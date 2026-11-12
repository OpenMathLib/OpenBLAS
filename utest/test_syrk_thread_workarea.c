/*****************************************************************************
Copyright (c) 2026, The OpenBLAS Project
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.
   3. Neither the name of the OpenBLAS project nor the names of
      its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********************************************************************************/

/* Issue #6224: every thread of the threaded SYRK/HERK driver
 * (driver/level3/level3_syrk_threaded.c) packs op(A) for its whole row range
 * and one K block into the sb part of its work area.  Nothing limited that to
 * the space left in the BUFFER_SIZE work area, so for large n it wrote past
 * the end of the work area into whatever was mapped next.
 *
 * This test calls the threaded driver directly and gives the calling thread
 * (job 0, which owns the widest row range of the lower-triangular partition)
 * a work area of its own, laid out like the one interface/syrk.c gets from
 * blas_memory_alloc(): BUFFER_SIZE bytes followed by WA_TAIL bytes that must
 * stay zero.  n is the smallest size at which the old driver writes
 * WA_MARGIN bytes past the end.  With k = 1 and ldc = 0 (every column of C
 * aliases the same n elements, so the values left in C are meaningless; the
 * packing of op(A), which overflowed, does not depend on C) the call costs
 * about n*n/2 multiply-adds and a few MB of memory although n is 10^4..10^5.
 * It is skipped where it would take more than OPENBLAS_UTEST_TIME_LIMIT
 * seconds (default 10, 0 = no limit), e.g. under qemu.
 *
 * It needs the internal driver symbols, i.e. a static libopenblas
 * (utest/CMakeLists.txt adds it only with BUILD_STATIC_LIBS). */

#include "openblas_utest.h"
#include <cblas.h>
#include <stdio.h>
#include <string.h>

#if defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3) && (MAX_CPU_NUMBER >= 2) && \
    (defined(BUILD_DOUBLE) || defined(BUILD_COMPLEX16))

#if !defined(_WIN32) || defined(__CYGWIN__)
#include <sys/time.h>
#endif

#define WA_MARGIN    (64 * 1024)   /* the old driver writes at least this far past the end */
#define WA_TAIL      (4 << 20)     /* zero memory behind the work area that is checked */
#define WA_ALIGN     (64 * 1024)   /* blas_memory_alloc() returns at least page-aligned areas */
#define WA_MAX_MADDS 4e10          /* never attempt more multiply-adds than this */

typedef int (*wa_driver_t)(blas_arg_t *, BLASLONG *, BLASLONG *, double *, double *, BLASLONG);

typedef struct {
  const char *name;
  wa_driver_t driver;            /* ?syrk_thread_LT */
  BLASLONG p, q, mn, cs;         /* GEMM_P, GEMM_Q, GEMM_UNROLL_MN, COMPSIZE */
  BLASLONG sb_off;               /* sb - start of work area, as in interface/syrk.c */
  BLASLONG n;                    /* size used */
  double madds;                  /* real multiply-adds of one call */
} wa_case_t;

static double wa_now(void)
{
#if defined(_WIN32) && !defined(__CYGWIN__)
  LARGE_INTEGER f, c;
  QueryPerformanceFrequency(&f); QueryPerformanceCounter(&c);
  return (double)c.QuadPart / (double)f.QuadPart;
#else
  /* not clock_gettime(): with glibc before 2.17 (manylinux1) that needs -lrt */
  struct timeval t;
  gettimeofday(&t, NULL);
  return (double)t.tv_sec + 1e-6 * (double)t.tv_usec;
#endif
}

/* End (bytes from the start of the work area) of what job 0 of the old
   driver writes: lower partition of level3_syrk_threaded.c with 2 threads,
   job 0 owns rows [0, w); with k = 1 its panel 0 (div_n columns) is at sb
   and panel 1 (w - div_n columns) at sb + GEMM_Q * div_n elements.  Keep
   this in sync with the partition and the panel layout of the driver. */
static double wa_old_end(const wa_case_t *t, BLASLONG n)
{
  BLASLONG mask = t->mn - 1, w, div_n;
  double dnum = (double)n * (double)n / 2.0, e = (double)t->cs * sizeof(double);
  w = (((BLASLONG)((sqrt(dnum) - 0.0) + mask)) / (mask + 1)) * (mask + 1);
  if (w > n || w < mask) w = n;
  div_n = (((w + 1) / 2 + mask) / (mask + 1)) * (mask + 1);
  if (w <= div_n) return (double)t->sb_off + (double)w * e;
  return (double)t->sb_off + ((double)t->q * (double)div_n + (double)(w - div_n)) * e;
}

/* Returns 0 if the parameters look inconsistent, e.g. because this file was
   compiled with other flags than the library that change the layout of the
   gotoblas table (BUILD_BFLOAT16/BUILD_HFLOAT16, SMALL_MATRIX_OPT, ...),
   from which DYNAMIC_ARCH builds read them. */
static int wa_setup(wa_case_t *t)
{
  BLASLONG n, i;
  t->n = 0; t->madds = 1e300;
  if (t->p < 1 || t->p > 65536 || t->q < 1 || t->q > 65536 || t->mn < 1 || t->mn > 1024) return 0;
  t->sb_off = (BLASLONG)GEMM_OFFSET_A + (BLASLONG)GEMM_OFFSET_B
            + (BLASLONG)(((size_t)t->p * t->q * t->cs * sizeof(double) + GEMM_ALIGN) & ~(size_t)GEMM_ALIGN);
  if (t->sb_off < 0 || t->sb_off >= (BLASLONG)BUFFER_SIZE / 2) return 0;
  n = (BLASLONG)(2.0 * sqrt(2.0) * ((double)BUFFER_SIZE - t->sb_off) / ((double)(t->q + 1) * t->cs * sizeof(double)) * 0.9);
  for (i = 0; i < 1000000 && wa_old_end(t, n) < (double)BUFFER_SIZE + WA_MARGIN; i++) n++;
  if (wa_old_end(t, n) < (double)BUFFER_SIZE + WA_MARGIN || wa_old_end(t, n) > (double)BUFFER_SIZE + WA_TAIL / 2) return 0;
  t->n = n;
  t->madds = (double)n * n / 2.0 * t->cs * t->cs;
  return 1;
}

static double wa_run(const wa_case_t *t, BLASLONG n, double *a, double *c, char *area)
{
  blas_arg_t args;
  double alpha[2] = {1.0, 0.0}, t0;
  memset(&args, 0, sizeof(args));
  args.n = n; args.k = 1;
  args.a = a; args.lda = 1;        /* op(A) = A^T, A is 1 x n */
  args.c = c; args.ldc = 0;        /* all columns alias c[0 .. n-1] */
  args.alpha = alpha; args.beta = NULL;
  args.nthreads = 2;
  t0 = wa_now();
  t->driver(&args, NULL, NULL, (double *)(area + GEMM_OFFSET_A), (double *)(area + t->sb_off), 0);
  return wa_now() - t0;
}

#ifdef USE_OPENMP
/* omp_get_thread_limit() is OpenMP 3.0; the OpenMP 2.0 runtime of MSVC does not have it */
static int wa_omp_thread_limit(void)
{
#if defined(_OPENMP) && (_OPENMP >= 200805)
  return omp_get_thread_limit();
#else
  return 2;
#endif
}
#endif

CTEST(syrk_thread, work_area_bound)
{
  wa_case_t cases[2], *t = NULL;
  int ncases = 0, used = 0, i, run;
#ifndef USE_OPENMP
  int saved;
#endif
  size_t o;
  double limit = 10.0, est = 0.0, secs = 0.0, *a, *c;
  long past = -1;
  char *raw, *area, *env;

#ifdef BUILD_DOUBLE
  cases[ncases].name = "dsyrk"; cases[ncases].driver = (wa_driver_t)dsyrk_thread_LT;
  cases[ncases].p = DGEMM_P; cases[ncases].q = DGEMM_Q; cases[ncases].mn = DGEMM_UNROLL_MN; cases[ncases].cs = 1;
  (void)wa_setup(&cases[ncases]); ncases++;
#endif
#ifdef BUILD_COMPLEX16
  cases[ncases].name = "zsyrk"; cases[ncases].driver = (wa_driver_t)zsyrk_thread_LT;
  cases[ncases].p = ZGEMM_P; cases[ncases].q = ZGEMM_Q; cases[ncases].mn = ZGEMM_UNROLL_MN; cases[ncases].cs = 2;
  (void)wa_setup(&cases[ncases]); ncases++;
#endif
  for (i = 0; i < ncases; i++) if (cases[i].n > 0 && (!t || cases[i].madds < t->madds)) t = &cases[i];
  if (!t) {
    CTEST_LOG("skipped: implausible GEMM parameters (built with other flags than the library?)");
    return;
  }

  env = getenv("OPENBLAS_UTEST_TIME_LIMIT");
  if (env) limit = atof(env);
  if (limit > 0 && t->madds > WA_MAX_MADDS) {
    CTEST_LOG("skipped: %s with n=%ld needs %.1e multiply-adds", t->name, (long)t->n, t->madds);
    return;
  }

  /* the two jobs wait for each other, so both must really run at the same time */
  if (openblas_get_num_procs() < 2) {
    CTEST_LOG("skipped: fewer than two CPUs");
    return;
  }
#ifdef USE_OPENMP
  /* no openblas_set_num_threads() here: in OpenMP builds that would make OpenBLAS ignore
     omp_set_num_threads() for the rest of the process */
  if (omp_get_max_threads() < 2 || wa_omp_thread_limit() < 2 || omp_in_parallel() || omp_get_dynamic()) {
    CTEST_LOG("skipped: two OpenMP threads are not available");
    return;
  }
#else
  saved = openblas_get_num_threads();
  openblas_set_num_threads(2);
  if (openblas_get_num_threads() < 2) {
    openblas_set_num_threads(saved);
    CTEST_LOG("skipped: two threads are not available");
    return;
  }
#endif

  raw = (char *)calloc(1, (size_t)BUFFER_SIZE + WA_TAIL + WA_ALIGN);   /* zero, mostly never touched */
  a = (double *)malloc((size_t)t->n * t->cs * sizeof(double));
  c = (double *)malloc((size_t)t->n * t->cs * sizeof(double));
  if (!raw || !a || !c) {
    free(raw); free(a); free(c);
#ifndef USE_OPENMP
    openblas_set_num_threads(saved);
#endif
    CTEST_LOG("skipped: not enough memory");
    return;
  }
  area = (char *)(((BLASULONG)raw + WA_ALIGN - 1) & ~(BLASULONG)(WA_ALIGN - 1));
  for (o = 0; o < (size_t)(t->n * t->cs); o++) a[o] = (double)(1 + o % 3);   /* no zeros */

  /* 1/64 of the work, far inside the work area (overestimates natively, where
     small sizes run less efficiently, by up to about 3x) */
  est = 64.0 * wa_run(t, t->n / 8, a, c, area);
  run = !(limit > 0 && est > limit);
  if (run) {
    memset(area + t->sb_off, 0, 64);
    secs = wa_run(t, t->n, a, c, area);
    for (o = 0; o < WA_TAIL; o++)
      if (area[(size_t)BUFFER_SIZE + o]) { past = (long)o; break; }
    for (i = 0; i < 64; i++) used |= area[t->sb_off + i] != 0;
  }
  free(raw); free(a); free(c);
#ifndef USE_OPENMP
  openblas_set_num_threads(saved);
#endif

  if (!run) {
    CTEST_LOG("skipped: %s with n=%ld would take about %.0f s", t->name, (long)t->n, est);
    return;
  }
  CTEST_LOG("%s n=%ld k=1, 2 threads, BUFFER_SIZE=%lu: %.2f s", t->name, (long)t->n, (unsigned long)BUFFER_SIZE, secs);
  if (past >= 0)
    CTEST_LOG("the threaded driver wrote past the end of its work area (first changed byte at +%ld)", past);
  ASSERT_EQUAL(-1, past);   /* nothing written behind the BUFFER_SIZE work area */
  ASSERT_TRUE(used);        /* and the work area was used at all */
}

#endif
