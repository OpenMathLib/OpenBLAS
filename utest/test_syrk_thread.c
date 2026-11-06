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

/*
 * Exactness tests for the threaded SYRK/HERK driver (driver/level3/level3_syrk_threaded.c)
 * and for threaded ?POTRF, which uses the same driver for its trailing updates (#6224).
 *
 * All data are small integers, so every product and partial sum is exact in single and double
 * precision and the threaded result must equal the exact reference.  The strict other triangle,
 * the padding rows (ldc > n), and guard zones before and after C and A must be bit-for-bit
 * unchanged.  The sizes exercise several K blocks (including the halving of the last two), the
 * sqrt-shaped row partition with widths that are not multiples of GEMM_UNROLL_MN, ranges with a
 * single DIVIDE_RATE chunk or chunks of unequal length, fewer ranges than threads, and thread
 * ranges wider than GEMM_P.
 */

#include "openblas_utest.h"
#include <stddef.h>
#include <string.h>
#include <cblas.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif

/* The threaded driver is not built with USE_SIMPLE_THREADED_LEVEL3 (always set on mips64). */
#if defined(SMP) && !defined(USE_SIMPLE_THREADED_LEVEL3)

#define ST_GUARD 4096
#define ST_SENT  0xA7

/* ---- thread control ------------------------------------------------------------------- */

static int st_saved, st_tmax;

static void st_begin(void)
{
  int procs = openblas_get_num_procs();
#ifdef USE_OPENMP
  st_saved = omp_get_max_threads();
  st_tmax = st_saved;
  /* a smaller team than requested would leave the jobs of the driver waiting for each other */
  if (omp_get_dynamic() || omp_in_parallel()) st_tmax = 1;
#else
  st_saved = openblas_get_num_threads();
  st_tmax = MAX_CPU_NUMBER;
#endif
  /* the threads of the level-3 drivers wait for each other: do not oversubscribe */
  if (st_tmax > procs) st_tmax = procs;
  if (st_tmax > 8) st_tmax = 8;
  if (st_tmax < 2) CTEST_LOG("fewer than two threads available: no threaded case runs");
}

/* Use T threads; returns 0 if that is not possible in this build or on this machine. */
static int st_set(int T)
{
  if (T > st_tmax) return 0;
#ifdef USE_OPENMP
  /* not openblas_set_num_threads(): that would make OpenBLAS ignore omp_set_num_threads() for the
     rest of the process, and the OpenMP server can only run as many jobs as its team has threads */
  omp_set_num_threads(T);
  return T;
#else
  openblas_set_num_threads(T);
  return openblas_get_num_threads() == T ? T : 0;
#endif
}

static void st_end(void)
{
#ifdef USE_OPENMP
  omp_set_num_threads(st_saved);
#else
  openblas_set_num_threads(st_saved);
#endif
}

/* ---- data ----------------------------------------------------------------------------- */

static unsigned int st_seed;

static int st_rand(void)        /* integers in [-3, 3] */
{
  st_seed = st_seed * 1103515245u + 12345u;
  return (int)((st_seed >> 16) % 7u) - 3;
}

typedef struct { char name[3]; int dbl, cplx, herk; } st_prec;

/* the routines built in this configuration */
static const st_prec st_routines[] = {
#ifdef BUILD_SINGLE
  {"s", 0, 0, 0},
#endif
#ifdef BUILD_DOUBLE
  {"d", 1, 0, 0},
#endif
#ifdef BUILD_COMPLEX
  {"c", 0, 1, 0}, {"ch", 0, 1, 1},
#endif
#ifdef BUILD_COMPLEX16
  {"z", 1, 1, 0}, {"zh", 1, 1, 1},
#endif
  {"", 0, 0, 0}
};

static size_t st_esz(const st_prec *p) { return (size_t)(p->dbl ? 8 : 4) * (p->cplx ? 2 : 1); }

static double st_get(const st_prec *p, const unsigned char *b, size_t i)    /* i counts reals */
{
  return p->dbl ? ((const double *)b)[i] : (double)((const float *)b)[i];
}

static void st_put(const st_prec *p, unsigned char *b, size_t i, double v)
{
  if (p->dbl) ((double *)b)[i] = v; else ((float *)b)[i] = (float)v;
}

/* bytes + guard zones on both sides; returns the start of the usable area */
static unsigned char *st_alloc(size_t bytes)
{
  unsigned char *m = (unsigned char *)malloc(bytes + 2 * ST_GUARD);
  if (m == NULL) return NULL;
  memset(m, ST_SENT, bytes + 2 * ST_GUARD);
  return m + ST_GUARD;
}

static void st_free(unsigned char *b) { if (b) free(b - ST_GUARD); }

static int st_guards_ok(const unsigned char *b, size_t bytes)
{
  size_t i;
  for (i = 0; i < ST_GUARD; i++)
    if (b[(ptrdiff_t)i - ST_GUARD] != ST_SENT || b[bytes + i] != ST_SENT) return 0;
  return 1;
}

static long st_nfail;

#define ST_FAIL(...) do { if (st_nfail++ < 10) fprintf(stderr, __VA_ARGS__); } while (0)

/* ---- one SYRK/HERK case ----------------------------------------------------------------- */

static void st_call(const st_prec *p, char uplo, char trans, blasint n, blasint k, const double *alpha,
                    unsigned char *a, blasint lda, const double *beta, unsigned char *c, blasint ldc)
{
  float fa[2], fb[2];
  double da[2], db[2];
  fa[0] = (float)alpha[0]; fa[1] = (float)alpha[1]; fb[0] = (float)beta[0]; fb[1] = (float)beta[1];
  da[0] = alpha[0]; da[1] = alpha[1]; db[0] = beta[0]; db[1] = beta[1];
  (void)fa; (void)fb; (void)da; (void)db;
#ifdef BUILD_SINGLE
  if (!p->dbl && !p->cplx) BLASFUNC(ssyrk)(&uplo, &trans, &n, &k, fa, (float *)a, &lda, fb, (float *)c, &ldc);
#endif
#ifdef BUILD_DOUBLE
  if (p->dbl && !p->cplx) BLASFUNC(dsyrk)(&uplo, &trans, &n, &k, da, (double *)a, &lda, db, (double *)c, &ldc);
#endif
#ifdef BUILD_COMPLEX
  if (!p->dbl && p->cplx && !p->herk) BLASFUNC(csyrk)(&uplo, &trans, &n, &k, fa, (float *)a, &lda, fb, (float *)c, &ldc);
  if (!p->dbl && p->herk) BLASFUNC(cherk)(&uplo, &trans, &n, &k, fa, (float *)a, &lda, fb, (float *)c, &ldc);
#endif
#ifdef BUILD_COMPLEX16
  if (p->dbl && p->cplx && !p->herk) BLASFUNC(zsyrk)(&uplo, &trans, &n, &k, da, (double *)a, &lda, db, (double *)c, &ldc);
  if (p->dbl && p->herk) BLASFUNC(zherk)(&uplo, &trans, &n, &k, da, (double *)a, &lda, db, (double *)c, &ldc);
#endif
}

/* alpha and beta for variant v: beta == 1 (no scaling pass), beta == 0, and general integers;
   complex alpha and beta for ?SYRK, real ones for ?HERK */
static void st_scalars(const st_prec *p, int v, double *alpha, double *beta)
{
  static const double ar[4] = {1, 2, -1, 3}, ai[4] = {0, -1, 2, 1}, br[4] = {1, -2, 0, 2}, bi[4] = {0, 1, 0, -1};
  alpha[0] = ar[v & 3]; beta[0] = br[v & 3];
  alpha[1] = (p->cplx && !p->herk) ? ai[v & 3] : 0;
  beta[1]  = (p->cplx && !p->herk) ? bi[v & 3] : 0;
}

/* C := alpha*op(A)*op(A)**T (or **H) + beta*C with T threads; returns 0 if T threads are not available */
static int st_syrk_case(const st_prec *p, char uplo, char trans, blasint n, blasint k, int T,
                        const double *alpha, const double *beta, unsigned int seed)
{
  int w = p->cplx ? 2 : 1, upper = (uplo == 'U'), notrans = (trans == 'N');
  blasint rows = notrans ? n : k, cols = notrans ? k : n;
  blasint lda = rows + 3, ldc = n + 2, i, j, l;
  size_t es = st_esz(p), abytes = (size_t)lda * cols * es, cbytes = (size_t)ldc * n * es;
  unsigned char *a, *c, *a0, *c0;
  int *xr, *xi;
  long before = st_nfail;

  if (st_set(T) == 0) return 0;
  a = st_alloc(abytes); c = st_alloc(cbytes);
  a0 = (unsigned char *)malloc(abytes); c0 = (unsigned char *)malloc(cbytes);
  xr = (int *)malloc(sizeof(int) * (size_t)n * k); xi = (int *)malloc(sizeof(int) * (size_t)n * k);
  if (!a || !c || !a0 || !c0 || !xr || !xi) { ST_FAIL("syrk_thread: out of memory\n"); goto done; }

  st_seed = seed;
  /* A: integers in [-3,3]; the padding rows hold 1000 and must not be used */
  for (j = 0; j < cols; j++)
    for (i = 0; i < lda; i++)
      for (l = 0; l < w; l++) st_put(p, a, ((size_t)i + (size_t)j * lda) * w + l, i < rows ? st_rand() : 1000.0);
  memcpy(a0, a, abytes);
  /* integer copy of op(A), row i contiguous, for the reference */
  for (i = 0; i < n; i++)
    for (l = 0; l < k; l++) {
      size_t e = notrans ? ((size_t)i + (size_t)l * lda) * w : ((size_t)l + (size_t)i * lda) * w;
      xr[(size_t)i * k + l] = (int)st_get(p, a, e);
      xi[(size_t)i * k + l] = p->cplx ? (int)st_get(p, a, e + 1) : 0;
    }
  /* C: integers in the referenced triangle (real diagonal for HERK), other values elsewhere */
  for (j = 0; j < n; j++)
    for (i = 0; i < ldc; i++)
      for (l = 0; l < w; l++) {
        int tri = i < n && (upper ? i <= j : i >= j);
        double v = tri ? st_rand() : (i < n ? -77777.0 : 55555.0);
        if (p->herk && i == j && l == 1) v = 0.0;
        st_put(p, c, ((size_t)i + (size_t)j * ldc) * w + l, v);
      }
  memcpy(c0, c, cbytes);

  st_call(p, uplo, trans, n, k, alpha, a, lda, beta, c, ldc);

  if (!st_guards_ok(a, abytes) || memcmp(a, a0, abytes)) ST_FAIL("  A or its guard zones modified\n");
  if (!st_guards_ok(c, cbytes)) ST_FAIL("  guard zone around C modified\n");
  for (j = 0; j < n; j++) {
    size_t col = (size_t)j * ldc * es;
    /* strict other triangle and padding rows of column j: bit-for-bit unchanged */
    if (upper ? memcmp(c + col + (j + 1) * es, c0 + col + (j + 1) * es, (ldc - j - 1) * es) != 0
              : (memcmp(c + col, c0 + col, j * es) != 0 || memcmp(c + col + n * es, c0 + col + n * es, (ldc - n) * es) != 0))
      ST_FAIL("  column %ld modified outside the referenced triangle\n", (long)j);
    for (i = upper ? 0 : j; i < (upper ? j + 1 : n); i++) {
      size_t e = ((size_t)i + (size_t)j * ldc) * w;
      const int *x_r = xr + (size_t)i * k, *x_i = xi + (size_t)i * k, *y_r = xr + (size_t)j * k, *y_i = xi + (size_t)j * k;
      int s1 = 0, s2 = 0, s3 = 0, s4 = 0;
      double sr, si, er, ei, cr = st_get(p, c0, e), ci = p->cplx ? st_get(p, c0, e + 1) : 0.0;
      for (l = 0; l < k; l++) { s1 += x_r[l] * y_r[l]; s2 += x_i[l] * y_i[l]; s3 += x_r[l] * y_i[l]; s4 += x_i[l] * y_r[l]; }
      if (!p->herk) { sr = s1 - s2; si = s3 + s4; }       /* x * y          */
      else if (notrans) { sr = s1 + s2; si = s4 - s3; }   /* x * conj(y)    */
      else { sr = s1 + s2; si = s3 - s4; }                /* conj(x) * y    */
      if (p->herk) {
        er = alpha[0] * sr + beta[0] * cr;
        ei = (i == j) ? 0.0 : alpha[0] * si + beta[0] * ci;
      } else if (p->cplx) {
        er = alpha[0] * sr - alpha[1] * si + beta[0] * cr - beta[1] * ci;
        ei = alpha[0] * si + alpha[1] * sr + beta[0] * ci + beta[1] * cr;
      } else {
        er = alpha[0] * sr + beta[0] * cr; ei = 0.0;
      }
      if (st_get(p, c, e) != er || (p->cplx && st_get(p, c, e + 1) != ei))
        ST_FAIL("  C(%ld,%ld) = (%g,%g), expected (%g,%g)\n", (long)i, (long)j,
                st_get(p, c, e), p->cplx ? st_get(p, c, e + 1) : 0.0, er, ei);
    }
  }
done:
  if (st_nfail > before)
    fprintf(stderr, "FAILED: %s%s uplo=%c trans=%c n=%ld k=%ld threads=%d alpha=(%g,%g) beta=(%g,%g)\n",
            p->name, p->herk ? "erk" : "syrk", uplo, trans, (long)n, (long)k, T, alpha[0], alpha[1], beta[0], beta[1]);
  st_free(a); st_free(c); free(a0); free(c0); free(xr); free(xi);
  return T;
}

static char st_trans(const st_prec *p, int t) { return t ? (p->herk ? 'C' : 'T') : 'N'; }

/* Every (precision, uplo, trans) variant with 2 threads and several K blocks.  The K values
   exceed GEMM_Q of the common cores (most of them twice), so the K loop runs full blocks and
   halves the last two; the k = 1100 cases cover GEMM_Q up to 1026 (POWER9, Apple M1, ...).
   n = 39 (lower), 65 (upper, complex) and 129 (upper, real) give two thread ranges for
   GEMM_UNROLL_MN up to 16, 16 and 32. */
CTEST(syrk_thread, k_blocks)
{
  const st_prec *p;
  int v;
  st_nfail = 0;
  st_begin();
  for (p = st_routines; p->name[0]; p++) {
    unsigned int r = (unsigned int)(p - st_routines);
    blasint kq = p->cplx ? (p->dbl ? 260 : 300) : (p->dbl ? 600 : 700);
    blasint nu = p->cplx ? 65 : 129;
    double alpha[2], beta[2];
    for (v = 0; v < 4; v++) {
      char uplo = (v & 2) ? 'L' : 'U';
      st_scalars(p, v, alpha, beta);
      st_syrk_case(p, uplo, st_trans(p, v & 1), uplo == 'L' ? 39 : nu, kq, 2, alpha, beta, 100u + 10u * v + r);
    }
    /* k = 1100: lower for the real types, alternating for the complex ones */
    st_scalars(p, v + r, alpha, beta);
    if (!p->cplx || (r & 1) == ((r >> 1) & 1))
      st_syrk_case(p, 'L', st_trans(p, r & 1), 39, 1100, 2, alpha, beta, 200u + r);
    else
      st_syrk_case(p, 'U', st_trans(p, r & 1), 65, 1100, 2, alpha, beta, 200u + r);
  }
  st_end();
  ASSERT_EQUAL(0, st_nfail);
}

/* Row partition sweep with 2, 3, 4 and 8 threads (as far as available).  (n+1)*n*k is about twice
   the interface threshold for threading and n >= 16 * threads.  n = 131, 139, 133 and 131 give at
   least two ranges for GEMM_UNROLL_MN up to 32 (both triangles); n = 55 and 71 give fewer ranges
   than threads, down to a single range.  Depending on GEMM_UNROLL_MN, the widths include ones that
   are not multiples of it, ranges with only one DIVIDE_RATE chunk and chunks of unequal length.
   The alpha = 0 cases run only the beta scaling of each thread's rows. */
CTEST(syrk_thread, partitions)
{
  static const int Ts[4] = {2, 3, 4, 8};
  static const blasint nbig[4] = {131, 139, 133, 131}, nsmall[4] = {0, 55, 71, 0};
  const st_prec *p;
  int t, u;
  st_nfail = 0;
  st_begin();
  for (t = 0; t < 4; t++) {
    int T = Ts[t], pass;
    for (pass = 0; pass < 2; pass++) {
      blasint n = pass ? nsmall[t] : nbig[t];
      blasint kr = (blasint)(2.0 * 439777.0 / ((double)n * (n + 1))) + 1;   /* real:    2 x 4 x 109944 */
      blasint kc = (blasint)(2.0 * 59297.0 / ((double)n * (n + 1))) + 1;    /* complex: 2 x 4 x 14824  */
      if (n == 0) continue;
      for (p = st_routines; p->name[0]; p++) {
        unsigned int r = (unsigned int)(p - st_routines);
        for (u = 0; u < 2; u++) {
          double alpha[2], beta[2];
          if (pass && u != (int)(r & 1)) continue;      /* small n: one triangle per routine */
          st_scalars(p, T + 2 * u + pass, alpha, beta);
          st_syrk_case(p, u ? 'L' : 'U', st_trans(p, (T + u + pass) & 1), n, p->cplx ? kc : kr, T, alpha, beta,
                       1000u + 100u * T + 20u * pass + 10u * u + r);
        }
      }
    }
    for (p = st_routines; p->name[0]; p++)
      if (p == st_routines || (p->herk && p->dbl)) {
        double alpha[2] = {0, 0}, beta[2] = {-2, 0};
        blasint n = nbig[t];
        st_syrk_case(p, (T & 1) ? 'U' : 'L', 'N', n, (blasint)(2.0 * 439777.0 / ((double)n * (n + 1))) + 1, T,
                     alpha, beta, 3000u + T);
      }
  }
  st_end();
  ASSERT_EQUAL(0, st_nfail);
}

/* Thread ranges wider than GEMM_P (n = 727 with 2 threads gives about 514 rows to the widest
   range), so each thread packs several row panels per K block and keeps the other threads'
   panels until its last one; k = 1 and 2 exercise the shortest K blocks. */
CTEST(syrk_thread, wide_ranges)
{
  const st_prec *p;
  st_nfail = 0;
  st_begin();
  for (p = st_routines; p->name[0]; p++) {
    double alpha[2], beta[2];
    char uplo, trans;
    blasint k;
    if (!p->cplx && !p->dbl)     { uplo = 'L'; trans = 'N'; k = 3; }
    else if (!p->cplx)           { uplo = 'U'; trans = 'T'; k = 2; }
    else if (!p->dbl && !p->herk) { uplo = 'U'; trans = 'N'; k = 1; }
    else if (p->dbl && p->herk)  { uplo = 'L'; trans = 'C'; k = 1; }
    else continue;                /* csyrk and zherk cover the complex drivers */
#if defined(ARCH_LOONGARCH64)
    /* LA464 (DYNAMIC_ARCH): DGEMM_P = 112 is not a multiple of DGEMM_UNROLL_MN, which breaks
       DSYRK row panels of GEMM_P rows with or without threads (a pre-existing bug outside the
       threaded driver) */
    if (p->dbl && !p->cplx) continue;
#endif
    st_scalars(p, (int)(p - st_routines) + 1, alpha, beta);
    st_syrk_case(p, uplo, trans, 727, k, 2, alpha, beta, 4000u + (unsigned)(p - st_routines));
  }
  st_end();
  ASSERT_EQUAL(0, st_nfail);
}

#if (BUFFER_SIZE <= (16 << 20)) && defined(BUILD_COMPLEX16)
/* #6224 through the public interface, where that is affordable: with a work area of 16 MiB or
   less (e.g. 32-bit x86) the widest of two thread ranges (about n/sqrt(2) rows) no longer fits
   one K block of GEMM_Q into the work area at n = 6600 (scaled with BUFFER_SIZE), so the driver
   has to shorten the K block (Nehalem: 209 instead of 256, K blocks of 209, 152 and 151).
   Without that it wrote several MB past the end of the work area (SIGSEGV or wrong results).  A(i,l) = u_i * v_l with small integers, so
   C(i,j) = u_i * conj(u_j) * sum_l v_l**2 exactly; every element of C is checked. */
static int st_small(unsigned int x, int m)          /* integer in [-(m/2), m/2] */
{
  x ^= x >> 16; x *= 0x7feb352dU; x ^= x >> 15; x *= 0x846ca68bU; x ^= x >> 16;
  return (int)(x % (unsigned int)m) - m / 2;
}

CTEST(syrk_thread, small_work_area)
{
  blasint n = (blasint)(6600L * (long)(BUFFER_SIZE >> 20) / 16), k = 512, i, j, l;
  double alpha = 1.0, beta = 0.0, sv = 0.0, *a, *c, *u, *v;
  char uplo = 'L', trans = 'N';
  long bad = 0;

  st_begin();
  if (st_set(2) == 0) { st_end(); CTEST_LOG("skipped: two threads are not available"); return; }
  a = (double *)malloc(sizeof(double) * 2 * (size_t)n * k);
  c = (double *)calloc((size_t)n * n * 2, sizeof(double));
  u = (double *)malloc(sizeof(double) * 2 * (size_t)n);
  v = (double *)malloc(sizeof(double) * (size_t)k);
  if (!a || !c || !u || !v) {
    free(a); free(c); free(u); free(v); st_end();
    CTEST_LOG("skipped: not enough memory");
    return;
  }
  for (i = 0; i < n; i++) {
    u[2 * i] = st_small(2u * (unsigned int)i + 1u, 81);
    u[2 * i + 1] = st_small(2u * (unsigned int)i + 2u, 61);
    if (u[2 * i] == 0 && u[2 * i + 1] == 0) u[2 * i] = 1;
  }
  for (l = 0; l < k; l++) { v[l] = st_small(100003u + (unsigned int)l, 7); if (v[l] == 0) v[l] = 1; sv += v[l] * v[l]; }
  for (l = 0; l < k; l++)
    for (i = 0; i < n; i++) {
      a[2 * ((size_t)l * n + i)] = u[2 * i] * v[l];
      a[2 * ((size_t)l * n + i) + 1] = u[2 * i + 1] * v[l];
    }

  BLASFUNC(zherk)(&uplo, &trans, &n, &k, &alpha, a, &n, &beta, c, &n);
  st_end();

  for (j = 0; j < n; j++)
    for (i = 0; i < n; i++) {
      const double *x = c + 2 * ((size_t)j * n + i);
      double er = 0.0, ei = 0.0;
      if (i >= j) {
        er = (u[2 * i] * u[2 * j] + u[2 * i + 1] * u[2 * j + 1]) * sv;
        if (i != j) ei = (u[2 * i + 1] * u[2 * j] - u[2 * i] * u[2 * j + 1]) * sv;
      }
      if (x[0] != er || x[1] != ei)
        if (bad++ < 10) fprintf(stderr, "  C(%ld,%ld) = (%g,%g), expected (%g,%g)\n", (long)i, (long)j, x[0], x[1], er, ei);
    }
  free(a); free(c); free(u); free(v);
  if (bad) fprintf(stderr, "FAILED: zherk uplo=L trans=N n=%ld k=%ld threads=2: %ld wrong elements\n", (long)n, (long)k, bad);
  ASSERT_EQUAL(0, bad);
}
#endif /* BUFFER_SIZE <= 16 MiB */

#ifndef NO_LAPACK
/* ---- threaded ?POTRF ------------------------------------------------------------------- */

/* A = L*L**T (L**H) with integer L, unit diagonal and entries in [-2,2]: the factorisation is
   exact, so ?POTRF must return exactly L (or U = L**T / L**H) and leave the other triangle,
   the padding row and the guard zones unchanged. */
static void st_potrf_case(const st_prec *p, char uplo, blasint n, int T, unsigned int seed)
{
  int w = p->cplx ? 2 : 1, upper = (uplo == 'U');
  blasint lda = n + 1, i, j, l, info = -1;
  size_t es = st_esz(p), bytes = (size_t)lda * n * es;
  double *L = (double *)calloc((size_t)n * n * 2, sizeof(double));
  unsigned char *a = st_alloc(bytes), *a0 = (unsigned char *)malloc(bytes);
  long before = st_nfail;

  if (st_set(T) == 0) { free(L); st_free(a); free(a0); return; }
  if (!L || !a || !a0) { ST_FAIL("potrf_thread: out of memory\n"); goto done; }
  st_seed = seed;
  for (j = 0; j < n; j++)
    for (i = j; i < n; i++) {
      double *x = L + 2 * ((size_t)i + (size_t)j * n);
      if (i == j) { x[0] = 1.0; x[1] = 0.0; }
      else { x[0] = st_rand() % 3; x[1] = p->cplx ? st_rand() % 3 : 0; }
    }
  for (j = 0; j < n; j++)
    for (i = 0; i < lda; i++) {
      size_t e = ((size_t)i + (size_t)j * lda) * w;
      blasint r = upper ? j : i, s = upper ? i : j;      /* (r, s): position in the lower triangle */
      double sr = 0, si = 0;
      if (!(i < n && (upper ? i <= j : i >= j))) {
        st_put(p, a, e, i < n ? -77777.0 : 55555.0);
        if (p->cplx) st_put(p, a, e + 1, 3.0);
        continue;
      }
      for (l = 0; l <= s; l++) {
        const double *x = L + 2 * ((size_t)r + (size_t)l * n), *y = L + 2 * ((size_t)s + (size_t)l * n);
        sr += x[0] * y[0] + x[1] * y[1];                /* L(r,l) * conj(L(s,l)) */
        si += x[1] * y[0] - x[0] * y[1];
      }
      st_put(p, a, e, sr);
      if (p->cplx) st_put(p, a, e + 1, upper ? -si : si);
    }
  memcpy(a0, a, bytes);

#ifdef BUILD_SINGLE
  if (!p->dbl && !p->cplx) BLASFUNC(spotrf)(&uplo, &n, (float *)a, &lda, &info);
#endif
#ifdef BUILD_DOUBLE
  if (p->dbl && !p->cplx) BLASFUNC(dpotrf)(&uplo, &n, (double *)a, &lda, &info);
#endif
#ifdef BUILD_COMPLEX
  if (!p->dbl && p->cplx) BLASFUNC(cpotrf)(&uplo, &n, (float *)a, &lda, &info);
#endif
#ifdef BUILD_COMPLEX16
  if (p->dbl && p->cplx) BLASFUNC(zpotrf)(&uplo, &n, (double *)a, &lda, &info);
#endif

  if (info != 0) ST_FAIL("  info = %ld\n", (long)info);
  if (!st_guards_ok(a, bytes)) ST_FAIL("  guard zone around A modified\n");
  for (j = 0; j < n; j++) {
    size_t col = (size_t)j * lda * es;
    if (upper ? memcmp(a + col + (j + 1) * es, a0 + col + (j + 1) * es, (lda - j - 1) * es) != 0
              : (memcmp(a + col, a0 + col, j * es) != 0 || memcmp(a + col + n * es, a0 + col + n * es, (lda - n) * es) != 0))
      ST_FAIL("  column %ld modified outside the referenced triangle\n", (long)j);
    for (i = upper ? 0 : j; i < (upper ? j + 1 : n); i++) {
      size_t e = ((size_t)i + (size_t)j * lda) * w;
      const double *x = upper ? L + 2 * ((size_t)j + (size_t)i * n) : L + 2 * ((size_t)i + (size_t)j * n);
      double er = x[0], ei = upper ? -x[1] : x[1];
      if (st_get(p, a, e) != er || (p->cplx && st_get(p, a, e + 1) != ei))
        ST_FAIL("  A(%ld,%ld) = (%g,%g), expected (%g,%g)\n", (long)i, (long)j,
                st_get(p, a, e), p->cplx ? st_get(p, a, e + 1) : 0.0, er, ei);
    }
  }
done:
  if (st_nfail > before) fprintf(stderr, "FAILED: %spotrf uplo=%c n=%ld threads=%d\n", p->name, uplo, (long)n, T);
  free(L); st_free(a); free(a0);
}

/* The sizes make the interface use the threaded driver (real: n >= 128 (s) / 64 (d) per thread,
   complex: n >= 64) and give trailing updates of at least 16 rows per thread, so that they go
   through the threaded SYRK/HERK driver (HERK_THREAD_LN for 'L', HERK_THREAD_UC for 'U'). */
CTEST(potrf_thread, exact)
{
  const st_prec *p;
  st_nfail = 0;
  st_begin();
  for (p = st_routines; p->name[0]; p++) {
    int T4 = st_tmax < 4 ? st_tmax : 4, T3 = st_tmax < 3 ? st_tmax : 3;
    if (p->herk || T3 < 2) continue;
    if (!p->cplx && !p->dbl) {
      st_potrf_case(p, 'L', 128 * 2 + 9, 2, 5000u);
    } else if (!p->cplx) {
      st_potrf_case(p, 'U', 64 * T3 + 13, T3, 5001u);
      st_potrf_case(p, 'L', 64 * T3 + 13, T3, 5002u);
    } else {
      st_potrf_case(p, 'U', 141, T4, 5003u + p->dbl);
      st_potrf_case(p, 'L', 141, T4, 5005u + p->dbl);
    }
  }
  st_end();
  ASSERT_EQUAL(0, st_nfail);
}
#endif /* NO_LAPACK */

#endif /* SMP && !USE_SIMPLE_THREADED_LEVEL3 */
