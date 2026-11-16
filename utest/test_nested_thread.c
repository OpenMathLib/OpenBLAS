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
 * Regression tests for nested parallel work inside thread-server jobs.
 *
 * Threaded level-2 drivers call level-1 kernels (DOT, COPY) from inside each
 * worker job.  On several targets those level-1 kernels are themselves
 * threaded above 10000 elements via num_cpu_avail(1):
 *   - DDOT / ZDOT on x86_64 (Haswell, Zen, SkylakeX, Sandybridge, ...) and
 *     SDOT / DDOT / CDOT / ZDOT on arm64 (Neoverse, SVE, A64FX, Vortex, TX2);
 *   - SCOPY / DCOPY / CCOPY / ZCOPY on arm64 (copy_thunderx2t99.c).
 *
 * Without an in-job check in num_cpu_avail(), workers re-entered exec_blas()
 * and deadlocked in exec_blas_async() once all pool threads were occupied by
 * the outer call.
 */

#include "openblas_utest.h"
#include <cblas.h>
#include <stdlib.h>
#include <string.h>

#if !defined(_WIN32)
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif
#ifdef USE_OPENMP
#include <omp.h>
#endif

#ifdef SMP

/* Run fn(nthreads) in a child process under a watchdog when using the POSIX
   pthreads server so that a deadlock fails cleanly instead of hanging the test
   suite.  OpenMP and native Windows builds run fn() in-process. */
static __inline int run_nested_case(int (*fn)(int))
{
#if !defined(_WIN32)
  if (openblas_get_parallel() == OPENBLAS_THREAD) {
    pid_t pid;
    int status, i, limit = 20;
    const char *env = getenv("OPENBLAS_UTEST_TIMEOUT");
    int saved = openblas_get_num_threads();
    int procs = openblas_get_num_procs();

    if (env != NULL && atoi(env) > 0) limit = atoi(env);
    if (saved < 2 && procs < 2) return fn(1);

    fflush(stdout);
    fflush(stderr);
    pid = fork();
    if (pid < 0) return 0;
    if (pid == 0) {
      int t = saved >= 2 ? saved : 2;
      if (saved < 2) openblas_set_num_threads(t);
      exit(fn(t));
    }
    for (i = 0; i < limit * 100; i++) {
      if (waitpid(pid, &status, WNOHANG) == pid) {
        if (WIFEXITED(status)) return WEXITSTATUS(status);
        return 128 + WTERMSIG(status);
      }
      usleep(10000);
    }
    kill(pid, SIGKILL);
    waitpid(pid, &status, 0);
    return -1;
  }
#endif
  {
    int t = openblas_get_num_threads();
    if (t < 1) t = 1;
    if (t > 16) t = 16;
    return fn(t);
  }
}

#ifdef BUILD_DOUBLE

static int nested_dgbmv_dot(int nthreads)
{
  blasint kl = 10000, ku = 0, lda = kl + 1;
  blasint n = (blasint)(nthreads < 4 ? 16 : 4 * nthreads);
  blasint m = kl + n + 10001, inc = 1, j;
  size_t i, na = (size_t)lda * (size_t)n;
  double alpha = 1.0, beta = 0.0;
  double *a, *x, *y;
  char trans = 'T';
  int bad = 0;

  a = (double *)malloc(sizeof(double) * na);
  x = (double *)malloc(sizeof(double) * (size_t)m);
  y = (double *)calloc((size_t)n, sizeof(double));
  if (a == NULL || x == NULL || y == NULL) {
    free(a); free(x); free(y);
    return 0;
  }
  for (i = 0; i < na; i++) a[i] = 1.0;
  for (i = 0; i < (size_t)m; i++) x[i] = 1.0;

  BLASFUNC(dgbmv)(&trans, &m, &n, &kl, &ku, &alpha, a, &lda, x, &inc, &beta, y, &inc);

  for (j = 0; j < n; j++)
    if (y[j] != (double)(kl + 1)) bad = 1;

  free(a);
  free(x);
  free(y);
  return bad;
}

CTEST(nested_thread, dgbmv_dot)
{
  ASSERT_EQUAL(0, run_nested_case(nested_dgbmv_dot));
}

static int nested_dger_strided_copy(int nthreads)
{
  blasint m = 10001, n = (blasint)(nthreads < 2 ? 8 : 4 * nthreads);
  blasint incx = 2, incy = 1, lda = m, i, j;
  double alpha = 1.0;
  double *x, *y, *a;
  int bad = 0;

  x = (double *)malloc(sizeof(double) * (size_t)m * (size_t)incx);
  y = (double *)malloc(sizeof(double) * (size_t)n);
  a = (double *)calloc((size_t)lda * (size_t)n, sizeof(double));
  if (x == NULL || y == NULL || a == NULL) {
    free(x); free(y); free(a);
    return 0;
  }
  for (i = 0; i < m * incx; i++) x[i] = (double)(i % 7) - 3.0;
  for (j = 0; j < n; j++) y[j] = (double)(j + 1);

  BLASFUNC(dger)(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);

  for (j = 0; j < n; j++)
    for (i = 0; i < m; i++)
      if (a[i + (size_t)j * (size_t)lda] != x[i * incx] * y[j]) bad = 1;

  free(x);
  free(y);
  free(a);
  return bad;
}

CTEST(nested_thread, dger_strided_copy)
{
  ASSERT_EQUAL(0, run_nested_case(nested_dger_strided_copy));
}

static int nested_dtbmv_strided_copy(int nthreads)
{
  blasint n = 10001, k = 0, lda = 1, incx = 2, i;
  double *a, *x;
  char uplo = 'U', trans = 'N', diag = 'N';
  int bad = 0;
  (void)nthreads;

  a = (double *)malloc(sizeof(double) * (size_t)n);
  x = (double *)malloc(sizeof(double) * (size_t)n * (size_t)incx);
  if (a == NULL || x == NULL) {
    free(a); free(x);
    return 0;
  }
  for (i = 0; i < n; i++) a[i] = 2.0;
  for (i = 0; i < n * incx; i++) x[i] = (double)(i % 5) - 2.0;

  BLASFUNC(dtbmv)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);

  for (i = 0; i < n; i++) {
    double orig = (double)((i * incx) % 5) - 2.0;
    if (x[i * incx] != 2.0 * orig) bad = 1;
  }
  free(a);
  free(x);
  return bad;
}

CTEST(nested_thread, dtbmv_strided_copy)
{
  ASSERT_EQUAL(0, run_nested_case(nested_dtbmv_strided_copy));
}

#endif /* BUILD_DOUBLE */

#ifdef BUILD_SINGLE

static int nested_sgbmv_dot(int nthreads)
{
  blasint kl = 22001, ku = 0, lda = kl + 1;
  blasint n = (blasint)(nthreads < 4 ? 16 : 4 * nthreads);
  blasint m = kl + n + 1, inc = 1, j;
  size_t i, na = (size_t)lda * (size_t)n;
  float alpha = 1.0f, beta = 0.0f;
  float *a, *x, *y;
  char trans = 'T';
  int bad = 0;

  a = (float *)malloc(sizeof(float) * na);
  x = (float *)malloc(sizeof(float) * (size_t)m);
  y = (float *)calloc((size_t)n, sizeof(float));
  if (a == NULL || x == NULL || y == NULL) {
    free(a); free(x); free(y);
    return 0;
  }
  for (i = 0; i < na; i++) a[i] = 1.0f;
  for (i = 0; i < (size_t)m; i++) x[i] = 1.0f;

  BLASFUNC(sgbmv)(&trans, &m, &n, &kl, &ku, &alpha, a, &lda, x, &inc, &beta, y, &inc);

  for (j = 0; j < n; j++)
    if (y[j] != (float)(kl + 1)) bad = 1;

  free(a);
  free(x);
  free(y);
  return bad;
}

CTEST(nested_thread, sgbmv_dot)
{
  ASSERT_EQUAL(0, run_nested_case(nested_sgbmv_dot));
}

#endif /* BUILD_SINGLE */

#ifdef BUILD_COMPLEX16

static int nested_zgbmv_zdot(int nthreads)
{
  blasint kl = 10000, ku = 0, lda = kl + 1;
  blasint n = (blasint)(nthreads < 4 ? 16 : 4 * nthreads);
  blasint m = kl + n + 1, inc = 1, j;
  size_t i, na = (size_t)lda * (size_t)n;
  double alpha[2] = {1.0, 0.0}, beta[2] = {0.0, 0.0};
  double *a, *x, *y;
  char trans = 'C';
  int bad = 0;

  a = (double *)malloc(sizeof(double) * 2 * na);
  x = (double *)malloc(sizeof(double) * 2 * (size_t)m);
  y = (double *)calloc(2 * (size_t)n, sizeof(double));
  if (a == NULL || x == NULL || y == NULL) {
    free(a); free(x); free(y);
    return 0;
  }
  for (i = 0; i < na; i++) { a[2 * i] = 1.0; a[2 * i + 1] = 1.0; }
  for (i = 0; i < (size_t)m; i++) { x[2 * i] = 1.0; x[2 * i + 1] = -2.0; }

  BLASFUNC(zgbmv)(&trans, &m, &n, &kl, &ku, alpha, a, &lda, x, &inc, beta, y, &inc);

  /* conj(1 + i) * (1 - 2i) = (1 - i) * (1 - 2i) = -1 - 3i */
  for (j = 0; j < n; j++) {
    if (y[2 * j] != -(double)(kl + 1) || y[2 * j + 1] != -3.0 * (double)(kl + 1))
      bad = 1;
  }

  free(a);
  free(x);
  free(y);
  return bad;
}

CTEST(nested_thread, zgbmv_zdot)
{
  ASSERT_EQUAL(0, run_nested_case(nested_zgbmv_zdot));
}

#endif /* BUILD_COMPLEX16 */

#endif /* SMP */
