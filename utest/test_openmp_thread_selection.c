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
   3. Neither the name of the OpenBLAS project nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#include <cblas.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int callback_calls;
static int callback_numjobs = 1;

static void test_threads_callback(int sync, openblas_dojob_callback dojob,
                                  int numjobs, size_t jobdata_elsize,
                                  void *jobdata, int dojob_data)
{
  int i;

  (void)sync;
  callback_calls++;
  callback_numjobs = numjobs;
  if (omp_in_parallel()) {
    /*
     * Use tasks from the existing team. DGEMM jobs exchange packed blocks,
     * so they must be allowed to make progress concurrently.
     */
    for (i = 0; i < numjobs; i++) {
#pragma omp task firstprivate(i, dojob, jobdata_elsize, jobdata, dojob_data)
      dojob(i, (char *)jobdata + (size_t)i * jobdata_elsize, dojob_data);
    }
#pragma omp taskwait
  } else {
#pragma omp parallel for num_threads(numjobs) schedule(static, 1) \
    firstprivate(dojob, jobdata_elsize, jobdata, dojob_data)
    for (i = 0; i < numjobs; i++)
      dojob(i, (char *)jobdata + (size_t)i * jobdata_elsize, dojob_data);
  }
}

static int check_dgemm(const char *label, int expected_jobs, blasint n,
                       const double *a, const double *b, double *c)
{
  size_t i;
  size_t elements = (size_t)n * (size_t)n;
  int failed = 0;

  memset(c, 0, elements * sizeof(*c));
  callback_calls = 0;
  callback_numjobs = 1;

  printf("Running %s\n", label);
  fflush(stdout);

  /*
   * TN avoids architecture-specific small-matrix kernels that would bypass
   * the threaded DGEMM dispatcher at this deliberately small test size.
   */
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, n, n, n,
              1.0, a, n, b, n, 0.0, c, n);

  if (callback_numjobs != expected_jobs) {
    fprintf(stderr, "%s: expected %d logical job(s), observed %d\n",
            label, expected_jobs, callback_numjobs);
    failed = 1;
  }
  if (expected_jobs > 1 && callback_calls == 0) {
    fprintf(stderr, "%s: threaded callback was not invoked\n", label);
    failed = 1;
  }

  for (i = 0; i < elements; i++) {
    if (c[i] != (double)n) {
      fprintf(stderr,
              "%s: incorrect DGEMM result at element %zu: %.17g != %d\n",
              label, i, c[i], (int)n);
      failed = 1;
      break;
    }
  }

  return failed;
}

int main(void)
{
  const double two_thread_work =
      2.0 * 65536.0 * (double)GEMM_MULTITHREAD_THRESHOLD;
  blasint n = 1;
  size_t elements;
  double *a;
  double *b;
  double *c;
  int failed = 0;
  int nested_failed = 0;
  int outer_team_size = 0;
  size_t i;

  /*
   * This must be the first OpenBLAS call. The callback records the logical
   * job count and runs the jobs without creating a nested OpenMP team.
   */
  openblas_set_threads_callback_function(test_threads_callback);

  omp_set_dynamic(0);
  while ((double)n * (double)n * (double)n <= two_thread_work)
    n++;
  n++; /* Keep the test just beyond the threshold boundary. */

  elements = (size_t)n * (size_t)n;
  a = malloc(elements * sizeof(*a));
  b = malloc(elements * sizeof(*b));
  c = malloc(elements * sizeof(*c));
  if (a == NULL || b == NULL || c == NULL) {
    fprintf(stderr, "failed to allocate DGEMM test matrices\n");
    free(c);
    free(b);
    free(a);
    return 1;
  }

  for (i = 0; i < elements; i++)
    a[i] = b[i] = 1.0;

  openblas_set_num_threads(1);
  failed |= check_dgemm("global count 1 outside OpenMP", 1, n, a, b, c);

  openblas_set_num_threads(2);
  failed |= check_dgemm("global count 2 outside OpenMP", 2, n, a, b, c);

#pragma omp parallel num_threads(2) shared(outer_team_size, nested_failed)
  {
#pragma omp single
    {
      outer_team_size = omp_get_num_threads();
      if (outer_team_size == 2) {
#pragma omp task
        nested_failed =
            check_dgemm("global count 2 inside OpenMP task", 1, n, a, b, c);
#pragma omp taskwait
      }
    }
  }

  if (outer_team_size == 2)
    failed |= nested_failed;
  else
    printf("SKIP: OpenMP runtime formed an outer team of %d thread(s)\n",
           outer_team_size);

  failed |= check_dgemm("global count 2 after OpenMP region", 2, n, a, b, c);

  free(c);
  free(b);
  free(a);

  if (failed == 0)
    printf("OpenMP thread selection test passed (DGEMM dimension %d)\n",
           (int)n);
  return failed;
}
