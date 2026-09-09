/*
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
      derived from this software without specific prior written permission.

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
*/

#include "cases.h"
#include "common.h"
#include "ref.h"
#include "tol.h"

static void check_sgemv_one(enum CBLAS_TRANSPOSE trans, int m, int n, int incx,
                            int incy) {
  int lda = m;
  int lenx = (trans == CblasNoTrans) ? n : m;
  int leny = (trans == CblasNoTrans) ? m : n;
  int absx = incx == 0 ? 1 : (incx > 0 ? incx : -incx);
  int absy = incy == 0 ? 1 : (incy > 0 ? incy : -incy);
  int nx = incx == 0 ? 1 : lenx * absx;
  int ny = incy == 0 ? 1 : leny * absy;
  float *A = xmalloc((size_t)lda * (size_t)n * sizeof(float));
  float *x = xmalloc((size_t)nx * sizeof(float));
  float *y = xmalloc((size_t)ny * sizeof(float));
  float *yr = xmalloc((size_t)ny * sizeof(float));
  fill_f32(A, lda * n, 3);
  fill_f32(x, nx, 4);
  fill_f32(y, ny, 5);
  memcpy(yr, y, (size_t)ny * sizeof(float));
  float alpha = 1.1f, beta = 0.7f;
  ref_sgemv(trans, m, n, alpha, A, lda, x, incx, beta, yr, incy);
  cblas_sgemv(CblasColMajor, trans, m, n, alpha, A, lda, x, incx, beta, y,
              incy);
  float maxe, maxv, tol = tol_s_l2(m > n ? m : n);
  char msg[160];
  snprintf(msg, sizeof(msg), "sgemv t=%d m=%d n=%d incx=%d incy=%d", (int)trans,
           m, n, incx, incy);
  if (!close_f32(y, yr, ny, tol, &maxe, &maxv))
    fail_f32(msg, maxe, maxv, tol);
  else
    pass_one();
  free(A);
  free(x);
  free(y);
  free(yr);
}

static void check_dgemv_one(enum CBLAS_TRANSPOSE trans, int m, int n, int incx,
                            int incy) {
  int lda = m;
  int lenx = (trans == CblasNoTrans) ? n : m;
  int leny = (trans == CblasNoTrans) ? m : n;
  int absx = incx == 0 ? 1 : (incx > 0 ? incx : -incx);
  int absy = incy == 0 ? 1 : (incy > 0 ? incy : -incy);
  int nx = incx == 0 ? 1 : lenx * absx;
  int ny = incy == 0 ? 1 : leny * absy;
  double *A = xmalloc((size_t)lda * (size_t)n * sizeof(double));
  double *x = xmalloc((size_t)nx * sizeof(double));
  double *y = xmalloc((size_t)ny * sizeof(double));
  double *yr = xmalloc((size_t)ny * sizeof(double));
  fill_f64(A, lda * n, 3);
  fill_f64(x, nx, 4);
  fill_f64(y, ny, 5);
  memcpy(yr, y, (size_t)ny * sizeof(double));
  double alpha = 1.1, beta = 0.7;
  ref_dgemv(trans, m, n, alpha, A, lda, x, incx, beta, yr, incy);
  cblas_dgemv(CblasColMajor, trans, m, n, alpha, A, lda, x, incx, beta, y,
              incy);
  double maxe, maxv, tol = tol_d_l2(m > n ? m : n);
  char msg[160];
  snprintf(msg, sizeof(msg), "dgemv t=%d m=%d n=%d incx=%d incy=%d", (int)trans,
           m, n, incx, incy);
  if (!close_f64(y, yr, ny, tol, &maxe, &maxv))
    fail_f64(msg, maxe, maxv, tol);
  else
    pass_one();
  free(A);
  free(x);
  free(y);
  free(yr);
}

void check_l2(void) {
  printf("==> L2 GEMV\n");
  enum CBLAS_TRANSPOSE tr[2] = {CblasNoTrans, CblasTrans};
  for (int s = 0; s < NS_L2; s++) {
    int n = SIZES_L2[s];
    /* Square and a few rectangular shapes. */
    int ms[3] = {n, n + (n > 1 ? 1 : 0), n > 2 ? n - 1 : n};
    int ns[3] = {n, n > 2 ? n - 1 : n, n + (n > 1 ? 2 : 0)};
    for (int t = 0; t < 2; t++) {
      for (int k = 0; k < 3; k++) {
        int m = ms[k];
        int nn = ns[k];
        if (m < 1 || nn < 1)
          continue;
        check_sgemv_one(tr[t], m, nn, 1, 1);
        check_dgemv_one(tr[t], m, nn, 1, 1);
        if (n <= 36) {
          check_sgemv_one(tr[t], m, nn, 2, 3);
          check_dgemv_one(tr[t], m, nn, 2, 3);
        }
      }
    }
  }
}
