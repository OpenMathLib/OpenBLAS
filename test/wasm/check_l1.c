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

static void check_saxpy_one(int n, int incx, int incy) {
  int absx = incx == 0 ? 1 : (incx > 0 ? incx : -incx);
  int absy = incy == 0 ? 1 : (incy > 0 ? incy : -incy);
  int nx = (n <= 0) ? 1 : (incx == 0 ? 1 : n * absx);
  int ny = (n <= 0) ? 1 : (incy == 0 ? 1 : n * absy);
  float *x = xmalloc((size_t)nx * sizeof(float));
  float *y = xmalloc((size_t)ny * sizeof(float));
  float *yr = xmalloc((size_t)ny * sizeof(float));
  fill_f32(x, nx, 1);
  fill_f32(y, ny, 2);
  memcpy(yr, y, (size_t)ny * sizeof(float));
  float alpha = 1.1f;
  ref_saxpy(n, alpha, x, incx, yr, incy);
  cblas_saxpy(n, alpha, x, incx, y, incy);
  float maxe, maxv;
  char msg[128];
  snprintf(msg, sizeof(msg), "saxpy n=%d incx=%d incy=%d", n, incx, incy);
  if (!close_f32(y, yr, ny, TOL_S_L1, &maxe, &maxv))
    fail_f32(msg, maxe, maxv, TOL_S_L1);
  else
    pass_one();
  free(x);
  free(y);
  free(yr);
}

static void check_daxpy_one(int n, int incx, int incy) {
  int absx = incx == 0 ? 1 : (incx > 0 ? incx : -incx);
  int absy = incy == 0 ? 1 : (incy > 0 ? incy : -incy);
  int nx = (n <= 0) ? 1 : (incx == 0 ? 1 : n * absx);
  int ny = (n <= 0) ? 1 : (incy == 0 ? 1 : n * absy);
  double *x = xmalloc((size_t)nx * sizeof(double));
  double *y = xmalloc((size_t)ny * sizeof(double));
  double *yr = xmalloc((size_t)ny * sizeof(double));
  fill_f64(x, nx, 1);
  fill_f64(y, ny, 2);
  memcpy(yr, y, (size_t)ny * sizeof(double));
  double alpha = 1.1;
  ref_daxpy(n, alpha, x, incx, yr, incy);
  cblas_daxpy(n, alpha, x, incx, y, incy);
  double maxe, maxv;
  char msg[128];
  snprintf(msg, sizeof(msg), "daxpy n=%d incx=%d incy=%d", n, incx, incy);
  if (!close_f64(y, yr, ny, TOL_D_L1, &maxe, &maxv))
    fail_f64(msg, maxe, maxv, TOL_D_L1);
  else
    pass_one();
  free(x);
  free(y);
  free(yr);
}

void check_l1(void) {
  printf("==> L1 AXPY\n");
  for (int s = 0; s < NS_L1; s++) {
    int n = SIZES_L1[s];
    check_saxpy_one(n, 1, 1);
    check_daxpy_one(n, 1, 1);
    for (int i = 0; i < NINCS; i++) {
      int inc = INCS[i];
      if (inc == 1)
        continue;
      check_saxpy_one(n, inc, inc);
      check_daxpy_one(n, inc, inc);
    }
  }
}
