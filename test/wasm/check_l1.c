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

/*
 * Level-1 CBLAS checks vs the scalar oracle in ref_l1.c / ref.c.
 *
 * Deep SAXPY/DAXPY grids stress tile remainders and strides (SIZES_L1).
 * The remaining Level-1 ops run on the compact SIZES_FULL grid.
 *
 * One static test_* function per public CBLAS interface: allocate, fill,
 * call the scalar reference, call OpenBLAS, then expect_close_* / expect_eq_size.
 */

#include "cases.h"
#include "common.h"
#include "ref.h"
#include "tol.h"

/* Storage length for an n-vector with stride |inc| (inc==0 → length 1). */
static int vec_len(int n, int inc) {
  int abs_inc = inc == 0 ? 1 : (inc > 0 ? inc : -inc);
  if (n <= 0)
    return 1;
  return inc == 0 ? 1 : n * abs_inc;
}

/* ------------------------------------------------------------------------- */
/* Deep AXPY grids (tile remainders + non-unit strides)                      */
/* ------------------------------------------------------------------------- */

/* saxpy: y := alpha * x + y */
static void test_saxpy(int n, int incx, int incy) {
  int nx = vec_len(n, incx);
  int ny = vec_len(n, incy);
  float *x = xmalloc((size_t)nx * sizeof(float));
  float *y = xmalloc((size_t)ny * sizeof(float));
  float *yr = xmalloc((size_t)ny * sizeof(float));
  float alpha = 1.1f;
  char msg[128];

  fill_f32(x, nx, 1);
  fill_f32(y, ny, 2);
  memcpy(yr, y, (size_t)ny * sizeof(float));

  ref_saxpy(n, alpha, x, incx, yr, incy);
  cblas_saxpy(n, alpha, x, incx, y, incy);

  snprintf(msg, sizeof(msg), "saxpy n=%d incx=%d incy=%d", n, incx, incy);
  expect_close_f32(msg, y, yr, ny, TOL_S_L1);

  free(x);
  free(y);
  free(yr);
}

/* daxpy: y := alpha * x + y */
static void test_daxpy(int n, int incx, int incy) {
  int nx = vec_len(n, incx);
  int ny = vec_len(n, incy);
  double *x = xmalloc((size_t)nx * sizeof(double));
  double *y = xmalloc((size_t)ny * sizeof(double));
  double *yr = xmalloc((size_t)ny * sizeof(double));
  double alpha = 1.1;
  char msg[128];

  fill_f64(x, nx, 1);
  fill_f64(y, ny, 2);
  memcpy(yr, y, (size_t)ny * sizeof(double));

  ref_daxpy(n, alpha, x, incx, yr, incy);
  cblas_daxpy(n, alpha, x, incx, y, incy);

  snprintf(msg, sizeof(msg), "daxpy n=%d incx=%d incy=%d", n, incx, incy);
  expect_close_f64(msg, y, yr, ny, TOL_D_L1);

  free(x);
  free(y);
  free(yr);
}

/* ------------------------------------------------------------------------- */
/* Real single-precision Level-1 (compact grid, unit / stride-2)             */
/* ------------------------------------------------------------------------- */

/* sswap: swap x and y */
static void test_sswap(int n, int inc) {
  int z = n * inc;
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float *xr = xmalloc((size_t)z * sizeof(float));
  float *yr = xmalloc((size_t)z * sizeof(float));

  fill_f32(x, z, 20);
  fill_f32(y, z, 21);
  memcpy(xr, x, (size_t)z * sizeof(float));
  memcpy(yr, y, (size_t)z * sizeof(float));

  ref_l1_swap(n, xr, inc, yr, inc, 1, 0);
  cblas_sswap(n, x, inc, y, inc);

  expect_close_f32("sswap", x, xr, z, TOL_S_L1);
  expect_close_f32("sswap-y", y, yr, z, TOL_S_L1);

  free(x);
  free(y);
  free(xr);
  free(yr);
}

/* scopy: y := x */
static void test_scopy(int n, int inc) {
  int z = n * inc;
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float *yr = xmalloc((size_t)z * sizeof(float));

  fill_f32(x, z, 20);
  fill_f32(y, z, 22);
  memcpy(yr, y, (size_t)z * sizeof(float));

  ref_l1_copy(n, x, inc, yr, inc, 1, 0);
  cblas_scopy(n, x, inc, y, inc);

  expect_close_f32("scopy", y, yr, z, TOL_S_L1);

  free(x);
  free(y);
  free(yr);
}

/* sscal: x := alpha * x */
static void test_sscal(int n, int inc) {
  int z = n * inc;
  float *x = xmalloc((size_t)z * sizeof(float));
  float *xr = xmalloc((size_t)z * sizeof(float));
  float a = -0.7f;

  fill_f32(x, z, 20);
  memcpy(xr, x, (size_t)z * sizeof(float));

  ref_l1_scal(n, &a, xr, inc, 1, 0, 1);
  cblas_sscal(n, a, x, inc);

  expect_close_f32("sscal", x, xr, z, TOL_S_L1);

  free(x);
  free(xr);
}

/* saxpy on the compact grid (in addition to the deep grid above) */
static void test_saxpy_compact(int n, int inc) {
  int z = n * inc;
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float *yr = xmalloc((size_t)z * sizeof(float));
  float a = -0.7f;

  fill_f32(x, z, 23);
  fill_f32(y, z, 24);
  memcpy(yr, y, (size_t)z * sizeof(float));

  ref_l1_axpy(n, &a, x, inc, yr, inc, 1, 0);
  cblas_saxpy(n, a, x, inc, y, inc);

  expect_close_f32("saxpy-compact", y, yr, z, TOL_S_L1);

  free(x);
  free(y);
  free(yr);
}

/* sdot: dot product x^T y */
static void test_sdot(int n, int inc) {
  int z = n * inc;
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float dr = 0.0f;
  float dg;

  fill_f32(x, z, 23);
  fill_f32(y, z, 24);

  ref_l1_dot(n, x, inc, y, inc, &dr, 1, 0, 0);
  dg = cblas_sdot(n, x, inc, y, inc);

  expect_close_f32("sdot", &dg, &dr, 1, TOL_S_L1 * (n + 1));

  free(x);
  free(y);
}

/* snrm2: Euclidean norm of x */
static void test_snrm2(int n, int inc) {
  int z = n * inc;
  float *x = xmalloc((size_t)z * sizeof(float));
  float nr;
  float ng;

  fill_f32(x, z, 23);
  nr = (float)ref_l1_nrm2(n, x, inc, 1, 0);
  ng = cblas_snrm2(n, x, inc);

  expect_close_f32("snrm2", &ng, &nr, 1, TOL_S_L1 * (n + 1));

  free(x);
}

/* sasum: sum of absolute values of x */
static void test_sasum(int n, int inc) {
  int z = n * inc;
  float *x = xmalloc((size_t)z * sizeof(float));
  float sr;
  float sg;

  fill_f32(x, z, 23);
  sr = (float)ref_l1_asum(n, x, inc, 1, 0);
  sg = cblas_sasum(n, x, inc);

  expect_close_f32("sasum", &sg, &sr, 1, TOL_S_L1 * (n + 1));

  free(x);
}

/* isamax: 1-based index of max |x_i| */
static void test_isamax(int n, int inc) {
  int z = n * inc;
  float *x = xmalloc((size_t)z * sizeof(float));
  size_t ir;
  size_t ig;

  fill_f32(x, z, 23);
  ir = ref_l1_iamax(n, x, inc, 1, 0);
  ig = (size_t)cblas_isamax(n, x, inc);

  expect_eq_size("isamax", ig, ir);

  free(x);
}

/* ------------------------------------------------------------------------- */
/* Real double-precision Level-1                                             */
/* ------------------------------------------------------------------------- */

/* dswap: swap x and y */
static void test_dswap(int n, int inc) {
  int z = n * inc;
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double *xr = xmalloc((size_t)z * sizeof(double));
  double *yr = xmalloc((size_t)z * sizeof(double));

  fill_f64(x, z, 20);
  fill_f64(y, z, 21);
  memcpy(xr, x, (size_t)z * sizeof(double));
  memcpy(yr, y, (size_t)z * sizeof(double));

  ref_l1_swap(n, xr, inc, yr, inc, 1, 1);
  cblas_dswap(n, x, inc, y, inc);

  expect_close_f64("dswap", x, xr, z, TOL_D_L1);
  expect_close_f64("dswap-y", y, yr, z, TOL_D_L1);

  free(x);
  free(y);
  free(xr);
  free(yr);
}

/* dcopy: y := x */
static void test_dcopy(int n, int inc) {
  int z = n * inc;
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double *yr = xmalloc((size_t)z * sizeof(double));

  fill_f64(x, z, 20);
  fill_f64(y, z, 22);
  memcpy(yr, y, (size_t)z * sizeof(double));

  ref_l1_copy(n, x, inc, yr, inc, 1, 1);
  cblas_dcopy(n, x, inc, y, inc);

  expect_close_f64("dcopy", y, yr, z, TOL_D_L1);

  free(x);
  free(y);
  free(yr);
}

/* dscal: x := alpha * x */
static void test_dscal(int n, int inc) {
  int z = n * inc;
  double *x = xmalloc((size_t)z * sizeof(double));
  double *xr = xmalloc((size_t)z * sizeof(double));
  double a = -0.7;

  fill_f64(x, z, 20);
  memcpy(xr, x, (size_t)z * sizeof(double));

  ref_l1_scal(n, &a, xr, inc, 1, 1, 1);
  cblas_dscal(n, a, x, inc);

  expect_close_f64("dscal", x, xr, z, TOL_D_L1);

  free(x);
  free(xr);
}

/* daxpy on the compact grid */
static void test_daxpy_compact(int n, int inc) {
  int z = n * inc;
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double *yr = xmalloc((size_t)z * sizeof(double));
  double a = -0.7;

  fill_f64(x, z, 23);
  fill_f64(y, z, 24);
  memcpy(yr, y, (size_t)z * sizeof(double));

  ref_l1_axpy(n, &a, x, inc, yr, inc, 1, 1);
  cblas_daxpy(n, a, x, inc, y, inc);

  expect_close_f64("daxpy-compact", y, yr, z, TOL_D_L1);

  free(x);
  free(y);
  free(yr);
}

/* ddot: dot product x^T y */
static void test_ddot(int n, int inc) {
  int z = n * inc;
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double dr = 0.0;
  double dg;

  fill_f64(x, z, 23);
  fill_f64(y, z, 24);

  ref_l1_dot(n, x, inc, y, inc, &dr, 1, 1, 0);
  dg = cblas_ddot(n, x, inc, y, inc);

  expect_close_f64("ddot", &dg, &dr, 1, TOL_D_L1 * (n + 1));

  free(x);
  free(y);
}

/* dnrm2: Euclidean norm of x */
static void test_dnrm2(int n, int inc) {
  int z = n * inc;
  double *x = xmalloc((size_t)z * sizeof(double));
  double nr;
  double ng;

  fill_f64(x, z, 23);
  nr = ref_l1_nrm2(n, x, inc, 1, 1);
  ng = cblas_dnrm2(n, x, inc);

  expect_close_f64("dnrm2", &ng, &nr, 1, TOL_D_L1 * (n + 1));

  free(x);
}

/* dasum: sum of absolute values of x */
static void test_dasum(int n, int inc) {
  int z = n * inc;
  double *x = xmalloc((size_t)z * sizeof(double));
  double sr;
  double sg;

  fill_f64(x, z, 23);
  sr = ref_l1_asum(n, x, inc, 1, 1);
  sg = cblas_dasum(n, x, inc);

  expect_close_f64("dasum", &sg, &sr, 1, TOL_D_L1 * (n + 1));

  free(x);
}

/* idamax: 1-based index of max |x_i| */
static void test_idamax(int n, int inc) {
  int z = n * inc;
  double *x = xmalloc((size_t)z * sizeof(double));
  size_t ir;
  size_t ig;

  fill_f64(x, z, 23);
  ir = ref_l1_iamax(n, x, inc, 1, 1);
  ig = (size_t)cblas_idamax(n, x, inc);

  expect_eq_size("idamax", ig, ir);

  free(x);
}

/* ------------------------------------------------------------------------- */
/* Complex single-precision Level-1 (interleaved re,im)                      */
/* ------------------------------------------------------------------------- */

/* cswap: swap x and y */
static void test_cswap(int n, int inc) {
  int z = 2 * n * inc;
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float *xr = xmalloc((size_t)z * sizeof(float));
  float *yr = xmalloc((size_t)z * sizeof(float));

  fill_c32(x, n * inc, 30);
  fill_c32(y, n * inc, 31);
  memcpy(xr, x, (size_t)z * sizeof(float));
  memcpy(yr, y, (size_t)z * sizeof(float));

  ref_l1_swap(n, xr, inc, yr, inc, 2, 0);
  cblas_cswap(n, x, inc, y, inc);

  expect_close_f32("cswap", x, xr, z, TOL_S_L1);
  expect_close_f32("cswap-y", y, yr, z, TOL_S_L1);

  free(x);
  free(y);
  free(xr);
  free(yr);
}

/* ccopy: y := x */
static void test_ccopy(int n, int inc) {
  int z = 2 * n * inc;
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float *yr = xmalloc((size_t)z * sizeof(float));

  fill_c32(x, n * inc, 30);
  fill_c32(y, n * inc, 32);
  memcpy(yr, y, (size_t)z * sizeof(float));

  ref_l1_copy(n, x, inc, yr, inc, 2, 0);
  cblas_ccopy(n, x, inc, y, inc);

  expect_close_f32("ccopy", y, yr, z, TOL_S_L1);

  free(x);
  free(y);
  free(yr);
}

/* cscal: x := alpha * x  (complex alpha) */
static void test_cscal(int n, int inc) {
  int z = 2 * n * inc;
  float *x = xmalloc((size_t)z * sizeof(float));
  float *xr = xmalloc((size_t)z * sizeof(float));
  float a[2] = {-0.7f, 0.2f};

  fill_c32(x, n * inc, 30);
  memcpy(xr, x, (size_t)z * sizeof(float));

  ref_l1_scal(n, a, xr, inc, 2, 0, 0);
  cblas_cscal(n, a, x, inc);

  expect_close_f32("cscal", x, xr, z, TOL_S_L1);

  free(x);
  free(xr);
}

/* csscal: x := alpha * x  (real alpha) */
static void test_csscal(int n, int inc) {
  int z = 2 * n * inc;
  float *x = xmalloc((size_t)z * sizeof(float));
  float *xr = xmalloc((size_t)z * sizeof(float));
  float a[2] = {-0.7f, 0.2f};

  fill_c32(x, n * inc, 32);
  memcpy(xr, x, (size_t)z * sizeof(float));

  ref_l1_scal(n, a, xr, inc, 2, 0, 1);
  cblas_csscal(n, a[0], x, inc);

  expect_close_f32("csscal", x, xr, z, TOL_S_L1);

  free(x);
  free(xr);
}

/* caxpy: y := alpha * x + y */
static void test_caxpy(int n, int inc) {
  int z = 2 * n * inc;
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float *yr = xmalloc((size_t)z * sizeof(float));
  float a[2] = {-0.7f, 0.2f};

  fill_c32(x, n * inc, 30);
  fill_c32(y, n * inc, 33);
  memcpy(yr, y, (size_t)z * sizeof(float));

  ref_l1_axpy(n, a, x, inc, yr, inc, 2, 0);
  cblas_caxpy(n, a, x, inc, y, inc);

  expect_close_f32("caxpy", y, yr, z, TOL_S_L1);

  free(x);
  free(y);
  free(yr);
}

/* cdotu_sub: unconjugated complex dot product */
static void test_cdotu_sub(int n, int inc) {
  int z = 2 * n * inc;
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float ru[2] = {0.0f, 0.0f};
  float gu[2] = {0.0f, 0.0f};

  fill_c32(x, n * inc, 30);
  fill_c32(y, n * inc, 33);

  ref_l1_dot(n, x, inc, y, inc, ru, 2, 0, 0);
  cblas_cdotu_sub(n, x, inc, y, inc, gu);

  expect_close_f32("cdotu", gu, ru, 2, TOL_S_L1 * (n + 1));

  free(x);
  free(y);
}

/* cdotc_sub: conjugated complex dot product */
static void test_cdotc_sub(int n, int inc) {
  int z = 2 * n * inc;
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float rc[2] = {0.0f, 0.0f};
  float gc[2] = {0.0f, 0.0f};

  fill_c32(x, n * inc, 30);
  fill_c32(y, n * inc, 33);

  ref_l1_dot(n, x, inc, y, inc, rc, 2, 0, 1);
  cblas_cdotc_sub(n, x, inc, y, inc, gc);

  expect_close_f32("cdotc", gc, rc, 2, TOL_S_L1 * (n + 1));

  free(x);
  free(y);
}

/* scnrm2: Euclidean norm of a complex vector */
static void test_scnrm2(int n, int inc) {
  int z = 2 * n * inc;
  float *x = xmalloc((size_t)z * sizeof(float));
  float nr;
  float ng;

  fill_c32(x, n * inc, 30);
  nr = (float)ref_l1_nrm2(n, x, inc, 2, 0);
  ng = cblas_scnrm2(n, x, inc);

  expect_close_f32("scnrm2", &ng, &nr, 1, TOL_S_L1 * (n + 1));

  free(x);
}

/* scasum: sum of |re| + |im| over a complex vector */
static void test_scasum(int n, int inc) {
  int z = 2 * n * inc;
  float *x = xmalloc((size_t)z * sizeof(float));
  float sr;
  float sg;

  fill_c32(x, n * inc, 30);
  sr = (float)ref_l1_asum(n, x, inc, 2, 0);
  sg = cblas_scasum(n, x, inc);

  expect_close_f32("scasum", &sg, &sr, 1, TOL_S_L1 * (n + 1));

  free(x);
}

/* icamax: 1-based index of max |z_i| */
static void test_icamax(int n, int inc) {
  int z = 2 * n * inc;
  float *x = xmalloc((size_t)z * sizeof(float));
  size_t ir;
  size_t ig;

  fill_c32(x, n * inc, 30);
  ir = ref_l1_iamax(n, x, inc, 2, 0);
  ig = (size_t)cblas_icamax(n, x, inc);

  expect_eq_size("icamax", ig, ir);

  free(x);
}

/* ------------------------------------------------------------------------- */
/* Complex double-precision Level-1                                          */
/* ------------------------------------------------------------------------- */

/* zswap: swap x and y */
static void test_zswap(int n, int inc) {
  int z = 2 * n * inc;
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double *xr = xmalloc((size_t)z * sizeof(double));
  double *yr = xmalloc((size_t)z * sizeof(double));

  fill_c64(x, n * inc, 40);
  fill_c64(y, n * inc, 41);
  memcpy(xr, x, (size_t)z * sizeof(double));
  memcpy(yr, y, (size_t)z * sizeof(double));

  ref_l1_swap(n, xr, inc, yr, inc, 2, 1);
  cblas_zswap(n, x, inc, y, inc);

  expect_close_f64("zswap", x, xr, z, TOL_D_L1);
  expect_close_f64("zswap-y", y, yr, z, TOL_D_L1);

  free(x);
  free(y);
  free(xr);
  free(yr);
}

/* zcopy: y := x */
static void test_zcopy(int n, int inc) {
  int z = 2 * n * inc;
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double *yr = xmalloc((size_t)z * sizeof(double));

  fill_c64(x, n * inc, 40);
  fill_c64(y, n * inc, 41);
  memcpy(yr, y, (size_t)z * sizeof(double));

  ref_l1_copy(n, x, inc, yr, inc, 2, 1);
  cblas_zcopy(n, x, inc, y, inc);

  expect_close_f64("zcopy", y, yr, z, TOL_D_L1);

  free(x);
  free(y);
  free(yr);
}

/* zscal: x := alpha * x  (complex alpha) */
static void test_zscal(int n, int inc) {
  int z = 2 * n * inc;
  double *x = xmalloc((size_t)z * sizeof(double));
  double *xr = xmalloc((size_t)z * sizeof(double));
  double a[2] = {-0.7, 0.2};

  fill_c64(x, n * inc, 40);
  memcpy(xr, x, (size_t)z * sizeof(double));

  ref_l1_scal(n, a, xr, inc, 2, 1, 0);
  cblas_zscal(n, a, x, inc);

  expect_close_f64("zscal", x, xr, z, TOL_D_L1);

  free(x);
  free(xr);
}

/* zdscal: x := alpha * x  (real alpha) */
static void test_zdscal(int n, int inc) {
  int z = 2 * n * inc;
  double *x = xmalloc((size_t)z * sizeof(double));
  double *xr = xmalloc((size_t)z * sizeof(double));
  double a[2] = {-0.7, 0.2};

  fill_c64(x, n * inc, 42);
  memcpy(xr, x, (size_t)z * sizeof(double));

  ref_l1_scal(n, a, xr, inc, 2, 1, 1);
  cblas_zdscal(n, a[0], x, inc);

  expect_close_f64("zdscal", x, xr, z, TOL_D_L1);

  free(x);
  free(xr);
}

/* zaxpy: y := alpha * x + y */
static void test_zaxpy(int n, int inc) {
  int z = 2 * n * inc;
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double *yr = xmalloc((size_t)z * sizeof(double));
  double a[2] = {-0.7, 0.2};

  fill_c64(x, n * inc, 40);
  fill_c64(y, n * inc, 43);
  memcpy(yr, y, (size_t)z * sizeof(double));

  ref_l1_axpy(n, a, x, inc, yr, inc, 2, 1);
  cblas_zaxpy(n, a, x, inc, y, inc);

  expect_close_f64("zaxpy", y, yr, z, TOL_D_L1);

  free(x);
  free(y);
  free(yr);
}

/* zdotu_sub: unconjugated complex dot product */
static void test_zdotu_sub(int n, int inc) {
  int z = 2 * n * inc;
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double ru[2] = {0.0, 0.0};
  double gu[2] = {0.0, 0.0};

  fill_c64(x, n * inc, 40);
  fill_c64(y, n * inc, 43);

  ref_l1_dot(n, x, inc, y, inc, ru, 2, 1, 0);
  cblas_zdotu_sub(n, x, inc, y, inc, gu);

  expect_close_f64("zdotu", gu, ru, 2, TOL_D_L1 * (n + 1));

  free(x);
  free(y);
}

/* zdotc_sub: conjugated complex dot product */
static void test_zdotc_sub(int n, int inc) {
  int z = 2 * n * inc;
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double rc[2] = {0.0, 0.0};
  double gc[2] = {0.0, 0.0};

  fill_c64(x, n * inc, 40);
  fill_c64(y, n * inc, 43);

  ref_l1_dot(n, x, inc, y, inc, rc, 2, 1, 1);
  cblas_zdotc_sub(n, x, inc, y, inc, gc);

  expect_close_f64("zdotc", gc, rc, 2, TOL_D_L1 * (n + 1));

  free(x);
  free(y);
}

/* dznrm2: Euclidean norm of a complex vector */
static void test_dznrm2(int n, int inc) {
  int z = 2 * n * inc;
  double *x = xmalloc((size_t)z * sizeof(double));
  double nr;
  double ng;

  fill_c64(x, n * inc, 40);
  nr = ref_l1_nrm2(n, x, inc, 2, 1);
  ng = cblas_dznrm2(n, x, inc);

  expect_close_f64("dznrm2", &ng, &nr, 1, TOL_D_L1 * (n + 1));

  free(x);
}

/* dzasum: sum of |re| + |im| over a complex vector */
static void test_dzasum(int n, int inc) {
  int z = 2 * n * inc;
  double *x = xmalloc((size_t)z * sizeof(double));
  double sr;
  double sg;

  fill_c64(x, n * inc, 40);
  sr = ref_l1_asum(n, x, inc, 2, 1);
  sg = cblas_dzasum(n, x, inc);

  expect_close_f64("dzasum", &sg, &sr, 1, TOL_D_L1 * (n + 1));

  free(x);
}

/* izamax: 1-based index of max |z_i| */
static void test_izamax(int n, int inc) {
  int z = 2 * n * inc;
  double *x = xmalloc((size_t)z * sizeof(double));
  size_t ir;
  size_t ig;

  fill_c64(x, n * inc, 40);
  ir = ref_l1_iamax(n, x, inc, 2, 1);
  ig = (size_t)cblas_izamax(n, x, inc);

  expect_eq_size("izamax", ig, ir);

  free(x);
}

/* ------------------------------------------------------------------------- */
/* Rotations                                                                 */
/* ------------------------------------------------------------------------- */

/* srot: plane rotation of (x, y) */
static void test_srot(int n) {
  int k = n > 32 ? 32 : n;
  float x[64], y[64], xr[64], yr[64];

  fill_f32(x, k, 50);
  fill_f32(y, k, 51);
  memcpy(xr, x, (size_t)k * sizeof(float));
  memcpy(yr, y, (size_t)k * sizeof(float));

  ref_l1_rot(k, xr, 1, yr, 1, 0.8, 0.6, 1, 0);
  cblas_srot(k, x, 1, y, 1, 0.8f, 0.6f);

  expect_close_f32("srot", x, xr, k, TOL_S_L1);
  expect_close_f32("srot-y", y, yr, k, TOL_S_L1);
}

/* drot: plane rotation of (x, y) */
static void test_drot(int n) {
  int k = n > 32 ? 32 : n;
  double x[64], y[64], xr[64], yr[64];

  fill_f64(x, k, 50);
  fill_f64(y, k, 51);
  memcpy(xr, x, (size_t)k * sizeof(double));
  memcpy(yr, y, (size_t)k * sizeof(double));

  ref_l1_rot(k, xr, 1, yr, 1, 0.8, 0.6, 1, 1);
  cblas_drot(k, x, 1, y, 1, 0.8, 0.6);

  expect_close_f64("drot", x, xr, k, TOL_D_L1);
  expect_close_f64("drot-y", y, yr, k, TOL_D_L1);
}

/* srotm: modified Givens rotation of (x, y) */
static void test_srotm(int n) {
  int k = n > 32 ? 32 : n;
  float x[64], y[64], xr[64], yr[64];
  float p[5] = {-1.0f, 0.8f, -0.2f, 0.3f, 1.1f};

  fill_f32(x, k, 52);
  fill_f32(y, k, 53);
  memcpy(xr, x, (size_t)k * sizeof(float));
  memcpy(yr, y, (size_t)k * sizeof(float));

  ref_l1_rotm(k, xr, 1, yr, 1, p, 0);
  cblas_srotm(k, x, 1, y, 1, p);

  expect_close_f32("srotm", x, xr, k, TOL_S_L1);
  expect_close_f32("srotm-y", y, yr, k, TOL_S_L1);
}

/* drotm: modified Givens rotation of (x, y) */
static void test_drotm(int n) {
  int k = n > 32 ? 32 : n;
  double x[64], y[64], xr[64], yr[64];
  double p[5] = {-1.0, 0.8, -0.2, 0.3, 1.1};

  fill_f64(x, k, 52);
  fill_f64(y, k, 53);
  memcpy(xr, x, (size_t)k * sizeof(double));
  memcpy(yr, y, (size_t)k * sizeof(double));

  ref_l1_rotm(k, xr, 1, yr, 1, p, 1);
  cblas_drotm(k, x, 1, y, 1, p);

  expect_close_f64("drotm", x, xr, k, TOL_D_L1);
  expect_close_f64("drotm-y", y, yr, k, TOL_D_L1);
}

/* csrot: real plane rotation applied to complex vectors */
static void test_csrot(int n) {
  int k = n > 32 ? 32 : n;
  float x[128], y[128], xr[128], yr[128];

  fill_c32(x, k, 54);
  fill_c32(y, k, 55);
  memcpy(xr, x, (size_t)2 * k * sizeof(float));
  memcpy(yr, y, (size_t)2 * k * sizeof(float));

  ref_l1_rot(k, xr, 1, yr, 1, 0.8, 0.6, 2, 0);
  cblas_csrot(k, x, 1, y, 1, 0.8f, 0.6f);

  expect_close_f32("csrot", x, xr, 2 * k, TOL_S_L1);
  expect_close_f32("csrot-y", y, yr, 2 * k, TOL_S_L1);
}

/* zdrot: real plane rotation applied to complex vectors */
static void test_zdrot(int n) {
  int k = n > 32 ? 32 : n;
  double x[128], y[128], xr[128], yr[128];

  fill_c64(x, k, 54);
  fill_c64(y, k, 55);
  memcpy(xr, x, (size_t)2 * k * sizeof(double));
  memcpy(yr, y, (size_t)2 * k * sizeof(double));

  ref_l1_rot(k, xr, 1, yr, 1, 0.8, 0.6, 2, 1);
  cblas_zdrot(k, x, 1, y, 1, 0.8, 0.6);

  expect_close_f64("zdrot", x, xr, 2 * k, TOL_D_L1);
  expect_close_f64("zdrot-y", y, yr, 2 * k, TOL_D_L1);
}

/* srotg: generate a real Givens rotation (residual check) */
static void test_srotg(void) {
  float a = 3.0f, b = 4.0f, c = 0.0f, s = 0.0f;
  float e;
  float zero = 0.0f;

  cblas_srotg(&a, &b, &c, &s);
  e = (c * 3.0f + s * 4.0f) - a;
  expect_close_f32("srotg", &e, &zero, 1, TOL_S_L1 * 8);
}

/* drotg: generate a real Givens rotation (residual check) */
static void test_drotg(void) {
  double a = 3.0, b = 4.0, c = 0.0, s = 0.0;
  double e;
  double zero = 0.0;

  cblas_drotg(&a, &b, &c, &s);
  e = (c * 3.0 + s * 4.0) - a;
  expect_close_f64("drotg", &e, &zero, 1, TOL_D_L1 * 8);
}

/* crotg: generate a complex Givens rotation (residual check) */
static void test_crotg(void) {
  float a[2] = {3.0f, 1.0f}, b[2] = {2.0f, -1.0f}, s[2] = {0.0f, 0.0f};
  float c = 0.0f;
  float e0, e1, e;
  float zero = 0.0f;

  cblas_crotg(a, b, &c, s);
  e0 = -(s[0] * 3.0f + s[1] * 1.0f) + c * 2.0f;
  e1 = -(-s[1] * 3.0f + s[0] * 1.0f) - c;
  e = fabsf(e0) + fabsf(e1);
  expect_close_f32("crotg", &e, &zero, 1, TOL_S_L1 * 16);
}

/* zrotg: generate a complex Givens rotation (residual check) */
static void test_zrotg(void) {
  double a[2] = {3.0, 1.0}, b[2] = {2.0, -1.0}, s[2] = {0.0, 0.0};
  double c = 0.0;
  double e0, e1, e;
  double zero = 0.0;

  cblas_zrotg(a, b, &c, s);
  e0 = -(s[0] * 3.0 + s[1] * 1.0) + c * 2.0;
  e1 = -(-s[1] * 3.0 + s[0] * 1.0) - c;
  e = fabs(e0) + fabs(e1);
  expect_close_f64("zrotg", &e, &zero, 1, TOL_D_L1 * 16);
}

/* srotmg: generate a modified Givens transformation */
static void test_srotmg(void) {
  float d1 = 1.0f, d2 = 2.0f, b1 = 3.0f, b2 = 4.0f, p[5] = {9, 9, 9, 9, 9};
  float rd1 = 1.0f, rd2 = 2.0f, rb1 = 3.0f, rp[5] = {9, 9, 9, 9, 9};
  float got[8], ref[8];

  ref_l1_rotmg(&rd1, &rd2, &rb1, &b2, rp, 0);
  cblas_srotmg(&d1, &d2, &b1, b2, p);

  got[0] = d1;
  got[1] = d2;
  got[2] = b1;
  got[3] = p[0];
  got[4] = p[1];
  got[5] = p[2];
  got[6] = p[3];
  got[7] = p[4];
  ref[0] = rd1;
  ref[1] = rd2;
  ref[2] = rb1;
  ref[3] = rp[0];
  ref[4] = rp[1];
  ref[5] = rp[2];
  ref[6] = rp[3];
  ref[7] = rp[4];

  expect_close_f32("srotmg", got, ref, 8, TOL_S_L1 * 16);
}

/* drotmg: generate a modified Givens transformation */
static void test_drotmg(void) {
  double d1 = 1.0, d2 = 2.0, b1 = 3.0, b2 = 4.0, p[5] = {9, 9, 9, 9, 9};
  double rd1 = 1.0, rd2 = 2.0, rb1 = 3.0, rp[5] = {9, 9, 9, 9, 9};
  double got[8], ref[8];

  ref_l1_rotmg(&rd1, &rd2, &rb1, &b2, rp, 1);
  cblas_drotmg(&d1, &d2, &b1, b2, p);

  got[0] = d1;
  got[1] = d2;
  got[2] = b1;
  got[3] = p[0];
  got[4] = p[1];
  got[5] = p[2];
  got[6] = p[3];
  got[7] = p[4];
  ref[0] = rd1;
  ref[1] = rd2;
  ref[2] = rb1;
  ref[3] = rp[0];
  ref[4] = rp[1];
  ref[5] = rp[2];
  ref[6] = rp[3];
  ref[7] = rp[4];

  expect_close_f64("drotmg", got, ref, 8, TOL_D_L1 * 16);
}

/* ------------------------------------------------------------------------- */
/* Drivers                                                                   */
/* ------------------------------------------------------------------------- */

static void run_deep_axpy(void) {
  printf("==> L1 AXPY (deep remainder / stride grid)\n");
  for (int s = 0; s < NS_L1; s++) {
    int n = SIZES_L1[s];
    test_saxpy(n, 1, 1);
    test_daxpy(n, 1, 1);
    for (int i = 0; i < NINCS; i++) {
      int inc = INCS[i];
      if (inc == 1)
        continue;
      test_saxpy(n, inc, inc);
      test_daxpy(n, inc, inc);
    }
  }
}

static void run_compact_l1(void) {
  printf("==> L1 compact S/D/C/Z + rotations\n");

  test_srotg();
  test_drotg();
  test_crotg();
  test_zrotg();
  test_srotmg();
  test_drotmg();

  for (int i = 0; i < NS_FULL; i++) {
    int n = SIZES_FULL[i];
    int incs[2] = {1, 2};

    for (int j = 0; j < 2; j++) {
      int inc = incs[j];

      test_sswap(n, inc);
      test_scopy(n, inc);
      test_sscal(n, inc);
      test_saxpy_compact(n, inc);
      test_sdot(n, inc);
      test_snrm2(n, inc);
      test_sasum(n, inc);
      test_isamax(n, inc);

      test_dswap(n, inc);
      test_dcopy(n, inc);
      test_dscal(n, inc);
      test_daxpy_compact(n, inc);
      test_ddot(n, inc);
      test_dnrm2(n, inc);
      test_dasum(n, inc);
      test_idamax(n, inc);

      test_cswap(n, inc);
      test_ccopy(n, inc);
      test_cscal(n, inc);
      test_csscal(n, inc);
      test_caxpy(n, inc);
      test_cdotu_sub(n, inc);
      test_cdotc_sub(n, inc);
      test_scnrm2(n, inc);
      test_scasum(n, inc);
      test_icamax(n, inc);

      test_zswap(n, inc);
      test_zcopy(n, inc);
      test_zscal(n, inc);
      test_zdscal(n, inc);
      test_zaxpy(n, inc);
      test_zdotu_sub(n, inc);
      test_zdotc_sub(n, inc);
      test_dznrm2(n, inc);
      test_dzasum(n, inc);
      test_izamax(n, inc);
    }

    test_srot(n);
    test_drot(n);
    test_srotm(n);
    test_drotm(n);
    test_csrot(n);
    test_zdrot(n);
  }
}

void check_l1(void) {
  run_deep_axpy();
  run_compact_l1();
}
