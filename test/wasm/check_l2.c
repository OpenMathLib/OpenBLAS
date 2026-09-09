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
 * Level-2 CBLAS checks vs the scalar oracle in ref_l2.c / ref.c.
 *
 * A deep SGEMV/DGEMV grid stresses tile remainders, rectangular shapes,
 * transposes and non-unit strides (SIZES_L2). Every other Level-2 op --
 * dense, banded and packed, for s/d/c/z -- runs on the compact SIZES_FULL
 * grid, where the matrices stay small enough to build by hand.
 *
 * One static test_* function per public CBLAS interface: allocate, fill,
 * call the scalar reference, call OpenBLAS, then expect_close_*.
 *
 * Matrices are column major throughout. Complex data is interleaved re,im,
 * so a complex n-vector occupies 2*n floats/doubles and the element (i, j)
 * of an n-by-n matrix starts at index 2*(i + j*n).
 */

#include "cases.h"
#include "common.h"
#include "ref.h"
#include "tol.h"

/* ------------------------------------------------------------------------- */
/* Shared helpers                                                            */
/* ------------------------------------------------------------------------- */

/* Storage length for an n-vector with stride |inc| (inc==0 → length 1). */
static int vec_len(int n, int inc) {
  int abs_inc = inc == 0 ? 1 : (inc > 0 ? inc : -inc);

  if (n <= 0)
    return 1;
  return inc == 0 ? 1 : n * abs_inc;
}

/* Level-2 tolerance grows with n; tag the failure message with the size. */
static void expect_l2_f32(const char *op, int n, const float *got,
                          const float *ref, int z) {
  char msg[96];

  snprintf(msg, sizeof(msg), "%s n=%d", op, n);
  expect_close_f32(msg, got, ref, z, tol_s_l2(n));
}

static void expect_l2_f64(const char *op, int n, const double *got,
                          const double *ref, int z) {
  char msg[96];

  snprintf(msg, sizeof(msg), "%s n=%d", op, n);
  expect_close_f64(msg, got, ref, z, tol_d_l2(n));
}

/* Half-bandwidth used by the banded tests: 3, or n-1 for the tiny sizes. */
static int band_width(int n) { return n > 3 ? 3 : n - 1; }

/* Index of element (i, j), i >= j, in a column-major packed lower triangle. */
static int packed_lower(int i, int j, int n) {
  return i + (2 * n - j - 1) * j / 2;
}

/* ------------------------------------------------------------------------- */
/* Deep GEMV grids (shapes, transposes, strides)                             */
/* ------------------------------------------------------------------------- */

/* sgemv: y := alpha*op(A)*x + beta*y */
static void test_sgemv(enum CBLAS_TRANSPOSE trans, int m, int n, int incx,
                       int incy) {
  int lda = m;
  int lenx = (trans == CblasNoTrans) ? n : m;
  int leny = (trans == CblasNoTrans) ? m : n;
  int nx = vec_len(lenx, incx);
  int ny = vec_len(leny, incy);
  float *A = xmalloc((size_t)lda * (size_t)n * sizeof(float));
  float *x = xmalloc((size_t)nx * sizeof(float));
  float *y = xmalloc((size_t)ny * sizeof(float));
  float *yr = xmalloc((size_t)ny * sizeof(float));
  float alpha = 1.1f;
  float beta = 0.7f;
  char msg[160];

  fill_f32(A, lda * n, 3);
  fill_f32(x, nx, 4);
  fill_f32(y, ny, 5);
  memcpy(yr, y, (size_t)ny * sizeof(float));

  ref_sgemv(trans, m, n, alpha, A, lda, x, incx, beta, yr, incy);
  cblas_sgemv(CblasColMajor, trans, m, n, alpha, A, lda, x, incx, beta, y,
              incy);

  snprintf(msg, sizeof(msg), "sgemv t=%d m=%d n=%d incx=%d incy=%d", (int)trans,
           m, n, incx, incy);
  expect_close_f32(msg, y, yr, ny, tol_s_l2(m > n ? m : n));

  free(A);
  free(x);
  free(y);
  free(yr);
}

/* dgemv: y := alpha*op(A)*x + beta*y */
static void test_dgemv(enum CBLAS_TRANSPOSE trans, int m, int n, int incx,
                       int incy) {
  int lda = m;
  int lenx = (trans == CblasNoTrans) ? n : m;
  int leny = (trans == CblasNoTrans) ? m : n;
  int nx = vec_len(lenx, incx);
  int ny = vec_len(leny, incy);
  double *A = xmalloc((size_t)lda * (size_t)n * sizeof(double));
  double *x = xmalloc((size_t)nx * sizeof(double));
  double *y = xmalloc((size_t)ny * sizeof(double));
  double *yr = xmalloc((size_t)ny * sizeof(double));
  double alpha = 1.1;
  double beta = 0.7;
  char msg[160];

  fill_f64(A, lda * n, 3);
  fill_f64(x, nx, 4);
  fill_f64(y, ny, 5);
  memcpy(yr, y, (size_t)ny * sizeof(double));

  ref_dgemv(trans, m, n, alpha, A, lda, x, incx, beta, yr, incy);
  cblas_dgemv(CblasColMajor, trans, m, n, alpha, A, lda, x, incx, beta, y,
              incy);

  snprintf(msg, sizeof(msg), "dgemv t=%d m=%d n=%d incx=%d incy=%d", (int)trans,
           m, n, incx, incy);
  expect_close_f64(msg, y, yr, ny, tol_d_l2(m > n ? m : n));

  free(A);
  free(x);
  free(y);
  free(yr);
}

/* ------------------------------------------------------------------------- */
/* Real single-precision matrix builders                                     */
/* ------------------------------------------------------------------------- */

/* Mirror the lower triangle of A into the strict upper one. */
static void symmetrize_f32(float *A, int n) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < j; i++)
      A[i + j * n] = A[j + i * n];
}

/* A pseudo-random symmetric matrix. */
static void make_sym_f32(float *A, int n, int seed) {
  fill_f32(A, n * n, seed);
  symmetrize_f32(A, n);
}

/* Drop the strict upper triangle, keeping the lower triangular part. */
static void zero_upper_f32(float *A, int n) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < j; i++)
      A[i + j * n] = 0.0f;
}

/* Copy the lower triangle of A into packed storage. */
static void pack_lower_f32(float *ap, const float *A, int n) {
  for (int j = 0; j < n; j++)
    for (int i = j; i < n; i++)
      ap[packed_lower(i, j, n)] = A[i + j * n];
}

/* A general band matrix, both dense (A) and in (kl=ku=k) band storage (ab),
 * where band element (i, j) lives at ab[k + i - j + j*ldab], ldab = 2*k + 1. */
static void make_gband_f32(float *A, float *ab, int n, int k) {
  int ldab = 2 * k + 1;

  memset(A, 0, (size_t)n * n * sizeof(float));
  memset(ab, 0, (size_t)ldab * n * sizeof(float));
  for (int j = 0; j < n; j++)
    for (int i = j - k; i <= j + k; i++) {
      float v;

      if (i < 0 || i >= n)
        continue;
      v = 0.15f + 0.01f * (float)(i + 2 * j);
      A[i + j * n] = v;
      ab[k + i - j + j * ldab] = v;
    }
}

/* A symmetric band matrix, both dense (A) and in lower band storage (ab),
 * where band element (i, j) lives at ab[i - j + j*ldab], ldab = k + 1.
 * The diagonal is 2, so the lower triangle is safe to invert in TBSV. */
static void make_sband_f32(float *A, float *ab, int n, int k) {
  int ldab = k + 1;

  memset(A, 0, (size_t)n * n * sizeof(float));
  memset(ab, 0, (size_t)ldab * n * sizeof(float));
  for (int j = 0; j < n; j++)
    for (int i = j; i <= j + k && i < n; i++) {
      float v = (i == j) ? 2.0f : 0.1f + 0.01f * (float)(i + j);

      A[i + j * n] = v;
      A[j + i * n] = v;
      ab[i - j + j * ldab] = v;
    }
}

/* Packed lower triangle of alpha*x*x^T, the SPR reference. */
static void ref_spr_f32(float *ap, const float *x, int n, float alpha) {
  for (int j = 0; j < n; j++)
    for (int i = j; i < n; i++)
      ap[packed_lower(i, j, n)] = alpha * x[i] * x[j];
}

/* ------------------------------------------------------------------------- */
/* Real single-precision Level-2 (compact grid)                              */
/* ------------------------------------------------------------------------- */

/* sgemv: y := alpha*A*x + beta*y */
static void test_sgemv_compact(int n) {
  float *A = xmalloc((size_t)n * n * sizeof(float));
  float *x = xmalloc((size_t)n * sizeof(float));
  float *y = xmalloc((size_t)n * sizeof(float));
  float *yr = xmalloc((size_t)n * sizeof(float));
  float alpha = 0.8f;
  float beta = -0.3f;

  make_sym_f32(A, n, 60);
  fill_f32(x, n, 61);
  fill_f32(y, n, 62);
  memcpy(yr, y, (size_t)n * sizeof(float));

  ref_l2_mv(n, A, x, yr, &alpha, &beta, 1, 0, 0);
  cblas_sgemv(CblasColMajor, CblasNoTrans, n, n, alpha, A, n, x, 1, beta, y, 1);

  expect_l2_f32("sgemv", n, y, yr, n);

  free(A);
  free(x);
  free(y);
  free(yr);
}

/* ssymv: y := alpha*A*x + beta*y, symmetric A read from its lower triangle */
static void test_ssymv(int n) {
  float *A = xmalloc((size_t)n * n * sizeof(float));
  float *x = xmalloc((size_t)n * sizeof(float));
  float *y = xmalloc((size_t)n * sizeof(float));
  float *yr = xmalloc((size_t)n * sizeof(float));
  float alpha = 0.8f;
  float beta = -0.3f;

  make_sym_f32(A, n, 60);
  fill_f32(x, n, 61);
  fill_f32(y, n, 63);
  memcpy(yr, y, (size_t)n * sizeof(float));

  ref_l2_mv(n, A, x, yr, &alpha, &beta, 1, 0, 0);
  cblas_ssymv(CblasColMajor, CblasLower, n, alpha, A, n, x, 1, beta, y, 1);

  expect_l2_f32("ssymv", n, y, yr, n);

  free(A);
  free(x);
  free(y);
  free(yr);
}

/* strmv: x := A*x, A lower triangular with a non-unit diagonal */
static void test_strmv(int n) {
  float *A = xmalloc((size_t)n * n * sizeof(float));
  float *x = xmalloc((size_t)n * sizeof(float));
  float *y = xmalloc((size_t)n * sizeof(float));
  float *yr = xmalloc((size_t)n * sizeof(float));
  float one = 1.0f;
  float zero = 0.0f;

  make_tri_f32(A, n, n, CblasLower, 0);
  fill_f32(x, n, 61);
  memcpy(y, x, (size_t)n * sizeof(float));
  memset(yr, 0, (size_t)n * sizeof(float));

  ref_l2_mv(n, A, x, yr, &one, &zero, 1, 0, 0);
  cblas_strmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, A, n, y,
              1);

  expect_l2_f32("strmv", n, y, yr, n);

  free(A);
  free(x);
  free(y);
  free(yr);
}

/* strsv: solve A*y = A*x on the same triangle, so y must come back as x */
static void test_strsv(int n) {
  float *A = xmalloc((size_t)n * n * sizeof(float));
  float *x = xmalloc((size_t)n * sizeof(float));
  float *y = xmalloc((size_t)n * sizeof(float));
  float one = 1.0f;
  float zero = 0.0f;

  make_tri_f32(A, n, n, CblasLower, 0);
  fill_f32(x, n, 61);
  memset(y, 0, (size_t)n * sizeof(float));

  ref_l2_mv(n, A, x, y, &one, &zero, 1, 0, 0);
  cblas_strsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, A, n, y,
              1);

  expect_l2_f32("strsv", n, y, x, n);

  free(A);
  free(x);
  free(y);
}

/* strmv: x := A^T*x, A upper triangular with an implicit unit diagonal.
 * The reference takes the transpose through its "conjugate" flag, which for
 * real data is a plain transpose. */
static void test_strmv_upper_unit(int n) {
  float *A = xmalloc((size_t)n * n * sizeof(float));
  float *x = xmalloc((size_t)n * sizeof(float));
  float *y = xmalloc((size_t)n * sizeof(float));
  float *yr = xmalloc((size_t)n * sizeof(float));
  float one = 1.0f;
  float zero = 0.0f;

  make_tri_f32(A, n, n, CblasUpper, 1);
  fill_f32(x, n, 61);
  memcpy(y, x, (size_t)n * sizeof(float));
  memset(yr, 0, (size_t)n * sizeof(float));

  ref_l2_mv(n, A, x, yr, &one, &zero, 1, 0, 1);
  cblas_strmv(CblasColMajor, CblasUpper, CblasTrans, CblasUnit, n, A, n, y, 1);

  expect_l2_f32("strmv-upper-unit", n, y, yr, n);

  free(A);
  free(x);
  free(y);
  free(yr);
}

/* strsv: solve A^T*y = A^T*x on the upper unit triangle */
static void test_strsv_upper_unit(int n) {
  float *A = xmalloc((size_t)n * n * sizeof(float));
  float *x = xmalloc((size_t)n * sizeof(float));
  float *y = xmalloc((size_t)n * sizeof(float));
  float one = 1.0f;
  float zero = 0.0f;

  make_tri_f32(A, n, n, CblasUpper, 1);
  fill_f32(x, n, 61);
  memset(y, 0, (size_t)n * sizeof(float));

  ref_l2_mv(n, A, x, y, &one, &zero, 1, 0, 1);
  cblas_strsv(CblasColMajor, CblasUpper, CblasTrans, CblasUnit, n, A, n, y, 1);

  expect_l2_f32("strsv-upper-unit", n, y, x, n);

  free(A);
  free(x);
  free(y);
}

/* sger: A := alpha*x*y^T + A, starting from a zero A */
static void test_sger(int n) {
  int nn = n * n;
  float *A = xmalloc((size_t)nn * sizeof(float));
  float *Ar = xmalloc((size_t)nn * sizeof(float));
  float *x = xmalloc((size_t)n * sizeof(float));
  float *y = xmalloc((size_t)n * sizeof(float));
  float alpha = 0.8f;

  fill_f32(x, n, 61);
  fill_f32(y, n, 68);
  memset(A, 0, (size_t)nn * sizeof(float));
  memset(Ar, 0, (size_t)nn * sizeof(float));

  ref_l2_rank(n, Ar, x, y, &alpha, 1, 0, 0, 0, 0);
  cblas_sger(CblasColMajor, n, n, alpha, x, 1, y, 1, A, n);

  expect_l2_f32("sger", n, A, Ar, nn);

  free(A);
  free(Ar);
  free(x);
  free(y);
}

/* ssyr: A := alpha*x*x^T + A, lower triangle only */
static void test_ssyr(int n) {
  int nn = n * n;
  float *A = xmalloc((size_t)nn * sizeof(float));
  float *Ar = xmalloc((size_t)nn * sizeof(float));
  float *x = xmalloc((size_t)n * sizeof(float));
  float alpha = 0.8f;

  fill_f32(x, n, 61);
  memset(A, 0, (size_t)nn * sizeof(float));
  memset(Ar, 0, (size_t)nn * sizeof(float));

  ref_l2_rank(n, Ar, x, x, &alpha, 1, 0, 0, 1, 0);
  cblas_ssyr(CblasColMajor, CblasLower, n, alpha, x, 1, A, n);

  expect_l2_f32("ssyr", n, A, Ar, nn);

  free(A);
  free(Ar);
  free(x);
}

/* ssyr2: A := alpha*(x*y^T + y*x^T) + A, lower triangle only.
 * A starts from a rank-1 update so the update is checked as an accumulation. */
static void test_ssyr2(int n) {
  int nn = n * n;
  float *A = xmalloc((size_t)nn * sizeof(float));
  float *Ar = xmalloc((size_t)nn * sizeof(float));
  float *x = xmalloc((size_t)n * sizeof(float));
  float *y = xmalloc((size_t)n * sizeof(float));
  float alpha = 0.8f;

  fill_f32(x, n, 61);
  fill_f32(y, n, 68);
  memset(Ar, 0, (size_t)nn * sizeof(float));
  ref_l2_rank(n, Ar, x, x, &alpha, 1, 0, 0, 1, 0);
  memcpy(A, Ar, (size_t)nn * sizeof(float));

  ref_l2_rank(n, Ar, x, y, &alpha, 1, 0, 0, 1, 0);
  ref_l2_rank(n, Ar, y, x, &alpha, 1, 0, 0, 1, 0);
  cblas_ssyr2(CblasColMajor, CblasLower, n, alpha, x, 1, y, 1, A, n);

  expect_l2_f32("ssyr2", n, A, Ar, nn);

  free(A);
  free(Ar);
  free(x);
  free(y);
}

/* sgbmv: y := alpha*A*x + beta*y with A in general band storage */
static void test_sgbmv(int n) {
  int k = band_width(n);
  int ldab = 2 * k + 1;
  float *A = xmalloc((size_t)n * n * sizeof(float));
  float *ab = xmalloc((size_t)ldab * n * sizeof(float));
  float *x = xmalloc((size_t)n * sizeof(float));
  float *y = xmalloc((size_t)n * sizeof(float));
  float *yr = xmalloc((size_t)n * sizeof(float));
  float alpha = 0.8f;
  float beta = -0.3f;

  make_gband_f32(A, ab, n, k);
  fill_f32(x, n, 61);
  fill_f32(y, n, 64);
  memcpy(yr, y, (size_t)n * sizeof(float));

  ref_l2_mv(n, A, x, yr, &alpha, &beta, 1, 0, 0);
  cblas_sgbmv(CblasColMajor, CblasNoTrans, n, n, k, k, alpha, ab, ldab, x, 1,
              beta, y, 1);

  expect_l2_f32("sgbmv", n, y, yr, n);

  free(A);
  free(ab);
  free(x);
  free(y);
  free(yr);
}

/* ssbmv: y := alpha*A*x + beta*y with A in symmetric band storage */
static void test_ssbmv(int n) {
  int k = band_width(n);
  int ldab = k + 1;
  float *A = xmalloc((size_t)n * n * sizeof(float));
  float *ab = xmalloc((size_t)ldab * n * sizeof(float));
  float *x = xmalloc((size_t)n * sizeof(float));
  float *y = xmalloc((size_t)n * sizeof(float));
  float *yr = xmalloc((size_t)n * sizeof(float));
  float alpha = 0.8f;
  float beta = -0.3f;

  make_sband_f32(A, ab, n, k);
  fill_f32(x, n, 61);
  fill_f32(y, n, 65);
  memcpy(yr, y, (size_t)n * sizeof(float));

  ref_l2_mv(n, A, x, yr, &alpha, &beta, 1, 0, 0);
  cblas_ssbmv(CblasColMajor, CblasLower, n, k, alpha, ab, ldab, x, 1, beta, y,
              1);

  expect_l2_f32("ssbmv", n, y, yr, n);

  free(A);
  free(ab);
  free(x);
  free(y);
  free(yr);
}

/* stbmv: x := A*x with A the lower triangle of a band matrix */
static void test_stbmv(int n) {
  int k = band_width(n);
  int ldab = k + 1;
  float *A = xmalloc((size_t)n * n * sizeof(float));
  float *ab = xmalloc((size_t)ldab * n * sizeof(float));
  float *x = xmalloc((size_t)n * sizeof(float));
  float *y = xmalloc((size_t)n * sizeof(float));
  float *yr = xmalloc((size_t)n * sizeof(float));
  float one = 1.0f;
  float zero = 0.0f;

  make_sband_f32(A, ab, n, k);
  zero_upper_f32(A, n);
  fill_f32(x, n, 61);
  memcpy(y, x, (size_t)n * sizeof(float));
  memset(yr, 0, (size_t)n * sizeof(float));

  ref_l2_mv(n, A, x, yr, &one, &zero, 1, 0, 0);
  cblas_stbmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, k, ab,
              ldab, y, 1);

  expect_l2_f32("stbmv", n, y, yr, n);

  free(A);
  free(ab);
  free(x);
  free(y);
  free(yr);
}

/* stbsv: solve A*y = A*x on the same triangular band */
static void test_stbsv(int n) {
  int k = band_width(n);
  int ldab = k + 1;
  float *A = xmalloc((size_t)n * n * sizeof(float));
  float *ab = xmalloc((size_t)ldab * n * sizeof(float));
  float *x = xmalloc((size_t)n * sizeof(float));
  float *y = xmalloc((size_t)n * sizeof(float));
  float one = 1.0f;
  float zero = 0.0f;

  make_sband_f32(A, ab, n, k);
  zero_upper_f32(A, n);
  fill_f32(x, n, 61);
  memset(y, 0, (size_t)n * sizeof(float));

  ref_l2_mv(n, A, x, y, &one, &zero, 1, 0, 0);
  cblas_stbsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, k, ab,
              ldab, y, 1);

  expect_l2_f32("stbsv", n, y, x, n);

  free(A);
  free(ab);
  free(x);
  free(y);
}

/* sspmv: y := alpha*A*x + beta*y with symmetric A in packed storage */
static void test_sspmv(int n) {
  int np = n * (n + 1) / 2;
  float *A = xmalloc((size_t)n * n * sizeof(float));
  float *ap = xmalloc((size_t)np * sizeof(float));
  float *x = xmalloc((size_t)n * sizeof(float));
  float *y = xmalloc((size_t)n * sizeof(float));
  float *yr = xmalloc((size_t)n * sizeof(float));
  float alpha = 0.8f;
  float beta = -0.3f;

  /* Pack the lower triangle, then mirror it so the dense reference sees the
   * same symmetric matrix. */
  make_tri_f32(A, n, n, CblasLower, 0);
  pack_lower_f32(ap, A, n);
  symmetrize_f32(A, n);

  fill_f32(x, n, 61);
  fill_f32(y, n, 66);
  memcpy(yr, y, (size_t)n * sizeof(float));

  ref_l2_mv(n, A, x, yr, &alpha, &beta, 1, 0, 0);
  cblas_sspmv(CblasColMajor, CblasLower, n, alpha, ap, x, 1, beta, y, 1);

  expect_l2_f32("sspmv", n, y, yr, n);

  free(A);
  free(ap);
  free(x);
  free(y);
  free(yr);
}

/* stpmv: x := A*x with A lower triangular in packed storage */
static void test_stpmv(int n) {
  int np = n * (n + 1) / 2;
  float *A = xmalloc((size_t)n * n * sizeof(float));
  float *ap = xmalloc((size_t)np * sizeof(float));
  float *x = xmalloc((size_t)n * sizeof(float));
  float *y = xmalloc((size_t)n * sizeof(float));
  float *yr = xmalloc((size_t)n * sizeof(float));
  float one = 1.0f;
  float zero = 0.0f;

  make_tri_f32(A, n, n, CblasLower, 0);
  pack_lower_f32(ap, A, n);
  fill_f32(x, n, 61);
  memcpy(y, x, (size_t)n * sizeof(float));
  memset(yr, 0, (size_t)n * sizeof(float));

  ref_l2_mv(n, A, x, yr, &one, &zero, 1, 0, 0);
  cblas_stpmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, ap, y,
              1);

  expect_l2_f32("stpmv", n, y, yr, n);

  free(A);
  free(ap);
  free(x);
  free(y);
  free(yr);
}

/* stpsv: solve A*y = A*x on the same packed triangle */
static void test_stpsv(int n) {
  int np = n * (n + 1) / 2;
  float *A = xmalloc((size_t)n * n * sizeof(float));
  float *ap = xmalloc((size_t)np * sizeof(float));
  float *x = xmalloc((size_t)n * sizeof(float));
  float *y = xmalloc((size_t)n * sizeof(float));
  float one = 1.0f;
  float zero = 0.0f;

  make_tri_f32(A, n, n, CblasLower, 0);
  pack_lower_f32(ap, A, n);
  fill_f32(x, n, 61);
  memset(y, 0, (size_t)n * sizeof(float));

  ref_l2_mv(n, A, x, y, &one, &zero, 1, 0, 0);
  cblas_stpsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, ap, y,
              1);

  expect_l2_f32("stpsv", n, y, x, n);

  free(A);
  free(ap);
  free(x);
  free(y);
}

/* sspr: packed A := alpha*x*x^T + A, starting from a zero A */
static void test_sspr(int n) {
  int np = n * (n + 1) / 2;
  float *ap = xmalloc((size_t)np * sizeof(float));
  float *apr = xmalloc((size_t)np * sizeof(float));
  float *x = xmalloc((size_t)n * sizeof(float));
  float alpha = 0.8f;

  fill_f32(x, n, 61);
  memset(ap, 0, (size_t)np * sizeof(float));

  cblas_sspr(CblasColMajor, CblasLower, n, alpha, x, 1, ap);
  ref_spr_f32(apr, x, n, alpha);

  expect_l2_f32("sspr", n, ap, apr, np);

  free(ap);
  free(apr);
  free(x);
}

/* sspr2: packed A := alpha*(x*y^T + y*x^T) + A.
 * A starts from the SPR update so this is checked as an accumulation. */
static void test_sspr2(int n) {
  int np = n * (n + 1) / 2;
  float *ap = xmalloc((size_t)np * sizeof(float));
  float *apr = xmalloc((size_t)np * sizeof(float));
  float *x = xmalloc((size_t)n * sizeof(float));
  float *y = xmalloc((size_t)n * sizeof(float));
  float alpha = 0.8f;

  fill_f32(x, n, 61);
  fill_f32(y, n, 67);
  ref_spr_f32(apr, x, n, alpha);
  memcpy(ap, apr, (size_t)np * sizeof(float));

  cblas_sspr2(CblasColMajor, CblasLower, n, alpha, x, 1, y, 1, ap);
  for (int j = 0; j < n; j++)
    for (int i = j; i < n; i++)
      apr[packed_lower(i, j, n)] += alpha * (x[i] * y[j] + y[i] * x[j]);

  expect_l2_f32("sspr2", n, ap, apr, np);

  free(ap);
  free(apr);
  free(x);
  free(y);
}

/* ------------------------------------------------------------------------- */
/* Real double-precision matrix builders                                     */
/* ------------------------------------------------------------------------- */

/* Mirror the lower triangle of A into the strict upper one. */
static void symmetrize_f64(double *A, int n) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < j; i++)
      A[i + j * n] = A[j + i * n];
}

/* A pseudo-random symmetric matrix. */
static void make_sym_f64(double *A, int n, int seed) {
  fill_f64(A, n * n, seed);
  symmetrize_f64(A, n);
}

/* Drop the strict upper triangle, keeping the lower triangular part. */
static void zero_upper_f64(double *A, int n) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < j; i++)
      A[i + j * n] = 0.0;
}

/* Copy the lower triangle of A into packed storage. */
static void pack_lower_f64(double *ap, const double *A, int n) {
  for (int j = 0; j < n; j++)
    for (int i = j; i < n; i++)
      ap[packed_lower(i, j, n)] = A[i + j * n];
}

/* A general band matrix, both dense (A) and in (kl=ku=k) band storage (ab). */
static void make_gband_f64(double *A, double *ab, int n, int k) {
  int ldab = 2 * k + 1;

  memset(A, 0, (size_t)n * n * sizeof(double));
  memset(ab, 0, (size_t)ldab * n * sizeof(double));
  for (int j = 0; j < n; j++)
    for (int i = j - k; i <= j + k; i++) {
      double v;

      if (i < 0 || i >= n)
        continue;
      v = 0.15 + 0.01 * (double)(i + 2 * j);
      A[i + j * n] = v;
      ab[k + i - j + j * ldab] = v;
    }
}

/* A symmetric band matrix, both dense (A) and in lower band storage (ab). */
static void make_sband_f64(double *A, double *ab, int n, int k) {
  int ldab = k + 1;

  memset(A, 0, (size_t)n * n * sizeof(double));
  memset(ab, 0, (size_t)ldab * n * sizeof(double));
  for (int j = 0; j < n; j++)
    for (int i = j; i <= j + k && i < n; i++) {
      double v = (i == j) ? 2.0 : 0.1 + 0.01 * (double)(i + j);

      A[i + j * n] = v;
      A[j + i * n] = v;
      ab[i - j + j * ldab] = v;
    }
}

/* Packed lower triangle of alpha*x*x^T, the SPR reference. */
static void ref_spr_f64(double *ap, const double *x, int n, double alpha) {
  for (int j = 0; j < n; j++)
    for (int i = j; i < n; i++)
      ap[packed_lower(i, j, n)] = alpha * x[i] * x[j];
}

/* ------------------------------------------------------------------------- */
/* Real double-precision Level-2 (compact grid)                              */
/* ------------------------------------------------------------------------- */

/* dgemv: y := alpha*A*x + beta*y */
static void test_dgemv_compact(int n) {
  double *A = xmalloc((size_t)n * n * sizeof(double));
  double *x = xmalloc((size_t)n * sizeof(double));
  double *y = xmalloc((size_t)n * sizeof(double));
  double *yr = xmalloc((size_t)n * sizeof(double));
  double alpha = 0.8;
  double beta = -0.3;

  make_sym_f64(A, n, 60);
  fill_f64(x, n, 61);
  fill_f64(y, n, 62);
  memcpy(yr, y, (size_t)n * sizeof(double));

  ref_l2_mv(n, A, x, yr, &alpha, &beta, 1, 1, 0);
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, alpha, A, n, x, 1, beta, y, 1);

  expect_l2_f64("dgemv", n, y, yr, n);

  free(A);
  free(x);
  free(y);
  free(yr);
}

/* dsymv: y := alpha*A*x + beta*y, symmetric A read from its lower triangle */
static void test_dsymv(int n) {
  double *A = xmalloc((size_t)n * n * sizeof(double));
  double *x = xmalloc((size_t)n * sizeof(double));
  double *y = xmalloc((size_t)n * sizeof(double));
  double *yr = xmalloc((size_t)n * sizeof(double));
  double alpha = 0.8;
  double beta = -0.3;

  make_sym_f64(A, n, 60);
  fill_f64(x, n, 61);
  fill_f64(y, n, 63);
  memcpy(yr, y, (size_t)n * sizeof(double));

  ref_l2_mv(n, A, x, yr, &alpha, &beta, 1, 1, 0);
  cblas_dsymv(CblasColMajor, CblasLower, n, alpha, A, n, x, 1, beta, y, 1);

  expect_l2_f64("dsymv", n, y, yr, n);

  free(A);
  free(x);
  free(y);
  free(yr);
}

/* dtrmv: x := A*x, A lower triangular with a non-unit diagonal */
static void test_dtrmv(int n) {
  double *A = xmalloc((size_t)n * n * sizeof(double));
  double *x = xmalloc((size_t)n * sizeof(double));
  double *y = xmalloc((size_t)n * sizeof(double));
  double *yr = xmalloc((size_t)n * sizeof(double));
  double one = 1.0;
  double zero = 0.0;

  make_tri_f64(A, n, n, CblasLower, 0);
  fill_f64(x, n, 61);
  memcpy(y, x, (size_t)n * sizeof(double));
  memset(yr, 0, (size_t)n * sizeof(double));

  ref_l2_mv(n, A, x, yr, &one, &zero, 1, 1, 0);
  cblas_dtrmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, A, n, y,
              1);

  expect_l2_f64("dtrmv", n, y, yr, n);

  free(A);
  free(x);
  free(y);
  free(yr);
}

/* dtrsv: solve A*y = A*x on the same triangle, so y must come back as x */
static void test_dtrsv(int n) {
  double *A = xmalloc((size_t)n * n * sizeof(double));
  double *x = xmalloc((size_t)n * sizeof(double));
  double *y = xmalloc((size_t)n * sizeof(double));
  double one = 1.0;
  double zero = 0.0;

  make_tri_f64(A, n, n, CblasLower, 0);
  fill_f64(x, n, 61);
  memset(y, 0, (size_t)n * sizeof(double));

  ref_l2_mv(n, A, x, y, &one, &zero, 1, 1, 0);
  cblas_dtrsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, A, n, y,
              1);

  expect_l2_f64("dtrsv", n, y, x, n);

  free(A);
  free(x);
  free(y);
}

/* dtrmv: x := A^T*x, A upper triangular with an implicit unit diagonal */
static void test_dtrmv_upper_unit(int n) {
  double *A = xmalloc((size_t)n * n * sizeof(double));
  double *x = xmalloc((size_t)n * sizeof(double));
  double *y = xmalloc((size_t)n * sizeof(double));
  double *yr = xmalloc((size_t)n * sizeof(double));
  double one = 1.0;
  double zero = 0.0;

  make_tri_f64(A, n, n, CblasUpper, 1);
  fill_f64(x, n, 61);
  memcpy(y, x, (size_t)n * sizeof(double));
  memset(yr, 0, (size_t)n * sizeof(double));

  ref_l2_mv(n, A, x, yr, &one, &zero, 1, 1, 1);
  cblas_dtrmv(CblasColMajor, CblasUpper, CblasTrans, CblasUnit, n, A, n, y, 1);

  expect_l2_f64("dtrmv-upper-unit", n, y, yr, n);

  free(A);
  free(x);
  free(y);
  free(yr);
}

/* dtrsv: solve A^T*y = A^T*x on the upper unit triangle */
static void test_dtrsv_upper_unit(int n) {
  double *A = xmalloc((size_t)n * n * sizeof(double));
  double *x = xmalloc((size_t)n * sizeof(double));
  double *y = xmalloc((size_t)n * sizeof(double));
  double one = 1.0;
  double zero = 0.0;

  make_tri_f64(A, n, n, CblasUpper, 1);
  fill_f64(x, n, 61);
  memset(y, 0, (size_t)n * sizeof(double));

  ref_l2_mv(n, A, x, y, &one, &zero, 1, 1, 1);
  cblas_dtrsv(CblasColMajor, CblasUpper, CblasTrans, CblasUnit, n, A, n, y, 1);

  expect_l2_f64("dtrsv-upper-unit", n, y, x, n);

  free(A);
  free(x);
  free(y);
}

/* dger: A := alpha*x*y^T + A, starting from a zero A */
static void test_dger(int n) {
  int nn = n * n;
  double *A = xmalloc((size_t)nn * sizeof(double));
  double *Ar = xmalloc((size_t)nn * sizeof(double));
  double *x = xmalloc((size_t)n * sizeof(double));
  double *y = xmalloc((size_t)n * sizeof(double));
  double alpha = 0.8;

  fill_f64(x, n, 61);
  fill_f64(y, n, 68);
  memset(A, 0, (size_t)nn * sizeof(double));
  memset(Ar, 0, (size_t)nn * sizeof(double));

  ref_l2_rank(n, Ar, x, y, &alpha, 1, 1, 0, 0, 0);
  cblas_dger(CblasColMajor, n, n, alpha, x, 1, y, 1, A, n);

  expect_l2_f64("dger", n, A, Ar, nn);

  free(A);
  free(Ar);
  free(x);
  free(y);
}

/* dsyr: A := alpha*x*x^T + A, lower triangle only */
static void test_dsyr(int n) {
  int nn = n * n;
  double *A = xmalloc((size_t)nn * sizeof(double));
  double *Ar = xmalloc((size_t)nn * sizeof(double));
  double *x = xmalloc((size_t)n * sizeof(double));
  double alpha = 0.8;

  fill_f64(x, n, 61);
  memset(A, 0, (size_t)nn * sizeof(double));
  memset(Ar, 0, (size_t)nn * sizeof(double));

  ref_l2_rank(n, Ar, x, x, &alpha, 1, 1, 0, 1, 0);
  cblas_dsyr(CblasColMajor, CblasLower, n, alpha, x, 1, A, n);

  expect_l2_f64("dsyr", n, A, Ar, nn);

  free(A);
  free(Ar);
  free(x);
}

/* dsyr2: A := alpha*(x*y^T + y*x^T) + A, lower triangle only.
 * A starts from a rank-1 update so the update is checked as an accumulation. */
static void test_dsyr2(int n) {
  int nn = n * n;
  double *A = xmalloc((size_t)nn * sizeof(double));
  double *Ar = xmalloc((size_t)nn * sizeof(double));
  double *x = xmalloc((size_t)n * sizeof(double));
  double *y = xmalloc((size_t)n * sizeof(double));
  double alpha = 0.8;

  fill_f64(x, n, 61);
  fill_f64(y, n, 68);
  memset(Ar, 0, (size_t)nn * sizeof(double));
  ref_l2_rank(n, Ar, x, x, &alpha, 1, 1, 0, 1, 0);
  memcpy(A, Ar, (size_t)nn * sizeof(double));

  ref_l2_rank(n, Ar, x, y, &alpha, 1, 1, 0, 1, 0);
  ref_l2_rank(n, Ar, y, x, &alpha, 1, 1, 0, 1, 0);
  cblas_dsyr2(CblasColMajor, CblasLower, n, alpha, x, 1, y, 1, A, n);

  expect_l2_f64("dsyr2", n, A, Ar, nn);

  free(A);
  free(Ar);
  free(x);
  free(y);
}

/* dgbmv: y := alpha*A*x + beta*y with A in general band storage */
static void test_dgbmv(int n) {
  int k = band_width(n);
  int ldab = 2 * k + 1;
  double *A = xmalloc((size_t)n * n * sizeof(double));
  double *ab = xmalloc((size_t)ldab * n * sizeof(double));
  double *x = xmalloc((size_t)n * sizeof(double));
  double *y = xmalloc((size_t)n * sizeof(double));
  double *yr = xmalloc((size_t)n * sizeof(double));
  double alpha = 0.8;
  double beta = -0.3;

  make_gband_f64(A, ab, n, k);
  fill_f64(x, n, 61);
  fill_f64(y, n, 64);
  memcpy(yr, y, (size_t)n * sizeof(double));

  ref_l2_mv(n, A, x, yr, &alpha, &beta, 1, 1, 0);
  cblas_dgbmv(CblasColMajor, CblasNoTrans, n, n, k, k, alpha, ab, ldab, x, 1,
              beta, y, 1);

  expect_l2_f64("dgbmv", n, y, yr, n);

  free(A);
  free(ab);
  free(x);
  free(y);
  free(yr);
}

/* dsbmv: y := alpha*A*x + beta*y with A in symmetric band storage */
static void test_dsbmv(int n) {
  int k = band_width(n);
  int ldab = k + 1;
  double *A = xmalloc((size_t)n * n * sizeof(double));
  double *ab = xmalloc((size_t)ldab * n * sizeof(double));
  double *x = xmalloc((size_t)n * sizeof(double));
  double *y = xmalloc((size_t)n * sizeof(double));
  double *yr = xmalloc((size_t)n * sizeof(double));
  double alpha = 0.8;
  double beta = -0.3;

  make_sband_f64(A, ab, n, k);
  fill_f64(x, n, 61);
  fill_f64(y, n, 65);
  memcpy(yr, y, (size_t)n * sizeof(double));

  ref_l2_mv(n, A, x, yr, &alpha, &beta, 1, 1, 0);
  cblas_dsbmv(CblasColMajor, CblasLower, n, k, alpha, ab, ldab, x, 1, beta, y,
              1);

  expect_l2_f64("dsbmv", n, y, yr, n);

  free(A);
  free(ab);
  free(x);
  free(y);
  free(yr);
}

/* dtbmv: x := A*x with A the lower triangle of a band matrix */
static void test_dtbmv(int n) {
  int k = band_width(n);
  int ldab = k + 1;
  double *A = xmalloc((size_t)n * n * sizeof(double));
  double *ab = xmalloc((size_t)ldab * n * sizeof(double));
  double *x = xmalloc((size_t)n * sizeof(double));
  double *y = xmalloc((size_t)n * sizeof(double));
  double *yr = xmalloc((size_t)n * sizeof(double));
  double one = 1.0;
  double zero = 0.0;

  make_sband_f64(A, ab, n, k);
  zero_upper_f64(A, n);
  fill_f64(x, n, 61);
  memcpy(y, x, (size_t)n * sizeof(double));
  memset(yr, 0, (size_t)n * sizeof(double));

  ref_l2_mv(n, A, x, yr, &one, &zero, 1, 1, 0);
  cblas_dtbmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, k, ab,
              ldab, y, 1);

  expect_l2_f64("dtbmv", n, y, yr, n);

  free(A);
  free(ab);
  free(x);
  free(y);
  free(yr);
}

/* dtbsv: solve A*y = A*x on the same triangular band */
static void test_dtbsv(int n) {
  int k = band_width(n);
  int ldab = k + 1;
  double *A = xmalloc((size_t)n * n * sizeof(double));
  double *ab = xmalloc((size_t)ldab * n * sizeof(double));
  double *x = xmalloc((size_t)n * sizeof(double));
  double *y = xmalloc((size_t)n * sizeof(double));
  double one = 1.0;
  double zero = 0.0;

  make_sband_f64(A, ab, n, k);
  zero_upper_f64(A, n);
  fill_f64(x, n, 61);
  memset(y, 0, (size_t)n * sizeof(double));

  ref_l2_mv(n, A, x, y, &one, &zero, 1, 1, 0);
  cblas_dtbsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, k, ab,
              ldab, y, 1);

  expect_l2_f64("dtbsv", n, y, x, n);

  free(A);
  free(ab);
  free(x);
  free(y);
}

/* dspmv: y := alpha*A*x + beta*y with symmetric A in packed storage */
static void test_dspmv(int n) {
  int np = n * (n + 1) / 2;
  double *A = xmalloc((size_t)n * n * sizeof(double));
  double *ap = xmalloc((size_t)np * sizeof(double));
  double *x = xmalloc((size_t)n * sizeof(double));
  double *y = xmalloc((size_t)n * sizeof(double));
  double *yr = xmalloc((size_t)n * sizeof(double));
  double alpha = 0.8;
  double beta = -0.3;

  make_tri_f64(A, n, n, CblasLower, 0);
  pack_lower_f64(ap, A, n);
  symmetrize_f64(A, n);

  fill_f64(x, n, 61);
  fill_f64(y, n, 66);
  memcpy(yr, y, (size_t)n * sizeof(double));

  ref_l2_mv(n, A, x, yr, &alpha, &beta, 1, 1, 0);
  cblas_dspmv(CblasColMajor, CblasLower, n, alpha, ap, x, 1, beta, y, 1);

  expect_l2_f64("dspmv", n, y, yr, n);

  free(A);
  free(ap);
  free(x);
  free(y);
  free(yr);
}

/* dtpmv: x := A*x with A lower triangular in packed storage */
static void test_dtpmv(int n) {
  int np = n * (n + 1) / 2;
  double *A = xmalloc((size_t)n * n * sizeof(double));
  double *ap = xmalloc((size_t)np * sizeof(double));
  double *x = xmalloc((size_t)n * sizeof(double));
  double *y = xmalloc((size_t)n * sizeof(double));
  double *yr = xmalloc((size_t)n * sizeof(double));
  double one = 1.0;
  double zero = 0.0;

  make_tri_f64(A, n, n, CblasLower, 0);
  pack_lower_f64(ap, A, n);
  fill_f64(x, n, 61);
  memcpy(y, x, (size_t)n * sizeof(double));
  memset(yr, 0, (size_t)n * sizeof(double));

  ref_l2_mv(n, A, x, yr, &one, &zero, 1, 1, 0);
  cblas_dtpmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, ap, y,
              1);

  expect_l2_f64("dtpmv", n, y, yr, n);

  free(A);
  free(ap);
  free(x);
  free(y);
  free(yr);
}

/* dtpsv: solve A*y = A*x on the same packed triangle */
static void test_dtpsv(int n) {
  int np = n * (n + 1) / 2;
  double *A = xmalloc((size_t)n * n * sizeof(double));
  double *ap = xmalloc((size_t)np * sizeof(double));
  double *x = xmalloc((size_t)n * sizeof(double));
  double *y = xmalloc((size_t)n * sizeof(double));
  double one = 1.0;
  double zero = 0.0;

  make_tri_f64(A, n, n, CblasLower, 0);
  pack_lower_f64(ap, A, n);
  fill_f64(x, n, 61);
  memset(y, 0, (size_t)n * sizeof(double));

  ref_l2_mv(n, A, x, y, &one, &zero, 1, 1, 0);
  cblas_dtpsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, ap, y,
              1);

  expect_l2_f64("dtpsv", n, y, x, n);

  free(A);
  free(ap);
  free(x);
  free(y);
}

/* dspr: packed A := alpha*x*x^T + A, starting from a zero A */
static void test_dspr(int n) {
  int np = n * (n + 1) / 2;
  double *ap = xmalloc((size_t)np * sizeof(double));
  double *apr = xmalloc((size_t)np * sizeof(double));
  double *x = xmalloc((size_t)n * sizeof(double));
  double alpha = 0.8;

  fill_f64(x, n, 61);
  memset(ap, 0, (size_t)np * sizeof(double));

  cblas_dspr(CblasColMajor, CblasLower, n, alpha, x, 1, ap);
  ref_spr_f64(apr, x, n, alpha);

  expect_l2_f64("dspr", n, ap, apr, np);

  free(ap);
  free(apr);
  free(x);
}

/* dspr2: packed A := alpha*(x*y^T + y*x^T) + A.
 * A starts from the SPR update so this is checked as an accumulation. */
static void test_dspr2(int n) {
  int np = n * (n + 1) / 2;
  double *ap = xmalloc((size_t)np * sizeof(double));
  double *apr = xmalloc((size_t)np * sizeof(double));
  double *x = xmalloc((size_t)n * sizeof(double));
  double *y = xmalloc((size_t)n * sizeof(double));
  double alpha = 0.8;

  fill_f64(x, n, 61);
  fill_f64(y, n, 67);
  ref_spr_f64(apr, x, n, alpha);
  memcpy(ap, apr, (size_t)np * sizeof(double));

  cblas_dspr2(CblasColMajor, CblasLower, n, alpha, x, 1, y, 1, ap);
  for (int j = 0; j < n; j++)
    for (int i = j; i < n; i++)
      apr[packed_lower(i, j, n)] += alpha * (x[i] * y[j] + y[i] * x[j]);

  expect_l2_f64("dspr2", n, ap, apr, np);

  free(ap);
  free(apr);
  free(x);
  free(y);
}

/* ------------------------------------------------------------------------- */
/* Complex single-precision matrix builders                                  */
/* ------------------------------------------------------------------------- */

/* A triangular matrix with a real diagonal, well conditioned for TRSV. */
static void make_tri_c32(float *A, int n, enum CBLAS_UPLO uplo, int unit) {
  memset(A, 0, (size_t)2 * n * n * sizeof(float));
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      int in_triangle =
          (uplo == CblasLower && i >= j) || (uplo == CblasUpper && i <= j);

      if (!in_triangle)
        continue;
      A[2 * (i + j * n)] = i == j ? (unit ? 1.0f : 2.0f + 0.03f * (float)i)
                                  : 0.08f * (float)(1 + (i + j) % 4);
      A[2 * (i + j * n) + 1] = i == j ? 0.0f : 0.04f * (float)(i - j);
    }
}

/* Mirror the lower triangle of A into the strict upper one, conjugated. */
static void hermitize_c32(float *A, int n) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < j; i++) {
      A[2 * (i + j * n)] = A[2 * (j + i * n)];
      A[2 * (i + j * n) + 1] = -A[2 * (j + i * n) + 1];
    }
}

/* A pseudo-random Hermitian matrix (real diagonal, conjugate off-diagonal). */
static void make_herm_c32(float *A, int n, int seed) {
  fill_c32(A, n * n, seed);
  for (int j = 0; j < n; j++)
    A[2 * (j + j * n) + 1] = 0.0f;
  hermitize_c32(A, n);
}

/* Drop the strict upper triangle, keeping the lower triangular part. */
static void zero_upper_c32(float *A, int n) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < j; i++) {
      A[2 * (i + j * n)] = 0.0f;
      A[2 * (i + j * n) + 1] = 0.0f;
    }
}

/* Copy the lower triangle of A into packed storage. */
static void pack_lower_c32(float *ap, const float *A, int n) {
  for (int j = 0; j < n; j++)
    for (int i = j; i < n; i++) {
      int q = 2 * packed_lower(i, j, n);

      ap[q] = A[2 * (i + j * n)];
      ap[q + 1] = A[2 * (i + j * n) + 1];
    }
}

/* A general band matrix, both dense (A) and in (kl=ku=k) band storage (ab). */
static void make_gband_c32(float *A, float *ab, int n, int k) {
  int ldab = 2 * k + 1;

  memset(A, 0, (size_t)2 * n * n * sizeof(float));
  memset(ab, 0, (size_t)2 * ldab * n * sizeof(float));
  for (int j = 0; j < n; j++)
    for (int i = j - k; i <= j + k; i++) {
      int ia;
      int ib;

      if (i < 0 || i >= n)
        continue;
      ia = 2 * (i + j * n);
      ib = 2 * (k + i - j + j * ldab);
      A[ia] = 0.2f + 0.01f * (float)(i + j);
      A[ia + 1] = 0.03f * (float)(i - j);
      ab[ib] = A[ia];
      ab[ib + 1] = A[ia + 1];
    }
}

/* A Hermitian band matrix, both dense (A) and in lower band storage (ab).
 * The diagonal is real and equal to 2, so the lower triangle is safe to
 * invert in TBSV. */
static void make_hband_c32(float *A, float *ab, int n, int k) {
  int ldab = k + 1;

  memset(A, 0, (size_t)2 * n * n * sizeof(float));
  memset(ab, 0, (size_t)2 * ldab * n * sizeof(float));
  for (int j = 0; j < n; j++)
    for (int i = j; i <= j + k && i < n; i++) {
      int ia = 2 * (i + j * n);
      int it = 2 * (j + i * n);
      int ib = 2 * (i - j + j * ldab);

      A[ia] = (i == j) ? 2.0f : 0.1f + 0.01f * (float)(i + j);
      A[ia + 1] = (i == j) ? 0.0f : 0.03f * (float)(i - j);
      A[it] = A[ia];
      A[it + 1] = -A[ia + 1];
      ab[ib] = A[ia];
      ab[ib + 1] = A[ia + 1];
    }
}

/* Packed lower triangle of alpha*x*x^H, the HPR reference. The diagonal of a
 * Hermitian rank update is real, so its imaginary part is forced to zero. */
static void ref_hpr_c32(float *ap, const float *x, int n, float alpha) {
  for (int j = 0; j < n; j++)
    for (int i = j; i < n; i++) {
      int q = 2 * packed_lower(i, j, n);
      float xir = x[2 * i];
      float xii = x[2 * i + 1];
      float xjr = x[2 * j];
      float xji = x[2 * j + 1];

      ap[q] = alpha * (xir * xjr + xii * xji);
      ap[q + 1] = (i == j) ? 0.0f : alpha * (xii * xjr - xir * xji);
    }
}

/* ------------------------------------------------------------------------- */
/* Complex single-precision Level-2 (compact grid)                           */
/* ------------------------------------------------------------------------- */

/* cgemv: y := alpha*A*x + beta*y */
static void test_cgemv(int n) {
  int z = 2 * n;
  float *A = xmalloc((size_t)2 * n * n * sizeof(float));
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float *yr = xmalloc((size_t)z * sizeof(float));
  float alpha[2] = {0.8f, -0.2f};
  float beta[2] = {-0.3f, 0.1f};

  make_herm_c32(A, n, 70);
  fill_c32(x, n, 71);
  fill_c32(y, n, 72);
  memcpy(yr, y, (size_t)z * sizeof(float));

  ref_l2_mv(n, A, x, yr, alpha, beta, 2, 0, 0);
  cblas_cgemv(CblasColMajor, CblasNoTrans, n, n, alpha, A, n, x, 1, beta, y, 1);

  expect_l2_f32("cgemv", n, y, yr, z);

  free(A);
  free(x);
  free(y);
  free(yr);
}

/* cgemv: y := alpha*A^H*x + beta*y */
static void test_cgemv_conjtrans(int n) {
  int z = 2 * n;
  float *A = xmalloc((size_t)2 * n * n * sizeof(float));
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float *yr = xmalloc((size_t)z * sizeof(float));
  float alpha[2] = {0.8f, -0.2f};
  float beta[2] = {-0.3f, 0.1f};

  make_herm_c32(A, n, 70);
  fill_c32(x, n, 71);
  fill_c32(y, n, 77);
  memcpy(yr, y, (size_t)z * sizeof(float));

  ref_l2_mv(n, A, x, yr, alpha, beta, 2, 0, 1);
  cblas_cgemv(CblasColMajor, CblasConjTrans, n, n, alpha, A, n, x, 1, beta, y,
              1);

  expect_l2_f32("cgemv-conjtrans", n, y, yr, z);

  free(A);
  free(x);
  free(y);
  free(yr);
}

/* chemv: y := alpha*A*x + beta*y, Hermitian A read from its lower triangle */
static void test_chemv(int n) {
  int z = 2 * n;
  float *A = xmalloc((size_t)2 * n * n * sizeof(float));
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float *yr = xmalloc((size_t)z * sizeof(float));
  float alpha[2] = {0.8f, -0.2f};
  float beta[2] = {-0.3f, 0.1f};

  make_herm_c32(A, n, 70);
  fill_c32(x, n, 71);
  fill_c32(y, n, 73);
  memcpy(yr, y, (size_t)z * sizeof(float));

  ref_l2_mv(n, A, x, yr, alpha, beta, 2, 0, 0);
  cblas_chemv(CblasColMajor, CblasLower, n, alpha, A, n, x, 1, beta, y, 1);

  expect_l2_f32("chemv", n, y, yr, z);

  free(A);
  free(x);
  free(y);
  free(yr);
}

/* ctrmv: x := A^H*x, A lower triangular with a non-unit diagonal */
static void test_ctrmv(int n) {
  int z = 2 * n;
  float *A = xmalloc((size_t)2 * n * n * sizeof(float));
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float *yr = xmalloc((size_t)z * sizeof(float));
  float one[2] = {1.0f, 0.0f};
  float zero[2] = {0.0f, 0.0f};

  make_tri_c32(A, n, CblasLower, 0);
  fill_c32(x, n, 71);
  memcpy(y, x, (size_t)z * sizeof(float));
  memset(yr, 0, (size_t)z * sizeof(float));

  ref_l2_mv(n, A, x, yr, one, zero, 2, 0, 1);
  cblas_ctrmv(CblasColMajor, CblasLower, CblasConjTrans, CblasNonUnit, n, A, n,
              y, 1);

  expect_l2_f32("ctrmv", n, y, yr, z);

  free(A);
  free(x);
  free(y);
  free(yr);
}

/* ctrsv: solve A^H*y = A^H*x on the same triangle */
static void test_ctrsv(int n) {
  int z = 2 * n;
  float *A = xmalloc((size_t)2 * n * n * sizeof(float));
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float one[2] = {1.0f, 0.0f};
  float zero[2] = {0.0f, 0.0f};

  make_tri_c32(A, n, CblasLower, 0);
  fill_c32(x, n, 71);
  memset(y, 0, (size_t)z * sizeof(float));

  ref_l2_mv(n, A, x, y, one, zero, 2, 0, 1);
  cblas_ctrsv(CblasColMajor, CblasLower, CblasConjTrans, CblasNonUnit, n, A, n,
              y, 1);

  expect_l2_f32("ctrsv", n, y, x, z);

  free(A);
  free(x);
  free(y);
}

/* ctrmv: x := A^H*x, A upper triangular with an implicit unit diagonal */
static void test_ctrmv_upper_unit(int n) {
  int z = 2 * n;
  float *A = xmalloc((size_t)2 * n * n * sizeof(float));
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float *yr = xmalloc((size_t)z * sizeof(float));
  float one[2] = {1.0f, 0.0f};
  float zero[2] = {0.0f, 0.0f};

  make_tri_c32(A, n, CblasUpper, 1);
  fill_c32(x, n, 71);
  memcpy(y, x, (size_t)z * sizeof(float));
  memset(yr, 0, (size_t)z * sizeof(float));

  ref_l2_mv(n, A, x, yr, one, zero, 2, 0, 1);
  cblas_ctrmv(CblasColMajor, CblasUpper, CblasConjTrans, CblasUnit, n, A, n, y,
              1);

  expect_l2_f32("ctrmv-upper-unit", n, y, yr, z);

  free(A);
  free(x);
  free(y);
  free(yr);
}

/* ctrsv: solve A^H*y = A^H*x on the upper unit triangle */
static void test_ctrsv_upper_unit(int n) {
  int z = 2 * n;
  float *A = xmalloc((size_t)2 * n * n * sizeof(float));
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float one[2] = {1.0f, 0.0f};
  float zero[2] = {0.0f, 0.0f};

  make_tri_c32(A, n, CblasUpper, 1);
  fill_c32(x, n, 71);
  memset(y, 0, (size_t)z * sizeof(float));

  ref_l2_mv(n, A, x, y, one, zero, 2, 0, 1);
  cblas_ctrsv(CblasColMajor, CblasUpper, CblasConjTrans, CblasUnit, n, A, n, y,
              1);

  expect_l2_f32("ctrsv-upper-unit", n, y, x, z);

  free(A);
  free(x);
  free(y);
}

/* cgeru: A := alpha*x*y^T + A, starting from a zero A */
static void test_cgeru(int n) {
  int nz = 2 * n * n;
  float *A = xmalloc((size_t)nz * sizeof(float));
  float *Ar = xmalloc((size_t)nz * sizeof(float));
  float *x = xmalloc((size_t)2 * n * sizeof(float));
  float *y = xmalloc((size_t)2 * n * sizeof(float));
  float alpha[2] = {0.8f, -0.2f};

  fill_c32(x, n, 71);
  fill_c32(y, n, 78);
  memset(A, 0, (size_t)nz * sizeof(float));
  memset(Ar, 0, (size_t)nz * sizeof(float));

  ref_l2_rank(n, Ar, x, y, alpha, 2, 0, 0, 0, 0);
  cblas_cgeru(CblasColMajor, n, n, alpha, x, 1, y, 1, A, n);

  expect_l2_f32("cgeru", n, A, Ar, nz);

  free(A);
  free(Ar);
  free(x);
  free(y);
}

/* cgerc: A := alpha*x*y^H + A, starting from a zero A */
static void test_cgerc(int n) {
  int nz = 2 * n * n;
  float *A = xmalloc((size_t)nz * sizeof(float));
  float *Ar = xmalloc((size_t)nz * sizeof(float));
  float *x = xmalloc((size_t)2 * n * sizeof(float));
  float *y = xmalloc((size_t)2 * n * sizeof(float));
  float alpha[2] = {0.8f, -0.2f};

  fill_c32(x, n, 71);
  fill_c32(y, n, 78);
  memset(A, 0, (size_t)nz * sizeof(float));
  memset(Ar, 0, (size_t)nz * sizeof(float));

  ref_l2_rank(n, Ar, x, y, alpha, 2, 0, 1, 0, 0);
  cblas_cgerc(CblasColMajor, n, n, alpha, x, 1, y, 1, A, n);

  expect_l2_f32("cgerc", n, A, Ar, nz);

  free(A);
  free(Ar);
  free(x);
  free(y);
}

/* cher: A := alpha*x*x^H + A with a real alpha, lower triangle only */
static void test_cher(int n) {
  int nz = 2 * n * n;
  float *A = xmalloc((size_t)nz * sizeof(float));
  float *Ar = xmalloc((size_t)nz * sizeof(float));
  float *x = xmalloc((size_t)2 * n * sizeof(float));
  float alpha = 0.6f;
  float alpha_c[2] = {0.6f, 0.0f};

  fill_c32(x, n, 71);
  memset(A, 0, (size_t)nz * sizeof(float));
  memset(Ar, 0, (size_t)nz * sizeof(float));

  ref_l2_rank(n, Ar, x, x, alpha_c, 2, 0, 1, 1, 1);
  cblas_cher(CblasColMajor, CblasLower, n, alpha, x, 1, A, n);

  expect_l2_f32("cher", n, A, Ar, nz);

  free(A);
  free(Ar);
  free(x);
}

/* cher2: A := alpha*x*y^H + conj(alpha)*y*x^H + A, lower triangle only.
 * A starts from a Hermitian rank-1 update so this is checked as an
 * accumulation. Note the second term carries the conjugate of alpha. */
static void test_cher2(int n) {
  int nz = 2 * n * n;
  float *A = xmalloc((size_t)nz * sizeof(float));
  float *Ar = xmalloc((size_t)nz * sizeof(float));
  float *x = xmalloc((size_t)2 * n * sizeof(float));
  float *y = xmalloc((size_t)2 * n * sizeof(float));
  float alpha[2] = {0.8f, -0.2f};
  float alpha_conj[2] = {0.8f, 0.2f};
  float base[2] = {0.6f, 0.0f};

  fill_c32(x, n, 71);
  fill_c32(y, n, 78);
  memset(Ar, 0, (size_t)nz * sizeof(float));
  ref_l2_rank(n, Ar, x, x, base, 2, 0, 1, 1, 1);
  memcpy(A, Ar, (size_t)nz * sizeof(float));

  ref_l2_rank(n, Ar, x, y, alpha, 2, 0, 1, 1, 1);
  ref_l2_rank(n, Ar, y, x, alpha_conj, 2, 0, 1, 1, 1);
  cblas_cher2(CblasColMajor, CblasLower, n, alpha, x, 1, y, 1, A, n);

  expect_l2_f32("cher2", n, A, Ar, nz);

  free(A);
  free(Ar);
  free(x);
  free(y);
}

/* cgbmv: y := alpha*A*x + beta*y with A in general band storage */
static void test_cgbmv(int n) {
  int z = 2 * n;
  int k = band_width(n);
  int ldab = 2 * k + 1;
  float *A = xmalloc((size_t)2 * n * n * sizeof(float));
  float *ab = xmalloc((size_t)2 * ldab * n * sizeof(float));
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float *yr = xmalloc((size_t)z * sizeof(float));
  float alpha[2] = {0.8f, -0.2f};
  float beta[2] = {-0.3f, 0.1f};

  make_gband_c32(A, ab, n, k);
  fill_c32(x, n, 71);
  fill_c32(y, n, 74);
  memcpy(yr, y, (size_t)z * sizeof(float));

  ref_l2_mv(n, A, x, yr, alpha, beta, 2, 0, 0);
  cblas_cgbmv(CblasColMajor, CblasNoTrans, n, n, k, k, alpha, ab, ldab, x, 1,
              beta, y, 1);

  expect_l2_f32("cgbmv", n, y, yr, z);

  free(A);
  free(ab);
  free(x);
  free(y);
  free(yr);
}

/* chbmv: y := alpha*A*x + beta*y with A in Hermitian band storage */
static void test_chbmv(int n) {
  int z = 2 * n;
  int k = band_width(n);
  int ldab = k + 1;
  float *A = xmalloc((size_t)2 * n * n * sizeof(float));
  float *ab = xmalloc((size_t)2 * ldab * n * sizeof(float));
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float *yr = xmalloc((size_t)z * sizeof(float));
  float alpha[2] = {0.8f, -0.2f};
  float beta[2] = {-0.3f, 0.1f};

  make_hband_c32(A, ab, n, k);
  fill_c32(x, n, 71);
  fill_c32(y, n, 75);
  memcpy(yr, y, (size_t)z * sizeof(float));

  ref_l2_mv(n, A, x, yr, alpha, beta, 2, 0, 0);
  cblas_chbmv(CblasColMajor, CblasLower, n, k, alpha, ab, ldab, x, 1, beta, y,
              1);

  expect_l2_f32("chbmv", n, y, yr, z);

  free(A);
  free(ab);
  free(x);
  free(y);
  free(yr);
}

/* ctbmv: x := A*x with A the lower triangle of a band matrix */
static void test_ctbmv(int n) {
  int z = 2 * n;
  int k = band_width(n);
  int ldab = k + 1;
  float *A = xmalloc((size_t)2 * n * n * sizeof(float));
  float *ab = xmalloc((size_t)2 * ldab * n * sizeof(float));
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float *yr = xmalloc((size_t)z * sizeof(float));
  float one[2] = {1.0f, 0.0f};
  float zero[2] = {0.0f, 0.0f};

  make_hband_c32(A, ab, n, k);
  zero_upper_c32(A, n);
  fill_c32(x, n, 71);
  memcpy(y, x, (size_t)z * sizeof(float));
  memset(yr, 0, (size_t)z * sizeof(float));

  ref_l2_mv(n, A, x, yr, one, zero, 2, 0, 0);
  cblas_ctbmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, k, ab,
              ldab, y, 1);

  expect_l2_f32("ctbmv", n, y, yr, z);

  free(A);
  free(ab);
  free(x);
  free(y);
  free(yr);
}

/* ctbsv: solve A*y = A*x on the same triangular band */
static void test_ctbsv(int n) {
  int z = 2 * n;
  int k = band_width(n);
  int ldab = k + 1;
  float *A = xmalloc((size_t)2 * n * n * sizeof(float));
  float *ab = xmalloc((size_t)2 * ldab * n * sizeof(float));
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float one[2] = {1.0f, 0.0f};
  float zero[2] = {0.0f, 0.0f};

  make_hband_c32(A, ab, n, k);
  zero_upper_c32(A, n);
  fill_c32(x, n, 71);
  memset(y, 0, (size_t)z * sizeof(float));

  ref_l2_mv(n, A, x, y, one, zero, 2, 0, 0);
  cblas_ctbsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, k, ab,
              ldab, y, 1);

  expect_l2_f32("ctbsv", n, y, x, z);

  free(A);
  free(ab);
  free(x);
  free(y);
}

/* chpmv: y := alpha*A*x + beta*y with Hermitian A in packed storage */
static void test_chpmv(int n) {
  int z = 2 * n;
  int np = n * (n + 1) / 2;
  float *A = xmalloc((size_t)2 * n * n * sizeof(float));
  float *ap = xmalloc((size_t)2 * np * sizeof(float));
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float *yr = xmalloc((size_t)z * sizeof(float));
  float alpha[2] = {0.8f, -0.2f};
  float beta[2] = {-0.3f, 0.1f};

  /* Pack the lower triangle, then mirror it conjugated so the dense reference
   * sees the same Hermitian matrix. */
  make_tri_c32(A, n, CblasLower, 0);
  pack_lower_c32(ap, A, n);
  hermitize_c32(A, n);

  fill_c32(x, n, 71);
  fill_c32(y, n, 76);
  memcpy(yr, y, (size_t)z * sizeof(float));

  ref_l2_mv(n, A, x, yr, alpha, beta, 2, 0, 0);
  cblas_chpmv(CblasColMajor, CblasLower, n, alpha, ap, x, 1, beta, y, 1);

  expect_l2_f32("chpmv", n, y, yr, z);

  free(A);
  free(ap);
  free(x);
  free(y);
  free(yr);
}

/* ctpmv: x := A*x with A lower triangular in packed storage */
static void test_ctpmv(int n) {
  int z = 2 * n;
  int np = n * (n + 1) / 2;
  float *A = xmalloc((size_t)2 * n * n * sizeof(float));
  float *ap = xmalloc((size_t)2 * np * sizeof(float));
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float *yr = xmalloc((size_t)z * sizeof(float));
  float one[2] = {1.0f, 0.0f};
  float zero[2] = {0.0f, 0.0f};

  make_tri_c32(A, n, CblasLower, 0);
  pack_lower_c32(ap, A, n);
  fill_c32(x, n, 71);
  memcpy(y, x, (size_t)z * sizeof(float));
  memset(yr, 0, (size_t)z * sizeof(float));

  ref_l2_mv(n, A, x, yr, one, zero, 2, 0, 0);
  cblas_ctpmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, ap, y,
              1);

  expect_l2_f32("ctpmv", n, y, yr, z);

  free(A);
  free(ap);
  free(x);
  free(y);
  free(yr);
}

/* ctpsv: solve A*y = A*x on the same packed triangle */
static void test_ctpsv(int n) {
  int z = 2 * n;
  int np = n * (n + 1) / 2;
  float *A = xmalloc((size_t)2 * n * n * sizeof(float));
  float *ap = xmalloc((size_t)2 * np * sizeof(float));
  float *x = xmalloc((size_t)z * sizeof(float));
  float *y = xmalloc((size_t)z * sizeof(float));
  float one[2] = {1.0f, 0.0f};
  float zero[2] = {0.0f, 0.0f};

  make_tri_c32(A, n, CblasLower, 0);
  pack_lower_c32(ap, A, n);
  fill_c32(x, n, 71);
  memset(y, 0, (size_t)z * sizeof(float));

  ref_l2_mv(n, A, x, y, one, zero, 2, 0, 0);
  cblas_ctpsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, ap, y,
              1);

  expect_l2_f32("ctpsv", n, y, x, z);

  free(A);
  free(ap);
  free(x);
  free(y);
}

/* chpr: packed A := alpha*x*x^H + A with a real alpha, from a zero A */
static void test_chpr(int n) {
  int np = n * (n + 1) / 2;
  float *ap = xmalloc((size_t)2 * np * sizeof(float));
  float *apr = xmalloc((size_t)2 * np * sizeof(float));
  float *x = xmalloc((size_t)2 * n * sizeof(float));
  float alpha = 0.6f;

  fill_c32(x, n, 71);
  memset(ap, 0, (size_t)2 * np * sizeof(float));

  cblas_chpr(CblasColMajor, CblasLower, n, alpha, x, 1, ap);
  ref_hpr_c32(apr, x, n, alpha);

  expect_l2_f32("chpr", n, ap, apr, 2 * np);

  free(ap);
  free(apr);
  free(x);
}

/* chpr2: packed A := alpha*x*y^H + conj(alpha)*y*x^H + A.
 * A starts from the HPR update so this is checked as an accumulation. */
static void test_chpr2(int n) {
  int np = n * (n + 1) / 2;
  float *ap = xmalloc((size_t)2 * np * sizeof(float));
  float *apr = xmalloc((size_t)2 * np * sizeof(float));
  float *x = xmalloc((size_t)2 * n * sizeof(float));
  float *y = xmalloc((size_t)2 * n * sizeof(float));
  float alpha[2] = {0.8f, -0.2f};

  fill_c32(x, n, 71);
  fill_c32(y, n, 79);
  ref_hpr_c32(apr, x, n, 0.6f);
  memcpy(ap, apr, (size_t)2 * np * sizeof(float));

  cblas_chpr2(CblasColMajor, CblasLower, n, alpha, x, 1, y, 1, ap);
  for (int j = 0; j < n; j++)
    for (int i = j; i < n; i++) {
      int q = 2 * packed_lower(i, j, n);
      /* p1 = x_i * conj(y_j) and p2 = y_i * conj(x_j) */
      float p1r = x[2 * i] * y[2 * j] + x[2 * i + 1] * y[2 * j + 1];
      float p1i = x[2 * i + 1] * y[2 * j] - x[2 * i] * y[2 * j + 1];
      float p2r = y[2 * i] * x[2 * j] + y[2 * i + 1] * x[2 * j + 1];
      float p2i = y[2 * i + 1] * x[2 * j] - y[2 * i] * x[2 * j + 1];

      /* alpha*p1 + conj(alpha)*p2; the diagonal stays real. */
      apr[q] +=
          alpha[0] * p1r - alpha[1] * p1i + alpha[0] * p2r + alpha[1] * p2i;
      apr[q + 1] = (i == j) ? 0.0f
                            : apr[q + 1] + alpha[0] * p1i + alpha[1] * p1r +
                                  alpha[0] * p2i - alpha[1] * p2r;
    }

  expect_l2_f32("chpr2", n, ap, apr, 2 * np);

  free(ap);
  free(apr);
  free(x);
  free(y);
}

/* ------------------------------------------------------------------------- */
/* Complex double-precision matrix builders                                  */
/* ------------------------------------------------------------------------- */

/* A triangular matrix with a real diagonal, well conditioned for TRSV. */
static void make_tri_c64(double *A, int n, enum CBLAS_UPLO uplo, int unit) {
  memset(A, 0, (size_t)2 * n * n * sizeof(double));
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      int in_triangle =
          (uplo == CblasLower && i >= j) || (uplo == CblasUpper && i <= j);

      if (!in_triangle)
        continue;
      A[2 * (i + j * n)] = i == j ? (unit ? 1.0 : 2.0 + 0.03 * (double)i)
                                  : 0.08 * (double)(1 + (i + j) % 4);
      A[2 * (i + j * n) + 1] = i == j ? 0.0 : 0.04 * (double)(i - j);
    }
}

/* Mirror the lower triangle of A into the strict upper one, conjugated. */
static void hermitize_c64(double *A, int n) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < j; i++) {
      A[2 * (i + j * n)] = A[2 * (j + i * n)];
      A[2 * (i + j * n) + 1] = -A[2 * (j + i * n) + 1];
    }
}

/* A pseudo-random Hermitian matrix (real diagonal, conjugate off-diagonal). */
static void make_herm_c64(double *A, int n, int seed) {
  fill_c64(A, n * n, seed);
  for (int j = 0; j < n; j++)
    A[2 * (j + j * n) + 1] = 0.0;
  hermitize_c64(A, n);
}

/* Drop the strict upper triangle, keeping the lower triangular part. */
static void zero_upper_c64(double *A, int n) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < j; i++) {
      A[2 * (i + j * n)] = 0.0;
      A[2 * (i + j * n) + 1] = 0.0;
    }
}

/* Copy the lower triangle of A into packed storage. */
static void pack_lower_c64(double *ap, const double *A, int n) {
  for (int j = 0; j < n; j++)
    for (int i = j; i < n; i++) {
      int q = 2 * packed_lower(i, j, n);

      ap[q] = A[2 * (i + j * n)];
      ap[q + 1] = A[2 * (i + j * n) + 1];
    }
}

/* A general band matrix, both dense (A) and in (kl=ku=k) band storage (ab). */
static void make_gband_c64(double *A, double *ab, int n, int k) {
  int ldab = 2 * k + 1;

  memset(A, 0, (size_t)2 * n * n * sizeof(double));
  memset(ab, 0, (size_t)2 * ldab * n * sizeof(double));
  for (int j = 0; j < n; j++)
    for (int i = j - k; i <= j + k; i++) {
      int ia;
      int ib;

      if (i < 0 || i >= n)
        continue;
      ia = 2 * (i + j * n);
      ib = 2 * (k + i - j + j * ldab);
      A[ia] = 0.2 + 0.01 * (double)(i + j);
      A[ia + 1] = 0.03 * (double)(i - j);
      ab[ib] = A[ia];
      ab[ib + 1] = A[ia + 1];
    }
}

/* A Hermitian band matrix, both dense (A) and in lower band storage (ab). */
static void make_hband_c64(double *A, double *ab, int n, int k) {
  int ldab = k + 1;

  memset(A, 0, (size_t)2 * n * n * sizeof(double));
  memset(ab, 0, (size_t)2 * ldab * n * sizeof(double));
  for (int j = 0; j < n; j++)
    for (int i = j; i <= j + k && i < n; i++) {
      int ia = 2 * (i + j * n);
      int it = 2 * (j + i * n);
      int ib = 2 * (i - j + j * ldab);

      A[ia] = (i == j) ? 2.0 : 0.1 + 0.01 * (double)(i + j);
      A[ia + 1] = (i == j) ? 0.0 : 0.03 * (double)(i - j);
      A[it] = A[ia];
      A[it + 1] = -A[ia + 1];
      ab[ib] = A[ia];
      ab[ib + 1] = A[ia + 1];
    }
}

/* Packed lower triangle of alpha*x*x^H, the HPR reference. */
static void ref_hpr_c64(double *ap, const double *x, int n, double alpha) {
  for (int j = 0; j < n; j++)
    for (int i = j; i < n; i++) {
      int q = 2 * packed_lower(i, j, n);
      double xir = x[2 * i];
      double xii = x[2 * i + 1];
      double xjr = x[2 * j];
      double xji = x[2 * j + 1];

      ap[q] = alpha * (xir * xjr + xii * xji);
      ap[q + 1] = (i == j) ? 0.0 : alpha * (xii * xjr - xir * xji);
    }
}

/* ------------------------------------------------------------------------- */
/* Complex double-precision Level-2 (compact grid)                           */
/* ------------------------------------------------------------------------- */

/* zgemv: y := alpha*A*x + beta*y */
static void test_zgemv(int n) {
  int z = 2 * n;
  double *A = xmalloc((size_t)2 * n * n * sizeof(double));
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double *yr = xmalloc((size_t)z * sizeof(double));
  double alpha[2] = {0.8, -0.2};
  double beta[2] = {-0.3, 0.1};

  make_herm_c64(A, n, 70);
  fill_c64(x, n, 71);
  fill_c64(y, n, 72);
  memcpy(yr, y, (size_t)z * sizeof(double));

  ref_l2_mv(n, A, x, yr, alpha, beta, 2, 1, 0);
  cblas_zgemv(CblasColMajor, CblasNoTrans, n, n, alpha, A, n, x, 1, beta, y, 1);

  expect_l2_f64("zgemv", n, y, yr, z);

  free(A);
  free(x);
  free(y);
  free(yr);
}

/* zgemv: y := alpha*A^H*x + beta*y */
static void test_zgemv_conjtrans(int n) {
  int z = 2 * n;
  double *A = xmalloc((size_t)2 * n * n * sizeof(double));
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double *yr = xmalloc((size_t)z * sizeof(double));
  double alpha[2] = {0.8, -0.2};
  double beta[2] = {-0.3, 0.1};

  make_herm_c64(A, n, 70);
  fill_c64(x, n, 71);
  fill_c64(y, n, 77);
  memcpy(yr, y, (size_t)z * sizeof(double));

  ref_l2_mv(n, A, x, yr, alpha, beta, 2, 1, 1);
  cblas_zgemv(CblasColMajor, CblasConjTrans, n, n, alpha, A, n, x, 1, beta, y,
              1);

  expect_l2_f64("zgemv-conjtrans", n, y, yr, z);

  free(A);
  free(x);
  free(y);
  free(yr);
}

/* zhemv: y := alpha*A*x + beta*y, Hermitian A read from its lower triangle */
static void test_zhemv(int n) {
  int z = 2 * n;
  double *A = xmalloc((size_t)2 * n * n * sizeof(double));
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double *yr = xmalloc((size_t)z * sizeof(double));
  double alpha[2] = {0.8, -0.2};
  double beta[2] = {-0.3, 0.1};

  make_herm_c64(A, n, 70);
  fill_c64(x, n, 71);
  fill_c64(y, n, 73);
  memcpy(yr, y, (size_t)z * sizeof(double));

  ref_l2_mv(n, A, x, yr, alpha, beta, 2, 1, 0);
  cblas_zhemv(CblasColMajor, CblasLower, n, alpha, A, n, x, 1, beta, y, 1);

  expect_l2_f64("zhemv", n, y, yr, z);

  free(A);
  free(x);
  free(y);
  free(yr);
}

/* ztrmv: x := A^H*x, A lower triangular with a non-unit diagonal */
static void test_ztrmv(int n) {
  int z = 2 * n;
  double *A = xmalloc((size_t)2 * n * n * sizeof(double));
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double *yr = xmalloc((size_t)z * sizeof(double));
  double one[2] = {1.0, 0.0};
  double zero[2] = {0.0, 0.0};

  make_tri_c64(A, n, CblasLower, 0);
  fill_c64(x, n, 71);
  memcpy(y, x, (size_t)z * sizeof(double));
  memset(yr, 0, (size_t)z * sizeof(double));

  ref_l2_mv(n, A, x, yr, one, zero, 2, 1, 1);
  cblas_ztrmv(CblasColMajor, CblasLower, CblasConjTrans, CblasNonUnit, n, A, n,
              y, 1);

  expect_l2_f64("ztrmv", n, y, yr, z);

  free(A);
  free(x);
  free(y);
  free(yr);
}

/* ztrsv: solve A^H*y = A^H*x on the same triangle */
static void test_ztrsv(int n) {
  int z = 2 * n;
  double *A = xmalloc((size_t)2 * n * n * sizeof(double));
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double one[2] = {1.0, 0.0};
  double zero[2] = {0.0, 0.0};

  make_tri_c64(A, n, CblasLower, 0);
  fill_c64(x, n, 71);
  memset(y, 0, (size_t)z * sizeof(double));

  ref_l2_mv(n, A, x, y, one, zero, 2, 1, 1);
  cblas_ztrsv(CblasColMajor, CblasLower, CblasConjTrans, CblasNonUnit, n, A, n,
              y, 1);

  expect_l2_f64("ztrsv", n, y, x, z);

  free(A);
  free(x);
  free(y);
}

/* ztrmv: x := A^H*x, A upper triangular with an implicit unit diagonal */
static void test_ztrmv_upper_unit(int n) {
  int z = 2 * n;
  double *A = xmalloc((size_t)2 * n * n * sizeof(double));
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double *yr = xmalloc((size_t)z * sizeof(double));
  double one[2] = {1.0, 0.0};
  double zero[2] = {0.0, 0.0};

  make_tri_c64(A, n, CblasUpper, 1);
  fill_c64(x, n, 71);
  memcpy(y, x, (size_t)z * sizeof(double));
  memset(yr, 0, (size_t)z * sizeof(double));

  ref_l2_mv(n, A, x, yr, one, zero, 2, 1, 1);
  cblas_ztrmv(CblasColMajor, CblasUpper, CblasConjTrans, CblasUnit, n, A, n, y,
              1);

  expect_l2_f64("ztrmv-upper-unit", n, y, yr, z);

  free(A);
  free(x);
  free(y);
  free(yr);
}

/* ztrsv: solve A^H*y = A^H*x on the upper unit triangle */
static void test_ztrsv_upper_unit(int n) {
  int z = 2 * n;
  double *A = xmalloc((size_t)2 * n * n * sizeof(double));
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double one[2] = {1.0, 0.0};
  double zero[2] = {0.0, 0.0};

  make_tri_c64(A, n, CblasUpper, 1);
  fill_c64(x, n, 71);
  memset(y, 0, (size_t)z * sizeof(double));

  ref_l2_mv(n, A, x, y, one, zero, 2, 1, 1);
  cblas_ztrsv(CblasColMajor, CblasUpper, CblasConjTrans, CblasUnit, n, A, n, y,
              1);

  expect_l2_f64("ztrsv-upper-unit", n, y, x, z);

  free(A);
  free(x);
  free(y);
}

/* zgeru: A := alpha*x*y^T + A, starting from a zero A */
static void test_zgeru(int n) {
  int nz = 2 * n * n;
  double *A = xmalloc((size_t)nz * sizeof(double));
  double *Ar = xmalloc((size_t)nz * sizeof(double));
  double *x = xmalloc((size_t)2 * n * sizeof(double));
  double *y = xmalloc((size_t)2 * n * sizeof(double));
  double alpha[2] = {0.8, -0.2};

  fill_c64(x, n, 71);
  fill_c64(y, n, 78);
  memset(A, 0, (size_t)nz * sizeof(double));
  memset(Ar, 0, (size_t)nz * sizeof(double));

  ref_l2_rank(n, Ar, x, y, alpha, 2, 1, 0, 0, 0);
  cblas_zgeru(CblasColMajor, n, n, alpha, x, 1, y, 1, A, n);

  expect_l2_f64("zgeru", n, A, Ar, nz);

  free(A);
  free(Ar);
  free(x);
  free(y);
}

/* zgerc: A := alpha*x*y^H + A, starting from a zero A */
static void test_zgerc(int n) {
  int nz = 2 * n * n;
  double *A = xmalloc((size_t)nz * sizeof(double));
  double *Ar = xmalloc((size_t)nz * sizeof(double));
  double *x = xmalloc((size_t)2 * n * sizeof(double));
  double *y = xmalloc((size_t)2 * n * sizeof(double));
  double alpha[2] = {0.8, -0.2};

  fill_c64(x, n, 71);
  fill_c64(y, n, 78);
  memset(A, 0, (size_t)nz * sizeof(double));
  memset(Ar, 0, (size_t)nz * sizeof(double));

  ref_l2_rank(n, Ar, x, y, alpha, 2, 1, 1, 0, 0);
  cblas_zgerc(CblasColMajor, n, n, alpha, x, 1, y, 1, A, n);

  expect_l2_f64("zgerc", n, A, Ar, nz);

  free(A);
  free(Ar);
  free(x);
  free(y);
}

/* zher: A := alpha*x*x^H + A with a real alpha, lower triangle only */
static void test_zher(int n) {
  int nz = 2 * n * n;
  double *A = xmalloc((size_t)nz * sizeof(double));
  double *Ar = xmalloc((size_t)nz * sizeof(double));
  double *x = xmalloc((size_t)2 * n * sizeof(double));
  double alpha = 0.6;
  double alpha_c[2] = {0.6, 0.0};

  fill_c64(x, n, 71);
  memset(A, 0, (size_t)nz * sizeof(double));
  memset(Ar, 0, (size_t)nz * sizeof(double));

  ref_l2_rank(n, Ar, x, x, alpha_c, 2, 1, 1, 1, 1);
  cblas_zher(CblasColMajor, CblasLower, n, alpha, x, 1, A, n);

  expect_l2_f64("zher", n, A, Ar, nz);

  free(A);
  free(Ar);
  free(x);
}

/* zher2: A := alpha*x*y^H + conj(alpha)*y*x^H + A, lower triangle only.
 * A starts from a Hermitian rank-1 update so this is checked as an
 * accumulation. Note the second term carries the conjugate of alpha. */
static void test_zher2(int n) {
  int nz = 2 * n * n;
  double *A = xmalloc((size_t)nz * sizeof(double));
  double *Ar = xmalloc((size_t)nz * sizeof(double));
  double *x = xmalloc((size_t)2 * n * sizeof(double));
  double *y = xmalloc((size_t)2 * n * sizeof(double));
  double alpha[2] = {0.8, -0.2};
  double alpha_conj[2] = {0.8, 0.2};
  double base[2] = {0.6, 0.0};

  fill_c64(x, n, 71);
  fill_c64(y, n, 78);
  memset(Ar, 0, (size_t)nz * sizeof(double));
  ref_l2_rank(n, Ar, x, x, base, 2, 1, 1, 1, 1);
  memcpy(A, Ar, (size_t)nz * sizeof(double));

  ref_l2_rank(n, Ar, x, y, alpha, 2, 1, 1, 1, 1);
  ref_l2_rank(n, Ar, y, x, alpha_conj, 2, 1, 1, 1, 1);
  cblas_zher2(CblasColMajor, CblasLower, n, alpha, x, 1, y, 1, A, n);

  expect_l2_f64("zher2", n, A, Ar, nz);

  free(A);
  free(Ar);
  free(x);
  free(y);
}

/* zgbmv: y := alpha*A*x + beta*y with A in general band storage */
static void test_zgbmv(int n) {
  int z = 2 * n;
  int k = band_width(n);
  int ldab = 2 * k + 1;
  double *A = xmalloc((size_t)2 * n * n * sizeof(double));
  double *ab = xmalloc((size_t)2 * ldab * n * sizeof(double));
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double *yr = xmalloc((size_t)z * sizeof(double));
  double alpha[2] = {0.8, -0.2};
  double beta[2] = {-0.3, 0.1};

  make_gband_c64(A, ab, n, k);
  fill_c64(x, n, 71);
  fill_c64(y, n, 74);
  memcpy(yr, y, (size_t)z * sizeof(double));

  ref_l2_mv(n, A, x, yr, alpha, beta, 2, 1, 0);
  cblas_zgbmv(CblasColMajor, CblasNoTrans, n, n, k, k, alpha, ab, ldab, x, 1,
              beta, y, 1);

  expect_l2_f64("zgbmv", n, y, yr, z);

  free(A);
  free(ab);
  free(x);
  free(y);
  free(yr);
}

/* zhbmv: y := alpha*A*x + beta*y with A in Hermitian band storage */
static void test_zhbmv(int n) {
  int z = 2 * n;
  int k = band_width(n);
  int ldab = k + 1;
  double *A = xmalloc((size_t)2 * n * n * sizeof(double));
  double *ab = xmalloc((size_t)2 * ldab * n * sizeof(double));
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double *yr = xmalloc((size_t)z * sizeof(double));
  double alpha[2] = {0.8, -0.2};
  double beta[2] = {-0.3, 0.1};

  make_hband_c64(A, ab, n, k);
  fill_c64(x, n, 71);
  fill_c64(y, n, 75);
  memcpy(yr, y, (size_t)z * sizeof(double));

  ref_l2_mv(n, A, x, yr, alpha, beta, 2, 1, 0);
  cblas_zhbmv(CblasColMajor, CblasLower, n, k, alpha, ab, ldab, x, 1, beta, y,
              1);

  expect_l2_f64("zhbmv", n, y, yr, z);

  free(A);
  free(ab);
  free(x);
  free(y);
  free(yr);
}

/* ztbmv: x := A*x with A the lower triangle of a band matrix */
static void test_ztbmv(int n) {
  int z = 2 * n;
  int k = band_width(n);
  int ldab = k + 1;
  double *A = xmalloc((size_t)2 * n * n * sizeof(double));
  double *ab = xmalloc((size_t)2 * ldab * n * sizeof(double));
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double *yr = xmalloc((size_t)z * sizeof(double));
  double one[2] = {1.0, 0.0};
  double zero[2] = {0.0, 0.0};

  make_hband_c64(A, ab, n, k);
  zero_upper_c64(A, n);
  fill_c64(x, n, 71);
  memcpy(y, x, (size_t)z * sizeof(double));
  memset(yr, 0, (size_t)z * sizeof(double));

  ref_l2_mv(n, A, x, yr, one, zero, 2, 1, 0);
  cblas_ztbmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, k, ab,
              ldab, y, 1);

  expect_l2_f64("ztbmv", n, y, yr, z);

  free(A);
  free(ab);
  free(x);
  free(y);
  free(yr);
}

/* ztbsv: solve A*y = A*x on the same triangular band */
static void test_ztbsv(int n) {
  int z = 2 * n;
  int k = band_width(n);
  int ldab = k + 1;
  double *A = xmalloc((size_t)2 * n * n * sizeof(double));
  double *ab = xmalloc((size_t)2 * ldab * n * sizeof(double));
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double one[2] = {1.0, 0.0};
  double zero[2] = {0.0, 0.0};

  make_hband_c64(A, ab, n, k);
  zero_upper_c64(A, n);
  fill_c64(x, n, 71);
  memset(y, 0, (size_t)z * sizeof(double));

  ref_l2_mv(n, A, x, y, one, zero, 2, 1, 0);
  cblas_ztbsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, k, ab,
              ldab, y, 1);

  expect_l2_f64("ztbsv", n, y, x, z);

  free(A);
  free(ab);
  free(x);
  free(y);
}

/* zhpmv: y := alpha*A*x + beta*y with Hermitian A in packed storage */
static void test_zhpmv(int n) {
  int z = 2 * n;
  int np = n * (n + 1) / 2;
  double *A = xmalloc((size_t)2 * n * n * sizeof(double));
  double *ap = xmalloc((size_t)2 * np * sizeof(double));
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double *yr = xmalloc((size_t)z * sizeof(double));
  double alpha[2] = {0.8, -0.2};
  double beta[2] = {-0.3, 0.1};

  make_tri_c64(A, n, CblasLower, 0);
  pack_lower_c64(ap, A, n);
  hermitize_c64(A, n);

  fill_c64(x, n, 71);
  fill_c64(y, n, 76);
  memcpy(yr, y, (size_t)z * sizeof(double));

  ref_l2_mv(n, A, x, yr, alpha, beta, 2, 1, 0);
  cblas_zhpmv(CblasColMajor, CblasLower, n, alpha, ap, x, 1, beta, y, 1);

  expect_l2_f64("zhpmv", n, y, yr, z);

  free(A);
  free(ap);
  free(x);
  free(y);
  free(yr);
}

/* ztpmv: x := A*x with A lower triangular in packed storage */
static void test_ztpmv(int n) {
  int z = 2 * n;
  int np = n * (n + 1) / 2;
  double *A = xmalloc((size_t)2 * n * n * sizeof(double));
  double *ap = xmalloc((size_t)2 * np * sizeof(double));
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double *yr = xmalloc((size_t)z * sizeof(double));
  double one[2] = {1.0, 0.0};
  double zero[2] = {0.0, 0.0};

  make_tri_c64(A, n, CblasLower, 0);
  pack_lower_c64(ap, A, n);
  fill_c64(x, n, 71);
  memcpy(y, x, (size_t)z * sizeof(double));
  memset(yr, 0, (size_t)z * sizeof(double));

  ref_l2_mv(n, A, x, yr, one, zero, 2, 1, 0);
  cblas_ztpmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, ap, y,
              1);

  expect_l2_f64("ztpmv", n, y, yr, z);

  free(A);
  free(ap);
  free(x);
  free(y);
  free(yr);
}

/* ztpsv: solve A*y = A*x on the same packed triangle */
static void test_ztpsv(int n) {
  int z = 2 * n;
  int np = n * (n + 1) / 2;
  double *A = xmalloc((size_t)2 * n * n * sizeof(double));
  double *ap = xmalloc((size_t)2 * np * sizeof(double));
  double *x = xmalloc((size_t)z * sizeof(double));
  double *y = xmalloc((size_t)z * sizeof(double));
  double one[2] = {1.0, 0.0};
  double zero[2] = {0.0, 0.0};

  make_tri_c64(A, n, CblasLower, 0);
  pack_lower_c64(ap, A, n);
  fill_c64(x, n, 71);
  memset(y, 0, (size_t)z * sizeof(double));

  ref_l2_mv(n, A, x, y, one, zero, 2, 1, 0);
  cblas_ztpsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, ap, y,
              1);

  expect_l2_f64("ztpsv", n, y, x, z);

  free(A);
  free(ap);
  free(x);
  free(y);
}

/* zhpr: packed A := alpha*x*x^H + A with a real alpha, from a zero A */
static void test_zhpr(int n) {
  int np = n * (n + 1) / 2;
  double *ap = xmalloc((size_t)2 * np * sizeof(double));
  double *apr = xmalloc((size_t)2 * np * sizeof(double));
  double *x = xmalloc((size_t)2 * n * sizeof(double));
  double alpha = 0.6;

  fill_c64(x, n, 71);
  memset(ap, 0, (size_t)2 * np * sizeof(double));

  cblas_zhpr(CblasColMajor, CblasLower, n, alpha, x, 1, ap);
  ref_hpr_c64(apr, x, n, alpha);

  expect_l2_f64("zhpr", n, ap, apr, 2 * np);

  free(ap);
  free(apr);
  free(x);
}

/* zhpr2: packed A := alpha*x*y^H + conj(alpha)*y*x^H + A.
 * A starts from the HPR update so this is checked as an accumulation. */
static void test_zhpr2(int n) {
  int np = n * (n + 1) / 2;
  double *ap = xmalloc((size_t)2 * np * sizeof(double));
  double *apr = xmalloc((size_t)2 * np * sizeof(double));
  double *x = xmalloc((size_t)2 * n * sizeof(double));
  double *y = xmalloc((size_t)2 * n * sizeof(double));
  double alpha[2] = {0.8, -0.2};

  fill_c64(x, n, 71);
  fill_c64(y, n, 79);
  ref_hpr_c64(apr, x, n, 0.6);
  memcpy(ap, apr, (size_t)2 * np * sizeof(double));

  cblas_zhpr2(CblasColMajor, CblasLower, n, alpha, x, 1, y, 1, ap);
  for (int j = 0; j < n; j++)
    for (int i = j; i < n; i++) {
      int q = 2 * packed_lower(i, j, n);
      /* p1 = x_i * conj(y_j) and p2 = y_i * conj(x_j) */
      double p1r = x[2 * i] * y[2 * j] + x[2 * i + 1] * y[2 * j + 1];
      double p1i = x[2 * i + 1] * y[2 * j] - x[2 * i] * y[2 * j + 1];
      double p2r = y[2 * i] * x[2 * j] + y[2 * i + 1] * x[2 * j + 1];
      double p2i = y[2 * i + 1] * x[2 * j] - y[2 * i] * x[2 * j + 1];

      /* alpha*p1 + conj(alpha)*p2; the diagonal stays real. */
      apr[q] +=
          alpha[0] * p1r - alpha[1] * p1i + alpha[0] * p2r + alpha[1] * p2i;
      apr[q + 1] = (i == j) ? 0.0
                            : apr[q + 1] + alpha[0] * p1i + alpha[1] * p1r +
                                  alpha[0] * p2i - alpha[1] * p2r;
    }

  expect_l2_f64("zhpr2", n, ap, apr, 2 * np);

  free(ap);
  free(apr);
  free(x);
  free(y);
}

/* ------------------------------------------------------------------------- */
/* Drivers                                                                   */
/* ------------------------------------------------------------------------- */

static void run_deep_gemv(void) {
  printf("==> L2 GEMV (deep shape / stride grid)\n");
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

        test_sgemv(tr[t], m, nn, 1, 1);
        test_dgemv(tr[t], m, nn, 1, 1);

        if (n <= 36) {
          test_sgemv(tr[t], m, nn, 2, 3);
          test_dgemv(tr[t], m, nn, 2, 3);
        }
      }
    }
  }
}

static void run_compact_l2(void) {
  printf("==> L2 dense, banded, and packed S/D/C/Z families\n");

  for (int i = 0; i < NS_FULL; i++) {
    int n = SIZES_FULL[i];

    test_sgemv_compact(n);
    test_ssymv(n);
    test_strmv(n);
    test_strsv(n);
    test_strmv_upper_unit(n);
    test_strsv_upper_unit(n);
    test_sger(n);
    test_ssyr(n);
    test_ssyr2(n);
    test_sgbmv(n);
    test_ssbmv(n);
    test_stbmv(n);
    test_stbsv(n);
    test_sspmv(n);
    test_stpmv(n);
    test_stpsv(n);
    test_sspr(n);
    test_sspr2(n);

    test_dgemv_compact(n);
    test_dsymv(n);
    test_dtrmv(n);
    test_dtrsv(n);
    test_dtrmv_upper_unit(n);
    test_dtrsv_upper_unit(n);
    test_dger(n);
    test_dsyr(n);
    test_dsyr2(n);
    test_dgbmv(n);
    test_dsbmv(n);
    test_dtbmv(n);
    test_dtbsv(n);
    test_dspmv(n);
    test_dtpmv(n);
    test_dtpsv(n);
    test_dspr(n);
    test_dspr2(n);

    test_cgemv(n);
    test_cgemv_conjtrans(n);
    test_chemv(n);
    test_ctrmv(n);
    test_ctrsv(n);
    test_ctrmv_upper_unit(n);
    test_ctrsv_upper_unit(n);
    test_cgeru(n);
    test_cgerc(n);
    test_cher(n);
    test_cher2(n);
    test_cgbmv(n);
    test_chbmv(n);
    test_ctbmv(n);
    test_ctbsv(n);
    test_chpmv(n);
    test_ctpmv(n);
    test_ctpsv(n);
    test_chpr(n);
    test_chpr2(n);

    test_zgemv(n);
    test_zgemv_conjtrans(n);
    test_zhemv(n);
    test_ztrmv(n);
    test_ztrsv(n);
    test_ztrmv_upper_unit(n);
    test_ztrsv_upper_unit(n);
    test_zgeru(n);
    test_zgerc(n);
    test_zher(n);
    test_zher2(n);
    test_zgbmv(n);
    test_zhbmv(n);
    test_ztbmv(n);
    test_ztbsv(n);
    test_zhpmv(n);
    test_ztpmv(n);
    test_ztpsv(n);
    test_zhpr(n);
    test_zhpr2(n);
  }
}

void check_l2(void) {
  run_deep_gemv();
  run_compact_l2();
}
