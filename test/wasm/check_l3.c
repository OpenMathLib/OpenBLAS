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
 * Level-3 CBLAS checks vs the scalar oracle in ref_l3.c.
 *
 * Level-3 kernels are the blocked ones, so the interesting failures live at
 * block and micro-tile remainders. Two grids are used:
 *
 *   - Deep grids (SIZES_L3 / SIZES_CZ) for GEMM, SYRK, TRMM and TRSM. These
 *     sweep every transpose / uplo / side combination plus a few rectangular
 *     shapes, and reach sizes large enough to exercise the packing path.
 *   - A compact grid (SIZES_FULL) for the remaining standard Level-3 calls
 *     (SYMM, HEMM, complex SYRK/HERK, SYR2K/HER2K and the unit-diagonal
 *     TRMM/TRSM variants), which only need broad remainder sampling.
 *
 * One static test_* function per public CBLAS interface, or per distinct case
 * where a single interface has several meaningful variants (left/right side,
 * unit/non-unit diagonal). Each one allocates, fills, calls the scalar
 * reference, calls OpenBLAS, then expect_close_* / expect_l3_*.
 *
 * All matrices are column major with the leading dimension equal to the
 * number of rows unless the shape forces otherwise, and complex matrices are
 * stored as interleaved re,im pairs.
 */

#include "cases.h"
#include "common.h"
#include "ref.h"
#include "tol.h"

/* ------------------------------------------------------------------------- */
/* Shared helpers                                                            */
/* ------------------------------------------------------------------------- */

/* Compact-grid comparison: the SIZES_FULL families all accumulate over n
 * terms, and the rank-2k ones sum two such products per element, so they get
 * twice the plain Level-3 budget. */
static void expect_l3_f32(const char *op, int n, const float *got,
                          const float *ref, int z) {
  char msg[96];

  snprintf(msg, sizeof(msg), "%s n=%d", op, n);
  expect_close_f32(msg, got, ref, z, tol_s_l3(n) * 2.0f);
}

static void expect_l3_f64(const char *op, int n, const double *got,
                          const double *ref, int z) {
  char msg[96];

  snprintf(msg, sizeof(msg), "%s n=%d", op, n);
  expect_close_f64(msg, got, ref, z, tol_d_l3(n) * 2.0);
}

/* The lower-triangle forms of SYRK / HERK / SYR2K / HER2K leave the strictly
 * upper triangle of C untouched, so zero it in both operands before
 * comparing. `width` is 1 for real and 2 for interleaved complex matrices. */
static void lower32(float *got, float *ref, int n, int width) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < j; i++)
      for (int q = 0; q < width; q++) {
        got[width * (i + j * n) + q] = 0.0f;
        ref[width * (i + j * n) + q] = 0.0f;
      }
}

static void lower64(double *got, double *ref, int n, int width) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < j; i++)
      for (int q = 0; q < width; q++) {
        got[width * (i + j * n) + q] = 0.0;
        ref[width * (i + j * n) + q] = 0.0;
      }
}

/* Complex counterpart of make_tri_f32: a dense n-by-n triangular matrix with
 * the other triangle explicitly zeroed, so the same buffer can be handed to
 * the reference GEMM. The diagonal is real and well away from zero to keep
 * TRSM well conditioned; `unit` materialises the implicit unit diagonal. */
static void tri_c32(float *A, int n, enum CBLAS_UPLO uplo, int unit) {
  memset(A, 0, (size_t)2 * n * n * sizeof(float));
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      int in_triangle = (uplo == CblasLower) ? (i >= j) : (i <= j);
      if (!in_triangle)
        continue;
      if (i == j) {
        A[2 * (i + j * n)] = unit ? 1.0f : 2.2f + 0.02f * (float)i;
        A[2 * (i + j * n) + 1] = 0.0f;
      } else {
        A[2 * (i + j * n)] = 0.07f * (float)(1 + (i + j) % 5);
        A[2 * (i + j * n) + 1] = 0.03f * (float)(i - j);
      }
    }
}

static void tri_c64(double *A, int n, enum CBLAS_UPLO uplo, int unit) {
  memset(A, 0, (size_t)2 * n * n * sizeof(double));
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      int in_triangle = (uplo == CblasLower) ? (i >= j) : (i <= j);
      if (!in_triangle)
        continue;
      if (i == j) {
        A[2 * (i + j * n)] = unit ? 1.0 : 2.2 + 0.02 * (double)i;
        A[2 * (i + j * n) + 1] = 0.0;
      } else {
        A[2 * (i + j * n)] = 0.07 * (double)(1 + (i + j) % 5);
        A[2 * (i + j * n) + 1] = 0.03 * (double)(i - j);
      }
    }
}

/* a := alpha * a over n interleaved complex entries. */
static void scale_c32(float *a, int n, const float *alpha) {
  for (int i = 0; i < n; i++) {
    float re = a[2 * i];
    float im = a[2 * i + 1];
    a[2 * i] = alpha[0] * re - alpha[1] * im;
    a[2 * i + 1] = alpha[0] * im + alpha[1] * re;
  }
}

static void scale_c64(double *a, int n, const double *alpha) {
  for (int i = 0; i < n; i++) {
    double re = a[2 * i];
    double im = a[2 * i + 1];
    a[2 * i] = alpha[0] * re - alpha[1] * im;
    a[2 * i + 1] = alpha[0] * im + alpha[1] * re;
  }
}

/* ------------------------------------------------------------------------- */
/* Deep GEMM grid (all transpose combinations, square and rectangular)       */
/* ------------------------------------------------------------------------- */

/* sgemm: C := alpha * op(A) * op(B) + beta * C */
static void test_sgemm(enum CBLAS_TRANSPOSE ta, enum CBLAS_TRANSPOSE tb, int m,
                       int n, int k) {
  int lda = (ta == CblasNoTrans) ? m : k;
  int ldb = (tb == CblasNoTrans) ? k : n;
  int ldc = m;
  int as = lda * ((ta == CblasNoTrans) ? k : m);
  int bs = ldb * ((tb == CblasNoTrans) ? n : k);
  float *A = xmalloc((size_t)as * sizeof(float));
  float *B = xmalloc((size_t)bs * sizeof(float));
  float *C = xmalloc((size_t)ldc * (size_t)n * sizeof(float));
  float *Cr = xmalloc((size_t)ldc * (size_t)n * sizeof(float));
  float alpha = 1.1f;
  float beta = 0.7f;
  char msg[160];

  fill_f32(A, as, 1);
  fill_f32(B, bs, 2);
  fill_f32(C, ldc * n, 3);
  memcpy(Cr, C, (size_t)ldc * (size_t)n * sizeof(float));

  ref_sgemm(ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, Cr, ldc);
  cblas_sgemm(CblasColMajor, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C,
              ldc);

  snprintf(msg, sizeof(msg), "sgemm ta=%d tb=%d m=%d n=%d k=%d", (int)ta,
           (int)tb, m, n, k);
  expect_close_f32(msg, C, Cr, ldc * n, tol_s_l3(k));

  free(A);
  free(B);
  free(C);
  free(Cr);
}

/* dgemm: C := alpha * op(A) * op(B) + beta * C */
static void test_dgemm(enum CBLAS_TRANSPOSE ta, enum CBLAS_TRANSPOSE tb, int m,
                       int n, int k) {
  int lda = (ta == CblasNoTrans) ? m : k;
  int ldb = (tb == CblasNoTrans) ? k : n;
  int ldc = m;
  int as = lda * ((ta == CblasNoTrans) ? k : m);
  int bs = ldb * ((tb == CblasNoTrans) ? n : k);
  double *A = xmalloc((size_t)as * sizeof(double));
  double *B = xmalloc((size_t)bs * sizeof(double));
  double *C = xmalloc((size_t)ldc * (size_t)n * sizeof(double));
  double *Cr = xmalloc((size_t)ldc * (size_t)n * sizeof(double));
  double alpha = 1.1;
  double beta = 0.7;
  char msg[160];

  fill_f64(A, as, 1);
  fill_f64(B, bs, 2);
  fill_f64(C, ldc * n, 3);
  memcpy(Cr, C, (size_t)ldc * (size_t)n * sizeof(double));

  ref_dgemm(ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, Cr, ldc);
  cblas_dgemm(CblasColMajor, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C,
              ldc);

  snprintf(msg, sizeof(msg), "dgemm ta=%d tb=%d m=%d n=%d k=%d", (int)ta,
           (int)tb, m, n, k);
  expect_close_f64(msg, C, Cr, ldc * n, tol_d_l3(k));

  free(A);
  free(B);
  free(C);
  free(Cr);
}

/* cgemm: C := alpha * op(A) * op(B) + beta * C, op also covering conj-trans */
static void test_cgemm(enum CBLAS_TRANSPOSE ta, enum CBLAS_TRANSPOSE tb, int m,
                       int n, int k) {
  int lda = (ta == CblasNoTrans) ? m : k;
  int ldb = (tb == CblasNoTrans) ? k : n;
  int ldc = m;
  int as = lda * ((ta == CblasNoTrans) ? k : m);
  int bs = ldb * ((tb == CblasNoTrans) ? n : k);
  float *A = xmalloc((size_t)as * 2 * sizeof(float));
  float *B = xmalloc((size_t)bs * 2 * sizeof(float));
  float *C = xmalloc((size_t)ldc * (size_t)n * 2 * sizeof(float));
  float *Cr = xmalloc((size_t)ldc * (size_t)n * 2 * sizeof(float));
  float alpha[2] = {1.1f, -0.3f};
  float beta[2] = {0.7f, 0.2f};
  char msg[160];

  fill_c32(A, as, 1);
  fill_c32(B, bs, 2);
  fill_c32(C, ldc * n, 3);
  memcpy(Cr, C, (size_t)ldc * (size_t)n * 2 * sizeof(float));

  ref_cgemm(ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, Cr, ldc);
  cblas_cgemm(CblasColMajor, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C,
              ldc);

  snprintf(msg, sizeof(msg), "cgemm ta=%d tb=%d m=%d n=%d k=%d", (int)ta,
           (int)tb, m, n, k);
  /* A complex multiply-add is four real products, hence the extra factor. */
  expect_close_f32(msg, C, Cr, ldc * n * 2, tol_s_l3(k) * 2.0f);

  free(A);
  free(B);
  free(C);
  free(Cr);
}

/* zgemm: C := alpha * op(A) * op(B) + beta * C, op also covering conj-trans */
static void test_zgemm(enum CBLAS_TRANSPOSE ta, enum CBLAS_TRANSPOSE tb, int m,
                       int n, int k) {
  int lda = (ta == CblasNoTrans) ? m : k;
  int ldb = (tb == CblasNoTrans) ? k : n;
  int ldc = m;
  int as = lda * ((ta == CblasNoTrans) ? k : m);
  int bs = ldb * ((tb == CblasNoTrans) ? n : k);
  double *A = xmalloc((size_t)as * 2 * sizeof(double));
  double *B = xmalloc((size_t)bs * 2 * sizeof(double));
  double *C = xmalloc((size_t)ldc * (size_t)n * 2 * sizeof(double));
  double *Cr = xmalloc((size_t)ldc * (size_t)n * 2 * sizeof(double));
  double alpha[2] = {1.1, -0.3};
  double beta[2] = {0.7, 0.2};
  char msg[160];

  fill_c64(A, as, 1);
  fill_c64(B, bs, 2);
  fill_c64(C, ldc * n, 3);
  memcpy(Cr, C, (size_t)ldc * (size_t)n * 2 * sizeof(double));

  ref_zgemm(ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, Cr, ldc);
  cblas_zgemm(CblasColMajor, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C,
              ldc);

  snprintf(msg, sizeof(msg), "zgemm ta=%d tb=%d m=%d n=%d k=%d", (int)ta,
           (int)tb, m, n, k);
  expect_close_f64(msg, C, Cr, ldc * n * 2, tol_d_l3(k) * 2.0);

  free(A);
  free(B);
  free(C);
  free(Cr);
}

/* ------------------------------------------------------------------------- */
/* Deep SYRK grid                                                            */
/* ------------------------------------------------------------------------- */

/* ssyrk: C := alpha * A * A^T + beta * C (or A^T * A), one triangle only */
static void test_ssyrk(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n,
                       int k) {
  int lda = (trans == CblasNoTrans) ? n : k;
  int as = lda * ((trans == CblasNoTrans) ? k : n);
  float *A = xmalloc((size_t)as * sizeof(float));
  float *C = xmalloc((size_t)n * (size_t)n * sizeof(float));
  float *Cr = xmalloc((size_t)n * (size_t)n * sizeof(float));
  float alpha = 1.1f;
  float beta = 0.7f;
  char msg[160];

  fill_f32(A, as, 6);
  fill_f32(C, n * n, 7);
  memcpy(Cr, C, (size_t)n * (size_t)n * sizeof(float));

  ref_ssyrk(uplo, trans, n, k, alpha, A, lda, beta, Cr, n);
  cblas_ssyrk(CblasColMajor, uplo, trans, n, k, alpha, A, lda, beta, C, n);

  /* Only compare the triangle that SYRK writes. */
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      int keep = (uplo == CblasUpper) ? (i <= j) : (i >= j);
      if (!keep) {
        C[i + j * n] = 0.0f;
        Cr[i + j * n] = 0.0f;
      }
    }

  snprintf(msg, sizeof(msg), "ssyrk uplo=%d t=%d n=%d k=%d", (int)uplo,
           (int)trans, n, k);
  expect_close_f32(msg, C, Cr, n * n, tol_s_l3(k));

  free(A);
  free(C);
  free(Cr);
}

/* dsyrk: C := alpha * A * A^T + beta * C (or A^T * A), one triangle only */
static void test_dsyrk(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n,
                       int k) {
  int lda = (trans == CblasNoTrans) ? n : k;
  int as = lda * ((trans == CblasNoTrans) ? k : n);
  double *A = xmalloc((size_t)as * sizeof(double));
  double *C = xmalloc((size_t)n * (size_t)n * sizeof(double));
  double *Cr = xmalloc((size_t)n * (size_t)n * sizeof(double));
  double alpha = 1.1;
  double beta = 0.7;
  char msg[160];

  fill_f64(A, as, 6);
  fill_f64(C, n * n, 7);
  memcpy(Cr, C, (size_t)n * (size_t)n * sizeof(double));

  ref_dsyrk(uplo, trans, n, k, alpha, A, lda, beta, Cr, n);
  cblas_dsyrk(CblasColMajor, uplo, trans, n, k, alpha, A, lda, beta, C, n);

  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      int keep = (uplo == CblasUpper) ? (i <= j) : (i >= j);
      if (!keep) {
        C[i + j * n] = 0.0;
        Cr[i + j * n] = 0.0;
      }
    }

  snprintf(msg, sizeof(msg), "dsyrk uplo=%d t=%d n=%d k=%d", (int)uplo,
           (int)trans, n, k);
  expect_close_f64(msg, C, Cr, n * n, tol_d_l3(k));

  free(A);
  free(C);
  free(Cr);
}

/* ------------------------------------------------------------------------- */
/* Deep TRMM / TRSM grids                                                    */
/* ------------------------------------------------------------------------- */

/* strmm: B := alpha * op(A) * B (left) or alpha * B * op(A) (right) */
static void test_strmm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
                       enum CBLAS_TRANSPOSE t, int m, int n) {
  int ka = (side == CblasLeft) ? m : n;
  float *A = xmalloc((size_t)ka * (size_t)ka * sizeof(float));
  float *B = xmalloc((size_t)m * (size_t)n * sizeof(float));
  float *Br = xmalloc((size_t)m * (size_t)n * sizeof(float));
  float alpha = 1.1f;
  char msg[160];

  make_tri_f32(A, ka, ka, uplo, 0);
  fill_f32(B, m * n, 8);
  memcpy(Br, B, (size_t)m * (size_t)n * sizeof(float));

  ref_strmm(side, uplo, t, CblasNonUnit, m, n, alpha, A, ka, Br, m);
  cblas_strmm(CblasColMajor, side, uplo, t, CblasNonUnit, m, n, alpha, A, ka, B,
              m);

  snprintf(msg, sizeof(msg), "strmm side=%d uplo=%d t=%d m=%d n=%d", (int)side,
           (int)uplo, (int)t, m, n);
  expect_close_f32(msg, B, Br, m * n, tol_s_l3(ka));

  free(A);
  free(B);
  free(Br);
}

/* dtrmm: B := alpha * op(A) * B (left) or alpha * B * op(A) (right) */
static void test_dtrmm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
                       enum CBLAS_TRANSPOSE t, int m, int n) {
  int ka = (side == CblasLeft) ? m : n;
  double *A = xmalloc((size_t)ka * (size_t)ka * sizeof(double));
  double *B = xmalloc((size_t)m * (size_t)n * sizeof(double));
  double *Br = xmalloc((size_t)m * (size_t)n * sizeof(double));
  double alpha = 1.1;
  char msg[160];

  make_tri_f64(A, ka, ka, uplo, 0);
  fill_f64(B, m * n, 8);
  memcpy(Br, B, (size_t)m * (size_t)n * sizeof(double));

  ref_dtrmm(side, uplo, t, CblasNonUnit, m, n, alpha, A, ka, Br, m);
  cblas_dtrmm(CblasColMajor, side, uplo, t, CblasNonUnit, m, n, alpha, A, ka, B,
              m);

  snprintf(msg, sizeof(msg), "dtrmm side=%d uplo=%d t=%d m=%d n=%d", (int)side,
           (int)uplo, (int)t, m, n);
  expect_close_f64(msg, B, Br, m * n, tol_d_l3(ka));

  free(A);
  free(B);
  free(Br);
}

/* strsm: solve op(A) X = B (left) or X op(A) = B (right), then check the
 * reconstruction B ~= op(A) X rather than comparing against a reference
 * solve, which would only restate the same rounding. */
static void test_strsm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
                       enum CBLAS_TRANSPOSE t, int m, int n) {
  int ka = (side == CblasLeft) ? m : n;
  float *A = xmalloc((size_t)ka * (size_t)ka * sizeof(float));
  float *Ad = xmalloc((size_t)ka * (size_t)ka * sizeof(float));
  float *B0 = xmalloc((size_t)m * (size_t)n * sizeof(float));
  float *X = xmalloc((size_t)m * (size_t)n * sizeof(float));
  float *Bhat = xmalloc((size_t)m * (size_t)n * sizeof(float));
  float alpha = 1.0f;
  char msg[160];

  make_tri_f32(A, ka, ka, uplo, 0);
  fill_f32(B0, m * n, 9);
  memcpy(X, B0, (size_t)m * (size_t)n * sizeof(float));

  cblas_strsm(CblasColMajor, side, uplo, t, CblasNonUnit, m, n, alpha, A, ka, X,
              m);

  /* Dense copy of the referenced triangle so the reference GEMM can multiply
   * it as an ordinary matrix. */
  memset(Ad, 0, (size_t)ka * (size_t)ka * sizeof(float));
  for (int j = 0; j < ka; j++)
    for (int i = 0; i < ka; i++) {
      int keep = (uplo == CblasLower) ? (i >= j) : (i <= j);
      if (keep)
        Ad[i + j * ka] = A[i + j * ka];
    }

  if (side == CblasLeft)
    ref_sgemm(t, CblasNoTrans, m, n, m, 1.0f, Ad, ka, X, m, 0.0f, Bhat, m);
  else
    ref_sgemm(CblasNoTrans, t, m, n, n, 1.0f, X, m, Ad, ka, 0.0f, Bhat, m);

  snprintf(msg, sizeof(msg), "strsm side=%d uplo=%d t=%d m=%d n=%d", (int)side,
           (int)uplo, (int)t, m, n);
  expect_close_f32(msg, Bhat, B0, m * n, tol_s_l3(ka) * 2.0f);

  free(A);
  free(Ad);
  free(B0);
  free(X);
  free(Bhat);
}

/* dtrsm: solve op(A) X = B (left) or X op(A) = B (right), checked by
 * reconstructing B ~= op(A) X. */
static void test_dtrsm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
                       enum CBLAS_TRANSPOSE t, int m, int n) {
  int ka = (side == CblasLeft) ? m : n;
  double *A = xmalloc((size_t)ka * (size_t)ka * sizeof(double));
  double *Ad = xmalloc((size_t)ka * (size_t)ka * sizeof(double));
  double *B0 = xmalloc((size_t)m * (size_t)n * sizeof(double));
  double *X = xmalloc((size_t)m * (size_t)n * sizeof(double));
  double *Bhat = xmalloc((size_t)m * (size_t)n * sizeof(double));
  double alpha = 1.0;
  char msg[160];

  make_tri_f64(A, ka, ka, uplo, 0);
  fill_f64(B0, m * n, 9);
  memcpy(X, B0, (size_t)m * (size_t)n * sizeof(double));

  cblas_dtrsm(CblasColMajor, side, uplo, t, CblasNonUnit, m, n, alpha, A, ka, X,
              m);

  memset(Ad, 0, (size_t)ka * (size_t)ka * sizeof(double));
  for (int j = 0; j < ka; j++)
    for (int i = 0; i < ka; i++) {
      int keep = (uplo == CblasLower) ? (i >= j) : (i <= j);
      if (keep)
        Ad[i + j * ka] = A[i + j * ka];
    }

  if (side == CblasLeft)
    ref_dgemm(t, CblasNoTrans, m, n, m, 1.0, Ad, ka, X, m, 0.0, Bhat, m);
  else
    ref_dgemm(CblasNoTrans, t, m, n, n, 1.0, X, m, Ad, ka, 0.0, Bhat, m);

  snprintf(msg, sizeof(msg), "dtrsm side=%d uplo=%d t=%d m=%d n=%d", (int)side,
           (int)uplo, (int)t, m, n);
  expect_close_f64(msg, Bhat, B0, m * n, tol_d_l3(ka) * 2.0);

  free(A);
  free(Ad);
  free(B0);
  free(X);
  free(Bhat);
}

/* ------------------------------------------------------------------------- */
/* Compact grid: real single precision                                       */
/* ------------------------------------------------------------------------- */

/* ssymm (left, lower): C := alpha * A * B + beta * C with A symmetric */
static void test_ssymm_left(int n) {
  int nn = n * n;
  float *A = xmalloc((size_t)nn * sizeof(float));
  float *B = xmalloc((size_t)nn * sizeof(float));
  float *C = xmalloc((size_t)nn * sizeof(float));
  float *R = xmalloc((size_t)nn * sizeof(float));
  float alpha = 0.8f;
  float beta = -0.3f;

  fill_f32(A, nn, 80);
  fill_f32(B, nn, 81);
  fill_f32(C, nn, 82);
  memcpy(R, C, (size_t)nn * sizeof(float));

  ref_ssymm(CblasLeft, CblasLower, n, n, alpha, A, n, B, n, beta, R, n);
  cblas_ssymm(CblasColMajor, CblasLeft, CblasLower, n, n, alpha, A, n, B, n,
              beta, C, n);

  expect_l3_f32("ssymm-left", n, C, R, nn);

  free(A);
  free(B);
  free(C);
  free(R);
}

/* ssymm (right, upper): C := alpha * B * A + beta * C with A symmetric */
static void test_ssymm_right(int n) {
  int nn = n * n;
  float *A = xmalloc((size_t)nn * sizeof(float));
  float *B = xmalloc((size_t)nn * sizeof(float));
  float *C = xmalloc((size_t)nn * sizeof(float));
  float *R = xmalloc((size_t)nn * sizeof(float));
  float alpha = 0.8f;
  float beta = -0.3f;

  fill_f32(A, nn, 80);
  fill_f32(B, nn, 81);
  fill_f32(C, nn, 83);
  memcpy(R, C, (size_t)nn * sizeof(float));

  ref_ssymm(CblasRight, CblasUpper, n, n, alpha, A, n, B, n, beta, R, n);
  cblas_ssymm(CblasColMajor, CblasRight, CblasUpper, n, n, alpha, A, n, B, n,
              beta, C, n);

  expect_l3_f32("ssymm-right", n, C, R, nn);

  free(A);
  free(B);
  free(C);
  free(R);
}

/* ssyr2k (lower): C := alpha * A * B^T + alpha * B * A^T + beta * C */
static void test_ssyr2k(int n) {
  int nn = n * n;
  float *A = xmalloc((size_t)nn * sizeof(float));
  float *B = xmalloc((size_t)nn * sizeof(float));
  float *C = xmalloc((size_t)nn * sizeof(float));
  float *R = xmalloc((size_t)nn * sizeof(float));
  float alpha = 0.8f;
  float beta = -0.3f;

  fill_f32(A, nn, 84);
  fill_f32(B, nn, 85);
  fill_f32(C, nn, 86);
  memcpy(R, C, (size_t)nn * sizeof(float));

  ref_ssyr2k(CblasLower, CblasNoTrans, n, n, alpha, A, n, B, n, beta, R, n);
  cblas_ssyr2k(CblasColMajor, CblasLower, CblasNoTrans, n, n, alpha, A, n, B, n,
               beta, C, n);

  lower32(C, R, n, 1);
  expect_l3_f32("ssyr2k", n, C, R, nn);

  free(A);
  free(B);
  free(C);
  free(R);
}

/* strmm (left, lower, non-unit): B := alpha * A * B */
static void test_strmm_full(int n) {
  int nn = n * n;
  float *A = xmalloc((size_t)nn * sizeof(float));
  float *B = xmalloc((size_t)nn * sizeof(float));
  float *R = xmalloc((size_t)nn * sizeof(float));
  float alpha = 0.8f;

  make_tri_f32(A, n, n, CblasLower, 0);
  fill_f32(B, nn, 87);
  memcpy(R, B, (size_t)nn * sizeof(float));

  ref_strmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, n, n, alpha, A,
            n, R, n);
  cblas_strmm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
              n, n, alpha, A, n, B, n);

  expect_l3_f32("strmm-full", n, B, R, nn);

  free(A);
  free(B);
  free(R);
}

/* strsm (left, lower, non-unit): solving A X = alpha * A * B recovers
 * X = alpha * B, so the expected result is known in closed form. */
static void test_strsm_full(int n) {
  int nn = n * n;
  float *A = xmalloc((size_t)nn * sizeof(float));
  float *B = xmalloc((size_t)nn * sizeof(float));
  float *R = xmalloc((size_t)nn * sizeof(float));
  float alpha = 0.8f;

  make_tri_f32(A, n, n, CblasLower, 0);
  fill_f32(B, nn, 87);

  /* Right-hand side alpha * A * B, built with the scalar reference. */
  ref_sgemm(CblasNoTrans, CblasNoTrans, n, n, n, alpha, A, n, B, n, 0.0f, R, n);
  memcpy(B, R, (size_t)nn * sizeof(float));

  cblas_strsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
              n, n, 1.0f, A, n, B, n);

  fill_f32(R, nn, 87);
  for (int i = 0; i < nn; i++)
    R[i] *= alpha;

  expect_l3_f32("strsm-full", n, B, R, nn);

  free(A);
  free(B);
  free(R);
}

/* strmm (right, upper, transposed, unit diagonal): B := alpha * B * A^T */
static void test_strmm_right_unit(int n) {
  int nn = n * n;
  float *A = xmalloc((size_t)nn * sizeof(float));
  float *B = xmalloc((size_t)nn * sizeof(float));
  float *R = xmalloc((size_t)nn * sizeof(float));
  float alpha = 0.8f;

  make_tri_f32(A, n, n, CblasUpper, 1);
  fill_f32(B, nn, 88);
  memcpy(R, B, (size_t)nn * sizeof(float));

  ref_strmm(CblasRight, CblasUpper, CblasTrans, CblasUnit, n, n, alpha, A, n, R,
            n);
  cblas_strmm(CblasColMajor, CblasRight, CblasUpper, CblasTrans, CblasUnit, n,
              n, alpha, A, n, B, n);

  expect_l3_f32("strmm-right-unit", n, B, R, nn);

  free(A);
  free(B);
  free(R);
}

/* strsm (right, upper, transposed, unit diagonal): solving X A^T =
 * alpha * B * A^T recovers X = alpha * B. */
static void test_strsm_right_unit(int n) {
  int nn = n * n;
  float *A = xmalloc((size_t)nn * sizeof(float));
  float *B = xmalloc((size_t)nn * sizeof(float));
  float *R = xmalloc((size_t)nn * sizeof(float));
  float alpha = 0.8f;

  make_tri_f32(A, n, n, CblasUpper, 1);
  fill_f32(B, nn, 88);

  ref_sgemm(CblasNoTrans, CblasTrans, n, n, n, alpha, B, n, A, n, 0.0f, R, n);
  memcpy(B, R, (size_t)nn * sizeof(float));

  cblas_strsm(CblasColMajor, CblasRight, CblasUpper, CblasTrans, CblasUnit, n,
              n, 1.0f, A, n, B, n);

  fill_f32(R, nn, 88);
  for (int i = 0; i < nn; i++)
    R[i] *= alpha;

  expect_l3_f32("strsm-right-unit", n, B, R, nn);

  free(A);
  free(B);
  free(R);
}

/* ------------------------------------------------------------------------- */
/* Compact grid: real double precision                                       */
/* ------------------------------------------------------------------------- */

/* dsymm (left, lower): C := alpha * A * B + beta * C with A symmetric */
static void test_dsymm_left(int n) {
  int nn = n * n;
  double *A = xmalloc((size_t)nn * sizeof(double));
  double *B = xmalloc((size_t)nn * sizeof(double));
  double *C = xmalloc((size_t)nn * sizeof(double));
  double *R = xmalloc((size_t)nn * sizeof(double));
  double alpha = 0.8;
  double beta = -0.3;

  fill_f64(A, nn, 80);
  fill_f64(B, nn, 81);
  fill_f64(C, nn, 82);
  memcpy(R, C, (size_t)nn * sizeof(double));

  ref_dsymm(CblasLeft, CblasLower, n, n, alpha, A, n, B, n, beta, R, n);
  cblas_dsymm(CblasColMajor, CblasLeft, CblasLower, n, n, alpha, A, n, B, n,
              beta, C, n);

  expect_l3_f64("dsymm-left", n, C, R, nn);

  free(A);
  free(B);
  free(C);
  free(R);
}

/* dsymm (right, upper): C := alpha * B * A + beta * C with A symmetric */
static void test_dsymm_right(int n) {
  int nn = n * n;
  double *A = xmalloc((size_t)nn * sizeof(double));
  double *B = xmalloc((size_t)nn * sizeof(double));
  double *C = xmalloc((size_t)nn * sizeof(double));
  double *R = xmalloc((size_t)nn * sizeof(double));
  double alpha = 0.8;
  double beta = -0.3;

  fill_f64(A, nn, 80);
  fill_f64(B, nn, 81);
  fill_f64(C, nn, 83);
  memcpy(R, C, (size_t)nn * sizeof(double));

  ref_dsymm(CblasRight, CblasUpper, n, n, alpha, A, n, B, n, beta, R, n);
  cblas_dsymm(CblasColMajor, CblasRight, CblasUpper, n, n, alpha, A, n, B, n,
              beta, C, n);

  expect_l3_f64("dsymm-right", n, C, R, nn);

  free(A);
  free(B);
  free(C);
  free(R);
}

/* dsyr2k (lower): C := alpha * A * B^T + alpha * B * A^T + beta * C */
static void test_dsyr2k(int n) {
  int nn = n * n;
  double *A = xmalloc((size_t)nn * sizeof(double));
  double *B = xmalloc((size_t)nn * sizeof(double));
  double *C = xmalloc((size_t)nn * sizeof(double));
  double *R = xmalloc((size_t)nn * sizeof(double));
  double alpha = 0.8;
  double beta = -0.3;

  fill_f64(A, nn, 84);
  fill_f64(B, nn, 85);
  fill_f64(C, nn, 86);
  memcpy(R, C, (size_t)nn * sizeof(double));

  ref_dsyr2k(CblasLower, CblasNoTrans, n, n, alpha, A, n, B, n, beta, R, n);
  cblas_dsyr2k(CblasColMajor, CblasLower, CblasNoTrans, n, n, alpha, A, n, B, n,
               beta, C, n);

  lower64(C, R, n, 1);
  expect_l3_f64("dsyr2k", n, C, R, nn);

  free(A);
  free(B);
  free(C);
  free(R);
}

/* dtrmm (left, lower, non-unit): B := alpha * A * B */
static void test_dtrmm_full(int n) {
  int nn = n * n;
  double *A = xmalloc((size_t)nn * sizeof(double));
  double *B = xmalloc((size_t)nn * sizeof(double));
  double *R = xmalloc((size_t)nn * sizeof(double));
  double alpha = 0.8;

  make_tri_f64(A, n, n, CblasLower, 0);
  fill_f64(B, nn, 87);
  memcpy(R, B, (size_t)nn * sizeof(double));

  ref_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, n, n, alpha, A,
            n, R, n);
  cblas_dtrmm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
              n, n, alpha, A, n, B, n);

  expect_l3_f64("dtrmm-full", n, B, R, nn);

  free(A);
  free(B);
  free(R);
}

/* dtrsm (left, lower, non-unit): solving A X = alpha * A * B recovers
 * X = alpha * B. */
static void test_dtrsm_full(int n) {
  int nn = n * n;
  double *A = xmalloc((size_t)nn * sizeof(double));
  double *B = xmalloc((size_t)nn * sizeof(double));
  double *R = xmalloc((size_t)nn * sizeof(double));
  double alpha = 0.8;

  make_tri_f64(A, n, n, CblasLower, 0);
  fill_f64(B, nn, 87);

  ref_dgemm(CblasNoTrans, CblasNoTrans, n, n, n, alpha, A, n, B, n, 0.0, R, n);
  memcpy(B, R, (size_t)nn * sizeof(double));

  cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
              n, n, 1.0, A, n, B, n);

  fill_f64(R, nn, 87);
  for (int i = 0; i < nn; i++)
    R[i] *= alpha;

  expect_l3_f64("dtrsm-full", n, B, R, nn);

  free(A);
  free(B);
  free(R);
}

/* dtrmm (right, upper, transposed, unit diagonal): B := alpha * B * A^T */
static void test_dtrmm_right_unit(int n) {
  int nn = n * n;
  double *A = xmalloc((size_t)nn * sizeof(double));
  double *B = xmalloc((size_t)nn * sizeof(double));
  double *R = xmalloc((size_t)nn * sizeof(double));
  double alpha = 0.8;

  make_tri_f64(A, n, n, CblasUpper, 1);
  fill_f64(B, nn, 88);
  memcpy(R, B, (size_t)nn * sizeof(double));

  ref_dtrmm(CblasRight, CblasUpper, CblasTrans, CblasUnit, n, n, alpha, A, n, R,
            n);
  cblas_dtrmm(CblasColMajor, CblasRight, CblasUpper, CblasTrans, CblasUnit, n,
              n, alpha, A, n, B, n);

  expect_l3_f64("dtrmm-right-unit", n, B, R, nn);

  free(A);
  free(B);
  free(R);
}

/* dtrsm (right, upper, transposed, unit diagonal): solving X A^T =
 * alpha * B * A^T recovers X = alpha * B. */
static void test_dtrsm_right_unit(int n) {
  int nn = n * n;
  double *A = xmalloc((size_t)nn * sizeof(double));
  double *B = xmalloc((size_t)nn * sizeof(double));
  double *R = xmalloc((size_t)nn * sizeof(double));
  double alpha = 0.8;

  make_tri_f64(A, n, n, CblasUpper, 1);
  fill_f64(B, nn, 88);

  ref_dgemm(CblasNoTrans, CblasTrans, n, n, n, alpha, B, n, A, n, 0.0, R, n);
  memcpy(B, R, (size_t)nn * sizeof(double));

  cblas_dtrsm(CblasColMajor, CblasRight, CblasUpper, CblasTrans, CblasUnit, n,
              n, 1.0, A, n, B, n);

  fill_f64(R, nn, 88);
  for (int i = 0; i < nn; i++)
    R[i] *= alpha;

  expect_l3_f64("dtrsm-right-unit", n, B, R, nn);

  free(A);
  free(B);
  free(R);
}

/* ------------------------------------------------------------------------- */
/* Compact grid: complex single precision                                    */
/* ------------------------------------------------------------------------- */

/* csymm (left, lower): C := alpha * A * B + beta * C with A complex
 * symmetric (A = A^T, not Hermitian) */
static void test_csymm(int n) {
  int nn = n * n;
  int z = 2 * nn;
  float *A = xmalloc((size_t)z * sizeof(float));
  float *B = xmalloc((size_t)z * sizeof(float));
  float *C = xmalloc((size_t)z * sizeof(float));
  float *R = xmalloc((size_t)z * sizeof(float));
  float alpha[2] = {0.8f, -0.2f};
  float beta[2] = {-0.3f, 0.1f};

  fill_c32(A, nn, 90);
  fill_c32(B, nn, 91);
  fill_c32(C, nn, 92);
  memcpy(R, C, (size_t)z * sizeof(float));

  ref_csymm(CblasLeft, CblasLower, n, n, alpha, A, n, B, n, beta, R, n);
  cblas_csymm(CblasColMajor, CblasLeft, CblasLower, n, n, alpha, A, n, B, n,
              beta, C, n);

  expect_l3_f32("csymm", n, C, R, z);

  free(A);
  free(B);
  free(C);
  free(R);
}

/* chemm (left, lower): C := alpha * A * B + beta * C with A Hermitian */
static void test_chemm(int n) {
  int nn = n * n;
  int z = 2 * nn;
  float *A = xmalloc((size_t)z * sizeof(float));
  float *B = xmalloc((size_t)z * sizeof(float));
  float *C = xmalloc((size_t)z * sizeof(float));
  float *R = xmalloc((size_t)z * sizeof(float));
  float alpha[2] = {0.8f, -0.2f};
  float beta[2] = {-0.3f, 0.1f};

  fill_c32(A, nn, 90);
  fill_c32(B, nn, 91);
  fill_c32(C, nn, 93);
  memcpy(R, C, (size_t)z * sizeof(float));

  ref_chemm(CblasLeft, CblasLower, n, n, alpha, A, n, B, n, beta, R, n);
  cblas_chemm(CblasColMajor, CblasLeft, CblasLower, n, n, alpha, A, n, B, n,
              beta, C, n);

  expect_l3_f32("chemm", n, C, R, z);

  free(A);
  free(B);
  free(C);
  free(R);
}

/* csyrk (lower): C := alpha * A * A^T + beta * C */
static void test_csyrk(int n) {
  int nn = n * n;
  int z = 2 * nn;
  float *A = xmalloc((size_t)z * sizeof(float));
  float *C = xmalloc((size_t)z * sizeof(float));
  float *R = xmalloc((size_t)z * sizeof(float));
  float alpha[2] = {0.8f, -0.2f};
  float beta[2] = {-0.3f, 0.1f};

  fill_c32(A, nn, 94);
  fill_c32(C, nn, 95);
  memcpy(R, C, (size_t)z * sizeof(float));

  ref_csyrk(CblasLower, CblasNoTrans, n, n, alpha, A, n, beta, R, n);
  cblas_csyrk(CblasColMajor, CblasLower, CblasNoTrans, n, n, alpha, A, n, beta,
              C, n);

  lower32(C, R, n, 2);
  expect_l3_f32("csyrk", n, C, R, z);

  free(A);
  free(C);
  free(R);
}

/* cherk (lower): C := alpha * A * A^H + beta * C with real alpha and beta */
static void test_cherk(int n) {
  int nn = n * n;
  int z = 2 * nn;
  float *A = xmalloc((size_t)z * sizeof(float));
  float *C = xmalloc((size_t)z * sizeof(float));
  float *R = xmalloc((size_t)z * sizeof(float));
  float alpha = 0.7f;
  float beta = -0.2f;

  fill_c32(A, nn, 94);
  fill_c32(C, nn, 96);
  /* A Hermitian C has a real diagonal, and HERK keeps it that way. */
  for (int i = 0; i < n; i++)
    C[2 * (i + i * n) + 1] = 0.0f;
  memcpy(R, C, (size_t)z * sizeof(float));

  ref_cherk(CblasLower, CblasNoTrans, n, n, alpha, A, n, beta, R, n);
  cblas_cherk(CblasColMajor, CblasLower, CblasNoTrans, n, n, alpha, A, n, beta,
              C, n);

  lower32(C, R, n, 2);
  expect_l3_f32("cherk", n, C, R, z);

  free(A);
  free(C);
  free(R);
}

/* csyr2k (lower): C := alpha * A * B^T + alpha * B * A^T + beta * C */
static void test_csyr2k(int n) {
  int nn = n * n;
  int z = 2 * nn;
  float *A = xmalloc((size_t)z * sizeof(float));
  float *B = xmalloc((size_t)z * sizeof(float));
  float *C = xmalloc((size_t)z * sizeof(float));
  float *R = xmalloc((size_t)z * sizeof(float));
  float alpha[2] = {0.8f, -0.2f};
  float beta[2] = {-0.3f, 0.1f};

  fill_c32(A, nn, 94);
  fill_c32(B, nn, 97);
  fill_c32(C, nn, 98);
  memcpy(R, C, (size_t)z * sizeof(float));

  ref_csyr2k(CblasLower, CblasNoTrans, n, n, alpha, A, n, B, n, beta, R, n);
  cblas_csyr2k(CblasColMajor, CblasLower, CblasNoTrans, n, n, alpha, A, n, B, n,
               beta, C, n);

  lower32(C, R, n, 2);
  expect_l3_f32("csyr2k", n, C, R, z);

  free(A);
  free(B);
  free(C);
  free(R);
}

/* cher2k (lower): C := alpha * A * B^H + conj(alpha) * B * A^H + beta * C
 * with real beta */
static void test_cher2k(int n) {
  int nn = n * n;
  int z = 2 * nn;
  float *A = xmalloc((size_t)z * sizeof(float));
  float *B = xmalloc((size_t)z * sizeof(float));
  float *C = xmalloc((size_t)z * sizeof(float));
  float *R = xmalloc((size_t)z * sizeof(float));
  float alpha[2] = {0.8f, -0.2f};
  float beta = -0.2f;

  fill_c32(A, nn, 94);
  fill_c32(B, nn, 97);
  fill_c32(C, nn, 101);
  for (int i = 0; i < n; i++)
    C[2 * (i + i * n) + 1] = 0.0f;
  memcpy(R, C, (size_t)z * sizeof(float));

  ref_cher2k(CblasLower, CblasNoTrans, n, n, alpha, A, n, B, n, beta, R, n);
  cblas_cher2k(CblasColMajor, CblasLower, CblasNoTrans, n, n, alpha, A, n, B, n,
               beta, C, n);

  lower32(C, R, n, 2);
  expect_l3_f32("cher2k", n, C, R, z);

  free(A);
  free(B);
  free(C);
  free(R);
}

/* ctrmm (left, lower, conj-transposed, non-unit): B := alpha * A^H * B */
static void test_ctrmm(int n) {
  int nn = n * n;
  int z = 2 * nn;
  float *A = xmalloc((size_t)z * sizeof(float));
  float *B = xmalloc((size_t)z * sizeof(float));
  float *R = xmalloc((size_t)z * sizeof(float));
  float alpha[2] = {0.8f, -0.2f};
  float zero[2] = {0.0f, 0.0f};

  tri_c32(A, n, CblasLower, 0);
  fill_c32(B, nn, 99);

  ref_cgemm(CblasConjTrans, CblasNoTrans, n, n, n, alpha, A, n, B, n, zero, R,
            n);
  cblas_ctrmm(CblasColMajor, CblasLeft, CblasLower, CblasConjTrans,
              CblasNonUnit, n, n, alpha, A, n, B, n);

  expect_l3_f32("ctrmm", n, B, R, z);

  free(A);
  free(B);
  free(R);
}

/* ctrsm (left, lower, conj-transposed, non-unit): solving A^H X =
 * alpha * A^H * B recovers X = alpha * B. */
static void test_ctrsm(int n) {
  int nn = n * n;
  int z = 2 * nn;
  float *A = xmalloc((size_t)z * sizeof(float));
  float *B = xmalloc((size_t)z * sizeof(float));
  float *R = xmalloc((size_t)z * sizeof(float));
  float alpha[2] = {0.8f, -0.2f};
  float one[2] = {1.0f, 0.0f};
  float zero[2] = {0.0f, 0.0f};

  tri_c32(A, n, CblasLower, 0);
  fill_c32(B, nn, 99);

  ref_cgemm(CblasConjTrans, CblasNoTrans, n, n, n, alpha, A, n, B, n, zero, R,
            n);
  memcpy(B, R, (size_t)z * sizeof(float));

  cblas_ctrsm(CblasColMajor, CblasLeft, CblasLower, CblasConjTrans,
              CblasNonUnit, n, n, one, A, n, B, n);

  fill_c32(R, nn, 99);
  scale_c32(R, nn, alpha);

  expect_l3_f32("ctrsm", n, B, R, z);

  free(A);
  free(B);
  free(R);
}

/* ctrmm (right, upper, transposed, unit diagonal): B := alpha * B * A^T */
static void test_ctrmm_right_unit(int n) {
  int nn = n * n;
  int z = 2 * nn;
  float *A = xmalloc((size_t)z * sizeof(float));
  float *B = xmalloc((size_t)z * sizeof(float));
  float *R = xmalloc((size_t)z * sizeof(float));
  float alpha[2] = {0.8f, -0.2f};
  float zero[2] = {0.0f, 0.0f};

  tri_c32(A, n, CblasUpper, 1);
  fill_c32(B, nn, 100);

  ref_cgemm(CblasNoTrans, CblasTrans, n, n, n, alpha, B, n, A, n, zero, R, n);
  cblas_ctrmm(CblasColMajor, CblasRight, CblasUpper, CblasTrans, CblasUnit, n,
              n, alpha, A, n, B, n);

  expect_l3_f32("ctrmm-right-unit", n, B, R, z);

  free(A);
  free(B);
  free(R);
}

/* ctrsm (right, upper, transposed, unit diagonal): solving X A^T =
 * alpha * B * A^T recovers X = alpha * B. */
static void test_ctrsm_right_unit(int n) {
  int nn = n * n;
  int z = 2 * nn;
  float *A = xmalloc((size_t)z * sizeof(float));
  float *B = xmalloc((size_t)z * sizeof(float));
  float *R = xmalloc((size_t)z * sizeof(float));
  float alpha[2] = {0.8f, -0.2f};
  float one[2] = {1.0f, 0.0f};
  float zero[2] = {0.0f, 0.0f};

  tri_c32(A, n, CblasUpper, 1);
  fill_c32(B, nn, 100);

  ref_cgemm(CblasNoTrans, CblasTrans, n, n, n, alpha, B, n, A, n, zero, R, n);
  memcpy(B, R, (size_t)z * sizeof(float));

  cblas_ctrsm(CblasColMajor, CblasRight, CblasUpper, CblasTrans, CblasUnit, n,
              n, one, A, n, B, n);

  fill_c32(R, nn, 100);
  scale_c32(R, nn, alpha);

  expect_l3_f32("ctrsm-right-unit", n, B, R, z);

  free(A);
  free(B);
  free(R);
}

/* ------------------------------------------------------------------------- */
/* Compact grid: complex double precision                                    */
/* ------------------------------------------------------------------------- */

/* zsymm (left, lower): C := alpha * A * B + beta * C with A complex
 * symmetric (A = A^T, not Hermitian) */
static void test_zsymm(int n) {
  int nn = n * n;
  int z = 2 * nn;
  double *A = xmalloc((size_t)z * sizeof(double));
  double *B = xmalloc((size_t)z * sizeof(double));
  double *C = xmalloc((size_t)z * sizeof(double));
  double *R = xmalloc((size_t)z * sizeof(double));
  double alpha[2] = {0.8, -0.2};
  double beta[2] = {-0.3, 0.1};

  fill_c64(A, nn, 90);
  fill_c64(B, nn, 91);
  fill_c64(C, nn, 92);
  memcpy(R, C, (size_t)z * sizeof(double));

  ref_zsymm(CblasLeft, CblasLower, n, n, alpha, A, n, B, n, beta, R, n);
  cblas_zsymm(CblasColMajor, CblasLeft, CblasLower, n, n, alpha, A, n, B, n,
              beta, C, n);

  expect_l3_f64("zsymm", n, C, R, z);

  free(A);
  free(B);
  free(C);
  free(R);
}

/* zhemm (left, lower): C := alpha * A * B + beta * C with A Hermitian */
static void test_zhemm(int n) {
  int nn = n * n;
  int z = 2 * nn;
  double *A = xmalloc((size_t)z * sizeof(double));
  double *B = xmalloc((size_t)z * sizeof(double));
  double *C = xmalloc((size_t)z * sizeof(double));
  double *R = xmalloc((size_t)z * sizeof(double));
  double alpha[2] = {0.8, -0.2};
  double beta[2] = {-0.3, 0.1};

  fill_c64(A, nn, 90);
  fill_c64(B, nn, 91);
  fill_c64(C, nn, 93);
  memcpy(R, C, (size_t)z * sizeof(double));

  ref_zhemm(CblasLeft, CblasLower, n, n, alpha, A, n, B, n, beta, R, n);
  cblas_zhemm(CblasColMajor, CblasLeft, CblasLower, n, n, alpha, A, n, B, n,
              beta, C, n);

  expect_l3_f64("zhemm", n, C, R, z);

  free(A);
  free(B);
  free(C);
  free(R);
}

/* zsyrk (lower): C := alpha * A * A^T + beta * C */
static void test_zsyrk(int n) {
  int nn = n * n;
  int z = 2 * nn;
  double *A = xmalloc((size_t)z * sizeof(double));
  double *C = xmalloc((size_t)z * sizeof(double));
  double *R = xmalloc((size_t)z * sizeof(double));
  double alpha[2] = {0.8, -0.2};
  double beta[2] = {-0.3, 0.1};

  fill_c64(A, nn, 94);
  fill_c64(C, nn, 95);
  memcpy(R, C, (size_t)z * sizeof(double));

  ref_zsyrk(CblasLower, CblasNoTrans, n, n, alpha, A, n, beta, R, n);
  cblas_zsyrk(CblasColMajor, CblasLower, CblasNoTrans, n, n, alpha, A, n, beta,
              C, n);

  lower64(C, R, n, 2);
  expect_l3_f64("zsyrk", n, C, R, z);

  free(A);
  free(C);
  free(R);
}

/* zherk (lower): C := alpha * A * A^H + beta * C with real alpha and beta */
static void test_zherk(int n) {
  int nn = n * n;
  int z = 2 * nn;
  double *A = xmalloc((size_t)z * sizeof(double));
  double *C = xmalloc((size_t)z * sizeof(double));
  double *R = xmalloc((size_t)z * sizeof(double));
  double alpha = 0.7;
  double beta = -0.2;

  fill_c64(A, nn, 94);
  fill_c64(C, nn, 96);
  for (int i = 0; i < n; i++)
    C[2 * (i + i * n) + 1] = 0.0;
  memcpy(R, C, (size_t)z * sizeof(double));

  ref_zherk(CblasLower, CblasNoTrans, n, n, alpha, A, n, beta, R, n);
  cblas_zherk(CblasColMajor, CblasLower, CblasNoTrans, n, n, alpha, A, n, beta,
              C, n);

  lower64(C, R, n, 2);
  expect_l3_f64("zherk", n, C, R, z);

  free(A);
  free(C);
  free(R);
}

/* zsyr2k (lower): C := alpha * A * B^T + alpha * B * A^T + beta * C */
static void test_zsyr2k(int n) {
  int nn = n * n;
  int z = 2 * nn;
  double *A = xmalloc((size_t)z * sizeof(double));
  double *B = xmalloc((size_t)z * sizeof(double));
  double *C = xmalloc((size_t)z * sizeof(double));
  double *R = xmalloc((size_t)z * sizeof(double));
  double alpha[2] = {0.8, -0.2};
  double beta[2] = {-0.3, 0.1};

  fill_c64(A, nn, 94);
  fill_c64(B, nn, 97);
  fill_c64(C, nn, 98);
  memcpy(R, C, (size_t)z * sizeof(double));

  ref_zsyr2k(CblasLower, CblasNoTrans, n, n, alpha, A, n, B, n, beta, R, n);
  cblas_zsyr2k(CblasColMajor, CblasLower, CblasNoTrans, n, n, alpha, A, n, B, n,
               beta, C, n);

  lower64(C, R, n, 2);
  expect_l3_f64("zsyr2k", n, C, R, z);

  free(A);
  free(B);
  free(C);
  free(R);
}

/* zher2k (lower): C := alpha * A * B^H + conj(alpha) * B * A^H + beta * C
 * with real beta */
static void test_zher2k(int n) {
  int nn = n * n;
  int z = 2 * nn;
  double *A = xmalloc((size_t)z * sizeof(double));
  double *B = xmalloc((size_t)z * sizeof(double));
  double *C = xmalloc((size_t)z * sizeof(double));
  double *R = xmalloc((size_t)z * sizeof(double));
  double alpha[2] = {0.8, -0.2};
  double beta = -0.2;

  fill_c64(A, nn, 94);
  fill_c64(B, nn, 97);
  fill_c64(C, nn, 101);
  for (int i = 0; i < n; i++)
    C[2 * (i + i * n) + 1] = 0.0;
  memcpy(R, C, (size_t)z * sizeof(double));

  ref_zher2k(CblasLower, CblasNoTrans, n, n, alpha, A, n, B, n, beta, R, n);
  cblas_zher2k(CblasColMajor, CblasLower, CblasNoTrans, n, n, alpha, A, n, B, n,
               beta, C, n);

  lower64(C, R, n, 2);
  expect_l3_f64("zher2k", n, C, R, z);

  free(A);
  free(B);
  free(C);
  free(R);
}

/* ztrmm (left, lower, conj-transposed, non-unit): B := alpha * A^H * B */
static void test_ztrmm(int n) {
  int nn = n * n;
  int z = 2 * nn;
  double *A = xmalloc((size_t)z * sizeof(double));
  double *B = xmalloc((size_t)z * sizeof(double));
  double *R = xmalloc((size_t)z * sizeof(double));
  double alpha[2] = {0.8, -0.2};
  double zero[2] = {0.0, 0.0};

  tri_c64(A, n, CblasLower, 0);
  fill_c64(B, nn, 99);

  ref_zgemm(CblasConjTrans, CblasNoTrans, n, n, n, alpha, A, n, B, n, zero, R,
            n);
  cblas_ztrmm(CblasColMajor, CblasLeft, CblasLower, CblasConjTrans,
              CblasNonUnit, n, n, alpha, A, n, B, n);

  expect_l3_f64("ztrmm", n, B, R, z);

  free(A);
  free(B);
  free(R);
}

/* ztrsm (left, lower, conj-transposed, non-unit): solving A^H X =
 * alpha * A^H * B recovers X = alpha * B. */
static void test_ztrsm(int n) {
  int nn = n * n;
  int z = 2 * nn;
  double *A = xmalloc((size_t)z * sizeof(double));
  double *B = xmalloc((size_t)z * sizeof(double));
  double *R = xmalloc((size_t)z * sizeof(double));
  double alpha[2] = {0.8, -0.2};
  double one[2] = {1.0, 0.0};
  double zero[2] = {0.0, 0.0};

  tri_c64(A, n, CblasLower, 0);
  fill_c64(B, nn, 99);

  ref_zgemm(CblasConjTrans, CblasNoTrans, n, n, n, alpha, A, n, B, n, zero, R,
            n);
  memcpy(B, R, (size_t)z * sizeof(double));

  cblas_ztrsm(CblasColMajor, CblasLeft, CblasLower, CblasConjTrans,
              CblasNonUnit, n, n, one, A, n, B, n);

  fill_c64(R, nn, 99);
  scale_c64(R, nn, alpha);

  expect_l3_f64("ztrsm", n, B, R, z);

  free(A);
  free(B);
  free(R);
}

/* ztrmm (right, upper, transposed, unit diagonal): B := alpha * B * A^T */
static void test_ztrmm_right_unit(int n) {
  int nn = n * n;
  int z = 2 * nn;
  double *A = xmalloc((size_t)z * sizeof(double));
  double *B = xmalloc((size_t)z * sizeof(double));
  double *R = xmalloc((size_t)z * sizeof(double));
  double alpha[2] = {0.8, -0.2};
  double zero[2] = {0.0, 0.0};

  tri_c64(A, n, CblasUpper, 1);
  fill_c64(B, nn, 100);

  ref_zgemm(CblasNoTrans, CblasTrans, n, n, n, alpha, B, n, A, n, zero, R, n);
  cblas_ztrmm(CblasColMajor, CblasRight, CblasUpper, CblasTrans, CblasUnit, n,
              n, alpha, A, n, B, n);

  expect_l3_f64("ztrmm-right-unit", n, B, R, z);

  free(A);
  free(B);
  free(R);
}

/* ztrsm (right, upper, transposed, unit diagonal): solving X A^T =
 * alpha * B * A^T recovers X = alpha * B. */
static void test_ztrsm_right_unit(int n) {
  int nn = n * n;
  int z = 2 * nn;
  double *A = xmalloc((size_t)z * sizeof(double));
  double *B = xmalloc((size_t)z * sizeof(double));
  double *R = xmalloc((size_t)z * sizeof(double));
  double alpha[2] = {0.8, -0.2};
  double one[2] = {1.0, 0.0};
  double zero[2] = {0.0, 0.0};

  tri_c64(A, n, CblasUpper, 1);
  fill_c64(B, nn, 100);

  ref_zgemm(CblasNoTrans, CblasTrans, n, n, n, alpha, B, n, A, n, zero, R, n);
  memcpy(B, R, (size_t)z * sizeof(double));

  cblas_ztrsm(CblasColMajor, CblasRight, CblasUpper, CblasTrans, CblasUnit, n,
              n, one, A, n, B, n);

  fill_c64(R, nn, 100);
  scale_c64(R, nn, alpha);

  expect_l3_f64("ztrsm-right-unit", n, B, R, z);

  free(A);
  free(B);
  free(R);
}

/* ------------------------------------------------------------------------- */
/* Drivers                                                                   */
/* ------------------------------------------------------------------------- */

static void run_deep_gemm(void) {
  enum CBLAS_TRANSPOSE tr[2] = {CblasNoTrans, CblasTrans};

  printf("==> L3 GEMM\n");
  for (int s = 0; s < NS_L3; s++) {
    int n = SIZES_L3[s];
    use_fill_case(s);
    for (int ia = 0; ia < 2; ia++)
      for (int ib = 0; ib < 2; ib++) {
        test_sgemm(tr[ia], tr[ib], n, n, n);
        test_dgemm(tr[ia], tr[ib], n, n, n);
        /* Small sizes also get a non-square shape with a distinct k. */
        if (n <= 36) {
          int m2 = n + 1;
          int n2 = n > 1 ? n - 1 : n;
          int k2 = n + 2;
          test_sgemm(tr[ia], tr[ib], m2, n2, k2);
          test_dgemm(tr[ia], tr[ib], m2, n2, k2);
        }
      }
  }
}

static void run_deep_complex_gemm(void) {
  enum CBLAS_TRANSPOSE ctr[3] = {CblasNoTrans, CblasTrans, CblasConjTrans};

  printf("==> L3 CGEMM/ZGEMM\n");
  for (int s = 0; s < NS_CZ; s++) {
    int n = SIZES_CZ[s];
    use_fill_case(s);
    for (int ia = 0; ia < 3; ia++)
      for (int ib = 0; ib < 3; ib++) {
        test_cgemm(ctr[ia], ctr[ib], n, n, n);
        test_zgemm(ctr[ia], ctr[ib], n, n, n);
      }
  }
}

static void run_deep_syrk(void) {
  enum CBLAS_TRANSPOSE tr[2] = {CblasNoTrans, CblasTrans};
  enum CBLAS_UPLO uplos[2] = {CblasLower, CblasUpper};

  printf("==> L3 SYRK\n");
  for (int s = 0; s < NS_L3; s++) {
    int n = SIZES_L3[s];
    use_fill_case(s);
    int k = n <= 36 ? n + 3 : n;
    for (int u = 0; u < 2; u++)
      for (int t = 0; t < 2; t++) {
        test_ssyrk(uplos[u], tr[t], n, k);
        test_dsyrk(uplos[u], tr[t], n, k);
      }
  }
}

static void run_deep_trmm(void) {
  enum CBLAS_TRANSPOSE tr[2] = {CblasNoTrans, CblasTrans};
  enum CBLAS_SIDE sides[2] = {CblasLeft, CblasRight};
  enum CBLAS_UPLO uplos[2] = {CblasLower, CblasUpper};

  printf("==> L3 TRMM\n");
  for (int s = 0; s < NS_L3; s++) {
    int n = SIZES_L3[s];
    use_fill_case(s);
    for (int si = 0; si < 2; si++)
      for (int u = 0; u < 2; u++)
        for (int t = 0; t < 2; t++) {
          test_strmm(sides[si], uplos[u], tr[t], n, n);
          test_dtrmm(sides[si], uplos[u], tr[t], n, n);
          if (n <= 36) {
            int m2 = n + 1;
            int n2 = n > 1 ? n - 1 : 1;
            test_strmm(sides[si], uplos[u], tr[t], m2, n2);
            test_dtrmm(sides[si], uplos[u], tr[t], m2, n2);
          }
        }
  }
}

static void run_deep_trsm(void) {
  enum CBLAS_TRANSPOSE tr[2] = {CblasNoTrans, CblasTrans};
  enum CBLAS_SIDE sides[2] = {CblasLeft, CblasRight};
  enum CBLAS_UPLO uplos[2] = {CblasLower, CblasUpper};

  printf("==> L3 TRSM\n");
  for (int s = 0; s < NS_L3; s++) {
    int n = SIZES_L3[s];
    use_fill_case(s);
    for (int si = 0; si < 2; si++)
      for (int u = 0; u < 2; u++)
        for (int t = 0; t < 2; t++) {
          test_strsm(sides[si], uplos[u], tr[t], n, n);
          test_dtrsm(sides[si], uplos[u], tr[t], n, n);
        }
  }
}

static void run_compact_l3(void) {
  printf("==> L3 SYMM/HEMM/SYRK/HERK/SYR2K/HER2K/TRMM/TRSM (compact grid)\n");
  for (int i = 0; i < NS_FULL; i++) {
    int n = SIZES_FULL[i];
    use_fill_case(i);

    test_ssymm_left(n);
    test_ssymm_right(n);
    test_ssyr2k(n);
    test_strmm_full(n);
    test_strsm_full(n);
    test_strmm_right_unit(n);
    test_strsm_right_unit(n);

    test_dsymm_left(n);
    test_dsymm_right(n);
    test_dsyr2k(n);
    test_dtrmm_full(n);
    test_dtrsm_full(n);
    test_dtrmm_right_unit(n);
    test_dtrsm_right_unit(n);

    test_csymm(n);
    test_chemm(n);
    test_csyrk(n);
    test_cherk(n);
    test_csyr2k(n);
    test_cher2k(n);
    test_ctrmm(n);
    test_ctrsm(n);
    test_ctrmm_right_unit(n);
    test_ctrsm_right_unit(n);

    test_zsymm(n);
    test_zhemm(n);
    test_zsyrk(n);
    test_zherk(n);
    test_zsyr2k(n);
    test_zher2k(n);
    test_ztrmm(n);
    test_ztrsm(n);
    test_ztrmm_right_unit(n);
    test_ztrsm_right_unit(n);
  }
}

void check_l3(void) {
  run_deep_gemm();
  run_deep_complex_gemm();
  run_deep_syrk();
  run_deep_trmm();
  run_deep_trsm();
  run_compact_l3();
}
