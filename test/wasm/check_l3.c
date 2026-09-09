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

static void check_sgemm_one(enum CBLAS_TRANSPOSE ta, enum CBLAS_TRANSPOSE tb,
                            int m, int n, int k) {
  int lda = (ta == CblasNoTrans) ? m : k;
  int ldb = (tb == CblasNoTrans) ? k : n;
  int ldc = m;
  int as = lda * ((ta == CblasNoTrans) ? k : m);
  int bs = ldb * ((tb == CblasNoTrans) ? n : k);
  float *A = xmalloc((size_t)as * sizeof(float));
  float *B = xmalloc((size_t)bs * sizeof(float));
  float *C = xmalloc((size_t)ldc * (size_t)n * sizeof(float));
  float *Cr = xmalloc((size_t)ldc * (size_t)n * sizeof(float));
  fill_f32(A, as, 1);
  fill_f32(B, bs, 2);
  fill_f32(C, ldc * n, 3);
  memcpy(Cr, C, (size_t)ldc * (size_t)n * sizeof(float));
  float alpha = 1.1f, beta = 0.7f;
  ref_sgemm(ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, Cr, ldc);
  cblas_sgemm(CblasColMajor, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C,
              ldc);
  float maxe, maxv, tol = tol_s_l3(k > 0 ? k : 1);
  char msg[160];
  snprintf(msg, sizeof(msg), "sgemm ta=%d tb=%d m=%d n=%d k=%d", (int)ta,
           (int)tb, m, n, k);
  if (!close_f32(C, Cr, ldc * n, tol, &maxe, &maxv))
    fail_f32(msg, maxe, maxv, tol);
  else
    pass_one();
  free(A);
  free(B);
  free(C);
  free(Cr);
}

static void check_dgemm_one(enum CBLAS_TRANSPOSE ta, enum CBLAS_TRANSPOSE tb,
                            int m, int n, int k) {
  int lda = (ta == CblasNoTrans) ? m : k;
  int ldb = (tb == CblasNoTrans) ? k : n;
  int ldc = m;
  int as = lda * ((ta == CblasNoTrans) ? k : m);
  int bs = ldb * ((tb == CblasNoTrans) ? n : k);
  double *A = xmalloc((size_t)as * sizeof(double));
  double *B = xmalloc((size_t)bs * sizeof(double));
  double *C = xmalloc((size_t)ldc * (size_t)n * sizeof(double));
  double *Cr = xmalloc((size_t)ldc * (size_t)n * sizeof(double));
  fill_f64(A, as, 1);
  fill_f64(B, bs, 2);
  fill_f64(C, ldc * n, 3);
  memcpy(Cr, C, (size_t)ldc * (size_t)n * sizeof(double));
  double alpha = 1.1, beta = 0.7;
  ref_dgemm(ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, Cr, ldc);
  cblas_dgemm(CblasColMajor, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C,
              ldc);
  double maxe, maxv, tol = tol_d_l3(k > 0 ? k : 1);
  char msg[160];
  snprintf(msg, sizeof(msg), "dgemm ta=%d tb=%d m=%d n=%d k=%d", (int)ta,
           (int)tb, m, n, k);
  if (!close_f64(C, Cr, ldc * n, tol, &maxe, &maxv))
    fail_f64(msg, maxe, maxv, tol);
  else
    pass_one();
  free(A);
  free(B);
  free(C);
  free(Cr);
}

static void check_cgemm_one(enum CBLAS_TRANSPOSE ta, enum CBLAS_TRANSPOSE tb,
                            int m, int n, int k) {
  int lda = (ta == CblasNoTrans) ? m : k;
  int ldb = (tb == CblasNoTrans) ? k : n;
  int ldc = m;
  int as = lda * ((ta == CblasNoTrans) ? k : m);
  int bs = ldb * ((tb == CblasNoTrans) ? n : k);
  float *A = xmalloc((size_t)as * 2 * sizeof(float));
  float *B = xmalloc((size_t)bs * 2 * sizeof(float));
  float *C = xmalloc((size_t)ldc * (size_t)n * 2 * sizeof(float));
  float *Cr = xmalloc((size_t)ldc * (size_t)n * 2 * sizeof(float));
  fill_c32(A, as, 1);
  fill_c32(B, bs, 2);
  fill_c32(C, ldc * n, 3);
  memcpy(Cr, C, (size_t)ldc * (size_t)n * 2 * sizeof(float));
  float alpha[2] = {1.1f, -0.3f};
  float beta[2] = {0.7f, 0.2f};
  ref_cgemm(ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, Cr, ldc);
  cblas_cgemm(CblasColMajor, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C,
              ldc);
  float maxe, maxv, tol = tol_s_l3(k > 0 ? k : 1) * 2.0f;
  char msg[160];
  snprintf(msg, sizeof(msg), "cgemm ta=%d tb=%d m=%d n=%d k=%d", (int)ta,
           (int)tb, m, n, k);
  if (!close_f32(C, Cr, ldc * n * 2, tol, &maxe, &maxv))
    fail_f32(msg, maxe, maxv, tol);
  else
    pass_one();
  free(A);
  free(B);
  free(C);
  free(Cr);
}

static void check_zgemm_one(enum CBLAS_TRANSPOSE ta, enum CBLAS_TRANSPOSE tb,
                            int m, int n, int k) {
  int lda = (ta == CblasNoTrans) ? m : k;
  int ldb = (tb == CblasNoTrans) ? k : n;
  int ldc = m;
  int as = lda * ((ta == CblasNoTrans) ? k : m);
  int bs = ldb * ((tb == CblasNoTrans) ? n : k);
  double *A = xmalloc((size_t)as * 2 * sizeof(double));
  double *B = xmalloc((size_t)bs * 2 * sizeof(double));
  double *C = xmalloc((size_t)ldc * (size_t)n * 2 * sizeof(double));
  double *Cr = xmalloc((size_t)ldc * (size_t)n * 2 * sizeof(double));
  fill_c64(A, as, 1);
  fill_c64(B, bs, 2);
  fill_c64(C, ldc * n, 3);
  memcpy(Cr, C, (size_t)ldc * (size_t)n * 2 * sizeof(double));
  double alpha[2] = {1.1, -0.3};
  double beta[2] = {0.7, 0.2};
  ref_zgemm(ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, Cr, ldc);
  cblas_zgemm(CblasColMajor, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C,
              ldc);
  double maxe, maxv, tol = tol_d_l3(k > 0 ? k : 1) * 2.0;
  char msg[160];
  snprintf(msg, sizeof(msg), "zgemm ta=%d tb=%d m=%d n=%d k=%d", (int)ta,
           (int)tb, m, n, k);
  if (!close_f64(C, Cr, ldc * n * 2, tol, &maxe, &maxv))
    fail_f64(msg, maxe, maxv, tol);
  else
    pass_one();
  free(A);
  free(B);
  free(C);
  free(Cr);
}

static void check_ssyrk_one(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans,
                            int n, int k) {
  int lda = (trans == CblasNoTrans) ? n : k;
  int as = lda * ((trans == CblasNoTrans) ? k : n);
  float *A = xmalloc((size_t)as * sizeof(float));
  float *C = xmalloc((size_t)n * (size_t)n * sizeof(float));
  float *Cr = xmalloc((size_t)n * (size_t)n * sizeof(float));
  fill_f32(A, as, 6);
  fill_f32(C, n * n, 7);
  memcpy(Cr, C, (size_t)n * (size_t)n * sizeof(float));
  float alpha = 1.1f, beta = 0.7f;
  ref_ssyrk(uplo, trans, n, k, alpha, A, lda, beta, Cr, n);
  cblas_ssyrk(CblasColMajor, uplo, trans, n, k, alpha, A, lda, beta, C, n);
  /* Only compare triangle that SYRK writes. */
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      int keep = (uplo == CblasUpper) ? (i <= j) : (i >= j);
      if (!keep) {
        C[i + j * n] = Cr[i + j * n] = 0.0f;
      }
    }
  float maxe, maxv, tol = tol_s_l3(k > 0 ? k : 1);
  char msg[160];
  snprintf(msg, sizeof(msg), "ssyrk uplo=%d t=%d n=%d k=%d", (int)uplo,
           (int)trans, n, k);
  if (!close_f32(C, Cr, n * n, tol, &maxe, &maxv))
    fail_f32(msg, maxe, maxv, tol);
  else
    pass_one();
  free(A);
  free(C);
  free(Cr);
}

static void check_dsyrk_one(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans,
                            int n, int k) {
  int lda = (trans == CblasNoTrans) ? n : k;
  int as = lda * ((trans == CblasNoTrans) ? k : n);
  double *A = xmalloc((size_t)as * sizeof(double));
  double *C = xmalloc((size_t)n * (size_t)n * sizeof(double));
  double *Cr = xmalloc((size_t)n * (size_t)n * sizeof(double));
  fill_f64(A, as, 6);
  fill_f64(C, n * n, 7);
  memcpy(Cr, C, (size_t)n * (size_t)n * sizeof(double));
  double alpha = 1.1, beta = 0.7;
  ref_dsyrk(uplo, trans, n, k, alpha, A, lda, beta, Cr, n);
  cblas_dsyrk(CblasColMajor, uplo, trans, n, k, alpha, A, lda, beta, C, n);
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      int keep = (uplo == CblasUpper) ? (i <= j) : (i >= j);
      if (!keep) {
        C[i + j * n] = Cr[i + j * n] = 0.0;
      }
    }
  double maxe, maxv, tol = tol_d_l3(k > 0 ? k : 1);
  char msg[160];
  snprintf(msg, sizeof(msg), "dsyrk uplo=%d t=%d n=%d k=%d", (int)uplo,
           (int)trans, n, k);
  if (!close_f64(C, Cr, n * n, tol, &maxe, &maxv))
    fail_f64(msg, maxe, maxv, tol);
  else
    pass_one();
  free(A);
  free(C);
  free(Cr);
}

static void check_strmm_one(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
                            enum CBLAS_TRANSPOSE t, int m, int n) {
  int ka = (side == CblasLeft) ? m : n;
  float *A = xmalloc((size_t)ka * (size_t)ka * sizeof(float));
  float *B = xmalloc((size_t)m * (size_t)n * sizeof(float));
  float *Br = xmalloc((size_t)m * (size_t)n * sizeof(float));
  make_tri_f32(A, ka, ka, uplo, 0);
  fill_f32(B, m * n, 8);
  memcpy(Br, B, (size_t)m * (size_t)n * sizeof(float));
  float alpha = 1.1f;
  ref_strmm(side, uplo, t, CblasNonUnit, m, n, alpha, A, ka, Br, m);
  cblas_strmm(CblasColMajor, side, uplo, t, CblasNonUnit, m, n, alpha, A, ka, B,
              m);
  float maxe, maxv, tol = tol_s_l3(ka);
  char msg[160];
  snprintf(msg, sizeof(msg), "strmm side=%d uplo=%d t=%d m=%d n=%d", (int)side,
           (int)uplo, (int)t, m, n);
  if (!close_f32(B, Br, m * n, tol, &maxe, &maxv))
    fail_f32(msg, maxe, maxv, tol);
  else
    pass_one();
  free(A);
  free(B);
  free(Br);
}

static void check_dtrmm_one(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
                            enum CBLAS_TRANSPOSE t, int m, int n) {
  int ka = (side == CblasLeft) ? m : n;
  double *A = xmalloc((size_t)ka * (size_t)ka * sizeof(double));
  double *B = xmalloc((size_t)m * (size_t)n * sizeof(double));
  double *Br = xmalloc((size_t)m * (size_t)n * sizeof(double));
  make_tri_f64(A, ka, ka, uplo, 0);
  fill_f64(B, m * n, 8);
  memcpy(Br, B, (size_t)m * (size_t)n * sizeof(double));
  double alpha = 1.1;
  ref_dtrmm(side, uplo, t, CblasNonUnit, m, n, alpha, A, ka, Br, m);
  cblas_dtrmm(CblasColMajor, side, uplo, t, CblasNonUnit, m, n, alpha, A, ka, B,
              m);
  double maxe, maxv, tol = tol_d_l3(ka);
  char msg[160];
  snprintf(msg, sizeof(msg), "dtrmm side=%d uplo=%d t=%d m=%d n=%d", (int)side,
           (int)uplo, (int)t, m, n);
  if (!close_f64(B, Br, m * n, tol, &maxe, &maxv))
    fail_f64(msg, maxe, maxv, tol);
  else
    pass_one();
  free(A);
  free(B);
  free(Br);
}

/* Verify TRSM by reconstructing B ≈ op(A) X (or X op(A)). */
static void check_strsm_one(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
                            enum CBLAS_TRANSPOSE t, int m, int n) {
  int ka = (side == CblasLeft) ? m : n;
  float *A = xmalloc((size_t)ka * (size_t)ka * sizeof(float));
  float *B0 = xmalloc((size_t)m * (size_t)n * sizeof(float));
  float *X = xmalloc((size_t)m * (size_t)n * sizeof(float));
  float *Bhat = xmalloc((size_t)m * (size_t)n * sizeof(float));
  make_tri_f32(A, ka, ka, uplo, 0);
  fill_f32(B0, m * n, 9);
  memcpy(X, B0, (size_t)m * (size_t)n * sizeof(float));
  float alpha = 1.0f;
  cblas_strsm(CblasColMajor, side, uplo, t, CblasNonUnit, m, n, alpha, A, ka, X,
              m);
  /* Bhat = op(A) X or X op(A) via dense gemm with triangular zeroed A copied. */
  float *Ad = xmalloc((size_t)ka * (size_t)ka * sizeof(float));
  memset(Ad, 0, (size_t)ka * (size_t)ka * sizeof(float));
  for (int j = 0; j < ka; j++)
    for (int i = 0; i < ka; i++) {
      int keep = (uplo == CblasLower) ? (i >= j) : (i <= j);
      if (keep)
        Ad[i + j * ka] = A[i + j * ka];
    }
  memset(Bhat, 0, (size_t)m * (size_t)n * sizeof(float));
  if (side == CblasLeft)
    ref_sgemm(t, CblasNoTrans, m, n, m, 1.0f, Ad, ka, X, m, 0.0f, Bhat, m);
  else
    ref_sgemm(CblasNoTrans, t, m, n, n, 1.0f, X, m, Ad, ka, 0.0f, Bhat, m);
  float maxe, maxv, tol = tol_s_l3(ka) * 2.0f;
  char msg[160];
  snprintf(msg, sizeof(msg), "strsm side=%d uplo=%d t=%d m=%d n=%d", (int)side,
           (int)uplo, (int)t, m, n);
  if (!close_f32(Bhat, B0, m * n, tol, &maxe, &maxv))
    fail_f32(msg, maxe, maxv, tol);
  else
    pass_one();
  free(A);
  free(Ad);
  free(B0);
  free(X);
  free(Bhat);
}

static void check_dtrsm_one(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
                            enum CBLAS_TRANSPOSE t, int m, int n) {
  int ka = (side == CblasLeft) ? m : n;
  double *A = xmalloc((size_t)ka * (size_t)ka * sizeof(double));
  double *B0 = xmalloc((size_t)m * (size_t)n * sizeof(double));
  double *X = xmalloc((size_t)m * (size_t)n * sizeof(double));
  double *Bhat = xmalloc((size_t)m * (size_t)n * sizeof(double));
  make_tri_f64(A, ka, ka, uplo, 0);
  fill_f64(B0, m * n, 9);
  memcpy(X, B0, (size_t)m * (size_t)n * sizeof(double));
  double alpha = 1.0;
  cblas_dtrsm(CblasColMajor, side, uplo, t, CblasNonUnit, m, n, alpha, A, ka, X,
              m);
  double *Ad = xmalloc((size_t)ka * (size_t)ka * sizeof(double));
  memset(Ad, 0, (size_t)ka * (size_t)ka * sizeof(double));
  for (int j = 0; j < ka; j++)
    for (int i = 0; i < ka; i++) {
      int keep = (uplo == CblasLower) ? (i >= j) : (i <= j);
      if (keep)
        Ad[i + j * ka] = A[i + j * ka];
    }
  memset(Bhat, 0, (size_t)m * (size_t)n * sizeof(double));
  if (side == CblasLeft)
    ref_dgemm(t, CblasNoTrans, m, n, m, 1.0, Ad, ka, X, m, 0.0, Bhat, m);
  else
    ref_dgemm(CblasNoTrans, t, m, n, n, 1.0, X, m, Ad, ka, 0.0, Bhat, m);
  double maxe, maxv, tol = tol_d_l3(ka) * 2.0;
  char msg[160];
  snprintf(msg, sizeof(msg), "dtrsm side=%d uplo=%d t=%d m=%d n=%d", (int)side,
           (int)uplo, (int)t, m, n);
  if (!close_f64(Bhat, B0, m * n, tol, &maxe, &maxv))
    fail_f64(msg, maxe, maxv, tol);
  else
    pass_one();
  free(A);
  free(Ad);
  free(B0);
  free(X);
  free(Bhat);
}

void check_l3(void) {
  enum CBLAS_TRANSPOSE tr[2] = {CblasNoTrans, CblasTrans};
  enum CBLAS_TRANSPOSE ctr[3] = {CblasNoTrans, CblasTrans, CblasConjTrans};
  enum CBLAS_SIDE sides[2] = {CblasLeft, CblasRight};
  enum CBLAS_UPLO uplos[2] = {CblasLower, CblasUpper};

  printf("==> L3 GEMM\n");
  for (int s = 0; s < NS_L3; s++) {
    int n = SIZES_L3[s];
    for (int ia = 0; ia < 2; ia++)
      for (int ib = 0; ib < 2; ib++) {
        check_sgemm_one(tr[ia], tr[ib], n, n, n);
        check_dgemm_one(tr[ia], tr[ib], n, n, n);
        if (n <= 36) {
          int m2 = n + 1;
          int n2 = n > 1 ? n - 1 : n;
          int k2 = n + 2;
          check_sgemm_one(tr[ia], tr[ib], m2, n2, k2);
          check_dgemm_one(tr[ia], tr[ib], m2, n2, k2);
        }
      }
  }

  printf("==> L3 CGEMM/ZGEMM\n");
  for (int s = 0; s < NS_CZ; s++) {
    int n = SIZES_CZ[s];
    for (int ia = 0; ia < 3; ia++)
      for (int ib = 0; ib < 3; ib++) {
        check_cgemm_one(ctr[ia], ctr[ib], n, n, n);
        check_zgemm_one(ctr[ia], ctr[ib], n, n, n);
      }
  }

  printf("==> L3 SYRK\n");
  for (int s = 0; s < NS_L3; s++) {
    int n = SIZES_L3[s];
    int k = n <= 36 ? n + 3 : n;
    for (int u = 0; u < 2; u++)
      for (int t = 0; t < 2; t++) {
        check_ssyrk_one(uplos[u], tr[t], n, k);
        check_dsyrk_one(uplos[u], tr[t], n, k);
      }
  }

  printf("==> L3 TRMM\n");
  for (int s = 0; s < NS_L3; s++) {
    int n = SIZES_L3[s];
    for (int si = 0; si < 2; si++)
      for (int u = 0; u < 2; u++)
        for (int t = 0; t < 2; t++) {
          check_strmm_one(sides[si], uplos[u], tr[t], n, n);
          check_dtrmm_one(sides[si], uplos[u], tr[t], n, n);
          if (n <= 36) {
            int m2 = n + 1;
            int n2 = n > 1 ? n - 1 : 1;
            check_strmm_one(sides[si], uplos[u], tr[t], m2, n2);
            check_dtrmm_one(sides[si], uplos[u], tr[t], m2, n2);
          }
        }
  }

  printf("==> L3 TRSM\n");
  for (int s = 0; s < NS_L3; s++) {
    int n = SIZES_L3[s];
    for (int si = 0; si < 2; si++)
      for (int u = 0; u < 2; u++)
        for (int t = 0; t < 2; t++) {
          check_strsm_one(sides[si], uplos[u], tr[t], n, n);
          check_dtrsm_one(sides[si], uplos[u], tr[t], n, n);
        }
  }
}
