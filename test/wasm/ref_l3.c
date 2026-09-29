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
 * Scalar IEEE Level-3 oracle for test/wasm.
 *
 * Plain triple loops using only IEEE `*` and `+`, so the result is
 * independent of the SIMD, FMA and blocking choices made by the kernels under
 * test. Matrices are column major and complex ones are interleaved re,im.
 *
 * The symmetric, Hermitian and triangular routines read only the triangle
 * that UPLO declares as stored, and the rank-k / rank-2k routines write only
 * that triangle, so callers do not have to mirror anything before or after.
 */

#include "ref.h"

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------------- */
/* Shared helpers                                                            */
/* ------------------------------------------------------------------------- */

/* Does (i,j) fall in the triangle that `uplo` declares stored? */
static int in_uplo(enum CBLAS_UPLO uplo, int i, int j) {
  return (uplo == CblasUpper) ? (i <= j) : (i >= j);
}

static void cmul_f32(float ar, float ai, float br, float bi, float *cr,
                     float *ci) {
  *cr = ar * br - ai * bi;
  *ci = ar * bi + ai * br;
}

static void cmul_f64(double ar, double ai, double br, double bi, double *cr,
                     double *ci) {
  *cr = ar * br - ai * bi;
  *ci = ar * bi + ai * br;
}

/* c := alpha * s + beta * c on one interleaved complex element. */
static void caccum_f32(const float *alpha, float sr, float si,
                       const float *beta, float *c) {
  float ar, ai;

  cmul_f32(alpha[0], alpha[1], sr, si, &ar, &ai);
  if (beta[0] == 0.0f && beta[1] == 0.0f) {
    c[0] = ar;
    c[1] = ai;
  } else {
    float br, bi;
    cmul_f32(beta[0], beta[1], c[0], c[1], &br, &bi);
    c[0] = ar + br;
    c[1] = ai + bi;
  }
}

static void caccum_f64(const double *alpha, double sr, double si,
                       const double *beta, double *c) {
  double ar, ai;

  cmul_f64(alpha[0], alpha[1], sr, si, &ar, &ai);
  if (beta[0] == 0.0 && beta[1] == 0.0) {
    c[0] = ar;
    c[1] = ai;
  } else {
    double br, bi;
    cmul_f64(beta[0], beta[1], c[0], c[1], &br, &bi);
    c[0] = ar + br;
    c[1] = ai + bi;
  }
}

static void cget_f32(const float *A, int i, int j, int lda,
                     enum CBLAS_TRANSPOSE t, float *r, float *im) {
  int idx;
  if (t == CblasNoTrans)
    idx = 2 * (i + j * lda);
  else if (t == CblasTrans)
    idx = 2 * (j + i * lda);
  else /* ConjTrans */
    idx = 2 * (j + i * lda);
  *r = A[idx];
  *im = A[idx + 1];
  if (t == CblasConjTrans)
    *im = -*im;
}

static void cget_f64(const double *A, int i, int j, int lda,
                     enum CBLAS_TRANSPOSE t, double *r, double *im) {
  int idx;
  if (t == CblasNoTrans)
    idx = 2 * (i + j * lda);
  else if (t == CblasTrans)
    idx = 2 * (j + i * lda);
  else
    idx = 2 * (j + i * lda);
  *r = A[idx];
  *im = A[idx + 1];
  if (t == CblasConjTrans)
    *im = -*im;
}

/* ------------------------------------------------------------------------- */
/* GEMM                                                                      */
/* ------------------------------------------------------------------------- */

void ref_sgemm(enum CBLAS_TRANSPOSE ta, enum CBLAS_TRANSPOSE tb, int m, int n,
               int k, float alpha, const float *A, int lda, const float *B,
               int ldb, float beta, float *C, int ldc) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++) {
      float s = 0.0f;
      for (int p = 0; p < k; p++) {
        float av = (ta == CblasNoTrans) ? A[i + p * lda] : A[p + i * lda];
        float bv = (tb == CblasNoTrans) ? B[p + j * ldb] : B[j + p * ldb];
        s += av * bv;
      }
      if (beta == 0.0f)
        C[i + j * ldc] = alpha * s;
      else
        C[i + j * ldc] = alpha * s + beta * C[i + j * ldc];
    }
}

void ref_dgemm(enum CBLAS_TRANSPOSE ta, enum CBLAS_TRANSPOSE tb, int m, int n,
               int k, double alpha, const double *A, int lda, const double *B,
               int ldb, double beta, double *C, int ldc) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++) {
      double s = 0.0;
      for (int p = 0; p < k; p++) {
        double av = (ta == CblasNoTrans) ? A[i + p * lda] : A[p + i * lda];
        double bv = (tb == CblasNoTrans) ? B[p + j * ldb] : B[j + p * ldb];
        s += av * bv;
      }
      if (beta == 0.0)
        C[i + j * ldc] = alpha * s;
      else
        C[i + j * ldc] = alpha * s + beta * C[i + j * ldc];
    }
}

void ref_cgemm(enum CBLAS_TRANSPOSE ta, enum CBLAS_TRANSPOSE tb, int m, int n,
               int k, const float *alpha, const float *A, int lda,
               const float *B, int ldb, const float *beta, float *C, int ldc) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++) {
      float sr = 0.0f, si = 0.0f;
      for (int p = 0; p < k; p++) {
        float ar, ai, br, bi, pr, pi;
        cget_f32(A, i, p, lda, ta, &ar, &ai);
        cget_f32(B, p, j, ldb, tb, &br, &bi);
        cmul_f32(ar, ai, br, bi, &pr, &pi);
        sr += pr;
        si += pi;
      }
      caccum_f32(alpha, sr, si, beta, &C[2 * (i + j * ldc)]);
    }
}

void ref_zgemm(enum CBLAS_TRANSPOSE ta, enum CBLAS_TRANSPOSE tb, int m, int n,
               int k, const double *alpha, const double *A, int lda,
               const double *B, int ldb, const double *beta, double *C,
               int ldc) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++) {
      double sr = 0.0, si = 0.0;
      for (int p = 0; p < k; p++) {
        double ar, ai, br, bi, pr, pi;
        cget_f64(A, i, p, lda, ta, &ar, &ai);
        cget_f64(B, p, j, ldb, tb, &br, &bi);
        cmul_f64(ar, ai, br, bi, &pr, &pi);
        sr += pr;
        si += pi;
      }
      caccum_f64(alpha, sr, si, beta, &C[2 * (i + j * ldc)]);
    }
}

/* ------------------------------------------------------------------------- */
/* SYMM / HEMM                                                               */
/* ------------------------------------------------------------------------- */

/* A[i,j] of a symmetric matrix, taken from the stored triangle only. */
static float sym_f32(const float *A, int lda, enum CBLAS_UPLO uplo, int i,
                     int j) {
  return in_uplo(uplo, i, j) ? A[i + j * lda] : A[j + i * lda];
}

static double sym_f64(const double *A, int lda, enum CBLAS_UPLO uplo, int i,
                      int j) {
  return in_uplo(uplo, i, j) ? A[i + j * lda] : A[j + i * lda];
}

/* Complex symmetric (A = A^T, no conjugation) element read. */
static void csym_f32(const float *A, int lda, enum CBLAS_UPLO uplo, int i,
                     int j, float *re, float *im) {
  int idx = in_uplo(uplo, i, j) ? 2 * (i + j * lda) : 2 * (j + i * lda);
  *re = A[idx];
  *im = A[idx + 1];
}

static void csym_f64(const double *A, int lda, enum CBLAS_UPLO uplo, int i,
                     int j, double *re, double *im) {
  int idx = in_uplo(uplo, i, j) ? 2 * (i + j * lda) : 2 * (j + i * lda);
  *re = A[idx];
  *im = A[idx + 1];
}

/* Hermitian (A = A^H) element read: the diagonal is real by definition and
 * the unstored triangle is the conjugate of the stored one. */
static void herm_f32(const float *A, int lda, enum CBLAS_UPLO uplo, int i,
                     int j, float *re, float *im) {
  if (i == j) {
    *re = A[2 * (i + i * lda)];
    *im = 0.0f;
  } else if (in_uplo(uplo, i, j)) {
    *re = A[2 * (i + j * lda)];
    *im = A[2 * (i + j * lda) + 1];
  } else {
    *re = A[2 * (j + i * lda)];
    *im = -A[2 * (j + i * lda) + 1];
  }
}

static void herm_f64(const double *A, int lda, enum CBLAS_UPLO uplo, int i,
                     int j, double *re, double *im) {
  if (i == j) {
    *re = A[2 * (i + i * lda)];
    *im = 0.0;
  } else if (in_uplo(uplo, i, j)) {
    *re = A[2 * (i + j * lda)];
    *im = A[2 * (i + j * lda) + 1];
  } else {
    *re = A[2 * (j + i * lda)];
    *im = -A[2 * (j + i * lda) + 1];
  }
}

void ref_ssymm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo, int m, int n,
               float alpha, const float *A, int lda, const float *B, int ldb,
               float beta, float *C, int ldc) {
  int ka = (side == CblasLeft) ? m : n;

  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++) {
      float s = 0.0f;
      for (int p = 0; p < ka; p++) {
        if (side == CblasLeft)
          s += sym_f32(A, lda, uplo, i, p) * B[p + j * ldb];
        else
          s += B[i + p * ldb] * sym_f32(A, lda, uplo, p, j);
      }
      if (beta == 0.0f)
        C[i + j * ldc] = alpha * s;
      else
        C[i + j * ldc] = alpha * s + beta * C[i + j * ldc];
    }
}

void ref_dsymm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo, int m, int n,
               double alpha, const double *A, int lda, const double *B,
               int ldb, double beta, double *C, int ldc) {
  int ka = (side == CblasLeft) ? m : n;

  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++) {
      double s = 0.0;
      for (int p = 0; p < ka; p++) {
        if (side == CblasLeft)
          s += sym_f64(A, lda, uplo, i, p) * B[p + j * ldb];
        else
          s += B[i + p * ldb] * sym_f64(A, lda, uplo, p, j);
      }
      if (beta == 0.0)
        C[i + j * ldc] = alpha * s;
      else
        C[i + j * ldc] = alpha * s + beta * C[i + j * ldc];
    }
}

void ref_csymm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo, int m, int n,
               const float *alpha, const float *A, int lda, const float *B,
               int ldb, const float *beta, float *C, int ldc) {
  int ka = (side == CblasLeft) ? m : n;

  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++) {
      float sr = 0.0f, si = 0.0f;
      for (int p = 0; p < ka; p++) {
        float ar, ai, br, bi, pr, pi;
        int bidx;
        if (side == CblasLeft) {
          csym_f32(A, lda, uplo, i, p, &ar, &ai);
          bidx = 2 * (p + j * ldb);
        } else {
          csym_f32(A, lda, uplo, p, j, &ar, &ai);
          bidx = 2 * (i + p * ldb);
        }
        br = B[bidx];
        bi = B[bidx + 1];
        cmul_f32(ar, ai, br, bi, &pr, &pi);
        sr += pr;
        si += pi;
      }
      caccum_f32(alpha, sr, si, beta, &C[2 * (i + j * ldc)]);
    }
}

void ref_zsymm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo, int m, int n,
               const double *alpha, const double *A, int lda, const double *B,
               int ldb, const double *beta, double *C, int ldc) {
  int ka = (side == CblasLeft) ? m : n;

  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++) {
      double sr = 0.0, si = 0.0;
      for (int p = 0; p < ka; p++) {
        double ar, ai, br, bi, pr, pi;
        int bidx;
        if (side == CblasLeft) {
          csym_f64(A, lda, uplo, i, p, &ar, &ai);
          bidx = 2 * (p + j * ldb);
        } else {
          csym_f64(A, lda, uplo, p, j, &ar, &ai);
          bidx = 2 * (i + p * ldb);
        }
        br = B[bidx];
        bi = B[bidx + 1];
        cmul_f64(ar, ai, br, bi, &pr, &pi);
        sr += pr;
        si += pi;
      }
      caccum_f64(alpha, sr, si, beta, &C[2 * (i + j * ldc)]);
    }
}

void ref_chemm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo, int m, int n,
               const float *alpha, const float *A, int lda, const float *B,
               int ldb, const float *beta, float *C, int ldc) {
  int ka = (side == CblasLeft) ? m : n;

  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++) {
      float sr = 0.0f, si = 0.0f;
      for (int p = 0; p < ka; p++) {
        float ar, ai, br, bi, pr, pi;
        int bidx;
        if (side == CblasLeft) {
          herm_f32(A, lda, uplo, i, p, &ar, &ai);
          bidx = 2 * (p + j * ldb);
        } else {
          herm_f32(A, lda, uplo, p, j, &ar, &ai);
          bidx = 2 * (i + p * ldb);
        }
        br = B[bidx];
        bi = B[bidx + 1];
        cmul_f32(ar, ai, br, bi, &pr, &pi);
        sr += pr;
        si += pi;
      }
      caccum_f32(alpha, sr, si, beta, &C[2 * (i + j * ldc)]);
    }
}

void ref_zhemm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo, int m, int n,
               const double *alpha, const double *A, int lda, const double *B,
               int ldb, const double *beta, double *C, int ldc) {
  int ka = (side == CblasLeft) ? m : n;

  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++) {
      double sr = 0.0, si = 0.0;
      for (int p = 0; p < ka; p++) {
        double ar, ai, br, bi, pr, pi;
        int bidx;
        if (side == CblasLeft) {
          herm_f64(A, lda, uplo, i, p, &ar, &ai);
          bidx = 2 * (p + j * ldb);
        } else {
          herm_f64(A, lda, uplo, p, j, &ar, &ai);
          bidx = 2 * (i + p * ldb);
        }
        br = B[bidx];
        bi = B[bidx + 1];
        cmul_f64(ar, ai, br, bi, &pr, &pi);
        sr += pr;
        si += pi;
      }
      caccum_f64(alpha, sr, si, beta, &C[2 * (i + j * ldc)]);
    }
}

/* ------------------------------------------------------------------------- */
/* SYRK / HERK                                                               */
/* ------------------------------------------------------------------------- */

void ref_ssyrk(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
               float alpha, const float *A, int lda, float beta, float *C,
               int ldc) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      if (!in_uplo(uplo, i, j))
        continue;
      float s = 0.0f;
      for (int p = 0; p < k; p++) {
        float ai = (trans == CblasNoTrans) ? A[i + p * lda] : A[p + i * lda];
        float aj = (trans == CblasNoTrans) ? A[j + p * lda] : A[p + j * lda];
        s += ai * aj;
      }
      if (beta == 0.0f)
        C[i + j * ldc] = alpha * s;
      else
        C[i + j * ldc] = alpha * s + beta * C[i + j * ldc];
    }
}

void ref_dsyrk(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
               double alpha, const double *A, int lda, double beta, double *C,
               int ldc) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      if (!in_uplo(uplo, i, j))
        continue;
      double s = 0.0;
      for (int p = 0; p < k; p++) {
        double ai = (trans == CblasNoTrans) ? A[i + p * lda] : A[p + i * lda];
        double aj = (trans == CblasNoTrans) ? A[j + p * lda] : A[p + j * lda];
        s += ai * aj;
      }
      if (beta == 0.0)
        C[i + j * ldc] = alpha * s;
      else
        C[i + j * ldc] = alpha * s + beta * C[i + j * ldc];
    }
}

void ref_csyrk(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
               const float *alpha, const float *A, int lda, const float *beta,
               float *C, int ldc) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      if (!in_uplo(uplo, i, j))
        continue;
      float sr = 0.0f, si = 0.0f;
      for (int p = 0; p < k; p++) {
        int ia = (trans == CblasNoTrans) ? 2 * (i + p * lda) : 2 * (p + i * lda);
        int ja = (trans == CblasNoTrans) ? 2 * (j + p * lda) : 2 * (p + j * lda);
        float pr, pi;
        cmul_f32(A[ia], A[ia + 1], A[ja], A[ja + 1], &pr, &pi);
        sr += pr;
        si += pi;
      }
      caccum_f32(alpha, sr, si, beta, &C[2 * (i + j * ldc)]);
    }
}

void ref_zsyrk(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
               const double *alpha, const double *A, int lda,
               const double *beta, double *C, int ldc) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      if (!in_uplo(uplo, i, j))
        continue;
      double sr = 0.0, si = 0.0;
      for (int p = 0; p < k; p++) {
        int ia = (trans == CblasNoTrans) ? 2 * (i + p * lda) : 2 * (p + i * lda);
        int ja = (trans == CblasNoTrans) ? 2 * (j + p * lda) : 2 * (p + j * lda);
        double pr, pi;
        cmul_f64(A[ia], A[ia + 1], A[ja], A[ja + 1], &pr, &pi);
        sr += pr;
        si += pi;
      }
      caccum_f64(alpha, sr, si, beta, &C[2 * (i + j * ldc)]);
    }
}

void ref_cherk(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
               float alpha, const float *A, int lda, float beta, float *C,
               int ldc) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      if (!in_uplo(uplo, i, j))
        continue;
      float sr = 0.0f, si = 0.0f;
      for (int p = 0; p < k; p++) {
        float xr, xi, yr, yi, pr, pi;
        if (trans == CblasNoTrans) {
          /* A * A^H */
          xr = A[2 * (i + p * lda)];
          xi = A[2 * (i + p * lda) + 1];
          yr = A[2 * (j + p * lda)];
          yi = -A[2 * (j + p * lda) + 1];
        } else {
          /* A^H * A */
          xr = A[2 * (p + i * lda)];
          xi = -A[2 * (p + i * lda) + 1];
          yr = A[2 * (p + j * lda)];
          yi = A[2 * (p + j * lda) + 1];
        }
        cmul_f32(xr, xi, yr, yi, &pr, &pi);
        sr += pr;
        si += pi;
      }
      int cidx = 2 * (i + j * ldc);
      if (beta == 0.0f) {
        C[cidx] = alpha * sr;
        C[cidx + 1] = alpha * si;
      } else {
        C[cidx] = alpha * sr + beta * C[cidx];
        C[cidx + 1] = alpha * si + beta * C[cidx + 1];
      }
      if (i == j)
        C[cidx + 1] = 0.0f;
    }
}

void ref_zherk(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
               double alpha, const double *A, int lda, double beta, double *C,
               int ldc) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      if (!in_uplo(uplo, i, j))
        continue;
      double sr = 0.0, si = 0.0;
      for (int p = 0; p < k; p++) {
        double xr, xi, yr, yi, pr, pi;
        if (trans == CblasNoTrans) {
          xr = A[2 * (i + p * lda)];
          xi = A[2 * (i + p * lda) + 1];
          yr = A[2 * (j + p * lda)];
          yi = -A[2 * (j + p * lda) + 1];
        } else {
          xr = A[2 * (p + i * lda)];
          xi = -A[2 * (p + i * lda) + 1];
          yr = A[2 * (p + j * lda)];
          yi = A[2 * (p + j * lda) + 1];
        }
        cmul_f64(xr, xi, yr, yi, &pr, &pi);
        sr += pr;
        si += pi;
      }
      int cidx = 2 * (i + j * ldc);
      if (beta == 0.0) {
        C[cidx] = alpha * sr;
        C[cidx + 1] = alpha * si;
      } else {
        C[cidx] = alpha * sr + beta * C[cidx];
        C[cidx + 1] = alpha * si + beta * C[cidx + 1];
      }
      if (i == j)
        C[cidx + 1] = 0.0;
    }
}

/* ------------------------------------------------------------------------- */
/* SYR2K / HER2K                                                             */
/* ------------------------------------------------------------------------- */

void ref_ssyr2k(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
                float alpha, const float *A, int lda, const float *B, int ldb,
                float beta, float *C, int ldc) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      if (!in_uplo(uplo, i, j))
        continue;
      float s = 0.0f;
      for (int p = 0; p < k; p++) {
        int no = (trans == CblasNoTrans);
        float aip = no ? A[i + p * lda] : A[p + i * lda];
        float ajp = no ? A[j + p * lda] : A[p + j * lda];
        float bip = no ? B[i + p * ldb] : B[p + i * ldb];
        float bjp = no ? B[j + p * ldb] : B[p + j * ldb];
        s += aip * bjp + bip * ajp;
      }
      if (beta == 0.0f)
        C[i + j * ldc] = alpha * s;
      else
        C[i + j * ldc] = alpha * s + beta * C[i + j * ldc];
    }
}

void ref_dsyr2k(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
                double alpha, const double *A, int lda, const double *B,
                int ldb, double beta, double *C, int ldc) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      if (!in_uplo(uplo, i, j))
        continue;
      double s = 0.0;
      for (int p = 0; p < k; p++) {
        int no = (trans == CblasNoTrans);
        double aip = no ? A[i + p * lda] : A[p + i * lda];
        double ajp = no ? A[j + p * lda] : A[p + j * lda];
        double bip = no ? B[i + p * ldb] : B[p + i * ldb];
        double bjp = no ? B[j + p * ldb] : B[p + j * ldb];
        s += aip * bjp + bip * ajp;
      }
      if (beta == 0.0)
        C[i + j * ldc] = alpha * s;
      else
        C[i + j * ldc] = alpha * s + beta * C[i + j * ldc];
    }
}

void ref_csyr2k(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
                const float *alpha, const float *A, int lda, const float *B,
                int ldb, const float *beta, float *C, int ldc) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      if (!in_uplo(uplo, i, j))
        continue;
      float sr = 0.0f, si = 0.0f;
      for (int p = 0; p < k; p++) {
        int no = (trans == CblasNoTrans);
        int ia = no ? 2 * (i + p * lda) : 2 * (p + i * lda);
        int ja = no ? 2 * (j + p * lda) : 2 * (p + j * lda);
        int ib = no ? 2 * (i + p * ldb) : 2 * (p + i * ldb);
        int jb = no ? 2 * (j + p * ldb) : 2 * (p + j * ldb);
        float pr, pi;
        cmul_f32(A[ia], A[ia + 1], B[jb], B[jb + 1], &pr, &pi);
        sr += pr;
        si += pi;
        cmul_f32(B[ib], B[ib + 1], A[ja], A[ja + 1], &pr, &pi);
        sr += pr;
        si += pi;
      }
      caccum_f32(alpha, sr, si, beta, &C[2 * (i + j * ldc)]);
    }
}

void ref_zsyr2k(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
                const double *alpha, const double *A, int lda, const double *B,
                int ldb, const double *beta, double *C, int ldc) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      if (!in_uplo(uplo, i, j))
        continue;
      double sr = 0.0, si = 0.0;
      for (int p = 0; p < k; p++) {
        int no = (trans == CblasNoTrans);
        int ia = no ? 2 * (i + p * lda) : 2 * (p + i * lda);
        int ja = no ? 2 * (j + p * lda) : 2 * (p + j * lda);
        int ib = no ? 2 * (i + p * ldb) : 2 * (p + i * ldb);
        int jb = no ? 2 * (j + p * ldb) : 2 * (p + j * ldb);
        double pr, pi;
        cmul_f64(A[ia], A[ia + 1], B[jb], B[jb + 1], &pr, &pi);
        sr += pr;
        si += pi;
        cmul_f64(B[ib], B[ib + 1], A[ja], A[ja + 1], &pr, &pi);
        sr += pr;
        si += pi;
      }
      caccum_f64(alpha, sr, si, beta, &C[2 * (i + j * ldc)]);
    }
}

void ref_cher2k(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
                const float *alpha, const float *A, int lda, const float *B,
                int ldb, float beta, float *C, int ldc) {
  float conj_alpha[2] = {alpha[0], -alpha[1]};

  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      if (!in_uplo(uplo, i, j))
        continue;
      /* s1 accumulates the alpha term, s2 the conj(alpha) one. */
      float s1r = 0.0f, s1i = 0.0f, s2r = 0.0f, s2i = 0.0f;
      for (int p = 0; p < k; p++) {
        float xr, xi, yr, yi, pr, pi;
        if (trans == CblasNoTrans) {
          /* A * B^H and B * A^H */
          xr = A[2 * (i + p * lda)];
          xi = A[2 * (i + p * lda) + 1];
          yr = B[2 * (j + p * ldb)];
          yi = -B[2 * (j + p * ldb) + 1];
          cmul_f32(xr, xi, yr, yi, &pr, &pi);
          s1r += pr;
          s1i += pi;
          xr = B[2 * (i + p * ldb)];
          xi = B[2 * (i + p * ldb) + 1];
          yr = A[2 * (j + p * lda)];
          yi = -A[2 * (j + p * lda) + 1];
        } else {
          /* A^H * B and B^H * A */
          xr = A[2 * (p + i * lda)];
          xi = -A[2 * (p + i * lda) + 1];
          yr = B[2 * (p + j * ldb)];
          yi = B[2 * (p + j * ldb) + 1];
          cmul_f32(xr, xi, yr, yi, &pr, &pi);
          s1r += pr;
          s1i += pi;
          xr = B[2 * (p + i * ldb)];
          xi = -B[2 * (p + i * ldb) + 1];
          yr = A[2 * (p + j * lda)];
          yi = A[2 * (p + j * lda) + 1];
        }
        cmul_f32(xr, xi, yr, yi, &pr, &pi);
        s2r += pr;
        s2i += pi;
      }
      float t1r, t1i, t2r, t2i;
      cmul_f32(alpha[0], alpha[1], s1r, s1i, &t1r, &t1i);
      cmul_f32(conj_alpha[0], conj_alpha[1], s2r, s2i, &t2r, &t2i);
      int cidx = 2 * (i + j * ldc);
      if (beta == 0.0f) {
        C[cidx] = t1r + t2r;
        C[cidx + 1] = t1i + t2i;
      } else {
        C[cidx] = t1r + t2r + beta * C[cidx];
        C[cidx + 1] = t1i + t2i + beta * C[cidx + 1];
      }
      if (i == j)
        C[cidx + 1] = 0.0f;
    }
}

void ref_zher2k(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
                const double *alpha, const double *A, int lda,
                const double *B, int ldb, double beta, double *C, int ldc) {
  double conj_alpha[2] = {alpha[0], -alpha[1]};

  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      if (!in_uplo(uplo, i, j))
        continue;
      double s1r = 0.0, s1i = 0.0, s2r = 0.0, s2i = 0.0;
      for (int p = 0; p < k; p++) {
        double xr, xi, yr, yi, pr, pi;
        if (trans == CblasNoTrans) {
          xr = A[2 * (i + p * lda)];
          xi = A[2 * (i + p * lda) + 1];
          yr = B[2 * (j + p * ldb)];
          yi = -B[2 * (j + p * ldb) + 1];
          cmul_f64(xr, xi, yr, yi, &pr, &pi);
          s1r += pr;
          s1i += pi;
          xr = B[2 * (i + p * ldb)];
          xi = B[2 * (i + p * ldb) + 1];
          yr = A[2 * (j + p * lda)];
          yi = -A[2 * (j + p * lda) + 1];
        } else {
          xr = A[2 * (p + i * lda)];
          xi = -A[2 * (p + i * lda) + 1];
          yr = B[2 * (p + j * ldb)];
          yi = B[2 * (p + j * ldb) + 1];
          cmul_f64(xr, xi, yr, yi, &pr, &pi);
          s1r += pr;
          s1i += pi;
          xr = B[2 * (p + i * ldb)];
          xi = -B[2 * (p + i * ldb) + 1];
          yr = A[2 * (p + j * lda)];
          yi = A[2 * (p + j * lda) + 1];
        }
        cmul_f64(xr, xi, yr, yi, &pr, &pi);
        s2r += pr;
        s2i += pi;
      }
      double t1r, t1i, t2r, t2i;
      cmul_f64(alpha[0], alpha[1], s1r, s1i, &t1r, &t1i);
      cmul_f64(conj_alpha[0], conj_alpha[1], s2r, s2i, &t2r, &t2i);
      int cidx = 2 * (i + j * ldc);
      if (beta == 0.0) {
        C[cidx] = t1r + t2r;
        C[cidx + 1] = t1i + t2i;
      } else {
        C[cidx] = t1r + t2r + beta * C[cidx];
        C[cidx + 1] = t1i + t2i + beta * C[cidx + 1];
      }
      if (i == j)
        C[cidx + 1] = 0.0;
    }
}

/* ------------------------------------------------------------------------- */
/* TRMM / TRSM                                                               */
/* ------------------------------------------------------------------------- */

static float tri_a_f32(const float *A, int lda, enum CBLAS_UPLO uplo,
                       enum CBLAS_TRANSPOSE t, enum CBLAS_DIAG diag, int i,
                       int j) {
  int ii = i, jj = j;
  if (t != CblasNoTrans) {
    ii = j;
    jj = i;
  }
  if (ii == jj)
    return (diag == CblasUnit) ? 1.0f : A[ii + jj * lda];
  if (uplo == CblasLower) {
    if (ii > jj)
      return A[ii + jj * lda];
    return 0.0f;
  }
  if (ii < jj)
    return A[ii + jj * lda];
  return 0.0f;
}

static double tri_a_f64(const double *A, int lda, enum CBLAS_UPLO uplo,
                        enum CBLAS_TRANSPOSE t, enum CBLAS_DIAG diag, int i,
                        int j) {
  int ii = i, jj = j;
  if (t != CblasNoTrans) {
    ii = j;
    jj = i;
  }
  if (ii == jj)
    return (diag == CblasUnit) ? 1.0 : A[ii + jj * lda];
  if (uplo == CblasLower) {
    if (ii > jj)
      return A[ii + jj * lda];
    return 0.0;
  }
  if (ii < jj)
    return A[ii + jj * lda];
  return 0.0;
}

void ref_strmm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
               enum CBLAS_TRANSPOSE trans, enum CBLAS_DIAG diag, int m, int n,
               float alpha, const float *A, int lda, float *B, int ldb) {
  float *T = (float *)malloc((size_t)m * (size_t)n * sizeof(float));
  if (!T)
    exit(1);
  memcpy(T, B, (size_t)m * (size_t)n * sizeof(float));
  if (side == CblasLeft) {
    int ka = m;
    for (int j = 0; j < n; j++)
      for (int i = 0; i < m; i++) {
        float s = 0.0f;
        for (int p = 0; p < ka; p++)
          s += tri_a_f32(A, lda, uplo, trans, diag, i, p) * T[p + j * ldb];
        B[i + j * ldb] = alpha * s;
      }
  } else {
    int ka = n;
    for (int j = 0; j < n; j++)
      for (int i = 0; i < m; i++) {
        float s = 0.0f;
        for (int p = 0; p < ka; p++)
          s += T[i + p * ldb] * tri_a_f32(A, lda, uplo, trans, diag, p, j);
        B[i + j * ldb] = alpha * s;
      }
  }
  free(T);
}

void ref_dtrmm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
               enum CBLAS_TRANSPOSE trans, enum CBLAS_DIAG diag, int m, int n,
               double alpha, const double *A, int lda, double *B, int ldb) {
  double *T = (double *)malloc((size_t)m * (size_t)n * sizeof(double));
  if (!T)
    exit(1);
  memcpy(T, B, (size_t)m * (size_t)n * sizeof(double));
  if (side == CblasLeft) {
    int ka = m;
    for (int j = 0; j < n; j++)
      for (int i = 0; i < m; i++) {
        double s = 0.0;
        for (int p = 0; p < ka; p++)
          s += tri_a_f64(A, lda, uplo, trans, diag, i, p) * T[p + j * ldb];
        B[i + j * ldb] = alpha * s;
      }
  } else {
    int ka = n;
    for (int j = 0; j < n; j++)
      for (int i = 0; i < m; i++) {
        double s = 0.0;
        for (int p = 0; p < ka; p++)
          s += T[i + p * ldb] * tri_a_f64(A, lda, uplo, trans, diag, p, j);
        B[i + j * ldb] = alpha * s;
      }
  }
  free(T);
}

/* TRSM by forward / back substitution on the triangle that op(A) exposes;
 * transposing swaps which of the two directions applies. */

static void trsm_left_f32(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans,
                          enum CBLAS_DIAG diag, int m, int n, float alpha,
                          const float *A, int lda, float *B, int ldb) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++)
      B[i + j * ldb] *= alpha;

  int lower = (uplo == CblasLower) ^ (trans != CblasNoTrans);
  if (lower) {
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < m; i++) {
        float s = B[i + j * ldb];
        for (int k = 0; k < i; k++)
          s -= tri_a_f32(A, lda, uplo, trans, diag, i, k) * B[k + j * ldb];
        float d = tri_a_f32(A, lda, uplo, trans, diag, i, i);
        B[i + j * ldb] = s / d;
      }
    }
  } else {
    for (int j = 0; j < n; j++) {
      for (int i = m - 1; i >= 0; i--) {
        float s = B[i + j * ldb];
        for (int k = i + 1; k < m; k++)
          s -= tri_a_f32(A, lda, uplo, trans, diag, i, k) * B[k + j * ldb];
        float d = tri_a_f32(A, lda, uplo, trans, diag, i, i);
        B[i + j * ldb] = s / d;
      }
    }
  }
}

static void trsm_right_f32(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans,
                           enum CBLAS_DIAG diag, int m, int n, float alpha,
                           const float *A, int lda, float *B, int ldb) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++)
      B[i + j * ldb] *= alpha;

  int lower = (uplo == CblasLower) ^ (trans != CblasNoTrans);
  if (lower) {
    for (int j = n - 1; j >= 0; j--) {
      for (int i = 0; i < m; i++) {
        float s = B[i + j * ldb];
        for (int k = j + 1; k < n; k++)
          s -= B[i + k * ldb] * tri_a_f32(A, lda, uplo, trans, diag, k, j);
        float d = tri_a_f32(A, lda, uplo, trans, diag, j, j);
        B[i + j * ldb] = s / d;
      }
    }
  } else {
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < m; i++) {
        float s = B[i + j * ldb];
        for (int k = 0; k < j; k++)
          s -= B[i + k * ldb] * tri_a_f32(A, lda, uplo, trans, diag, k, j);
        float d = tri_a_f32(A, lda, uplo, trans, diag, j, j);
        B[i + j * ldb] = s / d;
      }
    }
  }
}

static void trsm_left_f64(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans,
                          enum CBLAS_DIAG diag, int m, int n, double alpha,
                          const double *A, int lda, double *B, int ldb) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++)
      B[i + j * ldb] *= alpha;

  int lower = (uplo == CblasLower) ^ (trans != CblasNoTrans);
  if (lower) {
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < m; i++) {
        double s = B[i + j * ldb];
        for (int k = 0; k < i; k++)
          s -= tri_a_f64(A, lda, uplo, trans, diag, i, k) * B[k + j * ldb];
        double d = tri_a_f64(A, lda, uplo, trans, diag, i, i);
        B[i + j * ldb] = s / d;
      }
    }
  } else {
    for (int j = 0; j < n; j++) {
      for (int i = m - 1; i >= 0; i--) {
        double s = B[i + j * ldb];
        for (int k = i + 1; k < m; k++)
          s -= tri_a_f64(A, lda, uplo, trans, diag, i, k) * B[k + j * ldb];
        double d = tri_a_f64(A, lda, uplo, trans, diag, i, i);
        B[i + j * ldb] = s / d;
      }
    }
  }
}

static void trsm_right_f64(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans,
                           enum CBLAS_DIAG diag, int m, int n, double alpha,
                           const double *A, int lda, double *B, int ldb) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++)
      B[i + j * ldb] *= alpha;

  int lower = (uplo == CblasLower) ^ (trans != CblasNoTrans);
  if (lower) {
    for (int j = n - 1; j >= 0; j--) {
      for (int i = 0; i < m; i++) {
        double s = B[i + j * ldb];
        for (int k = j + 1; k < n; k++)
          s -= B[i + k * ldb] * tri_a_f64(A, lda, uplo, trans, diag, k, j);
        double d = tri_a_f64(A, lda, uplo, trans, diag, j, j);
        B[i + j * ldb] = s / d;
      }
    }
  } else {
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < m; i++) {
        double s = B[i + j * ldb];
        for (int k = 0; k < j; k++)
          s -= B[i + k * ldb] * tri_a_f64(A, lda, uplo, trans, diag, k, j);
        double d = tri_a_f64(A, lda, uplo, trans, diag, j, j);
        B[i + j * ldb] = s / d;
      }
    }
  }
}

void ref_strsm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
               enum CBLAS_TRANSPOSE trans, enum CBLAS_DIAG diag, int m, int n,
               float alpha, const float *A, int lda, float *B, int ldb) {
  if (side == CblasLeft)
    trsm_left_f32(uplo, trans, diag, m, n, alpha, A, lda, B, ldb);
  else
    trsm_right_f32(uplo, trans, diag, m, n, alpha, A, lda, B, ldb);
}

void ref_dtrsm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
               enum CBLAS_TRANSPOSE trans, enum CBLAS_DIAG diag, int m, int n,
               double alpha, const double *A, int lda, double *B, int ldb) {
  if (side == CblasLeft)
    trsm_left_f64(uplo, trans, diag, m, n, alpha, A, lda, B, ldb);
  else
    trsm_right_f64(uplo, trans, diag, m, n, alpha, A, lda, B, ldb);
}
