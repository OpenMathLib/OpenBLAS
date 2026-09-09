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

#include "ref.h"

#include <stdlib.h>
#include <string.h>

void ref_saxpy(int n, float alpha, const float *x, int incx, float *y,
               int incy) {
  if (n <= 0)
    return;
  if (incx == 0 && incy == 0) {
    y[0] += (float)n * alpha * x[0];
    return;
  }
  for (int i = 0; i < n; i++)
    y[i * (incy == 0 ? 0 : incy)] += alpha * x[i * (incx == 0 ? 0 : incx)];
}

void ref_daxpy(int n, double alpha, const double *x, int incx, double *y,
               int incy) {
  if (n <= 0)
    return;
  if (incx == 0 && incy == 0) {
    y[0] += (double)n * alpha * x[0];
    return;
  }
  for (int i = 0; i < n; i++)
    y[i * (incy == 0 ? 0 : incy)] += alpha * x[i * (incx == 0 ? 0 : incx)];
}

void ref_sgemv(enum CBLAS_TRANSPOSE trans, int m, int n, float alpha,
               const float *A, int lda, const float *x, int incx, float beta,
               float *y, int incy) {
  int leny = (trans == CblasNoTrans) ? m : n;
  int lenx = (trans == CblasNoTrans) ? n : m;
  for (int i = 0; i < leny; i++) {
    int yi = i * (incy == 0 ? 0 : incy);
    float s = 0.0f;
    for (int j = 0; j < lenx; j++) {
      int xj = j * (incx == 0 ? 0 : incx);
      float a = (trans == CblasNoTrans) ? A[i + j * lda] : A[j + i * lda];
      s += a * x[xj];
    }
    if (beta == 0.0f)
      y[yi] = alpha * s;
    else
      y[yi] = alpha * s + beta * y[yi];
  }
}

void ref_dgemv(enum CBLAS_TRANSPOSE trans, int m, int n, double alpha,
               const double *A, int lda, const double *x, int incx, double beta,
               double *y, int incy) {
  int leny = (trans == CblasNoTrans) ? m : n;
  int lenx = (trans == CblasNoTrans) ? n : m;
  for (int i = 0; i < leny; i++) {
    int yi = i * (incy == 0 ? 0 : incy);
    double s = 0.0;
    for (int j = 0; j < lenx; j++) {
      int xj = j * (incx == 0 ? 0 : incx);
      double a = (trans == CblasNoTrans) ? A[i + j * lda] : A[j + i * lda];
      s += a * x[xj];
    }
    if (beta == 0.0)
      y[yi] = alpha * s;
    else
      y[yi] = alpha * s + beta * y[yi];
  }
}

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
      float cr, ci;
      cmul_f32(alpha[0], alpha[1], sr, si, &cr, &ci);
      int cidx = 2 * (i + j * ldc);
      if (beta[0] == 0.0f && beta[1] == 0.0f) {
        C[cidx] = cr;
        C[cidx + 1] = ci;
      } else {
        float br, bi;
        cmul_f32(beta[0], beta[1], C[cidx], C[cidx + 1], &br, &bi);
        C[cidx] = cr + br;
        C[cidx + 1] = ci + bi;
      }
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
      double cr, ci;
      cmul_f64(alpha[0], alpha[1], sr, si, &cr, &ci);
      int cidx = 2 * (i + j * ldc);
      if (beta[0] == 0.0 && beta[1] == 0.0) {
        C[cidx] = cr;
        C[cidx + 1] = ci;
      } else {
        double br, bi;
        cmul_f64(beta[0], beta[1], C[cidx], C[cidx + 1], &br, &bi);
        C[cidx] = cr + br;
        C[cidx + 1] = ci + bi;
      }
    }
}

void ref_ssyrk(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
               float alpha, const float *A, int lda, float beta, float *C,
               int ldc) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      int write = (uplo == CblasUpper) ? (i <= j) : (i >= j);
      if (!write)
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
      int write = (uplo == CblasUpper) ? (i <= j) : (i >= j);
      if (!write)
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

void make_tri_f32(float *A, int n, int lda, enum CBLAS_UPLO uplo, int unit) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      float v = 0.0f;
      if (i == j)
        v = unit ? 1.0f : (1.0f + 0.1f * (float)((i % 5) + 1));
      else if (uplo == CblasLower && i > j)
        v = 0.1f * (float)((i + j) % 5 + 1);
      else if (uplo == CblasUpper && i < j)
        v = 0.1f * (float)((i + j) % 5 + 1);
      A[i + j * lda] = v;
    }
}

void make_tri_f64(double *A, int n, int lda, enum CBLAS_UPLO uplo, int unit) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      double v = 0.0;
      if (i == j)
        v = unit ? 1.0 : (1.0 + 0.1 * (double)((i % 5) + 1));
      else if (uplo == CblasLower && i > j)
        v = 0.1 * (double)((i + j) % 5 + 1);
      else if (uplo == CblasUpper && i < j)
        v = 0.1 * (double)((i + j) % 5 + 1);
      A[i + j * lda] = v;
    }
}

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

/* TRSM via reconstruct: not used for reference solution; implement forward/back
 * sub for Left Lower NoTrans NonUnit as primary, else use gemm-style inverse
 * via iterative substitution on triangular. */

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
