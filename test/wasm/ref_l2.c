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

/* Typed GEMV oracles for the deep shape / transpose / stride grid in
 * check_l2.c. */

void ref_sgemv(enum CBLAS_TRANSPOSE trans, int m, int n, float alpha,
               const float *A, int lda, const float *x, int incx, float beta,
               float *y, int incy) {
  int leny = (trans == CblasNoTrans) ? m : n;
  int lenx = (trans == CblasNoTrans) ? n : m;
  for (int i = 0; i < leny; i++) {
    int yi = i * incy;
    float s = 0.0f;
    for (int j = 0; j < lenx; j++) {
      int xj = j * incx;
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
    int yi = i * incy;
    double s = 0.0;
    for (int j = 0; j < lenx; j++) {
      int xj = j * incx;
      double a = (trans == CblasNoTrans) ? A[i + j * lda] : A[j + i * lda];
      s += a * x[xj];
    }
    if (beta == 0.0)
      y[yi] = alpha * s;
    else
      y[yi] = alpha * s + beta * y[yi];
  }
}

static double vget(const void *p, int i, int is_double) {
  return is_double ? ((const double *)p)[i] : ((const float *)p)[i];
}

static void vset(void *p, int i, int is_double, double x) {
  if (is_double)
    ((double *)p)[i] = x;
  else
    ((float *)p)[i] = (float)x;
}

/* Complex multiply: (ar + i*ai) * (br + i*bi). Real callers pass ai=bi=0. */
static void cmul(double ar, double ai, double br, double bi, double *cr,
                 double *ci) {
  *cr = ar * br - ai * bi;
  *ci = ar * bi + ai * br;
}

/*
 * Dense matrix-vector product for the compact L2 suite:
 *   y := alpha * op(A) * x + beta * y
 * with A n-by-n, unit strides, stored column-major.
 *
 * width      1 = real, 2 = complex (interleaved re,im)
 * conjugate  0 = NoTrans (A[i,j]), 1 = ConjTrans (conj(A[j,i]))
 */
void ref_l2_mv(int n, const void *a, const void *x, void *y,
               const void *alpha, const void *beta, int width, int is_double,
               int conjugate) {
  double alpha_r = vget(alpha, 0, is_double);
  double alpha_i = (width == 2) ? vget(alpha, 1, is_double) : 0.0;
  double beta_r = vget(beta, 0, is_double);
  double beta_i = (width == 2) ? vget(beta, 1, is_double) : 0.0;

  for (int i = 0; i < n; i++) {
    double sum_r = 0.0;
    double sum_i = 0.0;

    for (int j = 0; j < n; j++) {
      int ia = width * (conjugate ? (j + i * n) : (i + j * n));
      int ix = width * j;
      double a_r = vget(a, ia, is_double);
      double a_i = (width == 2) ? vget(a, ia + 1, is_double) : 0.0;
      double x_r = vget(x, ix, is_double);
      double x_i = (width == 2) ? vget(x, ix + 1, is_double) : 0.0;
      double prod_r, prod_i;

      if (conjugate)
        a_i = -a_i;

      cmul(a_r, a_i, x_r, x_i, &prod_r, &prod_i);
      sum_r += prod_r;
      sum_i += prod_i;
    }

    int iy = width * i;
    double y_r = vget(y, iy, is_double);
    double y_i = (width == 2) ? vget(y, iy + 1, is_double) : 0.0;
    double ax_r, ax_i, by_r, by_i;

    cmul(alpha_r, alpha_i, sum_r, sum_i, &ax_r, &ax_i);
    cmul(beta_r, beta_i, y_r, y_i, &by_r, &by_i);
    vset(y, iy, is_double, ax_r + by_r);
    if (width == 2)
      vset(y, iy + 1, is_double, ax_i + by_i);
  }
}

/*
 * Rank-1 update for the compact L2 suite:
 *   A := A + alpha * x * y^T   (or y^H if conjugate_y)
 *
 * If symmetric or hermitian is set, only the lower triangle (i >= j) is
 * written. Hermitian updates force a zero imaginary part on the diagonal.
 */
void ref_l2_rank(int n, void *a, const void *x, const void *y,
                 const void *alpha, int width, int is_double, int conjugate_y,
                 int symmetric, int hermitian) {
  double alpha_r = vget(alpha, 0, is_double);
  double alpha_i = (width == 2) ? vget(alpha, 1, is_double) : 0.0;

  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
      int ix, iy, ia;
      double x_r, x_i, y_r, y_i;
      double outer_r, outer_i, scaled_r, scaled_i;

      if ((symmetric || hermitian) && i < j)
        continue;

      ix = width * i;
      iy = width * j;
      ia = width * (i + j * n);

      x_r = vget(x, ix, is_double);
      x_i = (width == 2) ? vget(x, ix + 1, is_double) : 0.0;
      y_r = vget(y, iy, is_double);
      y_i = (width == 2) ? vget(y, iy + 1, is_double) : 0.0;
      if (conjugate_y)
        y_i = -y_i;

      cmul(x_r, x_i, y_r, y_i, &outer_r, &outer_i);
      cmul(alpha_r, alpha_i, outer_r, outer_i, &scaled_r, &scaled_i);

      vset(a, ia, is_double, vget(a, ia, is_double) + scaled_r);
      if (width == 2) {
        double imag = (hermitian && i == j)
                          ? 0.0
                          : vget(a, ia + 1, is_double) + scaled_i;
        vset(a, ia + 1, is_double, imag);
      }
    }
  }
}
