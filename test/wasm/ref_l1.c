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

#include <math.h>

/* Typed AXPY oracles for the deep remainder / stride grid in check_l1.c. */

void ref_saxpy(int n, float alpha, const float *x, int incx, float *y,
               int incy) {
  if (n <= 0)
    return;
  if (incx == 0 && incy == 0) {
    y[0] += (float)n * alpha * x[0];
    return;
  }
  for (int i = 0; i < n; i++)
    y[i * incy] += alpha * x[i * incx];
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
    y[i * incy] += alpha * x[i * incx];
}

static double getv(const void *p, int i, int is_double) {
  return is_double ? ((const double *)p)[i] : ((const float *)p)[i];
}

static void setv(void *p, int i, int is_double, double v) {
  if (is_double)
    ((double *)p)[i] = v;
  else
    ((float *)p)[i] = (float)v;
}

void ref_l1_swap(int n, void *x, int incx, void *y, int incy, int width,
                 int is_double) {
  for (int i = 0; i < n; i++)
    for (int q = 0; q < width; q++) {
      int ix = width * i * incx + q, iy = width * i * incy + q;
      double t = getv(x, ix, is_double);
      setv(x, ix, is_double, getv(y, iy, is_double));
      setv(y, iy, is_double, t);
    }
}

void ref_l1_copy(int n, const void *x, int incx, void *y, int incy, int width,
                 int is_double) {
  for (int i = 0; i < n; i++)
    for (int q = 0; q < width; q++)
      setv(y, width * i * incy + q, is_double,
           getv(x, width * i * incx + q, is_double));
}

void ref_l1_scal(int n, const void *alpha, void *x, int incx, int width,
                 int is_double, int real_alpha) {
  double ar = getv(alpha, 0, is_double);
  double ai = real_alpha ? 0.0 : getv(alpha, 1, is_double);
  for (int i = 0; i < n; i++) {
    int k = width * i * incx;
    double xr = getv(x, k, is_double);
    if (width == 1)
      setv(x, k, is_double, ar * xr);
    else {
      double xi = getv(x, k + 1, is_double);
      setv(x, k, is_double, ar * xr - ai * xi);
      setv(x, k + 1, is_double, ar * xi + ai * xr);
    }
  }
}

void ref_l1_axpy(int n, const void *alpha, const void *x, int incx, void *y,
                 int incy, int width, int is_double) {
  double ar = getv(alpha, 0, is_double);
  double ai = width == 1 ? 0.0 : getv(alpha, 1, is_double);
  for (int i = 0; i < n; i++) {
    int ix = width * i * incx, iy = width * i * incy;
    double xr = getv(x, ix, is_double), yr = getv(y, iy, is_double);
    if (width == 1)
      setv(y, iy, is_double, yr + ar * xr);
    else {
      double xi = getv(x, ix + 1, is_double);
      double yi = getv(y, iy + 1, is_double);
      setv(y, iy, is_double, yr + ar * xr - ai * xi);
      setv(y, iy + 1, is_double, yi + ar * xi + ai * xr);
    }
  }
}

void ref_l1_dot(int n, const void *x, int incx, const void *y, int incy,
                void *out, int width, int is_double, int conjugate) {
  double sr = 0.0, si = 0.0;
  for (int i = 0; i < n; i++) {
    int ix = width * i * incx, iy = width * i * incy;
    double xr = getv(x, ix, is_double), yr = getv(y, iy, is_double);
    if (width == 1)
      sr += xr * yr;
    else {
      double xi = getv(x, ix + 1, is_double);
      double yi = getv(y, iy + 1, is_double);
      if (conjugate)
        xi = -xi;
      sr += xr * yr - xi * yi;
      si += xr * yi + xi * yr;
    }
  }
  setv(out, 0, is_double, sr);
  if (width == 2)
    setv(out, 1, is_double, si);
}

double ref_l1_nrm2(int n, const void *x, int incx, int width, int is_double) {
  double scale = 0.0, ssq = 1.0;
  for (int i = 0; i < n; i++)
    for (int q = 0; q < width; q++) {
      double a = fabs(getv(x, width * i * incx + q, is_double));
      if (a != 0.0) {
        if (scale < a) {
          double r = scale / a;
          ssq = 1.0 + ssq * r * r;
          scale = a;
        } else {
          double r = a / scale;
          ssq += r * r;
        }
      }
    }
  return scale == 0.0 ? 0.0 : scale * sqrt(ssq);
}

double ref_l1_asum(int n, const void *x, int incx, int width, int is_double) {
  double sum = 0.0;
  for (int i = 0; i < n; i++)
    for (int q = 0; q < width; q++)
      sum += fabs(getv(x, width * i * incx + q, is_double));
  return sum;
}

size_t ref_l1_iamax(int n, const void *x, int incx, int width, int is_double) {
  size_t best = 0;
  double vmax = -1.0;
  for (int i = 0; i < n; i++) {
    double v = 0.0;
    for (int q = 0; q < width; q++)
      v += fabs(getv(x, width * i * incx + q, is_double));
    if (v > vmax) {
      vmax = v;
      best = (size_t)i;
    }
  }
  return n > 0 ? best : 0;
}

void ref_l1_rot(int n, void *x, int incx, void *y, int incy, double c,
                double s, int width, int is_double) {
  for (int i = 0; i < n; i++)
    for (int q = 0; q < width; q++) {
      int ix = width * i * incx + q, iy = width * i * incy + q;
      double xv = getv(x, ix, is_double), yv = getv(y, iy, is_double);
      setv(x, ix, is_double, c * xv + s * yv);
      setv(y, iy, is_double, c * yv - s * xv);
    }
}

void ref_l1_rotm(int n, void *x, int incx, void *y, int incy,
                 const void *param, int is_double) {
  double flag = getv(param, 0, is_double);
  if (flag == -2.0)
    return;
  double h11 = (flag < 0.0 || flag == 0.0) ? getv(param, 1, is_double) : 1.0;
  double h21 = (flag < 0.0 || flag > 0.0) ? getv(param, 2, is_double) : -1.0;
  double h12 = (flag < 0.0 || flag > 0.0) ? getv(param, 3, is_double) : 1.0;
  double h22 = (flag < 0.0 || flag == 0.0) ? getv(param, 4, is_double) : 1.0;
  for (int i = 0; i < n; i++) {
    int ix = i * incx, iy = i * incy;
    double w = getv(x, ix, is_double), z = getv(y, iy, is_double);
    setv(x, ix, is_double, w * h11 + z * h12);
    setv(y, iy, is_double, w * h21 + z * h22);
  }
}

/*
 * Modified Givens generator (BLAS ROTMG).
 *
 * Updates d1, d2, b1 and writes param[0..4]:
 *   flag = param[0] selects which H entries are stored:
 *     -2  identity (early exit; no H written)
 *     -1  full H: h11, h21, h12, h22 in param[1..4]
 *      0  off-diagonals only: h21, h12 (diags implied 1)
 *      1  diagonals only: h11, h22 (off-diags implied -1 / 1)
 *   See Lawson et al. / Netlib drotmg for the case split on |q1| vs |q2|.
 */
void ref_l1_rotmg(void *d1p, void *d2p, void *b1p, const void *b2p,
                  void *param, int is_double) {
  double d1 = getv(d1p, 0, is_double);
  double d2 = getv(d2p, 0, is_double);
  double b1 = getv(b1p, 0, is_double);
  double b2 = getv(b2p, 0, is_double);
  double flag = -1.0;
  double h11 = 0.0, h12 = 0.0, h21 = 0.0, h22 = 0.0;

  if (d1 < 0.0) {
    /* Negative d1: zero the state and return a full (zero) H. */
    flag = -1.0;
    d1 = 0.0;
    d2 = 0.0;
    b1 = 0.0;
  } else {
    double p2 = d2 * b2;
    if (p2 == 0.0) {
      /* No second component: leave d1/d2/b1 unchanged, flag = -2. */
      setv(param, 0, is_double, -2.0);
      return;
    }

    double p1 = d1 * b1;
    double q1 = p1 * b1;
    double q2 = p2 * b2;

    if (fabs(q1) > fabs(q2)) {
      /* Prefer scaling that keeps |h21|,|h12| from exploding. */
      h21 = -b2 / b1;
      h12 = p2 / p1;
      double u = 1.0 - h12 * h21;
      if (u <= 0.0) {
        flag = -1.0;
        d1 = 0.0;
        d2 = 0.0;
        b1 = 0.0;
      } else {
        flag = 0.0;
        d1 /= u;
        d2 /= u;
        b1 *= u;
      }
    } else if (q2 < 0.0) {
      flag = -1.0;
      d1 = 0.0;
      d2 = 0.0;
      b1 = 0.0;
    } else {
      h11 = p1 / p2;
      h22 = b1 / b2;
      double u = 1.0 + h11 * h22;
      double tmp = d2 / u;
      d2 = d1 / u;
      d1 = tmp;
      b1 = b2 * u;
      flag = 1.0;
    }
  }

  setv(d1p, 0, is_double, d1);
  setv(d2p, 0, is_double, d2);
  setv(b1p, 0, is_double, b1);
  setv(param, 0, is_double, flag);

  if (flag < 0.0) {
    setv(param, 1, is_double, h11);
    setv(param, 2, is_double, h21);
    setv(param, 3, is_double, h12);
    setv(param, 4, is_double, h22);
  } else if (flag == 0.0) {
    setv(param, 2, is_double, h21);
    setv(param, 3, is_double, h12);
  } else if (flag == 1.0) {
    setv(param, 1, is_double, h11);
    setv(param, 4, is_double, h22);
  }
}
