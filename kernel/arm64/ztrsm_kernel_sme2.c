/***************************************************************************
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
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

/* Complex TRSM kernels for the SME2 complex GEMM kernel: the update in GEMM_KERNEL, the solve in plain or NEON code. */

#include "common.h"
#include <arm_neon.h>

#if defined(LN) || defined(LT)
#define ZT_LEFT 1
#endif
#if defined(LN) || defined(RT)
#define ZT_BACK 1
#endif
#ifdef CONJ
#define ZT_CJ 1
#else
#define ZT_CJ 0
#endif
#ifdef ZT_LEFT
#define ZT_GEMM (ZT_CJ ? GEMM_KERNEL_L : GEMM_KERNEL_N)
#else
#define ZT_GEMM (ZT_CJ ? GEMM_KERNEL_R : GEMM_KERNEL_N)
#endif

#ifdef DOUBLE
typedef float64x2_t zt_v;
#define ZT_N 1
#define ZT_LD(p) vld1q_f64(p)
#define ZT_ST(p, v) vst1q_f64(p, v)
#define ZT_DUP(r, i) ((float64x2_t){r, i})
#define ZT_ZERO() vdupq_n_f64(0)
#define ZT_C0(c, n, m) vcmlaq_f64(c, n, m)
#define ZT_C90(c, n, m) vcmlaq_rot90_f64(c, n, m)
#else
typedef float32x4_t zt_v;
#define ZT_N 2
#define ZT_LD(p) vld1q_f32(p)
#define ZT_ST(p, v) vst1q_f32(p, v)
#define ZT_DUP(r, i) ((float32x4_t){r, i, r, i})
#define ZT_ZERO() vdupq_n_f32(0)
#define ZT_C0(c, n, m) vcmlaq_f32(c, n, m)
#define ZT_C90(c, n, m) vcmlaq_rot90_f32(c, n, m)
#endif

/* Panel width at j0 of n as the generic copies lay it out: full panels of u, then 8, 4, 2, 1. */
static inline int zt_w(BLASLONG n, BLASLONG j0, int u) {
  BLASLONG r = n - j0;
  if (r >= u) return u;
  int w = 1;
  while (2 * w <= r) w *= 2;
  return w;
}

#ifndef ZT_LEFT
/* c[0..len) += v[0..len) s for complex numbers. */
static inline void zt_axpy(int len, FLOAT *c, const FLOAT *v, FLOAT sr, FLOAT si) {
  const zt_v m = ZT_DUP(sr, si);
  int r = 0;
  for (; r + ZT_N <= len; r += ZT_N) {
    const zt_v n = ZT_LD(v + 2 * r);
    ZT_ST(c + 2 * r, ZT_C90(ZT_C0(ZT_LD(c + 2 * r), n, m), n, m));
  }
  for (; r < len; ++r) {
    const FLOAT vr = v[2 * r], vi = v[2 * r + 1];
    c[2 * r] += vr * sr - vi * si;
    c[2 * r + 1] += vr * si + vi * sr;
  }
}

/* x = c = c s for len complex numbers. */
static inline void zt_scale(int len, FLOAT *c, FLOAT *x, FLOAT sr, FLOAT si) {
  const zt_v m = ZT_DUP(sr, si);
  int r = 0;
  for (; r + ZT_N <= len; r += ZT_N) {
    const zt_v n = ZT_LD(c + 2 * r);
    const zt_v v = ZT_C90(ZT_C0(ZT_ZERO(), n, m), n, m);
    ZT_ST(c + 2 * r, v);
    ZT_ST(x + 2 * r, v);
  }
  for (; r < len; ++r) {
    const FLOAT cr = c[2 * r], ci = c[2 * r + 1];
    x[2 * r] = c[2 * r] = cr * sr - ci * si;
    x[2 * r + 1] = c[2 * r + 1] = cr * si + ci * sr;
  }
}
#endif

#ifdef ZT_LEFT
/* Left solve on real and imaginary planes of the block with rows contiguous: plain vectorizable loops along wn. */
static inline void zt_solve_l(int wm, const int wn, const FLOAT *t, FLOAT *x, FLOAT *c, BLASLONG ldc) {
  const FLOAT sg = ZT_CJ ? -1 : 1;
  FLOAT pr[GEMM_DEFAULT_UNROLL_M * GEMM_DEFAULT_UNROLL_N], pi[GEMM_DEFAULT_UNROLL_M * GEMM_DEFAULT_UNROLL_N];
  for (int j = 0; j < wn; ++j)
    for (int i = 0; i < wm; ++i) {
      pr[i * wn + j] = c[2 * (i + j * ldc)];
      pi[i * wn + j] = c[2 * (i + j * ldc) + 1];
    }
  for (int u = 0; u < wm; ++u) {
#ifdef ZT_BACK
    const int i = wm - 1 - u, r0 = 0, r1 = i;
#else
    const int i = u, r0 = i + 1, r1 = wm;
#endif
    const FLOAT *ti = t + 2 * (BLASLONG)i * wm;
    const FLOAT dr = ti[2 * i], di = sg * ti[2 * i + 1];
    FLOAT *restrict xr = pr + i * wn, *restrict xi = pi + i * wn, *restrict xo = x + 2 * (BLASLONG)i * wn;
    for (int j = 0; j < wn; ++j) {
      const FLOAT a = xr[j], b = xi[j];
      xr[j] = a * dr - b * di;
      xi[j] = a * di + b * dr;
      xo[2 * j] = xr[j];
      xo[2 * j + 1] = xi[j];
    }
    for (int r = r0; r < r1; ++r) {
      const FLOAT tr = ti[2 * r], tm = sg * ti[2 * r + 1];
      FLOAT *restrict yr = pr + r * wn, *restrict yi = pi + r * wn;
      for (int j = 0; j < wn; ++j) {
        yr[j] -= tr * xr[j] - tm * xi[j];
        yi[j] -= tr * xi[j] + tm * xr[j];
      }
    }
  }
  for (int j = 0; j < wn; ++j)
    for (int i = 0; i < wm; ++i) {
      c[2 * (i + j * ldc)] = pr[i * wn + j];
      c[2 * (i + j * ldc) + 1] = pi[i * wn + j];
    }
}
#endif

/* Solve of a wm x wn block against the triangle t (inverted diagonal); the solution also goes to x. */
static void zt_solve(int wm, int wn, const FLOAT *t, FLOAT *x, FLOAT *c, BLASLONG ldc) {
#ifdef ZT_LEFT
  if (wn == GEMM_DEFAULT_UNROLL_N) zt_solve_l(wm, GEMM_DEFAULT_UNROLL_N, t, x, c, ldc);
  else zt_solve_l(wm, wn, t, x, c, ldc);
#else
  const FLOAT sg = ZT_CJ ? -1 : 1;
  for (int u = 0; u < wn; ++u) {
#ifdef ZT_BACK
    const int i = wn - 1 - u, r0 = 0, r1 = i;
#else
    const int i = u, r0 = i + 1, r1 = wn;
#endif
    const FLOAT *ti = t + 2 * (BLASLONG)i * wn;
    FLOAT *xi = x + 2 * (BLASLONG)i * wm;
    zt_scale(wm, c + 2 * i * ldc, xi, ti[2 * i], sg * ti[2 * i + 1]);
    for (int r = r0; r < r1; ++r) zt_axpy(wm, c + 2 * r * ldc, xi, -ti[2 * r], -sg * ti[2 * r + 1]);
  }
#endif
}

/* C block -= au bu over ku steps, then the solve. */
static void zt_block(BLASLONG ku, FLOAT *au, FLOAT *bu, const FLOAT *t, FLOAT *x, int wm, int wn, FLOAT *c,
                     BLASLONG ldc) {
  if (ku > 0) ZT_GEMM(wm, wn, ku, -ONE, ZERO, au, bu, c, ldc);
  zt_solve(wm, wn, t, x, c, ldc);
}

int CNAME(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT dummy1, FLOAT dummy2, FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc,
          BLASLONG offset) {
  const int um = GEMM_UNROLL_M, un = GEMM_UNROLL_N;
  if (m <= 0 || n <= 0) return 0;
#ifdef ZT_LEFT
  for (BLASLONG j0 = 0; j0 < n;) {
    const int wn = zt_w(n, j0, un);
    FLOAT *bj = b + 2 * j0 * k, *cj = c + 2 * j0 * ldc;
#ifdef LT
    BLASLONG kk = offset;
    for (BLASLONG i0 = 0; i0 < m;) {
      const int wm = zt_w(m, i0, um);
      FLOAT *aa = a + 2 * i0 * k;
      zt_block(kk, aa, bj, aa + 2 * kk * wm, bj + 2 * kk * wn, wm, wn, cj + 2 * i0, ldc);
      kk += wm;
      i0 += wm;
    }
#else
    /* A panels from the last one back: the tails (1, 2, 4, ... rows), then the full panels */
    BLASLONG kk = m + offset, i1 = m;
    for (int wm = 1; i1 > 0; wm = wm < um ? 2 * wm : um) {
      if (wm < um && !(m & wm)) continue;
      const BLASLONG i0 = i1 - wm;
      FLOAT *aa = a + 2 * i0 * k;
      zt_block(k - kk, aa + 2 * wm * kk, bj + 2 * wn * kk, aa + 2 * (kk - wm) * wm, bj + 2 * (kk - wm) * wn, wm, wn,
               cj + 2 * i0, ldc);
      kk -= wm;
      i1 = i0;
    }
#endif
    j0 += wn;
  }
#else
#ifdef RN
  BLASLONG kk = -offset;
  for (BLASLONG j0 = 0; j0 < n;) {
    const int wn = zt_w(n, j0, un);
    FLOAT *bj = b + 2 * j0 * k;
    for (BLASLONG i0 = 0; i0 < m;) {
      const int wm = zt_w(m, i0, um);
      FLOAT *aa = a + 2 * i0 * k;
      zt_block(kk, aa, bj, bj + 2 * kk * wn, aa + 2 * kk * wm, wm, wn, c + 2 * (i0 + j0 * ldc), ldc);
      i0 += wm;
    }
    kk += wn;
    j0 += wn;
  }
#else
  /* B panels from the last one back: the tails (1, 2, 4, ... columns), then the full panels */
  BLASLONG kk = n - offset, j1 = n;
  for (int wn = 1; j1 > 0; wn = wn < un ? 2 * wn : un) {
    if (wn < un && !(n & wn)) continue;
    const BLASLONG j0 = j1 - wn;
    FLOAT *bj = b + 2 * j0 * k;
    for (BLASLONG i0 = 0; i0 < m;) {
      const int wm = zt_w(m, i0, um);
      FLOAT *aa = a + 2 * i0 * k;
      zt_block(k - kk, aa + 2 * wm * kk, bj + 2 * wn * kk, bj + 2 * (kk - wn) * wn, aa + 2 * (kk - wn) * wm, wm, wn,
               c + 2 * (i0 + j0 * ldc), ldc);
      i0 += wm;
    }
    kk -= wn;
    j1 = j0;
  }
#endif
#endif
  return 0;
}
