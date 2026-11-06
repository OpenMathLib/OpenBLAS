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

/* Copies into the A and B panels of the SME2 kernels: m rows of w <= 64 columns, element (i, c) to b[i * w + c]. */

#ifndef SME2_GEMM_COPY_H
#define SME2_GEMM_COPY_H

#include <arm_neon.h>
#include <string.h>
#include "sme2_gemm_detect.h"

#define S2C_W GEMM_DEFAULT_UNROLL_M

static inline BLASLONG s2c_clip(BLASLONG x, BLASLONG lo, BLASLONG hi) { return x < lo ? lo : (x > hi ? hi : x); }

#ifdef HAVE_SME2_GEMM
/* Streaming-mode copies: the SME unit moves data faster than NEON. Below S2C_SME elements NEON is faster. */
#define S2C_SME 512
#ifdef DOUBLE
#define S2C_VL 8
#define S2C_WL(a, b) svwhilelt_b64_s64(a, b)
#define S2C_LD svld1_f64
#define S2C_ST svst1_f64
#define S2C_DUP svdup_n_f64
#else
#define S2C_VL 16
#define S2C_WL(a, b) svwhilelt_b32_s64(a, b)
#define S2C_LD svld1_f32
#define S2C_ST svst1_f32
#define S2C_DUP svdup_n_f32
#endif

static inline __attribute__((always_inline)) void s2c_rows_s(BLASLONG m, BLASLONG w, const FLOAT *a, BLASLONG lda,
                                                              FLOAT *b) __arm_streaming {
  for (BLASLONG c = 0; c < w; c += S2C_VL) {
    const svbool_t p = S2C_WL(c, w);
    BLASLONG i = 0;
    for (; i + 4 <= m; i += 4) {
      const FLOAT *q = a + i * lda + c;
      FLOAT *d = b + i * w + c;
      S2C_ST(p, d, S2C_LD(p, q));
      S2C_ST(p, d + w, S2C_LD(p, q + lda));
      S2C_ST(p, d + 2 * w, S2C_LD(p, q + 2 * lda));
      S2C_ST(p, d + 3 * w, S2C_LD(p, q + 3 * lda));
    }
    for (; i < m; ++i) S2C_ST(p, b + i * w + c, S2C_LD(p, a + i * lda + c));
  }
}

/* Tile numbers must be literals; with a constant t the switches fold away. */
static inline __attribute__((always_inline)) void s2c_ldh(int t, uint32_t s, svbool_t p, const FLOAT *q) __arm_streaming
    __arm_inout("za") {
#ifdef DOUBLE
  switch (t) {
    case 0: svld1_hor_za64(0, s, p, q); break;
    case 1: svld1_hor_za64(1, s, p, q); break;
    case 2: svld1_hor_za64(2, s, p, q); break;
    case 3: svld1_hor_za64(3, s, p, q); break;
    case 4: svld1_hor_za64(4, s, p, q); break;
    case 5: svld1_hor_za64(5, s, p, q); break;
    case 6: svld1_hor_za64(6, s, p, q); break;
    default: svld1_hor_za64(7, s, p, q); break;
  }
#else
  switch (t) {
    case 0: svld1_hor_za32(0, s, p, q); break;
    case 1: svld1_hor_za32(1, s, p, q); break;
    case 2: svld1_hor_za32(2, s, p, q); break;
    default: svld1_hor_za32(3, s, p, q); break;
  }
#endif
}

static inline __attribute__((always_inline)) void s2c_stv(int t, uint32_t s, svbool_t p, FLOAT *q) __arm_streaming
    __arm_inout("za") {
#ifdef DOUBLE
  switch (t) {
    case 0: svst1_ver_za64(0, s, p, q); break;
    case 1: svst1_ver_za64(1, s, p, q); break;
    case 2: svst1_ver_za64(2, s, p, q); break;
    case 3: svst1_ver_za64(3, s, p, q); break;
    case 4: svst1_ver_za64(4, s, p, q); break;
    case 5: svst1_ver_za64(5, s, p, q); break;
    case 6: svst1_ver_za64(6, s, p, q); break;
    default: svst1_ver_za64(7, s, p, q); break;
  }
#else
  switch (t) {
    case 0: svst1_ver_za32(0, s, p, q); break;
    case 1: svst1_ver_za32(1, s, p, q); break;
    case 2: svst1_ver_za32(2, s, p, q); break;
    default: svst1_ver_za32(3, s, p, q); break;
  }
#endif
}

/* Columns c go into horizontal slice c % VL of tile c / VL and leave transposed through vertical slices. */
static inline __attribute__((always_inline)) void s2c_cols_s(BLASLONG m, BLASLONG w, const FLOAT *a, BLASLONG lda,
                                                              FLOAT *b) __arm_streaming __arm_inout("za") {
  BLASLONG i = 0;
  if (w <= S2C_VL) /* B panels: two tiles in turn, so that loads into one overlap stores from the other */
    for (; i + 2 * S2C_VL <= m; i += 2 * S2C_VL) {
      const svbool_t pt = S2C_WL(0, S2C_VL), pw = S2C_WL(0, w);
      for (int c = 0; c < w; ++c) {
        s2c_ldh(0, c, pt, a + i + c * lda);
        s2c_ldh(1, c, pt, a + i + S2C_VL + c * lda);
      }
      for (int r = 0; r < S2C_VL; ++r) s2c_stv(0, r, pw, b + (i + r) * w);
      for (int r = 0; r < S2C_VL; ++r) s2c_stv(1, r, pw, b + (i + S2C_VL + r) * w);
    }
  for (; i < m; i += S2C_VL) {
    const svbool_t pm = S2C_WL(i, m);
    const int rows = m - i < S2C_VL ? (int)(m - i) : S2C_VL;
    _Pragma("GCC unroll 8") for (int t = 0; t < S2C_W / S2C_VL; ++t) {
      if (t * S2C_VL >= w) break;
      for (int s = 0; s < S2C_VL && t * S2C_VL + s < w; ++s) s2c_ldh(t, s, pm, a + i + (t * S2C_VL + s) * lda);
    }
    for (int r = 0; r < rows; ++r) {
      FLOAT *d = b + (i + r) * w;
      _Pragma("GCC unroll 8") for (int t = 0; t < S2C_W / S2C_VL; ++t) {
        if (t * S2C_VL >= w) break;
        s2c_stv(t, r, S2C_WL(t * S2C_VL, w), d + t * S2C_VL);
      }
    }
  }
}

__arm_locally_streaming static __attribute__((noinline)) void s2c_rows_sme(BLASLONG m, BLASLONG w, const FLOAT *a,
                                                                            BLASLONG lda, FLOAT *b) {
  s2c_rows_s(m, w, a, lda, b);
}

__arm_new("za") __arm_locally_streaming static __attribute__((noinline)) void s2c_cols_sme(BLASLONG m, BLASLONG w,
                                                                                           const FLOAT *a, BLASLONG lda,
                                                                                           FLOAT *b) {
  s2c_cols_s(m, w, a, lda, b);
}
#endif

/* Element (i, c) at a[c + i * lda]. */
static inline void s2c_rows(BLASLONG m, BLASLONG w, const FLOAT *a, BLASLONG lda, FLOAT *b) {
#ifdef HAVE_SME2_GEMM
  if (m * w >= S2C_SME && s2_usable()) {
    s2c_rows_sme(m, w, a, lda, b);
    return;
  }
#endif
  for (BLASLONG i = 0; i < m; ++i) memcpy(b + i * w, a + i * lda, w * sizeof(FLOAT));
}

/* Element (i, c) at a[i + c * lda], transposed in registers 4 x 4 (fp32) or 2 x 2 (fp64) at a time. */
static inline void s2c_cols(BLASLONG m, BLASLONG w, const FLOAT *a, BLASLONG lda, FLOAT *b) {
#ifdef HAVE_SME2_GEMM
  if (m * w >= S2C_SME && s2_usable()) {
    s2c_cols_sme(m, w, a, lda, b);
    return;
  }
#endif
  BLASLONG c = 0;
#ifdef DOUBLE
  for (; c + 2 <= w; c += 2) {
    const FLOAT *a0 = a + c * lda, *a1 = a0 + lda;
    BLASLONG i = 0;
    for (; i + 2 <= m; i += 2) {
      const float64x2_t x0 = vld1q_f64(a0 + i), x1 = vld1q_f64(a1 + i);
      vst1q_f64(b + i * w + c, vtrn1q_f64(x0, x1));
      vst1q_f64(b + (i + 1) * w + c, vtrn2q_f64(x0, x1));
    }
    for (; i < m; ++i) {
      b[i * w + c] = a0[i];
      b[i * w + c + 1] = a1[i];
    }
  }
#else
  for (; c + 4 <= w; c += 4) {
    const FLOAT *a0 = a + c * lda, *a1 = a0 + lda, *a2 = a1 + lda, *a3 = a2 + lda;
    BLASLONG i = 0;
    for (; i + 4 <= m; i += 4) {
      const float32x4_t x0 = vld1q_f32(a0 + i), x1 = vld1q_f32(a1 + i), x2 = vld1q_f32(a2 + i),
                        x3 = vld1q_f32(a3 + i);
      const float64x2_t t0 = vreinterpretq_f64_f32(vtrn1q_f32(x0, x1)), t1 = vreinterpretq_f64_f32(vtrn2q_f32(x0, x1)),
                        t2 = vreinterpretq_f64_f32(vtrn1q_f32(x2, x3)), t3 = vreinterpretq_f64_f32(vtrn2q_f32(x2, x3));
      vst1q_f32(b + i * w + c, vreinterpretq_f32_f64(vtrn1q_f64(t0, t2)));
      vst1q_f32(b + (i + 1) * w + c, vreinterpretq_f32_f64(vtrn1q_f64(t1, t3)));
      vst1q_f32(b + (i + 2) * w + c, vreinterpretq_f32_f64(vtrn2q_f64(t0, t2)));
      vst1q_f32(b + (i + 3) * w + c, vreinterpretq_f32_f64(vtrn2q_f64(t1, t3)));
    }
    for (; i < m; ++i) {
      b[i * w + c] = a0[i];
      b[i * w + c + 1] = a1[i];
      b[i * w + c + 2] = a2[i];
      b[i * w + c + 3] = a3[i];
    }
  }
#endif
  for (; c < w; ++c)
    for (BLASLONG i = 0; i < m; ++i) b[i * w + c] = a[i + c * lda];
}

/* B panels of GEMM_DEFAULT_UNROLL_N columns, then 8/4/2/1 (the layout of generic/gemm_[nt]copy_16.c). */
#ifdef HAVE_SME2_GEMM
__arm_new("za") __arm_locally_streaming static __attribute__((noinline)) void s2c_b_panels_sme(int rows, BLASLONG m,
                                                                                               BLASLONG n,
                                                                                               const FLOAT *a,
                                                                                               BLASLONG lda, FLOAT *b) {
  for (BLASLONG j = 0, w = GEMM_DEFAULT_UNROLL_N; j < n; j += w, b += m * w) {
    while (w > n - j) w >>= 1;
    if (rows) s2c_rows_s(m, w, a + j, lda, b);
    else s2c_cols_s(m, w, a + j * lda, lda, b);
  }
}
#endif

static inline void s2c_b_panels(int rows, BLASLONG m, BLASLONG n, const FLOAT *a, BLASLONG lda, FLOAT *b) {
#ifdef HAVE_SME2_GEMM
  if (m * n >= S2C_SME && s2_usable()) {
    s2c_b_panels_sme(rows, m, n, a, lda, b);
    return;
  }
#endif
  for (BLASLONG j = 0, w = GEMM_DEFAULT_UNROLL_N; j < n; j += w, b += m * w) {
    while (w > n - j) w >>= 1;
    if (rows) s2c_rows(m, w, a + j, lda, b);
    else s2c_cols(m, w, a + j * lda, lda, b);
  }
}

static inline void s2c_zero(BLASLONG m, BLASLONG w, FLOAT *b) {
  if (m > 0) memset(b, 0, (size_t)(m * w) * sizeof(FLOAT));
}

/* Rows i1..i2 of a diagonal block, diagonal in column d0 + i: zeros left (zl) or right of it, 1 there with unit. */
static inline void s2c_tri(BLASLONG i1, BLASLONG i2, BLASLONG w, BLASLONG d0, int zl, int unit, FLOAT *b) {
  for (BLASLONG i = i1; i < i2; ++i) {
    const BLASLONG d = d0 + i;
    FLOAT *bi = b + i * w;
    for (BLASLONG c = zl ? 0 : d + 1; c < (zl ? d : w); ++c) bi[c] = ZERO;
    if (unit) bi[d] = ONE;
  }
}

#ifdef HAVE_SME2_GEMM
static inline __attribute__((always_inline)) void s2c_zero_s(BLASLONG m, BLASLONG w, FLOAT *b) __arm_streaming {
  for (BLASLONG k = 0; k < m * w; k += S2C_VL) S2C_ST(S2C_WL(k, m * w), b + k, S2C_DUP(0));
}

static inline __attribute__((always_inline)) void s2c_tri_s(BLASLONG i1, BLASLONG i2, BLASLONG w, BLASLONG d0, int zl,
                                                            int unit, FLOAT *b) __arm_streaming {
  for (BLASLONG i = i1; i < i2; ++i) {
    const BLASLONG d = d0 + i;
    for (BLASLONG c = 0; c < w; c += S2C_VL) {
      const svbool_t pw = S2C_WL(c, w), pd = S2C_WL(c, d), pd1 = S2C_WL(c, d + 1);
      S2C_ST(zl ? pd : svbic_b_z(pw, pw, pd1), b + i * w + c, S2C_DUP(0));
      if (unit) S2C_ST(svbic_b_z(pw, pd1, pd), b + i * w + c, S2C_DUP(1));
    }
  }
}
#endif

#endif
