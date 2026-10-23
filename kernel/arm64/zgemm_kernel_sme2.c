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

/* SME2 complex GEMM kernel of the level-3 driver, the TRMM kernel with TRMMKERNEL; plain loops without SME2. */

#include "common.h"
#include "sme2_gemm_detect.h"

/* Conjugation of A (bit 1) and of B (bit 0), as the generic complex kernels. */
#if defined(NN) || defined(NT) || defined(TN) || defined(TT)
#define S2Z_CJ 0
#elif defined(NR) || defined(NC) || defined(TR) || defined(TC)
#define S2Z_CJ 1
#elif defined(RN) || defined(RT) || defined(CN) || defined(CT)
#define S2Z_CJ 2
#else
#define S2Z_CJ 3
#endif

#ifdef HAVE_SME2_GEMM
#define S2Z_SC __arm_streaming_compatible
#else
#define S2Z_SC
#endif

/* Width of the panel at j0 of n as the generic copies lay it out: full panels of u, then 8, 4, 2, 1. */
static inline int s2z_w(BLASLONG n, BLASLONG j0, int u) S2Z_SC {
  BLASLONG r = n - j0;
  if (r >= u) return u;
  int w = 1;
  while (2 * w <= r) w *= 2;
  return w;
}

/* Depth range of the block at (i0, j0): the TRMM kernel skips the zero part of the triangle as trmmkernel_4x4.c. */
static inline void s2z_krange(BLASLONG i0, BLASLONG j0, int wm, int wn, BLASLONG k, BLASLONG offset, BLASLONG *ks,
                              BLASLONG *kl) S2Z_SC {
#ifdef TRMMKERNEL
#ifdef LEFT
  const BLASLONG off = offset + i0;
#else
  const BLASLONG off = j0 - offset;
#endif
#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
  BLASLONG s = off, e = k;
#elif defined(LEFT)
  BLASLONG s = 0, e = off + wm;
#else
  BLASLONG s = 0, e = off + wn;
#endif
  if (s < 0) s = 0;
  if (s > k) s = k;
  if (e > k) e = k;
  if (e < s) e = s;
  *ks = s;
  *kl = e - s;
#else
  *ks = 0;
  *kl = k;
#endif
}

#ifdef HAVE_SME2_GEMM
#include "sme2_gemm_tile.h"

/* A panels of 2 (fp32) or 4 (fp64) vectors, B panels of 2 vectors; complex numbers stay interleaved. */
#define S2Z_TA (2 * GEMM_DEFAULT_UNROLL_M / S2_VL)
#define S2Z_TB (2 * GEMM_DEFAULT_UNROLL_N / S2_VL)
#if S2Z_TB != 2 || S2Z_TA * S2Z_TB != S2_NT
#error "the SME2 complex kernel needs GEMM_UNROLL_N of one streaming vector and UNROLL_M x UNROLL_N to fill ZA"
#endif

#pragma clang attribute push(__attribute__((target("sme2"))), apply_to = function)

#ifdef DOUBLE
#define S2Z_CADD(x, y, r) svcadd_f64_x(svptrue_b64(), x, y, r)
#define S2Z_CMLA(c, x, y, r) svcmla_f64_x(svptrue_b64(), c, x, y, r)
#define S2Z_DUPQ(r, i) svdupq_n_f64(r, i)
#else
#define S2Z_CADD(x, y, r) svcadd_f32_x(svptrue_b32(), x, y, r)
#define S2Z_CMLA(c, x, y, r) svcmla_f32_x(svptrue_b32(), c, x, y, r)
#define S2Z_DUPQ(r, i) svdupq_n_f32(r, i, r, i)
#endif

/* Tile ib * TA + ia += B vector ib (rows) x A vector ia; ZA rows 2q, 2q + 1: A re b_q, A im b_q (interleaved). */
S2_INL void s2z_mac(const int TA, const int FULL, BLASLONG kb, const FLOAT *a, int ma, const FLOAT *b,
                    int nb) S2_S S2_ZA {
  const svbool_t pt = S2_PT();
  BLASLONG l = 0;
  if (FULL) {
    for (; l + 2 <= kb; l += 2) {
      const s2_v4 bb = s2_ld4(b + l * 2 * S2_VL);
      s2_v4 a0, a1;
      if (TA == 4) {
        a0 = s2_ld4(a + l * 4 * S2_VL);
        a1 = s2_ld4(a + (l + 1) * 4 * S2_VL);
      } else {
        a0 = s2_ld4(a + l * 2 * S2_VL);
        a1 = a0;
      }
      S2_UNROLL for (int U = 0; U < 2; ++U) {
        S2_UNROLL for (int ib = 0; ib < 2; ++ib) {
          const s2_v bv = s2_get4(bb, 2 * U + ib);
          S2_UNROLL for (int ia = 0; ia < TA; ++ia) {
            const s2_v av = TA == 4 ? s2_get4(U ? a1 : a0, ia) : s2_get4(a0, 2 * U + ia);
            s2_mopa(ib * TA + ia, pt, bv, av);
          }
        }
      }
    }
  }
  const svbool_t pb0 = S2_PW(0, nb), pb1 = S2_PW(S2_VL, nb);
  const int two = nb > S2_VL;
  for (; l < kb; ++l) {
    const FLOAT *ar = a + l * ma, *br = b + l * nb;
    const s2_v b0 = S2_LD1(pb0, br), b1 = two ? S2_LD1(pb1, br + S2_VL) : b0;
    S2_UNROLL for (int ia = 0; ia < TA; ++ia) {
      const s2_v av = S2_LD1(S2_PW(ia * S2_VL, ma), ar + ia * S2_VL);
      s2_mopa(ia, pt, b0, av);
      if (two) s2_mopa(TA + ia, pt, b1, av);
    }
  }
}

/* C (+)= alpha T for the two ZA rows x, y of one column of op(A) op(B) (T = x + i y with the conjugations). */
S2_INL void s2z_put(FLOAT *c, svbool_t pg, s2_v x, s2_v y, s2_v al, const int add) S2_S {
#if S2Z_CJ == 0 || S2Z_CJ == 3
  const s2_v t = S2Z_CADD(x, y, 90);
#else
  const s2_v t = S2Z_CADD(x, y, 270);
#endif
  s2_v v = add ? S2_LD1(pg, c) : S2_ZEROV();
  v = S2Z_CMLA(v, t, al, 0);
#if S2Z_CJ < 2
  v = S2Z_CMLA(v, t, al, 90);
#else
  v = S2Z_CMLA(v, t, al, 270);
#endif
  S2_ST1(pg, c, v);
}

/* C (+)= alpha op(A) op(B) for wm x wn complex numbers from the ZA layout of s2z_mac; FULL: whole block. */
S2_INL void s2z_store(const int TA, const int FULL, FLOAT *C, BLASLONG ldc, int wm, int wn, const s2_v al,
                      const int add) S2_S S2_ZA {
  const svbool_t pt = S2_PT();
  if (FULL) { /* four ZA rows (two columns) per move, all tiles at once for independent chains */
    S2_UNROLL for (int ib = 0; ib < 2; ++ib)
      S2_UNROLL for (int r = 0; r < S2_VL; r += 4) {
        FLOAT *c = C + (BLASLONG)(ib * S2_VL + r) * ldc;
        S2_UNROLL for (int ia = 0; ia < TA; ++ia) {
          const s2_v4 z = s2_rdh4(ib * TA + ia, r);
          s2z_put(c + ia * S2_VL, pt, s2_get4(z, 0), s2_get4(z, 1), al, add);
          s2z_put(c + 2 * ldc + ia * S2_VL, pt, s2_get4(z, 2), s2_get4(z, 3), al, add);
        }
      }
    return;
  }
  const int ma = 2 * wm;
  for (int q = 0; q < wn; ++q) {
    const int ib = 2 * q / S2_VL, r = 2 * q % S2_VL;
    S2_UNROLL for (int ia = 0; ia < TA; ++ia)
      s2z_put(C + 2 * q * ldc + ia * S2_VL, S2_PW(ia * S2_VL, ma), s2_rdh(ib * TA + ia, r), s2_rdh(ib * TA + ia, r + 1),
              al, add);
  }
}

typedef void (*s2z_bfn)(BLASLONG, const FLOAT *, int, const FLOAT *, int, FLOAT *, BLASLONG, FLOAT, FLOAT) S2_S S2_ZA;

#ifdef TRMMKERNEL
#define S2Z_ADD 0
#else
#define S2Z_ADD 1
#endif

/* Block of wm x wn complex numbers of C, TA A vectors, FULL: both panels at full width. */
#define S2Z_DEF(TA, FULL)                                                                                        \
  S2_NOINL void s2z_blk_##TA##_##FULL(BLASLONG kb, const FLOAT *a, int wm, const FLOAT *b, int wn, FLOAT *C,      \
                                      BLASLONG ldc, FLOAT ar, FLOAT ai) S2_S S2_ZA {                             \
    svzero_za();                                                                                                 \
    s2z_mac(TA, FULL, kb, a, 2 * wm, b, 2 * wn);                                                                 \
    s2z_store(TA, FULL, C, ldc, wm, wn, S2Z_DUPQ(ar, ai), S2Z_ADD);                                               \
  }
S2Z_DEF(1, 0)
S2Z_DEF(2, 0)
#ifdef DOUBLE
S2Z_DEF(4, 0)
S2Z_DEF(4, 1)
static const s2z_bfn s2z_blk[4] = {s2z_blk_1_0, s2z_blk_2_0, s2z_blk_4_0, s2z_blk_4_0};
#define S2Z_FULL s2z_blk_4_1
#else
S2Z_DEF(2, 1)
static const s2z_bfn s2z_blk[2] = {s2z_blk_1_0, s2z_blk_2_0};
#define S2Z_FULL s2z_blk_2_1
#endif

__arm_new("za") __arm_locally_streaming static void s2z_drive(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT ar,
                                                              FLOAT ai, const FLOAT *sa, const FLOAT *sb, FLOAT *C,
                                                              BLASLONG ldc, BLASLONG offset) {
  for (BLASLONG j0 = 0; j0 < n;) {
    const int wn = s2z_w(n, j0, GEMM_DEFAULT_UNROLL_N);
    for (BLASLONG i0 = 0; i0 < m;) {
      const int wm = s2z_w(m, i0, GEMM_DEFAULT_UNROLL_M);
      BLASLONG ks, kl;
      s2z_krange(i0, j0, wm, wn, k, offset, &ks, &kl);
      const FLOAT *a = sa + 2 * (i0 * k + ks * wm), *b = sb + 2 * (j0 * k + ks * wn);
      FLOAT *c = C + 2 * (i0 + j0 * ldc);
      if (wm == GEMM_DEFAULT_UNROLL_M && wn == GEMM_DEFAULT_UNROLL_N)
        S2Z_FULL(kl, a, wm, b, wn, c, ldc, ar, ai);
      else
        s2z_blk[(2 * wm + S2_VL - 1) / S2_VL - 1](kl, a, wm, b, wn, c, ldc, ar, ai);
      i0 += wm;
    }
    j0 += wn;
  }
}

#pragma clang attribute pop
#endif

static void s2z_loops(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alr, FLOAT ali, const FLOAT *sa, const FLOAT *sb,
                      FLOAT *C, BLASLONG ldc, BLASLONG offset) {
  const FLOAT sa_i = (S2Z_CJ & 2) ? -1 : 1, sb_i = (S2Z_CJ & 1) ? -1 : 1;
  for (BLASLONG j0 = 0; j0 < n;) {
    const int wn = s2z_w(n, j0, GEMM_DEFAULT_UNROLL_N);
    for (BLASLONG i0 = 0; i0 < m;) {
      const int wm = s2z_w(m, i0, GEMM_DEFAULT_UNROLL_M);
      BLASLONG ks, kl;
      s2z_krange(i0, j0, wm, wn, k, offset, &ks, &kl);
      const FLOAT *a = sa + 2 * (i0 * k + ks * wm), *b = sb + 2 * (j0 * k + ks * wn);
      for (int j = 0; j < wn; ++j)
        for (int i = 0; i < wm; ++i) {
          FLOAT sr = 0, si = 0;
          for (BLASLONG l = 0; l < kl; ++l) {
            const FLOAT xr = a[2 * (l * wm + i)], xi = sa_i * a[2 * (l * wm + i) + 1];
            const FLOAT yr = b[2 * (l * wn + j)], yi = sb_i * b[2 * (l * wn + j) + 1];
            sr += xr * yr - xi * yi;
            si += xr * yi + xi * yr;
          }
          FLOAT *c = C + 2 * (i0 + i + (j0 + j) * ldc);
#ifdef TRMMKERNEL
          c[0] = alr * sr - ali * si;
          c[1] = alr * si + ali * sr;
#else
          c[0] += alr * sr - ali * si;
          c[1] += alr * si + ali * sr;
#endif
        }
      i0 += wm;
    }
    j0 += wn;
  }
}

int CNAME(BLASLONG bm, BLASLONG bn, BLASLONG bk, FLOAT alphar, FLOAT alphai, FLOAT *ba, FLOAT *bb, FLOAT *C,
          BLASLONG ldc
#ifdef TRMMKERNEL
          , BLASLONG offset
#endif
) {
#ifndef TRMMKERNEL
  const BLASLONG offset = 0;
  if (bk <= 0) return 0;
#endif
  if (bm <= 0 || bn <= 0) return 0;
#ifdef HAVE_SME2_GEMM
  if (s2_usable()) {
    s2z_drive(bm, bn, bk, alphar, alphai, ba, bb, C, ldc, offset);
    return 0;
  }
#endif
  s2z_loops(bm, bn, bk, alphar, alphai, ba, bb, C, ldc, offset);
  return 0;
}
