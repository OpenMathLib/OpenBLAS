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

/* SME2 GEMM kernel of the level-3 driver, the TRMM kernel with TRMMKERNEL; plain loops without SME2. */

#include "common.h"
#include "sme2_gemm_detect.h"
#include "sme2_gemm_packed.h"

/* Depth range of the block at (i0, j0): the TRMM kernel skips the zero part of the triangle as trmmkernel_4x4.c. */
static inline void s2p_krange(BLASLONG i0, BLASLONG j0, int wm, int wn, BLASLONG k, BLASLONG offset, BLASLONG *ks,
                              BLASLONG *kl) S2P_SC {
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
#pragma clang attribute push(__attribute__((target("sme2"))), apply_to = function)

typedef void (*s2p_bfn)(BLASLONG, const FLOAT *, int, const FLOAT *, int, FLOAT *, BLASLONG, FLOAT, int) S2_S S2_ZA;

/* Block of wm rows and wn columns of C; cm 0: C += A B, 1: C -= A B, 2: C += alpha A B, 3: C = alpha A B. */
#define S2P_DEF(TN, BF, SUB)                                                                                     \
  S2_NOINL void s2p_blk_##TN##_##BF##_##SUB(BLASLONG kb, const FLOAT *a, int wm, const FLOAT *b, int wn, FLOAT *C, \
                                            BLASLONG ldc, FLOAT alpha, int cm) S2_S S2_ZA {                      \
    if (cm < 2) s2_c_load(1, TN, C, ldc, wn, wm, S2_LOAD, (FLOAT)1);                                             \
    else svzero_za();                                                                                            \
    s2p_mac(TN, BF, SUB, kb, a, wm, b, wn);                                                                      \
    if (cm < 2) s2_c_store(1, TN, C, ldc, wn, wm);                                                               \
    else s2_c_store_alpha(TN, C, ldc, wn, wm, alpha, cm == 2);                                                   \
  }
#define S2P_DEF4(TN) S2P_DEF(TN, 0, 0) S2P_DEF(TN, 0, 1) S2P_DEF(TN, 1, 0) S2P_DEF(TN, 1, 1)
#define S2P_ROW(TN) {{s2p_blk_##TN##_0_0, s2p_blk_##TN##_0_1}, {s2p_blk_##TN##_1_0, s2p_blk_##TN##_1_1}}
S2P_DEF4(1)
S2P_DEF4(2)
S2P_DEF4(3)
S2P_DEF4(4)
#ifdef DOUBLE
S2P_DEF4(5)
S2P_DEF4(6)
S2P_DEF4(7)
S2P_DEF4(8)
static const s2p_bfn s2p_blk[S2_NT][2][2] = {S2P_ROW(1), S2P_ROW(2), S2P_ROW(3), S2P_ROW(4),
                                             S2P_ROW(5), S2P_ROW(6), S2P_ROW(7), S2P_ROW(8)};
#else
static const s2p_bfn s2p_blk[S2_NT][2][2] = {S2P_ROW(1), S2P_ROW(2), S2P_ROW(3), S2P_ROW(4)};
#endif

__arm_new("za") __arm_locally_streaming static void s2p_drive(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha,
                                                              const FLOAT *sa, const FLOAT *sb, FLOAT *C,
                                                              BLASLONG ldc, BLASLONG offset) {
#ifdef TRMMKERNEL
  const int cm = 3;
#else
  const int cm = alpha == (FLOAT)1 ? 0 : (alpha == (FLOAT)-1 ? 1 : 2);
#endif
  for (BLASLONG j0 = 0; j0 < n;) {
    const int wn = s2p_wn(n, j0);
    for (BLASLONG i0 = 0; i0 < m; i0 += S2_NR) {
      const int wm = m - i0 < S2_NR ? (int)(m - i0) : S2_NR;
      BLASLONG ks, kl;
      s2p_krange(i0, j0, wm, wn, k, offset, &ks, &kl);
      s2p_blk[(wm + S2_VL - 1) / S2_VL - 1][wn == S2_VL][cm == 1](kl, sa + i0 * k + ks * wm, wm,
                                                                  sb + j0 * k + ks * wn, wn, C + i0 + j0 * ldc,
                                                                  ldc, alpha, cm);
    }
    j0 += wn;
  }
}

#pragma clang attribute pop
#endif

static void s2p_loops(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, const FLOAT *sa, const FLOAT *sb, FLOAT *C,
                      BLASLONG ldc, BLASLONG offset) {
  for (BLASLONG j0 = 0; j0 < n;) {
    const int wn = s2p_wn(n, j0);
    for (BLASLONG i0 = 0; i0 < m; i0 += S2P_UM) {
      const int wm = m - i0 < S2P_UM ? (int)(m - i0) : S2P_UM;
      BLASLONG ks, kl;
      s2p_krange(i0, j0, wm, wn, k, offset, &ks, &kl);
      const FLOAT *a = sa + i0 * k + ks * wm, *b = sb + j0 * k + ks * wn;
      for (int j = 0; j < wn; ++j)
        for (int i = 0; i < wm; ++i) {
          FLOAT s = 0;
          for (BLASLONG l = 0; l < kl; ++l) s += a[l * wm + i] * b[l * wn + j];
#ifdef TRMMKERNEL
          C[i0 + i + (j0 + j) * ldc] = alpha * s;
#else
          C[i0 + i + (j0 + j) * ldc] += alpha * s;
#endif
        }
    }
    j0 += wn;
  }
}

int CNAME(BLASLONG bm, BLASLONG bn, BLASLONG bk, FLOAT alpha, FLOAT *ba, FLOAT *bb, FLOAT *C, BLASLONG ldc
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
    s2p_drive(bm, bn, bk, alpha, ba, bb, C, ldc, offset);
    return 0;
  }
#endif
  s2p_loops(bm, bn, bk, alpha, ba, bb, C, ldc, offset);
  return 0;
}
