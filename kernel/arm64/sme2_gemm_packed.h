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

/* A panels: wm <= 64 rows at a[l * wm + i], the last one the remainder; B panels: VL columns, then 8, 4, 2, 1. */

#ifndef SME2_GEMM_PACKED_H
#define SME2_GEMM_PACKED_H

#define S2P_UM GEMM_DEFAULT_UNROLL_M
#define S2P_UN GEMM_DEFAULT_UNROLL_N

#ifdef HAVE_SME2_GEMM
#define S2P_SC __arm_streaming_compatible
#else
#define S2P_SC
#endif

/* Width of the B panel at column j0 of n: full panels, then the largest power of two left. */
static inline int s2p_wn(BLASLONG n, BLASLONG j0) S2P_SC {
  BLASLONG r = n - j0;
  if (r >= S2P_UN) return S2P_UN;
  int w = 1;
  while (2 * w <= r) w *= 2;
  return w;
}

#ifdef HAVE_SME2_GEMM
#include "sme2_gemm_tile.h"

#if S2P_UM != S2_NR || S2P_UN != S2_VL
#error "the SME2 packed kernels need GEMM_UNROLL_M 64 and GEMM_UNROLL_N of one streaming vector"
#endif

#pragma clang attribute push(__attribute__((target("sme2"))), apply_to = function)

S2_INL void s2p_op(const int SUB, int t, svbool_t p, s2_v a, s2_v b) S2_S S2_ZA {
  if (SUB) s2_mops(t, p, p, a, b);
  else s2_mopa(t, p, a, b);
}

/* ZA row j (column j of C), tiles 0..TN-1 (its rows): += (-= with SUB) A panel times B panel over kb steps. */
S2_INL void s2p_mac(const int TN, const int BFULL, const int SUB, BLASLONG kb, const FLOAT *a, int wm, const FLOAT *b,
                    int wn) S2_S S2_ZA {
  const svbool_t pt = S2_PT();
  BLASLONG l = 0;
  if (BFULL)
    for (; l + 4 <= kb; l += 4) {
      const s2_v4 bb = s2_ld4(b + l * S2_VL);
      S2_UNROLL for (int U = 0; U < 4; ++U) {
        const FLOAT *ar = a + (l + U) * wm;
        const s2_v4 a0 = s2_ld4n(ar, wm);
        const s2_v4 a1 = TN > 4 ? s2_ld4n(ar + 4 * S2_VL, wm - 4 * S2_VL) : a0;
        S2_UNROLL for (int Q = 0; Q < TN; ++Q) s2p_op(SUB, Q, pt, s2_get4(bb, U), s2_pick(a0, a1, Q));
      }
    }
  const svbool_t pb = S2_PW(0, wn);
  for (; l < kb; ++l) {
    const s2_v bv = S2_LD1(pb, b + l * wn);
    const FLOAT *ar = a + l * wm;
    const s2_v4 a0 = s2_ld4n(ar, wm);
    const s2_v4 a1 = TN > 4 ? s2_ld4n(ar + 4 * S2_VL, wm - 4 * S2_VL) : a0;
    S2_UNROLL for (int Q = 0; Q < TN; ++Q) s2p_op(SUB, Q, pt, bv, s2_pick(a0, a1, Q));
  }
}

#pragma clang attribute pop

#endif
#endif
