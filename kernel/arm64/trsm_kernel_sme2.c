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

/* SME2 TRSM kernels: each C block takes the update, then the solve in ZA; plain loops without SME2. */

#include "common.h"
#include "sme2_gemm_detect.h"
#include "sme2_gemm_packed.h"

#if defined(LN) || defined(LT)
#define S2T_LEFT 1
#endif
#if defined(LN) || defined(RT)
#define S2T_BACK 1
#endif

/* C block -= au bu over ku steps, then solved against the triangle t; the solution also goes to x. */
static void s2t_block_ref(BLASLONG ku, const FLOAT *au, const FLOAT *bu, const FLOAT *t, FLOAT *x, int wm, int wn,
                          FLOAT *c, BLASLONG ldc) {
  for (int j = 0; j < wn; ++j)
    for (int i = 0; i < wm; ++i) {
      FLOAT s = 0;
      for (BLASLONG l = 0; l < ku; ++l) s += au[l * wm + i] * bu[l * wn + j];
      c[i + j * ldc] -= s;
    }
#ifdef S2T_LEFT
  for (int u = 0; u < wm; ++u) {
#ifdef S2T_BACK
    const int s = wm - 1 - u, r0 = 0, r1 = s;
#else
    const int s = u, r0 = s + 1, r1 = wm;
#endif
    const FLOAT *tr = t + (BLASLONG)s * wm;
    for (int j = 0; j < wn; ++j) {
      const FLOAT v = c[s + j * ldc] * tr[s];
      x[s * wn + j] = v;
      c[s + j * ldc] = v;
      for (int r = r0; r < r1; ++r) c[r + j * ldc] -= v * tr[r];
    }
  }
#else
  for (int u = 0; u < wn; ++u) {
#ifdef S2T_BACK
    const int s = wn - 1 - u, r0 = 0, r1 = s;
#else
    const int s = u, r0 = s + 1, r1 = wn;
#endif
    const FLOAT *tr = t + (BLASLONG)s * wn;
    for (int i = 0; i < wm; ++i) {
      const FLOAT v = c[i + s * ldc] * tr[s];
      x[s * wm + i] = v;
      c[i + s * ldc] = v;
      for (int r = r0; r < r1; ++r) c[i + r * ldc] -= v * tr[r];
    }
  }
#endif
}

#ifdef HAVE_SME2_GEMM
#pragma clang attribute push(__attribute__((target("sme2"))), apply_to = function)

/* The solve in ZA; tile rows are columns of C. */
S2_INL void s2t_solve(const int TN, const FLOAT *t, FLOAT *x, int wm, int wn) S2_S S2_ZA {
  const svbool_t pt = S2_PT();
#ifdef S2T_LEFT
  const svbool_t pj = S2_PW(0, wn);
  for (int u = 0; u < wm; ++u) {
#ifdef S2T_BACK
    const int s = wm - 1 - u;
#else
    const int s = u;
#endif
    const int q = s / S2_VL;
    const FLOAT *tr = t + (BLASLONG)s * wm;
    const s2_v v = S2_MUL(s2_rdv(q, s % S2_VL), tr[s]);
    s2_wrv(q, s % S2_VL, v);
    S2_ST1(pj, x + (BLASLONG)s * wn, v);
    S2_UNROLL for (int Q = 0; Q < TN; ++Q) {
#ifdef S2T_BACK
      if (Q > q) break;
      const svbool_t pm = S2_PW(Q * S2_VL, s);
#else
      if (Q < q) continue;
      const svbool_t pm = svbic_b_z(pt, S2_PW(Q * S2_VL, wm), S2_PW(Q * S2_VL, s + 1));
#endif
      s2_mops(Q, pt, pm, v, S2_LD1(pm, tr + Q * S2_VL));
    }
  }
#else
  for (int u = 0; u < wn; ++u) {
#ifdef S2T_BACK
    const int s = wn - 1 - u;
    const svbool_t pr = S2_PW(0, s);
#else
    const int s = u;
    const svbool_t pr = svbic_b_z(pt, S2_PW(0, wn), S2_PW(0, s + 1));
#endif
    const FLOAT *tr = t + (BLASLONG)s * wn;
    const s2_v tv = S2_LD1(pr, tr);
    S2_UNROLL for (int Q = 0; Q < TN; ++Q) {
      const s2_v v = S2_MUL(s2_rdh(Q, s), tr[s]);
      s2_wrh(Q, s, v);
      S2_ST1(S2_PW(Q * S2_VL, wm), x + (BLASLONG)s * wm + Q * S2_VL, v);
      s2_mops(Q, pr, pt, tv, v);
    }
  }
#endif
}

typedef void (*s2t_bfn)(BLASLONG, const FLOAT *, const FLOAT *, const FLOAT *, FLOAT *, int, int, FLOAT *, BLASLONG);

#define S2T_DEF(TN)                                                                                              \
  __arm_new("za") __arm_locally_streaming static void s2t_blk_##TN(BLASLONG ku, const FLOAT *au, const FLOAT *bu,  \
                                                                   const FLOAT *t, FLOAT *x, int wm, int wn,     \
                                                                   FLOAT *c, BLASLONG ldc) {                     \
    s2_c_load(1, TN, c, ldc, wn, wm, S2_LOAD, (FLOAT)1);                                                         \
    if (wn == S2_VL) s2p_mac(TN, 1, 1, ku, au, wm, bu, wn);                                                      \
    else s2p_mac(TN, 0, 1, ku, au, wm, bu, wn);                                                                  \
    s2t_solve(TN, t, x, wm, wn);                                                                                 \
    s2_c_store(1, TN, c, ldc, wn, wm);                                                                           \
  }
S2T_DEF(1)
S2T_DEF(2)
S2T_DEF(3)
S2T_DEF(4)
#ifdef DOUBLE
S2T_DEF(5)
S2T_DEF(6)
S2T_DEF(7)
S2T_DEF(8)
static const s2t_bfn s2t_blk[S2_NT] = {s2t_blk_1, s2t_blk_2, s2t_blk_3, s2t_blk_4,
                                       s2t_blk_5, s2t_blk_6, s2t_blk_7, s2t_blk_8};
#else
static const s2t_bfn s2t_blk[S2_NT] = {s2t_blk_1, s2t_blk_2, s2t_blk_3, s2t_blk_4};
#endif

#pragma clang attribute pop
#endif

static void s2t_block(int sme, BLASLONG ku, const FLOAT *au, const FLOAT *bu, const FLOAT *t, FLOAT *x, int wm, int wn,
                      FLOAT *c, BLASLONG ldc) {
  if (ku < 0) ku = 0;
#ifdef HAVE_SME2_GEMM
  if (sme) {
    s2t_blk[(wm + S2_VL - 1) / S2_VL - 1](ku, au, bu, t, x, wm, wn, c, ldc);
    return;
  }
#endif
  s2t_block_ref(ku, au, bu, t, x, wm, wn, c, ldc);
}

int CNAME(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT dummy1, FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc,
          BLASLONG offset) {
  if (m <= 0 || n <= 0) return 0;
  int sme = 0;
#ifdef HAVE_SME2_GEMM
  sme = s2_usable();
#endif
#if defined(LT) || defined(LN)
  for (BLASLONG j0 = 0; j0 < n;) {
    const int wn = s2p_wn(n, j0);
    FLOAT *bj = b + j0 * k, *cj = c + j0 * ldc;
#ifdef LT
    BLASLONG kk = offset;
    for (BLASLONG i0 = 0; i0 < m; i0 += S2P_UM) {
      const int wm = m - i0 < S2P_UM ? (int)(m - i0) : S2P_UM;
      const FLOAT *aa = a + i0 * k;
      s2t_block(sme, kk, aa, bj, aa + kk * wm, bj + kk * wn, wm, wn, cj + i0, ldc);
      kk += wm;
    }
#else
    /* A panels from the last one back: the one of m % 64 rows, then the full ones */
    BLASLONG kk = m + offset;
    for (BLASLONG i0 = m % S2P_UM ? m - m % S2P_UM : m - S2P_UM; i0 >= 0; i0 -= S2P_UM) {
      const int wm = m - i0 < S2P_UM ? (int)(m - i0) : S2P_UM;
      const FLOAT *aa = a + i0 * k;
      s2t_block(sme, k - kk, aa + wm * kk, bj + wn * kk, aa + (kk - wm) * wm, bj + (kk - wm) * wn, wm, wn,
                cj + i0, ldc);
      kk -= wm;
    }
#endif
    j0 += wn;
  }
#else
#ifdef RN
  BLASLONG kk = -offset;
  for (BLASLONG j0 = 0; j0 < n;) {
    const int wn = s2p_wn(n, j0);
    FLOAT *bj = b + j0 * k;
    for (BLASLONG i0 = 0; i0 < m; i0 += S2P_UM) {
      const int wm = m - i0 < S2P_UM ? (int)(m - i0) : S2P_UM;
      FLOAT *aa = a + i0 * k;
      s2t_block(sme, kk, aa, bj, bj + kk * wn, aa + kk * wm, wm, wn, c + i0 + j0 * ldc, ldc);
    }
    kk += wn;
    j0 += wn;
  }
#else
  /* B panels from the last one back: the tails (1, 2, 4, ... columns), then the full panels */
  BLASLONG kk = n - offset, j1 = n;
  for (int wn = 1; j1 > 0; wn = wn < S2P_UN ? 2 * wn : S2P_UN) {
    if (wn < S2P_UN && !(n & wn)) continue;
    const BLASLONG j0 = j1 - wn;
    FLOAT *bj = b + j0 * k;
    for (BLASLONG i0 = 0; i0 < m; i0 += S2P_UM) {
      const int wm = m - i0 < S2P_UM ? (int)(m - i0) : S2P_UM;
      FLOAT *aa = a + i0 * k;
      s2t_block(sme, k - kk, aa + wm * kk, bj + wn * kk, bj + (kk - wn) * wn, aa + (kk - wn) * wm, wm, wn,
                c + i0 + j0 * ldc, ldc);
    }
    kk -= wn;
    j1 = j0;
  }
#endif
#endif
  return 0;
}
