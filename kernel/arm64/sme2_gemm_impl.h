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

/*
 * SME2 real GEMM for a 512-bit streaming vector length, included by sme_sgemm_kernel.c and sme_dgemm_kernel.c.
 *
 * The design follows C. Deng, W. Yang, J. Fang, D. Dong, "Demystifying ARM SME to Optimize General Matrix
 * Multiplications", arXiv:2512.21473, with the changes measured in MTGEMM-A (https://github.com/tesch1/mtgemm-a):
 *  - Goto loop nest on a row-major core problem; column-major C = op(A) op(B) is solved as C^T = op(B)^T op(A)^T.
 *  - mc, nc, kc from the paper's analytical model (TLB bound on kc, L2 bound, maximal compute-to-memory ratio).
 *  - A is packed into panels of VL rows with the transposition done in ZA (rows in, columns out).
 *  - B is packed "online": the first micro-kernel call of a block reads B from the source and writes the packed
 *    panel as a side effect, so there is no separate packing pass.
 *  - Micro-kernel of 1 x NT tiles (fp32: 16 x 64, fp64: 8 x 64) with SME2 multi-vector loads, a 2 x (NT/2) kernel
 *    for column tails of half a panel or more, and an edge kernel with all tiles along M.
 *  - Four-row ZA moves (MOVA vg4) for the C tiles and for the transposition of A, core-side L2 prefetch.
 *  - One thread per SME unit (Apple: per performance cluster) when the problem is large enough.
 *  - Problems of fewer than S2_SMALL_WORK multiply-adds run in plain loops, below the cost of streaming mode.
 */

#include "sme2_gemm_tile.h"
#include "sme2_gemm_detect.h"
#include <stdlib.h>
#include <string.h>
#if defined(__APPLE__)
#include <sys/sysctl.h>
#endif

/* The file is built with the SME flags of its target (no SME2); the functions here enable SME2 themselves. */
#pragma clang attribute push(__attribute__((target("sme2"))), apply_to = function)

#define S2_CH (S2_NT * S2_VL) /* columns of A moved through ZA per transposition pass */
#define S2_TN2 (S2_NT / 2)    /* tiles along N of the half-width tail kernel */
#define S2_W2 (S2_TN2 * S2_VL)
#define S2_PF_ROWS_B 8        /* B rows ahead for the core prefetch of B strips */

/* Blocking model parameters (Apple M4: SME load bandwidth holds up to 8 MB, 16 KB pages). */
#ifndef S2_L2_BYTES
#define S2_L2_BYTES (8L << 20)
#endif
#ifndef S2_PAGE_BYTES
#define S2_PAGE_BYTES 16384L
#endif
#ifndef S2_TLB_ENTRIES
#define S2_TLB_ENTRIES 160
#endif

/* A row of B for TN >= 4 (up to eight vectors): from the packed panel, or from the source and written to it. */
S2_INL void s2_b_row(const int TN, const int ONLINE, const FLOAT *src, FLOAT *dst, int ncols, s2_v4 *b0,
                     s2_v4 *b1) S2_S {
  if (ONLINE) {
    *b0 = ncols >= 4 * S2_VL ? s2_ld4(src) : s2_ld4n(src, s2_clip(ncols, 4 * S2_VL));
    s2_st4(dst, *b0);
    if (TN == 8) {
      *b1 = ncols >= 8 * S2_VL ? s2_ld4(src + 4 * S2_VL) : s2_ld4n(src + 4 * S2_VL, s2_clip(ncols - 4 * S2_VL, 4 * S2_VL));
      s2_st4(dst + 4 * S2_VL, *b1);
    } else {
      *b1 = *b0;
    }
  } else {
    *b0 = s2_ld4(dst);
    *b1 = TN == 8 ? s2_ld4(dst + 4 * S2_VL) : *b0;
  }
}

/*
 * Micro-kernel (paper Alg. 1): TM x TN tiles. Ar holds TM packed panels of VL rows (stride astr), Br is the
 * kb x (TN*VL) packed panel of B. ONLINE: B comes from the source Bs and is written to Br as a side effect.
 */
S2_INL void s2_kernel(const int TM, const int TN, const int ONLINE, int kb, const FLOAT *Ar, BLASLONG astr, FLOAT *Br,
                      const FLOAT *Bs, BLASLONG ldb, int ncols, FLOAT *C, BLASLONG ldc, int mrows, int cmode,
                      FLOAT beta, int pfb) S2_S S2_ZA {
  const int NR = TN * S2_VL;
  s2_c_load(TM, TN, C, ldc, mrows, ncols, cmode, beta);
  const svbool_t pt = S2_PT();
  int k = 0;
  if (TN >= 4) {
    for (; k + 4 <= kb; k += 4) {
      if (ONLINE && pfb && k + S2_PF_ROWS_B + 4 <= kb)
        for (int u = 0; u < 4; ++u) s2_pf_l2(Bs + (BLASLONG)(k + S2_PF_ROWS_B + u) * ldb, ncols * (int)sizeof(FLOAT));
      const s2_v4 a0 = s2_ld4(Ar + (BLASLONG)k * S2_VL);
      const s2_v4 a1 = TM > 1 ? s2_ld4(Ar + astr + (BLASLONG)k * S2_VL) : a0;
      S2_UNROLL for (int U = 0; U < 4; ++U) {
        s2_v4 b0, b1;
        s2_b_row(TN, ONLINE, Bs + (BLASLONG)(k + U) * ldb, Br + (BLASLONG)(k + U) * NR, ncols, &b0, &b1);
        S2_UNROLL for (int P = 0; P < TM; ++P) {
          const s2_v a = s2_get4(P == 0 ? a0 : a1, U);
          S2_UNROLL for (int Q = 0; Q < TN; ++Q) s2_mopa(P * TN + Q, pt, a, s2_pick(b0, b1, Q));
        }
      }
    }
    for (; k < kb; ++k) {
      s2_v4 b0, b1;
      s2_b_row(TN, ONLINE, Bs + (BLASLONG)k * ldb, Br + (BLASLONG)k * NR, ncols, &b0, &b1);
      S2_UNROLL for (int P = 0; P < TM; ++P) {
        const s2_v a = S2_LD1(pt, Ar + P * astr + (BLASLONG)k * S2_VL);
        S2_UNROLL for (int Q = 0; Q < TN; ++Q) s2_mopa(P * TN + Q, pt, a, s2_pick(b0, b1, Q));
      }
    }
  } else {
    for (; k + 4 <= kb; k += 4) {
      if (ONLINE && pfb && k + S2_PF_ROWS_B + 4 <= kb)
        for (int u = 0; u < 4; ++u) s2_pf_l2(Bs + (BLASLONG)(k + S2_PF_ROWS_B + u) * ldb, ncols * (int)sizeof(FLOAT));
      s2_v4 bb0, bb1;
      FLOAT *dst = Br + (BLASLONG)k * NR;
      if (ONLINE) {
        const svbool_t p0 = S2_PW(0, ncols), p1 = S2_PW(S2_VL, ncols);
        const FLOAT *s = Bs + (BLASLONG)k * ldb;
        if (TN == 1) {
          bb0 = S2_CREATE4(S2_LD1(p0, s), S2_LD1(p0, s + ldb), S2_LD1(p0, s + 2 * ldb), S2_LD1(p0, s + 3 * ldb));
          bb1 = bb0;
        } else {
          bb0 = S2_CREATE4(S2_LD1(p0, s), S2_LD1(p1, s + S2_VL), S2_LD1(p0, s + ldb), S2_LD1(p1, s + ldb + S2_VL));
          bb1 = S2_CREATE4(S2_LD1(p0, s + 2 * ldb), S2_LD1(p1, s + 2 * ldb + S2_VL), S2_LD1(p0, s + 3 * ldb),
                           S2_LD1(p1, s + 3 * ldb + S2_VL));
        }
        s2_st4(dst, bb0);
        if (TN == 2) s2_st4(dst + 4 * S2_VL, bb1);
      } else {
        bb0 = s2_ld4(dst);
        bb1 = TN == 2 ? s2_ld4(dst + 4 * S2_VL) : bb0;
      }
      /* Up to four A panels at a time, FMOPAs ordered so that consecutive ones hit different tiles. */
      S2_UNROLL for (int H = 0; H < (TM + 3) / 4; ++H) {
        const int P0 = H * 4, NP = TM - P0 < 4 ? TM - P0 : 4;
        const s2_v4 a0 = s2_ld4(Ar + P0 * astr + (BLASLONG)k * S2_VL);
        const s2_v4 a1 = NP > 1 ? s2_ld4(Ar + (P0 + 1) * astr + (BLASLONG)k * S2_VL) : a0;
        const s2_v4 a2 = NP > 2 ? s2_ld4(Ar + (P0 + 2) * astr + (BLASLONG)k * S2_VL) : a0;
        const s2_v4 a3 = NP > 3 ? s2_ld4(Ar + (P0 + 3) * astr + (BLASLONG)k * S2_VL) : a0;
        S2_UNROLL for (int U = 0; U < 4; ++U) {
          S2_UNROLL for (int PP = 0; PP < NP; ++PP) {
            const s2_v a = s2_get4(PP == 0 ? a0 : PP == 1 ? a1 : PP == 2 ? a2 : a3, U);
            S2_UNROLL for (int Q = 0; Q < TN; ++Q) s2_mopa((P0 + PP) * TN + Q, pt, a, s2_pick(bb0, bb1, U * TN + Q));
          }
        }
      }
    }
    for (; k < kb; ++k) {
      FLOAT *dst = Br + (BLASLONG)k * NR;
      s2_v b0, b1;
      if (ONLINE) {
        const FLOAT *s = Bs + (BLASLONG)k * ldb;
        b0 = S2_LD1(S2_PW(0, ncols), s);
        S2_ST1(pt, dst, b0);
        if (TN == 2) {
          b1 = S2_LD1(S2_PW(S2_VL, ncols), s + S2_VL);
          S2_ST1(pt, dst + S2_VL, b1);
        } else {
          b1 = b0;
        }
      } else {
        b0 = S2_LD1(pt, dst);
        b1 = TN == 2 ? S2_LD1(pt, dst + S2_VL) : b0;
      }
      S2_UNROLL for (int P = 0; P < TM; ++P) {
        const s2_v a = S2_LD1(pt, Ar + P * astr + (BLASLONG)k * S2_VL);
        s2_mopa(P * TN, pt, a, b0);
        if (TN == 2) s2_mopa(P * TN + 1, pt, a, b1);
      }
    }
  }
  s2_c_store(TM, TN, C, ldc, mrows, ncols);
}

/* One noinline function per kernel shape; the table index is TM - 1. */
typedef void (*s2_kfn)(int, const FLOAT *, BLASLONG, FLOAT *, const FLOAT *, BLASLONG, int, FLOAT *, BLASLONG, int,
                       int, FLOAT, int) S2_S S2_ZA;
#define S2_DEF_KERNEL(TM, TN, ON)                                                                                \
  S2_NOINL void s2_k_##TM##_##TN##_##ON(int kb, const FLOAT *Ar, BLASLONG astr, FLOAT *Br, const FLOAT *Bs,       \
                                        BLASLONG ldb, int ncols, FLOAT *C, BLASLONG ldc, int mrows, int cmode,    \
                                        FLOAT beta, int pfb) S2_S S2_ZA {                                         \
    s2_kernel(TM, TN, ON, kb, Ar, astr, Br, Bs, ldb, ncols, C, ldc, mrows, cmode, beta, pfb);                     \
  }
#define S2_DEF_BOTH(TM, TN) S2_DEF_KERNEL(TM, TN, 0) S2_DEF_KERNEL(TM, TN, 1)
#ifdef DOUBLE
S2_DEF_BOTH(1, 8)
S2_DEF_BOTH(1, 4)
S2_DEF_BOTH(2, 4)
#else
S2_DEF_BOTH(1, 4)
S2_DEF_BOTH(1, 2)
S2_DEF_BOTH(2, 2)
#endif
S2_DEF_BOTH(1, 1)
S2_DEF_BOTH(2, 1)
S2_DEF_BOTH(3, 1)
S2_DEF_BOTH(4, 1)
#ifdef DOUBLE
S2_DEF_BOTH(5, 1)
S2_DEF_BOTH(6, 1)
S2_DEF_BOTH(7, 1)
S2_DEF_BOTH(8, 1)
#endif

/*
 * Transposition in ZA (paper Fig. 6): `rows` (<= VL) rows of kb contiguous elements (row stride ld) go in through
 * horizontal slices of all tiles, CH columns at a time, and leave through vertical slices: dst[k * w + r].
 * w == VL gives the A panels; w > VL writes one VL-column group of a wider B panel.
 */
S2_INL void s2_pack_t(int rows, int kb, const FLOAT *src, BLASLONG ld, FLOAT *dst, const int w, FLOAT alpha,
                      const int scale, int pf) S2_S S2_ZA {
  const svbool_t pt = S2_PT();
  if (rows <= 0) {
    for (int k = 0; k < kb; ++k) S2_ST1(pt, dst + (BLASLONG)k * w, S2_ZEROV());
    return;
  }
  for (int k0 = 0; k0 < kb; k0 += S2_CH) {
    const int kn = kb - k0 < S2_CH ? kb - k0 : S2_CH;
    if (rows < S2_VL) svzero_za();
    int r = 0;
    if (kn >= S2_CH)
      for (; r + 4 <= rows; r += 4) {
        const FLOAT *q = src + (BLASLONG)r * ld + k0;
        if (pf && k0 + S2_CH < kb)
          for (int u = 0; u < 4; ++u) s2_pf_l2(q + u * ld + S2_CH, S2_CH * (int)sizeof(FLOAT));
        s2_rows_in4(q, ld, r);
      }
    for (; r < rows; ++r) {
      const FLOAT *q = src + (BLASLONG)r * ld + k0;
      if (pf && k0 + S2_CH < kb) s2_pf_l2(q + S2_CH, S2_CH * (int)sizeof(FLOAT));
      S2_UNROLL for (int G = 0; G < S2_NT / 4; ++G) {
        s2_v4 v = kn >= S2_CH ? s2_ld4(q + G * 4 * S2_VL) : s2_ld4n(q + G * 4 * S2_VL, s2_clip(kn - G * 4 * S2_VL, 4 * S2_VL));
        S2_UNROLL for (int Q = 0; Q < 4; ++Q) {
          s2_wrh(G * 4 + Q, r, s2_get4(v, Q));
        }
      }
    }
    S2_UNROLL for (int t = 0; t < S2_NT; ++t) {
      const int kt = kn - t * S2_VL;
      for (int c = 0; c < S2_VL && c < kt; c += 4) {
        s2_v4 v = s2_rdv4(t, c);
        if (scale) /* alpha on the way out, so that the four-row path serves scaled packing too */
          v = S2_CREATE4(S2_MUL(S2_GET4(v, 0), alpha), S2_MUL(S2_GET4(v, 1), alpha), S2_MUL(S2_GET4(v, 2), alpha),
                         S2_MUL(S2_GET4(v, 3), alpha));
        FLOAT *d = dst + (BLASLONG)(k0 + t * S2_VL + c) * w;
        if (w == S2_VL) {
          if (kt - c >= 4) s2_st4(d, v);
          else s2_st4n(d, v, (int64_t)(kt - c) * S2_VL);
        } else {
          S2_UNROLL for (int u = 0; u < 4; ++u)
            if (c + u < kt) S2_ST1(pt, d + (BLASLONG)u * w, s2_get4(v, u));
        }
      }
    }
  }
}

/* Specializations (as MTGEMM-A's templates): A panels, A panels scaled by alpha, and B panel groups. */
S2_NOINL void s2_pack_a_rows(int rows, int kb, const FLOAT *src, BLASLONG ld, FLOAT *dst) S2_S S2_ZA {
  s2_pack_t(rows, kb, src, ld, dst, S2_VL, (FLOAT)1, 0, 1);
}
S2_NOINL void s2_pack_a_rows_scaled(int rows, int kb, const FLOAT *src, BLASLONG ld, FLOAT *dst, FLOAT alpha) S2_S S2_ZA {
  s2_pack_t(rows, kb, src, ld, dst, S2_VL, alpha, 1, 1);
}
S2_NOINL void s2_pack_b_group(int rows, int kb, const FLOAT *src, BLASLONG ld, FLOAT *dst, int w) S2_S S2_ZA {
  if (w == S2_NR) s2_pack_t(rows, kb, src, ld, dst, S2_NR, (FLOAT)1, 0, 1);
  else if (w == S2_W2) s2_pack_t(rows, kb, src, ld, dst, S2_W2, (FLOAT)1, 0, 1);
  else s2_pack_t(rows, kb, src, ld, dst, S2_VL, (FLOAT)1, 0, 1);
}

/*
 * A with contiguous columns (element (r, k) at r + k * ld) into panels of VL rows: a copy, four panels and four
 * depth steps at a time (four x4 loads, one x4 store per panel), since single dependent loads run at the SME
 * unit's load latency.
 */
S2_INL void s2_pack_a_cols_t(int mb, int kb, const FLOAT *A, BLASLONG ld, FLOAT alpha, FLOAT *Ac, const int scale)
    S2_S {
  const svbool_t pt = S2_PT();
  for (int p0 = 0; p0 < mb; p0 += 4 * S2_VL) {
    const int64_t rows = s2_clip(mb - p0, 4 * S2_VL);
    const int np = (int)((rows + S2_VL - 1) / S2_VL);
    const svcount_t pc = S2_CW(0, rows);
    const FLOAT *s = A + p0;
    FLOAT *d = Ac + (BLASLONG)p0 * kb;
    int k = 0;
    for (; k + 4 <= kb; k += 4) {
      if (k + 8 + 4 <= kb)
        for (int u = 0; u < 4; ++u) s2_pf_l2(s + (BLASLONG)(k + 8 + u) * ld, (int)rows * (int)sizeof(FLOAT));
      s2_v4 x0 = S2_LD4(pc, s + (BLASLONG)k * ld), x1 = S2_LD4(pc, s + (BLASLONG)(k + 1) * ld);
      s2_v4 x2 = S2_LD4(pc, s + (BLASLONG)(k + 2) * ld), x3 = S2_LD4(pc, s + (BLASLONG)(k + 3) * ld);
      S2_UNROLL for (int q = 0; q < 4; ++q) {
        if (q >= np) break;
        s2_v4 v = S2_CREATE4(s2_get4(x0, q), s2_get4(x1, q), s2_get4(x2, q), s2_get4(x3, q));
        if (scale)
          v = S2_CREATE4(S2_MUL(S2_GET4(v, 0), alpha), S2_MUL(S2_GET4(v, 1), alpha), S2_MUL(S2_GET4(v, 2), alpha),
                         S2_MUL(S2_GET4(v, 3), alpha));
        s2_st4(d + (BLASLONG)q * S2_VL * kb + (BLASLONG)k * S2_VL, v);
      }
    }
    for (; k < kb; ++k) {
      const s2_v4 x = S2_LD4(pc, s + (BLASLONG)k * ld);
      S2_UNROLL for (int q = 0; q < 4; ++q) {
        if (q >= np) break;
        s2_v v = s2_get4(x, q);
        if (scale) v = S2_MUL(v, alpha);
        S2_ST1(pt, d + (BLASLONG)q * S2_VL * kb + (BLASLONG)k * S2_VL, v);
      }
    }
  }
}
S2_NOINL void s2_pack_a_cols(int mb, int kb, const FLOAT *A, BLASLONG ld, FLOAT alpha, FLOAT *Ac)
    S2_S {
  if (alpha != (FLOAT)1) s2_pack_a_cols_t(mb, kb, A, ld, alpha, Ac, 1);
  else s2_pack_a_cols_t(mb, kb, A, ld, alpha, Ac, 0);
}

/* B panel of `ncols` (<= w) columns from a B with contiguous columns (element (k, n) at k + n * ld). */
S2_INL void s2_pack_b_cols(int kb, int ncols, int w, const FLOAT *Bs, BLASLONG ld, FLOAT *Bd) S2_S S2_ZA {
  for (int g = 0; g < w; g += S2_VL)
    s2_pack_b_group((int)s2_clip(ncols - g, S2_VL), kb, Bs + (BLASLONG)g * ld, ld, Bd + g, w);
}

/* ---- blocking model (paper eqs. 1-3) ---- */

typedef struct { int mc, nc, kc; } s2_blocking;

S2_INL int s2_round_up(int x, int m) { return (x + m - 1) / m * m; }

static int s2_model_kc_max(int mr, int nr) {
  int best = 16;
  for (int kc = 16; kc <= 1 << 16; kc += 16) {
    const int ta = (int)(((long)mr * kc * sizeof(FLOAT) + S2_PAGE_BYTES - 1) / S2_PAGE_BYTES) + 1;
    const int tb = (int)(((long)nr * kc * sizeof(FLOAT) + S2_PAGE_BYTES - 1) / S2_PAGE_BYTES) + 1;
    if (ta + 2 * tb + mr < S2_TLB_ENTRIES) best = kc;
    else break;
  }
  return best;
}

static int s2_balance(int total, int block, int step) {
  if (block >= total) return s2_round_up(total, step);
  const int nblk = (total + block - 1) / block;
  const int b = s2_round_up((total + nblk - 1) / nblk, step);
  return b < block ? b : block;
}

static s2_blocking s2_model(int M, int N, int K, int mr, int nr) {
  const long budget = S2_L2_BYTES / (long)sizeof(FLOAT);
  int kcap = s2_model_kc_max(mr, nr);
  if (kcap > s2_round_up(K, 16)) kcap = s2_round_up(K, 16);
  const int mcap = s2_round_up(M, mr), ncap = s2_round_up(N, nr);
  double best = -1;
  s2_blocking r = {mr, nr, 16};
  for (int kc = 16; kc <= kcap; kc += 16)
    for (int mc = mr; mc <= mcap; mc += mr) {
      const long room = budget - (long)mc * kc;
      if (room <= 0) break;
      long nc = room / (2L * kc + 2L * mc);
      nc = nc / nr * nr;
      if (nc > ncap) nc = ncap;
      if (nc < nr) break;
      const double cmr = 2.0 * mc * nc * kc / ((double)mc * kc + (double)kc * nc + 2.0 * mc * nc);
      if (cmr > best * (1 + 1e-9)) {
        best = cmr;
        r.mc = mc;
        r.nc = (int)nc;
        r.kc = kc;
      }
    }
  r.kc = s2_balance(K, r.kc, 16);
  r.mc = s2_balance(M, r.mc, mr);
  r.nc = s2_balance(N, r.nc, nr);
  return r;
}

#if defined(_MSC_VER) && !defined(__clang__)
#define S2_TLS __declspec(thread)
#else
#define S2_TLS _Thread_local
#endif

/* The search takes up to tens of microseconds for tall problems, so each thread remembers its last shapes. */
static s2_blocking s2_blocking_for(int M, int N, int K) {
  static S2_TLS struct { int M, N, K; s2_blocking b; } cache[8];
  static S2_TLS int next;
  for (int i = 0; i < 8; ++i)
    if (cache[i].M == M && cache[i].N == N && cache[i].K == K) return cache[i].b;
  s2_blocking b = s2_model(M, N, K, S2_VL, S2_NR);
  b.mc = s2_round_up(b.mc < S2_VL ? S2_VL : b.mc, S2_VL);
  b.nc = s2_round_up(b.nc < S2_NR ? S2_NR : b.nc, S2_NR);
  if (b.kc < 1) b.kc = 1;
  cache[next].M = M;
  cache[next].N = N;
  cache[next].K = K;
  cache[next].b = b;
  next = (next + 1) & 7;
  return b;
}

/* ---- driver ---- */

/*
 * Row-major core problem C = alpha A B + beta C, C of M x N (row stride ldc).
 * a_cols: A has contiguous columns (element (i, k) at i + k * lda), else contiguous rows (i * lda + k).
 * b_cols: B has contiguous columns (element (k, j) at k + j * ldb), else contiguous rows (k * ldb + j).
 */
typedef struct {
  int M, N, K;
  FLOAT alpha, beta;
  const FLOAT *A, *B;
  FLOAT *C;
  BLASLONG lda, ldb, ldc;
  int a_cols, b_cols;
} s2_job;

static const s2_kfn s2_k_main[2] = {
#ifdef DOUBLE
    s2_k_1_8_0, s2_k_1_8_1};
static const s2_kfn s2_k_mid[2][2] = {{s2_k_1_4_0, s2_k_2_4_0}, {s2_k_1_4_1, s2_k_2_4_1}};
static const s2_kfn s2_k_edge[2][S2_NT] = {
    {s2_k_1_1_0, s2_k_2_1_0, s2_k_3_1_0, s2_k_4_1_0, s2_k_5_1_0, s2_k_6_1_0, s2_k_7_1_0, s2_k_8_1_0},
    {s2_k_1_1_1, s2_k_2_1_1, s2_k_3_1_1, s2_k_4_1_1, s2_k_5_1_1, s2_k_6_1_1, s2_k_7_1_1, s2_k_8_1_1}};
#else
    s2_k_1_4_0, s2_k_1_4_1};
static const s2_kfn s2_k_mid[2][2] = {{s2_k_1_2_0, s2_k_2_2_0}, {s2_k_1_2_1, s2_k_2_2_1}};
static const s2_kfn s2_k_edge[2][S2_NT] = {{s2_k_1_1_0, s2_k_2_1_0, s2_k_3_1_0, s2_k_4_1_0},
                                           {s2_k_1_1_1, s2_k_2_1_1, s2_k_3_1_1, s2_k_4_1_1}};
#endif

__arm_new("za") __arm_locally_streaming static void s2_drive(const s2_job *jb, s2_blocking blk, FLOAT *Ac, FLOAT *Bc) {
  const int M = jb->M, N = jb->N, K = jb->K, mc = blk.mc, nc = blk.nc, kc = blk.kc;
  const BLASLONG lda = jb->lda, ldb = jb->ldb, ldc = jb->ldc;
  const int online = !jb->b_cols;
  for (int i = 0; i < M; i += mc) {
    const int mb = M - i < mc ? M - i : mc;
    for (int k = 0; k < K; k += kc) {
      const int kb = K - k < kc ? K - k : kc;
      if (jb->a_cols) {
        s2_pack_a_cols(mb, kb, jb->A + i + (BLASLONG)k * lda, lda, jb->alpha, Ac);
      } else {
        const FLOAT *Ap = jb->A + (BLASLONG)i * lda + k;
        for (int p0 = 0; p0 < mb; p0 += S2_VL) {
          const int rows = mb - p0 < S2_VL ? mb - p0 : S2_VL;
          if (jb->alpha != (FLOAT)1)
            s2_pack_a_rows_scaled(rows, kb, Ap + (BLASLONG)p0 * lda, lda, Ac + (BLASLONG)p0 * kb, jb->alpha);
          else
            s2_pack_a_rows(rows, kb, Ap + (BLASLONG)p0 * lda, lda, Ac + (BLASLONG)p0 * kb);
        }
      }
      const int cmode = k > 0 ? S2_LOAD : (jb->beta == (FLOAT)0 ? S2_ZERO : (jb->beta == (FLOAT)1 ? S2_LOAD : S2_SCALE));
      for (int j = 0; j < N; j += nc) {
        const int nb = N - j < nc ? N - j : nc, nmain = nb / S2_NR * S2_NR;
        const FLOAT *Bs = jb->b_cols ? jb->B + (BLASLONG)j * ldb + k : jb->B + (BLASLONG)k * ldb + j;
        FLOAT *Cb = jb->C + (BLASLONG)i * ldc + j;
        /* Column tails wider than one vector go to the half-width kernel with two tile rows (predicated when
           narrower than it), the last VL or fewer columns to the edge kernel. */
        const int ecol = nb - nmain > S2_VL ? nmain + (nb - nmain - 1) / S2_W2 * S2_W2 + ((nb - nmain - 1) % S2_W2 >= S2_VL ? S2_W2 : 0) : nmain;
        if (!online) {
          for (int jj = 0; jj < nmain; jj += S2_NR)
            s2_pack_b_cols(kb, S2_NR, S2_NR, Bs + (BLASLONG)jj * ldb, ldb, Bc + (BLASLONG)jj * kb);
          for (int jj = nmain; jj < ecol; jj += S2_W2)
            s2_pack_b_cols(kb, nb - jj < S2_W2 ? nb - jj : S2_W2, S2_W2, Bs + (BLASLONG)jj * ldb, ldb, Bc + (BLASLONG)jj * kb);
          for (int jj = ecol; jj < nb; jj += S2_VL)
            s2_pack_b_cols(kb, nb - jj < S2_VL ? nb - jj : S2_VL, S2_VL, Bs + (BLASLONG)jj * ldb, ldb, Bc + (BLASLONG)jj * kb);
        }
        for (int ii = 0; ii < mb; ii += S2_VL) {
          const int rows = mb - ii < S2_VL ? mb - ii : S2_VL;
          for (int jj = 0; jj < nmain; jj += S2_NR) {
            FLOAT *Cp = Cb + (BLASLONG)ii * ldc + jj;
            FLOAT *Br = Bc + (BLASLONG)jj * kb;
            const FLOAT *Cn = jj + S2_NR < nmain ? Cp + S2_NR : Cb + (BLASLONG)(ii + S2_VL) * ldc;
            if (jj + S2_NR < nmain || ii + S2_VL < mb) /* next C tile towards L2 while this one computes */
              for (int r = 0; r < rows; ++r) s2_pf_l2(Cn + (BLASLONG)r * ldc, S2_NR * (int)sizeof(FLOAT));
            const int on = online && ii == 0;
            s2_k_main[on](kb, Ac + (BLASLONG)ii * kb, (BLASLONG)S2_VL * kb, Br, on ? Bs + jj : NULL, ldb, S2_NR, Cp,
                          ldc, rows, cmode, jb->beta, 1);
          }
        }
        for (int jj = nmain; jj < ecol; jj += S2_W2) {
          const int cols = nb - jj < S2_W2 ? nb - jj : S2_W2;
          for (int ii = 0; ii < mb; ii += 2 * S2_VL) {
            const int rows = mb - ii < 2 * S2_VL ? mb - ii : 2 * S2_VL, on = online && ii == 0;
            s2_k_mid[on][(rows + S2_VL - 1) / S2_VL - 1](kb, Ac + (BLASLONG)ii * kb, (BLASLONG)S2_VL * kb,
                                                        Bc + (BLASLONG)jj * kb, on ? Bs + jj : NULL, ldb, cols,
                                                        Cb + (BLASLONG)ii * ldc + jj, ldc, rows, cmode, jb->beta, 1);
          }
        }
        if (ecol < nb) /* edge kernel: all tiles stacked along M, one VL-wide column panel */
          for (int ii = 0; ii < mb; ii += S2_NT * S2_VL) {
            const int rows = mb - ii < S2_NT * S2_VL ? mb - ii : S2_NT * S2_VL;
            for (int jj = ecol; jj < nb; jj += S2_VL) {
              const int cols = nb - jj < S2_VL ? nb - jj : S2_VL, on = online && ii == 0;
              s2_k_edge[on][(rows + S2_VL - 1) / S2_VL - 1](kb, Ac + (BLASLONG)ii * kb, (BLASLONG)S2_VL * kb,
                                                           Bc + (BLASLONG)jj * kb, on ? Bs + jj : NULL, ldb, cols,
                                                           Cb + (BLASLONG)ii * ldc + jj, ldc, rows, cmode, jb->beta, 1);
            }
          }
      }
    }
  }
}

#define S2_THREAD_WORK (1L << 20) /* largest per-thread packing buffer; larger problems use blas_memory_alloc */

/* Bytes of packed A and B for blocking b (the model keeps this near S2_L2_BYTES, well inside BUFFER_SIZE). */
static size_t s2_work_bytes(s2_blocking b, size_t *a_bytes) {
  *a_bytes = ((size_t)s2_round_up(b.mc, 4 * S2_VL) * b.kc * sizeof(FLOAT) + 127) & ~(size_t)127;
  return *a_bytes + (size_t)b.nc * b.kc * sizeof(FLOAT) + 128;
}

/*
 * Per-thread packing buffer for problems up to S2_THREAD_WORK bytes, grown on demand and freed when the thread
 * exits: the per-call cost of blas_memory_alloc, and of packing on the stack next to the core's own data, is
 * significant against a small GEMM (fp32 20^3 column-major: 568 ns on the stack, 336 ns here).
 */
static pthread_key_t s2_work_key;
static pthread_once_t s2_work_once = PTHREAD_ONCE_INIT;
static void s2_work_free(void *p) { free(p); }
static void s2_work_key_init(void) { pthread_key_create(&s2_work_key, s2_work_free); }
static void *s2_thread_work(size_t bytes) {
  pthread_once(&s2_work_once, s2_work_key_init);
  size_t *p = (size_t *)pthread_getspecific(s2_work_key);
  if (p == NULL || p[0] < bytes) {
    size_t cap = p ? p[0] : 0;
    while (cap < bytes) cap = cap ? 2 * cap : 65536;
    free(p);
    void *q = NULL;
    if (posix_memalign(&q, 16384, cap + 128)) q = NULL;
    p = (size_t *)q;
    if (p) p[0] = cap;
    pthread_setspecific(s2_work_key, p);
    if (!p) return NULL;
  }
  return (char *)p + 128;
}

/* Runs one job on a workspace of at least s2_work_bytes() bytes (NULL: a per-thread or blas_memory_alloc buffer). */
static void s2_run(const s2_job *jb, void *work) {
  if (jb->M <= 0 || jb->N <= 0) return;
  const s2_blocking b = s2_blocking_for(jb->M, jb->N, jb->K);
  size_t a_bytes;
  const size_t bytes = s2_work_bytes(b, &a_bytes);
  void *buffer = NULL, *heap = NULL;
  if (work == NULL && bytes <= S2_THREAD_WORK) work = s2_thread_work(bytes);
  if (work == NULL && (buffer = blas_memory_alloc(0)) != NULL) work = (char *)buffer + GEMM_OFFSET_A;
  if (work == NULL && (work = heap = malloc(bytes + 128)) == NULL) return;
  char *base = (char *)(((uintptr_t)work + 127) & ~(uintptr_t)127);
  s2_drive(jb, b, (FLOAT *)base, (FLOAT *)(base + a_bytes));
  if (buffer) blas_memory_free(buffer);
  free(heap);
}

/* ---- runtime checks, threads, entry ---- */

/*
 * SME units that threads can use: Apple M4 and M5 have one per cluster of the fastest cores. M5 Pro and M5 Max,
 * whose second tier is "Performance" rather than "Efficiency", have two per "Super" cluster (third-party
 * measurement, not verified here). Elsewhere assume one.
 */
static int s2_units(void) {
  static int n = 0;
  if (n == 0) {
    int u = 1;
#if defined(__APPLE__)
    int cpus = 0, per = 0;
    char name[32] = "";
    size_t len = sizeof(int);
    if (sysctlbyname("hw.perflevel0.physicalcpu", &cpus, &len, NULL, 0) == 0 &&
        (len = sizeof(int), sysctlbyname("hw.perflevel0.cpusperl2", &per, &len, NULL, 0) == 0) && per > 0 &&
        cpus / per > 1)
      u = cpus / per;
    len = sizeof(name) - 1;
    if (sysctlbyname("hw.perflevel1.name", name, &len, NULL, 0) == 0 && strcmp(name, "Performance") == 0) u *= 2;
#endif
    n = u;
  }
  return n;
}

#ifdef SMP
static int s2_routine(blas_arg_t *args, BLASLONG *range_m, BLASLONG *range_n, FLOAT *sa, FLOAT *sb, BLASLONG pos) {
  s2_run((const s2_job *)args->common, sa);
  return 0;
}
#endif

#define S2_MAX_PARTS 4
#ifndef S2_SMALL_WORK
#define S2_SMALL_WORK 3000. /* multiply-adds below which s2_small is faster (Apple M4, both precisions) */
#endif

/* Column-major GEMM in plain loops (vectorized for NEON) for problems smaller than the cost of streaming mode. */
static void s2_small(int trans_a, int trans_b, BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, const FLOAT *a,
                     BLASLONG lda, const FLOAT *b, BLASLONG ldb, FLOAT beta, FLOAT *c, BLASLONG ldc) {
  const BLASLONG sbl = trans_b ? ldb : 1, sbj = trans_b ? 1 : ldb; /* stride of B along l and along j */
  for (BLASLONG j = 0; j < n; ++j) {
    FLOAT *cj = c + j * ldc;
    if (beta == (FLOAT)0)
      for (BLASLONG i = 0; i < m; ++i) cj[i] = 0;
    else if (beta != (FLOAT)1)
      for (BLASLONG i = 0; i < m; ++i) cj[i] *= beta;
  }
  if (trans_a) {
    for (BLASLONG j = 0; j < n; ++j)
      for (BLASLONG i = 0; i < m; ++i) {
        const FLOAT *ai = a + i * lda, *bj = b + j * sbj;
        FLOAT s = 0;
        for (BLASLONG l = 0; l < k; ++l) s += ai[l] * bj[l * sbl];
        c[i + j * ldc] += alpha * s;
      }
    return;
  }
  BLASLONG j = 0;
  for (; j + 4 <= n; j += 4) {
    FLOAT *c0 = c + j * ldc, *c1 = c0 + ldc, *c2 = c1 + ldc, *c3 = c2 + ldc;
    const FLOAT *bj = b + j * sbj;
    for (BLASLONG l = 0; l < k; ++l) {
      const FLOAT *al = a + l * lda, *bl = bj + l * sbl;
      const FLOAT t0 = alpha * bl[0], t1 = alpha * bl[sbj], t2 = alpha * bl[2 * sbj], t3 = alpha * bl[3 * sbj];
      for (BLASLONG i = 0; i < m; ++i) {
        const FLOAT x = al[i];
        c0[i] += t0 * x;
        c1[i] += t1 * x;
        c2[i] += t2 * x;
        c3[i] += t3 * x;
      }
    }
  }
  for (; j < n; ++j) {
    FLOAT *cj = c + j * ldc;
    for (BLASLONG l = 0; l < k; ++l) {
      const FLOAT t = alpha * b[j * sbj + l * sbl], *al = a + l * lda;
      for (BLASLONG i = 0; i < m; ++i) cj[i] += t * al[i];
    }
  }
}

/* Column-major C = alpha op(A) op(B) + beta C, solved as the row-major C^T = alpha op(B)^T op(A)^T + beta C^T. */
static void s2_gemm(int trans_a, int trans_b, BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, const FLOAT *a,
                    BLASLONG lda, const FLOAT *b, BLASLONG ldb, FLOAT beta, FLOAT *c, BLASLONG ldc) {
  if (m <= 0 || n <= 0) return;
  if ((double)m * n * k < S2_SMALL_WORK) {
    s2_small(trans_a, trans_b, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    return;
  }
  if (k <= 0 || alpha == (FLOAT)0) {
    for (BLASLONG j = 0; j < n; ++j)
      for (BLASLONG i = 0; i < m; ++i) c[i + j * ldc] = beta == (FLOAT)0 ? (FLOAT)0 : beta * c[i + j * ldc];
    return;
  }
  s2_job jb;
  jb.M = (int)n;
  jb.N = (int)m;
  jb.K = (int)k;
  jb.alpha = alpha;
  jb.beta = beta;
  jb.A = b;
  jb.lda = ldb;
  jb.a_cols = trans_b;
  jb.B = a;
  jb.ldb = lda;
  jb.b_cols = trans_a;
  jb.C = c;
  jb.ldc = ldc;

  int parts = 1;
#ifdef SMP
  /* One thread per SME unit from 2^22 multiply-adds, each part at least four panels wide along the split. */
  const int big = jb.M > jb.N ? jb.M : jb.N;
  if ((double)jb.M * jb.N * jb.K >= (double)(1 << 22)) {
    parts = s2_units();
    const int avail = num_cpu_avail(3);
    if (parts > avail) parts = avail;
    if (parts > big / (4 * S2_VL)) parts = big / (4 * S2_VL);
    if (parts > S2_MAX_PARTS) parts = S2_MAX_PARTS;
    if (parts < 1) parts = 1;
  }
  if (parts > 1) {
    s2_job pj[S2_MAX_PARTS];
    blas_arg_t args[S2_MAX_PARTS];
    blas_queue_t queue[S2_MAX_PARTS];
    const int split_m = jb.M >= jb.N;
    const int total = split_m ? jb.M : jb.N;
    const int step = s2_round_up((total + parts - 1) / parts, 4 * S2_VL);
    int used = 0;
    for (int p = 0, off = 0; p < parts && off < total; ++p, off += step) {
      const int len = total - off < step ? total - off : step;
      pj[p] = jb;
      if (split_m) {
        pj[p].M = len;
        pj[p].A = jb.a_cols ? jb.A + off : jb.A + (BLASLONG)off * jb.lda;
        pj[p].C = jb.C + (BLASLONG)off * jb.ldc;
      } else {
        pj[p].N = len;
        pj[p].B = jb.b_cols ? jb.B + (BLASLONG)off * jb.ldb : jb.B + off;
        pj[p].C = jb.C + off;
      }
      memset(&args[p], 0, sizeof(args[p]));
      args[p].common = &pj[p];
      memset(&queue[p], 0, sizeof(queue[p]));
#ifdef DOUBLE
      queue[p].mode = BLAS_DOUBLE | BLAS_REAL;
#else
      queue[p].mode = BLAS_SINGLE | BLAS_REAL;
#endif
      queue[p].routine = (void *)s2_routine;
      queue[p].args = &args[p];
      queue[p].next = &queue[p + 1];
      used = p + 1;
    }
    queue[used - 1].next = NULL;
    void *buffer = blas_memory_alloc(0); /* the calling thread runs queue[0] on the workspace it passes */
    if (buffer != NULL) {
      queue[0].sa = (char *)buffer + GEMM_OFFSET_A;
      exec_blas(used, queue);
      blas_memory_free(buffer);
      return;
    }
  }
#endif
  s2_run(&jb, NULL);
}

#pragma clang attribute pop
