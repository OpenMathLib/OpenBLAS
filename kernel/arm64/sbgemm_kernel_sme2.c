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

/* A panels of the SME2 kernels (sme2_gemm_packed.h) from op(A) with contiguous columns. */
/*
 * SME2 bf16 GEMM kernel for a 512-bit streaming vector length, the bfloat16 counterpart of
 * gemm_kernel_sme2.c: C += alpha * A * B on the packed panels of the level-3 driver.
 *
 * BFMOPA ZA32 with .H operands accumulates 2-wide dot products, so one instruction covers
 * two k steps. The sbgemm copies pack the panels k-pair interleaved (element (row r, k) of
 * a block of w rows at [pp * 2w + 2r + (k & 1)]) with the padding zeroed, and the kernel
 * needs one plain load per operand tile and k pair, no zipping. C accumulates in four ZA32
 * tiles of 16x16 fp32: tile rows are the wn columns of the block, tile columns the wm rows,
 * with the C load/store helpers of sme2_gemm_tile.h.
 *
 * A panels come in chunks of SBGEMM_UNROLL_M rows padded to the full width, B panels in
 * power-of-two groups of SBGEMM_UNROLL_N columns. Without SME2 at run time, or in the
 * compiler, plain loops run on the same layout.
 */

#include "common.h"
#include <string.h>

#include "sme2_gemm_detect.h"

#define S2B_UM GEMM_DEFAULT_UNROLL_M
#define S2B_UN GEMM_DEFAULT_UNROLL_N

#if S2B_UM != 64 || S2B_UN != 16
#error "the SME2 sbgemm kernel needs SBGEMM_UNROLL_M 64 and SBGEMM_UNROLL_N 16"
#endif

/* Width of the B panel group at column j0 of n: full groups, then the largest power of two left. */
static inline int s2b_wn(BLASLONG n, BLASLONG j0) {
  BLASLONG r = n - j0;
  if (r >= S2B_UN) return S2B_UN;
  int w = 1;
  while (2 * w <= r) w *= 2;
  return w;
}

static inline FLOAT s2b_cvt(IFLOAT h) {
  unsigned int u = (unsigned int)h << 16;
  FLOAT f;
  memcpy(&f, &u, sizeof(f));
  return f;
}

#ifdef HAVE_SME2_GEMM
#include "sme2_gemm_tile.h"

#if S2_VL != 16 || S2_NR != 64 || S2_NT != 4
#error "the SME2 sbgemm kernel needs a 512-bit streaming vector length"
#endif

/* The file is built with the SME flags of its target (no SME2); the functions here enable SME2 themselves. */
#pragma clang attribute push(__attribute__((target("sme2"))), apply_to = function)

S2_INL void s2b_op(const int SUB, int t, svbool_t p, svbfloat16_t b, svbfloat16_t a) S2_S S2_ZA {
  switch (t) {
    case 0: SUB ? svmops_za32_bf16_m(0, p, p, b, a) : svmopa_za32_bf16_m(0, p, p, b, a); break;
    case 1: SUB ? svmops_za32_bf16_m(1, p, p, b, a) : svmopa_za32_bf16_m(1, p, p, b, a); break;
    case 2: SUB ? svmops_za32_bf16_m(2, p, p, b, a) : svmopa_za32_bf16_m(2, p, p, b, a); break;
    default: SUB ? svmops_za32_bf16_m(3, p, p, b, a) : svmopa_za32_bf16_m(3, p, p, b, a); break;
  }
}

/* ZA += or -= the A panel times the B panel (wn columns) over the k-pair interleaved panels.
   The A loads read the zero padding of the full-width chunk; the B load is predicated to the
   group's actual 2*wn lanes, so lanes of a narrower group read as zero. */
S2_NOINL void s2b_mac(const int SUB, BLASLONG k, const IFLOAT *a, const IFLOAT *b, BLASLONG wn) S2_S S2_ZA {
  const svbool_t pt = svptrue_b16(); /* BFMOPA governs per 16-bit element */
  const svbool_t pb = svwhilelt_b16((int64_t)0, (int64_t)2 * wn);
  const BLASLONG pairs = (k + 1) >> 1;
  for (BLASLONG l = 0; l < pairs; ++l) {
    const svbfloat16_t bv = svreinterpret_bf16_u16(svld1_u16(pb, (const uint16_t *)b + l * 2 * wn));
    s2b_op(SUB, 0, pt, bv, svreinterpret_bf16_u16(svld1_u16(pt, (const uint16_t *)a + l * 2 * S2B_UM)));
    s2b_op(SUB, 1, pt, bv, svreinterpret_bf16_u16(svld1_u16(pt, (const uint16_t *)a + l * 2 * S2B_UM + 32)));
    s2b_op(SUB, 2, pt, bv, svreinterpret_bf16_u16(svld1_u16(pt, (const uint16_t *)a + l * 2 * S2B_UM + 64)));
    s2b_op(SUB, 3, pt, bv, svreinterpret_bf16_u16(svld1_u16(pt, (const uint16_t *)a + l * 2 * S2B_UM + 96)));
  }
}

/* Block of wm rows and wn columns of C; sgn 1 for C += alpha AB with alpha = +-1 (the AB part
   added or subtracted in ZA), sgn 0 for a general alpha scaled on the way out of ZA. */
S2_NOINL void s2b_blk(BLASLONG k, const IFLOAT *a, BLASLONG wm, const IFLOAT *b, BLASLONG wn, FLOAT *C,
                             BLASLONG ldc, FLOAT alpha, int sgn) S2_S S2_ZA {
  if (sgn) {
    s2_c_load(1, S2_NT, C, ldc, (int)wn, (int)wm, S2_LOAD, (FLOAT)1);
    s2b_mac(alpha == (FLOAT)-1, k, a, b, wn);
    s2_c_store(1, S2_NT, C, ldc, (int)wn, (int)wm);
  } else {
    svzero_za();
    s2b_mac(0, k, a, b, wn);
    s2_c_store_alpha(S2_NT, C, ldc, (int)wn, (int)wm, alpha, 1);
  }
}

__arm_new("za") __arm_locally_streaming static void s2b_drive(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha,
                                                              const IFLOAT *sa, const IFLOAT *sb, FLOAT *C,
                                                              BLASLONG ldc) {
  const BLASLONG pairs = (k + 1) >> 1;
  const int sgn = alpha == (FLOAT)1 || alpha == (FLOAT)-1;
  for (BLASLONG i0 = 0; i0 < m; i0 += S2B_UM) {
    const BLASLONG wm = m - i0 < S2B_UM ? m - i0 : S2B_UM;
    const IFLOAT *a = sa + i0 * pairs * 2;
    for (BLASLONG j0 = 0; j0 < n;) {
      const BLASLONG wn = n - j0 >= S2B_UN ? S2B_UN : s2b_wn(n, j0);
      s2b_blk(k, a, wm, sb + j0 * pairs * 2, wn, C + i0 + j0 * ldc, ldc, alpha, sgn);
      j0 += wn;
    }
  }
}

#pragma clang attribute pop
#endif

static void s2b_loops(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, const IFLOAT *sa, const IFLOAT *sb, FLOAT *C,
                      BLASLONG ldc) {
  const BLASLONG pairs = (k + 1) >> 1;
  for (BLASLONG i0 = 0; i0 < m; i0 += S2B_UM) {
    const BLASLONG wm = m - i0 < S2B_UM ? m - i0 : S2B_UM;
    const IFLOAT *a = sa + i0 * pairs * 2;
    for (BLASLONG j0 = 0; j0 < n;) {
      const BLASLONG wn = n - j0 >= S2B_UN ? S2B_UN : s2b_wn(n, j0);
      const IFLOAT *b = sb + j0 * pairs * 2;
      for (BLASLONG j = 0; j < wn; ++j)
        for (BLASLONG i = 0; i < wm; ++i) {
          FLOAT s = 0;
          for (BLASLONG p = 0; p < pairs; ++p)
            s += s2b_cvt(a[p * 2 * S2B_UM + 2 * i]) * s2b_cvt(b[p * 2 * wn + 2 * j]) +
                 s2b_cvt(a[p * 2 * S2B_UM + 2 * i + 1]) * s2b_cvt(b[p * 2 * wn + 2 * j + 1]);
          C[i0 + i + (j0 + j) * ldc] += alpha * s;
        }
      j0 += wn;
    }
  }
}

int CNAME(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, IFLOAT *ba, IFLOAT *bb, FLOAT *C, BLASLONG ldc) {
  if (m <= 0 || n <= 0 || k <= 0) return 0;
#ifdef HAVE_SME2_GEMM
  if (s2_usable()) {
    s2b_drive(m, n, k, alpha, ba, bb, C, ldc);
    return 0;
  }
#endif
  s2b_loops(m, n, k, alpha, ba, bb, C, ldc);
  return 0;
}
