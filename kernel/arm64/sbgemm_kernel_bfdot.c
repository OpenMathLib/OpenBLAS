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
 * NEON BFDOT sbgemm kernel for Apple cores (VORTEX): C += alpha * A * B on the packed panels
 * of the level-3 driver, the bfloat16 counterpart of the NEON sgemm kernels.
 *
 * One BFDOT does a two-wide dot product per lane, so the kernel zips rows l and l+1 of the
 * standard [k][m] and [k][n] panels into the lane layout (i, k & 1): the B pair of a column
 * is duplicated across the lanes, the A pairs cover four rows each, and acc[n][h] holds four
 * consecutive rows of C column n - one quad load/store per (column, row group). A missing
 * second row of an odd k is a zero vector, and blocks narrower than eight rows or columns run
 * in scalar code.
 *
 * Without FEAT_BF16 at run time, or in the compiler, plain loops run on the same layout.
 */

#include "common.h"
#include <arm_neon.h>
#include <string.h>
#if defined(__APPLE__)
#include <sys/sysctl.h>
#endif

#define S2BD_UM GEMM_DEFAULT_UNROLL_M
#define S2BD_UN GEMM_DEFAULT_UNROLL_N

#if S2BD_UM != 8 || S2BD_UN != 8
#error "the BFDOT sbgemm kernel needs SBGEMM_UNROLL_M 8 and SBGEMM_UNROLL_N 8"
#endif

/* Group width at offset j0 of a count of n: full groups, then the largest power of two left.
   With SBGEMM_UNROLL_M = SBGEMM_UNROLL_N the driver routes the A side through the O copies
   as well, so both panels come in this grouping. */
static inline int s2bd_w(BLASLONG n, BLASLONG j0) {
  BLASLONG r = n - j0;
  if (r >= S2BD_UN) return S2BD_UN;
  int w = 1;
  while (2 * w <= r) w *= 2;
  return w;
}

static inline FLOAT s2bd_cvt(IFLOAT h) {
  unsigned int u = (unsigned int)h << 16;
  FLOAT f;
  memcpy(&f, &u, sizeof(f));
  return f;
}

/* FEAT_BF16 at run time; VORTEX is the Apple target, so sysctl decides. */
static inline int s2bd_usable(void) {
  static int ok = -1;
  if (ok < 0) {
#if defined(__APPLE__)
    int v = 0;
    size_t len = sizeof(v);
    ok = sysctlbyname("hw.optional.arm.FEAT_BF16", &v, &len, NULL, 0) == 0 && v;
#else
    ok = 0;
#endif
  }
  return ok;
}

/* Scalar block of wm rows and wn columns, for tails and for the plain-loop fallback. */
static void s2bd_block(BLASLONG k, const IFLOAT *a, BLASLONG wm, const IFLOAT *b, BLASLONG wn, FLOAT *C,
                       BLASLONG ldc, FLOAT alpha) {
  for (BLASLONG j = 0; j < wn; ++j)
    for (BLASLONG i = 0; i < wm; ++i) {
      FLOAT s = 0;
      for (BLASLONG l = 0; l < k; ++l) s += s2bd_cvt(a[l * wm + i]) * s2bd_cvt(b[l * wn + j]);
      C[i + j * ldc] += alpha * s;
    }
}

#if defined(__clang__)
/* The file is built with the NEON flags of its target; the functions here add FEAT_BF16. */
#pragma clang attribute push(__attribute__((target("bf16"))), apply_to = function)

/* Full 8x8 block of C in sixteen quad accumulators, one quad per (column, row group). */
static void s2bd_mac(BLASLONG k, const IFLOAT *a, const IFLOAT *b, FLOAT *C, BLASLONG ldc, FLOAT alpha) {
  float32x4_t acc[8][2];
  for (int n = 0; n < 8; ++n) acc[n][0] = acc[n][1] = vdupq_n_f32(0);
  static const uint16_t zero[8]; /* the missing row of an odd k */
  for (BLASLONG l = 0; l < k; l += 2) {
    const uint16x8_t ak0 = vld1q_u16((const uint16_t *)a + l * 8);
    const uint16x8_t ak1 = vld1q_u16(l + 1 < k ? (const uint16_t *)a + (l + 1) * 8 : zero);
    const uint16x8_t bk0 = vld1q_u16((const uint16_t *)b + l * 8);
    const uint16x8_t bk1 = vld1q_u16(l + 1 < k ? (const uint16_t *)b + (l + 1) * 8 : zero);
    const uint16x8_t az0 = vzip1q_u16(ak0, ak1); /* (m0,k0),(m0,k1),...,(m3,k1) */
    const uint16x8_t az1 = vzip2q_u16(ak0, ak1); /* (m4,k0),...,(m7,k1) */
    const uint32x4_t bp0 = vreinterpretq_u32_u16(vzip1q_u16(bk0, bk1)); /* (n0,k0),(n0,k1),... */
    const uint32x4_t bp1 = vreinterpretq_u32_u16(vzip2q_u16(bk0, bk1));
/* vgetq_lane needs a literal lane, hence one macro expansion per column. */
#define S2BD_PAIR(N, BP, L)                                                                    \
  do {                                                                                         \
    const uint16x8_t bn = vreinterpretq_u16_u32(vdupq_n_u32(vgetq_lane_u32(BP, L)));          \
    acc[N][0] = vbfdotq_f32(acc[N][0], az0, bn);                                               \
    acc[N][1] = vbfdotq_f32(acc[N][1], az1, bn);                                               \
  } while (0);
    S2BD_PAIR(0, bp0, 0) S2BD_PAIR(1, bp0, 1) S2BD_PAIR(2, bp0, 2) S2BD_PAIR(3, bp0, 3)
    S2BD_PAIR(4, bp1, 0) S2BD_PAIR(5, bp1, 1) S2BD_PAIR(6, bp1, 2) S2BD_PAIR(7, bp1, 3)
#undef S2BD_PAIR
  }
  for (int n = 0; n < 8; ++n)
    for (int h = 0; h < 2; ++h) {
      FLOAT *p = C + n * ldc + 4 * h;
      vst1q_f32(p, vfmaq_n_f32(vld1q_f32(p), acc[n][h], alpha));
    }
}

static void s2bd_drive(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, const IFLOAT *sa, const IFLOAT *sb,
                       FLOAT *C, BLASLONG ldc) {
  for (BLASLONG i0 = 0; i0 < m;) {
    const BLASLONG wm = s2bd_w(m, i0);
    const IFLOAT *a = sa + i0 * k;
    for (BLASLONG j0 = 0; j0 < n;) {
      const BLASLONG wn = s2bd_w(n, j0);
      FLOAT *cb = C + i0 + j0 * ldc;
      if (wm == S2BD_UM && wn == S2BD_UN)
        s2bd_mac(k, a, sb + j0 * k, cb, ldc, alpha);
      else
        s2bd_block(k, a, wm, sb + j0 * k, wn, cb, ldc, alpha);
      j0 += wn;
    }
    i0 += wm;
  }
}

#pragma clang attribute pop
#endif

int CNAME(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, IFLOAT *ba, IFLOAT *bb, FLOAT *C, BLASLONG ldc) {
  if (m <= 0 || n <= 0 || k <= 0) return 0;
#if defined(__clang__)
  if (s2bd_usable()) {
    s2bd_drive(m, n, k, alpha, ba, bb, C, ldc);
    return 0;
  }
#endif
  for (BLASLONG i0 = 0; i0 < m;) {
    const BLASLONG wm = s2bd_w(m, i0);
    for (BLASLONG j0 = 0; j0 < n;) {
      const BLASLONG wn = s2bd_w(n, j0);
      s2bd_block(k, ba + i0 * k, wm, bb + j0 * k, wn, C + i0 + j0 * ldc, ldc, alpha);
      j0 += wn;
    }
    i0 += wm;
  }
  return 0;
}
