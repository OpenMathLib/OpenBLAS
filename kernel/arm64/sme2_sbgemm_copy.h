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

/* bfloat16 packing of the SME2 and BFDOT sbgemm kernels: k-pair interleaved panels.
   BFMOPA/BFDOT consume two k steps per instruction, so the panels hold element
   (row r, k) of a block of w rows at [pp * 2w + 2r + (k & 1)], one k pair per
   load. The odd k half of the last pair is zeroed by the copies; the A-side
   copies additionally zero-pad the ragged chunk to the full chunk width. The
   B-side copies pack their power-of-two groups at the actual width: consumers
   must not read beyond 2*w lanes of a w-column group (the SME2 kernel
   predicates its B load to them, the BFDOT kernel runs its vector path on full
   groups only), and C stores are predicated to the valid rows and columns. */

#ifndef SME2_SBGEMM_COPY_H
#define SME2_SBGEMM_COPY_H

#include <arm_neon.h>
#include <string.h>

#define S2BC_W GEMM_DEFAULT_UNROLL_M
#define S2BC_WN GEMM_DEFAULT_UNROLL_N

/* Width of the group at offset j0 of a count of n: full groups, then the largest power of two left. */
static inline BLASLONG s2bc_wn(BLASLONG n, BLASLONG j0) {
  BLASLONG r = n - j0;
  if (r >= S2BC_WN) return S2BC_WN;
  BLASLONG w = 1;
  while (2 * w <= r) w *= 2;
  return w;
}

/* Block of m k-steps by w rows into b (w u32 per k pair); element (k, r) at
   a[k + r * lda], so pairs of k are contiguous in the source. Rows beyond rows are zero.
   The main loop is a 4x4 u32 transpose: each quad load takes four pairs of one source
   row, each store writes one pair of four rows, so the pairs land in the panel order. */
static inline void s2bc_pairs_l(BLASLONG m, BLASLONG w, BLASLONG rows, const IFLOAT *a, BLASLONG lda, IFLOAT *b) {
  const unsigned pairs = (unsigned)((m + 1) >> 1), full = (unsigned)(m >> 1);
  unsigned *bo = (unsigned *)b;
  memset(bo, 0, (size_t)pairs * w * sizeof(unsigned));
  BLASLONG r = 0;
  for (; r + 4 <= rows; r += 4) {
    const uint16_t *s0 = (const uint16_t *)a + r * lda, *s1 = s0 + lda, *s2 = s1 + lda, *s3 = s2 + lda;
    unsigned p = 0;
    for (; p + 4 <= full; p += 4) {
      const uint32x4_t v0 = vreinterpretq_u32_u16(vld1q_u16(s0 + 2 * p));
      const uint32x4_t v1 = vreinterpretq_u32_u16(vld1q_u16(s1 + 2 * p));
      const uint32x4_t v2 = vreinterpretq_u32_u16(vld1q_u16(s2 + 2 * p));
      const uint32x4_t v3 = vreinterpretq_u32_u16(vld1q_u16(s3 + 2 * p));
      const uint32x4_t t0 = vtrn1q_u32(v0, v1), t1 = vtrn2q_u32(v0, v1);
      const uint32x4_t t2 = vtrn1q_u32(v2, v3), t3 = vtrn2q_u32(v2, v3);
      vst1q_u32(bo + (size_t)(p + 0) * w + r, vcombine_u32(vget_low_u32(t0), vget_low_u32(t2)));
      vst1q_u32(bo + (size_t)(p + 1) * w + r, vcombine_u32(vget_low_u32(t1), vget_low_u32(t3)));
      vst1q_u32(bo + (size_t)(p + 2) * w + r, vcombine_u32(vget_high_u32(t0), vget_high_u32(t2)));
      vst1q_u32(bo + (size_t)(p + 3) * w + r, vcombine_u32(vget_high_u32(t1), vget_high_u32(t3)));
    }
    for (; p < full; ++p)
      for (BLASLONG q = 0; q < 4; ++q) {
        const IFLOAT *s = a + (r + q) * lda + 2 * p;
        bo[(size_t)p * w + r + q] = (unsigned)s[0] | ((unsigned)s[1] << 16);
      }
  }
  for (; r < rows; ++r)
    for (unsigned p = 0; p < full; ++p) {
      const IFLOAT *s = a + r * lda + 2 * p;
      bo[(size_t)p * w + r] = (unsigned)s[0] | ((unsigned)s[1] << 16);
    }
  if (m & 1)
    for (BLASLONG q = 0; q < rows; ++q) bo[(size_t)(pairs - 1) * w + q] = (unsigned)a[q * lda + m - 1];
}

/* Same, but element (k, r) at a[r + k * lda]: along the rows both halves of a k pair
   form contiguous streams, so the main loop zips one quad of each into the panel order. */
static inline void s2bc_pairs_c(BLASLONG m, BLASLONG w, BLASLONG rows, const IFLOAT *a, BLASLONG lda, IFLOAT *b) {
  const unsigned pairs = (unsigned)((m + 1) >> 1), full = (unsigned)(m >> 1);
  unsigned *bo = (unsigned *)b;
  memset(bo, 0, (size_t)pairs * w * sizeof(unsigned));
  for (unsigned p = 0; p < full; ++p) {
    const uint16_t *se = (const uint16_t *)a + (BLASLONG)(2 * p) * lda;
    const uint16_t *so = se + lda;
    BLASLONG r = 0;
    for (; r + 8 <= rows; r += 8) {
      const uint16x8_t e = vld1q_u16(se + r), o = vld1q_u16(so + r);
      vst1q_u32(bo + (size_t)p * w + r, vreinterpretq_u32_u16(vzip1q_u16(e, o)));
      vst1q_u32(bo + (size_t)p * w + r + 4, vreinterpretq_u32_u16(vzip2q_u16(e, o)));
    }
    for (; r < rows; ++r) bo[(size_t)p * w + r] = (unsigned)se[r] | ((unsigned)so[r] << 16);
  }
  if (m & 1) {
    const uint16_t *se = (const uint16_t *)a + (BLASLONG)(m - 1) * lda;
    for (BLASLONG r = 0; r < rows; ++r) bo[(size_t)(pairs - 1) * w + r] = (unsigned)se[r];
  }
}

#endif
