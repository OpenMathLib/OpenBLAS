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
   a[k + r * lda], so pairs of k are contiguous in the source. Rows beyond rows are zero. */
static inline void s2bc_pairs_l(BLASLONG m, BLASLONG w, BLASLONG rows, const IFLOAT *a, BLASLONG lda, IFLOAT *b) {
  const unsigned pairs = (unsigned)((m + 1) >> 1);
  unsigned *bo = (unsigned *)b;
  for (unsigned p = 0; p < pairs; ++p)
    for (BLASLONG r = 0; r < w; ++r) {
      unsigned v = 0;
      if (r < rows && 2 * p < m) {
        const IFLOAT *s = a + r * lda + 2 * p;
        v = (unsigned)s[0] | (2 * p + 1 < m ? (unsigned)s[1] << 16 : 0);
      }
      bo[(size_t)p * w + r] = v;
    }
}

/* Same, but element (k, r) at a[r + k * lda]: pairs of k one lda apart in the source. */
static inline void s2bc_pairs_c(BLASLONG m, BLASLONG w, BLASLONG rows, const IFLOAT *a, BLASLONG lda, IFLOAT *b) {
  const unsigned pairs = (unsigned)((m + 1) >> 1);
  unsigned *bo = (unsigned *)b;
  for (unsigned p = 0; p < pairs; ++p)
    for (BLASLONG r = 0; r < w; ++r) {
      unsigned v = 0;
      if (r < rows) {
        if (2 * p < m) v = (unsigned)a[r + (BLASLONG)(2 * p) * lda];
        if (2 * p + 1 < m) v |= (unsigned)a[r + (BLASLONG)(2 * p + 1) * lda] << 16;
      }
      bo[(size_t)p * w + r] = v;
    }
}

#endif
