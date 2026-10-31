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
/* bfloat16 A panels of the SME2 sbgemm kernel (sbgemm_kernel_sme2.c), the same
   [k][m] layout in chunks of GEMM_UNROLL_M rows as gemm_ncopy_sme2.c. */

#ifndef SME2_SBGEMM_COPY_H
#define SME2_SBGEMM_COPY_H

#include <string.h>

#define S2BC_W GEMM_DEFAULT_UNROLL_M
#define S2BC_WN GEMM_DEFAULT_UNROLL_N

/* Element (i, c) at a[c + i * lda], bfloat16: a plain transpose for now. */
static inline void s2bc_cols(BLASLONG m, BLASLONG w, const IFLOAT *a, BLASLONG lda, IFLOAT *b) {
  for (BLASLONG i = 0; i < m; ++i)
    for (BLASLONG c = 0; c < w; ++c) b[i * w + c] = a[i + c * lda];
}

/* Element (i, c) at a[i * lda + c], bfloat16: contiguous rows. */
static inline void s2bc_rows(BLASLONG m, BLASLONG w, const IFLOAT *a, BLASLONG lda, IFLOAT *b) {
  for (BLASLONG i = 0; i < m; ++i) memcpy(b + i * w, a + i * lda, (size_t)w * sizeof(IFLOAT));
}

/* Width of the group at column j0 of n for the B copies: full groups of S2BC_N, then powers of two. */
static inline BLASLONG s2bc_wn(BLASLONG n, BLASLONG j0) {
  BLASLONG r = n - j0;
  if (r >= S2BC_WN) return S2BC_WN;
  BLASLONG w = 1;
  while (2 * w <= r) w *= 2;
  return w;
}

#endif
