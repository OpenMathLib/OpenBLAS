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

/* A panels of the SME2 kernels from the upper (LOWER: lower) triangle of a symmetric matrix, as symm_ucopy_2.c. */

#include "common.h"
#include "sme2_gemm_copy.h"

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, BLASLONG posX, BLASLONG posY, FLOAT *b) {
  for (BLASLONG j0 = 0; j0 < n; j0 += S2C_W) {
    const BLASLONG w = n - j0 < S2C_W ? n - j0 : S2C_W, lo = posX + j0;
    const FLOAT *an = a + posY + lo * lda, *at = a + lo + posY * lda; /* (i, c) at an[i + c * lda], at[c + i * lda] */
    const BLASLONG i1 = s2c_clip(lo - posY, 0, m), i2 = s2c_clip(lo + w - posY, i1, m);
#ifndef LOWER
    s2c_cols(i1, w, an, lda, b);
    s2c_rows(m - i2, w, at + i2 * lda, lda, b + i2 * w);
#else
    s2c_rows(i1, w, at, lda, b);
    s2c_cols(m - i2, w, an + i2, lda, b + i2 * w);
#endif
    for (BLASLONG i = i1; i < i2; ++i)
      for (BLASLONG c = 0; c < w; ++c) {
        const BLASLONG r = posY + i, col = lo + c;
#ifndef LOWER
        b[i * w + c] = r < col ? an[i + c * lda] : at[c + i * lda];
#else
        b[i * w + c] = r > col ? an[i + c * lda] : at[c + i * lda];
#endif
      }
    b += m * w;
  }
  return 0;
}
