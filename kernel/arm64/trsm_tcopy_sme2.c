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

/* A panels of the SME2 TRSM kernels from a transposed upper (LOWER: lower) triangular matrix, as trsm_utcopy_2.c. */

#include "common.h"
#include "sme2_gemm_copy.h"

#ifndef LOWER
#define KEEP(r, c) ((r) > (c))
#else
#define KEEP(r, c) ((r) < (c))
#endif
#ifndef UNIT
#define INV(x) (ONE / (x))
#else
#define INV(x) (ONE)
#endif

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, BLASLONG offset, FLOAT *b) {
  for (BLASLONG j0 = 0; j0 < n; j0 += S2C_W) {
    const BLASLONG w = n - j0 < S2C_W ? n - j0 : S2C_W, lo = offset + j0;
    const FLOAT *ap = a + j0;
    const BLASLONG i1 = s2c_clip(lo, 0, m), i2 = s2c_clip(lo + w, i1, m);
#ifndef LOWER
    s2c_zero(i1, w, b);
    s2c_rows(m - i2, w, ap + i2 * lda, lda, b + i2 * w);
#else
    s2c_rows(i1, w, ap, lda, b);
    s2c_zero(m - i2, w, b + i2 * w);
#endif
    for (BLASLONG i = i1; i < i2; ++i)
      for (BLASLONG c = 0; c < w; ++c) {
        const BLASLONG r = i, col = lo + c;
        b[i * w + c] = r == col ? INV(ap[c + i * lda]) : (KEEP(r, col) ? ap[c + i * lda] : ZERO);
      }
    b += m * w;
  }
  return 0;
}
