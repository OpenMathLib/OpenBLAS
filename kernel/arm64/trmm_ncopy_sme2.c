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

/* A panels of the SME2 kernels from an upper (LOWER: lower) triangular matrix, as trmm_uncopy_2.c. */

#include "common.h"
#include "sme2_gemm_copy.h"

#ifndef LOWER
#define ZL 1
#else
#define ZL 0
#endif
#ifdef UNIT
#define UNITD 1
#else
#define UNITD 0
#endif

#if ZL
#define PART(COPY, ZERO)   \
  COPY(i2, w, ap, lda, b); \
  ZERO(m - i2, w, b + i2 * w)
#else
#define PART(COPY, ZERO) \
  ZERO(i1, w, b);        \
  COPY(m - i1, w, ap + i1, lda, b + i1 * w)
#endif
/* Panels of w columns: the rows on the kept side of the diagonal copied, the others zeroed, then the diagonal block. */
#define PANELS(COPY, ZERO, TRI)                                                                         \
  for (BLASLONG j0 = 0, w; j0 < n; j0 += w, b += m * w) {                                               \
    w = n - j0 < S2C_W ? n - j0 : S2C_W;                                                                \
    const BLASLONG lo = posY + j0, i1 = s2c_clip(lo - posX, 0, m), i2 = s2c_clip(lo + w - posX, i1, m); \
    const FLOAT *ap = a + posX + lo * lda;                                                              \
    PART(COPY, ZERO);                                                                                   \
    TRI(i1, i2, w, posX - lo, ZL, UNITD, b);                                                            \
  }

#ifdef HAVE_SME2_GEMM
__arm_new("za") __arm_locally_streaming static
    __attribute__((noinline)) void panels_sme(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, BLASLONG posX,
                                              BLASLONG posY, FLOAT *b) {
  PANELS(s2c_cols_s, s2c_zero_s, s2c_tri_s)
}
#endif

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, BLASLONG posX, BLASLONG posY, FLOAT *b) {
#ifdef HAVE_SME2_GEMM
  if (m * n >= S2C_SME && s2_usable()) {
    panels_sme(m, n, a, lda, posX, posY, b);
    return 0;
  }
#endif
  PANELS(s2c_cols, s2c_zero, s2c_tri)
  return 0;
}
