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

/* Row interchanges of LASWP on n columns, the rows k1..k2 copied into panels of 16 columns, then 8, 4, 2, 1. */

#include "common.h"

int CNAME(BLASLONG n, BLASLONG k1, BLASLONG k2, FLOAT *a, BLASLONG lda, blasint *ipiv, FLOAT *buffer) {
  if (n <= 0 || k1 > k2) return 0;
  a--;
  ipiv += k1 - 1;
  for (BLASLONG j0 = 0; j0 < n;) {
    BLASLONG w = 16;
    while (w > n - j0) w /= 2;
    for (BLASLONG i = k1; i <= k2; ++i) {
      const BLASLONG ip = ipiv[i - k1];
      for (BLASLONG c = 0; c < w; ++c) {
        FLOAT *col = a + (j0 + c) * lda;
        const FLOAT x = col[i], y = col[ip];
        if (ip != i) {
          col[ip] = x;
          buffer[c] = y;
        } else {
          buffer[c] = x;
        }
      }
      buffer += w;
    }
    j0 += w;
  }
  return 0;
}
