/***************************************************************************
Copyright (c) 2021, The OpenBLAS Project
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
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#include <stdio.h>
#include <immintrin.h>
#include "common.h"

int CNAME(BLASLONG m, BLASLONG n, IFLOAT *a, BLASLONG lda, IFLOAT *b){
  BLASLONG i, j;

  IFLOAT *boffset;

  boffset   = b;

  BLASLONG n32 = n & ~31;
  BLASLONG m4 = m & ~3;
  BLASLONG m2 = m & ~1;

  uint32_t permute_table = {
    0, 0x10|0, 1, 0x10|1, 2, 0x10|2, 3, 0x10|3, 4, 0x10|4, 5, 0x10|5, 6, 0x10|6, 7, 0x10, 7,
    8, 0x10|8, 9, 0x10|9, 10, 0x10|10, 11, 0x10|11, 12, 0x10|12, 13, 0x10|13, 14, 0x10|14, 15, 0x10|15,
  };

  __m512i idx_lo = _mm512_loadu_si512(permute_table);
  __m512i idx_hi = _mm512_loadu_si512(permute_table + 16);

  for (j = 0; j < n32; j += 32) {
    for (i = 0; i < m4; i += 4) {
      /* bf16 fma need special memory layout:
       * for memory layout like below:
       *     a00, a01, a02, a03, a04, a05 ....
       *     a10, a11, a12, a13, a14, a15 ....
       * need to copy as:
       *     a00, a10, a01, a11, a02, a12, a03, a13, ...
       */
      __m512i a0 = _mm512_loadu_si512(&a[(i + 0)*lda + j]);
      __m512i a1 = _mm512_loadu_si512(&a[(i + 1)*lda + j]);
      __m512i a2 = _mm512_loadu_si512(&a[(i + 2)*lda + j]);
      __m512i a3 = _mm512_loadu_si512(&a[(i + 3)*lda + j]);

      __m512i a00 = _mm512_unpacklo_epi16(a0, a1);
      __m512i a01 = _mm512_unpackhi_epi16(a0, a1);
      __m512i a10 = _mm512_unpacklo_epi16(a2, a3);
      __m512i a11 = _mm512_unpackhi_epi16(a2, a3);

      a0 = _mm512_permutex2var_epi32(a00, idx_lo, a01);
      a1 = _mm512_permutex2var_epi32(a00, idx_hi, a01);
      a2 = _mm512_permutex2var_epi32(a10, idx_lo, a11);
      a3 = _mm512_permutex2var_epi32(a10, idx_hi, a11);

      _mm512_storeu_si512(boffset, a0);
      _mm512_storeu_si512(boffset + 32, a1);
      _mm512_storeu_si512(boffset + 64, a2);
      _mm512_storeu_si512(boffset + 96, a3);
      boffset += 128;
    }
    for (; i < m2; i += 2) {
      __m512i a0 = _mm512_loadu_si512(&a[(i + 0)*lda + j]);
      __m512i a1 = _mm512_loadu_si512(&a[(i + 1)*lda + j]);

      __m512i a00 = _mm512_unpacklo_epi16(a0, a1);
      __m512i a01 = _mm512_unpackhi_epi16(a0, a1);

      a0 = _mm512_permutex2var_epi32(a00, idx_lo, a01);
      a1 = _mm512_permutex2var_epi32(a00, idx_hi, a01);

      _mm512_storeu_si512(boffset, a0);
      _mm512_storeu_si512(boffset + 32, a1);
      boffset += 64;
    }
    for (; i < m; i++) {
      /* just copy the only remains row */
      __m512i a0 = _mm512_loadu_si512(&a[(i + 0)*lda + j]);
      _mm512_storeu_si512(boffset, a0);
      boffset += 32;
    }
  }
  if (j < n) {
    uint32_t remains = n - j;
    __mmask32 r_mask = (1UL << remains) - 1;
    if (remains > 16) {
      __mmask16 w_mask = (1UL << (remains - 16)) - 1;
      for (i = 0; i < m2; i += 2) {
        __m512i a0 = _mm512_maskz_loadu_epi16(r_mask, &a[(i + 0)*lda + j]);
        __m512i a1 = _mm512_maskz_loadu_epi16(r_mask, &a[(i + 1)*lda + j]);

        __m512i a00 = _mm512_unpacklo_epi16(a0, a1);
        __m512i a01 = _mm512_unpackhi_epi16(a0, a1);

        a0 = _mm512_permutex2var_epi32(a00, idx_lo, a01);
        a1 = _mm512_permutex2var_epi32(a00, idx_hi, a01);

        _mm512_storeu_si512(boffset, a0);
        _mm512_mask_storeu_epi32(boffset + 32, w_mask, a1);
        boffset += 2 * remains;
      }
    } else {
      __mmask16 w_mask = (1UL << remains ) - 1;
      for (i = 0; i < m2; i += 2) {
        __m512i a0 = _mm512_maskz_loadu_epi16(r_mask, &a[(i + 0)*lda + j]);
        __m512i a1 = _mm512_maskz_loadu_epi16(r_mask, &a[(i + 1)*lda + j]);

        __m512i a00 = _mm512_unpacklo_epi16(a0, a1);
        __m512i a01 = _mm512_unpackhi_epi16(a0, a1);

        a0 = _mm512_permutex2var_epi32(a00, idx_lo, a01);
        _mm512_mask_storeu_epi32(boffset, w_mask, a0);
        boffset += 2 * remains;
      }
    }
    for (; i < m; i++) {
        __m512i a0 = _mm512_maskz_loadu_epi16(r_mask, &a[(i + 0)*lda + j]);
        _mm512_mask_storeu_epi16(boffset, r_mask, a0);
        boffset += remains;
    }
  }
}
