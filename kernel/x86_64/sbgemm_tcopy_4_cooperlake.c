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

	IFLOAT *boffset0, *boffset1;

	boffset0   = b;

	BLASLONG n8 = n & ~7;
	BLASLONG m4 = m & ~3;
	BLASLONG m2 = m & ~1;

	for (j = 0; j < n8; j += 8) {
		boffset1 = boffset0 + m * 4;
		for (i = 0; i < m4; i +=4) {
			__m128i a0 = _mm_loadu_si128((void *)&a[(i + 0)*lda + j]);
			__m128i a1 = _mm_loadu_si128((void *)&a[(i + 1)*lda + j]);
			__m128i a2 = _mm_loadu_si128((void *)&a[(i + 2)*lda + j]);
			__m128i a3 = _mm_loadu_si128((void *)&a[(i + 3)*lda + j]);
			__m128i a00 = _mm_unpacklo_epi16(a0, a1);
			__m128i a01 = _mm_unpackhi_epi16(a0, a1);
			__m128i a10 = _mm_unpacklo_epi16(a2, a3);
			__m128i a11 = _mm_unpackhi_epi16(a2, a3);
			_mm_storeu_si128((void *)(boffset0 + 0), a00);
			_mm_storeu_si128((void *)(boffset0 + 8), a10);
			_mm_storeu_si128((void *)(boffset1 + 0), a01);
			_mm_storeu_si128((void *)(boffset1 + 8), a11);
			boffset0 += 16;
			boffset1 += 16;
		}
		for (; i < m2; i+= 2) {
			__m128i a0 = _mm_loadu_si128((void *)&a[(i + 0)*lda + j]);
			__m128i a1 = _mm_loadu_si128((void *)&a[(i + 1)*lda + j]);
			__m128i a00 = _mm_unpacklo_epi16(a0, a1);
			__m128i a01 = _mm_unpackhi_epi16(a0, a1);
			_mm_storeu_si128((void *)(boffset0 + 0), a00);
			_mm_storeu_si128((void *)(boffset1 + 0), a01);
			boffset0 += 8;
			boffset1 += 8;
		}
		for (; i < m; i++) {
			__m128d a0 = _mm_loadu_pd((void *)&a[(i + 0)*lda + j]);
			_mm_store_sd((void *)boffset0, a0);
			_mm_store_sd((void *)boffset1, _mm_permute_pd(a0, 0x1));
			boffset0 += 4;
			boffset1 += 4;
		}
		boffset0 = boffset1;
	}
	if (j < n) {
		uint32_t remains = n - j;
		__mmask8 r_mask = (1UL << remains) - 1;
		if (remains > 4) {
			boffset1 = boffset0 + m * 4;
			uint32_t tail1 = remains - 4;
			__mmask8 w_mask1 = (1UL << tail1) - 1;
			for (i = 0; i < m2; i += 2) {
				__m128i a0 = _mm_maskz_loadu_epi16(r_mask, &a[(i + 0)*lda + j]);
				__m128i a1 = _mm_maskz_loadu_epi16(r_mask, &a[(i + 1)*lda + j]);
				__m128i a00 = _mm_unpacklo_epi16(a0, a1);
				__m128i a01 = _mm_unpackhi_epi16(a0, a1);
				_mm_storeu_si128((void *)boffset0, a00);
				_mm_mask_storeu_epi32((void *)boffset1, w_mask1, a01);
				boffset0 += 8;
				boffset1 += 2 * tail1;
			}
			for (; i < m; i++) {
				__m128i a0 = _mm_maskz_loadu_epi16(r_mask, &a[(i + 0)*lda + j]);
				_mm_store_sd((void *)boffset0, (__m128d) a0);
				_mm_mask_storeu_epi16((void *)boffset1, w_mask1, (__m128i) _mm_permute_pd((__m128d) a0, 0x1));
				boffset0 += 4;
				boffset1 += tail1;
			}
		} else {
			for (i = 0; i < m2; i += 2) {
				__m128i a0 = _mm_maskz_loadu_epi16(r_mask, &a[(i + 0)*lda + j]);
				__m128i a1 = _mm_maskz_loadu_epi16(r_mask, &a[(i + 1)*lda + j]);
				__m128i a00 = _mm_unpacklo_epi16(a0, a1);
				_mm_mask_storeu_epi32((void *)boffset0, r_mask, a00);
				boffset0 += 2 * remains;
			}
			for (; i < m; i++) {
				__m128i a0 = _mm_maskz_loadu_epi16(r_mask, &a[(i + 0)*lda + j]);
				_mm_mask_storeu_epi16((void *)boffset0, r_mask, a0);
			}
		}
	}
	return 0;
}
