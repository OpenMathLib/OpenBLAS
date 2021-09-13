/***************************************************************************
 * Copyright (c) 2021, The OpenBLAS Project
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in
 * the documentation and/or other materials provided with the
 * distribution.
 * 3. Neither the name of the OpenBLAS project nor the names of
 * its contributors may be used to endorse or promote products
 * derived from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * *****************************************************************************/

#include <immintrin.h>
#include "common.h"

#define COPY_32(N) _mm512_storeu_si512(boffset + 32 * N, _mm512_loadu_si512(aoffset##N + i))
#define MASK_COPY_32(N) _mm512_mask_storeu_epi16(boffset + tail_m * N, mmask, _mm512_maskz_loadu_epi16(mmask, aoffset##N + i))
#define COPY_ODD_TAIL(N) *(boffset + N) = *(aoffset##N + i);


int CNAME(BLASLONG m, BLASLONG n, IFLOAT *a, BLASLONG lda, IFLOAT *b) {
	BLASLONG i, j;
	IFLOAT *aoffset, *boffset;
	IFLOAT *aoffset0, *aoffset1, *aoffset2, *aoffset3;
	IFLOAT *aoffset4, *aoffset5, *aoffset6, *aoffset7;
	IFLOAT *aoffset8, *aoffset9, *aoffset10, *aoffset11;
	IFLOAT *aoffset12, *aoffset13, *aoffset14, *aoffset15;

	aoffset = a;
	boffset = b;

	BLASLONG n16 = n & ~15;
	BLASLONG m32 = m & ~31;
	BLASLONG m2 = m & ~1;

	for (j = 0; j < n16; j += 16) {
		aoffset0  = aoffset;
		aoffset1  = aoffset0 + lda;
		aoffset2  = aoffset1 + lda;
		aoffset3  = aoffset2 + lda;
		aoffset4  = aoffset3 + lda;
		aoffset5  = aoffset4 + lda;
		aoffset6  = aoffset5 + lda;
		aoffset7  = aoffset6 + lda;
		aoffset8  = aoffset7 + lda;
		aoffset9  = aoffset8 + lda;
		aoffset10 = aoffset9 + lda;
		aoffset11 = aoffset10 + lda;
		aoffset12 = aoffset11 + lda;
		aoffset13 = aoffset12 + lda;
		aoffset14 = aoffset13 + lda;
		aoffset15 = aoffset14 + lda;
		aoffset += 16 * lda;
		for (i = 0; i < m32; i += 32) {
			COPY_32(0); COPY_32(1); COPY_32(2); COPY_32(3);
			COPY_32(4); COPY_32(5); COPY_32(6); COPY_32(7);
			COPY_32(8); COPY_32(9); COPY_32(10); COPY_32(11);
			COPY_32(12); COPY_32(13); COPY_32(14); COPY_32(15);
			boffset += 32 * 16;
		}
		if (i < m2) {
			int tail_m = m2 - i;
			__mmask32 mmask = (1UL << tail_m) - 1;
			MASK_COPY_32(0); MASK_COPY_32(1); MASK_COPY_32(2); MASK_COPY_32(3);
			MASK_COPY_32(4); MASK_COPY_32(5); MASK_COPY_32(6); MASK_COPY_32(7);
			MASK_COPY_32(8); MASK_COPY_32(9); MASK_COPY_32(10); MASK_COPY_32(11);
			MASK_COPY_32(12); MASK_COPY_32(13); MASK_COPY_32(14); MASK_COPY_32(15);
			i = m2;
			boffset += tail_m * 16;
		}
		if (i < m) {
			/* the tail odd k should put alone */
			COPY_ODD_TAIL(0); COPY_ODD_TAIL(1); COPY_ODD_TAIL(2); COPY_ODD_TAIL(3);
			COPY_ODD_TAIL(4); COPY_ODD_TAIL(5); COPY_ODD_TAIL(6); COPY_ODD_TAIL(7);
			COPY_ODD_TAIL(8); COPY_ODD_TAIL(9); COPY_ODD_TAIL(10); COPY_ODD_TAIL(11);
			COPY_ODD_TAIL(12); COPY_ODD_TAIL(13); COPY_ODD_TAIL(14); COPY_ODD_TAIL(15);
			boffset += 16;
		}
	}
	if (j < n) {
		int remain_n = n - j;
		aoffset0  = aoffset;
		aoffset1  = aoffset0 + lda;
		aoffset2  = aoffset1 + lda;
		aoffset3  = aoffset2 + lda;
		aoffset4  = aoffset3 + lda;
		aoffset5  = aoffset4 + lda;
		aoffset6  = aoffset5 + lda;
		aoffset7  = aoffset6 + lda;
		aoffset8  = aoffset7 + lda;
		aoffset9  = aoffset8 + lda;
		aoffset10 = aoffset9 + lda;
		aoffset11 = aoffset10 + lda;
		aoffset12 = aoffset11 + lda;
		aoffset13 = aoffset12 + lda;
		aoffset14 = aoffset13 + lda;
		aoffset15 = aoffset14 + lda;
		for (i = 0; i < m32; i += 32) {
			switch(remain_n) {
				case 15: COPY_32(14);
				case 14: COPY_32(13);
				case 13: COPY_32(12);
				case 12: COPY_32(11);
				case 11: COPY_32(10);
				case 10: COPY_32(9);
				case  9: COPY_32(8);
				case  8: COPY_32(7);
				case  7: COPY_32(6);
				case  6: COPY_32(5);
				case  5: COPY_32(4);
				case  4: COPY_32(3);
				case  3: COPY_32(2);
				case  2: COPY_32(1);
				case  1: COPY_32(0);
			}
			boffset += 32 * remain_n;
		}
		if (i < m2) {
			int tail_m = m2 - i;
			__mmask32 mmask = (1UL << tail_m) - 1;
			switch(remain_n) {
				case 15: MASK_COPY_32(14);
				case 14: MASK_COPY_32(13);
				case 13: MASK_COPY_32(12);
				case 12: MASK_COPY_32(11);
				case 11: MASK_COPY_32(10);
				case 10: MASK_COPY_32(9);
				case  9: MASK_COPY_32(8);
				case  8: MASK_COPY_32(7);
				case  7: MASK_COPY_32(6);
				case  6: MASK_COPY_32(5);
				case  5: MASK_COPY_32(4);
				case  4: MASK_COPY_32(3);
				case  3: MASK_COPY_32(2);
				case  2: MASK_COPY_32(1);
				case  1: MASK_COPY_32(0);
			}
			i = m2;
			boffset += tail_m * remain_n;
		}
		if (i < m) {
			switch(remain_n) {
				case 15: COPY_ODD_TAIL(14);
				case 14: COPY_ODD_TAIL(13);
				case 13: COPY_ODD_TAIL(12);
				case 12: COPY_ODD_TAIL(11);
				case 11: COPY_ODD_TAIL(10);
				case 10: COPY_ODD_TAIL(9);
				case  9: COPY_ODD_TAIL(8);
				case  8: COPY_ODD_TAIL(7);
				case  7: COPY_ODD_TAIL(6);
				case  6: COPY_ODD_TAIL(5);
				case  5: COPY_ODD_TAIL(4);
				case  4: COPY_ODD_TAIL(3);
				case  3: COPY_ODD_TAIL(2);
				case  2: COPY_ODD_TAIL(1);
				case  1: COPY_ODD_TAIL(0);
			}
		}
	}
	return 0;
}
