/***************************************************************************
Copyright (c) 2014, The OpenBLAS Project
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

#include <immintrin.h>

#define HAVE_KERNEL_4x4 1

static void dsymv_kernel_4x4(BLASLONG from, BLASLONG to, FLOAT **a, FLOAT *x, FLOAT *y, FLOAT *temp1, FLOAT *temp2)
{


	__m256d temp2_0, temp2_1, temp2_2, temp2_3; // temp2_0 temp2_1 temp2_2 temp2_3
	__m256d temp1_0, temp1_1, temp1_2, temp1_3;

	temp2_0 = _mm256_setzero_pd();
	temp2_1 = _mm256_setzero_pd();
	temp2_2 = _mm256_setzero_pd();
	temp2_3 = _mm256_setzero_pd();

	temp1_0 = _mm256_broadcastsd_pd(_mm_load_sd(&temp1[0]));
	temp1_1 = _mm256_broadcastsd_pd(_mm_load_sd(&temp1[1]));
	temp1_2 = _mm256_broadcastsd_pd(_mm_load_sd(&temp1[2]));
	temp1_3 = _mm256_broadcastsd_pd(_mm_load_sd(&temp1[3]));

#ifdef __AVX512CD__
	__m512d temp2_05, temp2_15, temp2_25, temp2_35; // temp2_0 temp2_1 temp2_2 temp2_3
	__m512d temp1_05, temp1_15, temp1_25, temp1_35;
	BLASLONG to2;
	int delta;

	temp2_05 = _mm512_setzero_pd();
	temp2_15 = _mm512_setzero_pd();
	temp2_25 = _mm512_setzero_pd();
	temp2_35 = _mm512_setzero_pd();

	temp1_05 = _mm512_broadcastsd_pd(_mm_load_sd(&temp1[0]));
	temp1_15 = _mm512_broadcastsd_pd(_mm_load_sd(&temp1[1]));
	temp1_25 = _mm512_broadcastsd_pd(_mm_load_sd(&temp1[2]));
	temp1_35 = _mm512_broadcastsd_pd(_mm_load_sd(&temp1[3]));

	delta = (to - from) & ~7;
	to2 = from + delta;


	for (; from < to2; from += 8) {
		__m512d _x, _y;
		__m512d a0, a1, a2, a3;

		_y = _mm512_loadu_pd(&y[from]);
		_x = _mm512_loadu_pd(&x[from]);

		a0 = _mm512_loadu_pd(&a[0][from]);
		a1 = _mm512_loadu_pd(&a[1][from]);
		a2 = _mm512_loadu_pd(&a[2][from]);
		a3 = _mm512_loadu_pd(&a[3][from]);

		_y += temp1_05 * a0 + temp1_15 * a1 + temp1_25 * a2 + temp1_35 * a3;

		temp2_05 += _x * a0;
		temp2_15 += _x * a1;
		temp2_25 += _x * a2;
		temp2_35 += _x * a3;

		_mm512_storeu_pd(&y[from], _y);

	};

	temp2_0 = _mm256_add_pd(_mm512_extractf64x4_pd(temp2_05, 0), _mm512_extractf64x4_pd(temp2_05, 1));
	temp2_1 = _mm256_add_pd(_mm512_extractf64x4_pd(temp2_15, 0), _mm512_extractf64x4_pd(temp2_15, 1));
	temp2_2 = _mm256_add_pd(_mm512_extractf64x4_pd(temp2_25, 0), _mm512_extractf64x4_pd(temp2_25, 1));
	temp2_3 = _mm256_add_pd(_mm512_extractf64x4_pd(temp2_35, 0), _mm512_extractf64x4_pd(temp2_35, 1));

#endif

	for (; from != to; from += 4) {
		__m256d _x, _y;
		__m256d a0, a1, a2, a3;

		_y = _mm256_loadu_pd(&y[from]);
		_x = _mm256_loadu_pd(&x[from]);

		a0 = _mm256_loadu_pd(&a[0][from]);
		a1 = _mm256_loadu_pd(&a[1][from]);
		a2 = _mm256_loadu_pd(&a[2][from]);
		a3 = _mm256_loadu_pd(&a[3][from]);

		_y += temp1_0 * a0 + temp1_1 * a1 + temp1_2 * a2 + temp1_3 * a3;

		temp2_0 += _x * a0;
		temp2_1 += _x * a1;
		temp2_2 += _x * a2;
		temp2_3 += _x * a3;

		_mm256_storeu_pd(&y[from], _y);

	};

	__m128d xmm0, xmm1, xmm2, xmm3;


	xmm0 = _mm_add_pd(_mm256_extractf128_pd(temp2_0, 0), _mm256_extractf128_pd(temp2_0, 1));
	xmm1 = _mm_add_pd(_mm256_extractf128_pd(temp2_1, 0), _mm256_extractf128_pd(temp2_1, 1));
	xmm2 = _mm_add_pd(_mm256_extractf128_pd(temp2_2, 0), _mm256_extractf128_pd(temp2_2, 1));
	xmm3 = _mm_add_pd(_mm256_extractf128_pd(temp2_3, 0), _mm256_extractf128_pd(temp2_3, 1));

	xmm0 = _mm_hadd_pd(xmm0, xmm0);
	xmm1 = _mm_hadd_pd(xmm1, xmm1);
	xmm2 = _mm_hadd_pd(xmm2, xmm2);
	xmm3 = _mm_hadd_pd(xmm3, xmm3);


	temp2[0] += xmm0[0];
	temp2[1] += xmm1[0];
	temp2[2] += xmm2[0];
	temp2[3] += xmm3[0];
} 


