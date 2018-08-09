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

#ifndef __AVX512CD__
#pragma GCC target("avx2,fma")
#endif

#ifdef __AVX2__

#define HAVE_KERNEL_16 1

#include <immintrin.h>

static void sdot_kernel_16( BLASLONG n, FLOAT *x, FLOAT *y, FLOAT *dot)

{
	int i = 0;
	__m256 accum_0, accum_1, accum_2, accum_3;

	accum_0 = _mm256_setzero_ps();
	accum_1 = _mm256_setzero_ps();
	accum_2 = _mm256_setzero_ps();
	accum_3 = _mm256_setzero_ps();

#ifdef __AVX512CD__
	__m512 accum_05, accum_15, accum_25, accum_35;
	int n64;
	n64 = n & (~63);

	accum_05 = _mm512_setzero_ps();
	accum_15 = _mm512_setzero_ps();
	accum_25 = _mm512_setzero_ps();
	accum_35 = _mm512_setzero_ps();

	for (; i < n64; i += 64) {
		accum_05 += _mm512_loadu_ps(&x[i+ 0]) * _mm512_loadu_ps(&y[i+ 0]);
		accum_15 += _mm512_loadu_ps(&x[i+16]) * _mm512_loadu_ps(&y[i+16]);
		accum_25 += _mm512_loadu_ps(&x[i+32]) * _mm512_loadu_ps(&y[i+32]);
		accum_35 += _mm512_loadu_ps(&x[i+48]) * _mm512_loadu_ps(&y[i+48]);
	}

	/*
	 * we need to fold our 512 bit wide accumulator vectors into 256 bit wide vectors so that the AVX2 code
	 * below can continue using the intermediate results in its loop
	 */
	accum_0 = _mm256_add_ps(_mm512_extractf32x8_ps(accum_05, 0), _mm512_extractf32x8_ps(accum_05, 1));
	accum_1 = _mm256_add_ps(_mm512_extractf32x8_ps(accum_15, 0), _mm512_extractf32x8_ps(accum_15, 1));
	accum_2 = _mm256_add_ps(_mm512_extractf32x8_ps(accum_25, 0), _mm512_extractf32x8_ps(accum_25, 1));
	accum_3 = _mm256_add_ps(_mm512_extractf32x8_ps(accum_35, 0), _mm512_extractf32x8_ps(accum_35, 1));

#endif
	for (; i < n; i += 32) {
		accum_0 += _mm256_loadu_ps(&x[i+ 0]) * _mm256_loadu_ps(&y[i+ 0]);
		accum_1 += _mm256_loadu_ps(&x[i+ 8]) * _mm256_loadu_ps(&y[i+ 8]);
		accum_2 += _mm256_loadu_ps(&x[i+16]) * _mm256_loadu_ps(&y[i+16]);
		accum_3 += _mm256_loadu_ps(&x[i+24]) * _mm256_loadu_ps(&y[i+24]);
	}

	/* we now have the partial sums of the dot product in the 4 accumulation vectors, time to consolidate */

	accum_0 = accum_0 + accum_1 + accum_2 + accum_3;

	__m128 half_accum0;

	/* Add upper half to lower half of each of the 256 bit vector to get a 128 bit vector */
	half_accum0 = _mm_add_ps(_mm256_extractf128_ps(accum_0, 0), _mm256_extractf128_ps(accum_0, 1));

	/* in 128 bit land there is a hadd operation to do the rest of the element-wise sum in one go */
	half_accum0 = _mm_hadd_ps(half_accum0, half_accum0);
	half_accum0 = _mm_hadd_ps(half_accum0, half_accum0);

	*dot = half_accum0[0];
}

#endif
