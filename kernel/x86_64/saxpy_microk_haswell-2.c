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

static void saxpy_kernel_16( BLASLONG n, FLOAT *x, FLOAT *y, FLOAT *alpha)
{
	BLASLONG i = 0;

	__m256 __alpha;

	__alpha =  _mm256_broadcastss_ps(_mm_load_ss(alpha));

	for (; i < n; i+= 32) {
		__m256 y0, y8, y16, y24;

		y0  = _mm256_loadu_ps(&y[i +  0]);
		y8  = _mm256_loadu_ps(&y[i +  8]);
		y16 = _mm256_loadu_ps(&y[i + 16]);
		y24 = _mm256_loadu_ps(&y[i + 24]);

		y0  += __alpha * _mm256_loadu_ps(&x[i +  0]);
		y8  += __alpha * _mm256_loadu_ps(&x[i +  8]);
		y16 += __alpha * _mm256_loadu_ps(&x[i + 16]);
		y24 += __alpha * _mm256_loadu_ps(&x[i + 24]);

		_mm256_storeu_ps(&y[i +  0], y0);
		_mm256_storeu_ps(&y[i +  8], y8);
		_mm256_storeu_ps(&y[i + 16], y16);
		_mm256_storeu_ps(&y[i + 24], y24);
	}
}
#endif

