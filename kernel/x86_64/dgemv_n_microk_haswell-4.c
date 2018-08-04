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



#define HAVE_KERNEL_4x4 1

#include <immintrin.h>

static void dgemv_kernel_4x4( BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y, FLOAT *alpha)
{

	int i = 0;

	__m256d x0, x1, x2, x3;
	__m256d __alpha;

	x0 = _mm256_broadcastsd_pd(_mm_load_sd(&x[0]));
	x1 = _mm256_broadcastsd_pd(_mm_load_sd(&x[1]));
	x2 = _mm256_broadcastsd_pd(_mm_load_sd(&x[2]));
	x3 = _mm256_broadcastsd_pd(_mm_load_sd(&x[3]));

	__alpha = _mm256_broadcastsd_pd(_mm_load_sd(alpha));


	for (i = 0; i < n; i+= 4) {
		__m256d tempY;
		__m256d sum;

		sum = _mm256_add_pd(
				_mm256_add_pd(
					_mm256_mul_pd(_mm256_loadu_pd(&ap[0][i]), x0),
					_mm256_mul_pd(_mm256_loadu_pd(&ap[1][i]), x1)),
				_mm256_add_pd(
					_mm256_mul_pd(_mm256_loadu_pd(&ap[2][i]), x2),
					_mm256_mul_pd(_mm256_loadu_pd(&ap[3][i]), x3))
			);

		tempY = _mm256_loadu_pd(&y[i]);
		tempY = _mm256_add_pd(tempY, _mm256_mul_pd(sum, __alpha));
		_mm256_storeu_pd(&y[i], tempY);
	}

} 


#define HAVE_KERNEL_4x2

static void dgemv_kernel_4x2( BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y, FLOAT *alpha) __attribute__ ((noinline));

static void dgemv_kernel_4x2( BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y, FLOAT *alpha)
{

	BLASLONG register i = 0;

	__asm__  __volatile__
	(
	"vbroadcastsd    (%2), %%ymm12	 \n\t"	// x0 
	"vbroadcastsd   8(%2), %%ymm13	 \n\t"	// x1 

	"vmovups	(%4,%0,8), %%ymm0	 \n\t"
	"vmovups	(%5,%0,8), %%ymm1	 \n\t"

	"vbroadcastsd    (%6), %%ymm6 	 \n\t"	// alpha 

        "addq		$4 , %0	  	 	      \n\t"
	"subq	        $4 , %1			      \n\t"		
	"jz		2f		      \n\t"

	"1:				 \n\t"

	"vmulpd        %%ymm0 , %%ymm12, %%ymm4      \n\t" 
	"vmulpd        %%ymm1 , %%ymm13, %%ymm5      \n\t" 
	"vmovups	(%4,%0,8), %%ymm0	 \n\t"
	"vmovups	(%5,%0,8), %%ymm1	 \n\t"

	"vmovups	-32(%3,%0,8), %%ymm8	       \n\t"	// 4 * y
	"vaddpd		 %%ymm4 , %%ymm5 , %%ymm4      \n\t"
	"vfmadd231pd     %%ymm6 , %%ymm4 , %%ymm8      \n\t"

	"vmovups         %%ymm8,   -32(%3,%0,8)	      \n\t"	// 4 * y

        "addq		$4 , %0	  	 	      \n\t"
	"subq	        $4 , %1			      \n\t"		
	"jnz		1b		      \n\t"
	

	"2:				 \n\t"

	"vmulpd        %%ymm0 , %%ymm12, %%ymm4      \n\t" 
	"vmulpd        %%ymm1 , %%ymm13, %%ymm5      \n\t" 


	"vmovups	-32(%3,%0,8), %%ymm8	       \n\t"	// 4 * y
	"vaddpd		 %%ymm4 , %%ymm5 , %%ymm4      \n\t"
	"vfmadd231pd     %%ymm6 , %%ymm4 , %%ymm8      \n\t"

	"vmovups  %%ymm8,   -32(%3,%0,8)	      \n\t"	// 4 * y


	"vzeroupper			              \n\t"


	:
          "+r" (i),	// 0	
	  "+r" (n)  	// 1
	:
          "r" (x),      // 2
          "r" (y),      // 3
          "r" (ap[0]),  // 4
          "r" (ap[1]),  // 5
          "r" (alpha)   // 6
	: "cc", 
	  "%xmm0", "%xmm1", 
	  "%xmm4", "%xmm5", 
	  "%xmm6", 
	  "%xmm8", 
	  "%xmm12", "%xmm13",
	  "memory"
	);
}
