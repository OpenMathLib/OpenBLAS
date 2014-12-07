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

#define HAVE_KERNEL_8 1
static void caxpy_kernel_8( BLASLONG n, FLOAT *x, FLOAT *y , FLOAT *alpha) __attribute__ ((noinline));

static void caxpy_kernel_8( BLASLONG n, FLOAT *x, FLOAT *y, FLOAT *alpha)
{


	BLASLONG register i = 0;

	__asm__  __volatile__
	(
	"vbroadcastss		(%4), %%xmm0		    \n\t"  // real part of alpha
	"vbroadcastss	       4(%4), %%xmm1		    \n\t"  // imag part of alpha

	".align 16				            \n\t"
	"1:				            \n\t"

	"prefetcht0  768(%2,%0,4)                           \n\t"
	"vmovups        (%2,%0,4), %%xmm5                   \n\t" // 2 complex values from x
	"vmovups      16(%2,%0,4), %%xmm7                   \n\t" // 2 complex values from x
	"vmovups      32(%2,%0,4), %%xmm9                   \n\t" // 2 complex values from x
	"vmovups      48(%2,%0,4), %%xmm11                  \n\t" // 2 complex values from x
	"prefetcht0  768(%3,%0,4)                           \n\t"

#if !defined(CONJ)
	"vfmaddps    (%3,%0,4), %%xmm0 , %%xmm5, %%xmm12    \n\t"
	"vpermilps	$0xb1 , %%xmm5 , %%xmm4 	    \n\t"  // exchange real and imag part
	"vmulps         %%xmm1, %%xmm4 , %%xmm4             \n\t"

	"vfmaddps  16(%3,%0,4), %%xmm0 , %%xmm7, %%xmm13    \n\t"
	"vpermilps	$0xb1 , %%xmm7 , %%xmm6 	    \n\t"  // exchange real and imag part
	"vmulps         %%xmm1, %%xmm6 , %%xmm6             \n\t"

	"vfmaddps  32(%3,%0,4), %%xmm0 , %%xmm9, %%xmm14    \n\t"
	"vpermilps	$0xb1 , %%xmm9 , %%xmm8 	    \n\t"  // exchange real and imag part
	"vmulps         %%xmm1, %%xmm8 , %%xmm8             \n\t"

	"vfmaddps  48(%3,%0,4), %%xmm0 , %%xmm11,%%xmm15    \n\t"
	"vpermilps	$0xb1 , %%xmm11, %%xmm10 	    \n\t"  // exchange real and imag part
	"vmulps         %%xmm1, %%xmm10, %%xmm10            \n\t"

        "vaddsubps      %%xmm4, %%xmm12, %%xmm12            \n\t"
        "vaddsubps      %%xmm6, %%xmm13, %%xmm13            \n\t"
        "vaddsubps      %%xmm8, %%xmm14, %%xmm14            \n\t"
        "vaddsubps      %%xmm10,%%xmm15, %%xmm15            \n\t"

#else

	"vmulps		%%xmm0, %%xmm5, %%xmm4		    \n\t" // a_r*x_r, a_r*x_i
	"vmulps		%%xmm1, %%xmm5, %%xmm5 		    \n\t" // a_i*x_r, a_i*x_i
	"vmulps		%%xmm0, %%xmm7, %%xmm6		    \n\t" // a_r*x_r, a_r*x_i
	"vmulps		%%xmm1, %%xmm7, %%xmm7 		    \n\t" // a_i*x_r, a_i*x_i
	"vmulps		%%xmm0, %%xmm9, %%xmm8		    \n\t" // a_r*x_r, a_r*x_i
	"vmulps		%%xmm1, %%xmm9, %%xmm9 		    \n\t" // a_i*x_r, a_i*x_i
	"vmulps		%%xmm0, %%xmm11, %%xmm10	    \n\t" // a_r*x_r, a_r*x_i
	"vmulps		%%xmm1, %%xmm11, %%xmm11	    \n\t" // a_i*x_r, a_i*x_i

	"vpermilps	$0xb1 , %%xmm4 , %%xmm4 	    \n\t"  // exchange real and imag part
	"vaddsubps	%%xmm4 ,%%xmm5 , %%xmm4 	    \n\t"
	"vpermilps	$0xb1 , %%xmm4 , %%xmm4 	    \n\t"  // exchange real and imag part

	"vpermilps	$0xb1 , %%xmm6 , %%xmm6 	    \n\t"  // exchange real and imag part
	"vaddsubps	%%xmm6 ,%%xmm7 , %%xmm6 	    \n\t"
	"vpermilps	$0xb1 , %%xmm6 , %%xmm6 	    \n\t"  // exchange real and imag part

	"vpermilps	$0xb1 , %%xmm8 , %%xmm8 	    \n\t"  // exchange real and imag part
	"vaddsubps	%%xmm8 ,%%xmm9 , %%xmm8 	    \n\t"
	"vpermilps	$0xb1 , %%xmm8 , %%xmm8 	    \n\t"  // exchange real and imag part

	"vpermilps	$0xb1 , %%xmm10, %%xmm10	    \n\t"  // exchange real and imag part
	"vaddsubps	%%xmm10,%%xmm11, %%xmm10	    \n\t"
	"vpermilps	$0xb1 , %%xmm10, %%xmm10	    \n\t"  // exchange real and imag part

	"vaddps	     (%3,%0,4) ,%%xmm4 , %%xmm12	    \n\t"
	"vaddps	   16(%3,%0,4) ,%%xmm6 , %%xmm13	    \n\t"
	"vaddps	   32(%3,%0,4) ,%%xmm8 , %%xmm14	    \n\t"
	"vaddps	   48(%3,%0,4) ,%%xmm10, %%xmm15	    \n\t"


#endif

	"vmovups	%%xmm12,   (%3,%0,4)		    \n\t"
	"vmovups	%%xmm13, 16(%3,%0,4)		    \n\t"
	"vmovups	%%xmm14, 32(%3,%0,4)		    \n\t"
	"vmovups	%%xmm15, 48(%3,%0,4)		    \n\t"

	"addq		$16, %0	  	 	             \n\t"
	"subq	        $8 , %1			             \n\t"		
	"jnz		1b		             \n\t"

	:
        : 
          "r" (i),	// 0	
	  "r" (n),  	// 1
          "r" (x),      // 2
          "r" (y),      // 3
          "r" (alpha)   // 4
	: "cc", 
	  "%xmm0", "%xmm1",
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	);

} 


