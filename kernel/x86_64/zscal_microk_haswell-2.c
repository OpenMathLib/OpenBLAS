/***************************************************************************
Copyright (c) 2014-2015, The OpenBLAS Project
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


#define HAVE_KERNEL 1
#define HAVE_KERNEL_8 1

static void zscal_kernel( BLASLONG n, FLOAT *alpha, FLOAT *x) __attribute__ ((noinline));

static void zscal_kernel( BLASLONG n, FLOAT *alpha, FLOAT *x)
{


	__asm__  __volatile__
	(
	"vbroadcastsd		(%2), %%ymm0		    \n\t"  // da_r	
	"vbroadcastsd          8(%2), %%ymm1		    \n\t"  // da_i 	

	"cmpq	        $8 , %0			            \n\t"
	"jb	3f					    \n\t"

	"addq	$128, %1				    \n\t"

	"vmovups	-128(%1), %%ymm4		    \n\t"
	"vmovups	 -96(%1), %%ymm5		    \n\t"
	"vmovups	 -64(%1), %%ymm6		    \n\t"
	"vmovups	 -32(%1), %%ymm7		    \n\t"

	"vpermilpd	$0x05 , %%ymm4, %%ymm12		    \n\t"
	"vpermilpd	$0x05 , %%ymm5, %%ymm13		    \n\t"
	"vpermilpd	$0x05 , %%ymm6, %%ymm14 	    \n\t"
	"vpermilpd	$0x05 , %%ymm7, %%ymm15		    \n\t"

	"subq	        $8 , %0			            \n\t"		
	"cmpq	        $8 , %0			            \n\t"
	"jb	2f					    \n\t"

	".p2align 4				            \n\t"
	"1:				            	    \n\t"

	//"prefetcht0     128(%1)				    \n\t"
	// ".align 2				            \n\t"

	"vmulpd		%%ymm1, %%ymm12, %%ymm8		    \n\t" // da_i*x1 , da_i *x0
	"vmulpd		%%ymm1, %%ymm13, %%ymm9		    \n\t"
	"vmulpd		%%ymm1, %%ymm14, %%ymm10	    \n\t"
	"vmulpd		%%ymm1, %%ymm15, %%ymm11	    \n\t"

	"vfmaddsub231pd	%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	"vmovups	   0(%1), %%ymm4		    \n\t"
	"vfmaddsub231pd	%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	"vmovups	  32(%1), %%ymm5		    \n\t"
	"vfmaddsub231pd	%%ymm0, %%ymm6 , %%ymm10	    \n\t"
	"vmovups	  64(%1), %%ymm6		    \n\t"
	"vfmaddsub231pd	%%ymm0, %%ymm7 , %%ymm11	    \n\t"
	"vmovups	  96(%1), %%ymm7		    \n\t"

	"vmovups	%%ymm8 , -128(%1)		    \n\t"
	"vpermilpd	$0x05 , %%ymm4, %%ymm12		    \n\t"
	"vmovups	%%ymm9 ,  -96(%1)		    \n\t"
	"vpermilpd	$0x05 , %%ymm5, %%ymm13		    \n\t"
	"vmovups	%%ymm10,  -64(%1)		    \n\t"
	"vpermilpd	$0x05 , %%ymm6, %%ymm14 	    \n\t"
	"vmovups	%%ymm11,  -32(%1)		    \n\t"
	"vpermilpd	$0x05 , %%ymm7, %%ymm15		    \n\t"

	"addq		$128 ,%1  	 	            \n\t"
	"subq	        $8 , %0			            \n\t"		
	"cmpq	        $8 , %0			            \n\t"
	"jae		1b		             	    \n\t"

	"2:				            	    \n\t"


	"vmulpd		%%ymm1, %%ymm12, %%ymm8		    \n\t" // da_i*x1 , da_i *x0
	"vmulpd		%%ymm1, %%ymm13, %%ymm9		    \n\t"
	"vmulpd		%%ymm1, %%ymm14, %%ymm10	    \n\t"
	"vmulpd		%%ymm1, %%ymm15, %%ymm11	    \n\t"

	"vfmaddsub231pd	%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	"vfmaddsub231pd	%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	"vfmaddsub231pd	%%ymm0, %%ymm6 , %%ymm10	    \n\t"
	"vfmaddsub231pd	%%ymm0, %%ymm7 , %%ymm11	    \n\t"

	"vmovups	%%ymm8 , -128(%1)		    \n\t"
	"vmovups	%%ymm9 ,  -96(%1)		    \n\t"
	"vmovups	%%ymm10,  -64(%1)		    \n\t"
	"vmovups	%%ymm11,  -32(%1)		    \n\t"

	"testq		$7, %0			            \n\t"
	"jz		6f			            \n\t"

	"3:				            	    \n\t"
	"testq		$4, %0			            \n\t"
	"jz		4f			            \n\t"
	"vmovups	   0(%1), %%ymm4		    \n\t"
	"vmovups	  32(%1), %%ymm5		    \n\t"
	"vpermilpd	$0x05 , %%ymm4 , %%ymm12	    \n\t"
	"vpermilpd	$0x05 , %%ymm5 , %%ymm13	    \n\t"
	"vmulpd		%%ymm1, %%ymm12, %%ymm8		    \n\t" // da_i*x1 , da_i *x0
	"vmulpd		%%ymm1, %%ymm13, %%ymm9		    \n\t"
	"vfmaddsub231pd	%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	"vfmaddsub231pd	%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	"vmovups	%%ymm8,    0(%1)		    \n\t"
	"vmovups	%%ymm9,   32(%1)		    \n\t"
	"addq		$64, %1  	 	            \n\t"

	"4:				            	    \n\t"
	"testq		$2, %0			            \n\t"
	"jz		5f			            \n\t"
	"vmovups	   0(%1), %%ymm4		    \n\t"
	"vpermilpd	$0x05 , %%ymm4 , %%ymm12	    \n\t"
	"vmulpd		%%ymm1, %%ymm12, %%ymm8		    \n\t" // da_i*x1 , da_i *x0
	"vfmaddsub231pd	%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	"vmovups	%%ymm8,    0(%1)		    \n\t"
	"addq		$32, %1  	 	            \n\t"

	"5:				            	    \n\t"
	"testq		$1, %0			            \n\t"
	"jz		6f			            \n\t"
	"vmovups	   0(%1), %%xmm4		    \n\t"
	"vpermilpd	$0x01 , %%xmm4 , %%xmm12	    \n\t"
	"vmulpd		%%xmm1, %%xmm12, %%xmm8		    \n\t" // da_i*x1 , da_i *x0
	"vfmaddsub231pd	%%xmm0, %%xmm4 , %%xmm8		    \n\t" // da_r*x0 , da_r *x1
	"vmovups	%%xmm8,    0(%1)		    \n\t"

	"6:				            	    \n\t"
	"vzeroupper					    \n\t"

	:
	  "+r" (n),  	// 0
          "+r" (x)      // 1
	:
          "r" (alpha)   // 2
	: "cc",
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	);

} 

static void zscal_kernel_8_zero_r( BLASLONG n, FLOAT *alpha, FLOAT *x) __attribute__ ((noinline));

static void zscal_kernel_8_zero_r( BLASLONG n, FLOAT *alpha, FLOAT *x)
{


	__asm__  __volatile__
	(
	"vxorpd	           %%ymm0, %%ymm0, %%ymm0	    \n\t"
	"vbroadcastsd          8(%2), %%ymm1		    \n\t"  // da_i 	

	"addq	$128, %1				    \n\t"

	"vmovups	-128(%1), %%ymm4		    \n\t"
	"vmovups	 -96(%1), %%ymm5		    \n\t"
	"vmovups	 -64(%1), %%ymm6		    \n\t"
	"vmovups	 -32(%1), %%ymm7		    \n\t"

	"vpermilpd	$0x05 , %%ymm4, %%ymm12		    \n\t"
	"vpermilpd	$0x05 , %%ymm5, %%ymm13		    \n\t"
	"vpermilpd	$0x05 , %%ymm6, %%ymm14 	    \n\t"
	"vpermilpd	$0x05 , %%ymm7, %%ymm15		    \n\t"

	"subq	        $8 , %0			            \n\t"		
	"jz	2f					    \n\t"

	".p2align 4				            \n\t"
	"1:				            	    \n\t"

	//"prefetcht0     128(%1)				    \n\t"
	// ".align 2				            \n\t"

	"vmovups	   0(%1), %%ymm4		    \n\t"
	"vmovups	  32(%1), %%ymm5		    \n\t"
	"vmovups	  64(%1), %%ymm6		    \n\t"
	"vmovups	  96(%1), %%ymm7		    \n\t"

	"vmulpd		%%ymm1, %%ymm12, %%ymm12	    \n\t" // da_i*x1 , da_i *x0
	"vaddsubpd	%%ymm12 , %%ymm0 , %%ymm8	    \n\t"
	"vmulpd		%%ymm1, %%ymm13, %%ymm13	    \n\t" 
	"vaddsubpd	%%ymm13 , %%ymm0 , %%ymm9	    \n\t"
	"vmulpd		%%ymm1, %%ymm14, %%ymm14	    \n\t" 
	"vaddsubpd	%%ymm14 , %%ymm0 , %%ymm10	    \n\t"
	"vmulpd		%%ymm1, %%ymm15, %%ymm15	    \n\t" 
	"vaddsubpd	%%ymm15 , %%ymm0 , %%ymm11	    \n\t"

	"vmovups	%%ymm8 , -128(%1)		    \n\t"
	"vpermilpd	$0x05 , %%ymm4, %%ymm12		    \n\t"
	"vmovups	%%ymm9 ,  -96(%1)		    \n\t"
	"vpermilpd	$0x05 , %%ymm5, %%ymm13		    \n\t"
	"vmovups	%%ymm10,  -64(%1)		    \n\t"
	"vpermilpd	$0x05 , %%ymm6, %%ymm14 	    \n\t"
	"vmovups	%%ymm11,  -32(%1)		    \n\t"
	"vpermilpd	$0x05 , %%ymm7, %%ymm15		    \n\t"

	"addq		$128 ,%1  	 	            \n\t"
	"subq	        $8 , %0			            \n\t"		
	"jnz		1b		             	    \n\t"

	"2:				            	    \n\t"

	"vmulpd		%%ymm1, %%ymm12, %%ymm12	    \n\t" // da_i*x1 , da_i *x0
	"vaddsubpd	%%ymm12 , %%ymm0 , %%ymm8	    \n\t"
	"vmulpd		%%ymm1, %%ymm13, %%ymm13	    \n\t" 
	"vaddsubpd	%%ymm13 , %%ymm0 , %%ymm9	    \n\t"
	"vmulpd		%%ymm1, %%ymm14, %%ymm14	    \n\t" 
	"vaddsubpd	%%ymm14 , %%ymm0 , %%ymm10	    \n\t"
	"vmulpd		%%ymm1, %%ymm15, %%ymm15	    \n\t" 
	"vaddsubpd	%%ymm15 , %%ymm0 , %%ymm11	    \n\t"

	"vmovups	%%ymm8 , -128(%1)		    \n\t"
	"vmovups	%%ymm9 ,  -96(%1)		    \n\t"
	"vmovups	%%ymm10,  -64(%1)		    \n\t"
	"vmovups	%%ymm11,  -32(%1)		    \n\t"

	"vzeroupper					    \n\t"

	:
	  "+r" (n),  	// 0
          "+r" (x)      // 1
	:
          "r" (alpha)   // 2
	: "cc",
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	);

} 



static void zscal_kernel_8_zero_i( BLASLONG n, FLOAT *alpha, FLOAT *x) __attribute__ ((noinline));

static void zscal_kernel_8_zero_i( BLASLONG n, FLOAT *alpha, FLOAT *x)
{


	__asm__  __volatile__
	(
	"vbroadcastsd		(%2), %%ymm0		    \n\t"  // da_r	

	"addq	$128, %1				    \n\t"

	"vmovups	-128(%1), %%ymm4		    \n\t"
	"vmovups	 -96(%1), %%ymm5		    \n\t"
	"vmovups	 -64(%1), %%ymm6		    \n\t"
	"vmovups	 -32(%1), %%ymm7		    \n\t"


	"subq	        $8 , %0			            \n\t"		
	"jz	2f					    \n\t"

	".p2align 4				            \n\t"
	"1:				            	    \n\t"

	//"prefetcht0     128(%1)				    \n\t"
	// ".align 2				            \n\t"

	"vmulpd		%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	"vmovups	   0(%1), %%ymm4		    \n\t"
	"vmulpd		%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	"vmovups	  32(%1), %%ymm5		    \n\t"
	"vmulpd		%%ymm0, %%ymm6 , %%ymm10	    \n\t" 
	"vmovups	  64(%1), %%ymm6		    \n\t"
	"vmulpd		%%ymm0, %%ymm7 , %%ymm11	    \n\t" 
	"vmovups	  96(%1), %%ymm7		    \n\t"

	"vmovups	%%ymm8 , -128(%1)		    \n\t"
	"vmovups	%%ymm9 ,  -96(%1)		    \n\t"
	"vmovups	%%ymm10,  -64(%1)		    \n\t"
	"vmovups	%%ymm11,  -32(%1)		    \n\t"

	"addq		$128 ,%1  	 	            \n\t"
	"subq	        $8 , %0			            \n\t"		
	"jnz		1b		             	    \n\t"

	"2:				            	    \n\t"


	"vmulpd		%%ymm0, %%ymm4 , %%ymm8		    \n\t" // da_r*x0 , da_r *x1
	"vmulpd		%%ymm0, %%ymm5 , %%ymm9		    \n\t"
	"vmulpd		%%ymm0, %%ymm6 , %%ymm10	    \n\t" 
	"vmulpd		%%ymm0, %%ymm7 , %%ymm11	    \n\t" 

	"vmovups	%%ymm8 , -128(%1)		    \n\t"
	"vmovups	%%ymm9 ,  -96(%1)		    \n\t"
	"vmovups	%%ymm10,  -64(%1)		    \n\t"
	"vmovups	%%ymm11,  -32(%1)		    \n\t"

	"vzeroupper					    \n\t"

	:
	  "+r" (n),  	// 0
          "+r" (x)      // 1
	:
          "r" (alpha)   // 2
	: "cc",
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	);

} 


static void zscal_kernel_8_zero( BLASLONG n, FLOAT *alpha, FLOAT *x) __attribute__ ((noinline));

static void zscal_kernel_8_zero( BLASLONG n, FLOAT *alpha, FLOAT *x)
{


	__asm__  __volatile__
	(
	"vxorpd	           %%ymm0, %%ymm0, %%ymm0	    \n\t"

	"addq	$128, %1				    \n\t"

	".p2align 4				            \n\t"
	"1:				            	    \n\t"

	//"prefetcht0     128(%1)				    \n\t"
	// ".align 2				            \n\t"

	"vmovups	%%ymm0 , -128(%1)		    \n\t"
	"vmovups	%%ymm0 ,  -96(%1)		    \n\t"
	"vmovups	%%ymm0 ,  -64(%1)		    \n\t"
	"vmovups	%%ymm0 ,  -32(%1)		    \n\t"

	"addq		$128 ,%1  	 	            \n\t"
	"subq	        $8 , %0			            \n\t"		
	"jnz		1b		             	    \n\t"

	"vzeroupper					    \n\t"

	:
	  "+r" (n),  	// 0
          "+r" (x)      // 1
	:
          "r" (alpha)   // 2
	: "cc",
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	);

} 



