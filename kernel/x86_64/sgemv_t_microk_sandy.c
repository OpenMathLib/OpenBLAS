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

static void  sgemv_kernel_16( long n, float alpha, float *a, long lda, float *x, float *y)
{

	//n = n / 16;

	__asm__ __volatile__
	(
	"movq 	        %0, 	 %%rax\n\t"		// n -> rax
	"vmovss         %1, 	 %%xmm1\n\t"		// alpha -> xmm1
	"movq		%2,	 %%rsi\n\t"		// adress of a -> rsi
	"movq	        %3,	 %%rcx\n\t"		// value of lda > rcx
	"movq		%4,	 %%rdi\n\t"		// adress of x -> rdi
	"movq		%5,	 %%rdx\n\t"		// adress of y -> rdx

	"leaq	(, %%rcx,4), %%rcx       \n\t"		// scale lda by size of float
	"leaq	(%%rsi,%%rcx,1), %%r8    \n\t"		// pointer to next line

	"vxorps		%%xmm12, %%xmm12, %%xmm12\n\t"	// set to zero
	"vxorps		%%xmm13, %%xmm13, %%xmm13\n\t"	// set to zero
	"vxorps		%%xmm14, %%xmm14, %%xmm14\n\t"	// set to zero
	"vxorps		%%xmm15, %%xmm15, %%xmm15\n\t"	// set to zero

	"sarq		$4, %%rax		 \n\t"	// n = n / 16

	".align 16				 \n\t"
	".L01LOOP%=:				 \n\t"
	// "prefetcht0	512(%%rsi)		 \n\t"
	"prefetcht0	(%%r8)		 	 \n\t" //prefetch next line of a
	"vmovups	(%%rsi), %%xmm4		 \n\t"
	"vmovups     4*4(%%rsi), %%xmm5		 \n\t"
	"vmovups     8*4(%%rsi), %%xmm6		 \n\t"
	"vmovups    12*4(%%rsi), %%xmm7		 \n\t"

	"vmulps      0*4(%%rdi), %%xmm4, %%xmm8 \n\t" // multiply a and c and add to temp
	"vmulps      4*4(%%rdi), %%xmm5, %%xmm9 \n\t" // multiply a and c and add to temp
	"vmulps      8*4(%%rdi), %%xmm6, %%xmm10\n\t" // multiply a and c and add to temp
	"vmulps     12*4(%%rdi), %%xmm7, %%xmm11\n\t" // multiply a and c and add to temp

	"vaddps		%%xmm12, %%xmm8 , %%xmm12\n\t"	
	"vaddps		%%xmm13, %%xmm9 , %%xmm13\n\t"	
	"vaddps		%%xmm14, %%xmm10, %%xmm14\n\t"	
	"vaddps		%%xmm15, %%xmm11, %%xmm15\n\t"	

        "addq		$16*4    ,   %%r8	 \n\t"  // increment prefetch pointer 
        "addq		$16*4    ,   %%rsi	 \n\t"  // increment pointer of a 
        "addq		$16*4    ,   %%rdi	 \n\t"  // increment pointer of c 
	"dec		%%rax			 \n\t"  // n = n -1
	"jnz		.L01LOOP%=		 \n\t"

	"vaddps		%%xmm12, %%xmm14, %%xmm12\n\t"	
	"vaddps		%%xmm13, %%xmm15, %%xmm13\n\t"	
	"vaddps		%%xmm12, %%xmm13, %%xmm12\n\t"	
	"vhaddps	%%xmm12, %%xmm12, %%xmm12\n\t"	
	"vhaddps	%%xmm12, %%xmm12, %%xmm12\n\t"	

	"vmulss		%%xmm12, %%xmm1, %%xmm12 \n\t"
	"vaddss	       (%%rdx), %%xmm12, %%xmm12\n\t"
	"vmovss		%%xmm12, (%%rdx)	 \n\t"  // store temp -> y

	:
        :
          "m" (n),	// 0	
	  "m" (alpha),  // 1
	  "m" (a),      // 2
          "m" (lda),    // 3
          "m" (x),      // 4
          "m" (y)       // 5
	: "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "cc",
	  "%xmm0", "%xmm1", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7",
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11",
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	);

} 



