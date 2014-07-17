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


static void sgemv_kernel_32( long n, float alpha, float *a, long lda, float *x, float *y)
{


	float *pre = a + lda*4*2;

	__asm __volatile
	(
	"movq 	        %0, 	 %%rax\n\t"		// n -> rax
	"vbroadcastss   %1, 	 %%ymm1\n\t"		// alpha -> ymm1
	"movq		%2,	 %%rsi\n\t"		// adress of a -> rsi
	"movq	        %3,	 %%rcx\n\t"		// value of lda > rcx
	"movq		%4,	 %%rdi\n\t"		// adress of x -> rdi
	"movq		%5,	 %%rdx\n\t"		// adress of y -> rdx
	"movq		%6,	 %%r8\n\t"		// address for prefetch
	"prefetcht0	(%%r8)\n\t"			// Prefetch
	"prefetcht0   64(%%r8)\n\t"			// Prefetch

	"vxorps		%%ymm12, %%ymm12, %%ymm12\n\t"	// set to zero
	"vxorps		%%ymm13, %%ymm13, %%ymm13\n\t"	// set to zero
	"vxorps		%%ymm14, %%ymm14, %%ymm14\n\t"	// set to zero
	"vxorps		%%ymm15, %%ymm15, %%ymm15\n\t"	// set to zero

	".L01LOOP%=:				 \n\t"
	"vbroadcastss	(%%rdi),   %%ymm0	 \n\t"	// load values of c
        "addq		$4     ,   %%rdi	 \n\t"  // increment pointer of c 

	"leaq     (%%r8 , %%rcx, 4), %%r8	 \n\t"	// add lda to pointer for prefetch
	"prefetcht0	(%%r8)\n\t"			// Prefetch
	"prefetcht0   64(%%r8)\n\t"			// Prefetch

	"vfmaddps %%ymm12,   0*4(%%rsi), %%ymm0, %%ymm12\n\t" // multiply a and c and add to temp
	"vfmaddps %%ymm13,   8*4(%%rsi), %%ymm0, %%ymm13\n\t" // multiply a and c and add to temp
	"vfmaddps %%ymm14,  16*4(%%rsi), %%ymm0, %%ymm14\n\t" // multiply a and c and add to temp
	"vfmaddps %%ymm15,  24*4(%%rsi), %%ymm0, %%ymm15\n\t" // multiply a and c and add to temp

	"leaq     (%%rsi, %%rcx, 4), %%rsi	 \n\t"	// add lda to pointer of a

	"dec		%%rax			 \n\t"  // n = n -1
	"jnz		.L01LOOP%=		 \n\t"

	"vmulps		%%ymm12, %%ymm1,  %%ymm12\n\t"  // scale by alpha
	"vmulps		%%ymm13, %%ymm1,  %%ymm13\n\t"  // scale by alpha
	"vmulps		%%ymm14, %%ymm1,  %%ymm14\n\t"  // scale by alpha
	"vmulps		%%ymm15, %%ymm1,  %%ymm15\n\t"  // scale by alpha

	"vmovups	%%ymm12,     (%%rdx)	 \n\t"  // store temp -> y
	"vmovups	%%ymm13,  8*4(%%rdx)	 \n\t"  // store temp -> y
	"vmovups	%%ymm14, 16*4(%%rdx)	 \n\t"  // store temp -> y
	"vmovups	%%ymm15, 24*4(%%rdx)	 \n\t"  // store temp -> y

	:
        :
          "m" (n),	// 0	
	  "m" (alpha),  // 1
	  "m" (a),      // 2
          "m" (lda),    // 3
          "m" (x),      // 4
          "m" (y),      // 5
	  "m" (pre)	// 6
	: "rax", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11",
	  "xmm0" , "xmm1", 
	  "xmm12", "xmm13", "xmm14", "xmm15",
	  "memory"
	);

} 

static void sgemv_kernel_16( long n, float alpha, float *a, long lda, float *x, float *y)
{

	float *pre = a + lda*4*3;

	__asm __volatile
	(
	"movq 	        %0, 	 %%rax\n\t"		// n -> rax
	"vbroadcastss   %1, 	 %%ymm1\n\t"		// alpha -> ymm1
	"movq		%2,	 %%rsi\n\t"		// adress of a -> rsi
	"movq	        %3,	 %%rcx\n\t"		// value of lda > rcx
	"movq		%4,	 %%rdi\n\t"		// adress of x -> rdi
	"movq		%5,	 %%rdx\n\t"		// adress of y -> rdx
	"movq		%6,	 %%r8\n\t"		// address for prefetch
	"prefetcht0	(%%r8)\n\t"			// Prefetch

	"vxorps		%%ymm12, %%ymm12, %%ymm12\n\t"	// set to zero
	"vxorps		%%ymm13, %%ymm13, %%ymm13\n\t"	// set to zero

	".L01LOOP%=:				 \n\t"
	"vbroadcastss	(%%rdi),   %%ymm0	 \n\t"	// load values of c
        "addq		$4     ,   %%rdi	 \n\t"  // increment pointer of c 

	"leaq     (%%r8 , %%rcx, 4), %%r8	 \n\t"	// add lda to pointer for prefetch
	"prefetcht0	(%%r8)\n\t"			// Prefetch

	"vfmaddps %%ymm12,   0*4(%%rsi), %%ymm0, %%ymm12\n\t" // multiply a and c and add to temp
	"vfmaddps %%ymm13,   8*4(%%rsi), %%ymm0, %%ymm13\n\t" // multiply a and c and add to temp

	"leaq     (%%rsi, %%rcx, 4), %%rsi	 \n\t"	// add lda to pointer of a

	"dec		%%rax			 \n\t"  // n = n -1
	"jnz		.L01LOOP%=		 \n\t"

	"vmulps		%%ymm12, %%ymm1,  %%ymm12\n\t"  // scale by alpha
	"vmulps		%%ymm13, %%ymm1,  %%ymm13\n\t"  // scale by alpha

	"vmovups	%%ymm12,     (%%rdx)	 \n\t"  // store temp -> y
	"vmovups	%%ymm13,  8*4(%%rdx)	 \n\t"  // store temp -> y

	:
        :
          "m" (n),	// 0	
	  "m" (alpha),  // 1
	  "m" (a),      // 2
          "m" (lda),    // 3
          "m" (x),      // 4
          "m" (y),      // 5
	  "m" (pre)	// 6
	: "rax", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11",
	  "xmm0" , "xmm1", 
	  "xmm12", "xmm13", "xmm14", "xmm15",
	  "memory"
	);

} 


static void sgemv_kernel_8( long n, float alpha, float *a, long lda, float *x, float *y)
{


	__asm __volatile
	(
	"movq 	        %0, 	 %%rax\n\t"		// n -> rax
	"vbroadcastss   %1, 	 %%ymm1\n\t"		// alpha -> ymm1
	"movq		%2,	 %%rsi\n\t"		// adress of a -> rsi
	"movq	        %3,	 %%rcx\n\t"		// value of lda > rcx
	"movq		%4,	 %%rdi\n\t"		// adress of x -> rdi
	"movq		%5,	 %%rdx\n\t"		// adress of y -> rdx

	"vxorps		%%ymm12, %%ymm12, %%ymm12\n\t"	// set to zero

	".L01LOOP%=:				 \n\t"
	"vbroadcastss	(%%rdi),   %%ymm0	 \n\t"	// load values of c
        "addq		$4     ,   %%rdi	 \n\t"  // increment pointer of c 

	"vfmaddps %%ymm12,   0*4(%%rsi), %%ymm0, %%ymm12\n\t" // multiply a and c and add to temp

	"leaq     (%%rsi, %%rcx, 4), %%rsi	 \n\t"	// add lda to pointer of a

	"dec		%%rax			 \n\t"  // n = n -1
	"jnz		.L01LOOP%=		 \n\t"

	"vmulps		%%ymm12, %%ymm1,  %%ymm12\n\t"  // scale by alpha

	"vmovups	%%ymm12,     (%%rdx)	 \n\t"  // store temp -> y

	:
        :
          "m" (n),	// 0	
	  "m" (alpha),  // 1
	  "m" (a),      // 2
          "m" (lda),    // 3
          "m" (x),      // 4
          "m" (y)       // 5
	: "rax", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11",
	  "xmm0" , "xmm1", 
	  "xmm12", "xmm13", "xmm14", "xmm15",
	  "memory"
	);

} 


static void sgemv_kernel_4( long n, float alpha, float *a, long lda, float *x, float *y)
{


	__asm __volatile
	(
	"movq 	        %0, 	 %%rax\n\t"		// n -> rax
	"vbroadcastss   %1, 	 %%xmm1\n\t"		// alpha -> xmm1
	"movq		%2,	 %%rsi\n\t"		// adress of a -> rsi
	"movq	        %3,	 %%rcx\n\t"		// value of lda > rcx
	"movq		%4,	 %%rdi\n\t"		// adress of x -> rdi
	"movq		%5,	 %%rdx\n\t"		// adress of y -> rdx

	"vxorps		%%xmm12, %%xmm12, %%xmm12\n\t"	// set to zero

	".L01LOOP%=:				 \n\t"
	"vbroadcastss	(%%rdi),   %%xmm0	 \n\t"	// load values of c
        "addq		$4     ,   %%rdi	 \n\t"  // increment pointer of c 

	"vfmaddps %%xmm12,   0*4(%%rsi), %%xmm0, %%xmm12\n\t" // multiply a and c and add to temp

	"leaq     (%%rsi, %%rcx, 4), %%rsi	 \n\t"	// add lda to pointer of a

	"dec		%%rax			 \n\t"  // n = n -1
	"jnz		.L01LOOP%=		 \n\t"

	"vmulps		%%xmm12, %%xmm1,  %%xmm12\n\t"  // scale by alpha

	"vmovups	%%xmm12,     (%%rdx)	 \n\t"  // store temp -> y

	:
        :
          "m" (n),	// 0	
	  "m" (alpha),  // 1
	  "m" (a),      // 2
          "m" (lda),    // 3
          "m" (x),      // 4
          "m" (y)       // 5
	: "rax", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11",
	  "xmm0" , "xmm1", 
	  "xmm12", "xmm13", "xmm14", "xmm15",
	  "memory"
	);

} 

static void sgemv_kernel_2( long n, float alpha, float *a, long lda, float *x, float *y)
{


	__asm __volatile
	(
	"movq 	        %0, 	 %%rax\n\t"		// n -> rax
	"vmovss         %1, 	 %%xmm1\n\t"		// alpha -> xmm1
	"movq		%2,	 %%rsi\n\t"		// adress of a -> rsi
	"movq	        %3,	 %%rcx\n\t"		// value of lda > rcx
	"movq		%4,	 %%rdi\n\t"		// adress of x -> rdi
	"movq		%5,	 %%rdx\n\t"		// adress of y -> rdx

	"vxorps		%%xmm12, %%xmm12, %%xmm12\n\t"	// set to zero
	"vxorps		%%xmm13, %%xmm13, %%xmm13\n\t"	// set to zero

	".L01LOOP%=:				 \n\t"
	"vmovss      	(%%rdi),   %%xmm0	 \n\t"	// load values of c
        "addq		$4     ,   %%rdi	 \n\t"  // increment pointer of c 

	"vfmaddss %%xmm12,   0*4(%%rsi), %%xmm0, %%xmm12\n\t" // multiply a and c and add to temp
	"vfmaddss %%xmm13,   1*4(%%rsi), %%xmm0, %%xmm13\n\t" // multiply a and c and add to temp

	"leaq     (%%rsi, %%rcx, 4), %%rsi	 \n\t"	// add lda to pointer of a

	"dec		%%rax			 \n\t"  // n = n -1
	"jnz		.L01LOOP%=		 \n\t"

	"vmulss		%%xmm12, %%xmm1,  %%xmm12\n\t"  // scale by alpha
	"vmulss		%%xmm13, %%xmm1,  %%xmm13\n\t"  // scale by alpha

	"vmovss 	%%xmm12,     (%%rdx)	 \n\t"  // store temp -> y
	"vmovss 	%%xmm13,    4(%%rdx)	 \n\t"  // store temp -> y

	:
        :
          "m" (n),	// 0	
	  "m" (alpha),  // 1
	  "m" (a),      // 2
          "m" (lda),    // 3
          "m" (x),      // 4
          "m" (y)       // 5
	: "rax", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11",
	  "xmm0" , "xmm1", 
	  "xmm12", "xmm13", "xmm14", "xmm15",
	  "memory"
	);

} 



static void sgemv_kernel_1( long n, float alpha, float *a, long lda, float *x, float *y)
{


	__asm __volatile
	(
	"movq 	        %0, 	 %%rax\n\t"		// n -> rax
	"vmovss         %1, 	 %%xmm1\n\t"		// alpha -> xmm1
	"movq		%2,	 %%rsi\n\t"		// adress of a -> rsi
	"movq	        %3,	 %%rcx\n\t"		// value of lda > rcx
	"movq		%4,	 %%rdi\n\t"		// adress of x -> rdi
	"movq		%5,	 %%rdx\n\t"		// adress of y -> rdx

	"vxorps		%%xmm12, %%xmm12, %%xmm12\n\t"	// set to zero

	".L01LOOP%=:				 \n\t"
	"vmovss      	(%%rdi),   %%xmm0	 \n\t"	// load values of c
        "addq		$4     ,   %%rdi	 \n\t"  // increment pointer of c 

	"vfmaddss %%xmm12,   0*4(%%rsi), %%xmm0, %%xmm12\n\t" // multiply a and c and add to temp

	"leaq     (%%rsi, %%rcx, 4), %%rsi	 \n\t"	// add lda to pointer of a

	"dec		%%rax			 \n\t"  // n = n -1
	"jnz		.L01LOOP%=		 \n\t"

	"vmulss		%%xmm12, %%xmm1,  %%xmm12\n\t"  // scale by alpha

	"vmovss 	%%xmm12,     (%%rdx)	 \n\t"  // store temp -> y

	:
        :
          "m" (n),	// 0	
	  "m" (alpha),  // 1
	  "m" (a),      // 2
          "m" (lda),    // 3
          "m" (x),      // 4
          "m" (y)       // 5
	: "rax", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11",
	  "xmm0" , "xmm1", 
	  "xmm12", "xmm13", "xmm14", "xmm15",
	  "memory"
	);

} 


