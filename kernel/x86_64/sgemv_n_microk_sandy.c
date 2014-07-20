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

static void  sgemv_kernel_64( long n, float alpha, float *a, long lda, float *x, float *y)
{


	float *pre = a + lda*2;

	__asm__  __volatile__
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

	"vxorps		%%ymm8 , %%ymm8 , %%ymm8 \n\t"	// set to zero
	"vxorps		%%ymm9 , %%ymm9 , %%ymm9 \n\t"	// set to zero
	"vxorps		%%ymm10, %%ymm10, %%ymm10\n\t"	// set to zero
	"vxorps		%%ymm11, %%ymm11, %%ymm11\n\t"	// set to zero
	"vxorps		%%ymm12, %%ymm12, %%ymm12\n\t"	// set to zero
	"vxorps		%%ymm13, %%ymm13, %%ymm13\n\t"	// set to zero
	"vxorps		%%ymm14, %%ymm14, %%ymm14\n\t"	// set to zero
	"vxorps		%%ymm15, %%ymm15, %%ymm15\n\t"	// set to zero
	".align 16				 \n\t"
	".L01LOOP%=:				 \n\t"
	"vbroadcastss	(%%rdi),   %%ymm0	 \n\t"	// load values of c
	"nop					 \n\t"
	"leaq     (%%r8 , %%rcx, 4), %%r8	 \n\t"	// add lda to pointer for prefetch

	"prefetcht0	(%%r8)\n\t"			// Prefetch
	"vmulps   0*4(%%rsi), %%ymm0, %%ymm4 \n\t" // multiply a and c and add to temp
	"vmulps   8*4(%%rsi), %%ymm0, %%ymm5 \n\t" // multiply a and c and add to temp
	"prefetcht0   64(%%r8)\n\t"			// Prefetch
	"vmulps  16*4(%%rsi), %%ymm0, %%ymm6 \n\t" // multiply a and c and add to temp
	"vmulps  24*4(%%rsi), %%ymm0, %%ymm7 \n\t" // multiply a and c and add to temp

	"vaddps %%ymm8 , %%ymm4, %%ymm8 \n\t" // multiply a and c and add to temp
	"vaddps %%ymm9 , %%ymm5, %%ymm9 \n\t" // multiply a and c and add to temp
	"prefetcht0  128(%%r8)\n\t"			// Prefetch
	"vaddps %%ymm10, %%ymm6, %%ymm10\n\t" // multiply a and c and add to temp
	"vaddps %%ymm11, %%ymm7, %%ymm11\n\t" // multiply a and c and add to temp

	"prefetcht0  192(%%r8)\n\t"			// Prefetch
	"vmulps  32*4(%%rsi), %%ymm0, %%ymm4 \n\t" // multiply a and c and add to temp
	"vmulps  40*4(%%rsi), %%ymm0, %%ymm5 \n\t" // multiply a and c and add to temp
	"vmulps  48*4(%%rsi), %%ymm0, %%ymm6 \n\t" // multiply a and c and add to temp
	"vmulps  56*4(%%rsi), %%ymm0, %%ymm7 \n\t" // multiply a and c and add to temp

	"vaddps %%ymm12, %%ymm4, %%ymm12\n\t" // multiply a and c and add to temp
	"vaddps %%ymm13, %%ymm5, %%ymm13\n\t" // multiply a and c and add to temp
	"vaddps %%ymm14, %%ymm6, %%ymm14\n\t" // multiply a and c and add to temp
	"vaddps %%ymm15, %%ymm7, %%ymm15\n\t" // multiply a and c and add to temp

        "addq		$4     ,   %%rdi	 \n\t"  // increment pointer of c 
	"leaq     (%%rsi, %%rcx, 4), %%rsi	 \n\t"	// add lda to pointer of a

	"dec		%%rax			 \n\t"  // n = n -1
	"jnz		.L01LOOP%=		 \n\t"

	"vmulps		%%ymm8 , %%ymm1,  %%ymm8 \n\t"  // scale by alpha
	"vmulps		%%ymm9 , %%ymm1,  %%ymm9 \n\t"  // scale by alpha
	"vmulps		%%ymm10, %%ymm1,  %%ymm10\n\t"  // scale by alpha
	"vmulps		%%ymm11, %%ymm1,  %%ymm11\n\t"  // scale by alpha
	"vmulps		%%ymm12, %%ymm1,  %%ymm12\n\t"  // scale by alpha
	"vmulps		%%ymm13, %%ymm1,  %%ymm13\n\t"  // scale by alpha
	"vmulps		%%ymm14, %%ymm1,  %%ymm14\n\t"  // scale by alpha
	"vmulps		%%ymm15, %%ymm1,  %%ymm15\n\t"  // scale by alpha

	"vmovups	%%ymm8 ,     (%%rdx)	 \n\t"  // store temp -> y
	"vmovups	%%ymm9 ,  8*4(%%rdx)	 \n\t"  // store temp -> y
	"vmovups	%%ymm10, 16*4(%%rdx)	 \n\t"  // store temp -> y
	"vmovups	%%ymm11, 24*4(%%rdx)	 \n\t"  // store temp -> y
	"vmovups	%%ymm12, 32*4(%%rdx)	 \n\t"  // store temp -> y
	"vmovups	%%ymm13, 40*4(%%rdx)	 \n\t"  // store temp -> y
	"vmovups	%%ymm14, 48*4(%%rdx)	 \n\t"  // store temp -> y
	"vmovups	%%ymm15, 56*4(%%rdx)	 \n\t"  // store temp -> y

	:
        :
          "m" (n),	// 0	
	  "m" (alpha),  // 1
	  "m" (a),      // 2
          "m" (lda),    // 3
          "m" (x),      // 4
          "m" (y),      // 5
	  "m" (pre)	// 6
	: "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "cc",
	  "%xmm0", "%xmm1", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7",
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11",
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	);

} 



static void  sgemv_kernel_32( long n, float alpha, float *a, long lda, float *x, float *y)
{


	float *pre = a + lda*3;

	__asm__  __volatile__
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

	"vxorps		%%ymm8 , %%ymm8 , %%ymm8 \n\t"	// set to zero
	"vxorps		%%ymm9 , %%ymm9 , %%ymm9 \n\t"	// set to zero
	"vxorps		%%ymm10, %%ymm10, %%ymm10\n\t"	// set to zero
	"vxorps		%%ymm11, %%ymm11, %%ymm11\n\t"	// set to zero
	".align 16				 \n\t"
	".L01LOOP%=:				 \n\t"
	"vbroadcastss	(%%rdi),   %%ymm0	 \n\t"	// load values of c
	"nop					 \n\t"
	"leaq     (%%r8 , %%rcx, 4), %%r8	 \n\t"	// add lda to pointer for prefetch

	"prefetcht0	(%%r8)\n\t"			// Prefetch
	"prefetcht0   64(%%r8)\n\t"			// Prefetch

	"vmulps   0*4(%%rsi), %%ymm0, %%ymm4 \n\t" // multiply a and c and add to temp
	"vmulps   8*4(%%rsi), %%ymm0, %%ymm5 \n\t" // multiply a and c and add to temp
	"vmulps  16*4(%%rsi), %%ymm0, %%ymm6 \n\t" // multiply a and c and add to temp
	"vmulps  24*4(%%rsi), %%ymm0, %%ymm7 \n\t" // multiply a and c and add to temp

	"vaddps %%ymm8 , %%ymm4, %%ymm8 \n\t" // multiply a and c and add to temp
	"vaddps %%ymm9 , %%ymm5, %%ymm9 \n\t" // multiply a and c and add to temp
	"vaddps %%ymm10, %%ymm6, %%ymm10\n\t" // multiply a and c and add to temp
	"vaddps %%ymm11, %%ymm7, %%ymm11\n\t" // multiply a and c and add to temp



        "addq		$4     ,   %%rdi	 \n\t"  // increment pointer of c 
	"leaq     (%%rsi, %%rcx, 4), %%rsi	 \n\t"	// add lda to pointer of a

	"dec		%%rax			 \n\t"  // n = n -1
	"jnz		.L01LOOP%=		 \n\t"

	"vmulps		%%ymm8 , %%ymm1,  %%ymm8 \n\t"  // scale by alpha
	"vmulps		%%ymm9 , %%ymm1,  %%ymm9 \n\t"  // scale by alpha
	"vmulps		%%ymm10, %%ymm1,  %%ymm10\n\t"  // scale by alpha
	"vmulps		%%ymm11, %%ymm1,  %%ymm11\n\t"  // scale by alpha

	"vmovups	%%ymm8 ,     (%%rdx)	 \n\t"  // store temp -> y
	"vmovups	%%ymm9 ,  8*4(%%rdx)	 \n\t"  // store temp -> y
	"vmovups	%%ymm10, 16*4(%%rdx)	 \n\t"  // store temp -> y
	"vmovups	%%ymm11, 24*4(%%rdx)	 \n\t"  // store temp -> y

	:
        :
          "m" (n),	// 0	
	  "m" (alpha),  // 1
	  "m" (a),      // 2
          "m" (lda),    // 3
          "m" (x),      // 4
          "m" (y),      // 5
	  "m" (pre)	// 6
	: "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "cc",
	  "%xmm0", "%xmm1", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7",
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11",
	  "memory"
	);



} 

static void  sgemv_kernel_16( long n, float alpha, float *a, long lda, float *x, float *y)
{

	float *pre = a + lda*3;
	
	__asm__  __volatile__
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

	"vxorps		%%ymm8 , %%ymm8 , %%ymm8 \n\t"	// set to zero
	"vxorps		%%ymm9 , %%ymm9 , %%ymm9 \n\t"	// set to zero
	".align 16				 \n\t"
	".L01LOOP%=:				 \n\t"
	"vbroadcastss	(%%rdi),   %%ymm0	 \n\t"	// load values of c
	"nop					 \n\t"
	"leaq     (%%r8 , %%rcx, 4), %%r8	 \n\t"	// add lda to pointer for prefetch

	"prefetcht0	(%%r8)\n\t"			// Prefetch

	"vmulps   0*4(%%rsi), %%ymm0, %%ymm4 \n\t" // multiply a and c and add to temp
	"vmulps   8*4(%%rsi), %%ymm0, %%ymm5 \n\t" // multiply a and c and add to temp

	"vaddps %%ymm8 , %%ymm4, %%ymm8 \n\t" // multiply a and c and add to temp
	"vaddps %%ymm9 , %%ymm5, %%ymm9 \n\t" // multiply a and c and add to temp

        "addq		$4     ,   %%rdi	 \n\t"  // increment pointer of c 
	"leaq     (%%rsi, %%rcx, 4), %%rsi	 \n\t"	// add lda to pointer of a

	"dec		%%rax			 \n\t"  // n = n -1
	"jnz		.L01LOOP%=		 \n\t"

	"vmulps		%%ymm8 , %%ymm1,  %%ymm8 \n\t"  // scale by alpha
	"vmulps		%%ymm9 , %%ymm1,  %%ymm9 \n\t"  // scale by alpha

	"vmovups	%%ymm8 ,     (%%rdx)	 \n\t"  // store temp -> y
	"vmovups	%%ymm9 ,  8*4(%%rdx)	 \n\t"  // store temp -> y

	:
        :
          "m" (n),	// 0	
	  "m" (alpha),  // 1
	  "m" (a),      // 2
          "m" (lda),    // 3
          "m" (x),      // 4
          "m" (y),      // 5
	  "m" (pre)	// 6
	: "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "cc",
	  "%xmm0", "%xmm1", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7",
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11",
	  "memory"
	);


} 


static void  sgemv_kernel_8( long n, float alpha, float *a, long lda, float *x, float *y)
{
	
	__asm__  __volatile__
	(
	"movq 	        %0, 	 %%rax\n\t"		// n -> rax
	"vbroadcastss   %1, 	 %%ymm1\n\t"		// alpha -> ymm1
	"movq		%2,	 %%rsi\n\t"		// adress of a -> rsi
	"movq	        %3,	 %%rcx\n\t"		// value of lda > rcx
	"movq		%4,	 %%rdi\n\t"		// adress of x -> rdi
	"movq		%5,	 %%rdx\n\t"		// adress of y -> rdx

	"vxorps		%%ymm8 , %%ymm8 , %%ymm8 \n\t"	// set to zero
	".align 16				 \n\t"
	".L01LOOP%=:				 \n\t"
	"vbroadcastss	(%%rdi),   %%ymm0	 \n\t"	// load values of c

	"vmulps   0*4(%%rsi), %%ymm0, %%ymm4 \n\t" // multiply a and c and add to temp
	"vaddps %%ymm8 , %%ymm4, %%ymm8 \n\t" // multiply a and c and add to temp

        "addq		$4     ,   %%rdi	 \n\t"  // increment pointer of c 
	"leaq     (%%rsi, %%rcx, 4), %%rsi	 \n\t"	// add lda to pointer of a

	"dec		%%rax			 \n\t"  // n = n -1
	"jnz		.L01LOOP%=		 \n\t"

	"vmulps		%%ymm8 , %%ymm1,  %%ymm8 \n\t"  // scale by alpha
	"vmovups	%%ymm8 ,     (%%rdx)	 \n\t"  // store temp -> y

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
	  "memory"
	);


} 


static void  sgemv_kernel_4( long n, float alpha, float *a, long lda, float *x, float *y)
{


	__asm__ __volatile__
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

	"vmulps   0*4(%%rsi), %%xmm0, %%xmm4 \n\t" 		// multiply a and c and add to temp
	"vaddps %%xmm12, %%xmm4, %%xmm12     \n\t" 		// multiply a and c and add to temp

        "addq		$4     ,   %%rdi	 \n\t"  // increment pointer of c 
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
	: "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", 
	  "%xmm0", "%xmm1", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11",
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	);

} 

static void  sgemv_kernel_2( long n, float alpha, float *a, long lda, float *x, float *y)
{


	__asm__ __volatile__
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

	"vmulps   0*4(%%rsi), %%xmm0, %%xmm4 \n\t" 		// multiply a and c and add to temp
	"vmulps   1*4(%%rsi), %%xmm0, %%xmm5 \n\t" 		// multiply a and c and add to temp

	"vaddps %%xmm12, %%xmm4, %%xmm12     \n\t" 		// multiply a and c and add to temp
	"vaddps %%xmm13, %%xmm5, %%xmm13     \n\t" 		// multiply a and c and add to temp

        "addq		$4     ,   %%rdi	 \n\t"  // increment pointer of c 
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
	: "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", 
	  "%xmm0", "%xmm1", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11",
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	);

} 



static void  sgemv_kernel_1( long n, float alpha, float *a, long lda, float *x, float *y)
{


	__asm__ __volatile__
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

	"vmulss   0*4(%%rsi), %%xmm0, %%xmm4 \n\t" 		// multiply a and c and add to temp
	"vaddss %%xmm12, %%xmm4, %%xmm12     \n\t" 		// multiply a and c and add to temp

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
	: "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", 
	  "%xmm0", "%xmm1", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11",
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	);

} 


