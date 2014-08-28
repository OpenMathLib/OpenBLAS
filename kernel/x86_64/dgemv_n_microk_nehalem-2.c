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

#define HAVE_KERNEL_16x4 1
static void dgemv_kernel_16x4( BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y) __attribute__ ((noinline));

static void dgemv_kernel_16x4( BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y)
{

	BLASLONG register i = 0;

	__asm__  __volatile__
	(
	"movddup    (%2), %%xmm12	 \n\t"	// x0 
	"movddup   8(%2), %%xmm13	 \n\t"	// x1 
	"movddup  16(%2), %%xmm14	 \n\t"	// x2 
	"movddup  24(%2), %%xmm15	 \n\t"	// x3 

	".align 16				 \n\t"
	".L01LOOP%=:				 \n\t"
	"prefetcht0	 192(%3,%0,8)		 \n\t"
	"movups	       (%3,%0,8), %%xmm4	 \n\t"	// 2 * y
	"movups      16(%3,%0,8), %%xmm5	 \n\t"	// 2 * y
	"movups      32(%3,%0,8), %%xmm6	 \n\t"	// 2 * y
	"movups      48(%3,%0,8), %%xmm7	 \n\t"	// 2 * y
	"movups	       (%4,%0,8), %%xmm8	 \n\t"	// 2 * a
	"movups      16(%4,%0,8), %%xmm9	 \n\t"	// 2 * a
	"movups      32(%4,%0,8), %%xmm10        \n\t"	// 2 * a
	"movups      48(%4,%0,8), %%xmm11        \n\t"  // 2 * a

	"prefetcht0	 192(%4,%0,8)		       \n\t"
	"mulpd           %%xmm12 , %%xmm8	       \n\t" // a * x
	"mulpd           %%xmm12 , %%xmm9	       \n\t" // a * x
	"mulpd           %%xmm12 , %%xmm10             \n\t" // a * x
	"mulpd           %%xmm12 , %%xmm11             \n\t" // a * x
	"addpd		 %%xmm8  , %%xmm4              \n\t" // y += a * x
	"addpd		 %%xmm9  , %%xmm5              \n\t" // y += a * x
	"addpd		 %%xmm10 , %%xmm6              \n\t" // y += a * x
	"addpd		 %%xmm11 , %%xmm7              \n\t" // y += a * x

	"prefetcht0	 192(%5,%0,8)		       \n\t"
	"movups	       (%5,%0,8), %%xmm8	       \n\t" // 2 * a
	"movups      16(%5,%0,8), %%xmm9	       \n\t" // 2 * a
	"movups      32(%5,%0,8), %%xmm10              \n\t" // 2 * a
	"movups      48(%5,%0,8), %%xmm11              \n\t" // 2 * a
	"mulpd           %%xmm13 , %%xmm8	       \n\t" // a * x
	"mulpd           %%xmm13 , %%xmm9	       \n\t" // a * x
	"mulpd           %%xmm13 , %%xmm10             \n\t" // a * x
	"mulpd           %%xmm13 , %%xmm11             \n\t" // a * x
	"addpd		 %%xmm8  , %%xmm4              \n\t" // y += a * x
	"addpd		 %%xmm9  , %%xmm5              \n\t" // y += a * x
	"addpd		 %%xmm10 , %%xmm6              \n\t" // y += a * x
	"addpd		 %%xmm11 , %%xmm7              \n\t" // y += a * x

	"prefetcht0	 192(%6,%0,8)		       \n\t"
	"movups	       (%6,%0,8), %%xmm8	       \n\t" // 2 * a
	"movups      16(%6,%0,8), %%xmm9	       \n\t" // 2 * a
	"movups      32(%6,%0,8), %%xmm10              \n\t" // 2 * a
	"movups      48(%6,%0,8), %%xmm11              \n\t" // 2 * a
	"mulpd           %%xmm14 , %%xmm8	       \n\t" // a * x
	"mulpd           %%xmm14 , %%xmm9	       \n\t" // a * x
	"mulpd           %%xmm14 , %%xmm10             \n\t" // a * x
	"mulpd           %%xmm14 , %%xmm11             \n\t" // a * x
	"addpd		 %%xmm8  , %%xmm4              \n\t" // y += a * x
	"addpd		 %%xmm9  , %%xmm5              \n\t" // y += a * x
	"addpd		 %%xmm10 , %%xmm6              \n\t" // y += a * x
	"addpd		 %%xmm11 , %%xmm7              \n\t" // y += a * x

	"prefetcht0	 192(%7,%0,8)		       \n\t"
	"movups	       (%7,%0,8), %%xmm8	       \n\t" // 2 * a
	"movups      16(%7,%0,8), %%xmm9	       \n\t" // 2 * a
	"movups      32(%7,%0,8), %%xmm10              \n\t" // 2 * a
	"movups      48(%7,%0,8), %%xmm11              \n\t" // 2 * a
	"mulpd           %%xmm15 , %%xmm8	       \n\t" // a * x
	"mulpd           %%xmm15 , %%xmm9	       \n\t" // a * x
	"mulpd           %%xmm15 , %%xmm10             \n\t" // a * x
	"mulpd           %%xmm15 , %%xmm11             \n\t" // a * x
	"addpd		 %%xmm8  , %%xmm4              \n\t" // y += a * x
	"addpd		 %%xmm9  , %%xmm5              \n\t" // y += a * x
	"addpd		 %%xmm10 , %%xmm6              \n\t" // y += a * x
	"addpd		 %%xmm11 , %%xmm7              \n\t" // y += a * x

	"movups  %%xmm4,   (%3,%0,8)		      \n\t"	// 4 * y
	"movups  %%xmm5, 16(%3,%0,8)		      \n\t"	// 4 * y
	"movups  %%xmm6, 32(%3,%0,8)		      \n\t"	// 4 * y
	"movups  %%xmm7, 48(%3,%0,8)		      \n\t"	// 4 * y

        "addq		$8 , %0	  	 	      \n\t"
	"subq	        $8 , %1			      \n\t"		
	"jnz		.L01LOOP%=		      \n\t"

	:
        : 
          "r" (i),	// 0	
	  "r" (n),  	// 1
          "r" (x),      // 2
          "r" (y),      // 3
          "r" (ap[0]),  // 4
          "r" (ap[1]),  // 5
          "r" (ap[2]),  // 6
          "r" (ap[3])   // 7
	: "cc", 
	  "%xmm4", "%xmm5", 
	  "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", 
	  "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	);

} 


