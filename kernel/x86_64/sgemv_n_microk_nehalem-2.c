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
static void sgemv_kernel_16x4( long n, float **ap, float *x, float *y) __attribute__ ((noinline));

static void sgemv_kernel_16x4( long n, float **ap, float *x, float *y)
{

	long register i = 0;

	__asm__  __volatile__
	(
	"movss    (%2), %%xmm12	 \n\t"	// x0 
	"movss   4(%2), %%xmm13	 \n\t"	// x1 
	"movss   8(%2), %%xmm14	 \n\t"	// x2 
	"movss  12(%2), %%xmm15	 \n\t"	// x3 
	"shufps $0,  %%xmm12, %%xmm12\n\t"	
	"shufps $0,  %%xmm13, %%xmm13\n\t"	
	"shufps $0,  %%xmm14, %%xmm14\n\t"	
	"shufps $0,  %%xmm15, %%xmm15\n\t"	

	".align 16				 \n\t"
	".L01LOOP%=:				 \n\t"
	"movups	       (%3,%0,4), %%xmm4	 \n\t"	// 4 * y
	"movups      16(%3,%0,4), %%xmm5	 \n\t"	// 4 * y
	"movups      32(%3,%0,4), %%xmm6	 \n\t"	// 4 * y
	"movups      48(%3,%0,4), %%xmm7	 \n\t"	// 4 * y

	"prefetcht0	 192(%4,%0,4)		       \n\t"

	"movups             (%4,%0,4), %%xmm8          \n\t" 
	"movups           16(%4,%0,4), %%xmm9          \n\t" 
	"movups           32(%4,%0,4), %%xmm10         \n\t" 
	"movups           48(%4,%0,4), %%xmm11         \n\t" 
	"mulps		%%xmm12, %%xmm8		       \n\t"
	"addps		%%xmm8 , %%xmm4		       \n\t"
	"mulps		%%xmm12, %%xmm9		       \n\t"
	"addps		%%xmm9 , %%xmm5		       \n\t"
	"mulps		%%xmm12, %%xmm10	       \n\t"
	"addps		%%xmm10, %%xmm6		       \n\t"
	"mulps		%%xmm12, %%xmm11	       \n\t"
	"addps		%%xmm11, %%xmm7		       \n\t"

	"prefetcht0	 192(%5,%0,4)		       \n\t"

	"movups             (%5,%0,4), %%xmm8          \n\t" 
	"movups           16(%5,%0,4), %%xmm9          \n\t" 
	"movups           32(%5,%0,4), %%xmm10         \n\t" 
	"movups           48(%5,%0,4), %%xmm11         \n\t" 
	"mulps		%%xmm13, %%xmm8		       \n\t"
	"addps		%%xmm8 , %%xmm4		       \n\t"
	"mulps		%%xmm13, %%xmm9		       \n\t"
	"addps		%%xmm9 , %%xmm5		       \n\t"
	"mulps		%%xmm13, %%xmm10	       \n\t"
	"addps		%%xmm10, %%xmm6		       \n\t"
	"mulps		%%xmm13, %%xmm11	       \n\t"
	"addps		%%xmm11, %%xmm7		       \n\t"

	"prefetcht0	 192(%6,%0,4)		       \n\t"

	"movups             (%6,%0,4), %%xmm8          \n\t" 
	"movups           16(%6,%0,4), %%xmm9          \n\t" 
	"movups           32(%6,%0,4), %%xmm10         \n\t" 
	"movups           48(%6,%0,4), %%xmm11         \n\t" 
	"mulps		%%xmm14, %%xmm8		       \n\t"
	"addps		%%xmm8 , %%xmm4		       \n\t"
	"mulps		%%xmm14, %%xmm9		       \n\t"
	"addps		%%xmm9 , %%xmm5		       \n\t"
	"mulps		%%xmm14, %%xmm10	       \n\t"
	"addps		%%xmm10, %%xmm6		       \n\t"
	"mulps		%%xmm14, %%xmm11	       \n\t"
	"addps		%%xmm11, %%xmm7		       \n\t"

	"prefetcht0	 192(%7,%0,4)		       \n\t"

	"movups             (%7,%0,4), %%xmm8          \n\t" 
	"movups           16(%7,%0,4), %%xmm9          \n\t" 
	"movups           32(%7,%0,4), %%xmm10         \n\t" 
	"movups           48(%7,%0,4), %%xmm11         \n\t" 
	"mulps		%%xmm15, %%xmm8		       \n\t"
	"addps		%%xmm8 , %%xmm4		       \n\t"
	"mulps		%%xmm15, %%xmm9		       \n\t"
	"addps		%%xmm9 , %%xmm5		       \n\t"
	"mulps		%%xmm15, %%xmm10	       \n\t"
	"addps		%%xmm10, %%xmm6		       \n\t"
	"mulps		%%xmm15, %%xmm11	       \n\t"
	"addps		%%xmm11, %%xmm7		       \n\t"


	"movups  %%xmm4,   (%3,%0,4)		      \n\t"	// 4 * y
	"movups  %%xmm5, 16(%3,%0,4)		      \n\t"	// 4 * y
	"movups  %%xmm6, 32(%3,%0,4)		      \n\t"	// 4 * y
	"movups  %%xmm7, 48(%3,%0,4)		      \n\t"	// 4 * y

        "addq		$16, %0	  	 	      \n\t"
	"subq	        $16, %1			      \n\t"		
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
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11",
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	);

} 


