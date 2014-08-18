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

#define HAVE_KERNEL_8x2 1
static void dsymv_kernel_8x2( BLASLONG n, FLOAT *a1, FLOAT *a2, FLOAT *x, FLOAT *y, FLOAT *temp1, FLOAT *temp2) __attribute__ ((noinline));

static void dsymv_kernel_8x2(BLASLONG n, FLOAT *a0, FLOAT *a1, FLOAT *x, FLOAT *y, FLOAT *temp1, FLOAT *temp2)
{

	BLASLONG register i = 0;

	__asm__  __volatile__
	(
	"vxorpd		%%xmm0 , %%xmm0 , %%xmm0     \n\t"	// temp2[0]
	"vxorpd		%%xmm1 , %%xmm1 , %%xmm1     \n\t"	// temp2[1]
	"vmovddup	(%6),    %%xmm2	             \n\t"	// temp1[0]
	"vmovddup      8(%6),    %%xmm3	             \n\t"	// temp1[1]

	"xorq		%0,%0			     \n\t"

	".align 16				 \n\t"
	".L01LOOP%=:				 \n\t"

	"prefetcht0	 192(%4,%0,8)		       \n\t"
	"vmovups	(%4,%0,8), %%xmm4	 \n\t"	// 2 * a0
	"vmovups      16(%4,%0,8), %%xmm5	 \n\t"	// 2 * a0
	"prefetcht0	 192(%2,%0,8)		       \n\t"
	"vmovups	(%2,%0,8), %%xmm8	 \n\t"	// 2 * x
	"vmovups      16(%2,%0,8), %%xmm9	 \n\t"	// 2 * x
	"prefetcht0	 192(%3,%0,8)		       \n\t"
	"vmovups      32(%4,%0,8), %%xmm6	 \n\t"	// 2 * a0
	"vmovups      48(%4,%0,8), %%xmm7	 \n\t"	// 2 * a0
	"vmovups      32(%2,%0,8), %%xmm10	 \n\t"	// 2 * x
	"vmovups      48(%2,%0,8), %%xmm11	 \n\t"	// 2 * x

	"prefetcht0	 192(%5,%0,8)		            \n\t"
	"vfmaddpd    (%3,%0,8), %%xmm2 , %%xmm4 , %%xmm12   \n\t" // y += temp1 * a0
	"vfmaddpd      %%xmm0 , %%xmm8 , %%xmm4 , %%xmm0    \n\t" // temp2 += a0 * x
	"vfmaddpd  16(%3,%0,8), %%xmm2 , %%xmm5 , %%xmm13   \n\t" // y += temp1 * a0
	"vmovups	(%5,%0,8), %%xmm4	            \n\t"	// 2 * a1
	"vfmaddpd      %%xmm0 , %%xmm9 , %%xmm5 , %%xmm0    \n\t" // temp2 += a0 * x
	"vfmaddpd  32(%3,%0,8), %%xmm2 , %%xmm6 , %%xmm14   \n\t" // y += temp1 * a0
	"vmovups      16(%5,%0,8), %%xmm5	            \n\t"	// 2 * a1
	"vfmaddpd      %%xmm0 , %%xmm10, %%xmm6 , %%xmm0    \n\t" // temp2 += a0 * x
	"vfmaddpd  48(%3,%0,8), %%xmm2 , %%xmm7 , %%xmm15   \n\t" // y += temp1 * a0
	"vmovups      32(%5,%0,8), %%xmm6	            \n\t"	// 2 * a1
	"vfmaddpd      %%xmm0 , %%xmm11, %%xmm7 , %%xmm0    \n\t" // temp2 += a0 * x
	"vmovups      48(%5,%0,8), %%xmm7	            \n\t"	// 2 * a1

	"vfmaddpd %%xmm12, %%xmm3 , %%xmm4 , %%xmm12   \n\t" // y += temp1 * a1
	"vfmaddpd %%xmm13, %%xmm3 , %%xmm5 , %%xmm13   \n\t" // y += temp1 * a1
	"vmovups  %%xmm12,   (%3,%0,8)		       \n\t"	// 2 * y
	"vfmaddpd %%xmm14, %%xmm3 , %%xmm6 , %%xmm14   \n\t" // y += temp1 * a1
	"vmovups  %%xmm13, 16(%3,%0,8)		       \n\t"	// 2 * y
	"vfmaddpd %%xmm15, %%xmm3 , %%xmm7 , %%xmm15   \n\t" // y += temp1 * a1
	"vmovups  %%xmm14, 32(%3,%0,8)		       \n\t"	// 2 * y

	"vfmaddpd %%xmm1 , %%xmm8 , %%xmm4 , %%xmm1   \n\t" // temp2 += a1 * x
	"vfmaddpd %%xmm1 , %%xmm9 , %%xmm5 , %%xmm1   \n\t" // temp2 += a1 * x
	"vmovups  %%xmm15, 48(%3,%0,8)		      \n\t"	// 2 * y
	"vfmaddpd %%xmm1 , %%xmm10, %%xmm6 , %%xmm1   \n\t" // temp2 += a1 * x
	"vfmaddpd %%xmm1 , %%xmm11, %%xmm7 , %%xmm1   \n\t" // temp2 += a1 * x

        "addq		$8, %0	  	 	      \n\t"
	"subq	        $8, %1			      \n\t"		
	"jnz		.L01LOOP%=		      \n\t"

	"vhaddpd        %%xmm0, %%xmm0, %%xmm0  \n\t"
	"vhaddpd        %%xmm1, %%xmm1, %%xmm1  \n\t"
	"vmovsd         %%xmm0 , (%7)		\n\t"	// save temp2
	"vmovsd         %%xmm1 ,8(%7)		\n\t"	// save temp2

	:
        : 
          "r" (i),	// 0	
	  "r" (n),  	// 1
          "r" (x),      // 2
          "r" (y),      // 3
          "r" (a0),  // 4
          "r" (a1),  // 5
          "r" (temp1),  // 6
          "r" (temp2)   // 7
	: "cc", 
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	);

} 


