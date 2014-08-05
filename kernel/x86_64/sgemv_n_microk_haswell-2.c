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
	"vzeroupper			 \n\t"
	"vbroadcastss    (%2), %%ymm12	 \n\t"	// x0 
	"vbroadcastss   4(%2), %%ymm13	 \n\t"	// x1 
	"vbroadcastss   8(%2), %%ymm14	 \n\t"	// x2 
	"vbroadcastss  12(%2), %%ymm15	 \n\t"	// x3 

	".align 16				 \n\t"
	".L01LOOP%=:				 \n\t"
	"vmovups	(%3,%0,4), %%ymm4	 \n\t"	// 8 * y
	"vmovups      32(%3,%0,4), %%ymm5	 \n\t"	// 8 * y

	"prefetcht0	 192(%4,%0,4)		       \n\t"
	"vfmadd231ps   (%4,%0,4), %%ymm12, %%ymm4      \n\t" 
	"vfmadd231ps 32(%4,%0,4), %%ymm12, %%ymm5      \n\t" 
	"prefetcht0	 192(%5,%0,4)		       \n\t"
	"vfmadd231ps   (%5,%0,4), %%ymm13, %%ymm4      \n\t" 
	"vfmadd231ps 32(%5,%0,4), %%ymm13, %%ymm5      \n\t" 
	"prefetcht0	 192(%6,%0,4)		       \n\t"
	"vfmadd231ps   (%6,%0,4), %%ymm14, %%ymm4      \n\t" 
	"vfmadd231ps 32(%6,%0,4), %%ymm14, %%ymm5      \n\t" 
	"prefetcht0	 192(%7,%0,4)		       \n\t"
	"vfmadd231ps   (%7,%0,4), %%ymm15, %%ymm4      \n\t" 
	"vfmadd231ps 32(%7,%0,4), %%ymm15, %%ymm5      \n\t" 

	"vmovups  %%ymm4,   (%3,%0,4)		      \n\t"	// 8 * y
	"vmovups  %%ymm5, 32(%3,%0,4)		      \n\t"	// 8 * y

        "addq		$16, %0	  	 	      \n\t"
	"subq	        $16, %1			      \n\t"		
	"jnz		.L01LOOP%=		      \n\t"
	"vzeroupper			 \n\t"

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
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	);

} 


