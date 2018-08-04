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
static void dsymv_kernel_4x4( BLASLONG from, BLASLONG to, FLOAT **a, FLOAT *x, FLOAT *y, FLOAT *temp1, FLOAT *temp2) __attribute__ ((noinline));

static void dsymv_kernel_4x4(BLASLONG from, BLASLONG to, FLOAT **a, FLOAT *x, FLOAT *y, FLOAT *temp1, FLOAT *temp2)
{


	__asm__  __volatile__
	(
	"vzeroupper				     \n\t"
	"vxorpd		%%ymm0 , %%ymm0 , %%ymm0     \n\t"	// temp2[0]
	"vxorpd		%%ymm1 , %%ymm1 , %%ymm1     \n\t"	// temp2[1]
	"vxorpd		%%ymm2 , %%ymm2 , %%ymm2     \n\t"	// temp2[2]
	"vxorpd		%%ymm3 , %%ymm3 , %%ymm3     \n\t"	// temp2[3]
	"vbroadcastsd   (%[temp1]),    %%ymm4	             \n\t"	// temp1[0]
	"vbroadcastsd  8(%[temp1]),    %%ymm5	             \n\t"	// temp1[1]
	"vbroadcastsd 16(%[temp1]),    %%ymm6	             \n\t"	// temp1[1]
	"vbroadcastsd 24(%[temp1]),    %%ymm7	             \n\t"	// temp1[1]

	".p2align 4				     \n\t"
	"1:				     \n\t"

	"vmovups	(%[y],%[from],8), %%ymm9	           \n\t"  // 2 * y
	"vmovups	(%[x],%[from],8), %%ymm8	           \n\t"  // 2 * x

	"vmovups	(%[a0],%[from],8), %%ymm12	           \n\t"  // 2 * a
	"vmovups	(%[a1],%[from],8), %%ymm13	           \n\t"  // 2 * a
	"vmovups	(%[a2],%[from],8), %%ymm14	           \n\t"  // 2 * a
	"vmovups	(%[a3],%[from],8), %%ymm15	           \n\t"  // 2 * a

	"vfmadd231pd	%%ymm4, %%ymm12 , %%ymm9  \n\t"  // y     += temp1 * a
	"vfmadd231pd	%%ymm8, %%ymm12 , %%ymm0  \n\t"  // temp2 += x * a

	"vfmadd231pd	%%ymm5, %%ymm13 , %%ymm9  \n\t"  // y     += temp1 * a
	"vfmadd231pd	%%ymm8, %%ymm13 , %%ymm1  \n\t"  // temp2 += x * a

	"vfmadd231pd	%%ymm6, %%ymm14 , %%ymm9  \n\t"  // y     += temp1 * a
	"vfmadd231pd	%%ymm8, %%ymm14 , %%ymm2  \n\t"  // temp2 += x * a

	"vfmadd231pd	%%ymm7, %%ymm15 , %%ymm9  \n\t"  // y     += temp1 * a
	"vfmadd231pd	%%ymm8, %%ymm15 , %%ymm3  \n\t"  // temp2 += x * a
	"addq		$4 , %[from]	  	 	      \n\t"

	"vmovups	%%ymm9 ,  -32(%[y],%[from],8)		   \n\t"

	"cmpq		%[from] , %[to]			      \n\t"
	"jnz		1b		      \n\t"

	"vmovsd		  (%[temp2]), %%xmm4		      \n\t"
	"vmovsd		 8(%[temp2]), %%xmm5		      \n\t"
	"vmovsd		16(%[temp2]), %%xmm6		      \n\t"
	"vmovsd		24(%[temp2]), %%xmm7		      \n\t"

	"vextractf128 $0x01, %%ymm0 , %%xmm12	      \n\t"
	"vextractf128 $0x01, %%ymm1 , %%xmm13	      \n\t"
	"vextractf128 $0x01, %%ymm2 , %%xmm14	      \n\t"
	"vextractf128 $0x01, %%ymm3 , %%xmm15	      \n\t"

	"vaddpd	        %%xmm0, %%xmm12, %%xmm0	      \n\t"
	"vaddpd	        %%xmm1, %%xmm13, %%xmm1	      \n\t"
	"vaddpd	        %%xmm2, %%xmm14, %%xmm2	      \n\t"
	"vaddpd	        %%xmm3, %%xmm15, %%xmm3	      \n\t"

	"vhaddpd        %%xmm0, %%xmm0, %%xmm0  \n\t"
	"vhaddpd        %%xmm1, %%xmm1, %%xmm1  \n\t"
	"vhaddpd        %%xmm2, %%xmm2, %%xmm2  \n\t"
	"vhaddpd        %%xmm3, %%xmm3, %%xmm3  \n\t"

	"vaddsd		%%xmm4, %%xmm0, %%xmm0  \n\t"
	"vaddsd		%%xmm5, %%xmm1, %%xmm1  \n\t"
	"vaddsd		%%xmm6, %%xmm2, %%xmm2  \n\t"
	"vaddsd		%%xmm7, %%xmm3, %%xmm3  \n\t"

	"vmovsd         %%xmm0 ,  (%[temp2])		\n\t"	// save temp2
	"vmovsd         %%xmm1 , 8(%[temp2])		\n\t"	// save temp2
	"vmovsd         %%xmm2 ,16(%[temp2])		\n\t"	// save temp2
	"vmovsd         %%xmm3 ,24(%[temp2])		\n\t"	// save temp2
	"vzeroupper				     \n\t"

	:
        : 
          [from] "r" (from),	// 0	
	  [to] "r" (to),  	// 1
          [x] "r" (x),      // 2
          [y] "r" (y),      // 3
          [a0] "r" (a[0]),	// 4
          [a1] "r" (a[1]),	// 5
          [a2] "r" (a[2]),	// 6
          [a3] "r" (a[3]),	// 7
          [temp1] "r" (temp1),  // 8
          [temp2] "r" (temp2)   // 9
	: "cc", 
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	);

} 


