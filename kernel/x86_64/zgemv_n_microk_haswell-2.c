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
static void zgemv_kernel_16x4( BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y) __attribute__ ((noinline));

static void zgemv_kernel_16x4( BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y)
{

	BLASLONG register i = 0;

	__asm__  __volatile__
	(
	"vzeroupper			 \n\t"

	"vbroadcastsd	  (%2), %%ymm0                  \n\t"  // real part x0
	"vbroadcastsd	 8(%2), %%ymm1                  \n\t"  // imag part x0
	"vbroadcastsd	16(%2), %%ymm2                  \n\t"  // real part x1
	"vbroadcastsd	24(%2), %%ymm3                  \n\t"  // imag part x1
	"vbroadcastsd	32(%2), %%ymm4                  \n\t"  // real part x2
	"vbroadcastsd	40(%2), %%ymm5                  \n\t"  // imag part x2
	"vbroadcastsd	48(%2), %%ymm6                  \n\t"  // real part x3
	"vbroadcastsd	56(%2), %%ymm7                  \n\t"  // imag part x3


	".align 16				        \n\t"
	".L01LOOP%=:				        \n\t"
	"prefetcht0      192(%4,%0,8)			\n\t"
	"vmovups	(%4,%0,8), %%ymm8	        \n\t" // 2 complex values form a0
	"vmovups      32(%4,%0,8), %%ymm9	        \n\t" // 2 complex values form a0

	"vxorpd		%%ymm12, %%ymm12, %%ymm12	\n\t"
	"vxorpd		%%ymm13, %%ymm13, %%ymm13	\n\t"
	"vxorpd		%%ymm14, %%ymm14, %%ymm14	\n\t"
	"vxorpd		%%ymm15, %%ymm15, %%ymm15	\n\t"

	"prefetcht0      192(%5,%0,8)			\n\t"
	"vmovups	(%5,%0,8), %%ymm10              \n\t" // 2 complex values form a1
	"vmovups      32(%5,%0,8), %%ymm11              \n\t" // 2 complex values form a1

	"vfmadd231pd      %%ymm8 , %%ymm0, %%ymm12      \n\t" // a_r[0] * x_r , a_i[0] * x_r, a_r[1] * x_r, a_i[1] * x_r
	"vfmadd231pd      %%ymm8 , %%ymm1, %%ymm13      \n\t" // a_r[0] * x_i , a_i[0] * x_i, a_r[1] * x_i, a_i[1] * x_i
	"vfmadd231pd      %%ymm9 , %%ymm0, %%ymm14      \n\t" // a_r[2] * x_r , a_i[2] * x_r, a_r[3] * x_r, a_i[3] * x_r
	"vfmadd231pd      %%ymm9 , %%ymm1, %%ymm15      \n\t" // a_r[2] * x_i , a_i[2] * x_i, a_r[3] * x_i, a_i[3] * x_i

	"prefetcht0      192(%6,%0,8)			\n\t"
	"vmovups	(%6,%0,8), %%ymm8	        \n\t" // 2 complex values form a2
	"vmovups      32(%6,%0,8), %%ymm9	        \n\t" // 2 complex values form a2

	"vfmadd231pd      %%ymm10, %%ymm2, %%ymm12      \n\t" // a_r[0] * x_r , a_i[0] * x_r, a_r[1] * x_r, a_i[1] * x_r
	"vfmadd231pd      %%ymm10, %%ymm3, %%ymm13      \n\t" // a_r[0] * x_i , a_i[0] * x_i, a_r[1] * x_i, a_i[1] * x_i
	"vfmadd231pd      %%ymm11, %%ymm2, %%ymm14      \n\t" // a_r[2] * x_r , a_i[2] * x_r, a_r[3] * x_r, a_i[3] * x_r
	"vfmadd231pd      %%ymm11, %%ymm3, %%ymm15      \n\t" // a_r[2] * x_i , a_i[2] * x_i, a_r[3] * x_i, a_i[3] * x_i

	"prefetcht0      192(%7,%0,8)			\n\t"
	"vmovups	(%7,%0,8), %%ymm10              \n\t" // 2 complex values form a3
	"vmovups      32(%7,%0,8), %%ymm11              \n\t" // 2 complex values form a3

	"vfmadd231pd      %%ymm8 , %%ymm4, %%ymm12      \n\t" // a_r[0] * x_r , a_i[0] * x_r, a_r[1] * x_r, a_i[1] * x_r
	"vfmadd231pd      %%ymm8 , %%ymm5, %%ymm13      \n\t" // a_r[0] * x_i , a_i[0] * x_i, a_r[1] * x_i, a_i[1] * x_i
	"vfmadd231pd      %%ymm9 , %%ymm4, %%ymm14      \n\t" // a_r[2] * x_r , a_i[2] * x_r, a_r[3] * x_r, a_i[3] * x_r
	"vfmadd231pd      %%ymm9 , %%ymm5, %%ymm15      \n\t" // a_r[2] * x_i , a_i[2] * x_i, a_r[3] * x_i, a_i[3] * x_i

	"vfmadd231pd      %%ymm10, %%ymm6, %%ymm12      \n\t" // a_r[0] * x_r , a_i[0] * x_r, a_r[1] * x_r, a_i[1] * x_r
	"vfmadd231pd      %%ymm10, %%ymm7, %%ymm13      \n\t" // a_r[0] * x_i , a_i[0] * x_i, a_r[1] * x_i, a_i[1] * x_i
	"vfmadd231pd      %%ymm11, %%ymm6, %%ymm14      \n\t" // a_r[2] * x_r , a_i[2] * x_r, a_r[3] * x_r, a_i[3] * x_r
	"vfmadd231pd      %%ymm11, %%ymm7, %%ymm15      \n\t" // a_r[2] * x_i , a_i[2] * x_i, a_r[3] * x_i, a_i[3] * x_i


#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
        "vpermilpd      $0x5 , %%ymm13, %%ymm13               \n\t"
        "vpermilpd      $0x5 , %%ymm15, %%ymm15               \n\t"
        "vaddsubpd      %%ymm13, %%ymm12, %%ymm8              \n\t"
        "vaddsubpd      %%ymm15, %%ymm14, %%ymm9              \n\t"
#else
        "vpermilpd      $0x5 , %%ymm12, %%ymm12               \n\t"
        "vpermilpd      $0x5 , %%ymm14, %%ymm14               \n\t"
        "vaddsubpd      %%ymm13, %%ymm12, %%ymm8              \n\t"
        "vaddsubpd      %%ymm15, %%ymm14, %%ymm9              \n\t"
        "vpermilpd      $0x5 , %%ymm8 , %%ymm8                \n\t"
        "vpermilpd      $0x5 , %%ymm9 , %%ymm9                \n\t"
#endif

	"prefetcht0      192(%3,%0,8)			\n\t"
	"vmovups	  (%3,%0,8),  %%ymm12           \n\t"
	"vmovups	32(%3,%0,8),  %%ymm13           \n\t"

#if !defined(XCONJ)
        "vaddpd         %%ymm8, %%ymm12, %%ymm12              \n\t"
        "vaddpd         %%ymm9, %%ymm13, %%ymm13              \n\t"
#else
        "vaddsubpd              %%ymm12, %%ymm8, %%ymm12              \n\t"
        "vaddsubpd              %%ymm13, %%ymm9, %%ymm13              \n\t"
#endif


	"vmovups  %%ymm12,   (%3,%0,8)		        \n\t" // 2 complex values to y	
	"vmovups  %%ymm13, 32(%3,%0,8)		        \n\t"	

        "addq		$8 , %0	  	 	        \n\t"
	"subq	        $4 , %1			        \n\t"		
	"jnz		.L01LOOP%=		        \n\t"
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
	  "%xmm0", "%xmm1", "%xmm2", "%xmm3", 
	  "%xmm4", "%xmm5", "%xmm6", "%xmm7", 
	  "%xmm8", "%xmm9", "%xmm10", "%xmm11", 
	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
	  "memory"
	);

} 


