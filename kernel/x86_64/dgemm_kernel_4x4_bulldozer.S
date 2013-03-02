/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/*   1. Redistributions of source code must retain the above         */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer.                                                  */
/*                                                                   */
/*   2. Redistributions in binary form must reproduce the above      */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer in the documentation and/or other materials       */
/*      provided with the distribution.                              */
/*                                                                   */
/*    THIS  SOFTWARE IS PROVIDED  BY THE  UNIVERSITY OF  TEXAS AT    */
/*    AUSTIN  ``AS IS''  AND ANY  EXPRESS OR  IMPLIED WARRANTIES,    */
/*    INCLUDING, BUT  NOT LIMITED  TO, THE IMPLIED  WARRANTIES OF    */
/*    MERCHANTABILITY  AND FITNESS FOR  A PARTICULAR  PURPOSE ARE    */
/*    DISCLAIMED.  IN  NO EVENT SHALL THE UNIVERSITY  OF TEXAS AT    */
/*    AUSTIN OR CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT,    */
/*    INCIDENTAL,  SPECIAL, EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES    */
/*    (INCLUDING, BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE    */
/*    GOODS  OR  SERVICES; LOSS  OF  USE,  DATA,  OR PROFITS;  OR    */
/*    BUSINESS INTERRUPTION) HOWEVER CAUSED  AND ON ANY THEORY OF    */
/*    LIABILITY, WHETHER  IN CONTRACT, STRICT  LIABILITY, OR TORT    */
/*    (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY WAY OUT    */
/*    OF  THE  USE OF  THIS  SOFTWARE,  EVEN  IF ADVISED  OF  THE    */
/*    POSSIBILITY OF SUCH DAMAGE.                                    */
/*                                                                   */
/* The views and conclusions contained in the software and           */
/* documentation are those of the authors and should not be          */
/* interpreted as representing official policies, either expressed   */
/* or implied, of The University of Texas at Austin.                 */
/*********************************************************************/

#define ASSEMBLER
#include "common.h"
 
#define OLD_M	%rdi
#define OLD_N	%rsi
#define M	%r13
#define N	%r14
#define K	%rdx

#define A	%rcx
#define B	%r8
#define C	%r9
#define LDC	%r10
	
#define I	%r11
#define AO	%rdi
#define BO	%rsi
#define	CO1	%r15
#define CO2	%r12
#define BB	%rbp
#define	J	%rbx

#ifndef WINDOWS_ABI

#define STACKSIZE 96

#define ALPHA	 48(%rsp)
#define OFFSET	 56(%rsp)
#define KK	 64(%rsp)
#define KKK	 72(%rsp)

#else

#define STACKSIZE 256

#define OLD_A		40 + STACKSIZE(%rsp)
#define OLD_B		48 + STACKSIZE(%rsp)
#define OLD_C		56 + STACKSIZE(%rsp)
#define OLD_LDC		64 + STACKSIZE(%rsp)
#define OLD_OFFSET	72 + STACKSIZE(%rsp)

#define ALPHA	224(%rsp)
#define OFFSET	232(%rsp)
#define KK	240(%rsp)
#define KKK	248(%rsp)

#endif

#define movapd	movaps
#define movupd	movups

#define KERNEL1(xx) \
	vfmaddpd  %xmm8,%xmm1,%xmm0,%xmm8 ;\
	vmovaps	%xmm2, %xmm0 ;\
	vmovups -14 * SIZE(AO, %rax, 4),%xmm2 ;\
	vfmaddpd %xmm9,%xmm3,%xmm0,%xmm9 ;\
	vfmaddpd %xmm12,%xmm2,%xmm1,%xmm12 ;\
	vmovddup	-14 * SIZE(BO, %rax, 4), %xmm1 ;\
	vfmaddpd %xmm13,%xmm2,%xmm3,%xmm13 ;\
	vmovddup	-13 * SIZE(BO, %rax, 4), %xmm3 ;\
	vfmaddpd %xmm10,%xmm1,%xmm0,%xmm10 ;\
	vfmaddpd %xmm11,%xmm3,%xmm0,%xmm11 ;\
	vfmaddpd %xmm15,%xmm2, %xmm3,%xmm15 ;\
	vfmaddpd %xmm14,%xmm2, %xmm1,%xmm14 ;\
	vmovups	-12 * SIZE(AO, %rax, 4), %xmm0 ;\
	vmovddup	-12 * SIZE(BO, %rax, 4), %xmm1 ;\
 	vmovddup	-11 * SIZE(BO, %rax, 4), %xmm3 ;\
	vmovaps	%xmm0, %xmm2

#define KERNEL2(xx) \
	vfmaddpd %xmm8,%xmm1, %xmm0,%xmm8 ;\
	vmovaps	%xmm2, %xmm0 ;\
	vmovups -10 * SIZE(AO, %rax, 4),%xmm2 ;\
/*A*/	vmovups	  (AO, %rax, 4), %xmm6 ;\
	vfmaddpd %xmm12,%xmm2, %xmm1,%xmm12 ;\
	vfmaddpd %xmm9,%xmm3, %xmm0,%xmm9 ;\
	vmovddup	-10 * SIZE(BO, %rax, 4), %xmm1 ;\
	vfmaddpd %xmm13,%xmm2, %xmm3,%xmm13 ;\
	vmovddup	 -9 * SIZE(BO, %rax, 4), %xmm3 ;\
	vfmaddpd %xmm10,%xmm1, %xmm0,%xmm10 ;\
	vfmaddpd %xmm14,%xmm2, %xmm1,%xmm14 ;\
/**/	vmovddup	  (BO, %rax, 4), %xmm1 ;\
	vfmaddpd %xmm11,%xmm3, %xmm0,%xmm11 ;\
	vfmaddpd %xmm15,%xmm2, %xmm3,%xmm15 ;\
 	vmovddup	 -7 * SIZE(BO, %rax, 4), %xmm3 ;\
	vmovaps	%xmm4, %xmm2

#define KERNEL3(xx) \
	vfmaddpd %xmm8,%xmm5, %xmm4, %xmm8 ;\
	vmovaps	%xmm2, %xmm4 ;\
	vmovups -6 * SIZE(AO, %rax, 4),%xmm2 ;\
	vfmaddpd %xmm12,%xmm2, %xmm5,%xmm12 ;\
	vfmaddpd %xmm9,%xmm3, %xmm4,%xmm9 ;\
	vmovddup	 -6 * SIZE(BO, %rax, 4), %xmm5 ;\
	vfmaddpd %xmm13,%xmm2, %xmm3,%xmm13 ;\
	vmovddup	 -5 * SIZE(BO, %rax, 4), %xmm3 ;\
	vfmaddpd %xmm10,%xmm5, %xmm4,%xmm10 ;\
	vfmaddpd %xmm14,%xmm2, %xmm5,%xmm14 ;\
	vfmaddpd %xmm11,%xmm3, %xmm4, %xmm11 ;\
	vmovups	 -4 * SIZE(AO, %rax, 4), %xmm4 ;\
	vmovddup	 -4 * SIZE(BO, %rax, 4), %xmm5 ;\
	vfmaddpd %xmm15,%xmm2,%xmm3,%xmm15 ;\
 	vmovddup	 -3 * SIZE(BO, %rax, 4), %xmm3 ;\
	vmovaps	%xmm4, %xmm2

#define KERNEL4(xx) \
	vfmaddpd %xmm8,%xmm5, %xmm4,%xmm8 ;\
	vmovaps	%xmm2, %xmm4 ;\
	vmovups -2 * SIZE(AO, %rax, 4),%xmm2 ;\
	vfmaddpd %xmm12,%xmm2, %xmm5 ,%xmm12;\
/*A*/	vmovups	  8 * SIZE(AO, %rax, 4), %xmm7 ;\
	vmovddup	 -2 * SIZE(BO, %rax, 4), %xmm5 ;\
	vfmaddpd %xmm9,%xmm3, %xmm4,%xmm9 ;\
	vfmaddpd %xmm13,%xmm2, %xmm3,%xmm13 ;\
	vmovddup	 -1 * SIZE(BO, %rax, 4), %xmm3 ;\
	vfmaddpd %xmm10,%xmm5, %xmm4,%xmm10 ;\
	vfmaddpd %xmm14,%xmm2, %xmm5,%xmm14 ;\
/**/	vmovddup	  8 * SIZE(BO, %rax, 4), %xmm5 ;\
	vfmaddpd %xmm11,%xmm3, %xmm4,%xmm11 ;\
	vfmaddpd %xmm15,%xmm2, %xmm3,%xmm15 ;\
 	vmovddup	  1 * SIZE(BO, %rax, 4), %xmm3 ;\
	vmovaps	%xmm6, %xmm2

#define KERNEL5(xx) \
	vfmaddpd %xmm8,%xmm1, %xmm6,%xmm8 ;\
	vmovaps	%xmm2, %xmm6 ;\
	vmovups 2 * SIZE(AO, %rax, 4),%xmm2 ;\
	vfmaddpd %xmm12,%xmm2, %xmm1,%xmm12 ;\
	vmovddup	  2 * SIZE(BO, %rax, 4), %xmm1 ;\
	vfmaddpd %xmm9,%xmm3, %xmm6,%xmm9 ;\
	vfmaddpd %xmm13,%xmm2, %xmm3,%xmm13 ;\
	vmovddup	  3 * SIZE(BO, %rax, 4), %xmm3 ;\
	vfmaddpd %xmm10,%xmm1, %xmm6,%xmm10 ;\
	vfmaddpd %xmm14,%xmm2, %xmm1,%xmm14 ;\
	vfmaddpd %xmm11,%xmm3, %xmm6,%xmm11 ;\
	vmovups	  4 * SIZE(AO, %rax, 4), %xmm6 ;\
	vmovddup	  4 * SIZE(BO, %rax, 4), %xmm1 ;\
	vfmaddpd %xmm15,%xmm2, %xmm3,%xmm15 ;\
 	vmovddup	  5 * SIZE(BO, %rax, 4), %xmm3 ;\
	vmovaps	%xmm6, %xmm2

#define KERNEL6(xx) \
	vfmaddpd %xmm8,%xmm1, %xmm6,%xmm8 ;\
	vmovaps	%xmm2, %xmm6 ;\
	vmovups	6 * SIZE(AO, %rax, 4),%xmm2 ;\
	vfmaddpd %xmm12,%xmm2, %xmm1,%xmm12 ;\
/*A*/	vmovups	 16 * SIZE(AO, %rax, 4), %xmm0 ;\
	vmovddup	  6 * SIZE(BO, %rax, 4), %xmm1 ;\
	vfmaddpd %xmm9,%xmm3, %xmm6,%xmm9 ;\
	vfmaddpd %xmm13,%xmm2, %xmm3,%xmm13 ;\
	vmovddup	  7 * SIZE(BO, %rax, 4), %xmm3 ;\
	vfmaddpd %xmm10,%xmm1, %xmm6,%xmm10 ;\
	vfmaddpd %xmm14,%xmm2, %xmm1,%xmm14 ;\
/**/	vmovddup	 16 * SIZE(BO, %rax, 4), %xmm1 ;\
	vfmaddpd %xmm11,%xmm3, %xmm6,%xmm11 ;\
	vfmaddpd %xmm15,%xmm2,%xmm3,%xmm15 ;\
 	vmovddup	  9 * SIZE(BO, %rax, 4), %xmm3 ;\
	vmovaps	%xmm7, %xmm2

#define KERNEL7(xx) \
	vfmaddpd %xmm8,%xmm5, %xmm7,%xmm8 ;\
	vmovaps	%xmm2, %xmm7 ;\
	vmovups 10 * SIZE(AO, %rax, 4),%xmm2 ;\
	vfmaddpd %xmm12,%xmm2, %xmm5,%xmm12 ;\
	vmovddup	 10 * SIZE(BO, %rax, 4), %xmm5 ;\
	vfmaddpd %xmm9,%xmm3, %xmm7,%xmm9 ;\
	vfmaddpd %xmm13,%xmm2, %xmm3,%xmm13 ;\
	vmovddup	 11 * SIZE(BO, %rax, 4), %xmm3 ;\
	vfmaddpd %xmm10,%xmm5, %xmm7,%xmm10 ;\
	vfmaddpd %xmm14,%xmm2, %xmm5,%xmm14 ;\
	vfmaddpd %xmm11,%xmm3, %xmm7,%xmm11 ;\
	vmovups	 12 * SIZE(AO, %rax, 4), %xmm7 ;\
	vmovddup	 12 * SIZE(BO, %rax, 4), %xmm5 ;\
	vfmaddpd %xmm15,%xmm2, %xmm3,%xmm15 ;\
 	vmovddup	 13 * SIZE(BO, %rax, 4), %xmm3 ;\
	vmovaps	%xmm7, %xmm2

#define KERNEL8(xx) \
	vfmaddpd %xmm8,%xmm5, %xmm7,%xmm8 ;\
	vmovaps	%xmm2, %xmm7 ;\
	vmovups 14 * SIZE(AO, %rax, 4),%xmm2 ;\
/*A*/	vmovups	 24 * SIZE(AO, %rax, 4), %xmm4 ;\
	vfmaddpd %xmm12,%xmm2, %xmm5,%xmm12 ;\
	vfmaddpd %xmm9,%xmm3, %xmm7,%xmm9 ;\
	vmovddup	 14 * SIZE(BO, %rax, 4), %xmm5 ;\
	vfmaddpd %xmm13, %xmm2, %xmm3,%xmm13 ;\
	vmovddup	 15 * SIZE(BO, %rax, 4), %xmm3 ;\
	vfmaddpd %xmm10,%xmm5, %xmm7,%xmm10 ;\
	vfmaddpd %xmm14,%xmm2, %xmm5,%xmm14 ;\
/**/	vmovddup	 24 * SIZE(BO, %rax, 4), %xmm5 ;\
	vfmaddpd %xmm11,%xmm3, %xmm7,%xmm11 ;\
	vfmaddpd %xmm15,%xmm2, %xmm3,%xmm15 ;\
 	vmovddup	 17 * SIZE(BO, %rax, 4), %xmm3 ;\
	vmovaps	%xmm0, %xmm2 ;\
	addq	$8 * SIZE, %rax ;\

#define KERNEL_SUB1(xx) \
        vfmaddpd   %xmm8, %xmm1, %xmm0,%xmm8 ;\
        vmovapd  %xmm2, %xmm0 ;\
	vmovups  -14 * SIZE(AO),%xmm2 ;\
        vfmaddpd   %xmm12, %xmm2, %xmm1,%xmm12 ;\
        vmovddup -14 * SIZE(BO), %xmm1 ;\
        vfmaddpd   %xmm9, %xmm3, %xmm0,%xmm9 ;\
        vfmaddpd   %xmm13, %xmm2, %xmm3,%xmm13 ;\
        vmovddup -13 * SIZE(BO), %xmm3 ;\
        vfmaddpd   %xmm10, %xmm1, %xmm0,%xmm10 ;\
        vfmaddpd   %xmm14,%xmm2, %xmm1 ,%xmm14 ;\
        vfmaddpd   %xmm11, %xmm3, %xmm0,%xmm11 ;\
        vfmaddpd   %xmm15, %xmm2, %xmm3,%xmm15 ;\
        vmovups  -12 * SIZE(AO), %xmm0 ;\
        vmovddup -12 * SIZE(BO), %xmm1 ;\
        vmovddup -11 * SIZE(BO), %xmm3 ;\
        vmovapd  %xmm0, %xmm2


#define KERNEL_SUB2(xx) \
	vfmaddpd %xmm8,%xmm1, %xmm0,%xmm8 ;\
	vmovaps	%xmm2, %xmm0 ;\
	vmovups -10 * SIZE(AO),%xmm2 ;\
	vfmaddpd %xmm12,%xmm2, %xmm1,%xmm12 ;\
	vfmaddpd %xmm9,%xmm3, %xmm0,%xmm9 ;\
	vfmaddpd %xmm13,%xmm2, %xmm3,%xmm13 ;\
	vmovddup	-10 * SIZE(BO), %xmm1 ;\
	vmovddup	 -9 * SIZE(BO), %xmm3 ;\
	vfmaddpd %xmm10,%xmm1, %xmm0,%xmm10 ;\
	vfmaddpd %xmm14,%xmm2, %xmm1,%xmm14 ;\
	vfmaddpd %xmm11,%xmm3, %xmm0,%xmm11 ;\
	vfmaddpd %xmm15,%xmm2, %xmm3,%xmm15 ;\
	vmovups	  (AO), %xmm0 ;\
	vmovddup	  (BO), %xmm1 ;\
 	vmovddup	 -7 * SIZE(BO), %xmm3 ;\
	vmovaps	%xmm4, %xmm2

#define KERNEL_SUB3(xx) \
	vfmaddpd %xmm8,%xmm5, %xmm4,%xmm8 ;\
	vmovaps	%xmm2, %xmm4 ;\
	vmovups -6 * SIZE(AO),%xmm2 ;\
	vfmaddpd %xmm12,%xmm2, %xmm5,%xmm12 ;\
	vfmaddpd %xmm9,%xmm3, %xmm4,%xmm9 ;\
	vmovddup	 -6 * SIZE(BO), %xmm5 ;\
	vfmaddpd %xmm13,%xmm2, %xmm3,%xmm13 ;\
	vmovddup	 -5 * SIZE(BO), %xmm3 ;\
	vfmaddpd %xmm10,%xmm5, %xmm4,%xmm10 ;\
	vfmaddpd %xmm14,%xmm2, %xmm5,%xmm14 ;\
	vfmaddpd %xmm11,%xmm3, %xmm4,%xmm11 ;\
	vmovups	 -4 * SIZE(AO), %xmm4 ;\
	vmovddup	 -4 * SIZE(BO), %xmm5 ;\
	vfmaddpd %xmm15,%xmm2, %xmm3,%xmm15 ;\
 	vmovddup	 -3 * SIZE(BO), %xmm3 ;\
	vmovaps	%xmm4, %xmm2

#define KERNEL_SUB4(xx) \
	vfmaddpd %xmm8,%xmm5, %xmm4,%xmm8 ;\
	vmovaps	%xmm2, %xmm4 ;\
	vmovups -2 * SIZE(AO),%xmm2 ;\
	vfmaddpd %xmm12,%xmm2, %xmm5,%xmm12 ;\
	vfmaddpd %xmm13,%xmm2, %xmm3,%xmm13 ;\
	vmovddup	 -2 * SIZE(BO), %xmm5 ;\
	vfmaddpd %xmm9,%xmm3, %xmm4,%xmm9 ;\
	vmovddup	 -1 * SIZE(BO), %xmm3 ;\
	vfmaddpd %xmm10,%xmm5, %xmm4,%xmm10 ;\
	vfmaddpd %xmm14,%xmm2, %xmm5,%xmm14 ;\
	vfmaddpd %xmm11,%xmm3, %xmm4,%xmm11 ;\
	vfmaddpd %xmm15,%xmm2, %xmm3,%xmm15 ;\
 	vmovddup	  1 * SIZE(BO), %xmm3 ;\
	vmovaps	%xmm0, %xmm2

	PROLOGUE
	PROFCODE
	
	subq	$STACKSIZE, %rsp
	movq	%rbx,   (%rsp)
	movq	%rbp,  8(%rsp)
	movq	%r12, 16(%rsp)
	movq	%r13, 24(%rsp)
	movq	%r14, 32(%rsp)
	movq	%r15, 40(%rsp)

	vzeroupper

#ifdef WINDOWS_ABI
	movq	%rdi,    48(%rsp)
	movq	%rsi,    56(%rsp)
	movups	%xmm6,   64(%rsp)
	movups	%xmm7,   80(%rsp)
	movups	%xmm8,   96(%rsp)
	movups	%xmm9,  112(%rsp)
	movups	%xmm10, 128(%rsp)
	movups	%xmm11, 144(%rsp)
	movups	%xmm12, 160(%rsp)
	movups	%xmm13, 176(%rsp)
	movups	%xmm14, 192(%rsp)
	movups	%xmm15, 208(%rsp)

	movq	ARG1,      OLD_M
	movq	ARG2,      OLD_N
	movq	ARG3,      K
	movq	OLD_A,     A
	movq	OLD_B,     B
	movq	OLD_C,     C
	movq	OLD_LDC,   LDC
#ifdef TRMMKERNEL
	movsd	OLD_OFFSET, %xmm12
#endif
	vmovaps	%xmm3, %xmm0

#else
	movq	STACKSIZE +  8(%rsp), LDC
#ifdef TRMMKERNEL
	movsd	STACKSIZE + 16(%rsp), %xmm12
#endif

#endif

	movq	OLD_M, M
	movq	OLD_N, N

	subq	$-16 * SIZE, A
	subq	$-16 * SIZE, B

	vmovsd	 %xmm0, ALPHA

	salq	$BASE_SHIFT, LDC

#ifdef TRMMKERNEL
	vmovsd	%xmm12, OFFSET
	vmovsd	%xmm12, KK
#ifndef LEFT
	negq	KK
#endif	
#endif
	movq	N,  J
	sarq	$2, J		# j = (n >> 2)
	jle	.L40
	ALIGN_4

.L01:
	movq	C, CO1			# coffset1 = c
	leaq	(C, LDC, 2), CO2	# coffset2 = c + ldc
	
	leaq	(C, LDC, 4), C		# c += 4 * ldc

#if defined(TRMMKERNEL) && defined(LEFT)
	movq	OFFSET, %rax
	movq	%rax, KK
#endif	

	movq	A, AO		# aoffset = a

	movq	K, %rax
	salq	$BASE_SHIFT + 2, %rax
	leaq	 (B, %rax), BB

	movq	M,  I
	sarq	$2, I	# i = (m >> 2)
	jle	.L20
	ALIGN_4

	.align 16
.L11:
#if !defined(TRMMKERNEL) || \
	(defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
	(defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
	movq	B, BO
#else
	movq	KK, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 4), AO
	leaq	(B,  %rax, 4), BO
#endif	

	vzeroall
	prefetcht0 256(CO1)
	prefetcht0 320(CO1)
	prefetcht0 256(CO2)
	prefetcht0 320(CO2)
	vmovups	-16 * SIZE(AO), %xmm0
	vmovddup	-16 * SIZE(BO), %xmm1
 	vmovddup	-15 * SIZE(BO), %xmm3
	vmovups	 -8 * SIZE(AO), %xmm4
	vmovddup	 -8 * SIZE(BO), %xmm5

	vmovaps	%xmm0, %xmm2


#ifndef TRMMKERNEL
	movq	K, %rax
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	movq	K, %rax
	subq	KK, %rax
	movq	%rax, KKK	
#else
	movq	KK, %rax
#ifdef LEFT
	addq	$4, %rax
#else
	addq	$4, %rax
#endif
	movq	%rax, KKK
#endif

	andq	$-8, %rax
	salq	$BASE_SHIFT, %rax
	leaq	(AO, %rax, 4), AO
	leaq	(BO, %rax, 4), BO
	negq	%rax
	NOBRANCH
	je	.L15
	ALIGN_4

	.align 16
.L12:
	prefetcht0 (AO,%rax,4)
	prefetcht0 (BO,%rax,4)
	KERNEL1(16 *  0)
	KERNEL2(16 *  0)
	KERNEL3(16 *  0)
	KERNEL4(16 *  0)
	KERNEL5(16 *  0)
	KERNEL6(16 *  0)
	KERNEL7(16 *  0)
	KERNEL8(16 *  0)
	NOBRANCH
	je	.L15
	KERNEL1(16 *  0)
	KERNEL2(16 *  0)
	KERNEL3(16 *  0)
	KERNEL4(16 *  0)
	KERNEL5(16 *  0)
	KERNEL6(16 *  0)
	KERNEL7(16 *  0)
	KERNEL8(16 *  0)
	NOBRANCH
	je	.L15
	KERNEL1(16 *  0)
	KERNEL2(16 *  0)
	KERNEL3(16 *  0)
	KERNEL4(16 *  0)
	KERNEL5(16 *  0)
	KERNEL6(16 *  0)
	KERNEL7(16 *  0)
	KERNEL8(16 *  0)
	NOBRANCH
	je	.L15
	KERNEL1(16 *  0)
	KERNEL2(16 *  0)
	KERNEL3(16 *  0)
	KERNEL4(16 *  0)
	KERNEL5(16 *  0)
	KERNEL6(16 *  0)
	KERNEL7(16 *  0)
	KERNEL8(16 *  0)
	NOBRANCH
	je	.L15
	KERNEL1(16 *  0)
	KERNEL2(16 *  0)
	KERNEL3(16 *  0)
	KERNEL4(16 *  0)
	KERNEL5(16 *  0)
	KERNEL6(16 *  0)
	KERNEL7(16 *  0)
	KERNEL8(16 *  0)
	NOBRANCH
	je	.L15
	KERNEL1(16 *  0)
	KERNEL2(16 *  0)
	KERNEL3(16 *  0)
	KERNEL4(16 *  0)
	KERNEL5(16 *  0)
	KERNEL6(16 *  0)
	KERNEL7(16 *  0)
	KERNEL8(16 *  0)
	NOBRANCH
	je	.L15
	KERNEL1(16 *  0)
	KERNEL2(16 *  0)
	KERNEL3(16 *  0)
	KERNEL4(16 *  0)
	KERNEL5(16 *  0)
	KERNEL6(16 *  0)
	KERNEL7(16 *  0)
	KERNEL8(16 *  0)
	NOBRANCH
	je	.L15
	KERNEL1(16 *  0)
	KERNEL2(16 *  0)
	KERNEL3(16 *  0)
	KERNEL4(16 *  0)
	KERNEL5(16 *  0)
	KERNEL6(16 *  0)
	KERNEL7(16 *  0)
	KERNEL8(16 *  0)
	jl	.L12
	ALIGN_4

.L15:
	vmovddup	ALPHA, %xmm7

#ifndef TRMMKERNEL
	movq	K, %rax
#else
	movq	KKK, %rax
#endif
	testq	$4, %rax
	je .L16
	ALIGN_4

	KERNEL_SUB1(16 *  0)
	KERNEL_SUB2(16 *  0)
	KERNEL_SUB3(16 *  0)
	KERNEL_SUB4(16 *  0)

	subq	$-16 * SIZE, BO
	subq	$-16 * SIZE, AO
	ALIGN_4

.L16:
#ifndef TRMMKERNEL
	movq	K, %rax
#else
	movq	KKK, %rax
#endif
	andq	$3, %rax		# if (k & 1)
	je .L19

	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 4), AO
	leaq	(BO, %rax, 4), BO
	negq	%rax
	ALIGN_4

.L17:
	vfmaddpd %xmm8,%xmm1, %xmm0,%xmm8
	vmovaps	%xmm2, %xmm0
	vmovups -14 * SIZE(AO, %rax, 4),%xmm2
	vfmaddpd %xmm12,%xmm2, %xmm1,%xmm12
	vmovddup	-14 * SIZE(BO, %rax, 4), %xmm1
	vfmaddpd %xmm9,%xmm3, %xmm0,%xmm9
	vfmaddpd %xmm13,%xmm2, %xmm3,%xmm13
	vmovddup	-13 * SIZE(BO, %rax, 4), %xmm3
	vfmaddpd %xmm10,%xmm1, %xmm0,%xmm10
	vfmaddpd %xmm14,%xmm2, %xmm1,%xmm14
	vfmaddpd %xmm11,%xmm3, %xmm0,%xmm11
	vmovups	-12 * SIZE(AO, %rax, 4), %xmm0
	vmovddup	-12 * SIZE(BO, %rax, 4), %xmm1
	vfmaddpd %xmm15,%xmm2, %xmm3,%xmm15
 	vmovddup	-11 * SIZE(BO, %rax, 4), %xmm3
	vmovaps	%xmm0, %xmm2

	addq	$SIZE, %rax
	jl	.L17
	ALIGN_4

.L19:
	// prefetch	   -8 * SIZE(BB)
	subq		 $-16 * SIZE, BB

#ifndef TRMMKERNEL

	vfmaddpd (CO1),%xmm7, %xmm8,%xmm8
	vfmaddpd 2 * SIZE(CO1),%xmm7, %xmm12,%xmm12
	.align 2
	vfmaddpd (CO1, LDC),%xmm7, %xmm9,%xmm9
	vfmaddpd 2 * SIZE(CO1, LDC),%xmm7, %xmm13,%xmm13
	.align 2
	vfmaddpd (CO2),%xmm7, %xmm10,%xmm10
	vfmaddpd 2 * SIZE(CO2),%xmm7, %xmm14,%xmm14
	.align 2
	vfmaddpd (CO2, LDC),%xmm7, %xmm11,%xmm11
	vfmaddpd 2 * SIZE(CO2, LDC),%xmm7, %xmm15,%xmm15

#else
	vmulpd	%xmm7, %xmm8,%xmm8
	vmulpd	%xmm7, %xmm12,%xmm12
	vmulpd	%xmm7, %xmm9,%xmm9
	vmulpd	%xmm7, %xmm13,%xmm13
	vmulpd	%xmm7, %xmm10,%xmm10
	vmulpd	%xmm7, %xmm14,%xmm14
	vmulpd	%xmm7, %xmm11,%xmm11
	vmulpd	%xmm7, %xmm15,%xmm15

#endif

	.align 2
	vmovups	%xmm8,  (CO1)
	vmovups	%xmm12, 2 * SIZE(CO1)
	.align 2
	vmovups	%xmm9,  (CO1, LDC)
	vmovups	%xmm13, 2 * SIZE(CO1, LDC)
	.align 2
	vmovups	%xmm10, (CO2)
	vmovups	%xmm14, 2 * SIZE(CO2)
	.align 2
	vmovups	%xmm11, (CO2, LDC)
	vmovups	%xmm15, 2 * SIZE(CO2, LDC)

#if (defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
    (defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
	movq	K, %rax
	subq	KKK, %rax
	leaq	(,%rax, SIZE), %rax
	leaq	(AO, %rax, 4), AO
	leaq	(BO, %rax, 4), BO
#endif

#if defined(TRMMKERNEL) && defined(LEFT)
	addq	$4, KK
#endif

	addq	$4 * SIZE, CO1		# coffset += 4
	addq	$4 * SIZE, CO2		# coffset += 4
	decq	I			# i --
	BRANCH
	jg	.L11
	ALIGN_4	

.L20:
	testq	$3, M
	je	.L39

	testq	$2, M
	je	.L30
	ALIGN_4

.L21:
#if !defined(TRMMKERNEL) || \
	(defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
	(defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
	movq	B, BO
#else
	movq	KK, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 2), AO
	leaq	(B,  %rax, 4), BO
#endif	

	vmovups	-16 * SIZE(AO), %xmm0
	vxorps	%xmm8, %xmm8, %xmm8
	vmovups	-12 * SIZE(AO), %xmm2
	vxorps	%xmm9, %xmm9 ,%xmm9
	vmovddup	-16 * SIZE(BO), %xmm1
	vxorps	%xmm10, %xmm10, %xmm10
	vmovddup	-15 * SIZE(BO), %xmm5
	vxorps	%xmm11, %xmm11, %xmm11
	vmovddup	 -8 * SIZE(BO), %xmm3

#ifndef TRMMKERNEL
	movq	K, %rax
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	movq	K, %rax
	subq	KK, %rax
	movq	%rax, KKK	
#else
	movq	KK, %rax
#ifdef LEFT
	addq	$2, %rax
#else
	addq	$4, %rax
#endif
	movq	%rax, KKK
#endif

	andq	$-4, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 2), AO
	leaq	(BO, %rax, 4), BO
	negq	%rax
	NOBRANCH
	je	.L26
	ALIGN_4

.L22:
	vfmaddpd %xmm8,%xmm0, %xmm1,%xmm8
	vfmaddpd %xmm9,%xmm0, %xmm5,%xmm9
	vmovddup	-14 * SIZE(BO, %rax, 4), %xmm1
	vfmaddpd %xmm10,%xmm0, %xmm1,%xmm10
	vmovddup	-13 * SIZE(BO, %rax, 4), %xmm5
	vfmaddpd %xmm11,%xmm0, %xmm5,%xmm11
	vmovups	-14 * SIZE(AO, %rax, 2), %xmm0
	vmovddup	-12 * SIZE(BO, %rax, 4), %xmm1
	vmovddup	-11 * SIZE(BO, %rax, 4), %xmm5
	vfmaddpd %xmm8,%xmm0, %xmm1,%xmm8
	vfmaddpd %xmm9,%xmm0, %xmm5,%xmm9
	vmovddup	-10 * SIZE(BO, %rax, 4), %xmm1
	vmovddup	 -9 * SIZE(BO, %rax, 4), %xmm5
	vfmaddpd %xmm10,%xmm0, %xmm1,%xmm10
	vfmaddpd %xmm11,%xmm0, %xmm5,%xmm11
	vmovddup	  (BO, %rax, 4), %xmm1
	vmovddup	 -7 * SIZE(BO, %rax, 4), %xmm5
	vmovups	 -8 * SIZE(AO, %rax, 2), %xmm0
	vfmaddpd %xmm8,%xmm2, %xmm3,%xmm8
	vfmaddpd %xmm9,%xmm2, %xmm5,%xmm9
	vmovddup	 -6 * SIZE(BO, %rax, 4), %xmm3
	vmovddup	 -5 * SIZE(BO, %rax, 4), %xmm5
	vfmaddpd %xmm10,%xmm2, %xmm3,%xmm10
	vfmaddpd %xmm11,%xmm2, %xmm5,%xmm11
	vmovups	-10 * SIZE(AO, %rax, 2), %xmm2
	vmovddup	 -4 * SIZE(BO, %rax, 4), %xmm3
	vmovddup	 -3 * SIZE(BO, %rax, 4), %xmm5
	vfmaddpd %xmm8,%xmm2, %xmm3,%xmm8
	vfmaddpd %xmm9,%xmm2, %xmm5,%xmm9
	vmovddup	 -2 * SIZE(BO, %rax, 4), %xmm3
	vmovddup	 -1 * SIZE(BO, %rax, 4), %xmm5
	vfmaddpd %xmm10,%xmm2, %xmm3,%xmm10
	vfmaddpd %xmm11,%xmm2, %xmm5,%xmm11
	vmovddup	  8 * SIZE(BO, %rax, 4), %xmm3
	vmovups	 -4 * SIZE(AO, %rax, 2), %xmm2
	vmovddup	  1 * SIZE(BO, %rax, 4), %xmm5

	addq	$4 * SIZE, %rax
	BRANCH
	jl	.L22
	ALIGN_4

.L26:
	vmovddup	ALPHA, %xmm7

#ifndef TRMMKERNEL
	movq	K, %rax
#else
	movq	KKK, %rax
#endif
	andq	$3, %rax		# if (k & 1)
	je .L29

	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 2), AO
	leaq	(BO, %rax, 4), BO
	negq	%rax
	ALIGN_4

.L27:
	vfmaddpd %xmm8,%xmm0, %xmm1,%xmm8
	vmovddup	-14 * SIZE(BO, %rax, 4), %xmm1
	vfmaddpd %xmm9,%xmm0, %xmm5,%xmm9
	vmovddup	-13 * SIZE(BO, %rax, 4), %xmm5
	vfmaddpd %xmm10,%xmm0, %xmm1,%xmm10
	vmovddup	-12 * SIZE(BO, %rax, 4), %xmm1
	vfmaddpd %xmm11,%xmm0, %xmm5,%xmm11
	vmovups	-14 * SIZE(AO, %rax, 2), %xmm0
	vmovddup	-11 * SIZE(BO, %rax, 4), %xmm5

	addq	$SIZE, %rax
	jl	.L27
	ALIGN_4

.L29:
#ifndef TRMMKERNEL

	vfmaddpd (CO1),%xmm7, %xmm8,%xmm8
	vfmaddpd (CO1, LDC),%xmm7, %xmm9,%xmm9
	vfmaddpd (CO2),%xmm7, %xmm10,%xmm10
	vfmaddpd (CO2, LDC),%xmm7, %xmm11,%xmm11

#else
	vmulpd	%xmm7, %xmm8,%xmm8
	vmulpd	%xmm7, %xmm9,%xmm9
	vmulpd	%xmm7, %xmm10,%xmm10
	vmulpd	%xmm7, %xmm11,%xmm11

#endif

	vmovups	%xmm8,  (CO1)
	vmovups	%xmm9,  (CO1, LDC)

	vmovups	%xmm10, (CO2)
	vmovups	%xmm11, (CO2, LDC)

#if (defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
    (defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
	movq	K, %rax
	subq	KKK, %rax
	leaq	(,%rax, SIZE), %rax
	leaq	(AO, %rax, 2), AO
	leaq	(BO, %rax, 4), BO
#endif

#if defined(TRMMKERNEL) && defined(LEFT)
	addq	$2, KK
#endif

	addq	$2 * SIZE, CO1
	addq	$2 * SIZE, CO2
	ALIGN_4	

.L30:
	testq	$1, M
	je	.L39

#if !defined(TRMMKERNEL) || \
	(defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
	(defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
	movq	B, BO
#else
	movq	KK, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 1), AO
	leaq	(B,  %rax, 4), BO
#endif	

	vmovddup	-16 * SIZE(AO), %xmm0
	vxorps	%xmm8, %xmm8, %xmm8
	vmovddup	-14 * SIZE(AO), %xmm2
	vxorps	%xmm9, %xmm9, %xmm9
	vmovddup	-15 * SIZE(AO), %xmm4
	vxorps	%xmm10, %xmm10,%xmm10
	vmovups	-16 * SIZE(BO), %xmm1
	vxorps	%xmm11, %xmm11,%xmm11
	vmovups	 -8 * SIZE(BO), %xmm3

#ifndef TRMMKERNEL
	movq	K, %rax
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	movq	K, %rax
	subq	KK, %rax
	movq	%rax, KKK	
#else
	movq	KK, %rax
#ifdef LEFT
	addq	$1, %rax
#else
	addq	$4, %rax
#endif
	movq	%rax, KKK
#endif

	andq	$-4, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 1), AO
	leaq	(BO, %rax, 4), BO
	negq	%rax
	NOBRANCH
	je	.L36
	ALIGN_4

.L32:
	vfmaddpd %xmm8,%xmm0, %xmm1,%xmm8
	vfmaddpd %xmm9,-14 * SIZE(BO, %rax, 4), %xmm0,%xmm9
	vmovups	-12 * SIZE(BO, %rax, 4), %xmm1
	vmovddup -12 * SIZE(AO, %rax, 1), %xmm0
	vfmaddpd %xmm10,%xmm4, %xmm1,%xmm10
	vfmaddpd %xmm11,-10 * SIZE(BO, %rax, 4), %xmm4,%xmm11
	vmovups	  (BO, %rax, 4), %xmm1
	vmovddup	-11 * SIZE(AO, %rax, 1), %xmm4
	vfmaddpd %xmm8,%xmm2, %xmm3,%xmm8
	vfmaddpd %xmm9,-6 * SIZE(BO, %rax, 4), %xmm2,%xmm9
	vmovups	 -4 * SIZE(BO, %rax, 4), %xmm3
	vmovddup	-13 * SIZE(AO, %rax, 1), %xmm2
	vfmaddpd %xmm10,%xmm2, %xmm3,%xmm10
	vfmaddpd %xmm11,-2 * SIZE(BO, %rax, 4), %xmm2,%xmm11
	vmovups	  8 * SIZE(BO, %rax, 4), %xmm3
	vmovddup	-10 * SIZE(AO, %rax, 1), %xmm2

	addq	$4 * SIZE, %rax
	BRANCH
	jl	.L32
	ALIGN_4

.L36:
	vmovddup	ALPHA, %xmm7

#ifndef TRMMKERNEL
	movq	K, %rax
#else
	movq	KKK, %rax
#endif
	andq	$3, %rax		# if (k & 1)
	je .L38

	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 1), AO
	leaq	(BO, %rax, 4), BO
	negq	%rax
	ALIGN_4

.L37:
	vfmaddpd %xmm8,%xmm0, %xmm1,%xmm8
	vfmaddpd %xmm9,-14 * SIZE(BO, %rax, 4), %xmm0,%xmm9
	vmovups	-12 * SIZE(BO, %rax, 4), %xmm1
	vmovddup -15 * SIZE(AO, %rax, 1), %xmm0

	addq	$SIZE, %rax
	jl	.L37
	ALIGN_4

.L38:
	vaddpd	%xmm10, %xmm8,%xmm8
	vaddpd	%xmm11, %xmm9,%xmm9

#ifndef TRMMKERNEL
	vmovsd	(CO1), %xmm0
	vmovhpd	(CO1, LDC), %xmm0,%xmm0
	vmovsd	(CO2), %xmm1
	vmovhpd	(CO2, LDC), %xmm1,%xmm1


	vfmaddpd	%xmm0, %xmm7,%xmm8,%xmm8
	vfmaddpd	%xmm1, %xmm7,%xmm9,%xmm9
#else

	vmulpd	%xmm7, %xmm8,%xmm8
	vmulpd	%xmm7, %xmm9,%xmm9

#endif

	vmovsd	%xmm8,  (CO1)
	vmovhpd	%xmm8,  (CO1, LDC)
	vmovsd	%xmm9,  (CO2)
	vmovhpd	%xmm9,  (CO2, LDC)

#if (defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
    (defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
	movq	K, %rax
	subq	KKK, %rax
	leaq	(,%rax, SIZE), %rax
	leaq	(AO, %rax, 1), AO
	leaq	(BO, %rax, 4), BO
#endif

#if defined(TRMMKERNEL) && defined(LEFT)
	addq	$1, KK
#endif
	ALIGN_4	
	
.L39:
#if defined(TRMMKERNEL) && !defined(LEFT)
	addq	$4, KK
#endif
	
	movq	BO, B

	decq	J			# j --
	jg	.L01
	ALIGN_4

.L40:
	testq	$3, N
	je	.L999

	testq	$2, N
	je	.L80
	ALIGN_4

.L41:
#if defined(TRMMKERNEL) && defined(LEFT)
	movq	OFFSET, %rax
	movq	%rax, KK
#endif	

	movq	C, CO1			# coffset1 = c
	leaq	(C, LDC, 1), CO2	# coffset2 = c + ldc
	movq	A, AO		# aoffset = a

	movq	K, %rax
	salq	$BASE_SHIFT + 1, %rax
	leaq	(B, %rax), BB

	movq	M,  I
	sarq	$2, I	# i = (m >> 2)
	jle	.L60
	ALIGN_4

.L51:
#if !defined(TRMMKERNEL) || \
	(defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
	(defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
	movq	B, BO
#else
	movq	KK, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 4), AO
	leaq	(B,  %rax, 2), BO
#endif	

	vmovddup	-16 * SIZE(BO), %xmm1
	vmovddup	-15 * SIZE(BO), %xmm5
	vmovddup	-12 * SIZE(BO), %xmm3
	vxorps	%xmm8, %xmm8,%xmm8
	vxorps	%xmm9, %xmm9,%xmm9
	vxorps	%xmm12, %xmm12,%xmm12
	vxorps	%xmm13, %xmm13,%xmm13
	vmovups	-16 * SIZE(AO), %xmm0
	vmovups	 -8 * SIZE(AO), %xmm4
	// prefetcht0 256(CO1)
	// prefetcht0 320(CO1)
	// prefetcht0 256(CO2)
	// prefetcht0 320(CO2)
	// prefetchnta	  24 * SIZE(CO1)
	// prefetchnta	  32 * SIZE(CO1)
	// prefetchw      3 * SIZE(CO1)
	vmovups	%xmm0, %xmm2
	// prefetchw      3 * SIZE(CO2)
	// prefetchnta -16 * SIZE(BB)
	// prefetch	 -16 * SIZE(BB)
	subq		 $-8 * SIZE, BB

#ifndef TRMMKERNEL
	movq	K, %rax
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	movq	K, %rax
	subq	KK, %rax
	movq	%rax, KKK	
#else
	movq	KK, %rax
#ifdef LEFT
	addq	$4, %rax
#else
	addq	$2, %rax
#endif
	movq	%rax, KKK
#endif


	andq	$-4, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 4), AO
	leaq	(BO, %rax, 2), BO
	negq	%rax
	NOBRANCH
	je	.L56
	ALIGN_4

.L52:
	vfmaddpd %xmm8,%xmm1, %xmm0,%xmm8
	vfmaddpd %xmm9,%xmm5, %xmm2,%xmm9
	vmovups  -14 * SIZE(AO, %rax, 4),%xmm2
	vfmaddpd %xmm12,%xmm2, %xmm1,%xmm12
	vmovups	-12 * SIZE(AO, %rax, 4), %xmm0
	vmovddup	-14 * SIZE(BO, %rax, 2), %xmm1
	vfmaddpd %xmm13,%xmm2, %xmm5,%xmm13
	vmovddup	-13 * SIZE(BO, %rax, 2), %xmm5
	vmovups 	-10 * SIZE(AO, %rax, 4), %xmm2
	vfmaddpd %xmm8,%xmm1, %xmm0,%xmm8
	vfmaddpd %xmm12,%xmm2, %xmm1,%xmm12
	vfmaddpd %xmm9,%xmm5, %xmm0,%xmm9
	vmovups	  (AO, %rax, 4), %xmm0
	vmovddup	 -8 * SIZE(BO, %rax, 2), %xmm1
	vfmaddpd %xmm13,%xmm2, %xmm5,%xmm13
	vmovddup	-11 * SIZE(BO, %rax, 2), %xmm5
	vmovups	-6 * SIZE(AO, %rax, 4), %xmm2
	vfmaddpd %xmm8,%xmm3, %xmm4,%xmm8
	vfmaddpd %xmm12,%xmm2, %xmm3,%xmm12
	vfmaddpd %xmm9,%xmm5, %xmm4,%xmm9
	vmovups	 -4 * SIZE(AO, %rax, 4), %xmm4
	vmovddup	-10 * SIZE(BO, %rax, 2), %xmm3
	vfmaddpd %xmm13,%xmm2, %xmm5,%xmm13
	vmovddup	 -9 * SIZE(BO, %rax, 2), %xmm5
	vmovups		-2 * SIZE(AO, %rax, 4), %xmm2
	vfmaddpd %xmm8,%xmm3, %xmm4,%xmm8
	vfmaddpd %xmm12,%xmm2, %xmm3,%xmm12
	vfmaddpd %xmm9,%xmm5, %xmm4,%xmm9
	vfmaddpd %xmm13,%xmm2, %xmm5,%xmm13
	vmovups	  8 * SIZE(AO, %rax, 4), %xmm4
	vmovddup	 -4 * SIZE(BO, %rax, 2), %xmm3
	vmovddup	 -7 * SIZE(BO, %rax, 2), %xmm5
	vmovaps	%xmm0, %xmm2

	addq	$4 * SIZE, %rax
	BRANCH
	jl	.L52
	ALIGN_4

.L56:
	vmovddup	ALPHA, %xmm7

#ifndef TRMMKERNEL
	movq	K, %rax
#else
	movq	KKK, %rax
#endif
	andq	$3, %rax		# if (k & 1)
	je .L59

	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 4), AO
	leaq	(BO, %rax, 2), BO
	negq	%rax
	ALIGN_4

.L57:
	vfmaddpd %xmm8,%xmm1, %xmm0,%xmm8
	vfmaddpd %xmm9,%xmm5, %xmm2,%xmm9
	vmovups -14 * SIZE(AO, %rax, 4),%xmm2
	vfmaddpd %xmm12,%xmm2, %xmm1,%xmm12
	vmovups	-12 * SIZE(AO, %rax, 4), %xmm0
	vmovddup	-14 * SIZE(BO, %rax, 2), %xmm1
	vfmaddpd %xmm13,%xmm2, %xmm5,%xmm13
	vmovddup	-13 * SIZE(BO, %rax, 2), %xmm5
	vmovaps	%xmm0, %xmm2

	addq	$SIZE, %rax
	jl	.L57
	ALIGN_4

.L59:
#ifndef TRMMKERNEL
	vfmaddpd	(CO1),%xmm7, %xmm8, %xmm8
	vfmaddpd	2 * SIZE(CO1),%xmm7, %xmm12, %xmm12
	vfmaddpd	(CO2),%xmm7, %xmm9, %xmm9
	vfmaddpd	2 * SIZE(CO2),%xmm7, %xmm13, %xmm13

#else
	vmulpd	%xmm7, %xmm8,%xmm8
	vmulpd	%xmm7, %xmm9,%xmm9
	vmulpd	%xmm7, %xmm12,%xmm12
	vmulpd	%xmm7, %xmm13,%xmm13

#endif

	vmovups	%xmm8,  (CO1)
	vmovups	%xmm12, 2 * SIZE(CO1)

	vmovups	%xmm9,  (CO2)
	vmovups	%xmm13, 2 * SIZE(CO2)

#if (defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
    (defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
	movq	K, %rax
	subq	KKK, %rax
	leaq	(,%rax, SIZE), %rax
	leaq	(AO, %rax, 4), AO
	leaq	(BO, %rax, 2), BO
#endif

#if defined(TRMMKERNEL) && defined(LEFT)
	addq	$4, KK
#endif

	addq	$4 * SIZE, CO1		# coffset += 4
	addq	$4 * SIZE, CO2		# coffset += 4
	decq	I			# i --
	jg	.L51
	ALIGN_4	

.L60:
	testq	$2, M
	je	.L70
	ALIGN_4

.L61:
#if !defined(TRMMKERNEL) || \
	(defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
	(defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
	movq	B, BO
#else
	movq	KK, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 2), AO
	leaq	(B,  %rax, 2), BO
#endif	

	vmovups	-16 * SIZE(AO), %xmm0
	vxorps	%xmm8, %xmm8,%xmm8
	vmovups	-12 * SIZE(AO), %xmm2
	vxorps	%xmm9, %xmm9,%xmm9
	vmovddup	-16 * SIZE(BO), %xmm1
	vxorps	%xmm10, %xmm10,%xmm10
	vmovddup	-15 * SIZE(BO), %xmm3
	vxorps	%xmm11, %xmm11,%xmm11

#ifndef TRMMKERNEL
	movq	K, %rax
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	movq	K, %rax
	subq	KK, %rax
	movq	%rax, KKK	
#else
	movq	KK, %rax
#ifdef LEFT
	addq	$2, %rax
#else
	addq	$2, %rax
#endif
	movq	%rax, KKK
#endif

	andq	$-4, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 2), AO
	leaq	(BO, %rax, 2), BO
	negq	%rax
	NOBRANCH
	je	.L66
	ALIGN_4

.L62:
	vfmaddpd %xmm8,%xmm0, %xmm1,%xmm8
	vmovddup	-14 * SIZE(BO, %rax, 2), %xmm1
	vfmaddpd %xmm9,%xmm0, %xmm3,%xmm9
	vmovups	-14 * SIZE(AO, %rax, 2), %xmm0
	vmovddup	-13 * SIZE(BO, %rax, 2), %xmm3
	vfmaddpd %xmm10,%xmm0, %xmm1,%xmm10
	vmovddup	-12 * SIZE(BO, %rax, 2), %xmm1
	vfmaddpd %xmm11,%xmm0, %xmm3,%xmm11
	vmovups	 -8 * SIZE(AO, %rax, 2), %xmm0
	vmovddup	-11 * SIZE(BO, %rax, 2), %xmm3
	vfmaddpd %xmm8,%xmm2, %xmm1,%xmm8
	vmovddup	-10 * SIZE(BO, %rax, 2), %xmm1
	vfmaddpd %xmm9,%xmm2, %xmm3,%xmm9
	vmovups	-10 * SIZE(AO, %rax, 2), %xmm2
	vmovddup	 -9 * SIZE(BO, %rax, 2), %xmm3
	vfmaddpd %xmm10,%xmm2, %xmm1,%xmm10
	vmovddup	 -8 * SIZE(BO, %rax, 2), %xmm1
	vfmaddpd %xmm11,%xmm2, %xmm3,%xmm11
	vmovups	 -4 * SIZE(AO, %rax, 2), %xmm2
	vmovddup	 -7 * SIZE(BO, %rax, 2), %xmm3

	addq	$4 * SIZE, %rax
	BRANCH
	jl	.L62
	ALIGN_4

.L66:
	vmovddup	ALPHA, %xmm7

#ifndef TRMMKERNEL
	movq	K, %rax
#else
	movq	KKK, %rax
#endif
	andq	$3, %rax		# if (k & 1)
	je .L69

	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 2), AO
	leaq	(BO, %rax, 2), BO
	negq	%rax
	ALIGN_4

.L67:
	vfmaddpd %xmm8,%xmm0, %xmm1,%xmm8
	vmovddup	-14 * SIZE(BO, %rax, 2), %xmm1
	vfmaddpd %xmm9,%xmm0, %xmm3,%xmm9
	vmovups	-14 * SIZE(AO, %rax, 2), %xmm0
	vmovddup	-13 * SIZE(BO, %rax, 2), %xmm3

	addq	$SIZE, %rax
	jl	.L67
	ALIGN_4

.L69:
	vaddpd	%xmm10, %xmm8,%xmm8
	vaddpd	%xmm11, %xmm9,%xmm9
	
#ifndef TRMMKERNEL

	vfmaddpd	(CO1),%xmm7, %xmm8,%xmm8
	vfmaddpd	(CO2),%xmm7, %xmm9,%xmm9

#else

	vmulpd	%xmm7, %xmm8,%xmm8
	vmulpd	%xmm7, %xmm9,%xmm9

#endif

	vmovups	%xmm8,  (CO1)
	vmovups	%xmm9,  (CO2)

#if (defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
    (defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
	movq	K, %rax
	subq	KKK, %rax
	leaq	(,%rax, SIZE), %rax
	leaq	(AO, %rax, 2), AO
	leaq	(BO, %rax, 2), BO
#endif

#if defined(TRMMKERNEL) && defined(LEFT)
	addq	$2, KK
#endif

	addq	$2 * SIZE, CO1		# coffset += 4
	addq	$2 * SIZE, CO2		# coffset += 4
	ALIGN_4	

.L70:
	testq	$1, M
	je	.L79
	ALIGN_4

.L71:
#if !defined(TRMMKERNEL) || \
	(defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
	(defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
	movq	B, BO
#else
	movq	KK, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 1), AO
	leaq	(B,  %rax, 2), BO
#endif	

	vmovddup	-16 * SIZE(AO), %xmm0
	vxorps	%xmm8, %xmm8,%xmm8
	vmovddup	-15 * SIZE(AO), %xmm1
	vxorps	%xmm9, %xmm9,%xmm9
	vmovddup	-14 * SIZE(AO), %xmm2
	vxorps	%xmm10, %xmm10,%xmm10
	vmovddup	-13 * SIZE(AO), %xmm3
	vxorps	%xmm11, %xmm11,%xmm11

#ifndef TRMMKERNEL
	movq	K, %rax
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	movq	K, %rax
	subq	KK, %rax
	movq	%rax, KKK	
#else
	movq	KK, %rax
#ifdef LEFT
	addq	$1, %rax
#else
	addq	$2, %rax
#endif
	movq	%rax, KKK
#endif

	andq	$-4, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 1), AO
	leaq	(BO, %rax, 2), BO
	negq	%rax
	NOBRANCH
	je	.L76
	ALIGN_4

.L72:
	vfmaddpd %xmm8,-16 * SIZE(BO, %rax, 2), %xmm0,%xmm8
	vmovddup	-12 * SIZE(AO, %rax, 1), %xmm0

	vfmaddpd %xmm9,-14 * SIZE(BO, %rax, 2), %xmm1,%xmm9
	vmovddup	-11 * SIZE(AO, %rax, 1), %xmm1

	vfmaddpd %xmm10,-12 * SIZE(BO, %rax, 2), %xmm2,%xmm10
	vmovddup	-10 * SIZE(AO, %rax, 1), %xmm2

	vfmaddpd %xmm11,-10 * SIZE(BO, %rax, 2), %xmm3,%xmm11
	vmovddup	 -9 * SIZE(AO, %rax, 1), %xmm3

	addq	$4 * SIZE, %rax
	BRANCH
	jl	.L72
	ALIGN_4

.L76:
	vmovddup	ALPHA, %xmm7

#ifndef TRMMKERNEL
	movq	K, %rax
#else
	movq	KKK, %rax
#endif
	andq	$3, %rax		# if (k & 1)
	je .L78

	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 1), AO
	leaq	(BO, %rax, 2), BO
	negq	%rax
	ALIGN_4

.L77:
	vfmaddpd %xmm8,-16 * SIZE(BO, %rax, 2), %xmm0,%xmm8
	vmovddup	-15 * SIZE(AO, %rax, 1), %xmm0

	addq	$SIZE, %rax
	jl	.L77
	ALIGN_4

.L78:
	vaddpd	%xmm9,  %xmm8,%xmm8
	vaddpd	%xmm11, %xmm10,%xmm10
	vaddpd	%xmm10, %xmm8,%xmm8

#ifndef TRMMKERNEL
	vmovsd	(CO1), %xmm0
	vmovhpd	(CO2), %xmm0,%xmm0
#endif

	vmulpd	%xmm7, %xmm8,%xmm8

#ifndef TRMMKERNEL
	vaddpd	%xmm0, %xmm8,%xmm8
#endif

	vmovsd	%xmm8,  (CO1)
	vmovhpd	%xmm8,  (CO2)

#if (defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
    (defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
	movq	K, %rax
	subq	KKK, %rax
	leaq	(,%rax, SIZE), %rax
	leaq	(AO, %rax, 1), AO
	leaq	(BO, %rax, 2), BO
#endif

#if defined(TRMMKERNEL) && defined(LEFT)
	addq	$1, KK
#endif
	ALIGN_4	
	
.L79:
#if defined(TRMMKERNEL) && !defined(LEFT)
	addq	$2, KK
#endif

	movq	BO, B

	leaq	(C, LDC, 2), C
	ALIGN_4

.L80:
	testq	$1, N
	je	.L999
	ALIGN_4

.L81:
#if defined(TRMMKERNEL) && defined(LEFT)
	movq	OFFSET, %rax
	movq	%rax, KK
#endif	

	movq	C, CO1			# coffset1 = c
	movq	A, AO		# aoffset = a

	movq	M,  I
	sarq	$2, I	# i = (m >> 2)
	jle	.L100
	ALIGN_4

.L91:
#if !defined(TRMMKERNEL) || \
	(defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
	(defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
	movq	B, BO
#else
	movq	KK, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 4), AO
	leaq	(B,  %rax, 1), BO
#endif	

	vmovups	 -8 * SIZE(AO), %xmm2
	vxorps	%xmm8, %xmm8,%xmm8
	vmovups	-16 * SIZE(AO), %xmm0
	vxorps	%xmm9, %xmm9,%xmm9
	vmovddup	-16 * SIZE(BO), %xmm1
	vxorps	%xmm12, %xmm12,%xmm12
	vmovddup	-14 * SIZE(BO), %xmm3
	vxorps	%xmm13, %xmm13,%xmm13
	vmovddup	-15 * SIZE(BO), %xmm5

	// prefetchw      3 * SIZE(CO1)

#ifndef TRMMKERNEL
	movq	K, %rax
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	movq	K, %rax
	subq	KK, %rax
	movq	%rax, KKK	
#else
	movq	KK, %rax
#ifdef LEFT
	addq	$4, %rax
#else
	addq	$1, %rax
#endif
	movq	%rax, KKK
#endif

	andq	$-4, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 4), AO
	leaq	(BO, %rax, 1), BO
	negq	%rax
	NOBRANCH
	je	.L96
	ALIGN_4

.L92:
	vfmaddpd %xmm8,%xmm1, %xmm0,%xmm8
	vfmaddpd %xmm12,-14 * SIZE(AO, %rax, 4), %xmm1,%xmm12
	vmovapd	-12 * SIZE(AO, %rax, 4), %xmm0
	vmovddup	-12 * SIZE(BO, %rax, 1), %xmm1
	vfmaddpd %xmm9,%xmm5, %xmm0,%xmm9
	vfmaddpd %xmm13,-10 * SIZE(AO, %rax, 4), %xmm5,%xmm13
	vmovapd	  (AO, %rax, 4), %xmm0
	vmovddup	-13 * SIZE(BO, %rax, 1), %xmm5
	vfmaddpd %xmm8,%xmm3, %xmm2,%xmm8
	vfmaddpd %xmm12,-6 * SIZE(AO, %rax, 4), %xmm3,%xmm12
	vmovapd	 -4 * SIZE(AO, %rax, 4), %xmm2
	vmovddup	-10 * SIZE(BO, %rax, 1), %xmm3
	vfmaddpd %xmm9,%xmm5, %xmm2,%xmm9
	vfmaddpd %xmm13,-2 * SIZE(AO, %rax, 4), %xmm5,%xmm13
	vmovapd	  8 * SIZE(AO, %rax, 4), %xmm2
	vmovddup	-11 * SIZE(BO, %rax, 1), %xmm5

	addq	$4 * SIZE, %rax
	BRANCH
	jl	.L92
	ALIGN_4

.L96:
	vmovddup	ALPHA, %xmm7

#ifndef TRMMKERNEL
	movq	K, %rax
#else
	movq	KKK, %rax
#endif
	andq	$3, %rax		# if (k & 1)
	je .L99

	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 4), AO
	leaq	(BO, %rax, 1), BO
	negq	%rax
	ALIGN_4

.L97:
	vfmaddpd %xmm8,%xmm1, %xmm0,%xmm8
	vfmaddpd %xmm12,-14 * SIZE(AO, %rax, 4), %xmm1,%xmm12
	vmovups	-12 * SIZE(AO, %rax, 4), %xmm0
	vmovddup	-15 * SIZE(BO, %rax, 1), %xmm1

	addq	$SIZE, %rax
	jl	.L97
	ALIGN_4

.L99:
	vaddpd	%xmm9,  %xmm8,%xmm8
	vaddpd	%xmm13, %xmm12,%xmm12

#ifndef TRMMKERNEL

	vfmaddpd		(CO1),%xmm7, %xmm8,%xmm8
	vfmaddpd	2 * SIZE(CO1),%xmm7,%xmm12,%xmm12

#else
	vmulpd	%xmm7, %xmm8,%xmm8
	vmulpd	%xmm7, %xmm12,%xmm12

#endif

	vmovups	%xmm8,  (CO1)
	vmovups	%xmm12, 2 * SIZE(CO1)

#if (defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
    (defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
	movq	K, %rax
	subq	KKK, %rax
	leaq	(,%rax, SIZE), %rax
	leaq	(AO, %rax, 4), AO
	leaq	(BO, %rax, 1), BO
#endif

#if defined(TRMMKERNEL) && defined(LEFT)
	addq	$4, KK
#endif

	addq	$4 * SIZE, CO1		# coffset += 4
	decq	I			# i --
	jg	.L91
	ALIGN_4	

.L100:
	testq	$2, M
	je	.L110
	ALIGN_4

.L101:
#if !defined(TRMMKERNEL) || \
	(defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
	(defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
	movq	B, BO
#else
	movq	KK, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 2), AO
	leaq	(B,  %rax, 1), BO
#endif	

	vmovddup	-16 * SIZE(BO), %xmm0
	vxorps	%xmm8, %xmm8,%xmm8
	vmovddup	-15 * SIZE(BO), %xmm1
	vxorps	%xmm9, %xmm9,%xmm9
	vmovddup	-14 * SIZE(BO), %xmm2
	vxorps	%xmm10, %xmm10,%xmm10
	vmovddup	-13 * SIZE(BO), %xmm3
	vxorps	%xmm11, %xmm11,%xmm11

#ifndef TRMMKERNEL
	movq	K, %rax
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	movq	K, %rax
	subq	KK, %rax
	movq	%rax, KKK	
#else
	movq	KK, %rax
#ifdef LEFT
	addq	$2, %rax
#else
	addq	$1, %rax
#endif
	movq	%rax, KKK
#endif

	andq	$-4, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 2), AO
	leaq	(BO, %rax, 1), BO
	negq	%rax
	NOBRANCH
	je	.L106
	ALIGN_4

.L102:
	vfmaddpd %xmm8,-16 * SIZE(AO, %rax, 2), %xmm0,%xmm8
	vmovddup	-12 * SIZE(BO, %rax, 1), %xmm0

	vfmaddpd %xmm9,-14 * SIZE(AO, %rax, 2), %xmm1,%xmm9
	vmovddup	-11 * SIZE(BO, %rax, 1), %xmm1

	vfmaddpd %xmm10,-12 * SIZE(AO, %rax, 2), %xmm2,%xmm10
	vmovddup	-10 * SIZE(BO, %rax, 1), %xmm2

	vfmaddpd %xmm11,-10 * SIZE(AO, %rax, 2), %xmm3,%xmm11
	vmovddup	 -9 * SIZE(BO, %rax, 1), %xmm3

	addq	$4 * SIZE, %rax
	BRANCH
	jl	.L102
	ALIGN_4

.L106:
	vmovddup	ALPHA, %xmm7

#ifndef TRMMKERNEL
	movq	K, %rax
#else
	movq	KKK, %rax
#endif
	andq	$3, %rax		# if (k & 1)
	je .L109

	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 2), AO
	leaq	(BO, %rax, 1), BO
	negq	%rax
	ALIGN_4

.L107:
	vmovddup	-16 * SIZE(BO, %rax, 1), %xmm0
	vfmaddpd %xmm8,-16 * SIZE(AO, %rax, 2), %xmm0,%xmm8

	addq	$SIZE, %rax
	jl	.L107
	ALIGN_4

.L109:
	vaddpd	%xmm9,  %xmm8,%xmm8
	vaddpd	%xmm11, %xmm10,%xmm10
	vaddpd	%xmm10, %xmm8,%xmm8

#ifndef TRMMKERNEL

	vfmaddpd	(CO1),%xmm7, %xmm8,%xmm8
#else
	vmulpd	%xmm7, %xmm8,%xmm8

#endif

	vmovups	%xmm8,  (CO1)

#if (defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
    (defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
	movq	K, %rax
	subq	KKK, %rax
	leaq	(,%rax, SIZE), %rax
	leaq	(AO, %rax, 2), AO
	leaq	(BO, %rax, 1), BO
#endif

#if defined(TRMMKERNEL) && defined(LEFT)
	addq	$2, KK
#endif

	addq	$2 * SIZE, CO1		# coffset += 4

	ALIGN_4	

.L110:
	testq	$1, M
	je	.L999
	ALIGN_4

.L111:
#if !defined(TRMMKERNEL) || \
	(defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
	(defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))
	movq	B, BO
#else
	movq	KK, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 1), AO
	leaq	(B,  %rax, 1), BO
#endif	

	vmovups	-16 * SIZE(AO), %xmm0
	vxorps	%xmm8, %xmm8,%xmm8
	movups	-14 * SIZE(AO), %xmm1
	vxorps	%xmm9, %xmm9,%xmm9

#ifndef TRMMKERNEL
	movq	K, %rax
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	movq	K, %rax
	subq	KK, %rax
	movq	%rax, KKK	
#else
	movq	KK, %rax
#ifdef LEFT
	addq	$1, %rax
#else
	addq	$1, %rax
#endif
	movq	%rax, KKK
#endif

	andq	$-4, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 1), AO
	leaq	(BO, %rax, 1), BO
	negq	%rax
	NOBRANCH
	je	.L116
	ALIGN_4

.L112:
	vfmaddpd %xmm8,-16 * SIZE(BO, %rax, 1), %xmm0,%xmm8
	vmovups	-12 * SIZE(AO, %rax, 1), %xmm0

	vfmaddpd %xmm9,-14 * SIZE(BO, %rax, 1), %xmm1,%xmm9
	vmovups	-10 * SIZE(AO, %rax, 1), %xmm1

	addq	$4 * SIZE, %rax
	BRANCH
	jl	.L112
	ALIGN_4

.L116:
	vmovddup	ALPHA, %xmm7

#ifndef TRMMKERNEL
	movq	K, %rax
#else
	movq	KKK, %rax
#endif
	andq	$3, %rax		# if (k & 1)
	je .L118

	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 1), AO
	leaq	(BO, %rax, 1), BO
	negq	%rax
	ALIGN_4

.L117:
	vmulsd	-16 * SIZE(BO, %rax, 1), %xmm0,%xmm0
	vaddsd	%xmm0, %xmm8,%xmm8
	vmovsd	-15 * SIZE(AO, %rax, 1), %xmm0

	addq	$SIZE, %rax
	jl	.L117
	ALIGN_4

.L118:
	vaddpd	%xmm9, %xmm8,%xmm8
	vhaddpd	%xmm8, %xmm8,%xmm8

#ifndef TRMMKERNEL
	vmovsd	(CO1), %xmm0
#endif

	vmulsd	%xmm7, %xmm8,%xmm8

#ifndef TRMMKERNEL
	vaddsd	%xmm0, %xmm8,%xmm8
#endif

	vmovsd	%xmm8,  (CO1)
	ALIGN_4
	
.L999:
	movq	   (%rsp), %rbx
	movq	  8(%rsp), %rbp
	movq	 16(%rsp), %r12
	movq	 24(%rsp), %r13
	movq	 32(%rsp), %r14
	movq	 40(%rsp), %r15

#ifdef WINDOWS_ABI
	movq	 48(%rsp), %rdi
	movq	 56(%rsp), %rsi
	movups	 64(%rsp), %xmm6
	movups	 80(%rsp), %xmm7
	movups	 96(%rsp), %xmm8
	movups	112(%rsp), %xmm9
	movups	128(%rsp), %xmm10
	movups	144(%rsp), %xmm11
	movups	160(%rsp), %xmm12
	movups	176(%rsp), %xmm13
	movups	192(%rsp), %xmm14
	movups	208(%rsp), %xmm15
#endif

	addq	$STACKSIZE, %rsp
	ret

	EPILOGUE
