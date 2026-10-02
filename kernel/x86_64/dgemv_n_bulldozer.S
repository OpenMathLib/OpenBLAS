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
#include "l2param.h"

#undef ALIGNED_ACCESS

#define A_PRE 256

#define VMOVUPS_A1(OFF, ADDR, REGS) vmovups	OFF(ADDR), REGS
#define VMOVUPS_A2(OFF, ADDR, BASE, SCALE, REGS) vmovups	OFF(ADDR, BASE, SCALE), REGS
#define VMOVUPS_YL1(OFF, ADDR, REGS) vmovups	OFF(ADDR), REGS
#define VMOVUPS_YS1(OFF, ADDR, REGS) vmovups	REGS,  OFF(ADDR)

#if GEMV_UNROLL < 2
#undef  GEMV_UNROLL
#define GEMV_UNROLL 2
#endif

#ifndef WINDOWS_ABI

#define STACKSIZE	64

#define OLD_M	  %rdi
#define OLD_N	  %rsi
#define OLD_A	  %rcx
#define OLD_LDA	  %r8
#define STACK_INCX	 8 + STACKSIZE(%rsp)
#define STACK_Y		16 + STACKSIZE(%rsp)
#define STACK_INCY	24 + STACKSIZE(%rsp)
#define STACK_BUFFER	32 + STACKSIZE(%rsp)
#define ALPHA	48	      (%rsp)

#else

#define STACKSIZE	256

#define OLD_M	  %rcx
#define OLD_N	  %rdx
#define OLD_A		 40 + STACKSIZE(%rsp)
#define OLD_LDA		 48 + STACKSIZE(%rsp)
#define OLD_X		 56 + STACKSIZE(%rsp)
#define STACK_INCX	 64 + STACKSIZE(%rsp)
#define STACK_Y		 72 + STACKSIZE(%rsp)
#define STACK_INCY	 80 + STACKSIZE(%rsp)
#define STACK_BUFFER	 88 + STACKSIZE(%rsp)
#define ALPHA	224	       (%rsp)

#endif

#define LDA	%r8
#define X	%r9

#define INCX	%rsi
#define INCY	%rdi

#define M	%r10
#define N	%r11
#define A	%r12
#define Y	%r14
#define BUFFER	%r13

#define I	%rax
#define A1	%rbx
#define A2	%rcx
#define LDA3	%rdx
#define Y1	%rbp

#ifdef ALIGNED_ACCESS
#define MM	%r15
#else
#define MM	M
#endif

	PROLOGUE
	PROFCODE

	subq	$STACKSIZE, %rsp
	movq	%rbx,  0(%rsp)
	movq	%rbp,  8(%rsp)
	movq	%r12, 16(%rsp)
	movq	%r13, 24(%rsp)
	movq	%r14, 32(%rsp)
	movq	%r15, 40(%rsp)

#ifdef WINDOWS_ABI
	movq	%rdi,    48(%rsp)
	movq	%rsi,    56(%rsp)
	vmovups	%xmm6,   64(%rsp)
	vmovups	%xmm7,   80(%rsp)
	vmovups	%xmm8,   96(%rsp)
	vmovups	%xmm9,  112(%rsp)
	vmovups	%xmm10, 128(%rsp)
	vmovups	%xmm11, 144(%rsp)
	vmovups	%xmm12, 160(%rsp)
	vmovups	%xmm13, 176(%rsp)
	vmovups	%xmm14, 192(%rsp)
	vmovups	%xmm15, 208(%rsp)

	movq	OLD_M,	      M
	movq	OLD_N,        N
	movq	OLD_A,        A
	movq	OLD_LDA,      LDA
	movq	OLD_X,        X
#else
	movq	OLD_M,	      M
	movq	OLD_N,        N
	movq	OLD_A,        A
	movq	OLD_LDA,      LDA
#endif

	movq	STACK_INCX,   INCX
	movq	STACK_Y,      Y
	movq	STACK_INCY,   INCY
	movq	STACK_BUFFER, BUFFER

#ifndef WINDOWS_ABI
	vmovsd	 %xmm0, ALPHA
#else
	vmovsd	 %xmm3, ALPHA
#endif

	leaq	-1(INCY), %rax

	leaq	(,INCX, SIZE), INCX
	leaq	(,INCY, SIZE), INCY
	leaq	(,LDA,  SIZE), LDA

	leaq	(LDA, LDA, 2), LDA3

	subq	$-16 * SIZE, A

#ifdef ALIGNED_ACCESS
	leaq	-1 (M), MM
	testq	$SIZE, A
	cmoveq	M, MM
#endif

	testq	N, N		# if n <= 0 goto END
	jle	L(999)
	testq	M, M		# if n <= 0 goto END
	jle	L(999)

#if !defined(COPY_FORCE) && !defined(ALIGNED_ACCESS)
#ifndef  NOCOPY_UNALIGNED
	movq	Y,    Y1
	andq	$0xf, Y1
	orq	Y1, %rax
#endif
	testq	%rax, %rax
 	cmoveq	Y, BUFFER
	je	L(10)
#endif

	movq	BUFFER, Y1

	vxorpd	%xmm4, %xmm4, %xmm4

	movq	M,   %rax
	addq	$16, %rax
	sarq	$4,  %rax
	ALIGN_3

L(01):
	vmovups	%xmm4,  0 * SIZE(Y1)
	vmovups	%xmm4,  2 * SIZE(Y1)
	vmovups	%xmm4,  4 * SIZE(Y1)
	vmovups	%xmm4,  6 * SIZE(Y1)
	vmovups	%xmm4,  8 * SIZE(Y1)
	vmovups	%xmm4, 10 * SIZE(Y1)
	vmovups	%xmm4, 12 * SIZE(Y1)
	vmovups	%xmm4, 14 * SIZE(Y1)
	subq	$-16 * SIZE, Y1
	decq	%rax
	jg	L(01)
	ALIGN_3

L(10):

#ifdef ALIGNED_ACCESS
	leaq	SIZE(BUFFER), %rax
	testq	$SIZE, A
	cmovne	%rax, BUFFER

	testq	$SIZE, LDA
	jne	L(50)
#endif

#if GEMV_UNROLL >= 8

	cmpq	$8, N
	jl	L(20)
	ALIGN_3

L(11):
	subq	$8, N

	leaq	16 * SIZE(BUFFER), Y1
	movq	A,  A1
	leaq	(A,  LDA, 4), A2
	leaq	(A,  LDA, 8), A

	vmovddup	(X), %xmm8
	addq	INCX, X
	vmovddup	(X), %xmm9
	addq	INCX, X
	vmovddup	(X), %xmm10
	addq	INCX, X
	vmovddup	(X), %xmm11
	addq	INCX, X
	vmovddup	(X), %xmm12
	addq	INCX, X
	vmovddup	(X), %xmm13
	addq	INCX, X
	vmovddup	(X), %xmm14
	addq	INCX, X
	vmovddup	(X), %xmm15
	addq	INCX, X

	vmovddup	ALPHA, %xmm0

	vmulpd	%xmm0, %xmm8  , %xmm8
	vmulpd	%xmm0, %xmm9  , %xmm9
	vmulpd	%xmm0, %xmm10 , %xmm10
	vmulpd	%xmm0, %xmm11 , %xmm11
	vmulpd	%xmm0, %xmm12 , %xmm12
	vmulpd	%xmm0, %xmm13 , %xmm13
	vmulpd	%xmm0, %xmm14 , %xmm14
	vmulpd	%xmm0, %xmm15 , %xmm15

#ifdef ALIGNED_ACCESS
	testq	$SIZE, A
	je	L(1X)

	vmovsd	 -16 * SIZE(Y1),         %xmm0

	vmovsd	 -16 * SIZE(A1),         %xmm4
	vmovsd	 -16 * SIZE(A1, LDA),    %xmm5
	vmovsd	 -16 * SIZE(A1, LDA, 2), %xmm6
	vmovsd	 -16 * SIZE(A1, LDA3),   %xmm7


	vfmaddsd %xmm0 , %xmm8 , %xmm4 , %xmm0
	vfmaddsd %xmm0 , %xmm9 , %xmm5 , %xmm0
	vfmaddsd %xmm0 , %xmm10, %xmm6 , %xmm0
	vfmaddsd %xmm0 , %xmm11, %xmm7 , %xmm0

	vmovsd	 -16 * SIZE(A2),         %xmm4
	vmovsd	 -16 * SIZE(A2, LDA),    %xmm5
	vmovsd	 -16 * SIZE(A2, LDA, 2), %xmm6
	vmovsd	 -16 * SIZE(A2, LDA3),   %xmm7

	vfmaddsd %xmm0 , %xmm12, %xmm4 , %xmm0
	vfmaddsd %xmm0 , %xmm13, %xmm5 , %xmm0
	vfmaddsd %xmm0 , %xmm14, %xmm6 , %xmm0
	vfmaddsd %xmm0 , %xmm15, %xmm7 , %xmm0

	vmovsd	 %xmm0, -16 * SIZE(Y1)

	addq	 $SIZE, A1
	addq	 $SIZE, A2
	addq	 $SIZE, Y1
	ALIGN_3

L(1X):
#endif

	movq	MM,  I
	sarq	$3,  I
	jle	L(15)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	VMOVUPS_YL1(-14 * SIZE, Y1, %xmm1)
	VMOVUPS_YL1(-12 * SIZE, Y1, %xmm2)
	VMOVUPS_YL1(-10 * SIZE, Y1, %xmm3)

	decq	 I
	jle	 L(14)
	ALIGN_5

L(13):


	prefetchnta	A_PRE(A1)
	vfmaddpd %xmm0 , -16 * SIZE(A1) , %xmm8, %xmm0
	vfmaddpd %xmm1 , -14 * SIZE(A1) , %xmm8, %xmm1
	vfmaddpd %xmm2 , -12 * SIZE(A1) , %xmm8, %xmm2
	vfmaddpd %xmm3 , -10 * SIZE(A1) , %xmm8, %xmm3
	nop

	prefetchnta	A_PRE(A1,LDA,1)
	vfmaddpd %xmm0 , -16 * SIZE(A1, LDA, 1) , %xmm9 , %xmm0
	vfmaddpd %xmm1 , -14 * SIZE(A1, LDA, 1) , %xmm9 , %xmm1
	vfmaddpd %xmm2 , -12 * SIZE(A1, LDA, 1) , %xmm9 , %xmm2
	vfmaddpd %xmm3 , -10 * SIZE(A1, LDA, 1) , %xmm9 , %xmm3

	prefetchnta	A_PRE(A1,LDA,2)
	vfmaddpd %xmm0 , -16 * SIZE(A1, LDA, 2) , %xmm10, %xmm0
	vfmaddpd %xmm1 , -14 * SIZE(A1, LDA, 2) , %xmm10, %xmm1
	vfmaddpd %xmm2 , -12 * SIZE(A1, LDA, 2) , %xmm10, %xmm2
	vfmaddpd %xmm3 , -10 * SIZE(A1, LDA, 2) , %xmm10, %xmm3

	prefetchnta	A_PRE(A1,LDA3,1)
	vfmaddpd %xmm0 , -16 * SIZE(A1, LDA3, 1) , %xmm11, %xmm0
	vfmaddpd %xmm1 , -14 * SIZE(A1, LDA3, 1) , %xmm11, %xmm1
	vfmaddpd %xmm2 , -12 * SIZE(A1, LDA3, 1) , %xmm11, %xmm2
	vfmaddpd %xmm3 , -10 * SIZE(A1, LDA3, 1) , %xmm11, %xmm3

	prefetchnta	A_PRE(A2)
	vfmaddpd %xmm0 , -16 * SIZE(A2) , %xmm12, %xmm0
	vfmaddpd %xmm1 , -14 * SIZE(A2) , %xmm12, %xmm1
	vfmaddpd %xmm2 , -12 * SIZE(A2) , %xmm12, %xmm2
	vfmaddpd %xmm3 , -10 * SIZE(A2) , %xmm12, %xmm3
	nop

	prefetchnta	A_PRE(A2,LDA,1)
	vfmaddpd %xmm0 , -16 * SIZE(A2, LDA, 1) , %xmm13, %xmm0
	vfmaddpd %xmm1 , -14 * SIZE(A2, LDA, 1) , %xmm13, %xmm1
	vfmaddpd %xmm2 , -12 * SIZE(A2, LDA, 1) , %xmm13, %xmm2
	vfmaddpd %xmm3 , -10 * SIZE(A2, LDA, 1) , %xmm13, %xmm3

	prefetchnta	A_PRE(A2,LDA,2)
	vfmaddpd %xmm0 , -16 * SIZE(A2, LDA, 2) , %xmm14, %xmm0
	vfmaddpd %xmm1 , -14 * SIZE(A2, LDA, 2) , %xmm14, %xmm1
	vfmaddpd %xmm2 , -12 * SIZE(A2, LDA, 2) , %xmm14, %xmm2
	vfmaddpd %xmm3 , -10 * SIZE(A2, LDA, 2) , %xmm14, %xmm3

	prefetchnta	A_PRE(A2,LDA3,1)
	vfmaddpd %xmm0 , -16 * SIZE(A2, LDA3, 1) , %xmm15, %xmm0
	vfmaddpd %xmm1 , -14 * SIZE(A2, LDA3, 1) , %xmm15, %xmm1
	vfmaddpd %xmm2 , -12 * SIZE(A2, LDA3, 1) , %xmm15, %xmm2
	vfmaddpd %xmm3 , -10 * SIZE(A2, LDA3, 1) , %xmm15, %xmm3


	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm0)
	VMOVUPS_YS1(-14 * SIZE, Y1, %xmm1)
	VMOVUPS_YS1(-12 * SIZE, Y1, %xmm2)
	VMOVUPS_YS1(-10 * SIZE, Y1, %xmm3)


	VMOVUPS_YL1( -8 * SIZE, Y1, %xmm0)
	VMOVUPS_YL1( -6 * SIZE, Y1, %xmm1)
	prefetchnta	A_PRE(Y1)
	VMOVUPS_YL1( -4 * SIZE, Y1, %xmm2)
	VMOVUPS_YL1( -2 * SIZE, Y1, %xmm3)

	subq	 $-8 * SIZE, A1
	subq	 $-8 * SIZE, A2
	subq	 $-8 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(13)
	ALIGN_3

L(14):

	vfmaddpd %xmm0 , -16 * SIZE(A1) , %xmm8, %xmm0
	vfmaddpd %xmm1 , -14 * SIZE(A1) , %xmm8, %xmm1
	vfmaddpd %xmm2 , -12 * SIZE(A1) , %xmm8, %xmm2
	vfmaddpd %xmm3 , -10 * SIZE(A1) , %xmm8, %xmm3

	vfmaddpd %xmm0 , -16 * SIZE(A1, LDA, 1) , %xmm9 , %xmm0
	vfmaddpd %xmm1 , -14 * SIZE(A1, LDA, 1) , %xmm9 , %xmm1
	vfmaddpd %xmm2 , -12 * SIZE(A1, LDA, 1) , %xmm9 , %xmm2
	vfmaddpd %xmm3 , -10 * SIZE(A1, LDA, 1) , %xmm9 , %xmm3

	vfmaddpd %xmm0 , -16 * SIZE(A1, LDA, 2) , %xmm10, %xmm0
	vfmaddpd %xmm1 , -14 * SIZE(A1, LDA, 2) , %xmm10, %xmm1
	vfmaddpd %xmm2 , -12 * SIZE(A1, LDA, 2) , %xmm10, %xmm2
	vfmaddpd %xmm3 , -10 * SIZE(A1, LDA, 2) , %xmm10, %xmm3

	vfmaddpd %xmm0 , -16 * SIZE(A1, LDA3, 1) , %xmm11, %xmm0
	vfmaddpd %xmm1 , -14 * SIZE(A1, LDA3, 1) , %xmm11, %xmm1
	vfmaddpd %xmm2 , -12 * SIZE(A1, LDA3, 1) , %xmm11, %xmm2
	vfmaddpd %xmm3 , -10 * SIZE(A1, LDA3, 1) , %xmm11, %xmm3

	vfmaddpd %xmm0 , -16 * SIZE(A2) , %xmm12, %xmm0
	vfmaddpd %xmm1 , -14 * SIZE(A2) , %xmm12, %xmm1
	vfmaddpd %xmm2 , -12 * SIZE(A2) , %xmm12, %xmm2
	vfmaddpd %xmm3 , -10 * SIZE(A2) , %xmm12, %xmm3

	vfmaddpd %xmm0 , -16 * SIZE(A2, LDA, 1) , %xmm13, %xmm0
	vfmaddpd %xmm1 , -14 * SIZE(A2, LDA, 1) , %xmm13, %xmm1
	vfmaddpd %xmm2 , -12 * SIZE(A2, LDA, 1) , %xmm13, %xmm2
	vfmaddpd %xmm3 , -10 * SIZE(A2, LDA, 1) , %xmm13, %xmm3

	vfmaddpd %xmm0 , -16 * SIZE(A2, LDA, 2) , %xmm14, %xmm0
	vfmaddpd %xmm1 , -14 * SIZE(A2, LDA, 2) , %xmm14, %xmm1
	vfmaddpd %xmm2 , -12 * SIZE(A2, LDA, 2) , %xmm14, %xmm2
	vfmaddpd %xmm3 , -10 * SIZE(A2, LDA, 2) , %xmm14, %xmm3

	vfmaddpd %xmm0 , -16 * SIZE(A2, LDA3, 1) , %xmm15, %xmm0
	vfmaddpd %xmm1 , -14 * SIZE(A2, LDA3, 1) , %xmm15, %xmm1
	vfmaddpd %xmm2 , -12 * SIZE(A2, LDA3, 1) , %xmm15, %xmm2
	vfmaddpd %xmm3 , -10 * SIZE(A2, LDA3, 1) , %xmm15, %xmm3



	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm0)
	VMOVUPS_YS1(-14 * SIZE, Y1, %xmm1)

	VMOVUPS_YS1(-12 * SIZE, Y1, %xmm2)
	VMOVUPS_YS1(-10 * SIZE, Y1, %xmm3)

	subq	 $-8 * SIZE, A1
	subq	 $-8 * SIZE, A2
	subq	 $-8 * SIZE, Y1
	ALIGN_3

L(15):
	testq	$4, MM
	je	L(16)

	VMOVUPS_A1(-16 * SIZE, A1, %xmm4)
	VMOVUPS_A1(-14 * SIZE, A1, %xmm5)
	VMOVUPS_A2(-16 * SIZE, A1, LDA, 1, %xmm6)
	VMOVUPS_A2(-14 * SIZE, A1, LDA, 1, %xmm7)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	VMOVUPS_YL1(-14 * SIZE, Y1, %xmm1)

	vfmaddpd %xmm0 , %xmm8 , %xmm4 , %xmm0
	vfmaddpd %xmm1 , %xmm8 , %xmm5 , %xmm1
	vfmaddpd %xmm0 , %xmm9 , %xmm6 , %xmm0
	vfmaddpd %xmm1 , %xmm9 , %xmm7 , %xmm1

	VMOVUPS_A2(-16 * SIZE, A1, LDA, 2, %xmm4)
	VMOVUPS_A2(-14 * SIZE, A1, LDA, 2, %xmm5)
	VMOVUPS_A2(-16 * SIZE, A1, LDA3, 1, %xmm6)
	VMOVUPS_A2(-14 * SIZE, A1, LDA3, 1, %xmm7)

	vfmaddpd %xmm0 , %xmm10, %xmm4 , %xmm0
	vfmaddpd %xmm1 , %xmm10, %xmm5 , %xmm1
	vfmaddpd %xmm0 , %xmm11, %xmm6 , %xmm0
	vfmaddpd %xmm1 , %xmm11, %xmm7 , %xmm1

	VMOVUPS_A1(-16 * SIZE, A2, %xmm4)
	VMOVUPS_A1(-14 * SIZE, A2, %xmm5)
	VMOVUPS_A2(-16 * SIZE, A2, LDA, 1, %xmm6)
	VMOVUPS_A2(-14 * SIZE, A2, LDA, 1, %xmm7)

	vfmaddpd %xmm0 , %xmm12, %xmm4 , %xmm0
	vfmaddpd %xmm1 , %xmm12, %xmm5 , %xmm1
	vfmaddpd %xmm0 , %xmm13, %xmm6 , %xmm0
	vfmaddpd %xmm1 , %xmm13, %xmm7 , %xmm1

	VMOVUPS_A2(-16 * SIZE, A2, LDA, 2, %xmm4)
	VMOVUPS_A2(-14 * SIZE, A2, LDA, 2, %xmm5)
	VMOVUPS_A2(-16 * SIZE, A2, LDA3, 1, %xmm6)
	VMOVUPS_A2(-14 * SIZE, A2, LDA3, 1, %xmm7)

	vfmaddpd %xmm0 , %xmm14, %xmm4 , %xmm0
	vfmaddpd %xmm1 , %xmm14, %xmm5 , %xmm1
	vfmaddpd %xmm0 , %xmm15, %xmm6 , %xmm0
	vfmaddpd %xmm1 , %xmm15, %xmm7 , %xmm1


	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm0)
	VMOVUPS_YS1(-14 * SIZE, Y1, %xmm1)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, A2
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(16):
	testq	$2, MM
	je	L(17)

	VMOVUPS_A1(-16 * SIZE, A1, %xmm4)
	VMOVUPS_A2(-16 * SIZE, A1, LDA,  1, %xmm5)
	VMOVUPS_A2(-16 * SIZE, A1, LDA,  2, %xmm6)
	VMOVUPS_A2(-16 * SIZE, A1, LDA3, 1, %xmm7)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm0)

	vfmaddpd %xmm0 , %xmm8 , %xmm4 , %xmm0
	vfmaddpd %xmm0 , %xmm9 , %xmm5 , %xmm0
	vfmaddpd %xmm0 , %xmm10, %xmm6 , %xmm0
	vfmaddpd %xmm0 , %xmm11, %xmm7 , %xmm0

	VMOVUPS_A1(-16 * SIZE, A2, %xmm4)
	VMOVUPS_A2(-16 * SIZE, A2, LDA,  1, %xmm5)
	VMOVUPS_A2(-16 * SIZE, A2, LDA,  2, %xmm6)
	VMOVUPS_A2(-16 * SIZE, A2, LDA3, 1, %xmm7)

	vfmaddpd %xmm0 , %xmm12, %xmm4 , %xmm0
	vfmaddpd %xmm0 , %xmm13, %xmm5 , %xmm0
	vfmaddpd %xmm0 , %xmm14, %xmm6 , %xmm0
	vfmaddpd %xmm0 , %xmm15, %xmm7 , %xmm0

	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm0)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(17):
	testq	$1, MM
	je	L(18)

	vmovsd	 -16 * SIZE(A1),         %xmm4
	vmovsd	 -16 * SIZE(A1, LDA),    %xmm5
	vmovsd	 -16 * SIZE(A1, LDA, 2), %xmm6
	vmovsd	 -16 * SIZE(A1, LDA3),   %xmm7

	vmovsd	 -16 * SIZE(Y1),         %xmm0


	vfmaddsd %xmm0 , %xmm8 , %xmm4 , %xmm0
	vfmaddsd %xmm0 , %xmm9 , %xmm5 , %xmm0
	vfmaddsd %xmm0 , %xmm10, %xmm6 , %xmm0
	vfmaddsd %xmm0 , %xmm11, %xmm7 , %xmm0

	vmovsd	 -16 * SIZE(A2),         %xmm4
	vmovsd	 -16 * SIZE(A2, LDA),    %xmm5
	vmovsd	 -16 * SIZE(A2, LDA, 2), %xmm6
	vmovsd	 -16 * SIZE(A2, LDA3),   %xmm7

	vfmaddsd %xmm0 , %xmm12, %xmm4 , %xmm0
	vfmaddsd %xmm0 , %xmm13, %xmm5 , %xmm0
	vfmaddsd %xmm0 , %xmm14, %xmm6 , %xmm0
	vfmaddsd %xmm0 , %xmm15, %xmm7 , %xmm0


	vmovsd	 %xmm0, -16 * SIZE(Y1)
	ALIGN_3

L(18):
	cmpq	$8, N
	jge	L(11)
	ALIGN_3

L(20):
#endif

#if GEMV_UNROLL >= 4

	cmpq	$4, N
	jl	L(30)

#if GEMV_UNROLL == 4
	ALIGN_3

L(21):
#endif

	subq	$4, N

	leaq	16 * SIZE(BUFFER), Y1
	movq	A,  A1
	leaq	(A,  LDA, 2), A2
	leaq	(A,  LDA, 4), A

	vmovddup	(X), %xmm12
	addq	INCX, X
	vmovddup	(X), %xmm13
	addq	INCX, X
	vmovddup	(X), %xmm14
	addq	INCX, X
	vmovddup	(X), %xmm15
	addq	INCX, X

	vmovddup	ALPHA, %xmm0

	vmulpd	%xmm0, %xmm12 , %xmm12
	vmulpd	%xmm0, %xmm13 , %xmm13
	vmulpd	%xmm0, %xmm14 , %xmm14
	vmulpd	%xmm0, %xmm15 , %xmm15

#ifdef ALIGNED_ACCESS
	testq	$SIZE, A
	je	L(2X)

	vmovsd	 -16 * SIZE(A1),         %xmm4
	vmovsd	 -16 * SIZE(A1, LDA),    %xmm5
	vmovsd	 -16 * SIZE(A2),         %xmm6
	vmovsd	 -16 * SIZE(A2, LDA),    %xmm7

	vmovsd	 -16 * SIZE(Y1),         %xmm0

	vfmaddsd %xmm0 , %xmm12, %xmm4 , %xmm0
	vfmaddsd %xmm0 , %xmm13, %xmm5 , %xmm0
	vfmaddsd %xmm0 , %xmm14, %xmm6 , %xmm0
	vfmaddsd %xmm0 , %xmm15, %xmm7 , %xmm0

	vmovsd	 %xmm0, -16 * SIZE(Y1)

	addq	 $SIZE, A1
	addq	 $SIZE, A2
	addq	 $SIZE, Y1
	ALIGN_3

L(2X):
#endif

	movq	MM,  I
	sarq	$3,  I
	jle	L(25)

	VMOVUPS_A1(-16 * SIZE, A1, %xmm0)
	VMOVUPS_A1(-14 * SIZE, A1, %xmm1)
	VMOVUPS_A1(-12 * SIZE, A1, %xmm2)
	VMOVUPS_A1(-10 * SIZE, A1, %xmm3)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm8)
	VMOVUPS_YL1(-14 * SIZE, Y1, %xmm9)
	VMOVUPS_YL1(-12 * SIZE, Y1, %xmm10)
	VMOVUPS_YL1(-10 * SIZE, Y1, %xmm11)

	VMOVUPS_A2(-16 * SIZE, A1, LDA, 1, %xmm4)
	VMOVUPS_A2(-14 * SIZE, A1, LDA, 1, %xmm5)
	VMOVUPS_A2(-12 * SIZE, A1, LDA, 1, %xmm6)
	VMOVUPS_A2(-10 * SIZE, A1, LDA, 1, %xmm7)

	decq	 I
	jle	 L(24)
	ALIGN_3

L(23):



	vfmaddpd %xmm8 , %xmm12, %xmm0 , %xmm8
	vfmaddpd %xmm9 , %xmm12, %xmm1 , %xmm9
	vfmaddpd %xmm10, %xmm12, %xmm2 , %xmm10
	vfmaddpd %xmm11, %xmm12, %xmm3 , %xmm11

	VMOVUPS_A1(-16 * SIZE, A2, %xmm0)
	VMOVUPS_A1(-14 * SIZE, A2, %xmm1)
	prefetchnta	A_PRE(A2)
	VMOVUPS_A1(-12 * SIZE, A2, %xmm2)
	VMOVUPS_A1(-10 * SIZE, A2, %xmm3)

	vfmaddpd %xmm8 , %xmm13, %xmm4 , %xmm8
	vfmaddpd %xmm9 , %xmm13, %xmm5 , %xmm9
	vfmaddpd %xmm10, %xmm13, %xmm6 , %xmm10
	vfmaddpd %xmm11, %xmm13, %xmm7 , %xmm11

	VMOVUPS_A2(-16 * SIZE, A2, LDA, 1, %xmm4)
	VMOVUPS_A2(-14 * SIZE, A2, LDA, 1, %xmm5)
	prefetchnta	A_PRE(A2, LDA, 1)
	VMOVUPS_A2(-12 * SIZE, A2, LDA, 1, %xmm6)
	VMOVUPS_A2(-10 * SIZE, A2, LDA, 1, %xmm7)

	vfmaddpd %xmm8 , %xmm14, %xmm0 , %xmm8
	vfmaddpd %xmm9 , %xmm14, %xmm1 , %xmm9
	vfmaddpd %xmm10, %xmm14, %xmm2 , %xmm10
	vfmaddpd %xmm11, %xmm14, %xmm3 , %xmm11

	VMOVUPS_A1( -8 * SIZE, A1, %xmm0)
	VMOVUPS_A1( -6 * SIZE, A1, %xmm1)
	prefetchnta	A_PRE(A1)
	VMOVUPS_A1( -4 * SIZE, A1, %xmm2)
	VMOVUPS_A1( -2 * SIZE, A1, %xmm3)

	vfmaddpd %xmm8 , %xmm15, %xmm4 , %xmm8
	vfmaddpd %xmm9 , %xmm15, %xmm5 , %xmm9
	vfmaddpd %xmm10, %xmm15, %xmm6 , %xmm10
	vfmaddpd %xmm11, %xmm15, %xmm7 , %xmm11

	VMOVUPS_A2( -8 * SIZE, A1, LDA, 1, %xmm4)
	VMOVUPS_A2( -6 * SIZE, A1, LDA, 1, %xmm5)
	prefetchnta	A_PRE(A1, LDA, 1)
	VMOVUPS_A2( -4 * SIZE, A1, LDA, 1, %xmm6)
	VMOVUPS_A2( -2 * SIZE, A1, LDA, 1, %xmm7)


	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm8)
	VMOVUPS_YS1(-14 * SIZE, Y1, %xmm9)
	VMOVUPS_YS1(-12 * SIZE, Y1, %xmm10)
	VMOVUPS_YS1(-10 * SIZE, Y1, %xmm11)

	VMOVUPS_YL1( -8 * SIZE, Y1, %xmm8)
	VMOVUPS_YL1( -6 * SIZE, Y1, %xmm9)
	prefetchnta	A_PRE(Y1)
	VMOVUPS_YL1( -4 * SIZE, Y1, %xmm10)
	VMOVUPS_YL1( -2 * SIZE, Y1, %xmm11)

	subq	 $-8 * SIZE, A1
	subq	 $-8 * SIZE, A2
	subq	 $-8 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(23)
	ALIGN_3

L(24):

	vfmaddpd %xmm8 , %xmm12, %xmm0 , %xmm8
	vfmaddpd %xmm9 , %xmm12, %xmm1 , %xmm9
	vfmaddpd %xmm10, %xmm12, %xmm2 , %xmm10
	vfmaddpd %xmm11, %xmm12, %xmm3 , %xmm11

	VMOVUPS_A1(-16 * SIZE, A2, %xmm0)
	VMOVUPS_A1(-14 * SIZE, A2, %xmm1)
	VMOVUPS_A1(-12 * SIZE, A2, %xmm2)
	VMOVUPS_A1(-10 * SIZE, A2, %xmm3)

	vfmaddpd %xmm8 , %xmm13, %xmm4 , %xmm8
	vfmaddpd %xmm9 , %xmm13, %xmm5 , %xmm9
	vfmaddpd %xmm10, %xmm13, %xmm6 , %xmm10
	vfmaddpd %xmm11, %xmm13, %xmm7 , %xmm11

	VMOVUPS_A2(-16 * SIZE, A2, LDA, 1, %xmm4)
	VMOVUPS_A2(-14 * SIZE, A2, LDA, 1, %xmm5)
	VMOVUPS_A2(-12 * SIZE, A2, LDA, 1, %xmm6)
	VMOVUPS_A2(-10 * SIZE, A2, LDA, 1, %xmm7)

	vfmaddpd %xmm8 , %xmm14, %xmm0 , %xmm8
	vfmaddpd %xmm9 , %xmm14, %xmm1 , %xmm9
	vfmaddpd %xmm10, %xmm14, %xmm2 , %xmm10
	vfmaddpd %xmm11, %xmm14, %xmm3 , %xmm11

	vfmaddpd %xmm8 , %xmm15, %xmm4 , %xmm8
	vfmaddpd %xmm9 , %xmm15, %xmm5 , %xmm9
	vfmaddpd %xmm10, %xmm15, %xmm6 , %xmm10
	vfmaddpd %xmm11, %xmm15, %xmm7 , %xmm11

	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm8)
	VMOVUPS_YS1(-14 * SIZE, Y1, %xmm9)
	VMOVUPS_YS1(-12 * SIZE, Y1, %xmm10)
	VMOVUPS_YS1(-10 * SIZE, Y1, %xmm11)

	subq	 $-8 * SIZE, A1
	subq	 $-8 * SIZE, A2
	subq	 $-8 * SIZE, Y1
	ALIGN_3

L(25):
	testq	$4, MM
	je	L(26)

	VMOVUPS_A1(-16 * SIZE, A1, %xmm0)
	VMOVUPS_A1(-14 * SIZE, A1, %xmm1)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm8)
	VMOVUPS_YL1(-14 * SIZE, Y1, %xmm9)

	vfmaddpd %xmm8 , %xmm12, %xmm0 , %xmm8
	vfmaddpd %xmm9 , %xmm12, %xmm1 , %xmm9

	VMOVUPS_A2(-16 * SIZE, A1, LDA, 1, %xmm4)
	VMOVUPS_A2(-14 * SIZE, A1, LDA, 1, %xmm5)

	vfmaddpd %xmm8 , %xmm13, %xmm4 , %xmm8
	vfmaddpd %xmm9 , %xmm13, %xmm5 , %xmm9

	VMOVUPS_A1(-16 * SIZE, A2, %xmm0)
	VMOVUPS_A1(-14 * SIZE, A2, %xmm1)

	vfmaddpd %xmm8 , %xmm14, %xmm0 , %xmm8
	vfmaddpd %xmm9 , %xmm14, %xmm1 , %xmm9

	VMOVUPS_A2(-16 * SIZE, A2, LDA, 1, %xmm4)
	VMOVUPS_A2(-14 * SIZE, A2, LDA, 1, %xmm5)

	vfmaddpd %xmm8 , %xmm15, %xmm4 , %xmm8
	vfmaddpd %xmm9 , %xmm15, %xmm5 , %xmm9

	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm8)
	VMOVUPS_YS1(-14 * SIZE, Y1, %xmm9)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, A2
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(26):
	testq	$2, MM
	je	L(27)

	VMOVUPS_A1(-16 * SIZE, A1, %xmm8)
	VMOVUPS_A2(-16 * SIZE, A1, LDA, 1, %xmm9)
	VMOVUPS_A1(-16 * SIZE, A2, %xmm10)
	VMOVUPS_A2(-16 * SIZE, A2, LDA, 1, %xmm11)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm0)

	vfmaddpd %xmm0 , %xmm12, %xmm8 , %xmm0
	vfmaddpd %xmm0 , %xmm13, %xmm9 , %xmm0
	vfmaddpd %xmm0 , %xmm14, %xmm10, %xmm0
	vfmaddpd %xmm0 , %xmm15, %xmm11, %xmm0

	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm0)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(27):
	testq	$1, MM
#if GEMV_UNROLL == 4
	je	L(28)
#else
	je	L(30)
#endif

	vmovsd	 -16 * SIZE(Y1), %xmm0

	vmovsd	 -16 * SIZE(A1), %xmm8
	vmovsd	 -16 * SIZE(A1, LDA), %xmm9
	vmovsd	 -16 * SIZE(A2), %xmm10
	vmovsd	 -16 * SIZE(A2, LDA), %xmm11

	vfmaddsd %xmm0 , %xmm12, %xmm8 , %xmm0
	vfmaddsd %xmm0 , %xmm13, %xmm9 , %xmm0
	vfmaddsd %xmm0 , %xmm14, %xmm10, %xmm0
	vfmaddsd %xmm0 , %xmm15, %xmm11, %xmm0

	vmovsd	 %xmm0, -16 * SIZE(Y1)
	ALIGN_3

#if GEMV_UNROLL == 4
L(28):
	cmpq	$4, N
	jge	L(21)
	ALIGN_3

#endif

L(30):
#endif

#if GEMV_UNROLL >= 2

	cmpq	$2, N
	jl	L(40)

#if GEMV_UNROLL == 2
	ALIGN_3

L(31):
#endif

	subq	$2, N

	leaq	16 * SIZE(BUFFER), Y1
	movq	A,  A1
	leaq	(A,  LDA), A2
	leaq	(A,  LDA, 2), A

	vmovddup	(X), %xmm12
	addq	INCX, X
	vmovddup	(X), %xmm13
	addq	INCX, X

	vmovddup	ALPHA, %xmm0

	vmulpd	%xmm0, %xmm12 , %xmm12
	vmulpd	%xmm0, %xmm13 , %xmm13

#ifdef ALIGNED_ACCESS
	testq	$SIZE, A
	je	L(3X)

	vmovsd	 -16 * SIZE(A1), %xmm4
	vmovsd	 -16 * SIZE(A2), %xmm5

	vmovsd	 -16 * SIZE(Y1), %xmm0

	vfmaddsd %xmm0 , %xmm12, %xmm4 , %xmm0
	vfmaddsd %xmm0 , %xmm13, %xmm5 , %xmm0

	vmovsd	 %xmm0, -16 * SIZE(Y1)

	addq	 $SIZE, A1
	addq	 $SIZE, A2
	addq	 $SIZE, Y1
	ALIGN_3

L(3X):
#endif

	movq	MM,  I
	sarq	$3,  I
	jle	L(35)

	VMOVUPS_A1(-16 * SIZE, A1, %xmm0)
	VMOVUPS_A1(-14 * SIZE, A1, %xmm1)
	VMOVUPS_A1(-12 * SIZE, A1, %xmm2)
	VMOVUPS_A1(-10 * SIZE, A1, %xmm3)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm8)
	VMOVUPS_YL1(-14 * SIZE, Y1, %xmm9)
	VMOVUPS_YL1(-12 * SIZE, Y1, %xmm10)
	VMOVUPS_YL1(-10 * SIZE, Y1, %xmm11)

	VMOVUPS_A1(-16 * SIZE, A2, %xmm4)
	VMOVUPS_A1(-14 * SIZE, A2, %xmm5)
	VMOVUPS_A1(-12 * SIZE, A2, %xmm6)
	VMOVUPS_A1(-10 * SIZE, A2, %xmm7)

	decq	 I
	jle	 L(34)
	ALIGN_3

L(33):


	vfmaddpd %xmm8 , %xmm12 , %xmm0 , %xmm8
	vfmaddpd %xmm9 , %xmm12 , %xmm1 , %xmm9
	vfmaddpd %xmm10, %xmm12 , %xmm2 , %xmm10
	vfmaddpd %xmm11, %xmm12 , %xmm3 , %xmm11

	vmovups  -8 * SIZE(A1), %xmm0
	vmovups  -6 * SIZE(A1), %xmm1
	prefetchnta	A_PRE(A1)
	vmovups  -4 * SIZE(A1), %xmm2
	vmovups  -2 * SIZE(A1), %xmm3


	vfmaddpd %xmm8 , %xmm13 , %xmm4 , %xmm8
	vfmaddpd %xmm9 , %xmm13 , %xmm5 , %xmm9
	prefetchnta	A_PRE(A2)
	vfmaddpd %xmm10, %xmm13 , %xmm6 , %xmm10
	vfmaddpd %xmm11, %xmm13 , %xmm7 , %xmm11

	vmovups  -8 * SIZE(A2), %xmm4
	vmovups  -6 * SIZE(A2), %xmm5
	vmovups  -4 * SIZE(A2), %xmm6
	vmovups  -2 * SIZE(A2)  , %xmm7


	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm8)
	VMOVUPS_YS1(-14 * SIZE, Y1, %xmm9)
	VMOVUPS_YS1(-12 * SIZE, Y1, %xmm10)
	VMOVUPS_YS1(-10 * SIZE, Y1, %xmm11)

	VMOVUPS_YL1( -8 * SIZE, Y1, %xmm8)
	VMOVUPS_YL1( -6 * SIZE, Y1, %xmm9)
	prefetchnta	A_PRE(Y1)
	VMOVUPS_YL1( -4 * SIZE, Y1, %xmm10)
	VMOVUPS_YL1( -2 * SIZE, Y1, %xmm11)

	subq	 $-8 * SIZE, A1
	subq	 $-8 * SIZE, A2
	subq	 $-8 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(33)
	ALIGN_3

L(34):
	vfmaddpd %xmm8 , %xmm12 , %xmm0 , %xmm8
	vfmaddpd %xmm9 , %xmm12 , %xmm1 , %xmm9
	vfmaddpd %xmm10, %xmm12 , %xmm2 , %xmm10
	vfmaddpd %xmm11, %xmm12 , %xmm3 , %xmm11

	vfmaddpd %xmm8 , %xmm13 , %xmm4 , %xmm8
	vfmaddpd %xmm9 , %xmm13 , %xmm5 , %xmm9
	vfmaddpd %xmm10, %xmm13 , %xmm6 , %xmm10
	vfmaddpd %xmm11, %xmm13 , %xmm7 , %xmm11

	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm8)
	VMOVUPS_YS1(-14 * SIZE, Y1, %xmm9)
	VMOVUPS_YS1(-12 * SIZE, Y1, %xmm10)
	VMOVUPS_YS1(-10 * SIZE, Y1, %xmm11)

	subq	 $-8 * SIZE, A1
	subq	 $-8 * SIZE, A2
	subq	 $-8 * SIZE, Y1
	ALIGN_3

L(35):
	testq	$4, MM
	je	L(36)


	VMOVUPS_A1(-16 * SIZE, A1, %xmm0)
	VMOVUPS_A1(-14 * SIZE, A1, %xmm1)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm8)
	VMOVUPS_YL1(-14 * SIZE, Y1, %xmm9)

	vfmaddpd %xmm8 , %xmm12 , %xmm0 , %xmm8
	vfmaddpd %xmm9 , %xmm12 , %xmm1 , %xmm9

	VMOVUPS_A1(-16 * SIZE, A2, %xmm4)
	VMOVUPS_A1(-14 * SIZE, A2, %xmm5)

	vfmaddpd %xmm8 , %xmm13 , %xmm4 , %xmm8
	vfmaddpd %xmm9 , %xmm13 , %xmm5 , %xmm9

	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm8)
	VMOVUPS_YS1(-14 * SIZE, Y1, %xmm9)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, A2
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(36):
	testq	$2, MM
	je	L(37)

	VMOVUPS_A1(-16 * SIZE, A1, %xmm8)
	VMOVUPS_A1(-16 * SIZE, A2, %xmm9)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm0)

	vfmaddpd %xmm0 , %xmm12 , %xmm8 , %xmm0
	vfmaddpd %xmm0 , %xmm13 , %xmm9 , %xmm0

	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm0)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(37):
	testq	$1, MM
#if GEMV_UNROLL == 2
	je	L(38)
#else
	je	L(40)
#endif

	vmovsd	 -16 * SIZE(Y1), %xmm0

	vmovsd	 -16 * SIZE(A1), %xmm8
	vmovsd	 -16 * SIZE(A2), %xmm9

	vfmaddsd %xmm0 , %xmm12 , %xmm8 , %xmm0
	vfmaddsd %xmm0 , %xmm13 , %xmm9 , %xmm0

	vmovsd	 %xmm0, -16 * SIZE(Y1)
	ALIGN_3

#if GEMV_UNROLL == 2
L(38):
	cmpq	$2, N
	jge	L(31)
	ALIGN_3

#endif

L(40):
	cmpq	$1, N
	jl	L(900)
#endif

	leaq	16 * SIZE(BUFFER), Y1
	movq	A,  A1

	vmovddup	(X), %xmm12
	addq	INCX, X

	vmovddup	ALPHA, %xmm0

	vmulpd	%xmm0, %xmm12 , %xmm12

#ifdef ALIGNED_ACCESS
	testq	$SIZE, A
	je	L(4X)

	vmovsd	 -16 * SIZE(A1), %xmm4
	vmovsd	 -16 * SIZE(Y1), %xmm0

	vfmaddsd %xmm0 , %xmm12 , %xmm4 , %xmm0

	vmovsd	 %xmm0, -16 * SIZE(Y1)

	addq	 $SIZE, A1
	addq	 $SIZE, Y1
	ALIGN_3

L(4X):
#endif

	movq	MM,  I
	sarq	$3,  I
	jle	L(45)

	VMOVUPS_A1(-16 * SIZE, A1, %xmm0)
	VMOVUPS_A1(-14 * SIZE, A1, %xmm1)
	VMOVUPS_A1(-12 * SIZE, A1, %xmm2)
	VMOVUPS_A1(-10 * SIZE, A1, %xmm3)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm8)
	VMOVUPS_YL1(-14 * SIZE, Y1, %xmm9)
	VMOVUPS_YL1(-12 * SIZE, Y1, %xmm10)
	VMOVUPS_YL1(-10 * SIZE, Y1, %xmm11)

	decq	 I
	jle	 L(44)
	ALIGN_3

L(43):


	vfmaddpd %xmm8 , %xmm12 , %xmm0 , %xmm8
	vfmaddpd %xmm9 , %xmm12 , %xmm1 , %xmm9
	vfmaddpd %xmm10, %xmm12 , %xmm2 , %xmm10
	vfmaddpd %xmm11, %xmm12 , %xmm3 , %xmm11

	VMOVUPS_A1( -8 * SIZE, A1, %xmm0)
	VMOVUPS_A1( -6 * SIZE, A1, %xmm1)
	VMOVUPS_A1( -4 * SIZE, A1, %xmm2)
	VMOVUPS_A1( -2 * SIZE, A1, %xmm3)


	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm8)
	VMOVUPS_YS1(-14 * SIZE, Y1, %xmm9)
	VMOVUPS_YS1(-12 * SIZE, Y1, %xmm10)
	VMOVUPS_YS1(-10 * SIZE, Y1, %xmm11)

	VMOVUPS_YL1( -8 * SIZE, Y1, %xmm8)
	VMOVUPS_YL1( -6 * SIZE, Y1, %xmm9)
	VMOVUPS_YL1( -4 * SIZE, Y1, %xmm10)
	VMOVUPS_YL1( -2 * SIZE, Y1, %xmm11)

	subq	 $-8 * SIZE, A1
	subq	 $-8 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(43)
	ALIGN_3

L(44):

	vfmaddpd %xmm8 , %xmm12 , %xmm0 , %xmm8
	vfmaddpd %xmm9 , %xmm12 , %xmm1 , %xmm9
	vfmaddpd %xmm10, %xmm12 , %xmm2 , %xmm10
	vfmaddpd %xmm11, %xmm12 , %xmm3 , %xmm11

	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm8)
	VMOVUPS_YS1(-14 * SIZE, Y1, %xmm9)
	VMOVUPS_YS1(-12 * SIZE, Y1, %xmm10)
	VMOVUPS_YS1(-10 * SIZE, Y1, %xmm11)

	subq	 $-8 * SIZE, A1
	subq	 $-8 * SIZE, Y1
	ALIGN_3

L(45):
	testq	$4, MM
	je	L(46)

	VMOVUPS_A1(-16 * SIZE, A1, %xmm0)
	VMOVUPS_A1(-14 * SIZE, A1, %xmm1)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm8)
	VMOVUPS_YL1(-14 * SIZE, Y1, %xmm9)

	vfmaddpd %xmm8 , %xmm12 , %xmm0 , %xmm8
	vfmaddpd %xmm9 , %xmm12 , %xmm1 , %xmm9

	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm8)
	VMOVUPS_YS1(-14 * SIZE, Y1, %xmm9)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(46):
	testq	$2, MM
	je	L(47)

	VMOVUPS_A1(-16 * SIZE, A1, %xmm8)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm0)

	vfmaddpd %xmm0 , %xmm12 , %xmm8 , %xmm0

	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm0)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(47):
	testq	$1, MM
	je	L(900)

	vmovsd	 -16 * SIZE(Y1), %xmm0
	vmovsd	 -16 * SIZE(A1), %xmm8

	vfmaddsd %xmm0 , %xmm12 , %xmm8 , %xmm0

	vmovsd	 %xmm0, -16 * SIZE(Y1)
	ALIGN_3

#ifdef ALIGNED_ACCESS
	jmp	L(900)
	ALIGN_3

L(50):
#if GEMV_UNROLL >= 4

	cmpq	$4, N
	jl	L(60)
	ALIGN_3

L(51):

	subq	$4, N

	leaq	16 * SIZE(BUFFER), Y1
	movq	A,  A1
	leaq	(A,  LDA, 2), A2
	leaq	(A,  LDA, 4), A

	vmovddup	(X), %xmm12
	addq	INCX, X
	vmovddup	(X), %xmm13
	addq	INCX, X
	vmovddup	(X), %xmm14
	addq	INCX, X
	vmovddup	(X), %xmm15
	addq	INCX, X

	vmovddup	ALPHA, %xmm0

	vmulpd	%xmm0, %xmm12 , %xmm12
	vmulpd	%xmm0, %xmm13 , %xmm13
	vmulpd	%xmm0, %xmm14 , %xmm14
	vmulpd	%xmm0, %xmm15 , %xmm15

	testq	$SIZE, A
	je	L(5X)

	vmovsd	 -16 * SIZE(A1),         %xmm4
	vmovsd	 -16 * SIZE(A1, LDA),    %xmm5
	vmovsd	 -16 * SIZE(A2),         %xmm6
	vmovsd	 -16 * SIZE(A2, LDA),    %xmm7

	vmovsd	 -16 * SIZE(Y1),         %xmm0

	vfmaddsd %xmm0 , %xmm12 , %xmm4 , %xmm0
	vfmaddsd %xmm0 , %xmm13 , %xmm5 , %xmm0
	vfmaddsd %xmm0 , %xmm14 , %xmm6 , %xmm0
	vfmaddsd %xmm0 , %xmm15 , %xmm7 , %xmm0

	vmovsd	 %xmm0, -16 * SIZE(Y1)

	addq	 $SIZE, A1
	addq	 $SIZE, A2
	addq	 $SIZE, Y1
	ALIGN_3

L(5X):
	vmovhpd	-16 * SIZE(A1, LDA), %xmm8, %xmm8
	vmovhpd	-16 * SIZE(A2, LDA), %xmm9, %xmm9

	movq	MM,  I
	sarq	$3,  I
	jle	L(55)

	VMOVUPS_A1(-16 * SIZE, A1, %xmm4)
	VMOVUPS_A1(-14 * SIZE, A1, %xmm5)
	VMOVUPS_A1(-12 * SIZE, A1, %xmm6)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	VMOVUPS_YL1(-14 * SIZE, Y1, %xmm1)
	VMOVUPS_YL1(-12 * SIZE, Y1, %xmm2)
	VMOVUPS_YL1(-10 * SIZE, Y1, %xmm3)

	decq	 I
	jle	 L(54)
	ALIGN_3

L(53):



	vfmaddpd %xmm0 , %xmm12 , %xmm4 , %xmm0
        VMOVUPS_A1(-10 * SIZE, A1, %xmm7)
	vfmaddpd %xmm1 , %xmm12 , %xmm5 , %xmm1
        VMOVUPS_A2(-15 * SIZE, A1, LDA, 1, %xmm4)

	prefetchnta	A_PRE(A1, LDA, 1)
	vfmaddpd %xmm2 , %xmm12 , %xmm6 , %xmm2
        VMOVUPS_A2(-13 * SIZE, A1, LDA, 1, %xmm5)
	vfmaddpd %xmm3 , %xmm12 , %xmm7 , %xmm3
        VMOVUPS_A2(-11 * SIZE, A1, LDA, 1, %xmm6)


        vshufpd  $0x01, %xmm4, %xmm8, %xmm8
	vfmaddpd %xmm0 , %xmm13 , %xmm8 , %xmm0
        VMOVUPS_A2( -9 * SIZE, A1, LDA, 1, %xmm8)
        vshufpd  $0x01, %xmm5, %xmm4, %xmm4
	vfmaddpd %xmm1 , %xmm13 , %xmm4 , %xmm1
        VMOVUPS_A1(-16 * SIZE, A2, %xmm4)

	prefetchnta	A_PRE(A2)
        vshufpd  $0x01, %xmm6, %xmm5, %xmm5
	vfmaddpd %xmm2 , %xmm13 , %xmm5 , %xmm2
        VMOVUPS_A1(-14 * SIZE, A2, %xmm5)
        vshufpd  $0x01, %xmm8, %xmm6, %xmm6
	vfmaddpd %xmm3 , %xmm13 , %xmm6 , %xmm3
        VMOVUPS_A1(-12 * SIZE, A2, %xmm6)


	vfmaddpd %xmm0 , %xmm14 , %xmm4 , %xmm0
        VMOVUPS_A1(-10 * SIZE, A2, %xmm7)
	vfmaddpd %xmm1 , %xmm14 , %xmm5 , %xmm1
        VMOVUPS_A2(-15 * SIZE, A2, LDA, 1, %xmm4)

	prefetchnta	A_PRE(A2, LDA, 1)
	vfmaddpd %xmm2 , %xmm14 , %xmm6 , %xmm2
        VMOVUPS_A2(-13 * SIZE, A2, LDA, 1, %xmm5)
	vfmaddpd %xmm3 , %xmm14 , %xmm7 , %xmm3
        VMOVUPS_A2(-11 * SIZE, A2, LDA, 1, %xmm6)


        vshufpd  $0x01, %xmm4, %xmm9, %xmm9
	vfmaddpd %xmm0 , %xmm15 , %xmm9 , %xmm0
        VMOVUPS_A2( -9 * SIZE, A2, LDA, 1, %xmm9)
        vshufpd  $0x01, %xmm5, %xmm4, %xmm4
	vfmaddpd %xmm1 , %xmm15 , %xmm4 , %xmm1
        VMOVUPS_A1( -8 * SIZE, A1, %xmm4)

	prefetchnta	A_PRE(A1)
        vshufpd  $0x01, %xmm6, %xmm5, %xmm5
	vfmaddpd %xmm2 , %xmm15 , %xmm5 , %xmm2
        VMOVUPS_A1( -6 * SIZE, A1, %xmm5)
        vshufpd  $0x01, %xmm9, %xmm6, %xmm6
	vfmaddpd %xmm3 , %xmm15 , %xmm6 , %xmm3
        VMOVUPS_A1( -4 * SIZE, A1, %xmm6)


        VMOVUPS_YS1(-16 * SIZE, Y1, %xmm0)
        VMOVUPS_YS1(-14 * SIZE, Y1, %xmm1)
        VMOVUPS_YS1(-12 * SIZE, Y1, %xmm2)
        VMOVUPS_YS1(-10 * SIZE, Y1, %xmm3)

        VMOVUPS_YL1( -8 * SIZE, Y1, %xmm0)
        VMOVUPS_YL1( -6 * SIZE, Y1, %xmm1)
	prefetchnta	A_PRE(Y1)
        VMOVUPS_YL1( -4 * SIZE, Y1, %xmm2)
        VMOVUPS_YL1( -2 * SIZE, Y1, %xmm3)

        subq     $-8 * SIZE, A1
        subq     $-8 * SIZE, A2
        subq     $-8 * SIZE, Y1

        subq     $1, I
        BRANCH
        jg      L(53)
        ALIGN_3


L(54):

	vfmaddpd %xmm0 , %xmm12 , %xmm4 , %xmm0
	VMOVUPS_A1(-10 * SIZE, A1, %xmm7)
	vfmaddpd %xmm1 , %xmm12 , %xmm5 , %xmm1
	VMOVUPS_A2(-15 * SIZE, A1, LDA, 1, %xmm4)

	vfmaddpd %xmm2 , %xmm12 , %xmm6 , %xmm2
	VMOVUPS_A2(-13 * SIZE, A1, LDA, 1, %xmm5)
	vfmaddpd %xmm3 , %xmm12 , %xmm7 , %xmm3
	VMOVUPS_A2(-11 * SIZE, A1, LDA, 1, %xmm6)

	vshufpd	$0x01, %xmm4, %xmm8, %xmm8
	vfmaddpd %xmm0 , %xmm13 , %xmm8 , %xmm0
	VMOVUPS_A2( -9 * SIZE, A1, LDA, 1, %xmm8)
	vshufpd	$0x01, %xmm5, %xmm4, %xmm4
	vfmaddpd %xmm1 , %xmm13 , %xmm4 , %xmm1
	VMOVUPS_A1(-16 * SIZE, A2, %xmm4)

	vshufpd	$0x01, %xmm6, %xmm5, %xmm5
	vfmaddpd %xmm2 , %xmm13 , %xmm5 , %xmm2
	VMOVUPS_A1(-14 * SIZE, A2, %xmm5)
	vshufpd	$0x01, %xmm8, %xmm6, %xmm6
	vfmaddpd %xmm3 , %xmm13 , %xmm6 , %xmm3
	VMOVUPS_A1(-12 * SIZE, A2, %xmm6)

	vfmaddpd %xmm0 , %xmm14 , %xmm4 , %xmm0
	VMOVUPS_A1(-10 * SIZE, A2, %xmm7)
	vfmaddpd %xmm1 , %xmm14 , %xmm5 , %xmm1
	VMOVUPS_A2(-15 * SIZE, A2, LDA, 1, %xmm4)

	vfmaddpd %xmm2 , %xmm14 , %xmm6 , %xmm2
	VMOVUPS_A2(-13 * SIZE, A2, LDA, 1, %xmm5)
	vfmaddpd %xmm3 , %xmm14 , %xmm7 , %xmm3
	VMOVUPS_A2(-11 * SIZE, A2, LDA, 1, %xmm6)

	vshufpd	$0x01, %xmm4, %xmm9, %xmm9
	vfmaddpd %xmm0 , %xmm15 , %xmm9 , %xmm0
	VMOVUPS_A2( -9 * SIZE, A2, LDA, 1, %xmm9)

	vshufpd	$0x01, %xmm5, %xmm4, %xmm4
	vfmaddpd %xmm1 , %xmm15 , %xmm4 , %xmm1
	vshufpd	$0x01, %xmm6, %xmm5, %xmm5
	vfmaddpd %xmm2 , %xmm15 , %xmm5 , %xmm2
	vshufpd	$0x01, %xmm9, %xmm6, %xmm6
	vfmaddpd %xmm3 , %xmm15 , %xmm6 , %xmm3

	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm0)
	VMOVUPS_YS1(-14 * SIZE, Y1, %xmm1)
	VMOVUPS_YS1(-12 * SIZE, Y1, %xmm2)
	VMOVUPS_YS1(-10 * SIZE, Y1, %xmm3)

	subq	 $-8 * SIZE, A1
	subq	 $-8 * SIZE, A2
	subq	 $-8 * SIZE, Y1
	ALIGN_3

L(55):
	testq	$4, MM
	je	L(56)

	VMOVUPS_A1(-16 * SIZE, A1, %xmm4)
	VMOVUPS_A1(-14 * SIZE, A1, %xmm5)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	VMOVUPS_YL1(-14 * SIZE, Y1, %xmm1)

	vfmaddpd %xmm0 , %xmm12 , %xmm4 , %xmm0
	vfmaddpd %xmm1 , %xmm12 , %xmm5 , %xmm1

	VMOVUPS_A2(-15 * SIZE, A1, LDA, 1, %xmm6)
	VMOVUPS_A2(-13 * SIZE, A1, LDA, 1, %xmm7)

	vshufpd	$0x01, %xmm6, %xmm8, %xmm8
	vfmaddpd %xmm0 , %xmm13 , %xmm8 , %xmm0
	movaps	 %xmm7,  %xmm8
	vshufpd	$0x01, %xmm7, %xmm6, %xmm6
	vfmaddpd %xmm1 , %xmm13 , %xmm6 , %xmm1

	VMOVUPS_A1(-16 * SIZE, A2, %xmm4)
	VMOVUPS_A1(-14 * SIZE, A2, %xmm5)

	vfmaddpd %xmm0 , %xmm14 , %xmm4 , %xmm0
	vfmaddpd %xmm1 , %xmm14 , %xmm5 , %xmm1

	VMOVUPS_A2(-15 * SIZE, A2, LDA, 1, %xmm6)
	VMOVUPS_A2(-13 * SIZE, A2, LDA, 1, %xmm7)

	vshufpd	$0x01, %xmm6, %xmm9, %xmm9
	vfmaddpd %xmm0 , %xmm15 , %xmm9 , %xmm0
	movaps	 %xmm7,  %xmm9
	vshufpd	$0x01, %xmm7, %xmm6, %xmm6
	vfmaddpd %xmm1 , %xmm15 , %xmm6 , %xmm1

	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm0)
	VMOVUPS_YS1(-14 * SIZE, Y1, %xmm1)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, A2
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(56):
	testq	$2, MM
	je	L(57)

	VMOVUPS_A1(-16 * SIZE, A1, %xmm4)
	VMOVUPS_A2(-15 * SIZE, A1, LDA, 1, %xmm5)
	VMOVUPS_A1(-16 * SIZE, A2, %xmm6)
	VMOVUPS_A2(-15 * SIZE, A2, LDA, 1, %xmm7)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm0)

	vfmaddpd %xmm0 , %xmm12 , %xmm4 , %xmm0
	vshufpd	$0x01, %xmm5, %xmm8, %xmm8
	vfmaddpd %xmm0 , %xmm13 , %xmm8 , %xmm0
	movaps	 %xmm5,  %xmm8
	vfmaddpd %xmm0 , %xmm14 , %xmm6 , %xmm0
	vshufpd	$0x01, %xmm7, %xmm9, %xmm9
	vfmaddpd %xmm0 , %xmm15 , %xmm9 , %xmm0
	movaps	 %xmm7,  %xmm9

	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm0)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(57):
	testq	$1, MM
	je	L(58)

	vmovsd	 -16 * SIZE(Y1), %xmm0

	vmovsd	 -16 * SIZE(A1), %xmm4
	vshufpd	 $0x01, %xmm8, %xmm8, %xmm8
	vmovsd	 -16 * SIZE(A2), %xmm6
	vshufpd	 $0x01, %xmm9, %xmm9, %xmm9

	vfmaddsd %xmm0 , %xmm12 , %xmm4 , %xmm0
	vfmaddsd %xmm0 , %xmm13 , %xmm8 , %xmm0
	vfmaddsd %xmm0 , %xmm14 , %xmm6 , %xmm0
	vfmaddsd %xmm0 , %xmm15 , %xmm9 , %xmm0

	vmovsd	 %xmm0, -16 * SIZE(Y1)
	ALIGN_3

L(58):
	cmpq	$4, N
	jge	L(51)
	ALIGN_3

L(60):
#endif

#if GEMV_UNROLL >= 2

	cmpq	$2, N
	jl	L(70)

#if GEMV_UNROLL == 2
	ALIGN_3

L(61):
#endif

	subq	$2, N

	leaq	16 * SIZE(BUFFER), Y1
	movq	A,  A1
	leaq	(A,  LDA), A2
	leaq	(A,  LDA, 2), A

	vmovddup	(X), %xmm12
	addq	INCX, X
	vmovddup	(X), %xmm13
	addq	INCX, X

	vmovddup	ALPHA, %xmm0

	vmulpd	%xmm0, %xmm12 , %xmm12
	vmulpd	%xmm0, %xmm13 , %xmm13

	testq	$SIZE, A
	je	L(6X)

	vmovsd	 -16 * SIZE(A1), %xmm4
	vmovsd	 -16 * SIZE(A2), %xmm5

	vmovsd	 -16 * SIZE(Y1), %xmm0

	vfmaddsd %xmm0 , %xmm12 , %xmm4 , %xmm0
	vfmaddsd %xmm0 , %xmm13 , %xmm5 , %xmm0

	vmovsd	 %xmm0, -16 * SIZE(Y1)

	addq	 $SIZE, A1
	addq	 $SIZE, A2
	addq	 $SIZE, Y1
	ALIGN_3

L(6X):
	vmovhpd	-16 * SIZE(A2), %xmm8, %xmm8

	movq	MM,  I
	sarq	$3,  I
	jle	L(65)

	VMOVUPS_A1(-16 * SIZE, A1, %xmm4)
	VMOVUPS_A1(-14 * SIZE, A1, %xmm5)
	VMOVUPS_A1(-12 * SIZE, A1, %xmm6)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	VMOVUPS_YL1(-14 * SIZE, Y1, %xmm1)
	VMOVUPS_YL1(-12 * SIZE, Y1, %xmm2)
	VMOVUPS_YL1(-10 * SIZE, Y1, %xmm3)

	decq	 I
	jle	 L(64)
	ALIGN_3

L(63):

	vfmaddpd %xmm0 , %xmm12 , %xmm4 , %xmm0
	VMOVUPS_A1(-10 * SIZE, A1, %xmm7)
	vfmaddpd %xmm1 , %xmm12 , %xmm5 , %xmm1
	VMOVUPS_A1(-15 * SIZE, A2, %xmm4)

	prefetchnta	A_PRE(A2)
	vfmaddpd %xmm2 , %xmm12 , %xmm6 , %xmm2
	VMOVUPS_A1(-13 * SIZE, A2, %xmm5)
	vfmaddpd %xmm3 , %xmm12 , %xmm7 , %xmm3
	VMOVUPS_A1(-11 * SIZE, A2, %xmm6)


	vshufpd	$0x01, %xmm4, %xmm8, %xmm8
	vfmaddpd %xmm0 , %xmm13 , %xmm8 , %xmm0
	VMOVUPS_A1( -9 * SIZE, A2, %xmm8)
	vshufpd	$0x01, %xmm5, %xmm4, %xmm4
	vfmaddpd %xmm1 , %xmm13 , %xmm4 , %xmm1
	VMOVUPS_A1( -8 * SIZE, A1, %xmm4)

	prefetchnta	A_PRE(A1)
	vshufpd	$0x01, %xmm6, %xmm5, %xmm5
	vfmaddpd %xmm2 , %xmm13 , %xmm5 , %xmm2
	VMOVUPS_A1( -6 * SIZE, A1, %xmm5)
	vshufpd	$0x01, %xmm8, %xmm6, %xmm6
	vfmaddpd %xmm3 , %xmm13 , %xmm6 , %xmm3
	VMOVUPS_A1( -4 * SIZE, A1, %xmm6)


	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm0)
	VMOVUPS_YS1(-14 * SIZE, Y1, %xmm1)
	VMOVUPS_YS1(-12 * SIZE, Y1, %xmm2)
	VMOVUPS_YS1(-10 * SIZE, Y1, %xmm3)

	VMOVUPS_YL1( -8 * SIZE, Y1, %xmm0)
	VMOVUPS_YL1( -6 * SIZE, Y1, %xmm1)
	prefetchnta	A_PRE(Y1)
	VMOVUPS_YL1( -4 * SIZE, Y1, %xmm2)
	VMOVUPS_YL1( -2 * SIZE, Y1, %xmm3)

	subq	 $-8 * SIZE, A1
	subq	 $-8 * SIZE, A2
	subq	 $-8 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(63)
	ALIGN_3

L(64):

	vfmaddpd %xmm0 , %xmm12 , %xmm4 , %xmm0
	VMOVUPS_A1(-10 * SIZE, A1, %xmm7)
	vfmaddpd %xmm1 , %xmm12 , %xmm5 , %xmm1
	VMOVUPS_A1(-15 * SIZE, A2, %xmm4)

	vfmaddpd %xmm2 , %xmm12 , %xmm6 , %xmm2
	VMOVUPS_A1(-13 * SIZE, A2, %xmm5)
	vfmaddpd %xmm3 , %xmm12 , %xmm7 , %xmm3
	VMOVUPS_A1(-11 * SIZE, A2, %xmm6)

	vshufpd	$0x01, %xmm4, %xmm8, %xmm8
	vfmaddpd %xmm0 , %xmm13 , %xmm8 , %xmm0
	VMOVUPS_A1( -9 * SIZE, A2, %xmm8)
	vshufpd	$0x01, %xmm5, %xmm4, %xmm4
	vfmaddpd %xmm1 , %xmm13 , %xmm4 , %xmm1

	vshufpd	$0x01, %xmm6, %xmm5, %xmm5
	vfmaddpd %xmm2 , %xmm13 , %xmm5 , %xmm2
	vshufpd	$0x01, %xmm8, %xmm6, %xmm6
	vfmaddpd %xmm3 , %xmm13 , %xmm6 , %xmm3

	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm0)
	VMOVUPS_YS1(-14 * SIZE, Y1, %xmm1)
	VMOVUPS_YS1(-12 * SIZE, Y1, %xmm2)
	VMOVUPS_YS1(-10 * SIZE, Y1, %xmm3)

	subq	 $-8 * SIZE, A1
	subq	 $-8 * SIZE, A2
	subq	 $-8 * SIZE, Y1
	ALIGN_3

L(65):
	testq	$4, MM
	je	L(66)


	VMOVUPS_A1(-16 * SIZE, A1, %xmm4)
	VMOVUPS_A1(-14 * SIZE, A1, %xmm5)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	VMOVUPS_YL1(-14 * SIZE, Y1, %xmm1)

	vfmaddpd %xmm0 , %xmm12 , %xmm4 , %xmm0
	vfmaddpd %xmm1 , %xmm12 , %xmm5 , %xmm1

	VMOVUPS_A1(-15 * SIZE, A2, %xmm6)
	VMOVUPS_A1(-13 * SIZE, A2, %xmm7)

	vshufpd	$0x01, %xmm6, %xmm8, %xmm8
	vfmaddpd %xmm0 , %xmm13 , %xmm8 , %xmm0
	movaps	 %xmm7,  %xmm8
	vshufpd	$0x01, %xmm7, %xmm6, %xmm6
	vfmaddpd %xmm1 , %xmm13 , %xmm6 , %xmm1

	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm0)
	VMOVUPS_YS1(-14 * SIZE, Y1, %xmm1)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, A2
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(66):
	testq	$2, MM
	je	L(67)

	VMOVUPS_A1(-16 * SIZE, A1, %xmm4)
	VMOVUPS_A1(-15 * SIZE, A2, %xmm5)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm0)

	vfmaddpd %xmm0 , %xmm12 , %xmm4 , %xmm0
	vshufpd	$0x01, %xmm5, %xmm8, %xmm8
	vfmaddpd %xmm0 , %xmm13 , %xmm8 , %xmm0
	movaps	 %xmm5,  %xmm8

	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm0)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(67):
	testq	$1, MM
#if GEMV_UNROLL == 2
	je	L(68)
#else
	je	L(70)
#endif

	vmovsd	 -16 * SIZE(Y1), %xmm0

	vmovsd	 -16 * SIZE(A1), %xmm4
	vshufpd	 $0x01, %xmm8, %xmm8 , %xmm8

	vfmaddsd %xmm0 , %xmm12 , %xmm4 , %xmm0
	vfmaddsd %xmm0 , %xmm13 , %xmm8 , %xmm0

	vmovsd	 %xmm0, -16 * SIZE(Y1)
	ALIGN_3

#if GEMV_UNROLL == 2
L(68):
	cmpq	$2, N
	jge	L(61)
	ALIGN_3

#endif

L(70):
	cmpq	$1, N
	jl	L(900)

#endif

	leaq	16 * SIZE(BUFFER), Y1
	movq	A,  A1

	vmovddup	(X), %xmm12
	addq	INCX, X

	vmovddup	ALPHA, %xmm0

	vmulpd	%xmm0, %xmm12 , %xmm12

	testq	$SIZE, A
	je	L(7X)

	vmovsd	 -16 * SIZE(A1), %xmm4
	vmovsd	 -16 * SIZE(Y1), %xmm0

	vfmaddsd %xmm0 , %xmm12 , %xmm4 , %xmm0

	vmovsd	 %xmm0, -16 * SIZE(Y1)

	addq	 $SIZE, A1
	addq	 $SIZE, Y1
	ALIGN_3

L(7X):

	movq	MM,  I
	sarq	$3,  I
	jle	L(75)

	VMOVUPS_A1(-16 * SIZE, A1, %xmm0)
	VMOVUPS_A1(-14 * SIZE, A1, %xmm1)
	VMOVUPS_A1(-12 * SIZE, A1, %xmm2)
	VMOVUPS_A1(-10 * SIZE, A1, %xmm3)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm8)
	VMOVUPS_YL1(-14 * SIZE, Y1, %xmm9)
	VMOVUPS_YL1(-12 * SIZE, Y1, %xmm10)
	VMOVUPS_YL1(-10 * SIZE, Y1, %xmm11)

	decq	 I
	jle	 L(74)
	ALIGN_3

L(73):

	vfmaddpd %xmm8 , %xmm12 , %xmm0 , %xmm8
	VMOVUPS_A1( -8 * SIZE, A1, %xmm0)
	vfmaddpd %xmm9 , %xmm12 , %xmm1 , %xmm9
	VMOVUPS_A1( -6 * SIZE, A1, %xmm1)

	prefetchnta	A_PRE(A1)
	vfmaddpd %xmm10, %xmm12 , %xmm2 , %xmm10
	VMOVUPS_A1( -4 * SIZE, A1, %xmm2)
	vfmaddpd %xmm11, %xmm12 , %xmm3 , %xmm11
	VMOVUPS_A1( -2 * SIZE, A1, %xmm3)


	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm8)
	VMOVUPS_YS1(-14 * SIZE, Y1, %xmm9)
	VMOVUPS_YS1(-12 * SIZE, Y1, %xmm10)
	VMOVUPS_YS1(-10 * SIZE, Y1, %xmm11)

	VMOVUPS_YL1( -8 * SIZE, Y1, %xmm8)
	VMOVUPS_YL1( -6 * SIZE, Y1, %xmm9)
	prefetchnta	A_PRE(Y1)
	VMOVUPS_YL1( -4 * SIZE, Y1, %xmm10)
	VMOVUPS_YL1( -2 * SIZE, Y1, %xmm11)

	subq	 $-8 * SIZE, A1
	subq	 $-8 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(73)
	ALIGN_3

L(74):

	vfmaddpd %xmm8 , %xmm12 , %xmm0 , %xmm8
	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm8)
	vfmaddpd %xmm9 , %xmm12 , %xmm1 , %xmm9
	VMOVUPS_YS1(-14 * SIZE, Y1, %xmm9)
	vfmaddpd %xmm10, %xmm12 , %xmm2 , %xmm10
	VMOVUPS_YS1(-12 * SIZE, Y1, %xmm10)
	vfmaddpd %xmm11, %xmm12 , %xmm3 , %xmm11
	VMOVUPS_YS1(-10 * SIZE, Y1, %xmm11)

	subq	 $-8 * SIZE, A1
	subq	 $-8 * SIZE, Y1
	ALIGN_3

L(75):
	testq	$4, MM
	je	L(76)

	VMOVUPS_A1(-16 * SIZE, A1, %xmm0)
	VMOVUPS_A1(-14 * SIZE, A1, %xmm1)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm8)
	VMOVUPS_YL1(-14 * SIZE, Y1, %xmm9)

	vfmaddpd %xmm8 , %xmm12 , %xmm0 , %xmm8
	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm8)
	vfmaddpd %xmm9 , %xmm12 , %xmm1 , %xmm9
	VMOVUPS_YS1(-14 * SIZE, Y1, %xmm9)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(76):
	testq	$2, MM
	je	L(77)

	VMOVUPS_A1(-16 * SIZE, A1, %xmm8)

	VMOVUPS_YL1(-16 * SIZE, Y1, %xmm0)

	vfmaddpd %xmm0 , %xmm12 , %xmm8 , %xmm0

	VMOVUPS_YS1(-16 * SIZE, Y1, %xmm0)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(77):
	testq	$1, MM
	je	L(900)

	vmovsd	 -16 * SIZE(Y1), %xmm0
	vmovsd	 -16 * SIZE(A1), %xmm8

	vfmaddsd %xmm0 , %xmm12 , %xmm8 , %xmm0

	vmovsd	 %xmm0, -16 * SIZE(Y1)
#endif
	ALIGN_3


L(900):
#ifndef COPY_FORCE
	cmpq	Y, BUFFER
	je	L(999)
#endif

	cmpq	$SIZE, INCY
	jne	L(950)

	testq	$SIZE, Y
	je	L(910)

	vmovsd	(Y), %xmm0
	vaddsd	(BUFFER), %xmm0, %xmm0
	vmovsd	%xmm0, (Y)

	addq	$SIZE, Y
	addq	$SIZE, BUFFER

	decq	M
	jle	L(999)
	ALIGN_4

L(910):
	testq	$SIZE, BUFFER
	jne	L(920)

	movq	M,   %rax
	sarq	$3,  %rax
	jle	L(914)
	ALIGN_3

L(912):

	vmovups	0 * SIZE(Y), %xmm0
	vmovups	2 * SIZE(Y), %xmm1
	vmovups	4 * SIZE(Y), %xmm2
	vmovups	6 * SIZE(Y), %xmm3

	vmovups	0 * SIZE(BUFFER), %xmm4
	vmovups	2 * SIZE(BUFFER), %xmm5
	vmovups	4 * SIZE(BUFFER), %xmm6
	vmovups	6 * SIZE(BUFFER), %xmm7


	vaddpd	%xmm4, %xmm0, %xmm0
	vaddpd	%xmm5, %xmm1, %xmm1
	vaddpd	%xmm6, %xmm2, %xmm2
	vaddpd	%xmm7, %xmm3, %xmm3

	vmovups	%xmm0,  0 * SIZE(Y)
	vmovups	%xmm1,  2 * SIZE(Y)
	vmovups	%xmm2,  4 * SIZE(Y)
	vmovups	%xmm3,  6 * SIZE(Y)

	addq	$8 * SIZE, Y
	addq	$8 * SIZE, BUFFER

	decq	%rax
	jg	L(912)
	ALIGN_3

L(914):
	testq	$7, M
	jle	L(999)

	testq	$4, M
	jle	L(915)

	vmovups	0 * SIZE(Y), %xmm0
	vmovups	2 * SIZE(Y), %xmm1

	vmovups	0 * SIZE(BUFFER), %xmm4
	vmovups	2 * SIZE(BUFFER), %xmm5

	vaddpd	%xmm4, %xmm0, %xmm0
	vaddpd	%xmm5, %xmm1, %xmm1

	vmovups	%xmm0,  0 * SIZE(Y)
	vmovups	%xmm1,  2 * SIZE(Y)

	addq	$4 * SIZE, Y
	addq	$4 * SIZE, BUFFER
	ALIGN_3

L(915):
	testq	$2, M
	jle	L(916)

	vmovups	(Y), %xmm0

	vmovups	(BUFFER), %xmm4

	vaddpd	%xmm4, %xmm0, %xmm0

	vmovups	%xmm0,  (Y)

	addq	$2 * SIZE, Y
	addq	$2 * SIZE, BUFFER
	ALIGN_3

L(916):
	testq	$1, M
	jle	L(999)

	vmovsd	(Y), %xmm0

	vmovsd	0 * SIZE(BUFFER), %xmm4

	vaddsd	%xmm4, %xmm0, %xmm0

	vmovsd	%xmm0, (Y)
	ALIGN_3

	jmp	L(999)
	ALIGN_4

L(920):
	vmovups	-1 * SIZE(BUFFER), %xmm4

	movq	M,   %rax
	sarq	$3,  %rax
	jle	L(924)
	ALIGN_3

L(922):

	vmovups	 0 * SIZE(Y), %xmm0
	vmovups	 2 * SIZE(Y), %xmm1
	vmovups	 4 * SIZE(Y), %xmm2
	vmovups	 6 * SIZE(Y), %xmm3

	vmovups	1 * SIZE(BUFFER), %xmm5
	vmovups	3 * SIZE(BUFFER), %xmm6
	vmovups	5 * SIZE(BUFFER), %xmm7
	vmovups	7 * SIZE(BUFFER), %xmm8

	vshufpd	 $0x01, %xmm5, %xmm4, %xmm4
	vshufpd	 $0x01, %xmm6, %xmm5, %xmm5
	vshufpd	 $0x01, %xmm7, %xmm6, %xmm6
	vshufpd	 $0x01, %xmm8, %xmm7, %xmm7


	vaddpd	%xmm4, %xmm0, %xmm0
	vaddpd	%xmm5, %xmm1, %xmm1
	vaddpd	%xmm6, %xmm2, %xmm2
	vaddpd	%xmm7, %xmm3, %xmm3

	vmovups	%xmm0,  0 * SIZE(Y)
	vmovups	%xmm1,  2 * SIZE(Y)
	vmovups	%xmm2,  4 * SIZE(Y)
	vmovups	%xmm3,  6 * SIZE(Y)

	vmovups	%xmm8, %xmm4

	addq	$8 * SIZE, Y
	addq	$8 * SIZE, BUFFER

	decq	%rax
	jg	L(922)
	ALIGN_3

L(924):
	testq	$7, M
	jle	L(999)

	testq	$4, M
	jle	L(925)

	vmovups	 0 * SIZE(Y), %xmm0
	vmovups	 2 * SIZE(Y), %xmm1

	vmovups	1 * SIZE(BUFFER), %xmm5
	vmovups	3 * SIZE(BUFFER), %xmm6

	vshufpd	 $0x01, %xmm5, %xmm4, %xmm4
	vshufpd	 $0x01, %xmm6, %xmm5, %xmm5

	vaddpd	%xmm4, %xmm0, %xmm0
	vaddpd	%xmm5, %xmm1, %xmm1

	vmovups	%xmm0,  0 * SIZE(Y)
	vmovups	%xmm1,  2 * SIZE(Y)

	vmovups	%xmm6, %xmm4

	addq	$4 * SIZE, Y
	addq	$4 * SIZE, BUFFER
	ALIGN_3

L(925):
	testq	$2, M
	jle	L(926)

	vmovups	 (Y), %xmm0

	vmovups	 1 * SIZE(BUFFER), %xmm5

	vshufpd	 $0x01, %xmm5, %xmm4, %xmm4

	vaddpd	%xmm4, %xmm0, %xmm0

	vmovups	%xmm0,  (Y)

	movaps	%xmm5, %xmm4

	addq	$2 * SIZE, Y
	addq	$2 * SIZE, BUFFER
	ALIGN_3

L(926):
	testq	$1, M
	jle	L(999)

	vmovsd	(Y), %xmm0

	vshufpd	$0x01, %xmm4 ,%xmm4, %xmm4

	vaddsd	%xmm4, %xmm0, %xmm0

	vmovsd	%xmm0, (Y)
	ALIGN_3

	jmp	L(999)
	ALIGN_4

L(950):
	testq	$SIZE, BUFFER
	je	L(960)

	vmovsd	(Y), %xmm0
	vaddsd	(BUFFER), %xmm0, %xmm0
	vmovsd	%xmm0, (Y)

	addq	INCY, Y
	addq	$SIZE, BUFFER

	decq	M
	jle	L(999)
	ALIGN_4

L(960):
	movq	Y,  Y1

	movq	M,   %rax
	sarq	$3,  %rax
	jle	L(964)
	ALIGN_3

L(962):
	vmovsd	(Y), %xmm0
	addq	INCY, Y
	vmovhpd	(Y), %xmm0, %xmm0
	addq	INCY, Y

	vmovups	0 * SIZE(BUFFER), %xmm4

	vmovsd	(Y), %xmm1
	addq	INCY, Y
	vmovhpd	(Y), %xmm1, %xmm1
	addq	INCY, Y

	vmovups	2 * SIZE(BUFFER), %xmm5

	vmovsd	(Y), %xmm2
	addq	INCY, Y
	vmovhpd	(Y), %xmm2, %xmm2
	addq	INCY, Y

	vmovups	4 * SIZE(BUFFER), %xmm6

	vaddpd	%xmm4, %xmm0, %xmm0

	vmovsd	(Y), %xmm3
	addq	INCY, Y
	vmovhpd	(Y), %xmm3, %xmm3
	addq	INCY, Y

	vmovups	6 * SIZE(BUFFER), %xmm7

	vaddpd	%xmm5, %xmm1, %xmm1

	vmovsd	%xmm0, (Y1)
	addq	INCY, Y1
	vmovhpd	%xmm0, (Y1)
	addq	INCY, Y1

	vaddpd	%xmm6, %xmm2, %xmm2

	vmovsd	%xmm1, (Y1)
	addq	INCY, Y1
	vmovhpd	%xmm1, (Y1)
	addq	INCY, Y1

	vaddpd	%xmm7, %xmm3, %xmm3

	vmovsd	%xmm2, (Y1)
	addq	INCY, Y1
	vmovhpd	%xmm2, (Y1)
	addq	INCY, Y1
	vmovsd	%xmm3, (Y1)
	addq	INCY, Y1
	vmovhpd	%xmm3, (Y1)
	addq	INCY, Y1

	addq	$8 * SIZE, BUFFER
	decq	%rax
	jg	L(962)
	ALIGN_3

L(964):
	testq	$7, M
	jle	L(999)

	testq	$4, M
	jle	L(965)

	vmovsd	(Y), %xmm0
	addq	INCY, Y
	vmovhpd	(Y), %xmm0, %xmm0
	addq	INCY, Y

	vmovups	0 * SIZE(BUFFER), %xmm4

	vmovsd	(Y), %xmm1
	addq	INCY, Y
	vmovhpd	(Y), %xmm1, %xmm1
	addq	INCY, Y

	vmovups	2 * SIZE(BUFFER), %xmm5

	vaddpd	%xmm4, %xmm0, %xmm0
	vaddpd	%xmm5, %xmm1, %xmm1

	vmovsd	%xmm0, (Y1)
	addq	INCY, Y1
	vmovhpd	%xmm0, (Y1)
	addq	INCY, Y1
	vmovsd	%xmm1, (Y1)
	addq	INCY, Y1
	vmovhpd	%xmm1, (Y1)
	addq	INCY, Y1

	addq	$4 * SIZE, BUFFER
	ALIGN_3

L(965):
	testq	$2, M
	jle	L(966)

	vmovsd	(Y), %xmm0
	addq	INCY, Y
	vmovhpd	(Y),%xmm0, %xmm0
	addq	INCY, Y

	vmovups	0 * SIZE(BUFFER), %xmm4

	vaddpd	%xmm4, %xmm0, %xmm0

	vmovsd	%xmm0, (Y1)
	addq	INCY, Y1
	vmovhpd	%xmm0, (Y1)
	addq	INCY, Y1

	addq	$2 * SIZE, BUFFER
	ALIGN_3

L(966):
	testq	$1, M
	jle	L(999)

	vmovsd	(Y), %xmm0

	vmovsd	0 * SIZE(BUFFER), %xmm4

	vaddsd	%xmm4, %xmm0, %xmm0

	vmovsd	%xmm0, (Y1)
	ALIGN_3

L(999):
	movq	  0(%rsp), %rbx
	movq	  8(%rsp), %rbp
	movq	 16(%rsp), %r12
	movq	 24(%rsp), %r13
	movq	 32(%rsp), %r14
	movq	 40(%rsp), %r15

#ifdef WINDOWS_ABI
	movq	 48(%rsp), %rdi
	movq	 56(%rsp), %rsi
	vmovups	 64(%rsp), %xmm6
	vmovups	 80(%rsp), %xmm7
	vmovups	 96(%rsp), %xmm8
	vmovups	112(%rsp), %xmm9
	vmovups	128(%rsp), %xmm10
	vmovups	144(%rsp), %xmm11
	vmovups	160(%rsp), %xmm12
	vmovups	176(%rsp), %xmm13
	vmovups	192(%rsp), %xmm14
	vmovups	208(%rsp), %xmm15
#endif

	addq	$STACKSIZE, %rsp

	ret
	EPILOGUE
