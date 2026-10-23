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

#if GEMV_UNROLL < 4
#undef  GEMV_UNROLL
#define GEMV_UNROLL 4
#endif

#ifndef WINDOWS_ABI

#define STACKSIZE	128

#define OLD_M	  %rdi
#define OLD_N	  %rsi
#define OLD_A	  %rcx
#define OLD_LDA	  %r8
#define STACK_INCX	 8 + STACKSIZE(%rsp)
#define STACK_Y		16 + STACKSIZE(%rsp)
#define STACK_INCY	24 + STACKSIZE(%rsp)
#define STACK_BUFFER	32 + STACKSIZE(%rsp)
#define ALPHA	48	      (%rsp)
#define MMM	56(%rsp)
#define NN	64(%rsp)
#define AA	72(%rsp)
#define LDAX	80(%rsp)
#define XX	96(%rsp)
#else

#define STACKSIZE	288

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

#define MMM	232(%rsp)
#define NN	240(%rsp)
#define AA	248(%rsp)
#define LDAX 256(%rsp)
#define XX	264(%rsp)

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

#ifndef WINDOWS_ABI
	movss	 %xmm0, ALPHA
#else
	movss	 %xmm3, ALPHA
#endif


	movq	M,MMM
	movq	A,AA
	movq	N,NN
	movq	LDA,LDAX
	movq	X,XX
	movq	STACK_Y,      Y
L(0t):
	xorq	I,I
	addq	$1,I
	salq	$22,I
	subq	I,MMM
	movq	I,M
	jge	L(00t)

	movq	MMM,M
	addq	I,M
	jle	L(999x)

L(00t):
	movq	AA,A
	movq	NN,N
	movq	LDAX,LDA
	movq	XX,X

	movq	STACK_INCX,   INCX
	movq	STACK_INCY,   INCY
	movq	STACK_BUFFER, BUFFER

	leaq	(,INCX, SIZE), INCX
	leaq	(,INCY, SIZE), INCY
	leaq	(,LDA,  SIZE), LDA

	leaq	(LDA, LDA, 2), LDA3

#ifdef ALIGNED_ACCESS
	movq	M, MM
	testq	$4 * SIZE - 1, A
	je	L(0X)
	cmpq	$3, M
	jle	L(0X)

	movq	A, MM
	sarq	$BASE_SHIFT, MM
	andq	$3, MM
	subq	$4, MM
	addq	M,  MM

L(0X):
#endif

	testq	N, N		# if n <= 0 goto END
	jle	L(999)
	testq	M, M		# if n <= 0 goto END
	jle	L(999)

	subq	$-32 * SIZE, A

	movq	BUFFER, Y1

	pxor	%xmm0, %xmm0

	movq	M,   %rax
#ifdef ALIGNED_ACCESS
	addq	$19, %rax
#else
	addq	$16, %rax
#endif
	sarq	$4,  %rax
	ALIGN_3

L(01):
	movaps	%xmm0,  0 * SIZE(Y1)
	movaps	%xmm0,  4 * SIZE(Y1)
	movaps	%xmm0,  8 * SIZE(Y1)
	movaps	%xmm0, 12 * SIZE(Y1)
	addq	$16 * SIZE, Y1
	decq	%rax
	jg	L(01)
	ALIGN_3

L(10):
#ifdef ALIGNED_ACCESS
	movq	A,  %rax
	andq	$4 * SIZE - 1, %rax
	addq	%rax, BUFFER

	testq	$4 * SIZE - 1, LDA
	jne	L(100)
#endif

#if GEMV_UNROLL >= 8

	cmpq	$8, N
	jl	L(20)
	ALIGN_3

L(11):
	subq	$8, N

	leaq	32 * SIZE(BUFFER), Y1
	movq	A,  A1
	leaq	(A,  LDA, 4), A2
	leaq	(A,  LDA, 8), A

	movss	(X), %xmm8
	addq	INCX, X
	movss	(X), %xmm9
	addq	INCX, X
	movss	(X), %xmm10
	addq	INCX, X
	movss	(X), %xmm11
	addq	INCX, X
	movss	(X), %xmm12
	addq	INCX, X
	movss	(X), %xmm13
	addq	INCX, X
	movss	(X), %xmm14
	addq	INCX, X
	movss	(X), %xmm15
	addq	INCX, X

	movss	ALPHA, %xmm0

	mulss	%xmm0, %xmm8
	shufps	$0, %xmm8,  %xmm8
	mulss	%xmm0, %xmm9
	shufps	$0, %xmm9,  %xmm9
	mulss	%xmm0, %xmm10
	shufps	$0, %xmm10, %xmm10
	mulss	%xmm0, %xmm11
	shufps	$0, %xmm11, %xmm11

	mulss	%xmm0, %xmm12
	shufps	$0, %xmm12, %xmm12
	mulss	%xmm0, %xmm13
	shufps	$0, %xmm13, %xmm13
	mulss	%xmm0, %xmm14
	shufps	$0, %xmm14, %xmm14
	mulss	%xmm0, %xmm15
	shufps	$0, %xmm15, %xmm15

#ifdef ALIGNED_ACCESS
	cmpq	$3, M
	jle	L(17)

	testq	$SIZE, A1
	je	L(1X)

	movss	 -32 * SIZE(A1), %xmm4
	movss	 -32 * SIZE(A1, LDA,  1), %xmm5
	movss	 -32 * SIZE(A1, LDA,  2), %xmm6
	movss	 -32 * SIZE(A1, LDA3, 1), %xmm7

	movss	 -32 * SIZE(Y1), %xmm0

	mulss	 %xmm8,  %xmm4
	addss	 %xmm4,  %xmm0
	movss	 -32 * SIZE(A2), %xmm4
	mulss	 %xmm9,  %xmm5
	addss	 %xmm5,  %xmm0
	movss	 -32 * SIZE(A2, LDA,  1), %xmm5
	mulss	 %xmm10, %xmm6
	addss	 %xmm6,  %xmm0
	movss	 -32 * SIZE(A2, LDA,  2), %xmm6
	mulss	 %xmm11, %xmm7
	addss	 %xmm7,  %xmm0
	movss	 -32 * SIZE(A2, LDA3, 1), %xmm7

	mulss	 %xmm12, %xmm4
	addss	 %xmm4,  %xmm0
	mulss	 %xmm13, %xmm5
	addss	 %xmm5,  %xmm0
	mulss	 %xmm14, %xmm6
	addss	 %xmm6,  %xmm0
	mulss	 %xmm15, %xmm7
	addss	 %xmm7,  %xmm0

	movss	 %xmm0, -32 * SIZE(Y1)

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, Y1
	ALIGN_3

L(1X):
	testq	$2 * SIZE, A1
	je	L(1XX)

	movsd	 -32 * SIZE(A1), %xmm4
	movsd	 -32 * SIZE(A1, LDA,  1), %xmm5
	movsd	 -32 * SIZE(A1, LDA,  2), %xmm6
	movsd	 -32 * SIZE(A1, LDA3, 1), %xmm7

	movsd	 -32 * SIZE(Y1), %xmm0

	mulps	 %xmm8,  %xmm4
	addps	 %xmm4,  %xmm0
	movsd	 -32 * SIZE(A2), %xmm4
	mulps	 %xmm9,  %xmm5
	addps	 %xmm5,  %xmm0
	movsd	 -32 * SIZE(A2, LDA,  1), %xmm5
	mulps	 %xmm10, %xmm6
	addps	 %xmm6,  %xmm0
	movsd	 -32 * SIZE(A2, LDA,  2), %xmm6
	mulps	 %xmm11, %xmm7
	addps	 %xmm7,  %xmm0
	movsd	 -32 * SIZE(A2, LDA3, 1), %xmm7

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm0
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm0
	mulps	 %xmm15, %xmm7
	addps	 %xmm7,  %xmm0

	movlps	 %xmm0, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(1XX):
#endif
	movq	MM,  I
	sarq	$4,  I
	jle	L(15)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm4)
	MOVUPS_A1 (-28 * SIZE, A1, %xmm5)
	MOVUPS_A1 (-24 * SIZE, A1, %xmm6)
	MOVUPS_A1 (-20 * SIZE, A1, %xmm7)

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YL1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YL1(-20 * SIZE, Y1, %xmm3)

	decq	 I
	jle	 L(14)
	ALIGN_3

L(13):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A1)
#endif

	mulps	 %xmm8,  %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A2 (-32 * SIZE, A1, LDA,  1, %xmm4)
	mulps	 %xmm8,  %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A2 (-28 * SIZE, A1, LDA,  1, %xmm5)
	mulps	 %xmm8,  %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_A2 (-24 * SIZE, A1, LDA,  1, %xmm6)
	mulps	 %xmm8,  %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_A2 (-20 * SIZE, A1, LDA,  1, %xmm7)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A1, LDA)
#endif

	mulps	 %xmm9,  %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A2 (-32 * SIZE, A1, LDA,  2, %xmm4)
	mulps	 %xmm9,  %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A2 (-28 * SIZE, A1, LDA,  2, %xmm5)
	mulps	 %xmm9,  %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_A2 (-24 * SIZE, A1, LDA,  2, %xmm6)
	mulps	 %xmm9,  %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_A2 (-20 * SIZE, A1, LDA,  2, %xmm7)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A1, LDA, 2)
#endif

	mulps	 %xmm10, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A2 (-32 * SIZE, A1, LDA3, 1, %xmm4)
	mulps	 %xmm10, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A2 (-28 * SIZE, A1, LDA3, 1, %xmm5)
	mulps	 %xmm10, %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_A2 (-24 * SIZE, A1, LDA3, 1, %xmm6)
	mulps	 %xmm10, %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_A2 (-20 * SIZE, A1, LDA3, 1, %xmm7)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A1, LDA3)
#endif

	mulps	 %xmm11, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A1 (-32 * SIZE, A2, %xmm4)
	mulps	 %xmm11, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A1 (-28 * SIZE, A2, %xmm5)
	mulps	 %xmm11, %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_A1 (-24 * SIZE, A2, %xmm6)
	mulps	 %xmm11, %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_A1 (-20 * SIZE, A2, %xmm7)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A2)
#endif

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A2 (-32 * SIZE, A2, LDA,  1, %xmm4)
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A2 (-28 * SIZE, A2, LDA,  1, %xmm5)
	mulps	 %xmm12, %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_A2 (-24 * SIZE, A2, LDA,  1, %xmm6)
	mulps	 %xmm12, %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_A2 (-20 * SIZE, A2, LDA,  1, %xmm7)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A2, LDA)
#endif

	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A2 (-32 * SIZE, A2, LDA,  2, %xmm4)
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A2 (-28 * SIZE, A2, LDA,  2, %xmm5)
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_A2 (-24 * SIZE, A2, LDA,  2, %xmm6)
	mulps	 %xmm13, %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_A2 (-20 * SIZE, A2, LDA,  2, %xmm7)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A2, LDA, 2)
#endif

	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A2 (-32 * SIZE, A2, LDA3, 1, %xmm4)
	mulps	 %xmm14, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A2 (-28 * SIZE, A2, LDA3, 1, %xmm5)
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_A2 (-24 * SIZE, A2, LDA3, 1, %xmm6)
	mulps	 %xmm14, %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_A2 (-20 * SIZE, A2, LDA3, 1, %xmm7)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A2, LDA3)
#endif

	mulps	 %xmm15, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A1 (-16 * SIZE, A1, %xmm4)
	mulps	 %xmm15, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A1 (-12 * SIZE, A1, %xmm5)
	mulps	 %xmm15, %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_A1 ( -8 * SIZE, A1, %xmm6)
	mulps	 %xmm15, %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_A1 ( -4 * SIZE, A1, %xmm7)

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE) - 128 + PREOFFSET(Y1)
#endif

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	MOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-12 * SIZE, Y1, %xmm1)
	MOVUPS_YL1( -8 * SIZE, Y1, %xmm2)
	MOVUPS_YL1( -4 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(13)
	ALIGN_3

L(14):
	mulps	 %xmm8,  %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A2 (-32 * SIZE, A1, LDA,  1, %xmm4)
	mulps	 %xmm8,  %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A2 (-28 * SIZE, A1, LDA,  1, %xmm5)
	mulps	 %xmm8,  %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_A2 (-24 * SIZE, A1, LDA,  1, %xmm6)
	mulps	 %xmm8,  %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_A2 (-20 * SIZE, A1, LDA,  1, %xmm7)

	mulps	 %xmm9,  %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A2 (-32 * SIZE, A1, LDA,  2, %xmm4)
	mulps	 %xmm9,  %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A2 (-28 * SIZE, A1, LDA,  2, %xmm5)
	mulps	 %xmm9,  %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_A2 (-24 * SIZE, A1, LDA,  2, %xmm6)
	mulps	 %xmm9,  %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_A2 (-20 * SIZE, A1, LDA,  2, %xmm7)

	mulps	 %xmm10, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A2 (-32 * SIZE, A1, LDA3, 1, %xmm4)
	mulps	 %xmm10, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A2 (-28 * SIZE, A1, LDA3, 1, %xmm5)
	mulps	 %xmm10, %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_A2 (-24 * SIZE, A1, LDA3, 1, %xmm6)
	mulps	 %xmm10, %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_A2 (-20 * SIZE, A1, LDA3, 1, %xmm7)

	mulps	 %xmm11, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A1 (-32 * SIZE, A2, %xmm4)
	mulps	 %xmm11, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A1 (-28 * SIZE, A2, %xmm5)
	mulps	 %xmm11, %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_A1 (-24 * SIZE, A2, %xmm6)
	mulps	 %xmm11, %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_A1 (-20 * SIZE, A2, %xmm7)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A2 (-32 * SIZE, A2, LDA,  1, %xmm4)
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A2 (-28 * SIZE, A2, LDA,  1, %xmm5)
	mulps	 %xmm12, %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_A2 (-24 * SIZE, A2, LDA,  1, %xmm6)
	mulps	 %xmm12, %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_A2 (-20 * SIZE, A2, LDA,  1, %xmm7)

	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A2 (-32 * SIZE, A2, LDA,  2, %xmm4)
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A2 (-28 * SIZE, A2, LDA,  2, %xmm5)
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_A2 (-24 * SIZE, A2, LDA,  2, %xmm6)
	mulps	 %xmm13, %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_A2 (-20 * SIZE, A2, LDA,  2, %xmm7)

	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A2 (-32 * SIZE, A2, LDA3, 1, %xmm4)
	mulps	 %xmm14, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A2 (-28 * SIZE, A2, LDA3, 1, %xmm5)
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_A2 (-24 * SIZE, A2, LDA3, 1, %xmm6)
	mulps	 %xmm14, %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_A2 (-20 * SIZE, A2, LDA3, 1, %xmm7)

	mulps	 %xmm15, %xmm4
	addps	 %xmm4,  %xmm0
	mulps	 %xmm15, %xmm5
	addps	 %xmm5,  %xmm1
	mulps	 %xmm15, %xmm6
	addps	 %xmm6,  %xmm2
	mulps	 %xmm15, %xmm7
	addps	 %xmm7,  %xmm3

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1
	ALIGN_3

L(15):
	testq	$8, MM
	je	L(16)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm4)
	MOVUPS_A1 (-28 * SIZE, A1, %xmm5)

	MOVUPS_A2 (-32 * SIZE, A1, LDA,  1, %xmm6)
	MOVUPS_A2 (-28 * SIZE, A1, LDA,  1, %xmm7)

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)

	mulps	 %xmm8,  %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A2 (-32 * SIZE, A1, LDA,  2, %xmm4)
	mulps	 %xmm8,  %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A2 (-28 * SIZE, A1, LDA,  2, %xmm5)

	mulps	 %xmm9,  %xmm6
	addps	 %xmm6,  %xmm0
	MOVUPS_A2 (-32 * SIZE, A1, LDA3, 1, %xmm6)
	mulps	 %xmm9,  %xmm7
	addps	 %xmm7,  %xmm1
	MOVUPS_A2 (-28 * SIZE, A1, LDA3, 1, %xmm7)

	mulps	 %xmm10, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A1 (-32 * SIZE, A2, %xmm4)
	mulps	 %xmm10, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A1 (-28 * SIZE, A2, %xmm5)

	mulps	 %xmm11, %xmm6
	addps	 %xmm6,  %xmm0
	MOVUPS_A2 (-32 * SIZE, A2, LDA,  1, %xmm6)
	mulps	 %xmm11, %xmm7
	addps	 %xmm7,  %xmm1
	MOVUPS_A2 (-28 * SIZE, A2, LDA,  1, %xmm7)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A2 (-32 * SIZE, A2, LDA,  2, %xmm4)
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A2 (-28 * SIZE, A2, LDA,  2, %xmm5)

	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm0
	MOVUPS_A2 (-32 * SIZE, A2, LDA3, 1, %xmm6)
	mulps	 %xmm13, %xmm7
	addps	 %xmm7,  %xmm1
	MOVUPS_A2 (-28 * SIZE, A2, LDA3, 1, %xmm7)

	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm0
	mulps	 %xmm14, %xmm5
	addps	 %xmm5,  %xmm1

	mulps	 %xmm15, %xmm6
	addps	 %xmm6,  %xmm0
	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	mulps	 %xmm15, %xmm7
	addps	 %xmm7,  %xmm1
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)

	addq	 $8 * SIZE, A1
	addq	 $8 * SIZE, A2
	addq	 $8 * SIZE, Y1
	ALIGN_3

L(16):
	testq	$4, MM
	je	L(17)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm4)
	MOVUPS_A2 (-32 * SIZE, A1, LDA,  1, %xmm5)
	MOVUPS_A2 (-32 * SIZE, A1, LDA,  2, %xmm6)
	MOVUPS_A2 (-32 * SIZE, A1, LDA3, 1, %xmm7)

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)

	mulps	 %xmm8,  %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A1 (-32 * SIZE, A2, %xmm4)
	mulps	 %xmm9,  %xmm5
	addps	 %xmm5,  %xmm0
	MOVUPS_A2 (-32 * SIZE, A2, LDA,  1, %xmm5)
	mulps	 %xmm10, %xmm6
	addps	 %xmm6,  %xmm0
	MOVUPS_A2 (-32 * SIZE, A2, LDA,  2, %xmm6)
	mulps	 %xmm11, %xmm7
	addps	 %xmm7,  %xmm0
	MOVUPS_A2 (-32 * SIZE, A2, LDA3, 1, %xmm7)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm0
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm0
	mulps	 %xmm15, %xmm7
	addps	 %xmm7,  %xmm0

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, A2
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(17):
	testq	$2, MM
	je	L(18)

	movsd	 -32 * SIZE(A1), %xmm4
	movsd	 -32 * SIZE(A1, LDA,  1), %xmm5
	movsd	 -32 * SIZE(A1, LDA,  2), %xmm6
	movsd	 -32 * SIZE(A1, LDA3, 1), %xmm7

	movsd	 -32 * SIZE(Y1), %xmm0

	mulps	 %xmm8,  %xmm4
	addps	 %xmm4,  %xmm0
	movsd	 -32 * SIZE(A2), %xmm4
	mulps	 %xmm9,  %xmm5
	addps	 %xmm5,  %xmm0
	movsd	 -32 * SIZE(A2, LDA,  1), %xmm5
	mulps	 %xmm10, %xmm6
	addps	 %xmm6,  %xmm0
	movsd	 -32 * SIZE(A2, LDA,  2), %xmm6
	mulps	 %xmm11, %xmm7
	addps	 %xmm7,  %xmm0
	movsd	 -32 * SIZE(A2, LDA3, 1), %xmm7

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm0
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm0
	mulps	 %xmm15, %xmm7
	addps	 %xmm7,  %xmm0

	movlps	 %xmm0, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(18):
	testq	$1, MM
	je	L(19)

	movss	 -32 * SIZE(A1), %xmm4
	movss	 -32 * SIZE(A1, LDA,  1), %xmm5
	movss	 -32 * SIZE(A1, LDA,  2), %xmm6
	movss	 -32 * SIZE(A1, LDA3, 1), %xmm7

	movss	 -32 * SIZE(Y1), %xmm0

	mulss	 %xmm8,  %xmm4
	addss	 %xmm4,  %xmm0
	movss	 -32 * SIZE(A2), %xmm4
	mulss	 %xmm9,  %xmm5
	addss	 %xmm5,  %xmm0
	movss	 -32 * SIZE(A2, LDA,  1), %xmm5
	mulss	 %xmm10, %xmm6
	addss	 %xmm6,  %xmm0
	movss	 -32 * SIZE(A2, LDA,  2), %xmm6
	mulss	 %xmm11, %xmm7
	addss	 %xmm7,  %xmm0
	movss	 -32 * SIZE(A2, LDA3, 1), %xmm7

	mulss	 %xmm12, %xmm4
	addss	 %xmm4,  %xmm0
	mulss	 %xmm13, %xmm5
	addss	 %xmm5,  %xmm0
	mulss	 %xmm14, %xmm6
	addss	 %xmm6,  %xmm0
	mulss	 %xmm15, %xmm7
	addss	 %xmm7,  %xmm0

	movss	 %xmm0, -32 * SIZE(Y1)
	ALIGN_3

L(19):
	cmpq	$8, N
	jge	L(11)
	ALIGN_3

L(20):
#endif

	cmpq	$4, N
	jl	L(30)

#if GEMV_UNROLL == 4
	ALIGN_3

L(21):
#endif
	subq	$4, N

	leaq	32 * SIZE(BUFFER), Y1
	movq	A,  A1
	leaq	(A,  LDA, 2), A2
	leaq	(A,  LDA, 4), A

	movss	(X), %xmm12
	addq	INCX, X
	movss	(X), %xmm13
	addq	INCX, X
	movss	(X), %xmm14
	addq	INCX, X
	movss	(X), %xmm15
	addq	INCX, X

	movss	ALPHA, %xmm0

	mulss	%xmm0, %xmm12
	mulss	%xmm0, %xmm13
	mulss	%xmm0, %xmm14
	mulss	%xmm0, %xmm15

	shufps	$0, %xmm12, %xmm12
	shufps	$0, %xmm13, %xmm13
	shufps	$0, %xmm14, %xmm14
	shufps	$0, %xmm15, %xmm15

#ifdef ALIGNED_ACCESS
	cmpq	$3, M
	jle	L(27)

	testq	$SIZE, A1
	je	L(2X)

	movss	 -32 * SIZE(A1), %xmm0
	movss	 -32 * SIZE(A1, LDA), %xmm1
	movss	 -32 * SIZE(A2), %xmm2
	movss	 -32 * SIZE(A2, LDA), %xmm3

	movss	 -32 * SIZE(Y1), %xmm8

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8
	mulss	 %xmm13, %xmm1
	addss	 %xmm1,  %xmm8
	mulss	 %xmm14, %xmm2
	addss	 %xmm2, %xmm8
	mulss	 %xmm15, %xmm3
	addss	 %xmm3, %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, Y1
	ALIGN_3

L(2X):
	testq	$2 * SIZE, A1
	je	L(2XX)

	movsd	 -32 * SIZE(A1), %xmm0
	movsd	 -32 * SIZE(A1, LDA), %xmm1
	movsd	 -32 * SIZE(A2), %xmm2
	movsd	 -32 * SIZE(A2, LDA), %xmm3

	movsd	 -32 * SIZE(Y1), %xmm8

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8
	mulps	 %xmm13, %xmm1
	addps	 %xmm1,  %xmm8
	mulps	 %xmm14, %xmm2
	addps	 %xmm2, %xmm8
	mulps	 %xmm15, %xmm3
	addps	 %xmm3, %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(2XX):
#endif

	movq	MM,  I
	sarq	$4,  I
	jle	L(25)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm8)
	MOVUPS_A1 (-28 * SIZE, A1, %xmm9)
	MOVUPS_A1 (-24 * SIZE, A1, %xmm10)
	MOVUPS_A1 (-20 * SIZE, A1, %xmm11)

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YL1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YL1(-20 * SIZE, Y1, %xmm3)

	MOVUPS_A2(-32 * SIZE, A1, LDA, 1, %xmm4)
	MOVUPS_A2(-28 * SIZE, A1, LDA, 1, %xmm5)
	MOVUPS_A2(-24 * SIZE, A1, LDA, 1, %xmm6)
	MOVUPS_A2(-20 * SIZE, A1, LDA, 1, %xmm7)

	decq	 I
	jle	 L(24)
	ALIGN_3

L(23):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1)
#endif

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_A1(-32 * SIZE, A2, %xmm8)
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	MOVUPS_A1(-28 * SIZE, A2, %xmm9)
	mulps	 %xmm12, %xmm10
	addps	 %xmm10,  %xmm2
	MOVUPS_A1(-24 * SIZE, A2, %xmm10)
	mulps	 %xmm12, %xmm11
	addps	 %xmm11,  %xmm3
	MOVUPS_A1(-20 * SIZE, A2, %xmm11)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1, LDA)
#endif

	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A2(-32 * SIZE, A2, LDA, 1, %xmm4)
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A2(-28 * SIZE, A2, LDA, 1, %xmm5)
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_A2(-24 * SIZE, A2, LDA, 1, %xmm6)
	mulps	 %xmm13, %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_A2(-20 * SIZE, A2, LDA, 1, %xmm7)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2)
#endif

	mulps	 %xmm14, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_A1 (-16 * SIZE, A1, %xmm8)
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm1
	MOVUPS_A1 (-12 * SIZE, A1, %xmm9)
	mulps	 %xmm14, %xmm10
	addps	 %xmm10,  %xmm2
	MOVUPS_A1 ( -8 * SIZE, A1, %xmm10)
	mulps	 %xmm14, %xmm11
	addps	 %xmm11,  %xmm3
	MOVUPS_A1 ( -4 * SIZE, A1, %xmm11)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2, LDA)
#endif

	mulps	 %xmm15, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A2(-16 * SIZE, A1, LDA, 1, %xmm4)
	mulps	 %xmm15, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A2(-12 * SIZE, A1, LDA, 1, %xmm5)
	mulps	 %xmm15, %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_A2( -8 * SIZE, A1, LDA, 1, %xmm6)
	mulps	 %xmm15, %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_A2( -4 * SIZE, A1, LDA, 1, %xmm7)

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(Y1)
#endif

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	MOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-12 * SIZE, Y1, %xmm1)
	MOVUPS_YL1( -8 * SIZE, Y1, %xmm2)
	MOVUPS_YL1( -4 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(23)
	ALIGN_3

L(24):
	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_A1(-32 * SIZE, A2, %xmm8)
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	MOVUPS_A1(-28 * SIZE, A2, %xmm9)
	mulps	 %xmm12, %xmm10
	addps	 %xmm10,  %xmm2
	MOVUPS_A1(-24 * SIZE, A2, %xmm10)
	mulps	 %xmm12, %xmm11
	addps	 %xmm11,  %xmm3
	MOVUPS_A1(-20 * SIZE, A2, %xmm11)

	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A2(-32 * SIZE, A2, LDA, 1, %xmm4)
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A2(-28 * SIZE, A2, LDA, 1, %xmm5)
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_A2(-24 * SIZE, A2, LDA, 1, %xmm6)
	mulps	 %xmm13, %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_A2(-20 * SIZE, A2, LDA, 1, %xmm7)

	mulps	 %xmm14, %xmm8
	addps	 %xmm8,  %xmm0
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm1
	mulps	 %xmm14, %xmm10
	addps	 %xmm10,  %xmm2
	mulps	 %xmm14, %xmm11
	addps	 %xmm11,  %xmm3

	mulps	 %xmm15, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	mulps	 %xmm15, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	mulps	 %xmm15, %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	mulps	 %xmm15, %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1
	ALIGN_3

L(25):
	testq	$8, MM
	je	L(26)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm8)
	MOVUPS_A1 (-28 * SIZE, A1, %xmm9)

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)

	MOVUPS_A2(-32 * SIZE, A1, LDA, 1, %xmm4)
	MOVUPS_A2(-28 * SIZE, A1, LDA, 1, %xmm5)

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_A1(-32 * SIZE, A2, %xmm8)
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	MOVUPS_A1(-28 * SIZE, A2, %xmm9)

	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A2(-32 * SIZE, A2, LDA, 1, %xmm4)
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A2(-28 * SIZE, A2, LDA, 1, %xmm5)
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm2

	mulps	 %xmm14, %xmm8
	addps	 %xmm8,  %xmm0
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm1

	mulps	 %xmm15, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	mulps	 %xmm15, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)

	addq	 $8 * SIZE, A1
	addq	 $8 * SIZE, A2
	addq	 $8 * SIZE, Y1
	ALIGN_3

L(26):
	testq	$4, MM
	je	L(27)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm8)

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)

	MOVUPS_A2(-32 * SIZE, A1, LDA, 1, %xmm4)

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_A1(-32 * SIZE, A2, %xmm8)

	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A2(-32 * SIZE, A2, LDA, 1, %xmm4)

	mulps	 %xmm14, %xmm8
	addps	 %xmm8,  %xmm0

	mulps	 %xmm15, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, A2
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(27):
	testq	$2, MM
	je	L(28)

	movsd	 -32 * SIZE(A1), %xmm0
	movsd	 -32 * SIZE(A1, LDA), %xmm1
	movsd	 -32 * SIZE(A2), %xmm2
	movsd	 -32 * SIZE(A2, LDA), %xmm3

	movsd	 -32 * SIZE(Y1), %xmm8

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8
	mulps	 %xmm13, %xmm1
	addps	 %xmm1,  %xmm8
	mulps	 %xmm14, %xmm2
	addps	 %xmm2, %xmm8
	mulps	 %xmm15, %xmm3
	addps	 %xmm3, %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(28):
	testq	$1, MM
#if GEMV_UNROLL == 4
	je	L(29)
#else
	je	L(30)
#endif

	movss	 -32 * SIZE(Y1), %xmm8

	movss	 -32 * SIZE(A1), %xmm0
	movss	 -32 * SIZE(A1, LDA), %xmm1
	movss	 -32 * SIZE(A2), %xmm2
	movss	 -32 * SIZE(A2, LDA), %xmm3

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8
	mulss	 %xmm13, %xmm1
	addss	 %xmm1,  %xmm8
	mulss	 %xmm14, %xmm2
	addss	 %xmm2, %xmm8
	mulss	 %xmm15, %xmm3
	addss	 %xmm3, %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)
	ALIGN_3

#if GEMV_UNROLL == 4
L(29):
	cmpq	$4, N
	jge	L(21)
#endif
	ALIGN_3

L(30):
	testq	N, N
	jle	L(990)

	cmpq	$3, N
	jne	L(40)

	leaq	32 * SIZE(BUFFER), Y1
	movq	A,  A1
	leaq	(A,  LDA, 2), A2

	movss	(X), %xmm12
	addq	INCX, X
	movss	(X), %xmm13
	addq	INCX, X
	movss	(X), %xmm14
	addq	INCX, X

	movss	ALPHA, %xmm0

	mulss	%xmm0, %xmm12
	mulss	%xmm0, %xmm13
	mulss	%xmm0, %xmm14

	shufps	$0, %xmm12, %xmm12
	shufps	$0, %xmm13, %xmm13
	shufps	$0, %xmm14, %xmm14

#ifdef ALIGNED_ACCESS
	cmpq	$3, M
	jle	L(37)

	testq	$SIZE, A1
	je	L(3X)

	movss	 -32 * SIZE(A1), %xmm0
	movss	 -32 * SIZE(A1, LDA), %xmm1
	movss	 -32 * SIZE(A2), %xmm2

	movss	 -32 * SIZE(Y1), %xmm8

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8
	mulss	 %xmm13, %xmm1
	addss	 %xmm1,  %xmm8
	mulss	 %xmm14, %xmm2
	addss	 %xmm2, %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, Y1
	ALIGN_3

L(3X):
	testq	$2 * SIZE, A1
	je	L(3XX)

	movsd	 -32 * SIZE(A1), %xmm0
	movsd	 -32 * SIZE(A1, LDA), %xmm1
	movsd	 -32 * SIZE(A2), %xmm2

	movsd	 -32 * SIZE(Y1), %xmm8

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8
	mulps	 %xmm13, %xmm1
	addps	 %xmm1,  %xmm8
	mulps	 %xmm14, %xmm2
	addps	 %xmm2, %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(3XX):
#endif

	movq	MM,  I
	sarq	$4,  I
	jle	L(35)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm8)
	MOVUPS_A1 (-28 * SIZE, A1, %xmm9)
	MOVUPS_A1 (-24 * SIZE, A1, %xmm10)
	MOVUPS_A1 (-20 * SIZE, A1, %xmm11)

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YL1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YL1(-20 * SIZE, Y1, %xmm3)

	MOVUPS_A2(-32 * SIZE, A1, LDA, 1, %xmm4)
	MOVUPS_A2(-28 * SIZE, A1, LDA, 1, %xmm5)
	MOVUPS_A2(-24 * SIZE, A1, LDA, 1, %xmm6)
	MOVUPS_A2(-20 * SIZE, A1, LDA, 1, %xmm7)

	decq	 I
	jle	 L(34)
	ALIGN_3

L(33):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A1)
#endif

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_A1(-32 * SIZE, A2, %xmm8)
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	MOVUPS_A1(-28 * SIZE, A2, %xmm9)
	mulps	 %xmm12, %xmm10
	addps	 %xmm10,  %xmm2
	MOVUPS_A1(-24 * SIZE, A2, %xmm10)
	mulps	 %xmm12, %xmm11
	addps	 %xmm11,  %xmm3
	MOVUPS_A1(-20 * SIZE, A2, %xmm11)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A1, LDA)
#endif

	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A2(-16 * SIZE, A1, LDA, 1, %xmm4)
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A2(-12 * SIZE, A1, LDA, 1, %xmm5)
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_A2( -8 * SIZE, A1, LDA, 1, %xmm6)
	mulps	 %xmm13, %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_A2( -4 * SIZE, A1, LDA, 1, %xmm7)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A2)
#endif

	mulps	 %xmm14, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_A1 (-16 * SIZE, A1, %xmm8)
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm1
	MOVUPS_A1 (-12 * SIZE, A1, %xmm9)
	mulps	 %xmm14, %xmm10
	addps	 %xmm10,  %xmm2
	MOVUPS_A1 ( -8 * SIZE, A1, %xmm10)
	mulps	 %xmm14, %xmm11
	addps	 %xmm11,  %xmm3
	MOVUPS_A1 ( -4 * SIZE, A1, %xmm11)

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(Y1)
#endif

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	MOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-12 * SIZE, Y1, %xmm1)
	MOVUPS_YL1( -8 * SIZE, Y1, %xmm2)
	MOVUPS_YL1( -4 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(33)
	ALIGN_3

L(34):
	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_A1(-32 * SIZE, A2, %xmm8)
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	MOVUPS_A1(-28 * SIZE, A2, %xmm9)
	mulps	 %xmm12, %xmm10
	addps	 %xmm10,  %xmm2
	MOVUPS_A1(-24 * SIZE, A2, %xmm10)
	mulps	 %xmm12, %xmm11
	addps	 %xmm11,  %xmm3
	MOVUPS_A1(-20 * SIZE, A2, %xmm11)

	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm0
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm1
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm2
	mulps	 %xmm13, %xmm7
	addps	 %xmm7,  %xmm3

	mulps	 %xmm14, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm1
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	mulps	 %xmm14, %xmm10
	addps	 %xmm10,  %xmm2
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	mulps	 %xmm14, %xmm11
	addps	 %xmm11,  %xmm3
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1
	ALIGN_3

L(35):
	testq	$8, MM
	je	L(36)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm8)
	MOVUPS_A1 (-28 * SIZE, A1, %xmm9)

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)

	MOVUPS_A2(-32 * SIZE, A1, LDA, 1, %xmm4)
	MOVUPS_A2(-28 * SIZE, A1, LDA, 1, %xmm5)

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_A1(-32 * SIZE, A2, %xmm8)
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	MOVUPS_A1(-28 * SIZE, A2, %xmm9)

	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm0
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm1
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm2

	mulps	 %xmm14, %xmm8
	addps	 %xmm8,  %xmm0
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm1

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)

	addq	 $8 * SIZE, A1
	addq	 $8 * SIZE, A2
	addq	 $8 * SIZE, Y1
	ALIGN_3

L(36):
	testq	$4, MM
	je	L(37)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm8)

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)

	MOVUPS_A2(-32 * SIZE, A1, LDA, 1, %xmm4)

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_A1(-32 * SIZE, A2, %xmm8)

	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm0

	mulps	 %xmm14, %xmm8
	addps	 %xmm8,  %xmm0

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, A2
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(37):
	testq	$2, MM
	je	L(38)

	movsd	 -32 * SIZE(A1), %xmm0
	movsd	 -32 * SIZE(A1, LDA), %xmm1
	movsd	 -32 * SIZE(A2), %xmm2

	movsd	 -32 * SIZE(Y1), %xmm8

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8
	mulps	 %xmm13, %xmm1
	addps	 %xmm1,  %xmm8
	mulps	 %xmm14, %xmm2
	addps	 %xmm2, %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(38):
	testq	$1, MM
	je	L(990)

	movss	 -32 * SIZE(Y1), %xmm8

	movss	 -32 * SIZE(A1), %xmm0
	movss	 -32 * SIZE(A1, LDA), %xmm1
	movss	 -32 * SIZE(A2), %xmm2

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8
	mulss	 %xmm13, %xmm1
	addss	 %xmm1,  %xmm8
	mulss	 %xmm14, %xmm2
	addss	 %xmm2, %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)
	jmp	 L(990)
	ALIGN_3

L(40):
	cmpq	$2, N
	jne	L(50)

	leaq	32 * SIZE(BUFFER), Y1
	movq	A,  A1
	leaq	(A,  LDA, 1), A2
	leaq	(A,  LDA, 2), A

	movss	(X), %xmm12
	addq	INCX, X
	movss	(X), %xmm13
	addq	INCX, X

	movss	ALPHA, %xmm0

	mulss	%xmm0, %xmm12
	mulss	%xmm0, %xmm13

	shufps	$0, %xmm12, %xmm12
	shufps	$0, %xmm13, %xmm13

#ifdef ALIGNED_ACCESS
	cmpq	$3, M
	jle	L(47)

	testq	$SIZE, A1
	je	L(4X)

	movss	 -32 * SIZE(Y1), %xmm8

	movss	 -32 * SIZE(A1), %xmm0
	movss	 -32 * SIZE(A2), %xmm1

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8
	mulss	 %xmm13, %xmm1
	addss	 %xmm1,  %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, Y1
	ALIGN_3

L(4X):
	testq	$2 * SIZE, A1
	je	L(4XX)

	movsd	 -32 * SIZE(A1), %xmm0
	movsd	 -32 * SIZE(A2), %xmm1

	movsd	 -32 * SIZE(Y1), %xmm8

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8
	mulps	 %xmm13, %xmm1
	addps	 %xmm1,  %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(4XX):
#endif

	movq	MM,  I
	sarq	$4,  I
	jle	L(45)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm8)
	MOVUPS_A1 (-28 * SIZE, A1, %xmm9)
	MOVUPS_A1 (-24 * SIZE, A1, %xmm10)
	MOVUPS_A1 (-20 * SIZE, A1, %xmm11)

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YL1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YL1(-20 * SIZE, Y1, %xmm3)

	MOVUPS_A1(-32 * SIZE, A2, %xmm4)
	MOVUPS_A1(-28 * SIZE, A2, %xmm5)
	MOVUPS_A1(-24 * SIZE, A2, %xmm6)
	MOVUPS_A1(-20 * SIZE, A2, %xmm7)

	decq	 I
	jle	 L(44)
	ALIGN_3

L(43):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A1)
#endif

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_A1 (-16 * SIZE, A1, %xmm8)
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	MOVUPS_A1 (-12 * SIZE, A1, %xmm9)
	mulps	 %xmm12, %xmm10
	addps	 %xmm10,  %xmm2
	MOVUPS_A1 ( -8 * SIZE, A1, %xmm10)
	mulps	 %xmm12, %xmm11
	addps	 %xmm11,  %xmm3
	MOVUPS_A1 ( -4 * SIZE, A1, %xmm11)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A2)
#endif

	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_A1(-16 * SIZE, A2, %xmm4)
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_A1(-12 * SIZE, A2, %xmm5)
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_A1( -8 * SIZE, A2, %xmm6)
	mulps	 %xmm13, %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_A1( -4 * SIZE, A2, %xmm7)

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(Y1)
#endif

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	MOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-12 * SIZE, Y1, %xmm1)
	MOVUPS_YL1( -8 * SIZE, Y1, %xmm2)
	MOVUPS_YL1( -4 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(43)
	ALIGN_3

L(44):
	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	mulps	 %xmm12, %xmm10
	addps	 %xmm10,  %xmm2
	mulps	 %xmm12, %xmm11
	addps	 %xmm11,  %xmm3

	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm2
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	mulps	 %xmm13, %xmm7
	addps	 %xmm7,  %xmm3
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1
	ALIGN_3

L(45):
	testq	$8, MM
	je	L(46)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm8)
	MOVUPS_A1 (-28 * SIZE, A1, %xmm9)

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_A1(-32 * SIZE, A2, %xmm4)
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	MOVUPS_A1(-28 * SIZE, A2, %xmm5)

	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm1
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)

	addq	 $8 * SIZE, A1
	addq	 $8 * SIZE, A2
	addq	 $8 * SIZE, Y1
	ALIGN_3

L(46):
	testq	$4, MM
	je	L(47)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm8)

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_A1(-32 * SIZE, A2, %xmm4)
	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm0
	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, A2
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(47):
	testq	$2, MM
	je	L(48)

	movsd	 -32 * SIZE(A1), %xmm0
	movsd	 -32 * SIZE(A2), %xmm1

	movsd	 -32 * SIZE(Y1), %xmm8

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8
	mulps	 %xmm13, %xmm1
	addps	 %xmm1,  %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(48):
	testq	$1, MM
	je	L(990)

	movss	 -32 * SIZE(Y1), %xmm8

	movss	 -32 * SIZE(A1), %xmm0
	movss	 -32 * SIZE(A2), %xmm1

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8
	mulss	 %xmm13, %xmm1
	addss	 %xmm1,  %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)
	jmp	 L(990)
	ALIGN_3

L(50):
	cmpq	$1, N
	jne	L(990)

	leaq	32 * SIZE(BUFFER), Y1
	movq	A,  A1

	movss	(X), %xmm12

	mulss	ALPHA, %xmm12
	shufps	$0, %xmm12, %xmm12

#ifdef ALIGNED_ACCESS
	cmpq	$3, M
	jle	L(57)

	testq	$SIZE, A1
	je	L(5X)

	movss	 -32 * SIZE(Y1), %xmm8
	movss	 -32 * SIZE(A1), %xmm0

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, Y1
	ALIGN_3

L(5X):
	testq	$2 * SIZE, A1
	je	L(5XX)

	movsd	 -32 * SIZE(Y1), %xmm8
	movsd	 -32 * SIZE(A1), %xmm0

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(5XX):
#endif

	movq	MM,  I
	sarq	$4,  I
	jle	L(55)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm8)
	MOVUPS_A1 (-28 * SIZE, A1, %xmm9)
	MOVUPS_A1 (-24 * SIZE, A1, %xmm10)
	MOVUPS_A1 (-20 * SIZE, A1, %xmm11)

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YL1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YL1(-20 * SIZE, Y1, %xmm3)

	decq	 I
	jle	 L(54)
	ALIGN_3

L(53):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(A1)
#endif

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_A1 (-16 * SIZE, A1, %xmm8)
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	MOVUPS_A1 (-12 * SIZE, A1, %xmm9)
	mulps	 %xmm12, %xmm10
	addps	 %xmm10,  %xmm2
	MOVUPS_A1 ( -8 * SIZE, A1, %xmm10)
	mulps	 %xmm12, %xmm11
	addps	 %xmm11,  %xmm3
	MOVUPS_A1 ( -4 * SIZE, A1, %xmm11)

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(Y1)
#endif

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	MOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-12 * SIZE, Y1, %xmm1)
	MOVUPS_YL1( -8 * SIZE, Y1, %xmm2)
	MOVUPS_YL1( -4 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(53)
	ALIGN_3

L(54):
	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	mulps	 %xmm12, %xmm10
	addps	 %xmm10,  %xmm2
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	mulps	 %xmm12, %xmm11
	addps	 %xmm11,  %xmm3
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, Y1
	ALIGN_3

L(55):
	testq	$8, MM
	je	L(56)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm8)
	MOVUPS_A1 (-28 * SIZE, A1, %xmm9)

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)

	addq	 $8 * SIZE, A1
	addq	 $8 * SIZE, Y1
	ALIGN_3

L(56):
	testq	$4, MM
	je	L(57)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm8)
	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(57):
	testq	$2, MM
	je	L(58)

	movsd	 -32 * SIZE(Y1), %xmm8
	movsd	 -32 * SIZE(A1), %xmm0

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(58):
	testq	$1, MM
	je	L(990)

	movss	 -32 * SIZE(Y1), %xmm8
	movss	 -32 * SIZE(A1), %xmm0

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)

#ifdef ALIGNED_ACCESS
       jmp	L(990)
       ALIGN_3

L(100):
	testq	$2 * SIZE - 1, LDA
	jne	L(200)

	cmpq	$4, N
	jl	L(110)
	ALIGN_3

L(101):
	subq	$4, N

	leaq	32 * SIZE(BUFFER), Y1
	movq	A,  A1
	leaq	(A,  LDA, 2), A2
	leaq	(A,  LDA, 4), A

	movss	(X), %xmm12
	addq	INCX, X
	movss	(X), %xmm13
	addq	INCX, X
	movss	(X), %xmm14
	addq	INCX, X
	movss	(X), %xmm15
	addq	INCX, X

	movss	ALPHA, %xmm0

	mulss	%xmm0, %xmm12
	mulss	%xmm0, %xmm13
	mulss	%xmm0, %xmm14
	mulss	%xmm0, %xmm15

	shufps	$0, %xmm12, %xmm12
	shufps	$0, %xmm13, %xmm13
	shufps	$0, %xmm14, %xmm14
	shufps	$0, %xmm15, %xmm15

	cmpq	$3, M
	jle	L(107)

	testq	$SIZE, A1
	je	L(10X)

	movss	 -32 * SIZE(A1), %xmm0
	movss	 -32 * SIZE(A1, LDA), %xmm1
	movss	 -32 * SIZE(A2), %xmm2
	movss	 -32 * SIZE(A2, LDA), %xmm3

	movss	 -32 * SIZE(Y1), %xmm8

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8
	mulss	 %xmm13, %xmm1
	addss	 %xmm1,  %xmm8
	mulss	 %xmm14, %xmm2
	addss	 %xmm2, %xmm8
	mulss	 %xmm15, %xmm3
	addss	 %xmm3, %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, Y1
	ALIGN_3

L(10X):
	testq	$2 * SIZE, A1
	je	L(10XX)

	movsd	 -32 * SIZE(A1), %xmm0
	movsd	 -32 * SIZE(A1, LDA), %xmm1
	movsd	 -32 * SIZE(A2), %xmm2
	movsd	 -32 * SIZE(A2, LDA), %xmm3

	movsd	 -32 * SIZE(Y1), %xmm8

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8
	mulps	 %xmm13, %xmm1
	addps	 %xmm1,  %xmm8
	mulps	 %xmm14, %xmm2
	addps	 %xmm2, %xmm8
	mulps	 %xmm15, %xmm3
	addps	 %xmm3, %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(10XX):
	movhps	-32 * SIZE(A1, LDA), %xmm8
	movhps	-32 * SIZE(A2, LDA), %xmm9

	movq	MM,  I
	sarq	$4,  I
	jle	L(105)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -28 * SIZE(A1), %xmm5
	movaps	 -24 * SIZE(A1), %xmm6

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YL1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YL1(-20 * SIZE, Y1, %xmm3)

	decq	 I
	jle	 L(104)
	ALIGN_3

L(103):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1)
#endif

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A1), %xmm7
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -30 * SIZE(A1, LDA), %xmm4
	mulps	 %xmm12, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -26 * SIZE(A1, LDA), %xmm5
	mulps	 %xmm12, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	 -22 * SIZE(A1, LDA), %xmm6

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1, LDA)
#endif

	shufps	 $0x4e, %xmm4, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0
	movaps	 -18 * SIZE(A1, LDA), %xmm8
	shufps	 $0x4e, %xmm5, %xmm4
	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -32 * SIZE(A2), %xmm4
	shufps	 $0x4e, %xmm6, %xmm5
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -28 * SIZE(A2), %xmm5
	shufps	 $0x4e, %xmm8, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	 -24 * SIZE(A2), %xmm6

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2)
#endif

	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A2), %xmm7
	mulps	 %xmm14, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -30 * SIZE(A2, LDA), %xmm4
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -26 * SIZE(A2, LDA), %xmm5
	mulps	 %xmm14, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	 -22 * SIZE(A2, LDA), %xmm6

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2, LDA)
#endif

	shufps	 $0x4e, %xmm4, %xmm9
	mulps	 %xmm15, %xmm9
	addps	 %xmm9,  %xmm0
	movaps	 -18 * SIZE(A2, LDA), %xmm9
	shufps	 $0x4e, %xmm5, %xmm4
	mulps	 %xmm15, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -16 * SIZE(A1), %xmm4
	shufps	 $0x4e, %xmm6, %xmm5
	mulps	 %xmm15, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -12 * SIZE(A1), %xmm5
	shufps	 $0x4e, %xmm9, %xmm6
	mulps	 %xmm15, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	  -8 * SIZE(A1), %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(Y1)
#endif

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	MOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-12 * SIZE, Y1, %xmm1)
	MOVUPS_YL1( -8 * SIZE, Y1, %xmm2)
	MOVUPS_YL1( -4 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(103)
	ALIGN_3

L(104):
	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A1), %xmm7
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -30 * SIZE(A1, LDA), %xmm4
	mulps	 %xmm12, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -26 * SIZE(A1, LDA), %xmm5
	mulps	 %xmm12, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	 -22 * SIZE(A1, LDA), %xmm6

	shufps	 $0x4e, %xmm4, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0
	movaps	 -18 * SIZE(A1, LDA), %xmm8
	shufps	 $0x4e, %xmm5, %xmm4
	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -32 * SIZE(A2), %xmm4
	shufps	 $0x4e, %xmm6, %xmm5
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -28 * SIZE(A2), %xmm5
	shufps	 $0x4e, %xmm8, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	 -24 * SIZE(A2), %xmm6

	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A2), %xmm7
	mulps	 %xmm14, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -30 * SIZE(A2, LDA), %xmm4
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -26 * SIZE(A2, LDA), %xmm5
	mulps	 %xmm14, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	 -22 * SIZE(A2, LDA), %xmm6

	shufps	 $0x4e, %xmm4, %xmm9
	mulps	 %xmm15, %xmm9
	addps	 %xmm9,  %xmm0
	movaps	 -18 * SIZE(A2, LDA), %xmm9
	shufps	 $0x4e, %xmm5, %xmm4
	mulps	 %xmm15, %xmm4
	addps	 %xmm4,  %xmm1
	shufps	 $0x4e, %xmm6, %xmm5
	mulps	 %xmm15, %xmm5
	addps	 %xmm5,  %xmm2
	shufps	 $0x4e, %xmm9, %xmm6
	mulps	 %xmm15, %xmm6
	addps	 %xmm6,  %xmm3

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1
	ALIGN_3

L(105):
	testq	$8, MM
	je	L(106)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -28 * SIZE(A1), %xmm5
	movaps	 -30 * SIZE(A1, LDA), %xmm6
	movaps	 -26 * SIZE(A1, LDA), %xmm7

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -32 * SIZE(A2), %xmm4
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -28 * SIZE(A2), %xmm10

	shufps	 $0x4e, %xmm6, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0
	movaps	 -30 * SIZE(A2, LDA), %xmm11
	shufps	 $0x4e, %xmm7, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm1
	movaps	 %xmm7, %xmm8
	movaps	 -26 * SIZE(A2, LDA), %xmm7

	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm0
	mulps	 %xmm14, %xmm10
	addps	 %xmm10,  %xmm1

	shufps	 $0x4e, %xmm11, %xmm9
	mulps	 %xmm15, %xmm9
	addps	 %xmm9,  %xmm0
	shufps	 $0x4e, %xmm7, %xmm11
	mulps	 %xmm15, %xmm11
	addps	 %xmm11,  %xmm1
	movaps	 %xmm7, %xmm9

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)

	addq	 $8 * SIZE, A1
	addq	 $8 * SIZE, A2
	addq	 $8 * SIZE, Y1
	ALIGN_3

L(106):
	testq	$4, MM
	je	L(107)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -30 * SIZE(A1, LDA), %xmm5
	movaps	 -32 * SIZE(A2), %xmm6
	movaps	 -30 * SIZE(A2, LDA), %xmm7

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	shufps	 $0x4e, %xmm5, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm0
	shufps	 $0x4e, %xmm7, %xmm9
	mulps	 %xmm15, %xmm9
	addps	 %xmm9,  %xmm0

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, A2
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(107):
	testq	$2, MM
	je	L(108)

	movsd	 -32 * SIZE(A1), %xmm4
	movsd	 -32 * SIZE(A1, LDA), %xmm5
	movsd	 -32 * SIZE(A2), %xmm6
	movsd	 -32 * SIZE(A2, LDA), %xmm7

	movsd	 -32 * SIZE(Y1), %xmm0

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm0
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm0
	mulps	 %xmm15, %xmm7
	addps	 %xmm7,  %xmm0

	movlps	 %xmm0, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(108):
	testq	$1, MM
	je	L(109)

	movss	 -32 * SIZE(Y1), %xmm0

	movss	 -32 * SIZE(A1), %xmm4
	movss	 -32 * SIZE(A1, LDA), %xmm5
	movss	 -32 * SIZE(A2), %xmm6
	movss	 -32 * SIZE(A2, LDA), %xmm7

	mulss	 %xmm12, %xmm4
	addss	 %xmm4,  %xmm0
	mulss	 %xmm13, %xmm5
	addss	 %xmm5,  %xmm0
	mulss	 %xmm14, %xmm6
	addss	 %xmm6,  %xmm0
	mulss	 %xmm15, %xmm7
	addss	 %xmm7,  %xmm0

	movss	 %xmm0, -32 * SIZE(Y1)
	ALIGN_3

L(109):
	cmpq	$4, N
	jge	L(101)
	ALIGN_3

L(110):
	testq	N, N
	jle	L(990)

	cmpq	$3, N
	jne	L(120)

	leaq	32 * SIZE(BUFFER), Y1
	movq	A,  A1
	leaq	(A,  LDA, 2), A2
	leaq	(A,  LDA, 4), A

	movss	(X), %xmm12
	addq	INCX, X
	movss	(X), %xmm13
	addq	INCX, X
	movss	(X), %xmm14
	addq	INCX, X

	movss	ALPHA, %xmm0

	mulss	%xmm0, %xmm12
	mulss	%xmm0, %xmm13
	mulss	%xmm0, %xmm14

	shufps	$0, %xmm12, %xmm12
	shufps	$0, %xmm13, %xmm13
	shufps	$0, %xmm14, %xmm14

	cmpq	$3, M
	jle	L(117)

	testq	$SIZE, A1
	je	L(11X)

	movss	 -32 * SIZE(A1), %xmm0
	movss	 -32 * SIZE(A1, LDA), %xmm1
	movss	 -32 * SIZE(A2), %xmm2

	movss	 -32 * SIZE(Y1), %xmm8

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8
	mulss	 %xmm13, %xmm1
	addss	 %xmm1,  %xmm8
	mulss	 %xmm14, %xmm2
	addss	 %xmm2, %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, Y1
	ALIGN_3

L(11X):
	testq	$2 * SIZE, A1
	je	L(11XX)

	movsd	 -32 * SIZE(A1), %xmm0
	movsd	 -32 * SIZE(A1, LDA), %xmm1
	movsd	 -32 * SIZE(A2), %xmm2

	movsd	 -32 * SIZE(Y1), %xmm8

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8
	mulps	 %xmm13, %xmm1
	addps	 %xmm1,  %xmm8
	mulps	 %xmm14, %xmm2
	addps	 %xmm2, %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(11XX):
	movhps	-32 * SIZE(A1, LDA), %xmm8
	movhps	-32 * SIZE(A2, LDA), %xmm9

	movq	MM,  I
	sarq	$4,  I
	jle	L(115)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -28 * SIZE(A1), %xmm5
	movaps	 -24 * SIZE(A1), %xmm6

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YL1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YL1(-20 * SIZE, Y1, %xmm3)

	decq	 I
	jle	 L(114)
	ALIGN_3

L(113):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A1)
#endif

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A1), %xmm7
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -30 * SIZE(A1, LDA), %xmm4
	mulps	 %xmm12, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -26 * SIZE(A1, LDA), %xmm5
	mulps	 %xmm12, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	 -22 * SIZE(A1, LDA), %xmm6

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A1, LDA)
#endif

	shufps	 $0x4e, %xmm4, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0
	movaps	 -18 * SIZE(A1, LDA), %xmm8
	shufps	 $0x4e, %xmm5, %xmm4
	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -32 * SIZE(A2), %xmm4
	shufps	 $0x4e, %xmm6, %xmm5
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -28 * SIZE(A2), %xmm5
	shufps	 $0x4e, %xmm8, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	 -24 * SIZE(A2), %xmm6

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A2)
#endif

	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A2), %xmm7
	mulps	 %xmm14, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -16 * SIZE(A1), %xmm4
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -12 * SIZE(A1), %xmm5
	mulps	 %xmm14, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	  -8 * SIZE(A1), %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(Y1)
#endif

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	MOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-12 * SIZE, Y1, %xmm1)
	MOVUPS_YL1( -8 * SIZE, Y1, %xmm2)
	MOVUPS_YL1( -4 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(113)
	ALIGN_3

L(114):
	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A1), %xmm7
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -30 * SIZE(A1, LDA), %xmm4
	mulps	 %xmm12, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -26 * SIZE(A1, LDA), %xmm5
	mulps	 %xmm12, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	 -22 * SIZE(A1, LDA), %xmm6

	shufps	 $0x4e, %xmm4, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0
	movaps	 -18 * SIZE(A1, LDA), %xmm8
	shufps	 $0x4e, %xmm5, %xmm4
	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -32 * SIZE(A2), %xmm4
	shufps	 $0x4e, %xmm6, %xmm5
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -28 * SIZE(A2), %xmm5
	shufps	 $0x4e, %xmm8, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	 -24 * SIZE(A2), %xmm6

	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A2), %xmm7
	mulps	 %xmm14, %xmm5
	addps	 %xmm5,  %xmm1
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm2
	mulps	 %xmm14, %xmm7
	addps	 %xmm7,  %xmm3

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1
	ALIGN_3

L(115):
	testq	$8, MM
	je	L(116)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -28 * SIZE(A1), %xmm5
	movaps	 -30 * SIZE(A1, LDA), %xmm6
	movaps	 -26 * SIZE(A1, LDA), %xmm7

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -32 * SIZE(A2), %xmm4
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -28 * SIZE(A2), %xmm10

	shufps	 $0x4e, %xmm6, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0
	shufps	 $0x4e, %xmm7, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm1
	movaps	 %xmm7, %xmm8

	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm0
	mulps	 %xmm14, %xmm10
	addps	 %xmm10,  %xmm1

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)

	addq	 $8 * SIZE, A1
	addq	 $8 * SIZE, A2
	addq	 $8 * SIZE, Y1
	ALIGN_3

L(116):
	testq	$4, MM
	je	L(117)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -30 * SIZE(A1, LDA), %xmm5
	movaps	 -32 * SIZE(A2), %xmm6

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	shufps	 $0x4e, %xmm5, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm0

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, A2
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(117):
	testq	$2, MM
	je	L(118)

	movsd	 -32 * SIZE(A1), %xmm4
	movsd	 -32 * SIZE(A1, LDA), %xmm5
	movsd	 -32 * SIZE(A2), %xmm6

	movsd	 -32 * SIZE(Y1), %xmm0

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm0
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm0

	movlps	 %xmm0, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(118):
	testq	$1, MM
	je	L(990)

	movss	 -32 * SIZE(Y1), %xmm0

	movss	 -32 * SIZE(A1), %xmm4
	movss	 -32 * SIZE(A1, LDA), %xmm5
	movss	 -32 * SIZE(A2), %xmm6

	mulss	 %xmm12, %xmm4
	addss	 %xmm4,  %xmm0
	mulss	 %xmm13, %xmm5
	addss	 %xmm5,  %xmm0
	mulss	 %xmm14, %xmm6
	addss	 %xmm6,  %xmm0

	movss	 %xmm0, -32 * SIZE(Y1)
	jmp	 L(990)
	ALIGN_3

L(120):
	cmpq	$2, N
	jl	L(130)

	leaq	32 * SIZE(BUFFER), Y1
	movq	A,  A1
	leaq	(A,  LDA, 1), A2
	leaq	(A,  LDA, 2), A

	movss	(X), %xmm12
	addq	INCX, X
	movss	(X), %xmm13
	addq	INCX, X

	movss	ALPHA, %xmm0

	mulss	%xmm0, %xmm12
	mulss	%xmm0, %xmm13

	shufps	$0, %xmm12, %xmm12
	shufps	$0, %xmm13, %xmm13

	cmpq	$3, M
	jle	L(127)

	testq	$SIZE, A1
	je	L(12X)

	movss	 -32 * SIZE(Y1), %xmm8

	movss	 -32 * SIZE(A1), %xmm0
	movss	 -32 * SIZE(A2), %xmm1

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8
	mulss	 %xmm13, %xmm1
	addss	 %xmm1,  %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, Y1
	ALIGN_3

L(12X):
	testq	$2 * SIZE, A1
	je	L(12XX)

	movsd	 -32 * SIZE(A1), %xmm0
	movsd	 -32 * SIZE(A2), %xmm1

	movsd	 -32 * SIZE(Y1), %xmm8

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8
	mulps	 %xmm13, %xmm1
	addps	 %xmm1,  %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(12XX):
	movhps	-32 * SIZE(A2), %xmm8

	movq	MM,  I
	sarq	$4,  I
	jle	L(125)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -28 * SIZE(A1), %xmm5
	movaps	 -24 * SIZE(A1), %xmm6

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YL1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YL1(-20 * SIZE, Y1, %xmm3)

	decq	 I
	jle	 L(124)
	ALIGN_3

L(123):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A1)
#endif

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A1), %xmm7
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -30 * SIZE(A1, LDA), %xmm4
	mulps	 %xmm12, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -26 * SIZE(A1, LDA), %xmm5
	mulps	 %xmm12, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	 -22 * SIZE(A1, LDA), %xmm6

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A2)
#endif

	shufps	 $0x4e, %xmm4, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0
	movaps	 -18 * SIZE(A1, LDA), %xmm8
	shufps	 $0x4e, %xmm5, %xmm4
	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -16 * SIZE(A1), %xmm4
	shufps	 $0x4e, %xmm6, %xmm5
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -12 * SIZE(A1), %xmm5
	shufps	 $0x4e, %xmm8, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	  -8 * SIZE(A1), %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(Y1)
#endif

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	MOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-12 * SIZE, Y1, %xmm1)
	MOVUPS_YL1( -8 * SIZE, Y1, %xmm2)
	MOVUPS_YL1( -4 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(123)
	ALIGN_3

L(124):
	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A1), %xmm7
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -30 * SIZE(A1, LDA), %xmm4
	mulps	 %xmm12, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -26 * SIZE(A1, LDA), %xmm5
	mulps	 %xmm12, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	 -22 * SIZE(A1, LDA), %xmm6

	shufps	 $0x4e, %xmm4, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0
	movaps	 -18 * SIZE(A1, LDA), %xmm8
	shufps	 $0x4e, %xmm5, %xmm4
	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm1
	shufps	 $0x4e, %xmm6, %xmm5
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm2
	shufps	 $0x4e, %xmm8, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm3

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1
	ALIGN_3

L(125):
	testq	$8, MM
	je	L(126)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -28 * SIZE(A1), %xmm5
	movaps	 -30 * SIZE(A2), %xmm6
	movaps	 -26 * SIZE(A2), %xmm7

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1

	shufps	 $0x4e, %xmm6, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0
	shufps	 $0x4e, %xmm7, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm1
	movaps	 %xmm7, %xmm8

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)

	addq	 $8 * SIZE, A1
	addq	 $8 * SIZE, A2
	addq	 $8 * SIZE, Y1
	ALIGN_3

L(126):
	testq	$4, MM
	je	L(127)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -30 * SIZE(A2), %xmm5

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	shufps	 $0x4e, %xmm5, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, A2
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(127):
	testq	$2, MM
	je	L(128)

	movsd	 -32 * SIZE(A1), %xmm0
	movsd	 -32 * SIZE(A2), %xmm1

	movsd	 -32 * SIZE(Y1), %xmm8

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8
	mulps	 %xmm13, %xmm1
	addps	 %xmm1,  %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(128):
	testq	$1, MM
	je	L(990)

	movss	 -32 * SIZE(Y1), %xmm8

	movss	 -32 * SIZE(A1), %xmm0
	movss	 -32 * SIZE(A2), %xmm1

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8
	mulss	 %xmm13, %xmm1
	addss	 %xmm1,  %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)
	jmp	 L(990)
	ALIGN_3

L(130):
	cmpq	$1, N
	jne	L(990)

	leaq	32 * SIZE(BUFFER), Y1
	movq	A,  A1

	movss	(X), %xmm12

	mulss	ALPHA, %xmm12
	shufps	$0, %xmm12, %xmm12

	cmpq	$3, M
	jle	L(137)

	testq	$SIZE, A1
	je	L(13X)

	movss	 -32 * SIZE(Y1), %xmm8
	movss	 -32 * SIZE(A1), %xmm0

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, Y1
	ALIGN_3

L(13X):
	testq	$2 * SIZE, A1
	je	L(13XX)

	movsd	 -32 * SIZE(Y1), %xmm8
	movsd	 -32 * SIZE(A1), %xmm0

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(13XX):
	movq	MM,  I
	sarq	$4,  I
	jle	L(135)

	movaps	-32 * SIZE(A1), %xmm8
	movaps	-28 * SIZE(A1), %xmm9
	movaps	-24 * SIZE(A1), %xmm10
	movaps	-20 * SIZE(A1), %xmm11

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YL1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YL1(-20 * SIZE, Y1, %xmm3)

	decq	 I
	jle	 L(134)
	ALIGN_3

L(133):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(A1)
#endif

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	movaps	-16 * SIZE(A1), %xmm8
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	movaps	-12 * SIZE(A1), %xmm9
	mulps	 %xmm12, %xmm10
	addps	 %xmm10,  %xmm2
	movaps	 -8 * SIZE(A1), %xmm10
	mulps	 %xmm12, %xmm11
	addps	 %xmm11,  %xmm3
	movaps	 -4 * SIZE(A1), %xmm11

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(Y1)
#endif

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	MOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-12 * SIZE, Y1, %xmm1)
	MOVUPS_YL1( -8 * SIZE, Y1, %xmm2)
	MOVUPS_YL1( -4 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(133)
	ALIGN_3

L(134):
	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	mulps	 %xmm12, %xmm10
	addps	 %xmm10,  %xmm2
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	mulps	 %xmm12, %xmm11
	addps	 %xmm11,  %xmm3
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, Y1
	ALIGN_3

L(135):
	testq	$8, MM
	je	L(136)

	movaps	-32 * SIZE(A1), %xmm8
	movaps	-28 * SIZE(A1), %xmm9

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)

	addq	 $8 * SIZE, A1
	addq	 $8 * SIZE, Y1
	ALIGN_3

L(136):
	testq	$4, MM
	je	L(137)

	movaps	-32 * SIZE(A1), %xmm8
	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(137):
	testq	$2, MM
	je	L(138)

	movsd	 -32 * SIZE(Y1), %xmm8
	movsd	 -32 * SIZE(A1), %xmm0

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(138):
	testq	$1, MM
	je	L(990)

	movss	 -32 * SIZE(Y1), %xmm8
	movss	 -32 * SIZE(A1), %xmm0

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)
	jmp	L(990)
	ALIGN_3

L(200):
	testq	$2 * SIZE, LDA
	jne	L(300)

	cmpq	$4, N
	jl	L(210)
	ALIGN_3

L(201):
	subq	$4, N

	leaq	32 * SIZE(BUFFER), Y1
	movq	A,  A1
	leaq	(A,  LDA, 2), A2
	leaq	(A,  LDA, 4), A

	movss	(X), %xmm12
	addq	INCX, X
	movss	(X), %xmm13
	addq	INCX, X
	movss	(X), %xmm14
	addq	INCX, X
	movss	(X), %xmm15
	addq	INCX, X

	movss	ALPHA, %xmm0

	mulss	%xmm0, %xmm12
	mulss	%xmm0, %xmm13
	mulss	%xmm0, %xmm14
	mulss	%xmm0, %xmm15

	shufps	$0, %xmm12, %xmm12
	shufps	$0, %xmm13, %xmm13
	shufps	$0, %xmm14, %xmm14
	shufps	$0, %xmm15, %xmm15

	cmpq	$3, M
	jle	L(207)

	testq	$SIZE, A1
	je	L(20X)

	movss	 -32 * SIZE(A1), %xmm0
	movss	 -32 * SIZE(A1, LDA), %xmm1
	movss	 -32 * SIZE(A2), %xmm2
	movss	 -32 * SIZE(A2, LDA), %xmm3

	movss	 -32 * SIZE(Y1), %xmm8

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8
	mulss	 %xmm13, %xmm1
	addss	 %xmm1,  %xmm8
	mulss	 %xmm14, %xmm2
	addss	 %xmm2, %xmm8
	mulss	 %xmm15, %xmm3
	addss	 %xmm3, %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, Y1
	ALIGN_3

L(20X):
	testq	$2 * SIZE, A1
	je	L(20XX)

	movsd	 -32 * SIZE(A1), %xmm0
	movsd	 -32 * SIZE(A1, LDA), %xmm1
	movsd	 -32 * SIZE(A2), %xmm2
	movsd	 -32 * SIZE(A2, LDA), %xmm3

	movsd	 -32 * SIZE(Y1), %xmm8

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8
	mulps	 %xmm13, %xmm1
	addps	 %xmm1,  %xmm8
	mulps	 %xmm14, %xmm2
	addps	 %xmm2, %xmm8
	mulps	 %xmm15, %xmm3
	addps	 %xmm3, %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(20XX):
	movaps	-33 * SIZE(A1, LDA), %xmm8
	movaps	-34 * SIZE(A2), %xmm9
	movaps	-35 * SIZE(A2, LDA), %xmm10

	movq	MM,  I
	sarq	$4,  I
	jle	L(205)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -28 * SIZE(A1), %xmm5
	movaps	 -24 * SIZE(A1), %xmm6

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YL1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YL1(-20 * SIZE, Y1, %xmm3)

	decq	 I
	jle	 L(204)
	ALIGN_3

L(203):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1)
#endif

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A1), %xmm7
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -29 * SIZE(A1, LDA), %xmm4
	mulps	 %xmm12, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -25 * SIZE(A1, LDA), %xmm5
	mulps	 %xmm12, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	 -21 * SIZE(A1, LDA), %xmm6

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1, LDA)
#endif

	movss	 %xmm4, %xmm8
	shufps	 $0x39, %xmm8, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0
	movaps	 -17 * SIZE(A1, LDA), %xmm8
	movss	 %xmm5, %xmm4
	shufps	 $0x39, %xmm4, %xmm4
	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -30 * SIZE(A2), %xmm4
	movss	 %xmm6, %xmm5
	shufps	 $0x39, %xmm5, %xmm5
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -26 * SIZE(A2), %xmm5
	movss	 %xmm8, %xmm6
	shufps	 $0x39, %xmm6, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	 -22 * SIZE(A2), %xmm6

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2)
#endif

	shufps	 $0x4e, %xmm4, %xmm9
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm0
	movaps	 -18 * SIZE(A2), %xmm9
	shufps	 $0x4e, %xmm5, %xmm4
	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -31 * SIZE(A2, LDA), %xmm4
	shufps	 $0x4e, %xmm6, %xmm5
	mulps	 %xmm14, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -27 * SIZE(A2, LDA), %xmm5
	shufps	 $0x4e, %xmm9, %xmm6
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	 -23 * SIZE(A2, LDA), %xmm6

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2, LDA)
#endif

	movss	 %xmm4,  %xmm10
	shufps	 $0x93, %xmm4, %xmm10
	mulps	 %xmm15, %xmm10
	addps	 %xmm10, %xmm0
	movaps	 -19 * SIZE(A2, LDA), %xmm10
	movss	 %xmm5,  %xmm4
	shufps	 $0x93, %xmm5, %xmm4
	mulps	 %xmm15, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -16 * SIZE(A1), %xmm4
	movss	 %xmm6,  %xmm5
	shufps	 $0x93, %xmm6, %xmm5
	mulps	 %xmm15, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -12 * SIZE(A1), %xmm5
	movss	 %xmm10,  %xmm6
	shufps	 $0x93, %xmm10, %xmm6
	mulps	 %xmm15, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	  -8 * SIZE(A1), %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(Y1)
#endif

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	MOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-12 * SIZE, Y1, %xmm1)
	MOVUPS_YL1( -8 * SIZE, Y1, %xmm2)
	MOVUPS_YL1( -4 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(203)
	ALIGN_3

L(204):
	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A1), %xmm7
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -29 * SIZE(A1, LDA), %xmm4
	mulps	 %xmm12, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -25 * SIZE(A1, LDA), %xmm5
	mulps	 %xmm12, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	 -21 * SIZE(A1, LDA), %xmm6

	movss	 %xmm4, %xmm8
	shufps	 $0x39, %xmm8, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0
	movaps	 -17 * SIZE(A1, LDA), %xmm8
	movss	 %xmm5, %xmm4
	shufps	 $0x39, %xmm4, %xmm4
	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -30 * SIZE(A2), %xmm4

	movss	 %xmm6, %xmm5
	shufps	 $0x39, %xmm5, %xmm5
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -26 * SIZE(A2), %xmm5
	movss	 %xmm8, %xmm6
	shufps	 $0x39, %xmm6, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	 -22 * SIZE(A2), %xmm6

	shufps	 $0x4e, %xmm4, %xmm9
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm0
	movaps	 -18 * SIZE(A2), %xmm9
	shufps	 $0x4e, %xmm5, %xmm4
	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -31 * SIZE(A2, LDA), %xmm4
	shufps	 $0x4e, %xmm6, %xmm5
	mulps	 %xmm14, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -27 * SIZE(A2, LDA), %xmm5
	shufps	 $0x4e, %xmm9, %xmm6
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	 -23 * SIZE(A2, LDA), %xmm6

	movss	 %xmm4,  %xmm10
	shufps	 $0x93, %xmm4, %xmm10
	mulps	 %xmm15, %xmm10
	addps	 %xmm10, %xmm0
	movaps	 -19 * SIZE(A2, LDA), %xmm10
	movss	 %xmm5,  %xmm4
	shufps	 $0x93, %xmm5, %xmm4
	mulps	 %xmm15, %xmm4
	addps	 %xmm4,  %xmm1

	movss	 %xmm6,  %xmm5
	shufps	 $0x93, %xmm6, %xmm5
	mulps	 %xmm15, %xmm5
	addps	 %xmm5,  %xmm2
	movss	 %xmm10,  %xmm6
	shufps	 $0x93, %xmm10, %xmm6
	mulps	 %xmm15, %xmm6
	addps	 %xmm6,  %xmm3

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1
	ALIGN_3

L(205):
	testq	$8, MM
	je	L(206)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -28 * SIZE(A1), %xmm5
	movaps	 -29 * SIZE(A1, LDA), %xmm6

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -30 * SIZE(A2), %xmm4
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -26 * SIZE(A2), %xmm5

	movss	 %xmm6, %xmm8
	shufps	 $0x39, %xmm8, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0
	movaps	 -25 * SIZE(A1, LDA), %xmm8
	movss	 %xmm8, %xmm6
	shufps	 $0x39, %xmm6, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm1
	movaps	 -31 * SIZE(A2, LDA), %xmm6

	shufps	 $0x4e, %xmm4, %xmm9
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm0
	movaps	 -27 * SIZE(A2, LDA), %xmm7

	shufps	 $0x4e, %xmm5, %xmm4
	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 %xmm5, %xmm9

	movss	 %xmm6,  %xmm10
	shufps	 $0x93, %xmm6, %xmm10
	mulps	 %xmm15, %xmm10
	addps	 %xmm10, %xmm0

	movss	 %xmm7,  %xmm6
	shufps	 $0x93, %xmm7, %xmm6
	mulps	 %xmm15, %xmm6
	addps	 %xmm6,  %xmm1
	movaps	 %xmm7, %xmm10

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)

	addq	 $8 * SIZE, A1
	addq	 $8 * SIZE, A2
	addq	 $8 * SIZE, Y1
	ALIGN_3

L(206):
	testq	$4, MM
	je	L(207)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -29 * SIZE(A1, LDA), %xmm5
	movaps	 -30 * SIZE(A2), %xmm6
	movaps	 -31 * SIZE(A2, LDA), %xmm7

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movss	 %xmm5, %xmm8
	shufps	 $0x39, %xmm8, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0

	shufps	 $0x4e, %xmm6, %xmm9
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm0
	movss	 %xmm7,  %xmm10
	shufps	 $0x93, %xmm7, %xmm10
	mulps	 %xmm15, %xmm10
	addps	 %xmm10, %xmm0

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, A2
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(207):
	testq	$2, MM
	je	L(208)

	movsd	 -32 * SIZE(A1), %xmm4
	movsd	 -32 * SIZE(A1, LDA), %xmm5
	movsd	 -32 * SIZE(A2), %xmm6
	movsd	 -32 * SIZE(A2, LDA), %xmm7

	movsd	 -32 * SIZE(Y1), %xmm0

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm0
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm0
	mulps	 %xmm15, %xmm7
	addps	 %xmm7,  %xmm0

	movlps	 %xmm0, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(208):
	testq	$1, MM
	je	L(209)

	movss	 -32 * SIZE(Y1), %xmm0

	movss	 -32 * SIZE(A1), %xmm4
	movss	 -32 * SIZE(A1, LDA), %xmm5
	movss	 -32 * SIZE(A2), %xmm6
	movss	 -32 * SIZE(A2, LDA), %xmm7

	mulss	 %xmm12, %xmm4
	addss	 %xmm4,  %xmm0
	mulss	 %xmm13, %xmm5
	addss	 %xmm5,  %xmm0
	mulss	 %xmm14, %xmm6
	addss	 %xmm6,  %xmm0
	mulss	 %xmm15, %xmm7
	addss	 %xmm7,  %xmm0

	movss	 %xmm0, -32 * SIZE(Y1)
	ALIGN_3

L(209):
	cmpq	$4, N
	jge	L(201)
	ALIGN_3

L(210):
	cmpq	$3, N
	jne	L(220)

	leaq	32 * SIZE(BUFFER), Y1
	movq	A,  A1
	leaq	(A,  LDA, 2), A2
	leaq	(A,  LDA, 4), A

	movss	(X), %xmm12
	addq	INCX, X
	movss	(X), %xmm13
	addq	INCX, X
	movss	(X), %xmm14
	addq	INCX, X

	movss	ALPHA, %xmm0

	mulss	%xmm0, %xmm12
	mulss	%xmm0, %xmm13
	mulss	%xmm0, %xmm14

	shufps	$0, %xmm12, %xmm12
	shufps	$0, %xmm13, %xmm13
	shufps	$0, %xmm14, %xmm14

	cmpq	$3, M
	jle	L(217)

	testq	$SIZE, A1
	je	L(21X)

	movss	 -32 * SIZE(A1), %xmm0
	movss	 -32 * SIZE(A1, LDA), %xmm1
	movss	 -32 * SIZE(A2), %xmm2

	movss	 -32 * SIZE(Y1), %xmm8

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8
	mulss	 %xmm13, %xmm1
	addss	 %xmm1,  %xmm8
	mulss	 %xmm14, %xmm2
	addss	 %xmm2, %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, Y1
	ALIGN_3

L(21X):
	testq	$2 * SIZE, A1
	je	L(21XX)

	movsd	 -32 * SIZE(A1), %xmm0
	movsd	 -32 * SIZE(A1, LDA), %xmm1
	movsd	 -32 * SIZE(A2), %xmm2

	movsd	 -32 * SIZE(Y1), %xmm8

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8
	mulps	 %xmm13, %xmm1
	addps	 %xmm1,  %xmm8
	mulps	 %xmm14, %xmm2
	addps	 %xmm2, %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(21XX):
	movaps	-33 * SIZE(A1, LDA), %xmm8
	movaps	-34 * SIZE(A2), %xmm9
	movaps	-35 * SIZE(A2, LDA), %xmm10

	movq	MM,  I
	sarq	$4,  I
	jle	L(215)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -28 * SIZE(A1), %xmm5
	movaps	 -24 * SIZE(A1), %xmm6

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YL1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YL1(-20 * SIZE, Y1, %xmm3)

	decq	 I
	jle	 L(214)
	ALIGN_3

L(213):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A1)
#endif

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A1), %xmm7
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -29 * SIZE(A1, LDA), %xmm4
	mulps	 %xmm12, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -25 * SIZE(A1, LDA), %xmm5
	mulps	 %xmm12, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	 -21 * SIZE(A1, LDA), %xmm6

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A1, LDA)
#endif

	movss	 %xmm4, %xmm8
	shufps	 $0x39, %xmm8, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0
	movaps	 -17 * SIZE(A1, LDA), %xmm8
	movss	 %xmm5, %xmm4
	shufps	 $0x39, %xmm4, %xmm4
	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -30 * SIZE(A2), %xmm4
	movss	 %xmm6, %xmm5
	shufps	 $0x39, %xmm5, %xmm5
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -26 * SIZE(A2), %xmm5
	movss	 %xmm8, %xmm6
	shufps	 $0x39, %xmm6, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	 -22 * SIZE(A2), %xmm6

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A2)
#endif

	shufps	 $0x4e, %xmm4, %xmm9
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm0
	movaps	 -18 * SIZE(A2), %xmm9
	shufps	 $0x4e, %xmm5, %xmm4
	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -16 * SIZE(A1), %xmm4
	shufps	 $0x4e, %xmm6, %xmm5
	mulps	 %xmm14, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -12 * SIZE(A1), %xmm5
	shufps	 $0x4e, %xmm9, %xmm6
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	  -8 * SIZE(A1), %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(Y1)
#endif

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	MOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-12 * SIZE, Y1, %xmm1)
	MOVUPS_YL1( -8 * SIZE, Y1, %xmm2)
	MOVUPS_YL1( -4 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(213)
	ALIGN_3

L(214):
	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A1), %xmm7
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -29 * SIZE(A1, LDA), %xmm4
	mulps	 %xmm12, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -25 * SIZE(A1, LDA), %xmm5
	mulps	 %xmm12, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	 -21 * SIZE(A1, LDA), %xmm6

	movss	 %xmm4, %xmm8
	shufps	 $0x39, %xmm8, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0
	movaps	 -17 * SIZE(A1, LDA), %xmm8
	movss	 %xmm5, %xmm4
	shufps	 $0x39, %xmm4, %xmm4
	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -30 * SIZE(A2), %xmm4

	movss	 %xmm6, %xmm5
	shufps	 $0x39, %xmm5, %xmm5
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -26 * SIZE(A2), %xmm5
	movss	 %xmm8, %xmm6
	shufps	 $0x39, %xmm6, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	 -22 * SIZE(A2), %xmm6

	shufps	 $0x4e, %xmm4, %xmm9
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm0
	movaps	 -18 * SIZE(A2), %xmm9
	shufps	 $0x4e, %xmm5, %xmm4
	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm1
	shufps	 $0x4e, %xmm6, %xmm5
	mulps	 %xmm14, %xmm5
	addps	 %xmm5,  %xmm2
	shufps	 $0x4e, %xmm9, %xmm6
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm3

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1
	ALIGN_3

L(215):
	testq	$8, MM
	je	L(216)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -28 * SIZE(A1), %xmm5
	movaps	 -29 * SIZE(A1, LDA), %xmm6

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -30 * SIZE(A2), %xmm4
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -26 * SIZE(A2), %xmm5

	movss	 %xmm6, %xmm8
	shufps	 $0x39, %xmm8, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0
	movaps	 -25 * SIZE(A1, LDA), %xmm8
	movss	 %xmm8, %xmm6
	shufps	 $0x39, %xmm6, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm1

	shufps	 $0x4e, %xmm4, %xmm9
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm0

	shufps	 $0x4e, %xmm5, %xmm4
	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 %xmm5, %xmm9

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)

	addq	 $8 * SIZE, A1
	addq	 $8 * SIZE, A2
	addq	 $8 * SIZE, Y1
	ALIGN_3

L(216):
	testq	$4, MM
	je	L(217)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -29 * SIZE(A1, LDA), %xmm5
	movaps	 -30 * SIZE(A2), %xmm6

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movss	 %xmm5, %xmm8
	shufps	 $0x39, %xmm8, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0

	shufps	 $0x4e, %xmm6, %xmm9
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm0
	movss	 %xmm7,  %xmm10

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, A2
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(217):
	testq	$2, MM
	je	L(218)

	movsd	 -32 * SIZE(A1), %xmm4
	movsd	 -32 * SIZE(A1, LDA), %xmm5
	movsd	 -32 * SIZE(A2), %xmm6

	movsd	 -32 * SIZE(Y1), %xmm0

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm0
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm0

	movlps	 %xmm0, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(218):
	testq	$1, MM
	je	L(990)

	movss	 -32 * SIZE(Y1), %xmm0

	movss	 -32 * SIZE(A1), %xmm4
	movss	 -32 * SIZE(A1, LDA), %xmm5
	movss	 -32 * SIZE(A2), %xmm6

	mulss	 %xmm12, %xmm4
	addss	 %xmm4,  %xmm0
	mulss	 %xmm13, %xmm5
	addss	 %xmm5,  %xmm0
	mulss	 %xmm14, %xmm6
	addss	 %xmm6,  %xmm0

	movss	 %xmm0, -32 * SIZE(Y1)
	jmp	 L(990)
	ALIGN_4

L(220):
	testq	N, N
	jle	L(990)

	cmpq	$2, N
	jne	L(230)

	leaq	32 * SIZE(BUFFER), Y1
	movq	A,  A1
	leaq	(A,  LDA, 1), A2
	leaq	(A,  LDA, 2), A

	movss	(X), %xmm12
	addq	INCX, X
	movss	(X), %xmm13
	addq	INCX, X

	movss	ALPHA, %xmm0

	mulss	%xmm0, %xmm12
	mulss	%xmm0, %xmm13

	shufps	$0, %xmm12, %xmm12
	shufps	$0, %xmm13, %xmm13

	cmpq	$3, M
	jle	L(227)

	testq	$SIZE, A1
	je	L(22X)

	movss	 -32 * SIZE(Y1), %xmm9

	movss	 -32 * SIZE(A1), %xmm0
	movss	 -32 * SIZE(A2), %xmm1

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm9
	mulss	 %xmm13, %xmm1
	addss	 %xmm1,  %xmm9

	movss	 %xmm9, -32 * SIZE(Y1)

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, Y1
	ALIGN_3

L(22X):
	testq	$2 * SIZE, A1
	je	L(22XX)

	movsd	 -32 * SIZE(A1), %xmm0
	movsd	 -32 * SIZE(A2), %xmm1

	movsd	 -32 * SIZE(Y1), %xmm9

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm9
	mulps	 %xmm13, %xmm1
	addps	 %xmm1,  %xmm9

	movlps	 %xmm9, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(22XX):
	movaps	-33 * SIZE(A1, LDA), %xmm8

	movq	MM,  I
	sarq	$4,  I
	jle	L(225)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -28 * SIZE(A1), %xmm5
	movaps	 -24 * SIZE(A1), %xmm6

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YL1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YL1(-20 * SIZE, Y1, %xmm3)

	decq	 I
	jle	 L(224)
	ALIGN_3

L(223):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A1)
#endif

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A1), %xmm7
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -29 * SIZE(A2), %xmm4
	mulps	 %xmm12, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -25 * SIZE(A2), %xmm5
	mulps	 %xmm12, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	 -21 * SIZE(A2), %xmm6

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A2)
#endif

	movss	 %xmm4, %xmm8
	shufps	 $0x39, %xmm8, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0
	movaps	 -17 * SIZE(A2), %xmm8
	movss	 %xmm5, %xmm4
	shufps	 $0x39, %xmm4, %xmm4
	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -16 * SIZE(A1), %xmm4

	movss	 %xmm6, %xmm5
	shufps	 $0x39, %xmm5, %xmm5
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -12 * SIZE(A1), %xmm5
	movss	 %xmm8, %xmm6
	shufps	 $0x39, %xmm6, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	  -8 * SIZE(A1), %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(Y1)
#endif

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	MOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-12 * SIZE, Y1, %xmm1)
	MOVUPS_YL1( -8 * SIZE, Y1, %xmm2)
	MOVUPS_YL1( -4 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(223)
	ALIGN_3

L(224):
	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A1), %xmm7
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -29 * SIZE(A2), %xmm4
	mulps	 %xmm12, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -25 * SIZE(A2), %xmm5
	mulps	 %xmm12, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	 -21 * SIZE(A2), %xmm6

	movss	 %xmm4, %xmm8
	shufps	 $0x39, %xmm8, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0
	movaps	 -17 * SIZE(A2), %xmm8
	movss	 %xmm5, %xmm4
	shufps	 $0x39, %xmm4, %xmm4
	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm1

	movss	 %xmm6, %xmm5
	shufps	 $0x39, %xmm5, %xmm5
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm2
	movss	 %xmm8, %xmm6
	shufps	 $0x39, %xmm6, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm3

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1
	ALIGN_3

L(225):
	testq	$8, MM
	je	L(226)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -28 * SIZE(A1), %xmm5

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -29 * SIZE(A2), %xmm6
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -25 * SIZE(A2), %xmm7

	movss	 %xmm6, %xmm8
	shufps	 $0x39, %xmm8, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0

	movss	 %xmm7, %xmm6
	shufps	 $0x39, %xmm6, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm1
	movaps	 %xmm7, %xmm8

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)

	addq	 $8 * SIZE, A1
	addq	 $8 * SIZE, A2
	addq	 $8 * SIZE, Y1
	ALIGN_3

L(226):
	testq	$4, MM
	je	L(227)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -29 * SIZE(A2), %xmm5

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movss	 %xmm5, %xmm8
	shufps	 $0x39, %xmm8, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8,  %xmm0

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, A2
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(227):
	testq	$2, MM
	je	L(228)

	movsd	 -32 * SIZE(A1), %xmm0
	movsd	 -32 * SIZE(A2), %xmm1

	movsd	 -32 * SIZE(Y1), %xmm9

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm9
	mulps	 %xmm13, %xmm1
	addps	 %xmm1,  %xmm9

	movlps	 %xmm9, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(228):
	testq	$1, MM
	je	L(990)

	movss	 -32 * SIZE(Y1), %xmm9

	movss	 -32 * SIZE(A1), %xmm0
	movss	 -32 * SIZE(A2), %xmm1

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm9
	mulss	 %xmm13, %xmm1
	addss	 %xmm1,  %xmm9

	movss	 %xmm9, -32 * SIZE(Y1)
	jmp	 L(990)
	ALIGN_3

L(230):
	cmpq	$1, N
	jne	L(990)

	leaq	32 * SIZE(BUFFER), Y1
	movq	A,  A1

	movss	(X), %xmm12

	mulss	ALPHA, %xmm12
	shufps	$0, %xmm12, %xmm12

	cmpq	$3, M
	jle	L(237)

	testq	$SIZE, A1
	je	L(23X)

	movss	 -32 * SIZE(Y1), %xmm8
	movss	 -32 * SIZE(A1), %xmm0

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, Y1
	ALIGN_3

L(23X):
	testq	$2 * SIZE, A1
	je	L(23XX)

	movsd	 -32 * SIZE(Y1), %xmm8
	movsd	 -32 * SIZE(A1), %xmm0

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(23XX):
	testq	$2 * SIZE, A1
	jne	L(230)

	movq	MM,  I
	sarq	$4,  I
	jle	L(235)

	movaps	-32 * SIZE(A1), %xmm8
	movaps	-28 * SIZE(A1), %xmm9
	movaps	-24 * SIZE(A1), %xmm10
	movaps	-20 * SIZE(A1), %xmm11

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YL1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YL1(-20 * SIZE, Y1, %xmm3)

	decq	 I
	jle	 L(234)
	ALIGN_3

L(233):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(A1)
#endif

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	movaps	-16 * SIZE(A1), %xmm8
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	movaps	-12 * SIZE(A1), %xmm9
	mulps	 %xmm12, %xmm10
	addps	 %xmm10,  %xmm2
	movaps	 -8 * SIZE(A1), %xmm10
	mulps	 %xmm12, %xmm11
	addps	 %xmm11,  %xmm3
	movaps	 -4 * SIZE(A1), %xmm11

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(Y1)
#endif

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	MOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-12 * SIZE, Y1, %xmm1)
	MOVUPS_YL1( -8 * SIZE, Y1, %xmm2)
	MOVUPS_YL1( -4 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(233)
	ALIGN_3

L(234):
	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	mulps	 %xmm12, %xmm10
	addps	 %xmm10,  %xmm2
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	mulps	 %xmm12, %xmm11
	addps	 %xmm11,  %xmm3
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, Y1
	ALIGN_3

L(235):
	testq	$8, MM
	je	L(236)

	movaps	-32 * SIZE(A1), %xmm8
	movaps	-28 * SIZE(A1), %xmm9

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)

	addq	 $8 * SIZE, A1
	addq	 $8 * SIZE, Y1
	ALIGN_3

L(236):
	testq	$4, MM
	je	L(237)

	movaps	-32 * SIZE(A1), %xmm8
	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(237):
	testq	$2, MM
	je	L(238)

	movsd	 -32 * SIZE(Y1), %xmm8
	movsd	 -32 * SIZE(A1), %xmm0

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(238):
	testq	$1, MM
	je	L(990)

	movss	 -32 * SIZE(Y1), %xmm8
	movss	 -32 * SIZE(A1), %xmm0

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)
	jmp	L(990)
	ALIGN_4

L(300):
	cmpq	$4, N
	jl	L(310)
	ALIGN_3

L(301):
	subq	$4, N

	leaq	32 * SIZE(BUFFER), Y1
	movq	A,  A1
	leaq	(A,  LDA, 2), A2
	leaq	(A,  LDA, 4), A

	movss	(X), %xmm12
	addq	INCX, X
	movss	(X), %xmm13
	addq	INCX, X
	movss	(X), %xmm14
	addq	INCX, X
	movss	(X), %xmm15
	addq	INCX, X

	movss	ALPHA, %xmm0

	mulss	%xmm0, %xmm12
	mulss	%xmm0, %xmm13
	mulss	%xmm0, %xmm14
	mulss	%xmm0, %xmm15

	shufps	$0, %xmm12, %xmm12
	shufps	$0, %xmm13, %xmm13
	shufps	$0, %xmm14, %xmm14
	shufps	$0, %xmm15, %xmm15

	cmpq	$3, M
	jle	L(307)

	testq	$SIZE, A1
	je	L(30X)

	movss	 -32 * SIZE(A1), %xmm0
	movss	 -32 * SIZE(A1, LDA), %xmm1
	movss	 -32 * SIZE(A2), %xmm2
	movss	 -32 * SIZE(A2, LDA), %xmm3

	movss	 -32 * SIZE(Y1), %xmm8

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8
	mulss	 %xmm13, %xmm1
	addss	 %xmm1,  %xmm8
	mulss	 %xmm14, %xmm2
	addss	 %xmm2, %xmm8
	mulss	 %xmm15, %xmm3
	addss	 %xmm3, %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, Y1
	ALIGN_3

L(30X):
	testq	$2 * SIZE, A1
	je	L(30XX)

	movsd	 -32 * SIZE(A1), %xmm0
	movsd	 -32 * SIZE(A1, LDA), %xmm1
	movsd	 -32 * SIZE(A2), %xmm2
	movsd	 -32 * SIZE(A2, LDA), %xmm3

	movsd	 -32 * SIZE(Y1), %xmm8

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8
	mulps	 %xmm13, %xmm1
	addps	 %xmm1,  %xmm8
	mulps	 %xmm14, %xmm2
	addps	 %xmm2, %xmm8
	mulps	 %xmm15, %xmm3
	addps	 %xmm3, %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(30XX):
	movaps	-35 * SIZE(A1, LDA), %xmm8
	movaps	-34 * SIZE(A2), %xmm9
	movaps	-33 * SIZE(A2, LDA), %xmm10

	movq	MM,  I
	sarq	$4,  I
	jle	L(305)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -28 * SIZE(A1), %xmm5
	movaps	 -24 * SIZE(A1), %xmm6

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YL1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YL1(-20 * SIZE, Y1, %xmm3)

	decq	 I
	jle	 L(304)
	ALIGN_3

L(303):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1)
#endif

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A1), %xmm7
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -31 * SIZE(A1, LDA), %xmm4
	mulps	 %xmm12, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -27 * SIZE(A1, LDA), %xmm5
	mulps	 %xmm12, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	 -23 * SIZE(A1, LDA), %xmm6

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1, LDA)
#endif

	movss	 %xmm4,  %xmm8
	shufps	 $0x93, %xmm4, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8, %xmm0
	movaps	 -19 * SIZE(A1, LDA), %xmm8
	movss	 %xmm5,  %xmm4
	shufps	 $0x93, %xmm5, %xmm4
	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -30 * SIZE(A2), %xmm4
	movss	 %xmm6,  %xmm5
	shufps	 $0x93, %xmm6, %xmm5
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -26 * SIZE(A2), %xmm5
	movss	 %xmm8,  %xmm6
	shufps	 $0x93, %xmm8, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	 -22 * SIZE(A2), %xmm6

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2)
#endif

	shufps	 $0x4e, %xmm4, %xmm9
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm0
	movaps	 -18 * SIZE(A2), %xmm9
	shufps	 $0x4e, %xmm5, %xmm4
	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -29 * SIZE(A2, LDA), %xmm4
	shufps	 $0x4e, %xmm6, %xmm5
	mulps	 %xmm14, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -25 * SIZE(A2, LDA), %xmm5
	shufps	 $0x4e, %xmm9, %xmm6
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	 -21 * SIZE(A2, LDA), %xmm6

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2, LDA)
#endif

	movss	 %xmm4, %xmm10
	shufps	 $0x39, %xmm10, %xmm10
	mulps	 %xmm15, %xmm10
	addps	 %xmm10,  %xmm0
	movaps	 -17 * SIZE(A2, LDA), %xmm10
	movss	 %xmm5, %xmm4
	shufps	 $0x39, %xmm4, %xmm4
	mulps	 %xmm15, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -16 * SIZE(A1), %xmm4
	movss	 %xmm6, %xmm5
	shufps	 $0x39, %xmm5, %xmm5
	mulps	 %xmm15, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -12 * SIZE(A1), %xmm5
	movss	 %xmm10, %xmm6
	shufps	 $0x39, %xmm6, %xmm6
	mulps	 %xmm15, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	  -8 * SIZE(A1), %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(Y1)
#endif

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	MOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-12 * SIZE, Y1, %xmm1)
	MOVUPS_YL1( -8 * SIZE, Y1, %xmm2)
	MOVUPS_YL1( -4 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(303)
	ALIGN_3

L(304):
	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A1), %xmm7
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -31 * SIZE(A1, LDA), %xmm4

	mulps	 %xmm12, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -27 * SIZE(A1, LDA), %xmm5
	mulps	 %xmm12, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	 -23 * SIZE(A1, LDA), %xmm6

	movss	 %xmm4,  %xmm8
	shufps	 $0x93, %xmm4, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8, %xmm0
	movaps	 -19 * SIZE(A1, LDA), %xmm8
	movss	 %xmm5,  %xmm4
	shufps	 $0x93, %xmm5, %xmm4
	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -30 * SIZE(A2), %xmm4

	movss	 %xmm6,  %xmm5
	shufps	 $0x93, %xmm6, %xmm5
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -26 * SIZE(A2), %xmm5
	movss	 %xmm8,  %xmm6
	shufps	 $0x93, %xmm8, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	 -22 * SIZE(A2), %xmm6

	shufps	 $0x4e, %xmm4, %xmm9
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm0
	movaps	 -18 * SIZE(A2), %xmm9
	shufps	 $0x4e, %xmm5, %xmm4
	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -29 * SIZE(A2, LDA), %xmm4
	shufps	 $0x4e, %xmm6, %xmm5
	mulps	 %xmm14, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -25 * SIZE(A2, LDA), %xmm5
	shufps	 $0x4e, %xmm9, %xmm6
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	 -21 * SIZE(A2, LDA), %xmm6

	movss	 %xmm4, %xmm10
	shufps	 $0x39, %xmm10, %xmm10
	mulps	 %xmm15, %xmm10
	addps	 %xmm10,  %xmm0
	movaps	 -17 * SIZE(A2, LDA), %xmm10
	movss	 %xmm5, %xmm4
	shufps	 $0x39, %xmm4, %xmm4
	mulps	 %xmm15, %xmm4
	addps	 %xmm4,  %xmm1
	movss	 %xmm6, %xmm5
	shufps	 $0x39, %xmm5, %xmm5
	mulps	 %xmm15, %xmm5
	addps	 %xmm5,  %xmm2
	movss	 %xmm10, %xmm6
	shufps	 $0x39, %xmm6, %xmm6
	mulps	 %xmm15, %xmm6
	addps	 %xmm6,  %xmm3

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1
	ALIGN_3

L(305):
	testq	$8, MM
	je	L(306)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -28 * SIZE(A1), %xmm5

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -31 * SIZE(A1, LDA), %xmm6
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -27 * SIZE(A1, LDA), %xmm7

	movss	 %xmm6,  %xmm8
	shufps	 $0x93, %xmm6, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8, %xmm0
	movaps	 -30 * SIZE(A2), %xmm4
	movss	 %xmm7,  %xmm6
	shufps	 $0x93, %xmm7, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm1
	movaps	 %xmm7, %xmm8
	movaps	 -26 * SIZE(A2), %xmm5

	shufps	 $0x4e, %xmm4, %xmm9
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm0
	movaps	 -29 * SIZE(A2, LDA), %xmm6
	shufps	 $0x4e, %xmm5, %xmm4
	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 %xmm5, %xmm9
	movaps	 -25 * SIZE(A2, LDA), %xmm7

	movss	 %xmm6, %xmm10
	shufps	 $0x39, %xmm10, %xmm10
	mulps	 %xmm15, %xmm10
	addps	 %xmm10,  %xmm0
	movss	 %xmm7, %xmm6
	shufps	 $0x39, %xmm6, %xmm6
	mulps	 %xmm15, %xmm6
	addps	 %xmm6,  %xmm1
	movaps	 %xmm7, %xmm10

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)

	addq	 $8 * SIZE, A1
	addq	 $8 * SIZE, A2
	addq	 $8 * SIZE, Y1
	ALIGN_3

L(306):
	testq	$4, MM
	je	L(307)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -31 * SIZE(A1, LDA), %xmm5

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -30 * SIZE(A2), %xmm6
	movss	 %xmm5,  %xmm8
	shufps	 $0x93, %xmm5, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8, %xmm0
	movaps	 -29 * SIZE(A2, LDA), %xmm7

	shufps	 $0x4e, %xmm6, %xmm9
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm0
	movss	 %xmm7, %xmm10
	shufps	 $0x39, %xmm10, %xmm10
	mulps	 %xmm15, %xmm10
	addps	 %xmm10,  %xmm0

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, A2
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(307):
	testq	$2, MM
	je	L(308)

	movsd	 -32 * SIZE(A1), %xmm4
	movsd	 -32 * SIZE(A1, LDA), %xmm5
	movsd	 -32 * SIZE(A2), %xmm6
	movsd	 -32 * SIZE(A2, LDA), %xmm7

	movsd	 -32 * SIZE(Y1), %xmm0

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm0
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm0
	mulps	 %xmm15, %xmm7
	addps	 %xmm7,  %xmm0

	movlps	 %xmm0, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(308):
	testq	$1, MM
	je	L(309)

	movss	 -32 * SIZE(Y1), %xmm0

	movss	 -32 * SIZE(A1), %xmm4
	movss	 -32 * SIZE(A1, LDA), %xmm5
	movss	 -32 * SIZE(A2), %xmm6
	movss	 -32 * SIZE(A2, LDA), %xmm7

	mulss	 %xmm12, %xmm4
	addss	 %xmm4,  %xmm0
	mulss	 %xmm13, %xmm5
	addss	 %xmm5,  %xmm0
	mulss	 %xmm14, %xmm6
	addss	 %xmm6,  %xmm0
	mulss	 %xmm15, %xmm7
	addss	 %xmm7,  %xmm0

	movss	 %xmm0, -32 * SIZE(Y1)
	ALIGN_3

L(309):
	cmpq	$4, N
	jge	L(301)
	ALIGN_3

L(310):
	cmpq	$3, N
	jne	L(320)

	leaq	32 * SIZE(BUFFER), Y1
	movq	A,  A1
	leaq	(A,  LDA, 2), A2
	leaq	(A,  LDA, 4), A

	movss	(X), %xmm12
	addq	INCX, X
	movss	(X), %xmm13
	addq	INCX, X
	movss	(X), %xmm14
	addq	INCX, X

	movss	ALPHA, %xmm0

	mulss	%xmm0, %xmm12
	mulss	%xmm0, %xmm13
	mulss	%xmm0, %xmm14

	shufps	$0, %xmm12, %xmm12
	shufps	$0, %xmm13, %xmm13
	shufps	$0, %xmm14, %xmm14

	cmpq	$3, M
	jle	L(317)

	testq	$SIZE, A1
	je	L(31X)

	movss	 -32 * SIZE(A1), %xmm0
	movss	 -32 * SIZE(A1, LDA), %xmm1
	movss	 -32 * SIZE(A2), %xmm2

	movss	 -32 * SIZE(Y1), %xmm8

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8
	mulss	 %xmm13, %xmm1
	addss	 %xmm1,  %xmm8
	mulss	 %xmm14, %xmm2
	addss	 %xmm2, %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, Y1
	ALIGN_3

L(31X):
	testq	$2 * SIZE, A1
	je	L(31XX)

	movsd	 -32 * SIZE(A1), %xmm0
	movsd	 -32 * SIZE(A1, LDA), %xmm1
	movsd	 -32 * SIZE(A2), %xmm2

	movsd	 -32 * SIZE(Y1), %xmm8

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8
	mulps	 %xmm13, %xmm1
	addps	 %xmm1,  %xmm8
	mulps	 %xmm14, %xmm2
	addps	 %xmm2, %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(31XX):
	movaps	-35 * SIZE(A1, LDA), %xmm8
	movaps	-34 * SIZE(A2), %xmm9
	movaps	-33 * SIZE(A2, LDA), %xmm10

	movq	MM,  I
	sarq	$4,  I
	jle	L(315)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -28 * SIZE(A1), %xmm5
	movaps	 -24 * SIZE(A1), %xmm6

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YL1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YL1(-20 * SIZE, Y1, %xmm3)

	decq	 I
	jle	 L(314)
	ALIGN_3

L(313):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A1)
#endif

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A1), %xmm7
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -31 * SIZE(A1, LDA), %xmm4
	mulps	 %xmm12, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -27 * SIZE(A1, LDA), %xmm5
	mulps	 %xmm12, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	 -23 * SIZE(A1, LDA), %xmm6

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A1, LDA)
#endif

	movss	 %xmm4,  %xmm8
	shufps	 $0x93, %xmm4, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8, %xmm0
	movaps	 -19 * SIZE(A1, LDA), %xmm8
	movss	 %xmm5,  %xmm4
	shufps	 $0x93, %xmm5, %xmm4
	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -30 * SIZE(A2), %xmm4
	movss	 %xmm6,  %xmm5
	shufps	 $0x93, %xmm6, %xmm5
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -26 * SIZE(A2), %xmm5
	movss	 %xmm8,  %xmm6
	shufps	 $0x93, %xmm8, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	 -22 * SIZE(A2), %xmm6

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A2)
#endif

	shufps	 $0x4e, %xmm4, %xmm9
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm0
	movaps	 -18 * SIZE(A2), %xmm9
	shufps	 $0x4e, %xmm5, %xmm4
	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -16 * SIZE(A1), %xmm4
	shufps	 $0x4e, %xmm6, %xmm5
	mulps	 %xmm14, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -12 * SIZE(A1), %xmm5
	shufps	 $0x4e, %xmm9, %xmm6
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	  -8 * SIZE(A1), %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(Y1)
#endif

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	MOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-12 * SIZE, Y1, %xmm1)
	MOVUPS_YL1( -8 * SIZE, Y1, %xmm2)
	MOVUPS_YL1( -4 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(313)
	ALIGN_3

L(314):
	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A1), %xmm7
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -31 * SIZE(A1, LDA), %xmm4

	mulps	 %xmm12, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -27 * SIZE(A1, LDA), %xmm5
	mulps	 %xmm12, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	 -23 * SIZE(A1, LDA), %xmm6

	movss	 %xmm4,  %xmm8
	shufps	 $0x93, %xmm4, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8, %xmm0
	movaps	 -19 * SIZE(A1, LDA), %xmm8
	movss	 %xmm5,  %xmm4
	shufps	 $0x93, %xmm5, %xmm4
	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -30 * SIZE(A2), %xmm4

	movss	 %xmm6,  %xmm5
	shufps	 $0x93, %xmm6, %xmm5
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -26 * SIZE(A2), %xmm5
	movss	 %xmm8,  %xmm6
	shufps	 $0x93, %xmm8, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	 -22 * SIZE(A2), %xmm6

	shufps	 $0x4e, %xmm4, %xmm9
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm0
	movaps	 -18 * SIZE(A2), %xmm9
	shufps	 $0x4e, %xmm5, %xmm4
	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm1
	shufps	 $0x4e, %xmm6, %xmm5
	mulps	 %xmm14, %xmm5
	addps	 %xmm5,  %xmm2
	shufps	 $0x4e, %xmm9, %xmm6
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm3

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1
	ALIGN_3

L(315):
	testq	$8, MM
	je	L(316)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -28 * SIZE(A1), %xmm5

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -31 * SIZE(A1, LDA), %xmm6
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -27 * SIZE(A1, LDA), %xmm7

	movss	 %xmm6,  %xmm8
	shufps	 $0x93, %xmm6, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8, %xmm0
	movaps	 -30 * SIZE(A2), %xmm4
	movss	 %xmm7,  %xmm6
	shufps	 $0x93, %xmm7, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm1
	movaps	 %xmm7, %xmm8
	movaps	 -26 * SIZE(A2), %xmm5

	shufps	 $0x4e, %xmm4, %xmm9
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm0
	shufps	 $0x4e, %xmm5, %xmm4
	mulps	 %xmm14, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 %xmm5, %xmm9

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)

	addq	 $8 * SIZE, A1
	addq	 $8 * SIZE, A2
	addq	 $8 * SIZE, Y1
	ALIGN_3

L(316):
	testq	$4, MM
	je	L(317)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -31 * SIZE(A1, LDA), %xmm5

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -30 * SIZE(A2), %xmm6
	movss	 %xmm5,  %xmm8
	shufps	 $0x93, %xmm5, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8, %xmm0

	shufps	 $0x4e, %xmm6, %xmm9
	mulps	 %xmm14, %xmm9
	addps	 %xmm9,  %xmm0

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, A2
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(317):
	testq	$2, MM
	je	L(318)

	movsd	 -32 * SIZE(A1), %xmm4
	movsd	 -32 * SIZE(A1, LDA), %xmm5
	movsd	 -32 * SIZE(A2), %xmm6

	movsd	 -32 * SIZE(Y1), %xmm0

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm0
	mulps	 %xmm14, %xmm6
	addps	 %xmm6,  %xmm0

	movlps	 %xmm0, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(318):
	testq	$1, MM
	je	L(990)

	movss	 -32 * SIZE(Y1), %xmm0

	movss	 -32 * SIZE(A1), %xmm4
	movss	 -32 * SIZE(A1, LDA), %xmm5
	movss	 -32 * SIZE(A2), %xmm6

	mulss	 %xmm12, %xmm4
	addss	 %xmm4,  %xmm0
	mulss	 %xmm13, %xmm5
	addss	 %xmm5,  %xmm0
	mulss	 %xmm14, %xmm6
	addss	 %xmm6,  %xmm0

	movss	 %xmm0, -32 * SIZE(Y1)
	jmp	 L(990)
	ALIGN_3

L(320):
	cmpq	$2, N
	jne	L(330)

	leaq	32 * SIZE(BUFFER), Y1
	movq	A,  A1
	leaq	(A,  LDA, 1), A2

	movss	(X), %xmm12
	addq	INCX, X
	movss	(X), %xmm13
	addq	INCX, X

	movss	ALPHA, %xmm0

	mulss	%xmm0, %xmm12
	mulss	%xmm0, %xmm13

	shufps	$0, %xmm12, %xmm12
	shufps	$0, %xmm13, %xmm13

	cmpq	$3, M
	jle	L(327)

	testq	$SIZE, A1
	je	L(32X)

	movss	 -32 * SIZE(Y1), %xmm9

	movss	 -32 * SIZE(A1), %xmm0
	movss	 -32 * SIZE(A2), %xmm1

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm9
	mulss	 %xmm13, %xmm1
	addss	 %xmm1,  %xmm9

	movss	 %xmm9, -32 * SIZE(Y1)

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, Y1
	ALIGN_3

L(32X):
	testq	$2 * SIZE, A1
	je	L(32XX)

	movsd	 -32 * SIZE(A1), %xmm0
	movsd	 -32 * SIZE(A2), %xmm1

	movsd	 -32 * SIZE(Y1), %xmm8

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8
	mulps	 %xmm13, %xmm1
	addps	 %xmm1,  %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(32XX):
	movaps	-35 * SIZE(A1, LDA), %xmm8

	movq	MM,  I
	sarq	$4,  I
	jle	L(325)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -28 * SIZE(A1), %xmm5
	movaps	 -24 * SIZE(A1), %xmm6

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YL1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YL1(-20 * SIZE, Y1, %xmm3)

	decq	 I
	jle	 L(324)
	ALIGN_3

L(323):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A1)
#endif

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A1), %xmm7
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -31 * SIZE(A2), %xmm4
	mulps	 %xmm12, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -27 * SIZE(A2), %xmm5
	mulps	 %xmm12, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	 -23 * SIZE(A2), %xmm6

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A2)
#endif

	movss	 %xmm4,  %xmm8
	shufps	 $0x93, %xmm4, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8, %xmm0
	movaps	 -19 * SIZE(A2), %xmm8
	movss	 %xmm5,  %xmm4
	shufps	 $0x93, %xmm5, %xmm4
	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm1
	movaps	 -16 * SIZE(A1), %xmm4
	movss	 %xmm6,  %xmm5
	shufps	 $0x93, %xmm6, %xmm5
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm2
	movaps	 -12 * SIZE(A1), %xmm5
	movss	 %xmm8,  %xmm6
	shufps	 $0x93, %xmm8, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm3
	movaps	  -8 * SIZE(A1), %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(Y1)
#endif

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	MOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-12 * SIZE, Y1, %xmm1)
	MOVUPS_YL1( -8 * SIZE, Y1, %xmm2)
	MOVUPS_YL1( -4 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(323)
	ALIGN_3

L(324):
	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -20 * SIZE(A1), %xmm7
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -31 * SIZE(A2), %xmm4

	mulps	 %xmm12, %xmm6
	addps	 %xmm6,  %xmm2
	movaps	 -27 * SIZE(A2), %xmm5
	mulps	 %xmm12, %xmm7
	addps	 %xmm7,  %xmm3
	movaps	 -23 * SIZE(A2), %xmm6

	movss	 %xmm4,  %xmm8
	shufps	 $0x93, %xmm4, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8, %xmm0
	movaps	 -19 * SIZE(A2), %xmm8
	movss	 %xmm5,  %xmm4
	shufps	 $0x93, %xmm5, %xmm4
	mulps	 %xmm13, %xmm4
	addps	 %xmm4,  %xmm1

	movss	 %xmm6,  %xmm5
	shufps	 $0x93, %xmm6, %xmm5
	mulps	 %xmm13, %xmm5
	addps	 %xmm5,  %xmm2
	movss	 %xmm8,  %xmm6
	shufps	 $0x93, %xmm8, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm3

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, A2
	subq	 $-16 * SIZE, Y1
	ALIGN_3

L(325):
	testq	$8, MM
	je	L(326)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -28 * SIZE(A1), %xmm5

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movaps	 -31 * SIZE(A2), %xmm6
	mulps	 %xmm12, %xmm5
	addps	 %xmm5,  %xmm1
	movaps	 -27 * SIZE(A2), %xmm7

	movss	 %xmm6,  %xmm8
	shufps	 $0x93, %xmm6, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8, %xmm0
	movss	 %xmm7,  %xmm6
	shufps	 $0x93, %xmm7, %xmm6
	mulps	 %xmm13, %xmm6
	addps	 %xmm6,  %xmm1
	movaps	 %xmm7, %xmm8

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)

	addq	 $8 * SIZE, A1
	addq	 $8 * SIZE, A2
	addq	 $8 * SIZE, Y1
	ALIGN_3

L(326):
	testq	$4, MM
	je	L(327)

	movaps	 -32 * SIZE(A1), %xmm4
	movaps	 -31 * SIZE(A2), %xmm5

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)

	mulps	 %xmm12, %xmm4
	addps	 %xmm4,  %xmm0
	movss	 %xmm5,  %xmm8
	shufps	 $0x93, %xmm5, %xmm8
	mulps	 %xmm13, %xmm8
	addps	 %xmm8, %xmm0

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, A2
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(327):
	testq	$2, MM
	je	L(328)

	movsd	 -32 * SIZE(A1), %xmm0
	movsd	 -32 * SIZE(A2), %xmm1

	movsd	 -32 * SIZE(Y1), %xmm8

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8
	mulps	 %xmm13, %xmm1
	addps	 %xmm1,  %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(328):
	testq	$1, MM
	je	L(990)

	movss	 -32 * SIZE(Y1), %xmm8

	movss	 -32 * SIZE(A1), %xmm0
	movss	 -32 * SIZE(A2), %xmm1

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8
	mulss	 %xmm13, %xmm1
	addss	 %xmm1,  %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)
	jmp	 L(990)
	ALIGN_3

L(330):
	cmpq	$1, N
	jne	L(990)

	leaq	32 * SIZE(BUFFER), Y1
	movq	A,  A1

	movss	(X), %xmm12

	mulss	ALPHA, %xmm12
	shufps	$0, %xmm12, %xmm12

	cmpq	$3, M
	jle	L(337)

	testq	$SIZE, A1
	je	L(33X)

	movss	 -32 * SIZE(Y1), %xmm8
	movss	 -32 * SIZE(A1), %xmm0

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, Y1
	ALIGN_3

L(33X):
	testq	$2 * SIZE, A1
	je	L(33XX)

	movsd	 -32 * SIZE(Y1), %xmm8
	movsd	 -32 * SIZE(A1), %xmm0

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(33XX):
	movq	MM,  I
	sarq	$4,  I
	jle	L(335)

	movaps	-32 * SIZE(A1), %xmm8
	movaps	-28 * SIZE(A1), %xmm9
	movaps	-24 * SIZE(A1), %xmm10
	movaps	-20 * SIZE(A1), %xmm11

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YL1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YL1(-20 * SIZE, Y1, %xmm3)

	decq	 I
	jle	 L(334)
	ALIGN_3

L(333):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(A1)
#endif

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	movaps	-16 * SIZE(A1), %xmm8
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	movaps	-12 * SIZE(A1), %xmm9
	mulps	 %xmm12, %xmm10
	addps	 %xmm10,  %xmm2
	movaps	 -8 * SIZE(A1), %xmm10
	mulps	 %xmm12, %xmm11
	addps	 %xmm11,  %xmm3
	movaps	 -4 * SIZE(A1), %xmm11

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(Y1)
#endif

	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	MOVUPS_YL1(-16 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-12 * SIZE, Y1, %xmm1)
	MOVUPS_YL1( -8 * SIZE, Y1, %xmm2)
	MOVUPS_YL1( -4 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, Y1

	subq	 $1, I
	BRANCH
	jg	L(333)
	ALIGN_3

L(334):
	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)
	mulps	 %xmm12, %xmm10
	addps	 %xmm10,  %xmm2
	MOVUPS_YS1(-24 * SIZE, Y1, %xmm2)
	mulps	 %xmm12, %xmm11
	addps	 %xmm11,  %xmm3
	MOVUPS_YS1(-20 * SIZE, Y1, %xmm3)

	subq	 $-16 * SIZE, A1
	subq	 $-16 * SIZE, Y1
	ALIGN_3

L(335):
	testq	$8, MM
	je	L(336)

	movaps	-32 * SIZE(A1), %xmm8
	movaps	-28 * SIZE(A1), %xmm9

	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)
	MOVUPS_YL1(-28 * SIZE, Y1, %xmm1)

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)
	mulps	 %xmm12, %xmm9
	addps	 %xmm9,  %xmm1
	MOVUPS_YS1(-28 * SIZE, Y1, %xmm1)

	addq	 $8 * SIZE, A1
	addq	 $8 * SIZE, Y1
	ALIGN_3

L(336):
	testq	$4, MM
	je	L(337)

	movaps	-32 * SIZE(A1), %xmm8
	MOVUPS_YL1(-32 * SIZE, Y1, %xmm0)

	mulps	 %xmm12, %xmm8
	addps	 %xmm8,  %xmm0
	MOVUPS_YS1(-32 * SIZE, Y1, %xmm0)

	addq	 $4 * SIZE, A1
	addq	 $4 * SIZE, Y1
	ALIGN_3

L(337):
	testq	$2, MM
	je	L(338)

	movsd	 -32 * SIZE(Y1), %xmm8
	movsd	 -32 * SIZE(A1), %xmm0

	mulps	 %xmm12, %xmm0
	addps	 %xmm0,  %xmm8

	movlps	 %xmm8, -32 * SIZE(Y1)

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, Y1
	ALIGN_3

L(338):
	testq	$1, MM
	je	L(990)

	movss	 -32 * SIZE(Y1), %xmm8
	movss	 -32 * SIZE(A1), %xmm0

	mulss	 %xmm12, %xmm0
	addss	 %xmm0,  %xmm8

	movss	 %xmm8, -32 * SIZE(Y1)
	jmp	L(990)
#endif
	ALIGN_4


L(990):
	movq	Y,  Y1

	movq	M,  %rax
	sarq	$3,  %rax
	jle	L(994)
	ALIGN_3
L(992):
	movsd	0 * SIZE(BUFFER), %xmm0
	movhps	2 * SIZE(BUFFER), %xmm0
	movsd	4 * SIZE(BUFFER), %xmm4
	movhps	6 * SIZE(BUFFER), %xmm4

	pshufd	$0x01, %xmm0, %xmm1
	pshufd	$0x02, %xmm0, %xmm2
	pshufd	$0x03, %xmm0, %xmm3

	pshufd	$0x01, %xmm4, %xmm5
	pshufd	$0x02, %xmm4, %xmm6
	pshufd	$0x03, %xmm4, %xmm7

	addss	(Y), %xmm0
	addq	INCY, Y
	addss	(Y), %xmm1
	addq	INCY, Y
	addss	(Y), %xmm2
	addq	INCY, Y
	addss	(Y), %xmm3
	addq	INCY, Y
	addss	(Y), %xmm4
	addq	INCY, Y
	addss	(Y), %xmm5
	addq	INCY, Y
	addss	(Y), %xmm6
	addq	INCY, Y
	addss	(Y), %xmm7
	addq	INCY, Y

	movss	%xmm0, (Y1)
	addq	INCY, Y1
	movss	%xmm1, (Y1)
	addq	INCY, Y1
	movss	%xmm2, (Y1)
	addq	INCY, Y1
	movss	%xmm3, (Y1)
	addq	INCY, Y1
	movss	%xmm4, (Y1)
	addq	INCY, Y1
	movss	%xmm5, (Y1)
	addq	INCY, Y1
	movss	%xmm6, (Y1)
	addq	INCY, Y1
	movss	%xmm7, (Y1)
	addq	INCY, Y1

	addq	$8 * SIZE, BUFFER
	decq	%rax
	jg	L(992)
	ALIGN_3

L(994):
	testq	$7, M
	jle	L(999)

	testq	$4, M
	jle	L(995)

	movsd	0 * SIZE(BUFFER), %xmm0
	movhps	2 * SIZE(BUFFER), %xmm0

	pshufd	$0x01, %xmm0, %xmm1
	pshufd	$0x02, %xmm0, %xmm2
	pshufd	$0x03, %xmm0, %xmm3

	addss	(Y), %xmm0
	addq	INCY, Y
	addss	(Y), %xmm1
	addq	INCY, Y
	addss	(Y), %xmm2
	addq	INCY, Y
	addss	(Y), %xmm3
	addq	INCY, Y

	movss	%xmm0, (Y1)
	addq	INCY, Y1
	movss	%xmm1, (Y1)
	addq	INCY, Y1
	movss	%xmm2, (Y1)
	addq	INCY, Y1
	movss	%xmm3, (Y1)
	addq	INCY, Y1

	addq	$4 * SIZE, BUFFER
	ALIGN_3

L(995):
	testq	$2, M
	jle	L(996)

	movsd	(BUFFER), %xmm0

	pshufd	$0x01, %xmm0, %xmm1

	addss	(Y), %xmm0
	addq	INCY, Y
	addss	(Y), %xmm1
	addq	INCY, Y

	movss	%xmm0, (Y1)
	addq	INCY, Y1
	movss	%xmm1, (Y1)
	addq	INCY, Y1

	addq	$2 * SIZE, BUFFER
	ALIGN_3

L(996):
	testq	$1, M
	jle	L(999)

	movss	(BUFFER), %xmm0

	addss	(Y), %xmm0

	movss	%xmm0, (Y1)
	ALIGN_3

L(999):
	leaq	(,M,SIZE),%rax
	addq	%rax,AA
	jmp		L(0t)
	ALIGN_4

L(999x):
	movq	  0(%rsp), %rbx
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
