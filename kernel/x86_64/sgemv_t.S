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
#define MMM		56(%rsp)
#define NN		64(%rsp)
#define AA		72(%rsp)
#define LDAX	80(%rsp)

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
#define MMM	232(%rsp)
#define NN	240(%rsp)
#define AA	248(%rsp)
#define LDAX 256(%rsp)

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
#define X1	%rbp

#define Y1	INCX

#ifdef ALIGNED_ACCESS
#define MM	%r15
#else
#define MM	M
#endif

#define ALPHA	%xmm7

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

	movq	OLD_M,	      MMM
	movq	OLD_N,        NN
	movq	OLD_A,        X
	movq	X,	      AA
	movq	OLD_LDA,      X
	movq	X,	      LDAX
	movq	OLD_X,        X
#else
	movq	OLD_M,	      MMM
	movq	OLD_N,        NN
	movq	OLD_A,        AA
	movq	OLD_LDA,      LDAX
#endif
#ifndef WINDOWS_ABI
	pshufd	$0, %xmm0, ALPHA
#else
	pshufd	$0, %xmm3, ALPHA
#endif


L(0t):
	xorq	M,M
	addq	$1,M
	salq	$22,M
	subq	M,MMM
	jge		L(00t)
	ALIGN_4

	movq	MMM,%rax
	addq	M,%rax
	jle		L(999x)
	movq	%rax,M

L(00t):
	movq	LDAX,LDA
	movq	NN,N
	movq	AA,A
	movq	STACK_INCX,   INCX
	movq	STACK_Y,      Y
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

	testq	M, M
	jle	L(999)
	testq	N, N
	jle	L(999)
	ALIGN_4

	subq	$-32 * SIZE, A

#ifdef ALIGNED_ACCESS
	movq	A,  %rax
	andq	$4 * SIZE - 1, %rax
	addq	%rax, BUFFER
#endif

	movq	BUFFER, X1

	movq	M,  I
	sarq	$3, I
	jle	L(05)
	ALIGN_4

L(02):
	movss	(X), %xmm0
	addq	INCX, X
	movss	(X), %xmm1
	addq	INCX, X

	movss	(X), %xmm2
	addq	INCX, X
	movss	(X), %xmm3
	addq	INCX, X

	movss	(X), %xmm4
	addq	INCX, X
	movss	(X), %xmm5
	addq	INCX, X

	movss	(X), %xmm6
	addq	INCX, X
	movss	(X), %xmm8
	addq	INCX, X

	movss	%xmm0, 0 * SIZE(X1)
	movss	%xmm1, 1 * SIZE(X1)
	movss	%xmm2, 2 * SIZE(X1)
	movss	%xmm3, 3 * SIZE(X1)
	movss	%xmm4, 4 * SIZE(X1)
	movss	%xmm5, 5 * SIZE(X1)
	movss	%xmm6, 6 * SIZE(X1)
	movss	%xmm8, 7 * SIZE(X1)

	addq	$8 * SIZE, X1
	decq	I
	jg	L(02)
	ALIGN_4

L(05):
	movq	M,  I
	andq	$7, I
	jle	L(10)
	ALIGN_2

L(06):
	movss	(X), %xmm0
	addq	INCX, X
	movss	%xmm0, 0 * SIZE(X1)
	addq	$SIZE, X1
	decq	I
	jg	L(06)
	ALIGN_4

L(10):
	movq	Y, Y1

#ifdef ALIGNED_ACCESS
	testq	$4 * SIZE - 1, LDA
	jne	L(100)
#endif

#if GEMV_UNROLL >= 8

	cmpq	$8, N
	jl	L(20)
	ALIGN_3

L(11):
	subq	$8, N

	leaq	32 * SIZE(BUFFER), X1

	movq	A, A1
	leaq	(A1, LDA, 4), A2
	leaq	(A1, LDA, 8), A

	xorps	%xmm8, %xmm8
	xorps	%xmm9, %xmm9
	xorps	%xmm10, %xmm10
	xorps	%xmm11, %xmm11
	xorps	%xmm12, %xmm12
	xorps	%xmm13, %xmm13
	xorps	%xmm14, %xmm14
	xorps	%xmm15, %xmm15

#ifdef ALIGNED_ACCESS
	cmpq	$3, M
	jle	L(17)

	testq	$SIZE, A1
	je	L(1X)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A1, LDA, 1), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	movss	-32 * SIZE(A1, LDA, 2), %xmm2
	mulss	%xmm4, %xmm2
	addss	%xmm2, %xmm10
	movss	-32 * SIZE(A1, LDA3, 1), %xmm3
	mulss	%xmm4, %xmm3
	addss	%xmm3, %xmm11
	movss	-32 * SIZE(A2), %xmm0
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm12
	movss	-32 * SIZE(A2, LDA, 1), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm13
	movss	-32 * SIZE(A2, LDA, 2), %xmm2
	mulss	%xmm4, %xmm2
	addss	%xmm2, %xmm14
	movss	-32 * SIZE(A2, LDA3, 1), %xmm3
	mulss	%xmm4, %xmm3
	addss	%xmm3, %xmm15

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, X1
	ALIGN_3

L(1X):
	testq	$2 * SIZE, A1
	je	L(1XX)

#ifdef movsd
	xorps	%xmm0, %xmm0
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(A1), %xmm0
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A1, LDA, 1), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
#ifdef movsd
	xorps	%xmm2, %xmm2
#endif
	movsd	-32 * SIZE(A1, LDA, 2), %xmm2
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
#ifdef movsd
	xorps	%xmm3, %xmm3
#endif
	movsd	-32 * SIZE(A1, LDA3, 1), %xmm3
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm11
	movsd	-32 * SIZE(A2), %xmm0
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm12
	movsd	-32 * SIZE(A2, LDA, 1), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm13
	movsd	-32 * SIZE(A2, LDA, 2), %xmm2
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm14
	movsd	-32 * SIZE(A2, LDA3, 1), %xmm3
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm15

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, X1
	ALIGN_3

L(1XX):
#endif

	MOVUPS_XL1 (-32 * SIZE, X1, %xmm4)
	MOVUPS_XL1 (-28 * SIZE, X1, %xmm5)

#ifdef PREFETCHW
	PREFETCHW	8 * SIZE(Y1)
#endif

	movq	MM,  I
	sarq	$4,  I
	jle	L(15)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-32 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A2 (-32 * SIZE, A1, LDA, 2, %xmm2)
	MOVUPS_A2 (-32 * SIZE, A1, LDA3, 1, %xmm3)

	decq	I
	jle	L(13)
	ALIGN_4

L(12):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A1)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-32 * SIZE, A2, %xmm0)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-32 * SIZE, A2, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A2 (-32 * SIZE, A2, LDA, 2, %xmm2)
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm11
	MOVUPS_A2 (-32 * SIZE, A2, LDA3, 1, %xmm3)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A1, LDA, 1)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm12
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm13
	MOVUPS_A2 (-28 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm14
	MOVUPS_A2 (-28 * SIZE, A1, LDA, 2, %xmm2)
	mulps	%xmm4, %xmm3
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm3, %xmm15
	MOVUPS_A2 (-28 * SIZE, A1, LDA3, 1, %xmm3)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A1, LDA, 2)
#endif

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A2, %xmm0)
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-28 * SIZE, A2, LDA, 1, %xmm1)
	mulps	%xmm5, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A2 (-28 * SIZE, A2, LDA, 2, %xmm2)
	mulps	%xmm5, %xmm3
	addps	%xmm3, %xmm11
	MOVUPS_A2 (-28 * SIZE, A2, LDA3, 1, %xmm3)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A1, LDA3)
#endif

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm12
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm13
	MOVUPS_A2 (-24 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm5, %xmm2
	addps	%xmm2, %xmm14
	MOVUPS_A2 (-24 * SIZE, A1, LDA, 2, %xmm2)
	mulps	%xmm5, %xmm3
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm3, %xmm15
	MOVUPS_A2 (-24 * SIZE, A1, LDA3, 1, %xmm3)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A2)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A2, %xmm0)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-24 * SIZE, A2, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A2 (-24 * SIZE, A2, LDA, 2, %xmm2)
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm11
	MOVUPS_A2 (-24 * SIZE, A2, LDA3, 1, %xmm3)
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm12

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A2, LDA, 1)
#endif

	MOVUPS_A1 (-20 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm13
	MOVUPS_A2 (-20 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm14
	MOVUPS_A2 (-20 * SIZE, A1, LDA, 2, %xmm2)
	mulps	%xmm4, %xmm3
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm3, %xmm15
	MOVUPS_A2 (-20 * SIZE, A1, LDA3, 1, %xmm3)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A2, LDA, 2)
#endif

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-20 * SIZE, A2, %xmm0)
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-20 * SIZE, A2, LDA, 1, %xmm1)
	mulps	%xmm5, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A2 (-20 * SIZE, A2, LDA, 2, %xmm2)
	mulps	%xmm5, %xmm3
	addps	%xmm3, %xmm11
	MOVUPS_A2 (-20 * SIZE, A2, LDA3, 1, %xmm3)
	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm12

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A2, LDA3)
#endif

	MOVUPS_A1 (-16 * SIZE, A1, %xmm0)
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm13
	MOVUPS_A2 (-16 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm5, %xmm2
	addps	%xmm2, %xmm14

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(X1)
#endif

	MOVUPS_A2 (-16 * SIZE, A1, LDA, 2, %xmm2)
	mulps	%xmm5, %xmm3
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm3, %xmm15
	MOVUPS_A2 (-16 * SIZE, A1, LDA3, 1, %xmm3)

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1

	decq	I
	jg	L(12)
	ALIGN_4

L(13):
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-32 * SIZE, A2, %xmm0)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-32 * SIZE, A2, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A2 (-32 * SIZE, A2, LDA, 2, %xmm2)
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm11
	MOVUPS_A2 (-32 * SIZE, A2, LDA3, 1, %xmm3)
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm12
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm13
	MOVUPS_A2 (-28 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm14
	MOVUPS_A2 (-28 * SIZE, A1, LDA, 2, %xmm2)
	mulps	%xmm4, %xmm3
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm3, %xmm15
	MOVUPS_A2 (-28 * SIZE, A1, LDA3, 1, %xmm3)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A2, %xmm0)
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-28 * SIZE, A2, LDA, 1, %xmm1)
	mulps	%xmm5, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A2 (-28 * SIZE, A2, LDA, 2, %xmm2)
	mulps	%xmm5, %xmm3
	addps	%xmm3, %xmm11
	MOVUPS_A2 (-28 * SIZE, A2, LDA3, 1, %xmm3)
	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm12
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm13
	MOVUPS_A2 (-24 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm5, %xmm2
	addps	%xmm2, %xmm14
	MOVUPS_A2 (-24 * SIZE, A1, LDA, 2, %xmm2)
	mulps	%xmm5, %xmm3
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm3, %xmm15
	MOVUPS_A2 (-24 * SIZE, A1, LDA3, 1, %xmm3)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A2, %xmm0)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-24 * SIZE, A2, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A2 (-24 * SIZE, A2, LDA, 2, %xmm2)
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm11
	MOVUPS_A2 (-24 * SIZE, A2, LDA3, 1, %xmm3)
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm12
	MOVUPS_A1 (-20 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm13
	MOVUPS_A2 (-20 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm14
	MOVUPS_A2 (-20 * SIZE, A1, LDA, 2, %xmm2)
	mulps	%xmm4, %xmm3
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm3, %xmm15
	MOVUPS_A2 (-20 * SIZE, A1, LDA3, 1, %xmm3)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-20 * SIZE, A2, %xmm0)
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-20 * SIZE, A2, LDA, 1, %xmm1)
	mulps	%xmm5, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A2 (-20 * SIZE, A2, LDA, 2, %xmm2)
	mulps	%xmm5, %xmm3
	addps	%xmm3, %xmm11
	MOVUPS_A2 (-20 * SIZE, A2, LDA3, 1, %xmm3)
	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm12
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm13
	mulps	%xmm5, %xmm2
	addps	%xmm2, %xmm14
	mulps	%xmm5, %xmm3
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm3, %xmm15

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1
	ALIGN_4

L(15):
	testq	$8, MM
	jle	L(16)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-32 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A2 (-32 * SIZE, A1, LDA, 2, %xmm2)
	MOVUPS_A2 (-32 * SIZE, A1, LDA3, 1, %xmm3)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-32 * SIZE, A2, %xmm0)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-32 * SIZE, A2, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A2 (-32 * SIZE, A2, LDA, 2, %xmm2)
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm11
	MOVUPS_A2 (-32 * SIZE, A2, LDA3, 1, %xmm3)
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm12
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm13
	MOVUPS_A2 (-28 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm14
	MOVUPS_A2 (-28 * SIZE, A1, LDA, 2, %xmm2)
	mulps	%xmm4, %xmm3
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm3, %xmm15
	MOVUPS_A2 (-28 * SIZE, A1, LDA3, 1, %xmm3)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A2, %xmm0)
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-28 * SIZE, A2, LDA, 1, %xmm1)
	mulps	%xmm5, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A2 (-28 * SIZE, A2, LDA, 2, %xmm2)
	mulps	%xmm5, %xmm3
	addps	%xmm3, %xmm11
	MOVUPS_A2 (-28 * SIZE, A2, LDA3, 1, %xmm3)
	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm12
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm13
	mulps	%xmm5, %xmm2
	addps	%xmm2, %xmm14
	mulps	%xmm5, %xmm3
	addps	%xmm3, %xmm15

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1
	ALIGN_4

L(16):
	testq	$4, MM
	jle	L(17)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-32 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A2 (-32 * SIZE, A1, LDA, 2, %xmm2)
	MOVUPS_A2 (-32 * SIZE, A1, LDA3, 1, %xmm3)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-32 * SIZE, A2, %xmm0)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-32 * SIZE, A2, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A2 (-32 * SIZE, A2, LDA, 2, %xmm2)
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm11
	MOVUPS_A2 (-32 * SIZE, A2, LDA3, 1, %xmm3)
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm12
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm13
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm14
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm15

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, A2
	addq	$4 * SIZE, X1
	ALIGN_4

L(17):
	testq	$2, MM
	jle	L(18)

#ifdef movsd
	xorps	%xmm0, %xmm0
#endif
	movsd	-32 * SIZE(A1), %xmm0
#ifdef movsd
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A1, LDA, 1), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
#ifdef movsd
	xorps	%xmm2, %xmm2
#endif
	movsd	-32 * SIZE(A1, LDA, 2), %xmm2
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
#ifdef movsd
	xorps	%xmm3, %xmm3
#endif
	movsd	-32 * SIZE(A1, LDA3, 1), %xmm3
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm11
	movsd	-32 * SIZE(A2), %xmm0
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm12
	movsd	-32 * SIZE(A2, LDA, 1), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm13
	movsd	-32 * SIZE(A2, LDA, 2), %xmm2
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm14
	movsd	-32 * SIZE(A2, LDA3, 1), %xmm3
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm15

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, A2
	addq	$2 * SIZE, X1
	ALIGN_4

L(18):
	testq	$1, MM
	jle	L(19)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A1, LDA, 1), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	movss	-32 * SIZE(A1, LDA, 2), %xmm2
	mulss	%xmm4, %xmm2
	addss	%xmm2, %xmm10
	movss	-32 * SIZE(A1, LDA3, 1), %xmm3
	mulss	%xmm4, %xmm3
	addss	%xmm3, %xmm11
	movss	-32 * SIZE(A2), %xmm0
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm12
	movss	-32 * SIZE(A2, LDA, 1), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm13
	movss	-32 * SIZE(A2, LDA, 2), %xmm2
	mulss	%xmm4, %xmm2
	addss	%xmm2, %xmm14
	movss	-32 * SIZE(A2, LDA3, 1), %xmm3
	mulss	%xmm4, %xmm3
	addss	%xmm3, %xmm15
	ALIGN_4

L(19):
#ifdef HAVE_SSE3
	haddps	%xmm9, %xmm8
	haddps	%xmm11, %xmm10
	haddps	%xmm10, %xmm8

	pshufd	$0x1, %xmm8, %xmm9
	pshufd	$0x2, %xmm8, %xmm10
	pshufd	$0x3, %xmm8, %xmm11

	haddps	%xmm13, %xmm12
	haddps	%xmm15, %xmm14
	haddps	%xmm14, %xmm12

	pshufd	$0x1, %xmm12, %xmm13
	pshufd	$0x2, %xmm12, %xmm14
	pshufd	$0x3, %xmm12, %xmm15
#else
	movaps	%xmm8, %xmm0
	unpcklps %xmm9, %xmm8
	unpckhps %xmm9, %xmm0

	movaps	%xmm10, %xmm1
	unpcklps %xmm11, %xmm10
	unpckhps %xmm11, %xmm1

	movaps	%xmm8, %xmm9
	unpcklps %xmm10, %xmm8
	unpckhps %xmm10, %xmm9

	movaps	%xmm0, %xmm10
	unpcklps %xmm1, %xmm0
	unpckhps %xmm1, %xmm10

	addps	%xmm9, %xmm8
	addps	%xmm0, %xmm10
	addps	%xmm10, %xmm8

	pshufd	$0x2, %xmm8, %xmm9
	pshufd	$0x1, %xmm8, %xmm10
	pshufd	$0x3, %xmm8, %xmm11

	movaps	%xmm12, %xmm0
	unpcklps %xmm13, %xmm12
	unpckhps %xmm13, %xmm0

	movaps	%xmm14, %xmm1
	unpcklps %xmm15, %xmm14
	unpckhps %xmm15, %xmm1

	movaps	%xmm12, %xmm13
	unpcklps %xmm14, %xmm12
	unpckhps %xmm14, %xmm13

	movaps	%xmm0, %xmm14
	unpcklps %xmm1, %xmm0
	unpckhps %xmm1, %xmm14

	addps	%xmm13, %xmm12
	addps	%xmm0, %xmm14
	addps	%xmm14, %xmm12

	pshufd	$0x2, %xmm12, %xmm13
	pshufd	$0x1, %xmm12, %xmm14
	pshufd	$0x3, %xmm12, %xmm15
#endif

	mulss	ALPHA, %xmm8
	mulss	ALPHA, %xmm9
	mulss	ALPHA, %xmm10
	mulss	ALPHA, %xmm11
	mulss	ALPHA, %xmm12
	mulss	ALPHA, %xmm13
	mulss	ALPHA, %xmm14
	mulss	ALPHA, %xmm15

	addss	(Y), %xmm8
	addq	INCY, Y
	addss	(Y), %xmm9
	addq	INCY, Y
	addss	(Y), %xmm10
	addq	INCY, Y
	addss	(Y), %xmm11
	addq	INCY, Y
	addss	(Y), %xmm12
	addq	INCY, Y
	addss	(Y), %xmm13
	addq	INCY, Y
	addss	(Y), %xmm14
	addq	INCY, Y
	addss	(Y), %xmm15
	addq	INCY, Y

	movss	%xmm8,  (Y1)
	addq	INCY, Y1
	movss	%xmm9,  (Y1)
	addq	INCY, Y1
	movss	%xmm10, (Y1)
	addq	INCY, Y1
	movss	%xmm11, (Y1)
	addq	INCY, Y1
	movss	%xmm12, (Y1)
	addq	INCY, Y1
	movss	%xmm13, (Y1)
	addq	INCY, Y1
	movss	%xmm14, (Y1)
	addq	INCY, Y1
	movss	%xmm15, (Y1)
	addq	INCY, Y1

	cmpq	$8, N
	jge	L(11)
	ALIGN_4

L(20):
#endif

	cmpq	$4, N
	jl	L(30)

#if GEMV_UNROLL == 4
	ALIGN_3

L(21):
#endif
	subq	$4, N

	leaq	32 * SIZE(BUFFER), X1

	movq	A, A1
	leaq	(A1, LDA, 2), A2
	leaq	(A1, LDA, 4), A

	xorps	%xmm8, %xmm8
	xorps	%xmm9, %xmm9
	xorps	%xmm10, %xmm10
	xorps	%xmm11, %xmm11

#ifdef ALIGNED_ACCESS
	cmpq	$3, M
	jle	L(27)

	testq	$SIZE, A1
	je	L(2X)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A1, LDA), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	movss	-32 * SIZE(A2), %xmm2
	mulss	%xmm4, %xmm2
	addss	%xmm2, %xmm10
	movss	-32 * SIZE(A2, LDA), %xmm3
	mulss	%xmm4, %xmm3
	addss	%xmm3, %xmm11

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, X1
	ALIGN_3

L(2X):
	testq	$2 * SIZE, A1
	je	L(2XX)

#ifdef movsd
	xorps	%xmm0, %xmm0
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(A1), %xmm0
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A1, LDA), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
#ifdef movsd
	xorps	%xmm2, %xmm2
#endif
	movsd	-32 * SIZE(A2), %xmm2
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
#ifdef movsd
	xorps	%xmm3, %xmm3
#endif
	movsd	-32 * SIZE(A2, LDA), %xmm3
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm11

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, X1
	ALIGN_3

L(2XX):
#endif

	MOVUPS_XL1 (-32 * SIZE, X1, %xmm4)
	MOVUPS_XL1 (-28 * SIZE, X1, %xmm5)

#if (GEMV_UNROLL == 4) && defined(PREFETCHW)
	PREFETCHW	4 * SIZE(Y1)
#endif

	movq	MM,  I
	sarq	$4,  I
	jle	L(25)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-32 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A1 (-32 * SIZE, A2, %xmm2)
	MOVUPS_A2 (-32 * SIZE, A2, LDA, 1, %xmm3)

	MOVUPS_A1 (-28 * SIZE, A1, %xmm12)
	MOVUPS_A2 (-28 * SIZE, A1, LDA, 1, %xmm13)
	MOVUPS_A1 (-28 * SIZE, A2, %xmm14)
	MOVUPS_A2 (-28 * SIZE, A2, LDA, 1, %xmm15)

	decq	I
	jle	L(23)
	ALIGN_4

L(22):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-24 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-24 * SIZE, A2, %xmm2)
	mulps	%xmm4, %xmm3
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm3, %xmm11
	MOVUPS_A2 (-24 * SIZE, A2, LDA, 1, %xmm3)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1, LDA)
#endif

	mulps	%xmm5, %xmm12
	addps	%xmm12, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm12)
	mulps	%xmm5, %xmm13
	addps	%xmm13, %xmm9
	MOVUPS_A2 (-20 * SIZE, A1, LDA, 1, %xmm13)
	mulps	%xmm5, %xmm14
	addps	%xmm14, %xmm10
	MOVUPS_A1 (-20 * SIZE, A2, %xmm14)
	mulps	%xmm5, %xmm15
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm15, %xmm11
	MOVUPS_A2 (-20 * SIZE, A2, LDA, 1, %xmm15)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-16 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-16 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-16 * SIZE, A2, %xmm2)
	mulps	%xmm4, %xmm3
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm3, %xmm11
	MOVUPS_A2 (-16 * SIZE, A2, LDA, 1, %xmm3)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2, LDA)
#endif

	mulps	%xmm5, %xmm12
	addps	%xmm12, %xmm8
	MOVUPS_A1 (-12 * SIZE, A1, %xmm12)
	mulps	%xmm5, %xmm13
	addps	%xmm13, %xmm9
	MOVUPS_A2 (-12 * SIZE, A1, LDA, 1, %xmm13)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(X1)
#endif

	mulps	%xmm5, %xmm14
	addps	%xmm14, %xmm10
	MOVUPS_A1 (-12 * SIZE, A2, %xmm14)
	mulps	%xmm5, %xmm15
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm15, %xmm11
	MOVUPS_A2 (-12 * SIZE, A2, LDA, 1, %xmm15)

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1

	decq	I
	jg	L(22)
	ALIGN_4

L(23):
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-24 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-24 * SIZE, A2, %xmm2)
	mulps	%xmm4, %xmm3
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm3, %xmm11
	MOVUPS_A2 (-24 * SIZE, A2, LDA, 1, %xmm3)

	mulps	%xmm5, %xmm12
	addps	%xmm12, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm12)
	mulps	%xmm5, %xmm13
	addps	%xmm13, %xmm9
	MOVUPS_A2 (-20 * SIZE, A1, LDA, 1, %xmm13)
	mulps	%xmm5, %xmm14
	addps	%xmm14, %xmm10
	MOVUPS_A1 (-20 * SIZE, A2, %xmm14)
	mulps	%xmm5, %xmm15
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm15, %xmm11
	MOVUPS_A2 (-20 * SIZE, A2, LDA, 1, %xmm15)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
	mulps	%xmm4, %xmm3
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm3, %xmm11

	mulps	%xmm5, %xmm12
	addps	%xmm12, %xmm8
	mulps	%xmm5, %xmm13
	addps	%xmm13, %xmm9
	mulps	%xmm5, %xmm14
	addps	%xmm14, %xmm10
	mulps	%xmm5, %xmm15
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm15, %xmm11

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1
	ALIGN_4

L(25):
	testq	$8, MM
	jle	L(26)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A2 (-32 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A1 (-32 * SIZE, A2, %xmm2)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A2 (-32 * SIZE, A2, LDA, 1, %xmm3)
	mulps	%xmm4, %xmm3
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm3, %xmm11

	MOVUPS_A1 (-28 * SIZE, A1, %xmm12)
	mulps	%xmm5, %xmm12
	addps	%xmm12, %xmm8
	MOVUPS_A2 (-28 * SIZE, A1, LDA, 1, %xmm13)
	mulps	%xmm5, %xmm13
	addps	%xmm13, %xmm9
	MOVUPS_A1 (-28 * SIZE, A2, %xmm14)
	mulps	%xmm5, %xmm14
	addps	%xmm14, %xmm10
	MOVUPS_A2 (-28 * SIZE, A2, LDA, 1, %xmm15)
	mulps	%xmm5, %xmm15
	addps	%xmm15, %xmm11

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1
	ALIGN_4

L(26):
	testq	$4, MM
	jle	L(27)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A2 (-32 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A1 (-32 * SIZE, A2, %xmm2)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A2 (-32 * SIZE, A2, LDA, 1, %xmm3)
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm11

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, A2
	addq	$4 * SIZE, X1
	ALIGN_4

L(27):
	testq	$2, MM
	jle	L(28)

#ifdef movsd
	xorps	%xmm0, %xmm0
#endif
	movsd	-32 * SIZE(A1), %xmm0
#ifdef movsd
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A1, LDA), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
#ifdef movsd
	xorps	%xmm2, %xmm2
#endif
	movsd	-32 * SIZE(A2), %xmm2
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
#ifdef movsd
	xorps	%xmm3, %xmm3
#endif
	movsd	-32 * SIZE(A2, LDA), %xmm3
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm11
	shufps	$0xe, %xmm4, %xmm4

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, A2
	addq	$2 * SIZE, X1
	ALIGN_4

L(28):
	testq	$1, MM
	jle	L(29)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A1, LDA), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	movss	-32 * SIZE(A2), %xmm2
	mulss	%xmm4, %xmm2
	addss	%xmm2, %xmm10
	movss	-32 * SIZE(A2, LDA), %xmm3
	mulss	%xmm4, %xmm3
	addss	%xmm3, %xmm11
	ALIGN_4

L(29):
#ifdef HAVE_SSE3
	haddps	%xmm9, %xmm8
	haddps	%xmm11, %xmm10
	haddps	%xmm10, %xmm8

	pshufd	$0x1, %xmm8, %xmm9
	pshufd	$0x2, %xmm8, %xmm10
	pshufd	$0x3, %xmm8, %xmm11
#else
	movaps	%xmm8, %xmm0
	unpcklps %xmm9, %xmm8
	unpckhps %xmm9, %xmm0

	movaps	%xmm10, %xmm1
	unpcklps %xmm11, %xmm10
	unpckhps %xmm11, %xmm1

	movaps	%xmm8, %xmm9
	unpcklps %xmm10, %xmm8
	unpckhps %xmm10, %xmm9

	movaps	%xmm0, %xmm10
	unpcklps %xmm1, %xmm0
	unpckhps %xmm1, %xmm10

	addps	%xmm9, %xmm8
	addps	%xmm0, %xmm10
	addps	%xmm10, %xmm8

	pshufd	$0x2, %xmm8, %xmm9
	pshufd	$0x1, %xmm8, %xmm10
	pshufd	$0x3, %xmm8, %xmm11
#endif

	mulss	ALPHA, %xmm8
	mulss	ALPHA, %xmm9
	mulss	ALPHA, %xmm10
	mulss	ALPHA, %xmm11

	addss	(Y), %xmm8
	addq	INCY, Y
	addss	(Y), %xmm9
	addq	INCY, Y
	addss	(Y), %xmm10
	addq	INCY, Y
	addss	(Y), %xmm11
	addq	INCY, Y

	movss	%xmm8, (Y1)
	addq	INCY, Y1
	movss	%xmm9, (Y1)
	addq	INCY, Y1
	movss	%xmm10, (Y1)
	addq	INCY, Y1
	movss	%xmm11, (Y1)
	addq	INCY, Y1

#if GEMV_UNROLL == 4
	cmpq	$4, N
	jge	L(21)
#endif
	ALIGN_4

L(30):
	cmpq	$3, N
	jne	L(40)

	leaq	32 * SIZE(BUFFER), X1

	movq	A, A1
	leaq	(A1, LDA, 2), A2
	leaq	(A1, LDA, 4), A

	xorps	%xmm8, %xmm8
	xorps	%xmm9, %xmm9
	xorps	%xmm10, %xmm10

#ifdef ALIGNED_ACCESS
	cmpq	$3, M
	jle	L(37)

	testq	$SIZE, A1
	je	L(3X)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A1, LDA), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	movss	-32 * SIZE(A2), %xmm2
	mulss	%xmm4, %xmm2
	addss	%xmm2, %xmm10
	movss	-32 * SIZE(A2, LDA), %xmm3
	mulss	%xmm4, %xmm3
	addss	%xmm3, %xmm11

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, X1
	ALIGN_3

L(3X):
	testq	$2 * SIZE, A1
	je	L(3XX)

#ifdef movsd
	xorps	%xmm0, %xmm0
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(A1), %xmm0
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A1, LDA), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
#ifdef movsd
	xorps	%xmm2, %xmm2
#endif
	movsd	-32 * SIZE(A2), %xmm2
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
#ifdef movsd
	xorps	%xmm3, %xmm3
#endif
	movsd	-32 * SIZE(A2, LDA), %xmm3
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm11

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, X1
	ALIGN_3

L(3XX):
#endif

	MOVUPS_XL1 (-32 * SIZE, X1, %xmm4)
	MOVUPS_XL1 (-28 * SIZE, X1, %xmm5)

#if (GEMV_UNROLL == 4) && defined(PREFETCHW)
	PREFETCHW	4 * SIZE(Y1)
#endif

	movq	MM,  I
	sarq	$4,  I
	jle	L(35)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-32 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A1 (-32 * SIZE, A2, %xmm2)

	MOVUPS_A1 (-28 * SIZE, A1, %xmm12)
	MOVUPS_A2 (-28 * SIZE, A1, LDA, 1, %xmm13)
	MOVUPS_A1 (-28 * SIZE, A2, %xmm14)

	decq	I
	jle	L(33)
	ALIGN_4

L(32):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A1)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-24 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm2
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-24 * SIZE, A2, %xmm2)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A1, LDA)
#endif

	mulps	%xmm5, %xmm12
	addps	%xmm12, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm12)
	mulps	%xmm5, %xmm13
	addps	%xmm13, %xmm9
	MOVUPS_A2 (-20 * SIZE, A1, LDA, 1, %xmm13)
	mulps	%xmm5, %xmm14
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm14, %xmm10
	MOVUPS_A1 (-20 * SIZE, A2, %xmm14)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A2)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-16 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-16 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm2
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-16 * SIZE, A2, %xmm2)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(X1)
#endif

	mulps	%xmm5, %xmm12
	addps	%xmm12, %xmm8
	MOVUPS_A1 (-12 * SIZE, A1, %xmm12)
	mulps	%xmm5, %xmm13
	addps	%xmm13, %xmm9
	MOVUPS_A2 (-12 * SIZE, A1, LDA, 1, %xmm13)
	mulps	%xmm5, %xmm14
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm14, %xmm10
	MOVUPS_A1 (-12 * SIZE, A2, %xmm14)

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1

	decq	I
	jg	L(32)
	ALIGN_4

L(33):
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-24 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm2
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-24 * SIZE, A2, %xmm2)

	mulps	%xmm5, %xmm12
	addps	%xmm12, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm12)
	mulps	%xmm5, %xmm13
	addps	%xmm13, %xmm9
	MOVUPS_A2 (-20 * SIZE, A1, LDA, 1, %xmm13)
	mulps	%xmm5, %xmm14
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm14, %xmm10
	MOVUPS_A1 (-20 * SIZE, A2, %xmm14)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)

	mulps	%xmm5, %xmm12
	addps	%xmm12, %xmm8
	mulps	%xmm5, %xmm13
	addps	%xmm13, %xmm9
	mulps	%xmm5, %xmm14
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm14, %xmm10

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1
	ALIGN_4

L(35):
	testq	$8, MM
	jle	L(36)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A2 (-32 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A1 (-32 * SIZE, A2, %xmm2)
	mulps	%xmm4, %xmm2
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm2, %xmm10

	MOVUPS_A1 (-28 * SIZE, A1, %xmm12)
	mulps	%xmm5, %xmm12
	addps	%xmm12, %xmm8
	MOVUPS_A2 (-28 * SIZE, A1, LDA, 1, %xmm13)
	mulps	%xmm5, %xmm13
	addps	%xmm13, %xmm9
	MOVUPS_A1 (-28 * SIZE, A2, %xmm14)
	mulps	%xmm5, %xmm14
	addps	%xmm14, %xmm10

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1
	ALIGN_4

L(36):
	testq	$4, MM
	jle	L(37)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A2 (-32 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A1 (-32 * SIZE, A2, %xmm2)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, A2
	addq	$4 * SIZE, X1
	ALIGN_4

L(37):
	testq	$2, MM
	jle	L(38)

#ifdef movsd
	xorps	%xmm0, %xmm0
#endif
	movsd	-32 * SIZE(A1), %xmm0
#ifdef movsd
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A1, LDA), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
#ifdef movsd
	xorps	%xmm2, %xmm2
#endif
	movsd	-32 * SIZE(A2), %xmm2
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
#ifdef movsd
	xorps	%xmm3, %xmm3
#endif

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, A2
	addq	$2 * SIZE, X1
	ALIGN_4

L(38):
	testq	$1, MM
	jle	L(39)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A1, LDA), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	movss	-32 * SIZE(A2), %xmm2
	mulss	%xmm4, %xmm2
	addss	%xmm2, %xmm10
	ALIGN_4

L(39):
#ifdef HAVE_SSE3
	haddps	%xmm9, %xmm8
	haddps	%xmm11, %xmm10
	haddps	%xmm10, %xmm8

	pshufd	$0x1, %xmm8, %xmm9
	pshufd	$0x2, %xmm8, %xmm10
#else
	movaps	%xmm8, %xmm0
	unpcklps %xmm9, %xmm8
	unpckhps %xmm9, %xmm0

	movaps	%xmm10, %xmm1
	unpcklps %xmm11, %xmm10
	unpckhps %xmm11, %xmm1

	movaps	%xmm8, %xmm9
	unpcklps %xmm10, %xmm8
	unpckhps %xmm10, %xmm9

	movaps	%xmm0, %xmm10
	unpcklps %xmm1, %xmm0
	unpckhps %xmm1, %xmm10

	addps	%xmm9, %xmm8
	addps	%xmm0, %xmm10
	addps	%xmm10, %xmm8

	pshufd	$0x2, %xmm8, %xmm9
	pshufd	$0x1, %xmm8, %xmm10
#endif

	mulss	ALPHA, %xmm8
	mulss	ALPHA, %xmm9
	mulss	ALPHA, %xmm10

	addss	(Y), %xmm8
	addq	INCY, Y
	addss	(Y), %xmm9
	addq	INCY, Y
	addss	(Y), %xmm10
	addq	INCY, Y

	movss	%xmm8, (Y1)
	addq	INCY, Y1
	movss	%xmm9, (Y1)
	addq	INCY, Y1
	movss	%xmm10, (Y1)
	addq	INCY, Y1
	jmp	L(999)
	ALIGN_4

L(40):
	cmpq	$2, N
	jne	L(50)

	leaq	32 * SIZE(BUFFER), X1

	movq	A, A1
	leaq	(A1, LDA), A2
	leaq	(A1, LDA, 2), A

	xorps	%xmm8, %xmm8
	xorps	%xmm9, %xmm9

#ifdef ALIGNED_ACCESS
	cmpq	$3, M
	jle	L(47)

	testq	$SIZE, A1
	je	L(4X)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A2), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, X1
	ALIGN_3

L(4X):
	testq	$2 * SIZE, A1
	je	L(4XX)

#ifdef movsd
	xorps	%xmm0, %xmm0
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(A1), %xmm0
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A2), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, X1
	ALIGN_3

L(4XX):
#endif

	MOVUPS_XL1 (-32 * SIZE, X1, %xmm4)
	MOVUPS_XL1 (-28 * SIZE, X1, %xmm5)

	movq	MM,  I
	sarq	$4,  I
	jle	L(45)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A1 (-32 * SIZE, A2, %xmm1)
	MOVUPS_A1 (-28 * SIZE, A1, %xmm12)
	MOVUPS_A1 (-28 * SIZE, A2, %xmm13)

	decq	I
	jle	L(43)
	ALIGN_4

L(42):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A1)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm1
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm1, %xmm9
	MOVUPS_A1 (-24 * SIZE, A2, %xmm1)

	mulps	%xmm5, %xmm12
	addps	%xmm12, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm12)
	mulps	%xmm5, %xmm13
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm13, %xmm9
	MOVUPS_A1 (-20 * SIZE, A2, %xmm13)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A2)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-16 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm1
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm1, %xmm9
	MOVUPS_A1 (-16 * SIZE, A2, %xmm1)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(X1)
#endif

	mulps	%xmm5, %xmm12
	addps	%xmm12, %xmm8
	MOVUPS_A1 (-12 * SIZE, A1, %xmm12)
	mulps	%xmm5, %xmm13
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm13, %xmm9
	MOVUPS_A1 (-12 * SIZE, A2, %xmm13)

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1

	decq	I
	jg	L(42)
	ALIGN_4

L(43):
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm1
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm1, %xmm9
	MOVUPS_A1 (-24 * SIZE, A2, %xmm1)

	mulps	%xmm5, %xmm12
	addps	%xmm12, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm12)
	mulps	%xmm5, %xmm13
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm13, %xmm9
	MOVUPS_A1 (-20 * SIZE, A2, %xmm13)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	mulps	%xmm4, %xmm1
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm1, %xmm9

	mulps	%xmm5, %xmm12
	addps	%xmm12, %xmm8
	mulps	%xmm5, %xmm13
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm13, %xmm9

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1
	ALIGN_4

L(45):
	testq	$8, MM
	jle	L(46)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-32 * SIZE, A2, %xmm1)
	mulps	%xmm4, %xmm1
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm1, %xmm9

	MOVUPS_A1 (-28 * SIZE, A1, %xmm12)
	mulps	%xmm5, %xmm12
	addps	%xmm12, %xmm8
	MOVUPS_A1 (-28 * SIZE, A2, %xmm13)
	mulps	%xmm5, %xmm13
	addps	%xmm13, %xmm9

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1
	ALIGN_4

L(46):
	testq	$4, MM
	jle	L(47)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-32 * SIZE, A2, %xmm1)
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, A2
	addq	$4 * SIZE, X1
	ALIGN_4

L(47):
	testq	$2, MM
	jle	L(48)

#ifdef movsd
	xorps	%xmm0, %xmm0
#endif
	movsd	-32 * SIZE(A1), %xmm0
#ifdef movsd
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A2), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	shufps	$0xe, %xmm4, %xmm4

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, A2
	addq	$2 * SIZE, X1
	ALIGN_4

L(48):
	testq	$1, MM
	jle	L(49)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A2), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	ALIGN_4

L(49):
#ifdef HAVE_SSE3
	haddps	%xmm9, %xmm8
	haddps	%xmm8, %xmm8
#else
	movaps	%xmm8, %xmm10
	unpcklps %xmm9, %xmm8
	unpckhps %xmm9, %xmm10

	addps	%xmm10, %xmm8
	movhlps %xmm8, %xmm9
	addps	%xmm9, %xmm8
#endif

	pshufd	$0x1, %xmm8, %xmm9

	mulss	ALPHA, %xmm8
	mulss	ALPHA, %xmm9

	addss	(Y), %xmm8
	addq	INCY, Y
	addss	(Y), %xmm9
	addq	INCY, Y

	movss	%xmm8, (Y1)
	addq	INCY, Y1
	movss	%xmm9, (Y1)
	addq	INCY, Y1
	jmp	L(999)
	ALIGN_4

L(50):
	cmpq	$1, N
	jne	L(999)

	leaq	32 * SIZE(BUFFER), X1

	movq	A, A1

	xorps	%xmm8, %xmm8
	xorps	%xmm9, %xmm9

#ifdef ALIGNED_ACCESS
	cmpq	$3, M
	jle	L(57)

	testq	$SIZE, A1
	je	L(5X)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, X1
	ALIGN_3

L(5X):
	testq	$2 * SIZE, A1
	je	L(5XX)


#ifdef movsd
	xorps	%xmm0, %xmm0
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(A1), %xmm0
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	shufps	$0xe, %xmm4, %xmm4

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, X1
	ALIGN_3

L(5XX):
#endif

	MOVUPS_XL1 (-32 * SIZE, X1, %xmm4)
	MOVUPS_XL1 (-28 * SIZE, X1, %xmm5)

	movq	MM,  I
	sarq	$4,  I
	jle	L(55)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A1 (-28 * SIZE, A1, %xmm12)

	decq	I
	jle	L(53)
	ALIGN_4

L(52):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(A1)
#endif

	mulps	%xmm4, %xmm0
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)

	mulps	%xmm5, %xmm12
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-20 * SIZE, A1, %xmm12)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(X1)
#endif

	mulps	%xmm4, %xmm0
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-16 * SIZE, A1, %xmm0)

	mulps	%xmm5, %xmm12
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-12 * SIZE, A1, %xmm12)

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, X1

	decq	I
	jg	L(52)
	ALIGN_4

L(53):
	mulps	%xmm4, %xmm0
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)

	mulps	%xmm5, %xmm12
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-20 * SIZE, A1, %xmm12)

	mulps	%xmm4, %xmm0
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm0, %xmm8

	mulps	%xmm5, %xmm12
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm12, %xmm9

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, X1
	ALIGN_4

L(55):
	testq	$8, MM
	jle	L(56)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm0
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm0, %xmm8

	MOVUPS_A1 (-28 * SIZE, A1, %xmm12)
	mulps	%xmm5, %xmm12
	addps	%xmm12, %xmm9

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, X1
	ALIGN_4

L(56):
	testq	$4, MM
	jle	L(57)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, X1
	ALIGN_4

L(57):
	testq	$2, MM
	jle	L(58)

#ifdef movsd
	xorps	%xmm0, %xmm0
#endif
	movsd	-32 * SIZE(A1), %xmm0
#ifdef movsd
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	shufps	$0xe, %xmm4, %xmm4

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, X1
	ALIGN_4

L(58):
	testq	$1, MM
	jle	L(59)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	ALIGN_4

L(59):
	addps	%xmm9, %xmm8

#ifdef HAVE_SSE3
	haddps	%xmm8, %xmm8
	haddps	%xmm8, %xmm8
#else
	pshufd	$1, %xmm8, %xmm9
	pshufd	$2, %xmm8, %xmm10
	pshufd	$3, %xmm8, %xmm11

	addss	%xmm9, %xmm8
	addss	%xmm11, %xmm10
	addss	%xmm10, %xmm8
#endif

	mulss	ALPHA, %xmm8

	addss	(Y), %xmm8
	movss	%xmm8, (Y1)

#ifdef ALIGNED_ACCESS
	jmp	L(999)
	ALIGN_4

L(100):
	testq	$2 * SIZE - 1, LDA
	jne	L(200)

	cmpq	$4, N
	jl	L(110)
	ALIGN_3

L(101):
	subq	$4, N

	leaq	32 * SIZE(BUFFER), X1

	movq	A, A1
	leaq	(A1, LDA, 2), A2
	leaq	(A1, LDA, 4), A

	xorps	%xmm8, %xmm8
	xorps	%xmm9, %xmm9
	xorps	%xmm10, %xmm10
	xorps	%xmm11, %xmm11

	cmpq	$3, M
	jle	L(107)

	testq	$SIZE, A1
	je	L(10X)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A1, LDA), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	movss	-32 * SIZE(A2), %xmm2
	mulss	%xmm4, %xmm2
	addss	%xmm2, %xmm10
	movss	-32 * SIZE(A2, LDA), %xmm3
	mulss	%xmm4, %xmm3
	addss	%xmm3, %xmm11

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, X1
	ALIGN_3

L(10X):
	testq	$2 * SIZE, A1
	je	L(10XX)

#ifdef movsd
	xorps	%xmm0, %xmm0
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(A1), %xmm0
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A1, LDA), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
#ifdef movsd
	xorps	%xmm2, %xmm2
#endif
	movsd	-32 * SIZE(A2), %xmm2
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
#ifdef movsd
	xorps	%xmm3, %xmm3
#endif
	movsd	-32 * SIZE(A2, LDA), %xmm3
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm11

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, X1
	ALIGN_3

L(10XX):
	MOVUPS_A2  (-34 * SIZE, A1, LDA, 1, %xmm12)
	MOVUPS_A2  (-34 * SIZE, A2, LDA, 1, %xmm13)

	MOVUPS_XL1 (-32 * SIZE, X1, %xmm4)
	MOVUPS_XL1 (-28 * SIZE, X1, %xmm5)

#ifdef PREFETCHW
	PREFETCHW	4 * SIZE(Y1)
#endif

	movq	MM,  I
	sarq	$4,  I
	jle	L(105)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-30 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A1 (-32 * SIZE, A2, %xmm2)
	MOVUPS_A2 (-30 * SIZE, A2, LDA, 1, %xmm3)

	decq	I
	jle	L(103)
	ALIGN_4

L(102):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm1, %xmm12
	mulps	%xmm4, %xmm12
	addps	%xmm12, %xmm9
	MOVUPS_A2 (-26 * SIZE, A1, LDA, 1, %xmm12)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-28 * SIZE, A2, %xmm2)
	shufps	$0x4e, %xmm3, %xmm13
	mulps	%xmm4, %xmm13
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm13, %xmm11
	MOVUPS_A2 (-26 * SIZE, A2, LDA, 1, %xmm13)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1, LDA)
#endif

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm12, %xmm1
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-22 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm5, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-24 * SIZE, A2, %xmm2)
	shufps	$0x4e, %xmm13, %xmm3
	mulps	%xmm5, %xmm3
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm3, %xmm11
	MOVUPS_A2 (-22 * SIZE, A2, LDA, 1, %xmm3)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm1, %xmm12
	mulps	%xmm4, %xmm12
	addps	%xmm12, %xmm9
	MOVUPS_A2 (-18 * SIZE, A1, LDA, 1, %xmm12)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-20 * SIZE, A2, %xmm2)
	shufps	$0x4e, %xmm3, %xmm13
	mulps	%xmm4, %xmm13
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm13, %xmm11
	MOVUPS_A2 (-18 * SIZE, A2, LDA, 1, %xmm13)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2, LDA)
#endif

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-16 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm12, %xmm1
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-14 * SIZE, A1, LDA, 1, %xmm1)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(X1)
#endif

	mulps	%xmm5, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-16 * SIZE, A2, %xmm2)
	shufps	$0x4e, %xmm13, %xmm3
	mulps	%xmm5, %xmm3
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm3, %xmm11
	MOVUPS_A2 (-14 * SIZE, A2, LDA, 1, %xmm3)

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1

	decq	I
	jg	L(102)
	ALIGN_4

L(103):
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm1, %xmm12
	mulps	%xmm4, %xmm12
	addps	%xmm12, %xmm9
	MOVUPS_A2 (-26 * SIZE, A1, LDA, 1, %xmm12)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-28 * SIZE, A2, %xmm2)
	shufps	$0x4e, %xmm3, %xmm13
	mulps	%xmm4, %xmm13
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm13, %xmm11
	MOVUPS_A2 (-26 * SIZE, A2, LDA, 1, %xmm13)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm12, %xmm1
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-22 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm5, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-24 * SIZE, A2, %xmm2)
	shufps	$0x4e, %xmm13, %xmm3
	mulps	%xmm5, %xmm3
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm3, %xmm11
	MOVUPS_A2 (-22 * SIZE, A2, LDA, 1, %xmm3)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm1, %xmm12
	mulps	%xmm4, %xmm12
	addps	%xmm12, %xmm9
	MOVUPS_A2 (-18 * SIZE, A1, LDA, 1, %xmm12)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-20 * SIZE, A2, %xmm2)
	shufps	$0x4e, %xmm3, %xmm13
	mulps	%xmm4, %xmm13
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm13, %xmm11
	MOVUPS_A2 (-18 * SIZE, A2, LDA, 1, %xmm13)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	shufps	$0x4e, %xmm12, %xmm1
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm9
	mulps	%xmm5, %xmm2
	addps	%xmm2, %xmm10
	shufps	$0x4e, %xmm13, %xmm3
	mulps	%xmm5, %xmm3
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm3, %xmm11

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1
	ALIGN_4

L(105):
	testq	$8, MM
	jle	L(106)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-30 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A1 (-32 * SIZE, A2, %xmm2)
	MOVUPS_A2 (-30 * SIZE, A2, LDA, 1, %xmm3)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm1, %xmm12
	mulps	%xmm4, %xmm12
	addps	%xmm12, %xmm9
	MOVUPS_A2 (-26 * SIZE, A1, LDA, 1, %xmm12)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-28 * SIZE, A2, %xmm2)
	shufps	$0x4e, %xmm3, %xmm13
	mulps	%xmm4, %xmm13
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm13, %xmm11
	MOVUPS_A2 (-26 * SIZE, A2, LDA, 1, %xmm13)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	shufps	$0x4e, %xmm12, %xmm1
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm9
	mulps	%xmm5, %xmm2
	addps	%xmm2, %xmm10
	shufps	$0x4e, %xmm13, %xmm3
	mulps	%xmm5, %xmm3
	addps	%xmm3, %xmm11

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1
	ALIGN_4

L(106):
	testq	$4, MM
	jle	L(107)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A2 (-30 * SIZE, A1, LDA, 1, %xmm1)
	shufps	$0x4e, %xmm1, %xmm12
	mulps	%xmm4, %xmm12
	addps	%xmm12, %xmm9

	MOVUPS_A1 (-32 * SIZE, A2, %xmm2)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
	MOVUPS_A2 (-30 * SIZE, A2, LDA, 1, %xmm3)
	shufps	$0x4e, %xmm3, %xmm13
	mulps	%xmm4, %xmm13
	addps	%xmm13, %xmm11

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, A2
	addq	$4 * SIZE, X1
	ALIGN_4

L(107):
	testq	$2, MM
	jle	L(108)

#ifdef movsd
	xorps	%xmm0, %xmm0
#endif
	movsd	-32 * SIZE(A1), %xmm0
#ifdef movsd
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A1, LDA), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
#ifdef movsd
	xorps	%xmm2, %xmm2
#endif
	movsd	-32 * SIZE(A2), %xmm2
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
#ifdef movsd
	xorps	%xmm3, %xmm3
#endif
	movsd	-32 * SIZE(A2, LDA), %xmm3
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm11
	shufps	$0xe, %xmm4, %xmm4

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, A2
	addq	$2 * SIZE, X1
	ALIGN_4

L(108):
	testq	$1, MM
	jle	L(109)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A1, LDA), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	movss	-32 * SIZE(A2), %xmm2
	mulss	%xmm4, %xmm2
	addss	%xmm2, %xmm10
	movss	-32 * SIZE(A2, LDA), %xmm3
	mulss	%xmm4, %xmm3
	addss	%xmm3, %xmm11
	ALIGN_4

L(109):
#ifdef HAVE_SSE3
	haddps	%xmm9, %xmm8
	haddps	%xmm11, %xmm10
	haddps	%xmm10, %xmm8

	pshufd	$0x1, %xmm8, %xmm9
	pshufd	$0x2, %xmm8, %xmm10
	pshufd	$0x3, %xmm8, %xmm11
#else
	movaps	%xmm8, %xmm0
	unpcklps %xmm9, %xmm8
	unpckhps %xmm9, %xmm0

	movaps	%xmm10, %xmm1
	unpcklps %xmm11, %xmm10
	unpckhps %xmm11, %xmm1

	movaps	%xmm8, %xmm9
	unpcklps %xmm10, %xmm8
	unpckhps %xmm10, %xmm9

	movaps	%xmm0, %xmm10
	unpcklps %xmm1, %xmm0
	unpckhps %xmm1, %xmm10

	addps	%xmm9, %xmm8
	addps	%xmm0, %xmm10
	addps	%xmm10, %xmm8

	pshufd	$0x2, %xmm8, %xmm9
	pshufd	$0x1, %xmm8, %xmm10
	pshufd	$0x3, %xmm8, %xmm11
#endif

	mulss	ALPHA, %xmm8
	mulss	ALPHA, %xmm9
	mulss	ALPHA, %xmm10
	mulss	ALPHA, %xmm11

	addss	(Y), %xmm8
	addq	INCY, Y
	addss	(Y), %xmm9
	addq	INCY, Y
	addss	(Y), %xmm10
	addq	INCY, Y
	addss	(Y), %xmm11
	addq	INCY, Y

	movss	%xmm8, (Y1)
	addq	INCY, Y1
	movss	%xmm9, (Y1)
	addq	INCY, Y1
	movss	%xmm10, (Y1)
	addq	INCY, Y1
	movss	%xmm11, (Y1)
	addq	INCY, Y1

	cmpq	$4, N
	jge	L(101)
	ALIGN_4

L(110):
	cmpq	$3, N
	jne	L(120)

	leaq	32 * SIZE(BUFFER), X1

	movq	A, A1
	leaq	(A1, LDA, 2), A2
	leaq	(A1, LDA, 4), A

	xorps	%xmm8, %xmm8
	xorps	%xmm9, %xmm9
	xorps	%xmm10, %xmm10

	cmpq	$3, M
	jle	L(117)

	testq	$SIZE, A1
	je	L(11X)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A1, LDA), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	movss	-32 * SIZE(A2), %xmm2
	mulss	%xmm4, %xmm2
	addss	%xmm2, %xmm10

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, X1
	ALIGN_3

L(11X):
	testq	$2 * SIZE, A1
	je	L(11XX)

#ifdef movsd
	xorps	%xmm0, %xmm0
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(A1), %xmm0
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A1, LDA), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
#ifdef movsd
	xorps	%xmm2, %xmm2
#endif
	movsd	-32 * SIZE(A2), %xmm2
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, X1
	ALIGN_3

L(11XX):
	MOVUPS_A2  (-34 * SIZE, A1, LDA, 1, %xmm12)

	MOVUPS_XL1 (-32 * SIZE, X1, %xmm4)
	MOVUPS_XL1 (-28 * SIZE, X1, %xmm5)

	movq	MM,  I
	sarq	$4,  I
	jle	L(115)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-30 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A1 (-32 * SIZE, A2, %xmm2)

	decq	I
	jle	L(113)
	ALIGN_4

L(112):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A1)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm1, %xmm12
	mulps	%xmm4, %xmm12
	addps	%xmm12, %xmm9
	MOVUPS_A2 (-26 * SIZE, A1, LDA, 1, %xmm12)
	mulps	%xmm4, %xmm2
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-28 * SIZE, A2, %xmm2)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A1, LDA)
#endif

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm12, %xmm1
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-22 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm5, %xmm2
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-24 * SIZE, A2, %xmm2)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A2)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm1, %xmm12
	mulps	%xmm4, %xmm12
	addps	%xmm12, %xmm9
	MOVUPS_A2 (-18 * SIZE, A1, LDA, 1, %xmm12)
	mulps	%xmm4, %xmm2
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-20 * SIZE, A2, %xmm2)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(X1)
#endif

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-16 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm12, %xmm1
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-14 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm5, %xmm2
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-16 * SIZE, A2, %xmm2)

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1

	decq	I
	jg	L(112)
	ALIGN_4

L(113):
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm1, %xmm12
	mulps	%xmm4, %xmm12
	addps	%xmm12, %xmm9
	MOVUPS_A2 (-26 * SIZE, A1, LDA, 1, %xmm12)
	mulps	%xmm4, %xmm2
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-28 * SIZE, A2, %xmm2)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm12, %xmm1
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm9
	MOVUPS_A2 (-22 * SIZE, A1, LDA, 1, %xmm1)
	mulps	%xmm5, %xmm2
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-24 * SIZE, A2, %xmm2)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm1, %xmm12
	mulps	%xmm4, %xmm12
	addps	%xmm12, %xmm9
	MOVUPS_A2 (-18 * SIZE, A1, LDA, 1, %xmm12)
	mulps	%xmm4, %xmm2
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-20 * SIZE, A2, %xmm2)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	shufps	$0x4e, %xmm12, %xmm1
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm9
	mulps	%xmm5, %xmm2
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm2, %xmm10

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1
	ALIGN_4

L(115):
	testq	$8, MM
	jle	L(116)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-30 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A1 (-32 * SIZE, A2, %xmm2)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm1, %xmm12
	mulps	%xmm4, %xmm12
	addps	%xmm12, %xmm9
	MOVUPS_A2 (-26 * SIZE, A1, LDA, 1, %xmm12)
	mulps	%xmm4, %xmm2
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm2, %xmm10
	MOVUPS_A1 (-28 * SIZE, A2, %xmm2)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	shufps	$0x4e, %xmm12, %xmm1
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm9
	mulps	%xmm5, %xmm2
	addps	%xmm2, %xmm10

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1
	ALIGN_4

L(116):
	testq	$4, MM
	jle	L(117)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A2 (-30 * SIZE, A1, LDA, 1, %xmm1)
	shufps	$0x4e, %xmm1, %xmm12
	mulps	%xmm4, %xmm12
	addps	%xmm12, %xmm9

	MOVUPS_A1 (-32 * SIZE, A2, %xmm2)
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, A2
	addq	$4 * SIZE, X1
	ALIGN_4

L(117):
	testq	$2, MM
	jle	L(118)

#ifdef movsd
	xorps	%xmm0, %xmm0
#endif
	movsd	-32 * SIZE(A1), %xmm0
#ifdef movsd
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A1, LDA), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
#ifdef movsd
	xorps	%xmm2, %xmm2
#endif
	movsd	-32 * SIZE(A2), %xmm2
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, A2
	addq	$2 * SIZE, X1
	ALIGN_4

L(118):
	testq	$1, MM
	jle	L(119)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A1, LDA), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	movss	-32 * SIZE(A2), %xmm2
	mulss	%xmm4, %xmm2
	addss	%xmm2, %xmm10
	ALIGN_4

L(119):
#ifdef HAVE_SSE3
	haddps	%xmm9, %xmm8
	haddps	%xmm11, %xmm10
	haddps	%xmm10, %xmm8

	pshufd	$0x1, %xmm8, %xmm9
	pshufd	$0x2, %xmm8, %xmm10
#else
	movaps	%xmm8, %xmm0
	unpcklps %xmm9, %xmm8
	unpckhps %xmm9, %xmm0

	movaps	%xmm10, %xmm1
	unpcklps %xmm11, %xmm10
	unpckhps %xmm11, %xmm1

	movaps	%xmm8, %xmm9
	unpcklps %xmm10, %xmm8
	unpckhps %xmm10, %xmm9

	movaps	%xmm0, %xmm10
	unpcklps %xmm1, %xmm0
	unpckhps %xmm1, %xmm10

	addps	%xmm9, %xmm8
	addps	%xmm0, %xmm10
	addps	%xmm10, %xmm8

	pshufd	$0x2, %xmm8, %xmm9
	pshufd	$0x1, %xmm8, %xmm10
#endif

	mulss	ALPHA, %xmm8
	mulss	ALPHA, %xmm9
	mulss	ALPHA, %xmm10

	addss	(Y), %xmm8
	addq	INCY, Y
	addss	(Y), %xmm9
	addq	INCY, Y
	addss	(Y), %xmm10

	movss	%xmm8, (Y1)
	addq	INCY, Y1
	movss	%xmm9, (Y1)
	addq	INCY, Y1
	movss	%xmm10, (Y1)
	jmp	L(999)
	ALIGN_4

L(120):
	cmpq	$2, N
	jne	L(130)

	leaq	32 * SIZE(BUFFER), X1

	movq	A, A1
	leaq	(A1, LDA), A2
	leaq	(A1, LDA, 2), A

	xorps	%xmm8, %xmm8
	xorps	%xmm9, %xmm9

	cmpq	$3, M
	jle	L(127)

	testq	$SIZE, A1
	je	L(12X)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A2), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, X1
	ALIGN_3

L(12X):
	testq	$2 * SIZE, A1
	je	L(12XX)

#ifdef movsd
	xorps	%xmm0, %xmm0
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(A1), %xmm0
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A2), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, X1
	ALIGN_3

L(12XX):
	MOVUPS_A1  (-34 * SIZE, A2, %xmm12)

	MOVUPS_XL1 (-32 * SIZE, X1, %xmm4)
	MOVUPS_XL1 (-28 * SIZE, X1, %xmm5)

	movq	MM,  I
	sarq	$4,  I
	jle	L(125)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A1 (-30 * SIZE, A2, %xmm1)

	decq	I
	jle	L(123)
	ALIGN_4

L(122):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A1)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm1, %xmm12
	mulps	%xmm4, %xmm12
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-26 * SIZE, A2, %xmm12)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm12, %xmm1
	mulps	%xmm5, %xmm1
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm1, %xmm9
	MOVUPS_A1 (-22 * SIZE, A2, %xmm1)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A2)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm1, %xmm12
	mulps	%xmm4, %xmm12
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-18 * SIZE, A2, %xmm12)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(X1)
#endif

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-16 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm12, %xmm1
	mulps	%xmm5, %xmm1
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm1, %xmm9
	MOVUPS_A1 (-14 * SIZE, A2, %xmm1)

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1

	decq	I
	jg	L(122)
	ALIGN_4

L(123):
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm1, %xmm12
	mulps	%xmm4, %xmm12
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-26 * SIZE, A2, %xmm12)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm12, %xmm1
	mulps	%xmm5, %xmm1
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm1, %xmm9
	MOVUPS_A1 (-22 * SIZE, A2, %xmm1)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm1, %xmm12
	mulps	%xmm4, %xmm12
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-18 * SIZE, A2, %xmm12)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	shufps	$0x4e, %xmm12, %xmm1
	mulps	%xmm5, %xmm1
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm1, %xmm9

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1
	ALIGN_4

L(125):
	testq	$8, MM
	jle	L(126)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A1 (-30 * SIZE, A2, %xmm1)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	shufps	$0x4e, %xmm1, %xmm12
	mulps	%xmm4, %xmm12
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-26 * SIZE, A2, %xmm12)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	shufps	$0x4e, %xmm12, %xmm1
	mulps	%xmm5, %xmm1
	addps	%xmm1, %xmm9

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1
	ALIGN_4

L(126):
	testq	$4, MM
	jle	L(127)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-30 * SIZE, A2, %xmm1)
	shufps	$0x4e, %xmm1, %xmm12
	mulps	%xmm4, %xmm12
	addps	%xmm12, %xmm9

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, A2
	addq	$4 * SIZE, X1
	ALIGN_4

L(127):
	testq	$2, MM
	jle	L(128)

#ifdef movsd
	xorps	%xmm0, %xmm0
#endif
	movsd	-32 * SIZE(A1), %xmm0
#ifdef movsd
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A2), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	shufps	$0xe, %xmm4, %xmm4

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, A2
	addq	$2 * SIZE, X1
	ALIGN_4

L(128):
	testq	$1, MM
	jle	L(129)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A2), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	ALIGN_4

L(129):
#ifdef HAVE_SSE3
	haddps	%xmm9, %xmm8
	haddps	%xmm8, %xmm8
#else
	movaps	%xmm8, %xmm10
	unpcklps %xmm9, %xmm8
	unpckhps %xmm9, %xmm10

	addps	%xmm10, %xmm8
	movhlps %xmm8, %xmm9
	addps	%xmm9, %xmm8
#endif

	pshufd	$0x1, %xmm8, %xmm9

	mulss	ALPHA, %xmm8
	mulss	ALPHA, %xmm9

	addss	(Y), %xmm8
	addq	INCY, Y
	addss	(Y), %xmm9
	addq	INCY, Y

	movss	%xmm8, (Y1)
	addq	INCY, Y1
	movss	%xmm9, (Y1)
	addq	INCY, Y1
	jmp	L(999)
	ALIGN_4

L(130):
	cmpq	$1, N
	jne	L(999)

	leaq	32 * SIZE(BUFFER), X1

	movq	A, A1

	xorps	%xmm8, %xmm8
	xorps	%xmm9, %xmm9

	cmpq	$3, M
	jle	L(137)

	testq	$SIZE, A1
	je	L(13X)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, X1
	ALIGN_3

L(13X):
	testq	$2 * SIZE, A1
	je	L(13XX)


#ifdef movsd
	xorps	%xmm0, %xmm0
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(A1), %xmm0
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	shufps	$0xe, %xmm4, %xmm4

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, X1
	ALIGN_3

L(13XX):

	MOVUPS_XL1 (-32 * SIZE, X1, %xmm4)
	MOVUPS_XL1 (-28 * SIZE, X1, %xmm5)

	movq	MM,  I
	sarq	$4,  I
	jle	L(135)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A1 (-28 * SIZE, A1, %xmm12)

	decq	I
	jle	L(133)
	ALIGN_4

L(132):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(A1)
#endif

	mulps	%xmm4, %xmm0
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)

	mulps	%xmm5, %xmm12
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-20 * SIZE, A1, %xmm12)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(X1)
#endif

	mulps	%xmm4, %xmm0
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-16 * SIZE, A1, %xmm0)

	mulps	%xmm5, %xmm12
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-12 * SIZE, A1, %xmm12)

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, X1

	decq	I
	jg	L(132)
	ALIGN_4

L(133):
	mulps	%xmm4, %xmm0
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)

	mulps	%xmm5, %xmm12
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-20 * SIZE, A1, %xmm12)

	mulps	%xmm4, %xmm0
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm0, %xmm8

	mulps	%xmm5, %xmm12
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm12, %xmm9

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, X1
	ALIGN_4

L(135):
	testq	$8, MM
	jle	L(136)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm0
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm0, %xmm8

	MOVUPS_A1 (-28 * SIZE, A1, %xmm12)
	mulps	%xmm5, %xmm12
	addps	%xmm12, %xmm9

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, X1
	ALIGN_4

L(136):
	testq	$4, MM
	jle	L(137)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, X1
	ALIGN_4

L(137):
	testq	$2, MM
	jle	L(138)

#ifdef movsd
	xorps	%xmm0, %xmm0
#endif
	movsd	-32 * SIZE(A1), %xmm0
#ifdef movsd
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	shufps	$0xe, %xmm4, %xmm4

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, X1
	ALIGN_4

L(138):
	testq	$1, MM
	jle	L(139)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	ALIGN_4

L(139):
	addps	%xmm9, %xmm8

#ifdef HAVE_SSE3
	haddps	%xmm8, %xmm8
	haddps	%xmm8, %xmm8
#else
	pshufd	$1, %xmm8, %xmm9
	pshufd	$2, %xmm8, %xmm10
	pshufd	$3, %xmm8, %xmm11

	addss	%xmm9, %xmm8
	addss	%xmm11, %xmm10
	addss	%xmm10, %xmm8
#endif

	mulss	ALPHA, %xmm8

	addss	(Y), %xmm8
	movss	%xmm8, (Y1)
	jmp	L(999)
	ALIGN_4

L(200):
	testq	$2 * SIZE, LDA
	jne	L(300)

	cmpq	$4, N
	jl	L(210)
	ALIGN_3

L(201):
	subq	$4, N

	leaq	32 * SIZE(BUFFER), X1

	movq	A, A1
	leaq	(A1, LDA, 2), A2
	leaq	(A1, LDA, 4), A

	xorps	%xmm8, %xmm8
	xorps	%xmm9, %xmm9
	xorps	%xmm10, %xmm10
	xorps	%xmm11, %xmm11

	cmpq	$3, M
	jle	L(207)

	testq	$SIZE, A1
	je	L(20X)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A1, LDA), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	movss	-32 * SIZE(A2), %xmm2
	mulss	%xmm4, %xmm2
	addss	%xmm2, %xmm10
	movss	-32 * SIZE(A2, LDA), %xmm3
	mulss	%xmm4, %xmm3
	addss	%xmm3, %xmm11

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, X1
	ALIGN_3

L(20X):
	testq	$2 * SIZE, A1
	je	L(20XX)

#ifdef movsd
	xorps	%xmm0, %xmm0
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(A1), %xmm0
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A1, LDA), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
#ifdef movsd
	xorps	%xmm2, %xmm2
#endif
	movsd	-32 * SIZE(A2), %xmm2
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
#ifdef movsd
	xorps	%xmm3, %xmm3
#endif
	movsd	-32 * SIZE(A2, LDA), %xmm3
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm11

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, X1
	ALIGN_3

L(20XX):
	movaps	-33 * SIZE(A1, LDA), %xmm12
	movaps	-34 * SIZE(A2), %xmm13
	movaps	-35 * SIZE(A2, LDA), %xmm14

	MOVUPS_XL1 (-32 * SIZE, X1, %xmm4)
	MOVUPS_XL1 (-28 * SIZE, X1, %xmm5)

#ifdef PREFETCHW
	PREFETCHW	4 * SIZE(Y1)
#endif

	movq	MM,  I
	sarq	$4,  I
	jle	L(205)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-29 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A1 (-30 * SIZE, A2, %xmm2)
	MOVUPS_A2 (-31 * SIZE, A2, LDA, 1, %xmm3)

	decq	I
	jle	L(203)
	ALIGN_4

L(202):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x39,  %xmm12, %xmm12
	mulps	 %xmm4,  %xmm12
	addps	 %xmm12, %xmm9
	MOVUPS_A2 (-25 * SIZE, A1, LDA, 1, %xmm12)

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	addps	 %xmm13, %xmm10
	MOVUPS_A1 (-26 * SIZE, A2, %xmm13)
	movss	 %xmm3,  %xmm14
	shufps	 $0x93, %xmm3, %xmm14
	mulps	 %xmm4, %xmm14
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	 %xmm14, %xmm11
	MOVUPS_A2 (-27 * SIZE, A2, LDA, 1, %xmm14)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1, LDA)
#endif

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	movss	 %xmm12,  %xmm1
	shufps	 $0x39,  %xmm1, %xmm1
	mulps	 %xmm5,  %xmm1
	addps	 %xmm1, %xmm9
	MOVUPS_A2 (-21 * SIZE, A1, LDA, 1, %xmm1)

	shufps	 $0x4e, %xmm13, %xmm2
	mulps	 %xmm5, %xmm2
	addps	 %xmm2, %xmm10
	MOVUPS_A1 (-22 * SIZE, A2, %xmm2)
	movss	 %xmm14,  %xmm3
	shufps	 $0x93, %xmm14, %xmm3
	mulps	 %xmm5, %xmm3
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	 %xmm3, %xmm11
	MOVUPS_A2 (-23 * SIZE, A2, LDA, 1, %xmm3)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x39,  %xmm12, %xmm12
	mulps	 %xmm4,  %xmm12
	addps	 %xmm12, %xmm9
	MOVUPS_A2 (-17 * SIZE, A1, LDA, 1, %xmm12)

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	addps	 %xmm13, %xmm10
	MOVUPS_A1 (-18 * SIZE, A2, %xmm13)
	movss	 %xmm3,  %xmm14
	shufps	 $0x93, %xmm3, %xmm14
	mulps	 %xmm4, %xmm14
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	 %xmm14, %xmm11
	MOVUPS_A2 (-19 * SIZE, A2, LDA, 1, %xmm14)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2, LDA)
#endif

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-16 * SIZE, A1, %xmm0)
	movss	 %xmm12,  %xmm1
	shufps	 $0x39,  %xmm1, %xmm1
	mulps	 %xmm5,  %xmm1
	addps	 %xmm1, %xmm9
	MOVUPS_A2 (-13 * SIZE, A1, LDA, 1, %xmm1)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(X1)
#endif

	shufps	 $0x4e, %xmm13, %xmm2
	mulps	 %xmm5, %xmm2
	addps	 %xmm2, %xmm10
	MOVUPS_A1 (-14 * SIZE, A2, %xmm2)
	movss	 %xmm14,  %xmm3
	shufps	 $0x93, %xmm14, %xmm3
	mulps	 %xmm5, %xmm3
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	 %xmm3, %xmm11
	MOVUPS_A2 (-15 * SIZE, A2, LDA, 1, %xmm3)

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1

	decq	I
	jg	L(202)
	ALIGN_4

L(203):
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x39,  %xmm12, %xmm12
	mulps	 %xmm4,  %xmm12
	addps	 %xmm12, %xmm9
	MOVUPS_A2 (-25 * SIZE, A1, LDA, 1, %xmm12)

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	addps	 %xmm13, %xmm10
	MOVUPS_A1 (-26 * SIZE, A2, %xmm13)
	movss	 %xmm3,  %xmm14
	shufps	 $0x93, %xmm3, %xmm14
	mulps	 %xmm4, %xmm14
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	 %xmm14, %xmm11
	MOVUPS_A2 (-27 * SIZE, A2, LDA, 1, %xmm14)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	movss	 %xmm12,  %xmm1
	shufps	 $0x39,  %xmm1, %xmm1
	mulps	 %xmm5,  %xmm1
	addps	 %xmm1, %xmm9
	MOVUPS_A2 (-21 * SIZE, A1, LDA, 1, %xmm1)

	shufps	 $0x4e, %xmm13, %xmm2
	mulps	 %xmm5, %xmm2
	addps	 %xmm2, %xmm10
	MOVUPS_A1 (-22 * SIZE, A2, %xmm2)
	movss	 %xmm14,  %xmm3
	shufps	 $0x93, %xmm14, %xmm3
	mulps	 %xmm5, %xmm3
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	 %xmm3, %xmm11
	MOVUPS_A2 (-23 * SIZE, A2, LDA, 1, %xmm3)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x39,  %xmm12, %xmm12
	mulps	 %xmm4,  %xmm12
	addps	 %xmm12, %xmm9
	MOVUPS_A2 (-17 * SIZE, A1, LDA, 1, %xmm12)

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	addps	 %xmm13, %xmm10
	MOVUPS_A1 (-18 * SIZE, A2, %xmm13)
	movss	 %xmm3,  %xmm14
	shufps	 $0x93, %xmm3, %xmm14
	mulps	 %xmm4, %xmm14
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	 %xmm14, %xmm11
	MOVUPS_A2 (-19 * SIZE, A2, LDA, 1, %xmm14)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	movss	 %xmm12,  %xmm1
	shufps	 $0x39,  %xmm1, %xmm1
	mulps	 %xmm5,  %xmm1
	addps	 %xmm1, %xmm9

	shufps	 $0x4e, %xmm13, %xmm2
	mulps	 %xmm5, %xmm2
	addps	 %xmm2, %xmm10
	movss	 %xmm14,  %xmm3
	shufps	 $0x93, %xmm14, %xmm3
	mulps	 %xmm5, %xmm3
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	 %xmm3, %xmm11

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1
	ALIGN_4

L(205):
	testq	$8, MM
	jle	L(206)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-29 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A1 (-30 * SIZE, A2, %xmm2)
	MOVUPS_A2 (-31 * SIZE, A2, LDA, 1, %xmm3)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x39,  %xmm12, %xmm12
	mulps	 %xmm4,  %xmm12
	addps	 %xmm12, %xmm9
	MOVUPS_A2 (-25 * SIZE, A1, LDA, 1, %xmm12)

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	addps	 %xmm13, %xmm10
	MOVUPS_A1 (-26 * SIZE, A2, %xmm13)
	movss	 %xmm3,  %xmm14
	shufps	 $0x93, %xmm3, %xmm14
	mulps	 %xmm4, %xmm14
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	 %xmm14, %xmm11
	MOVUPS_A2 (-27 * SIZE, A2, LDA, 1, %xmm14)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	movss	 %xmm12,  %xmm1
	shufps	 $0x39,  %xmm1, %xmm1
	mulps	 %xmm5,  %xmm1
	addps	 %xmm1, %xmm9

	shufps	 $0x4e, %xmm13, %xmm2
	mulps	 %xmm5, %xmm2
	addps	 %xmm2, %xmm10
	movss	 %xmm14,  %xmm3
	shufps	 $0x93, %xmm14, %xmm3
	mulps	 %xmm5, %xmm3
	addps	 %xmm3, %xmm11

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1
	ALIGN_4

L(206):
	testq	$4, MM
	jle	L(207)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-29 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A1 (-30 * SIZE, A2, %xmm2)
	MOVUPS_A2 (-31 * SIZE, A2, LDA, 1, %xmm3)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	movss	 %xmm1,  %xmm12
	shufps	 $0x39,  %xmm12, %xmm12
	mulps	 %xmm4,  %xmm12
	addps	 %xmm12, %xmm9

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	addps	 %xmm13, %xmm10
	movss	 %xmm3,  %xmm14
	shufps	 $0x93, %xmm3, %xmm14
	mulps	 %xmm4, %xmm14
	addps	 %xmm14, %xmm11

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, A2
	addq	$4 * SIZE, X1
	ALIGN_4

L(207):
	testq	$2, MM
	jle	L(208)

#ifdef movsd
	xorps	%xmm0, %xmm0
#endif
	movsd	-32 * SIZE(A1), %xmm0
#ifdef movsd
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A1, LDA), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
#ifdef movsd
	xorps	%xmm2, %xmm2
#endif
	movsd	-32 * SIZE(A2), %xmm2
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
#ifdef movsd
	xorps	%xmm3, %xmm3
#endif
	movsd	-32 * SIZE(A2, LDA), %xmm3
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm11
	shufps	$0xe, %xmm4, %xmm4

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, A2
	addq	$2 * SIZE, X1
	ALIGN_4

L(208):
	testq	$1, MM
	jle	L(209)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A1, LDA), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	movss	-32 * SIZE(A2), %xmm2
	mulss	%xmm4, %xmm2
	addss	%xmm2, %xmm10
	movss	-32 * SIZE(A2, LDA), %xmm3
	mulss	%xmm4, %xmm3
	addss	%xmm3, %xmm11
	ALIGN_4

L(209):
#ifdef HAVE_SSE3
	haddps	%xmm9, %xmm8
	haddps	%xmm11, %xmm10
	haddps	%xmm10, %xmm8

	pshufd	$0x1, %xmm8, %xmm9
	pshufd	$0x2, %xmm8, %xmm10
	pshufd	$0x3, %xmm8, %xmm11
#else
	movaps	%xmm8, %xmm0
	unpcklps %xmm9, %xmm8
	unpckhps %xmm9, %xmm0

	movaps	%xmm10, %xmm1
	unpcklps %xmm11, %xmm10
	unpckhps %xmm11, %xmm1

	movaps	%xmm8, %xmm9
	unpcklps %xmm10, %xmm8
	unpckhps %xmm10, %xmm9

	movaps	%xmm0, %xmm10
	unpcklps %xmm1, %xmm0
	unpckhps %xmm1, %xmm10

	addps	%xmm9, %xmm8
	addps	%xmm0, %xmm10
	addps	%xmm10, %xmm8

	pshufd	$0x2, %xmm8, %xmm9
	pshufd	$0x1, %xmm8, %xmm10
	pshufd	$0x3, %xmm8, %xmm11
#endif

	mulss	ALPHA, %xmm8
	mulss	ALPHA, %xmm9
	mulss	ALPHA, %xmm10
	mulss	ALPHA, %xmm11

	addss	(Y), %xmm8
	addq	INCY, Y
	addss	(Y), %xmm9
	addq	INCY, Y
	addss	(Y), %xmm10
	addq	INCY, Y
	addss	(Y), %xmm11
	addq	INCY, Y

	movss	%xmm8, (Y1)
	addq	INCY, Y1
	movss	%xmm9, (Y1)
	addq	INCY, Y1
	movss	%xmm10, (Y1)
	addq	INCY, Y1
	movss	%xmm11, (Y1)
	addq	INCY, Y1

	cmpq	$4, N
	jge	L(201)
	ALIGN_4

L(210):
	cmpq	$3, N
	jne	L(220)

	leaq	32 * SIZE(BUFFER), X1

	movq	A, A1
	leaq	(A1, LDA, 2), A2
	leaq	(A1, LDA, 4), A

	xorps	%xmm8, %xmm8
	xorps	%xmm9, %xmm9
	xorps	%xmm10, %xmm10

	cmpq	$3, M
	jle	L(217)

	testq	$SIZE, A1
	je	L(21X)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A1, LDA), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	movss	-32 * SIZE(A2), %xmm2
	mulss	%xmm4, %xmm2
	addss	%xmm2, %xmm10

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, X1
	ALIGN_3

L(21X):
	testq	$2 * SIZE, A1
	je	L(21XX)

#ifdef movsd
	xorps	%xmm0, %xmm0
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(A1), %xmm0
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A1, LDA), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
#ifdef movsd
	xorps	%xmm2, %xmm2
#endif
	movsd	-32 * SIZE(A2), %xmm2
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, X1
	ALIGN_3

L(21XX):
	movaps	-33 * SIZE(A1, LDA), %xmm12
	movaps	-34 * SIZE(A2), %xmm13

	MOVUPS_XL1 (-32 * SIZE, X1, %xmm4)
	MOVUPS_XL1 (-28 * SIZE, X1, %xmm5)

#ifdef PREFETCHW
	PREFETCHW	4 * SIZE(Y1)
#endif

	movq	MM,  I
	sarq	$4,  I
	jle	L(215)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-29 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A1 (-30 * SIZE, A2, %xmm2)

	decq	I
	jle	L(213)
	ALIGN_4

L(212):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A1)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x39,  %xmm12, %xmm12
	mulps	 %xmm4,  %xmm12
	addps	 %xmm12, %xmm9
	MOVUPS_A2 (-25 * SIZE, A1, LDA, 1, %xmm12)

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	addps	 %xmm13, %xmm10
	MOVUPS_A1 (-26 * SIZE, A2, %xmm13)
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A1, LDA)
#endif

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	movss	 %xmm12,  %xmm1
	shufps	 $0x39,  %xmm1, %xmm1
	mulps	 %xmm5,  %xmm1
	addps	 %xmm1, %xmm9
	MOVUPS_A2 (-21 * SIZE, A1, LDA, 1, %xmm1)

	shufps	 $0x4e, %xmm13, %xmm2
	mulps	 %xmm5, %xmm2
	addps	 %xmm2, %xmm10
	MOVUPS_A1 (-22 * SIZE, A2, %xmm2)
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A2)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x39,  %xmm12, %xmm12
	mulps	 %xmm4,  %xmm12
	addps	 %xmm12, %xmm9
	MOVUPS_A2 (-17 * SIZE, A1, LDA, 1, %xmm12)
	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	 %xmm13, %xmm10
	MOVUPS_A1 (-18 * SIZE, A2, %xmm13)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(X1)
#endif

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-16 * SIZE, A1, %xmm0)
	movss	 %xmm12,  %xmm1
	shufps	 $0x39,  %xmm1, %xmm1
	mulps	 %xmm5,  %xmm1
	addps	 %xmm1, %xmm9
	MOVUPS_A2 (-13 * SIZE, A1, LDA, 1, %xmm1)
	shufps	 $0x4e, %xmm13, %xmm2
	mulps	 %xmm5, %xmm2
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	 %xmm2, %xmm10
	MOVUPS_A1 (-14 * SIZE, A2, %xmm2)

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1

	decq	I
	jg	L(212)
	ALIGN_4

L(213):
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x39,  %xmm12, %xmm12
	mulps	 %xmm4,  %xmm12
	addps	 %xmm12, %xmm9
	MOVUPS_A2 (-25 * SIZE, A1, LDA, 1, %xmm12)

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	addps	 %xmm13, %xmm10
	MOVUPS_A1 (-26 * SIZE, A2, %xmm13)
	movss	 %xmm3,  %xmm14
	shufps	 $0x93, %xmm3, %xmm14
	mulps	 %xmm4, %xmm14
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	 %xmm14, %xmm11
	MOVUPS_A2 (-27 * SIZE, A2, LDA, 1, %xmm14)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	movss	 %xmm12,  %xmm1
	shufps	 $0x39,  %xmm1, %xmm1
	mulps	 %xmm5,  %xmm1
	addps	 %xmm1, %xmm9
	MOVUPS_A2 (-21 * SIZE, A1, LDA, 1, %xmm1)

	shufps	 $0x4e, %xmm13, %xmm2
	mulps	 %xmm5, %xmm2
	addps	 %xmm2, %xmm10
	MOVUPS_A1 (-22 * SIZE, A2, %xmm2)
	movss	 %xmm14,  %xmm3
	shufps	 $0x93, %xmm14, %xmm3
	mulps	 %xmm5, %xmm3
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	 %xmm3, %xmm11
	MOVUPS_A2 (-23 * SIZE, A2, LDA, 1, %xmm3)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x39,  %xmm12, %xmm12
	mulps	 %xmm4,  %xmm12
	addps	 %xmm12, %xmm9
	MOVUPS_A2 (-17 * SIZE, A1, LDA, 1, %xmm12)

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	addps	 %xmm13, %xmm10
	MOVUPS_A1 (-18 * SIZE, A2, %xmm13)
	movss	 %xmm3,  %xmm14
	shufps	 $0x93, %xmm3, %xmm14
	mulps	 %xmm4, %xmm14
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	 %xmm14, %xmm11
	MOVUPS_A2 (-19 * SIZE, A2, LDA, 1, %xmm14)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	movss	 %xmm12,  %xmm1
	shufps	 $0x39,  %xmm1, %xmm1
	mulps	 %xmm5,  %xmm1
	addps	 %xmm1, %xmm9

	shufps	 $0x4e, %xmm13, %xmm2
	mulps	 %xmm5, %xmm2
	addps	 %xmm2, %xmm10
	movss	 %xmm14,  %xmm3
	shufps	 $0x93, %xmm14, %xmm3
	mulps	 %xmm5, %xmm3
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	 %xmm3, %xmm11

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1
	ALIGN_4

L(215):
	testq	$8, MM
	jle	L(216)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-29 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A1 (-30 * SIZE, A2, %xmm2)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x39,  %xmm12, %xmm12
	mulps	 %xmm4,  %xmm12
	addps	 %xmm12, %xmm9
	MOVUPS_A2 (-25 * SIZE, A1, LDA, 1, %xmm12)

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	 %xmm13, %xmm10
	MOVUPS_A1 (-26 * SIZE, A2, %xmm13)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	movss	 %xmm12,  %xmm1
	shufps	 $0x39,  %xmm1, %xmm1
	mulps	 %xmm5,  %xmm1
	addps	 %xmm1, %xmm9

	shufps	 $0x4e, %xmm13, %xmm2
	mulps	 %xmm5, %xmm2
	addps	 %xmm2, %xmm10

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1
	ALIGN_4

L(216):
	testq	$4, MM
	jle	L(217)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-29 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A1 (-30 * SIZE, A2, %xmm2)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	movss	 %xmm1,  %xmm12
	shufps	 $0x39,  %xmm12, %xmm12
	mulps	 %xmm4,  %xmm12
	addps	 %xmm12, %xmm9

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	addps	 %xmm13, %xmm10

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, A2
	addq	$4 * SIZE, X1
	ALIGN_4

L(217):
	testq	$2, MM
	jle	L(218)

#ifdef movsd
	xorps	%xmm0, %xmm0
#endif
	movsd	-32 * SIZE(A1), %xmm0
#ifdef movsd
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A1, LDA), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
#ifdef movsd
	xorps	%xmm2, %xmm2
#endif
	movsd	-32 * SIZE(A2), %xmm2
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, A2
	addq	$2 * SIZE, X1
	ALIGN_4

L(218):
	testq	$1, MM
	jle	L(219)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A1, LDA), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	movss	-32 * SIZE(A2), %xmm2
	mulss	%xmm4, %xmm2
	addss	%xmm2, %xmm10
	ALIGN_4

L(219):
#ifdef HAVE_SSE3
	haddps	%xmm9, %xmm8
	haddps	%xmm11, %xmm10
	haddps	%xmm10, %xmm8

	pshufd	$0x1, %xmm8, %xmm9
	pshufd	$0x2, %xmm8, %xmm10
#else
	movaps	%xmm8, %xmm0
	unpcklps %xmm9, %xmm8
	unpckhps %xmm9, %xmm0

	movaps	%xmm10, %xmm1
	unpcklps %xmm11, %xmm10
	unpckhps %xmm11, %xmm1

	movaps	%xmm8, %xmm9
	unpcklps %xmm10, %xmm8
	unpckhps %xmm10, %xmm9

	movaps	%xmm0, %xmm10
	unpcklps %xmm1, %xmm0
	unpckhps %xmm1, %xmm10

	addps	%xmm9, %xmm8
	addps	%xmm0, %xmm10
	addps	%xmm10, %xmm8

	pshufd	$0x2, %xmm8, %xmm9
	pshufd	$0x1, %xmm8, %xmm10
#endif

	mulss	ALPHA, %xmm8
	mulss	ALPHA, %xmm9
	mulss	ALPHA, %xmm10

	addss	(Y), %xmm8
	addq	INCY, Y
	addss	(Y), %xmm9
	addq	INCY, Y
	addss	(Y), %xmm10

	movss	%xmm8, (Y1)
	addq	INCY, Y1
	movss	%xmm9, (Y1)
	addq	INCY, Y1
	movss	%xmm10, (Y1)
	jmp	L(999)
	ALIGN_4

L(220):
	testq	N, N
	jle	L(999)

	cmpq	$2, N
	jne	L(230)

	leaq	32 * SIZE(BUFFER), X1

	movq	A, A1
	leaq	(A1, LDA), A2
	leaq	(A1, LDA, 2), A

	xorps	%xmm8, %xmm8
	xorps	%xmm9, %xmm9

	cmpq	$3, M
	jle	L(227)

	testq	$SIZE, A1
	je	L(22X)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A2), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, X1
	ALIGN_3

L(22X):
	testq	$2 * SIZE, A1
	je	L(22XX)

#ifdef movsd
	xorps	%xmm0, %xmm0
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(A1), %xmm0
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A2), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, X1
	ALIGN_3

L(22XX):
	movaps	-33 * SIZE(A2), %xmm12

	MOVUPS_XL1 (-32 * SIZE, X1, %xmm4)
	MOVUPS_XL1 (-28 * SIZE, X1, %xmm5)

	movq	MM,  I
	sarq	$4,  I
	jle	L(225)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A1 (-29 * SIZE, A2, %xmm1)

	decq	I
	jle	L(223)
	ALIGN_4

L(222):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A1)
#endif

	mulps	%xmm4,  %xmm0
	addps	%xmm0,  %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm2)
	movss	%xmm1,  %xmm12
	shufps	$0x39,  %xmm12, %xmm12
	mulps	%xmm4,  %xmm12
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-25 * SIZE, A2, %xmm12)

	mulps	%xmm5,  %xmm2
	addps	%xmm2,  %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	movss	%xmm12, %xmm1
	shufps	$0x39,  %xmm1, %xmm1
	mulps	%xmm5,  %xmm1
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm1,  %xmm9
	MOVUPS_A1 (-21 * SIZE, A2, %xmm1)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A2)
#endif

	mulps	%xmm4,  %xmm0
	addps	%xmm0,  %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm2)
	movss	%xmm1,  %xmm12
	shufps	$0x39,  %xmm12, %xmm12
	mulps	%xmm4,  %xmm12
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-17 * SIZE, A2, %xmm12)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(X1)
#endif

	mulps	%xmm5,  %xmm2
	addps	%xmm2,  %xmm8
	MOVUPS_A1 (-16 * SIZE, A1, %xmm0)
	movss	%xmm12, %xmm1
	shufps	$0x39,  %xmm1, %xmm1
	mulps	%xmm5,  %xmm1
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm1,  %xmm9
	MOVUPS_A1 (-13 * SIZE, A2, %xmm1)

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1

	decq	I
	jg	L(222)
	ALIGN_4

L(223):
	mulps	%xmm4,  %xmm0
	addps	%xmm0,  %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm2)
	movss	%xmm1,  %xmm12
	shufps	$0x39,  %xmm12, %xmm12
	mulps	%xmm4,  %xmm12
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-25 * SIZE, A2, %xmm12)

	mulps	%xmm5,  %xmm2
	addps	%xmm2,  %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	movss	%xmm12, %xmm1
	shufps	$0x39,  %xmm1, %xmm1
	mulps	%xmm5,  %xmm1
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm1,  %xmm9
	MOVUPS_A1 (-21 * SIZE, A2, %xmm1)

	mulps	%xmm4,  %xmm0
	addps	%xmm0,  %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm2)
	movss	%xmm1,  %xmm12
	shufps	$0x39,  %xmm12, %xmm12
	mulps	%xmm4,  %xmm12
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-17 * SIZE, A2, %xmm12)

	mulps	%xmm5,  %xmm2
	addps	%xmm2,  %xmm8
	movss	%xmm12, %xmm1
	shufps	$0x39,  %xmm1, %xmm1
	mulps	%xmm5,  %xmm1
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm1,  %xmm9

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1
	ALIGN_4

L(225):
	testq	$8, MM
	jle	L(226)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A1 (-29 * SIZE, A2, %xmm1)

	mulps	%xmm4,  %xmm0
	addps	%xmm0,  %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm2)
	movss	%xmm1,  %xmm12
	shufps	$0x39,  %xmm12, %xmm12
	mulps	%xmm4,  %xmm12
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-25 * SIZE, A2, %xmm12)

	mulps	%xmm5,  %xmm2
	addps	%xmm2,  %xmm8
	movss	%xmm12, %xmm1
	shufps	$0x39,  %xmm1, %xmm1
	mulps	%xmm5,  %xmm1
	addps	%xmm1,  %xmm9

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1
	ALIGN_4

L(226):
	testq	$4, MM
	jle	L(227)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A1 (-29 * SIZE, A2, %xmm1)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8

	movss	 %xmm1,  %xmm12
	shufps	 $0x39,  %xmm12, %xmm12
	mulps	 %xmm4,  %xmm12
	addps	 %xmm12, %xmm9

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, A2
	addq	$4 * SIZE, X1
	ALIGN_4

L(227):
	testq	$2, MM
	jle	L(228)

#ifdef movsd
	xorps	%xmm0, %xmm0
#endif
	movsd	-32 * SIZE(A1), %xmm0
#ifdef movsd
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A2), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	shufps	$0xe, %xmm4, %xmm4

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, A2
	addq	$2 * SIZE, X1
	ALIGN_4

L(228):
	testq	$1, MM
	jle	L(229)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A2), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	ALIGN_4

L(229):
#ifdef HAVE_SSE3
	haddps	%xmm9, %xmm8
	haddps	%xmm8, %xmm8
#else
	movaps	%xmm8, %xmm10
	unpcklps %xmm9, %xmm8
	unpckhps %xmm9, %xmm10

	addps	%xmm10, %xmm8
	movhlps %xmm8, %xmm9
	addps	%xmm9, %xmm8
#endif

	pshufd	$0x1, %xmm8, %xmm9

	mulss	ALPHA, %xmm8
	mulss	ALPHA, %xmm9

	addss	(Y), %xmm8
	addq	INCY, Y
	addss	(Y), %xmm9
	addq	INCY, Y

	movss	%xmm8, (Y1)
	addq	INCY, Y1
	movss	%xmm9, (Y1)
	addq	INCY, Y1
	jmp	L(999)
	ALIGN_4

L(230):
	cmpq	$1, N
	jne	L(999)

	leaq	32 * SIZE(BUFFER), X1

	movq	A, A1

	xorps	%xmm8, %xmm8
	xorps	%xmm9, %xmm9

	cmpq	$3, M
	jle	L(237)

	testq	$SIZE, A1
	je	L(23X)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, X1
	ALIGN_3

L(23X):
	testq	$2 * SIZE, A1
	je	L(23XX)

#ifdef movsd
	xorps	%xmm0, %xmm0
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(A1), %xmm0
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	shufps	$0xe, %xmm4, %xmm4

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, X1
	ALIGN_3

L(23XX):
	MOVUPS_XL1 (-32 * SIZE, X1, %xmm4)
	MOVUPS_XL1 (-28 * SIZE, X1, %xmm5)


	movq	MM,  I
	sarq	$4,  I
	jle	L(235)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A1 (-28 * SIZE, A1, %xmm12)

	decq	I
	jle	L(233)
	ALIGN_4

L(232):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(A1)
#endif

	mulps	%xmm4, %xmm0
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)

	mulps	%xmm5, %xmm12
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-20 * SIZE, A1, %xmm12)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(X1)
#endif

	mulps	%xmm4, %xmm0
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-16 * SIZE, A1, %xmm0)

	mulps	%xmm5, %xmm12
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-12 * SIZE, A1, %xmm12)

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, X1

	decq	I
	jg	L(232)
	ALIGN_4

L(233):
	mulps	%xmm4, %xmm0
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)

	mulps	%xmm5, %xmm12
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-20 * SIZE, A1, %xmm12)

	mulps	%xmm4, %xmm0
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm0, %xmm8

	mulps	%xmm5, %xmm12
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm12, %xmm9

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, X1
	ALIGN_4

L(235):
	testq	$8, MM
	jle	L(236)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm0
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm0, %xmm8

	MOVUPS_A1 (-28 * SIZE, A1, %xmm12)
	mulps	%xmm5, %xmm12
	addps	%xmm12, %xmm9

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, X1
	ALIGN_4

L(236):
	testq	$4, MM
	jle	L(237)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, X1
	ALIGN_4

L(237):
	testq	$2, MM
	jle	L(238)

#ifdef movsd
	xorps	%xmm0, %xmm0
#endif
	movsd	-32 * SIZE(A1), %xmm0
#ifdef movsd
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	shufps	$0xe, %xmm4, %xmm4

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, X1
	ALIGN_4

L(238):
	testq	$1, MM
	jle	L(239)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	ALIGN_4

L(239):
	addps	%xmm9, %xmm8

#ifdef HAVE_SSE3
	haddps	%xmm8, %xmm8
	haddps	%xmm8, %xmm8
#else
	pshufd	$1, %xmm8, %xmm9
	pshufd	$2, %xmm8, %xmm10
	pshufd	$3, %xmm8, %xmm11

	addss	%xmm9, %xmm8
	addss	%xmm11, %xmm10
	addss	%xmm10, %xmm8
#endif

	mulss	ALPHA, %xmm8

	addss	(Y), %xmm8
	movss	%xmm8, (Y1)
	jmp	L(999)
	ALIGN_4

L(300):
	cmpq	$4, N
	jl	L(310)
	ALIGN_3

L(301):
	subq	$4, N

	leaq	32 * SIZE(BUFFER), X1

	movq	A, A1
	leaq	(A1, LDA, 2), A2
	leaq	(A1, LDA, 4), A

	xorps	%xmm8, %xmm8
	xorps	%xmm9, %xmm9
	xorps	%xmm10, %xmm10
	xorps	%xmm11, %xmm11

	cmpq	$3, M
	jle	L(307)

	testq	$SIZE, A1
	je	L(30X)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A1, LDA), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	movss	-32 * SIZE(A2), %xmm2
	mulss	%xmm4, %xmm2
	addss	%xmm2, %xmm10
	movss	-32 * SIZE(A2, LDA), %xmm3
	mulss	%xmm4, %xmm3
	addss	%xmm3, %xmm11

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, X1
	ALIGN_3

L(30X):
	testq	$2 * SIZE, A1
	je	L(30XX)

#ifdef movsd
	xorps	%xmm0, %xmm0
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(A1), %xmm0
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A1, LDA), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
#ifdef movsd
	xorps	%xmm2, %xmm2
#endif
	movsd	-32 * SIZE(A2), %xmm2
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
#ifdef movsd
	xorps	%xmm3, %xmm3
#endif
	movsd	-32 * SIZE(A2, LDA), %xmm3
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm11

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, X1
	ALIGN_3

L(30XX):
	movaps	-35 * SIZE(A1, LDA), %xmm12
	movaps	-34 * SIZE(A2), %xmm13
	movaps	-33 * SIZE(A2, LDA), %xmm14

	MOVUPS_XL1 (-32 * SIZE, X1, %xmm4)
	MOVUPS_XL1 (-28 * SIZE, X1, %xmm5)

#ifdef PREFETCHW
	PREFETCHW	4 * SIZE(Y1)
#endif

	movq	MM,  I
	sarq	$4,  I
	jle	L(305)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-31 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A1 (-30 * SIZE, A2, %xmm2)
	MOVUPS_A2 (-29 * SIZE, A2, LDA, 1, %xmm3)

	decq	I
	jle	L(303)
	ALIGN_4

L(302):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1)
#endif

	mulps	 %xmm4, %xmm0
	addps	 %xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x93, %xmm1, %xmm12
	mulps	 %xmm4, %xmm12
	addps	 %xmm12, %xmm9
	MOVUPS_A2 (-27 * SIZE, A1, LDA, 1, %xmm12)

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	addps	 %xmm13, %xmm10
	MOVUPS_A1 (-26 * SIZE, A2, %xmm13)
	movss	 %xmm3,  %xmm14
	shufps	 $0x39,  %xmm14, %xmm14
	mulps	 %xmm4,  %xmm14
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	 %xmm14, %xmm11
	MOVUPS_A2 (-25 * SIZE, A2, LDA, 1, %xmm14)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1, LDA)
#endif

	mulps	 %xmm5, %xmm0
	addps	 %xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	movss	 %xmm12,  %xmm1
	shufps	 $0x93, %xmm12, %xmm1
	mulps	 %xmm5, %xmm1
	addps	 %xmm1, %xmm9
	MOVUPS_A2 (-23 * SIZE, A1, LDA, 1, %xmm1)

	shufps	 $0x4e, %xmm13, %xmm2
	mulps	 %xmm5, %xmm2
	addps	 %xmm2, %xmm10
	MOVUPS_A1 (-22 * SIZE, A2, %xmm2)
	movss	 %xmm14,  %xmm3
	shufps	 $0x39,  %xmm3, %xmm3
	mulps	 %xmm5,  %xmm3
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	 %xmm3, %xmm11
	MOVUPS_A2 (-21 * SIZE, A2, LDA, 1, %xmm3)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2)
#endif

	mulps	 %xmm4, %xmm0
	addps	 %xmm0, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x93, %xmm1, %xmm12
	mulps	 %xmm4, %xmm12
	addps	 %xmm12, %xmm9
	MOVUPS_A2 (-19 * SIZE, A1, LDA, 1, %xmm12)

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	addps	 %xmm13, %xmm10
	MOVUPS_A1 (-18 * SIZE, A2, %xmm13)
	movss	 %xmm3,  %xmm14
	shufps	 $0x39,  %xmm14, %xmm14
	mulps	 %xmm4,  %xmm14
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	 %xmm14, %xmm11
	MOVUPS_A2 (-17 * SIZE, A2, LDA, 1, %xmm14)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2, LDA)
#endif

	mulps	 %xmm5, %xmm0
	addps	 %xmm0, %xmm8
	MOVUPS_A1 (-16 * SIZE, A1, %xmm0)
	movss	 %xmm12,  %xmm1
	shufps	 $0x93, %xmm12, %xmm1
	mulps	 %xmm5, %xmm1
	addps	 %xmm1, %xmm9
	MOVUPS_A2 (-15 * SIZE, A1, LDA, 1, %xmm1)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(X1)
#endif

	shufps	 $0x4e, %xmm13, %xmm2
	mulps	 %xmm5, %xmm2
	addps	 %xmm2, %xmm10
	MOVUPS_A1 (-14 * SIZE, A2, %xmm2)
	movss	 %xmm14,  %xmm3
	shufps	 $0x39,  %xmm3, %xmm3
	mulps	 %xmm5,  %xmm3
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	 %xmm3, %xmm11
	MOVUPS_A2 (-13 * SIZE, A2, LDA, 1, %xmm3)

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1

	decq	I
	jg	L(302)
	ALIGN_4

L(303):
	mulps	 %xmm4, %xmm0
	addps	 %xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x93, %xmm1, %xmm12
	mulps	 %xmm4, %xmm12
	addps	 %xmm12, %xmm9
	MOVUPS_A2 (-27 * SIZE, A1, LDA, 1, %xmm12)

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	addps	 %xmm13, %xmm10
	MOVUPS_A1 (-26 * SIZE, A2, %xmm13)
	movss	 %xmm3,  %xmm14
	shufps	 $0x39,  %xmm14, %xmm14
	mulps	 %xmm4,  %xmm14
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	 %xmm14, %xmm11
	MOVUPS_A2 (-25 * SIZE, A2, LDA, 1, %xmm14)

	mulps	 %xmm5, %xmm0
	addps	 %xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	movss	 %xmm12,  %xmm1
	shufps	 $0x93, %xmm12, %xmm1
	mulps	 %xmm5, %xmm1
	addps	 %xmm1, %xmm9
	MOVUPS_A2 (-23 * SIZE, A1, LDA, 1, %xmm1)

	shufps	 $0x4e, %xmm13, %xmm2
	mulps	 %xmm5, %xmm2
	addps	 %xmm2, %xmm10
	MOVUPS_A1 (-22 * SIZE, A2, %xmm2)
	movss	 %xmm14,  %xmm3
	shufps	 $0x39,  %xmm3, %xmm3
	mulps	 %xmm5,  %xmm3
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	 %xmm3, %xmm11
	MOVUPS_A2 (-21 * SIZE, A2, LDA, 1, %xmm3)

	mulps	 %xmm4, %xmm0
	addps	 %xmm0, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x93, %xmm1, %xmm12
	mulps	 %xmm4, %xmm12
	addps	 %xmm12, %xmm9
	MOVUPS_A2 (-19 * SIZE, A1, LDA, 1, %xmm12)

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	addps	 %xmm13, %xmm10
	MOVUPS_A1 (-18 * SIZE, A2, %xmm13)
	movss	 %xmm3,  %xmm14
	shufps	 $0x39,  %xmm14, %xmm14
	mulps	 %xmm4,  %xmm14
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	 %xmm14, %xmm11
	MOVUPS_A2 (-17 * SIZE, A2, LDA, 1, %xmm14)

	mulps	 %xmm5, %xmm0
	addps	 %xmm0, %xmm8
	movss	 %xmm12,  %xmm1
	shufps	 $0x93, %xmm12, %xmm1
	mulps	 %xmm5, %xmm1
	addps	 %xmm1, %xmm9

	shufps	 $0x4e, %xmm13, %xmm2
	mulps	 %xmm5, %xmm2
	addps	 %xmm2, %xmm10
	movss	 %xmm14,  %xmm3
	shufps	 $0x39,  %xmm3, %xmm3
	mulps	 %xmm5,  %xmm3
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	 %xmm3, %xmm11

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1
	ALIGN_4

L(305):
	testq	$8, MM
	jle	L(306)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-31 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A1 (-30 * SIZE, A2, %xmm2)
	MOVUPS_A2 (-29 * SIZE, A2, LDA, 1, %xmm3)

	mulps	 %xmm4, %xmm0
	addps	 %xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x93, %xmm1, %xmm12
	mulps	 %xmm4, %xmm12
	addps	 %xmm12, %xmm9
	MOVUPS_A2 (-27 * SIZE, A1, LDA, 1, %xmm12)

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	addps	 %xmm13, %xmm10
	MOVUPS_A1 (-26 * SIZE, A2, %xmm13)
	movss	 %xmm3,  %xmm14
	shufps	 $0x39,  %xmm14, %xmm14
	mulps	 %xmm4,  %xmm14
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	 %xmm14, %xmm11
	MOVUPS_A2 (-25 * SIZE, A2, LDA, 1, %xmm14)

	mulps	 %xmm5, %xmm0
	addps	 %xmm0, %xmm8
	movss	 %xmm12,  %xmm1
	shufps	 $0x93, %xmm12, %xmm1
	mulps	 %xmm5, %xmm1
	addps	 %xmm1, %xmm9

	shufps	 $0x4e, %xmm13, %xmm2
	mulps	 %xmm5, %xmm2
	addps	 %xmm2, %xmm10
	movss	 %xmm14,  %xmm3
	shufps	 $0x39,  %xmm3, %xmm3
	mulps	 %xmm5,  %xmm3
	addps	 %xmm3, %xmm11

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1
	ALIGN_4

L(306):
	testq	$4, MM
	jle	L(307)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-31 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A1 (-30 * SIZE, A2, %xmm2)
	MOVUPS_A2 (-29 * SIZE, A2, LDA, 1, %xmm3)

	mulps	 %xmm4, %xmm0
	addps	 %xmm0, %xmm8
	movss	 %xmm1,  %xmm12
	shufps	 $0x93, %xmm1, %xmm12
	mulps	 %xmm4, %xmm12
	addps	 %xmm12, %xmm9

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	addps	 %xmm13, %xmm10
	movss	 %xmm3,  %xmm14
	shufps	 $0x39,  %xmm14, %xmm14
	mulps	 %xmm4,  %xmm14
	addps	 %xmm14, %xmm11

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, A2
	addq	$4 * SIZE, X1
	ALIGN_4

L(307):
	testq	$2, MM
	jle	L(308)

#ifdef movsd
	xorps	%xmm0, %xmm0
#endif
	movsd	-32 * SIZE(A1), %xmm0
#ifdef movsd
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A1, LDA), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
#ifdef movsd
	xorps	%xmm2, %xmm2
#endif
	movsd	-32 * SIZE(A2), %xmm2
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10
#ifdef movsd
	xorps	%xmm3, %xmm3
#endif
	movsd	-32 * SIZE(A2, LDA), %xmm3
	mulps	%xmm4, %xmm3
	addps	%xmm3, %xmm11
	shufps	$0xe, %xmm4, %xmm4

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, A2
	addq	$2 * SIZE, X1
	ALIGN_4

L(308):
	testq	$1, MM
	jle	L(309)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A1, LDA), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	movss	-32 * SIZE(A2), %xmm2
	mulss	%xmm4, %xmm2
	addss	%xmm2, %xmm10
	movss	-32 * SIZE(A2, LDA), %xmm3
	mulss	%xmm4, %xmm3
	addss	%xmm3, %xmm11
	ALIGN_4

L(309):
#ifdef HAVE_SSE3
	haddps	%xmm9, %xmm8
	haddps	%xmm11, %xmm10
	haddps	%xmm10, %xmm8

	pshufd	$0x1, %xmm8, %xmm9
	pshufd	$0x2, %xmm8, %xmm10
	pshufd	$0x3, %xmm8, %xmm11
#else
	movaps	%xmm8, %xmm0
	unpcklps %xmm9, %xmm8
	unpckhps %xmm9, %xmm0

	movaps	%xmm10, %xmm1
	unpcklps %xmm11, %xmm10
	unpckhps %xmm11, %xmm1

	movaps	%xmm8, %xmm9
	unpcklps %xmm10, %xmm8
	unpckhps %xmm10, %xmm9

	movaps	%xmm0, %xmm10
	unpcklps %xmm1, %xmm0
	unpckhps %xmm1, %xmm10

	addps	%xmm9, %xmm8
	addps	%xmm0, %xmm10
	addps	%xmm10, %xmm8

	pshufd	$0x2, %xmm8, %xmm9
	pshufd	$0x1, %xmm8, %xmm10
	pshufd	$0x3, %xmm8, %xmm11
#endif

	mulss	ALPHA, %xmm8
	mulss	ALPHA, %xmm9
	mulss	ALPHA, %xmm10
	mulss	ALPHA, %xmm11

	addss	(Y), %xmm8
	addq	INCY, Y
	addss	(Y), %xmm9
	addq	INCY, Y
	addss	(Y), %xmm10
	addq	INCY, Y
	addss	(Y), %xmm11
	addq	INCY, Y

	movss	%xmm8, (Y1)
	addq	INCY, Y1
	movss	%xmm9, (Y1)
	addq	INCY, Y1
	movss	%xmm10, (Y1)
	addq	INCY, Y1
	movss	%xmm11, (Y1)
	addq	INCY, Y1

	cmpq	$4, N
	jge	L(301)
	ALIGN_4

L(310):
	testq	N, N
	jle	L(999)

	cmpq	$3, N
	jne	L(320)

	leaq	32 * SIZE(BUFFER), X1

	movq	A, A1
	leaq	(A1, LDA, 2), A2
	leaq	(A1, LDA, 4), A

	xorps	%xmm8, %xmm8
	xorps	%xmm9, %xmm9
	xorps	%xmm10, %xmm10

	cmpq	$3, M
	jle	L(317)

	testq	$SIZE, A1
	je	L(31X)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A1, LDA), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	movss	-32 * SIZE(A2), %xmm2
	mulss	%xmm4, %xmm2
	addss	%xmm2, %xmm10

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, X1
	ALIGN_3

L(31X):
	testq	$2 * SIZE, A1
	je	L(31XX)

#ifdef movsd
	xorps	%xmm0, %xmm0
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(A1), %xmm0
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A1, LDA), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
#ifdef movsd
	xorps	%xmm2, %xmm2
#endif
	movsd	-32 * SIZE(A2), %xmm2
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, X1
	ALIGN_3

L(31XX):
	movaps	-35 * SIZE(A1, LDA), %xmm12
	movaps	-34 * SIZE(A2), %xmm13

	MOVUPS_XL1 (-32 * SIZE, X1, %xmm4)
	MOVUPS_XL1 (-28 * SIZE, X1, %xmm5)

	movq	MM,  I
	sarq	$4,  I
	jle	L(315)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-31 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A1 (-30 * SIZE, A2, %xmm2)

	decq	I
	jle	L(313)
	ALIGN_4

L(312):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A1)
#endif

	mulps	 %xmm4, %xmm0
	addps	 %xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x93, %xmm1, %xmm12
	mulps	 %xmm4, %xmm12
	addps	 %xmm12, %xmm9
	MOVUPS_A2 (-27 * SIZE, A1, LDA, 1, %xmm12)

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	 %xmm13, %xmm10
	MOVUPS_A1 (-26 * SIZE, A2, %xmm13)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A1, LDA)
#endif

	mulps	 %xmm5, %xmm0
	addps	 %xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	movss	 %xmm12,  %xmm1
	shufps	 $0x93, %xmm12, %xmm1
	mulps	 %xmm5, %xmm1
	addps	 %xmm1, %xmm9
	MOVUPS_A2 (-23 * SIZE, A1, LDA, 1, %xmm1)
	shufps	 $0x4e, %xmm13, %xmm2
	mulps	 %xmm5, %xmm2
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	 %xmm2, %xmm10
	MOVUPS_A1 (-22 * SIZE, A2, %xmm2)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(A2)
#endif

	mulps	 %xmm4, %xmm0
	addps	 %xmm0, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x93, %xmm1, %xmm12
	mulps	 %xmm4, %xmm12
	addps	 %xmm12, %xmm9
	MOVUPS_A2 (-19 * SIZE, A1, LDA, 1, %xmm12)
	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	 %xmm13, %xmm10
	MOVUPS_A1 (-18 * SIZE, A2, %xmm13)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 3 - 128 + PREOFFSET(X1)
#endif

	mulps	 %xmm5, %xmm0
	addps	 %xmm0, %xmm8
	MOVUPS_A1 (-16 * SIZE, A1, %xmm0)
	movss	 %xmm12,  %xmm1
	shufps	 $0x93, %xmm12, %xmm1
	mulps	 %xmm5, %xmm1
	addps	 %xmm1, %xmm9
	MOVUPS_A2 (-15 * SIZE, A1, LDA, 1, %xmm1)
	shufps	 $0x4e, %xmm13, %xmm2
	mulps	 %xmm5, %xmm2
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	 %xmm2, %xmm10
	MOVUPS_A1 (-14 * SIZE, A2, %xmm2)

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1

	decq	I
	jg	L(312)
	ALIGN_4

L(313):
	mulps	 %xmm4, %xmm0
	addps	 %xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x93, %xmm1, %xmm12
	mulps	 %xmm4, %xmm12
	addps	 %xmm12, %xmm9
	MOVUPS_A2 (-27 * SIZE, A1, LDA, 1, %xmm12)

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	addps	 %xmm13, %xmm10
	MOVUPS_A1 (-26 * SIZE, A2, %xmm13)
	movss	 %xmm3,  %xmm14
	shufps	 $0x39,  %xmm14, %xmm14
	mulps	 %xmm4,  %xmm14
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	 %xmm14, %xmm11
	MOVUPS_A2 (-25 * SIZE, A2, LDA, 1, %xmm14)

	mulps	 %xmm5, %xmm0
	addps	 %xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	movss	 %xmm12,  %xmm1
	shufps	 $0x93, %xmm12, %xmm1
	mulps	 %xmm5, %xmm1
	addps	 %xmm1, %xmm9
	MOVUPS_A2 (-23 * SIZE, A1, LDA, 1, %xmm1)

	shufps	 $0x4e, %xmm13, %xmm2
	mulps	 %xmm5, %xmm2
	addps	 %xmm2, %xmm10
	MOVUPS_A1 (-22 * SIZE, A2, %xmm2)
	movss	 %xmm14,  %xmm3
	shufps	 $0x39,  %xmm3, %xmm3
	mulps	 %xmm5,  %xmm3
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	 %xmm3, %xmm11
	MOVUPS_A2 (-21 * SIZE, A2, LDA, 1, %xmm3)

	mulps	 %xmm4, %xmm0
	addps	 %xmm0, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x93, %xmm1, %xmm12
	mulps	 %xmm4, %xmm12
	addps	 %xmm12, %xmm9
	MOVUPS_A2 (-19 * SIZE, A1, LDA, 1, %xmm12)

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	addps	 %xmm13, %xmm10
	MOVUPS_A1 (-18 * SIZE, A2, %xmm13)
	movss	 %xmm3,  %xmm14
	shufps	 $0x39,  %xmm14, %xmm14
	mulps	 %xmm4,  %xmm14
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	 %xmm14, %xmm11
	MOVUPS_A2 (-17 * SIZE, A2, LDA, 1, %xmm14)

	mulps	 %xmm5, %xmm0
	addps	 %xmm0, %xmm8
	movss	 %xmm12,  %xmm1
	shufps	 $0x93, %xmm12, %xmm1
	mulps	 %xmm5, %xmm1
	addps	 %xmm1, %xmm9

	shufps	 $0x4e, %xmm13, %xmm2
	mulps	 %xmm5, %xmm2
	addps	 %xmm2, %xmm10
	movss	 %xmm14,  %xmm3
	shufps	 $0x39,  %xmm3, %xmm3
	mulps	 %xmm5,  %xmm3
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	 %xmm3, %xmm11

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1
	ALIGN_4

L(315):
	testq	$8, MM
	jle	L(316)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-31 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A1 (-30 * SIZE, A2, %xmm2)

	mulps	 %xmm4, %xmm0
	addps	 %xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x93, %xmm1, %xmm12
	mulps	 %xmm4, %xmm12
	addps	 %xmm12, %xmm9
	MOVUPS_A2 (-27 * SIZE, A1, LDA, 1, %xmm12)

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	 %xmm13, %xmm10
	MOVUPS_A1 (-26 * SIZE, A2, %xmm13)

	mulps	 %xmm5, %xmm0
	addps	 %xmm0, %xmm8
	movss	 %xmm12,  %xmm1
	shufps	 $0x93, %xmm12, %xmm1
	mulps	 %xmm5, %xmm1
	addps	 %xmm1, %xmm9

	shufps	 $0x4e, %xmm13, %xmm2
	mulps	 %xmm5, %xmm2
	addps	 %xmm2, %xmm10

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1
	ALIGN_4

L(316):
	testq	$4, MM
	jle	L(317)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A2 (-31 * SIZE, A1, LDA, 1, %xmm1)
	MOVUPS_A1 (-30 * SIZE, A2, %xmm2)

	mulps	 %xmm4, %xmm0
	addps	 %xmm0, %xmm8
	movss	 %xmm1,  %xmm12
	shufps	 $0x93, %xmm1, %xmm12
	mulps	 %xmm4, %xmm12
	addps	 %xmm12, %xmm9

	shufps	 $0x4e, %xmm2, %xmm13
	mulps	 %xmm4, %xmm13
	addps	 %xmm13, %xmm10

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, A2
	addq	$4 * SIZE, X1
	ALIGN_4

L(317):
	testq	$2, MM
	jle	L(318)

#ifdef movsd
	xorps	%xmm0, %xmm0
#endif
	movsd	-32 * SIZE(A1), %xmm0
#ifdef movsd
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A1, LDA), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
#ifdef movsd
	xorps	%xmm2, %xmm2
#endif
	movsd	-32 * SIZE(A2), %xmm2
	mulps	%xmm4, %xmm2
	addps	%xmm2, %xmm10

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, A2
	addq	$2 * SIZE, X1
	ALIGN_4

L(318):
	testq	$1, MM
	jle	L(319)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A1, LDA), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	movss	-32 * SIZE(A2), %xmm2
	mulss	%xmm4, %xmm2
	addss	%xmm2, %xmm10
	ALIGN_4

L(319):
#ifdef HAVE_SSE3
	haddps	%xmm9, %xmm8
	haddps	%xmm11, %xmm10
	haddps	%xmm10, %xmm8

	pshufd	$0x1, %xmm8, %xmm9
	pshufd	$0x2, %xmm8, %xmm10
#else
	movaps	%xmm8, %xmm0
	unpcklps %xmm9, %xmm8
	unpckhps %xmm9, %xmm0

	movaps	%xmm10, %xmm1
	unpcklps %xmm11, %xmm10
	unpckhps %xmm11, %xmm1

	movaps	%xmm8, %xmm9
	unpcklps %xmm10, %xmm8
	unpckhps %xmm10, %xmm9

	movaps	%xmm0, %xmm10
	unpcklps %xmm1, %xmm0
	unpckhps %xmm1, %xmm10

	addps	%xmm9, %xmm8
	addps	%xmm0, %xmm10
	addps	%xmm10, %xmm8

	pshufd	$0x2, %xmm8, %xmm9
	pshufd	$0x1, %xmm8, %xmm10
#endif

	mulss	ALPHA, %xmm8
	mulss	ALPHA, %xmm9
	mulss	ALPHA, %xmm10

	addss	(Y), %xmm8
	addq	INCY, Y
	addss	(Y), %xmm9
	addq	INCY, Y
	addss	(Y), %xmm10

	movss	%xmm8, (Y1)
	addq	INCY, Y1
	movss	%xmm9, (Y1)
	addq	INCY, Y1
	movss	%xmm10, (Y1)
	jmp	L(999)
	ALIGN_3

L(320):
	cmpq	$2, N
	jne	L(330)

	leaq	32 * SIZE(BUFFER), X1

	movq	A, A1
	leaq	(A1, LDA), A2

	xorps	%xmm8, %xmm8
	xorps	%xmm9, %xmm9

	cmpq	$3, M
	jle	L(327)

	testq	$SIZE, A1
	je	L(32X)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A2), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, A2
	addq	 $1 * SIZE, X1
	ALIGN_3

L(32X):
	testq	$2 * SIZE, A1
	je	L(32XX)

#ifdef movsd
	xorps	%xmm0, %xmm0
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(A1), %xmm0
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A2), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, A2
	addq	 $2 * SIZE, X1
	ALIGN_3

L(32XX):
	movaps	-35 * SIZE(A2), %xmm12

	MOVUPS_XL1 (-32 * SIZE, X1, %xmm4)
	MOVUPS_XL1 (-28 * SIZE, X1, %xmm5)

	movq	MM,  I
	sarq	$4,  I
	jle	L(325)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A1 (-31 * SIZE, A2, %xmm1)

	decq	I
	jle	L(323)
	ALIGN_4

L(322):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A1)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x93, %xmm1, %xmm12
	mulps	 %xmm4, %xmm12
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	 %xmm12, %xmm9
	MOVUPS_A1 (-27 * SIZE, A2, %xmm12)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	movss	 %xmm12,  %xmm1
	shufps	 $0x93, %xmm12, %xmm1
	mulps	 %xmm5, %xmm1
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	 %xmm1, %xmm9
	MOVUPS_A1 (-23 * SIZE, A2, %xmm1)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A2)
#endif

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x93, %xmm1, %xmm12
	mulps	 %xmm4, %xmm12
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	 %xmm12, %xmm9
	MOVUPS_A1 (-19 * SIZE, A2, %xmm12)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(X1)
#endif

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-16 * SIZE, A1, %xmm0)
	movss	 %xmm12,  %xmm1
	shufps	 $0x93, %xmm12, %xmm1
	mulps	 %xmm5, %xmm1
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	 %xmm1, %xmm9
	MOVUPS_A1 (-15 * SIZE, A2, %xmm1)

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1

	decq	I
	jg	L(322)
	ALIGN_4

L(323):
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x93, %xmm1, %xmm12
	mulps	 %xmm4, %xmm12
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	 %xmm12, %xmm9
	MOVUPS_A1 (-27 * SIZE, A2, %xmm12)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)
	movss	 %xmm12,  %xmm1
	shufps	 $0x93, %xmm12, %xmm1
	mulps	 %xmm5, %xmm1
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	 %xmm1, %xmm9
	MOVUPS_A1 (-23 * SIZE, A2, %xmm1)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-20 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x93, %xmm1, %xmm12
	mulps	 %xmm4, %xmm12
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	 %xmm12, %xmm9
	MOVUPS_A1 (-19 * SIZE, A2, %xmm12)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	movss	 %xmm12,  %xmm1
	shufps	 $0x93, %xmm12, %xmm1
	mulps	 %xmm5, %xmm1
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	 %xmm1, %xmm9

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, A2
	addq	$16 * SIZE, X1
	ALIGN_4

L(325):
	testq	$8, MM
	jle	L(326)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A1 (-31 * SIZE, A2, %xmm1)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-28 * SIZE, A1, %xmm0)
	movss	 %xmm1,  %xmm12
	shufps	 $0x93, %xmm1, %xmm12
	mulps	 %xmm4, %xmm12
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	 %xmm12, %xmm9
	MOVUPS_A1 (-27 * SIZE, A2, %xmm12)

	mulps	%xmm5, %xmm0
	addps	%xmm0, %xmm8
	movss	 %xmm12,  %xmm1
	shufps	 $0x93, %xmm12, %xmm1
	mulps	 %xmm5, %xmm1
	addps	 %xmm1, %xmm9

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1
	ALIGN_4

L(326):
	testq	$4, MM
	jle	L(327)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A1 (-31 * SIZE, A2, %xmm1)

	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	movss	 %xmm1,  %xmm12
	shufps	 $0x93, %xmm1, %xmm12
	mulps	 %xmm4, %xmm12
	addps	 %xmm12, %xmm9

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, A2
	addq	$4 * SIZE, X1
	ALIGN_4

L(327):
	testq	$2, MM
	jle	L(328)

#ifdef movsd
	xorps	%xmm0, %xmm0
#endif
	movsd	-32 * SIZE(A1), %xmm0
#ifdef movsd
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
#ifdef movsd
	xorps	%xmm1, %xmm1
#endif
	movsd	-32 * SIZE(A2), %xmm1
	mulps	%xmm4, %xmm1
	addps	%xmm1, %xmm9
	shufps	$0xe, %xmm4, %xmm4

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, A2
	addq	$2 * SIZE, X1
	ALIGN_4

L(328):
	testq	$1, MM
	jle	L(329)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	movss	-32 * SIZE(A2), %xmm1
	mulss	%xmm4, %xmm1
	addss	%xmm1, %xmm9
	ALIGN_4

L(329):
#ifdef HAVE_SSE3
	haddps	%xmm9, %xmm8
	haddps	%xmm8, %xmm8
#else
	movaps	%xmm8, %xmm10
	unpcklps %xmm9, %xmm8
	unpckhps %xmm9, %xmm10

	addps	%xmm10, %xmm8
	movhlps %xmm8, %xmm9
	addps	%xmm9, %xmm8
#endif

	pshufd	$0x1, %xmm8, %xmm9

	mulss	ALPHA, %xmm8
	mulss	ALPHA, %xmm9

	addss	(Y), %xmm8
	addq	INCY, Y
	addss	(Y), %xmm9
	addq	INCY, Y

	movss	%xmm8, (Y1)
	addq	INCY, Y1
	movss	%xmm9, (Y1)
	addq	INCY, Y1
	jmp	L(999)
	ALIGN_4

L(330):
	cmpq	$1, N
	jne	L(999)

	leaq	32 * SIZE(BUFFER), X1

	movq	A, A1

	xorps	%xmm8, %xmm8
	xorps	%xmm9, %xmm9

	cmpq	$3, M
	jle	L(337)

	testq	$SIZE, A1
	je	L(33X)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8

	addq	 $1 * SIZE, A1
	addq	 $1 * SIZE, X1
	ALIGN_3

L(33X):
	testq	$2 * SIZE, A1
	je	L(33XX)

#ifdef movsd
	xorps	%xmm0, %xmm0
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(A1), %xmm0
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	shufps	$0xe, %xmm4, %xmm4

	addq	 $2 * SIZE, A1
	addq	 $2 * SIZE, X1
	ALIGN_3

L(33XX):

	MOVUPS_XL1 (-32 * SIZE, X1, %xmm4)
	MOVUPS_XL1 (-28 * SIZE, X1, %xmm5)

	movq	MM,  I
	sarq	$4,  I
	jle	L(335)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	MOVUPS_A1 (-28 * SIZE, A1, %xmm12)

	decq	I
	jle	L(333)
	ALIGN_4

L(332):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(A1)
#endif

	mulps	%xmm4, %xmm0
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)

	mulps	%xmm5, %xmm12
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-20 * SIZE, A1, %xmm12)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(X1)
#endif

	mulps	%xmm4, %xmm0
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-16 * SIZE, A1, %xmm0)

	mulps	%xmm5, %xmm12
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-12 * SIZE, A1, %xmm12)

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, X1

	decq	I
	jg	L(332)
	ALIGN_4

L(333):
	mulps	%xmm4, %xmm0
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm0, %xmm8
	MOVUPS_A1 (-24 * SIZE, A1, %xmm0)

	mulps	%xmm5, %xmm12
	MOVUPS_XL1 (-20 * SIZE, X1, %xmm5)
	addps	%xmm12, %xmm9
	MOVUPS_A1 (-20 * SIZE, A1, %xmm12)

	mulps	%xmm4, %xmm0
	MOVUPS_XL1 (-16 * SIZE, X1, %xmm4)
	addps	%xmm0, %xmm8

	mulps	%xmm5, %xmm12
	MOVUPS_XL1 (-12 * SIZE, X1, %xmm5)
	addps	%xmm12, %xmm9

	addq	$16 * SIZE, A1
	addq	$16 * SIZE, X1
	ALIGN_4

L(335):
	testq	$8, MM
	jle	L(336)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm0
	MOVUPS_XL1 (-24 * SIZE, X1, %xmm4)
	addps	%xmm0, %xmm8

	MOVUPS_A1 (-28 * SIZE, A1, %xmm12)
	mulps	%xmm5, %xmm12
	addps	%xmm12, %xmm9

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, X1
	ALIGN_4

L(336):
	testq	$4, MM
	jle	L(337)

	MOVUPS_A1 (-32 * SIZE, A1, %xmm0)
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, X1
	ALIGN_4

L(337):
	testq	$2, MM
	jle	L(338)

#ifdef movsd
	xorps	%xmm0, %xmm0
#endif
	movsd	-32 * SIZE(A1), %xmm0
#ifdef movsd
	xorps	%xmm4, %xmm4
#endif
	movsd	-32 * SIZE(X1), %xmm4
	mulps	%xmm4, %xmm0
	addps	%xmm0, %xmm8
	shufps	$0xe, %xmm4, %xmm4

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, X1
	ALIGN_4

L(338):
	testq	$1, MM
	jle	L(339)

	movss	-32 * SIZE(A1), %xmm0
	movss	-32 * SIZE(X1), %xmm4
	mulss	%xmm4, %xmm0
	addss	%xmm0, %xmm8
	ALIGN_4

L(339):
	addps	%xmm9, %xmm8

#ifdef HAVE_SSE3
	haddps	%xmm8, %xmm8
	haddps	%xmm8, %xmm8
#else
	pshufd	$1, %xmm8, %xmm9
	pshufd	$2, %xmm8, %xmm10
	pshufd	$3, %xmm8, %xmm11

	addss	%xmm9, %xmm8
	addss	%xmm11, %xmm10
	addss	%xmm10, %xmm8
#endif

	mulss	ALPHA, %xmm8

	addss	(Y), %xmm8
	movss	%xmm8, (Y1)

	jmp	L(999)
#endif
	ALIGN_4

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
	ALIGN_4

	EPILOGUE
