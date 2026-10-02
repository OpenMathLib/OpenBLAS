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

#if GEMV_UNROLL < 2
#undef  GEMV_UNROLL
#define GEMV_UNROLL 2
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
#define MMM	56(%rsp)
#define NN	64(%rsp)
#define AA	72(%rsp)
#define LDAX	80(%rsp)
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
//Temp variables for M,N,A,LDA
#define MMM	224(%rsp)
#define NN	232(%rsp)
#define AA	240(%rsp)
#define LDAX	248(%rsp)

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
#define X1	%r15

#ifdef ALIGNED_ACCESS
#define MM	INCX
#else
#define MM	M
#endif

#define ALPHA	%xmm15

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

	movq	M,	      MMM
	movq	N,            NN
	movq	A,            AA
	movq	LDA,	      LDAX

#else
	movq	OLD_M,	      MMM
	movq	OLD_N,        NN
	movq	OLD_A,        AA
	movq	OLD_LDA,      LDAX
#endif
#ifdef HAVE_SSE3
#ifndef WINDOWS_ABI
	movddup	%xmm0, ALPHA
#else
	movddup	%xmm3, ALPHA
#endif
#else
#ifndef WINDOWS_ABI
	movapd	%xmm0, ALPHA
#else
	movapd	%xmm3, ALPHA
#endif
	unpcklpd ALPHA, ALPHA
#endif



L(0x):
	xorq	M,M
	addq	$1,M
	salq	$21,M
	subq	M,MMM
	jge L(00)

	movq	MMM,%rax
	addq	M,%rax
	jle	L(999x)
	movq	%rax,M

L(00):
	movq	LDAX,LDA
	movq	NN,N
	movq	AA,A
	movq	STACK_INCX,   INCX
	movq	STACK_Y,      Y
	movq	STACK_INCY,   INCY
	movq	STACK_BUFFER, BUFFER

	leaq	-1(INCX), %rax

	leaq	(,LDA,  SIZE), LDA
	leaq	(,INCX, SIZE), INCX
	leaq	(,INCY, SIZE), INCY

	leaq	(LDA, LDA, 2), LDA3

	subq	$-16 * SIZE, A

	testq	M, M
	jle	L(999)
	testq	N, N
	jle	L(999)

	movq	BUFFER, X1

#ifdef ALIGNED_ACCESS
	testq	$SIZE, A
	je	L(01)

	movsd	(X), %xmm0
	addq	INCX, X

	movsd	%xmm0, 1 * SIZE(BUFFER)
	addq	$1 * SIZE, BUFFER
	addq	$2 * SIZE, X1
	decq	M
	jle	L(10)
	ALIGN_4

L(01):
#endif

	movq	M,  I
	sarq	$3, I
	jle	L(05)
	ALIGN_4

L(02):
	movsd	(X), %xmm0
	addq	INCX, X
	movhpd	(X), %xmm0
	addq	INCX, X

	movsd	(X), %xmm1
	addq	INCX, X
	movhpd	(X), %xmm1
	addq	INCX, X

	movsd	(X), %xmm2
	addq	INCX, X
	movhpd	(X), %xmm2
	addq	INCX, X

	movsd	(X), %xmm3
	addq	INCX, X
	movhpd	(X), %xmm3
	addq	INCX, X

	movapd	%xmm0, 0 * SIZE(X1)
	movapd	%xmm1, 2 * SIZE(X1)
	movapd	%xmm2, 4 * SIZE(X1)
	movapd	%xmm3, 6 * SIZE(X1)

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
	movsd	(X), %xmm0
	addq	INCX, X
	movsd	%xmm0, 0 * SIZE(X1)
	addq	$SIZE, X1
	decq	I
	jg	L(06)
	ALIGN_4

L(10):
	movq	Y, Y1

#ifdef ALIGNED_ACCESS
	testq	$SIZE, LDA
	jne	L(50)
#endif

#if GEMV_UNROLL >= 8
	cmpq	$8, N
	jl	L(20)
	ALIGN_3

L(11):
	subq	$8, N

	leaq	16 * SIZE(BUFFER), X1

	movq	A, A1
	leaq	(A1, LDA, 4), A2
	leaq	(A1, LDA, 8), A

	xorps	%xmm0, %xmm0
	xorps	%xmm1, %xmm1
	xorps	%xmm2, %xmm2
	xorps	%xmm3, %xmm3
	xorps	%xmm4, %xmm4
	xorps	%xmm5, %xmm5
	xorps	%xmm6, %xmm6
	xorps	%xmm7, %xmm7

#ifdef PREFETCHW
	PREFETCHW	7 * SIZE(Y1)
#endif

#ifdef ALIGNED_ACCESS
	testq	$SIZE, A
	je	L(1X)

	movsd	-16 * SIZE(X1), %xmm12

	movsd	-16 * SIZE(A1), %xmm8
	mulsd	%xmm12, %xmm8
	addsd	%xmm8, %xmm0
	movsd	-16 * SIZE(A1, LDA), %xmm9
	mulsd	%xmm12, %xmm9
	addsd	%xmm9, %xmm1
	movsd	-16 * SIZE(A1, LDA, 2), %xmm10
	mulsd	%xmm12, %xmm10
	addsd	%xmm10, %xmm2
	movsd	-16 * SIZE(A1, LDA3), %xmm11
	mulsd	%xmm12, %xmm11
	addsd	%xmm11, %xmm3
	movsd	-16 * SIZE(A2), %xmm8
	mulsd	%xmm12, %xmm8
	addsd	%xmm8, %xmm4
	movsd	-16 * SIZE(A2, LDA), %xmm9
	mulsd	%xmm12, %xmm9
	addsd	%xmm9, %xmm5
	movsd	-16 * SIZE(A2, LDA, 2), %xmm10
	mulsd	%xmm12, %xmm10
	addsd	%xmm10, %xmm6
	movsd	-16 * SIZE(A2, LDA3), %xmm11
	mulsd	%xmm12, %xmm11
	addsd	%xmm11, %xmm7

	addq	 $SIZE, A1
	addq	 $SIZE, A2
	addq	 $SIZE, X1
	ALIGN_3

L(1X):
#endif

	movq	M,  I
	sarq	$3,  I
	jle	L(15)

	MOVUPS_A1(-16 * SIZE, A1, %xmm8)
	MOVUPS_A2(-16 * SIZE, A1, LDA, 1, %xmm9)
	MOVUPS_A2(-16 * SIZE, A1, LDA, 2, %xmm10)
	MOVUPS_A2(-16 * SIZE, A1, LDA3, 1, %xmm11)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)
	MOVUPS_XL1(-14 * SIZE, X1, %xmm13)

	decq	I
	jle	L(13)
	ALIGN_4

L(12):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A1)
#endif

	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-16 * SIZE, A2, %xmm8)
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A2(-16 * SIZE, A2, LDA, 1, %xmm9)
	mulpd	%xmm12, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A2(-16 * SIZE, A2, LDA, 2, %xmm10)
	mulpd	%xmm12, %xmm11
	addpd	%xmm11, %xmm3
	MOVUPS_A2(-16 * SIZE, A2, LDA3, 1, %xmm11)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A1, LDA)
#endif

	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm4
	MOVUPS_A1(-14 * SIZE, A1, %xmm8)
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm5
	MOVUPS_A2(-14 * SIZE, A1, LDA, 1, %xmm9)
	mulpd	%xmm12, %xmm10
	addpd	%xmm10, %xmm6
	MOVUPS_A2(-14 * SIZE, A1, LDA, 2, %xmm10)
	mulpd	%xmm12, %xmm11
	MOVUPS_XL1(-12 * SIZE, X1, %xmm12)
	addpd	%xmm11, %xmm7
	MOVUPS_A2(-14 * SIZE, A1, LDA3, 1, %xmm11)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A1, LDA, 2)
#endif

	mulpd	%xmm13, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-14 * SIZE, A2, %xmm8)
	mulpd	%xmm13, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A2(-14 * SIZE, A2, LDA, 1, %xmm9)
	mulpd	%xmm13, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A2(-14 * SIZE, A2, LDA, 2, %xmm10)
	mulpd	%xmm13, %xmm11
	addpd	%xmm11, %xmm3
	MOVUPS_A2(-14 * SIZE, A2, LDA3, 1, %xmm11)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A1, LDA3)
#endif

	mulpd	%xmm13, %xmm8
	addpd	%xmm8, %xmm4
	MOVUPS_A1(-12 * SIZE, A1, %xmm8)
	mulpd	%xmm13, %xmm9
	addpd	%xmm9, %xmm5
	MOVUPS_A2(-12 * SIZE, A1, LDA, 1, %xmm9)
	mulpd	%xmm13, %xmm10
	addpd	%xmm10, %xmm6
	MOVUPS_A2(-12 * SIZE, A1, LDA, 2, %xmm10)
	mulpd	%xmm13, %xmm11
	MOVUPS_XL1(-10 * SIZE, X1, %xmm13)
	addpd	%xmm11, %xmm7
	MOVUPS_A2(-12 * SIZE, A1, LDA3, 1, %xmm11)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A2)
#endif

	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-12 * SIZE, A2, %xmm8)
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A2(-12 * SIZE, A2, LDA, 1, %xmm9)
	mulpd	%xmm12, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A2(-12 * SIZE, A2, LDA, 2, %xmm10)
	mulpd	%xmm12, %xmm11
	addpd	%xmm11, %xmm3
	MOVUPS_A2(-12 * SIZE, A2, LDA3, 1, %xmm11)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A2, LDA)
#endif

	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm4
	MOVUPS_A1(-10 * SIZE, A1, %xmm8)
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm5
	MOVUPS_A2(-10 * SIZE, A1, LDA, 1, %xmm9)
	mulpd	%xmm12, %xmm10
	addpd	%xmm10, %xmm6
	MOVUPS_A2(-10 * SIZE, A1, LDA, 2, %xmm10)
	mulpd	%xmm12, %xmm11
	MOVUPS_XL1(-8 * SIZE, X1, %xmm12)
	addpd	%xmm11, %xmm7
	MOVUPS_A2(-10 * SIZE, A1, LDA3, 1, %xmm11)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A2, LDA, 2)
#endif

	mulpd	%xmm13, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-10 * SIZE, A2, %xmm8)
	mulpd	%xmm13, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A2(-10 * SIZE, A2, LDA, 1, %xmm9)
	mulpd	%xmm13, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A2(-10 * SIZE, A2, LDA, 2, %xmm10)
	mulpd	%xmm13, %xmm11
	addpd	%xmm11, %xmm3
	MOVUPS_A2(-10 * SIZE, A2, LDA3, 1, %xmm11)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(A2, LDA3)
#endif

	mulpd	%xmm13, %xmm8
	addpd	%xmm8, %xmm4
	MOVUPS_A1(-8 * SIZE, A1, %xmm8)
	mulpd	%xmm13, %xmm9
	addpd	%xmm9, %xmm5
	MOVUPS_A2(-8 * SIZE, A1, LDA, 1, %xmm9)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) - 128 + PREOFFSET(X1)
#endif

	mulpd	%xmm13, %xmm10
	addpd	%xmm10, %xmm6
	MOVUPS_A2(-8 * SIZE, A1, LDA, 2, %xmm10)
	mulpd	%xmm13, %xmm11
	MOVUPS_XL1(-6 * SIZE, X1, %xmm13)
	addpd	%xmm11, %xmm7
	MOVUPS_A2(-8 * SIZE, A1, LDA3, 1, %xmm11)

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1

	decq	I
	jg	L(12)
	ALIGN_4

L(13):
	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-16 * SIZE, A2, %xmm8)
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A2(-16 * SIZE, A2, LDA, 1, %xmm9)
	mulpd	%xmm12, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A2(-16 * SIZE, A2, LDA, 2, %xmm10)
	mulpd	%xmm12, %xmm11
	addpd	%xmm11, %xmm3
	MOVUPS_A2(-16 * SIZE, A2, LDA3, 1, %xmm11)

	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm4
	MOVUPS_A1(-14 * SIZE, A1, %xmm8)
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm5
	MOVUPS_A2(-14 * SIZE, A1, LDA, 1, %xmm9)
	mulpd	%xmm12, %xmm10
	addpd	%xmm10, %xmm6
	MOVUPS_A2(-14 * SIZE, A1, LDA, 2, %xmm10)
	mulpd	%xmm12, %xmm11
	MOVUPS_XL1(-12 * SIZE, X1, %xmm12)
	addpd	%xmm11, %xmm7
	MOVUPS_A2(-14 * SIZE, A1, LDA3, 1, %xmm11)

	mulpd	%xmm13, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-14 * SIZE, A2, %xmm8)
	mulpd	%xmm13, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A2(-14 * SIZE, A2, LDA, 1, %xmm9)
	mulpd	%xmm13, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A2(-14 * SIZE, A2, LDA, 2, %xmm10)
	mulpd	%xmm13, %xmm11
	addpd	%xmm11, %xmm3
	MOVUPS_A2(-14 * SIZE, A2, LDA3, 1, %xmm11)

	mulpd	%xmm13, %xmm8
	addpd	%xmm8, %xmm4
	MOVUPS_A1(-12 * SIZE, A1, %xmm8)
	mulpd	%xmm13, %xmm9
	addpd	%xmm9, %xmm5
	MOVUPS_A2(-12 * SIZE, A1, LDA, 1, %xmm9)
	mulpd	%xmm13, %xmm10
	addpd	%xmm10, %xmm6
	MOVUPS_A2(-12 * SIZE, A1, LDA, 2, %xmm10)
	mulpd	%xmm13, %xmm11
	MOVUPS_XL1(-10 * SIZE, X1, %xmm13)
	addpd	%xmm11, %xmm7
	MOVUPS_A2(-12 * SIZE, A1, LDA3, 1, %xmm11)

	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-12 * SIZE, A2, %xmm8)
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A2(-12 * SIZE, A2, LDA, 1, %xmm9)
	mulpd	%xmm12, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A2(-12 * SIZE, A2, LDA, 2, %xmm10)
	mulpd	%xmm12, %xmm11
	addpd	%xmm11, %xmm3
	MOVUPS_A2(-12 * SIZE, A2, LDA3, 1, %xmm11)

	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm4
	MOVUPS_A1(-10 * SIZE, A1, %xmm8)
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm5
	MOVUPS_A2(-10 * SIZE, A1, LDA, 1, %xmm9)
	mulpd	%xmm12, %xmm10
	addpd	%xmm10, %xmm6
	MOVUPS_A2(-10 * SIZE, A1, LDA, 2, %xmm10)
	mulpd	%xmm12, %xmm11
	addpd	%xmm11, %xmm7
	MOVUPS_A2(-10 * SIZE, A1, LDA3, 1, %xmm11)

	mulpd	%xmm13, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-10 * SIZE, A2, %xmm8)
	mulpd	%xmm13, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A2(-10 * SIZE, A2, LDA, 1, %xmm9)
	mulpd	%xmm13, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A2(-10 * SIZE, A2, LDA, 2, %xmm10)
	mulpd	%xmm13, %xmm11
	addpd	%xmm11, %xmm3
	MOVUPS_A2(-10 * SIZE, A2, LDA3, 1, %xmm11)

	mulpd	%xmm13, %xmm8
	addpd	%xmm8, %xmm4
	mulpd	%xmm13, %xmm9
	addpd	%xmm9, %xmm5
	mulpd	%xmm13, %xmm10
	addpd	%xmm10, %xmm6
	mulpd	%xmm13, %xmm11
	addpd	%xmm11, %xmm7

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1
	ALIGN_4

L(15):
	testq	$4, M
	jle	L(16)

	MOVUPS_A1(-16 * SIZE, A1, %xmm8)
	MOVUPS_A2(-16 * SIZE, A1, LDA, 1, %xmm9)
	MOVUPS_A2(-16 * SIZE, A1, LDA, 2, %xmm10)
	MOVUPS_A2(-16 * SIZE, A1, LDA3, 1, %xmm11)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)
	MOVUPS_XL1(-14 * SIZE, X1, %xmm13)

	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-16 * SIZE, A2, %xmm8)
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A2(-16 * SIZE, A2, LDA, 1, %xmm9)
	mulpd	%xmm12, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A2(-16 * SIZE, A2, LDA, 2, %xmm10)
	mulpd	%xmm12, %xmm11
	addpd	%xmm11, %xmm3
	MOVUPS_A2(-16 * SIZE, A2, LDA3, 1, %xmm11)

	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm4
	MOVUPS_A1(-14 * SIZE, A1, %xmm8)
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm5
	MOVUPS_A2(-14 * SIZE, A1, LDA, 1, %xmm9)
	mulpd	%xmm12, %xmm10
	addpd	%xmm10, %xmm6
	MOVUPS_A2(-14 * SIZE, A1, LDA, 2, %xmm10)
	mulpd	%xmm12, %xmm11
	addpd	%xmm11, %xmm7
	MOVUPS_A2(-14 * SIZE, A1, LDA3, 1, %xmm11)

	mulpd	%xmm13, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-14 * SIZE, A2, %xmm8)
	mulpd	%xmm13, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A2(-14 * SIZE, A2, LDA, 1, %xmm9)
	mulpd	%xmm13, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A2(-14 * SIZE, A2, LDA, 2, %xmm10)
	mulpd	%xmm13, %xmm11
	addpd	%xmm11, %xmm3
	MOVUPS_A2(-14 * SIZE, A2, LDA3, 1, %xmm11)

	mulpd	%xmm13, %xmm8
	addpd	%xmm8, %xmm4
	mulpd	%xmm13, %xmm9
	addpd	%xmm9, %xmm5
	mulpd	%xmm13, %xmm10
	addpd	%xmm10, %xmm6
	mulpd	%xmm13, %xmm11
	addpd	%xmm11, %xmm7

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, A2
	addq	$4 * SIZE, X1
	ALIGN_4

L(16):
	testq	$2, M
	jle	L(17)

	MOVUPS_A1(-16 * SIZE, A1, %xmm8)
	MOVUPS_A2(-16 * SIZE, A1, LDA, 1, %xmm9)
	MOVUPS_A2(-16 * SIZE, A1, LDA, 2, %xmm10)
	MOVUPS_A2(-16 * SIZE, A1, LDA3, 1, %xmm11)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)

	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-16 * SIZE, A2, %xmm8)
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A2(-16 * SIZE, A2, LDA, 1, %xmm9)
	mulpd	%xmm12, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A2(-16 * SIZE, A2, LDA, 2, %xmm10)
	mulpd	%xmm12, %xmm11
	addpd	%xmm11, %xmm3
	MOVUPS_A2(-16 * SIZE, A2, LDA3, 1, %xmm11)

	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm4
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm5
	mulpd	%xmm12, %xmm10
	addpd	%xmm10, %xmm6
	mulpd	%xmm12, %xmm11
	addpd	%xmm11, %xmm7

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, A2
	addq	$2 * SIZE, X1
	ALIGN_4

L(17):
	testq	$1, M
	je	L(18)

	movsd	-16 * SIZE(X1), %xmm12

	movsd	-16 * SIZE(A1), %xmm8
	mulsd	%xmm12, %xmm8
	addsd	%xmm8, %xmm0
	movsd	-16 * SIZE(A1, LDA), %xmm9
	mulsd	%xmm12, %xmm9
	addsd	%xmm9, %xmm1
	movsd	-16 * SIZE(A1, LDA, 2), %xmm10
	mulsd	%xmm12, %xmm10
	addsd	%xmm10, %xmm2
	movsd	-16 * SIZE(A1, LDA3), %xmm11
	mulsd	%xmm12, %xmm11
	addsd	%xmm11, %xmm3
	movsd	-16 * SIZE(A2), %xmm8
	mulsd	%xmm12, %xmm8
	addsd	%xmm8, %xmm4
	movsd	-16 * SIZE(A2, LDA), %xmm9
	mulsd	%xmm12, %xmm9
	addsd	%xmm9, %xmm5
	movsd	-16 * SIZE(A2, LDA, 2), %xmm10
	mulsd	%xmm12, %xmm10
	addsd	%xmm10, %xmm6
	movsd	-16 * SIZE(A2, LDA3), %xmm11
	mulsd	%xmm12, %xmm11
	addsd	%xmm11, %xmm7
	ALIGN_4

L(18):
#ifdef HAVE_SSE3
	haddpd	%xmm1, %xmm0
	haddpd	%xmm3, %xmm2
	haddpd	%xmm5, %xmm4
	haddpd	%xmm7, %xmm6
#else
	movapd	%xmm0, %xmm8
	unpcklpd %xmm1, %xmm0
	unpckhpd %xmm1, %xmm8

	movapd	%xmm2, %xmm9
	unpcklpd %xmm3, %xmm2
	unpckhpd %xmm3, %xmm9

	movapd	%xmm4, %xmm10
	unpcklpd %xmm5, %xmm4
	unpckhpd %xmm5, %xmm10

	movapd	%xmm6, %xmm11
	unpcklpd %xmm7, %xmm6
	unpckhpd %xmm7, %xmm11

	addpd	%xmm8,  %xmm0
	addpd	%xmm9,  %xmm2
	addpd	%xmm10, %xmm4
	addpd	%xmm11, %xmm6
#endif

	mulpd	ALPHA, %xmm0
	mulpd	ALPHA, %xmm2
	mulpd	ALPHA, %xmm4
	mulpd	ALPHA, %xmm6

	cmpq	$SIZE, INCY
	jne	L(19)

	movsd	 0 * SIZE(Y), %xmm8
	movhpd	 1 * SIZE(Y), %xmm8
	movsd	 2 * SIZE(Y), %xmm9
	movhpd	 3 * SIZE(Y), %xmm9
	movsd	 4 * SIZE(Y), %xmm10
	movhpd	 5 * SIZE(Y), %xmm10
	movsd	 6 * SIZE(Y), %xmm11
	movhpd	 7 * SIZE(Y), %xmm11
	addq	$8 * SIZE, Y

	addpd	%xmm8, %xmm0
	addpd	%xmm9,  %xmm2
	addpd	%xmm10, %xmm4
	addpd	%xmm11, %xmm6

	movlpd	%xmm0,  0 * SIZE(Y1)
	movhpd	%xmm0,  1 * SIZE(Y1)
	movlpd	%xmm2,  2 * SIZE(Y1)
	movhpd	%xmm2,  3 * SIZE(Y1)
	movlpd	%xmm4,  4 * SIZE(Y1)
	movhpd	%xmm4,  5 * SIZE(Y1)
	movlpd	%xmm6,  6 * SIZE(Y1)
	movhpd	%xmm6,  7 * SIZE(Y1)
	addq	$8 * SIZE, Y1

	cmpq	$8, N
	jge	L(11)
	jmp	L(20)
	ALIGN_4

L(19):
	movsd	(Y), %xmm8
	addq	INCY, Y
	movhpd	(Y), %xmm8
	addq	INCY, Y
	movsd	(Y), %xmm9
	addq	INCY, Y
	movhpd	(Y), %xmm9
	addq	INCY, Y
	movsd	(Y), %xmm10
	addq	INCY, Y
	movhpd	(Y), %xmm10
	addq	INCY, Y
	movsd	(Y), %xmm11
	addq	INCY, Y
	movhpd	(Y), %xmm11
	addq	INCY, Y

	addpd	%xmm8,  %xmm0
	addpd	%xmm9,  %xmm2
	addpd	%xmm10, %xmm4
	addpd	%xmm11, %xmm6

	movlpd	%xmm0, (Y1)
	addq	INCY, Y1
	movhpd	%xmm0, (Y1)
	addq	INCY, Y1
	movlpd	%xmm2, (Y1)
	addq	INCY, Y1
	movhpd	%xmm2, (Y1)
	addq	INCY, Y1
	movlpd	%xmm4, (Y1)
	addq	INCY, Y1
	movhpd	%xmm4, (Y1)
	addq	INCY, Y1
	movlpd	%xmm6, (Y1)
	addq	INCY, Y1
	movhpd	%xmm6, (Y1)
	addq	INCY, Y1

	cmpq	$8, N
	jge	L(11)
	ALIGN_4

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

	leaq	16 * SIZE(BUFFER), X1

	movq	A, A1
	leaq	(A1, LDA, 2), A2
	leaq	(A1, LDA, 4), A

	xorps	%xmm0, %xmm0
	xorps	%xmm1, %xmm1
	xorps	%xmm2, %xmm2
	xorps	%xmm3, %xmm3

#if (GEMV_UNROLL == 4 ) && defined(PREFETCHW)
	PREFETCHW	3 * SIZE(Y1)
#endif

#ifdef ALIGNED_ACCESS
	testq	$SIZE, A
	je	L(2X)

	movsd	-16 * SIZE(X1), %xmm12

	movsd	-16 * SIZE(A1), %xmm8
	mulsd	%xmm12, %xmm8
	addsd	%xmm8, %xmm0
	movsd	-16 * SIZE(A1, LDA), %xmm9
	mulsd	%xmm12, %xmm9
	addsd	%xmm9, %xmm1
	movsd	-16 * SIZE(A2), %xmm10
	mulsd	%xmm12, %xmm10
	addsd	%xmm10, %xmm2
	movsd	-16 * SIZE(A2, LDA), %xmm11
	mulsd	%xmm12, %xmm11
	addsd	%xmm11, %xmm3

	addq	 $SIZE, A1
	addq	 $SIZE, A2
	addq	 $SIZE, X1
	ALIGN_3

L(2X):
#endif

	movq	M,  I
	sarq	$3,  I
	jle	L(25)

	MOVUPS_A1(-16 * SIZE, A1, %xmm8)
	MOVUPS_A2(-16 * SIZE, A1, LDA, 1, %xmm9)
	MOVUPS_A1(-16 * SIZE, A2, %xmm10)
	MOVUPS_A2(-16 * SIZE, A2, LDA, 1, %xmm11)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)
	MOVUPS_XL1(-14 * SIZE, X1, %xmm13)

	decq	I
	jle	L(23)
	ALIGN_4

L(22):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1)
#endif

	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-14 * SIZE, A1, %xmm8)
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A2(-14 * SIZE, A1, LDA, 1, %xmm9)
	mulpd	%xmm12, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A1(-14 * SIZE, A2, %xmm10)
	mulpd	%xmm12, %xmm11
	MOVUPS_XL1(-12 * SIZE, X1, %xmm12)
	addpd	%xmm11, %xmm3
	MOVUPS_A2(-14 * SIZE, A2, LDA, 1, %xmm11)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1, LDA)
#endif

	mulpd	%xmm13, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-12 * SIZE, A1, %xmm8)
	mulpd	%xmm13, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A2(-12 * SIZE, A1, LDA, 1, %xmm9)
	mulpd	%xmm13, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A1(-12 * SIZE, A2, %xmm10)
	mulpd	%xmm13, %xmm11
	MOVUPS_XL1(-10 * SIZE, X1, %xmm13)
	addpd	%xmm11, %xmm3
	MOVUPS_A2(-12 * SIZE, A2, LDA, 1, %xmm11)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2)
#endif

	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-10 * SIZE, A1, %xmm8)
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A2(-10 * SIZE, A1, LDA, 1, %xmm9)
	mulpd	%xmm12, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A1(-10 * SIZE, A2, %xmm10)
	mulpd	%xmm12, %xmm11
	MOVUPS_XL1( -8 * SIZE, X1, %xmm12)
	addpd	%xmm11, %xmm3
	MOVUPS_A2(-10 * SIZE, A2, LDA, 1, %xmm11)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2, LDA)
#endif

	mulpd	%xmm13, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1( -8 * SIZE, A1, %xmm8)
	mulpd	%xmm13, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A2( -8 * SIZE, A1, LDA, 1, %xmm9)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(X1)
#endif

	mulpd	%xmm13, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A1( -8 * SIZE, A2, %xmm10)
	mulpd	%xmm13, %xmm11
	MOVUPS_XL1( -6 * SIZE, X1, %xmm13)
	addpd	%xmm11, %xmm3
	MOVUPS_A2( -8 * SIZE, A2, LDA, 1, %xmm11)

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1

	decq	I
	jg	L(22)
	ALIGN_4

L(23):
	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-14 * SIZE, A1, %xmm8)
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A2(-14 * SIZE, A1, LDA, 1, %xmm9)
	mulpd	%xmm12, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A1(-14 * SIZE, A2, %xmm10)
	mulpd	%xmm12, %xmm11
	MOVUPS_XL1(-12 * SIZE, X1, %xmm12)
	addpd	%xmm11, %xmm3
	MOVUPS_A2(-14 * SIZE, A2, LDA, 1, %xmm11)

	mulpd	%xmm13, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-12 * SIZE, A1, %xmm8)
	mulpd	%xmm13, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A2(-12 * SIZE, A1, LDA, 1, %xmm9)
	mulpd	%xmm13, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A1(-12 * SIZE, A2, %xmm10)
	mulpd	%xmm13, %xmm11
	MOVUPS_XL1(-10 * SIZE, X1, %xmm13)
	addpd	%xmm11, %xmm3
	MOVUPS_A2(-12 * SIZE, A2, LDA, 1, %xmm11)

	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-10 * SIZE, A1, %xmm8)
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A2(-10 * SIZE, A1, LDA, 1, %xmm9)
	mulpd	%xmm12, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A1(-10 * SIZE, A2, %xmm10)
	mulpd	%xmm12, %xmm11
	addpd	%xmm11, %xmm3
	MOVUPS_A2(-10 * SIZE, A2, LDA, 1, %xmm11)

	mulpd	%xmm13, %xmm8
	addpd	%xmm8, %xmm0
	mulpd	%xmm13, %xmm9
	addpd	%xmm9, %xmm1
	mulpd	%xmm13, %xmm10
	addpd	%xmm10, %xmm2
	mulpd	%xmm13, %xmm11
	addpd	%xmm11, %xmm3

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1
	ALIGN_4

L(25):
	testq	$4, M
	jle	L(26)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)

	MOVUPS_A1(-16 * SIZE, A1, %xmm8)
	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A2(-16 * SIZE, A1, LDA, 1, %xmm9)
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A1(-16 * SIZE, A2, %xmm10)
	mulpd	%xmm12, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A2(-16 * SIZE, A2, LDA, 1, %xmm11)
	mulpd	%xmm12, %xmm11
	addpd	%xmm11, %xmm3

	MOVUPS_XL1(-14 * SIZE, X1, %xmm13)

	MOVUPS_A1(-14 * SIZE, A1, %xmm8)
	mulpd	%xmm13, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A2(-14 * SIZE, A1, LDA, 1, %xmm9)
	mulpd	%xmm13, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A1(-14 * SIZE, A2, %xmm10)
	mulpd	%xmm13, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A2(-14 * SIZE, A2, LDA, 1, %xmm11)
	mulpd	%xmm13, %xmm11
	addpd	%xmm11, %xmm3

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, A2
	addq	$4 * SIZE, X1
	ALIGN_4

L(26):
	testq	$2, M
	jle	L(27)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)

	MOVUPS_A1(-16 * SIZE, A1, %xmm8)
	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A2(-16 * SIZE, A1, LDA, 1, %xmm9)
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm1
	MOVUPS_A1(-16 * SIZE, A2, %xmm10)
	mulpd	%xmm12, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A2(-16 * SIZE, A2, LDA, 1, %xmm11)
	mulpd	%xmm12, %xmm11
	addpd	%xmm11, %xmm3

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, A2
	addq	$2 * SIZE, X1
	ALIGN_4

L(27):
	testq	$1, M
	je	L(28)

	movsd	-16 * SIZE(X1), %xmm12

	movsd	-16 * SIZE(A1), %xmm8
	mulsd	%xmm12, %xmm8
	addsd	%xmm8, %xmm0
	movsd	-16 * SIZE(A1, LDA), %xmm9
	mulsd	%xmm12, %xmm9
	addsd	%xmm9, %xmm1
	movsd	-16 * SIZE(A2), %xmm10
	mulsd	%xmm12, %xmm10
	addsd	%xmm10, %xmm2
	movsd	-16 * SIZE(A2, LDA), %xmm11
	mulsd	%xmm12, %xmm11
	addsd	%xmm11, %xmm3
	ALIGN_4

L(28):
#ifdef HAVE_SSE3
	haddpd	%xmm1, %xmm0
	haddpd	%xmm3, %xmm2
#else
	movapd	%xmm0, %xmm8
	unpcklpd %xmm1, %xmm0
	unpckhpd %xmm1, %xmm8

	movapd	%xmm2, %xmm9
	unpcklpd %xmm3, %xmm2
	unpckhpd %xmm3, %xmm9

	addpd	%xmm8, %xmm0
	addpd	%xmm9, %xmm2
#endif

	mulpd	ALPHA, %xmm0
	mulpd	ALPHA, %xmm2

	cmpq	$SIZE, INCY
	jne	L(29)

	movsd	 0 * SIZE(Y), %xmm4
	movhpd	 1 * SIZE(Y), %xmm4
	movsd	 2 * SIZE(Y), %xmm5
	movhpd	 3 * SIZE(Y), %xmm5
	addq	$4 * SIZE, Y

	addpd	%xmm4, %xmm0
	addpd	%xmm5, %xmm2

	movlpd	%xmm0,  0 * SIZE(Y1)
	movhpd	%xmm0,  1 * SIZE(Y1)
	movlpd	%xmm2,  2 * SIZE(Y1)
	movhpd	%xmm2,  3 * SIZE(Y1)
	addq	$4 * SIZE, Y1

#if GEMV_UNROLL == 4
	cmpq	$4, N
	jge	L(21)
#endif
	jmp	L(30)
	ALIGN_4

L(29):
	movsd	(Y), %xmm4
	addq	INCY, Y
	movhpd	(Y), %xmm4
	addq	INCY, Y
	movsd	(Y), %xmm5
	addq	INCY, Y
	movhpd	(Y), %xmm5
	addq	INCY, Y

	addpd	%xmm4, %xmm0
	addpd	%xmm5, %xmm2

	movlpd	%xmm0, (Y1)
	addq	INCY, Y1
	movhpd	%xmm0, (Y1)
	addq	INCY, Y1
	movlpd	%xmm2, (Y1)
	addq	INCY, Y1
	movhpd	%xmm2, (Y1)
	addq	INCY, Y1

#if GEMV_UNROLL == 4
	cmpq	$4, N
	jge	L(21)
#endif
	ALIGN_4

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

	leaq	16 * SIZE(BUFFER), X1

	movq	A, A1
	leaq	(A1, LDA), A2
	leaq	(A1, LDA, 2), A

	xorps	%xmm0, %xmm0
	xorps	%xmm1, %xmm1
	xorps	%xmm2, %xmm2
	xorps	%xmm3, %xmm3

#if (GEMV_UNROLL == 2 ) && defined(PREFETCHW)
	PREFETCHW	2 * SIZE(Y1)
#endif

#ifdef ALIGNED_ACCESS
	testq	$SIZE, A
	je	L(3X)

	movsd	-16 * SIZE(X1), %xmm12

	movsd	-16 * SIZE(A1), %xmm8
	mulsd	%xmm12, %xmm8
	addsd	%xmm8, %xmm0
	movsd	-16 * SIZE(A2), %xmm9
	mulsd	%xmm12, %xmm9
	addsd	%xmm9, %xmm1

	addq	 $SIZE, A1
	addq	 $SIZE, A2
	addq	 $SIZE, X1
	ALIGN_3

L(3X):
#endif

	movq	M,  I
	sarq	$3,  I
	jle	L(35)

	MOVUPS_A1(-16 * SIZE, A1, %xmm8)
	MOVUPS_A1(-16 * SIZE, A2, %xmm9)
	MOVUPS_A1(-14 * SIZE, A1, %xmm10)
	MOVUPS_A1(-14 * SIZE, A2, %xmm11)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)
	MOVUPS_XL1(-14 * SIZE, X1, %xmm13)

	decq	I
	jle	L(33)
	ALIGN_4

L(32):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A1)
#endif

	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-12 * SIZE, A1, %xmm8)
	mulpd	%xmm12, %xmm9
	MOVUPS_XL1(-12 * SIZE, X1, %xmm12)
	addpd	%xmm9, %xmm1
	MOVUPS_A1(-12 * SIZE, A2, %xmm9)

	mulpd	%xmm13, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A1(-10 * SIZE, A1, %xmm10)
	mulpd	%xmm13, %xmm11
	MOVUPS_XL1(-10 * SIZE, X1, %xmm13)
	addpd	%xmm11, %xmm3
	MOVUPS_A1(-10 * SIZE, A2, %xmm11)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A2)
#endif

	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1( -8 * SIZE, A1, %xmm8)
	mulpd	%xmm12, %xmm9
	MOVUPS_XL1( -8 * SIZE, X1, %xmm12)
	addpd	%xmm9, %xmm1
	MOVUPS_A1( -8 * SIZE, A2, %xmm9)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(X1)
#endif

	mulpd	%xmm13, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A1( -6 * SIZE, A1, %xmm10)
	mulpd	%xmm13, %xmm11
	MOVUPS_XL1( -6 * SIZE, X1, %xmm13)
	addpd	%xmm11, %xmm3
	MOVUPS_A1( -6 * SIZE, A2, %xmm11)

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1

	decq	I
	jg	L(32)
	ALIGN_4

L(33):
	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-12 * SIZE, A1, %xmm8)
	mulpd	%xmm12, %xmm9
	MOVUPS_XL1(-12 * SIZE, X1, %xmm12)
	addpd	%xmm9, %xmm1
	MOVUPS_A1(-12 * SIZE, A2, %xmm9)

	mulpd	%xmm13, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A1(-10 * SIZE, A1, %xmm10)
	mulpd	%xmm13, %xmm11
	MOVUPS_XL1(-10 * SIZE, X1, %xmm13)
	addpd	%xmm11, %xmm3
	MOVUPS_A1(-10 * SIZE, A2, %xmm11)

	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm0
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm1

	mulpd	%xmm13, %xmm10
	addpd	%xmm10, %xmm2
	mulpd	%xmm13, %xmm11
	addpd	%xmm11, %xmm3

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1
	ALIGN_4

L(35):
	testq	$4, M
	jle	L(36)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)

	MOVUPS_A1(-16 * SIZE, A1, %xmm8)
	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-16 * SIZE, A2, %xmm9)
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm1

	MOVUPS_XL1(-14 * SIZE, X1, %xmm13)

	MOVUPS_A1(-14 * SIZE, A1, %xmm10)
	mulpd	%xmm13, %xmm10
	addpd	%xmm10, %xmm2
	MOVUPS_A1(-14 * SIZE, A2, %xmm11)
	mulpd	%xmm13, %xmm11
	addpd	%xmm11, %xmm3

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, A2
	addq	$4 * SIZE, X1
	ALIGN_4

L(36):
	testq	$2, M
	jle	L(37)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)

	MOVUPS_A1(-16 * SIZE, A1, %xmm8)
	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm0
	MOVUPS_A1(-16 * SIZE, A2, %xmm9)
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm1

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, A2
	addq	$2 * SIZE, X1
	ALIGN_4

L(37):
	testq	$1, M
	je	L(38)

	movsd	-16 * SIZE(X1), %xmm12

	movsd	-16 * SIZE(A1), %xmm8
	mulsd	%xmm12, %xmm8
	addsd	%xmm8, %xmm0
	movsd	-16 * SIZE(A2), %xmm9
	mulsd	%xmm12, %xmm9
	addsd	%xmm9, %xmm1
	ALIGN_4

L(38):
	addpd	%xmm2, %xmm0
	addpd	%xmm3, %xmm1

#ifdef HAVE_SSE3
	haddpd	%xmm1, %xmm0
#else
	movapd	%xmm0, %xmm8
	unpcklpd %xmm1, %xmm0
	unpckhpd %xmm1, %xmm8

	addpd	%xmm8, %xmm0
#endif

	mulpd	ALPHA, %xmm0

	movsd	(Y), %xmm4
	addq	INCY, Y
	movhpd	(Y), %xmm4
	addq	INCY, Y

	addpd	%xmm4, %xmm0

	movlpd	%xmm0, (Y1)
	addq	INCY, Y1
	movhpd	%xmm0, (Y1)
	addq	INCY, Y1

#if GEMV_UNROLL == 2
	cmpq	$2, N
	jge	L(31)
#endif
	ALIGN_4

L(40):
	cmpq	$1, N
	jl	L(999)

#endif

	leaq	16 * SIZE(BUFFER), X1

	movq	A, A1

	xorps	%xmm0, %xmm0
	xorps	%xmm1, %xmm1
	xorps	%xmm2, %xmm2
	xorps	%xmm3, %xmm3

#ifdef ALIGNED_ACCESS
	testq	$SIZE, A
	je	L(4X)

	movsd	-16 * SIZE(X1), %xmm12
	movsd	-16 * SIZE(A1), %xmm8
	mulsd	%xmm12, %xmm8
	addsd	%xmm8, %xmm0

	addq	 $SIZE, A1
	addq	 $SIZE, X1
	ALIGN_3

L(4X):
#endif

	movq	M,  I
	sarq	$3,  I
	jle	L(45)

	MOVUPS_A1(-16 * SIZE, A1, %xmm8)
	MOVUPS_A1(-14 * SIZE, A1, %xmm9)
	MOVUPS_A1(-12 * SIZE, A1, %xmm10)
	MOVUPS_A1(-10 * SIZE, A1, %xmm11)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)
	MOVUPS_XL1(-14 * SIZE, X1, %xmm13)

	decq	I
	jle	L(43)
	ALIGN_4

L(42):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(A1)
#endif

	mulpd	%xmm12, %xmm8
	MOVUPS_XL1(-12 * SIZE, X1, %xmm12)
	addpd	%xmm8,  %xmm0
	MOVUPS_A1( -8 * SIZE, A1, %xmm8)
	mulpd	%xmm13, %xmm9
	MOVUPS_XL1(-10 * SIZE, X1, %xmm13)
	addpd	%xmm9,  %xmm2
	MOVUPS_A1( -6 * SIZE, A1, %xmm9)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(X1)
#endif

	mulpd	%xmm12, %xmm10
	MOVUPS_XL1( -8 * SIZE, X1, %xmm12)
	addpd	%xmm10, %xmm0
	MOVUPS_A1( -4 * SIZE, A1, %xmm10)
	mulpd	%xmm13, %xmm11
	MOVUPS_XL1( -6 * SIZE, X1, %xmm13)
	addpd	%xmm11, %xmm2
	MOVUPS_A1( -2 * SIZE, A1, %xmm11)

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, X1

	decq	I
	jg	L(42)
	ALIGN_4

L(43):
	mulpd	%xmm12, %xmm8
	MOVUPS_XL1(-12 * SIZE, X1, %xmm12)
	addpd	%xmm8,  %xmm0
	mulpd	%xmm13, %xmm9
	MOVUPS_XL1(-10 * SIZE, X1, %xmm13)
	addpd	%xmm9,  %xmm2

	mulpd	%xmm12, %xmm10
	addpd	%xmm10, %xmm0
	mulpd	%xmm13, %xmm11
	addpd	%xmm11, %xmm2

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, X1
	ALIGN_4

L(45):
	testq	$4, M
	jle	L(46)

	MOVUPS_A1(-16 * SIZE, A1, %xmm8)
	MOVUPS_A1(-14 * SIZE, A1, %xmm9)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)
	MOVUPS_XL1(-14 * SIZE, X1, %xmm13)

	mulpd	%xmm12, %xmm8
	addpd	%xmm8,  %xmm0
	mulpd	%xmm13, %xmm9
	addpd	%xmm9,  %xmm2

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, X1
	ALIGN_4

L(46):
	testq	$2, M
	jle	L(47)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)

	MOVUPS_A1(-16 * SIZE, A1, %xmm8)
	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm0

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, X1
	ALIGN_4

L(47):
	testq	$1, M
	je	L(48)

	movsd	-16 * SIZE(X1), %xmm12

	movsd	-16 * SIZE(A1), %xmm8
	mulsd	%xmm12, %xmm8
	addsd	%xmm8, %xmm0
	ALIGN_4

L(48):
	addpd	%xmm2, %xmm0
	addpd	%xmm3, %xmm1

	addpd	%xmm1, %xmm0

#ifdef HAVE_SSE3
	haddpd	%xmm1, %xmm0
#else
	movapd	%xmm0, %xmm8
	unpcklpd %xmm1, %xmm0
	unpckhpd %xmm1, %xmm8

	addsd	%xmm8, %xmm0
#endif

	mulsd	ALPHA, %xmm0

	movsd	(Y), %xmm4
	addq	INCY, Y

	addsd	%xmm4, %xmm0

	movlpd	%xmm0, (Y1)
	addq	INCY, Y1

#ifdef ALIGNED_ACCESS
	jmp	L(999)
	ALIGN_4

L(50):
#if GEMV_UNROLL >= 4

	cmpq	$4, N
	jl	L(60)
	ALIGN_3

L(51):
	subq	$4, N

	leaq	16 * SIZE(BUFFER), X1

	movq	A, A1
	leaq	(A1, LDA, 2), A2
	leaq	(A1, LDA, 4), A

	xorps	%xmm0, %xmm0
	xorps	%xmm1, %xmm1
	xorps	%xmm2, %xmm2
	xorps	%xmm3, %xmm3

#ifdef PREFETCHW
	PREFETCHW	3 * SIZE(Y1)
#endif

#ifdef ALIGNED_ACCESS
	testq	$SIZE, A
	je	L(5X)

	movsd	-16 * SIZE(X1), %xmm12

	movsd	-16 * SIZE(A1), %xmm4
	mulsd	%xmm12, %xmm4
	addsd	%xmm4, %xmm0
	movsd	-16 * SIZE(A1, LDA), %xmm5
	mulsd	%xmm12, %xmm5
	addsd	%xmm5, %xmm1
	movsd	-16 * SIZE(A2), %xmm6
	mulsd	%xmm12, %xmm6
	addsd	%xmm6, %xmm2
	movsd	-16 * SIZE(A2, LDA), %xmm7
	mulsd	%xmm12, %xmm7
	addsd	%xmm7, %xmm3

	addq	 $SIZE, A1
	addq	 $SIZE, A2
	addq	 $SIZE, X1
	ALIGN_3

L(5X):
#endif

	movhpd	-16 * SIZE(A1, LDA), %xmm8
	movhpd	-16 * SIZE(A2, LDA), %xmm9

	movq	M,  I
	sarq	$3,  I
	jle	L(55)

	MOVUPS_A1(-16 * SIZE, A1, %xmm4)
	MOVUPS_A2(-15 * SIZE, A1, LDA, 1, %xmm5)
	MOVUPS_A1(-16 * SIZE, A2, %xmm6)
	MOVUPS_A2(-15 * SIZE, A2, LDA, 1, %xmm7)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)
	MOVUPS_XL1(-14 * SIZE, X1, %xmm13)

	decq	I
	jle	L(53)
	ALIGN_4

L(52):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A1)
#endif

	mulpd	%xmm12, %xmm4
	addpd	%xmm4, %xmm0
	MOVUPS_A1(-14 * SIZE, A1, %xmm4)
	shufpd	$1, %xmm5, %xmm8
	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm1
	MOVUPS_A2(-13 * SIZE, A1, LDA, 1, %xmm8)

	mulpd	%xmm12, %xmm6
	addpd	%xmm6, %xmm2
	MOVUPS_A1(-14 * SIZE, A2, %xmm6)
	shufpd	$1, %xmm7, %xmm9
	mulpd	%xmm12, %xmm9
	MOVUPS_XL1(-12 * SIZE, X1, %xmm12)
	addpd	%xmm9, %xmm3
	MOVUPS_A2(-13 * SIZE, A2, LDA, 1, %xmm9)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET + 8(A1, LDA)
#endif
	mulpd	%xmm13, %xmm4
	addpd	%xmm4, %xmm0
	MOVUPS_A1(-12 * SIZE, A1, %xmm4)
	shufpd	$1, %xmm8, %xmm5
	mulpd	%xmm13, %xmm5
	addpd	%xmm5, %xmm1
	MOVUPS_A2(-11 * SIZE, A1, LDA, 1, %xmm5)

	mulpd	%xmm13, %xmm6
	addpd	%xmm6, %xmm2
	MOVUPS_A1(-12 * SIZE, A2, %xmm6)
	shufpd	$1, %xmm9, %xmm7
	mulpd	%xmm13, %xmm7
	MOVUPS_XL1(-10 * SIZE, X1, %xmm13)
	addpd	%xmm7, %xmm3
	MOVUPS_A2(-11 * SIZE, A2, LDA, 1, %xmm7)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET(A2)
#endif

	mulpd	%xmm12, %xmm4
	addpd	%xmm4, %xmm0
	MOVUPS_A1(-10 * SIZE, A1, %xmm4)
	shufpd	$1, %xmm5, %xmm8
	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm1
	MOVUPS_A2( -9 * SIZE, A1, LDA, 1, %xmm8)

	mulpd	%xmm12, %xmm6
	addpd	%xmm6, %xmm2
	MOVUPS_A1(-10 * SIZE, A2, %xmm6)
	shufpd	$1, %xmm7, %xmm9
	mulpd	%xmm12, %xmm9
	MOVUPS_XL1(-8 * SIZE, X1, %xmm12)
	addpd	%xmm9, %xmm3
	MOVUPS_A2( -9 * SIZE, A2, LDA, 1, %xmm9)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET + 8(A2, LDA)
#endif

	mulpd	%xmm13, %xmm4
	addpd	%xmm4, %xmm0
	MOVUPS_A1(-8 * SIZE, A1, %xmm4)
	shufpd	$1, %xmm8, %xmm5
	mulpd	%xmm13, %xmm5
	addpd	%xmm5, %xmm1
	MOVUPS_A2(-7 * SIZE, A1, LDA, 1, %xmm5)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 2 - 128 + PREOFFSET + 8(X1)
#endif

	mulpd	%xmm13, %xmm6
	addpd	%xmm6, %xmm2
	MOVUPS_A1(-8 * SIZE, A2, %xmm6)
	shufpd	$1, %xmm9, %xmm7
	mulpd	%xmm13, %xmm7
	MOVUPS_XL1(-6 * SIZE, X1, %xmm13)
	addpd	%xmm7, %xmm3
	MOVUPS_A2(-7 * SIZE, A2, LDA, 1, %xmm7)

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1

	decq	I
	jg	L(52)
	ALIGN_4

L(53):
	mulpd	%xmm12, %xmm4
	addpd	%xmm4, %xmm0
	MOVUPS_A1(-14 * SIZE, A1, %xmm4)
	shufpd	$1, %xmm5, %xmm8
	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm1
	MOVUPS_A2(-13 * SIZE, A1, LDA, 1, %xmm8)

	mulpd	%xmm12, %xmm6
	addpd	%xmm6, %xmm2
	MOVUPS_A1(-14 * SIZE, A2, %xmm6)
	shufpd	$1, %xmm7, %xmm9
	mulpd	%xmm12, %xmm9
	MOVUPS_XL1(-12 * SIZE, X1, %xmm12)
	addpd	%xmm9, %xmm3
	MOVUPS_A2(-13 * SIZE, A2, LDA, 1, %xmm9)

	mulpd	%xmm13, %xmm4
	addpd	%xmm4, %xmm0
	MOVUPS_A1(-12 * SIZE, A1, %xmm4)
	shufpd	$1, %xmm8, %xmm5
	mulpd	%xmm13, %xmm5
	addpd	%xmm5, %xmm1
	MOVUPS_A2(-11 * SIZE, A1, LDA, 1, %xmm5)
	mulpd	%xmm13, %xmm6
	addpd	%xmm6, %xmm2
	MOVUPS_A1(-12 * SIZE, A2, %xmm6)
	shufpd	$1, %xmm9, %xmm7
	mulpd	%xmm13, %xmm7
	MOVUPS_XL1(-10 * SIZE, X1, %xmm13)
	addpd	%xmm7, %xmm3
	MOVUPS_A2(-11 * SIZE, A2, LDA, 1, %xmm7)

	mulpd	%xmm12, %xmm4
	addpd	%xmm4, %xmm0
	MOVUPS_A1(-10 * SIZE, A1, %xmm4)
	shufpd	$1, %xmm5, %xmm8
	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm1
	MOVUPS_A2( -9 * SIZE, A1, LDA, 1, %xmm8)

	mulpd	%xmm12, %xmm6
	addpd	%xmm6, %xmm2
	MOVUPS_A1(-10 * SIZE, A2, %xmm6)
	shufpd	$1, %xmm7, %xmm9
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm3
	MOVUPS_A2( -9 * SIZE, A2, LDA, 1, %xmm9)

	mulpd	%xmm13, %xmm4
	addpd	%xmm4, %xmm0
	shufpd	$1, %xmm8, %xmm5
	mulpd	%xmm13, %xmm5
	addpd	%xmm5, %xmm1

	mulpd	%xmm13, %xmm6
	addpd	%xmm6, %xmm2
	shufpd	$1, %xmm9, %xmm7
	mulpd	%xmm13, %xmm7
	addpd	%xmm7, %xmm3

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1
	ALIGN_4

L(55):
	testq	$4, M
	jle	L(56)

	MOVUPS_A1(-16 * SIZE, A1, %xmm4)
	MOVUPS_A2(-15 * SIZE, A1, LDA, 1, %xmm5)
	MOVUPS_A1(-16 * SIZE, A2, %xmm6)
	MOVUPS_A2(-15 * SIZE, A2, LDA, 1, %xmm7)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)
	MOVUPS_XL1(-14 * SIZE, X1, %xmm13)

	mulpd	%xmm12, %xmm4
	addpd	%xmm4, %xmm0
	MOVUPS_A1(-14 * SIZE, A1, %xmm4)

	shufpd	$1, %xmm5, %xmm8
	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm1
	MOVUPS_A2(-13 * SIZE, A1, LDA, 1, %xmm8)

	mulpd	%xmm12, %xmm6
	addpd	%xmm6, %xmm2
	MOVUPS_A1(-14 * SIZE, A2, %xmm6)

	shufpd	$1, %xmm7, %xmm9
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm3
	MOVUPS_A2(-13 * SIZE, A2, LDA, 1, %xmm9)

	mulpd	%xmm13, %xmm4
	addpd	%xmm4, %xmm0
	shufpd	$1, %xmm8, %xmm5
	mulpd	%xmm13, %xmm5
	addpd	%xmm5, %xmm1
	mulpd	%xmm13, %xmm6
	addpd	%xmm6, %xmm2
	shufpd	$1, %xmm9, %xmm7
	mulpd	%xmm13, %xmm7
	addpd	%xmm7, %xmm3

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, A2
	addq	$4 * SIZE, X1
	ALIGN_4

L(56):
	testq	$2, M
	jle	L(57)

	MOVUPS_A1(-16 * SIZE, A1, %xmm4)
	MOVUPS_A2(-15 * SIZE, A1, LDA, 1, %xmm5)
	MOVUPS_A1(-16 * SIZE, A2, %xmm6)
	MOVUPS_A2(-15 * SIZE, A2, LDA, 1, %xmm7)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)

	mulpd	%xmm12, %xmm4
	addpd	%xmm4, %xmm0
	shufpd	$1, %xmm5, %xmm8
	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm1
	movaps	%xmm5, %xmm8
	mulpd	%xmm12, %xmm6
	addpd	%xmm6, %xmm2
	shufpd	$1, %xmm7, %xmm9
	mulpd	%xmm12, %xmm9
	addpd	%xmm9, %xmm3
	movaps	%xmm7, %xmm9

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, A2
	addq	$2 * SIZE, X1
	ALIGN_4

L(57):
	testq	$1, M
	je	L(58)

	movsd	-16 * SIZE(X1), %xmm12

	movsd	-16 * SIZE(A1), %xmm4
	mulsd	%xmm12, %xmm4
	addsd	%xmm4, %xmm0
	shufpd	$1, %xmm8, %xmm8
	mulsd	%xmm12, %xmm8
	addsd	%xmm8, %xmm1
	movsd	-16 * SIZE(A2), %xmm6
	mulsd	%xmm12, %xmm6
	addsd	%xmm6, %xmm2
	shufpd	$1, %xmm9, %xmm9
	mulsd	%xmm12, %xmm9
	addsd	%xmm9, %xmm3
	ALIGN_4

L(58):
#ifdef HAVE_SSE3
	haddpd	%xmm1, %xmm0
	haddpd	%xmm3, %xmm2
#else
	movapd	%xmm0, %xmm4
	unpcklpd %xmm1, %xmm0
	unpckhpd %xmm1, %xmm4

	movapd	%xmm2, %xmm5
	unpcklpd %xmm3, %xmm2
	unpckhpd %xmm3, %xmm5

	addpd	%xmm4, %xmm0
	addpd	%xmm5, %xmm2
#endif

	mulpd	ALPHA, %xmm0
	mulpd	ALPHA, %xmm2

	cmpq	$SIZE, INCY
	jne	L(59)

	movsd	 0 * SIZE(Y), %xmm4
	movhpd	 1 * SIZE(Y), %xmm4
	movsd	 2 * SIZE(Y), %xmm5
	movhpd	 3 * SIZE(Y), %xmm5
	addq	$4 * SIZE, Y

	addpd	%xmm4, %xmm0
	addpd	%xmm5, %xmm2

	movlpd	%xmm0,  0 * SIZE(Y1)
	movhpd	%xmm0,  1 * SIZE(Y1)
	movlpd	%xmm2,  2 * SIZE(Y1)
	movhpd	%xmm2,  3 * SIZE(Y1)
	addq	$4 * SIZE, Y1

	cmpq	$4, N
	jge	L(51)
	jmp	L(60)
	ALIGN_4

L(59):
	movsd	(Y), %xmm4
	addq	INCY, Y
	movhpd	(Y), %xmm4
	addq	INCY, Y
	movsd	(Y), %xmm5
	addq	INCY, Y
	movhpd	(Y), %xmm5
	addq	INCY, Y

	addpd	%xmm4, %xmm0
	addpd	%xmm5, %xmm2

	movlpd	%xmm0, (Y1)
	addq	INCY, Y1
	movhpd	%xmm0, (Y1)
	addq	INCY, Y1
	movlpd	%xmm2, (Y1)
	addq	INCY, Y1
	movhpd	%xmm2, (Y1)
	addq	INCY, Y1
	cmpq	$4, N
	jge	L(51)
	ALIGN_4

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

	leaq	16 * SIZE(BUFFER), X1

	movq	A, A1
	leaq	(A1, LDA), A2
	leaq	(A1, LDA, 2), A

	xorps	%xmm0, %xmm0
	xorps	%xmm1, %xmm1
	xorps	%xmm2, %xmm2
	xorps	%xmm3, %xmm3

#if (GEMV_UNROLL == 2 ) && defined(PREFETCHW)
	PREFETCHW	2 * SIZE(Y1)
#endif

#ifdef ALIGNED_ACCESS
	testq	$SIZE, A
	je	L(6X)

	movsd	-16 * SIZE(X1), %xmm12

	movsd	-16 * SIZE(A1), %xmm4
	mulsd	%xmm12, %xmm4
	addsd	%xmm4, %xmm0
	movsd	-16 * SIZE(A2), %xmm5
	mulsd	%xmm12, %xmm5
	addsd	%xmm5, %xmm1

	addq	 $SIZE, A1
	addq	 $SIZE, A2
	addq	 $SIZE, X1
	ALIGN_3

L(6X):
#endif

	movhpd	-16 * SIZE(A2), %xmm8

	movq	M,  I
	sarq	$3,  I
	jle	L(65)

	MOVUPS_A1(-16 * SIZE, A1, %xmm4)
	MOVUPS_A1(-15 * SIZE, A2, %xmm5)
	MOVUPS_A1(-14 * SIZE, A1, %xmm6)
	MOVUPS_A1(-13 * SIZE, A2, %xmm7)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)
	MOVUPS_XL1(-14 * SIZE, X1, %xmm13)

	decq	I
	jle	L(63)
	ALIGN_4

L(62):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET(A1)
#endif

	mulpd	%xmm12, %xmm4
	addpd	%xmm4, %xmm0
	MOVUPS_A1(-12 * SIZE, A1, %xmm4)
	shufpd	$1, %xmm5, %xmm8
	mulpd	%xmm12, %xmm8
	MOVUPS_XL1(-12 * SIZE, X1, %xmm12)
	addpd	%xmm8, %xmm1
	MOVUPS_A1(-11 * SIZE, A2, %xmm9)

	mulpd	%xmm13, %xmm6
	addpd	%xmm6, %xmm0
	MOVUPS_A1(-10 * SIZE, A1, %xmm6)
	shufpd	$1, %xmm7, %xmm5
	mulpd	%xmm13, %xmm5
	MOVUPS_XL1(-10 * SIZE, X1, %xmm13)
	addpd	%xmm5, %xmm1
	MOVUPS_A1( -9 * SIZE, A2, %xmm8)

#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET + 8(A2)
#endif

	mulpd	%xmm12, %xmm4
	addpd	%xmm4, %xmm0
	MOVUPS_A1(-8 * SIZE, A1, %xmm4)
	shufpd	$1, %xmm9, %xmm7
	mulpd	%xmm12, %xmm7
	MOVUPS_XL1(-8 * SIZE, X1, %xmm12)
	addpd	%xmm7, %xmm1
	MOVUPS_A1(-7 * SIZE, A2, %xmm5)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 4 - 128 + PREOFFSET + 8(X1)
#endif

	mulpd	%xmm13, %xmm6
	addpd	%xmm6, %xmm0
	MOVUPS_A1(-6 * SIZE, A1, %xmm6)
	shufpd	$1, %xmm8, %xmm9
	mulpd	%xmm13, %xmm9
	MOVUPS_XL1(-6 * SIZE, X1, %xmm13)
	addpd	%xmm9, %xmm1
	MOVUPS_A1(-5 * SIZE, A2, %xmm7)

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1

	decq	I
	jg	L(62)
	ALIGN_4

L(63):
	mulpd	%xmm12, %xmm4
	addpd	%xmm4, %xmm0
	MOVUPS_A1(-12 * SIZE, A1, %xmm4)
	shufpd	$1, %xmm5, %xmm8
	mulpd	%xmm12, %xmm8
	MOVUPS_XL1(-12 * SIZE, X1, %xmm12)
	addpd	%xmm8, %xmm1
	MOVUPS_A1(-11 * SIZE, A2, %xmm9)

	mulpd	%xmm13, %xmm6
	addpd	%xmm6, %xmm0
	MOVUPS_A1(-10 * SIZE, A1, %xmm6)
	shufpd	$1, %xmm7, %xmm5
	mulpd	%xmm13, %xmm5
	MOVUPS_XL1(-10 * SIZE, X1, %xmm13)
	addpd	%xmm5, %xmm1
	MOVUPS_A1( -9 * SIZE, A2, %xmm8)

	mulpd	%xmm12, %xmm4
	addpd	%xmm4, %xmm0
	shufpd	$1, %xmm9, %xmm7
	mulpd	%xmm12, %xmm7
	addpd	%xmm7, %xmm1

	mulpd	%xmm13, %xmm6
	addpd	%xmm6, %xmm0
	shufpd	$1, %xmm8, %xmm9
	mulpd	%xmm13, %xmm9
	addpd	%xmm9, %xmm1

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, A2
	addq	$8 * SIZE, X1
	ALIGN_4

L(65):
	testq	$4, M
	jle	L(66)

	MOVUPS_A1(-16 * SIZE, A1, %xmm4)
	MOVUPS_A1(-15 * SIZE, A2, %xmm5)
	MOVUPS_A1(-14 * SIZE, A1, %xmm6)
	MOVUPS_A1(-13 * SIZE, A2, %xmm7)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)
	MOVUPS_XL1(-14 * SIZE, X1, %xmm13)

	mulpd	%xmm12, %xmm4
	addpd	%xmm4, %xmm0
	shufpd	$1, %xmm5, %xmm8
	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm1

	mulpd	%xmm13, %xmm6
	addpd	%xmm6, %xmm0
	shufpd	$1, %xmm7, %xmm5
	movaps	%xmm7, %xmm8
	mulpd	%xmm13, %xmm5
	addpd	%xmm5, %xmm1

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, A2
	addq	$4 * SIZE, X1
	ALIGN_4

L(66):
	testq	$2, M
	jle	L(67)

	MOVUPS_A1(-16 * SIZE, A1, %xmm4)
	MOVUPS_A1(-15 * SIZE, A2, %xmm5)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)

	mulpd	%xmm12, %xmm4
	addpd	%xmm4, %xmm0
	shufpd	$1, %xmm5, %xmm8
	mulpd	%xmm12, %xmm8
	addpd	%xmm8, %xmm1
	movaps	%xmm5, %xmm8

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, A2
	addq	$2 * SIZE, X1
	ALIGN_4

L(67):
	testq	$1, M
	je	L(68)

	movsd	-16 * SIZE(X1), %xmm12

	movsd	-16 * SIZE(A1), %xmm4
	mulsd	%xmm12, %xmm4
	addsd	%xmm4, %xmm0
	shufpd	$1, %xmm8, %xmm8
	mulsd	%xmm12, %xmm8
	addsd	%xmm8, %xmm1
	ALIGN_4

L(68):
	addpd	%xmm2, %xmm0
	addpd	%xmm3, %xmm1

#ifdef HAVE_SSE3
	haddpd	%xmm1, %xmm0
#else
	movapd	%xmm0, %xmm4
	unpcklpd %xmm1, %xmm0
	unpckhpd %xmm1, %xmm4

	addpd	%xmm4, %xmm0
#endif

	mulpd	ALPHA, %xmm0

	movsd	(Y), %xmm4
	addq	INCY, Y
	movhpd	(Y), %xmm4
	addq	INCY, Y

	addpd	%xmm4, %xmm0

	movlpd	%xmm0, (Y1)
	addq	INCY, Y1
	movhpd	%xmm0, (Y1)
	addq	INCY, Y1

#if GEMV_UNROLL == 2
	cmpq	$2, N
	jge	L(61)
#endif
	ALIGN_4

L(70):
	cmpq	$1, N
	jl	L(999)

#endif

	leaq	16 * SIZE(BUFFER), X1

	movq	A, A1

	xorps	%xmm0, %xmm0
	xorps	%xmm1, %xmm1
	xorps	%xmm2, %xmm2
	xorps	%xmm3, %xmm3

#ifdef ALIGNED_ACCESS
	testq	$SIZE, A
	je	L(7X)

	movsd	-16 * SIZE(X1), %xmm12

	movsd	-16 * SIZE(A1), %xmm4
	mulsd	%xmm12, %xmm4
	addsd	%xmm4, %xmm0

	addq	 $SIZE, A1
	addq	 $SIZE, X1
	ALIGN_3

L(7X):
#endif
	movq	M,  I
	sarq	$3,  I
	jle	L(75)

	MOVUPS_A1(-16 * SIZE, A1, %xmm4)
	MOVUPS_A1(-14 * SIZE, A1, %xmm5)
	MOVUPS_A1(-12 * SIZE, A1, %xmm6)
	MOVUPS_A1(-10 * SIZE, A1, %xmm7)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)
	MOVUPS_XL1(-14 * SIZE, X1, %xmm13)

	decq	I
	jle	L(73)
	ALIGN_4

L(72):
#ifdef PREFETCH
	PREFETCH	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(A1)
#endif

	mulpd	%xmm12, %xmm4
	MOVUPS_XL1(-12 * SIZE, X1, %xmm12)
	addpd	%xmm4,  %xmm0
	MOVUPS_A1( -8 * SIZE, A1, %xmm4)
	mulpd	%xmm13, %xmm5
	MOVUPS_XL1(-10 * SIZE, X1, %xmm13)
	addpd	%xmm5,  %xmm2
	MOVUPS_A1( -6 * SIZE, A1, %xmm5)

#ifdef PREFETCHW
	PREFETCH	(PREFETCHSIZE) * 8 - 128 + PREOFFSET(X1)
#endif

	mulpd	%xmm12, %xmm6
	MOVUPS_XL1( -8 * SIZE, X1, %xmm12)
	addpd	%xmm6, %xmm0
	MOVUPS_A1( -4 * SIZE, A1, %xmm6)
	mulpd	%xmm13, %xmm7
	MOVUPS_XL1( -6 * SIZE, X1, %xmm13)
	addpd	%xmm7, %xmm2
	MOVUPS_A1( -2 * SIZE, A1, %xmm7)

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, X1

	decq	I
	jg	L(72)
	ALIGN_4

L(73):
	mulpd	%xmm12, %xmm4
	MOVUPS_XL1(-12 * SIZE, X1, %xmm12)
	addpd	%xmm4,  %xmm0
	mulpd	%xmm13, %xmm5
	MOVUPS_XL1(-10 * SIZE, X1, %xmm13)
	addpd	%xmm5,  %xmm2

	mulpd	%xmm12, %xmm6
	addpd	%xmm6, %xmm0
	mulpd	%xmm13, %xmm7
	addpd	%xmm7, %xmm2

	addq	$8 * SIZE, A1
	addq	$8 * SIZE, X1
	ALIGN_4

L(75):
	testq	$4, M
	jle	L(76)

	MOVUPS_A1(-16 * SIZE, A1, %xmm4)
	MOVUPS_A1(-14 * SIZE, A1, %xmm5)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)
	MOVUPS_XL1(-14 * SIZE, X1, %xmm13)

	mulpd	%xmm12, %xmm4
	addpd	%xmm4,  %xmm0
	mulpd	%xmm13, %xmm5
	addpd	%xmm5,  %xmm2

	addq	$4 * SIZE, A1
	addq	$4 * SIZE, X1
	ALIGN_4

L(76):
	testq	$2, M
	jle	L(77)

	MOVUPS_XL1(-16 * SIZE, X1, %xmm12)

	MOVUPS_A1(-16 * SIZE, A1, %xmm4)
	mulpd	%xmm12, %xmm4
	addpd	%xmm4, %xmm0

	addq	$2 * SIZE, A1
	addq	$2 * SIZE, X1
	ALIGN_4

L(77):
	testq	$1, M
	je	L(78)

	movsd	-16 * SIZE(X1), %xmm12

	movsd	-16 * SIZE(A1), %xmm4
	mulsd	%xmm12, %xmm4
	addsd	%xmm4, %xmm0
	ALIGN_4

L(78):
	addpd	%xmm2, %xmm0
	addpd	%xmm3, %xmm1

	addpd	%xmm1, %xmm0

#ifdef HAVE_SSE3
	haddpd	%xmm1, %xmm0
#else
	movapd	%xmm0, %xmm4
	unpcklpd %xmm1, %xmm0
	unpckhpd %xmm1, %xmm4

	addsd	%xmm4, %xmm0
#endif

	mulsd	ALPHA, %xmm0

	movsd	(Y), %xmm4
	addq	INCY, Y

	addsd	%xmm4, %xmm0

	movlpd	%xmm0, (Y1)
	addq	INCY, Y1
#endif
	ALIGN_4

L(999):
	leaq	(, M, SIZE), %rax
	addq %rax,AA
	jmp L(0x);
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
