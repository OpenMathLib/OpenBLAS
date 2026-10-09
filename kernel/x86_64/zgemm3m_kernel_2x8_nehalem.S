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
#define OLD_K	%rdx

#define M	%r13
#define N	%r14
#define K	%r15

#define A	%rcx
#define B	%r8
#define C	%r9
#define LDC	%r10

#define I	%r11
#define AO	%rdi
#define BO	%rsi
#define	CO1	%rbx
#define CO2	%rbp
#define BB	%r12

#define BX	%rdx

#ifndef WINDOWS_ABI

#define STACKSIZE 128

#define OLD_LDC		 8 + STACKSIZE(%rsp)
#define OLD_OFFSET	16 + STACKSIZE(%rsp)

#define ALPHA_R	   48(%rsp)
#define ALPHA_I	   56(%rsp)
#define J	   64(%rsp)
#define OFFSET	   72(%rsp)
#define KK	   80(%rsp)
#define KKK	   88(%rsp)

#else

#define STACKSIZE 512

#define OLD_ALPHA_I	40 + STACKSIZE(%rsp)
#define OLD_A		48 + STACKSIZE(%rsp)
#define OLD_B		56 + STACKSIZE(%rsp)
#define OLD_C		64 + STACKSIZE(%rsp)
#define OLD_LDC		72 + STACKSIZE(%rsp)
#define OLD_OFFSET	80 + STACKSIZE(%rsp)

#define ALPHA_R	  224(%rsp)
#define ALPHA_I	  232(%rsp)
#define J	  240(%rsp)
#define OFFSET	  248(%rsp)
#define KK	  256(%rsp)
#define KKK	  264(%rsp)

#endif

#define PREFETCHSIZE  (8 * 1 - 4)
#define PREFETCH     prefetcht0

#if defined(OS_LINUX) && defined(CORE_BARCELONA)
	.align 32768
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

	movq	ARG1,      OLD_M
	movq	ARG2,      OLD_N
	movq	ARG3,      OLD_K
	movq	OLD_A,     A
	movq	OLD_B,     B
	movq	OLD_C,     C
	movq	OLD_LDC,   LDC
#ifdef TRMMKERNEL
	movq	OLD_OFFSET, %r11
#endif
	movaps	%xmm3, %xmm0
	movsd	OLD_ALPHA_I, %xmm1
#else
	movq	OLD_LDC,   LDC
#ifdef TRMMKERNEL
	movq	OLD_OFFSET, %r11
#endif
#endif

	movlps	 %xmm0, ALPHA_R
	movlps	 %xmm1, ALPHA_I

	subq	$-16 * SIZE, A
	subq	$-16 * SIZE, B

	movq	OLD_M, M
	movq	OLD_N, N
	movq	OLD_K, K

	salq	$ZBASE_SHIFT, LDC

#ifdef TRMMKERNEL
	movq	%r11, OFFSET
#ifndef LEFT
	negq	%r11
#endif
	movq	%r11, KK
#endif

	movq	N,  J
	sarq	$3, J
	NOBRANCH
	jle	L(30)
	ALIGN_4

L(01):
#if defined(TRMMKERNEL) && defined(LEFT)
        movq    OFFSET, %rax
	movq    %rax, KK
#endif

	movq	C, CO1
	leaq	(C, LDC, 4), CO2
	movq	A, AO

	movq	K, %rax
	salq	$BASE_SHIFT + 3, %rax
	leaq	(B, %rax), BB

	movq	M,  I
	sarq	$1, I
	NOBRANCH
	jle	L(20)
	ALIGN_4

L(11):
#if !defined(TRMMKERNEL) || \
	(defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
	(defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))

	movq	B, BO
#else
	movq	B, BO

	movq	KK, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 2), AO
	leaq	(BO, %rax, 8), BO
#endif

	prefetcht0	 -16 * SIZE(BB)
	subq		 $-8 * SIZE, BB

	xorps	%xmm1, %xmm1
	movaps	-16 * SIZE(AO), %xmm0
	xorps	%xmm2, %xmm2
	xorps	%xmm3, %xmm3
	xorps	%xmm4, %xmm4

	leaq	(LDC, LDC, 2), %rax

	xorps	%xmm8,  %xmm8
	prefetcht0     3 * SIZE(CO1)
	xorps	%xmm9,  %xmm9
	prefetcht0     7 * SIZE(CO1, LDC,  1)
	xorps	%xmm10, %xmm10
	prefetcht0     3 * SIZE(CO1, LDC,  2)
	xorps	%xmm11, %xmm11
	prefetcht0     7 * SIZE(CO1, %rax, 1)

	xorps	%xmm12, %xmm12
	prefetcht0     3 * SIZE(CO2)
	xorps	%xmm13, %xmm13
	prefetcht0     7 * SIZE(CO2, LDC,  1)
	xorps	%xmm14, %xmm14
	prefetcht0     3 * SIZE(CO2, LDC,  2)
	xorps	%xmm15, %xmm15
	prefetcht0     7 * SIZE(CO2, %rax, 1)

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
	addq	$8, %rax
#endif
	movq	%rax, KKK
#endif
	sarq	$2, %rax
	NOBRANCH
	jle	L(15)
	ALIGN_3

L(12):
	PREFETCH	(PREFETCHSIZE +  0) * SIZE(AO)

	addpd	%xmm1, %xmm12
	movaps	-16 * SIZE(BO), %xmm6
	addpd	%xmm2, %xmm13
	pshufd	$0x4e, %xmm6, %xmm2
	mulpd	%xmm0, %xmm6
	mulpd	%xmm0, %xmm2

	addpd	%xmm3, %xmm14
	movaps	-14 * SIZE(BO), %xmm3
	addpd	%xmm4, %xmm15
	pshufd	$0x4e, %xmm3, %xmm4
	mulpd	%xmm0, %xmm3
	mulpd	%xmm0, %xmm4

	addpd	%xmm6, %xmm8
	movaps	-12 * SIZE(BO), %xmm6
	addpd	%xmm2, %xmm9
	pshufd	$0x4e, %xmm6, %xmm2
	mulpd	%xmm0, %xmm6
	mulpd	%xmm0, %xmm2

	addpd	%xmm3, %xmm10
	movaps	-10 * SIZE(BO), %xmm3
	addpd	%xmm4, %xmm11
	pshufd	$0x4e, %xmm3, %xmm4
	mulpd	%xmm0, %xmm3
	mulpd	%xmm0, %xmm4

	addpd	%xmm6, %xmm12
	movaps	 -8 * SIZE(BO), %xmm1
	addpd	%xmm2, %xmm13
	movaps	-14 * SIZE(AO), %xmm5
	pshufd	$0x4e, %xmm1, %xmm2
	mulpd	%xmm5, %xmm1
	mulpd	%xmm5, %xmm2

	addpd	%xmm3, %xmm14
	movaps	 -6 * SIZE(BO), %xmm3
	addpd	%xmm4, %xmm15
	pshufd	$0x4e, %xmm3, %xmm4
	mulpd	%xmm5, %xmm3
	mulpd	%xmm5, %xmm4

	addpd	%xmm1, %xmm8
	movaps	 -4 * SIZE(BO), %xmm1
	addpd	%xmm2, %xmm9
	pshufd	$0x4e, %xmm1, %xmm2
	mulpd	%xmm5, %xmm1
	mulpd	%xmm5, %xmm2

	addpd	%xmm3, %xmm10
	movaps	 -2 * SIZE(BO), %xmm3
	addpd	%xmm4, %xmm11
	pshufd	$0x4e, %xmm3, %xmm4
	movaps	-12 * SIZE(AO), %xmm0
	mulpd	%xmm5, %xmm3
	mulpd	%xmm5, %xmm4

	addpd	%xmm1, %xmm12
	movaps	  0 * SIZE(BO), %xmm6
	addpd	%xmm2, %xmm13
	pshufd	$0x4e, %xmm6, %xmm2
	mulpd	%xmm0, %xmm6
	mulpd	%xmm0, %xmm2

	addpd	%xmm3, %xmm14
	movaps	  2 * SIZE(BO), %xmm3
	addpd	%xmm4, %xmm15
	pshufd	$0x4e, %xmm3, %xmm4
	mulpd	%xmm0, %xmm3
	mulpd	%xmm0, %xmm4

	addpd	%xmm6, %xmm8
	movaps	  4 * SIZE(BO), %xmm6
	addpd	%xmm2, %xmm9
	pshufd	$0x4e, %xmm6, %xmm2
	mulpd	%xmm0, %xmm6
	mulpd	%xmm0, %xmm2

	addpd	%xmm3, %xmm10
	movaps	  6 * SIZE(BO), %xmm3
	addpd	%xmm4, %xmm11
	pshufd	$0x4e, %xmm3, %xmm4
	mulpd	%xmm0, %xmm3
	mulpd	%xmm0, %xmm4

	addpd	%xmm6, %xmm12
	movaps	  8 * SIZE(BO), %xmm1
	addpd	%xmm2, %xmm13
	movaps	-10 * SIZE(AO), %xmm5
	pshufd	$0x4e, %xmm1, %xmm2
	mulpd	%xmm5, %xmm1
	mulpd	%xmm5, %xmm2

	addpd	%xmm3, %xmm14
	movaps	 10 * SIZE(BO), %xmm3
	addpd	%xmm4, %xmm15
	pshufd	$0x4e, %xmm3, %xmm4
	mulpd	%xmm5, %xmm3
	mulpd	%xmm5, %xmm4

	addpd	%xmm1, %xmm8
	movaps	 12 * SIZE(BO), %xmm1
	addpd	%xmm2, %xmm9
	pshufd	$0x4e, %xmm1, %xmm2
	mulpd	%xmm5, %xmm1
	mulpd	%xmm5, %xmm2

	addpd	%xmm3, %xmm10
	movaps	 14 * SIZE(BO), %xmm3
	addpd	%xmm4, %xmm11
	pshufd	$0x4e, %xmm3, %xmm4
	movaps	 -8 * SIZE(AO), %xmm0
	mulpd	%xmm5, %xmm3
	mulpd	%xmm5, %xmm4

	addq	$32 * SIZE, BO
	subq	$-8 * SIZE, AO
	decq	%rax
	BRANCH
	jg	L(12)
	ALIGN_3

L(15):
	movups	ALPHA_R, %xmm7

#ifndef TRMMKERNEL
	movq	K, %rax
#else
	movq	KKK, %rax
#endif
	andq	$3, %rax		# if (k & 1)
	BRANCH
	je	L(18)
	ALIGN_3

L(16):
	addpd	%xmm1, %xmm12
	movaps	-16 * SIZE(BO), %xmm1
	addpd	%xmm2, %xmm13
	pshufd	$0x4e, %xmm1, %xmm2
	mulpd	%xmm0, %xmm1
	mulpd	%xmm0, %xmm2

	addpd	%xmm3, %xmm14
	movaps	-14 * SIZE(BO), %xmm3
	addpd	%xmm4, %xmm15
	pshufd	$0x4e, %xmm3, %xmm4
	mulpd	%xmm0, %xmm3
	mulpd	%xmm0, %xmm4

	addpd	%xmm1, %xmm8
	movaps	-12 * SIZE(BO), %xmm1
	addpd	%xmm2, %xmm9
	pshufd	$0x4e, %xmm1, %xmm2
	mulpd	%xmm0, %xmm1
	mulpd	%xmm0, %xmm2

	addpd	%xmm3, %xmm10
	movaps	-10 * SIZE(BO), %xmm3
	addpd	%xmm4, %xmm11
	pshufd	$0x4e, %xmm3, %xmm4
	mulpd	%xmm0, %xmm3
	mulpd	%xmm0, %xmm4

	movaps	-14 * SIZE(AO), %xmm0

	addq	$2 * SIZE, AO
	addq	$8 * SIZE, BO

	subq	$1, %rax
	BRANCH
	jg	L(16)
	ALIGN_4

L(18):
	addpd	%xmm1, %xmm12
	addpd	%xmm2, %xmm13
	addpd	%xmm3, %xmm14
	addpd	%xmm4, %xmm15

	movaps	%xmm8, %xmm0
	shufpd	$2, %xmm9, %xmm8
	shufpd	$2, %xmm0, %xmm9
	movaps	%xmm10, %xmm0
	shufpd	$2, %xmm11, %xmm10
	shufpd	$2, %xmm0, %xmm11

	movaps	%xmm12, %xmm0
	shufpd	$2, %xmm13, %xmm12
	shufpd	$2, %xmm0, %xmm13
	movaps	%xmm14, %xmm0
	shufpd	$2, %xmm15, %xmm14
	shufpd	$2, %xmm0, %xmm15

	leaq	(LDC, LDC, 2), %rax

	movsd	0 * SIZE(CO1), %xmm0
	movhps	1 * SIZE(CO1), %xmm0
	movsd	2 * SIZE(CO1), %xmm1
	movhps	3 * SIZE(CO1), %xmm1

	movsd	0 * SIZE(CO1, LDC), %xmm2
	movhps	1 * SIZE(CO1, LDC), %xmm2
	movsd	2 * SIZE(CO1, LDC), %xmm3
	movhps	3 * SIZE(CO1, LDC), %xmm3

	movddup	 %xmm8,  %xmm4
	mulpd	 %xmm7,  %xmm4
	addpd	 %xmm4,  %xmm0

	unpckhpd %xmm8,  %xmm8
	mulpd	 %xmm7,  %xmm8
	addpd	 %xmm8,  %xmm1

	movddup	 %xmm9,  %xmm5
	mulpd	 %xmm7,  %xmm5
	addpd	 %xmm5,  %xmm2

	unpckhpd %xmm9,  %xmm9
	mulpd	 %xmm7,  %xmm9
	addpd	 %xmm9,  %xmm3

	movlps	%xmm0,  0 * SIZE(CO1)
	movhps	%xmm0,  1 * SIZE(CO1)
	movlps	%xmm1,  2 * SIZE(CO1)
	movhps	%xmm1,  3 * SIZE(CO1)

	movlps	%xmm2,  0 * SIZE(CO1, LDC)
	movhps	%xmm2,  1 * SIZE(CO1, LDC)
	movlps	%xmm3,  2 * SIZE(CO1, LDC)
	movhps	%xmm3,  3 * SIZE(CO1, LDC)

	movsd	0 * SIZE(CO1, LDC, 2), %xmm0
	movhps	1 * SIZE(CO1, LDC, 2), %xmm0
	movsd	2 * SIZE(CO1, LDC, 2), %xmm1
	movhps	3 * SIZE(CO1, LDC, 2), %xmm1

	movsd	0 * SIZE(CO1, %rax), %xmm2
	movhps	1 * SIZE(CO1, %rax), %xmm2
	movsd	2 * SIZE(CO1, %rax), %xmm3
	movhps	3 * SIZE(CO1, %rax), %xmm3

	movddup	 %xmm10,  %xmm4
	mulpd	 %xmm7,  %xmm4
	addpd	 %xmm4,  %xmm0

	unpckhpd %xmm10, %xmm10
	mulpd	 %xmm7,  %xmm10
	addpd	 %xmm10, %xmm1

	movddup	 %xmm11, %xmm5
	mulpd	 %xmm7,  %xmm5
	addpd	 %xmm5,  %xmm2

	unpckhpd %xmm11, %xmm11
	mulpd	 %xmm7,  %xmm11
	addpd	 %xmm11, %xmm3

	movlps	%xmm0,  0 * SIZE(CO1, LDC, 2)
	movhps	%xmm0,  1 * SIZE(CO1, LDC, 2)
	movlps	%xmm1,  2 * SIZE(CO1, LDC, 2)
	movhps	%xmm1,  3 * SIZE(CO1, LDC, 2)

	movlps	%xmm2,  0 * SIZE(CO1, %rax)
	movhps	%xmm2,  1 * SIZE(CO1, %rax)
	movlps	%xmm3,  2 * SIZE(CO1, %rax)
	movhps	%xmm3,  3 * SIZE(CO1, %rax)

	movsd	0 * SIZE(CO2), %xmm0
	movhps	1 * SIZE(CO2), %xmm0
	movsd	2 * SIZE(CO2), %xmm1
	movhps	3 * SIZE(CO2), %xmm1

	movsd	0 * SIZE(CO2, LDC), %xmm2
	movhps	1 * SIZE(CO2, LDC), %xmm2
	movsd	2 * SIZE(CO2, LDC), %xmm3
	movhps	3 * SIZE(CO2, LDC), %xmm3

	movddup	 %xmm12, %xmm4
	mulpd	 %xmm7,  %xmm4
	addpd	 %xmm4,  %xmm0

	unpckhpd %xmm12, %xmm12
	mulpd	 %xmm7,  %xmm12
	addpd	 %xmm12, %xmm1

	movddup	 %xmm13, %xmm5
	mulpd	 %xmm7,  %xmm5
	addpd	 %xmm5,  %xmm2

	unpckhpd %xmm13, %xmm13
	mulpd	 %xmm7,  %xmm13
	addpd	 %xmm13, %xmm3

	movlps	%xmm0,  0 * SIZE(CO2)
	movhps	%xmm0,  1 * SIZE(CO2)
	movlps	%xmm1,  2 * SIZE(CO2)
	movhps	%xmm1,  3 * SIZE(CO2)

	movlps	%xmm2,  0 * SIZE(CO2, LDC)
	movhps	%xmm2,  1 * SIZE(CO2, LDC)
	movlps	%xmm3,  2 * SIZE(CO2, LDC)
	movhps	%xmm3,  3 * SIZE(CO2, LDC)

	movsd	0 * SIZE(CO2, LDC, 2), %xmm0
	movhps	1 * SIZE(CO2, LDC, 2), %xmm0
	movsd	2 * SIZE(CO2, LDC, 2), %xmm1
	movhps	3 * SIZE(CO2, LDC, 2), %xmm1

	movsd	0 * SIZE(CO2, %rax), %xmm2
	movhps	1 * SIZE(CO2, %rax), %xmm2
	movsd	2 * SIZE(CO2, %rax), %xmm3
	movhps	3 * SIZE(CO2, %rax), %xmm3

	movddup	 %xmm14, %xmm4
	mulpd	 %xmm7,  %xmm4
	addpd	 %xmm4,  %xmm0

	unpckhpd %xmm14, %xmm14
	mulpd	 %xmm7,  %xmm14
	addpd	 %xmm14, %xmm1

	movddup	 %xmm15, %xmm5
	mulpd	 %xmm7,  %xmm5
	addpd	 %xmm5,  %xmm2

	unpckhpd %xmm15, %xmm15
	mulpd	 %xmm7,  %xmm15
	addpd	 %xmm15, %xmm3

	movlps	%xmm0,  0 * SIZE(CO2, LDC, 2)
	movhps	%xmm0,  1 * SIZE(CO2, LDC, 2)
	movlps	%xmm1,  2 * SIZE(CO2, LDC, 2)
	movhps	%xmm1,  3 * SIZE(CO2, LDC, 2)

	movlps	%xmm2,  0 * SIZE(CO2, %rax)
	movhps	%xmm2,  1 * SIZE(CO2, %rax)
	movlps	%xmm3,  2 * SIZE(CO2, %rax)
	movhps	%xmm3,  3 * SIZE(CO2, %rax)

	addq	$4 * SIZE, CO1
	addq	$4 * SIZE, CO2

	decq	I
	BRANCH
	jg	L(11)
	ALIGN_4

L(20):
	testq	$1, M
	BRANCH
	jle	L(29)
	ALIGN_4

#if !defined(TRMMKERNEL) || \
	(defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
	(defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))

	movq	B, BO
#else
	movq	B, BO

	movq	KK, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 1), AO
	leaq	(BO, %rax, 8), BO
#endif

	movddup	-16 * SIZE(AO), %xmm0
	xorps	%xmm8,  %xmm8
	movaps	-16 * SIZE(BO), %xmm1
	xorps	%xmm9,  %xmm9
	xorps	%xmm10, %xmm10
	xorps	%xmm11, %xmm11

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
	addq	$8, %rax
#endif
	movq	%rax, KKK
#endif
	sarq	$2, %rax
	NOBRANCH
	jle	L(25)
	ALIGN_3

L(22):
	mulpd	%xmm0, %xmm1
	addpd	%xmm1, %xmm8
	movaps	-14 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	addpd	%xmm1, %xmm9
	movaps	-12 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	addpd	%xmm1, %xmm10
	movaps	-10 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	movddup	-15 * SIZE(AO), %xmm0
	addpd	%xmm1, %xmm11
	movaps	 -8 * SIZE(BO), %xmm1

	mulpd	%xmm0, %xmm1
	addpd	%xmm1, %xmm8
	movaps	 -6 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	addpd	%xmm1, %xmm9
	movaps	 -4 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	addpd	%xmm1, %xmm10
	movaps	 -2 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	movddup	-14 * SIZE(AO), %xmm0
	addpd	%xmm1, %xmm11
	movaps	  0 * SIZE(BO), %xmm1

	mulpd	%xmm0, %xmm1
	addpd	%xmm1, %xmm8
	movaps	  2 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	addpd	%xmm1, %xmm9
	movaps	  4 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	addpd	%xmm1, %xmm10
	movaps	  6 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	movddup	-13 * SIZE(AO), %xmm0
	addpd	%xmm1, %xmm11
	movaps	  8 * SIZE(BO), %xmm1

	mulpd	%xmm0, %xmm1
	addpd	%xmm1, %xmm8
	movaps	 10 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	addpd	%xmm1, %xmm9
	movaps	 12 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	addpd	%xmm1, %xmm10
	movaps	 14 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	movddup	-12 * SIZE(AO), %xmm0
	addpd	%xmm1, %xmm11
	movaps	 16 * SIZE(BO), %xmm1

	subq	$ -4 * SIZE, AO
	subq	$-32 * SIZE, BO

	subq	$1, %rax
	BRANCH
	jg	L(22)
	ALIGN_3

L(25):
	movups	ALPHA_R, %xmm7

#ifndef TRMMKERNEL
	movq	K, %rax
#else
	movq	KKK, %rax
#endif
	andq	$3, %rax		# if (k & 1)
	BRANCH
	je	L(28)
	ALIGN_3

L(26):
	mulpd	%xmm0, %xmm1
	addpd	%xmm1, %xmm8
	movaps	-14 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	addpd	%xmm1, %xmm9
	movaps	-12 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	addpd	%xmm1, %xmm10
	movaps	-10 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	movddup	-15 * SIZE(AO), %xmm0
	addpd	%xmm1, %xmm11
	movaps	 -8 * SIZE(BO), %xmm1

	addq	$1 * SIZE, AO
	addq	$8 * SIZE, BO

	subq	$1, %rax
	BRANCH
	jg	L(26)
	ALIGN_4

L(28):
	leaq	(LDC, LDC, 2), %rax

	movsd	0 * SIZE(CO1), %xmm0
	movhps	1 * SIZE(CO1), %xmm0
	movsd	0 * SIZE(CO1, LDC), %xmm1
	movhps	1 * SIZE(CO1, LDC), %xmm1

	movsd	0 * SIZE(CO1, LDC, 2), %xmm2
	movhps	1 * SIZE(CO1, LDC, 2), %xmm2
	movsd	0 * SIZE(CO1, %rax), %xmm3
	movhps	1 * SIZE(CO1, %rax), %xmm3

	movddup	 %xmm8,  %xmm4
	mulpd	 %xmm7,  %xmm4
	addpd	 %xmm4,  %xmm0

	unpckhpd %xmm8,  %xmm8
	mulpd	 %xmm7,  %xmm8
	addpd	 %xmm8,  %xmm1

	movddup	 %xmm9,  %xmm5
	mulpd	 %xmm7,  %xmm5
	addpd	 %xmm5,  %xmm2

	unpckhpd %xmm9,  %xmm9
	mulpd	 %xmm7,  %xmm9
	addpd	 %xmm9,  %xmm3

	movlps	%xmm0,  0 * SIZE(CO1)
	movhps	%xmm0,  1 * SIZE(CO1)
	movlps	%xmm1,  0 * SIZE(CO1, LDC)
	movhps	%xmm1,  1 * SIZE(CO1, LDC)

	movlps	%xmm2,  0 * SIZE(CO1, LDC, 2)
	movhps	%xmm2,  1 * SIZE(CO1, LDC, 2)
	movlps	%xmm3,  0 * SIZE(CO1, %rax)
	movhps	%xmm3,  1 * SIZE(CO1, %rax)

	movsd	0 * SIZE(CO2), %xmm0
	movhps	1 * SIZE(CO2), %xmm0
	movsd	0 * SIZE(CO2, LDC), %xmm1
	movhps	1 * SIZE(CO2, LDC), %xmm1

	movsd	0 * SIZE(CO2, LDC, 2), %xmm2
	movhps	1 * SIZE(CO2, LDC, 2), %xmm2
	movsd	0 * SIZE(CO2, %rax), %xmm3
	movhps	1 * SIZE(CO2, %rax), %xmm3

	movddup	 %xmm10, %xmm4
	mulpd	 %xmm7,  %xmm4
	addpd	 %xmm4,  %xmm0

	unpckhpd %xmm10, %xmm10
	mulpd	 %xmm7,  %xmm10
	addpd	 %xmm10, %xmm1

	movddup	 %xmm11, %xmm5
	mulpd	 %xmm7,  %xmm5
	addpd	 %xmm5,  %xmm2

	unpckhpd %xmm11, %xmm11
	mulpd	 %xmm7,  %xmm11
	addpd	 %xmm11, %xmm3

	movlps	%xmm0,  0 * SIZE(CO2)
	movhps	%xmm0,  1 * SIZE(CO2)
	movlps	%xmm1,  0 * SIZE(CO2, LDC)
	movhps	%xmm1,  1 * SIZE(CO2, LDC)

	movlps	%xmm2,  0 * SIZE(CO2, LDC, 2)
	movhps	%xmm2,  1 * SIZE(CO2, LDC, 2)
	movlps	%xmm3,  0 * SIZE(CO2, %rax)
	movhps	%xmm3,  1 * SIZE(CO2, %rax)
	ALIGN_4

L(29):
#if defined(TRMMKERNEL) && !defined(LEFT)
	addq	$8, KK
#endif

	movq	BO, B

	leaq	(C, LDC, 8), C

	subq	$1, J
	BRANCH
	jg	L(01)
	ALIGN_4

L(30):
	testq	$4, N
	jle	L(50)
	ALIGN_4

#if defined(TRMMKERNEL) && defined(LEFT)
        movq    OFFSET, %rax
	movq    %rax, KK
#endif

	movq	C, CO1
	leaq	(C, LDC, 2), CO2
	movq	A, AO

	movq	M,  I
	sarq	$1, I
	NOBRANCH
	jle	L(40)
	ALIGN_4

L(31):
#if !defined(TRMMKERNEL) || \
	(defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
	(defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))

	movq	B, BO
#else
	movq	B, BO

	movq	KK, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 2), AO
	leaq	(BO, %rax, 4), BO
#endif

	xorps	%xmm1, %xmm1
	movaps	-16 * SIZE(AO), %xmm0
	xorps	%xmm2, %xmm2
	xorps	%xmm3, %xmm3
	xorps	%xmm4, %xmm4

	xorps	%xmm8,  %xmm8
	prefetcht0     3 * SIZE(CO1)
	xorps	%xmm9,  %xmm9
	prefetcht0     7 * SIZE(CO1, LDC,  1)
	xorps	%xmm10, %xmm10
	prefetcht0     3 * SIZE(CO2)
	xorps	%xmm11, %xmm11
	prefetcht0     7 * SIZE(CO2, LDC,  1)

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
	sarq	$2, %rax
	NOBRANCH
	jle	L(35)
	ALIGN_3

L(32):
	PREFETCH	(PREFETCHSIZE +  0) * SIZE(AO)

	addpd	%xmm1, %xmm8
	movaps	-16 * SIZE(BO), %xmm1
	addpd	%xmm2, %xmm9
	pshufd	$0x4e, %xmm1, %xmm2
	mulpd	%xmm0, %xmm1
	mulpd	%xmm0, %xmm2

	addpd	%xmm3, %xmm10
	movaps	-14 * SIZE(BO), %xmm3
	addpd	%xmm4, %xmm11
	pshufd	$0x4e, %xmm3, %xmm4
	mulpd	%xmm0, %xmm3
	mulpd	%xmm0, %xmm4

	movaps	-14 * SIZE(AO), %xmm0

	addpd	%xmm1, %xmm8
	movaps	-12 * SIZE(BO), %xmm1
	addpd	%xmm2, %xmm9
	pshufd	$0x4e, %xmm1, %xmm2
	mulpd	%xmm0, %xmm1
	mulpd	%xmm0, %xmm2

	addpd	%xmm3, %xmm10
	movaps	-10 * SIZE(BO), %xmm3
	addpd	%xmm4, %xmm11
	pshufd	$0x4e, %xmm3, %xmm4
	mulpd	%xmm0, %xmm3
	mulpd	%xmm0, %xmm4

	movaps	-12 * SIZE(AO), %xmm0

	addpd	%xmm1, %xmm8
	movaps	 -8 * SIZE(BO), %xmm1
	addpd	%xmm2, %xmm9
	pshufd	$0x4e, %xmm1, %xmm2
	mulpd	%xmm0, %xmm1
	mulpd	%xmm0, %xmm2

	addpd	%xmm3, %xmm10
	movaps	 -6 * SIZE(BO), %xmm3
	addpd	%xmm4, %xmm11
	pshufd	$0x4e, %xmm3, %xmm4
	mulpd	%xmm0, %xmm3
	mulpd	%xmm0, %xmm4

	movaps	-10 * SIZE(AO), %xmm0

	addpd	%xmm1, %xmm8
	movaps	 -4 * SIZE(BO), %xmm1
	addpd	%xmm2, %xmm9
	pshufd	$0x4e, %xmm1, %xmm2
	mulpd	%xmm0, %xmm1
	mulpd	%xmm0, %xmm2

	addpd	%xmm3, %xmm10
	movaps	 -2 * SIZE(BO), %xmm3
	addpd	%xmm4, %xmm11
	pshufd	$0x4e, %xmm3, %xmm4
	mulpd	%xmm0, %xmm3
	mulpd	%xmm0, %xmm4

	movaps	 -8 * SIZE(AO), %xmm0

	subq	$-8  * SIZE, AO
	subq	$-16 * SIZE, BO

	subq	$1, %rax
	BRANCH
	jg	L(32)
	ALIGN_3

L(35):
	movups	ALPHA_R, %xmm7

#ifndef TRMMKERNEL
	movq	K, %rax
#else
	movq	KKK, %rax
#endif
	andq	$3, %rax		# if (k & 1)
	BRANCH
	je	L(38)
	ALIGN_3

L(36):
	addpd	%xmm1, %xmm8
	movaps	-16 * SIZE(BO), %xmm1
	addpd	%xmm2, %xmm9
	pshufd	$0x4e, %xmm1, %xmm2
	mulpd	%xmm0, %xmm1
	mulpd	%xmm0, %xmm2

	addpd	%xmm3, %xmm10
	movaps	-14 * SIZE(BO), %xmm3
	addpd	%xmm4, %xmm11
	pshufd	$0x4e, %xmm3, %xmm4
	mulpd	%xmm0, %xmm3
	mulpd	%xmm0, %xmm4

	movaps	-14 * SIZE(AO), %xmm0

	addq	$2 * SIZE, AO
	addq	$4 * SIZE, BO

	subq	$1, %rax
	BRANCH
	jg	L(36)
	ALIGN_4

L(38):
	addpd	%xmm1, %xmm8
	addpd	%xmm2, %xmm9
	addpd	%xmm3, %xmm10
	addpd	%xmm4, %xmm11

	movaps	%xmm8, %xmm0
	shufpd	$2, %xmm9, %xmm8
	shufpd	$2, %xmm0, %xmm9

	movaps	%xmm10, %xmm0
	shufpd	$2, %xmm11, %xmm10
	shufpd	$2, %xmm0,  %xmm11

	movsd	0 * SIZE(CO1), %xmm0
	movhps	1 * SIZE(CO1), %xmm0
	movsd	2 * SIZE(CO1), %xmm1
	movhps	3 * SIZE(CO1), %xmm1

	movsd	0 * SIZE(CO1, LDC), %xmm2
	movhps	1 * SIZE(CO1, LDC), %xmm2
	movsd	2 * SIZE(CO1, LDC), %xmm3
	movhps	3 * SIZE(CO1, LDC), %xmm3

	movddup	 %xmm8,  %xmm4
	mulpd	 %xmm7,  %xmm4
	addpd	 %xmm4,  %xmm0

	unpckhpd %xmm8,  %xmm8
	mulpd	 %xmm7,  %xmm8
	addpd	 %xmm8,  %xmm1

	movddup	 %xmm9,  %xmm5
	mulpd	 %xmm7,  %xmm5
	addpd	 %xmm5,  %xmm2

	unpckhpd %xmm9,  %xmm9
	mulpd	 %xmm7,  %xmm9
	addpd	 %xmm9,  %xmm3

	movlps	%xmm0,  0 * SIZE(CO1)
	movhps	%xmm0,  1 * SIZE(CO1)
	movlps	%xmm1,  2 * SIZE(CO1)
	movhps	%xmm1,  3 * SIZE(CO1)

	movlps	%xmm2,  0 * SIZE(CO1, LDC)
	movhps	%xmm2,  1 * SIZE(CO1, LDC)
	movlps	%xmm3,  2 * SIZE(CO1, LDC)
	movhps	%xmm3,  3 * SIZE(CO1, LDC)

	movsd	0 * SIZE(CO2), %xmm0
	movhps	1 * SIZE(CO2), %xmm0
	movsd	2 * SIZE(CO2), %xmm1
	movhps	3 * SIZE(CO2), %xmm1

	movsd	0 * SIZE(CO2, LDC), %xmm2
	movhps	1 * SIZE(CO2, LDC), %xmm2
	movsd	2 * SIZE(CO2, LDC), %xmm3
	movhps	3 * SIZE(CO2, LDC), %xmm3

	movddup	 %xmm10, %xmm4
	mulpd	 %xmm7,  %xmm4
	addpd	 %xmm4,  %xmm0

	unpckhpd %xmm10, %xmm10
	mulpd	 %xmm7,  %xmm10
	addpd	 %xmm10, %xmm1

	movddup	 %xmm11, %xmm5
	mulpd	 %xmm7,  %xmm5
	addpd	 %xmm5,  %xmm2

	unpckhpd %xmm11, %xmm11
	mulpd	 %xmm7,  %xmm11
	addpd	 %xmm11, %xmm3

	movlps	%xmm0,  0 * SIZE(CO2)
	movhps	%xmm0,  1 * SIZE(CO2)
	movlps	%xmm1,  2 * SIZE(CO2)
	movhps	%xmm1,  3 * SIZE(CO2)

	movlps	%xmm2,  0 * SIZE(CO2, LDC)
	movhps	%xmm2,  1 * SIZE(CO2, LDC)
	movlps	%xmm3,  2 * SIZE(CO2, LDC)
	movhps	%xmm3,  3 * SIZE(CO2, LDC)

	addq	$4 * SIZE, CO1
	addq	$4 * SIZE, CO2
	decq	I
	BRANCH
	jg	L(31)
	ALIGN_4

L(40):
	testq	$1, M
	BRANCH
	jle	L(49)
	ALIGN_4

#if !defined(TRMMKERNEL) || \
	(defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
	(defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))

	movq	B, BO
#else
	movq	B, BO

	movq	KK, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 1), AO
	leaq	(BO, %rax, 4), BO
#endif

	movddup	-16 * SIZE(AO), %xmm0
	xorps	%xmm8,  %xmm8
	movaps	-16 * SIZE(BO), %xmm1
	xorps	%xmm9,  %xmm9
	xorps	%xmm10, %xmm10
	xorps	%xmm11, %xmm11

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
	sarq	$2, %rax
	NOBRANCH
	jle	L(45)
	ALIGN_3

L(42):
	mulpd	%xmm0, %xmm1
	addpd	%xmm1, %xmm8
	movaps	-14 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	movddup	-15 * SIZE(AO), %xmm0
	addpd	%xmm1, %xmm9
	movaps	-12 * SIZE(BO), %xmm1

	mulpd	%xmm0, %xmm1
	addpd	%xmm1, %xmm10
	movaps	-10 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	movddup	-14 * SIZE(AO), %xmm0
	addpd	%xmm1, %xmm11
	movaps	 -8 * SIZE(BO), %xmm1

	mulpd	%xmm0, %xmm1
	addpd	%xmm1, %xmm8
	movaps	 -6 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	movddup	-13 * SIZE(AO), %xmm0
	addpd	%xmm1, %xmm9
	movaps	 -4 * SIZE(BO), %xmm1

	mulpd	%xmm0, %xmm1
	addpd	%xmm1, %xmm10
	movaps	 -2 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	movddup	-12 * SIZE(AO), %xmm0
	addpd	%xmm1, %xmm11
	movaps	  0 * SIZE(BO), %xmm1

	subq	$ -4 * SIZE, AO
	subq	$-16 * SIZE, BO

	subq	$1, %rax
	BRANCH
	jg	L(42)
	ALIGN_3

L(45):
	movups	ALPHA_R, %xmm7

#ifndef TRMMKERNEL
	movq	K, %rax
#else
	movq	KKK, %rax
#endif
	andq	$3, %rax		# if (k & 1)
	BRANCH
	je	L(48)
	ALIGN_3

L(46):
	mulpd	%xmm0, %xmm1
	addpd	%xmm1, %xmm8
	movaps	-14 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	movddup	-15 * SIZE(AO), %xmm0
	addpd	%xmm1, %xmm9
	movaps	-12 * SIZE(BO), %xmm1

	addq	$1 * SIZE, AO
	addq	$4 * SIZE, BO

	subq	$1, %rax
	BRANCH
	jg	L(46)
	ALIGN_4

L(48):
	addpd	%xmm10, %xmm8
	addpd	%xmm11, %xmm9

	movsd	0 * SIZE(CO1), %xmm0
	movhps	1 * SIZE(CO1), %xmm0
	movsd	0 * SIZE(CO1, LDC), %xmm1
	movhps	1 * SIZE(CO1, LDC), %xmm1

	movsd	0 * SIZE(CO2), %xmm2
	movhps	1 * SIZE(CO2), %xmm2
	movsd	0 * SIZE(CO2, LDC), %xmm3
	movhps	1 * SIZE(CO2, LDC), %xmm3

	movddup	 %xmm8,  %xmm4
	mulpd	 %xmm7,  %xmm4
	addpd	 %xmm4,  %xmm0

	unpckhpd %xmm8,  %xmm8
	mulpd	 %xmm7,  %xmm8
	addpd	 %xmm8,  %xmm1

	movddup	 %xmm9,  %xmm5
	mulpd	 %xmm7,  %xmm5
	addpd	 %xmm5,  %xmm2

	unpckhpd %xmm9,  %xmm9
	mulpd	 %xmm7,  %xmm9
	addpd	 %xmm9,  %xmm3

	movlps	%xmm0,  0 * SIZE(CO1)
	movhps	%xmm0,  1 * SIZE(CO1)
	movlps	%xmm1,  0 * SIZE(CO1, LDC)
	movhps	%xmm1,  1 * SIZE(CO1, LDC)

	movlps	%xmm2,  0 * SIZE(CO2)
	movhps	%xmm2,  1 * SIZE(CO2)
	movlps	%xmm3,  0 * SIZE(CO2, LDC)
	movhps	%xmm3,  1 * SIZE(CO2, LDC)
	ALIGN_4

L(49):
#if defined(TRMMKERNEL) && !defined(LEFT)
	addq	$4, KK
#endif

	movq	BO, B

	leaq	(C, LDC, 4), C
	ALIGN_4

L(50):
	testq	$2, N
	jle	L(70)
	ALIGN_4

#if defined(TRMMKERNEL) && defined(LEFT)
        movq    OFFSET, %rax
	movq    %rax, KK
#endif

	movq	C, CO1
	leaq	(C, LDC), CO2
	movq	A, AO

	movq	M,  I
	sarq	$1, I
	NOBRANCH
	jle	L(60)
	ALIGN_4

L(51):
#if !defined(TRMMKERNEL) || \
	(defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
	(defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))

	movq	B, BO
#else
	movq	B, BO

	movq	KK, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 2), AO
	leaq	(BO, %rax, 2), BO
#endif

	xorps	%xmm1, %xmm1
	movaps	-16 * SIZE(AO), %xmm0
	xorps	%xmm2, %xmm2

	xorps	%xmm8,  %xmm8
	prefetcht0     3 * SIZE(CO1)
	xorps	%xmm9,  %xmm9
	prefetcht0     7 * SIZE(CO2)
	xorps	%xmm10, %xmm10
	xorps	%xmm11, %xmm11

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
	sarq	$2, %rax
	NOBRANCH
	jle	L(55)
	ALIGN_3

L(52):
	PREFETCH	(PREFETCHSIZE +  0) * SIZE(AO)

	addpd	%xmm1, %xmm8
	movaps	-16 * SIZE(BO), %xmm1
	addpd	%xmm2, %xmm9
	pshufd	$0x4e, %xmm1, %xmm2
	mulpd	%xmm0, %xmm1
	mulpd	%xmm0, %xmm2
	movaps	-14 * SIZE(AO), %xmm0

	addpd	%xmm1, %xmm10
	movaps	-14 * SIZE(BO), %xmm1
	addpd	%xmm2, %xmm11
	pshufd	$0x4e, %xmm1, %xmm2
	mulpd	%xmm0, %xmm1
	mulpd	%xmm0, %xmm2
	movaps	-12 * SIZE(AO), %xmm0

	addpd	%xmm1, %xmm8
	movaps	-12 * SIZE(BO), %xmm1
	addpd	%xmm2, %xmm9
	pshufd	$0x4e, %xmm1, %xmm2
	mulpd	%xmm0, %xmm1
	mulpd	%xmm0, %xmm2
	movaps	-10 * SIZE(AO), %xmm0

	addpd	%xmm1, %xmm10
	movaps	-10 * SIZE(BO), %xmm1
	addpd	%xmm2, %xmm11
	pshufd	$0x4e, %xmm1, %xmm2
	mulpd	%xmm0, %xmm1
	mulpd	%xmm0, %xmm2
	movaps	 -8 * SIZE(AO), %xmm0

	subq	$-8 * SIZE, AO
	subq	$-8 * SIZE, BO

	subq	$1, %rax
	BRANCH
	jg	L(52)

	addpd	%xmm10, %xmm8
	addpd	%xmm11, %xmm9
	ALIGN_3

L(55):
	movups	ALPHA_R, %xmm7

#ifndef TRMMKERNEL
	movq	K, %rax
#else
	movq	KKK, %rax
#endif
	andq	$3, %rax		# if (k & 1)
	BRANCH
	je	L(58)
	ALIGN_3

L(56):
	addpd	%xmm1, %xmm8
	movaps	-16 * SIZE(BO), %xmm1
	addpd	%xmm2, %xmm9
	pshufd	$0x4e, %xmm1, %xmm2
	mulpd	%xmm0, %xmm1
	mulpd	%xmm0, %xmm2
	movaps	-14 * SIZE(AO), %xmm0

	addq	$2 * SIZE, AO
	addq	$2 * SIZE, BO

	subq	$1, %rax
	BRANCH
	jg	L(56)
	ALIGN_4

L(58):
	addpd	%xmm1, %xmm8
	addpd	%xmm2, %xmm9

	movaps	%xmm8, %xmm0
	shufpd	$2, %xmm9, %xmm8
	shufpd	$2, %xmm0, %xmm9

	movsd	0 * SIZE(CO1), %xmm0
	movhps	1 * SIZE(CO1), %xmm0
	movsd	2 * SIZE(CO1), %xmm1
	movhps	3 * SIZE(CO1), %xmm1

	movsd	0 * SIZE(CO2), %xmm2
	movhps	1 * SIZE(CO2), %xmm2
	movsd	2 * SIZE(CO2), %xmm3
	movhps	3 * SIZE(CO2), %xmm3

	movddup	 %xmm8,  %xmm4
	mulpd	 %xmm7,  %xmm4
	addpd	 %xmm4,  %xmm0

	unpckhpd %xmm8,  %xmm8
	mulpd	 %xmm7,  %xmm8
	addpd	 %xmm8,  %xmm1

	movddup	 %xmm9,  %xmm5
	mulpd	 %xmm7,  %xmm5
	addpd	 %xmm5,  %xmm2

	unpckhpd %xmm9,  %xmm9
	mulpd	 %xmm7,  %xmm9
	addpd	 %xmm9,  %xmm3

	movlps	%xmm0,  0 * SIZE(CO1)
	movhps	%xmm0,  1 * SIZE(CO1)
	movlps	%xmm1,  2 * SIZE(CO1)
	movhps	%xmm1,  3 * SIZE(CO1)

	movlps	%xmm2,  0 * SIZE(CO2)
	movhps	%xmm2,  1 * SIZE(CO2)
	movlps	%xmm3,  2 * SIZE(CO2)
	movhps	%xmm3,  3 * SIZE(CO2)

	addq	$4 * SIZE, CO1
	addq	$4 * SIZE, CO2
	decq	I
	BRANCH
	jg	L(51)
	ALIGN_4

L(60):
	testq	$1, M
	BRANCH
	jle	L(69)
	ALIGN_4

#if !defined(TRMMKERNEL) || \
	(defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
	(defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))

	movq	B, BO
#else
	movq	B, BO

	movq	KK, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 1), AO
	leaq	(BO, %rax, 2), BO
#endif

	movddup	-16 * SIZE(AO), %xmm0
	xorps	%xmm8,  %xmm8
	movaps	-16 * SIZE(BO), %xmm1
	xorps	%xmm9,  %xmm9

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
	sarq	$2, %rax
	NOBRANCH
	jle	L(65)
	ALIGN_3

L(62):
	mulpd	%xmm0, %xmm1
	movddup	-15 * SIZE(AO), %xmm0
	addpd	%xmm1, %xmm8
	movaps	-14 * SIZE(BO), %xmm1

	mulpd	%xmm0, %xmm1
	movddup	-14 * SIZE(AO), %xmm0
	addpd	%xmm1, %xmm9
	movaps	-12 * SIZE(BO), %xmm1

	mulpd	%xmm0, %xmm1
	movddup	-13 * SIZE(AO), %xmm0
	addpd	%xmm1, %xmm8
	movaps	-10 * SIZE(BO), %xmm1

	mulpd	%xmm0, %xmm1
	movddup	-12 * SIZE(AO), %xmm0
	addpd	%xmm1, %xmm9
	movaps	 -8 * SIZE(BO), %xmm1

	subq	$-4 * SIZE, AO
	subq	$-8 * SIZE, BO

	subq	$1, %rax
	BRANCH
	jg	L(62)
	ALIGN_3

L(65):
	movups	ALPHA_R, %xmm7

#ifndef TRMMKERNEL
	movq	K, %rax
#else
	movq	KKK, %rax
#endif
	andq	$3, %rax		# if (k & 1)
	BRANCH
	je	L(68)
	ALIGN_3

L(66):
	mulpd	%xmm0, %xmm1
	movddup	-15 * SIZE(AO), %xmm0
	addpd	%xmm1, %xmm8
	movaps	-14 * SIZE(BO), %xmm1

	addq	$1 * SIZE, AO
	addq	$2 * SIZE, BO

	subq	$1, %rax
	BRANCH
	jg	L(66)
	ALIGN_4

L(68):
	addpd	%xmm9, %xmm8

	movsd	0 * SIZE(CO1), %xmm0
	movhps	1 * SIZE(CO1), %xmm0
	movsd	0 * SIZE(CO2), %xmm1
	movhps	1 * SIZE(CO2), %xmm1

	movddup	 %xmm8,  %xmm4
	mulpd	 %xmm7,  %xmm4
	addpd	 %xmm4,  %xmm0

	unpckhpd %xmm8,  %xmm8
	mulpd	 %xmm7,  %xmm8
	addpd	 %xmm8,  %xmm1

	movlps	%xmm0,  0 * SIZE(CO1)
	movhps	%xmm0,  1 * SIZE(CO1)
	movlps	%xmm1,  0 * SIZE(CO2)
	movhps	%xmm1,  1 * SIZE(CO2)
	ALIGN_4

L(69):
#if defined(TRMMKERNEL) && !defined(LEFT)
	addq	$2, KK
#endif

	movq	BO, B

	leaq	(C, LDC, 2), C
	ALIGN_4

L(70):
	testq	$1, N
	jle	L(999)
	ALIGN_4

#if defined(TRMMKERNEL) && defined(LEFT)
        movq    OFFSET, %rax
	movq    %rax, KK
#endif

	movq	C, CO1
	movq	A, AO

	movq	M,  I
	sarq	$1, I
	NOBRANCH
	jle	L(80)
	ALIGN_4

L(71):
#if !defined(TRMMKERNEL) || \
	(defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
	(defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))

	movq	B, BO
#else
	movq	B, BO

	movq	KK, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 2), AO
	leaq	(BO, %rax, 1), BO
#endif

	xorps	%xmm1, %xmm1
	movaps	-16 * SIZE(AO), %xmm0
	xorps	%xmm2, %xmm2

	xorps	%xmm8,  %xmm8
	prefetcht0     3 * SIZE(CO1)
	xorps	%xmm9,  %xmm9
	xorps	%xmm10, %xmm10
	xorps	%xmm11, %xmm11

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
	sarq	$2, %rax
	NOBRANCH
	jle	L(75)
	ALIGN_3

L(72):
	PREFETCH	(PREFETCHSIZE +  0) * SIZE(AO)

	addpd	%xmm1, %xmm8
	movddup	-16 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	movaps	-14 * SIZE(AO), %xmm0

	addpd	%xmm1, %xmm9
	movddup	-15 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	movaps	-12 * SIZE(AO), %xmm0

	addpd	%xmm1, %xmm8
	movddup	-14 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	movaps	-10 * SIZE(AO), %xmm0

	addpd	%xmm1, %xmm9
	movddup	-13 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	movaps	 -8 * SIZE(AO), %xmm0

	subq	$-8 * SIZE, AO
	subq	$-4 * SIZE, BO

	subq	$1, %rax
	BRANCH
	jg	L(72)

	addpd	%xmm9, %xmm8
	ALIGN_3

L(75):
	movups	ALPHA_R, %xmm7

#ifndef TRMMKERNEL
	movq	K, %rax
#else
	movq	KKK, %rax
#endif
	andq	$3, %rax		# if (k & 1)
	BRANCH
	je	L(78)
	ALIGN_3

L(76):
	addpd	%xmm1, %xmm8
	movddup	-16 * SIZE(BO), %xmm1
	mulpd	%xmm0, %xmm1
	movaps	-14 * SIZE(AO), %xmm0

	addq	$2 * SIZE, AO
	addq	$1 * SIZE, BO

	subq	$1, %rax
	BRANCH
	jg	L(76)
	ALIGN_4

L(78):
	addpd	%xmm1, %xmm8

	movsd	0 * SIZE(CO1), %xmm0
	movhps	1 * SIZE(CO1), %xmm0
	movsd	2 * SIZE(CO1), %xmm1
	movhps	3 * SIZE(CO1), %xmm1

	movddup	 %xmm8,  %xmm4
	mulpd	 %xmm7,  %xmm4
	addpd	 %xmm4,  %xmm0

	unpckhpd %xmm8,  %xmm8
	mulpd	 %xmm7,  %xmm8
	addpd	 %xmm8,  %xmm1

	movlps	%xmm0,  0 * SIZE(CO1)
	movhps	%xmm0,  1 * SIZE(CO1)
	movlps	%xmm1,  2 * SIZE(CO1)
	movhps	%xmm1,  3 * SIZE(CO1)

	addq	$4 * SIZE, CO1
	decq	I
	BRANCH
	jg	L(71)
	ALIGN_4

L(80):
	testq	$1, M
	BRANCH
	jle	L(999)
	ALIGN_4

#if !defined(TRMMKERNEL) || \
	(defined(TRMMKERNEL) &&  defined(LEFT) &&  defined(TRANSA)) || \
	(defined(TRMMKERNEL) && !defined(LEFT) && !defined(TRANSA))

	movq	B, BO
#else
	movq	B, BO

	movq	KK, %rax
	leaq	(, %rax, SIZE), %rax
	leaq	(AO, %rax, 1), AO
	leaq	(BO, %rax, 1), BO
#endif

#ifndef TRMMKERNEL
	movaps	-16 * SIZE(AO), %xmm0
	xorps	%xmm8,  %xmm8
	movaps	-16 * SIZE(BO), %xmm1
	xorps	%xmm9,  %xmm9
#else
	movsd	-16 * SIZE(AO), %xmm0
	movhpd	-15 * SIZE(AO), %xmm0
	xorps	%xmm8,  %xmm8
	movsd	-16 * SIZE(BO), %xmm1
	movhpd	-15 * SIZE(BO), %xmm1
	xorps	%xmm9,  %xmm9
#endif

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
	sarq	$2, %rax
	NOBRANCH
	jle	L(85)
	ALIGN_3

L(82):
	mulpd	%xmm0, %xmm1
#ifndef TRMMKERNEL
	movapd	-14 * SIZE(AO), %xmm0
#else
	movsd	-14 * SIZE(AO), %xmm0
	movhpd	-13 * SIZE(AO), %xmm0
#endif
	addpd	%xmm1, %xmm8
#ifndef TRMMKERNEL
	movapd	-14 * SIZE(BO), %xmm1
#else
	movsd	-14 * SIZE(BO), %xmm1
	movhpd	-13 * SIZE(BO), %xmm1
#endif

	mulpd	%xmm0, %xmm1
#ifndef TRMMKERNEL
	movapd	-12 * SIZE(AO), %xmm0
#else
	movsd	-12 * SIZE(AO), %xmm0
	movhpd	-11 * SIZE(AO), %xmm0
#endif
	addpd	%xmm1, %xmm9
#ifndef TRMMKERNEL
	movapd	-12 * SIZE(BO), %xmm1
#else
	movsd	-12 * SIZE(BO), %xmm1
	movhpd	-11 * SIZE(BO), %xmm1
#endif

	subq	$-4 * SIZE, AO
	subq	$-4 * SIZE, BO

	subq	$1, %rax
	BRANCH
	jg	L(82)

	addpd	%xmm9, %xmm8
	ALIGN_3

L(85):
	movups	ALPHA_R, %xmm7

#ifndef TRMMKERNEL
	movq	K, %rax
#else
	movq	KKK, %rax
#endif
	andq	$3, %rax		# if (k & 1)
	BRANCH
	je	L(88)
	ALIGN_3

L(86):
	mulsd	%xmm0, %xmm1
	movsd	-15 * SIZE(AO), %xmm0
	addsd	%xmm1, %xmm8
	movsd	-15 * SIZE(BO), %xmm1

	addq	$1 * SIZE, AO
	addq	$1 * SIZE, BO

	subq	$1, %rax
	BRANCH
	jg	L(86)
	ALIGN_4

L(88):
	haddpd	%xmm8, %xmm8

	movsd	0 * SIZE(CO1), %xmm0
	movhps	1 * SIZE(CO1), %xmm0

	movddup	 %xmm8,  %xmm4
	mulpd	 %xmm7,  %xmm4
	addpd	 %xmm4,  %xmm0

	movlps	%xmm0,  0 * SIZE(CO1)
	movhps	%xmm0,  1 * SIZE(CO1)
	ALIGN_4

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
