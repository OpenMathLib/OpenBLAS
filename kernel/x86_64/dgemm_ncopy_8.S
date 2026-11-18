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

#ifdef NEHALEM
#define PREFETCHSIZE	12
#define PREFETCH	prefetcht0
#define PREFETCHW	prefetcht0
#endif

#ifdef SANDYBRIDGE
#define PREFETCHSIZE	12
#define PREFETCH	prefetcht0
#define PREFETCHW	prefetcht0
#endif

#ifndef MOVAPS
#define MOVAPS	movaps
#endif

#ifndef WINDOWS_ABI

#define M	ARG1	/* rdi */
#define N	ARG2	/* rsi */
#define A	ARG3	/* rdx */
#define LDA	ARG4	/* rcx */
#define B	ARG5	/* r8  */

#define AO1	%r9
#define AO2	%r10
#define LDA3	%r11
#define J	%r12
#define MM	%r13

#else

#define STACKSIZE 128

#define M	ARG1	/* rcx */
#define N	ARG2	/* rdx */
#define A	ARG3	/* r8  */
#define LDA	ARG4	/* r9  */
#define OLD_B		40 + 32 + STACKSIZE(%rsp)

#define B	%r15

#define AO1	%r10
#define AO2	%r11
#define LDA3	%r12
#define J	%r13
#define MM	%r14

#endif

#define I	%rax

	PROLOGUE
	PROFCODE

#ifdef WINDOWS_ABI
	pushq	%r15
	pushq	%r14
#endif
	pushq	%r13
	pushq	%r12

#ifdef WINDOWS_ABI
	subq	$STACKSIZE, %rsp

	movups	%xmm6,    0(%rsp)
	movups	%xmm7,   16(%rsp)
	movups	%xmm8,   32(%rsp)
	movups	%xmm9,   48(%rsp)
	movups	%xmm10,  64(%rsp)
	movups	%xmm11,  80(%rsp)
	movups	%xmm12,  96(%rsp)

	movq	OLD_B,     B
#endif

	leaq	(,LDA, SIZE), LDA
	leaq	(LDA, LDA, 2), LDA3
	subq	$-16 * SIZE, B

	movq	M, MM
	leaq	-1(M), %rax
	testq	$SIZE, A
	cmovne	%rax, MM

	testq	$SIZE, LDA
	jne	L(50)

	movq	N,  J
	sarq	$3, J
	jle	L(20)
	ALIGN_4

L(11):
	movq	A, AO1
	leaq	(A, LDA, 4), AO2
	leaq	(A, LDA, 8), A

	testq	$SIZE, A
	je	L(12)

	movsd	0 * SIZE(AO1),         %xmm0
	movsd	0 * SIZE(AO1, LDA),    %xmm1
	movsd	0 * SIZE(AO1, LDA, 2), %xmm2
	movsd	0 * SIZE(AO1, LDA3),   %xmm3

	movsd	0 * SIZE(AO2),         %xmm4
	movsd	0 * SIZE(AO2, LDA),    %xmm5
	movsd	0 * SIZE(AO2, LDA, 2), %xmm6
	movsd	0 * SIZE(AO2, LDA3),   %xmm7

	unpcklpd %xmm1, %xmm0
	unpcklpd %xmm3, %xmm2
	unpcklpd %xmm5, %xmm4
	unpcklpd %xmm7, %xmm6

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm2,  -14 * SIZE(B)
	movaps	%xmm4,  -12 * SIZE(B)
	movaps	%xmm6,  -10 * SIZE(B)

	addq	$1 * SIZE, AO1
	addq	$1 * SIZE, AO2
	subq	$-8 * SIZE, B
	ALIGN_3

L(12):
	movq	MM, I
	sarq	$3, I
	jle	L(14)
	ALIGN_4

L(13):
#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * SIZE(AO1)
#endif

	MOVAPS	0 * SIZE(AO1),         %xmm0
	MOVAPS	0 * SIZE(AO1, LDA),    %xmm1
	MOVAPS	0 * SIZE(AO1, LDA, 2), %xmm2
	MOVAPS	0 * SIZE(AO1, LDA3),   %xmm3

	movaps	 %xmm0, %xmm8
	unpcklpd %xmm1, %xmm0
	movaps	 %xmm2, %xmm9
	unpcklpd %xmm3, %xmm2

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * SIZE(AO1, LDA)
#endif

	MOVAPS	0 * SIZE(AO2),         %xmm4
	MOVAPS	0 * SIZE(AO2, LDA),    %xmm5
	MOVAPS	0 * SIZE(AO2, LDA, 2), %xmm6
	MOVAPS	0 * SIZE(AO2, LDA3),   %xmm7

	movaps	 %xmm4, %xmm10
	unpcklpd %xmm5, %xmm4
	movaps	 %xmm6, %xmm11
	unpcklpd %xmm7, %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 +  0) * SIZE(B)
#endif

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm2,  -14 * SIZE(B)
	movaps	%xmm4,  -12 * SIZE(B)
	movaps	%xmm6,  -10 * SIZE(B)

	unpckhpd %xmm1, %xmm8
	unpckhpd %xmm3, %xmm9
	unpckhpd %xmm5, %xmm10
	unpckhpd %xmm7, %xmm11

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 +  8) * SIZE(B)
#endif

	movaps	%xmm8,   -8 * SIZE(B)
	movaps	%xmm9,   -6 * SIZE(B)
	movaps	%xmm10,  -4 * SIZE(B)
	movaps	%xmm11,  -2 * SIZE(B)

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * SIZE(AO1, LDA, 2)
#endif

	MOVAPS	2 * SIZE(AO1),         %xmm0
	MOVAPS	2 * SIZE(AO1, LDA),    %xmm1
	MOVAPS	2 * SIZE(AO1, LDA, 2), %xmm2
	MOVAPS	2 * SIZE(AO1, LDA3),   %xmm3

	movaps	 %xmm0, %xmm8
	unpcklpd %xmm1, %xmm0
	movaps	 %xmm2, %xmm9
	unpcklpd %xmm3, %xmm2

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * SIZE(AO1, LDA3)
#endif

	MOVAPS	2 * SIZE(AO2),         %xmm4
	MOVAPS	2 * SIZE(AO2, LDA),    %xmm5
	MOVAPS	2 * SIZE(AO2, LDA, 2), %xmm6
	MOVAPS	2 * SIZE(AO2, LDA3),   %xmm7

	movaps	 %xmm4, %xmm10
	unpcklpd %xmm5, %xmm4
	movaps	 %xmm6, %xmm11
	unpcklpd %xmm7, %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 + 16) * SIZE(B)
#endif

	movaps	%xmm0,    0 * SIZE(B)
	movaps	%xmm2,    2 * SIZE(B)
	movaps	%xmm4,    4 * SIZE(B)
	movaps	%xmm6,    6 * SIZE(B)

	unpckhpd %xmm1, %xmm8
	unpckhpd %xmm3, %xmm9
	unpckhpd %xmm5, %xmm10
	unpckhpd %xmm7, %xmm11

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 + 24) * SIZE(B)
#endif

	movaps	%xmm8,    8 * SIZE(B)
	movaps	%xmm9,   10 * SIZE(B)
	movaps	%xmm10,  12 * SIZE(B)
	movaps	%xmm11,  14 * SIZE(B)

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * SIZE(AO2)
#endif

	MOVAPS	4 * SIZE(AO1),         %xmm0
	MOVAPS	4 * SIZE(AO1, LDA),    %xmm1
	MOVAPS	4 * SIZE(AO1, LDA, 2), %xmm2
	MOVAPS	4 * SIZE(AO1, LDA3),   %xmm3

	movaps	 %xmm0, %xmm8
	unpcklpd %xmm1, %xmm0
	movaps	 %xmm2, %xmm9
	unpcklpd %xmm3, %xmm2

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * SIZE(AO2, LDA)
#endif

	MOVAPS	4 * SIZE(AO2),         %xmm4
	MOVAPS	4 * SIZE(AO2, LDA),    %xmm5
	MOVAPS	4 * SIZE(AO2, LDA, 2), %xmm6
	MOVAPS	4 * SIZE(AO2, LDA3),   %xmm7

	movaps	 %xmm4, %xmm10
	unpcklpd %xmm5, %xmm4
	movaps	 %xmm6, %xmm11
	unpcklpd %xmm7, %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 + 32) * SIZE(B)
#endif

	movaps	%xmm0,   16 * SIZE(B)
	movaps	%xmm2,   18 * SIZE(B)
	movaps	%xmm4,   20 * SIZE(B)
	movaps	%xmm6,   22 * SIZE(B)

	unpckhpd %xmm1, %xmm8
	unpckhpd %xmm3, %xmm9
	unpckhpd %xmm5, %xmm10
	unpckhpd %xmm7, %xmm11

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 + 40) * SIZE(B)
#endif

	movaps	%xmm8,   24 * SIZE(B)
	movaps	%xmm9,   26 * SIZE(B)
	movaps	%xmm10,  28 * SIZE(B)
	movaps	%xmm11,  30 * SIZE(B)

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * SIZE(AO2, LDA, 2)
#endif

	MOVAPS	6 * SIZE(AO1),         %xmm0
	MOVAPS	6 * SIZE(AO1, LDA),    %xmm1
	MOVAPS	6 * SIZE(AO1, LDA, 2), %xmm2
	MOVAPS	6 * SIZE(AO1, LDA3),   %xmm3

	movaps	 %xmm0, %xmm8
	unpcklpd %xmm1, %xmm0
	movaps	 %xmm2, %xmm9
	unpcklpd %xmm3, %xmm2

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * SIZE(AO2, LDA3)
#endif

	MOVAPS	6 * SIZE(AO2),         %xmm4
	MOVAPS	6 * SIZE(AO2, LDA),    %xmm5
	MOVAPS	6 * SIZE(AO2, LDA, 2), %xmm6
	MOVAPS	6 * SIZE(AO2, LDA3),   %xmm7

	movaps	 %xmm4, %xmm10
	unpcklpd %xmm5, %xmm4
	movaps	 %xmm6, %xmm11
	unpcklpd %xmm7, %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 + 48) * SIZE(B)
#endif

	movaps	%xmm0,   32 * SIZE(B)
	movaps	%xmm2,   34 * SIZE(B)
	movaps	%xmm4,   36 * SIZE(B)
	movaps	%xmm6,   38 * SIZE(B)

	unpckhpd %xmm1, %xmm8
	unpckhpd %xmm3, %xmm9
	unpckhpd %xmm5, %xmm10
	unpckhpd %xmm7, %xmm11

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 + 56) * SIZE(B)
#endif

	movaps	%xmm8,   40 * SIZE(B)
	movaps	%xmm9,   42 * SIZE(B)
	movaps	%xmm10,  44 * SIZE(B)
	movaps	%xmm11,  46 * SIZE(B)

	addq	$8 * SIZE, AO1
	addq	$8 * SIZE, AO2
	subq	$-64 * SIZE, B

	decq	I
	jg	L(13)
	ALIGN_4

L(14):
	testq	$4, MM
	jle	L(16)

	MOVAPS	0 * SIZE(AO1),         %xmm0
	MOVAPS	0 * SIZE(AO1, LDA),    %xmm1
	MOVAPS	0 * SIZE(AO1, LDA, 2), %xmm2
	MOVAPS	0 * SIZE(AO1, LDA3),   %xmm3

	MOVAPS	0 * SIZE(AO2),         %xmm4
	MOVAPS	0 * SIZE(AO2, LDA),    %xmm5
	MOVAPS	0 * SIZE(AO2, LDA, 2), %xmm6
	MOVAPS	0 * SIZE(AO2, LDA3),   %xmm7

	movaps	 %xmm0, %xmm8
	unpcklpd %xmm1, %xmm0
	movaps	 %xmm2, %xmm9
	unpcklpd %xmm3, %xmm2

	movaps	 %xmm4, %xmm10
	unpcklpd %xmm5, %xmm4
	movaps	 %xmm6, %xmm11
	unpcklpd %xmm7, %xmm6

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm2,  -14 * SIZE(B)
	movaps	%xmm4,  -12 * SIZE(B)
	movaps	%xmm6,  -10 * SIZE(B)

	unpckhpd %xmm1, %xmm8
	unpckhpd %xmm3, %xmm9
	unpckhpd %xmm5, %xmm10
	unpckhpd %xmm7, %xmm11

	movaps	%xmm8,   -8 * SIZE(B)
	movaps	%xmm9,   -6 * SIZE(B)
	movaps	%xmm10,  -4 * SIZE(B)
	movaps	%xmm11,  -2 * SIZE(B)

	MOVAPS	2 * SIZE(AO1),         %xmm0
	MOVAPS	2 * SIZE(AO1, LDA),    %xmm1
	MOVAPS	2 * SIZE(AO1, LDA, 2), %xmm2
	MOVAPS	2 * SIZE(AO1, LDA3),   %xmm3

	MOVAPS	2 * SIZE(AO2),         %xmm4
	MOVAPS	2 * SIZE(AO2, LDA),    %xmm5
	MOVAPS	2 * SIZE(AO2, LDA, 2), %xmm6
	MOVAPS	2 * SIZE(AO2, LDA3),   %xmm7

	movaps	 %xmm0, %xmm8
	unpcklpd %xmm1, %xmm0
	movaps	 %xmm2, %xmm9
	unpcklpd %xmm3, %xmm2

	movaps	 %xmm4, %xmm10
	unpcklpd %xmm5, %xmm4
	movaps	 %xmm6, %xmm11
	unpcklpd %xmm7, %xmm6

	movaps	%xmm0,    0 * SIZE(B)
	movaps	%xmm2,    2 * SIZE(B)
	movaps	%xmm4,    4 * SIZE(B)
	movaps	%xmm6,    6 * SIZE(B)

	unpckhpd %xmm1, %xmm8
	unpckhpd %xmm3, %xmm9
	unpckhpd %xmm5, %xmm10
	unpckhpd %xmm7, %xmm11

	movaps	%xmm8,    8 * SIZE(B)
	movaps	%xmm9,   10 * SIZE(B)
	movaps	%xmm10,  12 * SIZE(B)
	movaps	%xmm11,  14 * SIZE(B)

	addq	$4 * SIZE, AO1
	addq	$4 * SIZE, AO2
	subq	$-32 * SIZE, B
	ALIGN_4

L(16):
	testq	$2, MM
	jle	L(18)

	MOVAPS	0 * SIZE(AO1),         %xmm0
	MOVAPS	0 * SIZE(AO1, LDA),    %xmm1
	MOVAPS	0 * SIZE(AO1, LDA, 2), %xmm2
	MOVAPS	0 * SIZE(AO1, LDA3),   %xmm3

	MOVAPS	0 * SIZE(AO2),         %xmm4
	MOVAPS	0 * SIZE(AO2, LDA),    %xmm5
	MOVAPS	0 * SIZE(AO2, LDA, 2), %xmm6
	MOVAPS	0 * SIZE(AO2, LDA3),   %xmm7

	movaps	 %xmm0, %xmm8
	unpcklpd %xmm1, %xmm0
	movaps	 %xmm2, %xmm9
	unpcklpd %xmm3, %xmm2

	movaps	 %xmm4, %xmm10
	unpcklpd %xmm5, %xmm4
	movaps	 %xmm6, %xmm11
	unpcklpd %xmm7, %xmm6

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm2,  -14 * SIZE(B)
	movaps	%xmm4,  -12 * SIZE(B)
	movaps	%xmm6,  -10 * SIZE(B)

	unpckhpd %xmm1, %xmm8
	unpckhpd %xmm3, %xmm9
	unpckhpd %xmm5, %xmm10
	unpckhpd %xmm7, %xmm11

	movaps	%xmm8,   -8 * SIZE(B)
	movaps	%xmm9,   -6 * SIZE(B)
	movaps	%xmm10,  -4 * SIZE(B)
	movaps	%xmm11,  -2 * SIZE(B)

	addq	$2 * SIZE, AO1
	addq	$2 * SIZE, AO2
	subq	$-16 * SIZE, B
	ALIGN_4

L(18):
	testq	$1, MM
	jle	L(19)

	movsd	0 * SIZE(AO1),         %xmm0
	movsd	0 * SIZE(AO1, LDA),    %xmm1
	movsd	0 * SIZE(AO1, LDA, 2), %xmm2
	movsd	0 * SIZE(AO1, LDA3),   %xmm3

	movsd	0 * SIZE(AO2),         %xmm4
	movsd	0 * SIZE(AO2, LDA),    %xmm5
	movsd	0 * SIZE(AO2, LDA, 2), %xmm6
	movsd	0 * SIZE(AO2, LDA3),   %xmm7

	unpcklpd %xmm1, %xmm0
	unpcklpd %xmm3, %xmm2
	unpcklpd %xmm5, %xmm4
	unpcklpd %xmm7, %xmm6

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm2,  -14 * SIZE(B)
	movaps	%xmm4,  -12 * SIZE(B)
	movaps	%xmm6,  -10 * SIZE(B)

	subq	$-8 * SIZE, B
	ALIGN_4

L(19):
	decq	J
	jg	L(11)
	ALIGN_4

L(20):
	testq	$4, N
	jle	L(30)

	movq	A, AO1
	leaq	(A, LDA, 2), AO2
	leaq	(A, LDA, 4), A

	testq	$SIZE, A
	je	L(22)

	movsd	0 * SIZE(AO1),      %xmm0
	movsd	0 * SIZE(AO1, LDA), %xmm1
	movsd	0 * SIZE(AO2),      %xmm2
	movsd	0 * SIZE(AO2, LDA), %xmm3

	unpcklpd %xmm1, %xmm0
	unpcklpd %xmm3, %xmm2

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm2,  -14 * SIZE(B)

	addq	$1 * SIZE, AO1
	addq	$1 * SIZE, AO2
	subq	$-4 * SIZE, B
	ALIGN_3

L(22):
	movq	MM, I
	sarq	$3, I
	jle	L(24)
	ALIGN_4

L(23):
#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * 2 * SIZE(AO1)
#endif

	MOVAPS	0 * SIZE(AO1),      %xmm0
	MOVAPS	0 * SIZE(AO1, LDA), %xmm1
	MOVAPS	0 * SIZE(AO2),      %xmm2
	MOVAPS	0 * SIZE(AO2, LDA), %xmm3

	movaps	 %xmm0, %xmm4
	unpcklpd %xmm1, %xmm0
	movaps	 %xmm2, %xmm6
	unpcklpd %xmm3, %xmm2

	unpckhpd %xmm1, %xmm4
	unpckhpd %xmm3, %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 +  0) * SIZE(B)
#endif

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm2,  -14 * SIZE(B)
	movaps	%xmm4,  -12 * SIZE(B)
	movaps	%xmm6,  -10 * SIZE(B)

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * 2 * SIZE(AO1, LDA)
#endif

	MOVAPS	2 * SIZE(AO1),      %xmm0
	MOVAPS	2 * SIZE(AO1, LDA), %xmm1
	MOVAPS	2 * SIZE(AO2),      %xmm2
	MOVAPS	2 * SIZE(AO2, LDA), %xmm3

	movaps	 %xmm0, %xmm4
	unpcklpd %xmm1, %xmm0
	movaps	 %xmm2, %xmm6
	unpcklpd %xmm3, %xmm2

	unpckhpd %xmm1, %xmm4
	unpckhpd %xmm3, %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 +  8) * SIZE(B)
#endif

	movaps	%xmm0,  -8 * SIZE(B)
	movaps	%xmm2,  -6 * SIZE(B)
	movaps	%xmm4,  -4 * SIZE(B)
	movaps	%xmm6,  -2 * SIZE(B)

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * 2 * SIZE(AO2)
#endif

	MOVAPS	4 * SIZE(AO1),      %xmm0
	MOVAPS	4 * SIZE(AO1, LDA), %xmm1
	MOVAPS	4 * SIZE(AO2),      %xmm2
	MOVAPS	4 * SIZE(AO2, LDA), %xmm3

	movaps	 %xmm0, %xmm4
	unpcklpd %xmm1, %xmm0
	movaps	 %xmm2, %xmm6
	unpcklpd %xmm3, %xmm2

	unpckhpd %xmm1, %xmm4
	unpckhpd %xmm3, %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 + 16) * SIZE(B)
#endif

	movaps	%xmm0,   0 * SIZE(B)
	movaps	%xmm2,   2 * SIZE(B)
	movaps	%xmm4,   4 * SIZE(B)
	movaps	%xmm6,   6 * SIZE(B)

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * 2 * SIZE(AO2, LDA)
#endif

	MOVAPS	6 * SIZE(AO1),      %xmm0
	MOVAPS	6 * SIZE(AO1, LDA), %xmm1
	MOVAPS	6 * SIZE(AO2),      %xmm2
	MOVAPS	6 * SIZE(AO2, LDA), %xmm3

	movaps	 %xmm0, %xmm4
	unpcklpd %xmm1, %xmm0
	movaps	 %xmm2, %xmm6
	unpcklpd %xmm3, %xmm2

	unpckhpd %xmm1, %xmm4
	unpckhpd %xmm3, %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 + 24) * SIZE(B)
#endif

	movaps	%xmm0,   8 * SIZE(B)
	movaps	%xmm2,  10 * SIZE(B)
	movaps	%xmm4,  12 * SIZE(B)
	movaps	%xmm6,  14 * SIZE(B)

	addq	$8 * SIZE, AO1
	addq	$8 * SIZE, AO2
	subq	$-32 * SIZE, B

	decq	I
	jg	L(23)
	ALIGN_4

L(24):
	testq	$4, MM
	jle	L(26)

	MOVAPS	0 * SIZE(AO1),      %xmm0
	MOVAPS	0 * SIZE(AO1, LDA), %xmm1
	MOVAPS	0 * SIZE(AO2),      %xmm2
	MOVAPS	0 * SIZE(AO2, LDA), %xmm3

	movaps	 %xmm0, %xmm4
	unpcklpd %xmm1, %xmm0
	movaps	 %xmm2, %xmm6
	unpcklpd %xmm3, %xmm2

	unpckhpd %xmm1, %xmm4
	unpckhpd %xmm3, %xmm6

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm2,  -14 * SIZE(B)
	movaps	%xmm4,  -12 * SIZE(B)
	movaps	%xmm6,  -10 * SIZE(B)

	MOVAPS	2 * SIZE(AO1),      %xmm0
	MOVAPS	2 * SIZE(AO1, LDA), %xmm1
	MOVAPS	2 * SIZE(AO2),      %xmm2
	MOVAPS	2 * SIZE(AO2, LDA), %xmm3

	movaps	 %xmm0, %xmm4
	unpcklpd %xmm1, %xmm0
	movaps	 %xmm2, %xmm6
	unpcklpd %xmm3, %xmm2

	unpckhpd %xmm1, %xmm4
	unpckhpd %xmm3, %xmm6

	movaps	%xmm0,  -8 * SIZE(B)
	movaps	%xmm2,  -6 * SIZE(B)
	movaps	%xmm4,  -4 * SIZE(B)
	movaps	%xmm6,  -2 * SIZE(B)

	addq	$4 * SIZE, AO1
	addq	$4 * SIZE, AO2
	subq	$-16 * SIZE, B
	ALIGN_4

L(26):
	testq	$2, MM
	jle	L(28)

	MOVAPS	0 * SIZE(AO1),      %xmm0
	MOVAPS	0 * SIZE(AO1, LDA), %xmm1
	MOVAPS	0 * SIZE(AO2),      %xmm2
	MOVAPS	0 * SIZE(AO2, LDA), %xmm3

	movaps	 %xmm0, %xmm4
	unpcklpd %xmm1, %xmm0
	movaps	 %xmm2, %xmm6
	unpcklpd %xmm3, %xmm2

	unpckhpd %xmm1, %xmm4
	unpckhpd %xmm3, %xmm6

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm2,  -14 * SIZE(B)
	movaps	%xmm4,  -12 * SIZE(B)
	movaps	%xmm6,  -10 * SIZE(B)

	addq	$2 * SIZE, AO1
	addq	$2 * SIZE, AO2
	subq	$-8 * SIZE, B
	ALIGN_4

L(28):
	testq	$1, MM
	jle	L(30)

	movsd	0 * SIZE(AO1),      %xmm0
	movsd	0 * SIZE(AO1, LDA), %xmm1
	movsd	0 * SIZE(AO2),      %xmm2
	movsd	0 * SIZE(AO2, LDA), %xmm3

	unpcklpd %xmm1, %xmm0
	unpcklpd %xmm3, %xmm2

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm2,  -14 * SIZE(B)
	subq	$-4 * SIZE, B
	ALIGN_4

L(30):
	testq	$2, N
	jle	L(40)

	movq	A, AO1
	leaq	(A, LDA), AO2
	leaq	(A, LDA, 2), A

	testq	$SIZE, A
	je	L(32)

	movsd	0 * SIZE(AO1),      %xmm0
	movsd	0 * SIZE(AO2),      %xmm1

	unpcklpd %xmm1, %xmm0

	movaps	%xmm0,  -16 * SIZE(B)

	addq	$1 * SIZE, AO1
	addq	$1 * SIZE, AO2
	subq	$-2 * SIZE, B
	ALIGN_3

L(32):
	movq	MM, I
	sarq	$3, I
	jle	L(34)
	ALIGN_4

L(33):
#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * 4 * SIZE(AO1)
#endif

	MOVAPS	0 * SIZE(AO1),      %xmm0
	MOVAPS	0 * SIZE(AO2),      %xmm1
	MOVAPS	2 * SIZE(AO1),      %xmm2
	MOVAPS	2 * SIZE(AO2),      %xmm3

	movaps	 %xmm0, %xmm4
	unpcklpd %xmm1, %xmm0
	movaps	 %xmm2, %xmm6
	unpcklpd %xmm3, %xmm2

	unpckhpd %xmm1, %xmm4
	unpckhpd %xmm3, %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 +  0) * SIZE(B)
#endif

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm4,  -14 * SIZE(B)
	movaps	%xmm2,  -12 * SIZE(B)
	movaps	%xmm6,  -10 * SIZE(B)

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * 4 * SIZE(AO2)
#endif

	MOVAPS	4 * SIZE(AO1), %xmm0
	MOVAPS	4 * SIZE(AO2), %xmm1
	MOVAPS	6 * SIZE(AO1), %xmm2
	MOVAPS	6 * SIZE(AO2), %xmm3

	movaps	 %xmm0, %xmm4
	unpcklpd %xmm1, %xmm0
	movaps	 %xmm2, %xmm6
	unpcklpd %xmm3, %xmm2

	unpckhpd %xmm1, %xmm4
	unpckhpd %xmm3, %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 +  8) * SIZE(B)
#endif

	movaps	%xmm0,  -8 * SIZE(B)
	movaps	%xmm4,  -6 * SIZE(B)
	movaps	%xmm2,  -4 * SIZE(B)
	movaps	%xmm6,  -2 * SIZE(B)

	addq	$8 * SIZE, AO1
	addq	$8 * SIZE, AO2
	subq	$-16 * SIZE, B

	decq	I
	jg	L(33)
	ALIGN_4

L(34):
	testq	$4, MM
	jle	L(36)

	MOVAPS	0 * SIZE(AO1),      %xmm0
	MOVAPS	0 * SIZE(AO2),      %xmm1
	MOVAPS	2 * SIZE(AO1),      %xmm2
	MOVAPS	2 * SIZE(AO2),      %xmm3

	movaps	 %xmm0, %xmm4
	unpcklpd %xmm1, %xmm0
	unpckhpd %xmm1, %xmm4

	movaps	 %xmm2, %xmm6
	unpcklpd %xmm3, %xmm2
	unpckhpd %xmm3, %xmm6

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm4,  -14 * SIZE(B)
	movaps	%xmm2,  -12 * SIZE(B)
	movaps	%xmm6,  -10 * SIZE(B)

	addq	$4 * SIZE, AO1
	addq	$4 * SIZE, AO2
	subq	$-8 * SIZE, B
	ALIGN_4

L(36):
	testq	$2, MM
	jle	L(38)

	MOVAPS	0 * SIZE(AO1),      %xmm0
	MOVAPS	0 * SIZE(AO2),      %xmm1

	movaps	 %xmm0, %xmm2
	unpcklpd %xmm1, %xmm0
	unpckhpd %xmm1, %xmm2

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm2,  -14 * SIZE(B)

	addq	$2 * SIZE, AO1
	addq	$2 * SIZE, AO2
	subq	$-4 * SIZE, B
	ALIGN_4

L(38):
	testq	$1, MM
	jle	L(40)

	movsd	0 * SIZE(AO1),      %xmm0
	movsd	0 * SIZE(AO2),      %xmm1

	unpcklpd %xmm1, %xmm0

	movaps	%xmm0,  -16 * SIZE(B)
	subq	$-2 * SIZE, B
	ALIGN_4

L(40):
	testq	$1, N
	jle	L(999)

	movq	A, AO1

	testq	$SIZE, A
	jne	L(45)

	movq	MM, I
	sarq	$3, I
	jle	L(42)
	ALIGN_4

L(41):
#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * 8 * SIZE(AO1)
#endif

	MOVAPS	0 * SIZE(AO1),      %xmm0
	MOVAPS	2 * SIZE(AO1),      %xmm1
	MOVAPS	4 * SIZE(AO1),      %xmm2
	MOVAPS	6 * SIZE(AO1),      %xmm3

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 +  0) * SIZE(B)
#endif

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm1,  -14 * SIZE(B)
	movaps	%xmm2,  -12 * SIZE(B)
	movaps	%xmm3,  -10 * SIZE(B)

	addq	$8 * SIZE, AO1
	subq	$-8 * SIZE, B

	decq	I
	jg	L(41)
	ALIGN_4

L(42):
	testq	$4, MM
	jle	L(43)

	MOVAPS	0 * SIZE(AO1),      %xmm0
	MOVAPS	2 * SIZE(AO1),      %xmm1

	movaps	%xmm0, -16 * SIZE(B)
	movaps	%xmm1, -14 * SIZE(B)

	addq	$4 * SIZE, AO1
	subq	$-4 * SIZE, B
	ALIGN_4

L(43):
	testq	$2, MM
	jle	L(44)

	MOVAPS	0 * SIZE(AO1),      %xmm0

	movaps	%xmm0,  -16 * SIZE(B)

	addq	$2 * SIZE, AO1
	subq	$-2 * SIZE, B
	ALIGN_4

L(44):
	testq	$1, MM
	jle	L(999)

	movsd	0 * SIZE(AO1),      %xmm0

	movlpd	%xmm0,  -16 * SIZE(B)
	jmp	L(999)
	ALIGN_4

L(45):
	MOVAPS	-1 * SIZE(AO1),      %xmm0

	movq	M, I
	sarq	$3, I
	jle	L(46)
	ALIGN_4

L(46):
#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * 8 * SIZE(AO1)
#endif

	MOVAPS	1 * SIZE(AO1),      %xmm1
	MOVAPS	3 * SIZE(AO1),      %xmm2
	MOVAPS	5 * SIZE(AO1),      %xmm3
	MOVAPS	7 * SIZE(AO1),      %xmm4

	shufpd	$1, %xmm1, %xmm0
	shufpd	$1, %xmm2, %xmm1
	shufpd	$1, %xmm3, %xmm2
	shufpd	$1, %xmm4, %xmm3

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 +  0) * SIZE(B)
#endif

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm1,  -14 * SIZE(B)
	movaps	%xmm2,  -12 * SIZE(B)
	movaps	%xmm3,  -10 * SIZE(B)

	movaps	%xmm4, %xmm0

	addq	$8 * SIZE, AO1
	subq	$-8 * SIZE, B

	decq	I
	jg	L(46)
	ALIGN_4

L(47):
	testq	$4, M
	jle	L(48)

	MOVAPS	1 * SIZE(AO1),      %xmm1
	MOVAPS	3 * SIZE(AO1),      %xmm2

	shufpd	$1, %xmm1, %xmm0
	shufpd	$1, %xmm2, %xmm1

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm1,  -14 * SIZE(B)

	movaps	%xmm2, %xmm0

	addq	$4 * SIZE, AO1
	addq	$4 * SIZE, B
	ALIGN_4

L(48):
	testq	$2, M
	jle	L(49)

	MOVAPS	1 * SIZE(AO1),      %xmm1

	shufpd	$1, %xmm1, %xmm0

	movaps	%xmm0,  -16 * SIZE(B)

	movaps	%xmm1, %xmm0

	addq	$2 * SIZE, AO1
	subq	$-2 * SIZE, B
	ALIGN_4

L(49):
	testq	$1, M
	jle	L(999)

	shufpd	$1, %xmm0, %xmm0

	movlpd	%xmm0,  -16 * SIZE(B)
	jmp	L(999)
	ALIGN_4

L(50):
	movq	N,  J
	sarq	$3, J
	jle	L(60)
	ALIGN_4

L(51):
	movq	A, AO1
	leaq	(A, LDA, 4), AO2
	leaq	(A, LDA, 8), A

	testq	$SIZE, A
	je	L(52)

	movsd	0 * SIZE(AO1),         %xmm0
	movsd	0 * SIZE(AO1, LDA),    %xmm1
	movsd	0 * SIZE(AO1, LDA, 2), %xmm2
	movsd	0 * SIZE(AO1, LDA3),   %xmm3
	movsd	0 * SIZE(AO2),         %xmm4
	movsd	0 * SIZE(AO2, LDA),    %xmm5
	movsd	0 * SIZE(AO2, LDA, 2), %xmm6
	movsd	0 * SIZE(AO2, LDA3),   %xmm7

	unpcklpd %xmm1, %xmm0
	unpcklpd %xmm3, %xmm2
	unpcklpd %xmm5, %xmm4
	unpcklpd %xmm7, %xmm6

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm2,  -14 * SIZE(B)
	movaps	%xmm4,  -12 * SIZE(B)
	movaps	%xmm6,  -10 * SIZE(B)

	addq	$1 * SIZE, AO1
	addq	$1 * SIZE, AO2
	subq	$-8 * SIZE, B
	ALIGN_3

L(52):
	MOVAPS	-1 * SIZE(AO1, LDA),  %xmm9
	MOVAPS	-1 * SIZE(AO1, LDA3), %xmm10
	MOVAPS	-1 * SIZE(AO2, LDA),  %xmm11
	MOVAPS	-1 * SIZE(AO2, LDA3), %xmm12

	movq	MM, I
	sarq	$3, I
	jle	L(54)
	ALIGN_4

L(53):
#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * SIZE(AO1)
#endif

	MOVAPS	0 * SIZE(AO1),         %xmm0
	MOVAPS	1 * SIZE(AO1, LDA),    %xmm1
	MOVAPS	0 * SIZE(AO1, LDA, 2), %xmm2
	MOVAPS	1 * SIZE(AO1, LDA3),   %xmm3

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * SIZE(AO1, LDA)
#endif

	MOVAPS	0 * SIZE(AO2),         %xmm4
	MOVAPS	1 * SIZE(AO2, LDA),    %xmm5
	MOVAPS	0 * SIZE(AO2, LDA, 2), %xmm6
	MOVAPS	1 * SIZE(AO2, LDA3),   %xmm7

	movsd	 %xmm0, %xmm9
	movsd	 %xmm2, %xmm10
	movsd	 %xmm4, %xmm11
	movsd	 %xmm6, %xmm12

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 +  0) * SIZE(B)
#endif

	movaps	%xmm9,  -16 * SIZE(B)
	movaps	%xmm10, -14 * SIZE(B)
	movaps	%xmm11, -12 * SIZE(B)
	movaps	%xmm12, -10 * SIZE(B)

	shufpd	 $1, %xmm1, %xmm0
	shufpd	 $1, %xmm3, %xmm2
	shufpd	 $1, %xmm5, %xmm4
	shufpd	 $1, %xmm7, %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 +  8) * SIZE(B)
#endif

	movaps	%xmm0,   -8 * SIZE(B)
	movaps	%xmm2,   -6 * SIZE(B)
	movaps	%xmm4,   -4 * SIZE(B)
	movaps	%xmm6,   -2 * SIZE(B)

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * SIZE(AO1, LDA, 2)
#endif

	MOVAPS	2 * SIZE(AO1),         %xmm0
	MOVAPS	3 * SIZE(AO1, LDA),    %xmm9
	MOVAPS	2 * SIZE(AO1, LDA, 2), %xmm2
	MOVAPS	3 * SIZE(AO1, LDA3),   %xmm10

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * SIZE(AO1, LDA3)
#endif

	MOVAPS	2 * SIZE(AO2),         %xmm4
	MOVAPS	3 * SIZE(AO2, LDA),    %xmm11
	MOVAPS	2 * SIZE(AO2, LDA, 2), %xmm6
	MOVAPS	3 * SIZE(AO2, LDA3),   %xmm12

	movsd	 %xmm0, %xmm1
	movsd	 %xmm2, %xmm3
	movsd	 %xmm4, %xmm5
	movsd	 %xmm6, %xmm7

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 + 16) * SIZE(B)
#endif

	movaps	%xmm1,    0 * SIZE(B)
	movaps	%xmm3,    2 * SIZE(B)
	movaps	%xmm5,    4 * SIZE(B)
	movaps	%xmm7,    6 * SIZE(B)

	shufpd	 $1, %xmm9,  %xmm0
	shufpd	 $1, %xmm10, %xmm2
	shufpd	 $1, %xmm11, %xmm4
	shufpd	 $1, %xmm12, %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 + 24) * SIZE(B)
#endif

	movaps	%xmm0,    8 * SIZE(B)
	movaps	%xmm2,   10 * SIZE(B)
	movaps	%xmm4,   12 * SIZE(B)
	movaps	%xmm6,   14 * SIZE(B)

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * SIZE(AO2)
#endif

	MOVAPS	4 * SIZE(AO1),         %xmm0
	MOVAPS	5 * SIZE(AO1, LDA),    %xmm1
	MOVAPS	4 * SIZE(AO1, LDA, 2), %xmm2
	MOVAPS	5 * SIZE(AO1, LDA3),   %xmm3

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * SIZE(AO2, LDA)
#endif

	MOVAPS	4 * SIZE(AO2),         %xmm4
	MOVAPS	5 * SIZE(AO2, LDA),    %xmm5
	MOVAPS	4 * SIZE(AO2, LDA, 2), %xmm6
	MOVAPS	5 * SIZE(AO2, LDA3),   %xmm7

	movsd	 %xmm0, %xmm9
	movsd	 %xmm2, %xmm10
	movsd	 %xmm4, %xmm11
	movsd	 %xmm6, %xmm12

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 + 32) * SIZE(B)
#endif

	movaps	%xmm9,   16 * SIZE(B)
	movaps	%xmm10,  18 * SIZE(B)
	movaps	%xmm11,  20 * SIZE(B)
	movaps	%xmm12,  22 * SIZE(B)

	shufpd	 $1, %xmm1, %xmm0
	shufpd	 $1, %xmm3, %xmm2
	shufpd	 $1, %xmm5, %xmm4
	shufpd	 $1, %xmm7, %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 4 +  8) * SIZE(B)
#endif

	movaps	%xmm0,   24 * SIZE(B)
	movaps	%xmm2,   26 * SIZE(B)
	movaps	%xmm4,   28 * SIZE(B)
	movaps	%xmm6,   30 * SIZE(B)

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * SIZE(AO2, LDA, 2)
#endif

	MOVAPS	6 * SIZE(AO1),         %xmm0
	MOVAPS	7 * SIZE(AO1, LDA),    %xmm9
	MOVAPS	6 * SIZE(AO1, LDA, 2), %xmm2
	MOVAPS	7 * SIZE(AO1, LDA3),   %xmm10

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * SIZE(AO2, LDA3)
#endif

	MOVAPS	6 * SIZE(AO2),         %xmm4
	MOVAPS	7 * SIZE(AO2, LDA),    %xmm11
	MOVAPS	6 * SIZE(AO2, LDA, 2), %xmm6
	MOVAPS	7 * SIZE(AO2, LDA3),   %xmm12

	movsd	 %xmm0, %xmm1
	movsd	 %xmm2, %xmm3
	movsd	 %xmm4, %xmm5
	movsd	 %xmm6, %xmm7

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 + 40) * SIZE(B)
#endif

	movaps	%xmm1,   32 * SIZE(B)
	movaps	%xmm3,   34 * SIZE(B)
	movaps	%xmm5,   36 * SIZE(B)
	movaps	%xmm7,   38 * SIZE(B)

	shufpd	 $1, %xmm9,  %xmm0
	shufpd	 $1, %xmm10, %xmm2
	shufpd	 $1, %xmm11, %xmm4
	shufpd	 $1, %xmm12, %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 + 48) * SIZE(B)
#endif
	movaps	%xmm0,   40 * SIZE(B)
	movaps	%xmm2,   42 * SIZE(B)
	movaps	%xmm4,   44 * SIZE(B)
	movaps	%xmm6,   46 * SIZE(B)

	addq	$8 * SIZE, AO1
	addq	$8 * SIZE, AO2
	subq	$-64 * SIZE, B

	decq	I
	jg	L(53)
	ALIGN_4

L(54):
	testq	$4, MM
	jle	L(56)

	MOVAPS	0 * SIZE(AO1),         %xmm0
	MOVAPS	1 * SIZE(AO1, LDA),    %xmm1
	MOVAPS	0 * SIZE(AO1, LDA, 2), %xmm2
	MOVAPS	1 * SIZE(AO1, LDA3),   %xmm3
	MOVAPS	0 * SIZE(AO2),         %xmm4
	MOVAPS	1 * SIZE(AO2, LDA),    %xmm5
	MOVAPS	0 * SIZE(AO2, LDA, 2), %xmm6
	MOVAPS	1 * SIZE(AO2, LDA3),   %xmm7

	movsd	 %xmm0, %xmm9
	movsd	 %xmm2, %xmm10
	movsd	 %xmm4, %xmm11
	movsd	 %xmm6, %xmm12

	movaps	%xmm9,  -16 * SIZE(B)
	movaps	%xmm10, -14 * SIZE(B)
	movaps	%xmm11, -12 * SIZE(B)
	movaps	%xmm12, -10 * SIZE(B)

	shufpd	 $1, %xmm1, %xmm0
	shufpd	 $1, %xmm3, %xmm2
	shufpd	 $1, %xmm5, %xmm4
	shufpd	 $1, %xmm7, %xmm6

	movaps	%xmm0,   -8 * SIZE(B)
	movaps	%xmm2,   -6 * SIZE(B)
	movaps	%xmm4,   -4 * SIZE(B)
	movaps	%xmm6,   -2 * SIZE(B)

	MOVAPS	2 * SIZE(AO1),         %xmm0
	MOVAPS	3 * SIZE(AO1, LDA),    %xmm9
	MOVAPS	2 * SIZE(AO1, LDA, 2), %xmm2
	MOVAPS	3 * SIZE(AO1, LDA3),   %xmm10
	MOVAPS	2 * SIZE(AO2),         %xmm4
	MOVAPS	3 * SIZE(AO2, LDA),    %xmm11
	MOVAPS	2 * SIZE(AO2, LDA, 2), %xmm6
	MOVAPS	3 * SIZE(AO2, LDA3),   %xmm12

	movsd	 %xmm0, %xmm1
	movsd	 %xmm2, %xmm3
	movsd	 %xmm4, %xmm5
	movsd	 %xmm6, %xmm7

	movaps	%xmm1,    0 * SIZE(B)
	movaps	%xmm3,    2 * SIZE(B)
	movaps	%xmm5,    4 * SIZE(B)
	movaps	%xmm7,    6 * SIZE(B)

	shufpd	 $1, %xmm9,  %xmm0
	shufpd	 $1, %xmm10, %xmm2
	shufpd	 $1, %xmm11, %xmm4
	shufpd	 $1, %xmm12, %xmm6

	movaps	%xmm0,    8 * SIZE(B)
	movaps	%xmm2,   10 * SIZE(B)
	movaps	%xmm4,   12 * SIZE(B)
	movaps	%xmm6,   14 * SIZE(B)

	addq	$4 * SIZE, AO1
	addq	$4 * SIZE, AO2
	subq	$-32 * SIZE, B
	ALIGN_4

L(56):
	testq	$2, MM
	jle	L(58)

	MOVAPS	0 * SIZE(AO1),         %xmm0
	MOVAPS	1 * SIZE(AO1, LDA),    %xmm1
	MOVAPS	0 * SIZE(AO1, LDA, 2), %xmm2
	MOVAPS	1 * SIZE(AO1, LDA3),   %xmm3
	MOVAPS	0 * SIZE(AO2),         %xmm4
	MOVAPS	1 * SIZE(AO2, LDA),    %xmm5
	MOVAPS	0 * SIZE(AO2, LDA, 2), %xmm6
	MOVAPS	1 * SIZE(AO2, LDA3),   %xmm7

	movsd	 %xmm0, %xmm9
	movsd	 %xmm2, %xmm10
	movsd	 %xmm4, %xmm11
	movsd	 %xmm6, %xmm12

	movaps	%xmm9,  -16 * SIZE(B)
	movaps	%xmm10, -14 * SIZE(B)
	movaps	%xmm11, -12 * SIZE(B)
	movaps	%xmm12, -10 * SIZE(B)

	shufpd	 $1, %xmm1, %xmm0
	shufpd	 $1, %xmm3, %xmm2
	shufpd	 $1, %xmm5, %xmm4
	shufpd	 $1, %xmm7, %xmm6

	movaps	%xmm0,   -8 * SIZE(B)
	movaps	%xmm2,   -6 * SIZE(B)
	movaps	%xmm4,   -4 * SIZE(B)
	movaps	%xmm6,   -2 * SIZE(B)

	addq	$2 * SIZE, AO1
	addq	$2 * SIZE, AO2
	subq	$-16 * SIZE, B
	ALIGN_4

L(58):
	testq	$1, MM
	jle	L(59)

	movsd	0 * SIZE(AO1),         %xmm0
	movsd	0 * SIZE(AO1, LDA),    %xmm1
	movsd	0 * SIZE(AO1, LDA, 2), %xmm2
	movsd	0 * SIZE(AO1, LDA3),   %xmm3
	movsd	0 * SIZE(AO2),         %xmm4
	movsd	0 * SIZE(AO2, LDA),    %xmm5
	movsd	0 * SIZE(AO2, LDA, 2), %xmm6
	movsd	0 * SIZE(AO2, LDA3),   %xmm7

	unpcklpd %xmm1, %xmm0
	unpcklpd %xmm3, %xmm2
	unpcklpd %xmm5, %xmm4
	unpcklpd %xmm7, %xmm6

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm2,  -14 * SIZE(B)
	movaps	%xmm4,  -12 * SIZE(B)
	movaps	%xmm6,  -10 * SIZE(B)

	subq	$-8 * SIZE, B
	ALIGN_4

L(59):
	decq	J
	jg	L(51)
	ALIGN_4

L(60):
	testq	$4, N
	jle	L(70)

	movq	A, AO1
	leaq	(A, LDA, 2), AO2
	leaq	(A, LDA, 4), A

	testq	$SIZE, A
	je	L(62)

	movsd	0 * SIZE(AO1),      %xmm0
	movsd	0 * SIZE(AO1, LDA), %xmm1
	movsd	0 * SIZE(AO2),      %xmm2
	movsd	0 * SIZE(AO2, LDA), %xmm3

	unpcklpd %xmm1, %xmm0
	unpcklpd %xmm3, %xmm2

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm2,  -14 * SIZE(B)

	addq	$1 * SIZE, AO1
	addq	$1 * SIZE, AO2
	subq	$-4 * SIZE, B
	ALIGN_3

L(62):
	movaps	-1 * SIZE(AO1, LDA), %xmm5
	movaps	-1 * SIZE(AO2, LDA), %xmm7

	movq	MM, I
	sarq	$3, I
	jle	L(64)
	ALIGN_4

L(63):
#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * 2 * SIZE(AO1)
#endif

	MOVAPS	0 * SIZE(AO1),      %xmm0
	MOVAPS	1 * SIZE(AO1, LDA), %xmm1
	MOVAPS	0 * SIZE(AO2),      %xmm2
	MOVAPS	1 * SIZE(AO2, LDA), %xmm3

	movsd	 %xmm0, %xmm5
	movsd	 %xmm2, %xmm7
	shufpd	 $1, %xmm1, %xmm0
	shufpd	 $1, %xmm3, %xmm2

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 +  0) * SIZE(B)
#endif

	movaps	%xmm5,  -16 * SIZE(B)
	movaps	%xmm7,  -14 * SIZE(B)
	movaps	%xmm0,  -12 * SIZE(B)
	movaps	%xmm2,  -10 * SIZE(B)

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * 2 * SIZE(AO1, LDA)
#endif

	MOVAPS	2 * SIZE(AO1),      %xmm0
	MOVAPS	3 * SIZE(AO1, LDA), %xmm5
	MOVAPS	2 * SIZE(AO2),      %xmm2
	MOVAPS	3 * SIZE(AO2, LDA), %xmm7

	movsd	 %xmm0, %xmm1
	movsd	 %xmm2, %xmm3
	shufpd	 $1, %xmm5, %xmm0
	shufpd	 $1, %xmm7, %xmm2

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 +  8) * SIZE(B)
#endif

	movaps	%xmm1,  -8 * SIZE(B)
	movaps	%xmm3,  -6 * SIZE(B)
	movaps	%xmm0,  -4 * SIZE(B)
	movaps	%xmm2,  -2 * SIZE(B)

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * 2 * SIZE(AO2)
#endif

	MOVAPS	4 * SIZE(AO1),      %xmm0
	MOVAPS	5 * SIZE(AO1, LDA), %xmm1
	MOVAPS	4 * SIZE(AO2),      %xmm2
	MOVAPS	5 * SIZE(AO2, LDA), %xmm3

	movsd	 %xmm0, %xmm5
	movsd	 %xmm2, %xmm7
	shufpd	 $1, %xmm1, %xmm0
	shufpd	 $1, %xmm3, %xmm2

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 + 16) * SIZE(B)
#endif

	movaps	%xmm5,    0 * SIZE(B)
	movaps	%xmm7,    2 * SIZE(B)
	movaps	%xmm0,    4 * SIZE(B)
	movaps	%xmm2,    6 * SIZE(B)

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * 2 * SIZE(AO2, LDA)
#endif

	MOVAPS	6 * SIZE(AO1),      %xmm0
	MOVAPS	7 * SIZE(AO1, LDA), %xmm5
	MOVAPS	6 * SIZE(AO2),      %xmm2
	MOVAPS	7 * SIZE(AO2, LDA), %xmm7

	movsd	 %xmm0, %xmm1
	movsd	 %xmm2, %xmm3
	shufpd	 $1, %xmm5, %xmm0
	shufpd	 $1, %xmm7, %xmm2

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 + 24) * SIZE(B)
#endif

	movaps	%xmm1,   8 * SIZE(B)
	movaps	%xmm3,  10 * SIZE(B)
	movaps	%xmm0,  12 * SIZE(B)
	movaps	%xmm2,  14 * SIZE(B)

	addq	$8 * SIZE, AO1
	addq	$8 * SIZE, AO2
	subq	$-32 * SIZE, B

	decq	I
	jg	L(63)
	ALIGN_4

L(64):
	testq	$4, MM
	jle	L(66)

	MOVAPS	0 * SIZE(AO1),      %xmm0
	MOVAPS	1 * SIZE(AO1, LDA), %xmm1
	MOVAPS	0 * SIZE(AO2),      %xmm2
	MOVAPS	1 * SIZE(AO2, LDA), %xmm3

	movsd	 %xmm0, %xmm5
	shufpd	 $1, %xmm1, %xmm0
	movsd	 %xmm2, %xmm7
	shufpd	 $1, %xmm3, %xmm2

	movaps	%xmm5,  -16 * SIZE(B)
	movaps	%xmm7,  -14 * SIZE(B)
	movaps	%xmm0,  -12 * SIZE(B)
	movaps	%xmm2,  -10 * SIZE(B)

	MOVAPS	2 * SIZE(AO1),      %xmm0
	MOVAPS	3 * SIZE(AO1, LDA), %xmm5
	MOVAPS	2 * SIZE(AO2),      %xmm2
	MOVAPS	3 * SIZE(AO2, LDA), %xmm7

	movsd	 %xmm0, %xmm1
	shufpd	 $1, %xmm5, %xmm0
	movsd	 %xmm2, %xmm3
	shufpd	 $1, %xmm7, %xmm2

	movaps	%xmm1,  -8 * SIZE(B)
	movaps	%xmm3,  -6 * SIZE(B)
	movaps	%xmm0,  -4 * SIZE(B)
	movaps	%xmm2,  -2 * SIZE(B)

	addq	$4 * SIZE, AO1
	addq	$4 * SIZE, AO2
	subq	$-16 * SIZE, B
	ALIGN_4

L(66):
	testq	$2, MM
	jle	L(68)

	MOVAPS	0 * SIZE(AO1),      %xmm0
	MOVAPS	1 * SIZE(AO1, LDA), %xmm1
	MOVAPS	0 * SIZE(AO2),      %xmm2
	MOVAPS	1 * SIZE(AO2, LDA), %xmm3

	movsd	 %xmm0, %xmm5
	movsd	 %xmm2, %xmm7
	shufpd	 $1, %xmm1, %xmm0
	shufpd	 $1, %xmm3, %xmm2

	movaps	%xmm5,  -16 * SIZE(B)
	movaps	%xmm7,  -14 * SIZE(B)
	movaps	%xmm0,  -12 * SIZE(B)
	movaps	%xmm2,  -10 * SIZE(B)

	addq	$2 * SIZE, AO1
	addq	$2 * SIZE, AO2
	subq	$-8 * SIZE, B
	ALIGN_4

L(68):
	testq	$1, MM
	jle	L(70)

	movsd	0 * SIZE(AO1),      %xmm0
	movsd	0 * SIZE(AO1, LDA), %xmm1
	movsd	0 * SIZE(AO2),      %xmm2
	movsd	0 * SIZE(AO2, LDA), %xmm3

	unpcklpd %xmm1, %xmm0
	unpcklpd %xmm3, %xmm2

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm2,  -14 * SIZE(B)
	subq	$-4 * SIZE, B
	ALIGN_4

L(70):
	testq	$2, N
	jle	L(80)

	movq	A, AO1
	leaq	(A, LDA), AO2
	leaq	(A, LDA, 2), A

	testq	$SIZE, A
	je	L(72)

	movsd	0 * SIZE(AO1), %xmm0
	movsd	0 * SIZE(AO2), %xmm1

	unpcklpd %xmm1, %xmm0

	movaps	%xmm0,  -16 * SIZE(B)

	addq	$1 * SIZE, AO1
	addq	$1 * SIZE, AO2
	subq	$-2 * SIZE, B
	ALIGN_3

L(72):
	MOVAPS	-1 * SIZE(AO2), %xmm5

	movq	MM, I
	sarq	$3, I
	jle	L(74)
	ALIGN_4

L(73):
#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * 4 * SIZE(AO1)
#endif

	MOVAPS	0 * SIZE(AO1), %xmm0
	MOVAPS	1 * SIZE(AO2), %xmm1
	MOVAPS	2 * SIZE(AO1), %xmm2
	MOVAPS	3 * SIZE(AO2), %xmm3

	movsd	 %xmm0, %xmm5
	shufpd	 $1, %xmm1, %xmm0
	movsd	 %xmm2, %xmm1
	shufpd	 $1, %xmm3, %xmm2

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 +  0) * SIZE(B)
#endif

	movaps	%xmm5,  -16 * SIZE(B)
	movaps	%xmm0,  -14 * SIZE(B)
	movaps	%xmm1,  -12 * SIZE(B)
	movaps	%xmm2,  -10 * SIZE(B)

#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * 4 * SIZE(AO2)
#endif

	MOVAPS	4 * SIZE(AO1), %xmm0
	MOVAPS	5 * SIZE(AO2), %xmm1
	MOVAPS	6 * SIZE(AO1), %xmm2
	MOVAPS	7 * SIZE(AO2), %xmm5

	movsd	 %xmm0, %xmm3
	shufpd	 $1, %xmm1, %xmm0
	movsd	 %xmm2, %xmm1
	shufpd	 $1, %xmm5, %xmm2

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 +  8) * SIZE(B)
#endif

	movaps	%xmm3,   -8 * SIZE(B)
	movaps	%xmm0,   -6 * SIZE(B)
	movaps	%xmm1,   -4 * SIZE(B)
	movaps	%xmm2,   -2 * SIZE(B)

	addq	$8 * SIZE, AO1
	addq	$8 * SIZE, AO2
	subq	$-16 * SIZE, B

	decq	I
	jg	L(73)
	ALIGN_4

L(74):
	testq	$4, MM
	jle	L(76)

	MOVAPS	0 * SIZE(AO1), %xmm0
	MOVAPS	1 * SIZE(AO2), %xmm1
	MOVAPS	2 * SIZE(AO1), %xmm2
	MOVAPS	3 * SIZE(AO2), %xmm3

	movsd	 %xmm0, %xmm5
	shufpd	 $1, %xmm1, %xmm0
	movsd	 %xmm2, %xmm1
	shufpd	 $1, %xmm3, %xmm2

	movaps	%xmm5,  -16 * SIZE(B)
	movaps	%xmm0,  -14 * SIZE(B)
	movaps	%xmm1,  -12 * SIZE(B)
	movaps	%xmm2,  -10 * SIZE(B)

	movaps	%xmm3, %xmm5

	addq	$4 * SIZE, AO1
	addq	$4 * SIZE, AO2
	subq	$-8 * SIZE, B
	ALIGN_4

L(76):
	testq	$2, MM
	jle	L(78)

	MOVAPS	0 * SIZE(AO1), %xmm0
	MOVAPS	1 * SIZE(AO2), %xmm1

	movsd	 %xmm0, %xmm5
	shufpd	 $1, %xmm1, %xmm0

	movaps	%xmm5,  -16 * SIZE(B)
	movaps	%xmm0,  -14 * SIZE(B)

	addq	$2 * SIZE, AO1
	addq	$2 * SIZE, AO2
	subq	$-4 * SIZE, B
	ALIGN_4

L(78):
	testq	$1, MM
	jle	L(80)

	movsd	0 * SIZE(AO1),      %xmm0
	movsd	0 * SIZE(AO2),      %xmm1

	unpcklpd %xmm1, %xmm0

	movaps	%xmm0,  -16 * SIZE(B)
	subq	$-2 * SIZE, B
	ALIGN_4

L(80):
	testq	$1, N
	jle	L(999)

	movq	A, AO1

	testq	$SIZE, A
	jne	L(85)

	movq	MM, I
	sarq	$3, I
	jle	L(82)
	ALIGN_4

L(81):
#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * 8 * SIZE(AO1)
#endif

	MOVAPS	0 * SIZE(AO1),      %xmm0
	MOVAPS	2 * SIZE(AO1),      %xmm2
	MOVAPS	4 * SIZE(AO1),      %xmm4
	MOVAPS	6 * SIZE(AO1),      %xmm6

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 +  0) * SIZE(B)
#endif

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm2,  -14 * SIZE(B)
	movaps	%xmm4,  -12 * SIZE(B)
	movaps	%xmm6,  -10 * SIZE(B)

	addq	$8 * SIZE, AO1
	subq	$-8 * SIZE, B

	decq	I
	jg	L(81)
	ALIGN_4

L(82):
	testq	$4, MM
	jle	L(83)

	MOVAPS	0 * SIZE(AO1),      %xmm0
	MOVAPS	2 * SIZE(AO1),      %xmm2

	movaps	%xmm0, -16 * SIZE(B)
	movaps	%xmm2, -14 * SIZE(B)

	addq	$4 * SIZE, AO1
	subq	$-4 * SIZE, B
	ALIGN_4

L(83):
	testq	$2, MM
	jle	L(84)

	MOVAPS	0 * SIZE(AO1),      %xmm0

	movaps	%xmm0,  -16 * SIZE(B)

	addq	$2 * SIZE, AO1
	subq	$-2 * SIZE, B
	ALIGN_4

L(84):
	testq	$1, MM
	jle	L(999)

	movsd	0 * SIZE(AO1),      %xmm0

	movlpd	%xmm0,  -16 * SIZE(B)
	jmp	L(999)
	ALIGN_4

L(85):
	MOVAPS	-1 * SIZE(AO1),      %xmm0

	movq	M, I
	sarq	$3, I
	jle	L(86)
	ALIGN_4

L(86):
#ifdef PREFETCH
	PREFETCH	PREFETCHSIZE * 8 * SIZE(AO1)
#endif

	MOVAPS	1 * SIZE(AO1),      %xmm1
	MOVAPS	3 * SIZE(AO1),      %xmm2
	MOVAPS	5 * SIZE(AO1),      %xmm3
	MOVAPS	7 * SIZE(AO1),      %xmm4

	shufpd	$1, %xmm1, %xmm0
	shufpd	$1, %xmm2, %xmm1
	shufpd	$1, %xmm3, %xmm2
	shufpd	$1, %xmm4, %xmm3

#ifdef PREFETCHW
	PREFETCHW	(PREFETCHSIZE * 8 +  0) * SIZE(B)
#endif

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm1,  -14 * SIZE(B)
	movaps	%xmm2,  -12 * SIZE(B)
	movaps	%xmm3,  -10 * SIZE(B)

	movaps	%xmm4, %xmm0

	addq	$8 * SIZE, AO1
	subq	$-8 * SIZE, B

	decq	I
	jg	L(86)
	ALIGN_4

L(87):
	testq	$4, M
	jle	L(88)

	MOVAPS	1 * SIZE(AO1),      %xmm1
	MOVAPS	3 * SIZE(AO1),      %xmm2

	shufpd	$1, %xmm1, %xmm0
	shufpd	$1, %xmm2, %xmm1

	movaps	%xmm0,  -16 * SIZE(B)
	movaps	%xmm1,  -14 * SIZE(B)

	movaps	%xmm2, %xmm0

	addq	$4 * SIZE, AO1
	addq	$4 * SIZE, B
	ALIGN_4

L(88):
	testq	$2, M
	jle	L(89)

	MOVAPS	1 * SIZE(AO1),      %xmm1

	shufpd	$1, %xmm1, %xmm0

	movaps	%xmm0,  -16 * SIZE(B)

	movaps	%xmm1, %xmm0

	addq	$2 * SIZE, AO1
	subq	$-2 * SIZE, B
	ALIGN_4

L(89):
	testq	$1, M
	jle	L(999)

	shufpd	$1, %xmm0, %xmm0

	movlpd	%xmm0,  -16 * SIZE(B)
	ALIGN_4

L(999):
#ifdef WINDOWS_ABI
	movups	  0(%rsp), %xmm6
	movups	 16(%rsp), %xmm7
	movups	 32(%rsp), %xmm8
	movups	 48(%rsp), %xmm9
	movups	 64(%rsp), %xmm10
	movups	 80(%rsp), %xmm11
	movups	 96(%rsp), %xmm12

	addq	$STACKSIZE, %rsp
#endif

	popq	%r12
	popq	%r13

#ifdef WINDOWS_ABI
	popq	%r14
	popq	%r15
#endif
	ret

	EPILOGUE
