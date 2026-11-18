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

#define A_PRE 256
#define B_PRE 128


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

	vmovups	%xmm6,    0(%rsp)
	vmovups	%xmm7,   16(%rsp)
	vmovups	%xmm8,   32(%rsp)
	vmovups	%xmm9,   48(%rsp)
	vmovups	%xmm10,  64(%rsp)
	vmovups	%xmm11,  80(%rsp)
	vmovups	%xmm12,  96(%rsp)

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

	vmovsd	0 * SIZE(AO1),         %xmm0
	vmovsd	0 * SIZE(AO1, LDA),    %xmm1
	vmovsd	0 * SIZE(AO1, LDA, 2), %xmm2
	vmovsd	0 * SIZE(AO1, LDA3),   %xmm3

	vmovsd	0 * SIZE(AO2),         %xmm4
	vmovsd	0 * SIZE(AO2, LDA),    %xmm5
	vmovsd	0 * SIZE(AO2, LDA, 2), %xmm6
	vmovsd	0 * SIZE(AO2, LDA3),   %xmm7

	vunpcklpd %xmm1, %xmm0 , %xmm0
	vunpcklpd %xmm3, %xmm2 , %xmm2
	vunpcklpd %xmm5, %xmm4 , %xmm4
	vunpcklpd %xmm7, %xmm6 , %xmm6

	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm2,  -14 * SIZE(B)
	vmovups	%xmm4,  -12 * SIZE(B)
	vmovups	%xmm6,  -10 * SIZE(B)

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

	prefetchnta	A_PRE(AO1)
	vmovups	0 * SIZE(AO1),         %xmm0
	vmovups	0 * SIZE(AO1, LDA),    %xmm1
	prefetchnta	A_PRE(AO1, LDA)
	vmovups	0 * SIZE(AO1, LDA, 2), %xmm2
	vmovups	0 * SIZE(AO1, LDA3),   %xmm3

	prefetchnta	A_PRE(AO1, LDA, 2)
	vmovups	0 * SIZE(AO1, LDA, 2), %xmm2
	vmovups	 %xmm0, %xmm8
	prefetchnta	A_PRE(AO1, LDA3)
	vunpcklpd %xmm1, %xmm0 , %xmm0
	vmovups	 %xmm2, %xmm9
	vunpcklpd %xmm3, %xmm2 , %xmm2


	prefetchnta	A_PRE(AO2)
	vmovups	0 * SIZE(AO2),         %xmm4
	vmovups	0 * SIZE(AO2, LDA),    %xmm5
	prefetchnta	A_PRE(AO2, LDA)
	vmovups	0 * SIZE(AO2, LDA, 2), %xmm6
	vmovups	0 * SIZE(AO2, LDA3),   %xmm7

	prefetchnta	A_PRE(AO2, LDA, 2)
	vmovups	 %xmm4, %xmm10
	vunpcklpd %xmm5, %xmm4 , %xmm4
	prefetchnta	A_PRE(AO2, LDA3)
	vmovups	 %xmm6, %xmm11
	vunpcklpd %xmm7, %xmm6 , %xmm6

	prefetchw	B_PRE(B)
	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm2,  -14 * SIZE(B)
	vmovups	%xmm4,  -12 * SIZE(B)
	vmovups	%xmm6,  -10 * SIZE(B)

	vunpckhpd %xmm1, %xmm8 , %xmm8
	vunpckhpd %xmm3, %xmm9 , %xmm9
	vunpckhpd %xmm5, %xmm10, %xmm10
	vunpckhpd %xmm7, %xmm11, %xmm11


	prefetchw	B_PRE+64(B)
	vmovups	%xmm8,   -8 * SIZE(B)
	vmovups	%xmm9,   -6 * SIZE(B)
	vmovups	%xmm10,  -4 * SIZE(B)
	vmovups	%xmm11,  -2 * SIZE(B)

/***********************************************************************************************/

	vmovups	2 * SIZE(AO1),         %xmm0
	vmovups	2 * SIZE(AO1, LDA),    %xmm1
	vmovups	2 * SIZE(AO1, LDA, 2), %xmm2
	vmovups	2 * SIZE(AO1, LDA3),   %xmm3

	vmovups	 %xmm0, %xmm8
	vunpcklpd %xmm1, %xmm0 , %xmm0
	vmovups	 %xmm2, %xmm9
	vunpcklpd %xmm3, %xmm2 , %xmm2


	vmovups	2 * SIZE(AO2),         %xmm4
	vmovups	2 * SIZE(AO2, LDA),    %xmm5
	vmovups	2 * SIZE(AO2, LDA, 2), %xmm6
	vmovups	2 * SIZE(AO2, LDA3),   %xmm7

	vmovups	 %xmm4, %xmm10
	vunpcklpd %xmm5, %xmm4 , %xmm4
	vmovups	 %xmm6, %xmm11
	vunpcklpd %xmm7, %xmm6 , %xmm6


	prefetchw	B_PRE+128(B)
	vmovups	%xmm0,    0 * SIZE(B)
	vmovups	%xmm2,    2 * SIZE(B)
	vmovups	%xmm4,    4 * SIZE(B)
	vmovups	%xmm6,    6 * SIZE(B)

	vunpckhpd %xmm1, %xmm8 , %xmm8
	vunpckhpd %xmm3, %xmm9 , %xmm9
	vunpckhpd %xmm5, %xmm10, %xmm10
	vunpckhpd %xmm7, %xmm11, %xmm11


	prefetchw	B_PRE+192(B)
	vmovups	%xmm8,    8 * SIZE(B)
	vmovups	%xmm9,   10 * SIZE(B)
	vmovups	%xmm10,  12 * SIZE(B)
	vmovups	%xmm11,  14 * SIZE(B)

/***********************************************************************************************/

	vmovups	4 * SIZE(AO1),         %xmm0
	vmovups	4 * SIZE(AO1, LDA),    %xmm1
	vmovups	4 * SIZE(AO1, LDA, 2), %xmm2
	vmovups	4 * SIZE(AO1, LDA3),   %xmm3

	vmovups	 %xmm0, %xmm8
	vunpcklpd %xmm1, %xmm0 , %xmm0
	vmovups	 %xmm2, %xmm9
	vunpcklpd %xmm3, %xmm2 , %xmm2


	vmovups	4 * SIZE(AO2),         %xmm4
	vmovups	4 * SIZE(AO2, LDA),    %xmm5
	vmovups	4 * SIZE(AO2, LDA, 2), %xmm6
	vmovups	4 * SIZE(AO2, LDA3),   %xmm7

	vmovups	 %xmm4, %xmm10
	vunpcklpd %xmm5, %xmm4 , %xmm4
	vmovups	 %xmm6, %xmm11
	vunpcklpd %xmm7, %xmm6 , %xmm6


	prefetchw	B_PRE+256(B)
	vmovups	%xmm0,   16 * SIZE(B)
	vmovups	%xmm2,   18 * SIZE(B)
	vmovups	%xmm4,   20 * SIZE(B)
	vmovups	%xmm6,   22 * SIZE(B)

	vunpckhpd %xmm1, %xmm8 , %xmm8
	vunpckhpd %xmm3, %xmm9 , %xmm9
	vunpckhpd %xmm5, %xmm10, %xmm10
	vunpckhpd %xmm7, %xmm11, %xmm11

	prefetchw	B_PRE+320(B)
	vmovups	%xmm8,   24 * SIZE(B)
	vmovups	%xmm9,   26 * SIZE(B)
	vmovups	%xmm10,  28 * SIZE(B)
	vmovups	%xmm11,  30 * SIZE(B)

/***********************************************************************************************/

	vmovups	6 * SIZE(AO1),         %xmm0
	vmovups	6 * SIZE(AO1, LDA),    %xmm1
	vmovups	6 * SIZE(AO1, LDA, 2), %xmm2
	vmovups	6 * SIZE(AO1, LDA3),   %xmm3

	vmovups	 %xmm0, %xmm8
	vunpcklpd %xmm1, %xmm0 , %xmm0
	vmovups	 %xmm2, %xmm9
	vunpcklpd %xmm3, %xmm2 , %xmm2


	vmovups	6 * SIZE(AO2),         %xmm4
	vmovups	6 * SIZE(AO2, LDA),    %xmm5
	vmovups	6 * SIZE(AO2, LDA, 2), %xmm6
	vmovups	6 * SIZE(AO2, LDA3),   %xmm7

	vmovups	 %xmm4, %xmm10
	vunpcklpd %xmm5, %xmm4 , %xmm4
	vmovups	 %xmm6, %xmm11
	vunpcklpd %xmm7, %xmm6 , %xmm6


	prefetchw	B_PRE+384(B)
	vmovups	%xmm0,   32 * SIZE(B)
	vmovups	%xmm2,   34 * SIZE(B)
	vmovups	%xmm4,   36 * SIZE(B)
	vmovups	%xmm6,   38 * SIZE(B)

	vunpckhpd %xmm1, %xmm8 , %xmm8
	vunpckhpd %xmm3, %xmm9 , %xmm9
	vunpckhpd %xmm5, %xmm10, %xmm10
	vunpckhpd %xmm7, %xmm11, %xmm11

	prefetchw	B_PRE+448(B)
	vmovups	%xmm8,   40 * SIZE(B)
	vmovups	%xmm9,   42 * SIZE(B)
	vmovups	%xmm10,  44 * SIZE(B)
	vmovups	%xmm11,  46 * SIZE(B)

	addq	$8 * SIZE, AO1
	addq	$8 * SIZE, AO2
	subq	$-64 * SIZE, B

	decq	I
	jg	L(13)
	ALIGN_4

L(14):
	testq	$4, MM
	jle	L(16)

	vmovups	0 * SIZE(AO1),         %xmm0
	vmovups	0 * SIZE(AO1, LDA),    %xmm1
	vmovups	0 * SIZE(AO1, LDA, 2), %xmm2
	vmovups	0 * SIZE(AO1, LDA3),   %xmm3

	vmovups	0 * SIZE(AO2),         %xmm4
	vmovups	0 * SIZE(AO2, LDA),    %xmm5
	vmovups	0 * SIZE(AO2, LDA, 2), %xmm6
	vmovups	0 * SIZE(AO2, LDA3),   %xmm7

	vmovups	 %xmm0, %xmm8
	vunpcklpd %xmm1, %xmm0 , %xmm0
	vmovups	 %xmm2, %xmm9
	vunpcklpd %xmm3, %xmm2 , %xmm2

	vmovups	 %xmm4, %xmm10
	vunpcklpd %xmm5, %xmm4 , %xmm4
	vmovups	 %xmm6, %xmm11
	vunpcklpd %xmm7, %xmm6 , %xmm6

	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm2,  -14 * SIZE(B)
	vmovups	%xmm4,  -12 * SIZE(B)
	vmovups	%xmm6,  -10 * SIZE(B)

	vunpckhpd %xmm1, %xmm8 , %xmm8
	vunpckhpd %xmm3, %xmm9 , %xmm9
	vunpckhpd %xmm5, %xmm10, %xmm10
	vunpckhpd %xmm7, %xmm11, %xmm11

	vmovups	%xmm8,   -8 * SIZE(B)
	vmovups	%xmm9,   -6 * SIZE(B)
	vmovups	%xmm10,  -4 * SIZE(B)
	vmovups	%xmm11,  -2 * SIZE(B)

	vmovups	2 * SIZE(AO1),         %xmm0
	vmovups	2 * SIZE(AO1, LDA),    %xmm1
	vmovups	2 * SIZE(AO1, LDA, 2), %xmm2
	vmovups	2 * SIZE(AO1, LDA3),   %xmm3

	vmovups	2 * SIZE(AO2),         %xmm4
	vmovups	2 * SIZE(AO2, LDA),    %xmm5
	vmovups	2 * SIZE(AO2, LDA, 2), %xmm6
	vmovups	2 * SIZE(AO2, LDA3),   %xmm7

	vmovups	 %xmm0, %xmm8
	vunpcklpd %xmm1, %xmm0 , %xmm0
	vmovups	 %xmm2, %xmm9
	vunpcklpd %xmm3, %xmm2 , %xmm2

	vmovups	 %xmm4, %xmm10
	vunpcklpd %xmm5, %xmm4 , %xmm4
	vmovups	 %xmm6, %xmm11
	vunpcklpd %xmm7, %xmm6 , %xmm6

	vmovups	%xmm0,    0 * SIZE(B)
	vmovups	%xmm2,    2 * SIZE(B)
	vmovups	%xmm4,    4 * SIZE(B)
	vmovups	%xmm6,    6 * SIZE(B)

	vunpckhpd %xmm1, %xmm8 , %xmm8
	vunpckhpd %xmm3, %xmm9 , %xmm9
	vunpckhpd %xmm5, %xmm10, %xmm10
	vunpckhpd %xmm7, %xmm11, %xmm11

	vmovups	%xmm8,    8 * SIZE(B)
	vmovups	%xmm9,   10 * SIZE(B)
	vmovups	%xmm10,  12 * SIZE(B)
	vmovups	%xmm11,  14 * SIZE(B)

	addq	$4 * SIZE, AO1
	addq	$4 * SIZE, AO2
	subq	$-32 * SIZE, B
	ALIGN_4

L(16):
	testq	$2, MM
	jle	L(18)

	vmovups	0 * SIZE(AO1),         %xmm0
	vmovups	0 * SIZE(AO1, LDA),    %xmm1
	vmovups	0 * SIZE(AO1, LDA, 2), %xmm2
	vmovups	0 * SIZE(AO1, LDA3),   %xmm3

	vmovups	0 * SIZE(AO2),         %xmm4
	vmovups	0 * SIZE(AO2, LDA),    %xmm5
	vmovups	0 * SIZE(AO2, LDA, 2), %xmm6
	vmovups	0 * SIZE(AO2, LDA3),   %xmm7

	vmovups	 %xmm0, %xmm8
	vunpcklpd %xmm1, %xmm0 , %xmm0
	vmovups	 %xmm2, %xmm9
	vunpcklpd %xmm3, %xmm2 , %xmm2

	vmovups	 %xmm4, %xmm10
	vunpcklpd %xmm5, %xmm4 , %xmm4
	vmovups	 %xmm6, %xmm11
	vunpcklpd %xmm7, %xmm6 , %xmm6

	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm2,  -14 * SIZE(B)
	vmovups	%xmm4,  -12 * SIZE(B)
	vmovups	%xmm6,  -10 * SIZE(B)

	vunpckhpd %xmm1, %xmm8 , %xmm8
	vunpckhpd %xmm3, %xmm9 , %xmm9
	vunpckhpd %xmm5, %xmm10, %xmm10
	vunpckhpd %xmm7, %xmm11, %xmm11

	vmovups	%xmm8,   -8 * SIZE(B)
	vmovups	%xmm9,   -6 * SIZE(B)
	vmovups	%xmm10,  -4 * SIZE(B)
	vmovups	%xmm11,  -2 * SIZE(B)

	addq	$2 * SIZE, AO1
	addq	$2 * SIZE, AO2
	subq	$-16 * SIZE, B
	ALIGN_4

L(18):
	testq	$1, MM
	jle	L(19)

	vmovsd	0 * SIZE(AO1),         %xmm0
	vmovsd	0 * SIZE(AO1, LDA),    %xmm1
	vmovsd	0 * SIZE(AO1, LDA, 2), %xmm2
	vmovsd	0 * SIZE(AO1, LDA3),   %xmm3

	vmovsd	0 * SIZE(AO2),         %xmm4
	vmovsd	0 * SIZE(AO2, LDA),    %xmm5
	vmovsd	0 * SIZE(AO2, LDA, 2), %xmm6
	vmovsd	0 * SIZE(AO2, LDA3),   %xmm7

	vunpcklpd %xmm1, %xmm0 , %xmm0
	vunpcklpd %xmm3, %xmm2 , %xmm2
	vunpcklpd %xmm5, %xmm4 , %xmm4
	vunpcklpd %xmm7, %xmm6 , %xmm6

	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm2,  -14 * SIZE(B)
	vmovups	%xmm4,  -12 * SIZE(B)
	vmovups	%xmm6,  -10 * SIZE(B)

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

	vmovsd	0 * SIZE(AO1),      %xmm0
	vmovsd	0 * SIZE(AO1, LDA), %xmm1
	vmovsd	0 * SIZE(AO2),      %xmm2
	vmovsd	0 * SIZE(AO2, LDA), %xmm3

	vunpcklpd %xmm1, %xmm0 , %xmm0
	vunpcklpd %xmm3, %xmm2 , %xmm2

	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm2,  -14 * SIZE(B)

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

	vmovups	0 * SIZE(AO1),      %xmm0
	vmovups	0 * SIZE(AO1, LDA), %xmm1
	vmovups	0 * SIZE(AO2),      %xmm2
	vmovups	0 * SIZE(AO2, LDA), %xmm3

	vmovups	 %xmm0, %xmm4
	vunpcklpd %xmm1, %xmm0 , %xmm0
	vmovups	 %xmm2, %xmm6
	vunpcklpd %xmm3, %xmm2 , %xmm2

	vunpckhpd %xmm1, %xmm4 , %xmm4
	vunpckhpd %xmm3, %xmm6 , %xmm6


	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm2,  -14 * SIZE(B)
	vmovups	%xmm4,  -12 * SIZE(B)
	vmovups	%xmm6,  -10 * SIZE(B)


	vmovups	2 * SIZE(AO1),      %xmm0
	vmovups	2 * SIZE(AO1, LDA), %xmm1
	vmovups	2 * SIZE(AO2),      %xmm2
	vmovups	2 * SIZE(AO2, LDA), %xmm3

	vmovups	 %xmm0, %xmm4
	vunpcklpd %xmm1, %xmm0 , %xmm0
	vmovups	 %xmm2, %xmm6
	vunpcklpd %xmm3, %xmm2 , %xmm2

	vunpckhpd %xmm1, %xmm4 , %xmm4
	vunpckhpd %xmm3, %xmm6 , %xmm6


	vmovups	%xmm0,  -8 * SIZE(B)
	vmovups	%xmm2,  -6 * SIZE(B)
	vmovups	%xmm4,  -4 * SIZE(B)
	vmovups	%xmm6,  -2 * SIZE(B)


	vmovups	4 * SIZE(AO1),      %xmm0
	vmovups	4 * SIZE(AO1, LDA), %xmm1
	vmovups	4 * SIZE(AO2),      %xmm2
	vmovups	4 * SIZE(AO2, LDA), %xmm3

	vmovups	 %xmm0, %xmm4
	vunpcklpd %xmm1, %xmm0 , %xmm0
	vmovups	 %xmm2, %xmm6
	vunpcklpd %xmm3, %xmm2 , %xmm2

	vunpckhpd %xmm1, %xmm4 , %xmm4
	vunpckhpd %xmm3, %xmm6 , %xmm6


	vmovups	%xmm0,   0 * SIZE(B)
	vmovups	%xmm2,   2 * SIZE(B)
	vmovups	%xmm4,   4 * SIZE(B)
	vmovups	%xmm6,   6 * SIZE(B)


	vmovups	6 * SIZE(AO1),      %xmm0
	vmovups	6 * SIZE(AO1, LDA), %xmm1
	vmovups	6 * SIZE(AO2),      %xmm2
	vmovups	6 * SIZE(AO2, LDA), %xmm3

	vmovups	 %xmm0, %xmm4
	vunpcklpd %xmm1, %xmm0 , %xmm0
	vmovups	 %xmm2, %xmm6
	vunpcklpd %xmm3, %xmm2 , %xmm2

	vunpckhpd %xmm1, %xmm4 , %xmm4
	vunpckhpd %xmm3, %xmm6 , %xmm6


	vmovups	%xmm0,   8 * SIZE(B)
	vmovups	%xmm2,  10 * SIZE(B)
	vmovups	%xmm4,  12 * SIZE(B)
	vmovups	%xmm6,  14 * SIZE(B)

	addq	$8 * SIZE, AO1
	addq	$8 * SIZE, AO2
	subq	$-32 * SIZE, B

	decq	I
	jg	L(23)
	ALIGN_4

L(24):
	testq	$4, MM
	jle	L(26)

	vmovups	0 * SIZE(AO1),      %xmm0
	vmovups	0 * SIZE(AO1, LDA), %xmm1
	vmovups	0 * SIZE(AO2),      %xmm2
	vmovups	0 * SIZE(AO2, LDA), %xmm3

	vmovups	 %xmm0, %xmm4
	vunpcklpd %xmm1, %xmm0 , %xmm0
	vmovups	 %xmm2, %xmm6
	vunpcklpd %xmm3, %xmm2 , %xmm2

	vunpckhpd %xmm1, %xmm4 , %xmm4
	vunpckhpd %xmm3, %xmm6 , %xmm6

	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm2,  -14 * SIZE(B)
	vmovups	%xmm4,  -12 * SIZE(B)
	vmovups	%xmm6,  -10 * SIZE(B)

	vmovups	2 * SIZE(AO1),      %xmm0
	vmovups	2 * SIZE(AO1, LDA), %xmm1
	vmovups	2 * SIZE(AO2),      %xmm2
	vmovups	2 * SIZE(AO2, LDA), %xmm3

	vmovups	 %xmm0, %xmm4
	vunpcklpd %xmm1, %xmm0 , %xmm0
	vmovups	 %xmm2, %xmm6
	vunpcklpd %xmm3, %xmm2 , %xmm2

	vunpckhpd %xmm1, %xmm4 , %xmm4
	vunpckhpd %xmm3, %xmm6 , %xmm6

	vmovups	%xmm0,  -8 * SIZE(B)
	vmovups	%xmm2,  -6 * SIZE(B)
	vmovups	%xmm4,  -4 * SIZE(B)
	vmovups	%xmm6,  -2 * SIZE(B)

	addq	$4 * SIZE, AO1
	addq	$4 * SIZE, AO2
	subq	$-16 * SIZE, B
	ALIGN_4

L(26):
	testq	$2, MM
	jle	L(28)

	vmovups	0 * SIZE(AO1),      %xmm0
	vmovups	0 * SIZE(AO1, LDA), %xmm1
	vmovups	0 * SIZE(AO2),      %xmm2
	vmovups	0 * SIZE(AO2, LDA), %xmm3

	vmovups	 %xmm0, %xmm4
	vunpcklpd %xmm1, %xmm0 , %xmm0
	vmovups	 %xmm2, %xmm6
	vunpcklpd %xmm3, %xmm2 , %xmm2

	vunpckhpd %xmm1, %xmm4 , %xmm4
	vunpckhpd %xmm3, %xmm6 , %xmm6

	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm2,  -14 * SIZE(B)
	vmovups	%xmm4,  -12 * SIZE(B)
	vmovups	%xmm6,  -10 * SIZE(B)

	addq	$2 * SIZE, AO1
	addq	$2 * SIZE, AO2
	subq	$-8 * SIZE, B
	ALIGN_4

L(28):
	testq	$1, MM
	jle	L(30)

	vmovsd	0 * SIZE(AO1),      %xmm0
	vmovsd	0 * SIZE(AO1, LDA), %xmm1
	vmovsd	0 * SIZE(AO2),      %xmm2
	vmovsd	0 * SIZE(AO2, LDA), %xmm3

	vunpcklpd %xmm1, %xmm0 , %xmm0
	vunpcklpd %xmm3, %xmm2 , %xmm2

	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm2,  -14 * SIZE(B)
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

	vmovsd	0 * SIZE(AO1),      %xmm0
	vmovsd	0 * SIZE(AO2),      %xmm1

	vunpcklpd %xmm1, %xmm0 , %xmm0

	vmovups	%xmm0,  -16 * SIZE(B)

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

	vmovups	0 * SIZE(AO1),      %xmm0
	vmovups	0 * SIZE(AO2),      %xmm1
	vmovups	2 * SIZE(AO1),      %xmm2
	vmovups	2 * SIZE(AO2),      %xmm3

	vmovups	 %xmm0, %xmm4
	vunpcklpd %xmm1, %xmm0 , %xmm0
	vmovups	 %xmm2, %xmm6
	vunpcklpd %xmm3, %xmm2 , %xmm2

	vunpckhpd %xmm1, %xmm4 , %xmm4
	vunpckhpd %xmm3, %xmm6 , %xmm6


	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm4,  -14 * SIZE(B)
	vmovups	%xmm2,  -12 * SIZE(B)
	vmovups	%xmm6,  -10 * SIZE(B)


	vmovups	4 * SIZE(AO1), %xmm0
	vmovups	4 * SIZE(AO2), %xmm1
	vmovups	6 * SIZE(AO1), %xmm2
	vmovups	6 * SIZE(AO2), %xmm3

	vmovups	 %xmm0, %xmm4
	vunpcklpd %xmm1, %xmm0 , %xmm0
	vmovups	 %xmm2, %xmm6
	vunpcklpd %xmm3, %xmm2 , %xmm2

	vunpckhpd %xmm1, %xmm4 , %xmm4
	vunpckhpd %xmm3, %xmm6 , %xmm6


	vmovups	%xmm0,  -8 * SIZE(B)
	vmovups	%xmm4,  -6 * SIZE(B)
	vmovups	%xmm2,  -4 * SIZE(B)
	vmovups	%xmm6,  -2 * SIZE(B)

	addq	$8 * SIZE, AO1
	addq	$8 * SIZE, AO2
	subq	$-16 * SIZE, B

	decq	I
	jg	L(33)
	ALIGN_4

L(34):
	testq	$4, MM
	jle	L(36)

	vmovups	0 * SIZE(AO1),      %xmm0
	vmovups	0 * SIZE(AO2),      %xmm1
	vmovups	2 * SIZE(AO1),      %xmm2
	vmovups	2 * SIZE(AO2),      %xmm3

	vmovups	 %xmm0, %xmm4
	vunpcklpd %xmm1, %xmm0 , %xmm0
	vunpckhpd %xmm1, %xmm4 , %xmm4

	vmovups	 %xmm2, %xmm6
	vunpcklpd %xmm3, %xmm2 , %xmm2
	vunpckhpd %xmm3, %xmm6 , %xmm6

	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm4,  -14 * SIZE(B)
	vmovups	%xmm2,  -12 * SIZE(B)
	vmovups	%xmm6,  -10 * SIZE(B)

	addq	$4 * SIZE, AO1
	addq	$4 * SIZE, AO2
	subq	$-8 * SIZE, B
	ALIGN_4

L(36):
	testq	$2, MM
	jle	L(38)

	vmovups	0 * SIZE(AO1),      %xmm0
	vmovups	0 * SIZE(AO2),      %xmm1

	vmovups	 %xmm0, %xmm2
	vunpcklpd %xmm1, %xmm0 , %xmm0
	vunpckhpd %xmm1, %xmm2 , %xmm2

	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm2,  -14 * SIZE(B)

	addq	$2 * SIZE, AO1
	addq	$2 * SIZE, AO2
	subq	$-4 * SIZE, B
	ALIGN_4

L(38):
	testq	$1, MM
	jle	L(40)

	vmovsd	0 * SIZE(AO1),      %xmm0
	vmovsd	0 * SIZE(AO2),      %xmm1

	vunpcklpd %xmm1, %xmm0 , %xmm0

	vmovups	%xmm0,  -16 * SIZE(B)
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

	vmovups	0 * SIZE(AO1),      %xmm0
	vmovups	2 * SIZE(AO1),      %xmm1
	vmovups	4 * SIZE(AO1),      %xmm2
	vmovups	6 * SIZE(AO1),      %xmm3


	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm1,  -14 * SIZE(B)
	vmovups	%xmm2,  -12 * SIZE(B)
	vmovups	%xmm3,  -10 * SIZE(B)

	addq	$8 * SIZE, AO1
	subq	$-8 * SIZE, B

	decq	I
	jg	L(41)
	ALIGN_4

L(42):
	testq	$4, MM
	jle	L(43)

	vmovups	0 * SIZE(AO1),      %xmm0
	vmovups	2 * SIZE(AO1),      %xmm1

	vmovups	%xmm0, -16 * SIZE(B)
	vmovups	%xmm1, -14 * SIZE(B)

	addq	$4 * SIZE, AO1
	subq	$-4 * SIZE, B
	ALIGN_4

L(43):
	testq	$2, MM
	jle	L(44)

	vmovups	0 * SIZE(AO1),      %xmm0

	vmovups	%xmm0,  -16 * SIZE(B)

	addq	$2 * SIZE, AO1
	subq	$-2 * SIZE, B
	ALIGN_4

L(44):
	testq	$1, MM
	jle	L(999)

	vmovsd	0 * SIZE(AO1),      %xmm0

	vmovlpd	%xmm0,  -16 * SIZE(B)
	jmp	L(999)
	ALIGN_4

L(45):
	vmovups	-1 * SIZE(AO1),      %xmm0

	movq	M, I
	sarq	$3, I
	jle	L(46)
	ALIGN_4

L(46):

	vmovups	1 * SIZE(AO1),      %xmm1
	vmovups	3 * SIZE(AO1),      %xmm2
	vmovups	5 * SIZE(AO1),      %xmm3
	vmovups	7 * SIZE(AO1),      %xmm4

	vshufpd	$1, %xmm1, %xmm0 , %xmm0
	vshufpd	$1, %xmm2, %xmm1 , %xmm1
	vshufpd	$1, %xmm3, %xmm2 , %xmm2
	vshufpd	$1, %xmm4, %xmm3 , %xmm3


	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm1,  -14 * SIZE(B)
	vmovups	%xmm2,  -12 * SIZE(B)
	vmovups	%xmm3,  -10 * SIZE(B)

	vmovups	%xmm4, %xmm0

	addq	$8 * SIZE, AO1
	subq	$-8 * SIZE, B

	decq	I
	jg	L(46)
	ALIGN_4

L(47):
	testq	$4, M
	jle	L(48)

	vmovups	1 * SIZE(AO1),      %xmm1
	vmovups	3 * SIZE(AO1),      %xmm2

	vshufpd	$1, %xmm1, %xmm0 , %xmm0
	vshufpd	$1, %xmm2, %xmm1 , %xmm1

	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm1,  -14 * SIZE(B)

	vmovups	%xmm2, %xmm0

	addq	$4 * SIZE, AO1
	addq	$4 * SIZE, B
	ALIGN_4

L(48):
	testq	$2, M
	jle	L(49)

	vmovups	1 * SIZE(AO1),      %xmm1

	vshufpd	$1, %xmm1, %xmm0 , %xmm0

	vmovups	%xmm0,  -16 * SIZE(B)

	vmovups	%xmm1, %xmm0

	addq	$2 * SIZE, AO1
	subq	$-2 * SIZE, B
	ALIGN_4

L(49):
	testq	$1, M
	jle	L(999)

	vshufpd	$1, %xmm0, %xmm0 , %xmm0

	vmovlpd	%xmm0,  -16 * SIZE(B)
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

	vmovsd	0 * SIZE(AO1),         %xmm0
	vmovsd	0 * SIZE(AO1, LDA),    %xmm1
	vmovsd	0 * SIZE(AO1, LDA, 2), %xmm2
	vmovsd	0 * SIZE(AO1, LDA3),   %xmm3
	vmovsd	0 * SIZE(AO2),         %xmm4
	vmovsd	0 * SIZE(AO2, LDA),    %xmm5
	vmovsd	0 * SIZE(AO2, LDA, 2), %xmm6
	vmovsd	0 * SIZE(AO2, LDA3),   %xmm7

	vunpcklpd %xmm1, %xmm0 , %xmm0
	vunpcklpd %xmm3, %xmm2 , %xmm2
	vunpcklpd %xmm5, %xmm4 , %xmm4
	vunpcklpd %xmm7, %xmm6 , %xmm6

	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm2,  -14 * SIZE(B)
	vmovups	%xmm4,  -12 * SIZE(B)
	vmovups	%xmm6,  -10 * SIZE(B)

	addq	$1 * SIZE, AO1
	addq	$1 * SIZE, AO2
	subq	$-8 * SIZE, B
	ALIGN_3

L(52):
	vmovups	-1 * SIZE(AO1, LDA),  %xmm9
	vmovups	-1 * SIZE(AO1, LDA3), %xmm10
	vmovups	-1 * SIZE(AO2, LDA),  %xmm11
	vmovups	-1 * SIZE(AO2, LDA3), %xmm12

	movq	MM, I
	sarq	$3, I
	jle	L(54)
	ALIGN_4

L(53):

	vmovups	0 * SIZE(AO1),         %xmm0
	vmovups	1 * SIZE(AO1, LDA),    %xmm1
	vmovups	0 * SIZE(AO1, LDA, 2), %xmm2
	vmovups	1 * SIZE(AO1, LDA3),   %xmm3


	vmovups	0 * SIZE(AO2),         %xmm4
	vmovups	1 * SIZE(AO2, LDA),    %xmm5
	vmovups	0 * SIZE(AO2, LDA, 2), %xmm6
	vmovups	1 * SIZE(AO2, LDA3),   %xmm7

	vmovsd %xmm0, %xmm9 , %xmm9
	vmovsd %xmm2, %xmm10, %xmm10
	vmovsd %xmm4, %xmm11, %xmm11
	vmovsd %xmm6, %xmm12, %xmm12


	vmovups	%xmm9,  -16 * SIZE(B)
	vmovups	%xmm10, -14 * SIZE(B)
	vmovups	%xmm11, -12 * SIZE(B)
	vmovups	%xmm12, -10 * SIZE(B)

	vshufpd	 $1, %xmm1, %xmm0 , %xmm0
	vshufpd	 $1, %xmm3, %xmm2 , %xmm2
	vshufpd	 $1, %xmm5, %xmm4 , %xmm4
	vshufpd	 $1, %xmm7, %xmm6 , %xmm6


	vmovups	%xmm0,   -8 * SIZE(B)
	vmovups	%xmm2,   -6 * SIZE(B)
	vmovups	%xmm4,   -4 * SIZE(B)
	vmovups	%xmm6,   -2 * SIZE(B)


	vmovups	2 * SIZE(AO1),         %xmm0
	vmovups	3 * SIZE(AO1, LDA),    %xmm9
	vmovups	2 * SIZE(AO1, LDA, 2), %xmm2
	vmovups	3 * SIZE(AO1, LDA3),   %xmm10


	vmovups	2 * SIZE(AO2),         %xmm4
	vmovups	3 * SIZE(AO2, LDA),    %xmm11
	vmovups	2 * SIZE(AO2, LDA, 2), %xmm6
	vmovups	3 * SIZE(AO2, LDA3),   %xmm12

	vmovsd	 %xmm0, %xmm1 , %xmm1
	vmovsd	 %xmm2, %xmm3 , %xmm3
	vmovsd	 %xmm4, %xmm5 , %xmm5
	vmovsd	 %xmm6, %xmm7 , %xmm7


	vmovups	%xmm1,    0 * SIZE(B)
	vmovups	%xmm3,    2 * SIZE(B)
	vmovups	%xmm5,    4 * SIZE(B)
	vmovups	%xmm7,    6 * SIZE(B)

	vshufpd	 $1, %xmm9,  %xmm0 , %xmm0
	vshufpd	 $1, %xmm10, %xmm2 , %xmm2
	vshufpd	 $1, %xmm11, %xmm4 , %xmm4
	vshufpd	 $1, %xmm12, %xmm6 , %xmm6


	vmovups	%xmm0,    8 * SIZE(B)
	vmovups	%xmm2,   10 * SIZE(B)
	vmovups	%xmm4,   12 * SIZE(B)
	vmovups	%xmm6,   14 * SIZE(B)


	vmovups	4 * SIZE(AO1),         %xmm0
	vmovups	5 * SIZE(AO1, LDA),    %xmm1
	vmovups	4 * SIZE(AO1, LDA, 2), %xmm2
	vmovups	5 * SIZE(AO1, LDA3),   %xmm3


	vmovups	4 * SIZE(AO2),         %xmm4
	vmovups	5 * SIZE(AO2, LDA),    %xmm5
	vmovups	4 * SIZE(AO2, LDA, 2), %xmm6
	vmovups	5 * SIZE(AO2, LDA3),   %xmm7

	vmovsd	 %xmm0, %xmm9 , %xmm9
	vmovsd	 %xmm2, %xmm10, %xmm10
	vmovsd	 %xmm4, %xmm11, %xmm11
	vmovsd	 %xmm6, %xmm12, %xmm12


	vmovups	%xmm9,   16 * SIZE(B)
	vmovups	%xmm10,  18 * SIZE(B)
	vmovups	%xmm11,  20 * SIZE(B)
	vmovups	%xmm12,  22 * SIZE(B)

	vshufpd	 $1, %xmm1, %xmm0 , %xmm0
	vshufpd	 $1, %xmm3, %xmm2 , %xmm2
	vshufpd	 $1, %xmm5, %xmm4 , %xmm4
	vshufpd	 $1, %xmm7, %xmm6 , %xmm6


	vmovups	%xmm0,   24 * SIZE(B)
	vmovups	%xmm2,   26 * SIZE(B)
	vmovups	%xmm4,   28 * SIZE(B)
	vmovups	%xmm6,   30 * SIZE(B)


	vmovups	6 * SIZE(AO1),         %xmm0
	vmovups	7 * SIZE(AO1, LDA),    %xmm9
	vmovups	6 * SIZE(AO1, LDA, 2), %xmm2
	vmovups	7 * SIZE(AO1, LDA3),   %xmm10


	vmovups	6 * SIZE(AO2),         %xmm4
	vmovups	7 * SIZE(AO2, LDA),    %xmm11
	vmovups	6 * SIZE(AO2, LDA, 2), %xmm6
	vmovups	7 * SIZE(AO2, LDA3),   %xmm12

	vmovsd	 %xmm0, %xmm1 , %xmm1
	vmovsd	 %xmm2, %xmm3 , %xmm3
	vmovsd	 %xmm4, %xmm5 , %xmm5
	vmovsd	 %xmm6, %xmm7 , %xmm7


	vmovups	%xmm1,   32 * SIZE(B)
	vmovups	%xmm3,   34 * SIZE(B)
	vmovups	%xmm5,   36 * SIZE(B)
	vmovups	%xmm7,   38 * SIZE(B)

	vshufpd	 $1, %xmm9,  %xmm0 , %xmm0
	vshufpd	 $1, %xmm10, %xmm2 , %xmm2
	vshufpd	 $1, %xmm11, %xmm4 , %xmm4
	vshufpd	 $1, %xmm12, %xmm6 , %xmm6

	vmovups	%xmm0,   40 * SIZE(B)
	vmovups	%xmm2,   42 * SIZE(B)
	vmovups	%xmm4,   44 * SIZE(B)
	vmovups	%xmm6,   46 * SIZE(B)

	addq	$8 * SIZE, AO1
	addq	$8 * SIZE, AO2
	subq	$-64 * SIZE, B

	decq	I
	jg	L(53)
	ALIGN_4

L(54):
	testq	$4, MM
	jle	L(56)

	vmovups	0 * SIZE(AO1),         %xmm0
	vmovups	1 * SIZE(AO1, LDA),    %xmm1
	vmovups	0 * SIZE(AO1, LDA, 2), %xmm2
	vmovups	1 * SIZE(AO1, LDA3),   %xmm3
	vmovups	0 * SIZE(AO2),         %xmm4
	vmovups	1 * SIZE(AO2, LDA),    %xmm5
	vmovups	0 * SIZE(AO2, LDA, 2), %xmm6
	vmovups	1 * SIZE(AO2, LDA3),   %xmm7

	vmovsd	 %xmm0, %xmm9 , %xmm9
	vmovsd	 %xmm2, %xmm10, %xmm10
	vmovsd	 %xmm4, %xmm11, %xmm11
	vmovsd	 %xmm6, %xmm12, %xmm12

	vmovups	%xmm9,  -16 * SIZE(B)
	vmovups	%xmm10, -14 * SIZE(B)
	vmovups	%xmm11, -12 * SIZE(B)
	vmovups	%xmm12, -10 * SIZE(B)

	vshufpd	 $1, %xmm1, %xmm0 , %xmm0
	vshufpd	 $1, %xmm3, %xmm2 , %xmm2
	vshufpd	 $1, %xmm5, %xmm4 , %xmm4
	vshufpd	 $1, %xmm7, %xmm6 , %xmm6

	vmovups	%xmm0,   -8 * SIZE(B)
	vmovups	%xmm2,   -6 * SIZE(B)
	vmovups	%xmm4,   -4 * SIZE(B)
	vmovups	%xmm6,   -2 * SIZE(B)

	vmovups	2 * SIZE(AO1),         %xmm0
	vmovups	3 * SIZE(AO1, LDA),    %xmm9
	vmovups	2 * SIZE(AO1, LDA, 2), %xmm2
	vmovups	3 * SIZE(AO1, LDA3),   %xmm10
	vmovups	2 * SIZE(AO2),         %xmm4
	vmovups	3 * SIZE(AO2, LDA),    %xmm11
	vmovups	2 * SIZE(AO2, LDA, 2), %xmm6
	vmovups	3 * SIZE(AO2, LDA3),   %xmm12

	vmovsd	 %xmm0, %xmm1 , %xmm1
	vmovsd	 %xmm2, %xmm3 , %xmm3
	vmovsd	 %xmm4, %xmm5 , %xmm5
	vmovsd	 %xmm6, %xmm7 , %xmm7

	vmovups	%xmm1,    0 * SIZE(B)
	vmovups	%xmm3,    2 * SIZE(B)
	vmovups	%xmm5,    4 * SIZE(B)
	vmovups	%xmm7,    6 * SIZE(B)

	vshufpd	 $1, %xmm9,  %xmm0 , %xmm0
	vshufpd	 $1, %xmm10, %xmm2 , %xmm2
	vshufpd	 $1, %xmm11, %xmm4 , %xmm4
	vshufpd	 $1, %xmm12, %xmm6 , %xmm6

	vmovups	%xmm0,    8 * SIZE(B)
	vmovups	%xmm2,   10 * SIZE(B)
	vmovups	%xmm4,   12 * SIZE(B)
	vmovups	%xmm6,   14 * SIZE(B)

	addq	$4 * SIZE, AO1
	addq	$4 * SIZE, AO2
	subq	$-32 * SIZE, B
	ALIGN_4

L(56):
	testq	$2, MM
	jle	L(58)

	vmovups	0 * SIZE(AO1),         %xmm0
	vmovups	1 * SIZE(AO1, LDA),    %xmm1
	vmovups	0 * SIZE(AO1, LDA, 2), %xmm2
	vmovups	1 * SIZE(AO1, LDA3),   %xmm3
	vmovups	0 * SIZE(AO2),         %xmm4
	vmovups	1 * SIZE(AO2, LDA),    %xmm5
	vmovups	0 * SIZE(AO2, LDA, 2), %xmm6
	vmovups	1 * SIZE(AO2, LDA3),   %xmm7

	vmovsd	 %xmm0, %xmm9 , %xmm9
	vmovsd	 %xmm2, %xmm10, %xmm10
	vmovsd	 %xmm4, %xmm11, %xmm11
	vmovsd	 %xmm6, %xmm12, %xmm12

	vmovups	%xmm9,  -16 * SIZE(B)
	vmovups	%xmm10, -14 * SIZE(B)
	vmovups	%xmm11, -12 * SIZE(B)
	vmovups	%xmm12, -10 * SIZE(B)

	vshufpd	 $1, %xmm1, %xmm0 , %xmm0
	vshufpd	 $1, %xmm3, %xmm2 , %xmm2
	vshufpd	 $1, %xmm5, %xmm4 , %xmm4
	vshufpd	 $1, %xmm7, %xmm6 , %xmm6

	vmovups	%xmm0,   -8 * SIZE(B)
	vmovups	%xmm2,   -6 * SIZE(B)
	vmovups	%xmm4,   -4 * SIZE(B)
	vmovups	%xmm6,   -2 * SIZE(B)

	addq	$2 * SIZE, AO1
	addq	$2 * SIZE, AO2
	subq	$-16 * SIZE, B
	ALIGN_4

L(58):
	testq	$1, MM
	jle	L(59)

	vmovsd	0 * SIZE(AO1),         %xmm0
	vmovsd	0 * SIZE(AO1, LDA),    %xmm1
	vmovsd	0 * SIZE(AO1, LDA, 2), %xmm2
	vmovsd	0 * SIZE(AO1, LDA3),   %xmm3
	vmovsd	0 * SIZE(AO2),         %xmm4
	vmovsd	0 * SIZE(AO2, LDA),    %xmm5
	vmovsd	0 * SIZE(AO2, LDA, 2), %xmm6
	vmovsd	0 * SIZE(AO2, LDA3),   %xmm7

	vunpcklpd %xmm1, %xmm0 , %xmm0
	vunpcklpd %xmm3, %xmm2 , %xmm2
	vunpcklpd %xmm5, %xmm4 , %xmm4
	vunpcklpd %xmm7, %xmm6 , %xmm6

	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm2,  -14 * SIZE(B)
	vmovups	%xmm4,  -12 * SIZE(B)
	vmovups	%xmm6,  -10 * SIZE(B)

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

	vmovsd	0 * SIZE(AO1),      %xmm0
	vmovsd	0 * SIZE(AO1, LDA), %xmm1
	vmovsd	0 * SIZE(AO2),      %xmm2
	vmovsd	0 * SIZE(AO2, LDA), %xmm3

	vunpcklpd %xmm1, %xmm0 , %xmm0
	vunpcklpd %xmm3, %xmm2 , %xmm2

	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm2,  -14 * SIZE(B)

	addq	$1 * SIZE, AO1
	addq	$1 * SIZE, AO2
	subq	$-4 * SIZE, B
	ALIGN_3

L(62):
	vmovups	-1 * SIZE(AO1, LDA), %xmm5
	vmovups	-1 * SIZE(AO2, LDA), %xmm7

	movq	MM, I
	sarq	$3, I
	jle	L(64)
	ALIGN_4

L(63):

	vmovups	0 * SIZE(AO1),      %xmm0
	vmovups	1 * SIZE(AO1, LDA), %xmm1
	vmovups	0 * SIZE(AO2),      %xmm2
	vmovups	1 * SIZE(AO2, LDA), %xmm3

	vmovsd	 %xmm0, %xmm5 , %xmm5
	vmovsd	 %xmm2, %xmm7 , %xmm7
	vshufpd	 $1, %xmm1, %xmm0 , %xmm0
	vshufpd	 $1, %xmm3, %xmm2 , %xmm2


	vmovups	%xmm5,  -16 * SIZE(B)
	vmovups	%xmm7,  -14 * SIZE(B)
	vmovups	%xmm0,  -12 * SIZE(B)
	vmovups	%xmm2,  -10 * SIZE(B)


	vmovups	2 * SIZE(AO1),      %xmm0
	vmovups	3 * SIZE(AO1, LDA), %xmm5
	vmovups	2 * SIZE(AO2),      %xmm2
	vmovups	3 * SIZE(AO2, LDA), %xmm7

	vmovsd	 %xmm0, %xmm1 , %xmm1
	vmovsd	 %xmm2, %xmm3 , %xmm3
	vshufpd	 $1, %xmm5, %xmm0 , %xmm0
	vshufpd	 $1, %xmm7, %xmm2 , %xmm2


	vmovups	%xmm1,  -8 * SIZE(B)
	vmovups	%xmm3,  -6 * SIZE(B)
	vmovups	%xmm0,  -4 * SIZE(B)
	vmovups	%xmm2,  -2 * SIZE(B)


	vmovups	4 * SIZE(AO1),      %xmm0
	vmovups	5 * SIZE(AO1, LDA), %xmm1
	vmovups	4 * SIZE(AO2),      %xmm2
	vmovups	5 * SIZE(AO2, LDA), %xmm3

	vmovsd	 %xmm0, %xmm5 , %xmm5
	vmovsd	 %xmm2, %xmm7 , %xmm7
	vshufpd	 $1, %xmm1, %xmm0 , %xmm0
	vshufpd	 $1, %xmm3, %xmm2 , %xmm2


	vmovups	%xmm5,    0 * SIZE(B)
	vmovups	%xmm7,    2 * SIZE(B)
	vmovups	%xmm0,    4 * SIZE(B)
	vmovups	%xmm2,    6 * SIZE(B)


	vmovups	6 * SIZE(AO1),      %xmm0
	vmovups	7 * SIZE(AO1, LDA), %xmm5
	vmovups	6 * SIZE(AO2),      %xmm2
	vmovups	7 * SIZE(AO2, LDA), %xmm7

	vmovsd	 %xmm0, %xmm1 , %xmm1
	vmovsd	 %xmm2, %xmm3 , %xmm3
	vshufpd	 $1, %xmm5, %xmm0 , %xmm0
	vshufpd	 $1, %xmm7, %xmm2 , %xmm2


	vmovups	%xmm1,   8 * SIZE(B)
	vmovups	%xmm3,  10 * SIZE(B)
	vmovups	%xmm0,  12 * SIZE(B)
	vmovups	%xmm2,  14 * SIZE(B)

	addq	$8 * SIZE, AO1
	addq	$8 * SIZE, AO2
	subq	$-32 * SIZE, B

	decq	I
	jg	L(63)
	ALIGN_4

L(64):
	testq	$4, MM
	jle	L(66)

	vmovups	0 * SIZE(AO1),      %xmm0
	vmovups	1 * SIZE(AO1, LDA), %xmm1
	vmovups	0 * SIZE(AO2),      %xmm2
	vmovups	1 * SIZE(AO2, LDA), %xmm3

	vmovsd	 %xmm0, %xmm5 , %xmm5
	vshufpd	 $1, %xmm1, %xmm0 , %xmm0
	vmovsd	 %xmm2, %xmm7 , %xmm7
	vshufpd	 $1, %xmm3, %xmm2 , %xmm2

	vmovups	%xmm5,  -16 * SIZE(B)
	vmovups	%xmm7,  -14 * SIZE(B)
	vmovups	%xmm0,  -12 * SIZE(B)
	vmovups	%xmm2,  -10 * SIZE(B)

	vmovups	2 * SIZE(AO1),      %xmm0
	vmovups	3 * SIZE(AO1, LDA), %xmm5
	vmovups	2 * SIZE(AO2),      %xmm2
	vmovups	3 * SIZE(AO2, LDA), %xmm7

	vmovsd	 %xmm0, %xmm1 , %xmm1
	vshufpd	 $1, %xmm5, %xmm0 , %xmm0
	vmovsd	 %xmm2, %xmm3 , %xmm3
	vshufpd	 $1, %xmm7, %xmm2 , %xmm2

	vmovups	%xmm1,  -8 * SIZE(B)
	vmovups	%xmm3,  -6 * SIZE(B)
	vmovups	%xmm0,  -4 * SIZE(B)
	vmovups	%xmm2,  -2 * SIZE(B)

	addq	$4 * SIZE, AO1
	addq	$4 * SIZE, AO2
	subq	$-16 * SIZE, B
	ALIGN_4

L(66):
	testq	$2, MM
	jle	L(68)

	vmovups	0 * SIZE(AO1),      %xmm0
	vmovups	1 * SIZE(AO1, LDA), %xmm1
	vmovups	0 * SIZE(AO2),      %xmm2
	vmovups	1 * SIZE(AO2, LDA), %xmm3

	vmovsd	 %xmm0, %xmm5 , %xmm5
	vmovsd	 %xmm2, %xmm7 , %xmm7
	vshufpd	 $1, %xmm1, %xmm0 , %xmm0
	vshufpd	 $1, %xmm3, %xmm2 , %xmm2

	vmovups	%xmm5,  -16 * SIZE(B)
	vmovups	%xmm7,  -14 * SIZE(B)
	vmovups	%xmm0,  -12 * SIZE(B)
	vmovups	%xmm2,  -10 * SIZE(B)

	addq	$2 * SIZE, AO1
	addq	$2 * SIZE, AO2
	subq	$-8 * SIZE, B
	ALIGN_4

L(68):
	testq	$1, MM
	jle	L(70)

	vmovsd	0 * SIZE(AO1),      %xmm0
	vmovsd	0 * SIZE(AO1, LDA), %xmm1
	vmovsd	0 * SIZE(AO2),      %xmm2
	vmovsd	0 * SIZE(AO2, LDA), %xmm3

	vunpcklpd %xmm1, %xmm0 , %xmm0
	vunpcklpd %xmm3, %xmm2 , %xmm2

	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm2,  -14 * SIZE(B)
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

	vmovsd	0 * SIZE(AO1), %xmm0
	vmovsd	0 * SIZE(AO2), %xmm1

	vunpcklpd %xmm1, %xmm0 , %xmm0

	vmovups	%xmm0,  -16 * SIZE(B)

	addq	$1 * SIZE, AO1
	addq	$1 * SIZE, AO2
	subq	$-2 * SIZE, B
	ALIGN_3

L(72):
	vmovups	-1 * SIZE(AO2), %xmm5

	movq	MM, I
	sarq	$3, I
	jle	L(74)
	ALIGN_4

L(73):

	vmovups	0 * SIZE(AO1), %xmm0
	vmovups	1 * SIZE(AO2), %xmm1
	vmovups	2 * SIZE(AO1), %xmm2
	vmovups	3 * SIZE(AO2), %xmm3

	vmovsd	 %xmm0, %xmm5 , %xmm5
	vshufpd	 $1, %xmm1, %xmm0 , %xmm0
	vmovsd	 %xmm2, %xmm1 , %xmm1
	vshufpd	 $1, %xmm3, %xmm2 , %xmm2


	vmovups	%xmm5,  -16 * SIZE(B)
	vmovups	%xmm0,  -14 * SIZE(B)
	vmovups	%xmm1,  -12 * SIZE(B)
	vmovups	%xmm2,  -10 * SIZE(B)


	vmovups	4 * SIZE(AO1), %xmm0
	vmovups	5 * SIZE(AO2), %xmm1
	vmovups	6 * SIZE(AO1), %xmm2
	vmovups	7 * SIZE(AO2), %xmm5

	vmovsd	 %xmm0, %xmm3 , %xmm3
	vshufpd	 $1, %xmm1, %xmm0 , %xmm0
	vmovsd	 %xmm2, %xmm1 , %xmm1
	vshufpd	 $1, %xmm5, %xmm2 , %xmm2


	vmovups	%xmm3,   -8 * SIZE(B)
	vmovups	%xmm0,   -6 * SIZE(B)
	vmovups	%xmm1,   -4 * SIZE(B)
	vmovups	%xmm2,   -2 * SIZE(B)

	addq	$8 * SIZE, AO1
	addq	$8 * SIZE, AO2
	subq	$-16 * SIZE, B

	decq	I
	jg	L(73)
	ALIGN_4

L(74):
	testq	$4, MM
	jle	L(76)

	vmovups	0 * SIZE(AO1), %xmm0
	vmovups	1 * SIZE(AO2), %xmm1
	vmovups	2 * SIZE(AO1), %xmm2
	vmovups	3 * SIZE(AO2), %xmm3

	vmovsd	 %xmm0, %xmm5 , %xmm5
	vshufpd	 $1, %xmm1, %xmm0 , %xmm0
	vmovsd	 %xmm2, %xmm1 , %xmm1
	vshufpd	 $1, %xmm3, %xmm2 , %xmm2

	vmovups	%xmm5,  -16 * SIZE(B)
	vmovups	%xmm0,  -14 * SIZE(B)
	vmovups	%xmm1,  -12 * SIZE(B)
	vmovups	%xmm2,  -10 * SIZE(B)

	vmovups	%xmm3, %xmm5

	addq	$4 * SIZE, AO1
	addq	$4 * SIZE, AO2
	subq	$-8 * SIZE, B
	ALIGN_4

L(76):
	testq	$2, MM
	jle	L(78)

	vmovups	0 * SIZE(AO1), %xmm0
	vmovups	1 * SIZE(AO2), %xmm1

	vmovsd	 %xmm0, %xmm5 , %xmm5
	vshufpd	 $1, %xmm1, %xmm0 , %xmm0

	vmovups	%xmm5,  -16 * SIZE(B)
	vmovups	%xmm0,  -14 * SIZE(B)

	addq	$2 * SIZE, AO1
	addq	$2 * SIZE, AO2
	subq	$-4 * SIZE, B
	ALIGN_4

L(78):
	testq	$1, MM
	jle	L(80)

	vmovsd	0 * SIZE(AO1),      %xmm0
	vmovsd	0 * SIZE(AO2),      %xmm1

	vunpcklpd %xmm1, %xmm0 , %xmm0

	vmovups	%xmm0,  -16 * SIZE(B)
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

	vmovups	0 * SIZE(AO1),      %xmm0
	vmovups	2 * SIZE(AO1),      %xmm2
	vmovups	4 * SIZE(AO1),      %xmm4
	vmovups	6 * SIZE(AO1),      %xmm6


	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm2,  -14 * SIZE(B)
	vmovups	%xmm4,  -12 * SIZE(B)
	vmovups	%xmm6,  -10 * SIZE(B)

	addq	$8 * SIZE, AO1
	subq	$-8 * SIZE, B

	decq	I
	jg	L(81)
	ALIGN_4

L(82):
	testq	$4, MM
	jle	L(83)

	vmovups	0 * SIZE(AO1),      %xmm0
	vmovups	2 * SIZE(AO1),      %xmm2

	vmovups	%xmm0, -16 * SIZE(B)
	vmovups	%xmm2, -14 * SIZE(B)

	addq	$4 * SIZE, AO1
	subq	$-4 * SIZE, B
	ALIGN_4

L(83):
	testq	$2, MM
	jle	L(84)

	vmovups	0 * SIZE(AO1),      %xmm0

	vmovups	%xmm0,  -16 * SIZE(B)

	addq	$2 * SIZE, AO1
	subq	$-2 * SIZE, B
	ALIGN_4

L(84):
	testq	$1, MM
	jle	L(999)

	vmovsd	0 * SIZE(AO1),      %xmm0

	vmovlpd	%xmm0,  -16 * SIZE(B)
	jmp	L(999)
	ALIGN_4

L(85):
	vmovups	-1 * SIZE(AO1),      %xmm0

	movq	M, I
	sarq	$3, I
	jle	L(86)
	ALIGN_4

L(86):

	vmovups	1 * SIZE(AO1),      %xmm1
	vmovups	3 * SIZE(AO1),      %xmm2
	vmovups	5 * SIZE(AO1),      %xmm3
	vmovups	7 * SIZE(AO1),      %xmm4

	vshufpd	$1, %xmm1, %xmm0 , %xmm0
	vshufpd	$1, %xmm2, %xmm1 , %xmm1
	vshufpd	$1, %xmm3, %xmm2 , %xmm2
	vshufpd	$1, %xmm4, %xmm3 , %xmm3


	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm1,  -14 * SIZE(B)
	vmovups	%xmm2,  -12 * SIZE(B)
	vmovups	%xmm3,  -10 * SIZE(B)

	vmovups	%xmm4, %xmm0

	addq	$8 * SIZE, AO1
	subq	$-8 * SIZE, B

	decq	I
	jg	L(86)
	ALIGN_4

L(87):
	testq	$4, M
	jle	L(88)

	vmovups	1 * SIZE(AO1),      %xmm1
	vmovups	3 * SIZE(AO1),      %xmm2

	vshufpd	$1, %xmm1, %xmm0 , %xmm0
	vshufpd	$1, %xmm2, %xmm1 , %xmm1

	vmovups	%xmm0,  -16 * SIZE(B)
	vmovups	%xmm1,  -14 * SIZE(B)

	vmovups	%xmm2, %xmm0

	addq	$4 * SIZE, AO1
	addq	$4 * SIZE, B
	ALIGN_4

L(88):
	testq	$2, M
	jle	L(89)

	vmovups	1 * SIZE(AO1),      %xmm1

	vshufpd	$1, %xmm1, %xmm0 , %xmm0

	vmovups	%xmm0,  -16 * SIZE(B)

	vmovups	%xmm1, %xmm0

	addq	$2 * SIZE, AO1
	subq	$-2 * SIZE, B
	ALIGN_4

L(89):
	testq	$1, M
	jle	L(999)

	vshufpd	$1, %xmm0, %xmm0 , %xmm0

	vmovlpd	%xmm0,  -16 * SIZE(B)
	ALIGN_4

L(999):
#ifdef WINDOWS_ABI
	vmovups	  0(%rsp), %xmm6
	vmovups	 16(%rsp), %xmm7
	vmovups	 32(%rsp), %xmm8
	vmovups	 48(%rsp), %xmm9
	vmovups	 64(%rsp), %xmm10
	vmovups	 80(%rsp), %xmm11
	vmovups	 96(%rsp), %xmm12

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
