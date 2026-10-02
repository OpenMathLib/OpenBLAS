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

#define M	ARG1
#define N	ARG2
#define K	ARG3
#define A	ARG4
#define B	ARG5
#define C	ARG6
#define LDC	%r10

#define I	%r12
#define J	%r13
#define AO	%r14
#define BO	%r15
#define	CO	%rbp

#define KK	%r11
#define AORIG	 48(%rsp)

#define STACKSIZE 64

#define ALPHA	 8 + STACKSIZE(%rsp)
#define OFFSET	32 + STACKSIZE(%rsp)

#ifdef OPTERON
#define PREFETCH	prefetch
#define PREFETCHW	prefetchw
#else
#define PREFETCH	prefetcht0
#define PREFETCHW	prefetcht0
#endif

#define PREFETCHSIZE (5 + 4 * 10)

	PROLOGUE
	PROFCODE

#ifdef WINDOWS_ABI
	emms
#endif

	subq	$STACKSIZE, %rsp
	movq	%rbx,  0(%rsp)
	movq	%rbp,  8(%rsp)
	movq	%r12, 16(%rsp)
	movq	%r13, 24(%rsp)
	movq	%r14, 32(%rsp)
	movq	%r15, 40(%rsp)

	movq	24 + STACKSIZE(%rsp), LDC

#if defined(TRMMKERNEL) && !defined(LEFT)
	movq	OFFSET, %rax
	negq	%rax
	movq	%rax, KK
#endif

	addq	$8 * SIZE, A
	addq	$8 * SIZE, B

	salq	$BASE_SHIFT, LDC

#ifdef LN
       movq	M, %rax
       salq	$BASE_SHIFT, %rax
       addq	%rax, C
       imulq	K, %rax
       addq	%rax, A
#endif

#ifdef RT
       movq	N, %rax
       salq	$BASE_SHIFT, %rax
       imulq	K, %rax
       addq	%rax, B

       movq	N,   %rax
       imulq	LDC, %rax
       addq	%rax, C
#endif

#ifdef RN
       movq	OFFSET, %rax
       negq	%rax
       movq	%rax, KK
#endif

#ifdef RT
       movq	N, %rax
       subq	OFFSET, %rax
       movq	%rax, KK
#endif

	movq	N,  %rax
	testq	$1, %rax
	je	L(30)

#if defined(LT) || defined(RN)
	movq	A, AO
#else
	movq	A, %rax
	movq	%rax, AORIG
#endif

#ifdef RT
	movq	K, %rax
	salq	$0 + BASE_SHIFT, %rax
	subq	%rax, B
#endif

#ifdef RT
	subq	LDC, C
#endif
	movq	C, CO
#ifndef RT
	addq	LDC, C
#endif

#ifdef LN
	movq	OFFSET, %rax
	addq	M, %rax
	movq	%rax, KK
#endif

#ifdef LT
	movq	OFFSET, %rax
	movq	%rax, KK
#endif

	movq	M,  I
	sarq	$1, I
	je	L(40)
	ALIGN_4

L(31):
#ifdef LN
       movq	K, %rax
       salq	$1 + BASE_SHIFT, %rax
       subq	%rax, AORIG
#endif

#if defined(LN) || defined(RT)
	movq	KK, %rax
	salq	$BASE_SHIFT, %rax
	movq	AORIG, AO
	leaq	(AO, %rax, 2), AO
	leaq	(B,  %rax, 1), BO
#else
	movq	B, BO
#endif

	fldz
	fldz

#if   defined(HAVE_3DNOW)
	prefetchw	2 * SIZE(CO)
#elif defined(HAVE_SSE)
	prefetchnta	2 * SIZE(CO)
#endif

#if defined(LT) || defined(RN)
	movq	KK, %rax
#else
	movq	K,  %rax
	subq	KK, %rax
#endif
	sarq	$2, %rax
 	je	L(35)
	ALIGN_4

L(32):
	PREFETCH	(PREFETCHSIZE + 0) * SIZE(AO)

	FLD	 -8 * SIZE(BO)
	FLD	 -8 * SIZE(AO)
	fmul	 %st(1), %st
	faddp	 %st, %st(2)

	FLD	 -7 * SIZE(AO)
	fmulp	 %st, %st(1)
	faddp	 %st, %st(2)

	FLD	 -7 * SIZE(BO)
	FLD	 -6 * SIZE(AO)
	fmul	 %st(1), %st
	faddp	 %st, %st(2)

	FLD	 -5 * SIZE(AO)
	fmulp	 %st, %st(1)
	faddp	 %st, %st(2)

	FLD	 -6 * SIZE(BO)
	FLD	 -4 * SIZE(AO)
	fmul	 %st(1), %st
	faddp	 %st, %st(2)

	FLD	 -3 * SIZE(AO)
	fmulp	 %st, %st(1)
	faddp	 %st, %st(2)

	FLD	 -5 * SIZE(BO)
	FLD	 -2 * SIZE(AO)
	fmul	 %st(1), %st
	faddp	 %st, %st(2)

	FLD	 -1 * SIZE(AO)
	fmulp	 %st, %st(1)
	faddp	 %st, %st(2)

	addq	$8 * SIZE,AO
	addq	$4 * SIZE,BO

	decq	%rax
	jne	L(32)
	ALIGN_4

L(35):
#if defined(LT) || defined(RN)
	movq	KK, %rax
#else
	movq	K,  %rax
	subq	KK, %rax
#endif
	and	$3,  %rax
	je	L(38)
	ALIGN_4

L(36):
	FLD	 -8 * SIZE(BO)

	FLD	 -8 * SIZE(AO)
	fmul	 %st(1), %st
	faddp	 %st, %st(2)

	FLD	 -7 * SIZE(AO)
	fmulp	 %st, %st(1)
	faddp	 %st, %st(2)

	addq	$2 * SIZE,AO
	addq	$1 * SIZE,BO

	decq	%rax
	jne	 L(36)
	ALIGN_4

L(38):
#if defined(LN) || defined(RT)
	movq	KK, %rax
#ifdef LN
	subq	$2, %rax
#else
	subq	$1, %rax
#endif

	salq	$BASE_SHIFT, %rax

	movq	AORIG, AO
	leaq	(AO, %rax, 2), AO
	leaq	(B,  %rax, 1), BO
#endif

#if defined(LN) || defined(LT)
	FLD	-8 * SIZE(BO)
	fsubp	%st, %st(1)
	FLD	-7 * SIZE(BO)
	fsubp	%st, %st(2)
#else
	FLD	-8 * SIZE(AO)
	fsubp	%st, %st(1)
	FLD	-7 * SIZE(AO)
	fsubp	%st, %st(2)
#endif

#ifdef LN
       FLD	-5 * SIZE(AO)
       fmulp	%st, %st(2)

       FLD	-6 * SIZE(AO)
       fmul	%st(2), %st

       fsubrp	%st, %st(1)
       FLD	-8 * SIZE(AO)
       fmulp	%st, %st(1)
#endif

#ifdef LT
       FLD	-8 * SIZE(AO)
       fmulp	%st, %st(1)

       FLD	-7 * SIZE(AO)
       fmul	%st(1), %st

       fsubrp	%st, %st(2)

       FLD	-5 * SIZE(AO)
       fmulp	%st, %st(2)
#endif

#ifdef RN
       FLD	-8 * SIZE(BO)
       fmul	%st, %st(1)
       fmulp	%st, %st(2)
#endif

#ifdef RT
       FLD	-8 * SIZE(BO)
       fmul	%st, %st(1)
       fmulp	%st, %st(2)
#endif

#ifdef LN
	subq	$2 * SIZE, CO
#endif

#if defined(LN) || defined(LT)
	fld	%st
	FST	-8 * SIZE(BO)
	fxch	%st(1)
	fld	%st
	FST	-7 * SIZE(BO)
#else
	fld	%st
	FST	-8 * SIZE(AO)
	fxch	%st(1)
	fld	%st
	FST	-7 * SIZE(AO)
#endif

	FST	1 * SIZE(CO)
	FST	0 * SIZE(CO)

#ifndef LN
	addq	$2 * SIZE, CO
#endif

#if defined(LT) || defined(RN)
	movq	K,  %rax
	subq	KK, %rax
	salq	$BASE_SHIFT, %rax
	leaq	(AO, %rax, 2), AO
	leaq	(BO, %rax, 1), BO
#endif

#ifdef LN
	subq	$2, KK
#endif

#ifdef LT
	addq	$2, KK
#endif

#ifdef RT
       movq	K, %rax
       salq	$1 + BASE_SHIFT, %rax
       addq	%rax, AORIG
#endif

	decq	I
	jne	L(31)
	ALIGN_4

L(40):
	movq	 M, %rax
	andq	$1, %rax
	je	L(49)
	ALIGN_4

L(41):
#ifdef LN
       movq	K, %rax
       salq	$0 + BASE_SHIFT, %rax
       subq	%rax, AORIG
#endif

#if defined(LN) || defined(RT)
	movq	KK, %rax
	salq	$BASE_SHIFT, %rax
	movq	AORIG, AO
	leaq	(AO, %rax, 1), AO
	leaq	(B,  %rax, 1), BO
#else
	movq	B, BO
#endif

	fldz

#if defined(LT) || defined(RN)
	movq	KK, %rax
#else
	movq	K,  %rax
	subq	KK, %rax
#endif
	sarq	$2, %rax
 	je	L(45)
	ALIGN_4

L(42):
	PREFETCH	(PREFETCHSIZE + 0) * SIZE(AO)

	FLD	 -8 * SIZE(AO)
	FLD	 -8 * SIZE(BO)
	fmulp	 %st, %st(1)
	faddp	 %st, %st(1)

	FLD	 -7 * SIZE(AO)
	FLD	 -7 * SIZE(BO)
	fmulp	 %st, %st(1)
	faddp	 %st, %st(1)

	FLD	 -6 * SIZE(AO)
	FLD	 -6 * SIZE(BO)
	fmulp	 %st, %st(1)
	faddp	 %st, %st(1)

	FLD	 -5 * SIZE(AO)
	FLD	 -5 * SIZE(BO)
	fmulp	 %st, %st(1)
	faddp	 %st, %st(1)

	addq	$4 * SIZE,AO
	addq	$4 * SIZE,BO

	decq	%rax
	jne	L(42)
	ALIGN_4

L(45):
#if defined(LT) || defined(RN)
	movq	KK, %rax
#else
	movq	K,  %rax
	subq	KK, %rax
#endif
	and	$3,  %rax
	je	L(48)
	ALIGN_4

L(46):
	FLD	 -8 * SIZE(AO)

	FLD	 -8 * SIZE(BO)
	fmulp	 %st, %st(1)
	faddp	 %st, %st(1)

	addq	$1 * SIZE,AO
	addq	$1 * SIZE,BO

	decq	%rax
	jne	 L(46)
	ALIGN_4

L(48):
#if defined(LN) || defined(RT)
	movq	KK, %rax
#ifdef LN
	subq	$1, %rax
#else
	subq	$1, %rax
#endif

	salq	$BASE_SHIFT, %rax

	movq	AORIG, AO
	leaq	(AO, %rax, 1), AO
	leaq	(B,  %rax, 1), BO
#endif

#if defined(LN) || defined(LT)
	FLD	-8 * SIZE(BO)
	fsubp	%st, %st(1)
#else
	FLD	-8 * SIZE(AO)
	fsubp	%st, %st(1)
#endif

#ifdef LN
       FLD	-8 * SIZE(AO)
       fmulp	%st, %st(1)
#endif

#ifdef LT
       FLD	-8 * SIZE(AO)
       fmulp	%st, %st(1)
#endif

#ifdef RN
       FLD	-8 * SIZE(BO)
       fmulp	%st, %st(1)
#endif

#ifdef RT
       FLD	-8 * SIZE(BO)
       fmulp	%st, %st(1)
#endif

#ifdef LN
	subq	$1 * SIZE, CO
#endif

#if defined(LN) || defined(LT)
	fld	%st
	FST	-8 * SIZE(BO)
#else
	fld	%st
	FST	-8 * SIZE(AO)
#endif

	FST	0 * SIZE(CO)

#ifndef LN
	addq	$1 * SIZE, CO
#endif

#if defined(LT) || defined(RN)
	movq	K,  %rax
	subq	KK, %rax
	salq	$BASE_SHIFT, %rax
	leaq	(AO, %rax, 1), AO
	leaq	(BO, %rax, 1), BO
#endif

#ifdef LN
	subq	$1, KK
#endif

#ifdef LT
	addq	$1, KK
#endif

#ifdef RT
       movq	K, %rax
       salq	$0 + BASE_SHIFT, %rax
       addq	%rax, AORIG
#endif
	ALIGN_4

L(49):
#ifdef LN
       movq	K, %rax
       salq	$BASE_SHIFT, %rax
       leaq	(B, %rax, 1), B
#endif

#if defined(LT) || defined(RN)
	movq	BO, B
#endif

#ifdef RN
	addq	$1, KK
#endif

#ifdef RT
	subq	$1, KK
#endif
	ALIGN_4

L(30):
	movq	N,   %rax
	sarq	$1,  %rax
	movq	%rax, J
	je	L(999)
	ALIGN_4

L(01):
#if defined(LT) || defined(RN)
	movq	A, AO
#else
	movq	A, %rax
	movq	%rax, AORIG
#endif

#ifdef RT
	movq	K, %rax
	salq	$1 + BASE_SHIFT, %rax
	subq	%rax, B
#endif

	lea	(, LDC, 2), %rax

#ifdef RT
	subq	%rax, C
#endif
	movq	C, CO
#ifndef RT
	addq	%rax, C
#endif

#ifdef LN
	movq	OFFSET, %rax
	addq	M, %rax
	movq	%rax, KK
#endif

#ifdef LT
	movq	OFFSET, %rax
	movq	%rax, KK
#endif

	movq	M,  I
	sarq	$1, I
	je	L(20)
	ALIGN_4

L(11):
#ifdef LN
       movq	K, %rax
       salq	$1 + BASE_SHIFT, %rax
       subq	%rax, AORIG
#endif

#if defined(LN) || defined(RT)
	movq	KK, %rax
	salq	$BASE_SHIFT, %rax
	movq	AORIG, AO
	leaq	(AO, %rax, 2), AO
	leaq	(B,  %rax, 2), BO
#else
	movq	B, BO
#endif

	fldz
	fldz
	fldz
	fldz

#if   defined(HAVE_3DNOW)
	prefetchw	2 * SIZE(CO)
 	prefetchw	2 * SIZE(CO, LDC, 1)
#elif defined(HAVE_SSE)
	prefetchnta	2 * SIZE(CO)
 	prefetchnta	2 * SIZE(CO, LDC, 1)
#endif

#if defined(LT) || defined(RN)
	movq	KK, %rax
#else
	movq	K,  %rax
	subq	KK, %rax
#endif
	sarq	$2, %rax
 	je	L(15)
	ALIGN_4

L(12):
	PREFETCH	(PREFETCHSIZE + 0) * SIZE(AO)

	FLD	 -8 * SIZE(AO)

	FLD	 -8 * SIZE(BO)
	fld	 %st(1)
	fmul	 %st(1), %st
	faddp	 %st, %st(3)

	FLD	 -7 * SIZE(BO)
	fmul	 %st, %st(2)

	FLD	 -7 * SIZE(AO)
	fmul	 %st, %st(2)
	fmulp	 %st, %st(1)

	faddp	 %st, %st(6)
	faddp	 %st, %st(4)
	faddp	 %st, %st(2)

	FLD	 -6 * SIZE(AO)

	FLD	 -6 * SIZE(BO)
	fld	 %st(1)
	fmul	 %st(1), %st
	faddp	 %st, %st(3)

	FLD	 -5 * SIZE(BO)
	fmul	 %st, %st(2)

	FLD	 -5 * SIZE(AO)
	fmul	 %st, %st(2)
	fmulp	 %st, %st(1)

	faddp	 %st, %st(6)
	faddp	 %st, %st(4)
	faddp	 %st, %st(2)

	PREFETCH	(PREFETCHSIZE + 4) * SIZE(AO)

	FLD	 -4 * SIZE(AO)

	FLD	 -4 * SIZE(BO)
	fld	 %st(1)
	fmul	 %st(1), %st
	faddp	 %st, %st(3)

	FLD	 -3 * SIZE(BO)
	fmul	 %st, %st(2)

	FLD	 -3 * SIZE(AO)
	fmul	 %st, %st(2)
	fmulp	 %st, %st(1)

	faddp	 %st, %st(6)
	faddp	 %st, %st(4)
	faddp	 %st, %st(2)

	FLD	 -2 * SIZE(AO)

	FLD	 -2 * SIZE(BO)
	fld	 %st(1)
	fmul	 %st(1), %st
	faddp	 %st, %st(3)

	FLD	 -1 * SIZE(BO)
	fmul	 %st, %st(2)

	FLD	 -1 * SIZE(AO)
	fmul	 %st, %st(2)
	fmulp	 %st, %st(1)

	faddp	 %st, %st(6)
	faddp	 %st, %st(4)
	faddp	 %st, %st(2)

	addq	$8 * SIZE,AO
	addq	$8 * SIZE,BO

	decq	%rax
	jne	L(12)
	ALIGN_4

L(15):
#if defined(LT) || defined(RN)
	movq	KK, %rax
#else
	movq	K,  %rax
	subq	KK, %rax
#endif
	and	$3,  %rax
	je	L(18)
	ALIGN_4

L(16):
	FLD	 -8 * SIZE(AO)

	FLD	 -8 * SIZE(BO)
	fld	 %st(1)
	fmul	 %st(1), %st
	faddp	 %st, %st(3)

	FLD	 -7 * SIZE(BO)
	fmul	 %st, %st(2)

	FLD	 -7 * SIZE(AO)
	fmul	 %st, %st(2)
	fmulp	 %st, %st(1)

	faddp	 %st, %st(6)
	faddp	 %st, %st(4)
	faddp	 %st, %st(2)

	addq	$2 * SIZE,AO
	addq	$2 * SIZE,BO

	decq	%rax
	jne	 L(16)
	ALIGN_4

L(18):
#if defined(LN) || defined(RT)
	movq	KK, %rax
#ifdef LN
	subq	$2, %rax
#else
	subq	$2, %rax
#endif

	salq	$BASE_SHIFT, %rax

	movq	AORIG, AO
	leaq	(AO, %rax, 2), AO
	leaq	(B,  %rax, 2), BO
#endif

#if defined(LN) || defined(LT)
	FLD	-8 * SIZE(BO)
	fsubp	%st, %st(1)
	FLD	-7 * SIZE(BO)
	fsubp	%st, %st(2)
	FLD	-6 * SIZE(BO)
	fsubp	%st, %st(3)
	FLD	-5 * SIZE(BO)
	fsubp	%st, %st(4)
#else
	FLD	-8 * SIZE(AO)
	fsubp	%st, %st(1)
	FLD	-7 * SIZE(AO)
	fsubp	%st, %st(3)
	FLD	-6 * SIZE(AO)
	fsubp	%st, %st(2)
	FLD	-5 * SIZE(AO)
	fsubp	%st, %st(4)
#endif

#ifdef LN
       FLD	-5 * SIZE(AO)
       fmul	%st, %st(3)
       fmulp	%st, %st(4)

       FLD	-6 * SIZE(AO)
       fmul	%st(3), %st
       FLD	-6 * SIZE(AO)
       fmul	%st(5), %st

       fsubrp	%st, %st(3)
       fsubrp	%st, %st(1)

       FLD	-8 * SIZE(AO)
       fmul	%st, %st(1)
       fmulp	%st, %st(2)
#endif

#ifdef LT
       FLD	-8 * SIZE(AO)
       fmul	%st, %st(1)
       fmulp	%st, %st(2)

       FLD	-7 * SIZE(AO)
       fmul	%st(1), %st
       FLD	-7 * SIZE(AO)
       fmul	%st(3), %st

       fsubrp	%st, %st(5)
       fsubrp	%st, %st(3)

       FLD	-5 * SIZE(AO)
       fmul	%st, %st(3)
       fmulp	%st, %st(4)
#endif

#ifdef RN
       FLD	-8 * SIZE(BO)
       fmul	%st, %st(1)
       fmulp	%st, %st(3)

       FLD	-7 * SIZE(BO)
       fmul	%st(1), %st
       FLD	-7 * SIZE(BO)
       fmul	%st(4), %st

       fsubrp	%st, %st(5)
       fsubrp	%st, %st(2)

       FLD	-5 * SIZE(BO)
       fmul	%st, %st(2)
       fmulp	%st, %st(4)
#endif

#ifdef RT
       FLD	-5 * SIZE(BO)
       fmul	%st, %st(2)
       fmulp	%st, %st(4)

       FLD	-6 * SIZE(BO)
       fmul	%st(2), %st
       FLD	-6 * SIZE(BO)
       fmul	%st(5), %st

       fsubrp	%st, %st(4)
       fsubrp	%st, %st(1)

       FLD	-8 * SIZE(BO)
       fmul	%st, %st(1)
       fmulp	%st, %st(3)
#endif

#ifdef LN
	subq	$2 * SIZE, CO
#endif

#if defined(LN) || defined(LT)
	fld	%st
	FST	-8 * SIZE(BO)
	fxch	%st(1)
	fld	%st
	FST	-7 * SIZE(BO)
	fxch	%st(2)
	fld	%st
	FST	-6 * SIZE(BO)
	fxch	%st(3)
	fld	%st
	FST	-5 * SIZE(BO)

	FST	1 * SIZE(CO, LDC)
	FST	0 * SIZE(CO)
	FST	0 * SIZE(CO, LDC)
	FST	1 * SIZE(CO)
#else
	fld	%st
	FST	-8 * SIZE(AO)
	fxch	%st(2)
	fld	%st
	FST	-7 * SIZE(AO)
	fxch	%st(1)
	fld	%st
	FST	-6 * SIZE(AO)
	fxch	%st(3)
	fld	%st
	FST	-5 * SIZE(AO)

	FST	1 * SIZE(CO, LDC)
	FST	1 * SIZE(CO)
	FST	0 * SIZE(CO)
	FST	0 * SIZE(CO, LDC)
#endif

#ifndef LN
	addq	$2 * SIZE, CO
#endif

#if defined(LT) || defined(RN)
	movq	K,  %rax
	subq	KK, %rax
	salq	$BASE_SHIFT, %rax
	leaq	(AO, %rax, 2), AO
	leaq	(BO, %rax, 2), BO
#endif

#ifdef LN
	subq	$2, KK
#endif

#ifdef LT
	addq	$2, KK
#endif

#ifdef RT
       movq	K, %rax
       salq	$1 + BASE_SHIFT, %rax
       addq	%rax, AORIG
#endif

	decq	I
	jne	L(11)
	ALIGN_4

L(20):
	movq	 M, %rax
	andq	$1, %rax
	je	L(29)
	ALIGN_4

L(21):
#ifdef LN
       movq	K, %rax
       salq	$0 + BASE_SHIFT, %rax
       subq	%rax, AORIG
#endif

#if defined(LN) || defined(RT)
	movq	KK, %rax
	salq	$BASE_SHIFT, %rax
	movq	AORIG, AO
	leaq	(AO, %rax, 1), AO
	leaq	(B,  %rax, 2), BO
#else
	movq	B, BO
#endif

	fldz
	fldz

#if defined(LT) || defined(RN)
	movq	KK, %rax
#else
	movq	K,  %rax
	subq	KK, %rax
#endif
	sarq	$2, %rax
 	je	L(25)
	ALIGN_4

L(22):
	PREFETCH	(PREFETCHSIZE + 0) * SIZE(AO)

	FLD	 -8 * SIZE(AO)

	FLD	 -8 * SIZE(BO)
	fmul	 %st(1), %st
	faddp	 %st, %st(2)

	FLD	 -7 * SIZE(BO)
	fmulp	 %st, %st(1)
	faddp	 %st, %st(2)

	FLD	 -7 * SIZE(AO)

	FLD	 -6 * SIZE(BO)
	fmul	 %st(1), %st
	faddp	 %st, %st(2)

	FLD	 -5 * SIZE(BO)
	fmulp	 %st, %st(1)
	faddp	 %st, %st(2)

	FLD	 -6 * SIZE(AO)

	FLD	 -4 * SIZE(BO)
	fmul	 %st(1), %st
	faddp	 %st, %st(2)

	FLD	 -3 * SIZE(BO)
	fmulp	 %st, %st(1)
	faddp	 %st, %st(2)

	FLD	 -5 * SIZE(AO)

	FLD	 -2 * SIZE(BO)
	fmul	 %st(1), %st
	faddp	 %st, %st(2)

	FLD	 -1 * SIZE(BO)
	fmulp	 %st, %st(1)
	faddp	 %st, %st(2)

	addq	$4 * SIZE,AO
	addq	$8 * SIZE,BO

	decq	%rax
	jne	L(22)
	ALIGN_4

L(25):
#if defined(LT) || defined(RN)
	movq	KK, %rax
#else
	movq	K,  %rax
	subq	KK, %rax
#endif
	and	$3,  %rax
	je	L(28)
	ALIGN_4

L(26):
	FLD	 -8 * SIZE(AO)

	FLD	 -8 * SIZE(BO)
	fmul	 %st(1), %st
	faddp	 %st, %st(2)

	FLD	 -7 * SIZE(BO)
	fmulp	 %st, %st(1)
	faddp	 %st, %st(2)

	addq	$1 * SIZE,AO
	addq	$2 * SIZE,BO

	decq	%rax
	jne	 L(26)
	ALIGN_4

L(28):
#if defined(LN) || defined(RT)
	movq	KK, %rax
#ifdef LN
	subq	$1, %rax
#else
	subq	$2, %rax
#endif

	salq	$BASE_SHIFT, %rax

	movq	AORIG, AO
	leaq	(AO, %rax, 1), AO
	leaq	(B,  %rax, 2), BO
#endif

#if defined(LN) || defined(LT)
	FLD	-8 * SIZE(BO)
	fsubp	%st, %st(1)
	FLD	-7 * SIZE(BO)
	fsubp	%st, %st(2)
#else
	FLD	-8 * SIZE(AO)
	fsubp	%st, %st(1)
	FLD	-7 * SIZE(AO)
	fsubp	%st, %st(2)
#endif

#if defined(LN) || defined(LT)
       FLD	-8 * SIZE(AO)
       fmul	%st, %st(1)
       fmulp	%st, %st(2)
#endif

#ifdef RN
       FLD	-8 * SIZE(BO)
       fmulp	%st, %st(1)

       FLD	-7 * SIZE(BO)
       fmul	%st(1), %st

       fsubrp	%st, %st(2)

       FLD	-5 * SIZE(BO)
       fmulp	%st, %st(2)
#endif

#ifdef RT
       FLD	-5 * SIZE(BO)
       fmulp	%st, %st(2)

       FLD	-6 * SIZE(BO)
       fmul	%st(2), %st

       fsubrp	%st, %st(1)

       FLD	-8 * SIZE(BO)
       fmulp	%st, %st(1)
#endif

#ifdef LN
	subq	$1 * SIZE, CO
#endif

#if defined(LN) || defined(LT)
	fld	%st
	FST	-8 * SIZE(BO)
	fxch	%st(1)
	fld	%st
	FST	-7 * SIZE(BO)
#else
	fld	%st
	FST	-8 * SIZE(AO)
	fxch	%st(1)
	fld	%st
	FST	-7 * SIZE(AO)
#endif

	FST	0 * SIZE(CO, LDC)
	FST	0 * SIZE(CO)

#ifndef LN
	addq	$1 * SIZE, CO
#endif

#if defined(LT) || defined(RN)
	movq	K,  %rax
	subq	KK, %rax
	salq	$BASE_SHIFT, %rax
	leaq	(AO, %rax, 1), AO
	leaq	(BO, %rax, 2), BO
#endif

#ifdef LN
	subq	$1, KK
#endif

#ifdef LT
	addq	$1, KK
#endif

#ifdef RT
       movq	K, %rax
       salq	$0 + BASE_SHIFT, %rax
       addq	%rax, AORIG
#endif
	ALIGN_4

L(29):
#ifdef LN
       movq	K, %rax
       salq	$BASE_SHIFT, %rax
       leaq	(B, %rax, 2), B
#endif

#if defined(LT) || defined(RN)
	movq	BO, B
#endif

#ifdef RN
	addq	$2, KK
#endif

#ifdef RT
	subq	$2, KK
#endif

	decq	J
	jne	L(01)
	ALIGN_4

L(999):
	movq	  0(%rsp), %rbx
	movq	  8(%rsp), %rbp
	movq	 16(%rsp), %r12
	movq	 24(%rsp), %r13
	movq	 32(%rsp), %r14
	movq	 40(%rsp), %r15
	addq	$STACKSIZE, %rsp
	ret

	EPILOGUE
