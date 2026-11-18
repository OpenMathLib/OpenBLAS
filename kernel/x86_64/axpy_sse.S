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

#ifndef WINDOWS_ABI
#define M	ARG1
#define X	ARG4
#define INCX	ARG5
#define Y	ARG6
#define INCY	ARG2
#else
#define M	ARG1
#define X	ARG2
#define INCX	ARG3
#define Y	ARG4
#define INCY	%r10
#endif

#define	YY	%r11
#define ALPHA	%xmm15

#include "l1param.h"

	PROLOGUE
	PROFCODE

#ifndef WINDOWS_ABI
#ifndef XDOUBLE
	movq	 8(%rsp), INCY
#else
	movq	24(%rsp), INCY
#endif
	movaps	%xmm0,  ALPHA
#else


	movq	40(%rsp), X
	movq	48(%rsp), INCX
	movq	56(%rsp), Y
	movq	64(%rsp), INCY
#endif

	SAVEREGISTERS

#ifdef WINDOWS_ABI
	movaps	%xmm3,  ALPHA
#endif

	shufps	$0, ALPHA, ALPHA

	leaq	(, INCX, SIZE), INCX
	leaq	(, INCY, SIZE), INCY

	testq	M, M
	jle	L(19)

	cmpq	$SIZE, INCX
	jne	L(50)
	cmpq	$SIZE, INCY
	jne	L(50)

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y

	cmpq	$3, M
	jle	L(16)

	testq	$SIZE, Y
	je	L(00)

	movss	-32 * SIZE(X), %xmm0
	mulss	ALPHA, %xmm0
	addss	-32 * SIZE(Y), %xmm0
	movss	%xmm0, -32 * SIZE(Y)
	addq	$1 * SIZE, X
	addq	$1 * SIZE, Y
	decq	M
	jle	L(19)
	ALIGN_3

L(00):
	testq	$SIZE * 2, Y
	je	L(10)

	movsd	-32 * SIZE(X), %xmm0
	movsd	-32 * SIZE(Y), %xmm4
	mulps	ALPHA, %xmm0
	addps	%xmm4, %xmm0
	movsd	%xmm0, -32 * SIZE(Y)

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	subq	$2, M
	jle	L(19)
	ALIGN_3

L(10):
	testq	$SIZE * 3, X
	jne	L(20)

	movq	M,  %rax
	sarq	$5, %rax
	jle	L(13)

	movaps	-32 * SIZE(X), %xmm0
	movaps	-28 * SIZE(X), %xmm1
	movaps	-24 * SIZE(X), %xmm2
	movaps	-20 * SIZE(X), %xmm3

	decq	%rax
	jle L(12)
	ALIGN_4

L(11):
	movaps	-16 * SIZE(X), %xmm4
	movaps	-12 * SIZE(X), %xmm5

#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	movaps	%xmm0, -32 * SIZE(Y)

	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1
	movaps	%xmm1, -28 * SIZE(Y)

	movaps	 -8 * SIZE(X), %xmm6
	movaps	 -4 * SIZE(X), %xmm7

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	mulps	ALPHA, %xmm2
	addps	-24 * SIZE(Y), %xmm2
	movaps	%xmm2, -24 * SIZE(Y)

	mulps	ALPHA, %xmm3
	addps	-20 * SIZE(Y), %xmm3
	movaps	%xmm3, -20 * SIZE(Y)

	movaps	 0 * SIZE(X), %xmm0
	movaps	 4 * SIZE(X), %xmm1

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

	mulps	ALPHA, %xmm4
	addps	-16 * SIZE(Y), %xmm4
	movaps	%xmm4, -16 * SIZE(Y)

	mulps	ALPHA, %xmm5
	addps	-12 * SIZE(Y), %xmm5
	movaps	%xmm5, -12 * SIZE(Y)

	movaps	 8 * SIZE(X), %xmm2
	movaps	12 * SIZE(X), %xmm3

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif

	mulps	ALPHA, %xmm6
	addps	 -8 * SIZE(Y), %xmm6
	movaps	%xmm6,  -8 * SIZE(Y)

	mulps	ALPHA, %xmm7
	addps	 -4 * SIZE(Y), %xmm7
	movaps	%xmm7,  -4 * SIZE(Y)

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	decq	%rax
	jg	L(11)
	ALIGN_3

L(12):
	movaps	-16 * SIZE(X), %xmm4
	movaps	-12 * SIZE(X), %xmm5

	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	movaps	%xmm0, -32 * SIZE(Y)

	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1
	movaps	%xmm1, -28 * SIZE(Y)

	movaps	 -8 * SIZE(X), %xmm6
	movaps	 -4 * SIZE(X), %xmm7

	mulps	ALPHA, %xmm2
	addps	-24 * SIZE(Y), %xmm2
	movaps	%xmm2, -24 * SIZE(Y)

	mulps	ALPHA, %xmm3
	addps	-20 * SIZE(Y), %xmm3
	movaps	%xmm3, -20 * SIZE(Y)

	mulps	ALPHA, %xmm4
	addps	-16 * SIZE(Y), %xmm4
	movaps	%xmm4, -16 * SIZE(Y)

	mulps	ALPHA, %xmm5
	addps	-12 * SIZE(Y), %xmm5
	movaps	%xmm5, -12 * SIZE(Y)

	mulps	ALPHA, %xmm6
	addps	 -8 * SIZE(Y), %xmm6
	movaps	%xmm6,  -8 * SIZE(Y)

	mulps	ALPHA, %xmm7
	addps	 -4 * SIZE(Y), %xmm7
	movaps	%xmm7,  -4 * SIZE(Y)

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

L(13):
	movq	M,  %rax
	andq	$16, %rax
	jle	L(14)
	ALIGN_3

	movaps	-32 * SIZE(X), %xmm0
	movaps	-28 * SIZE(X), %xmm1
	movaps	-24 * SIZE(X), %xmm2
	movaps	-20 * SIZE(X), %xmm3

	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1
	mulps	ALPHA, %xmm2
	addps	-24 * SIZE(Y), %xmm2
	mulps	ALPHA, %xmm3
	addps	-20 * SIZE(Y), %xmm3

	movaps	%xmm0, -32 * SIZE(Y)
	movaps	%xmm1, -28 * SIZE(Y)
	movaps	%xmm2, -24 * SIZE(Y)
	movaps	%xmm3, -20 * SIZE(Y)

	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
	ALIGN_3

L(14):
	movq	M,  %rax
	andq	$8, %rax
	jle	L(15)
	ALIGN_3

	movaps	-32 * SIZE(X), %xmm0
	movaps	-28 * SIZE(X), %xmm1

	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1

	movaps	%xmm0, -32 * SIZE(Y)
	movaps	%xmm1, -28 * SIZE(Y)

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

L(15):
	movq	M,  %rax
	andq	$4, %rax
	jle	L(16)
	ALIGN_3

	movaps	-32 * SIZE(X), %xmm0

	mulps	ALPHA, %xmm0

	addps	-32 * SIZE(Y), %xmm0

	movaps	%xmm0, -32 * SIZE(Y)

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

L(16):
	movq	M,  %rax
	andq	$2, %rax
	jle	L(17)
	ALIGN_3

	movsd	-32 * SIZE(X), %xmm0
	movsd	-32 * SIZE(Y), %xmm4

	mulps	ALPHA, %xmm0
	addps	%xmm4, %xmm0

	movsd	%xmm0, -32 * SIZE(Y)

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

L(17):
	movq	M,  %rax
	andq	$1, %rax
	jle	L(19)
	ALIGN_3

	movss	-32 * SIZE(X), %xmm0
	mulss	ALPHA, %xmm0
	addss	-32 * SIZE(Y), %xmm0

	movss	%xmm0, 	-32 * SIZE(Y)
	ALIGN_3

L(19):
	xorq	%rax,%rax

	RESTOREREGISTERS

	ret
	ALIGN_3

L(20):

#ifdef ALIGNED_ACCESS

	testq	$SIZE, X
	jne	L(30)

	movhps	-32 * SIZE(X), %xmm0

	movq	M,  %rax
	sarq	$5, %rax
	jle	L(23)

	movaps	-30 * SIZE(X), %xmm1
	movaps	-26 * SIZE(X), %xmm2
	movaps	-22 * SIZE(X), %xmm3
	movaps	-18 * SIZE(X), %xmm4

	decq	%rax
	jle L(22)
	ALIGN_4

L(21):
	movaps	-14 * SIZE(X), %xmm5
	movaps	-10 * SIZE(X), %xmm6

#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	SHUFPD_1 %xmm1, %xmm0
	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	movaps	%xmm0, -32 * SIZE(Y)

	SHUFPD_1 %xmm2, %xmm1
	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1
	movaps	%xmm1, -28 * SIZE(Y)

	movaps	 -6 * SIZE(X), %xmm7
	movaps	 -2 * SIZE(X), %xmm0

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	SHUFPD_1 %xmm3, %xmm2
	mulps	ALPHA, %xmm2
	addps	-24 * SIZE(Y), %xmm2
	movaps	%xmm2, -24 * SIZE(Y)

	SHUFPD_1 %xmm4, %xmm3
	mulps	ALPHA, %xmm3
	addps	-20 * SIZE(Y), %xmm3
	movaps	%xmm3, -20 * SIZE(Y)

	movaps	  2 * SIZE(X), %xmm1
	movaps	  6 * SIZE(X), %xmm2

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

	SHUFPD_1 %xmm5, %xmm4
	mulps	ALPHA, %xmm4
	addps	-16 * SIZE(Y), %xmm4
	movaps	%xmm4, -16 * SIZE(Y)

	SHUFPD_1 %xmm6, %xmm5
	mulps	ALPHA, %xmm5
	addps	-12 * SIZE(Y), %xmm5
	movaps	%xmm5, -12 * SIZE(Y)

	movaps	 10 * SIZE(X), %xmm3
	movaps	 14 * SIZE(X), %xmm4

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif

	SHUFPD_1 %xmm7, %xmm6
	mulps	ALPHA, %xmm6
	addps	 -8 * SIZE(Y), %xmm6
	movaps	%xmm6,  -8 * SIZE(Y)

	SHUFPD_1 %xmm0, %xmm7
	mulps	ALPHA, %xmm7
	addps	 -4 * SIZE(Y), %xmm7
	movaps	%xmm7,  -4 * SIZE(Y)

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	decq	%rax
	jg	L(21)
	ALIGN_3

L(22):
	movaps	-14 * SIZE(X), %xmm5
	movaps	-10 * SIZE(X), %xmm6

	SHUFPD_1 %xmm1, %xmm0
	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	movaps	%xmm0, -32 * SIZE(Y)

	SHUFPD_1 %xmm2, %xmm1
	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1
	movaps	%xmm1, -28 * SIZE(Y)

	movaps	 -6 * SIZE(X), %xmm7
	movaps	 -2 * SIZE(X), %xmm0

	SHUFPD_1 %xmm3, %xmm2
	mulps	ALPHA, %xmm2
	addps	-24 * SIZE(Y), %xmm2
	movaps	%xmm2, -24 * SIZE(Y)

	SHUFPD_1 %xmm4, %xmm3
	mulps	ALPHA, %xmm3
	addps	-20 * SIZE(Y), %xmm3
	movaps	%xmm3, -20 * SIZE(Y)

	SHUFPD_1 %xmm5, %xmm4
	mulps	ALPHA, %xmm4
	addps	-16 * SIZE(Y), %xmm4
	movaps	%xmm4, -16 * SIZE(Y)

	SHUFPD_1 %xmm6, %xmm5
	mulps	ALPHA, %xmm5
	addps	-12 * SIZE(Y), %xmm5
	movaps	%xmm5, -12 * SIZE(Y)

	SHUFPD_1 %xmm7, %xmm6
	mulps	ALPHA, %xmm6
	addps	 -8 * SIZE(Y), %xmm6
	movaps	%xmm6,  -8 * SIZE(Y)

	SHUFPD_1 %xmm0, %xmm7
	mulps	ALPHA, %xmm7
	addps	 -4 * SIZE(Y), %xmm7
	movaps	%xmm7,  -4 * SIZE(Y)

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

L(23):
	movq	M,  %rax
	andq	$16, %rax
	jle	L(24)
	ALIGN_3

	movaps	-30 * SIZE(X), %xmm1
	movaps	-26 * SIZE(X), %xmm2
	movaps	-22 * SIZE(X), %xmm3
	movaps	-18 * SIZE(X), %xmm4

	SHUFPD_1 %xmm1, %xmm0
	SHUFPD_1 %xmm2, %xmm1
	SHUFPD_1 %xmm3, %xmm2
	SHUFPD_1 %xmm4, %xmm3

	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1
	mulps	ALPHA, %xmm2
	addps	-24 * SIZE(Y), %xmm2
	mulps	ALPHA, %xmm3
	addps	-20 * SIZE(Y), %xmm3

	movaps	%xmm0, -32 * SIZE(Y)
	movaps	%xmm1, -28 * SIZE(Y)
	movaps	%xmm2, -24 * SIZE(Y)
	movaps	%xmm3, -20 * SIZE(Y)

	movaps	%xmm4, %xmm0

	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
	ALIGN_3

L(24):
	movq	M,  %rax
	andq	$8, %rax
	jle	L(25)
	ALIGN_3

	movaps	-30 * SIZE(X), %xmm1
	movaps	-26 * SIZE(X), %xmm2

	SHUFPD_1 %xmm1, %xmm0
	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	SHUFPD_1 %xmm2, %xmm1
	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1

	movaps	%xmm0, -32 * SIZE(Y)
	movaps	%xmm1, -28 * SIZE(Y)
	movaps	%xmm2, %xmm0

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

L(25):
	movq	M,  %rax
	andq	$4, %rax
	jle	L(26)
	ALIGN_3

	movaps	-30 * SIZE(X), %xmm1

	SHUFPD_1 %xmm1, %xmm0
	mulps	ALPHA, %xmm0

	addps	-32 * SIZE(Y), %xmm0

	movaps	%xmm0, -32 * SIZE(Y)

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

L(26):
	movq	M,  %rax
	andq	$2, %rax
	jle	L(27)
	ALIGN_3

	movsd	-32 * SIZE(X), %xmm0
	movsd	-32 * SIZE(Y), %xmm4

	mulps	ALPHA, %xmm0
	addps	%xmm4, %xmm0

	movsd	%xmm0, -32 * SIZE(Y)

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

L(27):
	movq	M,  %rax
	andq	$1, %rax
	jle	L(29)
	ALIGN_3

	movss	-32 * SIZE(X), %xmm0
	mulss	ALPHA, %xmm0
	addss	-32 * SIZE(Y), %xmm0

	movss	%xmm0, 	-32 * SIZE(Y)
	addq	$SIZE, Y
	ALIGN_3

L(29):
	xorq	%rax,%rax

	RESTOREREGISTERS

	ret
	ALIGN_3

L(30):
	testq	$2 * SIZE, X
	jne	L(40)

	movaps	-33 * SIZE(X), %xmm0

	movq	M,  %rax
	sarq	$5, %rax
	jle	L(33)

	movaps	-29 * SIZE(X), %xmm1
	movaps	-25 * SIZE(X), %xmm2
	movaps	-21 * SIZE(X), %xmm3
	movaps	-17 * SIZE(X), %xmm4

	decq	%rax
	jle L(32)
	ALIGN_4

L(31):
	movaps	-13 * SIZE(X), %xmm5
	movaps	 -9 * SIZE(X), %xmm6

#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	movss	%xmm1, %xmm0
	SHUFPS_39 %xmm0, %xmm0
	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	movaps	%xmm0, -32 * SIZE(Y)

	movss	%xmm2, %xmm1
	SHUFPS_39 %xmm1, %xmm1
	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1
	movaps	%xmm1, -28 * SIZE(Y)

	movaps	 -5 * SIZE(X), %xmm7
	movaps	 -1 * SIZE(X), %xmm0

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movss	%xmm3, %xmm2
	SHUFPS_39 %xmm2, %xmm2
	mulps	ALPHA, %xmm2
	addps	-24 * SIZE(Y), %xmm2
	movaps	%xmm2, -24 * SIZE(Y)

	movss	%xmm4, %xmm3
	SHUFPS_39 %xmm3, %xmm3
	mulps	ALPHA, %xmm3
	addps	-20 * SIZE(Y), %xmm3
	movaps	%xmm3, -20 * SIZE(Y)

	movaps	  3 * SIZE(X), %xmm1
	movaps	  7 * SIZE(X), %xmm2

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

	movss	%xmm5, %xmm4
	SHUFPS_39 %xmm4, %xmm4
	mulps	ALPHA, %xmm4
	addps	-16 * SIZE(Y), %xmm4
	movaps	%xmm4, -16 * SIZE(Y)

	movss	%xmm6, %xmm5
	SHUFPS_39 %xmm5, %xmm5
	mulps	ALPHA, %xmm5
	addps	-12 * SIZE(Y), %xmm5
	movaps	%xmm5, -12 * SIZE(Y)

	movaps	 11 * SIZE(X), %xmm3
	movaps	 15 * SIZE(X), %xmm4

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif

	movss	%xmm7, %xmm6
	SHUFPS_39 %xmm6, %xmm6
	mulps	ALPHA, %xmm6
	addps	-8 * SIZE(Y), %xmm6
	movaps	%xmm6,  -8 * SIZE(Y)

	movss	%xmm0, %xmm7
	SHUFPS_39 %xmm7, %xmm7
	mulps	ALPHA, %xmm7
	addps	-4 * SIZE(Y), %xmm7
	movaps	%xmm7,  -4 * SIZE(Y)

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	decq	%rax
	jg	L(31)
	ALIGN_3

L(32):
	movaps	-13 * SIZE(X), %xmm5
	movaps	 -9 * SIZE(X), %xmm6

	movss	%xmm1, %xmm0
	SHUFPS_39 %xmm0, %xmm0
	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	movaps	%xmm0, -32 * SIZE(Y)

	movss	%xmm2, %xmm1
	SHUFPS_39 %xmm1, %xmm1
	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1
	movaps	%xmm1, -28 * SIZE(Y)

	movaps	 -5 * SIZE(X), %xmm7
	movaps	 -1 * SIZE(X), %xmm0

	movss	%xmm3, %xmm2
	SHUFPS_39 %xmm2, %xmm2
	mulps	ALPHA, %xmm2
	addps	-24 * SIZE(Y), %xmm2
	movaps	%xmm2, -24 * SIZE(Y)

	movss	%xmm4, %xmm3
	SHUFPS_39 %xmm3, %xmm3
	mulps	ALPHA, %xmm3
	addps	-20 * SIZE(Y), %xmm3
	movaps	%xmm3, -20 * SIZE(Y)

	movss	%xmm5, %xmm4
	SHUFPS_39 %xmm4, %xmm4
	mulps	ALPHA, %xmm4
	addps	-16 * SIZE(Y), %xmm4
	movaps	%xmm4, -16 * SIZE(Y)

	movss	%xmm6, %xmm5
	SHUFPS_39 %xmm5, %xmm5
	mulps	ALPHA, %xmm5
	addps	-12 * SIZE(Y), %xmm5
	movaps	%xmm5, -12 * SIZE(Y)

	movss	%xmm7, %xmm6
	SHUFPS_39 %xmm6, %xmm6
	mulps	ALPHA, %xmm6
	addps	-8 * SIZE(Y), %xmm6
	movaps	%xmm6,  -8 * SIZE(Y)

	movss	%xmm0, %xmm7
	SHUFPS_39 %xmm7, %xmm7
	mulps	ALPHA, %xmm7
	addps	-4 * SIZE(Y), %xmm7
	movaps	%xmm7,  -4 * SIZE(Y)

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

L(33):
	movq	M,  %rax
	andq	$16, %rax
	jle	L(34)
	ALIGN_3

	movaps	-29 * SIZE(X), %xmm1
	movaps	-25 * SIZE(X), %xmm2
	movaps	-21 * SIZE(X), %xmm3
	movaps	-17 * SIZE(X), %xmm4

	movss	%xmm1, %xmm0
	SHUFPS_39 %xmm0, %xmm0
	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0

	movss	%xmm2, %xmm1
	SHUFPS_39 %xmm1, %xmm1
	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1

	movss	%xmm3, %xmm2
	SHUFPS_39 %xmm2, %xmm2
	mulps	ALPHA, %xmm2
	addps	-24 * SIZE(Y), %xmm2

	movss	%xmm4, %xmm3
	SHUFPS_39 %xmm3, %xmm3
	mulps	ALPHA, %xmm3
	addps	-20 * SIZE(Y), %xmm3

	movaps	%xmm0, -32 * SIZE(Y)
	movaps	%xmm1, -28 * SIZE(Y)
	movaps	%xmm2, -24 * SIZE(Y)
	movaps	%xmm3, -20 * SIZE(Y)

	movaps	%xmm4, %xmm0

	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
	ALIGN_3

L(34):
	movq	M,  %rax
	andq	$8, %rax
	jle	L(35)
	ALIGN_3

	movaps	-29 * SIZE(X), %xmm1
	movaps	-25 * SIZE(X), %xmm2

	movss	%xmm1, %xmm0
	SHUFPS_39 %xmm0, %xmm0
	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0

	movss	%xmm2, %xmm1
	SHUFPS_39 %xmm1, %xmm1
	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1

	movaps	%xmm0, -32 * SIZE(Y)
	movaps	%xmm1, -28 * SIZE(Y)
	movaps	%xmm2, %xmm0

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

L(35):
	movq	M,  %rax
	andq	$4, %rax
	jle	L(36)
	ALIGN_3

	movaps	-29 * SIZE(X), %xmm1

	movss	%xmm1, %xmm0
	SHUFPS_39 %xmm0, %xmm0
	mulps	ALPHA, %xmm0

	addps	-32 * SIZE(Y), %xmm0

	movaps	%xmm0, -32 * SIZE(Y)

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

L(36):
	movq	M,  %rax
	andq	$2, %rax
	jle	L(37)
	ALIGN_3

	movsd	-32 * SIZE(X), %xmm0
	movsd	-32 * SIZE(Y), %xmm4

	mulps	ALPHA, %xmm0
	addps	%xmm4, %xmm0

	movsd	%xmm0, -32 * SIZE(Y)

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

L(37):
	movq	M,  %rax
	andq	$1, %rax
	jle	L(39)
	ALIGN_3

	movss	-32 * SIZE(X), %xmm0
	mulss	ALPHA, %xmm0
	addss	-32 * SIZE(Y), %xmm0

	movss	%xmm0, 	-32 * SIZE(Y)
	addq	$SIZE, Y
	ALIGN_3

L(39):
	xorq	%rax,%rax

	RESTOREREGISTERS

	ret
	ALIGN_3

L(40):
	movaps	-35 * SIZE(X), %xmm0

	movq	M,  %rax
	sarq	$5, %rax
	jle	L(43)

	movaps	-31 * SIZE(X), %xmm1
	movaps	-27 * SIZE(X), %xmm2
	movaps	-23 * SIZE(X), %xmm3
	movaps	-19 * SIZE(X), %xmm4

	decq	%rax
	jle L(42)
	ALIGN_4

L(41):
	movaps	-15 * SIZE(X), %xmm5
	movaps	-11 * SIZE(X), %xmm6

#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	movss	%xmm1, %xmm0
	shufps	$0x93, %xmm1, %xmm0
	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	movaps	%xmm0, -32 * SIZE(Y)

	movss	%xmm2, %xmm1
	shufps	$0x93, %xmm2, %xmm1
	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1
	movaps	%xmm1, -28 * SIZE(Y)

	movaps	 -7 * SIZE(X), %xmm7
	movaps	 -3 * SIZE(X), %xmm0

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movss	%xmm3, %xmm2
	shufps	$0x93, %xmm3, %xmm2
	mulps	ALPHA, %xmm2
	addps	-24 * SIZE(Y), %xmm2
	movaps	%xmm2, -24 * SIZE(Y)

	movss	%xmm4, %xmm3
	shufps	$0x93, %xmm4, %xmm3
	mulps	ALPHA, %xmm3
	addps	-20 * SIZE(Y), %xmm3
	movaps	%xmm3, -20 * SIZE(Y)

	movaps	  1 * SIZE(X), %xmm1
	movaps	  5 * SIZE(X), %xmm2

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

	movss	%xmm5, %xmm4
	shufps	$0x93, %xmm5, %xmm4
	mulps	ALPHA, %xmm4
	addps	-16 * SIZE(Y), %xmm4
	movaps	%xmm4, -16 * SIZE(Y)

	movss	%xmm6, %xmm5
	shufps	$0x93, %xmm6, %xmm5
	mulps	ALPHA, %xmm5
	addps	-12 * SIZE(Y), %xmm5
	movaps	%xmm5, -12 * SIZE(Y)

	movaps	  9 * SIZE(X), %xmm3
	movaps	 13 * SIZE(X), %xmm4

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif

	movss	%xmm7, %xmm6
	shufps	$0x93, %xmm7, %xmm6
	mulps	ALPHA, %xmm6
	addps	 -8 * SIZE(Y), %xmm6
	movaps	%xmm6, -8 * SIZE(Y)

	movss	%xmm0, %xmm7
	shufps	$0x93, %xmm0, %xmm7
	mulps	ALPHA, %xmm7
	addps	 -4 * SIZE(Y), %xmm7
	movaps	%xmm7, -4 * SIZE(Y)

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	decq	%rax
	jg	L(41)
	ALIGN_3

L(42):
	movaps	-15 * SIZE(X), %xmm5
	movaps	-11 * SIZE(X), %xmm6

	movss	%xmm1, %xmm0
	shufps	$0x93, %xmm1, %xmm0
	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	movaps	%xmm0, -32 * SIZE(Y)

	movss	%xmm2, %xmm1
	shufps	$0x93, %xmm2, %xmm1
	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1
	movaps	%xmm1, -28 * SIZE(Y)

	movaps	 -7 * SIZE(X), %xmm7
	movaps	 -3 * SIZE(X), %xmm0

	movss	%xmm3, %xmm2
	shufps	$0x93, %xmm3, %xmm2
	mulps	ALPHA, %xmm2
	addps	-24 * SIZE(Y), %xmm2
	movaps	%xmm2, -24 * SIZE(Y)

	movss	%xmm4, %xmm3
	shufps	$0x93, %xmm4, %xmm3
	mulps	ALPHA, %xmm3
	addps	-20 * SIZE(Y), %xmm3
	movaps	%xmm3, -20 * SIZE(Y)

	movss	%xmm5, %xmm4
	shufps	$0x93, %xmm5, %xmm4
	mulps	ALPHA, %xmm4
	addps	-16 * SIZE(Y), %xmm4
	movaps	%xmm4, -16 * SIZE(Y)

	movss	%xmm6, %xmm5
	shufps	$0x93, %xmm6, %xmm5
	mulps	ALPHA, %xmm5
	addps	-12 * SIZE(Y), %xmm5
	movaps	%xmm5, -12 * SIZE(Y)

	movss	%xmm7, %xmm6
	shufps	$0x93, %xmm7, %xmm6
	mulps	ALPHA, %xmm6
	addps	 -8 * SIZE(Y), %xmm6
	movaps	%xmm6, -8 * SIZE(Y)

	movss	%xmm0, %xmm7
	shufps	$0x93, %xmm0, %xmm7
	mulps	ALPHA, %xmm7
	addps	 -4 * SIZE(Y), %xmm7
	movaps	%xmm7, -4 * SIZE(Y)

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

L(43):
	movq	M,  %rax
	andq	$16, %rax
	jle	L(44)
	ALIGN_3

	movaps	-31 * SIZE(X), %xmm1
	movaps	-27 * SIZE(X), %xmm2
	movaps	-23 * SIZE(X), %xmm3
	movaps	-19 * SIZE(X), %xmm4

	movss	%xmm1, %xmm0
	shufps	$0x93, %xmm1, %xmm0
	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	movss	%xmm2, %xmm1
	shufps	$0x93, %xmm2, %xmm1
	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1

	movss	%xmm3, %xmm2
	shufps	$0x93, %xmm3, %xmm2
	mulps	ALPHA, %xmm2
	addps	-24 * SIZE(Y), %xmm2
	movss	%xmm4, %xmm3
	shufps	$0x93, %xmm4, %xmm3
	mulps	ALPHA, %xmm3
	addps	-20 * SIZE(Y), %xmm3

	movaps	%xmm0, -32 * SIZE(Y)
	movaps	%xmm1, -28 * SIZE(Y)
	movaps	%xmm2, -24 * SIZE(Y)
	movaps	%xmm3, -20 * SIZE(Y)

	movaps	%xmm4, %xmm0

	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
	ALIGN_3

L(44):
	movq	M,  %rax
	andq	$8, %rax
	jle	L(45)
	ALIGN_3

	movaps	-31 * SIZE(X), %xmm1
	movaps	-27 * SIZE(X), %xmm2

	movss	%xmm1, %xmm0
	shufps	$0x93, %xmm1, %xmm0
	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	movss	%xmm2, %xmm1
	shufps	$0x93, %xmm2, %xmm1
	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1

	movaps	%xmm0, -32 * SIZE(Y)
	movaps	%xmm1, -28 * SIZE(Y)
	movaps	%xmm2, %xmm0

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

L(45):
	movq	M,  %rax
	andq	$4, %rax
	jle	L(46)
	ALIGN_3

	movaps	-31 * SIZE(X), %xmm1

	movss	%xmm1, %xmm0
	shufps	$0x93, %xmm1, %xmm0
	mulps	ALPHA, %xmm0

	addps	-32 * SIZE(Y), %xmm0

	movaps	%xmm0, -32 * SIZE(Y)

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

L(46):
	movq	M,  %rax
	andq	$2, %rax
	jle	L(47)
	ALIGN_3

	movsd	-32 * SIZE(X), %xmm0
	movsd	-32 * SIZE(Y), %xmm4

	mulps	ALPHA, %xmm0
	addps	%xmm4, %xmm0

	movsd	%xmm0, -32 * SIZE(Y)

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

L(47):
	movq	M,  %rax
	andq	$1, %rax
	jle	L(49)
	ALIGN_3

	movss	-32 * SIZE(X), %xmm0
	mulss	ALPHA, %xmm0
	addss	-32 * SIZE(Y), %xmm0

	movss	%xmm0, 	-32 * SIZE(Y)
	addq	$SIZE, Y
	ALIGN_3

L(49):
	xorq	%rax,%rax

	RESTOREREGISTERS

	ret

#else

	movq	M,  %rax
	sarq	$5, %rax
	jle	L(23)

	movsd	-32 * SIZE(X), %xmm0
	movhps	-30 * SIZE(X), %xmm0
	movsd	-28 * SIZE(X), %xmm1
	movhps	-26 * SIZE(X), %xmm1
	movsd	-24 * SIZE(X), %xmm2
	movhps	-22 * SIZE(X), %xmm2
	movsd	-20 * SIZE(X), %xmm3
	movhps	-18 * SIZE(X), %xmm3

	decq	%rax
	jle L(22)
	ALIGN_4

L(21):
	movsd	-16 * SIZE(X), %xmm4
	movhps	-14 * SIZE(X), %xmm4
	movsd	-12 * SIZE(X), %xmm5
	movhps	-10 * SIZE(X), %xmm5

#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	movaps	%xmm0, -32 * SIZE(Y)

	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1
	movaps	%xmm1, -28 * SIZE(Y)

	movsd	 -8 * SIZE(X), %xmm6
	movhps	 -6 * SIZE(X), %xmm6
	movsd	 -4 * SIZE(X), %xmm7
	movhps	 -2 * SIZE(X), %xmm7

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	mulps	ALPHA, %xmm2
	addps	-24 * SIZE(Y), %xmm2
	movaps	%xmm2, -24 * SIZE(Y)

	mulps	ALPHA, %xmm3
	addps	-20 * SIZE(Y), %xmm3
	movaps	%xmm3, -20 * SIZE(Y)

	movsd	  0 * SIZE(X), %xmm0
	movhps	  2 * SIZE(X), %xmm0
	movsd	  4 * SIZE(X), %xmm1
	movhps	  6 * SIZE(X), %xmm1

#if defined(PREFETCHW) && !defined(FETCH128)
	PREFETCHW (PREFETCHSIZE +  64) - PREOFFSET(Y)
#endif

	mulps	ALPHA, %xmm4
	addps	-16 * SIZE(Y), %xmm4
	movaps	%xmm4, -16 * SIZE(Y)

	mulps	ALPHA, %xmm5
	addps	-12 * SIZE(Y), %xmm5
	movaps	%xmm5, -12 * SIZE(Y)

	movsd	  8 * SIZE(X), %xmm2
	movhps	 10 * SIZE(X), %xmm2
	movsd	 12 * SIZE(X), %xmm3
	movhps	 14 * SIZE(X), %xmm3

#if defined(PREFETCH) && !defined(FETCH128)
	PREFETCH (PREFETCHSIZE +  64) - PREOFFSET(X)
#endif

	mulps	ALPHA, %xmm6
	addps	 -8 * SIZE(Y), %xmm6
	movaps	%xmm6,  -8 * SIZE(Y)

	mulps	ALPHA, %xmm7
	addps	 -4 * SIZE(Y), %xmm7
	movaps	%xmm7,  -4 * SIZE(Y)

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	decq	%rax
	jg	L(21)
	ALIGN_3

L(22):
	movsd	-16 * SIZE(X), %xmm4
	movhps	-14 * SIZE(X), %xmm4
	movsd	-12 * SIZE(X), %xmm5
	movhps	-10 * SIZE(X), %xmm5

	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	movaps	%xmm0, -32 * SIZE(Y)

	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1
	movaps	%xmm1, -28 * SIZE(Y)

	movsd	 -8 * SIZE(X), %xmm6
	movhps	 -6 * SIZE(X), %xmm6
	movsd	 -4 * SIZE(X), %xmm7
	movhps	 -2 * SIZE(X), %xmm7

	mulps	ALPHA, %xmm2
	addps	-24 * SIZE(Y), %xmm2
	movaps	%xmm2, -24 * SIZE(Y)

	mulps	ALPHA, %xmm3
	addps	-20 * SIZE(Y), %xmm3
	movaps	%xmm3, -20 * SIZE(Y)

	mulps	ALPHA, %xmm4
	addps	-16 * SIZE(Y), %xmm4
	movaps	%xmm4, -16 * SIZE(Y)

	mulps	ALPHA, %xmm5
	addps	-12 * SIZE(Y), %xmm5
	movaps	%xmm5, -12 * SIZE(Y)

	mulps	ALPHA, %xmm6
	addps	 -8 * SIZE(Y), %xmm6
	movaps	%xmm6,  -8 * SIZE(Y)

	mulps	ALPHA, %xmm7
	addps	 -4 * SIZE(Y), %xmm7
	movaps	%xmm7,  -4 * SIZE(Y)

	subq	$-32 * SIZE, X
	subq	$-32 * SIZE, Y
	ALIGN_3

L(23):
	movq	M,  %rax
	andq	$16, %rax
	jle	L(24)
	ALIGN_3

	movsd	-32 * SIZE(X), %xmm0
	movhps	-30 * SIZE(X), %xmm0
	movsd	-28 * SIZE(X), %xmm1
	movhps	-26 * SIZE(X), %xmm1

	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	movaps	%xmm0, -32 * SIZE(Y)
	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1
	movaps	%xmm1, -28 * SIZE(Y)

	movsd	-24 * SIZE(X), %xmm2
	movhps	-22 * SIZE(X), %xmm2
	movsd	-20 * SIZE(X), %xmm3
	movhps	-18 * SIZE(X), %xmm3

	mulps	ALPHA, %xmm2
	addps	-24 * SIZE(Y), %xmm2
	movaps	%xmm2, -24 * SIZE(Y)
	mulps	ALPHA, %xmm3
	addps	-20 * SIZE(Y), %xmm3
	movaps	%xmm3, -20 * SIZE(Y)

	addq	$16 * SIZE, X
	addq	$16 * SIZE, Y
	ALIGN_3

L(24):
	movq	M,  %rax
	andq	$8, %rax
	jle	L(25)
	ALIGN_3

	movsd	-32 * SIZE(X), %xmm0
	movhps	-30 * SIZE(X), %xmm0
	movsd	-28 * SIZE(X), %xmm1
	movhps	-26 * SIZE(X), %xmm1

	mulps	ALPHA, %xmm0
	addps	-32 * SIZE(Y), %xmm0
	mulps	ALPHA, %xmm1
	addps	-28 * SIZE(Y), %xmm1

	movaps	%xmm0, -32 * SIZE(Y)
	movaps	%xmm1, -28 * SIZE(Y)

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

L(25):
	movq	M,  %rax
	andq	$4, %rax
	jle	L(26)
	ALIGN_3

	movsd	-32 * SIZE(X), %xmm0
	movhps	-30 * SIZE(X), %xmm0

	mulps	ALPHA, %xmm0

	addps	-32 * SIZE(Y), %xmm0

	movaps	%xmm0, -32 * SIZE(Y)

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

L(26):
	movq	M,  %rax
	andq	$2, %rax
	jle	L(27)
	ALIGN_3

	movsd	-32 * SIZE(X), %xmm0
	movsd	-32 * SIZE(Y), %xmm4

	mulps	ALPHA, %xmm0
	addps	%xmm4, %xmm0

	movsd	%xmm0, -32 * SIZE(Y)

	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

L(27):
	movq	M,  %rax
	andq	$1, %rax
	jle	L(29)
	ALIGN_3

	movss	-32 * SIZE(X), %xmm0
	mulss	ALPHA, %xmm0
	addss	-32 * SIZE(Y), %xmm0

	movss	%xmm0, -32 * SIZE(Y)
	addq	$SIZE, Y
	ALIGN_3

L(29):
	xorq	%rax,%rax

	RESTOREREGISTERS

	ret
#endif
	ALIGN_3


L(50):
	movq	M, %rax
	movq	Y, YY
//If incx==0 || incy==0, avoid unloop.
	cmpq	$0, INCX
	je  L(56)
	cmpq	$0, INCY
	je  L(56)

	sarq	$3,   %rax
	jle	L(55)
	ALIGN_3

L(51):
	movss	(X), %xmm0
	addq	INCX, X
	mulss	ALPHA, %xmm0
	movss	(YY), %xmm6
	addq	INCY, YY
	addss	%xmm6, %xmm0

	movss	(X), %xmm1
	addq	INCX, X
	mulss	ALPHA, %xmm1
	movss	(YY), %xmm6
	addq	INCY, YY
	addss	%xmm6, %xmm1

	movss	(X), %xmm2
	addq	INCX, X
	mulss	ALPHA, %xmm2
	movss	(YY), %xmm6
	addq	INCY, YY
	addss	%xmm6, %xmm2

	movss	(X), %xmm3
	addq	INCX, X
	mulss	ALPHA, %xmm3
	movss	(YY), %xmm6
	addq	INCY, YY
	addss	%xmm6, %xmm3

	movss	%xmm0, (Y)
	addq	INCY, Y
	movss	%xmm1, (Y)
	addq	INCY, Y
	movss	%xmm2, (Y)
	addq	INCY, Y
	movss	%xmm3, (Y)
	addq	INCY, Y

	movss	(X), %xmm0
	addq	INCX, X
	mulss	ALPHA, %xmm0
	movss	(YY), %xmm6
	addq	INCY, YY
	addss	%xmm6, %xmm0

	movss	(X), %xmm1
	addq	INCX, X
	mulss	ALPHA, %xmm1
	movss	(YY), %xmm6
	addq	INCY, YY
	addss	%xmm6, %xmm1

	movss	(X), %xmm2
	addq	INCX, X
	mulss	ALPHA, %xmm2
	movss	(YY), %xmm6
	addq	INCY, YY
	addss	%xmm6, %xmm2

	movss	(X), %xmm3
	addq	INCX, X
	mulss	ALPHA, %xmm3
	movss	(YY), %xmm6
	addq	INCY, YY
	addss	%xmm6, %xmm3

	movss	%xmm0, (Y)
	addq	INCY, Y
	movss	%xmm1, (Y)
	addq	INCY, Y
	movss	%xmm2, (Y)
	addq	INCY, Y
	movss	%xmm3, (Y)
	addq	INCY, Y

	decq	%rax
	jg	L(51)
	ALIGN_3

L(55):
	movq	M, %rax
	andq	$7,   %rax
	jle	L(59)
	ALIGN_3

L(56):
	movss	(X), %xmm0
	addq	INCX, X
	mulss	ALPHA, %xmm0
	movss	(Y), %xmm6
	addss	%xmm6, %xmm0
	movss	%xmm0, (Y)
	addq	INCY, Y
	decq	%rax
	jg	L(56)
	ALIGN_3

L(59):
	xorq	%rax,%rax

	RESTOREREGISTERS

	ret
	ALIGN_3


	EPILOGUE
