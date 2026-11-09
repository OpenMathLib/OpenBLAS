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
#else
#define M	ARG1
#define X	ARG2
#define INCX	ARG3
#endif
#define FLAG    %r9

#define XX	%r10
#define I	%rax

#include "l1param.h"

	PROLOGUE
	PROFCODE

#ifdef WINDOWS_ABI
	movq	40(%rsp), X
	movq	48(%rsp), INCX
	movq	64(%rsp), FLAG
	movaps	%xmm3, %xmm0
#else
	movq    24(%rsp), FLAG
#endif

	SAVEREGISTERS

	testq	M, M
	jle	L(999)

	leaq	(, INCX, SIZE), INCX

	xorps	%xmm1, %xmm1
	comisd	%xmm0, %xmm1
	jne	L(100)		# Alpha != ZERO
	jp	L(100)		# For Alpha = NaN

	cmpq $1, FLAG
	je 	L(100)		# disable the Alpha=zero path as it does not handle x=inf or nan
/* Alpha == ZERO */
	cmpq	$SIZE, INCX
	jne	L(50)

/* INCX == 1 */
	testq	$15, X		# aligned for quad word?
	je	L(05)

	movsd	%xmm1, 0 * SIZE(X)
	addq	$SIZE, X
	decq	M
	jle	L(999)
	ALIGN_3
L(05):

/* Aligned Mode */
	movq	M,  I	# rcx = n
	sarq	$4, I
	jle	L(12)
	ALIGN_4

L(11):
#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movaps	%xmm1,  0 * SIZE(X)
	movaps	%xmm1,  2 * SIZE(X)
	movaps	%xmm1,  4 * SIZE(X)
	movaps	%xmm1,  6 * SIZE(X)

#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(X)
#endif

	movaps	%xmm1,  8 * SIZE(X)
	movaps	%xmm1, 10 * SIZE(X)
	movaps	%xmm1, 12 * SIZE(X)
	movaps	%xmm1, 14 * SIZE(X)

	addq	$16 * SIZE, X
	decq	I
	jg	L(11)
	ALIGN_4

L(12):
	testq	$15, M
	je	L(999)
	testq	$8, M
	je	L(13)

	movaps	%xmm1,  0 * SIZE(X)
	movaps	%xmm1,  2 * SIZE(X)
	movaps	%xmm1,  4 * SIZE(X)
	movaps	%xmm1,  6 * SIZE(X)
	addq	$8 * SIZE, X
	ALIGN_3

L(13):
	testq  $4, M
	je    L(14)

	movaps	%xmm1,  0 * SIZE(X)
	movaps	%xmm1,  2 * SIZE(X)
	addq	$4 * SIZE, X
	ALIGN_3

L(14):
	testq  $2, M
	je    L(15)

	movaps	%xmm1,  0 * SIZE(X)
	addq	$2 * SIZE, X
	ALIGN_3

L(15):
	testq  $1, M
	je    L(999)

	movsd	%xmm1,  0 * SIZE(X)
	jmp	L(999)
	ALIGN_4

L(50):
	movq	M,  I
	sarq	$3, I
	jle	L(52)
	ALIGN_4

L(51):
#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movsd	%xmm1, (X)
	addq	INCX, X
	movsd	%xmm1, (X)
	addq	INCX, X
	movsd	%xmm1, (X)
	addq	INCX, X
	movsd	%xmm1, (X)
	addq	INCX, X
	movsd	%xmm1, (X)
	addq	INCX, X
	movsd	%xmm1, (X)
	addq	INCX, X
	movsd	%xmm1, (X)
	addq	INCX, X
	movsd	%xmm1, (X)
	addq	INCX, X

	decq	I
	jg	L(51)
	ALIGN_4

L(52):
	testq	$7, M
	je	L(999)

	testq	$4, M
	je	L(53)

	movsd	%xmm1, (X)
	addq	INCX, X
	movsd	%xmm1, (X)
	addq	INCX, X
	movsd	%xmm1, (X)
	addq	INCX, X
	movsd	%xmm1, (X)
	addq	INCX, X
	ALIGN_3

L(53):
	testq	$2, M
	je	L(54)

	movsd	%xmm1, (X)
	addq	INCX, X
	movsd	%xmm1, (X)
	addq	INCX, X
	ALIGN_3

L(54):
	testq	$1, M
	je	L(999)

	movsd	%xmm1, (X)
	jmp	L(999)
	ALIGN_4

/* Alpha != ZERO */

L(100):
	unpcklpd %xmm0, %xmm0

	cmpq	$SIZE, INCX
	jne	L(150)

	testq	$SIZE, X
	je	L(105)

	movsd	0 * SIZE(X), %xmm1
	mulsd	%xmm0, %xmm1
	movsd	%xmm1, 0 * SIZE(X)
	addq	$SIZE, X
	decq	M
	jle	L(999)
	ALIGN_3
L(105):
	subq	$-16 * SIZE, X

	movq	M,  I		# rcx = n
	sarq	$4, I
	jle	L(113)

#if defined(BARCELONA) || defined(SHANGHAI) || defined(BULLDOZER)

	movaps  %xmm0, %xmm1
	mulpd	-16 * SIZE(X), %xmm1
	movaps  %xmm0, %xmm2
	mulpd	-14 * SIZE(X), %xmm2
	movaps  %xmm0, %xmm3
	mulpd	-12 * SIZE(X), %xmm3
	movaps  %xmm0, %xmm4
	mulpd	-10 * SIZE(X), %xmm4
	movaps  %xmm0, %xmm5
	mulpd	 -8 * SIZE(X), %xmm5
	movaps  %xmm0, %xmm6
	mulpd	 -6 * SIZE(X), %xmm6
	movaps  %xmm0, %xmm7
	mulpd	 -4 * SIZE(X), %xmm7
	movaps  %xmm0, %xmm8
	mulpd	 -2 * SIZE(X), %xmm8

	decq	I
	jle	L(112)
	ALIGN_4

L(111):
#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movaps	%xmm1, -16 * SIZE(X)
	movaps	%xmm2, -14 * SIZE(X)
	movaps	%xmm3, -12 * SIZE(X)
	movaps	%xmm4, -10 * SIZE(X)

	movaps  %xmm0, %xmm1
	mulpd	  0 * SIZE(X), %xmm1
	movaps  %xmm0, %xmm2
	mulpd	  2 * SIZE(X), %xmm2
	movaps  %xmm0, %xmm3
	mulpd	  4 * SIZE(X), %xmm3
	movaps  %xmm0, %xmm4
	mulpd	  6 * SIZE(X), %xmm4

#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(X)
#endif

	movaps	%xmm5,  -8 * SIZE(X)
	movaps	%xmm6,  -6 * SIZE(X)
	movaps	%xmm7,  -4 * SIZE(X)
	movaps	%xmm8,  -2 * SIZE(X)

	movaps  %xmm0, %xmm5
	mulpd	  8 * SIZE(X), %xmm5
	movaps  %xmm0, %xmm6
	mulpd	 10 * SIZE(X), %xmm6
	movaps  %xmm0, %xmm7
	mulpd	 12 * SIZE(X), %xmm7
	movaps  %xmm0, %xmm8
	mulpd	 14 * SIZE(X), %xmm8

	subq	$-16 * SIZE, X
	decq	I
	jg	L(111)
	ALIGN_4

L(112):
	movaps	%xmm1, -16 * SIZE(X)
	movaps	%xmm2, -14 * SIZE(X)
	movaps	%xmm3, -12 * SIZE(X)
	movaps	%xmm4, -10 * SIZE(X)
	movaps	%xmm5,  -8 * SIZE(X)
	movaps	%xmm6,  -6 * SIZE(X)
	movaps	%xmm7,  -4 * SIZE(X)
	movaps	%xmm8,  -2 * SIZE(X)

#else
	movaps	-16 * SIZE(X), %xmm1
	movaps	-14 * SIZE(X), %xmm2
	movaps	-12 * SIZE(X), %xmm3
	movaps	-10 * SIZE(X), %xmm4
	movaps	 -8 * SIZE(X), %xmm5
	movaps	 -6 * SIZE(X), %xmm6
	movaps	 -4 * SIZE(X), %xmm7
	movaps	 -2 * SIZE(X), %xmm8

	decq	I
	jle	L(112)
	ALIGN_4

L(111):
#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	mulpd	%xmm0, %xmm1
	movaps	%xmm1, -16 * SIZE(X)
	movaps	 0 * SIZE(X), %xmm1
	mulpd	%xmm0, %xmm2
	movaps	%xmm2, -14 * SIZE(X)
	movaps	 2 * SIZE(X), %xmm2

	mulpd	%xmm0, %xmm3
	movaps	%xmm3, -12 * SIZE(X)
	movaps	 4  * SIZE(X), %xmm3
	mulpd	%xmm0, %xmm4
	movaps	%xmm4, -10 * SIZE(X)
	movaps	 6 * SIZE(X), %xmm4

#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(X)
#endif

	mulpd	%xmm0, %xmm5
	movaps	%xmm5, -8 * SIZE(X)
	movaps	 8 * SIZE(X), %xmm5
	mulpd	%xmm0, %xmm6
	movaps	%xmm6, -6 * SIZE(X)
	movaps	10 * SIZE(X), %xmm6

	mulpd	%xmm0, %xmm7
	movaps	%xmm7, -4 * SIZE(X)
	movaps	12 * SIZE(X), %xmm7
	mulpd	%xmm0, %xmm8
	movaps	%xmm8, -2 * SIZE(X)
	movaps	14 * SIZE(X), %xmm8

	subq	$-16 * SIZE, X
	decq	I
	jg	L(111)
	ALIGN_4

L(112):
	mulpd	%xmm0, %xmm1
	movaps	%xmm1, -16 * SIZE(X)
	mulpd	%xmm0, %xmm2
	movaps	%xmm2, -14 * SIZE(X)
	mulpd	%xmm0, %xmm3
	movaps	%xmm3, -12 * SIZE(X)
	mulpd	%xmm0, %xmm4
	movaps	%xmm4, -10 * SIZE(X)

	mulpd	%xmm0, %xmm5
	movaps	%xmm5,  -8 * SIZE(X)
	mulpd	%xmm0, %xmm6
	movaps	%xmm6,  -6 * SIZE(X)
	mulpd	%xmm0, %xmm7
	movaps	%xmm7,  -4 * SIZE(X)
	mulpd	%xmm0, %xmm8
	movaps	%xmm8,  -2 * SIZE(X)
#endif

	subq	$-16 * SIZE, X
	ALIGN_3

L(113):
	testq	$15, M
	je	L(999)

	testq	$8, M
	je	L(114)

	movaps	-16 * SIZE(X), %xmm1
	movaps	-14 * SIZE(X), %xmm2
	movaps	-12 * SIZE(X), %xmm3
	movaps	-10 * SIZE(X), %xmm4

	mulpd	%xmm0, %xmm1
	movaps	%xmm1, -16 * SIZE(X)
	mulpd	%xmm0, %xmm2
	movaps	%xmm2, -14 * SIZE(X)
	mulpd	%xmm0, %xmm3
	movaps	%xmm3, -12 * SIZE(X)
	mulpd	%xmm0, %xmm4
	movaps	%xmm4, -10 * SIZE(X)
	addq	$8 * SIZE, X
	ALIGN_3

L(114):
	testq	$4, M
	je	L(115)

	movaps	-16 * SIZE(X), %xmm1
	movaps	-14 * SIZE(X), %xmm2

	mulpd	%xmm0, %xmm1
	movaps	%xmm1, -16 * SIZE(X)
	mulpd	%xmm0, %xmm2
	movaps	%xmm2, -14 * SIZE(X)
	addq	$4 * SIZE, X
	ALIGN_3

L(115):
	testq	$2, M
	je	L(116)

	movaps	-16 * SIZE(X), %xmm1
	mulpd	%xmm0, %xmm1
	movaps	%xmm1, -16 * SIZE(X)
	addq	$2 * SIZE, X
	ALIGN_3

L(116):
	testq	$1, M
	je	L(999)

	movsd	-16 * SIZE(X), %xmm1
	mulsd	%xmm0, %xmm1
	movsd	%xmm1, -16 * SIZE(X)
	jmp	L(999)
	ALIGN_3

/* incx != 1 */

L(150):
	movq	X, XX
	movq	M,  I		# rcx = n
	sarq	$3, I		# (n >> 3)
	jle	L(152)
	ALIGN_4

L(151):
#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE + 64) - PREOFFSET(X)
#endif
	movsd	(X), %xmm1
	addq	INCX, X
	movsd	(X), %xmm2
	addq	INCX, X
	movsd	(X), %xmm3
	addq	INCX, X
	movsd	(X), %xmm4
	addq	INCX, X
	movsd	(X), %xmm5
	addq	INCX, X
	movsd	(X), %xmm6
	addq	INCX, X
	movsd	(X), %xmm7
	addq	INCX, X
	movsd	(X), %xmm8
	addq	INCX, X

	mulsd	%xmm0, %xmm1
	mulsd	%xmm0, %xmm2
	mulsd	%xmm0, %xmm3
	mulsd	%xmm0, %xmm4
	mulsd	%xmm0, %xmm5
	mulsd	%xmm0, %xmm6
	mulsd	%xmm0, %xmm7
	mulsd	%xmm0, %xmm8

	movsd	%xmm1, (XX)
	addq	INCX, XX
	movsd	%xmm2, (XX)
	addq	INCX, XX
	movsd	%xmm3, (XX)
	addq	INCX, XX
	movsd	%xmm4, (XX)
	addq	INCX, XX
	movsd	%xmm5, (XX)
	addq	INCX, XX
	movsd	%xmm6, (XX)
	addq	INCX, XX
	movsd	%xmm7, (XX)
	addq	INCX, XX
	movsd	%xmm8, (XX)
	addq	INCX, XX
	decq	I
	jg	L(151)
	ALIGN_4

L(152):
	testq	$7, M
	je	L(999)

	testq	$4, M
	je	L(153)

	movsd	 (X), %xmm1
	addq	INCX, X
	movsd	 (X), %xmm2
	addq	INCX, X
	movsd	 (X), %xmm3
	addq	INCX, X
	movsd	 (X), %xmm4
	addq	INCX, X

	mulsd	%xmm0, %xmm1
	mulsd	%xmm0, %xmm2
	mulsd	%xmm0, %xmm3
	mulsd	%xmm0, %xmm4

	movsd	%xmm1,  (XX)
	addq	INCX, XX
	movsd	%xmm2,  (XX)
	addq	INCX, XX
	movsd	%xmm3,  (XX)
	addq	INCX, XX
	movsd	%xmm4,  (XX)
	addq	INCX, XX
	ALIGN_3

L(153):
	testq	$2, M
	je	L(154)

	movsd	 (X), %xmm1
	addq	INCX, X
	movsd	 (X), %xmm2
	addq	INCX, X

	mulsd	%xmm0, %xmm1
	mulsd	%xmm0, %xmm2

	movsd	%xmm1,  (XX)
	addq	INCX, XX
	movsd	%xmm2,  (XX)
	addq	INCX, XX
	ALIGN_3

L(154):
	testq	$1, M
	je	L(999)

	movsd	 (X), %xmm1
	mulsd	%xmm0, %xmm1
	movsd	%xmm1,  (X)
	ALIGN_4

L(999):
	xorq	%rax, %rax

	RESTOREREGISTERS

	ret

	EPILOGUE
