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
#define M	ARG1	/* rdi */
#define X	ARG2	/* rsi */
#define INCX	ARG3	/* rdx */

#define I	%rax

#include "l1param.h"

	PROLOGUE
	PROFCODE

	SAVEREGISTERS

	pxor	%xmm0, %xmm0
	testq	M, M
	jle	.L999
	pxor	%xmm1, %xmm1
	testq	INCX, INCX
	je	.L999

	pxor	%xmm2, %xmm2
	leaq	(, INCX, SIZE), INCX
	pxor	%xmm3, %xmm3
	cmpq	$SIZE, INCX
	jne	.L40

	testq	$SIZE, X
	je	.L05

#ifndef DOUBLE
	movss	0 * SIZE(X), %xmm4
	cvtss2sd %xmm4, %xmm6
#else
	movsd   0 * SIZE(X), %xmm6
#endif
	mulsd	 %xmm6, %xmm6
	addsd	 %xmm6, %xmm3
//#ifndef DOUBLE
	addq	INCX, X
//#else
//	addq	$1 *SIZE, X
//#endif
	decq	M
	jle	.L998
	ALIGN_3

.L05:
	movq	M,  I
	sarq	$3, I
	jle	.L14
#ifndef DOUBLE
	movsd	0 * SIZE(X), %xmm4
	movsd	2 * SIZE(X), %xmm5
	movsd	4 * SIZE(X), %xmm6
	movsd	6 * SIZE(X), %xmm7
	addq	$8 * SIZE, X
#else
	movsd	0 * SIZE(X), %xmm4
	movsd	1 * SIZE(X), %xmm5
	movsd	2 * SIZE(X), %xmm6
	movsd	3 * SIZE(X), %xmm7
	movsd	4 * SIZE(X), %xmm8
	movsd	5 * SIZE(X), %xmm9
	movsd	6 * SIZE(X), %xmm10
	movsd	7 * SIZE(X), %xmm11
	addq	$8 * SIZE, X
#endif
	decq	I
	jle	.L12
	ALIGN_3

.L10:
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif
#ifndef DOUBLE
	cvtps2pd %xmm4, %xmm8
	cvtps2pd %xmm5, %xmm9
	cvtps2pd %xmm6, %xmm10
	cvtps2pd %xmm7, %xmm11
#else
	mulpd %xmm4, %xmm4
	mulpd %xmm5, %xmm5
	mulpd %xmm6, %xmm6
	mulpd %xmm7, %xmm7
	addpd %xmm4, %xmm0
	addpd %xmm5, %xmm1
	addpd %xmm6, %xmm2
	addpd %xmm7, %xmm3
#endif
#ifndef DOUBLE
	movsd	0 * SIZE(X), %xmm4
	movsd	2 * SIZE(X), %xmm5
	movsd	4 * SIZE(X), %xmm6
	movsd	6 * SIZE(X), %xmm7
#endif
	mulpd	 %xmm8,  %xmm8
	mulpd	 %xmm9,  %xmm9
	mulpd	 %xmm10, %xmm10
	mulpd	 %xmm11, %xmm11

	addpd	 %xmm8,  %xmm0
	addpd	 %xmm9,  %xmm1
	addpd	 %xmm10, %xmm2
	addpd	 %xmm11, %xmm3

#ifdef DOUBLE
	movsd	0 * SIZE(X), %xmm4
	movsd	1 * SIZE(X), %xmm5
	movsd	2 * SIZE(X), %xmm6
	movsd	3 * SIZE(X), %xmm7
	movsd	4 * SIZE(X), %xmm8
	movsd	5 * SIZE(X), %xmm9
	movsd	6 * SIZE(X), %xmm10
	movsd	7 * SIZE(X), %xmm11
#endif

	addq	$8 * SIZE, X
	decq	I
	jg	.L10
	ALIGN_3

.L12:
#ifndef DOUBLE
	cvtps2pd %xmm4, %xmm8
	cvtps2pd %xmm5, %xmm9
	cvtps2pd %xmm6, %xmm10
	cvtps2pd %xmm7, %xmm11
#else
	mulpd %xmm4, %xmm4
	mulpd %xmm5, %xmm5
	mulpd %xmm6, %xmm6
	mulpd %xmm7, %xmm7
	addpd %xmm4, %xmm0
	addpd %xmm5, %xmm1
	addpd %xmm6, %xmm2
	addpd %xmm7, %xmm3
#endif
	mulpd	 %xmm8,  %xmm8
	mulpd	 %xmm9,  %xmm9
	mulpd	 %xmm10, %xmm10
	mulpd	 %xmm11, %xmm11

	addpd	 %xmm8,  %xmm0
	addpd	 %xmm9,  %xmm1
	addpd	 %xmm10, %xmm2
	addpd	 %xmm11, %xmm3
	ALIGN_3


.L14:
	testq	$4,  M
	je	.L15
#ifndef DOUBLE
	movsd	0 * SIZE(X), %xmm4
	movsd	2 * SIZE(X), %xmm5
	cvtps2pd %xmm4, %xmm6
	cvtps2pd %xmm5, %xmm7
#else
	movupd	0 * SIZE(X), %xmm4
	movupd	1 * SIZE(X), %xmm5
	movupd	2 * SIZE(X), %xmm6
	movupd	3 * SIZE(X), %xmm7
	mulpd	 %xmm4, %xmm4
	mulpd	 %xmm5, %xmm5
	addpd	 %xmm4, %xmm0
	addpd	 %xmm5, %xmm1
#endif
	mulpd	 %xmm6, %xmm6
	mulpd	 %xmm7, %xmm7
	addpd	 %xmm6, %xmm0
	addpd	 %xmm7, %xmm1
#ifndef DOUBLW
	addq	$4 * SIZE, X
#else
	addq	$4 * SIZE, X
#endif
	ALIGN_3

.L15:
	testq	$2,  M
	je	.L16

	movsd	0 * SIZE(X), %xmm4
#ifndef DOUBLE
	cvtps2pd %xmm4, %xmm6
#else
	movsd	1 * SIZE(X), %xmm6
	mulpd	 %xmm4, %xmm4
	addpd	 %xmm4, %xmm2
#endif
	mulpd	 %xmm6, %xmm6
	addpd	 %xmm6, %xmm2
#ifndef DOUBLE
	addq	$2 * SIZE, X
#else
	addq	$2 * SIZE, X
#endif
	ALIGN_3

.L16:
	testq	$1,  M
	je	.L998

#ifndef DOUBLE
	movss	0 * SIZE(X), %xmm4
	cvtss2sd %xmm4, %xmm6
#else
	movsd   0 * SIZE(X), %xmm6
#endif
	mulsd	 %xmm6, %xmm6
	addsd	 %xmm6, %xmm3
	jmp	.L998
	ALIGN_4

.L40:
	movq	M,  I
	sarq	$3, I
	jle	.L44
	ALIGN_4

.L41:
#ifndef DOUBLE
	movss	(X), %xmm4
	addq	INCX, X
	movss	(X), %xmm5
	addq	INCX, X
	movss	(X), %xmm6
	addq	INCX, X
	movss	(X), %xmm7
	addq	INCX, X
	movss	(X), %xmm8
	addq	INCX, X
	movss	(X), %xmm9
	addq	INCX, X
	movss	(X), %xmm10
	addq	INCX, X
	movss	(X), %xmm11
	addq	INCX, X

	cvtss2sd %xmm4,  %xmm4
	cvtss2sd %xmm5,  %xmm5
	cvtss2sd %xmm6,  %xmm6
	cvtss2sd %xmm7,  %xmm7
	cvtss2sd %xmm8,  %xmm8
	cvtss2sd %xmm9,  %xmm9
	cvtss2sd %xmm10, %xmm10
	cvtss2sd %xmm11, %xmm11
#else
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
	movsd	(X), %xmm9
	addq	INCX, X
	movsd	(X), %xmm10
	addq	INCX, X
	movsd	(X), %xmm11
	addq	INCX, X
#endif
	mulsd	 %xmm4,  %xmm4
	mulsd	 %xmm5,  %xmm5
	mulsd	 %xmm6,  %xmm6
	mulsd	 %xmm7,  %xmm7

	addsd	 %xmm4,  %xmm0
	addsd	 %xmm5,  %xmm1
	addsd	 %xmm6,  %xmm2
	addsd	 %xmm7,  %xmm3

	mulsd	 %xmm8,  %xmm8
	mulsd	 %xmm9,  %xmm9
	mulsd	 %xmm10, %xmm10
	mulsd	 %xmm11, %xmm11

	addsd	 %xmm8,  %xmm0
	addsd	 %xmm9,  %xmm1
	addsd	 %xmm10, %xmm2
	addsd	 %xmm11, %xmm3

	decq	I
	jg	.L41
	ALIGN_3

.L44:
	testq	$4,  M
	je	.L45
#ifndef DOUBLE
	movss	(X), %xmm4
	addq	INCX, X
	movss	(X), %xmm5
	addq	INCX, X
	movss	(X), %xmm6
	addq	INCX, X
	movss	(X), %xmm7
	addq	INCX, X

	cvtss2sd %xmm4, %xmm8
	cvtss2sd %xmm5, %xmm9
	cvtss2sd %xmm6, %xmm10
	cvtss2sd %xmm7, %xmm11
#else
	movsd	(X), %xmm8
	addq	INCX, X
	movsd	(X), %xmm9
	addq	INCX, X
	movsd	(X), %xmm10
	addq	INCX, X
	movsd	(X), %xmm11
	addq	INCX, X
#endif
	mulsd	 %xmm8,  %xmm8
	mulsd	 %xmm9,  %xmm9
	mulsd	 %xmm10, %xmm10
	mulsd	 %xmm11, %xmm11

	addsd	 %xmm8,  %xmm0
	addsd	 %xmm9,  %xmm1
	addsd	 %xmm10, %xmm2
	addsd	 %xmm11, %xmm3
	ALIGN_3

.L45:
	testq	$2,  M
	je	.L46
#ifndef DOUBLE
	movss	(X), %xmm4
	addq	INCX, X
	movss	(X), %xmm5
	addq	INCX, X

	cvtss2sd %xmm4, %xmm6
	cvtss2sd %xmm5, %xmm7
#else
	movsd	(X), %xmm6
	addq	INCX, X
	movsd	(X), %xmm7
	addq	INCX, X
#endif
	mulsd	 %xmm6, %xmm6
	mulsd	 %xmm7, %xmm7
	addsd	 %xmm6, %xmm1
	addsd	 %xmm7, %xmm2
	ALIGN_3

.L46:
	testq	$1,  M
	je	.L998
#ifndef DOUBLE
	movss	(X), %xmm4
	cvtss2sd %xmm4, %xmm6
#else
	movsd   (X), %xmm6
#endif
	mulsd	 %xmm6, %xmm6
	addsd	 %xmm6, %xmm3
	ALIGN_4

.L998:
	addpd	%xmm1, %xmm0
	addpd	%xmm3, %xmm2
	addpd	%xmm2, %xmm0

#ifndef HAVE_SSE3
	movapd	%xmm0, %xmm1
	unpckhpd	%xmm0, %xmm0
	addsd	%xmm1, %xmm0
#else
	haddpd	%xmm0, %xmm0
#endif

	ALIGN_4

.L999:
	sqrtsd	%xmm0,  %xmm0
#ifndef DOUBLE
	cvtsd2ss %xmm0, %xmm0
#endif
	RESTOREREGISTERS

	ret

	EPILOGUE

