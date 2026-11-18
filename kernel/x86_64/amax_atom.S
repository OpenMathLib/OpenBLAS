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

#ifdef USE_MIN
#define maxsd	minsd
#endif

#include "l1param.h"

	PROLOGUE
	PROFCODE

	SAVEREGISTERS

	xorps	%xmm0, %xmm0
	leaq	(, INCX, SIZE), INCX

	testq	M, M
	jle	L(999)

	testq	INCX, INCX
	jle	L(999)

#ifdef USE_ABS
	pcmpeqb	%xmm15, %xmm15
	psrlq	$1, %xmm15
#endif

	movsd	(X), %xmm0
	addq	INCX, X

#ifdef USE_ABS
	andps	 %xmm15, %xmm0
#endif
	decq	M
	jle	L(999)

	movaps	 %xmm0, %xmm1
	movaps	 %xmm0, %xmm2
	movaps	 %xmm0, %xmm3

	cmpq	$SIZE, INCX
	jne	L(20)

	movq	M,  I
	sarq	$3, I
	jle	L(15)

	movsd	0 * SIZE(X), %xmm4
	movsd	1 * SIZE(X), %xmm5
	movsd	2 * SIZE(X), %xmm6
	movsd	3 * SIZE(X), %xmm7

	movsd	4 * SIZE(X), %xmm8
	movsd	5 * SIZE(X), %xmm9
	movsd	6 * SIZE(X), %xmm10
	movsd	7 * SIZE(X), %xmm11

	decq	I
	jle	L(13)
	ALIGN_4

L(12):
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

#ifdef USE_ABS
	andps	%xmm15, %xmm4
#endif
	maxsd	%xmm4,  %xmm1
	movsd	 8 * SIZE(X), %xmm4

#ifdef USE_ABS
	andps	%xmm15, %xmm5
#endif
	maxsd	%xmm5,  %xmm2
	movsd	 9 * SIZE(X), %xmm5

#ifdef USE_ABS
	andps	%xmm15, %xmm6
#endif
	maxsd	%xmm6,  %xmm1
	movsd	10 * SIZE(X), %xmm6

#ifdef USE_ABS
	andps	%xmm15, %xmm7
#endif
	maxsd	%xmm7,  %xmm2
	movsd	11 * SIZE(X), %xmm7

#ifdef USE_ABS
	andps	%xmm15, %xmm8
#endif
	maxsd	%xmm8,  %xmm1
	movsd	12 * SIZE(X), %xmm8

#ifdef USE_ABS
	andps	%xmm15, %xmm9
#endif
	maxsd	%xmm9,  %xmm2
	movsd	13 * SIZE(X), %xmm9

#ifdef USE_ABS
	andps	%xmm15, %xmm10
#endif
	maxsd	%xmm10, %xmm1
	movsd	14 * SIZE(X), %xmm10

#ifdef USE_ABS
	andps	%xmm15, %xmm11
#endif
	maxsd	%xmm11, %xmm2
	movsd	15 * SIZE(X), %xmm11

	addq	$8 * SIZE, X
	decq	I
	jg	L(12)
	ALIGN_4

L(13):
#ifdef USE_ABS
	andps	%xmm15, %xmm4
#endif
	maxsd	%xmm4,  %xmm0
#ifdef USE_ABS
	andps	%xmm15, %xmm5
#endif
	maxsd	%xmm5,  %xmm1
#ifdef USE_ABS
	andps	%xmm15, %xmm6
#endif
	maxsd	%xmm6,  %xmm2
#ifdef USE_ABS
	andps	%xmm15, %xmm7
#endif
	maxsd	%xmm7,  %xmm3

#ifdef USE_ABS
	andps	%xmm15, %xmm8
#endif
	maxsd	%xmm8,  %xmm0
#ifdef USE_ABS
	andps	%xmm15, %xmm9
#endif
	maxsd	%xmm9,  %xmm1
#ifdef USE_ABS
	andps	%xmm15, %xmm10
#endif
	maxsd	%xmm10, %xmm2
#ifdef USE_ABS
	andps	%xmm15, %xmm11
#endif
	maxsd	%xmm11, %xmm3

	addq	$8 * SIZE, X
	ALIGN_4

L(15):
	testq	$4, M
	jle	L(17)

	movsd	0 * SIZE(X), %xmm4
	movsd	1 * SIZE(X), %xmm5
	movsd	2 * SIZE(X), %xmm6
	movsd	3 * SIZE(X), %xmm7

#ifdef USE_ABS
	andps	%xmm15, %xmm4
#endif
	maxsd	%xmm4,  %xmm0
#ifdef USE_ABS
	andps	%xmm15, %xmm5
#endif
	maxsd	%xmm5,  %xmm1
#ifdef USE_ABS
	andps	%xmm15, %xmm6
#endif
	maxsd	%xmm6,  %xmm2
#ifdef USE_ABS
	andps	%xmm15, %xmm7
#endif
	maxsd	%xmm7,  %xmm3

	addq	$4 * SIZE, X
	ALIGN_3

L(17):
	testq	$2, M
	jle	L(18)

	movsd	0 * SIZE(X), %xmm4
	movsd	1 * SIZE(X), %xmm5

#ifdef USE_ABS
	andps	 %xmm15, %xmm4
#endif
	maxsd	%xmm4, %xmm1
#ifdef USE_ABS
	andps	 %xmm15, %xmm5
#endif
	maxsd	%xmm5, %xmm2
	addq	$2 * SIZE, X
	ALIGN_3

L(18):
	testq	$1, M
	jle	L(998)

	movsd	0 * SIZE(X), %xmm4
#ifdef USE_ABS
	andps	 %xmm15, %xmm4
#endif
	maxsd	%xmm4, %xmm3
	jmp	L(998)
	ALIGN_3

L(20):
	movq	M,  I
	sarq	$3, I
	jle	L(25)

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

	decq	I
	jle	L(23)
	ALIGN_4

L(22):
#ifdef USE_ABS
	andps	%xmm15, %xmm4
#endif
	addq	INCX, X
	maxsd	%xmm4,  %xmm1
	movsd	(X), %xmm4

#ifdef USE_ABS
	andps	%xmm15, %xmm5
#endif
	addq	INCX, X
	maxsd	%xmm5,  %xmm2
	movsd	(X), %xmm5

#ifdef USE_ABS
	andps	%xmm15, %xmm6
#endif
	addq	INCX, X
	maxsd	%xmm6,  %xmm1
	movsd	(X), %xmm6

#ifdef USE_ABS
	andps	%xmm15, %xmm7
#endif
	addq	INCX, X
	maxsd	%xmm7,  %xmm2
	movsd	(X), %xmm7

#ifdef USE_ABS
	andps	%xmm15, %xmm8
#endif
	addq	INCX, X
	maxsd	%xmm8,  %xmm1
	movsd	(X), %xmm8

#ifdef USE_ABS
	andps	%xmm15, %xmm9
#endif
	addq	INCX, X
	maxsd	%xmm9,  %xmm2
	movsd	(X), %xmm9

#ifdef USE_ABS
	andps	%xmm15, %xmm10
#endif
	addq	INCX, X
	maxsd	%xmm10, %xmm1
	movsd	(X), %xmm10

#ifdef USE_ABS
	andps	%xmm15, %xmm11
#endif
	addq	INCX, X
	maxsd	%xmm11, %xmm2
	movsd	(X), %xmm11

	decq	I
	jg	L(22)
	ALIGN_4

L(23):
#ifdef USE_ABS
	andps	%xmm15, %xmm4
#endif
	addq	INCX, X
	maxsd	%xmm4,  %xmm0
#ifdef USE_ABS
	andps	%xmm15, %xmm5
#endif
	maxsd	%xmm5,  %xmm1
#ifdef USE_ABS
	andps	%xmm15, %xmm6
#endif
	maxsd	%xmm6,  %xmm2
#ifdef USE_ABS
	andps	%xmm15, %xmm7
#endif
	maxsd	%xmm7,  %xmm3

#ifdef USE_ABS
	andps	%xmm15, %xmm8
#endif
	maxsd	%xmm8,  %xmm0
#ifdef USE_ABS
	andps	%xmm15, %xmm9
#endif
	maxsd	%xmm9,  %xmm1
#ifdef USE_ABS
	andps	%xmm15, %xmm10
#endif
	maxsd	%xmm10, %xmm2
#ifdef USE_ABS
	andps	%xmm15, %xmm11
#endif
	maxsd	%xmm11, %xmm3
	ALIGN_4

L(25):
	testq	$4, M
	jle	L(27)

	movsd	(X), %xmm4
	addq	INCX, X
	movsd	(X), %xmm5
	addq	INCX, X
	movsd	(X), %xmm6
	addq	INCX, X
	movsd	(X), %xmm7
	addq	INCX, X

#ifdef USE_ABS
	andps	%xmm15, %xmm4
#endif
	maxsd	%xmm4,  %xmm0
#ifdef USE_ABS
	andps	%xmm15, %xmm5
#endif
	maxsd	%xmm5,  %xmm1
#ifdef USE_ABS
	andps	%xmm15, %xmm6
#endif
	maxsd	%xmm6,  %xmm2
#ifdef USE_ABS
	andps	%xmm15, %xmm7
#endif
	maxsd	%xmm7,  %xmm3
	ALIGN_3

L(27):
	testq	$2, M
	jle	L(28)

	movsd	(X), %xmm4
	addq	INCX, X
	movsd	(X), %xmm5
	addq	INCX, X

#ifdef USE_ABS
	andps	 %xmm15, %xmm4
#endif
	maxsd	%xmm4, %xmm1
#ifdef USE_ABS
	andps	 %xmm15, %xmm5
#endif
	maxsd	%xmm5, %xmm2
	ALIGN_3

L(28):
	testq	$1, M
	jle	L(998)

	movsd	(X), %xmm4
	addq	INCX, X
#ifdef USE_ABS
	andps	 %xmm15, %xmm4
#endif
	maxsd	%xmm4, %xmm3
	ALIGN_3

L(998):
	maxsd	%xmm1, %xmm0
	maxsd	%xmm3, %xmm2
	maxsd	%xmm2, %xmm0
	ALIGN_4

L(999):

	RESTOREREGISTERS

	ret

	EPILOGUE
