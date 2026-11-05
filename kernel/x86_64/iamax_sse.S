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

#define RET	%rax
#define I	ARG4
#define XX	%r10
#define	MM	%r11

#define MAXPS	maxps
#define MAXSS	maxss
#ifdef USE_MIN
#undef MAXPS
#undef MAXSS
#define MAXPS	minps
#define MAXSS	minss
#endif

.macro LOAD_AND_COMPARE_TO_MXX REG
	movss	0 * SIZE(X), \REG
	addq	INCX, X
#ifdef USE_ABS
	andps	 %xmm15, \REG
#endif
	cmpeqss	%xmm0, \REG
.endm

#include "l1param.h"

	PROLOGUE
	PROFCODE

	SAVEREGISTERS

	pxor	%xmm0, %xmm0		/* Return Value(Float) */
	xor	RET, RET		/* Return Value(Int)   */
	testq	M, M
	jle	L(999)
	leaq	(, INCX, SIZE), INCX
	testq	INCX, INCX
	jle	L(999)

	movq	M, MM
	movq	X, XX

#ifdef USE_ABS
	pcmpeqb	%xmm15, %xmm15
	psrld	$1, %xmm15		/* Generate USE_ABS */
#endif

	movss	(X), %xmm0
	addq	INCX, X
	decq	M
	shufps	$0, %xmm0, %xmm0
#ifdef USE_ABS
	andps	%xmm15, %xmm0
#endif
	movups	%xmm0, %xmm1
	movups	%xmm0, %xmm2
	movups	%xmm0, %xmm3		/* Generating "seed value" */
	cmpq	$SIZE, INCX
	jne	L(80)			/* Incx != 1 goto L80 */

/* Analigned Check */
	testq	$3, X			/* 00000011 */
	jne	L(30)			/* Purely Unaligned Mode */

	cmpq	$8, M
	jle	L(30)			/* if M <= 8 goto Unaligned mode */

	testq	$4, X			/* bit test 000100 */
	je	L(05)

	movss	0 * SIZE(X), %xmm4
#ifdef USE_ABS
	andps	%xmm15, %xmm4
#endif
	MAXSS	  %xmm4, %xmm0
	decq	M
	addq	$SIZE, X
	ALIGN_3

L(05):
	testq	$8, X
	je	L(06)

	movsd	0 * SIZE(X), %xmm4
	unpcklps  %xmm4, %xmm4
#ifdef USE_ABS
	andps	%xmm15, %xmm4
#endif
	MAXPS	  %xmm4, %xmm1
	subq	$2, M
	addq	$2 * SIZE, X
	ALIGN_3

L(06):
	movq	M,  I
	sarq	$4, I
	jle	L(15)
	ALIGN_4

L(11):
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movups	 0 * SIZE(X), %xmm4
#ifdef USE_ABS
	andps	%xmm15, %xmm4
#endif
	MAXPS	%xmm4, %xmm0

	movups	 4 * SIZE(X), %xmm5
#ifdef USE_ABS
	andps	%xmm15, %xmm5
#endif
	MAXPS	%xmm5, %xmm1

	movups	 8 * SIZE(X), %xmm6
#ifdef USE_ABS
	andps	%xmm15, %xmm6
#endif
	MAXPS	%xmm6, %xmm2

	movups	12 * SIZE(X), %xmm7
#ifdef USE_ABS
	andps	%xmm15, %xmm7
#endif
	MAXPS	%xmm7, %xmm3

	addq	$16 * SIZE, X
	decq	I
	jg	L(11)
	ALIGN_4

L(15):
	andq	$15,  M
	jle	L(20)

	testq	$8, M
	je	L(16)

	movups	0 * SIZE(X), %xmm4
#ifdef USE_ABS
	andps	%xmm15, %xmm4
#endif
	MAXPS	%xmm4, %xmm0

	movups	4 * SIZE(X), %xmm5
#ifdef USE_ABS
	andps	%xmm15, %xmm5
#endif
	MAXPS	%xmm5, %xmm1
	addq	$8 * SIZE, X
	ALIGN_3

L(16):
	testq	$4, M
	je	L(17)

	movups	0 * SIZE(X), %xmm6
#ifdef USE_ABS
	andps	%xmm15, %xmm6
#endif
	MAXPS	%xmm6, %xmm2
	addq	$4 * SIZE, X
	ALIGN_3

L(17):
	testq	$2, M
	je	L(18)

	movsd	0 * SIZE(X), %xmm7
	unpcklps %xmm7, %xmm7
#ifdef USE_ABS
	andps	%xmm15, %xmm7
#endif
	MAXPS	%xmm7, %xmm3
	addq	$2 * SIZE, X

L(18):
	testq	$1, M
	je	L(20)

	movss	0 * SIZE(X), %xmm4
#ifdef USE_ABS
	andps	%xmm15, %xmm4
#endif
	MAXSS	%xmm4, %xmm0
	ALIGN_3

L(20):
	movq	XX, X
	movq	MM, M

	MAXPS	%xmm1, %xmm0
	MAXPS	%xmm3, %xmm2
	MAXPS	%xmm2, %xmm0
	movups	%xmm0, %xmm1
	movhlps %xmm0, %xmm0
	MAXPS	%xmm1, %xmm0
	movups	%xmm0, %xmm1
	shufps	$1, %xmm0, %xmm0
	MAXSS	%xmm1, %xmm0
	shufps	$0, %xmm0, %xmm0

	testq	$4, X
	je	L(21)

	movss	0 * SIZE(X), %xmm1

	decq	M
	addq	$SIZE, X

#ifdef USE_ABS
	andps	 %xmm15, %xmm1
#endif
	incq	RET
	comiss	%xmm0, %xmm1
	je	L(999)
	ALIGN_3

L(21):
	testq	$8, X
	je	L(22)

	movss	0 * SIZE(X), %xmm1
	movss	1 * SIZE(X), %xmm2

	subq	$2, M
	addq	$2 * SIZE, X

#ifdef USE_ABS
	andps	 %xmm15, %xmm1
	andps	 %xmm15, %xmm2
#endif
	incq	RET
	comiss	%xmm0, %xmm1
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm2
	je	L(999)
	ALIGN_3

L(22):
	movq	M,  I
	sarq	$3, I
	jle	L(25)
	ALIGN_4

L(23):
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movups	0 * SIZE(X), %xmm1
#ifdef USE_ABS
	andps	 %xmm15, %xmm1
#endif
	cmpeqps	%xmm0, %xmm1

	movups	4 * SIZE(X), %xmm3
#ifdef USE_ABS
	andps	 %xmm15, %xmm3
#endif
	cmpeqps	%xmm0, %xmm3

	orps	%xmm3, %xmm1
#ifndef C_SUN
	movmskps %xmm1, %r11
#else
	.long	0xd9500f4c
#endif
	testq	 $15, %r11
	jne	 L(24)

	addq	$8 * SIZE, X
	addq	$8, RET
	decq	I
	jg	L(23)
	jmp	L(25)
	ALIGN_3

L(24):
	movss	0 * SIZE(X), %xmm1
	movss	1 * SIZE(X), %xmm2
	movss	2 * SIZE(X), %xmm3
	movss	3 * SIZE(X), %xmm4
	movss	4 * SIZE(X), %xmm5
	movss	5 * SIZE(X), %xmm6
	movss	6 * SIZE(X), %xmm7
	movss	7 * SIZE(X), %xmm8
#ifdef USE_ABS
	andps	 %xmm15, %xmm1
	andps	 %xmm15, %xmm2
	andps	 %xmm15, %xmm3
	andps	 %xmm15, %xmm4
	andps	 %xmm15, %xmm5
	andps	 %xmm15, %xmm6
	andps	 %xmm15, %xmm7
	andps	 %xmm15, %xmm8
#endif

	addq	$8 * SIZE, X

	incq	RET
	comiss	%xmm0, %xmm1
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm2
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm3
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm4
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm5
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm6
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm7
	je	L(999)
	incq	RET
	jmp	L(999)
	ALIGN_4

L(25):
	testq	$4, M
	je	L(26)

	movss	0 * SIZE(X), %xmm1
	movss	1 * SIZE(X), %xmm2
	movss	2 * SIZE(X), %xmm3
	movss	3 * SIZE(X), %xmm4
#ifdef USE_ABS
	andps	 %xmm15, %xmm1
	andps	 %xmm15, %xmm2
	andps	 %xmm15, %xmm3
	andps	 %xmm15, %xmm4
#endif
	addq	$4 * SIZE, X
	incq	RET
	comiss	%xmm0, %xmm1
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm2
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm3
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm4
	je	L(999)
	ALIGN_3

L(26):
	testq	$2, M
	je	L(27)

	movss	0 * SIZE(X), %xmm1
	movss	1 * SIZE(X), %xmm2
#ifdef USE_ABS
	andps	 %xmm15, %xmm1
	andps	 %xmm15, %xmm2
#endif
	addq	$2 * SIZE, X
	incq	RET
	comiss	%xmm0, %xmm1
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm2
	je	L(999)
	ALIGN_3

L(27):
	incq	RET
	jmp	L(999)
	ALIGN_3

/* Unaligned Mode */
L(30):
	movq	M,  I
	sarq	$4, I
	jle	L(35)
	ALIGN_4

L(31):
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movsd	 0 * SIZE(X), %xmm4
	movhps	 2 * SIZE(X), %xmm4
#ifdef USE_ABS
	andps	%xmm15, %xmm4
#endif
	MAXPS	%xmm4, %xmm0

	movsd	 4 * SIZE(X), %xmm5
	movhps	 6 * SIZE(X), %xmm5
#ifdef USE_ABS
	andps	%xmm15, %xmm5
#endif
	MAXPS	%xmm5, %xmm1

	movsd	 8 * SIZE(X), %xmm6
	movhps	10 * SIZE(X), %xmm6
#ifdef USE_ABS
	andps	%xmm15, %xmm6
#endif
	MAXPS	%xmm6, %xmm2

	movsd	12 * SIZE(X), %xmm7
	movhps	14 * SIZE(X), %xmm7
#ifdef USE_ABS
	andps	%xmm15, %xmm7
#endif
	MAXPS	%xmm7, %xmm3

	addq	$16 * SIZE, X
	decq	I
	jg	L(31)
	ALIGN_4

L(35):
	andq	$15,  M
	jle	L(40)

	testq	$8, M
	je	L(36)

	movsd	0 * SIZE(X), %xmm4
	movhps	2 * SIZE(X), %xmm4
#ifdef USE_ABS
	andps	%xmm15, %xmm4
#endif
	MAXPS	%xmm4, %xmm0

	movsd	4 * SIZE(X), %xmm5
	movhps	6 * SIZE(X), %xmm5
#ifdef USE_ABS
	andps	%xmm15, %xmm5
#endif
	MAXPS	%xmm5, %xmm1

	addq	$8 * SIZE, X
	ALIGN_3

L(36):
	testq	$4, M
	je	L(37)

	movsd	0 * SIZE(X), %xmm6
	movhps	2 * SIZE(X), %xmm6
#ifdef USE_ABS
	andps	%xmm15, %xmm6
#endif
	MAXPS	%xmm6, %xmm2
	addq	$4 * SIZE, X
	ALIGN_3

L(37):
	testq	$2, M
	je	L(38)

	movsd	0 * SIZE(X), %xmm7
	unpcklps %xmm7, %xmm7
#ifdef USE_ABS
	andps	%xmm15, %xmm7
#endif
	MAXPS	%xmm7, %xmm3
	addq	$2 * SIZE, X

L(38):
	testq	$1, M
	je	L(40)

	movss	0 * SIZE(X), %xmm4
#ifdef USE_ABS
	andps	%xmm15, %xmm4
#endif
	MAXSS	%xmm4, %xmm0
	jmp	L(40)
	ALIGN_4

L(40):
	movq	XX, X
	movq	MM, M

	MAXPS	%xmm1, %xmm0
	MAXPS	%xmm3, %xmm2
	MAXPS	%xmm2, %xmm0
	movups	%xmm0, %xmm1
	movhlps %xmm0, %xmm0
	MAXPS	%xmm1, %xmm0
	movups	%xmm0, %xmm1
	shufps	$1, %xmm0, %xmm0
	MAXSS	%xmm1, %xmm0
	shufps	$0, %xmm0, %xmm0

	movq	M,  I
	sarq	$3, I
	jle	L(45)
	ALIGN_4

L(43):
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movsd	0 * SIZE(X), %xmm1
	movhps	2 * SIZE(X), %xmm1
#ifdef USE_ABS
	andps	 %xmm15, %xmm1
#endif
	cmpeqps	%xmm0, %xmm1

	movsd	4 * SIZE(X), %xmm3
	movhps	6 * SIZE(X), %xmm3
#ifdef USE_ABS
	andps	 %xmm15, %xmm3
#endif
	cmpeqps	%xmm0, %xmm3

	orps	%xmm3, %xmm1
#ifndef C_SUN
	movmskps %xmm1, %r11
#else
	.long	0xd9500f4c
#endif
	testq	 $15, %r11
	jne	 L(44)

	addq	$8 * SIZE, X
	addq	$8, RET
	decq	I
	jg	L(43)
	jmp	L(45)
	ALIGN_3

L(44):
	movss	0 * SIZE(X), %xmm1
	movss	1 * SIZE(X), %xmm2
	movss	2 * SIZE(X), %xmm3
	movss	3 * SIZE(X), %xmm4
	movss	4 * SIZE(X), %xmm5
	movss	5 * SIZE(X), %xmm6
	movss	6 * SIZE(X), %xmm7
	movss	7 * SIZE(X), %xmm8
#ifdef USE_ABS
	andps	 %xmm15, %xmm1
	andps	 %xmm15, %xmm2
	andps	 %xmm15, %xmm3
	andps	 %xmm15, %xmm4
	andps	 %xmm15, %xmm5
	andps	 %xmm15, %xmm6
	andps	 %xmm15, %xmm7
	andps	 %xmm15, %xmm8
#endif

	addq	$8 * SIZE, X

	incq	RET
	comiss	%xmm0, %xmm1
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm2
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm3
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm4
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm5
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm6
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm7
	je	L(999)
	incq	RET
	jmp	L(999)
	ALIGN_4

L(45):
	testq	$4, M
	je	L(46)

	movss	0 * SIZE(X), %xmm1
	movss	1 * SIZE(X), %xmm2
	movss	2 * SIZE(X), %xmm3
	movss	3 * SIZE(X), %xmm4
#ifdef USE_ABS
	andps	 %xmm15, %xmm1
	andps	 %xmm15, %xmm2
	andps	 %xmm15, %xmm3
	andps	 %xmm15, %xmm4
#endif
	addq	$4 * SIZE, X
	incq	RET
	comiss	%xmm0, %xmm1
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm2
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm3
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm4
	je	L(999)
	ALIGN_3

L(46):
	testq	$2, M
	je	L(47)

	movss	0 * SIZE(X), %xmm1
	movss	1 * SIZE(X), %xmm2
#ifdef USE_ABS
	andps	 %xmm15, %xmm1
	andps	 %xmm15, %xmm2
#endif
	addq	$2 * SIZE, X
	incq	RET
	comiss	%xmm0, %xmm1
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm2
	je	L(999)
	ALIGN_3

L(47):
	incq	RET
	jmp	L(999)
	ALIGN_3

L(80):
	movq	M,  I
	sarq	$3, I
	jle	L(85)
	ALIGN_4

L(81):
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	movss	0 * SIZE(X), %xmm4
	addq	INCX, X
#ifdef USE_ABS
	andps	%xmm15, %xmm4
#endif
	MAXSS	%xmm4, %xmm0

	movss	0 * SIZE(X), %xmm5
	addq	INCX, X
#ifdef USE_ABS
	andps	%xmm15, %xmm5
#endif
	MAXSS	%xmm5, %xmm1

	movss	0 * SIZE(X), %xmm6
	addq	INCX, X
#ifdef USE_ABS
	andps	%xmm15, %xmm6
#endif
	MAXSS	%xmm6, %xmm2

	movss	0 * SIZE(X), %xmm7
	addq	INCX, X
#ifdef USE_ABS
	andps	%xmm15, %xmm7
#endif
	MAXSS	%xmm7, %xmm3

	movss	0 * SIZE(X), %xmm4
	addq	INCX, X
#ifdef USE_ABS
	andps	%xmm15, %xmm4
#endif
	MAXSS	%xmm4, %xmm0

	movss	0 * SIZE(X), %xmm5
	addq	INCX, X
#ifdef USE_ABS
	andps	%xmm15, %xmm5
#endif
	MAXSS	%xmm5, %xmm1

	movss	0 * SIZE(X), %xmm6
	addq	INCX, X
#ifdef USE_ABS
	andps	%xmm15, %xmm6
#endif
	MAXSS	%xmm6, %xmm2

	movss	0 * SIZE(X), %xmm7
	addq	INCX, X
#ifdef USE_ABS
	andps	%xmm15, %xmm7
#endif
	MAXSS	%xmm7, %xmm3

	decq	I
	jg	L(81)
	ALIGN_4

L(85):
	andq	$7,  M
	jle	L(90)

	testq	$4, M
	je	L(86)

	movss	0 * SIZE(X), %xmm4
	addq	INCX, X
#ifdef USE_ABS
	andps	%xmm15, %xmm4
#endif
	MAXSS	%xmm4, %xmm0

	movss	0 * SIZE(X), %xmm5
	addq	INCX, X
#ifdef USE_ABS
	andps	%xmm15, %xmm5
#endif
	MAXSS	%xmm5, %xmm1

	movss	0 * SIZE(X), %xmm6
	addq	INCX, X
#ifdef USE_ABS
	andps	%xmm15, %xmm6
#endif
	MAXSS	%xmm6, %xmm2

	movss	0 * SIZE(X), %xmm7
	addq	INCX, X
#ifdef USE_ABS
	andps	%xmm15, %xmm7
#endif
	MAXSS	%xmm7, %xmm3
	ALIGN_3

L(86):
	testq	$2, M
	je	L(87)

	movss	0 * SIZE(X), %xmm4
	addq	INCX, X
#ifdef USE_ABS
	andps	%xmm15, %xmm4
#endif
	MAXSS	%xmm4, %xmm0

	movss	0 * SIZE(X), %xmm5
	addq	INCX, X
#ifdef USE_ABS
	andps	%xmm15, %xmm5
#endif
	MAXSS	%xmm5, %xmm1
	ALIGN_3

L(87):
	testq	$1, M
	je	L(90)

	movss	0 * SIZE(X), %xmm6
	addq	INCX, X
#ifdef USE_ABS
	andps	%xmm15, %xmm6
#endif
	MAXSS	%xmm6, %xmm2
	ALIGN_4

L(90):
	movq	XX, X
	movq	MM, M

	MAXSS	%xmm1, %xmm0
	MAXSS	%xmm3, %xmm2
	MAXSS	%xmm2, %xmm0
	shufps	$0, %xmm0, %xmm0

	movq	M,  I
	sarq	$3, I
	jle	L(95)
	ALIGN_4

L(93):
  LOAD_AND_COMPARE_TO_MXX %xmm1
  LOAD_AND_COMPARE_TO_MXX %xmm2
  LOAD_AND_COMPARE_TO_MXX %xmm3
  LOAD_AND_COMPARE_TO_MXX %xmm4
  LOAD_AND_COMPARE_TO_MXX %xmm5
  LOAD_AND_COMPARE_TO_MXX %xmm6
  LOAD_AND_COMPARE_TO_MXX %xmm7
  LOAD_AND_COMPARE_TO_MXX %xmm8

	orps	%xmm2, %xmm1
	orps	%xmm4, %xmm3
	orps	%xmm6, %xmm5
	orps	%xmm8, %xmm7
	orps	%xmm3, %xmm1
	orps	%xmm7, %xmm5
	orps	%xmm5, %xmm1

#ifndef C_SUN
	movmskps %xmm1, %r11
#else
	.long	0xd9500f4c
#endif
	testq	 $15, %r11
	jne	 L(94)

	addq	$8, RET
	decq	I
	jg	L(93)
	jmp	L(95)
	ALIGN_3

L(94):
	subq	INCX, X
	movss	0 * SIZE(X), %xmm8
	subq	INCX, X
	movss	0 * SIZE(X), %xmm7
	subq	INCX, X
	movss	0 * SIZE(X), %xmm6
	subq	INCX, X
	movss	0 * SIZE(X), %xmm5
	subq	INCX, X
	movss	0 * SIZE(X), %xmm4
	subq	INCX, X
	movss	0 * SIZE(X), %xmm3
	subq	INCX, X
	movss	0 * SIZE(X), %xmm2
	subq	INCX, X
	movss	0 * SIZE(X), %xmm1
#ifdef USE_ABS
	andps	 %xmm15, %xmm1
	andps	 %xmm15, %xmm2
	andps	 %xmm15, %xmm3
	andps	 %xmm15, %xmm4
	andps	 %xmm15, %xmm5
	andps	 %xmm15, %xmm6
	andps	 %xmm15, %xmm7
	andps	 %xmm15, %xmm8
#endif
	incq	RET
	comiss	%xmm0, %xmm1
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm2
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm3
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm4
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm5
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm6
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm7
	je	L(999)
	incq	RET
	jmp	L(999)
	ALIGN_4

L(95):
	testq	$4, M
	je	L(96)

	movss	0 * SIZE(X), %xmm1
	addq	INCX, X
	movss	0 * SIZE(X), %xmm2
	addq	INCX, X
	movss	0 * SIZE(X), %xmm3
	addq	INCX, X
	movss	0 * SIZE(X), %xmm4
	addq	INCX, X

#ifdef USE_ABS
	andps	 %xmm15, %xmm1
	andps	 %xmm15, %xmm2
	andps	 %xmm15, %xmm3
	andps	 %xmm15, %xmm4
#endif
	incq	RET
	comiss	%xmm0, %xmm1
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm2
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm3
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm4
	je	L(999)
	ALIGN_3

L(96):
	testq	$2, M
	je	L(97)

	movss	0 * SIZE(X), %xmm1
	addq	INCX, X
	movss	0 * SIZE(X), %xmm2
	addq	INCX, X
#ifdef USE_ABS
	andps	 %xmm15, %xmm1
	andps	 %xmm15, %xmm2
#endif
	incq	RET
	comiss	%xmm0, %xmm1
	je	L(999)
	incq	RET
	comiss	%xmm0, %xmm2
	je	L(999)
	ALIGN_3

L(97):
	incq	RET
	ALIGN_3

L(999):
	RESTOREREGISTERS

	ret

	EPILOGUE

