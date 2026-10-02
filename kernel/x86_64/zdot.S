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

#define N	ARG1	/* rdi */
#define X	ARG2	/* rsi */
#define INCX	ARG3	/* rdx */
#define Y	ARG4	/* rcx */
#ifndef WINDOWS_ABI
#define INCY	ARG5	/* r8  */
#else
#define INCY	%r10
#endif

#include "l1param.h"

	PROLOGUE
	PROFCODE

#ifdef WINDOWS_ABI
	emms

	movq	40(%rsp), INCY
#endif

	testq	N, N
	jle	L(88)

	salq	$ZBASE_SHIFT, INCX
	salq	$ZBASE_SHIFT, INCY

	fldz
	fldz
	fldz
	fldz

	cmpq	$2 * SIZE, INCX
	jne	L(14)
	cmpq	$2 * SIZE, INCY
	jne	L(14)

	movq	N,  %rax
	sarq	$1, %rax
	jle	L(15)
	ALIGN_3

L(16):
#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	FLD	 0 * SIZE(X)

	FLD	 0 * SIZE(Y)
	fmul	%st(1)
	faddp	%st, %st(2)

	FLD	 1 * SIZE(Y)
	fmulp	 %st, %st(1)
	faddp	%st, %st(2)
	FLD	 1 * SIZE(X)

	FLD	 0 * SIZE(Y)
	fmul	%st(1)
	faddp	%st, %st(4)

	FLD	 1 * SIZE(Y)
	fmulp	 %st, %st(1)
	faddp	%st, %st(4)
	FLD	 2 * SIZE(X)

	FLD	 2 * SIZE(Y)
	fmul	%st(1)
	faddp	%st, %st(2)

#ifdef PREFETCH
	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	FLD	 3 * SIZE(Y)
	fmulp	 %st, %st(1)
	faddp	%st, %st(2)
	FLD	 3 * SIZE(X)

	FLD	 2 * SIZE(Y)
	fmul	%st(1)
	faddp	%st, %st(4)

	FLD	 3 * SIZE(Y)
	fmulp	 %st, %st(1)
	faddp	%st, %st(4)

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	decq	%rax
	jg	L(16)
	ALIGN_3

L(15):
	movq	N, %rax
	andq	$1,   %rax
	jle	L(27)
	ALIGN_3

L(22):
	FLD	 0 * SIZE(X)

	FLD	 0 * SIZE(Y)
	fmul	%st(1)
	faddp	%st, %st(2)

	FLD	 1 * SIZE(Y)
	fmulp	 %st, %st(1)
	faddp	%st, %st(2)
	FLD	 1 * SIZE(X)

	FLD	 0 * SIZE(Y)
	fmul	%st(1)
	faddp	%st, %st(4)

	FLD	 1 * SIZE(Y)
	fmulp	 %st, %st(1)
	faddp	%st, %st(4)
	jmp	L(27)
	ALIGN_3

L(14):
	movq	N, %rax
	sarq	$1,   %rax
	jle	L(30)
	ALIGN_3


L(31):
	FLD	 0 * SIZE(X)

	FLD	 0 * SIZE(Y)
	fmul	%st(1)
	faddp	%st, %st(2)

	FLD	 1 * SIZE(Y)
	fmulp	 %st, %st(1)
	faddp	%st, %st(2)
	FLD	 1 * SIZE(X)

	FLD	 0 * SIZE(Y)
	fmul	%st(1)
	faddp	%st, %st(4)

	FLD	 1 * SIZE(Y)
	fmulp	 %st, %st(1)
	faddp	%st, %st(4)
	addq	INCX, X

	FLD	 0 * SIZE(X)
	addq	INCY, Y

	FLD	 0 * SIZE(Y)
	fmul	%st(1)
	faddp	%st, %st(2)

	FLD	 1 * SIZE(Y)
	fmulp	 %st, %st(1)
	faddp	%st, %st(2)
	FLD	 1 * SIZE(X)

	FLD	 0 * SIZE(Y)
	fmul	%st(1)
	faddp	%st, %st(4)

	FLD	 1 * SIZE(Y)
	fmulp	 %st, %st(1)
	faddp	%st, %st(4)
	addq	INCX, X
	addq	INCY, Y

	decq	%rax
	jg	L(31)
	ALIGN_3

L(30):
	movq	N,  %rax
	andq	$1, %rax
	jle	L(27)
	ALIGN_3

L(37):
	FLD	 0 * SIZE(X)

	FLD	 0 * SIZE(Y)
	fmul	%st(1)
	faddp	%st, %st(2)

	FLD	 1 * SIZE(Y)
	fmulp	 %st, %st(1)
	faddp	%st, %st(2)
	FLD	 1 * SIZE(X)

	FLD	 0 * SIZE(Y)
	fmul	%st(1)
	faddp	%st, %st(4)

	FLD	 1 * SIZE(Y)
	fmulp	 %st, %st(1)
	faddp	%st, %st(4)
	ALIGN_3

L(27):
#ifndef CONJ
	fsubp	%st, %st(3)
	faddp	%st, %st(1)
	fxch	%st(1)
#else
	faddp	%st, %st(3)
	fsubp	%st, %st(1)
	fxch	%st(1)
#endif
	ret
	ALIGN_3

L(88):
	fldz
	fldz

	ret
	EPILOGUE
