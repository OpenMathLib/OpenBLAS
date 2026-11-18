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

#define STACK	12
#define ARGS     0

#define STACK_N		 4 + STACK + ARGS(%esp)
#define STACK_X		 8 + STACK + ARGS(%esp)
#define STACK_INCX	12 + STACK + ARGS(%esp)
#define STACK_Y		16 + STACK + ARGS(%esp)
#define STACK_INCY	20 + STACK + ARGS(%esp)

#define N	%ebx
#define X	%esi
#define INCX	%ecx
#define Y	%edi
#define INCY	%edx

#include "l1param.h"

	PROLOGUE

#ifdef WINDOWS_ABI
	emms
#endif

	pushl	%edi
	pushl	%esi
	pushl	%ebx

	PROFCODE

	movl	STACK_N,     N
	movl	STACK_X,     X
	movl	STACK_INCX,  INCX
	movl	STACK_Y,     Y
	movl	STACK_INCY,  INCY

	sall	$BASE_SHIFT, INCX
	sall	$BASE_SHIFT, INCY

	fldz
	fldz
	fldz
	fldz

	cmpl	$SIZE, INCX
	jne	L(14)
	cmpl	$SIZE, INCY
	jne	L(14)

	movl	N, %eax
	sarl	$2,   %eax
	jle	L(15)
	ALIGN_3

L(16):
#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(X)
#endif

	FLD	0 * SIZE(X)
	FLD	0 * SIZE(Y)
	fmulp	%st, %st(1)
	faddp	%st,%st(1)
	FLD	1 * SIZE(X)
	FLD	1 * SIZE(Y)
	fmulp	%st, %st(1)
	faddp	%st,%st(2)

#ifdef PREFETCHW
	PREFETCHW (PREFETCHSIZE +  0) - PREOFFSET(Y)
#endif

	FLD	2 * SIZE(X)
	FLD	2 * SIZE(Y)
	fmulp	%st, %st(1)
	faddp	%st,%st(3)
	FLD	3 * SIZE(X)
	FLD	3 * SIZE(Y)
	fmulp	%st, %st(1)
	faddp	%st,%st(4)
	addl	$4 * SIZE, X
	addl	$4 * SIZE, Y
	decl	%eax
	jg	L(16)
	ALIGN_3

L(15):
	movl	N, %eax
	andl	$3,   %eax
	jle	L(27)
	ALIGN_3

L(22):
	FLD	(X)
	addl	$SIZE, X
	FLD	(Y)
	fmulp	%st, %st(1)
	addl	$SIZE, Y
	faddp	%st,%st(1)
	decl	%eax
	jg	L(22)

	jmp	L(27)
	ALIGN_3

L(14):
	movl	N, %eax
	sarl	$2,   %eax
	jle	L(30)
	ALIGN_3

L(31):
	FLD	(X)
	addl	INCX, X
	FLD	(Y)
	fmulp	%st, %st(1)
	addl	INCY, Y
	faddp	%st,%st(1)

	FLD	(X)
	addl	INCX, X
	FLD	(Y)
	fmulp	%st, %st(1)
	addl	INCY, Y
	faddp	%st,%st(2)

	FLD	(X)
	addl	INCX, X
	FLD	(Y)
	fmulp	%st, %st(1)
	addl	INCY, Y
	faddp	%st,%st(3)

	FLD	(X)
	addl	INCX, X
	FLD	(Y)
	fmulp	%st, %st(1)
	addl	INCY, Y
	faddp	%st,%st(4)

	decl	%eax
	jg	L(31)
	ALIGN_3

L(30):
	movl	N, %eax
	andl	$3,   %eax
	jle	L(27)
	ALIGN_3

L(37):
	FLD	(X)
	addl	INCX, X
	FLD	(Y)
	fmulp	%st, %st(1)
	addl	INCY, Y
	faddp	%st, %st(1)
	decl	%eax
	jg	L(37)
	ALIGN_3

L(27):
	faddp	%st,%st(2)
	faddp	%st,%st(2)
	faddp	%st,%st(1)

	popl	%ebx
	popl	%esi
	popl	%edi
	ret

	EPILOGUE
