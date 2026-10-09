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

#define A_PRE 512

#include "l1param.h"

	PROLOGUE
	PROFCODE

#ifdef WINDOWS_ABI
	movq	40(%rsp), INCY
#endif

	SAVEREGISTERS

	leaq	(, INCX, SIZE), INCX
	leaq	(, INCY, SIZE), INCY

	vxorps	%xmm0, %xmm0 , %xmm0
	vxorps	%xmm1, %xmm1 , %xmm1
	vxorps	%xmm2, %xmm2 , %xmm2
	vxorps	%xmm3, %xmm3 , %xmm3

	cmpq	$0, N
	jle	L(999)

	cmpq	$SIZE, INCX
	jne	L(50)
	cmpq	$SIZE, INCY
	jne	L(50)

	subq	$-16 * SIZE, X
	subq	$-16 * SIZE, Y

	testq	$SIZE, Y
	je	L(10)

	vmovsd	-16 * SIZE(X), %xmm0
	vmulsd	-16 * SIZE(Y), %xmm0 , %xmm0
	addq	$1 * SIZE, X
	addq	$1 * SIZE, Y
	decq	N
	ALIGN_2

L(10):

	movq	N,  %rax
	sarq	$4, %rax
	jle	L(14)

	vmovups	-16 * SIZE(X), %xmm4
	vmovups	-14 * SIZE(X), %xmm5
	vmovups	-12 * SIZE(X), %xmm6
	vmovups	-10 * SIZE(X), %xmm7

	vmovups	 -8 * SIZE(X), %xmm8
	vmovups	 -6 * SIZE(X), %xmm9
	vmovups	 -4 * SIZE(X), %xmm10
	vmovups	 -2 * SIZE(X), %xmm11

	decq	%rax
	jle	L(12)

	ALIGN_3

L(11):
	prefetchnta     A_PRE(Y)

	vfmaddpd %xmm0 , -16 * SIZE(Y), %xmm4 , %xmm0
	vfmaddpd %xmm1 , -14 * SIZE(Y), %xmm5 , %xmm1
	prefetchnta     A_PRE(X)
	vfmaddpd %xmm2 , -12 * SIZE(Y), %xmm6 , %xmm2
	vfmaddpd %xmm3 , -10 * SIZE(Y), %xmm7 , %xmm3

	vmovups	  0 * SIZE(X), %xmm4
	vfmaddpd %xmm0 ,  -8 * SIZE(Y), %xmm8 , %xmm0
	vfmaddpd %xmm1 ,  -6 * SIZE(Y), %xmm9 , %xmm1
	vmovups	  2 * SIZE(X), %xmm5
	vmovups	  4 * SIZE(X), %xmm6
	vfmaddpd %xmm2 ,  -4 * SIZE(Y), %xmm10, %xmm2
	vfmaddpd %xmm3 ,  -2 * SIZE(Y), %xmm11, %xmm3
	vmovups	  6 * SIZE(X), %xmm7

	prefetchnta     A_PRE+64(Y)

	vmovups	  8 * SIZE(X), %xmm8
	vmovups	 10 * SIZE(X), %xmm9
	prefetchnta     A_PRE+64(X)
	vmovups	 12 * SIZE(X), %xmm10
	vmovups	 14 * SIZE(X), %xmm11

	subq	$-16 * SIZE, X
	subq	$-16 * SIZE, Y

	decq	%rax
	jg	L(11)
	ALIGN_3

L(12):

	vfmaddpd %xmm0 , -16 * SIZE(Y), %xmm4 , %xmm0
	vfmaddpd %xmm1 , -14 * SIZE(Y), %xmm5 , %xmm1
	vfmaddpd %xmm2 , -12 * SIZE(Y), %xmm6 , %xmm2
	vfmaddpd %xmm3 , -10 * SIZE(Y), %xmm7 , %xmm3

	vfmaddpd %xmm0 ,  -8 * SIZE(Y), %xmm8 , %xmm0
	vfmaddpd %xmm1 ,  -6 * SIZE(Y), %xmm9 , %xmm1
	vfmaddpd %xmm2 ,  -4 * SIZE(Y), %xmm10, %xmm2
	vfmaddpd %xmm3 ,  -2 * SIZE(Y), %xmm11, %xmm3

	subq	$-16 * SIZE, X
	subq	$-16 * SIZE, Y
	ALIGN_3

L(14):
	testq	$15, N
	jle	L(999)

	testq	$8, N
	jle	L(15)

	vmovups	-16 * SIZE(X), %xmm4
	vmovups	-14 * SIZE(X), %xmm5
	vmovups	-12 * SIZE(X), %xmm6
	vmovups	-10 * SIZE(X), %xmm7

	vfmaddpd %xmm0 , -16 * SIZE(Y), %xmm4 , %xmm0
	vfmaddpd %xmm1 , -14 * SIZE(Y), %xmm5 , %xmm1
	vfmaddpd %xmm2 , -12 * SIZE(Y), %xmm6 , %xmm2
	vfmaddpd %xmm3 , -10 * SIZE(Y), %xmm7 , %xmm3

	addq	$8 * SIZE, X
	addq	$8 * SIZE, Y
	ALIGN_3

L(15):
	testq	$4, N
	jle	L(16)

	vmovups	-16 * SIZE(X), %xmm4
	vmovups	-14 * SIZE(X), %xmm5

	vfmaddpd %xmm0 , -16 * SIZE(Y), %xmm4 , %xmm0
	vfmaddpd %xmm1 , -14 * SIZE(Y), %xmm5 , %xmm1

	addq	$4 * SIZE, X
	addq	$4 * SIZE, Y
	ALIGN_3

L(16):
	testq	$2, N
	jle	L(17)

	vmovups	-16 * SIZE(X), %xmm4
	vfmaddpd %xmm0 , -16 * SIZE(Y), %xmm4 , %xmm0


	addq	$2 * SIZE, X
	addq	$2 * SIZE, Y
	ALIGN_3

L(17):
	testq	$1, N
	jle	L(999)

	vmovsd	-16 * SIZE(X), %xmm4
	vmovsd	-16 * SIZE(Y), %xmm5
	vfmaddpd %xmm0, %xmm4 , %xmm5 , %xmm0
	jmp	L(999)
	ALIGN_3


L(50):
	movq	N,  %rax
	sarq	$3, %rax
	jle	L(55)
	ALIGN_3

L(53):


	vmovsd	0 * SIZE(X), %xmm4
	addq	INCX, X
	vmovsd	0 * SIZE(Y), %xmm8
	addq	INCY, Y
	vmovsd	0 * SIZE(X), %xmm5
	addq	INCX, X
	vmovsd	0 * SIZE(Y), %xmm9
	addq	INCY, Y

	vmovsd	0 * SIZE(X), %xmm6
	addq	INCX, X
	vmovsd	0 * SIZE(Y), %xmm10
	addq	INCY, Y
	vmovsd	0 * SIZE(X), %xmm7
	addq	INCX, X
	vmovsd	0 * SIZE(Y), %xmm11
	addq	INCY, Y

	vfmaddpd %xmm0 , %xmm4 , %xmm8 , %xmm0
	vfmaddpd %xmm1 , %xmm5 , %xmm9 , %xmm1
	vfmaddpd %xmm2 , %xmm6 , %xmm10, %xmm2
	vfmaddpd %xmm3 , %xmm7 , %xmm11, %xmm3


	vmovsd	0 * SIZE(X), %xmm4
	addq	INCX, X
	vmovsd	0 * SIZE(Y), %xmm8
	addq	INCY, Y
	vmovsd	0 * SIZE(X), %xmm5
	addq	INCX, X
	vmovsd	0 * SIZE(Y), %xmm9
	addq	INCY, Y

	vmovsd	0 * SIZE(X), %xmm6
	addq	INCX, X
	vmovsd	0 * SIZE(Y), %xmm10
	addq	INCY, Y
	vmovsd	0 * SIZE(X), %xmm7
	addq	INCX, X
	vmovsd	0 * SIZE(Y), %xmm11
	addq	INCY, Y

	vfmaddpd %xmm0 , %xmm4 , %xmm8 , %xmm0
	vfmaddpd %xmm1 , %xmm5 , %xmm9 , %xmm1
	vfmaddpd %xmm2 , %xmm6 , %xmm10, %xmm2
	vfmaddpd %xmm3 , %xmm7 , %xmm11, %xmm3

	decq	%rax
	jg	L(53)
	ALIGN_3

L(55):
	movq	N, %rax
	andq	$7,   %rax
	jle	L(999)
	ALIGN_3

L(56):
	vmovsd	0 * SIZE(X), %xmm4
	addq	INCX, X
	vmovsd	0 * SIZE(Y), %xmm8
	addq	INCY, Y

	vfmaddpd %xmm0 , %xmm4 , %xmm8 , %xmm0

	decq	%rax
	jg	L(56)
	ALIGN_3

L(999):
	vaddpd	%xmm1, %xmm0 , %xmm0
	vaddpd	%xmm3, %xmm2 , %xmm2
	vaddpd	%xmm2, %xmm0 , %xmm0

	vhaddpd	%xmm0, %xmm0 , %xmm0

	RESTOREREGISTERS

	ret

	EPILOGUE
