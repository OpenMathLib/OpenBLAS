/*****************************************************************************
Copyright (c) 2011, Lab of Parallel Software and Computational Science,ICSAS
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.
   3. Neither the name of the ISCAS nor the names of its contributors may 
      be used to endorse or promote products derived from this software 
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

**********************************************************************************/

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

		
#define PREFETCH_DISTANCE 1864
		
#define N	$4

#define X	$8
#define INCX	$9

#define Y	$10
#define INCY	$11

#define I	$2
#define TEMP	$3

#define YY	$5

#define ALPHA	$f15

#define a1	$f0
#define a2	$f1
#define a3	$f2
#define a4	$f3
#define a5	$f4
#define a6	$f5
#define a7	$f6
#define a8	$f7

#define a9	$f8
#define a10	$f9
#define a11	$f10
#define a12	$f11
#define a13	$f12
#define a14	$f13
#define a15	$f14
#define a16	$f17

#define t1	$f18
#define t2	$f19
#define t3	$f20
#define t4	$f21

#define b1	$f22
#define b2	$f23
#define b3	$f24
#define b4	$f25

#define b5	$f26
#define b6	$f27
#define b7	$f28
#define b8	$f29


#define A1	 0
#define A2	 1
#define A3	 2
#define A4	 3
#define A5	 4
#define A6	 5
#define A7	 6
#define A8	 7

#define A9	 8
#define A10	 9
#define A11	 10
#define A12	 11
#define A13	 12
#define A14	 13
#define A15	 14
#define A16	 17

#define T1	 18
#define T2	 19
#define T3	 20
#define T4	 21

#define B1	 22
#define B2	 23
#define B3	 24
#define B4	 25

#define B5	 26
#define B6	 27
#define B7	 28
#define B8	 29

#define X_BASE 8
#define Y_BASE 10
		
#define gsLQC1_(base,fq,ft,offset) .word (0x32<<26|base<<21|ft<<16|0x1<<15|offset<<6|0x1<<5|fq)
#define gsLQC1(base,fq,ft,offset) gsLQC1_((base), (fq), (ft), (offset))

#define gsSQC1_(base,fq,ft,offset) .word (0x3A<<26|base<<21|ft<<16|0x1<<15|offset<<6|0x1<<5|fq)
#define gsSQC1(base,fq,ft,offset) gsSQC1_((base), (fq), (ft), (offset))

	PROLOGUE
	
#ifndef __64BIT__
	daddiu	$sp, $sp, -40
	sdc1	$f20, 0($sp)
	sdc1	$f22, 8($sp)
	sdc1	$f24, 16($sp)
	sdc1	$f26, 24($sp)
	sdc1	$f28, 32($sp)
#else
	daddiu	$sp, $sp, -48
	sdc1	$f24, 0($sp)
	sdc1	$f25, 8($sp)
	sdc1	$f26, 16($sp)
	sdc1	$f27, 24($sp)
	sdc1	$f28, 32($sp)
	sdc1	$f29, 40($sp)
#endif


	
	li	TEMP, SIZE

	blez	N, .L999
	dsll	INCX, INCX, BASE_SHIFT

	bne	INCX, TEMP, .L20
	dsll	INCY, INCY, BASE_SHIFT

	bne	INCY, TEMP, .L20
	dsra	I, N, 4

	blez	I, .L15
	daddiu	I, I, -1

	gsLQC1(X_BASE,A2,A1,0*SIZE)
	gsLQC1(X_BASE,A4,A3,2*SIZE)
	gsLQC1(X_BASE,A6,A5,4*SIZE)
	gsLQC1(X_BASE,A8,A7,6*SIZE)

	gsLQC1(X_BASE,A10,A9,8*SIZE)
	gsLQC1(X_BASE,A12,A11,10*SIZE)
	gsLQC1(X_BASE,A14,A13,12*SIZE)
	gsLQC1(X_BASE,A16,A15,14*SIZE)

	gsLQC1(Y_BASE,B2,B1,0*SIZE)
	gsLQC1(Y_BASE,B4,B3,2*SIZE)
	gsLQC1(Y_BASE,B6,B5,4*SIZE)
	gsLQC1(Y_BASE,B8,B7,6*SIZE)
		
	blez	I, .L13
	NOP
	.align 5

.L12:
		
	MADD	t1, b1, ALPHA, a1	
	MADD	t2, b2, ALPHA, a2
	gsSQC1(Y_BASE, T2, T1, 0*SIZE)		
	gsLQC1(Y_BASE,B2,B1,8*SIZE)

	MADD	t3, b3, ALPHA, a3
	MADD	t4, b4, ALPHA, a4
	gsSQC1(Y_BASE, T4, T3, 2*SIZE)
	gsLQC1(Y_BASE,B4,B3,10*SIZE)

	PREFETCHD(PREFETCH_DISTANCE*SIZE(Y))
	PREFETCHD((PREFETCH_DISTANCE+4)*SIZE(Y))

	MADD	t1, b5, ALPHA, a5
	MADD	t2, b6, ALPHA, a6
	gsSQC1(Y_BASE, T2, T1, 4*SIZE)		
	gsLQC1(Y_BASE,B6,B5,12*SIZE)

	MADD	t3, b7, ALPHA, a7
	MADD	t4, b8, ALPHA, a8
	gsSQC1(Y_BASE, T4, T3, 6*SIZE)
	gsLQC1(Y_BASE,B8,B7,14*SIZE)

	PREFETCHD((PREFETCH_DISTANCE+8)*SIZE(Y))
	PREFETCHD((PREFETCH_DISTANCE+12)*SIZE(Y))

	MADD	t1, b1, ALPHA, a9	
	MADD	t2, b2, ALPHA, a10
	gsSQC1(Y_BASE, T2, T1, 8*SIZE)		
	gsLQC1(Y_BASE,B2,B1,16*SIZE)

	MADD	t3, b3, ALPHA, a11
	MADD	t4, b4, ALPHA, a12
	gsSQC1(Y_BASE, T4, T3, 10*SIZE)
	gsLQC1(Y_BASE,B4,B3,18*SIZE)

	PREFETCHD(PREFETCH_DISTANCE*SIZE(X))
	PREFETCHD((PREFETCH_DISTANCE+4)*SIZE(X))

	MADD	t1, b5, ALPHA, a13	
	MADD	t2, b6, ALPHA, a14
	gsSQC1(Y_BASE, T2, T1, 12*SIZE)		
	gsLQC1(Y_BASE,B6,B5,20*SIZE)

	MADD	t3, b7, ALPHA, a15
	MADD	t4, b8, ALPHA, a16
	gsSQC1(Y_BASE, T4, T3, 14*SIZE)
	gsLQC1(Y_BASE,B8,B7,22*SIZE)
		
	PREFETCHD((PREFETCH_DISTANCE+8)*SIZE(X))
	PREFETCHD((PREFETCH_DISTANCE+12)*SIZE(X))

	gsLQC1(X_BASE,A2,A1,16*SIZE)
	gsLQC1(X_BASE,A4,A3,18*SIZE)
	gsLQC1(X_BASE,A6,A5,20*SIZE)
	gsLQC1(X_BASE,A8,A7,22*SIZE)

	gsLQC1(X_BASE,A10,A9,24*SIZE)
	gsLQC1(X_BASE,A12,A11,26*SIZE)
	gsLQC1(X_BASE,A14,A13,28*SIZE)
	gsLQC1(X_BASE,A16,A15,30*SIZE)


	daddiu	I, I, -1
	daddiu	Y, Y, 16 * SIZE
		
	daddiu	X, X, 16 * SIZE
	bgtz	I, .L12

	.align 5

.L13:

	MADD	t1, b1, ALPHA, a1	
	MADD	t2, b2, ALPHA, a2
	gsSQC1(Y_BASE, T2, T1, 0*SIZE)		
	gsLQC1(Y_BASE,B2,B1,8*SIZE)

	MADD	t3, b3, ALPHA, a3
	MADD	t4, b4, ALPHA, a4
	gsSQC1(Y_BASE, T4, T3, 2*SIZE)
	gsLQC1(Y_BASE,B4,B3,10*SIZE)


	MADD	t1, b5, ALPHA, a5
	MADD	t2, b6, ALPHA, a6
	gsSQC1(Y_BASE, T2, T1, 4*SIZE)		
	gsLQC1(Y_BASE,B6,B5,12*SIZE)

	MADD	t3, b7, ALPHA, a7
	MADD	t4, b8, ALPHA, a8
	gsSQC1(Y_BASE, T4, T3, 6*SIZE)
	gsLQC1(Y_BASE,B8,B7,14*SIZE)


	MADD	t1, b1, ALPHA, a9	
	MADD	t2, b2, ALPHA, a10
	gsSQC1(Y_BASE, T2, T1, 8*SIZE)		


	MADD	t3, b3, ALPHA, a11
	MADD	t4, b4, ALPHA, a12
	gsSQC1(Y_BASE, T4, T3, 10*SIZE)


	MADD	t1, b5, ALPHA, a13	
	MADD	t2, b6, ALPHA, a14
	gsSQC1(Y_BASE, T2, T1, 12*SIZE)		


	MADD	t3, b7, ALPHA, a15
	MADD	t4, b8, ALPHA, a16
	gsSQC1(Y_BASE, T4, T3, 14*SIZE)


	daddiu	X, X, 16 * SIZE
	daddiu	Y, Y, 16 * SIZE
	.align 5

.L15:
	andi	I,  N, 15

	blez	I, .L999
	NOP
	.align	3

.L16:
	LD	a1,  0 * SIZE(X)
	LD	b1,  0 * SIZE(Y)

	daddiu	X, X, SIZE
	daddiu	Y, Y, SIZE

	MADD	t1, b1, ALPHA, a1
	daddiu	I, I, -1

	bgtz	I, .L16
	ST	t1, -1 * SIZE(Y)


#ifndef __64BIT__
	ldc1	$f20, 0($sp)
	ldc1	$f22, 8($sp)
	ldc1	$f24, 16($sp)
	ldc1	$f26, 24($sp)
	ldc1	$f28, 32($sp)
	daddiu	$sp, $sp, 40
#else
	ldc1	$f24, 0($sp)
	ldc1	$f25, 8($sp)
	ldc1	$f26, 16($sp)
	ldc1	$f27, 24($sp)
	ldc1	$f28, 32($sp)
	ldc1	$f29, 40($sp)
	daddiu	$sp, $sp, 48
#endif

	j	$31
	NOP
	.align 5

.L20:
	dsra	I, N, 3
	move	YY, Y

	blez	I, .L25
	daddiu	I, I, -1

	LD	a1,  0 * SIZE(X)
	daddu	X, X, INCX
	LD	b1,  0 * SIZE(Y)
	daddu	Y, Y, INCY
	LD	a2,  0 * SIZE(X)
	daddu	X, X, INCX
	LD	b2,  0 * SIZE(Y)
	daddu	Y, Y, INCY
	LD	a3,  0 * SIZE(X)
	daddu	X, X, INCX
	LD	b3,  0 * SIZE(Y)
	daddu	Y, Y, INCY
	LD	a4,  0 * SIZE(X)
	daddu	X, X, INCX
	LD	b4,  0 * SIZE(Y)
	daddu	Y, Y, INCY
	LD	a5,  0 * SIZE(X)
	daddu	X, X, INCX
	LD	b5,  0 * SIZE(Y)
	daddu	Y, Y, INCY
	LD	a6,  0 * SIZE(X)
	daddu	X, X, INCX
	LD	b6,  0 * SIZE(Y)
	daddu	Y, Y, INCY
	LD	a7,  0 * SIZE(X)
	daddu	X, X, INCX
	LD	b7,  0 * SIZE(Y)
	daddu	Y, Y, INCY
	LD	a8,  0 * SIZE(X)
	daddu	X, X, INCX
	LD	b8,  0 * SIZE(Y)
	daddu	Y, Y, INCY

	blez	I, .L23
	NOP
	.align 5

.L22:
	MADD	t1, b1, ALPHA, a1
	LD	a1,  0 * SIZE(X)
	LD	b1,  0 * SIZE(Y)
	daddu	X, X, INCX
	daddu	Y, Y, INCY

	MADD	t2, b2, ALPHA, a2
	LD	a2,  0 * SIZE(X)
	LD	b2,  0 * SIZE(Y)
	daddu	X, X, INCX
	daddu	Y, Y, INCY

	MADD	t3, b3, ALPHA, a3
	LD	a3,  0 * SIZE(X)
	LD	b3,  0 * SIZE(Y)
	daddu	X, X, INCX
	daddu	Y, Y, INCY

	MADD	t4, b4, ALPHA, a4
	LD	a4,  0 * SIZE(X)
	LD	b4,  0 * SIZE(Y)
	daddu	X, X, INCX
	daddu	Y, Y, INCY

	ST	t1,  0 * SIZE(YY)
	daddu	YY, YY, INCY
	MADD	t1, b5, ALPHA, a5

	LD	a5,  0 * SIZE(X)
	LD	b5,  0 * SIZE(Y)
	daddu	X, X, INCX
	daddu	Y, Y, INCY

	ST	t2,  0 * SIZE(YY)
	daddu	YY, YY, INCY
	MADD	t2, b6, ALPHA, a6

	LD	a6,  0 * SIZE(X)
	LD	b6,  0 * SIZE(Y)
	daddu	X, X, INCX
	daddu	Y, Y, INCY

	ST	t3,  0 * SIZE(YY)
	daddu	YY, YY, INCY
	MADD	t3, b7, ALPHA, a7

	LD	a7,  0 * SIZE(X)
	LD	b7,  0 * SIZE(Y)
	daddu	X, X, INCX
	daddu	Y, Y, INCY

	ST	t4,  0 * SIZE(YY)
	daddu	YY, YY, INCY
	MADD	t4, b8, ALPHA, a8

	LD	a8,  0 * SIZE(X)
	daddu	X, X, INCX

	LD	b8,  0 * SIZE(Y)
	daddu	Y, Y, INCY

	ST	t1,  0 * SIZE(YY)
	daddu	YY, YY, INCY
	ST	t2,  0 * SIZE(YY)
	daddu	YY, YY, INCY
	ST	t3,  0 * SIZE(YY)
	daddu	YY, YY, INCY
	ST	t4,  0 * SIZE(YY)
	daddiu	I, I, -1

	bgtz	I, .L22
	daddu	YY, YY, INCY
	.align 5

.L23:
	MADD	t1, b1, ALPHA, a1
	MADD	t2, b2, ALPHA, a2
	MADD	t3, b3, ALPHA, a3
	MADD	t4, b4, ALPHA, a4

	ST	t1,  0 * SIZE(YY)
	daddu	YY, YY, INCY
	MADD	t1, b5, ALPHA, a5

	ST	t2,  0 * SIZE(YY)
	daddu	YY, YY, INCY
	MADD	t2, b6, ALPHA, a6

	ST	t3,  0 * SIZE(YY)
	daddu	YY, YY, INCY
	MADD	t3, b7, ALPHA, a7

	ST	t4,  0 * SIZE(YY)
	daddu	YY, YY, INCY
	MADD	t4, b8, ALPHA, a8

	ST	t1,  0 * SIZE(YY)
	daddu	YY, YY, INCY
	ST	t2,  0 * SIZE(YY)
	daddu	YY, YY, INCY
	ST	t3,  0 * SIZE(YY)
	daddu	YY, YY, INCY
	ST	t4,  0 * SIZE(YY)
	daddu	YY, YY, INCY
	.align 5

.L25:
	andi	I,  N, 7

	blez	I, .L999
	NOP
	.align	3

.L26:
	LD	a1,  0 * SIZE(X)
	LD	b1,  0 * SIZE(Y)

	MADD	t1, b1, ALPHA, a1
	daddu	X, X, INCX

	ST	t1,  0 * SIZE(Y)
	daddiu	I, I, -1

	bgtz	I, .L26
	daddu	Y, Y, INCY
	.align 5

.L999:

#ifndef __64BIT__
	ldc1	$f20, 0($sp)
	ldc1	$f22, 8($sp)
	ldc1	$f24, 16($sp)
	ldc1	$f26, 24($sp)
	ldc1	$f28, 32($sp)
	daddiu	$sp, $sp, 40
#else
	ldc1	$f24, 0($sp)
	ldc1	$f25, 8($sp)
	ldc1	$f26, 16($sp)
	ldc1	$f27, 24($sp)
	ldc1	$f28, 32($sp)
	ldc1	$f29, 40($sp)
	daddiu	$sp, $sp, 48
#endif

	j	$31
	NOP

	EPILOGUE
