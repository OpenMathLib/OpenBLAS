/***************************************************************************
Copyright (c) 2023, The OpenBLAS Project
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
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/
#define ASSEMBLER

#include "common.h"

/**********************************************************************/
/* scalar registers */
#define N	$r4		/* vector length (also return value) */
#define X	$r5		/* "X" vector address                */
#define INCX	$r6		/* "X" stride (bytes, after scaling) */
#define I	$r12		/* vector loop counter               */
#define J	$r13		/* scalar remainder counter          */
#define Z	$r14		/* 1-based element index cursor      */
#define T0	$r15
#define INDEX	$r16		/* index of max value                */
#define X0	$r17		/* original X pointer                */
#define T1	$r18
#define TSZ	$r19		/* bytes per complex element         */

#define MAXF	$f20		/* running maximum ($vr20 low lane)  */
#define FA	$f10
#define FB	$f11
#define FC	$f12

/* vector registers (all in the caller-saved range, none >= $vr24) */
#define V1	$vr1
#define V2	$vr2
#define V3	$vr3
#define V4	$vr4
#define V5	$vr5
#define V6	$vr6
#define V7	$vr7
#define V8	$vr8
#define V9	$vr9
#define V10	$vr10
#define V11	$vr11
#define V12	$vr12
#define V13	$vr13
#define V14	$vr14
#define V15	$vr15
#define V16	$vr16
#define VT0	$vr17
#define VT1	$vr18
#define VT2	$vr19	/* extra raw-load scratch: never a MAG2 operand */

#ifdef DOUBLE
#define FLD	fld.d
#define FABS	fabs.d
#define FADD	fadd.d
#define FCLT	fcmp.clt.d
#define FCEQ	fcmp.ceq.d
#define FSEL	fsel
#define VMAX	vfmax.d
#define BLKSC	$f9	/* = $vr9 low lane; used in scalar ops only */

#else
#define FLD	fld.s
#define FABS	fabs.s
#define FADD	fadd.s
#define FCLT	fcmp.clt.s
#define FCEQ	fcmp.ceq.s
#define FSEL	fsel
#define VMAX	vfmax.s
#define BLKSC	$f1
#endif

/* 32 complex elements per vector iteration:
 *   double : 32 vld  (1 complex per 16B register, 6 load rounds)
 *   single : 16 vld  (2 complex per 16B register)                    */
#define BLK_ELEMS	32
#define REM_MASK	31

/**********************************************************************/
/* MAG2: compute cabs1 into A.
 *   double : A=[re_a,im_a], B=[re_b,im_b]  -> A=[c_a, c_b]
 *   single : A=[re0,im0,re1,im1]           -> A=[c0, c1, c0, c1]      */
#ifdef DOUBLE
	.macro	MAG2	A, B
	vpickev.d	VT0, \B, \A
	vpickod.d	VT1, \B, \A
	vbitclri.d	VT0, VT0, 63
	vbitclri.d	VT1, VT1, 63
	vfadd.d		\A, VT0, VT1
	.endm
#else
	.macro	MAG2	A
	vbitclri.w	\A, \A, 31
	vpickev.w	VT0, \A, \A
	vpickod.w	VT1, \A, \A
	vfadd.s		\A, VT0, VT1
	.endm
#endif

/* SCALAR_STEP: process the complex element at (X, Z) against the
 * running max.  Z must hold the 1-based index of this element; X is
 * advanced by INCX, Z is NOT modified. */
	.macro	SCALAR_STEP
	FLD	FA, X, 0
	FLD	FB, X, SIZE
	add.d	X, X, INCX
	FABS	FA, FA
	FABS	FB, FB
	FADD	FC, FA, FB
	FCLT	$fcc0, MAXF, FC
	FSEL	MAXF, MAXF, FC, $fcc0
	bceqz	$fcc0, 9f
	move	INDEX, Z
9:
	.endm

/**********************************************************************/
	PROLOGUE
	PROFCODE

	li.d	INDEX, 0
	bge	$r0, N, .L999
	bge	$r0, INCX, .L999

	slli.d	INCX, INCX, ZBASE_SHIFT		/* INCX in bytes */
	li.d	T0, 1
	slli.d	TSZ, T0, ZBASE_SHIFT		/* bytes per complex element */

	move	X0, X

	/* INIT: element 1 (1-based) is the initial maximum */
	addi.d	Z, $r0, 1			/* Z = 1 */
	move	INDEX, Z
	FLD	FA, X, 0
	FLD	FB, X, SIZE
	add.d	X, X, INCX
	FABS	FA, FA
	FABS	FB, FB
	FADD	MAXF, FA, FB
	addi.d	N, N, -1			/* remaining = n - 1 */
	bge	$r0, N, .L999r

	bne	INCX, TSZ, .L60			/* INCX != 1 -> scalar path */

	/* peel until X is 16-byte aligned (vld natural alignment) */
1:	andi	T0, X, 15
	beqz	T0, 2f
	addi.d	Z, Z, 1
	SCALAR_STEP
	addi.d	N, N, -1
	bge	$r0, N, .L999r
	b	1b
2:
	/* contiguous: vector block loop, 32 elements per iteration */
	srai.d	I, N, 5				/* full 32-element blocks */
	addi.d	Z, Z, 1				/* next element index     */
	beqz	I, .L25
	.align	3
.L10:
#ifdef DOUBLE
	/* phase 1: elements 1-16 -> magnitudes in V1,V3,...,V15 */
	vld	V1, X, 0
	vld	V2, X, 16
	vld	V3, X, 32
	vld	V4, X, 48
	vld	V5, X, 64
	vld	V6, X, 80
	vld	V7, X, 96
	vld	V8, X, 112
	vld	V9, X, 128
	vld	V10, X, 144
	vld	V11, X, 160
	vld	V12, X, 176
	vld	V13, X, 192
	vld	V14, X, 208
	vld	V15, X, 224
	vld	V16, X, 240
	MAG2	V1, V2
	MAG2	V3, V4
	MAG2	V5, V6
	MAG2	V7, V8
	MAG2	V9, V10
	MAG2	V11, V12
	MAG2	V13, V14
	MAG2	V15, V16
	/* phase 2: elements 17-20 -> V2, V6 */
	vld	V2, X, 256
	vld	V4, X, 272
	vld	V6, X, 288
	vld	V8, X, 304
	MAG2	V2, V4
	MAG2	V6, V8
	/* phase 3: elements 21-24 -> V4, V10 */
	vld	V4, X, 320
	vld	V8, X, 336
	vld	V10, X, 352
	vld	V12, X, 368
	MAG2	V4, V8
	MAG2	V10, V12
	/* phase 4: elements 25-28 -> V8, V14 */
	vld	V8, X, 384
	vld	V12, X, 400
	vld	V14, X, 416
	vld	V16, X, 432
	MAG2	V8, V12
	MAG2	V14, V16
	/* phase 5: elements 29-30 -> V12 */
	vld	V12, X, 448
	vld	V16, X, 464
	MAG2	V12, V16
	/* phase 6: elements 31-32 -> V16 */
	vld	V16, X, 480
	vld	VT2, X, 496
	addi.d	X, X, 512
	MAG2	V16, VT2
#else
	vld	V1, X, 0
	vld	V2, X, 16
	vld	V3, X, 32
	vld	V4, X, 48
	vld	V5, X, 64
	vld	V6, X, 80
	vld	V7, X, 96
	vld	V8, X, 112
	MAG2	V1
	MAG2	V2
	MAG2	V3
	MAG2	V4
	vld	V9, X, 128
	vld	V10, X, 144
	vld	V11, X, 160
	vld	V12, X, 176
	MAG2	V5
	MAG2	V6
	MAG2	V7
	MAG2	V8
	vld	V13, X, 192
	vld	V14, X, 208
	vld	V15, X, 224
	vld	V16, X, 240
	addi.d	X, X, 256
	MAG2	V9
	MAG2	V10
	MAG2	V11
	MAG2	V12
	MAG2	V13
	MAG2	V14
	MAG2	V15
	MAG2	V16
#endif
	/* tree-reduce 32 magnitudes (16 registers) to one block max */
	VMAX	V1, V1, V2
	VMAX	V3, V3, V4
	VMAX	V5, V5, V6
	VMAX	V7, V7, V8
	VMAX	V9, V9, V10
	VMAX	V11, V11, V12
	VMAX	V13, V13, V14
	VMAX	V15, V15, V16
	VMAX	V1, V1, V3
	VMAX	V5, V5, V7
	VMAX	V9, V9, V11
	VMAX	V13, V13, V15
	VMAX	V1, V1, V5
	VMAX	V9, V9, V13
	VMAX	V1, V1, V9			/* [mA, mB]       */
#ifdef DOUBLE
	vreplvei.d	VT0, V1, 1		/* [mB, mB]       */
	vfmax.d	$vr9, V1, VT0			/* block max      */
#else
	vreplvei.w	VT0, V1, 1
	vfmax.s	$vr1, V1, VT0			/* block max      */
#endif
	FCLT	$fcc0, MAXF, BLKSC		/* running max < block max ? */
	FSEL	MAXF, MAXF, BLKSC, $fcc0	/* strict greater -> update  */
	bceqz	$fcc0, .Lcont
	move	INDEX, Z			/* INDEX = block start       */
.Lcont:
	addi.d	Z, Z, BLK_ELEMS
	addi.d	I, I, -1
	bnez	I, .L10

	/* FALLTHROUGH: rescan the winning block, locate exact position */
.L15:
	addi.d	T1, INDEX, -1
	mul.d	T1, T1, TSZ
	add.d	X, X0, T1			/* X = winning block start   */
	move	T0, $r0
	.align	3
.L16:
	addi.d	T0, T0, 1
	addi.d	T1, T0, -BLK_ELEMS
	bge	T1, $r0, .L17			/* safety bound              */
	FLD	FA, X, 0
	FLD	FB, X, SIZE
	add.d	X, X, TSZ
	FABS	FA, FA
	FABS	FB, FB
	FADD	FC, FA, FB
	FCEQ	$fcc0, FC, MAXF			/* == global max ?           */
	bceqz	$fcc0, .L16			/* not yet -> keep scanning  */
.L17:
	addi.d	T0, T0, -1
	add.d	INDEX, INDEX, T0		/* INDEX += block offset     */
	addi.d	T1, Z, -1
	mul.d	T1, T1, TSZ
	add.d	X, X0, T1			/* X = first unconsumed elem */

	/* scalar remainder (elements after the winning block) */
.L25:
	andi	J, N, REM_MASK
	beqz	J, .L999r
	.align	3
.Lrem:
	SCALAR_STEP
	addi.d	Z, Z, 1
	addi.d	J, J, -1
	bnez	J, .Lrem
	b	.L999r

	/* INCX != 1: purely scalar, reuses the remainder loop */
	.align	3
.L60:
	addi.d	Z, Z, 1				/* next element index = 2    */
	move	J, N
	beqz	J, .L999r
	b	.Lrem

.L999r:
	move	N, INDEX
	jirl	$r0, $r1, 0x0		/* return (EPILOGUE is only .end)  */
.L999:
	move	N, $r0			/* n <= 0 or incx <= 0 -> 0 */
	jirl	$r0, $r1, 0x0
	EPILOGUE
