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
#define INDEX	$r16		/* index of min value                */
#define X0	$r17		/* original X pointer                */
#define T1	$r18
#define TSZ	$r19		/* bytes per complex element         */

#define MINF	$f20		/* running minimum                   */
#define FA	$f10
#define FB	$f11
#define FC	$f12

/* vector registers ($xrN aliases $vrN low half) */
#define V1	$xr1
#define V2	$xr2
#define V3	$xr3
#define V4	$xr4
#define V5	$xr5
#define V6	$xr6
#define V7	$xr7
#define V8	$xr8
#define V9	$xr9
#define V10	$xr10
#define V11	$xr11
#define V12	$xr12
#define V13	$xr13
#define V14	$xr14
#define V15	$xr15
#define V16	$xr16
#define RA	$xr17		/* raw-load pool: MAG2 operands only */
#define RB	$xr18
#define RC	$xr19
#define RD	$xr23	/* NOTE: $xr20 holds MINF ($f20) low lane -- never touch */
#define VT0	$xr21		/* MAG2/tree temps, never MAG2 operands */
#define VT1	$xr22

#ifdef DOUBLE
#define FLD	fld.d
#define FABS	fabs.d
#define FADD	fadd.d
#define FCLT	fcmp.clt.d
#define FCEQ	fcmp.ceq.d
#define FSEL	fsel
#define VADD	xvfadd.d
#define VMIN	xvfmin.d
#define PICK	xvpickev.d
#define PICKO	xvpickod.d
#define ABSV	xvbitclri.d
#define ABSI	63
#define PERMI	xvpermi.d
#define BLKSC	$f9	/* = $xr9 low lane; used in scalar ops only */

#else
#define FLD	fld.s
#define FABS	fabs.s
#define FADD	fadd.s
#define FCLT	fcmp.clt.s
#define FCEQ	fcmp.ceq.s
#define FSEL	fsel
#define VADD	xvfadd.s
#define VMIN	xvfmin.s
#define PICK	xvpickev.w
#define PICKO	xvpickod.w
#define ABSV	xvbitclri.w
#define ABSI	31
#define PERMI	xvpermi.d
#define BLKSC	$f1
#endif

/* 64 complex elements per vector iteration:
 *   double : 32 xvld (2 complex per 32B register, 16 MAG2)
 *   single : 16 xvld (4 complex per 32B register,  8 MAG2)          */
#define BLK_ELEMS	64
#define REM_MASK	63

/**********************************************************************/
/* MAG2: cabs1 of 4 elements (double) / 8 elements (single) into D.
 * D must differ from A, B, VT0, VT1.  Lane order is interleaved
 * ([a0,b0,a1,b1,...]) -- irrelevant, see header comment.             */
	.macro	MAG2	D, A, B
	PICK	VT0, \B, \A
	PICKO	VT1, \B, \A
	ABSV	VT0, VT0, ABSI
	ABSV	VT1, VT1, ABSI
	VADD	\D, VT0, VT1
	.endm

/* SCALAR_STEP: process the complex element at (X, Z) against the
 * running min.  Z must hold the 1-based index of this element; X is
 * advanced by INCX, Z is NOT modified. */
	.macro	SCALAR_STEP
	FLD	FA, X, 0
	FLD	FB, X, SIZE
	add.d	X, X, INCX
	FABS	FA, FA
	FABS	FB, FB
	FADD	FC, FA, FB
	FCLT	$fcc0, FC, MINF
	FSEL	MINF, MINF, FC, $fcc0
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

	/* INIT: element 1 (1-based) is the initial minimum */
	addi.d	Z, $r0, 1			/* Z = 1 */
	move	INDEX, Z
	FLD	FA, X, 0
	FLD	FB, X, SIZE
	add.d	X, X, INCX
	FABS	FA, FA
	FABS	FB, FB
	FADD	MINF, FA, FB
	addi.d	N, N, -1			/* remaining = n - 1 */
	bge	$r0, N, .L999r

	bne	INCX, TSZ, .L60			/* INCX != 1 -> scalar path */

	/* peel until X is 32-byte aligned (xvld natural alignment) */
1:	andi	T0, X, 31
	beqz	T0, 2f
	addi.d	Z, Z, 1
	SCALAR_STEP
	addi.d	N, N, -1
	bge	$r0, N, .L999r
	b	1b
2:
	/* contiguous: vector block loop, 64 elements per iteration */
	srai.d	I, N, 6				/* full 64-element blocks */
	addi.d	Z, Z, 1				/* next element index     */
	beqz	I, .L25
	.align	3
.L10:
#ifdef DOUBLE
	xvld	RA, X, 0
	xvld	RB, X, 32
	MAG2	V1, RA, RB
	xvld	RC, X, 64
	xvld	RD, X, 96
	MAG2	V2, RC, RD
	xvld	RA, X, 128
	xvld	RB, X, 160
	MAG2	V3, RA, RB
	xvld	RC, X, 192
	xvld	RD, X, 224
	MAG2	V4, RC, RD
	xvld	RA, X, 256
	xvld	RB, X, 288
	MAG2	V5, RA, RB
	xvld	RC, X, 320
	xvld	RD, X, 352
	MAG2	V6, RC, RD
	xvld	RA, X, 384
	xvld	RB, X, 416
	MAG2	V7, RA, RB
	xvld	RC, X, 448
	xvld	RD, X, 480
	MAG2	V8, RC, RD
	xvld	RA, X, 512
	xvld	RB, X, 544
	MAG2	V9, RA, RB
	xvld	RC, X, 576
	xvld	RD, X, 608
	MAG2	V10, RC, RD
	xvld	RA, X, 640
	xvld	RB, X, 672
	MAG2	V11, RA, RB
	xvld	RC, X, 704
	xvld	RD, X, 736
	MAG2	V12, RC, RD
	xvld	RA, X, 768
	xvld	RB, X, 800
	MAG2	V13, RA, RB
	xvld	RC, X, 832
	xvld	RD, X, 864
	MAG2	V14, RC, RD
	xvld	RA, X, 896
	xvld	RB, X, 928
	MAG2	V15, RA, RB
	xvld	RC, X, 960
	xvld	RD, X, 992
	addi.d	X, X, 1024
	MAG2	V16, RC, RD
#else
	xvld	RA, X, 0
	xvld	RB, X, 32
	MAG2	V1, RA, RB
	xvld	RC, X, 64
	xvld	RD, X, 96
	MAG2	V2, RC, RD
	xvld	RA, X, 128
	xvld	RB, X, 160
	MAG2	V3, RA, RB
	xvld	RC, X, 192
	xvld	RD, X, 224
	MAG2	V4, RC, RD
	xvld	RA, X, 256
	xvld	RB, X, 288
	MAG2	V5, RA, RB
	xvld	RC, X, 320
	xvld	RD, X, 352
	MAG2	V6, RC, RD
	xvld	RA, X, 384
	xvld	RB, X, 416
	MAG2	V7, RA, RB
	xvld	RC, X, 448
	xvld	RD, X, 480
	addi.d	X, X, 512
	MAG2	V8, RC, RD
#endif
	/* tree-reduce 64 magnitudes to one block min */
#ifdef DOUBLE
	VMIN	V1, V1, V2
	VMIN	V3, V3, V4
	VMIN	V5, V5, V6
	VMIN	V7, V7, V8
	VMIN	V9, V9, V10
	VMIN	V11, V11, V12
	VMIN	V13, V13, V14
	VMIN	V15, V15, V16
	VMIN	V1, V1, V3
	VMIN	V5, V5, V7
	VMIN	V9, V9, V11
	VMIN	V13, V13, V15
	VMIN	V1, V1, V5
	VMIN	V9, V9, V13
	VMIN	V1, V1, V9			/* [m0,m1,m2,m3] */
	/* horizontal: per-half pair-min, then cross-half via xvpermi.d.
	 * NOTE: xvshuf4i.* immediates are mis-encoded by binutils 2.42
	 * (0x55 produced all-zero lanes in QEMU) -- do not use it.     */
	PICK	VT0, V1, V1			/* [m0,m0 | m2,m2]  */
	PICKO	VT1, V1, V1			/* [m1,m1 | m3,m3]  */
	VMIN	V1, VT0, VT1			/* [m01,m01|m23,m23] */
	PERMI	VT0, V1, 10			/* d0 = m23         */
	VMIN	$xr9, V1, VT0			/* block min        */
#else
	VMIN	V1, V1, V2
	VMIN	V3, V3, V4
	VMIN	V5, V5, V6
	VMIN	V7, V7, V8
	VMIN	V1, V1, V3
	VMIN	V5, V5, V7
	VMIN	V1, V1, V5			/* [v0..v7] */
	PICK	VT0, V1, V1			/* [v0,v2,v0,v2|v4,v6,v4,v6] */
	PICKO	VT1, V1, V1			/* [v1,v3,v1,v3|v5,v7,v5,v7] */
	VMIN	V1, VT0, VT1			/* [s0,s1,s0,s1|s2,s3,s2,s3] */
	PICK	VT0, V1, V1			/* [s0,s0,..|s2,s2,..]       */
	PICKO	VT1, V1, V1			/* [s1,s1,..|s3,s3,..]       */
	VMIN	V1, VT0, VT1			/* [t0,..|t1,..]             */
	PERMI	VT0, V1, 10			/* d0 = t1                   */
	VMIN	$xr1, V1, VT0			/* block min                 */
#endif
	FCLT	$fcc0, BLKSC, MINF		/* block min < running min ? */
	FSEL	MINF, MINF, BLKSC, $fcc0	/* strict smaller -> update  */
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
	FCEQ	$fcc0, FC, MINF			/* == global min ?           */
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
