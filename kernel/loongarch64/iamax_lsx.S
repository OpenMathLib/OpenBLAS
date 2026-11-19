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
#define TSZ	$r19		/* bytes per real element            */

#define MAXF	$f20		/* running maximum                   */
#define FA	$f17		/* caller-saved and no vector uses $vr17 */

/* vector registers ($xrN aliases $vrN low half) */
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
#define VT0	$vr21		/* tree temps, never raw operands   */
#define VT1	$vr22
#define VT2	$vr23

#ifdef DOUBLE
#define FLD	fld.d
#define FABS	fabs.d
#define FADD	fadd.d
#define FCLT	fcmp.clt.d
#define FCEQ	fcmp.ceq.d
#define FSEL	fsel
#define VMAX	vfmax.d
#define ABSV	vbitclri.d
#define ABSI	63
#define REPL	vreplvei.d
#define BLKSC	$f9	/* = $vr9 low lane; used in scalar ops only */

#else
#define FLD	fld.s
#define FABS	fabs.s
#define FADD	fadd.s
#define FCLT	fcmp.clt.s
#define FCEQ	fcmp.ceq.s
#define FSEL	fsel
#define VMAX	vfmax.s
#define ABSV	vbitclri.w
#define ABSI	31
#define REPL	vreplvei.w
#define BLKSC	$f9	/* = $vr9 low lane; used in scalar ops only */
#endif

/* 32 real elements per vector iteration:
 *   double : 16 vld (2 per 16B register)  -> 16 magnitudes
 *   single :  8 vld (4 per 16B register)  ->  8 magnitudes          */
#define BLK_ELEMS	32
#define REM_MASK	31

/**********************************************************************/
/* SCALAR_STEP: process the element at (X, Z) against the running
 * max.  Z must hold the 1-based index of this element; X is advanced
 * by INCX, Z is NOT modified. */
	.macro	SCALAR_STEP
	FLD	FA, X, 0
	add.d	X, X, INCX
	FABS	FA, FA
	FCLT	$fcc0, MAXF, FA
	FSEL	MAXF, MAXF, FA, $fcc0
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

	slli.d	INCX, INCX, BASE_SHIFT		/* INCX in bytes */
	li.d	T0, 1
	slli.d	TSZ, T0, BASE_SHIFT		/* bytes per element */

	move	X0, X

	/* INIT: element 1 (1-based) is the initial maximum */
	addi.d	Z, $r0, 1			/* Z = 1 */
	move	INDEX, Z
	FLD	FA, X, 0
	add.d	X, X, INCX
	FABS	FA, FA
	fmov.d	MAXF, FA
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
	/* contiguous: vector block loop, 64 elements per iteration */
	srai.d	I, N, 5				/* full 32-element blocks */
	addi.d	Z, Z, 1				/* next element index     */
	beqz	I, .L25
	.align	3
.L10:
#ifdef DOUBLE
	/* 16 registers x 2 doubles = 32 elements */
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
	addi.d	X, X, 256
	ABSV	V1, V1, ABSI
	ABSV	V2, V2, ABSI
	ABSV	V3, V3, ABSI
	ABSV	V4, V4, ABSI
	ABSV	V5, V5, ABSI
	ABSV	V6, V6, ABSI
	ABSV	V7, V7, ABSI
	ABSV	V8, V8, ABSI
	ABSV	V9, V9, ABSI
	ABSV	V10, V10, ABSI
	ABSV	V11, V11, ABSI
	ABSV	V12, V12, ABSI
	ABSV	V13, V13, ABSI
	ABSV	V14, V14, ABSI
	ABSV	V15, V15, ABSI
	ABSV	V16, V16, ABSI
	/* tree-reduce 16 registers to one block max */
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
	VMAX	V1, V1, V9
#else
	/* 8 registers x 4 floats = 32 elements */
	vld	V1, X, 0
	vld	V2, X, 16
	vld	V3, X, 32
	vld	V4, X, 48
	vld	V5, X, 64
	vld	V6, X, 80
	vld	V7, X, 96
	vld	V8, X, 112
	addi.d	X, X, 128
	ABSV	V1, V1, ABSI
	ABSV	V2, V2, ABSI
	ABSV	V3, V3, ABSI
	ABSV	V4, V4, ABSI
	ABSV	V5, V5, ABSI
	ABSV	V6, V6, ABSI
	ABSV	V7, V7, ABSI
	ABSV	V8, V8, ABSI
	/* tree-reduce 8 registers to one block max */
	VMAX	V1, V1, V2
	VMAX	V3, V3, V4
	VMAX	V5, V5, V6
	VMAX	V7, V7, V8
	VMAX	V1, V1, V3
	VMAX	V5, V5, V7
	VMAX	V1, V1, V5
#endif
	/* horizontal: fold the lanes down to lane 0 (128-bit regs) */
#ifdef DOUBLE
	/* two lanes: replicate lane 1 and fold */
	REPL	VT0, V1, 1
	VMAX	$vr9, V1, VT0			/* block max in lane 0 */
#else
	/* Four lanes.  Every lane is replicated out of the ORIGINAL
	 * register before any vfmax touches it: a vfmax that reads the
	 * accumulator also writes all of its lanes, so folding lane 2 in
	 * first would destroy the lane 3 value still to be merged. */
	REPL	VT0, V1, 1
	REPL	VT1, V1, 2
	REPL	VT2, V1, 3
	VMAX	V1, V1, VT0
	VMAX	V1, V1, VT1
	VMAX	$vr9, V1, VT2			/* block max in lane 0 */
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
	add.d	X, X, TSZ
	FABS	FA, FA
	FCEQ	$fcc0, FA, MAXF			/* == global max ?           */
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
