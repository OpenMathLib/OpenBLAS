/*******************************************************************************
Copyright (c) 2017, The OpenBLAS Project
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
*******************************************************************************/

#define ASSEMBLER
#include "common.h"

#define	N	x0	/* vector length */
#define	X	x1	/* X vector address */
#define	INC_X	x2	/* X stride */
#define	Y	x3	/* Y vector address */
#define	INC_Y	x4	/* Y stride */
#define I	x5	/* loop variable */

/*******************************************************************************
* Macro definitions
*******************************************************************************/

#define REG0	xzr
#define DOTF	d0
#define TMPX	d16
#define LD1VX	{v16.d}[0]
#define TMPY	d24
#define LD1VY	{v24.d}[0]
#define SZ	8

/******************************************************************************/

.macro KERNEL_F1
	ldr	TMPX, [X]
	ldr	TMPY, [Y]
	add	X, X, #SZ
	add	Y, Y, #SZ
	fmadd	DOTF, TMPX, TMPY, DOTF
.endm

.macro KERNEL_F16
	ldp	q16, q17, [X]
	ldp	q24, q25, [Y]

	ldp	q18, q19, [X, #32]
	ldp	q26, q27, [Y, #32]

	fmla	v0.2d, v16.2d, v24.2d
	fmla	v1.2d, v17.2d, v25.2d

	ldp	q20, q21, [X, #64]
	ldp	q28, q29, [Y, #64]

	fmla	v2.2d, v18.2d, v26.2d
	fmla	v3.2d, v19.2d, v27.2d

	ldp	q22, q23, [X, #96]
	ldp	q30, q31, [Y, #96]

	add	Y, Y, #128
	add	X, X, #128

	fmla	v4.2d, v20.2d, v28.2d
	fmla	v5.2d, v21.2d, v29.2d

	PRFM	PLDL1KEEP, [X, #896]
	PRFM	PLDL1KEEP, [Y, #896]
	PRFM	PLDL1KEEP, [X, #896+64]
	PRFM	PLDL1KEEP, [Y, #896+64]

	fmla	v6.2d, v22.2d, v30.2d
	fmla	v7.2d, v23.2d, v31.2d
.endm

.macro KERNEL_F32
	KERNEL_F16
	KERNEL_F16
.endm

.macro KERNEL_F32_FINALIZE
	fadd	v0.2d, v0.2d, v1.2d
	fadd	v2.2d, v2.2d, v3.2d
	fadd	v4.2d, v4.2d, v5.2d
	fadd	v6.2d, v6.2d, v7.2d
	fadd	v0.2d, v0.2d, v2.2d
	fadd	v4.2d, v4.2d, v6.2d
	fadd	v0.2d, v0.2d, v4.2d
	faddp	DOTF, v0.2d
.endm

.macro INIT_S
	lsl	INC_X, INC_X, #3
	lsl	INC_Y, INC_Y, #3
.endm

.macro KERNEL_S1
	ld1	LD1VX, [X], INC_X
	ld1	LD1VY, [Y], INC_Y
	fmadd	DOTF, TMPX, TMPY, DOTF
.endm

/*******************************************************************************
* End of macro definitions
*******************************************************************************/

	PROLOGUE

	fmov	DOTF, REG0
	fmov	d1, REG0
	fmov	d2, REG0
	fmov	d3, REG0
	fmov	d4, REG0
	fmov	d5, REG0
	fmov	d6, REG0
	fmov	d7, REG0

	cmp	N, xzr
	ble	dot_kernel_L999

	cmp	INC_X, #1
	bne	dot_kernel_S_BEGIN
	cmp	INC_Y, #1
	bne	dot_kernel_S_BEGIN

dot_kernel_F_BEGIN:

	asr	I, N, #5
	cmp	I, xzr
	beq	dot_kernel_F1

dot_kernel_F32:

	KERNEL_F32

	subs	I, I, #1
	bne	dot_kernel_F32

	KERNEL_F32_FINALIZE

dot_kernel_F1:

	ands	I, N, #31
	ble	dot_kernel_L999

dot_kernel_F10:

	KERNEL_F1

	subs    I, I, #1
        bne     dot_kernel_F10

	ret

dot_kernel_S_BEGIN:

	INIT_S

	asr	I, N, #2
	cmp	I, xzr
	ble	dot_kernel_S1

dot_kernel_S4:

	KERNEL_S1
	KERNEL_S1
	KERNEL_S1
	KERNEL_S1

	subs	I, I, #1
	bne	dot_kernel_S4

dot_kernel_S1:

	ands	I, N, #3
	ble	dot_kernel_L999

dot_kernel_S10:

	KERNEL_S1

	subs    I, I, #1
        bne     dot_kernel_S10

dot_kernel_L999:

	ret

	EPILOGUE
