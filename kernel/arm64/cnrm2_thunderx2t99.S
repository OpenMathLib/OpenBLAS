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
#define I	x5	/* loop variable */

/*******************************************************************************
* Macro definitions
*******************************************************************************/

#define TMPF	d16
#define SSQ	s0
#define SSQD	d0

/******************************************************************************/

.macro INIT
	fmov	SSQD, xzr
	fmov	d1, xzr
	fmov	d2, xzr
	fmov	d3, xzr
	fmov	d4, xzr
	fmov	d5, xzr
	fmov	d6, xzr
	fmov	d7, xzr
.endm

.macro KERNEL_F1
	ldr	TMPF, [X]
	add	X, X, #8
	fcvtl	v16.2d, v16.2s
	fmla	v0.2d, v16.2d, v16.2d
.endm

.macro KERNEL_F16
	ldur	q16, [X]
	ldur	q18, [X, #16]
	ldur	q20, [X, #32]
	ldur	q22, [X, #48]
	ldur	q24, [X, #64]
	ldur	q26, [X, #80]
	ldur	q28, [X, #96]
	ldur	q30, [X, #112]

	add	X, X, #128

	fcvtl2	v17.2d, v16.4s
	fcvtl	v16.2d, v16.2s
	fcvtl2	v19.2d, v18.4s
	fcvtl	v18.2d, v18.2s
	fcvtl2	v21.2d, v20.4s
	fcvtl	v20.2d, v20.2s
	fcvtl2	v23.2d, v22.4s
	fcvtl	v22.2d, v22.2s
	fcvtl2	v25.2d, v24.4s
	fcvtl	v24.2d, v24.2s
	fcvtl2	v27.2d, v26.4s
	fcvtl	v26.2d, v26.2s
	fcvtl2	v29.2d, v28.4s
	fcvtl	v28.2d, v28.2s
	fcvtl2	v31.2d, v30.4s
	fcvtl	v30.2d, v30.2s

	fmla	v0.2d, v16.2d, v16.2d
	fmla	v1.2d, v17.2d, v17.2d
	fmla	v2.2d, v18.2d, v18.2d
	fmla	v3.2d, v19.2d, v19.2d
	fmla	v4.2d, v20.2d, v20.2d
	fmla	v5.2d, v21.2d, v21.2d
	fmla	v6.2d, v22.2d, v22.2d
	fmla	v7.2d, v23.2d, v23.2d

	fmla	v0.2d, v24.2d, v24.2d
	fmla	v1.2d, v25.2d, v25.2d
	fmla	v2.2d, v26.2d, v26.2d
	fmla	v3.2d, v27.2d, v27.2d
	fmla	v4.2d, v28.2d, v28.2d
	fmla	v5.2d, v29.2d, v29.2d
	fmla	v6.2d, v30.2d, v30.2d
	fmla	v7.2d, v31.2d, v31.2d

	prfm	PLDL1KEEP, [X, #1024]
	prfm	PLDL1KEEP, [X, #1024+64]
.endm

.macro KERNEL_F16_FINALIZE
	fadd	v0.2d, v0.2d, v1.2d
	fadd	v2.2d, v2.2d, v3.2d
	fadd	v4.2d, v4.2d, v5.2d
	fadd	v6.2d, v6.2d, v7.2d

	fadd	v0.2d, v0.2d, v2.2d
	fadd	v4.2d, v4.2d, v6.2d

	fadd	v0.2d, v0.2d, v4.2d
.endm

.macro KERNEL_FINALIZE
	faddp	SSQD, v0.2d
.endm

.macro INIT_S
	lsl	INC_X, INC_X, #3
.endm

.macro KERNEL_S1
	ldr	TMPF, [X]
	add	X, X, INC_X
	fcvtl	v16.2d, v16.2s
	fmla	v0.2d, v16.2d, v16.2d
.endm

/*******************************************************************************
* End of macro definitions
*******************************************************************************/

	PROLOGUE

	INIT

	cmp	N, xzr
	ble	nrm2_kernel_zero
	cmp	INC_X, xzr
	ble	nrm2_kernel_zero
	cmp	INC_X, #1
	bne	nrm2_kernel_S_BEGIN

nrm2_kernel_F_BEGIN:

	asr	I, N, #4
	cmp	I, xzr
	beq	nrm2_kernel_S_BEGIN

	.align 5
nrm2_kernel_F16:

	KERNEL_F16

	subs	I, I, #1
	bne	nrm2_kernel_F16

	KERNEL_F16_FINALIZE

nrm2_kernel_F1:

	ands	I, N, #15
	ble	nrm2_kernel_L999

nrm2_kernel_F10:

	KERNEL_F1

	subs    I, I, #1
	bne     nrm2_kernel_F10

	b	nrm2_kernel_L999

nrm2_kernel_S_BEGIN:

	INIT_S

	asr	I, N, #2
	cmp	I, xzr
	ble	nrm2_kernel_S1

nrm2_kernel_S4:

	KERNEL_S1
	KERNEL_S1
	KERNEL_S1
	KERNEL_S1

	subs	I, I, #1
	bne	nrm2_kernel_S4

nrm2_kernel_S1:

	ands	I, N, #3
	ble	nrm2_kernel_L999

nrm2_kernel_S10:

	KERNEL_S1

	subs    I, I, #1
	bne     nrm2_kernel_S10

nrm2_kernel_L999:
	KERNEL_FINALIZE
	fsqrt	SSQD, SSQD
	fcvt	SSQ, SSQD
	ret

nrm2_kernel_zero:
	fmov	SSQ, wzr

	ret

	EPILOGUE
