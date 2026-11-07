/*******************************************************************************
Copyright (c) 2015, The OpenBLAS Project
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

#define	N		x0
#define	X		x1
#define	INC_X		x2

#define I		x3

#if !defined(DOUBLE)
#define SSQF s0
#endif

#define SSQ		d0
#define SCALE		d1
#define REGZERO		d5
#define REGONE		d6

/*******************************************************************************
* Macro definitions
*******************************************************************************/

.macro KERNEL_F1
#if !defined(DOUBLE)
	ldr	s4, [X], #4
	fcvt	d4, s4
#else
	ldr	d4, [X], #8
#endif
	fcmp	d4, REGZERO
	beq	2f     /* KERNEL_F1_NEXT_\@ */
	fabs	d4, d4
	fcmp	SCALE, d4
	bge	1f     /* KERNEL_F1_SCALE_GE_X_\@ */
	fdiv	d2, SCALE, d4
	fmul	d2, d2, d2
	fmul	d3, SSQ, d2
	fadd	SSQ, REGONE, d3
	fmov	SCALE, d4
	b	2f     /* KERNEL_F1_NEXT_\@ */
1:                     /* KERNEL_F1_SCALE_GE_X_\@: */
	fdiv	d2, d4, SCALE
	fmla	SSQ, d2, v2.d[0]
2:                     /* KERNEL_F1_NEXT_\@: */
.endm

.macro KERNEL_S1
#if !defined(DOUBLE)
	ldr	s4, [X]
	fcvt	d4, s4
#else
	ldr	d4, [X]
#endif
	fcmp	d4, REGZERO
	beq	KERNEL_S1_NEXT
	fabs	d4, d4
	fcmp	SCALE, d4
	bge	KERNEL_S1_SCALE_GE_X
	fdiv	d2, SCALE, d4
	fmul	d2, d2, d2
	fmul	d3, SSQ, d2
	fadd	SSQ, REGONE, d3
	fmov	SCALE, d4
	b	KERNEL_S1_NEXT
KERNEL_S1_SCALE_GE_X:
	fdiv	d2, d4, SCALE
	fmla	SSQ, d2, v2.d[0]
KERNEL_S1_NEXT:
	add	X, X, INC_X
.endm

.macro KERNEL_F8
	KERNEL_F1
	KERNEL_F1
	KERNEL_F1
	KERNEL_F1
	KERNEL_F1
	KERNEL_F1
	KERNEL_F1
	KERNEL_F1
.endm

.macro INIT_S
#if !defined(DOUBLE)
	lsl	INC_X, INC_X, #2		// INC_X * SIZE
#else
	lsl	INC_X, INC_X, #3		// INC_X * SIZE
#endif
.endm

.macro INIT
	eor	v1.16b, v1.16b, v1.16b		// scale=0.0
	fmov	SSQ, #1.0
	fmov	REGONE, SSQ
	fmov	REGZERO, SCALE
.endm

/*******************************************************************************
* End of macro definitions
*******************************************************************************/

	PROLOGUE

	.align 5

	/* Preserve the source for rare exact range recovery. */
	mov	x4, X
	mov	x5, INC_X
	INIT

	cmp	N, #0
	ble	.Lnrm2_kernel_L999

	cmp	INC_X, #0
	beq	.Lnrm2_kernel_L999


	cmp	INC_X, #1
	bne	.Lnrm2_kernel_S_BEGIN

.Lnrm2_kernel_F_BEGIN:

	asr	I, N, #3				// I = N / 8
	cmp	I, xzr
	ble	.Lnrm2_kernel_F1

.Lnrm2_kernel_F8:

	KERNEL_F8

	subs	I, I, #1
	bne	.Lnrm2_kernel_F8

.Lnrm2_kernel_F1:

	ands	I, N, #7
	ble	.Lnrm2_kernel_L999


.Lnrm2_kernel_F10:

	KERNEL_F1

	subs	I, I, #1
	bne	.Lnrm2_kernel_F10

	b	.Lnrm2_kernel_L999

.Lnrm2_kernel_S_BEGIN:

	INIT_S

	mov	I, N

	.align 5

.Lnrm2_kernel_S10:

	KERNEL_S1

	subs	I, I, #1
	bne	.Lnrm2_kernel_S10


.Lnrm2_kernel_L999:
	/* With fewer than 2^64 real components, a norm near overflow must
	 * have scale >= 2^991 (double) or 2^95 (float). Exclude nonfinite
	 * scale so the existing special-value path remains separate.
	 */
	fmov	x6, SCALE
	ubfx	x6, x6, #52, #11
#if !defined(DOUBLE)
	sub	x6, x6, #1118
	cmp	x6, #929
#else
	sub	x6, x6, #2014
	cmp	x6, #33
#endif
	blo	.Lnrm2_range
	fsqrt	SSQ, SSQ
	fmul	SSQ, SCALE, SSQ
#if !defined(DOUBLE)
	fcvt SSQF, SSQ
#endif
	ret

.Lnrm2_range:
	mov	X, x4
	mov	INC_X, x5
	mov	w3, #1
#if defined(__APPLE__)
#if defined(DOUBLE)
	b	_openblas_dnrm2_range
#else
	b	_openblas_snrm2_range
#endif
#else
#if defined(DOUBLE)
	b	openblas_dnrm2_range
#else
	b	openblas_snrm2_range
#endif
#endif

	EPILOGUE

