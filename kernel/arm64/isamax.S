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

#define	N	x0	/* vector length */
#define	X	x1	/* X vector address */
#define	INC_X	x2	/* X stride */
#define INDEX	x3	/* index of max/min value */
#define Z	x4	/* vector index */
#define I	x5	/* loop variable */
#define X_COPY	x6	/* copy of X address */
#define MAXF_Z	x7

/*******************************************************************************
* Macro definitions
*******************************************************************************/

#define MAXF	s5
#define TMPF	s6
#define TMPVF	{v6.s}[0]
#define SZ	4

/******************************************************************************/

.macro INIT_F1
	ldr	MAXF, [X], #SZ
	mov	Z, #1
	mov	INDEX, Z
	fabs	MAXF, MAXF
.endm

.macro KERNEL_F1
	ldr	TMPF, [X], #SZ
	add	Z, Z, #1
	fabs	TMPF, TMPF
	fcmp	TMPF, MAXF
	fcsel	MAXF, MAXF, TMPF, le
	csel	INDEX, INDEX, Z, le
.endm

.macro INIT_F4
	ld1	{v0.4s}, [X], #16
	fabs	v0.4s, v0.4s
	fmaxv	MAXF, v0.4s
	mov	Z, #5
	mov	MAXF_Z, #1
.endm

.macro KERNEL_F4
	ld1	{v0.4s}, [X], #16
	fabs	v0.4s, v0.4s
	fmaxv	TMPF, v0.4s
	PRFM	PLDL1KEEP, [X, #512]
	fcmp	TMPF, MAXF
	fcsel	MAXF, MAXF, TMPF, le
	csel	MAXF_Z, MAXF_Z, Z, le
	add	Z, Z, #4
.endm


.macro KERNEL_F4_FINALIZE
	mov	INDEX, MAXF_Z
	sub	MAXF_Z, MAXF_Z, #1
	lsl	MAXF_Z, MAXF_Z, #2
	add	X_COPY, X_COPY, MAXF_Z
	ldr	TMPF, [X_COPY], #SZ
	fabs	TMPF, TMPF
	fcmp	TMPF, MAXF
	beq	KERNEL_F4_FINALIZE_DONE
	add	INDEX, INDEX, #1
	ldr	TMPF, [X_COPY], #SZ
	fabs	TMPF, TMPF
	fcmp	TMPF, MAXF
	beq	KERNEL_F4_FINALIZE_DONE
	add	INDEX, INDEX, #1
	ldr	TMPF, [X_COPY], #SZ
	fabs	TMPF, TMPF
	fcmp	TMPF, MAXF
	beq	KERNEL_F4_FINALIZE_DONE
	add	INDEX, INDEX, #1
KERNEL_F4_FINALIZE_DONE:
.endm


.macro INIT_S
	lsl	INC_X, INC_X, #2
	ld1	TMPVF, [X], INC_X
	mov	Z, #1
	mov	INDEX, Z
	fabs	MAXF, TMPF
.endm

.macro KERNEL_S1
	ld1	TMPVF, [X], INC_X
	add	Z, Z, #1
	fabs	TMPF, TMPF
	fcmp	TMPF, MAXF
	fcsel	MAXF, MAXF, TMPF, le
	csel	INDEX, INDEX, Z, le
.endm

/*******************************************************************************
* End of macro definitions
*******************************************************************************/

	PROLOGUE

	cmp	N, xzr
	ble	iamax_kernel_zero
	cmp	INC_X, xzr
	ble	iamax_kernel_zero

	PRFM	PLDL1KEEP, [X]
	mov	X_COPY, X

	cmp	INC_X, #1
	bne	iamax_kernel_S_BEGIN

iamax_kernel_F_BEGIN:
	asr	I, N, #2
	cmp	I, xzr
	beq	iamax_kernel_F1_INIT

	INIT_F4
	subs	I, I, #1
	beq	iamax_kernel_F4_FINALIZE

iamax_kernel_F4:
	KERNEL_F4
	subs	I, I, #1
	bne	iamax_kernel_F4

iamax_kernel_F4_FINALIZE:
	KERNEL_F4_FINALIZE

iamax_kernel_F1:
	ands	I, N, #3
	ble	iamax_kernel_L999

iamax_kernel_F10:
	KERNEL_F1
	subs	I, I, #1
	bne	iamax_kernel_F10
	b	iamax_kernel_L999

iamax_kernel_F1_INIT:
	INIT_F1
	subs	N, N, #1
	b	iamax_kernel_F1

iamax_kernel_S_BEGIN:
	INIT_S

	subs	N, N, #1
	ble	iamax_kernel_L999

	asr	I, N, #2
	cmp	I, xzr
	ble	iamax_kernel_S1

iamax_kernel_S4:
	KERNEL_S1
	KERNEL_S1
	KERNEL_S1
	KERNEL_S1

	subs	I, I, #1
	bne	iamax_kernel_S4

iamax_kernel_S1:
	ands	I, N, #3
	ble	iamax_kernel_L999

iamax_kernel_S10:
	KERNEL_S1
	subs	I, I, #1
	bne	iamax_kernel_S10

iamax_kernel_L999:
	mov	x0, INDEX
	ret

iamax_kernel_zero:
	mov	x0, xzr
	ret

	EPILOGUE
