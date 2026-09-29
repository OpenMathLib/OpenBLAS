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

/*                   X0          X1          X2          s0         X3        x4       x5           x6            x7 */
/*int CNAME(BLASLONG bm,BLASLONG bn,BLASLONG bk,FLOAT alpha0,FLOAT* ba,FLOAT* bb,FLOAT* C,BLASLONG ldc, BLASLONG offset) */

#define origM		x0
#define origN		x1
#define origK		x2
#define origPA		x3
#define origPB		x4
#define pC		x5
#define LDC		x6
#define offset		x7
#define counterL	x8
#define counterI	x9
#define counterJ	x10
#define pB		x11
#define pCRow0		x12
#define pCRow1		x13
#define pCRow2		x14
#define pA		x15
#define temp		x16
#define tempOffset	x19
#define tempK		x20

/*
 * Cortex-A72 DTRMM 6x8 microkernel (mr=6, nr=8).
 * Derived from dgemm_kernel_6x8_cortexa72.S with TRMM offset/K panel
 * stepping (same scheme as dtrmm_kernel_8x4.S). SAVE overwrites C.
 */
#define A_PRE_SIZE	2560
#define B_PRE_SIZE	448
#define C_PRE_SIZE	128

#define alpha0		d7
#define alphaV0		v7.d[0]

.macro INIT6x8
	fmov		d8, xzr
	fmov		d9, xzr
	fmov		d10, xzr
	fmov		d11, xzr
	fmov		d12, xzr
	fmov		d13, xzr
	fmov		d14, xzr
	fmov		d15, xzr
	fmov		d16, xzr
	fmov		d17, xzr
	fmov		d18, xzr
	fmov		d19, xzr
	fmov		d20, xzr
	fmov		d21, xzr
	fmov		d22, xzr
	fmov		d23, xzr
	fmov		d24, xzr
	fmov		d25, xzr
	fmov		d26, xzr
	fmov		d27, xzr
	fmov		d28, xzr
	fmov		d29, xzr
	fmov		d30, xzr
	fmov		d31, xzr
.endm

.macro KERNEL6x8_SUB
	ld1	{v0.2d, v1.2d, v2.2d}, [pA], #48
	ld1	{v3.2d, v4.2d, v5.2d, v6.2d}, [pB], #64
	prfm	PLDL1KEEP, [pA, #A_PRE_SIZE]
	prfm	PLDL1KEEP, [pB, #B_PRE_SIZE]
	fmla	v8.2d, v0.2d, v3.d[0]
	fmla	v9.2d, v1.2d, v3.d[0]
	fmla	v10.2d, v2.2d, v3.d[0]
	fmla	v11.2d, v0.2d, v3.d[1]
	fmla	v12.2d, v1.2d, v3.d[1]
	fmla	v13.2d, v2.2d, v3.d[1]
	fmla	v14.2d, v0.2d, v4.d[0]
	fmla	v15.2d, v1.2d, v4.d[0]
	fmla	v16.2d, v2.2d, v4.d[0]
	fmla	v17.2d, v0.2d, v4.d[1]
	fmla	v18.2d, v1.2d, v4.d[1]
	fmla	v19.2d, v2.2d, v4.d[1]
	fmla	v20.2d, v0.2d, v5.d[0]
	fmla	v21.2d, v1.2d, v5.d[0]
	fmla	v22.2d, v2.2d, v5.d[0]
	fmla	v23.2d, v0.2d, v5.d[1]
	fmla	v24.2d, v1.2d, v5.d[1]
	fmla	v25.2d, v2.2d, v5.d[1]
	fmla	v26.2d, v0.2d, v6.d[0]
	fmla	v27.2d, v1.2d, v6.d[0]
	fmla	v28.2d, v2.2d, v6.d[0]
	fmla	v29.2d, v0.2d, v6.d[1]
	fmla	v30.2d, v1.2d, v6.d[1]
	fmla	v31.2d, v2.2d, v6.d[1]
.endm

.macro SAVE6x8
	mov	pCRow1, pCRow0
	prfm	PLDL2KEEP, [pCRow1, #C_PRE_SIZE]
	fmul	v0.2d, v8.2d, alphaV0
	fmul	v1.2d, v9.2d, alphaV0
	fmul	v2.2d, v10.2d, alphaV0
	st1	{v0.2d, v1.2d, v2.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	prfm	PLDL2KEEP, [pCRow1, #C_PRE_SIZE]
	fmul	v0.2d, v11.2d, alphaV0
	fmul	v1.2d, v12.2d, alphaV0
	fmul	v2.2d, v13.2d, alphaV0
	st1	{v0.2d, v1.2d, v2.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	prfm	PLDL2KEEP, [pCRow1, #C_PRE_SIZE]
	fmul	v0.2d, v14.2d, alphaV0
	fmul	v1.2d, v15.2d, alphaV0
	fmul	v2.2d, v16.2d, alphaV0
	st1	{v0.2d, v1.2d, v2.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	prfm	PLDL2KEEP, [pCRow1, #C_PRE_SIZE]
	fmul	v0.2d, v17.2d, alphaV0
	fmul	v1.2d, v18.2d, alphaV0
	fmul	v2.2d, v19.2d, alphaV0
	st1	{v0.2d, v1.2d, v2.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	prfm	PLDL2KEEP, [pCRow1, #C_PRE_SIZE]
	fmul	v0.2d, v20.2d, alphaV0
	fmul	v1.2d, v21.2d, alphaV0
	fmul	v2.2d, v22.2d, alphaV0
	st1	{v0.2d, v1.2d, v2.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	prfm	PLDL2KEEP, [pCRow1, #C_PRE_SIZE]
	fmul	v0.2d, v23.2d, alphaV0
	fmul	v1.2d, v24.2d, alphaV0
	fmul	v2.2d, v25.2d, alphaV0
	st1	{v0.2d, v1.2d, v2.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	prfm	PLDL2KEEP, [pCRow1, #C_PRE_SIZE]
	fmul	v0.2d, v26.2d, alphaV0
	fmul	v1.2d, v27.2d, alphaV0
	fmul	v2.2d, v28.2d, alphaV0
	st1	{v0.2d, v1.2d, v2.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	prfm	PLDL2KEEP, [pCRow1, #C_PRE_SIZE]
	fmul	v0.2d, v29.2d, alphaV0
	fmul	v1.2d, v30.2d, alphaV0
	fmul	v2.2d, v31.2d, alphaV0
	st1	{v0.2d, v1.2d, v2.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	add	pCRow0, pCRow0, #48
.endm

.macro INIT4x8
	fmov	d16, xzr
	fmov	d17, xzr
	fmov	d18, xzr
	fmov	d19, xzr
	fmov	d20, xzr
	fmov	d21, xzr
	fmov	d22, xzr
	fmov	d23, xzr
	fmov	d24, xzr
	fmov	d25, xzr
	fmov	d26, xzr
	fmov	d27, xzr
	fmov	d28, xzr
	fmov	d29, xzr
	fmov	d30, xzr
	fmov	d31, xzr
.endm

.macro KERNEL4x8_SUB
	ld1	{v0.2d, v1.2d}, [pA], #32
	ld1	{v3.2d, v4.2d, v5.2d, v6.2d}, [pB], #64
	fmla	v16.2d, v0.2d, v3.d[0]
	fmla	v17.2d, v1.2d, v3.d[0]
	fmla	v18.2d, v0.2d, v3.d[1]
	fmla	v19.2d, v1.2d, v3.d[1]
	fmla	v20.2d, v0.2d, v4.d[0]
	fmla	v21.2d, v1.2d, v4.d[0]
	fmla	v22.2d, v0.2d, v4.d[1]
	fmla	v23.2d, v1.2d, v4.d[1]
	fmla	v24.2d, v0.2d, v5.d[0]
	fmla	v25.2d, v1.2d, v5.d[0]
	fmla	v26.2d, v0.2d, v5.d[1]
	fmla	v27.2d, v1.2d, v5.d[1]
	fmla	v28.2d, v0.2d, v6.d[0]
	fmla	v29.2d, v1.2d, v6.d[0]
	fmla	v30.2d, v0.2d, v6.d[1]
	fmla	v31.2d, v1.2d, v6.d[1]
.endm

.macro SAVE4x8
	mov	pCRow1, pCRow0
	fmul	v0.2d, v16.2d, alphaV0
	fmul	v1.2d, v17.2d, alphaV0
	st1	{v0.2d, v1.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v18.2d, alphaV0
	fmul	v1.2d, v19.2d, alphaV0
	st1	{v0.2d, v1.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v20.2d, alphaV0
	fmul	v1.2d, v21.2d, alphaV0
	st1	{v0.2d, v1.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v22.2d, alphaV0
	fmul	v1.2d, v23.2d, alphaV0
	st1	{v0.2d, v1.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v24.2d, alphaV0
	fmul	v1.2d, v25.2d, alphaV0
	st1	{v0.2d, v1.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v26.2d, alphaV0
	fmul	v1.2d, v27.2d, alphaV0
	st1	{v0.2d, v1.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v28.2d, alphaV0
	fmul	v1.2d, v29.2d, alphaV0
	st1	{v0.2d, v1.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v30.2d, alphaV0
	fmul	v1.2d, v31.2d, alphaV0
	st1	{v0.2d, v1.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	add	pCRow0, pCRow0, #32
.endm

.macro INIT2x8
	fmov	d16, xzr
	fmov	d18, xzr
	fmov	d20, xzr
	fmov	d22, xzr
	fmov	d24, xzr
	fmov	d26, xzr
	fmov	d28, xzr
	fmov	d30, xzr
.endm

.macro KERNEL2x8_SUB
	ld1	{v0.2d}, [pA], #16
	ld1	{v3.2d, v4.2d, v5.2d, v6.2d}, [pB], #64
	fmla	v16.2d, v0.2d, v3.d[0]
	fmla	v18.2d, v0.2d, v3.d[1]
	fmla	v20.2d, v0.2d, v4.d[0]
	fmla	v22.2d, v0.2d, v4.d[1]
	fmla	v24.2d, v0.2d, v5.d[0]
	fmla	v26.2d, v0.2d, v5.d[1]
	fmla	v28.2d, v0.2d, v6.d[0]
	fmla	v30.2d, v0.2d, v6.d[1]
.endm

.macro SAVE2x8
	mov	pCRow1, pCRow0
	fmul	v0.2d, v16.2d, alphaV0
	str	q0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v18.2d, alphaV0
	str	q0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v20.2d, alphaV0
	str	q0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v22.2d, alphaV0
	str	q0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v24.2d, alphaV0
	str	q0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v26.2d, alphaV0
	str	q0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v28.2d, alphaV0
	str	q0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v30.2d, alphaV0
	str	q0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	add	pCRow0, pCRow0, #16
.endm

.macro INIT1x8
	fmov	d16, xzr
	fmov	d17, xzr
	fmov	d18, xzr
	fmov	d19, xzr
	fmov	d20, xzr
	fmov	d21, xzr
	fmov	d22, xzr
	fmov	d23, xzr
.endm

.macro KERNEL1x8_SUB
	ldr	d0, [pA], #8
	ld1	{v3.2d, v4.2d, v5.2d, v6.2d}, [pB], #64
	fmadd	d16, d0, d3, d16
	fmul	d1, d0, v3.d[1]
	fadd	d17, d17, d1
	fmadd	d18, d0, d4, d18
	fmul	d1, d0, v4.d[1]
	fadd	d19, d19, d1
	fmadd	d20, d0, d5, d20
	fmul	d1, d0, v5.d[1]
	fadd	d21, d21, d1
	fmadd	d22, d0, d6, d22
	fmul	d1, d0, v6.d[1]
	fadd	d23, d23, d1
.endm

.macro SAVE1x8
	mov	pCRow1, pCRow0
	fmul	d0, d16, alpha0
	str	d0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	d0, d17, alpha0
	str	d0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	d0, d18, alpha0
	str	d0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	d0, d19, alpha0
	str	d0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	d0, d20, alpha0
	str	d0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	d0, d21, alpha0
	str	d0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	d0, d22, alpha0
	str	d0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	d0, d23, alpha0
	str	d0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	add	pCRow0, pCRow0, #8
.endm

.macro INIT6x4
	fmov	d8, xzr
	fmov	d9, xzr
	fmov	d10, xzr
	fmov	d11, xzr
	fmov	d12, xzr
	fmov	d13, xzr
	fmov	d14, xzr
	fmov	d15, xzr
	fmov	d16, xzr
	fmov	d17, xzr
	fmov	d18, xzr
	fmov	d19, xzr
.endm

.macro KERNEL6x4_SUB
	ld1	{v0.2d, v1.2d, v2.2d}, [pA], #48
	ld1	{v3.2d, v4.2d}, [pB], #32
	fmla	v8.2d, v0.2d, v3.d[0]
	fmla	v9.2d, v1.2d, v3.d[0]
	fmla	v10.2d, v2.2d, v3.d[0]
	fmla	v11.2d, v0.2d, v3.d[1]
	fmla	v12.2d, v1.2d, v3.d[1]
	fmla	v13.2d, v2.2d, v3.d[1]
	fmla	v14.2d, v0.2d, v4.d[0]
	fmla	v15.2d, v1.2d, v4.d[0]
	fmla	v16.2d, v2.2d, v4.d[0]
	fmla	v17.2d, v0.2d, v4.d[1]
	fmla	v18.2d, v1.2d, v4.d[1]
	fmla	v19.2d, v2.2d, v4.d[1]
.endm

.macro SAVE6x4
	mov	pCRow1, pCRow0
	fmul	v0.2d, v8.2d, alphaV0
	fmul	v1.2d, v9.2d, alphaV0
	fmul	v2.2d, v10.2d, alphaV0
	st1	{v0.2d, v1.2d, v2.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v11.2d, alphaV0
	fmul	v1.2d, v12.2d, alphaV0
	fmul	v2.2d, v13.2d, alphaV0
	st1	{v0.2d, v1.2d, v2.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v14.2d, alphaV0
	fmul	v1.2d, v15.2d, alphaV0
	fmul	v2.2d, v16.2d, alphaV0
	st1	{v0.2d, v1.2d, v2.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v17.2d, alphaV0
	fmul	v1.2d, v18.2d, alphaV0
	fmul	v2.2d, v19.2d, alphaV0
	st1	{v0.2d, v1.2d, v2.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	add	pCRow0, pCRow0, #48
.endm

.macro INIT6x2
	fmov	d8, xzr
	fmov	d9, xzr
	fmov	d10, xzr
	fmov	d11, xzr
	fmov	d12, xzr
	fmov	d13, xzr
.endm

.macro KERNEL6x2_SUB
	ld1	{v0.2d, v1.2d, v2.2d}, [pA], #48
	ld1	{v3.2d}, [pB], #16
	fmla	v8.2d, v0.2d, v3.d[0]
	fmla	v9.2d, v1.2d, v3.d[0]
	fmla	v10.2d, v2.2d, v3.d[0]
	fmla	v11.2d, v0.2d, v3.d[1]
	fmla	v12.2d, v1.2d, v3.d[1]
	fmla	v13.2d, v2.2d, v3.d[1]
.endm

.macro SAVE6x2
	mov	pCRow1, pCRow0
	fmul	v0.2d, v8.2d, alphaV0
	fmul	v1.2d, v9.2d, alphaV0
	fmul	v2.2d, v10.2d, alphaV0
	st1	{v0.2d, v1.2d, v2.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v11.2d, alphaV0
	fmul	v1.2d, v12.2d, alphaV0
	fmul	v2.2d, v13.2d, alphaV0
	st1	{v0.2d, v1.2d, v2.2d}, [pCRow1]
	add	pCRow0, pCRow0, #48
.endm

.macro INIT6x1
	fmov	d8, xzr
	fmov	d9, xzr
	fmov	d10, xzr
.endm

.macro KERNEL6x1_SUB
	ld1	{v0.2d, v1.2d, v2.2d}, [pA], #48
	ldr	d3, [pB], #8
	fmla	v8.2d, v0.2d, v3.d[0]
	fmla	v9.2d, v1.2d, v3.d[0]
	fmla	v10.2d, v2.2d, v3.d[0]
.endm

.macro SAVE6x1
	fmul	v0.2d, v8.2d, alphaV0
	fmul	v1.2d, v9.2d, alphaV0
	fmul	v2.2d, v10.2d, alphaV0
	st1	{v0.2d, v1.2d, v2.2d}, [pCRow0]
	add	pCRow0, pCRow0, #48
.endm


.macro INIT4x4
	fmov	d16, xzr
	fmov	d17, xzr
	fmov	d18, xzr
	fmov	d19, xzr
	fmov	d20, xzr
	fmov	d21, xzr
	fmov	d22, xzr
	fmov	d23, xzr
.endm
.macro KERNEL4x4_SUB
	ld1	{v0.2d, v1.2d}, [pA], #32
	ld1	{v3.2d, v4.2d}, [pB], #32
	fmla	v16.2d, v0.2d, v3.d[0]
	fmla	v17.2d, v1.2d, v3.d[0]
	fmla	v18.2d, v0.2d, v3.d[1]
	fmla	v19.2d, v1.2d, v3.d[1]
	fmla	v20.2d, v0.2d, v4.d[0]
	fmla	v21.2d, v1.2d, v4.d[0]
	fmla	v22.2d, v0.2d, v4.d[1]
	fmla	v23.2d, v1.2d, v4.d[1]
.endm
.macro SAVE4x4
	mov	pCRow1, pCRow0
	fmul	v0.2d, v16.2d, alphaV0
	fmul	v1.2d, v17.2d, alphaV0
	st1	{v0.2d, v1.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v18.2d, alphaV0
	fmul	v1.2d, v19.2d, alphaV0
	st1	{v0.2d, v1.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v20.2d, alphaV0
	fmul	v1.2d, v21.2d, alphaV0
	st1	{v0.2d, v1.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v22.2d, alphaV0
	fmul	v1.2d, v23.2d, alphaV0
	st1	{v0.2d, v1.2d}, [pCRow1]
	add	pCRow0, pCRow0, #32
.endm

.macro INIT2x4
	fmov	d16, xzr
	fmov	d18, xzr
	fmov	d20, xzr
	fmov	d22, xzr
.endm
.macro KERNEL2x4_SUB
	ld1	{v0.2d}, [pA], #16
	ld1	{v3.2d, v4.2d}, [pB], #32
	fmla	v16.2d, v0.2d, v3.d[0]
	fmla	v18.2d, v0.2d, v3.d[1]
	fmla	v20.2d, v0.2d, v4.d[0]
	fmla	v22.2d, v0.2d, v4.d[1]
.endm
.macro SAVE2x4
	mov	pCRow1, pCRow0
	fmul	v0.2d, v16.2d, alphaV0
	str	q0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v18.2d, alphaV0
	str	q0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v20.2d, alphaV0
	str	q0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v22.2d, alphaV0
	str	q0, [pCRow1]
	add	pCRow0, pCRow0, #16
.endm

.macro INIT1x4
	fmov	d16, xzr
	fmov	d17, xzr
	fmov	d18, xzr
	fmov	d19, xzr
.endm
.macro KERNEL1x4_SUB
	ldr	d0, [pA], #8
	ld1	{v3.2d, v4.2d}, [pB], #32
	fmadd	d16, d0, d3, d16
	fmul	d1, d0, v3.d[1]
	fadd	d17, d17, d1
	fmadd	d18, d0, d4, d18
	fmul	d1, d0, v4.d[1]
	fadd	d19, d19, d1
.endm
.macro SAVE1x4
	mov	pCRow1, pCRow0
	fmul	d0, d16, alpha0
	str	d0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	d0, d17, alpha0
	str	d0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	d0, d18, alpha0
	str	d0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	d0, d19, alpha0
	str	d0, [pCRow1]
	add	pCRow0, pCRow0, #8
.endm

.macro INIT4x2
	fmov	d16, xzr
	fmov	d17, xzr
	fmov	d18, xzr
	fmov	d19, xzr
.endm
.macro KERNEL4x2_SUB
	ld1	{v0.2d, v1.2d}, [pA], #32
	ld1	{v3.2d}, [pB], #16
	fmla	v16.2d, v0.2d, v3.d[0]
	fmla	v17.2d, v1.2d, v3.d[0]
	fmla	v18.2d, v0.2d, v3.d[1]
	fmla	v19.2d, v1.2d, v3.d[1]
.endm
.macro SAVE4x2
	mov	pCRow1, pCRow0
	fmul	v0.2d, v16.2d, alphaV0
	fmul	v1.2d, v17.2d, alphaV0
	st1	{v0.2d, v1.2d}, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v18.2d, alphaV0
	fmul	v1.2d, v19.2d, alphaV0
	st1	{v0.2d, v1.2d}, [pCRow1]
	add	pCRow0, pCRow0, #32
.endm

.macro INIT2x2
	fmov	d16, xzr
	fmov	d18, xzr
.endm
.macro KERNEL2x2_SUB
	ld1	{v0.2d}, [pA], #16
	ld1	{v3.2d}, [pB], #16
	fmla	v16.2d, v0.2d, v3.d[0]
	fmla	v18.2d, v0.2d, v3.d[1]
.endm
.macro SAVE2x2
	mov	pCRow1, pCRow0
	fmul	v0.2d, v16.2d, alphaV0
	str	q0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	v0.2d, v18.2d, alphaV0
	str	q0, [pCRow1]
	add	pCRow0, pCRow0, #16
.endm

.macro INIT1x2
	fmov	d16, xzr
	fmov	d17, xzr
.endm
.macro KERNEL1x2_SUB
	ldr	d0, [pA], #8
	ld1	{v3.2d}, [pB], #16
	fmadd	d16, d0, d3, d16
	fmul	d1, d0, v3.d[1]
	fadd	d17, d17, d1
.endm
.macro SAVE1x2
	mov	pCRow1, pCRow0
	fmul	d0, d16, alpha0
	str	d0, [pCRow1]
	add	pCRow1, pCRow1, LDC
	fmul	d0, d17, alpha0
	str	d0, [pCRow1]
	add	pCRow0, pCRow0, #8
.endm

.macro INIT4x1
	fmov	d16, xzr
	fmov	d17, xzr
.endm
.macro KERNEL4x1_SUB
	ld1	{v0.2d, v1.2d}, [pA], #32
	ldr	d3, [pB], #8
	fmla	v16.2d, v0.2d, v3.d[0]
	fmla	v17.2d, v1.2d, v3.d[0]
.endm
.macro SAVE4x1
	fmul	v0.2d, v16.2d, alphaV0
	fmul	v1.2d, v17.2d, alphaV0
	st1	{v0.2d, v1.2d}, [pCRow0]
	add	pCRow0, pCRow0, #32
.endm

.macro INIT2x1
	fmov	d16, xzr
.endm
.macro KERNEL2x1_SUB
	ld1	{v0.2d}, [pA], #16
	ldr	d3, [pB], #8
	fmla	v16.2d, v0.2d, v3.d[0]
.endm
.macro SAVE2x1
	fmul	v0.2d, v16.2d, alphaV0
	str	q0, [pCRow0]
	add	pCRow0, pCRow0, #16
.endm

.macro INIT1x1
	fmov	d16, xzr
.endm
.macro KERNEL1x1_SUB
	ldr	d0, [pA], #8
	ldr	d3, [pB], #8
	fmadd	d16, d0, d3, d16
.endm
.macro SAVE1x1
	fmul	d0, d16, alpha0
	str	d0, [pCRow0]
	add	pCRow0, pCRow0, #8
.endm

	PROLOGUE
	.align 5
	add	sp, sp, #-(11 * 16)
	stp	d8, d9, [sp, #(0 * 16)]
	stp	d10, d11, [sp, #(1 * 16)]
	stp	d12, d13, [sp, #(2 * 16)]
	stp	d14, d15, [sp, #(3 * 16)]
	stp	d16, d17, [sp, #(4 * 16)]
	stp	x18, x19, [sp, #(5 * 16)]
	stp	x20, x21, [sp, #(6 * 16)]
	stp	x22, x23, [sp, #(7 * 16)]
	stp	x24, x25, [sp, #(8 * 16)]
	stp	x26, x27, [sp, #(9 * 16)]
	str	x28, [sp, #(10 * 16)]

	fmov	alpha0, d0
	lsl	LDC, LDC, #3

#if !defined(LEFT)
	neg	tempOffset, offset
#endif

	mov	counterJ, origN
	asr	counterJ, counterJ, #3
	cbz	counterJ, .Ldtrmm6x8_L4_BEGIN

.Ldtrmm6x8_L8_BEGIN:
	mov	pCRow0, pC
	add	pC, pC, LDC, lsl #3
#if defined(LEFT)
	mov	tempOffset, offset
#endif
	mov	pA, origPA

	mov	temp, #6
	udiv	counterI, origM, temp
	cbz	counterI, .Ldtrmm6x8_L8_Mrem

.Ldtrmm6x8_L8_M6:
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
	#else
	mov	pB, origPB
	mov	temp, #48
	mul	temp, tempOffset, temp
	add	pA, pA, temp
	mov	temp, #64
	mul	temp, tempOffset, temp
	add	pB, pB, temp
	#endif
	#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
	#elif defined(LEFT)
	add	tempK, tempOffset, #6
	#else
	add	tempK, tempOffset, #8
	#endif

	INIT6x8
	mov	counterL, tempK
	cbz	counterL, .Ldtrmm6x8_L8_M6s
.Ldtrmm6x8_L8_M6l:
	KERNEL6x8_SUB
	subs	counterL, counterL, #1
	bgt	.Ldtrmm6x8_L8_M6l
.Ldtrmm6x8_L8_M6s:
	SAVE6x8
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
	#if defined(LEFT)
	sub	tempK, tempK, #6
	#else
	sub	tempK, tempK, #8
	#endif
	mov	temp, #48
	mul	temp, tempK, temp
	add	pA, pA, temp
	mov	temp, #64
	mul	temp, tempK, temp
	add	pB, pB, temp
	#endif
	#if defined(LEFT)
	add	tempOffset, tempOffset, #6
	#endif

	subs	counterI, counterI, #1
	bgt	.Ldtrmm6x8_L8_M6

.Ldtrmm6x8_L8_Mrem:
	mov	temp, #6
	udiv	counterI, origM, temp
	msub	counterI, counterI, temp, origM
	cmp	counterI, #4
	blt	.Ldtrmm6x8_L8_M2c
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
	#else
	mov	pB, origPB
	mov	temp, #32
	mul	temp, tempOffset, temp
	add	pA, pA, temp
	mov	temp, #64
	mul	temp, tempOffset, temp
	add	pB, pB, temp
	#endif
	#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
	#elif defined(LEFT)
	add	tempK, tempOffset, #4
	#else
	add	tempK, tempOffset, #8
	#endif

	INIT4x8
	mov	counterL, tempK
	cbz	counterL, .Ldtrmm6x8_L8_M4s
.Ldtrmm6x8_L8_M4l:
	KERNEL4x8_SUB
	subs	counterL, counterL, #1
	bgt	.Ldtrmm6x8_L8_M4l
.Ldtrmm6x8_L8_M4s:
	SAVE4x8
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
	#if defined(LEFT)
	sub	tempK, tempK, #4
	#else
	sub	tempK, tempK, #8
	#endif
	mov	temp, #32
	mul	temp, tempK, temp
	add	pA, pA, temp
	mov	temp, #64
	mul	temp, tempK, temp
	add	pB, pB, temp
	#endif
	#if defined(LEFT)
	add	tempOffset, tempOffset, #4
	#endif

	sub	counterI, counterI, #4
.Ldtrmm6x8_L8_M2c:
	cmp	counterI, #2
	blt	.Ldtrmm6x8_L8_M1c
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
	#else
	mov	pB, origPB
	mov	temp, #16
	mul	temp, tempOffset, temp
	add	pA, pA, temp
	mov	temp, #64
	mul	temp, tempOffset, temp
	add	pB, pB, temp
	#endif
	#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
	#elif defined(LEFT)
	add	tempK, tempOffset, #2
	#else
	add	tempK, tempOffset, #8
	#endif

	INIT2x8
	mov	counterL, tempK
	cbz	counterL, .Ldtrmm6x8_L8_M2s
.Ldtrmm6x8_L8_M2l:
	KERNEL2x8_SUB
	subs	counterL, counterL, #1
	bgt	.Ldtrmm6x8_L8_M2l
.Ldtrmm6x8_L8_M2s:
	SAVE2x8
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
	#if defined(LEFT)
	sub	tempK, tempK, #2
	#else
	sub	tempK, tempK, #8
	#endif
	mov	temp, #16
	mul	temp, tempK, temp
	add	pA, pA, temp
	mov	temp, #64
	mul	temp, tempK, temp
	add	pB, pB, temp
	#endif
	#if defined(LEFT)
	add	tempOffset, tempOffset, #2
	#endif

	sub	counterI, counterI, #2
.Ldtrmm6x8_L8_M1c:
	cbz	counterI, .Ldtrmm6x8_L8_END
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
	#else
	mov	pB, origPB
	mov	temp, #8
	mul	temp, tempOffset, temp
	add	pA, pA, temp
	mov	temp, #64
	mul	temp, tempOffset, temp
	add	pB, pB, temp
	#endif
	#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
	#elif defined(LEFT)
	add	tempK, tempOffset, #1
	#else
	add	tempK, tempOffset, #8
	#endif

	INIT1x8
	mov	counterL, tempK
	cbz	counterL, .Ldtrmm6x8_L8_M1s
.Ldtrmm6x8_L8_M1l:
	KERNEL1x8_SUB
	subs	counterL, counterL, #1
	bgt	.Ldtrmm6x8_L8_M1l
.Ldtrmm6x8_L8_M1s:
	SAVE1x8
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
	#if defined(LEFT)
	sub	tempK, tempK, #1
	#else
	sub	tempK, tempK, #8
	#endif
	mov	temp, #8
	mul	temp, tempK, temp
	add	pA, pA, temp
	mov	temp, #64
	mul	temp, tempK, temp
	add	pB, pB, temp
	#endif
	#if defined(LEFT)
	add	tempOffset, tempOffset, #1
	#endif

.Ldtrmm6x8_L8_END:
	lsl	temp, origK, #6
	add	origPB, origPB, temp
#if !defined(LEFT)
	add	tempOffset, tempOffset, #8
#endif
	subs	counterJ, counterJ, #1
	bgt	.Ldtrmm6x8_L8_BEGIN

.Ldtrmm6x8_L4_BEGIN:
	tst	origN, #7
	beq	.Ldtrmm6x8_L999
	tst	origN, #4
	beq	.Ldtrmm6x8_L2_BEGIN

	mov	pCRow0, pC
	add	pC, pC, LDC, lsl #2
#if defined(LEFT)
	mov	tempOffset, offset
#endif
	mov	pA, origPA

	mov	temp, #6
	udiv	counterI, origM, temp
	cbz	counterI, .Ldtrmm6x8_L4_Mrem
.Ldtrmm6x8_L4_M6:
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
	#else
	mov	pB, origPB
	mov	temp, #48
	mul	temp, tempOffset, temp
	add	pA, pA, temp
	mov	temp, #32
	mul	temp, tempOffset, temp
	add	pB, pB, temp
	#endif
	#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
	#elif defined(LEFT)
	add	tempK, tempOffset, #6
	#else
	add	tempK, tempOffset, #4
	#endif

	INIT6x4
	mov	counterL, tempK
	cbz	counterL, .Ldtrmm6x8_L4_M6s
.Ldtrmm6x8_L4_M6l:
	KERNEL6x4_SUB
	subs	counterL, counterL, #1
	bgt	.Ldtrmm6x8_L4_M6l
.Ldtrmm6x8_L4_M6s:
	SAVE6x4
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
	#if defined(LEFT)
	sub	tempK, tempK, #6
	#else
	sub	tempK, tempK, #4
	#endif
	mov	temp, #48
	mul	temp, tempK, temp
	add	pA, pA, temp
	mov	temp, #32
	mul	temp, tempK, temp
	add	pB, pB, temp
	#endif
	#if defined(LEFT)
	add	tempOffset, tempOffset, #6
	#endif

	subs	counterI, counterI, #1
	bgt	.Ldtrmm6x8_L4_M6
.Ldtrmm6x8_L4_Mrem:
	mov	temp, #6
	udiv	counterI, origM, temp
	msub	counterI, counterI, temp, origM
	cmp	counterI, #4
	blt	.Ldtrmm6x8_L4_M2c
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
	#else
	mov	pB, origPB
	mov	temp, #32
	mul	temp, tempOffset, temp
	add	pA, pA, temp
	mov	temp, #32
	mul	temp, tempOffset, temp
	add	pB, pB, temp
	#endif
	#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
	#elif defined(LEFT)
	add	tempK, tempOffset, #4
	#else
	add	tempK, tempOffset, #4
	#endif

	INIT4x4
	mov	counterL, tempK
	cbz	counterL, .Ldtrmm6x8_L4_M4s
.Ldtrmm6x8_L4_M4l:
	KERNEL4x4_SUB
	subs	counterL, counterL, #1
	bgt	.Ldtrmm6x8_L4_M4l
.Ldtrmm6x8_L4_M4s:
	SAVE4x4
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
	#if defined(LEFT)
	sub	tempK, tempK, #4
	#else
	sub	tempK, tempK, #4
	#endif
	mov	temp, #32
	mul	temp, tempK, temp
	add	pA, pA, temp
	mov	temp, #32
	mul	temp, tempK, temp
	add	pB, pB, temp
	#endif
	#if defined(LEFT)
	add	tempOffset, tempOffset, #4
	#endif

	sub	counterI, counterI, #4
.Ldtrmm6x8_L4_M2c:
	cmp	counterI, #2
	blt	.Ldtrmm6x8_L4_M1c
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
	#else
	mov	pB, origPB
	mov	temp, #16
	mul	temp, tempOffset, temp
	add	pA, pA, temp
	mov	temp, #32
	mul	temp, tempOffset, temp
	add	pB, pB, temp
	#endif
	#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
	#elif defined(LEFT)
	add	tempK, tempOffset, #2
	#else
	add	tempK, tempOffset, #4
	#endif

	INIT2x4
	mov	counterL, tempK
	cbz	counterL, .Ldtrmm6x8_L4_M2s
.Ldtrmm6x8_L4_M2l:
	KERNEL2x4_SUB
	subs	counterL, counterL, #1
	bgt	.Ldtrmm6x8_L4_M2l
.Ldtrmm6x8_L4_M2s:
	SAVE2x4
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
	#if defined(LEFT)
	sub	tempK, tempK, #2
	#else
	sub	tempK, tempK, #4
	#endif
	mov	temp, #16
	mul	temp, tempK, temp
	add	pA, pA, temp
	mov	temp, #32
	mul	temp, tempK, temp
	add	pB, pB, temp
	#endif
	#if defined(LEFT)
	add	tempOffset, tempOffset, #2
	#endif

	sub	counterI, counterI, #2
.Ldtrmm6x8_L4_M1c:
	cbz	counterI, .Ldtrmm6x8_L4_END
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
	#else
	mov	pB, origPB
	mov	temp, #8
	mul	temp, tempOffset, temp
	add	pA, pA, temp
	mov	temp, #32
	mul	temp, tempOffset, temp
	add	pB, pB, temp
	#endif
	#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
	#elif defined(LEFT)
	add	tempK, tempOffset, #1
	#else
	add	tempK, tempOffset, #4
	#endif

	INIT1x4
	mov	counterL, tempK
	cbz	counterL, .Ldtrmm6x8_L4_M1s
.Ldtrmm6x8_L4_M1l:
	KERNEL1x4_SUB
	subs	counterL, counterL, #1
	bgt	.Ldtrmm6x8_L4_M1l
.Ldtrmm6x8_L4_M1s:
	SAVE1x4
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
	#if defined(LEFT)
	sub	tempK, tempK, #1
	#else
	sub	tempK, tempK, #4
	#endif
	mov	temp, #8
	mul	temp, tempK, temp
	add	pA, pA, temp
	mov	temp, #32
	mul	temp, tempK, temp
	add	pB, pB, temp
	#endif
	#if defined(LEFT)
	add	tempOffset, tempOffset, #1
	#endif

.Ldtrmm6x8_L4_END:
	add	origPB, origPB, origK, lsl #5
#if !defined(LEFT)
	add	tempOffset, tempOffset, #4
#endif

.Ldtrmm6x8_L2_BEGIN:
	tst	origN, #2
	beq	.Ldtrmm6x8_L1_BEGIN
	mov	pCRow0, pC
	add	pC, pC, LDC, lsl #1
#if defined(LEFT)
	mov	tempOffset, offset
#endif
	mov	pA, origPA
	mov	temp, #6
	udiv	counterI, origM, temp
	cbz	counterI, .Ldtrmm6x8_L2_Mrem
.Ldtrmm6x8_L2_M6:
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
	#else
	mov	pB, origPB
	mov	temp, #48
	mul	temp, tempOffset, temp
	add	pA, pA, temp
	mov	temp, #16
	mul	temp, tempOffset, temp
	add	pB, pB, temp
	#endif
	#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
	#elif defined(LEFT)
	add	tempK, tempOffset, #6
	#else
	add	tempK, tempOffset, #2
	#endif

	INIT6x2
	mov	counterL, tempK
	cbz	counterL, .Ldtrmm6x8_L2_M6s
.Ldtrmm6x8_L2_M6l:
	KERNEL6x2_SUB
	subs	counterL, counterL, #1
	bgt	.Ldtrmm6x8_L2_M6l
.Ldtrmm6x8_L2_M6s:
	SAVE6x2
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
	#if defined(LEFT)
	sub	tempK, tempK, #6
	#else
	sub	tempK, tempK, #2
	#endif
	mov	temp, #48
	mul	temp, tempK, temp
	add	pA, pA, temp
	mov	temp, #16
	mul	temp, tempK, temp
	add	pB, pB, temp
	#endif
	#if defined(LEFT)
	add	tempOffset, tempOffset, #6
	#endif

	subs	counterI, counterI, #1
	bgt	.Ldtrmm6x8_L2_M6
.Ldtrmm6x8_L2_Mrem:
	mov	temp, #6
	udiv	counterI, origM, temp
	msub	counterI, counterI, temp, origM
	cmp	counterI, #4
	blt	.Ldtrmm6x8_L2_M2c
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
	#else
	mov	pB, origPB
	mov	temp, #32
	mul	temp, tempOffset, temp
	add	pA, pA, temp
	mov	temp, #16
	mul	temp, tempOffset, temp
	add	pB, pB, temp
	#endif
	#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
	#elif defined(LEFT)
	add	tempK, tempOffset, #4
	#else
	add	tempK, tempOffset, #2
	#endif

	INIT4x2
	mov	counterL, tempK
	cbz	counterL, .Ldtrmm6x8_L2_M4s
.Ldtrmm6x8_L2_M4l:
	KERNEL4x2_SUB
	subs	counterL, counterL, #1
	bgt	.Ldtrmm6x8_L2_M4l
.Ldtrmm6x8_L2_M4s:
	SAVE4x2
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
	#if defined(LEFT)
	sub	tempK, tempK, #4
	#else
	sub	tempK, tempK, #2
	#endif
	mov	temp, #32
	mul	temp, tempK, temp
	add	pA, pA, temp
	mov	temp, #16
	mul	temp, tempK, temp
	add	pB, pB, temp
	#endif
	#if defined(LEFT)
	add	tempOffset, tempOffset, #4
	#endif

	sub	counterI, counterI, #4
.Ldtrmm6x8_L2_M2c:
	cmp	counterI, #2
	blt	.Ldtrmm6x8_L2_M1c
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
	#else
	mov	pB, origPB
	mov	temp, #16
	mul	temp, tempOffset, temp
	add	pA, pA, temp
	mov	temp, #16
	mul	temp, tempOffset, temp
	add	pB, pB, temp
	#endif
	#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
	#elif defined(LEFT)
	add	tempK, tempOffset, #2
	#else
	add	tempK, tempOffset, #2
	#endif

	INIT2x2
	mov	counterL, tempK
	cbz	counterL, .Ldtrmm6x8_L2_M2s
.Ldtrmm6x8_L2_M2l:
	KERNEL2x2_SUB
	subs	counterL, counterL, #1
	bgt	.Ldtrmm6x8_L2_M2l
.Ldtrmm6x8_L2_M2s:
	SAVE2x2
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
	#if defined(LEFT)
	sub	tempK, tempK, #2
	#else
	sub	tempK, tempK, #2
	#endif
	mov	temp, #16
	mul	temp, tempK, temp
	add	pA, pA, temp
	mov	temp, #16
	mul	temp, tempK, temp
	add	pB, pB, temp
	#endif
	#if defined(LEFT)
	add	tempOffset, tempOffset, #2
	#endif

	sub	counterI, counterI, #2
.Ldtrmm6x8_L2_M1c:
	cbz	counterI, .Ldtrmm6x8_L2_END
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
	#else
	mov	pB, origPB
	mov	temp, #8
	mul	temp, tempOffset, temp
	add	pA, pA, temp
	mov	temp, #16
	mul	temp, tempOffset, temp
	add	pB, pB, temp
	#endif
	#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
	#elif defined(LEFT)
	add	tempK, tempOffset, #1
	#else
	add	tempK, tempOffset, #2
	#endif

	INIT1x2
	mov	counterL, tempK
	cbz	counterL, .Ldtrmm6x8_L2_M1s
.Ldtrmm6x8_L2_M1l:
	KERNEL1x2_SUB
	subs	counterL, counterL, #1
	bgt	.Ldtrmm6x8_L2_M1l
.Ldtrmm6x8_L2_M1s:
	SAVE1x2
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
	#if defined(LEFT)
	sub	tempK, tempK, #1
	#else
	sub	tempK, tempK, #2
	#endif
	mov	temp, #8
	mul	temp, tempK, temp
	add	pA, pA, temp
	mov	temp, #16
	mul	temp, tempK, temp
	add	pB, pB, temp
	#endif
	#if defined(LEFT)
	add	tempOffset, tempOffset, #1
	#endif

.Ldtrmm6x8_L2_END:
	add	origPB, origPB, origK, lsl #4
#if !defined(LEFT)
	add	tempOffset, tempOffset, #2
#endif

.Ldtrmm6x8_L1_BEGIN:
	tst	origN, #1
	beq	.Ldtrmm6x8_L999
	mov	pCRow0, pC
#if defined(LEFT)
	mov	tempOffset, offset
#endif
	mov	pA, origPA
	mov	temp, #6
	udiv	counterI, origM, temp
	cbz	counterI, .Ldtrmm6x8_L1_Mrem
.Ldtrmm6x8_L1_M6:
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
	#else
	mov	pB, origPB
	mov	temp, #48
	mul	temp, tempOffset, temp
	add	pA, pA, temp
	mov	temp, #8
	mul	temp, tempOffset, temp
	add	pB, pB, temp
	#endif
	#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
	#elif defined(LEFT)
	add	tempK, tempOffset, #6
	#else
	add	tempK, tempOffset, #1
	#endif

	INIT6x1
	mov	counterL, tempK
	cbz	counterL, .Ldtrmm6x8_L1_M6s
.Ldtrmm6x8_L1_M6l:
	KERNEL6x1_SUB
	subs	counterL, counterL, #1
	bgt	.Ldtrmm6x8_L1_M6l
.Ldtrmm6x8_L1_M6s:
	SAVE6x1
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
	#if defined(LEFT)
	sub	tempK, tempK, #6
	#else
	sub	tempK, tempK, #1
	#endif
	mov	temp, #48
	mul	temp, tempK, temp
	add	pA, pA, temp
	mov	temp, #8
	mul	temp, tempK, temp
	add	pB, pB, temp
	#endif
	#if defined(LEFT)
	add	tempOffset, tempOffset, #6
	#endif

	subs	counterI, counterI, #1
	bgt	.Ldtrmm6x8_L1_M6
.Ldtrmm6x8_L1_Mrem:
	mov	temp, #6
	udiv	counterI, origM, temp
	msub	counterI, counterI, temp, origM
	cmp	counterI, #4
	blt	.Ldtrmm6x8_L1_M2c
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
	#else
	mov	pB, origPB
	mov	temp, #32
	mul	temp, tempOffset, temp
	add	pA, pA, temp
	mov	temp, #8
	mul	temp, tempOffset, temp
	add	pB, pB, temp
	#endif
	#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
	#elif defined(LEFT)
	add	tempK, tempOffset, #4
	#else
	add	tempK, tempOffset, #1
	#endif

	INIT4x1
	mov	counterL, tempK
	cbz	counterL, .Ldtrmm6x8_L1_M4s
.Ldtrmm6x8_L1_M4l:
	KERNEL4x1_SUB
	subs	counterL, counterL, #1
	bgt	.Ldtrmm6x8_L1_M4l
.Ldtrmm6x8_L1_M4s:
	SAVE4x1
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
	#if defined(LEFT)
	sub	tempK, tempK, #4
	#else
	sub	tempK, tempK, #1
	#endif
	mov	temp, #32
	mul	temp, tempK, temp
	add	pA, pA, temp
	mov	temp, #8
	mul	temp, tempK, temp
	add	pB, pB, temp
	#endif
	#if defined(LEFT)
	add	tempOffset, tempOffset, #4
	#endif

	sub	counterI, counterI, #4
.Ldtrmm6x8_L1_M2c:
	cmp	counterI, #2
	blt	.Ldtrmm6x8_L1_M1c
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
	#else
	mov	pB, origPB
	mov	temp, #16
	mul	temp, tempOffset, temp
	add	pA, pA, temp
	mov	temp, #8
	mul	temp, tempOffset, temp
	add	pB, pB, temp
	#endif
	#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
	#elif defined(LEFT)
	add	tempK, tempOffset, #2
	#else
	add	tempK, tempOffset, #1
	#endif

	INIT2x1
	mov	counterL, tempK
	cbz	counterL, .Ldtrmm6x8_L1_M2s
.Ldtrmm6x8_L1_M2l:
	KERNEL2x1_SUB
	subs	counterL, counterL, #1
	bgt	.Ldtrmm6x8_L1_M2l
.Ldtrmm6x8_L1_M2s:
	SAVE2x1
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
	#if defined(LEFT)
	sub	tempK, tempK, #2
	#else
	sub	tempK, tempK, #1
	#endif
	mov	temp, #16
	mul	temp, tempK, temp
	add	pA, pA, temp
	mov	temp, #8
	mul	temp, tempK, temp
	add	pB, pB, temp
	#endif
	#if defined(LEFT)
	add	tempOffset, tempOffset, #2
	#endif

	sub	counterI, counterI, #2
.Ldtrmm6x8_L1_M1c:
	cbz	counterI, .Ldtrmm6x8_L999
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mov	pB, origPB
	#else
	mov	pB, origPB
	mov	temp, #8
	mul	temp, tempOffset, temp
	add	pA, pA, temp
	mov	temp, #8
	mul	temp, tempOffset, temp
	add	pB, pB, temp
	#endif
	#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub	tempK, origK, tempOffset
	#elif defined(LEFT)
	add	tempK, tempOffset, #1
	#else
	add	tempK, tempOffset, #1
	#endif

	INIT1x1
	mov	counterL, tempK
	cbz	counterL, .Ldtrmm6x8_L1_M1s
.Ldtrmm6x8_L1_M1l:
	KERNEL1x1_SUB
	subs	counterL, counterL, #1
	bgt	.Ldtrmm6x8_L1_M1l
.Ldtrmm6x8_L1_M1s:
	SAVE1x1
	#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub	tempK, origK, tempOffset
	#if defined(LEFT)
	sub	tempK, tempK, #1
	#else
	sub	tempK, tempK, #1
	#endif
	mov	temp, #8
	mul	temp, tempK, temp
	add	pA, pA, temp
	mov	temp, #8
	mul	temp, tempK, temp
	add	pB, pB, temp
	#endif
	#if defined(LEFT)
	add	tempOffset, tempOffset, #1
	#endif

.Ldtrmm6x8_L999:
	mov	x0, #0
	ldp	d8, d9, [sp, #(0 * 16)]
	ldp	d10, d11, [sp, #(1 * 16)]
	ldp	d12, d13, [sp, #(2 * 16)]
	ldp	d14, d15, [sp, #(3 * 16)]
	ldp	d16, d17, [sp, #(4 * 16)]
	ldp	x18, x19, [sp, #(5 * 16)]
	ldp	x20, x21, [sp, #(6 * 16)]
	ldp	x22, x23, [sp, #(7 * 16)]
	ldp	x24, x25, [sp, #(8 * 16)]
	ldp	x26, x27, [sp, #(9 * 16)]
	ldr	x28, [sp, #(10 * 16)]
	add	sp, sp, #(11 * 16)
	ret
	EPILOGUE
