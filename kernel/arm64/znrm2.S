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
#define I	x5	/* loop variable */

/*******************************************************************************
* Macro definitions
*******************************************************************************/

#if !defined(DOUBLE)
#define TMPF	s6
#define SSQ	s0
#define TMPVF	{v6.s}[0]
#define SZ	4
#else
#define TMPF	d6
#define SSQ	d0
#define TMPVF	{v6.d}[0]
#define SZ	8
#endif

/******************************************************************************/

.macro KERNEL_F1
#if !defined(DOUBLE)
	ld1	{v1.2s}, [X], #8
	fmul	v1.2s, v1.2s, v1.2s
	faddp	TMPF, v1.2s
	fadd	SSQ, SSQ, TMPF
#else
	ld1	{v1.2d}, [X], #16
	fmul	v1.2d, v1.2d, v1.2d
	faddp	TMPF, v1.2d
	fadd	SSQ, SSQ, TMPF
#endif
.endm

.macro KERNEL_F8
#if !defined(DOUBLE)
	ld1	{v1.4s, v2.4s}, [X], #32
	fmla	v0.4s, v1.4s, v1.4s
	fmla	v5.4s, v2.4s, v2.4s
	ld1	{v3.4s,v4.4s}, [X], #32
	fmla	v0.4s, v3.4s, v3.4s
	fmla	v5.4s, v4.4s, v4.4s
	PRFM	PLDL1KEEP, [X, #1024]
#else // DOUBLE
	ld1	{v1.2d, v2.2d}, [X], #32
	fmla	v0.2d, v1.2d, v1.2d
	fmla	v5.2d, v2.2d, v2.2d
	ld1	{v3.2d, v4.2d}, [X], #32
	fmla	v0.2d, v3.2d, v3.2d
	fmla	v5.2d, v4.2d, v4.2d

	ld1	{v16.2d, v17.2d}, [X], #32
	fmla	v0.2d, v16.2d, v16.2d
	fmla	v5.2d, v17.2d, v17.2d
	ld1	{v18.2d, v19.2d}, [X], #32
	fmla	v0.2d, v18.2d, v18.2d
	fmla	v5.2d, v19.2d, v19.2d
#endif
.endm

.macro nrm2_kernel_F8_FINALIZE
#if !defined(DOUBLE)
	fadd	v0.4s, v0.4s, v5.4s
	ext	v1.16b, v0.16b, v0.16b, #8
	fadd	v0.2s, v0.2s,  v1.2s
	faddp	SSQ, v0.2s
#else
	fadd	v0.2d, v0.2d, v5.2d
	faddp	SSQ, v0.2d
#endif
.endm

.macro INIT_S
#if !defined(DOUBLE)
	lsl	INC_X, INC_X, #3
	ld1	{v1.2s}, [X], INC_X
	fmul	v1.2s, v1.2s, v1.2s
	faddp	SSQ, v1.2s
#else
	lsl	INC_X, INC_X, #4
	ld1	{v1.2d}, [X], INC_X
	fmul	v1.2d, v1.2d, v1.2d
	faddp	SSQ, v1.2d
#endif
.endm

.macro KERNEL_S1
#if !defined(DOUBLE)
	ld1	{v1.2s}, [X], INC_X
	fmul	v1.2s, v1.2s, v1.2s
	faddp	TMPF, v1.2s
	fadd	SSQ, SSQ, TMPF
#else
	ld1	{v1.2d}, [X], INC_X
	fmul	v1.2d, v1.2d, v1.2d
	faddp	TMPF, v1.2d
	fadd	SSQ, SSQ, TMPF
#endif
.endm

/*******************************************************************************
* End of macro definitions
*******************************************************************************/

	PROLOGUE

#if !defined(DOUBLE)
	fmov	SSQ, wzr
	fmov	s5, SSQ
#else
	fmov	SSQ, xzr
	fmov	d5, SSQ
#endif

	cmp	N, xzr
	ble	nrm2_kernel_zero
	cmp	INC_X, xzr
	ble	nrm2_kernel_zero
	cmp	INC_X, #1
	bne	nrm2_kernel_S_BEGIN

nrm2_kernel_F_BEGIN:

	asr	I, N, #3
	cmp	I, xzr
	beq	nrm2_kernel_F1_INIT

nrm2_kernel_F8:

	KERNEL_F8

	subs	I, I, #1
	bne	nrm2_kernel_F8

	nrm2_kernel_F8_FINALIZE

nrm2_kernel_F1:

	ands	I, N, #7
	ble	nrm2_kernel_L999

nrm2_kernel_F10:

	KERNEL_F1

	subs    I, I, #1
        bne     nrm2_kernel_F10

	b	nrm2_kernel_L999

nrm2_kernel_F1_INIT:

	b	nrm2_kernel_F1

nrm2_kernel_S_BEGIN:

	INIT_S

	subs	N, N, #1
	ble	nrm2_kernel_L999

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
	fsqrt	SSQ, SSQ
	ret

nrm2_kernel_zero:
	ret

	EPILOGUE
