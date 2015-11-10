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
#define	X	x3	/* X vector address */
#define	INC_X	x4	/* X stride */
#define I	x5	/* loop variable */

/*******************************************************************************
* Macro definitions
*******************************************************************************/

#if !defined(DOUBLE)
#define	DA_R	s0	/* real scale input value */
#define	DA_I	s1	/* imaginary scale input value */
#else
#define	DA_R	d0	/* real scale input value */
#define	DA_I	d1	/* imaginary scale input value */
#endif

/******************************************************************************/

.macro INIT

#if !defined(DOUBLE)
	ins	v0.s[1], v0.s[0]		// v0 = DA_R, DA_R 
	fneg	s2, DA_I
	ins	v1.s[1], v2.s[0]		// v1 = -DA_I, DA_I 
	ext	v1.8b, v1.8b, v1.8b, #4		// v1 = DA_I, -DA_I
#else
	ins	v0.d[1], v0.d[0]		// v0 = DA_R, DA_R 
	fneg	d2, DA_I
	ins	v1.d[1], v2.d[0]		// v1 = DA_I, DA_I 
	ext	v1.16b, v1.16b, v1.16b, #8	// v1 = DA_I, -DA_I
#endif

.endm

.macro KERNEL_F1

#if !defined(DOUBLE)
	ld1	{v2.2s}, [X]			// X1, X0
	ext	v3.8b, v2.8b, v2.8b, #4		// X0, X1
	fmul	v2.2s, v2.2s, v0.2s		// DA_R*X1, DA_R*X0
	fmla	v2.2s, v3.2s, v1.2s		// DA_R*X1+DA_I*X0, DA_R*X0-DA_I*X1
	st1	{v2.2s}, [X], #8
#else
	ld1	{v2.2d}, [X]			// X1, X0
	ext	v3.16b, v2.16b, v2.16b, #8	// X0, X1
	fmul	v2.2d, v2.2d, v0.2d		// DA_R*X1, DA_R*X0
	fmla	v2.2d, v3.2d, v1.2d		// DA_R*X1+DA_I*X0, DA_R*X0-DA_I*X1
	st1	{v2.2d}, [X], #16
#endif

.endm

.macro KERNEL_INIT_F4

#if !defined(DOUBLE)
	// Replicate the lower 2 floats into the upper 2 slots
	ins	v0.d[1], v0.d[0]		// v0 = DA_R, DA_R, DA_R, DA_R
	ins	v1.d[1], v1.d[0]		// v1 = DA_I, DA_I, DA_I, DA_I
#endif

.endm

.macro KERNEL_F4

#if !defined(DOUBLE)
	ld1	{v2.4s,v3.4s}, [X]		// V2 = X[3], X[2], X[1], X[0]
						// V3 = X[7], X[6], X[5], X[4]

	ext	v6.8b, v2.8b, v2.8b, #4		// V6 =  -  ,  -  , X[0], X[1]
	ins	v6.s[2], v2.s[3]		// V6 =  -  , X[3], X[0], X[1]
	ins	v6.s[3], v2.s[2]		// V6 = X[2], X[3], X[0], X[1]
	fmul	v2.4s, v0.4s, v2.4s		// X'[ix]   += DA_R * X[ix]
						// X'[ix+1] += DA_R * X[ix+1]
	fmla	v2.4s, v1.4s, v6.4s		// X'[ix]   += -DA_I * X[ix+1]
						// X'[ix+1] += DA_I * X[ix]

	ext	v7.8b, v3.8b, v3.8b, #4		// V7 =  -  ,  -  , X[4], X[5]
	ins	v7.s[2], v3.s[3]		// V7 =  -  , X[7], X[4], X[5]
	ins	v7.s[3], v3.s[2]		// V7 = X[6], X[7], X[4], X[5]
	fmul	v3.4s, v0.4s, v3.4s		// X'[ix]   += DA_R * X[ix]
						// X'[ix+1] += DA_R * X[ix+1]
	fmla	v3.4s, v1.4s, v7.4s		// X'[ix]   += -DA_I * X[ix+1]
						// X'[ix+1] += DA_I * X[ix]

	st1	{v2.4s,v3.4s}, [X], #32
#else // DOUBLE
	ld1	{v2.2d,v3.2d,v4.2d,v5.2d}, [X]	// CX0, CX1, CX2, CX3
	ext	v20.16b, v2.16b, v2.16b, #8	// X[ix], X[ix+1]
	ext	v21.16b, v3.16b, v3.16b, #8	// X[ix], X[ix+1]
	ext	v22.16b, v4.16b, v4.16b, #8	// X[ix], X[ix+1]
	ext	v23.16b, v5.16b, v5.16b, #8	// X[ix], X[ix+1]

	fmul	v2.2d, v0.2d, v2.2d
	fmla	v2.2d, v1.2d, v20.2d

	fmul	v3.2d, v0.2d, v3.2d
	fmla	v3.2d, v1.2d, v21.2d
	st1	{v2.2d,v3.2d}, [X], #32

	fmul	v4.2d, v0.2d, v4.2d
	fmla	v4.2d, v1.2d, v22.2d

	fmul	v5.2d, v0.2d, v5.2d
	fmla	v5.2d, v1.2d, v23.2d
	st1	{v4.2d,v5.2d}, [X], #32
#endif
	PRFM	PLDL1KEEP, [X, #1024]
.endm

.macro INIT_S

#if !defined(DOUBLE)
	lsl	INC_X, INC_X, #3
#else
	lsl	INC_X, INC_X, #4
#endif

.endm

.macro KERNEL_S1

#if !defined(DOUBLE)
	ld1	{v2.2s}, [X]			// X1, X0
	ext	v3.8b, v2.8b, v2.8b, #4		// X0, X1
	fmul	v2.2s, v2.2s, v0.2s		// DA_R*X1, DA_R*X0
	fmla	v2.2s, v3.2s, v1.2s		// DA_R*X1+DA_I*X0, DA_R*X0-DA_I*X1
	st1	{v2.2s}, [X], INC_X
#else
	ld1	{v2.2d}, [X]			// X1, X0
	ext	v3.16b, v2.16b, v2.16b, #8	// X0, X1
	fmul	v2.2d, v2.2d, v0.2d		// DA_R*X1, DA_R*X0
	fmla	v2.2d, v3.2d, v1.2d		// DA_R*X1+DA_I*X0, DA_R*X0-DA_I*X1
	st1	{v2.2d}, [X], INC_X
#endif

.endm

/*******************************************************************************
* End of macro definitions
*******************************************************************************/

	PROLOGUE

	cmp	N, xzr
	ble	zscal_kernel_L999

	fcmp	DA_R, #0.0
	bne	zscal_kernel_1

	fcmp	DA_I, #0.0
	beq	zscal_kernel_zero

	// TODO: special case DA_R == 0 && DA_I != 0

zscal_kernel_1:

	// TODO: special case DA_R != 0 && DA_I == 0

	INIT

	cmp	INC_X, #1
	bne	zscal_kernel_S_BEGIN

zscal_kernel_F_BEGIN:

	asr	I, N, #2
	cmp	I, xzr
	beq	zscal_kernel_F1

	KERNEL_INIT_F4

zscal_kernel_F4:

	KERNEL_F4

	subs	I, I, #1
	bne	zscal_kernel_F4

zscal_kernel_F1:

	ands	I, N, #3
	ble	zscal_kernel_L999

zscal_kernel_F10:

	KERNEL_F1

	subs    I, I, #1
        bne     zscal_kernel_F10

	mov	w0, wzr
	ret

zscal_kernel_S_BEGIN:

	INIT_S

	asr	I, N, #2
	cmp	I, xzr
	ble	zscal_kernel_S1

zscal_kernel_S4:

	KERNEL_S1
	KERNEL_S1
	KERNEL_S1
	KERNEL_S1

	subs	I, I, #1
	bne	zscal_kernel_S4

zscal_kernel_S1:

	ands	I, N, #3
	ble	zscal_kernel_L999

zscal_kernel_S10:

	KERNEL_S1

	subs    I, I, #1
        bne     zscal_kernel_S10

zscal_kernel_L999:

	mov	w0, wzr
	ret

zscal_kernel_zero:

	INIT_S

zscal_kernel_Z1:

	stp	DA_R, DA_I, [X]
	add	X, X, INC_X
	subs    N, N, #1
        bne     zscal_kernel_Z1

	mov	w0, wzr
	ret

	EPILOGUE
