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

/*                   X0          X1          X2          s0         X3        x4       x5           x6 */
/*int CNAME(BLASLONG bm,BLASLONG bn,BLASLONG bk,FLOAT alpha0,FLOAT* ba,FLOAT* bb,FLOAT* C,BLASLONG ldc )*/

#define origM		x0
#define origN		x1
#define origK		x2
#define origPA		x3
#define origPB		x4
#define pC		x5
#define LDC		x6
#define temp		x7
#define counterL	x8
#define counterI	x9
#define counterJ	x10
#define pB		x11
#define pCRow0		x12
#define pCRow1		x13
#define pCRow2		x14
#define pCRow3		x15
#define pA		x16
#define alpha		x17

#define alpha0		d10
#define alphaZ		z10.d
#define alphaV0		v10.d[0]

#define A_PRE_SIZE	2560
#define B_PRE_SIZE	448
#define C_PRE_SIZE	128

// 00 origM
// 01 origN
// 02 origK
// 03 origPA
// 04 origPB
// 05 pC
// 06 origLDC -> LDC
// 07 temp
// 08 counterL
// 09 counterI
// 10 counterJ
// 11 pB
// 12 pCRow0
// 13 pCRow1
// 14 pCRow2
// 15 pCRow3
// 16 pA
// 17
// 18 must save
// 19 must save
// 20 must save
// 21 must save
// 22 must save
// 23 must save
// 24 must save
// 25 must save
// 26 must save
// 27 must save
// 28 must save
// 29 frame
// 30 link
// 31 sp

//v00 ALPHA -> pA0_0
//v01 pA0_1
//v02 pA0_2
//v03 pA0_3
//v04 pA0_4
//v05 pA0_5
//v06 pA0_6
//v07 pA0_7
//v08 must save pB0_0
//v09 must save pB0_1
//v10 must save pB0_2 --> ALPHA0
//v11 must save pB0_3
//v12 must save pB1_0
//v13 must save pB1_1
//v14 must save pB1_2
//v15 must save pB1_3
//v16 must save C0
//v17 must save C1
//v18 must save C2
//v19 must save C3
//v20 must save C4
//v21 must save C5
//v22 must save C6
//v23 must save C7

/*******************************************************************************
* Macro definitions
*******************************************************************************/

.macro INITv1x8
    dup         z16.d, #0
    dup         z17.d, #0
    dup         z18.d, #0
    dup         z19.d, #0
    dup         z20.d, #0
    dup         z21.d, #0
    dup         z22.d, #0
    dup         z23.d, #0
.endm

.macro KERNELv1x8_I
    ld1d  z0.d, p1/z, [pA] 
    ld1d  z1.d, p1/z, [pA, x18, lsl #3]   // next one
    //incb  pA, all, mul #2
	add	pA, pA, x18, lsl #4	// pA = pA + cnt_active * 2 * 8

    ld1rd  z8.d, p0/z,  [pB]
    ld1rd  z9.d, p0/z,  [pB, 8]
    ld1rd  z10.d, p0/z, [pB, 16]
    ld1rd  z11.d, p0/z, [pB, 24]
    ld1rd  z12.d, p0/z, [pB, 32]
    ld1rd  z13.d, p0/z, [pB, 40]
    ld1rd  z14.d, p0/z, [pB, 48]
    ld1rd  z15.d, p0/z, [pB, 56]

    add pB, pB, 64

    fmla z16.d, p1/m, z0.d, z8.d
    ld1rd  z8.d, p0/z,  [pB]
    fmla z17.d, p1/m, z0.d, z9.d
    ld1rd  z9.d, p0/z,  [pB, 8]
    fmla z18.d, p1/m, z0.d, z10.d
    ld1rd  z10.d, p0/z, [pB, 16]
    fmla z19.d, p1/m, z0.d, z11.d
    ld1rd  z11.d, p0/z, [pB, 24]
    fmla z20.d, p1/m, z0.d, z12.d
    ld1rd  z12.d, p0/z, [pB, 32]
    fmla z21.d, p1/m, z0.d, z13.d
    ld1rd  z13.d, p0/z, [pB, 40]
    fmla z22.d, p1/m, z0.d, z14.d
    ld1rd  z14.d, p0/z, [pB, 48]
    fmla z23.d, p1/m, z0.d, z15.d
    ld1rd  z15.d, p0/z, [pB, 56]

    add pB, pB, 64
.endm

.macro KERNELv1x8_M1
    ld1d  z1.d, p1/z, [pA] 
	add	pA, pA, x18, lsl #3	// pA = pA + cnt_active  * 8

    fmla z16.d, p1/m, z0.d, z8.d
    ld1rd  z8.d, p0/z,  [pB]
    fmla z17.d, p1/m, z0.d, z9.d
    ld1rd  z9.d, p0/z,  [pB, 8]
    fmla z18.d, p1/m, z0.d, z10.d
    ld1rd  z10.d, p0/z, [pB, 16]
    fmla z19.d, p1/m, z0.d, z11.d
    ld1rd  z11.d, p0/z, [pB, 24]
    fmla z20.d, p1/m, z0.d, z12.d
    ld1rd  z12.d, p0/z, [pB, 32]
    fmla z21.d, p1/m, z0.d, z13.d
    ld1rd  z13.d, p0/z, [pB, 40]
    fmla z22.d, p1/m, z0.d, z14.d
    ld1rd  z14.d, p0/z, [pB, 48]
    fmla z23.d, p1/m, z0.d, z15.d
    ld1rd  z15.d, p0/z, [pB, 56]

    add pB, pB, 64
.endm

.macro KERNELv1x8_M2
    ld1d  z0.d, p1/z, [pA] 
	add	pA, pA, x18, lsl #3	// pA = pA + cnt_active  * 8

    fmla z16.d, p1/m, z1.d, z8.d
    ld1rd  z8.d, p0/z,  [pB]
    fmla z17.d, p1/m, z1.d, z9.d
    ld1rd  z9.d, p0/z,  [pB, 8]
    fmla z18.d, p1/m, z1.d, z10.d
    ld1rd  z10.d, p0/z, [pB, 16]
    fmla z19.d, p1/m, z1.d, z11.d
    ld1rd  z11.d, p0/z, [pB, 24]
    fmla z20.d, p1/m, z1.d, z12.d
    ld1rd  z12.d, p0/z, [pB, 32]
    fmla z21.d, p1/m, z1.d, z13.d
    ld1rd  z13.d, p0/z, [pB, 40]
    fmla z22.d, p1/m, z1.d, z14.d
    ld1rd  z14.d, p0/z, [pB, 48]
    fmla z23.d, p1/m, z1.d, z15.d
    ld1rd  z15.d, p0/z, [pB, 56]

    add pB, pB, 64
.endm

.macro KERNELv1x8_E
    fmla z16.d, p1/m, z1.d, z8.d
    fmla z17.d, p1/m, z1.d, z9.d
    fmla z18.d, p1/m, z1.d, z10.d
    fmla z19.d, p1/m, z1.d, z11.d
    fmla z20.d, p1/m, z1.d, z12.d
    fmla z21.d, p1/m, z1.d, z13.d
    fmla z22.d, p1/m, z1.d, z14.d
    fmla z23.d, p1/m, z1.d, z15.d
.endm

.macro KERNELv1x8_SUB
    ld1d  z0.d, p1/z, [pA] 
	add	pA, pA, x18, lsl #3	// pA = pA + cnt_active  * 8

    ld1rd  z8.d, p0/z,  [pB]
    ld1rd  z9.d, p0/z,  [pB, 8]
    ld1rd  z10.d, p0/z, [pB, 16]
    ld1rd  z11.d, p0/z, [pB, 24]
    ld1rd  z12.d, p0/z, [pB, 32]
    ld1rd  z13.d, p0/z, [pB, 40]
    ld1rd  z14.d, p0/z, [pB, 48]
    ld1rd  z15.d, p0/z, [pB, 56]

    add pB, pB, 64

    fmla z16.d, p1/m, z0.d, z8.d
    fmla z17.d, p1/m, z0.d, z9.d
    fmla z18.d, p1/m, z0.d, z10.d
    fmla z19.d, p1/m, z0.d, z11.d
    fmla z20.d, p1/m, z0.d, z12.d
    fmla z21.d, p1/m, z0.d, z13.d
    fmla z22.d, p1/m, z0.d, z14.d
    fmla z23.d, p1/m, z0.d, z15.d

.endm

.macro SAVEv1x8
	dup	alphaZ, alpha

	prfm	PLDL2KEEP, [pCRow0, #C_PRE_SIZE]

	add	pCRow1, pCRow0, LDC
    ld1d  z24.d, p1/z, [pCRow0] 
    fmla z24.d, p1/m, z16.d, alphaZ
    st1d  z24.d, p1, [pCRow0]

	add	pCRow2, pCRow1, LDC
    ld1d  z25.d, p1/z, [pCRow1] 
    fmla z25.d, p1/m, z17.d, alphaZ
    st1d  z25.d, p1, [pCRow1]

	add	pCRow1, pCRow2, LDC
    ld1d  z26.d, p1/z, [pCRow2] 
    fmla z26.d, p1/m, z18.d, alphaZ
    st1d z26.d, p1, [pCRow2]

	add	pCRow2, pCRow1, LDC
    ld1d  z27.d, p1/z, [pCRow1] 
    fmla z27.d, p1/m, z19.d, alphaZ
    st1d  z27.d, p1, [pCRow1]

	add	pCRow1, pCRow2, LDC
    ld1d  z28.d, p1/z, [pCRow2] 
    fmla z28.d, p1/m, z20.d, alphaZ
    st1d  z28.d, p1, [pCRow2]

	add	pCRow2, pCRow1, LDC
    ld1d  z29.d, p1/z, [pCRow1] 
    fmla z29.d, p1/m, z21.d, alphaZ
    st1d  z29.d, p1, [pCRow1]

	add	pCRow1, pCRow2, LDC
    ld1d  z30.d, p1/z, [pCRow2] 
    fmla z30.d, p1/m, z22.d, alphaZ
    st1d  z30.d, p1, [pCRow2]

	add	pCRow2, pCRow1, LDC
    ld1d  z31.d, p1/z, [pCRow1] 
    fmla z31.d, p1/m, z23.d, alphaZ
    st1d  z31.d, p1, [pCRow1]

	add	pCRow0, pCRow0, x18, lsl #3	// pC = pC + cnt_active  * 8

.endm

/******************************************************************************/

.macro INITv1x4
    dup         z16.d, #0
    dup         z17.d, #0
    dup         z18.d, #0
    dup         z19.d, #0
.endm

.macro KERNELv1x4_SUB
    ld1d  z0.d, p1/z, [pA] 
	add	pA, pA, x18, lsl #3	// pA = pA + cnt_active  * 8

    ld1rd  z8.d, p0/z,  [pB]
    ld1rd  z9.d, p0/z,  [pB, 8]
    ld1rd  z10.d, p0/z, [pB, 16]
    ld1rd  z11.d, p0/z, [pB, 24]

    add pB, pB, 32

    fmla z16.d, p1/m, z0.d, z8.d
    fmla z17.d, p1/m, z0.d, z9.d
    fmla z18.d, p1/m, z0.d, z10.d
    fmla z19.d, p1/m, z0.d, z11.d

.endm

.macro SAVEv1x4
	dup	alphaZ, alpha

	prfm	PLDL2KEEP, [pCRow0, #C_PRE_SIZE]

	add	pCRow1, pCRow0, LDC
    ld1d  z24.d, p1/z, [pCRow0] 
    fmla z24.d, p1/m, z16.d, alphaZ
    st1d  z24.d, p1, [pCRow0]

	add	pCRow2, pCRow1, LDC
    ld1d  z25.d, p1/z, [pCRow1] 
    fmla z25.d, p1/m, z17.d, alphaZ
    st1d  z25.d, p1, [pCRow1]

	add	pCRow1, pCRow2, LDC
    ld1d  z26.d, p1/z, [pCRow2] 
    fmla z26.d, p1/m, z18.d, alphaZ
    st1d z26.d, p1, [pCRow2]

	add	pCRow2, pCRow1, LDC
    ld1d  z27.d, p1/z, [pCRow1] 
    fmla z27.d, p1/m, z19.d, alphaZ
    st1d  z27.d, p1, [pCRow1]

	add	pCRow0, pCRow0, x18, lsl #3	// pC = pC + cnt_active  * 8

.endm

/******************************************************************************/

.macro INITv1x2
    dup         z16.d, #0
    dup         z17.d, #0
.endm

.macro KERNELv1x2_SUB
    ld1d  z0.d, p1/z, [pA] 
	add	pA, pA, x18, lsl #3	// pA = pA + cnt_active  * 8

    ld1rd  z8.d, p0/z,  [pB]
    ld1rd  z9.d, p0/z,  [pB, 8]

    add pB, pB, 16

    fmla z16.d, p1/m, z0.d, z8.d
    fmla z17.d, p1/m, z0.d, z9.d

.endm

.macro SAVEv1x2
	dup	alphaZ, alpha

	prfm	PLDL2KEEP, [pCRow0, #C_PRE_SIZE]

	add	pCRow1, pCRow0, LDC
    ld1d  z24.d, p1/z, [pCRow0] 
    fmla z24.d, p1/m, z16.d, alphaZ
    st1d  z24.d, p1, [pCRow0]

	add	pCRow2, pCRow1, LDC
    ld1d  z25.d, p1/z, [pCRow1] 
    fmla z25.d, p1/m, z17.d, alphaZ
    st1d  z25.d, p1, [pCRow1]

	add	pCRow0, pCRow0, x18, lsl #3	// pC = pC + cnt_active  * 8

.endm

/******************************************************************************/

.macro INITv1x1
    dup         z16.d, #0
.endm

.macro KERNELv1x1_SUB
    ld1d  z0.d, p1/z, [pA] 
	add	pA, pA, x18, lsl #3	// pA = pA + cnt_active  * 8

    ld1rd  z8.d, p0/z,  [pB]

    add pB, pB, 8

    fmla z16.d, p1/m, z0.d, z8.d

.endm

.macro SAVEv1x1
	dup	alphaZ, alpha

	prfm	PLDL2KEEP, [pCRow0, #C_PRE_SIZE]

	add	pCRow1, pCRow0, LDC
    ld1d  z24.d, p1/z, [pCRow0] 
    fmla z24.d, p1/m, z16.d, alphaZ
    st1d  z24.d, p1, [pCRow0]


	add	pCRow0, pCRow0, x18, lsl #3	// pC = pC + cnt_active  * 8

.endm


/*******************************************************************************
* End of macro definitions
*******************************************************************************/

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

	prfm	PLDL1KEEP, [origPB]
	prfm	PLDL1KEEP, [origPA]

	fmov	alpha, d0

	lsl	LDC, LDC, #3			// ldc = ldc * 8
    ptrue p0.d                  // create true predicate 

	mov	pB, origPB

	mov	counterJ, origN
	asr 	counterJ, counterJ, #3		// J = J / 8
	cmp 	counterJ, #0
	ble	.Ldgemm_kernel_L4_BEGIN

/******************************************************************************/

	.align 5
.Ldgemm_kernel_L8_BEGIN:
	mov	pCRow0, pC

    add pC, pCRow0, LDC, lsl #3 // add 8 x LDC

	mov	pA, origPA			// pA = start of A array

.Ldgemm_kernel_L8_Mv1_BEGIN:

    mov counterI, #0
    whilelt p1.d, counterI, origM               //SVE instruction
    cntp x18, p0, p1.d
	/* mov	counterI, origM */
	/* asr 	counterI, counterI, #3		// counterI = counterI / 8 */
	/* cmp 	counterI, #0 */
	/* ble	.Ldgemm_kernel_L4_M4_BEGIN */

	.align 5
.Ldgemm_kernel_L8_Mv1_20:

	mov	pB, origPB
    INITv1x8                     // fill with zeros

	asr 	counterL , origK, #3		// L = K / 8
	cmp	counterL , #2			// is there at least 4 to do?
	blt	.Ldgemm_kernel_L8_Mv1_32

	KERNELv1x8_I
	KERNELv1x8_M2
	KERNELv1x8_M1
	KERNELv1x8_M2
	KERNELv1x8_M1
	KERNELv1x8_M2
	KERNELv1x8_M1
	KERNELv1x8_M2

	subs	counterL, counterL, #2		// subtract 2
	ble	.Ldgemm_kernel_L8_Mv1_22a

	.align 5
.Ldgemm_kernel_L8_Mv1_22:

	KERNELv1x8_M1
	KERNELv1x8_M2
	KERNELv1x8_M1
	KERNELv1x8_M2
	KERNELv1x8_M1
	KERNELv1x8_M2
	KERNELv1x8_M1
	KERNELv1x8_M2

	subs	counterL, counterL, #1
	bgt	.Ldgemm_kernel_L8_Mv1_22

	.align 5
.Ldgemm_kernel_L8_Mv1_22a:

	KERNELv1x8_M1
	KERNELv1x8_M2
	KERNELv1x8_M1
	KERNELv1x8_M2
	KERNELv1x8_M1
	KERNELv1x8_M2
	KERNELv1x8_M1
	KERNELv1x8_E

	b	 .Ldgemm_kernel_L8_Mv1_44

	.align 5
.Ldgemm_kernel_L8_Mv1_32:

	tst	counterL, #1
	ble	.Ldgemm_kernel_L8_Mv1_40

	KERNELv1x8_I
	KERNELv1x8_M2
	KERNELv1x8_M1
	KERNELv1x8_M2
	KERNELv1x8_M1
	KERNELv1x8_M2
	KERNELv1x8_M1
	KERNELv1x8_E


	b	.Ldgemm_kernel_L8_Mv1_44

.Ldgemm_kernel_L8_Mv1_40:

	INITv1x8

.Ldgemm_kernel_L8_Mv1_44:

	ands	counterL , origK, #7
	ble	.Ldgemm_kernel_L8_Mv1_100

	.align 5
.Ldgemm_kernel_L8_Mv1_46:

	KERNELv1x8_SUB

	subs	counterL, counterL, #1
	bne	.Ldgemm_kernel_L8_Mv1_46

.Ldgemm_kernel_L8_Mv1_100:
	prfm	PLDL1KEEP, [pA]
	prfm	PLDL1KEEP, [pA, #64]
	prfm	PLDL1KEEP, [origPB]

	SAVEv1x8

.Ldgemm_kernel_L8_Mv1_END:

    incd    counterI
    whilelt p1.d, counterI, origM             //SVE instruction
    cntp x18, p0, p1.d
    b.any   .Ldgemm_kernel_L8_Mv1_20   

.Ldgemm_kernel_L8_END:

	lsl	temp, origK, #6 
	add	origPB, origPB, temp		// B = B + K * 8 * 8

	subs	counterJ, counterJ , #1		// j--
	bgt	.Ldgemm_kernel_L8_BEGIN

/******************************************************************************/
/******************************************************************************/

	.align 5
.Ldgemm_kernel_L4_BEGIN:

	mov	counterJ , origN
	tst	counterJ , #4
	ble	.Ldgemm_kernel_L2_BEGIN


	mov	pCRow0, pC

    add pC, pCRow0, LDC, lsl #2 // add 4 x LDC

	mov	pA, origPA			// pA = start of A array

.Ldgemm_kernel_L4_Mv1_BEGIN:

    mov counterI, #0
    whilelt p1.d, counterI, origM               //SVE instruction
    cntp x18, p0, p1.d

	.align 5
.Ldgemm_kernel_L4_Mv1_20:

	mov	pB, origPB
    INITv1x4                     // fill with zeros

	asr 	counterL , origK, #3		// L = K / 8
	cmp	counterL , #0			// is there at least 4 to do?
	blt	.Ldgemm_kernel_L4_Mv1_44

	.align 5
.Ldgemm_kernel_L4_Mv1_22:

	KERNELv1x4_SUB
	KERNELv1x4_SUB
	KERNELv1x4_SUB
	KERNELv1x4_SUB
	KERNELv1x4_SUB
	KERNELv1x4_SUB
	KERNELv1x4_SUB
	KERNELv1x4_SUB

	subs	counterL, counterL, #1
	bgt	.Ldgemm_kernel_L4_Mv1_22

.Ldgemm_kernel_L4_Mv1_44:

	ands	counterL , origK, #7
	ble	.Ldgemm_kernel_L4_Mv1_100

	.align 5
.Ldgemm_kernel_L4_Mv1_46:

	KERNELv1x4_SUB

	subs	counterL, counterL, #1
	bne	.Ldgemm_kernel_L4_Mv1_46

.Ldgemm_kernel_L4_Mv1_100:
	prfm	PLDL1KEEP, [pA]
	prfm	PLDL1KEEP, [pA, #64]
	prfm	PLDL1KEEP, [origPB]

	SAVEv1x4

.Ldgemm_kernel_L4_Mv1_END:

    incd    counterI
    whilelt p1.d, counterI, origM             //SVE instruction
    cntp x18, p0, p1.d
    b.any   .Ldgemm_kernel_L4_Mv1_20   


.Ldgemm_kernel_L4_END:
	add	origPB, origPB, origK, lsl #5	// B = B + K * 4 * 8

/******************************************************************************/
/******************************************************************************/

	.align 5
.Ldgemm_kernel_L2_BEGIN:

	mov	counterJ , origN
	tst	counterJ , #2
	ble	.Ldgemm_kernel_L1_BEGIN

	mov	pCRow0, pC

    add pC, pCRow0, LDC, lsl #1 // add 2 x LDC

	mov	pA, origPA			// pA = start of A array

.Ldgemm_kernel_L2_Mv1_BEGIN:

    mov counterI, #0
    whilelt p1.d, counterI, origM               //SVE instruction
    cntp x18, p0, p1.d

	.align 5
.Ldgemm_kernel_L2_Mv1_20:

	mov	pB, origPB
    INITv1x2                     // fill with zeros

	asr 	counterL , origK, #3		// L = K / 8
	cmp	counterL , #0			// is there at least 4 to do?
	blt	.Ldgemm_kernel_L2_Mv1_44

	.align 5
.Ldgemm_kernel_L2_Mv1_22:

	KERNELv1x2_SUB
	KERNELv1x2_SUB
	KERNELv1x2_SUB
	KERNELv1x2_SUB
	KERNELv1x2_SUB
	KERNELv1x2_SUB
	KERNELv1x2_SUB
	KERNELv1x2_SUB

	subs	counterL, counterL, #1
	bgt	.Ldgemm_kernel_L2_Mv1_22

.Ldgemm_kernel_L2_Mv1_44:

	ands	counterL , origK, #7
	ble	.Ldgemm_kernel_L2_Mv1_100

	.align 5
.Ldgemm_kernel_L2_Mv1_46:

	KERNELv1x2_SUB

	subs	counterL, counterL, #1
	bne	.Ldgemm_kernel_L2_Mv1_46

.Ldgemm_kernel_L2_Mv1_100:
	prfm	PLDL1KEEP, [pA]
	prfm	PLDL1KEEP, [pA, #64]
	prfm	PLDL1KEEP, [origPB]

	SAVEv1x2

.Ldgemm_kernel_L2_Mv1_END:

    incd    counterI
    whilelt p1.d, counterI, origM             //SVE instruction
    cntp x18, p0, p1.d
    b.any   .Ldgemm_kernel_L2_Mv1_20   


.Ldgemm_kernel_L2_END:
	add	origPB, origPB, origK, lsl #4	// B = B + K * 2 * 8

/******************************************************************************/
/******************************************************************************/

	.align 5
.Ldgemm_kernel_L1_BEGIN:

	mov	counterJ , origN
	tst	counterJ , #1
	ble	.Ldgemm_kernel_L999 // done

	mov	pCRow0, pC

    add pC, pCRow0, LDC, lsl #1 // add 2 x LDC

	mov	pA, origPA			// pA = start of A array

.Ldgemm_kernel_L1_Mv1_BEGIN:

    mov counterI, #0
    whilelt p1.d, counterI, origM               //SVE instruction
    cntp x18, p0, p1.d

	.align 5
.Ldgemm_kernel_L1_Mv1_20:

	mov	pB, origPB
    INITv1x1                     // fill with zeros

	asr 	counterL , origK, #3		// L = K / 8
	cmp	counterL , #0			// is there at least 4 to do?
	blt	.Ldgemm_kernel_L1_Mv1_44

	.align 5
.Ldgemm_kernel_L1_Mv1_22:

	KERNELv1x1_SUB
	KERNELv1x1_SUB
	KERNELv1x1_SUB
	KERNELv1x1_SUB
	KERNELv1x1_SUB
	KERNELv1x1_SUB
	KERNELv1x1_SUB
	KERNELv1x1_SUB

	subs	counterL, counterL, #1
	bgt	.Ldgemm_kernel_L1_Mv1_22

.Ldgemm_kernel_L1_Mv1_44:

	ands	counterL , origK, #7
	ble	.Ldgemm_kernel_L1_Mv1_100

	.align 5
.Ldgemm_kernel_L1_Mv1_46:

	KERNELv1x1_SUB

	subs	counterL, counterL, #1
	bne	.Ldgemm_kernel_L1_Mv1_46

.Ldgemm_kernel_L1_Mv1_100:
	prfm	PLDL1KEEP, [pA]
	prfm	PLDL1KEEP, [pA, #64]
	prfm	PLDL1KEEP, [origPB]

	SAVEv1x1

.Ldgemm_kernel_L1_Mv1_END:

    incd    counterI
    whilelt p1.d, counterI, origM             //SVE instruction
    cntp x18, p0, p1.d
    b.any   .Ldgemm_kernel_L1_Mv1_20   


.Ldgemm_kernel_L1_END:

/******************************************************************************/

.Ldgemm_kernel_L999:
	mov	x0, #0				// set return value
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
	add	sp, sp, #(11*16)
	ret

	EPILOGUE

