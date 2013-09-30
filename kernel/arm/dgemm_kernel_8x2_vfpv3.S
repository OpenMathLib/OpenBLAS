/***************************************************************************
Copyright (c) 2013, The OpenBLAS Project
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

/**************************************************************************************
* 2013/09/30 Saar
* 	 BLASTEST 		: OK
* 	 CTEST			: OK
* 	 TEST			: OK
*
*
* 2013/09/30 Saar
*	UNROLL_N		2
*	UNROLL_M		8
*	DGEMM_P			64
*	DGEMM_Q			64
*	DGEMM_R			512
*	A_PRE			192
*	B_PRE			32
*	C_PRE			64
*
* 	Performance on Odroid U2:
*
*		1 Core:		1.42 GFLOPS	ATLAS: 1.58	GFLOPS
*		2 Cores:	2.81 GFLOPS	ATLAS: -	GFLOPS
*		3 Cores:	4.05 GFLOPS	ATLAS: -	GFLOPS
*		4 Cores:	5.40 GFLOPS	ATLAS: 3.88	GFLOPS
**************************************************************************************/

#define ASSEMBLER
#include "common.h"

#define STACKSIZE 252

#define	OLD_M	r0
#define	OLD_N	r1
#define	OLD_K	r2
#define	OLD_A	r3
#define OLD_ALPHA d0

/******************************************************
* [fp, #-128] - [fp, #-64] is reserved
* for store and restore of floating point
* registers
*******************************************************/

#define C	[fp, #-248 ]
#define LDC	[fp, #-252 ]
#define M	[fp, #-256 ]
#define N	[fp, #-260 ]
#define K	[fp, #-264 ]
#define A	[fp, #-268 ]

#define ALPHA	[fp, #-276 ]

#define B	[fp, #4 ]
#define OLD_C	[fp, #8 ]
#define OLD_LDC	[fp, #12 ]

#define I	r0
#define J	r1
#define L	r2

#define	AO	r5
#define	BO	r6

#define	CO1	r8
#define	CO2	r9

#define K1	r7
#define BC	r12

#define A_PRE	192
#define A_PRE1	224
#define B_PRE	32
#define C_PRE	64

/**************************************************************************************
* Macro definitions
**************************************************************************************/

.macro INIT8x2

	vsub.f64	d8 , d8 ,  d8
	vmov.f64		d9 , d8
	vmov.f64		d10, d8
	vmov.f64		d11, d8
	vmov.f64		d12, d8
	vmov.f64		d13, d8
	vmov.f64		d14, d8
	vmov.f64		d15, d8
	vmov.f64		d16, d8
	vmov.f64		d17, d8
	vmov.f64		d18, d8
	vmov.f64		d19, d8
	vmov.f64		d20, d8
	vmov.f64		d21, d8
	vmov.f64		d22, d8
	vmov.f64		d23, d8

.endm



.macro KERNEL8x2

	fldmiad	BO!, { d24 , d25}
	fldd	d0, [ AO ]
	fmacd	d8  , d0,  d24
	fldd	d1, [ AO , #8  ]
	fmacd	d16 , d0,  d25
	fldd	d2, [ AO , #16  ]
	fmacd	d9  , d1,  d24
	fmacd	d17 , d1,  d25
	fldd	d3, [ AO , #24  ]
	fmacd	d10 , d2,  d24
	fmacd	d18 , d2,  d25
	fldd	d4, [ AO , #32  ]
	fmacd	d11 , d3,  d24
        pld     [AO , #A_PRE]
	fmacd	d19 , d3,  d25
	fldd	d5, [ AO , #40 ]

	fmacd	d12 , d4,  d24
	fmacd	d20 , d4,  d25
	fldd	d6, [ AO , #48 ]
	fmacd	d13 , d5,  d24
	fmacd	d21 , d5,  d25

	fldd	d7, [ AO , #56  ]
	fmacd	d14 , d6,  d24
	fmacd	d22 , d6,  d25
        pld     [AO , #A_PRE+32]
	fmacd	d15 , d7,  d24
	add	AO, AO, #64
	fmacd	d23 , d7,  d25

.endm

.macro SAVE8x2

	vldr		d0, ALPHA

	fldd	d24, [CO1]
	fldd	d25, [CO1, #8 ]
		
	fmacd	d24, d0 , d8
	fldd	d8  , [CO2]	
	fldd	d26, [CO1, #16]
	fmacd	d25, d0 , d9
	fldd	d9  , [CO2, #8 ]	
	fldd	d27, [CO1, #24]
	fmacd	d26, d0 , d10
	fldd	d10 , [CO2, #16 ]	
	fldd	d28, [CO1, #32]
	fmacd	d27, d0 , d11
	fldd	d11 , [CO2, #24 ]	
	fldd	d29, [CO1, #40]
	fmacd	d28, d0 , d12
	fldd	d12 , [CO2, #32 ]	
	fldd	d30, [CO1, #48]
	fmacd	d29, d0 , d13
	fldd	d13 , [CO2, #40 ]	
	fldd	d31, [CO1, #56]
	fmacd	d30, d0 , d14
	fldd	d14 , [CO2, #48 ]	
	fmacd	d31, d0 , d15
	fldd	d15 , [CO2, #56 ]	

	
	fmacd	d8 , d0 , d16
	fstd	d24, [CO1]
	fmacd	d9 , d0 , d17
	fstd	d25, [CO1, #8 ]
	fstd	d8 , [CO2]
	fmacd	d10, d0 , d18
	fstd	d26, [CO1, #16 ]
	fstd	d9 , [CO2, #8  ]
	fmacd	d11, d0 , d19
	fstd	d27, [CO1, #24 ]
	fstd	d10, [CO2, #16 ]
	fmacd	d12, d0 , d20
	fstd	d28, [CO1, #32 ]
	fstd	d11, [CO2, #24 ]
	fmacd	d13, d0 , d21
	fstd	d29, [CO1, #40 ]
	fstd	d12, [CO2, #32 ]
	fmacd	d14, d0 , d22
	fstd	d30, [CO1, #48 ]
	fstd	d13, [CO2, #40 ]
	fmacd	d15, d0 , d23
	fstd	d31, [CO1, #56 ]
	fstd	d14, [CO2, #48 ]

	add	CO1, CO1, #64
	fstd	d15, [CO2, #56 ]
	add	CO2, CO2, #64


.endm

.macro SAVE8x2_BAD

        vldr            d0, ALPHA
        vldm            CO2, { d24, d25, d26 , d27 , d28 , d29 , d30 , d31 }

        vmul.f64        d8 , d0 , d8
        vmul.f64        d9 , d0 , d9
        vmul.f64        d10, d0 , d10
        vmul.f64        d11, d0 , d11
        vmul.f64        d12, d0 , d12
        vmul.f64        d13, d0 , d13
        vmul.f64        d14, d0 , d14
        vmul.f64        d15, d0 , d15

        vmul.f64        d16, d0 , d16
        vmul.f64        d17, d0 , d17
        vmul.f64        d18, d0 , d18
        vmul.f64        d19, d0 , d19
        vmul.f64        d20, d0 , d20
        vmul.f64        d21, d0 , d21
        vmul.f64        d22, d0 , d22
        vmul.f64        d23, d0 , d23

        vldm            CO1, { d0 , d1 , d2  , d3  , d4  , d5  , d6  , d7 }

        vadd.f64        d16, d16, d24
        vadd.f64        d17, d17, d25
        vadd.f64        d18, d18, d26
        vadd.f64        d19, d19, d27

        vadd.f64        d20, d20, d28
        vadd.f64        d21, d21, d29
        vadd.f64        d22, d22, d30
        vadd.f64        d23, d23, d31

        vadd.f64        d8 , d8 , d0
        vadd.f64        d9 , d9 , d1
        vadd.f64        d10, d10, d2
        vadd.f64        d11, d11, d3

        vadd.f64        d12, d12, d4
        vadd.f64        d13, d13, d5
        vadd.f64        d14, d14, d6
        vadd.f64        d15, d15, d7

        vstm            CO2!, { d16, d17, d18 , d19 , d20 , d21 , d22 , d23 }
        vstm            CO1!, { d8 , d9 , d10 , d11 , d12, d13  , d14 , d15 }

.endm



/*************************************************************************************/


.macro INIT4x2

	vsub.f64	d8 ,  d8 ,  d8
	vsub.f64	d9 ,  d9 ,  d9
	vsub.f64	d10,  d10,  d10
	vsub.f64	d11,  d11,  d11
	vsub.f64	d12,  d12,  d12
	vsub.f64	d13,  d13,  d13
	vsub.f64	d14,  d14,  d14
	vsub.f64	d15,  d15,  d15

.endm

.macro KERNEL4x2

	vldm		AO!, { d0, d1 , d2, d3 }
	vldm		BO!, { d4, d5 }

	vmul.f64	d6 , d0 , d4
	vmul.f64	d7 , d1 , d4
	vadd.f64	d8 , d8 , d6
	vadd.f64	d9 , d9 , d7

	vmul.f64	d6 , d2 , d4
	vmul.f64	d7 , d3 , d4
	vadd.f64	d10, d10, d6
	vadd.f64	d11, d11, d7

	vmul.f64	d6 , d0 , d5
	vmul.f64	d7 , d1 , d5
	vadd.f64	d12, d12, d6
	vadd.f64	d13, d13, d7

	vmul.f64	d6 , d2 , d5
	vmul.f64	d7 , d3 , d5
	vadd.f64	d14, d14, d6
	vadd.f64	d15, d15, d7

.endm

.macro SAVE4x2

	vldr		d0, ALPHA

	vmul.f64	d8 , d0 , d8
	vmul.f64	d9 , d0 , d9
	vmul.f64	d10, d0 , d10
	vmul.f64	d11, d0 , d11
	vmul.f64	d12, d0 , d12
	vmul.f64	d13, d0 , d13
	vmul.f64	d14, d0 , d14
	vmul.f64	d15, d0 , d15
	
	vldm		CO1, { d0, d1 , d2 , d3 }
	vldm		CO2, { d4, d5 , d6 , d7 }

	vadd.f64	d8 , d8 , d0
	vadd.f64	d9 , d9 , d1
	vadd.f64	d10, d10, d2
	vadd.f64	d11, d11, d3

	vadd.f64	d12, d12, d4
	vadd.f64	d13, d13, d5
	vadd.f64	d14, d14, d6
	vadd.f64	d15, d15, d7

	vstm		CO1!, { d8 , d9 , d10 , d11 }
	vstm		CO2!, { d12, d13 ,d14 , d15 }

.endm



/*************************************************************************************/

.macro INIT2x2

	vsub.f64	d8 ,  d8 ,  d8
	vsub.f64	d9 ,  d9 ,  d9
	vsub.f64	d12,  d12,  d12
	vsub.f64	d13,  d13,  d13

.endm

.macro KERNEL2x2

	vldm		AO!, { d0, d1 }
	vldm		BO!, { d4, d5 }

	vmul.f64	d6 , d0 , d4
	vmul.f64	d7 , d1 , d4
	vadd.f64	d8 , d8 , d6
	vadd.f64	d9 , d9 , d7

	vmul.f64	d6 , d0 , d5
	vmul.f64	d7 , d1 , d5
	vadd.f64	d12, d12, d6
	vadd.f64	d13, d13, d7

.endm

.macro SAVE2x2

	vldr		d0, ALPHA

	vmul.f64	d8 , d0 , d8
	vmul.f64	d9 , d0 , d9
	vmul.f64	d12, d0 , d12
	vmul.f64	d13, d0 , d13
	
	vldm		CO1, { d0, d1 }
	vldm		CO2, { d4, d5 }

	vadd.f64	d8 , d8 , d0
	vadd.f64	d9 , d9 , d1

	vadd.f64	d12, d12, d4
	vadd.f64	d13, d13, d5

	vstm		CO1!, { d8 , d9  }
	vstm		CO2!, { d12, d13 }

.endm

/*************************************************************************************/

.macro INIT1x2

	vsub.f64	d8 ,  d8 ,  d8
	vsub.f64	d12,  d12,  d12

.endm

.macro KERNEL1x2

	vldm		AO!, { d0 }
	vldm		BO!, { d4, d5 }

	vmul.f64	d6 , d0 , d4
	vadd.f64	d8 , d8 , d6

	vmul.f64	d6 , d0 , d5
	vadd.f64	d12, d12, d6

.endm

.macro SAVE1x2

	vldr		d0, ALPHA

	vmul.f64	d8 , d0 , d8
	vmul.f64	d12, d0 , d12
	
	vldm		CO1, { d0 }
	vldm		CO2, { d4 }

	vadd.f64	d8 , d8 , d0
	vadd.f64	d12, d12, d4

	vstm		CO1!, { d8 }
	vstm		CO2!, { d12}

.endm

/*************************************************************************************/

.macro INIT8x1

	vsub.f64	d8 ,  d8 ,  d8
	vsub.f64	d9 ,  d9 ,  d9
	vsub.f64	d10,  d10,  d10
	vsub.f64	d11,  d11,  d11
	vsub.f64	d12,  d12,  d12
	vsub.f64	d13,  d13,  d13
	vsub.f64	d14,  d14,  d14
	vsub.f64	d15,  d15,  d15

.endm

.macro KERNEL8x1

	vldm		AO!, { d0, d1 , d2, d3, d4 , d5 , d6 , d7 }
	vldm		BO!, { d24 }

	vmul.f64	d26 , d0 , d24
	vmul.f64	d27 , d1 , d24
	vadd.f64	d8  , d8 , d26
	vadd.f64	d9  , d9 , d27

	vmul.f64	d28 , d2 , d24
	vmul.f64	d29 , d3 , d24
	vadd.f64	d10 , d10, d28
	vadd.f64	d11 , d11, d29

	vmul.f64	d26 , d4 , d24
	vmul.f64	d27 , d5 , d24
	vadd.f64	d12 , d12, d26
	vadd.f64	d13 , d13, d27

	vmul.f64	d28 , d6 , d24
	vmul.f64	d29 , d7 , d24
	vadd.f64	d14 , d14, d28
	vadd.f64	d15 , d15, d29


.endm

.macro SAVE8x1

	vldr		d0, ALPHA

	vmul.f64	d8 , d0 , d8
	vmul.f64	d9 , d0 , d9
	vmul.f64	d10, d0 , d10
	vmul.f64	d11, d0 , d11
	vmul.f64	d12, d0 , d12
	vmul.f64	d13, d0 , d13
	vmul.f64	d14, d0 , d14
	vmul.f64	d15, d0 , d15
	
	vldm		CO1, { d0, d1 , d2 , d3 , d4 , d5 , d6 , d7 }

	vadd.f64	d8 , d8 , d0
	vadd.f64	d9 , d9 , d1
	vadd.f64	d10, d10, d2
	vadd.f64	d11, d11, d3

	vadd.f64	d12, d12, d4
	vadd.f64	d13, d13, d5
	vadd.f64	d14, d14, d6
	vadd.f64	d15, d15, d7

	vstm		CO1!, { d8 , d9 , d10 , d11 , d12, d13 ,d14 , d15 }

.endm


/*************************************************************************************/

.macro INIT4x1

	vsub.f64	d8 ,  d8 ,  d8
	vsub.f64	d9 ,  d9 ,  d9
	vsub.f64	d10,  d10,  d10
	vsub.f64	d11,  d11,  d11

.endm

.macro KERNEL4x1

	vldm		AO!, { d0, d1 , d2, d3 }
	vldm		BO!, { d4 }

	vmul.f64	d6 , d0 , d4
	vmul.f64	d7 , d1 , d4
	vadd.f64	d8 , d8 , d6
	vadd.f64	d9 , d9 , d7

	vmul.f64	d6 , d2 , d4
	vmul.f64	d7 , d3 , d4
	vadd.f64	d10, d10, d6
	vadd.f64	d11, d11, d7

.endm

.macro SAVE4x1

	vldr		d0, ALPHA

	vmul.f64	d8 , d0 , d8
	vmul.f64	d9 , d0 , d9
	vmul.f64	d10, d0 , d10
	vmul.f64	d11, d0 , d11
	
	vldm		CO1, { d0, d1 , d2 , d3 }

	vadd.f64	d8 , d8 , d0
	vadd.f64	d9 , d9 , d1
	vadd.f64	d10, d10, d2
	vadd.f64	d11, d11, d3

	vstm		CO1!, { d8 , d9 , d10 , d11 }

.endm

/*************************************************************************************/

.macro INIT2x1

	vsub.f64	d8 ,  d8 ,  d8
	vsub.f64	d9 ,  d9 ,  d9

.endm

.macro KERNEL2x1

	vldm		AO!, { d0, d1 }
	vldm		BO!, { d4 }

	vmul.f64	d6 , d0 , d4
	vmul.f64	d7 , d1 , d4
	vadd.f64	d8 , d8 , d6
	vadd.f64	d9 , d9 , d7

.endm

.macro SAVE2x1

	vldr		d0, ALPHA

	vmul.f64	d8 , d0 , d8
	vmul.f64	d9 , d0 , d9
	
	vldm		CO1, { d0, d1 }

	vadd.f64	d8 , d8 , d0
	vadd.f64	d9 , d9 , d1

	vstm		CO1!, { d8 , d9 }

.endm

/*************************************************************************************/

.macro INIT1x1

	vsub.f64	d8 ,  d8 ,  d8

.endm

.macro KERNEL1x1

	vldm		AO!, { d0 }
	vldm		BO!, { d4 }

	vmul.f64	d6 , d0 , d4
	vadd.f64	d8 , d8 , d6

.endm

.macro SAVE1x1

	vldr		d0, ALPHA

	vmul.f64	d8 , d0 , d8
	
	vldm		CO1, { d0 }

	vadd.f64	d8 , d8 , d0

	vstm		CO1!, { d8 }

.endm



/**************************************************************************************
* End of macro definitions
**************************************************************************************/

	PROLOGUE

	.align 5

	push	{r4 - r9, fp}
	add	fp, sp, #24
	sub	sp, sp, #STACKSIZE				// reserve stack

	str	OLD_M, M
	str	OLD_N, N
	str	OLD_K, K
	str	OLD_A, A
	vstr	OLD_ALPHA, ALPHA

	sub	r3, fp, #128
	vstm	r3, { d8 - d15} 				// store floating point registers

	ldr	r3, OLD_LDC
	lsl	r3, r3, #3					// ldc = ldc * 8
	str	r3, LDC

	ldr	r3, OLD_C
	str	r3, C

	ldr	K1, K
	ldr	BC, B

	ldr	J, N
	asrs	J, J, #1					// J = J / 2
	ble	_L1_BEGIN

_L2_BEGIN:
	
	ldr	CO1, C						// CO1 = C
	ldr	r4 , LDC
	add	CO2, CO1, r4					// CO2 = C + LDC
	add	r3 , CO2, r4					// C = CO2 + LDC
	str	r3 , C						// store C

	ldr	AO, A						// AO = A
        pld     [AO , #A_PRE-96]
        pld     [AO , #A_PRE-64]
        pld     [AO , #A_PRE-32]

_L2_M8_BEGIN:

	ldr	I, M
	asrs	I, I, #3					// I = I / 8
	ble	_L2_M4_BEGIN

_L2_M8_20:

	pld	[CO1, #C_PRE]
	pld	[CO1, #C_PRE+32]
	pld	[CO2, #C_PRE]
	pld	[CO2, #C_PRE+32]
	INIT8x2

	mov	BO, BC
	asrs	L , K1, #3					// L = L / 8
	ble	_L2_M8_40
	.align 5

_L2_M8_22:

        pld     [BO , #B_PRE]
	KERNEL8x2
	KERNEL8x2
        pld     [BO , #B_PRE]
	KERNEL8x2
	KERNEL8x2

        pld     [BO , #B_PRE]
	KERNEL8x2
	KERNEL8x2
        pld     [BO , #B_PRE]
	KERNEL8x2
	KERNEL8x2

	subs	L, L, #1
	bgt	_L2_M8_22
	

_L2_M8_40:
	
	ands	L , K1, #7					// L = L % 8
	ble	_L2_M8_100

_L2_M8_42:

	KERNEL8x2

	subs	L, L, #1
	bgt	_L2_M8_42
	
_L2_M8_100:

	SAVE8x2

_L2_M8_END:

	subs	I, I, #1
	bgt	_L2_M8_20


_L2_M4_BEGIN:

	ldr	I, M
	tst	I , #7
	ble	_L2_END

	tst	I , #4
	ble	_L2_M2_BEGIN

_L2_M4_20:

	INIT4x2

	mov	BO, BC
	asrs	L , K1, #3					// L = L / 8
	ble	_L2_M4_40
	.align 5

_L2_M4_22:

	KERNEL4x2
	KERNEL4x2
	KERNEL4x2
	KERNEL4x2

	KERNEL4x2
	KERNEL4x2
	KERNEL4x2
	KERNEL4x2

	subs	L, L, #1
	bgt	_L2_M4_22
	

_L2_M4_40:
	
	ands	L , K1, #7					// L = L % 8
	ble	_L2_M4_100

_L2_M4_42:

	KERNEL4x2

	subs	L, L, #1
	bgt	_L2_M4_42
	
_L2_M4_100:

	SAVE4x2

_L2_M4_END:



_L2_M2_BEGIN:

	tst	I, #2					// I = I / 2
	ble	_L2_M1_BEGIN

_L2_M2_20:

	INIT2x2

	mov	BO, BC
	asrs	L , K1, #2					// L = L / 4
	ble	_L2_M2_40

_L2_M2_22:

	KERNEL2x2
	KERNEL2x2
	KERNEL2x2
	KERNEL2x2

	subs	L, L, #1
	bgt	_L2_M2_22
	

_L2_M2_40:
	
	ands	L , K1, #3					// L = L % 4
	ble	_L2_M2_100

_L2_M2_42:

	KERNEL2x2

	subs	L, L, #1
	bgt	_L2_M2_42
	
_L2_M2_100:

	SAVE2x2

_L2_M2_END:


_L2_M1_BEGIN:

	tst	I, #1					// I = I % 2
	ble	_L2_END

_L2_M1_20:

	INIT1x2

	mov	BO, BC
	asrs	L , K1, #2					// L = L / 4
	ble	_L2_M1_40

_L2_M1_22:

	KERNEL1x2
	KERNEL1x2
	KERNEL1x2
	KERNEL1x2

	subs	L, L, #1
	bgt	_L2_M1_22
	

_L2_M1_40:
	
	ands	L , K1, #3					// L = L % 4
	ble	_L2_M1_100

_L2_M1_42:

	KERNEL1x2

	subs	L, L, #1
	bgt	_L2_M1_42
	
_L2_M1_100:

	SAVE1x2


_L2_END:

	mov	r3, BC
	mov	r4, K1
	lsl	r4, r4, #4					// k * 2 * 8
	add	r3, r3, r4					// B = B + K * 2 * 8
	mov	BC, r3
	
	subs	J , #1						// j--
	bgt	_L2_BEGIN


_L1_BEGIN:	

	ldr	J, N
	tst	J , #1					// J = J % 2
	ble	_L999

	ldr	CO1, C						// CO1 = C
	ldr	r4 , LDC
	add	r3 , CO1, r4					// C = CO1 + LDC
	str	r3 , C						// store C

	ldr	AO, A						// AO = A



_L1_M8_BEGIN:

	ldr	I, M
	asrs	I, I, #3					// I = I / 8
	ble	_L1_M4_BEGIN

_L1_M8_20:

	INIT8x1

	mov	BO, BC
	asrs	L , K1, #3					// L = L / 8
	ble	_L1_M8_40

_L1_M8_22:

	KERNEL8x1
	KERNEL8x1
	KERNEL8x1
	KERNEL8x1
	KERNEL8x1
	KERNEL8x1
	KERNEL8x1
	KERNEL8x1

	subs	L, L, #1
	bgt	_L1_M8_22
	

_L1_M8_40:
	
	ands	L , K1, #7					// L = L % 8
	ble	_L1_M8_100

_L1_M8_42:

	KERNEL8x1

	subs	L, L, #1
	bgt	_L1_M8_42
	
_L1_M8_100:

	SAVE8x1

_L1_M8_END:

	subs	I, I, #1
	bgt	_L1_M8_20




_L1_M4_BEGIN:

	ldr	I, M
	tst	I, #7					// I = I % 8
	ble	_L1_END	

	tst	I, #4					// I = I % 8
	ble	_L1_M2_BEGIN

_L1_M4_20:

	INIT4x1

	mov	BO, BC
	asrs	L , K1, #2					// L = L / 4
	ble	_L1_M4_40

_L1_M4_22:

	KERNEL4x1
	KERNEL4x1
	KERNEL4x1
	KERNEL4x1

	subs	L, L, #1
	bgt	_L1_M4_22
	

_L1_M4_40:
	
	ands	L , K1, #3					// L = L % 4
	ble	_L1_M4_100

_L1_M4_42:

	KERNEL4x1

	subs	L, L, #1
	bgt	_L1_M4_42
	
_L1_M4_100:

	SAVE4x1

_L1_M4_END:




_L1_M2_BEGIN:

	tst	I, #2					// I = I % 4
	ble	_L1_M1_BEGIN

_L1_M2_20:

	INIT2x1

	mov	BO, BC
	asrs	L , K1, #2					// L = L / 4
	ble	_L1_M2_40

_L1_M2_22:

	KERNEL2x1
	KERNEL2x1
	KERNEL2x1
	KERNEL2x1

	subs	L, L, #1
	bgt	_L1_M2_22
	

_L1_M2_40:
	
	ands	L , K1, #3					// L = L % 4
	ble	_L1_M2_100

_L1_M2_42:

	KERNEL2x1

	subs	L, L, #1
	bgt	_L1_M2_42
	
_L1_M2_100:

	SAVE2x1

_L1_M2_END:



_L1_M1_BEGIN:

	tst	I, #1					// I = I % 4
	ble	_L1_END

_L1_M1_20:

	INIT1x1

	mov	BO, BC
	asrs	L , K1, #2					// L = L / 4
	ble	_L1_M1_40

_L1_M1_22:

	KERNEL1x1
	KERNEL1x1
	KERNEL1x1
	KERNEL1x1

	subs	L, L, #1
	bgt	_L1_M1_22
	

_L1_M1_40:
	
	ands	L , K1, #3					// L = L % 4
	ble	_L1_M1_100

_L1_M1_42:

	KERNEL1x1

	subs	L, L, #1
	bgt	_L1_M1_42
	
_L1_M1_100:

	SAVE1x1


_L1_END:



_L999:

	sub	r3, fp, #128
	vldm	r3, { d8 - d15}					// restore floating point registers

	movs	r0, #0						// set return value
	sub	sp, fp, #24
	pop	{r4 - r9, fp}
	bx	lr

	EPILOGUE

