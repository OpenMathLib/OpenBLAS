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
* 2013/09/28 Saar
* 	 BLASTEST 		: OK
* 	 CTEST			: OK
* 	 TEST			: OK
*
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


#define KK      [fp, #-240 ]
#define KKK     [fp, #-244]
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
#define OFFSET  [fp, #16 ]

#define I	r0
#define J	r1
#define L	r2

#define	AO	r5
#define	BO	r6

#define	CO1	r8
#define	CO2	r9

#define K1	r7
#define BC	r12

#define A_PRE	128
#define A_PRE1	160
#define B_PRE	128
#define C_PRE	32

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

.macro KERNEL8x2_START

	vldm		BO!, { d24 , d25}
	vldm		AO!, { d0, d1 , d2, d3, d4 , d5 , d6 , d7 }

	vmul.f64	d26 , d0 , d24
	vmul.f64	d29 , d0 , d25
	vmul.f64	d27 , d1 , d24
	vmul.f64	d30 , d1 , d25
	vmul.f64	d28 , d2 , d24
	vmul.f64	d31 , d2 , d25

        pld     [AO , #A_PRE]

	vadd.f64	d8  , d8 , d26
	vadd.f64	d16 , d16, d29
	vadd.f64	d9  , d9 , d27
	vadd.f64	d17 , d17, d30
	vadd.f64	d10 , d10, d28
	vadd.f64	d18 , d18, d31

	vmul.f64	d26 , d3 , d24
	vmul.f64	d27 , d4 , d24
	vmul.f64	d28 , d5 , d24
	vmul.f64	d29 , d3 , d25
	vmul.f64	d30 , d4 , d25
	vmul.f64	d31 , d5 , d25

        pld     [AO , #A_PRE1]

	vadd.f64	d11 , d11, d26
	vadd.f64	d12 , d12, d27
	vadd.f64	d13 , d13, d28
	vadd.f64	d19 , d19, d29
	vadd.f64	d20 , d20, d30
	vadd.f64	d21 , d21, d31

	vmul.f64	d26 , d6 , d24
	vmul.f64	d27 , d7 , d24
	vmul.f64	d29 , d6 , d25
	vmul.f64	d30 , d7 , d25

	vldm		BO!, { d24 , d25}
	vldm		AO!, { d0, d1 , d2, d3, d4 , d5 , d6 , d7 }

	vadd.f64	d14 , d14, d26
	vadd.f64	d15 , d15, d27
	vadd.f64	d22 , d22, d29
	vadd.f64	d23 , d23, d30

.endm


.macro KERNEL8x2_M

	vmul.f64	d26 , d0 , d24
	vmul.f64	d29 , d0 , d25
	vmul.f64	d27 , d1 , d24
	vmul.f64	d30 , d1 , d25
	vmul.f64	d28 , d2 , d24
	vmul.f64	d31 , d2 , d25

        pld     [AO , #A_PRE]

	vadd.f64	d8  , d8 , d26
	vadd.f64	d16 , d16, d29
	vadd.f64	d9  , d9 , d27
	vadd.f64	d17 , d17, d30
	vadd.f64	d10 , d10, d28
	vadd.f64	d18 , d18, d31

	vmul.f64	d26 , d3 , d24
	vmul.f64	d27 , d4 , d24
	vmul.f64	d28 , d5 , d24
	vmul.f64	d29 , d3 , d25
	vmul.f64	d30 , d4 , d25
	vmul.f64	d31 , d5 , d25

        pld     [AO , #A_PRE1]

	vadd.f64	d11 , d11, d26
	vadd.f64	d12 , d12, d27
	vadd.f64	d13 , d13, d28
	vadd.f64	d19 , d19, d29
	vadd.f64	d20 , d20, d30
	vadd.f64	d21 , d21, d31

	vmul.f64	d26 , d6 , d24
	vmul.f64	d27 , d7 , d24
	vmul.f64	d29 , d6 , d25
	vmul.f64	d30 , d7 , d25

	vldm		BO!, { d24 , d25}
	vldm		AO!, { d0, d1 , d2, d3, d4 , d5 , d6 , d7 }

	vadd.f64	d14 , d14, d26
	vadd.f64	d15 , d15, d27
	vadd.f64	d22 , d22, d29
	vadd.f64	d23 , d23, d30

.endm


.macro KERNEL8x2_END

	vmul.f64	d26 , d0 , d24
	vmul.f64	d29 , d0 , d25
	vmul.f64	d27 , d1 , d24
	vmul.f64	d30 , d1 , d25
	vmul.f64	d28 , d2 , d24
	vmul.f64	d31 , d2 , d25

        pld     [AO , #A_PRE]

	vadd.f64	d8  , d8 , d26
	vadd.f64	d16 , d16, d29
	vadd.f64	d9  , d9 , d27
	vadd.f64	d17 , d17, d30
	vadd.f64	d10 , d10, d28
	vadd.f64	d18 , d18, d31

	vmul.f64	d26 , d3 , d24
	vmul.f64	d27 , d4 , d24
	vmul.f64	d28 , d5 , d24
	vmul.f64	d29 , d3 , d25
	vmul.f64	d30 , d4 , d25
	vmul.f64	d31 , d5 , d25

        pld     [AO , #A_PRE1]

	vadd.f64	d11 , d11, d26
	vadd.f64	d12 , d12, d27
	vadd.f64	d13 , d13, d28
	vadd.f64	d19 , d19, d29
	vadd.f64	d20 , d20, d30
	vadd.f64	d21 , d21, d31

	vmul.f64	d26 , d6 , d24
	vmul.f64	d27 , d7 , d24
	vmul.f64	d29 , d6 , d25
	vmul.f64	d30 , d7 , d25

	vadd.f64	d14 , d14, d26
	vadd.f64	d15 , d15, d27
	vadd.f64	d22 , d22, d29
	vadd.f64	d23 , d23, d30

.endm




.macro KERNEL8x2

	vldm		BO!, { d24 , d25}
	vldm		AO!, { d0, d1 , d2, d3, d4 , d5 , d6 , d7 }

	vmul.f64	d26 , d0 , d24
	vmul.f64	d29 , d0 , d25
	vmul.f64	d27 , d1 , d24
	vmul.f64	d30 , d1 , d25
	vmul.f64	d28 , d2 , d24
	vmul.f64	d31 , d2 , d25

        pld     [AO , #A_PRE]

	vadd.f64	d8  , d8 , d26
	vadd.f64	d16 , d16, d29
	vadd.f64	d9  , d9 , d27
	vadd.f64	d17 , d17, d30
	vadd.f64	d10 , d10, d28
	vadd.f64	d18 , d18, d31

	vmul.f64	d26 , d3 , d24
	vmul.f64	d27 , d4 , d24
	vmul.f64	d28 , d5 , d24
	vmul.f64	d29 , d3 , d25
	vmul.f64	d30 , d4 , d25
	vmul.f64	d31 , d5 , d25

        pld     [AO , #A_PRE1]

	vadd.f64	d11 , d11, d26
	vadd.f64	d12 , d12, d27
	vadd.f64	d13 , d13, d28
	vadd.f64	d19 , d19, d29
	vadd.f64	d20 , d20, d30
	vadd.f64	d21 , d21, d31

	vmul.f64	d26 , d6 , d24
	vmul.f64	d27 , d7 , d24
	vmul.f64	d29 , d6 , d25
	vmul.f64	d30 , d7 , d25

	vadd.f64	d14 , d14, d26
	vadd.f64	d15 , d15, d27
	vadd.f64	d22 , d22, d29
	vadd.f64	d23 , d23, d30

.endm

.macro SAVE8x2

	vldr		d0, ALPHA

	vmul.f64	d8 , d0 , d8
	vmul.f64	d9 , d0 , d9
	vmul.f64	d10, d0 , d10
	vmul.f64	d11, d0 , d11
	vmul.f64	d12, d0 , d12
	vmul.f64	d13, d0 , d13
	vmul.f64	d14, d0 , d14
	vmul.f64	d15, d0 , d15
	
	vmul.f64	d16, d0 , d16
	vmul.f64	d17, d0 , d17
	vmul.f64	d18, d0 , d18
	vmul.f64	d19, d0 , d19
	vmul.f64	d20, d0 , d20
	vmul.f64	d21, d0 , d21
	vmul.f64	d22, d0 , d22
	vmul.f64	d23, d0 , d23

	vstm		CO2!, { d16, d17, d18 , d19 , d20 , d21 , d22 , d23 }
	vstm		CO1!, { d8 , d9 , d10 , d11 , d12, d13  , d14 , d15 }

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

	ldr	BC, B

        ldr     r3, OFFSET
#ifndef LEFT
        neg     r3 , r3
#endif
        str     r3 , KK

	ldr	J, N
	asrs	J, J, #1					// J = J / 2
	ble	_L1_BEGIN

_L2_BEGIN:
	
	ldr	CO1, C						// CO1 = C
	ldr	r4 , LDC
	add	CO2, CO1, r4					// CO2 = C + LDC
	add	r3 , CO2, r4					// C = CO2 + LDC
	str	r3 , C						// store C

#if defined(LEFT)
        ldr     r3 , OFFSET
        str     r3 , KK
#endif

	ldr	AO, A						// AO = A

_L2_M8_BEGIN:

	ldr	I, M
	asrs	I, I, #3					// I = I / 8
	ble	_L2_M4_BEGIN

_L2_M8_20:

	pld	[CO1, #C_PRE]
	pld	[CO2, #C_PRE]

	INIT8x2

#if  (defined(LEFT) &&  defined(TRANSA)) || \
    (!defined(LEFT) && !defined(TRANSA))

	mov	BO, BC
#else
	mov	BO, BC
        ldr     r3 , KK
        lsls    r4 , r3 , #4                                    // 2 double values
        add     BO , BO , r4
        lsls    r4 , r3 , #6                                    // 8 double values
        add     AO , AO , r4

#endif

#ifndef TRMMKERNEL
        ldr     L , K
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
        ldr     L , K
        ldr     r3, KK
        sub     L , L, r3
        str     L , KKK
#else
        ldr     L , KK
#ifdef LEFT
        add     L , L , #8        // number of values in AO
#else
        add     L , L , #2        // number of values in BO
#endif
        str     L , KKK
#endif

	asrs	K1 , L, #3					// L = L / 8
	ble	_L2_M8_40
	.align 5

_L2_M8_22:

        pld     [BO , #B_PRE]
	KERNEL8x2_START
	KERNEL8x2_M
        pld     [BO , #B_PRE]
	KERNEL8x2_M
	KERNEL8x2_M

        pld     [BO , #B_PRE]
	KERNEL8x2_M
	KERNEL8x2_M
        pld     [BO , #B_PRE]
	KERNEL8x2_M
	KERNEL8x2_END

	subs	K1, K1, #1
	bgt	_L2_M8_22
	

_L2_M8_40:
	
	ands	K1 , L, #7					// L = L % 8
	ble	_L2_M8_100

_L2_M8_42:

	KERNEL8x2

	subs	K1, K1, #1
	bgt	_L2_M8_42
	
_L2_M8_100:

	SAVE8x2

#if  (defined(LEFT) &&  defined(TRANSA)) || \
    (!defined(LEFT) && !defined(TRANSA))
        ldr     r3 , K
        ldr     r4 , KKK
        sub     r3 , r3 , r4
        lsls    r4 , r3 , #4                    // 2 double values
        add     BO , BO , r4
        lsls    r4 , r3 , #6                    // 8 double values
        add     AO , AO , r4
#endif

#if defined(LEFT)
        ldr     r3 , KK
        add     r3 , r3 , #8                    // number of values in AO
        str     r3 , KK
#endif


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

#if  (defined(LEFT) &&  defined(TRANSA)) || \
    (!defined(LEFT) && !defined(TRANSA))

	mov	BO, BC
#else
	mov	BO, BC
        ldr     r3 , KK
        lsls    r4 , r3 , #4                                    // 2 double values
        add     BO , BO , r4
        lsls    r4 , r3 , #5                                    // 4 double values
        add     AO , AO , r4

#endif

#ifndef TRMMKERNEL
        ldr     L , K
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
        ldr     L , K
        ldr     r3, KK
        sub     L , L, r3
        str     L , KKK
#else
        ldr     L , KK
#ifdef LEFT
        add     L , L , #4        // number of values in AO
#else
        add     L , L , #2        // number of values in BO
#endif
        str     L , KKK
#endif

	asrs	K1, L, #3					// L = L / 8
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

	subs	K1, K1, #1
	bgt	_L2_M4_22
	

_L2_M4_40:
	
	ands	K1, L, #7					// L = L % 8
	ble	_L2_M4_100

_L2_M4_42:

	KERNEL4x2

	subs	K1, K1, #1
	bgt	_L2_M4_42
	
_L2_M4_100:

	SAVE4x2

#if  (defined(LEFT) &&  defined(TRANSA)) || \
    (!defined(LEFT) && !defined(TRANSA))
        ldr     r3 , K
        ldr     r4 , KKK
        sub     r3 , r3 , r4
        lsls    r4 , r3 , #4                    // 2 double values
        add     BO , BO , r4
        lsls    r4 , r3 , #5                    // 4 double values
        add     AO , AO , r4
#endif

#if defined(LEFT)
        ldr     r3 , KK
        add     r3 , r3 , #4                    // number of values in AO
        str     r3 , KK
#endif

_L2_M4_END:



_L2_M2_BEGIN:

	tst	I, #2					// I = I / 2
	ble	_L2_M1_BEGIN

_L2_M2_20:

	INIT2x2

#if  (defined(LEFT) &&  defined(TRANSA)) || \
    (!defined(LEFT) && !defined(TRANSA))

	mov	BO, BC
#else
	mov	BO, BC
        ldr     r3 , KK
        lsls    r4 , r3 , #4                                    // 2 double values
        add     BO , BO , r4
        lsls    r4 , r3 , #4                                    // 2 double values
        add     AO , AO , r4

#endif

#ifndef TRMMKERNEL
        ldr     L , K
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
        ldr     L , K
        ldr     r3, KK
        sub     L , L, r3
        str     L , KKK
#else
        ldr     L , KK
#ifdef LEFT
        add     L , L , #2        // number of values in AO
#else
        add     L , L , #2        // number of values in BO
#endif
        str     L , KKK
#endif

	asrs	K1, L, #2					// L = L / 4
	ble	_L2_M2_40

_L2_M2_22:

	KERNEL2x2
	KERNEL2x2
	KERNEL2x2
	KERNEL2x2

	subs	K1, K1, #1
	bgt	_L2_M2_22
	

_L2_M2_40:
	
	ands	K1, L, #3					// L = L % 4
	ble	_L2_M2_100

_L2_M2_42:

	KERNEL2x2

	subs	K1, K1, #1
	bgt	_L2_M2_42
	
_L2_M2_100:

	SAVE2x2

#if  (defined(LEFT) &&  defined(TRANSA)) || \
    (!defined(LEFT) && !defined(TRANSA))
        ldr     r3 , K
        ldr     r4 , KKK
        sub     r3 , r3 , r4
        lsls    r4 , r3 , #4                    // 2 double values
        add     BO , BO , r4
        lsls    r4 , r3 , #4                    // 2 double values
        add     AO , AO , r4
#endif

#if defined(LEFT)
        ldr     r3 , KK
        add     r3 , r3 , #2                    // number of values in AO
        str     r3 , KK
#endif


_L2_M2_END:


_L2_M1_BEGIN:

	tst	I, #1					// I = I % 2
	ble	_L2_END

_L2_M1_20:

	INIT1x2

#if  (defined(LEFT) &&  defined(TRANSA)) || \
    (!defined(LEFT) && !defined(TRANSA))

	mov	BO, BC
#else
	mov	BO, BC
        ldr     r3 , KK
        lsls    r4 , r3 , #4                                    // 2 double values
        add     BO , BO , r4
        lsls    r4 , r3 , #3                                    // 1 double value
        add     AO , AO , r4

#endif

#ifndef TRMMKERNEL
        ldr     L , K
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
        ldr     L , K
        ldr     r3, KK
        sub     L , L, r3
        str     L , KKK
#else
        ldr     L , KK
#ifdef LEFT
        add     L , L , #1        // number of values in AO
#else
        add     L , L , #2        // number of values in BO
#endif
        str     L , KKK
#endif

	asrs	K1, L, #2					// L = L / 4
	ble	_L2_M1_40

_L2_M1_22:

	KERNEL1x2
	KERNEL1x2
	KERNEL1x2
	KERNEL1x2

	subs	K1, K1, #1
	bgt	_L2_M1_22
	

_L2_M1_40:
	
	ands	K1, L, #3					// L = L % 4
	ble	_L2_M1_100

_L2_M1_42:

	KERNEL1x2

	subs	K1, K1, #1
	bgt	_L2_M1_42
	
_L2_M1_100:

	SAVE1x2

#if  (defined(LEFT) &&  defined(TRANSA)) || \
    (!defined(LEFT) && !defined(TRANSA))
        ldr     r3 , K
        ldr     r4 , KKK
        sub     r3 , r3 , r4
        lsls    r4 , r3 , #4                    // 2 double values
        add     BO , BO , r4
        lsls    r4 , r3 , #3                    // 1 double value
        add     AO , AO , r4
#endif

#if defined(LEFT)
        ldr     r3 , KK
        add     r3 , r3 , #1                    // number of values in AO
        str     r3 , KK
#endif



_L2_END:

	mov	r3, BC
	ldr	r4, K
	lsl	r4, r4, #4					// k * 2 * 8
	add	r3, r3, r4					// B = B + K * 2 * 8
	mov	BC, r3
	
#if !defined(LEFT)
        ldr     r3 , KK
        add     r3 , r3 , #2                                    // number of values in BO
        str     r3 , KK
#endif

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

#if defined(LEFT)
        ldr     r3 , OFFSET
        str     r3 , KK
#endif

	ldr	AO, A						// AO = A



_L1_M8_BEGIN:

	ldr	I, M
	asrs	I, I, #3					// I = I / 8
	ble	_L1_M4_BEGIN

_L1_M8_20:

	INIT8x1

#if  (defined(LEFT) &&  defined(TRANSA)) || \
    (!defined(LEFT) && !defined(TRANSA))

	mov	BO, BC
#else
	mov	BO, BC
        ldr     r3 , KK
        lsls    r4 , r3 , #3                                    // 1 double value
        add     BO , BO , r4
        lsls    r4 , r3 , #6                                    // 8 double values
        add     AO , AO , r4

#endif

#ifndef TRMMKERNEL
        ldr     L , K
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
        ldr     L , K
        ldr     r3, KK
        sub     L , L, r3
        str     L , KKK
#else
        ldr     L , KK
#ifdef LEFT
        add     L , L , #8        // number of values in AO
#else
        add     L , L , #1        // number of values in BO
#endif
        str     L , KKK
#endif

	asrs	K1, L, #3					// L = L / 8
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

	subs	K1, K1, #1
	bgt	_L1_M8_22
	

_L1_M8_40:
	
	ands	K1, L, #7					// L = L % 8
	ble	_L1_M8_100

_L1_M8_42:

	KERNEL8x1

	subs	K1, K1, #1
	bgt	_L1_M8_42
	
_L1_M8_100:

	SAVE8x1

#if  (defined(LEFT) &&  defined(TRANSA)) || \
    (!defined(LEFT) && !defined(TRANSA))
        ldr     r3 , K
        ldr     r4 , KKK
        sub     r3 , r3 , r4
        lsls    r4 , r3 , #3                    // 1 double value
        add     BO , BO , r4
        lsls    r4 , r3 , #6                    // 8 double values
        add     AO , AO , r4
#endif

#if defined(LEFT)
        ldr     r3 , KK
        add     r3 , r3 , #8                    // number of values in AO
        str     r3 , KK
#endif


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

#if  (defined(LEFT) &&  defined(TRANSA)) || \
    (!defined(LEFT) && !defined(TRANSA))

	mov	BO, BC
#else
	mov	BO, BC
        ldr     r3 , KK
        lsls    r4 , r3 , #3                                    // 1 double value
        add     BO , BO , r4
        lsls    r4 , r3 , #5                                    // 4 double values
        add     AO , AO , r4

#endif

#ifndef TRMMKERNEL
        ldr     L , K
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
        ldr     L , K
        ldr     r3, KK
        sub     L , L, r3
        str     L , KKK
#else
        ldr     L , KK
#ifdef LEFT
        add     L , L , #4        // number of values in AO
#else
        add     L , L , #1        // number of values in BO
#endif
        str     L , KKK
#endif

	asrs	K1, L, #2					// L = L / 4
	ble	_L1_M4_40

_L1_M4_22:

	KERNEL4x1
	KERNEL4x1
	KERNEL4x1
	KERNEL4x1

	subs	K1, K1, #1
	bgt	_L1_M4_22
	

_L1_M4_40:
	
	ands	K1, L, #3					// L = L % 4
	ble	_L1_M4_100

_L1_M4_42:

	KERNEL4x1

	subs	K1, K1, #1
	bgt	_L1_M4_42
	
_L1_M4_100:

	SAVE4x1

#if  (defined(LEFT) &&  defined(TRANSA)) || \
    (!defined(LEFT) && !defined(TRANSA))
        ldr     r3 , K
        ldr     r4 , KKK
        sub     r3 , r3 , r4
        lsls    r4 , r3 , #3                    // 1 double value
        add     BO , BO , r4
        lsls    r4 , r3 , #5                    // 4 double values
        add     AO , AO , r4
#endif

#if defined(LEFT)
        ldr     r3 , KK
        add     r3 , r3 , #4                    // number of values in AO
        str     r3 , KK
#endif


_L1_M4_END:




_L1_M2_BEGIN:

	tst	I, #2					// I = I % 4
	ble	_L1_M1_BEGIN

_L1_M2_20:

	INIT2x1

#if  (defined(LEFT) &&  defined(TRANSA)) || \
    (!defined(LEFT) && !defined(TRANSA))

	mov	BO, BC
#else
	mov	BO, BC
        ldr     r3 , KK
        lsls    r4 , r3 , #3                                    // 1 double value
        add     BO , BO , r4
        lsls    r4 , r3 , #4                                    // 2 double values
        add     AO , AO , r4

#endif

#ifndef TRMMKERNEL
        ldr     L , K
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
        ldr     L , K
        ldr     r3, KK
        sub     L , L, r3
        str     L , KKK
#else
        ldr     L , KK
#ifdef LEFT
        add     L , L , #2        // number of values in AO
#else
        add     L , L , #1        // number of values in BO
#endif
        str     L , KKK
#endif

	asrs	K1, L, #2					// L = L / 4
	ble	_L1_M2_40

_L1_M2_22:

	KERNEL2x1
	KERNEL2x1
	KERNEL2x1
	KERNEL2x1

	subs	K1, K1, #1
	bgt	_L1_M2_22
	

_L1_M2_40:
	
	ands	K1 , L, #3					// L = L % 4
	ble	_L1_M2_100

_L1_M2_42:

	KERNEL2x1

	subs	K1, K1, #1
	bgt	_L1_M2_42
	
_L1_M2_100:

	SAVE2x1

#if  (defined(LEFT) &&  defined(TRANSA)) || \
    (!defined(LEFT) && !defined(TRANSA))
        ldr     r3 , K
        ldr     r4 , KKK
        sub     r3 , r3 , r4
        lsls    r4 , r3 , #3                    // 1 double value
        add     BO , BO , r4
        lsls    r4 , r3 , #4                    // 2 double values
        add     AO , AO , r4
#endif

#if defined(LEFT)
        ldr     r3 , KK
        add     r3 , r3 , #2                    // number of values in AO
        str     r3 , KK
#endif


_L1_M2_END:



_L1_M1_BEGIN:

	tst	I, #1					// I = I % 4
	ble	_L1_END

_L1_M1_20:

	INIT1x1

#if  (defined(LEFT) &&  defined(TRANSA)) || \
    (!defined(LEFT) && !defined(TRANSA))

	mov	BO, BC
#else
	mov	BO, BC
        ldr     r3 , KK
        lsls    r4 , r3 , #3                                    // 1 double value
        add     BO , BO , r4
        lsls    r4 , r3 , #3                                    // 1 double value
        add     AO , AO , r4

#endif

#ifndef TRMMKERNEL
        ldr     L , K
#elif (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
        ldr     L , K
        ldr     r3, KK
        sub     L , L, r3
        str     L , KKK
#else
        ldr     L , KK
#ifdef LEFT
        add     L , L , #1        // number of values in AO
#else
        add     L , L , #1        // number of values in BO
#endif
        str     L , KKK
#endif

	asrs	K1, L, #2					// L = L / 4
	ble	_L1_M1_40

_L1_M1_22:

	KERNEL1x1
	KERNEL1x1
	KERNEL1x1
	KERNEL1x1

	subs	K1, K1, #1
	bgt	_L1_M1_22
	

_L1_M1_40:
	
	ands	K1 , L, #3					// L = L % 4
	ble	_L1_M1_100

_L1_M1_42:

	KERNEL1x1

	subs	K1, K1, #1
	bgt	_L1_M1_42
	
_L1_M1_100:

	SAVE1x1

#if  (defined(LEFT) &&  defined(TRANSA)) || \
    (!defined(LEFT) && !defined(TRANSA))
        ldr     r3 , K
        ldr     r4 , KKK
        sub     r3 , r3 , r4
        lsls    r4 , r3 , #3                    // 1 double value
        add     BO , BO , r4
        lsls    r4 , r3 , #3                    // 1 double value
        add     AO , AO , r4
#endif

#if defined(LEFT)
        ldr     r3 , KK
        add     r3 , r3 , #1                    // number of values in AO
        str     r3 , KK
#endif



_L1_END:



_L999:

	sub	r3, fp, #128
	vldm	r3, { d8 - d15}					// restore floating point registers

	movs	r0, #0						// set return value
	sub	sp, fp, #24
	pop	{r4 - r9, fp}
	bx	lr

	EPILOGUE

