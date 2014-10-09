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
* 2013/11/23 Saar
* 	 BLASTEST 		: OK
* 	 CTEST			: OK
* 	 TEST			: OK
*
*
* 2013/11/02 Saar
*	UNROLL_N		4
*	UNROLL_M		4
*	DGEMM_P			128
*	DGEMM_Q			240
*	DGEMM_R			12288
*	A_PRE			128
*	B_PRE			128
*	C_PRE			32
*
* Performance on Odroid U2:
*
* 3072x3072		1 Core:		2.62 GFLOPS	ATLAS: 2.69	GFLOPS
* 3072x3072		2 Cores:	5.23 GFLOPS	ATLAS: 5.27	GFLOPS
* 3072x3072		3 Cores:	7.78 GFLOPS	ATLAS: 7.87	GFLOPS
* 3072x3072		4 Cores:       10.10 GFLOPS	ATLAS: 9.98	GFLOPS
**************************************************************************************/

#define ASSEMBLER
#include "common.h"

/*                   X0          X1          X2          s0        X3        x4       x5           x6*/
/*int CNAME(BLASLONG bm,BLASLONG bn,BLASLONG bk,FLOAT alpha,FLOAT* ba,FLOAT* bb,FLOAT* C,BLASLONG ldc*/


#define origM x0
#define origN x1
#define origK x2
#define origPA x3
#define origPB x4
#define pC x5
#define LDC x6
#define offset x7
#define counterL x8
#define counterI x9
#define pB x10
#define counterJ x11
#define tempALPHA x12
#define pCRow0 x13
#define pCRow1 x14
#define pCRow2 x15
#define pA x16

// 00 origM
// 01 origN
// 02 origK
// 03 origPA
// 04 origPB
// 05 pC
// 06 origLDC -> LDC
// 07 offset
// 08 counterL
// 09 counterI
// 10 pB
// 11 counterJ
// 12 tempALPHA      
// 13 pCRow0
// 14 pCRow1
// 15 pCRow2
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

//v00 orig ALPHA -> a00
//v01 a01
//v02 a02
//v03 a03
//v04 a10
//v05 a11
//v06 a12
//v07 a13
//v08 must save b00
//v09 must save b01
//v10 must save b02
//v11 must save b03
//v12 must save b10
//v13 must save b11
//v14 must save b12
//v15 must save b13
//v16 must save  C00
//v17 must save  C01
//v18  C02
//v19  C03
//v20  C10
//v21  C11
//v22  C12
//v23  C13
//v24  C20
//v25  C21
//v26  C22
//v27  C23
//v28  C30
//v29  C31
//v30  C32
//v31  C33

//        add     sp,sp,#-(6*16)
//        stp     x18,x19,[sp,#(0*16)]
//        stp     x20,x21,[sp,#(1*16)]


/**************************************************************************************
* Macro definitions
**************************************************************************************/

.macro INIT4x4

	fsub     	v16.4s , v16.4s , v16.4s
	fsub     	v20.4s , v20.4s , v20.4s
	fsub     	v24.4s , v24.4s , v24.4s
	fsub     	v28.4s , v28.4s , v28.4s

.endm

.macro KERNEL4x4_I

        ld1     {v8.2s},[pB],#8
        ld1     {v10.2s},[pB],#8
        ld1     {v0.4s},[pA],#16

        fmulx   v16.4s, v0.4s, v8.4s[0]
        fmulx   v20.4s, v0.4s, v8.4s[1]
	fmulx	v24.4s, v0.4s, v10.4s[0]
	fmulx	v28.4s, v0.4s, v10.4s[1]

        ld1     {v12.2s},[pB],#8   // for next round
        ld1     {v14.2s},[pB],#8   // for next round
        ld1     {v4.4s},[pA],#16   // for next round


.endm


.macro KERNEL4x4_M2

	fmla  	v16.4s, v4.4s, v12.s[0]
	fmla  	v20.4s, v4.4s, v12.s[1]
	fmla 	v24.4s, v4.4s, v14.s[0]
	fmla  	v28.4s, v4.4s, v14.s[1]

        ld1     {v8.2s},[pB],#8
        ld1     {v10.2s},[pB],#8
        ld1     {v0.4s},[pA],#16

.endm


.macro KERNEL4x4_M1

	fmla 	v16.4s, v0.4s, v8.s[0]
	fmla 	v20.4s, v0.4s, v8.s[1]
	fmla 	v24.4s, v0.4s, v10.s[0]
	fmla 	v28.4s, v0.4s, v10.s[1]

        ld1     {v12.2s},[pB],#8
        ld1     {v14.2s},[pB],#8
        ld1     {v4.4s},[pA],#16

.endm



.macro KERNEL4x4_E

	fmla 	v16.4s, v4.4s, v12.s[0]
	fmla 	v20.4s, v4.4s, v12.s[1]
	fmla 	v24.4s, v4.4s, v14.s[0]
	fmla 	v28.4s, v4.4s, v14.s[1]

.endm




.macro KERNEL4x4_SUB

        ld1     {v8.2s},[pB],#8
        ld1     {v10.2s},[pB],#8
	ld1	{v0.4s} , [pA],#16

	fmla 	v16.4s, v0.4s, v8.s[0]
	fmla 	v20.4s, v0.4s, v8.s[1]
	fmla 	v24.4s, v0.4s, v10.s[0]
	fmla 	v28.4s, v0.4s, v10.s[1]

.endm




.macro SAVE4x4

	add	pCRow1, pCRow0, LDC    // create a second row pointer from the first row pointer
	mov	v0.d[0], tempALPHA

        ld1     {v8.4s},[pCRow0]   // load 4 values of C from first row
        fmla     v8.4s ,v16.4s,v0.s[0]
	st1 	{v8.4s},[pCRow0],#16 // store C from first row

        ld1     {v12.4s},[pCRow1]   // load 4 values of C from second row
        fmla     v12.4s ,v20.4s,v0.s[0]
	st1 	{v12.4s},[pCRow1] // store C from second row

	add	pCRow2, pCRow1, LDC        // Row2 points to third row 

        ld1     {v8.4s},[pCRow2]   // load 4 values of C from third row
        fmla     v8.4s ,v24.4s,v0.s[0]
	st1 	{v8.4s} ,[pCRow2]  // store C from third row

	add	pCRow1, pCRow2 , LDC // row1 points to fourth row

        ld1     {v12.4s},[pCRow1]   // load 4 values of C from fourth row
        fmla     v12.4s ,v28.4s,v0.s[0]
	st1     {v12.4s},[pCRow1]  // store fourth row

.endm

/******************************************************************************/

.macro INIT2x4

	fsub		s16 , s16 , s16
	fmov		s17, s16
	fmov		s20, s16
	fmov		s21, s16
	fmov		s24, s16
	fmov		s25, s16
	fmov		s28, s16
	fmov		s29, s16

.endm



.macro KERNEL2x4_SUB

	ldr	s8 , [ pB ]
	ldr	s9 , [ pB, #4 ]
	ldr	s10, [ pB, #8 ]
	ldr	s11, [ pB, #12 ]

	ldr	s0 , [ pA ]
	ldr	s1 , [ pA, #4 ]

	fmadd 	s16  , s0,  s8,	s16  
	fmadd 	s17  , s1,  s8,	s17  

	fmadd 	s20  , s0,  s9,	s20  
	fmadd 	s21  , s1,  s9,	s21  

	fmadd 	s24  , s0,  s10,	s24  
	fmadd 	s25  , s1,  s10,	s25  

	fmadd 	s28  , s0,  s11,	s28  
	fmadd 	s29  , s1,  s11,	s29  
	add	pA , pA, #8
	add	pB , pB, #16

.endm

            #define F1ST( op1, op2, op3) fmadd op1, op2, op3, op1
            #define L1ST( op1, op2, op3) ldr op1, [op2,  op3]

.macro SAVE2x4

	add	pCRow1 , pCRow0, LDC
	add	pCRow2  , pCRow1, LDC
	mov	v0.d[0], tempALPHA

	L1ST (	s8,pCRow0, #0)
	L1ST (	s9,pCRow0, #4 )

	F1ST ( 	s8 , s0 , s16)
	F1ST ( 	s9 , s0 , s17)

	str 	s8 , [pCRow0, #0]
	str 	s9 , [pCRow0, #4 ]

	ldr	s12, [pCRow1, #0]
	ldr	s13, [pCRow1, #4 ]

	F1ST ( 	s12, s0 , s20)
	F1ST ( 	s13, s0 , s21)

	str 	s12, [pCRow1, #0]
	str 	s13, [pCRow1, #4 ]

	L1ST (	s8,pCRow2 , #0)
	L1ST (	s9,pCRow2 , #4 )

	F1ST ( 	s8 , s0 , s24)
	F1ST ( 	s9 , s0 , s25)

	str 	s8 , [pCRow2 , #0]
	str 	s9 , [pCRow2 , #4 ]

	add	pCRow1, pCRow2 , LDC

	ldr	s12, [pCRow1, #0]
	ldr	s13, [pCRow1, #4 ]

	F1ST ( 	s12, s0 , s28)
	F1ST ( 	s13, s0 , s29)

	str 	s12, [pCRow1, #0]
	str 	s13, [pCRow1, #4 ]

	add	pCRow0, pCRow0, #8

.endm


/******************************************************************************/

.macro INIT1x4

	fsub		s16 , s16 , s16
	fmov		s20, s16
	fmov		s24, s16
	fmov		s28, s16

.endm



.macro KERNEL1x4_SUB

	ldr	s8 , [ pB ]
	ldr	s9 , [ pB, #4 ]
	ldr	s10, [ pB, #8 ]
	ldr	s11, [ pB, #12 ]

	ldr	s0 , [ pA ]

	fmadd 	s16  , s0,  s8,	s16  
	fmadd 	s20  , s0,  s9,	s20  
	fmadd 	s24  , s0,  s10,	s24  
	fmadd 	s28  , s0,  s11,	s28  

	add	pA , pA, #4
	add	pB , pB, #16

.endm

.macro SAVE1x4

	add	pCRow1 , pCRow0, LDC
	add	pCRow2  , pCRow1, LDC

	mov	v0.d[0], tempALPHA

	L1ST (	s8,pCRow0, #0)
	F1ST ( 	s8 , s0 , s16)
	str 	s8 , [pCRow0, #0]

	L1ST (	s12,pCRow1, #0)
	F1ST ( 	s12, s0 , s20)
	str 	s12, [pCRow1, #0]

	L1ST (	s8,pCRow2 , #0)
	F1ST ( 	s8 , s0 , s24)
	str 	s8 , [pCRow2 , #0]

	add	pCRow1, pCRow2 , LDC

	L1ST (	s12,pCRow1, #0)
	F1ST ( 	s12, s0 , s28)
	str 	s12, [pCRow1, #0]

	add	pCRow0, pCRow0, #4

.endm

/******************************************************************************/
/******************************************************************************/

.macro INIT4x2

	fsub		s16 , s16 , s16
	fmov		s17, s16
	fmov		s18, s16
	fmov		s19, s16
	fmov		s20, s16
	fmov		s21, s16
	fmov		s22, s16
	fmov		s23, s16

.endm



.macro KERNEL4x2_SUB

	ldr	s8 , [ pB ]
	ldr	s9 , [ pB, #4 ]

	ldr	s0 , [ pA ]
	ldr	s1 , [ pA, #4 ]
	ldr	s2 , [ pA, #8 ]
	ldr	s3 , [ pA, #12 ]

	fmadd 	s16  , s0,  s8,	s16  
	fmadd 	s17  , s1,  s8,	s17  
	fmadd 	s18  , s2,  s8,	s18  
	fmadd 	s19  , s3,  s8,	s19  

	fmadd 	s20  , s0,  s9,	s20  
	fmadd 	s21  , s1,  s9,	s21  
	fmadd 	s22  , s2,  s9,	s22  
	fmadd 	s23  , s3,  s9,	s23  

	add	pA , pA, #16
	add	pB , pB, #8

.endm

.macro SAVE4x2

	add	pCRow1 , pCRow0, LDC

	mov	v0.d[0], tempALPHA

	L1ST (	s8,pCRow0, #0)
	L1ST (	s9,pCRow0, #4 )
	L1ST (	s10,pCRow0, #8 )
	L1ST (	s11,pCRow0, #12 )

	F1ST ( 	s8 , s0 , s16)
	F1ST ( 	s9 , s0 , s17)
	F1ST ( 	s10, s0 , s18)
	F1ST ( 	s11, s0 , s19)

	str 	s8 , [pCRow0]
	str 	s9 , [pCRow0, #4 ]
	str 	s10, [pCRow0, #8 ]
	str 	s11, [pCRow0, #12 ]

	L1ST (	s12,pCRow1, #0)
	L1ST (	s13,pCRow1, #4 )
	L1ST (	s14,pCRow1, #8 )
	L1ST (	s15,pCRow1, #12 )

	F1ST ( 	s12, s0 , s20)
	F1ST ( 	s13, s0 , s21)
	F1ST ( 	s14, s0 , s22)
	F1ST ( 	s15, s0 , s23)

	str 	s12, [pCRow1]
	str 	s13, [pCRow1, #4 ]
	str 	s14, [pCRow1, #8 ]
	str 	s15, [pCRow1, #12 ]

	add	pCRow0, pCRow0, #16

.endm


/******************************************************************************/

.macro INIT2x2

	fsub		s16 , s16 , s16
	fmov		s17, s16
	fmov		s20, s16
	fmov		s21, s16

.endm



.macro KERNEL2x2_SUB

	ldr	s8 , [ pB ]
	ldr	s9 , [ pB, #4 ]

	ldr	s0 , [ pA ]
	ldr	s1 , [ pA, #4 ]

	fmadd 	s16  , s0,  s8,	s16  
	fmadd 	s17  , s1,  s8,	s17  

	fmadd 	s20  , s0,  s9,	s20  
	fmadd 	s21  , s1,  s9,	s21  

	add	pA , pA, #8
	add	pB , pB, #8

.endm

.macro SAVE2x2

	add	pCRow1 , pCRow0, LDC

	mov	v0.d[0], tempALPHA

	L1ST (	s8,pCRow0, #0 )
	L1ST (	s9,pCRow0, #4 )

	F1ST ( 	s8 , s0 , s16)
	F1ST ( 	s9 , s0 , s17)

	str 	s8 , [pCRow0]
	str 	s9 , [pCRow0, #4 ]

	L1ST (	s12,pCRow1, #0 )
	L1ST (	s13,pCRow1, #4 )

	F1ST ( 	s12, s0 , s20)
	F1ST ( 	s13, s0 , s21)

	str 	s12, [pCRow1]
	str 	s13, [pCRow1, #4 ]

	add	pCRow0, pCRow0, #8

.endm

/******************************************************************************/

.macro INIT1x2

	fsub		s16 , s16 , s16
	fmov		s20, s16

.endm



.macro KERNEL1x2_SUB

	ldr	s8 , [ pB ]
	ldr	s9 , [ pB, #4 ]

	ldr	s0 , [ pA ]
	fmadd 	s16  , s0,  s8,	s16  
	fmadd 	s20  , s0,  s9,	s20  

	add	pA , pA, #4
	add	pB , pB, #8

.endm

.macro SAVE1x2

	add	pCRow1 , pCRow0, LDC

	mov	v0.d[0], tempALPHA

	L1ST (	s8,pCRow0, #0)
	F1ST ( 	s8 , s0 , s16)
	str 	s8 , [pCRow0]

	L1ST (	s12,pCRow1, #0)
	F1ST ( 	s12, s0 , s20)
	str 	s12, [pCRow1]

	add	pCRow0, pCRow0, #4

.endm

/******************************************************************************/
/******************************************************************************/

.macro INIT4x1

	fsub		s16 , s16 , s16
	fmov		s17, s16
	fmov		s18, s16
	fmov		s19, s16

.endm



.macro KERNEL4x1_SUB

	ldr	s8 , [ pB ]

	ldr	s0 , [ pA ]
	ldr	s1 , [ pA, #4 ]
	ldr	s2 , [ pA, #8 ]
	ldr	s3 , [ pA, #12 ]

	fmadd 	s16  , s0,  s8,	s16  
	fmadd 	s17  , s1,  s8,	s17  
	fmadd 	s18  , s2,  s8,	s18  
	fmadd 	s19  , s3,  s8,	s19  

	add	pA , pA, #16
	add	pB , pB, #4

.endm

.macro SAVE4x1


	mov	v0.d[0], tempALPHA

	L1ST (	s8,pCRow0, #0 )
	L1ST (	s9,pCRow0, #4 )
	L1ST (	s10,pCRow0, #8 )
	L1ST (	s11,pCRow0, #12 )

	F1ST ( 	s8 , s0 , s16)
	F1ST ( 	s9 , s0 , s17)
	F1ST ( 	s10, s0 , s18)
	F1ST ( 	s11, s0 , s19)

	str 	s8 , [pCRow0]
	str 	s9 , [pCRow0, #4 ]
	str 	s10, [pCRow0, #8 ]
	str 	s11, [pCRow0, #12 ]

	add	pCRow0, pCRow0, #16

.endm




/******************************************************************************/

.macro INIT2x1

	fsub		s16 , s16 , s16
	fmov		s17, s16

.endm



.macro KERNEL2x1_SUB

	ldr	s8 , [ pB ]

	ldr	s0 , [ pA ]
	ldr	s1 , [ pA, #4 ]

	fmadd 	s16  , s0,  s8,	s16  
	fmadd 	s17  , s1,  s8,	s17  

	add	pA , pA, #8
	add	pB , pB, #4

.endm

.macro SAVE2x1


	mov	v0.d[0], tempALPHA

	L1ST (	s8,pCRow0, #0 )
	L1ST (	s9,pCRow0, #4 )

	F1ST ( 	s8 , s0 , s16)
	F1ST ( 	s9 , s0 , s17)

	str 	s8 , [pCRow0]
	str 	s9 , [pCRow0, #4 ]

	add	pCRow0, pCRow0, #8

.endm

/******************************************************************************/

.macro INIT1x1

	fsub		s16 , s16 , s16

.endm



.macro KERNEL1x1_SUB

	ldr	s8 , [ pB ]

	ldr	s0 , [ pA ]

	fmadd 	s16  , s0,  s8,	s16  

	add	pA , pA, #4
	add	pB , pB, #4

.endm

.macro SAVE1x1


	mov	v0.d[0], tempALPHA

	L1ST (	s8,pCRow0, #0 )
	F1ST ( 	s8 , s0 , s16)
	str 	s8 , [pCRow0]

	add	pCRow0, pCRow0, #4

.endm





/**************************************************************************************
* End of macro definitions
**************************************************************************************/

	PROLOGUE

	.align 5
        add     sp,sp,#-(5*16)
        stp     d8,d9,[sp,#(0*16)]
        stp     d10,d11,[sp,#(1*16)]
        stp     d12,d13,[sp,#(2*16)]
        stp     d14,d15,[sp,#(3*16)]
        stp     d16,d17,[sp,#(4*16)]

        mov     tempALPHA, v0.d[0]
	lsl	LDC, LDC, #2					// ldc = ldc * 4

	mov	pB, origPB

	mov	counterJ, origN
	asr 	counterJ, counterJ, #2					// J = J / 4
	cmp 	counterJ, #0
	ble	sgemm_kernel_L2_BEGIN

sgemm_kernel_L4_BEGIN:

	mov	pCRow0, pC						// pCRow0 = C
        add     pC,pC,LDC, lsl #2

	mov	pA, origPA						// pA = start of A array



sgemm_kernel_L4_M4_BEGIN:

	mov	counterI, origM
	asr 	counterI, counterI, #2					// counterI = counterI / 4
	cmp 	counterI, #0
	ble	sgemm_kernel_L4_M2_BEGIN

sgemm_kernel_L4_M4_20:

	mov	pB, origPB
	asr 	counterL , origK, #1					// L = K / 2
	cmp	counterL , #2                                           // is there at least 4 to do?
	blt	sgemm_kernel_L4_M4_32



	KERNEL4x4_I     //do one in the K
	KERNEL4x4_M2    //do another in the K

	subs	counterL, counterL, #2  // subtract 2, since one is always done at the tail
	ble	sgemm_kernel_L4_M4_22a
	.align 5

sgemm_kernel_L4_M4_22:

	KERNEL4x4_M1
	KERNEL4x4_M2

	subs	counterL, counterL, #1
	bgt	sgemm_kernel_L4_M4_22

sgemm_kernel_L4_M4_22a:

	KERNEL4x4_M1
	KERNEL4x4_E

	b	 sgemm_kernel_L4_M4_44

sgemm_kernel_L4_M4_32:   // less than 4 to do in the K direction

	tst	counterL, #1
	ble	sgemm_kernel_L4_M4_40

	KERNEL4x4_I

	KERNEL4x4_E

	b	 sgemm_kernel_L4_M4_44


sgemm_kernel_L4_M4_40:

	INIT4x4


sgemm_kernel_L4_M4_44:

	ands	counterL , origK, #1
	ble	sgemm_kernel_L4_M4_100

sgemm_kernel_L4_M4_46:

	KERNEL4x4_SUB

	subs	counterL, counterL, #1
	bne	sgemm_kernel_L4_M4_46

sgemm_kernel_L4_M4_100:

	SAVE4x4

sgemm_kernel_L4_M4_END:

	subs	counterI, counterI, #1
	bne	sgemm_kernel_L4_M4_20


sgemm_kernel_L4_M2_BEGIN:

	mov	counterI, origM
	tst	counterI , #3
	ble	sgemm_kernel_L4_END

	tst	counterI, #2					// counterI = counterI / 2
	ble	sgemm_kernel_L4_M1_BEGIN

sgemm_kernel_L4_M2_20:

	INIT2x4

	mov	pB, origPB
	asr 	counterL , origK, #3					// counterL = counterL / 8
	cmp	counterL , #0
	ble	sgemm_kernel_L4_M2_40

sgemm_kernel_L4_M2_22:

	KERNEL2x4_SUB
	KERNEL2x4_SUB
	KERNEL2x4_SUB
	KERNEL2x4_SUB

	KERNEL2x4_SUB
	KERNEL2x4_SUB
	KERNEL2x4_SUB
	KERNEL2x4_SUB

	subs	counterL, counterL, #1
	bgt	sgemm_kernel_L4_M2_22


sgemm_kernel_L4_M2_40:

	ands	counterL , origK, #7					// counterL = counterL % 8
	ble	sgemm_kernel_L4_M2_100

sgemm_kernel_L4_M2_42:

	KERNEL2x4_SUB

	subs	counterL, counterL, #1
	bgt	sgemm_kernel_L4_M2_42

sgemm_kernel_L4_M2_100:

	SAVE2x4

sgemm_kernel_L4_M2_END:


sgemm_kernel_L4_M1_BEGIN:

	tst	counterI, #1					// counterI = counterI % 2
	ble	sgemm_kernel_L4_END

sgemm_kernel_L4_M1_20:

	INIT1x4

	mov	pB, origPB
	asr 	counterL , origK, #3					// counterL = counterL / 8
	cmp	counterL , #0
	ble	sgemm_kernel_L4_M1_40

sgemm_kernel_L4_M1_22:
	KERNEL1x4_SUB
	KERNEL1x4_SUB
	KERNEL1x4_SUB
	KERNEL1x4_SUB

	KERNEL1x4_SUB
	KERNEL1x4_SUB
	KERNEL1x4_SUB
	KERNEL1x4_SUB

	subs	counterL, counterL, #1
	bgt	sgemm_kernel_L4_M1_22


sgemm_kernel_L4_M1_40:

	ands	counterL , origK, #7					// counterL = counterL % 8
	ble	sgemm_kernel_L4_M1_100

sgemm_kernel_L4_M1_42:

	KERNEL1x4_SUB

	subs	counterL, counterL, #1
	bgt	sgemm_kernel_L4_M1_42

sgemm_kernel_L4_M1_100:

	SAVE1x4


sgemm_kernel_L4_END:

	add	origPB, origPB, origK, lsl #4 // B = B + K * 4 * 4

	subs	counterJ, counterJ , #1						// j--
	bgt	sgemm_kernel_L4_BEGIN



/*********************************************************************************************/

sgemm_kernel_L2_BEGIN:   // less than 2 left in N direction

	mov	counterJ , origN
	tst	counterJ , #3
	ble	sgemm_kernel_L999   // error, N was less than 4?

	tst	counterJ , #2
	ble	sgemm_kernel_L1_BEGIN

	mov	pCRow0, pC						// pCRow0 = pC
	add	pC , pC, LDC, lsl #1

	mov	pA, origPA						// pA = A



sgemm_kernel_L2_M4_BEGIN:

	mov	counterI, origM
	asr 	counterI, counterI, #2					// counterI = counterI / 4
	cmp	counterI,#0
	ble	sgemm_kernel_L2_M2_BEGIN

sgemm_kernel_L2_M4_20:

	INIT4x2

	mov	pB, origPB
	asr	counterL , origK, #3					// counterL = counterL / 8
	cmp	counterL,#0
	ble	sgemm_kernel_L2_M4_40
	.align 5

sgemm_kernel_L2_M4_22:
	KERNEL4x2_SUB
	KERNEL4x2_SUB
	KERNEL4x2_SUB
	KERNEL4x2_SUB

	KERNEL4x2_SUB
	KERNEL4x2_SUB
	KERNEL4x2_SUB
	KERNEL4x2_SUB

	subs	counterL, counterL, #1
	bgt	sgemm_kernel_L2_M4_22


sgemm_kernel_L2_M4_40:

	ands	counterL , origK, #7					// counterL = counterL % 8
	ble	sgemm_kernel_L2_M4_100

sgemm_kernel_L2_M4_42:

	KERNEL4x2_SUB

	subs	counterL, counterL, #1
	bgt	sgemm_kernel_L2_M4_42

sgemm_kernel_L2_M4_100:

	SAVE4x2

sgemm_kernel_L2_M4_END:

	subs	counterI, counterI, #1
	bgt	sgemm_kernel_L2_M4_20


sgemm_kernel_L2_M2_BEGIN:

	mov	counterI, origM
	tst	counterI , #3
	ble	sgemm_kernel_L2_END

	tst	counterI, #2					// counterI = counterI / 2
	ble	sgemm_kernel_L2_M1_BEGIN

sgemm_kernel_L2_M2_20:

	INIT2x2

	mov	pB, origPB
	asr	counterL , origK, #3					// counterL = counterL / 8
        cmp	counterL,#0
	ble	sgemm_kernel_L2_M2_40

sgemm_kernel_L2_M2_22:

	KERNEL2x2_SUB
	KERNEL2x2_SUB
	KERNEL2x2_SUB
	KERNEL2x2_SUB

	KERNEL2x2_SUB
	KERNEL2x2_SUB
	KERNEL2x2_SUB
	KERNEL2x2_SUB

	subs	counterL, counterL, #1
	bgt	sgemm_kernel_L2_M2_22


sgemm_kernel_L2_M2_40:

	ands	counterL , origK, #7					// counterL = counterL % 8
	ble	sgemm_kernel_L2_M2_100

sgemm_kernel_L2_M2_42:

	KERNEL2x2_SUB

	subs	counterL, counterL, #1
	bgt	sgemm_kernel_L2_M2_42

sgemm_kernel_L2_M2_100:

	SAVE2x2

sgemm_kernel_L2_M2_END:


sgemm_kernel_L2_M1_BEGIN:

	tst	counterI, #1					// counterI = counterI % 2
	ble	sgemm_kernel_L2_END

sgemm_kernel_L2_M1_20:

	INIT1x2

	mov	pB, origPB
	asr 	counterL , origK, #3					// counterL = counterL / 8
        cmp     counterL, #0
	ble	sgemm_kernel_L2_M1_40

sgemm_kernel_L2_M1_22:
	KERNEL1x2_SUB
	KERNEL1x2_SUB
	KERNEL1x2_SUB
	KERNEL1x2_SUB

	KERNEL1x2_SUB
	KERNEL1x2_SUB
	KERNEL1x2_SUB
	KERNEL1x2_SUB

	subs	counterL, counterL, #1
	bgt	sgemm_kernel_L2_M1_22


sgemm_kernel_L2_M1_40:

	ands	counterL , origK, #7					// counterL = counterL % 8
	ble	sgemm_kernel_L2_M1_100

sgemm_kernel_L2_M1_42:

	KERNEL1x2_SUB

	subs	counterL, counterL, #1
	bgt	sgemm_kernel_L2_M1_42

sgemm_kernel_L2_M1_100:

	SAVE1x2


sgemm_kernel_L2_END:
	add	origPB, origPB, origK, lsl #3					// B = B + K * 2 * 4

/*********************************************************************************************/

sgemm_kernel_L1_BEGIN:

	mov	counterJ , origN
	tst	counterJ , #1
	ble	sgemm_kernel_L999 // done


	mov	pCRow0, pC						// pCRow0 = C
	add	pC , pCRow0 , LDC                                 // C01 is the current line, update pC to point to next

	mov	pA, origPA						// pA = A



sgemm_kernel_L1_M4_BEGIN:

	mov	counterI, origM
	asr 	counterI, counterI, #2					// counterI = counterI / 4
	cmp	counterI, #0
	ble	sgemm_kernel_L1_M2_BEGIN

sgemm_kernel_L1_M4_20:

	INIT4x1

	mov	pB, origPB
	asr	counterL , origK, #3					// counterL = counterL / 8
	cmp	counterL , #0
	ble	sgemm_kernel_L1_M4_40
	.align 5

sgemm_kernel_L1_M4_22:
	KERNEL4x1_SUB
	KERNEL4x1_SUB
	KERNEL4x1_SUB
	KERNEL4x1_SUB

	KERNEL4x1_SUB
	KERNEL4x1_SUB
	KERNEL4x1_SUB
	KERNEL4x1_SUB

	subs	counterL, counterL, #1
	bgt	sgemm_kernel_L1_M4_22


sgemm_kernel_L1_M4_40:

	ands	counterL , origK, #7					// counterL = counterL % 8
	ble	sgemm_kernel_L1_M4_100

sgemm_kernel_L1_M4_42:

	KERNEL4x1_SUB

	subs	counterL, counterL, #1
	bgt	sgemm_kernel_L1_M4_42

sgemm_kernel_L1_M4_100:

	SAVE4x1

sgemm_kernel_L1_M4_END:

	subs	counterI, counterI, #1
	bgt	sgemm_kernel_L1_M4_20


sgemm_kernel_L1_M2_BEGIN:

	mov	counterI, origM
	tst	counterI , #3
	ble	sgemm_kernel_L1_END

	tst	counterI, #2					// counterI = counterI / 2
	ble	sgemm_kernel_L1_M1_BEGIN

sgemm_kernel_L1_M2_20:

	INIT2x1

	mov	pB, origPB
	asr 	counterL , origK, #3					// counterL = counterL / 8
	cmp	counterL , #0
	ble	sgemm_kernel_L1_M2_40

sgemm_kernel_L1_M2_22:

	KERNEL2x1_SUB
	KERNEL2x1_SUB
	KERNEL2x1_SUB
	KERNEL2x1_SUB

	KERNEL2x1_SUB
	KERNEL2x1_SUB
	KERNEL2x1_SUB
	KERNEL2x1_SUB

	subs	counterL, counterL, #1
	bgt	sgemm_kernel_L1_M2_22


sgemm_kernel_L1_M2_40:

	ands	counterL , origK, #7					// counterL = counterL % 8
	ble	sgemm_kernel_L1_M2_100

sgemm_kernel_L1_M2_42:

	KERNEL2x1_SUB

	subs	counterL, counterL, #1
	bgt	sgemm_kernel_L1_M2_42

sgemm_kernel_L1_M2_100:

	SAVE2x1

sgemm_kernel_L1_M2_END:


sgemm_kernel_L1_M1_BEGIN:

	tst	counterI, #1					// counterI = counterI % 2
	ble	sgemm_kernel_L1_END

sgemm_kernel_L1_M1_20:

	INIT1x1

	mov	pB, origPB
	asr 	counterL , origK, #3					// counterL = counterL / 8
	cmp	counterL , #0
	ble	sgemm_kernel_L1_M1_40

sgemm_kernel_L1_M1_22:
	KERNEL1x1_SUB
	KERNEL1x1_SUB
	KERNEL1x1_SUB
	KERNEL1x1_SUB

	KERNEL1x1_SUB
	KERNEL1x1_SUB
	KERNEL1x1_SUB
	KERNEL1x1_SUB

	subs	counterL, counterL, #1
	bgt	sgemm_kernel_L1_M1_22


sgemm_kernel_L1_M1_40:

	ands	counterL , origK, #7					// counterL = counterL % 8
	ble	sgemm_kernel_L1_M1_100

sgemm_kernel_L1_M1_42:

	KERNEL1x1_SUB

	subs	counterL, counterL, #1
	bgt	sgemm_kernel_L1_M1_42

sgemm_kernel_L1_M1_100:

	SAVE1x1


sgemm_kernel_L1_END:


sgemm_kernel_L999:
	mov	x0, #0						// set return value
        ldp     d8,d9,[sp,#(0*16)]
        ldp     d10,d11,[sp,#(1*16)]
        ldp     d12,d13,[sp,#(2*16)]
        ldp     d14,d15,[sp,#(3*16)]
        ldp     d16,d17,[sp,#(4*16)]
        add     sp,sp,#(5*16)
	ret

	EPILOGUE

