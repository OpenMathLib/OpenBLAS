/***************************************************************************
Copyright (c) 2013-2019, The OpenBLAS Project
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
#define MY_ALIGN .align 3

	srawi.		J,	N,	1
	ble		ZGEMM_L2_END

ZGEMM_L2_BEGIN:

	mr		BO,	B
	mr		BBO,	BBUFFER
	srawi.		T1,	K,	2
	ble		ZGEMM_L2_COPYB1

ZGEMM_L2_COPYB8:

	addi		T2,	PRE, 128
	dcbt		BO,	PRE
	dcbtst		BBO,	PRE
	dcbtst		BBO,	T2
	ZCOPYB_8
	addic.		T1,	T1,	-1

	bgt		ZGEMM_L2_COPYB8

ZGEMM_L2_COPYB1:

	andi.		T1,	K,	3
	ble		ZGEMM_L2_COPYB_END

ZGEMM_L2_COPYB_LOOP:

	ZCOPYB_2
	addic.          T1,     T1,     -1

	bgt             ZGEMM_L2_COPYB_LOOP

ZGEMM_L2_COPYB_END:

	mr		CO,	C
	mr		AO,	A
	slwi		T1,	LDC	,	1
	add		C,	C,	T1
	srawi.		I,	M,	3
	ble		ZGEMM_L2x8_END

ZGEMM_L2x8_BEGIN:


	mr		BO,	BBUFFER
	mr T1, K
    addi T1,T1, -1
    srawi.		L,	T1,	5 /**(K-1) % 32x */ 
	ZERO2x8  
	ble		ZGEMM_L2x8_SUB0
 

ZGEMM_L2x8_LOOP_START:

    LOAD2x8 0 
    li T2, 1024
	li T3, 1024+512
	li T4, 2048
	li T5, 2048+512
	mtctr		L

	MY_ALIGN
ZGEMM_L2x8_LOOP:
 	dcbt		AO,	PRE
	dcbt		BO,	PRE
    KERNEL2x8_L 128,64,0,0
	KERNEL2x8_L 128,64,1,0
	dcbt		AO,	T2	
	KERNEL2x8_L 128,64,2,0
	KERNEL2x8_L 128,64,3,0 
	dcbt		AO,	T3
	dcbt		BO,	T2
    KERNEL2x8_L 128,64,4,0
	KERNEL2x8_L 128,64,5,0
	dcbt		AO,	T4	
	KERNEL2x8_L 128,64,6,0
	KERNEL2x8_L 128,64,7,0  
	dcbt		AO,	T5	
	dcbt		BO,	T3
    KERNEL2x8_L 128,64,8,0
	KERNEL2x8_L 128,64,9,0
	KERNEL2x8_L 128,64,10,0
	KERNEL2x8_L 128,64,11,0  
	dcbt		BO,	T4
    KERNEL2x8_L 128,64,12,0
	KERNEL2x8_L 128,64,13,0
	KERNEL2x8_L 128,64,14,0
	KERNEL2x8_L 128,64,15,1 		
	bdnz		ZGEMM_L2x8_LOOP
 	MY_ALIGN  
ZGEMM_L2x8_LOOP_END:
    END2x8  AO, BO, 128, 64   	 
 
	b		ZGEMM_L2x8_SUB1
 
ZGEMM_L2x8_SUB0:

	andi.		L,	K,	63
 
	b		ZGEMM_L2x8_SUB2

ZGEMM_L2x8_SUB1:

	andi.		L,	T1,	31
	ble		ZGEMM_L2x8_SAVE

ZGEMM_L2x8_SUB2:
    srawi.      T1,L, 3
    ble ZGEMM_L2x8_SUB2_4
    mtctr		T1
    MY_ALIGN
ZGEMM_L2x8_SUB2_LOOP:
	LOAD2x8 0 
    KERNEL2x8_L  128,64, 0,0
    KERNEL2x8_L  128,64, 1,0
    KERNEL2x8_L  128,64, 2,0
    KERNEL2x8_E  128,64, 3,1
    bdnz ZGEMM_L2x8_SUB2_LOOP 
    MY_ALIGN  
ZGEMM_L2x8_SUB2_4:
    andi.      T1,L, 4
    ble ZGEMM_L2x8_SUB2_2
	LOAD2x8 0 
    KERNEL2x8_L  128,64, 0,0
    KERNEL2x8_E  128,64, 1,1
    MY_ALIGN
ZGEMM_L2x8_SUB2_2:
    andi.      T1,L, 2
    ble ZGEMM_L2x8_SUB2_1
	LOAD2x8 0 
    KERNEL2x8_E  128,64, 0,1
    MY_ALIGN    
ZGEMM_L2x8_SUB2_1:
    andi.      T1,L, 1
    ble ZGEMM_L2x8_SAVE	
    KERNEL2x8      

/*	addic.		L,	L,	-1
	bgt		ZGEMM_L2x8_SUB2_1*/

ZGEMM_L2x8_SAVE:

	SAVE2x8

	addic.		I,	I,	-1
	bgt		ZGEMM_L2x8_BEGIN

ZGEMM_L2x8_END:

ZGEMM_L2x4_BEGIN:

	andi.		T2,	M,	7
	ble		ZGEMM_L2x1_END

	andi.		T1,	M,	4
	ble		ZGEMM_L2x4_END
	mr		BO,	BBUFFER
	mr T1, K
    addi T1,T1, -1
    srawi.		L,	T1, 4 /**(K-1) % 16x */ 
	ZERO2x4  
	ble		ZGEMM_L2x4_SUB0 

ZGEMM_L2x4_LOOP_START:
    LOAD2x4 0  
	mtctr		L

	MY_ALIGN
ZGEMM_L2x4_LOOP: 
    KERNEL2x4_L 64,64,0,0
	KERNEL2x4_L 64,64,1,0 	
	KERNEL2x4_L 64,64,2,0
	KERNEL2x4_L 64,64,3,0  
    KERNEL2x4_L 64,64,4,0
	KERNEL2x4_L 64,64,5,0 
	KERNEL2x4_L 64,64,6,0
	KERNEL2x4_L 64,64,7,1	
	bdnz		ZGEMM_L2x4_LOOP
 	MY_ALIGN  
ZGEMM_L2x4_LOOP_END:
    END2x4  AO, BO, 64, 64   	 
 
	b		ZGEMM_L2x4_SUB1
 
ZGEMM_L2x4_SUB0:

	andi.		L,	K,	31
 
	b		ZGEMM_L2x4_SUB2

ZGEMM_L2x4_SUB1:

	andi.		L,	T1,	15
	ble		ZGEMM_L2x4_SAVE

ZGEMM_L2x4_SUB2:
    srawi.      T1,L, 3
    ble ZGEMM_L2x4_SUB2_4
    mtctr		T1
    MY_ALIGN
ZGEMM_L2x4_SUB2_LOOP:
	LOAD2x4 0 
    KERNEL2x4_L  64,64, 0,0
    KERNEL2x4_L  64,64, 1,0
    KERNEL2x4_L  64,64, 2,0
    KERNEL2x4_E  64,64, 3,1
    bdnz ZGEMM_L2x4_SUB2_LOOP 
    MY_ALIGN  
ZGEMM_L2x4_SUB2_4:
    andi.      T1,L, 4
    ble ZGEMM_L2x4_SUB2_2
	LOAD2x4 0 
    KERNEL2x4_L  64,64, 0,0
    KERNEL2x4_E  64,64, 1,1
    MY_ALIGN
ZGEMM_L2x4_SUB2_2:
    andi.      T1,L, 2
    ble ZGEMM_L2x4_SUB2_1
	LOAD2x4 0 
    KERNEL2x4_E  64,64, 0,1
    MY_ALIGN    
ZGEMM_L2x4_SUB2_1:
    andi.      T1,L, 1
    ble ZGEMM_L2x4_SAVE	
    KERNEL2x4  

ZGEMM_L2x4_SAVE:

	SAVE2x4

ZGEMM_L2x4_END:

ZGEMM_L2x2_BEGIN:


	andi.		T1,	M,	2
	ble		ZGEMM_L2x2_END
	mr		BO,	BBUFFER
	mr T1, K
    addi T1,T1, -1
    srawi.		L,	T1, 4 /**(K-1) % 16x */ 
	ZERO2x2 
	ble		ZGEMM_L2x2_SUB0 

ZGEMM_L2x2_LOOP_START:
    LOAD2x2 0  
	mtctr		L

	MY_ALIGN
ZGEMM_L2x2_LOOP: 
    KERNEL2x2_L 32,64,0,0
	KERNEL2x2_L 32,64,1,0 	
	KERNEL2x2_L 32,64,2,0
	KERNEL2x2_L 32,64,3,0  
    KERNEL2x2_L 32,64,4,0
	KERNEL2x2_L 32,64,5,0 
	KERNEL2x2_L 32,64,6,0
	KERNEL2x2_L 32,64,7,1	
	bdnz		ZGEMM_L2x2_LOOP
 	MY_ALIGN  
ZGEMM_L2x2_LOOP_END:
    END2x2  AO, BO, 32, 64   	 
 
	b		ZGEMM_L2x2_SUB1
 
ZGEMM_L2x2_SUB0:

	andi.		L,	K,	31
 
	b		ZGEMM_L2x2_SUB2

ZGEMM_L2x2_SUB1:

	andi.		L,	T1,	15
	ble		ZGEMM_L2x2_SAVE

ZGEMM_L2x2_SUB2:
    srawi.      T1,L, 3
    ble ZGEMM_L2x2_SUB2_4
    mtctr		T1
    MY_ALIGN
ZGEMM_L2x2_SUB2_LOOP:
	LOAD2x2 0 
    KERNEL2x2_L  32,64, 0,0
    KERNEL2x2_L  32,64, 1,0
    KERNEL2x2_L  32,64, 2,0
    KERNEL2x2_E  32,64, 3,1
    bdnz ZGEMM_L2x2_SUB2_LOOP 
    MY_ALIGN  
ZGEMM_L2x2_SUB2_4:
    andi.      T1,L, 4
    ble ZGEMM_L2x2_SUB2_2
	LOAD2x2 0 
    KERNEL2x2_L  32,64, 0,0
    KERNEL2x2_E  32,64, 1,1
    MY_ALIGN
ZGEMM_L2x2_SUB2_2:
    andi.      T1,L, 2
    ble ZGEMM_L2x2_SUB2_1
	LOAD2x2 0 
    KERNEL2x2_E  32,64, 0,1
    MY_ALIGN    
ZGEMM_L2x2_SUB2_1:
    andi.      T1,L, 1
    ble ZGEMM_L2x2_SAVE	
    KERNEL2x2 
ZGEMM_L2x2_SAVE:

	SAVE2x2

ZGEMM_L2x2_END:

ZGEMM_L2x1_BEGIN:


	andi.		T1,	M,	1
	ble		ZGEMM_L2x1_END
	mr		BO,	BBUFFER
	mr T1, K
    addi T1,T1, -1
    srawi.		L,	T1, 4 /**(K-1) % 16x */ 
	ZERO2x1  
	ble		ZGEMM_L2x1_SUB0 

ZGEMM_L2x1_LOOP_START:

    LOAD2x1 0  
	mtctr		L

	MY_ALIGN
ZGEMM_L2x1_LOOP: 
    KERNEL2x1_L 16,64,0,0
	KERNEL2x1_L 16,64,1,0 	
	KERNEL2x1_L 16,64,2,0
	KERNEL2x1_L 16,64,3,0  
    KERNEL2x1_L 16,64,4,0
	KERNEL2x1_L 16,64,5,0 
	KERNEL2x1_L 16,64,6,0
	KERNEL2x1_L 16,64,7,1 		
	bdnz		ZGEMM_L2x1_LOOP
 	MY_ALIGN  
ZGEMM_L2x1_LOOP_END:
    END2x1  AO, BO, 16, 64   	 
 
	b		ZGEMM_L2x1_SUB1
 
ZGEMM_L2x1_SUB0:

	andi.		L,	K,	31
 
	b		ZGEMM_L2x1_SUB2

ZGEMM_L2x1_SUB1:

	andi.		L,	T1,	15
	ble		ZGEMM_L2x1_SAVE

ZGEMM_L2x1_SUB2:
    srawi.      T1,L, 3
    ble ZGEMM_L2x1_SUB2_4
    mtctr		T1
    MY_ALIGN
ZGEMM_L2x1_SUB2_LOOP:
	LOAD2x1 0 
    KERNEL2x1_L  16,64, 0,0
    KERNEL2x1_L  16,64, 1,0
    KERNEL2x1_L  16,64, 2,0
    KERNEL2x1_E  16,64, 3,1
    bdnz ZGEMM_L2x1_SUB2_LOOP 
    MY_ALIGN  
ZGEMM_L2x1_SUB2_4:
    andi.      T1,L, 4
    ble ZGEMM_L2x1_SUB2_2
	LOAD2x1 0 
    KERNEL2x1_L  16,64, 0,0
    KERNEL2x1_E  16,64, 1,1
    MY_ALIGN
ZGEMM_L2x1_SUB2_2:
    andi.      T1,L, 2
    ble ZGEMM_L2x1_SUB2_1
	LOAD2x1 0 
    KERNEL2x1_E  16,64, 0,1
    MY_ALIGN    
ZGEMM_L2x1_SUB2_1:
    andi.      T1,L, 1
    ble ZGEMM_L2x1_SAVE	
    KERNEL2x1 

ZGEMM_L2x1_SAVE:

	SAVE2x1

ZGEMM_L2x1_END:

	slwi		T1,	K,	5
	add		B,	B,	T1

	addic.		J,	J,	-1
	bgt		ZGEMM_L2_BEGIN

	andi.		T2,	N,	1
	ble		L999

ZGEMM_L2_END:

	b		ZGEMM_L1_BEGIN

L999_H1:

	b		L999

ZGEMM_L1_BEGIN:
	andi.		T1,	N,	1
	ble		ZGEMM_L1_END

	mr		BO,	B
	mr		BBO,	BBUFFER 
	srawi.		T1,	K,	3 /*this time K/8 */
	ble		ZGEMM_L1_COPYB1

ZGEMM_L1_COPYB8:

	addi		T2,	PRE, 128
	dcbt		BO,	PRE
	dcbtst		BBO,	PRE
	dcbtst		BBO,	T2
	ZCOPYB_8
	addic.		T1,	T1,	-1

	bgt		ZGEMM_L1_COPYB8

ZGEMM_L1_COPYB1:

	andi.		T1,	K,	7
	ble		ZGEMM_L1_COPYB_END

ZGEMM_L1_COPYB_LOOP:

	ZCOPYB_1
	addic.          T1,     T1,     -1

	bgt             ZGEMM_L1_COPYB_LOOP

ZGEMM_L1_COPYB_END:

	mr		CO,	C
	mr		AO,	A
	srawi.		I,	M,	3
	ble		ZGEMM_L1x8_END

ZGEMM_L1x8_BEGIN:


	mr		BO,	BBUFFER
	mr T1, K
    addi T1,T1, -1
    srawi.		L,	T1,	5 /**(K-1) % 32x */ 
	ZERO1x8  
	ble		ZGEMM_L1x8_SUB0
 

ZGEMM_L1x8_LOOP_START:

    LOAD1x8 0 
    li T2, 1024
	li T3, 1024+512
	li T4, 2048
	li T5, 2048+512
	mtctr		L

	MY_ALIGN
ZGEMM_L1x8_LOOP:
 	dcbt		AO,	PRE
	dcbt		BO,	PRE
    KERNEL1x8_L 128,32,0,0
	KERNEL1x8_L 128,32,1,0
	dcbt		AO,	T2	
	KERNEL1x8_L 128,32,2,0
	KERNEL1x8_L 128,32,3,0 
	dcbt		AO,	T3
	dcbt		BO,	T2
    KERNEL1x8_L 128,32,4,0
	KERNEL1x8_L 128,32,5,0
	dcbt		AO,	T4	
	KERNEL1x8_L 128,32,6,0
	KERNEL1x8_L 128,32,7,0  
	dcbt		AO,	T5	
	dcbt		BO,	T3
    KERNEL1x8_L 128,32,8,0
	KERNEL1x8_L 128,32,9,0
	KERNEL1x8_L 128,32,10,0
	KERNEL1x8_L 128,32,11,0  
	dcbt		BO,	T4
    KERNEL1x8_L 128,32,12,0
	KERNEL1x8_L 128,32,13,0
	KERNEL1x8_L 128,32,14,0
	KERNEL1x8_L 128,32,15,1 		
	bdnz		ZGEMM_L1x8_LOOP
 	MY_ALIGN  
ZGEMM_L1x8_LOOP_END:
    END1x8  AO, BO, 128, 32   	 
 
	b		ZGEMM_L1x8_SUB1
 
ZGEMM_L1x8_SUB0:

	andi.		L,	K,	63
 
	b		ZGEMM_L1x8_SUB2

ZGEMM_L1x8_SUB1:

	andi.		L,	T1,	31
	ble		ZGEMM_L1x8_SAVE

ZGEMM_L1x8_SUB2:
    srawi.      T1,L, 3
    ble ZGEMM_L1x8_SUB2_4
    mtctr		T1
    MY_ALIGN
ZGEMM_L1x8_SUB2_LOOP:
	LOAD1x8 0 
    KERNEL1x8_L  128,32, 0,0
    KERNEL1x8_L  128,32, 1,0
    KERNEL1x8_L  128,32, 2,0
    KERNEL1x8_E  128,32, 3,1
    bdnz ZGEMM_L1x8_SUB2_LOOP 
    MY_ALIGN  
ZGEMM_L1x8_SUB2_4:
    andi.      T1,L, 4
    ble ZGEMM_L1x8_SUB2_2
	LOAD1x8 0 
    KERNEL1x8_L  128,32, 0,0
    KERNEL1x8_E  128,32, 1,1
    MY_ALIGN
ZGEMM_L1x8_SUB2_2:
    andi.      T1,L, 2
    ble ZGEMM_L1x8_SUB2_1
	LOAD1x8 0 
    KERNEL1x8_E  128,32, 0,1
    MY_ALIGN    
ZGEMM_L1x8_SUB2_1:
    andi.      T1,L, 1
    ble ZGEMM_L1x8_SAVE	
    KERNEL1x8      

/*	addic.		L,	L,	-1
	bgt		ZGEMM_L1x8_SUB2_1*/

ZGEMM_L1x8_SAVE:

	SAVE1x8

	addic.		I,	I,	-1
	bgt		ZGEMM_L1x8_BEGIN

ZGEMM_L1x8_END:

ZGEMM_L1x4_BEGIN:

	andi.		T2,	M,	7
	ble		ZGEMM_L1x1_END

	andi.		T1,	M,	4
	ble		ZGEMM_L1x4_END
	mr		BO,	BBUFFER
	mr T1, K
    addi T1,T1, -1
    srawi.		L,	T1, 5 /**(K-1) % 16x */ 
	ZERO1x4  
	ble		ZGEMM_L1x4_SUB0 

ZGEMM_L1x4_LOOP_START:
    LOAD1x4 0  
	mtctr		L

	MY_ALIGN
ZGEMM_L1x4_LOOP: 
    KERNEL1x4_L 64,32,0,0
	KERNEL1x4_L 64,32,1,0 	
	KERNEL1x4_L 64,32,2,0
	KERNEL1x4_L 64,32,3,0  
    KERNEL1x4_L 64,32,4,0
	KERNEL1x4_L 64,32,5,0 
	KERNEL1x4_L 64,32,6,0
	KERNEL1x4_L 64,32,7,0   
    KERNEL1x4_L 64,32,8,0
	KERNEL1x4_L 64,32,9,0
	KERNEL1x4_L 64,32,10,0
	KERNEL1x4_L 64,32,11,0   
    KERNEL1x4_L 64,32,12,0
	KERNEL1x4_L 64,32,13,0
	KERNEL1x4_L 64,32,14,0
	KERNEL1x4_L 64,32,15,1 		
	bdnz		ZGEMM_L1x4_LOOP
 	MY_ALIGN  
ZGEMM_L1x4_LOOP_END:
    END1x4  AO, BO, 64, 32   	 
 
	b		ZGEMM_L1x4_SUB1
 
ZGEMM_L1x4_SUB0:

	andi.		L,	K,	63
 
	b		ZGEMM_L1x4_SUB2

ZGEMM_L1x4_SUB1:

	andi.		L,	T1,	31
	ble		ZGEMM_L1x4_SAVE

ZGEMM_L1x4_SUB2:
    srawi.      T1,L, 3
    ble ZGEMM_L1x4_SUB2_4
    mtctr		T1
    MY_ALIGN
ZGEMM_L1x4_SUB2_LOOP:
	LOAD1x4 0 
    KERNEL1x4_L  64,32, 0,0
    KERNEL1x4_L  64,32, 1,0
    KERNEL1x4_L  64,32, 2,0
    KERNEL1x4_E  64,32, 3,1
    bdnz ZGEMM_L1x4_SUB2_LOOP 
    MY_ALIGN  
ZGEMM_L1x4_SUB2_4:
    andi.      T1,L, 4
    ble ZGEMM_L1x4_SUB2_2
	LOAD1x4 0 
    KERNEL1x4_L  64,32, 0,0
    KERNEL1x4_E  64,32, 1,1
    MY_ALIGN
ZGEMM_L1x4_SUB2_2:
    andi.      T1,L, 2
    ble ZGEMM_L1x4_SUB2_1
	LOAD1x4 0 
    KERNEL1x4_E  64,32, 0,1
    MY_ALIGN    
ZGEMM_L1x4_SUB2_1:
    andi.      T1,L, 1
    ble ZGEMM_L1x4_SAVE	
    KERNEL1x4  

ZGEMM_L1x4_SAVE:

	SAVE1x4

ZGEMM_L1x4_END:

ZGEMM_L1x2_BEGIN:


	andi.		T1,	M,	2
	ble		ZGEMM_L1x2_END
	mr		BO,	BBUFFER
	mr T1, K
    addi T1,T1, -1
    srawi.		L,	T1, 5 /**(K-1) % 16x */ 
	ZERO1x2  
	ble		ZGEMM_L1x2_SUB0 

ZGEMM_L1x2_LOOP_START:
    LOAD1x2 0  
	mtctr		L

	MY_ALIGN
ZGEMM_L1x2_LOOP: 
    KERNEL1x2_L 32,32,0,0
	KERNEL1x2_L 32,32,1,0 	
	KERNEL1x2_L 32,32,2,0
	KERNEL1x2_L 32,32,3,0  
    KERNEL1x2_L 32,32,4,0
	KERNEL1x2_L 32,32,5,0 
	KERNEL1x2_L 32,32,6,0
	KERNEL1x2_L 32,32,7,0   
    KERNEL1x2_L 32,32,8,0
	KERNEL1x2_L 32,32,9,0
	KERNEL1x2_L 32,32,10,0
	KERNEL1x2_L 32,32,11,0   
    KERNEL1x2_L 32,32,12,0
	KERNEL1x2_L 32,32,13,0
	KERNEL1x2_L 32,32,14,0
	KERNEL1x2_L 32,32,15,1 		
	bdnz		ZGEMM_L1x2_LOOP
 	MY_ALIGN  
ZGEMM_L1x2_LOOP_END:
    END1x2  AO, BO, 32, 32   	 
 
	b		ZGEMM_L1x2_SUB1
 
ZGEMM_L1x2_SUB0:

	andi.		L,	K,	63
 
	b		ZGEMM_L1x2_SUB2

ZGEMM_L1x2_SUB1:

	andi.		L,	T1,	31
	ble		ZGEMM_L1x2_SAVE

ZGEMM_L1x2_SUB2:
    srawi.      T1,L, 3
    ble ZGEMM_L1x2_SUB2_4
    mtctr		T1
    MY_ALIGN
ZGEMM_L1x2_SUB2_LOOP:
	LOAD1x2 0 
    KERNEL1x2_L  32,32, 0,0
    KERNEL1x2_L  32,32, 1,0
    KERNEL1x2_L  32,32, 2,0
    KERNEL1x2_E  32,32, 3,1
    bdnz ZGEMM_L1x2_SUB2_LOOP 
    MY_ALIGN  
ZGEMM_L1x2_SUB2_4:
    andi.      T1,L, 4
    ble ZGEMM_L1x2_SUB2_2
	LOAD1x2 0 
    KERNEL1x2_L  32,32, 0,0
    KERNEL1x2_E  32,32, 1,1
    MY_ALIGN
ZGEMM_L1x2_SUB2_2:
    andi.      T1,L, 2
    ble ZGEMM_L1x2_SUB2_1
	LOAD1x2 0 
    KERNEL1x2_E  32,32, 0,1
    MY_ALIGN    
ZGEMM_L1x2_SUB2_1:
    andi.      T1,L, 1
    ble ZGEMM_L1x2_SAVE	
    KERNEL1x2 
ZGEMM_L1x2_SAVE:

	SAVE1x2

ZGEMM_L1x2_END:

ZGEMM_L1x1_BEGIN:


	andi.		T1,	M,	1
	ble		ZGEMM_L1x1_END
	mr		BO,	BBUFFER
	mr T1, K
    addi T1,T1, -1
    srawi.		L,	T1, 5 /**(K-1) % 16x */ 
	ZERO1x1  
	ble		ZGEMM_L1x1_SUB0 

ZGEMM_L1x1_LOOP_START:

    LOAD1x1 0  
	mtctr		L

	MY_ALIGN
ZGEMM_L1x1_LOOP: 
    KERNEL1x1_L 16,32,0,0
	KERNEL1x1_L 16,32,1,0 	
	KERNEL1x1_L 16,32,2,0
	KERNEL1x1_L 16,32,3,0  
    KERNEL1x1_L 16,32,4,0
	KERNEL1x1_L 16,32,5,0 
	KERNEL1x1_L 16,32,6,0
	KERNEL1x1_L 16,32,7,0   
    KERNEL1x1_L 16,32,8,0
	KERNEL1x1_L 16,32,9,0
	KERNEL1x1_L 16,32,10,0
	KERNEL1x1_L 16,32,11,0   
    KERNEL1x1_L 16,32,12,0
	KERNEL1x1_L 16,32,13,0
	KERNEL1x1_L 16,32,14,0
	KERNEL1x1_L 16,32,15,1 		
	bdnz		ZGEMM_L1x1_LOOP
 	MY_ALIGN  
ZGEMM_L1x1_LOOP_END:
    END1x1  AO, BO, 16, 32   	 
 
	b		ZGEMM_L1x1_SUB1
 
ZGEMM_L1x1_SUB0:

	andi.		L,	K,	63
 
	b		ZGEMM_L1x1_SUB2

ZGEMM_L1x1_SUB1:

	andi.		L,	T1,	31
	ble		ZGEMM_L1x1_SAVE

ZGEMM_L1x1_SUB2:
    srawi.      T1,L, 3
    ble ZGEMM_L1x1_SUB2_4
    mtctr		T1
    MY_ALIGN
ZGEMM_L1x1_SUB2_LOOP:
	LOAD1x1 0 
    KERNEL1x1_L  16,32, 0,0
    KERNEL1x1_L  16,32, 1,0
    KERNEL1x1_L  16,32, 2,0
    KERNEL1x1_E  16,32, 3,1
    bdnz ZGEMM_L1x1_SUB2_LOOP 
    MY_ALIGN  
ZGEMM_L1x1_SUB2_4:
    andi.      T1,L, 4
    ble ZGEMM_L1x1_SUB2_2
	LOAD1x1 0 
    KERNEL1x1_L  16,32, 0,0
    KERNEL1x1_E  16,32, 1,1
    MY_ALIGN
ZGEMM_L1x1_SUB2_2:
    andi.      T1,L, 2
    ble ZGEMM_L1x1_SUB2_1
	LOAD1x1 0 
    KERNEL1x1_E  16,32, 0,1
    MY_ALIGN    
ZGEMM_L1x1_SUB2_1:
    andi.      T1,L, 1
    ble ZGEMM_L1x1_SAVE	
    KERNEL1x1 

ZGEMM_L1x1_SAVE:

	SAVE1x1

ZGEMM_L1x1_END:

ZGEMM_L1_END:
