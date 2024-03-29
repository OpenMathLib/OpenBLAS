/***************************************************************************
Copyright (c) 2013-2016, The OpenBLAS Project
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
* 2016/04/21 Werner Saar (wernsaar@googlemail.com)
* 	 BLASTEST 		: OK
* 	 CTEST			: OK
* 	 TEST			: OK
*	 LAPACK-TEST		: OK
**************************************************************************************/

/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/*   1. Redistributions of source code must retain the above         */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer.                                                  */
/*                                                                   */
/*   2. Redistributions in binary form must reproduce the above      */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer in the documentation and/or other materials       */
/*      provided with the distribution.                              */
/*                                                                   */
/*    THIS  SOFTWARE IS PROVIDED  BY THE  UNIVERSITY OF  TEXAS AT    */
/*    AUSTIN  ``AS IS''  AND ANY  EXPRESS OR  IMPLIED WARRANTIES,    */
/*    INCLUDING, BUT  NOT LIMITED  TO, THE IMPLIED  WARRANTIES OF    */
/*    MERCHANTABILITY  AND FITNESS FOR  A PARTICULAR  PURPOSE ARE    */
/*    DISCLAIMED.  IN  NO EVENT SHALL THE UNIVERSITY  OF TEXAS AT    */
/*    AUSTIN OR CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT,    */
/*    INCIDENTAL,  SPECIAL, EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES    */
/*    (INCLUDING, BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE    */
/*    GOODS  OR  SERVICES; LOSS  OF  USE,  DATA,  OR PROFITS;  OR    */
/*    BUSINESS INTERRUPTION) HOWEVER CAUSED  AND ON ANY THEORY OF    */
/*    LIABILITY, WHETHER  IN CONTRACT, STRICT  LIABILITY, OR TORT    */
/*    (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY WAY OUT    */
/*    OF  THE  USE OF  THIS  SOFTWARE,  EVEN  IF ADVISED  OF  THE    */
/*    POSSIBILITY OF SUCH DAMAGE.                                    */
/*                                                                   */
/* The views and conclusions contained in the software and           */
/* documentation are those of the authors and should not be          */
/* interpreted as representing official policies, either expressed   */
/* or implied, of The University of Texas at Austin.                 */
/*********************************************************************/

#define ASSEMBLER
#include "common.h"
#include "def_vsx.h"

#define	M	r3
#define	N	r4
#define	A	r5
#define	LDA	r6
#define B	r7

#define A0	r8
#define A1	r9
#define A2	r10
#define A3	r11

#define J	r12

#define PREA	r14
#define PREB 	r15
#define BO	r16
#define B8	r17
#define B4	r18
#define B2	r19
#define B1	r20
#define o8	r21
#define T2	r22
#define I	r23
#define o16	r24
#define o32	r25
#define o48	r26
#define B16	r29
#define M16	r30
#define T1	r31

#define o0	0

#ifdef POWER10
#include "dgemm_tcopy_macros_16_power10.S"
#endif
#include "dgemm_tcopy_macros_16_power8.S"

#define STACKSIZE 144
 


	PROLOGUE
	PROFCODE

	addi	SP, SP, -STACKSIZE

	li	r0, 0

    std r14,0(SP) 
    std r15,8(SP) 
    std r16,16(SP) 
    std r17,24(SP) 
    std r18,32(SP) 
    std r19,40(SP) 
    std r20,48(SP) 
    std r21,56(SP) 
    std r22,64(SP) 
    std r23,72(SP) 
    std r24,80(SP) 
    std r25,88(SP) 
    std r26,96(SP) 
    std r27,104(SP) 
    std r28,112(SP) 
    std r29,120(SP) 
    std r30,128(SP) 
    std r31,136(SP) 

	cmpwi	cr0, M, 0
	ble-	L999
	cmpwi	cr0, N, 0
	ble-	L999

	slwi	LDA, LDA, BASE_SHIFT
	slwi	M16, M, 4 + BASE_SHIFT

	li	T1,	-16
	li	T2,	-8
	li	PREA,	-4
	li	PREB,	-2

	and	B8, N, T1
	and	B4, N, T2
	and	B2, N, PREA
	and	B1, N, PREB
	
	mullw	B8, B8, M
	mullw	B4, B4, M
	mullw	B2, B2, M
	mullw	B1, B1, M

	slwi	B8, B8, BASE_SHIFT
	slwi	B4, B4, BASE_SHIFT
	slwi	B2, B2, BASE_SHIFT
	slwi	B1, B1, BASE_SHIFT

	add	B8, B8, B
	add	B4, B4, B
	add	B2, B2, B
	add	B1, B1, B

	li	PREA,  384 
	addi	PREB,  M16, 128

	li	o8,	8
	li	o16,	16
	li	o32,	32
	li	o48,	48

#include "dgemm_tcopy_logic_16_power8.S"

L999:

	li	r3, 0

    ld r14,0(SP) 
    ld r15,8(SP) 
    ld r16,16(SP) 
    ld r17,24(SP) 
    ld r18,32(SP) 
    ld r19,40(SP) 
    ld r20,48(SP) 
    ld r21,56(SP) 
    ld r22,64(SP) 
    ld r23,72(SP) 
    ld r24,80(SP) 
    ld r25,88(SP) 
    ld r26,96(SP) 
    ld r27,104(SP) 
    ld r28,112(SP) 
    ld r29,120(SP) 
    ld r30,128(SP) 
    ld r31,136(SP) 

	addi	SP, SP, STACKSIZE
 
	blr
	EPILOGUE


