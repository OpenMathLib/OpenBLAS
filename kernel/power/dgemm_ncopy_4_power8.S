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
* 2016/04/28 Werner Saar (wernsaar@googlemail.com)
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
#define o64	r17
#define o80	r18
#define o96	r19
#define o112	r20
#define o8	r21
#define T2	r22
#define I	r23
#define o16	r24
#define o32	r25
#define o48	r26
#define NOTU1	r27
#define NOTU2	r30
#define T1	r31

#define o0	0

#include "dgemm_ncopy_macros_4_power8.S"

#define STACKSIZE 384


	PROLOGUE
	PROFCODE

	addi	SP, SP, -STACKSIZE
	li	r0, 0

        stfd    f14,    0(SP)
        stfd    f15,    8(SP)
        stfd    f16,   16(SP)
        stfd    f17,   24(SP)
        stfd    f18,   32(SP)
        stfd    f19,   40(SP)
        stfd    f20,   48(SP)
        stfd    f21,   56(SP)
        stfd    f22,   64(SP)
        stfd    f23,   72(SP)
        stfd    f24,   80(SP)
        stfd    f25,   88(SP)
        stfd    f26,   96(SP)
        stfd    f27,  104(SP)
        stfd    f28,  112(SP)
        stfd    f29,  120(SP)
        stfd    f30,  128(SP)
        stfd    f31,  136(SP)


        std     r31,  144(SP)
        std     r30,  152(SP)
        std     r29,  160(SP)
        std     r28,  168(SP)
        std     r27,  176(SP)
        std     r26,  184(SP)
        std     r25,  192(SP)
        std     r24,  200(SP)
        std     r23,  208(SP)
        std     r22,  216(SP)
        std     r21,  224(SP)
        std     r20,  232(SP)
        std     r19,  240(SP)
        std     r18,  248(SP)
        std     r17,  256(SP)
        std     r16,  264(SP)
        std     r15,  272(SP)
        std     r14,  280(SP)

	cmpwi	cr0, M, 0
	ble-	L999
	cmpwi	cr0, N, 0
	ble-	L999

	slwi	LDA, LDA, BASE_SHIFT

	li	PREA,  384
	li	PREB,  384

	li	o8,	8
	li	o16,	16
	li	o32,	32
	li	o48,	48
	li	o64,	64
	li	o80,	80
	li	o96,	96
	li	o112,	112

#include "dgemm_ncopy_logic_4_power8.S"

L999:

	li	r3, 0

        lfd     f14,    0(SP)
        lfd     f15,    8(SP)
        lfd     f16,   16(SP)
        lfd     f17,   24(SP)
        lfd     f18,   32(SP)
        lfd     f19,   40(SP)
        lfd     f20,   48(SP)
        lfd     f21,   56(SP)
        lfd     f22,   64(SP)
        lfd     f23,   72(SP)
        lfd     f24,   80(SP)
        lfd     f25,   88(SP)
        lfd     f26,   96(SP)
        lfd     f27,  104(SP)
        lfd     f28,  112(SP)
        lfd     f29,  120(SP)
        lfd     f30,  128(SP)
        lfd     f31,  136(SP)

        ld      r31,  144(SP)
        ld      r30,  152(SP)
        ld      r29,  160(SP)
        ld      r28,  168(SP)
        ld      r27,  176(SP)
        ld      r26,  184(SP)
        ld      r25,  192(SP)
        ld      r24,  200(SP)
        ld      r23,  208(SP)
        ld      r22,  216(SP)
        ld      r21,  224(SP)
        ld      r20,  232(SP)
        ld      r19,  240(SP)
        ld      r18,  248(SP)
        ld      r17,  256(SP)
        ld      r16,  264(SP)
        ld      r15,  272(SP)
        ld      r14,  280(SP)

	addi	SP, SP, STACKSIZE

	blr
	EPILOGUE


