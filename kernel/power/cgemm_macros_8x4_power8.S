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
* 2016/03/18 Werner Saar (wernsaar@googlemail.com)
* 	 BLASTEST 		: OK
* 	 CTEST			: OK
* 	 TEST			: OK
*	 LAPACK-TEST		: OK
**************************************************************************************/

#if   defined(NN) || defined(NT) || defined(TN) || defined(TT)

	#define	XSFADD_R1   xsaddsp
	#define	XSFADD_R2   xssubsp
	#define	XSFADD_I1   xsaddsp
	#define	XSFADD_I2   xsaddsp

#elif  defined(CN) || defined(CT) || defined(RN) || defined(RT)

	#define	XSFADD_R1   xsaddsp
	#define	XSFADD_R2   xsaddsp
	#define	XSFADD_I1   xssubsp
	#define	XSFADD_I2   xsaddsp

#elif  defined(NC) || defined(TC) || defined(NR) || defined(TR)

	#define	XSFADD_R1   xsaddsp
	#define	XSFADD_R2   xsaddsp
	#define	XSFADD_I1   xsaddsp
	#define	XSFADD_I2   xssubsp

#else             // CC || CR || RC || RR

	#define	XSFADD_R1   xsaddsp
	#define	XSFADD_R2   xssubsp
	#define	XSFADD_I1   xssubsp
	#define	XSFADD_I2   xssubsp

#endif

/**********************************************************************************************
* Macros for N=4 and M=8
**********************************************************************************************/

.macro LOAD4x8_1

	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3

	lxvw4x		vs2,	o32,	AO		// load a4, a5

	lxvw4x		vs3,	o48,	AO		// load a6, a7


	addi		AO,	AO,	64

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3

	lxvw4x		vs25,	o16,	BO		//  load b2, b3

	xxspltw		vs12,	vs25,	0
	xxspltw		vs13,	vs25,	1
	xxspltw		vs14,	vs25,	2
	xxspltw		vs15,	vs25,	3


	addi		BO,	BO,	32

.endm

.macro KERNEL4x8_I1


	lxvw4x		vs4,	o0,	AO		// load a0, a1

	lxvw4x		vs5,	o16,	AO		// load a2, a3

	lxvw4x		vs6,	o32,	AO		// load a4, a5

	lxvw4x		vs7,	o48,	AO		// load a6, a7


	addi		AO,	AO,	64

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs16,	vs24,	0
	xxspltw		vs17,	vs24,	1
	xxspltw		vs18,	vs24,	2
	xxspltw		vs19,	vs24,	3

	lxvw4x		vs25,	o16,	BO		//  load b2, b3

	xxspltw		vs20,	vs25,	0
	xxspltw		vs21,	vs25,	1
	xxspltw		vs22,	vs25,	2
	xxspltw		vs23,	vs25,	3


	addi		BO,	BO,	32


	xvmulsp		vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs36,	vs2,	vs8		// a2_r*b0_r, a2_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs37,	vs2,	vs9		// a2_r*b0_i, a2_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs38,	vs3,	vs8		// a3_r*b0_r, a3_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs39,	vs3,	vs9		// a3_r*b0_i, a3_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmulsp		vs40,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs41,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmulsp		vs42,	vs1,	vs10		// a1_r*b1_r, a1_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs43,	vs1,	vs11		// a1_r*b1_i, a1_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmulsp		vs44,	vs2,	vs10		// a2_r*b1_r, a2_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs45,	vs2,	vs11		// a2_r*b1_i, a2_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmulsp		vs46,	vs3,	vs10		// a3_r*b1_r, a3_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs47,	vs3,	vs11		// a3_r*b1_i, a3_i*b1_i, a1_r*b1_i, a1_i*b1_i

	xvmulsp		vs48,	vs0,	vs12		// a0_r*b2_r, a0_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmulsp		vs49,	vs0,	vs13		// a0_r*b2_i, a0_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmulsp		vs50,	vs1,	vs12		// a1_r*b2_r, a1_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmulsp		vs51,	vs1,	vs13		// a1_r*b2_i, a1_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmulsp		vs52,	vs2,	vs12		// a2_r*b2_r, a2_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmulsp		vs53,	vs2,	vs13		// a2_r*b2_i, a2_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmulsp		vs54,	vs3,	vs12		// a3_r*b2_r, a3_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmulsp		vs55,	vs3,	vs13		// a3_r*b2_i, a3_i*b2_i, a1_r*b2_i, a1_i*b2_i

	xvmulsp		vs56,	vs0,	vs14		// a0_r*b3_r, a0_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmulsp		vs57,	vs0,	vs15		// a0_r*b3_i, a0_i*b3_i, a1_r*b3_i, a1_i*b3_i
	xvmulsp		vs58,	vs1,	vs14		// a1_r*b3_r, a1_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmulsp		vs59,	vs1,	vs15		// a1_r*b3_i, a1_i*b3_i, a1_r*b3_i, a1_i*b3_i
	xvmulsp		vs60,	vs2,	vs14		// a2_r*b3_r, a2_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmulsp		vs61,	vs2,	vs15		// a2_r*b3_i, a2_i*b3_i, a1_r*b3_i, a1_i*b3_i
	xvmulsp		vs62,	vs3,	vs14		// a3_r*b3_r, a3_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmulsp		vs63,	vs3,	vs15		// a3_r*b3_i, a3_i*b3_i, a1_r*b3_i, a1_i*b3_i


.endm

.macro KERNEL4x8_1

	xvmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i

	lxvw4x		vs24,	o0,	BO		//  load b0, b1
	lxvw4x		vs4,	o0,	AO		// load a0, a1

	xvmaddasp	vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i

	lxvw4x		vs25,	o16,	BO		//  load b2, b3
	lxvw4x		vs5,	o16,	AO		// load a2, a3

	xvmaddasp	vs36,	vs2,	vs8		// a2_r*b0_r, a2_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs37,	vs2,	vs9		// a2_r*b0_i, a2_i*b0_i, a1_r*b0_i, a1_i*b0_i

	lxvw4x		vs6,	o32,	AO		// load a4, a5
	lxvw4x		vs7,	o48,	AO		// load a6, a7

	xvmaddasp	vs38,	vs3,	vs8		// a3_r*b0_r, a3_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs39,	vs3,	vs9		// a3_r*b0_i, a3_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs40,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs41,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs42,	vs1,	vs10		// a1_r*b1_r, a1_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs43,	vs1,	vs11		// a1_r*b1_i, a1_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs44,	vs2,	vs10		// a2_r*b1_r, a2_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs45,	vs2,	vs11		// a2_r*b1_i, a2_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs46,	vs3,	vs10		// a3_r*b1_r, a3_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs47,	vs3,	vs11		// a3_r*b1_i, a3_i*b1_i, a1_r*b1_i, a1_i*b1_i

	xvmaddasp	vs48,	vs0,	vs12		// a0_r*b2_r, a0_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs49,	vs0,	vs13		// a0_r*b2_i, a0_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmaddasp	vs50,	vs1,	vs12		// a1_r*b2_r, a1_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs51,	vs1,	vs13		// a1_r*b2_i, a1_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmaddasp	vs52,	vs2,	vs12		// a2_r*b2_r, a2_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs53,	vs2,	vs13		// a2_r*b2_i, a2_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmaddasp	vs54,	vs3,	vs12		// a3_r*b2_r, a3_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs55,	vs3,	vs13		// a3_r*b2_i, a3_i*b2_i, a1_r*b2_i, a1_i*b2_i

	xxspltw		vs16,	vs24,	0
	xxspltw		vs17,	vs24,	1
	xxspltw		vs18,	vs24,	2
	xxspltw		vs19,	vs24,	3

	xvmaddasp	vs56,	vs0,	vs14		// a0_r*b3_r, a0_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs57,	vs0,	vs15		// a0_r*b3_i, a0_i*b3_i, a1_r*b3_i, a1_i*b3_i
	xvmaddasp	vs58,	vs1,	vs14		// a1_r*b3_r, a1_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs59,	vs1,	vs15		// a1_r*b3_i, a1_i*b3_i, a1_r*b3_i, a1_i*b3_i
	addi		BO,	BO,	32
	xvmaddasp	vs60,	vs2,	vs14		// a2_r*b3_r, a2_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs61,	vs2,	vs15		// a2_r*b3_i, a2_i*b3_i, a1_r*b3_i, a1_i*b3_i
	addi		AO,	AO,	64
	xvmaddasp	vs62,	vs3,	vs14		// a3_r*b3_r, a3_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs63,	vs3,	vs15		// a3_r*b3_i, a3_i*b3_i, a1_r*b3_i, a1_i*b3_i

	xxspltw		vs20,	vs25,	0
	xxspltw		vs21,	vs25,	1
	xxspltw		vs22,	vs25,	2
	xxspltw		vs23,	vs25,	3

.endm

.macro KERNEL4x8_2

	xvmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r, a4_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs4,	vs17		// a4_r*b0_i, a4_i*b0_i, a1_r*b0_i, a1_i*b0_i

	lxvw4x		vs24,	o0,	BO		//  load b0, b1
	lxvw4x		vs0,	o0,	AO		// load a0, a1

	xvmaddasp	vs34,	vs5,	vs16		// a5_r*b0_r, a5_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs5,	vs17		// a5_r*b0_i, a5_i*b0_i, a1_r*b0_i, a1_i*b0_i

	lxvw4x		vs25,	o16,	BO		//  load b2, b3
	lxvw4x		vs1,	o16,	AO		// load a2, a3

	xvmaddasp	vs36,	vs6,	vs16		// a6_r*b0_r, a6_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs37,	vs6,	vs17		// a6_r*b0_i, a6_i*b0_i, a1_r*b0_i, a1_i*b0_i
	lxvw4x		vs2,	o32,	AO		// load a4, a5
	lxvw4x		vs3,	o48,	AO		// load a6, a7
	xvmaddasp	vs38,	vs7,	vs16		// a7_r*b0_r, a7_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs39,	vs7,	vs17		// a7_r*b0_i, a7_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs40,	vs4,	vs18		// a4_r*b1_r, a4_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs41,	vs4,	vs19		// a4_r*b1_i, a4_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs42,	vs5,	vs18		// a5_r*b1_r, a5_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs43,	vs5,	vs19		// a5_r*b1_i, a5_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs44,	vs6,	vs18		// a6_r*b1_r, a6_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs45,	vs6,	vs19		// a6_r*b1_i, a6_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs46,	vs7,	vs18		// a7_r*b1_r, a7_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs47,	vs7,	vs19		// a7_r*b1_i, a7_i*b1_i, a1_r*b1_i, a1_i*b1_i

	xvmaddasp	vs48,	vs4,	vs20		// a4_r*b2_r, a4_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs49,	vs4,	vs21		// a4_r*b2_i, a4_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmaddasp	vs50,	vs5,	vs20		// a5_r*b2_r, a5_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs51,	vs5,	vs21		// a5_r*b2_i, a5_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmaddasp	vs52,	vs6,	vs20		// a6_r*b2_r, a6_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs53,	vs6,	vs21		// a6_r*b2_i, a6_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmaddasp	vs54,	vs7,	vs20		// a7_r*b2_r, a7_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs55,	vs7,	vs21		// a7_r*b2_i, a7_i*b2_i, a1_r*b2_i, a1_i*b2_i

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3

	xvmaddasp	vs56,	vs4,	vs22		// a4_r*b3_r, a4_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs57,	vs4,	vs23		// a4_r*b3_i, a4_i*b3_i, a1_r*b3_i, a1_i*b3_i
	xvmaddasp	vs58,	vs5,	vs22		// a5_r*b3_r, a5_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs59,	vs5,	vs23		// a5_r*b3_i, a5_i*b3_i, a1_r*b3_i, a1_i*b3_i
	addi		AO,	AO,	64
	xvmaddasp	vs60,	vs6,	vs22		// a6_r*b3_r, a6_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs61,	vs6,	vs23		// a6_r*b3_i, a6_i*b3_i, a1_r*b3_i, a1_i*b3_i
	addi		BO,	BO,	32
	xvmaddasp	vs62,	vs7,	vs22		// a7_r*b3_r, a7_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs63,	vs7,	vs23		// a7_r*b3_i, a7_i*b3_i, a1_r*b3_i, a1_i*b3_i

	xxspltw		vs12,	vs25,	0
	xxspltw		vs13,	vs25,	1
	xxspltw		vs14,	vs25,	2
	xxspltw		vs15,	vs25,	3

.endm

.macro KERNEL4x8_E2


	xvmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r, a4_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs4,	vs17		// a4_r*b0_i, a4_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs5,	vs16		// a5_r*b0_r, a5_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs5,	vs17		// a5_r*b0_i, a5_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs36,	vs6,	vs16		// a6_r*b0_r, a6_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs37,	vs6,	vs17		// a6_r*b0_i, a6_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs38,	vs7,	vs16		// a7_r*b0_r, a7_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs39,	vs7,	vs17		// a7_r*b0_i, a7_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs40,	vs4,	vs18		// a4_r*b1_r, a4_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs41,	vs4,	vs19		// a4_r*b1_i, a4_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs42,	vs5,	vs18		// a5_r*b1_r, a5_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs43,	vs5,	vs19		// a5_r*b1_i, a5_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs44,	vs6,	vs18		// a6_r*b1_r, a6_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs45,	vs6,	vs19		// a6_r*b1_i, a6_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs46,	vs7,	vs18		// a7_r*b1_r, a7_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs47,	vs7,	vs19		// a7_r*b1_i, a7_i*b1_i, a1_r*b1_i, a1_i*b1_i

	xvmaddasp	vs48,	vs4,	vs20		// a4_r*b2_r, a4_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs49,	vs4,	vs21		// a4_r*b2_i, a4_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmaddasp	vs50,	vs5,	vs20		// a5_r*b2_r, a5_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs51,	vs5,	vs21		// a5_r*b2_i, a5_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmaddasp	vs52,	vs6,	vs20		// a6_r*b2_r, a6_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs53,	vs6,	vs21		// a6_r*b2_i, a6_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmaddasp	vs54,	vs7,	vs20		// a7_r*b2_r, a7_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs55,	vs7,	vs21		// a7_r*b2_i, a7_i*b2_i, a1_r*b2_i, a1_i*b2_i

	xvmaddasp	vs56,	vs4,	vs22		// a4_r*b3_r, a4_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs57,	vs4,	vs23		// a4_r*b3_i, a4_i*b3_i, a1_r*b3_i, a1_i*b3_i
	xvmaddasp	vs58,	vs5,	vs22		// a5_r*b3_r, a5_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs59,	vs5,	vs23		// a5_r*b3_i, a5_i*b3_i, a1_r*b3_i, a1_i*b3_i
	xvmaddasp	vs60,	vs6,	vs22		// a6_r*b3_r, a6_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs61,	vs6,	vs23		// a6_r*b3_i, a6_i*b3_i, a1_r*b3_i, a1_i*b3_i
	xvmaddasp	vs62,	vs7,	vs22		// a7_r*b3_r, a7_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs63,	vs7,	vs23		// a7_r*b3_i, a7_i*b3_i, a1_r*b3_i, a1_i*b3_i


.endm

.macro KERNEL4x8_SUBI1


	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3

	lxvw4x		vs2,	o32,	AO		// load a4, a5

	lxvw4x		vs3,	o48,	AO		// load a6, a7


	addi		AO,	AO,	64

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3

	lxvw4x		vs25,	o16,	BO		//  load b2, b3

	xxspltw		vs12,	vs25,	0
	xxspltw		vs13,	vs25,	1
	xxspltw		vs14,	vs25,	2
	xxspltw		vs15,	vs25,	3


	addi		BO,	BO,	32


	xvmulsp		vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs36,	vs2,	vs8		// a2_r*b0_r, a2_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs37,	vs2,	vs9		// a2_r*b0_i, a2_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs38,	vs3,	vs8		// a3_r*b0_r, a3_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs39,	vs3,	vs9		// a3_r*b0_i, a3_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmulsp		vs40,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs41,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmulsp		vs42,	vs1,	vs10		// a1_r*b1_r, a1_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs43,	vs1,	vs11		// a1_r*b1_i, a1_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmulsp		vs44,	vs2,	vs10		// a2_r*b1_r, a2_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs45,	vs2,	vs11		// a2_r*b1_i, a2_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmulsp		vs46,	vs3,	vs10		// a3_r*b1_r, a3_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs47,	vs3,	vs11		// a3_r*b1_i, a3_i*b1_i, a1_r*b1_i, a1_i*b1_i

	xvmulsp		vs48,	vs0,	vs12		// a0_r*b2_r, a0_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmulsp		vs49,	vs0,	vs13		// a0_r*b2_i, a0_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmulsp		vs50,	vs1,	vs12		// a1_r*b2_r, a1_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmulsp		vs51,	vs1,	vs13		// a1_r*b2_i, a1_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmulsp		vs52,	vs2,	vs12		// a2_r*b2_r, a2_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmulsp		vs53,	vs2,	vs13		// a2_r*b2_i, a2_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmulsp		vs54,	vs3,	vs12		// a3_r*b2_r, a3_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmulsp		vs55,	vs3,	vs13		// a3_r*b2_i, a3_i*b2_i, a1_r*b2_i, a1_i*b2_i

	xvmulsp		vs56,	vs0,	vs14		// a0_r*b3_r, a0_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmulsp		vs57,	vs0,	vs15		// a0_r*b3_i, a0_i*b3_i, a1_r*b3_i, a1_i*b3_i
	xvmulsp		vs58,	vs1,	vs14		// a1_r*b3_r, a1_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmulsp		vs59,	vs1,	vs15		// a1_r*b3_i, a1_i*b3_i, a1_r*b3_i, a1_i*b3_i
	xvmulsp		vs60,	vs2,	vs14		// a2_r*b3_r, a2_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmulsp		vs61,	vs2,	vs15		// a2_r*b3_i, a2_i*b3_i, a1_r*b3_i, a1_i*b3_i
	xvmulsp		vs62,	vs3,	vs14		// a3_r*b3_r, a3_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmulsp		vs63,	vs3,	vs15		// a3_r*b3_i, a3_i*b3_i, a1_r*b3_i, a1_i*b3_i


.endm

.macro KERNEL4x8_SUB1


	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3

	lxvw4x		vs2,	o32,	AO		// load a4, a5

	lxvw4x		vs3,	o48,	AO		// load a6, a7


	addi		AO,	AO,	64

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3

	lxvw4x		vs25,	o16,	BO		//  load b2, b3

	xxspltw		vs12,	vs25,	0
	xxspltw		vs13,	vs25,	1
	xxspltw		vs14,	vs25,	2
	xxspltw		vs15,	vs25,	3


	addi		BO,	BO,	32


	xvmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs36,	vs2,	vs8		// a2_r*b0_r, a2_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs37,	vs2,	vs9		// a2_r*b0_i, a2_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs38,	vs3,	vs8		// a3_r*b0_r, a3_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs39,	vs3,	vs9		// a3_r*b0_i, a3_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs40,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs41,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs42,	vs1,	vs10		// a1_r*b1_r, a1_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs43,	vs1,	vs11		// a1_r*b1_i, a1_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs44,	vs2,	vs10		// a2_r*b1_r, a2_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs45,	vs2,	vs11		// a2_r*b1_i, a2_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs46,	vs3,	vs10		// a3_r*b1_r, a3_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs47,	vs3,	vs11		// a3_r*b1_i, a3_i*b1_i, a1_r*b1_i, a1_i*b1_i

	xvmaddasp	vs48,	vs0,	vs12		// a0_r*b2_r, a0_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs49,	vs0,	vs13		// a0_r*b2_i, a0_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmaddasp	vs50,	vs1,	vs12		// a1_r*b2_r, a1_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs51,	vs1,	vs13		// a1_r*b2_i, a1_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmaddasp	vs52,	vs2,	vs12		// a2_r*b2_r, a2_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs53,	vs2,	vs13		// a2_r*b2_i, a2_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmaddasp	vs54,	vs3,	vs12		// a3_r*b2_r, a3_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs55,	vs3,	vs13		// a3_r*b2_i, a3_i*b2_i, a1_r*b2_i, a1_i*b2_i

	xvmaddasp	vs56,	vs0,	vs14		// a0_r*b3_r, a0_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs57,	vs0,	vs15		// a0_r*b3_i, a0_i*b3_i, a1_r*b3_i, a1_i*b3_i
	xvmaddasp	vs58,	vs1,	vs14		// a1_r*b3_r, a1_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs59,	vs1,	vs15		// a1_r*b3_i, a1_i*b3_i, a1_r*b3_i, a1_i*b3_i
	xvmaddasp	vs60,	vs2,	vs14		// a2_r*b3_r, a2_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs61,	vs2,	vs15		// a2_r*b3_i, a2_i*b3_i, a1_r*b3_i, a1_i*b3_i
	xvmaddasp	vs62,	vs3,	vs14		// a3_r*b3_r, a3_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs63,	vs3,	vs15		// a3_r*b3_i, a3_i*b3_i, a1_r*b3_i, a1_i*b3_i


.endm

.macro SAVE4x8

	mr		T1,	CO

// N=0

	mr		T2,	T1

// N=0 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs33,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=0 M=2

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs34,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs35,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=0 M=4

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs36,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs37,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=0 M=6

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs38,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs39,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC


// N=1

	mr		T2,	T1

// N=1 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs40,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs41,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=1 M=2

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs42,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs43,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=1 M=4

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs44,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs45,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=1 M=6

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs46,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs47,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC


// N=2

	mr		T2,	T1

// N=2 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs48,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs49,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=2 M=2

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs50,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs51,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=2 M=4

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs52,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs53,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=2 M=6

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs54,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs55,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC


// N=3

	mr		T2,	T1

// N=3 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs56,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs57,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=3 M=2

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs58,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs59,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=3 M=4

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs60,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs61,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=3 M=6

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs62,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs63,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC

	addi		CO,	CO,	64

.endm


/**********************************************************************************************
* Macros for N=4 and M=4
**********************************************************************************************/

.macro LOAD4x4_1

	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3


	addi		AO,	AO,	32

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3

	lxvw4x		vs25,	o16,	BO		//  load b2, b3

	xxspltw		vs12,	vs25,	0
	xxspltw		vs13,	vs25,	1
	xxspltw		vs14,	vs25,	2
	xxspltw		vs15,	vs25,	3


	addi		BO,	BO,	32

.endm

.macro KERNEL4x4_I1


	lxvw4x		vs4,	o0,	AO		// load a0, a1

	lxvw4x		vs5,	o16,	AO		// load a2, a3


	addi		AO,	AO,	32

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs16,	vs24,	0
	xxspltw		vs17,	vs24,	1
	xxspltw		vs18,	vs24,	2
	xxspltw		vs19,	vs24,	3

	lxvw4x		vs25,	o16,	BO		//  load b2, b3

	xxspltw		vs20,	vs25,	0
	xxspltw		vs21,	vs25,	1
	xxspltw		vs22,	vs25,	2
	xxspltw		vs23,	vs25,	3


	addi		BO,	BO,	32


	xvmulsp		vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmulsp		vs36,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs37,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmulsp		vs38,	vs1,	vs10		// a1_r*b1_r, a1_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs39,	vs1,	vs11		// a1_r*b1_i, a1_i*b1_i, a1_r*b1_i, a1_i*b1_i

	xvmulsp		vs40,	vs0,	vs12		// a0_r*b2_r, a0_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmulsp		vs41,	vs0,	vs13		// a0_r*b2_i, a0_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmulsp		vs42,	vs1,	vs12		// a1_r*b2_r, a1_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmulsp		vs43,	vs1,	vs13		// a1_r*b2_i, a1_i*b2_i, a1_r*b2_i, a1_i*b2_i

	xvmulsp		vs44,	vs0,	vs14		// a0_r*b3_r, a0_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmulsp		vs45,	vs0,	vs15		// a0_r*b3_i, a0_i*b3_i, a1_r*b3_i, a1_i*b3_i
	xvmulsp		vs46,	vs1,	vs14		// a1_r*b3_r, a1_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmulsp		vs47,	vs1,	vs15		// a1_r*b3_i, a1_i*b3_i, a1_r*b3_i, a1_i*b3_i


.endm

.macro KERNEL4x4_1


	lxvw4x		vs4,	o0,	AO		// load a0, a1

	lxvw4x		vs5,	o16,	AO		// load a2, a3


	addi		AO,	AO,	32

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs16,	vs24,	0
	xxspltw		vs17,	vs24,	1
	xxspltw		vs18,	vs24,	2
	xxspltw		vs19,	vs24,	3

	lxvw4x		vs25,	o16,	BO		//  load b2, b3

	xxspltw		vs20,	vs25,	0
	xxspltw		vs21,	vs25,	1
	xxspltw		vs22,	vs25,	2
	xxspltw		vs23,	vs25,	3


	addi		BO,	BO,	32


	xvmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs36,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs37,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs38,	vs1,	vs10		// a1_r*b1_r, a1_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs39,	vs1,	vs11		// a1_r*b1_i, a1_i*b1_i, a1_r*b1_i, a1_i*b1_i

	xvmaddasp	vs40,	vs0,	vs12		// a0_r*b2_r, a0_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs41,	vs0,	vs13		// a0_r*b2_i, a0_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmaddasp	vs42,	vs1,	vs12		// a1_r*b2_r, a1_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs43,	vs1,	vs13		// a1_r*b2_i, a1_i*b2_i, a1_r*b2_i, a1_i*b2_i

	xvmaddasp	vs44,	vs0,	vs14		// a0_r*b3_r, a0_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs45,	vs0,	vs15		// a0_r*b3_i, a0_i*b3_i, a1_r*b3_i, a1_i*b3_i
	xvmaddasp	vs46,	vs1,	vs14		// a1_r*b3_r, a1_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs47,	vs1,	vs15		// a1_r*b3_i, a1_i*b3_i, a1_r*b3_i, a1_i*b3_i


.endm

.macro KERNEL4x4_2


	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3


	addi		AO,	AO,	32

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3

	lxvw4x		vs25,	o16,	BO		//  load b2, b3

	xxspltw		vs12,	vs25,	0
	xxspltw		vs13,	vs25,	1
	xxspltw		vs14,	vs25,	2
	xxspltw		vs15,	vs25,	3


	addi		BO,	BO,	32


	xvmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r, a4_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs4,	vs17		// a4_r*b0_i, a4_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs5,	vs16		// a5_r*b0_r, a5_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs5,	vs17		// a5_r*b0_i, a5_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs36,	vs4,	vs18		// a4_r*b1_r, a4_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs37,	vs4,	vs19		// a4_r*b1_i, a4_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs38,	vs5,	vs18		// a5_r*b1_r, a5_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs39,	vs5,	vs19		// a5_r*b1_i, a5_i*b1_i, a1_r*b1_i, a1_i*b1_i

	xvmaddasp	vs40,	vs4,	vs20		// a4_r*b2_r, a4_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs41,	vs4,	vs21		// a4_r*b2_i, a4_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmaddasp	vs42,	vs5,	vs20		// a5_r*b2_r, a5_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs43,	vs5,	vs21		// a5_r*b2_i, a5_i*b2_i, a1_r*b2_i, a1_i*b2_i

	xvmaddasp	vs44,	vs4,	vs22		// a4_r*b3_r, a4_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs45,	vs4,	vs23		// a4_r*b3_i, a4_i*b3_i, a1_r*b3_i, a1_i*b3_i
	xvmaddasp	vs46,	vs5,	vs22		// a5_r*b3_r, a5_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs47,	vs5,	vs23		// a5_r*b3_i, a5_i*b3_i, a1_r*b3_i, a1_i*b3_i


.endm

.macro KERNEL4x4_E2


	xvmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r, a4_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs4,	vs17		// a4_r*b0_i, a4_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs5,	vs16		// a5_r*b0_r, a5_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs5,	vs17		// a5_r*b0_i, a5_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs36,	vs4,	vs18		// a4_r*b1_r, a4_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs37,	vs4,	vs19		// a4_r*b1_i, a4_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs38,	vs5,	vs18		// a5_r*b1_r, a5_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs39,	vs5,	vs19		// a5_r*b1_i, a5_i*b1_i, a1_r*b1_i, a1_i*b1_i

	xvmaddasp	vs40,	vs4,	vs20		// a4_r*b2_r, a4_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs41,	vs4,	vs21		// a4_r*b2_i, a4_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmaddasp	vs42,	vs5,	vs20		// a5_r*b2_r, a5_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs43,	vs5,	vs21		// a5_r*b2_i, a5_i*b2_i, a1_r*b2_i, a1_i*b2_i

	xvmaddasp	vs44,	vs4,	vs22		// a4_r*b3_r, a4_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs45,	vs4,	vs23		// a4_r*b3_i, a4_i*b3_i, a1_r*b3_i, a1_i*b3_i
	xvmaddasp	vs46,	vs5,	vs22		// a5_r*b3_r, a5_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs47,	vs5,	vs23		// a5_r*b3_i, a5_i*b3_i, a1_r*b3_i, a1_i*b3_i


.endm

.macro KERNEL4x4_SUBI1


	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3


	addi		AO,	AO,	32

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3

	lxvw4x		vs25,	o16,	BO		//  load b2, b3

	xxspltw		vs12,	vs25,	0
	xxspltw		vs13,	vs25,	1
	xxspltw		vs14,	vs25,	2
	xxspltw		vs15,	vs25,	3


	addi		BO,	BO,	32


	xvmulsp		vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmulsp		vs36,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs37,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmulsp		vs38,	vs1,	vs10		// a1_r*b1_r, a1_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs39,	vs1,	vs11		// a1_r*b1_i, a1_i*b1_i, a1_r*b1_i, a1_i*b1_i

	xvmulsp		vs40,	vs0,	vs12		// a0_r*b2_r, a0_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmulsp		vs41,	vs0,	vs13		// a0_r*b2_i, a0_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmulsp		vs42,	vs1,	vs12		// a1_r*b2_r, a1_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmulsp		vs43,	vs1,	vs13		// a1_r*b2_i, a1_i*b2_i, a1_r*b2_i, a1_i*b2_i

	xvmulsp		vs44,	vs0,	vs14		// a0_r*b3_r, a0_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmulsp		vs45,	vs0,	vs15		// a0_r*b3_i, a0_i*b3_i, a1_r*b3_i, a1_i*b3_i
	xvmulsp		vs46,	vs1,	vs14		// a1_r*b3_r, a1_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmulsp		vs47,	vs1,	vs15		// a1_r*b3_i, a1_i*b3_i, a1_r*b3_i, a1_i*b3_i


.endm

.macro KERNEL4x4_SUB1


	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3


	addi		AO,	AO,	32

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3

	lxvw4x		vs25,	o16,	BO		//  load b2, b3

	xxspltw		vs12,	vs25,	0
	xxspltw		vs13,	vs25,	1
	xxspltw		vs14,	vs25,	2
	xxspltw		vs15,	vs25,	3


	addi		BO,	BO,	32


	xvmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs36,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs37,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs38,	vs1,	vs10		// a1_r*b1_r, a1_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs39,	vs1,	vs11		// a1_r*b1_i, a1_i*b1_i, a1_r*b1_i, a1_i*b1_i

	xvmaddasp	vs40,	vs0,	vs12		// a0_r*b2_r, a0_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs41,	vs0,	vs13		// a0_r*b2_i, a0_i*b2_i, a1_r*b2_i, a1_i*b2_i
	xvmaddasp	vs42,	vs1,	vs12		// a1_r*b2_r, a1_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs43,	vs1,	vs13		// a1_r*b2_i, a1_i*b2_i, a1_r*b2_i, a1_i*b2_i

	xvmaddasp	vs44,	vs0,	vs14		// a0_r*b3_r, a0_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs45,	vs0,	vs15		// a0_r*b3_i, a0_i*b3_i, a1_r*b3_i, a1_i*b3_i
	xvmaddasp	vs46,	vs1,	vs14		// a1_r*b3_r, a1_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs47,	vs1,	vs15		// a1_r*b3_i, a1_i*b3_i, a1_r*b3_i, a1_i*b3_i


.endm

.macro SAVE4x4

	mr		T1,	CO

// N=0

	mr		T2,	T1

// N=0 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs33,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=0 M=2

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs34,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs35,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC


// N=1

	mr		T2,	T1

// N=1 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs36,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs37,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=1 M=2

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs38,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs39,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC


// N=2

	mr		T2,	T1

// N=2 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs40,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs41,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=2 M=2

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs42,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs43,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC


// N=3

	mr		T2,	T1

// N=3 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs44,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs45,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=3 M=2

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs46,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs47,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC

	addi		CO,	CO,	32

.endm


/**********************************************************************************************
* Macros for N=4 and M=2
**********************************************************************************************/

.macro LOAD4x2_1

	lxvw4x		vs0,	o0,	AO		// load a0, a1


	addi		AO,	AO,	16

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3

	lxvw4x		vs25,	o16,	BO		//  load b2, b3

	xxspltw		vs12,	vs25,	0
	xxspltw		vs13,	vs25,	1
	xxspltw		vs14,	vs25,	2
	xxspltw		vs15,	vs25,	3


	addi		BO,	BO,	32

.endm

.macro KERNEL4x2_I1


	lxvw4x		vs4,	o0,	AO		// load a0, a1


	addi		AO,	AO,	16

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs16,	vs24,	0
	xxspltw		vs17,	vs24,	1
	xxspltw		vs18,	vs24,	2
	xxspltw		vs19,	vs24,	3

	lxvw4x		vs25,	o16,	BO		//  load b2, b3

	xxspltw		vs20,	vs25,	0
	xxspltw		vs21,	vs25,	1
	xxspltw		vs22,	vs25,	2
	xxspltw		vs23,	vs25,	3


	addi		BO,	BO,	32


	xvmulsp		vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmulsp		vs34,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs35,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i

	xvmulsp		vs36,	vs0,	vs12		// a0_r*b2_r, a0_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmulsp		vs37,	vs0,	vs13		// a0_r*b2_i, a0_i*b2_i, a1_r*b2_i, a1_i*b2_i

	xvmulsp		vs38,	vs0,	vs14		// a0_r*b3_r, a0_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmulsp		vs39,	vs0,	vs15		// a0_r*b3_i, a0_i*b3_i, a1_r*b3_i, a1_i*b3_i


.endm

.macro KERNEL4x2_1


	lxvw4x		vs4,	o0,	AO		// load a0, a1


	addi		AO,	AO,	16

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs16,	vs24,	0
	xxspltw		vs17,	vs24,	1
	xxspltw		vs18,	vs24,	2
	xxspltw		vs19,	vs24,	3

	lxvw4x		vs25,	o16,	BO		//  load b2, b3

	xxspltw		vs20,	vs25,	0
	xxspltw		vs21,	vs25,	1
	xxspltw		vs22,	vs25,	2
	xxspltw		vs23,	vs25,	3


	addi		BO,	BO,	32


	xvmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs34,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs35,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i

	xvmaddasp	vs36,	vs0,	vs12		// a0_r*b2_r, a0_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs37,	vs0,	vs13		// a0_r*b2_i, a0_i*b2_i, a1_r*b2_i, a1_i*b2_i

	xvmaddasp	vs38,	vs0,	vs14		// a0_r*b3_r, a0_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs39,	vs0,	vs15		// a0_r*b3_i, a0_i*b3_i, a1_r*b3_i, a1_i*b3_i


.endm

.macro KERNEL4x2_2


	lxvw4x		vs0,	o0,	AO		// load a0, a1


	addi		AO,	AO,	16

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3

	lxvw4x		vs25,	o16,	BO		//  load b2, b3

	xxspltw		vs12,	vs25,	0
	xxspltw		vs13,	vs25,	1
	xxspltw		vs14,	vs25,	2
	xxspltw		vs15,	vs25,	3


	addi		BO,	BO,	32


	xvmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r, a4_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs4,	vs17		// a4_r*b0_i, a4_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs34,	vs4,	vs18		// a4_r*b1_r, a4_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs35,	vs4,	vs19		// a4_r*b1_i, a4_i*b1_i, a1_r*b1_i, a1_i*b1_i

	xvmaddasp	vs36,	vs4,	vs20		// a4_r*b2_r, a4_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs37,	vs4,	vs21		// a4_r*b2_i, a4_i*b2_i, a1_r*b2_i, a1_i*b2_i

	xvmaddasp	vs38,	vs4,	vs22		// a4_r*b3_r, a4_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs39,	vs4,	vs23		// a4_r*b3_i, a4_i*b3_i, a1_r*b3_i, a1_i*b3_i


.endm

.macro KERNEL4x2_E2


	xvmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r, a4_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs4,	vs17		// a4_r*b0_i, a4_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs34,	vs4,	vs18		// a4_r*b1_r, a4_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs35,	vs4,	vs19		// a4_r*b1_i, a4_i*b1_i, a1_r*b1_i, a1_i*b1_i

	xvmaddasp	vs36,	vs4,	vs20		// a4_r*b2_r, a4_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs37,	vs4,	vs21		// a4_r*b2_i, a4_i*b2_i, a1_r*b2_i, a1_i*b2_i

	xvmaddasp	vs38,	vs4,	vs22		// a4_r*b3_r, a4_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs39,	vs4,	vs23		// a4_r*b3_i, a4_i*b3_i, a1_r*b3_i, a1_i*b3_i


.endm

.macro KERNEL4x2_SUBI1


	lxvw4x		vs0,	o0,	AO		// load a0, a1


	addi		AO,	AO,	16

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3

	lxvw4x		vs25,	o16,	BO		//  load b2, b3

	xxspltw		vs12,	vs25,	0
	xxspltw		vs13,	vs25,	1
	xxspltw		vs14,	vs25,	2
	xxspltw		vs15,	vs25,	3


	addi		BO,	BO,	32


	xvmulsp		vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmulsp		vs34,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs35,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i

	xvmulsp		vs36,	vs0,	vs12		// a0_r*b2_r, a0_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmulsp		vs37,	vs0,	vs13		// a0_r*b2_i, a0_i*b2_i, a1_r*b2_i, a1_i*b2_i

	xvmulsp		vs38,	vs0,	vs14		// a0_r*b3_r, a0_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmulsp		vs39,	vs0,	vs15		// a0_r*b3_i, a0_i*b3_i, a1_r*b3_i, a1_i*b3_i


.endm

.macro KERNEL4x2_SUB1


	lxvw4x		vs0,	o0,	AO		// load a0, a1


	addi		AO,	AO,	16

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3

	lxvw4x		vs25,	o16,	BO		//  load b2, b3

	xxspltw		vs12,	vs25,	0
	xxspltw		vs13,	vs25,	1
	xxspltw		vs14,	vs25,	2
	xxspltw		vs15,	vs25,	3


	addi		BO,	BO,	32


	xvmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs34,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs35,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i

	xvmaddasp	vs36,	vs0,	vs12		// a0_r*b2_r, a0_i*b2_r, a1_r*b2_r, a1_i*b2_r
	xvmaddasp	vs37,	vs0,	vs13		// a0_r*b2_i, a0_i*b2_i, a1_r*b2_i, a1_i*b2_i

	xvmaddasp	vs38,	vs0,	vs14		// a0_r*b3_r, a0_i*b3_r, a1_r*b3_r, a1_i*b3_r
	xvmaddasp	vs39,	vs0,	vs15		// a0_r*b3_i, a0_i*b3_i, a1_r*b3_i, a1_i*b3_i


.endm

.macro SAVE4x2

	mr		T1,	CO

// N=0

	mr		T2,	T1

// N=0 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs33,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC


// N=1

	mr		T2,	T1

// N=1 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs34,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs35,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC


// N=2

	mr		T2,	T1

// N=2 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs36,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs37,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC


// N=3

	mr		T2,	T1

// N=3 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs38,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs39,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC

	addi		CO,	CO,	16

.endm


/**********************************************************************************************
* Macros for N=4 and M=1
**********************************************************************************************/

.macro LOAD4x1_1

	lxsspx		vs0,	o0,	AO		// load a0_r
	lxsspx		vs1,	o4,	AO		// load a0_i

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1		//  load b0_r
	lxsspx		vs9,	o4,	T1		//  load b0_i

	addi		T1,	T1,8

	lxsspx		vs10,	o0,	T1		//  load b1_r
	lxsspx		vs11,	o4,	T1		//  load b1_i

	addi		T1,	T1,8

	lxsspx		vs12,	o0,	T1		//  load b2_r
	lxsspx		vs13,	o4,	T1		//  load b2_i

	addi		T1,	T1,8

	lxsspx		vs14,	o0,	T1		//  load b3_r
	lxsspx		vs15,	o4,	T1		//  load b3_i

	addi		BO,	BO,	32

.endm

.macro KERNEL4x1_I1


	lxsspx		vs4,	o0,	AO		// load a0_r
	lxsspx		vs5,	o4,	AO		// load a0_i

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1		//  load b0_r
	lxsspx		vs17,	o4,	T1		//  load b0_i

	addi		T1,	T1,8

	lxsspx		vs18,	o0,	T1		//  load b1_r
	lxsspx		vs19,	o4,	T1		//  load b1_i

	addi		T1,	T1,8

	lxsspx		vs20,	o0,	T1		//  load b2_r
	lxsspx		vs21,	o4,	T1		//  load b2_i

	addi		T1,	T1,8

	lxsspx		vs22,	o0,	T1		//  load b3_r
	lxsspx		vs23,	o4,	T1		//  load b3_i

	addi		BO,	BO,	32


	xsmulsp		vs32,	vs0,	vs8		// a0_r*b0_r
	xsmulsp		vs33,	vs1,	vs9		// a0_i*b0_i
	xsmulsp		vs34,	vs0,	vs9		// a0_r*b0_i
	xsmulsp		vs35,	vs1,	vs8		// a0_i*b0_r

	xsmulsp		vs36,	vs0,	vs10		// a0_r*b1_r
	xsmulsp		vs37,	vs1,	vs11		// a0_i*b1_i
	xsmulsp		vs38,	vs0,	vs11		// a0_r*b1_i
	xsmulsp		vs39,	vs1,	vs10		// a0_i*b1_r

	xsmulsp		vs40,	vs0,	vs12		// a0_r*b2_r
	xsmulsp		vs41,	vs1,	vs13		// a0_i*b2_i
	xsmulsp		vs42,	vs0,	vs13		// a0_r*b2_i
	xsmulsp		vs43,	vs1,	vs12		// a0_i*b2_r

	xsmulsp		vs44,	vs0,	vs14		// a0_r*b3_r
	xsmulsp		vs45,	vs1,	vs15		// a0_i*b3_i
	xsmulsp		vs46,	vs0,	vs15		// a0_r*b3_i
	xsmulsp		vs47,	vs1,	vs14		// a0_i*b3_r


.endm

.macro KERNEL4x1_1


	lxsspx		vs4,	o0,	AO		// load a0_r
	lxsspx		vs5,	o4,	AO		// load a0_i

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1		//  load b0_r
	lxsspx		vs17,	o4,	T1		//  load b0_i

	addi		T1,	T1,8

	lxsspx		vs18,	o0,	T1		//  load b1_r
	lxsspx		vs19,	o4,	T1		//  load b1_i

	addi		T1,	T1,8

	lxsspx		vs20,	o0,	T1		//  load b2_r
	lxsspx		vs21,	o4,	T1		//  load b2_i

	addi		T1,	T1,8

	lxsspx		vs22,	o0,	T1		//  load b3_r
	lxsspx		vs23,	o4,	T1		//  load b3_i

	addi		BO,	BO,	32


	xsmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r
	xsmaddasp	vs33,	vs1,	vs9		// a0_i*b0_i
	xsmaddasp	vs34,	vs0,	vs9		// a0_r*b0_i
	xsmaddasp	vs35,	vs1,	vs8		// a0_i*b0_r

	xsmaddasp	vs36,	vs0,	vs10		// a0_r*b1_r
	xsmaddasp	vs37,	vs1,	vs11		// a0_i*b1_i
	xsmaddasp	vs38,	vs0,	vs11		// a0_r*b1_i
	xsmaddasp	vs39,	vs1,	vs10		// a0_i*b1_r

	xsmaddasp	vs40,	vs0,	vs12		// a0_r*b2_r
	xsmaddasp	vs41,	vs1,	vs13		// a0_i*b2_i
	xsmaddasp	vs42,	vs0,	vs13		// a0_r*b2_i
	xsmaddasp	vs43,	vs1,	vs12		// a0_i*b2_r

	xsmaddasp	vs44,	vs0,	vs14		// a0_r*b3_r
	xsmaddasp	vs45,	vs1,	vs15		// a0_i*b3_i
	xsmaddasp	vs46,	vs0,	vs15		// a0_r*b3_i
	xsmaddasp	vs47,	vs1,	vs14		// a0_i*b3_r


.endm

.macro KERNEL4x1_2


	lxsspx		vs0,	o0,	AO		// load a0_r
	lxsspx		vs1,	o4,	AO		// load a0_i

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1		//  load b0_r
	lxsspx		vs9,	o4,	T1		//  load b0_i

	addi		T1,	T1,8

	lxsspx		vs10,	o0,	T1		//  load b1_r
	lxsspx		vs11,	o4,	T1		//  load b1_i

	addi		T1,	T1,8

	lxsspx		vs12,	o0,	T1		//  load b2_r
	lxsspx		vs13,	o4,	T1		//  load b2_i

	addi		T1,	T1,8

	lxsspx		vs14,	o0,	T1		//  load b3_r
	lxsspx		vs15,	o4,	T1		//  load b3_i

	addi		BO,	BO,	32


	xsmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r
	xsmaddasp	vs33,	vs5,	vs17		// a4_i*b0_i
	xsmaddasp	vs34,	vs4,	vs17		// a4_r*b0_i
	xsmaddasp	vs35,	vs5,	vs16		// a4_i*b0_r

	xsmaddasp	vs36,	vs4,	vs18		// a4_r*b1_r
	xsmaddasp	vs37,	vs5,	vs19		// a4_i*b1_i
	xsmaddasp	vs38,	vs4,	vs19		// a4_r*b1_i
	xsmaddasp	vs39,	vs5,	vs18		// a4_i*b1_r

	xsmaddasp	vs40,	vs4,	vs20		// a4_r*b2_r
	xsmaddasp	vs41,	vs5,	vs21		// a4_i*b2_i
	xsmaddasp	vs42,	vs4,	vs21		// a4_r*b2_i
	xsmaddasp	vs43,	vs5,	vs20		// a4_i*b2_r

	xsmaddasp	vs44,	vs4,	vs22		// a4_r*b3_r
	xsmaddasp	vs45,	vs5,	vs23		// a4_i*b3_i
	xsmaddasp	vs46,	vs4,	vs23		// a4_r*b3_i
	xsmaddasp	vs47,	vs5,	vs22		// a4_i*b3_r


.endm

.macro KERNEL4x1_E2


	xsmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r
	xsmaddasp	vs33,	vs5,	vs17		// a4_i*b0_i
	xsmaddasp	vs34,	vs4,	vs17		// a4_r*b0_i
	xsmaddasp	vs35,	vs5,	vs16		// a4_i*b0_r

	xsmaddasp	vs36,	vs4,	vs18		// a4_r*b1_r
	xsmaddasp	vs37,	vs5,	vs19		// a4_i*b1_i
	xsmaddasp	vs38,	vs4,	vs19		// a4_r*b1_i
	xsmaddasp	vs39,	vs5,	vs18		// a4_i*b1_r

	xsmaddasp	vs40,	vs4,	vs20		// a4_r*b2_r
	xsmaddasp	vs41,	vs5,	vs21		// a4_i*b2_i
	xsmaddasp	vs42,	vs4,	vs21		// a4_r*b2_i
	xsmaddasp	vs43,	vs5,	vs20		// a4_i*b2_r

	xsmaddasp	vs44,	vs4,	vs22		// a4_r*b3_r
	xsmaddasp	vs45,	vs5,	vs23		// a4_i*b3_i
	xsmaddasp	vs46,	vs4,	vs23		// a4_r*b3_i
	xsmaddasp	vs47,	vs5,	vs22		// a4_i*b3_r


.endm

.macro KERNEL4x1_SUBI1


	lxsspx		vs0,	o0,	AO		// load a0_r
	lxsspx		vs1,	o4,	AO		// load a0_i

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1		//  load b0_r
	lxsspx		vs9,	o4,	T1		//  load b0_i

	addi		T1,	T1,8

	lxsspx		vs10,	o0,	T1		//  load b1_r
	lxsspx		vs11,	o4,	T1		//  load b1_i

	addi		T1,	T1,8

	lxsspx		vs12,	o0,	T1		//  load b2_r
	lxsspx		vs13,	o4,	T1		//  load b2_i

	addi		T1,	T1,8

	lxsspx		vs14,	o0,	T1		//  load b3_r
	lxsspx		vs15,	o4,	T1		//  load b3_i

	addi		BO,	BO,	32


	xsmulsp		vs32,	vs0,	vs8		// a0_r*b0_r
	xsmulsp		vs33,	vs1,	vs9		// a0_i*b0_i
	xsmulsp		vs34,	vs0,	vs9		// a0_r*b0_i
	xsmulsp		vs35,	vs1,	vs8		// a0_i*b0_r

	xsmulsp		vs36,	vs0,	vs10		// a0_r*b1_r
	xsmulsp		vs37,	vs1,	vs11		// a0_i*b1_i
	xsmulsp		vs38,	vs0,	vs11		// a0_r*b1_i
	xsmulsp		vs39,	vs1,	vs10		// a0_i*b1_r

	xsmulsp		vs40,	vs0,	vs12		// a0_r*b2_r
	xsmulsp		vs41,	vs1,	vs13		// a0_i*b2_i
	xsmulsp		vs42,	vs0,	vs13		// a0_r*b2_i
	xsmulsp		vs43,	vs1,	vs12		// a0_i*b2_r

	xsmulsp		vs44,	vs0,	vs14		// a0_r*b3_r
	xsmulsp		vs45,	vs1,	vs15		// a0_i*b3_i
	xsmulsp		vs46,	vs0,	vs15		// a0_r*b3_i
	xsmulsp		vs47,	vs1,	vs14		// a0_i*b3_r


.endm

.macro KERNEL4x1_SUB1


	lxsspx		vs0,	o0,	AO		// load a0_r
	lxsspx		vs1,	o4,	AO		// load a0_i

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1		//  load b0_r
	lxsspx		vs9,	o4,	T1		//  load b0_i

	addi		T1,	T1,8

	lxsspx		vs10,	o0,	T1		//  load b1_r
	lxsspx		vs11,	o4,	T1		//  load b1_i

	addi		T1,	T1,8

	lxsspx		vs12,	o0,	T1		//  load b2_r
	lxsspx		vs13,	o4,	T1		//  load b2_i

	addi		T1,	T1,8

	lxsspx		vs14,	o0,	T1		//  load b3_r
	lxsspx		vs15,	o4,	T1		//  load b3_i

	addi		BO,	BO,	32


	xsmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r
	xsmaddasp	vs33,	vs1,	vs9		// a0_i*b0_i
	xsmaddasp	vs34,	vs0,	vs9		// a0_r*b0_i
	xsmaddasp	vs35,	vs1,	vs8		// a0_i*b0_r

	xsmaddasp	vs36,	vs0,	vs10		// a0_r*b1_r
	xsmaddasp	vs37,	vs1,	vs11		// a0_i*b1_i
	xsmaddasp	vs38,	vs0,	vs11		// a0_r*b1_i
	xsmaddasp	vs39,	vs1,	vs10		// a0_i*b1_r

	xsmaddasp	vs40,	vs0,	vs12		// a0_r*b2_r
	xsmaddasp	vs41,	vs1,	vs13		// a0_i*b2_i
	xsmaddasp	vs42,	vs0,	vs13		// a0_r*b2_i
	xsmaddasp	vs43,	vs1,	vs12		// a0_i*b2_r

	xsmaddasp	vs44,	vs0,	vs14		// a0_r*b3_r
	xsmaddasp	vs45,	vs1,	vs15		// a0_i*b3_i
	xsmaddasp	vs46,	vs0,	vs15		// a0_r*b3_i
	xsmaddasp	vs47,	vs1,	vs14		// a0_i*b3_r


.endm

.macro SAVE4x1

	mr		T1,	CO

// N=0

	mr		T2,	T1

// N=0 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5

#ifndef TRMMKERNEL
	lxsspx		vs0,	o0,	T2	// load c0_r
	lxsspx		vs1,	o4,	T2	// load c0_i
#else
	xxlxor		vs0,	vs0,	vs0
	xxlxor		vs1,	vs1,	vs1
#endif

	XSFADD_R1	vs4,	vs4,	vs32		// add a0_r * b0_r
	XSFADD_I1	vs5,	vs5,	vs35		// add a0_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs33		// add a0_i * b0_i
	XSFADD_I2	vs5,	vs5,	vs34		// add a0_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsaddsp		vs0,	vs0,	vs20
	xsaddsp		vs1,	vs1,	vs21


	stxsspx		vs0,	o0,	T2	// store c0_r
	stxsspx		vs1,	o4,	T2	// store c0_i

	addi		T2,	T2,	8
	add		T1,	T1,	LDC


// N=1

	mr		T2,	T1

// N=1 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5

#ifndef TRMMKERNEL
	lxsspx		vs0,	o0,	T2	// load c0_r
	lxsspx		vs1,	o4,	T2	// load c0_i
#else
	xxlxor		vs0,	vs0,	vs0
	xxlxor		vs1,	vs1,	vs1
#endif

	XSFADD_R1	vs4,	vs4,	vs36		// add a0_r * b0_r
	XSFADD_I1	vs5,	vs5,	vs39		// add a0_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs37		// add a0_i * b0_i
	XSFADD_I2	vs5,	vs5,	vs38		// add a0_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsaddsp		vs0,	vs0,	vs20
	xsaddsp		vs1,	vs1,	vs21


	stxsspx		vs0,	o0,	T2	// store c0_r
	stxsspx		vs1,	o4,	T2	// store c0_i

	addi		T2,	T2,	8
	add		T1,	T1,	LDC


// N=2

	mr		T2,	T1

// N=2 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5

#ifndef TRMMKERNEL
	lxsspx		vs0,	o0,	T2	// load c0_r
	lxsspx		vs1,	o4,	T2	// load c0_i
#else
	xxlxor		vs0,	vs0,	vs0
	xxlxor		vs1,	vs1,	vs1
#endif

	XSFADD_R1	vs4,	vs4,	vs40		// add a0_r * b0_r
	XSFADD_I1	vs5,	vs5,	vs43		// add a0_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs41		// add a0_i * b0_i
	XSFADD_I2	vs5,	vs5,	vs42		// add a0_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsaddsp		vs0,	vs0,	vs20
	xsaddsp		vs1,	vs1,	vs21


	stxsspx		vs0,	o0,	T2	// store c0_r
	stxsspx		vs1,	o4,	T2	// store c0_i

	addi		T2,	T2,	8
	add		T1,	T1,	LDC


// N=3

	mr		T2,	T1

// N=3 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5

#ifndef TRMMKERNEL
	lxsspx		vs0,	o0,	T2	// load c0_r
	lxsspx		vs1,	o4,	T2	// load c0_i
#else
	xxlxor		vs0,	vs0,	vs0
	xxlxor		vs1,	vs1,	vs1
#endif

	XSFADD_R1	vs4,	vs4,	vs44		// add a0_r * b0_r
	XSFADD_I1	vs5,	vs5,	vs47		// add a0_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs45		// add a0_i * b0_i
	XSFADD_I2	vs5,	vs5,	vs46		// add a0_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsaddsp		vs0,	vs0,	vs20
	xsaddsp		vs1,	vs1,	vs21


	stxsspx		vs0,	o0,	T2	// store c0_r
	stxsspx		vs1,	o4,	T2	// store c0_i

	addi		T2,	T2,	8
	add		T1,	T1,	LDC

	addi		CO,	CO,	8

.endm


/**********************************************************************************************
* Macros for N=2 and M=8
**********************************************************************************************/

.macro LOAD2x8_1

	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3

	lxvw4x		vs2,	o32,	AO		// load a4, a5

	lxvw4x		vs3,	o48,	AO		// load a6, a7


	addi		AO,	AO,	64

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	16

.endm

.macro KERNEL2x8_I1


	lxvw4x		vs4,	o0,	AO		// load a0, a1

	lxvw4x		vs5,	o16,	AO		// load a2, a3

	lxvw4x		vs6,	o32,	AO		// load a4, a5

	lxvw4x		vs7,	o48,	AO		// load a6, a7


	addi		AO,	AO,	64

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs16,	vs24,	0
	xxspltw		vs17,	vs24,	1
	xxspltw		vs18,	vs24,	2
	xxspltw		vs19,	vs24,	3


	addi		BO,	BO,	16


	xvmulsp		vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs36,	vs2,	vs8		// a2_r*b0_r, a2_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs37,	vs2,	vs9		// a2_r*b0_i, a2_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs38,	vs3,	vs8		// a3_r*b0_r, a3_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs39,	vs3,	vs9		// a3_r*b0_i, a3_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmulsp		vs40,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs41,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmulsp		vs42,	vs1,	vs10		// a1_r*b1_r, a1_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs43,	vs1,	vs11		// a1_r*b1_i, a1_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmulsp		vs44,	vs2,	vs10		// a2_r*b1_r, a2_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs45,	vs2,	vs11		// a2_r*b1_i, a2_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmulsp		vs46,	vs3,	vs10		// a3_r*b1_r, a3_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs47,	vs3,	vs11		// a3_r*b1_i, a3_i*b1_i, a1_r*b1_i, a1_i*b1_i


.endm

.macro KERNEL2x8_1


	lxvw4x		vs4,	o0,	AO		// load a0, a1

	lxvw4x		vs5,	o16,	AO		// load a2, a3

	lxvw4x		vs6,	o32,	AO		// load a4, a5

	lxvw4x		vs7,	o48,	AO		// load a6, a7


	addi		AO,	AO,	64

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs16,	vs24,	0
	xxspltw		vs17,	vs24,	1
	xxspltw		vs18,	vs24,	2
	xxspltw		vs19,	vs24,	3


	addi		BO,	BO,	16


	xvmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs36,	vs2,	vs8		// a2_r*b0_r, a2_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs37,	vs2,	vs9		// a2_r*b0_i, a2_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs38,	vs3,	vs8		// a3_r*b0_r, a3_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs39,	vs3,	vs9		// a3_r*b0_i, a3_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs40,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs41,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs42,	vs1,	vs10		// a1_r*b1_r, a1_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs43,	vs1,	vs11		// a1_r*b1_i, a1_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs44,	vs2,	vs10		// a2_r*b1_r, a2_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs45,	vs2,	vs11		// a2_r*b1_i, a2_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs46,	vs3,	vs10		// a3_r*b1_r, a3_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs47,	vs3,	vs11		// a3_r*b1_i, a3_i*b1_i, a1_r*b1_i, a1_i*b1_i


.endm

.macro KERNEL2x8_2


	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3

	lxvw4x		vs2,	o32,	AO		// load a4, a5

	lxvw4x		vs3,	o48,	AO		// load a6, a7


	addi		AO,	AO,	64

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	16


	xvmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r, a4_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs4,	vs17		// a4_r*b0_i, a4_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs5,	vs16		// a5_r*b0_r, a5_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs5,	vs17		// a5_r*b0_i, a5_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs36,	vs6,	vs16		// a6_r*b0_r, a6_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs37,	vs6,	vs17		// a6_r*b0_i, a6_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs38,	vs7,	vs16		// a7_r*b0_r, a7_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs39,	vs7,	vs17		// a7_r*b0_i, a7_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs40,	vs4,	vs18		// a4_r*b1_r, a4_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs41,	vs4,	vs19		// a4_r*b1_i, a4_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs42,	vs5,	vs18		// a5_r*b1_r, a5_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs43,	vs5,	vs19		// a5_r*b1_i, a5_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs44,	vs6,	vs18		// a6_r*b1_r, a6_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs45,	vs6,	vs19		// a6_r*b1_i, a6_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs46,	vs7,	vs18		// a7_r*b1_r, a7_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs47,	vs7,	vs19		// a7_r*b1_i, a7_i*b1_i, a1_r*b1_i, a1_i*b1_i


.endm

.macro KERNEL2x8_E2


	xvmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r, a4_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs4,	vs17		// a4_r*b0_i, a4_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs5,	vs16		// a5_r*b0_r, a5_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs5,	vs17		// a5_r*b0_i, a5_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs36,	vs6,	vs16		// a6_r*b0_r, a6_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs37,	vs6,	vs17		// a6_r*b0_i, a6_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs38,	vs7,	vs16		// a7_r*b0_r, a7_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs39,	vs7,	vs17		// a7_r*b0_i, a7_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs40,	vs4,	vs18		// a4_r*b1_r, a4_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs41,	vs4,	vs19		// a4_r*b1_i, a4_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs42,	vs5,	vs18		// a5_r*b1_r, a5_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs43,	vs5,	vs19		// a5_r*b1_i, a5_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs44,	vs6,	vs18		// a6_r*b1_r, a6_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs45,	vs6,	vs19		// a6_r*b1_i, a6_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs46,	vs7,	vs18		// a7_r*b1_r, a7_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs47,	vs7,	vs19		// a7_r*b1_i, a7_i*b1_i, a1_r*b1_i, a1_i*b1_i


.endm

.macro KERNEL2x8_SUBI1


	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3

	lxvw4x		vs2,	o32,	AO		// load a4, a5

	lxvw4x		vs3,	o48,	AO		// load a6, a7


	addi		AO,	AO,	64

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	16


	xvmulsp		vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs36,	vs2,	vs8		// a2_r*b0_r, a2_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs37,	vs2,	vs9		// a2_r*b0_i, a2_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs38,	vs3,	vs8		// a3_r*b0_r, a3_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs39,	vs3,	vs9		// a3_r*b0_i, a3_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmulsp		vs40,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs41,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmulsp		vs42,	vs1,	vs10		// a1_r*b1_r, a1_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs43,	vs1,	vs11		// a1_r*b1_i, a1_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmulsp		vs44,	vs2,	vs10		// a2_r*b1_r, a2_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs45,	vs2,	vs11		// a2_r*b1_i, a2_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmulsp		vs46,	vs3,	vs10		// a3_r*b1_r, a3_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs47,	vs3,	vs11		// a3_r*b1_i, a3_i*b1_i, a1_r*b1_i, a1_i*b1_i


.endm

.macro KERNEL2x8_SUB1


	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3

	lxvw4x		vs2,	o32,	AO		// load a4, a5

	lxvw4x		vs3,	o48,	AO		// load a6, a7


	addi		AO,	AO,	64

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	16


	xvmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs36,	vs2,	vs8		// a2_r*b0_r, a2_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs37,	vs2,	vs9		// a2_r*b0_i, a2_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs38,	vs3,	vs8		// a3_r*b0_r, a3_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs39,	vs3,	vs9		// a3_r*b0_i, a3_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs40,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs41,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs42,	vs1,	vs10		// a1_r*b1_r, a1_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs43,	vs1,	vs11		// a1_r*b1_i, a1_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs44,	vs2,	vs10		// a2_r*b1_r, a2_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs45,	vs2,	vs11		// a2_r*b1_i, a2_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs46,	vs3,	vs10		// a3_r*b1_r, a3_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs47,	vs3,	vs11		// a3_r*b1_i, a3_i*b1_i, a1_r*b1_i, a1_i*b1_i


.endm

.macro SAVE2x8

	mr		T1,	CO

// N=0

	mr		T2,	T1

// N=0 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs33,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=0 M=2

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs34,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs35,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=0 M=4

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs36,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs37,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=0 M=6

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs38,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs39,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC


// N=1

	mr		T2,	T1

// N=1 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs40,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs41,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=1 M=2

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs42,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs43,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=1 M=4

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs44,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs45,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=1 M=6

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs46,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs47,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC

	addi		CO,	CO,	64

.endm


/**********************************************************************************************
* Macros for N=2 and M=4
**********************************************************************************************/

.macro LOAD2x4_1

	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3


	addi		AO,	AO,	32

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	16

.endm

.macro KERNEL2x4_I1


	lxvw4x		vs4,	o0,	AO		// load a0, a1

	lxvw4x		vs5,	o16,	AO		// load a2, a3


	addi		AO,	AO,	32

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs16,	vs24,	0
	xxspltw		vs17,	vs24,	1
	xxspltw		vs18,	vs24,	2
	xxspltw		vs19,	vs24,	3


	addi		BO,	BO,	16


	xvmulsp		vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmulsp		vs36,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs37,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmulsp		vs38,	vs1,	vs10		// a1_r*b1_r, a1_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs39,	vs1,	vs11		// a1_r*b1_i, a1_i*b1_i, a1_r*b1_i, a1_i*b1_i


.endm

.macro KERNEL2x4_1


	lxvw4x		vs4,	o0,	AO		// load a0, a1

	lxvw4x		vs5,	o16,	AO		// load a2, a3


	addi		AO,	AO,	32

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs16,	vs24,	0
	xxspltw		vs17,	vs24,	1
	xxspltw		vs18,	vs24,	2
	xxspltw		vs19,	vs24,	3


	addi		BO,	BO,	16


	xvmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs36,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs37,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs38,	vs1,	vs10		// a1_r*b1_r, a1_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs39,	vs1,	vs11		// a1_r*b1_i, a1_i*b1_i, a1_r*b1_i, a1_i*b1_i


.endm

.macro KERNEL2x4_2


	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3


	addi		AO,	AO,	32

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	16


	xvmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r, a4_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs4,	vs17		// a4_r*b0_i, a4_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs5,	vs16		// a5_r*b0_r, a5_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs5,	vs17		// a5_r*b0_i, a5_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs36,	vs4,	vs18		// a4_r*b1_r, a4_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs37,	vs4,	vs19		// a4_r*b1_i, a4_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs38,	vs5,	vs18		// a5_r*b1_r, a5_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs39,	vs5,	vs19		// a5_r*b1_i, a5_i*b1_i, a1_r*b1_i, a1_i*b1_i


.endm

.macro KERNEL2x4_E2


	xvmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r, a4_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs4,	vs17		// a4_r*b0_i, a4_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs5,	vs16		// a5_r*b0_r, a5_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs5,	vs17		// a5_r*b0_i, a5_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs36,	vs4,	vs18		// a4_r*b1_r, a4_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs37,	vs4,	vs19		// a4_r*b1_i, a4_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs38,	vs5,	vs18		// a5_r*b1_r, a5_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs39,	vs5,	vs19		// a5_r*b1_i, a5_i*b1_i, a1_r*b1_i, a1_i*b1_i


.endm

.macro KERNEL2x4_SUBI1


	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3


	addi		AO,	AO,	32

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	16


	xvmulsp		vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmulsp		vs36,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs37,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmulsp		vs38,	vs1,	vs10		// a1_r*b1_r, a1_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs39,	vs1,	vs11		// a1_r*b1_i, a1_i*b1_i, a1_r*b1_i, a1_i*b1_i


.endm

.macro KERNEL2x4_SUB1


	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3


	addi		AO,	AO,	32

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	16


	xvmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs36,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs37,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i
	xvmaddasp	vs38,	vs1,	vs10		// a1_r*b1_r, a1_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs39,	vs1,	vs11		// a1_r*b1_i, a1_i*b1_i, a1_r*b1_i, a1_i*b1_i


.endm

.macro SAVE2x4

	mr		T1,	CO

// N=0

	mr		T2,	T1

// N=0 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs33,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=0 M=2

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs34,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs35,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC


// N=1

	mr		T2,	T1

// N=1 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs36,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs37,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=1 M=2

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs38,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs39,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC

	addi		CO,	CO,	32

.endm


/**********************************************************************************************
* Macros for N=2 and M=2
**********************************************************************************************/

.macro LOAD2x2_1

	lxvw4x		vs0,	o0,	AO		// load a0, a1


	addi		AO,	AO,	16

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	16

.endm

.macro KERNEL2x2_I1


	lxvw4x		vs4,	o0,	AO		// load a0, a1


	addi		AO,	AO,	16

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs16,	vs24,	0
	xxspltw		vs17,	vs24,	1
	xxspltw		vs18,	vs24,	2
	xxspltw		vs19,	vs24,	3


	addi		BO,	BO,	16


	xvmulsp		vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmulsp		vs34,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs35,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i


.endm

.macro KERNEL2x2_1


	lxvw4x		vs4,	o0,	AO		// load a0, a1


	addi		AO,	AO,	16

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs16,	vs24,	0
	xxspltw		vs17,	vs24,	1
	xxspltw		vs18,	vs24,	2
	xxspltw		vs19,	vs24,	3


	addi		BO,	BO,	16


	xvmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs34,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs35,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i


.endm

.macro KERNEL2x2_2


	lxvw4x		vs0,	o0,	AO		// load a0, a1


	addi		AO,	AO,	16

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	16


	xvmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r, a4_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs4,	vs17		// a4_r*b0_i, a4_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs34,	vs4,	vs18		// a4_r*b1_r, a4_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs35,	vs4,	vs19		// a4_r*b1_i, a4_i*b1_i, a1_r*b1_i, a1_i*b1_i


.endm

.macro KERNEL2x2_E2


	xvmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r, a4_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs4,	vs17		// a4_r*b0_i, a4_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs34,	vs4,	vs18		// a4_r*b1_r, a4_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs35,	vs4,	vs19		// a4_r*b1_i, a4_i*b1_i, a1_r*b1_i, a1_i*b1_i


.endm

.macro KERNEL2x2_SUBI1


	lxvw4x		vs0,	o0,	AO		// load a0, a1


	addi		AO,	AO,	16

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	16


	xvmulsp		vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmulsp		vs34,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmulsp		vs35,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i


.endm

.macro KERNEL2x2_SUB1


	lxvw4x		vs0,	o0,	AO		// load a0, a1


	addi		AO,	AO,	16

	lxvw4x		vs24,	o0,	BO		//  load b0, b1

	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	16


	xvmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i

	xvmaddasp	vs34,	vs0,	vs10		// a0_r*b1_r, a0_i*b1_r, a1_r*b1_r, a1_i*b1_r
	xvmaddasp	vs35,	vs0,	vs11		// a0_r*b1_i, a0_i*b1_i, a1_r*b1_i, a1_i*b1_i


.endm

.macro SAVE2x2

	mr		T1,	CO

// N=0

	mr		T2,	T1

// N=0 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs33,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC


// N=1

	mr		T2,	T1

// N=1 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs34,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs35,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC

	addi		CO,	CO,	16

.endm


/**********************************************************************************************
* Macros for N=2 and M=1
**********************************************************************************************/

.macro LOAD2x1_1

	lxsspx		vs0,	o0,	AO		// load a0_r
	lxsspx		vs1,	o4,	AO		// load a0_i

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1		//  load b0_r
	lxsspx		vs9,	o4,	T1		//  load b0_i

	addi		T1,	T1,8

	lxsspx		vs10,	o0,	T1		//  load b1_r
	lxsspx		vs11,	o4,	T1		//  load b1_i

	addi		BO,	BO,	16

.endm

.macro KERNEL2x1_I1


	lxsspx		vs4,	o0,	AO		// load a0_r
	lxsspx		vs5,	o4,	AO		// load a0_i

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1		//  load b0_r
	lxsspx		vs17,	o4,	T1		//  load b0_i

	addi		T1,	T1,8

	lxsspx		vs18,	o0,	T1		//  load b1_r
	lxsspx		vs19,	o4,	T1		//  load b1_i

	addi		BO,	BO,	16


	xsmulsp		vs32,	vs0,	vs8		// a0_r*b0_r
	xsmulsp		vs33,	vs1,	vs9		// a0_i*b0_i
	xsmulsp		vs34,	vs0,	vs9		// a0_r*b0_i
	xsmulsp		vs35,	vs1,	vs8		// a0_i*b0_r

	xsmulsp		vs36,	vs0,	vs10		// a0_r*b1_r
	xsmulsp		vs37,	vs1,	vs11		// a0_i*b1_i
	xsmulsp		vs38,	vs0,	vs11		// a0_r*b1_i
	xsmulsp		vs39,	vs1,	vs10		// a0_i*b1_r


.endm

.macro KERNEL2x1_1


	lxsspx		vs4,	o0,	AO		// load a0_r
	lxsspx		vs5,	o4,	AO		// load a0_i

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1		//  load b0_r
	lxsspx		vs17,	o4,	T1		//  load b0_i

	addi		T1,	T1,8

	lxsspx		vs18,	o0,	T1		//  load b1_r
	lxsspx		vs19,	o4,	T1		//  load b1_i

	addi		BO,	BO,	16


	xsmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r
	xsmaddasp	vs33,	vs1,	vs9		// a0_i*b0_i
	xsmaddasp	vs34,	vs0,	vs9		// a0_r*b0_i
	xsmaddasp	vs35,	vs1,	vs8		// a0_i*b0_r

	xsmaddasp	vs36,	vs0,	vs10		// a0_r*b1_r
	xsmaddasp	vs37,	vs1,	vs11		// a0_i*b1_i
	xsmaddasp	vs38,	vs0,	vs11		// a0_r*b1_i
	xsmaddasp	vs39,	vs1,	vs10		// a0_i*b1_r


.endm

.macro KERNEL2x1_2


	lxsspx		vs0,	o0,	AO		// load a0_r
	lxsspx		vs1,	o4,	AO		// load a0_i

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1		//  load b0_r
	lxsspx		vs9,	o4,	T1		//  load b0_i

	addi		T1,	T1,8

	lxsspx		vs10,	o0,	T1		//  load b1_r
	lxsspx		vs11,	o4,	T1		//  load b1_i

	addi		BO,	BO,	16


	xsmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r
	xsmaddasp	vs33,	vs5,	vs17		// a4_i*b0_i
	xsmaddasp	vs34,	vs4,	vs17		// a4_r*b0_i
	xsmaddasp	vs35,	vs5,	vs16		// a4_i*b0_r

	xsmaddasp	vs36,	vs4,	vs18		// a4_r*b1_r
	xsmaddasp	vs37,	vs5,	vs19		// a4_i*b1_i
	xsmaddasp	vs38,	vs4,	vs19		// a4_r*b1_i
	xsmaddasp	vs39,	vs5,	vs18		// a4_i*b1_r


.endm

.macro KERNEL2x1_E2


	xsmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r
	xsmaddasp	vs33,	vs5,	vs17		// a4_i*b0_i
	xsmaddasp	vs34,	vs4,	vs17		// a4_r*b0_i
	xsmaddasp	vs35,	vs5,	vs16		// a4_i*b0_r

	xsmaddasp	vs36,	vs4,	vs18		// a4_r*b1_r
	xsmaddasp	vs37,	vs5,	vs19		// a4_i*b1_i
	xsmaddasp	vs38,	vs4,	vs19		// a4_r*b1_i
	xsmaddasp	vs39,	vs5,	vs18		// a4_i*b1_r


.endm

.macro KERNEL2x1_SUBI1


	lxsspx		vs0,	o0,	AO		// load a0_r
	lxsspx		vs1,	o4,	AO		// load a0_i

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1		//  load b0_r
	lxsspx		vs9,	o4,	T1		//  load b0_i

	addi		T1,	T1,8

	lxsspx		vs10,	o0,	T1		//  load b1_r
	lxsspx		vs11,	o4,	T1		//  load b1_i

	addi		BO,	BO,	16


	xsmulsp		vs32,	vs0,	vs8		// a0_r*b0_r
	xsmulsp		vs33,	vs1,	vs9		// a0_i*b0_i
	xsmulsp		vs34,	vs0,	vs9		// a0_r*b0_i
	xsmulsp		vs35,	vs1,	vs8		// a0_i*b0_r

	xsmulsp		vs36,	vs0,	vs10		// a0_r*b1_r
	xsmulsp		vs37,	vs1,	vs11		// a0_i*b1_i
	xsmulsp		vs38,	vs0,	vs11		// a0_r*b1_i
	xsmulsp		vs39,	vs1,	vs10		// a0_i*b1_r


.endm

.macro KERNEL2x1_SUB1


	lxsspx		vs0,	o0,	AO		// load a0_r
	lxsspx		vs1,	o4,	AO		// load a0_i

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1		//  load b0_r
	lxsspx		vs9,	o4,	T1		//  load b0_i

	addi		T1,	T1,8

	lxsspx		vs10,	o0,	T1		//  load b1_r
	lxsspx		vs11,	o4,	T1		//  load b1_i

	addi		BO,	BO,	16


	xsmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r
	xsmaddasp	vs33,	vs1,	vs9		// a0_i*b0_i
	xsmaddasp	vs34,	vs0,	vs9		// a0_r*b0_i
	xsmaddasp	vs35,	vs1,	vs8		// a0_i*b0_r

	xsmaddasp	vs36,	vs0,	vs10		// a0_r*b1_r
	xsmaddasp	vs37,	vs1,	vs11		// a0_i*b1_i
	xsmaddasp	vs38,	vs0,	vs11		// a0_r*b1_i
	xsmaddasp	vs39,	vs1,	vs10		// a0_i*b1_r


.endm

.macro SAVE2x1

	mr		T1,	CO

// N=0

	mr		T2,	T1

// N=0 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5

#ifndef TRMMKERNEL
	lxsspx		vs0,	o0,	T2	// load c0_r
	lxsspx		vs1,	o4,	T2	// load c0_i
#else
	xxlxor		vs0,	vs0,	vs0
	xxlxor		vs1,	vs1,	vs1
#endif

	XSFADD_R1	vs4,	vs4,	vs32		// add a0_r * b0_r
	XSFADD_I1	vs5,	vs5,	vs35		// add a0_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs33		// add a0_i * b0_i
	XSFADD_I2	vs5,	vs5,	vs34		// add a0_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsaddsp		vs0,	vs0,	vs20
	xsaddsp		vs1,	vs1,	vs21


	stxsspx		vs0,	o0,	T2	// store c0_r
	stxsspx		vs1,	o4,	T2	// store c0_i

	addi		T2,	T2,	8
	add		T1,	T1,	LDC


// N=1

	mr		T2,	T1

// N=1 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5

#ifndef TRMMKERNEL
	lxsspx		vs0,	o0,	T2	// load c0_r
	lxsspx		vs1,	o4,	T2	// load c0_i
#else
	xxlxor		vs0,	vs0,	vs0
	xxlxor		vs1,	vs1,	vs1
#endif

	XSFADD_R1	vs4,	vs4,	vs36		// add a0_r * b0_r
	XSFADD_I1	vs5,	vs5,	vs39		// add a0_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs37		// add a0_i * b0_i
	XSFADD_I2	vs5,	vs5,	vs38		// add a0_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsaddsp		vs0,	vs0,	vs20
	xsaddsp		vs1,	vs1,	vs21


	stxsspx		vs0,	o0,	T2	// store c0_r
	stxsspx		vs1,	o4,	T2	// store c0_i

	addi		T2,	T2,	8
	add		T1,	T1,	LDC

	addi		CO,	CO,	8

.endm


/**********************************************************************************************
* Macros for N=1 and M=8
**********************************************************************************************/

.macro LOAD1x8_1

	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3

	lxvw4x		vs2,	o32,	AO		// load a4, a5

	lxvw4x		vs3,	o48,	AO		// load a6, a7


	addi		AO,	AO,	64

	lxvw4x		vs24,	o0,	BO		//  load b0, b1



	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	8

.endm

.macro KERNEL1x8_I1


	lxvw4x		vs4,	o0,	AO		// load a0, a1

	lxvw4x		vs5,	o16,	AO		// load a2, a3

	lxvw4x		vs6,	o32,	AO		// load a4, a5

	lxvw4x		vs7,	o48,	AO		// load a6, a7


	addi		AO,	AO,	64

	lxvw4x		vs24,	o0,	BO		//  load b0, b1



	xxspltw		vs16,	vs24,	0
	xxspltw		vs17,	vs24,	1
	xxspltw		vs18,	vs24,	2
	xxspltw		vs19,	vs24,	3


	addi		BO,	BO,	8


	xvmulsp		vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs36,	vs2,	vs8		// a2_r*b0_r, a2_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs37,	vs2,	vs9		// a2_r*b0_i, a2_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs38,	vs3,	vs8		// a3_r*b0_r, a3_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs39,	vs3,	vs9		// a3_r*b0_i, a3_i*b0_i, a1_r*b0_i, a1_i*b0_i


.endm

.macro KERNEL1x8_1


	lxvw4x		vs4,	o0,	AO		// load a0, a1

	lxvw4x		vs5,	o16,	AO		// load a2, a3

	lxvw4x		vs6,	o32,	AO		// load a4, a5

	lxvw4x		vs7,	o48,	AO		// load a6, a7


	addi		AO,	AO,	64

	lxvw4x		vs24,	o0,	BO		//  load b0, b1



	xxspltw		vs16,	vs24,	0
	xxspltw		vs17,	vs24,	1
	xxspltw		vs18,	vs24,	2
	xxspltw		vs19,	vs24,	3


	addi		BO,	BO,	8


	xvmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs36,	vs2,	vs8		// a2_r*b0_r, a2_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs37,	vs2,	vs9		// a2_r*b0_i, a2_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs38,	vs3,	vs8		// a3_r*b0_r, a3_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs39,	vs3,	vs9		// a3_r*b0_i, a3_i*b0_i, a1_r*b0_i, a1_i*b0_i


.endm

.macro KERNEL1x8_2


	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3

	lxvw4x		vs2,	o32,	AO		// load a4, a5

	lxvw4x		vs3,	o48,	AO		// load a6, a7


	addi		AO,	AO,	64

	lxvw4x		vs24,	o0,	BO		//  load b0, b1



	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	8


	xvmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r, a4_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs4,	vs17		// a4_r*b0_i, a4_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs5,	vs16		// a5_r*b0_r, a5_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs5,	vs17		// a5_r*b0_i, a5_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs36,	vs6,	vs16		// a6_r*b0_r, a6_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs37,	vs6,	vs17		// a6_r*b0_i, a6_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs38,	vs7,	vs16		// a7_r*b0_r, a7_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs39,	vs7,	vs17		// a7_r*b0_i, a7_i*b0_i, a1_r*b0_i, a1_i*b0_i


.endm

.macro KERNEL1x8_E2


	xvmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r, a4_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs4,	vs17		// a4_r*b0_i, a4_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs5,	vs16		// a5_r*b0_r, a5_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs5,	vs17		// a5_r*b0_i, a5_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs36,	vs6,	vs16		// a6_r*b0_r, a6_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs37,	vs6,	vs17		// a6_r*b0_i, a6_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs38,	vs7,	vs16		// a7_r*b0_r, a7_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs39,	vs7,	vs17		// a7_r*b0_i, a7_i*b0_i, a1_r*b0_i, a1_i*b0_i


.endm

.macro KERNEL1x8_SUBI1


	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3

	lxvw4x		vs2,	o32,	AO		// load a4, a5

	lxvw4x		vs3,	o48,	AO		// load a6, a7


	addi		AO,	AO,	64

	lxvw4x		vs24,	o0,	BO		//  load b0, b1



	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	8


	xvmulsp		vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs36,	vs2,	vs8		// a2_r*b0_r, a2_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs37,	vs2,	vs9		// a2_r*b0_i, a2_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs38,	vs3,	vs8		// a3_r*b0_r, a3_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs39,	vs3,	vs9		// a3_r*b0_i, a3_i*b0_i, a1_r*b0_i, a1_i*b0_i


.endm

.macro KERNEL1x8_SUB1


	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3

	lxvw4x		vs2,	o32,	AO		// load a4, a5

	lxvw4x		vs3,	o48,	AO		// load a6, a7


	addi		AO,	AO,	64

	lxvw4x		vs24,	o0,	BO		//  load b0, b1



	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	8


	xvmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs36,	vs2,	vs8		// a2_r*b0_r, a2_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs37,	vs2,	vs9		// a2_r*b0_i, a2_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs38,	vs3,	vs8		// a3_r*b0_r, a3_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs39,	vs3,	vs9		// a3_r*b0_i, a3_i*b0_i, a1_r*b0_i, a1_i*b0_i


.endm

.macro SAVE1x8

	mr		T1,	CO

// N=0

	mr		T2,	T1

// N=0 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs33,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=0 M=2

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs34,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs35,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=0 M=4

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs36,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs37,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=0 M=6

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs38,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs39,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC

	addi		CO,	CO,	64

.endm


/**********************************************************************************************
* Macros for N=1 and M=4
**********************************************************************************************/

.macro LOAD1x4_1

	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3


	addi		AO,	AO,	32

	lxvw4x		vs24,	o0,	BO		//  load b0, b1



	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	8

.endm

.macro KERNEL1x4_I1


	lxvw4x		vs4,	o0,	AO		// load a0, a1

	lxvw4x		vs5,	o16,	AO		// load a2, a3


	addi		AO,	AO,	32

	lxvw4x		vs24,	o0,	BO		//  load b0, b1



	xxspltw		vs16,	vs24,	0
	xxspltw		vs17,	vs24,	1
	xxspltw		vs18,	vs24,	2
	xxspltw		vs19,	vs24,	3


	addi		BO,	BO,	8


	xvmulsp		vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i


.endm

.macro KERNEL1x4_1


	lxvw4x		vs4,	o0,	AO		// load a0, a1

	lxvw4x		vs5,	o16,	AO		// load a2, a3


	addi		AO,	AO,	32

	lxvw4x		vs24,	o0,	BO		//  load b0, b1



	xxspltw		vs16,	vs24,	0
	xxspltw		vs17,	vs24,	1
	xxspltw		vs18,	vs24,	2
	xxspltw		vs19,	vs24,	3


	addi		BO,	BO,	8


	xvmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i


.endm

.macro KERNEL1x4_2


	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3


	addi		AO,	AO,	32

	lxvw4x		vs24,	o0,	BO		//  load b0, b1



	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	8


	xvmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r, a4_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs4,	vs17		// a4_r*b0_i, a4_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs5,	vs16		// a5_r*b0_r, a5_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs5,	vs17		// a5_r*b0_i, a5_i*b0_i, a1_r*b0_i, a1_i*b0_i


.endm

.macro KERNEL1x4_E2


	xvmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r, a4_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs4,	vs17		// a4_r*b0_i, a4_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs5,	vs16		// a5_r*b0_r, a5_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs5,	vs17		// a5_r*b0_i, a5_i*b0_i, a1_r*b0_i, a1_i*b0_i


.endm

.macro KERNEL1x4_SUBI1


	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3


	addi		AO,	AO,	32

	lxvw4x		vs24,	o0,	BO		//  load b0, b1



	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	8


	xvmulsp		vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmulsp		vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i


.endm

.macro KERNEL1x4_SUB1


	lxvw4x		vs0,	o0,	AO		// load a0, a1

	lxvw4x		vs1,	o16,	AO		// load a2, a3


	addi		AO,	AO,	32

	lxvw4x		vs24,	o0,	BO		//  load b0, b1



	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	8


	xvmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i
	xvmaddasp	vs34,	vs1,	vs8		// a1_r*b0_r, a1_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs35,	vs1,	vs9		// a1_r*b0_i, a1_i*b0_i, a1_r*b0_i, a1_i*b0_i


.endm

.macro SAVE1x4

	mr		T1,	CO

// N=0

	mr		T2,	T1

// N=0 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs33,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16

// N=0 M=2

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs34,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs35,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC

	addi		CO,	CO,	32

.endm


/**********************************************************************************************
* Macros for N=1 and M=2
**********************************************************************************************/

.macro LOAD1x2_1

	lxvw4x		vs0,	o0,	AO		// load a0, a1


	addi		AO,	AO,	16

	lxvw4x		vs24,	o0,	BO		//  load b0, b1



	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	8

.endm

.macro KERNEL1x2_I1


	lxvw4x		vs4,	o0,	AO		// load a0, a1


	addi		AO,	AO,	16

	lxvw4x		vs24,	o0,	BO		//  load b0, b1



	xxspltw		vs16,	vs24,	0
	xxspltw		vs17,	vs24,	1
	xxspltw		vs18,	vs24,	2
	xxspltw		vs19,	vs24,	3


	addi		BO,	BO,	8


	xvmulsp		vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i


.endm

.macro KERNEL1x2_1


	lxvw4x		vs4,	o0,	AO		// load a0, a1


	addi		AO,	AO,	16

	lxvw4x		vs24,	o0,	BO		//  load b0, b1



	xxspltw		vs16,	vs24,	0
	xxspltw		vs17,	vs24,	1
	xxspltw		vs18,	vs24,	2
	xxspltw		vs19,	vs24,	3


	addi		BO,	BO,	8


	xvmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i


.endm

.macro KERNEL1x2_2


	lxvw4x		vs0,	o0,	AO		// load a0, a1


	addi		AO,	AO,	16

	lxvw4x		vs24,	o0,	BO		//  load b0, b1



	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	8


	xvmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r, a4_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs4,	vs17		// a4_r*b0_i, a4_i*b0_i, a1_r*b0_i, a1_i*b0_i


.endm

.macro KERNEL1x2_E2


	xvmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r, a4_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs4,	vs17		// a4_r*b0_i, a4_i*b0_i, a1_r*b0_i, a1_i*b0_i


.endm

.macro KERNEL1x2_SUBI1


	lxvw4x		vs0,	o0,	AO		// load a0, a1


	addi		AO,	AO,	16

	lxvw4x		vs24,	o0,	BO		//  load b0, b1



	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	8


	xvmulsp		vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmulsp		vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i


.endm

.macro KERNEL1x2_SUB1


	lxvw4x		vs0,	o0,	AO		// load a0, a1


	addi		AO,	AO,	16

	lxvw4x		vs24,	o0,	BO		//  load b0, b1



	xxspltw		vs8,	vs24,	0
	xxspltw		vs9,	vs24,	1
	xxspltw		vs10,	vs24,	2
	xxspltw		vs11,	vs24,	3


	addi		BO,	BO,	8


	xvmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r, a0_i*b0_r, a1_r*b0_r, a1_i*b0_r
	xvmaddasp	vs33,	vs0,	vs9		// a0_r*b0_i, a0_i*b0_i, a1_r*b0_i, a1_i*b0_i


.endm

.macro SAVE1x2

	mr		T1,	CO

// N=0

	mr		T2,	T1

// N=0 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5
	xxlxor		vs6,	vs6,	vs6
	xxlxor		vs7,	vs7,	vs7

#ifndef TRMMKERNEL
	lxvw4x		vs0,	o0,	T2	// c0, c1
#else
	xxlxor		vs0,	vs0,	vs0
#endif

	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs8,	o0,	TBUFFER
	lxsspx		vs9,	o4,	TBUFFER
	lxsspx		vs10,	o8,	TBUFFER
	lxsspx		vs11,	o12,	TBUFFER

	stxvw4x		vs33,	o0,	TBUFFER

	lxsspx		vs12,	o0,	TBUFFER
	lxsspx		vs13,	o4,	TBUFFER
	lxsspx		vs14,	o8,	TBUFFER
	lxsspx		vs15,	o12,	TBUFFER

	XSFADD_R1	vs4,	vs4,	vs8		// add a0_r * b0_r
	XSFADD_I2	vs5,	vs5,	vs12		// add a0_r * b0_i
	XSFADD_R1	vs6,	vs6,	vs10		// add a1_r * b0_r
	XSFADD_I2	vs7,	vs7,	vs14		// add a1_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs13		// add a0_i * b0_i
	XSFADD_I1	vs5,	vs5,	vs9 		// add a0_i * b0_r
	XSFADD_R2	vs6,	vs6,	vs15		// add a1_i * b0_i
	XSFADD_I1	vs7,	vs7,	vs11		// add a1_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsmulsp		vs16,	vs6,	alpha_r		// r1_r * alpha_r
	xsmulsp		vs17,	vs7,	alpha_i		// r1_i * alpha_i
	xsmulsp		vs18,	vs6,	alpha_i		// r1_r * alpha_i
	xsmulsp		vs19,	vs7,	alpha_r		// r1_i * alpha_r

	xssubsp		vs22,	vs16,	vs17		// r1_r * alpha_r - r1_i * alpha_i
	xsaddsp		vs23,	vs18,	vs19		// r1_r * alpha_i + r1_i * alpha_r

	stxsspx		vs20,	o0,	TBUFFER		// store r0_r
	stxsspx		vs21,	o4,	TBUFFER		// store r0_i
	stxsspx		vs22,	o8,	TBUFFER		// store r1_r
	stxsspx		vs23,	o12,	TBUFFER		// store r1_i
	lxvw4x		vs1,	o0,	TBUFFER			// load r0_r, r0_i, r1_r, r1_i
	xvaddsp		vs0,	vs0,	vs1


	stxvw4x		vs0,	o0,	T2	// c0, c1

	addi		T2,	T2,	16
	add		T1,	T1,	LDC

	addi		CO,	CO,	16

.endm


/**********************************************************************************************
* Macros for N=1 and M=1
**********************************************************************************************/

.macro LOAD1x1_1

	lxsspx		vs0,	o0,	AO		// load a0_r
	lxsspx		vs1,	o4,	AO		// load a0_i

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1		//  load b0_r
	lxsspx		vs9,	o4,	T1		//  load b0_i

	addi		BO,	BO,	8

.endm

.macro KERNEL1x1_I1


	lxsspx		vs4,	o0,	AO		// load a0_r
	lxsspx		vs5,	o4,	AO		// load a0_i

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1		//  load b0_r
	lxsspx		vs17,	o4,	T1		//  load b0_i

	addi		BO,	BO,	8


	xsmulsp		vs32,	vs0,	vs8		// a0_r*b0_r
	xsmulsp		vs33,	vs1,	vs9		// a0_i*b0_i
	xsmulsp		vs34,	vs0,	vs9		// a0_r*b0_i
	xsmulsp		vs35,	vs1,	vs8		// a0_i*b0_r


.endm

.macro KERNEL1x1_1


	lxsspx		vs4,	o0,	AO		// load a0_r
	lxsspx		vs5,	o4,	AO		// load a0_i

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1		//  load b0_r
	lxsspx		vs17,	o4,	T1		//  load b0_i

	addi		BO,	BO,	8


	xsmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r
	xsmaddasp	vs33,	vs1,	vs9		// a0_i*b0_i
	xsmaddasp	vs34,	vs0,	vs9		// a0_r*b0_i
	xsmaddasp	vs35,	vs1,	vs8		// a0_i*b0_r


.endm

.macro KERNEL1x1_2


	lxsspx		vs0,	o0,	AO		// load a0_r
	lxsspx		vs1,	o4,	AO		// load a0_i

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1		//  load b0_r
	lxsspx		vs9,	o4,	T1		//  load b0_i

	addi		BO,	BO,	8


	xsmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r
	xsmaddasp	vs33,	vs5,	vs17		// a4_i*b0_i
	xsmaddasp	vs34,	vs4,	vs17		// a4_r*b0_i
	xsmaddasp	vs35,	vs5,	vs16		// a4_i*b0_r


.endm

.macro KERNEL1x1_E2


	xsmaddasp	vs32,	vs4,	vs16		// a4_r*b0_r
	xsmaddasp	vs33,	vs5,	vs17		// a4_i*b0_i
	xsmaddasp	vs34,	vs4,	vs17		// a4_r*b0_i
	xsmaddasp	vs35,	vs5,	vs16		// a4_i*b0_r


.endm

.macro KERNEL1x1_SUBI1


	lxsspx		vs0,	o0,	AO		// load a0_r
	lxsspx		vs1,	o4,	AO		// load a0_i

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1		//  load b0_r
	lxsspx		vs9,	o4,	T1		//  load b0_i

	addi		BO,	BO,	8


	xsmulsp		vs32,	vs0,	vs8		// a0_r*b0_r
	xsmulsp		vs33,	vs1,	vs9		// a0_i*b0_i
	xsmulsp		vs34,	vs0,	vs9		// a0_r*b0_i
	xsmulsp		vs35,	vs1,	vs8		// a0_i*b0_r


.endm

.macro KERNEL1x1_SUB1


	lxsspx		vs0,	o0,	AO		// load a0_r
	lxsspx		vs1,	o4,	AO		// load a0_i

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1		//  load b0_r
	lxsspx		vs9,	o4,	T1		//  load b0_i

	addi		BO,	BO,	8


	xsmaddasp	vs32,	vs0,	vs8		// a0_r*b0_r
	xsmaddasp	vs33,	vs1,	vs9		// a0_i*b0_i
	xsmaddasp	vs34,	vs0,	vs9		// a0_r*b0_i
	xsmaddasp	vs35,	vs1,	vs8		// a0_i*b0_r


.endm

.macro SAVE1x1

	mr		T1,	CO

// N=0

	mr		T2,	T1

// N=0 M=0

	xxlxor		vs4,	vs4,	vs4
	xxlxor		vs5,	vs5,	vs5

#ifndef TRMMKERNEL
	lxsspx		vs0,	o0,	T2	// load c0_r
	lxsspx		vs1,	o4,	T2	// load c0_i
#else
	xxlxor		vs0,	vs0,	vs0
	xxlxor		vs1,	vs1,	vs1
#endif

	XSFADD_R1	vs4,	vs4,	vs32		// add a0_r * b0_r
	XSFADD_I1	vs5,	vs5,	vs35		// add a0_r * b0_i

	XSFADD_R2	vs4,	vs4,	vs33		// add a0_i * b0_i
	XSFADD_I2	vs5,	vs5,	vs34		// add a0_i * b0_r

	xsmulsp		vs16,	vs4,	alpha_r		// r0_r * alpha_r
	xsmulsp		vs17,	vs5,	alpha_i		// r0_i * alpha_i
	xsmulsp		vs18,	vs4,	alpha_i		// r0_r * alpha_i
	xsmulsp		vs19,	vs5,	alpha_r		// r0_i * alpha_r

	xssubsp		vs20,	vs16,	vs17		// r0_r * alpha_r - r0_i * alpha_i
	xsaddsp		vs21,	vs18,	vs19		// r0_r * alpha_i + r0_i * alpha_r

	xsaddsp		vs0,	vs0,	vs20
	xsaddsp		vs1,	vs1,	vs21


	stxsspx		vs0,	o0,	T2	// store c0_r
	stxsspx		vs1,	o4,	T2	// store c0_i

	addi		T2,	T2,	8
	add		T1,	T1,	LDC

	addi		CO,	CO,	8

.endm

