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

/**********************************************************************************************
* Macros for N=8 and M=16
**********************************************************************************************/

.macro LOAD8x16_1

	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO
	lxvw4x		vs2,	o32,	AO
	lxvw4x		vs3,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	lxvw4x		vs29,	o16,	BO

	xxspltw		vs12,	vs29,	0
	xxspltw		vs13,	vs29,	1
	xxspltw		vs14,	vs29,	2
	xxspltw		vs15,	vs29,	3

	addi		BO,	BO,	32

.endm

.macro KERNEL8x16_I1


	lxvw4x		vs4,	o0,	AO
	lxvw4x		vs5,	o16,	AO
	lxvw4x		vs6,	o32,	AO
	lxvw4x		vs7,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0
	xxspltw		vs17,	vs28,	1
	xxspltw		vs18,	vs28,	2
	xxspltw		vs19,	vs28,	3

	lxvw4x		vs29,	o16,	BO

	xxspltw		vs20,	vs29,	0
	xxspltw		vs21,	vs29,	1
	xxspltw		vs22,	vs29,	2
	xxspltw		vs23,	vs29,	3

	addi		BO,	BO,	32


	xvmulsp		vs32,	vs0,	vs8
	xvmulsp		vs33,	vs1,	vs8
	xvmulsp		vs34,	vs2,	vs8
	xvmulsp		vs35,	vs3,	vs8

	xvmulsp		vs36,	vs0,	vs9
	xvmulsp		vs37,	vs1,	vs9
	xvmulsp		vs38,	vs2,	vs9
	xvmulsp		vs39,	vs3,	vs9

	xvmulsp		vs40,	vs0,	vs10
	xvmulsp		vs41,	vs1,	vs10
	xvmulsp		vs42,	vs2,	vs10
	xvmulsp		vs43,	vs3,	vs10

	xvmulsp		vs44,	vs0,	vs11
	xvmulsp		vs45,	vs1,	vs11
	xvmulsp		vs46,	vs2,	vs11
	xvmulsp		vs47,	vs3,	vs11

	xvmulsp		vs48,	vs0,	vs12
	xvmulsp		vs49,	vs1,	vs12
	xvmulsp		vs50,	vs2,	vs12
	xvmulsp		vs51,	vs3,	vs12

	xvmulsp		vs52,	vs0,	vs13
	xvmulsp		vs53,	vs1,	vs13
	xvmulsp		vs54,	vs2,	vs13
	xvmulsp		vs55,	vs3,	vs13

	xvmulsp		vs56,	vs0,	vs14
	xvmulsp		vs57,	vs1,	vs14
	xvmulsp		vs58,	vs2,	vs14
	xvmulsp		vs59,	vs3,	vs14

	xvmulsp		vs60,	vs0,	vs15
	xvmulsp		vs61,	vs1,	vs15
	xvmulsp		vs62,	vs2,	vs15
	xvmulsp		vs63,	vs3,	vs15


.endm

.macro KERNEL8x16_1

	xvmaddasp	vs32,	vs0,	vs8
	xvmaddasp	vs33,	vs1,	vs8
	lxvw4x		vs28,	o0,	BO
	lxvw4x		vs4,	o0,	AO
	xvmaddasp	vs34,	vs2,	vs8
	xvmaddasp	vs35,	vs3,	vs8

	xvmaddasp	vs36,	vs0,	vs9
	xvmaddasp	vs37,	vs1,	vs9
	lxvw4x		vs29,	o16,	BO
	lxvw4x		vs5,	o16,	AO
	xvmaddasp	vs38,	vs2,	vs9
	xvmaddasp	vs39,	vs3,	vs9

	xvmaddasp	vs40,	vs0,	vs10
	xvmaddasp	vs41,	vs1,	vs10
	lxvw4x		vs6,	o32,	AO
	lxvw4x		vs7,	o48,	AO
	xvmaddasp	vs42,	vs2,	vs10
	xvmaddasp	vs43,	vs3,	vs10

	xxspltw		vs16,	vs28,	0
	xxspltw		vs17,	vs28,	1
	xxspltw		vs18,	vs28,	2
	xxspltw		vs19,	vs28,	3

	xvmaddasp	vs44,	vs0,	vs11
	xvmaddasp	vs45,	vs1,	vs11
	xvmaddasp	vs46,	vs2,	vs11
	xvmaddasp	vs47,	vs3,	vs11

	xxspltw		vs20,	vs29,	0
	xxspltw		vs21,	vs29,	1
	xxspltw		vs22,	vs29,	2
	xxspltw		vs23,	vs29,	3

	xvmaddasp	vs48,	vs0,	vs12
	xvmaddasp	vs49,	vs1,	vs12
	xvmaddasp	vs50,	vs2,	vs12
	xvmaddasp	vs51,	vs3,	vs12

	xvmaddasp	vs52,	vs0,	vs13
	xvmaddasp	vs53,	vs1,	vs13
	xvmaddasp	vs54,	vs2,	vs13
	xvmaddasp	vs55,	vs3,	vs13

	xvmaddasp	vs56,	vs0,	vs14
	xvmaddasp	vs57,	vs1,	vs14
	addi		AO,	AO,	64
	addi		BO,	BO,	32
	xvmaddasp	vs58,	vs2,	vs14
	xvmaddasp	vs59,	vs3,	vs14

	xvmaddasp	vs60,	vs0,	vs15
	xvmaddasp	vs61,	vs1,	vs15
	xvmaddasp	vs62,	vs2,	vs15
	xvmaddasp	vs63,	vs3,	vs15


.endm

.macro KERNEL8x16_2

	xvmaddasp	vs32,	vs4,	vs16
	xvmaddasp	vs33,	vs5,	vs16

	lxvw4x		vs28,	o0,	BO
	lxvw4x		vs0,	o0,	AO

	xvmaddasp	vs34,	vs6,	vs16
	xvmaddasp	vs35,	vs7,	vs16

	xvmaddasp	vs36,	vs4,	vs17
	xvmaddasp	vs37,	vs5,	vs17

	lxvw4x		vs29,	o16,	BO
	lxvw4x		vs1,	o16,	AO

	xvmaddasp	vs38,	vs6,	vs17
	xvmaddasp	vs39,	vs7,	vs17

	lxvw4x		vs2,	o32,	AO
	lxvw4x		vs3,	o48,	AO

	xvmaddasp	vs40,	vs4,	vs18
	xvmaddasp	vs41,	vs5,	vs18
	xvmaddasp	vs42,	vs6,	vs18
	xvmaddasp	vs43,	vs7,	vs18

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	xvmaddasp	vs44,	vs4,	vs19
	xvmaddasp	vs45,	vs5,	vs19
	xvmaddasp	vs46,	vs6,	vs19
	xvmaddasp	vs47,	vs7,	vs19

	xxspltw		vs12,	vs29,	0
	xxspltw		vs13,	vs29,	1
	xxspltw		vs14,	vs29,	2
	xxspltw		vs15,	vs29,	3

	xvmaddasp	vs48,	vs4,	vs20
	xvmaddasp	vs49,	vs5,	vs20
	xvmaddasp	vs50,	vs6,	vs20
	xvmaddasp	vs51,	vs7,	vs20

	xvmaddasp	vs52,	vs4,	vs21
	xvmaddasp	vs53,	vs5,	vs21
	xvmaddasp	vs54,	vs6,	vs21
	xvmaddasp	vs55,	vs7,	vs21

	xvmaddasp	vs56,	vs4,	vs22
	xvmaddasp	vs57,	vs5,	vs22
	xvmaddasp	vs58,	vs6,	vs22
	xvmaddasp	vs59,	vs7,	vs22

	xvmaddasp	vs60,	vs4,	vs23
	xvmaddasp	vs61,	vs5,	vs23
	addi		AO,	AO,	64
	addi		BO,	BO,	32
	xvmaddasp	vs62,	vs6,	vs23
	xvmaddasp	vs63,	vs7,	vs23


.endm

.macro KERNEL8x16_E2


	xvmaddasp	vs32,	vs4,	vs16
	xvmaddasp	vs33,	vs5,	vs16
	xvmaddasp	vs34,	vs6,	vs16
	xvmaddasp	vs35,	vs7,	vs16

	xvmaddasp	vs36,	vs4,	vs17
	xvmaddasp	vs37,	vs5,	vs17
	xvmaddasp	vs38,	vs6,	vs17
	xvmaddasp	vs39,	vs7,	vs17

	xvmaddasp	vs40,	vs4,	vs18
	xvmaddasp	vs41,	vs5,	vs18
	xvmaddasp	vs42,	vs6,	vs18
	xvmaddasp	vs43,	vs7,	vs18

	xvmaddasp	vs44,	vs4,	vs19
	xvmaddasp	vs45,	vs5,	vs19
	xvmaddasp	vs46,	vs6,	vs19
	xvmaddasp	vs47,	vs7,	vs19

	xvmaddasp	vs48,	vs4,	vs20
	xvmaddasp	vs49,	vs5,	vs20
	xvmaddasp	vs50,	vs6,	vs20
	xvmaddasp	vs51,	vs7,	vs20

	xvmaddasp	vs52,	vs4,	vs21
	xvmaddasp	vs53,	vs5,	vs21
	xvmaddasp	vs54,	vs6,	vs21
	xvmaddasp	vs55,	vs7,	vs21

	xvmaddasp	vs56,	vs4,	vs22
	xvmaddasp	vs57,	vs5,	vs22
	xvmaddasp	vs58,	vs6,	vs22
	xvmaddasp	vs59,	vs7,	vs22

	xvmaddasp	vs60,	vs4,	vs23
	xvmaddasp	vs61,	vs5,	vs23
	xvmaddasp	vs62,	vs6,	vs23
	xvmaddasp	vs63,	vs7,	vs23


.endm

.macro KERNEL8x16_SUBI1


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO
	lxvw4x		vs2,	o32,	AO
	lxvw4x		vs3,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	lxvw4x		vs29,	o16,	BO

	xxspltw		vs12,	vs29,	0
	xxspltw		vs13,	vs29,	1
	xxspltw		vs14,	vs29,	2
	xxspltw		vs15,	vs29,	3

	addi		BO,	BO,	32


	xvmulsp		vs32,	vs0,	vs8
	xvmulsp		vs33,	vs1,	vs8
	xvmulsp		vs34,	vs2,	vs8
	xvmulsp		vs35,	vs3,	vs8

	xvmulsp		vs36,	vs0,	vs9
	xvmulsp		vs37,	vs1,	vs9
	xvmulsp		vs38,	vs2,	vs9
	xvmulsp		vs39,	vs3,	vs9

	xvmulsp		vs40,	vs0,	vs10
	xvmulsp		vs41,	vs1,	vs10
	xvmulsp		vs42,	vs2,	vs10
	xvmulsp		vs43,	vs3,	vs10

	xvmulsp		vs44,	vs0,	vs11
	xvmulsp		vs45,	vs1,	vs11
	xvmulsp		vs46,	vs2,	vs11
	xvmulsp		vs47,	vs3,	vs11

	xvmulsp		vs48,	vs0,	vs12
	xvmulsp		vs49,	vs1,	vs12
	xvmulsp		vs50,	vs2,	vs12
	xvmulsp		vs51,	vs3,	vs12

	xvmulsp		vs52,	vs0,	vs13
	xvmulsp		vs53,	vs1,	vs13
	xvmulsp		vs54,	vs2,	vs13
	xvmulsp		vs55,	vs3,	vs13

	xvmulsp		vs56,	vs0,	vs14
	xvmulsp		vs57,	vs1,	vs14
	xvmulsp		vs58,	vs2,	vs14
	xvmulsp		vs59,	vs3,	vs14

	xvmulsp		vs60,	vs0,	vs15
	xvmulsp		vs61,	vs1,	vs15
	xvmulsp		vs62,	vs2,	vs15
	xvmulsp		vs63,	vs3,	vs15


.endm

.macro KERNEL8x16_SUB1


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO
	lxvw4x		vs2,	o32,	AO
	lxvw4x		vs3,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	lxvw4x		vs29,	o16,	BO

	xxspltw		vs12,	vs29,	0
	xxspltw		vs13,	vs29,	1
	xxspltw		vs14,	vs29,	2
	xxspltw		vs15,	vs29,	3

	addi		BO,	BO,	32


	xvmaddasp	vs32,	vs0,	vs8
	xvmaddasp	vs33,	vs1,	vs8
	xvmaddasp	vs34,	vs2,	vs8
	xvmaddasp	vs35,	vs3,	vs8

	xvmaddasp	vs36,	vs0,	vs9
	xvmaddasp	vs37,	vs1,	vs9
	xvmaddasp	vs38,	vs2,	vs9
	xvmaddasp	vs39,	vs3,	vs9

	xvmaddasp	vs40,	vs0,	vs10
	xvmaddasp	vs41,	vs1,	vs10
	xvmaddasp	vs42,	vs2,	vs10
	xvmaddasp	vs43,	vs3,	vs10

	xvmaddasp	vs44,	vs0,	vs11
	xvmaddasp	vs45,	vs1,	vs11
	xvmaddasp	vs46,	vs2,	vs11
	xvmaddasp	vs47,	vs3,	vs11

	xvmaddasp	vs48,	vs0,	vs12
	xvmaddasp	vs49,	vs1,	vs12
	xvmaddasp	vs50,	vs2,	vs12
	xvmaddasp	vs51,	vs3,	vs12

	xvmaddasp	vs52,	vs0,	vs13
	xvmaddasp	vs53,	vs1,	vs13
	xvmaddasp	vs54,	vs2,	vs13
	xvmaddasp	vs55,	vs3,	vs13

	xvmaddasp	vs56,	vs0,	vs14
	xvmaddasp	vs57,	vs1,	vs14
	xvmaddasp	vs58,	vs2,	vs14
	xvmaddasp	vs59,	vs3,	vs14

	xvmaddasp	vs60,	vs0,	vs15
	xvmaddasp	vs61,	vs1,	vs15
	xvmaddasp	vs62,	vs2,	vs15
	xvmaddasp	vs63,	vs3,	vs15


.endm

.macro SAVE8x16

	mr		T1,	CO

#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1
	lxvw4x		vs2,	o32,	T1
	lxvw4x		vs3,	o48,	T1

#endif


	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs33,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif

	stxvw4x		vs34,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs2,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs2,	vs2,	vs28
#endif

	stxvw4x		vs35,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs3,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs3,	vs3,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1
	stxvw4x		vs2,	o32,	T1
	stxvw4x		vs3,	o48,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1
	lxvw4x		vs2,	o32,	T1
	lxvw4x		vs3,	o48,	T1

#endif


	stxvw4x		vs36,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs37,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif

	stxvw4x		vs38,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs2,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs2,	vs2,	vs28
#endif

	stxvw4x		vs39,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs3,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs3,	vs3,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1
	stxvw4x		vs2,	o32,	T1
	stxvw4x		vs3,	o48,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1
	lxvw4x		vs2,	o32,	T1
	lxvw4x		vs3,	o48,	T1

#endif


	stxvw4x		vs40,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs41,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif

	stxvw4x		vs42,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs2,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs2,	vs2,	vs28
#endif

	stxvw4x		vs43,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs3,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs3,	vs3,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1
	stxvw4x		vs2,	o32,	T1
	stxvw4x		vs3,	o48,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1
	lxvw4x		vs2,	o32,	T1
	lxvw4x		vs3,	o48,	T1

#endif


	stxvw4x		vs44,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs45,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif

	stxvw4x		vs46,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs2,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs2,	vs2,	vs28
#endif

	stxvw4x		vs47,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs3,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs3,	vs3,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1
	stxvw4x		vs2,	o32,	T1
	stxvw4x		vs3,	o48,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1
	lxvw4x		vs2,	o32,	T1
	lxvw4x		vs3,	o48,	T1

#endif


	stxvw4x		vs48,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs49,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif

	stxvw4x		vs50,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs2,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs2,	vs2,	vs28
#endif

	stxvw4x		vs51,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs3,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs3,	vs3,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1
	stxvw4x		vs2,	o32,	T1
	stxvw4x		vs3,	o48,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1
	lxvw4x		vs2,	o32,	T1
	lxvw4x		vs3,	o48,	T1

#endif


	stxvw4x		vs52,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs53,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif

	stxvw4x		vs54,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs2,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs2,	vs2,	vs28
#endif

	stxvw4x		vs55,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs3,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs3,	vs3,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1
	stxvw4x		vs2,	o32,	T1
	stxvw4x		vs3,	o48,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1
	lxvw4x		vs2,	o32,	T1
	lxvw4x		vs3,	o48,	T1

#endif


	stxvw4x		vs56,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs57,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif

	stxvw4x		vs58,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs2,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs2,	vs2,	vs28
#endif

	stxvw4x		vs59,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs3,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs3,	vs3,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1
	stxvw4x		vs2,	o32,	T1
	stxvw4x		vs3,	o48,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1
	lxvw4x		vs2,	o32,	T1
	lxvw4x		vs3,	o48,	T1

#endif


	stxvw4x		vs60,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs61,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif

	stxvw4x		vs62,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs2,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs2,	vs2,	vs28
#endif

	stxvw4x		vs63,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs3,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs3,	vs3,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1
	stxvw4x		vs2,	o32,	T1
	stxvw4x		vs3,	o48,	T1

	add		T1,	T1,	LDC

	addi		CO,	CO,	64

.endm


/**********************************************************************************************
* Macros for N=8 and M=8
**********************************************************************************************/

.macro LOAD8x8_1

	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	lxvw4x		vs29,	o16,	BO

	xxspltw		vs12,	vs29,	0
	xxspltw		vs13,	vs29,	1
	xxspltw		vs14,	vs29,	2
	xxspltw		vs15,	vs29,	3

	addi		BO,	BO,	32

.endm

.macro KERNEL8x8_I1


	lxvw4x		vs4,	o0,	AO
	lxvw4x		vs5,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0
	xxspltw		vs17,	vs28,	1
	xxspltw		vs18,	vs28,	2
	xxspltw		vs19,	vs28,	3

	lxvw4x		vs29,	o16,	BO

	xxspltw		vs20,	vs29,	0
	xxspltw		vs21,	vs29,	1
	xxspltw		vs22,	vs29,	2
	xxspltw		vs23,	vs29,	3

	addi		BO,	BO,	32


	xvmulsp		vs32,	vs0,	vs8
	xvmulsp		vs33,	vs1,	vs8

	xvmulsp		vs34,	vs0,	vs9
	xvmulsp		vs35,	vs1,	vs9

	xvmulsp		vs36,	vs0,	vs10
	xvmulsp		vs37,	vs1,	vs10

	xvmulsp		vs38,	vs0,	vs11
	xvmulsp		vs39,	vs1,	vs11

	xvmulsp		vs40,	vs0,	vs12
	xvmulsp		vs41,	vs1,	vs12

	xvmulsp		vs42,	vs0,	vs13
	xvmulsp		vs43,	vs1,	vs13

	xvmulsp		vs44,	vs0,	vs14
	xvmulsp		vs45,	vs1,	vs14

	xvmulsp		vs46,	vs0,	vs15
	xvmulsp		vs47,	vs1,	vs15


.endm

.macro KERNEL8x8_1


	lxvw4x		vs4,	o0,	AO
	lxvw4x		vs5,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0
	xxspltw		vs17,	vs28,	1
	xxspltw		vs18,	vs28,	2
	xxspltw		vs19,	vs28,	3

	lxvw4x		vs29,	o16,	BO

	xxspltw		vs20,	vs29,	0
	xxspltw		vs21,	vs29,	1
	xxspltw		vs22,	vs29,	2
	xxspltw		vs23,	vs29,	3

	addi		BO,	BO,	32


	xvmaddasp	vs32,	vs0,	vs8
	xvmaddasp	vs33,	vs1,	vs8

	xvmaddasp	vs34,	vs0,	vs9
	xvmaddasp	vs35,	vs1,	vs9

	xvmaddasp	vs36,	vs0,	vs10
	xvmaddasp	vs37,	vs1,	vs10

	xvmaddasp	vs38,	vs0,	vs11
	xvmaddasp	vs39,	vs1,	vs11

	xvmaddasp	vs40,	vs0,	vs12
	xvmaddasp	vs41,	vs1,	vs12

	xvmaddasp	vs42,	vs0,	vs13
	xvmaddasp	vs43,	vs1,	vs13

	xvmaddasp	vs44,	vs0,	vs14
	xvmaddasp	vs45,	vs1,	vs14

	xvmaddasp	vs46,	vs0,	vs15
	xvmaddasp	vs47,	vs1,	vs15


.endm

.macro KERNEL8x8_2


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	lxvw4x		vs29,	o16,	BO

	xxspltw		vs12,	vs29,	0
	xxspltw		vs13,	vs29,	1
	xxspltw		vs14,	vs29,	2
	xxspltw		vs15,	vs29,	3

	addi		BO,	BO,	32


	xvmaddasp	vs32,	vs4,	vs16
	xvmaddasp	vs33,	vs5,	vs16

	xvmaddasp	vs34,	vs4,	vs17
	xvmaddasp	vs35,	vs5,	vs17

	xvmaddasp	vs36,	vs4,	vs18
	xvmaddasp	vs37,	vs5,	vs18

	xvmaddasp	vs38,	vs4,	vs19
	xvmaddasp	vs39,	vs5,	vs19

	xvmaddasp	vs40,	vs4,	vs20
	xvmaddasp	vs41,	vs5,	vs20

	xvmaddasp	vs42,	vs4,	vs21
	xvmaddasp	vs43,	vs5,	vs21

	xvmaddasp	vs44,	vs4,	vs22
	xvmaddasp	vs45,	vs5,	vs22

	xvmaddasp	vs46,	vs4,	vs23
	xvmaddasp	vs47,	vs5,	vs23


.endm

.macro KERNEL8x8_E2


	xvmaddasp	vs32,	vs4,	vs16
	xvmaddasp	vs33,	vs5,	vs16

	xvmaddasp	vs34,	vs4,	vs17
	xvmaddasp	vs35,	vs5,	vs17

	xvmaddasp	vs36,	vs4,	vs18
	xvmaddasp	vs37,	vs5,	vs18

	xvmaddasp	vs38,	vs4,	vs19
	xvmaddasp	vs39,	vs5,	vs19

	xvmaddasp	vs40,	vs4,	vs20
	xvmaddasp	vs41,	vs5,	vs20

	xvmaddasp	vs42,	vs4,	vs21
	xvmaddasp	vs43,	vs5,	vs21

	xvmaddasp	vs44,	vs4,	vs22
	xvmaddasp	vs45,	vs5,	vs22

	xvmaddasp	vs46,	vs4,	vs23
	xvmaddasp	vs47,	vs5,	vs23


.endm

.macro KERNEL8x8_SUBI1


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	lxvw4x		vs29,	o16,	BO

	xxspltw		vs12,	vs29,	0
	xxspltw		vs13,	vs29,	1
	xxspltw		vs14,	vs29,	2
	xxspltw		vs15,	vs29,	3

	addi		BO,	BO,	32


	xvmulsp		vs32,	vs0,	vs8
	xvmulsp		vs33,	vs1,	vs8

	xvmulsp		vs34,	vs0,	vs9
	xvmulsp		vs35,	vs1,	vs9

	xvmulsp		vs36,	vs0,	vs10
	xvmulsp		vs37,	vs1,	vs10

	xvmulsp		vs38,	vs0,	vs11
	xvmulsp		vs39,	vs1,	vs11

	xvmulsp		vs40,	vs0,	vs12
	xvmulsp		vs41,	vs1,	vs12

	xvmulsp		vs42,	vs0,	vs13
	xvmulsp		vs43,	vs1,	vs13

	xvmulsp		vs44,	vs0,	vs14
	xvmulsp		vs45,	vs1,	vs14

	xvmulsp		vs46,	vs0,	vs15
	xvmulsp		vs47,	vs1,	vs15


.endm

.macro KERNEL8x8_SUB1


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	lxvw4x		vs29,	o16,	BO

	xxspltw		vs12,	vs29,	0
	xxspltw		vs13,	vs29,	1
	xxspltw		vs14,	vs29,	2
	xxspltw		vs15,	vs29,	3

	addi		BO,	BO,	32


	xvmaddasp	vs32,	vs0,	vs8
	xvmaddasp	vs33,	vs1,	vs8

	xvmaddasp	vs34,	vs0,	vs9
	xvmaddasp	vs35,	vs1,	vs9

	xvmaddasp	vs36,	vs0,	vs10
	xvmaddasp	vs37,	vs1,	vs10

	xvmaddasp	vs38,	vs0,	vs11
	xvmaddasp	vs39,	vs1,	vs11

	xvmaddasp	vs40,	vs0,	vs12
	xvmaddasp	vs41,	vs1,	vs12

	xvmaddasp	vs42,	vs0,	vs13
	xvmaddasp	vs43,	vs1,	vs13

	xvmaddasp	vs44,	vs0,	vs14
	xvmaddasp	vs45,	vs1,	vs14

	xvmaddasp	vs46,	vs0,	vs15
	xvmaddasp	vs47,	vs1,	vs15


.endm

.macro SAVE8x8

	mr		T1,	CO

#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1

#endif


	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs33,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1

#endif


	stxvw4x		vs34,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs35,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1

#endif


	stxvw4x		vs36,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs37,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1

#endif


	stxvw4x		vs38,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs39,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1

#endif


	stxvw4x		vs40,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs41,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1

#endif


	stxvw4x		vs42,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs43,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1

#endif


	stxvw4x		vs44,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs45,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1

#endif


	stxvw4x		vs46,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs47,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1

	add		T1,	T1,	LDC

	addi		CO,	CO,	32

.endm


/**********************************************************************************************
* Macros for N=8 and M=4
**********************************************************************************************/

.macro LOAD8x4_1

	lxvw4x		vs0,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	lxvw4x		vs29,	o16,	BO

	xxspltw		vs12,	vs29,	0
	xxspltw		vs13,	vs29,	1
	xxspltw		vs14,	vs29,	2
	xxspltw		vs15,	vs29,	3

	addi		BO,	BO,	32

.endm

.macro KERNEL8x4_I1


	lxvw4x		vs4,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0
	xxspltw		vs17,	vs28,	1
	xxspltw		vs18,	vs28,	2
	xxspltw		vs19,	vs28,	3

	lxvw4x		vs29,	o16,	BO

	xxspltw		vs20,	vs29,	0
	xxspltw		vs21,	vs29,	1
	xxspltw		vs22,	vs29,	2
	xxspltw		vs23,	vs29,	3

	addi		BO,	BO,	32


	xvmulsp		vs32,	vs0,	vs8

	xvmulsp		vs33,	vs0,	vs9

	xvmulsp		vs34,	vs0,	vs10

	xvmulsp		vs35,	vs0,	vs11

	xvmulsp		vs36,	vs0,	vs12

	xvmulsp		vs37,	vs0,	vs13

	xvmulsp		vs38,	vs0,	vs14

	xvmulsp		vs39,	vs0,	vs15


.endm

.macro KERNEL8x4_1


	lxvw4x		vs4,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0
	xxspltw		vs17,	vs28,	1
	xxspltw		vs18,	vs28,	2
	xxspltw		vs19,	vs28,	3

	lxvw4x		vs29,	o16,	BO

	xxspltw		vs20,	vs29,	0
	xxspltw		vs21,	vs29,	1
	xxspltw		vs22,	vs29,	2
	xxspltw		vs23,	vs29,	3

	addi		BO,	BO,	32


	xvmaddasp	vs32,	vs0,	vs8

	xvmaddasp	vs33,	vs0,	vs9

	xvmaddasp	vs34,	vs0,	vs10

	xvmaddasp	vs35,	vs0,	vs11

	xvmaddasp	vs36,	vs0,	vs12

	xvmaddasp	vs37,	vs0,	vs13

	xvmaddasp	vs38,	vs0,	vs14

	xvmaddasp	vs39,	vs0,	vs15


.endm

.macro KERNEL8x4_2


	lxvw4x		vs0,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	lxvw4x		vs29,	o16,	BO

	xxspltw		vs12,	vs29,	0
	xxspltw		vs13,	vs29,	1
	xxspltw		vs14,	vs29,	2
	xxspltw		vs15,	vs29,	3

	addi		BO,	BO,	32


	xvmaddasp	vs32,	vs4,	vs16

	xvmaddasp	vs33,	vs4,	vs17

	xvmaddasp	vs34,	vs4,	vs18

	xvmaddasp	vs35,	vs4,	vs19

	xvmaddasp	vs36,	vs4,	vs20

	xvmaddasp	vs37,	vs4,	vs21

	xvmaddasp	vs38,	vs4,	vs22

	xvmaddasp	vs39,	vs4,	vs23


.endm

.macro KERNEL8x4_E2


	xvmaddasp	vs32,	vs4,	vs16

	xvmaddasp	vs33,	vs4,	vs17

	xvmaddasp	vs34,	vs4,	vs18

	xvmaddasp	vs35,	vs4,	vs19

	xvmaddasp	vs36,	vs4,	vs20

	xvmaddasp	vs37,	vs4,	vs21

	xvmaddasp	vs38,	vs4,	vs22

	xvmaddasp	vs39,	vs4,	vs23


.endm

.macro KERNEL8x4_SUBI1


	lxvw4x		vs0,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	lxvw4x		vs29,	o16,	BO

	xxspltw		vs12,	vs29,	0
	xxspltw		vs13,	vs29,	1
	xxspltw		vs14,	vs29,	2
	xxspltw		vs15,	vs29,	3

	addi		BO,	BO,	32


	xvmulsp		vs32,	vs0,	vs8

	xvmulsp		vs33,	vs0,	vs9

	xvmulsp		vs34,	vs0,	vs10

	xvmulsp		vs35,	vs0,	vs11

	xvmulsp		vs36,	vs0,	vs12

	xvmulsp		vs37,	vs0,	vs13

	xvmulsp		vs38,	vs0,	vs14

	xvmulsp		vs39,	vs0,	vs15


.endm

.macro KERNEL8x4_SUB1


	lxvw4x		vs0,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	lxvw4x		vs29,	o16,	BO

	xxspltw		vs12,	vs29,	0
	xxspltw		vs13,	vs29,	1
	xxspltw		vs14,	vs29,	2
	xxspltw		vs15,	vs29,	3

	addi		BO,	BO,	32


	xvmaddasp	vs32,	vs0,	vs8

	xvmaddasp	vs33,	vs0,	vs9

	xvmaddasp	vs34,	vs0,	vs10

	xvmaddasp	vs35,	vs0,	vs11

	xvmaddasp	vs36,	vs0,	vs12

	xvmaddasp	vs37,	vs0,	vs13

	xvmaddasp	vs38,	vs0,	vs14

	xvmaddasp	vs39,	vs0,	vs15


.endm

.macro SAVE8x4

	mr		T1,	CO

#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1

#endif


	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif




	stxvw4x		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1

#endif


	stxvw4x		vs33,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif




	stxvw4x		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1

#endif


	stxvw4x		vs34,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif




	stxvw4x		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1

#endif


	stxvw4x		vs35,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif




	stxvw4x		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1

#endif


	stxvw4x		vs36,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif




	stxvw4x		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1

#endif


	stxvw4x		vs37,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif




	stxvw4x		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1

#endif


	stxvw4x		vs38,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif




	stxvw4x		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1

#endif


	stxvw4x		vs39,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif




	stxvw4x		vs0,	o0,	T1

	add		T1,	T1,	LDC

	addi		CO,	CO,	16

.endm


/**********************************************************************************************
* Macros for N=8 and M=2
**********************************************************************************************/

.macro LOAD8x2_1

	lxsspx		vs0,	o0,	AO
	lxsspx		vs1,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1
	lxsspx		vs10,	o8,	T1
	lxsspx		vs11,	o12,	T1

	addi		T1,	T1,	16

	lxsspx		vs12,	o0,	T1
	lxsspx		vs13,	o4,	T1
	lxsspx		vs14,	o8,	T1
	lxsspx		vs15,	o12,	T1

	addi		BO,	BO,	32

.endm

.macro KERNEL8x2_I1


	lxsspx		vs4,	o0,	AO
	lxsspx		vs5,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1
	lxsspx		vs17,	o4,	T1
	lxsspx		vs18,	o8,	T1
	lxsspx		vs19,	o12,	T1

	addi		T1,	T1,	16

	lxsspx		vs20,	o0,	T1
	lxsspx		vs21,	o4,	T1
	lxsspx		vs22,	o8,	T1
	lxsspx		vs23,	o12,	T1

	addi		BO,	BO,	32


	xsmulsp		vs32,	vs0,	vs8
	xsmulsp		vs33,	vs1,	vs8

	xsmulsp		vs34,	vs0,	vs9
	xsmulsp		vs35,	vs1,	vs9

	xsmulsp		vs36,	vs0,	vs10
	xsmulsp		vs37,	vs1,	vs10

	xsmulsp		vs38,	vs0,	vs11
	xsmulsp		vs39,	vs1,	vs11

	xsmulsp		vs40,	vs0,	vs12
	xsmulsp		vs41,	vs1,	vs12

	xsmulsp		vs42,	vs0,	vs13
	xsmulsp		vs43,	vs1,	vs13

	xsmulsp		vs44,	vs0,	vs14
	xsmulsp		vs45,	vs1,	vs14

	xsmulsp		vs46,	vs0,	vs15
	xsmulsp		vs47,	vs1,	vs15


.endm

.macro KERNEL8x2_1


	lxsspx		vs4,	o0,	AO
	lxsspx		vs5,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1
	lxsspx		vs17,	o4,	T1
	lxsspx		vs18,	o8,	T1
	lxsspx		vs19,	o12,	T1

	addi		T1,	T1,	16

	lxsspx		vs20,	o0,	T1
	lxsspx		vs21,	o4,	T1
	lxsspx		vs22,	o8,	T1
	lxsspx		vs23,	o12,	T1

	addi		BO,	BO,	32


	xsmaddasp	vs32,	vs0,	vs8
	xsmaddasp	vs33,	vs1,	vs8

	xsmaddasp	vs34,	vs0,	vs9
	xsmaddasp	vs35,	vs1,	vs9

	xsmaddasp	vs36,	vs0,	vs10
	xsmaddasp	vs37,	vs1,	vs10

	xsmaddasp	vs38,	vs0,	vs11
	xsmaddasp	vs39,	vs1,	vs11

	xsmaddasp	vs40,	vs0,	vs12
	xsmaddasp	vs41,	vs1,	vs12

	xsmaddasp	vs42,	vs0,	vs13
	xsmaddasp	vs43,	vs1,	vs13

	xsmaddasp	vs44,	vs0,	vs14
	xsmaddasp	vs45,	vs1,	vs14

	xsmaddasp	vs46,	vs0,	vs15
	xsmaddasp	vs47,	vs1,	vs15


.endm

.macro KERNEL8x2_2


	lxsspx		vs0,	o0,	AO
	lxsspx		vs1,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1
	lxsspx		vs10,	o8,	T1
	lxsspx		vs11,	o12,	T1

	addi		T1,	T1,	16

	lxsspx		vs12,	o0,	T1
	lxsspx		vs13,	o4,	T1
	lxsspx		vs14,	o8,	T1
	lxsspx		vs15,	o12,	T1

	addi		BO,	BO,	32


	xsmaddasp	vs32,	vs4,	vs16
	xsmaddasp	vs33,	vs5,	vs16

	xsmaddasp	vs34,	vs4,	vs17
	xsmaddasp	vs35,	vs5,	vs17

	xsmaddasp	vs36,	vs4,	vs18
	xsmaddasp	vs37,	vs5,	vs18

	xsmaddasp	vs38,	vs4,	vs19
	xsmaddasp	vs39,	vs5,	vs19

	xsmaddasp	vs40,	vs4,	vs20
	xsmaddasp	vs41,	vs5,	vs20

	xsmaddasp	vs42,	vs4,	vs21
	xsmaddasp	vs43,	vs5,	vs21

	xsmaddasp	vs44,	vs4,	vs22
	xsmaddasp	vs45,	vs5,	vs22

	xsmaddasp	vs46,	vs4,	vs23
	xsmaddasp	vs47,	vs5,	vs23


.endm

.macro KERNEL8x2_E2


	xsmaddasp	vs32,	vs4,	vs16
	xsmaddasp	vs33,	vs5,	vs16

	xsmaddasp	vs34,	vs4,	vs17
	xsmaddasp	vs35,	vs5,	vs17

	xsmaddasp	vs36,	vs4,	vs18
	xsmaddasp	vs37,	vs5,	vs18

	xsmaddasp	vs38,	vs4,	vs19
	xsmaddasp	vs39,	vs5,	vs19

	xsmaddasp	vs40,	vs4,	vs20
	xsmaddasp	vs41,	vs5,	vs20

	xsmaddasp	vs42,	vs4,	vs21
	xsmaddasp	vs43,	vs5,	vs21

	xsmaddasp	vs44,	vs4,	vs22
	xsmaddasp	vs45,	vs5,	vs22

	xsmaddasp	vs46,	vs4,	vs23
	xsmaddasp	vs47,	vs5,	vs23


.endm

.macro KERNEL8x2_SUBI1


	lxsspx		vs0,	o0,	AO
	lxsspx		vs1,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1
	lxsspx		vs10,	o8,	T1
	lxsspx		vs11,	o12,	T1

	addi		T1,	T1,	16

	lxsspx		vs12,	o0,	T1
	lxsspx		vs13,	o4,	T1
	lxsspx		vs14,	o8,	T1
	lxsspx		vs15,	o12,	T1

	addi		BO,	BO,	32


	xsmulsp		vs32,	vs0,	vs8
	xsmulsp		vs33,	vs1,	vs8

	xsmulsp		vs34,	vs0,	vs9
	xsmulsp		vs35,	vs1,	vs9

	xsmulsp		vs36,	vs0,	vs10
	xsmulsp		vs37,	vs1,	vs10

	xsmulsp		vs38,	vs0,	vs11
	xsmulsp		vs39,	vs1,	vs11

	xsmulsp		vs40,	vs0,	vs12
	xsmulsp		vs41,	vs1,	vs12

	xsmulsp		vs42,	vs0,	vs13
	xsmulsp		vs43,	vs1,	vs13

	xsmulsp		vs44,	vs0,	vs14
	xsmulsp		vs45,	vs1,	vs14

	xsmulsp		vs46,	vs0,	vs15
	xsmulsp		vs47,	vs1,	vs15


.endm

.macro KERNEL8x2_SUB1


	lxsspx		vs0,	o0,	AO
	lxsspx		vs1,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1
	lxsspx		vs10,	o8,	T1
	lxsspx		vs11,	o12,	T1

	addi		T1,	T1,	16

	lxsspx		vs12,	o0,	T1
	lxsspx		vs13,	o4,	T1
	lxsspx		vs14,	o8,	T1
	lxsspx		vs15,	o12,	T1

	addi		BO,	BO,	32


	xsmaddasp	vs32,	vs0,	vs8
	xsmaddasp	vs33,	vs1,	vs8

	xsmaddasp	vs34,	vs0,	vs9
	xsmaddasp	vs35,	vs1,	vs9

	xsmaddasp	vs36,	vs0,	vs10
	xsmaddasp	vs37,	vs1,	vs10

	xsmaddasp	vs38,	vs0,	vs11
	xsmaddasp	vs39,	vs1,	vs11

	xsmaddasp	vs40,	vs0,	vs12
	xsmaddasp	vs41,	vs1,	vs12

	xsmaddasp	vs42,	vs0,	vs13
	xsmaddasp	vs43,	vs1,	vs13

	xsmaddasp	vs44,	vs0,	vs14
	xsmaddasp	vs45,	vs1,	vs14

	xsmaddasp	vs46,	vs0,	vs15
	xsmaddasp	vs47,	vs1,	vs15


.endm

.macro SAVE8x2

	mr		T1,	CO

#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1
	lxsspx		vs1,	o4,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs32,	alpha_r
	xsmulsp		vs1,	vs33,	alpha_r

#else

	xsmulsp		vs28,	vs32,	alpha_r
	xsaddsp		vs0,	vs0,	vs28
	xsmulsp		vs28,	vs33,	alpha_r
	xsaddsp		vs1,	vs1,	vs28

#endif

	stxsspx		vs0,	o0,	T1
	stxsspx		vs1,	o4,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1
	lxsspx		vs1,	o4,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs34,	alpha_r
	xsmulsp		vs1,	vs35,	alpha_r

#else

	xsmulsp		vs28,	vs34,	alpha_r
	xsaddsp		vs0,	vs0,	vs28
	xsmulsp		vs28,	vs35,	alpha_r
	xsaddsp		vs1,	vs1,	vs28

#endif

	stxsspx		vs0,	o0,	T1
	stxsspx		vs1,	o4,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1
	lxsspx		vs1,	o4,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs36,	alpha_r
	xsmulsp		vs1,	vs37,	alpha_r

#else

	xsmulsp		vs28,	vs36,	alpha_r
	xsaddsp		vs0,	vs0,	vs28
	xsmulsp		vs28,	vs37,	alpha_r
	xsaddsp		vs1,	vs1,	vs28

#endif

	stxsspx		vs0,	o0,	T1
	stxsspx		vs1,	o4,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1
	lxsspx		vs1,	o4,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs38,	alpha_r
	xsmulsp		vs1,	vs39,	alpha_r

#else

	xsmulsp		vs28,	vs38,	alpha_r
	xsaddsp		vs0,	vs0,	vs28
	xsmulsp		vs28,	vs39,	alpha_r
	xsaddsp		vs1,	vs1,	vs28

#endif

	stxsspx		vs0,	o0,	T1
	stxsspx		vs1,	o4,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1
	lxsspx		vs1,	o4,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs40,	alpha_r
	xsmulsp		vs1,	vs41,	alpha_r

#else

	xsmulsp		vs28,	vs40,	alpha_r
	xsaddsp		vs0,	vs0,	vs28
	xsmulsp		vs28,	vs41,	alpha_r
	xsaddsp		vs1,	vs1,	vs28

#endif

	stxsspx		vs0,	o0,	T1
	stxsspx		vs1,	o4,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1
	lxsspx		vs1,	o4,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs42,	alpha_r
	xsmulsp		vs1,	vs43,	alpha_r

#else

	xsmulsp		vs28,	vs42,	alpha_r
	xsaddsp		vs0,	vs0,	vs28
	xsmulsp		vs28,	vs43,	alpha_r
	xsaddsp		vs1,	vs1,	vs28

#endif

	stxsspx		vs0,	o0,	T1
	stxsspx		vs1,	o4,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1
	lxsspx		vs1,	o4,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs44,	alpha_r
	xsmulsp		vs1,	vs45,	alpha_r

#else

	xsmulsp		vs28,	vs44,	alpha_r
	xsaddsp		vs0,	vs0,	vs28
	xsmulsp		vs28,	vs45,	alpha_r
	xsaddsp		vs1,	vs1,	vs28

#endif

	stxsspx		vs0,	o0,	T1
	stxsspx		vs1,	o4,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1
	lxsspx		vs1,	o4,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs46,	alpha_r
	xsmulsp		vs1,	vs47,	alpha_r

#else

	xsmulsp		vs28,	vs46,	alpha_r
	xsaddsp		vs0,	vs0,	vs28
	xsmulsp		vs28,	vs47,	alpha_r
	xsaddsp		vs1,	vs1,	vs28

#endif

	stxsspx		vs0,	o0,	T1
	stxsspx		vs1,	o4,	T1

	add		T1,	T1,	LDC

	addi		CO,	CO,	8

.endm


/**********************************************************************************************
* Macros for N=8 and M=1
**********************************************************************************************/

.macro LOAD8x1_1

	lxsspx		vs0,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1
	lxsspx		vs10,	o8,	T1
	lxsspx		vs11,	o12,	T1

	addi		T1,	T1,	16

	lxsspx		vs12,	o0,	T1
	lxsspx		vs13,	o4,	T1
	lxsspx		vs14,	o8,	T1
	lxsspx		vs15,	o12,	T1

	addi		BO,	BO,	32

.endm

.macro KERNEL8x1_I1


	lxsspx		vs4,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1
	lxsspx		vs17,	o4,	T1
	lxsspx		vs18,	o8,	T1
	lxsspx		vs19,	o12,	T1

	addi		T1,	T1,	16

	lxsspx		vs20,	o0,	T1
	lxsspx		vs21,	o4,	T1
	lxsspx		vs22,	o8,	T1
	lxsspx		vs23,	o12,	T1

	addi		BO,	BO,	32


	xsmulsp		vs32,	vs0,	vs8

	xsmulsp		vs33,	vs0,	vs9

	xsmulsp		vs34,	vs0,	vs10

	xsmulsp		vs35,	vs0,	vs11

	xsmulsp		vs36,	vs0,	vs12

	xsmulsp		vs37,	vs0,	vs13

	xsmulsp		vs38,	vs0,	vs14

	xsmulsp		vs39,	vs0,	vs15


.endm

.macro KERNEL8x1_1


	lxsspx		vs4,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1
	lxsspx		vs17,	o4,	T1
	lxsspx		vs18,	o8,	T1
	lxsspx		vs19,	o12,	T1

	addi		T1,	T1,	16

	lxsspx		vs20,	o0,	T1
	lxsspx		vs21,	o4,	T1
	lxsspx		vs22,	o8,	T1
	lxsspx		vs23,	o12,	T1

	addi		BO,	BO,	32


	xsmaddasp	vs32,	vs0,	vs8

	xsmaddasp	vs33,	vs0,	vs9

	xsmaddasp	vs34,	vs0,	vs10

	xsmaddasp	vs35,	vs0,	vs11

	xsmaddasp	vs36,	vs0,	vs12

	xsmaddasp	vs37,	vs0,	vs13

	xsmaddasp	vs38,	vs0,	vs14

	xsmaddasp	vs39,	vs0,	vs15


.endm

.macro KERNEL8x1_2


	lxsspx		vs0,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1
	lxsspx		vs10,	o8,	T1
	lxsspx		vs11,	o12,	T1

	addi		T1,	T1,	16

	lxsspx		vs12,	o0,	T1
	lxsspx		vs13,	o4,	T1
	lxsspx		vs14,	o8,	T1
	lxsspx		vs15,	o12,	T1

	addi		BO,	BO,	32


	xsmaddasp	vs32,	vs4,	vs16

	xsmaddasp	vs33,	vs4,	vs17

	xsmaddasp	vs34,	vs4,	vs18

	xsmaddasp	vs35,	vs4,	vs19

	xsmaddasp	vs36,	vs4,	vs20

	xsmaddasp	vs37,	vs4,	vs21

	xsmaddasp	vs38,	vs4,	vs22

	xsmaddasp	vs39,	vs4,	vs23


.endm

.macro KERNEL8x1_E2


	xsmaddasp	vs32,	vs4,	vs16

	xsmaddasp	vs33,	vs4,	vs17

	xsmaddasp	vs34,	vs4,	vs18

	xsmaddasp	vs35,	vs4,	vs19

	xsmaddasp	vs36,	vs4,	vs20

	xsmaddasp	vs37,	vs4,	vs21

	xsmaddasp	vs38,	vs4,	vs22

	xsmaddasp	vs39,	vs4,	vs23


.endm

.macro KERNEL8x1_SUBI1


	lxsspx		vs0,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1
	lxsspx		vs10,	o8,	T1
	lxsspx		vs11,	o12,	T1

	addi		T1,	T1,	16

	lxsspx		vs12,	o0,	T1
	lxsspx		vs13,	o4,	T1
	lxsspx		vs14,	o8,	T1
	lxsspx		vs15,	o12,	T1

	addi		BO,	BO,	32


	xsmulsp		vs32,	vs0,	vs8

	xsmulsp		vs33,	vs0,	vs9

	xsmulsp		vs34,	vs0,	vs10

	xsmulsp		vs35,	vs0,	vs11

	xsmulsp		vs36,	vs0,	vs12

	xsmulsp		vs37,	vs0,	vs13

	xsmulsp		vs38,	vs0,	vs14

	xsmulsp		vs39,	vs0,	vs15


.endm

.macro KERNEL8x1_SUB1


	lxsspx		vs0,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1
	lxsspx		vs10,	o8,	T1
	lxsspx		vs11,	o12,	T1

	addi		T1,	T1,	16

	lxsspx		vs12,	o0,	T1
	lxsspx		vs13,	o4,	T1
	lxsspx		vs14,	o8,	T1
	lxsspx		vs15,	o12,	T1

	addi		BO,	BO,	32


	xsmaddasp	vs32,	vs0,	vs8

	xsmaddasp	vs33,	vs0,	vs9

	xsmaddasp	vs34,	vs0,	vs10

	xsmaddasp	vs35,	vs0,	vs11

	xsmaddasp	vs36,	vs0,	vs12

	xsmaddasp	vs37,	vs0,	vs13

	xsmaddasp	vs38,	vs0,	vs14

	xsmaddasp	vs39,	vs0,	vs15


.endm

.macro SAVE8x1

	mr		T1,	CO

#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs32,	alpha_r

#else

	xsmulsp		vs28,	vs32,	alpha_r
	xsaddsp		vs0,	vs0,	vs28

#endif

	stxsspx		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs33,	alpha_r

#else

	xsmulsp		vs28,	vs33,	alpha_r
	xsaddsp		vs0,	vs0,	vs28

#endif

	stxsspx		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs34,	alpha_r

#else

	xsmulsp		vs28,	vs34,	alpha_r
	xsaddsp		vs0,	vs0,	vs28

#endif

	stxsspx		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs35,	alpha_r

#else

	xsmulsp		vs28,	vs35,	alpha_r
	xsaddsp		vs0,	vs0,	vs28

#endif

	stxsspx		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs36,	alpha_r

#else

	xsmulsp		vs28,	vs36,	alpha_r
	xsaddsp		vs0,	vs0,	vs28

#endif

	stxsspx		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs37,	alpha_r

#else

	xsmulsp		vs28,	vs37,	alpha_r
	xsaddsp		vs0,	vs0,	vs28

#endif

	stxsspx		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs38,	alpha_r

#else

	xsmulsp		vs28,	vs38,	alpha_r
	xsaddsp		vs0,	vs0,	vs28

#endif

	stxsspx		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs39,	alpha_r

#else

	xsmulsp		vs28,	vs39,	alpha_r
	xsaddsp		vs0,	vs0,	vs28

#endif

	stxsspx		vs0,	o0,	T1

	add		T1,	T1,	LDC

	addi		CO,	CO,	4

.endm


/**********************************************************************************************
* Macros for N=4 and M=16
**********************************************************************************************/

.macro LOAD4x16_1

	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO
	lxvw4x		vs2,	o32,	AO
	lxvw4x		vs3,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	addi		BO,	BO,	16

.endm

.macro KERNEL4x16_I1


	lxvw4x		vs4,	o0,	AO
	lxvw4x		vs5,	o16,	AO
	lxvw4x		vs6,	o32,	AO
	lxvw4x		vs7,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0
	xxspltw		vs17,	vs28,	1
	xxspltw		vs18,	vs28,	2
	xxspltw		vs19,	vs28,	3

	addi		BO,	BO,	16


	xvmulsp		vs32,	vs0,	vs8
	xvmulsp		vs33,	vs1,	vs8
	xvmulsp		vs34,	vs2,	vs8
	xvmulsp		vs35,	vs3,	vs8

	xvmulsp		vs36,	vs0,	vs9
	xvmulsp		vs37,	vs1,	vs9
	xvmulsp		vs38,	vs2,	vs9
	xvmulsp		vs39,	vs3,	vs9

	xvmulsp		vs40,	vs0,	vs10
	xvmulsp		vs41,	vs1,	vs10
	xvmulsp		vs42,	vs2,	vs10
	xvmulsp		vs43,	vs3,	vs10

	xvmulsp		vs44,	vs0,	vs11
	xvmulsp		vs45,	vs1,	vs11
	xvmulsp		vs46,	vs2,	vs11
	xvmulsp		vs47,	vs3,	vs11


.endm

.macro KERNEL4x16_1


	lxvw4x		vs4,	o0,	AO
	lxvw4x		vs5,	o16,	AO
	lxvw4x		vs6,	o32,	AO
	lxvw4x		vs7,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0
	xxspltw		vs17,	vs28,	1
	xxspltw		vs18,	vs28,	2
	xxspltw		vs19,	vs28,	3

	addi		BO,	BO,	16


	xvmaddasp	vs32,	vs0,	vs8
	xvmaddasp	vs33,	vs1,	vs8
	xvmaddasp	vs34,	vs2,	vs8
	xvmaddasp	vs35,	vs3,	vs8

	xvmaddasp	vs36,	vs0,	vs9
	xvmaddasp	vs37,	vs1,	vs9
	xvmaddasp	vs38,	vs2,	vs9
	xvmaddasp	vs39,	vs3,	vs9

	xvmaddasp	vs40,	vs0,	vs10
	xvmaddasp	vs41,	vs1,	vs10
	xvmaddasp	vs42,	vs2,	vs10
	xvmaddasp	vs43,	vs3,	vs10

	xvmaddasp	vs44,	vs0,	vs11
	xvmaddasp	vs45,	vs1,	vs11
	xvmaddasp	vs46,	vs2,	vs11
	xvmaddasp	vs47,	vs3,	vs11


.endm

.macro KERNEL4x16_2


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO
	lxvw4x		vs2,	o32,	AO
	lxvw4x		vs3,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	addi		BO,	BO,	16


	xvmaddasp	vs32,	vs4,	vs16
	xvmaddasp	vs33,	vs5,	vs16
	xvmaddasp	vs34,	vs6,	vs16
	xvmaddasp	vs35,	vs7,	vs16

	xvmaddasp	vs36,	vs4,	vs17
	xvmaddasp	vs37,	vs5,	vs17
	xvmaddasp	vs38,	vs6,	vs17
	xvmaddasp	vs39,	vs7,	vs17

	xvmaddasp	vs40,	vs4,	vs18
	xvmaddasp	vs41,	vs5,	vs18
	xvmaddasp	vs42,	vs6,	vs18
	xvmaddasp	vs43,	vs7,	vs18

	xvmaddasp	vs44,	vs4,	vs19
	xvmaddasp	vs45,	vs5,	vs19
	xvmaddasp	vs46,	vs6,	vs19
	xvmaddasp	vs47,	vs7,	vs19


.endm

.macro KERNEL4x16_E2


	xvmaddasp	vs32,	vs4,	vs16
	xvmaddasp	vs33,	vs5,	vs16
	xvmaddasp	vs34,	vs6,	vs16
	xvmaddasp	vs35,	vs7,	vs16

	xvmaddasp	vs36,	vs4,	vs17
	xvmaddasp	vs37,	vs5,	vs17
	xvmaddasp	vs38,	vs6,	vs17
	xvmaddasp	vs39,	vs7,	vs17

	xvmaddasp	vs40,	vs4,	vs18
	xvmaddasp	vs41,	vs5,	vs18
	xvmaddasp	vs42,	vs6,	vs18
	xvmaddasp	vs43,	vs7,	vs18

	xvmaddasp	vs44,	vs4,	vs19
	xvmaddasp	vs45,	vs5,	vs19
	xvmaddasp	vs46,	vs6,	vs19
	xvmaddasp	vs47,	vs7,	vs19


.endm

.macro KERNEL4x16_SUBI1


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO
	lxvw4x		vs2,	o32,	AO
	lxvw4x		vs3,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	addi		BO,	BO,	16


	xvmulsp		vs32,	vs0,	vs8
	xvmulsp		vs33,	vs1,	vs8
	xvmulsp		vs34,	vs2,	vs8
	xvmulsp		vs35,	vs3,	vs8

	xvmulsp		vs36,	vs0,	vs9
	xvmulsp		vs37,	vs1,	vs9
	xvmulsp		vs38,	vs2,	vs9
	xvmulsp		vs39,	vs3,	vs9

	xvmulsp		vs40,	vs0,	vs10
	xvmulsp		vs41,	vs1,	vs10
	xvmulsp		vs42,	vs2,	vs10
	xvmulsp		vs43,	vs3,	vs10

	xvmulsp		vs44,	vs0,	vs11
	xvmulsp		vs45,	vs1,	vs11
	xvmulsp		vs46,	vs2,	vs11
	xvmulsp		vs47,	vs3,	vs11


.endm

.macro KERNEL4x16_SUB1


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO
	lxvw4x		vs2,	o32,	AO
	lxvw4x		vs3,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	addi		BO,	BO,	16


	xvmaddasp	vs32,	vs0,	vs8
	xvmaddasp	vs33,	vs1,	vs8
	xvmaddasp	vs34,	vs2,	vs8
	xvmaddasp	vs35,	vs3,	vs8

	xvmaddasp	vs36,	vs0,	vs9
	xvmaddasp	vs37,	vs1,	vs9
	xvmaddasp	vs38,	vs2,	vs9
	xvmaddasp	vs39,	vs3,	vs9

	xvmaddasp	vs40,	vs0,	vs10
	xvmaddasp	vs41,	vs1,	vs10
	xvmaddasp	vs42,	vs2,	vs10
	xvmaddasp	vs43,	vs3,	vs10

	xvmaddasp	vs44,	vs0,	vs11
	xvmaddasp	vs45,	vs1,	vs11
	xvmaddasp	vs46,	vs2,	vs11
	xvmaddasp	vs47,	vs3,	vs11


.endm

.macro SAVE4x16

	mr		T1,	CO

#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1
	lxvw4x		vs2,	o32,	T1
	lxvw4x		vs3,	o48,	T1

#endif


	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs33,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif

	stxvw4x		vs34,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs2,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs2,	vs2,	vs28
#endif

	stxvw4x		vs35,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs3,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs3,	vs3,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1
	stxvw4x		vs2,	o32,	T1
	stxvw4x		vs3,	o48,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1
	lxvw4x		vs2,	o32,	T1
	lxvw4x		vs3,	o48,	T1

#endif


	stxvw4x		vs36,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs37,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif

	stxvw4x		vs38,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs2,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs2,	vs2,	vs28
#endif

	stxvw4x		vs39,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs3,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs3,	vs3,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1
	stxvw4x		vs2,	o32,	T1
	stxvw4x		vs3,	o48,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1
	lxvw4x		vs2,	o32,	T1
	lxvw4x		vs3,	o48,	T1

#endif


	stxvw4x		vs40,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs41,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif

	stxvw4x		vs42,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs2,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs2,	vs2,	vs28
#endif

	stxvw4x		vs43,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs3,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs3,	vs3,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1
	stxvw4x		vs2,	o32,	T1
	stxvw4x		vs3,	o48,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1
	lxvw4x		vs2,	o32,	T1
	lxvw4x		vs3,	o48,	T1

#endif


	stxvw4x		vs44,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs45,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif

	stxvw4x		vs46,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs2,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs2,	vs2,	vs28
#endif

	stxvw4x		vs47,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs3,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs3,	vs3,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1
	stxvw4x		vs2,	o32,	T1
	stxvw4x		vs3,	o48,	T1

	add		T1,	T1,	LDC

	addi		CO,	CO,	64

.endm


/**********************************************************************************************
* Macros for N=4 and M=8
**********************************************************************************************/

.macro LOAD4x8_1

	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	addi		BO,	BO,	16

.endm

.macro KERNEL4x8_I1


	lxvw4x		vs4,	o0,	AO
	lxvw4x		vs5,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0
	xxspltw		vs17,	vs28,	1
	xxspltw		vs18,	vs28,	2
	xxspltw		vs19,	vs28,	3

	addi		BO,	BO,	16


	xvmulsp		vs32,	vs0,	vs8
	xvmulsp		vs33,	vs1,	vs8

	xvmulsp		vs34,	vs0,	vs9
	xvmulsp		vs35,	vs1,	vs9

	xvmulsp		vs36,	vs0,	vs10
	xvmulsp		vs37,	vs1,	vs10

	xvmulsp		vs38,	vs0,	vs11
	xvmulsp		vs39,	vs1,	vs11


.endm

.macro KERNEL4x8_1


	lxvw4x		vs4,	o0,	AO
	lxvw4x		vs5,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0
	xxspltw		vs17,	vs28,	1
	xxspltw		vs18,	vs28,	2
	xxspltw		vs19,	vs28,	3

	addi		BO,	BO,	16


	xvmaddasp	vs32,	vs0,	vs8
	xvmaddasp	vs33,	vs1,	vs8

	xvmaddasp	vs34,	vs0,	vs9
	xvmaddasp	vs35,	vs1,	vs9

	xvmaddasp	vs36,	vs0,	vs10
	xvmaddasp	vs37,	vs1,	vs10

	xvmaddasp	vs38,	vs0,	vs11
	xvmaddasp	vs39,	vs1,	vs11


.endm

.macro KERNEL4x8_2


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	addi		BO,	BO,	16


	xvmaddasp	vs32,	vs4,	vs16
	xvmaddasp	vs33,	vs5,	vs16

	xvmaddasp	vs34,	vs4,	vs17
	xvmaddasp	vs35,	vs5,	vs17

	xvmaddasp	vs36,	vs4,	vs18
	xvmaddasp	vs37,	vs5,	vs18

	xvmaddasp	vs38,	vs4,	vs19
	xvmaddasp	vs39,	vs5,	vs19


.endm

.macro KERNEL4x8_E2


	xvmaddasp	vs32,	vs4,	vs16
	xvmaddasp	vs33,	vs5,	vs16

	xvmaddasp	vs34,	vs4,	vs17
	xvmaddasp	vs35,	vs5,	vs17

	xvmaddasp	vs36,	vs4,	vs18
	xvmaddasp	vs37,	vs5,	vs18

	xvmaddasp	vs38,	vs4,	vs19
	xvmaddasp	vs39,	vs5,	vs19


.endm

.macro KERNEL4x8_SUBI1


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	addi		BO,	BO,	16


	xvmulsp		vs32,	vs0,	vs8
	xvmulsp		vs33,	vs1,	vs8

	xvmulsp		vs34,	vs0,	vs9
	xvmulsp		vs35,	vs1,	vs9

	xvmulsp		vs36,	vs0,	vs10
	xvmulsp		vs37,	vs1,	vs10

	xvmulsp		vs38,	vs0,	vs11
	xvmulsp		vs39,	vs1,	vs11


.endm

.macro KERNEL4x8_SUB1


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	addi		BO,	BO,	16


	xvmaddasp	vs32,	vs0,	vs8
	xvmaddasp	vs33,	vs1,	vs8

	xvmaddasp	vs34,	vs0,	vs9
	xvmaddasp	vs35,	vs1,	vs9

	xvmaddasp	vs36,	vs0,	vs10
	xvmaddasp	vs37,	vs1,	vs10

	xvmaddasp	vs38,	vs0,	vs11
	xvmaddasp	vs39,	vs1,	vs11


.endm

.macro SAVE4x8

	mr		T1,	CO

#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1

#endif


	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs33,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1

#endif


	stxvw4x		vs34,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs35,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1

#endif


	stxvw4x		vs36,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs37,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1

#endif


	stxvw4x		vs38,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs39,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1

	add		T1,	T1,	LDC

	addi		CO,	CO,	32

.endm


/**********************************************************************************************
* Macros for N=4 and M=4
**********************************************************************************************/

.macro LOAD4x4_1

	lxvw4x		vs0,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	addi		BO,	BO,	16

.endm

.macro KERNEL4x4_I1


	lxvw4x		vs4,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0
	xxspltw		vs17,	vs28,	1
	xxspltw		vs18,	vs28,	2
	xxspltw		vs19,	vs28,	3

	addi		BO,	BO,	16


	xvmulsp		vs32,	vs0,	vs8

	xvmulsp		vs33,	vs0,	vs9

	xvmulsp		vs34,	vs0,	vs10

	xvmulsp		vs35,	vs0,	vs11


.endm

.macro KERNEL4x4_1


	lxvw4x		vs4,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0
	xxspltw		vs17,	vs28,	1
	xxspltw		vs18,	vs28,	2
	xxspltw		vs19,	vs28,	3

	addi		BO,	BO,	16


	xvmaddasp	vs32,	vs0,	vs8

	xvmaddasp	vs33,	vs0,	vs9

	xvmaddasp	vs34,	vs0,	vs10

	xvmaddasp	vs35,	vs0,	vs11


.endm

.macro KERNEL4x4_2


	lxvw4x		vs0,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	addi		BO,	BO,	16


	xvmaddasp	vs32,	vs4,	vs16

	xvmaddasp	vs33,	vs4,	vs17

	xvmaddasp	vs34,	vs4,	vs18

	xvmaddasp	vs35,	vs4,	vs19


.endm

.macro KERNEL4x4_E2


	xvmaddasp	vs32,	vs4,	vs16

	xvmaddasp	vs33,	vs4,	vs17

	xvmaddasp	vs34,	vs4,	vs18

	xvmaddasp	vs35,	vs4,	vs19


.endm

.macro KERNEL4x4_SUBI1


	lxvw4x		vs0,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	addi		BO,	BO,	16


	xvmulsp		vs32,	vs0,	vs8

	xvmulsp		vs33,	vs0,	vs9

	xvmulsp		vs34,	vs0,	vs10

	xvmulsp		vs35,	vs0,	vs11


.endm

.macro KERNEL4x4_SUB1


	lxvw4x		vs0,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1
	xxspltw		vs10,	vs28,	2
	xxspltw		vs11,	vs28,	3

	addi		BO,	BO,	16


	xvmaddasp	vs32,	vs0,	vs8

	xvmaddasp	vs33,	vs0,	vs9

	xvmaddasp	vs34,	vs0,	vs10

	xvmaddasp	vs35,	vs0,	vs11


.endm

.macro SAVE4x4

	mr		T1,	CO

#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1

#endif


	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif




	stxvw4x		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1

#endif


	stxvw4x		vs33,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif




	stxvw4x		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1

#endif


	stxvw4x		vs34,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif




	stxvw4x		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1

#endif


	stxvw4x		vs35,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif




	stxvw4x		vs0,	o0,	T1

	add		T1,	T1,	LDC

	addi		CO,	CO,	16

.endm


/**********************************************************************************************
* Macros for N=4 and M=2
**********************************************************************************************/

.macro LOAD4x2_1

	lxsspx		vs0,	o0,	AO
	lxsspx		vs1,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1
	lxsspx		vs10,	o8,	T1
	lxsspx		vs11,	o12,	T1

	addi		BO,	BO,	16

.endm

.macro KERNEL4x2_I1


	lxsspx		vs4,	o0,	AO
	lxsspx		vs5,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1
	lxsspx		vs17,	o4,	T1
	lxsspx		vs18,	o8,	T1
	lxsspx		vs19,	o12,	T1

	addi		BO,	BO,	16


	xsmulsp		vs32,	vs0,	vs8
	xsmulsp		vs33,	vs1,	vs8

	xsmulsp		vs34,	vs0,	vs9
	xsmulsp		vs35,	vs1,	vs9

	xsmulsp		vs36,	vs0,	vs10
	xsmulsp		vs37,	vs1,	vs10

	xsmulsp		vs38,	vs0,	vs11
	xsmulsp		vs39,	vs1,	vs11


.endm

.macro KERNEL4x2_1


	lxsspx		vs4,	o0,	AO
	lxsspx		vs5,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1
	lxsspx		vs17,	o4,	T1
	lxsspx		vs18,	o8,	T1
	lxsspx		vs19,	o12,	T1

	addi		BO,	BO,	16


	xsmaddasp	vs32,	vs0,	vs8
	xsmaddasp	vs33,	vs1,	vs8

	xsmaddasp	vs34,	vs0,	vs9
	xsmaddasp	vs35,	vs1,	vs9

	xsmaddasp	vs36,	vs0,	vs10
	xsmaddasp	vs37,	vs1,	vs10

	xsmaddasp	vs38,	vs0,	vs11
	xsmaddasp	vs39,	vs1,	vs11


.endm

.macro KERNEL4x2_2


	lxsspx		vs0,	o0,	AO
	lxsspx		vs1,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1
	lxsspx		vs10,	o8,	T1
	lxsspx		vs11,	o12,	T1

	addi		BO,	BO,	16


	xsmaddasp	vs32,	vs4,	vs16
	xsmaddasp	vs33,	vs5,	vs16

	xsmaddasp	vs34,	vs4,	vs17
	xsmaddasp	vs35,	vs5,	vs17

	xsmaddasp	vs36,	vs4,	vs18
	xsmaddasp	vs37,	vs5,	vs18

	xsmaddasp	vs38,	vs4,	vs19
	xsmaddasp	vs39,	vs5,	vs19


.endm

.macro KERNEL4x2_E2


	xsmaddasp	vs32,	vs4,	vs16
	xsmaddasp	vs33,	vs5,	vs16

	xsmaddasp	vs34,	vs4,	vs17
	xsmaddasp	vs35,	vs5,	vs17

	xsmaddasp	vs36,	vs4,	vs18
	xsmaddasp	vs37,	vs5,	vs18

	xsmaddasp	vs38,	vs4,	vs19
	xsmaddasp	vs39,	vs5,	vs19


.endm

.macro KERNEL4x2_SUBI1


	lxsspx		vs0,	o0,	AO
	lxsspx		vs1,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1
	lxsspx		vs10,	o8,	T1
	lxsspx		vs11,	o12,	T1

	addi		BO,	BO,	16


	xsmulsp		vs32,	vs0,	vs8
	xsmulsp		vs33,	vs1,	vs8

	xsmulsp		vs34,	vs0,	vs9
	xsmulsp		vs35,	vs1,	vs9

	xsmulsp		vs36,	vs0,	vs10
	xsmulsp		vs37,	vs1,	vs10

	xsmulsp		vs38,	vs0,	vs11
	xsmulsp		vs39,	vs1,	vs11


.endm

.macro KERNEL4x2_SUB1


	lxsspx		vs0,	o0,	AO
	lxsspx		vs1,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1
	lxsspx		vs10,	o8,	T1
	lxsspx		vs11,	o12,	T1

	addi		BO,	BO,	16


	xsmaddasp	vs32,	vs0,	vs8
	xsmaddasp	vs33,	vs1,	vs8

	xsmaddasp	vs34,	vs0,	vs9
	xsmaddasp	vs35,	vs1,	vs9

	xsmaddasp	vs36,	vs0,	vs10
	xsmaddasp	vs37,	vs1,	vs10

	xsmaddasp	vs38,	vs0,	vs11
	xsmaddasp	vs39,	vs1,	vs11


.endm

.macro SAVE4x2

	mr		T1,	CO

#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1
	lxsspx		vs1,	o4,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs32,	alpha_r
	xsmulsp		vs1,	vs33,	alpha_r

#else

	xsmulsp		vs28,	vs32,	alpha_r
	xsaddsp		vs0,	vs0,	vs28
	xsmulsp		vs28,	vs33,	alpha_r
	xsaddsp		vs1,	vs1,	vs28

#endif

	stxsspx		vs0,	o0,	T1
	stxsspx		vs1,	o4,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1
	lxsspx		vs1,	o4,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs34,	alpha_r
	xsmulsp		vs1,	vs35,	alpha_r

#else

	xsmulsp		vs28,	vs34,	alpha_r
	xsaddsp		vs0,	vs0,	vs28
	xsmulsp		vs28,	vs35,	alpha_r
	xsaddsp		vs1,	vs1,	vs28

#endif

	stxsspx		vs0,	o0,	T1
	stxsspx		vs1,	o4,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1
	lxsspx		vs1,	o4,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs36,	alpha_r
	xsmulsp		vs1,	vs37,	alpha_r

#else

	xsmulsp		vs28,	vs36,	alpha_r
	xsaddsp		vs0,	vs0,	vs28
	xsmulsp		vs28,	vs37,	alpha_r
	xsaddsp		vs1,	vs1,	vs28

#endif

	stxsspx		vs0,	o0,	T1
	stxsspx		vs1,	o4,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1
	lxsspx		vs1,	o4,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs38,	alpha_r
	xsmulsp		vs1,	vs39,	alpha_r

#else

	xsmulsp		vs28,	vs38,	alpha_r
	xsaddsp		vs0,	vs0,	vs28
	xsmulsp		vs28,	vs39,	alpha_r
	xsaddsp		vs1,	vs1,	vs28

#endif

	stxsspx		vs0,	o0,	T1
	stxsspx		vs1,	o4,	T1

	add		T1,	T1,	LDC

	addi		CO,	CO,	8

.endm


/**********************************************************************************************
* Macros for N=4 and M=1
**********************************************************************************************/

.macro LOAD4x1_1

	lxsspx		vs0,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1
	lxsspx		vs10,	o8,	T1
	lxsspx		vs11,	o12,	T1

	addi		BO,	BO,	16

.endm

.macro KERNEL4x1_I1


	lxsspx		vs4,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1
	lxsspx		vs17,	o4,	T1
	lxsspx		vs18,	o8,	T1
	lxsspx		vs19,	o12,	T1

	addi		BO,	BO,	16


	xsmulsp		vs32,	vs0,	vs8

	xsmulsp		vs33,	vs0,	vs9

	xsmulsp		vs34,	vs0,	vs10

	xsmulsp		vs35,	vs0,	vs11


.endm

.macro KERNEL4x1_1


	lxsspx		vs4,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1
	lxsspx		vs17,	o4,	T1
	lxsspx		vs18,	o8,	T1
	lxsspx		vs19,	o12,	T1

	addi		BO,	BO,	16


	xsmaddasp	vs32,	vs0,	vs8

	xsmaddasp	vs33,	vs0,	vs9

	xsmaddasp	vs34,	vs0,	vs10

	xsmaddasp	vs35,	vs0,	vs11


.endm

.macro KERNEL4x1_2


	lxsspx		vs0,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1
	lxsspx		vs10,	o8,	T1
	lxsspx		vs11,	o12,	T1

	addi		BO,	BO,	16


	xsmaddasp	vs32,	vs4,	vs16

	xsmaddasp	vs33,	vs4,	vs17

	xsmaddasp	vs34,	vs4,	vs18

	xsmaddasp	vs35,	vs4,	vs19


.endm

.macro KERNEL4x1_E2


	xsmaddasp	vs32,	vs4,	vs16

	xsmaddasp	vs33,	vs4,	vs17

	xsmaddasp	vs34,	vs4,	vs18

	xsmaddasp	vs35,	vs4,	vs19


.endm

.macro KERNEL4x1_SUBI1


	lxsspx		vs0,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1
	lxsspx		vs10,	o8,	T1
	lxsspx		vs11,	o12,	T1

	addi		BO,	BO,	16


	xsmulsp		vs32,	vs0,	vs8

	xsmulsp		vs33,	vs0,	vs9

	xsmulsp		vs34,	vs0,	vs10

	xsmulsp		vs35,	vs0,	vs11


.endm

.macro KERNEL4x1_SUB1


	lxsspx		vs0,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1
	lxsspx		vs10,	o8,	T1
	lxsspx		vs11,	o12,	T1

	addi		BO,	BO,	16


	xsmaddasp	vs32,	vs0,	vs8

	xsmaddasp	vs33,	vs0,	vs9

	xsmaddasp	vs34,	vs0,	vs10

	xsmaddasp	vs35,	vs0,	vs11


.endm

.macro SAVE4x1

	mr		T1,	CO

#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs32,	alpha_r

#else

	xsmulsp		vs28,	vs32,	alpha_r
	xsaddsp		vs0,	vs0,	vs28

#endif

	stxsspx		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs33,	alpha_r

#else

	xsmulsp		vs28,	vs33,	alpha_r
	xsaddsp		vs0,	vs0,	vs28

#endif

	stxsspx		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs34,	alpha_r

#else

	xsmulsp		vs28,	vs34,	alpha_r
	xsaddsp		vs0,	vs0,	vs28

#endif

	stxsspx		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs35,	alpha_r

#else

	xsmulsp		vs28,	vs35,	alpha_r
	xsaddsp		vs0,	vs0,	vs28

#endif

	stxsspx		vs0,	o0,	T1

	add		T1,	T1,	LDC

	addi		CO,	CO,	4

.endm


/**********************************************************************************************
* Macros for N=2 and M=16
**********************************************************************************************/

.macro LOAD2x16_1

	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO
	lxvw4x		vs2,	o32,	AO
	lxvw4x		vs3,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1

	addi		BO,	BO,	8

.endm

.macro KERNEL2x16_I1


	lxvw4x		vs4,	o0,	AO
	lxvw4x		vs5,	o16,	AO
	lxvw4x		vs6,	o32,	AO
	lxvw4x		vs7,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0
	xxspltw		vs17,	vs28,	1

	addi		BO,	BO,	8


	xvmulsp		vs32,	vs0,	vs8
	xvmulsp		vs33,	vs1,	vs8
	xvmulsp		vs34,	vs2,	vs8
	xvmulsp		vs35,	vs3,	vs8

	xvmulsp		vs36,	vs0,	vs9
	xvmulsp		vs37,	vs1,	vs9
	xvmulsp		vs38,	vs2,	vs9
	xvmulsp		vs39,	vs3,	vs9


.endm

.macro KERNEL2x16_1


	lxvw4x		vs4,	o0,	AO
	lxvw4x		vs5,	o16,	AO
	lxvw4x		vs6,	o32,	AO
	lxvw4x		vs7,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0
	xxspltw		vs17,	vs28,	1

	addi		BO,	BO,	8


	xvmaddasp	vs32,	vs0,	vs8
	xvmaddasp	vs33,	vs1,	vs8
	xvmaddasp	vs34,	vs2,	vs8
	xvmaddasp	vs35,	vs3,	vs8

	xvmaddasp	vs36,	vs0,	vs9
	xvmaddasp	vs37,	vs1,	vs9
	xvmaddasp	vs38,	vs2,	vs9
	xvmaddasp	vs39,	vs3,	vs9


.endm

.macro KERNEL2x16_2


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO
	lxvw4x		vs2,	o32,	AO
	lxvw4x		vs3,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1

	addi		BO,	BO,	8


	xvmaddasp	vs32,	vs4,	vs16
	xvmaddasp	vs33,	vs5,	vs16
	xvmaddasp	vs34,	vs6,	vs16
	xvmaddasp	vs35,	vs7,	vs16

	xvmaddasp	vs36,	vs4,	vs17
	xvmaddasp	vs37,	vs5,	vs17
	xvmaddasp	vs38,	vs6,	vs17
	xvmaddasp	vs39,	vs7,	vs17


.endm

.macro KERNEL2x16_E2


	xvmaddasp	vs32,	vs4,	vs16
	xvmaddasp	vs33,	vs5,	vs16
	xvmaddasp	vs34,	vs6,	vs16
	xvmaddasp	vs35,	vs7,	vs16

	xvmaddasp	vs36,	vs4,	vs17
	xvmaddasp	vs37,	vs5,	vs17
	xvmaddasp	vs38,	vs6,	vs17
	xvmaddasp	vs39,	vs7,	vs17


.endm

.macro KERNEL2x16_SUBI1


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO
	lxvw4x		vs2,	o32,	AO
	lxvw4x		vs3,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1

	addi		BO,	BO,	8


	xvmulsp		vs32,	vs0,	vs8
	xvmulsp		vs33,	vs1,	vs8
	xvmulsp		vs34,	vs2,	vs8
	xvmulsp		vs35,	vs3,	vs8

	xvmulsp		vs36,	vs0,	vs9
	xvmulsp		vs37,	vs1,	vs9
	xvmulsp		vs38,	vs2,	vs9
	xvmulsp		vs39,	vs3,	vs9


.endm

.macro KERNEL2x16_SUB1


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO
	lxvw4x		vs2,	o32,	AO
	lxvw4x		vs3,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1

	addi		BO,	BO,	8


	xvmaddasp	vs32,	vs0,	vs8
	xvmaddasp	vs33,	vs1,	vs8
	xvmaddasp	vs34,	vs2,	vs8
	xvmaddasp	vs35,	vs3,	vs8

	xvmaddasp	vs36,	vs0,	vs9
	xvmaddasp	vs37,	vs1,	vs9
	xvmaddasp	vs38,	vs2,	vs9
	xvmaddasp	vs39,	vs3,	vs9


.endm

.macro SAVE2x16

	mr		T1,	CO

#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1
	lxvw4x		vs2,	o32,	T1
	lxvw4x		vs3,	o48,	T1

#endif


	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs33,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif

	stxvw4x		vs34,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs2,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs2,	vs2,	vs28
#endif

	stxvw4x		vs35,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs3,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs3,	vs3,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1
	stxvw4x		vs2,	o32,	T1
	stxvw4x		vs3,	o48,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1
	lxvw4x		vs2,	o32,	T1
	lxvw4x		vs3,	o48,	T1

#endif


	stxvw4x		vs36,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs37,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif

	stxvw4x		vs38,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs2,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs2,	vs2,	vs28
#endif

	stxvw4x		vs39,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs3,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs3,	vs3,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1
	stxvw4x		vs2,	o32,	T1
	stxvw4x		vs3,	o48,	T1

	add		T1,	T1,	LDC

	addi		CO,	CO,	64

.endm


/**********************************************************************************************
* Macros for N=2 and M=8
**********************************************************************************************/

.macro LOAD2x8_1

	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1

	addi		BO,	BO,	8

.endm

.macro KERNEL2x8_I1


	lxvw4x		vs4,	o0,	AO
	lxvw4x		vs5,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0
	xxspltw		vs17,	vs28,	1

	addi		BO,	BO,	8


	xvmulsp		vs32,	vs0,	vs8
	xvmulsp		vs33,	vs1,	vs8

	xvmulsp		vs34,	vs0,	vs9
	xvmulsp		vs35,	vs1,	vs9


.endm

.macro KERNEL2x8_1


	lxvw4x		vs4,	o0,	AO
	lxvw4x		vs5,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0
	xxspltw		vs17,	vs28,	1

	addi		BO,	BO,	8


	xvmaddasp	vs32,	vs0,	vs8
	xvmaddasp	vs33,	vs1,	vs8

	xvmaddasp	vs34,	vs0,	vs9
	xvmaddasp	vs35,	vs1,	vs9


.endm

.macro KERNEL2x8_2


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1

	addi		BO,	BO,	8


	xvmaddasp	vs32,	vs4,	vs16
	xvmaddasp	vs33,	vs5,	vs16

	xvmaddasp	vs34,	vs4,	vs17
	xvmaddasp	vs35,	vs5,	vs17


.endm

.macro KERNEL2x8_E2


	xvmaddasp	vs32,	vs4,	vs16
	xvmaddasp	vs33,	vs5,	vs16

	xvmaddasp	vs34,	vs4,	vs17
	xvmaddasp	vs35,	vs5,	vs17


.endm

.macro KERNEL2x8_SUBI1


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1

	addi		BO,	BO,	8


	xvmulsp		vs32,	vs0,	vs8
	xvmulsp		vs33,	vs1,	vs8

	xvmulsp		vs34,	vs0,	vs9
	xvmulsp		vs35,	vs1,	vs9


.endm

.macro KERNEL2x8_SUB1


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1

	addi		BO,	BO,	8


	xvmaddasp	vs32,	vs0,	vs8
	xvmaddasp	vs33,	vs1,	vs8

	xvmaddasp	vs34,	vs0,	vs9
	xvmaddasp	vs35,	vs1,	vs9


.endm

.macro SAVE2x8

	mr		T1,	CO

#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1

#endif


	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs33,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1

#endif


	stxvw4x		vs34,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs35,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1

	add		T1,	T1,	LDC

	addi		CO,	CO,	32

.endm


/**********************************************************************************************
* Macros for N=2 and M=4
**********************************************************************************************/

.macro LOAD2x4_1

	lxvw4x		vs0,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1

	addi		BO,	BO,	8

.endm

.macro KERNEL2x4_I1


	lxvw4x		vs4,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0
	xxspltw		vs17,	vs28,	1

	addi		BO,	BO,	8


	xvmulsp		vs32,	vs0,	vs8

	xvmulsp		vs33,	vs0,	vs9


.endm

.macro KERNEL2x4_1


	lxvw4x		vs4,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0
	xxspltw		vs17,	vs28,	1

	addi		BO,	BO,	8


	xvmaddasp	vs32,	vs0,	vs8

	xvmaddasp	vs33,	vs0,	vs9


.endm

.macro KERNEL2x4_2


	lxvw4x		vs0,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1

	addi		BO,	BO,	8


	xvmaddasp	vs32,	vs4,	vs16

	xvmaddasp	vs33,	vs4,	vs17


.endm

.macro KERNEL2x4_E2


	xvmaddasp	vs32,	vs4,	vs16

	xvmaddasp	vs33,	vs4,	vs17


.endm

.macro KERNEL2x4_SUBI1


	lxvw4x		vs0,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1

	addi		BO,	BO,	8


	xvmulsp		vs32,	vs0,	vs8

	xvmulsp		vs33,	vs0,	vs9


.endm

.macro KERNEL2x4_SUB1


	lxvw4x		vs0,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0
	xxspltw		vs9,	vs28,	1

	addi		BO,	BO,	8


	xvmaddasp	vs32,	vs0,	vs8

	xvmaddasp	vs33,	vs0,	vs9


.endm

.macro SAVE2x4

	mr		T1,	CO

#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1

#endif


	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif




	stxvw4x		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1

#endif


	stxvw4x		vs33,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif




	stxvw4x		vs0,	o0,	T1

	add		T1,	T1,	LDC

	addi		CO,	CO,	16

.endm


/**********************************************************************************************
* Macros for N=2 and M=2
**********************************************************************************************/

.macro LOAD2x2_1

	lxsspx		vs0,	o0,	AO
	lxsspx		vs1,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1

	addi		BO,	BO,	8

.endm

.macro KERNEL2x2_I1


	lxsspx		vs4,	o0,	AO
	lxsspx		vs5,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1
	lxsspx		vs17,	o4,	T1

	addi		BO,	BO,	8


	xsmulsp		vs32,	vs0,	vs8
	xsmulsp		vs33,	vs1,	vs8

	xsmulsp		vs34,	vs0,	vs9
	xsmulsp		vs35,	vs1,	vs9


.endm

.macro KERNEL2x2_1


	lxsspx		vs4,	o0,	AO
	lxsspx		vs5,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1
	lxsspx		vs17,	o4,	T1

	addi		BO,	BO,	8


	xsmaddasp	vs32,	vs0,	vs8
	xsmaddasp	vs33,	vs1,	vs8

	xsmaddasp	vs34,	vs0,	vs9
	xsmaddasp	vs35,	vs1,	vs9


.endm

.macro KERNEL2x2_2


	lxsspx		vs0,	o0,	AO
	lxsspx		vs1,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1

	addi		BO,	BO,	8


	xsmaddasp	vs32,	vs4,	vs16
	xsmaddasp	vs33,	vs5,	vs16

	xsmaddasp	vs34,	vs4,	vs17
	xsmaddasp	vs35,	vs5,	vs17


.endm

.macro KERNEL2x2_E2


	xsmaddasp	vs32,	vs4,	vs16
	xsmaddasp	vs33,	vs5,	vs16

	xsmaddasp	vs34,	vs4,	vs17
	xsmaddasp	vs35,	vs5,	vs17


.endm

.macro KERNEL2x2_SUBI1


	lxsspx		vs0,	o0,	AO
	lxsspx		vs1,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1

	addi		BO,	BO,	8


	xsmulsp		vs32,	vs0,	vs8
	xsmulsp		vs33,	vs1,	vs8

	xsmulsp		vs34,	vs0,	vs9
	xsmulsp		vs35,	vs1,	vs9


.endm

.macro KERNEL2x2_SUB1


	lxsspx		vs0,	o0,	AO
	lxsspx		vs1,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1

	addi		BO,	BO,	8


	xsmaddasp	vs32,	vs0,	vs8
	xsmaddasp	vs33,	vs1,	vs8

	xsmaddasp	vs34,	vs0,	vs9
	xsmaddasp	vs35,	vs1,	vs9


.endm

.macro SAVE2x2

	mr		T1,	CO

#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1
	lxsspx		vs1,	o4,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs32,	alpha_r
	xsmulsp		vs1,	vs33,	alpha_r

#else

	xsmulsp		vs28,	vs32,	alpha_r
	xsaddsp		vs0,	vs0,	vs28
	xsmulsp		vs28,	vs33,	alpha_r
	xsaddsp		vs1,	vs1,	vs28

#endif

	stxsspx		vs0,	o0,	T1
	stxsspx		vs1,	o4,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1
	lxsspx		vs1,	o4,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs34,	alpha_r
	xsmulsp		vs1,	vs35,	alpha_r

#else

	xsmulsp		vs28,	vs34,	alpha_r
	xsaddsp		vs0,	vs0,	vs28
	xsmulsp		vs28,	vs35,	alpha_r
	xsaddsp		vs1,	vs1,	vs28

#endif

	stxsspx		vs0,	o0,	T1
	stxsspx		vs1,	o4,	T1

	add		T1,	T1,	LDC

	addi		CO,	CO,	8

.endm


/**********************************************************************************************
* Macros for N=2 and M=1
**********************************************************************************************/

.macro LOAD2x1_1

	lxsspx		vs0,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1

	addi		BO,	BO,	8

.endm

.macro KERNEL2x1_I1


	lxsspx		vs4,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1
	lxsspx		vs17,	o4,	T1

	addi		BO,	BO,	8


	xsmulsp		vs32,	vs0,	vs8

	xsmulsp		vs33,	vs0,	vs9


.endm

.macro KERNEL2x1_1


	lxsspx		vs4,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1
	lxsspx		vs17,	o4,	T1

	addi		BO,	BO,	8


	xsmaddasp	vs32,	vs0,	vs8

	xsmaddasp	vs33,	vs0,	vs9


.endm

.macro KERNEL2x1_2


	lxsspx		vs0,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1

	addi		BO,	BO,	8


	xsmaddasp	vs32,	vs4,	vs16

	xsmaddasp	vs33,	vs4,	vs17


.endm

.macro KERNEL2x1_E2


	xsmaddasp	vs32,	vs4,	vs16

	xsmaddasp	vs33,	vs4,	vs17


.endm

.macro KERNEL2x1_SUBI1


	lxsspx		vs0,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1

	addi		BO,	BO,	8


	xsmulsp		vs32,	vs0,	vs8

	xsmulsp		vs33,	vs0,	vs9


.endm

.macro KERNEL2x1_SUB1


	lxsspx		vs0,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1
	lxsspx		vs9,	o4,	T1

	addi		BO,	BO,	8


	xsmaddasp	vs32,	vs0,	vs8

	xsmaddasp	vs33,	vs0,	vs9


.endm

.macro SAVE2x1

	mr		T1,	CO

#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs32,	alpha_r

#else

	xsmulsp		vs28,	vs32,	alpha_r
	xsaddsp		vs0,	vs0,	vs28

#endif

	stxsspx		vs0,	o0,	T1

	add		T1,	T1,	LDC


#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs33,	alpha_r

#else

	xsmulsp		vs28,	vs33,	alpha_r
	xsaddsp		vs0,	vs0,	vs28

#endif

	stxsspx		vs0,	o0,	T1

	add		T1,	T1,	LDC

	addi		CO,	CO,	4

.endm


/**********************************************************************************************
* Macros for N=1 and M=16
**********************************************************************************************/

.macro LOAD1x16_1

	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO
	lxvw4x		vs2,	o32,	AO
	lxvw4x		vs3,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0

	addi		BO,	BO,	4

.endm

.macro KERNEL1x16_I1


	lxvw4x		vs4,	o0,	AO
	lxvw4x		vs5,	o16,	AO
	lxvw4x		vs6,	o32,	AO
	lxvw4x		vs7,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0

	addi		BO,	BO,	4


	xvmulsp		vs32,	vs0,	vs8
	xvmulsp		vs33,	vs1,	vs8
	xvmulsp		vs34,	vs2,	vs8
	xvmulsp		vs35,	vs3,	vs8


.endm

.macro KERNEL1x16_1


	lxvw4x		vs4,	o0,	AO
	lxvw4x		vs5,	o16,	AO
	lxvw4x		vs6,	o32,	AO
	lxvw4x		vs7,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0

	addi		BO,	BO,	4


	xvmaddasp	vs32,	vs0,	vs8
	xvmaddasp	vs33,	vs1,	vs8
	xvmaddasp	vs34,	vs2,	vs8
	xvmaddasp	vs35,	vs3,	vs8


.endm

.macro KERNEL1x16_2


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO
	lxvw4x		vs2,	o32,	AO
	lxvw4x		vs3,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0

	addi		BO,	BO,	4


	xvmaddasp	vs32,	vs4,	vs16
	xvmaddasp	vs33,	vs5,	vs16
	xvmaddasp	vs34,	vs6,	vs16
	xvmaddasp	vs35,	vs7,	vs16


.endm

.macro KERNEL1x16_E2


	xvmaddasp	vs32,	vs4,	vs16
	xvmaddasp	vs33,	vs5,	vs16
	xvmaddasp	vs34,	vs6,	vs16
	xvmaddasp	vs35,	vs7,	vs16


.endm

.macro KERNEL1x16_SUBI1


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO
	lxvw4x		vs2,	o32,	AO
	lxvw4x		vs3,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0

	addi		BO,	BO,	4


	xvmulsp		vs32,	vs0,	vs8
	xvmulsp		vs33,	vs1,	vs8
	xvmulsp		vs34,	vs2,	vs8
	xvmulsp		vs35,	vs3,	vs8


.endm

.macro KERNEL1x16_SUB1


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO
	lxvw4x		vs2,	o32,	AO
	lxvw4x		vs3,	o48,	AO

	addi		AO,	AO,	64

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0

	addi		BO,	BO,	4


	xvmaddasp	vs32,	vs0,	vs8
	xvmaddasp	vs33,	vs1,	vs8
	xvmaddasp	vs34,	vs2,	vs8
	xvmaddasp	vs35,	vs3,	vs8


.endm

.macro SAVE1x16

	mr		T1,	CO

#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1
	lxvw4x		vs2,	o32,	T1
	lxvw4x		vs3,	o48,	T1

#endif


	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs33,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif

	stxvw4x		vs34,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs2,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs2,	vs2,	vs28
#endif

	stxvw4x		vs35,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs3,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs3,	vs3,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1
	stxvw4x		vs2,	o32,	T1
	stxvw4x		vs3,	o48,	T1

	add		T1,	T1,	LDC

	addi		CO,	CO,	64

.endm


/**********************************************************************************************
* Macros for N=1 and M=8
**********************************************************************************************/

.macro LOAD1x8_1

	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0

	addi		BO,	BO,	4

.endm

.macro KERNEL1x8_I1


	lxvw4x		vs4,	o0,	AO
	lxvw4x		vs5,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0

	addi		BO,	BO,	4


	xvmulsp		vs32,	vs0,	vs8
	xvmulsp		vs33,	vs1,	vs8


.endm

.macro KERNEL1x8_1


	lxvw4x		vs4,	o0,	AO
	lxvw4x		vs5,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0

	addi		BO,	BO,	4


	xvmaddasp	vs32,	vs0,	vs8
	xvmaddasp	vs33,	vs1,	vs8


.endm

.macro KERNEL1x8_2


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0

	addi		BO,	BO,	4


	xvmaddasp	vs32,	vs4,	vs16
	xvmaddasp	vs33,	vs5,	vs16


.endm

.macro KERNEL1x8_E2


	xvmaddasp	vs32,	vs4,	vs16
	xvmaddasp	vs33,	vs5,	vs16


.endm

.macro KERNEL1x8_SUBI1


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0

	addi		BO,	BO,	4


	xvmulsp		vs32,	vs0,	vs8
	xvmulsp		vs33,	vs1,	vs8


.endm

.macro KERNEL1x8_SUB1


	lxvw4x		vs0,	o0,	AO
	lxvw4x		vs1,	o16,	AO

	addi		AO,	AO,	32

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0

	addi		BO,	BO,	4


	xvmaddasp	vs32,	vs0,	vs8
	xvmaddasp	vs33,	vs1,	vs8


.endm

.macro SAVE1x8

	mr		T1,	CO

#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1
	lxvw4x		vs1,	o16,	T1

#endif


	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif

	stxvw4x		vs33,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs1,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs1,	vs1,	vs28
#endif




	stxvw4x		vs0,	o0,	T1
	stxvw4x		vs1,	o16,	T1

	add		T1,	T1,	LDC

	addi		CO,	CO,	32

.endm


/**********************************************************************************************
* Macros for N=1 and M=4
**********************************************************************************************/

.macro LOAD1x4_1

	lxvw4x		vs0,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0

	addi		BO,	BO,	4

.endm

.macro KERNEL1x4_I1


	lxvw4x		vs4,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0

	addi		BO,	BO,	4


	xvmulsp		vs32,	vs0,	vs8


.endm

.macro KERNEL1x4_1


	lxvw4x		vs4,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs16,	vs28,	0

	addi		BO,	BO,	4


	xvmaddasp	vs32,	vs0,	vs8


.endm

.macro KERNEL1x4_2


	lxvw4x		vs0,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0

	addi		BO,	BO,	4


	xvmaddasp	vs32,	vs4,	vs16


.endm

.macro KERNEL1x4_E2


	xvmaddasp	vs32,	vs4,	vs16


.endm

.macro KERNEL1x4_SUBI1


	lxvw4x		vs0,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0

	addi		BO,	BO,	4


	xvmulsp		vs32,	vs0,	vs8


.endm

.macro KERNEL1x4_SUB1


	lxvw4x		vs0,	o0,	AO

	addi		AO,	AO,	16

	lxvw4x		vs28,	o0,	BO

	xxspltw		vs8,	vs28,	0

	addi		BO,	BO,	4


	xvmaddasp	vs32,	vs0,	vs8


.endm

.macro SAVE1x4

	mr		T1,	CO

#ifndef TRMMKERNEL

	lxvw4x		vs0,	o0,	T1

#endif


	stxvw4x		vs32,	o0,	TBUFFER

	lxsspx		vs4,	o0,	TBUFFER
	lxsspx		vs5,	o4,	TBUFFER
	lxsspx		vs6,	o8,	TBUFFER
	lxsspx		vs7,	o12,	TBUFFER

	xsmulsp		vs4,	vs4,	alpha_r
	xsmulsp		vs5,	vs5,	alpha_r
	xsmulsp		vs6,	vs6,	alpha_r
	xsmulsp		vs7,	vs7,	alpha_r

	stxsspx		vs4,	o0,	TBUFFER
	stxsspx		vs5,	o4,	TBUFFER
	stxsspx		vs6,	o8,	TBUFFER
	stxsspx		vs7,	o12,	TBUFFER

#ifdef TRMMKERNEL
	lxvw4x		vs0,	o0,	TBUFFER
#else
	lxvw4x		vs28,	o0,	TBUFFER
	xvaddsp		vs0,	vs0,	vs28
#endif




	stxvw4x		vs0,	o0,	T1

	add		T1,	T1,	LDC

	addi		CO,	CO,	16

.endm


/**********************************************************************************************
* Macros for N=1 and M=2
**********************************************************************************************/

.macro LOAD1x2_1

	lxsspx		vs0,	o0,	AO
	lxsspx		vs1,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1

	addi		BO,	BO,	4

.endm

.macro KERNEL1x2_I1


	lxsspx		vs4,	o0,	AO
	lxsspx		vs5,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1

	addi		BO,	BO,	4


	xsmulsp		vs32,	vs0,	vs8
	xsmulsp		vs33,	vs1,	vs8


.endm

.macro KERNEL1x2_1


	lxsspx		vs4,	o0,	AO
	lxsspx		vs5,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1

	addi		BO,	BO,	4


	xsmaddasp	vs32,	vs0,	vs8
	xsmaddasp	vs33,	vs1,	vs8


.endm

.macro KERNEL1x2_2


	lxsspx		vs0,	o0,	AO
	lxsspx		vs1,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1

	addi		BO,	BO,	4


	xsmaddasp	vs32,	vs4,	vs16
	xsmaddasp	vs33,	vs5,	vs16


.endm

.macro KERNEL1x2_E2


	xsmaddasp	vs32,	vs4,	vs16
	xsmaddasp	vs33,	vs5,	vs16


.endm

.macro KERNEL1x2_SUBI1


	lxsspx		vs0,	o0,	AO
	lxsspx		vs1,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1

	addi		BO,	BO,	4


	xsmulsp		vs32,	vs0,	vs8
	xsmulsp		vs33,	vs1,	vs8


.endm

.macro KERNEL1x2_SUB1


	lxsspx		vs0,	o0,	AO
	lxsspx		vs1,	o4,	AO

	addi		AO,	AO,	8

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1

	addi		BO,	BO,	4


	xsmaddasp	vs32,	vs0,	vs8
	xsmaddasp	vs33,	vs1,	vs8


.endm

.macro SAVE1x2

	mr		T1,	CO

#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1
	lxsspx		vs1,	o4,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs32,	alpha_r
	xsmulsp		vs1,	vs33,	alpha_r

#else

	xsmulsp		vs28,	vs32,	alpha_r
	xsaddsp		vs0,	vs0,	vs28
	xsmulsp		vs28,	vs33,	alpha_r
	xsaddsp		vs1,	vs1,	vs28

#endif

	stxsspx		vs0,	o0,	T1
	stxsspx		vs1,	o4,	T1

	add		T1,	T1,	LDC

	addi		CO,	CO,	8

.endm


/**********************************************************************************************
* Macros for N=1 and M=1
**********************************************************************************************/

.macro LOAD1x1_1

	lxsspx		vs0,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1

	addi		BO,	BO,	4

.endm

.macro KERNEL1x1_I1


	lxsspx		vs4,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1

	addi		BO,	BO,	4


	xsmulsp		vs32,	vs0,	vs8


.endm

.macro KERNEL1x1_1


	lxsspx		vs4,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs16,	o0,	T1

	addi		BO,	BO,	4


	xsmaddasp	vs32,	vs0,	vs8


.endm

.macro KERNEL1x1_2


	lxsspx		vs0,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1

	addi		BO,	BO,	4


	xsmaddasp	vs32,	vs4,	vs16


.endm

.macro KERNEL1x1_E2


	xsmaddasp	vs32,	vs4,	vs16


.endm

.macro KERNEL1x1_SUBI1


	lxsspx		vs0,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1

	addi		BO,	BO,	4


	xsmulsp		vs32,	vs0,	vs8


.endm

.macro KERNEL1x1_SUB1


	lxsspx		vs0,	o0,	AO

	addi		AO,	AO,	4

	mr		T1,	BO

	lxsspx		vs8,	o0,	T1

	addi		BO,	BO,	4


	xsmaddasp	vs32,	vs0,	vs8


.endm

.macro SAVE1x1

	mr		T1,	CO

#ifndef TRMMKERNEL

	lxsspx		vs0,	o0,	T1

#endif

#ifdef TRMMKERNEL

	xsmulsp		vs0,	vs32,	alpha_r

#else

	xsmulsp		vs28,	vs32,	alpha_r
	xsaddsp		vs0,	vs0,	vs28

#endif

	stxsspx		vs0,	o0,	T1

	add		T1,	T1,	LDC

	addi		CO,	CO,	4

.endm

