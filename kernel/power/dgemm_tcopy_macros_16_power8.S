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


/**********************************************************************************************
* Macros for N=4 and M=16
**********************************************************************************************/

#ifndef POWER10
#if defined(_AIX)
define(`COPY_4x16', `
#else
.macro COPY_4x16
#endif

	lxvd2x		vs32,	o0,	A0
	lxvd2x		vs33,	o16,	A0
	lxvd2x		vs34,	o32,	A0
	lxvd2x		vs35,	o48,	A0
	addi		A0,	A0,	64

	lxvd2x		vs40,	o0,	A1
	lxvd2x		vs41,	o16,	A1
	lxvd2x		vs42,	o32,	A1
	lxvd2x		vs43,	o48,	A1
	addi		A1,	A1,	64

	lxvd2x		vs48,	o0,	A2
	lxvd2x		vs49,	o16,	A2
	lxvd2x		vs50,	o32,	A2
	lxvd2x		vs51,	o48,	A2
	addi		A2,	A2,	64

	lxvd2x		vs4,	o0,	A3
	lxvd2x		vs5,	o16,	A3
	lxvd2x		vs6,	o32,	A3
	lxvd2x		vs7,	o48,	A3
	addi		A3,	A3,	64

	lxvd2x		vs36,	o0,	A0
	lxvd2x		vs37,	o16,	A0
	lxvd2x		vs38,	o32,	A0
	lxvd2x		vs39,	o48,	A0
	addi		A0,	A0,	64

	lxvd2x		vs44,	o0,	A1
	lxvd2x		vs45,	o16,	A1
	lxvd2x		vs46,	o32,	A1
	lxvd2x		vs47,	o48,	A1
	addi		A1,	A1,	64

	lxvd2x		vs12,	o0,	A2
	lxvd2x		vs13,	o16,	A2
	lxvd2x		vs2,	o32,	A2
	lxvd2x		vs3,	o48,	A2
	addi		A2,	A2,	64

	lxvd2x		vs8,	o0,	A3
	lxvd2x		vs9,	o16,	A3
	lxvd2x		vs10,	o32,	A3
	lxvd2x		vs11,	o48,	A3
	addi		A3,	A3,	64

	mr		T1,	BO

	stxvd2x		vs32,	o0,	T1
	stxvd2x		vs33,	o16,	T1
	stxvd2x		vs34,	o32,	T1
	stxvd2x		vs35,	o48,	T1
	addi		T1,	T1,	64

	stxvd2x		vs36,	o0,	T1
	stxvd2x		vs37,	o16,	T1
	stxvd2x		vs38,	o32,	T1
	stxvd2x		vs39,	o48,	T1

	addi		T1,	T1,	64

	stxvd2x		vs40,	o0,	T1
	stxvd2x		vs41,	o16,	T1
	stxvd2x		vs42,	o32,	T1
	stxvd2x		vs43,	o48,	T1
	addi		T1,	T1,	64

	stxvd2x		vs44,	o0,	T1
	stxvd2x		vs45,	o16,	T1
	stxvd2x		vs46,	o32,	T1
	stxvd2x		vs47,	o48,	T1

	addi		T1,	T1,	64

	stxvd2x		vs48,	o0,	T1
	stxvd2x		vs49,	o16,	T1
	stxvd2x		vs50,	o32,	T1
	stxvd2x		vs51,	o48,	T1
	addi		T1,	T1,	64

	stxvd2x		vs12,	o0,	T1
	stxvd2x		vs13,	o16,	T1
	stxvd2x		vs2,	o32,	T1
	stxvd2x		vs3,	o48,	T1

	addi		T1,	T1,	64

	stxvd2x		vs4,	o0,	T1
	stxvd2x		vs5,	o16,	T1
	stxvd2x		vs6,	o32,	T1
	stxvd2x		vs7,	o48,	T1
	addi		T1,	T1,	64

	stxvd2x		vs8,	o0,	T1
	stxvd2x		vs9,	o16,	T1
	stxvd2x		vs10,	o32,	T1
	stxvd2x		vs11,	o48,	T1

#if defined(_AIX)
')
#else
.endm
#endif


/**********************************************************************************************
* Macros for N=4 and M=8
**********************************************************************************************/

#if defined(_AIX)
define(`COPY_4x8', `
#else
.macro COPY_4x8
#endif

	lxvd2x		vs32,	o0,	A0
	lxvd2x		vs33,	o16,	A0
	lxvd2x		vs34,	o32,	A0
	lxvd2x		vs35,	o48,	A0
	addi		A0,	A0,	64


	lxvd2x		vs36,	o0,	A1
	lxvd2x		vs37,	o16,	A1
	lxvd2x		vs38,	o32,	A1
	lxvd2x		vs39,	o48,	A1
	addi		A1,	A1,	64


	lxvd2x		vs40,	o0,	A2
	lxvd2x		vs41,	o16,	A2
	lxvd2x		vs42,	o32,	A2
	lxvd2x		vs43,	o48,	A2
	addi		A2,	A2,	64


	lxvd2x		vs44,	o0,	A3
	lxvd2x		vs45,	o16,	A3
	lxvd2x		vs46,	o32,	A3
	lxvd2x		vs47,	o48,	A3
	addi		A3,	A3,	64


	mr		T1,	BO

	stxvd2x		vs32,	o0,	T1
	stxvd2x		vs33,	o16,	T1
	stxvd2x		vs34,	o32,	T1
	stxvd2x		vs35,	o48,	T1

	addi		T1,	T1,	64

	stxvd2x		vs36,	o0,	T1
	stxvd2x		vs37,	o16,	T1
	stxvd2x		vs38,	o32,	T1
	stxvd2x		vs39,	o48,	T1

	addi		T1,	T1,	64

	stxvd2x		vs40,	o0,	T1
	stxvd2x		vs41,	o16,	T1
	stxvd2x		vs42,	o32,	T1
	stxvd2x		vs43,	o48,	T1

	addi		T1,	T1,	64

	stxvd2x		vs44,	o0,	T1
	stxvd2x		vs45,	o16,	T1
	stxvd2x		vs46,	o32,	T1
	stxvd2x		vs47,	o48,	T1

#if defined(_AIX)
')
#else
.endm
#endif


/**********************************************************************************************
* Macros for N=4 and M=4
**********************************************************************************************/

#if defined(_AIX)
define(`COPY_4x4', `
#else
.macro COPY_4x4
#endif

	lxvd2x		vs32,	o0,	A0
	lxvd2x		vs33,	o16,	A0
	addi		A0,	A0,	32


	lxvd2x		vs34,	o0,	A1
	lxvd2x		vs35,	o16,	A1
	addi		A1,	A1,	32


	lxvd2x		vs36,	o0,	A2
	lxvd2x		vs37,	o16,	A2
	addi		A2,	A2,	32


	lxvd2x		vs38,	o0,	A3
	lxvd2x		vs39,	o16,	A3
	addi		A3,	A3,	32


	mr		T1,	BO

	stxvd2x		vs32,	o0,	T1
	stxvd2x		vs33,	o16,	T1

	stxvd2x		vs34,	o32,	T1
	stxvd2x		vs35,	o48,	T1

	addi		T1,	T1,	64

	stxvd2x		vs36,	o0,	T1
	stxvd2x		vs37,	o16,	T1

	stxvd2x		vs38,	o32,	T1
	stxvd2x		vs39,	o48,	T1

#if defined(_AIX)
')
#else
.endm
#endif
#endif


/**********************************************************************************************
* Macros for N=4 and M=2
**********************************************************************************************/

#if defined(_AIX)
define(`COPY_4x2', `
#else
.macro COPY_4x2
#endif

	lxvd2x		vs32,	o0,	A0
	addi		A0,	A0,	16


	lxvd2x		vs33,	o0,	A1
	addi		A1,	A1,	16


	lxvd2x		vs34,	o0,	A2
	addi		A2,	A2,	16


	lxvd2x		vs35,	o0,	A3
	addi		A3,	A3,	16


	mr		T1,	BO

	stxvd2x		vs32,	o0,	T1

	stxvd2x		vs33,	o16,	T1

	stxvd2x		vs34,	o32,	T1

	stxvd2x		vs35,	o48,	T1

#if defined(_AIX)
')
#else
.endm
#endif


/**********************************************************************************************
* Macros for N=4 and M=1
**********************************************************************************************/

#if defined(_AIX)
define(`COPY_4x1', `
#else
.macro COPY_4x1
#endif

	lxsdx		vs32,	o0,	A0
	addi		A0,	A0,	8


	lxsdx		vs33,	o0,	A1
	addi		A1,	A1,	8


	lxsdx		vs34,	o0,	A2
	addi		A2,	A2,	8


	lxsdx		vs35,	o0,	A3
	addi		A3,	A3,	8


	mr		T1,	BO

	stxsdx		vs32,	o0,	T1

	stxsdx		vs33,	o8,	T1

	addi		T1,	T1,	16

	stxsdx		vs34,	o0,	T1

	stxsdx		vs35,	o8,	T1

#if defined(_AIX)
')
#else
.endm
#endif


/**********************************************************************************************
* Macros for N=2 and M=16
**********************************************************************************************/

#ifndef POWER10
#if defined(_AIX)
define(`COPY_2x16', `
#else
.macro COPY_2x16
#endif

	lxvd2x		vs32,	o0,	A0
	lxvd2x		vs33,	o16,	A0
	lxvd2x		vs34,	o32,	A0
	lxvd2x		vs35,	o48,	A0
	addi		A0,	A0,	64

	lxvd2x		vs36,	o0,	A0
	lxvd2x		vs37,	o16,	A0
	lxvd2x		vs38,	o32,	A0
	lxvd2x		vs39,	o48,	A0
	addi		A0,	A0,	64


	lxvd2x		vs40,	o0,	A1
	lxvd2x		vs41,	o16,	A1
	lxvd2x		vs42,	o32,	A1
	lxvd2x		vs43,	o48,	A1
	addi		A1,	A1,	64

	lxvd2x		vs44,	o0,	A1
	lxvd2x		vs45,	o16,	A1
	lxvd2x		vs46,	o32,	A1
	lxvd2x		vs47,	o48,	A1
	addi		A1,	A1,	64


	mr		T1,	BO

	stxvd2x		vs32,	o0,	T1
	stxvd2x		vs33,	o16,	T1
	stxvd2x		vs34,	o32,	T1
	stxvd2x		vs35,	o48,	T1
	addi		T1,	T1,	64

	stxvd2x		vs36,	o0,	T1
	stxvd2x		vs37,	o16,	T1
	stxvd2x		vs38,	o32,	T1
	stxvd2x		vs39,	o48,	T1

	addi		T1,	T1,	64

	stxvd2x		vs40,	o0,	T1
	stxvd2x		vs41,	o16,	T1
	stxvd2x		vs42,	o32,	T1
	stxvd2x		vs43,	o48,	T1
	addi		T1,	T1,	64

	stxvd2x		vs44,	o0,	T1
	stxvd2x		vs45,	o16,	T1
	stxvd2x		vs46,	o32,	T1
	stxvd2x		vs47,	o48,	T1

#if defined(_AIX)
')
#else
.endm
#endif


/**********************************************************************************************
* Macros for N=2 and M=8
**********************************************************************************************/

#if defined(_AIX)
define(`COPY_2x8', `
#else
.macro COPY_2x8
#endif

	lxvd2x		vs32,	o0,	A0
	lxvd2x		vs33,	o16,	A0
	lxvd2x		vs34,	o32,	A0
	lxvd2x		vs35,	o48,	A0
	addi		A0,	A0,	64


	lxvd2x		vs36,	o0,	A1
	lxvd2x		vs37,	o16,	A1
	lxvd2x		vs38,	o32,	A1
	lxvd2x		vs39,	o48,	A1
	addi		A1,	A1,	64


	mr		T1,	BO

	stxvd2x		vs32,	o0,	T1
	stxvd2x		vs33,	o16,	T1
	stxvd2x		vs34,	o32,	T1
	stxvd2x		vs35,	o48,	T1

	addi		T1,	T1,	64

	stxvd2x		vs36,	o0,	T1
	stxvd2x		vs37,	o16,	T1
	stxvd2x		vs38,	o32,	T1
	stxvd2x		vs39,	o48,	T1

#if defined(_AIX)
')
#else
.endm
#endif


/**********************************************************************************************
* Macros for N=2 and M=4
**********************************************************************************************/

#if defined(_AIX)
define(`COPY_2x4', `
#else
.macro COPY_2x4
#endif

	lxvd2x		vs32,	o0,	A0
	lxvd2x		vs33,	o16,	A0
	addi		A0,	A0,	32


	lxvd2x		vs34,	o0,	A1
	lxvd2x		vs35,	o16,	A1
	addi		A1,	A1,	32


	mr		T1,	BO

	stxvd2x		vs32,	o0,	T1
	stxvd2x		vs33,	o16,	T1

	stxvd2x		vs34,	o32,	T1
	stxvd2x		vs35,	o48,	T1

#if defined(_AIX)
')
#else
.endm
#endif
#endif


/**********************************************************************************************
* Macros for N=2 and M=2
**********************************************************************************************/

#if defined(_AIX)
define(`COPY_2x2', `
#else
.macro COPY_2x2
#endif

	lxvd2x		vs32,	o0,	A0
	addi		A0,	A0,	16


	lxvd2x		vs33,	o0,	A1
	addi		A1,	A1,	16


	mr		T1,	BO

	stxvd2x		vs32,	o0,	T1

	stxvd2x		vs33,	o16,	T1

#if defined(_AIX)
')
#else
.endm
#endif


/**********************************************************************************************
* Macros for N=2 and M=1
**********************************************************************************************/

#if defined(_AIX)
define(`COPY_2x1', `
#else
.macro COPY_2x1
#endif

	lxsdx		vs32,	o0,	A0
	addi		A0,	A0,	8


	lxsdx		vs33,	o0,	A1
	addi		A1,	A1,	8


	mr		T1,	BO

	stxsdx		vs32,	o0,	T1

	stxsdx		vs33,	o8,	T1

#if defined(_AIX)
')
#else
.endm
#endif


/**********************************************************************************************
* Macros for N=1 and M=16
**********************************************************************************************/

#ifndef POWER10
#if defined(_AIX)
define(`COPY_1x16', `
#else
.macro COPY_1x16
#endif

	lxvd2x		vs32,	o0,	A0
	lxvd2x		vs33,	o16,	A0
	lxvd2x		vs34,	o32,	A0
	lxvd2x		vs35,	o48,	A0
	addi		A0,	A0,	64

	lxvd2x		vs36,	o0,	A0
	lxvd2x		vs37,	o16,	A0
	lxvd2x		vs38,	o32,	A0
	lxvd2x		vs39,	o48,	A0
	addi		A0,	A0,	64


	mr		T1,	BO

	stxvd2x		vs32,	o0,	T1
	stxvd2x		vs33,	o16,	T1
	stxvd2x		vs34,	o32,	T1
	stxvd2x		vs35,	o48,	T1
	addi		T1,	T1,	64

	stxvd2x		vs36,	o0,	T1
	stxvd2x		vs37,	o16,	T1
	stxvd2x		vs38,	o32,	T1
	stxvd2x		vs39,	o48,	T1

#if defined(_AIX)
')
#else
.endm
#endif


/**********************************************************************************************
* Macros for N=1 and M=8
**********************************************************************************************/

#if defined(_AIX)
define(`COPY_1x8', `
#else
.macro COPY_1x8
#endif

	lxvd2x		vs32,	o0,	A0
	lxvd2x		vs33,	o16,	A0
	lxvd2x		vs34,	o32,	A0
	lxvd2x		vs35,	o48,	A0
	addi		A0,	A0,	64


	mr		T1,	BO

	stxvd2x		vs32,	o0,	T1
	stxvd2x		vs33,	o16,	T1
	stxvd2x		vs34,	o32,	T1
	stxvd2x		vs35,	o48,	T1

#if defined(_AIX)
')
#else
.endm
#endif


/**********************************************************************************************
* Macros for N=1 and M=4
**********************************************************************************************/

#if defined(_AIX)
define(`COPY_1x4', `
#else
.macro COPY_1x4
#endif

	lxvd2x		vs32,	o0,	A0
	lxvd2x		vs33,	o16,	A0
	addi		A0,	A0,	32


	mr		T1,	BO

	stxvd2x		vs32,	o0,	T1
	stxvd2x		vs33,	o16,	T1

#if defined(_AIX)
')
#else
.endm
#endif
#endif


/**********************************************************************************************
* Macros for N=1 and M=2
**********************************************************************************************/

#if defined(_AIX)
define(`COPY_1x2', `
#else
.macro COPY_1x2
#endif

	lxvd2x		vs32,	o0,	A0
	addi		A0,	A0,	16


	mr		T1,	BO

	stxvd2x		vs32,	o0,	T1

#if defined(_AIX)
')
#else
.endm
#endif


/**********************************************************************************************
* Macros for N=1 and M=1
**********************************************************************************************/

#if defined(_AIX)
define(`COPY_1x1', `
#else
.macro COPY_1x1
#endif

	lxsdx		vs32,	o0,	A0
	addi		A0,	A0,	8


	mr		T1,	BO

	stxsdx		vs32,	o0,	T1

#if defined(_AIX)
')
#else
.endm
#endif

