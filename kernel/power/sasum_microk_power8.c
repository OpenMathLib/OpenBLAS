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
* 2016/03/28 Werner Saar (wernsaar@googlemail.com)
* 	 BLASTEST 		: OK
* 	 CTEST			: OK
* 	 TEST			: OK
*	 LAPACK-TEST		: OK
**************************************************************************************/

#define HAVE_KERNEL_32 1
static void sasum_kernel_32( BLASLONG n, FLOAT *x, FLOAT *svec) __attribute__ ((noinline));

static void sasum_kernel_32( BLASLONG n, FLOAT *x, FLOAT *svec)
{
	BLASLONG o16 = 16;
	BLASLONG o32 = 32;
	BLASLONG o48 = 48;
	BLASLONG o64 = 64;
	BLASLONG o80 = 80;
	BLASLONG o96 = 96;
	BLASLONG o112 = 112;
	BLASLONG pre = 384;

	__asm__
	(
	"dcbt		%1, %3			\n\t"

	"xxlxor		32, 32,	32		\n\t"
	"xxlxor		33, 33,	33		\n\t"
	"xxlxor		34, 34,	34		\n\t"
	"xxlxor		35, 35,	35		\n\t"
	"xxlxor		36, 36,	36		\n\t"
	"xxlxor		37, 37,	37		\n\t"
	"xxlxor		38, 38,	38		\n\t"
	"xxlxor		39, 39,	39		\n\t"

	"lxvw4x		40, 0, %1		\n\t"
	"lxvw4x		41, %4, %1		\n\t"
	"lxvw4x		42, %5, %1		\n\t"
	"lxvw4x		43, %6, %1		\n\t"
	"lxvw4x		44, %7, %1		\n\t"
	"lxvw4x		45, %8, %1		\n\t"
	"lxvw4x		46, %9, %1		\n\t"
	"lxvw4x		47, %10, %1		\n\t"

	"addi		%1, %1, 128		\n\t"
	"addic.		%2, %2, -32		\n\t"
	"ble		2f			\n\t"

	".p2align 5				\n\t"
	"1:					\n\t"
	"dcbt		%1, %3			\n\t"

	"xvabssp	48, 40			\n\t"
	"xvabssp	49, 41			\n\t"
	"xvabssp	50, 42			\n\t"
	"xvabssp	51, 43			\n\t"

	"lxvw4x		40, 0, %1		\n\t"
	"lxvw4x		41, %4, %1		\n\t"

	"xvabssp	52, 44			\n\t"
	"xvabssp	53, 45			\n\t"

	"lxvw4x		42, %5, %1		\n\t"
	"lxvw4x		43, %6, %1		\n\t"

	"xvabssp	54, 46			\n\t"
	"xvabssp	55, 47			\n\t"

	"lxvw4x		44, %7, %1		\n\t"
	"lxvw4x		45, %8, %1		\n\t"

	"xvaddsp	32, 32, 48		\n\t"
	"xvaddsp	33, 33, 49		\n\t"

	"lxvw4x		46, %9, %1		\n\t"
	"lxvw4x		47, %10, %1		\n\t"

	"xvaddsp	34, 34, 50		\n\t"
	"xvaddsp	35, 35, 51		\n\t"
	"addi		%1, %1, 128		\n\t"
	"xvaddsp	36, 36, 52		\n\t"
	"xvaddsp	37, 37, 53		\n\t"
	"addic.		%2, %2, -32		\n\t"
	"xvaddsp	38, 38, 54		\n\t"
	"xvaddsp	39, 39, 55		\n\t"

	"bgt		1b			\n\t"

	"2:					\n\t"
	"xvabssp	48, 40			\n\t"
	"xvabssp	49, 41			\n\t"
	"xvabssp	50, 42			\n\t"
	"xvabssp	51, 43			\n\t"
	"xvabssp	52, 44			\n\t"
	"xvabssp	53, 45			\n\t"
	"xvabssp	54, 46			\n\t"
	"xvabssp	55, 47			\n\t"

	"xvaddsp	32, 32, 48		\n\t"
	"xvaddsp	33, 33, 49		\n\t"
	"xvaddsp	34, 34, 50		\n\t"
	"xvaddsp	35, 35, 51		\n\t"
	"xvaddsp	36, 36, 52		\n\t"
	"xvaddsp	37, 37, 53		\n\t"
	"xvaddsp	38, 38, 54		\n\t"
	"xvaddsp	39, 39, 55		\n\t"

	"xvaddsp	32, 32, 33		\n\t"
	"xvaddsp	34, 34, 35		\n\t"
	"xvaddsp	36, 36, 37		\n\t"
	"xvaddsp	38, 38, 39		\n\t"

	"xvaddsp	32, 32, 34		\n\t"
	"xvaddsp	36, 36, 38		\n\t"

	"xvaddsp	32, 32, 36		\n\t"

	"stxvw4x	32, %y0			\n\t"

	:
	  "=m" (*svec),	// 0
	  "+b" (x),	// 1
	  "+r" (n)	// 2
	:
	  "r" (pre),	// 3
	  "r" (o16),	// 4
	  "r" (o32),	// 5
	  "r" (o48),	// 6
	  "r" (o64),	// 7
	  "r" (o80),	// 8
	  "r" (o96),	// 9
	  "r" (o112)	// 10
	:
	  "cr0","32","33","34","35","36","37","38","39",
	  "40","41","42","43","44","45","46","47",
	  "48","49","50","51","52","53","54","55"
	);
} 


