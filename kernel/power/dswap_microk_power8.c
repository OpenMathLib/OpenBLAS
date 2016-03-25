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
* 2016/03/25 Werner Saar (wernsaar@googlemail.com)
* 	 BLASTEST 		: OK
* 	 CTEST			: OK
* 	 TEST			: OK
*	 LAPACK-TEST		: OK
**************************************************************************************/

#define HAVE_KERNEL_32 1

static void dswap_kernel_32( BLASLONG n, FLOAT *x, FLOAT *y) __attribute__ ((noinline));

static void dswap_kernel_32( BLASLONG n, FLOAT *x, FLOAT *y)
{


	BLASLONG i = n;
	BLASLONG o16 = 16;
	BLASLONG o32 = 32;
	BLASLONG o48 = 48;
	BLASLONG o64 = 64;
	BLASLONG o80 = 80;
	BLASLONG o96 = 96;
	BLASLONG o112 = 112;
	FLOAT *x1=x;
	FLOAT *y1=y;
	FLOAT *x2=x+1;
	FLOAT *y2=y+1;
	BLASLONG pre = 384;
	BLASLONG alpha=0;

	__asm__  __volatile__
	(

	"addi		%3, %3, -8			    \n\t"	
	"addi		%4, %4, -8			    \n\t"	

	".align 5				            \n\t"
	"1:				                    \n\t"

	"lxvd2x		32, 0, %2			    \n\t"
	"lxvd2x		33, %5, %2			    \n\t"
	"lxvd2x		34, %6, %2			    \n\t"
	"lxvd2x		35, %7, %2			    \n\t"
	"lxvd2x		36, %8, %2			    \n\t"
	"lxvd2x		37, %9, %2			    \n\t"
	"lxvd2x		38, %10, %2			    \n\t"
	"lxvd2x		39, %11, %2			    \n\t"

	"addi		%2, %2, 128			    \n\t"

	"lxvd2x		40, 0, %2			    \n\t"
	"lxvd2x		41, %5, %2			    \n\t"
	"lxvd2x		42, %6, %2			    \n\t"
	"lxvd2x		43, %7, %2			    \n\t"
	"lxvd2x		44, %8, %2			    \n\t"
	"lxvd2x		45, %9, %2			    \n\t"
	"lxvd2x		46, %10, %2			    \n\t"
	"lxvd2x		47, %11, %2			    \n\t"

	"addi		%2, %2, 128			    \n\t"

	"lxvd2x		48, 0, %1			    \n\t"
	"lxvd2x		49, %5, %1			    \n\t"
	"lxvd2x		50, %6, %1			    \n\t"
	"lxvd2x		51, %7, %1			    \n\t"
	"lxvd2x		52, %8, %1			    \n\t"
	"lxvd2x		53, %9, %1			    \n\t"
	"lxvd2x		54, %10, %1			    \n\t"
	"lxvd2x		55, %11, %1			    \n\t"

	"addi		%1, %1, 128			    \n\t"

	"lxvd2x		56, 0, %1			    \n\t"
	"lxvd2x		57, %5, %1			    \n\t"
	"lxvd2x		58, %6, %1			    \n\t"
	"lxvd2x		59, %7, %1			    \n\t"
	"lxvd2x		60, %8, %1			    \n\t"
	"lxvd2x		61, %9, %1			    \n\t"
	"lxvd2x		62, %10, %1			    \n\t"
	"lxvd2x		63, %11, %1			    \n\t"

	"addi		%1, %1, 128			    \n\t"

	"stxvd2x		32, 0, %3			    \n\t"
	"stxvd2x		33, %5, %3			    \n\t"
	"stxvd2x		34, %6, %3			    \n\t"
	"stxvd2x		35, %7, %3			    \n\t"
	"stxvd2x		36, %8, %3			    \n\t"
	"stxvd2x		37, %9, %3			    \n\t"
	"stxvd2x		38, %10, %3			    \n\t"
	"stxvd2x		39, %11, %3			    \n\t"

	"addi		%3, %3, 128			    \n\t"

	"stxvd2x		40, 0, %3			    \n\t"
	"stxvd2x		41, %5, %3			    \n\t"
	"stxvd2x		42, %6, %3			    \n\t"
	"stxvd2x		43, %7, %3			    \n\t"
	"stxvd2x		44, %8, %3			    \n\t"
	"stxvd2x		45, %9, %3			    \n\t"
	"stxvd2x		46, %10, %3			    \n\t"
	"stxvd2x		47, %11, %3			    \n\t"

	"addi		%3, %3, 128			    \n\t"

	"stxvd2x		48, 0, %4			    \n\t"
	"stxvd2x		49, %5, %4			    \n\t"
	"stxvd2x		50, %6, %4			    \n\t"
	"stxvd2x		51, %7, %4			    \n\t"
	"stxvd2x		52, %8, %4			    \n\t"
	"stxvd2x		53, %9, %4			    \n\t"
	"stxvd2x		54, %10, %4			    \n\t"
	"stxvd2x		55, %11, %4			    \n\t"

	"addi		%4, %4, 128			    \n\t"

	"stxvd2x		56, 0, %4			    \n\t"
	"stxvd2x		57, %5, %4			    \n\t"
	"stxvd2x		58, %6, %4			    \n\t"
	"stxvd2x		59, %7, %4			    \n\t"
	"stxvd2x		60, %8, %4			    \n\t"
	"stxvd2x		61, %9, %4			    \n\t"
	"stxvd2x		62, %10, %4			    \n\t"
	"stxvd2x		63, %11, %4			    \n\t"

	"addi		%4, %4, 128			    \n\t"

	"addic.		%0 , %0	, -32  	 	             \n\t"
	"bgt		1b		             	     \n\t"

	"2:						     \n\t"

	:
        : 
          "r" (i),	// 0	
	  "r" (y1),  	// 1
          "r" (x1),     // 2
          "r" (y2),     // 3
          "r" (x2),     // 4
	  "r" (o16),	// 5
	  "r" (o32),	// 6
	  "r" (o48),    // 7
          "r" (o64),    // 8
          "r" (o80),    // 9
          "r" (o96),    // 10
          "r" (o112)    // 11
	: "cr0", "%0", "%2" , "%1", "%3", "%4", "memory"
	);

} 


