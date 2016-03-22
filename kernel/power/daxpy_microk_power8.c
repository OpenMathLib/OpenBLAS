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
* 2016/03/22 Werner Saar (wernsaar@googlemail.com)
* 	 BLASTEST 		: OK
* 	 CTEST			: OK
* 	 TEST			: OK
*	 LAPACK-TEST		: OK
**************************************************************************************/


#define HAVE_KERNEL_8 1
static void daxpy_kernel_8( BLASLONG n, FLOAT *x, FLOAT *y , FLOAT *alpha) __attribute__ ((noinline));

static void daxpy_kernel_8( BLASLONG n, FLOAT *x, FLOAT *y, FLOAT *alpha)
{


	BLASLONG i = n;
	BLASLONG o16 = 16;
	BLASLONG o32 = 32;
	BLASLONG o48 = 48;
	FLOAT *x1=x;
	FLOAT *y1=y;
	FLOAT *y2=y+1;
	BLASLONG pre = 384;

	__asm__  __volatile__
	(

	"lxsdx		33, %5, %4			    \n\t"
	"xxspltd	32, 33, 0			    \n\t"
	"addi		%8, %8, -8			    \n\t"

	"dcbt		%2, %9				    \n\t"
	"dcbt		%3, %9				    \n\t"

	"lxvd2x		40, 0, %2			    \n\t"
	"lxvd2x		41, %5, %2			    \n\t"
	"lxvd2x		42, %6, %2			    \n\t"
	"lxvd2x		43, %7, %2			    \n\t"

	"lxvd2x		48, 0, %3			    \n\t"
	"lxvd2x		49, %5, %3			    \n\t"
	"lxvd2x		50, %6, %3			    \n\t"
	"lxvd2x		51, %7, %3			    \n\t"

	"addi		%2, %2, 64			    \n\t"
	"addi		%3, %3, 64			    \n\t"

	"lxvd2x		44, 0, %2			    \n\t"
	"lxvd2x		45, %5, %2			    \n\t"
	"lxvd2x		46, %6, %2			    \n\t"
	"lxvd2x		47, %7, %2			    \n\t"

	"lxvd2x		52, 0, %3			    \n\t"
	"lxvd2x		53, %5, %3			    \n\t"
	"lxvd2x		54, %6, %3			    \n\t"
	"lxvd2x		55, %7, %3			    \n\t"

	"addi		%2, %2, 64			    \n\t"
	"addi		%3, %3, 64			    \n\t"

	"addic.		%0 , %0	, -16  	 	             \n\t"
	"ble		2f		             	     \n\t"

	".align 5				            \n\t"
	"1:				                    \n\t"

	"dcbt		%2, %9				    \n\t"
	"dcbt		%3, %9				    \n\t"

	"xvmaddadp	48, 40, 32		    	    \n\t"
	"xvmaddadp	49, 41, 32		    	    \n\t"

	"lxvd2x		40, 0, %2			    \n\t"
	"lxvd2x		41, %5, %2			    \n\t"

	"stxvd2x	48,  0, %8			    \n\t"
	"stxvd2x	49, %5, %8			    \n\t"

	"xvmaddadp	50, 42, 32		    	    \n\t"
	"xvmaddadp	51, 43, 32		    	    \n\t"

	"lxvd2x		42, %6, %2			    \n\t"
	"lxvd2x		43, %7, %2			    \n\t"

	"stxvd2x	50, %6, %8			    \n\t"
	"stxvd2x	51, %7, %8			    \n\t"

	"lxvd2x		48, 0, %3			    \n\t"
	"lxvd2x		49, %5, %3			    \n\t"
	"lxvd2x		50, %6, %3			    \n\t"
	"lxvd2x		51, %7, %3			    \n\t"

	"addi		%2, %2, 64			    \n\t"
	"addi		%8, %8, 64			    \n\t"

	"xvmaddadp	52, 44, 32		    	    \n\t"
	"addi		%3, %3, 64			    \n\t"
	"xvmaddadp	53, 45, 32		    	    \n\t"

	"lxvd2x		44, 0, %2			    \n\t"
	"lxvd2x		45, %5, %2			    \n\t"

	"stxvd2x	52,  0, %8			    \n\t"
	"stxvd2x	53, %5, %8			    \n\t"

	"xvmaddadp	54, 46, 32		    	    \n\t"
	"xvmaddadp	55, 47, 32		    	    \n\t"

	"lxvd2x		46, %6, %2			    \n\t"
	"lxvd2x		47, %7, %2			    \n\t"

	"stxvd2x	54, %6, %8			    \n\t"
	"stxvd2x	55, %7, %8			    \n\t"

	"addi		%2, %2, 64			    \n\t"
	"addi		%8, %8, 64			    \n\t"

	"lxvd2x		52, 0, %3			    \n\t"
	"lxvd2x		53, %5, %3			    \n\t"
	"lxvd2x		54, %6, %3			    \n\t"
	"lxvd2x		55, %7, %3			    \n\t"

	"addi		%3, %3, 64			    \n\t"


	"addic.		%0 , %0	, -16  	 	             \n\t"
	"bgt		1b		             	     \n\t"

	"2:						     \n\t"

	
	"xvmaddadp	48, 40, 32		    	    \n\t"
	"xvmaddadp	49, 41, 32		    	    \n\t"
	"xvmaddadp	50, 42, 32		    	    \n\t"
	"xvmaddadp	51, 43, 32		    	    \n\t"

	"xvmaddadp	52, 44, 32		    	    \n\t"
	"xvmaddadp	53, 45, 32		    	    \n\t"
	"xvmaddadp	54, 46, 32		    	    \n\t"
	"xvmaddadp	55, 47, 32		    	    \n\t"

	"stxvd2x	48,  0, %8			    \n\t"
	"stxvd2x	49, %5, %8			    \n\t"
	"stxvd2x	50, %6, %8			    \n\t"
	"stxvd2x	51, %7, %8			    \n\t"

	"addi		%8, %8, 64			    \n\t"

	"stxvd2x	52,  0, %8			    \n\t"
	"stxvd2x	53, %5, %8			    \n\t"
	"stxvd2x	54, %6, %8			    \n\t"
	"stxvd2x	55, %7, %8			    \n\t"

	"addi		%8, %8, 64			    \n\t"

	:
        : 
          "r" (i),	// 0	
	  "r" (n),  	// 1
          "r" (x1),     // 2
          "r" (y1),     // 3
          "r" (alpha),    // 4
	  "r" (o16),	// 5
	  "r" (o32),	// 6
	  "r" (o48),    // 7
	  "r" (y2),     // 8
	  "r" (pre)	// 9
	: "cr0", "%0", "%2" , "%3", "%8", "memory"
	);

} 


