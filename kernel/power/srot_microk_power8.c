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
* 2016/03/27 Werner Saar (wernsaar@googlemail.com)
*
* I don't use fused multiply-add ( precision problems with lapack )
*
* 	 BLASTEST 		: OK
* 	 CTEST			: OK
* 	 TEST			: OK
*	 LAPACK-TEST		: OK
**************************************************************************************/

#define HAVE_KERNEL_16 1

static void srot_kernel_16( BLASLONG n, FLOAT *x, FLOAT *y, FLOAT *c, FLOAT *s) __attribute__ ((noinline));

static void srot_kernel_16( BLASLONG n, FLOAT *x, FLOAT *y, FLOAT *c, FLOAT *s)
{


	BLASLONG i = n;
	BLASLONG o16 = 16;
	BLASLONG o32 = 32;
	BLASLONG o48 = 48;
	FLOAT *x1=x;
	FLOAT *y1=y;
	FLOAT *x2=x+1;
	FLOAT *y2=y+1;

	__asm__  __volatile__
	(

        "lxvw4x         36 , 0, %3                          \n\t"	// load c
        "lxvw4x         37 , 0, %4                          \n\t"	// load s
	"addi		%8 , %8, -4			     \n\t"
	"addi		%9 , %9, -4			     \n\t"

	"lxvw4x		32, 0, %1			    \n\t"	// load x
	"lxvw4x		33, %5, %1			    \n\t"
	"lxvw4x		34, %6, %1			    \n\t"
	"lxvw4x		35, %7, %1			    \n\t"

	"lxvw4x		40, 0, %2			    \n\t"	// load y
	"lxvw4x		41, %5, %2			    \n\t"
	"lxvw4x		42, %6, %2			    \n\t"
	"lxvw4x		43, %7, %2			    \n\t"

	"addi		%1, %1, 64			    \n\t"
	"addi		%2, %2, 64			    \n\t"

	"addic.		%0 , %0	, -16  	 	             \n\t"
	"ble		2f		             	     \n\t"

	".align 5				            \n\t"
	"1:				                    \n\t"

	"xvmulsp	48, 32, 36		    	    \n\t"	// c * x
	"xvmulsp	49, 33, 36		    	    \n\t"
	"xvmulsp	50, 34, 36		    	    \n\t"
	"xvmulsp	51, 35, 36		    	    \n\t"

	"xvmulsp	56, 40, 36		    	    \n\t"	// c * y
	"xvmulsp	57, 41, 36		    	    \n\t"
	"xvmulsp	58, 42, 36		    	    \n\t"
	"xvmulsp	59, 43, 36		    	    \n\t"

	"xvmulsp	52, 32, 37		    	    \n\t"	// s * x
	"xvmulsp	53, 33, 37		    	    \n\t"

	"lxvw4x		32, 0, %1			    \n\t"	// load x
	"lxvw4x		33, %5, %1			    \n\t"

	"xvmulsp	54, 34, 37		    	    \n\t"
	"xvmulsp	55, 35, 37		    	    \n\t"

	"lxvw4x		34, %6, %1			    \n\t"
	"lxvw4x		35, %7, %1			    \n\t"

	"xvmulsp	60, 40, 37		    	    \n\t"	// s * y
	"xvmulsp	61, 41, 37		    	    \n\t"

	"lxvw4x		40, 0, %2			    \n\t"	// load y
	"lxvw4x		41, %5, %2			    \n\t"

	"xvmulsp	62, 42, 37		    	    \n\t"
	"xvmulsp	63, 43, 37		    	    \n\t"

	"lxvw4x		42, %6, %2			    \n\t"
	"lxvw4x		43, %7, %2			    \n\t"

	"xvaddsp	48, 48 , 60			    \n\t"	// c * x + s * y 
	"xvaddsp	49, 49 , 61			    \n\t"	// c * x + s * y 

	"addi		%1, %1, 64			    \n\t"
	"addi		%2, %2, 64			    \n\t"

	"xvaddsp	50, 50 , 62			    \n\t"	// c * x + s * y 
	"xvaddsp	51, 51 , 63			    \n\t"	// c * x + s * y 

	"xvsubsp	56, 56 , 52			    \n\t"	// c * y - s * x
	"xvsubsp	57, 57 , 53			    \n\t"	// c * y - s * x
	"xvsubsp	58, 58 , 54			    \n\t"	// c * y - s * x
	"xvsubsp	59, 59 , 55			    \n\t"	// c * y - s * x

	"stxvw4x	48, 0, %8			    \n\t"	// store x
	"stxvw4x	49, %5, %8			    \n\t"
	"stxvw4x	50, %6, %8			    \n\t"
	"stxvw4x	51, %7, %8			    \n\t"

	"stxvw4x	56, 0, %9			    \n\t"	// store y
	"stxvw4x	57, %5, %9			    \n\t"
	"stxvw4x	58, %6, %9			    \n\t"
	"stxvw4x	59, %7, %9			    \n\t"

	"addi		%8, %8, 64			    \n\t"
	"addi		%9, %9, 64			    \n\t"

	"addic.		%0 , %0	, -16  	 	             \n\t"
	"bgt		1b		             	     \n\t"

	"2:						     \n\t"

	"xvmulsp	48, 32, 36		    	    \n\t"	// c * x
	"xvmulsp	49, 33, 36		    	    \n\t"
	"xvmulsp	50, 34, 36		    	    \n\t"
	"xvmulsp	51, 35, 36		    	    \n\t"

	"xvmulsp	56, 40, 36		    	    \n\t"	// c * y
	"xvmulsp	57, 41, 36		    	    \n\t"
	"xvmulsp	58, 42, 36		    	    \n\t"
	"xvmulsp	59, 43, 36		    	    \n\t"

	"xvmulsp	52, 32, 37		    	    \n\t"	// s * x
	"xvmulsp	53, 33, 37		    	    \n\t"
	"xvmulsp	54, 34, 37		    	    \n\t"
	"xvmulsp	55, 35, 37		    	    \n\t"

	"xvmulsp	60, 40, 37		    	    \n\t"	// s * y
	"xvmulsp	61, 41, 37		    	    \n\t"
	"xvmulsp	62, 42, 37		    	    \n\t"
	"xvmulsp	63, 43, 37		    	    \n\t"

	"xvaddsp	48, 48 , 60			    \n\t"	// c * x + s * y 
	"xvaddsp	49, 49 , 61			    \n\t"	// c * x + s * y 
	"xvaddsp	50, 50 , 62			    \n\t"	// c * x + s * y 
	"xvaddsp	51, 51 , 63			    \n\t"	// c * x + s * y 

	"xvsubsp	56, 56 , 52			    \n\t"	// c * y - s * x
	"xvsubsp	57, 57 , 53			    \n\t"	// c * y - s * x
	"xvsubsp	58, 58 , 54			    \n\t"	// c * y - s * x
	"xvsubsp	59, 59 , 55			    \n\t"	// c * y - s * x

	"stxvw4x	48, 0, %8			    \n\t"	// store x
	"stxvw4x	49, %5, %8			    \n\t"
	"stxvw4x	50, %6, %8			    \n\t"
	"stxvw4x	51, %7, %8			    \n\t"

	"stxvw4x	56, 0, %9			    \n\t"	// store y
	"stxvw4x	57, %5, %9			    \n\t"
	"stxvw4x	58, %6, %9			    \n\t"
	"stxvw4x	59, %7, %9			    \n\t"



	:
        : 
          "r" (i),	// 0	
	  "r" (x1),  	// 1
          "r" (y1),     // 2
          "r" (c),      // 3
          "r" (s),      // 4
	  "r" (o16),	// 5
	  "r" (o32),	// 6
	  "r" (o48),    // 7
	  "r" (x2),     // 8
	  "r" (y2)      // 9
	: "cr0", "%0", "%1" , "%2", "%8", "%9", "memory"
	);

} 


