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
* 2016/03/21 Werner Saar (wernsaar@googlemail.com)
* 	 BLASTEST 		: OK
* 	 CTEST			: OK
* 	 TEST			: OK
*	 LAPACK-TEST		: OK
**************************************************************************************/

#define HAVE_KERNEL_8 1
static void zdot_kernel_8( BLASLONG n, FLOAT *x, FLOAT *y , FLOAT *dot) __attribute__ ((noinline));

static void zdot_kernel_8( BLASLONG n, FLOAT *x, FLOAT *y, FLOAT *dot)
{


	BLASLONG i = n;
	BLASLONG o16 = 16;
	BLASLONG o32 = 32;
	BLASLONG o48 = 48;
	FLOAT *x1=x;
	FLOAT *y1=y;
	BLASLONG pre = 384;

	__asm__  __volatile__
	(
	"xxlxor		32,32,32			    \n\t"
	"xxlxor		33,33,33			    \n\t"
	"xxlxor		34,34,34			    \n\t"
	"xxlxor		35,35,35			    \n\t"
	"xxlxor		36,36,36			    \n\t"
	"xxlxor		37,37,37			    \n\t"
	"xxlxor		38,38,38			    \n\t"
	"xxlxor		39,39,39			    \n\t"

	"dcbt		%2, %8				    \n\t"
	"dcbt		%3, %8				    \n\t"

	"lxvd2x		40, 0, %2			    \n\t"	// x0_r, x0_i
	"lxvd2x		48, 0, %3			    \n\t"	// y0_r, y0_i
	"lxvd2x		41, %5, %2			    \n\t"	// x1_r, x1_i
	"lxvd2x		49, %5, %3			    \n\t"	// y1_r, y1_i
	"lxvd2x		42, %6, %2			    \n\t"	// x2_r, x2_i
	"lxvd2x		50, %6, %3			    \n\t"	// y2_r, y2_i
	"lxvd2x		43, %7, %2			    \n\t"	// x3_r, x3_i
	"lxvd2x		51, %7, %3			    \n\t"	// y3_r, y3_i

	"xxswapd	52,48				    \n\t"	// y0_i, y0_r
	"xxswapd	53,49				    \n\t"	// y1_i, y1_r
	"xxswapd	54,50				    \n\t"	// y2_i, y2_r
	"xxswapd	55,51				    \n\t"	// y3_i, y3_r

	"addi		%2, %2, 64			    \n\t"
	"addi		%3, %3, 64			    \n\t"


	"lxvd2x		44, 0, %2			    \n\t"	// x0_r, x0_i
	"lxvd2x		56, 0, %3			    \n\t"	// y0_r, y0_i
	"lxvd2x		45, %5, %2			    \n\t"	// x1_r, x1_i
	"lxvd2x		57, %5, %3			    \n\t"	// y1_r, y1_i
	"lxvd2x		46, %6, %2			    \n\t"	// x2_r, x2_i
	"lxvd2x		58, %6, %3			    \n\t"	// y2_r, y2_i
	"lxvd2x		47, %7, %2			    \n\t"	// x3_r, x3_i
	"lxvd2x		59, %7, %3			    \n\t"	// y3_r, y3_i

	"xxswapd	60,56				    \n\t"	// y0_i, y0_r
	"xxswapd	61,57				    \n\t"	// y1_i, y1_r
	"xxswapd	62,58				    \n\t"	// y2_i, y2_r
	"xxswapd	63,59				    \n\t"	// y3_i, y3_r

	"addi		%2, %2, 64			    \n\t"
	"addi		%3, %3, 64			    \n\t"

	"addic.		%0 , %0	, -8  	 	             \n\t"
	"ble		2f		             	     \n\t"

	".align 5				            \n\t"
	"1:				                    \n\t"

	"dcbt		%2, %8				    \n\t"
	"dcbt		%3, %8				    \n\t"

	"xvmaddadp	32, 40, 48		    	    \n\t"	// x0_r * y0_r , x0_i * y0_i
	"lxvd2x		48, 0, %3			    \n\t"	// y0_r, y0_i
	"xvmaddadp	34, 41, 49		    	    \n\t"	// x1_r * y1_r , x1_i * y1_i
	"lxvd2x		49, %5, %3			    \n\t"	// y1_r, y1_i

	"xvmaddadp	36, 42, 50		    	    \n\t"	// x2_r * y2_r , x2_i * y2_i
	"lxvd2x		50, %6, %3			    \n\t"	// y2_r, y2_i
	"xvmaddadp	38, 43, 51		    	    \n\t"	// x3_r * y3_r , x3_i * y3_i
	"lxvd2x		51, %7, %3			    \n\t"	// y3_r, y3_i

	"xvmaddadp	33, 40, 52		    	    \n\t"	// x0_r * y0_i , x0_i * y0_r
	"lxvd2x		40, 0, %2			    \n\t"	// x0_r, x0_i
	"xvmaddadp	35, 41, 53		    	    \n\t"	// x1_r * y1_i , x1_i * y1_r
	"lxvd2x		41, %5, %2			    \n\t"	// x1_r, x1_i

	"xvmaddadp	37, 42, 54		    	    \n\t"	// x2_r * y2_i , x2_i * y2_r
	"lxvd2x		42, %6, %2			    \n\t"	// x2_r, x2_i
	"xvmaddadp	39, 43, 55		    	    \n\t"	// x3_r * y3_i , x3_i * y3_r
	"lxvd2x		43, %7, %2			    \n\t"	// x3_r, x3_i

	"xxswapd	52,48				    \n\t"	// y0_i, y0_r
	"xxswapd	53,49				    \n\t"	// y1_i, y1_r

	"addi		%2, %2, 64			    \n\t"
	"addi		%3, %3, 64			    \n\t"

	"xxswapd	54,50				    \n\t"	// y2_i, y2_r
	"xxswapd	55,51				    \n\t"	// y3_i, y3_r

	"xvmaddadp	32, 44, 56		    	    \n\t"	// x0_r * y0_r , x0_i * y0_i
	"lxvd2x		56, 0, %3			    \n\t"	// y0_r, y0_i
	"xvmaddadp	34, 45, 57		    	    \n\t"	// x1_r * y1_r , x1_i * y1_i
	"lxvd2x		57, %5, %3			    \n\t"	// y1_r, y1_i
	"xvmaddadp	36, 46, 58		    	    \n\t"	// x2_r * y2_r , x2_i * y2_i
	"lxvd2x		58, %6, %3			    \n\t"	// y2_r, y2_i
	"xvmaddadp	38, 47, 59		    	    \n\t"	// x3_r * y3_r , x3_i * y3_i
	"lxvd2x		59, %7, %3			    \n\t"	// y3_r, y3_i

	"xvmaddadp	33, 44, 60		    	    \n\t"	// x0_r * y0_i , x0_i * y0_r
	"lxvd2x		44, 0, %2			    \n\t"	// x0_r, x0_i
	"xvmaddadp	35, 45, 61		    	    \n\t"	// x1_r * y1_i , x1_i * y1_r
	"lxvd2x		45, %5, %2			    \n\t"	// x1_r, x1_i
	"xvmaddadp	37, 46, 62		    	    \n\t"	// x2_r * y2_i , x2_i * y2_r
	"lxvd2x		46, %6, %2			    \n\t"	// x2_r, x2_i
	"xvmaddadp	39, 47, 63		    	    \n\t"	// x3_r * y3_i , x3_i * y3_r
	"lxvd2x		47, %7, %2			    \n\t"	// x3_r, x3_i

	"xxswapd	60,56				    \n\t"	// y0_i, y0_r
	"xxswapd	61,57				    \n\t"	// y1_i, y1_r

	"addi		%2, %2, 64			    \n\t"
	"addi		%3, %3, 64			    \n\t"

	"xxswapd	62,58				    \n\t"	// y2_i, y2_r
	"xxswapd	63,59				    \n\t"	// y3_i, y3_r

	"addic.		%0 , %0	, -8  	 	             \n\t"
	"bgt		1b		             	     \n\t"

	"2:						     \n\t"

	"xvmaddadp	32, 40, 48		    	    \n\t"	// x0_r * y0_r , x0_i * y0_i
	"xvmaddadp	34, 41, 49		    	    \n\t"	// x1_r * y1_r , x1_i * y1_i
	"xvmaddadp	36, 42, 50		    	    \n\t"	// x2_r * y2_r , x2_i * y2_i
	"xvmaddadp	38, 43, 51		    	    \n\t"	// x3_r * y3_r , x3_i * y3_i

	"xvmaddadp	33, 40, 52		    	    \n\t"	// x0_r * y0_i , x0_i * y0_r
	"xvmaddadp	35, 41, 53		    	    \n\t"	// x1_r * y1_i , x1_i * y1_r
	"xvmaddadp	37, 42, 54		    	    \n\t"	// x2_r * y2_i , x2_i * y2_r
	"xvmaddadp	39, 43, 55		    	    \n\t"	// x3_r * y3_i , x3_i * y3_r

	"xvmaddadp	32, 44, 56		    	    \n\t"	// x0_r * y0_r , x0_i * y0_i
	"xvmaddadp	34, 45, 57		    	    \n\t"	// x1_r * y1_r , x1_i * y1_i
	"xvmaddadp	36, 46, 58		    	    \n\t"	// x2_r * y2_r , x2_i * y2_i
	"xvmaddadp	38, 47, 59		    	    \n\t"	// x3_r * y3_r , x3_i * y3_i

	"xvmaddadp	33, 44, 60		    	    \n\t"	// x0_r * y0_i , x0_i * y0_r
	"xvmaddadp	35, 45, 61		    	    \n\t"	// x1_r * y1_i , x1_i * y1_r
	"xvmaddadp	37, 46, 62		    	    \n\t"	// x2_r * y2_i , x2_i * y2_r
	"xvmaddadp	39, 47, 63		    	    \n\t"	// x3_r * y3_i , x3_i * y3_r


	"xvadddp	32, 32, 34		     \n\t"
	"xvadddp	36, 36, 38		     \n\t"

	"xvadddp	33, 33, 35		     \n\t"
	"xvadddp	37, 37, 39		     \n\t"

	"xvadddp	32, 32, 36		     \n\t"
	"xvadddp	33, 33, 37		     \n\t"

	"stxvd2x	32, 0, %4		     \n\t"
	"stxvd2x	33, %5, %4		     \n\t"

	:
        : 
          "r" (i),	// 0	
	  "r" (n),  	// 1
          "r" (x1),     // 2
          "r" (y1),     // 3
          "r" (dot),    // 4
	  "r" (o16),	// 5
	  "r" (o32),	// 6
	  "r" (o48),    // 7
	  "r" (pre)	// 8
	: "cr0", "%0", "%2" , "%3", "memory"
	);

} 


