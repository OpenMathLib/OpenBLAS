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
* 2016/03/30 Werner Saar (wernsaar@googlemail.com)
* 	 BLASTEST 		: OK
* 	 CTEST			: OK
* 	 TEST			: OK
*	 LAPACK-TEST		: OK
**************************************************************************************/

#define HAVE_KERNEL_4x4 1

static void dgemv_kernel_4x4(BLASLONG n, FLOAT **ap, FLOAT *xo, FLOAT *y, FLOAT *alpha) __attribute__ ((noinline));

static void dgemv_kernel_4x4(BLASLONG n, FLOAT **ap, FLOAT *xo, FLOAT *y, FLOAT *alpha)
{
        BLASLONG i=n;
	BLASLONG o8  = 8;
	BLASLONG o16 = 16;
	BLASLONG o24 = 24;
	BLASLONG pre = 384;

        FLOAT *a0,*a1,*a2,*a3;
	FLOAT *y1=y+1;
        FLOAT x[4]  __attribute__ ((aligned (16)));;
        a0 = ap[0]+1;
        a1 = ap[1]+1;
        a2 = ap[2]+1;
        a3 = ap[3]+1;

	x[0]=xo[0] * *alpha;
	x[1]=xo[1] * *alpha;
	x[2]=xo[2] * *alpha;
	x[3]=xo[3] * *alpha;


	__asm__  __volatile__
	(
	"lxvdsx		32, 0 , %1			    \n\t"	// x0
	"lxvdsx		33,%3 , %1			    \n\t"	// x1
	"lxvdsx		34,%4 , %1			    \n\t"	// x2
	"lxvdsx		35,%5 , %1			    \n\t"	// x3
	"addi		%2 , %2 , -8			    \n\t"
	"addi		%6 , %6 , -8			    \n\t"
	"addi		%7 , %7 , -8			    \n\t"
	"addi		%8 , %8 , -8			    \n\t"
	"addi		%9 , %9 , -8			    \n\t"
	
	"lxvd2x		48, 0, %6			    \n\t"	// a0[0], a0[1] 
	"lxvd2x		49,%4, %6			    \n\t"	// a0[2], a0[3] 

	"lxvd2x		50, 0, %7			    \n\t"	// a1[0], a1[1] 
	"lxvd2x		51,%4, %7			    \n\t"	// a1[2], a1[3] 

	"lxvd2x		52, 0, %8			    \n\t"	// a2[0], a2[1] 
	"lxvd2x		53,%4, %8			    \n\t"	// a2[2], a2[3] 

	"lxvd2x		54, 0, %9			    \n\t"	// a3[0], a3[1] 
	"lxvd2x		55,%4, %9			    \n\t"	// a3[2], a3[3] 

	"addi		%6, %6, 32			    \n\t"
	"addi		%7, %7, 32			    \n\t"
	"addi		%8, %8, 32			    \n\t"
	"addi		%9, %9, 32			    \n\t"

	"addic.		%0 , %0	, -4  	 	             \n\t"
	"ble		2f		             	     \n\t"

	".align 5				            \n\t"
	"1:				                    \n\t"

	"dcbt		%2, %10				    \n\t"

	"lxvd2x		40, 0, %2			    \n\t"	// y0, y1
	"lxvd2x		41,%4, %2			    \n\t"	// y2, y3
	
	"dcbt		%6, %10				    \n\t"
	"dcbt		%7, %10				    \n\t"
	"dcbt		%8, %10				    \n\t"
	"dcbt		%9, %10				    \n\t"

	"xvmaddadp	40, 48, 32			    \n\t"
	"xvmaddadp	41, 49, 32			    \n\t"

	"lxvd2x		48, 0, %6			    \n\t"	// a0[0], a0[1] 
	"lxvd2x		49,%4, %6			    \n\t"	// a0[2], a0[3] 

	"xvmaddadp	40, 50, 33			    \n\t"
	"addi		%6, %6, 32			    \n\t"
	"xvmaddadp	41, 51, 33			    \n\t"

	"lxvd2x		50, 0, %7			    \n\t"	// a1[0], a1[1] 
	"lxvd2x		51,%4, %7			    \n\t"	// a1[2], a1[3] 

	"xvmaddadp	40, 52, 34			    \n\t"
	"addi		%7, %7, 32			    \n\t"
	"xvmaddadp	41, 53, 34			    \n\t"

	"lxvd2x		52, 0, %8			    \n\t"	// a2[0], a2[1] 
	"lxvd2x		53,%4, %8			    \n\t"	// a2[2], a2[3] 

	"xvmaddadp	40, 54, 35			    \n\t"
	"addi		%8, %8, 32			    \n\t"
	"xvmaddadp	41, 55, 35			    \n\t"

	"stxvd2x	40, 0, %2			    \n\t"	// y0, y1
	"stxvd2x	41,%4, %2			    \n\t"	// y2, y3

	"lxvd2x		54, 0, %9			    \n\t"	// a3[0], a3[1] 
	"lxvd2x		55,%4, %9			    \n\t"	// a3[2], a3[3] 

	"addi		%9, %9, 32			    \n\t"
	"addi		%2, %2, 32			    \n\t"

	"addic.		%0 , %0	, -4  	 	             \n\t"
	"ble		2f		             	     \n\t"


	"lxvd2x		40, 0, %2			    \n\t"	// y0, y1
	"lxvd2x		41,%4, %2			    \n\t"	// y2, y3
	
	"xvmaddadp	40, 48, 32			    \n\t"
	"xvmaddadp	41, 49, 32			    \n\t"

	"lxvd2x		48, 0, %6			    \n\t"	// a0[0], a0[1] 
	"lxvd2x		49,%4, %6			    \n\t"	// a0[2], a0[3] 

	"xvmaddadp	40, 50, 33			    \n\t"
	"addi		%6, %6, 32			    \n\t"
	"xvmaddadp	41, 51, 33			    \n\t"

	"lxvd2x		50, 0, %7			    \n\t"	// a1[0], a1[1] 
	"lxvd2x		51,%4, %7			    \n\t"	// a1[2], a1[3] 

	"xvmaddadp	40, 52, 34			    \n\t"
	"addi		%7, %7, 32			    \n\t"
	"xvmaddadp	41, 53, 34			    \n\t"

	"lxvd2x		52, 0, %8			    \n\t"	// a2[0], a2[1] 
	"lxvd2x		53,%4, %8			    \n\t"	// a2[2], a2[3] 

	"xvmaddadp	40, 54, 35			    \n\t"
	"addi		%8, %8, 32			    \n\t"
	"xvmaddadp	41, 55, 35			    \n\t"

	"stxvd2x	40, 0, %2			    \n\t"	// y0, y1
	"stxvd2x	41,%4, %2			    \n\t"	// y2, y3

	"lxvd2x		54, 0, %9			    \n\t"	// a3[0], a3[1] 
	"lxvd2x		55,%4, %9			    \n\t"	// a3[2], a3[3] 

	"addi		%9, %9, 32			    \n\t"
	"addi		%2, %2, 32			    \n\t"

	"addic.		%0 , %0	, -4  	 	             \n\t"
	"ble		2f		             	     \n\t"


	"lxvd2x		40, 0, %2			    \n\t"	// y0, y1
	"lxvd2x		41,%4, %2			    \n\t"	// y2, y3
	
	"xvmaddadp	40, 48, 32			    \n\t"
	"xvmaddadp	41, 49, 32			    \n\t"

	"lxvd2x		48, 0, %6			    \n\t"	// a0[0], a0[1] 
	"lxvd2x		49,%4, %6			    \n\t"	// a0[2], a0[3] 

	"xvmaddadp	40, 50, 33			    \n\t"
	"addi		%6, %6, 32			    \n\t"
	"xvmaddadp	41, 51, 33			    \n\t"

	"lxvd2x		50, 0, %7			    \n\t"	// a1[0], a1[1] 
	"lxvd2x		51,%4, %7			    \n\t"	// a1[2], a1[3] 

	"xvmaddadp	40, 52, 34			    \n\t"
	"addi		%7, %7, 32			    \n\t"
	"xvmaddadp	41, 53, 34			    \n\t"

	"lxvd2x		52, 0, %8			    \n\t"	// a2[0], a2[1] 
	"lxvd2x		53,%4, %8			    \n\t"	// a2[2], a2[3] 

	"xvmaddadp	40, 54, 35			    \n\t"
	"addi		%8, %8, 32			    \n\t"
	"xvmaddadp	41, 55, 35			    \n\t"

	"stxvd2x	40, 0, %2			    \n\t"	// y0, y1
	"stxvd2x	41,%4, %2			    \n\t"	// y2, y3

	"lxvd2x		54, 0, %9			    \n\t"	// a3[0], a3[1] 
	"lxvd2x		55,%4, %9			    \n\t"	// a3[2], a3[3] 

	"addi		%9, %9, 32			    \n\t"
	"addi		%2, %2, 32			    \n\t"

	"addic.		%0 , %0	, -4  	 	             \n\t"
	"ble		2f		             	     \n\t"


	"lxvd2x		40, 0, %2			    \n\t"	// y0, y1
	"lxvd2x		41,%4, %2			    \n\t"	// y2, y3
	
	"xvmaddadp	40, 48, 32			    \n\t"
	"xvmaddadp	41, 49, 32			    \n\t"

	"lxvd2x		48, 0, %6			    \n\t"	// a0[0], a0[1] 
	"lxvd2x		49,%4, %6			    \n\t"	// a0[2], a0[3] 

	"xvmaddadp	40, 50, 33			    \n\t"
	"addi		%6, %6, 32			    \n\t"
	"xvmaddadp	41, 51, 33			    \n\t"

	"lxvd2x		50, 0, %7			    \n\t"	// a1[0], a1[1] 
	"lxvd2x		51,%4, %7			    \n\t"	// a1[2], a1[3] 

	"xvmaddadp	40, 52, 34			    \n\t"
	"addi		%7, %7, 32			    \n\t"
	"xvmaddadp	41, 53, 34			    \n\t"

	"lxvd2x		52, 0, %8			    \n\t"	// a2[0], a2[1] 
	"lxvd2x		53,%4, %8			    \n\t"	// a2[2], a2[3] 

	"xvmaddadp	40, 54, 35			    \n\t"
	"addi		%8, %8, 32			    \n\t"
	"xvmaddadp	41, 55, 35			    \n\t"

	"stxvd2x	40, 0, %2			    \n\t"	// y0, y1
	"stxvd2x	41,%4, %2			    \n\t"	// y2, y3

	"lxvd2x		54, 0, %9			    \n\t"	// a3[0], a3[1] 
	"lxvd2x		55,%4, %9			    \n\t"	// a3[2], a3[3] 

	"addi		%9, %9, 32			    \n\t"
	"addi		%2, %2, 32			    \n\t"

	"addic.		%0 , %0	, -4  	 	             \n\t"
	"bgt		1b		             	     \n\t"

	"2:						     \n\t"

	"lxvd2x		40, 0, %2			    \n\t"	// y0, y1
	"lxvd2x		41,%4, %2			    \n\t"	// y2, y3

	"xvmaddadp	40, 48, 32			    \n\t"
	"xvmaddadp	41, 49, 32			    \n\t"

	"xvmaddadp	40, 50, 33			    \n\t"
	"xvmaddadp	41, 51, 33			    \n\t"

	"xvmaddadp	40, 52, 34			    \n\t"
	"xvmaddadp	41, 53, 34			    \n\t"

	"xvmaddadp	40, 54, 35			    \n\t"
	"xvmaddadp	41, 55, 35			    \n\t"

	"stxvd2x	40, 0, %2			    \n\t"	// y0, y1
	"stxvd2x	41,%4, %2			    \n\t"	// y2, y3

	:
        : 
          "r" (i),	// 0	
          "r" (x),      // 1
          "r" (y1),     // 2
	  "r" (o8),	// 3
	  "r" (o16),	// 4
	  "r" (o24),	// 5
	  "r" (a0),	// 6
	  "r" (a1),	// 7
	  "r" (a2),	// 8
	  "r" (a3),	// 9
	  "r" (pre)	// 10
	: "cr0", "%0", "%2" , "%6", "%7", "%8", "%9", "memory"
	);

} 


