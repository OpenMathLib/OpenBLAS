/***************************************************************************
Copyright (c) 2013-2018,The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms,with or without
modification,are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice,this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice,this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES,INCLUDING,BUT NOT LIMITED TO,THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT,INDIRECT,INCIDENTAL,SPECIAL,EXEMPLARY,OR CONSEQUENTIAL
DAMAGES (INCLUDING,BUT NOT LIMITED TO,PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE,DATA,OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY,WHETHER IN CONTRACT,STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE,EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#include "common.h"

static double dsdot_kernel_32(BLASLONG n, FLOAT *x, FLOAT *y)
{
    double dot;

    __asm__ volatile (   
        "vzero %%v0                      \n\t"
        "srlg  %%r0,%1,5                 \n\t"
        "xgr   %%r1,%%r1                 \n\t"
        "0:                              \n\t"
        "pfd 1,1024(%%r1,%2)             \n\t"
        "pfd 2,1024(%%r1,%3)             \n\t"

        "vl  %%v16,0(%%r1,%2)            \n\t"
        "vl  %%v17,16(%%r1,%2)           \n\t"
        "vl  %%v18,32(%%r1,%2)           \n\t"
        "vl  %%v19,48(%%r1,%2)           \n\t"
        "vl  %%v20,64(%%r1,%2)           \n\t"
        "vl  %%v21,80(%%r1,%2)           \n\t"
        "vl  %%v22,96(%%r1,%2)           \n\t"
        "vl  %%v23,112(%%r1,%2)          \n\t"

        "vl  %%v24,0(%%r1,%3)            \n\t"
        "vfmsb   %%v16,%%v16,%%v24       \n\t"
        "vl  %%v25,16(%%r1,%3)           \n\t"
        "vfmsb   %%v17,%%v17,%%v25       \n\t"
        "vl  %%v26,32(%%r1,%3)           \n\t"
        "vfmsb   %%v18,%%v18,%%v26       \n\t"
        "vl  %%v27,48(%%r1,%3)           \n\t"
        "vfmsb   %%v19,%%v19,%%v27       \n\t"
        "vl  %%v28,64(%%r1,%3)           \n\t"
        "vfmsb   %%v20,%%v20,%%v28       \n\t"
        "vl  %%v29,80(%%r1,%3)           \n\t"
        "vfmsb   %%v21,%%v21,%%v29       \n\t"
        "vl  %%v30,96(%%r1,%3)           \n\t"
        "vfmsb   %%v22,%%v22,%%v30       \n\t"
        "vl  %%v31,112(%%r1,%3)          \n\t"
        "vfmsb   %%v23,%%v23,%%v31       \n\t"

        "vflls   %%v24,%%v16             \n\t"
        "vflls   %%v25,%%v17             \n\t"
        "vflls   %%v26,%%v18             \n\t"
        "vflls   %%v27,%%v19             \n\t"
        "vflls   %%v28,%%v20             \n\t"
        "vflls   %%v29,%%v21             \n\t"
        "vflls   %%v30,%%v22             \n\t"
        "vflls   %%v31,%%v23             \n\t"

        "veslg   %%v16,%%v16,32          \n\t"
        "veslg   %%v17,%%v17,32          \n\t"
        "veslg   %%v18,%%v18,32          \n\t"
        "veslg   %%v19,%%v19,32          \n\t"
        "veslg   %%v20,%%v20,32          \n\t"
        "veslg   %%v21,%%v21,32          \n\t"
        "veslg   %%v22,%%v22,32          \n\t"
        "veslg   %%v23,%%v23,32          \n\t"

        "vflls   %%v16,%%v16             \n\t"
        "vflls   %%v17,%%v17             \n\t"
        "vflls   %%v18,%%v18             \n\t"
        "vflls   %%v19,%%v19             \n\t"
        "vflls   %%v20,%%v20             \n\t"
        "vflls   %%v21,%%v21             \n\t"
        "vflls   %%v22,%%v22             \n\t"
        "vflls   %%v23,%%v23             \n\t"

        "vfadb   %%v16,%%v16,%%v24       \n\t"
        "vfadb   %%v17,%%v17,%%v25       \n\t"
        "vfadb   %%v18,%%v18,%%v26       \n\t"
        "vfadb   %%v19,%%v19,%%v27       \n\t"
        "vfadb   %%v20,%%v20,%%v28       \n\t"
        "vfadb   %%v21,%%v21,%%v29       \n\t"
        "vfadb   %%v22,%%v22,%%v30       \n\t"
        "vfadb   %%v23,%%v23,%%v31       \n\t"
        "vfadb   %%v16,%%v16,%%v20       \n\t"
        "vfadb   %%v17,%%v17,%%v21       \n\t"
        "vfadb   %%v18,%%v18,%%v22       \n\t"
        "vfadb   %%v19,%%v19,%%v23       \n\t"
        "vfadb   %%v16,%%v16,%%v18       \n\t"
        "vfadb   %%v17,%%v17,%%v19       \n\t"
        "vfadb   %%v16,%%v16,%%v17       \n\t"
        "vfadb   %%v0,%%v16,%%v0         \n\t"

        "agfi   %%r1,128                 \n\t"
        "brctg  %%r0,0b                  \n\t"
        "vrepg  %%v1,%%v0,1              \n\t"
        "adbr   %%f0,%%f1                \n\t"
        "ldr    %0,%%f0                      "
        :"=f"(dot)
        :"r"(n),"ZR"((const FLOAT (*)[n])x),"ZR"((const FLOAT (*)[n])y)
        :"memory","cc","r0","r1","v0","v1","v2","v3","v16","v17","v18","v19","v20","v21","v22","v23","v24","v25","v26","v27","v28","v29","v30","v31"
    );

    return dot;
}

double CNAME(BLASLONG n,FLOAT *x,BLASLONG inc_x,FLOAT *y,BLASLONG inc_y)
{
	BLASLONG i=0;
	BLASLONG ix=0,iy=0;

	double  dot = 0.0 ;

	if ( n <= 0 )  return(dot);

	if ( (inc_x == 1) && (inc_y == 1) )
	{

		BLASLONG n1 = n & -32;

		if ( n1 )
			dot = dsdot_kernel_32(n1,x,y);

		i = n1;
		while(i < n)
		{

			dot += y[i] * x[i] ;
			i++ ;

		}
		return(dot);


	}

	BLASLONG n1 = n & -2;

	while(i < n1)
	{

		dot += y[iy] * x[ix] + y[iy+inc_y] * x[ix+inc_x];
		ix  += inc_x*2 ;
		iy  += inc_y*2 ;
		i+=2 ;

	}

	while(i < n)
	{

		dot += y[iy] * x[ix] ;
		ix  += inc_x ;
		iy  += inc_y ;
		i++ ;

	}
	return(dot);

}


