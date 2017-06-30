/***************************************************************************
Copyright (c) 2017, The OpenBLAS Project
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

#include "common.h"
#include <stdbool.h>

int CNAME(BLASLONG bm,BLASLONG bn,BLASLONG bk,FLOAT alpha,FLOAT* ba,FLOAT* bb,FLOAT* C,BLASLONG ldc)
{

   BLASLONG i,j,k;
   FLOAT *C0,*C1,*ptrba,*ptrbb;

   FLOAT res0_0;
   FLOAT res0_1;
   FLOAT res0_2;
   FLOAT res0_3;

   FLOAT res1_0;
   FLOAT res1_1;
   FLOAT res1_2;
   FLOAT res1_3;

   FLOAT a0;
   FLOAT a1;

   FLOAT b0;
   FLOAT b1;

   for (j=0; j<(bn/2); j+=2)
   {
        C0 = C;
        C1 = C0+ldc;

        ptrba = ba;

        for (i=0; i<bm/4; i+=1)
	{
		ptrbb = bb;

		res0_0 = 0;
		res0_1 = 0;
		res0_2 = 0;
		res0_3 = 0;

		res1_0 = 0;
		res1_1 = 0;
		res1_2 = 0;
		res1_3 = 0;

		for (k=0; k<bk; k++)
                {
			b0 = ptrbb[0];
			b1 = ptrbb[1];

			a0 = ptrba[0];
			res0_0 += a0*b0;
			res1_0 += a0*b1;

			a1 = ptrba[1];
			res0_1 += a1*b0;
			res1_1 += a1*b1;

			a0 = ptrba[2];
			res0_2 += a0*b0;
			res1_2 += a0*b1;

			a1 = ptrba[3];
			res0_3 += a1*b0;
			res1_3 += a1*b1;

			ptrba = ptrba+4;
			ptrbb = ptrbb+2;
                }

		res0_0 *= alpha;
		res0_1 *= alpha;
		res0_2 *= alpha;
		res0_3 *= alpha;

		res1_0 *= alpha;
		res1_1 *= alpha;
		res1_2 *= alpha;
		res1_3 *= alpha;

		C0[0] += res0_0;
		C0[1] += res0_1;
		C0[2] += res0_2;
		C0[3] += res0_3;

		C1[0] += res1_0;
		C1[1] += res1_1;
		C1[2] += res1_2;
		C1[3] += res1_3;

		C0 = C0+4;
		C1 = C1+4;

	}

	if ( bm & 2 )
	{
		ptrbb = bb;

		res0_0 = 0;
		res0_1 = 0;

		res1_0 = 0;
		res1_1 = 0;

		for (k=0; k<bk; k++)
                {
			b0 = ptrbb[0];
			b1 = ptrbb[1];

			a0 = ptrba[0];
			res0_0 += a0*b0;
			res1_0 += a0*b1;

			a1 = ptrba[1];
			res0_1 += a1*b0;
			res1_1 += a1*b1;

			ptrba = ptrba+2;
			ptrbb = ptrbb+2;
                }

		res0_0 *= alpha;
		res0_1 *= alpha;

		res1_0 *= alpha;
		res1_1 *= alpha;

		C0[0] = res0_0;
		C0[1] = res0_1;

		C1[0] = res1_0;
		C1[1] = res1_1;

		C0 = C0+2;
		C1 = C1+2;

	}

	if ( bm & 1 )
	{
		ptrbb = bb;

		res0_0 = 0;

		res1_0 = 0;

		for (k=0; k<bk; k++)
                {
			b0 = ptrbb[0];
			b1 = ptrbb[1];

			a0 = ptrba[0];
			res0_0 += a0*b0;
			res1_0 += a0*b1;

			ptrba = ptrba+1;
			ptrbb = ptrbb+2;
                }

		res0_0 *= alpha;

		res1_0 *= alpha;

		C0[0] = res0_0;

		C1[0] = res1_0;

		C0 = C0+1;
		C1 = C1+1;

	}

        k = (bk<<1);
        bb = bb+k;
        i = (ldc<<1);
        C = C+i;
    }

   for (j=0; j<(bn&1); j+=1) // do the Mx1 loops
   {
        C0 = C;

        ptrba = ba;

        for (i=0; i<bm/4; i+=1) // do blocks of 4x1 loops
	{
		ptrbb = bb;

		res0_0 = 0;
		res0_1 = 0;
		res0_2 = 0;
		res0_3 = 0;

		for (k=0; k<bk; k++)
                {
			b0 = ptrbb[0];

			a0 = ptrba[0];
			res0_0 += a0*b0;

			a1 = ptrba[1];
			res0_1 += a1*b0;

			a0 = ptrba[2];
			res0_2 += a0*b0;

			a1 = ptrba[3];
			res0_3 += a1*b0;

			ptrba = ptrba+4;
			ptrbb = ptrbb+1;
                }

		res0_0 *= alpha;
		res0_1 *= alpha;
		res0_2 *= alpha;
		res0_3 *= alpha;

		C0[0] = res0_0;
		C0[1] = res0_1;
		C0[2] = res0_2;
		C0[3] = res0_3;

		C0 = C0+4;

	}

	if ( bm & 2 ) // do any 2x1 loop
	{
		ptrbb = bb;

		res0_0 = 0;
		res0_1 = 0;

		for (k=0; k<bk; k++)
                {
			b0 = ptrbb[0];

			a0 = ptrba[0];
			res0_0 += a0*b0;

			a1 = ptrba[1];
			res0_1 += a1*b0;

			ptrba = ptrba+2;
			ptrbb = ptrbb+1;
                }

		res0_0 *= alpha;
		res0_1 *= alpha;

		C0[0] = res0_0;
		C0[1] = res0_1;

		C0 = C0+2;

	}

	if ( bm & 1 ) // do any 1x1 loop
	{
		ptrbb = bb;

		res0_0 = 0;

		for (k=0; k<bk; k++)
                {
			b0 = ptrbb[0];

			a0 = ptrba[0];
			res0_0 += a0*b0;

			ptrba = ptrba+1;
			ptrbb = ptrbb+1;
                }

		res0_0 *= alpha;

		C0[0] = res0_0;

		C0 = C0+1;

	}
        k = (bk<<0);
        bb = bb+k;
        C = C+ldc;
   }
   return 0;
}
