/***************************************************************************
 * Copyright (c) 2025, The OpenBLAS Project
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in
 * the documentation and/or other materials provided with the
 * distribution.
 * 3. Neither the name of the OpenBLAS project nor the names of
 * its contributors may be used to endorse or promote products
 * derived from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * *****************************************************************************/

#include "common.h"
#include "bf16_macros.h"

int CNAME(BLASLONG bm,BLASLONG bn,BLASLONG bk,FLOAT alpha,IFLOAT* ba,IFLOAT* bb,FLOAT* C,BLASLONG ldc
#ifdef TRMMKERNEL
		,BLASLONG offset
#endif
		)
{
   BLASLONG i,j,k;
   FLOAT *C0,*C1;
   IFLOAT *ptrba,*ptrbb;
#ifdef BGEMM
   float res0,res1,res2,res3;
#else
   FLOAT res0,res1,res2,res3;
#endif
   IFLOAT load0,load1,load2,load3,load4,load5,load6,load7;
   for (j=0; j<bn/2; j+=1)
     {
        C0 = C;
        C1 = C0+ldc;
        ptrba = ba;
        for (i=0; i<bm/2; i+=1)
          {
             ptrbb = bb;
             res0 = 0;
             res1 = 0;
             res2 = 0;
             res3 = 0;
             for (k=0; k<bk/4; k+=1)
               {
                  load0 = ptrba[2*0+0];
                  load1 = ptrbb[2*0+0];
                  res0 = res0+BF16TOF32(load0)*BF16TOF32(load1);
                  load2 = ptrba[2*0+1];
                  res1 = res1+BF16TOF32(load2)*BF16TOF32(load1);
                  load3 = ptrbb[2*0+1];
                  res2 = res2+BF16TOF32(load0)*BF16TOF32(load3);
                  res3 = res3+BF16TOF32(load2)*BF16TOF32(load3);
                  load4 = ptrba[2*1+0];
                  load5 = ptrbb[2*1+0];
                  res0 = res0+BF16TOF32(load4)*BF16TOF32(load5);
                  load6 = ptrba[2*1+1];
                  res1 = res1+BF16TOF32(load6)*BF16TOF32(load5);
                  load7 = ptrbb[2*1+1];
                  res2 = res2+BF16TOF32(load4)*BF16TOF32(load7);
                  res3 = res3+BF16TOF32(load6)*BF16TOF32(load7);
                  load0 = ptrba[2*2+0];
                  load1 = ptrbb[2*2+0];
                  res0 = res0+BF16TOF32(load0)*BF16TOF32(load1);
                  load2 = ptrba[2*2+1];
                  res1 = res1+BF16TOF32(load2)*BF16TOF32(load1);
                  load3 = ptrbb[2*2+1];
                  res2 = res2+BF16TOF32(load0)*BF16TOF32(load3);
                  res3 = res3+BF16TOF32(load2)*BF16TOF32(load3);
                  load4 = ptrba[2*3+0];
                  load5 = ptrbb[2*3+0];
                  res0 = res0+BF16TOF32(load4)*BF16TOF32(load5);
                  load6 = ptrba[2*3+1];
                  res1 = res1+BF16TOF32(load6)*BF16TOF32(load5);
                  load7 = ptrbb[2*3+1];
                  res2 = res2+BF16TOF32(load4)*BF16TOF32(load7);
                  res3 = res3+BF16TOF32(load6)*BF16TOF32(load7);
                  ptrba = ptrba+8;
                  ptrbb = ptrbb+8;
               }
             for (k=0; k<(bk&3); k+=1)
               {
                  load0 = ptrba[2*0+0];
                  load1 = ptrbb[2*0+0];
                  res0 = res0+BF16TOF32(load0)*BF16TOF32(load1);
                  load2 = ptrba[2*0+1];
                  res1 = res1+BF16TOF32(load2)*BF16TOF32(load1);
                  load3 = ptrbb[2*0+1];
                  res2 = res2+BF16TOF32(load0)*BF16TOF32(load3);
                  res3 = res3+BF16TOF32(load2)*BF16TOF32(load3);
                  ptrba = ptrba+2;
                  ptrbb = ptrbb+2;
               }
             res0 = res0*ALPHA;
             C0[0] = F32TOBF16(BF16TOF32(C0[0])+res0);
             res1 = res1*ALPHA;
             C0[1] = F32TOBF16(BF16TOF32(C0[1])+res1);
             res2 = res2*ALPHA;
             C1[0] = F32TOBF16(BF16TOF32(C1[0])+res2);
             res3 = res3*ALPHA;
             C1[1] = F32TOBF16(BF16TOF32(C1[1])+res3);
             C0 = C0+2;
             C1 = C1+2;
          }
        for (i=0; i<(bm&1); i+=1)
          {
             ptrbb = bb;
             res0 = 0;
             res1 = 0;
             for (k=0; k<bk; k+=1)
               {
                  load0 = ptrba[0+0];
                  load1 = ptrbb[2*0+0];
                  res0 = res0+BF16TOF32(load0)*BF16TOF32(load1);
                  load2 = ptrbb[2*0+1];
                  res1 = res1+BF16TOF32(load0)*BF16TOF32(load2);
                  ptrba = ptrba+1;
                  ptrbb = ptrbb+2;
               }
             res0 = res0*ALPHA;
             C0[0] = F32TOBF16(BF16TOF32(C0[0])+res0);
             res1 = res1*ALPHA;
             C1[0] = F32TOBF16(BF16TOF32(C1[0])+res1);
             C0 = C0+1;
             C1 = C1+1;
          }
        k = (bk<<1);
        bb = bb+k;
        i = (ldc<<1);
        C = C+i;
     }
   for (j=0; j<(bn&1); j+=1)
     {
        C0 = C;
        ptrba = ba;
        for (i=0; i<bm/2; i+=1)
          {
             ptrbb = bb;
             res0 = 0;
             res1 = 0;
             for (k=0; k<bk; k+=1)
               {
                  load0 = ptrba[2*0+0];
                  load1 = ptrbb[0+0];
                  res0 = res0+BF16TOF32(load0)*BF16TOF32(load1);
                  load2 = ptrba[2*0+1];
                  res1 = res1+BF16TOF32(load2)*BF16TOF32(load1);
                  ptrba = ptrba+2;
                  ptrbb = ptrbb+1;
               }
             res0 = res0*ALPHA;
             C0[0] = F32TOBF16(BF16TOF32(C0[0])+res0);
             res1 = res1*ALPHA;
             C0[1] = F32TOBF16(BF16TOF32(C0[1])+res1);
             C0 = C0+2;
          }
        for (i=0; i<(bm&1); i+=1)
          {
             ptrbb = bb;
             res0 = 0;
             for (k=0; k<bk; k+=1)
               {
                  load0 = ptrba[0+0];
                  load1 = ptrbb[0+0];
                  res0 = res0+BF16TOF32(load0)*BF16TOF32(load1);
                  ptrba = ptrba+1;
                  ptrbb = ptrbb+1;
               }
             res0 = res0*ALPHA;
             C0[0] = F32TOBF16(BF16TOF32(C0[0])+res0);
             C0 = C0+1;
          }
        k = (bk<<0);
        bb = bb+k;
        C = C+ldc;
     }
   return 0;
}
