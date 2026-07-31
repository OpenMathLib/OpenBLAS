/***************************************************************************
 * Copyright (c) 2026, The OpenBLAS Project
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

/*
 * Portable C GEMM micro-kernel with a 4x4 register tile (16 accumulators).
 *
 * This is a wider companion to gemmkernel_2x2.c intended for in-order scalar
 * cores whose FP FMA has a multi-cycle latency but 1/cycle throughput (e.g.
 * SiFive U74: fmadd.d latency 7, repeat rate 1).  A 2x2 tile exposes only 4
 * independent accumulator chains, which is fewer than the FMA latency and
 * leaves the FP pipe stalled on the accumulator dependency.  A 4x4 tile keeps
 * 16 independent chains -- comfortably above the latency -- and lowers the
 * load:FMA ratio from 1:1 to 1:2, so the single load/store pipe stops being
 * the bottleneck.  RV64G has 32 FP registers, so 16 accumulators + 4 A + 4 B
 * fit without spilling.
 *
 * Packed-data contract (identical to the 2x2 kernel, verified against the
 * generic tcopy_4 / ncopy_4 copy routines): the A operand is packed by
 * tcopy_<UNROLL_M> into MR-row micro-panels [A(r0,k)..A(r3,k)] per k, and the
 * B operand by ncopy_<UNROLL_N> into NR-col micro-panels [B(k,c0)..B(k,c3)]
 * per k.  Both dimensions are decomposed as 4 / 2 / 1 sub-blocks at the edges.
 */

#include "common.h"

#include "conversion_macros.h"

#ifdef BGEMM
#define C_TO_F32 TO_F32
#else
#define C_TO_F32
#endif

int CNAME(BLASLONG bm,BLASLONG bn,BLASLONG bk,FLOAT alpha,IFLOAT* ba,IFLOAT* bb,FLOAT* C,BLASLONG ldc
#ifdef TRMMKERNEL
		,BLASLONG offset
#endif
		)
{
   BLASLONG i,j,k;
   FLOAT *C0,*C1,*C2,*C3;
   IFLOAT *ptrba,*ptrbb;
   FLOAT r0c0,r1c0,r2c0,r3c0;
   FLOAT r0c1,r1c1,r2c1,r3c1;
   FLOAT r0c2,r1c2,r2c2,r3c2;
   FLOAT r0c3,r1c3,r2c3,r3c3;
   IFLOAT a0,a1,a2,a3,b0,b1,b2,b3;

   /* ==================== N panels of 4 ==================== */
   for (j=0; j<bn/4; j+=1)
     {
        C0 = C;
        C1 = C0+ldc;
        C2 = C1+ldc;
        C3 = C2+ldc;
        ptrba = ba;

        /* ---- 4x4 : 4 rows x 4 cols, 16 accumulators ---- */
        for (i=0; i<bm/4; i+=1)
          {
             ptrbb = bb;
             r0c0=r1c0=r2c0=r3c0=0;
             r0c1=r1c1=r2c1=r3c1=0;
             r0c2=r1c2=r2c2=r3c2=0;
             r0c3=r1c3=r2c3=r3c3=0;
             for (k=0; k<bk; k+=1)
               {
                  b0=ptrbb[0]; b1=ptrbb[1]; b2=ptrbb[2]; b3=ptrbb[3];
                  a0=ptrba[0]; a1=ptrba[1]; a2=ptrba[2]; a3=ptrba[3];
                  r0c0+=TO_F32(a0)*TO_F32(b0); r1c0+=TO_F32(a1)*TO_F32(b0); r2c0+=TO_F32(a2)*TO_F32(b0); r3c0+=TO_F32(a3)*TO_F32(b0);
                  r0c1+=TO_F32(a0)*TO_F32(b1); r1c1+=TO_F32(a1)*TO_F32(b1); r2c1+=TO_F32(a2)*TO_F32(b1); r3c1+=TO_F32(a3)*TO_F32(b1);
                  r0c2+=TO_F32(a0)*TO_F32(b2); r1c2+=TO_F32(a1)*TO_F32(b2); r2c2+=TO_F32(a2)*TO_F32(b2); r3c2+=TO_F32(a3)*TO_F32(b2);
                  r0c3+=TO_F32(a0)*TO_F32(b3); r1c3+=TO_F32(a1)*TO_F32(b3); r2c3+=TO_F32(a2)*TO_F32(b3); r3c3+=TO_F32(a3)*TO_F32(b3);
                  ptrba+=4; ptrbb+=4;
               }
             C0[0]=TO_OUTPUT(C_TO_F32(C0[0])+r0c0*ALPHA); C0[1]=TO_OUTPUT(C_TO_F32(C0[1])+r1c0*ALPHA); C0[2]=TO_OUTPUT(C_TO_F32(C0[2])+r2c0*ALPHA); C0[3]=TO_OUTPUT(C_TO_F32(C0[3])+r3c0*ALPHA);
             C1[0]=TO_OUTPUT(C_TO_F32(C1[0])+r0c1*ALPHA); C1[1]=TO_OUTPUT(C_TO_F32(C1[1])+r1c1*ALPHA); C1[2]=TO_OUTPUT(C_TO_F32(C1[2])+r2c1*ALPHA); C1[3]=TO_OUTPUT(C_TO_F32(C1[3])+r3c1*ALPHA);
             C2[0]=TO_OUTPUT(C_TO_F32(C2[0])+r0c2*ALPHA); C2[1]=TO_OUTPUT(C_TO_F32(C2[1])+r1c2*ALPHA); C2[2]=TO_OUTPUT(C_TO_F32(C2[2])+r2c2*ALPHA); C2[3]=TO_OUTPUT(C_TO_F32(C2[3])+r3c2*ALPHA);
             C3[0]=TO_OUTPUT(C_TO_F32(C3[0])+r0c3*ALPHA); C3[1]=TO_OUTPUT(C_TO_F32(C3[1])+r1c3*ALPHA); C3[2]=TO_OUTPUT(C_TO_F32(C3[2])+r2c3*ALPHA); C3[3]=TO_OUTPUT(C_TO_F32(C3[3])+r3c3*ALPHA);
             C0+=4; C1+=4; C2+=4; C3+=4;
          }
        /* ---- 2x4 : 2 rows x 4 cols ---- */
        if (bm & 2)
          {
             ptrbb = bb;
             r0c0=r1c0=0; r0c1=r1c1=0; r0c2=r1c2=0; r0c3=r1c3=0;
             for (k=0; k<bk; k+=1)
               {
                  b0=ptrbb[0]; b1=ptrbb[1]; b2=ptrbb[2]; b3=ptrbb[3];
                  a0=ptrba[0]; a1=ptrba[1];
                  r0c0+=TO_F32(a0)*TO_F32(b0); r1c0+=TO_F32(a1)*TO_F32(b0);
                  r0c1+=TO_F32(a0)*TO_F32(b1); r1c1+=TO_F32(a1)*TO_F32(b1);
                  r0c2+=TO_F32(a0)*TO_F32(b2); r1c2+=TO_F32(a1)*TO_F32(b2);
                  r0c3+=TO_F32(a0)*TO_F32(b3); r1c3+=TO_F32(a1)*TO_F32(b3);
                  ptrba+=2; ptrbb+=4;
               }
             C0[0]=TO_OUTPUT(C_TO_F32(C0[0])+r0c0*ALPHA); C0[1]=TO_OUTPUT(C_TO_F32(C0[1])+r1c0*ALPHA);
             C1[0]=TO_OUTPUT(C_TO_F32(C1[0])+r0c1*ALPHA); C1[1]=TO_OUTPUT(C_TO_F32(C1[1])+r1c1*ALPHA);
             C2[0]=TO_OUTPUT(C_TO_F32(C2[0])+r0c2*ALPHA); C2[1]=TO_OUTPUT(C_TO_F32(C2[1])+r1c2*ALPHA);
             C3[0]=TO_OUTPUT(C_TO_F32(C3[0])+r0c3*ALPHA); C3[1]=TO_OUTPUT(C_TO_F32(C3[1])+r1c3*ALPHA);
             C0+=2; C1+=2; C2+=2; C3+=2;
          }
        /* ---- 1x4 : 1 row x 4 cols ---- */
        if (bm & 1)
          {
             ptrbb = bb;
             r0c0=0; r0c1=0; r0c2=0; r0c3=0;
             for (k=0; k<bk; k+=1)
               {
                  b0=ptrbb[0]; b1=ptrbb[1]; b2=ptrbb[2]; b3=ptrbb[3];
                  a0=ptrba[0];
                  r0c0+=TO_F32(a0)*TO_F32(b0);
                  r0c1+=TO_F32(a0)*TO_F32(b1);
                  r0c2+=TO_F32(a0)*TO_F32(b2);
                  r0c3+=TO_F32(a0)*TO_F32(b3);
                  ptrba+=1; ptrbb+=4;
               }
             C0[0]=TO_OUTPUT(C_TO_F32(C0[0])+r0c0*ALPHA);
             C1[0]=TO_OUTPUT(C_TO_F32(C1[0])+r0c1*ALPHA);
             C2[0]=TO_OUTPUT(C_TO_F32(C2[0])+r0c2*ALPHA);
             C3[0]=TO_OUTPUT(C_TO_F32(C3[0])+r0c3*ALPHA);
             C0+=1; C1+=1; C2+=1; C3+=1;
          }
        bb = bb + bk*4;
        C  = C  + ldc*4;
     }

   /* ==================== N panel of 2 ==================== */
   if (bn & 2)
     {
        C0 = C;
        C1 = C0+ldc;
        ptrba = ba;

        for (i=0; i<bm/4; i+=1)
          {
             ptrbb = bb;
             r0c0=r1c0=r2c0=r3c0=0;
             r0c1=r1c1=r2c1=r3c1=0;
             for (k=0; k<bk; k+=1)
               {
                  b0=ptrbb[0]; b1=ptrbb[1];
                  a0=ptrba[0]; a1=ptrba[1]; a2=ptrba[2]; a3=ptrba[3];
                  r0c0+=TO_F32(a0)*TO_F32(b0); r1c0+=TO_F32(a1)*TO_F32(b0); r2c0+=TO_F32(a2)*TO_F32(b0); r3c0+=TO_F32(a3)*TO_F32(b0);
                  r0c1+=TO_F32(a0)*TO_F32(b1); r1c1+=TO_F32(a1)*TO_F32(b1); r2c1+=TO_F32(a2)*TO_F32(b1); r3c1+=TO_F32(a3)*TO_F32(b1);
                  ptrba+=4; ptrbb+=2;
               }
             C0[0]=TO_OUTPUT(C_TO_F32(C0[0])+r0c0*ALPHA); C0[1]=TO_OUTPUT(C_TO_F32(C0[1])+r1c0*ALPHA); C0[2]=TO_OUTPUT(C_TO_F32(C0[2])+r2c0*ALPHA); C0[3]=TO_OUTPUT(C_TO_F32(C0[3])+r3c0*ALPHA);
             C1[0]=TO_OUTPUT(C_TO_F32(C1[0])+r0c1*ALPHA); C1[1]=TO_OUTPUT(C_TO_F32(C1[1])+r1c1*ALPHA); C1[2]=TO_OUTPUT(C_TO_F32(C1[2])+r2c1*ALPHA); C1[3]=TO_OUTPUT(C_TO_F32(C1[3])+r3c1*ALPHA);
             C0+=4; C1+=4;
          }
        if (bm & 2)
          {
             ptrbb = bb;
             r0c0=r1c0=0; r0c1=r1c1=0;
             for (k=0; k<bk; k+=1)
               {
                  b0=ptrbb[0]; b1=ptrbb[1];
                  a0=ptrba[0]; a1=ptrba[1];
                  r0c0+=TO_F32(a0)*TO_F32(b0); r1c0+=TO_F32(a1)*TO_F32(b0);
                  r0c1+=TO_F32(a0)*TO_F32(b1); r1c1+=TO_F32(a1)*TO_F32(b1);
                  ptrba+=2; ptrbb+=2;
               }
             C0[0]=TO_OUTPUT(C_TO_F32(C0[0])+r0c0*ALPHA); C0[1]=TO_OUTPUT(C_TO_F32(C0[1])+r1c0*ALPHA);
             C1[0]=TO_OUTPUT(C_TO_F32(C1[0])+r0c1*ALPHA); C1[1]=TO_OUTPUT(C_TO_F32(C1[1])+r1c1*ALPHA);
             C0+=2; C1+=2;
          }
        if (bm & 1)
          {
             ptrbb = bb;
             r0c0=0; r0c1=0;
             for (k=0; k<bk; k+=1)
               {
                  b0=ptrbb[0]; b1=ptrbb[1];
                  a0=ptrba[0];
                  r0c0+=TO_F32(a0)*TO_F32(b0);
                  r0c1+=TO_F32(a0)*TO_F32(b1);
                  ptrba+=1; ptrbb+=2;
               }
             C0[0]=TO_OUTPUT(C_TO_F32(C0[0])+r0c0*ALPHA);
             C1[0]=TO_OUTPUT(C_TO_F32(C1[0])+r0c1*ALPHA);
             C0+=1; C1+=1;
          }
        bb = bb + bk*2;
        C  = C  + ldc*2;
     }

   /* ==================== N panel of 1 ==================== */
   if (bn & 1)
     {
        C0 = C;
        ptrba = ba;

        for (i=0; i<bm/4; i+=1)
          {
             ptrbb = bb;
             r0c0=r1c0=r2c0=r3c0=0;
             for (k=0; k<bk; k+=1)
               {
                  b0=ptrbb[0];
                  a0=ptrba[0]; a1=ptrba[1]; a2=ptrba[2]; a3=ptrba[3];
                  r0c0+=TO_F32(a0)*TO_F32(b0); r1c0+=TO_F32(a1)*TO_F32(b0); r2c0+=TO_F32(a2)*TO_F32(b0); r3c0+=TO_F32(a3)*TO_F32(b0);
                  ptrba+=4; ptrbb+=1;
               }
             C0[0]=TO_OUTPUT(C_TO_F32(C0[0])+r0c0*ALPHA); C0[1]=TO_OUTPUT(C_TO_F32(C0[1])+r1c0*ALPHA); C0[2]=TO_OUTPUT(C_TO_F32(C0[2])+r2c0*ALPHA); C0[3]=TO_OUTPUT(C_TO_F32(C0[3])+r3c0*ALPHA);
             C0+=4;
          }
        if (bm & 2)
          {
             ptrbb = bb;
             r0c0=r1c0=0;
             for (k=0; k<bk; k+=1)
               {
                  b0=ptrbb[0];
                  a0=ptrba[0]; a1=ptrba[1];
                  r0c0+=TO_F32(a0)*TO_F32(b0); r1c0+=TO_F32(a1)*TO_F32(b0);
                  ptrba+=2; ptrbb+=1;
               }
             C0[0]=TO_OUTPUT(C_TO_F32(C0[0])+r0c0*ALPHA); C0[1]=TO_OUTPUT(C_TO_F32(C0[1])+r1c0*ALPHA);
             C0+=2;
          }
        if (bm & 1)
          {
             ptrbb = bb;
             r0c0=0;
             for (k=0; k<bk; k+=1)
               {
                  r0c0+=TO_F32(ptrba[0])*TO_F32(ptrbb[0]);
                  ptrba+=1; ptrbb+=1;
               }
             C0[0]=TO_OUTPUT(C_TO_F32(C0[0])+r0c0*ALPHA);
             C0+=1;
          }
        bb = bb + bk;
        C  = C  + ldc;
     }

   return 0;
}
