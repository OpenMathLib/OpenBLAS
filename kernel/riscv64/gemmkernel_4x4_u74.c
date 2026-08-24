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

/*
 * U74 hybrid variant.  The DGEMM fast path (bm, bn multiples of 4, even bk,
 * non-TRMM) is served by a hand-scheduled scalar assembly micro-kernel,
 * dgemm_u74_asm, embedded below via a top-level __asm__ (readable source in
 * kern_u74.S in this directory): a 4x4 register tile with full operand
 * double-buffering (P/Q ping-pong) and load-before-FMA issue ordering that
 * reaches the U74 FP-pipe peak (~16.5 cycles per 16 fmadd.d).  All other
 * shapes, odd bk, and the TRMM builds fall back to the portable C 4x4 code
 * below.  Measured ~+20% single-core DGEMM over the C kernel; four-core HPL vs
 * the tuned C kernel is +4% at N=10000 and +10.7% at N=27456 (5.99 GFLOPS,
 * ~50% of the 12 GF peak - the best clean figure).  Validated against the full
 * BLAS Level-3 suite (DGEMM 17,496 calls, 0 failures) and HPL (residual PASSED).
 */

#include "common.h"

#include "../generic/conversion_macros.h"

#ifdef BGEMM
#define C_TO_F32 TO_F32
#else
#define C_TO_F32
#endif


extern int dgemm_u74_asm(BLASLONG,BLASLONG,BLASLONG,FLOAT,IFLOAT*,IFLOAT*,FLOAT*,BLASLONG);

__asm__(
	"	.text\n"
	"	.p2align 4\n"
	"	.type	dgemm_u74_asm,@function\n"
	"dgemm_u74_asm:\n"
	"	addi	sp, sp, -112\n"
	"	fsd	fs0, 0(sp)\n"
	"	fsd	fs1, 8(sp)\n"
	"	fsd	fs2, 16(sp)\n"
	"	fsd	fs3, 24(sp)\n"
	"	fsd	fs4, 32(sp)\n"
	"	fsd	fs5, 40(sp)\n"
	"	fsd	fs6, 48(sp)\n"
	"	fsd	fs7, 56(sp)\n"
	"	fsd	fs8, 64(sp)\n"
	"	fsd	fs9, 72(sp)\n"
	"	fsd	fs10, 80(sp)\n"
	"	fsd	fs11, 88(sp)\n"
	"	fsd	fa0, 96(sp)\n"
	"	slli	t0, a6, 3\n"
	"	srli	t1, a0, 2\n"
	"	srli	t2, a1, 2\n"
	"	beqz	t2, .Lk4done\n"
	".Lk4j:\n"
	"	mv	a0, a5\n"
	"	add	a1, a0, t0\n"
	"	add	a6, a1, t0\n"
	"	add	a7, a6, t0\n"
	"	mv	t3, a3\n"
	"	mv	t5, t1\n"
	".Lk4i:\n"
	"	mv	t4, a4\n"
	"	fmv.d.x	ft0, zero\n"
	"	fmv.d.x	ft1, zero\n"
	"	fmv.d.x	ft2, zero\n"
	"	fmv.d.x	ft3, zero\n"
	"	fmv.d.x	ft4, zero\n"
	"	fmv.d.x	ft5, zero\n"
	"	fmv.d.x	ft6, zero\n"
	"	fmv.d.x	ft7, zero\n"
	"	fmv.d.x	ft8, zero\n"
	"	fmv.d.x	ft9, zero\n"
	"	fmv.d.x	ft10, zero\n"
	"	fmv.d.x	ft11, zero\n"
	"	fmv.d.x	fa2, zero\n"
	"	fmv.d.x	fa3, zero\n"
	"	fmv.d.x	fa4, zero\n"
	"	fmv.d.x	fa5, zero\n"
	"	fld	fa0, 0(t3)\n"
	"	fld	fa1, 8(t3)\n"
	"	fld	fa6, 16(t3)\n"
	"	fld	fa7, 24(t3)\n"
	"	fld	fs0, 0(t4)\n"
	"	fld	fs1, 8(t4)\n"
	"	fld	fs2, 16(t4)\n"
	"	fld	fs3, 24(t4)\n"
	"	addi	t3, t3, 32\n"
	"	addi	t4, t4, 32\n"
	"	srli	t6, a2, 1\n"
	"	addi	t6, t6, -1\n"
	"	beqz	t6, .Lk4epi\n"
	"	.p2align 4\n"
	".Lk4body:\n"
	"	fld	fs4, 0(t3)\n"
	"	fmadd.d	ft0, fa0, fs0, ft0\n"
	"	fld	fs5, 8(t3)\n"
	"	fmadd.d	ft1, fa1, fs0, ft1\n"
	"	fld	fs6, 16(t3)\n"
	"	fmadd.d	ft2, fa6, fs0, ft2\n"
	"	fld	fs7, 24(t3)\n"
	"	fmadd.d	ft3, fa7, fs0, ft3\n"
	"	fld	fs8, 0(t4)\n"
	"	fmadd.d	ft4, fa0, fs1, ft4\n"
	"	fld	fs9, 8(t4)\n"
	"	fmadd.d	ft5, fa1, fs1, ft5\n"
	"	fld	fs10, 16(t4)\n"
	"	fmadd.d	ft6, fa6, fs1, ft6\n"
	"	fld	fs11, 24(t4)\n"
	"	fmadd.d	ft7, fa7, fs1, ft7\n"
	"	addi	t3, t3, 32\n"
	"	fmadd.d	ft8, fa0, fs2, ft8\n"
	"	addi	t4, t4, 32\n"
	"	fmadd.d	ft9, fa1, fs2, ft9\n"
	"	fmadd.d	ft10, fa6, fs2, ft10\n"
	"	fmadd.d	ft11, fa7, fs2, ft11\n"
	"	fmadd.d	fa2, fa0, fs3, fa2\n"
	"	fmadd.d	fa3, fa1, fs3, fa3\n"
	"	fmadd.d	fa4, fa6, fs3, fa4\n"
	"	fmadd.d	fa5, fa7, fs3, fa5\n"
	"	fld	fa0, 0(t3)\n"
	"	fmadd.d	ft0, fs4, fs8, ft0\n"
	"	fld	fa1, 8(t3)\n"
	"	fmadd.d	ft1, fs5, fs8, ft1\n"
	"	fld	fa6, 16(t3)\n"
	"	fmadd.d	ft2, fs6, fs8, ft2\n"
	"	fld	fa7, 24(t3)\n"
	"	fmadd.d	ft3, fs7, fs8, ft3\n"
	"	fld	fs0, 0(t4)\n"
	"	fmadd.d	ft4, fs4, fs9, ft4\n"
	"	fld	fs1, 8(t4)\n"
	"	fmadd.d	ft5, fs5, fs9, ft5\n"
	"	fld	fs2, 16(t4)\n"
	"	fmadd.d	ft6, fs6, fs9, ft6\n"
	"	fld	fs3, 24(t4)\n"
	"	fmadd.d	ft7, fs7, fs9, ft7\n"
	"	addi	t3, t3, 32\n"
	"	fmadd.d	ft8, fs4, fs10, ft8\n"
	"	addi	t4, t4, 32\n"
	"	fmadd.d	ft9, fs5, fs10, ft9\n"
	"	fmadd.d	ft10, fs6, fs10, ft10\n"
	"	fmadd.d	ft11, fs7, fs10, ft11\n"
	"	fmadd.d	fa2, fs4, fs11, fa2\n"
	"	fmadd.d	fa3, fs5, fs11, fa3\n"
	"	fmadd.d	fa4, fs6, fs11, fa4\n"
	"	fmadd.d	fa5, fs7, fs11, fa5\n"
	"	addi	t6, t6, -1\n"
	"	bnez	t6, .Lk4body\n"
	".Lk4epi:\n"
	"	fld	fs4, 0(t3)\n"
	"	fmadd.d	ft0, fa0, fs0, ft0\n"
	"	fld	fs5, 8(t3)\n"
	"	fmadd.d	ft1, fa1, fs0, ft1\n"
	"	fld	fs6, 16(t3)\n"
	"	fmadd.d	ft2, fa6, fs0, ft2\n"
	"	fld	fs7, 24(t3)\n"
	"	fmadd.d	ft3, fa7, fs0, ft3\n"
	"	fld	fs8, 0(t4)\n"
	"	fmadd.d	ft4, fa0, fs1, ft4\n"
	"	fld	fs9, 8(t4)\n"
	"	fmadd.d	ft5, fa1, fs1, ft5\n"
	"	fld	fs10, 16(t4)\n"
	"	fmadd.d	ft6, fa6, fs1, ft6\n"
	"	fld	fs11, 24(t4)\n"
	"	fmadd.d	ft7, fa7, fs1, ft7\n"
	"	addi	t3, t3, 32\n"
	"	fmadd.d	ft8, fa0, fs2, ft8\n"
	"	addi	t4, t4, 32\n"
	"	fmadd.d	ft9, fa1, fs2, ft9\n"
	"	fmadd.d	ft10, fa6, fs2, ft10\n"
	"	fmadd.d	ft11, fa7, fs2, ft11\n"
	"	fmadd.d	fa2, fa0, fs3, fa2\n"
	"	fmadd.d	fa3, fa1, fs3, fa3\n"
	"	fmadd.d	fa4, fa6, fs3, fa4\n"
	"	fmadd.d	fa5, fa7, fs3, fa5\n"
	"	fmadd.d	ft0, fs4, fs8, ft0\n"
	"	fmadd.d	ft1, fs5, fs8, ft1\n"
	"	fmadd.d	ft2, fs6, fs8, ft2\n"
	"	fmadd.d	ft3, fs7, fs8, ft3\n"
	"	fmadd.d	ft4, fs4, fs9, ft4\n"
	"	fmadd.d	ft5, fs5, fs9, ft5\n"
	"	fmadd.d	ft6, fs6, fs9, ft6\n"
	"	fmadd.d	ft7, fs7, fs9, ft7\n"
	"	fmadd.d	ft8, fs4, fs10, ft8\n"
	"	fmadd.d	ft9, fs5, fs10, ft9\n"
	"	fmadd.d	ft10, fs6, fs10, ft10\n"
	"	fmadd.d	ft11, fs7, fs10, ft11\n"
	"	fmadd.d	fa2, fs4, fs11, fa2\n"
	"	fmadd.d	fa3, fs5, fs11, fa3\n"
	"	fmadd.d	fa4, fs6, fs11, fa4\n"
	"	fmadd.d	fa5, fs7, fs11, fa5\n"
	"	fld	fs4, 96(sp)\n"
	"	fld	fs0, 0(a0)\n"
	"	fmadd.d	fs0, ft0, fs4, fs0\n"
	"	fsd	fs0, 0(a0)\n"
	"	fld	fs1, 8(a0)\n"
	"	fmadd.d	fs1, ft1, fs4, fs1\n"
	"	fsd	fs1, 8(a0)\n"
	"	fld	fs0, 16(a0)\n"
	"	fmadd.d	fs0, ft2, fs4, fs0\n"
	"	fsd	fs0, 16(a0)\n"
	"	fld	fs1, 24(a0)\n"
	"	fmadd.d	fs1, ft3, fs4, fs1\n"
	"	fsd	fs1, 24(a0)\n"
	"	fld	fs0, 0(a1)\n"
	"	fmadd.d	fs0, ft4, fs4, fs0\n"
	"	fsd	fs0, 0(a1)\n"
	"	fld	fs1, 8(a1)\n"
	"	fmadd.d	fs1, ft5, fs4, fs1\n"
	"	fsd	fs1, 8(a1)\n"
	"	fld	fs0, 16(a1)\n"
	"	fmadd.d	fs0, ft6, fs4, fs0\n"
	"	fsd	fs0, 16(a1)\n"
	"	fld	fs1, 24(a1)\n"
	"	fmadd.d	fs1, ft7, fs4, fs1\n"
	"	fsd	fs1, 24(a1)\n"
	"	fld	fs0, 0(a6)\n"
	"	fmadd.d	fs0, ft8, fs4, fs0\n"
	"	fsd	fs0, 0(a6)\n"
	"	fld	fs1, 8(a6)\n"
	"	fmadd.d	fs1, ft9, fs4, fs1\n"
	"	fsd	fs1, 8(a6)\n"
	"	fld	fs0, 16(a6)\n"
	"	fmadd.d	fs0, ft10, fs4, fs0\n"
	"	fsd	fs0, 16(a6)\n"
	"	fld	fs1, 24(a6)\n"
	"	fmadd.d	fs1, ft11, fs4, fs1\n"
	"	fsd	fs1, 24(a6)\n"
	"	fld	fs0, 0(a7)\n"
	"	fmadd.d	fs0, fa2, fs4, fs0\n"
	"	fsd	fs0, 0(a7)\n"
	"	fld	fs1, 8(a7)\n"
	"	fmadd.d	fs1, fa3, fs4, fs1\n"
	"	fsd	fs1, 8(a7)\n"
	"	fld	fs0, 16(a7)\n"
	"	fmadd.d	fs0, fa4, fs4, fs0\n"
	"	fsd	fs0, 16(a7)\n"
	"	fld	fs1, 24(a7)\n"
	"	fmadd.d	fs1, fa5, fs4, fs1\n"
	"	fsd	fs1, 24(a7)\n"
	"	addi	a0, a0, 32\n"
	"	addi	a1, a1, 32\n"
	"	addi	a6, a6, 32\n"
	"	addi	a7, a7, 32\n"
	"	addi	t5, t5, -1\n"
	"	bnez	t5, .Lk4i\n"
	"	slli	t6, a2, 5\n"
	"	add	a4, a4, t6\n"
	"	slli	t6, t0, 2\n"
	"	add	a5, a5, t6\n"
	"	addi	t2, t2, -1\n"
	"	bnez	t2, .Lk4j\n"
	".Lk4done:\n"
	"	fld	fs0, 0(sp)\n"
	"	fld	fs1, 8(sp)\n"
	"	fld	fs2, 16(sp)\n"
	"	fld	fs3, 24(sp)\n"
	"	fld	fs4, 32(sp)\n"
	"	fld	fs5, 40(sp)\n"
	"	fld	fs6, 48(sp)\n"
	"	fld	fs7, 56(sp)\n"
	"	fld	fs8, 64(sp)\n"
	"	fld	fs9, 72(sp)\n"
	"	fld	fs10, 80(sp)\n"
	"	fld	fs11, 88(sp)\n"
	"	addi	sp, sp, 112\n"
	"	li	a0, 0\n"
	"	ret\n"
	"	.size	dgemm_u74_asm, .-dgemm_u74_asm\n"
);

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

   #if defined(DOUBLE) && !defined(TRMMKERNEL)
   if (bm > 0 && bn > 0 && bk > 0 && ((bm & 3) == 0) && ((bn & 3) == 0) && ((bk & 1) == 0))
      return dgemm_u74_asm(bm, bn, bk, alpha, ba, bb, C, ldc);
#endif

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
