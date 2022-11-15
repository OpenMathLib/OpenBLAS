/***************************************************************************
Copyright (c) 2022, The OpenBLAS Project
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

#if !defined(DOUBLE)
#define VSETVL(n) vsetvl_e32m1(n)
#define VSETVL_MAX vsetvlmax_e32m1()
#define VSETVL_MAX_M1 vsetvlmax_e32m1()
#define FLOAT_V_T vfloat32m1_t
#define FLOAT_V_T_M1 vfloat32m1_t
#define VLEV_FLOAT vle32_v_f32m1
#define VLSEG2_FLOAT vlseg2e32_v_f32m1
#define VLSEG4_FLOAT vlseg4e32_v_f32m1
#define VFMVVF_FLOAT vfmv_v_f_f32m1
#define VFMACCVF_FLOAT vfmacc_vf_f32m1
#define VFMACCVV_FLOAT vfmacc_vv_f32m1
#define VFREDSUMVS_FLOAT vfredusum_vs_f32m1_f32m1
#define VFMVVF_FLOAT_M1 vfmv_v_f_f32m1
#define VFMVFS_FLOAT_M1 vfmv_f_s_f32m1_f32
#else
#define VSETVL(n) vsetvl_e64m1(n)
#define VSETVL_MAX vsetvlmax_e64m1()
#define VSETVL_MAX_M1 vsetvlmax_e64m1()
#define FLOAT_V_T vfloat64m1_t
#define FLOAT_V_T_M1 vfloat64m1_t
#define VLEV_FLOAT vle64_v_f64m1
#define VLSEG2_FLOAT vlseg2e64_v_f64m1
#define VLSEG4_FLOAT vlseg4e64_v_f64m1
#define VFMVVF_FLOAT vfmv_v_f_f64m1
#define VFMACCVF_FLOAT vfmacc_vf_f64m1
#define VFMACCVV_FLOAT vfmacc_vv_f64m1
#define VFREDSUMVS_FLOAT vfredusum_vs_f64m1_f64m1
#define VFMVVF_FLOAT_M1 vfmv_v_f_f64m1
#define VFMVFS_FLOAT_M1 vfmv_f_s_f64m1_f64
#endif

// Optimizes the implementation in ../generic/gemm_kernel_2x2.c

int CNAME(BLASLONG bm, BLASLONG bn, BLASLONG bk, FLOAT alpha, IFLOAT* ba, IFLOAT* bb, FLOAT* C, BLASLONG ldc
#ifdef TRMMKERNEL
		,BLASLONG offset
#endif
		)
{
    BLASLONG i,j,k;
    FLOAT *C0,*C1,*C2,*C3;
    IFLOAT *ptrba,*ptrbb;

    //fprintf(stderr, "gemm_kernel_4x4 bm=%ld bn=%ld bk=%ld alpha=%f ldc=%ld\n", bm, bn, bk, alpha, ldc); // KU

    FLOAT_V_T va0, va1, va2, va3; 
    FLOAT_V_T vb0, vb1, vb2, vb3;
    FLOAT_V_T vres0, vres1, vres2, vres3, vres4, vres5, vres6, vres7;
    FLOAT_V_T vres8, vres9, vres10, vres11, vres12, vres13, vres14, vres15;
    FLOAT_V_T_M1 vsum0, vsum1, vsum2, vsum3;
    FLOAT_V_T_M1 v_z0;

    v_z0 = VFMVVF_FLOAT_M1(0, VSETVL_MAX_M1);
    size_t vlmax = VSETVL_MAX;
    size_t vl;

    for (j = bn/4; j > 0; j--) {
        C0 = C;
        C1 = C0 + ldc;
        C2 = C1 + ldc;
        C3 = C2 + ldc;
        ptrba = ba;

        for (i = bm/4; i > 0; i--) {
            ptrbb = bb;

            vres0 = VFMVVF_FLOAT(0.0, vlmax);
            vres1 = VFMVVF_FLOAT(0.0, vlmax);
            vres2 = VFMVVF_FLOAT(0.0, vlmax);
            vres3 = VFMVVF_FLOAT(0.0, vlmax);
            vres4 = VFMVVF_FLOAT(0.0, vlmax);
            vres5 = VFMVVF_FLOAT(0.0, vlmax);
            vres6 = VFMVVF_FLOAT(0.0, vlmax);
            vres7 = VFMVVF_FLOAT(0.0, vlmax);
            vres8 = VFMVVF_FLOAT(0.0, vlmax);
            vres9 = VFMVVF_FLOAT(0.0, vlmax);
            vres10 = VFMVVF_FLOAT(0.0, vlmax);
            vres11 = VFMVVF_FLOAT(0.0, vlmax);
            vres12 = VFMVVF_FLOAT(0.0, vlmax);
            vres13 = VFMVVF_FLOAT(0.0, vlmax);
            vres14 = VFMVVF_FLOAT(0.0, vlmax);
            vres15 = VFMVVF_FLOAT(0.0, vlmax);

            for (k = bk; k > 0; k -= vl) {
                vl = VSETVL(k);

                VLSEG4_FLOAT(&va0, &va1, &va2, &va3, ptrba, vl); 
                VLSEG4_FLOAT(&vb0, &vb1, &vb2, &vb3, ptrbb, vl); 

                vres0 = VFMACCVV_FLOAT(vres0, va0, vb0, vl);
                vres1 = VFMACCVV_FLOAT(vres1, va1, vb0, vl);
                vres2 = VFMACCVV_FLOAT(vres2, va0, vb1, vl);
                vres3 = VFMACCVV_FLOAT(vres3, va1, vb1, vl);

                vres4 = VFMACCVV_FLOAT(vres4, va0, vb2, vl);
                vres5 = VFMACCVV_FLOAT(vres5, va1, vb2, vl);
                vres6 = VFMACCVV_FLOAT(vres6, va0, vb3, vl);
                vres7 = VFMACCVV_FLOAT(vres7, va1, vb3, vl);

                vres8 = VFMACCVV_FLOAT(vres8, va2, vb0, vl);
                vres9 = VFMACCVV_FLOAT(vres9, va3, vb0, vl);
                vres10 = VFMACCVV_FLOAT(vres10, va2, vb1, vl);
                vres11 = VFMACCVV_FLOAT(vres11, va3, vb1, vl);

                vres12 = VFMACCVV_FLOAT(vres12, va2, vb2, vl);
                vres13 = VFMACCVV_FLOAT(vres13, va3, vb2, vl);
                vres14 = VFMACCVV_FLOAT(vres14, va2, vb3, vl);
                vres15 = VFMACCVV_FLOAT(vres15, va3, vb3, vl);

                ptrba += vl*4;
                ptrbb += vl*4;
            }

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres1, v_z0, vlmax);
            vsum2 = VFREDSUMVS_FLOAT(vsum2, vres8, v_z0, vlmax);
            vsum3 = VFREDSUMVS_FLOAT(vsum3, vres9, v_z0, vlmax);
            C0[0] += alpha * VFMVFS_FLOAT_M1(vsum0);
            C0[1] += alpha * VFMVFS_FLOAT_M1(vsum1);
            C0[2] += alpha * VFMVFS_FLOAT_M1(vsum2);
            C0[3] += alpha * VFMVFS_FLOAT_M1(vsum3);

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres2, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres3, v_z0, vlmax);
            vsum2 = VFREDSUMVS_FLOAT(vsum2, vres10, v_z0, vlmax);
            vsum3 = VFREDSUMVS_FLOAT(vsum3, vres11, v_z0, vlmax);
            C1[0] += alpha * VFMVFS_FLOAT_M1(vsum0);
            C1[1] += alpha * VFMVFS_FLOAT_M1(vsum1);
            C1[2] += alpha * VFMVFS_FLOAT_M1(vsum2);
            C1[3] += alpha * VFMVFS_FLOAT_M1(vsum3);

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres4, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres5, v_z0, vlmax);
            vsum2 = VFREDSUMVS_FLOAT(vsum2, vres12, v_z0, vlmax);
            vsum3 = VFREDSUMVS_FLOAT(vsum3, vres13, v_z0, vlmax);
            C2[0] += alpha * VFMVFS_FLOAT_M1(vsum0);
            C2[1] += alpha * VFMVFS_FLOAT_M1(vsum1);
            C2[2] += alpha * VFMVFS_FLOAT_M1(vsum2);
            C2[3] += alpha * VFMVFS_FLOAT_M1(vsum3);

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres6, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres7, v_z0, vlmax);
            vsum2 = VFREDSUMVS_FLOAT(vsum2, vres14, v_z0, vlmax);
            vsum3 = VFREDSUMVS_FLOAT(vsum3, vres15, v_z0, vlmax);
            C3[0] += alpha * VFMVFS_FLOAT_M1(vsum0);
            C3[1] += alpha * VFMVFS_FLOAT_M1(vsum1);
            C3[2] += alpha * VFMVFS_FLOAT_M1(vsum2);
            C3[3] += alpha * VFMVFS_FLOAT_M1(vsum3);

            C0 += 4;
            C1 += 4;
            C2 += 4;
            C3 += 4;
        }

        if(bm & 2) {
            ptrbb = bb;

            vres0 = VFMVVF_FLOAT(0.0, vlmax);
            vres1 = VFMVVF_FLOAT(0.0, vlmax);
            vres2 = VFMVVF_FLOAT(0.0, vlmax);
            vres3 = VFMVVF_FLOAT(0.0, vlmax);
            vres4 = VFMVVF_FLOAT(0.0, vlmax);
            vres5 = VFMVVF_FLOAT(0.0, vlmax);
            vres6 = VFMVVF_FLOAT(0.0, vlmax);
            vres7 = VFMVVF_FLOAT(0.0, vlmax);

            for (k = bk; k > 0; k -= vl) {
                vl = VSETVL(k);

                VLSEG2_FLOAT(&va0, &va1, ptrba, vl); 
                VLSEG4_FLOAT(&vb0, &vb1, &vb2, &vb3, ptrbb, vl); 

                vres0 = VFMACCVV_FLOAT(vres0, va0, vb0, vl);
                vres1 = VFMACCVV_FLOAT(vres1, va1, vb0, vl);
                vres2 = VFMACCVV_FLOAT(vres2, va0, vb1, vl);
                vres3 = VFMACCVV_FLOAT(vres3, va1, vb1, vl);

                vres4 = VFMACCVV_FLOAT(vres4, va0, vb2, vl);
                vres5 = VFMACCVV_FLOAT(vres5, va1, vb2, vl);
                vres6 = VFMACCVV_FLOAT(vres6, va0, vb3, vl);
                vres7 = VFMACCVV_FLOAT(vres7, va1, vb3, vl);

                ptrba += vl*2;
                ptrbb += vl*4;
            }

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres1, v_z0, vlmax);
            C0[0] += alpha * VFMVFS_FLOAT_M1(vsum0);
            C0[1] += alpha * VFMVFS_FLOAT_M1(vsum1);
 
            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres2, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres3, v_z0, vlmax);
            C1[0] += alpha * VFMVFS_FLOAT_M1(vsum0);
            C1[1] += alpha * VFMVFS_FLOAT_M1(vsum1);

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres4, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres5, v_z0, vlmax);
            C2[0] += alpha * VFMVFS_FLOAT_M1(vsum0);
            C2[1] += alpha * VFMVFS_FLOAT_M1(vsum1);

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres6, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres7, v_z0, vlmax);
            C3[0] += alpha * VFMVFS_FLOAT_M1(vsum0);
            C3[1] += alpha * VFMVFS_FLOAT_M1(vsum1);

            C0 += 2;
            C1 += 2;
            C2 += 2;
            C3 += 2;
        }

        if(bm & 1) {
            ptrbb = bb;

            vres0 = VFMVVF_FLOAT(0.0, vlmax);
            vres1 = VFMVVF_FLOAT(0.0, vlmax);
            vres2 = VFMVVF_FLOAT(0.0, vlmax);
            vres3 = VFMVVF_FLOAT(0.0, vlmax);

            for (k = bk; k > 0; k -= vl) {
                vl = VSETVL(k);

                va0 = VLEV_FLOAT(ptrba, vl);
                VLSEG4_FLOAT(&vb0, &vb1, &vb2, &vb3, ptrbb, vl); 

                vres0 = VFMACCVV_FLOAT(vres0, va0, vb0, vl);
                vres1 = VFMACCVV_FLOAT(vres1, va0, vb1, vl);
                vres2 = VFMACCVV_FLOAT(vres2, va0, vb2, vl);
                vres3 = VFMACCVV_FLOAT(vres3, va0, vb3, vl);

                ptrba += vl;
                ptrbb += vl*4;
            }

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres1, v_z0, vlmax);
            vsum2 = VFREDSUMVS_FLOAT(vsum2, vres2, v_z0, vlmax);
            vsum3 = VFREDSUMVS_FLOAT(vsum3, vres3, v_z0, vlmax);
            C0[0] += alpha * VFMVFS_FLOAT_M1(vsum0);
            C1[0] += alpha * VFMVFS_FLOAT_M1(vsum1);
            C2[0] += alpha * VFMVFS_FLOAT_M1(vsum2);
            C3[0] += alpha * VFMVFS_FLOAT_M1(vsum3);

            C0 += 1;
            C1 += 1;
            C2 += 1;
            C3 += 1;
        }

        bb += (bk<<2);
        C += (ldc<<2);
    }

    if(bn & 2) {

        C0 = C;
        C1 = C0 + ldc;
        ptrba = ba;

        for (i = bm/4; i > 0; i--) {
            ptrbb = bb;

            vres0 = VFMVVF_FLOAT(0.0, vlmax);
            vres1 = VFMVVF_FLOAT(0.0, vlmax);
            vres2 = VFMVVF_FLOAT(0.0, vlmax);
            vres3 = VFMVVF_FLOAT(0.0, vlmax);

            vres4 = VFMVVF_FLOAT(0.0, vlmax);
            vres5 = VFMVVF_FLOAT(0.0, vlmax);
            vres6 = VFMVVF_FLOAT(0.0, vlmax);
            vres7 = VFMVVF_FLOAT(0.0, vlmax);

            for (k = bk; k > 0; k -= vl) {
                vl = VSETVL(k);

                VLSEG4_FLOAT(&va0, &va1, &va2, &va3, ptrba, vl);
                VLSEG2_FLOAT(&vb0, &vb1, ptrbb, vl);

                vres0 = VFMACCVV_FLOAT(vres0, va0, vb0, vl);
                vres1 = VFMACCVV_FLOAT(vres1, va1, vb0, vl);
                vres2 = VFMACCVV_FLOAT(vres2, va2, vb0, vl);
                vres3 = VFMACCVV_FLOAT(vres3, va3, vb0, vl);

                vres4 = VFMACCVV_FLOAT(vres4, va0, vb1, vl);
                vres5 = VFMACCVV_FLOAT(vres5, va1, vb1, vl);
                vres6 = VFMACCVV_FLOAT(vres6, va2, vb1, vl);
                vres7 = VFMACCVV_FLOAT(vres7, va3, vb1, vl);

                ptrba += vl*4;
                ptrbb += vl*2;
            }

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres1, v_z0, vlmax);
            vsum2 = VFREDSUMVS_FLOAT(vsum2, vres2, v_z0, vlmax);
            vsum3 = VFREDSUMVS_FLOAT(vsum3, vres3, v_z0, vlmax);
            C0[0] += alpha * VFMVFS_FLOAT_M1(vsum0);
            C0[1] += alpha * VFMVFS_FLOAT_M1(vsum1);
            C0[2] += alpha * VFMVFS_FLOAT_M1(vsum2);
            C0[3] += alpha * VFMVFS_FLOAT_M1(vsum3);

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres4, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres5, v_z0, vlmax);
            vsum2 = VFREDSUMVS_FLOAT(vsum2, vres6, v_z0, vlmax);
            vsum3 = VFREDSUMVS_FLOAT(vsum3, vres7, v_z0, vlmax);
            C1[0] += alpha * VFMVFS_FLOAT_M1(vsum0);
            C1[1] += alpha * VFMVFS_FLOAT_M1(vsum1);
            C1[2] += alpha * VFMVFS_FLOAT_M1(vsum2);
            C1[3] += alpha * VFMVFS_FLOAT_M1(vsum3);

            C0 += 4;
            C1 += 4;
        }

        if(bm & 2) {
            ptrbb = bb;

            vres0 = VFMVVF_FLOAT(0.0, vlmax);
            vres1 = VFMVVF_FLOAT(0.0, vlmax);
            vres2 = VFMVVF_FLOAT(0.0, vlmax);
            vres3 = VFMVVF_FLOAT(0.0, vlmax);

            for (k = bk; k > 0; k -= vl) {
                vl = VSETVL(k);

                VLSEG2_FLOAT(&va0, &va1, ptrba, vl);
                VLSEG2_FLOAT(&vb0, &vb1, ptrbb, vl);

                vres0 = VFMACCVV_FLOAT(vres0, va0, vb0, vl);
                vres1 = VFMACCVV_FLOAT(vres1, va1, vb0, vl);
                vres2 = VFMACCVV_FLOAT(vres2, va0, vb1, vl);
                vres3 = VFMACCVV_FLOAT(vres3, va1, vb1, vl);

                ptrba += vl*2;
                ptrbb += vl*2;
            }

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres1, v_z0, vlmax);
            vsum2 = VFREDSUMVS_FLOAT(vsum2, vres2, v_z0, vlmax);
            vsum3 = VFREDSUMVS_FLOAT(vsum3, vres3, v_z0, vlmax);
            C0[0] += alpha * VFMVFS_FLOAT_M1(vsum0);
            C0[1] += alpha * VFMVFS_FLOAT_M1(vsum1);
            C1[0] += alpha * VFMVFS_FLOAT_M1(vsum2);
            C1[1] += alpha * VFMVFS_FLOAT_M1(vsum3);

            C0 += 2;
            C1 += 2;
        }

        if(bm & 1) {
            ptrbb = bb;

            vres0 = VFMVVF_FLOAT(0.0, vlmax);
            vres1 = VFMVVF_FLOAT(0.0, vlmax);

            for (k = bk; k > 0; k -= vl) {
                vl = VSETVL(k);

                va0 = VLEV_FLOAT(ptrba, vl);
                VLSEG2_FLOAT(&vb0, &vb1, ptrbb, vl);

                vres0 = VFMACCVV_FLOAT(vres0, va0, vb0, vl);
                vres1 = VFMACCVV_FLOAT(vres1, va0, vb1, vl);

                ptrba += vl;
                ptrbb += vl*2;
            }

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres1, v_z0, vlmax);
            C0[0] += alpha * VFMVFS_FLOAT_M1(vsum0);
            C1[0] += alpha * VFMVFS_FLOAT_M1(vsum1);

            C0 += 1;
            C1 += 1;
        }

        bb += (bk<<1);
        C += (ldc<<1);
    }

    if(bn & 1) {
        C0 = C;
        ptrba = ba;
        for (i = bm/4; i > 0; i--) {
            ptrbb = bb;

            vres0 = VFMVVF_FLOAT(0.0, vlmax);
            vres1 = VFMVVF_FLOAT(0.0, vlmax);
            vres2 = VFMVVF_FLOAT(0.0, vlmax);
            vres3 = VFMVVF_FLOAT(0.0, vlmax);

            for (k = bk; k > 0; k -= vl) {
                vl = VSETVL(k);

                VLSEG4_FLOAT(&va0, &va1, &va2, &va3, ptrba, vl);
                vb0 = VLEV_FLOAT(ptrbb, vl);

                vres0 = VFMACCVV_FLOAT(vres0, va0, vb0, vl);
                vres1 = VFMACCVV_FLOAT(vres1, va1, vb0, vl);
                vres2 = VFMACCVV_FLOAT(vres2, va2, vb0, vl);
                vres3 = VFMACCVV_FLOAT(vres3, va3, vb0, vl);

                ptrba += vl*4;
                ptrbb += vl;
            }

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres1, v_z0, vlmax);
            vsum2 = VFREDSUMVS_FLOAT(vsum2, vres2, v_z0, vlmax);
            vsum3 = VFREDSUMVS_FLOAT(vsum3, vres3, v_z0, vlmax);
            C0[0] += alpha * VFMVFS_FLOAT_M1(vsum0);
            C0[1] += alpha * VFMVFS_FLOAT_M1(vsum1);
            C0[2] += alpha * VFMVFS_FLOAT_M1(vsum2);
            C0[3] += alpha * VFMVFS_FLOAT_M1(vsum3);

            C0 += 4;
        }

        if(bm & 2) {
            ptrbb = bb;

            vres0 = VFMVVF_FLOAT(0.0, vlmax);
            vres1 = VFMVVF_FLOAT(0.0, vlmax);

            for (k = bk; k > 0; k -= vl) {
                vl = VSETVL(k);

                VLSEG2_FLOAT(&va0, &va1, ptrba, vl);
                vb0 = VLEV_FLOAT(ptrbb, vl);

                vres0 = VFMACCVV_FLOAT(vres0, va0, vb0, vl);
                vres1 = VFMACCVV_FLOAT(vres1, va1, vb0, vl);

                ptrba += vl*2;
                ptrbb += vl;
            }

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres1, v_z0, vlmax);
            C0[0] += alpha * VFMVFS_FLOAT_M1(vsum0);
            C0[1] += alpha * VFMVFS_FLOAT_M1(vsum1);

            C0 += 2;
        }

        if(bm & 1) {
            ptrbb = bb;

            vres0 = VFMVVF_FLOAT(0.0, vlmax);

            for (k = bk; k > 0; k -= vl) {
                vl = VSETVL(k);

                va0 = VLEV_FLOAT(ptrba, vl);
                vb0 = VLEV_FLOAT(ptrbb, vl);

                vres0 = VFMACCVV_FLOAT(vres0, va0, vb0, vl);

                ptrba += vl;
                ptrbb += vl;
            }

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres0, v_z0, vlmax);
            C0[0] += alpha * VFMVFS_FLOAT_M1(vsum0);

            C0 += 1;
        }

        bb += (bk<<0);
        C += ldc;
    }

    return 0;
}
