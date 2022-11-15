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
#include <stdbool.h>

#if !defined(DOUBLE)
#define VSETVL(n)           vsetvl_e32m2(n)
#define VSETVL_MAX          vsetvlmax_e32m2()
#define VSETVL_MAX_M1       vsetvlmax_e32m1()
#define FLOAT_V_T           vfloat32m2_t
#define FLOAT_V_T_M1        vfloat32m1_t
#define VLEV_FLOAT          vle32_v_f32m2
#define VLSEG4_FLOAT        vlseg4e32_v_f32m2
#define VLSEG2_FLOAT        vlseg2e32_v_f32m2
#define VFMVVF_FLOAT        vfmv_v_f_f32m2
#define VFMUL_FLOAT         vfmul_vv_f32m2
#define VFMACCVF_FLOAT      vfmacc_vf_f32m2
#define VFMACCVV_FLOAT      vfmacc_vv_f32m2
#define VFREDSUMVS_FLOAT    vfredusum_vs_f32m2_f32m1
#define VFMVVF_FLOAT_M1     vfmv_v_f_f32m1
#define VFMVFS_FLOAT_M1     vfmv_f_s_f32m1_f32
#else
#define VSETVL(n)           vsetvl_e64m2(n)
#define VSETVL_MAX          vsetvlmax_e64m2()
#define VSETVL_MAX_M1       vsetvlmax_e64m1()
#define FLOAT_V_T           vfloat64m2_t
#define FLOAT_V_T_M1        vfloat64m1_t
#define VLEV_FLOAT          vle64_v_f64m2
#define VLSEG4_FLOAT        vlseg4e64_v_f64m2
#define VLSEG2_FLOAT        vlseg2e64_v_f64m2
#define VFMVVF_FLOAT        vfmv_v_f_f64m2
#define VFMUL_FLOAT         vfmul_vv_f64m2
#define VFMACCVF_FLOAT      vfmacc_vf_f64m2
#define VFMACCVV_FLOAT      vfmacc_vv_f64m2
#define VFREDSUMVS_FLOAT    vfredusum_vs_f64m2_f64m1
#define VFMVVF_FLOAT_M1     vfmv_v_f_f64m1
#define VFMVFS_FLOAT_M1     vfmv_f_s_f64m1_f64
#endif


// Optimizes the implementation in ../generic/trmmkernel_4x4.c

int CNAME(BLASLONG bm,BLASLONG bn,BLASLONG bk,FLOAT alpha,FLOAT* ba,FLOAT* bb,FLOAT* C,BLASLONG ldc ,BLASLONG offset)
{

   BLASLONG i,j,k;
   FLOAT *C0,*C1,*C2,*C3,*ptrba,*ptrbb;

   FLOAT_V_T va0, va1, va2, va3, vb0, vb1, vb2, vb3;
   FLOAT_V_T_M1 vsum0, vsum1, vsum2, vsum3, v_z0;
   v_z0 = VFMVVF_FLOAT_M1(0, VSETVL_MAX_M1);
   size_t vl;
   size_t vlmax = VSETVL_MAX;

   FLOAT_V_T vres0_0;
   FLOAT_V_T vres0_1;
   FLOAT_V_T vres0_2;
   FLOAT_V_T vres0_3;

   FLOAT_V_T vres1_0;
   FLOAT_V_T vres1_1;
   FLOAT_V_T vres1_2;
   FLOAT_V_T vres1_3;

   FLOAT_V_T vres2_0;
   FLOAT_V_T vres2_1;
   FLOAT_V_T vres2_2;
   FLOAT_V_T vres2_3;

   FLOAT_V_T vres3_0;
   FLOAT_V_T vres3_1;
   FLOAT_V_T vres3_2;
   FLOAT_V_T vres3_3;

   BLASLONG off, temp;

   bool left;
   bool transposed;
   bool backwards;

#ifdef LEFT
   left = true;
#else
   left = false;
#endif

#ifdef TRANSA
   transposed = true;
#else
   transposed = false;
#endif

   backwards = left != transposed;

   if (!left) {
      off = -offset;
   }


   for (j=0; j<bn/4; j+=1) // do blocks of the Mx4 loops 
   {
        C0 = C;
        C1 = C0+ldc;
        C2 = C1+ldc;
        C3 = C2+ldc;


        if (left) {
            off = offset;
        }

        ptrba = ba;

        for (i=0; i<bm/4; i+=1) // do blocks of 4x4
        {

            ptrbb = bb;
            if (backwards)
            {
                ptrba += off*4; // number of values in A
                ptrbb += off*4; // number of values in B
            }

            vres0_0 = VFMVVF_FLOAT(0, vlmax);
            vres0_1 = VFMVVF_FLOAT(0, vlmax);
            vres0_2 = VFMVVF_FLOAT(0, vlmax);
            vres0_3 = VFMVVF_FLOAT(0, vlmax);
            
            vres1_0 = VFMVVF_FLOAT(0, vlmax);
            vres1_1 = VFMVVF_FLOAT(0, vlmax);
            vres1_2 = VFMVVF_FLOAT(0, vlmax);
            vres1_3 = VFMVVF_FLOAT(0, vlmax);
            
            vres2_0 = VFMVVF_FLOAT(0, vlmax);
            vres2_1 = VFMVVF_FLOAT(0, vlmax);
            vres2_2 = VFMVVF_FLOAT(0, vlmax);
            vres2_3 = VFMVVF_FLOAT(0, vlmax);
            
            vres3_0 = VFMVVF_FLOAT(0, vlmax);
            vres3_1 = VFMVVF_FLOAT(0, vlmax);
            vres3_2 = VFMVVF_FLOAT(0, vlmax);
            vres3_3 = VFMVVF_FLOAT(0, vlmax);

            temp = backwards ? bk-off :
                         left ? off + 4 : // number of values in A
                                off + 4;  // number of values in B

            for (k = temp; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                VLSEG4_FLOAT(&va0, &va1, &va2, &va3, ptrba, vl);
                VLSEG4_FLOAT(&vb0, &vb1, &vb2, &vb3, ptrbb, vl);
            
                vres0_0 = VFMACCVV_FLOAT(vres0_0, va0, vb0, vl);
                vres1_0 = VFMACCVV_FLOAT(vres1_0, va0, vb1, vl);
                vres2_0 = VFMACCVV_FLOAT(vres2_0, va0, vb2, vl);
                vres3_0 = VFMACCVV_FLOAT(vres3_0, va0, vb3, vl);
            
                vres0_1 = VFMACCVV_FLOAT(vres0_1, va1, vb0, vl);
                vres1_1 = VFMACCVV_FLOAT(vres1_1, va1, vb1, vl);
                vres2_1 = VFMACCVV_FLOAT(vres2_1, va1, vb2, vl);
                vres3_1 = VFMACCVV_FLOAT(vres3_1, va1, vb3, vl);

                vres0_2 = VFMACCVV_FLOAT(vres0_2, va2, vb0, vl);
                vres1_2 = VFMACCVV_FLOAT(vres1_2, va2, vb1, vl);
                vres2_2 = VFMACCVV_FLOAT(vres2_2, va2, vb2, vl);
                vres3_2 = VFMACCVV_FLOAT(vres3_2, va2, vb3, vl);

                vres0_3 = VFMACCVV_FLOAT(vres0_3, va3, vb0, vl);
                vres1_3 = VFMACCVV_FLOAT(vres1_3, va3, vb1, vl);
                vres2_3 = VFMACCVV_FLOAT(vres2_3, va3, vb2, vl);
                vres3_3 = VFMACCVV_FLOAT(vres3_3, va3, vb3, vl);

                ptrba += vl * 4;
                ptrbb += vl * 4;
            }

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres0_0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres0_1, v_z0, vlmax);
            vsum2 = VFREDSUMVS_FLOAT(vsum2, vres0_2, v_z0, vlmax);
            vsum3 = VFREDSUMVS_FLOAT(vsum3, vres0_3, v_z0, vlmax);
            C0[0] = alpha * VFMVFS_FLOAT_M1(vsum0);
            C0[1] = alpha * VFMVFS_FLOAT_M1(vsum1);
            C0[2] = alpha * VFMVFS_FLOAT_M1(vsum2);
            C0[3] = alpha * VFMVFS_FLOAT_M1(vsum3);

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres1_0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres1_1, v_z0, vlmax);
            vsum2 = VFREDSUMVS_FLOAT(vsum2, vres1_2, v_z0, vlmax);
            vsum3 = VFREDSUMVS_FLOAT(vsum3, vres1_3, v_z0, vlmax);
            C1[0] = alpha * VFMVFS_FLOAT_M1(vsum0);
            C1[1] = alpha * VFMVFS_FLOAT_M1(vsum1);
            C1[2] = alpha * VFMVFS_FLOAT_M1(vsum2);
            C1[3] = alpha * VFMVFS_FLOAT_M1(vsum3);

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres2_0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres2_1, v_z0, vlmax);
            vsum2 = VFREDSUMVS_FLOAT(vsum2, vres2_2, v_z0, vlmax);
            vsum3 = VFREDSUMVS_FLOAT(vsum3, vres2_3, v_z0, vlmax);
            C2[0] = alpha * VFMVFS_FLOAT_M1(vsum0);
            C2[1] = alpha * VFMVFS_FLOAT_M1(vsum1);
            C2[2] = alpha * VFMVFS_FLOAT_M1(vsum2);
            C2[3] = alpha * VFMVFS_FLOAT_M1(vsum3);

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres3_0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres3_1, v_z0, vlmax);
            vsum2 = VFREDSUMVS_FLOAT(vsum2, vres3_2, v_z0, vlmax);
            vsum3 = VFREDSUMVS_FLOAT(vsum3, vres3_3, v_z0, vlmax);
            C3[0] = alpha * VFMVFS_FLOAT_M1(vsum0);
            C3[1] = alpha * VFMVFS_FLOAT_M1(vsum1);
            C3[2] = alpha * VFMVFS_FLOAT_M1(vsum2);
            C3[3] = alpha * VFMVFS_FLOAT_M1(vsum3);

            if (!backwards) {
                temp = bk-off;
                temp = left ? temp - 4 : // number of values in A
                              temp - 4;  // number of values in B

                ptrba += temp*4; // number of values in A
                ptrbb += temp*4; // number of values in B
            }
#ifdef LEFT
            off += 4; // number of values in A
#endif

            C0 = C0+4;
            C1 = C1+4;
            C2 = C2+4;
            C3 = C3+4;

        }

        if ( bm & 2 ) // do any 2x4 loop
        {

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
            ptrbb = bb;
#else
            ptrba += off*2;
            ptrbb = bb + off*4;
#endif

            vres0_0 = VFMVVF_FLOAT(0, vlmax);
            vres0_1 = VFMVVF_FLOAT(0, vlmax);
            
            vres1_0 = VFMVVF_FLOAT(0, vlmax);
            vres1_1 = VFMVVF_FLOAT(0, vlmax);
            
            vres2_0 = VFMVVF_FLOAT(0, vlmax);
            vres2_1 = VFMVVF_FLOAT(0, vlmax);
            
            vres3_0 = VFMVVF_FLOAT(0, vlmax);
            vres3_1 = VFMVVF_FLOAT(0, vlmax);


#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
            temp = bk-off;
#elif defined(LEFT)
            temp = off+2;   // number of values in A
#else
            temp = off+4;   // number of values in B
#endif
            for (k = temp; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                VLSEG2_FLOAT(&va0, &va1, ptrba, vl);
                VLSEG4_FLOAT(&vb0, &vb1, &vb2, &vb3, ptrbb, vl);
            
                vres0_0 = VFMACCVV_FLOAT(vres0_0, va0, vb0, vl);
                vres1_0 = VFMACCVV_FLOAT(vres1_0, va0, vb1, vl);
                vres2_0 = VFMACCVV_FLOAT(vres2_0, va0, vb2, vl);
                vres3_0 = VFMACCVV_FLOAT(vres3_0, va0, vb3, vl);

                vres0_1 = VFMACCVV_FLOAT(vres0_1, va1, vb0, vl);
                vres1_1 = VFMACCVV_FLOAT(vres1_1, va1, vb1, vl);
                vres2_1 = VFMACCVV_FLOAT(vres2_1, va1, vb2, vl);
                vres3_1 = VFMACCVV_FLOAT(vres3_1, va1, vb3, vl);

                ptrba += vl * 2;
                ptrbb += vl * 4;
            }

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres0_0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres0_1, v_z0, vlmax);
            vsum2 = VFREDSUMVS_FLOAT(vsum2, vres1_0, v_z0, vlmax);
            vsum3 = VFREDSUMVS_FLOAT(vsum3, vres1_1, v_z0, vlmax);

            C0[0] = alpha * VFMVFS_FLOAT_M1(vsum0);
            C0[1] = alpha * VFMVFS_FLOAT_M1(vsum1);
            C1[0] = alpha * VFMVFS_FLOAT_M1(vsum2);
            C1[1] = alpha * VFMVFS_FLOAT_M1(vsum3);
            
            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres2_0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres2_1, v_z0, vlmax);
            vsum2 = VFREDSUMVS_FLOAT(vsum2, vres3_0, v_z0, vlmax);
            vsum3 = VFREDSUMVS_FLOAT(vsum3, vres3_1, v_z0, vlmax);

            C2[0] = alpha * VFMVFS_FLOAT_M1(vsum0);
            C2[1] = alpha * VFMVFS_FLOAT_M1(vsum1);
            C3[0] = alpha * VFMVFS_FLOAT_M1(vsum2);
            C3[1] = alpha * VFMVFS_FLOAT_M1(vsum3);


#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
            temp = bk - off;
#ifdef LEFT
            temp -= 2; // number of values in A
#else
            temp -= 4; // number of values in B
#endif
            ptrba += temp*2;
            ptrbb += temp*4;
#endif

#ifdef LEFT
            off += 2; // number of values in A
#endif

            C0 = C0+2;
            C1 = C1+2;
            C2 = C2+2;
            C3 = C3+2;

        }

        if ( bm & 1 ) // do any 1x4 loop
        {

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
            ptrbb = bb;
#else
            ptrba += off*1;
            ptrbb = bb + off*4;
#endif

            vres0_0 = VFMVVF_FLOAT(0, vlmax);
            vres1_0 = VFMVVF_FLOAT(0, vlmax);
            vres2_0 = VFMVVF_FLOAT(0, vlmax);
            vres3_0 = VFMVVF_FLOAT(0, vlmax);


#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
            temp = bk-off;
#elif defined(LEFT)
            temp = off+1;   // number of values in A
#else
            temp = off+4;   // number of values in B
#endif

            for (k = temp; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                va0 = VLEV_FLOAT(ptrba, vl);
                VLSEG4_FLOAT(&vb0, &vb1, &vb2, &vb3, ptrbb, vl);
            
                vres0_0 = VFMACCVV_FLOAT(vres0_0, va0, vb0, vl);
                vres1_0 = VFMACCVV_FLOAT(vres1_0, va0, vb1, vl);
                vres2_0 = VFMACCVV_FLOAT(vres2_0, va0, vb2, vl);
                vres3_0 = VFMACCVV_FLOAT(vres3_0, va0, vb3, vl);

                ptrba += vl;
                ptrbb += vl * 4;
            }

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres0_0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres1_0, v_z0, vlmax);
            vsum2 = VFREDSUMVS_FLOAT(vsum2, vres2_0, v_z0, vlmax);
            vsum3 = VFREDSUMVS_FLOAT(vsum3, vres3_0, v_z0, vlmax);

            C0[0] = alpha * VFMVFS_FLOAT_M1(vsum0);
            C1[0] = alpha * VFMVFS_FLOAT_M1(vsum1);
            C2[0] = alpha * VFMVFS_FLOAT_M1(vsum2);
            C3[0] = alpha * VFMVFS_FLOAT_M1(vsum3);

#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
            temp = bk - off;
#ifdef LEFT
            temp -= 1; // number of values in A
#else
            temp -= 4; // number of values in B
#endif
            ptrba += temp*1;
            ptrbb += temp*4;
#endif

#ifdef LEFT
            off += 1; // number of values in A
#endif

            C0 = C0+1;
            C1 = C1+1;
            C2 = C2+1;
            C3 = C3+1;

        }


#if defined(TRMMKERNEL) && !defined(LEFT)
        off += 4;
#endif

        k = (bk<<2);
        bb = bb+k;
        i = (ldc<<2);
        C = C+i;
    }

    for (j=0; j<(bn&2); j+=2) // do the Mx2 loops 
    {
        C0 = C;
        C1 = C0+ldc;

#if defined(TRMMKERNEL) && defined(LEFT)
        off = offset;
#endif

        ptrba = ba;

        for (i=0; i<bm/4; i+=1) // do blocks of 4x2
        {

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
            ptrbb = bb;
#else
            ptrba += off*4;
            ptrbb = bb + off*2;
#endif
            vres0_0 = VFMVVF_FLOAT(0, vlmax);
            vres0_1 = VFMVVF_FLOAT(0, vlmax);
            vres0_2 = VFMVVF_FLOAT(0, vlmax);
            vres0_3 = VFMVVF_FLOAT(0, vlmax);
            
            vres1_0 = VFMVVF_FLOAT(0, vlmax);
            vres1_1 = VFMVVF_FLOAT(0, vlmax);
            vres1_2 = VFMVVF_FLOAT(0, vlmax);
            vres1_3 = VFMVVF_FLOAT(0, vlmax);
        

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
            temp = bk-off;
#elif defined(LEFT)
            temp = off+4;   // number of values in A
#else
            temp = off+2;   // number of values in B
#endif

            for (k = temp; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                VLSEG4_FLOAT(&va0, &va1, &va2, &va3, ptrba, vl);
                VLSEG2_FLOAT(&vb0, &vb1, ptrbb, vl);
            
                vres0_0 = VFMACCVV_FLOAT(vres0_0, va0, vb0, vl);
                vres1_0 = VFMACCVV_FLOAT(vres1_0, va0, vb1, vl);
            
                vres0_1 = VFMACCVV_FLOAT(vres0_1, va1, vb0, vl);
                vres1_1 = VFMACCVV_FLOAT(vres1_1, va1, vb1, vl);
            
                vres0_2 = VFMACCVV_FLOAT(vres0_2, va2, vb0, vl);
                vres1_2 = VFMACCVV_FLOAT(vres1_2, va2, vb1, vl);
            
                vres0_3 = VFMACCVV_FLOAT(vres0_3, va3, vb0, vl);
                vres1_3 = VFMACCVV_FLOAT(vres1_3, va3, vb1, vl);
            
                ptrba += vl * 4;
                ptrbb += vl * 2;
            }

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres0_0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres0_1, v_z0, vlmax);
            vsum2 = VFREDSUMVS_FLOAT(vsum2, vres0_2, v_z0, vlmax);
            vsum3 = VFREDSUMVS_FLOAT(vsum3, vres0_3, v_z0, vlmax);
            C0[0] = alpha * VFMVFS_FLOAT_M1(vsum0);
            C0[1] = alpha * VFMVFS_FLOAT_M1(vsum1);
            C0[2] = alpha * VFMVFS_FLOAT_M1(vsum2);
            C0[3] = alpha * VFMVFS_FLOAT_M1(vsum3);

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres1_0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres1_1, v_z0, vlmax);
            vsum2 = VFREDSUMVS_FLOAT(vsum2, vres1_2, v_z0, vlmax);
            vsum3 = VFREDSUMVS_FLOAT(vsum3, vres1_3, v_z0, vlmax);
            C1[0] = alpha * VFMVFS_FLOAT_M1(vsum0);
            C1[1] = alpha * VFMVFS_FLOAT_M1(vsum1);
            C1[2] = alpha * VFMVFS_FLOAT_M1(vsum2);
            C1[3] = alpha * VFMVFS_FLOAT_M1(vsum3);

#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
            temp = bk - off;
#ifdef LEFT
            temp -= 4; // number of values in A
#else
            temp -= 2; // number of values in B
#endif
            ptrba += temp*4;
            ptrbb += temp*2;
#endif

#ifdef LEFT
            off += 4; // number of values in A
#endif

            C0 = C0+4;
            C1 = C1+4;

        }

        if ( bm & 2 ) // do any 2x2 loop
        {

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
            ptrbb = bb;
#else
            ptrba += off*2;
            ptrbb = bb + off*2;
#endif

            vres0_0 = VFMVVF_FLOAT(0, vlmax);
            vres0_1 = VFMVVF_FLOAT(0, vlmax);
            
            vres1_0 = VFMVVF_FLOAT(0, vlmax);
            vres1_1 = VFMVVF_FLOAT(0, vlmax);


#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
            temp = bk-off;
#elif defined(LEFT)
            temp = off+2;   // number of values in A
#else
            temp = off+2;   // number of values in B
#endif
            for (k = temp; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                VLSEG2_FLOAT(&va0, &va1, ptrba, vl);
                VLSEG2_FLOAT(&vb0, &vb1, ptrbb, vl);
            
                vres0_0 = VFMACCVV_FLOAT(vres0_0, va0, vb0, vl);
                vres1_0 = VFMACCVV_FLOAT(vres1_0, va0, vb1, vl);
            
                vres0_1 = VFMACCVV_FLOAT(vres0_1, va1, vb0, vl);
                vres1_1 = VFMACCVV_FLOAT(vres1_1, va1, vb1, vl);
            
                ptrba += vl * 2;
                ptrbb += vl * 2;
            }

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres0_0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres0_1, v_z0, vlmax);
            vsum2 = VFREDSUMVS_FLOAT(vsum2, vres1_0, v_z0, vlmax);
            vsum3 = VFREDSUMVS_FLOAT(vsum3, vres1_1, v_z0, vlmax);

            C0[0] = alpha * VFMVFS_FLOAT_M1(vsum0);
            C0[1] = alpha * VFMVFS_FLOAT_M1(vsum1);
            C1[0] = alpha * VFMVFS_FLOAT_M1(vsum2);
            C1[1] = alpha * VFMVFS_FLOAT_M1(vsum3);


#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
            temp = bk - off;
#ifdef LEFT
            temp -= 2; // number of values in A
#else
            temp -= 2; // number of values in B
#endif
            ptrba += temp*2;
            ptrbb += temp*2;
#endif

#ifdef LEFT
            off += 2; // number of values in A
#endif

            C0 = C0+2;
            C1 = C1+2;

        }

        if ( bm & 1 ) // do any 1x2 loop
        {

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
            ptrbb = bb;
#else
            ptrba += off*1;
            ptrbb = bb + off*2;
#endif


            vres0_0 = VFMVVF_FLOAT(0, vlmax);
            vres1_0 = VFMVVF_FLOAT(0, vlmax);


#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
            temp = bk-off;
#elif defined(LEFT)
            temp = off+1;   // number of values in A
#else
            temp = off+2;   // number of values in B
#endif

            for (k = temp; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                va0 = VLEV_FLOAT(ptrba, vl);
                VLSEG2_FLOAT(&vb0, &vb1, ptrbb, vl);
            
                vres0_0 = VFMACCVV_FLOAT(vres0_0, va0, vb0, vl);
                vres1_0 = VFMACCVV_FLOAT(vres1_0, va0, vb1, vl);
                    
                ptrba += vl;
                ptrbb += vl * 2;
            }

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres0_0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres1_0, v_z0, vlmax);
            C0[0] = alpha * VFMVFS_FLOAT_M1(vsum0);
            C1[0] = alpha * VFMVFS_FLOAT_M1(vsum1);


#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
            temp = bk - off;
#ifdef LEFT
            temp -= 1; // number of values in A
#else
            temp -= 2; // number of values in B
#endif
            ptrba += temp*1;
            ptrbb += temp*2;
#endif

#ifdef LEFT
            off += 1; // number of values in A
#endif

            C0 = C0+1;
            C1 = C1+1;

        }


#if defined(TRMMKERNEL) && !defined(LEFT)
        off += 2;
#endif

        k = (bk<<1);
        bb = bb+k;
        i = (ldc<<1);
        C = C+i;
    }

    for (j=0; j<(bn&1); j+=1) // do the Mx1 loops
    {
        C0 = C;

#if defined(TRMMKERNEL) &&  defined(LEFT)
    off = offset;
#endif

        ptrba = ba;

        for (i=0; i<bm/4; i+=1) // do blocks of 4x1 loops
        {

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
            ptrbb = bb;
#else
            ptrba += off*4;
            ptrbb = bb + off*1;
#endif

            vres0_0 = VFMVVF_FLOAT(0, vlmax);
            vres0_1 = VFMVVF_FLOAT(0, vlmax);
            vres0_2 = VFMVVF_FLOAT(0, vlmax);
            vres0_3 = VFMVVF_FLOAT(0, vlmax);


#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
            temp = bk-off;
#elif defined(LEFT)
            temp = off+4;   // number of values in A
#else
            temp = off+1;   // number of values in B
#endif

            for (k = temp; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                VLSEG4_FLOAT(&va0, &va1, &va2, &va3, ptrba, vl);
                vb0 = VLEV_FLOAT(ptrbb, vl);
            
                vres0_0 = VFMACCVV_FLOAT(vres0_0, va0, vb0, vl);
           
                vres0_1 = VFMACCVV_FLOAT(vres0_1, va1, vb0, vl);
            
                vres0_2 = VFMACCVV_FLOAT(vres0_2, va2, vb0, vl);
            
                vres0_3 = VFMACCVV_FLOAT(vres0_3, va3, vb0, vl);
            
                ptrba += vl * 4;
                ptrbb += vl;
            }

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres0_0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres0_1, v_z0, vlmax);
            vsum2 = VFREDSUMVS_FLOAT(vsum2, vres0_2, v_z0, vlmax);
            vsum3 = VFREDSUMVS_FLOAT(vsum3, vres0_3, v_z0, vlmax);
            C0[0] = alpha * VFMVFS_FLOAT_M1(vsum0);
            C0[1] = alpha * VFMVFS_FLOAT_M1(vsum1);
            C0[2] = alpha * VFMVFS_FLOAT_M1(vsum2);
            C0[3] = alpha * VFMVFS_FLOAT_M1(vsum3);


#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
            temp = bk - off;
#ifdef LEFT
            temp -= 4; // number of values in A
#else
            temp -= 1; // number of values in B
#endif
            ptrba += temp*4;
            ptrbb += temp*1;
#endif

#ifdef LEFT
            off += 4; // number of values in A
#endif

            C0 = C0+4;

        }

        if ( bm & 2 ) // do any 2x1 loop
        {

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
            ptrbb = bb;
#else
            ptrba += off*2;
            ptrbb = bb + off*1;
#endif

            vres0_0 = VFMVVF_FLOAT(0, vlmax);
            vres0_1 = VFMVVF_FLOAT(0, vlmax);

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
            temp = bk-off;
#elif defined(LEFT)
            temp = off+2;   // number of values in A
#else
            temp = off+1;   // number of values in B
#endif

            for (k = temp; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                VLSEG2_FLOAT(&va0, &va1, ptrba, vl);
                vb0 = VLEV_FLOAT(ptrbb, vl);
            
                vres0_0 = VFMACCVV_FLOAT(vres0_0, va0, vb0, vl);
           
                vres0_1 = VFMACCVV_FLOAT(vres0_1, va1, vb0, vl);
           
                ptrba += vl * 2;
                ptrbb += vl;
            }

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres0_0, v_z0, vlmax);
            vsum1 = VFREDSUMVS_FLOAT(vsum1, vres0_1, v_z0, vlmax);
            C0[0] = alpha * VFMVFS_FLOAT_M1(vsum0);
            C0[1] = alpha * VFMVFS_FLOAT_M1(vsum1);


#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
            temp = bk - off;
#ifdef LEFT
            temp -= 2; // number of values in A
#else
            temp -= 1; // number of values in B
#endif
            ptrba += temp*2;
            ptrbb += temp*1;
#endif

#ifdef LEFT
            off += 2; // number of values in A
#endif

            C0 = C0+2;

        }

        if ( bm & 1 ) // do any 1x1 loop
        {

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
            ptrbb = bb;
#else
            ptrba += off*1;
            ptrbb = bb + off*1;
#endif

            vres0_0 = VFMVVF_FLOAT(0, vlmax);

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
            temp = bk-off;
#elif defined(LEFT)
            temp = off+1;   // number of values in A
#else
            temp = off+1;   // number of values in B
#endif

            for (k = temp; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                va0 = VLEV_FLOAT(ptrba, vl);
                vb0 = VLEV_FLOAT(ptrbb, vl);
            
                vres0_0 = VFMACCVV_FLOAT(vres0_0, va0, vb0, vl);
                  
                ptrba += vl;
                ptrbb += vl;
            }

            vsum0 = VFREDSUMVS_FLOAT(vsum0, vres0_0, v_z0, vlmax);
            C0[0] = alpha * VFMVFS_FLOAT_M1(vsum0);


#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
            temp = bk - off;
#ifdef LEFT
            temp -= 1; // number of values in A
#else
            temp -= 1; // number of values in B
#endif
            ptrba += temp*1;
            ptrbb += temp*1;
#endif

#ifdef LEFT
            off += 1; // number of values in A
#endif

            C0 = C0+1;

        }

#if defined(TRMMKERNEL) && !defined(LEFT)
        off += 1;
#endif

        k = (bk<<0);
        bb = bb+k;
        C = C+ldc;
    }
    return 0;
}
