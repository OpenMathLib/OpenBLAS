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
#define VSETVL(n) vsetvl_e32m4(n)
#define VSETVL_MAX vsetvlmax_e32m4()
#define VSETVL_MAX_M1 vsetvlmax_e32m1()
#define FLOAT_V_T vfloat32m4_t
#define FLOAT_V_T_M1 vfloat32m1_t
#define VLEV_FLOAT vle32_v_f32m4
#define VLSEG_FLOAT vlseg2e32_v_f32m4
#define VFMVVF_FLOAT vfmv_v_f_f32m4
#define VFMACCVF_FLOAT vfmacc_vf_f32m4
#define VFMACCVV_FLOAT vfmacc_vv_f32m4
#define VFREDSUMVS_FLOAT vfredusum_vs_f32m4_f32m1
#define VFMVVF_FLOAT_M1 vfmv_v_f_f32m1
#define VFMVFS_FLOAT_M1 vfmv_f_s_f32m1_f32
#else
#define VSETVL(n) vsetvl_e64m4(n)
#define VSETVL_MAX vsetvlmax_e64m4()
#define VSETVL_MAX_M1 vsetvlmax_e64m1()
#define FLOAT_V_T vfloat64m4_t
#define FLOAT_V_T_M1 vfloat64m1_t
#define VLEV_FLOAT vle64_v_f64m4
#define VLSEG_FLOAT vlseg2e64_v_f64m4
#define VFMVVF_FLOAT vfmv_v_f_f64m4
#define VFMACCVF_FLOAT vfmacc_vf_f64m4
#define VFMACCVV_FLOAT vfmacc_vv_f64m4
#define VFREDSUMVS_FLOAT vfredusum_vs_f64m4_f64m1
#define VFMVVF_FLOAT_M1 vfmv_v_f_f64m1
#define VFMVFS_FLOAT_M1 vfmv_f_s_f64m1_f64
#endif


// Optimizes the implementation in ../generic/trmmkernel_2x2.c


int CNAME(BLASLONG bm,BLASLONG bn,BLASLONG bk,FLOAT alpha,FLOAT* ba,FLOAT* bb,FLOAT* C,BLASLONG ldc
#ifdef TRMMKERNEL
        ,BLASLONG offset
#endif
        )
{
   BLASLONG i,j,k;
   FLOAT *C0,*C1,*ptrba,*ptrbb;
   BLASLONG off, temp;

   FLOAT_V_T va0, va1, vb0, vb1;
   FLOAT_V_T vres0, vres1, vres2, vres3;
   FLOAT_V_T_M1 v_res, v_z0;
   v_z0 = VFMVVF_FLOAT_M1(0, VSETVL_MAX_M1);
   size_t vl;
   size_t vlmax = VSETVL_MAX;
   
#if defined(TRMMKERNEL) && !defined(LEFT)
   off = -offset;
#else
   off = 0;
#endif

   for (j = bn/2; j > 0; j--)
   {
        C0 = C;
        C1 = C0+ldc;
#if defined(TRMMKERNEL) && defined(LEFT)
        off = offset;
#endif
        ptrba = ba;

        for (i = bm/2; i > 0; i--)
        {
#if (defined(LEFT) &&  defined(TRANSA)) || \
              (!defined(LEFT) && !defined(TRANSA))
             ptrbb = bb;
#else
              ptrba += off*2;
              ptrbb = bb + off*2;
#endif

#if (defined(LEFT) && !defined(TRANSA)) || \
             (!defined(LEFT) && defined(TRANSA))
             temp = bk-off;
#elif defined(LEFT)
             temp = off+2;
#else
             temp = off+2;
#endif
            vres0 = VFMVVF_FLOAT(0.0, vlmax);
            vres1 = VFMVVF_FLOAT(0.0, vlmax);
            vres2 = VFMVVF_FLOAT(0.0, vlmax);
            vres3 = VFMVVF_FLOAT(0.0, vlmax);
            for (k = temp; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                VLSEG_FLOAT(&va0, &va1, ptrba, vl);
                VLSEG_FLOAT(&vb0, &vb1, ptrbb, vl);

                vres0 = VFMACCVV_FLOAT(vres0, va0, vb0, vl);
                vres1 = VFMACCVV_FLOAT(vres1, va1, vb0, vl);
                vres2 = VFMACCVV_FLOAT(vres2, va0, vb1, vl);
                vres3 = VFMACCVV_FLOAT(vres3, va1, vb1, vl);

                ptrba += vl * 2;
                ptrbb += vl * 2;
            }
            v_res = VFREDSUMVS_FLOAT(v_res, vres0, v_z0, vlmax);
            C0[0] = alpha * VFMVFS_FLOAT_M1(v_res);
            v_res = VFREDSUMVS_FLOAT(v_res, vres1, v_z0, vlmax);
            C0[1] = alpha * VFMVFS_FLOAT_M1(v_res);
            v_res = VFREDSUMVS_FLOAT(v_res, vres2, v_z0, vlmax);
            C1[0] = alpha * VFMVFS_FLOAT_M1(v_res);
            v_res = VFREDSUMVS_FLOAT(v_res, vres3, v_z0, vlmax);
            C1[1] = alpha * VFMVFS_FLOAT_M1(v_res);

#if ( defined(LEFT) && defined(TRANSA)) || \
             (!defined(LEFT) && !defined(TRANSA))
             temp = bk - off;
#ifdef LEFT
             temp -= 2;
#else
             temp -= 2;
#endif
             ptrba += temp*2;
             ptrbb += temp*2;
#endif
#ifdef LEFT
             off += 2;
#endif
             C0 = C0+2;
             C1 = C1+2;
        }

        if (bm & 1)
        {
#if (defined(LEFT) &&  defined(TRANSA)) ||(!defined(LEFT) && !defined(TRANSA))
             ptrbb = bb;
#else
             ptrba += off;
             ptrbb = bb+off*2;
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
             temp = bk-off;
#elif defined(LEFT)
             temp = off+1;
#else
             temp = off+2;
#endif
            vres0 = VFMVVF_FLOAT(0.0, vlmax);
            vres1 = VFMVVF_FLOAT(0.0, vlmax);

            for (k = temp; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                va0 = VLEV_FLOAT(ptrba, vl);
                VLSEG_FLOAT(&vb0, &vb1, ptrbb, vl);
                  
                vres0 = VFMACCVV_FLOAT(vres0, va0, vb0, vl);
                vres1 = VFMACCVV_FLOAT(vres1, va0, vb1, vl);

                ptrba += vl;
                ptrbb += vl * 2;

            }
            v_res = VFREDSUMVS_FLOAT(v_res, vres0, v_z0, vlmax);
            C0[0] = alpha * VFMVFS_FLOAT_M1(v_res);
            v_res = VFREDSUMVS_FLOAT(v_res, vres1, v_z0, vlmax);
            C1[0] = alpha * VFMVFS_FLOAT_M1(v_res);

#if ( defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
             temp = bk-off;
#ifdef LEFT
             temp -= 1;
#else
             temp -= 2;
#endif
             ptrba += temp;
             ptrbb += temp*2;
#endif
#ifdef LEFT
             off += 1;
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

   if (bn & 1)
   {
        C0 = C;
#if defined(TRMMKERNEL) &&  defined(LEFT)
        off = offset;
#endif
        ptrba = ba;

        for (i = bm/2; i > 0; i--)
        {
#if (defined(LEFT) &&  defined(TRANSA)) || \
              (!defined(LEFT) && !defined(TRANSA))
              ptrbb = bb;
#else
              ptrba += off*2;
              ptrbb = bb + off;
#endif


#if (defined(LEFT) && !defined(TRANSA)) || \
             (!defined(LEFT) && defined(TRANSA))
             temp = bk-off;
#elif defined(LEFT)
             temp = off+2;
#else
             temp = off+1;
#endif
             vres0 = VFMVVF_FLOAT(0.0, vlmax);
             vres1 = VFMVVF_FLOAT(0.0, vlmax);
             
             for (k = temp; k > 0; k -= vl)
             {
                 vl = VSETVL(k);
                 vb0 = VLEV_FLOAT(ptrbb, vl);
                 VLSEG_FLOAT(&va0, &va1, ptrba, vl);
                   
                 vres0 = VFMACCVV_FLOAT(vres0, vb0, va0, vl);
                 vres1 = VFMACCVV_FLOAT(vres1, vb0, va1, vl);
             
                 ptrba += vl * 2;
                 ptrbb += vl;
             
             }
             v_res = VFREDSUMVS_FLOAT(v_res, vres0, v_z0, vlmax);
             C0[0] = alpha * VFMVFS_FLOAT_M1(v_res);
             v_res = VFREDSUMVS_FLOAT(v_res, vres1, v_z0, vlmax);
             C0[1] = alpha * VFMVFS_FLOAT_M1(v_res);

#if ( defined(LEFT) &&  defined(TRANSA)) || \
             (!defined(LEFT) && !defined(TRANSA))
             temp = bk - off;
#ifdef LEFT
             temp -= 2;
#else
             temp -= 1;
#endif
             ptrba += temp*2;
             ptrbb += temp;
#endif
#ifdef LEFT
             off += 2;
#endif

             C0 = C0+2;
        }

        if (bm & 1)
        {
#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
             ptrbb = bb;
#else
             ptrba += off;
             ptrbb = bb+off;
#endif

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
             temp = bk-off;
#elif defined(LEFT)
             temp = off + 1;
#else
             temp = off + 1;
#endif
             vres0 = VFMVVF_FLOAT(0.0, vlmax);
             
             for (k = temp; k > 0; k -= vl)
             {
                 vl = VSETVL(k);
                 va0 = VLEV_FLOAT(ptrba, vl);
                 vb0 = VLEV_FLOAT(ptrbb, vl);
                   
                 vres0 = VFMACCVV_FLOAT(vres0, vb0, va0, vl);
                 ptrba += vl;
                 ptrbb += vl;
             }
             v_res = VFREDSUMVS_FLOAT(v_res, vres0, v_z0, vlmax);
             C0[0] = alpha * VFMVFS_FLOAT_M1(v_res);

#if ( defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
             temp = bk-off;
#ifdef LEFT
             temp -= 1;
#else
             temp -= 1;
#endif
             ptrba += temp;
             ptrbb += temp;
#endif
#ifdef LEFT
             off += 1;
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

