/***************************************************************************
Copyright (c) 2020, The OpenBLAS Project
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
#define VSETVL_MAX vsetvlmax_e32m1()
#define FLOAT_V_T vfloat32m4_t
#define FLOAT_V_T_M1 vfloat32m1_t
#define VFMVFS_FLOAT vfmv_f_s_f32m1_f32
#define VLEV_FLOAT vle32_v_f32m4
#define VLSEV_FLOAT vlse32_v_f32m4
#define VFREDSUM_FLOAT vfredusum_vs_f32m4_f32m1
#define VFMACCVV_FLOAT vfmacc_vv_f32m4
#define VFMVVF_FLOAT vfmv_v_f_f32m4
#define VFMVVF_FLOAT_M1 vfmv_v_f_f32m1
#define VFDOTVV_FLOAT vfdot_vv_f32m4
#define VFMULVV_FLOAT vfmul_vv_f32m4
#define VFMSACVV_FLOAT vfmsac_vv_f32m4
#define VFNMSACVV_FLOAT vfnmsac_vv_f32m4
#else
#define VSETVL(n) vsetvl_e64m4(n)
#define VSETVL_MAX vsetvlmax_e64m1()
#define FLOAT_V_T vfloat64m4_t
#define FLOAT_V_T_M1 vfloat64m1_t
#define VFMVFS_FLOAT vfmv_f_s_f64m1_f64
#define VLEV_FLOAT vle64_v_f64m4
#define VLSEV_FLOAT vlse64_v_f64m4
#define VFREDSUM_FLOAT vfredusum_vs_f64m4_f64m1
#define VFMACCVV_FLOAT vfmacc_vv_f64m4
#define VFMVVF_FLOAT vfmv_v_f_f64m4
#define VFMVVF_FLOAT_M1 vfmv_v_f_f64m1
#define VFDOTVV_FLOAT vfdot_vv_f64m4
#define VFMULVV_FLOAT vfmul_vv_f64m4
#define VFMSACVV_FLOAT vfmsac_vv_f64m4
#define VFNMSACVV_FLOAT vfnmsac_vv_f64m4
#endif

OPENBLAS_COMPLEX_FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y)
{
        BLASLONG i=0, j=0;
        BLASLONG ix=0,iy=0;
        FLOAT dot[2];
        OPENBLAS_COMPLEX_FLOAT result;

        dot[0]=0.0;
        dot[1]=0.0;

        CREAL(result) = 0.0;
        CIMAG(result) = 0.0;

        if ( n < 1 )  return(result);

        unsigned int gvl = 0;
        FLOAT_V_T_M1 v_res, v_z0;
        gvl = VSETVL_MAX;
        v_res = VFMVVF_FLOAT_M1(0, gvl);
        v_z0 = VFMVVF_FLOAT_M1(0, gvl);

        FLOAT_V_T vr0, vr1, vx0, vx1, vy0, vy1;
        gvl = VSETVL(n);
        vr0 = VFMVVF_FLOAT(0, gvl);
        vr1 = VFMVVF_FLOAT(0, gvl);
        BLASLONG stride_x = inc_x * 2 * sizeof(FLOAT);
        BLASLONG stride_y = inc_y * 2 * sizeof(FLOAT);
        BLASLONG inc_xv = inc_x * 2 * gvl;
        BLASLONG inc_yv = inc_y * 2 * gvl;

        for(i=0,j=0; i<n/gvl; i++){
                vx0 = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                vx1 = VLSEV_FLOAT(&x[ix+1], stride_x, gvl);
                vy0 = VLSEV_FLOAT(&y[iy], stride_y, gvl);
                vy1 = VLSEV_FLOAT(&y[iy+1], stride_y, gvl);

                vr0 = VFMACCVV_FLOAT(vr0, vx0, vy0, gvl);
                vr1 = VFMACCVV_FLOAT(vr1, vx0, vy1, gvl);
#if !defined(CONJ)
                vr0 = VFNMSACVV_FLOAT(vr0, vx1, vy1, gvl);
                vr1 = VFMACCVV_FLOAT(vr1, vx1, vy0, gvl);
#else
                vr0 = VFMACCVV_FLOAT(vr0, vx1, vy1, gvl);
                vr1 = VFNMSACVV_FLOAT(vr1, vx1, vy0, gvl);
#endif
                j += gvl;
                ix += inc_xv;
                iy += inc_yv;
        }
        v_res = VFREDSUM_FLOAT(v_res, vr0, v_z0, gvl);
        dot[0] += VFMVFS_FLOAT(v_res);
        v_res = VFREDSUM_FLOAT(v_res, vr1, v_z0, gvl);
        dot[1] += VFMVFS_FLOAT(v_res);
        //tail
        if(j < n){
                gvl = VSETVL(n-j);
                vx0 = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                vx1 = VLSEV_FLOAT(&x[ix+1], stride_x, gvl);
                vy0 = VLSEV_FLOAT(&y[iy], stride_y, gvl);
                vy1 = VLSEV_FLOAT(&y[iy+1], stride_y, gvl);

#if !defined(CONJ)
                vr0 = VFMULVV_FLOAT(vx1, vy1, gvl);
                vr0 = VFMSACVV_FLOAT(vr0, vx0, vy0, gvl);
                vr1 = VFMULVV_FLOAT(vx0, vy1, gvl);
                vr1 = VFMACCVV_FLOAT(vr1, vx1, vy0, gvl);
#else
                vr0 = VFMULVV_FLOAT(vx0, vy0, gvl);
                vr0 = VFMACCVV_FLOAT(vr0, vx1, vy1, gvl);
                vr1 = VFMULVV_FLOAT(vx1, vy0, gvl);
                vr1 = VFMSACVV_FLOAT(vr1, vx0, vy1, gvl);
#endif
                v_res = VFREDSUM_FLOAT(v_res, vr0, v_z0, gvl);
                dot[0] += VFMVFS_FLOAT(v_res);
                v_res = VFREDSUM_FLOAT(v_res, vr1, v_z0, gvl);
                dot[1] += VFMVFS_FLOAT(v_res);
        }
        CREAL(result) = dot[0];
        CIMAG(result) = dot[1];
        return(result);
}
