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
#define RVV_EFLOAT RVV_E32
#define RVV_M RVV_M4
#define FLOAT_V_T float32xm4_t
#define VLEV_FLOAT vlev_float32xm4
#define VLSEV_FLOAT vlsev_float32xm4
#define VSEV_FLOAT vsev_float32xm4
#define VSSEV_FLOAT vssev_float32xm4
#define VFREDSUM_FLOAT vfredsumvs_float32xm4
#define VFMACCVV_FLOAT vfmaccvv_float32xm4
#define VFMACCVF_FLOAT vfmaccvf_float32xm4
#define VFMVVF_FLOAT vfmvvf_float32xm4
#define VFDOTVV_FLOAT vfdotvv_float32xm4
#define VFMULVV_FLOAT vfmulvv_float32xm4
#else
#define RVV_EFLOAT RVV_E64
#define RVV_M RVV_M4
#define FLOAT_V_T float64xm4_t
#define VLEV_FLOAT vlev_float64xm4
#define VLSEV_FLOAT vlsev_float64xm4
#define VSEV_FLOAT vsev_float64xm4
#define VSSEV_FLOAT vssev_float64xm4
#define VFREDSUM_FLOAT vfredsumvs_float64xm4
#define VFMACCVV_FLOAT vfmaccvv_float64xm4
#define VFMACCVF_FLOAT vfmaccvf_float64xm4
#define VFMVVF_FLOAT vfmvvf_float64xm4
#define VFDOTVV_FLOAT vfdotvv_float64xm4
#define VFMULVV_FLOAT vfmulvv_float64xm4
#endif

int CNAME(BLASLONG m, BLASLONG offset, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *buffer)
{
        BLASLONG i, j, k;
        BLASLONG ix,iy;
        BLASLONG jx,jy;
        FLOAT temp1;
        FLOAT temp2;
        FLOAT *a_ptr = a;
        unsigned int gvl = 0;

        FLOAT_V_T va, vx, vy, vr;
        BLASLONG stride_x, stride_y, inc_xv, inc_yv;

        BLASLONG m1 = m - offset;
        if(inc_x == 1 && inc_y == 1){
                a_ptr += m1 * lda;
                for (j=m1; j<m; j++)
                {
                        temp1 = alpha * x[j];
                        temp2 = 0.0;
                        if(j > 0){
                                i = 0;
                                gvl = vsetvli(j, RVV_EFLOAT, RVV_M);
                                vr = VFMVVF_FLOAT(0, gvl);
                                for(k = 0; k < j / gvl; k++){
                                        vy = VLEV_FLOAT(&y[i], gvl);
                                        va = VLEV_FLOAT(&a_ptr[i], gvl);
                                        vy = VFMACCVF_FLOAT(vy, temp1, va, gvl);
                                        VSEV_FLOAT(&y[i], vy, gvl);

                                        vx = VLEV_FLOAT(&x[i], gvl);
                                        vr = VFMACCVV_FLOAT(vr, vx, va, gvl);

                                        i += gvl;
                                }
                                va = VFMVVF_FLOAT(0, gvl);
                                va = VFREDSUM_FLOAT(vr, va, gvl);
                                temp2 = va[0];
                                if(i < j){
                                        gvl = vsetvli(j-i, RVV_EFLOAT, RVV_M);
                                        vy = VLEV_FLOAT(&y[i], gvl);
                                        va = VLEV_FLOAT(&a_ptr[i], gvl);
                                        vy = VFMACCVF_FLOAT(vy, temp1, va, gvl);
                                        VSEV_FLOAT(&y[i], vy, gvl);

                                        vx = VLEV_FLOAT(&x[i], gvl);
                                        vr = VFMULVV_FLOAT(vx, va, gvl);
                                        va = VFMVVF_FLOAT(0, gvl);
                                        va = VFREDSUM_FLOAT(vr, va, gvl);
                                        temp2 += va[0];
                                }
                        }
                        y[j] += temp1 * a_ptr[j] + alpha * temp2;
                        a_ptr += lda;
                }
        }else if(inc_x == 1){
                jy = m1 * inc_y;
                a_ptr += m1 * lda;
                stride_y = inc_y * sizeof(FLOAT);
                for (j=m1; j<m; j++)
                {
                        temp1 = alpha * x[j];
                        temp2 = 0.0;
                        if(j > 0){
                                iy = 0;
                                i = 0;
                                gvl = vsetvli(j, RVV_EFLOAT, RVV_M);
                                inc_yv = inc_y * gvl;
                                vr = VFMVVF_FLOAT(0, gvl);
                                for(k = 0; k < j / gvl; k++){
                                        vy = VLSEV_FLOAT(&y[iy], stride_y, gvl);
                                        va = VLEV_FLOAT(&a_ptr[i], gvl);
                                        vy = VFMACCVF_FLOAT(vy, temp1, va, gvl);
                                        VSSEV_FLOAT(&y[iy], stride_y, vy, gvl);

                                        vx = VLEV_FLOAT(&x[i], gvl);
                                        vr = VFMACCVV_FLOAT(vr, vx, va, gvl);

                                        i += gvl;
                                        iy += inc_yv;
                                }
                                va = VFMVVF_FLOAT(0, gvl);
                                va = VFREDSUM_FLOAT(vr, va, gvl);
                                temp2 = va[0];
                                if(i < j){
                                        gvl = vsetvli(j-i, RVV_EFLOAT, RVV_M);
                                        vy = VLSEV_FLOAT(&y[iy], stride_y, gvl);
                                        va = VLEV_FLOAT(&a_ptr[i], gvl);
                                        vy = VFMACCVF_FLOAT(vy, temp1, va, gvl);
                                        VSSEV_FLOAT(&y[iy], stride_y, vy, gvl);

                                        vx = VLEV_FLOAT(&x[i], gvl);
                                        vr = VFMULVV_FLOAT(vx, va, gvl);
                                        va = VFMVVF_FLOAT(0, gvl);
                                        va = VFREDSUM_FLOAT(vr, va, gvl);
                                        temp2 += va[0];
                                }
                        }
                        y[jy] += temp1 * a_ptr[j] + alpha * temp2;
                        a_ptr += lda;
                        jy    += inc_y;
                }
        }else if(inc_y == 1){
                jx = m1 * inc_x;
                a_ptr += m1 * lda;
                stride_x = inc_x * sizeof(FLOAT);
                for (j=m1; j<m; j++)
                {
                        temp1 = alpha * x[jx];
                        temp2 = 0.0;
                        if(j > 0){
                                ix = 0;
                                i = 0;
                                gvl = vsetvli(j, RVV_EFLOAT, RVV_M);
                                inc_xv = inc_x * gvl;
                                vr = VFMVVF_FLOAT(0, gvl);
                                for(k = 0; k < j / gvl; k++){
                                        vy = VLEV_FLOAT(&y[i], gvl);
                                        va = VLEV_FLOAT(&a_ptr[i], gvl);
                                        vy = VFMACCVF_FLOAT(vy, temp1, va, gvl);
                                        VSEV_FLOAT(&y[i], vy, gvl);

                                        vx = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                                        vr = VFMACCVV_FLOAT(vr, vx, va, gvl);

                                        i += gvl;
                                        ix += inc_xv;
                                }
                                va = VFMVVF_FLOAT(0, gvl);
                                va = VFREDSUM_FLOAT(vr, va, gvl);
                                temp2 = va[0];
                                if(i < j){
                                        gvl = vsetvli(j-i, RVV_EFLOAT, RVV_M);
                                        vy = VLEV_FLOAT(&y[i], gvl);
                                        va = VLEV_FLOAT(&a_ptr[i], gvl);
                                        vy = VFMACCVF_FLOAT(vy, temp1, va, gvl);
                                        VSEV_FLOAT(&y[i], vy, gvl);

                                        vx = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                                        vr = VFMULVV_FLOAT(vx, va, gvl);
                                        va = VFMVVF_FLOAT(0, gvl);
                                        va = VFREDSUM_FLOAT(vr, va, gvl);
                                        temp2 += va[0];
                                }
                        }
                        y[j] += temp1 * a_ptr[j] + alpha * temp2;
                        a_ptr += lda;
                        jx    += inc_x;
                }
        }else{
                jx = m1 * inc_x;
                jy = m1 * inc_y;
                a_ptr += m1 * lda;
                stride_x = inc_x * sizeof(FLOAT);
                stride_y = inc_y * sizeof(FLOAT);
                for (j=m1; j<m; j++)
                {
                        temp1 = alpha * x[jx];
                        temp2 = 0.0;
                        if(j > 0){
                                ix = 0;
                                iy = 0;
                                i = 0;
                                gvl = vsetvli(j, RVV_EFLOAT, RVV_M);
                                inc_xv = inc_x * gvl;
                                inc_yv = inc_y * gvl;
                                vr = VFMVVF_FLOAT(0, gvl);
                                for(k = 0; k < j / gvl; k++){
                                        vy = VLSEV_FLOAT(&y[iy], stride_y, gvl);
                                        va = VLEV_FLOAT(&a_ptr[i], gvl);
                                        vy = VFMACCVF_FLOAT(vy, temp1, va, gvl);
                                        VSSEV_FLOAT(&y[iy], stride_y, vy, gvl);

                                        vx = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                                        vr = VFMACCVV_FLOAT(vr, vx, va, gvl);

                                        i += gvl;
                                        ix += inc_xv;
                                        iy += inc_yv;
                                }
                                va = VFMVVF_FLOAT(0, gvl);
                                va = VFREDSUM_FLOAT(vr, va, gvl);
                                temp2 = va[0];
                                if(i < j){
                                        gvl = vsetvli(j-i, RVV_EFLOAT, RVV_M);
                                        vy = VLSEV_FLOAT(&y[iy], stride_y, gvl);
                                        va = VLEV_FLOAT(&a_ptr[i], gvl);
                                        vy = VFMACCVF_FLOAT(vy, temp1, va, gvl);
                                        VSSEV_FLOAT(&y[iy], stride_y, vy, gvl);

                                        vx = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                                        vr = VFMULVV_FLOAT(vx, va, gvl);
                                        va = VFMVVF_FLOAT(0, gvl);
                                        va = VFREDSUM_FLOAT(vr, va, gvl);
                                        temp2 += va[0];
                                }
                        }
                        y[jy] += temp1 * a_ptr[j] + alpha * temp2;
                        a_ptr += lda;
                        jx    += inc_x;
                        jy    += inc_y;
                }
        }
        return(0);
}

