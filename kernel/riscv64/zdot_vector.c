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
#define VFREDSUM_FLOAT vfredsumvs_float32xm4
#define VFMACCVV_FLOAT vfmaccvv_float32xm4
#define VFMVVF_FLOAT vfmvvf_float32xm4
#define VFDOTVV_FLOAT vfdotvv_float32xm4
#define VFMULVV_FLOAT vfmulvv_float32xm4
#define VFMSACVV_FLOAT vfmsacvv_float32xm4
#define VFNMSACVV_FLOAT vfnmsacvv_float32xm4
#else
#define RVV_EFLOAT RVV_E64
#define RVV_M RVV_M4
#define FLOAT_V_T float64xm4_t
#define VLEV_FLOAT vlev_float64xm4
#define VLSEV_FLOAT vlsev_float64xm4
#define VFREDSUM_FLOAT vfredsumvs_float64xm4
#define VFMACCVV_FLOAT vfmaccvv_float64xm4
#define VFMVVF_FLOAT vfmvvf_float64xm4
#define VFDOTVV_FLOAT vfdotvv_float64xm4
#define VFMULVV_FLOAT vfmulvv_float64xm4
#define VFMSACVV_FLOAT vfmsacvv_float64xm4
#define VFNMSACVV_FLOAT vfnmsacvv_float64xm4
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

        FLOAT_V_T vr0, vr1, vx0, vx1, vy0, vy1;
        gvl = vsetvli(n, RVV_EFLOAT, RVV_M);
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
        vx0 = VFMVVF_FLOAT(0, gvl);
        vr0 = VFREDSUM_FLOAT(vr0, vx0, gvl);
        dot[0] += vr0[0];
        vr1 = VFREDSUM_FLOAT(vr1, vx0, gvl);
        dot[1] += vr1[0];
        //tail
        if(j < n){
                gvl = vsetvli(n-j, RVV_EFLOAT, RVV_M);
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
                vx0 = VFMVVF_FLOAT(0, gvl);
                vr0 = VFREDSUM_FLOAT(vr0, vx0, gvl);
                dot[0] += vr0[0];
                vr1 = VFREDSUM_FLOAT(vr1, vx0, gvl);
                dot[1] += vr1[0];
        }
        CREAL(result) = dot[0];
        CIMAG(result) = dot[1];
        return(result);
}
