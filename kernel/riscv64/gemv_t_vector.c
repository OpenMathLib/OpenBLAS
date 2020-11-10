/***************************************************************************
Copyright (c) 2013, The OpenBLAS Project
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
#endif

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *buffer)
{
	BLASLONG i = 0, j = 0, k = 0;
	BLASLONG ix = 0, iy = 0;
	FLOAT *a_ptr = a;
        FLOAT temp;

        FLOAT_V_T va, vr, vx;
        unsigned int gvl = 0;
        if(inc_x == 1){
                for(i = 0; i < n; i++){
                        gvl = vsetvli(m, RVV_EFLOAT, RVV_M);
                        j = 0;
                        vr = VFMVVF_FLOAT(0, gvl);
                        for(k = 0; k < m/gvl; k++){
                                va = VLEV_FLOAT(&a_ptr[j], gvl);
                                vx = VLEV_FLOAT(&x[j], gvl);
                                vr = VFMACCVV_FLOAT(vr, va, vx, gvl);
                                j += gvl;
                        }
                        va = VFMVVF_FLOAT(0, gvl);
                        va = VFREDSUM_FLOAT(vr, va, gvl);
                        temp = va[0];
                        if(j < m){
                                gvl = vsetvli(m-j, RVV_EFLOAT, RVV_M);
                                va = VLEV_FLOAT(&a_ptr[j], gvl);
                                vx = VLEV_FLOAT(&x[j], gvl);
                                vr = VFMULVV_FLOAT(va, vx, gvl);

                                va = VFMVVF_FLOAT(0, gvl);
                                va = VFREDSUM_FLOAT(vr, va, gvl);
                                temp += va[0];
                        }
                        y[iy] += alpha * temp;
                        iy += inc_y;
                        a_ptr += lda;
                }
        }else{
                gvl = vsetvli(m, RVV_EFLOAT, RVV_M);
                BLASLONG stride_x = inc_x * sizeof(FLOAT);
                BLASLONG inc_xv = inc_x * gvl;
                for(i = 0; i < n; i++){
                        gvl = vsetvli(m, RVV_EFLOAT, RVV_M);
                        j = 0;
                        ix = 0;
                        vr = VFMVVF_FLOAT(0, gvl);
                        for(k = 0; k < m/gvl; k++){
                                va = VLEV_FLOAT(&a_ptr[j], gvl);
                                vx = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                                vr = VFMACCVV_FLOAT(vr, va, vx, gvl);
                                j += gvl;
                                ix += inc_xv;
                        }
                        va = VFMVVF_FLOAT(0, gvl);
                        va = VFREDSUM_FLOAT(vr, va, gvl);
                        temp = va[0];
                        if(j < m){
                                gvl = vsetvli(m-j, RVV_EFLOAT, RVV_M);
                                va = VLEV_FLOAT(&a_ptr[j], gvl);
                                vx = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                                vr = VFMULVV_FLOAT(va, vx, gvl);

                                va = VFMVVF_FLOAT(0, gvl);
                                va = VFREDSUM_FLOAT(vr, va, gvl);
                                temp += va[0];
                        }
                        y[iy] += alpha * temp;
                        iy += inc_y;
                        a_ptr += lda;
                }
        }
	return(0);
}

