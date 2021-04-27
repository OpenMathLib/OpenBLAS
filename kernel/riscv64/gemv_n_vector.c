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
#define VFMACCVF_FLOAT vfmaccvf_float32xm4
#else
#define RVV_EFLOAT RVV_E64
#define RVV_M RVV_M4
#define FLOAT_V_T float64xm4_t
#define VLEV_FLOAT vlev_float64xm4
#define VLSEV_FLOAT vlsev_float64xm4
#define VSEV_FLOAT vsev_float64xm4
#define VSSEV_FLOAT vssev_float64xm4
#define VFMACCVF_FLOAT vfmaccvf_float64xm4
#endif

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *buffer)
{
	BLASLONG i = 0, j = 0, k = 0;
        BLASLONG ix = 0, iy = 0;

	if(n < 0)  return(0);
        FLOAT *a_ptr = a;
        FLOAT temp = 0.0;
        FLOAT_V_T va0, va1, vy0, vy1;
        unsigned int gvl = 0;
        if(inc_y == 1){
                gvl = vsetvli(m, RVV_EFLOAT, RVV_M);
                if(gvl <= m/2){
                        for(k=0,j=0; k<m/(2*gvl); k++){
                                a_ptr = a;
                                ix = 0;
                                vy0 = VLEV_FLOAT(&y[j], gvl);
                                vy1 = VLEV_FLOAT(&y[j+gvl], gvl);
                                for(i = 0; i < n; i++){
                                        temp = alpha * x[ix];
                                        va0 = VLEV_FLOAT(&a_ptr[j], gvl);
                                        vy0 = VFMACCVF_FLOAT(vy0, temp, va0, gvl);

                                        va1 = VLEV_FLOAT(&a_ptr[j+gvl], gvl);
                                        vy1 = VFMACCVF_FLOAT(vy1, temp, va1, gvl);
                                        a_ptr += lda;
                                        ix += inc_x;
                                }
                                VSEV_FLOAT(&y[j], vy0, gvl);
                                VSEV_FLOAT(&y[j+gvl], vy1, gvl);
                                j += gvl * 2;
                        }
                }
                //tail
                for(;j < m;){
                        gvl = vsetvli(m-j, RVV_EFLOAT, RVV_M);
                        a_ptr = a;
                        ix = 0;
                        vy0 = VLEV_FLOAT(&y[j], gvl);
                        for(i = 0; i < n; i++){
                                temp = alpha * x[ix];
                                va0 = VLEV_FLOAT(&a_ptr[j], gvl);
                                vy0 = VFMACCVF_FLOAT(vy0, temp, va0, gvl);

                                a_ptr += lda;
                                ix += inc_x;
                        }
                        VSEV_FLOAT(&y[j], vy0, gvl);
                        j += gvl;
                }
        }else{
                BLASLONG stride_y = inc_y * sizeof(FLOAT);
                gvl = vsetvli(m, RVV_EFLOAT, RVV_M);
                if(gvl <= m/2){
                        BLASLONG inc_yv = inc_y * gvl;
                        for(k=0,j=0; k<m/(2*gvl); k++){
                                a_ptr = a;
                                ix = 0;
                                vy0 = VLSEV_FLOAT(&y[iy], stride_y, gvl);
                                vy1 = VLSEV_FLOAT(&y[iy+inc_yv], stride_y, gvl);
                                for(i = 0; i < n; i++){
                                        temp = alpha * x[ix];
                                        va0 = VLEV_FLOAT(&a_ptr[j], gvl);
                                        vy0 = VFMACCVF_FLOAT(vy0, temp, va0, gvl);

                                        va1 = VLEV_FLOAT(&a_ptr[j+gvl], gvl);
                                        vy1 = VFMACCVF_FLOAT(vy1, temp, va1, gvl);
                                        a_ptr += lda;
                                        ix += inc_x;
                                }
                                VSSEV_FLOAT(&y[iy], stride_y, vy0, gvl);
                                VSSEV_FLOAT(&y[iy+inc_yv], stride_y, vy1, gvl);
                                j += gvl * 2;
                                iy += inc_yv * 2;
                        }
                }
                //tail
                for(;j < m;){
                        gvl = vsetvli(m-j, RVV_EFLOAT, RVV_M);
                        a_ptr = a;
                        ix = 0;
                        vy0 = VLSEV_FLOAT(&y[j*inc_y], stride_y, gvl);
                        for(i = 0; i < n; i++){
                                temp = alpha * x[ix];
                                va0 = VLEV_FLOAT(&a_ptr[j], gvl);
                                vy0 = VFMACCVF_FLOAT(vy0, temp, va0, gvl);

                                a_ptr += lda;
                                ix += inc_x;
                        }
                        VSSEV_FLOAT(&y[j*inc_y], stride_y, vy0, gvl);
                        j += gvl;
                }
        }
	return(0);
}


