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
#define VLSEV_FLOAT vlsev_float32xm4
#define VFREDSUM_FLOAT vfredsumvs_float32xm4
#define VFMACCVV_FLOAT vfmaccvv_float32xm4
#define VFNMSACVV_FLOAT vfnmsacvv_float32xm4
#define VFMVVF_FLOAT vfmvvf_float32xm4
#define VFMULVV_FLOAT vfmulvv_float32xm4
#else
#define RVV_EFLOAT RVV_E64
#define RVV_M RVV_M4
#define FLOAT_V_T float64xm4_t
#define VLSEV_FLOAT vlsev_float64xm4
#define VFREDSUM_FLOAT vfredsumvs_float64xm4
#define VFMACCVV_FLOAT vfmaccvv_float64xm4
#define VFNMSACVV_FLOAT vfnmsacvv_float64xm4
#define VFMVVF_FLOAT vfmvvf_float64xm4
#define VFMULVV_FLOAT vfmulvv_float64xm4
#endif

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha_r, FLOAT alpha_i, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *buffer)
{
	BLASLONG i = 0, j = 0, k = 0;
	BLASLONG ix = 0, iy = 0;
	FLOAT *a_ptr = a;
        FLOAT temp_r, temp_i;

        FLOAT_V_T va0, va1, vx0, vx1, vr, vi;
        unsigned int gvl = 0;
        BLASLONG stride_x = inc_x * sizeof(FLOAT) * 2;
        BLASLONG stride_a = sizeof(FLOAT) * 2;
        gvl = vsetvli(m, RVV_EFLOAT, RVV_M);
        BLASLONG inc_xv = inc_x * gvl * 2;
        BLASLONG inc_av = gvl * 2;
        BLASLONG inc_y2 = inc_y * 2;
        BLASLONG lda2 = lda * 2;
        for(i = 0; i < n; i++){
                gvl = vsetvli(m, RVV_EFLOAT, RVV_M);
                j = 0;
                ix = 0;
                vr = VFMVVF_FLOAT(0, gvl);
                vi = VFMVVF_FLOAT(0, gvl);
                for(k = 0; k < m/gvl; k++){
                        va0 = VLSEV_FLOAT(&a_ptr[j], stride_a, gvl);
                        va1 = VLSEV_FLOAT(&a_ptr[j+1], stride_a, gvl);
                        vx0 = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                        vx1 = VLSEV_FLOAT(&x[ix+1], stride_x, gvl);
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
                        vr = VFMACCVV_FLOAT(vr, va0, vx0, gvl);
                        vr = VFNMSACVV_FLOAT(vr, va1, vx1, gvl);
                        vi = VFMACCVV_FLOAT(vi, va0, vx1, gvl);
                        vi = VFMACCVV_FLOAT(vi, va1, vx0, gvl);
#else
                        vr = VFMACCVV_FLOAT(vr, va0, vx0, gvl);
                        vr = VFMACCVV_FLOAT(vr, va1, vx1, gvl);
                        vi = VFMACCVV_FLOAT(vi, va0, vx1, gvl);
                        vi = VFNMSACVV_FLOAT(vi, va1, vx0, gvl);

#endif
                        j += inc_av;
                        ix += inc_xv;
                }
                va0 = VFMVVF_FLOAT(0, gvl);
                vx0 = VFREDSUM_FLOAT(vr, va0, gvl);
                temp_r = vx0[0];
                vx1 = VFREDSUM_FLOAT(vi, va0, gvl);
                temp_i = vx1[0];
                if(j/2 < m){
                        gvl = vsetvli(m-j/2, RVV_EFLOAT, RVV_M);
                        va0 = VLSEV_FLOAT(&a_ptr[j], stride_a, gvl);
                        va1 = VLSEV_FLOAT(&a_ptr[j+1], stride_a, gvl);
                        vx0 = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                        vx1 = VLSEV_FLOAT(&x[ix+1], stride_x, gvl);
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
                        vr = VFMULVV_FLOAT(va0, vx0, gvl);
                        vr = VFNMSACVV_FLOAT(vr, va1, vx1, gvl);
                        vi = VFMULVV_FLOAT(va0, vx1, gvl);
                        vi = VFMACCVV_FLOAT(vi, va1, vx0, gvl);
#else
                        vr = VFMULVV_FLOAT(va0, vx0, gvl);
                        vr = VFMACCVV_FLOAT(vr, va1, vx1, gvl);
                        vi = VFMULVV_FLOAT(va0, vx1, gvl);
                        vi = VFNMSACVV_FLOAT(vi, va1, vx0, gvl);

#endif
                        va0 = VFMVVF_FLOAT(0, gvl);
                        vx0 = VFREDSUM_FLOAT(vr, va0, gvl);
                        temp_r += vx0[0];
                        vx1 = VFREDSUM_FLOAT(vi, va0, gvl);
                        temp_i += vx1[0];
                }
#if !defined(XCONJ)
                y[iy]   += alpha_r * temp_r - alpha_i * temp_i;
                y[iy+1] += alpha_r * temp_i + alpha_i * temp_r;
#else
                y[iy]   += alpha_r * temp_r + alpha_i * temp_i;
                y[iy+1] -= alpha_r * temp_i - alpha_i * temp_r;
#endif
                iy += inc_y2;
                a_ptr += lda2;
        }
	return(0);
}

