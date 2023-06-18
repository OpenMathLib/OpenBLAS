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
#define FLOAT_V_T vfloat32m4_t
#define VLEV_FLOAT vle32_v_f32m4
#define VLSEV_FLOAT vlse32_v_f32m4
#define VSEV_FLOAT vse32_v_f32m4
#define VSSEV_FLOAT vsse32_v_f32m4
#define VFMACCVF_FLOAT vfmacc_vf_f32m4
#define VFNMSACVF_FLOAT vfnmsac_vf_f32m4
#else
#define VSETVL(n) vsetvl_e64m4(n)
#define FLOAT_V_T vfloat64m4_t
#define VLEV_FLOAT vle64_v_f64m4
#define VLSEV_FLOAT vlse64_v_f64m4
#define VSEV_FLOAT vse64_v_f64m4
#define VSSEV_FLOAT vsse64_v_f64m4
#define VFMACCVF_FLOAT vfmacc_vf_f64m4
#define VFNMSACVF_FLOAT vfnmsac_vf_f64m4
#endif

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha_r, FLOAT alpha_i, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *buffer)
{
	BLASLONG i = 0, j = 0, k = 0;
        BLASLONG ix = 0, iy = 0;
        FLOAT *a_ptr = a;
        FLOAT temp_r = 0.0, temp_i = 0.0;
        FLOAT_V_T va0, va1, vy0, vy1;
        unsigned int gvl = 0;
        BLASLONG stride_a = sizeof(FLOAT) * 2;
        BLASLONG stride_y = inc_y * sizeof(FLOAT) * 2;
        gvl = VSETVL(m);
        BLASLONG inc_yv = inc_y * gvl * 2;
        BLASLONG inc_x2 = inc_x * 2;
        BLASLONG lda2 = lda * 2;
        for(k=0,j=0; k<m/gvl; k++){
                a_ptr = a;
                ix = 0;
                vy0 = VLSEV_FLOAT(&y[iy], stride_y, gvl);
                vy1 = VLSEV_FLOAT(&y[iy+1], stride_y, gvl);
                for(i = 0; i < n; i++){
#if !defined(XCONJ)
			temp_r = alpha_r * x[ix]   - alpha_i * x[ix+1];
			temp_i = alpha_r * x[ix+1] + alpha_i * x[ix];
#else
			temp_r = alpha_r * x[ix]   + alpha_i * x[ix+1];
			temp_i = alpha_r * x[ix+1] - alpha_i * x[ix];
#endif

                        va0 = VLSEV_FLOAT(&a_ptr[j], stride_a, gvl);
                        va1 = VLSEV_FLOAT(&a_ptr[j+1], stride_a, gvl);
#if !defined(CONJ)
#if !defined(XCONJ)
			vy0 = VFMACCVF_FLOAT(vy0, temp_r, va0, gvl);
			vy0 = VFNMSACVF_FLOAT(vy0, temp_i, va1, gvl);
			vy1 = VFMACCVF_FLOAT(vy1, temp_r, va1, gvl);
			vy1 = VFMACCVF_FLOAT(vy1, temp_i, va0, gvl);
#else

			vy0 = VFMACCVF_FLOAT(vy0, temp_r, va0, gvl);
			vy0 = VFMACCVF_FLOAT(vy0, temp_i, va1, gvl);
			vy1 = VFMACCVF_FLOAT(vy1, temp_r, va1, gvl);
			vy1 = VFNMSACVF_FLOAT(vy1, temp_i, va0, gvl);
#endif

#else

#if !defined(XCONJ)
			vy0 = VFMACCVF_FLOAT(vy0, temp_r, va0, gvl);
			vy0 = VFMACCVF_FLOAT(vy0, temp_i, va1, gvl);
			vy1 = VFNMSACVF_FLOAT(vy1, temp_r, va1, gvl);
			vy1 = VFMACCVF_FLOAT(vy1, temp_i, va0, gvl);
#else
			vy0 = VFMACCVF_FLOAT(vy0, temp_r, va0, gvl);
			vy0 = VFNMSACVF_FLOAT(vy0, temp_i, va1, gvl);
			vy1 = VFNMSACVF_FLOAT(vy1, temp_r, va1, gvl);
			vy1 = VFNMSACVF_FLOAT(vy1, temp_i, va0, gvl);
#endif

#endif
                        a_ptr += lda2;
                        ix += inc_x2;
                }
                VSSEV_FLOAT(&y[iy], stride_y, vy0, gvl);
                VSSEV_FLOAT(&y[iy+1], stride_y, vy1, gvl);
                j += gvl * 2;
                iy += inc_yv;
        }
        //tail
        if(j/2 < m){
                gvl = VSETVL(m-j/2);
                a_ptr = a;
                ix = 0;
                vy0 = VLSEV_FLOAT(&y[iy], stride_y, gvl);
                vy1 = VLSEV_FLOAT(&y[iy+1], stride_y, gvl);
                for(i = 0; i < n; i++){
#if !defined(XCONJ)
			temp_r = alpha_r * x[ix]   - alpha_i * x[ix+1];
			temp_i = alpha_r * x[ix+1] + alpha_i * x[ix];
#else
			temp_r = alpha_r * x[ix]   + alpha_i * x[ix+1];
			temp_i = alpha_r * x[ix+1] - alpha_i * x[ix];
#endif

                        va0 = VLSEV_FLOAT(&a_ptr[j], stride_a, gvl);
                        va1 = VLSEV_FLOAT(&a_ptr[j+1], stride_a, gvl);
#if !defined(CONJ)

#if !defined(XCONJ)
			vy0 = VFMACCVF_FLOAT(vy0, temp_r, va0, gvl);
			vy0 = VFNMSACVF_FLOAT(vy0, temp_i, va1, gvl);
			vy1 = VFMACCVF_FLOAT(vy1, temp_r, va1, gvl);
			vy1 = VFMACCVF_FLOAT(vy1, temp_i, va0, gvl);
#else

			vy0 = VFMACCVF_FLOAT(vy0, temp_r, va0, gvl);
			vy0 = VFMACCVF_FLOAT(vy0, temp_i, va1, gvl);
			vy1 = VFMACCVF_FLOAT(vy1, temp_r, va1, gvl);
			vy1 = VFNMSACVF_FLOAT(vy1, temp_i, va0, gvl);
#endif

#else

#if !defined(XCONJ)
			vy0 = VFMACCVF_FLOAT(vy0, temp_r, va0, gvl);
			vy0 = VFMACCVF_FLOAT(vy0, temp_i, va1, gvl);
			vy1 = VFNMSACVF_FLOAT(vy1, temp_r, va1, gvl);
			vy1 = VFMACCVF_FLOAT(vy1, temp_i, va0, gvl);
#else
			vy0 = VFMACCVF_FLOAT(vy0, temp_r, va0, gvl);
			vy0 = VFNMSACVF_FLOAT(vy0, temp_i, va1, gvl);
			vy1 = VFNMSACVF_FLOAT(vy1, temp_r, va1, gvl);
			vy1 = VFNMSACVF_FLOAT(vy1, temp_i, va0, gvl);
#endif

#endif
                        a_ptr += lda2;
                        ix += inc_x2;
                }
                VSSEV_FLOAT(&y[iy], stride_y, vy0, gvl);
                VSSEV_FLOAT(&y[iy+1], stride_y, vy1, gvl);
        }
	return(0);
}


