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
#define FLOAT_V_T vfloat32m4_t
#define VLEV_FLOAT vle32_v_f32m4
#define VLSEV_FLOAT vlse32_v_f32m4
#define VSEV_FLOAT vse32_v_f32m4
#define VSSEV_FLOAT vsse32_v_f32m4
#define VLSEG_FLOAT vlseg2e32_v_f32m4
#define VSSEG_FLOAT vsseg2e32_v_f32m4
#define VLSSEG_FLOAT vlsseg2e32_v_f32m4
#define VSSSEG_FLOAT vssseg2e32_v_f32m4
#define VFMACCVF_FLOAT vfmacc_vf_f32m4
#define VFNMSACVF_FLOAT vfnmsac_vf_f32m4
#else
#define VSETVL(n) vsetvl_e64m4(n)
#define FLOAT_V_T vfloat64m4_t
#define VLEV_FLOAT vle64_v_f64m4
#define VLSEV_FLOAT vlse64_v_f64m4
#define VSEV_FLOAT vse64_v_f64m4
#define VSSEV_FLOAT vsse64_v_f64m4
#define VLSEG_FLOAT vlseg2e64_v_f64m4
#define VSSEG_FLOAT vsseg2e64_v_f64m4
#define VLSSEG_FLOAT vlsseg2e64_v_f64m4
#define VSSSEG_FLOAT vssseg2e64_v_f64m4
#define VFMACCVF_FLOAT vfmacc_vf_f64m4
#define VFNMSACVF_FLOAT vfnmsac_vf_f64m4
#endif

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha_r, FLOAT alpha_i, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *buffer)
{
    BLASLONG i;
    BLASLONG ix;
    FLOAT *a_ptr;
    FLOAT temp_r, temp_i;
    FLOAT_V_T va0, va1, vy0, vy1;

    BLASLONG stride_y = inc_y * sizeof(FLOAT) * 2;

    BLASLONG inc_x2 = inc_x * 2;
    BLASLONG lda2 = lda * 2;
    if (inc_y == 1)
    {
        for (size_t vl; m > 0; m -= vl, a += vl*2, y += vl*2) {
            vl = VSETVL(m);
            a_ptr = a;
            ix = 0;
            VLSEG_FLOAT(&vy0, &vy1, y, vl);

            for(i = 0; i < n; i++){
#if !defined(XCONJ)
                temp_r = alpha_r * x[ix]   - alpha_i * x[ix+1];
                temp_i = alpha_r * x[ix+1] + alpha_i * x[ix];
#else
                temp_r = alpha_r * x[ix]   + alpha_i * x[ix+1];
                temp_i = alpha_r * x[ix+1] - alpha_i * x[ix];
#endif

                VLSEG_FLOAT(&va0, &va1, a_ptr, vl);
#if !defined(CONJ)
#if !defined(XCONJ)
                vy0 = VFMACCVF_FLOAT(vy0, temp_r, va0, vl);
                vy0 = VFNMSACVF_FLOAT(vy0, temp_i, va1, vl);
                vy1 = VFMACCVF_FLOAT(vy1, temp_r, va1, vl);
                vy1 = VFMACCVF_FLOAT(vy1, temp_i, va0, vl);
#else
                vy0 = VFMACCVF_FLOAT(vy0, temp_r, va0, vl);
                vy0 = VFMACCVF_FLOAT(vy0, temp_i, va1, vl);
                vy1 = VFMACCVF_FLOAT(vy1, temp_r, va1, vl);
                vy1 = VFNMSACVF_FLOAT(vy1, temp_i, va0, vl);
#endif
#else
#if !defined(XCONJ)
                vy0 = VFMACCVF_FLOAT(vy0, temp_r, va0, vl);
                vy0 = VFMACCVF_FLOAT(vy0, temp_i, va1, vl);
                vy1 = VFNMSACVF_FLOAT(vy1, temp_r, va1, vl);
                vy1 = VFMACCVF_FLOAT(vy1, temp_i, va0, vl);
#else
                vy0 = VFMACCVF_FLOAT(vy0, temp_r, va0, vl);
                vy0 = VFNMSACVF_FLOAT(vy0, temp_i, va1, vl);
                vy1 = VFNMSACVF_FLOAT(vy1, temp_r, va1, vl);
                vy1 = VFNMSACVF_FLOAT(vy1, temp_i, va0, vl);
#endif
#endif
                a_ptr += lda2;
                ix += inc_x2;
            }
            VSSEG_FLOAT(y, vy0, vy1, vl);
        }

    }
    else
    {
        for (size_t vl; m > 0; m -= vl, a += vl*2, y += vl*inc_y*2) {
            vl = VSETVL(m);
            a_ptr = a;
            ix = 0;
            VLSSEG_FLOAT(&vy0, &vy1, y, stride_y, vl);

            for(i = 0; i < n; i++){
#if !defined(XCONJ)
                temp_r = alpha_r * x[ix]   - alpha_i * x[ix+1];
                temp_i = alpha_r * x[ix+1] + alpha_i * x[ix];
#else
                temp_r = alpha_r * x[ix]   + alpha_i * x[ix+1];
                temp_i = alpha_r * x[ix+1] - alpha_i * x[ix];
#endif

                VLSEG_FLOAT(&va0, &va1, a_ptr, vl);
#if !defined(CONJ)
#if !defined(XCONJ)
                vy0 = VFMACCVF_FLOAT(vy0, temp_r, va0, vl);
                vy0 = VFNMSACVF_FLOAT(vy0, temp_i, va1, vl);
                vy1 = VFMACCVF_FLOAT(vy1, temp_r, va1, vl);
                vy1 = VFMACCVF_FLOAT(vy1, temp_i, va0, vl);
#else
                vy0 = VFMACCVF_FLOAT(vy0, temp_r, va0, vl);
                vy0 = VFMACCVF_FLOAT(vy0, temp_i, va1, vl);
                vy1 = VFMACCVF_FLOAT(vy1, temp_r, va1, vl);
                vy1 = VFNMSACVF_FLOAT(vy1, temp_i, va0, vl);
#endif
#else
#if !defined(XCONJ)
                vy0 = VFMACCVF_FLOAT(vy0, temp_r, va0, vl);
                vy0 = VFMACCVF_FLOAT(vy0, temp_i, va1, vl);
                vy1 = VFNMSACVF_FLOAT(vy1, temp_r, va1, vl);
                vy1 = VFMACCVF_FLOAT(vy1, temp_i, va0, vl);
#else
                vy0 = VFMACCVF_FLOAT(vy0, temp_r, va0, vl);
                vy0 = VFNMSACVF_FLOAT(vy0, temp_i, va1, vl);
                vy1 = VFNMSACVF_FLOAT(vy1, temp_r, va1, vl);
                vy1 = VFNMSACVF_FLOAT(vy1, temp_i, va0, vl);
#endif
#endif
                a_ptr += lda2;
                ix += inc_x2;
            }
            VSSSEG_FLOAT(y, stride_y, vy0, vy1, vl);
        }
    }
    return(0);
}
