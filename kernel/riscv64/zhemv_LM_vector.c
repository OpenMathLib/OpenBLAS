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
#define VSSEV_FLOAT vssev_float32xm4
#define VFREDSUM_FLOAT vfredsumvs_float32xm4
#define VFMACCVV_FLOAT vfmaccvv_float32xm4
#define VFMACCVF_FLOAT vfmaccvf_float32xm4
#define VFMVVF_FLOAT vfmvvf_float32xm4
#define VFMULVV_FLOAT vfmulvv_float32xm4
#define VFNMSACVF_FLOAT vfnmsacvf_float32xm4
#define VFNMSACVV_FLOAT vfnmsacvv_float32xm4
#else
#define RVV_EFLOAT RVV_E64
#define RVV_M RVV_M4
#define FLOAT_V_T float64xm4_t
#define VLSEV_FLOAT vlsev_float64xm4
#define VSSEV_FLOAT vssev_float64xm4
#define VFREDSUM_FLOAT vfredsumvs_float64xm4
#define VFMACCVV_FLOAT vfmaccvv_float64xm4
#define VFMACCVF_FLOAT vfmaccvf_float64xm4
#define VFMVVF_FLOAT vfmvvf_float64xm4
#define VFMULVV_FLOAT vfmulvv_float64xm4
#define VFNMSACVF_FLOAT vfnmsacvf_float64xm4
#define VFNMSACVV_FLOAT vfnmsacvv_float64xm4
#endif

int CNAME(BLASLONG m, BLASLONG offset, FLOAT alpha_r, FLOAT alpha_i, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG incx, FLOAT *y, BLASLONG incy, FLOAT *buffer){
        BLASLONG i, j, k;
        BLASLONG ix, iy, ia;
        BLASLONG jx, jy, ja;
        FLOAT temp_r1, temp_i1;
        FLOAT temp_r2, temp_i2;
        FLOAT *a_ptr = a;
        unsigned int gvl = 0;


        FLOAT_V_T va0, va1, vx0, vx1, vy0, vy1, vr0, vr1;
        BLASLONG stride_x, stride_y, stride_a, inc_xv, inc_yv, inc_av, len, lda2;

        BLASLONG inc_x2 = incx * 2;
        BLASLONG inc_y2 = incy * 2;
        stride_x = inc_x2 * sizeof(FLOAT);
        stride_y = inc_y2 * sizeof(FLOAT);
        stride_a = 2 * sizeof(FLOAT);
        lda2 = lda * 2;

        jx = 0;
        jy = 0;
        ja = 0;
        for(j = 0; j < offset; j++){
                temp_r1 = alpha_r * x[jx]   - alpha_i * x[jx+1];;
                temp_i1 = alpha_r * x[jx+1] + alpha_i * x[jx];
                temp_r2 = 0;
                temp_i2 = 0;
                y[jy] += temp_r1 * a_ptr[ja];
                y[jy+1] += temp_i1 * a_ptr[ja];
                ix = jx + inc_x2;
                iy = jy + inc_y2;
                ia = ja + 2;
                i = j + 1;
                len = m - i;
                if(len > 0){
                        gvl = vsetvli(len, RVV_EFLOAT, RVV_M);
                        inc_xv = incx * gvl * 2;
                        inc_yv = incy * gvl * 2;
                        inc_av = gvl * 2;
                        vr0 = VFMVVF_FLOAT(0, gvl);
                        vr1 = VFMVVF_FLOAT(0, gvl);
                        for(k = 0; k < len / gvl; k++){
                                va0 = VLSEV_FLOAT(&a_ptr[ia], stride_a, gvl);
                                va1 = VLSEV_FLOAT(&a_ptr[ia+1], stride_a, gvl);
                                vy0 = VLSEV_FLOAT(&y[iy], stride_y, gvl);
                                vy1 = VLSEV_FLOAT(&y[iy+1], stride_y, gvl);
#ifndef HEMVREV
                                vy0 = VFMACCVF_FLOAT(vy0, temp_r1, va0, gvl);
                                vy0 = VFNMSACVF_FLOAT(vy0, temp_i1, va1, gvl);
                                vy1 = VFMACCVF_FLOAT(vy1, temp_r1, va1, gvl);
                                vy1 = VFMACCVF_FLOAT(vy1, temp_i1, va0, gvl);
#else
                                vy0 = VFMACCVF_FLOAT(vy0, temp_r1, va0, gvl);
                                vy0 = VFMACCVF_FLOAT(vy0, temp_i1, va1, gvl);
                                vy1 = VFNMSACVF_FLOAT(vy1, temp_r1, va1, gvl);
                                vy1 = VFMACCVF_FLOAT(vy1, temp_i1, va0, gvl);
#endif
                                VSSEV_FLOAT(&y[iy], stride_y, vy0, gvl);
                                VSSEV_FLOAT(&y[iy+1], stride_y, vy1, gvl);

                                vx0 = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                                vx1 = VLSEV_FLOAT(&x[ix+1], stride_x, gvl);
#ifndef HEMVREV
                                vr0 = VFMACCVV_FLOAT(vr0, vx0, va0, gvl);
                                vr0 = VFMACCVV_FLOAT(vr0, vx1, va1, gvl);
                                vr1 = VFMACCVV_FLOAT(vr1, vx1, va0, gvl);
                                vr1 = VFNMSACVV_FLOAT(vr1, vx0, va1, gvl);
#else
                                vr0 = VFMACCVV_FLOAT(vr0, vx0, va0, gvl);
                                vr0 = VFNMSACVV_FLOAT(vr0, vx1, va1, gvl);
                                vr1 = VFMACCVV_FLOAT(vr1, vx1, va0, gvl);
                                vr1 = VFMACCVV_FLOAT(vr1, vx0, va1, gvl);

#endif
                                i += gvl;
                                ix += inc_xv;
                                iy += inc_yv;
                                ia += inc_av;
                        }
                        va0 = VFMVVF_FLOAT(0, gvl);
                        vx0 = VFREDSUM_FLOAT(vr0, va0, gvl);
                        temp_r2 = vx0[0];
                        vx1 = VFREDSUM_FLOAT(vr1, va0, gvl);
                        temp_i2 = vx1[0];
                        if(i < m){
				gvl = vsetvli(m-i, RVV_EFLOAT, RVV_M);
                                va0 = VLSEV_FLOAT(&a_ptr[ia], stride_a, gvl);
                                va1 = VLSEV_FLOAT(&a_ptr[ia+1], stride_a, gvl);
                                vy0 = VLSEV_FLOAT(&y[iy], stride_y, gvl);
                                vy1 = VLSEV_FLOAT(&y[iy+1], stride_y, gvl);
#ifndef HEMVREV
                                vy0 = VFMACCVF_FLOAT(vy0, temp_r1, va0, gvl);
                                vy0 = VFNMSACVF_FLOAT(vy0, temp_i1, va1, gvl);
                                vy1 = VFMACCVF_FLOAT(vy1, temp_r1, va1, gvl);
                                vy1 = VFMACCVF_FLOAT(vy1, temp_i1, va0, gvl);
#else
                                vy0 = VFMACCVF_FLOAT(vy0, temp_r1, va0, gvl);
                                vy0 = VFMACCVF_FLOAT(vy0, temp_i1, va1, gvl);
                                vy1 = VFNMSACVF_FLOAT(vy1, temp_r1, va1, gvl);
                                vy1 = VFMACCVF_FLOAT(vy1, temp_i1, va0, gvl);
#endif
                                VSSEV_FLOAT(&y[iy], stride_y, vy0, gvl);
                                VSSEV_FLOAT(&y[iy+1], stride_y, vy1, gvl);

                                vx0 = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                                vx1 = VLSEV_FLOAT(&x[ix+1], stride_x, gvl);
#ifndef HEMVREV
                                vr0 = VFMULVV_FLOAT(vx0, va0, gvl);
                                vr0 = VFMACCVV_FLOAT(vr0, vx1, va1, gvl);
                                vr1 = VFMULVV_FLOAT(vx1, va0, gvl);
                                vr1 = VFNMSACVV_FLOAT(vr1, vx0, va1, gvl);
#else
                                vr0 = VFMULVV_FLOAT(vx0, va0, gvl);
                                vr0 = VFNMSACVV_FLOAT(vr0, vx1, va1, gvl);
                                vr1 = VFMULVV_FLOAT(vx1, va0, gvl);
                                vr1 = VFMACCVV_FLOAT(vr1, vx0, va1, gvl);
#endif

                                va0 = VFMVVF_FLOAT(0, gvl);
                                vx0 = VFREDSUM_FLOAT(vr0, va0, gvl);
                                temp_r2 += vx0[0];
                                vx1 = VFREDSUM_FLOAT(vr1, va0, gvl);
                                temp_i2 += vx1[0];
                        }
                }
		y[jy] += alpha_r * temp_r2 - alpha_i * temp_i2;
		y[jy+1] += alpha_r * temp_i2 + alpha_i * temp_r2;
		jx    += inc_x2;
		jy    += inc_y2;
		ja    += 2;
		a_ptr += lda2;
        }
	return(0);
}
