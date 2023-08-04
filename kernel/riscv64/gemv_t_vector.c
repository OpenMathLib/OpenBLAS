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
#define VSETVL(n) vsetvl_e32m4(n)
#define VSETVL_MAX vsetvlmax_e32m1()
#define FLOAT_V_T vfloat32m4_t
#define FLOAT_V_T_M1 vfloat32m1_t
#define VFMVFS_FLOAT vfmv_f_s_f32m1_f32
#define VLEV_FLOAT vle32_v_f32m4
#define VLSEV_FLOAT vlse32_v_f32m4
#define VFREDSUM_FLOAT vfredosum_vs_f32m4_f32m1
#define VFMACCVV_FLOAT vfmacc_vv_f32m4
#define VFMVVF_FLOAT vfmv_v_f_f32m4
#define VFMVVF_FLOAT_M1 vfmv_v_f_f32m1
#define VFDOTVV_FLOAT vfdot_vv_f32m4
#define VFMULVV_FLOAT vfmul_vv_f32m4
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
#endif

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *buffer)
{
	    BLASLONG i = 0, j = 0, k = 0;
	    BLASLONG ix = 0, iy = 0;
	    FLOAT *a_ptr = a;
        FLOAT temp;

        FLOAT_V_T va, vr, vx;
        unsigned int gvl = 0;
        FLOAT_V_T_M1 v_res, v_z0;
        gvl = VSETVL_MAX;
        v_res = VFMVVF_FLOAT_M1(0, gvl);
        v_z0 = VFMVVF_FLOAT_M1(0, gvl);

        if(inc_x == 1){
                for(i = 0; i < n; i++){
                        gvl = VSETVL(m);
                        j = 0;
                        vr = VFMVVF_FLOAT(0, gvl);
                        for(k = 0; k < m/gvl; k++){
                                va = VLEV_FLOAT(&a_ptr[j], gvl);
                                vx = VLEV_FLOAT(&x[j], gvl);
                                vr = VFMACCVV_FLOAT(vr, va, vx, gvl);
                                j += gvl;
                        }
                        v_res = VFREDSUM_FLOAT(v_res, vr, v_z0, gvl);
                        temp = (FLOAT)VFMVFS_FLOAT(v_res);
                        if(j < m){
                                gvl = VSETVL(m-j);
                                va = VLEV_FLOAT(&a_ptr[j], gvl);
                                vx = VLEV_FLOAT(&x[j], gvl);
                                vr = VFMULVV_FLOAT(va, vx, gvl);

                                v_res = VFREDSUM_FLOAT(v_res, vr, v_z0, gvl);
                                temp += (FLOAT)VFMVFS_FLOAT(v_res);
                        }
                        y[iy] += alpha * temp;
                        iy += inc_y;
                        a_ptr += lda;
                }
        }else{
                BLASLONG stride_x = inc_x * sizeof(FLOAT);
                
                for(i = 0; i < n; i++){
                        gvl = VSETVL(m);
						BLASLONG inc_xv = inc_x * gvl;
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
                        v_res = VFREDSUM_FLOAT(v_res, vr, v_z0, gvl);
                        temp = (FLOAT)VFMVFS_FLOAT(v_res);
                        if(j < m){
                                gvl = VSETVL(m-j);
                                va = VLEV_FLOAT(&a_ptr[j], gvl);
                                vx = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                                vr = VFMULVV_FLOAT(va, vx, gvl);

                                v_res = VFREDSUM_FLOAT(v_res, vr, v_z0, gvl);
                                temp += (FLOAT)VFMVFS_FLOAT(v_res);
                        }
                        y[iy] += alpha * temp;
                        iy += inc_y;
                        a_ptr += lda;
                }
        }
	return(0);
}

