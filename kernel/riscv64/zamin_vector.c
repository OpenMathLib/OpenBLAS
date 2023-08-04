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
#include <math.h>
#include <float.h>

#if !defined(DOUBLE)
#define VSETVL(n) vsetvl_e32m8(n)
#define VSETVL_MAX vsetvlmax_e32m1()
#define FLOAT_V_T vfloat32m8_t
#define FLOAT_V_T_M1 vfloat32m1_t
#define VFMVFS_FLOAT vfmv_f_s_f32m1_f32
#define VLSEV_FLOAT vlse32_v_f32m8
#define VFREDMINVS_FLOAT vfredmin_vs_f32m8_f32m1
#define MASK_T vbool4_t
#define VMFLTVF_FLOAT vmflt_vf_f32m8_b4
#define VFMVVF_FLOAT vfmv_v_f_f32m8
#define VFMVVF_FLOAT_M1 vfmv_v_f_f32m1
#define VFRSUBVF_MASK_FLOAT vfrsub_vf_f32m8_m
#define VFMINVV_FLOAT vfmin_vv_f32m8
#define VFADDVV_FLOAT vfadd_vv_f32m8
#else
#define VSETVL(n) vsetvl_e64m8(n)
#define VSETVL_MAX vsetvlmax_e32m1()
#define FLOAT_V_T vfloat64m8_t
#define FLOAT_V_T_M1 vfloat64m1_t
#define VFMVFS_FLOAT vfmv_f_s_f64m1_f64
#define VLSEV_FLOAT vlse64_v_f64m8
#define VFREDMINVS_FLOAT vfredmin_vs_f64m8_f64m1
#define MASK_T vbool8_t
#define VMFLTVF_FLOAT vmflt_vf_f64m8_b8
#define VFMVVF_FLOAT vfmv_v_f_f64m8
#define VFMVVF_FLOAT_M1 vfmv_v_f_f64m1
#define VFRSUBVF_MASK_FLOAT vfrsub_vf_f64m8_m
#define VFMINVV_FLOAT vfmin_vv_f64m8
#define VFADDVV_FLOAT vfadd_vv_f64m8
#endif

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
	BLASLONG i=0, j=0;
	BLASLONG ix=0;
	if (n <= 0 || inc_x <= 0) return(0.0);
	FLOAT minf=FLT_MAX;
        unsigned int gvl = 0;
        FLOAT_V_T v0, v1, v_min;
        FLOAT_V_T_M1 v_res, v_max;
        gvl = VSETVL_MAX;
        v_res = VFMVVF_FLOAT_M1(0, gvl);
        v_max = VFMVVF_FLOAT_M1(FLT_MAX, gvl);

        MASK_T mask0, mask1;
        BLASLONG stride_x = inc_x * sizeof(FLOAT) * 2;
        gvl = VSETVL(n);
        v_min = VFMVVF_FLOAT(FLT_MAX, gvl);
        BLASLONG inc_xv = inc_x * gvl * 2;
        for(; i<n/gvl; i++){
                v0 = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                v1 = VLSEV_FLOAT(&x[ix+1], stride_x, gvl);
                mask0 = VMFLTVF_FLOAT(v0, 0, gvl);
                v0 = VFRSUBVF_MASK_FLOAT(mask0, v0, v0, 0, gvl);
                mask1 = VMFLTVF_FLOAT(v1, 0, gvl);
                v1 = VFRSUBVF_MASK_FLOAT(mask1, v1, v1, 0, gvl);

                v0 = VFADDVV_FLOAT(v0, v1, gvl);
                v_min = VFMINVV_FLOAT(v_min, v0, gvl);

                j += gvl;
                ix += inc_xv;
        }
        v_res = VFREDMINVS_FLOAT(v_res, v_min, v_max, gvl);
        minf = VFMVFS_FLOAT(v_res);

        if(j<n){
                gvl = VSETVL(n-j);
                v0 = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                v1 = VLSEV_FLOAT(&x[ix+1], stride_x, gvl);
                mask0 = VMFLTVF_FLOAT(v0, 0, gvl);
                v0 = VFRSUBVF_MASK_FLOAT(mask0, v0, v0, 0, gvl);
                mask1 = VMFLTVF_FLOAT(v1, 0, gvl);
                v1 = VFRSUBVF_MASK_FLOAT(mask1, v1, v1, 0, gvl);
                v1 = VFADDVV_FLOAT(v0, v1, gvl);
                v_res = VFREDMINVS_FLOAT(v_res, v1, v_max, gvl);
                if(VFMVFS_FLOAT(v_res) < minf)
                        minf = VFMVFS_FLOAT(v_res);
        }
        return(minf);
}
