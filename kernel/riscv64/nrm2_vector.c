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
#define VSETVL_MAX vsetvlmax_e32m1()
#define FLOAT_V_T vfloat32m4_t
#define VFMVFS_FLOATM4 vfmv_f_s_f32m4_f32
#define FLOAT_V_T_M1 vfloat32m1_t
#define VFMVFS_FLOAT vfmv_f_s_f32m1_f32
#define VLEV_FLOAT vle32_v_f32m4
#define VLSEV_FLOAT vlse32_v_f32m4
#define VFREDSUM_FLOAT vfredusum_vs_f32m4_f32m1
#define VFMACCVV_FLOAT vfmacc_vv_f32m4
#define VFMVVF_FLOAT vfmv_v_f_f32m4
#define VFMVVF_FLOAT_M1 vfmv_v_f_f32m1
#define VFDOTVV_FLOAT vfdot_vv_f32m4
#define ABS fabsf
#define MASK_T vbool8_t
#define VFRSUBVF_MASK_FLOAT vfrsub_vf_f32m4_m
#define VMFGTVF_FLOAT vmfgt_vf_f32m4_b8
#define VMFIRSTM vmfirst_m_b8
#define VFDIVVF_FLOAT vfdiv_vf_f32m4
#define VMFLTVF_FLOAT vmflt_vf_f32m4_b8
#define VFREDMAXVS_FLOAT vfredmax_vs_f32m4_f32m1
#else
#define VSETVL(n) vsetvl_e64m4(n)
#define VSETVL_MAX vsetvlmax_e64m1()
#define FLOAT_V_T vfloat64m4_t
#define VFMVFS_FLOATM4 vfmv_f_s_f64m4_f64
#define FLOAT_V_T_M1 vfloat64m1_t
#define VFMVFS_FLOAT vfmv_f_s_f64m1_f64
#define VLEV_FLOAT vle64_v_f64m4
#define VLSEV_FLOAT vlse64_v_f64m4
#define VFREDSUM_FLOAT vfredusum_vs_f64m4_f64m1
#define VFMACCVV_FLOAT vfmacc_vv_f64m4
#define VFMVVF_FLOAT vfmv_v_f_f64m4
#define VFMVVF_FLOAT_M1 vfmv_v_f_f64m1
#define VFDOTVV_FLOAT vfdot_vv_f64m4
#define ABS fabs
#define MASK_T vbool16_t
#define VFRSUBVF_MASK_FLOAT vfrsub_vf_f64m4_m
#define VMFGTVF_FLOAT vmfgt_vf_f64m4_b16
#define VMFIRSTM vmfirst_m_b16
#define VFDIVVF_FLOAT vfdiv_vf_f64m4
#define VMFLTVF_FLOAT vmflt_vf_f64m4_b16
#define VFREDMAXVS_FLOAT vfredmax_vs_f64m4_f64m1
#endif

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
	BLASLONG i=0, j=0;

	if ( n < 0 )  return(0.0);
        if(n == 1) return (ABS(x[0]));

        FLOAT_V_T vr, v0, v_zero;
        unsigned int gvl = 0;
        FLOAT_V_T_M1 v_res, v_z0;
        gvl = VSETVL_MAX;
        v_res = VFMVVF_FLOAT_M1(0, gvl);
        v_z0 = VFMVVF_FLOAT_M1(0, gvl);

        FLOAT scale = 0.0, ssq = 0.0;
        MASK_T mask;
        BLASLONG index = 0;
        if(inc_x == 1){
                gvl = VSETVL(n);
                vr = VFMVVF_FLOAT(0, gvl);
                v_zero = VFMVVF_FLOAT(0, gvl);
                for(i=0,j=0; i<n/gvl; i++){
                        v0 = VLEV_FLOAT(&x[j], gvl);
                        //fabs(vector)
                        mask = VMFLTVF_FLOAT(v0, 0, gvl);
                        v0 = VFRSUBVF_MASK_FLOAT(mask, v0, v0, 0, gvl);
                        //if scale change
                        mask = VMFGTVF_FLOAT(v0, scale, gvl);
                        index = VMFIRSTM(mask, gvl);
                        if(index == -1){//no elements greater than scale
                                if(scale != 0.0){
                                        v0 = VFDIVVF_FLOAT(v0, scale, gvl);
                                        vr = VFMACCVV_FLOAT(vr, v0, v0, gvl);
                                }
                        }else{//found greater element
                                //ssq in vector vr: vr[0]
                                v_res = VFREDSUM_FLOAT(v_res, vr, v_z0, gvl);
                                //total ssq before current vector
                                ssq += VFMVFS_FLOAT(v_res);
                                //find max
                                v_res = VFREDMAXVS_FLOAT(v_res, v0, v_z0, gvl);
                                //update ssq before max_index
                                ssq = ssq * (scale/VFMVFS_FLOAT(v_res))*(scale/VFMVFS_FLOAT(v_res));
                                //update scale
                                scale = VFMVFS_FLOAT(v_res);
                                //ssq in vector vr
                                v0 = VFDIVVF_FLOAT(v0, scale, gvl);
                                vr = VFMACCVV_FLOAT(v_zero, v0, v0, gvl);
                        }
                        j += gvl;
                }
                //ssq in vector vr: vr[0]
                v_res = VFREDSUM_FLOAT(v_res, vr, v_z0, gvl);
                //total ssq now
                ssq += VFMVFS_FLOAT(v_res);

                //tail
                if(j < n){
                        gvl = VSETVL(n-j);
                        v0 = VLEV_FLOAT(&x[j], gvl);
                        //fabs(vector)
                        mask = VMFLTVF_FLOAT(v0, 0, gvl);
                        v0 = VFRSUBVF_MASK_FLOAT(mask, v0, v0, 0, gvl);
                        //if scale change
                        mask = VMFGTVF_FLOAT(v0, scale, gvl);
                        index = VMFIRSTM(mask, gvl);
                        if(index == -1){//no elements greater than scale
                                if(scale != 0.0)
                                        v0 = VFDIVVF_FLOAT(v0, scale, gvl);
                        }else{//found greater element
                                //find max
                                v_res = VFREDMAXVS_FLOAT(v_res, v0, v_z0, gvl);
                                //update ssq before max_index
                                ssq = ssq * (scale/VFMVFS_FLOAT(v_res))*(scale/VFMVFS_FLOAT(v_res));
                                //update scale
                                scale = VFMVFS_FLOAT(v_res);
                                v0 = VFDIVVF_FLOAT(v0, scale, gvl);
                        }
                        vr = VFMACCVV_FLOAT(v_zero, v0, v0, gvl);
                        //ssq in vector vr: vr[0]
                        v_res = VFREDSUM_FLOAT(v_res, vr, v_z0, gvl);
                        //total ssq now
                        ssq += VFMVFS_FLOAT(v_res);
                }
        }else{
                gvl = VSETVL(n);
                vr = VFMVVF_FLOAT(0, gvl);
                v_zero = VFMVVF_FLOAT(0, gvl);
                unsigned int stride_x = inc_x * sizeof(FLOAT);
                int idx = 0, inc_v = inc_x * gvl;
                for(i=0,j=0; i<n/gvl; i++){
                        v0 = VLSEV_FLOAT(&x[idx], stride_x, gvl);
                        //fabs(vector)
                        mask = VMFLTVF_FLOAT(v0, 0, gvl);
                        v0 = VFRSUBVF_MASK_FLOAT(mask, v0, v0, 0, gvl);
                        //if scale change
                        mask = VMFGTVF_FLOAT(v0, scale, gvl);
                        index = VMFIRSTM(mask, gvl);
                        if(index == -1){//no elements greater than scale
                                if(scale != 0.0){
                                        v0 = VFDIVVF_FLOAT(v0, scale, gvl);
                                        vr = VFMACCVV_FLOAT(vr, v0, v0, gvl);
                                }
                        }else{//found greater element
                                //ssq in vector vr: vr[0]
                                v_res = VFREDSUM_FLOAT(v_res, vr, v_z0, gvl);
                                //total ssq before current vector
                                ssq += VFMVFS_FLOAT(v_res);
                                //find max
                                v_res = VFREDMAXVS_FLOAT(v_res, v0, v_z0, gvl);
                                //update ssq before max_index
                                ssq = ssq * (scale/VFMVFS_FLOAT(v_res))*(scale/VFMVFS_FLOAT(v_res));
                                //update scale
                                scale = VFMVFS_FLOAT(v_res);
                                //ssq in vector vr
                                v0 = VFDIVVF_FLOAT(v0, scale, gvl);
                                vr = VFMACCVV_FLOAT(v_zero, v0, v0, gvl);
                        }
                        j += gvl;
                        idx += inc_v;
                }
                //ssq in vector vr: vr[0]
                v_res = VFREDSUM_FLOAT(v_res, vr, v_z0, gvl);
                //total ssq now
                ssq += VFMVFS_FLOAT(v_res);

                //tail
                if(j < n){
                        gvl = VSETVL(n-j);
                        v0 = VLSEV_FLOAT(&x[idx], stride_x, gvl);
                        //fabs(vector)
                        mask = VMFLTVF_FLOAT(v0, 0, gvl);
                        v0 = VFRSUBVF_MASK_FLOAT(mask, v0, v0, 0, gvl);
                        //if scale change
                        mask = VMFGTVF_FLOAT(v0, scale, gvl);
                        index = VMFIRSTM(mask, gvl);
                        if(index == -1){//no elements greater than scale
                                if(scale != 0.0)
                                        v0 = VFDIVVF_FLOAT(v0, scale, gvl);
                        }else{//found greater element
                                //find max
                                v_res = VFREDMAXVS_FLOAT(v_res, v0, v_z0, gvl);
                                //update ssq before max_index
                                ssq = ssq * (scale/VFMVFS_FLOAT(v_res))*(scale/VFMVFS_FLOAT(v_res));
                                //update scale
                                scale = VFMVFS_FLOATM4(vr);
                                v0 = VFDIVVF_FLOAT(v0, scale, gvl);
                        }
                        vr = VFMACCVV_FLOAT(v_zero, v0, v0, gvl);
                        //ssq in vector vr: vr[0]
                        v_res = VFREDSUM_FLOAT(v_res, vr, v_z0, gvl);
                        //total ssq now
                        ssq += VFMVFS_FLOAT(v_res);
                }
        }
	return(scale * sqrt(ssq));
}


