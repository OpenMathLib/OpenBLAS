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

#if defined(DOUBLE)
#define VFMVFS_FLOAT vfmv_f_s_f64m1_f64
#define VSETVL(n) vsetvl_e64m8(n)
#define VSETVL_MAX vsetvlmax_e64m1()
#define FLOAT_V_T vfloat64m8_t
#define FLOAT_V_T_M1 vfloat64m1_t
#define VLEV_FLOAT vle64_v_f64m8
#define VLSEV_FLOAT vlse64_v_f64m8
#define VFREDMAXVS_FLOAT vfredmax_vs_f64m8_f64m1
#define MASK_T vbool8_t
#define VMFLTVF_FLOAT vmflt_vf_f64m8_b8
#define VMFLTVV_FLOAT vmflt_vv_f64m8_b8
#define VFMVVF_FLOAT vfmv_v_f_f64m8
#define VFMVVF_FLOAT_M1 vfmv_v_f_f64m1
#define VFRSUBVF_MASK_FLOAT vfrsub_vf_f64m8_m
#define VFMAXVV_FLOAT vfmax_vv_f64m8
#define VMFGEVF_FLOAT vmfge_vf_f64m8_b8
#define VMFIRSTM vmfirst_m_b8
#define UINT_V_T vuint64m8_t
#define VIDV_MASK_UINT vid_v_u64m8_m
#define VIDV_UINT vid_v_u64m8
#define VADDVX_MASK_UINT vadd_vx_u64m8_m
#define VADDVX_UINT vadd_vx_u64m8
#define VMVVX_UINT vmv_v_x_u64m8
#else
#define VFMVFS_FLOAT vfmv_f_s_f32m1_f32
#define VSETVL(n) vsetvl_e32m8(n)
#define VSETVL_MAX vsetvlmax_e32m1()
#define FLOAT_V_T vfloat32m8_t
#define FLOAT_V_T_M1 vfloat32m1_t
#define VLEV_FLOAT vle32_v_f32m8
#define VLSEV_FLOAT vlse32_v_f32m8
#define VFREDMAXVS_FLOAT vfredmax_vs_f32m8_f32m1
#define MASK_T vbool4_t
#define VMFLTVF_FLOAT vmflt_vf_f32m8_b4
#define VMFLTVV_FLOAT vmflt_vv_f32m8_b4
#define VFMVVF_FLOAT vfmv_v_f_f32m8
#define VFMVVF_FLOAT_M1 vfmv_v_f_f32m1
#define VFRSUBVF_MASK_FLOAT vfrsub_vf_f32m8_m
#define VFMAXVV_FLOAT vfmax_vv_f32m8
#define VMFGEVF_FLOAT vmfge_vf_f32m8_b4
#define VMFIRSTM vmfirst_m_b4
#define UINT_V_T vuint32m8_t
#define VIDV_MASK_UINT vid_v_u32m8_m
#define VIDV_UINT vid_v_u32m8
#define VADDVX_MASK_UINT vadd_vx_u32m8_m
#define VADDVX_UINT vadd_vx_u32m8
#define VMVVX_UINT vmv_v_x_u32m8
#endif


BLASLONG CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
	BLASLONG i=0, j=0;
	FLOAT maxf=0.0;
#ifdef DOUBLE
        BLASLONG max_index = 0;
#else
        unsigned int max_index = 0;
#endif
	if (n <= 0 || inc_x <= 0) return(max_index);

        FLOAT_V_T vx, v_max;
        UINT_V_T v_max_index;
        MASK_T mask;
        unsigned int gvl = 0;
        FLOAT_V_T_M1 v_res, v_z0;
        gvl = VSETVL_MAX;
        v_res = VFMVVF_FLOAT_M1(0, gvl);
        v_z0 = VFMVVF_FLOAT_M1(0, gvl);

        if(inc_x == 1){
                gvl = VSETVL(n);
                v_max_index = VMVVX_UINT(0, gvl);
                v_max = VFMVVF_FLOAT(-1, gvl);
                for(i=0,j=0; i < n/gvl; i++){
                        vx = VLEV_FLOAT(&x[j], gvl);
                        //fabs(vector)
                        mask = VMFLTVF_FLOAT(vx, 0, gvl);
                        vx = VFRSUBVF_MASK_FLOAT(mask, vx, vx, 0, gvl);

                        //index where element greater than v_max
                        mask = VMFLTVV_FLOAT(v_max, vx, gvl);
                        v_max_index = VIDV_MASK_UINT(mask, v_max_index, gvl);
                        v_max_index = VADDVX_MASK_UINT(mask, v_max_index, v_max_index, j,gvl);

                        //update v_max and start_index j
                        v_max = VFMAXVV_FLOAT(v_max, vx, gvl);
                        j += gvl;
                }
                v_res = VFREDMAXVS_FLOAT(v_res, v_max, v_z0, gvl);
                maxf = VFMVFS_FLOAT(v_res);
                mask = VMFGEVF_FLOAT(v_max, maxf, gvl);
                max_index = VMFIRSTM(mask,gvl);
#ifdef DOUBLE
		max_index = *((BLASLONG *)&v_max_index+max_index);
#else
		max_index = *((unsigned int *)&v_max_index+max_index);
#endif
                if(j < n){
                        gvl = VSETVL(n-j);
                        vx = VLEV_FLOAT(&x[j], gvl);
                        //fabs(vector)
                        mask = VMFLTVF_FLOAT(vx, 0, gvl);
                        v_max = VFRSUBVF_MASK_FLOAT(mask, vx, vx, 0, gvl);

                        v_res = VFREDMAXVS_FLOAT(v_res, v_max, v_z0, gvl);
                        FLOAT cur_maxf = VFMVFS_FLOAT(v_res);
                        if(cur_maxf > maxf){
                                //tail index
                                v_max_index = VIDV_UINT(gvl);
                                v_max_index = VADDVX_UINT(v_max_index, j, gvl);

                                mask = VMFGEVF_FLOAT(v_max, cur_maxf, gvl);
                                max_index = VMFIRSTM(mask,gvl);
#ifdef DOUBLE
                                max_index = *((BLASLONG*)&v_max_index+max_index);
#else
                                max_index = *((unsigned int*)&v_max_index+max_index);
#endif
                        }
                }
        }else{
                gvl = VSETVL(n);
                unsigned int stride_x = inc_x * sizeof(FLOAT);
                unsigned int idx = 0, inc_v = gvl * inc_x;

                v_max_index = VMVVX_UINT(0, gvl);
                v_max = VFMVVF_FLOAT(-1, gvl);
                for(i=0,j=0; i < n/gvl; i++){
                        vx = VLSEV_FLOAT(&x[idx], stride_x, gvl);
                        //fabs(vector)
                        mask = VMFLTVF_FLOAT(vx, 0, gvl);
                        vx = VFRSUBVF_MASK_FLOAT(mask, vx, vx, 0, gvl);

                        //index where element greater than v_max
                        mask = VMFLTVV_FLOAT(v_max, vx, gvl);
                        v_max_index = VIDV_MASK_UINT(mask, v_max_index, gvl);
                        v_max_index = VADDVX_MASK_UINT(mask, v_max_index, v_max_index, j, gvl);

                        //update v_max and start_index j
                        v_max = VFMAXVV_FLOAT(v_max, vx, gvl);
                        j += gvl;
                        idx += inc_v;
                }
                v_res = VFREDMAXVS_FLOAT(v_res, v_max, v_z0, gvl);
                maxf = VFMVFS_FLOAT(v_res);
                mask = VMFGEVF_FLOAT(v_max, maxf, gvl);
                max_index = VMFIRSTM(mask,gvl);
#ifdef DOUBLE
                max_index = *((BLASLONG*)&v_max_index+max_index);
#else
                max_index = *((unsigned int*)&v_max_index+max_index);
#endif
                if(j < n){
                        gvl = VSETVL(n-j);
                        vx = VLSEV_FLOAT(&x[idx], stride_x, gvl);
                        //fabs(vector)
                        mask = VMFLTVF_FLOAT(vx, 0, gvl);
                        v_max = VFRSUBVF_MASK_FLOAT(mask, vx, vx, 0, gvl);

                        v_res = VFREDMAXVS_FLOAT(v_res, v_max, v_z0, gvl);
                        FLOAT cur_maxf = VFMVFS_FLOAT(v_res);
                        if(cur_maxf > maxf){
                                //tail index
                                v_max_index = VIDV_UINT(gvl);
                                v_max_index = VADDVX_UINT(v_max_index, j, gvl);

                                mask = VMFGEVF_FLOAT(v_max, cur_maxf, gvl);
                                max_index = VMFIRSTM(mask,gvl);
#ifdef DOUBLE
                                max_index = *((BLASLONG*)&v_max_index+max_index);
#else
                                max_index = *((unsigned int*)&v_max_index+max_index);
#endif
                        }
                }
        }
	return(max_index+1);
}
