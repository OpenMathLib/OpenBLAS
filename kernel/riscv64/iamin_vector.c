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

#if defined(DOUBLE)

#define ABS fabs
#define RVV_EFLOAT RVV_E64
#define RVV_M RVV_M8
#define FLOAT_V_T float64xm8_t
#define VLEV_FLOAT vlev_float64xm8
#define VLSEV_FLOAT vlsev_float64xm8
#define VFREDMINVS_FLOAT vfredminvs_float64xm8
#define MASK_T e64xm8_t
#define VMFLTVF_FLOAT vmfltvf_e64xm8_float64xm8
#define VMFLTVV_FLOAT vmfltvv_e64xm8_float64xm8
#define VFMVVF_FLOAT vfmvvf_float64xm8
#define VFRSUBVF_MASK_FLOAT vfrsubvf_mask_float64xm8
#define VFMINVV_FLOAT vfminvv_float64xm8
#define VMFLEVF_FLOAT vmflevf_e64xm8_float64xm8
#define VMFIRSTM vmfirstm_e64xm8
#define UINT_V_T uint64xm8_t
#define VIDV_MASK_UINT vidv_mask_uint64xm8
#define VIDV_UINT vidv_uint64xm8
#define VADDVX_MASK_UINT vaddvx_mask_uint64xm8
#define VADDVX_UINT vaddvx_uint64xm8
#define VMVVX_UINT vmvvx_uint64xm8
#else

#define ABS fabsf
#define RVV_EFLOAT RVV_E32
#define RVV_M RVV_M8
#define FLOAT_V_T float32xm8_t
#define VLEV_FLOAT vlev_float32xm8
#define VLSEV_FLOAT vlsev_float32xm8
#define VFREDMINVS_FLOAT vfredminvs_float32xm8
#define MASK_T e32xm8_t
#define VMFLTVF_FLOAT vmfltvf_e32xm8_float32xm8
#define VMFLTVV_FLOAT vmfltvv_e32xm8_float32xm8
#define VFMVVF_FLOAT vfmvvf_float32xm8
#define VFRSUBVF_MASK_FLOAT vfrsubvf_mask_float32xm8
#define VFMINVV_FLOAT vfminvv_float32xm8
#define VMFLEVF_FLOAT vmflevf_e32xm8_float32xm8
#define VMFIRSTM vmfirstm_e32xm8
#define UINT_V_T uint32xm8_t
#define VIDV_MASK_UINT vidv_mask_uint32xm8
#define VIDV_UINT vidv_uint32xm8
#define VADDVX_MASK_UINT vaddvx_mask_uint32xm8
#define VADDVX_UINT vaddvx_uint32xm8
#define VMVVX_UINT vmvvx_uint32xm8
#endif


BLASLONG CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
	BLASLONG i=0, j=0;
	FLOAT minf=FLT_MAX;
        unsigned int min_index = 0;
	if (n <= 0 || inc_x <= 0) return(min_index);

        FLOAT_V_T vx, v_min;
        UINT_V_T v_min_index;
        MASK_T mask;
        unsigned int gvl = 0;
        if(inc_x == 1){
                gvl = vsetvli(n, RVV_EFLOAT, RVV_M);
                v_min = VFMVVF_FLOAT(FLT_MAX, gvl);
                v_min_index = VMVVX_UINT(0, gvl);
                for(i=0,j=0; i < n/gvl; i++){
                        vx = VLEV_FLOAT(&x[j], gvl);
                        //fabs(vector)
                        mask = VMFLTVF_FLOAT(vx, 0, gvl);
                        vx = VFRSUBVF_MASK_FLOAT(vx, vx, 0, mask, gvl);

                        //index where element less than v_min
                        mask = VMFLTVV_FLOAT(vx, v_min, gvl);
                        v_min_index = VIDV_MASK_UINT(v_min_index, mask, gvl);
                        v_min_index = VADDVX_MASK_UINT(v_min_index, v_min_index, j, mask, gvl);

                        //update v_min and start_index j
                        v_min = VFMINVV_FLOAT(v_min, vx, gvl);
                        j += gvl;
                }
                vx = VFMVVF_FLOAT(FLT_MAX, gvl);
                vx = VFREDMINVS_FLOAT(v_min, vx, gvl);
                minf = vx[0];
                mask = VMFLEVF_FLOAT(v_min, minf, gvl);
                min_index = VMFIRSTM(mask,gvl);
                min_index = v_min_index[min_index];

                if(j < n){
                        gvl = vsetvli(n-j, RVV_EFLOAT, RVV_M);
                        vx = VLEV_FLOAT(&x[j], gvl);
                        //fabs(vector)
                        mask = VMFLTVF_FLOAT(vx, 0, gvl);
                        v_min = VFRSUBVF_MASK_FLOAT(vx, vx, 0, mask, gvl);

                        vx = VFMVVF_FLOAT(FLT_MAX, gvl);
                        vx = VFREDMINVS_FLOAT(v_min, vx, gvl);
                        FLOAT cur_minf = vx[0];
                        if(cur_minf < minf){
                                //tail index
                                v_min_index = VIDV_UINT(gvl);
                                v_min_index = VADDVX_UINT(v_min_index, j, gvl);

                                mask = VMFLEVF_FLOAT(v_min, cur_minf, gvl);
                                min_index = VMFIRSTM(mask,gvl);
                                min_index = v_min_index[min_index];
                        }
                }
        }else{
                gvl = vsetvli(n, RVV_EFLOAT, RVV_M);
                unsigned int stride_x = inc_x * sizeof(FLOAT);
                unsigned int idx = 0, inc_v = gvl * inc_x;

                v_min = VFMVVF_FLOAT(FLT_MAX, gvl);
                v_min_index = VMVVX_UINT(0, gvl);
                for(i=0,j=0; i < n/gvl; i++){
                        vx = VLSEV_FLOAT(&x[idx], stride_x, gvl);
                        //fabs(vector)
                        mask = VMFLTVF_FLOAT(vx, 0, gvl);
                        vx = VFRSUBVF_MASK_FLOAT(vx, vx, 0, mask, gvl);

                        //index where element less than v_min
                        mask = VMFLTVV_FLOAT(vx, v_min, gvl);
                        v_min_index = VIDV_MASK_UINT(v_min_index, mask, gvl);
                        v_min_index = VADDVX_MASK_UINT(v_min_index, v_min_index, j, mask, gvl);

                        //update v_min and start_index j
                        v_min = VFMINVV_FLOAT(v_min, vx, gvl);
                        j += gvl;
                        idx += inc_v;
                }
                vx = VFMVVF_FLOAT(FLT_MAX, gvl);
                vx = VFREDMINVS_FLOAT(v_min, vx, gvl);
                minf = vx[0];
                mask = VMFLEVF_FLOAT(v_min, minf, gvl);
                min_index = VMFIRSTM(mask,gvl);
                min_index = v_min_index[min_index];

                if(j < n){
                        gvl = vsetvli(n-j, RVV_EFLOAT, RVV_M);
                        vx = VLSEV_FLOAT(&x[idx], stride_x, gvl);
                        //fabs(vector)
                        mask = VMFLTVF_FLOAT(vx, 0, gvl);
                        v_min = VFRSUBVF_MASK_FLOAT(vx, vx, 0, mask, gvl);

                        vx = VFMVVF_FLOAT(FLT_MAX, gvl);
                        vx = VFREDMINVS_FLOAT(v_min, vx, gvl);
                        FLOAT cur_minf = vx[0];
                        if(cur_minf < minf){
                                //tail index
                                v_min_index = VIDV_UINT(gvl);
                                v_min_index = VADDVX_UINT(v_min_index, j, gvl);

                                mask = VMFLEVF_FLOAT(v_min, cur_minf, gvl);
                                min_index = VMFIRSTM(mask,gvl);
                                min_index = v_min_index[min_index];
                        }
                }
        }
	return(min_index+1);
}


