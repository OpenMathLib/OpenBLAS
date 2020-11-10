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

#define ABS fabs
#define RVV_EFLOAT RVV_E64
#define RVV_M RVV_M8
#define FLOAT_V_T float64xm8_t
#define VLEV_FLOAT vlev_float64xm8
#define VLSEV_FLOAT vlsev_float64xm8
#define VFREDMAXVS_FLOAT vfredmaxvs_float64xm8
#define MASK_T e64xm8_t
#define VMFLTVF_FLOAT vmfltvf_e64xm8_float64xm8
#define VMFLTVV_FLOAT vmfltvv_e64xm8_float64xm8
#define VFMVVF_FLOAT vfmvvf_float64xm8
#define VFRSUBVF_MASK_FLOAT vfrsubvf_mask_float64xm8
#define VFMAXVV_FLOAT vfmaxvv_float64xm8
#define VMFGEVF_FLOAT vmfgevf_e64xm8_float64xm8
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
#define VFREDMAXVS_FLOAT vfredmaxvs_float32xm8
#define MASK_T e32xm8_t
#define VMFLTVF_FLOAT vmfltvf_e32xm8_float32xm8
#define VMFLTVV_FLOAT vmfltvv_e32xm8_float32xm8
#define VFMVVF_FLOAT vfmvvf_float32xm8
#define VFRSUBVF_MASK_FLOAT vfrsubvf_mask_float32xm8
#define VFMAXVV_FLOAT vfmaxvv_float32xm8
#define VMFGEVF_FLOAT vmfgevf_e32xm8_float32xm8
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
	FLOAT maxf=0.0;
        unsigned int max_index = 0;
	if (n <= 0 || inc_x <= 0) return(max_index);

        FLOAT_V_T vx, v_max;
        UINT_V_T v_max_index;
        MASK_T mask;
        unsigned int gvl = 0;
        if(inc_x == 1){
                gvl = vsetvli(n, RVV_EFLOAT, RVV_M);
                v_max_index = VMVVX_UINT(0, gvl);
                v_max = VFMVVF_FLOAT(-1, gvl);
                for(i=0,j=0; i < n/gvl; i++){
                        vx = VLEV_FLOAT(&x[j], gvl);
                        //fabs(vector)
                        mask = VMFLTVF_FLOAT(vx, 0, gvl);
                        vx = VFRSUBVF_MASK_FLOAT(vx, vx, 0, mask, gvl);

                        //index where element greater than v_max
                        mask = VMFLTVV_FLOAT(v_max, vx, gvl);
                        v_max_index = VIDV_MASK_UINT(v_max_index, mask, gvl);
                        v_max_index = VADDVX_MASK_UINT(v_max_index, v_max_index, j, mask, gvl);

                        //update v_max and start_index j
                        v_max = VFMAXVV_FLOAT(v_max, vx, gvl);
                        j += gvl;
                }
                vx = VFMVVF_FLOAT(0, gvl);
                vx = VFREDMAXVS_FLOAT(v_max, vx, gvl);
                maxf = vx[0];
                mask = VMFGEVF_FLOAT(v_max, maxf, gvl);
                max_index = VMFIRSTM(mask,gvl);
                max_index = v_max_index[max_index];

                if(j < n){
                        gvl = vsetvli(n-j, RVV_EFLOAT, RVV_M);
                        vx = VLEV_FLOAT(&x[j], gvl);
                        //fabs(vector)
                        mask = VMFLTVF_FLOAT(vx, 0, gvl);
                        v_max = VFRSUBVF_MASK_FLOAT(vx, vx, 0, mask, gvl);

                        vx = VFMVVF_FLOAT(0, gvl);
                        vx = VFREDMAXVS_FLOAT(v_max, vx, gvl);
                        FLOAT cur_maxf = vx[0];
                        if(cur_maxf > maxf){
                                //tail index
                                v_max_index = VIDV_UINT(gvl);
                                v_max_index = VADDVX_UINT(v_max_index, j, gvl);

                                mask = VMFGEVF_FLOAT(v_max, cur_maxf, gvl);
                                max_index = VMFIRSTM(mask,gvl);
                                max_index = v_max_index[max_index];
                        }
                }
        }else{
                gvl = vsetvli(n, RVV_EFLOAT, RVV_M);
                unsigned int stride_x = inc_x * sizeof(FLOAT);
                unsigned int idx = 0, inc_v = gvl * inc_x;

                v_max_index = VMVVX_UINT(0, gvl);
                v_max = VFMVVF_FLOAT(-1, gvl);
                for(i=0,j=0; i < n/gvl; i++){
                        vx = VLSEV_FLOAT(&x[idx], stride_x, gvl);
                        //fabs(vector)
                        mask = VMFLTVF_FLOAT(vx, 0, gvl);
                        vx = VFRSUBVF_MASK_FLOAT(vx, vx, 0, mask, gvl);

                        //index where element greater than v_max
                        mask = VMFLTVV_FLOAT(v_max, vx, gvl);
                        v_max_index = VIDV_MASK_UINT(v_max_index, mask, gvl);
                        v_max_index = VADDVX_MASK_UINT(v_max_index, v_max_index, j, mask, gvl);

                        //update v_max and start_index j
                        v_max = VFMAXVV_FLOAT(v_max, vx, gvl);
                        j += gvl;
                        idx += inc_v;
                }
                vx = VFMVVF_FLOAT(0, gvl);
                vx = VFREDMAXVS_FLOAT(v_max, vx, gvl);
                maxf = vx[0];
                mask = VMFGEVF_FLOAT(v_max, maxf, gvl);
                max_index = VMFIRSTM(mask,gvl);
                max_index = v_max_index[max_index];

                if(j < n){
                        gvl = vsetvli(n-j, RVV_EFLOAT, RVV_M);
                        vx = VLSEV_FLOAT(&x[idx], stride_x, gvl);
                        //fabs(vector)
                        mask = VMFLTVF_FLOAT(vx, 0, gvl);
                        v_max = VFRSUBVF_MASK_FLOAT(vx, vx, 0, mask, gvl);

                        vx = VFMVVF_FLOAT(0, gvl);
                        vx = VFREDMAXVS_FLOAT(v_max, vx, gvl);
                        FLOAT cur_maxf = vx[0];
                        if(cur_maxf > maxf){
                                //tail index
                                v_max_index = VIDV_UINT(gvl);
                                v_max_index = VADDVX_UINT(v_max_index, j, gvl);

                                mask = VMFGEVF_FLOAT(v_max, cur_maxf, gvl);
                                max_index = VMFIRSTM(mask,gvl);
                                max_index = v_max_index[max_index];
                        }
                }
        }
	return(max_index+1);
}


