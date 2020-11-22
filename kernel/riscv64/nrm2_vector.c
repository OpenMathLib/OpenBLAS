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
#define RVV_EFLOAT RVV_E32
#define RVV_M RVV_M4
#define FLOAT_V_T float32xm4_t
#define VLEV_FLOAT vlev_float32xm4
#define VLSEV_FLOAT vlsev_float32xm4
#define VFREDSUM_FLOAT vfredsumvs_float32xm4
#define VFMACCVV_FLOAT vfmaccvv_float32xm4
#define VFMVVF_FLOAT vfmvvf_float32xm4
#define VFDOTVV_FLOAT vfdotvv_float32xm4
#define ABS fabsf
#define MASK_T e32xm4_t
#define VFRSUBVF_MASK_FLOAT vfrsubvf_mask_float32xm4
#define VMFGTVF_FLOAT vmfgtvf_e32xm4_float32xm4
#define VMFIRSTM vmfirstm_e32xm4
#define VFDIVVF_FLOAT vfdivvf_float32xm4
#define VMFLTVF_FLOAT vmfltvf_e32xm4_float32xm4
#define VFREDMAXVS_FLOAT vfredmaxvs_float32xm4
#else
#define RVV_EFLOAT RVV_E64
#define RVV_M RVV_M4
#define FLOAT_V_T float64xm4_t
#define VLEV_FLOAT vlev_float64xm4
#define VLSEV_FLOAT vlsev_float64xm4
#define VFREDSUM_FLOAT vfredsumvs_float64xm4
#define VFMACCVV_FLOAT vfmaccvv_float64xm4
#define VFMVVF_FLOAT vfmvvf_float64xm4
#define VFDOTVV_FLOAT vfdotvv_float64xm4
#define ABS fabs
#define MASK_T e64xm4_t
#define VFRSUBVF_MASK_FLOAT vfrsubvf_mask_float64xm4
#define VMFGTVF_FLOAT vmfgtvf_e64xm4_float64xm4
#define VMFIRSTM vmfirstm_e64xm4
#define VFDIVVF_FLOAT vfdivvf_float64xm4
#define VMFLTVF_FLOAT vmfltvf_e64xm4_float64xm4
#define VFREDMAXVS_FLOAT vfredmaxvs_float64xm4
#endif

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
	BLASLONG i=0, j=0;

	if ( n < 0 )  return(0.0);
        if(n == 1) return (ABS(x[0]));

        FLOAT_V_T vr, v0, v_zero;
        unsigned int gvl = 0;
        FLOAT scale = 0.0, ssq = 0.0;
        MASK_T mask;
        BLASLONG index = 0;
        if(inc_x == 1){
                gvl = vsetvli(n, RVV_EFLOAT, RVV_M);
                vr = VFMVVF_FLOAT(0, gvl);
                v_zero = VFMVVF_FLOAT(0, gvl);
                for(i=0,j=0; i<n/gvl; i++){
                        v0 = VLEV_FLOAT(&x[j], gvl);
                        //fabs(vector)
                        mask = VMFLTVF_FLOAT(v0, 0, gvl);
                        v0 = VFRSUBVF_MASK_FLOAT(v0, v0, 0, mask, gvl);
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
                                vr = VFREDSUM_FLOAT(vr, v_zero, gvl);
                                //total ssq before current vector
                                ssq += vr[0];
                                //find max
                                vr = VFREDMAXVS_FLOAT(v0, v_zero, gvl);
                                //update ssq before max_index
                                ssq = ssq * (scale/vr[0])*(scale/vr[0]);
                                //update scale
                                scale = vr[0];
                                //ssq in vector vr
                                v0 = VFDIVVF_FLOAT(v0, scale, gvl);
                                vr = VFMACCVV_FLOAT(v_zero, v0, v0, gvl);
                        }
                        j += gvl;
                }
                //ssq in vector vr: vr[0]
                vr = VFREDSUM_FLOAT(vr, v_zero, gvl);
                //total ssq now
                ssq += vr[0];

                //tail
                if(j < n){
                        gvl = vsetvli(n-j, RVV_EFLOAT, RVV_M);
                        v0 = VLEV_FLOAT(&x[j], gvl);
                        //fabs(vector)
                        mask = VMFLTVF_FLOAT(v0, 0, gvl);
                        v0 = VFRSUBVF_MASK_FLOAT(v0, v0, 0, mask, gvl);
                        //if scale change
                        mask = VMFGTVF_FLOAT(v0, scale, gvl);
                        index = VMFIRSTM(mask, gvl);
                        if(index == -1){//no elements greater than scale
                                if(scale != 0.0)
                                        v0 = VFDIVVF_FLOAT(v0, scale, gvl);
                        }else{//found greater element
                                //find max
                                vr = VFREDMAXVS_FLOAT(v0, v_zero, gvl);
                                //update ssq before max_index
                                ssq = ssq * (scale/vr[0])*(scale/vr[0]);
                                //update scale
                                scale = vr[0];
                                v0 = VFDIVVF_FLOAT(v0, scale, gvl);
                        }
                        vr = VFMACCVV_FLOAT(v_zero, v0, v0, gvl);
                        //ssq in vector vr: vr[0]
                        vr = VFREDSUM_FLOAT(vr, v_zero, gvl);
                        //total ssq now
                        ssq += vr[0];
                }
        }else{
                gvl = vsetvli(n, RVV_EFLOAT, RVV_M);
                vr = VFMVVF_FLOAT(0, gvl);
                v_zero = VFMVVF_FLOAT(0, gvl);
                unsigned int stride_x = inc_x * sizeof(FLOAT);
                int idx = 0, inc_v = inc_x * gvl;
                for(i=0,j=0; i<n/gvl; i++){
                        v0 = VLSEV_FLOAT(&x[idx], stride_x, gvl);
                        //fabs(vector)
                        mask = VMFLTVF_FLOAT(v0, 0, gvl);
                        v0 = VFRSUBVF_MASK_FLOAT(v0, v0, 0, mask, gvl);
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
                                vr = VFREDSUM_FLOAT(vr, v_zero, gvl);
                                //total ssq before current vector
                                ssq += vr[0];
                                //find max
                                vr = VFREDMAXVS_FLOAT(v0, v_zero, gvl);
                                //update ssq before max_index
                                ssq = ssq * (scale/vr[0])*(scale/vr[0]);
                                //update scale
                                scale = vr[0];
                                //ssq in vector vr
                                v0 = VFDIVVF_FLOAT(v0, scale, gvl);
                                vr = VFMACCVV_FLOAT(v_zero, v0, v0, gvl);
                        }
                        j += gvl;
                        idx += inc_v;
                }
                //ssq in vector vr: vr[0]
                vr = VFREDSUM_FLOAT(vr, v_zero, gvl);
                //total ssq now
                ssq += vr[0];

                //tail
                if(j < n){
                        gvl = vsetvli(n-j, RVV_EFLOAT, RVV_M);
                        v0 = VLSEV_FLOAT(&x[idx], stride_x, gvl);
                        //fabs(vector)
                        mask = VMFLTVF_FLOAT(v0, 0, gvl);
                        v0 = VFRSUBVF_MASK_FLOAT(v0, v0, 0, mask, gvl);
                        //if scale change
                        mask = VMFGTVF_FLOAT(v0, scale, gvl);
                        index = VMFIRSTM(mask, gvl);
                        if(index == -1){//no elements greater than scale
                                if(scale != 0.0)
                                        v0 = VFDIVVF_FLOAT(v0, scale, gvl);
                        }else{//found greater element
                                //find max
                                vr = VFREDMAXVS_FLOAT(v0, v_zero, gvl);
                                //update ssq before max_index
                                ssq = ssq * (scale/vr[0])*(scale/vr[0]);
                                //update scale
                                scale = vr[0];
                                v0 = VFDIVVF_FLOAT(v0, scale, gvl);
                        }
                        vr = VFMACCVV_FLOAT(v_zero, v0, v0, gvl);
                        //ssq in vector vr: vr[0]
                        vr = VFREDSUM_FLOAT(vr, v_zero, gvl);
                        //total ssq now
                        ssq += vr[0];
                }
        }
	return(scale * sqrt(ssq));
}


