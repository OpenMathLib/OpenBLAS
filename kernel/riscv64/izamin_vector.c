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

#define RVV_EFLOAT RVV_E64
#define FLOAT_V_T float64xm8_t
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
#define VFADDVV_FLOAT vfaddvv_float64xm8
#define VMVVX_UINT vmvvx_uint64xm8
#else

#define ABS fabsf
#define RVV_EFLOAT RVV_E32
#define FLOAT_V_T float32xm8_t
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
#define VFADDVV_FLOAT vfaddvv_float32xm8
#define VMVVX_UINT vmvvx_uint32xm8
#endif

#define RVV_M RVV_M8

BLASLONG CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
	BLASLONG i=0, j=0;
	FLOAT minf=FLT_MAX;
        unsigned int min_index = 0;
	if (n <= 0 || inc_x <= 0) return(min_index);

        FLOAT_V_T vx0, vx1, v_min;
        UINT_V_T v_min_index;
        MASK_T mask0, mask1;
        unsigned int gvl = 0;
        gvl = vsetvli(n, RVV_EFLOAT, RVV_M);
        v_min_index = VMVVX_UINT(0, gvl);
        v_min = VFMVVF_FLOAT(FLT_MAX, gvl);
        BLASLONG stride_x = inc_x * 2 * sizeof(FLOAT);
        BLASLONG inc_xv = gvl * inc_x * 2;
        BLASLONG ix = 0;
        for(i=0,j=0; i < n/gvl; i++){
                vx0 = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                //fabs(vector)
                mask0 = VMFLTVF_FLOAT(vx0, 0, gvl);
                vx0 = VFRSUBVF_MASK_FLOAT(vx0, vx0, 0, mask0, gvl);
/*
#if defined(DOUBLE)
asm volatile(
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e64,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+v"(vx0)
        :"v"(mask0), "f"(zero), "r"(gvl)
        :"v0");
#else
asm volatile(
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e32,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+v"(vx0)
        :"v"(mask0), "f"(zero), "r"(gvl)
        :"v0");
#endif
*/
                vx1 = VLSEV_FLOAT(&x[ix+1], stride_x, gvl);
                //fabs(vector)
                mask1 = VMFLTVF_FLOAT(vx1, 0, gvl);
                vx1 = VFRSUBVF_MASK_FLOAT(vx1, vx1, 0, mask1, gvl);
/*
#if defined(DOUBLE)
asm volatile(
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e64,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+v"(vx1)
        :"v"(mask1), "f"(zero), "r"(gvl)
        :"v0");
#else
asm volatile(
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e32,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+v"(vx1)
        :"v"(mask1), "f"(zero), "r"(gvl)
        :"v0");
#endif
*/
                vx0 = VFADDVV_FLOAT(vx0, vx1, gvl);

                //index where element less than v_min
                mask0 = VMFLTVV_FLOAT(vx0, v_min, gvl);
                v_min_index = VIDV_MASK_UINT(v_min_index, mask0, gvl);
/*
#if defined(DOUBLE)
asm volatile(
        "vor.vv v0, %1, %1 \n\t"
        "vsetvli x0, %2, e64,m8 \n\t"
        "vid.v %0, v0.t \n\t"
        :"+v"(v_min_index)
        :"v"(mask0), "r"(gvl)
        :"v0");
#else
asm volatile(
        "vor.vv v0, %1, %1 \n\t"
        "vsetvli x0, %2, e32,m8 \n\t"
        "vid.v %0, v0.t \n\t"
        :"+v"(v_min_index)
        :"v"(mask0), "r"(gvl)
        :"v0");
#endif
*/
                v_min_index = VADDVX_MASK_UINT(v_min_index, v_min_index, j, mask0, gvl);

                //update v_min and start_index j
                v_min = VFMINVV_FLOAT(v_min, vx0, gvl);
                j += gvl;
                ix += inc_xv;
        }
        vx0 = VFMVVF_FLOAT(FLT_MAX, gvl);
        vx0 = VFREDMINVS_FLOAT(v_min, vx0, gvl);
        minf = vx0[0];
        mask0 = VMFLEVF_FLOAT(v_min, minf, gvl);
        min_index = VMFIRSTM(mask0,gvl);
        min_index = v_min_index[min_index];

        if(j < n){
                gvl = vsetvli(n-j, RVV_EFLOAT, RVV_M);
                v_min_index = VMVVX_UINT(0, gvl);
                vx0 = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                //fabs(vector)
                mask0 = VMFLTVF_FLOAT(vx0, 0, gvl);
                vx0 = VFRSUBVF_MASK_FLOAT(vx0, vx0, 0, mask0, gvl);
/*
#if defined(DOUBLE)
asm volatile(
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e64,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+v"(vx0)
        :"v"(mask0), "f"(zero), "r"(gvl)
        :"v0");
#else
asm volatile(
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e32,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+v"(vx0)
        :"v"(mask0), "f"(zero), "r"(gvl)
        :"v0");
#endif
*/
                vx1 = VLSEV_FLOAT(&x[ix+1], stride_x, gvl);
                //fabs(vector)
                mask1 = VMFLTVF_FLOAT(vx1, 0, gvl);
                vx1 = VFRSUBVF_MASK_FLOAT(vx1, vx1, 0, mask1, gvl);
/*
#if defined(DOUBLE)
asm volatile(
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e64,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+v"(vx1)
        :"v"(mask1), "f"(zero), "r"(gvl)
        :"v0");
#else
asm volatile(
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e32,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+v"(vx1)
        :"v"(mask1), "f"(zero), "r"(gvl)
        :"v0");
#endif
*/
                v_min = VFADDVV_FLOAT(vx0, vx1, gvl);
                vx0 = VFMVVF_FLOAT(FLT_MAX, gvl);
                vx0 = VFREDMINVS_FLOAT(v_min, vx0, gvl);
                FLOAT cur_minf = vx0[0];
                if(cur_minf < minf){
                        //tail index
                        v_min_index = VIDV_UINT(gvl);
                        v_min_index = VADDVX_UINT(v_min_index, j, gvl);

                        mask0 = VMFLEVF_FLOAT(v_min, cur_minf, gvl);
                        min_index = VMFIRSTM(mask0,gvl);
                        min_index = v_min_index[min_index];
                }
        }
	return(min_index+1);
}


