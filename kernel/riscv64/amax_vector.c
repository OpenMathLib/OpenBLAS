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

#if !defined(DOUBLE)
#define RVV_EFLOAT RVV_E32
#define RVV_M RVV_M8
#define FLOAT_V_T float32xm8_t
#define VLEV_FLOAT vlev_float32xm8
#define VLSEV_FLOAT vlsev_float32xm8
#define VFREDMAXVS_FLOAT vfredmaxvs_float32xm8
#define MASK_T e32xm8_t
#define VMFLTVF_FLOAT vmfltvf_e32xm8_float32xm8
#define VFMVVF_FLOAT vfmvvf_float32xm8
#define VFRSUBVF_MASK_FLOAT vfrsubvf_mask_float32xm8
#define VFMAXVV_FLOAT vfmaxvv_float32xm8
#else
#define RVV_EFLOAT RVV_E64
#define RVV_M RVV_M8
#define FLOAT_V_T float64xm8_t
#define VLEV_FLOAT vlev_float64xm8
#define VLSEV_FLOAT vlsev_float64xm8
#define VFREDMAXVS_FLOAT vfredmaxvs_float64xm8
#define MASK_T e64xm8_t
#define VMFLTVF_FLOAT vmfltvf_e64xm8_float64xm8
#define VFMVVF_FLOAT vfmvvf_float64xm8
#define VFRSUBVF_MASK_FLOAT vfrsubvf_mask_float64xm8
#define VFMAXVV_FLOAT vfmaxvv_float64xm8
#endif

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
	BLASLONG i=0, j=0;
	BLASLONG ix=0;
	FLOAT maxf=0.0;
	if (n <= 0 || inc_x <= 0) return(maxf);
        unsigned int gvl = 0;
        FLOAT_V_T v0, v1, v_max;

        MASK_T mask0, mask1;
        FLOAT zero = 0.0;
        if(inc_x == 1){
                gvl = vsetvli(n, RVV_EFLOAT, RVV_M);
                if(gvl <= n/2){
                        v_max = VFMVVF_FLOAT(0, gvl);
                        for(i=0,j=0; i<n/(gvl*2); i++){
                                v0 = VLEV_FLOAT(&x[j], gvl);
                                mask0 = VMFLTVF_FLOAT(v0, 0, gvl);
                                //v0 = VFRSUBVF_MASK_FLOAT(v0, 0, mask0, gvl);
#if defined(DOUBLE)
asm volatile(
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e64,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+v"(v0)
        :"v"(mask0), "f"(zero), "r"(gvl)
        :"v0");
#else
asm volatile(
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e32,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+v"(v0)
        :"v"(mask0), "f"(zero), "r"(gvl)
        :"v0");
#endif

                                v_max = VFMAXVV_FLOAT(v_max, v0, gvl);

                                v1 = VLEV_FLOAT(&x[j+gvl], gvl);
                                mask1 = VMFLTVF_FLOAT(v1, 0, gvl);
                                //v1 = VFRSUBVF_MASK_FLOAT(v1, 0, mask1, gvl);
#if defined(DOUBLE)
asm volatile(
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e64,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+v"(v1)
        :"v"(mask1), "f"(zero), "r"(gvl)
        :"v0");
#else
asm volatile(
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e32,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+v"(v1)
        :"v"(mask1), "f"(zero), "r"(gvl)
        :"v0");
#endif

                                v_max = VFMAXVV_FLOAT(v_max, v1, gvl);
                                j += gvl*2;
                        }
                        v0 = VFMVVF_FLOAT(0, gvl);
                        v0 = VFREDMAXVS_FLOAT(v_max, v0, gvl);
                        maxf = v0[0];
                }
                for(;j<n;){
                        gvl = vsetvli(n-j, RVV_EFLOAT, RVV_M);
                        v0 = VLEV_FLOAT(&x[j], gvl);
                        mask0 = VMFLTVF_FLOAT(v0, 0, gvl);
                        //v0 = VFRSUBVF_MASK_FLOAT(v0, 0, mask0, gvl);
#if defined(DOUBLE)
asm volatile(
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e64,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+v"(v0)
        :"v"(mask0), "f"(zero), "r"(gvl)
        :"v0");
#else
asm volatile(
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e32,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+v"(v0)
        :"v"(mask0), "f"(zero), "r"(gvl)
        :"v0");
#endif

                        v1 = VFMVVF_FLOAT(0, gvl);
                        v0 = VFREDMAXVS_FLOAT(v0, v1, gvl);
                        if(v0[0] > maxf)
                                maxf = v0[0];
                        j += gvl;
                }
        }else{
                gvl = vsetvli(n, RVV_EFLOAT, RVV_M);
                BLASLONG stride_x = inc_x * sizeof(FLOAT);
                if(gvl <= n/2){
                        BLASLONG inc_xv = inc_x * gvl;
                        v_max = VFMVVF_FLOAT(0, gvl);
                        for(i=0,j=0; i<n/(gvl*2); i++){
                                v0 = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                                mask0 = VMFLTVF_FLOAT(v0, 0, gvl);
                                //v0 = VFRSUBVF_MASK_FLOAT(v0, 0, mask0, gvl);
#if defined(DOUBLE)
asm volatile(
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e64,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+v"(v0)
        :"v"(mask0), "f"(zero), "r"(gvl)
        :"v0");
#else
asm volatile(
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e32,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+v"(v0)
        :"v"(mask0), "f"(zero), "r"(gvl)
        :"v0");
#endif

                                v_max = VFMAXVV_FLOAT(v_max, v0, gvl);

                                v1 = VLSEV_FLOAT(&x[ix+inc_xv], stride_x, gvl);
                                mask1 = VMFLTVF_FLOAT(v1, 0, gvl);
                                //v1 = VFRSUBVF_MASK_FLOAT(v1, 0, mask1, gvl);
#if defined(DOUBLE)
asm volatile(
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e64,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+v"(v1)
        :"v"(mask1), "f"(zero), "r"(gvl)
        :"v0");
#else
asm volatile(
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e32,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+v"(v1)
        :"v"(mask1), "f"(zero), "r"(gvl)
        :"v0");
#endif

                                v_max = VFMAXVV_FLOAT(v_max, v1, gvl);
                                j += gvl*2;
                                ix += inc_xv*2;
                        }
                        v0 = VFMVVF_FLOAT(0, gvl);
                        v0 = VFREDMAXVS_FLOAT(v_max, v0, gvl);
                        maxf = v0[0];
                }
                for(;j<n;){
                        gvl = vsetvli(n-j, RVV_EFLOAT, RVV_M);
                        v0 = VLSEV_FLOAT(&x[j*inc_x], stride_x, gvl);
                        mask0 = VMFLTVF_FLOAT(v0, 0, gvl);
                        //v0 = VFRSUBVF_MASK_FLOAT(v0, 0, mask0, gvl);
#if defined(DOUBLE)
asm volatile(
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e64,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+v"(v0)
        :"v"(mask0), "f"(zero), "r"(gvl)
        :"v0");
#else
asm volatile(
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e32,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+v"(v0)
        :"v"(mask0), "f"(zero), "r"(gvl)
        :"v0");
#endif

                        v1 = VFMVVF_FLOAT(0, gvl);
                        v0 = VFREDMAXVS_FLOAT(v0, v1, gvl);
                        if(v0[0] > maxf)
                                maxf = v0[0];
                        j += gvl;
                }
        }
	return(maxf);
}


