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
#define VLEV_FLOAT vle32_v_f32m8
#define VLSEV_FLOAT vlse32_v_f32m8
#define VFREDMINVS_FLOAT vfredmin_vs_f32m8_f32m1
#define MASK_T vbool4_t
#define VMFLTVF_FLOAT vmflt_vf_f32m8_b4
#define VFMVVF_FLOAT vfmv_v_f_f32m8
#define VFMVVF_FLOAT_M1 vfmv_v_f_f32m1
#define VFRSUBVF_MASK_FLOAT vfrsub_vf_f32m8_m
#define VFMINVV_FLOAT vfmin_vv_f32m8
#else
#define VSETVL(n) vsetvl_e64m8(n)
#define VSETVL_MAX vsetvlmax_e32m1()
#define FLOAT_V_T vfloat64m8_t
#define FLOAT_V_T_M1 vfloat64m1_t
#define VLEV_FLOAT vle64_v_f64m8
#define VLSEV_FLOAT vlse64_v_f64m8
#define VFREDMINVS_FLOAT vfredmin_vs_f64m8_f64m1
#define MASK_T vbool8_t
#define VMFLTVF_FLOAT vmflt_vf_f64m8_b8
#define VFMVVF_FLOAT vfmv_v_f_f64m8
#define VFMVVF_FLOAT_M1 vfmv_v_f_f64m1
#define VFRSUBVF_MASK_FLOAT vfrsub_vf_f64m8_m
#define VFMINVV_FLOAT vfmin_vv_f64m8
#endif

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
	BLASLONG i=0, j=0;
	if (n <= 0 || inc_x <= 0) return(0.0);
	FLOAT minf=FLT_MAX;
        unsigned int gvl = 0;
        FLOAT_V_T v0, v1, v_min;
        FLOAT_V_T_M1 v_res, v_max;
        gvl = VSETVL_MAX;
        v_res = VFMVVF_FLOAT_M1(0, gvl);
        v_max = VFMVVF_FLOAT_M1(FLT_MAX, gvl);

        MASK_T mask0, mask1;
	    FLOAT zero = 0.0;
        if(inc_x == 1){
                gvl = VSETVL(n);
                if(gvl <= n/2){
                        v_min = VFMVVF_FLOAT(FLT_MAX, gvl);
                        for(i=0,j=0; i<n/(gvl*2); i++){
                                v0 = VLEV_FLOAT(&x[j], gvl);
                                mask0 = VMFLTVF_FLOAT(v0, 0, gvl);
                                //v0 = VFRSUBVF_MASK_FLOAT(v0, 0, mask0, gvl);
#if defined(DOUBLE)
asm volatile(
        "vsetvli    zero, zero, e8, m1\n\t"
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e64,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+vd"(v0)
        :"vd"(mask0), "f"(zero), "r"(gvl)
        :"v0");
#else
asm volatile(
        "vsetvli    zero, zero, e8, m1\n\t"
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e32,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+vd"(v0)
        :"vd"(mask0), "f"(zero), "r"(gvl)
        :"v0");
#endif
                                v_min = VFMINVV_FLOAT(v_min, v0, gvl);

                                v1 = VLEV_FLOAT(&x[j+gvl], gvl);
                                mask1 = VMFLTVF_FLOAT(v1, 0, gvl);
                                //v1 = VFRSUBVF_MASK_FLOAT(v1, 0, mask1, gvl);
#if defined(DOUBLE)
asm volatile(
        "vsetvli    zero, zero, e8, m1\n\t"
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e64,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+vd"(v1)
        :"vd"(mask1), "f"(zero), "r"(gvl)
        :"v0");
#else
asm volatile(
        "vsetvli    zero, zero, e8, m1\n\t"
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e32,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+vd"(v1)
        :"vd"(mask1), "f"(zero), "r"(gvl)
        :"v0");
#endif

                                v_min = VFMINVV_FLOAT(v_min, v1, gvl);
                                j += gvl*2;
                        }
                        v_res = VFREDMINVS_FLOAT(v_res, v_min, v_max, gvl);
                        minf = *((FLOAT*)&v_res);
                }
                for(;j<n;){
                        gvl = VSETVL(n-j);
                        v0 = VLEV_FLOAT(&x[j], gvl);
                        mask0 = VMFLTVF_FLOAT(v0, 0, gvl);
                        //v0 = VFRSUBVF_MASK_FLOAT(v0, 0, mask0, gvl);
#if defined(DOUBLE)
asm volatile(
        "vsetvli    zero, zero, e8, m1\n\t"
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e64,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+vd"(v0)
        :"vd"(mask0), "f"(zero), "r"(gvl)
        :"v0");
#else
asm volatile(
        "vsetvli    zero, zero, e8, m1\n\t"
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e32,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+vd"(v0)
        :"vd"(mask0), "f"(zero), "r"(gvl)
        :"v0");
#endif
                        v_res = VFREDMINVS_FLOAT(v_res, v0, v_max, gvl);
                        if(*((FLOAT*)&v_res) < minf)
                                minf = *((FLOAT*)&v_res);
                        j += gvl;
                }
        }else{
                gvl = VSETVL(n);
                BLASLONG stride_x = inc_x * sizeof(FLOAT);
                if(gvl <= n/2){
                        BLASLONG idx = 0, inc_xv = inc_x * gvl;
                        v_min = VFMVVF_FLOAT(FLT_MAX, gvl);
                        for(i=0,j=0; i<n/(gvl*2); i++){
                                v0 = VLSEV_FLOAT(&x[idx], stride_x, gvl);
                                mask0 = VMFLTVF_FLOAT(v0, 0, gvl);
                                //v0 = VFRSUBVF_MASK_FLOAT(v0, 0, mask0, gvl);
#if defined(DOUBLE)
asm volatile(
        "vsetvli    zero, zero, e8, m1\n\t"
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e64,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+vd"(v0)
        :"vd"(mask0), "f"(zero), "r"(gvl)
        :"v0");
#else
asm volatile(
        "vsetvli    zero, zero, e8, m1\n\t"
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e32,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+vd"(v0)
        :"vd"(mask0), "f"(zero), "r"(gvl)
        :"v0");
#endif
                                v_min = VFMINVV_FLOAT(v_min, v0, gvl);

                                v1 = VLSEV_FLOAT(&x[idx+inc_xv], stride_x, gvl);
                                mask1 = VMFLTVF_FLOAT(v1, 0, gvl);
                                //v1 = VFRSUBVF_MASK_FLOAT(v1, 0, mask1, gvl);
#if defined(DOUBLE)
asm volatile(
        "vsetvli    zero, zero, e8, m1\n\t"
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e64,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+vd"(v1)
        :"vd"(mask1), "f"(zero), "r"(gvl)
        :"v0");
#else
asm volatile(
        "vsetvli    zero, zero, e8, m1\n\t"
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e32,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+vd"(v1)
        :"vd"(mask1), "f"(zero), "r"(gvl)
        :"v0");
#endif

                                v_min = VFMINVV_FLOAT(v_min, v1, gvl);
                                j += gvl*2;
                                idx += inc_xv*2;
                        }
                        v_res = VFREDMINVS_FLOAT(v_res, v_min, v_max, gvl);
                        minf = *((FLOAT*)&v_res);
                }
                for(;j<n;){
                        gvl = VSETVL(n-j);
                        v0 = VLSEV_FLOAT(&x[j*inc_x], stride_x, gvl);
                        mask0 = VMFLTVF_FLOAT(v0, 0, gvl);
                        //v0 = VFRSUBVF_MASK_FLOAT(v0, 0, mask0, gvl);
#if defined(DOUBLE)
asm volatile(
        "vsetvli    zero, zero, e8, m1\n\t"
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e64,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+vd"(v0)
        :"vd"(mask0), "f"(zero), "r"(gvl)
        :"v0");
#else
asm volatile(
        "vsetvli    zero, zero, e8, m1\n\t"
        "vor.vv     v0, %1, %1\n\t"
        "vsetvli    x0, %3, e32,m8 \n\t"
        "vfrsub.vf  %0, %0, %2, v0.t \n\t"
        :"+vd"(v0)
        :"vd"(mask0), "f"(zero), "r"(gvl)
        :"v0");
#endif
                        v_res = VFREDMINVS_FLOAT(v_res, v0, v_max, gvl);
                        if(*((FLOAT*)&v_res) < minf)
                                minf = *((FLOAT*)&v_res);
                        j += gvl;
                }
        }
	return(minf);
}


