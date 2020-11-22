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
#define VFREDSUMVS_FLOAT vfredsumvs_float32xm8
#define MASK_T e32xm8_t
#define VMFLTVF_FLOAT vmfltvf_e32xm8_float32xm8
#define VFMVVF_FLOAT vfmvvf_float32xm8
#define VFRSUBVF_MASK_FLOAT vfrsubvf_mask_float32xm8
#define VFADDVV_FLOAT vfaddvv_float32xm8
#else
#define RVV_EFLOAT RVV_E64
#define RVV_M RVV_M8
#define FLOAT_V_T float64xm8_t
#define VLEV_FLOAT vlev_float64xm8
#define VLSEV_FLOAT vlsev_float64xm8
#define VFREDSUMVS_FLOAT vfredsumvs_float64xm8
#define MASK_T e64xm8_t
#define VMFLTVF_FLOAT vmfltvf_e64xm8_float64xm8
#define VFMVVF_FLOAT vfmvvf_float64xm8
#define VFRSUBVF_MASK_FLOAT vfrsubvf_mask_float64xm8
#define VFADDVV_FLOAT vfaddvv_float64xm8
#endif
FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
	BLASLONG i=0, j=0;
	BLASLONG ix=0;
	FLOAT asumf=0.0;
	if (n <= 0 || inc_x <= 0) return(asumf);
        unsigned int gvl = 0;
        FLOAT_V_T v0, v1, v_zero,v_sum;

        MASK_T mask0, mask1;
        if(inc_x == 1){
                gvl = vsetvli(n, RVV_EFLOAT, RVV_M);
                v_zero = VFMVVF_FLOAT(0, gvl);
                if(gvl <= n/2){
                        v_sum = VFMVVF_FLOAT(0, gvl);
                        for(i=0,j=0; i<n/(gvl*2); i++){
                                v0 = VLEV_FLOAT(&x[j], gvl);
                                mask0 = VMFLTVF_FLOAT(v0, 0, gvl);
                                v0 = VFRSUBVF_MASK_FLOAT(v0, v0, 0, mask0, gvl);
                                v_sum = VFADDVV_FLOAT(v_sum, v0, gvl);

                                v1 = VLEV_FLOAT(&x[j+gvl], gvl);
                                mask1 = VMFLTVF_FLOAT(v1, 0, gvl);
                                v1 = VFRSUBVF_MASK_FLOAT(v1, v1, 0, mask1, gvl);
                                v_sum = VFADDVV_FLOAT(v_sum, v1, gvl);
                                j += gvl * 2;
                        }
                        v0 = VFREDSUMVS_FLOAT(v_sum, v_zero, gvl);
                        asumf += v0[0];
                }
                for(;j<n;){
                        gvl = vsetvli(n-j, RVV_EFLOAT, RVV_M);
                        v0 = VLEV_FLOAT(&x[j], gvl);
                        mask0 = VMFLTVF_FLOAT(v0, 0, gvl);
                        v0 = VFRSUBVF_MASK_FLOAT(v0, v0, 0, mask0, gvl);
                        v0 = VFREDSUMVS_FLOAT(v0, v_zero, gvl);
                        asumf += v0[0];
                        j += gvl;
                }
        }else{
                gvl = vsetvli(n, RVV_EFLOAT, RVV_M);
                unsigned int stride_x = inc_x * sizeof(FLOAT);
                v_zero = VFMVVF_FLOAT(0, gvl);
                if(gvl <= n/2){
                        v_sum = VFMVVF_FLOAT(0, gvl);
                        BLASLONG inc_xv = inc_x * gvl;
                        for(i=0,j=0; i<n/(gvl*2); i++){
                                v0 = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                                mask0 = VMFLTVF_FLOAT(v0, 0, gvl);
                                v0 = VFRSUBVF_MASK_FLOAT(v0, v0, 0, mask0, gvl);
                                v_sum = VFADDVV_FLOAT(v_sum, v0, gvl);

                                v1 = VLSEV_FLOAT(&x[ix+inc_xv], stride_x, gvl);
                                mask1 = VMFLTVF_FLOAT(v1, 0, gvl);
                                v1 = VFRSUBVF_MASK_FLOAT(v1, v1, 0, mask1, gvl);
                                v_sum = VFADDVV_FLOAT(v_sum, v1, gvl);
                                j += gvl * 2;
                                inc_xv += inc_xv * 2;
                        }
                        v0 = VFREDSUMVS_FLOAT(v_sum, v_zero, gvl);
                        asumf += v0[0];
                }
                for(;j<n;){
                        gvl = vsetvli(n-j, RVV_EFLOAT, RVV_M);
                        v0 = VLSEV_FLOAT(&x[j*inc_x], stride_x, gvl);
                        mask0 = VMFLTVF_FLOAT(v0, 0, gvl);
                        v0 = VFRSUBVF_MASK_FLOAT(v0, v0, 0, mask0, gvl);
                        v0 = VFREDSUMVS_FLOAT(v0, v_zero, gvl);
                        asumf += v0[0];
                        j += gvl;
                }
        }
	return(asumf);
}


