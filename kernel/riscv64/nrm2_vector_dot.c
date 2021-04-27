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
#define RVV_M RVV_M8
#define FLOAT_V_T float32xm8_t
#define VLEV_FLOAT vlev_float32xm8
#define VLSEV_FLOAT vlsev_float32xm8
#define VFREDSUM_FLOAT vfredsumvs_float32xm8
#define VFMACCVV_FLOAT vfmaccvv_float32xm8
#define VFMVVF_FLOAT vfmvvf_float32xm8
#define VFDOTVV_FLOAT vfdotvv_float32xm8
#define ABS fabsf
#else
#define RVV_EFLOAT RVV_E64
#define RVV_M RVV_M8
#define FLOAT_V_T float64xm8_t
#define VLEV_FLOAT vlev_float64xm8
#define VLSEV_FLOAT vlsev_float64xm8
#define VFREDSUM_FLOAT vfredsumvs_float64xm8
#define VFMACCVV_FLOAT vfmaccvv_float64xm8
#define VFMVVF_FLOAT vfmvvf_float64xm8
#define VFDOTVV_FLOAT vfdotvv_float64xm8
#define ABS fabs
#endif

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
	BLASLONG i=0, j=0;
	double len = 0.0 ;

	if ( n < 0 )  return(0.0);
        if(n == 1) return (ABS(x[0]));

        FLOAT_V_T vr, v0, v1;
        unsigned int gvl = 0;
        if(inc_x == 1){
                gvl = vsetvli(n, RVV_EFLOAT, RVV_M);
                if(gvl < n/2){
                        vr = VFMVVF_FLOAT(0, gvl);
                        for(i=0,j=0; i<n/(2*gvl); i++){
                                v0 = VLEV_FLOAT(&x[j], gvl);
                                vr = VFMACCVV_FLOAT(vr, v0, v0, gvl);
                                j += gvl;

                                v1 = VLEV_FLOAT(&x[j], gvl);
                                vr = VFMACCVV_FLOAT(vr, v1, v1, gvl);
                                j += gvl;
                        }
                        v0 = VFMVVF_FLOAT(0, gvl);
                        v0 = VFREDSUM_FLOAT(vr, v0, gvl);
                        len += v0[0];
                }
                //tail
                for(;j < n;){
                        gvl = vsetvli(n-j, RVV_EFLOAT, RVV_M);
                        v0 = VLEV_FLOAT(&x[j], gvl);
                        //v1 = 0
                        v1 = VFMVVF_FLOAT(0, gvl);
                        //vr = VFDOTVV_FLOAT(v0, v0, gvl);
                        vr = VFMACCVV_FLOAT(v1, v0, v0, gvl);
                        v0 = VFREDSUM_FLOAT(vr, v1, gvl);
                        len += v0[0];

                        j += gvl;
                }
        }else{
                gvl = vsetvli(n, RVV_EFLOAT, RVV_M);
                unsigned int stride_x = inc_x * sizeof(FLOAT);
                if(gvl < n/2){
                        vr = VFMVVF_FLOAT(0, gvl);
                        for(i=0,j=0; i<n/(2*gvl); i++){
                                v0 = VLSEV_FLOAT(&x[j*inc_x], stride_x, gvl);
                                vr = VFMACCVV_FLOAT(vr, v0, v0, gvl);
                                j += gvl;

                                v1 = VLSEV_FLOAT(&x[j*inc_x], stride_x, gvl);
                                vr = VFMACCVV_FLOAT(vr, v1, v1, gvl);
                                j += gvl;
                        }
                        v0 = VFMVVF_FLOAT(0, gvl);
                        v0 = VFREDSUM_FLOAT(vr, v0, gvl);
                        len += v0[0];
                }
                //tail
                for(;j < n;){
                        gvl = vsetvli(n-j, RVV_EFLOAT, RVV_M);
                        v0 = VLSEV_FLOAT(&x[j*inc_x], stride_x, gvl);
                        //v1 = 0
                        v1 = VFMVVF_FLOAT(0, gvl);
                        //vr = VFDOTVV_FLOAT(v0, v0, gvl);
                        vr = VFMACCVV_FLOAT(v1, v0, v0, gvl);
                        v0 = VFREDSUM_FLOAT(vr, v1, gvl);
                        len += v0[0];

                        j += gvl;
                }
        }
	return(sqrt(len));
}


