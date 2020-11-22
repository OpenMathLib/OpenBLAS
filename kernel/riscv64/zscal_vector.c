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
#define VLSEV_FLOAT vlsev_float32xm4
#define VSSEV_FLOAT vssev_float32xm4
#define VFMACCVF_FLOAT vfmaccvf_float32xm4
#define VFMULVF_FLOAT vfmulvf_float32xm4
#define VFNMSACVF_FLOAT vfnmsacvf_float32xm4
#else
#define RVV_EFLOAT RVV_E64
#define RVV_M RVV_M4
#define FLOAT_V_T float64xm4_t
#define VLSEV_FLOAT vlsev_float64xm4
#define VSSEV_FLOAT vssev_float64xm4
#define VFMACCVF_FLOAT vfmaccvf_float64xm4
#define VFMULVF_FLOAT vfmulvf_float64xm4
#define VFNMSACVF_FLOAT vfnmsacvf_float64xm4
#endif

int CNAME(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT da_r,FLOAT da_i, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *dummy, BLASLONG dummy2)
{
        BLASLONG i=0, j=0;
        BLASLONG ix=0;


        if((n <= 0) || (inc_x <= 0))
                return(0);

        unsigned int gvl = 0;
        FLOAT_V_T vt, v0, v1;
        if(da_r == 0.0 && da_i == 0.0){
                memset(&x[0], 0, n * 2 * sizeof(FLOAT));
        }else if(da_r == 0.0){
                gvl = vsetvli(n, RVV_EFLOAT, RVV_M);
                BLASLONG stride_x = inc_x * 2 * sizeof(FLOAT);
                BLASLONG inc_xv = inc_x * 2 * gvl;
                for(i=0,j=0; i < n/gvl; i++){
                        v0 = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                        v1 = VLSEV_FLOAT(&x[ix+1], stride_x, gvl);

                        vt = VFMULVF_FLOAT(v1, -da_i, gvl);
                        v1 = VFMULVF_FLOAT(v0, da_i, gvl);

                        VSSEV_FLOAT(&x[ix], stride_x, vt, gvl);
                        VSSEV_FLOAT(&x[ix+1], stride_x, v1, gvl);

                        j += gvl;
                        ix += inc_xv;
                }
                if(j < n){
                        gvl = vsetvli(n-j, RVV_EFLOAT, RVV_M);
                        v0 = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                        v1 = VLSEV_FLOAT(&x[ix+1], stride_x, gvl);

                        vt = VFMULVF_FLOAT(v1, -da_i, gvl);
                        v1 = VFMULVF_FLOAT(v0, da_i, gvl);

                        VSSEV_FLOAT(&x[ix], stride_x, vt, gvl);
                        VSSEV_FLOAT(&x[ix+1], stride_x, v1, gvl);
                }
        }else if(da_i == 0.0){
                gvl = vsetvli(n, RVV_EFLOAT, RVV_M);
                BLASLONG stride_x = inc_x * 2 * sizeof(FLOAT);
                BLASLONG inc_xv = inc_x * 2 * gvl;
                for(i=0,j=0; i < n/gvl; i++){
                        v0 = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                        v1 = VLSEV_FLOAT(&x[ix+1], stride_x, gvl);

                        vt = VFMULVF_FLOAT(v0, da_r, gvl);
                        v1 = VFMULVF_FLOAT(v1, da_r, gvl);

                        VSSEV_FLOAT(&x[ix], stride_x, vt, gvl);
                        VSSEV_FLOAT(&x[ix+1], stride_x, v1, gvl);

                        j += gvl;
                        ix += inc_xv;
                }
                if(j < n){
                        gvl = vsetvli(n-j, RVV_EFLOAT, RVV_M);
                        v0 = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                        v1 = VLSEV_FLOAT(&x[ix+1], stride_x, gvl);

                        vt = VFMULVF_FLOAT(v0, da_r, gvl);
                        v1 = VFMULVF_FLOAT(v1, da_r, gvl);

                        VSSEV_FLOAT(&x[ix], stride_x, vt, gvl);
                        VSSEV_FLOAT(&x[ix+1], stride_x, v1, gvl);
                }
        }else{
                gvl = vsetvli(n, RVV_EFLOAT, RVV_M);
                BLASLONG stride_x = inc_x * 2 * sizeof(FLOAT);
                BLASLONG inc_xv = inc_x * 2 * gvl;
                for(i=0,j=0; i < n/gvl; i++){
                        v0 = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                        v1 = VLSEV_FLOAT(&x[ix+1], stride_x, gvl);

                        vt = VFMULVF_FLOAT(v0, da_r, gvl);
                        vt = VFNMSACVF_FLOAT(vt, da_i, v1, gvl);
                        v1 = VFMULVF_FLOAT(v1, da_r, gvl);
                        v1 = VFMACCVF_FLOAT(v1, da_i, v0, gvl);

                        VSSEV_FLOAT(&x[ix], stride_x, vt, gvl);
                        VSSEV_FLOAT(&x[ix+1], stride_x, v1, gvl);

                        j += gvl;
                        ix += inc_xv;
                }
                if(j < n){
                        gvl = vsetvli(n-j, RVV_EFLOAT, RVV_M);
                        v0 = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                        v1 = VLSEV_FLOAT(&x[ix+1], stride_x, gvl);

                        vt = VFMULVF_FLOAT(v0, da_r, gvl);
                        vt = VFNMSACVF_FLOAT(vt, da_i, v1, gvl);
                        v1 = VFMULVF_FLOAT(v1, da_r, gvl);
                        v1 = VFMACCVF_FLOAT(v1, da_i, v0, gvl);

                        VSSEV_FLOAT(&x[ix], stride_x, vt, gvl);
                        VSSEV_FLOAT(&x[ix+1], stride_x, v1, gvl);
                }
        }
        return(0);
}
