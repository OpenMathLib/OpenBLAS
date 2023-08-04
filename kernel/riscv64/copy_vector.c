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
#define VSETVL(n) vsetvl_e32m8(n)
#define FLOAT_V_T vfloat32m8_t
#define VLEV_FLOAT vle32_v_f32m8
#define VLSEV_FLOAT vlse32_v_f32m8
#define VSEV_FLOAT vse32_v_f32m8
#define VSSEV_FLOAT vsse32_v_f32m8
#else
#define VSETVL(n) vsetvl_e64m8(n)
#define FLOAT_V_T vfloat64m8_t
#define VLEV_FLOAT vle64_v_f64m8
#define VLSEV_FLOAT vlse64_v_f64m8
#define VSEV_FLOAT vse64_v_f64m8
#define VSSEV_FLOAT vsse64_v_f64m8
#endif

int CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y)
{
	BLASLONG i=0, j=0;
	BLASLONG ix=0,iy=0;
	if(n < 0)  return(0);

        BLASLONG stride_x, stride_y;
        FLOAT_V_T v0, v1, v2, v3;
        unsigned int gvl = 0;

        if(inc_x == 1 && inc_y == 1){
                memcpy(&y[0], &x[0], n*sizeof(FLOAT));
        }else if (inc_y == 1){
                gvl = VSETVL(n);
                stride_x = inc_x * sizeof(FLOAT);
                if(gvl <= n/4){
                        BLASLONG inc_xv = inc_x * gvl;
                        BLASLONG gvl3 = gvl * 3;
                        BLASLONG inc_xv3 = inc_xv * 3;
                        for(i=0,j=0; i<n/(4*gvl); i++){
                                v0 = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                                VSEV_FLOAT(&y[j], v0, gvl);
                                v1 = VLSEV_FLOAT(&x[ix+inc_xv], stride_x, gvl);
                                VSEV_FLOAT(&y[j+gvl], v1, gvl);

                                v2 = VLSEV_FLOAT(&x[ix+inc_xv*2], stride_x, gvl);
                                VSEV_FLOAT(&y[j+gvl*2], v2, gvl);
                                v3 = VLSEV_FLOAT(&x[ix+inc_xv3], stride_x, gvl);
                                VSEV_FLOAT(&y[j+gvl3], v3, gvl);
                                j += gvl * 4;
                                ix += inc_xv * 4;
                        }
                }
                for(;j<n;){
                        gvl = VSETVL(n-j);
                        v0 = VLSEV_FLOAT(&x[j*inc_x], stride_x, gvl);
                        VSEV_FLOAT(&y[j], v0, gvl);
                        j += gvl;
                }
        }else if(inc_x == 1){
                gvl = VSETVL(n);
                stride_y = inc_y * sizeof(FLOAT);
                if(gvl <= n/4){
                        BLASLONG inc_yv = inc_y * gvl;
                        BLASLONG inc_yv3 = inc_yv * 3;
                        BLASLONG gvl3 = gvl * 3;
                        for(i=0,j=0; i<n/(4*gvl); i++){
                                v0 = VLEV_FLOAT(&x[j], gvl);
                                VSSEV_FLOAT(&y[iy], stride_y, v0, gvl);
                                v1 = VLEV_FLOAT(&x[j+gvl], gvl);
                                VSSEV_FLOAT(&y[iy+inc_yv], stride_y, v1, gvl);

                                v2 = VLEV_FLOAT(&x[j+gvl*2], gvl);
                                VSSEV_FLOAT(&y[iy+inc_yv*2], stride_y, v2, gvl);
                                v3 = VLEV_FLOAT(&x[j+gvl3], gvl);
                                VSSEV_FLOAT(&y[iy+inc_yv3], stride_y, v3, gvl);
                                j += gvl * 4;
                                iy += inc_yv * 4;
                        }
                }
                for(;j<n;){
                        gvl = VSETVL(n-j);
                        v0 = VLEV_FLOAT(&x[j], gvl);
                        VSSEV_FLOAT(&y[j*inc_y], stride_y, v0, gvl);
                        j += gvl;
                }

        }else{
                gvl = VSETVL(n);
                stride_x = inc_x * sizeof(FLOAT);
                stride_y = inc_y * sizeof(FLOAT);
                if(gvl <= n/4){
                        BLASLONG inc_xv = inc_x * gvl;
                        BLASLONG inc_yv = inc_y * gvl;
                        BLASLONG inc_xv3 = inc_xv * 3;
                        BLASLONG inc_yv3 = inc_yv * 3;
                        for(i=0,j=0; i<n/(4*gvl); i++){
                                v0 = VLSEV_FLOAT(&x[ix], stride_x, gvl);
                                VSSEV_FLOAT(&y[iy], stride_y, v0, gvl);
                                v1 = VLSEV_FLOAT(&x[ix+inc_xv], stride_x, gvl);
                                VSSEV_FLOAT(&y[iy+inc_yv], stride_y, v1, gvl);

                                v2 = VLSEV_FLOAT(&x[ix+inc_xv*2], stride_x, gvl);
                                VSSEV_FLOAT(&y[iy+inc_yv*2], stride_y, v2, gvl);
                                v3 = VLSEV_FLOAT(&x[ix+inc_xv3], stride_x, gvl);
                                VSSEV_FLOAT(&y[iy+inc_yv3], stride_y, v3, gvl);

                                j += gvl * 4;
                                ix += inc_xv * 4;
                                iy += inc_yv * 4;
                        }
                }
                for(;j<n;){
                        gvl = VSETVL(n-j);
                        v0 = VLSEV_FLOAT(&x[j*inc_x], stride_x, gvl);
                        VSSEV_FLOAT(&y[j*inc_y], stride_y, v0, gvl);
                        j += gvl;
                }
        }
	return(0);
}


