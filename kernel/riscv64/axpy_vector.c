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
#define VSETVL(n) vsetvl_e32m4(n)
#define FLOAT_V_T vfloat32m4_t
#define VLEV_FLOAT vle32_v_f32m4
#define VLSEV_FLOAT vlse32_v_f32m4
#define VSEV_FLOAT vse32_v_f32m4
#define VSSEV_FLOAT vsse32_v_f32m4
#define VFMACCVF_FLOAT vfmacc_vf_f32m4
#else
#define VSETVL(n) vsetvl_e64m4(n)
#define FLOAT_V_T vfloat64m4_t
#define VLEV_FLOAT vle64_v_f64m4
#define VLSEV_FLOAT vlse64_v_f64m4
#define VSEV_FLOAT vse64_v_f64m4
#define VSSEV_FLOAT vsse64_v_f64m4
#define VFMACCVF_FLOAT vfmacc_vf_f64m4
#endif

int CNAME(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT da, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *dummy, BLASLONG dummy2)
{
	BLASLONG i=0, j=0, jx=0, jy=0;
	unsigned int gvl = 0;
	FLOAT_V_T vx0, vx1;
	FLOAT_V_T vy0, vy1;
	BLASLONG stride_x, stride_y;

	if (n < 0)  return(0);
	if (da == 0.0) return(0);

	if (inc_x == 1 && inc_y == 1) {

		gvl = VSETVL(n);

		if (gvl <= n/2) {
			for (i = 0, j=0; i < n/(2*gvl); i++, j+=2*gvl) {
				vx0 = VLEV_FLOAT(&x[j], gvl);
				vy0 = VLEV_FLOAT(&y[j], gvl);
				vy0 = VFMACCVF_FLOAT(vy0, da, vx0, gvl);
				VSEV_FLOAT(&y[j], vy0, gvl);

				vx1 = VLEV_FLOAT(&x[j+gvl], gvl);
				vy1 = VLEV_FLOAT(&y[j+gvl], gvl);
				vy1 = VFMACCVF_FLOAT(vy1, da, vx1, gvl);
				VSEV_FLOAT(&y[j+gvl], vy1, gvl);
			}
		}
		//tail
		for (; j < n; ) {
			gvl = VSETVL(n - j);
			vx0 = VLEV_FLOAT(&x[j], gvl);
			vy0 = VLEV_FLOAT(&y[j], gvl);
			vy0 = VFMACCVF_FLOAT(vy0, da, vx0, gvl);
			VSEV_FLOAT(&y[j], vy0, gvl);

			j += gvl;
		}
	}else if (inc_y == 1) {
		stride_x = inc_x * sizeof(FLOAT);
                gvl = VSETVL(n);
                if(gvl <= n/2){
                        BLASLONG inc_xv = inc_x * gvl;
                        for(i=0,j=0; i<n/(2*gvl); i++){
			        vx0 = VLSEV_FLOAT(&x[jx], stride_x, gvl);
                                vy0 = VLEV_FLOAT(&y[j], gvl);
                                vy0 = VFMACCVF_FLOAT(vy0, da, vx0, gvl);
                                VSEV_FLOAT(&y[j], vy0, gvl);

			        vx1 = VLSEV_FLOAT(&x[jx+inc_xv], stride_x, gvl);
                                vy1 = VLEV_FLOAT(&y[j+gvl], gvl);
                                vy1 = VFMACCVF_FLOAT(vy1, da, vx1, gvl);
                                VSEV_FLOAT(&y[j+gvl], vy1, gvl);

                                j += gvl * 2;
                                jx += inc_xv * 2;
                        }
                }
		for (; j<n; ) {
			gvl = VSETVL(n - j);
			vx0 = VLSEV_FLOAT(&x[j*inc_x], stride_x, gvl);
			vy0 = VLEV_FLOAT(&y[j], gvl);
			vy0 = VFMACCVF_FLOAT(vy0, da, vx0, gvl);
			VSEV_FLOAT(&y[j], vy0, gvl);
			j += gvl;
		}
        }else if(inc_x == 1){
		stride_y = inc_y * sizeof(FLOAT);
                gvl = VSETVL(n);
                if(gvl <= n/2){
                        BLASLONG inc_yv = inc_y * gvl;
                        for(i=0,j=0; i<n/(2*gvl); i++){
			        vx0 = VLEV_FLOAT(&x[j], gvl);
                                vy0 = VLSEV_FLOAT(&y[jy], stride_y, gvl);
                                vy0 = VFMACCVF_FLOAT(vy0, da, vx0, gvl);
                                VSSEV_FLOAT(&y[jy], stride_y, vy0, gvl);

			        vx1 = VLEV_FLOAT(&x[j+gvl], gvl);
                                vy1 = VLSEV_FLOAT(&y[jy+inc_yv], stride_y, gvl);
                                vy1 = VFMACCVF_FLOAT(vy1, da, vx1, gvl);
                                VSSEV_FLOAT(&y[jy+inc_yv], stride_y, vy1, gvl);

                                j += gvl * 2;
                                jy += inc_yv * 2;
                        }
                }
		for (; j<n; ) {
			gvl = VSETVL(n - j);
			vx0 = VLEV_FLOAT(&x[j], gvl);
			vy0 = VLSEV_FLOAT(&y[j*inc_y], stride_y, gvl);
			vy0 = VFMACCVF_FLOAT(vy0, da, vx0, gvl);
			VSSEV_FLOAT(&y[j*inc_y], stride_y, vy0, gvl);
			j += gvl;
		}
	}else{
		stride_x = inc_x * sizeof(FLOAT);
		stride_y = inc_y * sizeof(FLOAT);
                gvl = VSETVL(n);
                if(gvl <= n/2){
                        BLASLONG inc_xv = inc_x * gvl;
                        BLASLONG inc_yv = inc_y * gvl;
                        for(i=0,j=0; i<n/(2*gvl); i++){
			        vx0 = VLSEV_FLOAT(&x[jx], stride_x, gvl);
                                vy0 = VLSEV_FLOAT(&y[jy], stride_y, gvl);
                                vy0 = VFMACCVF_FLOAT(vy0, da, vx0, gvl);
                                VSSEV_FLOAT(&y[jy], stride_y, vy0, gvl);

			        vx1 = VLSEV_FLOAT(&x[jx+inc_xv], stride_x, gvl);
                                vy1 = VLSEV_FLOAT(&y[jy+inc_yv], stride_y, gvl);
                                vy1 = VFMACCVF_FLOAT(vy1, da, vx1, gvl);
                                VSSEV_FLOAT(&y[jy+inc_yv], stride_y, vy1, gvl);

                                j += gvl * 2;
                                jx += inc_xv * 2;
                                jy += inc_yv * 2;
                        }
                }
		for (; j<n; ) {
			gvl = VSETVL(n - j);
			vx0 = VLSEV_FLOAT(&x[j*inc_x], stride_x, gvl);
			vy0 = VLSEV_FLOAT(&y[j*inc_y], stride_y, gvl);
			vy0 = VFMACCVF_FLOAT(vy0, da, vx0, gvl);
			VSSEV_FLOAT(&y[j*inc_y], stride_y, vy0, gvl);
			j += gvl;
		}
	}
	return(0);
}


