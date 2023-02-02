/***************************************************************************
Copyright (c) 2022, The OpenBLAS Project
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

#include <stdio.h>
#include "common.h"

#if !defined(DOUBLE)
#define VSETVL(n) vsetvl_e32m2(n)
#define FLOAT_V_T vfloat32m2_t
#define VLSSEG2_FLOAT vlsseg2e32_v_f32m2
#define VSSEG2_FLOAT vsseg2e32_v_f32m2
#define VSSEG2_FLOAT_M vsseg2e32_v_f32m2_m
#define VBOOL_T vbool16_t
#define UINT_V_T vuint32m2_t
#define VID_V_UINT vid_v_u32m2
#define VMSLTU_VX_UINT vmsltu_vx_u32m2_b16
#else
#define VSETVL(n) vsetvl_e64m2(n)
#define FLOAT_V_T vfloat64m2_t
#define VLSSEG2_FLOAT vlsseg2e64_v_f64m2
#define VSSEG2_FLOAT vsseg2e64_v_f64m2
#define VSSEG2_FLOAT_M vsseg2e64_v_f64m2_m
#define VBOOL_T     vbool32_t
#define UINT_V_T     vuint64m2_t
#define VID_V_UINT   vid_v_u64m2
#define VMSLTU_VX_UINT vmsltu_vx_u64m2_b32

#endif


int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, BLASLONG offset, FLOAT *b){

    //fprintf(stderr, "%s , %s, m = %4ld  n = %4ld  lda = %4ld offset = %4ld\n", __FILE__, __FUNCTION__, m, n, lda, offset); // Debug

    BLASLONG i, ii, jj, js;

    FLOAT *ao;

    jj = offset;

    BLASLONG stride_lda = sizeof(FLOAT)*lda*2;

    FLOAT_V_T va0, va1;
    VBOOL_T vbool_cmp;
    UINT_V_T vindex;
    size_t vl;

    for (js = n; js > 0; js -= vl)
    {
        vl = VSETVL(js);
        ao = a;

        ii = 0;
        for (i = 0; i < m;)
        {
            if (ii == jj) 
            {
                vindex  = VID_V_UINT(vl);
                for (unsigned int j = 0; j < vl; j++) 
                {
                    VLSSEG2_FLOAT(&va0, &va1, ao, stride_lda, vl);
                    vbool_cmp = VMSLTU_VX_UINT(vindex, j, vl);
                    VSSEG2_FLOAT_M(vbool_cmp, b, va0, va1, vl);

                    compinv((b + j * 2), *(ao + j * lda * 2), *(ao + j * lda * 2 + 1));
                    ao  += 2;
                    b   += vl * 2;
                }
                i += vl;
                ii += vl;
            }
            else
            {
                if (ii > jj)
                {
                    VLSSEG2_FLOAT(&va0, &va1, ao, stride_lda, vl);
                    VSSEG2_FLOAT(b, va0, va1, vl);
                }
                ao  += 2;
                b   += vl * 2;
                i++;
                ii++;
            }
        }

        a += vl * lda * 2;
        jj += vl;
    }

    return 0;
}
