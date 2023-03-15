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
#define VLEV_FLOAT vle32_v_f32m2
#define VSEV_FLOAT vse32_v_f32m2
#define VLSEV_FLOAT vlse32_v_f32m2
#define VLSEG2_FLOAT vlseg2e32_v_f32m2
#define VLSSEG2_FLOAT vlsseg2e32_v_f32m2
#define VSSEG2_FLOAT vsseg2e32_v_f32m2
#define VBOOL_T vbool16_t
#define UINT_V_T vint32m2_t
#define VID_V_UINT vid_v_i32m2
#define VMSGTU_VX_UINT vmsgt_vx_i32m2_b16
#define VMSEQ_VX_UINT vmseq_vx_i32m2_b16
#define VFMERGE_VFM_FLOAT  vfmerge_vfm_f32m2
#else
#define VSETVL(n) vsetvl_e64m2(n)
#define FLOAT_V_T vfloat64m2_t
#define VLEV_FLOAT vle64_v_f64m2
#define VSEV_FLOAT vse64_v_f64m2
#define VLSEV_FLOAT vlse64_v_f64m2
#define VLSEG2_FLOAT vlseg2e64_v_f64m2
#define VLSSEG2_FLOAT vlsseg2e64_v_f64m2
#define VSSEG2_FLOAT vsseg2e64_v_f64m2
#define VBOOL_T     vbool32_t
#define UINT_V_T     vuint64m2_t
#define VID_V_UINT   vid_v_u64m2
#define VMSGTU_VX_UINT vmsgtu_vx_u64m2_b32
#define VMSEQ_VX_UINT vmseq_vx_u64m2_b32
#define VFMERGE_VFM_FLOAT  vfmerge_vfm_f64m2
#endif

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, BLASLONG posX, BLASLONG posY, FLOAT *b){

    BLASLONG i, js, X;

    FLOAT *ao;

    BLASLONG stride_lda = sizeof(FLOAT)*lda*2;
    
    FLOAT_V_T va0, va1;

    size_t vl;
#ifdef UNIT
    VBOOL_T vbool_eq;
#endif

    VBOOL_T vbool_cmp;
    UINT_V_T vindex;

    for (js = n; js > 0; js -= vl)
    {
        vl = VSETVL(js);
        X = posX;

        if (posX <= posY) 
        {
            ao = a + posY * 2 + posX * lda * 2;
        } 
        else 
        {
            ao = a + posX * 2 + posY * lda * 2;
        }

        i = 0;
        do 
        {
            if (X > posY) 
            {
                VLSSEG2_FLOAT(&va0, &va1, ao, stride_lda, vl);
                VSSEG2_FLOAT(b, va0, va1, vl);

                ao  += 2;
                b   += vl * 2;

                X ++;
                i ++;
            } 
            else if (X < posY) 
            {
                ao  += lda * 2;
                b   += vl * 2;
                X ++;
                i ++;
            } 
            else 
            {
                vindex  = VID_V_UINT(vl);
                for (unsigned int j = 0; j < vl; j++) 
                {
                    VLSSEG2_FLOAT(&va0, &va1, ao, stride_lda, vl);
                    vbool_cmp = VMSGTU_VX_UINT(vindex, j, vl);
                    va0 = VFMERGE_VFM_FLOAT(vbool_cmp, va0, ZERO, vl);
                    va1 = VFMERGE_VFM_FLOAT(vbool_cmp, va1, ZERO, vl);
#ifdef UNIT
                    vbool_eq = VMSEQ_VX_UINT(vindex, j, vl);
                    va0 =  VFMERGE_VFM_FLOAT(vbool_eq, va0, ONE, vl);
                    va1 =  VFMERGE_VFM_FLOAT(vbool_eq, va1, ZERO, vl);
#endif
                    VSSEG2_FLOAT(b, va0, va1, vl);
                    ao  += 2;
                    b   += vl * 2;
                }

                X += vl;
                i += vl;
            }
        } while (i < m);

        posY += vl;
    }

    return 0;
}
