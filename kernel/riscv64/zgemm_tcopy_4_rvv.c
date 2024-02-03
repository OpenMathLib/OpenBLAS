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

#include "common.h"

#if !defined(DOUBLE)
#define VSETVL(n)               __riscv_vsetvl_e32m1(n)
#define FLOAT_V_T               vfloat32m1_t
#define FLOAT_VX2_T             vfloat32m1x2_t
#define FLOAT_VX4_T             vfloat32m1x4_t
#define FLOAT_VX8_T             vfloat32m1x8_t
#define VLEV_FLOAT              __riscv_vle32_v_f32m1
#define VSEV_FLOAT              __riscv_vse32_v_f32m1
#define VLSSEG2_FLOAT           __riscv_vlsseg2e32_v_f32m1x2
#define VLSSEG4_FLOAT           __riscv_vlsseg4e32_v_f32m1x4
#define VLSSEG8_FLOAT           __riscv_vlsseg8e32_v_f32m1x8
#define VSSEG2_FLOAT            __riscv_vsseg2e32_v_f32m1x2
#define VSSEG4_FLOAT            __riscv_vsseg4e32_v_f32m1x4
#define VSSEG8_FLOAT            __riscv_vsseg8e32_v_f32m1x8
#else
#define VSETVL(n)               __riscv_vsetvl_e64m1(n)
#define FLOAT_V_T               vfloat64m1_t
#define FLOAT_VX2_T             vfloat64m1x2_t
#define FLOAT_VX4_T             vfloat64m1x4_t
#define FLOAT_VX8_T             vfloat64m1x8_t
#define VLEV_FLOAT              __riscv_vle64_v_f64m1
#define VSEV_FLOAT              __riscv_vse64_v_f64m1
#define VLSSEG2_FLOAT           __riscv_vlsseg2e64_v_f64m1x2
#define VLSSEG4_FLOAT           __riscv_vlsseg4e64_v_f64m1x4
#define VLSSEG8_FLOAT           __riscv_vlsseg8e64_v_f64m1x8
#define VSSEG2_FLOAT            __riscv_vsseg2e64_v_f64m1x2
#define VSSEG4_FLOAT            __riscv_vsseg4e64_v_f64m1x4
#define VSSEG8_FLOAT            __riscv_vsseg8e64_v_f64m1x8
#endif

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, FLOAT *b){

    BLASLONG i, j;

    IFLOAT *aoffset;
    IFLOAT *aoffset1;

    IFLOAT *boffset, *boffset1, *boffset2, *boffset3;

    FLOAT_V_T v0;
    FLOAT_VX2_T vx2;
    FLOAT_VX4_T vx4;
    FLOAT_VX8_T vx8;

    size_t vl;

    //fprintf(stderr, "%s m=%ld n=%ld lda=%ld\n", __FUNCTION__, m, n, lda);

    aoffset   = a;
    boffset   = b;
    boffset2  = b + 2 * m  * (n & ~3);
    boffset3  = b + 2 * m  * (n & ~1);

    for(j = (m >> 2); j > 0; j--) {

        aoffset1  = aoffset;
        aoffset += 8 * lda;

        boffset1  = boffset;
        boffset  += 32;

        for(i = (n >> 2); i > 0; i--) {
            vl = 4;

            vx8 = VLSSEG8_FLOAT(aoffset1, lda * sizeof(FLOAT) * 2, vl);
            VSSEG8_FLOAT(boffset1, vx8, vl);

            aoffset1 += 8;
            boffset1 += m * 8;
        }

        if (n & 2) {
            vl = 4;

            vx4 = VLSSEG4_FLOAT(aoffset1, lda * sizeof(FLOAT) * 2, vl);
            VSSEG4_FLOAT(boffset2, vx4, vl);

            aoffset1 += 4;
            boffset2 += 16;
        }

        if (n & 1) {
            vl = 4;

            vx2 = VLSSEG2_FLOAT(aoffset1, lda * sizeof(FLOAT) * 2, vl);
            VSSEG2_FLOAT(boffset3, vx2, vl);

            aoffset1 += 2;
            boffset3 += 8;
        }
    }

    if (m & 2) {
        aoffset1  = aoffset;
        aoffset += 4 * lda;

        boffset1  = boffset;
        boffset  += 16;

        for(i = (n >> 2); i > 0; i--) {
            vl = 2;

            vx8 = VLSSEG8_FLOAT(aoffset1, lda * sizeof(FLOAT) * 2, vl);
            VSSEG8_FLOAT(boffset1, vx8, vl);

            aoffset1 += 8;
            boffset1 += m * 8;
        }

        if (n & 2) {
            vl = 2;

            vx4 = VLSSEG4_FLOAT(aoffset1, lda * sizeof(FLOAT) * 2, vl);
            VSSEG4_FLOAT(boffset2, vx4, vl);

            aoffset1 += 4;
            boffset2 += 8;
        }

        if (n & 1) {
            vl = 2;

            vx2 = VLSSEG2_FLOAT(aoffset1, lda * sizeof(FLOAT) * 2, vl);
            VSSEG2_FLOAT(boffset3, vx2, vl);

            //aoffset1 += 2;
            boffset3 += 4;
        }
    }

    if (m & 1) {
        aoffset1  = aoffset;
        boffset1  = boffset;

        for(i = (n >> 2); i > 0; i--) {
            vl = 8;

            v0 = VLEV_FLOAT(aoffset1, vl);
            VSEV_FLOAT(boffset1, v0, vl);

            aoffset1 += 8;
            boffset1 += 8 * m;
        }

        if (n & 2) {
            vl = 4;

            v0 = VLEV_FLOAT(aoffset1, vl);
            VSEV_FLOAT(boffset2, v0, vl);

            aoffset1 += 4;
            //boffset2 += 4;
        }

        if (n & 1) {
           *(boffset3) = *(aoffset1);
           *(boffset3 + 1) = *(aoffset1 + 1);
        }
    }

    return 0;
}
