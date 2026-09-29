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
#define FLOAT_VX2_T             vfloat32m1x2_t
#define FLOAT_VX4_T             vfloat32m1x4_t
#define FLOAT_VX8_T             vfloat32m1x8_t
#define VGET_VX2                __riscv_vget_v_f32m1x2_f32m1
#define VSET_VX2                __riscv_vset_v_f32m1_f32m1x2
#define VSET_VX4                __riscv_vset_v_f32m1_f32m1x4
#define VSET_VX8                __riscv_vset_v_f32m1_f32m1x8
#define VLSEG2_FLOAT            __riscv_vlseg2e32_v_f32m1x2
#define VSSEG2_FLOAT            __riscv_vsseg2e32_v_f32m1x2
#define VSSEG4_FLOAT            __riscv_vsseg4e32_v_f32m1x4
#define VSSEG8_FLOAT            __riscv_vsseg8e32_v_f32m1x8
#define VSSSEG8_FLOAT           __riscv_vssseg8e32_v_f32m1x8
#else
#define VSETVL(n)               __riscv_vsetvl_e64m1(n)
#define FLOAT_VX2_T             vfloat64m1x2_t
#define FLOAT_VX4_T             vfloat64m1x4_t
#define FLOAT_VX8_T             vfloat64m1x8_t
#define VGET_VX2                __riscv_vget_v_f64m1x2_f64m1
#define VSET_VX2                __riscv_vset_v_f64m1_f64m1x2
#define VSET_VX4                __riscv_vset_v_f64m1_f64m1x4
#define VSET_VX8                __riscv_vset_v_f64m1_f64m1x8
#define VLSEG2_FLOAT            __riscv_vlseg2e64_v_f64m1x2
#define VSSEG2_FLOAT            __riscv_vsseg2e64_v_f64m1x2
#define VSSEG4_FLOAT            __riscv_vsseg4e64_v_f64m1x4
#define VSSEG8_FLOAT            __riscv_vsseg8e64_v_f64m1x8
#define VSSSEG8_FLOAT           __riscv_vssseg8e64_v_f64m1x8
#endif

/* Optimizes the implementation in ../generic/zgemm_ncopy_8.c
 *
 * The generic packer emits, for every complex element k of an 8-column panel,
 * the 16 floats  [r1,i1,r2,i2,r3,i3,r4,i4,r5,i5,r6,i6,r7,i7,r8,i8]  and relies
 * on the scalar compiler to build one base+offset address per store.  The
 * vector form below loads each of the eight source rows as one strided
 * segment-2 (real/imaginary) load and writes the panel with two strided
 * segment-8 stores whose byte stride is 16*sizeof(FLOAT).  This removes the
 * per-store address generation (the eight add/addi per iteration reported by
 * perf annotate) while producing the exact same packed byte sequence that
 * cgemm_kernel_* consumes.  The n&4 / n&2 / n&1 remainders keep the original
 * contiguous segment-8 / segment-4 / segment-2 layout.
 */

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, FLOAT *b){
    BLASLONG i, j;

    FLOAT *aoffset;
    FLOAT *aoffset1, *aoffset2, *aoffset3, *aoffset4;
    FLOAT *aoffset5, *aoffset6, *aoffset7, *aoffset8;
    FLOAT *boffset;

    FLOAT_VX2_T v1x2, v2x2, v3x2, v4x2, v5x2, v6x2, v7x2, v8x2;
    FLOAT_VX4_T vxx4;
    FLOAT_VX8_T vxx8a, vxx8b;

    ptrdiff_t pstride = (ptrdiff_t)(16 * sizeof(FLOAT));

    size_t vl;

    aoffset = a;
    boffset = b;
    lda *= 2;

    for (j = (n >> 3); j > 0; j--) {
        aoffset1  = aoffset;
        aoffset2  = aoffset1 + lda;
        aoffset3  = aoffset2 + lda;
        aoffset4  = aoffset3 + lda;
        aoffset5  = aoffset4 + lda;
        aoffset6  = aoffset5 + lda;
        aoffset7  = aoffset6 + lda;
        aoffset8  = aoffset7 + lda;
        aoffset  += 8 * lda;

        for (i = m; i > 0; i -= vl) {
            vl = VSETVL(i);

            v1x2 = VLSEG2_FLOAT(aoffset1, vl);
            v2x2 = VLSEG2_FLOAT(aoffset2, vl);
            v3x2 = VLSEG2_FLOAT(aoffset3, vl);
            v4x2 = VLSEG2_FLOAT(aoffset4, vl);
            v5x2 = VLSEG2_FLOAT(aoffset5, vl);
            v6x2 = VLSEG2_FLOAT(aoffset6, vl);
            v7x2 = VLSEG2_FLOAT(aoffset7, vl);
            v8x2 = VLSEG2_FLOAT(aoffset8, vl);

            vxx8a = VSET_VX8(vxx8a, 0, VGET_VX2(v1x2, 0));
            vxx8a = VSET_VX8(vxx8a, 1, VGET_VX2(v1x2, 1));
            vxx8a = VSET_VX8(vxx8a, 2, VGET_VX2(v2x2, 0));
            vxx8a = VSET_VX8(vxx8a, 3, VGET_VX2(v2x2, 1));
            vxx8a = VSET_VX8(vxx8a, 4, VGET_VX2(v3x2, 0));
            vxx8a = VSET_VX8(vxx8a, 5, VGET_VX2(v3x2, 1));
            vxx8a = VSET_VX8(vxx8a, 6, VGET_VX2(v4x2, 0));
            vxx8a = VSET_VX8(vxx8a, 7, VGET_VX2(v4x2, 1));

            vxx8b = VSET_VX8(vxx8b, 0, VGET_VX2(v5x2, 0));
            vxx8b = VSET_VX8(vxx8b, 1, VGET_VX2(v5x2, 1));
            vxx8b = VSET_VX8(vxx8b, 2, VGET_VX2(v6x2, 0));
            vxx8b = VSET_VX8(vxx8b, 3, VGET_VX2(v6x2, 1));
            vxx8b = VSET_VX8(vxx8b, 4, VGET_VX2(v7x2, 0));
            vxx8b = VSET_VX8(vxx8b, 5, VGET_VX2(v7x2, 1));
            vxx8b = VSET_VX8(vxx8b, 6, VGET_VX2(v8x2, 0));
            vxx8b = VSET_VX8(vxx8b, 7, VGET_VX2(v8x2, 1));

            VSSSEG8_FLOAT(boffset, pstride, vxx8a, vl);
            VSSSEG8_FLOAT(boffset + 8, pstride, vxx8b, vl);

            aoffset1 += vl * 2;
            aoffset2 += vl * 2;
            aoffset3 += vl * 2;
            aoffset4 += vl * 2;
            aoffset5 += vl * 2;
            aoffset6 += vl * 2;
            aoffset7 += vl * 2;
            aoffset8 += vl * 2;
            boffset  += vl * 16;
        }
    }

    if (n & 4) {
        aoffset1  = aoffset;
        aoffset2  = aoffset1 + lda;
        aoffset3  = aoffset2 + lda;
        aoffset4  = aoffset3 + lda;
        aoffset  += 4 * lda;

        for (i = m; i > 0; i -= vl) {
            vl = VSETVL(i);

            v1x2 = VLSEG2_FLOAT(aoffset1, vl);
            v2x2 = VLSEG2_FLOAT(aoffset2, vl);
            v3x2 = VLSEG2_FLOAT(aoffset3, vl);
            v4x2 = VLSEG2_FLOAT(aoffset4, vl);

            vxx8a = VSET_VX8(vxx8a, 0, VGET_VX2(v1x2, 0));
            vxx8a = VSET_VX8(vxx8a, 1, VGET_VX2(v1x2, 1));
            vxx8a = VSET_VX8(vxx8a, 2, VGET_VX2(v2x2, 0));
            vxx8a = VSET_VX8(vxx8a, 3, VGET_VX2(v2x2, 1));
            vxx8a = VSET_VX8(vxx8a, 4, VGET_VX2(v3x2, 0));
            vxx8a = VSET_VX8(vxx8a, 5, VGET_VX2(v3x2, 1));
            vxx8a = VSET_VX8(vxx8a, 6, VGET_VX2(v4x2, 0));
            vxx8a = VSET_VX8(vxx8a, 7, VGET_VX2(v4x2, 1));

            VSSEG8_FLOAT(boffset, vxx8a, vl);

            aoffset1 += vl * 2;
            aoffset2 += vl * 2;
            aoffset3 += vl * 2;
            aoffset4 += vl * 2;
            boffset  += vl * 8;
        }
    }

    if (n & 2) {
        aoffset1  = aoffset;
        aoffset2  = aoffset1 + lda;
        aoffset  += 2 * lda;

        for (i = m; i > 0; i -= vl) {
            vl = VSETVL(i);

            v1x2 = VLSEG2_FLOAT(aoffset1, vl);
            v2x2 = VLSEG2_FLOAT(aoffset2, vl);

            vxx4 = VSET_VX4(vxx4, 0, VGET_VX2(v1x2, 0));
            vxx4 = VSET_VX4(vxx4, 1, VGET_VX2(v1x2, 1));
            vxx4 = VSET_VX4(vxx4, 2, VGET_VX2(v2x2, 0));
            vxx4 = VSET_VX4(vxx4, 3, VGET_VX2(v2x2, 1));

            VSSEG4_FLOAT(boffset, vxx4, vl);

            aoffset1 += vl * 2;
            aoffset2 += vl * 2;
            boffset  += vl * 4;
        }
    }

    if (n & 1) {
        aoffset1  = aoffset;
        aoffset  += lda;

        for (i = m; i > 0; i -= vl) {
            vl = VSETVL(i);

            v1x2 = VLSEG2_FLOAT(aoffset1, vl);

            VSSEG2_FLOAT(boffset, v1x2, vl);

            aoffset1 += vl * 2;
            boffset  += vl * 2;
        }
    }

    return 0;
}