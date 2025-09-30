/***************************************************************************
Copyright (c) 2025, The OpenBLAS Project
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
#define VLSEG_FLOAT             __riscv_vlse32_v_f32m1
#define VLSSEG2_FLOAT           __riscv_vlsseg2e32_v_f32m1x2
#define VLSSEG4_FLOAT           __riscv_vlsseg4e32_v_f32m1x4
#define VLSSEG8_FLOAT           __riscv_vlsseg8e32_v_f32m1x8
#define VGET_VX2                __riscv_vget_v_f32m1x2_f32m1
#define VGET_VX4                __riscv_vget_v_f32m1x4_f32m1
#define VGET_VX8                __riscv_vget_v_f32m1x8_f32m1
#define VSET_VX2                __riscv_vset_v_f32m1_f32m1x2
#define VSET_VX4                __riscv_vset_v_f32m1_f32m1x4
#define VSET_VX8                __riscv_vset_v_f32m1_f32m1x8
#define VLEV_FLOAT              __riscv_vle32_v_f32m1
#define VSEV_FLOAT              __riscv_vse32_v_f32m1
#define VSSEG2_FLOAT            __riscv_vsseg2e32_v_f32m1x2
#define VSSEG4_FLOAT            __riscv_vsseg4e32_v_f32m1x4
#define VSSEG8_FLOAT            __riscv_vsseg8e32_v_f32m1x8
#else
#define VSETVL(n)               __riscv_vsetvl_e64m1(n)
#define FLOAT_V_T               vfloat64m1_t
#define FLOAT_VX2_T             vfloat64m1x2_t
#define FLOAT_VX4_T             vfloat64m1x4_t
#define FLOAT_VX8_T             vfloat64m1x8_t
#define VLSEG_FLOAT             __riscv_vlse64_v_f64m1
#define VLSSEG2_FLOAT           __riscv_vlsseg2e64_v_f64m1x2
#define VLSSEG4_FLOAT           __riscv_vlsseg4e64_v_f64m1x4
#define VLSSEG8_FLOAT           __riscv_vlsseg8e64_v_f64m1x8
#define VGET_VX2                __riscv_vget_v_f64m1x2_f64m1
#define VGET_VX4                __riscv_vget_v_f64m1x4_f64m1
#define VGET_VX8                __riscv_vget_v_f64m1x8_f64m1
#define VSET_VX2                __riscv_vset_v_f64m1_f64m1x2
#define VSET_VX4                __riscv_vset_v_f64m1_f64m1x4
#define VSET_VX8                __riscv_vset_v_f64m1_f64m1x8
#define VLEV_FLOAT              __riscv_vle64_v_f64m1
#define VSEV_FLOAT              __riscv_vse64_v_f64m1
#define VSSEG2_FLOAT            __riscv_vsseg2e64_v_f64m1x2
#define VSSEG4_FLOAT            __riscv_vsseg4e64_v_f64m1x4
#define VSSEG8_FLOAT            __riscv_vsseg8e64_v_f64m1x8
#endif

// Optimizes the implementation in ../generic/gemm_ncopy_16.c

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, FLOAT *b)
{
    BLASLONG i, j;

    FLOAT *a_offset;
    FLOAT *a_offset1, *a_offset2, *a_offset3, *a_offset4;
    FLOAT *a_offset5, *a_offset6, *a_offset7, *a_offset8;
    FLOAT *b_offset;

    FLOAT_V_T v1, v2, v3, v4, v5, v6, v7, v8;
    FLOAT_V_T v9, v10, v11, v12, v13, v14, v15, v16;
    FLOAT_VX2_T vx2, vx21;
    FLOAT_VX4_T vx4, vx41;
    FLOAT_VX8_T vx8, vx81;

    size_t vl;

    //fprintf(stderr, "gemm_ncopy_16 m=%ld n=%ld lda=%ld\n", m, n, lda);

    a_offset = a;
    b_offset = b;

    j = (n >> 4);
    if (j) {
        vl = VSETVL(8);

        do {
            a_offset1  = a_offset;
            a_offset2  = a_offset1  + lda * 8;
            a_offset += 16 * lda;

            i = m >> 3;
            if (i) {
                do {
                    vx8 = VLSSEG8_FLOAT(a_offset1, lda * sizeof(FLOAT), vl);
                    vx81 = VLSSEG8_FLOAT(a_offset2, lda * sizeof(FLOAT), vl);

                    v1 = VGET_VX8(vx8, 0);
                    v2 = VGET_VX8(vx8, 1);
                    v3 = VGET_VX8(vx8, 2);
                    v4 = VGET_VX8(vx8, 3);
                    v5 = VGET_VX8(vx8, 4);
                    v6 = VGET_VX8(vx8, 5);
                    v7 = VGET_VX8(vx8, 6);
                    v8 = VGET_VX8(vx8, 7);
                    v9 = VGET_VX8(vx81, 0);
                    v10 = VGET_VX8(vx81, 1);
                    v11 = VGET_VX8(vx81, 2);
                    v12 = VGET_VX8(vx81, 3);
                    v13 = VGET_VX8(vx81, 4);
                    v14 = VGET_VX8(vx81, 5);
                    v15 = VGET_VX8(vx81, 6);
                    v16 = VGET_VX8(vx81, 7);

                    VSEV_FLOAT(b_offset, v1, vl);
                    VSEV_FLOAT(b_offset + 8, v9, vl);
                    VSEV_FLOAT(b_offset + 16, v2, vl);
                    VSEV_FLOAT(b_offset + 24, v10, vl);
                    VSEV_FLOAT(b_offset + 32, v3, vl);
                    VSEV_FLOAT(b_offset + 40, v11, vl);
                    VSEV_FLOAT(b_offset + 48, v4, vl);
                    VSEV_FLOAT(b_offset + 56, v12, vl);
                    VSEV_FLOAT(b_offset + 64, v5, vl);
                    VSEV_FLOAT(b_offset + 72, v13, vl);
                    VSEV_FLOAT(b_offset + 80, v6, vl);
                    VSEV_FLOAT(b_offset + 88, v14, vl);
                    VSEV_FLOAT(b_offset + 96, v7, vl);
                    VSEV_FLOAT(b_offset + 104, v15, vl);
                    VSEV_FLOAT(b_offset + 112, v8, vl);
                    VSEV_FLOAT(b_offset + 120, v16, vl);

                    a_offset1 += 8;
                    a_offset2 += 8;
                    b_offset += 128;
                } while (--i);
            }

            if (m & 4) {
                vx4 = VLSSEG4_FLOAT(a_offset1, lda * sizeof(FLOAT), vl);
                vx41 = VLSSEG4_FLOAT(a_offset2, lda * sizeof(FLOAT), vl);

                v1 = VGET_VX4(vx4, 0);
                v2 = VGET_VX4(vx4, 1);
                v3 = VGET_VX4(vx4, 2);
                v4 = VGET_VX4(vx4, 3);
                v5 = VGET_VX4(vx41, 0);
                v6 = VGET_VX4(vx41, 1);
                v7 = VGET_VX4(vx41, 2);
                v8 = VGET_VX4(vx41, 3);

                VSEV_FLOAT(b_offset, v1, vl);
                VSEV_FLOAT(b_offset + 8, v5, vl);
                VSEV_FLOAT(b_offset + 16, v2, vl);
                VSEV_FLOAT(b_offset + 24, v6, vl);
                VSEV_FLOAT(b_offset + 32, v3, vl);
                VSEV_FLOAT(b_offset + 40, v7, vl);
                VSEV_FLOAT(b_offset + 48, v4, vl);
                VSEV_FLOAT(b_offset + 56, v8, vl);

                a_offset1 += 4;
                a_offset2 += 4;
                b_offset += 64;
            }

            if (m & 2) {
                vx2 = VLSSEG2_FLOAT(a_offset1, lda * sizeof(FLOAT), vl);
                vx21 = VLSSEG2_FLOAT(a_offset2, lda * sizeof(FLOAT), vl);

                v1 = VGET_VX2(vx2, 0);
                v2 = VGET_VX2(vx2, 1);
                v3 = VGET_VX2(vx21, 0);
                v4 = VGET_VX2(vx21, 1);

                VSEV_FLOAT(b_offset, v1, vl);
                VSEV_FLOAT(b_offset + 8, v3, vl);
                VSEV_FLOAT(b_offset + 16, v2, vl);
                VSEV_FLOAT(b_offset + 24, v4, vl);

                a_offset1 += 2;
                a_offset2 += 2;
                b_offset += 32;
            }

            if (m & 1) {
                v1 = VLSEG_FLOAT(a_offset1, lda * sizeof(FLOAT), vl);
                v2 = VLSEG_FLOAT(a_offset2, lda * sizeof(FLOAT), vl);

                VSEV_FLOAT(b_offset, v1, vl);
                VSEV_FLOAT(b_offset + 8, v2, vl);

                b_offset += 16;
            }
        } while (--j);
    }

    if (n & 8) {
        a_offset1  = a_offset;
        a_offset2  = a_offset1 + lda;
        a_offset3  = a_offset2 + lda;
        a_offset4  = a_offset3 + lda;
        a_offset5  = a_offset4 + lda;
        a_offset6  = a_offset5 + lda;
        a_offset7  = a_offset6 + lda;
        a_offset8  = a_offset7 + lda;
        a_offset += 8 * lda;

        for(i = m; i > 0; i -= vl) {
            vl = VSETVL(i);

            v1 = VLEV_FLOAT(a_offset1, vl);
            v2 = VLEV_FLOAT(a_offset2, vl);
            v3 = VLEV_FLOAT(a_offset3, vl);
            v4 = VLEV_FLOAT(a_offset4, vl);
            v5 = VLEV_FLOAT(a_offset5, vl);
            v6 = VLEV_FLOAT(a_offset6, vl);
            v7 = VLEV_FLOAT(a_offset7, vl);
            v8 = VLEV_FLOAT(a_offset8, vl);

            vx8 = VSET_VX8(vx8, 0, v1);
            vx8 = VSET_VX8(vx8, 1, v2);
            vx8 = VSET_VX8(vx8, 2, v3);
            vx8 = VSET_VX8(vx8, 3, v4);
            vx8 = VSET_VX8(vx8, 4, v5);
            vx8 = VSET_VX8(vx8, 5, v6);
            vx8 = VSET_VX8(vx8, 6, v7);
            vx8 = VSET_VX8(vx8, 7, v8);

            VSSEG8_FLOAT(b_offset, vx8, vl);

            a_offset1 += vl;
            a_offset2 += vl;
            a_offset3 += vl;
            a_offset4 += vl;
            a_offset5 += vl;
            a_offset6 += vl;
            a_offset7 += vl;
            a_offset8 += vl;
            b_offset += vl*8;
        }
    }

    if (n & 4) {
        a_offset1  = a_offset;
        a_offset2  = a_offset1 + lda;
        a_offset3  = a_offset2 + lda;
        a_offset4  = a_offset3 + lda;
        a_offset += 4 * lda;

        for(i = m; i > 0; i -= vl) {
            vl = VSETVL(i);

            v1 = VLEV_FLOAT(a_offset1, vl);
            v2 = VLEV_FLOAT(a_offset2, vl);
            v3 = VLEV_FLOAT(a_offset3, vl);
            v4 = VLEV_FLOAT(a_offset4, vl);

            vx4 = VSET_VX4(vx4, 0, v1);
            vx4 = VSET_VX4(vx4, 1, v2);
            vx4 = VSET_VX4(vx4, 2, v3);
            vx4 = VSET_VX4(vx4, 3, v4);

            VSSEG4_FLOAT(b_offset, vx4, vl);

            a_offset1 += vl;
            a_offset2 += vl;
            a_offset3 += vl;
            a_offset4 += vl;
            b_offset += vl*4;
        }
    }

    if (n & 2) {
        a_offset1  = a_offset;
        a_offset2  = a_offset1 + lda;
        a_offset += 2 * lda;

        for(i = m; i > 0; i -= vl) {
            vl = VSETVL(i);

            v1 = VLEV_FLOAT(a_offset1, vl);
            v2 = VLEV_FLOAT(a_offset2, vl);

            vx2 = VSET_VX2(vx2, 0, v1);
            vx2 = VSET_VX2(vx2, 1, v2);

            VSSEG2_FLOAT(b_offset, vx2, vl);

            a_offset1 += vl;
            a_offset2 += vl;
            b_offset += vl*2;
        }
    }

    if (n & 1) {
        a_offset1  = a_offset;

        for(i = m; i > 0; i -= vl) {
            vl = VSETVL(i);

            v1 = VLEV_FLOAT(a_offset1, vl);

            VSEV_FLOAT(b_offset, v1, vl);

            a_offset1 += vl;
            b_offset += vl;
        }
    }

    return 0;
}
