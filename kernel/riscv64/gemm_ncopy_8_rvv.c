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
#define VSETVL(n) vsetvl_e32m1(n)
#define FLOAT_V_T vfloat32m1_t
#define VLEV_FLOAT vle32_v_f32m1
#define VSEV_FLOAT vse32_v_f32m1
#define VSSEG2_FLOAT vsseg2e32_v_f32m1
#define VSSEG4_FLOAT vsseg4e32_v_f32m1
#define VSSEG8_FLOAT vsseg8e32_v_f32m1
#else
#define VSETVL(n) vsetvl_e64m1(n)
#define FLOAT_V_T vfloat64m1_t
#define VLEV_FLOAT vle64_v_f64m1
#define VSEV_FLOAT vse64_v_f64m1
#define VSSEG2_FLOAT vsseg2e64_v_f64m1
#define VSSEG4_FLOAT vsseg4e64_v_f64m1
#define VSSEG8_FLOAT vsseg8e64_v_f64m1
#endif

// Optimizes the implementation in ../generic/gemm_ncopy_8.c

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, FLOAT *b)
{
    BLASLONG i, j;

    FLOAT *a_offset;
    FLOAT *a_offset1, *a_offset2, *a_offset3, *a_offset4;
    FLOAT *a_offset5, *a_offset6, *a_offset7, *a_offset8;
    FLOAT *b_offset;

    FLOAT_V_T v1, v2, v3, v4, v5, v6, v7, v8;
    size_t vl;

    //fprintf(stderr, "gemm_ncopy_8 m=%ld n=%ld lda=%ld\n", m, n, lda);

    a_offset = a;
    b_offset = b;

    for(j = (n >> 3); j > 0; j--) {
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

            VSSEG8_FLOAT(b_offset, v1, v2, v3, v4, v5, v6, v7, v8, vl);

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

            VSSEG4_FLOAT(b_offset, v1, v2, v3, v4, vl);

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

            VSSEG2_FLOAT(b_offset, v1, v2, vl);

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
