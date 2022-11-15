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
#define VSETVL(n) vsetvl_e32m2(n)
#define FLOAT_V_T vfloat32m2_t
#define VLSEG2_FLOAT vlseg2e32_v_f32m2
#define VSSSEG2_FLOAT vssseg2e32_v_f32m2
#define VSSSEG4_FLOAT vssseg4e32_v_f32m2
#else
#define VSETVL(n) vsetvl_e64m2(n)
#define FLOAT_V_T vfloat64m2_t
#define VLSEG2_FLOAT vlseg2e64_v_f64m2
#define VSSSEG2_FLOAT vssseg2e64_v_f64m2
#define VSSSEG4_FLOAT vssseg4e64_v_f64m2
#endif

// Optimizes the implementation in ../generic/gemm_tcopy_2.c

int CNAME(BLASLONG m, BLASLONG n, IFLOAT *a, BLASLONG lda, IFLOAT *b)
{
    BLASLONG i, j;
    IFLOAT *a_offset, *a_offset1, *a_offset2;
    IFLOAT *b_offset, *b_offset1, *b_offset2;
    FLOAT_V_T v1a, v1b, v2a, v2b;
    size_t vl;

    //fprintf(stderr, "gemm_tcopy_2 m=%ld n=%ld lda=%ld\n", m, n, lda); // KU

    a_offset = a;
    b_offset = b;
    b_offset2 = b + m * (n & ~1);

    for(i = (m >> 1); i > 0; i--) {

        a_offset1 = a_offset;
        a_offset2 = a_offset + lda;
        a_offset += 2 * lda;

        b_offset1 = b_offset;
        b_offset += 4;

        for(j = (n >> 1); j > 0; j -= vl) {
            vl = VSETVL(j);

            VLSEG2_FLOAT(&v1a, &v1b, a_offset1, vl);
            VLSEG2_FLOAT(&v2a, &v2b, a_offset2, vl);

            VSSSEG4_FLOAT(b_offset1, m*2*sizeof(FLOAT), v1a, v1b, v2a, v2b, vl);

            a_offset1 += vl * 2;
	        a_offset2 += vl * 2;
            b_offset1 += vl * m * 2;
        }

        if (n & 1) {
	        *(b_offset2 + 0) = *(a_offset1 + 0);
	        *(b_offset2 + 1) = *(a_offset2 + 0);
	        b_offset2 += 2;
        }
    }

    if (m & 1) {

        for(j = (n >> 1); j > 0; j -= vl) {
            vl = VSETVL(j);

            VLSEG2_FLOAT(&v1a, &v1b, a_offset, vl);

            VSSSEG2_FLOAT(b_offset, m*2*sizeof(FLOAT), v1a, v1b, vl);

            a_offset += vl * 2;
            b_offset += vl * m * 2;
        }

        if (n & 1){
            *(b_offset2 + 0) = *(a_offset + 0);
        }
    }

    return 0;
}
