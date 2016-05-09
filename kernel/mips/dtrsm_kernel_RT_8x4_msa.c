/*******************************************************************************
Copyright (c) 2016, The OpenBLAS Project
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
*******************************************************************************/

#include "common.h"
#include "macros_msa.h"

static void dsolve_8x4_rt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG bk)
{
    v2f64 src_c0, src_c1, src_c2, src_c3, src_c4, src_c5, src_c6, src_c7;
    v2f64 src_c8, src_c9, src_c10, src_c11, src_c12, src_c13, src_c14, src_c15;
    v2f64 src_b0, src_b4, src_b5, src_b8, src_b9, src_b10, src_b12, src_b13;
    v2f64 src_b14, src_b15;
    FLOAT *c_nxt1line = c + ldc;
    FLOAT *c_nxt2line = c + 2 * ldc;
    FLOAT *c_nxt3line = c + 3 * ldc;

    LD_DP4(c, 2, src_c0, src_c1, src_c2, src_c3);
    LD_DP4(c_nxt1line, 2, src_c4, src_c5, src_c6, src_c7);
    LD_DP4(c_nxt2line, 2, src_c8, src_c9, src_c10, src_c11);
    LD_DP4(c_nxt3line, 2, src_c12, src_c13, src_c14, src_c15);

    if (bk > 0)
    {
        BLASLONG i;
        FLOAT *pba = a, *pbb = b;
        v2f64 src_b, src_b0, src_b1, src_b2, src_b3;
        v2f64 src_a0, src_a1, src_a2, src_a3, src_a4, src_a5, src_a6, src_a7;

        LD_DP4(pba, 2, src_a0, src_a1, src_a2, src_a3);
        LD_DP2(pbb, 2, src_b0, src_b1);

        for (i = (bk - 1); i--;)
        {
            pba += 8;
            pbb += 4;

            LD_DP4(pba, 2, src_a4, src_a5, src_a6, src_a7);
            LD_DP2(pbb, 2, src_b2, src_b3);

            src_b = (v2f64) __msa_ilvr_d((v2i64) src_b0, (v2i64) src_b0);
            src_c0 -= src_a0 * src_b;
            src_c1 -= src_a1 * src_b;
            src_c2 -= src_a2 * src_b;
            src_c3 -= src_a3 * src_b;

            src_b = (v2f64) __msa_ilvl_d((v2i64) src_b0, (v2i64) src_b0);
            src_c4 -= src_a0 * src_b;
            src_c5 -= src_a1 * src_b;
            src_c6 -= src_a2 * src_b;
            src_c7 -= src_a3 * src_b;

            src_b = (v2f64) __msa_ilvr_d((v2i64) src_b1, (v2i64) src_b1);
            src_c8  -= src_a0 * src_b;
            src_c9  -= src_a1 * src_b;
            src_c10 -= src_a2 * src_b;
            src_c11 -= src_a3 * src_b;

            src_b = (v2f64) __msa_ilvl_d((v2i64) src_b1, (v2i64) src_b1);
            src_c12 -= src_a0 * src_b;
            src_c13 -= src_a1 * src_b;
            src_c14 -= src_a2 * src_b;
            src_c15 -= src_a3 * src_b;

            src_a0 = src_a4;
            src_a1 = src_a5;
            src_a2 = src_a6;
            src_a3 = src_a7;
            src_b0 = src_b2;
            src_b1 = src_b3;
        }

        src_b = (v2f64) __msa_ilvr_d((v2i64) src_b0, (v2i64) src_b0);
        src_c0 -= src_a0 * src_b;
        src_c1 -= src_a1 * src_b;
        src_c2 -= src_a2 * src_b;
        src_c3 -= src_a3 * src_b;

        src_b = (v2f64) __msa_ilvl_d((v2i64) src_b0, (v2i64) src_b0);
        src_c4 -= src_a0 * src_b;
        src_c5 -= src_a1 * src_b;
        src_c6 -= src_a2 * src_b;
        src_c7 -= src_a3 * src_b;

        src_b = (v2f64) __msa_ilvr_d((v2i64) src_b1, (v2i64) src_b1);
        src_c8  -= src_a0 * src_b;
        src_c9  -= src_a1 * src_b;
        src_c10 -= src_a2 * src_b;
        src_c11 -= src_a3 * src_b;

        src_b = (v2f64) __msa_ilvl_d((v2i64) src_b1, (v2i64) src_b1);
        src_c12 -= src_a0 * src_b;
        src_c13 -= src_a1 * src_b;
        src_c14 -= src_a2 * src_b;
        src_c15 -= src_a3 * src_b;
    }

    a -= 32;
    b -= 16;

    src_b12 = LD_DP(b + 12);
    src_b13 = (v2f64) __msa_splati_d((v2i64) src_b12, 1);
    src_b12 = (v2f64) __msa_splati_d((v2i64) src_b12, 0);
    src_b14 = LD_DP(b + 14);
    src_b15 = (v2f64) __msa_splati_d((v2i64) src_b14, 1);
    src_b14 = (v2f64) __msa_splati_d((v2i64) src_b14, 0);

    src_b8 = LD_DP(b + 8);
    src_b9 = (v2f64) __msa_splati_d((v2i64) src_b8, 1);
    src_b8 = (v2f64) __msa_splati_d((v2i64) src_b8, 0);
    src_b10 = __msa_cast_to_vector_double(*(b + 10));
    src_b10 = (v2f64) __msa_splati_d((v2i64) src_b10, 0);

    src_b0 = __msa_cast_to_vector_double(*(b + 0));
    src_b0 = (v2f64) __msa_splati_d((v2i64) src_b0, 0);
    src_b4 = LD_DP(b + 4);
    src_b5 = (v2f64) __msa_splati_d((v2i64) src_b4, 1);
    src_b4 = (v2f64) __msa_splati_d((v2i64) src_b4, 0);

    src_c12 *= src_b15;
    src_c13 *= src_b15;
    src_c14 *= src_b15;
    src_c15 *= src_b15;

    src_c8 -= src_c12 * src_b14;
    src_c9 -= src_c13 * src_b14;
    src_c10 -= src_c14 * src_b14;
    src_c11 -= src_c15 * src_b14;

    src_c8 *= src_b10;
    src_c9 *= src_b10;
    src_c10 *= src_b10;
    src_c11 *= src_b10;

    src_c4 -= src_c12 * src_b13;
    src_c5 -= src_c13 * src_b13;
    src_c6 -= src_c14 * src_b13;
    src_c7 -= src_c15 * src_b13;

    src_c4 -= src_c8 * src_b9;
    src_c5 -= src_c9 * src_b9;
    src_c6 -= src_c10 * src_b9;
    src_c7 -= src_c11 * src_b9;

    src_c4 *= src_b5;
    src_c5 *= src_b5;
    src_c6 *= src_b5;
    src_c7 *= src_b5;

    src_c0 -= src_c12 * src_b12;
    src_c1 -= src_c13 * src_b12;
    src_c2 -= src_c14 * src_b12;
    src_c3 -= src_c15 * src_b12;

    src_c0 -= src_c8 * src_b8;
    src_c1 -= src_c9 * src_b8;
    src_c2 -= src_c10 * src_b8;
    src_c3 -= src_c11 * src_b8;

    src_c0 -= src_c4 * src_b4;
    src_c1 -= src_c5 * src_b4;
    src_c2 -= src_c6 * src_b4;
    src_c3 -= src_c7 * src_b4;

    src_c0 *= src_b0;
    src_c1 *= src_b0;
    src_c2 *= src_b0;
    src_c3 *= src_b0;

    ST_DP4(src_c12, src_c13, src_c14, src_c15, c_nxt3line, 2);
    ST_DP4(src_c12, src_c13, src_c14, src_c15, a + 24, 2);
    ST_DP4(src_c8, src_c9, src_c10, src_c11, c_nxt2line, 2);
    ST_DP4(src_c8, src_c9, src_c10, src_c11, a + 16, 2);
    ST_DP4(src_c4, src_c5, src_c6, src_c7, c_nxt1line, 2);
    ST_DP4(src_c4, src_c5, src_c6, src_c7, a + 8, 2);
    ST_DP4(src_c0, src_c1, src_c2, src_c3, c, 2);
    ST_DP4(src_c0, src_c1, src_c2, src_c3, a, 2);
}

static void dsolve_8x2_rt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, int bk)
{
    v2f64 src_c0, src_c1, src_c2, src_c3, src_c4, src_c5, src_c6, src_c7;
    v2f64 src_b0, src_b2, src_b3;

    LD_DP4(c, 2, src_c0, src_c1, src_c2, src_c3);
    LD_DP4(c + ldc, 2, src_c4, src_c5, src_c6, src_c7);

    if (bk > 0)
    {
        v2f64 src_a0, src_a1, src_a2, src_a3, src_b, src_b0;

        LD_DP4(a + 16, 2, src_a0, src_a1, src_a2, src_a3);
        src_b0 = LD_DP(b + 4);

        src_b = (v2f64) __msa_ilvr_d((v2i64) src_b0, (v2i64) src_b0);
        src_c0 -= src_a0 * src_b;
        src_c1 -= src_a1 * src_b;
        src_c2 -= src_a2 * src_b;
        src_c3 -= src_a3 * src_b;

        src_b = (v2f64) __msa_ilvl_d((v2i64) src_b0, (v2i64) src_b0);
        src_c4 -= src_a0 * src_b;
        src_c5 -= src_a1 * src_b;
        src_c6 -= src_a2 * src_b;
        src_c7 -= src_a3 * src_b;
    }

    src_b0 = __msa_cast_to_vector_double(*(b + 0));
    src_b0 = (v2f64) __msa_splati_d((v2i64) src_b0, 0);
    src_b2 = LD_DP(b + 2);
    src_b3 = (v2f64) __msa_splati_d((v2i64) src_b2, 1);
    src_b2 = (v2f64) __msa_splati_d((v2i64) src_b2, 0);

    src_c4 *= src_b3;
    src_c5 *= src_b3;
    src_c6 *= src_b3;
    src_c7 *= src_b3;

    src_c0 -= src_c4 * src_b2;
    src_c1 -= src_c5 * src_b2;
    src_c2 -= src_c6 * src_b2;
    src_c3 -= src_c7 * src_b2;

    src_c0 *= src_b0;
    src_c1 *= src_b0;
    src_c2 *= src_b0;
    src_c3 *= src_b0;

    ST_DP4(src_c0, src_c1, src_c2, src_c3, c, 2);
    ST_DP4(src_c4, src_c5, src_c6, src_c7, c + ldc, 2);

    ST_DP4(src_c0, src_c1, src_c2, src_c3, a, 2);
    ST_DP4(src_c4, src_c5, src_c6, src_c7, a + 8, 2);
}

static void dsolve_8x1_rt_msa(FLOAT *a, FLOAT *b, FLOAT *c)
{
    v2f64 src_c0, src_c1, src_c2, src_c3;
    v2f64 src_b0;

    LD_DP4(c, 2, src_c0, src_c1, src_c2, src_c3);

    src_b0 = __msa_cast_to_vector_double(*b);
    src_b0 = (v2f64) __msa_splati_d((v2i64) src_b0, 0);

    src_c0 *= src_b0;
    src_c1 *= src_b0;
    src_c2 *= src_b0;
    src_c3 *= src_b0;

    ST_DP4(src_c0, src_c1, src_c2, src_c3, c, 2);
    ST_DP4(src_c0, src_c1, src_c2, src_c3, a, 2);
}

static void dsolve_4x4_rt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG bk)
{
    v2f64 src_c0, src_c1, src_c2, src_c3, src_c4, src_c5, src_c6, src_c7;
    v2f64 src_b0, src_b4, src_b5, src_b8, src_b9, src_b10, src_b12, src_b13;
    v2f64 src_b14, src_b15;

    LD_DP2(c, 2, src_c0, src_c1);
    LD_DP2(c + ldc, 2, src_c2, src_c3);
    LD_DP2(c + 2 * ldc, 2, src_c4, src_c5);
    LD_DP2(c + 3 * ldc, 2, src_c6, src_c7);

    if (bk > 0)
    {
        BLASLONG i;
        FLOAT *aa = a + 16, *bb = b + 16;
        v2f64 src_a0, src_a1, src_b, src_b0, src_b1;

        for (i = bk; i--;)
        {
            LD_DP2(aa, 2, src_a0, src_a1);
            LD_DP2(bb, 2, src_b0, src_b1);

            src_b = (v2f64) __msa_ilvr_d((v2i64) src_b0, (v2i64) src_b0);
            src_c0 -= src_a0 * src_b;
            src_c1 -= src_a1 * src_b;

            src_b = (v2f64) __msa_ilvl_d((v2i64) src_b0, (v2i64) src_b0);
            src_c2 -= src_a0 * src_b;
            src_c3 -= src_a1 * src_b;

            src_b = (v2f64) __msa_ilvr_d((v2i64) src_b1, (v2i64) src_b1);
            src_c4 -= src_a0 * src_b;
            src_c5 -= src_a1 * src_b;

            src_b = (v2f64) __msa_ilvl_d((v2i64) src_b1, (v2i64) src_b1);
            src_c6 -= src_a0 * src_b;
            src_c7 -= src_a1 * src_b;

            aa += 4;
            bb += 4;
        }
    }

    src_b12 = LD_DP(b + 12);
    src_b13 = (v2f64) __msa_splati_d((v2i64) src_b12, 1);
    src_b12 = (v2f64) __msa_splati_d((v2i64) src_b12, 0);
    src_b14 = LD_DP(b + 14);
    src_b15 = (v2f64) __msa_splati_d((v2i64) src_b14, 1);
    src_b14 = (v2f64) __msa_splati_d((v2i64) src_b14, 0);

    src_b8 = LD_DP(b + 8);
    src_b9 = (v2f64) __msa_splati_d((v2i64) src_b8, 1);
    src_b8 = (v2f64) __msa_splati_d((v2i64) src_b8, 0);
    src_b10 = __msa_cast_to_vector_double(*(b + 10));
    src_b10 = (v2f64) __msa_splati_d((v2i64) src_b10, 0);

    src_b0 = __msa_cast_to_vector_double(*(b + 0));
    src_b0 = (v2f64) __msa_splati_d((v2i64) src_b0, 0);
    src_b4 = LD_DP(b + 4);
    src_b5 = (v2f64) __msa_splati_d((v2i64) src_b4, 1);
    src_b4 = (v2f64) __msa_splati_d((v2i64) src_b4, 0);

    src_c6 *= src_b15;
    src_c7 *= src_b15;

    src_c4 -= src_c6 * src_b14;
    src_c5 -= src_c7 * src_b14;

    src_c4 *= src_b10;
    src_c5 *= src_b10;

    src_c2 -= src_c6 * src_b13;
    src_c3 -= src_c7 * src_b13;

    src_c2 -= src_c4 * src_b9;
    src_c3 -= src_c5 * src_b9;

    src_c2 *= src_b5;
    src_c3 *= src_b5;

    src_c0 -= src_c6 * src_b12;
    src_c1 -= src_c7 * src_b12;

    src_c0 -= src_c4 * src_b8;
    src_c1 -= src_c5 * src_b8;

    src_c0 -= src_c2 * src_b4;
    src_c1 -= src_c3 * src_b4;

    src_c0 *= src_b0;
    src_c1 *= src_b0;

    ST_DP2(src_c6, src_c7, c + 3 * ldc, 2);
    ST_DP2(src_c4, src_c5, c + 2 * ldc, 2);
    ST_DP2(src_c2, src_c3, c + ldc, 2);
    ST_DP2(src_c0, src_c1, c, 2);

    ST_DP4(src_c4, src_c5, src_c6, src_c7, a + 8, 2);
    ST_DP4(src_c0, src_c1, src_c2, src_c3, a, 2);
}

static void dsolve_4x2_rt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, int bk)
{
    v2f64 src_c0, src_c1, src_c2, src_c3, src_b0, src_b2, src_b3;

    LD_DP2(c, 2, src_c0, src_c1);
    LD_DP2(c + ldc, 2, src_c2, src_c3);

    if (bk > 0)
    {
        v2f64 src_a0, src_a1, src_b, src_b0;

        LD_DP2(a + 8, 2, src_a0, src_a1);
        src_b0 = LD_DP(b + 4);

        src_b = (v2f64) __msa_ilvr_d((v2i64) src_b0, (v2i64) src_b0);
        src_c0 -= src_a0 * src_b;
        src_c1 -= src_a1 * src_b;

        src_b = (v2f64) __msa_ilvl_d((v2i64) src_b0, (v2i64) src_b0);
        src_c2 -= src_a0 * src_b;
        src_c3 -= src_a1 * src_b;
    }

    src_b0 = __msa_cast_to_vector_double(*(b + 0));
    src_b0 = (v2f64) __msa_splati_d((v2i64) src_b0, 0);
    src_b2 = LD_DP(b + 2);
    src_b3 = (v2f64) __msa_splati_d((v2i64) src_b2, 1);
    src_b2 = (v2f64) __msa_splati_d((v2i64) src_b2, 0);

    src_c2 *= src_b3;
    src_c3 *= src_b3;

    src_c0 -= src_c2 * src_b2;
    src_c1 -= src_c3 * src_b2;

    src_c0 *= src_b0;
    src_c1 *= src_b0;

    ST_DP2(src_c0, src_c1, c, 2);
    ST_DP2(src_c2, src_c3, c + ldc, 2);

    ST_DP4(src_c0, src_c1, src_c2, src_c3, a, 2);
}

static void dsolve_4x1_rt_msa(FLOAT *a, FLOAT *b, FLOAT *c)
{
    FLOAT b0, c0, c1, c2, c3;

    b0 = *(b + 0);

    c0 = *(c + 0);
    c1 = *(c + 1);
    c2 = *(c + 2);
    c3 = *(c + 3);

    c0 *= b0;
    c1 *= b0;
    c2 *= b0;
    c3 *= b0;

    *(a + 0) = c0;
    *(a + 1) = c1;
    *(a + 2) = c2;
    *(a + 3) = c3;

    *(c + 0) = c0;
    *(c + 1) = c1;
    *(c + 2) = c2;
    *(c + 3) = c3;
}

static void dsolve_2x4_rt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG bk)
{
    FLOAT b0, b4, b5, b8, b9, b10, b12, b13, b14, b15;
    FLOAT c0, c1, c0_nxt1, c1_nxt1, c0_nxt2, c1_nxt2, c0_nxt3, c1_nxt3;

    c0 = *(c + 0);
    c1 = *(c + 1);
    c0_nxt1 = *(c + 0 + 1 * ldc);
    c1_nxt1 = *(c + 1 + 1 * ldc);
    c0_nxt2 = *(c + 0 + 2 * ldc);
    c1_nxt2 = *(c + 1 + 2 * ldc);
    c0_nxt3 = *(c + 0 + 3 * ldc);
    c1_nxt3 = *(c + 1 + 3 * ldc);

    if (bk > 0)
    {
        BLASLONG i;
        FLOAT *aa = a + 8, *bb = b + 16;
        FLOAT a0, a1, b0, b1, b2, b3;

        for (i = bk; i--;)
        {
            a0 = aa[0];
            a1 = aa[1];

            b0 = bb[0];
            c0 -= a0 * b0;
            c1 -= a1 * b0;

            b1 = bb[1];
            c0_nxt1 -= a0 * b1;
            c1_nxt1 -= a1 * b1;

            b2 = bb[2];
            c0_nxt2 -= a0 * b2;
            c1_nxt2 -= a1 * b2;

            b3 = bb[3];
            c0_nxt3 -= a0 * b3;
            c1_nxt3 -= a1 * b3;

            aa += 2;
            bb += 4;
        }
    }

    b0 = *b;
    b4 = *(b + 4);
    b5 = *(b + 5);
    b8 = *(b + 8);
    b9 = *(b + 9);
    b10 = *(b + 10);
    b12 = *(b + 12);
    b13 = *(b + 13);
    b14 = *(b + 14);
    b15 = *(b + 15);

    c0_nxt3 *= b15;
    c1_nxt3 *= b15;

    c0_nxt2 -= c0_nxt3 * b14;
    c1_nxt2 -= c1_nxt3 * b14;
    c0_nxt2 *= b10;
    c1_nxt2 *= b10;

    c0_nxt1 -= c0_nxt3 * b13;
    c1_nxt1 -= c1_nxt3 * b13;
    c0_nxt1 -= c0_nxt2 * b9;
    c1_nxt1 -= c1_nxt2 * b9;
    c0_nxt1 *= b5;
    c1_nxt1 *= b5;

    c0 -= c0_nxt3 * b12;
    c1 -= c1_nxt3 * b12;
    c0 -= c0_nxt2 * b8;
    c1 -= c1_nxt2 * b8;
    c0 -= c0_nxt1 * b4;
    c1 -= c1_nxt1 * b4;
    c0 *= b0;
    c1 *= b0;

    *(a + 0) = c0;
    *(a + 1) = c1;
    *(a + 2) = c0_nxt1;
    *(a + 3) = c1_nxt1;
    *(a + 4) = c0_nxt2;
    *(a + 5) = c1_nxt2;
    *(a + 6) = c0_nxt3;
    *(a + 7) = c1_nxt3;

    *(c + 0) = c0;
    *(c + 1) = c1;

    *(c + 0 + 1 * ldc) = c0_nxt1;
    *(c + 1 + 1 * ldc) = c1_nxt1;

    *(c + 0 + 2 * ldc) = c0_nxt2;
    *(c + 1 + 2 * ldc) = c1_nxt2;

    *(c + 0 + 3 * ldc) = c0_nxt3;
    *(c + 1 + 3 * ldc) = c1_nxt3;
}

static void dsolve_2x2_rt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, int bk)
{
    FLOAT b0, b2, b3;
    FLOAT c0, c1, c0_nxt, c1_nxt;

    c0 = *(c + 0);
    c1 = *(c + 1);

    c0_nxt = *(c + 0 + ldc);
    c1_nxt = *(c + 1 + ldc);

    if (bk > 0)
    {
        FLOAT a0, a1, b0, b1;

        a0 = a[4];
        a1 = a[5];

        b0 = b[4];
        c0 -= a0 * b0;
        c1 -= a1 * b0;

        b1 = b[5];
        c0_nxt -= a0 * b1;
        c1_nxt -= a1 * b1;
    }

    b3 = *(b + 3);
    b2 = *(b + 2);
    b0 = *b;

    c0_nxt *= b3;
    c1_nxt *= b3;

    c0 -= c0_nxt * b2;
    c0 *= b0;

    c1 -= c1_nxt * b2;
    c1 *= b0;

    *(a + 0) = c0;
    *(a + 1) = c1;
    *(a + 2) = c0_nxt;
    *(a + 3) = c1_nxt;

    *(c + 0) = c0;
    *(c + 1) = c1;
    *(c + 0 + ldc) = c0_nxt;
    *(c + 1 + ldc) = c1_nxt;
}

static void dsolve_2x1_rt_msa(FLOAT *a, FLOAT *b, FLOAT *c)
{
    FLOAT b0, c0, c1;

    c0 = *(c + 0);
    c1 = *(c + 1);

    b0 = *b;

    c0 *= b0;
    c1 *= b0;

    *(a + 0) = c0;
    *(a + 1) = c1;

    *(c + 0) = c0;
    *(c + 1) = c1;
}

static void dsolve_1x4_rt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG bk)
{
    FLOAT b0, b4, b5, b8, b9, b10, b12, b13, b14, b15;
    FLOAT c0, c0_nxt1, c0_nxt2, c0_nxt3;

    c0 = *(c + 0);
    c0_nxt1 = *(c + 1 * ldc);
    c0_nxt2 = *(c + 2 * ldc);
    c0_nxt3 = *(c + 3 * ldc);

    if (bk > 0)
    {
        BLASLONG i;
        FLOAT *aa = a + 4, *bb = b + 16;

        for (i = bk; i--;)
        {
            c0 -= aa[0] * bb[0];
            c0_nxt1 -= aa[0] * bb[1];
            c0_nxt2 -= aa[0] * bb[2];
            c0_nxt3 -= aa[0] * bb[3];

            aa += 1;
            bb += 4;
        }
    }

    b0 = *b;
    b4 = *(b + 4);
    b5 = *(b + 5);
    b8 = *(b + 8);
    b9 = *(b + 9);
    b10 = *(b + 10);
    b12 = *(b + 12);
    b13 = *(b + 13);
    b14 = *(b + 14);
    b15 = *(b + 15);

    c0_nxt3 *= b15;

    c0_nxt2 -= c0_nxt3 * b14;
    c0_nxt2 *= b10;

    c0_nxt1 -= c0_nxt3 * b13;
    c0_nxt1 -= c0_nxt2 * b9;
    c0_nxt1 *= b5;

    c0 -= c0_nxt3 * b12;
    c0 -= c0_nxt2 * b8;
    c0 -= c0_nxt1 * b4;
    c0 *= b0;

    *(a + 0) = c0;
    *(a + 1) = c0_nxt1;
    *(a + 2) = c0_nxt2;
    *(a + 3) = c0_nxt3;

    *(c) = c0;
    *(c + 1 * ldc) = c0_nxt1;
    *(c + 2 * ldc) = c0_nxt2;
    *(c + 3 * ldc) = c0_nxt3;
}

static void dsolve_1x2_rt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG bk)
{
    FLOAT b0, b2, b3, c0, c0_nxt;

    c0 = *(c + 0);
    c0_nxt = *(c + ldc);

    if (bk > 0)
    {
        c0 -= a[2] * b[4];
        c0_nxt -= a[2] * b[5];
    }

    b3 = *(b + 3);
    b2 = *(b + 2);
    b0 = *b;

    c0_nxt *= b3;

    c0 -= c0_nxt * b2;
    c0 *= b0;

    *(a + 0) = c0;
    *(a + 1) = c0_nxt;

    *(c + 0) = c0;
    *(c + ldc) = c0_nxt;
}

int CNAME(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT dummy1, FLOAT *a, FLOAT *b,
          FLOAT *c, BLASLONG ldc, BLASLONG offset)
{
    BLASLONG i, j, kk;
    FLOAT *aa, *cc, *bb;

    kk = n;
    c += n * ldc;
    b += n * k;

    if (n & 3)
    {
        if (n & 1)
        {
            aa = a;
            c -= ldc;
            b -= k;
            bb = b + (kk - 1);
            cc = c;

            for (i = (m >> 3); i--;)
            {
                dsolve_8x1_rt_msa(aa + 8 * kk - 8, bb, cc);

                aa += 8 * k;
                cc += 8;
            }

            if (m & 7)
            {
                if (m & 4)
                {
                    dsolve_4x1_rt_msa(aa + 4 * kk - 4, bb, cc);

                    aa += 4 * k;
                    cc += 4;
                }

                if (m & 2)
                {
                    dsolve_2x1_rt_msa(aa + 2 * kk - 2, bb, cc);

                    aa += 2 * k;
                    cc += 2;
                }

                if (m & 1)
                {
                    *cc *= *bb;
                    *(aa + kk - 1) = *cc;

                    aa += k;
                    cc += 1;
                }

            }

            kk -= 1;
        }

        if (n & 2)
        {
            aa = a;
            c -= 2 * ldc;
            b -= 2 * k;
            bb = b + 2 * kk;
            cc = c;

            for (i = (m >> 3); i--;)
            {
                dsolve_8x2_rt_msa(aa + 8 * kk - 16, bb - 4, cc, ldc, k - kk);

                aa += 8 * k;
                cc += 8;
            }

            if (m & 7)
            {
                if (m & 4)
                {
                    dsolve_4x2_rt_msa(aa + 4 * kk - 8, bb - 4, cc, ldc, k - kk);

                    aa += 4 * k;
                    cc += 4;
                }

                if (m & 2)
                {
                    dsolve_2x2_rt_msa(aa + 2 * kk - 4, bb - 4, cc, ldc, k - kk);

                    aa += 2 * k;
                    cc += 2;
                }

                if (m & 1)
                {
                    dsolve_1x2_rt_msa(aa + kk - 2, bb - 4, cc, ldc, k - kk);
                }
            }

            kk -= 2;
        }
    }

    for (j = (n >> 2); j--;)
    {
        aa  = a;
        b -= 4 * k;
        bb = b + 4 * kk;
        c -= 4 * ldc;
        cc = c;

        for (i = (m >> 3); i--;)
        {
            dsolve_8x4_rt_msa(aa + kk * 8, bb, cc, ldc, k - kk);

            aa += 8 * k;
            cc += 8;
        }

        if (m & 7)
        {
            if (m & 4)
            {
                dsolve_4x4_rt_msa(aa + kk * 4 - 16, bb - 16, cc, ldc, k - kk);

                aa += 4 * k;
                cc += 4;
            }

            if (m & 2)
            {
                dsolve_2x4_rt_msa(aa + kk * 2 - 8, bb - 16, cc, ldc, k - kk);

                aa += 2 * k;
                cc += 2;
            }

            if (m & 1)
            {
                dsolve_1x4_rt_msa(aa + kk - 4, bb - 16, cc, ldc, k - kk);

                aa += k;
                cc += 1;
            }
        }

        kk -= 4;
    }

    return 0;
}
