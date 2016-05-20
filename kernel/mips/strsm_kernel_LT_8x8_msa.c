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

static void ssolve_8x8_lt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG bk)
{
    v4f32 src_c0, src_c1, src_c2, src_c3, src_c4, src_c5, src_c6, src_c7;
    v4f32 src_c8, src_c9, src_c10, src_c11, src_c12, src_c13, src_c14, src_c15;
    v4f32 res_c0, res_c1, res_c2, res_c3, res_c4, res_c5, res_c6, res_c7;
    v4f32 res_c8, res_c9, res_c10, res_c11, res_c12, res_c13, res_c14, res_c15;
    v4f32 src_a0, src_a1, src_a2, src_a3, src_a4, src_a5, src_a6, src_a7;
    v4f32 src_a9, src_a10, src_a11, src_a12, src_a13, src_a14, src_a15, src_a18;
    v4f32 src_a19, src_a20, src_a21, src_a22, src_a23, src_a27, src_a28;
    v4f32 src_a29, src_a30, src_a31, src_a36, src_a37, src_a38, src_a39;
    v4f32 src_a45, src_a46, src_a47, src_a54, src_a55, src_a63, src_a;
    FLOAT *c_nxt1line = c + ldc;
    FLOAT *c_nxt2line = c + 2 * ldc;
    FLOAT *c_nxt3line = c + 3 * ldc;
    FLOAT *c_nxt4line = c + 4 * ldc;
    FLOAT *c_nxt5line = c + 5 * ldc;
    FLOAT *c_nxt6line = c + 6 * ldc;
    FLOAT *c_nxt7line = c + 7 * ldc;

    if (bk)
    {
        BLASLONG k;
        v4f32 src_b, src_b0, src_b1, src_b2, src_b3;
        v4f32 res0, res1, res2, res3, res4, res5, res6, res7;
        v4f32 res8, res9, res10, res11, res12, res13, res14, res15;

        LD_SP2(a, 4, src_a0, src_a1);

        src_b = LD_SP(b + 0);
        SPLATI_W4_SP(src_b, src_b0, src_b1, src_b2, src_b3);
        res0 = src_a0 * src_b0;
        res1 = src_a1 * src_b0;
        res2 = src_a0 * src_b1;
        res3 = src_a1 * src_b1;
        res4 = src_a0 * src_b2;
        res5 = src_a1 * src_b2;
        res6 = src_a0 * src_b3;
        res7 = src_a1 * src_b3;

        src_b = LD_SP(b + 4);
        SPLATI_W4_SP(src_b, src_b0, src_b1, src_b2, src_b3);
        res8 = src_a0 * src_b0;
        res9 = src_a1 * src_b0;
        res10 = src_a0 * src_b1;
        res11 = src_a1 * src_b1;
        res12 = src_a0 * src_b2;
        res13 = src_a1 * src_b2;
        res14 = src_a0 * src_b3;
        res15 = src_a1 * src_b3;

        a += 8;
        b += 8;

        for (k = (bk - 1); k--;)
        {
            LD_SP2(a, 4, src_a0, src_a1);

            src_b = LD_SP(b + 0);
            SPLATI_W4_SP(src_b, src_b0, src_b1, src_b2, src_b3);
            res0 += src_a0 * src_b0;
            res1 += src_a1 * src_b0;
            res2 += src_a0 * src_b1;
            res3 += src_a1 * src_b1;
            res4 += src_a0 * src_b2;
            res5 += src_a1 * src_b2;
            res6 += src_a0 * src_b3;
            res7 += src_a1 * src_b3;

            src_b = LD_SP(b + 4);
            SPLATI_W4_SP(src_b, src_b0, src_b1, src_b2, src_b3);
            res8 += src_a0 * src_b0;
            res9 += src_a1 * src_b0;
            res10 += src_a0 * src_b1;
            res11 += src_a1 * src_b1;
            res12 += src_a0 * src_b2;
            res13 += src_a1 * src_b2;
            res14 += src_a0 * src_b3;
            res15 += src_a1 * src_b3;

            a += 8;
            b += 8;
        }

        LD_SP2(c, 4, src_c0, src_c1);
        LD_SP2(c_nxt1line, 4, src_c2, src_c3);
        LD_SP2(c_nxt2line, 4, src_c4, src_c5);
        LD_SP2(c_nxt3line, 4, src_c6, src_c7);
        LD_SP2(c_nxt4line, 4, src_c8, src_c9);
        LD_SP2(c_nxt5line, 4, src_c10, src_c11);
        LD_SP2(c_nxt6line, 4, src_c12, src_c13);
        LD_SP2(c_nxt7line, 4, src_c14, src_c15);

        src_c0 -= res0;
        src_c1 -= res1;
        src_c2 -= res2;
        src_c3 -= res3;
        src_c4 -= res4;
        src_c5 -= res5;
        src_c6 -= res6;
        src_c7 -= res7;
        src_c8 -= res8;
        src_c9 -= res9;
        src_c10 -= res10;
        src_c11 -= res11;
        src_c12 -= res12;
        src_c13 -= res13;
        src_c14 -= res14;
        src_c15 -= res15;
    }
    else
    {
        LD_SP2(c, 4, src_c0, src_c1);
        LD_SP2(c_nxt1line, 4, src_c2, src_c3);
        LD_SP2(c_nxt2line, 4, src_c4, src_c5);
        LD_SP2(c_nxt3line, 4, src_c6, src_c7);
        LD_SP2(c_nxt4line, 4, src_c8, src_c9);
        LD_SP2(c_nxt5line, 4, src_c10, src_c11);
        LD_SP2(c_nxt6line, 4, src_c12, src_c13);
        LD_SP2(c_nxt7line, 4, src_c14, src_c15);
    }

    TRANSPOSE4x4_SP_SP(src_c0, src_c2, src_c4, src_c6,
                       res_c0, res_c1, res_c2, res_c3);
    TRANSPOSE4x4_SP_SP(src_c8, src_c10, src_c12, src_c14,
                       res_c8, res_c9, res_c10, res_c11);
    TRANSPOSE4x4_SP_SP(src_c1, src_c3, src_c5, src_c7,
                       res_c4, res_c5, res_c6, res_c7);
    TRANSPOSE4x4_SP_SP(src_c9, src_c11, src_c13, src_c15,
                       res_c12, res_c13, res_c14, res_c15);

    src_a = LD_SP(a + 0);
    SPLATI_W4_SP(src_a, src_a0, src_a1, src_a2, src_a3);
    src_a = LD_SP(a + 4);
    SPLATI_W4_SP(src_a, src_a4, src_a5, src_a6, src_a7);

    res_c0 *= src_a0;
    res_c8 *= src_a0;
    res_c1 -= res_c0 * src_a1;
    res_c9 -= res_c8 * src_a1;
    res_c2 -= res_c0 * src_a2;
    res_c10 -= res_c8 * src_a2;
    res_c3 -= res_c0 * src_a3;
    res_c11 -= res_c8 * src_a3;
    res_c4 -= res_c0 * src_a4;
    res_c12 -= res_c8 * src_a4;
    res_c5 -= res_c0 * src_a5;
    res_c13 -= res_c8 * src_a5;
    res_c6 -= res_c0 * src_a6;
    res_c14 -= res_c8 * src_a6;
    res_c7 -= res_c0 * src_a7;
    res_c15 -= res_c8 * src_a7;

    src_a = LD_SP(a + 9);
    SPLATI_W4_SP(src_a, src_a9, src_a10, src_a11, src_a12);
    src_a13 = LD_SP(a + 13);
    src_a15 = (v4f32) __msa_splati_w((v4i32) src_a13, 2);
    src_a14 = (v4f32) __msa_splati_w((v4i32) src_a13, 1);
    src_a13 = (v4f32) __msa_splati_w((v4i32) src_a13, 0);

    res_c1 *= src_a9;
    res_c9 *= src_a9;
    res_c2 -= res_c1 * src_a10;
    res_c10 -= res_c9 * src_a10;
    res_c3 -= res_c1 * src_a11;
    res_c11 -= res_c9 * src_a11;
    res_c4 -= res_c1 * src_a12;
    res_c12 -= res_c9 * src_a12;
    res_c5 -= res_c1 * src_a13;
    res_c13 -= res_c9 * src_a13;
    res_c6 -= res_c1 * src_a14;
    res_c14 -= res_c9 * src_a14;
    res_c7 -= res_c1 * src_a15;
    res_c15 -= res_c9 * src_a15;

    src_a = LD_SP(a + 18);
    SPLATI_W4_SP(src_a, src_a18, src_a19, src_a20, src_a21);
    src_a22 = LD_SP(a + 22);
    src_a23 = (v4f32) __msa_splati_w((v4i32) src_a22, 1);
    src_a22 = (v4f32) __msa_splati_w((v4i32) src_a22, 0);

    res_c2 *= src_a18;
    res_c10 *= src_a18;
    res_c3 -= res_c2 * src_a19;
    res_c11 -= res_c10 * src_a19;
    res_c4 -= res_c2 * src_a20;
    res_c12 -= res_c10 * src_a20;
    res_c5 -= res_c2 * src_a21;
    res_c13 -= res_c10 * src_a21;
    res_c6 -= res_c2 * src_a22;
    res_c14 -= res_c10 * src_a22;
    res_c7 -= res_c2 * src_a23;
    res_c15 -= res_c10 * src_a23;

    src_a = LD_SP(a + 27);
    SPLATI_W4_SP(src_a, src_a27, src_a28, src_a29, src_a30);
    src_a31 = __msa_cast_to_vector_float(*(a + 31));
    src_a31 = (v4f32) __msa_splati_w((v4i32) src_a31, 0);

    res_c3 *= src_a27;
    res_c11 *= src_a27;
    res_c4 -= res_c3 * src_a28;
    res_c12 -= res_c11 * src_a28;
    res_c5 -= res_c3 * src_a29;
    res_c13 -= res_c11 * src_a29;
    res_c6 -= res_c3 * src_a30;
    res_c14 -= res_c11 * src_a30;
    res_c7 -= res_c3 * src_a31;
    res_c15 -= res_c11 * src_a31;

    ST_SP4(res_c0, res_c8, res_c1, res_c9, b, 4);
    ST_SP4(res_c2, res_c10, res_c3, res_c11, b + 16, 4);

    TRANSPOSE4x4_SP_SP(res_c0, res_c1, res_c2, res_c3,
                       src_c0, src_c2, src_c4, src_c6);
    TRANSPOSE4x4_SP_SP(res_c8, res_c9, res_c10, res_c11,
                       src_c8, src_c10, src_c12, src_c14);

    ST_SP(src_c0, c);
    ST_SP(src_c2, c_nxt1line);
    ST_SP(src_c4, c_nxt2line);
    ST_SP(src_c6, c_nxt3line);
    ST_SP(src_c8, c_nxt4line);
    ST_SP(src_c10, c_nxt5line);
    ST_SP(src_c12, c_nxt6line);
    ST_SP(src_c14, c_nxt7line);

    src_a = LD_SP(a + 36);
    SPLATI_W4_SP(src_a, src_a36, src_a37, src_a38, src_a39);

    res_c4 *= src_a36;
    res_c12 *= src_a36;
    res_c5 -= res_c4 * src_a37;
    res_c13 -= res_c12 * src_a37;
    res_c6 -= res_c4 * src_a38;
    res_c14 -= res_c12 * src_a38;
    res_c7 -= res_c4 * src_a39;
    res_c15 -= res_c12 * src_a39;

    src_a45 = LD_SP(a + 45);
    src_a47 = (v4f32) __msa_splati_w((v4i32) src_a45, 2);
    src_a46 = (v4f32) __msa_splati_w((v4i32) src_a45, 1);
    src_a45 = (v4f32) __msa_splati_w((v4i32) src_a45, 0);

    res_c5 *= src_a45;
    res_c13 *= src_a45;
    res_c6 -= res_c5 * src_a46;
    res_c14 -= res_c13 * src_a46;
    res_c7 -= res_c5 * src_a47;
    res_c15 -= res_c13 * src_a47;

    src_a54 = __msa_cast_to_vector_float(*(a + 54));
    src_a54 = (v4f32) __msa_splati_w((v4i32) src_a54, 0);
    src_a55 = __msa_cast_to_vector_float(*(a + 55));
    src_a55 = (v4f32) __msa_splati_w((v4i32) src_a55, 0);
    src_a63 = __msa_cast_to_vector_float(*(a + 63));
    src_a63 = (v4f32) __msa_splati_w((v4i32) src_a63, 0);

    res_c6 *= src_a54;
    res_c14 *= src_a54;
    res_c7 -= res_c6 * src_a55;
    res_c15 -= res_c14 * src_a55;

    res_c7 *= src_a63;
    res_c15 *= src_a63;

    ST_SP4(res_c4, res_c12, res_c5, res_c13, b + 32, 4);
    ST_SP4(res_c6, res_c14, res_c7, res_c15, b + 48, 4);

    TRANSPOSE4x4_SP_SP(res_c4, res_c5, res_c6, res_c7,
                       src_c1, src_c3, src_c5, src_c7);
    TRANSPOSE4x4_SP_SP(res_c12, res_c13, res_c14, res_c15,
                       src_c9, src_c11, src_c13, src_c15);

    ST_SP(src_c1, c + 4);
    ST_SP(src_c3, c_nxt1line + 4);
    ST_SP(src_c5, c_nxt2line + 4);
    ST_SP(src_c7, c_nxt3line + 4);
    ST_SP(src_c9, c_nxt4line + 4);
    ST_SP(src_c11, c_nxt5line + 4);
    ST_SP(src_c13, c_nxt6line + 4);
    ST_SP(src_c15, c_nxt7line + 4);
}

static void ssolve_8x4_lt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG bk)
{
    v4f32 src_c0, src_c1, src_c2, src_c3, src_c4, src_c5, src_c6, src_c7;
    v4f32 res_c0, res_c1, res_c2, res_c3, res_c4, res_c5, res_c6, res_c7;
    v4f32 src_a0, src_a1, src_a2, src_a3, src_a4, src_a5, src_a6, src_a7;
    v4f32 src_a9, src_a10, src_a11, src_a12, src_a13, src_a14, src_a15, src_a18;
    v4f32 src_a19, src_a20, src_a21, src_a22, src_a23, src_a27, src_a28;
    v4f32 src_a29, src_a30, src_a31, src_a36, src_a37, src_a38, src_a39;
    v4f32 src_a45, src_a46, src_a47, src_a54, src_a55, src_a63, src_a;
    FLOAT *c_nxt1line = c + ldc;
    FLOAT *c_nxt2line = c + 2 * ldc;
    FLOAT *c_nxt3line = c + 3 * ldc;

    if (bk)
    {
        BLASLONG k;
        v4f32 src_b, src_b0, src_b1, src_b2, src_b3;
        v4f32 res0, res1, res2, res3, res4, res5, res6, res7;

        LD_SP2(a, 4, src_a0, src_a1);

        src_b = LD_SP(b + 0);
        SPLATI_W4_SP(src_b, src_b0, src_b1, src_b2, src_b3);
        res0 = src_a0 * src_b0;
        res1 = src_a1 * src_b0;
        res2 = src_a0 * src_b1;
        res3 = src_a1 * src_b1;
        res4 = src_a0 * src_b2;
        res5 = src_a1 * src_b2;
        res6 = src_a0 * src_b3;
        res7 = src_a1 * src_b3;

        a += 8;
        b += 4;

        for (k = (bk - 1); k--;)
        {
            LD_SP2(a, 4, src_a0, src_a1);

            src_b = LD_SP(b + 0);
            SPLATI_W4_SP(src_b, src_b0, src_b1, src_b2, src_b3);
            res0 += src_a0 * src_b0;
            res1 += src_a1 * src_b0;
            res2 += src_a0 * src_b1;
            res3 += src_a1 * src_b1;
            res4 += src_a0 * src_b2;
            res5 += src_a1 * src_b2;
            res6 += src_a0 * src_b3;
            res7 += src_a1 * src_b3;

            a += 8;
            b += 4;
        }

        LD_SP2(c, 4, src_c0, src_c1);
        LD_SP2(c_nxt1line, 4, src_c2, src_c3);
        LD_SP2(c_nxt2line, 4, src_c4, src_c5);
        LD_SP2(c_nxt3line, 4, src_c6, src_c7);

        src_c0 -= res0;
        src_c1 -= res1;
        src_c2 -= res2;
        src_c3 -= res3;
        src_c4 -= res4;
        src_c5 -= res5;
        src_c6 -= res6;
        src_c7 -= res7;
    }
    else
    {
        LD_SP2(c, 4, src_c0, src_c1);
        LD_SP2(c_nxt1line, 4, src_c2, src_c3);
        LD_SP2(c_nxt2line, 4, src_c4, src_c5);
        LD_SP2(c_nxt3line, 4, src_c6, src_c7);
    }

    TRANSPOSE4x4_SP_SP(src_c0, src_c2, src_c4, src_c6,
                       res_c0, res_c1, res_c2, res_c3);
    TRANSPOSE4x4_SP_SP(src_c1, src_c3, src_c5, src_c7,
                       res_c4, res_c5, res_c6, res_c7);

    src_a = LD_SP(a + 0);
    SPLATI_W4_SP(src_a, src_a0, src_a1, src_a2, src_a3);
    src_a = LD_SP(a + 4);
    SPLATI_W4_SP(src_a, src_a4, src_a5, src_a6, src_a7);

    res_c0 *= src_a0;
    res_c1 -= res_c0 * src_a1;
    res_c2 -= res_c0 * src_a2;
    res_c3 -= res_c0 * src_a3;
    res_c4 -= res_c0 * src_a4;
    res_c5 -= res_c0 * src_a5;
    res_c6 -= res_c0 * src_a6;
    res_c7 -= res_c0 * src_a7;

    src_a = LD_SP(a + 9);
    SPLATI_W4_SP(src_a, src_a9, src_a10, src_a11, src_a12);
    src_a13 = LD_SP(a + 13);
    src_a15 = (v4f32) __msa_splati_w((v4i32) src_a13, 2);
    src_a14 = (v4f32) __msa_splati_w((v4i32) src_a13, 1);
    src_a13 = (v4f32) __msa_splati_w((v4i32) src_a13, 0);

    res_c1 *= src_a9;
    res_c2 -= res_c1 * src_a10;
    res_c3 -= res_c1 * src_a11;
    res_c4 -= res_c1 * src_a12;
    res_c5 -= res_c1 * src_a13;
    res_c6 -= res_c1 * src_a14;
    res_c7 -= res_c1 * src_a15;

    src_a = LD_SP(a + 18);
    SPLATI_W4_SP(src_a, src_a18, src_a19, src_a20, src_a21);
    src_a22 = LD_SP(a + 22);
    src_a23 = (v4f32) __msa_splati_w((v4i32) src_a22, 1);
    src_a22 = (v4f32) __msa_splati_w((v4i32) src_a22, 0);

    res_c2 *= src_a18;
    res_c3 -= res_c2 * src_a19;
    res_c4 -= res_c2 * src_a20;
    res_c5 -= res_c2 * src_a21;
    res_c6 -= res_c2 * src_a22;
    res_c7 -= res_c2 * src_a23;

    src_a = LD_SP(a + 27);
    SPLATI_W4_SP(src_a, src_a27, src_a28, src_a29, src_a30);
    src_a31 = __msa_cast_to_vector_float(*(a + 31));
    src_a31 = (v4f32) __msa_splati_w((v4i32) src_a31, 0);

    res_c3 *= src_a27;
    res_c4 -= res_c3 * src_a28;
    res_c5 -= res_c3 * src_a29;
    res_c6 -= res_c3 * src_a30;
    res_c7 -= res_c3 * src_a31;

    src_a = LD_SP(a + 36);
    SPLATI_W4_SP(src_a, src_a36, src_a37, src_a38, src_a39);

    res_c4 *= src_a36;
    res_c5 -= res_c4 * src_a37;
    res_c6 -= res_c4 * src_a38;
    res_c7 -= res_c4 * src_a39;

    src_a45 = LD_SP(a + 45);
    src_a47 = (v4f32) __msa_splati_w((v4i32) src_a45, 2);
    src_a46 = (v4f32) __msa_splati_w((v4i32) src_a45, 1);
    src_a45 = (v4f32) __msa_splati_w((v4i32) src_a45, 0);

    res_c5 *= src_a45;
    res_c6 -= res_c5 * src_a46;
    res_c7 -= res_c5 * src_a47;

    src_a54 = __msa_cast_to_vector_float(*(a + 54));
    src_a54 = (v4f32) __msa_splati_w((v4i32) src_a54, 0);
    src_a55 = __msa_cast_to_vector_float(*(a + 55));
    src_a55 = (v4f32) __msa_splati_w((v4i32) src_a55, 0);
    src_a63 = __msa_cast_to_vector_float(*(a + 63));
    src_a63 = (v4f32) __msa_splati_w((v4i32) src_a63, 0);

    res_c6 *= src_a54;
    res_c7 -= res_c6 * src_a55;
    res_c7 *= src_a63;

    ST_SP4(res_c0, res_c1, res_c2, res_c3, b, 4);
    b += 16;
    ST_SP4(res_c4, res_c5, res_c6, res_c7, b, 4);

    TRANSPOSE4x4_SP_SP(res_c0, res_c1, res_c2, res_c3,
                       src_c0, src_c2, src_c4, src_c6);
    TRANSPOSE4x4_SP_SP(res_c4, res_c5, res_c6, res_c7,
                       src_c1, src_c3, src_c5, src_c7);

    ST_SP2(src_c0, src_c1, c, 4);
    ST_SP2(src_c2, src_c3, c_nxt1line, 4);
    ST_SP2(src_c4, src_c5, c_nxt2line, 4);
    ST_SP2(src_c6, src_c7, c_nxt3line, 4);
}

static void ssolve_8x2_lt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG bk)
{
    FLOAT a0, a1, a2, a3, a4, a5, a6, a7, a9, a10, a11, a12, a13, a14, a15, a18;
    FLOAT a19, a20, a21, a22, a23, a27, a28, a29, a30, a31, a36, a37, a38, a39;
    FLOAT a45, a46, a47, a54, a55, a63;
    FLOAT c0, c1, c2, c3, c4, c5, c6, c7;
    FLOAT c0_nxt, c1_nxt, c2_nxt, c3_nxt, c4_nxt, c5_nxt, c6_nxt, c7_nxt;

    c0 = *(c + 0);
    c1 = *(c + 1);
    c2 = *(c + 2);
    c3 = *(c + 3);
    c4 = *(c + 4);
    c5 = *(c + 5);
    c6 = *(c + 6);
    c7 = *(c + 7);
    c0_nxt = *(c + 0 + ldc);
    c1_nxt = *(c + 1 + ldc);
    c2_nxt = *(c + 2 + ldc);
    c3_nxt = *(c + 3 + ldc);
    c4_nxt = *(c + 4 + ldc);
    c5_nxt = *(c + 5 + ldc);
    c6_nxt = *(c + 6 + ldc);
    c7_nxt = *(c + 7 + ldc);

    if (bk > 0)
    {
        BLASLONG k;
        FLOAT res[16];

        res[0] = a[0] * b[0];
        res[1] = a[1] * b[0];
        res[2] = a[2] * b[0];
        res[3] = a[3] * b[0];
        res[4] = a[4] * b[0];
        res[5] = a[5] * b[0];
        res[6] = a[6] * b[0];
        res[7] = a[7] * b[0];
        res[8]  = a[0] * b[1];
        res[9]  = a[1] * b[1];
        res[10] = a[2] * b[1];
        res[11] = a[3] * b[1];
        res[12] = a[4] * b[1];
        res[13] = a[5] * b[1];
        res[14] = a[6] * b[1];
        res[15] = a[7] * b[1];

        for (k = (bk - 1); k--;)
        {
            a += 8;
            b += 2;

            res[0] += a[0] * b[0];
            res[1] += a[1] * b[0];
            res[2] += a[2] * b[0];
            res[3] += a[3] * b[0];
            res[4] += a[4] * b[0];
            res[5] += a[5] * b[0];
            res[6] += a[6] * b[0];
            res[7] += a[7] * b[0];
            res[8]  += a[0] * b[1];
            res[9]  += a[1] * b[1];
            res[10] += a[2] * b[1];
            res[11] += a[3] * b[1];
            res[12] += a[4] * b[1];
            res[13] += a[5] * b[1];
            res[14] += a[6] * b[1];
            res[15] += a[7] * b[1];
        }

        c0 -= res[0];
        c1 -= res[1];
        c2 -= res[2];
        c3 -= res[3];
        c4 -= res[4];
        c5 -= res[5];
        c6 -= res[6];
        c7 -= res[7];
        c0_nxt -= res[8];
        c1_nxt -= res[9];
        c2_nxt -= res[10];
        c3_nxt -= res[11];
        c4_nxt -= res[12];
        c5_nxt -= res[13];
        c6_nxt -= res[14];
        c7_nxt -= res[15];

        a += 8;
        b += 2;
    }

    a0 = *(a + 0);
    a1 = *(a + 1);
    a2 = *(a + 2);
    a3 = *(a + 3);
    a4 = *(a + 4);
    a5 = *(a + 5);
    a6 = *(a + 6);
    a7 = *(a + 7);
    a9 = *(a + 9);
    a10 = *(a + 10);
    a11 = *(a + 11);
    a12 = *(a + 12);
    a13 = *(a + 13);
    a14 = *(a + 14);
    a15 = *(a + 15);
    a18 = *(a + 18);
    a19 = *(a + 19);
    a20 = *(a + 20);
    a21 = *(a + 21);
    a22 = *(a + 22);
    a23 = *(a + 23);
    a27 = *(a + 27);
    a28 = *(a + 28);
    a29 = *(a + 29);
    a30 = *(a + 30);
    a31 = *(a + 31);
    a36 = *(a + 36);
    a37 = *(a + 37);
    a38 = *(a + 38);
    a39 = *(a + 39);
    a45 = *(a + 45);
    a46 = *(a + 46);
    a47 = *(a + 47);
    a54 = *(a + 54);
    a55 = *(a + 55);
    a63 = *(a + 63);

    c0 *= a0;
    c0_nxt *= a0;

    c1 -= c0 * a1;
    c1_nxt -= c0_nxt * a1;
    c1 *= a9;
    c1_nxt *= a9;

    c2 -= c0 * a2;
    c2_nxt -= c0_nxt * a2;
    c2 -= c1 * a10;
    c2_nxt -= c1_nxt * a10;
    c2 *= a18;
    c2_nxt *= a18;

    c3 -= c0 * a3;
    c3_nxt -= c0_nxt * a3;
    c3 -= c1 * a11;
    c3_nxt -= c1_nxt * a11;
    c3 -= c2 * a19;
    c3_nxt -= c2_nxt * a19;
    c3 *= a27;
    c3_nxt *= a27;

    c4 -= c0 * a4;
    c4_nxt -= c0_nxt * a4;
    c4 -= c1 * a12;
    c4_nxt -= c1_nxt * a12;
    c4 -= c2 * a20;
    c4_nxt -= c2_nxt * a20;
    c4 -= c3 * a28;
    c4_nxt -= c3_nxt * a28;
    c4 *= a36;
    c4_nxt *= a36;

    c5 -= c0 * a5;
    c5_nxt -= c0_nxt * a5;
    c5 -= c1 * a13;
    c5_nxt -= c1_nxt * a13;
    c5 -= c2 * a21;
    c5_nxt -= c2_nxt * a21;
    c5 -= c3 * a29;
    c5_nxt -= c3_nxt * a29;
    c5 -= c4 * a37;
    c5_nxt -= c4_nxt * a37;
    c5 *= a45;
    c5_nxt *= a45;

    c6 -= c0 * a6;
    c6_nxt -= c0_nxt * a6;
    c6 -= c1 * a14;
    c6_nxt -= c1_nxt * a14;
    c6 -= c2 * a22;
    c6_nxt -= c2_nxt * a22;
    c6 -= c3 * a30;
    c6_nxt -= c3_nxt * a30;
    c6 -= c4 * a38;
    c6_nxt -= c4_nxt * a38;
    c6 -= c5 * a46;
    c6_nxt -= c5_nxt * a46;
    c6 *= a54;
    c6_nxt *= a54;

    c7 -= c0 * a7;
    c7_nxt -= c0_nxt * a7;
    c7 -= c1 * a15;
    c7_nxt -= c1_nxt * a15;
    c7 -= c2 * a23;
    c7_nxt -= c2_nxt * a23;
    c7 -= c3 * a31;
    c7_nxt -= c3_nxt * a31;
    c7 -= c4 * a39;
    c7_nxt -= c4_nxt * a39;
    c7 -= c5 * a47;
    c7_nxt -= c5_nxt * a47;
    c7 -= c6 * a55;
    c7_nxt -= c6_nxt * a55;
    c7 *= a63;
    c7_nxt *= a63;

    *(c + 0) = c0;
    *(c + 1) = c1;
    *(c + 2) = c2;
    *(c + 3) = c3;
    *(c + 4) = c4;
    *(c + 5) = c5;
    *(c + 6) = c6;
    *(c + 7) = c7;
    *(c + 0 + ldc) = c0_nxt;
    *(c + 1 + ldc) = c1_nxt;
    *(c + 2 + ldc) = c2_nxt;
    *(c + 3 + ldc) = c3_nxt;
    *(c + 4 + ldc) = c4_nxt;
    *(c + 5 + ldc) = c5_nxt;
    *(c + 6 + ldc) = c6_nxt;
    *(c + 7 + ldc) = c7_nxt;

    *(b + 0) = c0;
    *(b + 1) = c0_nxt;
    *(b + 2) = c1;
    *(b + 3) = c1_nxt;
    *(b + 4) = c2;
    *(b + 5) = c2_nxt;
    *(b + 6) = c3;
    *(b + 7) = c3_nxt;
    *(b + 8) = c4;
    *(b + 9) = c4_nxt;
    *(b + 10) = c5;
    *(b + 11) = c5_nxt;
    *(b + 12) = c6;
    *(b + 13) = c6_nxt;
    *(b + 14) = c7;
    *(b + 15) = c7_nxt;
}

static void ssolve_8x1_lt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG bk)
{
    FLOAT a0, a1, a2, a3, a4, a5, a6, a7, a9, a10, a11, a12, a13, a14, a15, a18;
    FLOAT a19, a20, a21, a22, a23, a27, a28, a29, a30, a31, a36, a37, a38, a39;
    FLOAT a45, a46, a47, a54, a55, a63, c0, c1, c2, c3, c4, c5, c6, c7;

    c0 = *(c + 0);
    c1 = *(c + 1);
    c2 = *(c + 2);
    c3 = *(c + 3);
    c4 = *(c + 4);
    c5 = *(c + 5);
    c6 = *(c + 6);
    c7 = *(c + 7);

    if (bk > 0)
    {
        BLASLONG i;
        FLOAT a0, a1, a2, a3, a4, a5, a6, a7;

        a0 = a[0] * b[0];
        a1 = a[1] * b[0];
        a2 = a[2] * b[0];
        a3 = a[3] * b[0];
        a4 = a[4] * b[0];
        a5 = a[5] * b[0];
        a6 = a[6] * b[0];
        a7 = a[7] * b[0];

        for (i = (bk - 1); i--; )
        {
            a += 8;
            b += 1;

            a0 += a[0] * b[0];
            a1 += a[1] * b[0];
            a2 += a[2] * b[0];
            a3 += a[3] * b[0];
            a4 += a[4] * b[0];
            a5 += a[5] * b[0];
            a6 += a[6] * b[0];
            a7 += a[7] * b[0];
        }

        c0 -= a0;
        c1 -= a1;
        c2 -= a2;
        c3 -= a3;
        c4 -= a4;
        c5 -= a5;
        c6 -= a6;
        c7 -= a7;

        a += 8;
        b += 1;
    }

    a0 = *(a + 0);
    a1 = *(a + 1);
    a2 = *(a + 2);
    a3 = *(a + 3);
    a4 = *(a + 4);
    a5 = *(a + 5);
    a6 = *(a + 6);
    a7 = *(a + 7);
    a9 = *(a + 9);
    a10 = *(a + 10);
    a11 = *(a + 11);
    a12 = *(a + 12);
    a13 = *(a + 13);
    a14 = *(a + 14);
    a15 = *(a + 15);
    a18 = *(a + 18);
    a19 = *(a + 19);
    a20 = *(a + 20);
    a21 = *(a + 21);
    a22 = *(a + 22);
    a23 = *(a + 23);
    a27 = *(a + 27);
    a28 = *(a + 28);
    a29 = *(a + 29);
    a30 = *(a + 30);
    a31 = *(a + 31);
    a36 = *(a + 36);
    a37 = *(a + 37);
    a38 = *(a + 38);
    a39 = *(a + 39);
    a45 = *(a + 45);
    a46 = *(a + 46);
    a47 = *(a + 47);
    a54 = *(a + 54);
    a55 = *(a + 55);
    a63 = *(a + 63);

    c0 *= a0;

    c1 -= c0 * a1;
    c1 *= a9;

    c2 -= c0 * a2;
    c2 -= c1 * a10;
    c2 *= a18;

    c3 -= c0 * a3;
    c3 -= c1 * a11;
    c3 -= c2 * a19;
    c3 *= a27;

    c4 -= c0 * a4;
    c4 -= c1 * a12;
    c4 -= c2 * a20;
    c4 -= c3 * a28;
    c4 *= a36;

    c5 -= c0 * a5;
    c5 -= c1 * a13;
    c5 -= c2 * a21;
    c5 -= c3 * a29;
    c5 -= c4 * a37;
    c5 *= a45;

    c6 -= c0 * a6;
    c6 -= c1 * a14;
    c6 -= c2 * a22;
    c6 -= c3 * a30;
    c6 -= c4 * a38;
    c6 -= c5 * a46;
    c6 *= a54;

    c7 -= c0 * a7;
    c7 -= c1 * a15;
    c7 -= c2 * a23;
    c7 -= c3 * a31;
    c7 -= c4 * a39;
    c7 -= c5 * a47;
    c7 -= c6 * a55;
    c7 *= a63;

    *(c + 0) = c0;
    *(c + 1) = c1;
    *(c + 2) = c2;
    *(c + 3) = c3;
    *(c + 4) = c4;
    *(c + 5) = c5;
    *(c + 6) = c6;
    *(c + 7) = c7;

    *(b + 0) = c0;
    *(b + 1) = c1;
    *(b + 2) = c2;
    *(b + 3) = c3;
    *(b + 4) = c4;
    *(b + 5) = c5;
    *(b + 6) = c6;
    *(b + 7) = c7;
}

static void ssolve_4x8_lt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG bk)
{
    v4f32 src_c0, src_c1, src_c2, src_c3, src_c4, src_c5, src_c6, src_c7;
    v4f32 res_c0, res_c1, res_c2, res_c3, res_c4, res_c5, res_c6, res_c7;
    v4f32 src_a0, src_a1, src_a2, src_a3, src_a5, src_a6, src_a7;
    v4f32 src_a10, src_a11, src_a15, src_a;
    FLOAT *c_nxt1line = c + ldc;
    FLOAT *c_nxt2line = c + 2 * ldc;
    FLOAT *c_nxt3line = c + 3 * ldc;
    FLOAT *c_nxt4line = c + 4 * ldc;
    FLOAT *c_nxt5line = c + 5 * ldc;
    FLOAT *c_nxt6line = c + 6 * ldc;
    FLOAT *c_nxt7line = c + 7 * ldc;

    if (bk > 0)
    {
        BLASLONG k;
        v4f32 src_b, src_b0, src_b1, src_b2, src_b3;
        v4f32 res0, res1, res2, res3, res4, res5, res6, res7;

        src_a0 = LD_SP(a);

        src_b = LD_SP(b + 0);
        SPLATI_W4_SP(src_b, src_b0, src_b1, src_b2, src_b3);
        res0 = src_a0 * src_b0;
        res1 = src_a0 * src_b1;
        res2 = src_a0 * src_b2;
        res3 = src_a0 * src_b3;

        src_b = LD_SP(b + 4);
        SPLATI_W4_SP(src_b, src_b0, src_b1, src_b2, src_b3);
        res4 = src_a0 * src_b0;
        res5 = src_a0 * src_b1;
        res6 = src_a0 * src_b2;
        res7 = src_a0 * src_b3;

        a += 4;
        b += 8;

        for (k = (bk - 1); k--;)
        {
            src_a0 = LD_SP(a);

            src_b = LD_SP(b + 0);
            SPLATI_W4_SP(src_b, src_b0, src_b1, src_b2, src_b3);
            res0 += src_a0 * src_b0;
            res1 += src_a0 * src_b1;
            res2 += src_a0 * src_b2;
            res3 += src_a0 * src_b3;

            src_b = LD_SP(b + 4);
            SPLATI_W4_SP(src_b, src_b0, src_b1, src_b2, src_b3);
            res4 += src_a0 * src_b0;
            res5 += src_a0 * src_b1;
            res6 += src_a0 * src_b2;
            res7 += src_a0 * src_b3;

            a += 4;
            b += 8;
        }

        src_c0 = LD_SP(c);
        src_c1 = LD_SP(c_nxt1line);
        src_c2 = LD_SP(c_nxt2line);
        src_c3 = LD_SP(c_nxt3line);
        src_c4 = LD_SP(c_nxt4line);
        src_c5 = LD_SP(c_nxt5line);
        src_c6 = LD_SP(c_nxt6line);
        src_c7 = LD_SP(c_nxt7line);

        src_c0 -= res0;
        src_c1 -= res1;
        src_c2 -= res2;
        src_c3 -= res3;
        src_c4 -= res4;
        src_c5 -= res5;
        src_c6 -= res6;
        src_c7 -= res7;
    }
    else
    {
        src_c0 = LD_SP(c);
        src_c1 = LD_SP(c_nxt1line);
        src_c2 = LD_SP(c_nxt2line);
        src_c3 = LD_SP(c_nxt3line);
        src_c4 = LD_SP(c_nxt4line);
        src_c5 = LD_SP(c_nxt5line);
        src_c6 = LD_SP(c_nxt6line);
        src_c7 = LD_SP(c_nxt7line);
    }

    TRANSPOSE4x4_SP_SP(src_c0, src_c1, src_c2, src_c3,
                       res_c0, res_c1, res_c2, res_c3);
    TRANSPOSE4x4_SP_SP(src_c4, src_c5, src_c6, src_c7,
                       res_c4, res_c5, res_c6, res_c7);

    src_a = LD_SP(a + 0);
    SPLATI_W4_SP(src_a, src_a0, src_a1, src_a2, src_a3);
    src_a5 = LD_SP(a + 5);
    src_a7 = (v4f32) __msa_splati_w((v4i32) src_a5, 2);
    src_a6 = (v4f32) __msa_splati_w((v4i32) src_a5, 1);
    src_a5 = (v4f32) __msa_splati_w((v4i32) src_a5, 0);
    src_a10 = __msa_cast_to_vector_float(*(a + 10));
    src_a10 = (v4f32) __msa_splati_w((v4i32) src_a10, 0);
    src_a11 = __msa_cast_to_vector_float(*(a + 11));
    src_a11 = (v4f32) __msa_splati_w((v4i32) src_a11, 0);
    src_a15 = __msa_cast_to_vector_float(*(a + 15));
    src_a15 = (v4f32) __msa_splati_w((v4i32) src_a15, 0);

    res_c0 *= src_a0;
    res_c4 *= src_a0;
    res_c1 -= res_c0 * src_a1;
    res_c5 -= res_c4 * src_a1;
    res_c2 -= res_c0 * src_a2;
    res_c6 -= res_c4 * src_a2;
    res_c3 -= res_c0 * src_a3;
    res_c7 -= res_c4 * src_a3;

    res_c1 *= src_a5;
    res_c5 *= src_a5;
    res_c2 -= res_c1 * src_a6;
    res_c6 -= res_c5 * src_a6;
    res_c3 -= res_c1 * src_a7;
    res_c7 -= res_c5 * src_a7;

    res_c2 *= src_a10;
    res_c6 *= src_a10;
    res_c3 -= res_c2 * src_a11;
    res_c7 -= res_c6 * src_a11;

    res_c3 *= src_a15;
    res_c7 *= src_a15;

    ST_SP4(res_c0, res_c4, res_c1, res_c5, b, 4);
    ST_SP4(res_c2, res_c6, res_c3, res_c7, b + 16, 4);

    TRANSPOSE4x4_SP_SP(res_c0, res_c1, res_c2, res_c3,
                       src_c0, src_c1, src_c2, src_c3);
    TRANSPOSE4x4_SP_SP(res_c4, res_c5, res_c6, res_c7,
                       src_c4, src_c5, src_c6, src_c7);

    ST_SP(src_c0, c);
    ST_SP(src_c1, c_nxt1line);
    ST_SP(src_c2, c_nxt2line);
    ST_SP(src_c3, c_nxt3line);
    ST_SP(src_c4, c_nxt4line);
    ST_SP(src_c5, c_nxt5line);
    ST_SP(src_c6, c_nxt6line);
    ST_SP(src_c7, c_nxt7line);
}

static void ssolve_4x4_lt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG bk)
{
    v4f32 src_c0, src_c1, src_c2, src_c3, res_c0, res_c1, res_c2, res_c3;
    v4f32 src_a0, src_a1, src_a2, src_a3, src_a5, src_a6, src_a7;
    v4f32 src_a10, src_a11, src_a15, src_a;
    FLOAT *c_nxt1line = c + ldc;
    FLOAT *c_nxt2line = c + 2 * ldc;
    FLOAT *c_nxt3line = c + 3 * ldc;

    if (bk > 0)
    {
        BLASLONG k;
        v4f32 src_b, src_b0, src_b1, src_b2, src_b3;
        v4f32 res0, res1, res2, res3;

        src_a0 = LD_SP(a);

        src_b = LD_SP(b + 0);
        SPLATI_W4_SP(src_b, src_b0, src_b1, src_b2, src_b3);
        res0 = src_a0 * src_b0;
        res1 = src_a0 * src_b1;
        res2 = src_a0 * src_b2;
        res3 = src_a0 * src_b3;

        a += 4;
        b += 4;

        for (k = (bk - 1) >> 1; k--;)
        {
            src_a0 = LD_SP(a);

            src_b = LD_SP(b + 0);
            SPLATI_W4_SP(src_b, src_b0, src_b1, src_b2, src_b3);
            res0 += src_a0 * src_b0;
            res1 += src_a0 * src_b1;
            res2 += src_a0 * src_b2;
            res3 += src_a0 * src_b3;

            a += 4;
            b += 4;

            src_a0 = LD_SP(a);

            src_b = LD_SP(b + 0);
            SPLATI_W4_SP(src_b, src_b0, src_b1, src_b2, src_b3);
            res0 += src_a0 * src_b0;
            res1 += src_a0 * src_b1;
            res2 += src_a0 * src_b2;
            res3 += src_a0 * src_b3;

            a += 4;
            b += 4;
        }

        if ((bk - 1) & 1)
        {
            src_a0 = LD_SP(a);

            src_b = LD_SP(b + 0);
            SPLATI_W4_SP(src_b, src_b0, src_b1, src_b2, src_b3);
            res0 += src_a0 * src_b0;
            res1 += src_a0 * src_b1;
            res2 += src_a0 * src_b2;
            res3 += src_a0 * src_b3;

            a += 4;
            b += 4;
        }

        src_c0 = LD_SP(c);
        src_c1 = LD_SP(c_nxt1line);
        src_c2 = LD_SP(c_nxt2line);
        src_c3 = LD_SP(c_nxt3line);

        src_c0 -= res0;
        src_c1 -= res1;
        src_c2 -= res2;
        src_c3 -= res3;
    }
    else
    {
        src_c0 = LD_SP(c);
        src_c1 = LD_SP(c_nxt1line);
        src_c2 = LD_SP(c_nxt2line);
        src_c3 = LD_SP(c_nxt3line);
    }

    TRANSPOSE4x4_SP_SP(src_c0, src_c1, src_c2, src_c3,
                       res_c0, res_c1, res_c2, res_c3);

    src_a = LD_SP(a + 0);
    SPLATI_W4_SP(src_a, src_a0, src_a1, src_a2, src_a3);
    src_a5 = LD_SP(a + 5);
    src_a7 = (v4f32) __msa_splati_w((v4i32) src_a5, 2);
    src_a6 = (v4f32) __msa_splati_w((v4i32) src_a5, 1);
    src_a5 = (v4f32) __msa_splati_w((v4i32) src_a5, 0);
    src_a10 = __msa_cast_to_vector_float(*(a + 10));
    src_a10 = (v4f32) __msa_splati_w((v4i32) src_a10, 0);
    src_a11 = __msa_cast_to_vector_float(*(a + 11));
    src_a11 = (v4f32) __msa_splati_w((v4i32) src_a11, 0);
    src_a15 = __msa_cast_to_vector_float(*(a + 15));
    src_a15 = (v4f32) __msa_splati_w((v4i32) src_a15, 0);

    res_c0 *= src_a0;
    res_c1 -= res_c0 * src_a1;
    res_c2 -= res_c0 * src_a2;
    res_c3 -= res_c0 * src_a3;

    res_c1 *= src_a5;
    res_c2 -= res_c1 * src_a6;
    res_c3 -= res_c1 * src_a7;

    res_c2 *= src_a10;
    res_c3 -= res_c2 * src_a11;

    res_c3 *= src_a15;

    ST_SP4(res_c0, res_c1, res_c2, res_c3, b, 4);

    TRANSPOSE4x4_SP_SP(res_c0, res_c1, res_c2, res_c3,
                       src_c0, src_c1, src_c2, src_c3);

    ST_SP(src_c0, c);
    ST_SP(src_c1, c_nxt1line);
    ST_SP(src_c2, c_nxt2line);
    ST_SP(src_c3, c_nxt3line);
}

static void ssolve_4x2_lt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG bk)
{
    FLOAT c0, c1, c2, c3, c0_nxt, c1_nxt, c2_nxt, c3_nxt;
    FLOAT a0, a1, a2, a3, a5, a6, a7, a10, a11, a15;

    c0 = *(c + 0);
    c1 = *(c + 1);
    c2 = *(c + 2);
    c3 = *(c + 3);
    c0_nxt = *(c + 0 + ldc);
    c1_nxt = *(c + 1 + ldc);
    c2_nxt = *(c + 2 + ldc);
    c3_nxt = *(c + 3 + ldc);

    if (bk > 0)
    {
        BLASLONG k;
        FLOAT res[8];

        res[0] = a[0] * b[0];
        res[1] = a[1] * b[0];
        res[2] = a[2] * b[0];
        res[3] = a[3] * b[0];
        res[4] = a[0] * b[1];
        res[5] = a[1] * b[1];
        res[6] = a[2] * b[1];
        res[7] = a[3] * b[1];

        for (k = (bk - 1); k--;)
        {
            a += 4;
            b += 2;

            res[0] += a[0] * b[0];
            res[1] += a[1] * b[0];
            res[2] += a[2] * b[0];
            res[3] += a[3] * b[0];
            res[4] += a[0] * b[1];
            res[5] += a[1] * b[1];
            res[6] += a[2] * b[1];
            res[7] += a[3] * b[1];
        }

        c0 -= res[0];
        c1 -= res[1];
        c2 -= res[2];
        c3 -= res[3];
        c0_nxt -= res[4];
        c1_nxt -= res[5];
        c2_nxt -= res[6];
        c3_nxt -= res[7];

        a += 4;
        b += 2;
    }

    a0 = *(a + 0);
    a1 = *(a + 1);
    a2 = *(a + 2);
    a3 = *(a + 3);
    a5 = *(a + 5);
    a6 = *(a + 6);
    a7 = *(a + 7);
    a10 = *(a + 10);
    a11 = *(a + 11);
    a15 = *(a + 15);

    c0 *= a0;
    c0_nxt *= a0;

    c1 -= c0 * a1;
    c1_nxt -= c0_nxt * a1;

    c1 *= a5;
    c1_nxt *= a5;

    c2 -= c0 * a2;
    c2_nxt -= c0_nxt * a2;

    c2 -= c1 * a6;
    c2_nxt -= c1_nxt * a6;

    c2 *= a10;
    c2_nxt *= a10;

    c3 -= c0 * a3;
    c3_nxt -= c0_nxt * a3;

    c3 -= c1 * a7;
    c3_nxt -= c1_nxt * a7;

    c3 -= c2 * a11;
    c3_nxt -= c2_nxt * a11;

    c3 *= a15;
    c3_nxt *= a15;

    *(b + 0) = c0;
    *(b + 1) = c0_nxt;
    *(b + 2) = c1;
    *(b + 3) = c1_nxt;
    *(b + 4) = c2;
    *(b + 5) = c2_nxt;
    *(b + 6) = c3;
    *(b + 7) = c3_nxt;

    *(c + 0) = c0;
    *(c + 1) = c1;
    *(c + 2) = c2;
    *(c + 3) = c3;
    *(c + 0 + ldc) = c0_nxt;
    *(c + 1 + ldc) = c1_nxt;
    *(c + 2 + ldc) = c2_nxt;
    *(c + 3 + ldc) = c3_nxt;
}

static void ssolve_4x1_lt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG bk)
{
    FLOAT a0, a1, a2, a3, a5, a6, a7, a10, a11, a15, c0, c1, c2, c3;

    c0 = *(c + 0);
    c1 = *(c + 1);
    c2 = *(c + 2);
    c3 = *(c + 3);

    if (bk > 0)
    {
        BLASLONG k;
        FLOAT t0, t1, t2, t3;

        t0 = a[0] * b[0];
        t1 = a[1] * b[0];
        t2 = a[2] * b[0];
        t3 = a[3] * b[0];

        for (k = (bk - 1); k--;)
        {
            a += 4;
            b += 1;

            t0 += a[0] * b[0];
            t1 += a[1] * b[0];
            t2 += a[2] * b[0];
            t3 += a[3] * b[0];
        }

        c0 -= t0;
        c1 -= t1;
        c2 -= t2;
        c3 -= t3;

        a += 4;
        b += 1;
    }

    a0 = *(a + 0);
    a1 = *(a + 1);
    a2 = *(a + 2);
    a3 = *(a + 3);
    a5 = *(a + 5);
    a6 = *(a + 6);
    a7 = *(a + 7);
    a10 = *(a + 10);
    a11 = *(a + 11);
    a15 = *(a + 15);

    c0 *= a0;

    c1 -= c0 * a1;
    c1 *= a5;

    c2 -= c0 * a2;
    c2 -= c1 * a6;
    c2 *= a10;

    c3 -= c0 * a3;
    c3 -= c1 * a7;
    c3 -= c2 * a11;
    c3 *= a15;

    *(b + 0) = c0;
    *(b + 1) = c1;
    *(b + 2) = c2;
    *(b + 3) = c3;

    *(c + 0) = c0;
    *(c + 1) = c1;
    *(c + 2) = c2;
    *(c + 3) = c3;
}

static void ssolve_2x8_lt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG bk)
{
    FLOAT a0, a1, a3, c0, c1, c0_nxt1, c1_nxt1, c0_nxt2, c1_nxt2;
    FLOAT c0_nxt3, c1_nxt3, c0_nxt4, c1_nxt4, c0_nxt5, c1_nxt5;
    FLOAT c0_nxt6, c1_nxt6, c0_nxt7, c1_nxt7;

    c0 = *(c + 0);
    c1 = *(c + 1);
    c0_nxt1 = *(c + ldc);
    c1_nxt1 = *(c + 1 + ldc);
    c0_nxt2 = *(c + 2 * ldc);
    c1_nxt2 = *(c + 1 + 2 * ldc);
    c0_nxt3 = *(c + 3 * ldc);
    c1_nxt3 = *(c + 1 + 3 * ldc);
    c0_nxt4 = *(c + 4 * ldc);
    c1_nxt4 = *(c + 1 + 4 * ldc);
    c0_nxt5 = *(c + 5 * ldc);
    c1_nxt5 = *(c + 1 + 5 * ldc);
    c0_nxt6 = *(c + 6 * ldc);
    c1_nxt6 = *(c + 1 + 6 * ldc);
    c0_nxt7 = *(c + 7 * ldc);
    c1_nxt7 = *(c + 1 + 7 * ldc);

    if (bk > 0)
    {
        BLASLONG k;
        FLOAT res[16];

        res[0] = a[0] * b[0];
        res[1] = a[1] * b[0];
        res[2] = a[0] * b[1];
        res[3] = a[1] * b[1];
        res[4] = a[0] * b[2];
        res[5] = a[1] * b[2];
        res[6] = a[0] * b[3];
        res[7] = a[1] * b[3];
        res[8] = a[0] * b[4];
        res[9] = a[1] * b[4];
        res[10] = a[0] * b[5];
        res[11] = a[1] * b[5];
        res[12] = a[0] * b[6];
        res[13] = a[1] * b[6];
        res[14] = a[0] * b[7];
        res[15] = a[1] * b[7];

        for (k = (bk - 1); k--;)
        {
            a += 2;
            b += 8;

            res[0] += a[0] * b[0];
            res[1] += a[1] * b[0];
            res[2] += a[0] * b[1];
            res[3] += a[1] * b[1];
            res[4] += a[0] * b[2];
            res[5] += a[1] * b[2];
            res[6] += a[0] * b[3];
            res[7] += a[1] * b[3];
            res[8] += a[0] * b[4];
            res[9] += a[1] * b[4];
            res[10] += a[0] * b[5];
            res[11] += a[1] * b[5];
            res[12] += a[0] * b[6];
            res[13] += a[1] * b[6];
            res[14] += a[0] * b[7];
            res[15] += a[1] * b[7];
        }

        c0 -= res[0];
        c1 -= res[1];
        c0_nxt1 -= res[2];
        c1_nxt1 -= res[3];
        c0_nxt2 -= res[4];
        c1_nxt2 -= res[5];
        c0_nxt3 -= res[6];
        c1_nxt3 -= res[7];
        c0_nxt4 -= res[8];
        c1_nxt4 -= res[9];
        c0_nxt5 -= res[10];
        c1_nxt5 -= res[11];
        c0_nxt6 -= res[12];
        c1_nxt6 -= res[13];
        c0_nxt7 -= res[14];
        c1_nxt7 -= res[15];

        a += 2;
        b += 8;
    }

    a0 = *a;
    a1 = *(a + 1);
    a3 = *(a + 3);

    c0 = c0 * a0;
    c1 = (c1 - c0 * a1) * a3;

    c0_nxt1 = c0_nxt1 * a0;
    c1_nxt1 = (c1_nxt1 - c0_nxt1 * a1) * a3;

    c0_nxt2 = c0_nxt2 * a0;
    c1_nxt2 = (c1_nxt2 - c0_nxt2 * a1) * a3;

    c0_nxt3 = c0_nxt3 * a0;
    c1_nxt3 = (c1_nxt3 - c0_nxt3 * a1) * a3;

    c0_nxt4 = c0_nxt4 * a0;
    c1_nxt4 = (c1_nxt4 - c0_nxt4 * a1) * a3;

    c0_nxt5 = c0_nxt5 * a0;
    c1_nxt5 = (c1_nxt5 - c0_nxt5 * a1) * a3;

    c0_nxt6 = c0_nxt6 * a0;
    c1_nxt6 = (c1_nxt6 - c0_nxt6 * a1) * a3;

    c0_nxt7 = c0_nxt7 * a0;
    c1_nxt7 = (c1_nxt7 - c0_nxt7 * a1) * a3;

    *(b + 0) = c0;
    *(b + 1) = c0_nxt1;
    *(b + 2) = c0_nxt2;
    *(b + 3) = c0_nxt3;
    *(b + 4) = c0_nxt4;
    *(b + 5) = c0_nxt5;
    *(b + 6) = c0_nxt6;
    *(b + 7) = c0_nxt7;
    *(b + 8) = c1;
    *(b + 9) = c1_nxt1;
    *(b + 10) = c1_nxt2;
    *(b + 11) = c1_nxt3;
    *(b + 12) = c1_nxt4;
    *(b + 13) = c1_nxt5;
    *(b + 14) = c1_nxt6;
    *(b + 15) = c1_nxt7;

    *(c + 0) = c0;
    *(c + 1) = c1;
    *(c + 0 + ldc) = c0_nxt1;
    *(c + 1 + ldc) = c1_nxt1;
    *(c + 0 + 2 * ldc) = c0_nxt2;
    *(c + 1 + 2 * ldc) = c1_nxt2;
    *(c + 0 + 3 * ldc) = c0_nxt3;
    *(c + 1 + 3 * ldc) = c1_nxt3;
    *(c + 0 + 4 * ldc) = c0_nxt4;
    *(c + 1 + 4 * ldc) = c1_nxt4;
    *(c + 0 + 5 * ldc) = c0_nxt5;
    *(c + 1 + 5 * ldc) = c1_nxt5;
    *(c + 0 + 6 * ldc) = c0_nxt6;
    *(c + 1 + 6 * ldc) = c1_nxt6;
    *(c + 0 + 7 * ldc) = c0_nxt7;
    *(c + 1 + 7 * ldc) = c1_nxt7;
}

static void ssolve_2x4_lt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG bk)
{
    FLOAT a0, a1, a3, c0, c1, c0_nxt1, c1_nxt1;
    FLOAT c0_nxt2, c1_nxt2, c0_nxt3, c1_nxt3;

    c0 = *(c + 0);
    c1 = *(c + 1);
    c0_nxt1 = *(c + ldc);
    c1_nxt1 = *(c + 1 + ldc);
    c0_nxt2 = *(c + 2 * ldc);
    c1_nxt2 = *(c + 1 + 2 * ldc);
    c0_nxt3 = *(c + 3 * ldc);
    c1_nxt3 = *(c + 1 + 3 * ldc);

    if (bk > 0)
    {
        BLASLONG k;
        FLOAT res[8];

        res[0] = a[0] * b[0];
        res[1] = a[1] * b[0];
        res[2] = a[0] * b[1];
        res[3] = a[1] * b[1];
        res[4] = a[0] * b[2];
        res[5] = a[1] * b[2];
        res[6] = a[0] * b[3];
        res[7] = a[1] * b[3];

        for (k = (bk - 1); k--;)
        {
            a += 2;
            b += 4;

            res[0] += a[0] * b[0];
            res[1] += a[1] * b[0];
            res[2] += a[0] * b[1];
            res[3] += a[1] * b[1];
            res[4] += a[0] * b[2];
            res[5] += a[1] * b[2];
            res[6] += a[0] * b[3];
            res[7] += a[1] * b[3];
        }

        c0 -= res[0];
        c1 -= res[1];
        c0_nxt1 -= res[2];
        c1_nxt1 -= res[3];
        c0_nxt2 -= res[4];
        c1_nxt2 -= res[5];
        c0_nxt3 -= res[6];
        c1_nxt3 -= res[7];

        a += 2;
        b += 4;
    }

    a0 = *a;
    a1 = *(a + 1);
    a3 = *(a + 3);

    c0 *= a0;
    c0_nxt1 *= a0;
    c0_nxt2 *= a0;
    c0_nxt3 *= a0;

    c1 -= c0 * a1;
    c1_nxt1 -= c0_nxt1 * a1;
    c1_nxt2 -= c0_nxt2 * a1;
    c1_nxt3 -= c0_nxt3 * a1;
    c1 *= a3;
    c1_nxt1 *= a3;
    c1_nxt2 *= a3;
    c1_nxt3 *= a3;

    *(b + 0) = c0;
    *(b + 1) = c0_nxt1;
    *(b + 2) = c0_nxt2;
    *(b + 3) = c0_nxt3;
    *(b + 4) = c1;
    *(b + 5) = c1_nxt1;
    *(b + 6) = c1_nxt2;
    *(b + 7) = c1_nxt3;

    *(c + 0) = c0;
    *(c + 1) = c1;
    *(c + 0 + ldc) = c0_nxt1;
    *(c + 1 + ldc) = c1_nxt1;
    *(c + 0 + 2 * ldc) = c0_nxt2;
    *(c + 1 + 2 * ldc) = c1_nxt2;
    *(c + 0 + 3 * ldc) = c0_nxt3;
    *(c + 1 + 3 * ldc) = c1_nxt3;
}

static void ssolve_2x2_lt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG bk)
{
    FLOAT a0, a1, a3, c0, c1, c0_nxt, c1_nxt;

    c0 = *(c + 0);
    c1 = *(c + 1);
    c0_nxt = *(c + ldc);
    c1_nxt = *(c + 1 + ldc);

    if (bk > 0)
    {
        BLASLONG k;
        FLOAT res0, res1, res2, res3;

        res0 = a[0] * b[0];
        res1 = a[1] * b[0];
        res2 = a[0] * b[1];
        res3 = a[1] * b[1];

        for (k = (bk - 1); k--;)
        {
            a += 2;
            b += 2;

            res0 += a[0] * b[0];
            res1 += a[1] * b[0];
            res2 += a[0] * b[1];
            res3 += a[1] * b[1];
        }

        c0 -= res0;
        c1 -= res1;

        c0_nxt -= res2;
        c1_nxt -= res3;

        a += 2;
        b += 2;
    }

    a0 = *a;
    a1 = *(a + 1);
    a3 = *(a + 3);

    c0 *= a0;
    c0_nxt *= a0;
    c1 -= c0 * a1;
    c1_nxt -= c0_nxt * a1;
    c1 *= a3;
    c1_nxt *= a3;

    *(b + 0) = c0;
    *(b + 1) = c0_nxt;
    *(b + 2) = c1;
    *(b + 3) = c1_nxt;

    *(c + 0) = c0;
    *(c + 1) = c1;
    *(c + 0 + ldc) = c0_nxt;
    *(c + 1 + ldc) = c1_nxt;
}

static void ssolve_2x1_lt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG bk)
{
    FLOAT c0, c1;

    c0 = *(c + 0);
    c1 = *(c + 1);

    if (bk > 0)
    {
        BLASLONG k;
        FLOAT res0, res1;

        res0 = a[0] * b[0];
        res1 = a[1] * b[0];

        for (k = (bk - 1); k--;)
        {
            a += 2;
            b += 1;

            res0 += a[0] * b[0];
            res1 += a[1] * b[0];
        }

        c0 -= res0;
        c1 -= res1;

        a += 2;
        b += 1;
    }

    c0 *= *(a + 0);

    c1 -= c0 * *(a + 1);
    c1 *= *(a + 3);

    *(b + 0) = c0;
    *(b + 1) = c1;

    *(c + 0) = c0;
    *(c + 1) = c1;
}

static void ssolve_1x8_lt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG bk)
{
    if (bk > 0)
    {
        BLASLONG k;
        FLOAT c0, c1, c2, c3, c4, c5, c6, c7;

        c0 = a[0] * b[0];
        c1 = a[0] * b[1];
        c2 = a[0] * b[2];
        c3 = a[0] * b[3];
        c4 = a[0] * b[4];
        c5 = a[0] * b[5];
        c6 = a[0] * b[6];
        c7 = a[0] * b[7];

        for (k = (bk - 1); k--;)
        {
            a += 1;
            b += 8;

            c0 += a[0] * b[0];
            c1 += a[0] * b[1];
            c2 += a[0] * b[2];
            c3 += a[0] * b[3];
            c4 += a[0] * b[4];
            c5 += a[0] * b[5];
            c6 += a[0] * b[6];
            c7 += a[0] * b[7];
        }

        *(c + 0 * ldc) -= c0;
        *(c + 1 * ldc) -= c1;
        *(c + 2 * ldc) -= c2;
        *(c + 3 * ldc) -= c3;
        *(c + 4 * ldc) -= c4;
        *(c + 5 * ldc) -= c5;
        *(c + 6 * ldc) -= c6;
        *(c + 7 * ldc) -= c7;

        a += 1;
        b += 8;
    }

    *c *= *a;
    *(c + ldc) *= *a;
    *(c + 2 * ldc) *= *a;
    *(c + 3 * ldc) *= *a;
    *(c + 4 * ldc) *= *a;
    *(c + 5 * ldc) *= *a;
    *(c + 6 * ldc) *= *a;
    *(c + 7 * ldc) *= *a;

    *b = *c;
    *(b + 1) = *(c + ldc);
    *(b + 2) = *(c + 2 * ldc);
    *(b + 3) = *(c + 3 * ldc);
    *(b + 4) = *(c + 4 * ldc);
    *(b + 5) = *(c + 5 * ldc);
    *(b + 6) = *(c + 6 * ldc);
    *(b + 7) = *(c + 7 * ldc);
}

static void ssolve_1x4_lt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG bk)
{
    FLOAT c0, c1, c2, c3;

    c0 = *(c + 0 * ldc);
    c1 = *(c + 1 * ldc);
    c2 = *(c + 2 * ldc);
    c3 = *(c + 3 * ldc);

    if (bk > 0)
    {
        BLASLONG k;
        FLOAT res0, res1, res2, res3;

        res0 = a[0] * b[0];
        res1 = a[0] * b[1];
        res2 = a[0] * b[2];
        res3 = a[0] * b[3];

        for (k = (bk - 1); k--;)
        {
            a += 1;
            b += 4;

            res0 += a[0] * b[0];
            res1 += a[0] * b[1];
            res2 += a[0] * b[2];
            res3 += a[0] * b[3];
        }

        c0 -= res0;
        c1 -= res1;
        c2 -= res2;
        c3 -= res3;
        a += 1;
        b += 4;
    }

    c0 *= *a;
    c1 *= *a;
    c2 *= *a;
    c3 *= *a;

    *c = c0;
    *(c + ldc) = c1;
    *(c + 2 * ldc) = c2;
    *(c + 3 * ldc) = c3;

    *b = *c;
    *(b + 1) = *(c + ldc);
    *(b + 2) = *(c + 2 * ldc);
    *(b + 3) = *(c + 3 * ldc);
}

static void ssolve_1x2_lt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG bk)
{
    FLOAT c0, c1;

    c0 = *c;
    c1 = *(c + ldc);

    if (bk > 0)
    {
        BLASLONG k;
        FLOAT res0, res1;

        res0 = a[0] * b[0];
        res1 = a[0] * b[1];

        for (k = (bk - 1); k--;)
        {
            a += 1;
            b += 2;

            res0 += a[0] * b[0];
            res1 += a[0] * b[1];
        }

        c0 -= res0;
        c1 -= res1;

        a += 1;
        b += 2;
    }

    *c = c0 * *a;
    *(c + ldc) = c1 * *a;

    *b = *c;
    *(b + 1) = *(c + ldc);
}

static void ssolve_1x1_lt_msa(FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG bk)
{
    if (bk)
    {
        BLASLONG k;
        FLOAT res;

        res = a[0] * b[0];

        for (k = (bk - 1); k--;)
        {
            a++;
            b++;

            res += a[0] * b[0];
        }

        *c -= res;

        a++;
        b++;
    }

    *c *= *a;
    *b = *c;
}

int CNAME(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT dummy1, FLOAT *a, FLOAT *b,
          FLOAT *c, BLASLONG ldc, BLASLONG offset)
{
    FLOAT *aa, *cc;
    BLASLONG i, j, kk;

    for (j = (n >> 3); j--;)
    {
        kk = offset;
        aa = a;
        cc = c;

        for (i = (m >> 3); i--;)
        {
            ssolve_8x8_lt_msa(aa, b, cc, ldc, kk);

            aa += 8 * k;
            cc += 8;
            kk += 8;
        }

        if (m & 7)
        {
            if (m & 4)
            {
                ssolve_4x8_lt_msa(aa, b, cc, ldc, kk);

                aa += 4 * k;
                cc += 4;
                kk += 4;
            }

            if (m & 2)
            {
                ssolve_2x8_lt_msa(aa, b, cc, ldc, kk);

                aa += 2 * k;
                cc += 2;
                kk += 2;
            }

            if (m & 1)
            {
                ssolve_1x8_lt_msa(aa, b, cc, ldc, kk);

                aa += k;
                cc += 1;
                kk += 1;
            }
        }

        b += 8 * k;
        c += 8 * ldc;
    }

    if (n & 7)
    {
        if (n & 4)
        {
            kk = offset;
            aa = a;
            cc = c;

            for (i = (m >> 3); i--;)
            {
                ssolve_8x4_lt_msa(aa, b, cc, ldc, kk);

                aa += 8 * k;
                cc += 8;
                kk += 8;
            }

            if (m & 7)
            {
                if (m & 4)
                {
                    ssolve_4x4_lt_msa(aa, b, cc, ldc, kk);

                    aa += 4 * k;
                    cc += 4;
                    kk += 4;
                }

                if (m & 2)
                {
                    ssolve_2x4_lt_msa(aa, b, cc, ldc, kk);

                    aa += 2 * k;
                    cc += 2;
                    kk += 2;
                }

                if (m & 1)
                {
                    ssolve_1x4_lt_msa(aa, b, cc, ldc, kk);

                    aa += k;
                    cc += 1;
                    kk += 1;
                }
            }

            b += 4 * k;
            c += 4 * ldc;
        }

        if (n & 2)
        {
            kk = offset;
            aa = a;
            cc = c;

            for (i = (m >> 3); i--;)
            {
                ssolve_8x2_lt_msa(aa, b, cc, ldc, kk);

                aa += 8 * k;
                cc += 8;
                kk += 8;
            }

            if (m & 7)
            {
                if (m & 4)
                {
                    ssolve_4x2_lt_msa(aa, b, cc, ldc, kk);

                    aa += 4 * k;
                    cc += 4;
                    kk += 4;
                }

                if (m & 2)
                {
                    ssolve_2x2_lt_msa(aa, b, cc, ldc, kk);

                    aa += 2 * k;
                    cc += 2;
                    kk += 2;
                }

                if (m & 1)
                {
                    ssolve_1x2_lt_msa(aa, b, cc, ldc, kk);

                    aa += k;
                    cc += 1;
                    kk += 1;
                }
            }

            b += 2 * k;
            c += 2 * ldc;
        }

        if (n & 1)
        {
            kk = offset;
            aa = a;
            cc = c;

            for (i = (m >> 3); i--;)
            {
                ssolve_8x1_lt_msa(aa, b, cc, kk);

                aa += 8 * k;
                cc += 8;
                kk += 8;
            }

            if (m & 7)
            {
                if (m & 4)
                {
                    ssolve_4x1_lt_msa(aa, b, cc, kk);

                    aa += 4 * k;
                    cc += 4;
                    kk += 4;
                }

                if (m & 2)
                {
                    ssolve_2x1_lt_msa(aa, b, cc, kk);

                    aa += 2 * k;
                    cc += 2;
                    kk += 2;
                }

                if (m & 1)
                {
                    ssolve_1x1_lt_msa(aa, b, cc, kk);

                    aa += k;
                    cc += 1;
                    kk += 1;
                }
            }

            b += k;
            c += ldc;
        }
    }

    return 0;
}
