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

#undef OP0
#undef OP1
#undef OP2
#undef OP3
#undef OP4

#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
    #define OP0  -=
    #define OP1  +=
    #define OP2  +=
#else
    #define OP0  +=
    #define OP1  +=
    #define OP2  -=
#endif

#define ZGEMV_T_4x4()                        \
    LD_DP4(pa0 + k, 2, t0, t1, t2, t3);      \
    LD_DP4(pa1 + k, 2, t4, t5, t6, t7);      \
    LD_DP4(pa2 + k, 2, t8, t9, t10, t11);    \
    LD_DP4(pa3 + k, 2, t12, t13, t14, t15);  \
                                             \
    PCKEVOD_D2_DP(t1, t0, src0r, src0i);     \
    PCKEVOD_D2_DP(t3, t2, src1r, src1i);     \
    PCKEVOD_D2_DP(t5, t4, src2r, src2i);     \
    PCKEVOD_D2_DP(t7, t6, src3r, src3i);     \
    PCKEVOD_D2_DP(t9, t8, src4r, src4i);     \
    PCKEVOD_D2_DP(t11, t10, src5r, src5i);   \
    PCKEVOD_D2_DP(t13, t12, src6r, src6i);   \
    PCKEVOD_D2_DP(t15, t14, src7r, src7i);   \
                                             \
    tp0r += src0r * x0r;                     \
    tp0r += src1r * x1r;                     \
    tp0r OP0 src0i * x0i;                    \
    tp0r OP0 src1i * x1i;                    \
                                             \
    tp1r += src2r * x0r;                     \
    tp1r += src3r * x1r;                     \
    tp1r OP0 src2i * x0i;                    \
    tp1r OP0 src3i * x1i;                    \
                                             \
    tp2r += src4r * x0r;                     \
    tp2r += src5r * x1r;                     \
    tp2r OP0 src4i * x0i;                    \
    tp2r OP0 src5i * x1i;                    \
                                             \
    tp3r += src6r * x0r;                     \
    tp3r += src7r * x1r;                     \
    tp3r OP0 src6i * x0i;                    \
    tp3r OP0 src7i * x1i;                    \
                                             \
    tp0i OP1 src0r * x0i;                    \
    tp0i OP1 src1r * x1i;                    \
    tp0i OP2 src0i * x0r;                    \
    tp0i OP2 src1i * x1r;                    \
                                             \
    tp1i OP1 src2r * x0i;                    \
    tp1i OP1 src3r * x1i;                    \
    tp1i OP2 src2i * x0r;                    \
    tp1i OP2 src3i * x1r;                    \
                                             \
    tp2i OP1 src4r * x0i;                    \
    tp2i OP1 src5r * x1i;                    \
    tp2i OP2 src4i * x0r;                    \
    tp2i OP2 src5i * x1r;                    \
                                             \
    tp3i OP1 src6r * x0i;                    \
    tp3i OP1 src7r * x1i;                    \
    tp3i OP2 src6i * x0r;                    \
    tp3i OP2 src7i * x1r;                    \

#define ZGEMV_T_4x2()                     \
    LD_DP4(pa0 + k, 2, t0, t1, t2, t3);   \
    LD_DP4(pa1 + k, 2, t4, t5, t6, t7);   \
                                          \
    PCKEVOD_D2_DP(t1, t0, src0r, src0i);  \
    PCKEVOD_D2_DP(t3, t2, src1r, src1i);  \
    PCKEVOD_D2_DP(t5, t4, src2r, src2i);  \
    PCKEVOD_D2_DP(t7, t6, src3r, src3i);  \
                                          \
    tp0r += src0r * x0r;                  \
    tp0r += src1r * x1r;                  \
    tp0r OP0 src0i * x0i;                 \
    tp0r OP0 src1i * x1i;                 \
                                          \
    tp1r += src2r * x0r;                  \
    tp1r += src3r * x1r;                  \
    tp1r OP0 src2i * x0i;                 \
    tp1r OP0 src3i * x1i;                 \
                                          \
    tp0i OP1 src0r * x0i;                 \
    tp0i OP1 src1r * x1i;                 \
    tp0i OP2 src0i * x0r;                 \
    tp0i OP2 src1i * x1r;                 \
                                          \
    tp1i OP1 src2r * x0i;                 \
    tp1i OP1 src3r * x1i;                 \
    tp1i OP2 src2i * x0r;                 \
    tp1i OP2 src3i * x1r;                 \

#define ZGEMV_T_4x1()                     \
    LD_DP4(pa0 + k, 2, t0, t1, t2, t3);   \
                                          \
    PCKEVOD_D2_DP(t1, t0, src0r, src0i);  \
    PCKEVOD_D2_DP(t3, t2, src1r, src1i);  \
                                          \
    tp0r += src0r * x0r;                  \
    tp0r += src1r * x1r;                  \
    tp0r OP0 src0i * x0i;                 \
    tp0r OP0 src1i * x1i;                 \
                                          \
    tp0i OP1 src0r * x0i;                 \
    tp0i OP1 src1r * x1i;                 \
    tp0i OP2 src0i * x0r;                 \
    tp0i OP2 src1i * x1r;                 \

#define ZGEMV_T_2x4()                       \
    LD_DP2(pa0 + k, 2, t0, t1);             \
    LD_DP2(pa1 + k, 2, t4, t5);             \
    LD_DP2(pa2 + k, 2, t8, t9);             \
    LD_DP2(pa3 + k, 2, t12, t13);           \
                                            \
    PCKEVOD_D2_DP(t1, t0, src0r, src0i);    \
    PCKEVOD_D2_DP(t5, t4, src2r, src2i);    \
    PCKEVOD_D2_DP(t9, t8, src4r, src4i);    \
    PCKEVOD_D2_DP(t13, t12, src6r, src6i);  \
                                            \
    tp0r += src0r * x0r;                    \
    tp0r OP0 src0i * x0i;                   \
                                            \
    tp1r += src2r * x0r;                    \
    tp1r OP0 src2i * x0i;                   \
                                            \
    tp2r += src4r * x0r;                    \
    tp2r OP0 src4i * x0i;                   \
                                            \
    tp3r += src6r * x0r;                    \
    tp3r OP0 src6i * x0i;                   \
                                            \
    tp0i OP1 src0r * x0i;                   \
    tp0i OP2 src0i * x0r;                   \
                                            \
    tp1i OP1 src2r * x0i;                   \
    tp1i OP2 src2i * x0r;                   \
                                            \
    tp2i OP1 src4r * x0i;                   \
    tp2i OP2 src4i * x0r;                   \
                                            \
    tp3i OP1 src6r * x0i;                   \
    tp3i OP2 src6i * x0r;                   \

#define ZGEMV_T_2x2()                     \
    LD_DP2(pa0 + k, 2, t0, t1);           \
    LD_DP2(pa1 + k, 2, t4, t5);           \
                                          \
    PCKEVOD_D2_DP(t1, t0, src0r, src0i);  \
    PCKEVOD_D2_DP(t5, t4, src2r, src2i);  \
                                          \
    tp0r += src0r * x0r;                  \
    tp0r OP0 src0i * x0i;                 \
                                          \
    tp1r += src2r * x0r;                  \
    tp1r OP0 src2i * x0i;                 \
                                          \
    tp0i OP1 src0r * x0i;                 \
    tp0i OP2 src0i * x0r;                 \
                                          \
    tp1i OP1 src2r * x0i;                 \
    tp1i OP2 src2i * x0r;                 \

#define ZGEMV_T_2x1()                     \
    LD_DP2(pa0 + k, 2, t0, t1);           \
                                          \
    PCKEVOD_D2_DP(t1, t0, src0r, src0i);  \
                                          \
    tp0r += src0r * x0r;                  \
    tp0r OP0 src0i * x0i;                 \
                                          \
    tp0i OP1 src0r * x0i;                 \
    tp0i OP2 src0i * x0r;                 \

#define ZGEMV_T_1x4()                           \
    temp0r  += pa0[k + 0] * x[0 * inc_x2];      \
    temp0r OP0 pa0[k + 1] * x[0 * inc_x2 + 1];  \
    temp1r  += pa1[k + 0] * x[0 * inc_x2];      \
    temp1r OP0 pa1[k + 1] * x[0 * inc_x2 + 1];  \
    temp2r  += pa2[k + 0] * x[0 * inc_x2];      \
    temp2r OP0 pa2[k + 1] * x[0 * inc_x2 + 1];  \
    temp3r  += pa3[k + 0] * x[0 * inc_x2];      \
    temp3r OP0 pa3[k + 1] * x[0 * inc_x2 + 1];  \
                                                \
    temp0i OP1 pa0[k + 0] * x[0 * inc_x2 + 1];  \
    temp0i OP2 pa0[k + 1] * x[0 * inc_x2];      \
    temp1i OP1 pa1[k + 0] * x[0 * inc_x2 + 1];  \
    temp1i OP2 pa1[k + 1] * x[0 * inc_x2];      \
    temp2i OP1 pa2[k + 0] * x[0 * inc_x2 + 1];  \
    temp2i OP2 pa2[k + 1] * x[0 * inc_x2];      \
    temp3i OP1 pa3[k + 0] * x[0 * inc_x2 + 1];  \
    temp3i OP2 pa3[k + 1] * x[0 * inc_x2];      \

#define ZGEMV_T_1x2()                           \
    temp0r  += pa0[k + 0] * x[0 * inc_x2];      \
    temp0r OP0 pa0[k + 1] * x[0 * inc_x2 + 1];  \
    temp1r  += pa1[k + 0] * x[0 * inc_x2];      \
    temp1r OP0 pa1[k + 1] * x[0 * inc_x2 + 1];  \
                                                \
    temp0i OP1 pa0[k + 0] * x[0 * inc_x2 + 1];  \
    temp0i OP2 pa0[k + 1] * x[0 * inc_x2];      \
    temp1i OP1 pa1[k + 0] * x[0 * inc_x2 + 1];  \
    temp1i OP2 pa1[k + 1] * x[0 * inc_x2];      \

#define ZGEMV_T_1x1()                           \
    temp0r  += pa0[k + 0] * x[0 * inc_x2];      \
    temp0r OP0 pa0[k + 1] * x[0 * inc_x2 + 1];  \
                                                \
    temp0i OP1 pa0[k + 0] * x[0 * inc_x2 + 1];  \
    temp0i OP2 pa0[k + 1] * x[0 * inc_x2];      \

#define ZSCALE_STORE_Y4_GP()    \
    res0r = y[0 * inc_y2];      \
    res1r = y[1 * inc_y2];      \
    res2r = y[2 * inc_y2];      \
    res3r = y[3 * inc_y2];      \
                                \
    res0i = y[0 * inc_y2 + 1];  \
    res1i = y[1 * inc_y2 + 1];  \
    res2i = y[2 * inc_y2 + 1];  \
    res3i = y[3 * inc_y2 + 1];  \
                                \
    res0r  += alphar * temp0r;  \
    res0r OP0 alphai * temp0i;  \
    res1r  += alphar * temp1r;  \
    res1r OP0 alphai * temp1i;  \
    res2r  += alphar * temp2r;  \
    res2r OP0 alphai * temp2i;  \
    res3r  += alphar * temp3r;  \
    res3r OP0 alphai * temp3i;  \
                                \
    res0i OP1 alphar * temp0i;  \
    res0i OP2 alphai * temp0r;  \
    res1i OP1 alphar * temp1i;  \
    res1i OP2 alphai * temp1r;  \
    res2i OP1 alphar * temp2i;  \
    res2i OP2 alphai * temp2r;  \
    res3i OP1 alphar * temp3i;  \
    res3i OP2 alphai * temp3r;  \
                                \
    y[0 * inc_y2] = res0r;      \
    y[1 * inc_y2] = res1r;      \
    y[2 * inc_y2] = res2r;      \
    y[3 * inc_y2] = res3r;      \
                                \
    y[0 * inc_y2 + 1] = res0i;  \
    y[1 * inc_y2 + 1] = res1i;  \
    y[2 * inc_y2 + 1] = res2i;  \
    y[3 * inc_y2 + 1] = res3i;  \

#define ZSCALE_STORE_Y2_GP()    \
    res0r = y[0 * inc_y2];      \
    res1r = y[1 * inc_y2];      \
                                \
    res0i = y[0 * inc_y2 + 1];  \
    res1i = y[1 * inc_y2 + 1];  \
                                \
    res0r  += alphar * temp0r;  \
    res0r OP0 alphai * temp0i;  \
    res1r  += alphar * temp1r;  \
    res1r OP0 alphai * temp1i;  \
                                \
    res0i OP1 alphar * temp0i;  \
    res0i OP2 alphai * temp0r;  \
    res1i OP1 alphar * temp1i;  \
    res1i OP2 alphai * temp1r;  \
                                \
    y[0 * inc_y2] = res0r;      \
    y[1 * inc_y2] = res1r;      \
                                \
    y[0 * inc_y2 + 1] = res0i;  \
    y[1 * inc_y2 + 1] = res1i;  \

#define ZSCALE_STORE_Y1_GP()    \
    res0r = y[0 * inc_y2];      \
    res0i = y[0 * inc_y2 + 1];  \
                                \
    res0r  += alphar * temp0r;  \
    res0r OP0 alphai * temp0i;  \
                                \
    res0i OP1 alphar * temp0i;  \
    res0i OP2 alphai * temp0r;  \
                                \
    y[0 * inc_y2] = res0r;      \
    y[0 * inc_y2 + 1] = res0i;  \

#define ZLOAD_X4_VECTOR()             \
    LD_DP4(x, 2, x0, x1, x2, x3);     \
    PCKEVOD_D2_DP(x1, x0, x0r, x0i);  \
    PCKEVOD_D2_DP(x3, x2, x1r, x1i);  \

#define ZLOAD_X2_VECTOR()             \
    LD_DP2(x, 2, x0, x1);             \
    PCKEVOD_D2_DP(x1, x0, x0r, x0i);  \

#define ZLOAD_X4_GP()                                                                      \
    x0r = (v2f64) __msa_insert_d((v2i64) tp0r, 0, *((long long *) (x + 0 * inc_x2)));      \
    x0r = (v2f64) __msa_insert_d((v2i64) x0r,  1, *((long long *) (x + 1 * inc_x2)));      \
    x1r = (v2f64) __msa_insert_d((v2i64) tp0r, 0, *((long long *) (x + 2 * inc_x2)));      \
    x1r = (v2f64) __msa_insert_d((v2i64) x1r,  1, *((long long *) (x + 3 * inc_x2)));      \
    x0i = (v2f64) __msa_insert_d((v2i64) tp0r, 0, *((long long *) (x + 0 * inc_x2 + 1)));  \
    x0i = (v2f64) __msa_insert_d((v2i64) x0i,  1, *((long long *) (x + 1 * inc_x2 + 1)));  \
    x1i = (v2f64) __msa_insert_d((v2i64) tp0r, 0, *((long long *) (x + 2 * inc_x2 + 1)));  \
    x1i = (v2f64) __msa_insert_d((v2i64) x1i,  1, *((long long *) (x + 3 * inc_x2 + 1)));  \

#define ZLOAD_X2_GP()                                                                      \
    x0r = (v2f64) __msa_insert_d((v2i64) tp0r, 0, *((long long *) (x + 0 * inc_x2)));      \
    x0r = (v2f64) __msa_insert_d((v2i64) x0r,  1, *((long long *) (x + 1 * inc_x2)));      \
    x0i = (v2f64) __msa_insert_d((v2i64) tp0r, 0, *((long long *) (x + 0 * inc_x2 + 1)));  \
    x0i = (v2f64) __msa_insert_d((v2i64) x0i,  1, *((long long *) (x + 1 * inc_x2 + 1)));  \

#define ZGEMV_T_MSA()                      \
    for (j = (n >> 2); j--;)               \
    {                                      \
        tp0r = tp1r = tp2r = tp3r = zero;  \
        tp0i = tp1i = tp2i = tp3i = zero;  \
                                           \
        k = 0;                             \
        x = srcx_org;                      \
                                           \
        for (i = (m >> 2); i--;)           \
        {                                  \
            ZLOAD_X4();                    \
            ZGEMV_T_4x4();                 \
                                           \
            k += 2 * 4;                    \
            x += inc_x2 * 4;               \
        }                                  \
                                           \
        if (m & 2)                         \
        {                                  \
            ZLOAD_X2();                    \
            ZGEMV_T_2x4();                 \
                                           \
            k += 2 * 2;                    \
            x += inc_x2 * 2;               \
        }                                  \
                                           \
        temp0r = tp0r[0] + tp0r[1];        \
        temp1r = tp1r[0] + tp1r[1];        \
        temp2r = tp2r[0] + tp2r[1];        \
        temp3r = tp3r[0] + tp3r[1];        \
        temp0i = tp0i[0] + tp0i[1];        \
        temp1i = tp1i[0] + tp1i[1];        \
        temp2i = tp2i[0] + tp2i[1];        \
        temp3i = tp3i[0] + tp3i[1];        \
                                           \
        if (m & 1)                         \
        {                                  \
            ZGEMV_T_1x4();                 \
                                           \
            k += 2;                        \
            x += inc_x2;                   \
        }                                  \
                                           \
        ZSCALE_STORE_Y4_GP();              \
                                           \
        pa0 += 4 * lda2;                   \
        pa1 += 4 * lda2;                   \
        pa2 += 4 * lda2;                   \
        pa3 += 4 * lda2;                   \
        y += 4 * inc_y2;                   \
    }                                      \
                                           \
    if (n & 2)                             \
    {                                      \
        tp0r = tp1r = zero;                \
        tp0i = tp1i = zero;                \
                                           \
        k = 0;                             \
        x = srcx_org;                      \
                                           \
        for (i = (m >> 2); i--;)           \
        {                                  \
            ZLOAD_X4();                    \
            ZGEMV_T_4x2();                 \
                                           \
            k += 2 * 4;                    \
            x += inc_x2 * 4;               \
        }                                  \
                                           \
        if (m & 2)                         \
        {                                  \
            ZLOAD_X2();                    \
            ZGEMV_T_2x2();                 \
                                           \
            k += 2 * 2;                    \
            x += inc_x2 * 2;               \
        }                                  \
                                           \
        temp0r = tp0r[0] + tp0r[1];        \
        temp1r = tp1r[0] + tp1r[1];        \
        temp0i = tp0i[0] + tp0i[1];        \
        temp1i = tp1i[0] + tp1i[1];        \
                                           \
        if (m & 1)                         \
        {                                  \
            ZGEMV_T_1x2();                 \
                                           \
            k += 2;                        \
            x += inc_x2;                   \
        }                                  \
                                           \
        ZSCALE_STORE_Y2_GP();              \
                                           \
        pa0 += 2 * lda2;                   \
        pa1 += 2 * lda2;                   \
        y += 2 * inc_y2;                   \
    }                                      \
                                           \
    if (n & 1)                             \
    {                                      \
        tp0r = zero;                       \
        tp0i = zero;                       \
                                           \
        k = 0;                             \
        x = srcx_org;                      \
                                           \
        for (i = (m >> 2); i--;)           \
        {                                  \
            ZLOAD_X4();                    \
            ZGEMV_T_4x1();                 \
                                           \
            k += 2 * 4;                    \
            x += inc_x2 * 4;               \
        }                                  \
                                           \
        if (m & 2)                         \
        {                                  \
            ZLOAD_X2();                    \
            ZGEMV_T_2x1();                 \
                                           \
            k += 2 * 2;                    \
            x += inc_x2 * 2;               \
        }                                  \
                                           \
        temp0r = tp0r[0] + tp0r[1];        \
        temp0i = tp0i[0] + tp0i[1];        \
                                           \
        if (m & 1)                         \
        {                                  \
            ZGEMV_T_1x1();                 \
                                           \
            k += 2;                        \
            x += inc_x2;                   \
        }                                  \
                                           \
        ZSCALE_STORE_Y1_GP();              \
                                           \
        pa0  += lda2;                      \
        y += inc_y2;                       \
    }                                      \

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alphar, FLOAT alphai,
          FLOAT *A, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y,
          BLASLONG inc_y, FLOAT *buffer)
{
    BLASLONG i, j, k;
    BLASLONG inc_x2, inc_y2, lda2;
    FLOAT *pa0, *pa1, *pa2, *pa3;
    FLOAT *srcx_org = x;
    FLOAT temp0r, temp0i, temp2r, temp2i, temp1r, temp1i, temp3r, temp3i;
    FLOAT res0r, res0i, res2r, res2i, res1r, res1i, res3r, res3i;
    v2f64 zero = {0};
    v2f64 x0, x1, x2, x3, x0r, x1r, x0i, x1i;
    v2f64 t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
    v2f64 src0r, src1r, src2r, src3r, src4r, src5r, src6r, src7r;
    v2f64 src0i, src1i, src2i, src3i, src4i, src5i, src6i, src7i;
    v2f64 tp0r, tp1r, tp2r, tp3r, tp0i, tp1i, tp2i, tp3i;

    lda2 = 2 * lda;

    pa0 = A;
    pa1 = A + lda2;
    pa2 = A + 2 * lda2;
    pa3 = A + 3 * lda2;

    inc_x2 = 2 * inc_x;
    inc_y2 = 2 * inc_y;

    if (2 == inc_x2)
    {
        #define ZLOAD_X4  ZLOAD_X4_VECTOR
        #define ZLOAD_X2  ZLOAD_X2_VECTOR

        ZGEMV_T_MSA();

        #undef ZLOAD_X4
        #undef ZLOAD_X2
    }
    else
    {
        #define ZLOAD_X4  ZLOAD_X4_GP
        #define ZLOAD_X2  ZLOAD_X2_GP

        ZGEMV_T_MSA();

        #undef ZLOAD_X4
        #undef ZLOAD_X2
    }

    return(0);
}

#undef OP0
#undef OP1
#undef OP2
