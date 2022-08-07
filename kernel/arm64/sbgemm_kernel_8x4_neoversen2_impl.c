/***************************************************************************
 * Copyright (c) 2022, The OpenBLAS Project
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in
 * the documentation and/or other materials provided with the
 * distribution.
 * 3. Neither the name of the OpenBLAS project nor the names of
 * its contributors may be used to endorse or promote products
 * derived from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * *****************************************************************************/

#include <arm_sve.h>

#include "common.h"

#ifdef ALPHA_ONE
#define LOAD_C(M, N) \
    mc##M##N = svld1_gather_index(pg32, ptr_c0##N + 2 * M , off_vc);

#define LOAD_C_LOW(M, N) \
    mc##M##N = svld1_gather_index(pg32_low, ptr_c0##N + 2 * M, off_vc);

#define LOAD_C_EVEN(M, N) \
    mc##M##N = svld1_gather_index(pg32_even, ptr_c0##N + 2 * M, off_vc);

#define LOAD_C_FIRST(M, N) \
    mc##M##N = svld1_gather_index(pg32_first, ptr_c0##N + 2 * M, off_vc);

#define STORE_C(M, N) \
    svst1_scatter_index(pg32, ptr_c0##N + 2 * M, off_vc, mc##M##N);

#define STORE_C_LOW(M, N) \
    svst1_scatter_index(pg32_low, ptr_c0##N + 2 * M, off_vc, mc##M##N);

#define STORE_C_EVEN(M, N) \
    svst1_scatter_index(pg32_even, ptr_c0##N + 2 * M, off_vc, mc##M##N);

#define STORE_C_FIRST(M, N) \
    svst1_scatter_index(pg32_first, ptr_c0##N + 2 * M, off_vc, mc##M##N);

#else
#define LOAD_C(M, N)  \
    mc##M##N = svdup_f32(0);  \
    oc##M##N = svld1_gather_index(pg32, ptr_c0##N + 2 * M , off_vc);

#define LOAD_C_LOW(M, N)  \
    mc##M##N = svdup_f32(0);  \
    oc##M##N = svld1_gather_index(pg32_low, ptr_c0##N + 2 * M , off_vc);

#define LOAD_C_EVEN(M, N)  \
    mc##M##N = svdup_f32(0);  \
    oc##M##N = svld1_gather_index(pg32_even, ptr_c0##N + 2 * M , off_vc);

#define LOAD_C_FIRST(M, N)  \
    mc##M##N = svdup_f32(0);  \
    oc##M##N = svld1_gather_index(pg32_first, ptr_c0##N + 2 * M , off_vc);

#define STORE_C(M, N) \
    mc##M##N = svmad_z(pg32, svalpha, mc##M##N, oc##M##N);  \
    svst1_scatter_index(pg32, ptr_c0##N + 2 * M, off_vc, mc##M##N);

#define STORE_C_LOW(M, N) \
    mc##M##N = svmad_z(pg32_low, svalpha, mc##M##N, oc##M##N);  \
    svst1_scatter_index(pg32_low, ptr_c0##N + 2 * M, off_vc, mc##M##N);

#define STORE_C_EVEN(M, N) \
    mc##M##N = svmad_z(pg32_even, svalpha, mc##M##N, oc##M##N);  \
    svst1_scatter_index(pg32_even, ptr_c0##N + 2 * M, off_vc, mc##M##N);

#define STORE_C_FIRST(M, N) \
    mc##M##N = svmad_z(pg32_first, svalpha, mc##M##N, oc##M##N);  \
    svst1_scatter_index(pg32_first, ptr_c0##N + 2 * M, off_vc, mc##M##N);

#endif

#define LOAD_A(M) ma##M = svld1_bf16(pg16, ptr_a##M);

#define LOAD_B(N) mb##N = svld1_bf16(pg16, ptr_b##N);

#define MATMUL(M, N) mc##M##N = svbfmmla(mc##M##N, ma##M, mb##N);

#define LOAD_KREST_1(NAME, M)                                    \
    m##NAME##M = svdupq_bf16(*(ptr_##NAME##M), zero, zero, zero, \
                             *(ptr_##NAME##M + 1), zero, zero, zero);

#define LOAD_KREST_1_LOW(NAME, M)                                            \
    m##NAME##M = svdupq_bf16(*(ptr_##NAME##M), zero, zero, zero, zero, zero, \
                             zero, zero);

#define LOAD_KREST_2(NAME, M)                                           \
    m##NAME##M =                                                        \
        svdupq_bf16(*(ptr_##NAME##M), *(ptr_##NAME##M + 1), zero, zero, \
                    *(ptr_##NAME##M + 2), *(ptr_##NAME##M + 3), zero, zero);

#define LOAD_KREST_2_LOW(NAME, M)                                          \
    m##NAME##M = svdupq_bf16(*(ptr_##NAME##M), *(ptr_##NAME##M + 1), zero, \
                             zero, zero, zero, zero, zero);

#define LOAD_KREST_3(NAME, M)                                         \
    m##NAME##M =                                                      \
        svdupq_bf16(*(ptr_##NAME##M), *(ptr_##NAME##M + 1),           \
                    *(ptr_##NAME##M + 2), zero, *(ptr_##NAME##M + 3), \
                    *(ptr_##NAME##M + 4), *(ptr_##NAME##M + 5), zero);

#define LOAD_KREST_3_LOW(NAME, M)                           \
    m##NAME##M =                                            \
        svdupq_bf16(*(ptr_##NAME##M), *(ptr_##NAME##M + 1), \
                    *(ptr_##NAME##M + 2), zero, zero, zero, zero, zero);


#ifdef ALPHA_ONE
int sbgemm_kernel_neoversen2_alpha_one(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, IFLOAT * A, IFLOAT * B, FLOAT * C, BLASLONG ldc)
#else
int sbgemm_kernel_neoversen2_alpha(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, IFLOAT * A, IFLOAT * B, FLOAT * C, BLASLONG ldc)
#endif
{
    bfloat16_t *ptr_a = (bfloat16_t *)A;
    bfloat16_t *ptr_b = (bfloat16_t *)B;
    FLOAT *ptr_c = C;

    bfloat16_t *ptr_a0, *ptr_a1, *ptr_a2, *ptr_a3;
    bfloat16_t *ptr_b0, *ptr_b1;
    FLOAT *ptr_c00, *ptr_c01;

    svbfloat16_t ma0, ma1, ma2, ma3, mb0, mb1;
    svfloat32_t mc00, mc01, mc10, mc11, mc20, mc21, mc30, mc31;
#ifndef ALPHA_ONE
    svfloat32_t oc00, oc01, oc10, oc11, oc20, oc21, oc30, oc31;
#endif
    svbool_t pg16 = svptrue_b16();
    svbool_t pg16_low = svdupq_b16(1, 1, 1, 1, 0, 0, 0, 0);
    svbool_t pg32 = svptrue_b32();
    svbool_t pg32_low = svdupq_b32(1, 1, 0, 0);
    svbool_t pg32_even = svdupq_b32(1, 0, 1, 0);
    svbool_t pg32_first = svdupq_b32(1, 0, 0, 0);
    svfloat32_t svalpha = svdup_f32(alpha);
    bfloat16 tmp = 0;
    bfloat16_t zero = *((bfloat16_t *)&tmp);
    BLASLONG krest = k & 3;

    // 00 01 10 11
    svuint32_t off_vc = svdupq_u32(0, (uint32_t)ldc, 1, (uint32_t)ldc + 1);

    for (BLASLONG j = 0; j < n / 4; j++) {
        ptr_c00 = ptr_c;
        ptr_c01 = ptr_c + 2 * ldc;
        ptr_c += 4 * ldc;

        ptr_a = (bfloat16_t *)A;

        for (BLASLONG i = 0; i < m / 8; i++) {
            ptr_a0 = ptr_a;
            ptr_a1 = ptr_a0 + 2 * k;
            ptr_a2 = ptr_a1 + 2 * k;
            ptr_a3 = ptr_a2 + 2 * k;
            ptr_a += 8 * k;

            ptr_b0 = ptr_b;
            ptr_b1 = ptr_b0 + 2 * k;

            LOAD_C(0, 0); LOAD_C(0, 1);
            LOAD_C(1, 0); LOAD_C(1, 1);
            LOAD_C(2, 0); LOAD_C(2, 1);
            LOAD_C(3, 0); LOAD_C(3, 1);

            for (BLASLONG p = 0; p < k / 4; p++) {
                LOAD_A(0); LOAD_A(1); LOAD_A(2); LOAD_A(3);
                LOAD_B(0); LOAD_B(1);

                MATMUL(0, 0); MATMUL(0, 1);
                MATMUL(1, 0); MATMUL(1, 1);
                MATMUL(2, 0); MATMUL(2, 1);
                MATMUL(3, 0); MATMUL(3, 1);

                ptr_a0 += 8; ptr_a1 += 8; ptr_a2 += 8; ptr_a3 += 8;
                ptr_b0 += 8; ptr_b1 += 8;
            }

            if (krest) {
                if (krest == 1) {
                    LOAD_KREST_1(a, 0); LOAD_KREST_1(a, 1);
                    LOAD_KREST_1(a, 2); LOAD_KREST_1(a, 3);
                    LOAD_KREST_1(b, 0); LOAD_KREST_1(b, 1);
                } else if (krest == 2) {
                    LOAD_KREST_2(a, 0); LOAD_KREST_2(a, 1);
                    LOAD_KREST_2(a, 2); LOAD_KREST_2(a, 3);
                    LOAD_KREST_2(b, 0); LOAD_KREST_2(b, 1);
                } else if (krest == 3) {
                    LOAD_KREST_3(a, 0); LOAD_KREST_3(a, 1);
                    LOAD_KREST_3(a, 2); LOAD_KREST_3(a, 3);
                    LOAD_KREST_3(b, 0); LOAD_KREST_3(b, 1);
                }
                MATMUL(0, 0); MATMUL(0, 1);
                MATMUL(1, 0); MATMUL(1, 1);
                MATMUL(2, 0); MATMUL(2, 1);
                MATMUL(3, 0); MATMUL(3, 1);
            }

            STORE_C(0, 0); STORE_C(0, 1);
            STORE_C(1, 0); STORE_C(1, 1);
            STORE_C(2, 0); STORE_C(2, 1);
            STORE_C(3, 0); STORE_C(3, 1);

            ptr_c00 += 8; ptr_c01 += 8;
        }

        if (m & 4) {
            ptr_a0 = ptr_a;
            ptr_a1 = ptr_a0 + 2 * k;
            ptr_a += 4 * k;

            ptr_b0 = ptr_b;
            ptr_b1 = ptr_b0 + 2 * k;

            LOAD_C(0, 0); LOAD_C(0, 1);
            LOAD_C(1, 0); LOAD_C(1, 1);

            for (BLASLONG p = 0; p < k / 4; p++) {
                LOAD_A(0); LOAD_A(1);
                LOAD_B(0); LOAD_B(1);

                MATMUL(0, 0); MATMUL(0, 1);
                MATMUL(1, 0); MATMUL(1, 1);

                ptr_a0 += 8; ptr_a1 += 8;
                ptr_b0 += 8; ptr_b1 += 8;
            }

            if (krest) {
                if (krest == 1) {
                    LOAD_KREST_1(a, 0); LOAD_KREST_1(a, 1);
                    LOAD_KREST_1(b, 0); LOAD_KREST_1(b, 1);
                } else if (krest == 2) {
                    LOAD_KREST_2(a, 0); LOAD_KREST_2(a, 1);
                    LOAD_KREST_2(b, 0); LOAD_KREST_2(b, 1);
                } else if (krest == 3) {
                    LOAD_KREST_3(a, 0); LOAD_KREST_3(a, 1);
                    LOAD_KREST_3(b, 0); LOAD_KREST_3(b, 1);
                }
                MATMUL(0, 0); MATMUL(0, 1);
                MATMUL(1, 0); MATMUL(1, 1);
            }
            
            STORE_C(0, 0); STORE_C(0, 1);
            STORE_C(1, 0); STORE_C(1, 1);

            ptr_c00 += 4; ptr_c01 += 4;
        }

        if (m & 2) {
            ptr_a0 = ptr_a;
            ptr_a += 2 * k;

            ptr_b0 = ptr_b;
            ptr_b1 = ptr_b0 + 2 * k;

            LOAD_C(0, 0); LOAD_C(0, 1);

            for (BLASLONG p = 0; p < k / 4; p++) {
                LOAD_A(0);
                LOAD_B(0); LOAD_B(1);

                MATMUL(0, 0); MATMUL(0, 1);

                ptr_a0 += 8;
                ptr_b0 += 8; ptr_b1 += 8;
            }

            if (krest) {
                if (krest == 1) {
                    LOAD_KREST_1(a, 0);
                    LOAD_KREST_1(b, 0); LOAD_KREST_1(b, 1);
                } else if (krest == 2) {
                    LOAD_KREST_2(a, 0);
                    LOAD_KREST_2(b, 0); LOAD_KREST_2(b, 1);
                } else if (krest == 3) {
                    LOAD_KREST_3(a, 0);
                    LOAD_KREST_3(b, 0); LOAD_KREST_3(b, 1);
                }
                MATMUL(0, 0); MATMUL(0, 1);
            }
            STORE_C(0, 0); STORE_C(0, 1);
            ptr_c00 += 2; ptr_c01 += 2;
        }

        if (m & 1) {
            ptr_a0 = ptr_a;

            ptr_b0 = ptr_b;
            ptr_b1 = ptr_b0 + 2 * k;

            LOAD_C_LOW(0, 0); LOAD_C_LOW(0, 1);

            for (BLASLONG p = 0; p < k / 4; p++) {
                ma0 = svld1_bf16(pg16_low, ptr_a0);
                LOAD_B(0); LOAD_B(1);

                MATMUL(0, 0); MATMUL(0, 1);

                ptr_a0 += 4;
                ptr_b0 += 8;
                ptr_b1 += 8;
            }

            if (krest) {
                if (krest == 1) {
                    LOAD_KREST_1_LOW(a, 0);
                    LOAD_KREST_1(b, 0); LOAD_KREST_1(b, 1);
                } else if (krest == 2) {
                    LOAD_KREST_2_LOW(a, 0);
                    LOAD_KREST_2(b, 0); LOAD_KREST_2(b, 1);
                } else if (krest == 3) {
                    LOAD_KREST_3_LOW(a, 0);
                    LOAD_KREST_3(b, 0); LOAD_KREST_3(b, 1);
                }
                MATMUL(0, 0); MATMUL(0, 1);
            }
            STORE_C_LOW(0, 0); STORE_C_LOW(0, 1);
        }

        ptr_b += 4 * k;
    }

    if (n & 2) {
        ptr_c00 = ptr_c;
        ptr_c += 2 * ldc;

        ptr_a = (bfloat16_t *)A;

        for (BLASLONG i = 0; i < m / 8; i++) {
            ptr_a0 = ptr_a;
            ptr_a1 = ptr_a0 + 2 * k;
            ptr_a2 = ptr_a1 + 2 * k;
            ptr_a3 = ptr_a2 + 2 * k;
            ptr_a += 8 * k;

            ptr_b0 = ptr_b;

            LOAD_C(0, 0);
            LOAD_C(1, 0);
            LOAD_C(2, 0);
            LOAD_C(3, 0);

            for (BLASLONG p = 0; p < k / 4; p++) {
                LOAD_A(0); LOAD_A(1); LOAD_A(2); LOAD_A(3);
                LOAD_B(0);

                MATMUL(0, 0);
                MATMUL(1, 0);
                MATMUL(2, 0);
                MATMUL(3, 0);

                ptr_a0 += 8; ptr_a1 += 8; ptr_a2 += 8; ptr_a3 += 8;
                ptr_b0 += 8;
            }
            if (krest) {
                if (krest == 1) {
                    LOAD_KREST_1(a, 0); LOAD_KREST_1(a, 1);
                    LOAD_KREST_1(a, 2); LOAD_KREST_1(a, 3);
                    LOAD_KREST_1(b, 0);
                } else if (krest == 2) {
                    LOAD_KREST_2(a, 0); LOAD_KREST_2(a, 1);
                    LOAD_KREST_2(a, 2); LOAD_KREST_2(a, 3);
                    LOAD_KREST_2(b, 0);
                } else if (krest == 3) {
                    LOAD_KREST_3(a, 0); LOAD_KREST_3(a, 1);
                    LOAD_KREST_3(a, 2); LOAD_KREST_3(a, 3);
                    LOAD_KREST_3(b, 0);
                }
                MATMUL(0, 0);
                MATMUL(1, 0);
                MATMUL(2, 0);
                MATMUL(3, 0);
            }
            
            STORE_C(0, 0);
            STORE_C(1, 0);
            STORE_C(2, 0);
            STORE_C(3, 0);

            ptr_c00 += 8;
        }

        if (m & 4) {
            ptr_a0 = ptr_a;
            ptr_a1 = ptr_a0 + 2 * k;
            ptr_a += 4 * k;

            ptr_b0 = ptr_b;

            LOAD_C(0, 0);
            LOAD_C(1, 0);

            for (BLASLONG p = 0; p < k / 4; p++) {
                LOAD_A(0); LOAD_A(1);
                LOAD_B(0);

                MATMUL(0, 0);
                MATMUL(1, 0);

                ptr_a0 += 8; ptr_a1 += 8;
                ptr_b0 += 8;
            }
            if (krest) {
                if (krest == 1) {
                    LOAD_KREST_1(a, 0); LOAD_KREST_1(a, 1);
                    LOAD_KREST_1(b, 0);
                } else if (krest == 2) {
                    LOAD_KREST_2(a, 0); LOAD_KREST_2(a, 1);
                    LOAD_KREST_2(b, 0);
                } else if (krest == 3) {
                    LOAD_KREST_3(a, 0); LOAD_KREST_3(a, 1);
                    LOAD_KREST_3(b, 0);
                }
                MATMUL(0, 0);
                MATMUL(1, 0);
            }
            STORE_C(0, 0)
            STORE_C(1, 0)

            ptr_c00 += 4;
        }
        
        if (m & 2) {
            ptr_a0 = ptr_a;
            ptr_a += 2 * k;
            ptr_b0 = ptr_b;

            LOAD_C(0, 0);
            for (BLASLONG p = 0; p < k / 4; p++) {
                LOAD_A(0);
                LOAD_B(0);
                MATMUL(0, 0);
                ptr_a0 += 8;
                ptr_b0 += 8;
            }
            if (krest) {
                if (krest == 1) {
                    LOAD_KREST_1(a, 0);
                    LOAD_KREST_1(b, 0);
                } else if (krest == 2) {
                    LOAD_KREST_2(a, 0);
                    LOAD_KREST_2(b, 0);
                } else if (krest == 3) {
                    LOAD_KREST_3(a, 0);
                    LOAD_KREST_3(b, 0);
                }
                MATMUL(0, 0);
            }
            STORE_C(0, 0);
            ptr_c00 += 2;
        }

        if (m & 1) {
            ptr_a0 = ptr_a;

            ptr_b0 = ptr_b;

            LOAD_C(0, 0);

            for (BLASLONG p = 0; p < k / 4; p++) {
                ma0 = svld1_bf16(pg16_low, ptr_a0);
                LOAD_B(0);
                MATMUL(0, 0);
                ptr_a0 += 4;
                ptr_b0 += 8;
            }
            if (krest) {
                if (krest == 1) {
                    LOAD_KREST_1_LOW(a, 0);
                    LOAD_KREST_1(b, 0);
                } else if (krest == 2) {
                    LOAD_KREST_2_LOW(a, 0);
                    LOAD_KREST_2(b, 0);
                } else if (krest == 3) {
                    LOAD_KREST_3_LOW(a, 0);
                    LOAD_KREST_3(b, 0);
                }
                MATMUL(0, 0);
            }
            STORE_C_LOW(0, 0);
        }
        
        ptr_b += 2 * k;
    }

    if (n & 1) {
        ptr_c00 = ptr_c;
        ptr_a = (bfloat16_t *) A;

        for (BLASLONG i = 0; i < m / 8; i++) {
            ptr_a0 = ptr_a;
            ptr_a1 = ptr_a0 + 2 * k;
            ptr_a2 = ptr_a1 + 2 * k;
            ptr_a3 = ptr_a2 + 2 * k;
            ptr_a += 8 * k;

            ptr_b0 = ptr_b;

            LOAD_C_EVEN(0, 0);
            LOAD_C_EVEN(1, 0);
            LOAD_C_EVEN(2, 0);
            LOAD_C_EVEN(3, 0);

            for (BLASLONG p = 0; p < k / 4; p++) {
                LOAD_A(0); LOAD_A(1); LOAD_A(2); LOAD_A(3);
                mb0 = svld1_bf16(pg16_low, ptr_b0);

                MATMUL(0, 0);
                MATMUL(1, 0);
                MATMUL(2, 0);
                MATMUL(3, 0);

                ptr_a0 += 8; ptr_a1 += 8; ptr_a2 += 8; ptr_a3 += 8;
                ptr_b0 += 4;
            }
            if (krest) {
                if (krest == 1) {
                    LOAD_KREST_1(a, 0); LOAD_KREST_1(a, 1);
                    LOAD_KREST_1(a, 2); LOAD_KREST_1(a, 3);
                    LOAD_KREST_1_LOW(b, 0);
                } else if (krest == 2) {
                    LOAD_KREST_2(a, 0); LOAD_KREST_2(a, 1);
                    LOAD_KREST_2(a, 2); LOAD_KREST_2(a, 3);
                    LOAD_KREST_2_LOW(b, 0);
                } else if (krest == 3) {
                    LOAD_KREST_3(a, 0); LOAD_KREST_3(a, 1);
                    LOAD_KREST_3(a, 2); LOAD_KREST_3(a, 3);
                    LOAD_KREST_3_LOW(b, 0);
                }
                MATMUL(0, 0);
                MATMUL(1, 0);
                MATMUL(2, 0);
                MATMUL(3, 0);
            }
            STORE_C_EVEN(0, 0)
            STORE_C_EVEN(1, 0);
            STORE_C_EVEN(2, 0);
            STORE_C_EVEN(3, 0);

            ptr_c00 += 8;
        }
        
        if (m & 4) {
            ptr_a0 = ptr_a;
            ptr_a1 = ptr_a0 + 2 * k;
            ptr_a += 4 * k;

            ptr_b0 = ptr_b;

            LOAD_C_EVEN(0, 0);
            LOAD_C_EVEN(1, 0);

            for (BLASLONG p = 0; p < k / 4; p++) {
                LOAD_A(0); LOAD_A(1);
                mb0 = svld1_bf16(pg16_low, ptr_b0);

                MATMUL(0, 0);
                MATMUL(1, 0);

                ptr_a0 += 8; ptr_a1 += 8;
                ptr_b0 += 4;
            }
            if (krest) {
                if (krest == 1) {
                    LOAD_KREST_1(a, 0); LOAD_KREST_1(a, 1);
                    LOAD_KREST_1_LOW(b, 0);
                } else if (krest == 2) {
                    LOAD_KREST_2(a, 0); LOAD_KREST_2(a, 1);
                    LOAD_KREST_2_LOW(b, 0);
                } else if (krest == 3) {
                    LOAD_KREST_3(a, 0); LOAD_KREST_3(a, 1);
                    LOAD_KREST_3_LOW(b, 0);
                }
                MATMUL(0, 0);
                MATMUL(1, 0);
            }
            STORE_C_EVEN(0, 0)
            STORE_C_EVEN(1, 0)

            ptr_c00 += 4;
        }

        if (m & 2) {
            ptr_a0 = ptr_a;
            ptr_a += 2 * k;

            ptr_b0 = ptr_b;
            
            LOAD_C_EVEN(0, 0);

            for (BLASLONG p = 0; p < k / 4; p++) {
                LOAD_A(0);
                mb0 = svld1_bf16(pg16_low, ptr_b0);

                MATMUL(0, 0);

                ptr_a0 += 8;
                ptr_b0 += 4;
            }
            if (krest) {
                if (krest == 1) {
                    LOAD_KREST_1(a, 0);
                    LOAD_KREST_1_LOW(b, 0);
                } else if (krest == 2) {
                    LOAD_KREST_2(a, 0);
                    LOAD_KREST_2_LOW(b, 0);
                } else if (krest == 3) {
                    LOAD_KREST_3(a, 0);
                    LOAD_KREST_3_LOW(b, 0);
                }
                MATMUL(0, 0);
            }
            STORE_C_EVEN(0, 0);
            ptr_c00 += 2;
        }
        if (m & 1) {
            ptr_a0 = ptr_a;
            ptr_b0 = ptr_b;
            LOAD_C_FIRST(0, 0);
            for (BLASLONG p = 0; p < k / 4; p++) {
                ma0 = svld1_bf16(pg16_low, ptr_a0);
                mb0 = svld1_bf16(pg16_low, ptr_b0);

                MATMUL(0, 0);

                ptr_a0 += 4;
                ptr_b0 += 4;
            }
            if (krest) {
                if (krest == 1) {
                    LOAD_KREST_1_LOW(a, 0);
                    LOAD_KREST_1_LOW(b, 0);
                } else if (krest == 2) {
                    LOAD_KREST_2_LOW(a, 0);
                    LOAD_KREST_2_LOW(b, 0);
                } else if (krest == 3) {
                    LOAD_KREST_3_LOW(a, 0);
                    LOAD_KREST_3_LOW(b, 0);
                }
                MATMUL(0, 0);
            }
            STORE_C_FIRST(0, 0);
        }
    }

    return 0;
}