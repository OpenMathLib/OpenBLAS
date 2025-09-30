/***************************************************************************
Copyright (c) 2025 The OpenBLAS Project
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

#include <arm_sve.h>
#include <arm_neon.h>

#define UPDATE_PTRSx2 \
    a_ptr1 = a_ptr0 + lda;

#define UPDATE_PTRSx4 \
    UPDATE_PTRSx2 \
    a_ptr2 = a_ptr1 + lda; \
    a_ptr3 = a_ptr2 + lda;

#define UPDATE_PTRSx8 \
    UPDATE_PTRSx4 \
    a_ptr4 = a_ptr3 + lda; \
    a_ptr5 = a_ptr4 + lda; \
    a_ptr6 = a_ptr5 + lda; \
    a_ptr7 = a_ptr6 + lda;

#define LANESx2(MACRO, offset) \
    MACRO(offset, 0) \
    MACRO(offset, 1)

#define LANESx4(MACRO, offset) \
    LANESx2(MACRO, offset) \
    MACRO(offset, 2) \
    MACRO(offset, 3)

#define LANESx8(MACRO, offset) \
    LANESx4(MACRO, offset) \
    MACRO(offset, 4) \
    MACRO(offset, 5) \
    MACRO(offset, 6) \
    MACRO(offset, 7)

#define LOAD_A_VEC(offset, vec) \
    svbfloat16_t a_vec ## offset ## vec = svld1(pg_full, &a_ptr ## vec[i] + offset * sve_size_bf16);

#define UPDATE_ACCUMULATORS_FROM_LANE(offset, lane) \
    acc ## offset ## 0 = svbfmlalb_lane(acc ## offset ## 0, a_vec ## offset ## lane, x_vec, lane); \
    acc ## offset ## 1 = svbfmlalt_lane(acc ## offset ## 1, a_vec ## offset ## lane, x_vec, lane);

#define INIT_ACCUMULATORS(offset) \
    svfloat32_t acc ## offset ## 0 = svdup_f32(0.0); \
    svfloat32_t acc ## offset ## 1 = svdup_f32(0.0);

#define UPDATE_ACCUMULATORS(offset) \
    acc ## offset ## 0 = svbfmlalb(acc ## offset ## 0, a_vec ## offset ## 0, x_vec); \
    acc ## offset ## 1 = svbfmlalt(acc ## offset ## 1, a_vec ## offset ## 0, x_vec);

#define STORE_ACCUMULATORS(offset) \
    svbfloat16_t acc ## offset ## 0_bf16 = svcvt_bf16_x(pg_full, acc ## offset ## 0); \
    svbfloat16_t acc ## offset ## 1_bf16 = svcvt_bf16_x(pg_full, acc ## offset ## 1); \
    svbfloat16_t combined ## offset = svtrn1(acc ## offset ## 0_bf16, acc ## offset ## 1_bf16); \
    svst1(pg_full, &y[i] + offset * sve_size_bf16, combined ## offset);

#define ALPHA_OP(offset) \
    acc ## offset ## 0 = svmul_x(pg_full, acc ## offset ## 0, svalpha); \
    acc ## offset ## 1 = svmul_x(pg_full, acc ## offset ## 1, svalpha);

#define BETA_OP(offset) \
    svbfloat16_t y_vec ## offset = svld1(pg_full, &y[i] + offset * sve_size_bf16); \
    acc ## offset ## 0 = svbfmlalb(acc ## offset ## 0, svbeta16, y_vec ## offset); \
    acc ## offset ## 1 = svbfmlalt(acc ## offset ## 1, svbeta16, y_vec ## offset);

int CNAME(BLASLONG m, BLASLONG n, FLOAT alpha, IFLOAT *a_in, BLASLONG lda, IFLOAT *x_in, BLASLONG inc_x, FLOAT beta, FLOAT *y_in, BLASLONG inc_y) 
{
    BLASLONG ix;
    bfloat16_t *a_ptr0, *a_ptr1, *a_ptr2, *a_ptr3, *a_ptr4, *a_ptr5, *a_ptr6, *a_ptr7;
    BLASLONG sve_size_bf16 = svcnth();
    BLASLONG sve_size2_bf16 = sve_size_bf16 * 2;
    BLASLONG sve_size3_bf16 = sve_size_bf16 * 3;
    svbool_t pg_full = svptrue_b16();
    svbool_t pg_tail = svwhilelt_b16_s64(0, m % sve_size_bf16);

    BLASLONG n8 = n & -8;
    BLASLONG n4 = n & -4;
    BLASLONG n2 = n & -2;
    
    bfloat16_t *a = (bfloat16_t*)a_in;
    bfloat16_t *x = (bfloat16_t*)x_in;
    bfloat16_t *y = (bfloat16_t*)y_in;

    bfloat16_t alpha_bf16, beta_bf16;
    memcpy(&alpha_bf16, &alpha, sizeof(bfloat16_t));
    memcpy(&beta_bf16, &beta, sizeof(bfloat16_t));
    float beta_f32 = vcvtah_f32_bf16(beta_bf16);
    float alpha_f32 = vcvtah_f32_bf16(alpha_bf16);
    svfloat32_t svalpha = svdup_f32(vcvtah_f32_bf16(alpha_bf16));
    svbfloat16_t svbeta16 = svdup_bf16(beta_bf16);
    
    BLASLONG i = 0;
    if (inc_y == 1) {
        for (; (i + sve_size3_bf16 - 1) < m; i+= sve_size3_bf16) {
            INIT_ACCUMULATORS(0);
            INIT_ACCUMULATORS(1);
            INIT_ACCUMULATORS(2);

            BLASLONG j = 0;
            ix = 0;
            a_ptr0 = a;
            UPDATE_PTRSx8;

            if (inc_x == 1) {
                for (; j < n4; j+= 4) {
                    uint64_t* x_u64 = (uint64_t*)&x[ix];
                    svbfloat16_t x_vec = svreinterpret_bf16_u64(svdup_u64(*x_u64));

                    LANESx4(LOAD_A_VEC, 0);
                    LANESx4(LOAD_A_VEC, 1);
                    LANESx4(LOAD_A_VEC, 2);
                    LANESx4(UPDATE_ACCUMULATORS_FROM_LANE, 0);
                    LANESx4(UPDATE_ACCUMULATORS_FROM_LANE, 1);
                    LANESx4(UPDATE_ACCUMULATORS_FROM_LANE, 2);

                    ix += 4;

                    a_ptr0 += 4 * lda;
                    UPDATE_PTRSx4;
                }

                for (; j < n2; j+= 2) {
                    uint32_t* x_u32 = (uint32_t*)&x[ix];
                    svbfloat16_t x_vec = svreinterpret_bf16_u32(svdup_u32(*x_u32));

                    LANESx2(LOAD_A_VEC, 0);
                    LANESx2(LOAD_A_VEC, 1);
                    LANESx2(LOAD_A_VEC, 2);
                    LANESx2(UPDATE_ACCUMULATORS_FROM_LANE, 0);
                    LANESx2(UPDATE_ACCUMULATORS_FROM_LANE, 1);
                    LANESx2(UPDATE_ACCUMULATORS_FROM_LANE, 2);

                    ix += 2;

                    a_ptr0 += 2 * lda;
                    UPDATE_PTRSx2;
                }
            }

            for (; j < n; j++) {
                svbfloat16_t x_vec = svdup_bf16(x[ix]);
                LOAD_A_VEC(0, 0);
                LOAD_A_VEC(1, 0);
                LOAD_A_VEC(2, 0);
                UPDATE_ACCUMULATORS(0);
                UPDATE_ACCUMULATORS(1);
                UPDATE_ACCUMULATORS(2);

                ix += inc_x;
                a_ptr0 += lda;
            }

            if (alpha_f32 != ONE) {
                ALPHA_OP(0);
                ALPHA_OP(1);
                ALPHA_OP(2);
            }
            if (beta_f32 != ZERO) {
                BETA_OP(0);
                BETA_OP(1);
                BETA_OP(2);
            }

            STORE_ACCUMULATORS(0);
            STORE_ACCUMULATORS(1);
            STORE_ACCUMULATORS(2);
        }

        for (; (i + sve_size_bf16 - 1) < m; i+= sve_size_bf16) {
            INIT_ACCUMULATORS(0);

            BLASLONG j = 0;
            ix = 0;
            a_ptr0 = a;
            UPDATE_PTRSx8;

            if (inc_x == 1) {
                for (; j < n8; j+= 8) {
                    svbfloat16_t x_vec = svld1rq(pg_full, &x[ix]);
                    LANESx8(LOAD_A_VEC, 0);
                    LANESx8(UPDATE_ACCUMULATORS_FROM_LANE, 0);

                    ix += 8;

                    a_ptr0 += 8 * lda;
                    UPDATE_PTRSx8;
                }

                for (; j < n4; j+= 4) {
                    uint64_t* x_u64 = (uint64_t*)&x[ix];
                    svbfloat16_t x_vec = svreinterpret_bf16_u64(svdup_u64(*x_u64));

                    LANESx4(LOAD_A_VEC, 0);
                    LANESx4(UPDATE_ACCUMULATORS_FROM_LANE, 0);

                    ix += 4;

                    a_ptr0 += 4 * lda;
                    UPDATE_PTRSx4;
                }

                for (; j < n2; j+= 2) {
                    uint32_t* x_u32 = (uint32_t*)&x[ix];
                    svbfloat16_t x_vec = svreinterpret_bf16_u32(svdup_u32(*x_u32));

                    LANESx2(LOAD_A_VEC, 0);
                    LANESx2(UPDATE_ACCUMULATORS_FROM_LANE, 0);

                    ix += 2;

                    a_ptr0 += 2 * lda;
                    UPDATE_PTRSx2;
                }
            }

            for (; j < n; j++) {
                svbfloat16_t x_vec = svdup_bf16(x[ix]);
                LOAD_A_VEC(0, 0);
                UPDATE_ACCUMULATORS(0);

                ix += inc_x;
                a_ptr0 += lda;
            }

            if (alpha_f32 != ONE) {
                ALPHA_OP(0);
            }
            if (beta_f32 != ZERO) {
                BETA_OP(0);
            }

            STORE_ACCUMULATORS(0);
        }

        if (i < m) {
            svfloat32_t acc0 = svdup_f32(0.0);
            svfloat32_t acc1 = svdup_f32(0.0);

            ix = 0;
            a_ptr0 = a;
            for (BLASLONG j = 0; j < n; j++) {
                svbfloat16_t x_vec0 = svdup_bf16(x[ix]);
                svbfloat16_t a_vec0 = svld1(pg_tail, &a_ptr0[i]);

                acc0 = svbfmlalb(acc0, a_vec0, x_vec0);
                acc1 = svbfmlalt(acc1, a_vec0, x_vec0);

                ix += inc_x;
                a_ptr0 += lda;
            }

            if (alpha_f32 != ONE) {
                acc0 = svmul_x(pg_full, acc0, svalpha);
                acc1 = svmul_x(pg_full, acc1, svalpha);
            }
            if (beta_f32 != ZERO) {
                svbfloat16_t y_vec = svld1(pg_tail, &y[i]);
                acc0 = svbfmlalb(acc0, svbeta16, y_vec);
                acc1 = svbfmlalt(acc1, svbeta16, y_vec);
            }

            svbfloat16_t acc0_bf16 = svcvt_bf16_x(pg_full, acc0);
            svbfloat16_t acc1_bf16 = svcvt_bf16_x(pg_full, acc1);
            svbfloat16_t combined = svtrn1(acc0_bf16, acc1_bf16);
            svst1(pg_tail, &y[i], combined);
        }

        return 0;
    }

    // Scalar fallback
    BLASLONG iy = 0;
    for (; i < m; i++) {
        float temp = 0.0;

        ix = 0;
        a_ptr0 = a;
        for (BLASLONG j = 0; j < n; j++) {
            temp += vcvtah_f32_bf16(a_ptr0[i]) * vcvtah_f32_bf16(x[ix]);
            ix += inc_x;
            a_ptr0 += lda;
        }

        if (beta_f32 == ZERO) {
            y[iy] = vcvth_bf16_f32(alpha_f32 * temp);
        } else {
            y[iy] = vcvth_bf16_f32(alpha_f32 * temp + beta_f32 * vcvtah_f32_bf16(y[iy]));
        }

        iy += inc_y;
    }

    return (0);
}
