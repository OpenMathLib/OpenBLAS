/***************************************************************************
 * Copyright (c) 2025, The OpenBLAS Project
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
#include <arm_neon.h>

#include "common.h"

#ifdef BGEMM

#ifdef ALPHA_ONE

#define TO16 vcvth_bf16_f32
#define TO32 vcvtah_f32_bf16

#define UPDATE_C(PG, PTR, DST, SRC)                                            \
  do {                                                                         \
    DST = svreinterpret_f32_u32(svld1uh_u32((pghalf), (uint16_t*)PTR));        \
    DST = svadd_z((PG), SRC, DST);                                             \
    svtmp16 = svcvt_bf16_f32_z((PG), DST);                                     \
    svtmp16 = svuzp1_bf16(svtmp16, svtmp16);                                   \
    svst1_bf16((pghalf), (PTR), svtmp16);                                      \
  } while (0);
#define UPDATE_C2(ptr, tmp, vector)                             \
    *(ptr) = TO16(vector[0] + TO32(*ptr));                      \
    *(ptr + 1) = TO16(vector[1] + TO32(*(ptr + 1)));
#define UPDATE_C1(ptr, value) *ptr = TO16(TO32(*ptr) + (value))

#else

#define UPDATE_C(PG, PTR, DST, SRC)                                            \
  do {                                                                         \
    DST = svreinterpret_f32_u32(svld1uh_u32((pghalf), (uint16_t*)PTR));        \
    DST = svmad_z((PG), svalpha, SRC, DST);                                    \
    svtmp16 = svcvt_bf16_f32_z((PG), DST);                                     \
    svtmp16 = svuzp1_bf16(svtmp16, svtmp16);                                   \
    svst1_bf16((pghalf), (PTR), svtmp16);                                      \
  } while (0);
#define UPDATE_C2(ptr, tmp, vector)                                 \
    *(ptr) = TO16(vector[0] * alpha + TO32(*ptr));                  \
    *(ptr + 1) = TO16(vector[1] * alpha + TO32(*(ptr + 1)));
#define UPDATE_C1(ptr, value) *ptr = TO16(TO32(*ptr) + (value) * alpha)

#endif

#else

#ifdef ALPHA_ONE

#define UPDATE_C(PG, PTR, DST, SRC)                                            \
  do {                                                                         \
    DST = svld1_f32((PG), (PTR));                                              \
    DST = svadd_z((PG), SRC, DST);                                             \
    svst1_f32((PG), (PTR), DST);                                               \
  } while (0);
#define UPDATE_C2(ptr, tmp, vector)                     \
    tmp = vld1_f32(ptr);                                \
    tmp = vadd_f32(vector, tmp);                        \
    vst1_f32(ptr, tmp);
#define UPDATE_C1(ptr, value) *ptr = *ptr + (value)

#else

#define UPDATE_C(PG, PTR, DST, SRC)                                            \
  do {                                                                         \
    DST = svld1_f32((PG), (PTR));                                              \
    DST = svmad_z((PG), svalpha, SRC, DST);                                    \
    svst1_f32((PG), (PTR), DST);                                               \
  } while (0);
#define UPDATE_C2(ptr, tmp, vector)                 \
    tmp = vld1_f32(ptr);                            \
    tmp = vmla_n_f32(tmp, vector, alpha);           \
    vst1_f32(ptr, tmp);
#define UPDATE_C1(ptr, value) *ptr = *ptr + (value) * alpha

#endif

#endif

#ifdef BGEMM
#define OUTPUT_FLOAT bfloat16_t
#else
#define OUTPUT_FLOAT float
#endif

#ifdef ALPHA_ONE
static int bgemm_kernel_neoversev1_alpha_one(BLASLONG m, BLASLONG n, BLASLONG k,
                                       float alpha, IFLOAT *AA, IFLOAT *BB,
                                       FLOAT *CC, BLASLONG ldc)
#else
static int bgemm_kernel_neoversev1_alpha(BLASLONG m, BLASLONG n, BLASLONG k,
                                   float alpha, IFLOAT *AA, IFLOAT *BB, FLOAT *CC,
                                   BLASLONG ldc)
#endif
{
    BLASLONG pad_k = (k + 3) & ~3;

#ifndef ALPHA_ONE
    svfloat32_t svalpha = svdup_f32(alpha);
#endif

    bfloat16_t *ptr_a = (bfloat16_t *)AA;
    bfloat16_t *ptr_b = (bfloat16_t *)BB;
    OUTPUT_FLOAT *ptr_c =(OUTPUT_FLOAT*)CC;

    bfloat16_t *ptr_a0;
    bfloat16_t *ptr_b0;
    OUTPUT_FLOAT *ptr_c0, *ptr_c1, *ptr_c2, *ptr_c3;
    svfloat32_t tmp0, tmp1, tmp2, tmp3;
#ifdef BGEMM
    svbfloat16_t svtmp16;
#else
    float32x2_t tmp4, tmp5, tmp6, tmp7;
#endif

    const int sve_size_bf16 = svcnth();
    const int num_accumulators = sve_size_bf16 >> 1;

    svbool_t pgtrue = svptrue_b16();
#ifdef BGEMM
    // For BF16 load/store we use half the vector size
    svbool_t pghalf = svwhilelt_b16(0, num_accumulators);
#endif

    // N values are 4x2 packed matrices
    int n_step = 0;
    const int n2 = n & -2;
    const int n4 = n & -4;

    // For 256-bit this would be 8
    const int m_acc = (m & -num_accumulators);
    const int m2 = m & -2;

    for (; n_step < n4; n_step += 4) {
        ptr_a = (bfloat16_t *)AA;
        ptr_c0 = ptr_c;
        ptr_c1 = ptr_c0 + ldc;
        ptr_c2 = ptr_c1 + ldc;
        ptr_c3 = ptr_c2 + ldc;
        ptr_c += 4 * ldc;

        int m_step = 0;
        for (; m_step < m_acc; m_step += num_accumulators) {
            svfloat32_t acc0 = svdup_f32(0);
            svfloat32_t acc1 = svdup_f32(0);
            svfloat32_t acc2 = svdup_f32(0);
            svfloat32_t acc3 = svdup_f32(0);

            ptr_a0 = ptr_a;
            ptr_b0 = ptr_b;
            ptr_a += num_accumulators * pad_k;

            // Load entire 2VL block
            for (BLASLONG p = 0; p < pad_k; p += 4) {
                svbfloat16_t ma0 = svld1_bf16(pgtrue, ptr_a0);
                svbfloat16_t ma1 = svld1_bf16(pgtrue, ptr_a0 + sve_size_bf16);
                svbfloat16_t mb0 = svld1rq_bf16(pgtrue, ptr_b0);
                svbfloat16_t mb1 = svld1rq_bf16(pgtrue, ptr_b0 + 8);

                acc0 = svbfmmla_f32(acc0, mb0, ma0);
                acc1 = svbfmmla_f32(acc1, mb0, ma1);
                acc2 = svbfmmla_f32(acc2, mb1, ma0);
                acc3 = svbfmmla_f32(acc3, mb1, ma1);

                ptr_a0 += sve_size_bf16 * 2;
                ptr_b0 += 16;
            }

            svfloat32_t out0 = svreinterpret_f32_u64(svuzp1_u64(svreinterpret_u64_f32(acc0), svreinterpret_u64_f32(acc1)));
            svfloat32_t out1 = svreinterpret_f32_u64(svuzp2_u64(svreinterpret_u64_f32(acc0), svreinterpret_u64_f32(acc1)));

            svfloat32_t out2 = svreinterpret_f32_u64(svuzp1_u64(svreinterpret_u64_f32(acc2), svreinterpret_u64_f32(acc3)));
            svfloat32_t out3 = svreinterpret_f32_u64(svuzp2_u64(svreinterpret_u64_f32(acc2), svreinterpret_u64_f32(acc3)));
            
            UPDATE_C(pgtrue, ptr_c0, tmp0, out0);
            UPDATE_C(pgtrue, ptr_c1, tmp1, out1);
            UPDATE_C(pgtrue, ptr_c2, tmp2, out2);
            UPDATE_C(pgtrue, ptr_c3, tmp3, out3);

            ptr_c0 += num_accumulators;
            ptr_c1 += num_accumulators;
            ptr_c2 += num_accumulators;
            ptr_c3 += num_accumulators;
        }

        for (; m_step < m2; m_step += 2) {
            float32x4_t acc0 = {0,0,0,0};
            float32x4_t acc1 = {0,0,0,0};

            ptr_a0 = ptr_a;
            ptr_b0 = ptr_b;
            ptr_a += 2 * pad_k;

            for (BLASLONG p = 0; p < pad_k; p += 4) {
                bfloat16x8_t ma0 = vld1q_bf16(ptr_a0);
                bfloat16x8_t mb0 = vld1q_bf16(ptr_b0);
                bfloat16x8_t mb1 = vld1q_bf16(ptr_b0 + 8);

                acc0 = vbfmmlaq_f32(acc0, mb0, ma0);
                acc1 = vbfmmlaq_f32(acc1, mb1, ma0);

                ptr_a0 += 8;
                ptr_b0 += 16;
            }

            UPDATE_C2(ptr_c0, tmp4, vget_low_f32(acc0));
            UPDATE_C2(ptr_c1, tmp5, vget_high_f32(acc0));
            UPDATE_C2(ptr_c2, tmp6, vget_low_f32(acc1));
            UPDATE_C2(ptr_c3, tmp7, vget_high_f32(acc1));
            
            ptr_c0 += 2;
            ptr_c1 += 2;
            ptr_c2 += 2;
            ptr_c3 += 2;
        }
        
        // Final row is always a contiguous single row
        if (m & 1) {
            ptr_a0 = ptr_a;
            ptr_b0 = ptr_b;
            float32x4_t acc0 = {0,0,0,0};
            float32x4_t acc1 = {0,0,0,0};

            for (BLASLONG p = 0; p < pad_k; p += 4) {
                /// Same A value can be used for both B values
                bfloat16x8_t ma0 = vreinterpretq_bf16_u64(vdupq_n_u64(
                    *((uint64_t*)ptr_a0)
                ));
                bfloat16x8_t mb0 = vld1q_bf16(ptr_b0);
                bfloat16x8_t mb1 = vld1q_bf16(ptr_b0 + 8);

                acc0 = vbfmmlaq_f32(acc0, mb0, ma0);
                acc1 = vbfmmlaq_f32(acc1, mb1, ma0);

                ptr_a0 += 4;
                ptr_b0 += 16;
            }

            UPDATE_C1(ptr_c0, acc0[1]);
            UPDATE_C1(ptr_c1, acc0[3]);
            UPDATE_C1(ptr_c2, acc1[1]);
            UPDATE_C1(ptr_c3, acc1[3]);
        }

        ptr_b += 4 * pad_k;
    }

    for (; n_step < n2; n_step += 2) {
        ptr_a = (bfloat16_t *)AA;
        ptr_c0 = ptr_c;
        ptr_c1 = ptr_c0 + ldc;
        ptr_c += 2 * ldc;

        // Sets of two are contiguously packed so yay
        int m_step = 0;
        for (; m_step < m_acc; m_step += num_accumulators) {
            svfloat32_t acc0 = svdup_f32(0);
            svfloat32_t acc1 = svdup_f32(0);
            ptr_a0 = ptr_a;
            ptr_b0 = ptr_b;
            ptr_a += num_accumulators * pad_k;

            // Load entire 8x4 block
            for (BLASLONG p = 0; p < pad_k; p += 4) {
                svbfloat16_t ma0 = svld1_bf16(pgtrue, ptr_a0);
                svbfloat16_t ma1 = svld1_bf16(pgtrue, ptr_a0 + sve_size_bf16);
                svbfloat16_t mb0 = svld1rq_bf16(pgtrue, ptr_b0);

                acc0 = svbfmmla_f32(acc0, mb0, ma0);
                acc1 = svbfmmla_f32(acc1, mb0, ma1);

                ptr_a0 += sve_size_bf16 * 2;
                ptr_b0 += 8;
            }

            svfloat32_t out0 = svreinterpret_f32_u64(svuzp1_u64(svreinterpret_u64_f32(acc0), svreinterpret_u64_f32(acc1)));
            svfloat32_t out1 = svreinterpret_f32_u64(svuzp2_u64(svreinterpret_u64_f32(acc0), svreinterpret_u64_f32(acc1)));

            UPDATE_C(pgtrue, ptr_c0, tmp0, out0);
            UPDATE_C(pgtrue, ptr_c1, tmp1, out1);

            ptr_c0 += num_accumulators;
            ptr_c1 += num_accumulators;
        }

        for (; m_step < m2; m_step += 2) {
            float32x4_t acc = {0,0,0,0};

            ptr_a0 = ptr_a;
            ptr_b0 = ptr_b;
            ptr_a += 2 * pad_k;

            for (BLASLONG p = 0; p < pad_k; p += 4) {
                bfloat16x8_t ma0 = vld1q_bf16(ptr_a0);
                bfloat16x8_t mb0 = vld1q_bf16(ptr_b0);

                acc = vbfmmlaq_f32(acc, mb0, ma0);

                ptr_a0 += 8;
                ptr_b0 += 8;
            }

            UPDATE_C2(ptr_c0, tmp4, vget_low_f32(acc));
            UPDATE_C2(ptr_c1, tmp5, vget_high_f32(acc));
            
            ptr_c0 += 2;
            ptr_c1 += 2;
        }

        // Final row is always a contiguous single row
        if (m & 1) {
            ptr_a0 = ptr_a;
            ptr_b0 = ptr_b;
            float32x4_t acc = {0,0,0,0};

            for (BLASLONG p = 0; p < pad_k; p += 4) {
                /// Same A value can be used for both B values
                bfloat16x8_t ma0 = vreinterpretq_bf16_u64(vdupq_n_u64(
                    *((uint64_t*)ptr_a0)
                ));
                bfloat16x8_t mb0 = vld1q_bf16(ptr_b0);

                acc = vbfmmlaq_f32(acc, mb0, ma0);

                ptr_a0 += 4;
                ptr_b0 += 8;
            }

            UPDATE_C1(ptr_c0, acc[0]);
            UPDATE_C1(ptr_c1, acc[2]);
        }

        ptr_b += 2 * pad_k;
    }

    if (n & 1) {
        ptr_a = (bfloat16_t *)AA;
        ptr_c0 = ptr_c;

        int m_step = 0;
        for (; m_step < m_acc; m_step += num_accumulators) {
            ptr_a0 = ptr_a;
            ptr_b0 = ptr_b;
            ptr_a += num_accumulators * pad_k;

            svfloat32_t acc0 = svdup_f32(0);
            svfloat32_t acc1 = svdup_f32(0);

            // Load entire 8x4 block
            for (BLASLONG p = 0; p < pad_k; p += 4) {
                uint64_t* ptr_b0_u64 = (uint64_t*)ptr_b0;

                svbfloat16_t ma0 = svld1_bf16(pgtrue, ptr_a0);
                svbfloat16_t ma1 = svld1_bf16(pgtrue, ptr_a0 + sve_size_bf16);
                svbfloat16_t mb0 = svreinterpret_bf16_u64(svdup_u64(*ptr_b0_u64));

                acc0 = svbfmmla_f32(acc0, mb0, ma0);
                acc1 = svbfmmla_f32(acc1, mb0, ma1);

                ptr_a0 += sve_size_bf16 * 2;
                ptr_b0 += 4;
            }

            svfloat32_t out0 = svreinterpret_f32_u64(svuzp1_u64(svreinterpret_u64_f32(acc0), svreinterpret_u64_f32(acc1)));
            UPDATE_C(pgtrue, ptr_c0, tmp0, out0);

            ptr_c0 += num_accumulators;
        }

        for (; m_step < m2; m_step += 2) {
            float32x4_t acc = {0, 0, 0, 0};

            ptr_a0 = ptr_a;
            ptr_b0 = ptr_b;
            ptr_a += 2 * pad_k;

            for (BLASLONG p = 0; p < pad_k; p += 4) {
                bfloat16x8_t ma0 = vld1q_bf16(ptr_a0);
                bfloat16x8_t mb0 = vcombine_bf16(vld1_bf16(ptr_b0), vdup_n_bf16(vcvth_bf16_f32(0.0f)));

                acc = vbfmmlaq_f32(acc, mb0, ma0);

                ptr_a0 += 8;
                ptr_b0 += 4;
            }

            UPDATE_C2(ptr_c0, tmp4, vget_low_f32(acc));
            
            ptr_c0 += 2;
        }

        if (m & 1) {
            ptr_a0 = ptr_a;
            ptr_b0 = ptr_b;
            float32x2_t acc = {0,0};

            for (BLASLONG p = 0; p < pad_k; p += 4) {
                bfloat16x4_t ma0 = vld1_bf16(ptr_a0);
                bfloat16x4_t mb0 = vld1_bf16(ptr_b0);

                acc = vbfdot_f32(acc, ma0, mb0);

                ptr_a0 += 4;
                ptr_b0 += 4;
            }

            UPDATE_C1(ptr_c0, acc[0] + acc[1]);
        }
    }

    return 0;
}
