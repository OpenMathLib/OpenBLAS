/*

AUTOGENERATED KERNEL
Script: ./kernel/riscv64/generate_kernel.py
Settings:
 LMUL=2
 M=8
 M_tail_scalar_from=2
 N=8
 __riscv_='__riscv_'
 complex=False
 conjugate=False
 cpu='zvl128b'
 force_acc_double=False
 index_type='BLASLONG'
 op='trmm'
 param_precision='float'
 reg_width_bits=128
 tail_policy=''
 trace=False

Derived:
 ELEN_ACC=32
 ELEN_PARAM=32
 LMUL_ACC=2
 VFMACC='__riscv_vfmacc_vf_f32m2'
 VFMUL='__riscv_vfmul_vf_f32m2'
 VLEV='__riscv_vle32_v_f32m2'
 VLSEV='__riscv_vlse32_v_f32m2'
 VMACC_TO_ACC='__riscv_vfmacc_vf_f32m2'
 VMUL_TO_ACC='__riscv_vfmul_vf_f32m2'
 VSETVL='__riscv_vsetvl_e32m2'
 VSEV='__riscv_vse32_v_f32m2'
 VSSEV='__riscv_vsse32_v_f32m2'
 acc_vector_t='vfloat32m2_t'
 output='strmm_kernel_8x8_zvl128b.c'
 param_scalar_t='float'
 param_vector_t='vfloat32m2_t'

*/

#include "common.h"

#if defined(LEFT) != defined(TRANSA)
#define BACKWARDS
#endif

int CNAME(BLASLONG M, BLASLONG N, BLASLONG K, FLOAT alpha, FLOAT *A, FLOAT *B, FLOAT *C, BLASLONG ldc, BLASLONG offset)

{
    BLASLONG gvl = 0;
    BLASLONG m_top = 0;
    BLASLONG n_top = 0;

    // -- MAIN PASS

    for (BLASLONG j = 0; j < N / 8; j += 1) {
        m_top = 0;
        BLASLONG gvl = __riscv_vsetvl_e32m2(8);

        for (BLASLONG i = 0; i < M / 8; i += 1) {
            BLASLONG ai = m_top * K;
            BLASLONG bi = n_top * K;
            BLASLONG pass_K = K;
#ifdef LEFT
            BLASLONG off = offset + m_top;
#else
            BLASLONG off = -offset + n_top;
#endif
#ifdef BACKWARDS
            ai += off * 8;
            bi += off * 8;
            pass_K -= off;
#else
#ifdef LEFT
            pass_K = off + 8;
#else
            pass_K = off + 8;
#endif
#endif
            float B0 = B[bi + 0];
            float B1 = B[bi + 1];
            float B2 = B[bi + 2];
            float B3 = B[bi + 3];
            float B4 = B[bi + 4];
            float B5 = B[bi + 5];
            float B6 = B[bi + 6];
            float B7 = B[bi + 7];
            bi += 8;

            vfloat32m2_t A0 = __riscv_vle32_v_f32m2(&A[ai + 0 * gvl], gvl);
            ai += 8;

            vfloat32m2_t result0 = __riscv_vfmul_vf_f32m2(A0, B0, gvl);
            vfloat32m2_t result1 = __riscv_vfmul_vf_f32m2(A0, B1, gvl);
            vfloat32m2_t result2 = __riscv_vfmul_vf_f32m2(A0, B2, gvl);
            vfloat32m2_t result3 = __riscv_vfmul_vf_f32m2(A0, B3, gvl);
            vfloat32m2_t result4 = __riscv_vfmul_vf_f32m2(A0, B4, gvl);
            vfloat32m2_t result5 = __riscv_vfmul_vf_f32m2(A0, B5, gvl);
            vfloat32m2_t result6 = __riscv_vfmul_vf_f32m2(A0, B6, gvl);
            vfloat32m2_t result7 = __riscv_vfmul_vf_f32m2(A0, B7, gvl);

            for (BLASLONG k = 1; k < pass_K; k++) {
                B0 = B[bi + 0];
                B1 = B[bi + 1];
                B2 = B[bi + 2];
                B3 = B[bi + 3];
                B4 = B[bi + 4];
                B5 = B[bi + 5];
                B6 = B[bi + 6];
                B7 = B[bi + 7];
                bi += 8;

                A0 = __riscv_vle32_v_f32m2(&A[ai + 0 * gvl], gvl);
                ai += 8;

                result0 = __riscv_vfmacc_vf_f32m2(result0, B0, A0, gvl);
                result1 = __riscv_vfmacc_vf_f32m2(result1, B1, A0, gvl);
                result2 = __riscv_vfmacc_vf_f32m2(result2, B2, A0, gvl);
                result3 = __riscv_vfmacc_vf_f32m2(result3, B3, A0, gvl);
                result4 = __riscv_vfmacc_vf_f32m2(result4, B4, A0, gvl);
                result5 = __riscv_vfmacc_vf_f32m2(result5, B5, A0, gvl);
                result6 = __riscv_vfmacc_vf_f32m2(result6, B6, A0, gvl);
                result7 = __riscv_vfmacc_vf_f32m2(result7, B7, A0, gvl);
            }

            BLASLONG ci = n_top * ldc + m_top;

            vfloat32m2_t c0 = __riscv_vfmul_vf_f32m2(result0, alpha, gvl);
            vfloat32m2_t c1 = __riscv_vfmul_vf_f32m2(result1, alpha, gvl);
            vfloat32m2_t c2 = __riscv_vfmul_vf_f32m2(result2, alpha, gvl);
            vfloat32m2_t c3 = __riscv_vfmul_vf_f32m2(result3, alpha, gvl);
            vfloat32m2_t c4 = __riscv_vfmul_vf_f32m2(result4, alpha, gvl);
            vfloat32m2_t c5 = __riscv_vfmul_vf_f32m2(result5, alpha, gvl);
            vfloat32m2_t c6 = __riscv_vfmul_vf_f32m2(result6, alpha, gvl);
            vfloat32m2_t c7 = __riscv_vfmul_vf_f32m2(result7, alpha, gvl);
            __riscv_vse32_v_f32m2(&C[ci], c0, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c1, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c2, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c3, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c4, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c5, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c6, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c7, gvl);
            m_top += 8;
        }

        // -- tails for main pass

        if (M & 4) {
            gvl = __riscv_vsetvl_e32m2(4);

            BLASLONG ai = m_top * K;
            BLASLONG bi = n_top * K;
            BLASLONG pass_K = K;
#ifdef LEFT
            BLASLONG off = offset + m_top;
#else
            BLASLONG off = -offset + n_top;
#endif
#ifdef BACKWARDS
            ai += off * 4;
            bi += off * 8;
            pass_K -= off;
#else
#ifdef LEFT
            pass_K = off + 4;
#else
            pass_K = off + 8;
#endif
#endif
            float B0 = B[bi + 0];
            float B1 = B[bi + 1];
            float B2 = B[bi + 2];
            float B3 = B[bi + 3];
            float B4 = B[bi + 4];
            float B5 = B[bi + 5];
            float B6 = B[bi + 6];
            float B7 = B[bi + 7];
            bi += 8;

            vfloat32m2_t A0 = __riscv_vle32_v_f32m2(&A[ai + 0 * gvl], gvl);
            ai += 4;

            vfloat32m2_t result0 = __riscv_vfmul_vf_f32m2(A0, B0, gvl);
            vfloat32m2_t result1 = __riscv_vfmul_vf_f32m2(A0, B1, gvl);
            vfloat32m2_t result2 = __riscv_vfmul_vf_f32m2(A0, B2, gvl);
            vfloat32m2_t result3 = __riscv_vfmul_vf_f32m2(A0, B3, gvl);
            vfloat32m2_t result4 = __riscv_vfmul_vf_f32m2(A0, B4, gvl);
            vfloat32m2_t result5 = __riscv_vfmul_vf_f32m2(A0, B5, gvl);
            vfloat32m2_t result6 = __riscv_vfmul_vf_f32m2(A0, B6, gvl);
            vfloat32m2_t result7 = __riscv_vfmul_vf_f32m2(A0, B7, gvl);

            for (BLASLONG k = 1; k < pass_K; k++) {
                B0 = B[bi + 0];
                B1 = B[bi + 1];
                B2 = B[bi + 2];
                B3 = B[bi + 3];
                B4 = B[bi + 4];
                B5 = B[bi + 5];
                B6 = B[bi + 6];
                B7 = B[bi + 7];
                bi += 8;

                A0 = __riscv_vle32_v_f32m2(&A[ai + 0 * gvl], gvl);
                ai += 4;

                result0 = __riscv_vfmacc_vf_f32m2(result0, B0, A0, gvl);
                result1 = __riscv_vfmacc_vf_f32m2(result1, B1, A0, gvl);
                result2 = __riscv_vfmacc_vf_f32m2(result2, B2, A0, gvl);
                result3 = __riscv_vfmacc_vf_f32m2(result3, B3, A0, gvl);
                result4 = __riscv_vfmacc_vf_f32m2(result4, B4, A0, gvl);
                result5 = __riscv_vfmacc_vf_f32m2(result5, B5, A0, gvl);
                result6 = __riscv_vfmacc_vf_f32m2(result6, B6, A0, gvl);
                result7 = __riscv_vfmacc_vf_f32m2(result7, B7, A0, gvl);
            }

            BLASLONG ci = n_top * ldc + m_top;

            vfloat32m2_t c0 = __riscv_vfmul_vf_f32m2(result0, alpha, gvl);
            vfloat32m2_t c1 = __riscv_vfmul_vf_f32m2(result1, alpha, gvl);
            vfloat32m2_t c2 = __riscv_vfmul_vf_f32m2(result2, alpha, gvl);
            vfloat32m2_t c3 = __riscv_vfmul_vf_f32m2(result3, alpha, gvl);
            vfloat32m2_t c4 = __riscv_vfmul_vf_f32m2(result4, alpha, gvl);
            vfloat32m2_t c5 = __riscv_vfmul_vf_f32m2(result5, alpha, gvl);
            vfloat32m2_t c6 = __riscv_vfmul_vf_f32m2(result6, alpha, gvl);
            vfloat32m2_t c7 = __riscv_vfmul_vf_f32m2(result7, alpha, gvl);
            __riscv_vse32_v_f32m2(&C[ci], c0, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c1, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c2, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c3, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c4, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c5, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c6, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c7, gvl);
            m_top += 4;
        }

        if (M & 2) {
            float result0 = 0;
            float result1 = 0;
            float result2 = 0;
            float result3 = 0;
            float result4 = 0;
            float result5 = 0;
            float result6 = 0;
            float result7 = 0;
            float result8 = 0;
            float result9 = 0;
            float result10 = 0;
            float result11 = 0;
            float result12 = 0;
            float result13 = 0;
            float result14 = 0;
            float result15 = 0;
            BLASLONG ai = m_top * K;
            BLASLONG bi = n_top * K;
            BLASLONG pass_K = K;
#ifdef LEFT
            BLASLONG off = offset + m_top;
#else
            BLASLONG off = -offset + n_top;
#endif
#ifdef BACKWARDS
            ai += off * 2;
            bi += off * 8;
            pass_K -= off;
#else
#ifdef LEFT
            pass_K = off + 2;
#else
            pass_K = off + 8;
#endif
#endif

            for (BLASLONG k = 0; k < pass_K; k++) {
                result0 += A[ai + 0] * B[bi + 0];
                result1 += A[ai + 1] * B[bi + 0];
                result2 += A[ai + 0] * B[bi + 1];
                result3 += A[ai + 1] * B[bi + 1];
                result4 += A[ai + 0] * B[bi + 2];
                result5 += A[ai + 1] * B[bi + 2];
                result6 += A[ai + 0] * B[bi + 3];
                result7 += A[ai + 1] * B[bi + 3];
                result8 += A[ai + 0] * B[bi + 4];
                result9 += A[ai + 1] * B[bi + 4];
                result10 += A[ai + 0] * B[bi + 5];
                result11 += A[ai + 1] * B[bi + 5];
                result12 += A[ai + 0] * B[bi + 6];
                result13 += A[ai + 1] * B[bi + 6];
                result14 += A[ai + 0] * B[bi + 7];
                result15 += A[ai + 1] * B[bi + 7];
                ai += 2;
                bi += 8;
            }

            BLASLONG ci = n_top * ldc + m_top;
            C[ci + 0 * ldc + 0] = alpha * result0;
            C[ci + 0 * ldc + 1] = alpha * result1;
            C[ci + 1 * ldc + 0] = alpha * result2;
            C[ci + 1 * ldc + 1] = alpha * result3;
            C[ci + 2 * ldc + 0] = alpha * result4;
            C[ci + 2 * ldc + 1] = alpha * result5;
            C[ci + 3 * ldc + 0] = alpha * result6;
            C[ci + 3 * ldc + 1] = alpha * result7;
            C[ci + 4 * ldc + 0] = alpha * result8;
            C[ci + 4 * ldc + 1] = alpha * result9;
            C[ci + 5 * ldc + 0] = alpha * result10;
            C[ci + 5 * ldc + 1] = alpha * result11;
            C[ci + 6 * ldc + 0] = alpha * result12;
            C[ci + 6 * ldc + 1] = alpha * result13;
            C[ci + 7 * ldc + 0] = alpha * result14;
            C[ci + 7 * ldc + 1] = alpha * result15;
            m_top += 2;
        }

        if (M & 1) {
            float result0 = 0;
            float result1 = 0;
            float result2 = 0;
            float result3 = 0;
            float result4 = 0;
            float result5 = 0;
            float result6 = 0;
            float result7 = 0;
            BLASLONG ai = m_top * K;
            BLASLONG bi = n_top * K;
            BLASLONG pass_K = K;
#ifdef LEFT
            BLASLONG off = offset + m_top;
#else
            BLASLONG off = -offset + n_top;
#endif
#ifdef BACKWARDS
            ai += off * 1;
            bi += off * 8;
            pass_K -= off;
#else
#ifdef LEFT
            pass_K = off + 1;
#else
            pass_K = off + 8;
#endif
#endif

            for (BLASLONG k = 0; k < pass_K; k++) {
                result0 += A[ai + 0] * B[bi + 0];
                result1 += A[ai + 0] * B[bi + 1];
                result2 += A[ai + 0] * B[bi + 2];
                result3 += A[ai + 0] * B[bi + 3];
                result4 += A[ai + 0] * B[bi + 4];
                result5 += A[ai + 0] * B[bi + 5];
                result6 += A[ai + 0] * B[bi + 6];
                result7 += A[ai + 0] * B[bi + 7];
                ai += 1;
                bi += 8;
            }

            BLASLONG ci = n_top * ldc + m_top;
            C[ci + 0 * ldc + 0] = alpha * result0;
            C[ci + 1 * ldc + 0] = alpha * result1;
            C[ci + 2 * ldc + 0] = alpha * result2;
            C[ci + 3 * ldc + 0] = alpha * result3;
            C[ci + 4 * ldc + 0] = alpha * result4;
            C[ci + 5 * ldc + 0] = alpha * result5;
            C[ci + 6 * ldc + 0] = alpha * result6;
            C[ci + 7 * ldc + 0] = alpha * result7;
            m_top += 1;
        }

        n_top += 8;
    }

    // -- tails for N=4

    if (N & 4) {
        gvl = __riscv_vsetvl_e32m2(8);
        m_top = 0;

        for (BLASLONG i = 0; i < M / 8; i += 1) {
            BLASLONG ai = m_top * K;
            BLASLONG bi = n_top * K;
            BLASLONG pass_K = K;
#ifdef LEFT
            BLASLONG off = offset + m_top;
#else
            BLASLONG off = -offset + n_top;
#endif
#ifdef BACKWARDS
            ai += off * 8;
            bi += off * 4;
            pass_K -= off;
#else
#ifdef LEFT
            pass_K = off + 8;
#else
            pass_K = off + 4;
#endif
#endif
            float B0 = B[bi + 0];
            float B1 = B[bi + 1];
            float B2 = B[bi + 2];
            float B3 = B[bi + 3];
            bi += 4;

            vfloat32m2_t A0 = __riscv_vle32_v_f32m2(&A[ai + 0 * gvl], gvl);
            ai += 8;

            vfloat32m2_t result0 = __riscv_vfmul_vf_f32m2(A0, B0, gvl);
            vfloat32m2_t result1 = __riscv_vfmul_vf_f32m2(A0, B1, gvl);
            vfloat32m2_t result2 = __riscv_vfmul_vf_f32m2(A0, B2, gvl);
            vfloat32m2_t result3 = __riscv_vfmul_vf_f32m2(A0, B3, gvl);

            for (BLASLONG k = 1; k < pass_K; k++) {
                B0 = B[bi + 0];
                B1 = B[bi + 1];
                B2 = B[bi + 2];
                B3 = B[bi + 3];
                bi += 4;

                A0 = __riscv_vle32_v_f32m2(&A[ai + 0 * gvl], gvl);
                ai += 8;

                result0 = __riscv_vfmacc_vf_f32m2(result0, B0, A0, gvl);
                result1 = __riscv_vfmacc_vf_f32m2(result1, B1, A0, gvl);
                result2 = __riscv_vfmacc_vf_f32m2(result2, B2, A0, gvl);
                result3 = __riscv_vfmacc_vf_f32m2(result3, B3, A0, gvl);
            }

            BLASLONG ci = n_top * ldc + m_top;

            vfloat32m2_t c0 = __riscv_vfmul_vf_f32m2(result0, alpha, gvl);
            vfloat32m2_t c1 = __riscv_vfmul_vf_f32m2(result1, alpha, gvl);
            vfloat32m2_t c2 = __riscv_vfmul_vf_f32m2(result2, alpha, gvl);
            vfloat32m2_t c3 = __riscv_vfmul_vf_f32m2(result3, alpha, gvl);
            __riscv_vse32_v_f32m2(&C[ci], c0, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c1, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c2, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c3, gvl);
            m_top += 8;
        }

        if (M & 4) {
            gvl = __riscv_vsetvl_e32m2(4);

            BLASLONG ai = m_top * K;
            BLASLONG bi = n_top * K;
            BLASLONG pass_K = K;
#ifdef LEFT
            BLASLONG off = offset + m_top;
#else
            BLASLONG off = -offset + n_top;
#endif
#ifdef BACKWARDS
            ai += off * 4;
            bi += off * 4;
            pass_K -= off;
#else
#ifdef LEFT
            pass_K = off + 4;
#else
            pass_K = off + 4;
#endif
#endif
            float B0 = B[bi + 0];
            float B1 = B[bi + 1];
            float B2 = B[bi + 2];
            float B3 = B[bi + 3];
            bi += 4;

            vfloat32m2_t A0 = __riscv_vle32_v_f32m2(&A[ai + 0 * gvl], gvl);
            ai += 4;

            vfloat32m2_t result0 = __riscv_vfmul_vf_f32m2(A0, B0, gvl);
            vfloat32m2_t result1 = __riscv_vfmul_vf_f32m2(A0, B1, gvl);
            vfloat32m2_t result2 = __riscv_vfmul_vf_f32m2(A0, B2, gvl);
            vfloat32m2_t result3 = __riscv_vfmul_vf_f32m2(A0, B3, gvl);

            for (BLASLONG k = 1; k < pass_K; k++) {
                B0 = B[bi + 0];
                B1 = B[bi + 1];
                B2 = B[bi + 2];
                B3 = B[bi + 3];
                bi += 4;

                A0 = __riscv_vle32_v_f32m2(&A[ai + 0 * gvl], gvl);
                ai += 4;

                result0 = __riscv_vfmacc_vf_f32m2(result0, B0, A0, gvl);
                result1 = __riscv_vfmacc_vf_f32m2(result1, B1, A0, gvl);
                result2 = __riscv_vfmacc_vf_f32m2(result2, B2, A0, gvl);
                result3 = __riscv_vfmacc_vf_f32m2(result3, B3, A0, gvl);
            }

            BLASLONG ci = n_top * ldc + m_top;

            vfloat32m2_t c0 = __riscv_vfmul_vf_f32m2(result0, alpha, gvl);
            vfloat32m2_t c1 = __riscv_vfmul_vf_f32m2(result1, alpha, gvl);
            vfloat32m2_t c2 = __riscv_vfmul_vf_f32m2(result2, alpha, gvl);
            vfloat32m2_t c3 = __riscv_vfmul_vf_f32m2(result3, alpha, gvl);
            __riscv_vse32_v_f32m2(&C[ci], c0, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c1, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c2, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c3, gvl);
            m_top += 4;
        }

        if (M & 2) {
            float result0 = 0;
            float result1 = 0;
            float result2 = 0;
            float result3 = 0;
            float result4 = 0;
            float result5 = 0;
            float result6 = 0;
            float result7 = 0;
            BLASLONG ai = m_top * K;
            BLASLONG bi = n_top * K;
            BLASLONG pass_K = K;
#ifdef LEFT
            BLASLONG off = offset + m_top;
#else
            BLASLONG off = -offset + n_top;
#endif
#ifdef BACKWARDS
            ai += off * 2;
            bi += off * 4;
            pass_K -= off;
#else
#ifdef LEFT
            pass_K = off + 2;
#else
            pass_K = off + 4;
#endif
#endif

            for (BLASLONG k = 0; k < pass_K; k++) {
                result0 += A[ai + 0] * B[bi + 0];
                result1 += A[ai + 1] * B[bi + 0];
                result2 += A[ai + 0] * B[bi + 1];
                result3 += A[ai + 1] * B[bi + 1];
                result4 += A[ai + 0] * B[bi + 2];
                result5 += A[ai + 1] * B[bi + 2];
                result6 += A[ai + 0] * B[bi + 3];
                result7 += A[ai + 1] * B[bi + 3];
                ai += 2;
                bi += 4;
            }

            BLASLONG ci = n_top * ldc + m_top;
            C[ci + 0 * ldc + 0] = alpha * result0;
            C[ci + 0 * ldc + 1] = alpha * result1;
            C[ci + 1 * ldc + 0] = alpha * result2;
            C[ci + 1 * ldc + 1] = alpha * result3;
            C[ci + 2 * ldc + 0] = alpha * result4;
            C[ci + 2 * ldc + 1] = alpha * result5;
            C[ci + 3 * ldc + 0] = alpha * result6;
            C[ci + 3 * ldc + 1] = alpha * result7;
            m_top += 2;
        }

        if (M & 1) {
            float result0 = 0;
            float result1 = 0;
            float result2 = 0;
            float result3 = 0;
            BLASLONG ai = m_top * K;
            BLASLONG bi = n_top * K;
            BLASLONG pass_K = K;
#ifdef LEFT
            BLASLONG off = offset + m_top;
#else
            BLASLONG off = -offset + n_top;
#endif
#ifdef BACKWARDS
            ai += off * 1;
            bi += off * 4;
            pass_K -= off;
#else
#ifdef LEFT
            pass_K = off + 1;
#else
            pass_K = off + 4;
#endif
#endif

            for (BLASLONG k = 0; k < pass_K; k++) {
                result0 += A[ai + 0] * B[bi + 0];
                result1 += A[ai + 0] * B[bi + 1];
                result2 += A[ai + 0] * B[bi + 2];
                result3 += A[ai + 0] * B[bi + 3];
                ai += 1;
                bi += 4;
            }

            BLASLONG ci = n_top * ldc + m_top;
            C[ci + 0 * ldc + 0] = alpha * result0;
            C[ci + 1 * ldc + 0] = alpha * result1;
            C[ci + 2 * ldc + 0] = alpha * result2;
            C[ci + 3 * ldc + 0] = alpha * result3;
            m_top += 1;
        }

        n_top += 4;
    }

    // -- tails for N=2

    if (N & 2) {
        gvl = __riscv_vsetvl_e32m2(8);
        m_top = 0;

        for (BLASLONG i = 0; i < M / 8; i += 1) {
            BLASLONG ai = m_top * K;
            BLASLONG bi = n_top * K;
            BLASLONG pass_K = K;
#ifdef LEFT
            BLASLONG off = offset + m_top;
#else
            BLASLONG off = -offset + n_top;
#endif
#ifdef BACKWARDS
            ai += off * 8;
            bi += off * 2;
            pass_K -= off;
#else
#ifdef LEFT
            pass_K = off + 8;
#else
            pass_K = off + 2;
#endif
#endif
            float B0 = B[bi + 0];
            float B1 = B[bi + 1];
            bi += 2;

            vfloat32m2_t A0 = __riscv_vle32_v_f32m2(&A[ai + 0 * gvl], gvl);
            ai += 8;

            vfloat32m2_t result0 = __riscv_vfmul_vf_f32m2(A0, B0, gvl);
            vfloat32m2_t result1 = __riscv_vfmul_vf_f32m2(A0, B1, gvl);

            for (BLASLONG k = 1; k < pass_K; k++) {
                B0 = B[bi + 0];
                B1 = B[bi + 1];
                bi += 2;

                A0 = __riscv_vle32_v_f32m2(&A[ai + 0 * gvl], gvl);
                ai += 8;

                result0 = __riscv_vfmacc_vf_f32m2(result0, B0, A0, gvl);
                result1 = __riscv_vfmacc_vf_f32m2(result1, B1, A0, gvl);
            }

            BLASLONG ci = n_top * ldc + m_top;

            vfloat32m2_t c0 = __riscv_vfmul_vf_f32m2(result0, alpha, gvl);
            vfloat32m2_t c1 = __riscv_vfmul_vf_f32m2(result1, alpha, gvl);
            __riscv_vse32_v_f32m2(&C[ci], c0, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c1, gvl);
            m_top += 8;
        }

        if (M & 4) {
            gvl = __riscv_vsetvl_e32m2(4);

            BLASLONG ai = m_top * K;
            BLASLONG bi = n_top * K;
            BLASLONG pass_K = K;
#ifdef LEFT
            BLASLONG off = offset + m_top;
#else
            BLASLONG off = -offset + n_top;
#endif
#ifdef BACKWARDS
            ai += off * 4;
            bi += off * 2;
            pass_K -= off;
#else
#ifdef LEFT
            pass_K = off + 4;
#else
            pass_K = off + 2;
#endif
#endif
            float B0 = B[bi + 0];
            float B1 = B[bi + 1];
            bi += 2;

            vfloat32m2_t A0 = __riscv_vle32_v_f32m2(&A[ai + 0 * gvl], gvl);
            ai += 4;

            vfloat32m2_t result0 = __riscv_vfmul_vf_f32m2(A0, B0, gvl);
            vfloat32m2_t result1 = __riscv_vfmul_vf_f32m2(A0, B1, gvl);

            for (BLASLONG k = 1; k < pass_K; k++) {
                B0 = B[bi + 0];
                B1 = B[bi + 1];
                bi += 2;

                A0 = __riscv_vle32_v_f32m2(&A[ai + 0 * gvl], gvl);
                ai += 4;

                result0 = __riscv_vfmacc_vf_f32m2(result0, B0, A0, gvl);
                result1 = __riscv_vfmacc_vf_f32m2(result1, B1, A0, gvl);
            }

            BLASLONG ci = n_top * ldc + m_top;

            vfloat32m2_t c0 = __riscv_vfmul_vf_f32m2(result0, alpha, gvl);
            vfloat32m2_t c1 = __riscv_vfmul_vf_f32m2(result1, alpha, gvl);
            __riscv_vse32_v_f32m2(&C[ci], c0, gvl);
            ci += ldc - gvl * 0;
            __riscv_vse32_v_f32m2(&C[ci], c1, gvl);
            m_top += 4;
        }

        if (M & 2) {
            float result0 = 0;
            float result1 = 0;
            float result2 = 0;
            float result3 = 0;
            BLASLONG ai = m_top * K;
            BLASLONG bi = n_top * K;
            BLASLONG pass_K = K;
#ifdef LEFT
            BLASLONG off = offset + m_top;
#else
            BLASLONG off = -offset + n_top;
#endif
#ifdef BACKWARDS
            ai += off * 2;
            bi += off * 2;
            pass_K -= off;
#else
#ifdef LEFT
            pass_K = off + 2;
#else
            pass_K = off + 2;
#endif
#endif

            for (BLASLONG k = 0; k < pass_K; k++) {
                result0 += A[ai + 0] * B[bi + 0];
                result1 += A[ai + 1] * B[bi + 0];
                result2 += A[ai + 0] * B[bi + 1];
                result3 += A[ai + 1] * B[bi + 1];
                ai += 2;
                bi += 2;
            }

            BLASLONG ci = n_top * ldc + m_top;
            C[ci + 0 * ldc + 0] = alpha * result0;
            C[ci + 0 * ldc + 1] = alpha * result1;
            C[ci + 1 * ldc + 0] = alpha * result2;
            C[ci + 1 * ldc + 1] = alpha * result3;
            m_top += 2;
        }

        if (M & 1) {
            float result0 = 0;
            float result1 = 0;
            BLASLONG ai = m_top * K;
            BLASLONG bi = n_top * K;
            BLASLONG pass_K = K;
#ifdef LEFT
            BLASLONG off = offset + m_top;
#else
            BLASLONG off = -offset + n_top;
#endif
#ifdef BACKWARDS
            ai += off * 1;
            bi += off * 2;
            pass_K -= off;
#else
#ifdef LEFT
            pass_K = off + 1;
#else
            pass_K = off + 2;
#endif
#endif

            for (BLASLONG k = 0; k < pass_K; k++) {
                result0 += A[ai + 0] * B[bi + 0];
                result1 += A[ai + 0] * B[bi + 1];
                ai += 1;
                bi += 2;
            }

            BLASLONG ci = n_top * ldc + m_top;
            C[ci + 0 * ldc + 0] = alpha * result0;
            C[ci + 1 * ldc + 0] = alpha * result1;
            m_top += 1;
        }

        n_top += 2;
    }

    // -- tails for N=1

    if (N & 1) {
        gvl = __riscv_vsetvl_e32m2(8);
        m_top = 0;

        for (BLASLONG i = 0; i < M / 8; i += 1) {
            BLASLONG ai = m_top * K;
            BLASLONG bi = n_top * K;
            BLASLONG pass_K = K;
#ifdef LEFT
            BLASLONG off = offset + m_top;
#else
            BLASLONG off = -offset + n_top;
#endif
#ifdef BACKWARDS
            ai += off * 8;
            bi += off * 1;
            pass_K -= off;
#else
#ifdef LEFT
            pass_K = off + 8;
#else
            pass_K = off + 1;
#endif
#endif
            float B0 = B[bi + 0];
            bi += 1;

            vfloat32m2_t A0 = __riscv_vle32_v_f32m2(&A[ai + 0 * gvl], gvl);
            ai += 8;

            vfloat32m2_t result0 = __riscv_vfmul_vf_f32m2(A0, B0, gvl);

            for (BLASLONG k = 1; k < pass_K; k++) {
                B0 = B[bi + 0];
                bi += 1;

                A0 = __riscv_vle32_v_f32m2(&A[ai + 0 * gvl], gvl);
                ai += 8;

                result0 = __riscv_vfmacc_vf_f32m2(result0, B0, A0, gvl);
            }

            BLASLONG ci = n_top * ldc + m_top;

            vfloat32m2_t c0 = __riscv_vfmul_vf_f32m2(result0, alpha, gvl);
            __riscv_vse32_v_f32m2(&C[ci], c0, gvl);
            m_top += 8;
        }

        if (M & 4) {
            gvl = __riscv_vsetvl_e32m2(4);

            BLASLONG ai = m_top * K;
            BLASLONG bi = n_top * K;
            BLASLONG pass_K = K;
#ifdef LEFT
            BLASLONG off = offset + m_top;
#else
            BLASLONG off = -offset + n_top;
#endif
#ifdef BACKWARDS
            ai += off * 4;
            bi += off * 1;
            pass_K -= off;
#else
#ifdef LEFT
            pass_K = off + 4;
#else
            pass_K = off + 1;
#endif
#endif
            float B0 = B[bi + 0];
            bi += 1;

            vfloat32m2_t A0 = __riscv_vle32_v_f32m2(&A[ai + 0 * gvl], gvl);
            ai += 4;

            vfloat32m2_t result0 = __riscv_vfmul_vf_f32m2(A0, B0, gvl);

            for (BLASLONG k = 1; k < pass_K; k++) {
                B0 = B[bi + 0];
                bi += 1;

                A0 = __riscv_vle32_v_f32m2(&A[ai + 0 * gvl], gvl);
                ai += 4;

                result0 = __riscv_vfmacc_vf_f32m2(result0, B0, A0, gvl);
            }

            BLASLONG ci = n_top * ldc + m_top;

            vfloat32m2_t c0 = __riscv_vfmul_vf_f32m2(result0, alpha, gvl);
            __riscv_vse32_v_f32m2(&C[ci], c0, gvl);
            m_top += 4;
        }

        if (M & 2) {
            float result0 = 0;
            float result1 = 0;
            BLASLONG ai = m_top * K;
            BLASLONG bi = n_top * K;
            BLASLONG pass_K = K;
#ifdef LEFT
            BLASLONG off = offset + m_top;
#else
            BLASLONG off = -offset + n_top;
#endif
#ifdef BACKWARDS
            ai += off * 2;
            bi += off * 1;
            pass_K -= off;
#else
#ifdef LEFT
            pass_K = off + 2;
#else
            pass_K = off + 1;
#endif
#endif

            for (BLASLONG k = 0; k < pass_K; k++) {
                result0 += A[ai + 0] * B[bi + 0];
                result1 += A[ai + 1] * B[bi + 0];
                ai += 2;
                bi += 1;
            }

            BLASLONG ci = n_top * ldc + m_top;
            C[ci + 0 * ldc + 0] = alpha * result0;
            C[ci + 0 * ldc + 1] = alpha * result1;
            m_top += 2;
        }

        if (M & 1) {
            float result0 = 0;
            BLASLONG ai = m_top * K;
            BLASLONG bi = n_top * K;
            BLASLONG pass_K = K;
#ifdef LEFT
            BLASLONG off = offset + m_top;
#else
            BLASLONG off = -offset + n_top;
#endif
#ifdef BACKWARDS
            ai += off * 1;
            bi += off * 1;
            pass_K -= off;
#else
#ifdef LEFT
            pass_K = off + 1;
#else
            pass_K = off + 1;
#endif
#endif

            for (BLASLONG k = 0; k < pass_K; k++) {
                result0 += A[ai + 0] * B[bi + 0];
                ai += 1;
                bi += 1;
            }

            BLASLONG ci = n_top * ldc + m_top;
            C[ci + 0 * ldc + 0] = alpha * result0;
            m_top += 1;
        }

        n_top += 1;
    }

    return 0;
}
