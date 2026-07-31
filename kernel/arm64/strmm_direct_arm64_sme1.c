/*
 Copyright (c) Qualcomm Technologies, Inc. and/or its subsidiaries.
 SPDX-License-Identifier: BSD-3-Clause-Clear
*/

#include "common.h"
#include <stdlib.h>
#include <inttypes.h>
//#include "sme_abi.h"
#if defined(HAVE_SME)

#if defined(__ARM_FEATURE_SME) && defined(__clang__) && __clang_major__ >= 16
#include <arm_sme.h>
#endif

#if !defined(TRANSA)
/* Function Definitions */
static uint64_t sve_cntw() {
    uint64_t cnt;
    asm volatile(
        "rdsvl  %[res], #1\n"
        "lsr    %[res], %[res], #2\n"
        : [res] "=r" (cnt) ::
    );
    return cnt;
}
#endif

#if defined(__ARM_FEATURE_SME) && defined(__ARM_FEATURE_LOCALLY_STREAMING) && defined(__clang__) && __clang_major__ >= 16

#if !defined(TRANSA)
// Transpose 1SVL x N panel of A
__attribute__((always_inline))
inline static void transpose_panel_lower(const float *restrict a, float *restrict b,
                uint64_t rows, uint64_t cols,
                uint64_t a_step, uint64_t rows_index)
__arm_out("za") __arm_streaming {
    // for Lower Trangular Matrix
    uint64_t svl = svcntw();
    uint64_t col_batch = svl;

    svzero_za();
    uint64_t last_rows_index = rows_index + rows - 1;
    for (uint64_t k = 0; k < cols; k += col_batch) {
        if (last_rows_index < k) {
            // Early exit: if all rows are above the diagonal, no valid elements remain
            break;
        }
        // Load to horizontal slices
        for (uint64_t row = 0; row < rows; row++) {
            svbool_t pg_row = svwhilelt_b32_u64(k, MIN(rows_index + row + 1, cols));
            svld1_hor_za32(0, row, pg_row, &a[row * a_step + k]);
        }

        // Save from vertical slices
        col_batch = MIN(col_batch, cols - k);
        for (uint64_t col = 0; col < col_batch; col++) {
            svst1_ver_za32(0, col, svptrue_b32(), &b[(col + k) * svl]);
        }
    }
}

__attribute__((always_inline))
inline static void transpose_panel_upper(const float *restrict a, float *restrict b,
                uint64_t rows, uint64_t cols,
                uint64_t a_step, uint64_t rows_index)
__arm_out("za") __arm_streaming {
    // for Upper Trangular Matrix
    uint64_t svl = svcntw();
    uint64_t col_batch = svl;

    svzero_za();
    // Start from column k = rows_index to ensure we only process the upper triangle (k >= rows_index)
    for (uint64_t k = rows_index; k < cols; k += col_batch) {
        // Load to horizontal slices
        for (uint64_t row = 0; row < rows; row++) {
            svbool_t pg_row = svwhilelt_b32_u64(k, cols);
            svld1_hor_za32(0, row, pg_row, &a[row * a_step + k]);
        }

        // Save from vertical slices
        col_batch = MIN(col_batch, cols - k);
        for (uint64_t col = 0, real_col = k; col < col_batch; col++, real_col++) {
            // Only the upper triangular part of the matrix is stored.
            svbool_t pg_col = svwhilelt_b32_u64(rows_index, real_col + 1);
            svst1_ver_za32(0, col, pg_col, &b[(col + k) * svl]);
        }
    }
}

__arm_new("za") __arm_locally_streaming
static void strmm_direct_sme1_preprocess(uint64_t nbr, uint64_t nbc,
                const float *restrict a, float *restrict a_mod) {
    const uint64_t num_rows = nbr;
    uint64_t row_batch = svcntsw();
    for (uint64_t row_idx = 0; row_idx < num_rows; row_idx += row_batch) {
        // Transpose 1SVL x N panel of A
        row_batch = MIN(row_batch, num_rows - row_idx);
#if !defined(UPPER)
        transpose_panel_lower(&a[row_idx * nbc], &a_mod[row_idx * nbc], row_batch, nbc, nbc, row_idx);
#else
        transpose_panel_upper(&a[row_idx * nbc], &a_mod[row_idx * nbc], row_batch, nbc, nbc, row_idx);
#endif
    }
}
#endif

// Outer product kernel.
// Computes a 2SVL x 2SVL block of C, utilizing all four FP32 tiles of ZA.
__attribute__((always_inline)) inline void
kernel_2x2(const float *A, const float *B, float *C, size_t shared_dim,
           size_t ldc, size_t block_rows, size_t block_cols, float alpha, uint64_t row_idx)
    __arm_out("za") __arm_streaming {
    const uint64_t svl = svcntw();
    size_t ldb = ldc;
    // Predicate set-up
    svbool_t pg = svptrue_b32();
    svbool_t pg_a_0 = svwhilelt_b32_u64(0, block_rows);
    svbool_t pg_a_1 = svwhilelt_b32_u64(svl, block_rows);

#if (!defined(TRANSA) && !defined(UPPER)) || (defined(TRANSA) && defined(UPPER))
#define pg_a_0_full pg_a_0
#define pg_a_1_full pg_a_1
#endif
    svbool_t pg_b_0 = svwhilelt_b32_u64(0, block_cols);
    svbool_t pg_b_1 = svwhilelt_b32_u64(svl, block_cols);

#define pg_c_0 pg_b_0
#define pg_c_1 pg_b_1

    svzero_za();
    svfloat32_t alpha_vec = svdup_f32(alpha);
    // Iterate through shared dimension (K)
#if (!defined(TRANSA) && defined(UPPER)) || (defined(TRANSA) && !defined(UPPER))
    for (size_t k = row_idx, valid_index = 1; k < shared_dim; k++,valid_index++) {
        pg_a_0 = svwhilelt_b32_u64(0, MIN(valid_index, block_rows));
        pg_a_1 = svwhilelt_b32_u64(svl, MIN(valid_index, block_rows));
#else
    for (size_t k = 0; k < MIN(row_idx + block_rows, shared_dim); k++) {
        // If k exceeds row_idx, mask out rows before (k - row_idx)
        // This ensures only valid rows are included for lower triangular logic.
        if (k > row_idx) {
            pg_a_0 = svnot_b_z(pg_a_0_full, svwhilelt_b32_u64(0, k - row_idx));
            pg_a_1 = svnot_b_z(pg_a_1_full, svwhilelt_b32_u64(svl, k - row_idx));
        }
#endif

#if !defined(TRANSA)
        // Load column of A
        svfloat32_t col_a_0 = svld1(pg_a_0, &A[k * svl]);
        svfloat32_t col_a_1 = svld1(pg_a_1, &A[(k + shared_dim) * svl]);
#else
        svfloat32_t col_a_0 = svld1(pg_a_0, &A[k * shared_dim]);
        svfloat32_t col_a_1 = svld1(pg_a_1, &A[k * shared_dim + svl]);
#endif
        col_a_0 = svmul_x(pg_a_0, alpha_vec, col_a_0);
        col_a_1 = svmul_x(pg_a_1, alpha_vec, col_a_1);
        // Load row of B
        svfloat32_t row_b_0 = svld1(pg_b_0, &B[k * ldb]);
        svfloat32_t row_b_1 = svld1(pg_b_1, &B[k * ldb + svl]);
        // Perform outer product
        svmopa_za32_m(/*tile*/0, pg_a_0, pg, col_a_0, row_b_0);
        svmopa_za32_m(/*tile*/1, pg_a_0, pg, col_a_0, row_b_1);
        svmopa_za32_m(/*tile*/2, pg_a_1, pg, col_a_1, row_b_0);
        svmopa_za32_m(/*tile*/3, pg_a_1, pg, col_a_1, row_b_1);
    }

    // Store to C from ZA
    for (size_t i = 0; i < MIN(svl, block_rows); i++) {
      svst1_hor_za32(/*tile*/0, /*slice*/i, pg_c_0, &C[i * ldc]);
      svst1_hor_za32(/*tile*/1, /*slice*/i, pg_c_1, &C[i * ldc + svl]);
    }
    for (size_t i = svl; i < block_rows; i++) {
      svst1_hor_za32(/*tile*/2, /*slice*/i, pg_c_0, &C[i * ldc]);
      svst1_hor_za32(/*tile*/3, /*slice*/i, pg_c_1, &C[i * ldc + svl]);
    }

}

__arm_new("za") __arm_locally_streaming
static inline void strmm_direct_alpha_sme1_2VLx2VL(uint64_t m, uint64_t k, uint64_t n, const float* alpha,\
                                   const float *ba, float *restrict bb) {
    const uint64_t num_rows = m;
    const uint64_t num_cols = n;

    const float *restrict a_ptr = ba;
    const float *restrict b_ptr = bb;
    float *restrict c_ptr = bb;

    const uint64_t svl = svcntw();
    const uint64_t svl_x2 = 2*svl;
    const uint64_t ldc = n;


    uint64_t row_idx = 0;
#if (!defined(TRANSA) && defined(UPPER)) || (defined(TRANSA) && !defined(UPPER))
    // 2x2 loop
    uint64_t row_batch = svl_x2;
    // Block over rows of C (panels of A)
    for (; row_idx < num_rows; row_idx += row_batch) {
        row_batch = MIN(row_batch, num_rows - row_idx);
#else
    // Calculate the remainder of num_rows divided by 2VL to determine tail tile size
    uint64_t row_batch = num_rows % svl_x2;
    // If there's no remainder, use full tile size (2VL) for initial batch
    if (row_batch == 0) row_batch = svl_x2;
    // Loop from bottom to top, processing rows in batches
    for (uint64_t index = num_rows; index > 0; index -= row_batch, row_batch = svl_x2) {
        // Compute the starting row index for the current batch
        row_idx = index - row_batch;
#endif
        uint64_t col_idx = 0;
        uint64_t col_batch = svl_x2;
        // Block over column dimension of C
        for (; col_idx < num_cols; col_idx += col_batch) {
            col_batch = MIN(col_batch, num_cols - col_idx);
#if !defined(TRANSA)
            kernel_2x2(&a_ptr[row_idx * k], &b_ptr[col_idx],
                    &c_ptr[row_idx * ldc + col_idx], k,
                    ldc, row_batch, col_batch, *alpha, row_idx);
#else
            kernel_2x2(&a_ptr[row_idx], &b_ptr[col_idx],
                &c_ptr[row_idx * ldc + col_idx], k,
                ldc, row_batch, col_batch, *alpha, row_idx);
#endif
        }
    }

    return;
}

#else
#if !defined(TRANSA)
static void strmm_direct_sme1_preprocess(uint64_t nbr, uint64_t nbc,
                                         const float *restrict a, float *restrict a_mod) {}
#endif
static void strmm_direct_alpha_sme1_2VLx2VL(uint64_t m, uint64_t k, uint64_t n, const float* alpha,\
                                            const float *ba, float *restrict bb){}
#endif

void CNAME (BLASLONG M, BLASLONG N, float alpha, float * __restrict A,\
            BLASLONG strideA, float * __restrict B, BLASLONG strideB){
    if (alpha == 0.0f) {
        for (BLASLONG i = 0; i < M; i++) {
            float *row = B + i * strideB;
            for (BLASLONG j = 0; j < N; j++) {
                row[j] = 0.0f;
            }
        }
        return;
    }

#if !defined(TRANSA)
    uint64_t m_mod, vl_elms;

    vl_elms = sve_cntw();

    m_mod = (((uint64_t)M + vl_elms - 1) / vl_elms) * vl_elms;

    float *A_mod = (float *) malloc(m_mod*M*sizeof(float));
    strmm_direct_sme1_preprocess(M, M, A, A_mod);
    /* Calculate B = alpha*A*B*/
    strmm_direct_alpha_sme1_2VLx2VL(M, M, N, &alpha, A_mod, B);
    free(A_mod);
#else
    strmm_direct_alpha_sme1_2VLx2VL(M, M, N, &alpha, A, B);
#endif
}

#else
void CNAME (BLASLONG M, BLASLONG N, float alpha, float * __restrict A,\
            BLASLONG strideA, float * __restrict B, BLASLONG strideB){
            }

#endif
