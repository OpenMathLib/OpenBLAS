/*
 Copyright (c) Qualcomm Technologies, Inc. and/or its subsidiaries.
 SPDX-License-Identifier: BSD-3-Clause-Clear
*/

#include "common.h"
#include <stdlib.h>
#include <inttypes.h>
#if defined(HAVE_SME)

#if defined(DYNAMIC_ARCH)
#define COMBINE(a,b) a ## b
#define COMBINE2(a,b) COMBINE(a,b)
#define SGEMM_PREPROCESS_BASE sgemm_direct_sme1_preprocess
#define SGEMM_PREPROCESS COMBINE2(SGEMM_PREPROCESS_BASE,TS)
#define SGEMM_DIRECT2X2_BASE sgemm_direct_alpha_beta_sme1_2VLx2VL
#define SGEMM_DIRECT2X2 COMBINE2(SGEMM_DIRECT2X2_BASE,TS)
#else
#define SGEMM_PREPROCESS sgemm_direct_sme1_preprocess
#define SGEMM_DIRECT2X2 sgemm_direct_alpha_beta_sme1_2VLx2VL
#endif

#if defined(__ARM_FEATURE_SME) && defined(__clang__) && __clang_major__ >= 16
#include <arm_sme.h>
#endif

/* Function prototypes */
extern void SGEMM_PREPROCESS (uint64_t nbr, uint64_t nbc,\

                                  const float * restrict a, float *  a_mod) ;

/* Function Definitions */
#if !defined(TRANSA)
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
// Outer product kernel.
// Computes a 2SVL x 2SVL block of C, utilizing all four FP32 tiles of ZA.
__attribute__((always_inline)) inline void
kernel_2x2(const float *A, const float *B, float *C, size_t shared_dim,
           size_t ldc, size_t block_rows, size_t block_cols, float alpha,
           float beta, uint64_t row_idx, uint64_t col_idx)
    __arm_out("za") __arm_streaming {

  const uint64_t svl = svcntw();
#if defined(TRANSA)
  size_t ldb = ldc;
#endif
  // Predicate set-up
  svbool_t pg = svptrue_b32();
  svbool_t pg_a_0 = svwhilelt_b32_u64(0, block_rows);
  svbool_t pg_a_1 = svwhilelt_b32_u64(svl, block_rows);

  svbool_t pg_b_0 = svwhilelt_b32_u64(0, block_cols);
  svbool_t pg_b_1 = svwhilelt_b32_u64(svl, block_cols);

#define pg_c_0 pg_b_0
#define pg_c_1 pg_b_1

  svzero_za();
  // beta == 0 must not read C; ZA is already initialized to zero.
  if (beta != 0.0f) {
    svfloat32_t beta_vec = svdup_f32(beta);

    // Load C to ZA
    for (size_t i = 0; i < MIN(svl, block_rows); i++) {
      svfloat32_t row_c_0 = svld1(pg_c_0, &C[i * ldc]);
      row_c_0 = svmul_x(pg, beta_vec, row_c_0);
      svwrite_hor_za32_f32_m(/*tile*/0, /*slice*/i, pg_c_0, row_c_0);

      svfloat32_t row_c_1 = svld1(pg_c_1, &C[i * ldc + svl]);
      row_c_1 = svmul_x(pg, beta_vec, row_c_1);
      svwrite_hor_za32_f32_m(/*tile*/1, /*slice*/i, pg_c_1, row_c_1);
    }
    for (size_t i = svl; i < block_rows; i++) {
      svfloat32_t row_c_0 = svld1(pg_c_0, &C[i * ldc]);
      row_c_0 = svmul_x(pg, beta_vec, row_c_0);
      svwrite_hor_za32_f32_m(/*tile*/2, /*slice*/i - svl, pg_c_0, row_c_0);

      svfloat32_t row_c_1 = svld1(pg_c_1, &C[i * ldc + svl]);
      row_c_1 = svmul_x(pg, beta_vec, row_c_1);
      svwrite_hor_za32_f32_m(/*tile*/3, /*slice*/i - svl, pg_c_1, row_c_1);
    }
  }

  svfloat32_t alpha_vec = svdup_f32(alpha);
  // Iterate through shared dimension (K)
  for (size_t k = 0; k < shared_dim; k++) {
#if !defined(TRANSA)
    // Load column of A
    svfloat32_t col_a_0 = svld1(pg_a_0, &A[k * svl]);
    col_a_0 = svmul_x(pg, alpha_vec, col_a_0);
    svfloat32_t col_a_1 = svld1(pg_a_1, &A[(k + shared_dim) * svl]);
    col_a_1 = svmul_x(pg, alpha_vec, col_a_1);

    // Load row of A**T
    svfloat32_t row_b_0 = svld1(pg_b_0, &B[k * svl]);
    svfloat32_t row_b_1 = svld1(pg_b_1, &B[(k + shared_dim) * svl]);
#else
    // Load column of A**T
    svfloat32_t col_a_0 = svld1(pg_a_0, &A[k * ldb]);
    col_a_0 = svmul_x(pg, alpha_vec, col_a_0);

    svfloat32_t col_a_1 = svld1(pg_a_1, &A[k * ldb + svl]);
    col_a_1 = svmul_x(pg, alpha_vec, col_a_1);

    // Load row of A
    svfloat32_t row_b_0 = svld1(pg_b_0, &B[k * ldb]);
    svfloat32_t row_b_1 = svld1(pg_b_1, &B[k * ldb + svl]);
#endif
    // Perform outer product
    svmopa_za32_m(/*tile*/0, pg, pg, col_a_0, row_b_0);
    svmopa_za32_m(/*tile*/1, pg, pg, col_a_0, row_b_1);
    svmopa_za32_m(/*tile*/2, pg, pg, col_a_1, row_b_0);
    svmopa_za32_m(/*tile*/3, pg, pg, col_a_1, row_b_1);
  }

#if defined(UPPER)
#define pg_c_0_full pg_c_0
#define pg_c_1_full pg_c_1

  bool need_update_pg_b = true;
  size_t last_invalid_index = col_idx - row_idx;
  // For Upper, If col_idx - row_idx >= 2*svl, we don't need to update the predicate due to all elements above the digonal
  if (col_idx - row_idx >= 2*svl) {
      need_update_pg_b = false;
  }
  // Store to C from ZA
  for (size_t i = 0; i < MIN(svl, block_rows); i++, last_invalid_index++) {
    if (need_update_pg_b) {
       pg_c_0 = svnot_b_z(pg_c_0_full, svwhilelt_b32_u64(0, last_invalid_index));
       pg_c_1 = svnot_b_z(pg_c_1_full, svwhilelt_b32_u64(svl, last_invalid_index));
    }

    svst1_hor_za32(/*tile*/0, /*slice*/i, pg_c_0, &C[i * ldc]);
    svst1_hor_za32(/*tile*/1, /*slice*/i, pg_c_1, &C[i * ldc + svl]);
  }
  for (size_t i = svl; i < block_rows; i++,last_invalid_index++) {
      if (need_update_pg_b) {
        pg_c_0 = svnot_b_z(pg_c_0_full, svwhilelt_b32_u64(0, last_invalid_index));
        pg_c_1 = svnot_b_z(pg_c_1_full, svwhilelt_b32_u64(svl, last_invalid_index));
      }
      svst1_hor_za32(/*tile*/2, /*slice*/i - svl, pg_c_0, &C[i * ldc]);
      svst1_hor_za32(/*tile*/3, /*slice*/i - svl, pg_c_1, &C[i * ldc + svl]);
  }
#else
  // Store to C from ZA
  size_t valid_index = row_idx - col_idx + 1;
  for (size_t i = 0; i < MIN(svl, block_rows); i++, valid_index++) {
    pg_c_0 = svwhilelt_b32_u64(0, MIN(valid_index, block_cols));
    pg_c_1 = svwhilelt_b32_u64(svl, MIN(valid_index, block_cols));
    svst1_hor_za32(/*tile*/0, /*slice*/i, pg_c_0, &C[i * ldc]);
    svst1_hor_za32(/*tile*/1, /*slice*/i, pg_c_1, &C[i * ldc + svl]);
  }
  for (size_t i = svl; i < block_rows; i++, valid_index++) {
    pg_c_0 = svwhilelt_b32_u64(0, MIN(valid_index, block_cols));
    pg_c_1 = svwhilelt_b32_u64(svl, MIN(valid_index, block_cols));
    svst1_hor_za32(/*tile*/2, /*slice*/i - svl, pg_c_0, &C[i * ldc]);
    svst1_hor_za32(/*tile*/3, /*slice*/i - svl, pg_c_1, &C[i * ldc + svl]);
  }
#endif
}

__arm_new("za") __arm_locally_streaming
static void ssyrk_direct_sme1_2VLx2VL(uint64_t n, uint64_t k, const float* alpha,\
                                   const float *ba, const float* beta, float *restrict bc) {
  const uint64_t num_rows = n;
  const uint64_t num_cols = n;

  const float *restrict a_ptr = ba;
  const float *restrict b_ptr = ba;
  float *restrict c_ptr = bc;

  const uint64_t svl = svcntw();
  const uint64_t ldc = n;

  // Block over rows of C (panels of A)
  uint64_t row_idx = 0;

  // 2x2 loop
  uint64_t row_batch = 2*svl;

  // Block over row dimension of C
  for (; row_idx < num_rows; row_idx += row_batch) {
    row_batch = MIN(row_batch, num_rows - row_idx);
    uint64_t col_batch = 2*svl;
#if defined(UPPER)
    // for UPLO is upper, Start from column col_idx = rows_index to ensure we only process the upper triangle (col_idx >= rows_index)
    for (uint64_t col_idx = row_idx; col_idx < num_cols; col_idx += col_batch) {
       col_batch = MIN(col_batch, num_cols - col_idx);
#else
    // for UPLO is lower, we only process the lower triangle part (col_idx <= row_idxx)
    for (uint64_t col_idx = 0; col_idx < num_cols && col_idx <= row_idx; col_idx += col_batch) {
#endif
      col_batch = MIN(col_batch, num_cols - col_idx);
#if !defined(TRANSA)
      kernel_2x2(&a_ptr[row_idx * k], &b_ptr[col_idx * k],
                 &c_ptr[row_idx * ldc + col_idx], k,
				 ldc, row_batch, col_batch, *alpha, *beta, row_idx, col_idx);
#else
      kernel_2x2(&a_ptr[row_idx], &b_ptr[col_idx],
                &c_ptr[row_idx * ldc + col_idx], k,
                ldc, row_batch, col_batch, *alpha, *beta, row_idx, col_idx);
#endif

    }
  }
  return;
}

#else
static void ssyrk_direct_sme1_2VLx2VL(uint64_t n, uint64_t k, const float* alpha,\
                                   const float *ba, const float* beta, float *restrict bc){}
#endif

void CNAME (BLASLONG N, BLASLONG K, float alpha, float * __restrict A,\
            BLASLONG strideA, float beta, float * __restrict C, BLASLONG strideC){
        if (alpha == 0.0f || K == 0) {
                if (beta == 1.0f)
                        return;
                ssyrk_direct_sme1_2VLx2VL(N, 0, &alpha, A, &beta, C);
                return;
        }

#if !defined(TRANSA)
        uint64_t n_mod, vl_elms;

        vl_elms = sve_cntw();

        n_mod = (((uint64_t)N + vl_elms - 1) / vl_elms) * vl_elms;

        float *A_mod = (float *) malloc(n_mod*K*sizeof(float));

	    /* Prevent compiler optimization by reading from memory instead
	     * of reading directly from vector (z) registers.
	     * */
        asm volatile("" : : :"p0", "p1", "p2", "p3", "p4", "p5", "p6", "p7",
                         "p8", "p9", "p10", "p11", "p12", "p13", "p14", "p15",
                         "z0", "z1", "z2", "z3", "z4", "z5", "z6", "z7",
                         "z8", "z9", "z10", "z11", "z12", "z13", "z14", "z15",
                         "z16", "z17", "z18", "z19", "z20", "z21", "z22", "z23",
                         "z24", "z25", "z26", "z27", "z28", "z29", "z30", "z31");

        /* Pre-process the left matrix to make it suitable for
           matrix sum of outer-product calculation
         */
        SGEMM_PREPROCESS (N, K, A, A_mod);
        asm volatile("" : : :"p0", "p1", "p2", "p3", "p4", "p5", "p6", "p7",
                         "p8", "p9", "p10", "p11", "p12", "p13", "p14", "p15",
                         "z0", "z1", "z2", "z3", "z4", "z5", "z6", "z7",
                         "z8", "z9", "z10", "z11", "z12", "z13", "z14", "z15",
                         "z16", "z17", "z18", "z19", "z20", "z21", "z22", "z23",
                         "z24", "z25", "z26", "z27", "z28", "z29", "z30", "z31");
        ssyrk_direct_sme1_2VLx2VL(N, K, &alpha, A_mod, &beta, C);
        free(A_mod);
#else
        ssyrk_direct_sme1_2VLx2VL(N, K, &alpha, A, &beta, C);
#endif

}

#else

void CNAME (BLASLONG N, BLASLONG K, float alpha, float * __restrict A,\
            BLASLONG strideA, float beta, float * __restrict C, BLASLONG strideC){
fprintf(stderr,"empty ssyrk_direct kernel should never be called\n");
}

#endif
