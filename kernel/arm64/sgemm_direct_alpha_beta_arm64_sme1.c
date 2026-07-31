/*
 Copyright (c) Qualcomm Technologies, Inc. and/or its subsidiaries.
 SPDX-License-Identifier: BSD-3-Clause-Clear
*/

#include "common.h"
#include <stdlib.h>
#include <inttypes.h>

#if defined(DYNAMIC_ARCH)
#define COMBINE(a,b) a ## b
#define COMBINE2(a,b) COMBINE(a,b)
#define SME1_PREPROCESS_BASE sgemm_direct_sme1_preprocess
#define SME1_PREPROCESS COMBINE2(SME1_PREPROCESS_BASE,TS)
#define SME1_KERNEL2X2_BASE sgemm_direct_alpha_beta_sme1_2VLx2VL
#define SME1_KERNEL2X2 COMBINE2(SME1_KERNEL2X2_BASE,TS)
#else
#define SME1_PREPROCESS sgemm_direct_sme1_preprocess
#define SME1_KERNEL2X2  sgemm_direct_alpha_beta_sme1_2VLx2VL
#endif

/* Function prototypes */
extern void SME1_PREPROCESS(uint64_t nbr, uint64_t nbc,\
                                  const float * restrict a, float *  a_mod);

#if defined(HAVE_SME)
#include "sme_abi.h"

#if defined(__ARM_FEATURE_SME) && defined(__clang__) && __clang_major__ >= 16
#include <arm_sme.h>
#endif

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

#if defined(__ARM_FEATURE_SME) && defined(__ARM_FEATURE_LOCALLY_STREAMING) && defined(__clang__) && __clang_major__ >= 16
// Outer product kernel.
// Computes a 2SVL x 2SVL block of C, utilizing all four FP32 tiles of ZA.
__attribute__((always_inline)) inline void
kernel_2x2(const float *A, const float *B, float *C, size_t shared_dim,
           size_t ldc, size_t block_rows, size_t block_cols, float alpha, float beta)
    __arm_out("za") __arm_streaming {

  const uint64_t svl = svcntw();
  size_t ldb = ldc;
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
      svwrite_hor_za32_f32_m(/*tile*/2, /*slice*/i, pg_c_0, row_c_0);

      svfloat32_t row_c_1 = svld1(pg_c_1, &C[i * ldc + svl]);
      row_c_1 = svmul_x(pg, beta_vec, row_c_1);
      svwrite_hor_za32_f32_m(/*tile*/3, /*slice*/i, pg_c_1, row_c_1);
    }
  }

  svfloat32_t alpha_vec = svdup_f32(alpha);
  // Iterate through shared dimension (K)
  for (size_t k = 0; k < shared_dim; k++) {
    // Load column of A
    svfloat32_t col_a_0 = svld1(pg_a_0, &A[k * svl]);
	col_a_0 = svmul_x(pg, alpha_vec, col_a_0);
    svfloat32_t col_a_1 = svld1(pg_a_1, &A[(k + shared_dim) * svl]);
    col_a_1 = svmul_x(pg, alpha_vec, col_a_1);
    // Load row of B
    svfloat32_t row_b_0 = svld1(pg_b_0, &B[k * ldb]);
    svfloat32_t row_b_1 = svld1(pg_b_1, &B[k * ldb + svl]);
    // Perform outer product
    svmopa_za32_m(/*tile*/0, pg, pg, col_a_0, row_b_0);
    svmopa_za32_m(/*tile*/1, pg, pg, col_a_0, row_b_1);
    svmopa_za32_m(/*tile*/2, pg, pg, col_a_1, row_b_0);
    svmopa_za32_m(/*tile*/3, pg, pg, col_a_1, row_b_1);
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
return;
}

__arm_new("za") __arm_locally_streaming
void SME1_KERNEL2X2(uint64_t m, uint64_t k, uint64_t n, const float* alpha,\
                                   const float *ba, const float *restrict bb, const float* beta,\
                                   float *restrict C) {

  const uint64_t num_rows = m;
  const uint64_t num_cols = n;

  const float *restrict a_ptr = ba;
  const float *restrict b_ptr = bb;
  float *restrict c_ptr = C;

  const uint64_t svl = svcntw();
  const uint64_t ldc = n;

  // Block over rows of C (panels of A)
  uint64_t row_idx = 0;

  // 2x2 loop
  uint64_t row_batch = 2*svl;

  // Block over row dimension of C
  for (; row_idx < num_rows; row_idx += row_batch) {
    row_batch = MIN(row_batch, num_rows - row_idx);

    uint64_t col_idx = 0;
    uint64_t col_batch = 2*svl;

    // Block over column dimension of C
    for (; col_idx < num_cols; col_idx += col_batch) {
      col_batch = MIN(col_batch, num_cols - col_idx);

      kernel_2x2(&a_ptr[row_idx * k], &b_ptr[col_idx],
                 &c_ptr[row_idx * ldc + col_idx], k,
				 ldc, row_batch, col_batch, *alpha, *beta);
    }
  }
  return;
}

#else
void SME1_KERNEL2X2(uint64_t m, uint64_t k, uint64_t n, const float* alpha,\
                                   const float *ba, const float *restrict bb, const float* beta,\
                                   float *restrict C){fprintf(stderr,"empty sgemm_alpha_beta2x2 should never get called!!!\n");}
#endif

/*void sgemm_kernel_direct (BLASLONG M, BLASLONG N, BLASLONG K,\
       float * __restrict A, BLASLONG strideA, float * __restrict B,\
       BLASLONG strideB , float * __restrict R, BLASLONG strideR)
*/
void CNAME (BLASLONG M, BLASLONG N, BLASLONG K, float alpha, float * __restrict A,\
            BLASLONG strideA, float * __restrict B, BLASLONG strideB ,\
            float beta, float * __restrict R, BLASLONG strideR){
        if (alpha == 0.0f || K == 0) {
                if (beta == 1.0f)
                        return;
                SME1_KERNEL2X2(M, 0, N, &alpha, A, B, &beta, R);
                return;
        }

        uint64_t m_mod, vl_elms;

        vl_elms = sve_cntw();

        m_mod = (((uint64_t)M + vl_elms - 1) / vl_elms) * vl_elms;

        float *A_mod = (float *) malloc(m_mod*K*sizeof(float));

	    /* Prevent compiler optimization by reading from memory instead
	     * of reading directly from vector (z) registers.
	     * */
        asm volatile("" : : :"p0", "p1", "p2", "p3", "p4", "p5", "p6", "p7",
                         "p8", "p9", "p10", "p11", "p12", "p13", "p14", "p15", "d8", "d9", "d10", "d11", "d12", "d13", "d14", "d15",
                         "z0", "z1", "z2", "z3", "z4", "z5", "z6", "z7",
                         "z8", "z9", "z10", "z11", "z12", "z13", "z14", "z15",
                         "z16", "z17", "z18", "z19", "z20", "z21", "z22", "z23",
                         "z24", "z25", "z26", "z27", "z28", "z29", "z30", "z31");

        /* Pre-process the left matrix to make it suitable for
           matrix sum of outer-product calculation
         */

        SME1_PREPROCESS(M, K, A, A_mod);

        asm volatile("" : : :"p0", "p1", "p2", "p3", "p4", "p5", "p6", "p7",
                         "p8", "p9", "p10", "p11", "p12", "p13", "p14", "p15","d8", "d9", "d10", "d11", "d12", "d13", "d14", "d15",
                         "z0", "z1", "z2", "z3", "z4", "z5", "z6", "z7",
                         "z8", "z9", "z10", "z11", "z12", "z13", "z14", "z15",
                         "z16", "z17", "z18", "z19", "z20", "z21", "z22", "z23",
                         "z24", "z25", "z26", "z27", "z28", "z29", "z30", "z31");

        /* Calculate C = alpha*A*B + beta*C */

        SME1_KERNEL2X2(M, K, N, &alpha, A_mod, B, &beta, R);

        free(A_mod);
}

#else

void CNAME (BLASLONG M, BLASLONG N, BLASLONG K, float alpha, float * __restrict A,\
            BLASLONG strideA, float * __restrict B, BLASLONG strideB ,\
            float beta, float * __restrict R, BLASLONG strideR){fprintf(stderr,"empty sgemm_direct_alpha_beta should not be called!!!\n");}
#endif


