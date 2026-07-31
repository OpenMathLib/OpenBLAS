/*
 Copyright (c) Qualcomm Technologies, Inc. and/or its subsidiaries.
 SPDX-License-Identifier: BSD-3-Clause-Clear
*/

#include "common.h"
#include <stdlib.h>
#include <inttypes.h>
// #include "sme_abi.h"
#if defined(HAVE_SME)

#if defined(__ARM_FEATURE_SME) && defined(__clang__) && __clang_major__ >= 16
#include <arm_sme.h>
#endif

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

/* Function prototypes */
extern void SGEMM_PREPROCESS(uint64_t nbr, uint64_t nbc,\
                                  const float * restrict a, float *  a_mod);

extern void SGEMM_DIRECT2X2(uint64_t m, uint64_t k, uint64_t n, const float* alpha,\
                                   const float *ba, const float *restrict bb, const float* beta,\
                                   float *restrict C);
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

#if defined(UPPER)
__arm_new("za") __arm_locally_streaming
static void ssymm_direct_sme1_preprocessLU(uint64_t nbr, uint64_t nbc,
                const float *restrict a, float *restrict a_mod)
{
  // const uint64_t num_rows = nbr;
  // const uint64_t num_cols = nbc;
  const uint64_t svl = svcntw();
  uint64_t row_batch = svl;

  const float *restrict pSrc;
  float *restrict pDst;
  for (uint64_t row_idx = 0; row_idx < nbr; row_idx += row_batch)
  {
    row_batch = MIN(row_batch, nbr - row_idx);

    // Fill in the lower triangle and Transpose 1SVL x N panel of A
    uint64_t col_batch = svl;

    for (uint64_t col_idx = 0; col_idx < nbc; col_idx += col_batch)
    {
      svzero_za();

      if (col_idx == row_idx)
      {
        pSrc = &a[(row_idx)*nbc + col_idx];
        pDst = &a_mod[(col_idx)*svl + row_idx * nbc];
        // Load horizontal slices, filling lower elements
        const svbool_t pg_row = svwhilelt_b32_u64(col_idx, nbc);
        for (int64_t row = row_batch - 1; row >= 0; row--)
        {
          svld1_hor_za32(0, row, pg_row, &pSrc[row * nbc]);
          svld1_ver_za32(0, row, pg_row, &pSrc[row * nbc]);
        }
        // Save vertical slices
        col_batch = MIN(col_batch, nbc - col_idx);
        for (uint64_t col = 0; col < col_batch; col++)
        {
          svst1_ver_za32(0, col, svptrue_b32(), &pDst[col * svl]);
        }
      }
      else if (col_idx > row_idx)
      {
        pSrc = &a[(row_idx)*nbc + col_idx];
        pDst = &a_mod[(col_idx)*svl + row_idx * nbc];
        // Load horizontal slices
        const svbool_t pg_row = svwhilelt_b32_u64(col_idx, nbc);
        for (uint64_t row = 0; row < row_batch; row++)
        {
          svld1_hor_za32(0, row, pg_row, &pSrc[row * nbc]);
        }
        // Save vertical slices
        col_batch = MIN(col_batch, nbc - col_idx);
        for (uint64_t col = 0; col < col_batch; col++)
        {
          svst1_ver_za32(0, col, svptrue_b32(), &pDst[col * svl]);
        }
      }
      else if (col_idx < row_idx)
      {
        pSrc = &a[row_idx + col_idx * nbc];
        pDst = &a_mod[(col_idx)*svl + row_idx * nbc];
        // Load horizontal slices
        const svbool_t pg_row = svwhilelt_b32_u64(row_idx, nbc);
        for (uint64_t row = 0; row < svl; row++)
        {
          svld1_hor_za32(0, row, pg_row, &pSrc[row * nbc]);
        }
        // Save vertical slices
        col_batch = MIN(col_batch, nbc - col_idx);
        for (uint64_t col = 0; col < svl; col++)
        {
          svst1_hor_za32(0, col, svptrue_b32(), &pDst[col * svl]);
        }
      }
    }
  }
}
#endif

//
#if defined(LOWER)
__arm_new("za") __arm_locally_streaming
static void ssymm_direct_sme1_preprocessLL(uint64_t nbr, uint64_t nbc,
                const float *restrict a, float *restrict a_mod)
{
  // const uint64_t num_rows = nbr;
  const uint64_t svl = svcntw();
  uint64_t row_batch = svl;

  const float *restrict pSrc;
  float *restrict pDst;
  for (uint64_t row_idx = 0; row_idx < nbr; row_idx += row_batch)
  {
    row_batch = MIN(row_batch, nbr - row_idx);

    // Fill in the upper triangle and Transpose 1SVL x N panel of A
    uint64_t col_batch = svl;

    for (uint64_t col_idx = 0; col_idx < nbc; col_idx += col_batch)
    {
      svzero_za();

      if (col_idx == row_idx)
      {
        pSrc = &a[(row_idx)*nbc + col_idx];
        pDst = &a_mod[(col_idx)*svl + row_idx * nbc];
        // Load horizontal slices, filling upper elements
        const svbool_t pg_row = svwhilelt_b32_u64(col_idx, nbc);
        for (uint64_t row = 0; row < row_batch; row++)
        {
          svld1_hor_za32(0, row, pg_row, &pSrc[row * nbc]);
          svld1_ver_za32(0, row, pg_row, &pSrc[row * nbc]);
        }
        // Save vertical slices
        col_batch = MIN(col_batch, nbc - col_idx);
        for (uint64_t col = 0; col < col_batch; col++)
        {
          svst1_ver_za32(0, col, svptrue_b32(), &pDst[col * svl]);
        }
      }
      else if (col_idx > row_idx)
      {
        pSrc = &a[row_idx + col_idx * nbc];
        pDst = &a_mod[(col_idx)*svl + row_idx * nbc];
        // Load horizontal slices
        const svbool_t pg_row = svptrue_b32();
        for (uint64_t row = 0; row < row_batch; row++)
        {
          svld1_hor_za32(0, row, pg_row, &pSrc[row * nbc]);
        }
        // Save vertical slices
        col_batch = MIN(col_batch, nbc - col_idx);
        for (uint64_t col = 0; col < col_batch; col++)
        {
          svst1_hor_za32(0, col, svptrue_b32(), &pDst[col * svl]);
        }
      }
      else if (col_idx < row_idx)
      {
        pSrc = &a[(row_idx)*nbc + col_idx];
        pDst = &a_mod[(col_idx)*svl + row_idx * nbc];
        // Load horizontal slices
        const svbool_t pg_row = svwhilelt_b32_u64(col_idx, nbc);
        for (uint64_t row = 0; row < row_batch; row++)
        {
          svld1_hor_za32(0, row, pg_row, &pSrc[row * nbc]);
        }
        // Save vertical slices
        col_batch = MIN(col_batch, nbc - col_idx);
        for (uint64_t col = 0; col < col_batch; col++)
        {
          svst1_ver_za32(0, col, svptrue_b32(), &pDst[col * svl]);
        }
      }
    }
  }
}
#endif
#else
#if defined(UPPER)
static void ssymm_direct_sme1_preprocessLU(uint64_t nbr, uint64_t nbc,
                const float *restrict a, float *restrict a_mod){}
#endif
#if defined(LOWER)
static void ssymm_direct_sme1_preprocessLL(uint64_t nbr, uint64_t nbc,
                const float *restrict a, float *restrict a_mod){}
#endif
#endif

//
void CNAME(BLASLONG M, BLASLONG N, float alpha, float *__restrict A,
           BLASLONG strideA, float *__restrict B, BLASLONG strideB,
           float beta, float *__restrict R, BLASLONG strideR)
{
  if (alpha == 0.0f) {
    if (beta == 1.0f)
      return;
    SGEMM_DIRECT2X2(M, 0, N, &alpha, A, B, &beta, R);
    return;
  }

  uint64_t vl_elms = sve_cntw(); // vl_elem = 16
  uint64_t m_mod = (((uint64_t)M + vl_elms - 1) / vl_elms) * vl_elms;

  /* Pre-process the left matrix to make it suitable for
     matrix sum of outer-product calculation
   */
  float *A_mod = (float *)malloc(m_mod * M * sizeof(float));

#if defined(UPPER)
  ssymm_direct_sme1_preprocessLU(M, M, A, A_mod);
#elif defined(LOWER)
  ssymm_direct_sme1_preprocessLL(M, M, A, A_mod);
#endif

  /* Calculate C = alpha*A*B + beta*C */
  SGEMM_DIRECT2X2(M, M, N, &alpha, A_mod, B, &beta, R);
  free(A_mod);
}

#else

void CNAME (BLASLONG M, BLASLONG N, float alpha, float * __restrict A,\
            BLASLONG strideA, float * __restrict B, BLASLONG strideB ,\
            float beta, float * __restrict R, BLASLONG strideR){}

#endif
