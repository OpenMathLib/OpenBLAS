/*****************************************************************************
Copyright (c) 2025, The OpenBLAS Project
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
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

/**
 * Unit tests for BGEMM: BFloat16-in, BFloat16-out GEMM.
 *
 * Strategy: compute the same operation via SBGEMM (BF16-in, float32-out) as
 * a trusted reference, then widen the BGEMM BF16 output back to float32 and
 * compare element-wise.  Because both kernels use the same MMA accumulation
 * the only difference is the store path; relative tolerance is set to allow
 * for the one extra BF16 rounding on the output.
 *
 * BF16 helper: bfloat16_to_float() converts a stored bfloat16 value to
 * float32 for comparison without depending on the BFLOAT16CONVERSION macro.
 */

#include "utest/openblas_utest.h"
#include "common.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifdef BUILD_BFLOAT16

#define DATASIZE 100

/* -----------------------------------------------------------------------
 * BF16 helpers
 * --------------------------------------------------------------------- */

/* Widen a bfloat16 bit-pattern to float32 by zero-extending the low 16 bits.
 * This matches the OpenBLAS BF16TOF32 macro on both endiannesses because the
 * bfloat16 typedef stores the *high* 16 bits of a float32 regardless of host
 * endianness. */
static inline float bf16_to_float(bfloat16 b)
{
    float f = 0.0f;
    unsigned short *pf = (unsigned short *)&f;
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    pf[0] = (unsigned short)b;
#else
    pf[1] = (unsigned short)b;
#endif
    return f;
}

/* Truncate float32 to bfloat16 (round-to-nearest, ties-to-even omitted for
 * simplicity — matching what sbstobf16_ produces from a float already
 * rounded to BF16 representable values in the test data). */
static inline bfloat16 float_to_bf16(float f)
{
    unsigned short *p = (unsigned short *)&f;
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    return (bfloat16)p[0];
#else
    return (bfloat16)p[1];
#endif
}

/* -----------------------------------------------------------------------
 * Test data (static to avoid stack pressure for DATASIZE=100)
 * --------------------------------------------------------------------- */

struct DATA_BGEMM {
    bfloat16 a_bf16[DATASIZE * DATASIZE];   /* BF16 input A */
    bfloat16 b_bf16[DATASIZE * DATASIZE];   /* BF16 input B */
    bfloat16 c_bgemm[DATASIZE * DATASIZE];  /* BGEMM BF16 output  */
    float    c_sbgemm[DATASIZE * DATASIZE]; /* SBGEMM float32 reference */
    bfloat16 c_init[DATASIZE * DATASIZE];   /* saved initial C (BF16) */
    float    c_init_f[DATASIZE * DATASIZE]; /* saved initial C (float32) */
};

static struct DATA_BGEMM data_bgemm;

/* -----------------------------------------------------------------------
 * Core check helper
 *
 * 1. Fill A, B, C with random BF16 values.
 * 2. Run SBGEMM (trusted, float32 C) with float32 alpha/beta.
 * 3. Run BGEMM  (tested,  BF16   C) with BF16   alpha/beta.
 * 4. Widen BGEMM output to float32 and diff against SBGEMM reference.
 * Returns the max absolute difference (normalised by |ref| + 1).
 * --------------------------------------------------------------------- */
static float check_bgemm(char transa, char transb,
                          blasint m, blasint n, blasint k,
                          float alpha_f, float beta_f,
                          blasint lda, blasint ldb, blasint ldc)
{
    blasint i, j;
    float max_rdiff = 0.0f;

    /* Convert scalar parameters to BF16 for BGEMM */
    bfloat16 alpha_bf16 = float_to_bf16(alpha_f);
    bfloat16 beta_bf16  = float_to_bf16(beta_f);

    /* Widen back to float so both kernels use the same representable value */
    float alpha_rep = bf16_to_float(alpha_bf16);
    float beta_rep  = bf16_to_float(beta_bf16);

    /* Determine A and B sizes based on transpose flags */
    blasint a_rows = (transa == 'N' || transa == 'n') ? m : k;
    blasint a_cols = (transa == 'N' || transa == 'n') ? k : m;
    blasint b_rows = (transb == 'N' || transb == 'n') ? k : n;
    blasint b_cols = (transb == 'N' || transb == 'n') ? n : k;
    (void)b_rows; /* unused, kept for clarity */
    (void)a_rows;

    /* Fill A and B with random BF16 values in [0.5, 1.5] */
    for (i = 0; i < a_cols * lda; i++)
        data_bgemm.a_bf16[i] = float_to_bf16(0.5f + (float)rand() / (float)RAND_MAX);
    for (i = 0; i < b_cols * ldb; i++)
        data_bgemm.b_bf16[i] = float_to_bf16(0.5f + (float)rand() / (float)RAND_MAX);

    /* Fill C with random BF16 values, saving copies for both kernels */
    for (i = 0; i < n * ldc; i++) {
        data_bgemm.c_init[i]  = float_to_bf16(0.5f + (float)rand() / (float)RAND_MAX);
        data_bgemm.c_init_f[i] = bf16_to_float(data_bgemm.c_init[i]);
    }

    /* Copy initial C into both output buffers */
    memcpy(data_bgemm.c_bgemm,  data_bgemm.c_init,   n * ldc * sizeof(bfloat16));
    memcpy(data_bgemm.c_sbgemm, data_bgemm.c_init_f,  n * ldc * sizeof(float));

    /* Reference: SBGEMM (BF16 A/B, float32 C) */
    BLASFUNC(sbgemm)(&transa, &transb, &m, &n, &k,
                     &alpha_rep,
                     data_bgemm.a_bf16, &lda,
                     data_bgemm.b_bf16, &ldb,
                     &beta_rep,
                     data_bgemm.c_sbgemm, &ldc);

    /* Tested: BGEMM (BF16 A/B/C/alpha/beta) */
    BLASFUNC(bgemm)(&transa, &transb, &m, &n, &k,
                    &alpha_bf16,
                    data_bgemm.a_bf16, &lda,
                    data_bgemm.b_bf16, &ldb,
                    &beta_bf16,
                    data_bgemm.c_bgemm, &ldc);

    /* Compare: widen BGEMM BF16 output to float32, diff vs SBGEMM float32 */
    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            float ref  = data_bgemm.c_sbgemm[j * ldc + i];
            float test = bf16_to_float(data_bgemm.c_bgemm[j * ldc + i]);
            float rdiff = fabsf(test - ref) / (fabsf(ref) + 1.0f);
            if (rdiff > max_rdiff)
                max_rdiff = rdiff;
        }
    }
    return max_rdiff;
}

/*
 * Tolerance: SBGEMM accumulates in float32, stores float32.
 * BGEMM does the same accumulation then converts float32 → BF16 on store,
 * introducing at most 1 BF16 ULP (≈ 2^-7 ≈ 0.0078).  Use 0.01 to match
 * the tolerance used in compare_sgemm_sbgemm.c.
 */
#define BGEMM_TOL 0.01f

/* =======================================================================
 * Test cases: NN, TN, NT, TT × a few representative sizes + edge cases
 * ===================================================================== */

/* ----- NN ----- */
CTEST(bgemm, nn_M50_N50_K50)
{
    float rdiff = check_bgemm('N', 'N', 50, 50, 50,
                               1.5f, 0.0f, 50, 50, 50);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

CTEST(bgemm, nn_M100_N100_K100)
{
    float rdiff = check_bgemm('N', 'N', 100, 100, 100,
                               1.0f, 0.0f, 100, 100, 100);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

/* Hits the 16x8 primary tile exactly (kernel UNROLL_M=16, UNROLL_N=8) */
CTEST(bgemm, nn_M16_N8_K16)
{
    float rdiff = check_bgemm('N', 'N', 16, 8, 16,
                               1.0f, 0.0f, 16, 16, 16);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

/* m and n not multiples of tile — exercises all remainder paths */
CTEST(bgemm, nn_remainder_M17_N9_K7)
{
    float rdiff = check_bgemm('N', 'N', 17, 9, 7,
                               1.0f, 0.0f, 17, 7, 17);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

/* m=1 and n=1 — exercises the scalar (m&1) and (n&1) remainder paths */
CTEST(bgemm, nn_M1_N1_K8)
{
    float rdiff = check_bgemm('N', 'N', 1, 1, 8,
                               1.0f, 0.0f, 1, 8, 1);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

/* ----- TN ----- */
CTEST(bgemm, tn_M50_N50_K50)
{
    float rdiff = check_bgemm('T', 'N', 50, 50, 50,
                               1.5f, 0.0f, 50, 50, 50);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

CTEST(bgemm, tn_M100_N50_K50)
{
    float rdiff = check_bgemm('T', 'N', 100, 50, 50,
                               1.0f, 0.0f, 50, 50, 100);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

/* ----- NT ----- */
CTEST(bgemm, nt_M50_N50_K100)
{
    float rdiff = check_bgemm('N', 'T', 50, 50, 100,
                               1.0f, 0.0f, 50, 50, 50);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

/* ----- TT ----- */
CTEST(bgemm, tt_M50_N50_K50)
{
    float rdiff = check_bgemm('T', 'T', 50, 50, 50,
                               1.5f, 2.0f, 50, 50, 50);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

/* ----- alpha = 0: C should equal beta * C_init ----- */
CTEST(bgemm, alpha_zero_nn)
{
    float rdiff = check_bgemm('N', 'N', 50, 50, 50,
                               0.0f, 2.0f, 50, 50, 50);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

/* ----- beta = 0: C should be fresh accumulation, no add from old C ----- */
CTEST(bgemm, beta_zero_nn)
{
    float rdiff = check_bgemm('N', 'N', 50, 50, 50,
                               1.0f, 0.0f, 50, 50, 50);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

/* ----- beta = 1: accumulate into C, no scaling on existing C ----- */
CTEST(bgemm, beta_one_nn)
{
    float rdiff = check_bgemm('N', 'N', 50, 50, 50,
                               2.0f, 1.0f, 50, 50, 50);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

/* ----- odd k: exercises the k%2==1 remainder path in the kernel ----- */
CTEST(bgemm, nn_odd_k_M32_N8_K7)
{
    float rdiff = check_bgemm('N', 'N', 32, 8, 7,
                               1.0f, 0.0f, 32, 7, 32);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

CTEST(bgemm, nn_odd_k_M16_N8_K1)
{
    float rdiff = check_bgemm('N', 'N', 16, 8, 1,
                               1.0f, 0.0f, 16, 1, 16);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

/* ----- Large n>8 sweep (n&4, n&2, n&1 paths) ----- */
CTEST(bgemm, nn_n4_path_M32_N4_K16)
{
    float rdiff = check_bgemm('N', 'N', 32, 4, 16,
                               1.0f, 0.0f, 32, 16, 32);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

CTEST(bgemm, nn_n2_path_M32_N2_K16)
{
    float rdiff = check_bgemm('N', 'N', 32, 2, 16,
                               1.0f, 0.0f, 32, 16, 32);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

CTEST(bgemm, nn_n1_path_M32_N1_K16)
{
    float rdiff = check_bgemm('N', 'N', 32, 1, 16,
                               1.0f, 0.0f, 32, 16, 32);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

/* ----- m remainder paths: m&8, m&4, m&2, m&1 within the n>=8 loop ----- */
CTEST(bgemm, nn_m8_path_M8_N8_K16)
{
    float rdiff = check_bgemm('N', 'N', 8, 8, 16,
                               1.0f, 0.0f, 8, 16, 8);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

CTEST(bgemm, nn_m4_path_M4_N8_K16)
{
    float rdiff = check_bgemm('N', 'N', 4, 8, 16,
                               1.0f, 0.0f, 4, 16, 4);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

CTEST(bgemm, nn_m2_path_M2_N8_K16)
{
    float rdiff = check_bgemm('N', 'N', 2, 8, 16,
                               1.0f, 0.0f, 2, 16, 2);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

CTEST(bgemm, nn_m1_path_M1_N8_K16)
{
    float rdiff = check_bgemm('N', 'N', 1, 8, 16,
                               1.0f, 0.0f, 1, 16, 1);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

/* ----- m=32 sweep hits the m>=32 fast path (2x16 tile) ----- */
CTEST(bgemm, nn_m32_N8_K16)
{
    float rdiff = check_bgemm('N', 'N', 32, 8, 16,
                               1.0f, 0.0f, 32, 16, 32);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

/* ----- alpha and beta both BF16-representable non-trivial values ----- */
CTEST(bgemm, nn_alpha1p5_beta2_M50_N50_K50)
{
    /* alpha=1.5, beta=2.0 are exactly BF16-representable */
    float rdiff = check_bgemm('N', 'N', 50, 50, 50,
                               1.5f, 2.0f, 50, 50, 50);
    ASSERT_DBL_NEAR_TOL(0.0f, rdiff, BGEMM_TOL);
}

#endif /* BUILD_BFLOAT16 */
