/*****************************************************************************
Copyright (c) 2023, The OpenBLAS Project
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
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

**********************************************************************************/

#include "utest/openblas_utest.h"
#include <cblas.h>
#include "common.h"

#define DATASIZE 100
#define INCREMENT 2

struct DATA_CGEMM {
	float a_test[DATASIZE * DATASIZE * 2];
    float a_verify[DATASIZE * DATASIZE * 2];
	float b_test[DATASIZE * DATASIZE * 2];
    float b_verify[DATASIZE * DATASIZE * 2];
    float c_test[DATASIZE * DATASIZE * 2];
	float c_verify[DATASIZE * DATASIZE * 2];
};

#ifdef BUILD_COMPLEX
static struct DATA_CGEMM data_cgemm;

/**
 * Test cgemm with the conjugate matrices by conjugating and not transposed matrices
 * and comparing it with the non-conjugate cgemm.
 *
 * param transa specifies op(A), the transposition (conjugation) operation applied to A
 * param transb specifies op(B), the transposition (conjugation) operation applied to B
 * param m specifies the number of rows of the matrix op(A) and of the matrix C
 * param n specifies the number of columns of the matrix op(B) and the number of columns of the matrix C
 * param k specifies the number of columns of the matrix op(A) and the number of rows of the matrix op(B)
 * param alpha - scaling factor for the matrix-matrix product
 * param lda - leading dimension of matrix A
 * param ldb - leading dimension of matrix B
 * param beta - scaling factor for matrix C
 * param ldc - leading dimension of matrix C
 * return norm of difference
 */
static float check_cgemm(char transa, char transb, blasint m, blasint n, blasint k, 
                         float *alpha, blasint lda, blasint ldb, float *beta, blasint ldc)
{
	blasint i;
	float alpha_conj[] = {1.0f, 0.0f}; 
	char transa_verify = transa;
    char transb_verify = transb;
    char cc[2]="C", cr[2]="R";

    blasint arows = k, acols = m;
    blasint brows = n, bcols = k;

    if (transa == 'T' || transa == 'C'){
        arows = m; acols = k;
    }

    if (transb == 'T' || transb == 'C'){
        brows = k; bcols = n;
    }

	srand_generate(data_cgemm.a_test, arows * lda * 2);
	srand_generate(data_cgemm.b_test, brows * ldb * 2);
    srand_generate(data_cgemm.c_test, n * ldc * 2);

	for (i = 0; i < arows * lda * 2; i++)
		data_cgemm.a_verify[i] = data_cgemm.a_test[i];

	for (i = 0; i < brows * ldb * 2; i++)
		data_cgemm.b_verify[i] = data_cgemm.b_test[i];

    for (i = 0; i < n * ldc * 2; i++)
		data_cgemm.c_verify[i] = data_cgemm.c_test[i];

	if (transa == 'R'){
		BLASFUNC(cimatcopy)(cc, cr, &arows, &acols, alpha_conj, data_cgemm.a_verify, &lda, &lda);
		transa_verify = 'N';
	}

    if (transb == 'R'){
		BLASFUNC(cimatcopy)(cc, cr, &brows, &bcols, alpha_conj, data_cgemm.b_verify, &ldb, &ldb);
		transb_verify = 'N';
	}

	BLASFUNC(cgemm)(&transa_verify, &transb_verify, &m, &n, &k, alpha, data_cgemm.a_verify, &lda,
	 				data_cgemm.b_verify, &ldb, beta, data_cgemm.c_verify, &ldc);

	BLASFUNC(cgemm)(&transa, &transb, &m, &n, &k, alpha, data_cgemm.a_test, &lda,
	 				data_cgemm.b_test, &ldb, beta, data_cgemm.c_test, &ldc);

	return smatrix_difference(data_cgemm.c_test, data_cgemm.c_verify, m, n, ldc*2);
}

/**
 * Test cgemm with the conjugate matrices by conjugating and not transposed matrices 
 * and comparing it with the non-conjugate cgemm.
 * Test with the following options:
 *
 * matrix A is conjugate and transposed
 * matrix B is conjugate and not transposed
 */
CTEST(cgemm, conjtransa_conjnotransb)
{
	blasint n = DATASIZE, m = DATASIZE, k = DATASIZE;
    blasint lda = DATASIZE, ldb = DATASIZE, ldc = DATASIZE;
	char transa = 'C';
	char transb = 'R';
	float alpha[] = {-2.0, 1.0f};
    float beta[] = {1.0f, -1.0f};

	float norm = check_cgemm(transa, transb, m, n, k, alpha, lda, ldb, beta, ldc);

	ASSERT_DBL_NEAR_TOL(0.0f, norm, SINGLE_EPS);
}

/**
 * Test cgemm with the conjugate matrices by conjugating and not transposed matrices 
 * and comparing it with the non-conjugate cgemm.
 * Test with the following options:
 *
 * matrix A is not conjugate and not transposed
 * matrix B is conjugate and not transposed
 */
CTEST(cgemm, notransa_conjnotransb)
{
	blasint n = DATASIZE, m = DATASIZE, k = DATASIZE;
    blasint lda = DATASIZE, ldb = DATASIZE, ldc = DATASIZE;
	char transa = 'N';
	char transb = 'R';
	float alpha[] = {-2.0, 1.0f};
    float beta[] = {1.0f, -1.0f};

	float norm = check_cgemm(transa, transb, m, n, k, alpha, lda, ldb, beta, ldc);

	ASSERT_DBL_NEAR_TOL(0.0f, norm, SINGLE_EPS);
}

/**
 * Test cgemm with the conjugate matrices by conjugating and not transposed matrices 
 * and comparing it with the non-conjugate cgemm.
 * Test with the following options:
 *
 * matrix A is conjugate and not transposed
 * matrix B is conjugate and transposed
 */
CTEST(cgemm, conjnotransa_conjtransb)
{
	blasint n = DATASIZE, m = DATASIZE, k = DATASIZE;
    blasint lda = DATASIZE, ldb = DATASIZE, ldc = DATASIZE;
	char transa = 'R';
	char transb = 'C';
	float alpha[] = {-2.0, 1.0f};
    float beta[] = {1.0f, -1.0f};

	float norm = check_cgemm(transa, transb, m, n, k, alpha, lda, ldb, beta, ldc);

	ASSERT_DBL_NEAR_TOL(0.0f, norm, SINGLE_EPS);
}

/**
 * Test cgemm with the conjugate matrices by conjugating and not transposed matrices 
 * and comparing it with the non-conjugate cgemm.
 * Test with the following options:
 *
 * matrix A is conjugate and not transposed
 * matrix B is not conjugate and not transposed
 */
CTEST(cgemm, conjnotransa_notransb)
{
	blasint n = DATASIZE, m = DATASIZE, k = DATASIZE;
    blasint lda = DATASIZE, ldb = DATASIZE, ldc = DATASIZE;
	char transa = 'R';
	char transb = 'N';
	float alpha[] = {-2.0, 1.0f};
    float beta[] = {1.0f, -1.0f};

	float norm = check_cgemm(transa, transb, m, n, k, alpha, lda, ldb, beta, ldc);

	ASSERT_DBL_NEAR_TOL(0.0f, norm, SINGLE_EPS);
}

/**
 * Test cgemm with the conjugate matrices by conjugating and not transposed matrices 
 * and comparing it with the non-conjugate cgemm.
 * Test with the following options:
 *
 * matrix A is conjugate and not transposed
 * matrix B is conjugate and not transposed
 */
CTEST(cgemm, conjnotransa_conjnotransb)
{
	blasint n = DATASIZE, m = DATASIZE, k = DATASIZE;
    blasint lda = DATASIZE, ldb = DATASIZE, ldc = DATASIZE;
	char transa = 'R';
	char transb = 'R';
	float alpha[] = {-2.0, 1.0f};
    float beta[] = {1.0f, -1.0f};

	float norm = check_cgemm(transa, transb, m, n, k, alpha, lda, ldb, beta, ldc);

	ASSERT_DBL_NEAR_TOL(0.0f, norm, SINGLE_EPS);
}

/**
 * Test cgemm with the conjugate matrices by conjugating and not transposed matrices 
 * and comparing it with the non-conjugate cgemm.
 * Test with the following options:
 *
 * matrix A is conjugate and not transposed
 * matrix B is transposed
 */
CTEST(cgemm, conjnotransa_transb)
{
	blasint n = DATASIZE, m = DATASIZE, k = DATASIZE;
    blasint lda = DATASIZE, ldb = DATASIZE, ldc = DATASIZE;
	char transa = 'R';
	char transb = 'T';
	float alpha[] = {-2.0, 1.0f};
    float beta[] = {1.0f, -1.0f};

	float norm = check_cgemm(transa, transb, m, n, k, alpha, lda, ldb, beta, ldc);

	ASSERT_DBL_NEAR_TOL(0.0f, norm, SINGLE_EPS);
}

/**
 * Test cgemm with the conjugate matrices by conjugating and not transposed matrices 
 * and comparing it with the non-conjugate cgemm.
 * Test with the following options:
 *
 * matrix A is transposed
 * matrix B is conjugate and not transposed
 */
CTEST(cgemm, transa_conjnotransb)
{
	blasint n = DATASIZE, m = DATASIZE, k = DATASIZE;
    blasint lda = DATASIZE, ldb = DATASIZE, ldc = DATASIZE;
	char transa = 'T';
	char transb = 'R';
	float alpha[] = {-2.0, 1.0f};
    float beta[] = {1.0f, -1.0f};

	float norm = check_cgemm(transa, transb, m, n, k, alpha, lda, ldb, beta, ldc);

	ASSERT_DBL_NEAR_TOL(0.0f, norm, SINGLE_EPS);
}

#ifndef NO_CBLAS
CTEST(cgemm, strided_batch_complex_scalars)
{
    float a[] = {1.0f, 2.0f};
    float b[] = {3.0f, 4.0f};
    float c[] = {5.0f, 6.0f};
    float alpha[] = {2.0f, 3.0f};
    float beta[] = {4.0f, -1.0f};

    cblas_cgemm_batch_strided(CblasColMajor, CblasNoTrans, CblasNoTrans,
                              1, 1, 1, alpha, a, 1, 2, b, 1, 2,
                              beta, c, 1, 2, 1);

    ASSERT_DBL_NEAR_TOL(-14.0, c[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(24.0, c[1], SINGLE_EPS);

    beta[0] = beta[1] = 0.0;
    c[0] = 5.0;
    c[1] = 6.0;
    cblas_cgemm_batch_strided(CblasColMajor, CblasNoTrans, CblasNoTrans,
                              1, 1, 1, alpha, a, 1, 2, b, 1, 2,
                              beta, c, 1, 2, 1);

    ASSERT_DBL_NEAR_TOL(-40.0, c[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(5.0, c[1], SINGLE_EPS);
}

CTEST(cgemm, batch_validation_and_lowercase_trans)
{
    float a[] = {1.0f, 2.0f};
    float b[] = {3.0f, 4.0f};
    float c[] = {0.0f, 0.0f};
    float alpha[] = {1.0f, 0.0f};
    float beta[] = {0.0f, 0.0f};
    const float *a_arr[2] = {a, a};
    const float *b_arr[2] = {b, b};
    float *c_arr[2] = {c, c};
    CBLAS_TRANSPOSE ta_arr[2] = {CblasNoTrans, CblasNoTrans};
    CBLAS_TRANSPOSE tb_arr[2] = {CblasNoTrans, CblasNoTrans};
    blasint m_arr[2] = {1, 1}, n_arr[2] = {1, 1}, k_arr[2] = {1, 1};
    blasint ld_arr[2] = {1, 1};
    blasint bad_gsize[2] = {2, -2};
    char ta_f = 'r', tb_f = 'n';
    blasint one = 1;

    set_xerbla("CGEMM_BATCH ", 0);
    cblas_cgemm_batch((CBLAS_ORDER)INVALID, ta_arr, tb_arr, m_arr, n_arr, k_arr,
                      alpha, (const void **)a_arr, ld_arr, (const void **)b_arr, ld_arr,
                      beta, (void **)c_arr, ld_arr, 1, &one);
    ASSERT_EQUAL(TRUE, check_error());

    set_xerbla("CGEMM_BATCH ", 14);
    cblas_cgemm_batch(CblasColMajor, ta_arr, tb_arr, m_arr, n_arr, k_arr,
                      alpha, (const void **)a_arr, ld_arr, (const void **)b_arr, ld_arr,
                      beta, (void **)c_arr, ld_arr, -1, &one);
    ASSERT_EQUAL(TRUE, check_error());

    set_xerbla("CGEMM_BATCH ", 15);
    cblas_cgemm_batch(CblasColMajor, ta_arr, tb_arr, m_arr, n_arr, k_arr,
                      alpha, (const void **)a_arr, ld_arr, (const void **)b_arr, ld_arr,
                      beta, (void **)c_arr, ld_arr, 2, bad_gsize);
    ASSERT_EQUAL(TRUE, check_error());

    set_xerbla("CGEMM_BATCH_STRIDED ", 11);
    cblas_cgemm_batch_strided(CblasColMajor, CblasNoTrans, CblasNoTrans,
                              2, 2, 2, alpha, a, 2, 4, b, 1, 4, beta, c, 2, 4, 1);
    ASSERT_EQUAL(TRUE, check_error());

    set_xerbla("CGEMM_BATCH_STRIDED ", 15);
    cblas_cgemm_batch_strided(CblasColMajor, CblasNoTrans, CblasNoTrans,
                              2, 2, 2, alpha, a, 2, 4, b, 2, 4, beta, c, 1, 4, 1);
    ASSERT_EQUAL(TRUE, check_error());

    set_xerbla("CGEMM_BATCH_STRIDED ", 17);
    cblas_cgemm_batch_strided(CblasColMajor, CblasNoTrans, CblasNoTrans,
                              1, 1, 1, alpha, a, 1, 2, b, 1, 2, beta, c, 1, 2, -1);
    ASSERT_EQUAL(TRUE, check_error());

    /* conj(1 + 2i) * (3 + 4i) = (1 - 2i)(3 + 4i) = 11 - 2i */
    BLASFUNC(cgemm_batch)(&ta_f, &tb_f, &one, &one, &one, alpha,
                          (float **)a_arr, &one, (float **)b_arr, &one,
                          beta, c_arr, &one, &one, &one);
    ASSERT_DBL_NEAR_TOL(11.0, c[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(-2.0, c[1], SINGLE_EPS);

    c[0] = 0.0f;
    c[1] = 0.0f;
    BLASFUNC(cgemm_batch_strided)(&ta_f, &tb_f, &one, &one, &one, alpha,
                                  a, &one, &one, b, &one, &one,
                                  beta, c, &one, &one, &one);
    ASSERT_DBL_NEAR_TOL(11.0, c[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(-2.0, c[1], SINGLE_EPS);
}
#endif
#endif

#if defined(BUILD_BFLOAT16) && !defined(NO_CBLAS)
CTEST(cgemm, sbgemm_batch_and_sbgemv_order)
{
    /* 0x3f80 is 1.0f in bfloat16, 0x4000 is 2.0f in bfloat16 */
    bfloat16 a[2] = {0x3f80, 0x4000};
    bfloat16 b[2] = {0x4000, 0x3f80};
    float c[2] = {0.0f, 0.0f};
    const bfloat16 *a_arr[1] = {a};
    const bfloat16 *b_arr[1] = {b};
    float *c_arr[1] = {c};
    CBLAS_TRANSPOSE ta = CblasNoTrans, tb = CblasNoTrans;
    blasint m = 1, n = 1, k = 2, ld_a = 1, ld_b = 2, ld_c = 1, gsize = 1;
    float alpha = 1.0f, beta = 0.0f;

    cblas_sbgemm_batch(CblasColMajor, &ta, &tb, &m, &n, &k,
                       &alpha, a_arr, &ld_a, b_arr, &ld_b,
                       &beta, c_arr, &ld_c, 1, &gsize);
    /* 1.0 * 2.0 + 2.0 * 1.0 = 4.0 */
    ASSERT_DBL_NEAR_TOL(4.0, c[0], SINGLE_EPS);

    c[0] = 0.0f;
    cblas_sbgemm_batch_strided(CblasColMajor, CblasNoTrans, CblasNoTrans,
                               1, 1, 2, 1.0f, a, 1, 2, b, 2, 2, 0.0f, c, 1, 1, 1);
    ASSERT_DBL_NEAR_TOL(4.0, c[0], SINGLE_EPS);

    set_xerbla("SBGEMM_BATCH ", 14);
    cblas_sbgemm_batch(CblasColMajor, &ta, &tb, &m, &n, &k,
                       &alpha, a_arr, &ld_a, b_arr, &ld_b,
                       &beta, c_arr, &ld_c, -1, &gsize);
    ASSERT_EQUAL(TRUE, check_error());

    set_xerbla("SBGEMV ", 0);
    cblas_sbgemv((CBLAS_ORDER)INVALID, CblasNoTrans, 1, 1, 1.0f, a, 1, b, 1, 0.0f, c, 1);
    ASSERT_EQUAL(TRUE, check_error());
}
#endif
