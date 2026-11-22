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

struct DATA_ZGEMM {
	double a_test[DATASIZE * DATASIZE * 2];
    double a_verify[DATASIZE * DATASIZE * 2];
	double b_test[DATASIZE * DATASIZE * 2];
    double b_verify[DATASIZE * DATASIZE * 2];
    double c_test[DATASIZE * DATASIZE * 2];
	double c_verify[DATASIZE * DATASIZE * 2];
};

#ifdef BUILD_COMPLEX16
static struct DATA_ZGEMM data_zgemm;

/**
 * Test zgemm with the conjugate matrices by conjugating and not transposed matrices
 * and comparing it with the non-conjugate zgemm.
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
static double check_zgemm(char transa, char transb, blasint m, blasint n, blasint k, 
                         double *alpha, blasint lda, blasint ldb, double *beta, blasint ldc)
{
	blasint i;
	double alpha_conj[] = {1.0, 0.0}; 
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

	drand_generate(data_zgemm.a_test, arows * lda * 2);
	drand_generate(data_zgemm.b_test, brows * ldb * 2);
    drand_generate(data_zgemm.c_test, n * ldc * 2);

	for (i = 0; i < arows * lda * 2; i++)
		data_zgemm.a_verify[i] = data_zgemm.a_test[i];

	for (i = 0; i < brows * ldb * 2; i++)
		data_zgemm.b_verify[i] = data_zgemm.b_test[i];

    for (i = 0; i < n * ldc * 2; i++)
		data_zgemm.c_verify[i] = data_zgemm.c_test[i];

	if (transa == 'R'){
		BLASFUNC(zimatcopy)(cc, cr, &arows, &acols, alpha_conj, data_zgemm.a_verify, &lda, &lda);
		transa_verify = 'N';
	}

    if (transb == 'R'){
		BLASFUNC(zimatcopy)(cc, cr, &brows, &bcols, alpha_conj, data_zgemm.b_verify, &ldb, &ldb);
	    transb_verify = 'N';
	}

	BLASFUNC(zgemm)(&transa_verify, &transb_verify, &m, &n, &k, alpha, data_zgemm.a_verify, &lda,
	 				data_zgemm.b_verify, &ldb, beta, data_zgemm.c_verify, &ldc);

	BLASFUNC(zgemm)(&transa, &transb, &m, &n, &k, alpha, data_zgemm.a_test, &lda,
	 				data_zgemm.b_test, &ldb, beta, data_zgemm.c_test, &ldc);

	return dmatrix_difference(data_zgemm.c_test, data_zgemm.c_verify, m, n, ldc*2);
}

/**
 * Test zgemm with the conjugate matrices by conjugating and not transposed matrices 
 * and comparing it with the non-conjugate zgemm.
 * Test with the following options:
 *
 * matrix A is conjugate and transposed
 * matrix B is conjugate and not transposed
 */
CTEST(zgemm, conjtransa_conjnotransb)
{
	blasint n = DATASIZE, m = DATASIZE, k = DATASIZE;
    blasint lda = DATASIZE, ldb = DATASIZE, ldc = DATASIZE;
	char transa = 'C';
	char transb = 'R';
	double alpha[] = {-2.0, 1.0};
    double beta[] = {1.0, -1.0};

	double norm = check_zgemm(transa, transb, m, n, k, alpha, lda, ldb, beta, ldc);

	ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS);
}

/**
 * Test zgemm with the conjugate matrices by conjugating and not transposed matrices 
 * and comparing it with the non-conjugate zgemm.
 * Test with the following options:
 *
 * matrix A is not conjugate and not transposed
 * matrix B is conjugate and not transposed
 */
CTEST(zgemm, notransa_conjnotransb)
{
	blasint n = DATASIZE, m = DATASIZE, k = DATASIZE;
    blasint lda = DATASIZE, ldb = DATASIZE, ldc = DATASIZE;
	char transa = 'N';
	char transb = 'R';
	double alpha[] = {-2.0, 1.0};
    double beta[] = {1.0, -1.0};

	double norm = check_zgemm(transa, transb, m, n, k, alpha, lda, ldb, beta, ldc);

	ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS);
}

/**
 * Test zgemm with the conjugate matrices by conjugating and not transposed matrices 
 * and comparing it with the non-conjugate zgemm.
 * Test with the following options:
 *
 * matrix A is conjugate and not transposed
 * matrix B is conjugate and transposed
 */
CTEST(zgemm, conjnotransa_conjtransb)
{
	blasint n = DATASIZE, m = DATASIZE, k = DATASIZE;
    blasint lda = DATASIZE, ldb = DATASIZE, ldc = DATASIZE;
	char transa = 'R';
	char transb = 'C';
	double alpha[] = {-2.0, 1.0};
    double beta[] = {1.0, -1.0};

	double norm = check_zgemm(transa, transb, m, n, k, alpha, lda, ldb, beta, ldc);

	ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS);
}

/**
 * Test zgemm with the conjugate matrices by conjugating and not transposed matrices 
 * and comparing it with the non-conjugate zgemm.
 * Test with the following options:
 *
 * matrix A is conjugate and not transposed
 * matrix B is not conjugate and not transposed
 */
CTEST(zgemm, conjnotransa_notransb)
{
	blasint n = DATASIZE, m = DATASIZE, k = DATASIZE;
    blasint lda = DATASIZE, ldb = DATASIZE, ldc = DATASIZE;
	char transa = 'R';
	char transb = 'N';
	double alpha[] = {-2.0, 1.0};
    double beta[] = {1.0, -1.0};

	double norm = check_zgemm(transa, transb, m, n, k, alpha, lda, ldb, beta, ldc);

	ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS);
}

/**
 * Test zgemm with the conjugate matrices by conjugating and not transposed matrices 
 * and comparing it with the non-conjugate zgemm.
 * Test with the following options:
 *
 * matrix A is conjugate and not transposed
 * matrix B is conjugate and not transposed
 */
CTEST(zgemm, conjnotransa_conjnotransb)
{
	blasint n = DATASIZE, m = DATASIZE, k = DATASIZE;
    blasint lda = DATASIZE, ldb = DATASIZE, ldc = DATASIZE;
	char transa = 'R';
	char transb = 'R';
	double alpha[] = {-2.0, 1.0};
    double beta[] = {1.0, -1.0};

	double norm = check_zgemm(transa, transb, m, n, k, alpha, lda, ldb, beta, ldc);

	ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS);
}

/**
 * Test zgemm with the conjugate matrices by conjugating and not transposed matrices 
 * and comparing it with the non-conjugate zgemm.
 * Test with the following options:
 *
 * matrix A is conjugate and not transposed
 * matrix B is transposed
 */
CTEST(zgemm, conjnotransa_transb)
{
	blasint n = DATASIZE, m = DATASIZE, k = DATASIZE;
    blasint lda = DATASIZE, ldb = DATASIZE, ldc = DATASIZE;
	char transa = 'R';
	char transb = 'T';
	double alpha[] = {-2.0, 1.0};
    double beta[] = {1.0, -1.0};

	double norm = check_zgemm(transa, transb, m, n, k, alpha, lda, ldb, beta, ldc);

	ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS);
}

/**
 * Test zgemm with the conjugate matrices by conjugating and not transposed matrices 
 * and comparing it with the non-conjugate zgemm.
 * Test with the following options:
 *
 * matrix A is transposed
 * matrix B is conjugate and not transposed
 */
CTEST(zgemm, transa_conjnotransb)
{
	blasint n = DATASIZE, m = DATASIZE, k = DATASIZE;
    blasint lda = DATASIZE, ldb = DATASIZE, ldc = DATASIZE;
	char transa = 'T';
	char transb = 'R';
	double alpha[] = {-2.0, 1.0};
    double beta[] = {1.0, -1.0};

	double norm = check_zgemm(transa, transb, m, n, k, alpha, lda, ldb, beta, ldc);

	ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS);
}

#ifndef NO_CBLAS
CTEST(zgemm, strided_batch_complex_scalars)
{
    double a[] = {1.0, 2.0};
    double b[] = {3.0, 4.0};
    double c[] = {5.0, 6.0};
    double alpha[] = {2.0, 3.0};
    double beta[] = {4.0, -1.0};

    cblas_zgemm_batch_strided(CblasColMajor, CblasNoTrans, CblasNoTrans,
                              1, 1, 1, alpha, a, 1, 2, b, 1, 2,
                              beta, c, 1, 2, 1);

    ASSERT_DBL_NEAR_TOL(-14.0, c[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(24.0, c[1], DOUBLE_EPS);

    beta[0] = beta[1] = 0.0;
    c[0] = 5.0;
    c[1] = 6.0;
    cblas_zgemm_batch_strided(CblasColMajor, CblasNoTrans, CblasNoTrans,
                              1, 1, 1, alpha, a, 1, 2, b, 1, 2,
                              beta, c, 1, 2, 1);

    ASSERT_DBL_NEAR_TOL(-40.0, c[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(5.0, c[1], DOUBLE_EPS);
}

CTEST(zgemm, batch_validation_and_lowercase_trans)
{
    double a[] = {1.0, 2.0};
    double b[] = {3.0, 4.0};
    double c[] = {0.0, 0.0};
    double alpha[] = {1.0, 0.0};
    double beta[] = {0.0, 0.0};
    const double *a_arr[2] = {a, a};
    const double *b_arr[2] = {b, b};
    double *c_arr[2] = {c, c};
    CBLAS_TRANSPOSE ta_arr[2] = {CblasNoTrans, CblasNoTrans};
    CBLAS_TRANSPOSE tb_arr[2] = {CblasNoTrans, CblasNoTrans};
    blasint m_arr[2] = {1, 1}, n_arr[2] = {1, 1}, k_arr[2] = {1, 1};
    blasint ld_arr[2] = {1, 1};
    blasint bad_gsize[2] = {2, -2};
    char ta_f = 'r', tb_f = 'n';
    blasint one = 1;

    set_xerbla("ZGEMM_BATCH ", 0);
    cblas_zgemm_batch((CBLAS_ORDER)INVALID, ta_arr, tb_arr, m_arr, n_arr, k_arr,
                      alpha, (const void **)a_arr, ld_arr, (const void **)b_arr, ld_arr,
                      beta, (void **)c_arr, ld_arr, 1, &one);
    ASSERT_EQUAL(TRUE, check_error());

    set_xerbla("ZGEMM_BATCH ", 14);
    cblas_zgemm_batch(CblasColMajor, ta_arr, tb_arr, m_arr, n_arr, k_arr,
                      alpha, (const void **)a_arr, ld_arr, (const void **)b_arr, ld_arr,
                      beta, (void **)c_arr, ld_arr, -1, &one);
    ASSERT_EQUAL(TRUE, check_error());

    set_xerbla("ZGEMM_BATCH ", 15);
    cblas_zgemm_batch(CblasColMajor, ta_arr, tb_arr, m_arr, n_arr, k_arr,
                      alpha, (const void **)a_arr, ld_arr, (const void **)b_arr, ld_arr,
                      beta, (void **)c_arr, ld_arr, 2, bad_gsize);
    ASSERT_EQUAL(TRUE, check_error());

    set_xerbla("ZGEMM_BATCH_STRIDED ", 11);
    cblas_zgemm_batch_strided(CblasColMajor, CblasNoTrans, CblasNoTrans,
                              2, 2, 2, alpha, a, 2, 4, b, 1, 4, beta, c, 2, 4, 1);
    ASSERT_EQUAL(TRUE, check_error());

    set_xerbla("ZGEMM_BATCH_STRIDED ", 15);
    cblas_zgemm_batch_strided(CblasColMajor, CblasNoTrans, CblasNoTrans,
                              2, 2, 2, alpha, a, 2, 4, b, 2, 4, beta, c, 1, 4, 1);
    ASSERT_EQUAL(TRUE, check_error());

    set_xerbla("ZGEMM_BATCH_STRIDED ", 17);
    cblas_zgemm_batch_strided(CblasColMajor, CblasNoTrans, CblasNoTrans,
                              1, 1, 1, alpha, a, 1, 2, b, 1, 2, beta, c, 1, 2, -1);
    ASSERT_EQUAL(TRUE, check_error());

    /* conj(1 + 2i) * (3 + 4i) = (1 - 2i)(3 + 4i) = 11 - 2i */
    BLASFUNC(zgemm_batch)(&ta_f, &tb_f, &one, &one, &one, alpha,
                          (double **)a_arr, &one, (double **)b_arr, &one,
                          beta, c_arr, &one, &one, &one);
    ASSERT_DBL_NEAR_TOL(11.0, c[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(-2.0, c[1], DOUBLE_EPS);

    c[0] = 0.0;
    c[1] = 0.0;
    BLASFUNC(zgemm_batch_strided)(&ta_f, &tb_f, &one, &one, &one, alpha,
                                  a, &one, &one, b, &one, &one,
                                  beta, c, &one, &one, &one);
    ASSERT_DBL_NEAR_TOL(11.0, c[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(-2.0, c[1], DOUBLE_EPS);
}
#endif
#endif
