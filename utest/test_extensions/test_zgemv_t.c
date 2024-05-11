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
#include "common.h"

#define N 100
#define M 100
#define INCREMENT 2

struct DATA_ZGEMV_T {
    double a_test[N * M * 2];
    double a_verify[N * M * 2];
    double y_test[M * INCREMENT * 2];
    double y_verify[M * INCREMENT * 2];
    double x_test[M * INCREMENT * 2];
    double x_verify[M * INCREMENT * 2];
};

// DOUBLE_EPS_ZGEMV = MAX_VAL * NUMBER OF OPERATIONS * DBL_EPSILON
// DOUBLE_EPS_ZGEMV = 5.0 * O(100 * 100) * 2.2e-16 = 1e-11
#define DOUBLE_EPS_ZGEMV 1e-11

#ifdef BUILD_COMPLEX16
static struct DATA_ZGEMV_T data_zgemv_t;

/**
 * Find product of matrix-vector multiplication
 *
 * param n specifies number of columns of A
 * param m specifies number of rows of A and size of vector x
 * param lda specifies leading dimension of A
 * param inc_x specifies increment of vector x
 */
static void matrix_vector_product(blasint n, blasint m, blasint lda, blasint inc_x)
{
    blasint i;
    blasint one=1;
    double *a_ptr = data_zgemv_t.a_verify;
    double *x_ptr = data_zgemv_t.x_test;
    double *x_res = data_zgemv_t.x_verify;

    openblas_complex_double result;

    for (i = 0; i < n * inc_x; i += inc_x)
    {
        result = BLASFUNC(zdotu)(&lda, a_ptr, &one, x_ptr, &inc_x);
        x_res[0] = CREAL(result);
        x_res[1] = CIMAG(result);
        a_ptr += lda * 2;
        x_res += 2 * inc_x;
    }
}

/**
 * Test zgemv by comparing it against zomatcopy, zaxpby and
 * reference func matrix_vector_product
 *
 * zomatcopy perform operation: op(A)
 * matrix_vector_product perform operation: A*x
 * zaxpby perform operation: alpha*x + beta*y
 *
 * param api specifies tested api (C or Fortran)
 * param order specifies row or column major order
 * param trans specifies op(A), the transposition operation
 * applied to the matrix A
 * param m specifies number of rows of A
 * param n specifies number of columns of A
 * param alpha specifies scalar alpha
 * param lda specifies leading dimension of the matrix A
 * param inc_x specifies increment for vector x
 * param beta specifies scalar beta
 * param inc_y specifies increment for vector y
 * return norm of difference between zgemv and result of reference funcs
 */
static double check_zgemv(char api, char order, char trans, blasint m, blasint n, double *alpha,
                          blasint lda, blasint inc_x, double *beta, blasint inc_y)
{
    blasint i;

    enum CBLAS_ORDER corder;
    enum CBLAS_TRANSPOSE ctrans;

    // Transpose parameters for zomatcopy
    // zgemv_t perform operation on transposed matrix, no need to transpose a_verify
    char trans_copy;
    char ctrans_copy;

    // Param alpha for zomatcopy, scale on alpha perform zaxpby
    double alpha_one[] = {1.0, 0.0};

    memset(data_zgemv_t.x_verify, 0.0, m * inc_x * 2 * sizeof(double));

    // Fill matrix A, vectors x, y
    drand_generate(data_zgemv_t.a_test, lda * n * 2);
    drand_generate(data_zgemv_t.x_test, m * inc_x * 2);
    drand_generate(data_zgemv_t.y_test, m * inc_y * 2);

    // Copy vector y for reference funcs
    for (i = 0; i < m * inc_y * 2; i++)
    {
        data_zgemv_t.y_verify[i] = data_zgemv_t.y_test[i];
    }

    if (api == 'F') {
        if (trans == 'T') trans_copy = 'N';
        if (trans == 'C') trans_copy = 'R';
        if (trans == 'U') trans_copy = 'R';
        if (trans == 'D') trans_copy = 'N';

        // Perform operation: op(A)
        BLASFUNC(zomatcopy)(&order, &trans_copy, &m, &n, alpha_one, 
                            data_zgemv_t.a_test, &lda, data_zgemv_t.a_verify, &lda);

        // Find A*x
        matrix_vector_product(n, m, lda, inc_x);

        // Find conj(x)
        if (trans == 'U' || trans == 'D')
        {
            zconjugate_vector(m, inc_x, data_zgemv_t.x_verify);
        }

        // Find alpha*x+beta*y
        BLASFUNC(zaxpby)(&n, alpha, data_zgemv_t.x_verify, &inc_x, beta, 
                         data_zgemv_t.y_verify, &inc_y);

        BLASFUNC(zgemv)(&trans, &m, &n, alpha, data_zgemv_t.a_test, &lda, 
                        data_zgemv_t.x_test, &inc_x, beta, data_zgemv_t.y_test, &inc_y);
    }
#ifndef NO_CBLAS
    else {
        if (order == 'C') corder = CblasColMajor;
        if (order == 'R') corder = CblasRowMajor;
        if (trans == 'T') {ctrans = CblasTrans; ctrans_copy = (corder == CblasRowMajor) ? CblasTrans : CblasNoTrans;}
        if (trans == 'N') {ctrans = CblasNoTrans; ctrans_copy = (corder == CblasRowMajor) ? CblasNoTrans : CblasTrans;}
        if (trans == 'C') {ctrans = CblasConjTrans; ctrans_copy = (corder == CblasRowMajor) ? CblasConjTrans : CblasConjNoTrans;}
        if (trans == 'R') {ctrans = CblasConjNoTrans; ctrans_copy = (corder == CblasRowMajor) ? CblasConjNoTrans : CblasConjTrans;}

        // Perform operation: op(A)
        cblas_zomatcopy(corder, ctrans_copy, m, n, alpha_one, data_zgemv_t.a_test, lda, data_zgemv_t.a_verify, lda);

        // Find A*x
        matrix_vector_product(n, m, lda, inc_x);

        // Find alpha*x+beta*y
        cblas_zaxpby(n, alpha, data_zgemv_t.x_verify, inc_x, beta, data_zgemv_t.y_verify, inc_y);

        cblas_zgemv(corder, ctrans, m, n, alpha, data_zgemv_t.a_test,
                    lda, data_zgemv_t.x_test, inc_x, beta, data_zgemv_t.y_test, inc_y);
    }
#endif

    // Find the differences between output vector caculated by zgemv and reference funcs
    for (i = 0; i < m * inc_y * 2; i++)
        data_zgemv_t.y_test[i] -= data_zgemv_t.y_verify[i];

    // Find the norm of differences
    return BLASFUNC(dznrm2)(&m, data_zgemv_t.y_test, &inc_y);
}

/**
 * Check if error function was called with expected function name
 * and param info
 *
 * param order specifies row or column major order
 * param trans specifies op(A), the transposition operation
 * applied to the matrix A
 * param m specifies number of rows of A
 * param n specifies number of columns of A
 * param lda specifies leading dimension of the matrix A
 * param inc_x specifies increment for vector x
 * param inc_y specifies increment for vector y
 * param expected_info - expected invalid parameter number
 * return TRUE if everything is ok, otherwise FALSE
 */
static int check_badargs(char order, char trans, blasint m, blasint n,
                         blasint lda, blasint inc_x, blasint inc_y, int expected_info)
{
    double alpha[] = {1.0, 1.0};
    double a[] = {1.0, 1.0};
    double x[] = {1.0, 1.0};
    double beta[] = {1.0, 1.0};
    double y[] = {1.0, 1.0};

    set_xerbla("ZGEMV ", expected_info);

    BLASFUNC(zgemv)(&trans, &m, &n, alpha, a, &lda, x, 
                    &inc_x, beta, y, &inc_y);

    return check_error();
}
#ifndef NO_CBLAS
/**
 * C API specific function
 * Check if error function was called with expected function name
 * and param info
 *
 * param order specifies row or column major order
 * param trans specifies op(A), the transposition operation
 * applied to the matrix A
 * param m specifies number of rows of A
 * param n specifies number of columns of A
 * param lda specifies leading dimension of the matrix A
 * param inc_x specifies increment for vector x
 * param inc_y specifies increment for vector y
 * param expected_info - expected invalid parameter number
 * return TRUE if everything is ok, otherwise FALSE
 */
static int c_api_check_badargs(CBLAS_ORDER corder, CBLAS_TRANSPOSE ctrans, blasint m, blasint n,
                               blasint lda, blasint inc_x, blasint inc_y, int expected_info)
{
    double alpha[] = {1.0, 1.0};
    double a[] = {1.0, 1.0};
    double x[] = {1.0, 1.0};
    double beta[] = {1.0, 1.0};
    double y[] = {1.0, 1.0};

    set_xerbla("ZGEMV ", expected_info);

    cblas_zgemv(corder, ctrans, m, n, alpha, a, lda, x, inc_x, beta, y, inc_y);

    return check_error();
}

/**
 * Fortran API specific test
 * Test zgemv by comparing it against reference
 * with the following options:
 *
 * Column Major
 * Transposition
 * Square matrix
 * inc x = 2, inc y = 1
 * alpha_r = 1.0, alpha_i = 1.0
 * beta_r = 2.0, beta_i = 2.0
 */
CTEST(zgemv, colmajor_trans_col_100_row_100_inc_x_1_y_1)
{
    blasint m = 100, n = 100;
    blasint lda = 100;
    char order = 'C';
    char trans = 'T';

    double alpha[] = {2.0, 1.0};
    double beta[] = {1.0, 2.0};

    blasint inc_x = 1;
    blasint inc_y = 1;

    double norm = check_zgemv('F', order, trans, m, n, alpha, lda,
                              inc_x, beta, inc_y);

    ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS_ZGEMV);
}

/**
 * Fortran API specific test
 * Test zgemv by comparing it against reference
 * with the following options:
 *
 * Column Major
 * Transposition
 * Square matrix
 * inc x = 2, inc y = 1
 * alpha_r = 1.0, alpha_i = 1.0
 * beta_r = 2.0, beta_i = 2.0
 */
CTEST(zgemv, colmajor_trans_col_100_row_100_inc_x_2_y_1)
{
    blasint m = 100, n = 100;
    blasint lda = 100;
    char order = 'C';
    char trans = 'T';

    double alpha[] = {1.0, 1.0};
    double beta[] = {2.0, 2.0};

    blasint inc_x = 2;
    blasint inc_y = 1;

    double norm = check_zgemv('F', order, trans, m, n, alpha, lda,
                              inc_x, beta, inc_y);

    ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS_ZGEMV);
}

/**
 * Fortran API specific test
 * Test zgemv by comparing it against reference
 * with the following options:
 *
 * Column Major
 * Transposition and conjugate A
 * Square matrix
 * inc x = 1, inc y = 1
 * alpha_r = 2.0, alpha_i = 1.0
 * beta_r = 2.0, beta_i = 1.0
 */
CTEST(zgemv, colmajor_conjtrans_col_100_row_100_inc_x_1_y_1)
{
    blasint m = 100, n = 100;
    blasint lda = 100;
    char order = 'C';
    char trans = 'C';

    double alpha[] = {2.0, 1.0};
    double beta[] = {2.0, 1.0};

    blasint inc_x = 1;
    blasint inc_y = 1;

    double norm = check_zgemv('F', order, trans, m, n, alpha, lda,
                              inc_x, beta, inc_y);

    ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS_ZGEMV);
}

/**
 * Fortran API specific test
 * Test zgemv by comparing it against reference
 * with the following options:
 *
 * Column Major
 * Transposition and conjugate A
 * Square matrix
 * inc x = 1, inc y = 2
 * alpha_r = 2.0, alpha_i = 1.0
 * beta_r = 2.0, beta_i = 1.0
 */
CTEST(zgemv, colmajor_conjtrans_col_100_row_100_inc_x_1_y_2)
{
    blasint m = 100, n = 100;
    blasint lda = 100;
    char order = 'C';
    char trans = 'C';

    double alpha[] = {1.0, 1.0};
    double beta[] = {1.0, 1.0};

    blasint inc_x = 1;
    blasint inc_y = 2;

    double norm = check_zgemv('F', order, trans, m, n, alpha, lda,
                              inc_x, beta, inc_y);

    ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS_ZGEMV);
}

/**
 * Fortran API specific test
 * Test zgemv by comparing it against reference
 * with the following options:
 *
 * Column Major
 * Transposition and x conjugate
 * Square matrix
 * inc x = 1, inc y = 1
 * alpha_r = 2.0, alpha_i = 1.0
 * beta_r = 2.0, beta_i = 1.0
 */
CTEST(zgemv, colmajor_trans_x_conj_col_100_row_100_inc_x_1_y_1)
{
    blasint m = 100, n = 100;
    blasint lda = 100;
    char order = 'C';
    char trans = 'U';

    double alpha[] = {1.0, 1.0};
    double beta[] = {1.0, 1.0};

    blasint inc_x = 1;
    blasint inc_y = 1;

    double norm = check_zgemv('F', order, trans, m, n, alpha, lda,
                              inc_x, beta, inc_y);

    ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS_ZGEMV);
}

/**
 * Fortran API specific test
 * Test zgemv by comparing it against reference
 * with the following options:
 *
 * Column Major
 * Transposition and x conjugate
 * Square matrix
 * inc x = 2, inc y = 2
 * alpha_r = 1.0, alpha_i = 2.0
 * beta_r = 1.0, beta_i = 1.0
 */
CTEST(zgemv, colmajor_trans_x_conj_col_100_row_100_inc_x_2_y_2)
{
    blasint m = 100, n = 100;
    blasint lda = 100;
    char order = 'C';
    char trans = 'U';

    double alpha[] = {1.0, 2.0};
    double beta[] = {1.0, 1.0};

    blasint inc_x = 2;
    blasint inc_y = 2;

    double norm = check_zgemv('F', order, trans, m, n, alpha, lda,
                              inc_x, beta, inc_y);

    ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS_ZGEMV);
}

/**
 * Fortran API specific test
 * Test zgemv by comparing it against reference
 * with the following options:
 *
 * Column Major
 * Transposition, conjugate A, conjugate x
 * Square matrix
 * inc x = 2, inc y = 2
 * alpha_r = 2.0, alpha_i = 1.0
 * beta_r = 1.0, beta_i = 2.0
 */
CTEST(zgemv, colmajor_conjtrans_x_conj_col_100_row_100_inc_x_1_y_2)
{
    blasint m = 100, n = 100;
    blasint lda = 100;
    char order = 'C';
    char trans = 'D';

    double alpha[] = {2.0, 1.0};
    double beta[] = {1.0, 2.0};

    blasint inc_x = 1;
    blasint inc_y = 2;

    double norm = check_zgemv('F', order, trans, m, n, alpha, lda,
                              inc_x, beta, inc_y);

    ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS_ZGEMV);
}

/**
 * C API specific test
 * Test zgemv by comparing it against reference
 * with the following options:
 *
 * Column Major
 * Transposition, conjugate A, conjugate x
 * Square matrix
 * inc x = 2, inc y = 1
 * alpha_r = 2.0, alpha_i = 1.0
 * beta_r = 1.0, beta_i = 2.0
 */
CTEST(zgemv, c_api_colmajor_trans_col_100_row_100_inc_x_1_y_1)
{
    blasint m = 100, n = 100;
    blasint lda = 100;
    char order = 'C';
    char trans = 'T';

    double alpha[] = {2.0, 1.0};
    double beta[] = {1.0, 2.0};

    blasint inc_x = 1;
    blasint inc_y = 1;

    double norm = check_zgemv('C', order, trans, m, n, alpha, lda,
                              inc_x, beta, inc_y);

    ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS_ZGEMV);
}

/**
 * C API specific test
 * Test zgemv by comparing it against reference
 * with the following options:
 *
 * Column Major
 * Transposition and conjugate A
 * Square matrix
 * inc x = 2, inc y = 1
 * alpha_r = 1.0, alpha_i = 1.0
 * beta_r = 1.0, beta_i = 2.0
 */
CTEST(zgemv, c_api_colmajor_conjtrans_col_100_row_100_inc_x_1_y_1)
{
    blasint m = 100, n = 100;
    blasint lda = 100;
    char order = 'C';
    char trans = 'C';

    double alpha[] = {1.0, 1.0};
    double beta[] = {1.0, 2.0};

    blasint inc_x = 1;
    blasint inc_y = 1;

    double norm = check_zgemv('C', order, trans, m, n, alpha, lda,
                              inc_x, beta, inc_y);

    ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS_ZGEMV);
}

/**
 * C API specific test
 * Test zgemv by comparing it against reference
 * with the following options:
 *
 * Column Major
 * Transposition and conjugate A
 * Square matrix
 * inc x = 1, inc y = 2
 * alpha_r = 1.0, alpha_i = 1.0
 * beta_r = 1.0, beta_i = 2.0
 */
CTEST(zgemv, c_api_colmajor_conjtrans_col_100_row_100_inc_x_1_y_2)
{
    blasint m = 100, n = 100;
    blasint lda = 100;
    char order = 'C';
    char trans = 'C';

    double alpha[] = {1.0, 1.0};
    double beta[] = {1.0, 2.0};

    blasint inc_x = 1;
    blasint inc_y = 2;

    double norm = check_zgemv('C', order, trans, m, n, alpha, lda,
                              inc_x, beta, inc_y);

    ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS_ZGEMV);
}

/**
 * C API specific test
 * Test zgemv by comparing it against reference
 * with the following options:
 *
 * Row Major
 * Square matrix
 * inc x = 1, inc y = 1
 * alpha_r = 2.0, alpha_i = 1.0
 * beta_r = 1.0, beta_i = 1.0
 */
CTEST(zgemv, c_api_rowmajor_notrans_col_100_row_100_inc_x_1_y_1)
{
    blasint m = 100, n = 100;
    blasint lda = 100;
    char order = 'R';
    char trans = 'N';

    double alpha[] = {2.0, 1.0};
    double beta[] = {1.0, 1.0};

    blasint inc_x = 1;
    blasint inc_y = 1;

    double norm = check_zgemv('C', order, trans, m, n, alpha, lda,
                              inc_x, beta, inc_y);

    ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS_ZGEMV);
}

/**
 * C API specific test
 * Test zgemv by comparing it against reference
 * with the following options:
 *
 * Row Major
 * No trans
 * Square matrix
 * inc x = 2, inc y = 2
 * alpha_r = 1.0, alpha_i = 1.0
 * beta_r = 3.0, beta_i = 2.0
 */
CTEST(zgemv, c_api_rowmajor_notrans_col_100_row_100_inc_x_2_y_2)
{
    blasint m = 100, n = 100;
    blasint lda = 100;
    char order = 'R';
    char trans = 'N';

    double alpha[] = {1.0, 1.0};
    double beta[] = {3.0, 1.0};

    blasint inc_x = 2;
    blasint inc_y = 2;

    double norm = check_zgemv('C', order, trans, m, n, alpha, lda,
                              inc_x, beta, inc_y);

    ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS_ZGEMV);
}

/**
 * C API specific test
 * Test zgemv by comparing it against reference
 * with the following options:
 *
 * Column Major
 * Conjugate
 * Square matrix
 * inc x = 1, inc y = 1
 * alpha_r = 1.0, alpha_i = 3.0
 * beta_r = 1.0, beta_i = 2.5
 */
CTEST(zgemv, c_api_rowmajor_conj_col_100_row_100_inc_x_1_y_1)
{
    blasint m = 100, n = 100;
    blasint lda = 100;
    char order = 'R';
    char trans = 'R';

    double alpha[] = {1.0, 3.0};
    double beta[] = {1.0, 2.5};

    blasint inc_x = 1;
    blasint inc_y = 1;

    double norm = check_zgemv('C', order, trans, m, n, alpha, lda,
                              inc_x, beta, inc_y);

    ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS_ZGEMV);
}

/**
 * C API specific test
 * Test zgemv by comparing it against reference
 * with the following options:
 *
 * Row Major
 * Conjugate
 * Square matrix
 * inc x = 2, inc y = 1
 * alpha_r = 1.0, alpha_i = 1.0
 * beta_r = 1.0, beta_i = 1.5
 */
CTEST(zgemv, c_api_rowmajor_conj_col_100_row_100_inc_x_2_y_1)
{
    blasint m = 100, n = 100;
    blasint lda = 100;
    char order = 'R';
    char trans = 'R';

    double alpha[] = {1.0, 1.0};
    double beta[] = {1.0, 1.5};

    blasint inc_x = 2;
    blasint inc_y = 1;

    double norm = check_zgemv('C', order, trans, m, n, alpha, lda,
                              inc_x, beta, inc_y);

    ASSERT_DBL_NEAR_TOL(0.0, norm, DOUBLE_EPS_ZGEMV);
}

/**
 * Fortran API specific test
 * Test error function for an invalid param inc_y.
 * Must be positive
 *
 * Column major
 */
CTEST(zgemv, xerbla_invalid_inc_y)
{
    char order = 'C';
    char trans = 'T';

    blasint m = 1, n = 1;
    blasint lda = 1;

    blasint inc_x = 1;
    blasint inc_y = 0;

    int expected_info = 11;

    int passed = check_badargs(order, trans, m, n, lda, inc_x, inc_y, expected_info);
    ASSERT_EQUAL(TRUE, passed);
}

/**
 * C API specific test
 * Test error function for an invalid param inc_y.
 * Must be positive
 *
 * Column major
 */
CTEST(zgemv, c_api_xerbla_invalid_inc_y_col_major)
{
    enum CBLAS_ORDER corder = CblasColMajor;
    enum CBLAS_TRANSPOSE ctrans = CblasTrans;

    blasint m = 1, n = 1;
    blasint lda = 1;

    blasint inc_x = 1;
    blasint inc_y = 0;

    int expected_info = 11;

    int passed = c_api_check_badargs(corder, ctrans, m, n, lda, inc_x, inc_y, expected_info);
    ASSERT_EQUAL(TRUE, passed);
}

/**
 * C API specific test
 * Test error function for an invalid param inc_y.
 * Must be positive
 *
 * Row major
 */
CTEST(zgemv, c_api_xerbla_invalid_inc_y_row_major)
{
    enum CBLAS_ORDER corder = CblasRowMajor;
    enum CBLAS_TRANSPOSE ctrans = CblasNoTrans;

    blasint m = 1, n = 1;
    blasint lda = 1;

    blasint inc_x = 1;
    blasint inc_y = 0;

    int expected_info = 11;

    int passed = c_api_check_badargs(corder, ctrans, m, n, lda, inc_x, inc_y, expected_info);
    ASSERT_EQUAL(TRUE, passed);
}

/**
 * Fortran API specific test
 * Test error function for an invalid param inc_x.
 * Must be positive
 *
 * Column major
 */
CTEST(zgemv, xerbla_invalid_inc_x)
{
    char order = 'C';
    char trans = 'T';
    blasint m = 1, n = 1;
    blasint lda = 1;

    blasint inc_x = 0;
    blasint inc_y = 1;

    int expected_info = 8;

    int passed = check_badargs(order, trans, m, n, lda, inc_x, inc_y, expected_info);
    ASSERT_EQUAL(TRUE, passed);
}

/**
 * C API specific test
 * Test error function for an invalid param inc_x.
 * Must be positive
 *
 * Column major
 */
CTEST(zgemv, c_api_xerbla_invalid_inc_x_col_major)
{
    enum CBLAS_ORDER corder = CblasColMajor;
    enum CBLAS_TRANSPOSE ctrans = CblasTrans;

    blasint m = 1, n = 1;
    blasint lda = 1;

    blasint inc_x = 0;
    blasint inc_y = 1;

    int expected_info = 8;

    int passed = c_api_check_badargs(corder, ctrans, m, n, lda, inc_x, inc_y, expected_info);
    ASSERT_EQUAL(TRUE, passed);
}

/**
 * C API specific test
 * Test error function for an invalid param inc_x.
 * Must be positive
 *
 * Row major
 */
CTEST(zgemv, c_api_xerbla_invalid_inc_x_row_major)
{
    enum CBLAS_ORDER corder = CblasColMajor;
    enum CBLAS_TRANSPOSE ctrans = CblasTrans;

    blasint m = 1, n = 1;
    blasint lda = 1;

    blasint inc_x = 0;
    blasint inc_y = 1;

    int expected_info = 8;

    int passed = c_api_check_badargs(corder, ctrans, m, n, lda, inc_x, inc_y, expected_info);
    ASSERT_EQUAL(TRUE, passed);
}

/**
 * Fortran API specific test
 * Test error function for an invalid param n.
 * Must be positive.
 *
 * Column major
 */
CTEST(zgemv, xerbla_invalid_n)
{
    char order = 'C';
    char trans = 'T';

    blasint m = 1, n = INVALID;
    blasint lda = 1;

    blasint inc_x = 1;
    blasint inc_y = 1;

    int expected_info = 3;

    int passed = check_badargs(order, trans, m, n, lda, inc_x, inc_y, expected_info);
    ASSERT_EQUAL(TRUE, passed);
}

/**
 * C API specific test
 * Test error function for an invalid param n.
 * Must be positive.
 *
 * Column major
 */
CTEST(zgemv, c_api_xerbla_invalid_n_col_major)
{
    enum CBLAS_ORDER corder = CblasColMajor;
    enum CBLAS_TRANSPOSE ctrans = CblasTrans;

    blasint m = 1, n = INVALID;
    blasint lda = 1;

    blasint inc_x = 1;
    blasint inc_y = 1;

    int expected_info = 3;

    int passed = c_api_check_badargs(corder, ctrans, m, n, lda, inc_x, inc_y, expected_info);
    ASSERT_EQUAL(TRUE, passed);
}

/**
 * C API specific test
 * Test error function for an invalid param n.
 * Must be positive.
 *
 * Row major
 */
CTEST(zgemv, c_api_xerbla_invalid_n_row_major)
{
    enum CBLAS_ORDER corder = CblasRowMajor;
    enum CBLAS_TRANSPOSE ctrans = CblasNoTrans;

    blasint m = INVALID, n = 1;
    blasint lda = 1;

    blasint inc_x = 1;
    blasint inc_y = 1;

    int expected_info = 3;

    int passed = c_api_check_badargs(corder, ctrans, m, n, lda, inc_x, inc_y, expected_info);
    ASSERT_EQUAL(TRUE, passed);
}

/**
 * Fortran API specific test
 * Test error function for an invalid param m.
 * Must be positive.
 *
 * Column major
 */
CTEST(zgemv, xerbla_invalid_m)
{
    char order = 'C';
    char trans = 'T';

    blasint m = INVALID, n = 1;
    blasint lda = 1;

    blasint inc_x = 1;
    blasint inc_y = 1;

    int expected_info = 2;

    int passed = check_badargs(order, trans, m, n, lda, inc_x, inc_y, expected_info);
    ASSERT_EQUAL(TRUE, passed);
}

/**
 * C API specific test
 * Test error function for an invalid param m.
 * Must be positive.
 *
 * Column major
 */
CTEST(zgemv, c_api_xerbla_invalid_m_col_major)
{
    enum CBLAS_ORDER corder = CblasColMajor;
    enum CBLAS_TRANSPOSE ctrans = CblasTrans;

    blasint m = INVALID, n = 1;
    blasint lda = 1;

    blasint inc_x = 1;
    blasint inc_y = 1;

    int expected_info = 2;

    int passed = c_api_check_badargs(corder, ctrans, m, n, lda, inc_x, inc_y, expected_info);
    ASSERT_EQUAL(TRUE, passed);
}

/**
 * C API specific test
 * Test error function for an invalid param m.
 * Must be positive.
 *
 * Row major
 */
CTEST(zgemv, c_api_xerbla_invalid_m_row_major)
{
    enum CBLAS_ORDER corder = CblasRowMajor;
    enum CBLAS_TRANSPOSE ctrans = CblasNoTrans;

    blasint m = 1, n = INVALID;
    blasint lda = 1;

    blasint inc_x = 1;
    blasint inc_y = 1;

    int expected_info = 2;

    int passed = c_api_check_badargs(corder, ctrans, m, n, lda, inc_x, inc_y, expected_info);
    ASSERT_EQUAL(TRUE, passed);
}

/**
 * Fortran API specific test
 * Test error function for an invalid param lda.
 * lda must be at least n.
 *
 * Column major
 */
CTEST(zgemv, xerbla_invalid_lda)
{
    char order = 'C';
    char trans = 'T';

    blasint m = 1, n = 1;
    blasint lda = INVALID;

    blasint inc_x = 1;
    blasint inc_y = 1;

    int expected_info = 6;

    int passed = check_badargs(order, trans, m, n, lda, inc_x, inc_y, expected_info);
    ASSERT_EQUAL(TRUE, passed);
}

/**
 * C API specific test
 * Test error function for an invalid param lda.
 * If matrices are stored using col major layout,
 * lda must be at least m.
 *
 * Column major
 */
CTEST(zgemv, c_api_xerbla_invalid_lda_col_major)
{
    enum CBLAS_ORDER corder = CblasColMajor;
    enum CBLAS_TRANSPOSE ctrans = CblasTrans;

    blasint m = 1, n = 1;
    blasint lda = INVALID;

    blasint inc_x = 1;
    blasint inc_y = 1;

    int expected_info = 6;

    int passed = c_api_check_badargs(corder, ctrans, m, n, lda, inc_x, inc_y, expected_info);
    ASSERT_EQUAL(TRUE, passed);
}

/**
 * C API specific test
 * Test error function for an invalid param lda.
 * If matrices are stored using col major layout,
 * lda must be at least n.
 *
 * Column major
 */
CTEST(zgemv, c_api_xerbla_invalid_lda_row_major)
{
    enum CBLAS_ORDER corder = CblasRowMajor;
    enum CBLAS_TRANSPOSE ctrans = CblasNoTrans;

    blasint m = 1, n = 1;
    blasint lda = INVALID;

    blasint inc_x = 1;
    blasint inc_y = 1;

    int expected_info = 6;

    int passed = c_api_check_badargs(corder, ctrans, m, n, lda, inc_x, inc_y, expected_info);
    ASSERT_EQUAL(TRUE, passed);
}

/**
 * Fortran API specific test
 * Test error function for an invalid param trans.
 *
 * Column major
 */
CTEST(zgemv, xerbla_invalid_trans)
{
    char order = 'C';
    char trans = 'Z';

    blasint m = 1, n = 1;
    blasint lda = 1;

    blasint inc_x = 1;
    blasint inc_y = 1;

    int expected_info = 1;

    int passed = check_badargs(order, trans, m, n, lda, inc_x, inc_y, expected_info);
    ASSERT_EQUAL(TRUE, passed);
}

/**
 * C API specific test
 * Test error function for an invalid param trans.
 *
 * Column major
 */
CTEST(zgemv, c_api_xerbla_invalid_trans_col_major)
{
    enum CBLAS_ORDER corder = CblasColMajor;
    enum CBLAS_TRANSPOSE ctrans = INVALID;

    blasint m = 1, n = 1;
    blasint lda = 1;

    blasint inc_x = 1;
    blasint inc_y = 1;

    int expected_info = 1;

    int passed = c_api_check_badargs(corder, ctrans, m, n, lda, inc_x, inc_y, expected_info);
    ASSERT_EQUAL(TRUE, passed);
}

/**
 * C API specific test
 * Test error function for an invalid param trans.
 *
 * Row major
 */
CTEST(zgemv, c_api_xerbla_invalid_trans_row_major)
{
    enum CBLAS_ORDER corder = CblasRowMajor;
    enum CBLAS_TRANSPOSE ctrans = INVALID;

    blasint m = 1, n = 1;
    blasint lda = 1;

    blasint inc_x = 1;
    blasint inc_y = 1;

    int expected_info = 1;

    int passed = c_api_check_badargs(corder, ctrans, m, n, lda, inc_x, inc_y, expected_info);
    ASSERT_EQUAL(TRUE, passed);
}

/**
 * C API specific test
 * Test error function for an invalid param order.
 */
CTEST(zgemv, c_api_xerbla_invalid_order_col_major)
{
    enum CBLAS_ORDER corder = INVALID;
    enum CBLAS_TRANSPOSE ctrans = CblasTrans;

    blasint m = 1, n = 1;
    blasint lda = 1;

    blasint inc_x = 1;
    blasint inc_y = 1;

    int expected_info = 0;

    int passed = c_api_check_badargs(corder, ctrans, m, n, lda, inc_x, inc_y, expected_info);
    ASSERT_EQUAL(TRUE, passed);
}
#endif
#endif
