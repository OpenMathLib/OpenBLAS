#include "openblas_utest.h"
#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>

#define LAPACK_ROW_MAJOR               101
blasint LAPACKE_dgesvd( blasint matrix_layout, char jobu, char jobvt,
                           blasint m, blasint n, double* a,
                           blasint lda, double* s, double* u, blasint ldu,
                           double* vt, blasint ldvt, double* superb );
                                                                                 

#define DATASIZE 100

double s[DATASIZE];
double u[DATASIZE*DATASIZE];
double vt[DATASIZE*DATASIZE];
double X[DATASIZE*DATASIZE];
double superb[DATASIZE];
double tmp[DATASIZE*DATASIZE];
double m[DATASIZE*DATASIZE];

CTEST(kernel_regress,skx_avx)
{
#ifdef BUILD_DOUBLE
    double norm;
    int i, j, info;
    srand(0);
	for (i = 0; i < DATASIZE*DATASIZE; i++) {
        m[i] = (rand()+0.0)/RAND_MAX * 10;
        tmp[i] = m[i];
    }

    info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'A', 'A', DATASIZE, DATASIZE, m, DATASIZE,
                        s, u, DATASIZE, vt, DATASIZE, superb);

	for (i = 0; i < DATASIZE; i++) {
	    for (j = 0; j < DATASIZE; j++) {
            u[i*DATASIZE+j] = u[i*DATASIZE+j]*s[j];
        }
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                DATASIZE, DATASIZE, DATASIZE, 1, u, DATASIZE, vt, DATASIZE, 0, X, DATASIZE);

	for (i = 0; i < DATASIZE*DATASIZE; i++) {
        X[i] = X[i] - tmp[i];
    }
    
    norm = cblas_dnrm2(DATASIZE*DATASIZE, X, 1);
    ASSERT_DBL_NEAR_TOL(0.0, norm, 1e-10);
#endif
}

#ifdef BUILD_DOUBLE
CTEST(kernel_regress, dgesv_small_thread_clamp)
{
    blasint n = 25, nrhs = 2, lda = 25, ldb = 25, info = 0;
    blasint ipiv[25];
    double a[25 * 25], b[25 * 2];
    int t_idx, i, j;
    int counts[2] = {1, 2};
#if defined(SMP) && !defined(USE_OPENMP)
    int saved = openblas_get_num_threads();
    int max_procs = openblas_get_num_procs();
#endif

    for (t_idx = 0; t_idx < 2; t_idx++) {
#if defined(SMP) && !defined(USE_OPENMP)
        if (counts[t_idx] > max_procs)
            continue;
        openblas_set_num_threads(counts[t_idx]);
#else
        (void)counts;
        if (t_idx > 0)
            break;
#endif
        for (j = 0; j < n; j++) {
            for (i = 0; i < n; i++) {
                a[j * lda + i] = (i == j) ? 50.0 : 1.0 / (double)(i + j + 2);
            }
            b[j] = 1.0;
            b[ldb + j] = (double)(j + 1);
        }
        info = -1;
        BLASFUNC(dgesv)(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
        ASSERT_EQUAL(0, info);
        for (j = 0; j < n; j++) {
            ASSERT_TRUE(b[j] == b[j]);
            ASSERT_TRUE(b[ldb + j] == b[ldb + j]);
        }
    }
#if defined(SMP) && !defined(USE_OPENMP)
    openblas_set_num_threads(saved);
#endif
}

CTEST(kernel_regress, dlauum_blocked)
{
    blasint n = 160, lda = 160, info = 0;
    double *u_mat = (double *)calloc((size_t)(n * lda), sizeof(double));
    double *l_mat = (double *)calloc((size_t)(n * lda), sizeof(double));
    char uplo_u = 'U', uplo_l = 'L';
    int i, j;

    ASSERT_NOT_NULL(u_mat);
    ASSERT_NOT_NULL(l_mat);
    for (j = 0; j < n; j++) {
        for (i = 0; i <= j; i++) {
            u_mat[j * lda + i] = (i == j) ? 2.0 : 0.01;
            l_mat[i * lda + j] = (i == j) ? 2.0 : 0.01;
        }
    }
    BLASFUNC(dlauum)(&uplo_u, &n, u_mat, &lda, &info);
    ASSERT_EQUAL(0, info);
    ASSERT_DBL_NEAR_TOL(4.0 + (double)(n - 1) * 0.0001, u_mat[0], 1e-10);

    info = -1;
    BLASFUNC(dlauum)(&uplo_l, &n, l_mat, &lda, &info);
    ASSERT_EQUAL(0, info);
    ASSERT_DBL_NEAR_TOL(4.0 + (double)(n - 1) * 0.0001, l_mat[0], 1e-10);

    free(u_mat);
    free(l_mat);
}
#endif

#ifdef BUILD_COMPLEX16
CTEST(kernel_regress, zgesv_small_thread_clamp)
{
    blasint n = 25, nrhs = 1, lda = 25, ldb = 25, info = 0;
    blasint ipiv[25];
    double a[25 * 25 * 2], b[25 * 2];
    int t_idx, i, j;
    int counts[2] = {1, 2};
#if defined(SMP) && !defined(USE_OPENMP)
    int saved = openblas_get_num_threads();
    int max_procs = openblas_get_num_procs();
#endif

    for (t_idx = 0; t_idx < 2; t_idx++) {
#if defined(SMP) && !defined(USE_OPENMP)
        if (counts[t_idx] > max_procs)
            continue;
        openblas_set_num_threads(counts[t_idx]);
#else
        (void)counts;
        if (t_idx > 0)
            break;
#endif
        for (j = 0; j < n; j++) {
            for (i = 0; i < n; i++) {
                a[(j * lda + i) * 2] = (i == j) ? 50.0 : 0.1;
                a[(j * lda + i) * 2 + 1] = (i == j) ? 1.0 : -0.05;
            }
            b[j * 2] = 1.0;
            b[j * 2 + 1] = 0.5;
        }
        info = -1;
        BLASFUNC(zgesv)(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
        ASSERT_EQUAL(0, info);
        for (j = 0; j < n; j++) {
            ASSERT_TRUE(b[j * 2] == b[j * 2]);
            ASSERT_TRUE(b[j * 2 + 1] == b[j * 2 + 1]);
        }
    }
#if defined(SMP) && !defined(USE_OPENMP)
    openblas_set_num_threads(saved);
#endif
}
#endif
