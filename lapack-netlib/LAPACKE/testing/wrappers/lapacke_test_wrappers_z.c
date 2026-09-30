/*****************************************************************************
 * LAPACKE test wrappers (double precision complex, spike subset)
 *
 * Each wrapper exports a Fortran-callable symbol <NAME>_TEST with the exact
 * argument list of the corresponding LAPACK routine and forwards to
 * LAPACKE; see lapacke_test_wrappers.h for the wrapper scheme and the
 * (layout, layer) build combinations.
 *****************************************************************************/

#include <stdio.h>

#include "lapacke_test_wrappers.h"
#define API_SUFFIX(x) x
static const int layout = LAPACKE_TEST_LAYOUT;

/*****************************************************************************
 * ZGETRF( M, N, A, LDA, IPIV, INFO )
 *****************************************************************************/
//#define ZGETRF_TEST LAPACK_GLOBAL_SUFFIX(zgetrf_test, ZGETRF_TEST)
void ZGETRF_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZGETRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgetrf)(layout, *m, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zgetrf_work)(layout, *m, *n, a_r, lda_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZGETRF", ret);
}

/*****************************************************************************
 * ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 *****************************************************************************/
//#define ZGETRS_TEST LAPACK_GLOBAL_SUFFIX(zgetrs_test, ZGETRS_TEST)
void ZGETRS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_double *a, const lapack_int *lda,
                 const lapack_int *ipiv, lapack_complex_double *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_complex_double *b_r = b;
    lapack_int lda_r = *lda;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZGETRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgetrs)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                     ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zgetrs_work)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                          ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZGETRS", ret);
}

/*****************************************************************************
 * ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
 *****************************************************************************/
//#define ZGETRI_TEST LAPACK_GLOBAL_SUFFIX(zgetri_test, ZGETRI_TEST)
void ZGETRI_TEST(const lapack_int *n, lapack_complex_double *a,
                 const lapack_int *lda, const lapack_int *ipiv,
                 lapack_complex_double *work, const lapack_int *lwork,
                 lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    if (*lwork == -1) {
        ret = API_SUFFIX(LAPACKE_zgetri_work)(LAPACK_COL_MAJOR, *n, a, *lda,
                                              ipiv, work, *lwork);
        *info = lapacke_test_info("ZGETRI", ret);
        return;
    }
#endif

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZGETRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgetri)(layout, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_zgetri_work)(layout, *n, a_r, lda_r, ipiv, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZGETRI", ret);
}

/*****************************************************************************
 * ZPOTRF( UPLO, N, A, LDA, INFO )
 *****************************************************************************/
//#define ZPOTRF_TEST LAPACK_GLOBAL_SUFFIX(zpotrf_test, ZPOTRF_TEST)
void ZPOTRF_TEST(const char *uplo, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZPOTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpotrf)(layout, *uplo, *n, a_r, lda_r);
#else
    ret = API_SUFFIX(LAPACKE_zpotrf_work)(layout, *uplo, *n, a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zpo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZPOTRF", ret);
}

/*****************************************************************************
 * ZPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
 *****************************************************************************/
//#define ZPOTRS_TEST LAPACK_GLOBAL_SUFFIX(zpotrs_test, ZPOTRS_TEST)
void ZPOTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_double *a, const lapack_int *lda,
                 lapack_complex_double *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_double *a_r = a;
    lapack_complex_double *b_r = b;
    lapack_int lda_r = *lda;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZPOTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zpotrs)(layout, *uplo, *n, *nrhs, a_r, lda_r, b_r,
                                     ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zpotrs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZPOTRS", ret);
}

/*****************************************************************************
 * ZGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 *****************************************************************************/
//#define ZGEQRF_TEST LAPACK_GLOBAL_SUFFIX(zgeqrf_test, ZGEQRF_TEST)
void ZGEQRF_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_double *a, const lapack_int *lda,
                 lapack_complex_double *tau, lapack_complex_double *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    if (*lwork == -1) {
        ret = API_SUFFIX(LAPACKE_zgeqrf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("ZGEQRF", ret);
        return;
    }
#endif

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZGEQRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgeqrf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_zgeqrf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZGEQRF", ret);
}

/*****************************************************************************
 * ZUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 *****************************************************************************/
//#define ZUNGQR_TEST LAPACK_GLOBAL_SUFFIX(zungqr_test, ZUNGQR_TEST)
void ZUNGQR_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 lapack_complex_double *a, const lapack_int *lda,
                 const lapack_complex_double *tau, lapack_complex_double *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    if (*lwork == -1) {
        ret = API_SUFFIX(LAPACKE_zungqr_work)(LAPACK_COL_MAJOR, *m, *n, *k, a,
                                              *lda, tau, work, *lwork);
        *info = lapacke_test_info("ZUNGQR", ret);
        return;
    }
#endif

    lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("ZUNGQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zungqr)(layout, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_zungqr_work)(layout, *m, *n, *k, a_r, lda_r, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("ZUNGQR", ret);
}

/*****************************************************************************
 * ZUNMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
 *****************************************************************************/
//#define ZUNMQR_TEST LAPACK_GLOBAL_SUFFIX(zunmqr_test, ZUNMQR_TEST)
void ZUNMQR_TEST(const char *side, const char *trans, const lapack_int *m,
                 const lapack_int *n, const lapack_int *k,
                 const lapack_complex_double *a, const lapack_int *lda,
                 const lapack_complex_double *tau, lapack_complex_double *c,
                 const lapack_int *ldc, lapack_complex_double *work,
                 const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN side_len, FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    if (*lwork == -1) {
        ret = API_SUFFIX(LAPACKE_zunmqr_work)(LAPACK_COL_MAJOR, *side, *trans,
                                              *m, *n, *k, a, *lda, tau, c, *ldc,
                                              work, *lwork);
        *info = lapacke_test_info("ZUNMQR", ret);
        return;
    }
#endif

    const lapack_complex_double *a_r = a;
    lapack_complex_double *c_r = c;
    lapack_int lda_r = *lda;
    lapack_int ldc_r = *ldc;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int r = API_SUFFIX(LAPACKE_lsame)(*side, 'l') ? *m : *n;
    a_r = lapacke_test_zge_cm_to_rm(r, *k, a, *lda, &lda_r);
    c_r = lapacke_test_zge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    if (a_r == NULL || c_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("ZUNMQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zunmqr)(layout, *side, *trans, *m, *n, *k, a_r,
                                     lda_r, tau, c_r, ldc_r);
#else
    ret = API_SUFFIX(LAPACKE_zunmqr_work)(layout, *side, *trans, *m, *n, *k,
                                          a_r, lda_r, tau, c_r, ldc_r, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("ZUNMQR", ret);
}

/*****************************************************************************
 * ZGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
 *****************************************************************************/
//#define ZGELS_TEST LAPACK_GLOBAL_SUFFIX(zgels_test, ZGELS_TEST)
void ZGELS_TEST(const char *trans, const lapack_int *m, const lapack_int *n,
                const lapack_int *nrhs, lapack_complex_double *a,
                const lapack_int *lda, lapack_complex_double *b,
                const lapack_int *ldb, lapack_complex_double *work,
                const lapack_int *lwork, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                ,
                FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    if (*lwork == -1) {
        ret = API_SUFFIX(LAPACKE_zgels_work)(LAPACK_COL_MAJOR, *trans, *m, *n,
                                             *nrhs, a, *lda, b, *ldb, work,
                                             *lwork);
        *info = lapacke_test_info("ZGELS", ret);
        return;
    }
#endif

    lapack_complex_double *a_r = a;
    lapack_complex_double *b_r = b;
    lapack_int lda_r = *lda;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int brows = MAX(*m, *n);
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_zge_cm_to_rm(brows, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("ZGELS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_zgels)(layout, *trans, *m, *n, *nrhs, a_r, lda_r,
                                    b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_zgels_work)(layout, *trans, *m, *n, *nrhs, a_r,
                                         lda_r, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_zge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_zge_rm_to_cm(brows, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("ZGELS", ret);
}

/*****************************************************************************
 * DOUBLE PRECISION FUNCTION ZLANGE( NORM, M, N, A, LDA, WORK )
 *****************************************************************************/
//#define ZLANGE_TEST LAPACK_GLOBAL_SUFFIX(zlange_test, ZLANGE_TEST)
double ZLANGE_TEST(const char *norm, const lapack_int *m, const lapack_int *n,
                   const lapack_complex_double *a, const lapack_int *lda,
                   double *work
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN norm_len
#endif
)
{
    double res = 0;
    const lapack_complex_double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_zge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("ZLANGE", &info);
        return (double)info;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    res = API_SUFFIX(LAPACKE_zlange)(layout, *norm, *m, *n, a_r, lda_r);
#else
    res = API_SUFFIX(LAPACKE_zlange_work)(layout, *norm, *m, *n, a_r, lda_r,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    return res;
}
