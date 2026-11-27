/*****************************************************************************
 * LAPACKE test wrappers (double precision, spike subset)
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
 * DGETRF( M, N, A, LDA, IPIV, INFO )
 *****************************************************************************/
//#define DGETRF_TEST LAPACK_GLOBAL_SUFFIX(dgetrf_test, DGETRF_TEST)
void DGETRF_TEST(const lapack_int *m, const lapack_int *n, double *a,
                 const lapack_int *lda, lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DGETRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgetrf)(layout, *m, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_dgetrf_work)(layout, *m, *n, a_r, lda_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DGETRF", ret);
}

/*****************************************************************************
 * DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 *****************************************************************************/
//#define DGETRS_TEST LAPACK_GLOBAL_SUFFIX(dgetrs_test, DGETRS_TEST)
void DGETRS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
                 const double *a, const lapack_int *lda, const lapack_int *ipiv,
                 double *b, const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const double *a_r = a;
    double *b_r = b;
    lapack_int lda_r = *lda;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DGETRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgetrs)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                     ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dgetrs_work)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                          ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DGETRS", ret);
}

/*****************************************************************************
 * DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
 *****************************************************************************/
//#define DGETRI_TEST LAPACK_GLOBAL_SUFFIX(dgetri_test, DGETRI_TEST)
void DGETRI_TEST(const lapack_int *n, double *a, const lapack_int *lda,
                 const lapack_int *ipiv, double *work, const lapack_int *lwork,
                 lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    if (*lwork == -1) {
        ret = API_SUFFIX(LAPACKE_dgetri_work)(LAPACK_COL_MAJOR, *n, a, *lda,
                                              ipiv, work, *lwork);
        *info = lapacke_test_info("DGETRI", ret);
        return;
    }
#endif

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DGETRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgetri)(layout, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_dgetri_work)(layout, *n, a_r, lda_r, ipiv, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DGETRI", ret);
}

/*****************************************************************************
 * DPOTRF( UPLO, N, A, LDA, INFO )
 *****************************************************************************/
//#define DPOTRF_TEST LAPACK_GLOBAL_SUFFIX(dpotrf_test, DPOTRF_TEST)
void DPOTRF_TEST(const char *uplo, const lapack_int *n, double *a,
                 const lapack_int *lda, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DPOTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpotrf)(layout, *uplo, *n, a_r, lda_r);
#else
    ret = API_SUFFIX(LAPACKE_dpotrf_work)(layout, *uplo, *n, a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dpo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DPOTRF", ret);
}

/*****************************************************************************
 * DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
 *****************************************************************************/
//#define DPOTRS_TEST LAPACK_GLOBAL_SUFFIX(dpotrs_test, DPOTRS_TEST)
void DPOTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const double *a, const lapack_int *lda, double *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const double *a_r = a;
    double *b_r = b;
    lapack_int lda_r = *lda;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DPOTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dpotrs)(layout, *uplo, *n, *nrhs, a_r, lda_r, b_r,
                                     ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dpotrs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DPOTRS", ret);
}

/*****************************************************************************
 * DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 *****************************************************************************/
//#define DGEQRF_TEST LAPACK_GLOBAL_SUFFIX(dgeqrf_test, DGEQRF_TEST)
void DGEQRF_TEST(const lapack_int *m, const lapack_int *n, double *a,
                 const lapack_int *lda, double *tau, double *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    if (*lwork == -1) {
        ret = API_SUFFIX(LAPACKE_dgeqrf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("DGEQRF", ret);
        return;
    }
#endif

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DGEQRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgeqrf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_dgeqrf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DGEQRF", ret);
}

/*****************************************************************************
 * DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 *****************************************************************************/
//#define DORGQR_TEST LAPACK_GLOBAL_SUFFIX(dorgqr_test, DORGQR_TEST)
void DORGQR_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 double *a, const lapack_int *lda, const double *tau,
                 double *work, const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    if (*lwork == -1) {
        ret = API_SUFFIX(LAPACKE_dorgqr_work)(LAPACK_COL_MAJOR, *m, *n, *k, a,
                                              *lda, tau, work, *lwork);
        *info = lapacke_test_info("DORGQR", ret);
        return;
    }
#endif

    double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("DORGQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dorgqr)(layout, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_dorgqr_work)(layout, *m, *n, *k, a_r, lda_r, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("DORGQR", ret);
}

/*****************************************************************************
 * DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
 *****************************************************************************/
//#define DORMQR_TEST LAPACK_GLOBAL_SUFFIX(dormqr_test, DORMQR_TEST)
void DORMQR_TEST(const char *side, const char *trans, const lapack_int *m,
                 const lapack_int *n, const lapack_int *k, const double *a,
                 const lapack_int *lda, const double *tau, double *c,
                 const lapack_int *ldc, double *work, const lapack_int *lwork,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN side_len, FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    if (*lwork == -1) {
        ret = API_SUFFIX(LAPACKE_dormqr_work)(LAPACK_COL_MAJOR, *side, *trans,
                                              *m, *n, *k, a, *lda, tau, c, *ldc,
                                              work, *lwork);
        *info = lapacke_test_info("DORMQR", ret);
        return;
    }
#endif

    const double *a_r = a;
    double *c_r = c;
    lapack_int lda_r = *lda;
    lapack_int ldc_r = *ldc;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int r = API_SUFFIX(LAPACKE_lsame)(*side, 'l') ? *m : *n;
    a_r = lapacke_test_dge_cm_to_rm(r, *k, a, *lda, &lda_r);
    c_r = lapacke_test_dge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    if (a_r == NULL || c_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("DORMQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dormqr)(layout, *side, *trans, *m, *n, *k, a_r,
                                     lda_r, tau, c_r, ldc_r);
#else
    ret = API_SUFFIX(LAPACKE_dormqr_work)(layout, *side, *trans, *m, *n, *k,
                                          a_r, lda_r, tau, c_r, ldc_r, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("DORMQR", ret);
}

/*****************************************************************************
 * DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
 *****************************************************************************/
//#define DGELS_TEST LAPACK_GLOBAL_SUFFIX(dgels_test, DGELS_TEST)
void DGELS_TEST(const char *trans, const lapack_int *m, const lapack_int *n,
                const lapack_int *nrhs, double *a, const lapack_int *lda,
                double *b, const lapack_int *ldb, double *work,
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
        ret = API_SUFFIX(LAPACKE_dgels_work)(LAPACK_COL_MAJOR, *trans, *m, *n,
                                             *nrhs, a, *lda, b, *ldb, work,
                                             *lwork);
        *info = lapacke_test_info("DGELS", ret);
        return;
    }
#endif

    double *a_r = a;
    double *b_r = b;
    lapack_int lda_r = *lda;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int brows = MAX(*m, *n);
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_dge_cm_to_rm(brows, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("DGELS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_dgels)(layout, *trans, *m, *n, *nrhs, a_r, lda_r,
                                    b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_dgels_work)(layout, *trans, *m, *n, *nrhs, a_r,
                                         lda_r, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_dge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_dge_rm_to_cm(brows, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("DGELS", ret);
}

/*****************************************************************************
 * DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
 *****************************************************************************/
//#define DLANGE_TEST LAPACK_GLOBAL_SUFFIX(dlange_test, DLANGE_TEST)
double DLANGE_TEST(const char *norm, const lapack_int *m, const lapack_int *n,
                   const double *a, const lapack_int *lda, double *work
#ifdef LAPACK_FORTRAN_STRLEN_END
                   ,
                   FORTRAN_STRLEN norm_len
#endif
)
{
    double res = 0;
    const double *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_dge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("DLANGE", &info);
        return (double)info;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    res = API_SUFFIX(LAPACKE_dlange)(layout, *norm, *m, *n, a_r, lda_r);
#else
    res = API_SUFFIX(LAPACKE_dlange_work)(layout, *norm, *m, *n, a_r, lda_r,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    return res;
}
