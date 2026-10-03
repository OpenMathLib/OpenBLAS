/*****************************************************************************
 * LAPACKE test wrappers (single precision, spike subset)
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
 * SGETRF( M, N, A, LDA, IPIV, INFO )
 *****************************************************************************/
//#define SGETRF_TEST LAPACK_GLOBAL_SUFFIX(sgetrf_test, SGETRF_TEST)
void SGETRF_TEST(const lapack_int *m, const lapack_int *n, float *a,
                 const lapack_int *lda, lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SGETRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgetrf)(layout, *m, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_sgetrf_work)(layout, *m, *n, a_r, lda_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SGETRF", ret);
}

/*****************************************************************************
 * SGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 *****************************************************************************/
//#define SGETRS_TEST LAPACK_GLOBAL_SUFFIX(sgetrs_test, SGETRS_TEST)
void SGETRS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
                 const float *a, const lapack_int *lda, const lapack_int *ipiv,
                 float *b, const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const float *a_r = a;
    float *b_r = b;
    lapack_int lda_r = *lda;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SGETRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgetrs)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                     ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_sgetrs_work)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                          ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SGETRS", ret);
}

/*****************************************************************************
 * SGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
 *****************************************************************************/
//#define SGETRI_TEST LAPACK_GLOBAL_SUFFIX(sgetri_test, SGETRI_TEST)
void SGETRI_TEST(const lapack_int *n, float *a, const lapack_int *lda,
                 const lapack_int *ipiv, float *work, const lapack_int *lwork,
                 lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    if (*lwork == -1) {
        ret = API_SUFFIX(LAPACKE_sgetri_work)(LAPACK_COL_MAJOR, *n, a, *lda,
                                              ipiv, work, *lwork);
        *info = lapacke_test_info("SGETRI", ret);
        return;
    }
#endif

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SGETRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgetri)(layout, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_sgetri_work)(layout, *n, a_r, lda_r, ipiv, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SGETRI", ret);
}

/*****************************************************************************
 * SPOTRF( UPLO, N, A, LDA, INFO )
 *****************************************************************************/
//#define SPOTRF_TEST LAPACK_GLOBAL_SUFFIX(spotrf_test, SPOTRF_TEST)
void SPOTRF_TEST(const char *uplo, const lapack_int *n, float *a,
                 const lapack_int *lda, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_spo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SPOTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spotrf)(layout, *uplo, *n, a_r, lda_r);
#else
    ret = API_SUFFIX(LAPACKE_spotrf_work)(layout, *uplo, *n, a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_spo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SPOTRF", ret);
}

/*****************************************************************************
 * SPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
 *****************************************************************************/
//#define SPOTRS_TEST LAPACK_GLOBAL_SUFFIX(spotrs_test, SPOTRS_TEST)
void SPOTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const float *a, const lapack_int *lda, float *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const float *a_r = a;
    float *b_r = b;
    lapack_int lda_r = *lda;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_spo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SPOTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_spotrs)(layout, *uplo, *n, *nrhs, a_r, lda_r, b_r,
                                     ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_spotrs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SPOTRS", ret);
}

/*****************************************************************************
 * SGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 *****************************************************************************/
//#define SGEQRF_TEST LAPACK_GLOBAL_SUFFIX(sgeqrf_test, SGEQRF_TEST)
void SGEQRF_TEST(const lapack_int *m, const lapack_int *n, float *a,
                 const lapack_int *lda, float *tau, float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    if (*lwork == -1) {
        ret = API_SUFFIX(LAPACKE_sgeqrf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("SGEQRF", ret);
        return;
    }
#endif

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SGEQRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgeqrf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_sgeqrf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SGEQRF", ret);
}

/*****************************************************************************
 * SORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 *****************************************************************************/
//#define SORGQR_TEST LAPACK_GLOBAL_SUFFIX(sorgqr_test, SORGQR_TEST)
void SORGQR_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 float *a, const lapack_int *lda, const float *tau, float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    if (*lwork == -1) {
        ret = API_SUFFIX(LAPACKE_sorgqr_work)(LAPACK_COL_MAJOR, *m, *n, *k, a,
                                              *lda, tau, work, *lwork);
        *info = lapacke_test_info("SORGQR", ret);
        return;
    }
#endif

    float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("SORGQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sorgqr)(layout, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_sorgqr_work)(layout, *m, *n, *k, a_r, lda_r, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("SORGQR", ret);
}

/*****************************************************************************
 * SORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
 *****************************************************************************/
//#define SORMQR_TEST LAPACK_GLOBAL_SUFFIX(sormqr_test, SORMQR_TEST)
void SORMQR_TEST(const char *side, const char *trans, const lapack_int *m,
                 const lapack_int *n, const lapack_int *k, const float *a,
                 const lapack_int *lda, const float *tau, float *c,
                 const lapack_int *ldc, float *work, const lapack_int *lwork,
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
        ret = API_SUFFIX(LAPACKE_sormqr_work)(LAPACK_COL_MAJOR, *side, *trans,
                                              *m, *n, *k, a, *lda, tau, c, *ldc,
                                              work, *lwork);
        *info = lapacke_test_info("SORMQR", ret);
        return;
    }
#endif

    const float *a_r = a;
    float *c_r = c;
    lapack_int lda_r = *lda;
    lapack_int ldc_r = *ldc;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int r = API_SUFFIX(LAPACKE_lsame)(*side, 'l') ? *m : *n;
    a_r = lapacke_test_sge_cm_to_rm(r, *k, a, *lda, &lda_r);
    c_r = lapacke_test_sge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    if (a_r == NULL || c_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("SORMQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sormqr)(layout, *side, *trans, *m, *n, *k, a_r,
                                     lda_r, tau, c_r, ldc_r);
#else
    ret = API_SUFFIX(LAPACKE_sormqr_work)(layout, *side, *trans, *m, *n, *k,
                                          a_r, lda_r, tau, c_r, ldc_r, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("SORMQR", ret);
}

/*****************************************************************************
 * SGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
 *****************************************************************************/
//#define SGELS_TEST LAPACK_GLOBAL_SUFFIX(sgels_test, SGELS_TEST)
void SGELS_TEST(const char *trans, const lapack_int *m, const lapack_int *n,
                const lapack_int *nrhs, float *a, const lapack_int *lda,
                float *b, const lapack_int *ldb, float *work,
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
        ret = API_SUFFIX(LAPACKE_sgels_work)(LAPACK_COL_MAJOR, *trans, *m, *n,
                                             *nrhs, a, *lda, b, *ldb, work,
                                             *lwork);
        *info = lapacke_test_info("SGELS", ret);
        return;
    }
#endif

    float *a_r = a;
    float *b_r = b;
    lapack_int lda_r = *lda;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int brows = MAX(*m, *n);
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_sge_cm_to_rm(brows, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("SGELS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_sgels)(layout, *trans, *m, *n, *nrhs, a_r, lda_r,
                                    b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_sgels_work)(layout, *trans, *m, *n, *nrhs, a_r,
                                         lda_r, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_sge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_sge_rm_to_cm(brows, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("SGELS", ret);
}

/*****************************************************************************
 * REAL FUNCTION SLANGE( NORM, M, N, A, LDA, WORK )
 *****************************************************************************/
//#define SLANGE_TEST LAPACK_GLOBAL_SUFFIX(slange_test, SLANGE_TEST)
float SLANGE_TEST(const char *norm, const lapack_int *m, const lapack_int *n,
                  const float *a, const lapack_int *lda, float *work
#ifdef LAPACK_FORTRAN_STRLEN_END
                  ,
                  FORTRAN_STRLEN norm_len
#endif
)
{
    float res = 0;
    const float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_sge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("SLANGE", &info);
        return (float)info;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    res = API_SUFFIX(LAPACKE_slange)(layout, *norm, *m, *n, a_r, lda_r);
#else
    res = API_SUFFIX(LAPACKE_slange_work)(layout, *norm, *m, *n, a_r, lda_r,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    return res;
}
