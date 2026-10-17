/*****************************************************************************
 * LAPACKE test wrappers (single precision complex, spike subset)
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
 * CGETRF( M, N, A, LDA, IPIV, INFO )
 *****************************************************************************/
//#define CGETRF_TEST LAPACK_GLOBAL_SUFFIX(cgetrf_test, CGETRF_TEST)
void CGETRF_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_float *a, const lapack_int *lda,
                 lapack_int *ipiv, lapack_int *info)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CGETRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgetrf)(layout, *m, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_cgetrf_work)(layout, *m, *n, a_r, lda_r, ipiv);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CGETRF", ret);
}

/*****************************************************************************
 * CGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
 *****************************************************************************/
//#define CGETRS_TEST LAPACK_GLOBAL_SUFFIX(cgetrs_test, CGETRS_TEST)
void CGETRS_TEST(const char *trans, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_float *a, const lapack_int *lda,
                 const lapack_int *ipiv, lapack_complex_float *b,
                 const lapack_int *ldb, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN trans_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_complex_float *b_r = b;
    lapack_int lda_r = *lda;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CGETRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgetrs)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                     ipiv, b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_cgetrs_work)(layout, *trans, *n, *nrhs, a_r, lda_r,
                                          ipiv, b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CGETRS", ret);
}

/*****************************************************************************
 * CGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
 *****************************************************************************/
//#define CGETRI_TEST LAPACK_GLOBAL_SUFFIX(cgetri_test, CGETRI_TEST)
void CGETRI_TEST(const lapack_int *n, lapack_complex_float *a,
                 const lapack_int *lda, const lapack_int *ipiv,
                 lapack_complex_float *work, const lapack_int *lwork,
                 lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    if (*lwork == -1) {
        ret = API_SUFFIX(LAPACKE_cgetri_work)(LAPACK_COL_MAJOR, *n, a, *lda,
                                              ipiv, work, *lwork);
        *info = lapacke_test_info("CGETRI", ret);
        return;
    }
#endif

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*n, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CGETRI", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgetri)(layout, *n, a_r, lda_r, ipiv);
#else
    ret = API_SUFFIX(LAPACKE_cgetri_work)(layout, *n, a_r, lda_r, ipiv, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CGETRI", ret);
}

/*****************************************************************************
 * CPOTRF( UPLO, N, A, LDA, INFO )
 *****************************************************************************/
//#define CPOTRF_TEST LAPACK_GLOBAL_SUFFIX(cpotrf_test, CPOTRF_TEST)
void CPOTRF_TEST(const char *uplo, const lapack_int *n, lapack_complex_float *a,
                 const lapack_int *lda, lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CPOTRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpotrf)(layout, *uplo, *n, a_r, lda_r);
#else
    ret = API_SUFFIX(LAPACKE_cpotrf_work)(layout, *uplo, *n, a_r, lda_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cpo_rm_to_cm(*uplo, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CPOTRF", ret);
}

/*****************************************************************************
 * CPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
 *****************************************************************************/
//#define CPOTRS_TEST LAPACK_GLOBAL_SUFFIX(cpotrs_test, CPOTRS_TEST)
void CPOTRS_TEST(const char *uplo, const lapack_int *n, const lapack_int *nrhs,
                 const lapack_complex_float *a, const lapack_int *lda,
                 lapack_complex_float *b, const lapack_int *ldb,
                 lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                 ,
                 FORTRAN_STRLEN uplo_len
#endif
)
{
    lapack_int ret = 0;
    const lapack_complex_float *a_r = a;
    lapack_complex_float *b_r = b;
    lapack_int lda_r = *lda;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cpo_cm_to_rm(*uplo, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(*n, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CPOTRS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cpotrs)(layout, *uplo, *n, *nrhs, a_r, lda_r, b_r,
                                     ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_cpotrs_work)(layout, *uplo, *n, *nrhs, a_r, lda_r,
                                          b_r, ldb_r);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*n, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CPOTRS", ret);
}

/*****************************************************************************
 * CGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 *****************************************************************************/
//#define CGEQRF_TEST LAPACK_GLOBAL_SUFFIX(cgeqrf_test, CGEQRF_TEST)
void CGEQRF_TEST(const lapack_int *m, const lapack_int *n,
                 lapack_complex_float *a, const lapack_int *lda,
                 lapack_complex_float *tau, lapack_complex_float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    if (*lwork == -1) {
        ret = API_SUFFIX(LAPACKE_cgeqrf_work)(LAPACK_COL_MAJOR, *m, *n, a, *lda,
                                              tau, work, *lwork);
        *info = lapacke_test_info("CGEQRF", ret);
        return;
    }
#endif

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CGEQRF", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgeqrf)(layout, *m, *n, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_cgeqrf_work)(layout, *m, *n, a_r, lda_r, tau, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CGEQRF", ret);
}

/*****************************************************************************
 * CUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 *****************************************************************************/
//#define CUNGQR_TEST LAPACK_GLOBAL_SUFFIX(cungqr_test, CUNGQR_TEST)
void CUNGQR_TEST(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                 lapack_complex_float *a, const lapack_int *lda,
                 const lapack_complex_float *tau, lapack_complex_float *work,
                 const lapack_int *lwork, lapack_int *info)
{
    lapack_int ret = 0;
#if LAPACKE_TEST_HIGH_LEVEL
    if (*lwork == -1) {
        ret = API_SUFFIX(LAPACKE_cungqr_work)(LAPACK_COL_MAJOR, *m, *n, *k, a,
                                              *lda, tau, work, *lwork);
        *info = lapacke_test_info("CUNGQR", ret);
        return;
    }
#endif

    lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapacke_test_report_alloc_failure("CUNGQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cungqr)(layout, *m, *n, *k, a_r, lda_r, tau);
#else
    ret = API_SUFFIX(LAPACKE_cungqr_work)(layout, *m, *n, *k, a_r, lda_r, tau,
                                          work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    LAPACKE_free(a_r);
#endif
    *info = lapacke_test_info("CUNGQR", ret);
}

/*****************************************************************************
 * CUNMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
 *****************************************************************************/
//#define CUNMQR_TEST LAPACK_GLOBAL_SUFFIX(cunmqr_test, CUNMQR_TEST)
void CUNMQR_TEST(const char *side, const char *trans, const lapack_int *m,
                 const lapack_int *n, const lapack_int *k,
                 const lapack_complex_float *a, const lapack_int *lda,
                 const lapack_complex_float *tau, lapack_complex_float *c,
                 const lapack_int *ldc, lapack_complex_float *work,
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
        ret = API_SUFFIX(LAPACKE_cunmqr_work)(LAPACK_COL_MAJOR, *side, *trans,
                                              *m, *n, *k, a, *lda, tau, c, *ldc,
                                              work, *lwork);
        *info = lapacke_test_info("CUNMQR", ret);
        return;
    }
#endif

    const lapack_complex_float *a_r = a;
    lapack_complex_float *c_r = c;
    lapack_int lda_r = *lda;
    lapack_int ldc_r = *ldc;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int r = API_SUFFIX(LAPACKE_lsame)(*side, 'l') ? *m : *n;
    a_r = lapacke_test_cge_cm_to_rm(r, *k, a, *lda, &lda_r);
    c_r = lapacke_test_cge_cm_to_rm(*m, *n, c, *ldc, &ldc_r);
    if (a_r == NULL || c_r == NULL) {
        LAPACKE_free((void *)a_r);
        LAPACKE_free(c_r);
        lapacke_test_report_alloc_failure("CUNMQR", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cunmqr)(layout, *side, *trans, *m, *n, *k, a_r,
                                     lda_r, tau, c_r, ldc_r);
#else
    ret = API_SUFFIX(LAPACKE_cunmqr_work)(layout, *side, *trans, *m, *n, *k,
                                          a_r, lda_r, tau, c_r, ldc_r, work,
                                          *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, c_r, ldc_r, c, *ldc);
    LAPACKE_free((void *)a_r);
    LAPACKE_free(c_r);
#endif
    *info = lapacke_test_info("CUNMQR", ret);
}

/*****************************************************************************
 * CGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
 *****************************************************************************/
//#define CGELS_TEST LAPACK_GLOBAL_SUFFIX(cgels_test, CGELS_TEST)
void CGELS_TEST(const char *trans, const lapack_int *m, const lapack_int *n,
                const lapack_int *nrhs, lapack_complex_float *a,
                const lapack_int *lda, lapack_complex_float *b,
                const lapack_int *ldb, lapack_complex_float *work,
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
        ret = API_SUFFIX(LAPACKE_cgels_work)(LAPACK_COL_MAJOR, *trans, *m, *n,
                                             *nrhs, a, *lda, b, *ldb, work,
                                             *lwork);
        *info = lapacke_test_info("CGELS", ret);
        return;
    }
#endif

    lapack_complex_float *a_r = a;
    lapack_complex_float *b_r = b;
    lapack_int lda_r = *lda;
    lapack_int ldb_r = *ldb;
#if LAPACKE_TEST_ROW_MAJOR
    const lapack_int brows = MAX(*m, *n);
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    b_r = lapacke_test_cge_cm_to_rm(brows, *nrhs, b, *ldb, &ldb_r);
    if (a_r == NULL || b_r == NULL) {
        LAPACKE_free(a_r);
        LAPACKE_free(b_r);
        lapacke_test_report_alloc_failure("CGELS", info);
        return;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    ret = API_SUFFIX(LAPACKE_cgels)(layout, *trans, *m, *n, *nrhs, a_r, lda_r,
                                    b_r, ldb_r);
#else
    ret = API_SUFFIX(LAPACKE_cgels_work)(layout, *trans, *m, *n, *nrhs, a_r,
                                         lda_r, b_r, ldb_r, work, *lwork);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    lapacke_test_cge_rm_to_cm(*m, *n, a_r, lda_r, a, *lda);
    lapacke_test_cge_rm_to_cm(brows, *nrhs, b_r, ldb_r, b, *ldb);
    LAPACKE_free(a_r);
    LAPACKE_free(b_r);
#endif
    *info = lapacke_test_info("CGELS", ret);
}

/*****************************************************************************
 * REAL FUNCTION CLANGE( NORM, M, N, A, LDA, WORK )
 *****************************************************************************/
//#define CLANGE_TEST LAPACK_GLOBAL_SUFFIX(clange_test, CLANGE_TEST)
float CLANGE_TEST(const char *norm, const lapack_int *m, const lapack_int *n,
                  const lapack_complex_float *a, const lapack_int *lda,
                  float *work
#ifdef LAPACK_FORTRAN_STRLEN_END
                  ,
                  FORTRAN_STRLEN norm_len
#endif
)
{
    float res = 0;
    const lapack_complex_float *a_r = a;
    lapack_int lda_r = *lda;
#if LAPACKE_TEST_ROW_MAJOR
    a_r = lapacke_test_cge_cm_to_rm(*m, *n, a, *lda, &lda_r);
    if (a_r == NULL) {
        lapack_int info;
        lapacke_test_report_alloc_failure("CLANGE", &info);
        return (float)info;
    }
#endif

#if LAPACKE_TEST_HIGH_LEVEL
    res = API_SUFFIX(LAPACKE_clange)(layout, *norm, *m, *n, a_r, lda_r);
#else
    res = API_SUFFIX(LAPACKE_clange_work)(layout, *norm, *m, *n, a_r, lda_r,
                                          work);
#endif

#if LAPACKE_TEST_ROW_MAJOR
    LAPACKE_free((void *)a_r);
#endif
    return res;
}
