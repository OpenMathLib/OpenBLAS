/*****************************************************************************
 * Shared helpers for the LAPACKE test wrapper libraries.
 * See lapacke_test_wrappers.h.
 *****************************************************************************/

#include <stdio.h>
#include <string.h>
#include "lapacke_test_wrappers.h"
#define API_SUFFIX(x) x
void xerbla (const char*, const lapack_int*, int);
/* The testing XERBLA provided by the test programs (non-halting). */
#if 0
#define fortran_xerbla LAPACK_GLOBAL_SUFFIX(xerbla, XERBLA)
void fortran_xerbla(const char *srname, const lapack_int *info
#ifdef LAPACK_FORTRAN_STRLEN_END
                    ,
                    FORTRAN_STRLEN srname_len
#endif
);
#endif
/**
 * \brief Map a LAPACKE info return value back to Fortran numbering and report
 * argument errors through the testing XERBLA.
 *
 * LAPACKE returns negative info shifted by one because of the leading
 * matrix_layout argument. This undoes that shift and calls the (interposed,
 * non-halting) testing XERBLA so that errors detected by LAPACKE itself
 * (rather than by the underlying Fortran routine) are also reported to the
 * CHKXER machinery. Memory errors and an unexpected matrix_layout error are
 * passed through with a diagnostic on stderr.
 *
 * \param[in] srname Fortran routine name reported to XERBLA.
 * \param[in] ret    Return value of the LAPACKE call.
 * \return The info value in Fortran numbering.
 */
lapack_int lapacke_test_info(const char *srname, lapack_int ret)
{
    lapack_int info = ret;
    lapack_int pos;
    if (ret == LAPACK_WORK_MEMORY_ERROR ||
        ret == LAPACK_TRANSPOSE_MEMORY_ERROR) {
        fprintf(stderr,
                "LAPACKE test wrapper %s: memory allocation error"
                " (info = %d)\n",
                srname, (int)ret);
        return ret;
    }
    if (ret < 0) {
        if (ret == -1) {
            /* Invalid matrix_layout cannot originate from these wrappers */
            fprintf(stderr,
                    "LAPACKE test wrapper %s: unexpected"
                    " matrix_layout error\n",
                    srname);
            return ret;
        }
        info = ret + 1; /* undo the matrix_layout shift */
        pos = -info;
        //fortran_xerbla(srname, &pos
        xerbla(srname, &pos
#ifdef LAPACK_FORTRAN_STRLEN_END
                       ,
                       (FORTRAN_STRLEN)strlen(srname)
#endif
        );
    }
    return info;
}

#if LAPACKE_TEST_LAYOUT == LAPACK_ROW_MAJOR

/**
 * \brief Allocate a row-major shadow copy of an m-by-n column-major general
 * matrix.
 *
 * \param[in]  m   Number of rows.
 * \param[in]  n   Number of columns.
 * \param[in]  a   Column-major source matrix.
 * \param[in]  lda Leading dimension of a.
 * \param[out] ldr Leading dimension of the shadow copy.
 * \return The shadow copy (release with LAPACKE_free), or NULL if the
 *         allocation failed.
 */
double *lapacke_test_dge_cm_to_rm(lapack_int m, lapack_int n, const double *a,
                                  lapack_int lda, lapack_int *ldr)
{
    double *r;
    *ldr = MAX(1, n);
    r = (double *)LAPACKE_malloc(sizeof(double) * MAX(1, m) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_dge_trans)(LAPACK_COL_MAJOR, m, n, a, lda, r, *ldr);
    }
    return r;
}

/**
 * \brief Copy a row-major shadow buffer back into a column-major general
 * matrix.
 *
 * \param[in]  m   Number of rows.
 * \param[in]  n   Number of columns.
 * \param[in]  r   Row-major shadow buffer.
 * \param[in]  ldr Leading dimension of r.
 * \param[out] a   Column-major destination matrix.
 * \param[in]  lda Leading dimension of a.
 */
void lapacke_test_dge_rm_to_cm(lapack_int m, lapack_int n, const double *r,
                               lapack_int ldr, double *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_dge_trans)(LAPACK_ROW_MAJOR, m, n, r, ldr, a, lda);
}

/**
 * \brief Allocate a row-major shadow copy of the uplo triangle of a
 * column-major symmetric positive definite matrix.
 *
 * Mirrors LAPACKE's own dpo_trans usage: only the uplo triangle is copied.
 *
 * \param[in]  uplo 'U' or 'L' triangle to copy.
 * \param[in]  n    Matrix order.
 * \param[in]  a    Column-major source matrix.
 * \param[in]  lda  Leading dimension of a.
 * \param[out] ldr  Leading dimension of the shadow copy.
 * \return The shadow copy (release with LAPACKE_free), or NULL if the
 *         allocation failed.
 */
double *lapacke_test_dpo_cm_to_rm(char uplo, lapack_int n, const double *a,
                                  lapack_int lda, lapack_int *ldr)
{
    double *r;
    *ldr = MAX(1, n);
    r = (double *)LAPACKE_malloc(sizeof(double) * MAX(1, n) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_dpo_trans)(LAPACK_COL_MAJOR, uplo, n, a, lda, r,
                                      *ldr);
    }
    return r;
}

/**
 * \brief Copy the uplo triangle of a row-major shadow buffer back into a
 * column-major symmetric positive definite matrix.
 *
 * \param[in]  uplo 'U' or 'L' triangle to copy.
 * \param[in]  n    Matrix order.
 * \param[in]  r    Row-major shadow buffer.
 * \param[in]  ldr  Leading dimension of r.
 * \param[out] a    Column-major destination matrix.
 * \param[in]  lda  Leading dimension of a.
 */
void lapacke_test_dpo_rm_to_cm(char uplo, lapack_int n, const double *r,
                               lapack_int ldr, double *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_dpo_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ldr, a, lda);
}

/**
 * \brief Allocate a row-major shadow copy of the uplo triangle of a
 * column-major symmetric matrix.
 *
 * Mirrors LAPACKE's own dsy_trans usage: only the uplo triangle is copied.
 *
 * \param[in]  uplo 'U' or 'L' triangle to copy.
 * \param[in]  n    Matrix order.
 * \param[in]  a    Column-major source matrix.
 * \param[in]  lda  Leading dimension of a.
 * \param[out] ldr  Leading dimension of the shadow copy.
 * \return The shadow copy (release with LAPACKE_free), or NULL if the
 *         allocation failed.
 */
double *lapacke_test_dsy_cm_to_rm(char uplo, lapack_int n, const double *a,
                                  lapack_int lda, lapack_int *ldr)
{
    double *r;
    *ldr = MAX(1, n);
    r = (double *)LAPACKE_malloc(sizeof(double) * MAX(1, n) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_dsy_trans)(LAPACK_COL_MAJOR, uplo, n, a, lda, r,
                                      *ldr);
    }
    return r;
}

/**
 * \brief Copy the uplo triangle of a row-major shadow buffer back into a
 * column-major symmetric matrix.
 *
 * \param[in]  uplo 'U' or 'L' triangle to copy.
 * \param[in]  n    Matrix order.
 * \param[in]  r    Row-major shadow buffer.
 * \param[in]  ldr  Leading dimension of r.
 * \param[out] a    Column-major destination matrix.
 * \param[in]  lda  Leading dimension of a.
 */
void lapacke_test_dsy_rm_to_cm(char uplo, lapack_int n, const double *r,
                               lapack_int ldr, double *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_dsy_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ldr, a, lda);
}

/**
 * \brief Allocate a row-major shadow copy of the referenced triangle of a
 * column-major triangular matrix.
 *
 * Mirrors LAPACKE's own dtr_trans usage: only the uplo triangle is copied,
 * without the diagonal for unit-diagonal (diag = 'U') matrices.
 *
 * \param[in]  uplo 'U' or 'L' triangle to copy.
 * \param[in]  diag 'N' for a full triangle, 'U' to skip the (unit)
 *                  diagonal.
 * \param[in]  n    Matrix order.
 * \param[in]  a    Column-major source matrix.
 * \param[in]  lda  Leading dimension of a.
 * \param[out] ldr  Leading dimension of the shadow copy.
 * \return The shadow copy (release with LAPACKE_free), or NULL if the
 *         allocation failed.
 */
double *lapacke_test_dtr_cm_to_rm(char uplo, char diag, lapack_int n,
                                  const double *a, lapack_int lda,
                                  lapack_int *ldr)
{
    double *r;
    *ldr = MAX(1, n);
    r = (double *)LAPACKE_malloc(sizeof(double) * MAX(1, n) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_dtr_trans)(LAPACK_COL_MAJOR, uplo, diag, n, a, lda,
                                      r, *ldr);
    }
    return r;
}

/**
 * \brief Copy the referenced triangle of a row-major shadow buffer back into a
 * column-major triangular matrix.
 *
 * \param[in]  uplo 'U' or 'L' triangle to copy.
 * \param[in]  diag 'N' for a full triangle, 'U' to skip the (unit)
 *                  diagonal.
 * \param[in]  n    Matrix order.
 * \param[in]  r    Row-major shadow buffer.
 * \param[in]  ldr  Leading dimension of r.
 * \param[out] a    Column-major destination matrix.
 * \param[in]  lda  Leading dimension of a.
 */
void lapacke_test_dtr_rm_to_cm(char uplo, char diag, lapack_int n,
                               const double *r, lapack_int ldr, double *a,
                               lapack_int lda)
{
    API_SUFFIX(LAPACKE_dtr_trans)(LAPACK_ROW_MAJOR, uplo, diag, n, r, ldr, a,
                                  lda);
}

/**
 * \brief Allocate a row-major shadow copy of an m-by-n column-major general
 * matrix.
 *
 * \param[in]  m   Number of rows.
 * \param[in]  n   Number of columns.
 * \param[in]  a   Column-major source matrix.
 * \param[in]  lda Leading dimension of a.
 * \param[out] ldr Leading dimension of the shadow copy.
 * \return The shadow copy (release with LAPACKE_free), or NULL if the
 *         allocation failed.
 */
float *lapacke_test_sge_cm_to_rm(lapack_int m, lapack_int n, const float *a,
                                 lapack_int lda, lapack_int *ldr)
{
    float *r;
    *ldr = MAX(1, n);
    r = (float *)LAPACKE_malloc(sizeof(float) * MAX(1, m) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_sge_trans)(LAPACK_COL_MAJOR, m, n, a, lda, r, *ldr);
    }
    return r;
}

/**
 * \brief Copy a row-major shadow buffer back into a column-major general
 * matrix.
 *
 * \param[in]  m   Number of rows.
 * \param[in]  n   Number of columns.
 * \param[in]  r   Row-major shadow buffer.
 * \param[in]  ldr Leading dimension of r.
 * \param[out] a   Column-major destination matrix.
 * \param[in]  lda Leading dimension of a.
 */
void lapacke_test_sge_rm_to_cm(lapack_int m, lapack_int n, const float *r,
                               lapack_int ldr, float *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_sge_trans)(LAPACK_ROW_MAJOR, m, n, r, ldr, a, lda);
}

/**
 * \brief Allocate a row-major shadow copy of the uplo triangle of a
 * column-major symmetric positive definite matrix.
 *
 * Mirrors LAPACKE's own dpo_trans usage: only the uplo triangle is copied.
 *
 * \param[in]  uplo 'U' or 'L' triangle to copy.
 * \param[in]  n    Matrix order.
 * \param[in]  a    Column-major source matrix.
 * \param[in]  lda  Leading dimension of a.
 * \param[out] ldr  Leading dimension of the shadow copy.
 * \return The shadow copy (release with LAPACKE_free), or NULL if the
 *         allocation failed.
 */
float *lapacke_test_spo_cm_to_rm(char uplo, lapack_int n, const float *a,
                                 lapack_int lda, lapack_int *ldr)
{
    float *r;
    *ldr = MAX(1, n);
    r = (float *)LAPACKE_malloc(sizeof(float) * MAX(1, n) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_spo_trans)(LAPACK_COL_MAJOR, uplo, n, a, lda, r,
                                      *ldr);
    }
    return r;
}

/**
 * \brief Copy the uplo triangle of a row-major shadow buffer back into a
 * column-major symmetric positive definite matrix.
 *
 * \param[in]  uplo 'U' or 'L' triangle to copy.
 * \param[in]  n    Matrix order.
 * \param[in]  r    Row-major shadow buffer.
 * \param[in]  ldr  Leading dimension of r.
 * \param[out] a    Column-major destination matrix.
 * \param[in]  lda  Leading dimension of a.
 */
void lapacke_test_spo_rm_to_cm(char uplo, lapack_int n, const float *r,
                               lapack_int ldr, float *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_spo_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ldr, a, lda);
}

/**
 * \brief Allocate a row-major shadow copy of an m-by-n column-major general
 * matrix.
 *
 * \param[in]  m   Number of rows.
 * \param[in]  n   Number of columns.
 * \param[in]  a   Column-major source matrix.
 * \param[in]  lda Leading dimension of a.
 * \param[out] ldr Leading dimension of the shadow copy.
 * \return The shadow copy (release with LAPACKE_free), or NULL if the
 *         allocation failed.
 */
lapack_complex_float *lapacke_test_cge_cm_to_rm(lapack_int m, lapack_int n,
                                                const lapack_complex_float *a,
                                                lapack_int lda, lapack_int *ldr)
{
    lapack_complex_float *r;
    *ldr = MAX(1, n);
    r = (lapack_complex_float *)LAPACKE_malloc(sizeof(lapack_complex_float) *
                                               MAX(1, m) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_cge_trans)(LAPACK_COL_MAJOR, m, n, a, lda, r, *ldr);
    }
    return r;
}

/**
 * \brief Copy a row-major shadow buffer back into a column-major general
 * matrix.
 *
 * \param[in]  m   Number of rows.
 * \param[in]  n   Number of columns.
 * \param[in]  r   Row-major shadow buffer.
 * \param[in]  ldr Leading dimension of r.
 * \param[out] a   Column-major destination matrix.
 * \param[in]  lda Leading dimension of a.
 */
void lapacke_test_cge_rm_to_cm(lapack_int m, lapack_int n,
                               const lapack_complex_float *r, lapack_int ldr,
                               lapack_complex_float *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_cge_trans)(LAPACK_ROW_MAJOR, m, n, r, ldr, a, lda);
}

/**
 * \brief Allocate a row-major shadow copy of the uplo triangle of a
 * column-major symmetric positive definite matrix.
 *
 * Mirrors LAPACKE's own dpo_trans usage: only the uplo triangle is copied.
 *
 * \param[in]  uplo 'U' or 'L' triangle to copy.
 * \param[in]  n    Matrix order.
 * \param[in]  a    Column-major source matrix.
 * \param[in]  lda  Leading dimension of a.
 * \param[out] ldr  Leading dimension of the shadow copy.
 * \return The shadow copy (release with LAPACKE_free), or NULL if the
 *         allocation failed.
 */
lapack_complex_float *lapacke_test_cpo_cm_to_rm(char uplo, lapack_int n,
                                                const lapack_complex_float *a,
                                                lapack_int lda, lapack_int *ldr)
{
    lapack_complex_float *r;
    *ldr = MAX(1, n);
    r = (lapack_complex_float *)LAPACKE_malloc(sizeof(lapack_complex_float) *
                                               MAX(1, n) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_cpo_trans)(LAPACK_COL_MAJOR, uplo, n, a, lda, r,
                                      *ldr);
    }
    return r;
}

/**
 * \brief Copy the uplo triangle of a row-major shadow buffer back into a
 * column-major symmetric positive definite matrix.
 *
 * \param[in]  uplo 'U' or 'L' triangle to copy.
 * \param[in]  n    Matrix order.
 * \param[in]  r    Row-major shadow buffer.
 * \param[in]  ldr  Leading dimension of r.
 * \param[out] a    Column-major destination matrix.
 * \param[in]  lda  Leading dimension of a.
 */
void lapacke_test_cpo_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_float *r, lapack_int ldr,
                               lapack_complex_float *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_cpo_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ldr, a, lda);
}

/**
 * \brief Allocate a row-major shadow copy of an m-by-n column-major general
 * matrix.
 *
 * \param[in]  m   Number of rows.
 * \param[in]  n   Number of columns.
 * \param[in]  a   Column-major source matrix.
 * \param[in]  lda Leading dimension of a.
 * \param[out] ldr Leading dimension of the shadow copy.
 * \return The shadow copy (release with LAPACKE_free), or NULL if the
 *         allocation failed.
 */
lapack_complex_double *lapacke_test_zge_cm_to_rm(lapack_int m, lapack_int n,
                                                 const lapack_complex_double *a,
                                                 lapack_int lda,
                                                 lapack_int *ldr)
{
    lapack_complex_double *r;
    *ldr = MAX(1, n);
    r = (lapack_complex_double *)LAPACKE_malloc(sizeof(lapack_complex_double) *
                                                MAX(1, m) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_zge_trans)(LAPACK_COL_MAJOR, m, n, a, lda, r, *ldr);
    }
    return r;
}

/**
 * \brief Copy a row-major shadow buffer back into a column-major general
 * matrix.
 *
 * \param[in]  m   Number of rows.
 * \param[in]  n   Number of columns.
 * \param[in]  r   Row-major shadow buffer.
 * \param[in]  ldr Leading dimension of r.
 * \param[out] a   Column-major destination matrix.
 * \param[in]  lda Leading dimension of a.
 */
void lapacke_test_zge_rm_to_cm(lapack_int m, lapack_int n,
                               const lapack_complex_double *r, lapack_int ldr,
                               lapack_complex_double *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_zge_trans)(LAPACK_ROW_MAJOR, m, n, r, ldr, a, lda);
}

/**
 * \brief Allocate a row-major shadow copy of the uplo triangle of a
 * column-major symmetric positive definite matrix.
 *
 * Mirrors LAPACKE's own dpo_trans usage: only the uplo triangle is copied.
 *
 * \param[in]  uplo 'U' or 'L' triangle to copy.
 * \param[in]  n    Matrix order.
 * \param[in]  a    Column-major source matrix.
 * \param[in]  lda  Leading dimension of a.
 * \param[out] ldr  Leading dimension of the shadow copy.
 * \return The shadow copy (release with LAPACKE_free), or NULL if the
 *         allocation failed.
 */
lapack_complex_double *lapacke_test_zpo_cm_to_rm(char uplo, lapack_int n,
                                                 const lapack_complex_double *a,
                                                 lapack_int lda,
                                                 lapack_int *ldr)
{
    lapack_complex_double *r;
    *ldr = MAX(1, n);
    r = (lapack_complex_double *)LAPACKE_malloc(sizeof(lapack_complex_double) *
                                                MAX(1, n) * (*ldr));
    if (r != NULL) {
        API_SUFFIX(LAPACKE_zpo_trans)(LAPACK_COL_MAJOR, uplo, n, a, lda, r,
                                      *ldr);
    }
    return r;
}

/**
 * \brief Copy the uplo triangle of a row-major shadow buffer back into a
 * column-major symmetric positive definite matrix.
 *
 * \param[in]  uplo 'U' or 'L' triangle to copy.
 * \param[in]  n    Matrix order.
 * \param[in]  r    Row-major shadow buffer.
 * \param[in]  ldr  Leading dimension of r.
 * \param[out] a    Column-major destination matrix.
 * \param[in]  lda  Leading dimension of a.
 */
void lapacke_test_zpo_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_double *r, lapack_int ldr,
                               lapack_complex_double *a, lapack_int lda)
{
    API_SUFFIX(LAPACKE_zpo_trans)(LAPACK_ROW_MAJOR, uplo, n, r, ldr, a, lda);
}

/**
 * \brief Report a failed shadow buffer allocation and set info to
 * LAPACK_TRANSPOSE_MEMORY_ERROR.
 *
 * Prints a diagnostic on stderr and sets info to
 * LAPACK_TRANSPOSE_MEMORY_ERROR, matching what LAPACKE itself returns when
 * its transposition buffers cannot be allocated.
 *
 * \param[in]  srname Routine name for the diagnostic.
 * \param[out] info   The wrapper's info result.
 */
void lapacke_test_report_alloc_failure(const char *srname, lapack_int *info)
{
    *info = LAPACK_TRANSPOSE_MEMORY_ERROR;
    fprintf(stderr,
            "LAPACKE test wrapper %s: shadow buffer allocation"
            " failed\n",
            srname);
}

#endif /* LAPACK_ROW_MAJOR */
