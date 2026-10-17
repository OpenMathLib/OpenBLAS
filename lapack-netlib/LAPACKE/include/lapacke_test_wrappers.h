/*****************************************************************************
 * Shared helpers for the LAPACKE test wrapper libraries.
 *
 * The wrappers (one source file per precision, lapacke_test_wrappers_<x>.c)
 * export Fortran-callable symbols <NAME>_TEST with the exact argument list
 * of the corresponding LAPACK routine. The LAPACK testing programs are
 * rewritten by the extended-API source generator so that calls to
 * allowlisted routines resolve to these wrappers, which forward to LAPACKE.
 * This way the unmodified numerical and error-exit checks of the Fortran
 * test suite exercise LAPACKE.
 *
 * The wrapper libraries are compiled once per (layout, layer) combination:
 *
 *   LAPACKE_TEST_LAYOUT=LAPACK_COL_MAJOR -- faithful passthrough
 *   LAPACKE_TEST_LAYOUT=LAPACK_ROW_MAJOR -- caller data is transposed into
 *     row-major shadow buffers so that LAPACKE's row-major conversion
 *     machinery is exercised; results are transposed back afterwards.
 *
 *   LAPACKE_TEST_LAYER=LAPACKE_TEST_LAYER_WORK -- call the middle-level
 *     (_work) interface; the caller's WORK/LWORK arguments pass through,
 *     which keeps workspace queries, minimal-workspace tests and LWORK
 *     error exits faithful to the Fortran routine.
 *   LAPACKE_TEST_LAYER=LAPACKE_TEST_LAYER_HIGH -- call the high-level
 *     interface, which performs NaN checks and allocates its own workspace.
 *     The high-level interface has no workspace-query notion, so LWORK=-1
 *     calls are serviced through the middle level (in column-major: only
 *     the dimensions determine the result, and the caller's leading
 *     dimensions are column-major-valid). The caller's WORK/LWORK are
 *     otherwise ignored, so LWORK error exits cannot fire; suite runs for
 *     this layer must disable the error-exit tests.
 *
 * The helpers declared here are implemented in lapacke_test_wrappers.c,
 * which is compiled into each wrapper library.
 *****************************************************************************/

#ifndef LAPACKE_TEST_WRAPPERS_H
#define LAPACKE_TEST_WRAPPERS_H

#include "lapacke.h"
#include "lapacke_utils.h"

#ifndef LAPACKE_TEST_LAYOUT
#error "LAPACKE_TEST_LAYOUT must be defined as LAPACK_{COL,ROW}_MAJOR"
#endif

/** LAPACKE_TEST_LAYER value selecting the middle-level (_work) API. */
#define LAPACKE_TEST_LAYER_WORK 1
/** LAPACKE_TEST_LAYER value selecting the high-level API. */
#define LAPACKE_TEST_LAYER_HIGH 2

#ifndef LAPACKE_TEST_LAYER
#error "LAPACKE_TEST_LAYER must be defined as LAPACKE_TEST_LAYER_{WORK,HIGH}"
#endif

/** Nonzero when the wrappers should call the high-level interface. */
#define LAPACKE_TEST_HIGH_LEVEL (LAPACKE_TEST_LAYER == LAPACKE_TEST_LAYER_HIGH)
/** Nonzero when the wrappers should use row-major layout. */
#define LAPACKE_TEST_ROW_MAJOR (LAPACKE_TEST_LAYOUT == LAPACK_ROW_MAJOR)

/** Map a LAPACKE info return value back to Fortran numbering and report
 *  argument errors through the testing XERBLA. */
lapack_int lapacke_test_info(const char *srname, lapack_int ret);

#if LAPACKE_TEST_LAYOUT == LAPACK_ROW_MAJOR

/* One cm_to_rm/rm_to_cm pair per conventional-storage matrix type,
 * mirroring LAPACKE's own <x><type>_trans usage. Band and packed storage
 * have their own buffer geometry and get their helpers with the wrappers
 * that first need them. */

/** Allocate a row-major shadow copy of an m-by-n column-major general
 *  matrix. */
double *lapacke_test_dge_cm_to_rm(lapack_int m, lapack_int n, const double *a,
                                  lapack_int lda, lapack_int *ldr);

/** Copy a row-major shadow buffer back into a column-major general
 *  matrix. */
void lapacke_test_dge_rm_to_cm(lapack_int m, lapack_int n, const double *r,
                               lapack_int ldr, double *a, lapack_int lda);

/** Allocate a row-major shadow copy of the uplo triangle of a column-major
 *  symmetric positive definite matrix. */
double *lapacke_test_dpo_cm_to_rm(char uplo, lapack_int n, const double *a,
                                  lapack_int lda, lapack_int *ldr);

/** Copy the uplo triangle of a row-major shadow buffer back into a
 *  column-major symmetric positive definite matrix. */
void lapacke_test_dpo_rm_to_cm(char uplo, lapack_int n, const double *r,
                               lapack_int ldr, double *a, lapack_int lda);

/** Allocate a row-major shadow copy of the uplo triangle of a column-major
 *  symmetric matrix. */
double *lapacke_test_dsy_cm_to_rm(char uplo, lapack_int n, const double *a,
                                  lapack_int lda, lapack_int *ldr);

/** Copy the uplo triangle of a row-major shadow buffer back into a
 *  column-major symmetric matrix. */
void lapacke_test_dsy_rm_to_cm(char uplo, lapack_int n, const double *r,
                               lapack_int ldr, double *a, lapack_int lda);

/** Allocate a row-major shadow copy of the referenced triangle of a
 *  column-major triangular matrix. */
double *lapacke_test_dtr_cm_to_rm(char uplo, char diag, lapack_int n,
                                  const double *a, lapack_int lda,
                                  lapack_int *ldr);

/** Copy the referenced triangle of a row-major shadow buffer back into a
 *  column-major triangular matrix. */
void lapacke_test_dtr_rm_to_cm(char uplo, char diag, lapack_int n,
                               const double *r, lapack_int ldr, double *a,
                               lapack_int lda);

/** Allocate a row-major shadow copy of an m-by-n column-major general
 *  matrix. */
float *lapacke_test_sge_cm_to_rm(lapack_int m, lapack_int n, const float *a,
                                 lapack_int lda, lapack_int *ldr);

/** Copy a row-major shadow buffer back into a column-major general
 *  matrix. */
void lapacke_test_sge_rm_to_cm(lapack_int m, lapack_int n, const float *r,
                               lapack_int ldr, float *a, lapack_int lda);

/** Allocate a row-major shadow copy of the uplo triangle of a column-major
 *  symmetric positive definite matrix. */
float *lapacke_test_spo_cm_to_rm(char uplo, lapack_int n, const float *a,
                                 lapack_int lda, lapack_int *ldr);

/** Copy the uplo triangle of a row-major shadow buffer back into a
 *  column-major symmetric positive definite matrix. */
void lapacke_test_spo_rm_to_cm(char uplo, lapack_int n, const float *r,
                               lapack_int ldr, float *a, lapack_int lda);

/** Allocate a row-major shadow copy of an m-by-n column-major general
 *  matrix. */
lapack_complex_float *lapacke_test_cge_cm_to_rm(lapack_int m, lapack_int n,
                                                const lapack_complex_float *a,
                                                lapack_int lda,
                                                lapack_int *ldr);

/** Copy a row-major shadow buffer back into a column-major general
 *  matrix. */
void lapacke_test_cge_rm_to_cm(lapack_int m, lapack_int n,
                               const lapack_complex_float *r, lapack_int ldr,
                               lapack_complex_float *a, lapack_int lda);

/** Allocate a row-major shadow copy of the uplo triangle of a column-major
 *  symmetric positive definite matrix. */
lapack_complex_float *lapacke_test_cpo_cm_to_rm(char uplo, lapack_int n,
                                                const lapack_complex_float *a,
                                                lapack_int lda,
                                                lapack_int *ldr);

/** Copy the uplo triangle of a row-major shadow buffer back into a
 *  column-major symmetric positive definite matrix. */
void lapacke_test_cpo_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_float *r, lapack_int ldr,
                               lapack_complex_float *a, lapack_int lda);

/** Allocate a row-major shadow copy of an m-by-n column-major general
 *  matrix. */
lapack_complex_double *lapacke_test_zge_cm_to_rm(lapack_int m, lapack_int n,
                                                 const lapack_complex_double *a,
                                                 lapack_int lda,
                                                 lapack_int *ldr);

/** Copy a row-major shadow buffer back into a column-major general
 *  matrix. */
void lapacke_test_zge_rm_to_cm(lapack_int m, lapack_int n,
                               const lapack_complex_double *r, lapack_int ldr,
                               lapack_complex_double *a, lapack_int lda);

/** Allocate a row-major shadow copy of the uplo triangle of a column-major
 *  symmetric positive definite matrix. */
lapack_complex_double *lapacke_test_zpo_cm_to_rm(char uplo, lapack_int n,
                                                 const lapack_complex_double *a,
                                                 lapack_int lda,
                                                 lapack_int *ldr);

/** Copy the uplo triangle of a row-major shadow buffer back into a
 *  column-major symmetric positive definite matrix. */
void lapacke_test_zpo_rm_to_cm(char uplo, lapack_int n,
                               const lapack_complex_double *r, lapack_int ldr,
                               lapack_complex_double *a, lapack_int lda);

/** Report a failed shadow buffer allocation and set info to
 *  LAPACK_TRANSPOSE_MEMORY_ERROR. */
void lapacke_test_report_alloc_failure(const char *srname, lapack_int *info);

#endif /* LAPACK_ROW_MAJOR */
#endif /* LAPACKE_TEST_WRAPPERS_H */
