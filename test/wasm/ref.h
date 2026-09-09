/*
Copyright (c) 2026, The OpenBLAS Project
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
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef TEST_WASM_REF_H
#define TEST_WASM_REF_H

#include <stddef.h>

#include "cblas.h"

/* Level 1 typed (ref_l1.c) — deep AXPY grid. */
void ref_saxpy(int n, float alpha, const float *x, int incx, float *y,
               int incy);
void ref_daxpy(int n, double alpha, const double *x, int incx, double *y,
               int incy);

/* Level 2 typed (ref_l2.c) — deep GEMV grid. */
void ref_sgemv(enum CBLAS_TRANSPOSE trans, int m, int n, float alpha,
               const float *A, int lda, const float *x, int incx, float beta,
               float *y, int incy);
void ref_dgemv(enum CBLAS_TRANSPOSE trans, int m, int n, double alpha,
               const double *A, int lda, const double *x, int incx, double beta,
               double *y, int incy);

/* Level 3 (ref_l3.c). Symmetric, Hermitian and triangular routines read only
 * the triangle named by `uplo`; rank-k / rank-2k write only that triangle. */

void ref_sgemm(enum CBLAS_TRANSPOSE ta, enum CBLAS_TRANSPOSE tb, int m, int n,
               int k, float alpha, const float *A, int lda, const float *B,
               int ldb, float beta, float *C, int ldc);
void ref_dgemm(enum CBLAS_TRANSPOSE ta, enum CBLAS_TRANSPOSE tb, int m, int n,
               int k, double alpha, const double *A, int lda, const double *B,
               int ldb, double beta, double *C, int ldc);
void ref_cgemm(enum CBLAS_TRANSPOSE ta, enum CBLAS_TRANSPOSE tb, int m, int n,
               int k, const float *alpha, const float *A, int lda,
               const float *B, int ldb, const float *beta, float *C, int ldc);
void ref_zgemm(enum CBLAS_TRANSPOSE ta, enum CBLAS_TRANSPOSE tb, int m, int n,
               int k, const double *alpha, const double *A, int lda,
               const double *B, int ldb, const double *beta, double *C,
               int ldc);

void ref_ssymm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo, int m, int n,
               float alpha, const float *A, int lda, const float *B, int ldb,
               float beta, float *C, int ldc);
void ref_dsymm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo, int m, int n,
               double alpha, const double *A, int lda, const double *B,
               int ldb, double beta, double *C, int ldc);
void ref_csymm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo, int m, int n,
               const float *alpha, const float *A, int lda, const float *B,
               int ldb, const float *beta, float *C, int ldc);
void ref_zsymm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo, int m, int n,
               const double *alpha, const double *A, int lda, const double *B,
               int ldb, const double *beta, double *C, int ldc);

void ref_chemm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo, int m, int n,
               const float *alpha, const float *A, int lda, const float *B,
               int ldb, const float *beta, float *C, int ldc);
void ref_zhemm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo, int m, int n,
               const double *alpha, const double *A, int lda, const double *B,
               int ldb, const double *beta, double *C, int ldc);

void ref_ssyrk(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
               float alpha, const float *A, int lda, float beta, float *C,
               int ldc);
void ref_dsyrk(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
               double alpha, const double *A, int lda, double beta, double *C,
               int ldc);
void ref_csyrk(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
               const float *alpha, const float *A, int lda, const float *beta,
               float *C, int ldc);
void ref_zsyrk(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
               const double *alpha, const double *A, int lda,
               const double *beta, double *C, int ldc);

void ref_cherk(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
               float alpha, const float *A, int lda, float beta, float *C,
               int ldc);
void ref_zherk(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
               double alpha, const double *A, int lda, double beta, double *C,
               int ldc);

void ref_ssyr2k(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
                float alpha, const float *A, int lda, const float *B, int ldb,
                float beta, float *C, int ldc);
void ref_dsyr2k(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
                double alpha, const double *A, int lda, const double *B,
                int ldb, double beta, double *C, int ldc);
void ref_csyr2k(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
                const float *alpha, const float *A, int lda, const float *B,
                int ldb, const float *beta, float *C, int ldc);
void ref_zsyr2k(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
                const double *alpha, const double *A, int lda,
                const double *B, int ldb, const double *beta, double *C,
                int ldc);

void ref_cher2k(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
                const float *alpha, const float *A, int lda, const float *B,
                int ldb, float beta, float *C, int ldc);
void ref_zher2k(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
                const double *alpha, const double *A, int lda,
                const double *B, int ldb, double beta, double *C, int ldc);

void ref_strmm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
               enum CBLAS_TRANSPOSE trans, enum CBLAS_DIAG diag, int m, int n,
               float alpha, const float *A, int lda, float *B, int ldb);
void ref_dtrmm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
               enum CBLAS_TRANSPOSE trans, enum CBLAS_DIAG diag, int m, int n,
               double alpha, const double *A, int lda, double *B, int ldb);

void ref_strsm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
               enum CBLAS_TRANSPOSE trans, enum CBLAS_DIAG diag, int m, int n,
               float alpha, const float *A, int lda, float *B, int ldb);
void ref_dtrsm(enum CBLAS_SIDE side, enum CBLAS_UPLO uplo,
               enum CBLAS_TRANSPOSE trans, enum CBLAS_DIAG diag, int m, int n,
               double alpha, const double *A, int lda, double *B, int ldb);

void ref_l1_swap(int n, void *x, int incx, void *y, int incy, int width,
                 int is_double);
void ref_l1_copy(int n, const void *x, int incx, void *y, int incy, int width,
                 int is_double);
void ref_l1_scal(int n, const void *alpha, void *x, int incx, int width,
                 int is_double, int real_alpha);
void ref_l1_axpy(int n, const void *alpha, const void *x, int incx, void *y,
                 int incy, int width, int is_double);
void ref_l1_dot(int n, const void *x, int incx, const void *y, int incy,
                void *out, int width, int is_double, int conjugate);
double ref_l1_nrm2(int n, const void *x, int incx, int width, int is_double);
double ref_l1_asum(int n, const void *x, int incx, int width, int is_double);
size_t ref_l1_iamax(int n, const void *x, int incx, int width, int is_double);
void ref_l1_rot(int n, void *x, int incx, void *y, int incy, double c,
                double s, int width, int is_double);
void ref_l1_rotm(int n, void *x, int incx, void *y, int incy,
                 const void *param, int is_double);
void ref_l1_rotmg(void *d1, void *d2, void *b1, const void *b2, void *param,
                  int is_double);

void ref_l2_mv(int n, const void *a, const void *x, void *y,
               const void *alpha, const void *beta, int width, int is_double,
               int conjugate);
void ref_l2_rank(int n, void *a, const void *x, const void *y,
                 const void *alpha, int width, int is_double, int conjugate_y,
                 int symmetric, int hermitian);

#endif /* TEST_WASM_REF_H */
