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

#include "cblas.h"

void ref_saxpy(int n, float alpha, const float *x, int incx, float *y,
               int incy);
void ref_daxpy(int n, double alpha, const double *x, int incx, double *y,
               int incy);

void ref_sgemv(enum CBLAS_TRANSPOSE trans, int m, int n, float alpha,
               const float *A, int lda, const float *x, int incx, float beta,
               float *y, int incy);
void ref_dgemv(enum CBLAS_TRANSPOSE trans, int m, int n, double alpha,
               const double *A, int lda, const double *x, int incx, double beta,
               double *y, int incy);

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

void ref_ssyrk(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
               float alpha, const float *A, int lda, float beta, float *C,
               int ldc);
void ref_dsyrk(enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE trans, int n, int k,
               double alpha, const double *A, int lda, double beta, double *C,
               int ldc);

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

void make_tri_f32(float *A, int n, int lda, enum CBLAS_UPLO uplo, int unit);
void make_tri_f64(double *A, int n, int lda, enum CBLAS_UPLO uplo, int unit);

#endif /* TEST_WASM_REF_H */
