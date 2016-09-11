#ifndef CBLAS_H
#define CBLAS_H

#include <stddef.h>
#include "common.h"

#ifdef __cplusplus
extern "C" {
	/* Assume C declarations for C++ */
#endif  /* __cplusplus */

/*Set the number of threads on runtime.*/
void openblas_set_num_threads(int num_threads);
void goto_set_num_threads(int num_threads);

/*Get the number of threads on runtime.*/
int openblas_get_num_threads(void);

/*Get the number of physical processors (cores).*/
int openblas_get_num_procs(void);

/*Get the build configure on runtime.*/
char* openblas_get_config(void);

/*Get the CPU corename on runtime.*/
char* openblas_get_corename(void);

/* Get the parallelization type which is used by OpenBLAS */
int openblas_get_parallel(void);
/* OpenBLAS is compiled for sequential use  */
#define OPENBLAS_SEQUENTIAL  0
/* OpenBLAS is compiled using normal threading model */
#define OPENBLAS_THREAD  1
/* OpenBLAS is compiled using OpenMP threading model */
#define OPENBLAS_OPENMP 2


/*
 * Since all of GotoBlas was written without const,
 * we disable it at build time.
 */
#ifndef OPENBLAS_CONST
# define OPENBLAS_CONST const
#endif


#define CBLAS_INDEX size_t

typedef enum CBLAS_ORDER     {CblasRowMajor=101, CblasColMajor=102} CBLAS_ORDER;
typedef enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113, CblasConjNoTrans=114} CBLAS_TRANSPOSE;
typedef enum CBLAS_UPLO      {CblasUpper=121, CblasLower=122} CBLAS_UPLO;
typedef enum CBLAS_DIAG      {CblasNonUnit=131, CblasUnit=132} CBLAS_DIAG;
typedef enum CBLAS_SIDE      {CblasLeft=141, CblasRight=142} CBLAS_SIDE;

float  cblas_sdsdot(blasint n, float alpha, OPENBLAS_CONST float *x, blasint incx, OPENBLAS_CONST float *y, blasint incy);
double cblas_dsdot (blasint n, OPENBLAS_CONST float *x, blasint incx, OPENBLAS_CONST float *y, blasint incy);
float  cblas_sdot(blasint n, OPENBLAS_CONST float  *x, blasint incx, OPENBLAS_CONST float  *y, blasint incy);
double cblas_ddot(blasint n, OPENBLAS_CONST double *x, blasint incx, OPENBLAS_CONST double *y, blasint incy);

openblas_complex_float  cblas_cdotu(blasint n, OPENBLAS_CONST float  *x, blasint incx, OPENBLAS_CONST float  *y, blasint incy);
openblas_complex_float  cblas_cdotc(blasint n, OPENBLAS_CONST float  *x, blasint incx, OPENBLAS_CONST float  *y, blasint incy);
openblas_complex_double cblas_zdotu(blasint n, OPENBLAS_CONST double *x, blasint incx, OPENBLAS_CONST double *y, blasint incy);
openblas_complex_double cblas_zdotc(blasint n, OPENBLAS_CONST double *x, blasint incx, OPENBLAS_CONST double *y, blasint incy);

void  cblas_cdotu_sub(blasint n, OPENBLAS_CONST float  *x, blasint incx, OPENBLAS_CONST float  *y, blasint incy, openblas_complex_float  *ret);
void  cblas_cdotc_sub(blasint n, OPENBLAS_CONST float  *x, blasint incx, OPENBLAS_CONST float  *y, blasint incy, openblas_complex_float  *ret);
void  cblas_zdotu_sub(blasint n, OPENBLAS_CONST double *x, blasint incx, OPENBLAS_CONST double *y, blasint incy, openblas_complex_double *ret);
void  cblas_zdotc_sub(blasint n, OPENBLAS_CONST double *x, blasint incx, OPENBLAS_CONST double *y, blasint incy, openblas_complex_double *ret);

float  cblas_sasum (blasint n, OPENBLAS_CONST float  *x, blasint incx);
double cblas_dasum (blasint n, OPENBLAS_CONST double *x, blasint incx);
float  cblas_scasum(blasint n, OPENBLAS_CONST float  *x, blasint incx);
double cblas_dzasum(blasint n, OPENBLAS_CONST double *x, blasint incx);

float  cblas_snrm2 (blasint N, OPENBLAS_CONST float  *X, blasint incX);
double cblas_dnrm2 (blasint N, OPENBLAS_CONST double *X, blasint incX);
float  cblas_scnrm2(blasint N, OPENBLAS_CONST float  *X, blasint incX);
double cblas_dznrm2(blasint N, OPENBLAS_CONST double *X, blasint incX);

CBLAS_INDEX cblas_isamax(blasint n, OPENBLAS_CONST float  *x, blasint incx);
CBLAS_INDEX cblas_idamax(blasint n, OPENBLAS_CONST double *x, blasint incx);
CBLAS_INDEX cblas_icamax(blasint n, OPENBLAS_CONST float  *x, blasint incx);
CBLAS_INDEX cblas_izamax(blasint n, OPENBLAS_CONST double *x, blasint incx);

void cblas_saxpy(blasint n, float alpha, OPENBLAS_CONST float *x, blasint incx, float *y, blasint incy);
void cblas_daxpy(blasint n, double alpha, OPENBLAS_CONST double *x, blasint incx, double *y, blasint incy);
void cblas_caxpy(blasint n, OPENBLAS_CONST float *alpha, OPENBLAS_CONST float *x, blasint incx, float *y, blasint incy);
void cblas_zaxpy(blasint n, OPENBLAS_CONST double *alpha, OPENBLAS_CONST double *x, blasint incx, double *y, blasint incy);

void cblas_scopy(blasint n, OPENBLAS_CONST float *x, blasint incx, float *y, blasint incy);
void cblas_dcopy(blasint n, OPENBLAS_CONST double *x, blasint incx, double *y, blasint incy);
void cblas_ccopy(blasint n, OPENBLAS_CONST float *x, blasint incx, float *y, blasint incy);
void cblas_zcopy(blasint n, OPENBLAS_CONST double *x, blasint incx, double *y, blasint incy);

void cblas_sswap(blasint n, float *x, blasint incx, float *y, blasint incy);
void cblas_dswap(blasint n, double *x, blasint incx, double *y, blasint incy);
void cblas_cswap(blasint n, float *x, blasint incx, float *y, blasint incy);
void cblas_zswap(blasint n, double *x, blasint incx, double *y, blasint incy);

void cblas_srot(blasint N, float *X, blasint incX, float *Y, blasint incY, float c, float s);
void cblas_drot(blasint N, double *X, blasint incX, double *Y, blasint incY, double c, double  s);

void cblas_srotg(float *a, float *b, float *c, float *s);
void cblas_drotg(double *a, double *b, double *c, double *s);

void cblas_srotm(blasint N, float *X, blasint incX, float *Y, blasint incY, OPENBLAS_CONST float *P);
void cblas_drotm(blasint N, double *X, blasint incX, double *Y, blasint incY, OPENBLAS_CONST double *P);

void cblas_srotmg(float *d1, float *d2, float *b1, float b2, float *P);
void cblas_drotmg(double *d1, double *d2, double *b1, double b2, double *P);

void cblas_sscal(blasint N, float alpha, float *X, blasint incX);
void cblas_dscal(blasint N, double alpha, double *X, blasint incX);
void cblas_cscal(blasint N, OPENBLAS_CONST float *alpha, float *X, blasint incX);
void cblas_zscal(blasint N, OPENBLAS_CONST double *alpha, double *X, blasint incX);
void cblas_csscal(blasint N, float alpha, float *X, blasint incX);
void cblas_zdscal(blasint N, double alpha, double *X, blasint incX);

void cblas_sgemv(enum CBLAS_ORDER order,  enum CBLAS_TRANSPOSE trans,  blasint m, blasint n,
		 float alpha, OPENBLAS_CONST float  *a, blasint lda,  OPENBLAS_CONST float  *x, blasint incx,  float beta,  float  *y, blasint incy);
void cblas_dgemv(enum CBLAS_ORDER order,  enum CBLAS_TRANSPOSE trans,  blasint m, blasint n,
		 double alpha, OPENBLAS_CONST double  *a, blasint lda,  OPENBLAS_CONST double  *x, blasint incx,  double beta,  double  *y, blasint incy);
void cblas_cgemv(enum CBLAS_ORDER order,  enum CBLAS_TRANSPOSE trans,  blasint m, blasint n,
		 OPENBLAS_CONST float *alpha, OPENBLAS_CONST float  *a, blasint lda,  OPENBLAS_CONST float  *x, blasint incx,  OPENBLAS_CONST float *beta,  float  *y, blasint incy);
void cblas_zgemv(enum CBLAS_ORDER order,  enum CBLAS_TRANSPOSE trans,  blasint m, blasint n,
		 OPENBLAS_CONST double *alpha, OPENBLAS_CONST double  *a, blasint lda,  OPENBLAS_CONST double  *x, blasint incx,  OPENBLAS_CONST double *beta,  double  *y, blasint incy);

void cblas_sger (enum CBLAS_ORDER order, blasint M, blasint N, float   alpha, OPENBLAS_CONST float  *X, blasint incX, OPENBLAS_CONST float  *Y, blasint incY, float  *A, blasint lda);
void cblas_dger (enum CBLAS_ORDER order, blasint M, blasint N, double  alpha, OPENBLAS_CONST double *X, blasint incX, OPENBLAS_CONST double *Y, blasint incY, double *A, blasint lda);
void cblas_cgeru(enum CBLAS_ORDER order, blasint M, blasint N, OPENBLAS_CONST float  *alpha, OPENBLAS_CONST float  *X, blasint incX, OPENBLAS_CONST float  *Y, blasint incY, float  *A, blasint lda);
void cblas_cgerc(enum CBLAS_ORDER order, blasint M, blasint N, OPENBLAS_CONST float  *alpha, OPENBLAS_CONST float  *X, blasint incX, OPENBLAS_CONST float  *Y, blasint incY, float  *A, blasint lda);
void cblas_zgeru(enum CBLAS_ORDER order, blasint M, blasint N, OPENBLAS_CONST double *alpha, OPENBLAS_CONST double *X, blasint incX, OPENBLAS_CONST double *Y, blasint incY, double *A, blasint lda);
void cblas_zgerc(enum CBLAS_ORDER order, blasint M, blasint N, OPENBLAS_CONST double *alpha, OPENBLAS_CONST double *X, blasint incX, OPENBLAS_CONST double *Y, blasint incY, double *A, blasint lda);

void cblas_strsv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag, blasint N, OPENBLAS_CONST float *A, blasint lda, float *X, blasint incX);
void cblas_dtrsv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag, blasint N, OPENBLAS_CONST double *A, blasint lda, double *X, blasint incX);
void cblas_ctrsv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag, blasint N, OPENBLAS_CONST float *A, blasint lda, float *X, blasint incX);
void cblas_ztrsv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag, blasint N, OPENBLAS_CONST double *A, blasint lda, double *X, blasint incX);

void cblas_strmv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag, blasint N, OPENBLAS_CONST float *A, blasint lda, float *X, blasint incX);
void cblas_dtrmv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag, blasint N, OPENBLAS_CONST double *A, blasint lda, double *X, blasint incX);
void cblas_ctrmv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag, blasint N, OPENBLAS_CONST float *A, blasint lda, float *X, blasint incX);
void cblas_ztrmv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag, blasint N, OPENBLAS_CONST double *A, blasint lda, double *X, blasint incX);

void cblas_ssyr(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, float alpha, OPENBLAS_CONST float *X, blasint incX, float *A, blasint lda);
void cblas_dsyr(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, double alpha, OPENBLAS_CONST double *X, blasint incX, double *A, blasint lda);
void cblas_cher(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, float alpha, OPENBLAS_CONST float *X, blasint incX, float *A, blasint lda);
void cblas_zher(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, double alpha, OPENBLAS_CONST double *X, blasint incX, double *A, blasint lda);

void cblas_ssyr2(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo,blasint N, float alpha, OPENBLAS_CONST float *X,
                blasint incX, OPENBLAS_CONST float *Y, blasint incY, float *A, blasint lda);
void cblas_dsyr2(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, double alpha, OPENBLAS_CONST double *X,
                blasint incX, OPENBLAS_CONST double *Y, blasint incY, double *A, blasint lda);
void cblas_cher2(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, OPENBLAS_CONST float *alpha, OPENBLAS_CONST float *X, blasint incX,
                OPENBLAS_CONST float *Y, blasint incY, float *A, blasint lda);
void cblas_zher2(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, OPENBLAS_CONST double *alpha, OPENBLAS_CONST double *X, blasint incX,
                OPENBLAS_CONST double *Y, blasint incY, double *A, blasint lda);

void cblas_sgbmv(enum CBLAS_ORDER order, enum CBLAS_TRANSPOSE TransA, blasint M, blasint N,
                 blasint KL, blasint KU, float alpha, OPENBLAS_CONST float *A, blasint lda, OPENBLAS_CONST float *X, blasint incX, float beta, float *Y, blasint incY);
void cblas_dgbmv(enum CBLAS_ORDER order, enum CBLAS_TRANSPOSE TransA, blasint M, blasint N,
                 blasint KL, blasint KU, double alpha, OPENBLAS_CONST double *A, blasint lda, OPENBLAS_CONST double *X, blasint incX, double beta, double *Y, blasint incY);
void cblas_cgbmv(enum CBLAS_ORDER order, enum CBLAS_TRANSPOSE TransA, blasint M, blasint N,
                 blasint KL, blasint KU, OPENBLAS_CONST float *alpha, OPENBLAS_CONST float *A, blasint lda, OPENBLAS_CONST float *X, blasint incX, OPENBLAS_CONST float *beta, float *Y, blasint incY);
void cblas_zgbmv(enum CBLAS_ORDER order, enum CBLAS_TRANSPOSE TransA, blasint M, blasint N,
                 blasint KL, blasint KU, OPENBLAS_CONST double *alpha, OPENBLAS_CONST double *A, blasint lda, OPENBLAS_CONST double *X, blasint incX, OPENBLAS_CONST double *beta, double *Y, blasint incY);

void cblas_ssbmv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, blasint K, float alpha, OPENBLAS_CONST float *A,
                 blasint lda, OPENBLAS_CONST float *X, blasint incX, float beta, float *Y, blasint incY);
void cblas_dsbmv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, blasint K, double alpha, OPENBLAS_CONST double *A,
                 blasint lda, OPENBLAS_CONST double *X, blasint incX, double beta, double *Y, blasint incY);


void cblas_stbmv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag,
                 blasint N, blasint K, OPENBLAS_CONST float *A, blasint lda, float *X, blasint incX);
void cblas_dtbmv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag,
                 blasint N, blasint K, OPENBLAS_CONST double *A, blasint lda, double *X, blasint incX);
void cblas_ctbmv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag,
                 blasint N, blasint K, OPENBLAS_CONST float *A, blasint lda, float *X, blasint incX);
void cblas_ztbmv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag,
                 blasint N, blasint K, OPENBLAS_CONST double *A, blasint lda, double *X, blasint incX);

void cblas_stbsv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag,
                 blasint N, blasint K, OPENBLAS_CONST float *A, blasint lda, float *X, blasint incX);
void cblas_dtbsv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag,
                 blasint N, blasint K, OPENBLAS_CONST double *A, blasint lda, double *X, blasint incX);
void cblas_ctbsv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag,
                 blasint N, blasint K, OPENBLAS_CONST float *A, blasint lda, float *X, blasint incX);
void cblas_ztbsv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag,
                 blasint N, blasint K, OPENBLAS_CONST double *A, blasint lda, double *X, blasint incX);

void cblas_stpmv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag,
                 blasint N, OPENBLAS_CONST float *Ap, float *X, blasint incX);
void cblas_dtpmv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag,
                 blasint N, OPENBLAS_CONST double *Ap, double *X, blasint incX);
void cblas_ctpmv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag,
                 blasint N, OPENBLAS_CONST float *Ap, float *X, blasint incX);
void cblas_ztpmv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag,
                 blasint N, OPENBLAS_CONST double *Ap, double *X, blasint incX);

void cblas_stpsv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag,
                 blasint N, OPENBLAS_CONST float *Ap, float *X, blasint incX);
void cblas_dtpsv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag,
                 blasint N, OPENBLAS_CONST double *Ap, double *X, blasint incX);
void cblas_ctpsv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag,
                 blasint N, OPENBLAS_CONST float *Ap, float *X, blasint incX);
void cblas_ztpsv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA, enum CBLAS_DIAG Diag,
                 blasint N, OPENBLAS_CONST double *Ap, double *X, blasint incX);

void cblas_ssymv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, float alpha, OPENBLAS_CONST float *A,
                 blasint lda, OPENBLAS_CONST float *X, blasint incX, float beta, float *Y, blasint incY);
void cblas_dsymv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, double alpha, OPENBLAS_CONST double *A,
                 blasint lda, OPENBLAS_CONST double *X, blasint incX, double beta, double *Y, blasint incY);
void cblas_chemv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, OPENBLAS_CONST float *alpha, OPENBLAS_CONST float *A,
                 blasint lda, OPENBLAS_CONST float *X, blasint incX, OPENBLAS_CONST float *beta, float *Y, blasint incY);
void cblas_zhemv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, OPENBLAS_CONST double *alpha, OPENBLAS_CONST double *A,
                 blasint lda, OPENBLAS_CONST double *X, blasint incX, OPENBLAS_CONST double *beta, double *Y, blasint incY);


void cblas_sspmv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, float alpha, OPENBLAS_CONST float *Ap,
                 OPENBLAS_CONST float *X, blasint incX, float beta, float *Y, blasint incY);
void cblas_dspmv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, double alpha, OPENBLAS_CONST double *Ap,
                 OPENBLAS_CONST double *X, blasint incX, double beta, double *Y, blasint incY);

void cblas_sspr(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, float alpha, OPENBLAS_CONST float *X, blasint incX, float *Ap);
void cblas_dspr(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, double alpha, OPENBLAS_CONST double *X, blasint incX, double *Ap);

void cblas_chpr(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, float alpha, OPENBLAS_CONST float *X, blasint incX, float *A);
void cblas_zhpr(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, double alpha, OPENBLAS_CONST double *X,blasint incX, double *A);

void cblas_sspr2(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, float alpha, OPENBLAS_CONST float *X, blasint incX, OPENBLAS_CONST float *Y, blasint incY, float *A);
void cblas_dspr2(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, double alpha, OPENBLAS_CONST double *X, blasint incX, OPENBLAS_CONST double *Y, blasint incY, double *A);
void cblas_chpr2(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, OPENBLAS_CONST float *alpha, OPENBLAS_CONST float *X, blasint incX, OPENBLAS_CONST float *Y, blasint incY, float *Ap);
void cblas_zhpr2(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, OPENBLAS_CONST double *alpha, OPENBLAS_CONST double *X, blasint incX, OPENBLAS_CONST double *Y, blasint incY, double *Ap);

void cblas_chbmv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, blasint K,
		 OPENBLAS_CONST float *alpha, OPENBLAS_CONST float *A, blasint lda, OPENBLAS_CONST float *X, blasint incX, OPENBLAS_CONST float *beta, float *Y, blasint incY);
void cblas_zhbmv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N, blasint K,
		 OPENBLAS_CONST double *alpha, OPENBLAS_CONST double *A, blasint lda, OPENBLAS_CONST double *X, blasint incX, OPENBLAS_CONST double *beta, double *Y, blasint incY);

void cblas_chpmv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N,
		 OPENBLAS_CONST float *alpha, OPENBLAS_CONST float *Ap, OPENBLAS_CONST float *X, blasint incX, OPENBLAS_CONST float *beta, float *Y, blasint incY);
void cblas_zhpmv(enum CBLAS_ORDER order, enum CBLAS_UPLO Uplo, blasint N,
		 OPENBLAS_CONST double *alpha, OPENBLAS_CONST double *Ap, OPENBLAS_CONST double *X, blasint incX, OPENBLAS_CONST double *beta, double *Y, blasint incY);

void cblas_sgemm(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, blasint M, blasint N, blasint K,
		 float alpha, OPENBLAS_CONST float *A, blasint lda, OPENBLAS_CONST float *B, blasint ldb, float beta, float *C, blasint ldc);
void cblas_dgemm(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, blasint M, blasint N, blasint K,
		 double alpha, OPENBLAS_CONST double *A, blasint lda, OPENBLAS_CONST double *B, blasint ldb, double beta, double *C, blasint ldc);
void cblas_cgemm(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, blasint M, blasint N, blasint K,
		 OPENBLAS_CONST float *alpha, OPENBLAS_CONST float *A, blasint lda, OPENBLAS_CONST float *B, blasint ldb, OPENBLAS_CONST float *beta, float *C, blasint ldc);
void cblas_cgemm3m(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, blasint M, blasint N, blasint K,
		 OPENBLAS_CONST float *alpha, OPENBLAS_CONST float *A, blasint lda, OPENBLAS_CONST float *B, blasint ldb, OPENBLAS_CONST float *beta, float *C, blasint ldc);
void cblas_zgemm(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, blasint M, blasint N, blasint K,
		 OPENBLAS_CONST double *alpha, OPENBLAS_CONST double *A, blasint lda, OPENBLAS_CONST double *B, blasint ldb, OPENBLAS_CONST double *beta, double *C, blasint ldc);
void cblas_zgemm3m(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, blasint M, blasint N, blasint K,
		 OPENBLAS_CONST double *alpha, OPENBLAS_CONST double *A, blasint lda, OPENBLAS_CONST double *B, blasint ldb, OPENBLAS_CONST double *beta, double *C, blasint ldc);


void cblas_ssymm(enum CBLAS_ORDER Order, enum CBLAS_SIDE Side, enum CBLAS_UPLO Uplo, blasint M, blasint N,
                 float alpha, OPENBLAS_CONST float *A, blasint lda, OPENBLAS_CONST float *B, blasint ldb, float beta, float *C, blasint ldc);
void cblas_dsymm(enum CBLAS_ORDER Order, enum CBLAS_SIDE Side, enum CBLAS_UPLO Uplo, blasint M, blasint N,
                 double alpha, OPENBLAS_CONST double *A, blasint lda, OPENBLAS_CONST double *B, blasint ldb, double beta, double *C, blasint ldc);
void cblas_csymm(enum CBLAS_ORDER Order, enum CBLAS_SIDE Side, enum CBLAS_UPLO Uplo, blasint M, blasint N,
                 OPENBLAS_CONST float *alpha, OPENBLAS_CONST float *A, blasint lda, OPENBLAS_CONST float *B, blasint ldb, OPENBLAS_CONST float *beta, float *C, blasint ldc);
void cblas_zsymm(enum CBLAS_ORDER Order, enum CBLAS_SIDE Side, enum CBLAS_UPLO Uplo, blasint M, blasint N,
                 OPENBLAS_CONST double *alpha, OPENBLAS_CONST double *A, blasint lda, OPENBLAS_CONST double *B, blasint ldb, OPENBLAS_CONST double *beta, double *C, blasint ldc);

void cblas_ssyrk(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE Trans,
		 blasint N, blasint K, float alpha, OPENBLAS_CONST float *A, blasint lda, float beta, float *C, blasint ldc);
void cblas_dsyrk(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE Trans,
		 blasint N, blasint K, double alpha, OPENBLAS_CONST double *A, blasint lda, double beta, double *C, blasint ldc);
void cblas_csyrk(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE Trans,
		 blasint N, blasint K, OPENBLAS_CONST float *alpha, OPENBLAS_CONST float *A, blasint lda, OPENBLAS_CONST float *beta, float *C, blasint ldc);
void cblas_zsyrk(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE Trans,
		 blasint N, blasint K, OPENBLAS_CONST double *alpha, OPENBLAS_CONST double *A, blasint lda, OPENBLAS_CONST double *beta, double *C, blasint ldc);

void cblas_ssyr2k(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE Trans,
		  blasint N, blasint K, float alpha, OPENBLAS_CONST float *A, blasint lda, OPENBLAS_CONST float *B, blasint ldb, float beta, float *C, blasint ldc);
void cblas_dsyr2k(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE Trans,
		  blasint N, blasint K, double alpha, OPENBLAS_CONST double *A, blasint lda, OPENBLAS_CONST double *B, blasint ldb, double beta, double *C, blasint ldc);
void cblas_csyr2k(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE Trans,
		  blasint N, blasint K, OPENBLAS_CONST float *alpha, OPENBLAS_CONST float *A, blasint lda, OPENBLAS_CONST float *B, blasint ldb, OPENBLAS_CONST float *beta, float *C, blasint ldc);
void cblas_zsyr2k(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE Trans,
		  blasint N, blasint K, OPENBLAS_CONST double *alpha, OPENBLAS_CONST double *A, blasint lda, OPENBLAS_CONST double *B, blasint ldb, OPENBLAS_CONST double *beta, double *C, blasint ldc);

void cblas_strmm(enum CBLAS_ORDER Order, enum CBLAS_SIDE Side, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA,
                 enum CBLAS_DIAG Diag, blasint M, blasint N, float alpha, OPENBLAS_CONST float *A, blasint lda, float *B, blasint ldb);
void cblas_dtrmm(enum CBLAS_ORDER Order, enum CBLAS_SIDE Side, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA,
                 enum CBLAS_DIAG Diag, blasint M, blasint N, double alpha, OPENBLAS_CONST double *A, blasint lda, double *B, blasint ldb);
void cblas_ctrmm(enum CBLAS_ORDER Order, enum CBLAS_SIDE Side, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA,
                 enum CBLAS_DIAG Diag, blasint M, blasint N, OPENBLAS_CONST float *alpha, OPENBLAS_CONST float *A, blasint lda, float *B, blasint ldb);
void cblas_ztrmm(enum CBLAS_ORDER Order, enum CBLAS_SIDE Side, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA,
                 enum CBLAS_DIAG Diag, blasint M, blasint N, OPENBLAS_CONST double *alpha, OPENBLAS_CONST double *A, blasint lda, double *B, blasint ldb);

void cblas_strsm(enum CBLAS_ORDER Order, enum CBLAS_SIDE Side, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA,
                 enum CBLAS_DIAG Diag, blasint M, blasint N, float alpha, OPENBLAS_CONST float *A, blasint lda, float *B, blasint ldb);
void cblas_dtrsm(enum CBLAS_ORDER Order, enum CBLAS_SIDE Side, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA,
                 enum CBLAS_DIAG Diag, blasint M, blasint N, double alpha, OPENBLAS_CONST double *A, blasint lda, double *B, blasint ldb);
void cblas_ctrsm(enum CBLAS_ORDER Order, enum CBLAS_SIDE Side, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA,
                 enum CBLAS_DIAG Diag, blasint M, blasint N, OPENBLAS_CONST float *alpha, OPENBLAS_CONST float *A, blasint lda, float *B, blasint ldb);
void cblas_ztrsm(enum CBLAS_ORDER Order, enum CBLAS_SIDE Side, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA,
                 enum CBLAS_DIAG Diag, blasint M, blasint N, OPENBLAS_CONST double *alpha, OPENBLAS_CONST double *A, blasint lda, double *B, blasint ldb);

void cblas_chemm(enum CBLAS_ORDER Order, enum CBLAS_SIDE Side, enum CBLAS_UPLO Uplo, blasint M, blasint N,
                 OPENBLAS_CONST float *alpha, OPENBLAS_CONST float *A, blasint lda, OPENBLAS_CONST float *B, blasint ldb, OPENBLAS_CONST float *beta, float *C, blasint ldc);
void cblas_zhemm(enum CBLAS_ORDER Order, enum CBLAS_SIDE Side, enum CBLAS_UPLO Uplo, blasint M, blasint N,
                 OPENBLAS_CONST double *alpha, OPENBLAS_CONST double *A, blasint lda, OPENBLAS_CONST double *B, blasint ldb, OPENBLAS_CONST double *beta, double *C, blasint ldc);

void cblas_cherk(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE Trans, blasint N, blasint K,
                 float alpha, OPENBLAS_CONST float *A, blasint lda, float beta, float *C, blasint ldc);
void cblas_zherk(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE Trans, blasint N, blasint K,
                 double alpha, OPENBLAS_CONST double *A, blasint lda, double beta, double *C, blasint ldc);

void cblas_cher2k(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE Trans, blasint N, blasint K,
                  OPENBLAS_CONST float *alpha, OPENBLAS_CONST float *A, blasint lda, OPENBLAS_CONST float *B, blasint ldb, float beta, float *C, blasint ldc);
void cblas_zher2k(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE Trans, blasint N, blasint K,
                  OPENBLAS_CONST double *alpha, OPENBLAS_CONST double *A, blasint lda, OPENBLAS_CONST double *B, blasint ldb, double beta, double *C, blasint ldc);

void cblas_xerbla(blasint p, char *rout, char *form, ...);

/*** BLAS extensions ***/

void cblas_saxpby(blasint n, float alpha, OPENBLAS_CONST float *x, blasint incx,float beta, float *y, blasint incy);

void cblas_daxpby(blasint n, double alpha, OPENBLAS_CONST double *x, blasint incx,double beta, double *y, blasint incy);

void cblas_caxpby(blasint n, OPENBLAS_CONST float *alpha, OPENBLAS_CONST float *x, blasint incx,OPENBLAS_CONST float *beta, float *y, blasint incy);

void cblas_zaxpby(blasint n, OPENBLAS_CONST double *alpha, OPENBLAS_CONST double *x, blasint incx,OPENBLAS_CONST double *beta, double *y, blasint incy);

void cblas_somatcopy(enum CBLAS_ORDER CORDER, enum CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, float calpha, OPENBLAS_CONST float *a,
		     blasint clda, float *b, blasint cldb);
void cblas_domatcopy(enum CBLAS_ORDER CORDER, enum CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, double calpha, OPENBLAS_CONST double *a,
		     blasint clda, double *b, blasint cldb);
void cblas_comatcopy(enum CBLAS_ORDER CORDER, enum CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, OPENBLAS_CONST float* calpha, OPENBLAS_CONST float* a,
		     blasint clda, float*b, blasint cldb);
void cblas_zomatcopy(enum CBLAS_ORDER CORDER, enum CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, OPENBLAS_CONST double* calpha, OPENBLAS_CONST double* a,
		     blasint clda,  double *b, blasint cldb);

void cblas_simatcopy(enum CBLAS_ORDER CORDER, enum CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, float calpha, float *a,
		     blasint clda, blasint cldb);
void cblas_dimatcopy(enum CBLAS_ORDER CORDER, enum CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, double calpha, double *a,
		     blasint clda, blasint cldb);
void cblas_cimatcopy(enum CBLAS_ORDER CORDER, enum CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, OPENBLAS_CONST float* calpha, float* a,
		     blasint clda, blasint cldb);
void cblas_zimatcopy(enum CBLAS_ORDER CORDER, enum CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, OPENBLAS_CONST double* calpha, double* a,
		     blasint clda, blasint cldb);

void cblas_sgeadd(enum CBLAS_ORDER CORDER,blasint crows, blasint ccols, float calpha, float *a, blasint clda, float cbeta,
		  float *c, blasint cldc);
void cblas_dgeadd(enum CBLAS_ORDER CORDER,blasint crows, blasint ccols, double calpha, double *a, blasint clda, double cbeta,
		  double *c, blasint cldc);
void cblas_cgeadd(enum CBLAS_ORDER CORDER,blasint crows, blasint ccols, OPENBLAS_CONST float *calpha, float *a, blasint clda, OPENBLAS_CONST float *cbeta,
		  float *c, blasint cldc);
void cblas_zgeadd(enum CBLAS_ORDER CORDER,blasint crows, blasint ccols, OPENBLAS_CONST double *calpha, double *a, blasint clda, OPENBLAS_CONST double *cbeta,
		  double *c, blasint cldc);


#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
