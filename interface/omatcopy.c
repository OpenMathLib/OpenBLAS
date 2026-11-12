/***************************************************************************
Copyright (c) 2014, The OpenBLAS Project
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
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

/***********************************************************
 * 2014/06/09 Saar
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#if defined(DOUBLE)
#define ERROR_NAME "DOMATCOPY"
#else
#define ERROR_NAME "SOMATCOPY"
#endif

#define BlasRowMajor 0
#define BlasColMajor 1
#define BlasNoTrans  0
#define BlasTrans    1

#ifdef SMP
#ifndef GEMM_MULTITHREAD_THRESHOLD
#define GEMM_MULTITHREAD_THRESHOLD 4
#endif

static int omatcopy_rn(BLASLONG m, BLASLONG n, BLASLONG dummy0, FLOAT alpha, FLOAT *a, BLASLONG lda,
		       FLOAT *b, BLASLONG ldb, FLOAT *dummy1, BLASLONG dummy2, void *dummy3)
{
	return OMATCOPY_K_RN(m, n, alpha, a, lda, b, ldb);
}

static int omatcopy_rt(BLASLONG m, BLASLONG n, BLASLONG dummy0, FLOAT alpha, FLOAT *a, BLASLONG lda,
		       FLOAT *b, BLASLONG ldb, FLOAT *dummy1, BLASLONG dummy2, void *dummy3)
{
	return OMATCOPY_K_RT(m, n, alpha, a, lda, b, ldb);
}

static int omatcopy_cn(BLASLONG m, BLASLONG n, BLASLONG dummy0, FLOAT alpha, FLOAT *a, BLASLONG lda,
		       FLOAT *b, BLASLONG ldb, FLOAT *dummy1, BLASLONG dummy2, void *dummy3)
{
	return OMATCOPY_K_CN(n, m, alpha, a, lda, b, ldb);
}

static int omatcopy_ct(BLASLONG m, BLASLONG n, BLASLONG dummy0, FLOAT alpha, FLOAT *a, BLASLONG lda,
		       FLOAT *b, BLASLONG ldb, FLOAT *dummy1, BLASLONG dummy2, void *dummy3)
{
	return OMATCOPY_K_CT(n, m, alpha, a, lda, b, ldb);
}
#endif

OPENBLAS_EXPORT

#ifndef CBLAS 
void NAME( char* ORDER, char* TRANS, blasint *rows, blasint *cols, FLOAT *alpha, FLOAT *a, blasint *lda, FLOAT *b, blasint *ldb)
{

	char Order, Trans;
	int order=-1,trans=-1;
	blasint info = -1;
#ifdef SMP
	int nthreads;
#endif

	Order = *ORDER;
	Trans = *TRANS;

	TOUPPER(Order);
	TOUPPER(Trans);

	if ( Order == 'C' ) order = BlasColMajor;
	if ( Order == 'R' ) order = BlasRowMajor;
	if ( Trans == 'N' ) trans = BlasNoTrans;
	if ( Trans == 'R' ) trans = BlasNoTrans;
	if ( Trans == 'T' ) trans = BlasTrans;
	if ( Trans == 'C' ) trans = BlasTrans;
#else 
void CNAME(enum CBLAS_ORDER CORDER, enum CBLAS_TRANSPOSE CTRANS, blasint crows, blasint ccols, FLOAT calpha, FLOAT *a, blasint clda, FLOAT *b, blasint cldb)
{
	blasint *rows, *cols, *lda, *ldb; 
	FLOAT   *alpha; 
	int order=-1,trans=-1;
	blasint info = -1;
#ifdef SMP
	int nthreads;
#endif

	if ( CORDER == CblasColMajor ) order = BlasColMajor; 
	if ( CORDER == CblasRowMajor ) order = BlasRowMajor; 

	if ( CTRANS == CblasNoTrans || CTRANS == CblasConjNoTrans ) trans = BlasNoTrans; 
	if ( CTRANS == CblasTrans   || CTRANS == CblasConjTrans   ) trans = BlasTrans; 

	rows = &crows; 
	cols = &ccols; 
	lda  = &clda; 
	ldb  = &cldb; 
	alpha = &calpha; 

#endif
	if ( order == BlasColMajor)
	{
        	if ( trans == BlasNoTrans  &&  *ldb < MAX(1,*rows) ) info = 9;
        	if ( trans == BlasTrans    &&  *ldb < MAX(1,*cols) ) info = 9;
	}
	if ( order == BlasRowMajor)
	{
        	if ( trans == BlasNoTrans  &&  *ldb < MAX(1,*cols) ) info = 9;
        	if ( trans == BlasTrans    &&  *ldb < MAX(1,*rows) ) info = 9;
	}

	if ( order == BlasColMajor &&  *lda < MAX(1,*rows) ) info = 7;
	if ( order == BlasRowMajor &&  *lda < MAX(1,*cols) ) info = 7;
	if ( *cols < 0 ) info = 4;
	if ( *rows < 0 ) info = 3;
	if ( trans < 0 ) info = 2;
	if ( order < 0 ) info = 1;

	if (info >= 0) {
    		BLASFUNC(xerbla)(ERROR_NAME, &info, sizeof(ERROR_NAME));
    		return;
  	}

	if ((*rows == 0) || (*cols == 0)) return;

#ifdef SMP
	nthreads = num_cpu_avail(1);
	/* Small copies stay on one thread */
	if ((double)*rows * (double)*cols < 65536. * GEMM_MULTITHREAD_THRESHOLD) nthreads = 1;
	/* In the transposed copies each thread writes a narrow band of every row
	   (or column) of B. Keep the bands at least 32 elements (256 bytes of double,
	   128 of float) wide, otherwise the threads keep writing into the same cache
	   lines. */
	if (trans == BlasTrans) {
		blasint split = (order == BlasRowMajor) ? *rows : *cols;
		if (nthreads > split / 32) nthreads = MAX(1, split / 32);
	}

	if (nthreads > 1) {
#ifdef DOUBLE
		int mode = BLAS_DOUBLE | BLAS_REAL;
#else
		int mode = BLAS_SINGLE | BLAS_REAL;
#endif
		/* Each thread gets a block of rows (row-major) or columns (column-major)
		   of A, so the column-major wrappers swap m and n. For the transposed
		   copies a block of B then starts width elements further along, not width
		   rows or columns, which is what BLAS_TRANSB_T does. */
		if (order == BlasRowMajor) {
			if (trans == BlasNoTrans)
				blas_level1_thread(mode, *rows, *cols, 0, alpha, a, *lda, b, *ldb, NULL, 0,
						   (int (*)(void))omatcopy_rn, nthreads);
			else
				blas_level1_thread(mode | BLAS_TRANSB_T, *rows, *cols, 0, alpha, a, *lda, b, *ldb, NULL, 0,
						   (int (*)(void))omatcopy_rt, nthreads);
		} else {
			if (trans == BlasNoTrans)
				blas_level1_thread(mode, *cols, *rows, 0, alpha, a, *lda, b, *ldb, NULL, 0,
						   (int (*)(void))omatcopy_cn, nthreads);
			else
				blas_level1_thread(mode | BLAS_TRANSB_T, *cols, *rows, 0, alpha, a, *lda, b, *ldb, NULL, 0,
						   (int (*)(void))omatcopy_ct, nthreads);
		}
		return;
	}
#endif

	if ( order == BlasColMajor )
	{
		if ( trans == BlasNoTrans )
		{
			OMATCOPY_K_CN(*rows, *cols, *alpha, a, *lda, b, *ldb );
		}
		else
		{
			OMATCOPY_K_CT(*rows, *cols, *alpha, a, *lda, b, *ldb );
		}
	}
	else
	{
		if ( trans == BlasNoTrans )
		{
			OMATCOPY_K_RN(*rows, *cols, *alpha, a, *lda, b, *ldb );
		}
		else
		{
			OMATCOPY_K_RT(*rows, *cols, *alpha, a, *lda, b, *ldb );
		}
	}

	return;

}


