/*****************************************************************************
Copyright (c) 2011-2014, The OpenBLAS Project
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
      derived from this software without specific prior written 
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

**********************************************************************************/

#include "openblas_utest.h"

#ifdef BUILD_COMPLEX16
CTEST( zdotu,zdotu_n_1)
{
	blasint N=1,incX=1,incY=1;
	double x1[]={1.0,1.0};
	double y1[]={1.0,2.0};
	
	openblas_complex_double result1=openblas_make_complex_double(0.0,0.0);
        openblas_complex_double result2=openblas_make_complex_double(-1.0000,3.0000);
#ifdef RETURN_BY_STACK
	BLASFUNC(zdotu)(&result1,&N,x1,&incX,y1,&incY);
#else
	result1=BLASFUNC(zdotu)(&N,x1,&incX,y1,&incY);
#endif
	
#ifdef OPENBLAS_COMPLEX_STRUCT
	ASSERT_DBL_NEAR_TOL(result2.real, result1.real, DOUBLE_EPS);
	ASSERT_DBL_NEAR_TOL(result2.imag, result1.imag, DOUBLE_EPS);
#else
	ASSERT_DBL_NEAR_TOL(creal(result2), creal(result1), DOUBLE_EPS);
	ASSERT_DBL_NEAR_TOL(cimag(result2), cimag(result1), DOUBLE_EPS);
#endif
	
}

CTEST(zdotu, zdotu_offset_1)
{
	blasint N=1,incX=1,incY=1;
	double x1[]={1.0,2.0,3.0,4.0};
	double y1[]={5.0,6.0,7.0,8.0};
	
	openblas_complex_double result1=openblas_make_complex_double(0.0,0.0);
        openblas_complex_double result2=openblas_make_complex_double(-9.0,32.0);
#ifdef RETURN_BY_STACK
	BLASFUNC(zdotu)(&result1,&N,x1+1,&incX,y1+1,&incY);
#else
	result1=BLASFUNC(zdotu)(&N,x1+1,&incX,y1+1,&incY);
#endif
	
#ifdef OPENBLAS_COMPLEX_STRUCT
	ASSERT_DBL_NEAR_TOL(result2.real, result1.real, DOUBLE_EPS);
	ASSERT_DBL_NEAR_TOL(result2.imag, result1.imag, DOUBLE_EPS);
#else
	ASSERT_DBL_NEAR_TOL(creal(result2), creal(result1), DOUBLE_EPS);
	ASSERT_DBL_NEAR_TOL(cimag(result2), cimag(result1), DOUBLE_EPS);
#endif

}

CTEST(zdotu, complex16_l1_threaded)
{
	blasint i;
	blasint N = 20005, inc = 1;
	blasint idx;
	double asum, nrm2;
	double *x = (double *)calloc((size_t)(N * 2), sizeof(double));
	double *y = (double *)calloc((size_t)(N * 2), sizeof(double));
	openblas_complex_double dot_res = openblas_make_complex_double(0.0, 0.0);

	ASSERT_NOT_NULL(x);
	ASSERT_NOT_NULL(y);
	for (i = 0; i < N; i++) {
		x[2 * i] = 3.0;
		x[2 * i + 1] = 4.0;
		y[2 * i] = 1.0;
		y[2 * i + 1] = 0.0;
	}

#ifdef RETURN_BY_STACK
	BLASFUNC(zdotu)(&dot_res, &N, x, &inc, y, &inc);
#else
	dot_res = BLASFUNC(zdotu)(&N, x, &inc, y, &inc);
#endif
#ifdef OPENBLAS_COMPLEX_STRUCT
	ASSERT_DBL_NEAR_TOL(3.0 * (double)N, dot_res.real, DOUBLE_EPS);
	ASSERT_DBL_NEAR_TOL(4.0 * (double)N, dot_res.imag, DOUBLE_EPS);
#else
	ASSERT_DBL_NEAR_TOL(3.0 * (double)N, creal(dot_res), DOUBLE_EPS);
	ASSERT_DBL_NEAR_TOL(4.0 * (double)N, cimag(dot_res), DOUBLE_EPS);
#endif

	asum = BLASFUNC(dzasum)(&N, x, &inc);
	ASSERT_DBL_NEAR_TOL(7.0 * (double)N, asum, DOUBLE_EPS * (double)N);

	nrm2 = BLASFUNC(dznrm2)(&N, x, &inc);
	ASSERT_DBL_NEAR_TOL(5.0 * sqrt((double)N), nrm2, 1e-10);

	x[2 * (N - 3)] = 10.0;
	x[2 * (N - 3) + 1] = -12.0;
	idx = BLASFUNC(izamax)(&N, x, &inc);
	ASSERT_EQUAL(N - 2, idx);

	BLASFUNC(zcopy)(&N, x, &inc, y, &inc);
	ASSERT_DBL_NEAR_TOL(10.0, y[2 * (N - 3)], DOUBLE_EPS);
	ASSERT_DBL_NEAR_TOL(-12.0, y[2 * (N - 3) + 1], DOUBLE_EPS);

	free(x);
	free(y);
}
#endif

#ifdef BUILD_COMPLEX
CTEST(zdotu, complex_l1_threaded)
{
	blasint i;
	blasint N = 20005, inc = 1;
	blasint idx;
	float asum, nrm2;
	float *x = (float *)calloc((size_t)(N * 2), sizeof(float));
	float *y = (float *)calloc((size_t)(N * 2), sizeof(float));
	openblas_complex_float dot_res = openblas_make_complex_float(0.0f, 0.0f);

	ASSERT_NOT_NULL(x);
	ASSERT_NOT_NULL(y);
	for (i = 0; i < N; i++) {
		x[2 * i] = 0.5f;
		x[2 * i + 1] = 0.5f;
		y[2 * i] = 1.0f;
		y[2 * i + 1] = 0.0f;
	}

#ifdef RETURN_BY_STACK
	BLASFUNC(cdotu)(&dot_res, &N, x, &inc, y, &inc);
#else
	dot_res = BLASFUNC(cdotu)(&N, x, &inc, y, &inc);
#endif
#ifdef OPENBLAS_COMPLEX_STRUCT
	ASSERT_DBL_NEAR_TOL(0.5 * (double)N, (double)dot_res.real, 1e-1);
	ASSERT_DBL_NEAR_TOL(0.5 * (double)N, (double)dot_res.imag, 1e-1);
#else
	ASSERT_DBL_NEAR_TOL(0.5 * (double)N, (double)crealf(dot_res), 1e-1);
	ASSERT_DBL_NEAR_TOL(0.5 * (double)N, (double)cimagf(dot_res), 1e-1);
#endif

	asum = BLASFUNC(scasum)(&N, x, &inc);
	ASSERT_DBL_NEAR_TOL((double)N, (double)asum, 1e-1);

	nrm2 = BLASFUNC(scnrm2)(&N, x, &inc);
	ASSERT_DBL_NEAR_TOL(sqrt(0.5 * (double)N), (double)nrm2, 1e-2);

	x[2 * (N - 3)] = 10.0f;
	x[2 * (N - 3) + 1] = -12.0f;
	idx = BLASFUNC(icamax)(&N, x, &inc);
	ASSERT_EQUAL(N - 2, idx);

	BLASFUNC(ccopy)(&N, x, &inc, y, &inc);
	ASSERT_DBL_NEAR_TOL(10.0, (double)y[2 * (N - 3)], SINGLE_EPS);
	ASSERT_DBL_NEAR_TOL(-12.0, (double)y[2 * (N - 3) + 1], SINGLE_EPS);

	free(x);
	free(y);
}
#endif

