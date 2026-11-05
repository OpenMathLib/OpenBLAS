/*****************************************************************************
Copyright (c) 2011-2022, The OpenBLAS Project
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
#include <math.h>
#include "openblas_utest.h"
#if defined(BUILD_DOUBLE)

#ifndef INFINITY
#define INFINITY HUGE_VAL
#endif

CTEST(dnrm2,dnrm2_inf)
{
	int i;
	double x[29];
	blasint incx=1;
	blasint n=28;
	double res1=0.0f, res2=INFINITY;

	for (i=0;i<n;i++)x[i]=0.0f;
	x[10]=-INFINITY;
	res1=BLASFUNC(dnrm2)(&n, x, &incx);
	ASSERT_DBL_NEAR_TOL(res2, res1, DOUBLE_EPS);

}
CTEST(dnrm2,dnrm2_tiny)
{
	int i;
	double x[29];
	blasint incx=1;
	blasint n=28;
	double res1=0.0f, res2=0.0f;

	for (i=0;i<n;i++)x[i]=7.457008414e-310;
	res1=BLASFUNC(dnrm2)(&n, x, &incx);
	ASSERT_DBL_NEAR_TOL(res2, res1, DOUBLE_EPS);
}
CTEST(dnrm2,dnrm2_neg_incx)
{
	int i;
	double x[5];
	blasint incx=-1;
	blasint n=5;
	double res1, res2;

	for (i=0;i<n;i++)x[i]=10.0;
	res1=BLASFUNC(dnrm2)(&n, x, &incx);
	res2 = sqrt(500.0);
	ASSERT_DBL_NEAR_TOL(res2, res1, DOUBLE_EPS);
}
#endif

/* Negative increments are supported by the public NRM2 interfaces. Keep
   complex RVV byte strides signed when the wrapper starts at the last item. */
#if !defined(NO_CBLAS) && (defined(BUILD_COMPLEX) || defined(BUILD_COMPLEX16))
#include "cblas.h"
#endif

#if defined(BUILD_COMPLEX)
static void check_cnrm2_stride(int c_api, blasint increment) {
    blasint n = 3;
    float x[10] = {3.0f, 4.0f, 0.0f, 0.0f, 0.0f,
                   0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    float result;
    (void)c_api;
#ifndef NO_CBLAS
    if (c_api) result = cblas_scnrm2(n, x, increment);
    else
#endif
    result = BLASFUNC(scnrm2)(&n, x, &increment);
    ASSERT_TRUE(isfinite(result));
    ASSERT_DBL_NEAR_TOL(5.0, result, SINGLE_EPS);
}
CTEST(cnrm2, negative_unit_fortran) { check_cnrm2_stride(0, -1); }
CTEST(cnrm2, negative_nonunit_fortran) { check_cnrm2_stride(0, -2); }
CTEST(cnrm2, positive_nonunit_fortran) { check_cnrm2_stride(0, 2); }
#ifndef NO_CBLAS
CTEST(cnrm2, negative_unit_cblas) { check_cnrm2_stride(1, -1); }
CTEST(cnrm2, negative_nonunit_cblas) { check_cnrm2_stride(1, -2); }
CTEST(cnrm2, positive_nonunit_cblas) { check_cnrm2_stride(1, 2); }
#endif
#endif

#if defined(BUILD_COMPLEX16)
static void check_znrm2_stride(int c_api, blasint increment) {
    blasint n = 3;
    double x[10] = {3.0, 4.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0};
    double result;
    (void)c_api;
#ifndef NO_CBLAS
    if (c_api) result = cblas_dznrm2(n, x, increment);
    else
#endif
    result = BLASFUNC(dznrm2)(&n, x, &increment);
    ASSERT_TRUE(isfinite(result));
    ASSERT_DBL_NEAR_TOL(5.0, result, DOUBLE_EPS);
}
CTEST(znrm2, negative_unit_fortran) { check_znrm2_stride(0, -1); }
CTEST(znrm2, negative_nonunit_fortran) { check_znrm2_stride(0, -2); }
CTEST(znrm2, positive_nonunit_fortran) { check_znrm2_stride(0, 2); }
#ifndef NO_CBLAS
CTEST(znrm2, negative_unit_cblas) { check_znrm2_stride(1, -1); }
CTEST(znrm2, negative_nonunit_cblas) { check_znrm2_stride(1, -2); }
CTEST(znrm2, positive_nonunit_cblas) { check_znrm2_stride(1, 2); }
#endif
#endif
