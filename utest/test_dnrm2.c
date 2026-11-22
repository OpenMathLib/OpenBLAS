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

#if (defined(__x86_64__) || (defined(__i386__) && !defined(OPENBLAS_GENERIC))) && (defined(__GNUC__) || defined(__clang__))
CTEST(dnrm2,dnrm2_x87_precision_control)
{
	int i;
	blasint n = 4 + 4096;
	blasint n_pair = 2;
	blasint n_zero = 0;
	blasint incx = 1;
	blasint inc_zero = 0;
	double x[4 + 4096];
	double x_pair[2];
	unsigned short cw_orig = 0, cw_64 = 0, cw_53 = 0, cw_24 = 0, cw_after = 0;
	double res_64, res_53, res_24, res_pair;

	for (i = 0; i < 4; i++) x[i] = 1.0;
	for (i = 4; i < n; i++) x[i] = ldexp(1.0, -27);
	x_pair[0] = 1.0;
	x_pair[1] = ldexp(1.0, -15);

	__asm__ __volatile__("fnstcw %0" : "=m"(cw_orig));

	cw_64 = (unsigned short)(cw_orig | 0x0300u);
	__asm__ __volatile__("fldcw %0" : : "m"(cw_64));
	res_64 = BLASFUNC(dnrm2)(&n, x, &incx);

	/* 53-bit significand (Windows MSVC default 0x027f) must match 64-bit mode */
	cw_53 = (unsigned short)((cw_orig & ~0x0300u) | 0x0200u);
	__asm__ __volatile__("fldcw %0" : : "m"(cw_53));
	res_53 = BLASFUNC(dnrm2)(&n, x, &incx);
	__asm__ __volatile__("fnstcw %0" : "=m"(cw_after));
	__asm__ __volatile__("fldcw %0" : : "m"(cw_orig));
	ASSERT_EQUAL(cw_53, cw_after);
	ASSERT_DBL_NEAR_TOL(res_64, res_53, DOUBLE_EPS);

	/* 24-bit significand (0x007f) and early-return paths */
	cw_24 = (unsigned short)(cw_orig & ~0x0300u);
	__asm__ __volatile__("fldcw %0" : : "m"(cw_24));
	res_24 = BLASFUNC(dnrm2)(&n, x, &incx);
	__asm__ __volatile__("fnstcw %0" : "=m"(cw_after));
	ASSERT_EQUAL(cw_24, cw_after);
	ASSERT_DBL_NEAR_TOL(res_64, res_24, DOUBLE_EPS);

	res_pair = BLASFUNC(dnrm2)(&n_pair, x_pair, &incx);
	__asm__ __volatile__("fnstcw %0" : "=m"(cw_after));
	ASSERT_EQUAL(cw_24, cw_after);
	ASSERT_DBL_NEAR_TOL(sqrt(1.0 + ldexp(1.0, -30)), res_pair, DOUBLE_EPS);

	res_24 = BLASFUNC(dnrm2)(&n_zero, x, &incx);
	__asm__ __volatile__("fnstcw %0" : "=m"(cw_after));
	ASSERT_EQUAL(cw_24, cw_after);
	ASSERT_DBL_NEAR_TOL(0.0, res_24, DOUBLE_EPS);

	res_24 = BLASFUNC(dnrm2)(&n, x, &inc_zero);
	__asm__ __volatile__("fnstcw %0" : "=m"(cw_after));
	__asm__ __volatile__("fldcw %0" : : "m"(cw_orig));
	ASSERT_EQUAL(cw_24, cw_after);
}
#endif
#endif

#if defined(BUILD_COMPLEX16)
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
CTEST(dnrm2,dznrm2_x87_precision_control)
{
	int i;
	blasint n = 2 + 2048;
	blasint incx = 1;
	double x[(2 + 2048) * 2];
	unsigned short cw_orig = 0, cw_64 = 0, cw_53 = 0, cw_after = 0;
	double res_64, res_53;

	for (i = 0; i < 4; i++) x[i] = 1.0;
	for (i = 4; i < n * 2; i++) x[i] = ldexp(1.0, -27);

	__asm__ __volatile__("fnstcw %0" : "=m"(cw_orig));

	cw_64 = (unsigned short)(cw_orig | 0x0300u);
	__asm__ __volatile__("fldcw %0" : : "m"(cw_64));
	res_64 = BLASFUNC(dznrm2)(&n, x, &incx);

	cw_53 = (unsigned short)((cw_orig & ~0x0300u) | 0x0200u);
	__asm__ __volatile__("fldcw %0" : : "m"(cw_53));
	res_53 = BLASFUNC(dznrm2)(&n, x, &incx);
	__asm__ __volatile__("fnstcw %0" : "=m"(cw_after));
	__asm__ __volatile__("fldcw %0" : : "m"(cw_orig));

	ASSERT_EQUAL(cw_53, cw_after);
	ASSERT_DBL_NEAR_TOL(res_64, res_53, DOUBLE_EPS);
}
#endif
#endif
