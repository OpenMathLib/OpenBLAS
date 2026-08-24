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
#if defined(BUILD_SINGLE) && defined(BUILD_DOUBLE)
CTEST(dsdot,dsdot_n_1)
{
	float x= 0.172555164F;
	float y= -0.0138700781F;
	blasint incx=1;
	blasint incy=1;
	blasint n=1;

	double res1=0.0f, res2=-0.00239335360107;

	res1=BLASFUNC(dsdot)(&n, &x, &incx, &y, &incy);
	ASSERT_DBL_NEAR_TOL(res2, res1, DOUBLE_EPS);

}
#endif
#ifdef ARMV8
#if defined(BUILD_SINGLE)
CTEST(sdot,sdot_n_1)
{
    static float x[64], y[64];
    for (int i = 0; i < 64; i++) { x[i] = 1.0f; y[i] = 1.0f; }
    blasint n = 64, inc = 1;
    float junk[4] = {1e6f, 1e6f, 1e6f, 1e6f};
    __asm__ volatile("ld1 {v0.4s}, [%0]" :: "r"(junk) : "v0");
    float r = BLASFUNC(sdot)(&n, x, &inc, y, &inc);
    ASSERT_DBL_NEAR_TOL(64.,r, DOUBLE_EPS);
}
#endif
#if defined(BUILD_DOUBLE)
CTEST(ddot,ddot_n_1)
{
    static double x[64], y[64];
    for (int i = 0; i < 64; i++) { x[i] = 1.0f; y[i] = 1.0f; }
    blasint n = 64, inc = 1;
    double junk[4] = {1e6f, 1e6f, 1e6f, 1e6f};
    __asm__ volatile("ld1 {v0.4s}, [%0]" :: "r"(junk) : "v0");
    double r = BLASFUNC(ddot)(&n, x, &inc, y, &inc);
    ASSERT_DBL_NEAR_TOL(64.,r, DOUBLE_EPS);
}
#endif
#endif
