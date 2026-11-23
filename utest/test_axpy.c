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

#ifdef BUILD_DOUBLE
CTEST(axpy,daxpy_inc_0)
{
	blasint i;
	blasint N=8,incX=0,incY=0;
	double a=0.25;
	double x1[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	double y1[]={2.0,4.0,6.0,8.0,2.0,4.0,6.0,8.0};

	double x2[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	double y2[]={4.0,4.0,6.0,8.0,2.0,4.0,6.0,8.0};

	//OpenBLAS
	BLASFUNC(daxpy)(&N,&a,x1,&incX,y1,&incY);

	for(i=0; i<N; i++){
		ASSERT_DBL_NEAR_TOL(x2[i], x1[i], DOUBLE_EPS);
		ASSERT_DBL_NEAR_TOL(y2[i], y1[i], DOUBLE_EPS);
	}
}
#endif

#ifdef BUILD_COMPLEX16
CTEST(axpy,zaxpy_inc_0)
{
	blasint i;
	blasint N=4,incX=0,incY=0;
	double a[2]={0.25,0.5};
	double x1[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	double y1[]={2.0,4.0,6.0,8.0,2.0,4.0,6.0,8.0};
	double x2[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	double y2[]={-3.0,9.0,6.0,8.0,2.0,4.0,6.0,8.0};

	//OpenBLAS
	BLASFUNC(zaxpy)(&N,a,x1,&incX,y1,&incY);

	for(i=0; i<2*N; i++){
		ASSERT_DBL_NEAR_TOL(x2[i], x1[i], DOUBLE_EPS);
		ASSERT_DBL_NEAR_TOL(y2[i], y1[i], DOUBLE_EPS);
	}
}

CTEST(axpy,zaxpy_incx_0)
{
        blasint i;
        blasint N=4,incX=0,incY=1;
        double a[2]={0.25,0.5};
        double x1[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
        double y1[]={2.0,4.0,6.0,8.0,2.0,4.0,6.0,8.0};
        double x2[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
        double y2[]={0.75,5.25,4.75,9.25,0.75,5.25,4.75,9.25};

        //OpenBLAS
        BLASFUNC(zaxpy)(&N,a,x1,&incX,y1,&incY);

        for(i=0; i<2*N; i++){
		ASSERT_DBL_NEAR_TOL(x2[i], x1[i], DOUBLE_EPS);
                ASSERT_DBL_NEAR_TOL(y2[i], y1[i], DOUBLE_EPS);
        }
}

#endif

#ifdef BUILD_SINGLE
CTEST(axpy,saxpy_inc_0)
{
	blasint i;
	blasint N=8,incX=0,incY=0;
	float a=0.25;
	float x1[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	float y1[]={2.0,4.0,6.0,8.0,2.0,4.0,6.0,8.0};
	float x2[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	float y2[]={4.0,4.0,6.0,8.0,2.0,4.0,6.0,8.0};

	//OpenBLAS
	BLASFUNC(saxpy)(&N,&a,x1,&incX,y1,&incY);

	for(i=0; i<N; i++){
		ASSERT_DBL_NEAR_TOL(x2[i], x1[i], DOUBLE_EPS);
		ASSERT_DBL_NEAR_TOL(y2[i], y1[i], DOUBLE_EPS);
	}
}
#endif

#ifdef BUILD_COMPLEX
CTEST(axpy,caxpy_inc_0)
{
	blasint i;
	blasint N=4,incX=0,incY=0;
	float a[2]={0.25,0.5};
	float x1[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	float y1[]={2.0,4.0,6.0,8.0,2.0,4.0,6.0,8.0};
	float x2[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	float y2[]={-3.0,9.0,6.0,8.0,2.0,4.0,6.0,8.0};

	//OpenBLAS
	BLASFUNC(caxpy)(&N,a,x1,&incX,y1,&incY);

	for(i=0; i<2*N; i++){
		ASSERT_DBL_NEAR_TOL(x2[i], x1[i], DOUBLE_EPS);
		ASSERT_DBL_NEAR_TOL(y2[i], y1[i], DOUBLE_EPS);
	}
}

CTEST(axpy,caxpy_incx_0)
{
	blasint i;
	blasint N=4,incX=0,incY=1;
	float a[2]={0.25,0.5};
	float x1[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	float y1[]={2.0,4.0,6.0,8.0,2.0,4.0,6.0,8.0};
        double x2[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
        double y2[]={0.75,5.25,4.75,9.25,0.75,5.25,4.75,9.25};

	//OpenBLAS
	BLASFUNC(caxpy)(&N,a,x1,&incX,y1,&incY);

	for(i=0; i<2*N; i++){
		ASSERT_DBL_NEAR_TOL(x2[i], x1[i], DOUBLE_EPS);
		ASSERT_DBL_NEAR_TOL(y2[i], y1[i], DOUBLE_EPS);
	}
}
#endif



#ifdef BUILD_DOUBLE
CTEST(axpy,daxpy_zero_output_stride)
{
    const blasint lengths[] = {1, 2, 9, 17, 33};
    const blasint strides[] = {1, 2, -1, -2};
    const double seeds[] = {-2, 0, 2};
    const double alphas[] = {0, 1, -1, 0.5};
    int ni, si, yi, ai, k;
    for (ni = 0; ni < 5; ++ni)
        for (si = 0; si < 4; ++si)
            for (yi = 0; yi < 3; ++yi)
                for (ai = 0; ai < 4; ++ai) {
                    blasint n = lengths[ni], incx = strides[si], incy = 0;
                    double alpha = alphas[ai], x[65], y[3] = {99, 0, 99};
                    double expected = seeds[yi] + n * alpha;
                    y[1] = seeds[yi];
                    for (k = 0; k < 65; ++k) x[k] = 1;
                    BLASFUNC(daxpy)(&n, &alpha, x, &incx, y + 1, &incy);
                    ASSERT_DBL_NEAR_TOL(expected, y[1], DOUBLE_EPS);
                    ASSERT_DBL_NEAR_TOL(99, y[0], DOUBLE_EPS);
                    ASSERT_DBL_NEAR_TOL(99, y[2], DOUBLE_EPS);
                    for (k = 0; k < 65; ++k)
                        ASSERT_DBL_NEAR_TOL(1, x[k], DOUBLE_EPS);
                }
}
#endif

#ifdef BUILD_SINGLE
CTEST(axpy,saxpy_zero_output_stride)
{
    const blasint lengths[] = {1, 2, 9, 17, 33};
    const blasint strides[] = {1, 2, -1, -2};
    const float seeds[] = {-2, 0, 2};
    const float alphas[] = {0, 1, -1, 0.5};
    int ni, si, yi, ai, k;
    for (ni = 0; ni < 5; ++ni)
        for (si = 0; si < 4; ++si)
            for (yi = 0; yi < 3; ++yi)
                for (ai = 0; ai < 4; ++ai) {
                    blasint n = lengths[ni], incx = strides[si], incy = 0;
                    float alpha = alphas[ai], x[65], y[3] = {99, 0, 99};
                    float expected = seeds[yi] + n * alpha;
                    y[1] = seeds[yi];
                    for (k = 0; k < 65; ++k) x[k] = 1;
                    BLASFUNC(saxpy)(&n, &alpha, x, &incx, y + 1, &incy);
                    ASSERT_DBL_NEAR_TOL(expected, y[1], DOUBLE_EPS);
                    ASSERT_DBL_NEAR_TOL(99, y[0], DOUBLE_EPS);
                    ASSERT_DBL_NEAR_TOL(99, y[2], DOUBLE_EPS);
                    for (k = 0; k < 65; ++k)
                        ASSERT_DBL_NEAR_TOL(1, x[k], DOUBLE_EPS);
                }
}
#endif

#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
#ifdef BUILD_COMPLEX
CTEST(axpy,caxpy_zero_output_stride)
{
    const blasint lengths[] = {1, 2, 9, 17, 33};
    const blasint strides[] = {1, 2, -1, -2};
    int ni, si, k;
    for (ni = 0; ni < 5; ++ni)
        for (si = 0; si < 4; ++si) {
            blasint n = lengths[ni], incx = strides[si], incy = 0;
            float alpha[2] = {1.0f, 2.0f};
            float x[132], y[6] = {99.0f, 99.0f, 0.5f, -1.5f, 99.0f, 99.0f};
            for (k = 0; k < 66; ++k) {
                x[2 * k] = 3.0f;
                x[2 * k + 1] = 4.0f;
            }
            BLASFUNC(caxpy)(&n, alpha, x, &incx, y + 2, &incy);
            /* (1 + 2i) * (3 + 4i) = -5 + 10i */
            ASSERT_DBL_NEAR_TOL(0.5 - 5.0 * (double)n, (double)y[2], SINGLE_EPS);
            ASSERT_DBL_NEAR_TOL(-1.5 + 10.0 * (double)n, (double)y[3], SINGLE_EPS);
            ASSERT_DBL_NEAR_TOL(99.0, (double)y[0], SINGLE_EPS);
            ASSERT_DBL_NEAR_TOL(99.0, (double)y[1], SINGLE_EPS);
            ASSERT_DBL_NEAR_TOL(99.0, (double)y[4], SINGLE_EPS);
            ASSERT_DBL_NEAR_TOL(99.0, (double)y[5], SINGLE_EPS);
        }
}
#endif

#ifdef BUILD_COMPLEX16
CTEST(axpy,zaxpy_zero_output_stride)
{
    const blasint lengths[] = {1, 2, 9, 17, 33};
    const blasint strides[] = {1, 2, -1, -2};
    int ni, si, k;
    for (ni = 0; ni < 5; ++ni)
        for (si = 0; si < 4; ++si) {
            blasint n = lengths[ni], incx = strides[si], incy = 0;
            double alpha[2] = {1.0, 2.0};
            double x[132], y[6] = {99.0, 99.0, 0.5, -1.5, 99.0, 99.0};
            for (k = 0; k < 66; ++k) {
                x[2 * k] = 3.0;
                x[2 * k + 1] = 4.0;
            }
            BLASFUNC(zaxpy)(&n, alpha, x, &incx, y + 2, &incy);
            ASSERT_DBL_NEAR_TOL(0.5 - 5.0 * (double)n, y[2], DOUBLE_EPS);
            ASSERT_DBL_NEAR_TOL(-1.5 + 10.0 * (double)n, y[3], DOUBLE_EPS);
            ASSERT_DBL_NEAR_TOL(99.0, y[0], DOUBLE_EPS);
            ASSERT_DBL_NEAR_TOL(99.0, y[1], DOUBLE_EPS);
            ASSERT_DBL_NEAR_TOL(99.0, y[4], DOUBLE_EPS);
            ASSERT_DBL_NEAR_TOL(99.0, y[5], DOUBLE_EPS);
        }
}
#endif
#endif
