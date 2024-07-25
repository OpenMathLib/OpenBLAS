#include "openblas_utest.h"
#include <cblas.h>

#ifndef NAN
#define NAN 0.0/0.0
#endif
#ifndef INFINITY
#define INFINITY 1.0/0.0
#endif

#ifdef BUILD_SINGLE

CTEST(sgemv, 0_nan_inf)
{
    blasint N = 17;
    blasint incX = 1;
    blasint incY = 1;
    float alpha = 0.0;
    float beta = 0.0;
    char  trans = 'N';
    float A[N * N];
    float X[N];
    float Y[N];

    memset(A, 0, sizeof(A));
    memset(X, 0, sizeof(X));
    for (int i = 0; i < (N - 1); i += 2)
    {
        Y[i]     = NAN;
        Y[i + 1] = INFINITY;
    }
    Y[N - 1] = NAN;
    BLASFUNC(sgemv)(&trans, &N, &N, &alpha, A, &N, X, &incX, &beta, Y, &incY);
    for (int i = 0; i < N; i ++)
        ASSERT_TRUE(Y[i] == 0.0);
}

CTEST(sgemv, 0_nan_inf_incy_2)
{
    blasint N  = 17;
    blasint Ny = 33;
    blasint incX = 1;
    blasint incY = 2;
    float alpha = 0.0;
    float beta = 0.0;
    char  trans = 'N';
    float A[N * N];
    float X[N];
    float Y[Ny];
    float *ay = Y;

    memset(A, 0, sizeof(A));
    memset(X, 0, sizeof(X));
    memset(Y, 0, sizeof(Y));
    for (int i = 0; i < (N - 1); i += 2)
    {
        ay[0]   = NAN;
        ay     += 2;
        ay[0]   = INFINITY;
        ay     += 2;
    }
    Y[Ny - 1] = NAN;
    BLASFUNC(sgemv)(&trans, &N, &N, &alpha, A, &N, X, &incX, &beta, Y, &incY);
    for (int i = 0; i < Ny; i ++)
        ASSERT_TRUE(Y[i] == 0.0);
}

#endif

#ifdef BUILD_DOUBLE
CTEST(dgemv, 0_nan_inf)
{
    blasint N = 17;
    blasint incX = 1;
    blasint incY = 1;
    double alpha = 0.0;
    double beta = 0.0;
    char  trans = 'N';
    double A[N * N];
    double X[N];
    double Y[N];

    memset(A, 0, sizeof(A));
    memset(X, 0, sizeof(X));
    for (int i = 0; i < (N - 1); i += 2)
    {
        Y[i]     = NAN;
        Y[i + 1] = INFINITY;
    }
    Y[N - 1] = NAN;
    BLASFUNC(dgemv)(&trans, &N, &N, &alpha, A, &N, X, &incX, &beta, Y, &incY);
    for (int i = 0; i < N; i ++)
        ASSERT_TRUE(Y[i] == 0.0);
}

CTEST(dgemv, 0_nan_inf_incy_2)
{
    blasint N  = 17;
    blasint Ny = 33;
    blasint incX = 1;
    blasint incY = 2;
    double alpha = 0.0;
    double beta = 0.0;
    char  trans = 'N';
    double A[N * N];
    double X[N];
    double Y[Ny];
    double *ay = Y;

    memset(A, 0, sizeof(A));
    memset(X, 0, sizeof(X));
    memset(Y, 0, sizeof(Y));
    for (int i = 0; i < (N - 1); i += 2)
    {
        ay[0]   = NAN;
        ay     += 2;
        ay[0]   = INFINITY;
        ay     += 2;
    }
    Y[Ny - 1] = NAN;
    BLASFUNC(dgemv)(&trans, &N, &N, &alpha, A, &N, X, &incX, &beta, Y, &incY);
    for (int i = 0; i < Ny; i ++)
        ASSERT_TRUE(Y[i] == 0.0);
}

#endif
