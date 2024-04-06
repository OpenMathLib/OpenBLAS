#include "openblas_utest.h"
#include <cblas.h>
#ifdef BUILD_COMPLEX16

#ifndef NAN
#define NAN 0.0/0.0
#endif
#ifndef INFINITY
#define INFINITY 1.0/0.0
#endif

CTEST(zscal, i_nan)
{
    blasint N=9;
    blasint incX=1;
    double i[] = {0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1 };
    double nan[] = {NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0};
    BLASFUNC(zscal)(&N, i, nan, &incX);
    ASSERT_TRUE(isnan(nan[0]));
    ASSERT_TRUE(isnan(nan[1]));
    ASSERT_TRUE(isnan(nan[16]));
    ASSERT_TRUE(isnan(nan[17]));
}

CTEST(zscal, i_nan_inc_2)
{
    blasint N=9;
    blasint incX=1;
    double i[] = {0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1 };
    double nan[] = {NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0,
                    NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0};
    BLASFUNC(zscal)(&N, i, nan, &incX);
    ASSERT_TRUE(isnan(nan[0]));
    ASSERT_TRUE(isnan(nan[1]));
    ASSERT_TRUE(isnan(nan[16]));
    ASSERT_TRUE(isnan(nan[17]));
}

CTEST(zscal, nan_i)
{
    blasint N=9;
    blasint incX=1;
    double i[] = {0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1 };
    double nan[] = {NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0};
    BLASFUNC(zscal)(&N, nan, i, &incX);
    ASSERT_TRUE(isnan(i[0]));
    ASSERT_TRUE(isnan(i[1]));
    ASSERT_TRUE(isnan(i[16]));
    ASSERT_TRUE(isnan(i[17]));
}

CTEST(zscal, nan_i_inc_2)
{
    blasint N=9;
    blasint incX=1;
    double i[] = {0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1,
                  0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1 };
    double nan[] = {NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0};
    BLASFUNC(zscal)(&N, nan, i, &incX);
    ASSERT_TRUE(isnan(i[0]));
    ASSERT_TRUE(isnan(i[1]));
    ASSERT_TRUE(isnan(i[16]));
    ASSERT_TRUE(isnan(i[17]));
}

CTEST(zscal, i_inf)
{
    blasint N=9;
    blasint incX=1;
    double i[] = {0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1 };
    double inf[] = {INFINITY, 0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0};
    BLASFUNC(zscal)(&N, i, inf, &incX);
    ASSERT_TRUE(isnan(inf[0]));
    ASSERT_TRUE(isinf(inf[1]));
    ASSERT_TRUE(isnan(inf[16]));
    ASSERT_TRUE(isinf(inf[17]));
}

CTEST(zscal, i_inf_inc_2)
{
    blasint N=9;
    blasint incX=2;
    double i[] = {0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1 };
    double inf[] = {INFINITY, 0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0,
                    INFINITY, 0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0};
    BLASFUNC(zscal)(&N, i, inf, &incX);
    ASSERT_TRUE(isnan(inf[0]));
    ASSERT_TRUE(isinf(inf[1]));
    ASSERT_TRUE(isnan(inf[16]));
    ASSERT_TRUE(isinf(inf[17]));
}

CTEST(zscal, inf_i)
{
    blasint N=9;
    blasint incX=1;
    double i[] = {0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1 };
    double inf[] = {INFINITY, 0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0};
    BLASFUNC(zscal)(&N, inf, i, &incX);
    ASSERT_TRUE(isnan(i[0]));
    ASSERT_TRUE(isinf(i[1]));
    ASSERT_TRUE(isnan(i[16]));
    ASSERT_TRUE(isinf(i[17]));
}

CTEST(zscal, inf_i_inc_2)
{
    blasint N=9;
    blasint incX=2;
    double i[] = {0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1,
                  0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1 };
    double inf[] = {INFINITY, 0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0};
    BLASFUNC(zscal)(&N, inf, i, &incX);
    ASSERT_TRUE(isnan(i[0]));
    ASSERT_TRUE(isinf(i[1]));
    ASSERT_TRUE(isnan(i[16]));
    ASSERT_TRUE(isinf(i[17]));
}

#endif
