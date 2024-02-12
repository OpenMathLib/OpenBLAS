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
    double i[] = {0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1 };
    double nan[] = {NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0};
    cblas_zscal(9, i, &nan, 1);
    ASSERT_TRUE(isnan(nan[0]));
    ASSERT_TRUE(isnan(nan[1]));
    ASSERT_TRUE(isnan(nan[16]));
    ASSERT_TRUE(isnan(nan[17]));
}

CTEST(zscal, i_nan_inc_2)
{
    double i[] = {0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1 };
    double nan[] = {NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0,
                    NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0};
    cblas_zscal(9, i, &nan, 2);
    ASSERT_TRUE(isnan(nan[0]));
    ASSERT_TRUE(isnan(nan[1]));
    ASSERT_TRUE(isnan(nan[16]));
    ASSERT_TRUE(isnan(nan[17]));
}

CTEST(zscal, nan_i)
{
    double i[] = {0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1 };
    double nan[] = {NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0};
    cblas_zscal(9, &nan, &i, 1);
    ASSERT_TRUE(isnan(i[0]));
    ASSERT_TRUE(isnan(i[1]));
    ASSERT_TRUE(isnan(i[16]));
    ASSERT_TRUE(isnan(i[17]));
}

CTEST(zscal, nan_i_inc_2)
{
    double i[] = {0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1,
                  0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1 };
    double nan[] = {NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0, NAN,0};
    cblas_zscal(9, &nan, &i, 2);
    ASSERT_TRUE(isnan(i[0]));
    ASSERT_TRUE(isnan(i[1]));
    ASSERT_TRUE(isnan(i[16]));
    ASSERT_TRUE(isnan(i[17]));
}

CTEST(zscal, i_inf)
{
    double i[] = {0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1 };
    double inf[] = {INFINITY, 0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0};
    cblas_zscal(9, i, &inf, 1);
    ASSERT_TRUE(isnan(inf[0]));
    ASSERT_TRUE(isinf(inf[1]));
    ASSERT_TRUE(isnan(inf[16]));
    ASSERT_TRUE(isinf(inf[17]));
}

CTEST(zscal, i_inf_inc_2)
{
    double i[] = {0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1 };
    double inf[] = {INFINITY, 0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0,
                    INFINITY, 0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0};
    cblas_zscal(9, i, &inf, 2);
    ASSERT_TRUE(isnan(inf[0]));
    ASSERT_TRUE(isinf(inf[1]));
    ASSERT_TRUE(isnan(inf[16]));
    ASSERT_TRUE(isinf(inf[17]));
}

CTEST(zscal, inf_i)
{
    double i[] = {0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1 };
    double inf[] = {INFINITY, 0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0};
    cblas_zscal(9, &inf, &i, 1);
    ASSERT_TRUE(isnan(i[0]));
    ASSERT_TRUE(isinf(i[1]));
    ASSERT_TRUE(isnan(i[16]));
    ASSERT_TRUE(isinf(i[17]));
}

CTEST(zscal, inf_i_inc_2)
{
    double i[] = {0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1,
                  0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1, 0,1 };
    double inf[] = {INFINITY, 0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0, INFINITY,0};
    cblas_zscal(9, &inf, &i, 2);
    ASSERT_TRUE(isnan(i[0]));
    ASSERT_TRUE(isinf(i[1]));
    ASSERT_TRUE(isnan(i[16]));
    ASSERT_TRUE(isinf(i[17]));
}

#endif
