/*****************************************************************************
Copyright (c) 2023, The OpenBLAS Project
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

#include "utest/openblas_utest.h"
#include <cblas.h>

#ifdef BUILD_COMPLEX

/**
 * Fortran API specific test
 * Test crotg by comparing it against pre-calculated values
 */
CTEST(crotg, zero_a)
{
    float sa[2] = {0.0f, 0.0f};
    float sb[2] = {1.0f, 1.0f};
    float ss[2];
    float sc;
    BLASFUNC(crotg)(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.0f, sc, SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.70711f, ss[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(-0.70711f, ss[1], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.41421f, sa[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0f, sa[1], SINGLE_EPS);
}

/**
 * Fortran API specific test
 * Test crotg by comparing it against pre-calculated values
 */
CTEST(crotg, zero_b)
{
    float sa[2] = {1.0f, 1.0f};
    float sb[2] = {0.0f, 0.0f};
    float ss[2];
    float sc;
    BLASFUNC(crotg)(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(1.0f, sc, SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0f, ss[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0f, ss[1], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.0f, sa[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.0f, sa[1], SINGLE_EPS);
}

/**
 * Fortran API specific test
 * Test crotg by comparing it against pre-calculated values
 */
CTEST(crotg, zero_real)
{
    float sa[2] = {0.0f, 1.0f};
    float sb[2] = {0.0f, 1.0f};
    float ss[2];
    float sc;
    BLASFUNC(crotg)(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.70711f, sc, SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.70711f, ss[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0f, ss[1], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0f, sa[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.41421f, sa[1], SINGLE_EPS);
}

/**
 * Fortran API specific test
 * Test crotg by comparing it against pre-calculated values
 */
CTEST(crotg, positive_real_positive_img)
{
    float sa[2] = {3.0f, 4.0f};
    float sb[2] = {4.0f, 6.0f};
    float ss[2];
    float sc;
    BLASFUNC(crotg)(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.5698f, sc, SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.82052f, ss[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(-0.04558f, ss[1], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(5.26498f, sa[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(7.01997f, sa[1], SINGLE_EPS);
}

/**
 * Test the scaled path with a nonzero imaginary part of b
 */
CTEST(crotg, scaled_complex_b)
{
    float scale = (float)ldexp(1.0, 70);
    float expected = (float)(1.0 / sqrt(3.0));
    float sa[2] = {scale, 0.0f};
    float sb[2] = {scale, scale};
    float ss[2];
    float sc;
    BLASFUNC(crotg)(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(expected, sc, SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(expected, ss[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(-expected, ss[1], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(sqrt(3.0), sa[0] / scale, SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0f, sa[1] / scale, SINGLE_EPS);
}

/**
 * Test the scaled path when a is much smaller than b
 */
CTEST(crotg, scaled_separated_magnitudes)
{
    float scale_a = (float)ldexp(1.0, 30);
    float scale_b = (float)ldexp(1.0, 80);
    float sa[2] = {scale_a, scale_a};
    float sb[2] = {scale_b, scale_b};
    float ss[2];
    float sc;
    BLASFUNC(crotg)(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(1.0f, sc / (scale_a / scale_b), SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.0f, ss[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0f, ss[1], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.0f, sa[0] / scale_b, SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.0f, sa[1] / scale_b, SINGLE_EPS);
}

/**
 * Fortran API specific test
 * Test crotg by comparing it against pre-calculated values
 */
CTEST(crotg, negative_real_positive_img)
{
    float sa[2] = {-3.0f, 4.0f};
    float sb[2] = {-4.0f, 6.0f};
    float ss[2];
    float sc;
    BLASFUNC(crotg)(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.5698f, sc, SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.82052f, ss[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.04558f, ss[1], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(-5.26498f, sa[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(7.01997f, sa[1], SINGLE_EPS);
}

/**
 * Fortran API specific test
 * Test crotg by comparing it against pre-calculated values
 */
CTEST(crotg, positive_real_negative_img)
{
    float sa[2] = {3.0f, -4.0f};
    float sb[2] = {4.0f, -6.0f};
    float ss[2];
    float sc;
    BLASFUNC(crotg)(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.5698f, sc, SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.82052f, ss[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.04558f, ss[1], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(5.26498f, sa[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(-7.01997f, sa[1], SINGLE_EPS);
}

/**
 * Fortran API specific test
 * Test crotg by comparing it against pre-calculated values
 */
CTEST(crotg, negative_real_negative_img)
{
    float sa[2] = {-3.0f, -4.0f};
    float sb[2] = {-4.0f, -6.0f};
    float ss[2];
    float sc;
    BLASFUNC(crotg)(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.5698f, sc, SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.82052f, ss[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(-0.04558f, ss[1], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(-5.26498f, sa[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(-7.01997f, sa[1], SINGLE_EPS);
}
#ifndef NO_CBLAS
/**
 * C API specific test
 * Test crotg by comparing it against pre-calculated values
 */
CTEST(crotg, c_api_zero_a)
{
    float sa[2] = {0.0f, 0.0f};
    float sb[2] = {1.0f, 1.0f};
    float ss[2];
    float sc;
    cblas_crotg(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.0f, sc, SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.70711f, ss[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(-0.70711f, ss[1], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.41421f, sa[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0f, sa[1], SINGLE_EPS);
}

/**
 * C API specific test
 * Test crotg by comparing it against pre-calculated values
 */
CTEST(crotg, c_api_zero_b)
{
    float sa[2] = {1.0f, 1.0f};
    float sb[2] = {0.0f, 0.0f};
    float ss[2];
    float sc;
    cblas_crotg(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(1.0f, sc, SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0f, ss[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0f, ss[1], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.0f, sa[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.0f, sa[1], SINGLE_EPS);
}

/**
 * C API specific test
 * Test crotg by comparing it against pre-calculated values
 */
CTEST(crotg, c_api_zero_real)
{
    float sa[2] = {0.0f, 1.0f};
    float sb[2] = {0.0f, 1.0f};
    float ss[2];
    float sc;
    cblas_crotg(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.70711f, sc, SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.70711f, ss[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0f, ss[1], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0f, sa[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.41421f, sa[1], SINGLE_EPS);
}

/**
 * C API specific test
 * Test crotg by comparing it against pre-calculated values
 */
CTEST(crotg, c_api_positive_real_positive_img)
{
    float sa[2] = {3.0f, 4.0f};
    float sb[2] = {4.0f, 6.0f};
    float ss[2];
    float sc;
    cblas_crotg(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.5698f, sc, SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.82052f, ss[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(-0.04558f, ss[1], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(5.26498f, sa[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(7.01997f, sa[1], SINGLE_EPS);
}

/**
 * Test the scaled path with a nonzero imaginary part of b
 */
CTEST(crotg, c_api_scaled_complex_b)
{
    float scale = (float)ldexp(1.0, 70);
    float expected = (float)(1.0 / sqrt(3.0));
    float sa[2] = {scale, 0.0f};
    float sb[2] = {scale, scale};
    float ss[2];
    float sc;
    cblas_crotg(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(expected, sc, SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(expected, ss[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(-expected, ss[1], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(sqrt(3.0), sa[0] / scale, SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0f, sa[1] / scale, SINGLE_EPS);
}

/**
 * Test the scaled path when a is much smaller than b
 */
CTEST(crotg, c_api_scaled_separated_magnitudes)
{
    float scale_a = (float)ldexp(1.0, 30);
    float scale_b = (float)ldexp(1.0, 80);
    float sa[2] = {scale_a, scale_a};
    float sb[2] = {scale_b, scale_b};
    float ss[2];
    float sc;
    cblas_crotg(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(1.0f, sc / (scale_a / scale_b), SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.0f, ss[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0f, ss[1], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.0f, sa[0] / scale_b, SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.0f, sa[1] / scale_b, SINGLE_EPS);
}

/**
 * C API specific test
 * Test crotg by comparing it against pre-calculated values
 */
CTEST(crotg, c_api_negative_real_positive_img)
{
    float sa[2] = {-3.0f, 4.0f};
    float sb[2] = {-4.0f, 6.0f};
    float ss[2];
    float sc;
    cblas_crotg(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.5698f, sc, SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.82052f, ss[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.04558f, ss[1], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(-5.26498f, sa[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(7.01997f, sa[1], SINGLE_EPS);
}

/**
 * C API specific test
 * Test crotg by comparing it against pre-calculated values
 */
CTEST(crotg, c_api_positive_real_negative_img)
{
    float sa[2] = {3.0f, -4.0f};
    float sb[2] = {4.0f, -6.0f};
    float ss[2];
    float sc;
    cblas_crotg(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.5698f, sc, SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.82052f, ss[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.04558f, ss[1], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(5.26498f, sa[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(-7.01997f, sa[1], SINGLE_EPS);
}

/**
 * C API specific test
 * Test crotg by comparing it against pre-calculated values
 */
CTEST(crotg, c_api_negative_real_negative_img)
{
    float sa[2] = {-3.0f, -4.0f};
    float sb[2] = {-4.0f, -6.0f};
    float ss[2];
    float sc;
    cblas_crotg(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.5698f, sc, SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.82052f, ss[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(-0.04558f, ss[1], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(-5.26498f, sa[0], SINGLE_EPS);
    ASSERT_DBL_NEAR_TOL(-7.01997f, sa[1], SINGLE_EPS);
}
#endif
#endif
