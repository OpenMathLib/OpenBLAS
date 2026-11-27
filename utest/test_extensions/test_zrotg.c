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

#ifdef BUILD_COMPLEX16

/**
 * Fortran API specific test
 * Test zrotg by comparing it against pre-calculated values
 */
CTEST(zrotg, zero_a)
{
    double sa[2] = {0.0, 0.0};
    double sb[2] = {1.0, 1.0};
    double ss[2];
    double sc;
    BLASFUNC(zrotg)(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.0, sc, DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.70710678118655, ss[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(-0.70710678118655, ss[1], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.41421356237310, sa[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0, sa[1], DOUBLE_EPS);
}

/**
 * Fortran API specific tests
 * Test zrotg by comparing it against pre-calculated values
 */
CTEST(zrotg, zero_b)
{
    double sa[2] = {1.0, 1.0};
    double sb[2] = {0.0, 0.0};
    double ss[2];
    double sc;
    BLASFUNC(zrotg)(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(1.0, sc, DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0, ss[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0, ss[1], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.0, sa[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.0, sa[1], DOUBLE_EPS);
}

/**
 * Fortran API specific test
 * Test zrotg by comparing it against pre-calculated values
 */
CTEST(zrotg, zero_real)
{
    double sa[2] = {0.0, 1.0};
    double sb[2] = {0.0, 1.0};
    double ss[2];
    double sc;
    BLASFUNC(zrotg)(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.70710678118654, sc, DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.70710678118654, ss[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0, ss[1], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0, sa[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.41421356237309, sa[1], DOUBLE_EPS);
}

/**
 * Fortran API specific test
 * Test zrotg by comparing it against pre-calculated values
 */
CTEST(zrotg, positive_real_positive_img)
{
    double sa[2] = {3.0, 4.0};
    double sb[2] = {4.0, 6.0};
    double ss[2];
    double sc;
    BLASFUNC(zrotg)(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.56980288229818, sc, DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.82051615050939, ss[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(-0.04558423058385, ss[1], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(5.26497863243527, sa[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(7.01997150991369, sa[1], DOUBLE_EPS);
}

/**
 * Test the scaled path with a nonzero imaginary part of b
 */
CTEST(zrotg, scaled_complex_b)
{
    double scale = ldexp(1.0, 512);
    double expected = 1.0 / sqrt(3.0);
    double sa[2] = {scale, 0.0};
    double sb[2] = {scale, scale};
    double ss[2];
    double sc;
    BLASFUNC(zrotg)(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(expected, sc, DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(expected, ss[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(-expected, ss[1], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(sqrt(3.0), sa[0] / scale, DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0, sa[1] / scale, DOUBLE_EPS);
}

/**
 * Test the scaled path when a is much smaller than b
 */
CTEST(zrotg, scaled_separated_magnitudes)
{
    double scale_a = ldexp(1.0, 300);
    double scale_b = ldexp(1.0, 600);
    double sa[2] = {scale_a, scale_a};
    double sb[2] = {scale_b, scale_b};
    double ss[2];
    double sc;
    BLASFUNC(zrotg)(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(1.0, sc / (scale_a / scale_b), DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.0, ss[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0, ss[1], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.0, sa[0] / scale_b, DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.0, sa[1] / scale_b, DOUBLE_EPS);
}

/**
 * Fortran API specific test
 * Test zrotg by comparing it against pre-calculated values
 */
CTEST(zrotg, negative_real_positive_img)
{
    double sa[2] = {-3.0, 4.0};
    double sb[2] = {-4.0, 6.0};
    double ss[2];
    double sc;
    BLASFUNC(zrotg)(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.56980288229818, sc, DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.82051615050939, ss[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.04558423058385, ss[1], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(-5.26497863243527, sa[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(7.01997150991369, sa[1], DOUBLE_EPS);
}

/**
 * Fortran API specific test
 * Test zrotg by comparing it against pre-calculated values
 */
CTEST(zrotg, positive_real_negative_img)
{
    double sa[2] = {3.0, -4.0};
    double sb[2] = {4.0, -6.0};
    double ss[2];
    double sc;
    BLASFUNC(zrotg)(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.56980288229818, sc, DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.82051615050939, ss[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.04558423058385, ss[1], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(5.26497863243527, sa[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(-7.01997150991369, sa[1], DOUBLE_EPS);
}

/**
 * Fortran API specific test
 * Test zrotg by comparing it against pre-calculated values
 */
CTEST(zrotg, negative_real_negative_img)
{
    double sa[2] = {-3.0, -4.0};
    double sb[2] = {-4.0, -6.0};
    double ss[2];
    double sc;
    BLASFUNC(zrotg)(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.56980288229818, sc, DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.82051615050939, ss[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(-0.04558423058385, ss[1], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(-5.26497863243527, sa[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(-7.01997150991369, sa[1], DOUBLE_EPS);
}

/**
 * Fortran API specific test
 * Test zrotg with a complex a much smaller than b, where the cosine
 * is computed in the small-cosine branch
 */
CTEST(zrotg, unbalanced_complex_a)
{
    double sa[2] = {ldexp(1.0, -200), ldexp(1.0, -200)};
    double sb[2] = {ldexp(1.0, 400), 0.0};
    double ss[2];
    double sc;
    BLASFUNC(zrotg)(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(1.41421356237310, ldexp(sc, 600), DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.70710678118655, ss[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.70710678118655, ss[1], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.41421356237310, ldexp(sa[0], -399), DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.41421356237310, ldexp(sa[1], -399), DOUBLE_EPS);
}

#ifndef NO_CBLAS
/**
 * C API specific test
 * Test zrotg by comparing it against pre-calculated values
 */
CTEST(zrotg, c_api_zero_a)
{
    double sa[2] = {0.0, 0.0};
    double sb[2] = {1.0, 1.0};
    double ss[2];
    double sc;
    cblas_zrotg(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.0, sc, DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.70710678118655, ss[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(-0.70710678118655, ss[1], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.41421356237310, sa[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0, sa[1], DOUBLE_EPS);
}

/**
 * C API specific test
 * Test zrotg by comparing it against pre-calculated values
 */
CTEST(zrotg, c_api_zero_b)
{
    double sa[2] = {1.0, 1.0};
    double sb[2] = {0.0, 0.0};
    double ss[2];
    double sc;
    cblas_zrotg(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(1.0, sc, DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0, ss[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0, ss[1], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.0, sa[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.0, sa[1], DOUBLE_EPS);
}

/**
 * C API specific test
 * Test zrotg by comparing it against pre-calculated values
 */
CTEST(zrotg, c_api_zero_real)
{
    double sa[2] = {0.0, 1.0};
    double sb[2] = {0.0, 1.0};
    double ss[2];
    double sc;
    cblas_zrotg(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.70710678118654, sc, DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.70710678118654, ss[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0, ss[1], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0, sa[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.41421356237309, sa[1], DOUBLE_EPS);
}

/**
 * C API specific test
 * Test zrotg by comparing it against pre-calculated values
 */
CTEST(zrotg, c_api_positive_real_positive_img)
{
    double sa[2] = {3.0, 4.0};
    double sb[2] = {4.0, 6.0};
    double ss[2];
    double sc;
    cblas_zrotg(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.56980288229818, sc, DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.82051615050939, ss[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(-0.04558423058385, ss[1], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(5.26497863243527, sa[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(7.01997150991369, sa[1], DOUBLE_EPS);
}

/**
 * Test the scaled path with a nonzero imaginary part of b
 */
CTEST(zrotg, c_api_scaled_complex_b)
{
    double scale = ldexp(1.0, 512);
    double expected = 1.0 / sqrt(3.0);
    double sa[2] = {scale, 0.0};
    double sb[2] = {scale, scale};
    double ss[2];
    double sc;
    cblas_zrotg(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(expected, sc, DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(expected, ss[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(-expected, ss[1], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(sqrt(3.0), sa[0] / scale, DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0, sa[1] / scale, DOUBLE_EPS);
}

/**
 * Test the scaled path when a is much smaller than b
 */
CTEST(zrotg, c_api_scaled_separated_magnitudes)
{
    double scale_a = ldexp(1.0, 300);
    double scale_b = ldexp(1.0, 600);
    double sa[2] = {scale_a, scale_a};
    double sb[2] = {scale_b, scale_b};
    double ss[2];
    double sc;
    cblas_zrotg(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(1.0, sc / (scale_a / scale_b), DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.0, ss[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.0, ss[1], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.0, sa[0] / scale_b, DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.0, sa[1] / scale_b, DOUBLE_EPS);
}

/**
 * C API specific test
 * Test zrotg by comparing it against pre-calculated values
 */
CTEST(zrotg, c_api_negative_real_positive_img)
{
    double sa[2] = {-3.0, 4.0};
    double sb[2] = {-4.0, 6.0};
    double ss[2];
    double sc;
    cblas_zrotg(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.56980288229818, sc, DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.82051615050939, ss[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.04558423058385, ss[1], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(-5.26497863243527, sa[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(7.01997150991369, sa[1], DOUBLE_EPS);
}

/**
 * C API specific test
 * Test zrotg by comparing it against pre-calculated values
 */
CTEST(zrotg, c_api_positive_real_negative_img)
{
    double sa[2] = {3.0, -4.0};
    double sb[2] = {4.0, -6.0};
    double ss[2];
    double sc;
    cblas_zrotg(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.56980288229818, sc, DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.82051615050939, ss[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.04558423058385, ss[1], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(5.26497863243527, sa[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(-7.01997150991369, sa[1], DOUBLE_EPS);
}

/**
 * C API specific test
 * Test zrotg by comparing it against pre-calculated values
 */
CTEST(zrotg, c_api_negative_real_negative_img)
{
    double sa[2] = {-3.0, -4.0};
    double sb[2] = {-4.0, -6.0};
    double ss[2];
    double sc;
    cblas_zrotg(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(0.56980288229818, sc, DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.82051615050939, ss[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(-0.04558423058385, ss[1], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(-5.26497863243527, sa[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(-7.01997150991369, sa[1], DOUBLE_EPS);
}

/**
 * C API specific test
 * Test zrotg with a complex a much smaller than b, where the cosine
 * is computed in the small-cosine branch
 */
CTEST(zrotg, c_api_unbalanced_complex_a)
{
    double sa[2] = {ldexp(1.0, -200), ldexp(1.0, -200)};
    double sb[2] = {ldexp(1.0, 400), 0.0};
    double ss[2];
    double sc;
    cblas_zrotg(sa, sb, &sc, ss);
    ASSERT_DBL_NEAR_TOL(1.41421356237310, ldexp(sc, 600), DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.70710678118655, ss[0], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(0.70710678118655, ss[1], DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.41421356237310, ldexp(sa[0], -399), DOUBLE_EPS);
    ASSERT_DBL_NEAR_TOL(1.41421356237310, ldexp(sa[1], -399), DOUBLE_EPS);
}
#endif
#endif
