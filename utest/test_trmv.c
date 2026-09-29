/*****************************************************************************
Copyright (c) 2026, The OpenBLAS Project
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
   3. Neither the name of OpenBLAS nor the names of its contributors may
      be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

**********************************************************************************/

#include "openblas_utest.h"
#include <cblas.h>

#if defined(BUILD_SINGLE) && !defined(NO_CBLAS)

/*
 * Regression test for equal unit strides in the ARM64 ASIMD dot kernel
 * under LTO.
 */
CTEST(strmv, small_rowmajor_lower_notrans_unit_stride)
{
    const float lower[9] = {
        2.0f, 0.0f, 0.0f,
        3.0f, 4.0f, 0.0f,
        5.0f, 6.0f, 7.0f
    };
    const float expected[3] = {22.0f, 85.0f, 252.0f};
    float vector[3] = {11.0f, 13.0f, 17.0f};
    int i;

    cblas_strmv(CblasRowMajor, CblasLower, CblasNoTrans,
                CblasNonUnit, 3, lower, 3, vector, 1);

    for (i = 0; i < 3; i++)
        ASSERT_DBL_NEAR_TOL(expected[i], vector[i], SINGLE_EPS);
}

#endif
