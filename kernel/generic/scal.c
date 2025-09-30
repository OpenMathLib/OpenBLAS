/***************************************************************************
Copyright (c) 2013, 2025 The OpenBLAS Project
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
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#include "common.h"

int CNAME(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT da, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *dummy, BLASLONG dummy2)
{
    BLASLONG i = 0, j = 0;
#if defined(BFLOAT16)
    float x_float, da_float;
    SBF16TOS_K(1, &da, 1, &da_float, 1);
#else
    float x_float;
    float da_float = da;
#endif

    if ((n <= 0) || (inc_x <= 0))
        return (0);

    if (dummy2 == 0)
    {
        while (j < n)
        {

            if (da_float == 0.0)
                x_float = 0.0;
            else
            {
#if defined(BFLOAT16)
                SBF16TOS_K(1, &x[i], 1, &x_float, 1);
#else
                float x_float = x[i];
#endif
                x_float = da_float * x_float;
            }

#if defined(BFLOAT16)
            SBSTOBF16_K(1, &x_float, 1, &x[i], 1);
#else
            x[i] = x_float;
#endif

            i += inc_x;
            j++;
        }
    }
    else
    {

        while (j < n)
        {
#if defined(BFLOAT16)
            SBF16TOS_K(1, &x[i], 1, &x_float, 1);
#else
            float x_float = x[i];
#endif
            if (da == 0.0)
                if (!isnan(x_float) && !isinf(x_float))
                {
                    x_float = 0.0;
                }
                else
                {
                    x_float = NAN;
                }
            else
            {
                x_float = da_float * x_float;
            }

#if defined(BFLOAT16)
            SBSTOBF16_K(1, &x_float, 1, &x[i], 1);
#else
            x[i] = x_float;
#endif

            i += inc_x;
            j++;
        }
    }
    return 0;
}
