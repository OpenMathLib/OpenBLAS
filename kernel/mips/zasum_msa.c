/*******************************************************************************
Copyright (c) 2016, The OpenBLAS Project
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
*******************************************************************************/

#include "common.h"
#include <math.h>
#include "macros_msa.h"

#define AND_VEC_D(in)   ((v2f64) ((v2i64) in & and_vec))

#define PROCESS_ZD(inc_val)                           \
    if (n > 8)                                        \
    {                                                 \
        n -= 8;                                       \
                                                      \
        LD_DP8_INC(x, inc_val, src0, src1, src2,      \
                   src3, src4, src5, src6, src7);     \
                                                      \
        sum_abs0 = AND_VEC_D(src0);                   \
        sum_abs1 = AND_VEC_D(src1);                   \
        sum_abs2 = AND_VEC_D(src2);                   \
        sum_abs3 = AND_VEC_D(src3);                   \
        sum_abs0 += AND_VEC_D(src4);                  \
        sum_abs1 += AND_VEC_D(src5);                  \
        sum_abs2 += AND_VEC_D(src6);                  \
        sum_abs3 += AND_VEC_D(src7);                  \
    }                                                 \
    else                                              \
    {                                                 \
        sum_abs0 = zero_v;                            \
        sum_abs1 = zero_v;                            \
        sum_abs2 = zero_v;                            \
        sum_abs3 = zero_v;                            \
    }                                                 \
                                                      \
    for (i = (n >> 3); i--;)                          \
    {                                                 \
        LD_DP8_INC(x, inc_val, src0, src1, src2,      \
                   src3, src4, src5, src6, src7);     \
                                                      \
        sum_abs0 += AND_VEC_D(src0);                  \
        sum_abs1 += AND_VEC_D(src1);                  \
        sum_abs2 += AND_VEC_D(src2);                  \
        sum_abs3 += AND_VEC_D(src3);                  \
        sum_abs0 += AND_VEC_D(src4);                  \
        sum_abs1 += AND_VEC_D(src5);                  \
        sum_abs2 += AND_VEC_D(src6);                  \
        sum_abs3 += AND_VEC_D(src7);                  \
    }                                                 \
                                                      \
    if (n & 7)                                        \
    {                                                 \
        if ((n & 4) && (n & 2) && (n & 1))            \
        {                                             \
            LD_DP7_INC(x, inc_val, src0, src1, src2,  \
                       src3, src4, src5, src6);       \
                                                      \
            sum_abs0 += AND_VEC_D(src0);              \
            sum_abs1 += AND_VEC_D(src1);              \
            sum_abs2 += AND_VEC_D(src2);              \
            sum_abs3 += AND_VEC_D(src3);              \
            sum_abs0 += AND_VEC_D(src4);              \
            sum_abs1 += AND_VEC_D(src5);              \
            sum_abs2 += AND_VEC_D(src6);              \
        }                                             \
        else if ((n & 4) && (n & 2))                  \
        {                                             \
            LD_DP6_INC(x, inc_val, src0, src1, src2,  \
                       src3, src4, src5);             \
                                                      \
            sum_abs0 += AND_VEC_D(src0);              \
            sum_abs1 += AND_VEC_D(src1);              \
            sum_abs2 += AND_VEC_D(src2);              \
            sum_abs3 += AND_VEC_D(src3);              \
            sum_abs0 += AND_VEC_D(src4);              \
            sum_abs1 += AND_VEC_D(src5);              \
        }                                             \
        else if ((n & 4) && (n & 1))                  \
        {                                             \
            LD_DP5_INC(x, inc_val, src0, src1, src2,  \
                       src3, src4);                   \
                                                      \
            sum_abs0 += AND_VEC_D(src0);              \
            sum_abs1 += AND_VEC_D(src1);              \
            sum_abs2 += AND_VEC_D(src2);              \
            sum_abs3 += AND_VEC_D(src3);              \
            sum_abs0 += AND_VEC_D(src4);              \
        }                                             \
        else if ((n & 2) && (n & 1))                  \
        {                                             \
            LD_DP3_INC(x, inc_val, src0, src1, src2); \
                                                      \
            sum_abs0 += AND_VEC_D(src0);              \
            sum_abs1 += AND_VEC_D(src1);              \
            sum_abs2 += AND_VEC_D(src2);              \
        }                                             \
        else if (n & 4)                               \
        {                                             \
            LD_DP4_INC(x, inc_val, src0, src1, src2,  \
                       src3);                         \
                                                      \
            sum_abs0 += AND_VEC_D(src0);              \
            sum_abs1 += AND_VEC_D(src1);              \
            sum_abs2 += AND_VEC_D(src2);              \
            sum_abs3 += AND_VEC_D(src3);              \
        }                                             \
        else if (n & 2)                               \
        {                                             \
            LD_DP2_INC(x, inc_val, src0, src1);       \
                                                      \
            sum_abs0 += AND_VEC_D(src0);              \
            sum_abs1 += AND_VEC_D(src1);              \
        }                                             \
        else if (n & 1)                               \
        {                                             \
            src0 = LD_DP(x);                          \
                                                      \
            sum_abs0 += AND_VEC_D(src0);              \
        }                                             \
    }                                                 \
                                                      \
    sum_abs0 += sum_abs1 + sum_abs2 + sum_abs3;       \
    sumf = sum_abs0[0] + sum_abs0[1];

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
    BLASLONG i;
    FLOAT sumf = 0.0;
    v2f64 src0, src1, src2, src3, src4, src5, src6, src7;
    v2f64 sum_abs0, sum_abs1, sum_abs2, sum_abs3;
    v2f64 zero_v = {0};
    v2i64 and_vec = {0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF};

    if (n <= 0 || inc_x <= 0) return (sumf);

    if (1 == inc_x)
    {
        PROCESS_ZD(2);
    }
    else
    {
        inc_x *= 2;
        PROCESS_ZD(inc_x);
    }

    return (sumf);
}
