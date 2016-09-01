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

#define AND_VEC_W(in)   ((v4f32) ((v4i32) in & and_vec))

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
    BLASLONG i = 0;
    FLOAT data0, data1, data2, sumf = 0.0;
    v4f32 src0, src1, src2, src3, src4, src5, src6, src7;
    v4f32 sum_abs0, sum_abs1, sum_abs2, sum_abs3;
    v4f32 zero_v = {0};
    v4i32 and_vec = {0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF};

    if (n <= 0 || inc_x <= 0) return (sumf);

    if (1 == inc_x)
    {
        if (n > 31)
        {
            n -= 32;

            LD_SP8_INC(x, 4, src0, src1, src2, src3, src4, src5, src6, src7);

            sum_abs0 = AND_VEC_W(src0);
            sum_abs1 = AND_VEC_W(src1);
            sum_abs2 = AND_VEC_W(src2);
            sum_abs3 = AND_VEC_W(src3);
            sum_abs0 += AND_VEC_W(src4);
            sum_abs1 += AND_VEC_W(src5);
            sum_abs2 += AND_VEC_W(src6);
            sum_abs3 += AND_VEC_W(src7);
        }
        else
        {
            sum_abs0 = zero_v;
            sum_abs1 = zero_v;
            sum_abs2 = zero_v;
            sum_abs3 = zero_v;
        }

        for (i = 0; i < (n >> 5); i++)
        {
            LD_SP8_INC(x, 4, src0, src1, src2, src3, src4, src5, src6, src7);

            sum_abs0 += AND_VEC_W(src0);
            sum_abs1 += AND_VEC_W(src1);
            sum_abs2 += AND_VEC_W(src2);
            sum_abs3 += AND_VEC_W(src3);
            sum_abs0 += AND_VEC_W(src4);
            sum_abs1 += AND_VEC_W(src5);
            sum_abs2 += AND_VEC_W(src6);
            sum_abs3 += AND_VEC_W(src7);
        }

        if (n & 31)
        {
            if ((n & 16) && (n & 8) && (n & 4))
            {
                LD_SP7_INC(x, 4, src0, src1, src2, src3, src4, src5, src6);

                sum_abs0 += AND_VEC_W(src0);
                sum_abs1 += AND_VEC_W(src1);
                sum_abs2 += AND_VEC_W(src2);
                sum_abs3 += AND_VEC_W(src3);
                sum_abs0 += AND_VEC_W(src4);
                sum_abs1 += AND_VEC_W(src5);
                sum_abs2 += AND_VEC_W(src6);

                sum_abs0 = sum_abs0 + sum_abs1 + sum_abs2 + sum_abs3;

                sumf += sum_abs0[0];
                sumf += sum_abs0[1];
                sumf += sum_abs0[2];
                sumf += sum_abs0[3];
            }
            else if ((n & 16) && (n & 8))
            {
                LD_SP6_INC(x, 4, src0, src1, src2, src3, src4, src5);

                sum_abs0 += AND_VEC_W(src0);
                sum_abs1 += AND_VEC_W(src1);
                sum_abs2 += AND_VEC_W(src2);
                sum_abs3 += AND_VEC_W(src3);
                sum_abs0 += AND_VEC_W(src4);
                sum_abs1 += AND_VEC_W(src5);

                sum_abs0 = sum_abs0 + sum_abs1 + sum_abs2 + sum_abs3;

                sumf += sum_abs0[0];
                sumf += sum_abs0[1];
                sumf += sum_abs0[2];
                sumf += sum_abs0[3];
            }
            else if ((n & 16) && (n & 4))
            {
                LD_SP5_INC(x, 4, src0, src1, src2, src3, src4);

                sum_abs0 += AND_VEC_W(src0);
                sum_abs1 += AND_VEC_W(src1);
                sum_abs2 += AND_VEC_W(src2);
                sum_abs3 += AND_VEC_W(src3);
                sum_abs0 += AND_VEC_W(src4);

                sum_abs0 = sum_abs0 + sum_abs1 + sum_abs2 + sum_abs3;

                sumf += sum_abs0[0];
                sumf += sum_abs0[1];
                sumf += sum_abs0[2];
                sumf += sum_abs0[3];
            }
            else if ((n & 8) && (n & 4))
            {
                LD_SP3_INC(x, 4, src0, src1, src2);

                sum_abs0 += AND_VEC_W(src0);
                sum_abs1 += AND_VEC_W(src1);
                sum_abs2 += AND_VEC_W(src2);

                sum_abs0 = sum_abs0 + sum_abs1 + sum_abs2 + sum_abs3;

                sumf += sum_abs0[0];
                sumf += sum_abs0[1];
                sumf += sum_abs0[2];
                sumf += sum_abs0[3];
            }
            else if (n & 16)
            {
                LD_SP4_INC(x, 4, src0, src1, src2, src3);

                sum_abs0 += AND_VEC_W(src0);
                sum_abs1 += AND_VEC_W(src1);
                sum_abs2 += AND_VEC_W(src2);
                sum_abs3 += AND_VEC_W(src3);

                sum_abs0 = sum_abs0 + sum_abs1 + sum_abs2 + sum_abs3;

                sumf += sum_abs0[0];
                sumf += sum_abs0[1];
                sumf += sum_abs0[2];
                sumf += sum_abs0[3];
            }
            else if (n & 8)
            {
                LD_SP2_INC(x, 4, src0, src1);

                sum_abs0 += AND_VEC_W(src0);
                sum_abs1 += AND_VEC_W(src1);

                sum_abs0 = sum_abs0 + sum_abs1 + sum_abs2 + sum_abs3;

                sumf += sum_abs0[0];
                sumf += sum_abs0[1];
                sumf += sum_abs0[2];
                sumf += sum_abs0[3];
            }
            else if (n & 4)
            {
                src0 = LD_SP(x); x += 4;

                sum_abs0 += AND_VEC_W(src0);

                sum_abs0 = sum_abs0 + sum_abs1 + sum_abs2 + sum_abs3;

                sumf += sum_abs0[0];
                sumf += sum_abs0[1];
                sumf += sum_abs0[2];
                sumf += sum_abs0[3];
            }
            else
            {
                sum_abs0 = sum_abs0 + sum_abs1 + sum_abs2 + sum_abs3;

                sumf += sum_abs0[0];
                sumf += sum_abs0[1];
                sumf += sum_abs0[2];
                sumf += sum_abs0[3];
            }

            if (n & 2)
            {
                sumf += fabsf(*(x + 0));
                sumf += fabsf(*(x + 1));
                x += 2;
            }

            if (n & 1)
            {
                sumf += fabsf(*(x + 0));
            }
        }
        else
        {
            sum_abs0 = sum_abs0 + sum_abs1 + sum_abs2 + sum_abs3;

            sumf += sum_abs0[0];
            sumf += sum_abs0[1];
            sumf += sum_abs0[2];
            sumf += sum_abs0[3];
        }
    }
    else
    {
        if (n > 8)
        {
            n -= 8;

            src0 = (v4f32) __msa_insert_w((v4i32) zero_v, 0, *((int *) x));
            x += inc_x;
            src0 = (v4f32) __msa_insert_w((v4i32) src0, 1, *((int *) x));
            x += inc_x;
            src0 = (v4f32) __msa_insert_w((v4i32) src0, 2, *((int *) x));
            x += inc_x;
            src0 = (v4f32) __msa_insert_w((v4i32) src0, 3, *((int *) x));
            x += inc_x;
            src4 = (v4f32) __msa_insert_w((v4i32) zero_v, 0, *((int *) x));
            x += inc_x;
            src4 = (v4f32) __msa_insert_w((v4i32) src4, 1, *((int *) x));
            x += inc_x;
            src4 = (v4f32) __msa_insert_w((v4i32) src4, 2, *((int *) x));
            x += inc_x;
            src4 = (v4f32) __msa_insert_w((v4i32) src4, 3, *((int *) x));
            x += inc_x;

            sum_abs0 = AND_VEC_W(src0);
            sum_abs1 = AND_VEC_W(src4);
        }
        else
        {
            sum_abs0 = zero_v;
            sum_abs1 = zero_v;
        }

        for (i = (n >> 3); i--;)
        {
            src0 = (v4f32) __msa_insert_w((v4i32) zero_v, 0, *((int *) x));
            x += inc_x;
            src0 = (v4f32) __msa_insert_w((v4i32) src0, 1, *((int *) x));
            x += inc_x;
            src0 = (v4f32) __msa_insert_w((v4i32) src0, 2, *((int *) x));
            x += inc_x;
            src0 = (v4f32) __msa_insert_w((v4i32) src0, 3, *((int *) x));
            x += inc_x;
            src4 = (v4f32) __msa_insert_w((v4i32) zero_v, 0, *((int *) x));
            x += inc_x;
            src4 = (v4f32) __msa_insert_w((v4i32) src4, 1, *((int *) x));
            x += inc_x;
            src4 = (v4f32) __msa_insert_w((v4i32) src4, 2, *((int *) x));
            x += inc_x;
            src4 = (v4f32) __msa_insert_w((v4i32) src4, 3, *((int *) x));
            x += inc_x;

            sum_abs0 += AND_VEC_W(src0);
            sum_abs1 += AND_VEC_W(src4);
        }

        if (n & 4)
        {
            src0 = (v4f32) __msa_insert_w((v4i32) zero_v, 0, *((int *) x));
            x += inc_x;
            src0 = (v4f32) __msa_insert_w((v4i32) src0, 1, *((int *) x));
            x += inc_x;
            src0 = (v4f32) __msa_insert_w((v4i32) src0, 2, *((int *) x));
            x += inc_x;
            src0 = (v4f32) __msa_insert_w((v4i32) src0, 3, *((int *) x));
            x += inc_x;

            sum_abs0 += AND_VEC_W(src0);
        }

        sum_abs0 += sum_abs1;

        sumf += sum_abs0[0];
        sumf += sum_abs0[1];
        sumf += sum_abs0[2];
        sumf += sum_abs0[3];

        if ((n & 2) && (n & 1))
        {
            data0 = fabsf(*x); x += inc_x;
            data1 = fabsf(*x); x += inc_x;
            data2 = fabsf(*x);

            sumf += data0;
            sumf += data1;
            sumf += data2;
        }
        else if (n & 2)
        {
            data0 = fabsf(*x); x += inc_x;
            data1 = fabsf(*x);

            sumf += data0;
            sumf += data1;
        }
        else if (n & 1)
        {
            data0 = fabsf(*x);

            sumf += data0;
        }
    }

    return (sumf);
}
