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
#include "macros_msa.h"

/* return float, x,y float */
#if defined(DSDOT)
double CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y)
#else
FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y)
#endif
{
    BLASLONG i = 0;
    double dot = 0.0;
    float x0, x1, x2, x3, y0, y1, y2, y3;
    v4f32 vx0, vx1, vx2, vx3, vx4, vx5, vx6, vx7;
    v4f32 vy0, vy1, vy2, vy3, vy4, vy5, vy6, vy7;
    v4f32 dot0 = {0, 0, 0, 0};

    if (n < 0) return (dot);

    if ((1 == inc_x) && (1 == inc_y))
    {
        for (i = (n >> 5); i--;)
        {
			LD_SP8_INC(x, 4, vx0, vx1, vx2, vx3, vx4, vx5, vx6, vx7);
			LD_SP8_INC(y, 4, vy0, vy1, vy2, vy3, vy4, vy5, vy6, vy7);

            dot0 += (vy0 * vx0);
            dot0 += (vy1 * vx1);
            dot0 += (vy2 * vx2);
            dot0 += (vy3 * vx3);
            dot0 += (vy4 * vx4);
            dot0 += (vy5 * vx5);
            dot0 += (vy6 * vx6);
            dot0 += (vy7 * vx7);
        }

        if (n & 31)
        {
            if ((n & 16) && (n & 8) && (n & 4))
            {
                LD_SP7_INC(x, 4, vx0, vx1, vx2, vx3, vx4, vx5, vx6);
                LD_SP7_INC(y, 4, vy0, vy1, vy2, vy3, vy4, vy5, vy6);

                dot0 += (vy0 * vx0);
                dot0 += (vy1 * vx1);
                dot0 += (vy2 * vx2);
                dot0 += (vy3 * vx3);
                dot0 += (vy4 * vx4);
                dot0 += (vy5 * vx5);
                dot0 += (vy6 * vx6);
            }
            else if ((n & 16) && (n & 8))
            {
                LD_SP6_INC(x, 4, vx0, vx1, vx2, vx3, vx4, vx5);
                LD_SP6_INC(y, 4, vy0, vy1, vy2, vy3, vy4, vy5);

                dot0 += (vy0 * vx0);
                dot0 += (vy1 * vx1);
                dot0 += (vy2 * vx2);
                dot0 += (vy3 * vx3);
                dot0 += (vy4 * vx4);
                dot0 += (vy5 * vx5);
            }
            else if ((n & 16) && (n & 4))
            {
                LD_SP5_INC(x, 4, vx0, vx1, vx2, vx3, vx4);
                LD_SP5_INC(y, 4, vy0, vy1, vy2, vy3, vy4);

                dot0 += (vy0 * vx0);
                dot0 += (vy1 * vx1);
                dot0 += (vy2 * vx2);
                dot0 += (vy3 * vx3);
                dot0 += (vy4 * vx4);
            }
            else if ((n & 8) && (n & 4))
            {
                LD_SP3_INC(x, 4, vx0, vx1, vx2);
                LD_SP3_INC(y, 4, vy0, vy1, vy2);

                dot0 += (vy0 * vx0);
                dot0 += (vy1 * vx1);
                dot0 += (vy2 * vx2);
            }
            else if (n & 16)
            {
				LD_SP4_INC(x, 4, vx0, vx1, vx2, vx3);
				LD_SP4_INC(y, 4, vy0, vy1, vy2, vy3);

                dot0 += (vy0 * vx0);
                dot0 += (vy1 * vx1);
                dot0 += (vy2 * vx2);
                dot0 += (vy3 * vx3);
            }
            else if (n & 8)
            {
				LD_SP2_INC(x, 4, vx0, vx1);
				LD_SP2_INC(y, 4, vy0, vy1);

                dot0 += (vy0 * vx0);
                dot0 += (vy1 * vx1);
            }
            else if (n & 4)
            {
                vx0 = LD_SP(x); x += 4;
                vy0 = LD_SP(y); y += 4;

                dot0 += (vy0 * vx0);
            }

            if ((n & 2) && (n & 1))
            {
                LD_GP3_INC(x, 1, x0, x1, x2);
                LD_GP3_INC(y, 1, y0, y1, y2);

                dot += (y0 * x0);
                dot += (y1 * x1);
                dot += (y2 * x2);
            }
            else if (n & 2)
            {
                LD_GP2_INC(x, 1, x0, x1);
                LD_GP2_INC(y, 1, y0, y1);

                dot += (y0 * x0);
                dot += (y1 * x1);
            }
            else if (n & 1)
            {
                x0 = *x;
                y0 = *y;

                dot += (y0 * x0);
            }
        }

        dot += dot0[0];
        dot += dot0[1];
        dot += dot0[2];
        dot += dot0[3];
    }
    else
    {
        for (i = (n >> 2); i--;)
        {
            LD_GP4_INC(x, inc_x, x0, x1, x2, x3);
            LD_GP4_INC(y, inc_y, y0, y1, y2, y3);

            dot += (y0 * x0);
            dot += (y1 * x1);
            dot += (y2 * x2);
            dot += (y3 * x3);
        }

        if ((n & 2) && (n & 1))
        {
            LD_GP3_INC(x, inc_x, x0, x1, x2);
            LD_GP3_INC(y, inc_y, y0, y1, y2);

            dot += (y0 * x0);
            dot += (y1 * x1);
            dot += (y2 * x2);
        }
        else if (n & 2)
        {
            LD_GP2_INC(x, inc_x, x0, x1);
            LD_GP2_INC(y, inc_y, y0, y1);

            dot += (y0 * x0);
            dot += (y1 * x1);
        }
        else if (n & 1)
        {
            x0 = *x;
            y0 = *y;

            dot += (y0 * x0);
        }
    }

    return (dot);
}
