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

#if !defined(CONJ)
    #define OP1     -=
    #define OP2     +=
    #define OP3     -
    #define OP4     +
#else
    #define OP1     +=
    #define OP2     -=
    #define OP3     +
    #define OP4     -
#endif

OPENBLAS_COMPLEX_FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y)
{
    BLASLONG i = 0;
    FLOAT dot[2];
    BLASLONG inc_x2, inc_y2;
    v2f64 vx0, vx1, vx2, vx3, vx4, vx5, vx6, vx7;
    v2f64 vy0, vy1, vy2, vy3, vy4, vy5, vy6, vy7;
    v2f64 vx0r, vx0i, vx1r, vx1i, vx2r, vx2i, vx3r, vx3i;
    v2f64 vy0r, vy0i, vy1r, vy1i, vy2r, vy2i, vy3r, vy3i;
    v2f64 dot0 = {0, 0};
    v2f64 dot1 = {0, 0};
    v2f64 dot2 = {0, 0};
    v2f64 dot3 = {0, 0};
    v2f64 dot4 = {0, 0};
    v2f64 dot5 = {0, 0};
    v2f64 dot6 = {0, 0};
    v2f64 dot7 = {0, 0};
    v2f64 zero = {0, 0};
    OPENBLAS_COMPLEX_FLOAT result;

    dot[0] = 0.0;
    dot[1] = 0.0;

    CREAL(result) = 0.0;
    CIMAG(result) = 0.0;

    if (n < 1) return (result);

    inc_x2 = 2 * inc_x;
    inc_y2 = 2 * inc_y;


#ifdef ENABLE_PREFETCH
    if ((1 == inc_x) && (1 == inc_y))
    {
        double *x_pref, *y_pref;
        BLASLONG pref_offset;

        pref_offset = (BLASLONG)x & (L1_DATA_LINESIZE - 1);
        if (pref_offset > 0)
        {
            pref_offset = L1_DATA_LINESIZE - pref_offset;
        }
        pref_offset = pref_offset / sizeof(double);
        x_pref = x + pref_offset + 32;

        pref_offset = (BLASLONG)y & (L1_DATA_LINESIZE - 1);
        if (pref_offset > 0)
        {
            pref_offset = L1_DATA_LINESIZE - pref_offset;
        }
        pref_offset = pref_offset / sizeof(double);
        y_pref = y + pref_offset + 32;

        for (i = (n >> 3); i--;)
        {
            __asm__ __volatile__(
                "pref   0,   0(%[x_pref])\n\t"
                "pref   0,  32(%[x_pref])\n\t"
                "pref   0,  64(%[x_pref])\n\t"
                "pref   0,  96(%[x_pref])\n\t"
                "pref   0,   0(%[y_pref])\n\t"
                "pref   0,  32(%[y_pref])\n\t"
                "pref   0,  64(%[y_pref])\n\t"
                "pref   0,  96(%[y_pref])\n\t"

                : : [x_pref] "r" (x_pref), [y_pref] "r" (y_pref)
            );

            x_pref += 16;
            y_pref += 16;

            LD_DP8_INC(x, 2, vx0, vx1, vx2, vx3, vx4, vx5, vx6, vx7);
            LD_DP8_INC(y, 2, vy0, vy1, vy2, vy3, vy4, vy5, vy6, vy7);

            PCKEVOD_D2_DP(vx1, vx0, vx0r, vx0i);
            PCKEVOD_D2_DP(vx3, vx2, vx1r, vx1i);
            PCKEVOD_D2_DP(vx5, vx4, vx2r, vx2i);
            PCKEVOD_D2_DP(vx7, vx6, vx3r, vx3i);

            PCKEVOD_D2_DP(vy1, vy0, vy0r, vy0i);
            PCKEVOD_D2_DP(vy3, vy2, vy1r, vy1i);
            PCKEVOD_D2_DP(vy5, vy4, vy2r, vy2i);
            PCKEVOD_D2_DP(vy7, vy6, vy3r, vy3i);

            dot0 += (vx0r * vy0r);
            dot0 OP1 (vx0i * vy0i);
            dot1 OP2 (vx0i * vy0r);
            dot1 += (vx0r * vy0i);

            dot2 += (vx1r * vy1r);
            dot2 OP1 (vx1i * vy1i);
            dot3 OP2 (vx1i * vy1r);
            dot3 += (vx1r * vy1i);

            dot4 += (vx2r * vy2r);
            dot4 OP1 (vx2i * vy2i);
            dot5 OP2 (vx2i * vy2r);
            dot5 += (vx2r * vy2i);

            dot6 += (vx3r * vy3r);
            dot6 OP1 (vx3i * vy3i);
            dot7 OP2 (vx3i * vy3r);
            dot7 += (vx3r * vy3i);
        }
    }
    else
#endif
    for (i = (n >> 3); i--;)
    {
        LD_DP8_INC(x, inc_x2, vx0, vx1, vx2, vx3, vx4, vx5, vx6, vx7);
        LD_DP8_INC(y, inc_y2, vy0, vy1, vy2, vy3, vy4, vy5, vy6, vy7);

        PCKEVOD_D2_DP(vx1, vx0, vx0r, vx0i);
        PCKEVOD_D2_DP(vx3, vx2, vx1r, vx1i);
        PCKEVOD_D2_DP(vx5, vx4, vx2r, vx2i);
        PCKEVOD_D2_DP(vx7, vx6, vx3r, vx3i);

        PCKEVOD_D2_DP(vy1, vy0, vy0r, vy0i);
        PCKEVOD_D2_DP(vy3, vy2, vy1r, vy1i);
        PCKEVOD_D2_DP(vy5, vy4, vy2r, vy2i);
        PCKEVOD_D2_DP(vy7, vy6, vy3r, vy3i);

        dot0 += (vx0r * vy0r);
        dot0 OP1 (vx0i * vy0i);
        dot1 OP2 (vx0i * vy0r);
        dot1 += (vx0r * vy0i);

        dot2 += (vx1r * vy1r);
        dot2 OP1 (vx1i * vy1i);
        dot3 OP2 (vx1i * vy1r);
        dot3 += (vx1r * vy1i);

        dot4 += (vx2r * vy2r);
        dot4 OP1 (vx2i * vy2i);
        dot5 OP2 (vx2i * vy2r);
        dot5 += (vx2r * vy2i);

        dot6 += (vx3r * vy3r);
        dot6 OP1 (vx3i * vy3i);
        dot7 OP2 (vx3i * vy3r);
        dot7 += (vx3r * vy3i);
    }

    if (n & 7)
    {
        if (n & 4)
        {
            LD_DP4_INC(x, inc_x2, vx0, vx1, vx2, vx3);
            LD_DP4_INC(y, inc_y2, vy0, vy1, vy2, vy3);

            PCKEVOD_D2_DP(vx1, vx0, vx0r, vx0i);
            PCKEVOD_D2_DP(vx3, vx2, vx1r, vx1i);

            PCKEVOD_D2_DP(vy1, vy0, vy0r, vy0i);
            PCKEVOD_D2_DP(vy3, vy2, vy1r, vy1i);

            dot0 += (vx0r * vy0r);
            dot0 OP1 (vx0i * vy0i);
            dot1 OP2 (vx0i * vy0r);
            dot1 += (vx0r * vy0i);

            dot2 += (vx1r * vy1r);
            dot2 OP1 (vx1i * vy1i);
            dot3 OP2 (vx1i * vy1r);
            dot3 += (vx1r * vy1i);
        }

        if (n & 2)
        {
            LD_DP2_INC(x, inc_x2, vx0, vx1);
            LD_DP2_INC(y, inc_y2, vy0, vy1);
            PCKEVOD_D2_DP(vx1, vx0, vx0r, vx0i);
            PCKEVOD_D2_DP(vy1, vy0, vy0r, vy0i);

            dot0 += (vx0r * vy0r);
            dot0 OP1 (vx0i * vy0i);
            dot1 OP2 (vx0i * vy0r);
            dot1 += (vx0r * vy0i);
        }

        if (n & 1)
        {
            vx0 = LD_DP(x);
            vy0 = LD_DP(y);
            PCKEVOD_D2_DP(zero, vx0, vx0r, vx0i);
            PCKEVOD_D2_DP(zero, vy0, vy0r, vy0i);

            dot0 += (vx0r * vy0r);
            dot0 OP1 (vx0i * vy0i);
            dot1 OP2 (vx0i * vy0r);
            dot1 += (vx0r * vy0i);
        }
    }

    dot0 += dot2 + dot4 + dot6;
    dot1 += dot3 + dot5 + dot7;

    dot[0] += (dot0[0] + dot0[1]);
    dot[1] += (dot1[0] + dot1[1]);

    CREAL(result) = dot[0];
    CIMAG(result) = dot[1];

    return (result);
}
