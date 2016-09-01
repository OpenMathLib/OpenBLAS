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
	#define OP2		+=
	#define OP3		-
	#define OP4		+
#else
	#define OP2		-=
	#define OP3		+
	#define OP4		-
#endif

#define DOT16_KERNEL(OPR0, OPR1)  \
	dot0 += (vx0r * vy0r);		  \
	dot0 OPR0## = (vx0i * vy0i);  \
	dot1 OPR1## = (vx0i * vy0r);  \
	dot1 += (vx0r * vy0i); 	      \
								  \
	dot0 += (vx1r * vy1r);        \
	dot0 OPR0## = (vx1i * vy1i);  \
	dot1 OPR1## = (vx1i * vy1r);  \
	dot1 += (vx1r * vy1i); 	      \
								  \
	dot0 += (vx2r * vy2r);        \
	dot0 OPR0## = (vx2i * vy2i);  \
	dot1 OPR1## = (vx2i * vy2r);  \
	dot1 += (vx2r * vy2i); 	      \
								  \
	dot0 += (vx3r * vy3r);        \
	dot0 OPR0## = (vx3i * vy3i);  \
	dot1 OPR1## = (vx3i * vy3r);  \
	dot1 += (vx3r * vy3i);

#define DOT12_KERNEL(OPR0, OPR1)  \
	dot0 += (vx0r * vy0r);		  \
	dot0 OPR0## = (vx0i * vy0i);  \
	dot1 OPR1## = (vx0i * vy0r);  \
	dot1 += (vx0r * vy0i); 	      \
								  \
	dot0 += (vx1r * vy1r);        \
	dot0 OPR0## = (vx1i * vy1i);  \
	dot1 OPR1## = (vx1i * vy1r);  \
	dot1 += (vx1r * vy1i);		  \
								  \
	dot0 += (vx2r * vy2r);        \
	dot0 OPR0## = (vx2i * vy2i);  \
	dot1 OPR1## = (vx2i * vy2r);  \
	dot1 += (vx2r * vy2i);

#define DOT8_KERNEL(OPR0, OPR1)   \
	dot0 += (vx0r * vy0r);		  \
	dot0 OPR0## = (vx0i * vy0i);  \
	dot1 OPR1## = (vx0i * vy0r);  \
	dot1 += (vx0r * vy0i); 	      \
								  \
	dot0 += (vx1r * vy1r);        \
	dot0 OPR0## = (vx1i * vy1i);  \
	dot1 OPR1## = (vx1i * vy1r);  \
	dot1 += (vx1r * vy1i);

#define DOT4_KERNEL(OPR0, OPR1)   \
	dot0 += (vx0r * vy0r);		  \
	dot0 OPR0## = (vx0i * vy0i);  \
	dot1 OPR1## = (vx0i * vy0r);  \
	dot1 += (vx0r * vy0i);

/* return float, x,y float */
/* cdotc -  CONJ */
/* cdotu - !CONJ */
#ifndef _MSC_VER
#include <complex.h>
FLOAT _Complex CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y)
#else
OPENBLAS_COMPLEX_FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y)
#endif
{
    BLASLONG i = 0;
    FLOAT dot[2];
    BLASLONG inc_x2;
    BLASLONG inc_y2;
    FLOAT x0, x1, x2, x3, x4, x5, x6, x7;
    FLOAT y0, y1, y2, y3, y4, y5, y6, y7;
    v4f32 vx0, vx1, vx2, vx3, vx4, vx5, vx6, vx7;
    v4f32 vy0, vy1, vy2, vy3, vy4, vy5, vy6, vy7;
	v4f32 vx0r, vx0i, vx1r, vx1i, vx2r, vx2i, vx3r, vx3i;
	v4f32 vy0r, vy0i, vy1r, vy1i, vy2r, vy2i, vy3r, vy3i;
    v4f32 dot0 = {0, 0, 0, 0};
    v4f32 dot1 = {0, 0, 0, 0};
    openblas_complex_float result;

    dot[0] = 0.0;
    dot[1] = 0.0;

    __real__(result) = 0.0;
    __imag__(result) = 0.0;

    if ( n < 1 ) return(result);

    if ((1 == inc_x) && (1 == inc_y))
    {
        for (i = (n >> 4); i--;)
        {
			LD_SP8_INC(x, 4, vx0, vx1, vx2, vx3, vx4, vx5, vx6, vx7);
			LD_SP8_INC(y, 4, vy0, vy1, vy2, vy3, vy4, vy5, vy6, vy7);

			PCKEVOD_W2_SP(vx1, vx0, vx0r, vx0i);
			PCKEVOD_W2_SP(vx3, vx2, vx1r, vx1i);
			PCKEVOD_W2_SP(vx5, vx4, vx2r, vx2i);
			PCKEVOD_W2_SP(vx7, vx6, vx3r, vx3i);

			PCKEVOD_W2_SP(vy1, vy0, vy0r, vy0i);
			PCKEVOD_W2_SP(vy3, vy2, vy1r, vy1i);
			PCKEVOD_W2_SP(vy5, vy4, vy2r, vy2i);
			PCKEVOD_W2_SP(vy7, vy6, vy3r, vy3i);

		#if !defined(CONJ)
			DOT16_KERNEL(-, +);
		#else
			DOT16_KERNEL(+, -);
		#endif
        }

        if (n & 15)
        {
            if ((n & 8) && (n & 4))
            {
				LD_SP4_INC(x, 4, vx0, vx1, vx2, vx3);
				LD_SP4_INC(y, 4, vy0, vy1, vy2, vy3);
				LD_SP2_INC(x, 4, vx4, vx5);
				LD_SP2_INC(y, 4, vy4, vy5);

				PCKEVOD_W2_SP(vx1, vx0, vx0r, vx0i);
				PCKEVOD_W2_SP(vx3, vx2, vx1r, vx1i);
				PCKEVOD_W2_SP(vx5, vx4, vx2r, vx2i);

				PCKEVOD_W2_SP(vy1, vy0, vy0r, vy0i);
				PCKEVOD_W2_SP(vy3, vy2, vy1r, vy1i);
				PCKEVOD_W2_SP(vy5, vy4, vy2r, vy2i);

			#if !defined(CONJ)
				DOT12_KERNEL(-, +);
			#else
				DOT12_KERNEL(+, -);
			#endif
            }
            else if (n & 8)
            {
				LD_SP4_INC(x, 4, vx0, vx1, vx2, vx3);
				LD_SP4_INC(y, 4, vy0, vy1, vy2, vy3);

				PCKEVOD_W2_SP(vx1, vx0, vx0r, vx0i);
				PCKEVOD_W2_SP(vx3, vx2, vx1r, vx1i);

				PCKEVOD_W2_SP(vy1, vy0, vy0r, vy0i);
				PCKEVOD_W2_SP(vy3, vy2, vy1r, vy1i);

			#if !defined(CONJ)
				DOT8_KERNEL(-, +);
			#else
				DOT8_KERNEL(+, -);
			#endif
            }
			else if (n & 4)
            {
				LD_SP2_INC(x, 4, vx0, vx1);
				LD_SP2_INC(y, 4, vy0, vy1);
				PCKEVOD_W2_SP(vx1, vx0, vx0r, vx0i);
				PCKEVOD_W2_SP(vy1, vy0, vy0r, vy0i);

			#if !defined(CONJ)
				DOT4_KERNEL(-, +);
			#else
				DOT4_KERNEL(+, -);
			#endif
            }

			if ((n & 2) && (n & 1))
			{
                LD_GP6_INC(x, 1, x0, x1, x2, x3, x4, x5);
                LD_GP6_INC(y, 1, y0, y1, y2, y3, y4, y5);

				dot[0] += ( x0 * y0 OP3 x1 * y1 );
				dot[1] OP2 ( x1 * y0 OP4 x0 * y1 );

				dot[0] += ( x2 * y2 OP3 x3 * y3 );
				dot[1] OP2 ( x3 * y2 OP4 x2 * y3 );

				dot[0] += ( x4 * y4 OP3 x5 * y5 );
				dot[1] OP2 ( x5 * y4 OP4 x4 * y5 );
			}
			else if (n & 2)
			{
                LD_GP4_INC(x, 1, x0, x1, x2, x3);
                LD_GP4_INC(y, 1, y0, y1, y2, y3);

				dot[0] += ( x0 * y0 OP3 x1 * y1 );
				dot[1] OP2 ( x1 * y0 OP4 x0 * y1 );

				dot[0] += ( x2 * y2 OP3 x3 * y3 );
				dot[1] OP2 ( x3 * y2 OP4 x2 * y3 );
			}
			else if (n & 1)
			{
                LD_GP2_INC(x, 1, x0, x1);
                LD_GP2_INC(y, 1, y0, y1);

				dot[0] += ( x0 * y0 OP3 x1 * y1 );
				dot[1] OP2 ( x1 * y0 OP4 x0 * y1 );
			}
        }

		dot[0] += (dot0[0] + dot0[1] + dot0[2] + dot0[3]);
		dot[1] += (dot1[0] + dot1[1] + dot1[2] + dot1[3]);
	}
	else
	{
		inc_x2 = 2 * inc_x;
		inc_y2 = 2 * inc_y;

		for (i = (n >> 2); i--;)
		{
			x0 = *x;
			x1 = *(x + 1);
			x += inc_x2;
			x2 = *x;
			x3 = *(x + 1);
			x += inc_x2;
			x4 = *x;
			x5 = *(x + 1);
			x += inc_x2;
			x6 = *x;
			x7 = *(x + 1);
			x += inc_x2;

			y0 = *y;
			y1 = *(y + 1);
			y += inc_y2;
			y2 = *y;
			y3 = *(y + 1);
			y += inc_y2;
			y4 = *y;
			y5 = *(y + 1);
			y += inc_y2;
			y6 = *y;
			y7 = *(y + 1);
			y += inc_y2;

			dot[0] += ( x0 * y0 OP3 x1 * y1 );
			dot[1] OP2 ( x1 * y0 OP4 x0 * y1 );

			dot[0] += ( x2 * y2 OP3 x3 * y3 );
			dot[1] OP2 ( x3 * y2 OP4 x2 * y3 );

			dot[0] += ( x4 * y4 OP3 x5 * y5 );
			dot[1] OP2 ( x5 * y4 OP4 x4 * y5 );

			dot[0] += ( x6 * y6 OP3 x7 * y7 );
			dot[1] OP2 ( x7 * y6 OP4 x6 * y7 );
		}

		if ((n & 2) && (n & 1))
		{
			x0 = *x;
			x1 = *(x + 1);
			x += inc_x2;
			x2 = *x;
			x3 = *(x + 1);
			x += inc_x2;
			x4 = *x;
			x5 = *(x + 1);
			x += inc_x2;

			y0 = *y;
			y1 = *(y + 1);
			y += inc_y2;
			y2 = *y;
			y3 = *(y + 1);
			y += inc_y2;
			y4 = *y;
			y5 = *(y + 1);
			y += inc_y2;

			dot[0] += ( x0 * y0 OP3 x1 * y1 );
			dot[1] OP2 ( x1 * y0 OP4 x0 * y1 );

			dot[0] += ( x2 * y2 OP3 x3 * y3 );
			dot[1] OP2 ( x3 * y2 OP4 x2 * y3 );

			dot[0] += ( x4 * y4 OP3 x5 * y5 );
			dot[1] OP2 ( x5 * y4 OP4 x4 * y5 );
		}
		else if (n & 2)
		{
			x0 = *x;
			x1 = *(x + 1);
			x += inc_x2;
			x2 = *x;
			x3 = *(x + 1);
			x += inc_x2;

			y0 = *y;
			y1 = *(y + 1);
			y += inc_y2;
			y2 = *y;
			y3 = *(y + 1);
			y += inc_y2;

			dot[0] += ( x0 * y0 OP3 x1 * y1 );
			dot[1] OP2 ( x1 * y0 OP4 x0 * y1 );

			dot[0] += ( x2 * y2 OP3 x3 * y3 );
			dot[1] OP2 ( x3 * y2 OP4 x2 * y3 );
		}
		else if (n & 1)
		{
			x0 = *x;
			x1 = *(x + 1);
			x += inc_x2;

			y0 = *y;
			y1 = *(y + 1);
			y += inc_y2;

			dot[0] += ( x0 * y0 OP3 x1 * y1 );
			dot[1] OP2 ( x1 * y0 OP4 x0 * y1 );
		}
	}

    __real__(result) = dot[0];
    __imag__(result) = dot[1];

    return(result);
}
