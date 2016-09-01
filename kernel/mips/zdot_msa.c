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

/* return double, x,y double */
/* zdotc -  CONJ */
/* zdotu - !CONJ */
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
    v2f64 vx0, vx1, vx2, vx3, vx4, vx5, vx6, vx7;
    v2f64 vy0, vy1, vy2, vy3, vy4, vy5, vy6, vy7;
	v2f64 vx0r, vx0i, vx1r, vx1i, vx2r, vx2i, vx3r, vx3i;
	v2f64 vy0r, vy0i, vy1r, vy1i, vy2r, vy2i, vy3r, vy3i;
    v2f64 dot0 = {0, 0};
    v2f64 dot1 = {0, 0};
    v2f64 zero = {0, 0};
    openblas_complex_double result;

    dot[0] = 0.0;
    dot[1] = 0.0;

    __real__(result) = 0.0;
    __imag__(result) = 0.0;

    if ( n < 1 ) return(result);

    inc_x2 = 2 * inc_x;
    inc_y2 = 2 * inc_y;

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

	#if !defined(CONJ)
		DOT16_KERNEL(-, +);
	#else
		DOT16_KERNEL(+, -);
	#endif
	}

	if (n & 7)
	{
		if ((n & 4) && (n & 2))
		{
			LD_DP4_INC(x, inc_x2, vx0, vx1, vx2, vx3);
			LD_DP4_INC(y, inc_y2, vy0, vy1, vy2, vy3);
			LD_DP2_INC(x, inc_x2, vx4, vx5);
			LD_DP2_INC(y, inc_y2, vy4, vy5);

			PCKEVOD_D2_DP(vx1, vx0, vx0r, vx0i);
			PCKEVOD_D2_DP(vx3, vx2, vx1r, vx1i);
			PCKEVOD_D2_DP(vx5, vx4, vx2r, vx2i);

			PCKEVOD_D2_DP(vy1, vy0, vy0r, vy0i);
			PCKEVOD_D2_DP(vy3, vy2, vy1r, vy1i);
			PCKEVOD_D2_DP(vy5, vy4, vy2r, vy2i);

		#if !defined(CONJ)
			DOT12_KERNEL(-, +);
		#else
			DOT12_KERNEL(+, -);
		#endif
		}
		else if (n & 4)
		{
			LD_DP4_INC(x, inc_x2, vx0, vx1, vx2, vx3);
			LD_DP4_INC(y, inc_y2, vy0, vy1, vy2, vy3);

			PCKEVOD_D2_DP(vx1, vx0, vx0r, vx0i);
			PCKEVOD_D2_DP(vx3, vx2, vx1r, vx1i);

			PCKEVOD_D2_DP(vy1, vy0, vy0r, vy0i);
			PCKEVOD_D2_DP(vy3, vy2, vy1r, vy1i);

		#if !defined(CONJ)
			DOT8_KERNEL(-, +);
		#else
			DOT8_KERNEL(+, -);
		#endif
		}
		else if (n & 2)
		{
			LD_DP2_INC(x, inc_x2, vx0, vx1);
			LD_DP2_INC(y, inc_y2, vy0, vy1);
			PCKEVOD_D2_DP(vx1, vx0, vx0r, vx0i);
			PCKEVOD_D2_DP(vy1, vy0, vy0r, vy0i);

		#if !defined(CONJ)
			DOT4_KERNEL(-, +);
		#else
			DOT4_KERNEL(+, -);
		#endif
		}

		if (n & 1)
		{
			vx0 = LD_DP(x);
			vy0 = LD_DP(y);
			PCKEVOD_D2_DP(zero, vx0, vx0r, vx0i);
			PCKEVOD_D2_DP(zero, vy0, vy0r, vy0i);

		#if !defined(CONJ)
			DOT4_KERNEL(-, +);
		#else
			DOT4_KERNEL(+, -);
		#endif
		}
	}

	dot[0] += (dot0[0] + dot0[1]);
	dot[1] += (dot1[0] + dot1[1]);

    __real__(result) = dot[0];
    __imag__(result) = dot[1];

    return(result);
}
