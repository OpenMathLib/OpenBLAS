/*************************************************************************************
 * *
 * *    Copyright (C) 2012-2013  Werner Saar (wernsaar@googlemail.com)
 * *
 * *    This program is free software: you can redistribute it and/or modify
 * *    it under the terms of the GNU General Public License as published by
 * *    the Free Software Foundation, either version 3 of the License, or
 * *    (at your option) any later version.
 * *
 * *    This program is distributed in the hope that it will be useful,
 * *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 * *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * *    GNU General Public License for more details.
 * *
 * *    You should have received a copy of the GNU General Public License
 * *    along with this program.  If not, see http://www.gnu.org/licenses/.
 * *
 * **************************************************************************************/

/**************************************************************************************
 * * 2013/09/15 Saar
 * *	 BLASTEST float		: OK
 * * 	 BLASTEST double	: OK
 * 	 CTEST			: OK
 * 	 TEST			: OK
 * *
 * **************************************************************************************/


#include "common.h"

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha_r, FLOAT alpha_i, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *buffer)
{
	BLASLONG i;
	BLASLONG ix,iy;
	BLASLONG j;
	FLOAT *a_ptr;
	FLOAT temp_r,temp_i;
	BLASLONG inc_x2,inc_y2;
	BLASLONG lda2;
	BLASLONG i2;

	if( alpha_r == 0.0 && alpha_i == 0.0 ) return(0);

	lda2 = 2*lda;

	inc_x2 = 2 * inc_x;
	inc_y2 = 2 * inc_y;

	ix = 0;
	a_ptr = a;

#if !defined(CONJ)
	for (j=0; j<n; j++)
	{

#if !defined(XCONJ)
		temp_r = alpha_r * x[ix]   - alpha_i * x[ix+1];
		temp_i = alpha_r * x[ix+1] + alpha_i * x[ix];
#else
		temp_r = alpha_r * x[ix]   + alpha_i * x[ix+1];
		temp_i = alpha_r * x[ix+1] - alpha_i * x[ix];
#endif
		iy = 0;
		for (i=0; i<m; i++)
		{
			i2 = 2*i;
#if !defined(XCONJ)
			y[iy]   += temp_r * a_ptr[i2]   - temp_i * a_ptr[i2+1];
			y[iy+1] += temp_r * a_ptr[i2+1] + temp_i * a_ptr[i2];
#else
			y[iy]   += temp_r * a_ptr[i2]   + temp_i * a_ptr[i2+1];
			y[iy+1] += temp_r * a_ptr[i2+1] - temp_i * a_ptr[i2];
#endif

			iy += inc_y2;
		}
		a_ptr += lda2;
		ix    += inc_x2;
	}

#else
	for (j=0; j<n; j++)
	{

#if !defined(XCONJ)
		temp_r = alpha_r * x[ix]   - alpha_i * x[ix+1];
		temp_i = alpha_r * x[ix+1] + alpha_i * x[ix];
#else
		temp_r = alpha_r * x[ix]   + alpha_i * x[ix+1];
		temp_i = alpha_r * x[ix+1] - alpha_i * x[ix];
#endif
		iy = 0;
		for (i=0; i<m; i++)
		{
			i2 = 2*i;
#if !defined(XCONJ)
			y[iy]   += temp_r * a_ptr[i2]   + temp_i * a_ptr[i2+1];
			y[iy+1] -= temp_r * a_ptr[i2+1] - temp_i * a_ptr[i2];
#else
			y[iy]   += temp_r * a_ptr[i2]   - temp_i * a_ptr[i2+1];
			y[iy+1] -= temp_r * a_ptr[i2+1] + temp_i * a_ptr[i2];
#endif

			iy += inc_y2;
		}
		a_ptr += lda2;
		ix    += inc_x2;
	}


#endif

	return(0);
}
	

