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
 * * 2013/09/14 Saar
 * *	 BLASTEST float		: OK
 * * 	 BLASTEST double	: OK
 * 	 CTEST			: OK
 * 	 TEST			: OK
 * *
 * **************************************************************************************/


#include "common.h"

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *buffer)
{
	BLASLONG i;
	BLASLONG ix,iy;
	BLASLONG j;
	FLOAT *a_ptr;
	FLOAT temp;

	iy = 0;
	a_ptr = a;

	for (j=0; j<n; j++)
	{
		temp = 0.0;
		ix = 0;
		for (i=0; i<m; i++)
		{
			temp += a_ptr[i] * x[ix];
			ix    += inc_x;
		}
		y[iy] += alpha * temp;
		iy += inc_y;
		a_ptr += lda;
	}

}
	

