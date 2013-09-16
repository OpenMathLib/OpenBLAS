/*************************************************************************************
*
*    Copyright (C) 2012-2013  Werner Saar (wernsaar@googlemail.com)
*
*    This program is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program.  If not, see http://www.gnu.org/licenses/.
*
**************************************************************************************/

/**************************************************************************************
* 2013/09/15 Saar
*	 BLASTEST float		: OK
* 	 BLASTEST double	: OK
* 	 CTEST			: OK
* 	 TEST			: OK
*
**************************************************************************************/


#include "common.h"

int CNAME(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT da_r, FLOAT da_i, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *dummy, BLASLONG dummy2)
{
	BLASLONG i=0;
	BLASLONG ix,iy;

	if ( n < 0     )  return(0);
	if ( da_r == 0.0 && da_i == 0.0 ) return(0);

	ix = 0;
	iy = 0;

	BLASLONG inc_x2 = 2 * inc_x;
	BLASLONG inc_y2 = 2 * inc_y;

	while(i < n)
	{
#if !defined(CONJ)
		y[iy]   += ( da_r * x[ix]   - da_i * x[ix+1] ) ;
		y[iy+1] += ( da_r * x[ix+1] + da_i * x[ix]   ) ;
#else
		y[iy]   += ( da_r * x[ix]   + da_i * x[ix+1] ) ;
		y[iy+1] -= ( da_r * x[ix+1] - da_i * x[ix]   ) ;
#endif
		ix += inc_x2 ;
		iy += inc_y2 ;
		i++ ;

	}
	return(0);

}
	

