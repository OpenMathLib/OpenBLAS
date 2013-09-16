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
* 2013/09/14 Saar
*	 BLASTEST float		: OK
* 	 BLASTEST double	: OK
* 	 CTEST			: OK
* 	 TEST			: OK
*
**************************************************************************************/

#include "common.h"

int CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT c, FLOAT s)
{
	BLASLONG i=0;
	BLASLONG ix=0,iy=0;
	FLOAT temp[2];

	if ( n <= 0     )  return(0);

	BLASLONG inc_x2 = 2 * inc_x ;
	BLASLONG inc_y2 = 2 * inc_y ;

	while(i < n)
	{
		temp[0]   = c*x[ix]   + s*y[iy] ;
		temp[1]   = c*x[ix+1] + s*y[iy+1] ;
		y[iy]     = c*y[iy]   - s*x[ix] ;
		y[iy+1]   = c*y[iy+1] - s*x[ix+1] ;
		x[ix]     = temp[0] ;
		x[ix+1]   = temp[1] ;

		ix += inc_x2 ;
		iy += inc_y2 ;
		i++ ;

	}
	return(0);

}
	

