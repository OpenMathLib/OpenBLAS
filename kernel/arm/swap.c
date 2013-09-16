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
* 2013/08/20 Saar
*	 BLASTEST float		OK
* 	 BLASTEST double	OK
*
**************************************************************************************/

#include "common.h"
#include <stdio.h>

int CNAME(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT dummy3, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *dummy, BLASLONG dummy2)
{
	BLASLONG i=0;
	BLASLONG ix=0,iy=0;
	FLOAT temp;

	if ( n < 0     )  return(0);

	while(i < n)
	{

		temp  = x[ix] ;
		x[ix] = y[iy] ;
		y[iy] = temp ;

		ix += inc_x ;
		iy += inc_y ;
		i++ ;

	}
	return(0);

}
	

