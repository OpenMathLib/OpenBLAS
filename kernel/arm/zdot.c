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
*	 BLASTEST float		: FAIL
* 	 BLASTEST double	: FAIL
* 	 CTEST			: OK
* 	 TEST			: OK
*
**************************************************************************************/

#include "common.h"
#include <complex.h>

FLOAT _Complex CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y)
{
	BLASLONG i=0;
	BLASLONG ix=0,iy=0;
	FLOAT dot[2];
	FLOAT _Complex result;
	
	dot[0]=0.0;
	dot[1]=0.0;

	__real__ result = 0.0 ;
	__imag__ result = 0.0 ;

	if ( n < 0 )  return(result);

	BLASLONG inc_x2 = 2 * inc_x ;
	BLASLONG inc_y2 = 2 * inc_y ;

	while(i < n)
	{
#if !defined(CONJ)
		dot[0] += ( x[ix]   * y[iy] - x[ix+1] * y[iy+1] ) ;
		dot[1] += ( x[ix+1] * y[iy] + x[ix]   * y[iy+1] ) ;
#else
		dot[0] += ( x[ix]   * y[iy] + x[ix+1] * y[iy+1] ) ;
		dot[1] -= ( x[ix+1] * y[iy] - x[ix]   * y[iy+1] ) ;
#endif
		ix  += inc_x2 ;
		iy  += inc_y2 ;
		i++ ;

	}
	__real__ result = dot[0];
	__imag__ result = dot[1];
	return(result);

}
	

