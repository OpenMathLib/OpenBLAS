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
* 	 CTEST			: NoTest
* 	 TEST			: NoTest
*
**************************************************************************************/

#include "common.h"
#include <math.h>

#if defined(DOUBLE)

#define ABS fabs

#else

#define ABS fabsf

#endif

#define CABS1(x,i)	ABS(x[i])+ABS(x[i+1])

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
	BLASLONG i=0;
	BLASLONG ix=0;
	FLOAT minf[2];
	BLASLONG min=0;
	BLASLONG inc_x2;

	if (n < 0 || inc_x < 1 ) return(0.0);

	inc_x2 = 2 * inc_x;

	minf[0] = ABS(x[ix]);
	minf[1] = ABS(x[ix+1]);

	while(i < n)
	{
		if( CABS1(x,ix) < CABS1(minf,0) ) 
		{
			min = i;
			minf[0] = ABS(x[ix]);
			minf[1] = ABS(x[ix+1]);
		}
		ix += inc_x2;
		i++;
	}
	return(CABS1(minf,0));
}
	

