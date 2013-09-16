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


FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
	BLASLONG i=0;
	BLASLONG ix=0;
	FLOAT minf=0.0;

	if (n < 0 || inc_x < 1 ) return(minf);

	minf=ABS(x[0]);

	while(i < n)
	{
		if( ABS(x[ix]) < ABS(minf) ) 
		{
			minf = ABS(x[ix]);
		}
		ix += inc_x;
		i++;
	}
	return(minf);
}
	

