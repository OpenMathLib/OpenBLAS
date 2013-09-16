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
*	 BLASTEST float		: NoTest
* 	 BLASTEST double	: NoTest
* 	 CTEST			: NoTest
* 	 TEST			: NoTest
*
**************************************************************************************/

#include "common.h"
#include <math.h>


FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
	BLASLONG i=0;
	BLASLONG ix=0;
	FLOAT maxf=0.0;

	if (n < 0 || inc_x < 1 ) return(maxf);

	maxf=x[0];

	while(i < n)
	{
		if( x[ix] > maxf ) 
		{
			maxf = x[ix];
		}
		ix += inc_x;
		i++;
	}
	return(maxf);
}
	

