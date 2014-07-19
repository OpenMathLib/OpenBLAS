/***************************************************************************
Copyright (c) 2014, The OpenBLAS Project
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
*****************************************************************************/


#include "common.h"

#if defined(BULLDOZER) || defined(PILEDRIVER)
#include "sgemv_t_microk_bulldozer.c"
#endif

static void copy_x(BLASLONG n, FLOAT *src, FLOAT *dest, BLASLONG inc_src)
{
	BLASLONG i;
	for ( i=0; i<n; i++ )
	{
		*dest = *src;
		dest++;
		src += inc_src;
	}
}

static void  sgemv_kernel_1( BLASLONG n, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *x, FLOAT *y)
{

	FLOAT register temp0 = 0.0;
	BLASLONG i;
	for ( i=0; i<n ; i++)
	{
		temp0 += a[i] * x[i];
	}
	temp0 *= alpha ;
	*y += temp0;
}




int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *buffer)
{
	BLASLONG i;
	BLASLONG j;
	FLOAT *a_ptr;
	FLOAT *x_ptr;
	FLOAT *y_ptr;
	FLOAT *a_ptrl;
	BLASLONG m1;
	BLASLONG register m2;
	FLOAT *xbuffer;
	xbuffer = buffer;
	BLASLONG register Mblock;

	m1 = m / 1024 ;
	m2 = m % 1024 ;

	x_ptr = x;
	a_ptr = a;

	for (j=0; j<m1; j++)
	{

		if ( inc_x == 1 )
			xbuffer = x_ptr;
		else
			copy_x(1024,x_ptr,xbuffer,inc_x);

		y_ptr = y;
		a_ptrl = a_ptr;

		for(i = 0; i<n; i++ )
		{
			sgemv_kernel_16(1024,alpha,a_ptrl,lda,xbuffer,y_ptr);
			y_ptr += inc_y;
			a_ptrl += lda;
		}
		a_ptr += 1024;	
		x_ptr += 1024 * inc_x;
	}

	if ( m2 == 0 ) return(0);

	Mblock = 512;
	while ( Mblock >= 16 )
	{
	  if ( m2 & Mblock)
	  {

		if ( inc_x == 1 )
			xbuffer = x_ptr;
		else
			copy_x(Mblock,x_ptr,xbuffer,inc_x);

		y_ptr = y;
		a_ptrl = a_ptr;

		for(i = 0; i<n; i++ )
		{
			sgemv_kernel_16(Mblock,alpha,a_ptrl,lda,xbuffer,y_ptr);
			y_ptr += inc_y;
			a_ptrl += lda;
		}
		a_ptr += Mblock;	
		x_ptr += Mblock * inc_x;


 	  }
	  Mblock /= 2;

	}

        if ( m2 & Mblock)
	{

		if ( inc_x == 1 )
			xbuffer = x_ptr;
		else
			copy_x(Mblock,x_ptr,xbuffer,inc_x);

		y_ptr = y;
		a_ptrl = a_ptr;

		for(i = 0; i<n; i++ )
		{
			sgemv_kernel_1(Mblock,alpha,a_ptrl,lda,xbuffer,y_ptr);
			y_ptr += inc_y;
			a_ptrl += lda;
		}
		a_ptr += Mblock;	
		x_ptr += Mblock * inc_x;


 	}
	Mblock /= 2;


        if ( m2 & Mblock)
	{

		if ( inc_x == 1 )
			xbuffer = x_ptr;
		else
			copy_x(Mblock,x_ptr,xbuffer,inc_x);

		y_ptr = y;
		a_ptrl = a_ptr;

		for(i = 0; i<n; i++ )
		{
			sgemv_kernel_1(Mblock,alpha,a_ptrl,lda,xbuffer,y_ptr);
			y_ptr += inc_y;
			a_ptrl += lda;
		}
		a_ptr += Mblock;	
		x_ptr += Mblock * inc_x;


 	}
	Mblock /= 2;

        if ( m2 & Mblock)
	{

		if ( inc_x == 1 )
			xbuffer = x_ptr;
		else
			copy_x(Mblock,x_ptr,xbuffer,inc_x);

		y_ptr = y;
		a_ptrl = a_ptr;

		for(i = 0; i<n; i++ )
		{
			sgemv_kernel_1(Mblock,alpha,a_ptrl,lda,xbuffer,y_ptr);
			y_ptr += inc_y;
			a_ptrl += lda;
		}
		a_ptr += Mblock;	
		x_ptr += Mblock * inc_x;


 	}
	Mblock /= 2;

        if ( m2 & Mblock)
	{

		xbuffer = x_ptr;

		y_ptr = y;
		a_ptrl = a_ptr;

		for(i = 0; i<n; i++ )
		{
			sgemv_kernel_1(Mblock,alpha,a_ptrl,lda,xbuffer,y_ptr);
			y_ptr += inc_y;
			a_ptrl += lda;
		}


 	}

	return(0);
}


