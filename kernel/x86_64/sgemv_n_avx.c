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
#include "sgemv_n_microk_bulldozer.c"
#elif defined(HASWELL)
#include "sgemv_n_microk_haswell.c"
#else
#include "sgemv_n_microk_sandy.c"
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

static void add_y(BLASLONG n, FLOAT *src, FLOAT *dest, BLASLONG inc_dest)
{
	BLASLONG i;
	for ( i=0; i<n; i++ )
	{
		*dest += *src;
		src++;
		dest += inc_dest;
	}
}

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *buffer)
{
	BLASLONG i;
	BLASLONG j;
	FLOAT *a_ptr;
	FLOAT *x_ptr;
	FLOAT *y_ptr;
	BLASLONG n1;
	BLASLONG m1;
	BLASLONG register m2;
	BLASLONG register n2;
	FLOAT *xbuffer,*ybuffer;
	xbuffer = buffer;
	ybuffer = xbuffer + 2048 + 256;
	
	n1 = n / 512 ;
	n2 = n % 512 ;

	m1 = m / 64;
	m2 = m % 64;

	y_ptr = y;
	x_ptr = x;

	for (j=0; j<n1; j++)
	{

		if ( inc_x == 1 )
			xbuffer = x_ptr;
		else
			copy_x(512,x_ptr,xbuffer,inc_x);

		a_ptr = a + j * 512 * lda;
		y_ptr = y;

		for(i = 0; i<m1; i++ )
		{
			sgemv_kernel_64(512,alpha,a_ptr,lda,xbuffer,ybuffer);
			add_y(64,ybuffer,y_ptr,inc_y);
			y_ptr += 64 * inc_y;
			a_ptr += 64;			

		}

		if ( m2 & 32 )
		{
			sgemv_kernel_32(512,alpha,a_ptr,lda,xbuffer,ybuffer);
			add_y(32,ybuffer,y_ptr,inc_y);
			y_ptr += 32 * inc_y;
			a_ptr += 32;			

		}

		if ( m2 & 16 )
		{
			sgemv_kernel_16(512,alpha,a_ptr,lda,xbuffer,ybuffer);
			add_y(16,ybuffer,y_ptr,inc_y);
			y_ptr += 16 * inc_y;
			a_ptr += 16;			
		}
		if ( m2 & 8 )
		{
			sgemv_kernel_8(512,alpha,a_ptr,lda,xbuffer,ybuffer);
			add_y(8,ybuffer,y_ptr,inc_y);
			y_ptr += 8 * inc_y;
			a_ptr += 8;			
		}
		if ( m2 & 4 )
		{
			sgemv_kernel_4(512,alpha,a_ptr,lda,xbuffer,ybuffer);
			add_y(4,ybuffer,y_ptr,inc_y);
			y_ptr += 4 * inc_y;
			a_ptr += 4;			
		}
		if ( m2 & 2 )
		{
			sgemv_kernel_2(512,alpha,a_ptr,lda,xbuffer,ybuffer);
			add_y(2,ybuffer,y_ptr,inc_y);
			y_ptr += 2 * inc_y;
			a_ptr += 2;			
		}
		if ( m2 & 1 )
		{
			sgemv_kernel_1(512,alpha,a_ptr,lda,xbuffer,ybuffer);
			add_y(1,ybuffer,y_ptr,inc_y);
		}
		x_ptr += 512 * inc_x;

	}

	if ( n2 > 0 )
	{

		if ( inc_x == 1 )
			xbuffer = x_ptr;
		else
			copy_x(n2,x_ptr,xbuffer,inc_x);

		a_ptr = a + n1 * 512 * lda;
		y_ptr = y;

		for(i = 0; i<m1; i++ )
		{
			sgemv_kernel_64(n2,alpha,a_ptr,lda,xbuffer,ybuffer);
			add_y(64,ybuffer,y_ptr,inc_y);
			y_ptr += 64 * inc_y;
			a_ptr += 64;			

		}

		if ( m2 & 32 )
		{
			sgemv_kernel_32(n2,alpha,a_ptr,lda,xbuffer,ybuffer);
			add_y(32,ybuffer,y_ptr,inc_y);
			y_ptr += 32 * inc_y;
			a_ptr += 32;			

		}
		if ( m2 & 16 )
		{
			sgemv_kernel_16(n2,alpha,a_ptr,lda,xbuffer,ybuffer);
			add_y(16,ybuffer,y_ptr,inc_y);
			y_ptr += 16 * inc_y;
			a_ptr += 16;			
		}
		if ( m2 & 8 )
		{
			sgemv_kernel_8(n2,alpha,a_ptr,lda,xbuffer,ybuffer);
			add_y(8,ybuffer,y_ptr,inc_y);
			y_ptr += 8 * inc_y;
			a_ptr += 8;			
		}
		if ( m2 & 4 )
		{
			sgemv_kernel_4(n2,alpha,a_ptr,lda,xbuffer,ybuffer);
			add_y(4,ybuffer,y_ptr,inc_y);
			y_ptr += 4 * inc_y;
			a_ptr += 4;			
		}
		if ( m2 & 2 )
		{
			sgemv_kernel_2(n2,alpha,a_ptr,lda,xbuffer,ybuffer);
			add_y(2,ybuffer,y_ptr,inc_y);
			y_ptr += 2 * inc_y;
			a_ptr += 2;			
		}
		if ( m2 & 1 )
		{
			sgemv_kernel_1(n2,alpha,a_ptr,lda,xbuffer,ybuffer);
			add_y(1,ybuffer,y_ptr,inc_y);
		}


	}
	return(0);
}


