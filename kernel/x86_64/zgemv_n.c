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

#if defined(HASWELL)
#include "zgemv_n_microk_haswell-2.c"
#endif


#define NBMAX 1024

#ifndef HAVE_KERNEL_16x4

static void zgemv_kernel_16x4(BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y)
{
	BLASLONG i;
	FLOAT *a0,*a1,*a2,*a3;
	a0 = ap[0];
	a1 = ap[1];
	a2 = ap[2];
	a3 = ap[3];

	for ( i=0; i< 2*n; i+=2 )
	{
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
		y[i]   += a0[i]*x[0] - a0[i+1] * x[1];
		y[i+1] += a0[i]*x[1] + a0[i+1] * x[0];
		y[i]   += a1[i]*x[2] - a1[i+1] * x[3];
		y[i+1] += a1[i]*x[3] + a1[i+1] * x[2];
		y[i]   += a2[i]*x[4] - a2[i+1] * x[5];
		y[i+1] += a2[i]*x[5] + a2[i+1] * x[4];
		y[i]   += a3[i]*x[6] - a3[i+1] * x[7];
		y[i+1] += a3[i]*x[7] + a3[i+1] * x[6];
#else 
		y[i]   += a0[i]*x[0] + a0[i+1] * x[1];
		y[i+1] += a0[i]*x[1] - a0[i+1] * x[0];
		y[i]   += a1[i]*x[2] + a1[i+1] * x[3];
		y[i+1] += a1[i]*x[3] - a1[i+1] * x[2];
		y[i]   += a2[i]*x[4] + a2[i+1] * x[5];
		y[i+1] += a2[i]*x[5] - a2[i+1] * x[4];
		y[i]   += a3[i]*x[6] + a3[i+1] * x[7];
		y[i+1] += a3[i]*x[7] - a3[i+1] * x[6];
#endif
	}
}
	
#endif

static void zgemv_kernel_16x1(BLASLONG n, FLOAT *ap, FLOAT *x, FLOAT *y)
{
	BLASLONG i;
	FLOAT *a0;
	a0 = ap;

	for ( i=0; i< 2*n; i+=2 )
	{
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
		y[i]   += a0[i]*x[0] - a0[i+1] * x[1];
		y[i+1] += a0[i]*x[1] + a0[i+1] * x[0];
#else 
		y[i]   += a0[i]*x[0] + a0[i+1] * x[1];
		y[i+1] += a0[i]*x[1] - a0[i+1] * x[0];
#endif

	}
}
	

static void zero_y(BLASLONG n, FLOAT *dest)
{
	BLASLONG i;
	for ( i=0; i<2*n; i++ )
	{
		*dest = 0.0;
		dest++;
	}
}



static void add_y(BLASLONG n, FLOAT *src, FLOAT *dest, BLASLONG inc_dest,FLOAT alpha_r, FLOAT alpha_i)
{
	BLASLONG i;
	FLOAT temp_r;
	FLOAT temp_i;
	for ( i=0; i<n; i++ )
	{
#if !defined(XCONJ) 
		temp_r = alpha_r * src[0] - alpha_i * src[1];
		temp_i = alpha_r * src[1] + alpha_i * src[0];
#else
		temp_r =  alpha_r * src[0] + alpha_i * src[1];
		temp_i = -alpha_r * src[1] + alpha_i * src[0];
#endif

		*dest += temp_r;
		*(dest+1) += temp_i;

		src+=2;
		dest += inc_dest;
	}
}

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha_r,FLOAT alpha_i, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *buffer)
{
	BLASLONG i;
	BLASLONG j;
	FLOAT *a_ptr;
	FLOAT *x_ptr;
	FLOAT *y_ptr;
	FLOAT *ap[4];
	BLASLONG n1;
	BLASLONG m1;
	BLASLONG m2;
	BLASLONG n2;
	FLOAT xbuffer[8],*ybuffer;

	ybuffer = buffer;
	
	inc_x *= 2;
	inc_y *= 2;
	lda   *= 2;

	n1 = n / 4 ;
	n2 = n % 4 ;
	
	m1 = m - ( m % 16 );
	m2 = (m % NBMAX) - (m % 16) ;
	
	y_ptr = y;

	BLASLONG NB = NBMAX;

	while ( NB == NBMAX )
	{
		
		m1 -= NB;
		if ( m1 < 0)
		{
			if ( m2 == 0 ) break;	
			NB = m2;
		}
		
		a_ptr = a;
		x_ptr = x;
		zero_y(NB,ybuffer);
		for( i = 0; i < n1 ; i++)
		{

			xbuffer[0] = x_ptr[0];
			xbuffer[1] = x_ptr[1];
			x_ptr += inc_x;	
			xbuffer[2] = x_ptr[0];
			xbuffer[3] = x_ptr[1];
			x_ptr += inc_x;	
			xbuffer[4] = x_ptr[0];
			xbuffer[5] = x_ptr[1];
			x_ptr += inc_x;	
			xbuffer[6] = x_ptr[0];
			xbuffer[7] = x_ptr[1];
			x_ptr += inc_x;	

			ap[0] = a_ptr;
			ap[1] = a_ptr + lda;
			ap[2] = ap[1] + lda;
			ap[3] = ap[2] + lda;
			zgemv_kernel_16x4(NB,ap,xbuffer,ybuffer);
			a_ptr += 4 * lda;
		}

		for( i = 0; i < n2 ; i++)
		{
			xbuffer[0] = x_ptr[0];
			xbuffer[1] = x_ptr[1];
			x_ptr += inc_x;	
			zgemv_kernel_16x1(NB,a_ptr,xbuffer,ybuffer);
			a_ptr += 1 * lda;

		}
		add_y(NB,ybuffer,y_ptr,inc_y,alpha_r,alpha_i);
		a     += 2 * NB;
		y_ptr += NB * inc_y;
	}

	j=0;
	while ( j < (m % 16))
	{
		a_ptr = a;
		x_ptr = x;
		FLOAT temp_r = 0.0;
		FLOAT temp_i = 0.0;
		for( i = 0; i < n; i++ )
		{
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
			temp_r += a_ptr[0] * x_ptr[0] - a_ptr[1] * x_ptr[1];
			temp_i += a_ptr[0] * x_ptr[1] + a_ptr[1] * x_ptr[0];
#else
			temp_r += a_ptr[0] * x_ptr[0] + a_ptr[1] * x_ptr[1];
			temp_i += a_ptr[0] * x_ptr[1] - a_ptr[1] * x_ptr[0];
#endif

			a_ptr += lda;
			x_ptr += inc_x;
		}

#if !defined(XCONJ) 
		y_ptr[0] += alpha_r * temp_r - alpha_i * temp_i;
		y_ptr[1] += alpha_r * temp_i + alpha_i * temp_r;
#else
		y_ptr[0] += alpha_r * temp_r + alpha_i * temp_i;
		y_ptr[1] -= alpha_r * temp_i - alpha_i * temp_r;
#endif
		y_ptr += inc_y;
		a+=2;
		j++;
	}
	return(0);
}


