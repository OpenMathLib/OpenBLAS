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


#if defined(BULLDOZER)
#include "zgemv_t_microk_bulldozer-2.c"
#elif defined(HASWELL)
#include "zgemv_t_microk_haswell-2.c"
#endif


#define NBMAX 1028

#ifndef HAVE_KERNEL_16x4

static void zgemv_kernel_16x4(BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y)
{
	BLASLONG i;
	FLOAT *a0,*a1,*a2,*a3;
	a0 = ap[0];
	a1 = ap[1];
	a2 = ap[2];
	a3 = ap[3];
	FLOAT temp_r0 = 0.0;
	FLOAT temp_r1 = 0.0;
	FLOAT temp_r2 = 0.0;
	FLOAT temp_r3 = 0.0;
	FLOAT temp_i0 = 0.0;
	FLOAT temp_i1 = 0.0;
	FLOAT temp_i2 = 0.0;
	FLOAT temp_i3 = 0.0;


	for ( i=0; i< 2*n; i+=2 )
	{
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
		temp_r0 += a0[i]*x[i]   - a0[i+1]*x[i+1];		
		temp_i0 += a0[i]*x[i+1] + a0[i+1]*x[i];		
		temp_r1 += a1[i]*x[i]   - a1[i+1]*x[i+1];		
		temp_i1 += a1[i]*x[i+1] + a1[i+1]*x[i];		
		temp_r2 += a2[i]*x[i]   - a2[i+1]*x[i+1];		
		temp_i2 += a2[i]*x[i+1] + a2[i+1]*x[i];		
		temp_r3 += a3[i]*x[i]   - a3[i+1]*x[i+1];		
		temp_i3 += a3[i]*x[i+1] + a3[i+1]*x[i];		
#else
		temp_r0 += a0[i]*x[i]   + a0[i+1]*x[i+1];		
		temp_i0 += a0[i]*x[i+1] - a0[i+1]*x[i];		
		temp_r1 += a1[i]*x[i]   + a1[i+1]*x[i+1];		
		temp_i1 += a1[i]*x[i+1] - a1[i+1]*x[i];		
		temp_r2 += a2[i]*x[i]   + a2[i+1]*x[i+1];		
		temp_i2 += a2[i]*x[i+1] - a2[i+1]*x[i];		
		temp_r3 += a3[i]*x[i]   + a3[i+1]*x[i+1];		
		temp_i3 += a3[i]*x[i+1] - a3[i+1]*x[i];		
#endif
	}
	y[0] = temp_r0;
	y[1] = temp_i0;
	y[2] = temp_r1;
	y[3] = temp_i1;
	y[4] = temp_r2;
	y[5] = temp_i2;
	y[6] = temp_r3;
	y[7] = temp_i3;
}
	
#endif

static void zgemv_kernel_16x1(BLASLONG n, FLOAT *ap, FLOAT *x, FLOAT *y)
{
	BLASLONG i;
	FLOAT *a0;
	a0 = ap;
	FLOAT temp_r = 0.0;
	FLOAT temp_i = 0.0;

	for ( i=0; i< 2*n; i+=2 )
	{
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
		temp_r += a0[i]*x[i]   - a0[i+1]*x[i+1];		
		temp_i += a0[i]*x[i+1] + a0[i+1]*x[i];		
#else
		temp_r += a0[i]*x[i]   + a0[i+1]*x[i+1];		
		temp_i += a0[i]*x[i+1] - a0[i+1]*x[i];		
#endif
	}
	*y      = temp_r;
	*(y+1)  = temp_i;
}
	
static void copy_x(BLASLONG n, FLOAT *src, FLOAT *dest, BLASLONG inc_src)
{
        BLASLONG i;
        for ( i=0; i<n; i++ )
        {
                *dest     = *src;
                *(dest+1) = *(src+1);
                dest+=2;
                src += inc_src;
        }
}


int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha_r, FLOAT alpha_i, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *buffer)
{
	BLASLONG i;
	BLASLONG j;
	FLOAT *a_ptr;
	FLOAT *x_ptr;
	FLOAT *y_ptr;
	FLOAT *ap[8];
	BLASLONG n1;
	BLASLONG m1;
	BLASLONG m2;
	BLASLONG n2;
	FLOAT ybuffer[8],*xbuffer;

        inc_x *= 2;
        inc_y *= 2;
        lda   *= 2;

	xbuffer = buffer;
	
	n1 = n / 4 ;
	n2 = n % 4 ;
	
	m1 = m - ( m % 16 );
	m2 = (m % NBMAX) - (m % 16) ;
	

	BLASLONG NB = NBMAX;

	while ( NB == NBMAX )
	{
		
		m1 -= NB;
		if ( m1 < 0)
		{
			if ( m2 == 0 ) break;	
			NB = m2;
		}
		
		y_ptr = y;
		a_ptr = a;
		x_ptr = x;
		copy_x(NB,x_ptr,xbuffer,inc_x);
		for( i = 0; i < n1 ; i++)
		{
			ap[0] = a_ptr;
			ap[1] = a_ptr + lda;
			ap[2] = ap[1] + lda;
			ap[3] = ap[2] + lda;
			zgemv_kernel_16x4(NB,ap,xbuffer,ybuffer);
			a_ptr += 4 * lda;

#if !defined(XCONJ)
			y_ptr[0] += alpha_r * ybuffer[0] - alpha_i * ybuffer[1];
			y_ptr[1] += alpha_r * ybuffer[1] + alpha_i * ybuffer[0];
			y_ptr  += inc_y;
			y_ptr[0] += alpha_r * ybuffer[2] - alpha_i * ybuffer[3];
			y_ptr[1] += alpha_r * ybuffer[3] + alpha_i * ybuffer[2];
			y_ptr  += inc_y;
			y_ptr[0] += alpha_r * ybuffer[4] - alpha_i * ybuffer[5];
			y_ptr[1] += alpha_r * ybuffer[5] + alpha_i * ybuffer[4];
			y_ptr  += inc_y;
			y_ptr[0] += alpha_r * ybuffer[6] - alpha_i * ybuffer[7];
			y_ptr[1] += alpha_r * ybuffer[7] + alpha_i * ybuffer[6];
			y_ptr  += inc_y;
#else
			y_ptr[0] += alpha_r * ybuffer[0] + alpha_i * ybuffer[1];
			y_ptr[1] -= alpha_r * ybuffer[1] - alpha_i * ybuffer[0];
			y_ptr  += inc_y;
			y_ptr[0] += alpha_r * ybuffer[2] + alpha_i * ybuffer[3];
			y_ptr[1] -= alpha_r * ybuffer[3] - alpha_i * ybuffer[2];
			y_ptr  += inc_y;
			y_ptr[0] += alpha_r * ybuffer[4] + alpha_i * ybuffer[5];
			y_ptr[1] -= alpha_r * ybuffer[5] - alpha_i * ybuffer[4];
			y_ptr  += inc_y;
			y_ptr[0] += alpha_r * ybuffer[6] + alpha_i * ybuffer[7];
			y_ptr[1] -= alpha_r * ybuffer[7] - alpha_i * ybuffer[6];
			y_ptr  += inc_y;
#endif
		}

		for( i = 0; i < n2 ; i++)
		{
			zgemv_kernel_16x1(NB,a_ptr,xbuffer,ybuffer);
			a_ptr += 1 * lda;

#if !defined(XCONJ)
			y_ptr[0] += alpha_r * ybuffer[0] - alpha_i * ybuffer[1];
			y_ptr[1] += alpha_r * ybuffer[1] + alpha_i * ybuffer[0];
			y_ptr  += inc_y;
#else
			y_ptr[0] += alpha_r * ybuffer[0] + alpha_i * ybuffer[1];
			y_ptr[1] -= alpha_r * ybuffer[1] - alpha_i * ybuffer[0];
			y_ptr  += inc_y;
#endif

		}
		a += 2* NB;
		x += NB * inc_x;	
	}

	BLASLONG m3 = m % 16;
	if ( m3 == 0 ) return(0);

	x_ptr = x;
	copy_x(m3,x_ptr,xbuffer,inc_x);
	j=0;
	a_ptr = a;
	y_ptr = y;
	while ( j < n)
	{
		FLOAT temp_r = 0.0;
		FLOAT temp_i = 0.0;
		for( i = 0; i < m3*2; i+=2 )
		{
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
			temp_r += a_ptr[i] * xbuffer[i]   - a_ptr[i+1] * xbuffer[i+1];
			temp_i += a_ptr[i] * xbuffer[i+1] + a_ptr[i+1] * xbuffer[i];
#else
			temp_r += a_ptr[i] * xbuffer[i]   + a_ptr[i+1] * xbuffer[i+1];
			temp_i += a_ptr[i] * xbuffer[i+1] - a_ptr[i+1] * xbuffer[i];
#endif
		}
		a_ptr += lda;

#if !defined(XCONJ) 
                y_ptr[0] += alpha_r * temp_r - alpha_i * temp_i;
                y_ptr[1] += alpha_r * temp_i + alpha_i * temp_r;
#else
                y_ptr[0] += alpha_r * temp_r + alpha_i * temp_i;
                y_ptr[1] -= alpha_r * temp_i - alpha_i * temp_r;
#endif

		y_ptr += inc_y;
		j++;
	}
	return(0);
}


