/***************************************************************************
Copyright (c) 2020, The OpenBLAS Project
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
#include <altivec.h>

typedef __vector unsigned char  vec_t;
typedef FLOAT v4sf_t __attribute__ ((vector_size (16)));
typedef __vector_pair          __attribute__((aligned(8))) vecp_t;

#include "dgemv_n_microk_power10.c"

#define MMA(X, APTR, ACC) \
        rX = (vec_t *) & X; \
        rowA = *((vecp_t*)((void*)&APTR)); \
        __builtin_mma_xvf64gerpp (ACC, rowA, rX[0]);

#define SAVE(ACC, Z) \
        rowC = (v4sf_t *) &y[Z]; \
        __builtin_mma_disassemble_acc ((void *)result, ACC); \
        result[0][1] = result[1][0]; \
        result[2][1] = result[3][0]; \
        rowC[0] += valpha * result[0]; \
        rowC[1] += valpha * result[2];

void
dgemv_kernel_4x128 (BLASLONG n, FLOAT * a_ptr, BLASLONG lda, FLOAT * xo,
                    FLOAT * y, FLOAT alpha)
{
  BLASLONG i, j, tmp;
  FLOAT *a0 = a_ptr;
  FLOAT *x1 = xo;
  vector double valpha = { alpha, alpha };
  v4sf_t *rowC;
  __vector_quad acc0, acc1, acc2, acc3, acc4, acc5, acc6, acc7;
  v4sf_t result[4];
  vecp_t rowA;
  vec_t *rX;
  tmp = (n / 32) * 32;
  for (i = 0; i < tmp; i += 32)
    {
      xo = x1;
      a0 = a_ptr;
      __builtin_mma_xxsetaccz (&acc0);
      __builtin_mma_xxsetaccz (&acc1);
      __builtin_mma_xxsetaccz (&acc2);
      __builtin_mma_xxsetaccz (&acc3);
      __builtin_mma_xxsetaccz (&acc4);
      __builtin_mma_xxsetaccz (&acc5);
      __builtin_mma_xxsetaccz (&acc6);
      __builtin_mma_xxsetaccz (&acc7);
      for (j = 0; j < 32; j++)
        {
          __builtin_prefetch (xo+j);
          __builtin_prefetch (a0+i+j+lda);
          MMA (xo[j], a0[i + 0 + j * lda], &acc0);
          MMA (xo[j], a0[i + 4 + j * lda], &acc1);
          MMA (xo[j], a0[i + 8 + j * lda], &acc2);
          MMA (xo[j], a0[i + 12 + j * lda], &acc3);
          MMA (xo[j], a0[i + 16 + j * lda], &acc4);
          MMA (xo[j], a0[i + 20 + j * lda], &acc5);
          MMA (xo[j], a0[i + 24 + j * lda], &acc6);
          MMA (xo[j], a0[i + 28 + j * lda], &acc7);
        }
      xo += 32;
      a0 += lda << 5;
      for (j = 0; j < 32; j++)
        {
          __builtin_prefetch (xo+j);
          __builtin_prefetch (a0+i+j+lda);
          MMA (xo[j], a0[i + 0 + j * lda], &acc0);
          MMA (xo[j], a0[i + 4 + j * lda], &acc1);
          MMA (xo[j], a0[i + 8 + j * lda], &acc2);
          MMA (xo[j], a0[i + 12 + j * lda], &acc3);
          MMA (xo[j], a0[i + 16 + j * lda], &acc4);
          MMA (xo[j], a0[i + 20 + j * lda], &acc5);
          MMA (xo[j], a0[i + 24 + j * lda], &acc6);
          MMA (xo[j], a0[i + 28 + j * lda], &acc7);
        }
      xo += 32;
      a0 += lda << 5;
      for (j = 0; j < 32; j++)
        {
          __builtin_prefetch (xo+j);
          __builtin_prefetch (a0+i+j+lda);
          MMA (xo[j], a0[i + 0 + j * lda], &acc0);
          MMA (xo[j], a0[i + 4 + j * lda], &acc1);
          MMA (xo[j], a0[i + 8 + j * lda], &acc2);
          MMA (xo[j], a0[i + 12 + j * lda], &acc3);
          MMA (xo[j], a0[i + 16 + j * lda], &acc4);
          MMA (xo[j], a0[i + 20 + j * lda], &acc5);
          MMA (xo[j], a0[i + 24 + j * lda], &acc6);
          MMA (xo[j], a0[i + 28 + j * lda], &acc7);
        }
      xo += 32;
      a0 += lda << 5;
      for (j = 0; j < 32; j++)
        {
          __builtin_prefetch (xo+j);
          __builtin_prefetch (a0+i+j+lda);
          MMA (xo[j], a0[i + 0 + j * lda], &acc0);
          MMA (xo[j], a0[i + 4 + j * lda], &acc1);
          MMA (xo[j], a0[i + 8 + j * lda], &acc2);
          MMA (xo[j], a0[i + 12 + j * lda], &acc3);
          MMA (xo[j], a0[i + 16 + j * lda], &acc4);
          MMA (xo[j], a0[i + 20 + j * lda], &acc5);
          MMA (xo[j], a0[i + 24 + j * lda], &acc6);
          MMA (xo[j], a0[i + 28 + j * lda], &acc7);
        }
      xo += 32;
      a0 += lda << 5;
      SAVE (&acc0, i + 0);
      SAVE (&acc1, i + 4);
      SAVE (&acc2, i + 8);
      SAVE (&acc3, i + 12);
      SAVE (&acc4, i + 16);
      SAVE (&acc5, i + 20);
      SAVE (&acc6, i + 24);
      SAVE (&acc7, i + 28);

    }
  for (i = tmp; i < n; i += 4)
    {
      xo = x1;
      a0 = a_ptr;
      __builtin_mma_xxsetaccz (&acc0);
      for (j = 0; j < 32; j++)
        {
          __builtin_prefetch (xo+j);
          __builtin_prefetch (a0+i+j+lda);
          MMA (xo[j], a0[i + j * lda], &acc0);
        }
      xo += 32;
      a0 += lda << 5;
      for (j = 0; j < 32; j++)
        {
          __builtin_prefetch (xo+j);
          __builtin_prefetch (a0+i+j+lda);
          MMA (xo[j], a0[i + j * lda], &acc0);
        }
      xo += 32;
      a0 += lda << 5;
      for (j = 0; j < 32; j++)
        {
          __builtin_prefetch (xo+j);
          __builtin_prefetch (a0+i+j+lda);
          MMA (xo[j], a0[i + j * lda], &acc0);
        }
      xo += 32;
      a0 += lda << 5;
      for (j = 0; j < 32; j++)
        {
          __builtin_prefetch (xo+j);
          __builtin_prefetch (a0+i+j+lda);
          MMA (xo[j], a0[i + j * lda], &acc0);
        }
      xo += 32;
      a0 += lda << 5;
      SAVE (&acc0, i);
    }
}


#define NBMAX 4096

#ifndef HAVE_KERNEL_4x4

static void dgemv_kernel_4x4(BLASLONG n, FLOAT *a_ptr, BLASLONG lda, FLOAT *xo, FLOAT *y, FLOAT alpha)
{
	BLASLONG i;
	FLOAT x[4]  __attribute__ ((aligned (16)));;
	FLOAT *a0 = a_ptr;
	FLOAT *a1 = a0 + lda;
	FLOAT *a2 = a1 + lda;
	FLOAT *a3 = a2 + lda;


	for ( i=0; i<4; i++)
		x[i] = xo[i] * alpha;

	for ( i=0; i< n; i+=4 )
	{
		y[i] += a0[i]*x[0] + a1[i]*x[1] + a2[i]*x[2] + a3[i]*x[3];		
		y[i+1] += a0[i+1]*x[0] + a1[i+1]*x[1] + a2[i+1]*x[2] + a3[i+1]*x[3];		
		y[i+2] += a0[i+2]*x[0] + a1[i+2]*x[1] + a2[i+2]*x[2] + a3[i+2]*x[3];		
		y[i+3] += a0[i+3]*x[0] + a1[i+3]*x[1] + a2[i+3]*x[2] + a3[i+3]*x[3];		
	}
}

#endif

#ifndef HAVE_KERNEL_4x2

static void dgemv_kernel_4x2(BLASLONG n, FLOAT *a0, FLOAT *a1, FLOAT *xo, FLOAT *y, FLOAT alpha)
{
	BLASLONG i;
	FLOAT x[4]  __attribute__ ((aligned (16)));;

	for ( i=0; i<2; i++)
		x[i] = xo[i] * alpha;

	for ( i=0; i< n; i+=4 )
	{
		y[i] += a0[i]*x[0] + a1[i]*x[1];		
		y[i+1] += a0[i+1]*x[0] + a1[i+1]*x[1];		
		y[i+2] += a0[i+2]*x[0] + a1[i+2]*x[1];		
		y[i+3] += a0[i+3]*x[0] + a1[i+3]*x[1];		
	}
}


#endif

#ifndef HAVE_KERNEL_4x1

static void dgemv_kernel_4x1(BLASLONG n, FLOAT *a0, FLOAT *xo, FLOAT *y, FLOAT alpha)
{
	BLASLONG i;
	FLOAT x[4]  __attribute__ ((aligned (16)));;

	for ( i=0; i<1; i++)
		x[i] = xo[i] * alpha;

	for ( i=0; i< n; i+=4 )
	{
		y[i] += a0[i]*x[0];		
		y[i+1] += a0[i+1]*x[0];		
		y[i+2] += a0[i+2]*x[0];		
		y[i+3] += a0[i+3]*x[0];		
	}
}


#endif


static void add_y(BLASLONG n, FLOAT *src, FLOAT *dest, BLASLONG inc_dest)
{
	BLASLONG i;
	if ( inc_dest != 1 )
	{
		for ( i=0; i<n; i++ )
		{
			*dest += *src;
			src++;
			dest += inc_dest;
		}
		return;
	}

}

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *buffer)
{

	BLASLONG i;
	FLOAT *a_ptr;
	FLOAT *x_ptr;
	FLOAT *y_ptr;
	BLASLONG n1;
	BLASLONG m1;
	BLASLONG m2;
	BLASLONG m3;
	BLASLONG n2;
	BLASLONG lda4 =  lda << 2;
	BLASLONG lda128 = lda << 7;

	FLOAT xbuffer[8] __attribute__ ((aligned (16)));
	FLOAT *ybuffer;

        if ( m < 1 ) return(0);
        if ( n < 1 ) return(0);

	ybuffer = buffer;
	BLASLONG n128 = n >> 7;
	n1 = (n - (n128 * 128)) >> 2;
	n2 = (n - (n128 * 128)) & 3;

        m3 = m & 3  ;
        m1 = m & -4 ;
        m2 = (m & (NBMAX-1)) - m3 ;

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
		
		if ( inc_y != 1 )
			memset(ybuffer,0,NB*8);
		else
			ybuffer = y_ptr;

		if ( inc_x == 1 )
		{

			for( i = 0; i < n128 ; i++)
			{
				dgemv_kernel_4x128(NB,a_ptr,lda,x_ptr,ybuffer,alpha);
				a_ptr += lda128;
				x_ptr += 128;
			}

			for( i = 0; i < n1 ; i++)
			{
				dgemv_kernel_4x4(NB,a_ptr,lda,x_ptr,ybuffer,alpha);
				a_ptr += lda4;
				x_ptr += 4;	
			}

			if ( n2 & 2 )
			{
				dgemv_kernel_4x2(NB,a_ptr,a_ptr+lda,x_ptr,ybuffer,alpha);
				a_ptr += lda*2;
				x_ptr += 2;	
			}


			if ( n2 & 1 )
			{
				dgemv_kernel_4x1(NB,a_ptr,x_ptr,ybuffer,alpha);
				a_ptr += lda;
				x_ptr += 1;	

			}


		}
		else
		{
			for( i = 0; i < n128 ; i++)
			{
	                        FLOAT xbuffer[128] __attribute__ ((aligned (16)));
				BLASLONG j;
				for ( j = 0; j < 128 ; j++)
				{
					xbuffer[j] = x_ptr[0];
				        x_ptr += inc_x;
				}
				dgemv_kernel_4x128(NB,a_ptr,lda,xbuffer,ybuffer,alpha);
				a_ptr += lda128;
			}

			for( i = 0; i < n1 ; i++)
			{
				xbuffer[0] = x_ptr[0];
				x_ptr += inc_x;	
				xbuffer[1] =  x_ptr[0];
				x_ptr += inc_x;	
				xbuffer[2] =  x_ptr[0];
				x_ptr += inc_x;	
				xbuffer[3] = x_ptr[0];
				x_ptr += inc_x;	
				dgemv_kernel_4x4(NB,a_ptr,lda,xbuffer,ybuffer,alpha);
				a_ptr += lda4;
			}

			for( i = 0; i < n2 ; i++)
			{
				xbuffer[0] = x_ptr[0];
				x_ptr += inc_x;	
				dgemv_kernel_4x1(NB,a_ptr,xbuffer,ybuffer,alpha);
				a_ptr += lda;

			}

		}

		a     += NB;
		if ( inc_y != 1 )
		{
			add_y(NB,ybuffer,y_ptr,inc_y);
			y_ptr += NB * inc_y;
		}
		else
			y_ptr += NB ;

	}

	if ( m3 == 0 ) return(0);

	if ( m3 == 3 )
	{
		a_ptr = a;
		x_ptr = x;
		FLOAT temp0 = 0.0;
		FLOAT temp1 = 0.0;
		FLOAT temp2 = 0.0;
		if ( lda == 3 && inc_x ==1 )
		{

			for( i = 0; i < ( n & -4 ); i+=4 )
			{

				temp0 += a_ptr[0] * x_ptr[0] + a_ptr[3] * x_ptr[1];
				temp1 += a_ptr[1] * x_ptr[0] + a_ptr[4] * x_ptr[1];
				temp2 += a_ptr[2] * x_ptr[0] + a_ptr[5] * x_ptr[1];

				temp0 += a_ptr[6] * x_ptr[2] + a_ptr[9]  * x_ptr[3];
				temp1 += a_ptr[7] * x_ptr[2] + a_ptr[10] * x_ptr[3];
				temp2 += a_ptr[8] * x_ptr[2] + a_ptr[11] * x_ptr[3];

				a_ptr += 12;
				x_ptr += 4;
			}

			for( ; i < n; i++ )
			{
				temp0 += a_ptr[0] * x_ptr[0];
				temp1 += a_ptr[1] * x_ptr[0];
				temp2 += a_ptr[2] * x_ptr[0];
				a_ptr += 3;
				x_ptr ++;
			}

		}
		else
		{

			for( i = 0; i < n; i++ )
			{
				temp0 += a_ptr[0] * x_ptr[0];
				temp1 += a_ptr[1] * x_ptr[0];
				temp2 += a_ptr[2] * x_ptr[0];
				a_ptr += lda;
				x_ptr += inc_x;


			}

		}
		y_ptr[0] += alpha * temp0;
		y_ptr += inc_y;
		y_ptr[0] += alpha * temp1;
		y_ptr += inc_y;
		y_ptr[0] += alpha * temp2;
		return(0);
	}


	if ( m3 == 2 )
	{
		a_ptr = a;
		x_ptr = x;
		FLOAT temp0 = 0.0;
		FLOAT temp1 = 0.0;
		if ( lda == 2 && inc_x ==1 )
		{

			for( i = 0; i < (n & -4) ; i+=4 )
			{
				temp0 += a_ptr[0] * x_ptr[0] + a_ptr[2] * x_ptr[1];
				temp1 += a_ptr[1] * x_ptr[0] + a_ptr[3] * x_ptr[1];
				temp0 += a_ptr[4] * x_ptr[2] + a_ptr[6] * x_ptr[3];
				temp1 += a_ptr[5] * x_ptr[2] + a_ptr[7] * x_ptr[3];
				a_ptr += 8;
				x_ptr += 4;

			}


			for( ; i < n; i++ )
			{
				temp0 += a_ptr[0]   * x_ptr[0];
				temp1 += a_ptr[1]   * x_ptr[0];
				a_ptr += 2;
				x_ptr ++;
			}

		}
		else
		{

			for( i = 0; i < n; i++ )
			{
				temp0 += a_ptr[0] * x_ptr[0];
				temp1 += a_ptr[1] * x_ptr[0];
				a_ptr += lda;
				x_ptr += inc_x;


			}

		}
		y_ptr[0] += alpha * temp0;
		y_ptr += inc_y;
		y_ptr[0] += alpha * temp1;
		return(0);
	}

	if ( m3 == 1 )
	{
		a_ptr = a;
		x_ptr = x;
		FLOAT temp = 0.0;
		if ( lda == 1 && inc_x ==1 )
		{

			for( i = 0; i < (n & -4); i+=4 )
			{
				temp += a_ptr[i] * x_ptr[i] + a_ptr[i+1] * x_ptr[i+1] + a_ptr[i+2] * x_ptr[i+2] + a_ptr[i+3] * x_ptr[i+3];
	
			}

			for( ; i < n; i++ )
			{
				temp += a_ptr[i] * x_ptr[i];
			}

		}
		else
		{

			for( i = 0; i < n; i++ )
			{
				temp += a_ptr[0] * x_ptr[0];
				a_ptr += lda;
				x_ptr += inc_x;
			}

		}
		y_ptr[0] += alpha * temp;
		return(0);
	}


	return(0);
}


