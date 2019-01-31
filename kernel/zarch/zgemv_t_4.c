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

#define NBMAX 1024

static void zgemv_kernel_4x4(BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y, FLOAT *alpha)
{
    __asm__ volatile (
		"vzero  %%v16                      \n\t"
		"vzero  %%v17                      \n\t"
		"vzero  %%v18                      \n\t"
		"vzero  %%v19                      \n\t"
        "xgr   %%r1,%%r1                   \n\t"
        "srlg  %%r0,%0,1                   \n\t"
        "0:                                \n\t"
        "pfd 1,1024(%%r1,%1)               \n\t"
        "pfd 1,1024(%%r1,%2)               \n\t"
        "pfd 1,1024(%%r1,%3)               \n\t"
        "pfd 1,1024(%%r1,%4)               \n\t"
		"pfd 1,1024(%%r1,%5)               \n\t"

		"vl     %%v20,0(%%r1,%5)           \n\t"
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
        "vleg   %%v21,8(%%r1,%5),0         \n\t"
        "wflcdb %%v21,%%v21                \n\t"
        "vleg   %%v21,0(%%r1,%5),1         \n\t"
#else
        "vleg   %%v21,0(%%r1,%5),1         \n\t"
        "vflcdb %%v21,%%v21                \n\t"
        "vleg   %%v21,8(%%r1,%5),0         \n\t"
#endif

        "vlrepg %%v24,0(%%r1,%1)           \n\t"
        "vlrepg %%v25,8(%%r1,%1)           \n\t"
		"vlrepg %%v26,0(%%r1,%2)           \n\t"
        "vlrepg %%v27,8(%%r1,%2)           \n\t"
        
        "vfmadb   %%v16,%%v24,%%v20,%%v16  \n\t"
        "vfmadb   %%v16,%%v25,%%v21,%%v16  \n\t"
        "vfmadb   %%v17,%%v26,%%v20,%%v17  \n\t"
        "vfmadb   %%v17,%%v27,%%v21,%%v17  \n\t"

        "vlrepg %%v28,0(%%r1,%3)           \n\t"
		"vlrepg %%v29,8(%%r1,%3)           \n\t"
        "vlrepg %%v30,0(%%r1,%4)           \n\t"
        "vlrepg %%v31,8(%%r1,%4)           \n\t"
        
        "vfmadb   %%v18,%%v28,%%v20,%%v18  \n\t"
        "vfmadb   %%v18,%%v29,%%v21,%%v18  \n\t"
        "vfmadb   %%v19,%%v30,%%v20,%%v19  \n\t"
        "vfmadb   %%v19,%%v31,%%v21,%%v19  \n\t"

		"vl     %%v22,16(%%r1,%5)          \n\t"
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
		"vleg   %%v23,24(%%r1,%5),0        \n\t"
        "wflcdb %%v23,%%v23                \n\t"
        "vleg   %%v23,16(%%r1,%5),1        \n\t"
#else
		"vleg   %%v23,16(%%r1,%5),1        \n\t"
        "vflcdb %%v23,%%v23                \n\t"
        "vleg   %%v23,24(%%r1,%5),0        \n\t"
#endif

        "vlrepg %%v24,16(%%r1,%1)          \n\t"
        "vlrepg %%v25,24(%%r1,%1)          \n\t"
		"vlrepg %%v26,16(%%r1,%2)          \n\t"
        "vlrepg %%v27,24(%%r1,%2)          \n\t"
        
        "vfmadb   %%v16,%%v24,%%v22,%%v16  \n\t"
        "vfmadb   %%v16,%%v25,%%v23,%%v16  \n\t"
        "vfmadb   %%v17,%%v26,%%v22,%%v17  \n\t"
        "vfmadb   %%v17,%%v27,%%v23,%%v17  \n\t"

        "vlrepg %%v28,16(%%r1,%3)          \n\t"
		"vlrepg %%v29,24(%%r1,%3)          \n\t"
        "vlrepg %%v30,16(%%r1,%4)          \n\t"
        "vlrepg %%v31,24(%%r1,%4)          \n\t"
        
        "vfmadb   %%v18,%%v28,%%v22,%%v18  \n\t"
        "vfmadb   %%v18,%%v29,%%v23,%%v18  \n\t"
        "vfmadb   %%v19,%%v30,%%v22,%%v19  \n\t"
        "vfmadb   %%v19,%%v31,%%v23,%%v19  \n\t"

        "agfi   %%r1,32                    \n\t"
        "brctg  %%r0,0b                    \n\t"

		"vpdi %%v20,%%v16,%%v16,4          \n\t"
        "vpdi %%v21,%%v17,%%v17,4          \n\t"
        "vpdi %%v22,%%v18,%%v18,4          \n\t"
        "vpdi %%v23,%%v19,%%v19,4          \n\t"
#if !defined(XCONJ)
		"vlrepg %%v24,0(%7)                \n\t"
		"vleg   %%v25,8(%7),0              \n\t"
        "wflcdb %%v25,%%v25                \n\t"
        "vleg   %%v25,8(%7),1              \n\t"
#else
		"vleg   %%v24,0(%7),1              \n\t"
        "vflcdb %%v24,%%v24                \n\t"
        "vleg   %%v24,0(%7),0              \n\t"
		"vlrepg %%v25,8(%7)                \n\t"
#endif
		"vl  %%v26,0(%6)                   \n\t"
		"vl  %%v27,16(%6)                  \n\t"
		"vl  %%v28,32(%6)                  \n\t"
		"vl  %%v29,48(%6)                  \n\t"
		"vfmadb   %%v26,%%v16,%%v24,%%v26  \n\t"
        "vfmadb   %%v26,%%v20,%%v25,%%v26  \n\t"
		"vfmadb   %%v27,%%v17,%%v24,%%v27  \n\t"
        "vfmadb   %%v27,%%v21,%%v25,%%v27  \n\t"
		"vfmadb   %%v28,%%v18,%%v24,%%v28  \n\t"
        "vfmadb   %%v28,%%v22,%%v25,%%v28  \n\t"
		"vfmadb   %%v29,%%v19,%%v24,%%v29  \n\t"
        "vfmadb   %%v29,%%v23,%%v25,%%v29  \n\t"
		"vst  %%v26,0(%6)                  \n\t"
		"vst  %%v27,16(%6)                 \n\t"
		"vst  %%v28,32(%6)                 \n\t"
		"vst  %%v29,48(%6)                     "
        :
        :"r"(n),"ZR"((const FLOAT (*)[n * 2])ap[0]),"ZR"((const FLOAT (*)[n * 2])ap[1]),"ZR"((const FLOAT (*)[n * 2])ap[2]),"ZR"((const FLOAT (*)[n * 2])ap[3]),"ZR"((const FLOAT (*)[n * 2])x),"ZQ"((FLOAT (*)[8])y),"ZQ"((const FLOAT (*)[2])alpha)
        :"memory","cc","r0","r1","v16","v17","v18","v19","v20","v21","v22","v23","v24","v25","v26","v27","v28","v29","v30","v31"
    );
}

static void zgemv_kernel_4x2(BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y, FLOAT *alpha)
{
    __asm__ volatile (
		"vzero  %%v16                      \n\t"
		"vzero  %%v17                      \n\t"
        "xgr   %%r1,%%r1                   \n\t"
        "srlg  %%r0,%0,1                   \n\t"
        "0:                                \n\t"
        "pfd 1,1024(%%r1,%1)               \n\t"
        "pfd 1,1024(%%r1,%2)               \n\t"
        "pfd 1,1024(%%r1,%3)               \n\t"

		"vl     %%v18,0(%%r1,%3)           \n\t"
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
        "vleg   %%v19,8(%%r1,%3),0         \n\t"
        "wflcdb %%v19,%%v19                \n\t"
        "vleg   %%v19,0(%%r1,%3),1         \n\t"
#else
        "vleg   %%v19,0(%%r1,%3),1         \n\t"
        "vflcdb %%v19,%%v19                \n\t"
        "vleg   %%v19,8(%%r1,%3),0         \n\t"
#endif

        "vlrepg %%v20,0(%%r1,%1)           \n\t"
        "vlrepg %%v21,8(%%r1,%1)           \n\t"
		"vlrepg %%v22,0(%%r1,%2)           \n\t"
        "vlrepg %%v23,8(%%r1,%2)           \n\t"
        
        "vfmadb   %%v16,%%v20,%%v18,%%v16  \n\t"
        "vfmadb   %%v16,%%v21,%%v19,%%v16  \n\t"
        "vfmadb   %%v17,%%v22,%%v18,%%v17  \n\t"
        "vfmadb   %%v17,%%v23,%%v19,%%v17  \n\t"

		"vl     %%v18,16(%%r1,%3)           \n\t"
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
        "vleg   %%v19,24(%%r1,%3),0         \n\t"
        "wflcdb %%v19,%%v19                \n\t"
        "vleg   %%v19,16(%%r1,%3),1         \n\t"
#else
        "vleg   %%v19,16(%%r1,%3),1         \n\t"
        "vflcdb %%v19,%%v19                \n\t"
        "vleg   %%v19,24(%%r1,%3),0         \n\t"
#endif

        "vlrepg %%v20,16(%%r1,%1)           \n\t"
        "vlrepg %%v21,24(%%r1,%1)           \n\t"
		"vlrepg %%v22,16(%%r1,%2)           \n\t"
        "vlrepg %%v23,24(%%r1,%2)           \n\t"
        
        "vfmadb   %%v16,%%v20,%%v18,%%v16  \n\t"
        "vfmadb   %%v16,%%v21,%%v19,%%v16  \n\t"
        "vfmadb   %%v17,%%v22,%%v18,%%v17  \n\t"
        "vfmadb   %%v17,%%v23,%%v19,%%v17  \n\t"

        "agfi   %%r1,32                    \n\t"
        "brctg  %%r0,0b                    \n\t"

		"vpdi %%v18,%%v16,%%v16,4          \n\t"
        "vpdi %%v19,%%v17,%%v17,4          \n\t"
#if !defined(XCONJ)
		"vlrepg %%v20,0(%5)                \n\t"
		"vleg   %%v21,8(%5),0              \n\t"
        "wflcdb %%v21,%%v21                \n\t"
        "vleg   %%v21,8(%5),1              \n\t"
#else
		"vleg   %%v20,0(%5),1              \n\t"
        "vflcdb %%v20,%%v20                \n\t"
        "vleg   %%v20,0(%5),0              \n\t"
		"vlrepg %%v21,8(%5)                \n\t"
#endif
		"vl  %%v22,0(%4)                   \n\t"
		"vl  %%v23,16(%4)                  \n\t"
		"vfmadb   %%v22,%%v16,%%v20,%%v22  \n\t"
        "vfmadb   %%v22,%%v18,%%v21,%%v22  \n\t"
		"vfmadb   %%v23,%%v17,%%v20,%%v23  \n\t"
        "vfmadb   %%v23,%%v19,%%v21,%%v23  \n\t"
		"vst  %%v22,0(%4)                  \n\t"
		"vst  %%v23,16(%4)                 \n\t"
        :
        :"r"(n),"ZR"((const FLOAT (*)[n * 2])ap[0]),"ZR"((const FLOAT (*)[n * 2])ap[1]),"ZR"((const FLOAT (*)[n * 2])x),"ZQ"((FLOAT (*)[4])y),"ZQ"((const FLOAT (*)[2])alpha)
        :"memory","cc","r0","r1","v16","v17","v18","v19","v20","v21","v22","v23"
    );
}

static void zgemv_kernel_4x1(BLASLONG n, FLOAT *ap, FLOAT *x, FLOAT *y, FLOAT *alpha)
{
    __asm__ volatile (
		"vzero  %%v16                      \n\t"
        "xgr   %%r1,%%r1                   \n\t"
        "srlg  %%r0,%0,1                   \n\t"
        "0:                                \n\t"
        "pfd 1,1024(%%r1,%1)               \n\t"
        "pfd 1,1024(%%r1,%2)               \n\t"

		"vl     %%v17,0(%%r1,%2)           \n\t"
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
        "vleg   %%v18,8(%%r1,%2),0         \n\t"
        "wflcdb %%v18,%%v18                \n\t"
        "vleg   %%v18,0(%%r1,%2),1         \n\t"
#else
        "vleg   %%v18,0(%%r1,%2),1         \n\t"
        "vflcdb %%v18,%%v18                \n\t"
        "vleg   %%v18,8(%%r1,%2),0         \n\t"
#endif

        "vlrepg %%v19,0(%%r1,%1)           \n\t"
        "vlrepg %%v20,8(%%r1,%1)           \n\t"
        
        "vfmadb   %%v16,%%v19,%%v17,%%v16  \n\t"
        "vfmadb   %%v16,%%v20,%%v18,%%v16  \n\t"

		"vl     %%v17,16(%%r1,%2)           \n\t"
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
        "vleg   %%v18,24(%%r1,%2),0         \n\t"
        "wflcdb %%v18,%%v18                \n\t"
        "vleg   %%v18,16(%%r1,%2),1         \n\t"
#else
        "vleg   %%v18,16(%%r1,%2),1         \n\t"
        "vflcdb %%v18,%%v18                \n\t"
        "vleg   %%v18,24(%%r1,%2),0         \n\t"
#endif

        "vlrepg %%v19,16(%%r1,%1)           \n\t"
        "vlrepg %%v20,24(%%r1,%1)           \n\t"
        
        "vfmadb   %%v16,%%v19,%%v17,%%v16  \n\t"
        "vfmadb   %%v16,%%v20,%%v18,%%v16  \n\t"

        "agfi   %%r1,32                    \n\t"
        "brctg  %%r0,0b                    \n\t"

		"vpdi %%v17,%%v16,%%v16,4          \n\t"
#if !defined(XCONJ)
		"vlrepg %%v18,0(%4)                \n\t"
		"vleg   %%v19,8(%4),0              \n\t"
        "wflcdb %%v19,%%v19                \n\t"
        "vleg   %%v19,8(%4),1              \n\t"
#else
		"vleg   %%v18,0(%4),1              \n\t"
        "vflcdb %%v18,%%v18                \n\t"
        "vleg   %%v18,0(%4),0              \n\t"
		"vlrepg %%v19,8(%4)                \n\t"
#endif
		"vl  %%v20,0(%3)                   \n\t"
		"vfmadb   %%v20,%%v16,%%v18,%%v20  \n\t"
        "vfmadb   %%v20,%%v17,%%v19,%%v20  \n\t"
		"vst  %%v20,0(%3)                  \n\t"
        :
        :"r"(n),"ZR"((const FLOAT (*)[n * 2])ap),"ZR"((const FLOAT (*)[n * 2])x),"ZQ"((FLOAT (*)[2])y),"ZQ"((const FLOAT (*)[2])alpha)
        :"memory","cc","r0","r1","v16","v17","v18","v19","v20"
    );
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
	BLASLONG m3;
	BLASLONG n2;
	BLASLONG lda4;
	FLOAT ybuffer[8],*xbuffer;
	FLOAT alpha[2];

        if ( m < 1 ) return(0);
        if ( n < 1 ) return(0);

        inc_x <<= 1;
        inc_y <<= 1;
        lda   <<= 1;
	lda4    = lda << 2;

	xbuffer = buffer;
	
	n1 = n  >> 2 ;
	n2 = n  &  3 ;
	
	m3 = m & 3 ;
	m1 = m - m3;
	m2 = (m & (NBMAX-1)) - m3 ;
	
	alpha[0] = alpha_r;
	alpha[1] = alpha_i;

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
		ap[0] = a_ptr;
		ap[1] = a_ptr + lda;
		ap[2] = ap[1] + lda;
		ap[3] = ap[2] + lda;
		if ( inc_x != 2 )
			copy_x(NB,x_ptr,xbuffer,inc_x);
		else
			xbuffer = x_ptr;
		
		if ( inc_y == 2 )
		{

			for( i = 0; i < n1 ; i++)
			{
				zgemv_kernel_4x4(NB,ap,xbuffer,y_ptr,alpha);
				ap[0] += lda4;
				ap[1] += lda4;
				ap[2] += lda4;
				ap[3] += lda4;
				a_ptr += lda4;
				y_ptr += 8;
				
			}

			if ( n2 & 2 )
			{
				zgemv_kernel_4x2(NB,ap,xbuffer,y_ptr,alpha);
				a_ptr += lda * 2;
				y_ptr += 4;

			}

			if ( n2 & 1 )
			{
				zgemv_kernel_4x1(NB,a_ptr,xbuffer,y_ptr,alpha);
				/* a_ptr += lda;
				y_ptr += 2; */

			}

		}
		else
		{

			for( i = 0; i < n1 ; i++)
			{
				memset(ybuffer,0,sizeof(ybuffer));
				zgemv_kernel_4x4(NB,ap,xbuffer,ybuffer,alpha);
				ap[0] += lda4;
				ap[1] += lda4;
				ap[2] += lda4;
				ap[3] += lda4;
				a_ptr += lda4;

				y_ptr[0] += ybuffer[0];
				y_ptr[1] += ybuffer[1];
				y_ptr  += inc_y;
				y_ptr[0] += ybuffer[2];
				y_ptr[1] += ybuffer[3];
				y_ptr  += inc_y;
				y_ptr[0] += ybuffer[4];
				y_ptr[1] += ybuffer[5];
				y_ptr  += inc_y;
				y_ptr[0] += ybuffer[6];
				y_ptr[1] += ybuffer[7];
				y_ptr  += inc_y;

			}

			for( i = 0; i < n2 ; i++)
			{
				memset(ybuffer,0,sizeof(ybuffer));
				zgemv_kernel_4x1(NB,a_ptr,xbuffer,ybuffer,alpha);
				a_ptr += lda;
				y_ptr[0] += ybuffer[0];
				y_ptr[1] += ybuffer[1];
				y_ptr  += inc_y;

			}

		}
		a += 2 * NB;
		x += NB * inc_x;	
	}



	if ( m3 == 0 ) return(0);

        x_ptr = x;
        j=0;
        a_ptr = a;
        y_ptr = y;

	if ( m3 == 3 )
	{

                FLOAT temp_r ;
                FLOAT temp_i ;
		FLOAT x0 = x_ptr[0];
		FLOAT x1 = x_ptr[1];
		x_ptr += inc_x;
		FLOAT x2 = x_ptr[0];
		FLOAT x3 = x_ptr[1];
		x_ptr += inc_x;
		FLOAT x4 = x_ptr[0];
		FLOAT x5 = x_ptr[1];
	        while ( j < n)
        	{
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
                       	temp_r  = a_ptr[0] * x0 - a_ptr[1] * x1; 
                       	temp_i  = a_ptr[0] * x1 + a_ptr[1] * x0; 
                       	temp_r += a_ptr[2] * x2 - a_ptr[3] * x3; 
                       	temp_i += a_ptr[2] * x3 + a_ptr[3] * x2; 
                       	temp_r += a_ptr[4] * x4 - a_ptr[5] * x5;
                       	temp_i += a_ptr[4] * x5 + a_ptr[5] * x4;
#else

                       	temp_r  = a_ptr[0] * x0 + a_ptr[1] * x1; 
                       	temp_i  = a_ptr[0] * x1 - a_ptr[1] * x0; 
                       	temp_r += a_ptr[2] * x2 + a_ptr[3] * x3; 
                       	temp_i += a_ptr[2] * x3 - a_ptr[3] * x2; 
                       	temp_r += a_ptr[4] * x4 + a_ptr[5] * x5;
                       	temp_i += a_ptr[4] * x5 - a_ptr[5] * x4;
#endif

#if !defined(XCONJ) 
                	y_ptr[0] += alpha_r * temp_r - alpha_i * temp_i;
                	y_ptr[1] += alpha_r * temp_i + alpha_i * temp_r;
#else
                	y_ptr[0] += alpha_r * temp_r + alpha_i * temp_i;
                	y_ptr[1] -= alpha_r * temp_i - alpha_i * temp_r;
#endif

                	a_ptr += lda;
                	y_ptr += inc_y;
                	j++;
        	}
        	return(0);
	}


	if ( m3 == 2 )
	{

                FLOAT temp_r ;
                FLOAT temp_i ;
                FLOAT temp_r1 ;
                FLOAT temp_i1 ;
		FLOAT x0 = x_ptr[0];
		FLOAT x1 = x_ptr[1];
		x_ptr += inc_x;
		FLOAT x2 = x_ptr[0];
		FLOAT x3 = x_ptr[1];
		FLOAT ar = alpha[0];
		FLOAT ai = alpha[1];

	        while ( j < ( n & -2 ))
        	{
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
                       	temp_r  = a_ptr[0] * x0 - a_ptr[1] * x1; 
                       	temp_i  = a_ptr[0] * x1 + a_ptr[1] * x0; 
                       	temp_r += a_ptr[2] * x2 - a_ptr[3] * x3; 
                       	temp_i += a_ptr[2] * x3 + a_ptr[3] * x2; 
                	a_ptr += lda;
                       	temp_r1  = a_ptr[0] * x0 - a_ptr[1] * x1; 
                       	temp_i1  = a_ptr[0] * x1 + a_ptr[1] * x0; 
                       	temp_r1 += a_ptr[2] * x2 - a_ptr[3] * x3; 
                       	temp_i1 += a_ptr[2] * x3 + a_ptr[3] * x2; 
#else

                       	temp_r  = a_ptr[0] * x0 + a_ptr[1] * x1; 
                       	temp_i  = a_ptr[0] * x1 - a_ptr[1] * x0; 
                       	temp_r += a_ptr[2] * x2 + a_ptr[3] * x3; 
                       	temp_i += a_ptr[2] * x3 - a_ptr[3] * x2; 
                	a_ptr += lda;
                       	temp_r1  = a_ptr[0] * x0 + a_ptr[1] * x1; 
                       	temp_i1  = a_ptr[0] * x1 - a_ptr[1] * x0; 
                       	temp_r1 += a_ptr[2] * x2 + a_ptr[3] * x3; 
                       	temp_i1 += a_ptr[2] * x3 - a_ptr[3] * x2; 
#endif

#if !defined(XCONJ) 
                	y_ptr[0] += ar * temp_r - ai * temp_i;
                	y_ptr[1] += ar * temp_i + ai * temp_r;
                	y_ptr += inc_y;
                	y_ptr[0] += ar * temp_r1 - ai * temp_i1;
                	y_ptr[1] += ar * temp_i1 + ai * temp_r1;
#else
                	y_ptr[0] += ar * temp_r + ai * temp_i;
                	y_ptr[1] -= ar * temp_i - ai * temp_r;
                	y_ptr += inc_y;
                	y_ptr[0] += ar * temp_r1 + ai * temp_i1;
                	y_ptr[1] -= ar * temp_i1 - ai * temp_r1;
#endif

                	a_ptr += lda;
                	y_ptr += inc_y;
                	j+=2;
        	}


	        while ( j < n)
        	{
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
                       	temp_r  = a_ptr[0] * x0 - a_ptr[1] * x1; 
                       	temp_i  = a_ptr[0] * x1 + a_ptr[1] * x0; 
                       	temp_r += a_ptr[2] * x2 - a_ptr[3] * x3; 
                       	temp_i += a_ptr[2] * x3 + a_ptr[3] * x2; 
#else

                       	temp_r  = a_ptr[0] * x0 + a_ptr[1] * x1; 
                       	temp_i  = a_ptr[0] * x1 - a_ptr[1] * x0; 
                       	temp_r += a_ptr[2] * x2 + a_ptr[3] * x3; 
                       	temp_i += a_ptr[2] * x3 - a_ptr[3] * x2; 
#endif

#if !defined(XCONJ) 
                	y_ptr[0] += ar * temp_r - ai * temp_i;
                	y_ptr[1] += ar * temp_i + ai * temp_r;
#else
                	y_ptr[0] += ar * temp_r + ai * temp_i;
                	y_ptr[1] -= ar * temp_i - ai * temp_r;
#endif

                	a_ptr += lda;
                	y_ptr += inc_y;
                	j++;
        	}

        	return(0);
	}


	if ( m3 == 1 )
	{

                FLOAT temp_r ;
                FLOAT temp_i ;
                FLOAT temp_r1 ;
                FLOAT temp_i1 ;
		FLOAT x0 = x_ptr[0];
		FLOAT x1 = x_ptr[1];
		FLOAT ar = alpha[0];
		FLOAT ai = alpha[1];

	        while ( j < ( n & -2 ))
        	{
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
                       	temp_r  = a_ptr[0] * x0 - a_ptr[1] * x1; 
                       	temp_i  = a_ptr[0] * x1 + a_ptr[1] * x0; 
                	a_ptr += lda;
                       	temp_r1  = a_ptr[0] * x0 - a_ptr[1] * x1; 
                       	temp_i1  = a_ptr[0] * x1 + a_ptr[1] * x0; 
#else

                       	temp_r  = a_ptr[0] * x0 + a_ptr[1] * x1; 
                       	temp_i  = a_ptr[0] * x1 - a_ptr[1] * x0; 
                	a_ptr += lda;
                       	temp_r1  = a_ptr[0] * x0 + a_ptr[1] * x1; 
                       	temp_i1  = a_ptr[0] * x1 - a_ptr[1] * x0; 
#endif

#if !defined(XCONJ) 
                	y_ptr[0] += ar * temp_r - ai * temp_i;
                	y_ptr[1] += ar * temp_i + ai * temp_r;
                	y_ptr += inc_y;
                	y_ptr[0] += ar * temp_r1 - ai * temp_i1;
                	y_ptr[1] += ar * temp_i1 + ai * temp_r1;
#else
                	y_ptr[0] += ar * temp_r + ai * temp_i;
                	y_ptr[1] -= ar * temp_i - ai * temp_r;
                	y_ptr += inc_y;
                	y_ptr[0] += ar * temp_r1 + ai * temp_i1;
                	y_ptr[1] -= ar * temp_i1 - ai * temp_r1;
#endif

                	a_ptr += lda;
                	y_ptr += inc_y;
                	j+=2;
        	}

	        while ( j < n)
        	{
#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
                       	temp_r  = a_ptr[0] * x0 - a_ptr[1] * x1; 
                       	temp_i  = a_ptr[0] * x1 + a_ptr[1] * x0; 
#else

                       	temp_r  = a_ptr[0] * x0 + a_ptr[1] * x1; 
                       	temp_i  = a_ptr[0] * x1 - a_ptr[1] * x0; 
#endif

#if !defined(XCONJ) 
                	y_ptr[0] += ar * temp_r - ai * temp_i;
                	y_ptr[1] += ar * temp_i + ai * temp_r;
#else
                	y_ptr[0] += ar * temp_r + ai * temp_i;
                	y_ptr[1] -= ar * temp_i - ai * temp_r;
#endif

                	a_ptr += lda;
                	y_ptr += inc_y;
                	j++;
        	}
        	return(0);
	}

	return(0);
}
