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

#include <stdio.h>
#include <stdlib.h>
#ifdef __CYGWIN32__
#include <sys/time.h>
#endif
#include "common.h"

#define SINGLE_EPS 1e-04
#define DOUBLE_EPS 1e-13

int assert_dbl_near(double exp, double real, double tol) {
    double diff = exp - real;
    double absdiff = diff;
    /* avoid using fabs and linking with a math lib */
    if(diff < 0) {
      absdiff *= -1;
    }
    if (absdiff > tol) {
        return 0;
    }
    return 1;
}

#ifdef COMPLEX
int zgemv_n_c(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha_r, FLOAT alpha_i, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y)
{
	BLASLONG i;
	BLASLONG ix,iy;
	BLASLONG j;
	FLOAT *a_ptr;
	FLOAT temp_r,temp_i;
	BLASLONG inc_x2,inc_y2;
	BLASLONG lda2;
	BLASLONG i2;

	lda2 = 2*lda;

	ix = 0;
	a_ptr = a;

	if ( inc_x == 1 && inc_y == 1 )
	{

	   for (j=0; j<n; j++)
	   {

#if !defined(XCONJ)
		temp_r = alpha_r * x[ix]   - alpha_i * x[ix+1];
		temp_i = alpha_r * x[ix+1] + alpha_i * x[ix];
#else
		temp_r = alpha_r * x[ix]   + alpha_i * x[ix+1];
		temp_i = alpha_r * x[ix+1] - alpha_i * x[ix];
#endif
		iy = 0;
		i2=0;

		for (i=0; i<m; i++)
		{
#if !defined(CONJ)

#if !defined(XCONJ)
			printf("\nParO: %f %f %f %f\n", a_ptr[i2], a_ptr[i2+1], temp_r, temp_i);
			y[iy]   += temp_r * a_ptr[i2]   - temp_i * a_ptr[i2+1];
			y[iy+1] += temp_r * a_ptr[i2+1] + temp_i * a_ptr[i2];
#else
			y[iy]   += temp_r * a_ptr[i2]   + temp_i * a_ptr[i2+1];
			y[iy+1] += temp_r * a_ptr[i2+1] - temp_i * a_ptr[i2];
#endif

#else

#if !defined(XCONJ)
			y[iy]   += temp_r * a_ptr[i2]   + temp_i * a_ptr[i2+1];
			y[iy+1] -= temp_r * a_ptr[i2+1] - temp_i * a_ptr[i2];
#else
			y[iy]   += temp_r * a_ptr[i2]   - temp_i * a_ptr[i2+1];
			y[iy+1] -= temp_r * a_ptr[i2+1] + temp_i * a_ptr[i2];
#endif

#endif
			i2 += 2;
			iy += 2;
		}
		a_ptr += lda2;
		ix    += 2;
	   }

	   return(0);

	}


	inc_x2 = 2 * inc_x;
	inc_y2 = 2 * inc_y;

	for (j=0; j<n; j++)
	{

#if !defined(XCONJ)
		temp_r = alpha_r * x[ix]   - alpha_i * x[ix+1];
		temp_i = alpha_r * x[ix+1] + alpha_i * x[ix];
#else
		temp_r = alpha_r * x[ix]   + alpha_i * x[ix+1];
		temp_i = alpha_r * x[ix+1] - alpha_i * x[ix];
#endif
		iy = 0;
		i2=0;

		for (i=0; i<m; i++)
		{
#if !defined(CONJ)

#if !defined(XCONJ)
			y[iy]   += temp_r * a_ptr[i2]   - temp_i * a_ptr[i2+1];
			y[iy+1] += temp_r * a_ptr[i2+1] + temp_i * a_ptr[i2];
#else
			y[iy]   += temp_r * a_ptr[i2]   + temp_i * a_ptr[i2+1];
			y[iy+1] += temp_r * a_ptr[i2+1] - temp_i * a_ptr[i2];
#endif

#else

#if !defined(XCONJ)
			y[iy]   += temp_r * a_ptr[i2]   + temp_i * a_ptr[i2+1];
			y[iy+1] -= temp_r * a_ptr[i2+1] - temp_i * a_ptr[i2];
#else
			y[iy]   += temp_r * a_ptr[i2]   - temp_i * a_ptr[i2+1];
			y[iy+1] -= temp_r * a_ptr[i2+1] + temp_i * a_ptr[i2];
#endif

#endif
			i2 += 2;
			iy += inc_y2;
		}
		a_ptr += lda2;
		ix    += inc_x2;
	}


	return(0);
}
int zgemv_t_c(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha_r, FLOAT alpha_i, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y)
{
	BLASLONG i;
	BLASLONG ix,iy;
	BLASLONG j;
	FLOAT *a_ptr;
	FLOAT temp_r,temp_i;
	BLASLONG inc_x2,inc_y2;
	BLASLONG lda2;
	BLASLONG i2;

	lda2 = 2*lda;

	iy = 0;
	a_ptr = a;

	if ( inc_x == 1 && inc_y == 1 )
	{

	   for (j=0; j<n; j++)
	   {
		temp_r = 0.0;
		temp_i = 0.0;
		ix = 0;
		i2=0;

		for (i=0; i<m; i++)
		{

#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
			temp_r += a_ptr[i2] * x[ix]   - a_ptr[i2+1] * x[ix+1];
			temp_i += a_ptr[i2] * x[ix+1] + a_ptr[i2+1] * x[ix];
#else
			temp_r += a_ptr[i2] * x[ix]   + a_ptr[i2+1] * x[ix+1];
			temp_i += a_ptr[i2] * x[ix+1] - a_ptr[i2+1] * x[ix];
#endif

			i2 += 2;
			ix += 2;
		}

#if !defined(XCONJ)
		y[iy]   += alpha_r * temp_r - alpha_i * temp_i;
		y[iy+1] += alpha_r * temp_i + alpha_i * temp_r;
#else
		y[iy]   += alpha_r * temp_r + alpha_i * temp_i;
		y[iy+1] -= alpha_r * temp_i - alpha_i * temp_r;
#endif

		a_ptr += lda2;
		iy    += 2;
	   }

	   return(0);

	}


	inc_x2 = 2 * inc_x;
	inc_y2 = 2 * inc_y;

	for (j=0; j<n; j++)
	{
		temp_r = 0.0;
		temp_i = 0.0;
		ix = 0;
		i2=0;

		for (i=0; i<m; i++)
		{

#if ( !defined(CONJ) && !defined(XCONJ) ) || ( defined(CONJ) && defined(XCONJ) )
			temp_r += a_ptr[i2] * x[ix]   - a_ptr[i2+1] * x[ix+1];
			temp_i += a_ptr[i2] * x[ix+1] + a_ptr[i2+1] * x[ix];
#else
			temp_r += a_ptr[i2] * x[ix]   + a_ptr[i2+1] * x[ix+1];
			temp_i += a_ptr[i2] * x[ix+1] - a_ptr[i2+1] * x[ix];
#endif

			i2 += 2;
			ix += inc_x2;
		}

#if !defined(XCONJ)
		y[iy]   += alpha_r * temp_r - alpha_i * temp_i;
		y[iy+1] += alpha_r * temp_i + alpha_i * temp_r;
#else
		y[iy]   += alpha_r * temp_r + alpha_i * temp_i;
		y[iy+1] -= alpha_r * temp_i - alpha_i * temp_r;
#endif

		a_ptr += lda2;
		iy    += inc_y2;
	}

	return(0);

}
#else
int gemv_n_c(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y)
{
	BLASLONG i;
	BLASLONG ix,iy;
	BLASLONG j;
	FLOAT *a_ptr;
	FLOAT temp;

	ix = 0;
	a_ptr = a;

	for (j=0; j<n; j++)
	{
		temp = alpha * x[ix];
		iy = 0;
		for (i=0; i<m; i++)
		{
			y[iy] += temp * a_ptr[i];
			iy += inc_y;
		}
		a_ptr += lda;
		ix    += inc_x;
	}

	return(0);
}
int gemv_t_c(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y)
{
	BLASLONG i;
	BLASLONG ix,iy;
	BLASLONG j;
	FLOAT *a_ptr;
	FLOAT temp;

	iy = 0;
	a_ptr = a;

	for (j=0; j<n; j++)
	{
		temp = 0.0;
		ix = 0;
		for (i=0; i<m; i++)
		{
			temp += a_ptr[i] * x[ix];
			ix    += inc_x;
		}
		y[iy] += alpha * temp;
		iy += inc_y;
		a_ptr += lda;
	}
  
	return(0);
}
#endif

#undef GEMV
#ifndef COMPLEX
#ifdef DOUBLE
#define GEMV   BLASFUNC(dgemv)
#else
#define GEMV   BLASFUNC(sgemv)
#endif
#else
#ifdef DOUBLE
#define GEMV   BLASFUNC(zgemv)
#else
#define GEMV   BLASFUNC(cgemv)
#endif
#endif

#if defined(__WIN32__) || defined(__WIN64__)

#ifndef DELTA_EPOCH_IN_MICROSECS
#define DELTA_EPOCH_IN_MICROSECS 11644473600000000ULL
#endif

int gettimeofday(struct timeval *tv, void *tz){

  FILETIME ft;
  unsigned __int64 tmpres = 0;
  static int tzflag;

  if (NULL != tv)
    {
      GetSystemTimeAsFileTime(&ft);

      tmpres |= ft.dwHighDateTime;
      tmpres <<= 32;
      tmpres |= ft.dwLowDateTime;

      /*converting file time to unix epoch*/
      tmpres /= 10;  /*convert into microseconds*/
      tmpres -= DELTA_EPOCH_IN_MICROSECS;
      tv->tv_sec = (long)(tmpres / 1000000UL);
      tv->tv_usec = (long)(tmpres % 1000000UL);
    }

  return 0;
}

#endif

#if !defined(__WIN32__) && !defined(__WIN64__) && !defined(__CYGWIN32__) && 0

static void *huge_malloc(BLASLONG size){
  int shmid;
  void *address;

#ifndef SHM_HUGETLB
#define SHM_HUGETLB 04000
#endif

  if ((shmid =shmget(IPC_PRIVATE,
		     (size + HUGE_PAGESIZE) & ~(HUGE_PAGESIZE - 1),
		     SHM_HUGETLB | IPC_CREAT |0600)) < 0) {
    printf( "Memory allocation failed(shmget).\n");
    exit(1);
  }

  address = shmat(shmid, NULL, SHM_RND);

  if ((BLASLONG)address == -1){
    printf( "Memory allocation failed(shmat).\n");
    exit(1);
  }

  shmctl(shmid, IPC_RMID, 0);

  return address;
}

#define malloc huge_malloc

#endif

int main(int argc, char *argv[]){

  FLOAT *a, *x, *y, *y_c;
  FLOAT alpha[] = {1.0, 1.0};
  FLOAT beta [] = {1.0, 1.0};
  char trans='N';
  blasint m, i, j;
  blasint inc_x=1,inc_y=1;
  blasint n=0;
  int has_param_n = 0;
  int has_param_m = 0;
  int loops = 1;
  int l;
  char *p;

  int from =   1;
  int to   = 200;
  int step =   1;

  struct timeval start, stop;
  double time1,timeg,timeg_c;

  blasint iy;
  int test = 1;

  argc--;argv++;

  if (argc > 0) { from     = atol(*argv);		argc--; argv++;}
  if (argc > 0) { to       = MAX(atol(*argv), from);	argc--; argv++;}
  if (argc > 0) { step     = atol(*argv);		argc--; argv++;}


  int tomax = to;

  if ((p = getenv("OPENBLAS_LOOPS")))  loops = atoi(p);
  if ((p = getenv("OPENBLAS_INCX")))   inc_x = atoi(p);
  if ((p = getenv("OPENBLAS_INCY")))   inc_y = atoi(p);
  if ((p = getenv("OPENBLAS_TRANS")))  trans=*p;
  if ((p = getenv("OPENBLAS_PARAM_N"))) {
	  n = atoi(p);
	  if ((n>0)) has_param_n = 1;
  	  if ( n > tomax ) tomax = n;
  }
  if ( has_param_n == 0 )
  	if ((p = getenv("OPENBLAS_PARAM_M"))) {
		  m = atoi(p);
		  if ((m>0)) has_param_m = 1;
  	  	  if ( m > tomax ) tomax = m;
  	}



  fprintf(stderr, "From : %3d  To : %3d Step = %3d Trans = '%c' Inc_x = %d Inc_y = %d Loops = %d\n", from, to, step,trans,inc_x,inc_y,loops);

  if (( a = (FLOAT *)malloc(sizeof(FLOAT) * tomax * tomax * COMPSIZE)) == NULL){
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }

  if (( x = (FLOAT *)malloc(sizeof(FLOAT) * tomax * abs(inc_x) * COMPSIZE)) == NULL){
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }

  if (( y = (FLOAT *)malloc(sizeof(FLOAT) * tomax * abs(inc_y) * COMPSIZE)) == NULL){
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }

  if (( y_c = (FLOAT *)malloc(sizeof(FLOAT) * tomax * abs(inc_y) * COMPSIZE)) == NULL){
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }

#ifdef linux
  srandom(getpid());
#endif

  fprintf(stderr, "    SIZE            Flops           Time          CTime        Test\n");

  if (has_param_m == 0)
  {

  	for(m = from; m <= to; m += step)
  	{
   		timeg=0;
      timeg_c=0;
   		if ( has_param_n == 0 ) n = m;
   		fprintf(stderr, " %6dx%d :", (int)m,(int)n);
   		for(j = 0; j < m; j++){
      			for(i = 0; i < n * COMPSIZE; i++){
				a[i + j * m * COMPSIZE] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
      			}
   		}

    		for (l=0; l<loops; l++)
    		{

   			for(i = 0; i < n * COMPSIZE * abs(inc_x); i++){
				x[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
   			}

   			for(i = 0; i < n * COMPSIZE * abs(inc_y); i++){
				y[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
        y_c[i]= y[i];
   			}
    			gettimeofday( &start, (struct timezone *)0);
    			GEMV (&trans, &m, &n, alpha, a, &m, x, &inc_x, beta, y, &inc_y );
    			gettimeofday( &stop, (struct timezone *)0);
    			time1 = (double)(stop.tv_sec - start.tv_sec) + (double)((stop.tv_usec - start.tv_usec)) * 1.e-6;
			    timeg += time1;

          gettimeofday( &start, (struct timezone *)0);
#ifdef COMPLEX
          if (trans == 'N')
            zgemv_n_c(m, n, 0, alpha[0], alpha[1], a, m, x, inc_x, y_c, inc_y);
          else
            zgemv_t_c(m, n, 0, alpha[0], alpha[1], a, m, x, inc_x, y_c, inc_y);
#else
          if (trans == 'N')
            gemv_n_c(m, n, 0, *alpha, a, m, x, inc_x, y_c, inc_y);
          else
            gemv_t_c(m, n, 0, *alpha, a, m, x, inc_x, y_c, inc_y);
#endif
          gettimeofday( &stop, (struct timezone *)0);
          time1 = (double)(stop.tv_sec - start.tv_sec) + (double)((stop.tv_usec - start.tv_usec)) * 1.e-6;
          timeg_c += time1;

          iy = 0;
#ifdef COMPLEX
          for (i = 0; i < m * 2; i++)
#else
          for (i = 0; i < m; i++)
#endif
          {
            test &= assert_dbl_near(y[iy], y_c[iy], SINGLE_EPS);
            iy += inc_y;
          }

    		}

    		timeg /= loops;
        timeg_c /= loops;

    		fprintf(stderr, "%10.2f MFlops %10.6f sec %10.6f sec    %s\n", 2. * (double)m / timeg * 1.e-6, timeg, timeg_c, test ? "PASS" : "FAILD");

  	}
  }
  else
  {

  	for(n = from; n <= to; n += step)
  	{
   		timeg=0;
      timeg_c=0;
   		fprintf(stderr, " %6dx%d :", (int)m,(int)n);
   		for(j = 0; j < m; j++){
      			for(i = 0; i < n * COMPSIZE; i++){
				a[i + j * m * COMPSIZE] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
      			}
   		}

    		for (l=0; l<loops; l++)
    		{

   			for(i = 0; i < n * COMPSIZE * abs(inc_x); i++){
				x[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
   			}

   			for(i = 0; i < n * COMPSIZE * abs(inc_y); i++){
				y[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
        y_c[i]= y[i];
   			}
    			gettimeofday( &start, (struct timezone *)0);
    			GEMV (&trans, &m, &n, alpha, a, &m, x, &inc_x, beta, y, &inc_y );
    			gettimeofday( &stop, (struct timezone *)0);
    			time1 = (double)(stop.tv_sec - start.tv_sec) + (double)((stop.tv_usec - start.tv_usec)) * 1.e-6;
			    timeg += time1;

          gettimeofday( &start, (struct timezone *)0);
#ifdef COMPLEX
          if (trans == 'N')
            zgemv_n_c(m, n, 0, alpha[0], alpha[1], a, m, x, inc_x, y_c, inc_y);
          else
            zgemv_t_c(m, n, 0, alpha[0], alpha[1], a, m, x, inc_x, y_c, inc_y);
#else
          if (trans == 'N')
            gemv_n_c(m, n, 0, *alpha, a, m, x, inc_x, y_c, inc_y);
          else
            gemv_t_c(m, n, 0, *alpha, a, m, x, inc_x, y_c, inc_y);
#endif
          gettimeofday( &stop, (struct timezone *)0);
          time1 = (double)(stop.tv_sec - start.tv_sec) + (double)((stop.tv_usec - start.tv_usec)) * 1.e-6;
          timeg_c += time1;

          iy = 0;
#ifdef COMPLEX
          for (i = 0; i < m * 2; i++)
#else
          for (i = 0; i < m; i++)
#endif
          {
            test &= assert_dbl_near(y[iy], y_c[iy], SINGLE_EPS);
            iy += inc_y;
          }

    		}

    		timeg /= loops;
        timeg_c /= loops;

    		fprintf(stderr, "%10.2f MFlops %10.6f sec %10.6f sec    %s\n", 2. * (double)m / timeg * 1.e-6, timeg, timeg_c, test ? "PASS" : "FAILD");

  	}
  }

  return 0;
}

// void main(int argc, char *argv[]) __attribute__((weak, alias("MAIN__")));
