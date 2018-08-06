/***************************************************************************
Copyright (c) 2014, The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above swapright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above swapright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE SWAPRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
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
int zswap_c(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT dummy3, FLOAT dummy4, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *dummy, BLASLONG dummy2)
{
	BLASLONG i=0;
	BLASLONG ix=0,iy=0;
	FLOAT temp[2];
	BLASLONG inc_x2;
	BLASLONG inc_y2;

	if ( n < 0     )  return(0);

	inc_x2 = 2 * inc_x;
	inc_y2 = 2 * inc_y;

	while(i < n)
	{

		temp[0]  = x[ix]   ;
		temp[1]  = x[ix+1] ;
		x[ix]    = y[iy]   ;
		x[ix+1]  = y[iy+1] ;
		y[iy]    = temp[0] ;
		y[iy+1]  = temp[1] ;

		ix += inc_x2 ;
		iy += inc_y2 ;
		i++ ;

	}
	return(0);
}
#else
int swap_c(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT dummy3, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *dummy, BLASLONG dummy2)
{
	BLASLONG i=0;
	BLASLONG ix=0,iy=0;
	FLOAT temp;

	if ( n < 0     )  return(0);

	while(i < n)
	{

		temp  = x[ix] ;
		x[ix] = y[iy] ;
		y[iy] = temp ;

		ix += inc_x ;
		iy += inc_y ;
		i++ ;

	}
	return(0);
}
#endif

#undef SWAP
#ifdef COMPLEX
#ifdef DOUBLE
#define SWAP   BLASFUNC(zswap)
#else
#define SWAP   BLASFUNC(cswap)
#endif
#else
#ifdef DOUBLE
#define SWAP   BLASFUNC(dswap)
#else
#define SWAP   BLASFUNC(sswap)
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

  FLOAT *x, *y, *x_c, *y_c;
  blasint m, i;
  blasint inc_x=1,inc_y=1;
  int loops = 1;
  int l;
  char *p;

  int from =   1;
  int to   = 200;
  int step =   1;

  struct timeval start, stop;
  double time1,timeg,timeg_c;

  blasint ix,iy;
  int test = 1;

  argc--;argv++;

  if (argc > 0) { from     = atol(*argv);		argc--; argv++;}
  if (argc > 0) { to       = MAX(atol(*argv), from);	argc--; argv++;}
  if (argc > 0) { step     = atol(*argv);		argc--; argv++;}

  if ((p = getenv("OPENBLAS_LOOPS")))  loops = atoi(p);
  if ((p = getenv("OPENBLAS_INCX")))   inc_x = atoi(p);
  if ((p = getenv("OPENBLAS_INCY")))   inc_y = atoi(p);

  fprintf(stderr, "From : %3d  To : %3d Step = %3d Inc_x = %d Inc_y = %d Loops = %d\n", from, to, step,inc_x,inc_y,loops);

  if (( x = (FLOAT *)malloc(sizeof(FLOAT) * to * abs(inc_x) * COMPSIZE)) == NULL){
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }

  if (( y = (FLOAT *)malloc(sizeof(FLOAT) * to * abs(inc_y) * COMPSIZE)) == NULL){
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }

  if (( x_c = (FLOAT *)malloc(sizeof(FLOAT) * to * abs(inc_x) * COMPSIZE)) == NULL){
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }

  if (( y_c = (FLOAT *)malloc(sizeof(FLOAT) * to * abs(inc_y) * COMPSIZE)) == NULL){
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }

#ifdef linux
  srandom(getpid());
#endif

  fprintf(stderr, "    SIZE            Flops           Time          CTime        Test\n");

  for(m = from; m <= to; m += step)
  {

   timeg=0;
   timeg_c=0;

   fprintf(stderr, " %6d :", (int)m);


   for (l=0; l<loops; l++)
   {

   	for(i = 0; i < m * COMPSIZE * abs(inc_x); i++){
			x[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
      x_c[i] = x[i];
   	}

   	for(i = 0; i < m * COMPSIZE * abs(inc_y); i++){
			y[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
      y_c[i] = y[i];
   	}
    	gettimeofday( &start, (struct timezone *)0);
    	SWAP (&m, x, &inc_x, y, &inc_y );
    	gettimeofday( &stop, (struct timezone *)0);
    	time1 = (double)(stop.tv_sec - start.tv_sec) + (double)((stop.tv_usec - start.tv_usec)) * 1.e-6;
	    timeg += time1;

      gettimeofday( &start, (struct timezone *)0);
#ifdef COMPLEX
      zswap_c(m, 0, 0, 0, 0, x_c, inc_x, y_c, inc_y, NULL, 0);
#else
      swap_c(m, 0, 0, 0, x_c, inc_x, y_c, inc_y, NULL, 0);
#endif
    	gettimeofday( &stop, (struct timezone *)0);
    	time1 = (double)(stop.tv_sec - start.tv_sec) + (double)((stop.tv_usec - start.tv_usec)) * 1.e-6;
	    timeg_c += time1;

      ix = 0;
      iy = 0;
#ifdef COMPLEX
      for (i = 0; i < m * 2; i++)
#else
      for (i = 0; i < m; i++)
#endif
      {
        test &= assert_dbl_near(x[ix], x_c[ix], SINGLE_EPS);
        test &= assert_dbl_near(y[ix], y_c[ix], SINGLE_EPS);
        ix += inc_x;
        iy += inc_y;
      }
    }

    timeg /= loops;
    timeg_c /= loops;

#ifdef COMPLEX
    fprintf(stderr, "%10.2f MFlops %10.6f sec %10.6f sec    %s\n", 6. * (double)m / timeg * 1.e-6, timeg, timeg_c, test ? "PASS" : "FAILD");
#else
    fprintf(stderr, "%10.2f MFlops %10.6f sec %10.6f sec    %s\n", 1. * (double)m / timeg * 1.e-6, timeg, timeg_c, test ? "PASS" : "FAILD");
#endif

  }

  return 0;
}

// void main(int argc, char *argv[]) __attribute__((weak, alias("MAIN__")));
