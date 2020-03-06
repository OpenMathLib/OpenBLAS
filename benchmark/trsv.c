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
#include <time.h>
#include "common.h"


#undef GEMV
#undef TRSV

#ifndef COMPLEX

#ifdef DOUBLE
#define TRSV   BLASFUNC(dtrsv)
#else
#define TRSV   BLASFUNC(strsv)
#endif

#else

#ifdef DOUBLE
#define TRSV   BLASFUNC(ztrsv)
#else
#define TRSV   BLASFUNC(ctrsv)
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

  FLOAT *a, *x;
  blasint n = 0, i, j;
  blasint inc_x=1;
  int loops = 1;
  int l;
  char *p;

  int from =   1;
  int to   = 200;
  int step =   1;

  struct timespec time_start, time_end;
  time_t seconds = 0;

  double time1,timeg;
  long long nanos = 0;

  argc--;argv++;

  if (argc > 0) { from     = atol(*argv);		argc--; argv++;}
  if (argc > 0) { to       = MAX(atol(*argv), from);	argc--; argv++;}
  if (argc > 0) { step     = atol(*argv);		argc--; argv++;}

  char uplo ='L';
  char transa = 'N';
  char diag ='U';

  if ((p = getenv("OPENBLAS_LOOPS")))  loops = atoi(p);
  if ((p = getenv("OPENBLAS_INCX")))   inc_x = atoi(p);
  if ((p = getenv("OPENBLAS_TRANSA")))  transa=*p;
  if ((p = getenv("OPENBLAS_DIAG")))  diag=*p;
  if ((p = getenv("OPENBLAS_UPLO")))  uplo=*p;

  fprintf(stderr, "From : %3d  To : %3d Step = %3d Transa = '%c' Inc_x = %d uplo=%c diag=%c loop = %d\n", from, to, step,transa,inc_x,
          uplo,diag,loops);


#ifdef linux
  srandom(getpid());
#endif

  fprintf(stderr, "   SIZE       Flops\n");
  fprintf(stderr, "============================================\n");

  for(n = from; n <= to; n += step)
  {
      timeg=0;
      if (( a = (FLOAT *)malloc(sizeof(FLOAT) * n * n * COMPSIZE)) == NULL){
          fprintf(stderr,"Out of Memory!!\n");exit(1);
      }

      if (( x = (FLOAT *)malloc(sizeof(FLOAT) * n * abs(inc_x) * COMPSIZE)) == NULL){
          fprintf(stderr,"Out of Memory!!\n");exit(1);
      }

      for(j = 0; j < n; j++){
          for(i = 0; i < n * COMPSIZE; i++){
              a[i + j * n * COMPSIZE] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
          }
      }

      for(i = 0; i < n * COMPSIZE * abs(inc_x); i++){
          x[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
      }

      for(l =0;l< loops;l++){

          clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&time_start);

          TRSV(&uplo,&transa,&diag,&n,a,&n,x,&inc_x);

          clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&time_end);
          nanos = time_end.tv_nsec - time_start.tv_nsec;
          seconds = time_end.tv_sec - time_start.tv_sec;

          time1 = seconds + nanos /1.e9;
          timeg += time1;
      }


      timeg /= loops;
      long long muls = n*(n+1)/2.0;
      long long adds = (n - 1.0)*n/2.0;

      fprintf(stderr, "%10d   %10.2f MFlops %10.6f sec\n", n,(muls+adds) / timeg * 1.e-6, timeg);
      if(a != NULL){
        free(a);
      }

      if( x != NULL){
        free(x);
      }

  }

  return 0;
}

