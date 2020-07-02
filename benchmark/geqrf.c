/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/*   1. Redistributions of source code must retain the above         */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer.                                                  */
/*                                                                   */
/*   2. Redistributions in binary form must reproduce the above      */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer in the documentation and/or other materials       */
/*      provided with the distribution.                              */
/*                                                                   */
/*    THIS  SOFTWARE IS PROVIDED  BY THE  UNIVERSITY OF  TEXAS AT    */
/*    AUSTIN  ``AS IS''  AND ANY  EXPRESS OR  IMPLIED WARRANTIES,    */
/*    INCLUDING, BUT  NOT LIMITED  TO, THE IMPLIED  WARRANTIES OF    */
/*    MERCHANTABILITY  AND FITNESS FOR  A PARTICULAR  PURPOSE ARE    */
/*    DISCLAIMED.  IN  NO EVENT SHALL THE UNIVERSITY  OF TEXAS AT    */
/*    AUSTIN OR CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT,    */
/*    INCIDENTAL,  SPECIAL, EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES    */
/*    (INCLUDING, BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE    */
/*    GOODS  OR  SERVICES; LOSS  OF  USE,  DATA,  OR PROFITS;  OR    */
/*    BUSINESS INTERRUPTION) HOWEVER CAUSED  AND ON ANY THEORY OF    */
/*    LIABILITY, WHETHER  IN CONTRACT, STRICT  LIABILITY, OR TORT    */
/*    (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY WAY OUT    */
/*    OF  THE  USE OF  THIS  SOFTWARE,  EVEN  IF ADVISED  OF  THE    */
/*    POSSIBILITY OF SUCH DAMAGE.                                    */
/*                                                                   */
/* The views and conclusions contained in the software and           */
/* documentation are those of the authors and should not be          */
/* interpreted as representing official policies, either expressed   */
/* or implied, of The University of Texas at Austin.                 */
/*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#ifdef __CYGWIN32__
#include <sys/time.h>
#endif
#include "common.h"

#undef GEQRF

#ifndef COMPLEX
#ifdef XDOUBLE
#define GEQRF BLASFUNC(qgeqrf)
#elif defined(DOUBLE)
#define GEQRF BLASFUNC(dgeqrf)
#else
#define GEQRF BLASFUNC(sgeqrf)
#endif
#else
#ifdef XDOUBLE
#define GEQRF BLASFUNC(xgeqrf)
#elif defined(DOUBLE)
#define GEQRF BLASFUNC(zgeqrf)
#else
#define GEQRF BLASFUNC(cgeqrf)
#endif
#endif

extern void GEQRF(blasint *m, blasint *n, FLOAT *a, blasint *lda, FLOAT *tau, FLOAT *work, blasint *lwork, blasint *info);

#if defined(__WIN32__) || defined(__WIN64__)

#ifndef DELTA_EPOCH_IN_MICROSECS
#define DELTA_EPOCH_IN_MICROSECS 11644473600000000ULL
#endif

int gettimeofday(struct timeval *tv, void *tz)
{

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
    tmpres /= 10; /*convert into microseconds*/
    tmpres -= DELTA_EPOCH_IN_MICROSECS;
    tv->tv_sec = (long)(tmpres / 1000000UL);
    tv->tv_usec = (long)(tmpres % 1000000UL);
  }

  return 0;
}

#endif

#if !defined(__WIN32__) && !defined(__WIN64__) && !defined(__CYGWIN32__) && 0

static void *huge_malloc(BLASLONG size)
{
  int shmid;
  void *address;

#ifndef SHM_HUGETLB
#define SHM_HUGETLB 04000
#endif

  if ((shmid = shmget(IPC_PRIVATE,
                      (size + HUGE_PAGESIZE) & ~(HUGE_PAGESIZE - 1),
                      SHM_HUGETLB | IPC_CREAT | 0600)) < 0)
  {
    printf("Memory allocation failed(shmget).\n");
    exit(1);
  }

  address = shmat(shmid, NULL, SHM_RND);

  if ((BLASLONG)address == -1)
  {
    printf("Memory allocation failed(shmat).\n");
    exit(1);
  }

  shmctl(shmid, IPC_RMID, 0);

  return address;
}

#define malloc huge_malloc

#endif

int main(int argc, char *argv[])
{

  FLOAT *a, *work;
  FLOAT wkopt[4];
  FLOAT *tau;
  blasint m, n, i = 0, j, info, lwork;
  int loops = 1;
  int l;
  int has_param_m = 0;
  int has_param_n = 0;
  char *p;
#ifndef COMPLEX
  int addfac = 1;
  int mulfac = 1;
#else
  int addfac = 2;
  int mulfac = 6;
#endif
  int from = 1;
  int to = 200;
  int step = 1;
  double flops = 1.0, mults = 0.0, adds = 0.0;

  if ((p = getenv("OPENBLAS_LOOPS")))
    loops = atoi(p);

  struct timeval start, stop;
  double time1, timeg;

  argc--;
  argv++;

  if (argc > 0)
  {
    from = atol(*argv);
    argc--;
    argv++;
  }
  if (argc > 0)
  {
    to = MAX(atol(*argv), from);
    argc--;
    argv++;
  }
  if (argc > 0)
  {
    step = atol(*argv);
    argc--;
    argv++;
  }

  if ((p = getenv("OPENBLAS_PARAM_M")))
  {
    m = atoi(p);
    has_param_m = 1;
  }
  else
  {
    m = to;
  }
  if ((p = getenv("OPENBLAS_PARAM_N")))
  {
    n = atoi(p);
    has_param_n = 1;
  }
  else
  {
    n = to;
  }

  fprintf(stderr, "From : %3d  To : %3d Step = %3d\n", from, to, step);

  if ((a = (FLOAT *)malloc(sizeof(FLOAT) * m * n * COMPSIZE)) == NULL)
  {
    fprintf(stderr, "Out of Memory!!\n");
    exit(1);
  }
  to = MAX(m, n);
  if ((tau = (FLOAT *)malloc(sizeof(FLOAT) * to * COMPSIZE)) == NULL)
  {
    fprintf(stderr, "Out of Memory!!\n");
    exit(1);
  }

  for (j = 0; j < m; j++)
  {
    for (i = 0; i < n * COMPSIZE; i++)
    {
      a[(long)i + (long)j * (long)n * COMPSIZE] = ((FLOAT)rand() / (FLOAT)RAND_MAX) - 0.5;
    }
  }

  lwork = -1;

  GEQRF(&m, &n, a, &m, tau, wkopt, &lwork, &info);

  lwork = (blasint)wkopt[0];
  if ((work = (FLOAT *)malloc(sizeof(FLOAT) * lwork * COMPSIZE)) == NULL)
  {
    fprintf(stderr, "Out of Memory!!\n");
    exit(1);
  }

#ifdef linux
  srandom(getpid());
#endif
  fprintf(stderr, "   SIZE           FLops           Time          Lwork\n");
  for (i = from; i <= to; i += step)
  {
    if (!has_param_m)
    {
      m = i;
    }
    if (!has_param_n)
    {
      n = i;
    }

    fprintf(stderr, " %6d : ", (int)m);

    gettimeofday(&start, (struct timezone *)0);
    for (l = 0; l < loops; l++)
    {
      lwork = -1;
      GEQRF(&m, &n, a, &m, tau, wkopt, &lwork, &info);

      lwork = (blasint)wkopt[0];
      GEQRF(&m, &n, a, &m, tau, work, &lwork, &info);
    }

    gettimeofday(&stop, (struct timezone *)0);

    if (info)
    {
      fprintf(stderr, "the %d-th parameter had an illegal value .. \n", info);
      exit(1);
    }

    time1 = (double)(stop.tv_sec - start.tv_sec) + (double)((stop.tv_usec - start.tv_usec)) * 1.e-6;
    timeg = time1 / loops;

    if (m >= n)
    {
      mults = (double)n * (((23. / 6.) + (double)m + (double)n / 2.) + (double)n * ((double)m - (double)n / 3.));
      adds = (double)n * ((5. / 6.) + (double)n * (1. / 2. + ((double)m - (double)n / 3.)));
    }
    else
    {
      mults = (double)m * (((23. / 6.) + 2. * (double)n - (double)m / 2.) + (double)m * ((double)n - (double)m / 3.));
      adds = (double)m * ((5. / 6.) + (double)n - (double)m / 2. + (double)m * ((double)n - (double)m / 3.));
    }
    flops = (double)mulfac * mults + (double)addfac * adds;
    fprintf(stderr,
            " %10.2f MFlops : %10.6f Sec : %d\n",
            flops / timeg * 1.e-6, timeg, lwork);
  }

  return 0;
}

// void main(int argc, char *argv[]) __attribute__((weak, alias("MAIN__")));
