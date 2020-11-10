#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#ifdef __CYGWIN32__
#include <sys/time.h>
#endif
#include "common.h"

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

#if defined(__WIN32__) || defined(__WIN64__) || !defined(_POSIX_TIMERS)
  struct timeval start, stop;
#else
  struct timespec start = { 0, 0 }, stop = { 0, 0 };
#endif

double getsec()
{
#if defined(__WIN32__) || defined(__WIN64__) || !defined(_POSIX_TIMERS)
    return (double)(stop.tv_sec - start.tv_sec) + (double)((stop.tv_usec - start.tv_usec)) * 1.e-6;
#else
    return (double)(stop.tv_sec - start.tv_sec) + (double)((stop.tv_nsec - start.tv_nsec)) * 1.e-9;
#endif
}

void begin() {
#if defined(__WIN32__) || defined(__WIN64__) || !defined(_POSIX_TIMERS)
    gettimeofday( &start, (struct timezone *)0);
#else
    clock_gettime(CLOCK_REALTIME, &start);
#endif
}

void end() {
#if defined(__WIN32__) || defined(__WIN64__) || !defined(_POSIX_TIMERS)
    gettimeofday( &stop, (struct timezone *)0);
#else
    clock_gettime(CLOCK_REALTIME, &stop);
#endif
}