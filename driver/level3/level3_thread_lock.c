/*********************************************************************/
/* Copyright 2026 The OpenBLAS Project.                              */
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
/*    THIS  SOFTWARE IS PROVIDED  BY THE OPENBLAS PROJECT ``AS IS''  */
/*    AND ANY  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     */
/*    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND      */
/*    FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT   */
/*    SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE LIABLE FOR ANY   */
/*    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR           */
/*    CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE OF    */
/*    THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH      */
/*    DAMAGE.                                                        */
/*                                                                   */
/*********************************************************************/

#include "common.h"

#ifdef USE_OPENMP

static omp_lock_t level3_lock, critical_section_lock;
static volatile BLASULONG init_lock = 0;
static _Atomic BLASULONG omp_lock_initialized = 0;
static volatile BLASULONG parallel_section_left = MAX_PARALLEL_NUMBER;

static void blas_level3_thread_lock_init(void)
{
  while (omp_lock_initialized == 0) {
    blas_lock(&init_lock);
    if (omp_lock_initialized == 0) {
      omp_init_lock(&level3_lock);
      omp_init_lock(&critical_section_lock);
      WMB;
      omp_lock_initialized = 1;
    }
    blas_unlock(&init_lock);
  }
}

void blas_level3_thread_enter(void)
{
  blas_level3_thread_lock_init();

  omp_set_lock(&level3_lock);
  omp_set_lock(&critical_section_lock);

  parallel_section_left--;

  if (parallel_section_left != 0)
    omp_unset_lock(&level3_lock);

  omp_unset_lock(&critical_section_lock);
}

void blas_level3_thread_leave(void)
{
  omp_set_lock(&critical_section_lock);

  parallel_section_left++;

  if (parallel_section_left == 1)
    omp_unset_lock(&level3_lock);

  omp_unset_lock(&critical_section_lock);
}

#elif defined(OS_WINDOWS)

static CRITICAL_SECTION level3_lock;
static volatile BLASULONG init_lock = 0;
static volatile BLASULONG level3_lock_initialized = 0;

static void blas_level3_thread_lock_init(void)
{
  while (level3_lock_initialized == 0) {
    blas_lock(&init_lock);
    if (level3_lock_initialized == 0) {
      InitializeCriticalSection((PCRITICAL_SECTION)&level3_lock);
      WMB;
      level3_lock_initialized = 1;
    }
    blas_unlock(&init_lock);
  }
}

void blas_level3_thread_enter(void)
{
  blas_level3_thread_lock_init();
  EnterCriticalSection((PCRITICAL_SECTION)&level3_lock);
}

void blas_level3_thread_leave(void)
{
  LeaveCriticalSection((PCRITICAL_SECTION)&level3_lock);
}

#else

static pthread_mutex_t level3_lock = PTHREAD_MUTEX_INITIALIZER;

void blas_level3_thread_enter(void)
{
  pthread_mutex_lock(&level3_lock);
}

void blas_level3_thread_leave(void)
{
  pthread_mutex_unlock(&level3_lock);
}

#endif
