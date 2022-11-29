/***************************************************************************
Copyright (c) 2017, The OpenBLAS Project
Copyright (c) 2022, Arm Ltd
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

// Some compilers will report feature support for SVE without the appropriate
// header available
#ifdef HAVE_SVE
#if defined __has_include 
#if __has_include(<arm_sve.h>) && __ARM_FEATURE_SVE
#define USE_SVE
#endif 
#endif
#endif

#ifdef USE_SVE
#include "dot_kernel_sve.c"
#endif
#include "dot_kernel_asimd.c"

#if defined(SMP)
extern int blas_level1_thread_with_return_value(int mode, BLASLONG m, BLASLONG n,
	BLASLONG k, void *alpha, void *a, BLASLONG lda, void *b, BLASLONG ldb,
	void *c, BLASLONG ldc, int (*function)(), int nthreads);
#endif

static RETURN_TYPE dot_compute(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y)
{
	RETURN_TYPE  dot = 0.0 ;

	if ( n <= 0 ) return dot;

#ifdef USE_SVE
	if (inc_x == 1 && inc_y == 1) {
		return dot_kernel_sve(n, x, y);
	}
#endif

	return dot_kernel_asimd(n, x, inc_x, y, inc_y);
}

#if defined(SMP)
static int dot_thread_function(BLASLONG n, BLASLONG dummy0,
	BLASLONG dummy1, FLOAT dummy2, FLOAT *x, BLASLONG inc_x, FLOAT *y,
	BLASLONG inc_y, FLOAT *result, BLASLONG dummy3)
{
	*(RETURN_TYPE *)result = dot_compute(n, x, inc_x, y, inc_y);

	return 0;
}
#endif

RETURN_TYPE CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y)
{
#if defined(SMP)
	int nthreads;
	FLOAT dummy_alpha;
#endif
	RETURN_TYPE dot = 0.0;

#if defined(SMP)
	if (inc_x == 0 || inc_y == 0 || n <= 10000)
		nthreads = 1;
	else
		nthreads = num_cpu_avail(1);

	if (nthreads == 1) {
		dot = dot_compute(n, x, inc_x, y, inc_y);
	} else {
		int mode, i;
		char result[MAX_CPU_NUMBER * sizeof(double) * 2];
		RETURN_TYPE *ptr;

#if !defined(DOUBLE)
		mode = BLAS_SINGLE  | BLAS_REAL;
#else
		mode = BLAS_DOUBLE  | BLAS_REAL;
#endif

		blas_level1_thread_with_return_value(mode, n, 0, 0, &dummy_alpha,
				   x, inc_x, y, inc_y, result, 0,
				   ( void *)dot_thread_function, nthreads);

		ptr = (RETURN_TYPE *)result;
		for (i = 0; i < nthreads; i++) {
			dot = dot + (*ptr);
			ptr = (RETURN_TYPE *)(((char *)ptr) + sizeof(double) * 2);
		}
	}
#else
	dot = dot_compute(n, x, inc_x, y, inc_y);
#endif

	return dot;
}
