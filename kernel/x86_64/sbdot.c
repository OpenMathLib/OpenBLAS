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

#if defined(COOPERLAKE) || defined(SAPPHIRERAPIDS)
#include "sbdot_microk_cooperlake.c"
#endif

/* the inputs are converted to fp32 in blocks of SBDOT_BLOCK elements held in
   a stack array; these SDOTU_K calls are too short to start threads of their
   own (kernel/arm64/dot.c threads above 10000 elements) */
#define SBDOT_BLOCK 1024

#ifdef HAVE_SBDOT_ACCL_KERNEL
/* keep the generic path out of line, so that the AVX512-BF16 path does not
   set up its stack frame (HAVE_SBDOT_ACCL_KERNEL implies GCC, clang or
   NVHPC) */
#define SBDOT_GENERIC_NOINLINE __attribute__((noinline))
#else
#define SBDOT_GENERIC_NOINLINE
#endif

/* same results as SBF16TOS_K (bf16to.c): denormals become signed zeros, NaNs
   are forced to quiet NaNs, all other values are converted exactly */
static inline void sbdot_bf16tos(BLASLONG n, bfloat16 *in, BLASLONG inc, float *out)
{
    union { uint32_t u; float f; } v;
    uint32_t b, e;
    BLASLONG i;

    for (i = 0; i < n; i++) {
        b = in[i * inc];
        e = b & 0x7f80u;
        v.u = b << 16;
        v.u = (e == 0) ? (v.u & 0x80000000u) : v.u;
        v.u = (e == 0x7f80u && (b & 0x007fu)) ? (v.u | 0x00400000u) : v.u;
        out[i] = v.f;
    }
}

static SBDOT_GENERIC_NOINLINE float sbdot_generic(BLASLONG n, bfloat16 *x, BLASLONG inc_x, bfloat16 *y, BLASLONG inc_y)
{
    float d = 0.0;
    /* both blocks share one stack array, aligned to 64 bytes by hand (16
       spare floats): this file is also built by MSVC, which does not accept
       __attribute__((aligned)) */
    float buffer[2 * SBDOT_BLOCK + 16];
    float *x_fp32 = (float *)(((BLASLONG)buffer + 63) & ~63);
    float *y_fp32 = x_fp32 + SBDOT_BLOCK;
    float r;
    BLASLONG i, m;

    for (i = 0; i < n; i += SBDOT_BLOCK) {
        m = n - i;
        if (m > SBDOT_BLOCK) m = SBDOT_BLOCK;

        sbdot_bf16tos(m, x + i * inc_x, inc_x, x_fp32);
        sbdot_bf16tos(m, y + i * inc_y, inc_y, y_fp32);

        r = SDOTU_K(m, x_fp32, 1, y_fp32, 1);
        /* the first block is taken as is, so that n <= SBDOT_BLOCK returns
           the SDOTU_K result bit for bit */
        d = (i == 0) ? r : d + r;
    }

    return d;
}

static float sbdot_compute(BLASLONG n, bfloat16 *x, BLASLONG inc_x, bfloat16 *y, BLASLONG inc_y)
{
#ifdef HAVE_SBDOT_ACCL_KERNEL
    if ((inc_x == 1) && (inc_y == 1)) {
        return sbdot_accl_kernel(n, x, y);
    }
#endif

    return sbdot_generic(n, x, inc_x, y, inc_y);
}

#if defined(SMP)
static int sbdot_thread_func(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, bfloat16 dummy2,
		            bfloat16 *x, BLASLONG inc_x, bfloat16 *y, BLASLONG inc_y,
			    float *result, BLASLONG dummy3)
{
    *(float *)result = sbdot_compute(n, x, inc_x, y, inc_y);
    return 0;
}

extern int blas_level1_thread_with_return_value(int mode, BLASLONG m, BLASLONG n, BLASLONG k, void *alpha,
                            void *a, BLASLONG lda, void *b, BLASLONG ldb, void *c, BLASLONG ldc,
			    int (*function)(), int nthreads);
#endif

float CNAME(BLASLONG n, bfloat16 *x, BLASLONG inc_x, bfloat16 *y, BLASLONG inc_y)
{
    float dot_result = 0.0;

    if (n <= 0)  return 0.0;

#if defined(SMP)
    int nthreads;
    int thread_thres = 40960;
    bfloat16 dummy_alpha;
#endif

#if defined(SMP)
    if (inc_x == 0 || inc_y == 0 || n <= thread_thres)
        nthreads = 1;
    else
        nthreads = num_cpu_avail(1);

    int best_threads = (int) (n/(float)thread_thres + 0.5);

    if (best_threads < nthreads) {
        nthreads = best_threads;
    }

    if (nthreads <= 1) {
        dot_result = sbdot_compute(n, x, inc_x, y, inc_y);
    } else {
        char thread_result[MAX_CPU_NUMBER * sizeof(double) * 2];
        int mode = BLAS_BFLOAT16 | BLAS_REAL;
        blas_level1_thread_with_return_value(mode, n, 0, 0, &dummy_alpha,
                                             x, inc_x, y, inc_y, thread_result, 0,
                                             (void *)sbdot_thread_func, nthreads);
        float * ptr = (float *)thread_result;
        for (int i = 0; i < nthreads; i++) {
            dot_result += (*ptr);
            ptr = (float *)(((char *)ptr) + sizeof(double) * 2);
        }
    }
#else
    dot_result = sbdot_compute(n, x, inc_x, y, inc_y);
#endif

    return dot_result;
}
