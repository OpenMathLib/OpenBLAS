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

#if defined (COOPERLAKE) || defined (SAPPHIRERAPIDS)
#include "sbgemv_n_microk_cooperlake.c"
#endif

#define ALIGN64_ALLOC(alloc_size, TYPE, ptr_align, ptr)   \
    ptr = (TYPE *) malloc(sizeof(TYPE)*alloc_size + 63); \
    ptr_align = ((int)(((uintptr_t)ptr & (uintptr_t)0x3F))!=0) ? (TYPE *)((char *)ptr + (64 - (int)((uintptr_t)ptr & (uintptr_t)0x3F))) : ptr

#define ALIGN64_FREE(ptr) \
    free(ptr)

/* bfloat16 -> float exactly as SBF16TOS_K (kernel/x86_64/bf16to.c) converts it:
 * denormals become signed zero, NaNs are made quiet */
static inline float sbgemv_bf16tos(bfloat16 v)
{
    uint32_t u = ((uint32_t)v) << 16;
    float f;

    if ((u & 0x7f800000u) == 0)
        u &= 0x80000000u;
    else if (((u & 0x7f800000u) == 0x7f800000u) && (u & 0x007f0000u))
        u |= 0x00400000u;
    memcpy(&f, &u, sizeof(f));
    return f;
}

#define SBGEMV_N_BLOCK 64

/* Computes directly from the bfloat16 data and the strided x/y, without allocating copies of A, x
 * or y.  Every y[i] is accumulated over the columns in increasing order, as before. */
static void sbgemv_kernel_n_generic(BLASLONG m, BLASLONG n, float alpha, bfloat16 *a, BLASLONG lda, bfloat16 *x, BLASLONG incx, float beta, float *y, BLASLONG incy)
{
    float accum[SBGEMV_N_BLOCK];

    for (BLASLONG i0 = 0; i0 < m; i0 += SBGEMV_N_BLOCK) {
        BLASLONG mb = MIN(m - i0, SBGEMV_N_BLOCK);

        for (BLASLONG i = 0; i < mb; i++) accum[i] = 0.0;

        for (BLASLONG j = 0; j < n; j++) {
            float tmp_x = sbgemv_bf16tos(x[j * incx]);
            bfloat16 *a_j = a + j * lda + i0;
            for (BLASLONG i = 0; i < mb; i++) {
                accum[i] += sbgemv_bf16tos(a_j[i]) * tmp_x;
            }
        }

        for (BLASLONG i = 0; i < mb; i++) {
            if (beta == ZERO) {
                y[(i0 + i) * incy] = alpha * accum[i];
            } else {
                y[(i0 + i) * incy] = alpha * accum[i] + beta * y[(i0 + i) * incy];
            }
        }
    }
}

#ifdef HAVE_SBGEMV_N_ACCL_KERNEL
static void bf16_compress_vector(BLASLONG n, bfloat16 * src, bfloat16 * target, BLASLONG inc)
{
    for(BLASLONG i=0; i<n; i++) {
        target[i] = src[i*inc];
    }
}

static void fp32_compress_vector(BLASLONG n, float * src, float * target, BLASLONG inc)
{
    for(BLASLONG i=0; i<n; i++) {
        target[i] = src[i*inc];
    }
}

static void fp32_expand_vector(BLASLONG n, float * src, float * target, BLASLONG inc)
{
    for(BLASLONG i=0; i<n; i++) {
        target[i*inc] = src[i];
    }
}
#endif

int CNAME(BLASLONG m, BLASLONG n, float alpha, bfloat16 *a, BLASLONG lda, bfloat16 *x, BLASLONG incx, float beta, float * y, BLASLONG incy)
{
    if ( m < 1 || n < 1) return(0);

#ifdef HAVE_SBGEMV_N_ACCL_KERNEL
    bfloat16 * xbuffer_align = x;
    float    * ybuffer_align = y;

    bfloat16 * xbuffer = NULL;
    float    * ybuffer = NULL;

    if (incx != 1) {
        ALIGN64_ALLOC(n, bfloat16, xbuffer_align, xbuffer);
    }

    if (incy != 1) {
        ALIGN64_ALLOC(m, float, ybuffer_align, ybuffer);
    }

    if ((incx != 1 && xbuffer == NULL) || (incy != 1 && ybuffer == NULL)) {
        /* out of memory for the contiguous copies of x/y: use the strided C code */
        ALIGN64_FREE(xbuffer);
        ALIGN64_FREE(ybuffer);
        sbgemv_kernel_n_generic(m, n, alpha, a, lda, x, incx, beta, y, incy);
        return(0);
    }

    if (incx != 1) {
        bf16_compress_vector(n, x, xbuffer_align, incx);
    }

    if (incy != 1) {
        if (beta != ZERO) {
            fp32_compress_vector(m, y, ybuffer_align, incy);
        }
    }

    sbgemv_kernel_n(m, n, alpha, a, lda, xbuffer_align, beta, ybuffer_align);

    if (incy != 1) {
        fp32_expand_vector(m, ybuffer_align, y, incy);
        ALIGN64_FREE(ybuffer);
    }

    if (incx != 1) {
        ALIGN64_FREE(xbuffer);
    }
#else
    sbgemv_kernel_n_generic(m, n, alpha, a, lda, x, incx, beta, y, incy);
#endif

	return(0);
}
