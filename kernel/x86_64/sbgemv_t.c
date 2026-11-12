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
#include "sbgemv_t_microk_cooperlake.c"
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

/* Computes directly from the bfloat16 data and the strided x/y, without allocating copies of A, x
 * or y (m and n are already swapped: y has m elements, x has n).  Every y[i] is accumulated over
 * j in increasing order, as before; four outputs share each converted x[j]. */
static void sbgemv_kernel_t_generic(BLASLONG m, BLASLONG n, float alpha, bfloat16 *a, BLASLONG lda, bfloat16 *x, BLASLONG incx, float beta, float *y, BLASLONG incy)
{
    BLASLONG i = 0;

    for (; i + 4 <= m; i += 4) {
        bfloat16 *a0 = a + (i + 0) * lda;
        bfloat16 *a1 = a + (i + 1) * lda;
        bfloat16 *a2 = a + (i + 2) * lda;
        bfloat16 *a3 = a + (i + 3) * lda;
        float accum0 = 0.0, accum1 = 0.0, accum2 = 0.0, accum3 = 0.0;

        for (BLASLONG j = 0; j < n; j++) {
            float tmp_x = sbgemv_bf16tos(x[j * incx]);
            accum0 += sbgemv_bf16tos(a0[j]) * tmp_x;
            accum1 += sbgemv_bf16tos(a1[j]) * tmp_x;
            accum2 += sbgemv_bf16tos(a2[j]) * tmp_x;
            accum3 += sbgemv_bf16tos(a3[j]) * tmp_x;
        }
        if (beta == ZERO) {
            y[(i + 0) * incy] = alpha * accum0;
            y[(i + 1) * incy] = alpha * accum1;
            y[(i + 2) * incy] = alpha * accum2;
            y[(i + 3) * incy] = alpha * accum3;
        } else {
            y[(i + 0) * incy] = alpha * accum0 + beta * y[(i + 0) * incy];
            y[(i + 1) * incy] = alpha * accum1 + beta * y[(i + 1) * incy];
            y[(i + 2) * incy] = alpha * accum2 + beta * y[(i + 2) * incy];
            y[(i + 3) * incy] = alpha * accum3 + beta * y[(i + 3) * incy];
        }
    }

    for (; i < m; i++) {
        bfloat16 *a0 = a + i * lda;
        float accum = 0.0;

        for (BLASLONG j = 0; j < n; j++) {
            accum += sbgemv_bf16tos(a0[j]) * sbgemv_bf16tos(x[j * incx]);
        }
        if (beta == ZERO) {
            y[i * incy] = alpha * accum;
        } else {
            y[i * incy] = alpha * accum + beta * y[i * incy];
        }
    }
}

#ifdef HAVE_SBGEMV_T_ACCL_KERNEL
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

    // Switch m and n
    BLASLONG t = m;
    m = n;
    n = t;

#ifdef HAVE_SBGEMV_T_ACCL_KERNEL
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
        sbgemv_kernel_t_generic(m, n, alpha, a, lda, x, incx, beta, y, incy);
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

    sbgemv_kernel_t(m, n, alpha, a, lda, xbuffer_align, beta, ybuffer_align);

    if (incy != 1) {
        fp32_expand_vector(m, ybuffer_align, y, incy);
        ALIGN64_FREE(ybuffer);
    }

    if (incx != 1) {
        ALIGN64_FREE(xbuffer);
    }
#else
    sbgemv_kernel_t_generic(m, n, alpha, a, lda, x, incx, beta, y, incy);
#endif

	return(0);
}
