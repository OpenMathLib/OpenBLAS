/***************************************************************************
Copyright (c) 2026, The OpenBLAS Project
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
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

/*
 * WASM SIMD128 GEMV_N: y += alpha * A * x  (AXPY each column of A into y).
 *
 * Compiled twice: SGEMV_N (-UDOUBLE) and DGEMV_N (-DDOUBLE).
 *
 * Four columns share one streaming of y. Inner loop unrolls four v128 lanes
 * (16 floats / 8 doubles). IEEE mul+add (not relaxed madd).
 * Non-unit inc_y stays scalar (no scatter).
 */

#include "common.h"

#if defined(__wasm_simd128__)
#include <wasm_simd128.h>

#ifdef DOUBLE
#define GEMV_VLEN 2
#define GEMV_SPLAT wasm_f64x2_splat
#define GEMV_MUL wasm_f64x2_mul
#define GEMV_ADD wasm_f64x2_add
#else
#define GEMV_VLEN 4
#define GEMV_SPLAT wasm_f32x4_splat
#define GEMV_MUL wasm_f32x4_mul
#define GEMV_ADD wasm_f32x4_add
#endif

#define GEMV_UNROLL 4
#define GEMV_CHUNK (GEMV_VLEN * GEMV_UNROLL)
#define GEMV_LOAD(p) wasm_v128_load((const void *)(p))
#define GEMV_STORE(p, v) wasm_v128_store((void *)(p), (v))
#define GEMV_MADD(y, a, x) GEMV_ADD((y), GEMV_MUL((a), (x)))

static void gemv_n_axpy4(BLASLONG m, const FLOAT *a0, const FLOAT *a1,
                         const FLOAT *a2, const FLOAT *a3, FLOAT *y, FLOAT t0,
                         FLOAT t1, FLOAT t2, FLOAT t3) {
  const v128_t v0 = GEMV_SPLAT(t0);
  const v128_t v1 = GEMV_SPLAT(t1);
  const v128_t v2 = GEMV_SPLAT(t2);
  const v128_t v3 = GEMV_SPLAT(t3);
  BLASLONG i = 0;
  const BLASLONG n_main = m & ~(BLASLONG)(GEMV_CHUNK - 1);

  for (; i < n_main; i += GEMV_CHUNK) {
    v128_t acc;
    acc = GEMV_LOAD(y + i + 0 * GEMV_VLEN);
    acc = GEMV_MADD(acc, v0, GEMV_LOAD(a0 + i + 0 * GEMV_VLEN));
    acc = GEMV_MADD(acc, v1, GEMV_LOAD(a1 + i + 0 * GEMV_VLEN));
    acc = GEMV_MADD(acc, v2, GEMV_LOAD(a2 + i + 0 * GEMV_VLEN));
    acc = GEMV_MADD(acc, v3, GEMV_LOAD(a3 + i + 0 * GEMV_VLEN));
    GEMV_STORE(y + i + 0 * GEMV_VLEN, acc);

    acc = GEMV_LOAD(y + i + 1 * GEMV_VLEN);
    acc = GEMV_MADD(acc, v0, GEMV_LOAD(a0 + i + 1 * GEMV_VLEN));
    acc = GEMV_MADD(acc, v1, GEMV_LOAD(a1 + i + 1 * GEMV_VLEN));
    acc = GEMV_MADD(acc, v2, GEMV_LOAD(a2 + i + 1 * GEMV_VLEN));
    acc = GEMV_MADD(acc, v3, GEMV_LOAD(a3 + i + 1 * GEMV_VLEN));
    GEMV_STORE(y + i + 1 * GEMV_VLEN, acc);

    acc = GEMV_LOAD(y + i + 2 * GEMV_VLEN);
    acc = GEMV_MADD(acc, v0, GEMV_LOAD(a0 + i + 2 * GEMV_VLEN));
    acc = GEMV_MADD(acc, v1, GEMV_LOAD(a1 + i + 2 * GEMV_VLEN));
    acc = GEMV_MADD(acc, v2, GEMV_LOAD(a2 + i + 2 * GEMV_VLEN));
    acc = GEMV_MADD(acc, v3, GEMV_LOAD(a3 + i + 2 * GEMV_VLEN));
    GEMV_STORE(y + i + 2 * GEMV_VLEN, acc);

    acc = GEMV_LOAD(y + i + 3 * GEMV_VLEN);
    acc = GEMV_MADD(acc, v0, GEMV_LOAD(a0 + i + 3 * GEMV_VLEN));
    acc = GEMV_MADD(acc, v1, GEMV_LOAD(a1 + i + 3 * GEMV_VLEN));
    acc = GEMV_MADD(acc, v2, GEMV_LOAD(a2 + i + 3 * GEMV_VLEN));
    acc = GEMV_MADD(acc, v3, GEMV_LOAD(a3 + i + 3 * GEMV_VLEN));
    GEMV_STORE(y + i + 3 * GEMV_VLEN, acc);
  }
  for (; i + GEMV_VLEN <= m; i += GEMV_VLEN) {
    v128_t acc = GEMV_LOAD(y + i);
    acc = GEMV_MADD(acc, v0, GEMV_LOAD(a0 + i));
    acc = GEMV_MADD(acc, v1, GEMV_LOAD(a1 + i));
    acc = GEMV_MADD(acc, v2, GEMV_LOAD(a2 + i));
    acc = GEMV_MADD(acc, v3, GEMV_LOAD(a3 + i));
    GEMV_STORE(y + i, acc);
  }
  for (; i < m; i++)
    y[i] += t0 * a0[i] + t1 * a1[i] + t2 * a2[i] + t3 * a3[i];
}

static void gemv_n_axpy1(BLASLONG m, const FLOAT *a, FLOAT *y, FLOAT t) {
  const v128_t vt = GEMV_SPLAT(t);
  BLASLONG i = 0;
  const BLASLONG n_main = m & ~(BLASLONG)(GEMV_CHUNK - 1);

  for (; i < n_main; i += GEMV_CHUNK) {
    v128_t acc;
    acc = GEMV_LOAD(y + i + 0 * GEMV_VLEN);
    GEMV_STORE(y + i + 0 * GEMV_VLEN, GEMV_MADD(acc, vt, GEMV_LOAD(a + i + 0 * GEMV_VLEN)));
    acc = GEMV_LOAD(y + i + 1 * GEMV_VLEN);
    GEMV_STORE(y + i + 1 * GEMV_VLEN, GEMV_MADD(acc, vt, GEMV_LOAD(a + i + 1 * GEMV_VLEN)));
    acc = GEMV_LOAD(y + i + 2 * GEMV_VLEN);
    GEMV_STORE(y + i + 2 * GEMV_VLEN, GEMV_MADD(acc, vt, GEMV_LOAD(a + i + 2 * GEMV_VLEN)));
    acc = GEMV_LOAD(y + i + 3 * GEMV_VLEN);
    GEMV_STORE(y + i + 3 * GEMV_VLEN, GEMV_MADD(acc, vt, GEMV_LOAD(a + i + 3 * GEMV_VLEN)));
  }
  for (; i + GEMV_VLEN <= m; i += GEMV_VLEN) {
    v128_t acc = GEMV_LOAD(y + i);
    GEMV_STORE(y + i, GEMV_MADD(acc, vt, GEMV_LOAD(a + i)));
  }
  for (; i < m; i++)
    y[i] += t * a[i];
}

#undef GEMV_MADD
#undef GEMV_STORE
#undef GEMV_LOAD
#undef GEMV_CHUNK
#undef GEMV_UNROLL
#undef GEMV_ADD
#undef GEMV_MUL
#undef GEMV_SPLAT
#undef GEMV_VLEN
#endif

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha, FLOAT *a,
          BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y,
          FLOAT *buffer) {
  BLASLONG i, j, ix, iy;
  FLOAT *a_ptr;
  FLOAT temp;

  (void)dummy1;
  (void)buffer;

  if (m < 1 || n < 1)
    return 0;
  if (alpha == (FLOAT)0.0)
    return 0;

#if defined(__wasm_simd128__)
  if (inc_y == 1) {
    ix = 0;
    j = 0;
    for (; j + 4 <= n; j += 4) {
      gemv_n_axpy4(m, a + j * lda, a + (j + 1) * lda, a + (j + 2) * lda,
                   a + (j + 3) * lda, y, alpha * x[ix],
                   alpha * x[ix + inc_x], alpha * x[ix + 2 * inc_x],
                   alpha * x[ix + 3 * inc_x]);
      ix += 4 * inc_x;
    }
    for (; j < n; j++) {
      gemv_n_axpy1(m, a + j * lda, y, alpha * x[ix]);
      ix += inc_x;
    }
    return 0;
  }
#endif

  ix = 0;
  a_ptr = a;
  for (j = 0; j < n; j++) {
    temp = alpha * x[ix];
    iy = 0;
    for (i = 0; i < m; i++) {
      y[iy] += temp * a_ptr[i];
      iy += inc_y;
    }
    a_ptr += lda;
    ix += inc_x;
  }
  return 0;
}
