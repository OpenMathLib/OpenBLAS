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
 * WASM SIMD128 GEMV_T: y += alpha * A^T * x  (dot each column of A with x).
 *
 * Compiled twice: SGEMV_T (-UDOUBLE) and DGEMV_T (-DDOUBLE).
 *
 * A single-column reduction is the ~4 GFLOPS SDOT ceiling on wasm128.
 * This kernel keeps several independent column accumulators and loads x
 * once per inner step:
 *
 *   float:  8 columns, two f32x4 steps (8 elements) per inner iter
 *   double: 4 columns, two f64x2 steps (4 elements) per inner iter
 *
 * IEEE mul+add (not relaxed madd). Horizontal add only after the inner
 * loop. Non-unit inc_x stays scalar (no gather).
 */

#include "common.h"

#if defined(__wasm_simd128__)
#include <wasm_simd128.h>

#ifdef DOUBLE
#define GEMV_VLEN 2
#define GEMV_MUL wasm_f64x2_mul
#define GEMV_ADD wasm_f64x2_add
#define GEMV_ZERO wasm_f64x2_const(0.0, 0.0)
#else
#define GEMV_VLEN 4
#define GEMV_MUL wasm_f32x4_mul
#define GEMV_ADD wasm_f32x4_add
#define GEMV_ZERO wasm_f32x4_const(0.0f, 0.0f, 0.0f, 0.0f)
#endif

#define GEMV_LOAD(p) wasm_v128_load((const void *)(p))
#define GEMV_MADD(acc, a, x) GEMV_ADD((acc), GEMV_MUL((a), (x)))
#define GEMV_STEP (2 * GEMV_VLEN)

#ifdef DOUBLE
static inline FLOAT gemv_hsum(v128_t v) {
  return wasm_f64x2_extract_lane(v, 0) + wasm_f64x2_extract_lane(v, 1);
}
#else
static inline FLOAT gemv_hsum(v128_t v) {
  v128_t s = wasm_f32x4_add(v, wasm_i32x4_shuffle(v, v, 2, 3, 0, 1));
  s = wasm_f32x4_add(s, wasm_i32x4_shuffle(s, s, 1, 0, 3, 2));
  return wasm_f32x4_extract_lane(s, 0);
}
#endif

#define GEMV_DOT2(acc, ap, xp0, xp1, i) \
  do { \
    (acc) = GEMV_MADD((acc), GEMV_LOAD((ap) + (i)), (xp0)); \
    (acc) = GEMV_MADD((acc), GEMV_LOAD((ap) + (i) + GEMV_VLEN), (xp1)); \
  } while (0)

static void gemv_t_1col(BLASLONG m, const FLOAT *a, const FLOAT *x, FLOAT *y,
                        FLOAT alpha) {
  v128_t acc0 = GEMV_ZERO;
  v128_t acc1 = GEMV_ZERO;
  BLASLONG i = 0;
  const BLASLONG n2 = m & ~(BLASLONG)(GEMV_STEP - 1);
  for (; i < n2; i += GEMV_STEP) {
    v128_t x0 = GEMV_LOAD(x + i);
    v128_t x1 = GEMV_LOAD(x + i + GEMV_VLEN);
    acc0 = GEMV_MADD(acc0, GEMV_LOAD(a + i), x0);
    acc1 = GEMV_MADD(acc1, GEMV_LOAD(a + i + GEMV_VLEN), x1);
  }
  FLOAT t = gemv_hsum(GEMV_ADD(acc0, acc1));
  for (; i < m; i++)
    t += a[i] * x[i];
  *y += alpha * t;
}

static void gemv_t_4col(BLASLONG m, const FLOAT *a0, const FLOAT *a1,
                        const FLOAT *a2, const FLOAT *a3, const FLOAT *x,
                        FLOAT *y, FLOAT alpha, BLASLONG inc_y) {
  v128_t acc0 = GEMV_ZERO, acc1 = GEMV_ZERO, acc2 = GEMV_ZERO, acc3 = GEMV_ZERO;
  BLASLONG i = 0;
  const BLASLONG n2 = m & ~(BLASLONG)(GEMV_STEP - 1);
  for (; i < n2; i += GEMV_STEP) {
    v128_t x0 = GEMV_LOAD(x + i);
    v128_t x1 = GEMV_LOAD(x + i + GEMV_VLEN);
    GEMV_DOT2(acc0, a0, x0, x1, i);
    GEMV_DOT2(acc1, a1, x0, x1, i);
    GEMV_DOT2(acc2, a2, x0, x1, i);
    GEMV_DOT2(acc3, a3, x0, x1, i);
  }
  FLOAT t0 = gemv_hsum(acc0);
  FLOAT t1 = gemv_hsum(acc1);
  FLOAT t2 = gemv_hsum(acc2);
  FLOAT t3 = gemv_hsum(acc3);
  for (; i < m; i++) {
    t0 += a0[i] * x[i];
    t1 += a1[i] * x[i];
    t2 += a2[i] * x[i];
    t3 += a3[i] * x[i];
  }
  y[0] += alpha * t0;
  y[inc_y] += alpha * t1;
  y[2 * inc_y] += alpha * t2;
  y[3 * inc_y] += alpha * t3;
}

#ifndef DOUBLE
static void gemv_t_8col(BLASLONG m, const FLOAT *a0, const FLOAT *a1,
                        const FLOAT *a2, const FLOAT *a3, const FLOAT *a4,
                        const FLOAT *a5, const FLOAT *a6, const FLOAT *a7,
                        const FLOAT *x, FLOAT *y, FLOAT alpha, BLASLONG inc_y) {
  v128_t acc0 = GEMV_ZERO, acc1 = GEMV_ZERO, acc2 = GEMV_ZERO, acc3 = GEMV_ZERO;
  v128_t acc4 = GEMV_ZERO, acc5 = GEMV_ZERO, acc6 = GEMV_ZERO, acc7 = GEMV_ZERO;
  BLASLONG i = 0;
  const BLASLONG n2 = m & ~(BLASLONG)(GEMV_STEP - 1);
  for (; i < n2; i += GEMV_STEP) {
    v128_t x0 = GEMV_LOAD(x + i);
    v128_t x1 = GEMV_LOAD(x + i + GEMV_VLEN);
    GEMV_DOT2(acc0, a0, x0, x1, i);
    GEMV_DOT2(acc1, a1, x0, x1, i);
    GEMV_DOT2(acc2, a2, x0, x1, i);
    GEMV_DOT2(acc3, a3, x0, x1, i);
    GEMV_DOT2(acc4, a4, x0, x1, i);
    GEMV_DOT2(acc5, a5, x0, x1, i);
    GEMV_DOT2(acc6, a6, x0, x1, i);
    GEMV_DOT2(acc7, a7, x0, x1, i);
  }
  FLOAT t0 = gemv_hsum(acc0);
  FLOAT t1 = gemv_hsum(acc1);
  FLOAT t2 = gemv_hsum(acc2);
  FLOAT t3 = gemv_hsum(acc3);
  FLOAT t4 = gemv_hsum(acc4);
  FLOAT t5 = gemv_hsum(acc5);
  FLOAT t6 = gemv_hsum(acc6);
  FLOAT t7 = gemv_hsum(acc7);
  for (; i < m; i++) {
    t0 += a0[i] * x[i];
    t1 += a1[i] * x[i];
    t2 += a2[i] * x[i];
    t3 += a3[i] * x[i];
    t4 += a4[i] * x[i];
    t5 += a5[i] * x[i];
    t6 += a6[i] * x[i];
    t7 += a7[i] * x[i];
  }
  y[0] += alpha * t0;
  y[inc_y] += alpha * t1;
  y[2 * inc_y] += alpha * t2;
  y[3 * inc_y] += alpha * t3;
  y[4 * inc_y] += alpha * t4;
  y[5 * inc_y] += alpha * t5;
  y[6 * inc_y] += alpha * t6;
  y[7 * inc_y] += alpha * t7;
}
#endif

#undef GEMV_DOT2
#undef GEMV_STEP
#undef GEMV_MADD
#undef GEMV_LOAD
#undef GEMV_ZERO
#undef GEMV_ADD
#undef GEMV_MUL
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
  if (inc_x == 1) {
    j = 0;
    iy = 0;
#ifndef DOUBLE
    for (; j + 8 <= n; j += 8) {
      gemv_t_8col(m, a + j * lda, a + (j + 1) * lda, a + (j + 2) * lda,
                  a + (j + 3) * lda, a + (j + 4) * lda, a + (j + 5) * lda,
                  a + (j + 6) * lda, a + (j + 7) * lda, x, y + iy, alpha,
                  inc_y);
      iy += 8 * inc_y;
    }
#endif
    for (; j + 4 <= n; j += 4) {
      gemv_t_4col(m, a + j * lda, a + (j + 1) * lda, a + (j + 2) * lda,
                  a + (j + 3) * lda, x, y + iy, alpha, inc_y);
      iy += 4 * inc_y;
    }
    for (; j < n; j++) {
      gemv_t_1col(m, a + j * lda, x, y + iy, alpha);
      iy += inc_y;
    }
    return 0;
  }
#endif

  iy = 0;
  a_ptr = a;
  for (j = 0; j < n; j++) {
    temp = 0.0;
    ix = 0;
    for (i = 0; i < m; i++) {
      temp += a_ptr[i] * x[ix];
      ix += inc_x;
    }
    y[iy] += alpha * temp;
    iy += inc_y;
    a_ptr += lda;
  }
  return 0;
}
