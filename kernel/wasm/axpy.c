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
 * WASM SIMD128 AXPY: y += da * x
 *
 * Compiled twice: SAXPY (-UDOUBLE) and DAXPY (-DDOUBLE). Unit-stride uses
 * eight independent v128 lanes (32 floats / 16 doubles) so load/mul/add/store
 * can overlap; remainder is one vector then scalar. IEEE mul+add (not
 * relaxed madd): AXPY is checked to machine epsilon, and relaxed madd in
 * the generic V_SIMD path previously slowed L1.
 *
 * Non-unit stride stays scalar (WASM SIMD128 has no gather). inc==0 is the
 * scalar path so y[0] += n * da * x[0] still holds.
 */

#include "common.h"

#if defined(__wasm_simd128__)
#include <wasm_simd128.h>

#ifdef DOUBLE
#define AXPY_VLEN 2
#define AXPY_SPLAT wasm_f64x2_splat
#define AXPY_MUL wasm_f64x2_mul
#define AXPY_ADD wasm_f64x2_add
#else
#define AXPY_VLEN 4
#define AXPY_SPLAT wasm_f32x4_splat
#define AXPY_MUL wasm_f32x4_mul
#define AXPY_ADD wasm_f32x4_add
#endif

#define AXPY_UNROLL 8
#define AXPY_CHUNK (AXPY_VLEN * AXPY_UNROLL)

#define AXPY_LOAD(p) wasm_v128_load((const void *)(p))
#define AXPY_STORE(p, v) wasm_v128_store((void *)(p), (v))
#define AXPY_MADD(y, a, x) AXPY_ADD((y), AXPY_MUL((a), (x)))

static void axpy_kernel_unit(BLASLONG n, const FLOAT *x, FLOAT *y, FLOAT da) {
  const v128_t va = AXPY_SPLAT(da);
  BLASLONG i = 0;
  const BLASLONG n_main = n & ~(BLASLONG)(AXPY_CHUNK - 1);

  for (; i < n_main; i += AXPY_CHUNK) {
    v128_t x0 = AXPY_LOAD(x + i + 0 * AXPY_VLEN);
    v128_t x1 = AXPY_LOAD(x + i + 1 * AXPY_VLEN);
    v128_t x2 = AXPY_LOAD(x + i + 2 * AXPY_VLEN);
    v128_t x3 = AXPY_LOAD(x + i + 3 * AXPY_VLEN);
    v128_t x4 = AXPY_LOAD(x + i + 4 * AXPY_VLEN);
    v128_t x5 = AXPY_LOAD(x + i + 5 * AXPY_VLEN);
    v128_t x6 = AXPY_LOAD(x + i + 6 * AXPY_VLEN);
    v128_t x7 = AXPY_LOAD(x + i + 7 * AXPY_VLEN);

    v128_t y0 = AXPY_LOAD(y + i + 0 * AXPY_VLEN);
    v128_t y1 = AXPY_LOAD(y + i + 1 * AXPY_VLEN);
    v128_t y2 = AXPY_LOAD(y + i + 2 * AXPY_VLEN);
    v128_t y3 = AXPY_LOAD(y + i + 3 * AXPY_VLEN);
    v128_t y4 = AXPY_LOAD(y + i + 4 * AXPY_VLEN);
    v128_t y5 = AXPY_LOAD(y + i + 5 * AXPY_VLEN);
    v128_t y6 = AXPY_LOAD(y + i + 6 * AXPY_VLEN);
    v128_t y7 = AXPY_LOAD(y + i + 7 * AXPY_VLEN);

    AXPY_STORE(y + i + 0 * AXPY_VLEN, AXPY_MADD(y0, va, x0));
    AXPY_STORE(y + i + 1 * AXPY_VLEN, AXPY_MADD(y1, va, x1));
    AXPY_STORE(y + i + 2 * AXPY_VLEN, AXPY_MADD(y2, va, x2));
    AXPY_STORE(y + i + 3 * AXPY_VLEN, AXPY_MADD(y3, va, x3));
    AXPY_STORE(y + i + 4 * AXPY_VLEN, AXPY_MADD(y4, va, x4));
    AXPY_STORE(y + i + 5 * AXPY_VLEN, AXPY_MADD(y5, va, x5));
    AXPY_STORE(y + i + 6 * AXPY_VLEN, AXPY_MADD(y6, va, x6));
    AXPY_STORE(y + i + 7 * AXPY_VLEN, AXPY_MADD(y7, va, x7));
  }

  for (; i + AXPY_VLEN <= n; i += AXPY_VLEN) {
    v128_t yi = AXPY_LOAD(y + i);
    v128_t xi = AXPY_LOAD(x + i);
    AXPY_STORE(y + i, AXPY_MADD(yi, va, xi));
  }

  for (; i < n; i++)
    y[i] += da * x[i];
}

#undef AXPY_LOAD
#undef AXPY_STORE
#undef AXPY_MADD
#undef AXPY_CHUNK
#undef AXPY_UNROLL
#undef AXPY_VLEN
#undef AXPY_SPLAT
#undef AXPY_MUL
#undef AXPY_ADD
#endif

int CNAME(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT da, FLOAT *x,
          BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *dummy,
          BLASLONG dummy2) {
  BLASLONG i = 0;
  BLASLONG ix = 0, iy = 0;

  (void)dummy0;
  (void)dummy1;
  (void)dummy;
  (void)dummy2;

  if (n <= 0)
    return 0;
  if (da == 0.0)
    return 0;

  if ((inc_x == 1) && (inc_y == 1)) {
#if defined(__wasm_simd128__)
    axpy_kernel_unit(n, x, y, da);
#else
    while (i < n) {
      y[i] += da * x[i];
      i++;
    }
#endif
    return 0;
  }

  {
    BLASLONG n1 = n & ~(BLASLONG)3;
    while (i < n1) {
      FLOAT m1 = da * x[ix];
      FLOAT m2 = da * x[ix + inc_x];
      FLOAT m3 = da * x[ix + 2 * inc_x];
      FLOAT m4 = da * x[ix + 3 * inc_x];

      y[iy] += m1;
      y[iy + inc_y] += m2;
      y[iy + 2 * inc_y] += m3;
      y[iy + 3 * inc_y] += m4;

      ix += inc_x * 4;
      iy += inc_y * 4;
      i += 4;
    }
  }

  while (i < n) {
    y[iy] += da * x[ix];
    ix += inc_x;
    iy += inc_y;
    i++;
  }
  return 0;
}
