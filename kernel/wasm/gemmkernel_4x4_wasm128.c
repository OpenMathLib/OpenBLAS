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
 * WASM SIMD128 GEMM micro-kernel, 4x4 register tile.
 *
 * Packed-data contract matches generic gemmkernel_4x4.c / tcopy_4 / ncopy_4:
 * A panel is [A(r0,k)..A(r3,k)] per k; B panel is [B(k,c0)..B(k,c3)] per k.
 */

#include "common.h"
#include "../generic/conversion_macros.h"

#if defined(__wasm_simd128__)
#include <wasm_simd128.h>
#if defined(__wasm_relaxed_simd__)
#define MADD_F32(a, b, c) wasm_f32x4_relaxed_madd((a), (b), (c))
#define MADD_F64(a, b, c) wasm_f64x2_relaxed_madd((a), (b), (c))
#else
#define MADD_F32(a, b, c) wasm_f32x4_add((c), wasm_f32x4_mul((a), (b)))
#define MADD_F64(a, b, c) wasm_f64x2_add((c), wasm_f64x2_mul((a), (b)))
#endif
#endif

#ifdef BGEMM
#define C_TO_F32 TO_F32
#else
#define C_TO_F32
#endif

int CNAME(BLASLONG bm, BLASLONG bn, BLASLONG bk, FLOAT alpha, IFLOAT *ba,
          IFLOAT *bb, FLOAT *C, BLASLONG ldc
#ifdef TRMMKERNEL
          ,
          BLASLONG offset
#endif
) {
  BLASLONG i, j, k;
  FLOAT *C0, *C1, *C2, *C3;
  IFLOAT *ptrba, *ptrbb;
  FLOAT r0c0, r1c0, r2c0, r3c0;
  FLOAT r0c1, r1c1, r2c1, r3c1;
  FLOAT r0c2, r1c2, r2c2, r3c2;
  FLOAT r0c3, r1c3, r2c3, r3c3;
  IFLOAT a0, a1, a2, a3, b0, b1, b2, b3;

  for (j = 0; j < bn / 4; j += 1) {
    C0 = C;
    C1 = C0 + ldc;
    C2 = C1 + ldc;
    C3 = C2 + ldc;
    ptrba = ba;

    for (i = 0; i < bm / 4; i += 1) {
      ptrbb = bb;
#if defined(__wasm_simd128__) && !defined(BGEMM)
#ifndef DOUBLE
      {
        v128_t acc0 = wasm_f32x4_splat(0.0f);
        v128_t acc1 = wasm_f32x4_splat(0.0f);
        v128_t acc2 = wasm_f32x4_splat(0.0f);
        v128_t acc3 = wasm_f32x4_splat(0.0f);
        for (k = 0; k < bk; k += 1) {
          v128_t va = wasm_v128_load(ptrba);
          v128_t vb = wasm_v128_load(ptrbb);
          v128_t vb0 = wasm_i32x4_shuffle(vb, vb, 0, 0, 0, 0);
          v128_t vb1 = wasm_i32x4_shuffle(vb, vb, 1, 1, 1, 1);
          v128_t vb2 = wasm_i32x4_shuffle(vb, vb, 2, 2, 2, 2);
          v128_t vb3 = wasm_i32x4_shuffle(vb, vb, 3, 3, 3, 3);
          acc0 = MADD_F32(va, vb0, acc0);
          acc1 = MADD_F32(va, vb1, acc1);
          acc2 = MADD_F32(va, vb2, acc2);
          acc3 = MADD_F32(va, vb3, acc3);
          ptrba += 4;
          ptrbb += 4;
        }
        v128_t valpha = wasm_f32x4_splat(alpha);
        wasm_v128_store(C0, MADD_F32(acc0, valpha, wasm_v128_load(C0)));
        wasm_v128_store(C1, MADD_F32(acc1, valpha, wasm_v128_load(C1)));
        wasm_v128_store(C2, MADD_F32(acc2, valpha, wasm_v128_load(C2)));
        wasm_v128_store(C3, MADD_F32(acc3, valpha, wasm_v128_load(C3)));
      }
#else
      {
        v128_t a0l = wasm_f64x2_splat(0.0), a0h = wasm_f64x2_splat(0.0);
        v128_t a1l = wasm_f64x2_splat(0.0), a1h = wasm_f64x2_splat(0.0);
        v128_t a2l = wasm_f64x2_splat(0.0), a2h = wasm_f64x2_splat(0.0);
        v128_t a3l = wasm_f64x2_splat(0.0), a3h = wasm_f64x2_splat(0.0);
        for (k = 0; k < bk; k += 1) {
          v128_t va01 = wasm_v128_load(ptrba);
          v128_t va23 = wasm_v128_load(ptrba + 2);
          v128_t vb01 = wasm_v128_load(ptrbb);
          v128_t vb23 = wasm_v128_load(ptrbb + 2);
          v128_t b0v = wasm_i64x2_shuffle(vb01, vb01, 0, 0);
          v128_t b1v = wasm_i64x2_shuffle(vb01, vb01, 1, 1);
          v128_t b2v = wasm_i64x2_shuffle(vb23, vb23, 0, 0);
          v128_t b3v = wasm_i64x2_shuffle(vb23, vb23, 1, 1);
          a0l = MADD_F64(va01, b0v, a0l);
          a0h = MADD_F64(va23, b0v, a0h);
          a1l = MADD_F64(va01, b1v, a1l);
          a1h = MADD_F64(va23, b1v, a1h);
          a2l = MADD_F64(va01, b2v, a2l);
          a2h = MADD_F64(va23, b2v, a2h);
          a3l = MADD_F64(va01, b3v, a3l);
          a3h = MADD_F64(va23, b3v, a3h);
          ptrba += 4;
          ptrbb += 4;
        }
        v128_t valpha = wasm_f64x2_splat(alpha);
        wasm_v128_store(C0, MADD_F64(a0l, valpha, wasm_v128_load(C0)));
        wasm_v128_store(C0 + 2, MADD_F64(a0h, valpha, wasm_v128_load(C0 + 2)));
        wasm_v128_store(C1, MADD_F64(a1l, valpha, wasm_v128_load(C1)));
        wasm_v128_store(C1 + 2, MADD_F64(a1h, valpha, wasm_v128_load(C1 + 2)));
        wasm_v128_store(C2, MADD_F64(a2l, valpha, wasm_v128_load(C2)));
        wasm_v128_store(C2 + 2, MADD_F64(a2h, valpha, wasm_v128_load(C2 + 2)));
        wasm_v128_store(C3, MADD_F64(a3l, valpha, wasm_v128_load(C3)));
        wasm_v128_store(C3 + 2, MADD_F64(a3h, valpha, wasm_v128_load(C3 + 2)));
      }
#endif
#else
      r0c0 = r1c0 = r2c0 = r3c0 = 0;
      r0c1 = r1c1 = r2c1 = r3c1 = 0;
      r0c2 = r1c2 = r2c2 = r3c2 = 0;
      r0c3 = r1c3 = r2c3 = r3c3 = 0;
      for (k = 0; k < bk; k += 1) {
        b0 = ptrbb[0];
        b1 = ptrbb[1];
        b2 = ptrbb[2];
        b3 = ptrbb[3];
        a0 = ptrba[0];
        a1 = ptrba[1];
        a2 = ptrba[2];
        a3 = ptrba[3];
        r0c0 += TO_F32(a0) * TO_F32(b0);
        r1c0 += TO_F32(a1) * TO_F32(b0);
        r2c0 += TO_F32(a2) * TO_F32(b0);
        r3c0 += TO_F32(a3) * TO_F32(b0);
        r0c1 += TO_F32(a0) * TO_F32(b1);
        r1c1 += TO_F32(a1) * TO_F32(b1);
        r2c1 += TO_F32(a2) * TO_F32(b1);
        r3c1 += TO_F32(a3) * TO_F32(b1);
        r0c2 += TO_F32(a0) * TO_F32(b2);
        r1c2 += TO_F32(a1) * TO_F32(b2);
        r2c2 += TO_F32(a2) * TO_F32(b2);
        r3c2 += TO_F32(a3) * TO_F32(b2);
        r0c3 += TO_F32(a0) * TO_F32(b3);
        r1c3 += TO_F32(a1) * TO_F32(b3);
        r2c3 += TO_F32(a2) * TO_F32(b3);
        r3c3 += TO_F32(a3) * TO_F32(b3);
        ptrba += 4;
        ptrbb += 4;
      }
      C0[0] = TO_OUTPUT(C_TO_F32(C0[0]) + r0c0 * ALPHA);
      C0[1] = TO_OUTPUT(C_TO_F32(C0[1]) + r1c0 * ALPHA);
      C0[2] = TO_OUTPUT(C_TO_F32(C0[2]) + r2c0 * ALPHA);
      C0[3] = TO_OUTPUT(C_TO_F32(C0[3]) + r3c0 * ALPHA);
      C1[0] = TO_OUTPUT(C_TO_F32(C1[0]) + r0c1 * ALPHA);
      C1[1] = TO_OUTPUT(C_TO_F32(C1[1]) + r1c1 * ALPHA);
      C1[2] = TO_OUTPUT(C_TO_F32(C1[2]) + r2c1 * ALPHA);
      C1[3] = TO_OUTPUT(C_TO_F32(C1[3]) + r3c1 * ALPHA);
      C2[0] = TO_OUTPUT(C_TO_F32(C2[0]) + r0c2 * ALPHA);
      C2[1] = TO_OUTPUT(C_TO_F32(C2[1]) + r1c2 * ALPHA);
      C2[2] = TO_OUTPUT(C_TO_F32(C2[2]) + r2c2 * ALPHA);
      C2[3] = TO_OUTPUT(C_TO_F32(C2[3]) + r3c2 * ALPHA);
      C3[0] = TO_OUTPUT(C_TO_F32(C3[0]) + r0c3 * ALPHA);
      C3[1] = TO_OUTPUT(C_TO_F32(C3[1]) + r1c3 * ALPHA);
      C3[2] = TO_OUTPUT(C_TO_F32(C3[2]) + r2c3 * ALPHA);
      C3[3] = TO_OUTPUT(C_TO_F32(C3[3]) + r3c3 * ALPHA);
#endif
      C0 += 4;
      C1 += 4;
      C2 += 4;
      C3 += 4;
    }

    if (bm & 2) {
      ptrbb = bb;
      r0c0 = r1c0 = 0;
      r0c1 = r1c1 = 0;
      r0c2 = r1c2 = 0;
      r0c3 = r1c3 = 0;
      for (k = 0; k < bk; k += 1) {
        b0 = ptrbb[0];
        b1 = ptrbb[1];
        b2 = ptrbb[2];
        b3 = ptrbb[3];
        a0 = ptrba[0];
        a1 = ptrba[1];
        r0c0 += TO_F32(a0) * TO_F32(b0);
        r1c0 += TO_F32(a1) * TO_F32(b0);
        r0c1 += TO_F32(a0) * TO_F32(b1);
        r1c1 += TO_F32(a1) * TO_F32(b1);
        r0c2 += TO_F32(a0) * TO_F32(b2);
        r1c2 += TO_F32(a1) * TO_F32(b2);
        r0c3 += TO_F32(a0) * TO_F32(b3);
        r1c3 += TO_F32(a1) * TO_F32(b3);
        ptrba += 2;
        ptrbb += 4;
      }
      C0[0] = TO_OUTPUT(C_TO_F32(C0[0]) + r0c0 * ALPHA);
      C0[1] = TO_OUTPUT(C_TO_F32(C0[1]) + r1c0 * ALPHA);
      C1[0] = TO_OUTPUT(C_TO_F32(C1[0]) + r0c1 * ALPHA);
      C1[1] = TO_OUTPUT(C_TO_F32(C1[1]) + r1c1 * ALPHA);
      C2[0] = TO_OUTPUT(C_TO_F32(C2[0]) + r0c2 * ALPHA);
      C2[1] = TO_OUTPUT(C_TO_F32(C2[1]) + r1c2 * ALPHA);
      C3[0] = TO_OUTPUT(C_TO_F32(C3[0]) + r0c3 * ALPHA);
      C3[1] = TO_OUTPUT(C_TO_F32(C3[1]) + r1c3 * ALPHA);
      C0 += 2;
      C1 += 2;
      C2 += 2;
      C3 += 2;
    }
    if (bm & 1) {
      ptrbb = bb;
      r0c0 = r0c1 = r0c2 = r0c3 = 0;
      for (k = 0; k < bk; k += 1) {
        a0 = ptrba[0];
        r0c0 += TO_F32(a0) * TO_F32(ptrbb[0]);
        r0c1 += TO_F32(a0) * TO_F32(ptrbb[1]);
        r0c2 += TO_F32(a0) * TO_F32(ptrbb[2]);
        r0c3 += TO_F32(a0) * TO_F32(ptrbb[3]);
        ptrba += 1;
        ptrbb += 4;
      }
      C0[0] = TO_OUTPUT(C_TO_F32(C0[0]) + r0c0 * ALPHA);
      C1[0] = TO_OUTPUT(C_TO_F32(C1[0]) + r0c1 * ALPHA);
      C2[0] = TO_OUTPUT(C_TO_F32(C2[0]) + r0c2 * ALPHA);
      C3[0] = TO_OUTPUT(C_TO_F32(C3[0]) + r0c3 * ALPHA);
      C0 += 1;
      C1 += 1;
      C2 += 1;
      C3 += 1;
    }
    bb = bb + bk * 4;
    C = C + ldc * 4;
  }

  if (bn & 2) {
    C0 = C;
    C1 = C0 + ldc;
    ptrba = ba;
    for (i = 0; i < bm / 4; i += 1) {
      ptrbb = bb;
      r0c0 = r1c0 = r2c0 = r3c0 = 0;
      r0c1 = r1c1 = r2c1 = r3c1 = 0;
      for (k = 0; k < bk; k += 1) {
        b0 = ptrbb[0];
        b1 = ptrbb[1];
        a0 = ptrba[0];
        a1 = ptrba[1];
        a2 = ptrba[2];
        a3 = ptrba[3];
        r0c0 += TO_F32(a0) * TO_F32(b0);
        r1c0 += TO_F32(a1) * TO_F32(b0);
        r2c0 += TO_F32(a2) * TO_F32(b0);
        r3c0 += TO_F32(a3) * TO_F32(b0);
        r0c1 += TO_F32(a0) * TO_F32(b1);
        r1c1 += TO_F32(a1) * TO_F32(b1);
        r2c1 += TO_F32(a2) * TO_F32(b1);
        r3c1 += TO_F32(a3) * TO_F32(b1);
        ptrba += 4;
        ptrbb += 2;
      }
      C0[0] = TO_OUTPUT(C_TO_F32(C0[0]) + r0c0 * ALPHA);
      C0[1] = TO_OUTPUT(C_TO_F32(C0[1]) + r1c0 * ALPHA);
      C0[2] = TO_OUTPUT(C_TO_F32(C0[2]) + r2c0 * ALPHA);
      C0[3] = TO_OUTPUT(C_TO_F32(C0[3]) + r3c0 * ALPHA);
      C1[0] = TO_OUTPUT(C_TO_F32(C1[0]) + r0c1 * ALPHA);
      C1[1] = TO_OUTPUT(C_TO_F32(C1[1]) + r1c1 * ALPHA);
      C1[2] = TO_OUTPUT(C_TO_F32(C1[2]) + r2c1 * ALPHA);
      C1[3] = TO_OUTPUT(C_TO_F32(C1[3]) + r3c1 * ALPHA);
      C0 += 4;
      C1 += 4;
    }
    if (bm & 2) {
      ptrbb = bb;
      r0c0 = r1c0 = r0c1 = r1c1 = 0;
      for (k = 0; k < bk; k += 1) {
        b0 = ptrbb[0];
        b1 = ptrbb[1];
        a0 = ptrba[0];
        a1 = ptrba[1];
        r0c0 += TO_F32(a0) * TO_F32(b0);
        r1c0 += TO_F32(a1) * TO_F32(b0);
        r0c1 += TO_F32(a0) * TO_F32(b1);
        r1c1 += TO_F32(a1) * TO_F32(b1);
        ptrba += 2;
        ptrbb += 2;
      }
      C0[0] = TO_OUTPUT(C_TO_F32(C0[0]) + r0c0 * ALPHA);
      C0[1] = TO_OUTPUT(C_TO_F32(C0[1]) + r1c0 * ALPHA);
      C1[0] = TO_OUTPUT(C_TO_F32(C1[0]) + r0c1 * ALPHA);
      C1[1] = TO_OUTPUT(C_TO_F32(C1[1]) + r1c1 * ALPHA);
      C0 += 2;
      C1 += 2;
    }
    if (bm & 1) {
      ptrbb = bb;
      r0c0 = r0c1 = 0;
      for (k = 0; k < bk; k += 1) {
        a0 = ptrba[0];
        r0c0 += TO_F32(a0) * TO_F32(ptrbb[0]);
        r0c1 += TO_F32(a0) * TO_F32(ptrbb[1]);
        ptrba += 1;
        ptrbb += 2;
      }
      C0[0] = TO_OUTPUT(C_TO_F32(C0[0]) + r0c0 * ALPHA);
      C1[0] = TO_OUTPUT(C_TO_F32(C1[0]) + r0c1 * ALPHA);
      C0 += 1;
      C1 += 1;
    }
    bb = bb + bk * 2;
    C = C + ldc * 2;
  }

  if (bn & 1) {
    C0 = C;
    ptrba = ba;
    for (i = 0; i < bm / 4; i += 1) {
      ptrbb = bb;
      r0c0 = r1c0 = r2c0 = r3c0 = 0;
      for (k = 0; k < bk; k += 1) {
        b0 = ptrbb[0];
        a0 = ptrba[0];
        a1 = ptrba[1];
        a2 = ptrba[2];
        a3 = ptrba[3];
        r0c0 += TO_F32(a0) * TO_F32(b0);
        r1c0 += TO_F32(a1) * TO_F32(b0);
        r2c0 += TO_F32(a2) * TO_F32(b0);
        r3c0 += TO_F32(a3) * TO_F32(b0);
        ptrba += 4;
        ptrbb += 1;
      }
      C0[0] = TO_OUTPUT(C_TO_F32(C0[0]) + r0c0 * ALPHA);
      C0[1] = TO_OUTPUT(C_TO_F32(C0[1]) + r1c0 * ALPHA);
      C0[2] = TO_OUTPUT(C_TO_F32(C0[2]) + r2c0 * ALPHA);
      C0[3] = TO_OUTPUT(C_TO_F32(C0[3]) + r3c0 * ALPHA);
      C0 += 4;
    }
    if (bm & 2) {
      ptrbb = bb;
      r0c0 = r1c0 = 0;
      for (k = 0; k < bk; k += 1) {
        b0 = ptrbb[0];
        a0 = ptrba[0];
        a1 = ptrba[1];
        r0c0 += TO_F32(a0) * TO_F32(b0);
        r1c0 += TO_F32(a1) * TO_F32(b0);
        ptrba += 2;
        ptrbb += 1;
      }
      C0[0] = TO_OUTPUT(C_TO_F32(C0[0]) + r0c0 * ALPHA);
      C0[1] = TO_OUTPUT(C_TO_F32(C0[1]) + r1c0 * ALPHA);
      C0 += 2;
    }
    if (bm & 1) {
      ptrbb = bb;
      r0c0 = 0;
      for (k = 0; k < bk; k += 1) {
        r0c0 += TO_F32(ptrba[0]) * TO_F32(ptrbb[0]);
        ptrba += 1;
        ptrbb += 1;
      }
      C0[0] = TO_OUTPUT(C_TO_F32(C0[0]) + r0c0 * ALPHA);
      C0 += 1;
    }
  }

  return 0;
}
