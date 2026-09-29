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
 * WASM SIMD128 GEMM micro-kernel:
 *   single: 8x4 (adjacent f32x4 rows for V8 WasmRevecReducer)
 *   double: 4x8 (adjacent f64x2 rows; Haswell-like N unroll)
 *
 * Packed-data contract matches gemm_{n,t}copy_{8,4}:
 *   SGEMM A is 8-wide, B is 4-wide; DGEMM A is 4-wide, B is 8-wide.
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
#ifdef DOUBLE
  FLOAT *C4, *C5, *C6, *C7;
#endif
  IFLOAT *ptrba, *ptrbb;
  FLOAT r0c0, r1c0, r2c0, r3c0;
  FLOAT r0c1, r1c1, r2c1, r3c1;
  FLOAT r0c2, r1c2, r2c2, r3c2;
  FLOAT r0c3, r1c3, r2c3, r3c3;
  IFLOAT a0, a1, a2, a3;
#ifndef DOUBLE
  IFLOAT a4, a5, a6, a7;
#endif
  IFLOAT b0, b1, b2, b3;
#ifdef DOUBLE
  IFLOAT b4, b5, b6, b7;
#endif

#ifndef DOUBLE
  /* ---- SGEMM 8x4 ---- */
  for (j = 0; j < bn / 4; j += 1) {
    C0 = C;
    C1 = C0 + ldc;
    C2 = C1 + ldc;
    C3 = C2 + ldc;
    ptrba = ba;

    for (i = 0; i < bm / 8; i += 1) {
      ptrbb = bb;
#if defined(__wasm_simd128__) && !defined(BGEMM)
      {
        v128_t acc0_lo = wasm_f32x4_splat(0.0f);
        v128_t acc0_hi = wasm_f32x4_splat(0.0f);
        v128_t acc1_lo = wasm_f32x4_splat(0.0f);
        v128_t acc1_hi = wasm_f32x4_splat(0.0f);
        v128_t acc2_lo = wasm_f32x4_splat(0.0f);
        v128_t acc2_hi = wasm_f32x4_splat(0.0f);
        v128_t acc3_lo = wasm_f32x4_splat(0.0f);
        v128_t acc3_hi = wasm_f32x4_splat(0.0f);
        for (k = 0; k < bk; k += 1) {
          v128_t va_lo = wasm_v128_load(ptrba);
          v128_t va_hi = wasm_v128_load(ptrba + 4);
          v128_t vb0 = wasm_v128_load32_splat(ptrbb + 0);
          v128_t vb1 = wasm_v128_load32_splat(ptrbb + 1);
          v128_t vb2 = wasm_v128_load32_splat(ptrbb + 2);
          v128_t vb3 = wasm_v128_load32_splat(ptrbb + 3);
          acc0_lo = MADD_F32(va_lo, vb0, acc0_lo);
          acc0_hi = MADD_F32(va_hi, vb0, acc0_hi);
          acc1_lo = MADD_F32(va_lo, vb1, acc1_lo);
          acc1_hi = MADD_F32(va_hi, vb1, acc1_hi);
          acc2_lo = MADD_F32(va_lo, vb2, acc2_lo);
          acc2_hi = MADD_F32(va_hi, vb2, acc2_hi);
          acc3_lo = MADD_F32(va_lo, vb3, acc3_lo);
          acc3_hi = MADD_F32(va_hi, vb3, acc3_hi);
          ptrba += 8;
          ptrbb += 4;
        }
        v128_t valpha = wasm_f32x4_splat(alpha);
        wasm_v128_store(C0, MADD_F32(acc0_lo, valpha, wasm_v128_load(C0)));
        wasm_v128_store(C0 + 4, MADD_F32(acc0_hi, valpha, wasm_v128_load(C0 + 4)));
        wasm_v128_store(C1, MADD_F32(acc1_lo, valpha, wasm_v128_load(C1)));
        wasm_v128_store(C1 + 4, MADD_F32(acc1_hi, valpha, wasm_v128_load(C1 + 4)));
        wasm_v128_store(C2, MADD_F32(acc2_lo, valpha, wasm_v128_load(C2)));
        wasm_v128_store(C2 + 4, MADD_F32(acc2_hi, valpha, wasm_v128_load(C2 + 4)));
        wasm_v128_store(C3, MADD_F32(acc3_lo, valpha, wasm_v128_load(C3)));
        wasm_v128_store(C3 + 4, MADD_F32(acc3_hi, valpha, wasm_v128_load(C3 + 4)));
      }
#else
      {
        FLOAT r4c0, r5c0, r6c0, r7c0;
        FLOAT r4c1, r5c1, r6c1, r7c1;
        FLOAT r4c2, r5c2, r6c2, r7c2;
        FLOAT r4c3, r5c3, r6c3, r7c3;
        r0c0 = r1c0 = r2c0 = r3c0 = r4c0 = r5c0 = r6c0 = r7c0 = 0;
        r0c1 = r1c1 = r2c1 = r3c1 = r4c1 = r5c1 = r6c1 = r7c1 = 0;
        r0c2 = r1c2 = r2c2 = r3c2 = r4c2 = r5c2 = r6c2 = r7c2 = 0;
        r0c3 = r1c3 = r2c3 = r3c3 = r4c3 = r5c3 = r6c3 = r7c3 = 0;
        for (k = 0; k < bk; k += 1) {
          b0 = ptrbb[0];
          b1 = ptrbb[1];
          b2 = ptrbb[2];
          b3 = ptrbb[3];
          a0 = ptrba[0];
          a1 = ptrba[1];
          a2 = ptrba[2];
          a3 = ptrba[3];
          a4 = ptrba[4];
          a5 = ptrba[5];
          a6 = ptrba[6];
          a7 = ptrba[7];
          r0c0 += TO_F32(a0) * TO_F32(b0);
          r1c0 += TO_F32(a1) * TO_F32(b0);
          r2c0 += TO_F32(a2) * TO_F32(b0);
          r3c0 += TO_F32(a3) * TO_F32(b0);
          r4c0 += TO_F32(a4) * TO_F32(b0);
          r5c0 += TO_F32(a5) * TO_F32(b0);
          r6c0 += TO_F32(a6) * TO_F32(b0);
          r7c0 += TO_F32(a7) * TO_F32(b0);
          r0c1 += TO_F32(a0) * TO_F32(b1);
          r1c1 += TO_F32(a1) * TO_F32(b1);
          r2c1 += TO_F32(a2) * TO_F32(b1);
          r3c1 += TO_F32(a3) * TO_F32(b1);
          r4c1 += TO_F32(a4) * TO_F32(b1);
          r5c1 += TO_F32(a5) * TO_F32(b1);
          r6c1 += TO_F32(a6) * TO_F32(b1);
          r7c1 += TO_F32(a7) * TO_F32(b1);
          r0c2 += TO_F32(a0) * TO_F32(b2);
          r1c2 += TO_F32(a1) * TO_F32(b2);
          r2c2 += TO_F32(a2) * TO_F32(b2);
          r3c2 += TO_F32(a3) * TO_F32(b2);
          r4c2 += TO_F32(a4) * TO_F32(b2);
          r5c2 += TO_F32(a5) * TO_F32(b2);
          r6c2 += TO_F32(a6) * TO_F32(b2);
          r7c2 += TO_F32(a7) * TO_F32(b2);
          r0c3 += TO_F32(a0) * TO_F32(b3);
          r1c3 += TO_F32(a1) * TO_F32(b3);
          r2c3 += TO_F32(a2) * TO_F32(b3);
          r3c3 += TO_F32(a3) * TO_F32(b3);
          r4c3 += TO_F32(a4) * TO_F32(b3);
          r5c3 += TO_F32(a5) * TO_F32(b3);
          r6c3 += TO_F32(a6) * TO_F32(b3);
          r7c3 += TO_F32(a7) * TO_F32(b3);
          ptrba += 8;
          ptrbb += 4;
        }
        C0[0] = TO_OUTPUT(C_TO_F32(C0[0]) + r0c0 * ALPHA);
        C0[1] = TO_OUTPUT(C_TO_F32(C0[1]) + r1c0 * ALPHA);
        C0[2] = TO_OUTPUT(C_TO_F32(C0[2]) + r2c0 * ALPHA);
        C0[3] = TO_OUTPUT(C_TO_F32(C0[3]) + r3c0 * ALPHA);
        C0[4] = TO_OUTPUT(C_TO_F32(C0[4]) + r4c0 * ALPHA);
        C0[5] = TO_OUTPUT(C_TO_F32(C0[5]) + r5c0 * ALPHA);
        C0[6] = TO_OUTPUT(C_TO_F32(C0[6]) + r6c0 * ALPHA);
        C0[7] = TO_OUTPUT(C_TO_F32(C0[7]) + r7c0 * ALPHA);
        C1[0] = TO_OUTPUT(C_TO_F32(C1[0]) + r0c1 * ALPHA);
        C1[1] = TO_OUTPUT(C_TO_F32(C1[1]) + r1c1 * ALPHA);
        C1[2] = TO_OUTPUT(C_TO_F32(C1[2]) + r2c1 * ALPHA);
        C1[3] = TO_OUTPUT(C_TO_F32(C1[3]) + r3c1 * ALPHA);
        C1[4] = TO_OUTPUT(C_TO_F32(C1[4]) + r4c1 * ALPHA);
        C1[5] = TO_OUTPUT(C_TO_F32(C1[5]) + r5c1 * ALPHA);
        C1[6] = TO_OUTPUT(C_TO_F32(C1[6]) + r6c1 * ALPHA);
        C1[7] = TO_OUTPUT(C_TO_F32(C1[7]) + r7c1 * ALPHA);
        C2[0] = TO_OUTPUT(C_TO_F32(C2[0]) + r0c2 * ALPHA);
        C2[1] = TO_OUTPUT(C_TO_F32(C2[1]) + r1c2 * ALPHA);
        C2[2] = TO_OUTPUT(C_TO_F32(C2[2]) + r2c2 * ALPHA);
        C2[3] = TO_OUTPUT(C_TO_F32(C2[3]) + r3c2 * ALPHA);
        C2[4] = TO_OUTPUT(C_TO_F32(C2[4]) + r4c2 * ALPHA);
        C2[5] = TO_OUTPUT(C_TO_F32(C2[5]) + r5c2 * ALPHA);
        C2[6] = TO_OUTPUT(C_TO_F32(C2[6]) + r6c2 * ALPHA);
        C2[7] = TO_OUTPUT(C_TO_F32(C2[7]) + r7c2 * ALPHA);
        C3[0] = TO_OUTPUT(C_TO_F32(C3[0]) + r0c3 * ALPHA);
        C3[1] = TO_OUTPUT(C_TO_F32(C3[1]) + r1c3 * ALPHA);
        C3[2] = TO_OUTPUT(C_TO_F32(C3[2]) + r2c3 * ALPHA);
        C3[3] = TO_OUTPUT(C_TO_F32(C3[3]) + r3c3 * ALPHA);
        C3[4] = TO_OUTPUT(C_TO_F32(C3[4]) + r4c3 * ALPHA);
        C3[5] = TO_OUTPUT(C_TO_F32(C3[5]) + r5c3 * ALPHA);
        C3[6] = TO_OUTPUT(C_TO_F32(C3[6]) + r6c3 * ALPHA);
        C3[7] = TO_OUTPUT(C_TO_F32(C3[7]) + r7c3 * ALPHA);
      }
#endif
      C0 += 8;
      C1 += 8;
      C2 += 8;
      C3 += 8;
    }

    if (bm & 4) {
      ptrbb = bb;
#if defined(__wasm_simd128__) && !defined(BGEMM)
      {
        v128_t acc0 = wasm_f32x4_splat(0.0f);
        v128_t acc1 = wasm_f32x4_splat(0.0f);
        v128_t acc2 = wasm_f32x4_splat(0.0f);
        v128_t acc3 = wasm_f32x4_splat(0.0f);
        for (k = 0; k < bk; k += 1) {
          v128_t va = wasm_v128_load(ptrba);
          v128_t vb0 = wasm_v128_load32_splat(ptrbb + 0);
          v128_t vb1 = wasm_v128_load32_splat(ptrbb + 1);
          v128_t vb2 = wasm_v128_load32_splat(ptrbb + 2);
          v128_t vb3 = wasm_v128_load32_splat(ptrbb + 3);
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
    for (i = 0; i < bm / 8; i += 1) {
      ptrbb = bb;
      {
        FLOAT r4c0, r5c0, r6c0, r7c0;
        FLOAT r4c1, r5c1, r6c1, r7c1;
        r0c0 = r1c0 = r2c0 = r3c0 = r4c0 = r5c0 = r6c0 = r7c0 = 0;
        r0c1 = r1c1 = r2c1 = r3c1 = r4c1 = r5c1 = r6c1 = r7c1 = 0;
        for (k = 0; k < bk; k += 1) {
          b0 = ptrbb[0];
          b1 = ptrbb[1];
          a0 = ptrba[0];
          a1 = ptrba[1];
          a2 = ptrba[2];
          a3 = ptrba[3];
          a4 = ptrba[4];
          a5 = ptrba[5];
          a6 = ptrba[6];
          a7 = ptrba[7];
          r0c0 += TO_F32(a0) * TO_F32(b0);
          r1c0 += TO_F32(a1) * TO_F32(b0);
          r2c0 += TO_F32(a2) * TO_F32(b0);
          r3c0 += TO_F32(a3) * TO_F32(b0);
          r4c0 += TO_F32(a4) * TO_F32(b0);
          r5c0 += TO_F32(a5) * TO_F32(b0);
          r6c0 += TO_F32(a6) * TO_F32(b0);
          r7c0 += TO_F32(a7) * TO_F32(b0);
          r0c1 += TO_F32(a0) * TO_F32(b1);
          r1c1 += TO_F32(a1) * TO_F32(b1);
          r2c1 += TO_F32(a2) * TO_F32(b1);
          r3c1 += TO_F32(a3) * TO_F32(b1);
          r4c1 += TO_F32(a4) * TO_F32(b1);
          r5c1 += TO_F32(a5) * TO_F32(b1);
          r6c1 += TO_F32(a6) * TO_F32(b1);
          r7c1 += TO_F32(a7) * TO_F32(b1);
          ptrba += 8;
          ptrbb += 2;
        }
        C0[0] = TO_OUTPUT(C_TO_F32(C0[0]) + r0c0 * ALPHA);
        C0[1] = TO_OUTPUT(C_TO_F32(C0[1]) + r1c0 * ALPHA);
        C0[2] = TO_OUTPUT(C_TO_F32(C0[2]) + r2c0 * ALPHA);
        C0[3] = TO_OUTPUT(C_TO_F32(C0[3]) + r3c0 * ALPHA);
        C0[4] = TO_OUTPUT(C_TO_F32(C0[4]) + r4c0 * ALPHA);
        C0[5] = TO_OUTPUT(C_TO_F32(C0[5]) + r5c0 * ALPHA);
        C0[6] = TO_OUTPUT(C_TO_F32(C0[6]) + r6c0 * ALPHA);
        C0[7] = TO_OUTPUT(C_TO_F32(C0[7]) + r7c0 * ALPHA);
        C1[0] = TO_OUTPUT(C_TO_F32(C1[0]) + r0c1 * ALPHA);
        C1[1] = TO_OUTPUT(C_TO_F32(C1[1]) + r1c1 * ALPHA);
        C1[2] = TO_OUTPUT(C_TO_F32(C1[2]) + r2c1 * ALPHA);
        C1[3] = TO_OUTPUT(C_TO_F32(C1[3]) + r3c1 * ALPHA);
        C1[4] = TO_OUTPUT(C_TO_F32(C1[4]) + r4c1 * ALPHA);
        C1[5] = TO_OUTPUT(C_TO_F32(C1[5]) + r5c1 * ALPHA);
        C1[6] = TO_OUTPUT(C_TO_F32(C1[6]) + r6c1 * ALPHA);
        C1[7] = TO_OUTPUT(C_TO_F32(C1[7]) + r7c1 * ALPHA);
      }
      C0 += 8;
      C1 += 8;
    }
    if (bm & 4) {
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
    for (i = 0; i < bm / 8; i += 1) {
      ptrbb = bb;
      {
        FLOAT r4c0, r5c0, r6c0, r7c0;
        r0c0 = r1c0 = r2c0 = r3c0 = r4c0 = r5c0 = r6c0 = r7c0 = 0;
        for (k = 0; k < bk; k += 1) {
          b0 = ptrbb[0];
          a0 = ptrba[0];
          a1 = ptrba[1];
          a2 = ptrba[2];
          a3 = ptrba[3];
          a4 = ptrba[4];
          a5 = ptrba[5];
          a6 = ptrba[6];
          a7 = ptrba[7];
          r0c0 += TO_F32(a0) * TO_F32(b0);
          r1c0 += TO_F32(a1) * TO_F32(b0);
          r2c0 += TO_F32(a2) * TO_F32(b0);
          r3c0 += TO_F32(a3) * TO_F32(b0);
          r4c0 += TO_F32(a4) * TO_F32(b0);
          r5c0 += TO_F32(a5) * TO_F32(b0);
          r6c0 += TO_F32(a6) * TO_F32(b0);
          r7c0 += TO_F32(a7) * TO_F32(b0);
          ptrba += 8;
          ptrbb += 1;
        }
        C0[0] = TO_OUTPUT(C_TO_F32(C0[0]) + r0c0 * ALPHA);
        C0[1] = TO_OUTPUT(C_TO_F32(C0[1]) + r1c0 * ALPHA);
        C0[2] = TO_OUTPUT(C_TO_F32(C0[2]) + r2c0 * ALPHA);
        C0[3] = TO_OUTPUT(C_TO_F32(C0[3]) + r3c0 * ALPHA);
        C0[4] = TO_OUTPUT(C_TO_F32(C0[4]) + r4c0 * ALPHA);
        C0[5] = TO_OUTPUT(C_TO_F32(C0[5]) + r5c0 * ALPHA);
        C0[6] = TO_OUTPUT(C_TO_F32(C0[6]) + r6c0 * ALPHA);
        C0[7] = TO_OUTPUT(C_TO_F32(C0[7]) + r7c0 * ALPHA);
      }
      C0 += 8;
    }
    if (bm & 4) {
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

#else
  /* ---- DGEMM 4x8 ---- */
  for (j = 0; j < bn / 8; j += 1) {
    C0 = C;
    C1 = C0 + ldc;
    C2 = C1 + ldc;
    C3 = C2 + ldc;
    C4 = C3 + ldc;
    C5 = C4 + ldc;
    C6 = C5 + ldc;
    C7 = C6 + ldc;
    ptrba = ba;

    for (i = 0; i < bm / 4; i += 1) {
      ptrbb = bb;
#if defined(__wasm_simd128__) && !defined(BGEMM)
      {
        v128_t a0l = wasm_f64x2_splat(0.0), a0h = wasm_f64x2_splat(0.0);
        v128_t a1l = wasm_f64x2_splat(0.0), a1h = wasm_f64x2_splat(0.0);
        v128_t a2l = wasm_f64x2_splat(0.0), a2h = wasm_f64x2_splat(0.0);
        v128_t a3l = wasm_f64x2_splat(0.0), a3h = wasm_f64x2_splat(0.0);
        v128_t a4l = wasm_f64x2_splat(0.0), a4h = wasm_f64x2_splat(0.0);
        v128_t a5l = wasm_f64x2_splat(0.0), a5h = wasm_f64x2_splat(0.0);
        v128_t a6l = wasm_f64x2_splat(0.0), a6h = wasm_f64x2_splat(0.0);
        v128_t a7l = wasm_f64x2_splat(0.0), a7h = wasm_f64x2_splat(0.0);
        for (k = 0; k < bk; k += 1) {
          v128_t va01 = wasm_v128_load(ptrba);
          v128_t va23 = wasm_v128_load(ptrba + 2);
          v128_t b0v = wasm_v128_load64_splat(ptrbb + 0);
          v128_t b1v = wasm_v128_load64_splat(ptrbb + 1);
          v128_t b2v = wasm_v128_load64_splat(ptrbb + 2);
          v128_t b3v = wasm_v128_load64_splat(ptrbb + 3);
          v128_t b4v = wasm_v128_load64_splat(ptrbb + 4);
          v128_t b5v = wasm_v128_load64_splat(ptrbb + 5);
          v128_t b6v = wasm_v128_load64_splat(ptrbb + 6);
          v128_t b7v = wasm_v128_load64_splat(ptrbb + 7);
          a0l = MADD_F64(va01, b0v, a0l);
          a0h = MADD_F64(va23, b0v, a0h);
          a1l = MADD_F64(va01, b1v, a1l);
          a1h = MADD_F64(va23, b1v, a1h);
          a2l = MADD_F64(va01, b2v, a2l);
          a2h = MADD_F64(va23, b2v, a2h);
          a3l = MADD_F64(va01, b3v, a3l);
          a3h = MADD_F64(va23, b3v, a3h);
          a4l = MADD_F64(va01, b4v, a4l);
          a4h = MADD_F64(va23, b4v, a4h);
          a5l = MADD_F64(va01, b5v, a5l);
          a5h = MADD_F64(va23, b5v, a5h);
          a6l = MADD_F64(va01, b6v, a6l);
          a6h = MADD_F64(va23, b6v, a6h);
          a7l = MADD_F64(va01, b7v, a7l);
          a7h = MADD_F64(va23, b7v, a7h);
          ptrba += 4;
          ptrbb += 8;
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
        wasm_v128_store(C4, MADD_F64(a4l, valpha, wasm_v128_load(C4)));
        wasm_v128_store(C4 + 2, MADD_F64(a4h, valpha, wasm_v128_load(C4 + 2)));
        wasm_v128_store(C5, MADD_F64(a5l, valpha, wasm_v128_load(C5)));
        wasm_v128_store(C5 + 2, MADD_F64(a5h, valpha, wasm_v128_load(C5 + 2)));
        wasm_v128_store(C6, MADD_F64(a6l, valpha, wasm_v128_load(C6)));
        wasm_v128_store(C6 + 2, MADD_F64(a6h, valpha, wasm_v128_load(C6 + 2)));
        wasm_v128_store(C7, MADD_F64(a7l, valpha, wasm_v128_load(C7)));
        wasm_v128_store(C7 + 2, MADD_F64(a7h, valpha, wasm_v128_load(C7 + 2)));
      }
#else
      {
        FLOAT r0c4, r1c4, r2c4, r3c4;
        FLOAT r0c5, r1c5, r2c5, r3c5;
        FLOAT r0c6, r1c6, r2c6, r3c6;
        FLOAT r0c7, r1c7, r2c7, r3c7;
        r0c0 = r1c0 = r2c0 = r3c0 = 0;
        r0c1 = r1c1 = r2c1 = r3c1 = 0;
        r0c2 = r1c2 = r2c2 = r3c2 = 0;
        r0c3 = r1c3 = r2c3 = r3c3 = 0;
        r0c4 = r1c4 = r2c4 = r3c4 = 0;
        r0c5 = r1c5 = r2c5 = r3c5 = 0;
        r0c6 = r1c6 = r2c6 = r3c6 = 0;
        r0c7 = r1c7 = r2c7 = r3c7 = 0;
        for (k = 0; k < bk; k += 1) {
          b0 = ptrbb[0];
          b1 = ptrbb[1];
          b2 = ptrbb[2];
          b3 = ptrbb[3];
          b4 = ptrbb[4];
          b5 = ptrbb[5];
          b6 = ptrbb[6];
          b7 = ptrbb[7];
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
          r0c4 += TO_F32(a0) * TO_F32(b4);
          r1c4 += TO_F32(a1) * TO_F32(b4);
          r2c4 += TO_F32(a2) * TO_F32(b4);
          r3c4 += TO_F32(a3) * TO_F32(b4);
          r0c5 += TO_F32(a0) * TO_F32(b5);
          r1c5 += TO_F32(a1) * TO_F32(b5);
          r2c5 += TO_F32(a2) * TO_F32(b5);
          r3c5 += TO_F32(a3) * TO_F32(b5);
          r0c6 += TO_F32(a0) * TO_F32(b6);
          r1c6 += TO_F32(a1) * TO_F32(b6);
          r2c6 += TO_F32(a2) * TO_F32(b6);
          r3c6 += TO_F32(a3) * TO_F32(b6);
          r0c7 += TO_F32(a0) * TO_F32(b7);
          r1c7 += TO_F32(a1) * TO_F32(b7);
          r2c7 += TO_F32(a2) * TO_F32(b7);
          r3c7 += TO_F32(a3) * TO_F32(b7);
          ptrba += 4;
          ptrbb += 8;
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
        C4[0] = TO_OUTPUT(C_TO_F32(C4[0]) + r0c4 * ALPHA);
        C4[1] = TO_OUTPUT(C_TO_F32(C4[1]) + r1c4 * ALPHA);
        C4[2] = TO_OUTPUT(C_TO_F32(C4[2]) + r2c4 * ALPHA);
        C4[3] = TO_OUTPUT(C_TO_F32(C4[3]) + r3c4 * ALPHA);
        C5[0] = TO_OUTPUT(C_TO_F32(C5[0]) + r0c5 * ALPHA);
        C5[1] = TO_OUTPUT(C_TO_F32(C5[1]) + r1c5 * ALPHA);
        C5[2] = TO_OUTPUT(C_TO_F32(C5[2]) + r2c5 * ALPHA);
        C5[3] = TO_OUTPUT(C_TO_F32(C5[3]) + r3c5 * ALPHA);
        C6[0] = TO_OUTPUT(C_TO_F32(C6[0]) + r0c6 * ALPHA);
        C6[1] = TO_OUTPUT(C_TO_F32(C6[1]) + r1c6 * ALPHA);
        C6[2] = TO_OUTPUT(C_TO_F32(C6[2]) + r2c6 * ALPHA);
        C6[3] = TO_OUTPUT(C_TO_F32(C6[3]) + r3c6 * ALPHA);
        C7[0] = TO_OUTPUT(C_TO_F32(C7[0]) + r0c7 * ALPHA);
        C7[1] = TO_OUTPUT(C_TO_F32(C7[1]) + r1c7 * ALPHA);
        C7[2] = TO_OUTPUT(C_TO_F32(C7[2]) + r2c7 * ALPHA);
        C7[3] = TO_OUTPUT(C_TO_F32(C7[3]) + r3c7 * ALPHA);
      }
#endif
      C0 += 4;
      C1 += 4;
      C2 += 4;
      C3 += 4;
      C4 += 4;
      C5 += 4;
      C6 += 4;
      C7 += 4;
    }

    if (bm & 2) {
      ptrbb = bb;
      {
        FLOAT r0c4, r1c4, r0c5, r1c5, r0c6, r1c6, r0c7, r1c7;
        r0c0 = r1c0 = r0c1 = r1c1 = r0c2 = r1c2 = r0c3 = r1c3 = 0;
        r0c4 = r1c4 = r0c5 = r1c5 = r0c6 = r1c6 = r0c7 = r1c7 = 0;
        for (k = 0; k < bk; k += 1) {
          b0 = ptrbb[0];
          b1 = ptrbb[1];
          b2 = ptrbb[2];
          b3 = ptrbb[3];
          b4 = ptrbb[4];
          b5 = ptrbb[5];
          b6 = ptrbb[6];
          b7 = ptrbb[7];
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
          r0c4 += TO_F32(a0) * TO_F32(b4);
          r1c4 += TO_F32(a1) * TO_F32(b4);
          r0c5 += TO_F32(a0) * TO_F32(b5);
          r1c5 += TO_F32(a1) * TO_F32(b5);
          r0c6 += TO_F32(a0) * TO_F32(b6);
          r1c6 += TO_F32(a1) * TO_F32(b6);
          r0c7 += TO_F32(a0) * TO_F32(b7);
          r1c7 += TO_F32(a1) * TO_F32(b7);
          ptrba += 2;
          ptrbb += 8;
        }
        C0[0] = TO_OUTPUT(C_TO_F32(C0[0]) + r0c0 * ALPHA);
        C0[1] = TO_OUTPUT(C_TO_F32(C0[1]) + r1c0 * ALPHA);
        C1[0] = TO_OUTPUT(C_TO_F32(C1[0]) + r0c1 * ALPHA);
        C1[1] = TO_OUTPUT(C_TO_F32(C1[1]) + r1c1 * ALPHA);
        C2[0] = TO_OUTPUT(C_TO_F32(C2[0]) + r0c2 * ALPHA);
        C2[1] = TO_OUTPUT(C_TO_F32(C2[1]) + r1c2 * ALPHA);
        C3[0] = TO_OUTPUT(C_TO_F32(C3[0]) + r0c3 * ALPHA);
        C3[1] = TO_OUTPUT(C_TO_F32(C3[1]) + r1c3 * ALPHA);
        C4[0] = TO_OUTPUT(C_TO_F32(C4[0]) + r0c4 * ALPHA);
        C4[1] = TO_OUTPUT(C_TO_F32(C4[1]) + r1c4 * ALPHA);
        C5[0] = TO_OUTPUT(C_TO_F32(C5[0]) + r0c5 * ALPHA);
        C5[1] = TO_OUTPUT(C_TO_F32(C5[1]) + r1c5 * ALPHA);
        C6[0] = TO_OUTPUT(C_TO_F32(C6[0]) + r0c6 * ALPHA);
        C6[1] = TO_OUTPUT(C_TO_F32(C6[1]) + r1c6 * ALPHA);
        C7[0] = TO_OUTPUT(C_TO_F32(C7[0]) + r0c7 * ALPHA);
        C7[1] = TO_OUTPUT(C_TO_F32(C7[1]) + r1c7 * ALPHA);
      }
      C0 += 2;
      C1 += 2;
      C2 += 2;
      C3 += 2;
      C4 += 2;
      C5 += 2;
      C6 += 2;
      C7 += 2;
    }
    if (bm & 1) {
      ptrbb = bb;
      {
        FLOAT r0c4, r0c5, r0c6, r0c7;
        r0c0 = r0c1 = r0c2 = r0c3 = r0c4 = r0c5 = r0c6 = r0c7 = 0;
        for (k = 0; k < bk; k += 1) {
          a0 = ptrba[0];
          r0c0 += TO_F32(a0) * TO_F32(ptrbb[0]);
          r0c1 += TO_F32(a0) * TO_F32(ptrbb[1]);
          r0c2 += TO_F32(a0) * TO_F32(ptrbb[2]);
          r0c3 += TO_F32(a0) * TO_F32(ptrbb[3]);
          r0c4 += TO_F32(a0) * TO_F32(ptrbb[4]);
          r0c5 += TO_F32(a0) * TO_F32(ptrbb[5]);
          r0c6 += TO_F32(a0) * TO_F32(ptrbb[6]);
          r0c7 += TO_F32(a0) * TO_F32(ptrbb[7]);
          ptrba += 1;
          ptrbb += 8;
        }
        C0[0] = TO_OUTPUT(C_TO_F32(C0[0]) + r0c0 * ALPHA);
        C1[0] = TO_OUTPUT(C_TO_F32(C1[0]) + r0c1 * ALPHA);
        C2[0] = TO_OUTPUT(C_TO_F32(C2[0]) + r0c2 * ALPHA);
        C3[0] = TO_OUTPUT(C_TO_F32(C3[0]) + r0c3 * ALPHA);
        C4[0] = TO_OUTPUT(C_TO_F32(C4[0]) + r0c4 * ALPHA);
        C5[0] = TO_OUTPUT(C_TO_F32(C5[0]) + r0c5 * ALPHA);
        C6[0] = TO_OUTPUT(C_TO_F32(C6[0]) + r0c6 * ALPHA);
        C7[0] = TO_OUTPUT(C_TO_F32(C7[0]) + r0c7 * ALPHA);
      }
      C0 += 1;
      C1 += 1;
      C2 += 1;
      C3 += 1;
      C4 += 1;
      C5 += 1;
      C6 += 1;
      C7 += 1;
    }
    bb = bb + bk * 8;
    C = C + ldc * 8;
  }

  if (bn & 4) {
    C0 = C;
    C1 = C0 + ldc;
    C2 = C1 + ldc;
    C3 = C2 + ldc;
    ptrba = ba;
    for (i = 0; i < bm / 4; i += 1) {
      ptrbb = bb;
#if defined(__wasm_simd128__) && !defined(BGEMM)
      {
        v128_t a0l = wasm_f64x2_splat(0.0), a0h = wasm_f64x2_splat(0.0);
        v128_t a1l = wasm_f64x2_splat(0.0), a1h = wasm_f64x2_splat(0.0);
        v128_t a2l = wasm_f64x2_splat(0.0), a2h = wasm_f64x2_splat(0.0);
        v128_t a3l = wasm_f64x2_splat(0.0), a3h = wasm_f64x2_splat(0.0);
        for (k = 0; k < bk; k += 1) {
          v128_t va01 = wasm_v128_load(ptrba);
          v128_t va23 = wasm_v128_load(ptrba + 2);
          v128_t b0v = wasm_v128_load64_splat(ptrbb + 0);
          v128_t b1v = wasm_v128_load64_splat(ptrbb + 1);
          v128_t b2v = wasm_v128_load64_splat(ptrbb + 2);
          v128_t b3v = wasm_v128_load64_splat(ptrbb + 3);
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
      r0c0 = r1c0 = r0c1 = r1c1 = r0c2 = r1c2 = r0c3 = r1c3 = 0;
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
#endif

  return 0;
}
