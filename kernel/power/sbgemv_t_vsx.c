/***************************************************************************
Copyright (c) 2024, The OpenBLAS Project
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

#ifndef SBGEMV_T_VSX_C
#define SBGEMV_T_VSX_C

#include "sbgemv_common.c"

#ifndef BF16GEMV_T_X
#define BF16GEMV_T_X
#define BF16GEMV_T_8 BF16GEMV_T_VSX_8
#define BF16GEMV_T_4 BF16GEMV_T_VSX_4
#define BF16GEMV_T_2 BF16GEMV_T_VSX_2
#define BF16GEMV_T_1 BF16GEMV_T_VSX_1
#endif

#define USE_BFGEMV_8_T_VSX

static void BF16GEMV_T_VSX_1(BLASLONG n, BLASLONG lda, IFLOAT *ap, IFLOAT *x, FLOAT *y, FLOAT alpha, FLOAT beta)
{
  IFLOAT *a0;
  vec_bf16 *va0, *v_x;
  vec_f32 temp0 = { 0, 0, 0, 0 };
  vec_bf16 zero = { 0, 0, 0, 0, 0, 0, 0, 0 };
  vec_f32 inp[2];

  a0 = ap;
  va0 = (vec_bf16 *)a0;
  v_x = (vec_bf16 *)x;
  BLASLONG n8 = n / 8;
  BLASLONG i = 0;

  for (; i < n8; i++) {
    vec_load_vec2(v_x, i, inp, zero);

    temp0 += vec_load_mult(&va0[i], inp, zero);
  }

  n &= 7;
  if (n > 4) {
    vec_loadN_vec2(v_x, i, inp, n, zero);

    temp0 += vec_loadN_mult(&va0[i], inp, n, zero);
  } else if (n) {
    inp[0] = vec_loadNHi_vec(v_x, i, n, zero);

    temp0 += vec_loadNHi_mult(&va0[i], inp[0], n, zero);
  }

  y[0] = (alpha * (temp0[0] + temp0[1] + temp0[2] + temp0[3])) + (beta * y[0]);
}

static void BF16GEMV_T_VSX_2(BLASLONG n, BLASLONG lda, IFLOAT *ap, IFLOAT *x, FLOAT *y, FLOAT alpha, FLOAT beta)
{
  IFLOAT *a0, *a1;
  vec_bf16 *va0, *va1, *v_x;
  vec_f32 temp0 = { 0, 0, 0, 0 };
  vec_f32 temp1 = { 0, 0, 0, 0 };
  vec_bf16 zero = { 0, 0, 0, 0, 0, 0, 0, 0 };
  vec_f32 inp[2];

  a0 = ap;
  a1 = ap + lda;
  va0 = (vec_bf16 *)a0;
  va1 = (vec_bf16 *)a1;
  v_x = (vec_bf16 *)x;
  BLASLONG n8 = n / 8;
  BLASLONG i = 0;

  for (; i < n8; i++) {
    vec_load_vec2(v_x, i, inp, zero);

    temp0 += vec_load_mult(&va0[i], inp, zero);
    temp1 += vec_load_mult(&va1[i], inp, zero);
  }

  n &= 7;
  if (n > 4) {
    vec_loadN_vec2(v_x, i, inp, n, zero);

    temp0 += vec_loadN_mult(&va0[i], inp, n, zero);
    temp1 += vec_loadN_mult(&va1[i], inp, n, zero);
  } else if (n) {
    inp[0] = vec_loadNHi_vec(v_x, i, n, zero);

    temp0 += vec_loadNHi_mult(&va0[i], inp[0], n, zero);
    temp1 += vec_loadNHi_mult(&va1[i], inp[0], n, zero);
  }

  y[0] = (alpha * (temp0[0] + temp0[1] + temp0[2] + temp0[3])) + (beta * y[0]);
  y[1] = (alpha * (temp1[0] + temp1[1] + temp1[2] + temp1[3])) + (beta * y[1]);
}

static void BF16GEMV_T_VSX_4(BLASLONG n, BLASLONG lda, IFLOAT *ap, IFLOAT *x, FLOAT *y, FLOAT alpha, FLOAT beta)
{
  IFLOAT *a0, *a1, *a2, *a3;
  vec_bf16 *va0, *va1, *va2, *va3, *v_x;
  vec_f32 temp0 = { 0, 0, 0, 0 };
  vec_f32 temp1 = { 0, 0, 0, 0 };
  vec_f32 temp2 = { 0, 0, 0, 0 };
  vec_f32 temp3 = { 0, 0, 0, 0 };
  vec_bf16 zero = { 0, 0, 0, 0, 0, 0, 0, 0 };
  vec_f32 inp[2];

  a0 = ap;
  a1 = ap + lda;
  a2 = a1 + lda;
  a3 = a2 + lda;
  va0 = (vec_bf16 *)a0;
  va1 = (vec_bf16 *)a1;
  va2 = (vec_bf16 *)a2;
  va3 = (vec_bf16 *)a3;
  v_x = (vec_bf16 *)x;
  BLASLONG n8 = n / 8;
  BLASLONG i = 0;

  for (; i < n8; i++) {
    vec_load_vec2(v_x, i, inp, zero);

    temp0 += vec_load_mult(&va0[i], inp, zero);
    temp1 += vec_load_mult(&va1[i], inp, zero);
    temp2 += vec_load_mult(&va2[i], inp, zero);
    temp3 += vec_load_mult(&va3[i], inp, zero);
  }

  n &= 7;
  if (n > 4) {
    vec_loadN_vec2(v_x, i, inp, n, zero);

    temp0 += vec_loadN_mult(&va0[i], inp, n, zero);
    temp1 += vec_loadN_mult(&va1[i], inp, n, zero);
    temp2 += vec_loadN_mult(&va2[i], inp, n, zero);
    temp3 += vec_loadN_mult(&va3[i], inp, n, zero);
  } else if (n) {
    inp[0] = vec_loadNHi_vec(v_x, i, n, zero);

    temp0 += vec_loadNHi_mult(&va0[i], inp[0], n, zero);
    temp1 += vec_loadNHi_mult(&va1[i], inp[0], n, zero);
    temp2 += vec_loadNHi_mult(&va2[i], inp[0], n, zero);
    temp3 += vec_loadNHi_mult(&va3[i], inp[0], n, zero);
  }

  vec_f32 t0, t1, t2, t3;
  vec_f32 a = { alpha, alpha, alpha, alpha };
  vec_f32 b = { beta, beta, beta, beta };
  vec_f32 *v_y = (vec_f32 *) y;

  t0 = vec_mergeh(temp0, temp2);
  t1 = vec_mergel(temp0, temp2);
  t2 = vec_mergeh(temp1, temp3);
  t3 = vec_mergel(temp1, temp3);
  temp0 = vec_mergeh(t0, t2);
  temp1 = vec_mergel(t0, t2);
  temp2 = vec_mergeh(t1, t3);
  temp3 = vec_mergel(t1, t3);
  temp0 += temp1 + temp2 + temp3;

  v_y[0] = (a * temp0) + (b * v_y[0]);
}

#ifdef USE_BFGEMV_8_T_VSX
static void BF16GEMV_T_VSX_8(BLASLONG n, BLASLONG lda, IFLOAT *ap, IFLOAT *x, FLOAT *y, FLOAT alpha, FLOAT beta)
{
  IFLOAT *a0, *a1, *a2, *a3, *a4, *a5, *a6, *a7;
  vec_bf16 *va0, *va1, *va2, *va3, *va4, *va5, *va6, *va7, *v_x;
  vec_f32 temp0 = { 0, 0, 0, 0 };
  vec_f32 temp1 = { 0, 0, 0, 0 };
  vec_f32 temp2 = { 0, 0, 0, 0 };
  vec_f32 temp3 = { 0, 0, 0, 0 };
  vec_f32 temp4 = { 0, 0, 0, 0 };
  vec_f32 temp5 = { 0, 0, 0, 0 };
  vec_f32 temp6 = { 0, 0, 0, 0 };
  vec_f32 temp7 = { 0, 0, 0, 0 };
  vec_bf16 zero = { 0, 0, 0, 0, 0, 0, 0, 0 };
  vec_f32 inp[2];

  a0 = ap;
  a1 = ap + lda;
  a2 = a1 + lda;
  a3 = a2 + lda;
  a4 = a3 + lda;
  a5 = a4 + lda;
  a6 = a5 + lda;
  a7 = a6 + lda;
  va0 = (vec_bf16 *)a0;
  va1 = (vec_bf16 *)a1;
  va2 = (vec_bf16 *)a2;
  va3 = (vec_bf16 *)a3;
  va4 = (vec_bf16 *)a4;
  va5 = (vec_bf16 *)a5;
  va6 = (vec_bf16 *)a6;
  va7 = (vec_bf16 *)a7;
  v_x = (vec_bf16 *)x;
  BLASLONG n8 = n / 8;
  BLASLONG i = 0;

  for (; i < n8; i++) {
    vec_load_vec2(v_x, i, inp, zero);

    temp0 += vec_load_mult(&va0[i], inp, zero);
    temp1 += vec_load_mult(&va1[i], inp, zero);
    temp2 += vec_load_mult(&va2[i], inp, zero);
    temp3 += vec_load_mult(&va3[i], inp, zero);
    temp4 += vec_load_mult(&va4[i], inp, zero);
    temp5 += vec_load_mult(&va5[i], inp, zero);
    temp6 += vec_load_mult(&va6[i], inp, zero);
    temp7 += vec_load_mult(&va7[i], inp, zero);
  }

  n &= 7;
  if (n > 4) {
    vec_loadN_vec2(v_x, i, inp, n, zero);

    temp0 += vec_loadN_mult(&va0[i], inp, n, zero);
    temp1 += vec_loadN_mult(&va1[i], inp, n, zero);
    temp2 += vec_loadN_mult(&va2[i], inp, n, zero);
    temp3 += vec_loadN_mult(&va3[i], inp, n, zero);
    temp4 += vec_loadN_mult(&va4[i], inp, n, zero);
    temp5 += vec_loadN_mult(&va5[i], inp, n, zero);
    temp6 += vec_loadN_mult(&va6[i], inp, n, zero);
    temp7 += vec_loadN_mult(&va7[i], inp, n, zero);
  } else if (n) {
    inp[0] = vec_loadNHi_vec(v_x, i, n, zero);

    temp0 += vec_loadNHi_mult(&va0[i], inp[0], n, zero);
    temp1 += vec_loadNHi_mult(&va1[i], inp[0], n, zero);
    temp2 += vec_loadNHi_mult(&va2[i], inp[0], n, zero);
    temp3 += vec_loadNHi_mult(&va3[i], inp[0], n, zero);
    temp4 += vec_loadNHi_mult(&va4[i], inp[0], n, zero);
    temp5 += vec_loadNHi_mult(&va5[i], inp[0], n, zero);
    temp6 += vec_loadNHi_mult(&va6[i], inp[0], n, zero);
    temp7 += vec_loadNHi_mult(&va7[i], inp[0], n, zero);
  }

  vec_f32 t0, t1, t2, t3;
  vec_f32 a = { alpha, alpha, alpha, alpha };
  vec_f32 b = { beta, beta, beta, beta };
  vec_f32 *v_y = (vec_f32 *) y;

  t0 = vec_mergeh(temp0, temp2);
  t1 = vec_mergel(temp0, temp2);
  t2 = vec_mergeh(temp1, temp3);
  t3 = vec_mergel(temp1, temp3);
  temp0 = vec_mergeh(t0, t2);
  temp1 = vec_mergel(t0, t2);
  temp2 = vec_mergeh(t1, t3);
  temp3 = vec_mergel(t1, t3);
  temp0 += temp1 + temp2 + temp3;

  t0 = vec_mergeh(temp4, temp6);
  t1 = vec_mergel(temp4, temp6);
  t2 = vec_mergeh(temp5, temp7);
  t3 = vec_mergel(temp5, temp7);
  temp4 = vec_mergeh(t0, t2);
  temp5 = vec_mergel(t0, t2);
  temp6 = vec_mergeh(t1, t3);
  temp7 = vec_mergel(t1, t3);
  temp4 += temp5 + temp6 + temp7;

  vec_load_pair(inp, v_y);
  inp[0] = (a * temp0) + (b * inp[0]);
  inp[1] = (a * temp4) + (b * inp[1]);
  vec_store_pair(v_y, inp);
}
#endif

#include "sbgemv_t.c"
#endif

