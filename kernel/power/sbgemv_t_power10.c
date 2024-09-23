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

#ifndef SBGEMV_T_MMA_C
#define SBGEMV_T_MMA_C

#define USE_BFGEMV_T_MMA

#ifdef USE_BFGEMV_T_MMA
#include "sbgemv_common_power10.c"

#ifndef BF16GEMV_T_X
#define BF16GEMV_T_X
#define BF16GEMV_T_8 BF16GEMV_T_MMA_8
#define BF16GEMV_T_4 BF16GEMV_T_MMA_4
#define BF16GEMV_T_2 BF16GEMV_T_MMA_2
#define BF16GEMV_T_1 BF16GEMV_T_MMA_1
#endif

#define USE_BFGEMV_8_T_MMA

static void BF16GEMV_T_MMA_1(BLASLONG n, BLASLONG lda, IFLOAT *ap, IFLOAT *x, FLOAT *y, FLOAT alpha, FLOAT beta)
{
  IFLOAT *a0;
  vec_bf16 *va0, *v_x;
  __vector_quad temp0;
  vec_f32 temp00[4];
  vec_bf16 inp[2];

  __builtin_mma_xxsetaccz(&temp0);

  a0 = ap;
  va0 = (vec_bf16 *)a0;
  v_x = (vec_bf16 *)x;
  BLASLONG n8 = n / 8;
  BLASLONG i = 0;

  for (; i + 2 <= n8; i += 2) {
    vec_load_pair((vec_f32 *)inp, (vec_f32 *)&v_x[i]);

    vec_load_mult2_mma(&temp0, &va0[i + 0], inp);
  }

  if (n8 & 1) {
    inp[0] = (vec_bf16)vec_load_vec(&v_x[i]);

    vec_load_mult_mma(&temp0, &va0[i], inp[0]);

    i++;
  }

  n &= 7;
  if (n) {
    inp[0] = vec_loadN(&v_x[i], n);

    vec_loadN_mult_mma(&temp0, &va0[i], inp[0], n);
  }

  __builtin_mma_disassemble_acc((void*)temp00, &temp0);

  y[0] = (alpha * (temp00[0][0] + temp00[1][1] + temp00[2][2] + temp00[3][3])) + (beta * y[0]);
}

static void BF16GEMV_T_MMA_2(BLASLONG n, BLASLONG lda, IFLOAT *ap, IFLOAT *x, FLOAT *y, FLOAT alpha, FLOAT beta)
{
  IFLOAT *a0, *a1;
  vec_bf16 *va0, *va1, *v_x;
  __vector_quad temp0, temp1;
  vec_f32 temp00[4], temp01[4];
  vec_bf16 inp[2];

  __builtin_mma_xxsetaccz(&temp0);
  __builtin_mma_xxsetaccz(&temp1);

  a0 = ap;
  a1 = ap + lda;
  va0 = (vec_bf16 *)a0;
  va1 = (vec_bf16 *)a1;
  v_x = (vec_bf16 *)x;
  BLASLONG n8 = n / 8;
  BLASLONG i = 0;

  for (; i + 2 <= n8; i += 2) {
    vec_load_pair((vec_f32 *)inp, (vec_f32 *)&v_x[i]);

    vec_load_mult2_mma(&temp0, &va0[i + 0], inp);
    vec_load_mult2_mma(&temp1, &va1[i + 0], inp);
  }

  if (n8 & 1) {
    inp[0] = (vec_bf16)vec_load_vec(&v_x[i]);

    vec_load_mult_mma(&temp0, &va0[i], inp[0]);
    vec_load_mult_mma(&temp1, &va1[i], inp[0]);

    i++;
  }

  n &= 7;
  if (n) {
    inp[0] = vec_loadN(&v_x[i], n);

    vec_loadN_mult_mma(&temp0, &va0[i], inp[0], n);
    vec_loadN_mult_mma(&temp1, &va1[i], inp[0], n);
  }

  __builtin_mma_disassemble_acc((void*)temp00, &temp0);
  __builtin_mma_disassemble_acc((void*)temp01, &temp1);

  y[0] = (alpha * (temp00[0][0] + temp00[1][1] + temp00[2][2] + temp00[3][3])) + (beta * y[0]);
  y[1] = (alpha * (temp01[0][0] + temp01[1][1] + temp01[2][2] + temp01[3][3])) + (beta * y[1]);
}

static void BF16GEMV_T_MMA_4(BLASLONG n, BLASLONG lda, IFLOAT *ap, IFLOAT *x, FLOAT *y, FLOAT alpha, FLOAT beta)
{
  IFLOAT *a0, *a1, *a2, *a3;
  vec_bf16 *va0, *va1, *va2, *va3, *v_x;
  __vector_quad temp0, temp1, temp2, temp3;
  vec_f32 temp00[4], temp01[4], temp02[4], temp03[4];
  vec_bf16 inp[2];

  __builtin_mma_xxsetaccz(&temp0);
  __builtin_mma_xxsetaccz(&temp1);
  __builtin_mma_xxsetaccz(&temp2);
  __builtin_mma_xxsetaccz(&temp3);

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

  for (; i + 2 <= n8; i += 2) {
    vec_load_pair((vec_f32 *)inp, (vec_f32 *)&v_x[i]);

    vec_load_mult2_mma(&temp0, &va0[i + 0], inp);
    vec_load_mult2_mma(&temp1, &va1[i + 0], inp);
    vec_load_mult2_mma(&temp2, &va2[i + 0], inp);
    vec_load_mult2_mma(&temp3, &va3[i + 0], inp);
  }

  if (n8 & 1) {
    inp[0] = (vec_bf16)vec_load_vec(&v_x[i]);

    vec_load_mult_mma(&temp0, &va0[i], inp[0]);
    vec_load_mult_mma(&temp1, &va1[i], inp[0]);
    vec_load_mult_mma(&temp2, &va2[i], inp[0]);
    vec_load_mult_mma(&temp3, &va3[i], inp[0]);

    i++;
  }

  n &= 7;
  if (n) {
    inp[0] = vec_loadN(&v_x[i], n);

    vec_loadN_mult_mma(&temp0, &va0[i], inp[0], n);
    vec_loadN_mult_mma(&temp1, &va1[i], inp[0], n);
    vec_loadN_mult_mma(&temp2, &va2[i], inp[0], n);
    vec_loadN_mult_mma(&temp3, &va3[i], inp[0], n);
  }

  __builtin_mma_disassemble_acc((void*)temp00, &temp0);
  __builtin_mma_disassemble_acc((void*)temp01, &temp1);
  __builtin_mma_disassemble_acc((void*)temp02, &temp2);
  __builtin_mma_disassemble_acc((void*)temp03, &temp3);

  vec_f32 t0, t1, t2, t3, t4, t5, t6, t7;
  vec_f32 a = { alpha, alpha, alpha, alpha };
  vec_f32 b = { beta, beta, beta, beta };
  vec_f32 *v_y = (vec_f32 *) y;

  t0 = vec_mergeh(temp00[0], temp01[0]);
  t1 = vec_mergeh(temp02[0], temp03[0]);
  t2 = vec_mergeo(temp00[1], temp01[1]);
  t3 = vec_mergeo(temp02[1], temp03[1]);
  t4 = vec_mergel(temp00[2], temp01[2]);
  t5 = vec_mergel(temp02[2], temp03[2]);
  t6 = vec_mergeo(temp00[3], temp01[3]);
  t7 = vec_mergeo(temp02[3], temp03[3]);
  t0 = vec_xxpermdi(t0, t1, 0);
  t2 = vec_xxpermdi(t2, t3, 0);
  t4 = vec_xxpermdi(t4, t5, 0);
  t6 = vec_xxpermdi(t6, t7, 3);

  t0 += t2 + t4 + t6;

  v_y[0] = (a * t0) + (b * v_y[0]);
}

#ifdef USE_BFGEMV_8_T_MMA
static void BF16GEMV_T_MMA_8(BLASLONG n, BLASLONG lda, IFLOAT *ap, IFLOAT *x, FLOAT *y, FLOAT alpha, FLOAT beta)
{
  IFLOAT *a0, *a1, *a2, *a3, *a4, *a5, *a6, *a7;
  vec_bf16 *va0, *va1, *va2, *va3, *va4, *va5, *va6, *va7, *v_x;
  __vector_quad temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7;
  vec_f32 temp00[4], temp01[4], temp02[4], temp03[4], temp04[4], temp05[4], temp06[4], temp07[4];
  vec_bf16 inp[2];

  __builtin_mma_xxsetaccz(&temp0);
  __builtin_mma_xxsetaccz(&temp1);
  __builtin_mma_xxsetaccz(&temp2);
  __builtin_mma_xxsetaccz(&temp3);
  __builtin_mma_xxsetaccz(&temp4);
  __builtin_mma_xxsetaccz(&temp5);
  __builtin_mma_xxsetaccz(&temp6);
  __builtin_mma_xxsetaccz(&temp7);

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

  for (; i + 2 <= n8; i += 2) {
    vec_load_pair((vec_f32 *)inp, (vec_f32 *)&v_x[i]);

    vec_load_mult2_mma(&temp0, &va0[i + 0], inp);
    vec_load_mult2_mma(&temp1, &va1[i + 0], inp);
    vec_load_mult2_mma(&temp2, &va2[i + 0], inp);
    vec_load_mult2_mma(&temp3, &va3[i + 0], inp);
    vec_load_mult2_mma(&temp4, &va4[i + 0], inp);
    vec_load_mult2_mma(&temp5, &va5[i + 0], inp);
    vec_load_mult2_mma(&temp6, &va6[i + 0], inp);
    vec_load_mult2_mma(&temp7, &va7[i + 0], inp);
  }

  if (n8 & 1) {
    inp[0] = (vec_bf16)vec_load_vec(&v_x[i]);

    vec_load_mult_mma(&temp0, &va0[i], inp[0]);
    vec_load_mult_mma(&temp1, &va1[i], inp[0]);
    vec_load_mult_mma(&temp2, &va2[i], inp[0]);
    vec_load_mult_mma(&temp3, &va3[i], inp[0]);
    vec_load_mult_mma(&temp4, &va4[i], inp[0]);
    vec_load_mult_mma(&temp5, &va5[i], inp[0]);
    vec_load_mult_mma(&temp6, &va6[i], inp[0]);
    vec_load_mult_mma(&temp7, &va7[i], inp[0]);

    i++;
  }

  n &= 7;
  if (n) {
    inp[0] = vec_loadN(&v_x[i], n);

    vec_loadN_mult_mma(&temp0, &va0[i], inp[0], n);
    vec_loadN_mult_mma(&temp1, &va1[i], inp[0], n);
    vec_loadN_mult_mma(&temp2, &va2[i], inp[0], n);
    vec_loadN_mult_mma(&temp3, &va3[i], inp[0], n);
    vec_loadN_mult_mma(&temp4, &va4[i], inp[0], n);
    vec_loadN_mult_mma(&temp5, &va5[i], inp[0], n);
    vec_loadN_mult_mma(&temp6, &va6[i], inp[0], n);
    vec_loadN_mult_mma(&temp7, &va7[i], inp[0], n);
  }

  __builtin_mma_disassemble_acc((void*)temp00, &temp0);
  __builtin_mma_disassemble_acc((void*)temp01, &temp1);
  __builtin_mma_disassemble_acc((void*)temp02, &temp2);
  __builtin_mma_disassemble_acc((void*)temp03, &temp3);
  __builtin_mma_disassemble_acc((void*)temp04, &temp4);
  __builtin_mma_disassemble_acc((void*)temp05, &temp5);
  __builtin_mma_disassemble_acc((void*)temp06, &temp6);
  __builtin_mma_disassemble_acc((void*)temp07, &temp7);

  vec_f32 t0, t1, t2, t3, t4, t5, t6, t7, t10, t11, t12, t13, t14, t15, t16, t17;
  vec_f32 a = { alpha, alpha, alpha, alpha };
  vec_f32 b = { beta, beta, beta, beta };
  vec_f32 *v_y = (vec_f32 *) y;

  t0 = vec_mergeh(temp00[0], temp01[0]);
  t1 = vec_mergeh(temp02[0], temp03[0]);
  t2 = vec_mergeo(temp00[1], temp01[1]);
  t3 = vec_mergeo(temp02[1], temp03[1]);
  t4 = vec_mergel(temp00[2], temp01[2]);
  t5 = vec_mergel(temp02[2], temp03[2]);
  t6 = vec_mergeo(temp00[3], temp01[3]);
  t7 = vec_mergeo(temp02[3], temp03[3]);
  t0 = vec_xxpermdi(t0, t1, 0);
  t2 = vec_xxpermdi(t2, t3, 0);
  t4 = vec_xxpermdi(t4, t5, 0);
  t6 = vec_xxpermdi(t6, t7, 3);

  t0 += t2 + t4 + t6;

  t10 = vec_mergeh(temp04[0], temp05[0]);
  t11 = vec_mergeh(temp06[0], temp07[0]);
  t12 = vec_mergeo(temp04[1], temp05[1]);
  t13 = vec_mergeo(temp06[1], temp07[1]);
  t14 = vec_mergel(temp04[2], temp05[2]);
  t15 = vec_mergel(temp06[2], temp07[2]);
  t16 = vec_mergeo(temp04[3], temp05[3]);
  t17 = vec_mergeo(temp06[3], temp07[3]);
  t10 = vec_xxpermdi(t10, t11, 0);
  t12 = vec_xxpermdi(t12, t13, 0);
  t14 = vec_xxpermdi(t14, t15, 0);
  t16 = vec_xxpermdi(t16, t17, 3);

  t10 += t12 + t14 + t16;

  vec_f32 inp2[2];
  vec_load_pair(inp2, v_y);
  inp2[0] = (a * t0) + (b * inp2[0]);
  inp2[1] = (a * t10) + (b * inp2[1]);
  vec_store_pair(v_y, inp2);
}
#endif

#include "sbgemv_t.c"
#else
#include "sbgemv_t_vsx.c"
#endif
#endif

