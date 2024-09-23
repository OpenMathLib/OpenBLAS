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

#ifndef SBGEMV_COMMON_MMA_C
#define SBGEMV_COMMON_MMA_C
#include "sbgemv_common.c"

FORCEINLINE void vec_load_mult_mma(__vector_quad *out, vec_bf16 *in, vec_bf16 inp)
{
  vec_bf16 in0 = (vec_bf16)vec_load_vec(in);

  __builtin_mma_xvbf16ger2pp(out, (vec_uc8)in0, (vec_uc8)inp);
}

FORCEINLINE void vec_load_mult2_mma(__vector_quad *out, vec_bf16 *in, vec_bf16 *inp)
{
  vec_bf16 in0[2];

  vec_load_pair((vec_f32 *)in0, (vec_f32 *)in);

  __builtin_mma_xvbf16ger2pp(out, (vec_uc8)in0[0], (vec_uc8)inp[0]);
  __builtin_mma_xvbf16ger2pp(out, (vec_uc8)in0[1], (vec_uc8)inp[1]);
}

FORCEINLINE void vec_loadN_mult_mma(__vector_quad *out, vec_bf16 *in, vec_bf16 inp, BLASLONG n)
{
  vec_bf16 in0 = vec_loadN(in, n);

  __builtin_mma_xvbf16ger2pp(out, (vec_uc8)in0, (vec_uc8)inp);
}

FORCEINLINE void vec_mult1_mma(__vector_quad *out, vec_bf16 in0, vec_bf16 inp)
{
  vec_bf16 in00 = vec_mergeh(in0, in0);

  __builtin_mma_xvbf16ger2(out, (vec_uc8)inp, (vec_uc8)in00);
}

FORCEINLINE void vec_mult2_mma(__vector_quad *out, vec_bf16 in0, vec_bf16 inp)
{
  vec_bf16 in01 = vec_mergel(in0, in0);

  vec_mult1_mma(&out[0], in0, inp);

  __builtin_mma_xvbf16ger2(&out[1], (vec_uc8)inp, (vec_uc8)in01);
}

FORCEINLINE void vec_mult4_mma(__vector_quad *out, vec_bf16 *in0, vec_bf16 inp)
{
  vec_mult2_mma(out + 0, in0[0], inp);
  vec_mult2_mma(out + 2, in0[1], inp);
}

FORCEINLINE void vec_loadN_mult11_mma(__vector_quad *out, vec_bf16 *in, vec_bf16 inp, BLASLONG n)
{
  vec_bf16 in0 = vec_loadN(in, n);

  vec_mult1_mma(out, in0, inp);
}

FORCEINLINE void vec_loadN_mult12_mma(__vector_quad *out, vec_bf16 *in, vec_bf16 inp, BLASLONG n)
{
  vec_bf16 in0 = vec_loadN(in, n);

  vec_mult2_mma(out, in0, inp);
}

FORCEINLINE void vec_load_mult12_mma(__vector_quad *out, vec_bf16 *in, vec_bf16 inp)
{
  vec_bf16 in0 = (vec_bf16)vec_load_vec(in);

  vec_mult2_mma(out, in0, inp);
}

FORCEINLINE void vec_load_mult18_mma(__vector_quad *out, vec_bf16 *in, vec_bf16 inp)
{
  vec_bf16 in0[4];

  vec_load_pair((vec_f32 *)(in0 + 0), (vec_f32 *)(in + 0));
  vec_load_pair((vec_f32 *)(in0 + 2), (vec_f32 *)(in + 2));

  vec_mult4_mma(&out[0], in0 + 0, inp);
  vec_mult4_mma(&out[4], in0 + 2, inp);
}

FORCEINLINE void vec_reduce1_mma(__vector_quad *out, vec_f32 *temp, vec_f32 v_alpha, vec_f32 *vy0)
{
  __builtin_mma_disassemble_acc((void*)temp, &out[0]);

  vy0[0] += (temp[0] * v_alpha);
}

FORCEINLINE void vec_reduce2_mma(__vector_quad *out, vec_f32 *temp, vec_f32 v_alpha, vec_f32 *vy0)
{
  vec_reduce1_mma(&out[0], &temp[0], v_alpha, &vy0[0]);
  vec_reduce1_mma(&out[1], &temp[4], v_alpha, &vy0[1]);
}

FORCEINLINE void vec_reduce8_mma(__vector_quad *out, vec_f32 *temp, vec_f32 v_alpha, vec_f32 *vy0)
{
  vec_reduce2_mma(&out[0], &temp[0],  v_alpha, vy0 + 0);
  vec_reduce2_mma(&out[2], &temp[8],  v_alpha, vy0 + 2);
  vec_reduce2_mma(&out[4], &temp[16], v_alpha, vy0 + 4);
  vec_reduce2_mma(&out[6], &temp[24], v_alpha, vy0 + 6);
}

FORCEINLINE void vec_mult11a_mma(__vector_quad *out, vec_bf16 in0, vec_bf16 in1, vec_bf16 inp)
{
  vec_bf16 in00 = vec_mergeh(in0, in1);

  __builtin_mma_xvbf16ger2(out, (vec_uc8)inp, (vec_uc8)in00);
}

FORCEINLINE void vec_mult2a_mma(__vector_quad *out, vec_bf16 in0, vec_bf16 in1, vec_bf16 inp)
{
  vec_bf16 in01 = vec_mergel(in0, in1);

  vec_mult11a_mma(&out[0], in0, in1, inp);

  __builtin_mma_xvbf16ger2(&out[1], (vec_uc8)inp, (vec_uc8)in01);
}

FORCEINLINE void vec_mult4a_mma(__vector_quad *out, vec_bf16 *in0, vec_bf16 *in1, vec_bf16 inp)
{
  vec_mult2a_mma(out + 0, in0[0], in1[0], inp);
  vec_mult2a_mma(out + 2, in0[1], in1[1], inp);
}

FORCEINLINE void vec_loadN_mult11a_mma(__vector_quad *out, vec_bf16 *ina, vec_bf16 *inb, vec_bf16 inp, BLASLONG n)
{
  vec_bf16 in0 = vec_loadN(ina, n);
  vec_bf16 in1 = vec_loadN(inb, n);

  vec_mult11a_mma(out, in0, in1, inp);
}

FORCEINLINE void vec_load_mult22a_mma(__vector_quad *out, vec_bf16 *ina, vec_bf16 *inb, vec_bf16 inp)
{
  vec_bf16 in0 = (vec_bf16)vec_load_vec(ina);
  vec_bf16 in1 = (vec_bf16)vec_load_vec(inb);

  vec_mult2a_mma(out, in0, in1, inp);
}

FORCEINLINE void vec_load_mult28a_mma(__vector_quad *out, vec_bf16 *ina, vec_bf16 *inb, vec_bf16 inp)
{
  vec_bf16 in0[4], in1[4];

  vec_load_pair((vec_f32 *)(in0 + 0), (vec_f32 *)(ina + 0));
  vec_load_pair((vec_f32 *)(in1 + 0), (vec_f32 *)(inb + 0));
  vec_load_pair((vec_f32 *)(in0 + 2), (vec_f32 *)(ina + 2));
  vec_load_pair((vec_f32 *)(in1 + 2), (vec_f32 *)(inb + 2));

  vec_mult4a_mma(&out[0], in0 + 0, in1 + 0, inp);
  vec_mult4a_mma(&out[4], in0 + 2, in1 + 2, inp);
}

FORCEINLINE void vec_loadN_mult22a_mma(__vector_quad *out, vec_bf16 *ina, vec_bf16 *inb, vec_bf16 inp, BLASLONG n)
{
  vec_bf16 in0 = vec_loadN(ina, n);
  vec_bf16 in1 = vec_loadN(inb, n);

  vec_mult2a_mma(out, in0, in1, inp);
}

FORCEINLINE void vec_mult11b_mma(__vector_quad *out, vec_bf16 in0, vec_bf16 in1, vec_bf16 inp)
{
  vec_bf16 in00 = vec_mergeh(in0, in1);

  __builtin_mma_xvbf16ger2pp(out, (vec_uc8)inp, (vec_uc8)in00);
}

FORCEINLINE void vec_mult2b_mma(__vector_quad *out, vec_bf16 in0, vec_bf16 in1, vec_bf16 inp)
{
  vec_bf16 in01 = vec_mergel(in0, in1);

  vec_mult11b_mma(&out[0], in0, in1, inp);

  __builtin_mma_xvbf16ger2pp(&out[1], (vec_uc8)inp, (vec_uc8)in01);
}

FORCEINLINE void vec_mult4b_mma(__vector_quad *out, vec_bf16 *in0, vec_bf16 *in1, vec_bf16 inp)
{
  vec_mult2b_mma(out + 0, in0[0], in1[0], inp);
  vec_mult2b_mma(out + 2, in0[1], in1[1], inp);
}

FORCEINLINE void vec_loadN_mult11b_mma(__vector_quad *out, vec_bf16 *ina, vec_bf16 *inb, vec_bf16 inp, BLASLONG n)
{
  vec_bf16 in0 = vec_loadN(ina, n);
  vec_bf16 in1 = vec_loadN(inb, n);

  vec_mult11b_mma(out, in0, in1, inp);
}

FORCEINLINE void vec_load_mult22b_mma(__vector_quad *out, vec_bf16 *ina, vec_bf16 *inb, vec_bf16 inp)
{
  vec_bf16 in0 = (vec_bf16)vec_load_vec(ina);
  vec_bf16 in1 = (vec_bf16)vec_load_vec(inb);

  vec_mult2b_mma(out, in0, in1, inp);
}

FORCEINLINE void vec_load_mult28b_mma(__vector_quad *out, vec_bf16 *ina, vec_bf16 *inb, vec_bf16 inp)
{
  vec_bf16 in0[4], in1[4];

  vec_load_pair((vec_f32 *)(in0 + 0), (vec_f32 *)(ina + 0));
  vec_load_pair((vec_f32 *)(in1 + 0), (vec_f32 *)(inb + 0));
  vec_load_pair((vec_f32 *)(in0 + 2), (vec_f32 *)(ina + 2));
  vec_load_pair((vec_f32 *)(in1 + 2), (vec_f32 *)(inb + 2));

  vec_mult4b_mma(&out[0], in0 + 0, in1 + 0, inp);
  vec_mult4b_mma(&out[4], in0 + 2, in1 + 2, inp);
}

FORCEINLINE void vec_loadN_mult22b_mma(__vector_quad *out, vec_bf16 *ina, vec_bf16 *inb, vec_bf16 inp, BLASLONG n)
{
  vec_bf16 in0 = vec_loadN(ina, n);
  vec_bf16 in1 = vec_loadN(inb, n);

  vec_mult2b_mma(out, in0, in1, inp);
}

FORCEINLINE void vec_load4_pair(vec_f32 *vy0, vec_f32 *v_y)
{
  vec_load_pair(vy0 + 0, v_y + 0);
  vec_load_pair(vy0 + 2, v_y + 2);
  vec_load_pair(vy0 + 4, v_y + 4);
  vec_load_pair(vy0 + 6, v_y + 6);
}

FORCEINLINE void vec_store4_pair(vec_f32 *v_y, vec_f32 *vy0)
{
  vec_store_pair(v_y + 0, vy0 + 0);
  vec_store_pair(v_y + 2, vy0 + 2);
  vec_store_pair(v_y + 4, vy0 + 4);
  vec_store_pair(v_y + 6, vy0 + 6);
}

#endif
