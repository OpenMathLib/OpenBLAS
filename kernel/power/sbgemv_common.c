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

#ifndef SBGEMV_COMMON_C
#define SBGEMV_COMMON_C
#include "common.h"

#include <altivec.h>

#define FORCEINLINE      inline __attribute__((always_inline))

#ifdef __clang__
#define uint16_t         unsigned short
#define uint32_t         unsigned int
#define uint64_t         unsigned long long
#endif

#ifdef _ARCH_PWR10
#ifdef __has_builtin
#if !__has_builtin(__builtin_vsx_assemble_pair)
#define __builtin_vsx_assemble_pair __builtin_mma_assemble_pair
#endif
#if !__has_builtin(__builtin_vsx_disassemble_pair)
#define __builtin_vsx_disassemble_pair __builtin_mma_disassemble_pair
#endif
#endif

#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
#define __builtin_vsx_assemble_pair2(vp0, v0, v1) __builtin_vsx_assemble_pair(vp0, v1, v0)
#else
#define __builtin_vsx_assemble_pair2(vp0, v0, v1) __builtin_vsx_assemble_pair(vp0, v0, v1)
#endif

#define USE_VECTOR_PAIRS
#endif

typedef __vector IFLOAT        vec_bf16;
typedef __vector FLOAT         vec_f32;
typedef __vector unsigned char vec_uc8;

#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
#define BF16_HI(data, zero)     (vec_f32)vec_mergeh(data, zero)
#define BF16_LO(data, zero)     (vec_f32)vec_mergel(data, zero)
#else
#define BF16_HI(data, zero)     (vec_f32)vec_mergeh(zero, data)
#define BF16_LO(data, zero)     (vec_f32)vec_mergel(zero, data)
#endif

FORCEINLINE vec_uc8 vec_load_vec(void *src)
{
  return vec_xl(0, (unsigned char *)(src));
}

FORCEINLINE void vec_load_pair(vec_f32 *dst, vec_f32 *src)
{
#ifdef USE_VECTOR_PAIRS
  __vector_pair vy0p;
  vy0p = *(__vector_pair *)(src);
  __builtin_vsx_disassemble_pair((void *)(dst), &vy0p);
#else
  dst[0] = src[0];
  dst[1] = src[1];
#endif
}

FORCEINLINE void vec_store_pair(vec_f32 *dst, vec_f32 *src)
{
#ifdef USE_VECTOR_PAIRS
  __vector_pair vy0p;
  __builtin_vsx_assemble_pair2(&vy0p, (vec_uc8)src[1], (vec_uc8)src[0]);
  *(__vector_pair *)(dst) = vy0p;
#else
  dst[0] = src[0];
  dst[1] = src[1];
#endif
}

FORCEINLINE vec_bf16 vec_loadN(void *src, BLASLONG n)
{
  IFLOAT *src2 = (IFLOAT *)(src);
#ifdef _ARCH_PWR9
  return vec_xl_len(src2, n * sizeof(IFLOAT));
#else
  __attribute__((aligned(16))) IFLOAT data[sizeof(vec_bf16) / sizeof(IFLOAT)];
  memset(data, 0, sizeof(vec_bf16));
  if (n & 4) {
    memcpy(data, src2, sizeof(uint64_t));
  }
  if (n & 2) {
    BLASLONG n4 = n & 4;
    memcpy(data + n4, src2 + n4, sizeof(uint32_t));
  }
  if (n & 1) {
    BLASLONG n6 = n & 6;
    data[n6] = src2[n6];
  }
  return (vec_bf16)vec_load_vec(data);
#endif
}

FORCEINLINE vec_f32 vec_loadNHi(void *src, BLASLONG n, vec_bf16 zero)
{
  vec_bf16 data = vec_loadN(src, n);
  return BF16_HI(data, zero);
}

FORCEINLINE vec_f32 vec_loadN_f32(void *src, BLASLONG n)
{
#ifndef _ARCH_PWR9
  if (n & 4) {
    return (vec_f32)vec_load_vec(src);
  }
#endif
  return (vec_f32)vec_loadN(src, n * (sizeof(FLOAT) / sizeof(IFLOAT)));
}

FORCEINLINE void vec_loadN2_f32(vec_f32 *data, vec_f32 *src, BLASLONG n)
{
  data[0] = src[0];
  data[1] = vec_loadN_f32(&src[1], n);
}

FORCEINLINE void vec_storeN_f32(vec_f32 data, void *dst, BLASLONG n)
{
  FLOAT *dst2 = (FLOAT *)(dst);
#ifdef _ARCH_PWR9
  vec_xst_len(data, dst2, n * sizeof(FLOAT));
#else
  if (n & 4) {
    vec_xst(data, 0, dst2);
    return;
  }
  __attribute__((aligned(16))) FLOAT data2[sizeof(vec_f32) / sizeof(FLOAT)];
  vec_xst(data, 0, data2);
  if (n & 2) {
    memcpy(dst2, data2, sizeof(uint64_t));
  }
  if (n & 1) {
    BLASLONG n2 = n & 2;
    dst2[n2] = data2[n2];
  }
#endif
}

FORCEINLINE void vec_storeN2_f32(vec_f32 *data, vec_f32 *dst, BLASLONG n)
{
  dst[0] = data[0];
  vec_storeN_f32(data[1], &dst[1], n);
}

FORCEINLINE vec_f32 vec_mult(vec_f32 *inp, vec_bf16 in0, vec_bf16 zero)
{
  vec_f32 v_in00 = BF16_HI(in0, zero);
  vec_f32 v_in01 = BF16_LO(in0, zero);

  return (inp[0] * v_in00) + (inp[1] * v_in01);
}

FORCEINLINE vec_f32 vec_load_mult(vec_bf16 *in, vec_f32 *inp, vec_bf16 zero)
{
  vec_bf16 in0 = (vec_bf16)vec_load_vec(in);

  return vec_mult(inp, in0, zero);
}

FORCEINLINE void vec_load_vec2(vec_bf16 *in, BLASLONG i, vec_f32 *v_x0, vec_bf16 zero)
{
  vec_bf16 inp = (vec_bf16)vec_load_vec(&in[i]);

  v_x0[0] = BF16_HI(inp, zero);
  v_x0[1] = BF16_LO(inp, zero);
}

FORCEINLINE void vec_mult2(vec_f32 v_x0, vec_bf16 in0, vec_bf16 zero, vec_f32 *vy0)
{
  vec_f32 v_in00 = BF16_HI(in0, zero);
  vec_f32 v_in01 = BF16_LO(in0, zero);

  vy0[0] += (v_x0 * v_in00);
  vy0[1] += (v_x0 * v_in01);
}

FORCEINLINE void vec_load_mult2(vec_f32 v_x0, vec_bf16 *in, vec_bf16 zero, vec_f32 *vy0)
{
  vec_bf16 in0 = (vec_bf16)vec_load_vec(in);

  vec_mult2(v_x0, in0, zero, vy0);
}

FORCEINLINE vec_f32 vec_loadN_mult(vec_bf16 *in, vec_f32 *inp, BLASLONG n, vec_bf16 zero)
{
  vec_bf16 in0 = vec_loadN(in, n);

  return vec_mult(inp, in0, zero);
}

FORCEINLINE void vec_loadN_vec2(vec_bf16 *in, BLASLONG i, vec_f32 *v_x0, BLASLONG n, vec_bf16 zero)
{
  vec_bf16 inp = vec_loadN(&in[i], n);

  v_x0[0] = BF16_HI(inp, zero);
  v_x0[1] = BF16_LO(inp, zero);
}

FORCEINLINE void vec_loadN_mult2(vec_f32 v_x0, vec_bf16 *in, BLASLONG n, vec_bf16 zero, vec_f32 *vy0)
{
  vec_bf16 in0 = vec_loadN(in, n);

  vec_mult2(v_x0, in0, zero, vy0);
}

FORCEINLINE vec_f32 vec_loadNHi_mult(vec_bf16 *in, vec_f32 v_inp0, BLASLONG n, vec_bf16 zero)
{
  vec_f32 v_in00 = vec_loadNHi(in, n, zero);

  return (v_inp0 * v_in00);
}

FORCEINLINE vec_f32 vec_loadNHi_multi2(vec_f32 v_x0, vec_bf16 *in, BLASLONG n, vec_bf16 zero)
{
  vec_f32 v_in00 = vec_loadNHi(in, n, zero);

  return (v_x0 * v_in00);
}

FORCEINLINE vec_f32 vec_loadNHi_vec(vec_bf16 *in, BLASLONG i, BLASLONG n, vec_bf16 zero)
{
  return vec_loadNHi(&in[i], n, zero);
}

FORCEINLINE void copy_x(BLASLONG n, IFLOAT *src, IFLOAT *dest, BLASLONG inc_src)
{
  for (BLASLONG i = 0; i < n; i++) {
    *dest++ = *src;
    src += inc_src;
  }
}

FORCEINLINE void copy_y_beta(BLASLONG n, FLOAT *src, FLOAT *dest, BLASLONG inc_src, FLOAT beta)
{
  if (beta == 0) {
    memset(dest, 0, sizeof(FLOAT) * n);
  } else if (beta == 1) {
    for (BLASLONG i = 0; i < n; i++) {
      *dest++ = *src;
      src += inc_src;
    }
  } else {
    for (BLASLONG i = 0; i < n; i++) {
      *dest++ = *src * beta;
      src += inc_src;
    }
  }
}

FORCEINLINE void copy_y(BLASLONG n, FLOAT *src, FLOAT *dest, BLASLONG inc_src, FLOAT beta)
{
  if (beta == 0) {
    for (BLASLONG i = 0; i < n; i++) {
      *dest = *src++;
      dest += inc_src;
    }
  } else if (beta == 1) {
    for (BLASLONG i = 0; i < n; i++) {
      *dest += *src++;
      dest += inc_src;
    }
  } else {
    for (BLASLONG i = 0; i < n; i++) {
      *dest = *src++ + (beta * *dest);
      dest += inc_src;
    }
  }
}

FORCEINLINE void add_y(BLASLONG n, FLOAT *src, FLOAT *dest, BLASLONG inc_dest)
{
  for (BLASLONG i = 0; i < n; i++) {
    *dest = *src++;
    dest += inc_dest;
  }
}
#endif
