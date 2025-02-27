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
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#include "common.h"

#include <arm_neon.h>

#define A_ELEMENT_K(m, offset_k) A[(i + (m)) + (k + offset_k) * lda]
#define A_ELEMENT(m) A_ELEMENT_K(m, 0)

#define B_ELEMENT_K(n, offset_k) B[(k + offset_k) + (j + (n)) * ldb]
#define B_ELEMENT(n) B_ELEMENT_K(n, 0)

#define C_ELEMENT(m, n) C[(i + (m)) + (j + (n)) * ldc]
#define PACK_ELEMENT_K(n, offset_k) packed_b[(k + offset_k) * 4 + n]
#define PACK_ELEMENT(n) PACK_ELEMENT_K(n, 0)

// ASIMD
#define DECLARE_RESULT_VECTOR4(m, n)                                           \
  float32x4_t result##m##n = vdupq_n_f32(0.0);
#define DECLARE_RESULT(m, n) float32_t result##m##n = 0.0;
#define BROADCAST_LOAD_A4(m, offset_k)                                         \
  float32x4_t a##m##_k##offset_k = vld1q_dup_f32(&A_ELEMENT_K(m, offset_k));
#define LOAD_A1(m, offset_k)                                                   \
  float32_t a##m##_k##offset_k = A_ELEMENT_K(m, offset_k);
#define BROADCAST_LOAD_B4(n, offset_k)                                         \
  float32x4_t b##n##_k##offset_k = vld1q_dup_f32(&B_ELEMENT_K(n, offset_k));
#define LOAD_B1(n, offset_k)                                                   \
  float32_t b##n##_k##offset_k = B_ELEMENT_K(n, offset_k);
#define VECTOR_LOAD_A4(m, offset_k)                                            \
  float32x4_t a##m##_k##offset_k = vld1q_f32(&A_ELEMENT_K(m, offset_k));
#define VECTOR_LOAD_B_K4(n, offset_k)                                          \
  float32x4_t b##k##n##_k##offset_k = vld1q_f32(&B_ELEMENT_K(n, offset_k));
#define TRANSPOSE_B4_K4(                                                       \
  n0, n1, n2, n3, offset_k0, offset_k1, offset_k2, offset_k3)                  \
  float32x4_t b##t##n0##_k##offset_k0 =                                        \
    vzip1q_f32(b##k##n0##_k##offset_k0, b##k##n1##_k##offset_k0);              \
  float32x4_t b##t##n0##_k##offset_k1 =                                        \
    vzip2q_f32(b##k##n0##_k##offset_k0, b##k##n1##_k##offset_k0);              \
  float32x4_t b##t##n0##_k##offset_k2 =                                        \
    vzip1q_f32(b##k##n2##_k##offset_k0, b##k##n3##_k##offset_k0);              \
  float32x4_t b##t##n0##_k##offset_k3 =                                        \
    vzip2q_f32(b##k##n2##_k##offset_k0, b##k##n3##_k##offset_k0);              \
  float32x4_t b##n0##_k##offset_k0 = vreinterpretq_f32_f64(                    \
    vzip1q_f64(vreinterpretq_f64_f32(b##t##n0##_k##offset_k0),                 \
               vreinterpretq_f64_f32(b##t##n0##_k##offset_k2)));               \
  float32x4_t b##n0##_k##offset_k1 = vreinterpretq_f32_f64(                    \
    vzip2q_f64(vreinterpretq_f64_f32(b##t##n0##_k##offset_k0),                 \
               vreinterpretq_f64_f32(b##t##n0##_k##offset_k2)));               \
  float32x4_t b##n0##_k##offset_k2 = vreinterpretq_f32_f64(                    \
    vzip1q_f64(vreinterpretq_f64_f32(b##t##n0##_k##offset_k1),                 \
               vreinterpretq_f64_f32(b##t##n0##_k##offset_k3)));               \
  float32x4_t b##n0##_k##offset_k3 = vreinterpretq_f32_f64(                    \
    vzip2q_f64(vreinterpretq_f64_f32(b##t##n0##_k##offset_k1),                 \
               vreinterpretq_f64_f32(b##t##n0##_k##offset_k3)));

#define GATHER_LOAD_B4(n, offset_k)                                            \
  float32x4_t b##n##_k##offset_k = vdupq_n_f32(B_ELEMENT_K(n, offset_k));      \
  b##n##_k##offset_k =                                                         \
    vsetq_lane_f32(B_ELEMENT_K(n + 1, offset_k), b##n##_k##offset_k, 1);       \
  b##n##_k##offset_k =                                                         \
    vsetq_lane_f32(B_ELEMENT_K(n + 2, offset_k), b##n##_k##offset_k, 2);       \
  b##n##_k##offset_k =                                                         \
    vsetq_lane_f32(B_ELEMENT_K(n + 3, offset_k), b##n##_k##offset_k, 3);
#define VECTOR_UNPACK_B4(n, offset_k)                                          \
  float32x4_t b##n##_k##offset_k = vld1q_f32(&PACK_ELEMENT_K(n, offset_k));
#define VECTOR_PACK_B4(n, offset_k)                                            \
  vst1q_f32(&PACK_ELEMENT_K(n, offset_k), b##n##_k##offset_k);
#define PACK_B4(n, offset_k)                                                   \
  PACK_ELEMENT_K(n, offset_k) = vgetq_lane_f32(b##n##_k##offset_k, 0);
#define BROADCAST_UNPACK_B4(n, offset_k)                                       \
  float32x4_t b##n##_k##offset_k = vdupq_n_f32(PACK_ELEMENT_K(n, offset_k));
#define UPDATE_RESULT_VECTOR4(m, n, offset_k)                                  \
  result##m##n =                                                               \
    vfmaq_f32(result##m##n, a##m##_k##offset_k, b##n##_k##offset_k);
#define UPDATE_RESULT(m, n, offset_k)                                          \
  result##m##n = result##m##n + a##m##_k##offset_k * b##n##_k##offset_k;
#define UPDATE_RESULT_VECTOR4_LANE4(m, n, outer, lane, offset_k)               \
  result##m##n = vfmaq_laneq_f32(                                              \
    result##m##n, a##m##_k##offset_k, b##outer##_k##offset_k, lane);
#ifdef B0
#define VECTOR_STORE4(m, n)                                                    \
  vst1q_f32(&C_ELEMENT(m, n), vmulq_f32(result##m##n, vdupq_n_f32(alpha)));
#define STORE(m, n) C_ELEMENT(m, n) = alpha * result##m##n;
#define SCATTER_STORE4(m, n)                                                   \
  result##m##n = vmulq_f32(result##m##n, vdupq_n_f32(alpha));                  \
  C_ELEMENT(m, n + 0) = vgetq_lane_f32(result##m##n, 0);                       \
  C_ELEMENT(m, n + 1) = vgetq_lane_f32(result##m##n, 1);                       \
  C_ELEMENT(m, n + 2) = vgetq_lane_f32(result##m##n, 2);                       \
  C_ELEMENT(m, n + 3) = vgetq_lane_f32(result##m##n, 3);
#else
#define VECTOR_STORE4(m, n)                                                    \
  result##m##n = vmulq_f32(result##m##n, vdupq_n_f32(alpha));                  \
  result##m##n =                                                               \
    vfmaq_f32(result##m##n, vld1q_f32(&C_ELEMENT(m, n)), vdupq_n_f32(beta));   \
  vst1q_f32(&C_ELEMENT(m, n), result##m##n);
#define STORE(m, n)                                                            \
  C_ELEMENT(m, n) = C_ELEMENT(m, n) * beta + alpha * result##m##n;
#define SCATTER_STORE4(m, n)                                                   \
  result##m##n = vmulq_f32(result##m##n, vdupq_n_f32(alpha));                  \
  C_ELEMENT(m, n + 0) =                                                        \
    C_ELEMENT(m, n + 0) * beta + vgetq_lane_f32(result##m##n, 0);              \
  C_ELEMENT(m, n + 1) =                                                        \
    C_ELEMENT(m, n + 1) * beta + vgetq_lane_f32(result##m##n, 1);              \
  C_ELEMENT(m, n + 2) =                                                        \
    C_ELEMENT(m, n + 2) * beta + vgetq_lane_f32(result##m##n, 2);              \
  C_ELEMENT(m, n + 3) =                                                        \
    C_ELEMENT(m, n + 3) * beta + vgetq_lane_f32(result##m##n, 3);
#endif

#ifdef B0
int
CNAME(BLASLONG M,
      BLASLONG N,
      BLASLONG K,
      IFLOAT* A,
      BLASLONG lda,
      FLOAT alpha,
      IFLOAT* B,
      BLASLONG ldb,
      FLOAT* C,
      BLASLONG ldc)
#else
int
CNAME(BLASLONG M,
      BLASLONG N,
      BLASLONG K,
      IFLOAT* A,
      BLASLONG lda,
      FLOAT alpha,
      IFLOAT* B,
      BLASLONG ldb,
      FLOAT beta,
      FLOAT* C,
      BLASLONG ldc)
#endif
{
  BLASLONG i, j, k;
  BLASLONG n4 = N & ~3;
  BLASLONG n1 = N;
  BLASLONG m8 = M & ~7;
  BLASLONG m4 = M & ~3;
  BLASLONG m1 = M;
  BLASLONG k4 = K & ~3;
  BLASLONG k1 = K;
  int pack_b = M >= 8 && N >= 8 && K >= 8 ? 1 : 0;
  FLOAT* packed_b;
  if (pack_b)
    packed_b = (FLOAT*)malloc(K * 4 * sizeof(FLOAT));

  j = 0;
  for (; j < n4; j += 4) {

    i = 0;
    for (; i < m8; i += 8) {

      k = 0;
      DECLARE_RESULT_VECTOR4(0, 0);
      DECLARE_RESULT_VECTOR4(0, 1);
      DECLARE_RESULT_VECTOR4(0, 2);
      DECLARE_RESULT_VECTOR4(0, 3);
      DECLARE_RESULT_VECTOR4(4, 0);
      DECLARE_RESULT_VECTOR4(4, 1);
      DECLARE_RESULT_VECTOR4(4, 2);
      DECLARE_RESULT_VECTOR4(4, 3);

      if (pack_b) {
        if (i == 0) {
          for (; k < k4; k += 4) {

            VECTOR_LOAD_B_K4(0, 0);
            VECTOR_LOAD_B_K4(1, 0);
            VECTOR_LOAD_B_K4(2, 0);
            VECTOR_LOAD_B_K4(3, 0);
            TRANSPOSE_B4_K4(0, 1, 2, 3, 0, 1, 2, 3);
            VECTOR_PACK_B4(0, 0);
            VECTOR_PACK_B4(0, 1);
            VECTOR_PACK_B4(0, 2);
            VECTOR_PACK_B4(0, 3);
            VECTOR_LOAD_A4(0, 0);
            UPDATE_RESULT_VECTOR4_LANE4(0, 0, 0, 0, 0);
            UPDATE_RESULT_VECTOR4_LANE4(0, 1, 0, 1, 0);
            UPDATE_RESULT_VECTOR4_LANE4(0, 2, 0, 2, 0);
            UPDATE_RESULT_VECTOR4_LANE4(0, 3, 0, 3, 0);
            VECTOR_LOAD_A4(0, 1);
            UPDATE_RESULT_VECTOR4_LANE4(0, 0, 0, 0, 1);
            UPDATE_RESULT_VECTOR4_LANE4(0, 1, 0, 1, 1);
            UPDATE_RESULT_VECTOR4_LANE4(0, 2, 0, 2, 1);
            UPDATE_RESULT_VECTOR4_LANE4(0, 3, 0, 3, 1);
            VECTOR_LOAD_A4(0, 2);
            UPDATE_RESULT_VECTOR4_LANE4(0, 0, 0, 0, 2);
            UPDATE_RESULT_VECTOR4_LANE4(0, 1, 0, 1, 2);
            UPDATE_RESULT_VECTOR4_LANE4(0, 2, 0, 2, 2);
            UPDATE_RESULT_VECTOR4_LANE4(0, 3, 0, 3, 2);
            VECTOR_LOAD_A4(0, 3);
            UPDATE_RESULT_VECTOR4_LANE4(0, 0, 0, 0, 3);
            UPDATE_RESULT_VECTOR4_LANE4(0, 1, 0, 1, 3);
            UPDATE_RESULT_VECTOR4_LANE4(0, 2, 0, 2, 3);
            UPDATE_RESULT_VECTOR4_LANE4(0, 3, 0, 3, 3);
            VECTOR_LOAD_A4(4, 0);
            UPDATE_RESULT_VECTOR4_LANE4(4, 0, 0, 0, 0);
            UPDATE_RESULT_VECTOR4_LANE4(4, 1, 0, 1, 0);
            UPDATE_RESULT_VECTOR4_LANE4(4, 2, 0, 2, 0);
            UPDATE_RESULT_VECTOR4_LANE4(4, 3, 0, 3, 0);
            VECTOR_LOAD_A4(4, 1);
            UPDATE_RESULT_VECTOR4_LANE4(4, 0, 0, 0, 1);
            UPDATE_RESULT_VECTOR4_LANE4(4, 1, 0, 1, 1);
            UPDATE_RESULT_VECTOR4_LANE4(4, 2, 0, 2, 1);
            UPDATE_RESULT_VECTOR4_LANE4(4, 3, 0, 3, 1);
            VECTOR_LOAD_A4(4, 2);
            UPDATE_RESULT_VECTOR4_LANE4(4, 0, 0, 0, 2);
            UPDATE_RESULT_VECTOR4_LANE4(4, 1, 0, 1, 2);
            UPDATE_RESULT_VECTOR4_LANE4(4, 2, 0, 2, 2);
            UPDATE_RESULT_VECTOR4_LANE4(4, 3, 0, 3, 2);
            VECTOR_LOAD_A4(4, 3);
            UPDATE_RESULT_VECTOR4_LANE4(4, 0, 0, 0, 3);
            UPDATE_RESULT_VECTOR4_LANE4(4, 1, 0, 1, 3);
            UPDATE_RESULT_VECTOR4_LANE4(4, 2, 0, 2, 3);
            UPDATE_RESULT_VECTOR4_LANE4(4, 3, 0, 3, 3);
          }
          for (; k < k1; k++) {

            BROADCAST_LOAD_B4(0, 0);
            PACK_B4(0, 0);
            VECTOR_LOAD_A4(0, 0);
            UPDATE_RESULT_VECTOR4(0, 0, 0);
            BROADCAST_LOAD_B4(1, 0);
            PACK_B4(1, 0);
            UPDATE_RESULT_VECTOR4(0, 1, 0);
            VECTOR_LOAD_A4(4, 0);
            UPDATE_RESULT_VECTOR4(4, 0, 0);
            UPDATE_RESULT_VECTOR4(4, 1, 0);
            BROADCAST_LOAD_B4(2, 0);
            PACK_B4(2, 0);
            UPDATE_RESULT_VECTOR4(0, 2, 0);
            UPDATE_RESULT_VECTOR4(4, 2, 0);
            BROADCAST_LOAD_B4(3, 0);
            PACK_B4(3, 0);
            UPDATE_RESULT_VECTOR4(0, 3, 0);
            UPDATE_RESULT_VECTOR4(4, 3, 0);
          }
        } else {
          for (; k < k1; k++) {

            VECTOR_UNPACK_B4(0, 0);
            VECTOR_LOAD_A4(0, 0);
            UPDATE_RESULT_VECTOR4_LANE4(0, 0, 0, 0, 0);
            UPDATE_RESULT_VECTOR4_LANE4(0, 1, 0, 1, 0);
            UPDATE_RESULT_VECTOR4_LANE4(0, 2, 0, 2, 0);
            UPDATE_RESULT_VECTOR4_LANE4(0, 3, 0, 3, 0);
            VECTOR_LOAD_A4(4, 0);
            UPDATE_RESULT_VECTOR4_LANE4(4, 0, 0, 0, 0);
            UPDATE_RESULT_VECTOR4_LANE4(4, 1, 0, 1, 0);
            UPDATE_RESULT_VECTOR4_LANE4(4, 2, 0, 2, 0);
            UPDATE_RESULT_VECTOR4_LANE4(4, 3, 0, 3, 0);
          }
        }
      } else {
        for (; k < k4; k += 4) {

          VECTOR_LOAD_B_K4(0, 0);
          VECTOR_LOAD_B_K4(1, 0);
          VECTOR_LOAD_B_K4(2, 0);
          VECTOR_LOAD_B_K4(3, 0);
          TRANSPOSE_B4_K4(0, 1, 2, 3, 0, 1, 2, 3);
          VECTOR_LOAD_A4(0, 0);
          UPDATE_RESULT_VECTOR4_LANE4(0, 0, 0, 0, 0);
          UPDATE_RESULT_VECTOR4_LANE4(0, 1, 0, 1, 0);
          UPDATE_RESULT_VECTOR4_LANE4(0, 2, 0, 2, 0);
          UPDATE_RESULT_VECTOR4_LANE4(0, 3, 0, 3, 0);
          VECTOR_LOAD_A4(0, 1);
          UPDATE_RESULT_VECTOR4_LANE4(0, 0, 0, 0, 1);
          UPDATE_RESULT_VECTOR4_LANE4(0, 1, 0, 1, 1);
          UPDATE_RESULT_VECTOR4_LANE4(0, 2, 0, 2, 1);
          UPDATE_RESULT_VECTOR4_LANE4(0, 3, 0, 3, 1);
          VECTOR_LOAD_A4(0, 2);
          UPDATE_RESULT_VECTOR4_LANE4(0, 0, 0, 0, 2);
          UPDATE_RESULT_VECTOR4_LANE4(0, 1, 0, 1, 2);
          UPDATE_RESULT_VECTOR4_LANE4(0, 2, 0, 2, 2);
          UPDATE_RESULT_VECTOR4_LANE4(0, 3, 0, 3, 2);
          VECTOR_LOAD_A4(0, 3);
          UPDATE_RESULT_VECTOR4_LANE4(0, 0, 0, 0, 3);
          UPDATE_RESULT_VECTOR4_LANE4(0, 1, 0, 1, 3);
          UPDATE_RESULT_VECTOR4_LANE4(0, 2, 0, 2, 3);
          UPDATE_RESULT_VECTOR4_LANE4(0, 3, 0, 3, 3);
          VECTOR_LOAD_A4(4, 0);
          UPDATE_RESULT_VECTOR4_LANE4(4, 0, 0, 0, 0);
          UPDATE_RESULT_VECTOR4_LANE4(4, 1, 0, 1, 0);
          UPDATE_RESULT_VECTOR4_LANE4(4, 2, 0, 2, 0);
          UPDATE_RESULT_VECTOR4_LANE4(4, 3, 0, 3, 0);
          VECTOR_LOAD_A4(4, 1);
          UPDATE_RESULT_VECTOR4_LANE4(4, 0, 0, 0, 1);
          UPDATE_RESULT_VECTOR4_LANE4(4, 1, 0, 1, 1);
          UPDATE_RESULT_VECTOR4_LANE4(4, 2, 0, 2, 1);
          UPDATE_RESULT_VECTOR4_LANE4(4, 3, 0, 3, 1);
          VECTOR_LOAD_A4(4, 2);
          UPDATE_RESULT_VECTOR4_LANE4(4, 0, 0, 0, 2);
          UPDATE_RESULT_VECTOR4_LANE4(4, 1, 0, 1, 2);
          UPDATE_RESULT_VECTOR4_LANE4(4, 2, 0, 2, 2);
          UPDATE_RESULT_VECTOR4_LANE4(4, 3, 0, 3, 2);
          VECTOR_LOAD_A4(4, 3);
          UPDATE_RESULT_VECTOR4_LANE4(4, 0, 0, 0, 3);
          UPDATE_RESULT_VECTOR4_LANE4(4, 1, 0, 1, 3);
          UPDATE_RESULT_VECTOR4_LANE4(4, 2, 0, 2, 3);
          UPDATE_RESULT_VECTOR4_LANE4(4, 3, 0, 3, 3);
        }
        for (; k < k1; k++) {

          BROADCAST_LOAD_B4(0, 0);
          VECTOR_LOAD_A4(0, 0);
          UPDATE_RESULT_VECTOR4(0, 0, 0);
          BROADCAST_LOAD_B4(1, 0);
          UPDATE_RESULT_VECTOR4(0, 1, 0);
          VECTOR_LOAD_A4(4, 0);
          UPDATE_RESULT_VECTOR4(4, 0, 0);
          UPDATE_RESULT_VECTOR4(4, 1, 0);
          BROADCAST_LOAD_B4(2, 0);
          UPDATE_RESULT_VECTOR4(0, 2, 0);
          UPDATE_RESULT_VECTOR4(4, 2, 0);
          BROADCAST_LOAD_B4(3, 0);
          UPDATE_RESULT_VECTOR4(0, 3, 0);
          UPDATE_RESULT_VECTOR4(4, 3, 0);
        }
      }
      VECTOR_STORE4(0, 0);
      VECTOR_STORE4(0, 1);
      VECTOR_STORE4(0, 2);
      VECTOR_STORE4(0, 3);
      VECTOR_STORE4(4, 0);
      VECTOR_STORE4(4, 1);
      VECTOR_STORE4(4, 2);
      VECTOR_STORE4(4, 3);
    }
    for (; i < m4; i += 4) {

      k = 0;
      DECLARE_RESULT_VECTOR4(0, 0);
      DECLARE_RESULT_VECTOR4(0, 1);
      DECLARE_RESULT_VECTOR4(0, 2);
      DECLARE_RESULT_VECTOR4(0, 3);

      if (pack_b) {
        for (; k < k1; k++) {

          VECTOR_UNPACK_B4(0, 0);
          VECTOR_LOAD_A4(0, 0);
          UPDATE_RESULT_VECTOR4_LANE4(0, 0, 0, 0, 0);
          UPDATE_RESULT_VECTOR4_LANE4(0, 1, 0, 1, 0);
          UPDATE_RESULT_VECTOR4_LANE4(0, 2, 0, 2, 0);
          UPDATE_RESULT_VECTOR4_LANE4(0, 3, 0, 3, 0);
        }
      } else {
        for (; k < k4; k += 4) {

          VECTOR_LOAD_B_K4(0, 0);
          VECTOR_LOAD_B_K4(1, 0);
          VECTOR_LOAD_B_K4(2, 0);
          VECTOR_LOAD_B_K4(3, 0);
          TRANSPOSE_B4_K4(0, 1, 2, 3, 0, 1, 2, 3);
          VECTOR_LOAD_A4(0, 0);
          UPDATE_RESULT_VECTOR4_LANE4(0, 0, 0, 0, 0);
          UPDATE_RESULT_VECTOR4_LANE4(0, 1, 0, 1, 0);
          UPDATE_RESULT_VECTOR4_LANE4(0, 2, 0, 2, 0);
          UPDATE_RESULT_VECTOR4_LANE4(0, 3, 0, 3, 0);
          VECTOR_LOAD_A4(0, 1);
          UPDATE_RESULT_VECTOR4_LANE4(0, 0, 0, 0, 1);
          UPDATE_RESULT_VECTOR4_LANE4(0, 1, 0, 1, 1);
          UPDATE_RESULT_VECTOR4_LANE4(0, 2, 0, 2, 1);
          UPDATE_RESULT_VECTOR4_LANE4(0, 3, 0, 3, 1);
          VECTOR_LOAD_A4(0, 2);
          UPDATE_RESULT_VECTOR4_LANE4(0, 0, 0, 0, 2);
          UPDATE_RESULT_VECTOR4_LANE4(0, 1, 0, 1, 2);
          UPDATE_RESULT_VECTOR4_LANE4(0, 2, 0, 2, 2);
          UPDATE_RESULT_VECTOR4_LANE4(0, 3, 0, 3, 2);
          VECTOR_LOAD_A4(0, 3);
          UPDATE_RESULT_VECTOR4_LANE4(0, 0, 0, 0, 3);
          UPDATE_RESULT_VECTOR4_LANE4(0, 1, 0, 1, 3);
          UPDATE_RESULT_VECTOR4_LANE4(0, 2, 0, 2, 3);
          UPDATE_RESULT_VECTOR4_LANE4(0, 3, 0, 3, 3);
        }
        for (; k < k1; k++) {

          BROADCAST_LOAD_B4(0, 0);
          VECTOR_LOAD_A4(0, 0);
          UPDATE_RESULT_VECTOR4(0, 0, 0);
          BROADCAST_LOAD_B4(1, 0);
          UPDATE_RESULT_VECTOR4(0, 1, 0);
          BROADCAST_LOAD_B4(2, 0);
          UPDATE_RESULT_VECTOR4(0, 2, 0);
          BROADCAST_LOAD_B4(3, 0);
          UPDATE_RESULT_VECTOR4(0, 3, 0);
        }
      }
      VECTOR_STORE4(0, 0);
      VECTOR_STORE4(0, 1);
      VECTOR_STORE4(0, 2);
      VECTOR_STORE4(0, 3);
    }
    for (; i < m1; i++) {

      k = 0;
      DECLARE_RESULT_VECTOR4(0, 0);

      if (pack_b) {
        for (; k < k1; k++) {

          VECTOR_UNPACK_B4(0, 0);
          BROADCAST_LOAD_A4(0, 0);
          UPDATE_RESULT_VECTOR4(0, 0, 0);
        }
      } else {
        for (; k < k1; k++) {

          GATHER_LOAD_B4(0, 0);
          BROADCAST_LOAD_A4(0, 0);
          UPDATE_RESULT_VECTOR4(0, 0, 0);
        }
      }
      SCATTER_STORE4(0, 0);
    }
  }
  for (; j < n1; j++) {

    i = 0;
    for (; i < m8; i += 8) {

      k = 0;
      DECLARE_RESULT_VECTOR4(0, 0);
      DECLARE_RESULT_VECTOR4(4, 0);

      for (; k < k1; k++) {

        BROADCAST_LOAD_B4(0, 0);
        VECTOR_LOAD_A4(0, 0);
        UPDATE_RESULT_VECTOR4(0, 0, 0);
        VECTOR_LOAD_A4(4, 0);
        UPDATE_RESULT_VECTOR4(4, 0, 0);
      }
      VECTOR_STORE4(0, 0);
      VECTOR_STORE4(4, 0);
    }
    for (; i < m4; i += 4) {

      k = 0;
      DECLARE_RESULT_VECTOR4(0, 0);

      for (; k < k1; k++) {

        BROADCAST_LOAD_B4(0, 0);
        VECTOR_LOAD_A4(0, 0);
        UPDATE_RESULT_VECTOR4(0, 0, 0);
      }
      VECTOR_STORE4(0, 0);
    }
    for (; i < m1; i++) {

      k = 0;
      DECLARE_RESULT(0, 0);

      for (k = 0; k < K; k++) {
        LOAD_B1(0, 0);
        LOAD_A1(0, 0);
        UPDATE_RESULT(0, 0, 0);
      }
      STORE(0, 0);
    }
  }

  if (pack_b)
    free(packed_b);

  return 0;
}
