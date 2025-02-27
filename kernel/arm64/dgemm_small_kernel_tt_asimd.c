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

#define A_ELEMENT_K(m, offset_k) A[(i + (m)) * lda + (k + offset_k)]
#define A_ELEMENT(m) A_ELEMENT_K(m, 0)

#define B_ELEMENT_K(n, offset_k) B[(k + offset_k) * ldb + (j + (n))]
#define B_ELEMENT(n) B_ELEMENT_K(n, 0)

#define C_ELEMENT(m, n) C[(i + (m)) + (j + (n)) * ldc]
#define PACK_ELEMENT_K(m, offset_k) packed_a[(k + offset_k) * 4 + m]
#define PACK_ELEMENT(m) PACK_ELEMENT_K(m, 0)

// ASIMD
#define DECLARE_RESULT_VECTOR2(n, m)                                           \
  float64x2_t result##n##m = vdupq_n_f64(0.0);
#define DECLARE_RESULT(n, m) float64_t result##n##m = 0.0;
#define BROADCAST_LOAD_B2(n, offset_k)                                         \
  float64x2_t b##n##_k##offset_k = vld1q_dup_f64(&B_ELEMENT_K(n, offset_k));
#define LOAD_B1(n, offset_k)                                                   \
  float64_t b##n##_k##offset_k = B_ELEMENT_K(n, offset_k);
#define BROADCAST_LOAD_A2(m, offset_k)                                         \
  float64x2_t a##m##_k##offset_k = vld1q_dup_f64(&A_ELEMENT_K(m, offset_k));
#define LOAD_A1(m, offset_k)                                                   \
  float64_t a##m##_k##offset_k = A_ELEMENT_K(m, offset_k);
#define VECTOR_LOAD_B2(n, offset_k)                                            \
  float64x2_t b##n##_k##offset_k = vld1q_f64(&B_ELEMENT_K(n, offset_k));
#define VECTOR_LOAD_A_K2(m, offset_k)                                          \
  float64x2_t a##k##m##_k##offset_k = vld1q_f64(&A_ELEMENT_K(m, offset_k));
#define TRANSPOSE_A2_K2(m0, m1, offset_k0, offset_k1)                          \
  float64x2_t a##m0##_k##offset_k0 =                                           \
    vzip1q_f64(a##k##m0##_k##offset_k0, a##k##m1##_k##offset_k0);              \
  float64x2_t a##m0##_k##offset_k1 =                                           \
    vzip2q_f64(a##k##m0##_k##offset_k0, a##k##m1##_k##offset_k0);

#define GATHER_LOAD_A2(m, offset_k)                                            \
  float64x2_t a##m##_k##offset_k = vdupq_n_f64(A_ELEMENT_K(m, offset_k));      \
  a##m##_k##offset_k =                                                         \
    vsetq_lane_f64(A_ELEMENT_K(m + 1, offset_k), a##m##_k##offset_k, 1);
#define VECTOR_UNPACK_A2(m, offset_k)                                          \
  float64x2_t a##m##_k##offset_k = vld1q_f64(&PACK_ELEMENT_K(m, offset_k));
#define VECTOR_PACK_A2(m, offset_k)                                            \
  vst1q_f64(&PACK_ELEMENT_K(m, offset_k), a##m##_k##offset_k);
#define PACK_A2(m, offset_k)                                                   \
  PACK_ELEMENT_K(m, offset_k) = vgetq_lane_f64(a##m##_k##offset_k, 0);
#define BROADCAST_UNPACK_A2(m, offset_k)                                       \
  float64x2_t a##m##_k##offset_k = vdupq_n_f64(PACK_ELEMENT_K(m, offset_k));
#define UPDATE_RESULT_VECTOR2(n, m, offset_k)                                  \
  result##n##m =                                                               \
    vfmaq_f64(result##n##m, b##n##_k##offset_k, a##m##_k##offset_k);
#define UPDATE_RESULT(n, m, offset_k)                                          \
  result##n##m = result##n##m + b##n##_k##offset_k * a##m##_k##offset_k;
#define UPDATE_RESULT_VECTOR2_LANE2(n, m, outer, lane, offset_k)               \
  result##n##m = vfmaq_laneq_f64(                                              \
    result##n##m, b##n##_k##offset_k, a##outer##_k##offset_k, lane);
#ifdef B0
#define VECTOR_STORE2(n, m)                                                    \
  vst1q_f64(&C_ELEMENT(m, n), vmulq_f64(result##n##m, vdupq_n_f64(alpha)));
#define STORE(n, m) C_ELEMENT(m, n) = alpha * result##n##m;
#define SCATTER_STORE2(n, m)                                                   \
  result##n##m = vmulq_f64(result##n##m, vdupq_n_f64(alpha));                  \
  C_ELEMENT(m, n + 0) = vgetq_lane_f64(result##n##m, 0);                       \
  C_ELEMENT(m, n + 1) = vgetq_lane_f64(result##n##m, 1);
#else
#define VECTOR_STORE2(n, m)                                                    \
  result##n##m = vmulq_f64(result##n##m, vdupq_n_f64(alpha));                  \
  result##n##m =                                                               \
    vfmaq_f64(result##n##m, vld1q_f64(&C_ELEMENT(m, n)), vdupq_n_f64(beta));   \
  vst1q_f64(&C_ELEMENT(m, n), result##n##m);
#define STORE(n, m)                                                            \
  C_ELEMENT(m, n) = C_ELEMENT(m, n) * beta + alpha * result##n##m;
#define SCATTER_STORE2(n, m)                                                   \
  result##n##m = vmulq_f64(result##n##m, vdupq_n_f64(alpha));                  \
  C_ELEMENT(m, n + 0) =                                                        \
    C_ELEMENT(m, n + 0) * beta + vgetq_lane_f64(result##n##m, 0);              \
  C_ELEMENT(m, n + 1) =                                                        \
    C_ELEMENT(m, n + 1) * beta + vgetq_lane_f64(result##n##m, 1);
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
  BLASLONG m4 = M & ~3;
  BLASLONG m2 = M & ~1;
  BLASLONG m1 = M;
  BLASLONG n4 = N & ~3;
  BLASLONG n2 = N & ~1;
  BLASLONG n1 = N;
  BLASLONG k2 = K & ~1;
  BLASLONG k1 = K;
  int pack_a = M >= 8 && N >= 8 && K >= 8 ? 1 : 0;
  FLOAT* packed_a;
  if (pack_a)
    packed_a = (FLOAT*)malloc(K * 4 * sizeof(FLOAT));

  i = 0;
  for (; i < m4; i += 4) {

    j = 0;
    for (; j < n4; j += 4) {

      k = 0;
      DECLARE_RESULT_VECTOR2(0, 0);
      DECLARE_RESULT_VECTOR2(0, 1);
      DECLARE_RESULT_VECTOR2(0, 2);
      DECLARE_RESULT_VECTOR2(0, 3);
      DECLARE_RESULT_VECTOR2(2, 0);
      DECLARE_RESULT_VECTOR2(2, 1);
      DECLARE_RESULT_VECTOR2(2, 2);
      DECLARE_RESULT_VECTOR2(2, 3);

      if (pack_a) {
        if (j == 0) {
          for (; k < k2; k += 2) {

            VECTOR_LOAD_A_K2(0, 0);
            VECTOR_LOAD_A_K2(1, 0);
            TRANSPOSE_A2_K2(0, 1, 0, 1);
            VECTOR_PACK_A2(0, 0);
            VECTOR_PACK_A2(0, 1);
            VECTOR_LOAD_B2(0, 0);
            UPDATE_RESULT_VECTOR2_LANE2(0, 0, 0, 0, 0);
            UPDATE_RESULT_VECTOR2_LANE2(0, 1, 0, 1, 0);
            VECTOR_LOAD_B2(0, 1);
            UPDATE_RESULT_VECTOR2_LANE2(0, 0, 0, 0, 1);
            UPDATE_RESULT_VECTOR2_LANE2(0, 1, 0, 1, 1);
            VECTOR_LOAD_A_K2(2, 0);
            VECTOR_LOAD_A_K2(3, 0);
            TRANSPOSE_A2_K2(2, 3, 0, 1);
            VECTOR_PACK_A2(2, 0);
            VECTOR_PACK_A2(2, 1);
            UPDATE_RESULT_VECTOR2_LANE2(0, 2, 2, 0, 0);
            UPDATE_RESULT_VECTOR2_LANE2(0, 3, 2, 1, 0);
            UPDATE_RESULT_VECTOR2_LANE2(0, 2, 2, 0, 1);
            UPDATE_RESULT_VECTOR2_LANE2(0, 3, 2, 1, 1);
            VECTOR_LOAD_B2(2, 0);
            UPDATE_RESULT_VECTOR2_LANE2(2, 0, 0, 0, 0);
            UPDATE_RESULT_VECTOR2_LANE2(2, 1, 0, 1, 0);
            UPDATE_RESULT_VECTOR2_LANE2(2, 2, 2, 0, 0);
            UPDATE_RESULT_VECTOR2_LANE2(2, 3, 2, 1, 0);
            VECTOR_LOAD_B2(2, 1);
            UPDATE_RESULT_VECTOR2_LANE2(2, 0, 0, 0, 1);
            UPDATE_RESULT_VECTOR2_LANE2(2, 1, 0, 1, 1);
            UPDATE_RESULT_VECTOR2_LANE2(2, 2, 2, 0, 1);
            UPDATE_RESULT_VECTOR2_LANE2(2, 3, 2, 1, 1);
          }
          for (; k < k1; k++) {

            BROADCAST_LOAD_A2(0, 0);
            PACK_A2(0, 0);
            VECTOR_LOAD_B2(0, 0);
            UPDATE_RESULT_VECTOR2(0, 0, 0);
            BROADCAST_LOAD_A2(1, 0);
            PACK_A2(1, 0);
            UPDATE_RESULT_VECTOR2(0, 1, 0);
            VECTOR_LOAD_B2(2, 0);
            UPDATE_RESULT_VECTOR2(2, 0, 0);
            UPDATE_RESULT_VECTOR2(2, 1, 0);
            BROADCAST_LOAD_A2(2, 0);
            PACK_A2(2, 0);
            UPDATE_RESULT_VECTOR2(0, 2, 0);
            UPDATE_RESULT_VECTOR2(2, 2, 0);
            BROADCAST_LOAD_A2(3, 0);
            PACK_A2(3, 0);
            UPDATE_RESULT_VECTOR2(0, 3, 0);
            UPDATE_RESULT_VECTOR2(2, 3, 0);
          }
        } else {
          for (; k < k1; k++) {

            VECTOR_UNPACK_A2(0, 0);
            VECTOR_LOAD_B2(0, 0);
            UPDATE_RESULT_VECTOR2_LANE2(0, 0, 0, 0, 0);
            UPDATE_RESULT_VECTOR2_LANE2(0, 1, 0, 1, 0);
            VECTOR_UNPACK_A2(2, 0);
            UPDATE_RESULT_VECTOR2_LANE2(0, 2, 2, 0, 0);
            UPDATE_RESULT_VECTOR2_LANE2(0, 3, 2, 1, 0);
            VECTOR_LOAD_B2(2, 0);
            UPDATE_RESULT_VECTOR2_LANE2(2, 0, 0, 0, 0);
            UPDATE_RESULT_VECTOR2_LANE2(2, 1, 0, 1, 0);
            UPDATE_RESULT_VECTOR2_LANE2(2, 2, 2, 0, 0);
            UPDATE_RESULT_VECTOR2_LANE2(2, 3, 2, 1, 0);
          }
        }
      } else {
        for (; k < k2; k += 2) {

          VECTOR_LOAD_A_K2(0, 0);
          VECTOR_LOAD_A_K2(1, 0);
          TRANSPOSE_A2_K2(0, 1, 0, 1);
          VECTOR_LOAD_B2(0, 0);
          UPDATE_RESULT_VECTOR2_LANE2(0, 0, 0, 0, 0);
          UPDATE_RESULT_VECTOR2_LANE2(0, 1, 0, 1, 0);
          VECTOR_LOAD_B2(0, 1);
          UPDATE_RESULT_VECTOR2_LANE2(0, 0, 0, 0, 1);
          UPDATE_RESULT_VECTOR2_LANE2(0, 1, 0, 1, 1);
          VECTOR_LOAD_A_K2(2, 0);
          VECTOR_LOAD_A_K2(3, 0);
          TRANSPOSE_A2_K2(2, 3, 0, 1);
          UPDATE_RESULT_VECTOR2_LANE2(0, 2, 2, 0, 0);
          UPDATE_RESULT_VECTOR2_LANE2(0, 3, 2, 1, 0);
          UPDATE_RESULT_VECTOR2_LANE2(0, 2, 2, 0, 1);
          UPDATE_RESULT_VECTOR2_LANE2(0, 3, 2, 1, 1);
          VECTOR_LOAD_B2(2, 0);
          UPDATE_RESULT_VECTOR2_LANE2(2, 0, 0, 0, 0);
          UPDATE_RESULT_VECTOR2_LANE2(2, 1, 0, 1, 0);
          UPDATE_RESULT_VECTOR2_LANE2(2, 2, 2, 0, 0);
          UPDATE_RESULT_VECTOR2_LANE2(2, 3, 2, 1, 0);
          VECTOR_LOAD_B2(2, 1);
          UPDATE_RESULT_VECTOR2_LANE2(2, 0, 0, 0, 1);
          UPDATE_RESULT_VECTOR2_LANE2(2, 1, 0, 1, 1);
          UPDATE_RESULT_VECTOR2_LANE2(2, 2, 2, 0, 1);
          UPDATE_RESULT_VECTOR2_LANE2(2, 3, 2, 1, 1);
        }
        for (; k < k1; k++) {

          BROADCAST_LOAD_A2(0, 0);
          VECTOR_LOAD_B2(0, 0);
          UPDATE_RESULT_VECTOR2(0, 0, 0);
          BROADCAST_LOAD_A2(1, 0);
          UPDATE_RESULT_VECTOR2(0, 1, 0);
          VECTOR_LOAD_B2(2, 0);
          UPDATE_RESULT_VECTOR2(2, 0, 0);
          UPDATE_RESULT_VECTOR2(2, 1, 0);
          BROADCAST_LOAD_A2(2, 0);
          UPDATE_RESULT_VECTOR2(0, 2, 0);
          UPDATE_RESULT_VECTOR2(2, 2, 0);
          BROADCAST_LOAD_A2(3, 0);
          UPDATE_RESULT_VECTOR2(0, 3, 0);
          UPDATE_RESULT_VECTOR2(2, 3, 0);
        }
      }
      SCATTER_STORE2(0, 0);
      SCATTER_STORE2(0, 1);
      SCATTER_STORE2(0, 2);
      SCATTER_STORE2(0, 3);
      SCATTER_STORE2(2, 0);
      SCATTER_STORE2(2, 1);
      SCATTER_STORE2(2, 2);
      SCATTER_STORE2(2, 3);
    }
    for (; j < n2; j += 2) {

      k = 0;
      DECLARE_RESULT_VECTOR2(0, 0);
      DECLARE_RESULT_VECTOR2(0, 1);
      DECLARE_RESULT_VECTOR2(0, 2);
      DECLARE_RESULT_VECTOR2(0, 3);

      if (pack_a) {
        for (; k < k1; k++) {

          VECTOR_UNPACK_A2(0, 0);
          VECTOR_LOAD_B2(0, 0);
          UPDATE_RESULT_VECTOR2_LANE2(0, 0, 0, 0, 0);
          UPDATE_RESULT_VECTOR2_LANE2(0, 1, 0, 1, 0);
          VECTOR_UNPACK_A2(2, 0);
          UPDATE_RESULT_VECTOR2_LANE2(0, 2, 2, 0, 0);
          UPDATE_RESULT_VECTOR2_LANE2(0, 3, 2, 1, 0);
        }
      } else {
        for (; k < k2; k += 2) {

          VECTOR_LOAD_A_K2(0, 0);
          VECTOR_LOAD_A_K2(1, 0);
          TRANSPOSE_A2_K2(0, 1, 0, 1);
          VECTOR_LOAD_B2(0, 0);
          UPDATE_RESULT_VECTOR2_LANE2(0, 0, 0, 0, 0);
          UPDATE_RESULT_VECTOR2_LANE2(0, 1, 0, 1, 0);
          VECTOR_LOAD_B2(0, 1);
          UPDATE_RESULT_VECTOR2_LANE2(0, 0, 0, 0, 1);
          UPDATE_RESULT_VECTOR2_LANE2(0, 1, 0, 1, 1);
          VECTOR_LOAD_A_K2(2, 0);
          VECTOR_LOAD_A_K2(3, 0);
          TRANSPOSE_A2_K2(2, 3, 0, 1);
          UPDATE_RESULT_VECTOR2_LANE2(0, 2, 2, 0, 0);
          UPDATE_RESULT_VECTOR2_LANE2(0, 3, 2, 1, 0);
          UPDATE_RESULT_VECTOR2_LANE2(0, 2, 2, 0, 1);
          UPDATE_RESULT_VECTOR2_LANE2(0, 3, 2, 1, 1);
        }
        for (; k < k1; k++) {

          BROADCAST_LOAD_A2(0, 0);
          VECTOR_LOAD_B2(0, 0);
          UPDATE_RESULT_VECTOR2(0, 0, 0);
          BROADCAST_LOAD_A2(1, 0);
          UPDATE_RESULT_VECTOR2(0, 1, 0);
          BROADCAST_LOAD_A2(2, 0);
          UPDATE_RESULT_VECTOR2(0, 2, 0);
          BROADCAST_LOAD_A2(3, 0);
          UPDATE_RESULT_VECTOR2(0, 3, 0);
        }
      }
      SCATTER_STORE2(0, 0);
      SCATTER_STORE2(0, 1);
      SCATTER_STORE2(0, 2);
      SCATTER_STORE2(0, 3);
    }
    for (; j < n1; j++) {

      k = 0;
      DECLARE_RESULT_VECTOR2(0, 0);
      DECLARE_RESULT_VECTOR2(0, 2);

      if (pack_a) {
        for (; k < k1; k++) {

          VECTOR_UNPACK_A2(0, 0);
          BROADCAST_LOAD_B2(0, 0);
          UPDATE_RESULT_VECTOR2(0, 0, 0);
          VECTOR_UNPACK_A2(2, 0);
          UPDATE_RESULT_VECTOR2(0, 2, 0);
        }
      } else {
        for (; k < k1; k++) {

          GATHER_LOAD_A2(0, 0);
          BROADCAST_LOAD_B2(0, 0);
          UPDATE_RESULT_VECTOR2(0, 0, 0);
          GATHER_LOAD_A2(2, 0);
          UPDATE_RESULT_VECTOR2(0, 2, 0);
        }
      }
      VECTOR_STORE2(0, 0);
      VECTOR_STORE2(0, 2);
    }
  }
  for (; i < m2; i += 2) {

    j = 0;
    for (; j < n4; j += 4) {

      k = 0;
      DECLARE_RESULT_VECTOR2(0, 0);
      DECLARE_RESULT_VECTOR2(0, 1);
      DECLARE_RESULT_VECTOR2(2, 0);
      DECLARE_RESULT_VECTOR2(2, 1);

      for (; k < k2; k += 2) {

        VECTOR_LOAD_A_K2(0, 0);
        VECTOR_LOAD_A_K2(1, 0);
        TRANSPOSE_A2_K2(0, 1, 0, 1);
        VECTOR_LOAD_B2(0, 0);
        UPDATE_RESULT_VECTOR2_LANE2(0, 0, 0, 0, 0);
        UPDATE_RESULT_VECTOR2_LANE2(0, 1, 0, 1, 0);
        VECTOR_LOAD_B2(0, 1);
        UPDATE_RESULT_VECTOR2_LANE2(0, 0, 0, 0, 1);
        UPDATE_RESULT_VECTOR2_LANE2(0, 1, 0, 1, 1);
        VECTOR_LOAD_B2(2, 0);
        UPDATE_RESULT_VECTOR2_LANE2(2, 0, 0, 0, 0);
        UPDATE_RESULT_VECTOR2_LANE2(2, 1, 0, 1, 0);
        VECTOR_LOAD_B2(2, 1);
        UPDATE_RESULT_VECTOR2_LANE2(2, 0, 0, 0, 1);
        UPDATE_RESULT_VECTOR2_LANE2(2, 1, 0, 1, 1);
      }
      for (; k < k1; k++) {

        BROADCAST_LOAD_A2(0, 0);
        VECTOR_LOAD_B2(0, 0);
        UPDATE_RESULT_VECTOR2(0, 0, 0);
        BROADCAST_LOAD_A2(1, 0);
        UPDATE_RESULT_VECTOR2(0, 1, 0);
        VECTOR_LOAD_B2(2, 0);
        UPDATE_RESULT_VECTOR2(2, 0, 0);
        UPDATE_RESULT_VECTOR2(2, 1, 0);
      }
      SCATTER_STORE2(0, 0);
      SCATTER_STORE2(0, 1);
      SCATTER_STORE2(2, 0);
      SCATTER_STORE2(2, 1);
    }
    for (; j < n2; j += 2) {

      k = 0;
      DECLARE_RESULT_VECTOR2(0, 0);
      DECLARE_RESULT_VECTOR2(0, 1);

      for (; k < k2; k += 2) {

        VECTOR_LOAD_A_K2(0, 0);
        VECTOR_LOAD_A_K2(1, 0);
        TRANSPOSE_A2_K2(0, 1, 0, 1);
        VECTOR_LOAD_B2(0, 0);
        UPDATE_RESULT_VECTOR2_LANE2(0, 0, 0, 0, 0);
        UPDATE_RESULT_VECTOR2_LANE2(0, 1, 0, 1, 0);
        VECTOR_LOAD_B2(0, 1);
        UPDATE_RESULT_VECTOR2_LANE2(0, 0, 0, 0, 1);
        UPDATE_RESULT_VECTOR2_LANE2(0, 1, 0, 1, 1);
      }
      for (; k < k1; k++) {

        BROADCAST_LOAD_A2(0, 0);
        VECTOR_LOAD_B2(0, 0);
        UPDATE_RESULT_VECTOR2(0, 0, 0);
        BROADCAST_LOAD_A2(1, 0);
        UPDATE_RESULT_VECTOR2(0, 1, 0);
      }
      SCATTER_STORE2(0, 0);
      SCATTER_STORE2(0, 1);
    }
    for (; j < n1; j++) {

      k = 0;
      DECLARE_RESULT_VECTOR2(0, 0);

      for (; k < k1; k++) {

        GATHER_LOAD_A2(0, 0);
        BROADCAST_LOAD_B2(0, 0);
        UPDATE_RESULT_VECTOR2(0, 0, 0);
      }
      VECTOR_STORE2(0, 0);
    }
  }
  for (; i < m1; i++) {

    j = 0;
    for (; j < n4; j += 4) {

      k = 0;
      DECLARE_RESULT_VECTOR2(0, 0);
      DECLARE_RESULT_VECTOR2(2, 0);

      for (; k < k1; k++) {

        BROADCAST_LOAD_A2(0, 0);
        VECTOR_LOAD_B2(0, 0);
        UPDATE_RESULT_VECTOR2(0, 0, 0);
        VECTOR_LOAD_B2(2, 0);
        UPDATE_RESULT_VECTOR2(2, 0, 0);
      }
      SCATTER_STORE2(0, 0);
      SCATTER_STORE2(2, 0);
    }
    for (; j < n2; j += 2) {

      k = 0;
      DECLARE_RESULT_VECTOR2(0, 0);

      for (; k < k1; k++) {

        BROADCAST_LOAD_A2(0, 0);
        VECTOR_LOAD_B2(0, 0);
        UPDATE_RESULT_VECTOR2(0, 0, 0);
      }
      SCATTER_STORE2(0, 0);
    }
    for (; j < n1; j++) {

      k = 0;
      DECLARE_RESULT(0, 0);

      for (k = 0; k < K; k++) {
        LOAD_A1(0, 0);
        LOAD_B1(0, 0);
        UPDATE_RESULT(0, 0, 0);
      }
      STORE(0, 0);
    }
  }

  if (pack_a)
    free(packed_a);

  return 0;
}
