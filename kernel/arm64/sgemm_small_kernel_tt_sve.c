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
#include <arm_sve.h>
#ifdef __ARM_NEON_SVE_BRIDGE
#include <arm_neon_sve_bridge.h>
#else
#define svdup_neonq_f32(fixed_reg)                                             \
  ({                                                                           \
    svfloat32_t scalable_reg;                                                  \
    asm("mov %0.q, %q1" : "=w"(scalable_reg) : "w"(fixed_reg) :);              \
    scalable_reg;                                                              \
  })
#define svdup_neonq_f64(fixed_reg)                                             \
  ({                                                                           \
    svfloat64_t scalable_reg;                                                  \
    asm("mov %0.q, %q1" : "=w"(scalable_reg) : "w"(fixed_reg) :);              \
    scalable_reg;                                                              \
  })
#endif

#define A_ELEMENT_K(m, offset_k) A[(i + (m)) * lda + (k + offset_k)]
#define A_ELEMENT(m) A_ELEMENT_K(m, 0)

#define B_ELEMENT_K(n, offset_k) B[(k + offset_k) * ldb + (j + (n))]
#define B_ELEMENT(n) B_ELEMENT_K(n, 0)

#define C_ELEMENT(m, n) C[(i + (m)) + (j + (n)) * ldc]

#define PACK_ELEMENT_K(m, offset_k) packed_a[(k + offset_k) * v_size2 + m]
#define PACK_ELEMENT(m) PACK_ELEMENT_K(m, 0)

// ASIMD
#define DECLARE_RESULT_VECTOR4(m, n)                                           \
  float32x4_t result##m##n = vdupq_n_f32(0.0);
#define DECLARE_RESULT(m, n) float32_t result##m##n = 0.0;
#define BROADCAST_LOAD_A4(m, offset_k)                                         \
  float32x4_t a##m##_k##offset_k = vld1q_dup_f32(&A_ELEMENT_K(m, offset_k));
#define LOAD_A1(m, offset_k)                                                   \
  float32_t a##m##_k##offset_k = A_ELEMENT_K(m, offset_k);
#define VECTOR_LOAD_B4(n, offset_k)                                            \
  float32x4_t b##n##_k##offset_k = vld1q_f32(&B_ELEMENT_K(n, offset_k));
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
#define PACK_B0(n, offset_k)                                                   \
  PACK_ELEMENT_K(n, offset_k) = vget_lane_f32(b##n##_k##offset_k, 0);
#define UPDATE_RESULT_VECTOR4(m, n, offset_k)                                  \
  result##m##n =                                                               \
    vfmaq_f32(result##m##n, a##m##_k##offset_k, b##n##_k##offset_k);
#define UPDATE_RESULT(m, n, offset_k)                                          \
  result##m##n = result##m##n + a##m##_k##offset_k * b##n##_k##offset_k;
#ifdef B0
#define VECTOR_STORE4(m, n)                                                    \
  vst1q_f32(&C_ELEMENT(m, n), vmulq_f32(result##m##n, vdupq_n_f32(alpha)));
#define STORE(m, n) C_ELEMENT(m, n) = alpha * result##m##n;
#else
#define VECTOR_STORE4(m, n)                                                    \
  result##m##n = vmulq_f32(result##m##n, vdupq_n_f32(alpha));                  \
  result##m##n =                                                               \
    vfmaq_f32(result##m##n, vld1q_f32(&C_ELEMENT(m, n)), vdupq_n_f32(beta));   \
  vst1q_f32(&C_ELEMENT(m, n), result##m##n);
#define STORE(m, n)                                                            \
  C_ELEMENT(m, n) = C_ELEMENT(m, n) * beta + alpha * result##m##n;
#endif

// SVE
#define DECLARE_RESULT_VECTOR(m, n) svfloat32_t result##m##n = svdup_f32(0.0);
#define BROADCAST_LOAD_A(m, offset_k)                                          \
  svfloat32_t a##s##m##_k##offset_k = svdup_f32(A_ELEMENT_K(m, offset_k));
#define BROADCAST_LOAD_B(n, offset_k)                                          \
  svfloat32_t b##s##n##_k##offset_k = svdup_f32(B_ELEMENT_K(n, offset_k));

#define QUADWORD_LOAD_B(n, offset_k)                                           \
  svfloat32_t b##s##n##_k##offset_k =                                          \
    svld1rq(pg_true, &B_ELEMENT_K(n, offset_k));
#define GATHER_LOAD_A(pg, offset_m, m, offset_k)                                         \
  svfloat32_t a##s##m##_k##offset_k =                                          \
    svld1_gather_index(pg, &A_ELEMENT_K(offset_m, offset_k), lda_vec);
#define PACK_A(m, offset_k)                                                    \
  svst1(pg_first, &PACK_ELEMENT_K(m, offset_k), a##s##m##_k##offset_k);
#define VECTOR_PACK_A(offset_m, m, offset_k)                                             \
  svst1(pg_true, &PACK_ELEMENT_K(m* v_size, offset_k), a##s##m##_k##offset_k);
#define QUADWORD_PACK_A(m, offset_k)                                           \
  svst1(pg_quad, &PACK_ELEMENT_K(m, offset_k), a##s##m##_k##offset_k);
#define UNPACK_VECTOR_A(offset_m, m, offset_k)                                           \
  svfloat32_t a##s##m##_k##offset_k =                                          \
    svld1(pg_true, &PACK_ELEMENT_K(offset_m, offset_k));
#define UNPACK_BROADCAST_A(m, offset_k)                                        \
  svfloat32_t a##s##m##_k##offset_k = svdup_f32(PACK_ELEMENT_K(m, offset_k));
#define UNPACK_QUADWORD_A(m, offset_k)                                         \
  svfloat32_t a##s##m##_k##offset_k =                                          \
    svld1rq(pg_true, &PACK_ELEMENT_K(m, offset_k));
#define UPDATE_RESULT_VECTOR(pg, m, n, offset_k)                               \
  result##m##n =                                                               \
    svmla_m(pg, result##m##n, a##s##m##_k##offset_k, b##s##n##_k##offset_k);
#define UPDATE_RESULT_VECTOR_QUADWORD(m, n, outer, lane, offset_k)             \
  result##m##n = svmla_lane(                                                   \
    result##m##n, a##s##m##_k##offset_k, b##s##outer##_k##offset_k, lane);
#ifdef B0
#define VECTOR_STORE(pg, offset_m, m, n)                                                 \
  result##m##n = svmul_m(pg, result##m##n, alpha_vec);                         \
  svst1(pg, &C_ELEMENT(offset_m, n), result##m##n);
#define SCATTER_STORE(pg, m, n)                                                \
  result##m##n = svmul_m(pg, result##m##n, alpha_vec);                         \
  svst1_scatter_index(                                                         \
    pg, &C_ELEMENT(v_size* m, n), svindex_u32(0LL, ldc), result##m##n);
#else
#define VECTOR_STORE(pg, offset_m, m, n)                                                 \
  result##m##n = svmul_m(pg, result##m##n, alpha_vec);                         \
  result##m##n =                                                               \
    svmla_m(pg, result##m##n, svld1(pg, &C_ELEMENT(offset_m, n)), beta_vec); \
  svst1(pg, &C_ELEMENT(offset_m, n), result##m##n);
#define SCATTER_STORE(pg, m, n)                                                \
  result##m##n = svmul_m(pg, result##m##n, alpha_vec);                         \
  result##m##n = svmla_m(                                                      \
    pg,                                                                        \
    result##m##n,                                                              \
    svld1_gather_index(pg, &C_ELEMENT(v_size * m, n), svindex_u32(0LL, ldc)),  \
    beta_vec);                                                                 \
  svst1_scatter_index(                                                         \
    pg, &C_ELEMENT(v_size* m, n), svindex_u32(0LL, ldc), result##m##n);
#endif

#ifndef LIKELY
#ifdef __GNUC__
#define LIKELY(x) __builtin_expect(!!(x), 1)
#else
#define LIKELY(x) (x)
#endif
#endif
#ifndef UNLIKELY
#ifdef __GNUC__
#define UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define UNLIKELY(x) (x)
#endif
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
  const uint64_t v_size = svcntw();
  const uint64_t v_size2 = v_size * 2;
  const svbool_t pg_true = svptrue_b32();
  const svfloat32_t alpha_vec = svdup_f32(alpha);
#ifndef B0
  const svfloat32_t beta_vec = svdup_f32(beta);
#endif
  const svuint32_t lda_vec = svindex_u32(0LL, lda);

  const BLASLONG v_m2 = M & -v_size2;
  const BLASLONG v_m1 = M & -v_size;
  const BLASLONG n8 = N & -8;
  const BLASLONG n4 = N & -4;

  const int pack_a = M >= v_size2 && N >= 8 && K >= 8 ? 1 : 0;
  FLOAT* packed_a =
    (pack_a) ? packed_a = (FLOAT*)malloc(K * 2 * v_size * sizeof(FLOAT)) : NULL;

  BLASLONG i = 0;
  for (; i < v_m2; i += v_size2) {

    BLASLONG j = 0;
    for (; j < n8; j += 8) {

      BLASLONG k = 0;
      DECLARE_RESULT_VECTOR(0, 0);
      DECLARE_RESULT_VECTOR(0, 1);
      DECLARE_RESULT_VECTOR(0, 2);
      DECLARE_RESULT_VECTOR(0, 3);
      DECLARE_RESULT_VECTOR(0, 4);
      DECLARE_RESULT_VECTOR(0, 5);
      DECLARE_RESULT_VECTOR(0, 6);
      DECLARE_RESULT_VECTOR(0, 7);
      DECLARE_RESULT_VECTOR(1, 0);
      DECLARE_RESULT_VECTOR(1, 1);
      DECLARE_RESULT_VECTOR(1, 2);
      DECLARE_RESULT_VECTOR(1, 3);
      DECLARE_RESULT_VECTOR(1, 4);
      DECLARE_RESULT_VECTOR(1, 5);
      DECLARE_RESULT_VECTOR(1, 6);
      DECLARE_RESULT_VECTOR(1, 7);

      if (LIKELY(packed_a != NULL)) {
        if (j == 0) {
          for (; k < K; k++) {

            QUADWORD_LOAD_B(0, 0);
            GATHER_LOAD_A(pg_true, 0, 0, 0);
            VECTOR_PACK_A(0, 0, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(0, 0, 0, 0, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(0, 1, 0, 1, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(0, 2, 0, 2, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(0, 3, 0, 3, 0);
            QUADWORD_LOAD_B(4, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(0, 4, 4, 0, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(0, 5, 4, 1, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(0, 6, 4, 2, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(0, 7, 4, 3, 0);
            GATHER_LOAD_A(pg_true, v_size, 1, 0);
            VECTOR_PACK_A(v_size, 1, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(1, 0, 0, 0, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(1, 1, 0, 1, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(1, 2, 0, 2, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(1, 3, 0, 3, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(1, 4, 4, 0, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(1, 5, 4, 1, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(1, 6, 4, 2, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(1, 7, 4, 3, 0);
          }
        } else {
          for (; k < K; k++) {

            QUADWORD_LOAD_B(0, 0);
            UNPACK_VECTOR_A(0, 0, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(0, 0, 0, 0, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(0, 1, 0, 1, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(0, 2, 0, 2, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(0, 3, 0, 3, 0);
            QUADWORD_LOAD_B(4, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(0, 4, 4, 0, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(0, 5, 4, 1, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(0, 6, 4, 2, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(0, 7, 4, 3, 0);
            UNPACK_VECTOR_A(v_size, 1, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(1, 0, 0, 0, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(1, 1, 0, 1, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(1, 2, 0, 2, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(1, 3, 0, 3, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(1, 4, 4, 0, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(1, 5, 4, 1, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(1, 6, 4, 2, 0);
            UPDATE_RESULT_VECTOR_QUADWORD(1, 7, 4, 3, 0);
          }
        }
      } else {
        for (; k < K; k++) {

          QUADWORD_LOAD_B(0, 0);
          GATHER_LOAD_A(pg_true, 0, 0, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(0, 0, 0, 0, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(0, 1, 0, 1, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(0, 2, 0, 2, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(0, 3, 0, 3, 0);
          QUADWORD_LOAD_B(4, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(0, 4, 4, 0, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(0, 5, 4, 1, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(0, 6, 4, 2, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(0, 7, 4, 3, 0);
          GATHER_LOAD_A(pg_true, v_size, 1, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(1, 0, 0, 0, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(1, 1, 0, 1, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(1, 2, 0, 2, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(1, 3, 0, 3, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(1, 4, 4, 0, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(1, 5, 4, 1, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(1, 6, 4, 2, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(1, 7, 4, 3, 0);
        }
      }
      VECTOR_STORE(pg_true, 0, 0, 0);
      VECTOR_STORE(pg_true, 0, 0, 1);
      VECTOR_STORE(pg_true, 0, 0, 2);
      VECTOR_STORE(pg_true, 0, 0, 3);
      VECTOR_STORE(pg_true, 0, 0, 4);
      VECTOR_STORE(pg_true, 0, 0, 5);
      VECTOR_STORE(pg_true, 0, 0, 6);
      VECTOR_STORE(pg_true, 0, 0, 7);
      VECTOR_STORE(pg_true, v_size, 1, 0);
      VECTOR_STORE(pg_true, v_size, 1, 1);
      VECTOR_STORE(pg_true, v_size, 1, 2);
      VECTOR_STORE(pg_true, v_size, 1, 3);
      VECTOR_STORE(pg_true, v_size, 1, 4);
      VECTOR_STORE(pg_true, v_size, 1, 5);
      VECTOR_STORE(pg_true, v_size, 1, 6);
      VECTOR_STORE(pg_true, v_size, 1, 7);
    }
    for (; j < n4; j += 4) {

      BLASLONG k = 0;
      DECLARE_RESULT_VECTOR(0, 0);
      DECLARE_RESULT_VECTOR(0, 1);
      DECLARE_RESULT_VECTOR(0, 2);
      DECLARE_RESULT_VECTOR(0, 3);
      DECLARE_RESULT_VECTOR(1, 0);
      DECLARE_RESULT_VECTOR(1, 1);
      DECLARE_RESULT_VECTOR(1, 2);
      DECLARE_RESULT_VECTOR(1, 3);

      if (LIKELY(packed_a != NULL)) {
        for (; k < K; k++) {

          QUADWORD_LOAD_B(0, 0);
          UNPACK_VECTOR_A(0, 0, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(0, 0, 0, 0, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(0, 1, 0, 1, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(0, 2, 0, 2, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(0, 3, 0, 3, 0);
          UNPACK_VECTOR_A(v_size, 1, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(1, 0, 0, 0, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(1, 1, 0, 1, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(1, 2, 0, 2, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(1, 3, 0, 3, 0);
        }
      } else {
        for (; k < K; k++) {

          QUADWORD_LOAD_B(0, 0);
          GATHER_LOAD_A(pg_true, 0, 0, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(0, 0, 0, 0, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(0, 1, 0, 1, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(0, 2, 0, 2, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(0, 3, 0, 3, 0);
          GATHER_LOAD_A(pg_true, v_size, 1, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(1, 0, 0, 0, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(1, 1, 0, 1, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(1, 2, 0, 2, 0);
          UPDATE_RESULT_VECTOR_QUADWORD(1, 3, 0, 3, 0);
        }
      }
      VECTOR_STORE(pg_true, 0, 0, 0);
      VECTOR_STORE(pg_true, 0, 0, 1);
      VECTOR_STORE(pg_true, 0, 0, 2);
      VECTOR_STORE(pg_true, 0, 0, 3);
      VECTOR_STORE(pg_true, v_size, 1, 0);
      VECTOR_STORE(pg_true, v_size, 1, 1);
      VECTOR_STORE(pg_true, v_size, 1, 2);
      VECTOR_STORE(pg_true, v_size, 1, 3);
    }
    for (; j < N; j++) {

      BLASLONG k = 0;
      DECLARE_RESULT_VECTOR(0, 0);
      DECLARE_RESULT_VECTOR(1, 0);

      if (LIKELY(packed_a != NULL)) {
        for (; k < K; k++) {

          BROADCAST_LOAD_B(0, 0);
          UNPACK_VECTOR_A(0, 0, 0);
          UPDATE_RESULT_VECTOR(pg_true, 0, 0, 0);
          UNPACK_VECTOR_A(v_size, 1, 0);
          UPDATE_RESULT_VECTOR(pg_true, 1, 0, 0);
        }
      } else {
        for (; k < K; k++) {

          BROADCAST_LOAD_B(0, 0);
          GATHER_LOAD_A(pg_true, 0, 0, 0);
          UPDATE_RESULT_VECTOR(pg_true, 0, 0, 0);
          GATHER_LOAD_A(pg_true, v_size, 1, 0);
          UPDATE_RESULT_VECTOR(pg_true, 1, 0, 0);
        }
      }
      VECTOR_STORE(pg_true, 0, 0, 0);
      VECTOR_STORE(pg_true, v_size, 1, 0);
    }
  }
  for (; i < v_m1; i += v_size) {

    BLASLONG j = 0;
    for (; j < n8; j += 8) {

      BLASLONG k = 0;
      DECLARE_RESULT_VECTOR(0, 0);
      DECLARE_RESULT_VECTOR(0, 1);
      DECLARE_RESULT_VECTOR(0, 2);
      DECLARE_RESULT_VECTOR(0, 3);
      DECLARE_RESULT_VECTOR(0, 4);
      DECLARE_RESULT_VECTOR(0, 5);
      DECLARE_RESULT_VECTOR(0, 6);
      DECLARE_RESULT_VECTOR(0, 7);

      for (; k < K; k++) {

        QUADWORD_LOAD_B(0, 0);
        GATHER_LOAD_A(pg_true, 0, 0, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 0, 0, 0, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 1, 0, 1, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 2, 0, 2, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 3, 0, 3, 0);
        QUADWORD_LOAD_B(4, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 4, 4, 0, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 5, 4, 1, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 6, 4, 2, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 7, 4, 3, 0);
      }
      VECTOR_STORE(pg_true, 0, 0, 0);
      VECTOR_STORE(pg_true, 0, 0, 1);
      VECTOR_STORE(pg_true, 0, 0, 2);
      VECTOR_STORE(pg_true, 0, 0, 3);
      VECTOR_STORE(pg_true, 0, 0, 4);
      VECTOR_STORE(pg_true, 0, 0, 5);
      VECTOR_STORE(pg_true, 0, 0, 6);
      VECTOR_STORE(pg_true, 0, 0, 7);
    }
    for (; j < n4; j += 4) {

      BLASLONG k = 0;
      DECLARE_RESULT_VECTOR(0, 0);
      DECLARE_RESULT_VECTOR(0, 1);
      DECLARE_RESULT_VECTOR(0, 2);
      DECLARE_RESULT_VECTOR(0, 3);

      for (; k < K; k++) {

        QUADWORD_LOAD_B(0, 0);
        GATHER_LOAD_A(pg_true, 0, 0, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 0, 0, 0, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 1, 0, 1, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 2, 0, 2, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 3, 0, 3, 0);
      }
      VECTOR_STORE(pg_true, 0, 0, 0);
      VECTOR_STORE(pg_true, 0, 0, 1);
      VECTOR_STORE(pg_true, 0, 0, 2);
      VECTOR_STORE(pg_true, 0, 0, 3);
    }
    for (; j < N; j++) {

      BLASLONG k = 0;
      DECLARE_RESULT_VECTOR(0, 0);

      for (; k < K; k++) {

        BROADCAST_LOAD_B(0, 0);
        GATHER_LOAD_A(pg_true, 0, 0, 0);
        UPDATE_RESULT_VECTOR(pg_true, 0, 0, 0);
      }
      VECTOR_STORE(pg_true, 0, 0, 0);
    }
  }
  for (; i < M; i += v_size) {
    const svbool_t pg_tail = svwhilelt_b32((uint32_t)i, (uint32_t)(M));

    BLASLONG j = 0;
    for (; j < n8; j += 8) {

      BLASLONG k = 0;
      DECLARE_RESULT_VECTOR(0, 0);
      DECLARE_RESULT_VECTOR(0, 1);
      DECLARE_RESULT_VECTOR(0, 2);
      DECLARE_RESULT_VECTOR(0, 3);
      DECLARE_RESULT_VECTOR(0, 4);
      DECLARE_RESULT_VECTOR(0, 5);
      DECLARE_RESULT_VECTOR(0, 6);
      DECLARE_RESULT_VECTOR(0, 7);

      for (; k < K; k++) {

        QUADWORD_LOAD_B(0, 0);
        GATHER_LOAD_A(pg_tail, 0, 0, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 0, 0, 0, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 1, 0, 1, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 2, 0, 2, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 3, 0, 3, 0);
        QUADWORD_LOAD_B(4, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 4, 4, 0, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 5, 4, 1, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 6, 4, 2, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 7, 4, 3, 0);
      }
      VECTOR_STORE(pg_tail, 0, 0, 0);
      VECTOR_STORE(pg_tail, 0, 0, 1);
      VECTOR_STORE(pg_tail, 0, 0, 2);
      VECTOR_STORE(pg_tail, 0, 0, 3);
      VECTOR_STORE(pg_tail, 0, 0, 4);
      VECTOR_STORE(pg_tail, 0, 0, 5);
      VECTOR_STORE(pg_tail, 0, 0, 6);
      VECTOR_STORE(pg_tail, 0, 0, 7);
    }
    for (; j < n4; j += 4) {

      BLASLONG k = 0;
      DECLARE_RESULT_VECTOR(0, 0);
      DECLARE_RESULT_VECTOR(0, 1);
      DECLARE_RESULT_VECTOR(0, 2);
      DECLARE_RESULT_VECTOR(0, 3);

      for (; k < K; k++) {

        QUADWORD_LOAD_B(0, 0);
        GATHER_LOAD_A(pg_tail, 0, 0, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 0, 0, 0, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 1, 0, 1, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 2, 0, 2, 0);
        UPDATE_RESULT_VECTOR_QUADWORD(0, 3, 0, 3, 0);
      }
      VECTOR_STORE(pg_tail, 0, 0, 0);
      VECTOR_STORE(pg_tail, 0, 0, 1);
      VECTOR_STORE(pg_tail, 0, 0, 2);
      VECTOR_STORE(pg_tail, 0, 0, 3);
    }
    for (; j < N; j++) {

      BLASLONG k = 0;
      DECLARE_RESULT_VECTOR(0, 0);

      for (; k < K; k++) {

        BROADCAST_LOAD_B(0, 0);
        GATHER_LOAD_A(pg_tail, 0, 0, 0);
        UPDATE_RESULT_VECTOR(pg_tail, 0, 0, 0);
      }
      VECTOR_STORE(pg_tail, 0, 0, 0);
    }
  }

  if (pack_a)
    free(packed_a);

  return 0;
}
