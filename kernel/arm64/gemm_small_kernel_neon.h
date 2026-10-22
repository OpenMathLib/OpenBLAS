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
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/


/* Register-blocked NEON GEMM for small problems (SMALL_MATRIX_OPT): blocks of C in registers, no packing. */

#ifndef GEMM_SMALL_KERNEL_NEON_H
#define GEMM_SMALL_KERNEL_NEON_H

#include <arm_neon.h>
#ifdef DOUBLE
#define SK_NV 2
typedef float64x2_t sk_v;
#define SK_LD vld1q_f64
#define SK_ST vst1q_f64
#define SK_DUP vdupq_n_f64
#define SK_FMA_L(acc, a, b, l) vfmaq_laneq_f64(acc, a, b, l)
#define SK_MUL vmulq_f64
#define SK_FMA vfmaq_f64
#else
#define SK_NV 4
typedef float32x4_t sk_v;
#define SK_LD vld1q_f32
#define SK_ST vst1q_f32
#define SK_DUP vdupq_n_f32
#define SK_FMA_L(acc, a, b, l) vfmaq_laneq_f32(acc, a, b, l)
#define SK_MUL vmulq_f32
#define SK_FMA vfmaq_f32
#endif
#define SK_MR (2 * SK_NV)
#define SK_NR 4

#define SK_INL static inline __attribute__((always_inline))

SK_INL sk_v sk_fmal(sk_v acc, sk_v a, sk_v b, const int t) {
  switch (t) {
    case 0: return SK_FMA_L(acc, a, b, 0);
#if SK_NV > 2
    case 2: return SK_FMA_L(acc, a, b, 2);
    case 3: return SK_FMA_L(acc, a, b, 3);
#endif
    default: return SK_FMA_L(acc, a, b, 1);
  }
}

/* C block of MV vectors x NR columns at (a, b, c); the first skip rows of the block are left unchanged in C. */
/* As in the reference GEMM for A op(B), alpha scales op(B) before the products (no overflow of A B for a small alpha). */
SK_INL void sk_blk(const int TB, const int MV, const int NR, int skip, BLASLONG k, const FLOAT *a, BLASLONG lda,
                   const FLOAT *b, BLASLONG ldb, FLOAT alpha, FLOAT beta, FLOAT *c, BLASLONG ldc) {
  const sk_v va = SK_DUP(alpha), vb = SK_DUP(beta);
  sk_v acc[4][2];
  for (int j = 0; j < NR; ++j)
    for (int h = 0; h < MV; ++h) acc[j][h] = SK_DUP(0);
  BLASLONG l = 0;
  if (!TB) /* B(l, j) at b[l + j * ldb]: SK_NV depth steps per load of each column */
    for (; l + SK_NV <= k; l += SK_NV) {
      sk_v bv[4];
      for (int j = 0; j < NR; ++j) bv[j] = SK_MUL(SK_LD(b + l + j * ldb), va);
#pragma clang loop unroll(full)
      for (int t = 0; t < SK_NV; ++t) {
        sk_v av[2];
        for (int h = 0; h < MV; ++h) av[h] = SK_LD(a + (l + t) * lda + h * SK_NV);
#pragma clang loop unroll(full)
        for (int j = 0; j < NR; ++j)
          for (int h = 0; h < MV; ++h) acc[j][h] = sk_fmal(acc[j][h], av[h], bv[j], t);
      }
    }
  else if (NR == 4) /* B(l, j) at b[j + l * ldb]: the four columns of one depth step in 4 / SK_NV loads */
    for (; l < k; ++l) {
      sk_v av[2], bv[4 / SK_NV];
      for (int h = 0; h < MV; ++h) av[h] = SK_LD(a + l * lda + h * SK_NV);
      for (int q = 0; q < 4 / SK_NV; ++q) bv[q] = SK_MUL(SK_LD(b + l * ldb + q * SK_NV), va);
#pragma clang loop unroll(full)
      for (int j = 0; j < 4; ++j)
        for (int h = 0; h < MV; ++h) acc[j][h] = sk_fmal(acc[j][h], av[h], bv[j / SK_NV], j % SK_NV);
    }
  const BLASLONG sl = TB ? ldb : 1, sj = TB ? 1 : ldb;
  for (; l < k; ++l) {
    sk_v av[2];
    for (int h = 0; h < MV; ++h) av[h] = SK_LD(a + l * lda + h * SK_NV);
    for (int j = 0; j < NR; ++j) {
      const sk_v bj = SK_DUP(alpha * b[l * sl + j * sj]);
      for (int h = 0; h < MV; ++h) acc[j][h] = SK_FMA(acc[j][h], av[h], bj);
    }
  }
  for (int j = 0; j < NR; ++j) {
    FLOAT *cj = c + j * ldc;
    for (int h = 0; h < MV; ++h) {
      FLOAT *ch = cj + h * SK_NV;
      const sk_v r = beta == (FLOAT)0 ? acc[j][h] : SK_FMA(acc[j][h], SK_LD(ch), vb);
      if (skip == 0) SK_ST(ch, r);
      else {
        FLOAT tmp[SK_NV];
        SK_ST(tmp, r);
        for (int i = skip; i < SK_NV; ++i) ch[i] = tmp[i];
      }
    }
  }
}

/* Rows fewer than one vector: plain loops. */
static void sk_rows(int tb, BLASLONG m, BLASLONG n, BLASLONG k, const FLOAT *a, BLASLONG lda, const FLOAT *b,
                   BLASLONG ldb, FLOAT alpha, FLOAT beta, FLOAT *c, BLASLONG ldc) {
  const BLASLONG sl = tb ? ldb : 1, sj = tb ? 1 : ldb;
  for (BLASLONG j = 0; j < n; ++j)
    for (BLASLONG i = 0; i < m; ++i) {
      FLOAT s = 0;
      for (BLASLONG l = 0; l < k; ++l) s += a[i + l * lda] * (alpha * b[l * sl + j * sj]);
      c[i + j * ldc] = s + (beta == (FLOAT)0 ? (FLOAT)0 : beta * c[i + j * ldc]);
    }
}

/* One block column of NR columns: two-vector row blocks, one-vector blocks, and the row tail as an overlapped block. */
SK_INL void sk_col(const int TB, const int NR, BLASLONG m, BLASLONG k, const FLOAT *a, BLASLONG lda, const FLOAT *b, BLASLONG ldb,
                 FLOAT alpha, FLOAT beta, FLOAT *c, BLASLONG ldc) {
  BLASLONG i = 0;
  for (; i + 2 * SK_NV <= m; i += 2 * SK_NV) sk_blk(TB, 2, NR, 0, k, a + i, lda, b, ldb, alpha, beta, c + i, ldc);
  if (i + SK_NV <= m) {
    sk_blk(TB, 1, NR, 0, k, a + i, lda, b, ldb, alpha, beta, c + i, ldc);
    i += SK_NV;
  }
  if (i < m) {
    if (m >= SK_NV) sk_blk(TB, 1, NR, (int)(SK_NV - (m - i)), k, a + m - SK_NV, lda, b, ldb, alpha, beta, c + m - SK_NV, ldc);
    else sk_rows(TB, m, NR, k, a, lda, b, ldb, alpha, beta, c, ldc);
  }
}

/* C = alpha A op(B) + beta C with A column-major; tb: B transposed. */
static void sk_gemm_n(int tb, BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, const FLOAT *a, BLASLONG lda, const FLOAT *b,
                     BLASLONG ldb, FLOAT beta, FLOAT *c, BLASLONG ldc) {
  const BLASLONG sj = tb ? 1 : ldb;
  BLASLONG j = 0;
  if (tb)
    for (; j + 4 <= n; j += 4) sk_col(1, 4, m, k, a, lda, b + j * sj, ldb, alpha, beta, c + j * ldc, ldc);
  else
    for (; j + 4 <= n; j += 4) sk_col(0, 4, m, k, a, lda, b + j * sj, ldb, alpha, beta, c + j * ldc, ldc);
#define SK_TAIL(NR)                                                                                                   \
  if (tb) sk_col(1, NR, m, k, a, lda, b + j * sj, ldb, alpha, beta, c + j * ldc, ldc);                                \
  else sk_col(0, NR, m, k, a, lda, b + j * sj, ldb, alpha, beta, c + j * ldc, ldc);
  switch (n - j) {
    case 3: SK_TAIL(3) break;
    case 2: SK_TAIL(2) break;
    case 1: SK_TAIL(1) break;
  }
#undef SK_TAIL
}

#ifdef DOUBLE
#define SK_ADDV vaddvq_f64
#else
#define SK_ADDV vaddvq_f32
#endif

/* C block of R rows x NR columns for C = alpha A^T B + beta C: dot products along the contiguous depth. */
SK_INL void sk_dot(const int R, const int NR, BLASLONG k, const FLOAT *a, BLASLONG lda, const FLOAT *b, BLASLONG ldb,
                 FLOAT alpha, FLOAT beta, FLOAT *c, BLASLONG ldc) {
  sk_v acc[4][4];
  for (int r = 0; r < R; ++r)
    for (int j = 0; j < NR; ++j) acc[r][j] = SK_DUP(0);
  BLASLONG l = 0;
  for (; l + SK_NV <= k; l += SK_NV) {
    sk_v av[4], bv[4];
    for (int r = 0; r < R; ++r) av[r] = SK_LD(a + r * lda + l);
    for (int j = 0; j < NR; ++j) bv[j] = SK_LD(b + j * ldb + l);
    for (int r = 0; r < R; ++r)
      for (int j = 0; j < NR; ++j) acc[r][j] = SK_FMA(acc[r][j], av[r], bv[j]);
  }
  for (int r = 0; r < R; ++r)
    for (int j = 0; j < NR; ++j) {
      FLOAT s = SK_ADDV(acc[r][j]);
      for (BLASLONG q = l; q < k; ++q) s += a[r * lda + q] * b[j * ldb + q];
      FLOAT *cc = c + r + j * ldc;
      *cc = alpha * s + (beta == (FLOAT)0 ? (FLOAT)0 : beta * *cc);
    }
}

SK_INL void sk_dot_rows(const int NR, BLASLONG m, BLASLONG k, const FLOAT *a, BLASLONG lda, const FLOAT *b, BLASLONG ldb,
                      FLOAT alpha, FLOAT beta, FLOAT *c, BLASLONG ldc) {
  BLASLONG i = 0;
  for (; i + 4 <= m; i += 4) sk_dot(4, NR, k, a + i * lda, lda, b, ldb, alpha, beta, c + i, ldc);
  switch (m - i) {
    case 3: sk_dot(3, NR, k, a + i * lda, lda, b, ldb, alpha, beta, c + i, ldc); break;
    case 2: sk_dot(2, NR, k, a + i * lda, lda, b, ldb, alpha, beta, c + i, ldc); break;
    case 1: sk_dot(1, NR, k, a + i * lda, lda, b, ldb, alpha, beta, c + i, ldc); break;
  }
}

/* C = alpha A^T B + beta C. */
static void sk_gemm_tn(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, const FLOAT *a, BLASLONG lda, const FLOAT *b,
                      BLASLONG ldb, FLOAT beta, FLOAT *c, BLASLONG ldc) {
  BLASLONG j = 0;
  for (; j + 4 <= n; j += 4) sk_dot_rows(4, m, k, a, lda, b + j * ldb, ldb, alpha, beta, c + j * ldc, ldc);
  switch (n - j) {
    case 3: sk_dot_rows(3, m, k, a, lda, b + j * ldb, ldb, alpha, beta, c + j * ldc, ldc); break;
    case 2: sk_dot_rows(2, m, k, a, lda, b + j * ldb, ldb, alpha, beta, c + j * ldc, ldc); break;
    case 1: sk_dot_rows(1, m, k, a, lda, b + j * ldb, ldb, alpha, beta, c + j * ldc, ldc); break;
  }
}

/* C = alpha A^T B^T + beta C: plain loops. */
static void sk_gemm_tt(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, const FLOAT *a, BLASLONG lda, const FLOAT *b,
                       BLASLONG ldb, FLOAT beta, FLOAT *c, BLASLONG ldc) {
  for (BLASLONG j = 0; j < n; ++j)
    for (BLASLONG i = 0; i < m; ++i) {
      FLOAT s = 0;
      for (BLASLONG l = 0; l < k; ++l) s += a[l + i * lda] * b[j + l * ldb];
      FLOAT *cc = c + i + j * ldc;
      *cc = alpha * s + (beta == (FLOAT)0 ? (FLOAT)0 : beta * *cc);
    }
}

/* The small-matrix kernel interface; B0: beta = 0, C not read. */
#ifdef B0
#define SK_KERNEL(call)                                                                                              \
  int CNAME(BLASLONG M, BLASLONG N, BLASLONG K, IFLOAT *A, BLASLONG lda, FLOAT alpha, IFLOAT *B, BLASLONG ldb,      \
            FLOAT *C, BLASLONG ldc) {                                                                               \
    const FLOAT beta = ZERO;                                                                                         \
    call;                                                                                                            \
    return 0;                                                                                                        \
  }
#else
#define SK_KERNEL(call)                                                                                              \
  int CNAME(BLASLONG M, BLASLONG N, BLASLONG K, IFLOAT *A, BLASLONG lda, FLOAT alpha, IFLOAT *B, BLASLONG ldb,      \
            FLOAT beta, FLOAT *C, BLASLONG ldc) {                                                                    \
    call;                                                                                                            \
    return 0;                                                                                                        \
  }
#endif

#endif
