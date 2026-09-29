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

/*
 * SYMM, SYRK, SYR2K, TRMM and TRSM on arm64 SME targets, on top of the SME GEMM kernel (SME_[SD]GEMM_KERNEL).
 *
 * The level-3 driver runs these routines on the NEON GEMM kernel of the target. Here each routine splits its
 * symmetric or triangular dimension in two, recursively: the off-diagonal blocks are GEMM calls on the SME
 * kernel, and only blocks of at most S3_NB rows or columns remain. Those run as one GEMM on a full copy of
 * the block (SYMM, SYRK, SYR2K, TRMM), or, for TRSM, as substitution in C on blocks of at most S3_NB_TRSM, so
 * that every solve keeps the arithmetic of substitution (the generic TRSM kernels of these targets run at 10-16
 * GFLOPS on such thin blocks). All arguments are in the column-major form that the interface files build.
 *
 * Included by interface/{symm,syrk,syr2k,trsm}.c; they call the s3_*_hook functions at the end of this file,
 * which return 1 when they have handled the call.
 */

#if !defined(COMPLEX) && !defined(XDOUBLE) && !defined(BFLOAT16) && !defined(HFLOAT16) && defined(ARCH_ARM64) && \
    (defined(USE_SGEMM_KERNEL_DIRECT) || defined(DYNAMIC_ARCH))
#define SME_LEVEL3 1

#include <string.h>
#include <arm_neon.h>

#define S3_NB 64          /* largest diagonal block; also the split granularity */
#define S3_MIN_WORK 2e5   /* multiply-adds below which the level-3 driver keeps the call */
#define S3_BCHUNK 4096    /* columns (rows) of B per temporary copy in the TRMM base case */
#ifndef S3_NB_TRSM
#define S3_NB_TRSM 32     /* TRSM blocks solved by substitution */
#endif

#ifdef DYNAMIC_ARCH
extern char *gotoblas_corename(void);
#endif

/* Whether to take this path: large enough, and the SME GEMM kernel is what gemm.c calls on this core. */
static inline int s3_enabled(double work) {
  if (work < S3_MIN_WORK) return 0;
#if defined(DYNAMIC_ARCH)
  const char *c = gotoblas_corename();
  if (strcmp(c, "armv9sme") != 0
#if defined(__clang__)
      && strcmp(c, "vortexm4") != 0
#endif
  )
    return 0;
#endif
  return 1;
}

/* Column-major C = alpha op(A) op(B) + beta C on the SME kernel. */
static inline void s3_gemm(int ta, int tb, BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, FLOAT *a, BLASLONG lda,
                    FLOAT *b, BLASLONG ldb, FLOAT beta, FLOAT *c, BLASLONG ldc) {
  if (m <= 0 || n <= 0) return;
  char *TA = ta ? "T" : "N", *TB = tb ? "T" : "N";
#ifdef DOUBLE
  SME_DGEMM_KERNEL(TA, TB, m, n, k, &alpha, a, lda, b, ldb, &beta, c, ldc);
#else
  SME_SGEMM_KERNEL(TA, TB, m, n, k, &alpha, a, lda, b, ldb, &beta, c, ldc);
#endif
}

/* First part of a split of t > S3_NB: about half, a multiple of S3_NB. */
static inline BLASLONG s3_split(BLASLONG t) { return (t / 2 + S3_NB - 1) / S3_NB * S3_NB; }

/* Address of element (i, j) of op(A), where op(A) is A (tr 0) or its transpose (tr 1). */
#define S3_OP(a, lda, tr, i, j) ((tr) ? (a) + (j) + (BLASLONG)(i) * (lda) : (a) + (i) + (BLASLONG)(j) * (lda))

/* ---- SYMM: C = alpha A B + beta C (side 0) or alpha B A + beta C (side 1), A symmetric ---- */

/* The t x t diagonal block of a symmetric matrix with both triangles filled. */
static inline void s3_sym_full(int lower, BLASLONG t, FLOAT *a, BLASLONG lda, FLOAT *f) {
  for (BLASLONG j = 0; j < t; ++j)
    for (BLASLONG i = 0; i < t; ++i) f[i + j * t] = ((i >= j) == lower) ? a[i + j * lda] : a[j + i * lda];
}

static inline void s3_symm_rec(int side, int lower, BLASLONG m, BLASLONG n, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *b,
                        BLASLONG ldb, FLOAT beta, FLOAT *c, BLASLONG ldc) {
  const BLASLONG t = side ? n : m;
  if (t <= S3_NB) {
    FLOAT f[S3_NB * S3_NB];
    s3_sym_full(lower, t, a, lda, f);
    if (!side) s3_gemm(0, 0, m, n, m, alpha, f, t, b, ldb, beta, c, ldc);
    else s3_gemm(0, 0, m, n, n, alpha, b, ldb, f, t, beta, c, ldc);
    return;
  }
  const BLASLONG t1 = s3_split(t), t2 = t - t1;
  FLOAT *a22 = a + t1 + t1 * lda;
  /* A21 is stored for lower; for upper it is the transpose of the stored A12 */
  FLOAT *a21 = lower ? a + t1 : a + t1 * lda;
  const int tr21 = !lower;
  if (!side) { /* C1 = A11 B1 + A21^T B2, C2 = A21 B1 + A22 B2 */
    s3_symm_rec(0, lower, t1, n, alpha, a, lda, b, ldb, beta, c, ldc);
    s3_gemm(!tr21, 0, t1, n, t2, alpha, a21, lda, b + t1, ldb, 1, c, ldc);
    s3_symm_rec(0, lower, t2, n, alpha, a22, lda, b + t1, ldb, beta, c + t1, ldc);
    s3_gemm(tr21, 0, t2, n, t1, alpha, a21, lda, b, ldb, 1, c + t1, ldc);
  } else { /* C1 = B1 A11 + B2 A21, C2 = B1 A21^T + B2 A22 */
    s3_symm_rec(1, lower, m, t1, alpha, a, lda, b, ldb, beta, c, ldc);
    s3_gemm(0, tr21, m, t1, t2, alpha, b + t1 * ldb, ldb, a21, lda, 1, c, ldc);
    s3_symm_rec(1, lower, m, t2, alpha, a22, lda, b + t1 * ldb, ldb, beta, c + t1 * ldc, ldc);
    s3_gemm(0, !tr21, m, t2, t1, alpha, b, ldb, a21, lda, 1, c + t1 * ldc, ldc);
  }
}

/* ---- SYRK (two 0) and SYR2K (two 1) on the lower or upper triangle of the n x n C ---- */

/* Rows r0.. of op(X): X is n x k for tr 0, k x n for tr 1. */
#define S3_ROWS(x, ld, tr, r0) ((tr) ? (x) + (BLASLONG)(r0) * (ld) : (x) + (r0))

/* C = alpha (op(A) op(B)^T + [two] op(B) op(A)^T) + beta C; SYRK passes B = A. */
static inline void s3_syrk_rec(int two, int lower, int tr, BLASLONG n, BLASLONG k, FLOAT alpha, FLOAT *a, BLASLONG lda,
                        FLOAT *b, BLASLONG ldb, FLOAT beta, FLOAT *c, BLASLONG ldc) {
  if (n <= S3_NB) {
    FLOAT f[S3_NB * S3_NB];
    s3_gemm(tr, !tr, n, n, k, alpha, a, lda, b, ldb, 0, f, n);
    if (two) s3_gemm(tr, !tr, n, n, k, alpha, b, ldb, a, lda, 1, f, n);
    for (BLASLONG j = 0; j < n; ++j)
      for (BLASLONG i = lower ? j : 0; i < (lower ? n : j + 1); ++i)
        c[i + j * ldc] = f[i + j * n] + (beta == 0 ? 0 : beta * c[i + j * ldc]);
    return;
  }
  const BLASLONG n1 = s3_split(n), n2 = n - n1;
  FLOAT *a2 = S3_ROWS(a, lda, tr, n1), *b2 = S3_ROWS(b, ldb, tr, n1);
  s3_syrk_rec(two, lower, tr, n1, k, alpha, a, lda, b, ldb, beta, c, ldc);
  if (lower) { /* C21 = alpha (A2 B1^T + [two] B2 A1^T) + beta C21 */
    s3_gemm(tr, !tr, n2, n1, k, alpha, a2, lda, b, ldb, beta, c + n1, ldc);
    if (two) s3_gemm(tr, !tr, n2, n1, k, alpha, b2, ldb, a, lda, 1, c + n1, ldc);
  } else { /* C12 = alpha (A1 B2^T + [two] B1 A2^T) + beta C12 */
    s3_gemm(tr, !tr, n1, n2, k, alpha, a, lda, b2, ldb, beta, c + n1 * ldc, ldc);
    if (two) s3_gemm(tr, !tr, n1, n2, k, alpha, b, ldb, a2, lda, 1, c + n1 * ldc, ldc);
  }
  s3_syrk_rec(two, lower, tr, n2, k, alpha, a2, lda, b2, ldb, beta, c + n1 + n1 * ldc, ldc);
}

/* ---- TRMM: B = alpha op(A) B (side 0) or alpha B op(A) (side 1), A triangular ---- */

/* The t x t diagonal block of op(A) as a full matrix: zeros outside the triangle, ones on a unit diagonal. */
static inline void s3_tri_full(int lower_op, int tr, int unit, BLASLONG t, FLOAT *a, BLASLONG lda, FLOAT *f) {
  for (BLASLONG j = 0; j < t; ++j)
    for (BLASLONG i = 0; i < t; ++i)
      f[i + j * t] = i == j ? (unit ? 1 : *S3_OP(a, lda, tr, i, i))
                            : (((i > j) == lower_op) ? *S3_OP(a, lda, tr, i, j) : 0);
}

/* lower_op: op(A) is lower triangular; work holds S3_NB * S3_BCHUNK elements. */
static inline void s3_trmm_rec(int side, int lower_op, int tr, int unit, BLASLONG m, BLASLONG n, FLOAT alpha, FLOAT *a,
                        BLASLONG lda, FLOAT *b, BLASLONG ldb, FLOAT *work) {
  const BLASLONG t = side ? n : m;
  if (t <= S3_NB) {
    FLOAT f[S3_NB * S3_NB];
    s3_tri_full(lower_op, tr, unit, t, a, lda, f);
    if (!side) { /* B = alpha F B, S3_BCHUNK columns at a time through a copy */
      for (BLASLONG j0 = 0; j0 < n; j0 += S3_BCHUNK) {
        const BLASLONG nc = n - j0 < S3_BCHUNK ? n - j0 : S3_BCHUNK;
        for (BLASLONG j = 0; j < nc; ++j) memcpy(work + j * t, b + (j0 + j) * ldb, t * sizeof(FLOAT));
        s3_gemm(0, 0, t, nc, t, alpha, f, t, work, t, 0, b + j0 * ldb, ldb);
      }
    } else { /* B = alpha B F, S3_BCHUNK rows at a time through a copy */
      for (BLASLONG i0 = 0; i0 < m; i0 += S3_BCHUNK) {
        const BLASLONG mc = m - i0 < S3_BCHUNK ? m - i0 : S3_BCHUNK;
        for (BLASLONG j = 0; j < t; ++j) memcpy(work + j * mc, b + i0 + j * ldb, mc * sizeof(FLOAT));
        s3_gemm(0, 0, mc, t, t, alpha, work, mc, f, t, 0, b + i0, ldb);
      }
    }
    return;
  }
  const BLASLONG t1 = s3_split(t), t2 = t - t1;
  FLOAT *a22 = a + t1 + t1 * lda;
  FLOAT *a21 = S3_OP(a, lda, tr, t1, 0), *a12 = S3_OP(a, lda, tr, 0, t1);
  if (!side) {
    FLOAT *b2 = b + t1;
    if (lower_op) { /* B2 = A22 B2 + A21 B1, then B1 = A11 B1 */
      s3_trmm_rec(0, 1, tr, unit, t2, n, alpha, a22, lda, b2, ldb, work);
      s3_gemm(tr, 0, t2, n, t1, alpha, a21, lda, b, ldb, 1, b2, ldb);
      s3_trmm_rec(0, 1, tr, unit, t1, n, alpha, a, lda, b, ldb, work);
    } else { /* B1 = A11 B1 + A12 B2, then B2 = A22 B2 */
      s3_trmm_rec(0, 0, tr, unit, t1, n, alpha, a, lda, b, ldb, work);
      s3_gemm(tr, 0, t1, n, t2, alpha, a12, lda, b2, ldb, 1, b, ldb);
      s3_trmm_rec(0, 0, tr, unit, t2, n, alpha, a22, lda, b2, ldb, work);
    }
  } else {
    FLOAT *b2 = b + t1 * ldb;
    if (lower_op) { /* B1 = B1 A11 + B2 A21, then B2 = B2 A22 */
      s3_trmm_rec(1, 1, tr, unit, m, t1, alpha, a, lda, b, ldb, work);
      s3_gemm(0, tr, m, t1, t2, alpha, b2, ldb, a21, lda, 1, b, ldb);
      s3_trmm_rec(1, 1, tr, unit, m, t2, alpha, a22, lda, b2, ldb, work);
    } else { /* B2 = B2 A22 + B1 A12, then B1 = B1 A11 */
      s3_trmm_rec(1, 0, tr, unit, m, t2, alpha, a22, lda, b2, ldb, work);
      s3_gemm(0, tr, m, t2, t1, alpha, b, ldb, a12, lda, 1, b2, ldb);
      s3_trmm_rec(1, 0, tr, unit, m, t1, alpha, a, lda, b, ldb, work);
    }
  }
}

/* ---- TRSM: op(A) X = alpha B (side 0) or X op(A) = alpha B (side 1), X overwrites B ---- */

/*
 * Left-side base solve: x[p][c] = alpha B(row of solve step p, j0 + c) for 16 columns, and back. Rows of B are
 * strided, so blocks of rows x columns are moved with NEON transposes (4 x 4 fp32, 2 x 2 fp64) instead of one
 * strided scalar per element, which cost as much as the solve itself. Step p is row p (fwd) or t - 1 - p.
 */
#ifdef DOUBLE
#define S3_TB 2
#else
#define S3_TB 4
#endif
static inline void s3_tile_io(int load, int fwd, BLASLONG t, BLASLONG nc, FLOAT alpha, FLOAT *b, BLASLONG ldb,
                              FLOAT x[][16]) {
  const BLASLONG cfull = nc / S3_TB * S3_TB, pfull = t / S3_TB * S3_TB;
  for (BLASLONG c0 = 0; c0 < cfull; c0 += S3_TB)
    for (BLASLONG p0 = 0; p0 < pfull; p0 += S3_TB) {
      /* rows rb..rb+TB-1 of B hold steps p0..p0+TB-1 (ascending if fwd, descending otherwise) */
      const BLASLONG rb = fwd ? p0 : t - S3_TB - p0;
      FLOAT *col = b + rb + c0 * ldb;
#ifdef DOUBLE
      if (load) {
        const float64x2_t v0 = vld1q_f64(col), v1 = vld1q_f64(col + ldb), a = vdupq_n_f64(alpha);
        const float64x2_t r0 = vmulq_f64(vzip1q_f64(v0, v1), a), r1 = vmulq_f64(vzip2q_f64(v0, v1), a);
        vst1q_f64(&x[fwd ? p0 : p0 + 1][c0], r0);
        vst1q_f64(&x[fwd ? p0 + 1 : p0][c0], r1);
      } else {
        const float64x2_t r0 = vld1q_f64(&x[fwd ? p0 : p0 + 1][c0]), r1 = vld1q_f64(&x[fwd ? p0 + 1 : p0][c0]);
        vst1q_f64(col, vzip1q_f64(r0, r1));
        vst1q_f64(col + ldb, vzip2q_f64(r0, r1));
      }
#else
      float32x4_t v0, v1, v2, v3;
      if (load) {
        v0 = vld1q_f32(col); v1 = vld1q_f32(col + ldb); v2 = vld1q_f32(col + 2 * ldb); v3 = vld1q_f32(col + 3 * ldb);
      } else { /* rows of the tile, in B's row order */
        v0 = vld1q_f32(&x[fwd ? p0 : p0 + 3][c0]); v1 = vld1q_f32(&x[fwd ? p0 + 1 : p0 + 2][c0]);
        v2 = vld1q_f32(&x[fwd ? p0 + 2 : p0 + 1][c0]); v3 = vld1q_f32(&x[fwd ? p0 + 3 : p0][c0]);
      }
      const float32x4x2_t t01 = vtrnq_f32(v0, v1), t23 = vtrnq_f32(v2, v3);
      float32x4_t r0 = vcombine_f32(vget_low_f32(t01.val[0]), vget_low_f32(t23.val[0]));
      float32x4_t r1 = vcombine_f32(vget_low_f32(t01.val[1]), vget_low_f32(t23.val[1]));
      float32x4_t r2 = vcombine_f32(vget_high_f32(t01.val[0]), vget_high_f32(t23.val[0]));
      float32x4_t r3 = vcombine_f32(vget_high_f32(t01.val[1]), vget_high_f32(t23.val[1]));
      if (load) {
        r0 = vmulq_n_f32(r0, alpha); r1 = vmulq_n_f32(r1, alpha); r2 = vmulq_n_f32(r2, alpha); r3 = vmulq_n_f32(r3, alpha);
        vst1q_f32(&x[fwd ? p0 : p0 + 3][c0], r0); vst1q_f32(&x[fwd ? p0 + 1 : p0 + 2][c0], r1);
        vst1q_f32(&x[fwd ? p0 + 2 : p0 + 1][c0], r2); vst1q_f32(&x[fwd ? p0 + 3 : p0][c0], r3);
      } else {
        vst1q_f32(col, r0); vst1q_f32(col + ldb, r1); vst1q_f32(col + 2 * ldb, r2); vst1q_f32(col + 3 * ldb, r3);
      }
#endif
    }
  /* the rest element by element: columns beyond the last full group, and steps beyond the last full group */
  for (BLASLONG c = 0; c < 16; ++c)
    for (BLASLONG p = (c < cfull ? pfull : 0); p < t; ++p) {
      FLOAT *e = b + (fwd ? p : t - 1 - p) + c * ldb;
      if (load) x[p][c] = c < nc ? alpha * *e : 0;
      else if (c < nc) *e = x[p][c];
    }
}

/*
 * Substitution on a block of at most S3_NB_TRSM, in dot form: each result accumulates in registers over the
 * results solved before it, 16 columns (left side) or 16 rows (right side) at a time, so that each vector FMA
 * needs one load. Indices run in solve order; op(A) is first copied in that order. The diagonal is applied as a
 * reciprocal, as the OpenBLAS TRSM kernels do; alpha comes first.
 */
static inline void s3_trsm_base(int side, int lower_op, int tr, int unit, BLASLONG m, BLASLONG n, FLOAT alpha,
                                FLOAT *a, BLASLONG lda, FLOAT *b, BLASLONG ldb) {
  const BLASLONG t = side ? n : m;
  FLOAT rinv[S3_NB_TRSM], ta[S3_NB_TRSM * S3_NB_TRSM];
  /* solve order: left lower and right upper forward, the others backward */
  const int fwd = side ? !lower_op : lower_op;
#define S3_ORD(x) (fwd ? (x) : t - 1 - (x))
  for (BLASLONG p = 0; p < t; ++p) {
    const BLASLONG q = S3_ORD(p);
    rinv[p] = unit ? 1 : 1 / *S3_OP(a, lda, tr, q, q);
    for (BLASLONG r = 0; r < p; ++r) /* coefficient of solved result r in result p */
      ta[p * t + r] = side ? *S3_OP(a, lda, tr, S3_ORD(r), q) : *S3_OP(a, lda, tr, q, S3_ORD(r));
  }
  if (!side) {
    for (BLASLONG j0 = 0; j0 < n; j0 += 16) { /* 16 columns of B at a time, solved in a local tile */
      const BLASLONG nc = n - j0 < 16 ? n - j0 : 16;
      FLOAT x[S3_NB_TRSM][16];
      s3_tile_io(1, fwd, t, nc, alpha, b + j0 * ldb, ldb, x);
      for (BLASLONG p = 0; p < t; ++p) {
        FLOAT acc[16];
        for (int c = 0; c < 16; ++c) acc[c] = x[p][c];
        for (BLASLONG r = 0; r < p; ++r) {
          const FLOAT l = ta[p * t + r];
          for (int c = 0; c < 16; ++c) acc[c] -= l * x[r][c];
        }
        for (int c = 0; c < 16; ++c) x[p][c] = acc[c] * rinv[p];
      }
      s3_tile_io(0, fwd, t, nc, alpha, b + j0 * ldb, ldb, x);
    }
  } else {
    for (BLASLONG i0 = 0; i0 < m; i0 += 16) { /* 16 rows of B at a time, solved in a local tile */
      const BLASLONG mr = m - i0 < 16 ? m - i0 : 16;
      FLOAT x[S3_NB_TRSM][16];
      for (BLASLONG p = 0; p < t; ++p) {
        const FLOAT *bj = b + i0 + S3_ORD(p) * ldb;
        FLOAT acc[16];
        if (mr == 16)
          for (int r = 0; r < 16; ++r) acc[r] = alpha * bj[r];
        else
          for (int r = 0; r < 16; ++r) acc[r] = r < mr ? alpha * bj[r] : 0;
        for (BLASLONG q = 0; q < p; ++q) {
          const FLOAT c = ta[p * t + q];
          for (int r = 0; r < 16; ++r) acc[r] -= c * x[q][r];
        }
        for (int r = 0; r < 16; ++r) x[p][r] = acc[r] * rinv[p];
      }
      for (BLASLONG p = 0; p < t; ++p) {
        FLOAT *bj = b + i0 + S3_ORD(p) * ldb;
        for (BLASLONG r = 0; r < mr; ++r) bj[r] = x[p][r];
      }
    }
  }
#undef S3_ORD
}

static inline void s3_trsm_rec(int side, int lower_op, int tr, int unit, BLASLONG m, BLASLONG n, FLOAT alpha,
                               FLOAT *a, BLASLONG lda, FLOAT *b, BLASLONG ldb) {
  const BLASLONG t = side ? n : m;
  if (t <= S3_NB_TRSM) {
    s3_trsm_base(side, lower_op, tr, unit, m, n, alpha, a, lda, b, ldb);
    return;
  }
  const BLASLONG t1 = (t / 2 + S3_NB_TRSM - 1) / S3_NB_TRSM * S3_NB_TRSM, t2 = t - t1;
  FLOAT *a22 = a + t1 + t1 * lda;
  FLOAT *a21 = S3_OP(a, lda, tr, t1, 0), *a12 = S3_OP(a, lda, tr, 0, t1);
  if (!side) {
    FLOAT *b2 = b + t1;
    if (lower_op) { /* X1 = A11 \ alpha B1; B2 = alpha B2 - A21 X1; X2 = A22 \ B2 */
      s3_trsm_rec(0, 1, tr, unit, t1, n, alpha, a, lda, b, ldb);
      s3_gemm(tr, 0, t2, n, t1, -1, a21, lda, b, ldb, alpha, b2, ldb);
      s3_trsm_rec(0, 1, tr, unit, t2, n, 1, a22, lda, b2, ldb);
    } else { /* X2 = A22 \ alpha B2; B1 = alpha B1 - A12 X2; X1 = A11 \ B1 */
      s3_trsm_rec(0, 0, tr, unit, t2, n, alpha, a22, lda, b2, ldb);
      s3_gemm(tr, 0, t1, n, t2, -1, a12, lda, b2, ldb, alpha, b, ldb);
      s3_trsm_rec(0, 0, tr, unit, t1, n, 1, a, lda, b, ldb);
    }
  } else {
    FLOAT *b2 = b + t1 * ldb;
    if (lower_op) { /* X2 = alpha B2 / A22; B1 = alpha B1 - X2 A21; X1 = B1 / A11 */
      s3_trsm_rec(1, 1, tr, unit, m, t2, alpha, a22, lda, b2, ldb);
      s3_gemm(0, tr, m, t1, t2, -1, b2, ldb, a21, lda, alpha, b, ldb);
      s3_trsm_rec(1, 1, tr, unit, m, t1, 1, a, lda, b, ldb);
    } else { /* X1 = alpha B1 / A11; B2 = alpha B2 - X1 A12; X2 = B2 / A22 */
      s3_trsm_rec(1, 0, tr, unit, m, t1, alpha, a, lda, b, ldb);
      s3_gemm(0, tr, m, t2, t1, -1, b, ldb, a12, lda, alpha, b2, ldb);
      s3_trsm_rec(1, 0, tr, unit, m, t2, 1, a22, lda, b2, ldb);
    }
  }
}

/* ---- entry points for the interface files (arguments as they build them); 1: the call was handled here ---- */

static inline int s3_symm_hook(int side, int uplo, blas_arg_t *args) {
  if (!s3_enabled(side ? (double)args->m * args->n * args->n : (double)args->m * args->m * args->n)) return 0;
  if (!side)
    s3_symm_rec(0, uplo, args->m, args->n, *(FLOAT *)args->alpha, (FLOAT *)args->a, args->lda, (FLOAT *)args->b,
                args->ldb, *(FLOAT *)args->beta, (FLOAT *)args->c, args->ldc);
  else /* for the right side args->a is the general matrix and args->b the symmetric one */
    s3_symm_rec(1, uplo, args->m, args->n, *(FLOAT *)args->alpha, (FLOAT *)args->b, args->ldb, (FLOAT *)args->a,
                args->lda, *(FLOAT *)args->beta, (FLOAT *)args->c, args->ldc);
  return 1;
}

/* SYRK (two 0) or SYR2K (two 1) */
static inline int s3_syrk_hook(int two, int uplo, int trans, blas_arg_t *args) {
  if (args->k <= 0 || !s3_enabled((double)args->n * args->n * args->k * (two ? 1.0 : 0.5))) return 0;
  s3_syrk_rec(two, uplo, trans & 1, args->n, args->k, *(FLOAT *)args->alpha, (FLOAT *)args->a, args->lda,
              two ? (FLOAT *)args->b : (FLOAT *)args->a, two ? args->ldb : args->lda, *(FLOAT *)args->beta,
              (FLOAT *)args->c, args->ldc);
  return 1;
}

/* TRSM, or TRMM when interface/trsm.c is compiled with TRMM; alpha is in args->beta as trsm.c stores it */
static inline int s3_trxm_hook(int side, int uplo, int trans, int unit, blas_arg_t *args) {
  if ((side ? args->n : args->m) <= S3_NB_TRSM ||
      !s3_enabled(side ? (double)args->m * args->n * args->n : (double)args->m * args->m * args->n))
    return 0;
  const int lower_op = (uplo == 1) != (trans & 1);
#ifndef TRMM
  s3_trsm_rec(side, lower_op, trans & 1, unit == 0, args->m, args->n, *(FLOAT *)args->beta, (FLOAT *)args->a,
              args->lda, (FLOAT *)args->b, args->ldb);
#else
  FLOAT *work = (FLOAT *)blas_memory_alloc(0);
  if (!work) return 0;
  s3_trmm_rec(side, lower_op, trans & 1, unit == 0, args->m, args->n, *(FLOAT *)args->beta, (FLOAT *)args->a,
              args->lda, (FLOAT *)args->b, args->ldb, work);
  blas_memory_free(work);
#endif
  return 1;
}

#endif
