/***************************************************************************
Copyright (c) 2025, The OpenBLAS Project
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

/* Packed int8 GEMM (OpenBLAS extension, declared in cblas.h).
 *
 *   C(m x n) = op(A)(m x k) * op(B)(k x n)
 *
 * with every operand packed four int8 values per int32 word, value 0 in the
 * low byte (the layout dp4a consumes). K4 is the number of packed words per
 * row, so K = 4 * K4. See cblas.h for the full contract; this file only
 * implements it.
 *
 * Both CBLAS orders are handled by gathering op(A) and op(B) into the row-major
 * layout the HIP kernels consume, rather than by reinterpreting the caller's
 * leading dimensions. That keeps one meaning for lda/ldb whatever the order:
 * the distance between consecutive entries of a stored row in a row-major call
 * and of a stored column in a column-major one, which is what cblas.h
 * documents.
 *
 * The kernels compute
 *
 *     Ck[M x N] = Ak[M x K4] * Btk[N x K4]^T
 *
 * so the problem above maps onto them as M = n, N = m, Ak = op(B)^T and
 * Btk = op(A)^T. The result needs a temporary unless ldc == m, because the
 * kernel writes N x M row-major, which is m x n with stride m in either CBLAS
 * order only when ldc == m; otherwise it is written to a scratch buffer and
 * copied into place with the caller's stride.
 *
 * Both the gathered operands and the scratch result are per-call host
 * allocations, freed before returning, so the entry points are reentrant and
 * thread-safe like the rest of the CBLAS surface.
 *
 * Dispatch: in a build with USE_HIP=1 the tuned kernel runs the product; in any
 * other build, when OPENBLAS_RDNA2_GPU=0, or when the
 * kernel reports an error (no usable device, for instance), the portable CPU
 * path below computes the same product. The CPU path accumulates in 32-bit
 * two's-complement wraparound, bit-identical to v_dot4_i32_i8, so a GPU/CPU
 * comparison is exact rather than tolerance-based.
 */

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"

/* Status codes documented in cblas.h. */
#define GEMM_I8_OK        0
#define GEMM_I8_ERR_ARG  (-1)
#define GEMM_I8_ERR_ALLOC (-2)

#if defined(OPENBLAS_USE_HIP)

extern int rdna2_gemm_i8_host_tuned(int *C, const int *A, const int *Bt,
                                    int M, int N, int K4);
extern int rdna2_gemm_i8_scale_host_tuned(float *C, const int *A, const int *Bt,
                                          int M, int N, int K4,
                                          const float *sa, const float *sb,
                                          int zp_b);

/* The GPU path is on by default, matching the level-3 f32 bridge; the same
 * variable switches both off, and it is read on every call so one process can
 * compare the two paths. */
static int i8_gpu_enabled(void) {
  const char *e = getenv("OPENBLAS_RDNA2_GPU");

  return !(e && (e[0] == '0' || e[0] == 'n' || e[0] == 'N'));
}

#endif

/* Normalise a CBLAS transpose code to 0/1. Conjugation is a no-op on 8-bit
 * integers, so the conjugate codes are accepted as their plain equivalents.
 * Returns 0 when the code is not one of the four CBLAS values. */
static int i8_trans(enum CBLAS_TRANSPOSE t, int *ok) {
  *ok = 1;

  switch (t) {
    case CblasNoTrans:     return 0;
    case CblasTrans:       return 1;
    case CblasConjTrans:   return 1;
    case CblasConjNoTrans: return 0;
    default:               *ok = 0; return 0;
  }
}

/* Signed byte dot product of one packed word pair, wrapped to 32 bits.
 * Four multiplies, exactly the arithmetic v_dot4_i32_i8 performs. */
static uint32_t i8_dot4(uint32_t a, uint32_t b) {
  uint32_t acc = 0;
  int k;

  for (k = 0; k < 4; ++k) {
    const int x = (int)(int8_t)(unsigned char)((a >> (8 * k)) & 0xffu);
    const int y = (int)(int8_t)(unsigned char)((b >> (8 * k)) & 0xffu);

    acc += (uint32_t)(x * y);
  }

  return acc;
}

/* Sum of the four signed bytes of a packed word (row sums for the zero point). */
static int i8_sum4(uint32_t w) {
  int s = 0;
  int k;

  for (k = 0; k < 4; ++k)
    s += (int)(int8_t)(unsigned char)((w >> (8 * k)) & 0xffu);

  return s;
}

/* Position of element (r, c) of op(X) in the caller's buffer. ld is a row
 * stride in a row-major call and a column stride in a column-major one. */
static size_t i8_idx(int rm, int trans, blasint ld, blasint r, blasint c) {
  if (rm) return trans ? (size_t)c * (size_t)ld + (size_t)r
                       : (size_t)r * (size_t)ld + (size_t)c;
  return trans ? (size_t)c + (size_t)r * (size_t)ld
               : (size_t)r + (size_t)c * (size_t)ld;
}

/* C element store, honouring the caller's order. */
static size_t i8_cidx(blasint i, blasint j, blasint ldc, int rm) {
  return rm ? (size_t)i * (size_t)ldc + (size_t)j
            : (size_t)i + (size_t)j * (size_t)ldc;
}

/* Gather both operands into the kernels' layout. Returns 0 on allocation
 * failure, 1 otherwise. pb receives op(B)^T (n rows), pa op(A)^T (m rows). */
static int i8_gather(blasint m, blasint n, blasint k4, int ta, int tb, int rm,
                     const int *a, blasint lda, const int *b, blasint ldb,
                     int **pa_out, int **pb_out) {
  int *pb = (int *)malloc((size_t)n * (size_t)k4 * sizeof(int));
  int *pa = (int *)malloc((size_t)m * (size_t)k4 * sizeof(int));

  if (pb == NULL || pa == NULL) {
    free(pa);
    free(pb);
    return 0;
  }

  /* The kernels want op(A) and op(B)^T as row-major k4-wide buffers: pa holds
   * op(A)(i, p) and pb holds op(B)(p, j), i.e. row j of pb is column j of
   * op(B). */
  {
    blasint i, p;

    for (i = 0; i < m; ++i)
      for (p = 0; p < k4; ++p)
        pa[(size_t)i * (size_t)k4 + (size_t)p] = a[i8_idx(rm, ta, lda, i, p)];

    for (i = 0; i < n; ++i)
      for (p = 0; p < k4; ++p)
        pb[(size_t)i * (size_t)k4 + (size_t)p] = b[i8_idx(rm, tb, ldb, p, i)];
  }

  *pa_out = pa;
  *pb_out = pb;
  return 1;
}

/* C = op(A) * op(B), no scaling. Arguments are already validated. */
static int i8_gemm_do(blasint m, blasint n, blasint k4, int ta, int tb, int rm,
                      const int *a, blasint lda, const int *b, blasint ldb,
                      int *c, blasint ldc) {
  int *pa = NULL, *pb = NULL;
  int rc = GEMM_I8_OK;
  int done = 0;

  if (!i8_gather(m, n, k4, ta, tb, rm, a, lda, b, ldb, &pa, &pb))
    return GEMM_I8_ERR_ALLOC;

#if defined(OPENBLAS_USE_HIP)
  if (i8_gpu_enabled()) {
    int *tmp = NULL;
    int *dst = c;
    int ready = 1;

    if (ldc != m) {
      tmp = (int *)malloc((size_t)m * (size_t)n * sizeof(int));
      if (tmp == NULL) ready = 0;  /* fall through to the CPU path */
      else dst = tmp;
    }

    if (ready && rdna2_gemm_i8_host_tuned(dst, pb, pa, (int)n, (int)m, (int)k4) == 0) {
      if (tmp != NULL) {
        blasint i, j;

        for (j = 0; j < n; ++j)
          for (i = 0; i < m; ++i)
            c[i8_cidx(i, j, ldc, rm)] = tmp[(size_t)j * (size_t)m + (size_t)i];
      }
      done = 1;
    }
    free(tmp);
  }
#endif

  if (!done) {
    blasint i, j, p;

    for (j = 0; j < n; ++j) {
      const int *brow = pb + (size_t)j * (size_t)k4;

      for (i = 0; i < m; ++i) {
        const int *arow = pa + (size_t)i * (size_t)k4;
        uint32_t acc = 0;
        int32_t out;

        for (p = 0; p < k4; ++p)
          acc += i8_dot4((uint32_t)arow[p], (uint32_t)brow[p]);

        memcpy(&out, &acc, sizeof(out));
        c[i8_cidx(i, j, ldc, rm)] = (int)out;
      }
    }

    rc = GEMM_I8_OK;
  }

  free(pa);
  free(pb);
  return rc;
}

/* C = (op(A) * op(B) - zp_b * rowsum(op(B))) * sa * sb. Arguments validated. */
static int i8_gemm_scale_do(blasint m, blasint n, blasint k4, int ta, int tb,
                            int rm, const int *a, blasint lda, const int *b,
                            blasint ldb, float *c, blasint ldc, const float *sa,
                            const float *sb, int zp_b) {
  int *pa = NULL, *pb = NULL;
  int rc = GEMM_I8_OK;
  int done = 0;

  if (!i8_gather(m, n, k4, ta, tb, rm, a, lda, b, ldb, &pa, &pb))
    return GEMM_I8_ERR_ALLOC;

#if defined(OPENBLAS_USE_HIP)
  if (i8_gpu_enabled()) {
    float *tmp = NULL;
    float *dst = c;
    int ready = 1;

    if (ldc != m) {
      tmp = (float *)malloc((size_t)m * (size_t)n * sizeof(float));
      if (tmp == NULL) ready = 0;  /* fall through to the CPU path */
      else dst = tmp;
    }

    /* The kernel applies sa to the rows of its A operand, which is op(B)^T,
     * i.e. the columns of op(B): the API's sb. Its sb indexes the rows of its
     * Bt operand, which is op(A)^T, i.e. the rows of op(A): the API's sa.
     * Hence the swap. */
    if (ready && rdna2_gemm_i8_scale_host_tuned(dst, pb, pa, (int)n, (int)m,
                                                (int)k4, sb, sa, zp_b) == 0) {
      if (tmp != NULL) {
        blasint i, j;

        for (j = 0; j < n; ++j)
          for (i = 0; i < m; ++i)
            c[i8_cidx(i, j, ldc, rm)] =
                tmp[(size_t)j * (size_t)m + (size_t)i];
      }
      done = 1;
    }
    free(tmp);
  }
#endif

  if (!done) {
    blasint i, j, p;

    for (j = 0; j < n; ++j) {
      const int *brow = pb + (size_t)j * (size_t)k4;
      int64_t rowsum = 0;

      if (zp_b != 0) {
        for (p = 0; p < k4; ++p) rowsum += i8_sum4((uint32_t)brow[p]);
      }

      for (i = 0; i < m; ++i) {
        const int *arow = pa + (size_t)i * (size_t)k4;
        uint32_t acc = 0;
        int32_t acc_s;
        int64_t corr;
        float v;

        for (p = 0; p < k4; ++p)
          acc += i8_dot4((uint32_t)arow[p], (uint32_t)brow[p]);

        memcpy(&acc_s, &acc, sizeof(acc_s));
        corr = (int64_t)acc_s - (int64_t)zp_b * rowsum;
        v = (float)corr;
        if (sa != NULL) v *= sa[i];
        if (sb != NULL) v *= sb[j];

        c[i8_cidx(i, j, ldc, rm)] = v;
      }
    }

    rc = GEMM_I8_OK;
  }

  free(pa);
  free(pb);
  return rc;
}

/* Validate one entry point's dimensions and leading dimensions. */
static int i8_valid(blasint m, blasint n, blasint k4, int ta, int tb, int rm,
                    const int *a, blasint lda, const int *b, blasint ldb,
                    const void *c, blasint ldc) {
  blasint minlda, minldb, minldc;

  if (m < 0 || n < 0 || k4 < 0) return 0;
  if (m > INT_MAX || n > INT_MAX || k4 > INT_MAX) return 0;

  if (rm) {
    /* Row-major: ld is the row stride, so it must cover the stored columns. */
    minlda = ta ? m : k4;
    minldb = tb ? k4 : n;
    minldc = n;
  } else {
    /* Column-major: ld is the stored row count. */
    minlda = ta ? k4 : m;
    minldb = tb ? n : k4;
    minldc = m;
  }

  if (minlda < 1) minlda = 1;
  if (minldb < 1) minldb = 1;
  if (minldc < 1) minldc = 1;

  if (lda < minlda || ldb < minldb || ldc < minldc) return 0;

  return 1;
}

int cblas_gemm_i8(enum CBLAS_ORDER Order,
                  enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB,
                  blasint M, blasint N, blasint K4,
                  OPENBLAS_CONST int *A, blasint lda,
                  OPENBLAS_CONST int *B, blasint ldb,
                  int *C, blasint ldc) {
  int ta, tb, ok, rm;
  blasint i, j;

  ta = i8_trans(TransA, &ok);
  if (!ok) return GEMM_I8_ERR_ARG;
  tb = i8_trans(TransB, &ok);
  if (!ok) return GEMM_I8_ERR_ARG;

  if (Order == CblasRowMajor)      rm = 1;
  else if (Order == CblasColMajor) rm = 0;
  else return GEMM_I8_ERR_ARG;

  if (!i8_valid(M, N, K4, ta, tb, rm, A, lda, B, ldb, C, ldc))
    return GEMM_I8_ERR_ARG;

  /* Quick return before the data pointers are required: no element is
   * referenced when either extent is zero, as in BLAS. */
  if (M == 0 || N == 0) return GEMM_I8_OK;
  if (C == NULL) return GEMM_I8_ERR_ARG;
  if (K4 > 0 && (A == NULL || B == NULL)) return GEMM_I8_ERR_ARG;

  /* Empty product: C = 0, the only value a k = 0 product can have. */
  if (K4 == 0) {
    for (j = 0; j < N; ++j)
      for (i = 0; i < M; ++i) C[i8_cidx(i, j, ldc, rm)] = 0;
    return GEMM_I8_OK;
  }

  return i8_gemm_do(M, N, K4, ta, tb, rm, A, lda, B, ldb, C, ldc);
}

int cblas_gemm_i8_scale(enum CBLAS_ORDER Order,
                        enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB,
                        blasint M, blasint N, blasint K4,
                        OPENBLAS_CONST int *A, blasint lda,
                        OPENBLAS_CONST int *B, blasint ldb,
                        float *C, blasint ldc,
                        OPENBLAS_CONST float *sa, OPENBLAS_CONST float *sb,
                        int zp_b) {
  int ta, tb, ok, rm;
  blasint i, j;

  ta = i8_trans(TransA, &ok);
  if (!ok) return GEMM_I8_ERR_ARG;
  tb = i8_trans(TransB, &ok);
  if (!ok) return GEMM_I8_ERR_ARG;

  if (Order == CblasRowMajor)      rm = 1;
  else if (Order == CblasColMajor) rm = 0;
  else return GEMM_I8_ERR_ARG;

  if (!i8_valid(M, N, K4, ta, tb, rm, A, lda, B, ldb, C, ldc))
    return GEMM_I8_ERR_ARG;

  if (M == 0 || N == 0) return GEMM_I8_OK;
  if (C == NULL) return GEMM_I8_ERR_ARG;
  if (K4 > 0 && (A == NULL || B == NULL)) return GEMM_I8_ERR_ARG;

  if (K4 == 0) {
    for (j = 0; j < N; ++j)
      for (i = 0; i < M; ++i) C[i8_cidx(i, j, ldc, rm)] = 0.0f;
    return GEMM_I8_OK;
  }

  return i8_gemm_scale_do(M, N, K4, ta, tb, rm, A, lda, B, ldb, C, ldc,
                          sa, sb, zp_b);
}
