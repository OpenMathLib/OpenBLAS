/* Must match common.h: on Windows 64-bit OpenBLAS uses long long. */
#if defined(_WIN64)
typedef long long BLASLONG;
#else
typedef long BLASLONG;
#endif
typedef float FLOAT;

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#ifdef __cplusplus
extern "C" {
#endif

int rdna2_sgemm_f32_host_tuned(float *C, const float *A, const float *Bt,
                               int M, int N, int K);

#ifdef __cplusplus
}
#endif

/* Bridges between OpenBLAS' level-3 entry points and the hand-written
 * RDNA2/RDNA4 HIP kernels in gpu/rdna2/rdna2_kernels.hip.
 *
 * The HIP kernel computes, in row-major terms,
 *
 *     C[i][j] = sum_p A[i][p] * Bt[j][p]        A is M x K, Bt is N x K
 *
 * so its row-major M x N result, read with row stride N, is exactly a
 * column-major N x M matrix. Every routine below calls the kernel with
 * (M, N) = (n, m), which makes the kernel's output the column-major m x n
 * result OpenBLAS expects. Column-major B (k x n, ldb == k) is already the
 * kernel's N x K left operand; a column-major A (m x k) is contiguous in m and
 * is transposed once into a temporary.
 *
 * All four bridges funnel into rdna2_gemm_cm_core(), so every level-3 routine
 * runs the same GEMM path that the standalone harness verifies against the CPU.
 */

/* Below this the CPU path wins, because the host entry point allocates, copies
 * and (for column-major A) transposes on every call. Measured single-call
 * crossover is around 1500^3; see gpu/rdna2/bench_full.txt. */
#define RDNA2_GPU_MIN_MNK 2000000000.0

/* Largest dimension the host launchers accept. */
#define RDNA2_MAX_DIM 65535

/* Block size of the blocked triangular solve. */
#define RDNA2_TRSM_NB 128

/* Read on every call (not cached) so a process can compare both paths, and so a
 * long-running caller can switch the GPU off without restarting. */
static int rdna2_gpu_enabled(void) {
  const char *e = getenv("OPENBLAS_RDNA2_GPU");
  return !(e && (e[0] == '0' || e[0] == 'n' || e[0] == 'N'));
}

static double rdna2_min_mnk(void) {
  const char *e = getenv("OPENBLAS_RDNA2_GPU_MIN_MNK");
  return (e && e[0]) ? atof(e) : RDNA2_GPU_MIN_MNK;
}

static int rdna2_debug(void) {
  const char *e = getenv("OPENBLAS_RDNA2_GPU_DEBUG");
  return e && e[0] == '1';
}

static long rdna2_ms(clock_t since) {
  return (long)((clock() - since) * 1000.0 / (double)CLOCKS_PER_SEC);
}

/* Column-major -> row-major copy, 64x64 tiled:
 *
 *     dst[i*ldd + p] = src[i + p*lds]      i in [0,rows), p in [0,cols)
 *
 * With ldd == cols this is a plain transpose. The naive version stores with a
 * stride of cols*4 bytes, one cache line per element; at 2048^3 that alone cost
 * more than the kernel. Tiled, both the store block (64 rows of 256 B) and the
 * matching load block stay in cache.
 */
static void rdna2_cm_to_rm(const float *src, BLASLONG lds, float *dst,
                           BLASLONG ldd, BLASLONG rows, BLASLONG cols) {
  const BLASLONG T = 64;
  BLASLONG i0, p0;

  for (p0 = 0; p0 < cols; p0 += T) {
    const BLASLONG p1 = (p0 + T < cols) ? p0 + T : cols;

    for (i0 = 0; i0 < rows; i0 += T) {
      const BLASLONG i1 = (i0 + T < rows) ? i0 + T : rows;
      BLASLONG p, i;

      for (p = p0; p < p1; ++p) {
        const float *s = src + (size_t)p * (size_t)lds + (size_t)i0;
        float *d = dst + (size_t)i0 * (size_t)ldd + (size_t)p;

        for (i = i0; i < i1; ++i) {
          *d = *s++;
          d += ldd;
        }
      }
    }
  }
}

/* C(m x n) = alpha * A(m x k) * B(k x n) + beta * C, every operand column-major
 * with strides lda, ldb and ldc.
 *
 * Returns 1 when C was produced on the GPU, 0 when the caller must fall back to
 * the CPU path. There is deliberately no size gate here: each bridge gates once
 * for the whole problem, and the blocked triangular solve issues sub-block GEMMs
 * that are individually small but must still run.
 */
static int rdna2_gemm_cm_core(BLASLONG m, BLASLONG n, BLASLONG k,
                              double alpha, const float *a, BLASLONG lda,
                              const float *b, BLASLONG ldb,
                              double beta, float *c, BLASLONG ldc) {
  const float *left = b;   /* kernel's N x K (== n x k) row-major operand */
  const float *right = a;  /* kernel's M x K (== m x k) row-major operand */
  float *bt = NULL, *at = NULL, *tmp = NULL;
  int direct, ok;

  if (ldb != k) {
    BLASLONG i, j;

    bt = (float *)malloc((size_t)n * (size_t)k * sizeof(float));
    if (bt == NULL) return 0;
    /* The kernel's left operand is B^T stored row-major n x k, so
     * bt[j*k + p] = B(p, j) = b[p + j*ldb]. That is each column of B (its first
     * k entries) laid out as a contiguous row, i.e. a plain copy per column and
     * not the faithful column-major -> row-major copy used for the other
     * operand. When ldb == k the same bytes are produced, which is why B can be
     * passed through untouched in that case. */
    for (j = 0; j < n; ++j) {
      const float *s = b + (size_t)j * (size_t)ldb;
      float *d = bt + (size_t)j * (size_t)k;

      for (i = 0; i < k; ++i) d[i] = s[i];
    }
    left = bt;
  }

  /* The kernel's second operand must be A stored row-major with row stride k.
   * A column-major m x k matrix read with stride k is its transpose, not
   * itself, so this copy is required even when lda == k. (Only the left operand
   * gets that shortcut: a column-major k x n matrix with ldb == k already is the
   * row-major n x k matrix the kernel wants.) */
  at = (float *)malloc((size_t)m * (size_t)k * sizeof(float));
  if (at == NULL) {
    free(bt);
    return 0;
  }
  rdna2_cm_to_rm(a, lda, at, k, m, k);
  right = at;

  /* The kernel writes its n x m row-major result straight into a column-major
   * m x n matrix only when ldc == m and the output aliases neither input.
   * TRMM and TRSM update b in place, so they always take the temporary. */
  direct = (ldc == m) && (c != b) && (c != a) && (alpha == 1.0) && (beta == 0.0);

  if (!direct) {
    tmp = (float *)malloc((size_t)m * (size_t)n * sizeof(float));
    if (tmp == NULL) {
      free(bt);
      free(at);
      return 0;
    }
  }

  ok = rdna2_sgemm_f32_host_tuned(direct ? c : tmp, left, right,
                                  (int)n, (int)m, (int)k) == 0;

  if (ok && !direct) {
    BLASLONG i, j;

    for (j = 0; j < n; ++j) {
      for (i = 0; i < m; ++i) {
        const float d = tmp[(size_t)j * (size_t)m + (size_t)i];
        float *cp = c + (size_t)i + (size_t)j * (size_t)ldc;

        if (alpha == 1.0 && beta == 0.0)
          *cp = d;
        else if (beta == 0.0)
          *cp = (float)(alpha * (double)d);
        else
          *cp = (float)(alpha * (double)d + beta * (double)*cp);
      }
    }
  }

  free(bt);
  free(at);
  free(tmp);

  return ok;
}

/* Returns 1 when the product was computed on the GPU, 0 when the caller should
 * carry on with the generic CPU path. */
int rdna2_oblas_sgemm_cm(BLASLONG m, BLASLONG n, BLASLONG k, double alpha,
                         const float *a, BLASLONG lda, const float *b,
                         BLASLONG ldb, double beta, float *c, BLASLONG ldc,
                         int transa, int transb) {
  int handled;
  clock_t t0;

  if (!rdna2_gpu_enabled()) return 0;
  /* Only the untransposed layouts are wired up; the core handles any leading
   * dimension, so lda/ldb/ldc are not constrained here. Transposed operands
   * need a different gather and stay on the CPU path. */
  if (transa != 0 || transb != 0) return 0;
  if (m <= 0 || n <= 0 || k <= 0) return 0;
  if (m > RDNA2_MAX_DIM || n > RDNA2_MAX_DIM || k > RDNA2_MAX_DIM) return 0;
  if ((double)m * (double)n * (double)k < rdna2_min_mnk()) return 0;

  t0 = clock();
  handled = rdna2_gemm_cm_core(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

  if (rdna2_debug())
    fprintf(stderr, "[rdna2] sgemm %ldx%ldx%ld handled=%d %ldms\n",
            (long)m, (long)n, (long)k, handled, rdna2_ms(t0));

  return handled;
}

/* SSYMM: C := alpha*A*B + beta*C (side=left) or alpha*B*A + beta*C (side=right),
 * where the m x m (left) or n x n (right) matrix A is symmetric and only the
 * uplo triangle is stored.
 *
 * A is expanded into a dense symmetric buffer, after which the operation is an
 * ordinary GEMM. The mirror image is always read from the stored triangle, so
 * the triangle the caller did not fill in is never touched.
 *
 * side: 0 = left, 1 = right.  uplo: 0 = upper, 1 = lower.
 */
int rdna2_oblas_ssymm_cm(BLASLONG m, BLASLONG n, double alpha,
                         const float *a, BLASLONG lda,
                         const float *b, BLASLONG ldb,
                         double beta, float *c, BLASLONG ldc,
                         int side, int uplo) {
  BLASLONG kdim, i, j;
  float *sym;
  int ok;
  clock_t t0;

  if (!rdna2_gpu_enabled()) return 0;
  if (m <= 0 || n <= 0) return 0;
  if (m > RDNA2_MAX_DIM || n > RDNA2_MAX_DIM) return 0;

  kdim = side ? n : m;
  if (kdim > RDNA2_MAX_DIM) return 0;
  if ((double)m * (double)n * (double)kdim < rdna2_min_mnk()) return 0;
  if (lda < kdim || ldb < m || ldc < m) return 0;

  sym = (float *)malloc((size_t)kdim * (size_t)kdim * sizeof(float));
  if (sym == NULL) return 0;

  for (j = 0; j < kdim; ++j) {
    for (i = 0; i < kdim; ++i) {
      const int stored = (uplo == 0) ? (i <= j) : (i >= j);

      sym[(size_t)i + (size_t)j * (size_t)kdim] =
          stored ? a[(size_t)i + (size_t)j * (size_t)lda]
                 : a[(size_t)j + (size_t)i * (size_t)lda];
    }
  }

  t0 = clock();

  if (!side)
    ok = rdna2_gemm_cm_core(m, n, m, alpha, sym, kdim, b, ldb, beta, c, ldc);
  else
    ok = rdna2_gemm_cm_core(m, n, n, alpha, b, ldb, sym, kdim, beta, c, ldc);

  if (rdna2_debug())
    fprintf(stderr, "[rdna2] ssymm %ldx%ld side=%d uplo=%d handled=%d %ldms\n",
            (long)m, (long)n, side, uplo, ok, rdna2_ms(t0));

  free(sym);
  return ok;
}

/* TRMM: B := alpha*op(A)*B (side=left) or alpha*B*op(A) (side=right), with A
 * triangular and m x m for side=left, n x n for side=right.
 *
 * OpenBLAS' TRMM has no separate output, so B is both the input and the result
 * (column-major, stride ldb). A is expanded into a dense buffer that is zero
 * outside the stored triangle, which makes the product identical to the
 * triangular one while letting it go through the GEMM core. TRMM has no beta:
 * the result always overwrites B.
 *
 * trans & 1 selects op(A) = A^T. unit_diag != 0 means A has an implicit unit
 * diagonal.
 */
int rdna2_oblas_strmm_cm(BLASLONG m, BLASLONG n, double alpha,
                         const float *a, BLASLONG lda,
                         float *b, BLASLONG ldb,
                         int side, int uplo, int trans, int unit_diag) {
  BLASLONG kdim, i, j;
  const int tr = trans & 1;
  float *d;
  int ok;
  clock_t t0;

  if (!rdna2_gpu_enabled()) return 0;
  if (m <= 0 || n <= 0) return 0;
  if (m > RDNA2_MAX_DIM || n > RDNA2_MAX_DIM) return 0;

  kdim = side ? n : m;
  if (kdim > RDNA2_MAX_DIM) return 0;
  if ((double)m * (double)n * (double)kdim < rdna2_min_mnk()) return 0;
  if (lda < kdim || ldb < m) return 0;

  d = (float *)malloc((size_t)kdim * (size_t)kdim * sizeof(float));
  if (d == NULL) return 0;

  for (j = 0; j < kdim; ++j) {
    for (i = 0; i < kdim; ++i) {
      /* Source index in the stored A. */
      const BLASLONG si = tr ? j : i;
      const BLASLONG sj = tr ? i : j;
      const int stored = (uplo == 0) ? (si <= sj) : (si >= sj);
      float v = 0.0f;

      if (stored) v = a[(size_t)si + (size_t)sj * (size_t)lda];
      if (si == sj && unit_diag) v = 1.0f;

      d[(size_t)i + (size_t)j * (size_t)kdim] = v;
    }
  }

  t0 = clock();

  if (!side)
    ok = rdna2_gemm_cm_core(m, n, m, alpha, d, kdim, b, ldb, 0.0, b, ldb);
  else
    ok = rdna2_gemm_cm_core(m, n, n, alpha, b, ldb, d, kdim, 0.0, b, ldb);

  if (rdna2_debug())
    fprintf(stderr, "[rdna2] strmm %ldx%ld side=%d uplo=%d trans=%d unit=%d handled=%d %ldms\n",
            (long)m, (long)n, side, uplo, tr, unit_diag, ok, rdna2_ms(t0));

  free(d);
  return ok;
}

/* Solve T * X = Y in place: T is the dense k x k triangular factor with stride
 * ldt, X is k x n column-major with stride ldx, and Y is already in X.
 * uplo: 0 = upper, 1 = lower. Returns 1 on success.
 *
 * Blocked substitution: the diagonal block is solved on the host, and the
 * trailing rank-nb update is a GEMM, which is what puts the work on the GPU.
 */
static int rdna2_trsm_left(BLASLONG k, BLASLONG n, const float *t, BLASLONG ldt,
                           float *x, BLASLONG ldx, int uplo) {
  if (uplo == 0) {
    /* Upper: back substitution, blocks from the bottom up. */
    BLASLONG p1;

    for (p1 = k; p1 > 0; ) {
      const BLASLONG p0 = (p1 > RDNA2_TRSM_NB) ? p1 - RDNA2_TRSM_NB : 0;
      BLASLONG i, j, q;

      for (j = 0; j < n; ++j) {
        float *xj = x + (size_t)j * (size_t)ldx;

        for (i = p1 - 1; i >= p0; --i) {
          double s = (double)xj[i];

          for (q = i + 1; q < p1; ++q)
            s -= (double)t[(size_t)i + (size_t)q * (size_t)ldt] * (double)xj[q];

          xj[i] = (float)(s / (double)t[(size_t)i + (size_t)i * (size_t)ldt]);
        }
      }

      /* X[0:p0, :] -= T[0:p0, p0:p1] * X[p0:p1, :] */
      if (p0 > 0 &&
          !rdna2_gemm_cm_core(p0, n, p1 - p0, -1.0,
                              t + (size_t)p0 * (size_t)ldt, ldt,
                              x + p0, ldx, 1.0, x, ldx))
        return 0;

      p1 = p0;
    }
  } else {
    /* Lower: forward substitution, blocks from the top down. */
    BLASLONG p0;

    for (p0 = 0; p0 < k; ) {
      const BLASLONG p1 = (p0 + RDNA2_TRSM_NB < k) ? p0 + RDNA2_TRSM_NB : k;
      BLASLONG i, j, q;

      for (j = 0; j < n; ++j) {
        float *xj = x + (size_t)j * (size_t)ldx;

        for (i = p0; i < p1; ++i) {
          double s = (double)xj[i];

          for (q = p0; q < i; ++q)
            s -= (double)t[(size_t)i + (size_t)q * (size_t)ldt] * (double)xj[q];

          xj[i] = (float)(s / (double)t[(size_t)i + (size_t)i * (size_t)ldt]);
        }
      }

      /* X[p1:k, :] -= T[p1:k, p0:p1] * X[p0:p1, :] */
      if (p1 < k &&
          !rdna2_gemm_cm_core(k - p1, n, p1 - p0, -1.0,
                              t + p1 + (size_t)p0 * (size_t)ldt, ldt,
                              x + p0, ldx, 1.0, x + p1, ldx))
        return 0;

      p0 = p1;
    }
  }

  return 1;
}

/* TRSM: solve op(A)*X = alpha*B (side=left) or X*op(A) = alpha*B (side=right)
 * in place in b, with A triangular.
 *
 * side=left solves directly. side=right is solved as the transposed system
 * op(A)^T * X^T = alpha*B^T, so both cases share one left-looking blocked
 * solver; the two transposes of B cost a pass over b each.
 *
 * The solve runs in a scratch copy and b is only written once the whole solve
 * succeeded, so a late GPU failure leaves b intact for the CPU fallback.
 */
int rdna2_oblas_strsm_cm(BLASLONG m, BLASLONG n, double alpha,
                         const float *a, BLASLONG lda,
                         float *b, BLASLONG ldb,
                         int side, int uplo, int trans, int unit_diag) {
  BLASLONG kdim, i, j;
  const int tr = trans & 1;
  int tri;
  float *d = NULL, *x = NULL;
  int ok = 0;
  clock_t t0;

  if (!rdna2_gpu_enabled()) return 0;
  if (m <= 0 || n <= 0) return 0;
  if (m > RDNA2_MAX_DIM || n > RDNA2_MAX_DIM) return 0;

  kdim = side ? n : m;
  if (kdim > RDNA2_MAX_DIM) return 0;
  if ((double)m * (double)n * (double)kdim < rdna2_min_mnk()) return 0;
  if (lda < kdim || ldb < m) return 0;

  /* Triangle of the factor the solver actually uses: op(A) for side=left and
   * op(A)^T for side=right. Transposing flips the triangle, so side=left is
   * !uplo when trans is set, and side=right is its complement. */
  tri = side ? (tr ? uplo : !uplo) : (tr ? !uplo : uplo);

  d = (float *)malloc((size_t)kdim * (size_t)kdim * sizeof(float));
  if (d == NULL) return 0;

  for (j = 0; j < kdim; ++j) {
    for (i = 0; i < kdim; ++i) {
      /* Dense copy of T: op(A) for side=left, op(A)^T for side=right. */
      BLASLONG si, sj;
      int stored;
      float v = 0.0f;

      if (!side) {
        si = tr ? j : i;
        sj = tr ? i : j;
      } else {
        si = tr ? i : j;
        sj = tr ? j : i;
      }

      stored = (uplo == 0) ? (si <= sj) : (si >= sj);
      if (stored) v = a[(size_t)si + (size_t)sj * (size_t)lda];
      if (si == sj && unit_diag) v = 1.0f;

      d[(size_t)i + (size_t)j * (size_t)kdim] = v;
    }
  }

  t0 = clock();

  if (!side) {
    x = (float *)malloc((size_t)m * (size_t)n * sizeof(float));
    if (x != NULL) {
      for (j = 0; j < n; ++j)
        for (i = 0; i < m; ++i)
          x[(size_t)i + (size_t)j * (size_t)m] =
              (float)(alpha * (double)b[(size_t)i + (size_t)j * (size_t)ldb]);

      ok = rdna2_trsm_left(m, n, d, kdim, x, m, tri);

      if (ok)
        for (j = 0; j < n; ++j)
          for (i = 0; i < m; ++i)
            b[(size_t)i + (size_t)j * (size_t)ldb] =
                x[(size_t)i + (size_t)j * (size_t)m];
    }
  } else {
    /* op(A)^T * Y = alpha * B^T, Y is n x m. */
    x = (float *)malloc((size_t)n * (size_t)m * sizeof(float));
    if (x != NULL) {
      for (j = 0; j < m; ++j)
        for (i = 0; i < n; ++i)
          x[(size_t)i + (size_t)j * (size_t)n] =
              (float)(alpha * (double)b[(size_t)j + (size_t)i * (size_t)ldb]);

      ok = rdna2_trsm_left(n, m, d, kdim, x, n, tri);

      if (ok)
        for (j = 0; j < n; ++j)
          for (i = 0; i < m; ++i)
            b[(size_t)i + (size_t)j * (size_t)ldb] =
                x[(size_t)j + (size_t)i * (size_t)n];
    }
  }

  if (rdna2_debug())
    fprintf(stderr, "[rdna2] strsm %ldx%ld side=%d uplo=%d trans=%d unit=%d handled=%d %ldms\n",
            (long)m, (long)n, side, uplo, tr, unit_diag, ok, rdna2_ms(t0));

  free(d);
  free(x);
  return ok;
}
