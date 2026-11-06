/* Checks that the built OpenBLAS library really runs the level-3 GPU bridges.
 *
 * Links against libopenblas_rdna2-*.a as built with `TARGET=GFX1201 USE_HIP=1`
 * (or any build with USE_HIP=1) and checks two things:
 *
 *   SGEMM  calls cblas_sgemm twice on the same inputs: once with the HIP path
 *          enabled, once with OPENBLAS_RDNA2_GPU=0. The CPU result is the
 *          reference, so a mismatch means the GPU path is wrong.
 *
 *   SYRK   calls cblas_ssyrk / cblas_ssyr2k over a range of layouts and compares
 *   SYR2K  the stored triangle against a double precision reference computed
 *          here, and checks the two contracts the bridges must keep: the
 *          triangle the caller did not store is untouched, and so is the ldc
 *          padding. Small odd shapes run with the size gate forced open so the
 *          device sees partial tiles; one case per op runs at the default gate
 *          so the real dispatch decision is covered as well. Without a GPU the
 *          same cases still run, on the CPU kernels, and still check the
 *          triangle contract.
 *
 *   SYMM   cover every order/side/uplo/trans/diag combination against a double
 *   TRMM   precision reference (an LU solve for TRSM, so the reference shares
 *   TRSM   none of the triangle logic under test). The unstored triangle of the
 *          factor, and the stored diagonal of a unit-diagonal one, are filled
 *          with a value that no correct result can produce, so a bridge that
 *          reads them fails loudly; the ldb/ldc padding must be untouched.
 *
 * This build may use NUM_THREADS=1, so the CPU timings here are single-threaded
 * and only meaningful as a sanity check; the CPU-vs-GPU numbers in bench_full.txt
 * use a 32-thread OpenBLAS.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "cblas.h"

static double now_s(void) {
  return (double)clock() / (double)CLOCKS_PER_SEC;
}

#ifdef _WIN32
#define SETENV(k, v) _putenv_s((k), (v))
#else
#define SETENV(k, v) setenv((k), (v), 1)
#endif

/* ------------------------------------------------------------------------- */
/* SYRK / SYR2K                                                              */
/* ------------------------------------------------------------------------- */

/* Stored triangle, in mathematical terms for both orders: uplo 0 = upper
 * (i <= j), 1 = lower (i >= j). The CBLAS row-major entry points map the enum
 * onto the same triangle of the mathematical result, only the storage differs. */
static int sym_stored(int uplo, int i, int j) {
  return (uplo == 0) ? (i <= j) : (i >= j);
}

/* Element (i,j) of a stored factor: op(A) for tr == 0, op(A)^T for tr != 0.
 * Column-major: nrowa x cols with stride ld. Row-major: the row stride is ld. */
static double sym_op(const float *a, int ld, int order, int tr, int i, int p) {
  if (order == CblasColMajor)
    return tr ? (double)a[(size_t)p + (size_t)i * ld]
              : (double)a[(size_t)i + (size_t)p * ld];
  return tr ? (double)a[(size_t)p * ld + (size_t)i]
            : (double)a[(size_t)i * ld + (size_t)p];
}

/* A wrong triangle can be wrong everywhere at once, and one line per element
 * buries the log (a 1300x1300 clobber is 845k lines), so print the first few
 * and count the rest. */
#define SYM_REPORT_MAX 8
static void violation(long *count, const char *tag, int row, int col,
                      const char *what, double got, double want, double rel) {
  if ((*count)++ < SYM_REPORT_MAX) {
    if (what)
      printf("%s: %s element (%d,%d) was written\n", tag, what, row, col);
    else
      printf("%s: (%d,%d) got %.6f want %.6f (rel %.2e)\n", tag, row, col, got,
             want, rel);
  }
}

/* Element (i,j) of an m x n matrix with stride ld, in either store order. */
static double mat_get(const float *x, int ld, int order, int i, int j) {
  return (order == CblasColMajor) ? (double)x[(size_t)i + (size_t)j * ld]
                                  : (double)x[(size_t)i * ld + (size_t)j];
}

static void mat_put(float *x, int ld, int order, int i, int j, float v) {
  if (order == CblasColMajor)
    x[(size_t)i + (size_t)j * ld] = v;
  else
    x[(size_t)i * ld + (size_t)j] = v;
}

/* Bytes to allocate for an m x n matrix with row or column stride ld. */
static size_t mat_size(int ld, int order, int m, int n) {
  return (order == CblasColMajor) ? (size_t)ld * (size_t)n
                                  : (size_t)ld * (size_t)m;
}

/* Element (i,p) of a symmetric factor given only its uplo triangle. The mirror
 * image is read out of the stored triangle, never from the unused one. */
static double sym_el(const float *a, int ld, int order, int uplo, int i, int p) {
  const int stored = sym_stored(uplo, i, p);
  return mat_get(a, ld, order, stored ? i : p, stored ? p : i);
}

/* Element (i,p) of op(T) for a triangular factor T: op(T) is T (tr == 0) or
 * T^T, zero outside the stored triangle, with an implicit unit diagonal when
 * unit is set. */
static double tri_el(const float *a, int ld, int order, int uplo, int tr,
                     int unit, int i, int p) {
  const int si = tr ? p : i;
  const int sj = tr ? i : p;
  if (!sym_stored(uplo, si, sj)) return 0.0;
  if (si == sj && unit) return 1.0;
  return mat_get(a, ld, order, si, sj);
}

/* Runs one SYRK (two == 0) or SYR2K (two == 1) case and compares it against a
 * double precision reference. Returns 0 when the case passed. */
static int sym_case(int two, int order, int uplo, int trans, int n, int k,
                    double alpha, double beta) {
  const int ta = trans ? CblasTrans : CblasNoTrans;
  const CBLAS_UPLO up = uplo ? CblasLower : CblasUpper;
  const int ldc = n + 5;
  const int lda = ((order == CblasColMajor) ? (trans ? k : n) : (trans ? n : k)) + 3;
  const int ldb = lda;
  const int rows = trans ? k : n; /* stored rows of A: k x n when transposed */
  const int cols = trans ? n : k; /* stored columns */
  /* Column-major: lda strides the rows, so the allocation is lda * cols.
   * Row-major: lda is the row stride, so it is lda * rows. */
  const size_t asz = (order == CblasColMajor) ? (size_t)lda * (size_t)cols
                                              : (size_t)lda * (size_t)rows;
  const size_t csz = (size_t)ldc * (size_t)n;
  const float sentinel = -12345.5f;
  float *A = (float *)malloc(asz * sizeof(float));
  float *B = (float *)malloc(asz * sizeof(float));
  float *C = (float *)malloc(csz * sizeof(float));
  float *C0 = (float *)malloc(csz * sizeof(float));
  size_t i;
  int idx, bad = 0;
  long nviol = 0;

  if (!A || !B || !C || !C0) {
    printf("syrk/syr2k: allocation failed at n=%d k=%d\n", n, k);
    free(A);
    free(B);
    free(C);
    free(C0);
    return 1;
  }

  /* Deterministic operands, including the padding rows of every column. */
  for (i = 0; i < asz; ++i) {
    A[i] = (float)((double)((i * 7) % 23) - 11.0) / 9.0f;
    B[i] = (float)((double)((i * 11) % 19) - 9.0) / 7.0f;
  }
  /* C starts at a value that no result can produce, so an element the routine
   * never writes is detectable. */
  for (i = 0; i < csz; ++i) C0[i] = C[i] = (float)(i % 17) * 0.5f + sentinel;

  if (two)
    cblas_ssyr2k(order, up, (CBLAS_TRANSPOSE)ta, n, k, (float)alpha, A, lda, B,
                 ldb, (float)beta, C, ldc);
  else
    cblas_ssyrk(order, up, (CBLAS_TRANSPOSE)ta, n, k, (float)alpha, A, lda,
                (float)beta, C, ldc);

  /* Walk the whole allocation: in-matrix elements of the stored triangle must
   * match the reference, everything else (the other triangle and the ldc
   * padding) must still hold the value written before the call. */
  for (idx = 0; idx < (int)csz; ++idx) {
    const int row = (order == CblasColMajor) ? idx % ldc : idx / ldc;
    const int col = (order == CblasColMajor) ? idx / ldc : idx % ldc;

    if (row >= n || col >= n) {
      if (C[idx] != C0[idx]) {
        violation(&nviol, "syrk/syr2k", row, col, "padding", 0.0, 0.0, 0.0);
        bad = 1;
      }
      continue;
    }
    if (!sym_stored(uplo, row, col)) {
      if (C[idx] != C0[idx]) {
        violation(&nviol, "syrk/syr2k", row, col, "unstored triangle", 0.0, 0.0,
                  0.0);
        bad = 1;
      }
      continue;
    }
    {
      double ref = 0.0, got, err, scale;
      int p;

      for (p = 0; p < k; ++p) {
        const double a_ip = sym_op(A, lda, order, trans, row, p);
        const double a_jp = sym_op(A, lda, order, trans, col, p);

        if (two) {
          const double b_ip = sym_op(B, ldb, order, trans, row, p);
          const double b_jp = sym_op(B, ldb, order, trans, col, p);

          ref += a_ip * b_jp + b_ip * a_jp;
        } else {
          ref += a_ip * a_jp;
        }
      }
      ref = alpha * ref + beta * (double)C0[idx];
      got = (double)C[idx];
      err = fabs(got - ref);
      scale = fabs(ref) > 1.0 ? fabs(ref) : 1.0;
      if (err / scale > 5e-4) {
        violation(&nviol, "syrk/syr2k", row, col, NULL, got, ref, err / scale);
        bad = 1;
      }
    }
  }

  if (nviol > SYM_REPORT_MAX)
    printf("syrk/syr2k: ... %ld more violation(s) in this case\n",
           nviol - SYM_REPORT_MAX);

  if (bad)
    printf("  case: %s order=%d uplo=%d trans=%d n=%d k=%d alpha=%g beta=%g\n",
           two ? "syr2k" : "syrk", order, uplo, trans, n, k, alpha, beta);
  else
    printf("ok   %s order=%d uplo=%d trans=%d n=%d k=%d alpha=%g beta=%g\n",
           two ? "syr2k" : "syrk", order, uplo, trans, n, k, alpha, beta);

  free(A);
  free(B);
  free(C);
  free(C0);
  return bad;
}

/* ------------------------------------------------------------------------- */
/* SYMM / TRMM / TRSM                                                        */
/* ------------------------------------------------------------------------- */

/* Written into every element the operation must not look at (the unstored
 * triangle of a symmetric or triangular factor, a unit diagonal). A pipeline
 * that reads it turns this into errors of this order, far above the tolerances
 * below, without producing the NaNs that would also poison a correct build. */
#define UNUSED_FILL 1.0e6f

/* Everything outside the logical m x n block of a matrix with stride ld must
 * still hold what it was filled with. */
static int pad_intact(const char *tag, const float *x, const float *x0, int ld,
                      int order, int m, int n, long *nviol) {
  int i, j, bad = 0;

  if (order == CblasColMajor) {
    for (j = 0; j < n; ++j)
      for (i = m; i < ld; ++i)
        if (x[(size_t)i + (size_t)j * ld] != x0[(size_t)i + (size_t)j * ld]) {
          violation(nviol, tag, i, j, "padding", 0.0, 0.0, 0.0);
          bad = 1;
        }
  } else {
    for (i = 0; i < m; ++i)
      for (j = n; j < ld; ++j)
        if (x[(size_t)i * ld + (size_t)j] != x0[(size_t)i * ld + (size_t)j]) {
          violation(nviol, tag, i, j, "padding", 0.0, 0.0, 0.0);
          bad = 1;
        }
  }
  return bad;
}

/* Compares one result element against the reference and reports a mismatch. */
static void check_el(const char *tag, int i, int j, double got, double want,
                     double tol, long *nviol, int *bad) {
  const double err = fabs(got - want);
  const double scale = fabs(want) > 1.0 ? fabs(want) : 1.0;

  if (err / scale > tol) {
    violation(nviol, tag, i, j, NULL, got, want, err / scale);
    *bad = 1;
  }
}

/* C := alpha*A*B + beta*C (side=left) or alpha*B*A + beta*C (side=right), A
 * symmetric with only its uplo triangle stored. Checks every element of C, the
 * ldc padding, and that A's other triangle is never read. */
static int symm_case(int order, int side, int uplo, int m, int n, double alpha,
                     double beta) {
  const int kdim = side ? n : m;
  const int lda = kdim + 3;
  const int ldb = (m > n ? m : n) + 4;
  const int ldc = (m > n ? m : n) + 5;
  const float sentinel = -12345.5f;
  const size_t asz = mat_size(lda, order, kdim, kdim);
  const size_t bsz = mat_size(ldb, order, m, n);
  const size_t csz = mat_size(ldc, order, m, n);
  float *A = (float *)malloc(asz * sizeof(float));
  float *B = (float *)malloc(bsz * sizeof(float));
  float *C = (float *)malloc(csz * sizeof(float));
  float *C0 = (float *)malloc(csz * sizeof(float));
  size_t k;
  int i, j, bad = 0;
  long nviol = 0;

  if (!A || !B || !C || !C0) {
    printf("symm: allocation failed at m=%d n=%d\n", m, n);
    free(A);
    free(B);
    free(C);
    free(C0);
    return 1;
  }

  for (k = 0; k < asz; ++k) A[k] = UNUSED_FILL;
  for (k = 0; k < bsz; ++k) B[k] = sentinel;
  for (k = 0; k < csz; ++k) C0[k] = C[k] = sentinel;

  for (j = 0; j < kdim; ++j)
    for (i = 0; i < kdim; ++i)
      if (sym_stored(uplo, i, j))
        mat_put(A, lda, order, i, j,
                (float)((double)((i * 5 + j * 3) % 17) - 8.0) / 8.0f);
  for (j = 0; j < n; ++j)
    for (i = 0; i < m; ++i) {
      const float v = (float)((double)((i * 7 + j * 11) % 23) - 11.0) / 9.0f;

      mat_put(B, ldb, order, i, j, v);
      mat_put(C0, ldc, order, i, j, v * 0.5f + 0.25f);
    }
  memcpy(C, C0, csz * sizeof(float));

  cblas_ssymm(order, side ? CblasRight : CblasLeft,
              uplo ? CblasLower : CblasUpper, m, n, (float)alpha, A, lda, B, ldb,
              (float)beta, C, ldc);

  for (j = 0; j < n; ++j)
    for (i = 0; i < m; ++i) {
      double ref = 0.0;
      int p;

      for (p = 0; p < kdim; ++p)
        ref += (side ? mat_get(B, ldb, order, i, p) *
                           sym_el(A, lda, order, uplo, p, j)
                     : sym_el(A, lda, order, uplo, i, p) *
                           mat_get(B, ldb, order, p, j));
      ref = alpha * ref + beta * (double)mat_get(C0, ldc, order, i, j);
      check_el("symm", i, j, mat_get(C, ldc, order, i, j), ref, 5e-4, &nviol,
               &bad);
    }
  bad |= pad_intact("symm", C, C0, ldc, order, m, n, &nviol);

  if (nviol > SYM_REPORT_MAX)
    printf("symm: ... %ld more violation(s) in this case\n",
           nviol - SYM_REPORT_MAX);
  printf("%s symm order=%d side=%d uplo=%d m=%d n=%d alpha=%g beta=%g\n",
         bad ? "FAIL" : "ok  ", order, side, uplo, m, n, alpha, beta);

  free(A);
  free(B);
  free(C);
  free(C0);
  return bad;
}

/* B := alpha*op(A)*B (side=left) or alpha*B*op(A) (side=right), in place, A
 * triangular with only its uplo triangle stored. The unused triangle and the
 * stored diagonal of a unit-diagonal factor carry UNUSED_FILL, so reading
 * either one shows up as a large error. */
static int trmm_case(int order, int side, int uplo, int trans, int unit, int m,
                     int n, double alpha) {
  const int kdim = side ? n : m;
  const int lda = kdim + 3;
  const int ldb = (m > n ? m : n) + 4;
  const float sentinel = -12345.5f;
  const size_t asz = mat_size(lda, order, kdim, kdim);
  const size_t bsz = mat_size(ldb, order, m, n);
  float *A = (float *)malloc(asz * sizeof(float));
  float *B = (float *)malloc(bsz * sizeof(float));
  float *B0 = (float *)malloc(bsz * sizeof(float));
  size_t k;
  int i, j, bad = 0;
  long nviol = 0;

  if (!A || !B || !B0) {
    printf("trmm: allocation failed at m=%d n=%d\n", m, n);
    free(A);
    free(B);
    free(B0);
    return 1;
  }

  for (k = 0; k < asz; ++k) A[k] = UNUSED_FILL;
  for (j = 0; j < kdim; ++j)
    for (i = 0; i < kdim; ++i)
      if (sym_stored(uplo, i, j) && !(i == j && unit))
        mat_put(A, lda, order, i, j,
                (float)((double)((i * 5 + j * 3) % 17) - 8.0) / 16.0f);
  for (k = 0; k < bsz; ++k) B0[k] = B[k] = sentinel;
  for (j = 0; j < n; ++j)
    for (i = 0; i < m; ++i)
      mat_put(B0, ldb, order, i, j,
              (float)((double)((i * 7 + j * 11) % 23) - 11.0) / 11.0f);
  memcpy(B, B0, bsz * sizeof(float));

  cblas_strmm(order, side ? CblasRight : CblasLeft,
              uplo ? CblasLower : CblasUpper,
              trans ? CblasTrans : CblasNoTrans,
              unit ? CblasUnit : CblasNonUnit, m, n, (float)alpha, A, lda, B,
              ldb);

  for (j = 0; j < n; ++j)
    for (i = 0; i < m; ++i) {
      double ref = 0.0;
      int p;

      for (p = 0; p < kdim; ++p)
        ref += (side ? tri_el(A, lda, order, uplo, trans, unit, p, j) *
                           mat_get(B0, ldb, order, i, p)
                     : tri_el(A, lda, order, uplo, trans, unit, i, p) *
                           mat_get(B0, ldb, order, p, j));
      check_el("trmm", i, j, mat_get(B, ldb, order, i, j), alpha * ref, 5e-4,
               &nviol, &bad);
    }
  bad |= pad_intact("trmm", B, B0, ldb, order, m, n, &nviol);

  if (nviol > SYM_REPORT_MAX)
    printf("trmm: ... %ld more violation(s) in this case\n",
           nviol - SYM_REPORT_MAX);
  printf("%s trmm order=%d side=%d uplo=%d trans=%d diag=%d m=%d n=%d alpha=%g\n",
         bad ? "FAIL" : "ok  ", order, side, uplo, trans, unit, m, n, alpha);

  free(A);
  free(B);
  free(B0);
  return bad;
}

/* Solves M*Y = X in double precision with partial pivoting. M (n x n) is
 * destroyed, X is n x nrhs with row stride ldx. Returns 0 if M is singular. */
static int lu_solve(double *M, int n, double *X, int nrhs, int ldx) {
  int k, i, j, r;

  for (k = 0; k < n; ++k) {
    int piv = k;
    double best = fabs(M[(size_t)k + (size_t)k * n]);

    for (i = k + 1; i < n; ++i) {
      const double v = fabs(M[(size_t)i + (size_t)k * n]);

      if (v > best) {
        best = v;
        piv = i;
      }
    }
    if (best < 1e-9) return 0;
    if (piv != k) {
      for (j = 0; j < n; ++j) {
        const double t = M[(size_t)k + (size_t)j * n];

        M[(size_t)k + (size_t)j * n] = M[(size_t)piv + (size_t)j * n];
        M[(size_t)piv + (size_t)j * n] = t;
      }
      for (r = 0; r < nrhs; ++r) {
        const double t = X[(size_t)k + (size_t)r * ldx];

        X[(size_t)k + (size_t)r * ldx] = X[(size_t)piv + (size_t)r * ldx];
        X[(size_t)piv + (size_t)r * ldx] = t;
      }
    }
    for (i = k + 1; i < n; ++i) {
      const double f =
          M[(size_t)i + (size_t)k * n] / M[(size_t)k + (size_t)k * n];

      if (f == 0.0) continue;
      for (j = k + 1; j < n; ++j)
        M[(size_t)i + (size_t)j * n] -= f * M[(size_t)k + (size_t)j * n];
      for (r = 0; r < nrhs; ++r)
        X[(size_t)i + (size_t)r * ldx] -= f * X[(size_t)k + (size_t)r * ldx];
    }
  }

  for (r = 0; r < nrhs; ++r)
    for (i = n - 1; i >= 0; --i) {
      double s = X[(size_t)i + (size_t)r * ldx];

      for (j = i + 1; j < n; ++j)
        s -= M[(size_t)i + (size_t)j * n] * X[(size_t)j + (size_t)r * ldx];
      X[(size_t)i + (size_t)r * ldx] = s / M[(size_t)i + (size_t)i * n];
    }
  return 1;
}

/* Solves op(A)*X = alpha*B (side=left) or X*op(A) = alpha*B (side=right), in
 * place on B, A triangular. The reference solves a dense copy of op(A) with a
 * pivoting LU, so it does not share the triangle logic under test. */
static int trsm_case(int order, int side, int uplo, int trans, int unit, int m,
                     int n, double alpha) {
  const int kdim = side ? n : m;
  const int lda = kdim + 3;
  const int ldb = (m > n ? m : n) + 4;
  const int nrhs = side ? m : n;
  const float sentinel = -12345.5f;
  const size_t asz = mat_size(lda, order, kdim, kdim);
  const size_t bsz = mat_size(ldb, order, m, n);
  float *A = (float *)malloc(asz * sizeof(float));
  float *B = (float *)malloc(bsz * sizeof(float));
  float *B0 = (float *)malloc(bsz * sizeof(float));
  double *M = (double *)malloc((size_t)kdim * (size_t)kdim * sizeof(double));
  double *X = (double *)malloc((size_t)kdim * (size_t)nrhs * sizeof(double));
  size_t k;
  int i, j, bad = 0;
  long nviol = 0;

  if (!A || !B || !B0 || !M || !X) {
    printf("trsm: allocation failed at m=%d n=%d\n", m, n);
    free(A);
    free(B);
    free(B0);
    free(M);
    free(X);
    return 1;
  }

  for (k = 0; k < asz; ++k) A[k] = UNUSED_FILL;
  for (j = 0; j < kdim; ++j)
    for (i = 0; i < kdim; ++i)
      if (sym_stored(uplo, i, j) && !(i == j && unit))
        mat_put(A, lda, order, i, j,
                (i == j) ? (float)(1.0 + 0.5 * ((j % 3) + 1) / 4.0)
                         : (float)((double)((i * 5 + j * 3) % 7) - 3.0) / 32.0f);
  for (k = 0; k < bsz; ++k) B0[k] = B[k] = sentinel;
  for (j = 0; j < n; ++j)
    for (i = 0; i < m; ++i)
      mat_put(B0, ldb, order, i, j,
              (float)((double)((i * 7 + j * 11) % 23) - 11.0) / 11.0f);
  memcpy(B, B0, bsz * sizeof(float));

  cblas_strsm(order, side ? CblasRight : CblasLeft,
              uplo ? CblasLower : CblasUpper,
              trans ? CblasTrans : CblasNoTrans,
              unit ? CblasUnit : CblasNonUnit, m, n, (float)alpha, A, lda, B,
              ldb);

  /* side=right is solved as op(A)^T * Y = alpha*B^T, so the dense copy is
   * transposed there; side=left solves against op(A) directly. */
  for (j = 0; j < kdim; ++j)
    for (i = 0; i < kdim; ++i)
      M[(size_t)i + (size_t)j * kdim] = tri_el(
          A, lda, order, uplo, trans, unit, side ? j : i, side ? i : j);
  /* X := alpha*B, or alpha*B^T for side=right, so both sides are a plain
   * solve against a matrix the code under test built nothing of. */
  for (j = 0; j < n; ++j)
    for (i = 0; i < m; ++i) {
      const double v = alpha * mat_get(B0, ldb, order, i, j);

      if (side)
        X[(size_t)j + (size_t)i * kdim] = v;
      else
        X[(size_t)i + (size_t)j * kdim] = v;
    }

  if (!lu_solve(M, kdim, X, nrhs, kdim)) {
    printf("FAIL trsm: reference matrix is singular (m=%d n=%d)\n", m, n);
    bad = 1;
  } else {
    for (j = 0; j < n; ++j)
      for (i = 0; i < m; ++i) {
        const double ref = side ? X[(size_t)j + (size_t)i * kdim]
                                : X[(size_t)i + (size_t)j * kdim];

        check_el("trsm", i, j, mat_get(B, ldb, order, i, j), ref, 1e-3, &nviol,
                 &bad);
      }
    bad |= pad_intact("trsm", B, B0, ldb, order, m, n, &nviol);
  }

  if (nviol > SYM_REPORT_MAX)
    printf("trsm: ... %ld more violation(s) in this case\n",
           nviol - SYM_REPORT_MAX);
  printf("%s trsm order=%d side=%d uplo=%d trans=%d diag=%d m=%d n=%d alpha=%g\n",
         bad ? "FAIL" : "ok  ", order, side, uplo, trans, unit, m, n, alpha);

  free(A);
  free(B);
  free(B0);
  free(M);
  free(X);
  return bad;
}

int main(int argc, char **argv) {
  setvbuf(stdout, NULL, _IONBF, 0);
  const int n = (argc > 1) ? atoi(argv[1]) : 2048;
  const double flops = 2.0 * (double)n * (double)n * (double)n;

  float *A = (float *)malloc((size_t)n * n * sizeof(float));
  float *B = (float *)malloc((size_t)n * n * sizeof(float));
  float *Cgpu = (float *)malloc((size_t)n * n * sizeof(float));
  float *Ccpu = (float *)malloc((size_t)n * n * sizeof(float));
  if (!A || !B || !Cgpu || !Ccpu) {
    fprintf(stderr, "allocation failed\n");
    return 1;
  }

  srand(1234);
  for (size_t i = 0; i < (size_t)n * n; ++i) {
    A[i] = (float)((rand() % 2001) - 1000) / 1000.0f;
    B[i] = (float)((rand() % 2001) - 1000) / 1000.0f;
  }
  memset(Ccpu, 0, (size_t)n * n * sizeof(float));

  /* GPU path (default when OPENBLAS_RDNA2_GPU is unset). */
  SETENV("OPENBLAS_RDNA2_GPU", "1");
  memset(Cgpu, 0, (size_t)n * n * sizeof(float));
  double t0 = now_s();
  cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0f, A, n,
              B, n, 0.0f, Cgpu, n);
  double t_gpu = now_s() - t0;

  /* CPU path, same library, GPU disabled. */
  SETENV("OPENBLAS_RDNA2_GPU", "0");
  t0 = now_s();
  cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0f, A, n,
              B, n, 0.0f, Ccpu, n);
  double t_cpu = now_s() - t0;

  double max_abs = 0.0, sum_sq = 0.0, c_sq = 0.0;
  for (size_t i = 0; i < (size_t)n * n; ++i) {
    double d = (double)Cgpu[i] - (double)Ccpu[i];
    if (d < 0) d = -d;
    if (d > max_abs) max_abs = d;
    sum_sq += d * d;
    c_sq += (double)Ccpu[i] * (double)Ccpu[i];
  }
  const double rms_rel = (c_sq > 0) ? sqrt(sum_sq / c_sq) : 0.0;

  printf("cblas_sgemm %dx%dx%d via libopenblas_rdna2 (build: TARGET=GFX1201 USE_HIP=1)\n", n, n, n);
  printf("  GPU path : %8.2f ms  %8.1f GFLOPS\n", t_gpu * 1e3, flops / t_gpu / 1e9);
  printf("  CPU path : %8.2f ms  %8.1f GFLOPS   (1 thread in this build)\n",
         t_cpu * 1e3, flops / t_cpu / 1e9);
  printf("  GPU vs CPU: max |diff| %.3e   rms_rel %.3e   %s\n", max_abs, rms_rel,
         (max_abs < 1e-3) ? "MATCH" : "MISMATCH");
  int failures = (max_abs < 1e-3) ? 0 : 1;

  free(A);
  free(B);
  free(Cgpu);
  free(Ccpu);

  /* ----------------------------------------------------------------------- */
  /* SYRK / SYR2K: odd shapes first, at the default size gate, so the real    */
  /* dispatch decision is exercised end to end.                               */
  /* ----------------------------------------------------------------------- */
  {
    int order, uplo, trans;
    int sym_failures = 0, cases = 0;
    /* n^2*k = 2.2e9, just past the bridge's default 2e9 threshold. */
    const int gate_n = 1300;

    SETENV("OPENBLAS_RDNA2_GPU", "1");
    printf("\nsyrk/syr2k above the default size gate (n=k=%d)\n", gate_n);
    sym_failures += sym_case(0, CblasColMajor, 0, 0, gate_n, gate_n, 1.0, 0.0);
    cases++;
    sym_failures += sym_case(1, CblasColMajor, 0, 0, gate_n, gate_n, 1.0, 0.0);
    cases++;

    /* The gate is not restored afterwards: the remaining cases are small on
     * purpose and must reach the device, or the CPU kernels on a GPU-less
     * machine, so the triangle and padding contracts are checked at sizes the
     * GEMM crossover would otherwise filter out. */
    SETENV("OPENBLAS_RDNA2_GPU_MIN_MNK", "0");
    printf("syrk/syr2k edge shapes (size gate forced open)\n");
    for (order = 0; order <= 1; ++order)
      for (uplo = 0; uplo <= 1; ++uplo)
        for (trans = 0; trans <= 1; ++trans) {
          const int ord = order ? CblasRowMajor : CblasColMajor;

          sym_failures += sym_case(0, ord, uplo, trans, 97, 61, 1.0, 0.0);
          cases++;
          sym_failures += sym_case(1, ord, uplo, trans, 97, 61, 1.0, 0.0);
          cases++;
          /* alpha/beta not 1/0: the triangle merge on the host path. */
          sym_failures += sym_case(0, ord, uplo, trans, 53, 130, 0.5, -0.25);
          cases++;
          sym_failures += sym_case(1, ord, uplo, trans, 53, 130, 0.5, -0.25);
          cases++;
        }

    printf("syrk/syr2k: %d cases, %d failure(s)\n", cases, sym_failures);
    failures += sym_failures;
  }

  /* ----------------------------------------------------------------------- */
  /* SYMM / TRMM / TRSM: one case per op above the default size gate, then     */
  /* every layout combination at shapes the gate would otherwise filter out.   */
  /* ----------------------------------------------------------------------- */
  {
    int order, side, uplo, trans, diag;
    int cases, sec_failures;
    /* m*n*kdim = 2.2e9, just past the bridge's default 2e9 threshold. */
    const int gate_mn = 1300;

    SETENV("OPENBLAS_RDNA2_GPU", "1");
    SETENV("OPENBLAS_RDNA2_GPU_MIN_MNK", "2000000000");
    printf("\nsymm/trmm/trsm at the default size gate (m=n=%d)\n", gate_mn);
    cases = sec_failures = 0;
    sec_failures += symm_case(CblasColMajor, 0, 0, gate_mn, gate_mn, 1.0, 0.0);
    cases++;
    sec_failures += trmm_case(CblasColMajor, 0, 0, 0, 0, gate_mn, gate_mn, 1.0);
    cases++;
    sec_failures += trsm_case(CblasColMajor, 0, 0, 0, 0, gate_mn, gate_mn, 1.0);
    cases++;
    printf("gate cases: %d cases, %d failure(s)\n", cases, sec_failures);
    failures += sec_failures;

    /* Every remaining case is small on purpose, so it reaches the device, or
     * the CPU kernels on a GPU-less machine, at shapes a GEMM crossover would
     * never route there. */
    SETENV("OPENBLAS_RDNA2_GPU_MIN_MNK", "0");

    cases = sec_failures = 0;
    printf("symm edge shapes (size gate forced open)\n");
    for (order = 0; order <= 1; ++order)
      for (side = 0; side <= 1; ++side)
        for (uplo = 0; uplo <= 1; ++uplo) {
          const int ord = order ? CblasRowMajor : CblasColMajor;

          sec_failures += symm_case(ord, side, uplo, 97, 61, 1.0, 0.0);
          cases++;
          sec_failures += symm_case(ord, side, uplo, 53, 130, 0.5, -0.25);
          cases++;
        }
    printf("symm: %d cases, %d failure(s)\n", cases, sec_failures);
    failures += sec_failures;

    cases = sec_failures = 0;
    printf("trmm edge shapes (size gate forced open)\n");
    for (order = 0; order <= 1; ++order)
      for (side = 0; side <= 1; ++side)
        for (uplo = 0; uplo <= 1; ++uplo)
          for (trans = 0; trans <= 1; ++trans)
            for (diag = 0; diag <= 1; ++diag) {
              const int ord = order ? CblasRowMajor : CblasColMajor;

              sec_failures += trmm_case(ord, side, uplo, trans, diag, 97, 61, 1.0);
              cases++;
              sec_failures +=
                  trmm_case(ord, side, uplo, trans, diag, 53, 130, -0.75);
              cases++;
            }
    printf("trmm: %d cases, %d failure(s)\n", cases, sec_failures);
    failures += sec_failures;

    cases = sec_failures = 0;
    printf("trsm edge shapes (size gate forced open)\n");
    for (order = 0; order <= 1; ++order)
      for (side = 0; side <= 1; ++side)
        for (uplo = 0; uplo <= 1; ++uplo)
          for (trans = 0; trans <= 1; ++trans)
            for (diag = 0; diag <= 1; ++diag) {
              const int ord = order ? CblasRowMajor : CblasColMajor;

              sec_failures += trsm_case(ord, side, uplo, trans, diag, 97, 61, 1.0);
              cases++;
              sec_failures +=
                  trsm_case(ord, side, uplo, trans, diag, 53, 130, 0.75);
              cases++;
            }
    printf("trsm: %d cases, %d failure(s)\n", cases, sec_failures);
    failures += sec_failures;
  }

  if (failures == 0)
    printf("ALL LEVEL-3 BRIDGE CHECKS PASSED\n");
  else
    printf("LEVEL-3 BRIDGE TEST FAILED (%d failure(s))\n", failures);
  return failures ? 2 : 0;
}
