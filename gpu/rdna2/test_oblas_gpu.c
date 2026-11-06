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
static void sym_violation(long *count, int row, int col, const char *what,
                          double got, double want, double rel) {
  if ((*count)++ < SYM_REPORT_MAX) {
    if (what)
      printf("syrk/syr2k: %s element (%d,%d) was written\n", what, row, col);
    else
      printf("syrk/syr2k: (%d,%d) got %.6f want %.6f (rel %.2e)\n", row, col,
             got, want, rel);
  }
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
        sym_violation(&nviol, row, col, "padding", 0.0, 0.0, 0.0);
        bad = 1;
      }
      continue;
    }
    if (!sym_stored(uplo, row, col)) {
      if (C[idx] != C0[idx]) {
        sym_violation(&nviol, row, col, "unstored triangle", 0.0, 0.0, 0.0);
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
        sym_violation(&nviol, row, col, NULL, got, ref, err / scale);
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

  if (failures == 0)
    printf("ALL LEVEL-3 BRIDGE CHECKS PASSED\n");
  else
    printf("LEVEL-3 BRIDGE TEST FAILED (%d failure(s))\n", failures);
  return failures ? 2 : 0;
}
