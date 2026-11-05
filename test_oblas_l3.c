/* A/B correctness harness for the RDNA2 level-3 bridges.
 *
 * Each case is computed twice through the same library: once with the GPU path
 * forced off and once with it forced on, both with the size gate disabled so the
 * bridge is actually entered. The two results are then compared.
 *
 * The triangle of the symmetric/triangular operand that BLAS promises not to
 * read is filled with NaN, so a bridge that reads the wrong triangle fails
 * loudly instead of quietly returning a plausible number.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cblas.h"

#ifdef _WIN32
static void env_set(const char *k, const char *v) { _putenv_s(k, v); }
#else
static void env_set(const char *k, const char *v) { setenv(k, v, 1); }
#endif

static unsigned rng_state = 12345u;

static float rnd(void) {
  rng_state = rng_state * 1103515245u + 12345u;
  return (float)(((rng_state >> 8) & 0xFFFF) / 32768.0) - 1.0f;
}

/* Symmetric operand: write one triangle, NaN the other. */
static void fill_sym(float *a, int kdim, int lda, int uplo) {
  int i, j;

  for (j = 0; j < kdim; ++j) {
    for (i = 0; i < kdim; ++i) {
      const int stored = (uplo == 0) ? (i <= j) : (i >= j);
      a[i + j * lda] = stored ? rnd() : NAN;
    }
  }
}

/* Triangular operand: write one triangle, NaN the other. When dominant, the
 * diagonal is scaled up so the solve stays well conditioned; a unit diagonal is
 * stored as zero, which must be ignored. */
static void fill_tri(float *a, int kdim, int lda, int uplo, int unit_diag,
                     int dominant) {
  int i, j;

  for (j = 0; j < kdim; ++j) {
    for (i = 0; i < kdim; ++i) {
      const int stored = (uplo == 0) ? (i <= j) : (i >= j);

      if (!stored)
        a[i + j * lda] = NAN;
      else if (i == j)
        a[i + j * lda] = unit_diag ? 0.0f
                                   : (dominant ? (float)kdim + 2.0f : 1.5f);
      else
        a[i + j * lda] = dominant ? 0.25f * rnd() : rnd();
    }
  }
}

static void fill_rand(float *x, long n) {
  long i;

  for (i = 0; i < n; ++i) x[i] = rnd();
}

static int compare(const char *name, const float *ref, const float *got, long n) {
  long i;
  double maxabs = 0.0, maxrel = 0.0, scale = 0.0;

  for (i = 0; i < n; ++i) {
    const double d = fabs((double)got[i] - (double)ref[i]);
    const double s = fabs((double)ref[i]);

    if (isnan((double)got[i])) {
      printf("  %-28s FAIL: NaN at %ld\n", name, i);
      return 0;
    }
    if (d > maxabs) maxabs = d;
    if (s > scale) scale = s;
  }

  maxrel = (scale > 0.0) ? maxabs / scale : maxabs;
  printf("  %-28s max|diff|=%.3e  max_rel=%.3e  %s\n", name, maxabs, maxrel,
         (maxrel < 5e-5) ? "MATCH" : "MISMATCH");
  return maxrel < 5e-5;
}

/* ---------------------------------------------------------------- sgemm --- */
static int case_sgemm(int m, int n, int k) {
  float *a = malloc((size_t)m * k * sizeof(float));
  float *b = malloc((size_t)k * n * sizeof(float));
  float *ref = malloc((size_t)m * n * sizeof(float));
  float *got = malloc((size_t)m * n * sizeof(float));
  int ok;

  fill_rand(a, (long)m * k);
  fill_rand(b, (long)k * n);

  env_set("OPENBLAS_RDNA2_GPU", "0");
  memcpy(got, ref, (size_t)m * n * sizeof(float));
  for (int i = 0; i < m * n; ++i) got[i] = 0.0f;
  cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0f, a, m, b,
              k, 0.0f, got, m);
  memcpy(ref, got, (size_t)m * n * sizeof(float));

  env_set("OPENBLAS_RDNA2_GPU", "1");
  for (int i = 0; i < m * n; ++i) got[i] = 0.0f;
  cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0f, a, m, b,
              k, 0.0f, got, m);

  ok = compare("sgemm", ref, got, (long)m * n);

  free(a); free(b); free(ref); free(got);
  return ok;
}

/* ---------------------------------------------------------------- ssymm --- */
static void symm_once(int m, int n, int side, int uplo, float alpha, float beta,
                      const float *a, int lda, const float *b, int ldb, float *c,
                      int ldc) {
  cblas_ssymm(CblasColMajor, side ? CblasRight : CblasLeft,
              uplo ? CblasLower : CblasUpper, m, n, alpha, a, lda, b, ldb, beta,
              c, ldc);
}

static int case_symm(int m, int n, int side, int uplo, float alpha, float beta) {
  const int kdim = side ? n : m;
  float *a = malloc((size_t)kdim * kdim * sizeof(float));
  float *b = malloc((size_t)m * n * sizeof(float));
  float *ref = malloc((size_t)m * n * sizeof(float));
  float *got = malloc((size_t)m * n * sizeof(float));
  char name[64];
  int ok;

  fill_sym(a, kdim, kdim, uplo);
  fill_rand(b, (long)m * n);

  env_set("OPENBLAS_RDNA2_GPU", "0");
  memcpy(ref, b, (size_t)m * n * sizeof(float));
  symm_once(m, n, side, uplo, alpha, beta, a, kdim, b, m, ref, m);

  env_set("OPENBLAS_RDNA2_GPU", "1");
  memcpy(got, b, (size_t)m * n * sizeof(float));
  symm_once(m, n, side, uplo, alpha, beta, a, kdim, b, m, got, m);

  snprintf(name, sizeof(name), "ssymm %s/%s a=%.1f b=%.1f",
           side ? "R" : "L", uplo ? "L" : "U", (double)alpha, (double)beta);
  ok = compare(name, ref, got, (long)m * n);

  free(a); free(b); free(ref); free(got);
  return ok;
}

/* ---------------------------------------------------------------- strmm --- */
static void trmm_once(int m, int n, int side, int uplo, int trans, int unit,
                      float alpha, const float *a, int lda, float *b, int ldb) {
  cblas_strmm(CblasColMajor, side ? CblasRight : CblasLeft,
              uplo ? CblasLower : CblasUpper,
              trans ? CblasTrans : CblasNoTrans,
              unit ? CblasUnit : CblasNonUnit, m, n, alpha, a, lda, b, ldb);
}

static int case_trmm(int m, int n, int side, int uplo, int trans, int unit) {
  const int kdim = side ? n : m;
  const float alpha = 0.75f;
  float *a = malloc((size_t)kdim * kdim * sizeof(float));
  float *b0 = malloc((size_t)m * n * sizeof(float));
  float *ref = malloc((size_t)m * n * sizeof(float));
  float *got = malloc((size_t)m * n * sizeof(float));
  char name[64];
  int ok;

  fill_tri(a, kdim, kdim, uplo, unit, 0);
  fill_rand(b0, (long)m * n);

  env_set("OPENBLAS_RDNA2_GPU", "0");
  memcpy(ref, b0, (size_t)m * n * sizeof(float));
  trmm_once(m, n, side, uplo, trans, unit, alpha, a, kdim, ref, m);

  env_set("OPENBLAS_RDNA2_GPU", "1");
  memcpy(got, b0, (size_t)m * n * sizeof(float));
  trmm_once(m, n, side, uplo, trans, unit, alpha, a, kdim, got, m);

  snprintf(name, sizeof(name), "strmm %s/%s/%s/%s", side ? "R" : "L",
           uplo ? "L" : "U", trans ? "T" : "N", unit ? "U" : "N");
  ok = compare(name, ref, got, (long)m * n);

  free(a); free(b0); free(ref); free(got);
  return ok;
}

/* ---------------------------------------------------------------- strsm --- */
static void trsm_once(int m, int n, int side, int uplo, int trans, int unit,
                      float alpha, const float *a, int lda, float *b, int ldb) {
  cblas_strsm(CblasColMajor, side ? CblasRight : CblasLeft,
              uplo ? CblasLower : CblasUpper,
              trans ? CblasTrans : CblasNoTrans,
              unit ? CblasUnit : CblasNonUnit, m, n, alpha, a, lda, b, ldb);
}

static int case_trsm(int m, int n, int side, int uplo, int trans, int unit) {
  const int kdim = side ? n : m;
  const float alpha = 1.25f;
  float *a = malloc((size_t)kdim * kdim * sizeof(float));
  float *b0 = malloc((size_t)m * n * sizeof(float));
  float *ref = malloc((size_t)m * n * sizeof(float));
  float *got = malloc((size_t)m * n * sizeof(float));
  char name[64];
  int ok;

  fill_tri(a, kdim, kdim, uplo, unit, 1);
  fill_rand(b0, (long)m * n);

  env_set("OPENBLAS_RDNA2_GPU", "0");
  memcpy(ref, b0, (size_t)m * n * sizeof(float));
  trsm_once(m, n, side, uplo, trans, unit, alpha, a, kdim, ref, m);

  env_set("OPENBLAS_RDNA2_GPU", "1");
  memcpy(got, b0, (size_t)m * n * sizeof(float));
  trsm_once(m, n, side, uplo, trans, unit, alpha, a, kdim, got, m);

  snprintf(name, sizeof(name), "strsm %s/%s/%s/%s", side ? "R" : "L",
           uplo ? "L" : "U", trans ? "T" : "N", unit ? "U" : "N");
  ok = compare(name, ref, got, (long)m * n);

  free(a); free(b0); free(ref); free(got);
  return ok;
}

int main(int argc, char **argv) {
  const int big = (argc > 1) ? atoi(argv[1]) : 1;
  int pass = 1, i;

  /* Any size is exercised, so the measured crossover gate is disabled and the
   * bridge is always entered. */
  env_set("OPENBLAS_RDNA2_GPU_MIN_MNK", "1");

  printf("== sgemm ==\n");
  pass &= case_sgemm(256, 192, 128);
  if (big) pass &= case_sgemm(1024, 768, 512);

  printf("== ssymm ==\n");
  pass &= case_symm(200, 150, 0, 0, 1.0f, 0.0f);
  pass &= case_symm(200, 150, 0, 1, 0.75f, 0.25f);
  pass &= case_symm(180, 260, 1, 0, -0.5f, 1.0f);
  pass &= case_symm(180, 260, 1, 1, 1.0f, 0.0f);
  if (big) pass &= case_symm(400, 300, 0, 0, 1.0f, 0.0f);

  printf("== strmm ==\n");
  for (i = 0; i < 8; ++i) {
    const int side = i & 1, uplo = (i >> 1) & 1, trans = (i >> 2) & 1;
    pass &= case_trmm(200 + 40 * i, 150, side, uplo, trans, 0);
  }
  pass &= case_trmm(210, 160, 0, 0, 0, 1);
  pass &= case_trmm(210, 160, 1, 1, 1, 1);
  if (big) pass &= case_trmm(400, 300, 0, 1, 1, 0);

  printf("== strsm ==\n");
  for (i = 0; i < 8; ++i) {
    const int side = i & 1, uplo = (i >> 1) & 1, trans = (i >> 2) & 1;
    pass &= case_trsm(200 + 40 * i, 150, side, uplo, trans, 0);
  }
  pass &= case_trsm(210, 160, 0, 0, 0, 1);
  pass &= case_trsm(210, 160, 1, 1, 1, 1);
  if (big) pass &= case_trsm(400, 300, 0, 1, 1, 0);

  printf("\n%s\n", pass ? "ALL CASES MATCH" : "FAILURES PRESENT");
  return pass ? 0 : 1;
}
