/* Layer isolation for the RDNA2 level-3 hooks.
 *
 * For each case this computes:
 *   cpu  = cblas_* with the GPU path forced off, vs a naive scalar reference
 *   gpu  = cblas_* with the GPU path forced on,  vs the same naive reference
 *
 * If cpu is already wrong, the harness is misusing CBLAS (or the naive reference
 * is wrong) and the GPU result cannot be judged by it. If cpu is right and gpu
 * is wrong, the interface hook is passing bad arguments.
 *
 * The unused triangle holds zeros rather than NaN so both the CPU driver and the
 * bridge are exercised on well-defined input.
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

static unsigned rng = 4242u;
static float rnd(void) {
  rng = rng * 1103515245u + 12345u;
  return (float)(((rng >> 8) & 0xFFFF) / 32768.0) - 1.0f;
}

/* A with the stored triangle filled and the rest zeroed. */
static void fill_tri(float *a, int kdim, int lda, int uplo, int unit,
                     int dominant) {
  for (int j = 0; j < kdim; ++j)
    for (int i = 0; i < kdim; ++i) {
      const int stored = (uplo == 0) ? (i <= j) : (i >= j);
      if (!stored)
        a[i + j * lda] = 0.0f;
      else if (i == j)
        a[i + j * lda] = unit ? 1.0f : (dominant ? (float)kdim + 2.0f : 1.5f);
      else
        a[i + j * lda] = dominant ? 0.25f * rnd() : rnd();
    }
}

/* Dense op(A), and the triangle the dense matrix actually occupies. */
static void densify(const float *a, int kdim, int lda, int uplo, int tr, int unit,
                    float *d, int *tri) {
  for (int j = 0; j < kdim; ++j)
    for (int i = 0; i < kdim; ++i) {
      const int si = tr ? j : i;
      const int sj = tr ? i : j;
      const int stored = (uplo == 0) ? (si <= sj) : (si >= sj);
      float v = 0.0f;
      if (stored) v = a[si + sj * lda];
      if (si == sj && unit) v = 1.0f;
      d[i + j * kdim] = v;
    }
  /* Triangle occupied by the dense op(A). Transposing swaps the triangle, so
   * side=left (which solves with op(A)) uses this and side=right (which
   * solves with op(A)^T) uses its complement. */
  *tri = tr ? (uplo ? 0 : 1) : uplo;
}

static double err(const float *ref, const float *got, int n) {
  double maxabs = 0.0, scale = 0.0;
  for (int i = 0; i < n; ++i) {
    const double d = fabs((double)got[i] - (double)ref[i]);
    if (d > maxabs) maxabs = d;
    if (fabs((double)ref[i]) > scale) scale = fabs((double)ref[i]);
  }
  return (scale > 0.0) ? maxabs / scale : maxabs;
}

static void one_trmm(int m, int n, int side, int uplo, int tr, int unit) {
  const int kdim = side ? n : m;
  const double alpha = 0.75;
  float *a = malloc((size_t)kdim * kdim * sizeof(float));
  float *d = malloc((size_t)kdim * kdim * sizeof(float));
  float *b0 = malloc((size_t)m * n * sizeof(float));
  float *ref = malloc((size_t)m * n * sizeof(float));
  float *cpu = malloc((size_t)m * n * sizeof(float));
  float *gpu = malloc((size_t)m * n * sizeof(float));
  int tri;

  fill_tri(a, kdim, kdim, uplo, unit, 0);
  densify(a, kdim, kdim, uplo, tr, unit, d, &tri);
  for (int i = 0; i < m * n; ++i) b0[i] = rnd();

  /* naive: C = alpha*op(A)*B  or  alpha*B*op(A) */
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < m; ++i) {
      double s = 0.0;
      for (int p = 0; p < kdim; ++p)
        s += side ? (double)b0[i + p * m] * d[p + j * kdim]
                  : (double)d[i + p * kdim] * b0[p + j * m];
      ref[i + j * m] = (float)(alpha * s);
    }

  env_set("OPENBLAS_RDNA2_GPU", "0");
  memcpy(cpu, b0, (size_t)m * n * sizeof(float));
  cblas_strmm(CblasColMajor, side ? CblasRight : CblasLeft,
              uplo ? CblasLower : CblasUpper, tr ? CblasTrans : CblasNoTrans,
              unit ? CblasUnit : CblasNonUnit, m, n, (float)alpha, a, kdim, cpu,
              m);

  env_set("OPENBLAS_RDNA2_GPU", "1");
  memcpy(gpu, b0, (size_t)m * n * sizeof(float));
  cblas_strmm(CblasColMajor, side ? CblasRight : CblasLeft,
              uplo ? CblasLower : CblasUpper, tr ? CblasTrans : CblasNoTrans,
              unit ? CblasUnit : CblasNonUnit, m, n, (float)alpha, a, kdim, gpu,
              m);

  printf("trmm %3dx%-3d %s/%s/%s/%s  cpu_err=%.2e  gpu_err=%.2e  %s\n", m, n,
         side ? "R" : "L", uplo ? "L" : "U", tr ? "T" : "N", unit ? "U" : "N",
         err(ref, cpu, m * n), err(ref, gpu, m * n),
         (err(ref, gpu, m * n) < 1e-5) ? "ok" : "BAD");

  free(a); free(d); free(b0); free(ref); free(cpu); free(gpu);
}

static void one_trsm(int m, int n, int side, int uplo, int tr, int unit) {
  const int kdim = side ? n : m;
  const double alpha = 1.25;
  float *a = malloc((size_t)kdim * kdim * sizeof(float));
  float *d = malloc((size_t)kdim * kdim * sizeof(float));
  float *b0 = malloc((size_t)m * n * sizeof(float));
  float *ref = malloc((size_t)m * n * sizeof(float));
  float *cpu = malloc((size_t)m * n * sizeof(float));
  float *gpu = malloc((size_t)m * n * sizeof(float));
  int tri;

  fill_tri(a, kdim, kdim, uplo, unit, 1);
  densify(a, kdim, kdim, uplo, tr, unit, d, &tri);
  for (int i = 0; i < m * n; ++i) b0[i] = rnd();

  /* naive solve, working on the dense factor */
  if (!side) {
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < m; ++i) ref[i + j * m] = (float)(alpha * b0[i + j * m]);
    if (tri == 0) {
      for (int j = 0; j < n; ++j)
        for (int i = kdim - 1; i >= 0; --i) {
          double s = ref[i + j * m];
          for (int q = i + 1; q < kdim; ++q)
            s -= (double)d[i + q * kdim] * ref[q + j * m];
          ref[i + j * m] = (float)(s / d[i + i * kdim]);
        }
    } else {
      for (int j = 0; j < n; ++j)
        for (int i = 0; i < kdim; ++i) {
          double s = ref[i + j * m];
          for (int q = 0; q < i; ++q)
            s -= (double)d[i + q * kdim] * ref[q + j * m];
          ref[i + j * m] = (float)(s / d[i + i * kdim]);
        }
    }
  } else {
    /* X*op(A) = Y  <=>  op(A)^T X^T = Y^T */
    float *yd = malloc((size_t)kdim * m * sizeof(float));
    float *t = malloc((size_t)kdim * kdim * sizeof(float));
    const int triT = 1 - tri;  /* T = op(A)^T flips the triangle */
    for (int j = 0; j < m; ++j)
      for (int i = 0; i < kdim; ++i) yd[i + j * kdim] = (float)(alpha * b0[j + i * m]);
    for (int j = 0; j < kdim; ++j)
      for (int i = 0; i < kdim; ++i) t[i + j * kdim] = d[j + i * kdim];
    if (triT == 0) {
      for (int j = 0; j < m; ++j)
        for (int i = kdim - 1; i >= 0; --i) {
          double s = yd[i + j * kdim];
          for (int q = i + 1; q < kdim; ++q) s -= (double)t[i + q * kdim] * yd[q + j * kdim];
          yd[i + j * kdim] = (float)(s / t[i + i * kdim]);
        }
    } else {
      for (int j = 0; j < m; ++j)
        for (int i = 0; i < kdim; ++i) {
          double s = yd[i + j * kdim];
          for (int q = 0; q < i; ++q) s -= (double)t[i + q * kdim] * yd[q + j * kdim];
          yd[i + j * kdim] = (float)(s / t[i + i * kdim]);
        }
    }
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < m; ++i) ref[i + j * m] = yd[j + i * kdim];
    free(yd); free(t);
  }

  env_set("OPENBLAS_RDNA2_GPU", "0");
  memcpy(cpu, b0, (size_t)m * n * sizeof(float));
  cblas_strsm(CblasColMajor, side ? CblasRight : CblasLeft,
              uplo ? CblasLower : CblasUpper, tr ? CblasTrans : CblasNoTrans,
              unit ? CblasUnit : CblasNonUnit, m, n, (float)alpha, a, kdim, cpu,
              m);

  env_set("OPENBLAS_RDNA2_GPU", "1");
  memcpy(gpu, b0, (size_t)m * n * sizeof(float));
  cblas_strsm(CblasColMajor, side ? CblasRight : CblasLeft,
              uplo ? CblasLower : CblasUpper, tr ? CblasTrans : CblasNoTrans,
              unit ? CblasUnit : CblasNonUnit, m, n, (float)alpha, a, kdim, gpu,
              m);

  printf("trsm %3dx%-3d %s/%s/%s/%s  cpu_err=%.2e  gpu_err=%.2e  %s\n", m, n,
         side ? "R" : "L", uplo ? "L" : "U", tr ? "T" : "N", unit ? "U" : "N",
         err(ref, cpu, m * n), err(ref, gpu, m * n),
         (err(ref, gpu, m * n) < 1e-5) ? "ok" : "BAD");

  free(a); free(d); free(b0); free(ref); free(cpu); free(gpu);
}

int main(void) {
  const int sizes[3][2] = {{8, 5}, {300, 5}, {5, 300}};

  env_set("OPENBLAS_RDNA2_GPU_MIN_MNK", "1");

  for (int s = 0; s < 3; ++s)
    for (int side = 0; side < 2; ++side)
      for (int uplo = 0; uplo < 2; ++uplo)
        for (int tr = 0; tr < 2; ++tr) {
          one_trmm(sizes[s][0], sizes[s][1], side, uplo, tr, 0);
          one_trsm(sizes[s][0], sizes[s][1], side, uplo, tr, 0);
        }

  /* unit diagonal, small sizes only */
  for (int side = 0; side < 2; ++side)
    for (int uplo = 0; uplo < 2; ++uplo)
      for (int tr = 0; tr < 2; ++tr) {
        one_trmm(9, 6, side, uplo, tr, 1);
        one_trsm(9, 6, side, uplo, tr, 1);
      }

  return 0;
}
