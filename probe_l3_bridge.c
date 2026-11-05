/* Direct probe of the RDNA2 level-3 bridges, bypassing the BLAS dispatch.
 *
 * This separates "is the bridge arithmetic right" from "does the interface hook
 * pass the right arguments". Every routine is compared against a naive scalar
 * reference computed here, so a mismatch is unambiguous.
 *
 * The unused triangle is filled with a large sentinel value instead of NaN so a
 * wrong-triangle bug shows up as a large numeric error rather than a NaN.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _WIN32
static void env_set(const char *k, const char *v) { _putenv_s(k, v); }
#else
static void env_set(const char *k, const char *v) { setenv(k, v, 1); }
#endif

extern int rdna2_oblas_ssymm_cm(long long, long long, double, const float *,
                                long long, const float *, long long, double,
                                float *, long long, int, int);
extern int rdna2_oblas_strmm_cm(long long, long long, double, const float *,
                                long long, float *, long long, int, int, int, int);
extern int rdna2_oblas_strsm_cm(long long, long long, double, const float *,
                                long long, float *, long long, int, int, int, int);
extern int rdna2_oblas_sgemm_cm(long long, long long, long long, double,
                                const float *, long long, const float *,
                                long long, double, float *, long long, int, int);

static unsigned rng = 999u;
static float rnd(void) {
  rng = rng * 1103515245u + 12345u;
  return (float)(((rng >> 8) & 0xFFFF) / 32768.0) - 1.0f;
}

/* Dense full triangle view. sentinel fills the triangle that must not be read. */
static void fill_sym(float *a, int kdim, int lda, int uplo) {
  for (int j = 0; j < kdim; ++j)
    for (int i = 0; i < kdim; ++i)
      a[i + j * lda] = ((uplo == 0) ? (i <= j) : (i >= j)) ? rnd() : 1.0e6f;
}

static void fill_tri(float *a, int kdim, int lda, int uplo, int unit,
                     int dominant) {
  for (int j = 0; j < kdim; ++j)
    for (int i = 0; i < kdim; ++i) {
      const int stored = (uplo == 0) ? (i <= j) : (i >= j);
      if (!stored)
        a[i + j * lda] = 1.0e6f;
      else if (i == j)
        a[i + j * lda] = unit ? 0.0f : (dominant ? (float)kdim + 2.0f : 1.5f);
      else
        a[i + j * lda] = dominant ? 0.25f * rnd() : rnd();
    }
}

/* Dense symmetric expansion, for the reference product only. */
static void expand_sym(const float *a, int kdim, int lda, int uplo, float *out) {
  for (int j = 0; j < kdim; ++j)
    for (int i = 0; i < kdim; ++i)
      out[i + j * kdim] =
          ((uplo == 0) ? (i <= j) : (i >= j)) ? a[i + j * lda] : a[j + i * lda];
}

/* Dense op(A) expansion, for the reference product only. */
static void expand_op(const float *a, int kdim, int lda, int uplo, int tr,
                      int unit, float *out) {
  for (int j = 0; j < kdim; ++j) {
    for (int i = 0; i < kdim; ++i) {
      const int si = tr ? j : i;
      const int sj = tr ? i : j;
      const int stored = (uplo == 0) ? (si <= sj) : (si >= sj);
      float v = 0.0f;
      if (stored) v = a[si + sj * lda];
      if (si == sj && unit) v = 1.0f;
      out[i + j * kdim] = v;
    }
  }
}

static int cmp(const char *name, const float *ref, const float *got, int n) {
  double maxabs = 0.0, scale = 0.0;
  int bad = -1;

  for (int i = 0; i < n; ++i) {
    const double d = fabs((double)got[i] - (double)ref[i]);
    if (d > maxabs) { maxabs = d; bad = i; }
    if (fabs((double)ref[i]) > scale) scale = fabs((double)ref[i]);
  }
  const double rel = (scale > 0.0) ? maxabs / scale : maxabs;
  printf("%-34s max|diff|=%.4e rel=%.4e %s", name, maxabs, rel,
         (rel < 1e-5) ? "OK" : "BAD");
  if (rel >= 1e-5 && bad >= 0)
    printf("  (first at %d: ref=%.6f got=%.6f)", bad, ref[bad], got[bad]);
  printf("\n");
  return rel < 1e-5;
}

static int t_symm(int m, int n, int side, int uplo) {
  const int kdim = side ? n : m;
  const double alpha = 0.75, beta = 0.5;
  float *a = malloc((size_t)kdim * kdim * sizeof(float));
  float *d = malloc((size_t)kdim * kdim * sizeof(float));
  float *b = malloc((size_t)m * n * sizeof(float));
  float *c0 = malloc((size_t)m * n * sizeof(float));
  float *ref = malloc((size_t)m * n * sizeof(float));
  float *got = malloc((size_t)m * n * sizeof(float));
  char name[64];
  int handled, i;

  fill_sym(a, kdim, kdim, uplo);
  expand_sym(a, kdim, kdim, uplo, d);
  for (i = 0; i < m * n; ++i) b[i] = rnd();
  for (i = 0; i < m * n; ++i) c0[i] = rnd();

  /* reference */
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < m; ++i) {
      double s = 0.0;
      for (int p = 0; p < kdim; ++p)
        s += (side ? b[i + p * m] * d[p + j * kdim]
                   : d[i + p * kdim] * b[p + j * m]);
      ref[i + j * m] = (float)(alpha * s + beta * c0[i + j * m]);
    }

  memcpy(got, c0, (size_t)m * n * sizeof(float));
  handled = rdna2_oblas_ssymm_cm(m, n, alpha, a, kdim, b, m, beta, got, m,
                                 side, uplo);
  snprintf(name, sizeof(name), "ssymm %dx%d %s/%s h=%d", m, n,
           side ? "R" : "L", uplo ? "L" : "U", handled);
  i = cmp(name, ref, got, m * n);

  free(a); free(d); free(b); free(c0); free(ref); free(got);
  return i;
}

static int t_trmm(int m, int n, int side, int uplo, int tr, int unit) {
  const int kdim = side ? n : m;
  const double alpha = 0.75;
  float *a = malloc((size_t)kdim * kdim * sizeof(float));
  float *d = malloc((size_t)kdim * kdim * sizeof(float));
  float *b = malloc((size_t)m * n * sizeof(float));
  float *ref = malloc((size_t)m * n * sizeof(float));
  float *got = malloc((size_t)m * n * sizeof(float));
  char name[64];
  int handled, i;

  fill_tri(a, kdim, kdim, uplo, unit, 0);
  expand_op(a, kdim, kdim, uplo, tr, unit, d);
  for (i = 0; i < m * n; ++i) b[i] = rnd();

  /* reference: C = alpha * op(A) * B  (side=left) or alpha * B * op(A) */
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < m; ++i) {
      double s = 0.0;
      for (int p = 0; p < kdim; ++p)
        s += (side ? b[i + p * m] * d[p + j * kdim]
                   : d[i + p * kdim] * b[p + j * m]);
      ref[i + j * m] = (float)(alpha * s);
    }

  memcpy(got, b, (size_t)m * n * sizeof(float));
  handled = rdna2_oblas_strmm_cm(m, n, alpha, a, kdim, got, m, side, uplo, tr,
                                 unit);
  snprintf(name, sizeof(name), "strmm %dx%d %s/%s/%s/%s h=%d", m, n,
           side ? "R" : "L", uplo ? "L" : "U", tr ? "T" : "N", unit ? "U" : "N",
           handled);
  i = cmp(name, ref, got, m * n);

  free(a); free(d); free(b); free(ref); free(got);
  return i;
}

static int t_trsm(int m, int n, int side, int uplo, int tr, int unit) {
  const int kdim = side ? n : m;
  const double alpha = 1.25;
  float *a = malloc((size_t)kdim * kdim * sizeof(float));
  float *d = malloc((size_t)kdim * kdim * sizeof(float));
  float *b = malloc((size_t)m * n * sizeof(float));
  float *ref = malloc((size_t)m * n * sizeof(float));
  float *got = malloc((size_t)m * n * sizeof(float));
  char name[64];
  int handled, i;

  fill_tri(a, kdim, kdim, uplo, unit, 1);
  expand_op(a, kdim, kdim, uplo, tr, unit, d);
  for (i = 0; i < m * n; ++i) b[i] = rnd();

  /* reference: solve op(A)*X = alpha*B (left) or X*op(A) = alpha*B (right),
   * by plain substitution on the dense expansion. */
  memcpy(ref, b, (size_t)m * n * sizeof(float));
  for (i = 0; i < m * n; ++i) ref[i] = (float)(alpha * (double)ref[i]);

  if (!side) {
    /* The solve matrix is op(A) = A or A^T, and transposing flips the triangle,
     * so branch on the triangle of op(A), not on the stored uplo. */
    const int triD = tr ? !uplo : uplo;
    if (triD == 0) {
      for (int j = 0; j < n; ++j)
        for (int i = kdim - 1; i >= 0; --i) {
          double s = ref[i + j * m];
          for (int q = i + 1; q < kdim; ++q)
            s -= (double)d[i + q * kdim] * (double)ref[q + j * m];
          ref[i + j * m] = (float)(s / (double)d[i + i * kdim]);
        }
    } else {
      for (int j = 0; j < n; ++j)
        for (int i = 0; i < kdim; ++i) {
          double s = ref[i + j * m];
          for (int q = 0; q < i; ++q)
            s -= (double)d[i + q * kdim] * (double)ref[q + j * m];
          ref[i + j * m] = (float)(s / (double)d[i + i * kdim]);
        }
    }
  } else {
    /* X*op(A) = Y  <=>  op(A)^T X^T = Y^T */
    float *yd = malloc((size_t)kdim * m * sizeof(float));
    for (int j = 0; j < m; ++j)
      for (int i = 0; i < kdim; ++i) yd[i + j * kdim] = ref[j + i * m];
    /* T = op(A)^T */
    float *t = malloc((size_t)kdim * kdim * sizeof(float));
    for (int j = 0; j < kdim; ++j)
      for (int i = 0; i < kdim; ++i) t[i + j * kdim] = d[j + i * kdim];
    const int tri = (tr ? uplo : !uplo);
    if (tri == 0) {
      for (int j = 0; j < m; ++j)
        for (int i = kdim - 1; i >= 0; --i) {
          double s = yd[i + j * kdim];
          for (int q = i + 1; q < kdim; ++q)
            s -= (double)t[i + q * kdim] * (double)yd[q + j * kdim];
          yd[i + j * kdim] = (float)(s / (double)t[i + i * kdim]);
        }
    } else {
      for (int j = 0; j < m; ++j)
        for (int i = 0; i < kdim; ++i) {
          double s = yd[i + j * kdim];
          for (int q = 0; q < i; ++q)
            s -= (double)t[i + q * kdim] * (double)yd[q + j * kdim];
          yd[i + j * kdim] = (float)(s / (double)t[i + i * kdim]);
        }
    }
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < m; ++i) ref[i + j * m] = yd[j + i * kdim];
    free(yd); free(t);
  }

  memcpy(got, b, (size_t)m * n * sizeof(float));
  handled = rdna2_oblas_strsm_cm(m, n, alpha, a, kdim, got, m, side, uplo, tr,
                                 unit);
  snprintf(name, sizeof(name), "strsm %dx%d %s/%s/%s/%s h=%d", m, n,
           side ? "R" : "L", uplo ? "L" : "U", tr ? "T" : "N", unit ? "U" : "N",
           handled);
  i = cmp(name, ref, got, m * n);

  free(a); free(d); free(b); free(ref); free(got);
  return i;
}

static int t_sgemm(int m, int n, int k) {
  float *a = malloc((size_t)m * k * sizeof(float));
  float *b = malloc((size_t)k * n * sizeof(float));
  float *ref = malloc((size_t)m * n * sizeof(float));
  float *got = malloc((size_t)m * n * sizeof(float));
  char name[64];
  int handled, i;

  for (i = 0; i < m * k; ++i) a[i] = rnd();
  for (i = 0; i < k * n; ++i) b[i] = rnd();

  for (int j = 0; j < n; ++j)
    for (int i2 = 0; i2 < m; ++i2) {
      double s = 0.0;
      for (int p = 0; p < k; ++p) s += (double)a[i2 + p * k] * (double)b[p + j * k];
      ref[i2 + j * m] = (float)s;
    }

  for (i = 0; i < m * n; ++i) got[i] = 0.0f;
  /* lda == k, ldb == k, ldc == m: the layout the bridge accepts. */
  handled = rdna2_oblas_sgemm_cm(m, n, k, 1.0, a, k, b, k, 0.0, got, m, 0, 0);
  snprintf(name, sizeof(name), "sgemm %dx%dx%d h=%d", m, n, k, handled);
  i = cmp(name, ref, got, m * n);

  free(a); free(b); free(ref); free(got);
  return i;
}

int main(void) {
  int pass = 1;
  const int ms[] = {5, 9, 130, 140, 300};

  env_set("OPENBLAS_RDNA2_GPU_MIN_MNK", "1");

  pass &= t_sgemm(5, 3, 4);
  pass &= t_sgemm(64, 32, 16);
  pass &= t_sgemm(200, 100, 96);

  for (int k = 0; k < 5; ++k) {
    const int m = ms[k], n = (k < 2) ? 3 : 7;
    for (int side = 0; side < 2; ++side) {
      for (int uplo = 0; uplo < 2; ++uplo) {
        pass &= t_symm(m, n, side, uplo);
        for (int tr = 0; tr < 2; ++tr) {
          pass &= t_trmm(m, n, side, uplo, tr, 0);
          pass &= t_trsm(m, n, side, uplo, tr, 0);
          if (m < 10) {
            pass &= t_trmm(m, n, side, uplo, tr, 1);
            pass &= t_trsm(m, n, side, uplo, tr, 1);
          }
        }
      }
    }
  }

  printf("\n%s\n", pass ? "BRIDGE PROBE ALL OK" : "BRIDGE PROBE FAILURES");
  return pass ? 0 : 1;
}
