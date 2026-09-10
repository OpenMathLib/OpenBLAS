/*
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
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef TEST_WASM_COMMON_H
#define TEST_WASM_COMMON_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cblas.h"

extern int g_fail;
extern int g_pass;

/* Sign domain of generated values. Never includes 0. */
typedef enum {
  FILL_POS = 0,  /* R^+ = (0, +inf) */
  FILL_NEG = 1,  /* R^- = (-inf, 0) */
  FILL_REAL = 2, /* R \ {0}: mixed signs */
} FillDomain;

/* Magnitude spread. NEAR stays O(1) at the large end so the mixed
 * absolute/relative tolerance in expect_close_* still gates kernel bugs. */
typedef enum {
  FILL_NEAR = 0, /* close to 0: log-uniform in [1e-4, 1] (f32) / [1e-8, 1] (f64) */
  FILL_FAR = 1,  /* far from 0: log-uniform in [1e2, 1e4] (f32) / [1e4, 1e8] (f64) */
} FillSpread;

typedef struct {
  FillDomain domain;
  FillSpread spread;
} FillSpec;

/* Active fill spec used by fill_vec_* / fill_mat_* / fill_f32 and friends.
 * Drivers cycle FILL_CASES across the size grid (see use_fill_case). */
extern FillSpec g_fill;

static const FillSpec FILL_CASES[] = {
    {FILL_POS, FILL_NEAR}, {FILL_POS, FILL_FAR},   {FILL_NEG, FILL_NEAR},
    {FILL_NEG, FILL_FAR},  {FILL_REAL, FILL_NEAR}, {FILL_REAL, FILL_FAR},
};
static const int NFILL_CASES = (int)(sizeof(FILL_CASES) / sizeof(FILL_CASES[0]));

static inline void use_fill_case(int i) {
  g_fill = FILL_CASES[(i < 0 ? 0 : i) % NFILL_CASES];
}

static inline const char *fill_spec_name(void) {
  static const char *names[3][2] = {
      {"R+ near 0", "R+ far from 0"},
      {"R- near 0", "R- far from 0"},
      {"R near 0", "R far from 0"},
  };
  return names[g_fill.domain][g_fill.spread];
}

static inline void *xmalloc(size_t n) {
  void *p = malloc(n);
  if (!p) {
    fprintf(stderr, "oom (%zu)\n", n);
    exit(1);
  }
  return p;
}

static inline unsigned fill_hash(int i, int seed) {
  unsigned x = (unsigned)i * 0x9e3779b9u + (unsigned)seed * 0x85ebca6bu;
  x ^= x >> 16;
  x *= 0x7feb352du;
  x ^= x >> 15;
  x *= 0x846ca68bu;
  x ^= x >> 16;
  return x;
}

/* Deterministic (0, 1) from (index, seed). */
static inline double fill_unit(int i, int seed) {
  return ((double)(fill_hash(i, seed) % 1000003u) + 1.0) / 1000004.0;
}

static inline float fill_mag_f32(int i, int seed, FillSpread spread) {
  float t = (float)fill_unit(i, seed);
  float log_lo = (spread == FILL_NEAR) ? logf(1e-4f) : logf(1e2f);
  float log_hi = (spread == FILL_NEAR) ? logf(1.0f) : logf(1e4f);
  return expf(log_lo + t * (log_hi - log_lo));
}

static inline double fill_mag_f64(int i, int seed, FillSpread spread) {
  double t = fill_unit(i, seed);
  double log_lo = (spread == FILL_NEAR) ? log(1e-8) : log(1e4);
  double log_hi = (spread == FILL_NEAR) ? log(1.0) : log(1e8);
  return exp(log_lo + t * (log_hi - log_lo));
}

static inline float fill_signed_f32(float mag, unsigned h, FillDomain domain) {
  switch (domain) {
  case FILL_POS:
    return mag;
  case FILL_NEG:
    return -mag;
  default:
    return (h & 1u) ? mag : -mag;
  }
}

static inline double fill_signed_f64(double mag, unsigned h, FillDomain domain) {
  switch (domain) {
  case FILL_POS:
    return mag;
  case FILL_NEG:
    return -mag;
  default:
    return (h & 1u) ? mag : -mag;
  }
}

/* Fill an n-vector from g_fill (strictly + / strictly − / mixed, near or far). */
static inline void fill_vec_f32(float *x, int n, int seed) {
  for (int i = 0; i < n; i++) {
    unsigned h = fill_hash(i, seed);
    x[i] = fill_signed_f32(fill_mag_f32(i, seed, g_fill.spread), h,
                           g_fill.domain);
  }
}

static inline void fill_vec_f64(double *x, int n, int seed) {
  for (int i = 0; i < n; i++) {
    unsigned h = fill_hash(i, seed);
    x[i] = fill_signed_f64(fill_mag_f64(i, seed, g_fill.spread), h,
                           g_fill.domain);
  }
}

/* Column-major matrix: lda rows allocated, cols columns. Padding in the
 * leading dimension is filled too, matching the existing dense fixtures. */
static inline void fill_mat_f32(float *A, int rows, int cols, int lda,
                                int seed) {
  (void)rows;
  fill_vec_f32(A, lda * cols, seed);
}

static inline void fill_mat_f64(double *A, int rows, int cols, int lda,
                                int seed) {
  (void)rows;
  fill_vec_f64(A, lda * cols, seed);
}

static inline void fill_f32(float *a, int n, int seed) {
  fill_vec_f32(a, n, seed);
}

static inline void fill_f64(double *a, int n, int seed) {
  fill_vec_f64(a, n, seed);
}

/* Complex as interleaved re,im pairs (length 2*n floats/doubles). Each part
 * follows the same real-line domain and spread as fill_vec_*. */
static inline void fill_c32(float *a, int n, int seed) {
  fill_vec_f32(a, 2 * n, seed);
}

static inline void fill_c64(double *a, int n, int seed) {
  fill_vec_f64(a, 2 * n, seed);
}

static inline void fill_mat_c32(float *A, int rows, int cols, int lda,
                                int seed) {
  (void)rows;
  fill_vec_f32(A, 2 * lda * cols, seed);
}

static inline void fill_mat_c64(double *A, int rows, int cols, int lda,
                                int seed) {
  (void)rows;
  fill_vec_f64(A, 2 * lda * cols, seed);
}

/* Dense triangular fixture for TRMV/TRSV/TRMM/TRSM tests (real only).
 * Intentionally independent of g_fill: O(1) diagonally-dominant entries keep
 * triangular solves well-conditioned. */
static inline void make_tri_f32(float *A, int n, int lda, enum CBLAS_UPLO uplo,
                                int unit) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      float v = 0.0f;
      if (i == j)
        v = unit ? 1.0f : (1.0f + 0.1f * (float)((i % 5) + 1));
      else if (uplo == CblasLower && i > j)
        v = 0.1f * (float)((i + j) % 5 + 1);
      else if (uplo == CblasUpper && i < j)
        v = 0.1f * (float)((i + j) % 5 + 1);
      A[i + j * lda] = v;
    }
}

static inline void make_tri_f64(double *A, int n, int lda, enum CBLAS_UPLO uplo,
                                int unit) {
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      double v = 0.0;
      if (i == j)
        v = unit ? 1.0 : (1.0 + 0.1 * (double)((i % 5) + 1));
      else if (uplo == CblasLower && i > j)
        v = 0.1 * (double)((i + j) % 5 + 1);
      else if (uplo == CblasUpper && i < j)
        v = 0.1 * (double)((i + j) % 5 + 1);
      A[i + j * lda] = v;
    }
}

static inline int close_f32(const float *got, const float *ref, int n, float tol,
                            float *out_maxe, float *out_maxv) {
  float maxe = 0.0f, maxv = 0.0f;
  for (int i = 0; i < n; i++) {
    float e = fabsf(got[i] - ref[i]);
    float v = fabsf(ref[i]);
    if (e > maxe)
      maxe = e;
    if (v > maxv)
      maxv = v;
  }
  if (out_maxe)
    *out_maxe = maxe;
  if (out_maxv)
    *out_maxv = maxv;
  return maxe <= tol * (1.0f + maxv);
}

static inline int close_f64(const double *got, const double *ref, int n,
                            double tol, double *out_maxe, double *out_maxv) {
  double maxe = 0.0, maxv = 0.0;
  for (int i = 0; i < n; i++) {
    double e = fabs(got[i] - ref[i]);
    double v = fabs(ref[i]);
    if (e > maxe)
      maxe = e;
    if (v > maxv)
      maxv = v;
  }
  if (out_maxe)
    *out_maxe = maxe;
  if (out_maxv)
    *out_maxv = maxv;
  return maxe <= tol * (1.0 + maxv);
}

static inline void pass_one(void) { g_pass++; }

static inline void fail_f32(const char *msg, float maxe, float maxv, float tol) {
  fprintf(stderr, "FAIL %s [%s] maxe=%.6g maxv=%.6g tol=%.6g\n", msg,
          fill_spec_name(), maxe, maxv, tol);
  g_fail++;
}

static inline void fail_f64(const char *msg, double maxe, double maxv,
                            double tol) {
  fprintf(stderr, "FAIL %s [%s] maxe=%.6g maxv=%.6g tol=%.6g\n", msg,
          fill_spec_name(), maxe, maxv, tol);
  g_fail++;
}

/* Pass if got[] matches the scalar reference within relative tolerance. */
static inline void expect_close_f32(const char *name, const float *got,
                                    const float *ref, int n, float tol) {
  float maxe;
  float maxv;

  if (close_f32(got, ref, n, tol, &maxe, &maxv))
    pass_one();
  else
    fail_f32(name, maxe, maxv, tol);
}

/* Pass if got[] matches the scalar reference within relative tolerance. */
static inline void expect_close_f64(const char *name, const double *got,
                                    const double *ref, int n, double tol) {
  double maxe;
  double maxv;

  if (close_f64(got, ref, n, tol, &maxe, &maxv))
    pass_one();
  else
    fail_f64(name, maxe, maxv, tol);
}

/* Pass if an integer-sized result matches the scalar reference exactly. */
static inline void expect_eq_size(const char *name, size_t got, size_t ref) {
  if (got == ref) {
    pass_one();
    return;
  }
  fprintf(stderr, "FAIL %s [%s] got=%zu ref=%zu\n", name, fill_spec_name(), got,
          ref);
  g_fail++;
}

void check_l1(void);
void check_l2(void);
void check_l3(void);

#endif /* TEST_WASM_COMMON_H */
