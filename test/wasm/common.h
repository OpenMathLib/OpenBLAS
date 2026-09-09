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

static inline void *xmalloc(size_t n) {
  void *p = malloc(n);
  if (!p) {
    fprintf(stderr, "oom (%zu)\n", n);
    exit(1);
  }
  return p;
}

static inline void fill_f32(float *a, int n, int seed) {
  for (int i = 0; i < n; i++)
    a[i] = (float)((i * 17 + seed * 13) % 19) / 19.0f - 0.5f;
}

static inline void fill_f64(double *a, int n, int seed) {
  for (int i = 0; i < n; i++)
    a[i] = (double)((i * 17 + seed * 13) % 19) / 19.0 - 0.5;
}

/* Complex as interleaved re,im pairs (length 2*n floats/doubles). */
static inline void fill_c32(float *a, int n, int seed) {
  fill_f32(a, 2 * n, seed);
}

static inline void fill_c64(double *a, int n, int seed) {
  fill_f64(a, 2 * n, seed);
}

/* Dense triangular fixture for TRMV/TRSV/TRMM/TRSM tests (real only). */
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
  fprintf(stderr, "FAIL %s maxe=%.6g maxv=%.6g tol=%.6g\n", msg, maxe, maxv,
          tol);
  g_fail++;
}

static inline void fail_f64(const char *msg, double maxe, double maxv,
                            double tol) {
  fprintf(stderr, "FAIL %s maxe=%.6g maxv=%.6g tol=%.6g\n", msg, maxe, maxv,
          tol);
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
  fprintf(stderr, "FAIL %s got=%zu ref=%zu\n", name, got, ref);
  g_fail++;
}

void check_l1(void);
void check_l2(void);
void check_l3(void);

#endif /* TEST_WASM_COMMON_H */
