/*******************************************************************************
Copyright (c) 2025, The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the OpenBLAS project nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cblas.h"

static float h2f(hfloat16 h)
{
    unsigned x = h, s = (x >> 15) & 1, e = (x >> 10) & 0x1f, f = x & 0x3ff;
    if (e == 0)  return s ? -ldexpf((float)f, -24) : ldexpf((float)f, -24);
    if (e == 31) return f ? (float)NAN : (s ? -(float)INFINITY : (float)INFINITY);
    return ldexpf((float)(1024 + f), (int)e - 25) * (s ? -1.f : 1.f);
}

static hfloat16 f2h(float f)
{
    if (f != f)    return (hfloat16)0x7e00u;          /* NaN  */
    if (f == 0.f)  return 0;
    unsigned s = f < 0; f = fabsf(f);
    int e; float m = frexpf(f, &e); e += 14;
    if (e <= 0)    return (hfloat16)(s << 15);
    if (e >= 31)   return (hfloat16)((s << 15) | 0x7c00u);
    unsigned frac = (unsigned)((m * 2.f - 1.f) * 1024.f);
    return (hfloat16)((s << 15) | ((unsigned)e << 10) | (frac & 0x3ffu));
}

static void ref_mm(int m, int n, int k,
                   const hfloat16 *A, const hfloat16 *B, float *C)
{
    memset(C, 0, (size_t)m * n * sizeof(float));
    for (int j = 0; j < n; j++)
        for (int l = 0; l < k; l++) {
            float bv = h2f(B[l + j * k]);
            for (int i = 0; i < m; i++)
                C[i + j * m] += h2f(A[i + l * m]) * bv;
        }
}

/*
 * Usage: ./shgemm_test [M [N [K]]]
 *
 *   M, N, K   matrix dimensions (default 128 each)
 *
 * Runs exactly 50 iterations. After every iteration prints whether
 * the result matches the scalar reference.
 *
 * Examples:
 *   ./shgemm_test                    → 128×128×128, 50 iters
 *   ./shgemm_test 512                → 512×512×512, 50 iters
 *   ./shgemm_test 128 128 128        → 128×128×128, 50 iters
 *   ./shgemm_test 1024 2048 1024     → M=1024 N=2048 K=1024, 50 iters
 */
int main(int argc, char **argv)
{
    int M = (argc > 1) ? atoi(argv[1]) : 128;
    int N = (argc > 2) ? atoi(argv[2]) : M;
    int K = (argc > 3) ? atoi(argv[3]) : M;

    if (M <= 0 || N <= 0 || K <= 0) {
        fprintf(stderr,
                "Usage: %s [M [N [K]]]\n"
                "  M=%d  N=%d  K=%d\n",
                argv[0], M, N, K);
        return 1;
    }

#define ITERS 50

    printf("SHGEMM test: M=%d N=%d K=%d  iters=%d\n\n", M, N, K, ITERS);

    hfloat16 *A  = malloc((size_t)M * K * sizeof(hfloat16));
    hfloat16 *B  = malloc((size_t)K * N * sizeof(hfloat16));
    float    *CC = calloc((size_t)M * N, sizeof(float));
    float    *DD = calloc((size_t)M * N, sizeof(float));

    if (!A || !B || !CC || !DD) {
        fprintf(stderr, "malloc failed\n"); return 1;
    }

    /* Fill A and B with 1.0 in FP16; every C[i][j] must equal K exactly. */
    hfloat16 one = f2h(1.0f);
    for (int i = 0; i < M * K; i++) A[i] = one;
    for (int i = 0; i < K * N; i++) B[i] = one;

    /* Scalar reference — run once (O(M·N·K), slow for large sizes). */
    ref_mm(M, N, K, A, B, DD);

    float tol    = (float)K * 0.001f + 0.5f;
    int   total_pass = 0, total_fail = 0;

    for (int it = 1; it <= ITERS; it++) {
        memset(CC, 0, (size_t)M * N * sizeof(float));
        cblas_shgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                     M, N, K, 1.0f, A, M, B, K, 0.0f, CC, M);

        /* Compare against reference */
        int ok = 1, fi = -1;
        for (int i = 0; i < M * N; i++) {
            if (fabsf(CC[i] - DD[i]) > tol) { ok = 0; fi = i; break; }
        }

        if (ok) {
            printf("iter %2d/%d: MATRICES ARE SAME\n", it, ITERS);
            total_pass++;
        } else {
            printf("iter %2d/%d: MATRICES ARE DIFFERENT"
                   " [idx %d: got %.4f expected %.4f]\n",
                   it, ITERS, fi, CC[fi], DD[fi]);
            total_fail++;
        }
    }

    printf("\nTotal: %d SAME, %d DIFFERENT\n", total_pass, total_fail);

    free(A); free(B); free(CC); free(DD);
    return total_fail ? 1 : 0;
}

