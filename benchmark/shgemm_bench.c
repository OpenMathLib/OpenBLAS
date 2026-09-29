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
#include <time.h>
#include <math.h>
#include <cblas.h>

/* ------------------------------------------------------------------ */
/* hfloat16 is typedef uint16_t from cblas.h                          */
/* FP32 -> FP16 bit pattern                                           */
/* ------------------------------------------------------------------ */
static hfloat16 f2h(float f)
{
    if (f != f)    return (hfloat16)0x7e00u;
    if (f == 0.f)  return 0;
    unsigned s = (f < 0); f = fabsf(f);
    int e; float m = frexpf(f, &e); e += 14;
    if (e <= 0)    return (hfloat16)(s << 15);
    if (e >= 31)   return (hfloat16)((s << 15) | 0x7c00u);
    unsigned frac = (unsigned)((m * 2.f - 1.f) * 1024.f);
    return (hfloat16)((s << 15) | ((unsigned)e << 10) | (frac & 0x3ffu));
}

/* FP16 bit pattern -> FP32 */
static float h2f(hfloat16 h)
{
    unsigned x = h, s = (x>>15)&1, e = (x>>10)&0x1f, f = x&0x3ff;
    float v;
    if (e == 0)       v = ldexpf((float)f, -24);
    else if (e == 31) v = f ? (float)NAN : (float)INFINITY;
    else              v = ldexpf((float)(1024+f), (int)e-25);
    return s ? -v : v;
}

/* ------------------------------------------------------------------ */
/* Monotonic wall-clock in seconds                                     */
/* ------------------------------------------------------------------ */
static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ------------------------------------------------------------------ */
/* Fill FP32 array with reproducible pseudo-random values in [-0.5,0.5]*/
/* ------------------------------------------------------------------ */
static void fill_rand_f32(float *m, int n, unsigned int seed)
{
    unsigned int s = seed;
    for (int i = 0; i < n; i++) {
        s = s * 1664525u + 1013904223u;
        m[i] = (float)(s >> 8) / (float)(1 << 23) - 0.5f;
    }
}

/* ------------------------------------------------------------------ */
/* Convert FP32 array to FP16                                         */
/* ------------------------------------------------------------------ */
static void f32_to_f16(const float *src, hfloat16 *dst, int n)
{
    for (int i = 0; i < n; i++)
        dst[i] = f2h(src[i]);
}

/* ------------------------------------------------------------------ */
/* Max absolute difference between two FP32 arrays                    */
/* ------------------------------------------------------------------ */
static double max_diff(const float *a, const float *b, int n)
{
    double d = 0.0;
    for (int i = 0; i < n; i++) {
        double x = fabs((double)a[i] - (double)b[i]);
        if (x > d) d = x;
    }
    return d;
}

int main(int argc, char *argv[])
{
    int lo    = argc > 1 ? atoi(argv[1]) : 256;
    int hi    = argc > 2 ? atoi(argv[2]) : 4096;
    int step  = argc > 3 ? atoi(argv[3]) : 256;
    int iters = argc > 4 ? atoi(argv[4]) : 10;

    if (lo <= 0 || hi < lo || step <= 0 || iters <= 0) {
        fprintf(stderr,
                "Usage: %s [lo [hi [step [iters]]]]\n"
                "  lo=%d hi=%d step=%d iters=%d\n",
                argv[0], lo, hi, step, iters);
        return 1;
    }

    printf("\n%-6s  %-18s  %-14s  %-18s  %-14s  %-12s\n",
           "N", "SGEMM(GFLOPS)", "SGEMM_ms",
           "SHGEMM(GFLOPS)", "SHGEMM_ms", "MaxAbsDiff");
    printf("%-6s  %-18s  %-14s  %-18s  %-14s  %-12s\n",
           "------", "------------------", "--------------",
           "------------------", "--------------", "------------");

    for (int N = lo; N <= hi; N += step) {
        long long flops = 2LL * N * N * N;
        int       sz    = N * N;

        float      *Af32  = malloc((size_t)sz * sizeof(float));
        float      *Bf32  = malloc((size_t)sz * sizeof(float));
        float      *C_ref = calloc((size_t)sz,  sizeof(float));
        float      *C_shg = calloc((size_t)sz,  sizeof(float));
        hfloat16   *Af16  = malloc((size_t)sz * sizeof(hfloat16));
        hfloat16   *Bf16  = malloc((size_t)sz * sizeof(hfloat16));

        if (!Af32 || !Bf32 || !C_ref || !C_shg || !Af16 || !Bf16) {
            fprintf(stderr, "malloc failed for N=%d\n", N);
            return 1;
        }

        fill_rand_f32(Af32, sz, 42u);
        fill_rand_f32(Bf32, sz, 137u);
        f32_to_f16(Af32, Af16, sz);
        f32_to_f16(Bf32, Bf16, sz);

        /* ---- Warm up + time cblas_sgemm ---- */
        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    N, N, N, 1.0f, Af32, N, Bf32, N, 0.0f, C_ref, N);

        double t0 = now_sec();
        for (int it = 0; it < iters; it++)
            cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                        N, N, N, 1.0f, Af32, N, Bf32, N, 0.0f, C_ref, N);
        double sgemm_ms = ((now_sec() - t0) / iters) * 1000.0;
        double sgemm_gf = (double)flops / (sgemm_ms * 1e-3) / 1e9;

        /* ---- Warm up + time cblas_shgemm ---- */
        cblas_shgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                     N, N, N, 1.0f, Af16, N, Bf16, N, 0.0f, C_shg, N);

        t0 = now_sec();
        for (int it = 0; it < iters; it++)
            cblas_shgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                         N, N, N, 1.0f, Af16, N, Bf16, N, 0.0f, C_shg, N);
        double shg_ms = ((now_sec() - t0) / iters) * 1000.0;
        double shg_gf = (double)flops / (shg_ms * 1e-3) / 1e9;

        /* ---- Correctness check ---- */
        memset(C_ref, 0, (size_t)sz * sizeof(float));
        memset(C_shg, 0, (size_t)sz * sizeof(float));
        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    N, N, N, 1.0f, Af32, N, Bf32, N, 0.0f, C_ref, N);
        cblas_shgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                     N, N, N, 1.0f, Af16, N, Bf16, N, 0.0f, C_shg, N);

        double diff = max_diff(C_ref, C_shg, sz);

        printf("%-6d  %-18.2f  %-14.4f  %-18.2f  %-14.4f  %-12.2e\n",
               N, sgemm_gf, sgemm_ms, shg_gf, shg_ms, diff);

        free(Af32); free(Bf32); free(C_ref); free(C_shg);
        free(Af16); free(Bf16);
    }

    return 0;
}
