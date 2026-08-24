//#include <algorithm>
#include <arm_neon.h>
//#include <cstddef>
#include <stddef.h>
#include <stdalign.h>
#include <stdbool.h>
#include <stdio.h>
#include "common.h"
#ifndef stdmin
#define stdmin(a,b)   (a>b? b:a)
#endif

#define KERNEL_ALPHA 0
#define USE_VECTORIZED_PACKING 1


// ========================================================================
// OPTIMIZED CACHE BLOCKING PARAMETERS
// ========================================================================
#define MC  256
#define KC  512
#define NC  1024

// ========================================================================
// 16x16 PURE-SME MACRO-KERNEL (Fused Beta Scaling)
// ========================================================================
static inline void dgemm_sme_compute_16x16_tile(int current_K, const double *A_ptr, const double *B_ptr, double *C_ptr, ptrdiff_t ldc, int beta_mode, const double *beta_ptr)
{
    ptrdiff_t ldc_bytes = ldc * sizeof(double);

    asm volatile(
    "smstart\n\t"
    // Enable all 64-bit (double precision) lanes in predicate register p0
    "ptrue p0.d\n\t"

    // =========================================================
    // PHASE 1: INITIALIZE 'ZA' ACCUMULATOR
    // =========================================================
    "cmp %w[beta_mode], #0\n\t"
    "b.eq 10f\n\t" // Jump to zero_za
    "cmp %w[beta_mode], #1\n\t"
    "b.eq 11f\n\t" // Jump to scale_beta
    "b 12f\n\t"    // Jump to direct_load (default)

    // --- PATH 0: beta == 0.0 (Instant Zero) ---
    "10:\n\t"
    "zero {za}\n\t"
    "b 19f\n\t" // Done with Phase 1

    // --- PATH 1: beta != 1.0 (Vectorized Load & Scale) ---
    "11:\n\t"
    "ld1rd z31.d, p0/z, [%[beta_ptr]]\n\t" // Load and replicate beta into z31
    "mov w12, #0\n\t"
    "mov x13, %[c]\n\t"

    "110:\n\t"
    "ld1d z0.d, p0/z, [x13]\n\t"
    "fmul z0.d, p0/m, z0.d, z31.d\n\t"
    "mova za0v.d[w12, 0], p0/m, z0.d\n\t"

    "add x14, x13, #64\n\t"
    "ld1d z1.d, p0/z, [x14]\n\t"
    "fmul z1.d, p0/m, z1.d, z31.d\n\t"
    "mova za2v.d[w12, 0], p0/m, z1.d\n\t"

    "add x13, x13, %[ldc_bytes]\n\t"
    "add w12, w12, #1\n\t"
    "cmp w12, #8\n\t"
    "b.ne 110b\n\t"

    "mov w15, #0\n\t"
    "111:\n\t"
    "ld1d z0.d, p0/z, [x13]\n\t"
    "fmul z0.d, p0/m, z0.d, z31.d\n\t"
    "mova za1v.d[w15, 0], p0/m, z0.d\n\t"

    "add x14, x13, #64\n\t"
    "ld1d z1.d, p0/z, [x14]\n\t"
    "fmul z1.d, p0/m, z1.d, z31.d\n\t"
    "mova za3v.d[w15, 0], p0/m, z1.d\n\t"

    "add x13, x13, %[ldc_bytes]\n\t"
    "add w15, w15, #1\n\t"
    "cmp w15, #8\n\t"
    "b.ne 111b\n\t"
    "b 19f\n\t" // Done with Phase 1

    // --- PATH 2: k > 0 or beta == 1.0 (Standard Direct Load) ---
    "12:\n\t"
    "mov w12, #0\n\t"
    "mov x13, %[c]\n\t"

    "100:\n\t"
    "ld1d {za0v.d[w12, 0]}, p0/z, [x13]\n\t"
    "add x14, x13, #64\n\t"
    "ld1d {za2v.d[w12, 0]}, p0/z, [x14]\n\t"
    "add x13, x13, %[ldc_bytes]\n\t"
    "add w12, w12, #1\n\t"
    "cmp w12, #8\n\t"
    "b.ne 100b\n\t"

    "mov w15, #0\n\t"
    "101:\n\t"
    "ld1d {za1v.d[w15, 0]}, p0/z, [x13]\n\t"
    "add x14, x13, #64\n\t"
    "ld1d {za3v.d[w15, 0]}, p0/z, [x14]\n\t"
    "add x13, x13, %[ldc_bytes]\n\t"
    "add w15, w15, #1\n\t"
    "cmp w15, #8\n\t"
    "b.ne 101b\n\t"

    "19:\n\t" // Entry point for Phase 2

    // =========================================================
    // PHASE 2: COMPUTE (ZA += A * B)
    // =========================================================
    "mov w10, %w[k]\n\t"
    "cbz w10, 3f\n\t"

    "lsr w11, w10, #2\n\t"
    "cbz w11, 12f\n\t"

    "11:\n\t"
    "ld1d z0.d, p0/z, [%[a], #0, mul vl]\n\t"
    "ld1d z1.d, p0/z, [%[a], #1, mul vl]\n\t"
    "ld1d z2.d, p0/z, [%[b], #0, mul vl]\n\t"
    "ld1d z3.d, p0/z, [%[b], #1, mul vl]\n\t"

    "ld1d z4.d, p0/z, [%[a], #2, mul vl]\n\t"
    "ld1d z5.d, p0/z, [%[a], #3, mul vl]\n\t"
    "ld1d z6.d, p0/z, [%[b], #2, mul vl]\n\t"
    "ld1d z7.d, p0/z, [%[b], #3, mul vl]\n\t"

    "ld1d z8.d, p0/z, [%[a], #4, mul vl]\n\t"
    "ld1d z9.d, p0/z, [%[a], #5, mul vl]\n\t"
    "ld1d z10.d, p0/z, [%[b], #4, mul vl]\n\t"
    "ld1d z11.d, p0/z, [%[b], #5, mul vl]\n\t"

    "ld1d z12.d, p0/z, [%[a], #6, mul vl]\n\t"
    "ld1d z13.d, p0/z, [%[a], #7, mul vl]\n\t"
    "ld1d z14.d, p0/z, [%[b], #6, mul vl]\n\t"
    "ld1d z15.d, p0/z, [%[b], #7, mul vl]\n\t"

    "fmopa za0.d, p0/m, p0/m, z0.d, z2.d\n\t"
    "fmopa za1.d, p0/m, p0/m, z0.d, z3.d\n\t"
    "fmopa za2.d, p0/m, p0/m, z1.d, z2.d\n\t"
    "fmopa za3.d, p0/m, p0/m, z1.d, z3.d\n\t"

    "fmopa za0.d, p0/m, p0/m, z4.d, z6.d\n\t"
    "fmopa za1.d, p0/m, p0/m, z4.d, z7.d\n\t"
    "fmopa za2.d, p0/m, p0/m, z5.d, z6.d\n\t"
    "fmopa za3.d, p0/m, p0/m, z5.d, z7.d\n\t"

    "fmopa za0.d, p0/m, p0/m, z8.d, z10.d\n\t"
    "fmopa za1.d, p0/m, p0/m, z8.d, z11.d\n\t"
    "fmopa za2.d, p0/m, p0/m, z9.d, z10.d\n\t"
    "fmopa za3.d, p0/m, p0/m, z9.d, z11.d\n\t"

    "fmopa za0.d, p0/m, p0/m, z12.d, z14.d\n\t"
    "fmopa za1.d, p0/m, p0/m, z12.d, z15.d\n\t"
    "fmopa za2.d, p0/m, p0/m, z13.d, z14.d\n\t"
    "fmopa za3.d, p0/m, p0/m, z13.d, z15.d\n\t"

    "add %[a], %[a], #512\n\t"
    "add %[b], %[b], #512\n\t"
    "subs w11, w11, #1\n\t"
    "b.ne 11b\n\t"

    // --- Remainder Loop ---
    "12:\n\t"
    "and w10, w10, #3\n\t"
    "cbz w10, 3f\n\t"

    "13:\n\t"
    "ld1d z0.d, p0/z, [%[a], #0, mul vl]\n\t"
    "ld1d z1.d, p0/z, [%[a], #1, mul vl]\n\t"
    "ld1d z2.d, p0/z, [%[b], #0, mul vl]\n\t"
    "ld1d z3.d, p0/z, [%[b], #1, mul vl]\n\t"
    "fmopa za0.d, p0/m, p0/m, z0.d, z2.d\n\t"
    "fmopa za1.d, p0/m, p0/m, z0.d, z3.d\n\t"
    "fmopa za2.d, p0/m, p0/m, z1.d, z2.d\n\t"
    "fmopa za3.d, p0/m, p0/m, z1.d, z3.d\n\t"
    "add %[a], %[a], #128\n\t"
    "add %[b], %[b], #128\n\t"
    "subs w10, w10, #1\n\t"
    "b.ne 13b\n\t"

    // =========================================================
    // PHASE 3: STORE 'ZA' ACCUMULATOR BACK TO MATRIX 'C'
    // =========================================================
    "3:\n\t"
    "mov w12, #0\n\t"
    "mov x13, %[c]\n\t"

    "200:\n\t"
    "st1d {za0v.d[w12, 0]}, p0, [x13]\n\t"
    "add x14, x13, #64\n\t"
    "st1d {za2v.d[w12, 0]}, p0, [x14]\n\t"
    "add x13, x13, %[ldc_bytes]\n\t"
    "add w12, w12, #1\n\t"
    "cmp w12, #8\n\t"
    "b.ne 200b\n\t"

    "mov w15, #0\n\t"
    "201:\n\t"
    "st1d {za1v.d[w15, 0]}, p0, [x13]\n\t"
    "add x14, x13, #64\n\t"
    "st1d {za3v.d[w15, 0]}, p0, [x14]\n\t"
    "add x13, x13, %[ldc_bytes]\n\t"
    "add w15, w15, #1\n\t"
    "cmp w15, #8\n\t"
    "b.ne 201b\n\t"
    "smstop\n\t"
    "msr fpsr, xzr\n\t"
    : [a] "+&r"(A_ptr), [b] "+&r"(B_ptr)
    : [k] "r"(current_K), [c] "r"(C_ptr), [ldc_bytes] "r"(ldc_bytes), [beta_mode] "r"(beta_mode), [beta_ptr] "r"(beta_ptr)
    : "x10", "x11", "x12", "x13", "x14", "w15", "memory", "cc",
       "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9",
       "v10", "v11", "v12", "v13", "v14", "v15", "v16", "v17", "v18",
       "v19", "v20", "v21", "v22", "v23", "v24", "v25", "v26", "v27",
       "v28", "v29", "v30", "v31",
       "z0", "z1", "z2", "z3", "z4", "z5", "z6", "z7", "z8", "z9",
       "z10", "z11", "z12", "z13", "z14", "z15", "z16", "z17", "z18",
       "z19", "z20", "z21", "z22", "z23", "z24", "z25", "z26", "z27",
       "z28", "z29", "z30", "z31", "za", "p0", "p1", "p2", "p3", "p4",
       "p5", "p6", "p7", "p8", "p9", "p10", "p11", "p12", "p13", "p14", "p15");


}

// ========================================================================
// VERSION 0: C = alpha * A * B + beta * C (NN)
// ========================================================================
static void dgemm_sme_NN(int M, int N, int K, double alpha, const double *A, int lda, const double *B, int ldb, double beta, double *C, int ldc)
{
#if PROFILING
    CALI_CXX_MARK_FUNCTION;
#endif
    // Hardware bypass for skipped main loops
    if (alpha == 0.0 || K == 0) {
        if (beta != 1.0) {
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < M; ++i) {
                    C[i + j * ldc] = (beta == 0.0) ? 0.0 : C[i + j * ldc] * beta;
                }
            }
        }
        return;
    }

    alignas(256) /*thread_local*/ static double A_pack[MC * KC];
    alignas(256) /*thread_local*/ static double B_pack[KC * NC];

#pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < N; j += NC) {
#if PROFILING
        CALI_MARK_BEGIN("loop_j");
#endif
        int current_N = stdmin(NC, N - j);
        int N_pad = (current_N + 15) & ~15;

        for (int k = 0; k < K; k += KC) {
#if PROFILING
            CALI_MARK_BEGIN("loop_k");
#endif
            int current_K = stdmin(KC, K - k);
            ptrdiff_t panel_stride_B = 16 * (ptrdiff_t)current_K;

            // Route the beta scaling natively inside the SME Assembly
            int beta_mode = 2; // Default: K>0, directly load C
            if (k == 0) {
                if (beta == 0.0) {
                    beta_mode = 0;
                }
                else if (beta != 1.0) {
                    beta_mode = 1;
                }
            }

#if PROFILING
            CALI_MARK_BEGIN("pack_B");
#endif
#if USE_VECTORIZED_PACKING
            int N_main = current_N & ~15;
            for (int jj = 0; jj < N_main; jj += 16) {
                double *__restrict B_out = &B_pack[(jj / 16) * panel_stride_B];

                for (int kk = 0; kk < current_K; ++kk) {
                    const double *__restrict B_in_ptr = &B[(j + jj) * (ptrdiff_t)ldb + k + kk];
                    double *__restrict B_out_ptr = &B_out[kk * 16];

                    ptrdiff_t ldb_sz = (ptrdiff_t)ldb;

                    B_out_ptr[0] = B_in_ptr[0];
                    B_out_ptr[1] = B_in_ptr[1 * ldb_sz];
                    B_out_ptr[2] = B_in_ptr[2 * ldb_sz];
                    B_out_ptr[3] = B_in_ptr[3 * ldb_sz];
                    B_out_ptr[4] = B_in_ptr[4 * ldb_sz];
                    B_out_ptr[5] = B_in_ptr[5 * ldb_sz];
                    B_out_ptr[6] = B_in_ptr[6 * ldb_sz];
                    B_out_ptr[7] = B_in_ptr[7 * ldb_sz];
                    B_out_ptr[8] = B_in_ptr[8 * ldb_sz];
                    B_out_ptr[9] = B_in_ptr[9 * ldb_sz];
                    B_out_ptr[10] = B_in_ptr[10 * ldb_sz];
                    B_out_ptr[11] = B_in_ptr[11 * ldb_sz];
                    B_out_ptr[12] = B_in_ptr[12 * ldb_sz];
                    B_out_ptr[13] = B_in_ptr[13 * ldb_sz];
                    B_out_ptr[14] = B_in_ptr[14 * ldb_sz];
                    B_out_ptr[15] = B_in_ptr[15 * ldb_sz];
                }
            }
            if (N_main < current_N) {
                int jj = N_main;
                double *__restrict B_out = &B_pack[(jj / 16) * panel_stride_B];
                for (int bc = 0; bc < 16; ++bc) {
                    if (jj + bc < current_N) {
                        const double *__restrict B_col = &B[(j + jj + bc) * (ptrdiff_t)ldb + k];
                        for (int kk = 0; kk < current_K; ++kk) {
                            B_out[kk * 16 + bc] = B_col[kk];
                        }
                    }
                    else {
                        for (int kk = 0; kk < current_K; ++kk) {
                            B_out[kk * 16 + bc] = 0.0;
                        }
                    }
                }
            }
#else
            int N_main = current_N & ~15;
            for (int jj = 0; jj < N_main; jj += 16) {
                double *__restrict B_out = &B_pack[(jj / 16) * panel_stride_B];
                for (int bc = 0; bc < 16; ++bc) {
                    const double *__restrict B_col = &B[(j + jj + bc) * (ptrdiff_t)ldb + k];
                    for (int kk = 0; kk < current_K; ++kk) {
                        B_out[kk * 16 + bc] = B_col[kk];
                    }
                }
            }
            if (N_main < current_N) {
                int jj = N_main;
                double *__restrict B_out = &B_pack[(jj / 16) * panel_stride_B];
                for (int bc = 0; bc < 16; ++bc) {
                    if (jj + bc < current_N) {
                        const double *__restrict B_col = &B[(j + jj + bc) * (ptrdiff_t)ldb + k];
                        for (int kk = 0; kk < current_K; ++kk) {
                            B_out[kk * 16 + bc] = B_col[kk];
                        }
                    }
                    else {
                        for (int kk = 0; kk < current_K; ++kk) {
                            B_out[kk * 16 + bc] = 0.0;
                        }
                    }
                }
            }
#endif
#if PROFILING
            CALI_MARK_END("pack_B");
#endif

            for (int i = 0; i < M; i += MC) {
#if PROFILING
                CALI_MARK_BEGIN("loop_i");
#endif
                int current_M = stdmin(MC, M - i);
                int M_pad = (current_M + 15) & ~15;
                ptrdiff_t panel_stride_A = 16 * (ptrdiff_t)current_K;

#if PROFILING
                CALI_MARK_BEGIN("pack_A");
#endif
#if USE_VECTORIZED_PACKING
                float64x2_t valpha = vdupq_n_f64(alpha);

                int M_main = current_M & ~15;
                for (int ii = 0; ii < M_main; ii += 16) {
                    double *__restrict A_out = &A_pack[(ii / 16) * panel_stride_A];

                    for (int kk = 0; kk < current_K; ++kk) {
                        const double *__restrict A_col = &A[(k + kk) * (ptrdiff_t)lda + i + ii];
                        double *__restrict A_out_ptr = &A_out[kk * 16];

                        // Alpha is mathematically confined here (M*K operations exact)
                        float64x2_t v0 = vmulq_f64(vld1q_f64(&A_col[0]), valpha);
                        float64x2_t v1 = vmulq_f64(vld1q_f64(&A_col[2]), valpha);
                        float64x2_t v2 = vmulq_f64(vld1q_f64(&A_col[4]), valpha);
                        float64x2_t v3 = vmulq_f64(vld1q_f64(&A_col[6]), valpha);
                        float64x2_t v4 = vmulq_f64(vld1q_f64(&A_col[8]), valpha);
                        float64x2_t v5 = vmulq_f64(vld1q_f64(&A_col[10]), valpha);
                        float64x2_t v6 = vmulq_f64(vld1q_f64(&A_col[12]), valpha);
                        float64x2_t v7 = vmulq_f64(vld1q_f64(&A_col[14]), valpha);

                        vst1q_f64(&A_out_ptr[0], v0);
                        vst1q_f64(&A_out_ptr[2], v1);
                        vst1q_f64(&A_out_ptr[4], v2);
                        vst1q_f64(&A_out_ptr[6], v3);
                        vst1q_f64(&A_out_ptr[8], v4);
                        vst1q_f64(&A_out_ptr[10], v5);
                        vst1q_f64(&A_out_ptr[12], v6);
                        vst1q_f64(&A_out_ptr[14], v7);
                    }
                }
                if (M_main < current_M) {
                    int ii = M_main;
                    double *__restrict A_out = &A_pack[(ii / 16) * panel_stride_A];
                    for (int kk = 0; kk < current_K; ++kk) {
                        const double *__restrict A_col = &A[(k + kk) * (ptrdiff_t)lda + i + ii];
                        int br = 0;
                        for (; br < current_M - ii; ++br) {
                            A_out[kk * 16 + br] = A_col[br] * alpha;
                        }
                        for (; br < 16; ++br) {
                            A_out[kk * 16 + br] = 0.0;
                        }
                    }
                }
#else
                int M_main = current_M & ~15;
                for (int ii = 0; ii < M_main; ii += 16) {
                    double *__restrict A_out = &A_pack[(ii / 16) * panel_stride_A];
                    for (int kk = 0; kk < current_K; ++kk) {
                        const double *__restrict A_col = &A[(k + kk) * (ptrdiff_t)lda + i + ii];
                        for (int br = 0; br < 16; ++br) {
                            A_out[kk * 16 + br] = A_col[br] * alpha;
                        }
                    }
                }
                if (M_main < current_M) {
                    int ii = M_main;
                    double *__restrict A_out = &A_pack[(ii / 16) * panel_stride_A];
                    for (int kk = 0; kk < current_K; ++kk) {
                        const double *__restrict A_col = &A[(k + kk) * (ptrdiff_t)lda + i + ii];
                        int br = 0;
                        for (; br < current_M - ii; ++br) {
                            A_out[kk * 16 + br] = A_col[br] * alpha;
                        }
                        for (; br < 16; ++br) {
                            A_out[kk * 16 + br] = 0.0;
                        }
                    }
                }
#endif
#if PROFILING
                CALI_MARK_END("pack_A");
#endif

                for (int jj = 0; jj < N_pad; jj += 16) {
#if PROFILING
                    CALI_MARK_BEGIN("loop_jj");
#endif
                    int current_N_block = stdmin(16, current_N - jj);
                    double *B_ptr = &B_pack[(jj / 16) * panel_stride_B];

                    for (int ii = 0; ii < M_pad; ii += 16) {
#if PROFILING
                        CALI_MARK_BEGIN("loop_ii");
#endif
                        int current_M_block = stdmin(16, current_M - ii);
                        double *A_ptr = &A_pack[(ii / 16) * panel_stride_A];
                        double *C_ptr = &C[(i + ii) + (j + jj) * ldc];
#if PROFILING
                        CALI_MARK_BEGIN("kernel");
#endif
                        if (current_M_block == 16 && current_N_block == 16) {
                            dgemm_sme_compute_16x16_tile(current_K, A_ptr, B_ptr, C_ptr, ldc, beta_mode, &beta);
                        }
                        else {
                            alignas(256) double C_buffer[256];
                            // We strictly only require initialization if the kernel natively relies on reading C.
                            // If beta == 0 (beta_mode == 0), the kernel entirely skips reading the buffer and writes over it.
                            if (beta_mode != 0) {
                                for (int idx = 0; idx < 256; ++idx) {
                                    C_buffer[idx] = 0.0;
                                }
                                for (int bc = 0; bc < current_N_block; ++bc) {
                                    for (int br = 0; br < current_M_block; ++br) {
                                        C_buffer[br + bc * 16] = C_ptr[br + bc * ldc];
                                    }
                                }
                            }
                            dgemm_sme_compute_16x16_tile(current_K, A_ptr, B_ptr, C_buffer, 16, beta_mode, &beta);
                            for (int bc = 0; bc < current_N_block; ++bc) {
                                for (int br = 0; br < current_M_block; ++br) {
                                    C_ptr[br + bc * ldc] = C_buffer[br + bc * 16];
                                }
                            }
                        }
#if PROFILING
                        CALI_MARK_END("kernel");
#endif
#if PROFILING
                        CALI_MARK_END("loop_ii");
#endif
                    }
#if PROFILING
                    CALI_MARK_END("loop_jj");
#endif
                }
#if PROFILING
                CALI_MARK_END("loop_i");
#endif
            }
#if PROFILING
            CALI_MARK_END("loop_k");
#endif
        }
#if PROFILING
        CALI_MARK_END("loop_j");
#endif
    }
}

// ========================================================================
// VERSION 1: C = alpha * A^T * B + beta * C (TN)
// ========================================================================
static void dgemm_sme_TN(int M, int N, int K, double alpha, const double *A, int lda, const double *B, int ldb, double beta, double *C, int ldc)
{
    if (alpha == 0.0 || K == 0) {
        if (beta != 1.0) {
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < M; ++i) {
                    C[i + j * ldc] = (beta == 0.0) ? 0.0 : C[i + j * ldc] * beta;
                }
            }
        }
        return;
    }

    alignas(256) /*thread_local*/ static double A_pack[MC * KC];
    alignas(256) /*thread_local*/ static double B_pack[KC * NC];

#pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < N; j += NC) {
        int current_N = stdmin(NC, N - j);
        int N_pad = (current_N + 15) & ~15;

        for (int k = 0; k < K; k += KC) {
            int current_K = stdmin(KC, K - k);
            ptrdiff_t panel_stride_B = 16 * (ptrdiff_t)current_K;

            int beta_mode = 2;
            if (k == 0) {
                if (beta == 0.0) {
                    beta_mode = 0;
                }
                else if (beta != 1.0) {
                    beta_mode = 1;
                }
            }

#if USE_VECTORIZED_PACKING
            int N_main = current_N & ~15;
            for (int jj = 0; jj < N_main; jj += 16) {
                double *__restrict B_out = &B_pack[(jj / 16) * panel_stride_B];

                for (int kk = 0; kk < current_K; ++kk) {
                    const double *__restrict B_in_ptr = &B[(j + jj) * (ptrdiff_t)ldb + k + kk];
                    double *__restrict B_out_ptr = &B_out[kk * 16];

                    ptrdiff_t ldb_sz = (ptrdiff_t)ldb;

                    B_out_ptr[0] = B_in_ptr[0];
                    B_out_ptr[1] = B_in_ptr[1 * ldb_sz];
                    B_out_ptr[2] = B_in_ptr[2 * ldb_sz];
                    B_out_ptr[3] = B_in_ptr[3 * ldb_sz];
                    B_out_ptr[4] = B_in_ptr[4 * ldb_sz];
                    B_out_ptr[5] = B_in_ptr[5 * ldb_sz];
                    B_out_ptr[6] = B_in_ptr[6 * ldb_sz];
                    B_out_ptr[7] = B_in_ptr[7 * ldb_sz];
                    B_out_ptr[8] = B_in_ptr[8 * ldb_sz];
                    B_out_ptr[9] = B_in_ptr[9 * ldb_sz];
                    B_out_ptr[10] = B_in_ptr[10 * ldb_sz];
                    B_out_ptr[11] = B_in_ptr[11 * ldb_sz];
                    B_out_ptr[12] = B_in_ptr[12 * ldb_sz];
                    B_out_ptr[13] = B_in_ptr[13 * ldb_sz];
                    B_out_ptr[14] = B_in_ptr[14 * ldb_sz];
                    B_out_ptr[15] = B_in_ptr[15 * ldb_sz];
                }
            }
            if (N_main < current_N) {
                int jj = N_main;
                double *__restrict B_out = &B_pack[(jj / 16) * panel_stride_B];
                for (int bc = 0; bc < 16; ++bc) {
                    if (jj + bc < current_N) {
                        const double *__restrict B_col = &B[(j + jj + bc) * (ptrdiff_t)ldb + k];
                        for (int kk = 0; kk < current_K; ++kk) {
                            B_out[kk * 16 + bc] = B_col[kk];
                        }
                    }
                    else {
                        for (int kk = 0; kk < current_K; ++kk) {
                            B_out[kk * 16 + bc] = 0.0;
                        }
                    }
                }
            }
#else
            int N_main = current_N & ~15;
            for (int jj = 0; jj < N_main; jj += 16) {
                double *__restrict B_out = &B_pack[(jj / 16) * panel_stride_B];
                for (int bc = 0; bc < 16; ++bc) {
                    const double *__restrict B_col = &B[(j + jj + bc) * (ptrdiff_t)ldb + k];
                    for (int kk = 0; kk < current_K; ++kk) {
                        B_out[kk * 16 + bc] = B_col[kk];
                    }
                }
            }
            if (N_main < current_N) {
                int jj = N_main;
                double *__restrict B_out = &B_pack[(jj / 16) * panel_stride_B];
                for (int bc = 0; bc < 16; ++bc) {
                    if (jj + bc < current_N) {
                        const double *__restrict B_col = &B[(j + jj + bc) * (ptrdiff_t)ldb + k];
                        for (int kk = 0; kk < current_K; ++kk) {
                            B_out[kk * 16 + bc] = B_col[kk];
                        }
                    }
                    else {
                        for (int kk = 0; kk < current_K; ++kk) {
                            B_out[kk * 16 + bc] = 0.0;
                        }
                    }
                }
            }
#endif

            for (int i = 0; i < M; i += MC) {
                int current_M = stdmin(MC, M - i);
                int M_pad = (current_M + 15) & ~15;
                ptrdiff_t panel_stride_A = 16 * (ptrdiff_t)current_K;

#if USE_VECTORIZED_PACKING
                int M_main = current_M & ~15;
                for (int ii = 0; ii < M_main; ii += 16) {
                    double *__restrict A_out = &A_pack[(ii / 16) * panel_stride_A];
                    for (int kk = 0; kk < current_K; ++kk) {
                        const double *__restrict A_in_ptr = &A[(i + ii) * (ptrdiff_t)lda + k + kk];
                        double *__restrict A_out_ptr = &A_out[kk * 16];
                        ptrdiff_t lda_sz = (ptrdiff_t)lda;

                        A_out_ptr[0] = A_in_ptr[0] * alpha;
                        A_out_ptr[1] = A_in_ptr[1 * lda_sz] * alpha;
                        A_out_ptr[2] = A_in_ptr[2 * lda_sz] * alpha;
                        A_out_ptr[3] = A_in_ptr[3 * lda_sz] * alpha;
                        A_out_ptr[4] = A_in_ptr[4 * lda_sz] * alpha;
                        A_out_ptr[5] = A_in_ptr[5 * lda_sz] * alpha;
                        A_out_ptr[6] = A_in_ptr[6 * lda_sz] * alpha;
                        A_out_ptr[7] = A_in_ptr[7 * lda_sz] * alpha;
                        A_out_ptr[8] = A_in_ptr[8 * lda_sz] * alpha;
                        A_out_ptr[9] = A_in_ptr[9 * lda_sz] * alpha;
                        A_out_ptr[10] = A_in_ptr[10 * lda_sz] * alpha;
                        A_out_ptr[11] = A_in_ptr[11 * lda_sz] * alpha;
                        A_out_ptr[12] = A_in_ptr[12 * lda_sz] * alpha;
                        A_out_ptr[13] = A_in_ptr[13 * lda_sz] * alpha;
                        A_out_ptr[14] = A_in_ptr[14 * lda_sz] * alpha;
                        A_out_ptr[15] = A_in_ptr[15 * lda_sz] * alpha;
                    }
                }
                if (M_main < current_M) {
                    int ii = M_main;
                    double *__restrict A_out = &A_pack[(ii / 16) * panel_stride_A];
                    for (int br = 0; br < 16; ++br) {
                        if (ii + br < current_M) {
                            const double *__restrict A_row = &A[(i + ii + br) * (ptrdiff_t)lda + k];
                            for (int kk = 0; kk < current_K; ++kk) {
                                A_out[kk * 16 + br] = A_row[kk] * alpha;
                            }
                        }
                        else {
                            for (int kk = 0; kk < current_K; ++kk) {
                                A_out[kk * 16 + br] = 0.0;
                            }
                        }
                    }
                }
#else
                int M_main = current_M & ~15;
                for (int ii = 0; ii < M_main; ii += 16) {
                    double *__restrict A_out = &A_pack[(ii / 16) * panel_stride_A];
                    for (int br = 0; br < 16; ++br) {
                        const double *__restrict A_row = &A[(i + ii + br) * (ptrdiff_t)lda + k];
                        for (int kk = 0; kk < current_K; ++kk) {
                            A_out[kk * 16 + br] = A_row[kk] * alpha;
                        }
                    }
                }
                if (M_main < current_M) {
                    int ii = M_main;
                    double *__restrict A_out = &A_pack[(ii / 16) * panel_stride_A];
                    for (int br = 0; br < 16; ++br) {
                        if (ii + br < current_M) {
                            const double *__restrict A_row = &A[(i + ii + br) * (ptrdiff_t)lda + k];
                            for (int kk = 0; kk < current_K; ++kk) {
                                A_out[kk * 16 + br] = A_row[kk] * alpha;
                            }
                        }
                        else {
                            for (int kk = 0; kk < current_K; ++kk) {
                                A_out[kk * 16 + br] = 0.0;
                            }
                        }
                    }
                }
#endif

                for (int jj = 0; jj < N_pad; jj += 16) {
                    int current_N_block = stdmin(16, current_N - jj);
                    double *B_ptr = &B_pack[(jj / 16) * panel_stride_B];

                    for (int ii = 0; ii < M_pad; ii += 16) {
                        int current_M_block = stdmin(16, current_M - ii);
                        double *A_ptr = &A_pack[(ii / 16) * panel_stride_A];
                        double *C_ptr = &C[(i + ii) + (j + jj) * ldc];

                        if (current_M_block == 16 && current_N_block == 16) {
                            dgemm_sme_compute_16x16_tile(current_K, A_ptr, B_ptr, C_ptr, ldc, beta_mode, &beta);
                        }
                        else {
                            alignas(256) double C_buffer[256];
                            if (beta_mode != 0) {
                                for (int idx = 0; idx < 256; ++idx) {
                                    C_buffer[idx] = 0.0;
                                }
                                for (int bc = 0; bc < current_N_block; ++bc) {
                                    for (int br = 0; br < current_M_block; ++br) {
                                        C_buffer[br + bc * 16] = C_ptr[br + bc * ldc];
                                    }
                                }
                            }
                            dgemm_sme_compute_16x16_tile(current_K, A_ptr, B_ptr, C_buffer, 16, beta_mode, &beta);
                            for (int bc = 0; bc < current_N_block; ++bc) {
                                for (int br = 0; br < current_M_block; ++br) {
                                    C_ptr[br + bc * ldc] = C_buffer[br + bc * 16];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// ========================================================================
// VERSION 2: C = alpha * A * B^T + beta * C (NT)
// ========================================================================
static void dgemm_sme_NT(int M, int N, int K, double alpha, const double *A, int lda, const double *B, int ldb, double beta, double *C, int ldc)
{
    if (alpha == 0.0 || K == 0) {
        if (beta != 1.0) {
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < M; ++i) {
                    C[i + j * ldc] = (beta == 0.0) ? 0.0 : C[i + j * ldc] * beta;
                }
            }
        }
        return;
    }

    alignas(256) /*thread_local*/ static double A_pack[MC * KC];
    alignas(256) /*thread_local*/ static double B_pack[KC * NC];

#pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < N; j += NC) {
        int current_N = stdmin(NC, N - j);
        int N_pad = (current_N + 15) & ~15;

        for (int k = 0; k < K; k += KC) {
            int current_K = stdmin(KC, K - k);
            ptrdiff_t panel_stride_B = 16 * (ptrdiff_t)current_K;

            int beta_mode = 2;
            if (k == 0) {
                if (beta == 0.0) {
                    beta_mode = 0;
                }
                else if (beta != 1.0) {
                    beta_mode = 1;
                }
            }

#if USE_VECTORIZED_PACKING
            int N_main = current_N & ~15;
            for (int jj = 0; jj < N_main; jj += 16) {
                double *__restrict B_out = &B_pack[(jj / 16) * panel_stride_B];

                for (int kk = 0; kk < current_K; ++kk) {
                    const double *__restrict B_col = &B[(k + kk) * (ptrdiff_t)ldb + j + jj];
                    double *__restrict B_out_ptr = &B_out[kk * 16];

                    float64x2_t v0 = vld1q_f64(&B_col[0]);
                    float64x2_t v1 = vld1q_f64(&B_col[2]);
                    float64x2_t v2 = vld1q_f64(&B_col[4]);
                    float64x2_t v3 = vld1q_f64(&B_col[6]);
                    float64x2_t v4 = vld1q_f64(&B_col[8]);
                    float64x2_t v5 = vld1q_f64(&B_col[10]);
                    float64x2_t v6 = vld1q_f64(&B_col[12]);
                    float64x2_t v7 = vld1q_f64(&B_col[14]);

                    vst1q_f64(&B_out_ptr[0], v0);
                    vst1q_f64(&B_out_ptr[2], v1);
                    vst1q_f64(&B_out_ptr[4], v2);
                    vst1q_f64(&B_out_ptr[6], v3);
                    vst1q_f64(&B_out_ptr[8], v4);
                    vst1q_f64(&B_out_ptr[10], v5);
                    vst1q_f64(&B_out_ptr[12], v6);
                    vst1q_f64(&B_out_ptr[14], v7);
                }
            }
            if (N_main < current_N) {
                int jj = N_main;
                double *__restrict B_out = &B_pack[(jj / 16) * panel_stride_B];
                for (int kk = 0; kk < current_K; ++kk) {
                    const double *__restrict B_row = &B[(k + kk) * (ptrdiff_t)ldb + j + jj];
                    int bc = 0;
                    for (; bc < current_N - jj; ++bc) {
                        B_out[kk * 16 + bc] = B_row[bc];
                    }
                    for (; bc < 16; ++bc) {
                        B_out[kk * 16 + bc] = 0.0;
                    }
                }
            }
#else
            int N_main = current_N & ~15;
            for (int jj = 0; jj < N_main; jj += 16) {
                double *__restrict B_out = &B_pack[(jj / 16) * panel_stride_B];
                for (int kk = 0; kk < current_K; ++kk) {
                    const double *__restrict B_row = &B[(k + kk) * (ptrdiff_t)ldb + j + jj];
                    for (int bc = 0; bc < 16; ++bc) {
                        B_out[kk * 16 + bc] = B_row[bc];
                    }
                }
            }
            if (N_main < current_N) {
                int jj = N_main;
                double *__restrict B_out = &B_pack[(jj / 16) * panel_stride_B];
                for (int kk = 0; kk < current_K; ++kk) {
                    const double *__restrict B_row = &B[(k + kk) * (ptrdiff_t)ldb + j + jj];
                    int bc = 0;
                    for (; bc < current_N - jj; ++bc) {
                        B_out[kk * 16 + bc] = B_row[bc];
                    }
                    for (; bc < 16; ++bc) {
                        B_out[kk * 16 + bc] = 0.0;
                    }
                }
            }
#endif

            for (int i = 0; i < M; i += MC) {
                int current_M = stdmin(MC, M - i);
                int M_pad = (current_M + 15) & ~15;
                ptrdiff_t panel_stride_A = 16 * (ptrdiff_t)current_K;

#if USE_VECTORIZED_PACKING
                float64x2_t valpha = vdupq_n_f64(alpha);

                int M_main = current_M & ~15;
                for (int ii = 0; ii < M_main; ii += 16) {
                    double *__restrict A_out = &A_pack[(ii / 16) * panel_stride_A];

                    for (int kk = 0; kk < current_K; ++kk) {
                        const double *__restrict A_col = &A[(k + kk) * (ptrdiff_t)lda + i + ii];
                        double *__restrict A_out_ptr = &A_out[kk * 16];

                        float64x2_t v0 = vmulq_f64(vld1q_f64(&A_col[0]), valpha);
                        float64x2_t v1 = vmulq_f64(vld1q_f64(&A_col[2]), valpha);
                        float64x2_t v2 = vmulq_f64(vld1q_f64(&A_col[4]), valpha);
                        float64x2_t v3 = vmulq_f64(vld1q_f64(&A_col[6]), valpha);
                        float64x2_t v4 = vmulq_f64(vld1q_f64(&A_col[8]), valpha);
                        float64x2_t v5 = vmulq_f64(vld1q_f64(&A_col[10]), valpha);
                        float64x2_t v6 = vmulq_f64(vld1q_f64(&A_col[12]), valpha);
                        float64x2_t v7 = vmulq_f64(vld1q_f64(&A_col[14]), valpha);

                        vst1q_f64(&A_out_ptr[0], v0);
                        vst1q_f64(&A_out_ptr[2], v1);
                        vst1q_f64(&A_out_ptr[4], v2);
                        vst1q_f64(&A_out_ptr[6], v3);
                        vst1q_f64(&A_out_ptr[8], v4);
                        vst1q_f64(&A_out_ptr[10], v5);
                        vst1q_f64(&A_out_ptr[12], v6);
                        vst1q_f64(&A_out_ptr[14], v7);
                    }
                }
                if (M_main < current_M) {
                    int ii = M_main;
                    double *__restrict A_out = &A_pack[(ii / 16) * panel_stride_A];
                    for (int kk = 0; kk < current_K; ++kk) {
                        const double *__restrict A_col = &A[(k + kk) * (ptrdiff_t)lda + i + ii];
                        int br = 0;
                        for (; br < current_M - ii; ++br) {
                            A_out[kk * 16 + br] = A_col[br] * alpha;
                        }
                        for (; br < 16; ++br) {
                            A_out[kk * 16 + br] = 0.0;
                        }
                    }
                }
#else
                int M_main = current_M & ~15;
                for (int ii = 0; ii < M_main; ii += 16) {
                    double *__restrict A_out = &A_pack[(ii / 16) * panel_stride_A];
                    for (int kk = 0; kk < current_K; ++kk) {
                        const double *__restrict A_col = &A[(k + kk) * (ptrdiff_t)lda + i + ii];
                        for (int br = 0; br < 16; ++br) {
                            A_out[kk * 16 + br] = A_col[br] * alpha;
                        }
                    }
                }
                if (M_main < current_M) {
                    int ii = M_main;
                    double *__restrict A_out = &A_pack[(ii / 16) * panel_stride_A];
                    for (int kk = 0; kk < current_K; ++kk) {
                        const double *__restrict A_col = &A[(k + kk) * (ptrdiff_t)lda + i + ii];
                        int br = 0;
                        for (; br < current_M - ii; ++br) {
                            A_out[kk * 16 + br] = A_col[br] * alpha;
                        }
                        for (; br < 16; ++br) {
                            A_out[kk * 16 + br] = 0.0;
                        }
                    }
                }
#endif

                for (int jj = 0; jj < N_pad; jj += 16) {
                    int current_N_block = stdmin(16, current_N - jj);
                    double *B_ptr = &B_pack[(jj / 16) * panel_stride_B];

                    for (int ii = 0; ii < M_pad; ii += 16) {
                        int current_M_block = stdmin(16, current_M - ii);
                        double *A_ptr = &A_pack[(ii / 16) * panel_stride_A];
                        double *C_ptr = &C[(i + ii) + (j + jj) * ldc];

                        if (current_M_block == 16 && current_N_block == 16) {
                            dgemm_sme_compute_16x16_tile(current_K, A_ptr, B_ptr, C_ptr, ldc, beta_mode, &beta);
                        }
                        else {
                            alignas(256) double C_buffer[256];
                            if (beta_mode != 0) {
                                for (int idx = 0; idx < 256; ++idx) {
                                    C_buffer[idx] = 0.0;
                                }
                                for (int bc = 0; bc < current_N_block; ++bc) {
                                    for (int br = 0; br < current_M_block; ++br) {
                                        C_buffer[br + bc * 16] = C_ptr[br + bc * ldc];
                                    }
                                }
                            }
                            dgemm_sme_compute_16x16_tile(current_K, A_ptr, B_ptr, C_buffer, 16, beta_mode, &beta);
                            for (int bc = 0; bc < current_N_block; ++bc) {
                                for (int br = 0; br < current_M_block; ++br) {
                                    C_ptr[br + bc * ldc] = C_buffer[br + bc * 16];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// ========================================================================
// VERSION 3: C = alpha * A^T * B^T + beta * C (TT)
// ========================================================================
static void dgemm_sme_TT(int M, int N, int K, double alpha, const double *A, int lda, const double *B, int ldb, double beta, double *C, int ldc)
{
    if (alpha == 0.0 || K == 0) {
        if (beta != 1.0) {
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < M; ++i) {
                    C[i + j * ldc] = (beta == 0.0) ? 0.0 : C[i + j * ldc] * beta;
                }
            }
        }
        return;
    }

    alignas(256) /*thread_local*/ static double A_pack[MC * KC];
    alignas(256) /*thread_local*/ static double B_pack[KC * NC];

#pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < N; j += NC) {
        int current_N = stdmin(NC, N - j);
        int N_pad = (current_N + 15) & ~15;

        for (int k = 0; k < K; k += KC) {
            int current_K = stdmin(KC, K - k);
            ptrdiff_t panel_stride_B = 16 * (ptrdiff_t)current_K;

            int beta_mode = 2;
            if (k == 0) {
                if (beta == 0.0) {
                    beta_mode = 0;
                }
                else if (beta != 1.0) {
                    beta_mode = 1;
                }
            }

#if USE_VECTORIZED_PACKING
            int N_main = current_N & ~15;
            for (int jj = 0; jj < N_main; jj += 16) {
                double *__restrict B_out = &B_pack[(jj / 16) * panel_stride_B];

                for (int kk = 0; kk < current_K; ++kk) {
                    const double *__restrict B_col = &B[(k + kk) * (ptrdiff_t)ldb + j + jj];
                    double *__restrict B_out_ptr = &B_out[kk * 16];

                    vst1q_f64(&B_out_ptr[0], vld1q_f64(&B_col[0]));
                    vst1q_f64(&B_out_ptr[2], vld1q_f64(&B_col[2]));
                    vst1q_f64(&B_out_ptr[4], vld1q_f64(&B_col[4]));
                    vst1q_f64(&B_out_ptr[6], vld1q_f64(&B_col[6]));
                    vst1q_f64(&B_out_ptr[8], vld1q_f64(&B_col[8]));
                    vst1q_f64(&B_out_ptr[10], vld1q_f64(&B_col[10]));
                    vst1q_f64(&B_out_ptr[12], vld1q_f64(&B_col[12]));
                    vst1q_f64(&B_out_ptr[14], vld1q_f64(&B_col[14]));
                }
            }
            if (N_main < current_N) {
                int jj = N_main;
                double *__restrict B_out = &B_pack[(jj / 16) * panel_stride_B];
                for (int kk = 0; kk < current_K; ++kk) {
                    const double *__restrict B_row = &B[(k + kk) * (ptrdiff_t)ldb + j + jj];
                    int bc = 0;
                    for (; bc < current_N - jj; ++bc) {
                        B_out[kk * 16 + bc] = B_row[bc];
                    }
                    for (; bc < 16; ++bc) {
                        B_out[kk * 16 + bc] = 0.0;
                    }
                }
            }
#else
            int N_main = current_N & ~15;
            for (int jj = 0; jj < N_main; jj += 16) {
                double *__restrict B_out = &B_pack[(jj / 16) * panel_stride_B];
                for (int kk = 0; kk < current_K; ++kk) {
                    const double *__restrict B_row = &B[(k + kk) * (ptrdiff_t)ldb + j + jj];
                    for (int bc = 0; bc < 16; ++bc) {
                        B_out[kk * 16 + bc] = B_row[bc];
                    }
                }
            }
            if (N_main < current_N) {
                int jj = N_main;
                double *__restrict B_out = &B_pack[(jj / 16) * panel_stride_B];
                for (int kk = 0; kk < current_K; ++kk) {
                    const double *__restrict B_row = &B[(k + kk) * (ptrdiff_t)ldb + j + jj];
                    int bc = 0;
                    for (; bc < current_N - jj; ++bc) {
                        B_out[kk * 16 + bc] = B_row[bc];
                    }
                    for (; bc < 16; ++bc) {
                        B_out[kk * 16 + bc] = 0.0;
                    }
                }
            }
#endif

            for (int i = 0; i < M; i += MC) {
                int current_M = stdmin(MC, M - i);
                int M_pad = (current_M + 15) & ~15;
                ptrdiff_t panel_stride_A = 16 * (ptrdiff_t)current_K;

#if USE_VECTORIZED_PACKING
                int M_main = current_M & ~15;
                for (int ii = 0; ii < M_main; ii += 16) {
                    double *__restrict A_out = &A_pack[(ii / 16) * panel_stride_A];
                    for (int kk = 0; kk < current_K; ++kk) {
                        const double *__restrict A_in_ptr = &A[(i + ii) * (ptrdiff_t)lda + k + kk];
                        double *__restrict A_out_ptr = &A_out[kk * 16];
                        ptrdiff_t lda_sz = (ptrdiff_t)lda;

                        A_out_ptr[0] = A_in_ptr[0] * alpha;
                        A_out_ptr[1] = A_in_ptr[1 * lda_sz] * alpha;
                        A_out_ptr[2] = A_in_ptr[2 * lda_sz] * alpha;
                        A_out_ptr[3] = A_in_ptr[3 * lda_sz] * alpha;
                        A_out_ptr[4] = A_in_ptr[4 * lda_sz] * alpha;
                        A_out_ptr[5] = A_in_ptr[5 * lda_sz] * alpha;
                        A_out_ptr[6] = A_in_ptr[6 * lda_sz] * alpha;
                        A_out_ptr[7] = A_in_ptr[7 * lda_sz] * alpha;
                        A_out_ptr[8] = A_in_ptr[8 * lda_sz] * alpha;
                        A_out_ptr[9] = A_in_ptr[9 * lda_sz] * alpha;
                        A_out_ptr[10] = A_in_ptr[10 * lda_sz] * alpha;
                        A_out_ptr[11] = A_in_ptr[11 * lda_sz] * alpha;
                        A_out_ptr[12] = A_in_ptr[12 * lda_sz] * alpha;
                        A_out_ptr[13] = A_in_ptr[13 * lda_sz] * alpha;
                        A_out_ptr[14] = A_in_ptr[14 * lda_sz] * alpha;
                        A_out_ptr[15] = A_in_ptr[15 * lda_sz] * alpha;
                    }
                }
                if (M_main < current_M) {
                    int ii = M_main;
                    double *__restrict A_out = &A_pack[(ii / 16) * panel_stride_A];
                    for (int br = 0; br < 16; ++br) {
                        if (ii + br < current_M) {
                            const double *__restrict A_row = &A[(i + ii + br) * (ptrdiff_t)lda + k];
                            for (int kk = 0; kk < current_K; ++kk) {
                                A_out[kk * 16 + br] = A_row[kk] * alpha;
                            }
                        }
                        else {
                            for (int kk = 0; kk < current_K; ++kk) {
                                A_out[kk * 16 + br] = 0.0;
                            }
                        }
                    }
                }
#else
                int M_main = current_M & ~15;
                for (int ii = 0; ii < M_main; ii += 16) {
                    double *__restrict A_out = &A_pack[(ii / 16) * panel_stride_A];
                    for (int br = 0; br < 16; ++br) {
                        const double *__restrict A_row = &A[(i + ii + br) * (ptrdiff_t)lda + k];
                        for (int kk = 0; kk < current_K; ++kk) {
                            A_out[kk * 16 + br] = A_row[kk] * alpha;
                        }
                    }
                }
                if (M_main < current_M) {
                    int ii = M_main;
                    double *__restrict A_out = &A_pack[(ii / 16) * panel_stride_A];
                    for (int br = 0; br < 16; ++br) {
                        if (ii + br < current_M) {
                            const double *__restrict A_row = &A[(i + ii + br) * (ptrdiff_t)lda + k];
                            for (int kk = 0; kk < current_K; ++kk) {
                                A_out[kk * 16 + br] = A_row[kk] * alpha;
                            }
                        }
                        else {
                            for (int kk = 0; kk < current_K; ++kk) {
                                A_out[kk * 16 + br] = 0.0;
                            }
                        }
                    }
                }
#endif

                for (int jj = 0; jj < N_pad; jj += 16) {
                    int current_N_block = stdmin(16, current_N - jj);
                    double *B_ptr = &B_pack[(jj / 16) * panel_stride_B];

                    for (int ii = 0; ii < M_pad; ii += 16) {
                        int current_M_block = stdmin(16, current_M - ii);
                        double *A_ptr = &A_pack[(ii / 16) * panel_stride_A];
                        double *C_ptr = &C[(i + ii) + (j + jj) * ldc];

                        if (current_M_block == 16 && current_N_block == 16) {
                            dgemm_sme_compute_16x16_tile(current_K, A_ptr, B_ptr, C_ptr, ldc, beta_mode, &beta);
                        }
                        else {
                            alignas(256) double C_buffer[256];
                            if (beta_mode != 0) {
                                for (int idx = 0; idx < 256; ++idx) {
                                    C_buffer[idx] = 0.0;
                                }
                                for (int bc = 0; bc < current_N_block; ++bc) {
                                    for (int br = 0; br < current_M_block; ++br) {
                                        C_buffer[br + bc * 16] = C_ptr[br + bc * ldc];
                                    }
                                }
                            }
                            dgemm_sme_compute_16x16_tile(current_K, A_ptr, B_ptr, C_buffer, 16, beta_mode, &beta);
                            for (int bc = 0; bc < current_N_block; ++bc) {
                                for (int br = 0; br < current_M_block; ++br) {
                                    C_ptr[br + bc * ldc] = C_buffer[br + bc * 16];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// ========================================================================
// WRAPPER: BLAS ABI Compatible dgemm
// ========================================================================
void CNAME(const char *transa, const char *transb, const BLASLONG m, const BLASLONG n, const BLASLONG k, const double *alpha, const double *a, const BLASLONG lda, const double *b, const BLASLONG ldb, const double *beta, double *c, const BLASLONG ldc)
{
    bool trans_a = (*transa == 'T' || *transa == 't' || *transa == 'C' || *transa == 'c');
    bool trans_b = (*transb == 'T' || *transb == 't' || *transb == 'C' || *transb == 'c');
    if (!trans_a && !trans_b) {
        dgemm_sme_NN(m, n, k, *alpha, a, lda, b, ldb, *beta, c, ldc);
    }
    else if (trans_a && !trans_b) {
        dgemm_sme_TN(m, n, k, *alpha, a, lda, b, ldb, *beta, c, ldc);
    }
    else if (!trans_a && trans_b) {
        dgemm_sme_NT(m, n, k, *alpha, a, lda, b, ldb, *beta, c, ldc);
    }
    else {
        dgemm_sme_TT(m, n, k, *alpha, a, lda, b, ldb, *beta, c, ldc);
    }
}
