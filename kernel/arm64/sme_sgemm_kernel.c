//#include <algorithm>
#include <arm_neon.h>
//#include <cstddef>
#include <stddef.h>
#include <stdalign.h>
#include <stdbool.h>
#include "common.h"
#ifndef stdmin
#define stdmin(a,b)   (a>b? b:a)
#endif

#define MC 512
#define KC 1024
#define NC 2048

static inline void sgemm_sme_compute_32x32_tile(blasint current_K, const float *A_ptr, const float *B_ptr, float *C_ptr, size_t ldc, blasint beta_mode, const float *beta_ptr)
{
    size_t ldc_bytes = ldc * sizeof(float);

    asm volatile("smstart\n\t"
		 "ptrue p0.s\n\t"
                 "cmp %w[beta_mode], #0\n\t"
                 "b.eq 10f\n\t"
                 "cmp %w[beta_mode], #1\n\t"
                 "b.eq 11f\n\t"
                 "b 12f\n\t"

                 "10:\n\t"
                 "zero {za}\n\t"
                 "b 19f\n\t"

                 "11:\n\t"
                 "ld1rw z31.s, p0/z, [%[beta_ptr]]\n\t"
                 "mov w12, #0\n\t"
                 "mov x13, %[c]\n\t"
                 "110:\n\t"
                 "ld1w z0.s, p0/z, [x13]\n\t"
                 "add x14, x13, #64\n\t"
                 "ld1w z1.s, p0/z, [x14]\n\t"
                 "fmul z0.s, p0/m, z0.s, z31.s\n\t"
                 "fmul z1.s, p0/m, z1.s, z31.s\n\t"
                 "mova za0v.s[w12, 0], p0/m, z0.s\n\t"
                 "mova za2v.s[w12, 0], p0/m, z1.s\n\t" // FIXED: za2 (Bottom-Left)
                 "add x13, x13, %[ldc_bytes]\n\t"
                 "add w12, w12, #1\n\t"
                 "cmp w12, #16\n\t"
                 "b.ne 110b\n\t"

                 "mov w15, #0\n\t"
                 "111:\n\t"
                 "ld1w z0.s, p0/z, [x13]\n\t"
                 "add x14, x13, #64\n\t"
                 "ld1w z1.s, p0/z, [x14]\n\t"
                 "fmul z0.s, p0/m, z0.s, z31.s\n\t"
                 "fmul z1.s, p0/m, z1.s, z31.s\n\t"
                 "mova za1v.s[w15, 0], p0/m, z0.s\n\t" // FIXED: za1 (Top-Right)
                 "mova za3v.s[w15, 0], p0/m, z1.s\n\t"
                 "add x13, x13, %[ldc_bytes]\n\t"
                 "add w15, w15, #1\n\t"
                 "cmp w15, #16\n\t"
                 "b.ne 111b\n\t"
                 "b 19f\n\t"

                 "12:\n\t"
                 "mov w12, #0\n\t"
                 "mov x13, %[c]\n\t"
                 "100:\n\t"
                 "ld1w {za0v.s[w12, 0]}, p0/z, [x13]\n\t"
                 "add x14, x13, #64\n\t"
                 "ld1w {za2v.s[w12, 0]}, p0/z, [x14]\n\t" // FIXED: za2
                 "add x13, x13, %[ldc_bytes]\n\t"
                 "add w12, w12, #1\n\t"
                 "cmp w12, #16\n\t"
                 "b.ne 100b\n\t"

                 "mov w15, #0\n\t"
                 "101:\n\t"
                 "ld1w {za1v.s[w15, 0]}, p0/z, [x13]\n\t" // FIXED: za1
                 "add x14, x13, #64\n\t"
                 "ld1w {za3v.s[w15, 0]}, p0/z, [x14]\n\t"
                 "add x13, x13, %[ldc_bytes]\n\t"
                 "add w15, w15, #1\n\t"
                 "cmp w15, #16\n\t"
                 "b.ne 101b\n\t"

                 "19:\n\t"
                 "mov w10, %w[k]\n\t"
                 "cbz w10, 3f\n\t"
                 "11:\n\t"
                 "ld1w z0.s, p0/z, [%[a], #0, mul vl]\n\t"
                 "ld1w z1.s, p0/z, [%[a], #1, mul vl]\n\t"
                 "ld1w z2.s, p0/z, [%[b], #0, mul vl]\n\t"
                 "ld1w z3.s, p0/z, [%[b], #1, mul vl]\n\t"
                 "fmopa za0.s, p0/m, p0/m, z0.s, z2.s\n\t"
                 "fmopa za1.s, p0/m, p0/m, z0.s, z3.s\n\t"
                 "fmopa za2.s, p0/m, p0/m, z1.s, z2.s\n\t"
                 "fmopa za3.s, p0/m, p0/m, z1.s, z3.s\n\t"
                 "add %[a], %[a], #128\n\t"
                 "add %[b], %[b], #128\n\t"
                 "subs w10, w10, #1\n\t"
                 "b.ne 11b\n\t"

                 "3:\n\t"
                 "mov w12, #0\n\t"
                 "mov x13, %[c]\n\t"
                 "200:\n\t"
                 "st1w {za0v.s[w12, 0]}, p0, [x13]\n\t"
                 "add x14, x13, #64\n\t"
                 "st1w {za2v.s[w12, 0]}, p0, [x14]\n\t" // FIXED: za2
                 "add x13, x13, %[ldc_bytes]\n\t"
                 "add w12, w12, #1\n\t"
                 "cmp w12, #16\n\t"
                 "b.ne 200b\n\t"

                 "mov w15, #0\n\t"
                 "201:\n\t"
                 "st1w {za1v.s[w15, 0]}, p0, [x13]\n\t" // FIXED: za1
                 "add x14, x13, #64\n\t"
                 "st1w {za3v.s[w15, 0]}, p0, [x14]\n\t"
                 "add x13, x13, %[ldc_bytes]\n\t"
                 "add w15, w15, #1\n\t"
                 "cmp w15, #16\n\t"
                 "b.ne 201b\n\t"
		 "smstop\n\t"
		"msr fpsr, xzr\n\t"
    : [a] "+&r"(A_ptr), [b] "+&r"(B_ptr)
    : [k] "r"(current_K), [c] "r"(C_ptr), [ldc_bytes] "r"(ldc_bytes), [beta_mode] "r"(beta_mode), [beta_ptr] "r"(beta_ptr)
    : "x10", "x12", "x13", "x14", "x15", "memory", "cc", 
                         "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7",
                         "v8", "v9", "v10", "v11", "v12", "v13", "v14", "v15",
                         "v16", "v17", "v18", "v19", "v20", "v21", "v22", "v23",
                         "v24", "v25", "v26", "v27", "v28", "v29", "v30", "v31",
                         "z0", "z1", "z2", "z3", "z4", "z5", "z6", "z7",
                         "z8", "z9", "z10", "z11", "z12", "z13", "z14", "z15",
                         "z16", "z17", "z18", "z19", "z20", "z21", "z22", "z23",
                         "z24", "z25", "z26", "z27", "z28", "z29", "z30", "z31",
      			 "p0", "p1", "p2", "p3", "p4", "p5", "p6", "p7",
                         "p8", "p9", "p10", "p11", "p12", "p13", "p14", "p15","za");


}

static void sgemm_sme_NN(blasint M, blasint N, blasint K, float alpha, const float *A, blasint lda, const float *B, blasint ldb, float beta, float *C, blasint ldc)
{
    if (alpha == 0.0f || K == 0) {
        if (beta != 1.0f) {
            for (blasint j = 0; j < N; ++j) {
                for (blasint i = 0; i < M; ++i) {
                    C[i + j * ldc] = (beta == 0.0f) ? 0.0f : C[i + j * ldc] * beta;
                }
            }
        }
        return;
    }

    alignas(256) /*thread_local*/ static float A_pack[MC * KC];
    alignas(256) /*thread_local*/ static float B_pack[KC * NC];

#pragma omp parallel for schedule(dynamic, 1)
    for (blasint j = 0; j < N; j += NC) {
        blasint current_N = stdmin(NC, N - j);
        blasint N_pad = (current_N + 31) & ~31;
        for (blasint k = 0; k < K; k += KC) {
            blasint current_K = stdmin(KC, K - k);
            size_t panel_stride_B = 32 * (size_t)current_K;
            int beta_mode = (k == 0) ? ((beta == 0.0f) ? 0 : (beta != 1.0f ? 1 : 2)) : 2;

            blasint N_main = current_N & ~31;
            for (blasint jj = 0; jj < N_main; jj += 32) {
                float *__restrict B_out = &B_pack[(jj / 32) * panel_stride_B];
                for (blasint kk = 0; kk < current_K; ++kk) {
                    const float *__restrict B_in = &B[(j + jj) * (size_t)ldb + k + kk];
                    for (int bc = 0; bc < 32; ++bc) {
                        B_out[kk * 32 + bc] = B_in[bc * ldb];
                    }
                }
            }
            if (N_main < current_N) {
                float *__restrict B_out = &B_pack[(N_main / 32) * panel_stride_B];
                for (int bc = 0; bc < 32; ++bc) {
                    if (N_main + bc < current_N) {
                        const float *__restrict B_col = &B[(j + N_main + bc) * (size_t)ldb + k];
                        for (blasint kk = 0; kk < current_K; ++kk) {
                            B_out[kk * 32 + bc] = B_col[kk];
                        }
                    }
                    else {
                        for (blasint kk = 0; kk < current_K; ++kk) {
                            B_out[kk * 32 + bc] = 0.0f;
                        }
                    }
                }
            }

            for (blasint i = 0; i < M; i += MC) {
                blasint current_M = stdmin(MC, M - i);
                blasint M_pad = (current_M + 31) & ~31;
                size_t panel_stride_A = 32 * (size_t)current_K;
                float32x4_t valpha = vdupq_n_f32(alpha);

                blasint M_main = current_M & ~31;
                for (blasint ii = 0; ii < M_main; ii += 32) {
                    float *__restrict A_out = &A_pack[(ii / 32) * panel_stride_A];
                    for (blasint kk = 0; kk < current_K; ++kk) {
                        const float *__restrict A_col = &A[(k + kk) * (size_t)lda + i + ii];
                        float *__restrict A_out_ptr = &A_out[kk * 32];
                        for (int v = 0; v < 8; ++v) {
                            vst1q_f32(&A_out_ptr[v * 4], vmulq_f32(vld1q_f32(&A_col[v * 4]), valpha));
                        }
                    }
                }
                if (M_main < current_M) {
                    float *__restrict A_out = &A_pack[(M_main / 32) * panel_stride_A];
                    for (blasint kk = 0; kk < current_K; ++kk) {
                        const float *__restrict A_col = &A[(k + kk) * (size_t)lda + i + M_main];
                        for (blasint br = 0; br < current_M - M_main; ++br) {
                            A_out[kk * 32 + br] = A_col[br] * alpha;
                        }
                        for (int br = current_M - M_main; br < 32; ++br) {
                            A_out[kk * 32 + br] = 0.0f;
                        }
                    }
                }

                for (blasint jj = 0; jj < N_pad; jj += 32) {
                    blasint current_N_block = stdmin(32, current_N - jj);
                    float *B_ptr = &B_pack[(jj / 32) * panel_stride_B];
                    for (blasint ii = 0; ii < M_pad; ii += 32) {
                        int current_M_block = stdmin(32, current_M - ii);
                        float *A_ptr = &A_pack[(ii / 32) * panel_stride_A];
                        float *C_ptr = &C[(i + ii) + (j + jj) * ldc];

                        if (current_M_block == 32 && current_N_block == 32) {
                            sgemm_sme_compute_32x32_tile(current_K, A_ptr, B_ptr, C_ptr, ldc, beta_mode, &beta);
                        }
                        else {
                            alignas(256) float C_buffer[1024];
                            if (beta_mode != 0) {
                                for (int idx = 0; idx < 1024; ++idx) {
                                    C_buffer[idx] = 0.0f;
                                }
                                for (blasint bc = 0; bc < current_N_block; ++bc) {
                                    for (blasint br = 0; br < current_M_block; ++br) {
                                        C_buffer[br + bc * 32] = C_ptr[br + bc * ldc];
                                    }
                                }
                            }
                            sgemm_sme_compute_32x32_tile(current_K, A_ptr, B_ptr, C_buffer, 32, beta_mode, &beta);
                            for (blasint bc = 0; bc < current_N_block; ++bc) {
                                for (blasint br = 0; br < current_M_block; ++br) {
                                    C_ptr[br + bc * ldc] = C_buffer[br + bc * 32];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

static void sgemm_sme_TN(int M, int N, int K, float alpha, const float *A, int lda, const float *B, int ldb, float beta, float *C, int ldc)
{
    if (alpha == 0.0f || K == 0) {
        if (beta != 1.0f) {
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < M; ++i) {
                    C[i + j * ldc] = (beta == 0.0f) ? 0.0f : C[i + j * ldc] * beta;
                }
            }
        }
        return;
    }

    alignas(256) /*thread_local*/ static float A_pack[MC * KC];
    alignas(256) /*thread_local*/ static float B_pack[KC * NC];

#pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < N; j += NC) {
        int current_N = stdmin(NC, N - j);
        int N_pad = (current_N + 31) & ~31;
        for (int k = 0; k < K; k += KC) {
            int current_K = stdmin(KC, K - k);
            size_t panel_stride_B = 32 * (size_t)current_K;
            int beta_mode = (k == 0) ? ((beta == 0.0f) ? 0 : (beta != 1.0f ? 1 : 2)) : 2;

            int N_main = current_N & ~31;
            for (int jj = 0; jj < N_main; jj += 32) {
                float *__restrict B_out = &B_pack[(jj / 32) * panel_stride_B];
                for (int kk = 0; kk < current_K; ++kk) {
                    const float *__restrict B_in = &B[(j + jj) * (size_t)ldb + k + kk];
                    for (int bc = 0; bc < 32; ++bc) {
                        B_out[kk * 32 + bc] = B_in[bc * ldb];
                    }
                }
            }
            if (N_main < current_N) {
                float *__restrict B_out = &B_pack[(N_main / 32) * panel_stride_B];
                for (int bc = 0; bc < 32; ++bc) {
                    if (N_main + bc < current_N) {
                        const float *__restrict B_col = &B[(j + N_main + bc) * (size_t)ldb + k];
                        for (int kk = 0; kk < current_K; ++kk) {
                            B_out[kk * 32 + bc] = B_col[kk];
                        }
                    }
                    else {
                        for (int kk = 0; kk < current_K; ++kk) {
                            B_out[kk * 32 + bc] = 0.0f;
                        }
                    }
                }
            }

            for (int i = 0; i < M; i += MC) {
                int current_M = stdmin(MC, M - i);
                int M_pad = (current_M + 31) & ~31;
                size_t panel_stride_A = 32 * (size_t)current_K;

                int M_main = current_M & ~31;
                for (int ii = 0; ii < M_main; ii += 32) {
                    float *__restrict A_out = &A_pack[(ii / 32) * panel_stride_A];
                    for (int kk = 0; kk < current_K; ++kk) {
                        const float *__restrict A_in = &A[(i + ii) * (size_t)lda + k + kk];
                        for (int br = 0; br < 32; ++br) {
                            A_out[kk * 32 + br] = A_in[br * lda] * alpha;
                        }
                    }
                }
                if (M_main < current_M) {
                    float *__restrict A_out = &A_pack[(M_main / 32) * panel_stride_A];
                    for (int br = 0; br < 32; ++br) {
                        if (M_main + br < current_M) {
                            const float *__restrict A_row = &A[(i + M_main + br) * (size_t)lda + k];
                            for (int kk = 0; kk < current_K; ++kk) {
                                A_out[kk * 32 + br] = A_row[kk] * alpha;
                            }
                        }
                        else {
                            for (int kk = 0; kk < current_K; ++kk) {
                                A_out[kk * 32 + br] = 0.0f;
                            }
                        }
                    }
                }

                for (int jj = 0; jj < N_pad; jj += 32) {
                    int current_N_block = stdmin(32, current_N - jj);
                    float *B_ptr = &B_pack[(jj / 32) * panel_stride_B];
                    for (int ii = 0; ii < M_pad; ii += 32) {
                        int current_M_block = stdmin(32, current_M - ii);
                        float *A_ptr = &A_pack[(ii / 32) * panel_stride_A];
                        float *C_ptr = &C[(i + ii) + (j + jj) * ldc];

                        if (current_M_block == 32 && current_N_block == 32) {
                            sgemm_sme_compute_32x32_tile(current_K, A_ptr, B_ptr, C_ptr, ldc, beta_mode, &beta);
                        }
                        else {
                            alignas(256) float C_buffer[1024];
                            if (beta_mode != 0) {
                                for (int idx = 0; idx < 1024; ++idx) {
                                    C_buffer[idx] = 0.0f;
                                }
                                for (int bc = 0; bc < current_N_block; ++bc) {
                                    for (int br = 0; br < current_M_block; ++br) {
                                        C_buffer[br + bc * 32] = C_ptr[br + bc * ldc];
                                    }
                                }
                            }
                            sgemm_sme_compute_32x32_tile(current_K, A_ptr, B_ptr, C_buffer, 32, beta_mode, &beta);
                            for (int bc = 0; bc < current_N_block; ++bc) {
                                for (int br = 0; br < current_M_block; ++br) {
                                    C_ptr[br + bc * ldc] = C_buffer[br + bc * 32];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

static void sgemm_sme_NT(int M, int N, int K, float alpha, const float *A, int lda, const float *B, int ldb, float beta, float *C, int ldc)
{
    if (alpha == 0.0f || K == 0) {
        if (beta != 1.0f) {
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < M; ++i) {
                    C[i + j * ldc] = (beta == 0.0f) ? 0.0f : C[i + j * ldc] * beta;
                }
            }
        }
        return;
    }

    alignas(256) /*thread_local*/ static float A_pack[MC * KC];
    alignas(256) /*thread_local*/ static float B_pack[KC * NC];

#pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < N; j += NC) {
        int current_N = stdmin(NC, N - j);
        int N_pad = (current_N + 31) & ~31;
        for (int k = 0; k < K; k += KC) {
            int current_K = stdmin(KC, K - k);
            size_t panel_stride_B = 32 * (size_t)current_K;
            int beta_mode = (k == 0) ? ((beta == 0.0f) ? 0 : (beta != 1.0f ? 1 : 2)) : 2;

            int N_main = current_N & ~31;
            for (int jj = 0; jj < N_main; jj += 32) {
                float *__restrict B_out = &B_pack[(jj / 32) * panel_stride_B];
                for (int kk = 0; kk < current_K; ++kk) {
                    const float *__restrict B_col = &B[(k + kk) * (size_t)ldb + j + jj];
                    float *__restrict B_out_ptr = &B_out[kk * 32];
                    for (int v = 0; v < 8; ++v) {
                        vst1q_f32(&B_out_ptr[v * 4], vld1q_f32(&B_col[v * 4]));
                    }
                }
            }
            if (N_main < current_N) {
                float *__restrict B_out = &B_pack[(N_main / 32) * panel_stride_B];
                for (int kk = 0; kk < current_K; ++kk) {
                    const float *__restrict B_row = &B[(k + kk) * (size_t)ldb + j + N_main];
                    for (int bc = 0; bc < current_N - N_main; ++bc) {
                        B_out[kk * 32 + bc] = B_row[bc];
                    }
                    for (int bc = current_N - N_main; bc < 32; ++bc) {
                        B_out[kk * 32 + bc] = 0.0f;
                    }
                }
            }

            for (int i = 0; i < M; i += MC) {
                int current_M = stdmin(MC, M - i);
                int M_pad = (current_M + 31) & ~31;
                size_t panel_stride_A = 32 * (size_t)current_K;
                float32x4_t valpha = vdupq_n_f32(alpha);

                int M_main = current_M & ~31;
                for (int ii = 0; ii < M_main; ii += 32) {
                    float *__restrict A_out = &A_pack[(ii / 32) * panel_stride_A];
                    for (int kk = 0; kk < current_K; ++kk) {
                        const float *__restrict A_col = &A[(k + kk) * (size_t)lda + i + ii];
                        float *__restrict A_out_ptr = &A_out[kk * 32];
                        for (int v = 0; v < 8; ++v) {
                            vst1q_f32(&A_out_ptr[v * 4], vmulq_f32(vld1q_f32(&A_col[v * 4]), valpha));
                        }
                    }
                }
                if (M_main < current_M) {
                    float *__restrict A_out = &A_pack[(M_main / 32) * panel_stride_A];
                    for (int kk = 0; kk < current_K; ++kk) {
                        const float *__restrict A_col = &A[(k + kk) * (size_t)lda + i + M_main];
                        for (int br = 0; br < current_M - M_main; ++br) {
                            A_out[kk * 32 + br] = A_col[br] * alpha;
                        }
                        for (int br = current_M - M_main; br < 32; ++br) {
                            A_out[kk * 32 + br] = 0.0f;
                        }
                    }
                }

                for (int jj = 0; jj < N_pad; jj += 32) {
                    int current_N_block = stdmin(32, current_N - jj);
                    float *B_ptr = &B_pack[(jj / 32) * panel_stride_B];
                    for (int ii = 0; ii < M_pad; ii += 32) {
                        int current_M_block = stdmin(32, current_M - ii);
                        float *A_ptr = &A_pack[(ii / 32) * panel_stride_A];
                        float *C_ptr = &C[(i + ii) + (j + jj) * ldc];

                        if (current_M_block == 32 && current_N_block == 32) {
                            sgemm_sme_compute_32x32_tile(current_K, A_ptr, B_ptr, C_ptr, ldc, beta_mode, &beta);
                        }
                        else {
                            alignas(256) float C_buffer[1024];
                            if (beta_mode != 0) {
                                for (int idx = 0; idx < 1024; ++idx) {
                                    C_buffer[idx] = 0.0f;
                                }
                                for (int bc = 0; bc < current_N_block; ++bc) {
                                    for (int br = 0; br < current_M_block; ++br) {
                                        C_buffer[br + bc * 32] = C_ptr[br + bc * ldc];
                                    }
                                }
                            }
                            sgemm_sme_compute_32x32_tile(current_K, A_ptr, B_ptr, C_buffer, 32, beta_mode, &beta);
                            for (int bc = 0; bc < current_N_block; ++bc) {
                                for (int br = 0; br < current_M_block; ++br) {
                                    C_ptr[br + bc * ldc] = C_buffer[br + bc * 32];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

static void sgemm_sme_TT(int M, int N, int K, float alpha, const float *A, int lda, const float *B, int ldb, float beta, float *C, int ldc)
{
    if (alpha == 0.0f || K == 0) {
        if (beta != 1.0f) {
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < M; ++i) {
                    C[i + j * ldc] = (beta == 0.0f) ? 0.0f : C[i + j * ldc] * beta;
                }
            }
        }
        return;
    }

    alignas(256) /*thread_local*/ static float A_pack[MC * KC];
    alignas(256) /*thread_local*/ static float B_pack[KC * NC];

#pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < N; j += NC) {
        int current_N = stdmin(NC, N - j);
        int N_pad = (current_N + 31) & ~31;
        for (int k = 0; k < K; k += KC) {
            int current_K = stdmin(KC, K - k);
            size_t panel_stride_B = 32 * (size_t)current_K;
            int beta_mode = (k == 0) ? ((beta == 0.0f) ? 0 : (beta != 1.0f ? 1 : 2)) : 2;

            int N_main = current_N & ~31;
            for (int jj = 0; jj < N_main; jj += 32) {
                float *__restrict B_out = &B_pack[(jj / 32) * panel_stride_B];
                for (int kk = 0; kk < current_K; ++kk) {
                    const float *__restrict B_col = &B[(k + kk) * (size_t)ldb + j + jj];
                    float *__restrict B_out_ptr = &B_out[kk * 32];
                    for (int v = 0; v < 8; ++v) {
                        vst1q_f32(&B_out_ptr[v * 4], vld1q_f32(&B_col[v * 4]));
                    }
                }
            }
            if (N_main < current_N) {
                float *__restrict B_out = &B_pack[(N_main / 32) * panel_stride_B];
                for (int kk = 0; kk < current_K; ++kk) {
                    const float *__restrict B_row = &B[(k + kk) * (size_t)ldb + j + N_main];
                    for (int bc = 0; bc < current_N - N_main; ++bc) {
                        B_out[kk * 32 + bc] = B_row[bc];
                    }
                    for (int bc = current_N - N_main; bc < 32; ++bc) {
                        B_out[kk * 32 + bc] = 0.0f;
                    }
                }
            }

            for (int i = 0; i < M; i += MC) {
                int current_M = stdmin(MC, M - i);
                int M_pad = (current_M + 31) & ~31;
                size_t panel_stride_A = 32 * (size_t)current_K;

                int M_main = current_M & ~31;
                for (int ii = 0; ii < M_main; ii += 32) {
                    float *__restrict A_out = &A_pack[(ii / 32) * panel_stride_A];
                    for (int kk = 0; kk < current_K; ++kk) {
                        const float *__restrict A_in = &A[(i + ii) * (size_t)lda + k + kk];
                        for (int br = 0; br < 32; ++br) {
                            A_out[kk * 32 + br] = A_in[br * lda] * alpha;
                        }
                    }
                }
                if (M_main < current_M) {
                    float *__restrict A_out = &A_pack[(M_main / 32) * panel_stride_A];
                    for (int br = 0; br < 32; ++br) {
                        if (M_main + br < current_M) {
                            const float *__restrict A_row = &A[(i + M_main + br) * (size_t)lda + k];
                            for (int kk = 0; kk < current_K; ++kk) {
                                A_out[kk * 32 + br] = A_row[kk] * alpha;
                            }
                        }
                        else {
                            for (int kk = 0; kk < current_K; ++kk) {
                                A_out[kk * 32 + br] = 0.0f;
                            }
                        }
                    }
                }

                for (int jj = 0; jj < N_pad; jj += 32) {
                    int current_N_block = stdmin(32, current_N - jj);
                    float *B_ptr = &B_pack[(jj / 32) * panel_stride_B];
                    for (int ii = 0; ii < M_pad; ii += 32) {
                        int current_M_block = stdmin(32, current_M - ii);
                        float *A_ptr = &A_pack[(ii / 32) * panel_stride_A];
                        float *C_ptr = &C[(i + ii) + (j + jj) * ldc];

                        if (current_M_block == 32 && current_N_block == 32) {
                            sgemm_sme_compute_32x32_tile(current_K, A_ptr, B_ptr, C_ptr, ldc, beta_mode, &beta);
                        }
                        else {
                            alignas(256) float C_buffer[1024];
                            if (beta_mode != 0) {
                                for (int idx = 0; idx < 1024; ++idx) {
                                    C_buffer[idx] = 0.0f;
                                }
                                for (int bc = 0; bc < current_N_block; ++bc) {
                                    for (int br = 0; br < current_M_block; ++br) {
                                        C_buffer[br + bc * 32] = C_ptr[br + bc * ldc];
                                    }
                                }
                            }
                            sgemm_sme_compute_32x32_tile(current_K, A_ptr, B_ptr, C_buffer, 32, beta_mode, &beta);
                            for (int bc = 0; bc < current_N_block; ++bc) {
                                for (int br = 0; br < current_M_block; ++br) {
                                    C_ptr[br + bc * ldc] = C_buffer[br + bc * 32];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void CNAME(char *transa,  char *transb, BLASLONG m,  BLASLONG n,  BLASLONG k,  float *alpha,  float *a,  BLASLONG lda,  float *b,  BLASLONG ldb,  float *beta, float *c,  BLASLONG ldc)
{
    bool trans_a = (*transa == 'T' || *transa == 't' || *transa == 'C' || *transa == 'c');
    bool trans_b = (*transb == 'T' || *transb == 't' || *transb == 'C' || *transb == 'c');
    if (!trans_a && !trans_b) {
        sgemm_sme_NN(m, n, k, *alpha, a, lda, b, ldb, *beta, c, ldc);
    }
    else if (trans_a && !trans_b) {
        sgemm_sme_TN(m, n, k, *alpha, a, lda, b, ldb, *beta, c, ldc);
    }
    else if (!trans_a && trans_b) {
        sgemm_sme_NT(m, n, k, *alpha, a, lda, b, ldb, *beta, c, ldc);
    }
    else {
        sgemm_sme_TT(m, n, k, *alpha, a, lda, b, ldb, *beta, c, ldc);
    }
}
