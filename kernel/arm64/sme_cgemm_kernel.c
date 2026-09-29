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
typedef float _Complex cfloat;
static cfloat CMUL(cfloat a, cfloat b,bool conja, bool conjb) {
float ra=creal(a);
float rb=creal(b);
float ia=conja ? -cimag(a) : cimag(a);
float ib=conjb ? -cimag(b) : cimag(b);
float r1=ra*rb;
float r2=ia*ib;
float r=r1-r2;
float i=(ra+ia)*(rb+ib)-r1-r2;
cfloat res={r,i};
return res;
}
#define KERNEL_ALPHA 0
#define USE_VECTORIZED_PACKING 1


static cfloat czero={0.,0.};
static cfloat cone={1.,0.};

#define MC 256
#define KC 2048
#define NC 1024

static inline void cgemm_sme_compute_16x16_tile(blasint current_K, const float *A_ptr, const float *B_ptr, cfloat *C_ptr, size_t ldc, int beta_mode, const cfloat *beta_ptr)
{
    size_t ldc_bytes = ldc * sizeof(cfloat);

    asm volatile("smstart\n\t"
		 "ptrue p0.s\n\t"
                 "zero {za}\n\t"
                 "cmp %w[beta_mode], #0\n\t"
                 "b.eq 19f\n\t"

                 // FIX 1: Correct SVE mnemonics for loading 32-bit floats (Real/Imag)
                 "cmp %w[beta_mode], #1\n\t"
                 "b.ne 10f\n\t"
                 "ld1rw z30.s, p0/z, [%[beta_ptr]]\n\t" // Load beta.real()
                 "add x15, %[beta_ptr], #4\n\t"         // 4-byte offset for float
                 "ld1rw z31.s, p0/z, [x15]\n\t"         // Load beta.imag()
                 "10:\n\t"

                 "mov w12, #0\n\t"
                 "mov x13, %[c]\n\t"
                 "100:\n\t"
                 "ld1w z0.s, p0/z, [x13]\n\t"
                 "add x14, x13, #64\n\t"
                 "ld1w z1.s, p0/z, [x14]\n\t"
                 "uzp1 z2.s, z0.s, z1.s\n\t" // z2 = C_re
                 "uzp2 z3.s, z0.s, z1.s\n\t" // z3 = C_im

                 "cmp %w[beta_mode], #1\n\t"
                 "b.ne 101f\n\t"

                 // FIX 2: Fully implemented Complex Beta Multiplication using movprfx
                 "movprfx z4, z2\n\t"
                 "fmul z4.s, p0/m, z4.s, z30.s\n\t" // z4 = Cre * Bre
                 "movprfx z5, z3\n\t"
                 "fmul z5.s, p0/m, z5.s, z31.s\n\t" // z5 = Cim * Bim
                 "movprfx z6, z4\n\t"
                 "fsub z6.s, p0/m, z6.s, z5.s\n\t" // z6 = Cre' (Cre*Bre - Cim*Bim)

                 "movprfx z8, z2\n\t"
                 "fmul z8.s, p0/m, z8.s, z31.s\n\t" // z8 = Cre * Bim
                 "movprfx z9, z3\n\t"
                 "fmul z9.s, p0/m, z9.s, z30.s\n\t" // z9 = Cim * Bre
                 "movprfx z7, z8\n\t"
                 "fadd z7.s, p0/m, z7.s, z9.s\n\t" // z7 = Cim' (Cre*Bim + Cim*Bre)

                 "mova za0v.s[w12, 0], p0/m, z6.s\n\t"
                 "mova za2v.s[w12, 0], p0/m, z7.s\n\t"
                 "b 102f\n\t"

                 "101:\n\t" // Fallback: beta == 1.0
                 "mova za0v.s[w12, 0], p0/m, z2.s\n\t"
                 "mova za2v.s[w12, 0], p0/m, z3.s\n\t"

                 "102:\n\t"
                 "add x13, x13, %[ldc_bytes]\n\t"
                 "add w12, w12, #1\n\t"
                 "cmp w12, #16\n\t"
                 "b.ne 100b\n\t"

                 "19:\n\t"
                 "mov w10, %w[k]\n\t"
                 "cbz w10, 3f\n\t"
                 "11:\n\t"
                 "ld1w z0.s, p0/z, [%[a], #0, mul vl]\n\t"
                 "ld1w z1.s, p0/z, [%[a], #1, mul vl]\n\t"
                 "ld1w z2.s, p0/z, [%[b], #0, mul vl]\n\t"
                 "ld1w z3.s, p0/z, [%[b], #1, mul vl]\n\t"
                 "fmopa za0.s, p0/m, p0/m, z0.s, z2.s\n\t"
                 "fmopa za1.s, p0/m, p0/m, z1.s, z3.s\n\t"
                 "fmopa za2.s, p0/m, p0/m, z0.s, z3.s\n\t"
                 "fmopa za3.s, p0/m, p0/m, z1.s, z2.s\n\t"
                 "add %[a], %[a], #128\n\t"
                 "add %[b], %[b], #128\n\t"
                 "subs w10, w10, #1\n\t"
                 "b.ne 11b\n\t"

                 "3:\n\t"
                 "mov w12, #0\n\t"
                 "mov x13, %[c]\n\t"
                 "200:\n\t"
                 "mova z0.s, p0/m, za0v.s[w12, 0]\n\t"
                 "mova z1.s, p0/m, za1v.s[w12, 0]\n\t"
                 "mova z2.s, p0/m, za2v.s[w12, 0]\n\t"
                 "mova z3.s, p0/m, za3v.s[w12, 0]\n\t"

                 // FIX 3: Non-destructive SVE arithmetic on store
                 "movprfx z4, z0\n\t"
                 "fsub z4.s, p0/m, z4.s, z1.s\n\t"
                 "movprfx z5, z2\n\t"
                 "fadd z5.s, p0/m, z5.s, z3.s\n\t"

                 "zip1 z6.s, z4.s, z5.s\n\t"
                 "zip2 z7.s, z4.s, z5.s\n\t"
                 "st1w z6.s, p0, [x13]\n\t"
                 "add x14, x13, #64\n\t"
                 "st1w z7.s, p0, [x14]\n\t"
                 "add x13, x13, %[ldc_bytes]\n\t"
                 "add w12, w12, #1\n\t"
                 "cmp w12, #16\n\t"
                 "b.ne 200b\n\t"
		 "smstop\n\t"
                 "msr fpsr, xzr\n\t"
    : [a] "+&r"(A_ptr), [b] "+&r"(B_ptr)
    : [k] "r"(current_K), [c] "r"(C_ptr), [ldc_bytes] "r"(ldc_bytes), [beta_mode] "r"(beta_mode), [beta_ptr] "r"(beta_ptr)
    // Updated Clobber list to cover x15 and z8-z9, z30-z31
    :  "x10", "x12", "x13", "x14", "x15", "cc", "memory",
       "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7",
       "v8", "v9", "v10", "v11", "v12", "v13", "v14", "v15",
       "v16", "v17", "v18", "v19", "v20", "v21", "v22", "v23",
       "v24", "v25", "v26", "v27", "v28", "v29", "v30", "v31",
       "z0", "z1", "z2", "z3", "z4", "z5", "z6", "z7",
       "z8", "z9", "z10", "z11", "z12", "z13", "z14", "z15",
       "z16", "z17", "z18", "z19", "z20", "z21", "z22", "z23",
       "z24", "z25", "z26", "z27", "z28", "z29", "z30", "z31", "za",
       "p0", "p1", "p2", "p3", "p4", "p5", "p6", "p7",
       "p8", "p9", "p10", "p11", "p12", "p13", "p14", "p15");


}

static void cgemm_sme_NN(int M, int N, int K, const cfloat alpha, const cfloat *A, int lda, const cfloat *b, int ldb, const cfloat beta, cfloat *C, int ldc, bool conja, bool conjb)
{

    if (alpha == czero || K == 0) {
        if (beta == czero) {
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < M; ++i) {
                    C[i + j * ldc] = czero;
                }
            }
        } else {
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < M; ++i) {
                    C[i + j * ldc] *= beta;
                }
            }
	}
        return;
    }



    alignas(256) /*thread_local*/ static float A_pack[MC * KC * 2];
    alignas(256) /*thread_local*/ static float B_pack[KC * NC * 2];
#pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < N; j += NC) {
        int current_N = stdmin(NC, N - j);
        int N_pad = (current_N + 15) & ~15;
        for (int k = 0; k < K; k += KC) {
            int current_K = stdmin(KC, K - k);
            size_t panel_stride_B = 32 * (size_t)current_K;
            int beta_mode = (k == 0) ? ((beta==czero) ? 0 : (beta != cone ? 1 : 2)) : 2;

            for (int jj = 0; jj < current_N; jj += 16) {
                float *__restrict B_out = &B_pack[(jj / 16) * panel_stride_B];
                for (int kk = 0; kk < current_K; ++kk) {
                    for (int bc = 0; bc < 16; ++bc) {
                        cfloat val = (jj + bc < current_N) ? b[(j + jj + bc) * ldb + k + kk] : czero;
                        B_out[kk * 32 + bc] = creal(val);
                        B_out[kk * 32 + 16 + bc] = conjb ? -cimag(val) : cimag(val);
                    }
                }
            }
            for (int i = 0; i < M; i += MC) {
                int current_M = stdmin(MC, M - i);
                int M_pad = (current_M + 15) & ~15;
                size_t panel_stride_A = 32 * (size_t)current_K;
                for (int ii = 0; ii < current_M; ii += 16) {
                    float *__restrict A_out = &A_pack[(ii / 16) * panel_stride_A];
                    for (int kk = 0; kk < current_K; ++kk) {
                        for (int br = 0; br < 16; ++br) {
                            cfloat val = (ii + br < current_M) ? CMUL(A[(k + kk) * lda + i + ii + br] , alpha, conja, 0) : czero;
                            A_out[kk * 32 + br] = creal(val);
                            A_out[kk * 32 + 16 + br] = cimag(val);
                        }
                    }
                }
                for (int jj = 0; jj < N_pad; jj += 16) {
                    int current_N_block = stdmin(16, current_N - jj);
                    float *B_ptr = &B_pack[(jj / 16) * panel_stride_B];
                    for (int ii = 0; ii < M_pad; ii += 16) {
                        int current_M_block = stdmin(16, current_M - ii);
                        float *A_ptr = &A_pack[(ii / 16) * panel_stride_A];
                        cfloat *C_ptr = &C[(i + ii) + (j + jj) * ldc];
                        if (current_M_block == 16 && current_N_block == 16) {
                            cgemm_sme_compute_16x16_tile(current_K, A_ptr, B_ptr, C_ptr, ldc, beta_mode, &beta);
                        }
                        else {
                            alignas(256) cfloat C_buffer[256];
                            if (beta_mode == 0) {
                                for (int idx = 0; idx < 256; ++idx) {
                                    C_buffer[idx] =  czero;
                                }
                            }
                                for (int bc = 0; bc < current_N_block; ++bc) {
                                    for (int br = 0; br < current_M_block; ++br) {
                                        C_buffer[br + bc * 16] = C_ptr[br + bc * ldc];
                                    }
                                }
                            cgemm_sme_compute_16x16_tile(current_K, A_ptr, B_ptr, C_buffer, 16, beta_mode, &beta);
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

static void cgemm_sme_TN(int M, int N, int K, const cfloat alpha, const cfloat *A, int lda, const cfloat *b, int ldb, const cfloat beta, cfloat *C, int ldc, bool conja, bool conjb)
{
    if (alpha == czero || K == 0) {
        if (beta == czero) {
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < M; ++i) {
                    C[i + j * ldc] = czero;
                }
            }
        } else {
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < M; ++i) {
                    C[i + j * ldc] = CMUL(C[i+j*ldc],beta,0,0);
                }
            }
	}
        return;
    }
    alignas(256) /*thread_local*/ static float A_pack[MC * KC * 2];
    alignas(256) /*thread_local*/ static float B_pack[KC * NC * 2];
#pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < N; j += NC) {
        int current_N = stdmin(NC, N - j);
        int N_pad = (current_N + 15) & ~15;
        for (int k = 0; k < K; k += KC) {
            int current_K = stdmin(KC, K - k);
            size_t panel_stride_B = 32 * (size_t)current_K;
            int beta_mode = (k == 0) ? ((beta==czero) ? 0 : (beta != cone ? 1 : 2)) : 2;

            for (int jj = 0; jj < current_N; jj += 16) {
                float *__restrict B_out = &B_pack[(jj / 16) * panel_stride_B];
                for (int kk = 0; kk < current_K; ++kk) {
                    for (int bc = 0; bc < 16; ++bc) {
                        cfloat val = (jj + bc < current_N) ? b[(j + jj + bc) * ldb + k + kk] : czero;
                        B_out[kk * 32 + bc] = creal(val);
                        B_out[kk * 32 + 16 + bc] = conjb ? -cimag(val) : cimag(val);
                    }
                }
            }
            for (int i = 0; i < M; i += MC) {
                int current_M = stdmin(MC, M - i);
                int M_pad = (current_M + 15) & ~15;
                size_t panel_stride_A = 32 * (size_t)current_K;
                for (int ii = 0; ii < current_M; ii += 16) {
                    float *__restrict A_out = &A_pack[(ii / 16) * panel_stride_A];
                    for (int kk = 0; kk < current_K; ++kk) {
                        for (int br = 0; br < 16; ++br) {
                            cfloat val = (ii + br < current_M) ? CMUL(A[(i + ii + br) * lda + k + kk],alpha,conja,0) : czero;
                            A_out[kk * 32 + br] = creal(val);
                            A_out[kk * 32 + 16 + br] = cimag(val);
                        }
                    }
                }
                for (int jj = 0; jj < N_pad; jj += 16) {
                    int current_N_block = stdmin(16, current_N - jj);
                    float *B_ptr = &B_pack[(jj / 16) * panel_stride_B];
                    for (int ii = 0; ii < M_pad; ii += 16) {
                        int current_M_block = stdmin(16, current_M - ii);
                        float *A_ptr = &A_pack[(ii / 16) * panel_stride_A];
                        cfloat *C_ptr = &C[(i + ii) + (j + jj) * ldc];
                        if (current_M_block == 16 && current_N_block == 16) {
                            cgemm_sme_compute_16x16_tile(current_K, A_ptr, B_ptr, C_ptr, ldc, beta_mode, &beta);
                        }
                        else {
                            alignas(256) cfloat C_buffer[256];
                                for (int idx = 0; idx < 256; ++idx) {
                                    C_buffer[idx] =  czero;
                                }
                            if (beta_mode != 0) {
                                for (int bc = 0; bc < current_N_block; ++bc) {
                                    for (int br = 0; br < current_M_block; ++br) {
                                        C_buffer[br + bc * 16] = C_ptr[br + bc * ldc];
                                    }
                                }
                            
                            }
                            cgemm_sme_compute_16x16_tile(current_K, A_ptr, B_ptr, C_buffer, 16, beta_mode, &beta);
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

static void cgemm_sme_NT(int M, int N, int K, const cfloat alpha, const cfloat *A, int lda, const cfloat *b, int ldb, const cfloat beta, cfloat *C, int ldc, bool conja, bool conjb)
{
    if (alpha == czero || K == 0) {
        if (beta == czero) {
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < M; ++i) {
                    C[i + j * ldc] = czero;
                }
            }
        } else {
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < M; ++i) {
                    C[i + j * ldc] *= beta;
                }
            }
	}
        return;
    }
    alignas(256) /*thread_local*/ static float A_pack[MC * KC * 2];
    alignas(256) /*thread_local*/ static float B_pack[KC * NC * 2];
#pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < N; j += NC) {
        int current_N = stdmin(NC, N - j);
        int N_pad = (current_N + 15) & ~15;
        for (int k = 0; k < K; k += KC) {
            int current_K = stdmin(KC, K - k);
            size_t panel_stride_B = 32 * (size_t)current_K;
            int beta_mode = (k == 0) ? ((beta==czero) ? 0 : (beta != cone ? 1 : 2)) : 2;

            for (int jj = 0; jj < current_N; jj += 16) {
                float *__restrict B_out = &B_pack[(jj / 16) * panel_stride_B];
                for (int kk = 0; kk < current_K; kk++) {
                    for (int bc = 0; bc < 16; ++bc) {
                        cfloat val = (jj + bc < current_N) ? b[(k + kk) * ldb + j + jj + bc] : czero ;
                        B_out[kk * 32 + bc] = creal(val);
                        B_out[kk * 32 + 16 + bc] = (conjb) ? -cimag(val) : cimag(val);
                    }
                }
            }
            for (int i = 0; i < M; i += MC) {
                int current_M = stdmin(MC, M - i);
                int M_pad = (current_M + 15) & ~15;
                size_t panel_stride_A = 32 * (size_t)current_K;
                for (int ii = 0; ii < current_M; ii += 16) {
                    float *__restrict A_out = &A_pack[(ii / 16) * panel_stride_A];
                    for (int kk = 0; kk < current_K; kk++) {
                        for (int br = 0; br < 16; ++br) {
                            cfloat val = (ii + br < current_M) ? CMUL(A[(k + kk) * lda + i + ii + br],alpha,conja,0) :  czero;
                            A_out[kk * 32 + br] = creal(val);
                            A_out[kk * 32 + 16 + br] = cimag(val);
                        }
                    }
                }
                for (int jj = 0; jj < N_pad; jj += 16) {
                    int current_N_block = stdmin(16, current_N - jj);
                    float *B_ptr = &B_pack[(jj / 16) * panel_stride_B];
                    for (int ii = 0; ii < M_pad; ii += 16) {
                        int current_M_block = stdmin(16, current_M - ii);
                        float *A_ptr = &A_pack[(ii / 16) * panel_stride_A];
                        cfloat *C_ptr = &C[(i + ii) + (j + jj) * ldc];
                        if (current_M_block == 16 && current_N_block == 16) {
                            cgemm_sme_compute_16x16_tile(current_K, A_ptr, B_ptr, C_ptr, ldc, beta_mode, &beta);
                        }
                        else {
                            alignas(256) cfloat C_buffer[256];
                            if (beta_mode != 0) {
                                for (int idx = 0; idx < 256; ++idx) {
                                    C_buffer[idx] =  czero;
                                }
                            }
                                for (int bc = 0; bc < current_N_block; ++bc) {
                                    for (int br = 0; br < current_M_block; ++br) {
                                        C_buffer[br + bc * 16] = C_ptr[br + bc * ldc];
                                    }
                                }
                        
                            cgemm_sme_compute_16x16_tile(current_K, A_ptr, B_ptr, C_buffer, 16, beta_mode, &beta);
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

static void cgemm_sme_TT(int M, int N, int K, const cfloat alpha, const cfloat *A, int lda, const cfloat *b, int ldb, const cfloat beta, cfloat *C, int ldc, bool conja, bool conjb)
{
    if (alpha == czero || K == 0) {
        if (beta == czero) {
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < M; ++i) {
                    C[i + j * ldc] = czero;
                }
            }
        } else {
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < M; ++i) {
                    C[i + j * ldc] *= beta;
                }
            }
	}
        return;
    }
    alignas(256) /*thread_local*/ static float A_pack[MC * KC * 2];
    alignas(256) /*thread_local*/ static float B_pack[KC * NC * 2];
#pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < N; j += NC) {
        int current_N = stdmin(NC, N - j);
        int N_pad = (current_N + 15) & ~15;
        for (int k = 0; k < K; k += KC) {
            int current_K = stdmin(KC, K - k);
            size_t panel_stride_B = 32 * (size_t)current_K;
            int beta_mode = (k == 0) ? ((beta == czero) ? 0 : (beta != cone ? 1 : 2)) : 2;

            for (int jj = 0; jj < current_N; jj += 16) {
                float *__restrict B_out = &B_pack[(jj / 16) * panel_stride_B];
                for (int kk = 0; kk < current_K; kk++) {
                    for (int bc = 0; bc < 16; ++bc) {
                        cfloat val = (jj + bc < current_N) ? b[(k + kk) * ldb + j + jj + bc] : czero;
                        B_out[kk * 32 + bc] = creal(val);
                        B_out[kk * 32 + 16 + bc] = (conjb) ? -cimag(val) : cimag(val);
                    }
                }
            }
            for (int i = 0; i < M; i += MC) {
                int current_M = stdmin(MC, M - i);
                int M_pad = (current_M + 15) & ~15;
                size_t panel_stride_A = 32 * (size_t)current_K;
                for (int ii = 0; ii < current_M; ii += 16) {
                    float *__restrict A_out = &A_pack[(ii / 16) * panel_stride_A];
                    for (int kk = 0; kk < current_K; kk++) {
                        for (int br = 0; br < 16; ++br) {
                            cfloat val = (ii + br < current_M) ? CMUL( A[(i + ii + br) * lda + k + kk],alpha,conja,0) : czero;
                            A_out[kk * 32 + br] = creal(val);
                            A_out[kk * 32 + 16 + br] = cimag(val);
                        }
                    }
                }
                for (int jj = 0; jj < N_pad; jj += 16) {
                    int current_N_block = stdmin(16, current_N - jj);
                    float *B_ptr = &B_pack[(jj / 16) * panel_stride_B];
                    for (int ii = 0; ii < M_pad; ii += 16) {
                        int current_M_block = stdmin(16, current_M - ii);
                        float *A_ptr = &A_pack[(ii / 16) * panel_stride_A];
                        cfloat *C_ptr = &C[(i + ii) + (j + jj) * ldc];
                        if (current_M_block == 16 && current_N_block == 16) {
                            cgemm_sme_compute_16x16_tile(current_K, A_ptr, B_ptr, C_ptr, ldc, beta_mode, &beta);
                        }
                        else {
                            alignas(256) cfloat C_buffer[256];
                            if (beta_mode != 0) {
                                for (int idx = 0; idx < 256; ++idx) {
                                    C_buffer[idx] =  czero;
                                }
			    }
                                for (int bc = 0; bc < current_N_block; ++bc) {
                                    for (int br = 0; br < current_M_block; ++br) {
                                        C_buffer[br + bc * 16] = C_ptr[br + bc * ldc];
                                    }
                                }
                            
                            cgemm_sme_compute_16x16_tile(current_K, A_ptr, B_ptr, C_buffer, 16, beta_mode, &beta);
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

void CNAME(const char *transa, const char *transb, const BLASLONG m, const BLASLONG n, const BLASLONG k, const float alpha_r, const float alpha_i, const float *a, const BLASLONG lda, const float *b, const BLASLONG ldb, const float beta_r, const float beta_i, float *c, const BLASLONG ldc)
{
cfloat alpha={alpha_r,alpha_i};
cfloat beta={beta_r,beta_i};
    bool trans_a = (*transa == 'T' || *transa == 't' || *transa == 'C' || *transa == 'c');
    bool trans_b = (*transb == 'T' || *transb == 't' || *transb == 'C' || *transb == 'c');
    if (!trans_a && !trans_b) {
    bool conja=(*transa == 'R' || *transa == 'r');
    bool conjb=(*transb == 'R' || *transb == 'r');
        cgemm_sme_NN(m, n, k, alpha, (const cfloat*) a, lda, (const cfloat*) b, ldb, beta, (cfloat*)c, ldc, conja, conjb);
    }
    else if (trans_a && !trans_b) {
    bool conja=(*transa == 'C' || *transa == 'c');
    bool conjb=(*transb == 'R' || *transb == 'r');
        cgemm_sme_TN(m, n, k, alpha, (const cfloat*) a, lda, (const cfloat*) b, ldb, beta, (cfloat*)c, ldc, conja, conjb);
    }
    else if (!trans_a && trans_b) {
    bool conja=(*transa == 'R' || *transa == 'r');
    bool conjb=(*transb == 'C' || *transb == 'c');
        cgemm_sme_NT(m, n, k, alpha, (const cfloat*) a, lda, (const cfloat*) b, ldb, beta, (cfloat*)c, ldc, conja, conjb);
    }
    else {
    bool conja=(*transa == 'C' || *transa == 'c');
    bool conjb=(*transb == 'C' || *transb == 'c');
        cgemm_sme_TT(m, n, k, alpha, (const cfloat*) a, lda, (const cfloat*) b, ldb, beta, (cfloat*)c, ldc, conja, conjb);
    }
}
