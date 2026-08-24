#include <arm_neon.h>
#include <stddef.h>
#include <stdalign.h>
#include <stdbool.h>
#include "common.h"
#ifndef stdmin
#define stdmin(a,b)   (a>b? b:a)
#endif

typedef double _Complex zdouble;

static zdouble CDMUL(zdouble a, zdouble b,bool conja, bool conjb) {
double ra=creal(a);
double rb=creal(b);
double ia=conja ? -cimag(a) : cimag(a);
double ib=conjb ? -cimag(b) : cimag(b);
double r1=ra*rb;
double r2=ia*ib;
double r=r1-r2;
double i=(ra+ia)*(rb+ib)-r1-r2;
zdouble res={r,i};
return res;
}


static zdouble cdzero={0.,0.};
static zdouble cdone={1.,0.};


#define MC 128
#define KC 1024
#define NC 512

static inline void zgemm_sme_compute_8x8_tile(int current_K, const double *A_ptr, const double *B_ptr, zdouble *C_ptr, size_t ldc, int beta_mode, const zdouble *beta_ptr)
{
    size_t ldc_bytes = ldc * sizeof(zdouble);

    asm volatile("smstart\n\t"
		 "ptrue p0.d\n\t"
                 "zero {za}\n\t"
                 "cmp %w[beta_mode], #0\n\t"
                 "b.eq 19f\n\t"

                 // FIX: Load Beta Real and Imaginary components and replicate them
                 "cmp %w[beta_mode], #1\n\t"
                 "b.ne 10f\n\t"
                 "ld1rd z30.d, p0/z, [%[beta_ptr]]\n\t" // Load beta.real()
                 "add x15, %[beta_ptr], #8\n\t"
                 "ld1rd z31.d, p0/z, [x15]\n\t" // Load beta.imag()
                 "10:\n\t"

                 "mov w12, #0\n\t"
                 "mov x13, %[c]\n\t"
                 "100:\n\t"
                 "ld1d z0.d, p0/z, [x13]\n\t"
                 "add x14, x13, #64\n\t"
                 "ld1d z1.d, p0/z, [x14]\n\t"
                 "uzp1 z2.d, z0.d, z1.d\n\t" // z2 = C_re
                 "uzp2 z3.d, z0.d, z1.d\n\t" // z3 = C_im

                 "cmp %w[beta_mode], #1\n\t"
                 "b.ne 101f\n\t"

                 // FIX: Fully implemented Complex Beta Multiplication using movprfx
                 "movprfx z4, z2\n\t"
                 "fmul z4.d, p0/m, z4.d, z30.d\n\t" // z4 = Cre * Bre
                 "movprfx z5, z3\n\t"
                 "fmul z5.d, p0/m, z5.d, z31.d\n\t" // z5 = Cim * Bim
                 "movprfx z6, z4\n\t"
                 "fsub z6.d, p0/m, z6.d, z5.d\n\t" // z6 = Cre' (Cre*Bre - Cim*Bim)

                 "movprfx z8, z2\n\t"
                 "fmul z8.d, p0/m, z8.d, z31.d\n\t" // z8 = Cre * Bim
                 "movprfx z9, z3\n\t"
                 "fmul z9.d, p0/m, z9.d, z30.d\n\t" // z9 = Cim * Bre
                 "movprfx z7, z8\n\t"
                 "fadd z7.d, p0/m, z7.d, z9.d\n\t" // z7 = Cim' (Cre*Bim + Cim*Bre)

                 "mova za0v.d[w12, 0], p0/m, z6.d\n\t"
                 "mova za2v.d[w12, 0], p0/m, z7.d\n\t"
                 "b 102f\n\t"

                 "101:\n\t" // Fallback: beta == 1.0
                 "mova za0v.d[w12, 0], p0/m, z2.d\n\t"
                 "mova za2v.d[w12, 0], p0/m, z3.d\n\t"

                 "102:\n\t"
                 "add x13, x13, %[ldc_bytes]\n\t"
                 "add w12, w12, #1\n\t"
                 "cmp w12, #8\n\t"
                 "b.ne 100b\n\t"

                 "19:\n\t"
                 "mov w10, %w[k]\n\t"
                 "cbz w10, 3f\n\t"
                 "11:\n\t"
                 "ld1d z0.d, p0/z, [%[a], #0, mul vl]\n\t"
                 "ld1d z1.d, p0/z, [%[a], #1, mul vl]\n\t"
                 "ld1d z2.d, p0/z, [%[b], #0, mul vl]\n\t"
                 "ld1d z3.d, p0/z, [%[b], #1, mul vl]\n\t"
                 "fmopa za0.d, p0/m, p0/m, z0.d, z2.d\n\t"
                 "fmopa za1.d, p0/m, p0/m, z1.d, z3.d\n\t"
                 "fmopa za2.d, p0/m, p0/m, z0.d, z3.d\n\t"
                 "fmopa za3.d, p0/m, p0/m, z1.d, z2.d\n\t"
                 "add %[a], %[a], #128\n\t"
                 "add %[b], %[b], #128\n\t"
                 "subs w10, w10, #1\n\t"
                 "b.ne 11b\n\t"

                 "3:\n\t"
                 "mov w12, #0\n\t"
                 "mov x13, %[c]\n\t"
                 "200:\n\t"
                 "mova z0.d, p0/m, za0v.d[w12, 0]\n\t"
                 "mova z1.d, p0/m, za1v.d[w12, 0]\n\t"
                 "mova z2.d, p0/m, za2v.d[w12, 0]\n\t"
                 "mova z3.d, p0/m, za3v.d[w12, 0]\n\t"

                 // FIX: Non-destructive SVE arithmetic
                 "movprfx z4, z0\n\t"
                 "fsub z4.d, p0/m, z4.d, z1.d\n\t"
                 "movprfx z5, z2\n\t"
                 "fadd z5.d, p0/m, z5.d, z3.d\n\t"

                 "zip1 z6.d, z4.d, z5.d\n\t"
                 "zip2 z7.d, z4.d, z5.d\n\t"
                 "st1d z6.d, p0, [x13]\n\t"
                 "add x14, x13, #64\n\t"
                 "st1d z7.d, p0, [x14]\n\t"
                 "add x13, x13, %[ldc_bytes]\n\t"
                 "add w12, w12, #1\n\t"
                 "cmp w12, #8\n\t"
                 "b.ne 200b\n\t"
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
      "z24", "z25", "z26", "z27", "z28", "z29", "z30", "z31", "za",
      "p0", "p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8",
      "p9", "p10", "p11", "p12", "p13", "p14", "p15");

}

static void zgemm_sme_NN(int M, int N, int K, const zdouble alpha, const zdouble *A, int lda, const zdouble *b, int ldb, const zdouble beta, zdouble *C, int ldc, bool conja, bool conjb)
{

    if (alpha == cdzero || K == 0) {
        if (beta != cdone) {
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < M; ++i) {
                    C[i + j * ldc] = (beta == cdzero) ? 0.0f : C[i + j * ldc] * beta;
                }
            }
        }
        return;
    }

    alignas(256) /*thread_local*/ static double A_pack[MC * KC * 2];
    alignas(256) /*thread_local*/ static double B_pack[KC * NC * 2];
#pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < N; j += NC) {
        int current_N = stdmin(NC, N - j);
        int N_pad = (current_N + 7) & ~7;
        for (int k = 0; k < K; k += KC) {
            int current_K = stdmin(KC, K - k);
            size_t panel_stride_B = 16 * (size_t)current_K;
            int beta_mode = (k == 0) ? ((beta == cdzero) ? 0 : (beta != cdone ? 1 : 2)) : 2;

            for (int jj = 0; jj < current_N; jj += 8) {
                double *__restrict B_out = &B_pack[(jj / 8) * panel_stride_B];
                for (int kk = 0; kk < current_K; ++kk) {
                    for (int bc = 0; bc < 8; ++bc) {
                        zdouble val = (jj + bc < current_N) ? b[(j + jj + bc) * ldb + k + kk] : cdzero;
                        B_out[kk * 16 + bc] = creal(val);
                        B_out[kk * 16 + 8 + bc] = conjb ? -cimag(val) : cimag(val);
                    }
                }
            }
            for (int i = 0; i < M; i += MC) {
                int current_M = stdmin(MC, M - i);
                int M_pad = (current_M + 7) & ~7;
                size_t panel_stride_A = 16 * (size_t)current_K;
                for (int ii = 0; ii < current_M; ii += 8) {
                    double *__restrict A_out = &A_pack[(ii / 8) * panel_stride_A];
                    for (int kk = 0; kk < current_K; ++kk) {
                        for (int br = 0; br < 8; ++br) {
                            zdouble val = (ii + br < current_M) ? CDMUL(A[(k + kk) * lda + i + ii + br] ,alpha,conja,0) : cdzero;
                            A_out[kk * 16 + br] = creal(val);
                            A_out[kk * 16 + 8 + br] = cimag(val);
                        }
                    }
                }
                for (int jj = 0; jj < N_pad; jj += 8) {
                    int current_N_block = stdmin(8, current_N - jj);
                    double *B_ptr = &B_pack[(jj / 8) * panel_stride_B];
                    for (int ii = 0; ii < M_pad; ii += 8) {
                        int current_M_block = stdmin(8, current_M - ii);
                        double *A_ptr = &A_pack[(ii / 8) * panel_stride_A];
                        zdouble *C_ptr = &C[(i + ii) + (j + jj) * ldc];
                        if (current_M_block == 8 && current_N_block == 8) {
                            zgemm_sme_compute_8x8_tile(current_K, A_ptr, B_ptr, C_ptr, ldc, beta_mode, &beta);
                        }
                        else {
                            alignas(256) zdouble C_buffer[64];
                            
                                for (int idx = 0; idx < 64; ++idx) {
                                    C_buffer[idx] = cdzero;
                                }
                            if (beta_mode != 0) {
                                for (int bc = 0; bc < current_N_block; ++bc) {
                                    for (int br = 0; br < current_M_block; ++br) {
                                        C_buffer[br + bc * 8] = C_ptr[br + bc * ldc];
                                    }
                                }
                            }
                            zgemm_sme_compute_8x8_tile(current_K, A_ptr, B_ptr, C_buffer, 8, beta_mode, &beta);
                            for (int bc = 0; bc < current_N_block; ++bc) {
                                for (int br = 0; br < current_M_block; ++br) {
                                    C_ptr[br + bc * ldc] = C_buffer[br + bc * 8];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

static void zgemm_sme_TN(int M, int N, int K, const zdouble alpha, const zdouble *A, int lda, const zdouble *b, int ldb, const zdouble beta, zdouble *C, int ldc, bool conja, bool conjb)
{
    if (alpha == cdzero || K == 0) {
        if (beta != cdone) {
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < M; ++i) {
                    C[i + j * ldc] = (beta == cdzero) ? 0.0f : C[i + j * ldc] * beta;
                }
            }
        }
        return;
    }

    alignas(256) /*thread_local*/ static double A_pack[MC * KC * 2];
    alignas(256) /*thread_local*/ static double B_pack[KC * NC * 2];
#pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < N; j += NC) {
        int current_N = stdmin(NC, N - j);
        int N_pad = (current_N + 7) & ~7;
        for (int k = 0; k < K; k += KC) {
            int current_K = stdmin(KC, K - k);
            size_t panel_stride_B = 16 * (size_t)current_K;
            int beta_mode = (k == 0) ? ((beta == cdzero) ? 0 : (beta != cdone ? 1 : 2)) : 2;

            for (int jj = 0; jj < current_N; jj += 8) {
                double *__restrict B_out = &B_pack[(jj / 8) * panel_stride_B];
                for (int kk = 0; kk < current_K; ++kk) {
                    for (int bc = 0; bc < 8; ++bc) {
                        zdouble val = (jj + bc < current_N) ? b[(j + jj + bc) * ldb + k + kk] : cdzero;
                        B_out[kk * 16 + bc] = creal(val);
                        B_out[kk * 16 + 8 + bc] = conjb ? -cimag(val) : cimag(val);
                    }
                }
            }
            for (int i = 0; i < M; i += MC) {
                int current_M = stdmin(MC, M - i);
                int M_pad = (current_M + 7) & ~7;
                size_t panel_stride_A = 16 * (size_t)current_K;
                for (int ii = 0; ii < current_M; ii += 8) {
                    double *__restrict A_out = &A_pack[(ii / 8) * panel_stride_A];
                    for (int kk = 0; kk < current_K; ++kk) {
                        for (int br = 0; br < 8; ++br) {
                            zdouble val = (ii + br < current_M) ? CDMUL(A[(i + ii + br) * lda + k + kk],alpha,conja,0) : cdzero;
                            A_out[kk * 16 + br] = creal(val);
                            A_out[kk * 16 + 8 + br] = cimag(val);
                        }
                    }
                }
                for (int jj = 0; jj < N_pad; jj += 8) {
                    int current_N_block = stdmin(8, current_N - jj);
                    double *B_ptr = &B_pack[(jj / 8) * panel_stride_B];
                    for (int ii = 0; ii < M_pad; ii += 8) {
                        int current_M_block = stdmin(8, current_M - ii);
                        double *A_ptr = &A_pack[(ii / 8) * panel_stride_A];
                        zdouble *C_ptr = &C[(i + ii) + (j + jj) * ldc];
                        if (current_M_block == 8 && current_N_block == 8) {
                            zgemm_sme_compute_8x8_tile(current_K, A_ptr, B_ptr, C_ptr, ldc, beta_mode, &beta);
                        }
                        else {
                            alignas(256) zdouble C_buffer[64];
                            
                                for (int idx = 0; idx < 64; ++idx) {
                                    C_buffer[idx] = cdzero;
                                }
                            if (beta_mode != 0) {
                                for (int bc = 0; bc < current_N_block; ++bc) {
                                    for (int br = 0; br < current_M_block; ++br) {
                                        C_buffer[br + bc * 8] = C_ptr[br + bc * ldc];
                                    }
                                }
                            }
                            zgemm_sme_compute_8x8_tile(current_K, A_ptr, B_ptr, C_buffer, 8, beta_mode, &beta);
                            for (int bc = 0; bc < current_N_block; ++bc) {
                                for (int br = 0; br < current_M_block; ++br) {
                                    C_ptr[br + bc * ldc] = C_buffer[br + bc * 8];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

static void zgemm_sme_NT(blasint M, blasint N, blasint K, const zdouble alpha, const zdouble *A, blasint lda, const zdouble *b, blasint ldb, const zdouble beta, zdouble *C, blasint ldc, bool conja, bool conjb)
{
#if 0
    if (alpha == cdzero || K == 0) {
        if (beta != cdone) {
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < M; ++i) {
                    C[i + j * ldc] = (beta == cdzero) ? 0.0f : C[i + j * ldc] * beta;
                }
            }
        }
        return;
    }
#endif
    if (alpha == cdzero || K == 0) {
        if (beta == cdzero) {
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < M; ++i) {
                    C[i + j * ldc] = cdzero;
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


    alignas(256) /*thread_local*/ static double A_pack[MC * KC * 2];
    alignas(256) /*thread_local*/ static double B_pack[KC * NC * 2 *2];
#pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < N; j += NC) {
        int current_N = stdmin(NC, N - j);
        int N_pad = (current_N + 7) & ~7;
        for (int k = 0; k < K; k += KC) {
            int current_K = stdmin(KC, K - k);
            size_t panel_stride_B = 16 * (size_t)current_K;
            int beta_mode = (k == 0) ? ((beta == cdzero) ? 0 : (beta != cdone ? 1 : 2)) : 2;

            for (int jj = 0; jj < current_N; jj += 8) {
                double *__restrict B_out = &B_pack[(jj / 8) * panel_stride_B];
                for (int kk = 0; kk < current_K; ++kk) {
                    for (int bc = 0; bc < 8; ++bc) {
                        zdouble val = (jj + bc < current_N) ? b[(k + kk) * ldb + j + jj + bc] : cdzero;
                        B_out[kk * 16 + bc] = creal(val);
                        B_out[kk * 16 + 8 + bc] = conjb ? -cimag(val) : cimag(val);
                    }
                }
            }
            for (int i = 0; i < M; i += MC) {
                int current_M = stdmin(MC, M - i);
                int M_pad = (current_M + 7) & ~7;
                size_t panel_stride_A = 16 * (size_t)current_K;
                for (int ii = 0; ii < current_M; ii += 8) {
                    double *__restrict A_out = &A_pack[(ii / 8) * panel_stride_A];
                    for (int kk = 0; kk < current_K; ++kk) {
                        for (int br = 0; br < 8; ++br) {
                            zdouble val = (ii + br < current_M) ? CDMUL(A[(k + kk) * lda + i + ii + br] ,alpha,conja,0) : cdzero;
                            A_out[kk * 16 + br] = creal(val);
                            A_out[kk * 16 + 8 + br] = cimag(val);
                        }
                    }
                }
                for (int jj = 0; jj < N_pad; jj += 8) {
                    int current_N_block = stdmin(8, current_N - jj);
                    double *B_ptr = &B_pack[(jj / 8) * panel_stride_B];
                    for (int ii = 0; ii < M_pad; ii += 8) {
                        int current_M_block = stdmin(8, current_M - ii);
                        double *A_ptr = &A_pack[(ii / 8) * panel_stride_A];
                        zdouble *C_ptr = &C[(i + ii) + (j + jj) * ldc];
                        if (current_M_block == 8 && current_N_block == 8) {
                            zgemm_sme_compute_8x8_tile(current_K, A_ptr, B_ptr, C_ptr, ldc, beta_mode, &beta);
                        }
                        else {
                            alignas(256) zdouble C_buffer[64];
                            
                                for (int idx = 0; idx < 64; ++idx) {
                                    C_buffer[idx] = cdzero;
                                }
                            if (beta_mode != 0) {
                                for (int bc = 0; bc < current_N_block; ++bc) {
                                    for (int br = 0; br < current_M_block; ++br) {
                                        C_buffer[br + bc * 8] = C_ptr[br + bc * ldc];
                                    }
                                }
                            }
                            zgemm_sme_compute_8x8_tile(current_K, A_ptr, B_ptr, C_buffer, 8, beta_mode, &beta);
                            for (int bc = 0; bc < current_N_block; ++bc) {
                                for (int br = 0; br < current_M_block; ++br) {
                                    C_ptr[br + bc * ldc] = C_buffer[br + bc * 8];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

static void zgemm_sme_TT(int M, int N, int K, const zdouble alpha, const zdouble *A, int lda, const zdouble *b, int ldb, const zdouble beta, zdouble *C, int ldc, bool conja, bool conjb)
{
    if (alpha == cdzero || K == 0) {
        if (beta != cdone) {
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < M; ++i) {
                    C[i + j * ldc] = (beta == cdzero) ? 0.0f : C[i + j * ldc] * beta;
                }
            }
        }
        return;
    }

    alignas(256) /*thread_local*/ static double A_pack[MC * KC * 2];
    alignas(256) /*thread_local*/ static double B_pack[KC * NC * 2];
#pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < N; j += NC) {
        int current_N = stdmin(NC, N - j);
        int N_pad = (current_N + 7) & ~7;
        for (int k = 0; k < K; k += KC) {
            int current_K = stdmin(KC, K - k);
            size_t panel_stride_B = 16 * (size_t)current_K;
            int beta_mode = (k == 0) ? ((beta == cdzero) ? 0 : (beta != cdone ? 1 : 2)) : 2;

            for (int jj = 0; jj < current_N; jj += 8) {
                double *__restrict B_out = &B_pack[(jj / 8) * panel_stride_B];
                for (int kk = 0; kk < current_K; ++kk) {
                    for (int bc = 0; bc < 8; ++bc) {
                        zdouble val = (jj + bc < current_N) ? b[(k + kk) * ldb + j + jj + bc] : cdzero;
                        B_out[kk * 16 + bc] = creal(val);
                        B_out[kk * 16 + 8 + bc] = conjb ? -cimag(val) : cimag(val);
                    }
                }
            }
            for (int i = 0; i < M; i += MC) {
                int current_M = stdmin(MC, M - i);
                int M_pad = (current_M + 7) & ~7;
                size_t panel_stride_A = 16 * (size_t)current_K;
                for (int ii = 0; ii < current_M; ii += 8) {
                    double *__restrict A_out = &A_pack[(ii / 8) * panel_stride_A];
                    for (int kk = 0; kk < current_K; ++kk) {
                        for (int br = 0; br < 8; ++br) {
                            zdouble val = (ii + br < current_M) ? CDMUL(A[(i + ii + br) * lda + k + kk],alpha,conja,0) : cdzero;
                            A_out[kk * 16 + br] = creal(val);
                            A_out[kk * 16 + 8 + br] = cimag(val);
                        }
                    }
                }
                for (int jj = 0; jj < N_pad; jj += 8) {
                    int current_N_block = stdmin(8, current_N - jj);
                    double *B_ptr = &B_pack[(jj / 8) * panel_stride_B];
                    for (int ii = 0; ii < M_pad; ii += 8) {
                        int current_M_block = stdmin(8, current_M - ii);
                        double *A_ptr = &A_pack[(ii / 8) * panel_stride_A];
                        zdouble *C_ptr = &C[(i + ii) + (j + jj) * ldc];
                        if (current_M_block == 8 && current_N_block == 8) {
                            zgemm_sme_compute_8x8_tile(current_K, A_ptr, B_ptr, C_ptr, ldc, beta_mode, &beta);
                        }
                        else {
                            alignas(256) zdouble C_buffer[64];
                            
                                for (int idx = 0; idx < 64; ++idx) {
                                    C_buffer[idx] = cdzero;
                                }
                            if (beta_mode != 0) {
                                for (int bc = 0; bc < current_N_block; ++bc) {
                                    for (int br = 0; br < current_M_block; ++br) {
                                        C_buffer[br + bc * 8] = C_ptr[br + bc * ldc];
                                    }
                                }
                            }
                            zgemm_sme_compute_8x8_tile(current_K, A_ptr, B_ptr, C_buffer, 8, beta_mode, &beta);
                            for (int bc = 0; bc < current_N_block; ++bc) {
                                for (int br = 0; br < current_M_block; ++br) {
                                    C_ptr[br + bc * ldc] = C_buffer[br + bc * 8];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void CNAME(const char *transa, const char *transb, const BLASLONG m, const BLASLONG n, const BLASLONG k, const double alpha_r, const double alpha_i, const double *a, const BLASLONG lda, const double *b, const BLASLONG ldb, const double beta_r, const double beta_i, double *c, const BLASLONG ldc)
{
zdouble alpha={alpha_r,alpha_i};
zdouble beta={beta_r,beta_i};

    bool trans_a = (*transa == 'T' || *transa == 't' || *transa == 'C' || *transa == 'c');
    bool trans_b = (*transb == 'T' || *transb == 't' || *transb == 'C' || *transb == 'c');
    if (!trans_a && !trans_b) {
	bool conja = (*transa == 'R' || *transa == 'r');
	bool conjb = (*transb == 'R' || *transb == 'r');
        zgemm_sme_NN(m, n, k, alpha, (const zdouble*)a, lda, (const zdouble*)b, ldb, beta, (zdouble*)c, ldc, conja, conjb);
    }
    else if (trans_a && !trans_b) {
	bool conja = (*transa == 'C' || *transa == 'c');
	bool conjb = (*transb == 'R' || *transb == 'r');
        zgemm_sme_TN(m, n, k, alpha, (const zdouble*)a, lda, (const zdouble*)b, ldb, beta, (zdouble*)c, ldc, conja, conjb);
    }
    else if (!trans_a && trans_b) {
	bool conja = (*transa == 'R' || *transa == 'r');
	bool conjb = (*transb == 'C' || *transb == 'c');
        zgemm_sme_NT(m, n, k, alpha, (const zdouble*)a, lda, (const zdouble*)b, ldb, beta, (zdouble*)c, ldc, conja, conjb);
    }
    else {
	bool conja = (*transa == 'C' || *transa == 'c');
	bool conjb = (*transb == 'C' || *transb == 'c');
        zgemm_sme_TT(m, n, k, alpha, (const zdouble*)a, lda, (const zdouble*)b, ldb, beta, (zdouble*)c, ldc, conja, conjb);
    }
}

