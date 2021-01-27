#include "sbgemm.h"

#include <immintrin.h>
// Walk around those intrinsics that missed by compiler
#define MM256_LOADU_EPI16(addr)   \
            _mm256_maskz_loadu_epi16(~0, (addr))
#define MM256_STOREU_EPI16(addr, reg)  \
            _mm256_mask_storeu_epi16((addr), ~0, (reg))

#include <stdio.h>
void print_block(BLASLONG m, BLASLONG n, bfloat16 * mat)
{
    printf("---- BLOCK %ld x %ld ----\n", m, n);
    for (BLASLONG i=0; i<m; i++) {
        for (BLASLONG j=0; j<n; j++) {
            printf("%-4X  ", *(mat + i*n +j));
        }
        printf("\n");
    }
    printf("---- End of BLOCK ----\n");
}

void COL_MAJOR_INCOPY_KERNEL_Kx32(BLASLONG k, bfloat16 * A, BLASLONG lda, bfloat16 * block_A)
{
    BLASLONG tag_k_2x = k & (~1);

    __m512i array512_0, array512_1, array512_2, array512_3;

    BLASLONG idx_src_base0, idx_src_base1;
    BLASLONG idx_target_base0, idx_target_base1;

    BLASLONG LDA_2x = 2*lda;
    BLASLONG BF16_BLOCK_T_M_2x = 2*32;
    idx_src_base0 = 0;
    idx_src_base1 = lda;
    idx_target_base0 = 0;
    idx_target_base1 = 32;
    for (BLASLONG idx_k = 0; idx_k < tag_k_2x; idx_k += 2) {
        array512_0 = _mm512_loadu_si512(&A[idx_src_base0]);
        array512_1 = _mm512_loadu_si512(&A[idx_src_base1]);
        array512_2 = _mm512_unpacklo_epi16(array512_0, array512_1);
        array512_3 = _mm512_unpackhi_epi16(array512_0, array512_1);
        _mm512_storeu_si512(&block_A[idx_target_base0], array512_2);
        _mm512_storeu_si512(&block_A[idx_target_base1], array512_3);

        idx_src_base0 += LDA_2x;
        idx_src_base1 += LDA_2x;
        idx_target_base0 += BF16_BLOCK_T_M_2x;
        idx_target_base1 += BF16_BLOCK_T_M_2x;
    }

    if (tag_k_2x != k) {
        __m512i ZERO512 = _mm512_setzero_si512();
        array512_0 = _mm512_loadu_si512(&A[idx_src_base0]);
        array512_2 = _mm512_unpacklo_epi16(array512_0, ZERO512);
        array512_3 = _mm512_unpackhi_epi16(array512_0, ZERO512);
        _mm512_storeu_si512(&block_A[idx_target_base0], array512_2);
        _mm512_storeu_si512(&block_A[idx_target_base1], array512_3);
   }

#ifdef DEBUG_PROFILE
   print_block(BF16_BLOCK_THRES_K, BF16_BLOCK_THRES_M, block_A);
#endif
}

void COL_MAJOR_INCOPY_KERNEL_Kx32m(BLASLONG k, BLASLONG m, bfloat16 * A, BLASLONG lda, bfloat16 * block_A)
{
    BLASLONG tag_k_2x = k & (~1);
    unsigned int tail_mask_value = (((unsigned int)0xffffffff) >> (32-m));
    __mmask32 tail_mask = *((__mmask32*) &tail_mask_value);

    __m512i array512_0, array512_1, array512_2, array512_3;

    BLASLONG idx_src_base0, idx_src_base1;
    BLASLONG idx_target_base0, idx_target_base1;

    BLASLONG LDA_2x = 2*lda;
    BLASLONG BF16_BLOCK_T_M_2x = 2*32;
    idx_src_base0 = 0;
    idx_src_base1 = lda;
    idx_target_base0 = 0;
    idx_target_base1 = 32;
    for (BLASLONG idx_k = 0; idx_k < tag_k_2x; idx_k += 2) {
        array512_0 = _mm512_maskz_loadu_epi16(tail_mask, &A[idx_src_base0]);
        array512_1 = _mm512_maskz_loadu_epi16(tail_mask, &A[idx_src_base1]);
        array512_2 = _mm512_unpacklo_epi16(array512_0, array512_1);
        array512_3 = _mm512_unpackhi_epi16(array512_0, array512_1);
        _mm512_storeu_si512(&block_A[idx_target_base0], array512_2);
        _mm512_storeu_si512(&block_A[idx_target_base1], array512_3);

        idx_src_base0 += LDA_2x;
        idx_src_base1 += LDA_2x;
        idx_target_base0 += BF16_BLOCK_T_M_2x;
        idx_target_base1 += BF16_BLOCK_T_M_2x;
    }

    if (tag_k_2x != k) {
        __m512i ZERO512 = _mm512_setzero_si512();
        array512_0 = _mm512_maskz_loadu_epi16(tail_mask, &A[idx_src_base0]);
        array512_2 = _mm512_unpacklo_epi16(array512_0, ZERO512);
        array512_3 = _mm512_unpackhi_epi16(array512_0, ZERO512);
        _mm512_storeu_si512(&block_A[idx_target_base0], array512_2);
        _mm512_storeu_si512(&block_A[idx_target_base1], array512_3);
    }

#ifdef DEBUG_PROFILE
   print_block(BF16_BLOCK_THRES_K, BF16_BLOCK_THRES_M, block_A);
#endif
}

void COL_MAJOR_INCOPY_KERNEL_Kx16(BLASLONG k, BLASLONG m, bfloat16 * A, BLASLONG lda, bfloat16 * block_A)
{
    BLASLONG tag_k_2x = k & (~1);

    __m256i array256_0, array256_1, array256_2, array256_3;

    BLASLONG idx_src_base0, idx_src_base1;
    BLASLONG idx_target_base0;

    BLASLONG LDA_2x = 2*lda;
    idx_src_base0 = 0;
    idx_src_base1 = lda;
    idx_target_base0 = 0;
    for (BLASLONG idx_k = 0; idx_k < tag_k_2x; idx_k += 2) {
        array256_0 = MM256_LOADU_EPI16(&A[idx_src_base0]);
        array256_1 = MM256_LOADU_EPI16(&A[idx_src_base1]);
        array256_2 = _mm256_unpacklo_epi16(array256_0, array256_1);
        array256_3 = _mm256_unpackhi_epi16(array256_0, array256_1);
        // Store in one row of block_B
        MM256_STOREU_EPI16(&block_A[idx_target_base0],      array256_2);
        MM256_STOREU_EPI16(&block_A[idx_target_base0 + 16], array256_3);

        idx_src_base0 += LDA_2x;
        idx_src_base1 += LDA_2x;
        idx_target_base0 += 32;
    }

    if (tag_k_2x != k) {
        __m256i ZERO256 = _mm256_setzero_si256();
        array256_0 = MM256_LOADU_EPI16(&A[idx_src_base0]);
        array256_2 = _mm256_unpacklo_epi16(array256_0, ZERO256);
        array256_3 = _mm256_unpackhi_epi16(array256_0, ZERO256);
        // Store in one row of block_B
        MM256_STOREU_EPI16(&block_A[idx_target_base0],      array256_2);
        MM256_STOREU_EPI16(&block_A[idx_target_base0 + 16], array256_3);
    }

#ifdef DEBUG_PROFILE
   print_block(BF16_BLOCK_THRES_K, BF16_BLOCK_THRES_M, block_A);
#endif
}

void COL_MAJOR_INCOPY_KERNEL_Kx16m(BLASLONG k, BLASLONG m, bfloat16 * A, BLASLONG lda, bfloat16 * block_A)
{
    BLASLONG tag_k_2x = k & (~1);
    unsigned short tail_mask_value = (((unsigned short)0xffff) >> (16-m));
    __mmask16 tail_mask = *((__mmask16*) &tail_mask_value);

    __m256i array256_0, array256_1, array256_2, array256_3;

    BLASLONG idx_src_base0, idx_src_base1;
    BLASLONG idx_target_base0;

    BLASLONG LDA_2x = 2*lda;
    idx_src_base0 = 0;
    idx_src_base1 = lda;
    idx_target_base0 = 0;
    for (BLASLONG idx_k = 0; idx_k < tag_k_2x; idx_k += 2) {
        array256_0 = _mm256_maskz_loadu_epi16(tail_mask, &A[idx_src_base0]);
        array256_1 = _mm256_maskz_loadu_epi16(tail_mask, &A[idx_src_base1]);
        array256_2 = _mm256_unpacklo_epi16(array256_0, array256_1);
        array256_3 = _mm256_unpackhi_epi16(array256_0, array256_1);
        // Store in one row of block_B
        MM256_STOREU_EPI16(&block_A[idx_target_base0],      array256_2);
        MM256_STOREU_EPI16(&block_A[idx_target_base0 + 16], array256_3);

        idx_src_base0 += LDA_2x;
        idx_src_base1 += LDA_2x;
        idx_target_base0 += 32;
    }

    if (tag_k_2x != k) {
        __m256i ZERO256 = _mm256_setzero_si256();
        array256_0 = _mm256_maskz_loadu_epi16(tail_mask, &A[idx_src_base0]);
        array256_2 = _mm256_unpacklo_epi16(array256_0, ZERO256);
        array256_3 = _mm256_unpackhi_epi16(array256_0, ZERO256);
        // Store in one row of block_B
        MM256_STOREU_EPI16(&block_A[idx_target_base0],      array256_2);
        MM256_STOREU_EPI16(&block_A[idx_target_base0 + 16], array256_3);
    }

#ifdef DEBUG_PROFILE
   print_block(BF16_BLOCK_THRES_K, BF16_BLOCK_THRES_M, block_A);
#endif
}

void COL_MAJOR_ONCOPY_KERNEL_8x32(BLASLONG k, bfloat16 * B, BLASLONG ldb, bfloat16 * block_B)
{
    BLASLONG tag_k_32x = k & (~31);
    BLASLONG idx_src_base0, idx_src_base1, idx_src_base2, idx_src_base3, idx_src_base4, idx_src_base5, idx_src_base6, idx_src_base7;
    BLASLONG idx_target_base0;

    idx_src_base0 = 0;
    idx_src_base1 = 1*ldb;
    idx_src_base2 = 2*ldb;
    idx_src_base3 = 3*ldb;
    idx_src_base4 = 4*ldb;
    idx_src_base5 = 5*ldb;
    idx_src_base6 = 6*ldb;
    idx_src_base7 = 7*ldb;
    idx_target_base0 = 0;

    for (BLASLONG idx_k = 0; idx_k < tag_k_32x; idx_k += 32) {
        _mm512_storeu_si512(&block_B[idx_target_base0+ 32*0], _mm512_loadu_si512(&B[idx_src_base0+idx_k]));
        _mm512_storeu_si512(&block_B[idx_target_base0+ 32*1], _mm512_loadu_si512(&B[idx_src_base1+idx_k]));
        _mm512_storeu_si512(&block_B[idx_target_base0+ 32*2], _mm512_loadu_si512(&B[idx_src_base2+idx_k]));
        _mm512_storeu_si512(&block_B[idx_target_base0+ 32*3], _mm512_loadu_si512(&B[idx_src_base3+idx_k]));
        _mm512_storeu_si512(&block_B[idx_target_base0+ 32*4], _mm512_loadu_si512(&B[idx_src_base4+idx_k]));
        _mm512_storeu_si512(&block_B[idx_target_base0+ 32*5], _mm512_loadu_si512(&B[idx_src_base5+idx_k]));
        _mm512_storeu_si512(&block_B[idx_target_base0+ 32*6], _mm512_loadu_si512(&B[idx_src_base6+idx_k]));
        _mm512_storeu_si512(&block_B[idx_target_base0+ 32*7], _mm512_loadu_si512(&B[idx_src_base7+idx_k]));
        idx_target_base0 += 32*8;
    }

    if (tag_k_32x != k) {
        unsigned int tail_mask_value = (((unsigned int)0xffffffff) >> (32-(k-tag_k_32x)));
        __mmask32 tail_mask = *((__mmask32*) &tail_mask_value);
        _mm512_storeu_si512(&block_B[idx_target_base0+ 32*0], _mm512_maskz_loadu_epi16(tail_mask, &B[idx_src_base0+tag_k_32x]));
        _mm512_storeu_si512(&block_B[idx_target_base0+ 32*1], _mm512_maskz_loadu_epi16(tail_mask, &B[idx_src_base1+tag_k_32x]));
        _mm512_storeu_si512(&block_B[idx_target_base0+ 32*2], _mm512_maskz_loadu_epi16(tail_mask, &B[idx_src_base2+tag_k_32x]));
        _mm512_storeu_si512(&block_B[idx_target_base0+ 32*3], _mm512_maskz_loadu_epi16(tail_mask, &B[idx_src_base3+tag_k_32x]));
        _mm512_storeu_si512(&block_B[idx_target_base0+ 32*4], _mm512_maskz_loadu_epi16(tail_mask, &B[idx_src_base4+tag_k_32x]));
        _mm512_storeu_si512(&block_B[idx_target_base0+ 32*5], _mm512_maskz_loadu_epi16(tail_mask, &B[idx_src_base5+tag_k_32x]));
        _mm512_storeu_si512(&block_B[idx_target_base0+ 32*6], _mm512_maskz_loadu_epi16(tail_mask, &B[idx_src_base6+tag_k_32x]));
        _mm512_storeu_si512(&block_B[idx_target_base0+ 32*7], _mm512_maskz_loadu_epi16(tail_mask, &B[idx_src_base7+tag_k_32x]));
    }

#ifdef DEBUG_PROFILE
   print_block(BF16_BLOCK_THRES_N, BF16_BLOCK_THRES_K, block_B);
#endif
}

void COL_MAJOR_ONCOPY_KERNEL_Nx32(BLASLONG n, BLASLONG k, bfloat16 * B, BLASLONG ldb, bfloat16 * block_B)
{
    BLASLONG tag_k_32x = k & (~31);
    BLASLONG tag_n_2x  = n & (~1);
    BLASLONG idx_src_base0;
    BLASLONG idx_target_base0;

    BLASLONG LDB_2x = 2*ldb;

    idx_target_base0 = 0;

    for (BLASLONG idx_k = 0; idx_k < tag_k_32x; idx_k += 32) {
        idx_src_base0 = 0;
        for (BLASLONG idx_n = 0; idx_n < tag_n_2x; idx_n += 2) {
            _mm512_storeu_si512(&block_B[idx_target_base0+ 32*0], _mm512_loadu_si512(&B[idx_src_base0       + idx_k]));
            _mm512_storeu_si512(&block_B[idx_target_base0+ 32*1], _mm512_loadu_si512(&B[idx_src_base0 + ldb + idx_k]));
            idx_src_base0 += LDB_2x;
            idx_target_base0 += 64;
        }
        
        if (tag_n_2x != n) {
            _mm512_storeu_si512(&block_B[idx_target_base0],  _mm512_loadu_si512(&B[idx_src_base0       + idx_k]));
            idx_target_base0 += 32;
        }
    }

    if (tag_k_32x != k) {
        unsigned int tail_mask_value = (((unsigned int)0xffffffff) >> (32-(k-tag_k_32x)));
        __mmask32 tail_mask = *((__mmask32*) &tail_mask_value);
        idx_src_base0 = 0;
        for (BLASLONG idx_n = 0; idx_n < tag_n_2x; idx_n += 2) {
            _mm512_storeu_si512(&block_B[idx_target_base0+ 32*0], _mm512_maskz_loadu_epi16(tail_mask, &B[idx_src_base0       + tag_k_32x]));
            _mm512_storeu_si512(&block_B[idx_target_base0+ 32*1], _mm512_maskz_loadu_epi16(tail_mask, &B[idx_src_base0 + ldb + tag_k_32x]));
            idx_src_base0 += LDB_2x;
            idx_target_base0 += 64;
        }
        
        if (tag_n_2x != n) {
            _mm512_storeu_si512(&block_B[idx_target_base0],  _mm512_maskz_loadu_epi16(tail_mask, &B[idx_src_base0       + tag_k_32x]));
        }
    }

#ifdef DEBUG_PROFILE
   print_block(BF16_BLOCK_THRES_N, BF16_BLOCK_THRES_K, block_B);
#endif
}

// Scale matrix C while beta is not ZERO or ONE
void sbgemm_scal_operation(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N, OPENBLAS_CONST float beta, float *C, OPENBLAS_CONST blasint ldc)
{
    BLASLONG tag_n_Nx = N & (~3);
    BLASLONG tag_n_Mx = M & (~15);

    BLASLONG LDC4x = ldc*4;
    BLASLONG idx_base_0 = 0;
    BLASLONG idx_base_1 = ldc;
    BLASLONG idx_base_2 = ldc*2;
    BLASLONG idx_base_3 = ldc*3;

    unsigned short tail_mask_value = (((unsigned short)0xffff) >> (16-M+tag_n_Mx));
    __mmask16 tail_mask = *((__mmask16*) &tail_mask_value);

    __m512 array_512_0, array_512_1, array_512_2, array_512_3;

    __m512  BETAVECTOR  = _mm512_set1_ps(beta);

    if (Order == CblasColMajor) {
        for (BLASLONG idx_n = 0; idx_n < tag_n_Nx; idx_n += 4) {
            for (BLASLONG idx_m = 0; idx_m < tag_n_Mx; idx_m += 16) {
                array_512_0 = _mm512_loadu_ps(&C[idx_base_0+idx_m]);
                array_512_1 = _mm512_loadu_ps(&C[idx_base_1+idx_m]);
                array_512_2 = _mm512_loadu_ps(&C[idx_base_2+idx_m]);
                array_512_3 = _mm512_loadu_ps(&C[idx_base_3+idx_m]);

                array_512_0 = _mm512_mul_ps(BETAVECTOR, array_512_0);
                array_512_1 = _mm512_mul_ps(BETAVECTOR, array_512_1);
                array_512_2 = _mm512_mul_ps(BETAVECTOR, array_512_2);
                array_512_3 = _mm512_mul_ps(BETAVECTOR, array_512_3);

                _mm512_storeu_ps(&C[idx_base_0+idx_m], array_512_0);
                _mm512_storeu_ps(&C[idx_base_1+idx_m], array_512_1);
                _mm512_storeu_ps(&C[idx_base_2+idx_m], array_512_2);
                _mm512_storeu_ps(&C[idx_base_3+idx_m], array_512_3);
            }

            if (tag_n_Mx != M) {
                array_512_0 = _mm512_maskz_loadu_ps(tail_mask, &C[idx_base_0+tag_n_Mx]);
                array_512_1 = _mm512_maskz_loadu_ps(tail_mask, &C[idx_base_1+tag_n_Mx]);
                array_512_2 = _mm512_maskz_loadu_ps(tail_mask, &C[idx_base_2+tag_n_Mx]);
                array_512_3 = _mm512_maskz_loadu_ps(tail_mask, &C[idx_base_3+tag_n_Mx]);

                array_512_0 = _mm512_mul_ps(BETAVECTOR, array_512_0);
                array_512_1 = _mm512_mul_ps(BETAVECTOR, array_512_1);
                array_512_2 = _mm512_mul_ps(BETAVECTOR, array_512_2);
                array_512_3 = _mm512_mul_ps(BETAVECTOR, array_512_3);

                _mm512_mask_storeu_ps(&C[idx_base_0+tag_n_Mx], tail_mask, array_512_0);
                _mm512_mask_storeu_ps(&C[idx_base_1+tag_n_Mx], tail_mask, array_512_1);
                _mm512_mask_storeu_ps(&C[idx_base_2+tag_n_Mx], tail_mask, array_512_2);
                _mm512_mask_storeu_ps(&C[idx_base_3+tag_n_Mx], tail_mask, array_512_3);
            }

            idx_base_0 += LDC4x;
            idx_base_1 += LDC4x;
            idx_base_2 += LDC4x;
            idx_base_3 += LDC4x;
        }

        if (tag_n_Nx != N) {
            for (BLASLONG idx_n = tag_n_Nx; idx_n < N; idx_n++) {
                for (BLASLONG idx_m = 0; idx_m < tag_n_Mx; idx_m += 16) {
                    array_512_0 = _mm512_loadu_ps(&C[idx_base_0+idx_m]);
                    array_512_0 = _mm512_mul_ps(BETAVECTOR, array_512_0);
                    _mm512_storeu_ps(&C[idx_base_0+idx_m], array_512_0);
                }

                if (tag_n_Mx != M) {
                    array_512_0 = _mm512_maskz_loadu_ps(tail_mask, &C[idx_base_0+tag_n_Mx]);
                    array_512_0 = _mm512_mul_ps(BETAVECTOR, array_512_0);
                    _mm512_mask_storeu_ps(&C[idx_base_0+tag_n_Mx], tail_mask, array_512_0);
                }
                idx_base_0 += ldc;
            }
        }
    } else {

    }
}

// Scale matrix C while beta is not ZERO or ONE
void sbgemm_zero_operation(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N, float *C, OPENBLAS_CONST blasint ldc)
{
    BLASLONG tag_n_Nx = N & (~3);
    BLASLONG tag_n_Mx = M & (~15);

    BLASLONG LDC4x = ldc*4;
    BLASLONG idx_base_0 = 0;
    BLASLONG idx_base_1 = ldc;
    BLASLONG idx_base_2 = ldc*2;
    BLASLONG idx_base_3 = ldc*3;

    unsigned short tail_mask_value = (((unsigned short)0xffff) >> (16-M+tag_n_Mx));
    __mmask16 tail_mask = *((__mmask16*) &tail_mask_value);

    __m512  ZEROVECTOR  = _mm512_setzero_ps();

    if (Order == CblasColMajor) {
        for (BLASLONG idx_n = 0; idx_n < tag_n_Nx; idx_n += 4) {
            for (BLASLONG idx_m = 0; idx_m < tag_n_Mx; idx_m += 16) {
                _mm512_storeu_ps(&C[idx_base_0+idx_m], ZEROVECTOR);
                _mm512_storeu_ps(&C[idx_base_1+idx_m], ZEROVECTOR);
                _mm512_storeu_ps(&C[idx_base_2+idx_m], ZEROVECTOR);
                _mm512_storeu_ps(&C[idx_base_3+idx_m], ZEROVECTOR);
            }

            if (tag_n_Mx != M) {
                _mm512_mask_storeu_ps(&C[idx_base_0+tag_n_Mx], tail_mask, ZEROVECTOR);
                _mm512_mask_storeu_ps(&C[idx_base_1+tag_n_Mx], tail_mask, ZEROVECTOR);
                _mm512_mask_storeu_ps(&C[idx_base_2+tag_n_Mx], tail_mask, ZEROVECTOR);
                _mm512_mask_storeu_ps(&C[idx_base_3+tag_n_Mx], tail_mask, ZEROVECTOR);
            }

            idx_base_0 += LDC4x;
            idx_base_1 += LDC4x;
            idx_base_2 += LDC4x;
            idx_base_3 += LDC4x;
        }

        if (tag_n_Nx != N) {
            for (BLASLONG idx_n = tag_n_Nx; idx_n < N; idx_n++) {
                for (BLASLONG idx_m = 0; idx_m < tag_n_Mx; idx_m += 16) {
                    _mm512_storeu_ps(&C[idx_base_0+idx_m], ZEROVECTOR);
                }

                if (tag_n_Mx != M) {
                    _mm512_mask_storeu_ps(&C[idx_base_0+tag_n_Mx], tail_mask, ZEROVECTOR);
                }
                idx_base_0 += ldc;
            }
        }
    } else {

    }
}