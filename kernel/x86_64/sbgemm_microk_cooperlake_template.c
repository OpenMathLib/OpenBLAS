#include "sbgemm.h"
#include "bf16_common_macros.h"
#include <immintrin.h>

#undef STORE16_COMPLETE_RESULT
#undef STORE16_MASK_COMPLETE_RESULT
#undef SBGEMM_BLOCK_KERNEL_32x8x32
#undef SBGEMM_BLOCK_KERNEL_16x8x32
#undef SBGEMM_BLOCK_KERNEL_32xNx32
#undef SBGEMM_BLOCK_KERNEL_16xNx32
#undef SBGEMM_BLOCKING_KERNEL_2

#ifndef ONE_ALPHA      // ALPHA is not ONE
    #define STORE16_COMPLETE_RESULT       STORE16_COMPLETE_RESULT_ALPHA_ONE
    #define STORE16_MASK_COMPLETE_RESULT  STORE16_MASK_COMPLETE_RESULT_ALPHA_ONE
    #define SBGEMM_BLOCK_KERNEL_32x8x32   sbgemm_block_kernel_32x8x32_alpha
    #define SBGEMM_BLOCK_KERNEL_16x8x32   sbgemm_block_kernel_16x8x32_alpha
    #define SBGEMM_BLOCK_KERNEL_32xNx32   sbgemm_block_kernel_32xNx32_alpha
    #define SBGEMM_BLOCK_KERNEL_16xNx32   sbgemm_block_kernel_16xNx32_alpha
    #define SBGEMM_BLOCKING_KERNEL_2      sbgemm_blocking_kernel_2_alpha
#else                  // ALPHA is ONE
    #define STORE16_COMPLETE_RESULT       STORE16_COMPLETE_RESULT_ONE_ONE
    #define STORE16_MASK_COMPLETE_RESULT  STORE16_MASK_COMPLETE_RESULT_ONE_ONE
    #define SBGEMM_BLOCK_KERNEL_32x8x32   sbgemm_block_kernel_32x8x32_one
    #define SBGEMM_BLOCK_KERNEL_16x8x32   sbgemm_block_kernel_16x8x32_one
    #define SBGEMM_BLOCK_KERNEL_32xNx32   sbgemm_block_kernel_32xNx32_one
    #define SBGEMM_BLOCK_KERNEL_16xNx32   sbgemm_block_kernel_16xNx32_one
    #define SBGEMM_BLOCKING_KERNEL_2      sbgemm_blocking_kernel_2_one
#endif


// SBGEMM Kernel for 16<M<=32, N=8, K can be any number, but the processing will take 32 as a base
#ifndef ONE_ALPHA      // ALPHA is not ONE
void sbgemm_block_kernel_32x8x32_alpha(BLASLONG m, BLASLONG k, float alpha, bfloat16 *A, bfloat16 *B, float *C, int ldc)
#else                  // ALPHA is ONE
void sbgemm_block_kernel_32x8x32_one(BLASLONG m, BLASLONG k, float alpha, bfloat16 *A, bfloat16 *B, float *C, int ldc)
#endif
{
    int  SHUFFLE_MAGIC_NO = 0x39;
    BLASLONG tag_k_32x = k & (~31);
    BLASLONG idxA_base = 0;
    BLASLONG idxB_base = 0;
    BLASLONG width = 32;

#ifndef ONE_ALPHA
    __m512  ALPHAVECTOR = _mm512_set1_ps(alpha);
#endif

    __m512i arrayA_512_0, arrayA_512_1;
    __m512i arrayB_512_0, arrayB_512_1, arrayB_512_2, arrayB_512_3, arrayB_512_4, arrayB_512_5, arrayB_512_6, arrayB_512_7;
    __m512  result_512_0, result_512_1, result_512_2, result_512_3, result_512_4, result_512_5, result_512_6, result_512_7,
            result_512_8, result_512_9, result_512_10, result_512_11, result_512_12, result_512_13, result_512_14, result_512_15;
    __m512  result_512_tmp_0, result_512_tmp_1, result_512_tmp_2, result_512_tmp_3;

    __m512i M512_EPI32_8      = _mm512_set1_epi32(8);
    __m512i shuffle_idx_base0 = _mm512_set_epi32(23, 22, 21, 20, 7, 6, 5, 4, 19, 18, 17, 16, 3, 2, 1, 0);
    __m512i shuffle_idx_base1 = _mm512_add_epi32(shuffle_idx_base0, M512_EPI32_8);

    result_512_0  = _mm512_setzero_ps();
    result_512_1  = _mm512_setzero_ps();
    result_512_2  = _mm512_setzero_ps();
    result_512_3  = _mm512_setzero_ps();
    result_512_4  = _mm512_setzero_ps();
    result_512_5  = _mm512_setzero_ps();
    result_512_6  = _mm512_setzero_ps();
    result_512_7  = _mm512_setzero_ps();
    result_512_8  = _mm512_setzero_ps();
    result_512_9  = _mm512_setzero_ps();
    result_512_10 = _mm512_setzero_ps();
    result_512_11 = _mm512_setzero_ps();
    result_512_12 = _mm512_setzero_ps();
    result_512_13 = _mm512_setzero_ps();
    result_512_14 = _mm512_setzero_ps();
    result_512_15 = _mm512_setzero_ps();

    for (BLASLONG idx_k = 0; idx_k < k; idx_k += 32) {
        // Load B with unroll 8
        idxB_base = idx_k << 3;
        arrayB_512_0 = _mm512_loadu_si512(&B[idxB_base + 32*0]);
        arrayB_512_1 = _mm512_loadu_si512(&B[idxB_base + 32*1]);
        arrayB_512_2 = _mm512_loadu_si512(&B[idxB_base + 32*2]);
        arrayB_512_3 = _mm512_loadu_si512(&B[idxB_base + 32*3]);
        arrayB_512_4 = _mm512_loadu_si512(&B[idxB_base + 32*4]);
        arrayB_512_5 = _mm512_loadu_si512(&B[idxB_base + 32*5]);
        arrayB_512_6 = _mm512_loadu_si512(&B[idxB_base + 32*6]);
        arrayB_512_7 = _mm512_loadu_si512(&B[idxB_base + 32*7]);

        if (idx_k == tag_k_32x) {width = k - tag_k_32x;}

        for (BLASLONG idx = 0; idx < width;) {
            // Each two rows are a group for 32-pair bf16 elements
            idxA_base = idx << 5;
            arrayA_512_0 = _mm512_loadu_si512(&A[idxA_base]);
            arrayA_512_1 = _mm512_loadu_si512(&A[idxA_base + 32]);

            result_512_0  = _mm512_dpbf16_ps(result_512_0,  (__m512bh) arrayA_512_0, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_0)));
            result_512_1  = _mm512_dpbf16_ps(result_512_1,  (__m512bh) arrayA_512_0, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_1)));
            result_512_2  = _mm512_dpbf16_ps(result_512_2,  (__m512bh) arrayA_512_0, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_2)));
            result_512_3  = _mm512_dpbf16_ps(result_512_3,  (__m512bh) arrayA_512_0, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_3)));
            result_512_4  = _mm512_dpbf16_ps(result_512_4,  (__m512bh) arrayA_512_0, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_4)));
            result_512_5  = _mm512_dpbf16_ps(result_512_5,  (__m512bh) arrayA_512_0, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_5)));
            result_512_6  = _mm512_dpbf16_ps(result_512_6,  (__m512bh) arrayA_512_0, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_6)));
            result_512_7  = _mm512_dpbf16_ps(result_512_7,  (__m512bh) arrayA_512_0, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_7)));
            result_512_8  = _mm512_dpbf16_ps(result_512_8,  (__m512bh) arrayA_512_1, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_0)));
            result_512_9  = _mm512_dpbf16_ps(result_512_9,  (__m512bh) arrayA_512_1, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_1)));
            result_512_10 = _mm512_dpbf16_ps(result_512_10, (__m512bh) arrayA_512_1, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_2)));
            result_512_11 = _mm512_dpbf16_ps(result_512_11, (__m512bh) arrayA_512_1, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_3)));
            result_512_12 = _mm512_dpbf16_ps(result_512_12, (__m512bh) arrayA_512_1, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_4)));
            result_512_13 = _mm512_dpbf16_ps(result_512_13, (__m512bh) arrayA_512_1, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_5)));
            result_512_14 = _mm512_dpbf16_ps(result_512_14, (__m512bh) arrayA_512_1, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_6)));
            result_512_15 = _mm512_dpbf16_ps(result_512_15, (__m512bh) arrayA_512_1, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_7)));

            arrayB_512_0 = _mm512_shuffle_epi32(arrayB_512_0, SHUFFLE_MAGIC_NO);
            arrayB_512_1 = _mm512_shuffle_epi32(arrayB_512_1, SHUFFLE_MAGIC_NO);
            arrayB_512_2 = _mm512_shuffle_epi32(arrayB_512_2, SHUFFLE_MAGIC_NO);
            arrayB_512_3 = _mm512_shuffle_epi32(arrayB_512_3, SHUFFLE_MAGIC_NO);
            arrayB_512_4 = _mm512_shuffle_epi32(arrayB_512_4, SHUFFLE_MAGIC_NO);
            arrayB_512_5 = _mm512_shuffle_epi32(arrayB_512_5, SHUFFLE_MAGIC_NO);
            arrayB_512_6 = _mm512_shuffle_epi32(arrayB_512_6, SHUFFLE_MAGIC_NO);
            arrayB_512_7 = _mm512_shuffle_epi32(arrayB_512_7, SHUFFLE_MAGIC_NO);

            idx += 2;
            // Every 4 loops we need to switch to next 128 bits of arrayB registers
            if ((idx & (~7)) == idx) {
                arrayB_512_0 = _mm512_shuffle_i32x4(arrayB_512_0, arrayB_512_0, SHUFFLE_MAGIC_NO);
                arrayB_512_1 = _mm512_shuffle_i32x4(arrayB_512_1, arrayB_512_1, SHUFFLE_MAGIC_NO);
                arrayB_512_2 = _mm512_shuffle_i32x4(arrayB_512_2, arrayB_512_2, SHUFFLE_MAGIC_NO);
                arrayB_512_3 = _mm512_shuffle_i32x4(arrayB_512_3, arrayB_512_3, SHUFFLE_MAGIC_NO);
                arrayB_512_4 = _mm512_shuffle_i32x4(arrayB_512_4, arrayB_512_4, SHUFFLE_MAGIC_NO);
                arrayB_512_5 = _mm512_shuffle_i32x4(arrayB_512_5, arrayB_512_5, SHUFFLE_MAGIC_NO);
                arrayB_512_6 = _mm512_shuffle_i32x4(arrayB_512_6, arrayB_512_6, SHUFFLE_MAGIC_NO);
                arrayB_512_7 = _mm512_shuffle_i32x4(arrayB_512_7, arrayB_512_7, SHUFFLE_MAGIC_NO);
            }
        }
    }

    if (m != 32) {
        unsigned short tail_mask_value = (((unsigned short)0xffff) >> (32-m));
        __mmask16 tail_mask = *((__mmask16*) &tail_mask_value);
        result_512_tmp_0 = _mm512_permutex2var_ps(result_512_0, shuffle_idx_base0, result_512_8);
        result_512_tmp_1 = _mm512_permutex2var_ps(result_512_0, shuffle_idx_base1, result_512_8);
        result_512_tmp_2 = _mm512_permutex2var_ps(result_512_1, shuffle_idx_base0, result_512_9);
        result_512_tmp_3 = _mm512_permutex2var_ps(result_512_1, shuffle_idx_base1, result_512_9);
        STORE16_COMPLETE_RESULT(result_512_tmp_0, (&C[ldc*0]))
        STORE16_MASK_COMPLETE_RESULT(result_512_tmp_1, (&C[ldc*0+16]), tail_mask)
        STORE16_COMPLETE_RESULT(result_512_tmp_2, (&C[ldc*1]))
        STORE16_MASK_COMPLETE_RESULT(result_512_tmp_3, (&C[ldc*1+16]), tail_mask)
        result_512_tmp_0 = _mm512_permutex2var_ps(result_512_2, shuffle_idx_base0, result_512_10);
        result_512_tmp_1 = _mm512_permutex2var_ps(result_512_2, shuffle_idx_base1, result_512_10);
        result_512_tmp_2 = _mm512_permutex2var_ps(result_512_3, shuffle_idx_base0, result_512_11);
        result_512_tmp_3 = _mm512_permutex2var_ps(result_512_3, shuffle_idx_base1, result_512_11);
        STORE16_COMPLETE_RESULT(result_512_tmp_0, (&C[ldc*2]))
        STORE16_MASK_COMPLETE_RESULT(result_512_tmp_1, (&C[ldc*2+16]), tail_mask)
        STORE16_COMPLETE_RESULT(result_512_tmp_2, (&C[ldc*3]))
        STORE16_MASK_COMPLETE_RESULT(result_512_tmp_3, (&C[ldc*3+16]), tail_mask)
        result_512_tmp_0 = _mm512_permutex2var_ps(result_512_4, shuffle_idx_base0, result_512_12);
        result_512_tmp_1 = _mm512_permutex2var_ps(result_512_4, shuffle_idx_base1, result_512_12);
        result_512_tmp_2 = _mm512_permutex2var_ps(result_512_5, shuffle_idx_base0, result_512_13);
        result_512_tmp_3 = _mm512_permutex2var_ps(result_512_5, shuffle_idx_base1, result_512_13);
        STORE16_COMPLETE_RESULT(result_512_tmp_0, (&C[ldc*4]))
        STORE16_MASK_COMPLETE_RESULT(result_512_tmp_1, (&C[ldc*4+16]), tail_mask)
        STORE16_COMPLETE_RESULT(result_512_tmp_2, (&C[ldc*5]))
        STORE16_MASK_COMPLETE_RESULT(result_512_tmp_3, (&C[ldc*5+16]), tail_mask)
        result_512_tmp_0 = _mm512_permutex2var_ps(result_512_6, shuffle_idx_base0, result_512_14);
        result_512_tmp_1 = _mm512_permutex2var_ps(result_512_6, shuffle_idx_base1, result_512_14);
        result_512_tmp_2 = _mm512_permutex2var_ps(result_512_7, shuffle_idx_base0, result_512_15);
        result_512_tmp_3 = _mm512_permutex2var_ps(result_512_7, shuffle_idx_base1, result_512_15);
        STORE16_COMPLETE_RESULT(result_512_tmp_0, (&C[ldc*6]))
        STORE16_MASK_COMPLETE_RESULT(result_512_tmp_1, (&C[ldc*6+16]), tail_mask)
        STORE16_COMPLETE_RESULT(result_512_tmp_2, (&C[ldc*7]))
        STORE16_MASK_COMPLETE_RESULT(result_512_tmp_3, (&C[ldc*7+16]), tail_mask)
    } else {
        result_512_tmp_0 = _mm512_permutex2var_ps(result_512_0, shuffle_idx_base0, result_512_8);
        result_512_tmp_1 = _mm512_permutex2var_ps(result_512_0, shuffle_idx_base1, result_512_8);
        result_512_tmp_2 = _mm512_permutex2var_ps(result_512_1, shuffle_idx_base0, result_512_9);
        result_512_tmp_3 = _mm512_permutex2var_ps(result_512_1, shuffle_idx_base1, result_512_9);
        STORE16_COMPLETE_RESULT(result_512_tmp_0, (&C[ldc*0]))
        STORE16_COMPLETE_RESULT(result_512_tmp_1, (&C[ldc*0+16]))
        STORE16_COMPLETE_RESULT(result_512_tmp_2, (&C[ldc*1]))
        STORE16_COMPLETE_RESULT(result_512_tmp_3, (&C[ldc*1+16]))
        result_512_tmp_0 = _mm512_permutex2var_ps(result_512_2, shuffle_idx_base0, result_512_10);
        result_512_tmp_1 = _mm512_permutex2var_ps(result_512_2, shuffle_idx_base1, result_512_10);
        result_512_tmp_2 = _mm512_permutex2var_ps(result_512_3, shuffle_idx_base0, result_512_11);
        result_512_tmp_3 = _mm512_permutex2var_ps(result_512_3, shuffle_idx_base1, result_512_11);
        STORE16_COMPLETE_RESULT(result_512_tmp_0, (&C[ldc*2]))
        STORE16_COMPLETE_RESULT(result_512_tmp_1, (&C[ldc*2+16]))
        STORE16_COMPLETE_RESULT(result_512_tmp_2, (&C[ldc*3]))
        STORE16_COMPLETE_RESULT(result_512_tmp_3, (&C[ldc*3+16]))
        result_512_tmp_0 = _mm512_permutex2var_ps(result_512_4, shuffle_idx_base0, result_512_12);
        result_512_tmp_1 = _mm512_permutex2var_ps(result_512_4, shuffle_idx_base1, result_512_12);
        result_512_tmp_2 = _mm512_permutex2var_ps(result_512_5, shuffle_idx_base0, result_512_13);
        result_512_tmp_3 = _mm512_permutex2var_ps(result_512_5, shuffle_idx_base1, result_512_13);
        STORE16_COMPLETE_RESULT(result_512_tmp_0, (&C[ldc*4]))
        STORE16_COMPLETE_RESULT(result_512_tmp_1, (&C[ldc*4+16]))
        STORE16_COMPLETE_RESULT(result_512_tmp_2, (&C[ldc*5]))
        STORE16_COMPLETE_RESULT(result_512_tmp_3, (&C[ldc*5+16]))
        result_512_tmp_0 = _mm512_permutex2var_ps(result_512_6, shuffle_idx_base0, result_512_14);
        result_512_tmp_1 = _mm512_permutex2var_ps(result_512_6, shuffle_idx_base1, result_512_14);
        result_512_tmp_2 = _mm512_permutex2var_ps(result_512_7, shuffle_idx_base0, result_512_15);
        result_512_tmp_3 = _mm512_permutex2var_ps(result_512_7, shuffle_idx_base1, result_512_15);
        STORE16_COMPLETE_RESULT(result_512_tmp_0, (&C[ldc*6]))
        STORE16_COMPLETE_RESULT(result_512_tmp_1, (&C[ldc*6+16]))
        STORE16_COMPLETE_RESULT(result_512_tmp_2, (&C[ldc*7]))
        STORE16_COMPLETE_RESULT(result_512_tmp_3, (&C[ldc*7+16]))
    }
}

// SBGEMM Kernel for M<=16, N=8, K can be any number, but the processing will take 32 as a base
#ifndef ONE_ALPHA      // ALPHA is not ONE
void sbgemm_block_kernel_16x8x32_alpha(BLASLONG m, BLASLONG k, float alpha, bfloat16 *A, bfloat16 *B, float *C, int ldc)
#else                  // ALPHA is ONE
void sbgemm_block_kernel_16x8x32_one(BLASLONG m, BLASLONG k, float alpha, bfloat16 *A, bfloat16 *B, float *C, int ldc)
#endif
{
    int  SHUFFLE_MAGIC_NO = 0x39;
    BLASLONG tag_k_32x = k & (~31);
    BLASLONG idxB_base = 0;
    BLASLONG width = 32;

#ifndef ONE_ALPHA
    __m512  ALPHAVECTOR = _mm512_set1_ps(alpha);
#endif

    __m512i arrayA_512_0;
    __m512i arrayB_512_0, arrayB_512_1, arrayB_512_2, arrayB_512_3, arrayB_512_4, arrayB_512_5, arrayB_512_6, arrayB_512_7;
    __m512  result_512_0, result_512_1, result_512_2, result_512_3, result_512_4, result_512_5, result_512_6, result_512_7;

    result_512_0  = _mm512_setzero_ps();
    result_512_1  = _mm512_setzero_ps();
    result_512_2  = _mm512_setzero_ps();
    result_512_3  = _mm512_setzero_ps();
    result_512_4  = _mm512_setzero_ps();
    result_512_5  = _mm512_setzero_ps();
    result_512_6  = _mm512_setzero_ps();
    result_512_7  = _mm512_setzero_ps();

    for (BLASLONG idx_k = 0; idx_k < k; idx_k += 32) {
        // Load B with unroll 8
        idxB_base = idx_k << 3;
        arrayB_512_0 = _mm512_loadu_si512(&B[idxB_base + 32*0]);
        arrayB_512_1 = _mm512_loadu_si512(&B[idxB_base + 32*1]);
        arrayB_512_2 = _mm512_loadu_si512(&B[idxB_base + 32*2]);
        arrayB_512_3 = _mm512_loadu_si512(&B[idxB_base + 32*3]);
        arrayB_512_4 = _mm512_loadu_si512(&B[idxB_base + 32*4]);
        arrayB_512_5 = _mm512_loadu_si512(&B[idxB_base + 32*5]);
        arrayB_512_6 = _mm512_loadu_si512(&B[idxB_base + 32*6]);
        arrayB_512_7 = _mm512_loadu_si512(&B[idxB_base + 32*7]);

        if (idx_k == tag_k_32x) {width = k - tag_k_32x;}

        for (BLASLONG idx = 0; idx < width;) {
            // Each two rows are a group for 32-pair bf16 elements
            // Load two rows into a 512 register
            arrayA_512_0 = _mm512_loadu_si512(&A[idx<<4]);

            result_512_0  = _mm512_dpbf16_ps(result_512_0,  (__m512bh) arrayA_512_0, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_0)));
            result_512_1  = _mm512_dpbf16_ps(result_512_1,  (__m512bh) arrayA_512_0, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_1)));
            result_512_2  = _mm512_dpbf16_ps(result_512_2,  (__m512bh) arrayA_512_0, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_2)));
            result_512_3  = _mm512_dpbf16_ps(result_512_3,  (__m512bh) arrayA_512_0, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_3)));
            result_512_4  = _mm512_dpbf16_ps(result_512_4,  (__m512bh) arrayA_512_0, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_4)));
            result_512_5  = _mm512_dpbf16_ps(result_512_5,  (__m512bh) arrayA_512_0, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_5)));
            result_512_6  = _mm512_dpbf16_ps(result_512_6,  (__m512bh) arrayA_512_0, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_6)));
            result_512_7  = _mm512_dpbf16_ps(result_512_7,  (__m512bh) arrayA_512_0, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512_7)));

            arrayB_512_0 = _mm512_shuffle_epi32(arrayB_512_0, SHUFFLE_MAGIC_NO);
            arrayB_512_1 = _mm512_shuffle_epi32(arrayB_512_1, SHUFFLE_MAGIC_NO);
            arrayB_512_2 = _mm512_shuffle_epi32(arrayB_512_2, SHUFFLE_MAGIC_NO);
            arrayB_512_3 = _mm512_shuffle_epi32(arrayB_512_3, SHUFFLE_MAGIC_NO);
            arrayB_512_4 = _mm512_shuffle_epi32(arrayB_512_4, SHUFFLE_MAGIC_NO);
            arrayB_512_5 = _mm512_shuffle_epi32(arrayB_512_5, SHUFFLE_MAGIC_NO);
            arrayB_512_6 = _mm512_shuffle_epi32(arrayB_512_6, SHUFFLE_MAGIC_NO);
            arrayB_512_7 = _mm512_shuffle_epi32(arrayB_512_7, SHUFFLE_MAGIC_NO);

            idx += 2;
            // Every 4 loops we need to switch to next 128 bits of arrayB registers
            if ((idx & (~7)) == idx) {
                arrayB_512_0 = _mm512_shuffle_i32x4(arrayB_512_0, arrayB_512_0, SHUFFLE_MAGIC_NO);
                arrayB_512_1 = _mm512_shuffle_i32x4(arrayB_512_1, arrayB_512_1, SHUFFLE_MAGIC_NO);
                arrayB_512_2 = _mm512_shuffle_i32x4(arrayB_512_2, arrayB_512_2, SHUFFLE_MAGIC_NO);
                arrayB_512_3 = _mm512_shuffle_i32x4(arrayB_512_3, arrayB_512_3, SHUFFLE_MAGIC_NO);
                arrayB_512_4 = _mm512_shuffle_i32x4(arrayB_512_4, arrayB_512_4, SHUFFLE_MAGIC_NO);
                arrayB_512_5 = _mm512_shuffle_i32x4(arrayB_512_5, arrayB_512_5, SHUFFLE_MAGIC_NO);
                arrayB_512_6 = _mm512_shuffle_i32x4(arrayB_512_6, arrayB_512_6, SHUFFLE_MAGIC_NO);
                arrayB_512_7 = _mm512_shuffle_i32x4(arrayB_512_7, arrayB_512_7, SHUFFLE_MAGIC_NO);
            }
        }
    }

    if (m != 16) {
        unsigned short tail_mask_value = (((unsigned short)0xffff) >> (16-m));
        __mmask16 tail_mask = *((__mmask16*) &tail_mask_value);

        result_512_0 = _mm512_shuffle_f32x4(result_512_0, result_512_0, 0xd8);
        result_512_1 = _mm512_shuffle_f32x4(result_512_1, result_512_1, 0xd8);
        result_512_2 = _mm512_shuffle_f32x4(result_512_2, result_512_2, 0xd8);
        result_512_3 = _mm512_shuffle_f32x4(result_512_3, result_512_3, 0xd8);
        STORE16_MASK_COMPLETE_RESULT(result_512_0, (&C[ldc*0]), tail_mask)
        STORE16_MASK_COMPLETE_RESULT(result_512_1, (&C[ldc*1]), tail_mask)
        STORE16_MASK_COMPLETE_RESULT(result_512_2, (&C[ldc*2]), tail_mask)
        STORE16_MASK_COMPLETE_RESULT(result_512_3, (&C[ldc*3]), tail_mask)
        result_512_4 = _mm512_shuffle_f32x4(result_512_4, result_512_4, 0xd8);
        result_512_5 = _mm512_shuffle_f32x4(result_512_5, result_512_5, 0xd8);
        result_512_6 = _mm512_shuffle_f32x4(result_512_6, result_512_6, 0xd8);
        result_512_7 = _mm512_shuffle_f32x4(result_512_7, result_512_7, 0xd8);
        STORE16_MASK_COMPLETE_RESULT(result_512_4, (&C[ldc*4]), tail_mask)
        STORE16_MASK_COMPLETE_RESULT(result_512_5, (&C[ldc*5]), tail_mask)
        STORE16_MASK_COMPLETE_RESULT(result_512_6, (&C[ldc*6]), tail_mask)
        STORE16_MASK_COMPLETE_RESULT(result_512_7, (&C[ldc*7]), tail_mask)
    } else {
        result_512_0 = _mm512_shuffle_f32x4(result_512_0, result_512_0, 0xd8);
        result_512_1 = _mm512_shuffle_f32x4(result_512_1, result_512_1, 0xd8);
        result_512_2 = _mm512_shuffle_f32x4(result_512_2, result_512_2, 0xd8);
        result_512_3 = _mm512_shuffle_f32x4(result_512_3, result_512_3, 0xd8);
        STORE16_COMPLETE_RESULT(result_512_0, (&C[ldc*0]))
        STORE16_COMPLETE_RESULT(result_512_1, (&C[ldc*1]))
        STORE16_COMPLETE_RESULT(result_512_2, (&C[ldc*2]))
        STORE16_COMPLETE_RESULT(result_512_3, (&C[ldc*3]))
        result_512_4 = _mm512_shuffle_f32x4(result_512_4, result_512_4, 0xd8);
        result_512_5 = _mm512_shuffle_f32x4(result_512_5, result_512_5, 0xd8);
        result_512_6 = _mm512_shuffle_f32x4(result_512_6, result_512_6, 0xd8);
        result_512_7 = _mm512_shuffle_f32x4(result_512_7, result_512_7, 0xd8);
        STORE16_COMPLETE_RESULT(result_512_4, (&C[ldc*4]))
        STORE16_COMPLETE_RESULT(result_512_5, (&C[ldc*5]))
        STORE16_COMPLETE_RESULT(result_512_6, (&C[ldc*6]))
        STORE16_COMPLETE_RESULT(result_512_7, (&C[ldc*7]))
    }
}

// SBGEMM Kernel for 16<M<=32, N<8, K can be any number, but the processing will take 32 as a base
#ifndef ONE_ALPHA      // ALPHA is not ONE
void sbgemm_block_kernel_32xNx32_alpha(BLASLONG m, BLASLONG n, BLASLONG k, float alpha, bfloat16 *A, bfloat16 *B, float *C, int ldc)
#else                  // ALPHA is ONE
void sbgemm_block_kernel_32xNx32_one(BLASLONG m, BLASLONG n, BLASLONG k, float alpha, bfloat16 *A, bfloat16 *B, float *C, int ldc)
#endif
{
    int  SHUFFLE_MAGIC_NO = 0x39;
    BLASLONG tag_k_32x = k & (~31);
    BLASLONG idxA_base = 0;
    BLASLONG idxB_base = 0;
    BLASLONG width = 32;

#ifndef ONE_ALPHA
    __m512  ALPHAVECTOR = _mm512_set1_ps(alpha);
#endif

    __m512i arrayA_512[2];
    __m512i arrayB_512[8];
    __m512  result_512[16];
    __m512  result_512_tmp_0, result_512_tmp_1;

    __m512i M512_EPI32_8      = _mm512_set1_epi32(8);
    __m512i shuffle_idx_base0 = _mm512_set_epi32(23, 22, 21, 20, 7, 6, 5, 4, 19, 18, 17, 16, 3, 2, 1, 0);
    __m512i shuffle_idx_base1 = _mm512_add_epi32(shuffle_idx_base0, M512_EPI32_8);

    for (int i = 0; i < 15; i += 2) {
        result_512[i]    = _mm512_setzero_ps();
        result_512[i+1]  = _mm512_setzero_ps();
    }

    for (BLASLONG idx_k = 0; idx_k < k; idx_k += 32) {
        // Load B with unroll n
        for (int i = 0; i < n; i ++) {
            arrayB_512[i] = _mm512_loadu_si512(&B[idxB_base]);
            idxB_base += 32;
        }

        if (idx_k == tag_k_32x) {width = k - tag_k_32x;}

        for (BLASLONG idx = 0; idx < width;) {
            // Each two rows are a group for 32-pair bf16 elements
            idxA_base = idx << 5;
            arrayA_512[0] = _mm512_loadu_si512(&A[idxA_base]);
            arrayA_512[1] = _mm512_loadu_si512(&A[idxA_base + 32]);

            for (int i = 0; i < n; i++) {
                result_512[i]   = _mm512_dpbf16_ps(result_512[i]  , (__m512bh) arrayA_512[0], (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512[i])));
                result_512[i+8] = _mm512_dpbf16_ps(result_512[i+8], (__m512bh) arrayA_512[1], (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512[i])));
                arrayB_512[i]   = _mm512_shuffle_epi32(arrayB_512[i], SHUFFLE_MAGIC_NO);
            }

            idx += 2;
            // Every 4 loops we need to switch to next 128 bits of arrayB registers
            if ((idx & (~7)) == idx) {
                for (int i = 0; i < n; i++) {
                    arrayB_512[i] = _mm512_shuffle_i32x4(arrayB_512[i], arrayB_512[i], SHUFFLE_MAGIC_NO);
                }
            }
        }
    }

    if (m != 32) {
        unsigned short tail_mask_value = (((unsigned short)0xffff) >> (32-m));
        __mmask16 tail_mask = *((__mmask16*) &tail_mask_value);
        for (int i = 0; i < n; i++) {
            result_512_tmp_0 = _mm512_permutex2var_ps(result_512[i], shuffle_idx_base0, result_512[i+8]);
            result_512_tmp_1 = _mm512_permutex2var_ps(result_512[i], shuffle_idx_base1, result_512[i+8]);
            STORE16_COMPLETE_RESULT(result_512_tmp_0, (&C[ldc*i]))
            STORE16_MASK_COMPLETE_RESULT(result_512_tmp_1, (&C[ldc*i+16]), tail_mask)
        }
    } else {
        for (int i = 0; i < n; i++) {
            result_512_tmp_0 = _mm512_permutex2var_ps(result_512[i], shuffle_idx_base0, result_512[i+8]);
            result_512_tmp_1 = _mm512_permutex2var_ps(result_512[i], shuffle_idx_base1, result_512[i+8]);
            STORE16_COMPLETE_RESULT(result_512_tmp_0, (&C[ldc*i]))
            STORE16_COMPLETE_RESULT(result_512_tmp_1, (&C[ldc*i+16]))
        }
    }
}

// SBGEMM Kernel for 16<=M, N<8, K can be any number, but the processing will take 32 as a base
#ifndef ONE_ALPHA      // ALPHA is not ONE
void sbgemm_block_kernel_16xNx32_alpha(BLASLONG m, BLASLONG n, BLASLONG k, float alpha, bfloat16 *A, bfloat16 *B, float *C, int ldc)
#else                  // ALPHA is ONE
void sbgemm_block_kernel_16xNx32_one(BLASLONG m, BLASLONG n, BLASLONG k, float alpha, bfloat16 *A, bfloat16 *B, float *C, int ldc)
#endif
{
    int  SHUFFLE_MAGIC_NO = 0x39;
    BLASLONG tag_k_32x = k & (~31);
    BLASLONG idxB_base = 0;
    BLASLONG width = 32;

#ifndef ONE_ALPHA
    __m512  ALPHAVECTOR = _mm512_set1_ps(alpha);
#endif

    __m512i arrayA_512;
    __m512i arrayB_512[8];
    __m512  result_512[8];

    for (int i = 0; i < 8; i += 2) {
        result_512[i]    = _mm512_setzero_ps();
        result_512[i+1]  = _mm512_setzero_ps();
    }

    for (BLASLONG idx_k = 0; idx_k < k; idx_k += 32) {
        // Load B with unroll n
        for (int i = 0; i < n; i ++) {
            arrayB_512[i] = _mm512_loadu_si512(&B[idxB_base]);
            idxB_base += 32;
        }

        if (idx_k == tag_k_32x) {width = k - tag_k_32x;}

        for (BLASLONG idx = 0; idx < width;) {
            // Each two rows are a group for 32-pair bf16 elements
            // Load two rows into a 512 register
            arrayA_512 = _mm512_loadu_si512(&A[idx<<4]);

            for (int i = 0; i < n; i ++) {
                result_512[i]  = _mm512_dpbf16_ps(result_512[i],  (__m512bh) arrayA_512, (__m512bh) _mm512_broadcastd_epi32(_mm512_castsi512_si128(arrayB_512[i])));
                arrayB_512[i] = _mm512_shuffle_epi32(arrayB_512[i], SHUFFLE_MAGIC_NO);
            }

            idx += 2;
            // Every 4 loops we need to switch to next 128 bits of arrayB registers
            if ((idx & (~7)) == idx) {
                for (int i = 0; i < n; i++) {
                    arrayB_512[i] = _mm512_shuffle_i32x4(arrayB_512[i], arrayB_512[i], SHUFFLE_MAGIC_NO);
                }
            }
        }
    }

    if (m != 16) {
        unsigned short tail_mask_value = (((unsigned short)0xffff) >> (16-m));
        __mmask16 tail_mask = *((__mmask16*) &tail_mask_value);
        for (int i = 0; i < n; i++) {
            result_512[i] = _mm512_shuffle_f32x4(result_512[i], result_512[i], 0xd8);
            STORE16_MASK_COMPLETE_RESULT(result_512[i], (&C[ldc*i]), tail_mask)
        }
    } else {
        for (int i = 0; i < n; i++) {
            result_512[i] = _mm512_shuffle_f32x4(result_512[i], result_512[i], 0xd8);
            STORE16_COMPLETE_RESULT(result_512[i], (&C[ldc*i]))
        }
    }
}
#ifndef ONE_ALPHA      // ALPHA is not ONE
void sbgemm_blocking_kernel_2_alpha(blasint M, blasint N, blasint K, float alpha, bfloat16 *A, blasint lda, bfloat16 *B, blasint ldb, float *C, blasint ldc, bfloat16 * block_A, bfloat16 * block_B)
#else                  // ALPHA is ONE
void sbgemm_blocking_kernel_2_one(blasint M, blasint N, blasint K, float alpha, bfloat16 *A, blasint lda, bfloat16 *B, blasint ldb, float *C, blasint ldc, bfloat16 * block_A, bfloat16 * block_B)
#endif
{
    BLASLONG m_step, n_step, k_step, k_step_round32;
    BLASLONG tag_m_Nx = M & (~(BF16_BLOCK_THRES_M-1));

    BLASLONG n_from, n_to;
    BLASLONG tag_n_Nx;

    n_from = 0;
    n_to = (BF16_BLOCK_THRES_N > N) ? N : BF16_BLOCK_THRES_N;
    tag_n_Nx = n_to & (~(BF16_BLOCK_STEP_N-1));

    k_step = (K > BF16_BLOCK_THRES_K) ? BF16_BLOCK_THRES_K : K;
    k_step_round32 = k_step & (~31);
    k_step_round32 = (k_step > k_step_round32) ? (k_step_round32 + 32) : k_step_round32;

    if (M >= BF16_BLOCK_THRES_M) {
        while (n_from < N) {
            for (BLASLONG idx_k = 0; idx_k < K;) {
                // Use Kx32 kernel when BF16_BLOCK_THRES_M==32, Kx16 kernel when BF16_BLOCK_THRES_M==16, ...
                COL_MAJOR_INCOPY_KERNEL_Kx32(k_step, &A(idx_k, 0), lda, block_A);
                // TODO: MT
                for (BLASLONG idx_n = n_from; idx_n < tag_n_Nx; idx_n += BF16_BLOCK_STEP_N) {
                    // Use 8x32 kernel when BF16_BLOCK_THRES_N==8, 4x32 kernel when BF16_BLOCK_THRES_N==4, ...
                    COL_MAJOR_ONCOPY_KERNEL_8x32(k_step, &B(idx_n, idx_k), ldb, block_B + (idx_n-n_from)*k_step_round32);
                    SBGEMM_BLOCK_KERNEL_32x8x32(32, k_step, alpha, block_A, block_B + (idx_n-n_from)*k_step_round32, &C(idx_n, 0), ldc);
                }

                if (tag_n_Nx != n_to) {
                    n_step = n_to - tag_n_Nx;
                    COL_MAJOR_ONCOPY_KERNEL_Nx32(n_step, k_step, &B(tag_n_Nx, idx_k), ldb, block_B + (tag_n_Nx-n_from)*k_step_round32);
                    SBGEMM_BLOCK_KERNEL_32xNx32(32, n_step, k_step, alpha, block_A, block_B + (tag_n_Nx-n_from)*k_step_round32, &C(tag_n_Nx, 0), ldc);
                }

                for (BLASLONG idx_m = BF16_BLOCK_THRES_M; idx_m < tag_m_Nx; idx_m += BF16_BLOCK_THRES_M) {
                    COL_MAJOR_INCOPY_KERNEL_Kx32(k_step, &A(idx_k, idx_m), lda, block_A);
                    for (BLASLONG idx_n = n_from; idx_n < tag_n_Nx; idx_n += BF16_BLOCK_STEP_N) {
                        SBGEMM_BLOCK_KERNEL_32x8x32(32, k_step, alpha, block_A, block_B + (idx_n-n_from)*k_step_round32, &C(idx_n, idx_m), ldc);
                    }

                    if (tag_n_Nx != n_to) {
                        n_step = n_to - tag_n_Nx;
                        SBGEMM_BLOCK_KERNEL_32xNx32(32, n_step, k_step, alpha, block_A, block_B + (tag_n_Nx-n_from)*k_step_round32, &C(tag_n_Nx, idx_m), ldc);
                    }
                }

                if (tag_m_Nx != M) {
                    m_step = M - tag_m_Nx;
                    if (m_step > 16) {
                        COL_MAJOR_INCOPY_KERNEL_Kx32m(k_step, m_step, &A(idx_k, tag_m_Nx), lda, block_A);
                        for (BLASLONG idx_n = n_from; idx_n < tag_n_Nx; idx_n += BF16_BLOCK_STEP_N) {
                            SBGEMM_BLOCK_KERNEL_32x8x32(m_step, k_step, alpha, block_A, block_B + (idx_n-n_from)*k_step_round32, &C(idx_n, tag_m_Nx), ldc);
                        }

                        if (tag_n_Nx != n_to) {
                            n_step = n_to - tag_n_Nx;
                            SBGEMM_BLOCK_KERNEL_32xNx32(m_step, n_step, k_step, alpha, block_A, block_B + (tag_n_Nx-n_from)*k_step_round32, &C(tag_n_Nx, tag_m_Nx), ldc);
                        }
                    } else if (m_step == 16) {
                        COL_MAJOR_INCOPY_KERNEL_Kx16(k_step, m_step, &A(idx_k, tag_m_Nx), lda, block_A);
                        for (BLASLONG idx_n = n_from; idx_n < tag_n_Nx; idx_n += BF16_BLOCK_STEP_N) {
                            SBGEMM_BLOCK_KERNEL_16x8x32(m_step, k_step, alpha, block_A, block_B + (idx_n-n_from)*k_step_round32, &C(idx_n, tag_m_Nx), ldc);
                        }

                        if (tag_n_Nx != n_to) {
                            n_step = n_to - tag_n_Nx;
                            SBGEMM_BLOCK_KERNEL_16xNx32(m_step, n_step, k_step, alpha, block_A, block_B + (tag_n_Nx-n_from)*k_step_round32, &C(tag_n_Nx, tag_m_Nx), ldc);
                        }
                    } else {
                        COL_MAJOR_INCOPY_KERNEL_Kx16m(k_step, m_step, &A(idx_k, tag_m_Nx), lda, block_A);
                        for (BLASLONG idx_n = n_from; idx_n < tag_n_Nx; idx_n += BF16_BLOCK_STEP_N) {
                            SBGEMM_BLOCK_KERNEL_16x8x32(m_step, k_step, alpha, block_A, block_B + (idx_n-n_from)*k_step_round32, &C(idx_n, tag_m_Nx), ldc);
                        }

                        if (tag_n_Nx != n_to) {
                            n_step = n_to - tag_n_Nx;
                            SBGEMM_BLOCK_KERNEL_16xNx32(m_step, n_step, k_step, alpha, block_A, block_B + (tag_n_Nx-n_from)*k_step_round32, &C(tag_n_Nx, tag_m_Nx), ldc);
                        }
                    }
                }

                idx_k += k_step;
                k_step = K - idx_k;
                k_step = (k_step > BF16_BLOCK_THRES_K) ? BF16_BLOCK_THRES_K : k_step;
                k_step_round32 = k_step & (~31);
                k_step_round32 = (k_step > k_step_round32) ? (k_step_round32 + 32) : k_step_round32;
            }

            n_from = n_to;
            n_to += BF16_BLOCK_THRES_N;
            n_to = (n_to > N) ? N : n_to;
            tag_n_Nx = n_to & (~(BF16_BLOCK_STEP_N-1));           
        }
    } else {
        m_step = M - tag_m_Nx;
        while (n_from < N) {
            for (BLASLONG idx_k = 0; idx_k < K;) {
                // Use Kx32 kernel when BF16_BLOCK_THRES_M==32, Kx16 kernel when BF16_BLOCK_THRES_M==16, ...
                COL_MAJOR_INCOPY_KERNEL_Kx32m(k_step, m_step, &A(idx_k, 0), lda, block_A);
                // TODO: MT
                for (BLASLONG idx_n = n_from; idx_n < tag_n_Nx; idx_n += BF16_BLOCK_STEP_N) {
                    // Use 8x32 kernel when BF16_BLOCK_THRES_N==8, 4x32 kernel when BF16_BLOCK_THRES_N==4, ...
                    COL_MAJOR_ONCOPY_KERNEL_8x32(k_step, &B(idx_n, idx_k), ldb, block_B + (idx_n-n_from)*k_step_round32);
                    SBGEMM_BLOCK_KERNEL_32x8x32(m_step, k_step, alpha, block_A, block_B + (idx_n-n_from)*k_step_round32, &C(idx_n, 0), ldc);
                }

                if (tag_n_Nx != n_to) {
                    n_step = n_to - tag_n_Nx;
                    COL_MAJOR_ONCOPY_KERNEL_Nx32(n_step, k_step, &B(tag_n_Nx, idx_k), ldb, block_B + (tag_n_Nx-n_from)*k_step_round32);
                    SBGEMM_BLOCK_KERNEL_32xNx32(m_step, n_step, k_step, alpha, block_A, block_B + (tag_n_Nx-n_from)*k_step_round32, &C(tag_n_Nx, 0), ldc);
                }

                idx_k += k_step;
                k_step = K - idx_k;
                k_step = (k_step > BF16_BLOCK_THRES_K) ? BF16_BLOCK_THRES_K : k_step;
                k_step_round32 = k_step & (~31);
                k_step_round32 = (k_step > k_step_round32) ? (k_step_round32 + 32) : k_step_round32;
            }
            n_from = n_to;
            n_to += BF16_BLOCK_THRES_N;
            n_to = (n_to > N) ? N : n_to;
            tag_n_Nx = n_to & (~(BF16_BLOCK_STEP_N-1));
        }
    }
}

#ifndef ONE_ALPHA      // ALPHA is not ONE
void sbgemm_internal_kernel_alpha(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K,
		 OPENBLAS_CONST float alpha, OPENBLAS_CONST bfloat16 *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST bfloat16 *B, OPENBLAS_CONST blasint ldb, float *C, OPENBLAS_CONST blasint ldc)
#else                  // ALPHA is ONE
void sbgemm_internal_kernel_one(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K,
		 OPENBLAS_CONST float alpha, OPENBLAS_CONST bfloat16 *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST bfloat16 *B, OPENBLAS_CONST blasint ldb, float *C, OPENBLAS_CONST blasint ldc)
#endif
{
    bfloat16 block_A[BF16_BLOCK_THRES_K * BF16_BLOCK_THRES_M];
    bfloat16 block_B[BF16_BLOCK_THRES_N * BF16_BLOCK_THRES_K];

    // TODO: assume no trans for both A and B, to complement these scenarios later
    if (Order == CblasColMajor) {
        SBGEMM_BLOCKING_KERNEL_2(M, N, K, alpha, A, lda, B, ldb, C, ldc, block_A, block_B);
    } else {
        
    }
}