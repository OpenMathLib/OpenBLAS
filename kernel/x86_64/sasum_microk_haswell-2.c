#if (( defined(__GNUC__)  && __GNUC__   > 6 && defined(__AVX2__)) || (defined(__clang__) && __clang_major__ >= 6))

#define HAVE_KERNEL_32 1

#include <immintrin.h>

static FLOAT sasum_kernel_32(BLASLONG n, FLOAT *x1)
{
    BLASLONG i = 0;
    __m256 accum_0, accum_1, accum_2, accum_3;

    accum_0 = _mm256_setzero_ps();
    accum_1 = _mm256_setzero_ps();
    accum_2 = _mm256_setzero_ps();
    accum_3 = _mm256_setzero_ps();

    __m256i abs_mask = _mm256_set1_epi32(0x7fffffff);
    for (; i < n; i += 32) {
        accum_0 += (__m256)_mm256_and_si256(_mm256_loadu_si256(&x1[i+ 0]), abs_mask);
        accum_1 += (__m256)_mm256_and_si256(_mm256_loadu_si256(&x1[i+ 8]), abs_mask);
        accum_2 += (__m256)_mm256_and_si256(_mm256_loadu_si256(&x1[i+16]), abs_mask);
        accum_3 += (__m256)_mm256_and_si256(_mm256_loadu_si256(&x1[i+24]), abs_mask);
    }

    accum_0 = accum_0 + accum_1 + accum_2 + accum_3;

    __m128 half_accum0;
    half_accum0 = _mm_add_ps(_mm256_extractf128_ps(accum_0, 0), _mm256_extractf128_ps(accum_0, 1));

    half_accum0 = _mm_hadd_ps(half_accum0, half_accum0);
    half_accum0 = _mm_hadd_ps(half_accum0, half_accum0);

    return half_accum0[0];

}
#endif
