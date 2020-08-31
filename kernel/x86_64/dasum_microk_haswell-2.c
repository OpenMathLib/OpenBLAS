#if (( defined(__GNUC__)  && __GNUC__   > 6 ) || (defined(__clang__) && __clang_major__ >= 6)) && defined(__AVX2__)

#define HAVE_KERNEL_16 1

#include <immintrin.h>

static FLOAT dasum_kernel_16(BLASLONG n, FLOAT *x1)
{
    BLASLONG i = 0;
    __m256d accum_0, accum_1, accum_2, accum_3;

    accum_0 = _mm256_setzero_pd();
    accum_1 = _mm256_setzero_pd();
    accum_2 = _mm256_setzero_pd();
    accum_3 = _mm256_setzero_pd();

     __m256i abs_mask = _mm256_set1_epi64x(0x7fffffffffffffff);
    for (; i < n; i += 16) {
        accum_0 += (__m256d)_mm256_and_si256(_mm256_loadu_si256(&x1[i+ 0]), abs_mask);
        accum_1 += (__m256d)_mm256_and_si256(_mm256_loadu_si256(&x1[i+ 4]), abs_mask);
        accum_2 += (__m256d)_mm256_and_si256(_mm256_loadu_si256(&x1[i+ 8]), abs_mask);
        accum_3 += (__m256d)_mm256_and_si256(_mm256_loadu_si256(&x1[i+12]), abs_mask);
    }

    accum_0 = accum_0 + accum_1 + accum_2 + accum_3;

    __m128d half_accum0;
    half_accum0 = _mm_add_pd(_mm256_extractf128_pd(accum_0, 0), _mm256_extractf128_pd(accum_0, 1));

    half_accum0 = _mm_hadd_pd(half_accum0, half_accum0);

    return half_accum0[0];

}
#endif
