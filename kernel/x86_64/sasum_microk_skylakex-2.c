/* need a new enough GCC for avx512 support */
#if (( defined(__GNUC__)  && __GNUC__   > 6 && defined(__AVX2__)) || (defined(__clang__) && __clang_major__ >= 6))

#if defined(__AVX512CD__)
#define HAVE_KERNEL_32 1

#include <immintrin.h>

static FLOAT sasum_kernel_32(BLASLONG n, FLOAT *x1)
{
    BLASLONG i = 0;

    __m512 accum_0, accum_1;

    accum_0 = _mm512_setzero_ps();
    accum_1 = _mm512_setzero_ps();

    for (; i < n; i += 32) {
        accum_0 += _mm512_abs_ps(_mm512_loadu_ps(&x1[i+ 0]));
        accum_1 += _mm512_abs_ps(_mm512_loadu_ps(&x1[i+ 16]));
    }

    accum_0 += accum_1;
    return _mm512_reduce_add_ps(accum_0);
}
#endif
#endif
