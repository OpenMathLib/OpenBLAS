/* need a new enough GCC for avx512 support */
#if (( defined(__GNUC__)  && __GNUC__   > 6 && defined(__AVX2__)) || (defined(__clang__) && __clang_major__ >= 6))

#if defined(__AVX512CD__)
#define HAVE_KERNEL_16 1

#include <immintrin.h>

static FLOAT dasum_kernel_16(BLASLONG n, FLOAT *x1)
{
    BLASLONG i = 0;

    __m512d accum_0, accum_1;

    accum_0 = _mm512_setzero_pd();
    accum_1 = _mm512_setzero_pd();

    for (; i < n; i += 16) {
        accum_0 += _mm512_abs_pd(_mm512_loadu_pd(&x1[i+ 0]));
        accum_1 += _mm512_abs_pd(_mm512_loadu_pd(&x1[i+ 8]));
    }

    accum_0 += accum_1;
    return _mm512_reduce_add_pd(accum_0);
}
#endif
#endif
