/* need a new enough GCC for avx512 support */
#ifdef __NVCOMPILER
#define NVCOMPVERS ( __NVCOMPILER_MAJOR__ * 100 + __NVCOMPILER_MINOR__ )
#endif
#if (( defined(__GNUC__)  && __GNUC__   > 6 && defined(__AVX512CD__)) || (defined(__clang__) && __clang_major__ >= 9)) || (defined(__NVCOMPILER) && NVCOMPVERS >= 2203 )

#define HAVE_SASUM_KERNEL 1

#ifndef ABS_K
#define ABS_K(a) ((a) > 0 ? (a) : (-(a)))
#endif

#include <immintrin.h>
#include <stdint.h>

static FLOAT sasum_kernel(BLASLONG n, FLOAT *x1)
{
    BLASLONG i = 0;
    FLOAT sumf = 0.0;

    BLASLONG tail_index_SSE = n&(~7);
    BLASLONG tail_index_AVX512 = n&(~255);

    /* Deliberately no peeling to an alignment boundary: peeling makes the
       grouping of the sum into accumulators - and hence the rounding of the
       result - depend on the buffer address. Unaligned loads keep the result
       a function of n alone, and cost the same as aligned loads when the
       address happens to be aligned. */
    if (n >= 256) {
        __m512 accum_0, accum_1, accum_2, accum_3;
        accum_0 = _mm512_setzero_ps();
        accum_1 = _mm512_setzero_ps();
        accum_2 = _mm512_setzero_ps();
        accum_3 = _mm512_setzero_ps();

        for (i = 0; i < tail_index_AVX512; i += 64) {
            accum_0 += _mm512_abs_ps(_mm512_loadu_ps(&x1[i + 0]));
            accum_1 += _mm512_abs_ps(_mm512_loadu_ps(&x1[i +16]));
            accum_2 += _mm512_abs_ps(_mm512_loadu_ps(&x1[i +32]));
            accum_3 += _mm512_abs_ps(_mm512_loadu_ps(&x1[i +48]));
        }

        accum_0 = accum_0 + accum_1 + accum_2 + accum_3;
        sumf += _mm512_reduce_add_ps(accum_0);
    }

    if (n >= 8) {
        __m128 accum_20, accum_21;
        accum_20 = _mm_setzero_ps();
        accum_21 = _mm_setzero_ps();

        __m128i abs_mask2 = _mm_set1_epi32(0x7fffffff);
        for (i = tail_index_AVX512; i < tail_index_SSE; i += 8) {
            accum_20 += (__m128)_mm_and_si128(_mm_loadu_si128((__m128i*)&x1[i + 0]), abs_mask2);
            accum_21 += (__m128)_mm_and_si128(_mm_loadu_si128((__m128i*)&x1[i + 4]), abs_mask2);
        }
        
        accum_20 += accum_21;
        accum_20 = _mm_hadd_ps(accum_20, accum_20);
        accum_20 = _mm_hadd_ps(accum_20, accum_20);
        
        sumf += accum_20[0];
    }

    for (i = tail_index_SSE; i < n; i++) {
        sumf += ABS_K(x1[i]);
    }

    return sumf;
}
#endif
