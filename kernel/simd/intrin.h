#ifndef _INTRIN_H_
#define _INTRIN_H_

#ifdef __cplusplus
extern "C" {
#endif
// include head
/** SSE **/
#ifdef HAVE_SSE
#include <xmmintrin.h>
#endif
/** SSE2 **/
#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
/** SSE3 **/
#ifdef HAVE_SSE3
#include <pmmintrin.h>
#endif
/** SSSE3 **/
#ifdef HAVE_SSSE3
#include <tmmintrin.h>
#endif
/** SSE41 **/
#ifdef HAVE_SSE4_1
#include <smmintrin.h>
#endif

/** AVX **/
#ifdef HAVE_AVX
#include <immintrin.h>
#endif

// distribute
#if defined(HAVE_AVX512VL) || defined(HAVE_AVX512BF16)
#include "intrin_avx512.h"
#elif defined(HAVE_AVX2)
#include "intrin_avx.h"
#elif defined(HAVE_SSE2)
#include "intrin_sse.h"
#endif

#ifndef V_SIMD
    #define V_SIMD 0
    #define V_SIMD_F64 0
#endif

#ifdef __cplusplus
}
#endif
#endif // _INTRIN_H_
