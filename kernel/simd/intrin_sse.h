#define V_SIMD 128
#define V_SIMD_F64 1
/***************************
 * Data Type
 ***************************/
#ifdef DOUBLE
typedef __m128d  v_f32;
#else
typedef __m128  v_f32;
#endif

#define v_nlanes_f32 4
/***************************
 * Arithmetic
 ***************************/
#ifdef DOUBLE
#define v_add_f32 _mm_add_pd
#define v_mul_f32 _mm_mul_pd
#else
#define v_add_f32 _mm_add_ps
#define v_mul_f32 _mm_mul_ps
#endif
#ifdef HAVE_FMA3
    // multiply and add, a*b + c
#ifdef DOUBLE
    #define v_muladd_f32 _mm_fmadd_pd
#else
	#define v_muladd_f32 _mm_fmadd_ps
#endif
#elif defined(HAVE_FMA4)
    // multiply and add, a*b + c
    #ifdef DOUBLE
        #define v_muladd_f32 _mm_macc_pd
    #else
	#define v_muladd_f32 _mm_macc_ps
    #endif
#else
    // multiply and add, a*b + c
    BLAS_FINLINE v_f32 v_muladd_f32(v_f32 a, v_f32 b, v_f32 c)
    { return v_add_f32(v_mul_f32(a, b), c); }
#endif // HAVE_FMA3

// Horizontal add: Calculates the sum of all vector elements.
#ifdef DOUBLE
BLAS_FINLINE double v_sum_f32(__m128d a)
{
#ifdef HAVE_SSE3
    __m128d sum_halves = _mm_hadd_pd(a, a);
    return _mm_cvtsd_f64(_mm_hadd_pd(sum_halves, sum_halves));
#else
    __m128d t1 = _mm_movehl_pd(a, a);
    __m128d t2 = _mm_add_pd(a, t1);
    __m128d t3 = _mm_shuffle_pd(t2, t2, 1);
    __m128d t4 = _mm_add_ss(t2, t3);
    return _mm_cvtsd_f64(t4);
#endif
}
#else
// Horizontal add: Calculates the sum of all vector elements.
BLAS_FINLINE float v_sum_f32(__m128 a)
{
#ifdef HAVE_SSE3
    __m128 sum_halves = _mm_hadd_ps(a, a);
    return _mm_cvtss_f32(_mm_hadd_ps(sum_halves, sum_halves));
#else
    __m128 t1 = _mm_movehl_ps(a, a);
    __m128 t2 = _mm_add_ps(a, t1);
    __m128 t3 = _mm_shuffle_ps(t2, t2, 1);
    __m128 t4 = _mm_add_ss(t2, t3);
    return _mm_cvtss_f32(t4);
#endif
}
#endif
/***************************
 * memory
 ***************************/
// unaligned load
#ifdef DOUBLE
#define v_loadu_f32 _mm_loadu_pd
#define v_storeu_f32 _mm_storeu_pd
#define v_setall_f32(VAL) _mm_set1_pd(VAL)
#define v_zero_f32 _mm_setzero_pd
#else
#define v_loadu_f32 _mm_loadu_ps
#define v_storeu_f32 _mm_storeu_ps
#define v_setall_f32(VAL) _mm_set1_ps(VAL)
#define v_zero_f32 _mm_setzero_ps
#endif
