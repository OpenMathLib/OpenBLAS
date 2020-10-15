#define V_SIMD 128
#define V_SIMD_F64 1
/***************************
 * Data Type
 ***************************/
typedef __m128  v_f32;
#define v_nlanes_f32 4
/***************************
 * Arithmetic
 ***************************/
#define v_add_f32 _mm_add_ps
#define v_mul_f32 _mm_mul_ps
#ifdef HAVE_FMA3
    // multiply and add, a*b + c
    #define v_muladd_f32 _mm_fmadd_ps
#elif defined(HAVE_FMA4)
    // multiply and add, a*b + c
    #define v_muladd_f32 _mm_macc_ps
#else
    // multiply and add, a*b + c
    BLAS_FINLINE v_f32 v_muladd_f32(v_f32 a, v_f32 b, v_f32 c)
    { return v_add_f32(v_mul_f32(a, b), c); }
#endif // HAVE_FMA3

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
/***************************
 * memory
 ***************************/
// unaligned load
#define v_loadu_f32 _mm_loadu_ps
#define v_storeu_f32 _mm_storeu_ps
#define v_setall_f32(VAL) _mm_set1_ps(VAL)
#define v_zero_f32 _mm_setzero_ps