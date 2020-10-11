#define V_SIMD 512
#define V_SIMD_F64 1
/***************************
 * Data Type
 ***************************/
typedef __m512  v_f32;
#define v_nlanes_f32 16
/***************************
 * Arithmetic
 ***************************/
#define v_add_f32 _mm512_add_ps
#define v_mul_f32 _mm512_mul_ps
// multiply and add, a*b + c
#define v_muladd_f32 _mm512_fmadd_ps

BLAS_FINLINE float v_sum_f32(v_f32 a)
{
    __m512 h64 = _mm512_shuffle_f32x4(a, a, _MM_SHUFFLE(3, 2, 3, 2));
    __m512 sum32 = _mm512_add_ps(a, h64);
    __m512 h32 = _mm512_shuffle_f32x4(sum32, sum32, _MM_SHUFFLE(1, 0, 3, 2));
    __m512 sum16 = _mm512_add_ps(sum32, h32);
    __m512 h16 = _mm512_permute_ps(sum16, _MM_SHUFFLE(1, 0, 3, 2));
    __m512 sum8 = _mm512_add_ps(sum16, h16);
    __m512 h4 = _mm512_permute_ps(sum8, _MM_SHUFFLE(2, 3, 0, 1));
    __m512 sum4 = _mm512_add_ps(sum8, h4);
    return _mm_cvtss_f32(_mm512_castps512_ps128(sum4));
}
/***************************
 * memory
 ***************************/
// unaligned load
#define v_loadu_f32(PTR) _mm512_loadu_ps((const __m512*)(PTR))
#define v_storeu_f32 _mm512_storeu_ps
#define v_setall_f32(VAL) _mm512_set1_ps(VAL)
#define v_zero_f32 _mm512_setzero_ps
