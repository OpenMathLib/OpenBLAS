#define V_SIMD 128
#define V_SIMD_F64 1
/*
Data Type
*/
typedef __m128  v_f32;
#define v_nlanes_f32 4
/*
arithmetic
*/
#define v_add_f32 _mm_add_ps
#define v_mul_f32 _mm_mul_ps
/*
memory
*/
// unaligned load
#define v_loadu_f32 _mm_loadu_ps
#define v_storeu_f32 _mm_storeu_ps
#define v_setall_f32(VAL) _mm_set1_ps(VAL)