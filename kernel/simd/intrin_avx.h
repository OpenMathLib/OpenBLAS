#define V_SIMD 256
#define V_SIMD_F64 1
/*
Data Type
*/
typedef __m256  v_f32;
#define v_nlanes_f32 8
/*
arithmetic
*/
#define v_add_f32 _mm256_add_ps
#define v_mul_f32 _mm256_mul_ps
/*
memory
*/
// unaligned load
#define v_loadu_f32 _mm256_loadu_ps
#define v_storeu_f32 _mm256_storeu_ps
#define v_setall_f32(VAL) _mm256_set1_ps(VAL)