#define V_SIMD 512
#define V_SIMD_F64 1
/*
Data Type
*/
typedef __m512  v_f32;
#define v_nlanes_f32 16
/*
arithmetic
*/
#define v_add_f32 _mm512_add_ps
#define v_mul_f32 _mm512_mul_ps
// multiply and add, a*b + c
#define v_muladd_f32 _mm512_fmadd_ps
/*
memory
*/
// unaligned load
#define v_loadu_f32(PTR) _mm512_loadu_ps((const __m512*)(PTR))
#define v_storeu_f32 _mm512_storeu_ps
#define v_setall_f32(VAL) _mm512_set1_ps(VAL)
