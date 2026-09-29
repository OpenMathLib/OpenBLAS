#include <wasm_simd128.h>

#define V_SIMD 128
#define V_SIMD_F64 1

/***************************
 * Data Type
 ***************************/
typedef v128_t v_f32;
typedef v128_t v_f64;
#define v_nlanes_f32 4
#define v_nlanes_f64 2

/***************************
 * Arithmetic
 ***************************/
#define v_add_f32 wasm_f32x4_add
#define v_add_f64 wasm_f64x2_add
#define v_sub_f32 wasm_f32x4_sub
#define v_sub_f64 wasm_f64x2_sub
#define v_mul_f32 wasm_f32x4_mul
#define v_mul_f64 wasm_f64x2_mul

#ifdef __wasm_relaxed_simd__
BLAS_FINLINE v_f32 v_muladd_f32(v_f32 a, v_f32 b, v_f32 c)
{ return wasm_f32x4_relaxed_madd(a, b, c); }

BLAS_FINLINE v_f64 v_muladd_f64(v_f64 a, v_f64 b, v_f64 c)
{ return wasm_f64x2_relaxed_madd(a, b, c); }
#else
BLAS_FINLINE v_f32 v_muladd_f32(v_f32 a, v_f32 b, v_f32 c)
{ return v_add_f32(v_mul_f32(a, b), c); }

BLAS_FINLINE v_f64 v_muladd_f64(v_f64 a, v_f64 b, v_f64 c)
{ return v_add_f64(v_mul_f64(a, b), c); }
#endif

BLAS_FINLINE v_f32 v_mulsub_f32(v_f32 a, v_f32 b, v_f32 c)
{ return v_sub_f32(v_mul_f32(a, b), c); }

BLAS_FINLINE v_f64 v_mulsub_f64(v_f64 a, v_f64 b, v_f64 c)
{ return v_sub_f64(v_mul_f64(a, b), c); }

BLAS_FINLINE v_f64 v_cvt_f32_f64_lo(v_f32 a)
{ return wasm_f64x2_promote_low_f32x4(a); }

BLAS_FINLINE v_f64 v_cvt_f32_f64_hi(v_f32 a)
{
    v128_t hi = wasm_i32x4_shuffle(a, a, 2, 3, 0, 1);
    return wasm_f64x2_promote_low_f32x4(hi);
}

/***************************
 * reduction
 ***************************/
BLAS_FINLINE float v_sum_f32(v_f32 a)
{
    return wasm_f32x4_extract_lane(a, 0)
         + wasm_f32x4_extract_lane(a, 1)
         + wasm_f32x4_extract_lane(a, 2)
         + wasm_f32x4_extract_lane(a, 3);
}

BLAS_FINLINE double v_sum_f64(v_f64 a)
{
    return wasm_f64x2_extract_lane(a, 0)
         + wasm_f64x2_extract_lane(a, 1);
}

/***************************
 * memory
 ***************************/
#define v_loadu_f32(a) wasm_v128_load((const float*)a)
#define v_loadu_f64(a) wasm_v128_load((const double*)a)
#define v_storeu_f32(a, v) wasm_v128_store((float*)a, v)
#define v_storeu_f64(a, v) wasm_v128_store((double*)a, v)
#define v_setall_f32(VAL) wasm_f32x4_splat(VAL)
#define v_setall_f64(VAL) wasm_f64x2_splat(VAL)
#define v_zero_f32() wasm_f32x4_splat(0.0f)
#define v_zero_f64() wasm_f64x2_splat(0.0)
