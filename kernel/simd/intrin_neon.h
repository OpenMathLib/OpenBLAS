#define V_SIMD 128
#ifdef __aarch64__
    #define V_SIMD_F64 1
#else
    #define V_SIMD_F64 0
#endif
/***************************
 * Data Type
 ***************************/
typedef float32x4_t v_f32;
#define v_nlanes_f32 4
/***************************
 * Arithmetic
 ***************************/
#define v_add_f32 vaddq_f32
#define v_mul_f32 vmulq_f32

// FUSED F32
#ifdef HAVE_VFPV4 // FMA
    // multiply and add, a*b + c
    BLAS_FINLINE v_f32 v_muladd_f32(v_f32 a, v_f32 b, v_f32 c)
    { return vfmaq_f32(c, a, b); }
#else
    // multiply and add, a*b + c
    BLAS_FINLINE v_f32 v_muladd_f32(v_f32 a, v_f32 b, v_f32 c)
    { return vmlaq_f32(c, a, b); }
#endif

// Horizontal add: Calculates the sum of all vector elements.
BLAS_FINLINE float v_sum_f32(float32x4_t a)
{
    float32x2_t r = vadd_f32(vget_high_f32(a), vget_low_f32(a));
    return vget_lane_f32(vpadd_f32(r, r), 0);
}
/***************************
 * memory
 ***************************/
// unaligned load
#define v_loadu_f32(a) vld1q_f32((const float*)a)
#define v_storeu_f32 vst1q_f32
#define v_setall_f32(VAL) vdupq_n_f32(VAL)
#define v_zero_f32() vdupq_n_f32(0.0f)