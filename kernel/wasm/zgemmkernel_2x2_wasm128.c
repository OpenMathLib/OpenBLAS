#include "common.h"
#if defined(__wasm_simd128__)
#include <wasm_simd128.h>
#if defined(__wasm_relaxed_simd__)
#define MADD_F32(a, b, c) wasm_f32x4_relaxed_madd((a), (b), (c))
#define MADD_F64(a, b, c) wasm_f64x2_relaxed_madd((a), (b), (c))
#else
#define MADD_F32(a, b, c) wasm_f32x4_add((c), wasm_f32x4_mul((a), (b)))
#define MADD_F64(a, b, c) wasm_f64x2_add((c), wasm_f64x2_mul((a), (b)))
#endif
#endif
/********************************
  SIMD128 2x2 complex GEMM (CGEMM/ZGEMM). Remainders stay scalar.
  ADD1 a*c
*********************************/
int CNAME(BLASLONG bm,BLASLONG bn,BLASLONG bk,FLOAT alphar,FLOAT alphai,FLOAT* ba,FLOAT* bb,FLOAT* C,BLASLONG ldc
#ifdef	TRMMKERNEL
		, BLASLONG offset
#endif
		)
{
   BLASLONG i,j,k;
   FLOAT *C0,*C1,*ptrba,*ptrbb;
   FLOAT res0,res1,res2,res3,res4,res5,res6,res7,load0,load1,load2,load3,load4,load5,load6,load7,load8,load9,load10,load11,load12,load13,load14,load15;
   for (j=0; j<bn/2; j+=1)
     {
        C0 = C;
        C1 = C0+2*ldc;
        ptrba = ba;
        for (i=0; i<bm/2; i+=1)
          {
             ptrbb = bb;
             res0 = 0;
             res1 = 0;
             res2 = 0;
             res3 = 0;
             res4 = 0;
             res5 = 0;
             res6 = 0;
             res7 = 0;

#if defined(__wasm_simd128__) && !defined(TRMMKERNEL)
#ifndef DOUBLE
             {
                v128_t acc0 = wasm_f32x4_splat(0.0f);
                v128_t acc1 = wasm_f32x4_splat(0.0f);
                const float sign_even_mem[4] __attribute__((aligned(16))) = {-0.0f, 0.0f, -0.0f, 0.0f};
                const float sign_odd_mem[4] __attribute__((aligned(16))) = {0.0f, -0.0f, 0.0f, -0.0f};
                v128_t sign_even __attribute__((unused)) = wasm_v128_load(sign_even_mem);
                v128_t sign_odd __attribute__((unused)) = wasm_v128_load(sign_odd_mem);
                for (k = 0; k < bk; k += 1) {
                   v128_t va = wasm_v128_load(ptrba);
                   v128_t vb = wasm_v128_load(ptrbb);
                   v128_t va_m, vrot;
#if defined(NN) || defined(NT) || defined(TN) || defined(TT)
                   va_m = va;
                   vrot = wasm_v128_xor(wasm_i32x4_shuffle(va, va, 1, 0, 3, 2), sign_even);
#elif defined(NR) || defined(NC) || defined(TR) || defined(TC)
                   va_m = va;
                   vrot = wasm_v128_xor(wasm_i32x4_shuffle(va, va, 1, 0, 3, 2), sign_odd);
#elif defined(RN) || defined(RT) || defined(CN) || defined(CT)
                   va_m = wasm_v128_xor(va, sign_odd);
                   vrot = wasm_i32x4_shuffle(va, va, 1, 0, 3, 2);
#else
                   va_m = wasm_v128_xor(va, sign_odd);
                   vrot = wasm_f32x4_neg(wasm_i32x4_shuffle(va, va, 1, 0, 3, 2));
#endif
                   acc0 = MADD_F32(va_m, wasm_i32x4_shuffle(vb, vb, 0, 0, 0, 0), acc0);
                   acc0 = MADD_F32(vrot, wasm_i32x4_shuffle(vb, vb, 1, 1, 1, 1), acc0);
                   acc1 = MADD_F32(va_m, wasm_i32x4_shuffle(vb, vb, 2, 2, 2, 2), acc1);
                   acc1 = MADD_F32(vrot, wasm_i32x4_shuffle(vb, vb, 3, 3, 3, 3), acc1);
                   ptrba += 4;
                   ptrbb += 4;
                }
                res0 = wasm_f32x4_extract_lane(acc0, 0);
                res1 = wasm_f32x4_extract_lane(acc0, 1);
                res2 = wasm_f32x4_extract_lane(acc0, 2);
                res3 = wasm_f32x4_extract_lane(acc0, 3);
                res4 = wasm_f32x4_extract_lane(acc1, 0);
                res5 = wasm_f32x4_extract_lane(acc1, 1);
                res6 = wasm_f32x4_extract_lane(acc1, 2);
                res7 = wasm_f32x4_extract_lane(acc1, 3);
             }
#else
             {
                v128_t acc00 = wasm_f64x2_splat(0.0);
                v128_t acc10 = wasm_f64x2_splat(0.0);
                v128_t acc01 = wasm_f64x2_splat(0.0);
                v128_t acc11 = wasm_f64x2_splat(0.0);
                const double sign_re_mem[2] __attribute__((aligned(16))) = {-0.0, 0.0};
                const double sign_im_mem[2] __attribute__((aligned(16))) = {0.0, -0.0};
                v128_t sign_re __attribute__((unused)) = wasm_v128_load(sign_re_mem);
                v128_t sign_im __attribute__((unused)) = wasm_v128_load(sign_im_mem);
                for (k = 0; k < bk; k += 1) {
                   v128_t a0 = wasm_v128_load(ptrba);
                   v128_t a1 = wasm_v128_load(ptrba + 2);
                   v128_t b0 = wasm_v128_load(ptrbb);
                   v128_t b1 = wasm_v128_load(ptrbb + 2);
                   v128_t a0m, a1m, r0, r1;
#if defined(NN) || defined(NT) || defined(TN) || defined(TT)
                   a0m = a0; a1m = a1;
                   r0 = wasm_v128_xor(wasm_i64x2_shuffle(a0, a0, 1, 0), sign_re);
                   r1 = wasm_v128_xor(wasm_i64x2_shuffle(a1, a1, 1, 0), sign_re);
#elif defined(NR) || defined(NC) || defined(TR) || defined(TC)
                   a0m = a0; a1m = a1;
                   r0 = wasm_v128_xor(wasm_i64x2_shuffle(a0, a0, 1, 0), sign_im);
                   r1 = wasm_v128_xor(wasm_i64x2_shuffle(a1, a1, 1, 0), sign_im);
#elif defined(RN) || defined(RT) || defined(CN) || defined(CT)
                   a0m = wasm_v128_xor(a0, sign_im); a1m = wasm_v128_xor(a1, sign_im);
                   r0 = wasm_i64x2_shuffle(a0, a0, 1, 0);
                   r1 = wasm_i64x2_shuffle(a1, a1, 1, 0);
#else
                   a0m = wasm_v128_xor(a0, sign_im); a1m = wasm_v128_xor(a1, sign_im);
                   r0 = wasm_f64x2_neg(wasm_i64x2_shuffle(a0, a0, 1, 0));
                   r1 = wasm_f64x2_neg(wasm_i64x2_shuffle(a1, a1, 1, 0));
#endif
                   v128_t b0r = wasm_i64x2_shuffle(b0, b0, 0, 0);
                   v128_t b0i = wasm_i64x2_shuffle(b0, b0, 1, 1);
                   v128_t b1r = wasm_i64x2_shuffle(b1, b1, 0, 0);
                   v128_t b1i = wasm_i64x2_shuffle(b1, b1, 1, 1);
                   acc00 = MADD_F64(a0m, b0r, acc00);
                   acc00 = MADD_F64(r0, b0i, acc00);
                   acc10 = MADD_F64(a1m, b0r, acc10);
                   acc10 = MADD_F64(r1, b0i, acc10);
                   acc01 = MADD_F64(a0m, b1r, acc01);
                   acc01 = MADD_F64(r0, b1i, acc01);
                   acc11 = MADD_F64(a1m, b1r, acc11);
                   acc11 = MADD_F64(r1, b1i, acc11);
                   ptrba += 4;
                   ptrbb += 4;
                }
                res0 = wasm_f64x2_extract_lane(acc00, 0);
                res1 = wasm_f64x2_extract_lane(acc00, 1);
                res2 = wasm_f64x2_extract_lane(acc10, 0);
                res3 = wasm_f64x2_extract_lane(acc10, 1);
                res4 = wasm_f64x2_extract_lane(acc01, 0);
                res5 = wasm_f64x2_extract_lane(acc01, 1);
                res6 = wasm_f64x2_extract_lane(acc11, 0);
                res7 = wasm_f64x2_extract_lane(acc11, 1);
             }
#endif
#else
             for (k=0; k<bk/4; k+=1)
               {
#if   defined(NN) || defined(NT) || defined(TN) || defined(TT)
                  load0 = ptrba[4*0+0];
                  load1 = ptrbb[4*0+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[4*0+1];
                  res1 = res1+load2*load1;
                  load3 = ptrbb[4*0+1];
                  res0 = res0-load2*load3;
                  res1 = res1+load0*load3;
                  load4 = ptrba[4*0+2];
                  res2 = res2+load4*load1;
                  load5 = ptrba[4*0+3];
                  res3 = res3+load5*load1;
                  res2 = res2-load5*load3;
                  res3 = res3+load4*load3;
                  load6 = ptrbb[4*0+2];
                  res4 = res4+load0*load6;
                  res5 = res5+load2*load6;
                  load7 = ptrbb[4*0+3];
                  res4 = res4-load2*load7;
                  res5 = res5+load0*load7;
                  res6 = res6+load4*load6;
                  res7 = res7+load5*load6;
                  res6 = res6-load5*load7;
                  res7 = res7+load4*load7;
                  load8 = ptrba[4*1+0];
                  load9 = ptrbb[4*1+0];
                  res0 = res0+load8*load9;
                  load10 = ptrba[4*1+1];
                  res1 = res1+load10*load9;
                  load11 = ptrbb[4*1+1];
                  res0 = res0-load10*load11;
                  res1 = res1+load8*load11;
                  load12 = ptrba[4*1+2];
                  res2 = res2+load12*load9;
                  load13 = ptrba[4*1+3];
                  res3 = res3+load13*load9;
                  res2 = res2-load13*load11;
                  res3 = res3+load12*load11;
                  load14 = ptrbb[4*1+2];
                  res4 = res4+load8*load14;
                  res5 = res5+load10*load14;
                  load15 = ptrbb[4*1+3];
                  res4 = res4-load10*load15;
                  res5 = res5+load8*load15;
                  res6 = res6+load12*load14;
                  res7 = res7+load13*load14;
                  res6 = res6-load13*load15;
                  res7 = res7+load12*load15;
                  load0 = ptrba[4*2+0];
                  load1 = ptrbb[4*2+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[4*2+1];
                  res1 = res1+load2*load1;
                  load3 = ptrbb[4*2+1];
                  res0 = res0-load2*load3;
                  res1 = res1+load0*load3;
                  load4 = ptrba[4*2+2];
                  res2 = res2+load4*load1;
                  load5 = ptrba[4*2+3];
                  res3 = res3+load5*load1;
                  res2 = res2-load5*load3;
                  res3 = res3+load4*load3;
                  load6 = ptrbb[4*2+2];
                  res4 = res4+load0*load6;
                  res5 = res5+load2*load6;
                  load7 = ptrbb[4*2+3];
                  res4 = res4-load2*load7;
                  res5 = res5+load0*load7;
                  res6 = res6+load4*load6;
                  res7 = res7+load5*load6;
                  res6 = res6-load5*load7;
                  res7 = res7+load4*load7;
                  load8 = ptrba[4*3+0];
                  load9 = ptrbb[4*3+0];
                  res0 = res0+load8*load9;
                  load10 = ptrba[4*3+1];
                  res1 = res1+load10*load9;
                  load11 = ptrbb[4*3+1];
                  res0 = res0-load10*load11;
                  res1 = res1+load8*load11;
                  load12 = ptrba[4*3+2];
                  res2 = res2+load12*load9;
                  load13 = ptrba[4*3+3];
                  res3 = res3+load13*load9;
                  res2 = res2-load13*load11;
                  res3 = res3+load12*load11;
                  load14 = ptrbb[4*3+2];
                  res4 = res4+load8*load14;
                  res5 = res5+load10*load14;
                  load15 = ptrbb[4*3+3];
                  res4 = res4-load10*load15;
                  res5 = res5+load8*load15;
                  res6 = res6+load12*load14;
                  res7 = res7+load13*load14;
                  res6 = res6-load13*load15;
                  res7 = res7+load12*load15;
#endif
#if   defined(NR) || defined(NC) || defined(TR) || defined(TC)
                  load0 = ptrba[4*0+0];
                  load1 = ptrbb[4*0+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[4*0+1];
                  res1 = res1+load2*load1;
                  load3 = ptrbb[4*0+1];
                  res0 = res0+load2*load3;
                  res1 = res1-load0*load3;
                  load4 = ptrba[4*0+2];
                  res2 = res2+load4*load1;
                  load5 = ptrba[4*0+3];
                  res3 = res3+load5*load1;
                  res2 = res2+load5*load3;
                  res3 = res3-load4*load3;
                  load6 = ptrbb[4*0+2];
                  res4 = res4+load0*load6;
                  res5 = res5+load2*load6;
                  load7 = ptrbb[4*0+3];
                  res4 = res4+load2*load7;
                  res5 = res5-load0*load7;
                  res6 = res6+load4*load6;
                  res7 = res7+load5*load6;
                  res6 = res6+load5*load7;
                  res7 = res7-load4*load7;
                  load8 = ptrba[4*1+0];
                  load9 = ptrbb[4*1+0];
                  res0 = res0+load8*load9;
                  load10 = ptrba[4*1+1];
                  res1 = res1+load10*load9;
                  load11 = ptrbb[4*1+1];
                  res0 = res0+load10*load11;
                  res1 = res1-load8*load11;
                  load12 = ptrba[4*1+2];
                  res2 = res2+load12*load9;
                  load13 = ptrba[4*1+3];
                  res3 = res3+load13*load9;
                  res2 = res2+load13*load11;
                  res3 = res3-load12*load11;
                  load14 = ptrbb[4*1+2];
                  res4 = res4+load8*load14;
                  res5 = res5+load10*load14;
                  load15 = ptrbb[4*1+3];
                  res4 = res4+load10*load15;
                  res5 = res5-load8*load15;
                  res6 = res6+load12*load14;
                  res7 = res7+load13*load14;
                  res6 = res6+load13*load15;
                  res7 = res7-load12*load15;
                  load0 = ptrba[4*2+0];
                  load1 = ptrbb[4*2+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[4*2+1];
                  res1 = res1+load2*load1;
                  load3 = ptrbb[4*2+1];
                  res0 = res0+load2*load3;
                  res1 = res1-load0*load3;
                  load4 = ptrba[4*2+2];
                  res2 = res2+load4*load1;
                  load5 = ptrba[4*2+3];
                  res3 = res3+load5*load1;
                  res2 = res2+load5*load3;
                  res3 = res3-load4*load3;
                  load6 = ptrbb[4*2+2];
                  res4 = res4+load0*load6;
                  res5 = res5+load2*load6;
                  load7 = ptrbb[4*2+3];
                  res4 = res4+load2*load7;
                  res5 = res5-load0*load7;
                  res6 = res6+load4*load6;
                  res7 = res7+load5*load6;
                  res6 = res6+load5*load7;
                  res7 = res7-load4*load7;
                  load8 = ptrba[4*3+0];
                  load9 = ptrbb[4*3+0];
                  res0 = res0+load8*load9;
                  load10 = ptrba[4*3+1];
                  res1 = res1+load10*load9;
                  load11 = ptrbb[4*3+1];
                  res0 = res0+load10*load11;
                  res1 = res1-load8*load11;
                  load12 = ptrba[4*3+2];
                  res2 = res2+load12*load9;
                  load13 = ptrba[4*3+3];
                  res3 = res3+load13*load9;
                  res2 = res2+load13*load11;
                  res3 = res3-load12*load11;
                  load14 = ptrbb[4*3+2];
                  res4 = res4+load8*load14;
                  res5 = res5+load10*load14;
                  load15 = ptrbb[4*3+3];
                  res4 = res4+load10*load15;
                  res5 = res5-load8*load15;
                  res6 = res6+load12*load14;
                  res7 = res7+load13*load14;
                  res6 = res6+load13*load15;
                  res7 = res7-load12*load15;
#endif
#if   defined(RN) || defined(RT) || defined(CN) || defined(CT)
                  load0 = ptrba[4*0+0];
                  load1 = ptrbb[4*0+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[4*0+1];
                  res1 = res1-load2*load1;
                  load3 = ptrbb[4*0+1];
                  res0 = res0+load2*load3;
                  res1 = res1+load0*load3;
                  load4 = ptrba[4*0+2];
                  res2 = res2+load4*load1;
                  load5 = ptrba[4*0+3];
                  res3 = res3-load5*load1;
                  res2 = res2+load5*load3;
                  res3 = res3+load4*load3;
                  load6 = ptrbb[4*0+2];
                  res4 = res4+load0*load6;
                  res5 = res5-load2*load6;
                  load7 = ptrbb[4*0+3];
                  res4 = res4+load2*load7;
                  res5 = res5+load0*load7;
                  res6 = res6+load4*load6;
                  res7 = res7-load5*load6;
                  res6 = res6+load5*load7;
                  res7 = res7+load4*load7;
                  load8 = ptrba[4*1+0];
                  load9 = ptrbb[4*1+0];
                  res0 = res0+load8*load9;
                  load10 = ptrba[4*1+1];
                  res1 = res1-load10*load9;
                  load11 = ptrbb[4*1+1];
                  res0 = res0+load10*load11;
                  res1 = res1+load8*load11;
                  load12 = ptrba[4*1+2];
                  res2 = res2+load12*load9;
                  load13 = ptrba[4*1+3];
                  res3 = res3-load13*load9;
                  res2 = res2+load13*load11;
                  res3 = res3+load12*load11;
                  load14 = ptrbb[4*1+2];
                  res4 = res4+load8*load14;
                  res5 = res5-load10*load14;
                  load15 = ptrbb[4*1+3];
                  res4 = res4+load10*load15;
                  res5 = res5+load8*load15;
                  res6 = res6+load12*load14;
                  res7 = res7-load13*load14;
                  res6 = res6+load13*load15;
                  res7 = res7+load12*load15;
                  load0 = ptrba[4*2+0];
                  load1 = ptrbb[4*2+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[4*2+1];
                  res1 = res1-load2*load1;
                  load3 = ptrbb[4*2+1];
                  res0 = res0+load2*load3;
                  res1 = res1+load0*load3;
                  load4 = ptrba[4*2+2];
                  res2 = res2+load4*load1;
                  load5 = ptrba[4*2+3];
                  res3 = res3-load5*load1;
                  res2 = res2+load5*load3;
                  res3 = res3+load4*load3;
                  load6 = ptrbb[4*2+2];
                  res4 = res4+load0*load6;
                  res5 = res5-load2*load6;
                  load7 = ptrbb[4*2+3];
                  res4 = res4+load2*load7;
                  res5 = res5+load0*load7;
                  res6 = res6+load4*load6;
                  res7 = res7-load5*load6;
                  res6 = res6+load5*load7;
                  res7 = res7+load4*load7;
                  load8 = ptrba[4*3+0];
                  load9 = ptrbb[4*3+0];
                  res0 = res0+load8*load9;
                  load10 = ptrba[4*3+1];
                  res1 = res1-load10*load9;
                  load11 = ptrbb[4*3+1];
                  res0 = res0+load10*load11;
                  res1 = res1+load8*load11;
                  load12 = ptrba[4*3+2];
                  res2 = res2+load12*load9;
                  load13 = ptrba[4*3+3];
                  res3 = res3-load13*load9;
                  res2 = res2+load13*load11;
                  res3 = res3+load12*load11;
                  load14 = ptrbb[4*3+2];
                  res4 = res4+load8*load14;
                  res5 = res5-load10*load14;
                  load15 = ptrbb[4*3+3];
                  res4 = res4+load10*load15;
                  res5 = res5+load8*load15;
                  res6 = res6+load12*load14;
                  res7 = res7-load13*load14;
                  res6 = res6+load13*load15;
                  res7 = res7+load12*load15;
#endif
#if   defined(RR) || defined(RC) || defined(CR) || defined(CC)
                  load0 = ptrba[4*0+0];
                  load1 = ptrbb[4*0+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[4*0+1];
                  res1 = res1-load2*load1;
                  load3 = ptrbb[4*0+1];
                  res0 = res0-load2*load3;
                  res1 = res1-load0*load3;
                  load4 = ptrba[4*0+2];
                  res2 = res2+load4*load1;
                  load5 = ptrba[4*0+3];
                  res3 = res3-load5*load1;
                  res2 = res2-load5*load3;
                  res3 = res3-load4*load3;
                  load6 = ptrbb[4*0+2];
                  res4 = res4+load0*load6;
                  res5 = res5-load2*load6;
                  load7 = ptrbb[4*0+3];
                  res4 = res4-load2*load7;
                  res5 = res5-load0*load7;
                  res6 = res6+load4*load6;
                  res7 = res7-load5*load6;
                  res6 = res6-load5*load7;
                  res7 = res7-load4*load7;
                  load8 = ptrba[4*1+0];
                  load9 = ptrbb[4*1+0];
                  res0 = res0+load8*load9;
                  load10 = ptrba[4*1+1];
                  res1 = res1-load10*load9;
                  load11 = ptrbb[4*1+1];
                  res0 = res0-load10*load11;
                  res1 = res1-load8*load11;
                  load12 = ptrba[4*1+2];
                  res2 = res2+load12*load9;
                  load13 = ptrba[4*1+3];
                  res3 = res3-load13*load9;
                  res2 = res2-load13*load11;
                  res3 = res3-load12*load11;
                  load14 = ptrbb[4*1+2];
                  res4 = res4+load8*load14;
                  res5 = res5-load10*load14;
                  load15 = ptrbb[4*1+3];
                  res4 = res4-load10*load15;
                  res5 = res5-load8*load15;
                  res6 = res6+load12*load14;
                  res7 = res7-load13*load14;
                  res6 = res6-load13*load15;
                  res7 = res7-load12*load15;
                  load0 = ptrba[4*2+0];
                  load1 = ptrbb[4*2+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[4*2+1];
                  res1 = res1-load2*load1;
                  load3 = ptrbb[4*2+1];
                  res0 = res0-load2*load3;
                  res1 = res1-load0*load3;
                  load4 = ptrba[4*2+2];
                  res2 = res2+load4*load1;
                  load5 = ptrba[4*2+3];
                  res3 = res3-load5*load1;
                  res2 = res2-load5*load3;
                  res3 = res3-load4*load3;
                  load6 = ptrbb[4*2+2];
                  res4 = res4+load0*load6;
                  res5 = res5-load2*load6;
                  load7 = ptrbb[4*2+3];
                  res4 = res4-load2*load7;
                  res5 = res5-load0*load7;
                  res6 = res6+load4*load6;
                  res7 = res7-load5*load6;
                  res6 = res6-load5*load7;
                  res7 = res7-load4*load7;
                  load8 = ptrba[4*3+0];
                  load9 = ptrbb[4*3+0];
                  res0 = res0+load8*load9;
                  load10 = ptrba[4*3+1];
                  res1 = res1-load10*load9;
                  load11 = ptrbb[4*3+1];
                  res0 = res0-load10*load11;
                  res1 = res1-load8*load11;
                  load12 = ptrba[4*3+2];
                  res2 = res2+load12*load9;
                  load13 = ptrba[4*3+3];
                  res3 = res3-load13*load9;
                  res2 = res2-load13*load11;
                  res3 = res3-load12*load11;
                  load14 = ptrbb[4*3+2];
                  res4 = res4+load8*load14;
                  res5 = res5-load10*load14;
                  load15 = ptrbb[4*3+3];
                  res4 = res4-load10*load15;
                  res5 = res5-load8*load15;
                  res6 = res6+load12*load14;
                  res7 = res7-load13*load14;
                  res6 = res6-load13*load15;
                  res7 = res7-load12*load15;
#endif
                  ptrba = ptrba+16;
                  ptrbb = ptrbb+16;
               }
             for (k=0; k<(bk&3); k+=1)
               {
#if   defined(NN) || defined(NT) || defined(TN) || defined(TT)
                  load0 = ptrba[4*0+0];
                  load1 = ptrbb[4*0+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[4*0+1];
                  res1 = res1+load2*load1;
                  load3 = ptrbb[4*0+1];
                  res0 = res0-load2*load3;
                  res1 = res1+load0*load3;
                  load4 = ptrba[4*0+2];
                  res2 = res2+load4*load1;
                  load5 = ptrba[4*0+3];
                  res3 = res3+load5*load1;
                  res2 = res2-load5*load3;
                  res3 = res3+load4*load3;
                  load6 = ptrbb[4*0+2];
                  res4 = res4+load0*load6;
                  res5 = res5+load2*load6;
                  load7 = ptrbb[4*0+3];
                  res4 = res4-load2*load7;
                  res5 = res5+load0*load7;
                  res6 = res6+load4*load6;
                  res7 = res7+load5*load6;
                  res6 = res6-load5*load7;
                  res7 = res7+load4*load7;
#endif
#if   defined(NR) || defined(NC) || defined(TR) || defined(TC)
                  load0 = ptrba[4*0+0];
                  load1 = ptrbb[4*0+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[4*0+1];
                  res1 = res1+load2*load1;
                  load3 = ptrbb[4*0+1];
                  res0 = res0+load2*load3;
                  res1 = res1-load0*load3;
                  load4 = ptrba[4*0+2];
                  res2 = res2+load4*load1;
                  load5 = ptrba[4*0+3];
                  res3 = res3+load5*load1;
                  res2 = res2+load5*load3;
                  res3 = res3-load4*load3;
                  load6 = ptrbb[4*0+2];
                  res4 = res4+load0*load6;
                  res5 = res5+load2*load6;
                  load7 = ptrbb[4*0+3];
                  res4 = res4+load2*load7;
                  res5 = res5-load0*load7;
                  res6 = res6+load4*load6;
                  res7 = res7+load5*load6;
                  res6 = res6+load5*load7;
                  res7 = res7-load4*load7;
#endif
#if   defined(RN) || defined(RT) || defined(CN) || defined(CT)
                  load0 = ptrba[4*0+0];
                  load1 = ptrbb[4*0+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[4*0+1];
                  res1 = res1-load2*load1;
                  load3 = ptrbb[4*0+1];
                  res0 = res0+load2*load3;
                  res1 = res1+load0*load3;
                  load4 = ptrba[4*0+2];
                  res2 = res2+load4*load1;
                  load5 = ptrba[4*0+3];
                  res3 = res3-load5*load1;
                  res2 = res2+load5*load3;
                  res3 = res3+load4*load3;
                  load6 = ptrbb[4*0+2];
                  res4 = res4+load0*load6;
                  res5 = res5-load2*load6;
                  load7 = ptrbb[4*0+3];
                  res4 = res4+load2*load7;
                  res5 = res5+load0*load7;
                  res6 = res6+load4*load6;
                  res7 = res7-load5*load6;
                  res6 = res6+load5*load7;
                  res7 = res7+load4*load7;
#endif
#if   defined(RR) || defined(RC) || defined(CR) || defined(CC)
                  load0 = ptrba[4*0+0];
                  load1 = ptrbb[4*0+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[4*0+1];
                  res1 = res1-load2*load1;
                  load3 = ptrbb[4*0+1];
                  res0 = res0-load2*load3;
                  res1 = res1-load0*load3;
                  load4 = ptrba[4*0+2];
                  res2 = res2+load4*load1;
                  load5 = ptrba[4*0+3];
                  res3 = res3-load5*load1;
                  res2 = res2-load5*load3;
                  res3 = res3-load4*load3;
                  load6 = ptrbb[4*0+2];
                  res4 = res4+load0*load6;
                  res5 = res5-load2*load6;
                  load7 = ptrbb[4*0+3];
                  res4 = res4-load2*load7;
                  res5 = res5-load0*load7;
                  res6 = res6+load4*load6;
                  res7 = res7-load5*load6;
                  res6 = res6-load5*load7;
                  res7 = res7-load4*load7;
#endif
                  ptrba = ptrba+4;
                  ptrbb = ptrbb+4;
               }
#endif
             load0 = res0*alphar;
             C0[0] = C0[0]+load0;
             load1 = res1*alphar;
             C0[1] = C0[1]+load1;
             load0 = res1*alphai;
             C0[0] = C0[0]-load0;
             load1 = res0*alphai;
             C0[1] = C0[1]+load1;
             load2 = res2*alphar;
             C0[2] = C0[2]+load2;
             load3 = res3*alphar;
             C0[3] = C0[3]+load3;
             load2 = res3*alphai;
             C0[2] = C0[2]-load2;
             load3 = res2*alphai;
             C0[3] = C0[3]+load3;
             load4 = res4*alphar;
             C1[0] = C1[0]+load4;
             load5 = res5*alphar;
             C1[1] = C1[1]+load5;
             load4 = res5*alphai;
             C1[0] = C1[0]-load4;
             load5 = res4*alphai;
             C1[1] = C1[1]+load5;
             load6 = res6*alphar;
             C1[2] = C1[2]+load6;
             load7 = res7*alphar;
             C1[3] = C1[3]+load7;
             load6 = res7*alphai;
             C1[2] = C1[2]-load6;
             load7 = res6*alphai;
             C1[3] = C1[3]+load7;
             C0 = C0+4;
             C1 = C1+4;
          }
        for (i=0; i<(bm&1); i+=1)
          {
             ptrbb = bb;
             res0 = 0;
             res1 = 0;
             res2 = 0;
             res3 = 0;
             for (k=0; k<bk; k+=1)
               {
#if   defined(NN) || defined(NT) || defined(TN) || defined(TT)
                  load0 = ptrba[2*0+0];
                  load1 = ptrbb[4*0+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[2*0+1];
                  res1 = res1+load2*load1;
                  load3 = ptrbb[4*0+1];
                  res0 = res0-load2*load3;
                  res1 = res1+load0*load3;
                  load4 = ptrbb[4*0+2];
                  res2 = res2+load0*load4;
                  res3 = res3+load2*load4;
                  load5 = ptrbb[4*0+3];
                  res2 = res2-load2*load5;
                  res3 = res3+load0*load5;
#endif
#if   defined(NR) || defined(NC) || defined(TR) || defined(TC)
				  load0 = ptrba[2*0+0];
                  load1 = ptrbb[4*0+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[2*0+1];
                  res1 = res1+load2*load1;
                  load3 = ptrbb[4*0+1];
                  res0 = res0+load2*load3;
                  res1 = res1-load0*load3;
                  load4 = ptrbb[4*0+2];
                  res2 = res2+load0*load4;
                  res3 = res3+load2*load4;
                  load5 = ptrbb[4*0+3];
                  res2 = res2+load2*load5;
                  res3 = res3-load0*load5;
#endif
#if   defined(RN) || defined(RT) || defined(CN) || defined(CT)
                  load0 = ptrba[2*0+0];
                  load1 = ptrbb[4*0+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[2*0+1];
                  res1 = res1-load2*load1;
                  load3 = ptrbb[4*0+1];
                  res0 = res0+load2*load3;
                  res1 = res1+load0*load3;
                  load4 = ptrbb[4*0+2];
                  res2 = res2+load0*load4;
                  res3 = res3-load2*load4;
                  load5 = ptrbb[4*0+3];
                  res2 = res2+load2*load5;
                  res3 = res3+load0*load5;
#endif
#if   defined(RR) || defined(RC) || defined(CR) || defined(CC)
                  load0 = ptrba[2*0+0];
                  load1 = ptrbb[4*0+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[2*0+1];
                  res1 = res1-load2*load1;
                  load3 = ptrbb[4*0+1];
                  res0 = res0-load2*load3;
                  res1 = res1-load0*load3;
                  load4 = ptrbb[4*0+2];
                  res2 = res2+load0*load4;
                  res3 = res3-load2*load4;
                  load5 = ptrbb[4*0+3];
                  res2 = res2-load2*load5;
                  res3 = res3-load0*load5;
#endif
                  ptrba = ptrba+2;
                  ptrbb = ptrbb+4;
               }
             load0 = res0*alphar;
             C0[0] = C0[0]+load0;
             load1 = res1*alphar;
             C0[1] = C0[1]+load1;
             load0 = res1*alphai;
             C0[0] = C0[0]-load0;
             load1 = res0*alphai;
             C0[1] = C0[1]+load1;
             load2 = res2*alphar;
             C1[0] = C1[0]+load2;
             load3 = res3*alphar;
             C1[1] = C1[1]+load3;
             load2 = res3*alphai;
             C1[0] = C1[0]-load2;
             load3 = res2*alphai;
             C1[1] = C1[1]+load3;
             C0 = C0+2;
             C1 = C1+2;
          }
        k = (bk<<2);
        bb = bb+k;
        i = (ldc<<2);
        C = C+i;
     }
   for (j=0; j<(bn&1); j+=1)
     {
        C0 = C;
        ptrba = ba;
        for (i=0; i<bm/2; i+=1)
          {
             ptrbb = bb;
             res0 = 0;
             res1 = 0;
             res2 = 0;
             res3 = 0;
             for (k=0; k<bk; k+=1)
               {
#if   defined(NN) || defined(NT) || defined(TN) || defined(TT)
                  load0 = ptrba[4*0+0];
                  load1 = ptrbb[2*0+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[4*0+1];
                  res1 = res1+load2*load1;
                  load3 = ptrbb[2*0+1];
                  res0 = res0-load2*load3;
                  res1 = res1+load0*load3;
                  load4 = ptrba[4*0+2];
                  res2 = res2+load4*load1;
                  load5 = ptrba[4*0+3];
                  res3 = res3+load5*load1;
                  res2 = res2-load5*load3;
                  res3 = res3+load4*load3;
#endif
#if   defined(NR) || defined(NC) || defined(TR) || defined(TC)
                  load0 = ptrba[4*0+0];
                  load1 = ptrbb[2*0+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[4*0+1];
                  res1 = res1+load2*load1;
                  load3 = ptrbb[2*0+1];
                  res0 = res0+load2*load3;
                  res1 = res1-load0*load3;
                  load4 = ptrba[4*0+2];
                  res2 = res2+load4*load1;
                  load5 = ptrba[4*0+3];
                  res3 = res3+load5*load1;
                  res2 = res2+load5*load3;
                  res3 = res3-load4*load3;
#endif
#if   defined(RN) || defined(RT) || defined(CN) || defined(CT)
                  load0 = ptrba[4*0+0];
                  load1 = ptrbb[2*0+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[4*0+1];
                  res1 = res1-load2*load1;
                  load3 = ptrbb[2*0+1];
                  res0 = res0+load2*load3;
                  res1 = res1+load0*load3;
                  load4 = ptrba[4*0+2];
                  res2 = res2+load4*load1;
                  load5 = ptrba[4*0+3];
                  res3 = res3-load5*load1;
                  res2 = res2+load5*load3;
                  res3 = res3+load4*load3;
#endif
#if   defined(RR) || defined(RC) || defined(CR) || defined(CC)
                  load0 = ptrba[4*0+0];
                  load1 = ptrbb[2*0+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[4*0+1];
                  res1 = res1-load2*load1;
                  load3 = ptrbb[2*0+1];
                  res0 = res0-load2*load3;
                  res1 = res1-load0*load3;
                  load4 = ptrba[4*0+2];
                  res2 = res2+load4*load1;
                  load5 = ptrba[4*0+3];
                  res3 = res3-load5*load1;
                  res2 = res2-load5*load3;
                  res3 = res3-load4*load3;
#endif
                  ptrba = ptrba+4;
                  ptrbb = ptrbb+2;
               }
             load0 = res0*alphar;
             C0[0] = C0[0]+load0;
             load1 = res1*alphar;
             C0[1] = C0[1]+load1;
             load0 = res1*alphai;
             C0[0] = C0[0]-load0;
             load1 = res0*alphai;
             C0[1] = C0[1]+load1;
             load2 = res2*alphar;
             C0[2] = C0[2]+load2;
             load3 = res3*alphar;
             C0[3] = C0[3]+load3;
             load2 = res3*alphai;
             C0[2] = C0[2]-load2;
             load3 = res2*alphai;
             C0[3] = C0[3]+load3;
             C0 = C0+4;
          }
        for (i=0; i<(bm&1); i+=1)
          {
             ptrbb = bb;
             res0 = 0;
             res1 = 0;
             for (k=0; k<bk; k+=1)
               {
#if   defined(NN) || defined(NT) || defined(TN) || defined(TT)
                  load0 = ptrba[2*0+0];
                  load1 = ptrbb[2*0+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[2*0+1];
                  res1 = res1+load2*load1;
                  load3 = ptrbb[2*0+1];
                  res0 = res0-load2*load3;
                  res1 = res1+load0*load3;
#endif
#if   defined(NR) || defined(NC) || defined(TR) || defined(TC)
                  load0 = ptrba[2*0+0];
                  load1 = ptrbb[2*0+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[2*0+1];
                  res1 = res1+load2*load1;
                  load3 = ptrbb[2*0+1];
                  res0 = res0+load2*load3;
                  res1 = res1-load0*load3;
#endif
#if   defined(RN) || defined(RT) || defined(CN) || defined(CT)
                  load0 = ptrba[2*0+0];
                  load1 = ptrbb[2*0+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[2*0+1];
                  res1 = res1-load2*load1;
                  load3 = ptrbb[2*0+1];
                  res0 = res0+load2*load3;
                  res1 = res1+load0*load3;
#endif
#if   defined(RR) || defined(RC) || defined(CR) || defined(CC)
                  load0 = ptrba[2*0+0];
                  load1 = ptrbb[2*0+0];
                  res0 = res0+load0*load1;
                  load2 = ptrba[2*0+1];
                  res1 = res1-load2*load1;
                  load3 = ptrbb[2*0+1];
                  res0 = res0-load2*load3;
                  res1 = res1-load0*load3;
#endif
                  ptrba = ptrba+2;
                  ptrbb = ptrbb+2;
               }
             load0 = res0*alphar;
             C0[0] = C0[0]+load0;
             load1 = res1*alphar;
             C0[1] = C0[1]+load1;
             load0 = res1*alphai;
             C0[0] = C0[0]-load0;
             load1 = res0*alphai;
             C0[1] = C0[1]+load1;
             C0 = C0+2;
          }
        k = (bk<<1);
        bb = bb+k;
        i = (ldc<<1);
        C = C+i;
     }
   return 0;
}
