/***************************************************************************
Copyright (c) 2026, The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

/* SME2 types, ZA tile access and C block transfers for sme2_gemm_impl.h and the packed level-3 kernels. */

#ifndef SME2_GEMM_TILE_H
#define SME2_GEMM_TILE_H

#include <arm_sme.h>
#include <stdint.h>

/* The file is built with the SME flags of its target (no SME2); the functions here enable SME2 themselves. */
#pragma clang attribute push(__attribute__((target("sme2"))), apply_to = function)

#define S2_S __arm_streaming
#define S2_ZA __arm_inout("za")
#define S2_INL static inline __attribute__((always_inline))
#define S2_NOINL static __attribute__((noinline))
#define S2_UNROLL _Pragma("GCC unroll 8")

#ifdef DOUBLE
#define S2_VL 8
#define S2_NT 8
typedef svfloat64_t s2_v;
typedef svfloat64x4_t s2_v4;
#define S2_PT() svptrue_b64()
#define S2_PW(i, n) svwhilelt_b64_s64((int64_t)(i), (int64_t)(n))
#define S2_CT() svptrue_c64()
#define S2_CW(i, n) svwhilelt_c64_s64((int64_t)(i), (int64_t)(n), 4)
#define S2_LD1(p, q) svld1_f64(p, q)
#define S2_ST1(p, q, v) svst1_f64(p, q, v)
#define S2_LD4(c, q) svld1_f64_x4(c, q)
#define S2_ST4(c, q, v) svst1_f64_x4(c, q, v)
#define S2_GET4(v, i) svget4_f64(v, i)
#define S2_CREATE4(a, b, c, d) svcreate4_f64(a, b, c, d)
#define S2_ZEROV() svdup_n_f64(0)
#define S2_MUL(v, a) svmul_n_f64_x(svptrue_b64(), v, a)
#define S2_MOPA(t, p, a, b) svmopa_za64_f64_m(t, p, p, a, b)
#define S2_RDH(t, s) svread_hor_za64_f64_m(svundef_f64(), svptrue_b64(), t, s)
#define S2_WRH(t, s, v) svwrite_hor_za64_f64_m(t, s, svptrue_b64(), v)
#define S2_RDV4(t, s) svread_ver_za64_f64_vg4(t, s)
#define S2_RDH4(t, s) svread_hor_za64_f64_vg4(t, s)
#define S2_WRH4(t, s, v) svwrite_hor_za64_f64_vg4(t, s, v)
#define S2_MOPS(t, pn, pm, a, b) svmops_za64_f64_m(t, pn, pm, a, b)
#define S2_RDV(t, s) svread_ver_za64_f64_m(svundef_f64(), svptrue_b64(), t, s)
#define S2_WRV(t, s, v) svwrite_ver_za64_f64_m(t, s, svptrue_b64(), v)
#define S2_MLA(c, v, a) svmla_n_f64_x(svptrue_b64(), c, v, a)
#else
#define S2_VL 16
#define S2_NT 4
typedef svfloat32_t s2_v;
typedef svfloat32x4_t s2_v4;
#define S2_PT() svptrue_b32()
#define S2_PW(i, n) svwhilelt_b32_s64((int64_t)(i), (int64_t)(n))
#define S2_CT() svptrue_c32()
#define S2_CW(i, n) svwhilelt_c32_s64((int64_t)(i), (int64_t)(n), 4)
#define S2_LD1(p, q) svld1_f32(p, q)
#define S2_ST1(p, q, v) svst1_f32(p, q, v)
#define S2_LD4(c, q) svld1_f32_x4(c, q)
#define S2_ST4(c, q, v) svst1_f32_x4(c, q, v)
#define S2_GET4(v, i) svget4_f32(v, i)
#define S2_CREATE4(a, b, c, d) svcreate4_f32(a, b, c, d)
#define S2_ZEROV() svdup_n_f32(0)
#define S2_MUL(v, a) svmul_n_f32_x(svptrue_b32(), v, a)
#define S2_MOPA(t, p, a, b) svmopa_za32_f32_m(t, p, p, a, b)
#define S2_RDH(t, s) svread_hor_za32_f32_m(svundef_f32(), svptrue_b32(), t, s)
#define S2_WRH(t, s, v) svwrite_hor_za32_f32_m(t, s, svptrue_b32(), v)
#define S2_RDV4(t, s) svread_ver_za32_f32_vg4(t, s)
#define S2_RDH4(t, s) svread_hor_za32_f32_vg4(t, s)
#define S2_WRH4(t, s, v) svwrite_hor_za32_f32_vg4(t, s, v)
#define S2_MOPS(t, pn, pm, a, b) svmops_za32_f32_m(t, pn, pm, a, b)
#define S2_RDV(t, s) svread_ver_za32_f32_m(svundef_f32(), svptrue_b32(), t, s)
#define S2_WRV(t, s, v) svwrite_ver_za32_f32_m(t, s, svptrue_b32(), v)
#define S2_MLA(c, v, a) svmla_n_f32_x(svptrue_b32(), c, v, a)
#endif

#define S2_NR (S2_NT * S2_VL) /* columns of one row of NT tiles */

enum { S2_ZERO = 0, S2_LOAD = 1, S2_SCALE = 2 };

/* Tile numbers of the ZA intrinsics must be literals; with a constant t these switches fold away when inlined. */
S2_INL void s2_mopa(int t, svbool_t p, s2_v a, s2_v b) S2_S S2_ZA {
  switch (t) {
    case 0: S2_MOPA(0, p, a, b); break;
    case 1: S2_MOPA(1, p, a, b); break;
    case 2: S2_MOPA(2, p, a, b); break;
    case 3: S2_MOPA(3, p, a, b); break;
#ifdef DOUBLE
    case 4: S2_MOPA(4, p, a, b); break;
    case 5: S2_MOPA(5, p, a, b); break;
    case 6: S2_MOPA(6, p, a, b); break;
    case 7: S2_MOPA(7, p, a, b); break;
#endif
  }
}
S2_INL s2_v s2_rdh(int t, uint32_t s) S2_S S2_ZA {
  switch (t) {
    case 0: return S2_RDH(0, s);
    case 1: return S2_RDH(1, s);
    case 2: return S2_RDH(2, s);
#ifdef DOUBLE
    case 3: return S2_RDH(3, s);
    case 4: return S2_RDH(4, s);
    case 5: return S2_RDH(5, s);
    case 6: return S2_RDH(6, s);
    default: return S2_RDH(7, s);
#else
    default: return S2_RDH(3, s);
#endif
  }
}
S2_INL void s2_wrh(int t, uint32_t s, s2_v v) S2_S S2_ZA {
  switch (t) {
    case 0: S2_WRH(0, s, v); break;
    case 1: S2_WRH(1, s, v); break;
    case 2: S2_WRH(2, s, v); break;
    case 3: S2_WRH(3, s, v); break;
#ifdef DOUBLE
    case 4: S2_WRH(4, s, v); break;
    case 5: S2_WRH(5, s, v); break;
    case 6: S2_WRH(6, s, v); break;
    case 7: S2_WRH(7, s, v); break;
#endif
  }
}
/* Outer product subtracted from tile t, rows under pn and columns under pm. */
S2_INL void s2_mops(int t, svbool_t pn, svbool_t pm, s2_v a, s2_v b) S2_S S2_ZA {
  switch (t) {
    case 0: S2_MOPS(0, pn, pm, a, b); break;
    case 1: S2_MOPS(1, pn, pm, a, b); break;
    case 2: S2_MOPS(2, pn, pm, a, b); break;
    case 3: S2_MOPS(3, pn, pm, a, b); break;
#ifdef DOUBLE
    case 4: S2_MOPS(4, pn, pm, a, b); break;
    case 5: S2_MOPS(5, pn, pm, a, b); break;
    case 6: S2_MOPS(6, pn, pm, a, b); break;
    case 7: S2_MOPS(7, pn, pm, a, b); break;
#endif
  }
}
S2_INL s2_v s2_rdv(int t, uint32_t s) S2_S S2_ZA {
  switch (t) {
    case 0: return S2_RDV(0, s);
    case 1: return S2_RDV(1, s);
    case 2: return S2_RDV(2, s);
#ifdef DOUBLE
    case 3: return S2_RDV(3, s);
    case 4: return S2_RDV(4, s);
    case 5: return S2_RDV(5, s);
    case 6: return S2_RDV(6, s);
    default: return S2_RDV(7, s);
#else
    default: return S2_RDV(3, s);
#endif
  }
}
S2_INL void s2_wrv(int t, uint32_t s, s2_v v) S2_S S2_ZA {
  switch (t) {
    case 0: S2_WRV(0, s, v); break;
    case 1: S2_WRV(1, s, v); break;
    case 2: S2_WRV(2, s, v); break;
    case 3: S2_WRV(3, s, v); break;
#ifdef DOUBLE
    case 4: S2_WRV(4, s, v); break;
    case 5: S2_WRV(5, s, v); break;
    case 6: S2_WRV(6, s, v); break;
    case 7: S2_WRV(7, s, v); break;
#endif
  }
}
S2_INL s2_v4 s2_rdv4(int t, uint32_t s) S2_S S2_ZA {
  switch (t) {
    case 0: return S2_RDV4(0, s);
    case 1: return S2_RDV4(1, s);
    case 2: return S2_RDV4(2, s);
#ifdef DOUBLE
    case 3: return S2_RDV4(3, s);
    case 4: return S2_RDV4(4, s);
    case 5: return S2_RDV4(5, s);
    case 6: return S2_RDV4(6, s);
    default: return S2_RDV4(7, s);
#else
    default: return S2_RDV4(3, s);
#endif
  }
}
/* Four consecutive horizontal slices of tile t at once (MOVA vg4). */
S2_INL s2_v4 s2_rdh4(int t, uint32_t s) S2_S S2_ZA {
  switch (t) {
    case 0: return S2_RDH4(0, s);
    case 1: return S2_RDH4(1, s);
    case 2: return S2_RDH4(2, s);
#ifdef DOUBLE
    case 3: return S2_RDH4(3, s);
    case 4: return S2_RDH4(4, s);
    case 5: return S2_RDH4(5, s);
    case 6: return S2_RDH4(6, s);
    default: return S2_RDH4(7, s);
#else
    default: return S2_RDH4(3, s);
#endif
  }
}
S2_INL void s2_wrh4(int t, uint32_t s, s2_v4 v) S2_S S2_ZA {
  switch (t) {
    case 0: S2_WRH4(0, s, v); break;
    case 1: S2_WRH4(1, s, v); break;
    case 2: S2_WRH4(2, s, v); break;
    case 3: S2_WRH4(3, s, v); break;
#ifdef DOUBLE
    case 4: S2_WRH4(4, s, v); break;
    case 5: S2_WRH4(5, s, v); break;
    case 6: S2_WRH4(6, s, v); break;
    case 7: S2_WRH4(7, s, v); break;
#endif
  }
}
S2_INL s2_v s2_get4(s2_v4 v, int i) S2_S {
  switch (i) {
    case 0: return S2_GET4(v, 0);
    case 1: return S2_GET4(v, 1);
    case 2: return S2_GET4(v, 2);
    default: return S2_GET4(v, 3);
  }
}
S2_INL s2_v4 s2_scale4(s2_v4 v, FLOAT a) S2_S {
  return S2_CREATE4(S2_MUL(S2_GET4(v, 0), a), S2_MUL(S2_GET4(v, 1), a), S2_MUL(S2_GET4(v, 2), a), S2_MUL(S2_GET4(v, 3), a));
}
/* Vector i of the pair (x0, x1) of four-vector groups. */
S2_INL s2_v s2_pick(s2_v4 x0, s2_v4 x1, int i) S2_S { return i < 4 ? s2_get4(x0, i) : s2_get4(x1, i - 4); }

S2_INL int64_t s2_clip(int64_t n, int64_t hi) __arm_streaming_compatible { return n < 0 ? 0 : (n > hi ? hi : n); }
S2_INL s2_v4 s2_ld4(const FLOAT *p) S2_S { return S2_LD4(S2_CT(), p); }
S2_INL s2_v4 s2_ld4n(const FLOAT *p, int64_t n) S2_S { return S2_LD4(S2_CW(0, n), p); }
S2_INL void s2_st4(FLOAT *p, s2_v4 v) S2_S { S2_ST4(S2_CT(), p, v); }
S2_INL void s2_st4n(FLOAT *p, s2_v4 v, int64_t n) S2_S { S2_ST4(S2_CW(0, n), p, v); }

/* Core-side prefetch into L2 (prfm pldl2keep); the SME unit reads through L2. */
S2_INL void s2_pf_l2(const void *p, int bytes) __arm_streaming_compatible {
  for (int l = 0; l < bytes; l += 128) __builtin_prefetch((const char *)p + l, 0, 2);
}

/* Four rows x 64 columns into horizontal slices r..r+3 of every tile: strided-register x4 loads put the four rows
   of one tile into z(4t)..z(4t+3), so one MOVA vg4 per tile suffices. */
S2_INL void s2_rows_in4(const FLOAT *src, BLASLONG ld, uint32_t r) S2_S S2_ZA {
  const FLOAT *p1 = src + ld, *p2 = src + 2 * ld, *p3 = src + 3 * ld;
#ifdef DOUBLE
  __asm__ volatile(
      "ptrue pn8.d\n"
      "ld1d {z16.d, z20.d, z24.d, z28.d}, pn8/z, [%[a0]]\n"
      "ld1d {z17.d, z21.d, z25.d, z29.d}, pn8/z, [%[a1]]\n"
      "ld1d {z18.d, z22.d, z26.d, z30.d}, pn8/z, [%[a2]]\n"
      "ld1d {z19.d, z23.d, z27.d, z31.d}, pn8/z, [%[a3]]\n"
      "mova za0h.d[%w[r], 0:3], {z16.d - z19.d}\n"
      "mova za1h.d[%w[r], 0:3], {z20.d - z23.d}\n"
      "mova za2h.d[%w[r], 0:3], {z24.d - z27.d}\n"
      "mova za3h.d[%w[r], 0:3], {z28.d - z31.d}\n"
      "ld1d {z16.d, z20.d, z24.d, z28.d}, pn8/z, [%[a0], #4, mul vl]\n"
      "ld1d {z17.d, z21.d, z25.d, z29.d}, pn8/z, [%[a1], #4, mul vl]\n"
      "ld1d {z18.d, z22.d, z26.d, z30.d}, pn8/z, [%[a2], #4, mul vl]\n"
      "ld1d {z19.d, z23.d, z27.d, z31.d}, pn8/z, [%[a3], #4, mul vl]\n"
      "mova za4h.d[%w[r], 0:3], {z16.d - z19.d}\n"
      "mova za5h.d[%w[r], 0:3], {z20.d - z23.d}\n"
      "mova za6h.d[%w[r], 0:3], {z24.d - z27.d}\n"
      "mova za7h.d[%w[r], 0:3], {z28.d - z31.d}\n"
      :
      : [a0] "r"(src), [a1] "r"(p1), [a2] "r"(p2), [a3] "r"(p3), [r] "Ucj"(r)
      : "p8", "z16", "z17", "z18", "z19", "z20", "z21", "z22", "z23", "z24", "z25", "z26", "z27", "z28", "z29",
        "z30", "z31", "memory");
#else
  __asm__ volatile(
      "ptrue pn8.s\n"
      "ld1w {z16.s, z20.s, z24.s, z28.s}, pn8/z, [%[a0]]\n"
      "ld1w {z17.s, z21.s, z25.s, z29.s}, pn8/z, [%[a1]]\n"
      "ld1w {z18.s, z22.s, z26.s, z30.s}, pn8/z, [%[a2]]\n"
      "ld1w {z19.s, z23.s, z27.s, z31.s}, pn8/z, [%[a3]]\n"
      "mova za0h.s[%w[r], 0:3], {z16.s - z19.s}\n"
      "mova za1h.s[%w[r], 0:3], {z20.s - z23.s}\n"
      "mova za2h.s[%w[r], 0:3], {z24.s - z27.s}\n"
      "mova za3h.s[%w[r], 0:3], {z28.s - z31.s}\n"
      :
      : [a0] "r"(src), [a1] "r"(p1), [a2] "r"(p2), [a3] "r"(p3), [r] "Ucj"(r)
      : "p8", "z16", "z17", "z18", "z19", "z20", "z21", "z22", "z23", "z24", "z25", "z26", "z27", "z28", "z29",
        "z30", "z31", "memory");
#endif
}

/* Inverse of s2_rows_in4: horizontal slices r..r+3 of every tile to four rows of 64 columns. */
S2_INL void s2_rows_out4(FLOAT *dst, BLASLONG ld, uint32_t r) S2_S S2_ZA {
  FLOAT *p1 = dst + ld, *p2 = dst + 2 * ld, *p3 = dst + 3 * ld;
#ifdef DOUBLE
  __asm__ volatile(
      "ptrue pn8.d\n"
      "mova {z16.d - z19.d}, za0h.d[%w[r], 0:3]\n"
      "mova {z20.d - z23.d}, za1h.d[%w[r], 0:3]\n"
      "mova {z24.d - z27.d}, za2h.d[%w[r], 0:3]\n"
      "mova {z28.d - z31.d}, za3h.d[%w[r], 0:3]\n"
      "st1d {z16.d, z20.d, z24.d, z28.d}, pn8, [%[a0]]\n"
      "st1d {z17.d, z21.d, z25.d, z29.d}, pn8, [%[a1]]\n"
      "st1d {z18.d, z22.d, z26.d, z30.d}, pn8, [%[a2]]\n"
      "st1d {z19.d, z23.d, z27.d, z31.d}, pn8, [%[a3]]\n"
      "mova {z16.d - z19.d}, za4h.d[%w[r], 0:3]\n"
      "mova {z20.d - z23.d}, za5h.d[%w[r], 0:3]\n"
      "mova {z24.d - z27.d}, za6h.d[%w[r], 0:3]\n"
      "mova {z28.d - z31.d}, za7h.d[%w[r], 0:3]\n"
      "st1d {z16.d, z20.d, z24.d, z28.d}, pn8, [%[a0], #4, mul vl]\n"
      "st1d {z17.d, z21.d, z25.d, z29.d}, pn8, [%[a1], #4, mul vl]\n"
      "st1d {z18.d, z22.d, z26.d, z30.d}, pn8, [%[a2], #4, mul vl]\n"
      "st1d {z19.d, z23.d, z27.d, z31.d}, pn8, [%[a3], #4, mul vl]\n"
      :
      : [a0] "r"(dst), [a1] "r"(p1), [a2] "r"(p2), [a3] "r"(p3), [r] "Ucj"(r)
      : "p8", "z16", "z17", "z18", "z19", "z20", "z21", "z22", "z23", "z24", "z25", "z26", "z27", "z28", "z29",
        "z30", "z31", "memory");
#else
  __asm__ volatile(
      "ptrue pn8.s\n"
      "mova {z16.s - z19.s}, za0h.s[%w[r], 0:3]\n"
      "mova {z20.s - z23.s}, za1h.s[%w[r], 0:3]\n"
      "mova {z24.s - z27.s}, za2h.s[%w[r], 0:3]\n"
      "mova {z28.s - z31.s}, za3h.s[%w[r], 0:3]\n"
      "st1w {z16.s, z20.s, z24.s, z28.s}, pn8, [%[a0]]\n"
      "st1w {z17.s, z21.s, z25.s, z29.s}, pn8, [%[a1]]\n"
      "st1w {z18.s, z22.s, z26.s, z30.s}, pn8, [%[a2]]\n"
      "st1w {z19.s, z23.s, z27.s, z31.s}, pn8, [%[a3]]\n"
      :
      : [a0] "r"(dst), [a1] "r"(p1), [a2] "r"(p2), [a3] "r"(p3), [r] "Ucj"(r)
      : "p8", "z16", "z17", "z18", "z19", "z20", "z21", "z22", "z23", "z24", "z25", "z26", "z27", "z28", "z29",
        "z30", "z31", "memory");
#endif
}

#ifndef DOUBLE
/* 2 x 2 fp32 tiles (h = tile row): four C rows of 32 columns via strided x2 loads, one MOVA vg4 per tile. */
#define S2_X2_ARGS                                                                                               \
  : : [a0] "r"(src), [a1] "r"(p1), [a2] "r"(p2), [a3] "r"(p3), [r] "Ucj"(r)                                       \
  : "p8", "z16", "z17", "z18", "z19", "z24", "z25", "z26", "z27", "memory"
#define S2_X2_LD                                                                                                 \
  "ptrue pn8.s\n ld1w {z16.s, z24.s}, pn8/z, [%[a0]]\n ld1w {z17.s, z25.s}, pn8/z, [%[a1]]\n"                    \
  "ld1w {z18.s, z26.s}, pn8/z, [%[a2]]\n ld1w {z19.s, z27.s}, pn8/z, [%[a3]]\n"
#define S2_X2_ST                                                                                                 \
  "st1w {z16.s, z24.s}, pn8, [%[a0]]\n st1w {z17.s, z25.s}, pn8, [%[a1]]\n"                                     \
  "st1w {z18.s, z26.s}, pn8, [%[a2]]\n st1w {z19.s, z27.s}, pn8, [%[a3]]\n"
S2_INL void s2_rows_in4_x2(const float *src, BLASLONG ld, uint32_t r, int h) S2_S S2_ZA {
  const float *p1 = src + ld, *p2 = src + 2 * ld, *p3 = src + 3 * ld;
  if (h == 0)
    __asm__ volatile(S2_X2_LD "mova za0h.s[%w[r], 0:3], {z16.s - z19.s}\n mova za1h.s[%w[r], 0:3], {z24.s - z27.s}\n"
                     S2_X2_ARGS);
  else
    __asm__ volatile(S2_X2_LD "mova za2h.s[%w[r], 0:3], {z16.s - z19.s}\n mova za3h.s[%w[r], 0:3], {z24.s - z27.s}\n"
                     S2_X2_ARGS);
}
S2_INL void s2_rows_out4_x2(float *dst, BLASLONG ld, uint32_t r, int h) S2_S S2_ZA {
  const float *src = dst;
  float *p1 = dst + ld, *p2 = dst + 2 * ld, *p3 = dst + 3 * ld;
  if (h == 0)
    __asm__ volatile("ptrue pn8.s\n mova {z16.s - z19.s}, za0h.s[%w[r], 0:3]\n mova {z24.s - z27.s}, za1h.s[%w[r], 0:3]\n"
                     S2_X2_ST S2_X2_ARGS);
  else
    __asm__ volatile("ptrue pn8.s\n mova {z16.s - z19.s}, za2h.s[%w[r], 0:3]\n mova {z24.s - z27.s}, za3h.s[%w[r], 0:3]\n"
                     S2_X2_ST S2_X2_ARGS);
}
#endif

#ifdef DOUBLE
/* fp64 half-width block (h = tile row): four C rows of 32 columns into tiles 4h..4h+3, one MOVA vg4 per tile. */
S2_INL void s2_rows_in4_h(const double *src, BLASLONG ld, uint32_t r, int h) S2_S S2_ZA {
  const double *p1 = src + ld, *p2 = src + 2 * ld, *p3 = src + 3 * ld;
#define S2_H_LD                                                                                                  \
  "ptrue pn8.d\n ld1d {z16.d, z20.d, z24.d, z28.d}, pn8/z, [%[a0]]\n"                                            \
  "ld1d {z17.d, z21.d, z25.d, z29.d}, pn8/z, [%[a1]]\n ld1d {z18.d, z22.d, z26.d, z30.d}, pn8/z, [%[a2]]\n"     \
  "ld1d {z19.d, z23.d, z27.d, z31.d}, pn8/z, [%[a3]]\n"
#define S2_H_ARGS                                                                                                \
  : : [a0] "r"(src), [a1] "r"(p1), [a2] "r"(p2), [a3] "r"(p3), [r] "Ucj"(r)                                     \
  : "p8", "z16", "z17", "z18", "z19", "z20", "z21", "z22", "z23", "z24", "z25", "z26", "z27", "z28", "z29",   \
    "z30", "z31", "memory"
  if (h == 0)
    __asm__ volatile(S2_H_LD "mova za0h.d[%w[r], 0:3], {z16.d - z19.d}\n mova za1h.d[%w[r], 0:3], {z20.d - z23.d}\n"
                     "mova za2h.d[%w[r], 0:3], {z24.d - z27.d}\n mova za3h.d[%w[r], 0:3], {z28.d - z31.d}\n" S2_H_ARGS);
  else
    __asm__ volatile(S2_H_LD "mova za4h.d[%w[r], 0:3], {z16.d - z19.d}\n mova za5h.d[%w[r], 0:3], {z20.d - z23.d}\n"
                     "mova za6h.d[%w[r], 0:3], {z24.d - z27.d}\n mova za7h.d[%w[r], 0:3], {z28.d - z31.d}\n" S2_H_ARGS);
#undef S2_H_LD
#undef S2_H_ARGS
}
#endif

/* C block of TM x TN tiles into ZA: rows < mrows, columns < ncols; tile (p, q) holds rows p*VL.., columns q*VL.. */
S2_INL void s2_c_load(const int TM, const int TN, FLOAT *C, BLASLONG ldc, int mrows, int ncols, int cmode,
                      FLOAT beta) S2_S S2_ZA {
  if (cmode == S2_ZERO) { /* rows and columns outside the block hold stale values but are never stored */
    svzero_za();
    return;
  }
  int s0 = 0;
  if (TM == 1 && TN == S2_NT && cmode == S2_LOAD && ncols == TN * S2_VL)
    for (; s0 + 4 <= mrows; s0 += 4) s2_rows_in4(C + (BLASLONG)s0 * ldc, ldc, s0);
#ifndef DOUBLE
  if (TM == 2 && TN == 2 && cmode == S2_LOAD && ncols == 2 * S2_VL && mrows == 2 * S2_VL) {
    for (uint32_t r = 0; r < S2_VL; r += 4) s2_rows_in4_x2(C + (BLASLONG)r * ldc, ldc, r, 0);
    for (uint32_t r = 0; r < S2_VL; r += 4) s2_rows_in4_x2(C + (BLASLONG)(S2_VL + r) * ldc, ldc, r, 1);
    return;
  }
#else
  if (TN == 4 && cmode == S2_LOAD && ncols == 4 * S2_VL && mrows == TM * S2_VL) {
    S2_UNROLL for (int P = 0; P < TM; ++P)
      for (uint32_t r = 0; r < S2_VL; r += 4) s2_rows_in4_h(C + (BLASLONG)(P * S2_VL + r) * ldc, ldc, r, P);
    return;
  }
#endif
  const int sc = cmode == S2_SCALE;
  S2_UNROLL for (int P = 0; P < TM; ++P) {
    const int prow = mrows - P * S2_VL < S2_VL ? mrows - P * S2_VL : S2_VL;
    int s = P == 0 ? s0 : 0;
    for (; s < prow; ++s) {
      FLOAT *row = C + (BLASLONG)(P * S2_VL + s) * ldc;
      if (TN % 4 == 0) {
        S2_UNROLL for (int G = 0; G < TN / 4; ++G) {
          s2_v4 v = ncols >= (G + 1) * 4 * S2_VL ? s2_ld4(row + G * 4 * S2_VL)
                                                 : s2_ld4n(row + G * 4 * S2_VL, s2_clip(ncols - G * 4 * S2_VL, 4 * S2_VL));
          S2_UNROLL for (int Q = 0; Q < 4; ++Q) {
            s2_v x = s2_get4(v, Q);
            if (sc) x = S2_MUL(x, beta);
            s2_wrh(P * TN + G * 4 + Q, s, x);
          }
        }
      } else {
        S2_UNROLL for (int Q = 0; Q < TN; ++Q) {
          s2_v x = S2_LD1(S2_PW(Q * S2_VL, ncols), row + Q * S2_VL);
          if (sc) x = S2_MUL(x, beta);
          s2_wrh(P * TN + Q, s, x);
        }
      }
    }
  }
}

S2_INL void s2_c_store(const int TM, const int TN, FLOAT *C, BLASLONG ldc, int mrows, int ncols) S2_S S2_ZA {
  int s0 = 0;
  if (TM == 1 && TN == S2_NT && ncols == TN * S2_VL)
    for (; s0 + 4 <= mrows; s0 += 4) s2_rows_out4(C + (BLASLONG)s0 * ldc, ldc, s0);
#ifndef DOUBLE
  if (TM == 2 && TN == 2 && ncols == 2 * S2_VL && mrows == 2 * S2_VL) {
    for (uint32_t r = 0; r < S2_VL; r += 4) s2_rows_out4_x2(C + (BLASLONG)r * ldc, ldc, r, 0);
    for (uint32_t r = 0; r < S2_VL; r += 4) s2_rows_out4_x2(C + (BLASLONG)(S2_VL + r) * ldc, ldc, r, 1);
    return;
  }
#endif
  S2_UNROLL for (int P = 0; P < TM; ++P) {
    const int prow = mrows - P * S2_VL < S2_VL ? mrows - P * S2_VL : S2_VL;
    int s = P == 0 ? s0 : 0;
    for (; s + 4 <= prow; s += 4) { /* four rows per ZA move, stores predicated for partial widths */
      FLOAT *row = C + (BLASLONG)(P * S2_VL + s) * ldc;
      if (TN % 4 == 0) {
        S2_UNROLL for (int G = 0; G < TN / 4; ++G) {
          const int t = P * TN + G * 4;
          const s2_v4 q0 = s2_rdh4(t, s), q1 = s2_rdh4(t + 1, s), q2 = s2_rdh4(t + 2, s), q3 = s2_rdh4(t + 3, s);
          S2_UNROLL for (int u = 0; u < 4; ++u) {
            const s2_v4 v = S2_CREATE4(s2_get4(q0, u), s2_get4(q1, u), s2_get4(q2, u), s2_get4(q3, u));
            if (ncols >= (G + 1) * 4 * S2_VL) s2_st4(row + u * ldc + G * 4 * S2_VL, v);
            else s2_st4n(row + u * ldc + G * 4 * S2_VL, v, s2_clip(ncols - G * 4 * S2_VL, 4 * S2_VL));
          }
        }
      } else {
        S2_UNROLL for (int Q = 0; Q < TN; ++Q) {
          const svbool_t pg = S2_PW(Q * S2_VL, ncols);
          const s2_v4 q = s2_rdh4(P * TN + Q, s);
          S2_UNROLL for (int u = 0; u < 4; ++u) S2_ST1(pg, row + u * ldc + Q * S2_VL, s2_get4(q, u));
        }
      }
    }
    for (; s < prow; ++s) {
      FLOAT *row = C + (BLASLONG)(P * S2_VL + s) * ldc;
      if (TN % 4 == 0) {
        S2_UNROLL for (int G = 0; G < TN / 4; ++G) {
          const int t = P * TN + G * 4;
          s2_v4 v = S2_CREATE4(s2_rdh(t, s), s2_rdh(t + 1, s), s2_rdh(t + 2, s), s2_rdh(t + 3, s));
          if (ncols >= (G + 1) * 4 * S2_VL) s2_st4(row + G * 4 * S2_VL, v);
          else s2_st4n(row + G * 4 * S2_VL, v, s2_clip(ncols - G * 4 * S2_VL, 4 * S2_VL));
        }
      } else {
        S2_UNROLL for (int Q = 0; Q < TN; ++Q) S2_ST1(S2_PW(Q * S2_VL, ncols), row + Q * S2_VL, s2_rdh(P * TN + Q, s));
      }
    }
  }
}

/* C = C + alpha ZA (add) or C = alpha ZA from a row of TN tiles: rows < mrows, columns < ncols. */
S2_INL void s2_c_store_alpha(const int TN, FLOAT *C, BLASLONG ldc, int mrows, int ncols, FLOAT alpha, int add) S2_S
    S2_ZA {
  for (int s = 0; s < mrows; ++s) {
    FLOAT *row = C + (BLASLONG)s * ldc;
    S2_UNROLL for (int Q = 0; Q < TN; ++Q) {
      const svbool_t pg = S2_PW(Q * S2_VL, ncols);
      const s2_v z = s2_rdh(Q, s);
      S2_ST1(pg, row + Q * S2_VL, add ? S2_MLA(S2_LD1(pg, row + Q * S2_VL), z, alpha) : S2_MUL(z, alpha));
    }
  }
}

#pragma clang attribute pop

#endif
