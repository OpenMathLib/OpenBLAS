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


/* GEMV T, y += alpha A^T x: SME2 with 16 columns at a time accumulated in ZA; NEON for small problems and without SME2. */

#include "common.h"
#include "sme2_gemv.h"

/* Smallest m and m * n for the SME2 path, which takes the columns in blocks of 16. */
#ifndef GVT_SME2_MIN
#define GVT_SME2_MIN 4096
#endif
#ifndef GVT_SME2_MINM
#define GVT_SME2_MINM 48
#endif

/* y += alpha A^T x with a contiguous x, four columns at a time. */
static void gvt_neon(BLASLONG m, BLASLONG n, FLOAT alpha, const FLOAT *a, BLASLONG lda, const FLOAT *x, FLOAT *y,
                     BLASLONG incy) {
  BLASLONG j = 0;
  for (; j + 4 <= n; j += 4) {
    const FLOAT *a0 = a + j * lda, *a1 = a0 + lda, *a2 = a1 + lda, *a3 = a2 + lda;
    gv_v s0 = GV_ZERO(), s1 = s0, s2 = s0, s3 = s0, u0 = s0, u1 = s0, u2 = s0, u3 = s0;
    BLASLONG i = 0;
    for (; i + 2 * GV_NV <= m; i += 2 * GV_NV) {
      const gv_v x0 = GV_LD(x + i), x1 = GV_LD(x + i + GV_NV);
      s0 = GV_FMA(s0, GV_LD(a0 + i), x0);
      u0 = GV_FMA(u0, GV_LD(a0 + i + GV_NV), x1);
      s1 = GV_FMA(s1, GV_LD(a1 + i), x0);
      u1 = GV_FMA(u1, GV_LD(a1 + i + GV_NV), x1);
      s2 = GV_FMA(s2, GV_LD(a2 + i), x0);
      u2 = GV_FMA(u2, GV_LD(a2 + i + GV_NV), x1);
      s3 = GV_FMA(s3, GV_LD(a3 + i), x0);
      u3 = GV_FMA(u3, GV_LD(a3 + i + GV_NV), x1);
    }
    FLOAT r0 = GV_ADDV(s0, u0), r1 = GV_ADDV(s1, u1), r2 = GV_ADDV(s2, u2), r3 = GV_ADDV(s3, u3);
    for (; i < m; i++) {
      r0 += a0[i] * x[i];
      r1 += a1[i] * x[i];
      r2 += a2[i] * x[i];
      r3 += a3[i] * x[i];
    }
    y[j * incy] += alpha * r0;
    y[(j + 1) * incy] += alpha * r1;
    y[(j + 2) * incy] += alpha * r2;
    y[(j + 3) * incy] += alpha * r3;
  }
  for (; j < n; j++) {
    const FLOAT *aj = a + j * lda;
    gv_v s0 = GV_ZERO(), u0 = s0;
    BLASLONG i = 0;
    for (; i + 2 * GV_NV <= m; i += 2 * GV_NV) {
      s0 = GV_FMA(s0, GV_LD(aj + i), GV_LD(x + i));
      u0 = GV_FMA(u0, GV_LD(aj + i + GV_NV), GV_LD(x + i + GV_NV));
    }
    FLOAT r0 = GV_ADDV(s0, u0);
    for (; i < m; i++) r0 += aj[i] * x[i];
    y[j * incy] += alpha * r0;
  }
}

#ifdef HAVE_SME2_GEMM
#pragma clang attribute push(__attribute__((target("sme2"))), apply_to = function)

#ifdef DOUBLE
#define GVT_T "d"
#define GVT_LD "ld1d"
#define GVT_SH "3"
#else
#define GVT_T "s"
#define GVT_LD "ld1w"
#define GVT_SH "2"
#endif

/* Full block of 16 columns, ZA row s + 16 k accumulates column 4 s + k over rows i..i+GV_GL-1; strided loads feed the FMLA without moves. */
S2_INL void gvt_chunk(BLASLONG i, BLASLONG m, const FLOAT *x, const FLOAT *a, BLASLONG ldb) S2_S S2_ZA {
  __asm__ volatile(
      "whilelt pn8." GVT_T ", %[i], %[m], vlx4\n"
      GVT_LD " {z0." GVT_T "-z3." GVT_T "}, pn8/z, [%[x], %[i], lsl #" GVT_SH "]\n"
      "mov w8, #0\n"
      "mov x10, %[a]\n"
      "1:\n"
      "add x11, x10, %[ldb]\n"
      "add x12, x11, %[ldb]\n"
      "add x13, x12, %[ldb]\n"
      GVT_LD " {z16." GVT_T ", z20." GVT_T ", z24." GVT_T ", z28." GVT_T "}, pn8/z, [x10, %[i], lsl #" GVT_SH "]\n"
      GVT_LD " {z17." GVT_T ", z21." GVT_T ", z25." GVT_T ", z29." GVT_T "}, pn8/z, [x11, %[i], lsl #" GVT_SH "]\n"
      GVT_LD " {z18." GVT_T ", z22." GVT_T ", z26." GVT_T ", z30." GVT_T "}, pn8/z, [x12, %[i], lsl #" GVT_SH "]\n"
      GVT_LD " {z19." GVT_T ", z23." GVT_T ", z27." GVT_T ", z31." GVT_T "}, pn8/z, [x13, %[i], lsl #" GVT_SH "]\n"
      "add x10, x13, %[ldb]\n"
      "fmla za." GVT_T "[w8, 0, vgx4], {z16." GVT_T "-z19." GVT_T "}, z0." GVT_T "\n"
      "fmla za." GVT_T "[w8, 0, vgx4], {z20." GVT_T "-z23." GVT_T "}, z1." GVT_T "\n"
      "fmla za." GVT_T "[w8, 0, vgx4], {z24." GVT_T "-z27." GVT_T "}, z2." GVT_T "\n"
      "fmla za." GVT_T "[w8, 0, vgx4], {z28." GVT_T "-z31." GVT_T "}, z3." GVT_T "\n"
      "add w8, w8, #1\n"
      "cmp w8, #4\n"
      "b.lt 1b\n"
      :
      : [i] "r"(i), [m] "r"(m), [x] "r"(x), [a] "r"(a), [ldb] "r"(ldb)
      : "x8", "x10", "x11", "x12", "x13", "p8", "z0", "z1", "z2", "z3", "z16", "z17", "z18", "z19", "z20", "z21", "z22",
        "z23", "z24", "z25", "z26", "z27", "z28", "z29", "z30", "z31", "za", "cc", "memory");
}

#define GVT_SUM4(r) GV_ADD(GV_ADD(S2_GET4(r, 0), S2_GET4(r, 1)), GV_ADD(S2_GET4(r, 2), S2_GET4(r, 3)))
#define GVT_VSUM(t, s, acc) acc = GV_ADD(acc, GVT_SUM4(S2_RDV4(t, s)))

/* y[0..16) += alpha * sums of the 16 column vectors in ZA groups 0..3 (column 4 s + k in vector k of group s). */
S2_INL void gvt_reduce(FLOAT alpha, FLOAT *y) S2_S S2_ZA {
  const s2_v4 r0 = GV_RDZ4(0), r1 = GV_RDZ4(1), r2 = GV_RDZ4(2), r3 = GV_RDZ4(3);
  s2_v t0 = S2_ZEROV();
#ifdef DOUBLE
  s2_v t1 = S2_ZEROV();
  S2_WRH4(0, 0, r0);
  S2_WRH4(0, 4, r1);
  S2_WRH4(1, 0, r2);
  S2_WRH4(1, 4, r3);
  GVT_VSUM(0, 0, t0);
  GVT_VSUM(0, 4, t0);
  GVT_VSUM(1, 0, t1);
  GVT_VSUM(1, 4, t1);
  S2_ST1(S2_PT(), y + S2_VL, GV_MLAV(S2_LD1(S2_PT(), y + S2_VL), t1, alpha));
#else
  s2_v t1 = S2_ZEROV();
  S2_WRH4(0, 0, r0);
  S2_WRH4(0, 4, r1);
  S2_WRH4(0, 8, r2);
  S2_WRH4(0, 12, r3);
  GVT_VSUM(0, 0, t0);
  GVT_VSUM(0, 4, t1);
  GVT_VSUM(0, 8, t0);
  GVT_VSUM(0, 12, t1);
  t0 = GV_ADD(t0, t1);
#endif
  S2_ST1(S2_PT(), y, GV_MLAV(S2_LD1(S2_PT(), y), t0, alpha));
}

/* n is a multiple of 16. */
__arm_new("za") __arm_locally_streaming static void gvt_sme2(BLASLONG m, BLASLONG n, FLOAT alpha, const FLOAT *a,
                                                             BLASLONG lda, const FLOAT *x, FLOAT *y) {
  for (BLASLONG j0 = 0; j0 < n; j0 += 16) {
    svzero_za();
    for (BLASLONG i = 0; i < m; i += GV_GL) gvt_chunk(i, m, x, a + j0 * lda, lda * (BLASLONG)sizeof(FLOAT));
    gvt_reduce(alpha, y + j0);
  }
}

#pragma clang attribute pop
#endif

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x,
          FLOAT *y, BLASLONG inc_y, FLOAT *buffer) {
  if (m < 1 || n < 1 || alpha == ZERO) return 0;
  FLOAT *xx = x;
  if (inc_x != 1) {
    xx = buffer;
    gv_gather(m, x, inc_x, xx);
    buffer += m;
  }
#ifdef HAVE_SME2_GEMM
  const BLASLONG n16 = n & ~(BLASLONG)15;
  if (n16 && m >= GVT_SME2_MINM && m * n >= GVT_SME2_MIN && s2_usable()) {
    FLOAT *yy = y;
    if (inc_y != 1) {
      yy = buffer;
      gv_gather(n16, y, inc_y, yy);
    }
    gvt_sme2(m, n16, alpha, a, lda, xx, yy);
    if (inc_y != 1) gv_scatter(n16, yy, y, inc_y);
    a += n16 * lda;
    y += n16 * inc_y;
    n -= n16;
  }
#endif
  gvt_neon(m, n, alpha, a, lda, xx, y, inc_y);
  return 0;
}
