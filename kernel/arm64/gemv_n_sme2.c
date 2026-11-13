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


/* GEMV N, y += alpha A x: SME2 with up to 64 vectors of y updated in ZA; NEON for small problems and without SME2. */

#include "common.h"
#include "sme2_gemv.h"

/* Smallest m and m * n for the SME2 path, including the streaming mode switch. */
#ifndef GVN_SME2_MIN
#define GVN_SME2_MIN 1024
#endif
#ifndef GVN_SME2_MINM
#define GVN_SME2_MINM 16
#endif

#ifndef DOUBLE
/* fp32 small problems and CPUs without SME2: the NEON kernel of other arm64 targets. */
#pragma push_macro("CNAME")
#undef CNAME
#define CNAME sgemv_n_neon_vortexm4
#include "sgemv_n_neon.c"
#pragma pop_macro("CNAME")
#else
/* y (contiguous) += alpha A x, four columns at a time. */
static void gvn_neon(BLASLONG m, BLASLONG n, FLOAT alpha, const FLOAT *a, BLASLONG lda, const FLOAT *x, BLASLONG incx,
                     FLOAT *y) {
  BLASLONG j = 0;
  for (; j + 4 <= n; j += 4) {
    const FLOAT *a0 = a + j * lda, *a1 = a0 + lda, *a2 = a1 + lda, *a3 = a2 + lda;
    const FLOAT x0 = alpha * x[j * incx], x1 = alpha * x[(j + 1) * incx];
    const FLOAT x2 = alpha * x[(j + 2) * incx], x3 = alpha * x[(j + 3) * incx];
    BLASLONG i = 0;
    for (; i + 2 * GV_NV <= m; i += 2 * GV_NV) {
      gv_v y0 = GV_LD(y + i), y1 = GV_LD(y + i + GV_NV);
      y0 = GV_FMAN(y0, GV_LD(a0 + i), x0);
      y1 = GV_FMAN(y1, GV_LD(a0 + i + GV_NV), x0);
      y0 = GV_FMAN(y0, GV_LD(a1 + i), x1);
      y1 = GV_FMAN(y1, GV_LD(a1 + i + GV_NV), x1);
      y0 = GV_FMAN(y0, GV_LD(a2 + i), x2);
      y1 = GV_FMAN(y1, GV_LD(a2 + i + GV_NV), x2);
      y0 = GV_FMAN(y0, GV_LD(a3 + i), x3);
      y1 = GV_FMAN(y1, GV_LD(a3 + i + GV_NV), x3);
      GV_ST(y + i, y0);
      GV_ST(y + i + GV_NV, y1);
    }
    for (; i < m; i++) y[i] += a0[i] * x0 + a1[i] * x1 + a2[i] * x2 + a3[i] * x3;
  }
  for (; j < n; j++) {
    const FLOAT *aj = a + j * lda;
    const FLOAT xj = alpha * x[j * incx];
    for (BLASLONG i = 0; i < m; i++) y[i] += aj[i] * xj;
  }
}
#endif

#ifdef HAVE_SME2_GEMM
#pragma clang attribute push(__attribute__((target("sme2"))), apply_to = function)

/* Columns per 128-bit segment of x, used by the indexed FMLA. */
#ifdef DOUBLE
#define GVN_NQ 2
#define GVN_LANES(g, pc, aj, lda, xq, r) \
  GV_MLAL4(g, S2_LD4(pc, aj), xq, 0);    \
  if (r > 1) GV_MLAL4(g, S2_LD4(pc, (aj) + (lda)), xq, 1);
#else
#define GVN_NQ 4
#define GVN_LANES(g, pc, aj, lda, xq, r)                          \
  GV_MLAL4(g, S2_LD4(pc, aj), xq, 0);                             \
  if (r > 1) GV_MLAL4(g, S2_LD4(pc, (aj) + (lda)), xq, 1);         \
  if (r > 2) GV_MLAL4(g, S2_LD4(pc, (aj) + 2 * (lda)), xq, 2);     \
  if (r > 3) GV_MLAL4(g, S2_LD4(pc, (aj) + 3 * (lda)), xq, 3);
#endif

/* ZA groups 0..ng-1 += A (alpha x) for a block of mb rows, column after column as the reference BLAS; full: mb = GV_G * GV_GL. */
S2_INL void gvn_block(int ng, int full, BLASLONG mb, BLASLONG n, FLOAT alpha, const FLOAT *a, BLASLONG lda,
                      const FLOAT *x) S2_S S2_ZA {
  BLASLONG j = 0;
  for (; j + GVN_NQ <= n; j += GVN_NQ) {
    const s2_v xq = S2_MUL(GV_LDRQ(S2_PT(), x + j), alpha);
    const FLOAT *aj = a + j * lda;
    for (int g = 0; g < ng; g++) {
      const svcount_t pc = full ? S2_CT() : S2_CW(g * GV_GL, mb);
      GVN_LANES(g, pc, aj + g * GV_GL, lda, xq, GVN_NQ)
    }
  }
  if (j < n) {
    const int r = (int)(n - j);
    const s2_v xq = S2_MUL(GV_LDRQ(S2_PW(j, n), x + j), alpha);
    const FLOAT *aj = a + j * lda;
    for (int g = 0; g < ng; g++) {
      const svcount_t pc = full ? S2_CT() : S2_CW(g * GV_GL, mb);
      GVN_LANES(g, pc, aj + g * GV_GL, lda, xq, r)
    }
  }
}

__arm_new("za") __arm_locally_streaming static void gvn_sme2(BLASLONG m, BLASLONG n, FLOAT alpha, const FLOAT *a,
                                                             BLASLONG lda, const FLOAT *x, FLOAT *y) {
  for (BLASLONG i0 = 0; i0 < m; i0 += GV_G * GV_GL) {
    const BLASLONG mb = m - i0 < GV_G * GV_GL ? m - i0 : GV_G * GV_GL;
    const int ng = (int)((mb + GV_GL - 1) / GV_GL);
    for (int g = 0; g < ng; g++) GV_WRZ4(g, S2_LD4(S2_CW(g * GV_GL, mb), y + i0 + g * GV_GL));
    if (ng == GV_G && mb == GV_G * GV_GL)
      gvn_block(GV_G, 1, mb, n, alpha, a + i0, lda, x);
    else
      gvn_block(ng, 0, mb, n, alpha, a + i0, lda, x);
    for (int g = 0; g < ng; g++) S2_ST4(S2_CW(g * GV_GL, mb), y + i0 + g * GV_GL, GV_RDZ4(g));
  }
}

#pragma clang attribute pop
#endif

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x,
          FLOAT *y, BLASLONG inc_y, FLOAT *buffer) {
  if (m < 1 || n < 1 || alpha == ZERO) return 0;
  FLOAT *yy = y;
#ifdef HAVE_SME2_GEMM
  if (m >= GVN_SME2_MINM && m * n >= GVN_SME2_MIN && s2_usable()) {
    FLOAT *xx = x;
    if (inc_x != 1) {
      xx = buffer;
      gv_gather(n, x, inc_x, xx);
      buffer += n;
    }
    if (inc_y != 1) {
      yy = buffer;
      gv_gather(m, y, inc_y, yy);
    }
    gvn_sme2(m, n, alpha, a, lda, xx, yy);
    if (inc_y != 1) gv_scatter(m, yy, y, inc_y);
    return 0;
  }
#endif
#ifndef DOUBLE
  return sgemv_n_neon_vortexm4(m, n, dummy1, alpha, a, lda, x, inc_x, y, inc_y, buffer);
#else
  if (inc_y != 1) {
    yy = buffer;
    gv_gather(m, y, inc_y, yy);
  }
  gvn_neon(m, n, alpha, a, lda, x, inc_x, yy);
  if (inc_y != 1) gv_scatter(m, yy, y, inc_y);
  return 0;
#endif
}
