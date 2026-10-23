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


/* Shared parts of the SME2 GEMV kernels: NEON helpers for small problems and CPUs without SME2, ZA macros. */
#ifndef SME2_GEMV_H
#define SME2_GEMV_H

#include <arm_neon.h>
#include "sme2_gemm_detect.h"

#ifdef DOUBLE
#define GV_NV 2
typedef float64x2_t gv_v;
#define GV_LD(p) vld1q_f64(p)
#define GV_ST(p, v) vst1q_f64(p, v)
#define GV_FMA(c, a, b) vfmaq_f64(c, a, b)
#define GV_FMAN(c, a, s) vfmaq_n_f64(c, a, s)
#define GV_ZERO() vdupq_n_f64(0)
#define GV_ADDV(a, b) vaddvq_f64(vaddq_f64(a, b))
#else
#define GV_NV 4
typedef float32x4_t gv_v;
#define GV_LD(p) vld1q_f32(p)
#define GV_ST(p, v) vst1q_f32(p, v)
#define GV_FMA(c, a, b) vfmaq_f32(c, a, b)
#define GV_FMAN(c, a, s) vfmaq_n_f32(c, a, s)
#define GV_ZERO() vdupq_n_f32(0)
#define GV_ADDV(a, b) vaddvq_f32(vaddq_f32(a, b))
#endif

/* Copy n elements with stride inc (the pointer as the level-2 interface passes it) to a contiguous d, and back. */
static inline void gv_gather(BLASLONG n, const FLOAT *s, BLASLONG inc, FLOAT *d) {
  for (BLASLONG i = 0; i < n; i++) d[i] = s[i * inc];
}
static inline void gv_scatter(BLASLONG n, const FLOAT *s, FLOAT *d, BLASLONG inc) {
  for (BLASLONG i = 0; i < n; i++) d[i * inc] = s[i];
}

#ifdef HAVE_SME2_GEMM
#include "sme2_gemm_tile.h"

#ifdef DOUBLE
#define GV_RDZ4(s) svread_za64_f64_vg1x4(s)
#define GV_WRZ4(s, v) svwrite_za64_f64_vg1x4(s, v)
#define GV_MLAL4(s, a, z, l) svmla_lane_za64_f64_vg1x4(s, a, z, l)
#define GV_LDRQ(p, q) svld1rq_f64(p, q)
#define GV_ADD(a, b) svadd_f64_x(svptrue_b64(), a, b)
#define GV_MLAV(c, v, s) svmla_n_f64_x(svptrue_b64(), c, v, s)
#else
#define GV_RDZ4(s) svread_za32_f32_vg1x4(s)
#define GV_WRZ4(s, v) svwrite_za32_f32_vg1x4(s, v)
#define GV_MLAL4(s, a, z, l) svmla_lane_za32_f32_vg1x4(s, a, z, l)
#define GV_LDRQ(p, q) svld1rq_f32(p, q)
#define GV_ADD(a, b) svadd_f32_x(svptrue_b32(), a, b)
#define GV_MLAV(c, v, s) svmla_n_f32_x(svptrue_b32(), c, v, s)
#endif

/* ZA vg1x4 groups: 16 slice indices of 4 vectors each, 64 vectors in all. */
#define GV_G 16
#define GV_GL (4 * S2_VL)
#endif

#endif
