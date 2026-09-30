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

/* Copies into the A panels of the SME2 kernels: m rows of w <= 64 columns, element (i, c) to b[i * w + c]. */

#ifndef SME2_GEMM_COPY_H
#define SME2_GEMM_COPY_H

#include <arm_neon.h>
#include <string.h>

#define S2C_W GEMM_DEFAULT_UNROLL_M

static inline BLASLONG s2c_clip(BLASLONG x, BLASLONG lo, BLASLONG hi) { return x < lo ? lo : (x > hi ? hi : x); }

/* Element (i, c) at a[c + i * lda]. */
static inline void s2c_rows(BLASLONG m, BLASLONG w, const FLOAT *a, BLASLONG lda, FLOAT *b) {
  for (BLASLONG i = 0; i < m; ++i) memcpy(b + i * w, a + i * lda, w * sizeof(FLOAT));
}

/* Element (i, c) at a[i + c * lda], transposed in registers 4 x 4 (fp32) or 2 x 2 (fp64) at a time. */
static inline void s2c_cols(BLASLONG m, BLASLONG w, const FLOAT *a, BLASLONG lda, FLOAT *b) {
  BLASLONG c = 0;
#ifdef DOUBLE
  for (; c + 2 <= w; c += 2) {
    const FLOAT *a0 = a + c * lda, *a1 = a0 + lda;
    BLASLONG i = 0;
    for (; i + 2 <= m; i += 2) {
      const float64x2_t x0 = vld1q_f64(a0 + i), x1 = vld1q_f64(a1 + i);
      vst1q_f64(b + i * w + c, vtrn1q_f64(x0, x1));
      vst1q_f64(b + (i + 1) * w + c, vtrn2q_f64(x0, x1));
    }
    for (; i < m; ++i) {
      b[i * w + c] = a0[i];
      b[i * w + c + 1] = a1[i];
    }
  }
#else
  for (; c + 4 <= w; c += 4) {
    const FLOAT *a0 = a + c * lda, *a1 = a0 + lda, *a2 = a1 + lda, *a3 = a2 + lda;
    BLASLONG i = 0;
    for (; i + 4 <= m; i += 4) {
      const float32x4_t x0 = vld1q_f32(a0 + i), x1 = vld1q_f32(a1 + i), x2 = vld1q_f32(a2 + i),
                        x3 = vld1q_f32(a3 + i);
      const float64x2_t t0 = vreinterpretq_f64_f32(vtrn1q_f32(x0, x1)), t1 = vreinterpretq_f64_f32(vtrn2q_f32(x0, x1)),
                        t2 = vreinterpretq_f64_f32(vtrn1q_f32(x2, x3)), t3 = vreinterpretq_f64_f32(vtrn2q_f32(x2, x3));
      vst1q_f32(b + i * w + c, vreinterpretq_f32_f64(vtrn1q_f64(t0, t2)));
      vst1q_f32(b + (i + 1) * w + c, vreinterpretq_f32_f64(vtrn1q_f64(t1, t3)));
      vst1q_f32(b + (i + 2) * w + c, vreinterpretq_f32_f64(vtrn2q_f64(t0, t2)));
      vst1q_f32(b + (i + 3) * w + c, vreinterpretq_f32_f64(vtrn2q_f64(t1, t3)));
    }
    for (; i < m; ++i) {
      b[i * w + c] = a0[i];
      b[i * w + c + 1] = a1[i];
      b[i * w + c + 2] = a2[i];
      b[i * w + c + 3] = a3[i];
    }
  }
#endif
  for (; c < w; ++c)
    for (BLASLONG i = 0; i < m; ++i) b[i * w + c] = a[i + c * lda];
}

static inline void s2c_zero(BLASLONG m, BLASLONG w, FLOAT *b) {
  if (m > 0) memset(b, 0, (size_t)(m * w) * sizeof(FLOAT));
}

#endif
