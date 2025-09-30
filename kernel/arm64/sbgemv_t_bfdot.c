/***************************************************************************
Copyright (c) 2025, The OpenBLAS Project
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
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#include <arm_neon.h>
#include "common.h"

#ifdef BGEMM
#define INNER_FLOAT bfloat16_t
#define TO32(x) vcvtah_f32_bf16(x)
#define FROM32(x) vcvth_bf16_f32(x)
#else
#define INNER_FLOAT float
#define TO32(x) x
#define FROM32(x) x
#endif

int CNAME(BLASLONG m, BLASLONG n, FLOAT alpha, bfloat16 *a, BLASLONG lda, bfloat16 *x, BLASLONG incx, FLOAT beta, FLOAT *y_in, BLASLONG incy)
{
  if (m < 1 || n < 1) return(0);

  BLASLONG i;
  BLASLONG ix,iy;
  BLASLONG j;
  bfloat16_t *a_ptr;
  bfloat16_t *x_ptr;
  float temp, temp0, temp1, temp2, temp3;

#ifdef BGEMM
  bfloat16_t alpha_bf16, beta_bf16;
  memcpy(&alpha_bf16, &alpha, sizeof(bfloat16_t));
  memcpy(&beta_bf16, &beta, sizeof(bfloat16_t));
  float alpha_f32 = vcvtah_f32_bf16(alpha_bf16);
  float beta_f32 = vcvtah_f32_bf16(beta_bf16);
#else
  float alpha_f32 = alpha;
  float beta_f32 = beta;
#endif
  INNER_FLOAT *y = (INNER_FLOAT *)y_in;
  INNER_FLOAT *y_ptr;

  iy = 0;
  a_ptr = (bfloat16_t*)(a);
  x_ptr = (bfloat16_t*)(x);

  if (incx == 1) {
    BLASLONG width = n / 4;

    bfloat16_t *a0_ptr = a_ptr + lda * width * 0;
    bfloat16_t *a1_ptr = a_ptr + lda * width * 1;
    bfloat16_t *a2_ptr = a_ptr + lda * width * 2;
    bfloat16_t *a3_ptr = a_ptr + lda * width * 3;

    INNER_FLOAT *y0_ptr = y + incy * width * 0;
    INNER_FLOAT *y1_ptr = y + incy * width * 1;
    INNER_FLOAT *y2_ptr = y + incy * width * 2;
    INNER_FLOAT *y3_ptr = y + incy * width * 3;

    for (j = 0; j < width; j++) {
      float32x4_t temp0_vec = vdupq_n_f32(0.0f);
      float32x4_t temp1_vec = vdupq_n_f32(0.0f);
      float32x4_t temp2_vec = vdupq_n_f32(0.0f);
      float32x4_t temp3_vec = vdupq_n_f32(0.0f);

      i = 0;
      while (i + 7 < m) {
        bfloat16x8_t x_vec = vld1q_bf16(x_ptr + i);

        bfloat16x8_t a0_vec = vld1q_bf16(a0_ptr + i);
        bfloat16x8_t a1_vec = vld1q_bf16(a1_ptr + i);
        bfloat16x8_t a2_vec = vld1q_bf16(a2_ptr + i);
        bfloat16x8_t a3_vec = vld1q_bf16(a3_ptr + i);

        temp0_vec = vbfdotq_f32(temp0_vec, a0_vec, x_vec);
        temp1_vec = vbfdotq_f32(temp1_vec, a1_vec, x_vec);
        temp2_vec = vbfdotq_f32(temp2_vec, a2_vec, x_vec);
        temp3_vec = vbfdotq_f32(temp3_vec, a3_vec, x_vec);

        i += 8;
      }
      if (i + 3 < m) {
        float32x2_t t0 = vdup_n_f32(0.0f);
        float32x2_t t1 = vdup_n_f32(0.0f);
        float32x2_t t2 = vdup_n_f32(0.0f);
        float32x2_t t3 = vdup_n_f32(0.0f);

        bfloat16x4_t x_vec = vld1_bf16(x_ptr + i);

        bfloat16x4_t a0_vec = vld1_bf16(a0_ptr + i);
        bfloat16x4_t a1_vec = vld1_bf16(a1_ptr + i);
        bfloat16x4_t a2_vec = vld1_bf16(a2_ptr + i);
        bfloat16x4_t a3_vec = vld1_bf16(a3_ptr + i);

        t0 = vbfdot_f32(t0, a0_vec, x_vec);
        t1 = vbfdot_f32(t1, a1_vec, x_vec);
        t2 = vbfdot_f32(t2, a2_vec, x_vec);
        t3 = vbfdot_f32(t3, a3_vec, x_vec);

        float32x2_t temp0_vec_low = vget_low_f32(temp0_vec);
        float32x2_t temp1_vec_low = vget_low_f32(temp1_vec);
        float32x2_t temp2_vec_low = vget_low_f32(temp2_vec);
        float32x2_t temp3_vec_low = vget_low_f32(temp3_vec);

        temp0_vec = vcombine_f32(vadd_f32(t0, temp0_vec_low), vget_high_f32(temp0_vec));
        temp1_vec = vcombine_f32(vadd_f32(t1, temp1_vec_low), vget_high_f32(temp1_vec));
        temp2_vec = vcombine_f32(vadd_f32(t2, temp2_vec_low), vget_high_f32(temp2_vec));
        temp3_vec = vcombine_f32(vadd_f32(t3, temp3_vec_low), vget_high_f32(temp3_vec));

        i += 4;
      }

      if (beta_f32 == 0.0f) {
        temp0 = alpha_f32 * vaddvq_f32(temp0_vec);
        temp1 = alpha_f32 * vaddvq_f32(temp1_vec);
        temp2 = alpha_f32 * vaddvq_f32(temp2_vec);
        temp3 = alpha_f32 * vaddvq_f32(temp3_vec);
      } else {
        temp0 = alpha_f32 * vaddvq_f32(temp0_vec) + beta_f32 * TO32(y0_ptr[iy]);
        temp1 = alpha_f32 * vaddvq_f32(temp1_vec) + beta_f32 * TO32(y1_ptr[iy]);
        temp2 = alpha_f32 * vaddvq_f32(temp2_vec) + beta_f32 * TO32(y2_ptr[iy]);
        temp3 = alpha_f32 * vaddvq_f32(temp3_vec) + beta_f32 * TO32(y3_ptr[iy]);
      }

      for (; i < m; ++i) {
        temp0 = temp0 + alpha_f32 * vcvtah_f32_bf16(a0_ptr[i]) * vcvtah_f32_bf16(x_ptr[i]);
        temp1 = temp1 + alpha_f32 * vcvtah_f32_bf16(a1_ptr[i]) * vcvtah_f32_bf16(x_ptr[i]);
        temp2 = temp2 + alpha_f32 * vcvtah_f32_bf16(a2_ptr[i]) * vcvtah_f32_bf16(x_ptr[i]);
        temp3 = temp3 + alpha_f32 * vcvtah_f32_bf16(a3_ptr[i]) * vcvtah_f32_bf16(x_ptr[i]);
      }

      y0_ptr[iy] = FROM32(temp0);
      y1_ptr[iy] = FROM32(temp1);
      y2_ptr[iy] = FROM32(temp2);
      y3_ptr[iy] = FROM32(temp3);

      iy += incy;

      a0_ptr += lda;
      a1_ptr += lda;
      a2_ptr += lda;
      a3_ptr += lda;
    }

    a_ptr = a3_ptr;
    y_ptr = y3_ptr;
    for (j = width * 4; j < n; j++) {
      float32x4_t temp0_vec = vdupq_n_f32(0.0f);
      i = 0;
      while (i + 7 < m) {
        bfloat16x8_t x_vec = vld1q_bf16(x_ptr + i);
        bfloat16x8_t a0_vec = vld1q_bf16(a_ptr + i);
        temp0_vec = vbfdotq_f32(temp0_vec, a0_vec, x_vec);

        i += 8;
      }
      if (i + 3 < m) {
        float32x2_t t0 = vdup_n_f32(0.0f);
        bfloat16x4_t x_vec = vld1_bf16(x_ptr + i);
        bfloat16x4_t a0_vec = vld1_bf16(a_ptr + i);

        t0 = vbfdot_f32(t0, a0_vec, x_vec);
        float32x2_t temp0_vec_low = vget_low_f32(temp0_vec);
        temp0_vec = vcombine_f32(vadd_f32(t0, temp0_vec_low), vget_high_f32(temp0_vec));

        i += 4;
      }

      if (beta_f32 == 0.0f) {
        temp = alpha_f32 * vaddvq_f32(temp0_vec);
      } else {
        temp = alpha_f32 * vaddvq_f32(temp0_vec) + beta_f32 * TO32(y_ptr[iy]);
      }
      for (; i < m; ++i) {
        temp += alpha_f32 * vcvtah_f32_bf16(a_ptr[i]) * vcvtah_f32_bf16(x_ptr[i]);
      }
      y_ptr[iy] = FROM32(temp);

      iy += incy;

      a_ptr += lda;
    }
    return(0);
  }

  for (j = 0; j < n; j++) {
    temp = 0.0;
    ix = 0;
    for (i = 0; i < m; i++) {
      temp += vcvtah_f32_bf16(a_ptr[i]) * vcvtah_f32_bf16(x_ptr[ix]);
      ix += incx;
    }
    if (beta_f32 == 0.0f) {
      y[iy] = FROM32(alpha_f32 * temp);
    }
    else {
      y[iy] = FROM32(alpha_f32 * temp + beta_f32 * TO32(y[iy]));
    }
    iy += incy;
    a_ptr += lda;
  }
  return (0);
}
