/***************************************************************************
 * Copyright (c) 2025, The OpenBLAS Project
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in
 * the documentation and/or other materials provided with the
 * distribution.
 * 3. Neither the name of the OpenBLAS project nor the names of
 * its contributors may be used to endorse or promote products
 * derived from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * *****************************************************************************/

#include "common.h"

#include <arm_neon.h>

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT beta_in, IFLOAT *dummy2,
          BLASLONG dummy3, IFLOAT *dummy4, BLASLONG dummy5, FLOAT *c,
          BLASLONG ldc) {
    BLASLONG i, j;
    BLASLONG chunk, remain;
    
    bfloat16_t *ptr_c, *ptr_c0;

    bfloat16x8_t x0, z0;
    float32x4_t y0, y1;
      
    float x;
    bfloat16_t z;

    bfloat16_t zero_bf16 = vcvth_bf16_f32(0.0f); 
    bfloat16x8_t zeros = vdupq_n_bf16(zero_bf16);

    bfloat16_t beta_bf16;
    memcpy(&beta_bf16, &beta_in, sizeof(bfloat16_t));
    float beta = vcvtah_f32_bf16(beta_bf16);
    float32x4_t beta_neon = vdupq_n_f32(beta);

    ptr_c = (bfloat16_t *)c;
    
    chunk = m >> 3;
    remain = m & 7;

    if (beta == 0.0f){
      for (j = 0; j < n; j ++){
        ptr_c0 = ptr_c;
        ptr_c += ldc;

        for (i = 0; i < chunk; i ++){
          vst1q_bf16(ptr_c0, zeros);
          ptr_c0 += 8;
        }

        for (i = 0; i < remain; i ++){
          ptr_c0[0] = zero_bf16;
          ptr_c0 ++;
        }
      }
    } else {
      for (j = 0; j < n; j ++){
        ptr_c0 = ptr_c;
        ptr_c += ldc;

        for (i = 0; i < chunk; i ++){
          x0 = vld1q_bf16(ptr_c0);

          y0 = vcvtq_low_f32_bf16(x0);
          y1 = vcvtq_high_f32_bf16(x0);

          y0 = vmulq_f32(y0, beta_neon);
          y1 = vmulq_f32(y1, beta_neon);

          z0 = vcvtq_low_bf16_f32(y0);
          z0 = vcvtq_high_bf16_f32(z0, y1);
          
          vst1q_bf16(ptr_c0, z0);

          ptr_c0 += 8;
        }

        for (i = 0; i < remain; i ++){
          x = vcvtah_f32_bf16(ptr_c0[0]);
          z = vcvth_bf16_f32(x * beta);

          ptr_c0[0] = z;
          ptr_c0 ++;
        }
      } 
    }
    return 0;
};
