/***************************************************************************
Copyright (c) 2024, The OpenBLAS Project
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

#include <arm_sve.h>
#include "common.h"

#ifdef DOUBLE
#define SV_COUNT svcntd
#define SV_TYPE svfloat64_t
#define SV_TRUE svptrue_b64
#define SV_WHILE svwhilelt_b64
#define SV_DUP svdup_f64
#else
#define SV_COUNT svcntw
#define SV_TYPE svfloat32_t
#define SV_TRUE svptrue_b32
#define SV_WHILE svwhilelt_b32
#define SV_DUP svdup_f32
#endif

int CNAME(BLASLONG M, BLASLONG N, BLASLONG dummy1, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *buffer)
{
  const uint64_t v_size = SV_COUNT();
  const uint64_t v_size2 = v_size * 2;
  const svbool_t pg_true = SV_TRUE();
  const BLASLONG n4 = N & -4;
  const BLASLONG n2 = N & -2;
  const BLASLONG v_m1 = M & -v_size;
  const BLASLONG v_m2 = M & -v_size2;

  BLASLONG iy = 0;

  if (inc_x == 1) {
    BLASLONG j = 0;

    for (; j < n4; j += 4) {
      SV_TYPE temp_vec1 = SV_DUP(0.0);
      SV_TYPE temp_vec2 = SV_DUP(0.0);
      SV_TYPE temp_vec3 = SV_DUP(0.0);
      SV_TYPE temp_vec4 = SV_DUP(0.0);
      BLASLONG i = 0;
      
      for (; i < v_m1; i += v_size) {
        SV_TYPE a_vec1 = svld1(pg_true, a + i);
        SV_TYPE a_vec2 = svld1(pg_true, a + i + lda);
        SV_TYPE a_vec3 = svld1(pg_true, a + i + lda * 2);
        SV_TYPE a_vec4 = svld1(pg_true, a + i + lda * 3);
        SV_TYPE x_vec = svld1(pg_true, x + i);
        temp_vec1 = svmla_x(pg_true, temp_vec1, a_vec1, x_vec);
        temp_vec2 = svmla_x(pg_true, temp_vec2, a_vec2, x_vec);
        temp_vec3 = svmla_x(pg_true, temp_vec3, a_vec3, x_vec);
        temp_vec4 = svmla_x(pg_true, temp_vec4, a_vec4, x_vec);
      }

      for (; i < M; i += v_size) {
        svbool_t pg = SV_WHILE(i, M);
        SV_TYPE a_vec1 = svld1(pg, a + i);
        SV_TYPE a_vec2 = svld1(pg, a + i + lda);
        SV_TYPE a_vec3 = svld1(pg, a + i + lda * 2);
        SV_TYPE a_vec4 = svld1(pg, a + i + lda * 3);
        SV_TYPE x_vec = svld1(pg, x + i);
        temp_vec1 = svmla_x(pg, temp_vec1, a_vec1, x_vec);
        temp_vec2 = svmla_x(pg, temp_vec2, a_vec2, x_vec);
        temp_vec3 = svmla_x(pg, temp_vec3, a_vec3, x_vec);
        temp_vec4 = svmla_x(pg, temp_vec4, a_vec4, x_vec);
      }

      FLOAT temp1 = svaddv(pg_true, temp_vec1);
      FLOAT temp2 = svaddv(pg_true, temp_vec2);
      FLOAT temp3 = svaddv(pg_true, temp_vec3);
      FLOAT temp4 = svaddv(pg_true, temp_vec4);

      y[iy] += alpha * temp1;
      y[iy + inc_y] += alpha * temp2;
      y[iy + inc_y * 2] += alpha * temp3;
      y[iy + inc_y * 3] += alpha * temp4;

      iy += inc_y * 4;
      a += lda * 4;
    }

    for (; j < n2; j += 2) {
      SV_TYPE temp_vec1 = SV_DUP(0.0);
      SV_TYPE temp_vec2 = SV_DUP(0.0);
      BLASLONG i = 0;
      for (; i < v_m1; i += v_size) {
        SV_TYPE a_vec1 = svld1(pg_true, a + i);
        SV_TYPE a_vec2 = svld1(pg_true, a + lda + i);
        SV_TYPE x_vec = svld1(pg_true, x + i);
        temp_vec1 = svmla_x(pg_true, temp_vec1, a_vec1, x_vec);
        temp_vec2 = svmla_x(pg_true, temp_vec2, a_vec2, x_vec);
      }
      for (; i < M; i += v_size) {
        svbool_t pg = SV_WHILE(i, M);
        SV_TYPE a_vec1 = svld1(pg, a + i);
        SV_TYPE a_vec2 = svld1(pg, a + lda + i);
        SV_TYPE x_vec = svld1(pg, x + i);
        temp_vec1 = svmla_x(pg_true, temp_vec1, a_vec1, x_vec);
        temp_vec2 = svmla_x(pg_true, temp_vec2, a_vec2, x_vec);
      }
      FLOAT temp1 = svaddv(pg_true, temp_vec1);
      y[iy] += alpha * temp1;
      FLOAT temp2 = svaddv(pg_true, temp_vec2);
      y[iy + inc_y] += alpha * temp2;
      iy += inc_y + inc_y;
      a += lda + lda;
    }

    for (; j < N; j++) {
      SV_TYPE temp_vec = SV_DUP(0.0);
      BLASLONG i = 0;
      for (; i < v_m2; i += v_size2) {
        SV_TYPE a_vec1 = svld1(pg_true, a + i);
        SV_TYPE a_vec2 = svld1(pg_true, a + i + v_size);
        SV_TYPE x_vec1 = svld1(pg_true, x + i);
        SV_TYPE x_vec2 = svld1(pg_true, x + i + v_size);
        temp_vec = svadd_x(pg_true, temp_vec, 
          svadd_x(pg_true,
            svmul_x(pg_true, a_vec1, x_vec1),
            svmul_x(pg_true, a_vec2, x_vec2)
          )
        );
      }
      for (; i < v_m1; i += v_size) {
        SV_TYPE a_vec = svld1(pg_true, a + i);
        SV_TYPE x_vec = svld1(pg_true, x + i);
        temp_vec = svmla_x(pg_true, temp_vec, a_vec, x_vec);
      }
      for (; i < M; i += v_size) {
        svbool_t pg = SV_WHILE(i, M);
        SV_TYPE a_vec = svld1(pg, a + i);
        SV_TYPE x_vec = svld1(pg, x + i);
        temp_vec = svmla_x(pg, temp_vec, a_vec, x_vec);
      }
      FLOAT temp = svaddv(pg_true, temp_vec);
      y[iy] += alpha * temp;
      iy += inc_y;
      a += lda;
    }
    return(0);
  }

  BLASLONG j = 0;
  for (j = 0; j < N; j++) {
    FLOAT temp = 0.0;
    BLASLONG ix = 0;
    BLASLONG i = 0;
    for (i = 0; i < M; i++) {
      temp += a[i] * x[ix];
      ix += inc_x;
    }
    y[iy] += alpha * temp;
    iy += inc_y;
    a += lda;
  }

  return (0);
}
