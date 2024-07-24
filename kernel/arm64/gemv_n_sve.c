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
#ifndef DOUBLE
  const BLASLONG n8 = N & -8;
#endif
  const BLASLONG n4 = N & -4;
#ifdef DOUBLE
  const BLASLONG n2 = N & -2;
#endif
  const BLASLONG v_m1 = M & -v_size;
  const BLASLONG v_m2 = M & -v_size2;

  BLASLONG ix = 0;

  if (inc_y == 1) {
    BLASLONG j = 0;
    if (inc_x == 1) {
#ifndef DOUBLE
      for (; j < n8; j += 8) {
        SV_TYPE temp_vec1 = svmul_x(pg_true, svld1rq(pg_true, &x[ix]), alpha);
        SV_TYPE temp_vec2 = svmul_x(pg_true, svld1rq(pg_true, &x[ix + 4]), alpha);

        BLASLONG i = 0;
        for (; i < v_m1; i += v_size) {
          SV_TYPE a_vec1 = svld1(pg_true, a + i);
          SV_TYPE a_vec2 = svld1(pg_true, a + i + lda);
          SV_TYPE a_vec3 = svld1(pg_true, a + i + lda * 2);
          SV_TYPE a_vec4 = svld1(pg_true, a + i + lda * 3);
          SV_TYPE a_vec5 = svld1(pg_true, a + i + lda * 4);
          SV_TYPE a_vec6 = svld1(pg_true, a + i + lda * 5);
          SV_TYPE a_vec7 = svld1(pg_true, a + i + lda * 6);
          SV_TYPE a_vec8 = svld1(pg_true, a + i + lda * 7);
          SV_TYPE y_vec = svld1(pg_true, y + i);
          y_vec = svmla_lane(y_vec, a_vec1, temp_vec1, 0);
          y_vec = svmla_lane(y_vec, a_vec2, temp_vec1, 1);
          y_vec = svmla_lane(y_vec, a_vec3, temp_vec1, 2);
          y_vec = svmla_lane(y_vec, a_vec4, temp_vec1, 3);
          y_vec = svmla_lane(y_vec, a_vec5, temp_vec2, 0);
          y_vec = svmla_lane(y_vec, a_vec6, temp_vec2, 1);
          y_vec = svmla_lane(y_vec, a_vec7, temp_vec2, 2);
          y_vec = svmla_lane(y_vec, a_vec8, temp_vec2, 3);
          svst1(pg_true, y + i, y_vec);
        }

        for (; i < M; i += v_size) {
          svbool_t pg = SV_WHILE(i, M);
          SV_TYPE a_vec1 = svld1(pg, a + i);
          SV_TYPE a_vec2 = svld1(pg, a + i + lda);
          SV_TYPE a_vec3 = svld1(pg, a + i + lda * 2);
          SV_TYPE a_vec4 = svld1(pg, a + i + lda * 3);
          SV_TYPE a_vec5 = svld1(pg, a + i + lda * 4);
          SV_TYPE a_vec6 = svld1(pg, a + i + lda * 5);
          SV_TYPE a_vec7 = svld1(pg, a + i + lda * 6);
          SV_TYPE a_vec8 = svld1(pg, a + i + lda * 7);

          SV_TYPE y_vec = svld1(pg, y + i);
          y_vec = svmla_lane(y_vec, a_vec1, temp_vec1, 0);
          y_vec = svmla_lane(y_vec, a_vec2, temp_vec1, 1);
          y_vec = svmla_lane(y_vec, a_vec3, temp_vec1, 2);
          y_vec = svmla_lane(y_vec, a_vec4, temp_vec1, 3);
          y_vec = svmla_lane(y_vec, a_vec5, temp_vec2, 0);
          y_vec = svmla_lane(y_vec, a_vec6, temp_vec2, 1);
          y_vec = svmla_lane(y_vec, a_vec7, temp_vec2, 2);
          y_vec = svmla_lane(y_vec, a_vec8, temp_vec2, 3);
          svst1(pg, y + i, y_vec);
        }

        a += lda * 8;
        ix += 8;
      }
      for (; j < n4; j += 4) {
        SV_TYPE temp_vec1 = svmul_x(pg_true, svld1rq(pg_true, &x[ix]), alpha);

        BLASLONG i = 0;
        for (; i < v_m1; i += v_size) {
          SV_TYPE a_vec1 = svld1(pg_true, a + i);
          SV_TYPE a_vec2 = svld1(pg_true, a + i + lda);
          SV_TYPE a_vec3 = svld1(pg_true, a + i + lda * 2);
          SV_TYPE a_vec4 = svld1(pg_true, a + i + lda * 3);
          SV_TYPE y_vec = svld1(pg_true, y + i);
          y_vec = svmla_lane(y_vec, a_vec1, temp_vec1, 0);
          y_vec = svmla_lane(y_vec, a_vec2, temp_vec1, 1);
          y_vec = svmla_lane(y_vec, a_vec3, temp_vec1, 2);
          y_vec = svmla_lane(y_vec, a_vec4, temp_vec1, 3);
          svst1(pg_true, y + i, y_vec);
        }

        for (; i < M; i += v_size) {
          svbool_t pg = SV_WHILE(i, M);
          SV_TYPE a_vec1 = svld1(pg, a + i);
          SV_TYPE a_vec2 = svld1(pg, a + i + lda);
          SV_TYPE a_vec3 = svld1(pg, a + i + lda * 2);
          SV_TYPE a_vec4 = svld1(pg, a + i + lda * 3);

          SV_TYPE y_vec = svld1(pg, y + i);
          y_vec = svmla_lane(y_vec, a_vec1, temp_vec1, 0);
          y_vec = svmla_lane(y_vec, a_vec2, temp_vec1, 1);
          y_vec = svmla_lane(y_vec, a_vec3, temp_vec1, 2);
          y_vec = svmla_lane(y_vec, a_vec4, temp_vec1, 3);
          svst1(pg, y + i, y_vec);
        }

        a += lda * 4;
        ix += 4;
      }
#else
      for (; j < n4; j += 4) {
        SV_TYPE temp_vec1 = svmul_x(pg_true, svld1rq(pg_true, &x[ix]), alpha);
        SV_TYPE temp_vec2 = svmul_x(pg_true, svld1rq(pg_true, &x[ix + 2]), alpha);

        BLASLONG i = 0;
        for (; i < v_m1; i += v_size) {
          SV_TYPE a_vec1 = svld1(pg_true, a + i);
          SV_TYPE a_vec2 = svld1(pg_true, a + i + lda);
          SV_TYPE a_vec3 = svld1(pg_true, a + i + lda * 2);
          SV_TYPE a_vec4 = svld1(pg_true, a + i + lda * 3);
          SV_TYPE y_vec = svld1(pg_true, y + i);
          y_vec = svmla_lane(y_vec, a_vec1, temp_vec1, 0);
          y_vec = svmla_lane(y_vec, a_vec2, temp_vec1, 1);
          y_vec = svmla_lane(y_vec, a_vec3, temp_vec2, 0);
          y_vec = svmla_lane(y_vec, a_vec4, temp_vec2, 1);
          svst1(pg_true, y + i, y_vec);
        }
        for (; i < M; i += v_size) {
          svbool_t pg = SV_WHILE(i, M);
          SV_TYPE a_vec1 = svld1(pg, a + i);
          SV_TYPE a_vec2 = svld1(pg, a + i + lda);
          SV_TYPE a_vec3 = svld1(pg, a + i + lda * 2);
          SV_TYPE a_vec4 = svld1(pg, a + i + lda * 3);
          SV_TYPE y_vec = svld1(pg, y + i);
          y_vec = svmla_lane(y_vec, a_vec1, temp_vec1, 0);
          y_vec = svmla_lane(y_vec, a_vec2, temp_vec1, 1);
          y_vec = svmla_lane(y_vec, a_vec3, temp_vec2, 0);
          y_vec = svmla_lane(y_vec, a_vec4, temp_vec2, 1);
          svst1(pg, y + i, y_vec);
        }

        a += lda * 4;
        ix += 4;
      }
      for (; j < n2; j += 2) {
        SV_TYPE temp_vec1 = svmul_x(pg_true, svld1rq(pg_true, &x[ix]), alpha);

        BLASLONG i = 0;
        for (; i < v_m1; i += v_size) {
          SV_TYPE a_vec1 = svld1(pg_true, a + i);
          SV_TYPE a_vec2 = svld1(pg_true, a + i + lda);
          SV_TYPE y_vec = svld1(pg_true, y + i);
          y_vec = svmla_lane(y_vec, a_vec1, temp_vec1, 0);
          y_vec = svmla_lane(y_vec, a_vec2, temp_vec1, 1);
          svst1(pg_true, y + i, y_vec);
        }
        for (; i < M; i += v_size) {
          svbool_t pg = SV_WHILE(i, M);
          SV_TYPE a_vec1 = svld1(pg, a + i);
          SV_TYPE a_vec2 = svld1(pg, a + i + lda);
          SV_TYPE y_vec = svld1(pg, y + i);
          y_vec = svmla_lane(y_vec, a_vec1, temp_vec1, 0);
          y_vec = svmla_lane(y_vec, a_vec2, temp_vec1, 1);
          svst1(pg, y + i, y_vec);
        }

        a += lda * 2;
        ix += 2;
      }
#endif
    }

    for (; j < N; j++) {
      SV_TYPE temp_vec1 = SV_DUP(alpha * x[ix]);
      SV_TYPE temp_vec2 = temp_vec1;
      
      BLASLONG i = 0;
      for (; i < v_m1; i += v_size) {
        SV_TYPE a_vec = svld1(pg_true, a + i);
        SV_TYPE y_vec = svld1(pg_true, y + i);
        y_vec = svmla_x(pg_true, y_vec, temp_vec1, a_vec);
        svst1(pg_true, y + i, y_vec);
      }
      for (; i < M; i += v_size) {
        svbool_t pg = SV_WHILE(i, M);
        SV_TYPE a_vec = svld1(pg, a + i);
        SV_TYPE y_vec = svld1(pg, y + i);
        y_vec = svmla_x(pg, y_vec, temp_vec1, a_vec);
        svst1(pg, y + i, y_vec);
      }
      a += lda;
      ix += inc_x;
    }
    return(0);
  }

  for (BLASLONG j = 0; j < N; j++) {
    FLOAT temp = alpha * x[ix];
    BLASLONG iy = 0;
    for (BLASLONG i = 0; i < M; i++) {
      y[iy] += temp * a[i];
      iy += inc_y;
    }
    a += lda;
    ix += inc_x;
  }
  return (0);
}
