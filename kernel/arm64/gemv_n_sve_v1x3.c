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

#include <arm_sve.h>

#include "common.h"

#ifdef DOUBLE
#define SV_COUNT svcntd
#define SV_TYPE svfloat64_t
#define SV_TRUE svptrue_b64
#define SV_WHILE svwhilelt_b64_s64
#define SV_DUP svdup_f64
#else
#define SV_COUNT svcntw
#define SV_TYPE svfloat32_t
#define SV_TRUE svptrue_b32
#define SV_WHILE svwhilelt_b32_s64
#define SV_DUP svdup_f32
#endif

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha, FLOAT *a,
          BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y,
          FLOAT *buffer)
{
  BLASLONG i, j;
  BLASLONG ix = 0;
  BLASLONG iy;
  FLOAT *a_ptr = a;
  FLOAT temp;

  if (inc_y == 1) {
    BLASLONG width = n / 3; // Only process full 3-column blocks
    BLASLONG sve_size = SV_COUNT();
    svbool_t pg_full = SV_TRUE();
    svbool_t pg_tail = SV_WHILE(0, m % sve_size);

    FLOAT *a0_ptr = a_ptr + lda * width * 0;
    FLOAT *a1_ptr = a_ptr + lda * width * 1;
    FLOAT *a2_ptr = a_ptr + lda * width * 2;

    FLOAT *x0_ptr = x + inc_x * width * 0;
    FLOAT *x1_ptr = x + inc_x * width * 1;
    FLOAT *x2_ptr = x + inc_x * width * 2;

    for (j = 0; j < width; j++) {
      SV_TYPE temp0_vec = SV_DUP(alpha * x0_ptr[ix]);
      SV_TYPE temp1_vec = SV_DUP(alpha * x1_ptr[ix]);
      SV_TYPE temp2_vec = SV_DUP(alpha * x2_ptr[ix]);

      i = 0;
      while ((i + sve_size - 1) < m) {
        SV_TYPE y0_vec = svld1(pg_full, y + i);

        SV_TYPE a00_vec = svld1(pg_full, a0_ptr + i);
        SV_TYPE a01_vec = svld1(pg_full, a1_ptr + i);
        SV_TYPE a02_vec = svld1(pg_full, a2_ptr + i);

        y0_vec = svmla_x(pg_full, y0_vec, temp0_vec, a00_vec);
        y0_vec = svmla_x(pg_full, y0_vec, temp1_vec, a01_vec);
        y0_vec = svmla_x(pg_full, y0_vec, temp2_vec, a02_vec);

        svst1(pg_full, y + i, y0_vec);
        i += sve_size;
      }

      if (i < m) {
        SV_TYPE y0_vec = svld1(pg_tail, y + i);

        SV_TYPE a00_vec = svld1(pg_tail, a0_ptr + i);
        SV_TYPE a01_vec = svld1(pg_tail, a1_ptr + i);
        SV_TYPE a02_vec = svld1(pg_tail, a2_ptr + i);

        y0_vec = svmla_m(pg_tail, y0_vec, temp0_vec, a00_vec);
        y0_vec = svmla_m(pg_tail, y0_vec, temp1_vec, a01_vec);
        y0_vec = svmla_m(pg_tail, y0_vec, temp2_vec, a02_vec);

        svst1(pg_tail, y + i, y0_vec);
      }
      a0_ptr += lda;
      a1_ptr += lda;
      a2_ptr += lda;
      ix += inc_x;
    }
    // Handle remaining n % 3 columns
    for (j = width * 3; j < n; j++) {
      FLOAT *a_col = a + j * lda;
      temp = alpha * x[j * inc_x];
      SV_TYPE temp_vec = SV_DUP(temp);

      i = 0;
      while ((i + sve_size - 1) < m) {
        SV_TYPE y_vec = svld1(pg_full, y + i);

        SV_TYPE a_vec = svld1(pg_full, a_col + i);

        y_vec = svmla_x(pg_full, y_vec, temp_vec, a_vec);

        svst1(pg_full, y + i, y_vec);
        i += sve_size;
      }
      if (i < m) {
        SV_TYPE y_vec = svld1(pg_tail, y + i);

        SV_TYPE a_vec = svld1(pg_tail, a_col + i);

        y_vec = svmla_m(pg_tail, y_vec, temp_vec, a_vec);

        svst1(pg_tail, y + i, y_vec);
      }
    }
    return(0);
  }

  // Fallback scalar loop
  for (j = 0; j < n; j++) {
    temp = alpha * x[ix];
    iy = 0;
    for (i = 0; i < m; i++) {
      y[iy] += temp * a_ptr[i];
      iy += inc_y;
    }
    a_ptr += lda;
    ix += inc_x;
  }
  return (0);
}
