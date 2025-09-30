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
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
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

static FLOAT dot_kernel_sve(BLASLONG n, FLOAT* x, FLOAT* y)
{
  SV_TYPE temp0_vec = SV_DUP(0.0);
  SV_TYPE temp1_vec = SV_DUP(0.0);
  SV_TYPE temp2_vec = SV_DUP(0.0);
  SV_TYPE temp3_vec = SV_DUP(0.0);
  SV_TYPE temp4_vec = SV_DUP(0.0);
  SV_TYPE temp5_vec = SV_DUP(0.0);
  SV_TYPE temp6_vec = SV_DUP(0.0);
  SV_TYPE temp7_vec = SV_DUP(0.0);

  BLASLONG i = 0;
  BLASLONG sve_size = SV_COUNT();

  while ((i + sve_size * 8 - 1) < n) {
    FLOAT *x0_ptr = x + i;
    SV_TYPE x0_vec = svld1_vnum(SV_TRUE(), x0_ptr, 0);
    SV_TYPE x1_vec = svld1_vnum(SV_TRUE(), x0_ptr, 1);
    SV_TYPE x2_vec = svld1_vnum(SV_TRUE(), x0_ptr, 2);
    SV_TYPE x3_vec = svld1_vnum(SV_TRUE(), x0_ptr, 3);
    SV_TYPE x4_vec = svld1_vnum(SV_TRUE(), x0_ptr, 4);
    SV_TYPE x5_vec = svld1_vnum(SV_TRUE(), x0_ptr, 5);
    SV_TYPE x6_vec = svld1_vnum(SV_TRUE(), x0_ptr, 6);
    SV_TYPE x7_vec = svld1_vnum(SV_TRUE(), x0_ptr, 7);

    FLOAT *y0_ptr = y + i;
    SV_TYPE y0_vec = svld1_vnum(SV_TRUE(), y0_ptr, 0);
    SV_TYPE y1_vec = svld1_vnum(SV_TRUE(), y0_ptr, 1);
    SV_TYPE y2_vec = svld1_vnum(SV_TRUE(), y0_ptr, 2);
    SV_TYPE y3_vec = svld1_vnum(SV_TRUE(), y0_ptr, 3);
    SV_TYPE y4_vec = svld1_vnum(SV_TRUE(), y0_ptr, 4);
    SV_TYPE y5_vec = svld1_vnum(SV_TRUE(), y0_ptr, 5);
    SV_TYPE y6_vec = svld1_vnum(SV_TRUE(), y0_ptr, 6);
    SV_TYPE y7_vec = svld1_vnum(SV_TRUE(), y0_ptr, 7);

    temp0_vec = svmla_x(SV_TRUE(), temp0_vec, x0_vec, y0_vec);
    temp1_vec = svmla_x(SV_TRUE(), temp1_vec, x1_vec, y1_vec);
    temp2_vec = svmla_x(SV_TRUE(), temp2_vec, x2_vec, y2_vec);
    temp3_vec = svmla_x(SV_TRUE(), temp3_vec, x3_vec, y3_vec);
    temp4_vec = svmla_x(SV_TRUE(), temp4_vec, x4_vec, y4_vec);
    temp5_vec = svmla_x(SV_TRUE(), temp5_vec, x5_vec, y5_vec);
    temp6_vec = svmla_x(SV_TRUE(), temp6_vec, x6_vec, y6_vec);
    temp7_vec = svmla_x(SV_TRUE(), temp7_vec, x7_vec, y7_vec);

    i += sve_size * 8;
  }

  if (i < n) {
    svbool_t pg0 = SV_WHILE(i + sve_size * 0, n);
    svbool_t pg1 = SV_WHILE(i + sve_size * 1, n);
    svbool_t pg2 = SV_WHILE(i + sve_size * 2, n);
    svbool_t pg3 = SV_WHILE(i + sve_size * 3, n);
    svbool_t pg4 = SV_WHILE(i + sve_size * 4, n);
    svbool_t pg5 = SV_WHILE(i + sve_size * 5, n);
    svbool_t pg6 = SV_WHILE(i + sve_size * 6, n);
    svbool_t pg7 = SV_WHILE(i + sve_size * 7, n);

    FLOAT *x0_ptr = x + i;
    SV_TYPE x0_vec = svld1_vnum(pg0, x0_ptr, 0);
    SV_TYPE x1_vec = svld1_vnum(pg1, x0_ptr, 1);
    SV_TYPE x2_vec = svld1_vnum(pg2, x0_ptr, 2);
    SV_TYPE x3_vec = svld1_vnum(pg3, x0_ptr, 3);
    SV_TYPE x4_vec = svld1_vnum(pg4, x0_ptr, 4);
    SV_TYPE x5_vec = svld1_vnum(pg5, x0_ptr, 5);
    SV_TYPE x6_vec = svld1_vnum(pg6, x0_ptr, 6);
    SV_TYPE x7_vec = svld1_vnum(pg7, x0_ptr, 7);

    FLOAT *y0_ptr = y + i;
    SV_TYPE y0_vec = svld1_vnum(pg0, y0_ptr, 0);
    SV_TYPE y1_vec = svld1_vnum(pg1, y0_ptr, 1);
    SV_TYPE y2_vec = svld1_vnum(pg2, y0_ptr, 2);
    SV_TYPE y3_vec = svld1_vnum(pg3, y0_ptr, 3);
    SV_TYPE y4_vec = svld1_vnum(pg4, y0_ptr, 4);
    SV_TYPE y5_vec = svld1_vnum(pg5, y0_ptr, 5);
    SV_TYPE y6_vec = svld1_vnum(pg6, y0_ptr, 6);
    SV_TYPE y7_vec = svld1_vnum(pg7, y0_ptr, 7);

    temp0_vec = svmla_m(pg0, temp0_vec, x0_vec, y0_vec);
    temp1_vec = svmla_m(pg1, temp1_vec, x1_vec, y1_vec);
    temp2_vec = svmla_m(pg2, temp2_vec, x2_vec, y2_vec);
    temp3_vec = svmla_m(pg3, temp3_vec, x3_vec, y3_vec);
    temp4_vec = svmla_m(pg4, temp4_vec, x4_vec, y4_vec);
    temp5_vec = svmla_m(pg5, temp5_vec, x5_vec, y5_vec);
    temp6_vec = svmla_m(pg6, temp6_vec, x6_vec, y6_vec);
    temp7_vec = svmla_m(pg7, temp7_vec, x7_vec, y7_vec);
  }

  temp0_vec = svadd_x(SV_TRUE(), temp0_vec, temp1_vec);
  temp2_vec = svadd_x(SV_TRUE(), temp2_vec, temp3_vec);
  temp4_vec = svadd_x(SV_TRUE(), temp4_vec, temp5_vec);
  temp6_vec = svadd_x(SV_TRUE(), temp6_vec, temp7_vec);
  temp0_vec = svadd_x(SV_TRUE(), temp0_vec, temp2_vec);
  temp4_vec = svadd_x(SV_TRUE(), temp4_vec, temp6_vec);
  temp0_vec = svadd_x(SV_TRUE(), temp0_vec, temp4_vec);

  return svaddv(SV_TRUE(), temp0_vec);
}
