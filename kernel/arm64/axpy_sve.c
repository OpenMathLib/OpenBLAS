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
#define SV_TYPE svfloat64_t
#define SV_COUNT svcntd
#define SV_DUP svdup_f64
#define SV_WHILE svwhilelt_b64_s64
#define SV_TRUE svptrue_b64
#else
#define SV_TYPE svfloat32_t
#define SV_COUNT svcntw
#define SV_DUP svdup_f32
#define SV_WHILE svwhilelt_b32_s64
#define SV_TRUE svptrue_b32
#endif

int CNAME(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT da, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *dummy, BLASLONG dummy2) {
  BLASLONG i = 0;
  BLASLONG ix = 0, iy = 0;
  BLASLONG sve_size = SV_COUNT();

  if (n < 0) return (0);
  if (da == 0.0) return (0);

  if (inc_x == 1 && inc_y == 1) {
    SV_TYPE da_vec = SV_DUP(da);
    for (i = 0; i + sve_size - 1 < n; i += sve_size) {
      SV_TYPE x_vec = svld1(SV_TRUE(), &x[i]);
      SV_TYPE y_vec = svld1(SV_TRUE(), &y[i]);
      y_vec = svmla_x(SV_TRUE(), y_vec, da_vec, x_vec);
      svst1(SV_TRUE(), &y[i], y_vec);
    }

    if (i < n) {
      svbool_t pg = SV_WHILE(i, n);
      SV_TYPE x_vec = svld1(pg, &x[i]);
      SV_TYPE y_vec = svld1(pg, &y[i]);
      y_vec = svmla_x(pg, y_vec, da_vec, x_vec);
      svst1(pg, &y[i], y_vec);
    }
    return (0);
  }

  while (i < n) {
    y[iy] += da * x[ix];
    ix += inc_x;
    iy += inc_y;
    i++;
  }

  return (0);
}
