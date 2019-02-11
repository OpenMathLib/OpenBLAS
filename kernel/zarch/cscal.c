/***************************************************************************
Copyright (c) 2013-2019, The OpenBLAS Project
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
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#include "common.h"

static void cscal_kernel_16(BLASLONG n, FLOAT *alpha, FLOAT *x) {
  __asm__("vlrepf %%v0,0(%[alpha])\n\t"
    "vlef   %%v1,4(%[alpha]),0\n\t"
    "vlef   %%v1,4(%[alpha]),2\n\t"
    "vflcsb %%v1,%%v1\n\t"
    "vlef   %%v1,4(%[alpha]),1\n\t"
    "vlef   %%v1,4(%[alpha]),3\n\t"
    "srlg %[n],%[n],4\n\t"
    "xgr   %%r1,%%r1\n\t"
    "0:\n\t"
    "pfd 2, 1024(%%r1,%[x])\n\t"
    "vl   %%v16,0(%%r1,%[x])\n\t"
    "vl   %%v17,16(%%r1,%[x])\n\t"
    "vl   %%v18,32(%%r1,%[x])\n\t"
    "vl   %%v19,48(%%r1,%[x])\n\t"
    "vl   %%v20,64(%%r1,%[x])\n\t"
    "vl   %%v21,80(%%r1,%[x])\n\t"
    "vl   %%v22,96(%%r1,%[x])\n\t"
    "vl   %%v23,112(%%r1,%[x])\n\t"
    "verllg   %%v24,%%v16,32\n\t"
    "verllg   %%v25,%%v17,32\n\t"
    "verllg   %%v26,%%v18,32\n\t"
    "verllg   %%v27,%%v19,32\n\t"
    "verllg   %%v28,%%v20,32\n\t"
    "verllg   %%v29,%%v21,32\n\t"
    "verllg   %%v30,%%v22,32\n\t"
    "verllg   %%v31,%%v23,32\n\t"
    "vfmsb %%v16,%%v16,%%v0\n\t"
    "vfmsb %%v17,%%v17,%%v0\n\t"
    "vfmsb %%v18,%%v18,%%v0\n\t"
    "vfmsb %%v19,%%v19,%%v0\n\t"
    "vfmsb %%v20,%%v20,%%v0\n\t"
    "vfmsb %%v21,%%v21,%%v0\n\t"
    "vfmsb %%v22,%%v22,%%v0\n\t"
    "vfmsb %%v23,%%v23,%%v0\n\t"
    "vfmasb %%v16,%%v24,%%v1,%%v16\n\t"
    "vfmasb %%v17,%%v25,%%v1,%%v17\n\t"
    "vfmasb %%v18,%%v26,%%v1,%%v18\n\t"
    "vfmasb %%v19,%%v27,%%v1,%%v19\n\t"
    "vfmasb %%v20,%%v28,%%v1,%%v20\n\t"
    "vfmasb %%v21,%%v29,%%v1,%%v21\n\t"
    "vfmasb %%v22,%%v30,%%v1,%%v22\n\t"
    "vfmasb %%v23,%%v31,%%v1,%%v23\n\t"
    "vst %%v16,0(%%r1,%[x])\n\t"
    "vst %%v17,16(%%r1,%[x])\n\t"
    "vst %%v18,32(%%r1,%[x])\n\t"
    "vst %%v19,48(%%r1,%[x])\n\t"
    "vst %%v20,64(%%r1,%[x])\n\t"
    "vst %%v21,80(%%r1,%[x])\n\t"
    "vst %%v22,96(%%r1,%[x])\n\t"
    "vst %%v23,112(%%r1,%[x])\n\t"
    "agfi  %%r1,128\n\t"
    "brctg %[n],0b"
    : "+m"(*(struct { FLOAT x[n * 2]; } *) x),[n] "+&r"(n)
    : [x] "a"(x), "m"(*(const struct { FLOAT x[2]; } *) alpha),
       [alpha] "a"(alpha)
    : "cc", "r1", "v0", "v1", "v16", "v17", "v18", "v19", "v20", "v21",
       "v22", "v23", "v24", "v25", "v26", "v27", "v28", "v29", "v30",
       "v31");
}

static void cscal_kernel_16_zero_r(BLASLONG n, FLOAT *alpha, FLOAT *x) {
  __asm__("vlef   %%v0,4(%[alpha]),0\n\t"
    "vlef   %%v0,4(%[alpha]),2\n\t"
    "vflcsb %%v0,%%v0\n\t"
    "vlef   %%v0,4(%[alpha]),1\n\t"
    "vlef   %%v0,4(%[alpha]),3\n\t"
    "srlg %[n],%[n],4\n\t"
    "xgr   %%r1,%%r1\n\t"
    "0:\n\t"
    "pfd 2, 1024(%%r1,%[x])\n\t"
    "vl   %%v16,0(%%r1,%[x])\n\t"
    "vl   %%v17,16(%%r1,%[x])\n\t"
    "vl   %%v18,32(%%r1,%[x])\n\t"
    "vl   %%v19,48(%%r1,%[x])\n\t"
    "vl   %%v20,64(%%r1,%[x])\n\t"
    "vl   %%v21,80(%%r1,%[x])\n\t"
    "vl   %%v22,96(%%r1,%[x])\n\t"
    "vl   %%v23,112(%%r1,%[x])\n\t"
    "verllg   %%v16,%%v16,32\n\t"
    "verllg   %%v17,%%v17,32\n\t"
    "verllg   %%v18,%%v18,32\n\t"
    "verllg   %%v19,%%v19,32\n\t"
    "verllg   %%v20,%%v20,32\n\t"
    "verllg   %%v21,%%v21,32\n\t"
    "verllg   %%v22,%%v22,32\n\t"
    "verllg   %%v23,%%v23,32\n\t"
    "vfmsb %%v16,%%v16,%%v0\n\t"
    "vfmsb %%v17,%%v17,%%v0\n\t"
    "vfmsb %%v18,%%v18,%%v0\n\t"
    "vfmsb %%v19,%%v19,%%v0\n\t"
    "vfmsb %%v20,%%v20,%%v0\n\t"
    "vfmsb %%v21,%%v21,%%v0\n\t"
    "vfmsb %%v22,%%v22,%%v0\n\t"
    "vfmsb %%v23,%%v23,%%v0\n\t"
    "vst %%v16,0(%%r1,%[x])\n\t"
    "vst %%v17,16(%%r1,%[x])\n\t"
    "vst %%v18,32(%%r1,%[x])\n\t"
    "vst %%v19,48(%%r1,%[x])\n\t"
    "vst %%v20,64(%%r1,%[x])\n\t"
    "vst %%v21,80(%%r1,%[x])\n\t"
    "vst %%v22,96(%%r1,%[x])\n\t"
    "vst %%v23,112(%%r1,%[x])\n\t"
    "agfi  %%r1,128\n\t"
    "brctg %[n],0b"
    : "+m"(*(struct { FLOAT x[n * 2]; } *) x),[n] "+&r"(n)
    : [x] "a"(x), "m"(*(const struct { FLOAT x[2]; } *) alpha),
       [alpha] "a"(alpha)
    : "cc", "r1", "v0", "v16", "v17", "v18", "v19", "v20", "v21", "v22",
       "v23");
}

static void cscal_kernel_16_zero_i(BLASLONG n, FLOAT *alpha, FLOAT *x) {
  __asm__("vlrepf %%v0,0(%[alpha])\n\t"
    "srlg %[n],%[n],4\n\t"
    "xgr   %%r1,%%r1\n\t"
    "0:\n\t"
    "pfd 2, 1024(%%r1,%[x])\n\t"
    "vl   %%v16,0(%%r1,%[x])\n\t"
    "vl   %%v17,16(%%r1,%[x])\n\t"
    "vl   %%v18,32(%%r1,%[x])\n\t"
    "vl   %%v19,48(%%r1,%[x])\n\t"
    "vl   %%v20,64(%%r1,%[x])\n\t"
    "vl   %%v21,80(%%r1,%[x])\n\t"
    "vl   %%v22,96(%%r1,%[x])\n\t"
    "vl   %%v23,112(%%r1,%[x])\n\t"
    "vfmsb %%v16,%%v16,%%v0\n\t"
    "vfmsb %%v17,%%v17,%%v0\n\t"
    "vfmsb %%v18,%%v18,%%v0\n\t"
    "vfmsb %%v19,%%v19,%%v0\n\t"
    "vfmsb %%v20,%%v20,%%v0\n\t"
    "vfmsb %%v21,%%v21,%%v0\n\t"
    "vfmsb %%v22,%%v22,%%v0\n\t"
    "vfmsb %%v23,%%v23,%%v0\n\t"
    "vst %%v16,0(%%r1,%[x])\n\t"
    "vst %%v17,16(%%r1,%[x])\n\t"
    "vst %%v18,32(%%r1,%[x])\n\t"
    "vst %%v19,48(%%r1,%[x])\n\t"
    "vst %%v20,64(%%r1,%[x])\n\t"
    "vst %%v21,80(%%r1,%[x])\n\t"
    "vst %%v22,96(%%r1,%[x])\n\t"
    "vst %%v23,112(%%r1,%[x])\n\t"
    "agfi  %%r1,128\n\t"
    "brctg %[n],0b"
    : "+m"(*(struct { FLOAT x[n * 2]; } *) x),[n] "+&r"(n)
    : [x] "a"(x), "m"(*(const struct { FLOAT x[2]; } *) alpha),
       [alpha] "a"(alpha)
    : "cc", "r1", "v0", "v16", "v17", "v18", "v19", "v20", "v21", "v22",
       "v23");
}

static void cscal_kernel_16_zero(BLASLONG n, FLOAT *x) {
  __asm__("vzero %%v0\n\t"
    "srlg %[n],%[n],4\n\t"
    "xgr   %%r1,%%r1\n\t"
    "0:\n\t"
    "pfd 2, 1024(%%r1,%[x])\n\t"
    "vst  %%v0,0(%%r1,%[x])\n\t"
    "vst  %%v0,16(%%r1,%[x])\n\t"
    "vst  %%v0,32(%%r1,%[x])\n\t"
    "vst  %%v0,48(%%r1,%[x])\n\t"
    "vst  %%v0,64(%%r1,%[x])\n\t"
    "vst  %%v0,80(%%r1,%[x])\n\t"
    "vst  %%v0,96(%%r1,%[x])\n\t"
    "vst  %%v0,112(%%r1,%[x])\n\t"
    "agfi  %%r1,128\n\t"
    "brctg %[n],0b"
    : "=m"(*(struct { FLOAT x[n * 2]; } *) x),[n] "+&r"(n)
    : [x] "a"(x)
    : "cc", "r1", "v0");
}

static void cscal_kernel_inc_8(BLASLONG n, FLOAT *alpha, FLOAT *x,
                               BLASLONG inc_x) {
  BLASLONG i;
  BLASLONG inc_x2 = 2 * inc_x;
  BLASLONG inc_x3 = inc_x2 + inc_x;
  FLOAT t0, t1, t2, t3;
  FLOAT da_r = alpha[0];
  FLOAT da_i = alpha[1];

  for (i = 0; i < n; i += 4) {
    t0 = da_r * x[0] - da_i * x[1];
    t1 = da_r * x[inc_x] - da_i * x[inc_x + 1];
    t2 = da_r * x[inc_x2] - da_i * x[inc_x2 + 1];
    t3 = da_r * x[inc_x3] - da_i * x[inc_x3 + 1];

    x[1] = da_i * x[0] + da_r * x[1];
    x[inc_x + 1] = da_i * x[inc_x] + da_r * x[inc_x + 1];
    x[inc_x2 + 1] = da_i * x[inc_x2] + da_r * x[inc_x2 + 1];
    x[inc_x3 + 1] = da_i * x[inc_x3] + da_r * x[inc_x3 + 1];

    x[0] = t0;
    x[inc_x] = t1;
    x[inc_x2] = t2;
    x[inc_x3] = t3;

    x += 4 * inc_x;
  }
}

int CNAME(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT da_r, FLOAT da_i,
          FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *dummy,
          BLASLONG dummy2) {
  BLASLONG i = 0, j = 0;
  FLOAT temp0;
  FLOAT temp1;
  FLOAT alpha[2] __attribute__ ((aligned(16)));

  if (inc_x != 1) {
    inc_x <<= 1;

    if (da_r == 0.0) {

      BLASLONG n1 = n & -2;

      if (da_i == 0.0) {

        while (j < n1) {

          x[i] = 0.0;
          x[i + 1] = 0.0;
          x[i + inc_x] = 0.0;
          x[i + 1 + inc_x] = 0.0;
          i += 2 * inc_x;
          j += 2;

        }

        while (j < n) {

          x[i] = 0.0;
          x[i + 1] = 0.0;
          i += inc_x;
          j++;

        }

      } else {

        while (j < n1) {

          temp0 = -da_i * x[i + 1];
          x[i + 1] = da_i * x[i];
          x[i] = temp0;
          temp1 = -da_i * x[i + 1 + inc_x];
          x[i + 1 + inc_x] = da_i * x[i + inc_x];
          x[i + inc_x] = temp1;
          i += 2 * inc_x;
          j += 2;

        }

        while (j < n) {

          temp0 = -da_i * x[i + 1];
          x[i + 1] = da_i * x[i];
          x[i] = temp0;
          i += inc_x;
          j++;

        }

      }

    } else {

      if (da_i == 0.0) {
        BLASLONG n1 = n & -2;

        while (j < n1) {

          temp0 = da_r * x[i];
          x[i + 1] = da_r * x[i + 1];
          x[i] = temp0;
          temp1 = da_r * x[i + inc_x];
          x[i + 1 + inc_x] = da_r * x[i + 1 + inc_x];
          x[i + inc_x] = temp1;
          i += 2 * inc_x;
          j += 2;

        }

        while (j < n) {

          temp0 = da_r * x[i];
          x[i + 1] = da_r * x[i + 1];
          x[i] = temp0;
          i += inc_x;
          j++;

        }

      } else {

        BLASLONG n1 = n & -8;
        if (n1 > 0) {
          alpha[0] = da_r;
          alpha[1] = da_i;
          cscal_kernel_inc_8(n1, alpha, x, inc_x);
          j = n1;
          i = n1 * inc_x;
        }

        while (j < n) {

          temp0 = da_r * x[i] - da_i * x[i + 1];
          x[i + 1] = da_r * x[i + 1] + da_i * x[i];
          x[i] = temp0;
          i += inc_x;
          j++;

        }

      }

    }

    return (0);
  }

  BLASLONG n1 = n & -16;
  if (n1 > 0) {

    alpha[0] = da_r;
    alpha[1] = da_i;

    if (da_r == 0.0)
      if (da_i == 0)
        cscal_kernel_16_zero(n1, x);
      else
        cscal_kernel_16_zero_r(n1, alpha, x);
    else if (da_i == 0)
      cscal_kernel_16_zero_i(n1, alpha, x);
    else
      cscal_kernel_16(n1, alpha, x);

    i = n1 << 1;
    j = n1;
  }

  if (da_r == 0.0) {

    if (da_i == 0.0) {

      while (j < n) {

        x[i] = 0.0;
        x[i + 1] = 0.0;
        i += 2;
        j++;

      }

    } else {

      while (j < n) {

        temp0 = -da_i * x[i + 1];
        x[i + 1] = da_i * x[i];
        x[i] = temp0;
        i += 2;
        j++;

      }

    }

  } else {

    if (da_i == 0.0) {

      while (j < n) {

        temp0 = da_r * x[i];
        x[i + 1] = da_r * x[i + 1];
        x[i] = temp0;
        i += 2;
        j++;

      }

    } else {

      while (j < n) {

        temp0 = da_r * x[i] - da_i * x[i + 1];
        x[i + 1] = da_r * x[i + 1] + da_i * x[i];
        x[i] = temp0;
        i += 2;
        j++;

      }

    }

  }

  return (0);
}
