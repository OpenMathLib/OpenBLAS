/***************************************************************************
Copyright (c) 2017, The OpenBLAS Project
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

static void zaxpy_kernel_8(BLASLONG n, FLOAT *x, FLOAT *y, FLOAT *alpha)
{
    __asm__ volatile(
#if !defined(CONJ)
        "vlrepg %%v0,0(%3)              \n\t"
        "vleg   %%v1,8(%3),0            \n\t"
        "wflcdb %%v1,%%v1               \n\t"
        "vleg   %%v1,8(%3),1            \n\t"
#else        
        "vleg   %%v0,0(%3),1            \n\t"
        "vflcdb %%v0,%%v0               \n\t"
        "vleg   %%v0,0(%3),0            \n\t"
        "vlrepg %%v1,8(%3)              \n\t"
#endif
        "srlg %%r0,%0,3                 \n\t"
        "xgr  %%r1,%%r1                 \n\t"
        "0:                             \n\t"
        "pfd 1, 1024(%%r1,%1)           \n\t"
        "pfd 2, 1024(%%r1,%2)           \n\t"

        "vl   %%v16,0(%%r1,%1)          \n\t"
        "vl   %%v17,16(%%r1,%1)         \n\t"
        "vl   %%v18,32(%%r1,%1)         \n\t"
        "vl   %%v19,48(%%r1,%1)         \n\t"
        "vl   %%v20,0(%%r1,%2)          \n\t"
        "vl   %%v21,16(%%r1,%2)         \n\t"
        "vl   %%v22,32(%%r1,%2)         \n\t"
        "vl   %%v23,48(%%r1,%2)         \n\t"
        "vpdi %%v24,%%v16,%%v16,4       \n\t"
        "vpdi %%v25,%%v17,%%v17,4       \n\t"
        "vpdi %%v26,%%v18,%%v18,4       \n\t"
        "vpdi %%v27,%%v19,%%v19,4       \n\t"

        "vfmadb %%v28,%%v16,%%v0,%%v20  \n\t"
        "vfmadb %%v29,%%v17,%%v0,%%v21  \n\t"
        "vfmadb %%v30,%%v18,%%v0,%%v22  \n\t"
        "vfmadb %%v31,%%v19,%%v0,%%v23  \n\t"

        "vfmadb %%v28,%%v24,%%v1,%%v28  \n\t"
        "vfmadb %%v29,%%v25,%%v1,%%v29  \n\t"
        "vfmadb %%v30,%%v26,%%v1,%%v30  \n\t"
        "vfmadb %%v31,%%v27,%%v1,%%v31  \n\t"

        "vst %%v28,0(%%r1,%2)           \n\t"
        "vst %%v29,16(%%r1,%2)          \n\t"
        "vst %%v30,32(%%r1,%2)          \n\t"
        "vst %%v31,48(%%r1,%2)          \n\t"

        "vl   %%v16,64(%%r1,%1)         \n\t"
        "vl   %%v17,80(%%r1,%1)         \n\t"
        "vl   %%v18,96(%%r1,%1)         \n\t"
        "vl   %%v19,112(%%r1,%1)        \n\t"
        "vl   %%v20,64(%%r1,%2)         \n\t"
        "vl   %%v21,80(%%r1,%2)         \n\t"
        "vl   %%v22,96(%%r1,%2)         \n\t"
        "vl   %%v23,112(%%r1,%2)        \n\t"
        "vpdi %%v24,%%v16,%%v16,4       \n\t"
        "vpdi %%v25,%%v17,%%v17,4       \n\t"
        "vpdi %%v26,%%v18,%%v18,4       \n\t"
        "vpdi %%v27,%%v19,%%v19,4       \n\t"

        "vfmadb %%v28,%%v16,%%v0,%%v20  \n\t"
        "vfmadb %%v29,%%v17,%%v0,%%v21  \n\t"
        "vfmadb %%v30,%%v18,%%v0,%%v22  \n\t"
        "vfmadb %%v31,%%v19,%%v0,%%v23  \n\t"

        "vfmadb %%v28,%%v24,%%v1,%%v28  \n\t"
        "vfmadb %%v29,%%v25,%%v1,%%v29  \n\t"
        "vfmadb %%v30,%%v26,%%v1,%%v30  \n\t"
        "vfmadb %%v31,%%v27,%%v1,%%v31  \n\t"

        "vst %%v28,64(%%r1,%2)          \n\t"
        "vst %%v29,80(%%r1,%2)          \n\t"
        "vst %%v30,96(%%r1,%2)          \n\t"
        "vst %%v31,112(%%r1,%2)         \n\t"

        "agfi  %%r1,128                 \n\t"
        "brctg %%r0,0b                      "
        :
        :"r"(n),"ZR"((const FLOAT (*)[n * 2])x),"ZR"((FLOAT (*)[n * 2])y),"ZQ"((const FLOAT (*)[2])alpha)
        :"memory","cc","r0","r1","v0","v1","v16","v17","v18","v19","v20","v21","v22","v23","v24","v25","v26","v27","v28","v29","v30","v31"
    );
}

int CNAME(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT da_r, FLOAT da_i, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *dummy, BLASLONG dummy2) {
    BLASLONG i = 0;
    BLASLONG ix = 0, iy = 0;
    FLOAT da[2] __attribute__ ((aligned(16)));

    if (n <= 0) return (0);

    if ((inc_x == 1) && (inc_y == 1)) {

        BLASLONG n1 = n & -8;

        if (n1) {
            da[0] = da_r;
            da[1] = da_i;
            zaxpy_kernel_8(n1, x, y, da);
            ix = 2 * n1;
        }
        i = n1;
        while (i < n) {
#if !defined(CONJ)
            y[ix] += (da_r * x[ix] - da_i * x[ix + 1]);
            y[ix + 1] += (da_r * x[ix + 1] + da_i * x[ix]);
#else
            y[ix] += (da_r * x[ix] + da_i * x[ix + 1]);
            y[ix + 1] -= (da_r * x[ix + 1] - da_i * x[ix]);
#endif
            i++;
            ix += 2;

        }
        return (0);


    }

    inc_x *= 2;
    inc_y *= 2;

    while (i < n) {

#if !defined(CONJ)
        y[iy] += (da_r * x[ix] - da_i * x[ix + 1]);
        y[iy + 1] += (da_r * x[ix + 1] + da_i * x[ix]);
#else
        y[iy] += (da_r * x[ix] + da_i * x[ix + 1]);
        y[iy + 1] -= (da_r * x[ix + 1] - da_i * x[ix]);
#endif
        ix += inc_x;
        iy += inc_y;
        i++;

    }
    return (0);

}


