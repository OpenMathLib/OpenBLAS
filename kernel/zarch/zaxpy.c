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

static void  zaxpy_kernel_8(BLASLONG n, FLOAT *x, FLOAT *y, FLOAT da_r,FLOAT da_i) {
    __asm__ ("pfd   1, 0(%[x_tmp]) \n\t"
            "pfd    2, 0(%[y_tmp]) \n\t"
            "lgdr   %%r1,%[alpha_r]    \n\t"
            "vlvgp  %%v28,%%r1,%%r1    \n\t"
            "lgdr   %%r1,%[alpha_i]    \n\t"
            "vlvgp  %%v29,%%r1,%%r1    \n\t"
            "sllg   %[tmp],%[tmp],4    \n\t"
            "xgr    %%r1,%%r1          \n\t"         

            ".align 16 \n\t"
            "1:     \n\t"
            "pfd    1, 256(%%r1,%[x_tmp]) \n\t"
            "pfd    2, 256(%%r1,%[y_tmp]) \n\t"
            "vleg   %%v16 ,  0(%%r1,%[y_tmp]),0 \n\t"
            "vleg   %%v17 ,  8(%%r1,%[y_tmp]),0 \n\t"
            "vleg   %%v16 , 16(%%r1,%[y_tmp]),1 \n\t"
            "vleg   %%v17 , 24(%%r1,%[y_tmp]),1 \n\t"
            "vleg   %%v18 , 32(%%r1,%[y_tmp]),0 \n\t"
            "vleg   %%v19 , 40(%%r1,%[y_tmp]),0 \n\t"
            "vleg   %%v18 , 48(%%r1,%[y_tmp]),1 \n\t"
            "vleg   %%v19 , 56(%%r1,%[y_tmp]),1 \n\t"
            "vleg   %%v24 ,  0(%%r1,%[x_tmp]),0 \n\t"
            "vleg   %%v25 ,  8(%%r1,%[x_tmp]),0 \n\t"
            "vleg   %%v24 , 16(%%r1,%[x_tmp]),1 \n\t"
            "vleg   %%v25 , 24(%%r1,%[x_tmp]),1 \n\t"
            "vleg   %%v26 , 32(%%r1,%[x_tmp]),0 \n\t"
            "vleg   %%v27 , 40(%%r1,%[x_tmp]),0 \n\t"
            "vleg   %%v26 , 48(%%r1,%[x_tmp]),1 \n\t"
            "vleg   %%v27 , 56(%%r1,%[x_tmp]),1 \n\t"
#if !defined(CONJ)
            "vfmsdb %%v16,  %%v25, %%v29,%%v16  \n\t"
            "vfmadb %%v17,  %%v24, %%v29, %%v17 \n\t"
            "vfmsdb %%v18,  %%v27, %%v29, %%v18 \n\t"
            "vfmadb %%v19,  %%v26, %%v29, %%v19 \n\t"

            "vfmsdb %%v16, %%v24, %%v28 ,%%v16  \n\t"
            "vfmadb %%v17, %%v25, %%v28, %%v17  \n\t"
            "vfmsdb %%v18, %%v26, %%v28, %%v18  \n\t"
            "vfmadb %%v19, %%v27, %%v28, %%v19  \n\t"
#else
            "vfmadb %%v16, %%v25, %%v29, %%v16  \n\t"
            "vfmsdb %%v17, %%v25, %%v28, %%v17  \n\t"
            "vfmadb %%v18, %%v27, %%v29, %%v18  \n\t"
            "vfmsdb %%v19, %%v27, %%v28, %%v19  \n\t"
            "vfmadb %%v16, %%v24, %%v28, %%v16  \n\t"
            "vfmsdb %%v17, %%v24, %%v29, %%v17  \n\t"
            "vfmadb %%v18, %%v26, %%v28, %%v18  \n\t"
            "vfmsdb %%v19, %%v26, %%v29, %%v19  \n\t"

#endif 
            "vsteg %%v16 ,  0(%%r1,%[y_tmp]),0  \n\t"
            "vsteg %%v17 ,  8(%%r1,%[y_tmp]),0  \n\t"
            "vsteg %%v16 , 16(%%r1,%[y_tmp]),1  \n\t"
            "vsteg %%v17 , 24(%%r1,%[y_tmp]),1  \n\t"

            "vsteg %%v18 , 32(%%r1,%[y_tmp]),0  \n\t"
            "vsteg %%v19 , 40(%%r1,%[y_tmp]),0  \n\t"
            "vsteg %%v18 , 48(%%r1,%[y_tmp]),1  \n\t"
            "vsteg %%v19 , 56(%%r1,%[y_tmp]),1  \n\t"

            "vleg %%v20 , 64(%%r1,%[y_tmp]),0   \n\t"
            "vleg %%v21 , 72(%%r1,%[y_tmp]),0   \n\t"
            "vleg %%v20 , 80(%%r1,%[y_tmp]),1   \n\t"
            "vleg %%v21 , 88(%%r1,%[y_tmp]),1   \n\t"

            "vleg %%v22 ,  96(%%r1,%[y_tmp]),0  \n\t"
            "vleg %%v23 , 104(%%r1,%[y_tmp]),0  \n\t"
            "vleg %%v22 , 112(%%r1,%[y_tmp]),1  \n\t"
            "vleg %%v23 , 120(%%r1,%[y_tmp]),1  \n\t"

            "vleg %%v24 , 64(%%r1,%[x_tmp]),0   \n\t"
            "vleg %%v25 , 72(%%r1,%[x_tmp]),0   \n\t"
            "vleg %%v24 , 80(%%r1,%[x_tmp]),1   \n\t"
            "vleg %%v25 , 88(%%r1,%[x_tmp]),1   \n\t"

            "vleg %%v26 ,  96(%%r1,%[x_tmp]),0  \n\t"
            "vleg %%v27 , 104(%%r1,%[x_tmp]),0  \n\t"
            "vleg %%v26 , 112(%%r1,%[x_tmp]),1  \n\t"
            "vleg %%v27 , 120(%%r1,%[x_tmp]),1  \n\t"
#if !defined(CONJ)
            "vfmsdb %%v20,  %%v25, %%v29,%%v20   \n\t"
            "vfmadb %%v21,  %%v24, %%v29, %%v21  \n\t"
            "vfmsdb %%v22,  %%v27, %%v29, %%v22  \n\t"
            "vfmadb %%v23,  %%v26, %%v29, %%v23  \n\t"

            "vfmsdb %%v20, %%v24, %%v28 ,%%v20   \n\t"
            "vfmadb %%v21, %%v25, %%v28, %%v21   \n\t"
            "vfmsdb %%v22, %%v26, %%v28, %%v22   \n\t"
            "vfmadb %%v23, %%v27, %%v28, %%v23   \n\t"
#else
            "vfmadb %%v20, %%v25, %%v29, %%v20   \n\t"
            "vfmsdb %%v21, %%v25, %%v28, %%v21   \n\t"
            "vfmadb %%v22, %%v27, %%v29, %%v22   \n\t"
            "vfmsdb %%v23, %%v27, %%v28, %%v23   \n\t"
            "vfmadb %%v20, %%v24, %%v28, %%v20   \n\t"
            "vfmsdb %%v21, %%v24, %%v29, %%v21   \n\t"
            "vfmadb %%v22, %%v26, %%v28, %%v22   \n\t"
            "vfmsdb %%v23, %%v26, %%v29, %%v23   \n\t"
#endif 
            "vsteg %%v20 , 64(%%r1,%[y_tmp]),0   \n\t"
            "vsteg %%v21 , 72(%%r1,%[y_tmp]),0   \n\t"
            "vsteg %%v20 , 80(%%r1,%[y_tmp]),1   \n\t"
            "vsteg %%v21 , 88(%%r1,%[y_tmp]),1   \n\t"

            "vsteg %%v22 ,  96(%%r1,%[y_tmp]),0  \n\t"
            "vsteg %%v23 , 104(%%r1,%[y_tmp]),0  \n\t"
            "vsteg %%v22 , 112(%%r1,%[y_tmp]),1  \n\t"
            "vsteg %%v23 , 120(%%r1,%[y_tmp]),1  \n\t"

            "la     %%r1,128(%%r1) \n\t"
            "clgrjl %%r1,%[tmp],1b        \n\t"   
            : [mem_y] "+m" (*(double (*)[2*n])y),[tmp]"+&r"(n)
            : [mem_x] "m" (*(const double (*)[2*n])x), [x_tmp] "a"(x), [y_tmp] "a"(y), [alpha_r] "f"(da_r),[alpha_i] "f"(da_i)
            : "cc",  "r1","v16",
            "v17","v18","v19","v20","v21","v22","v23","v24","v25","v26","v27","v28","v29"
            );

}

int CNAME(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT da_r, FLOAT da_i, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *dummy, BLASLONG dummy2) {
    BLASLONG i = 0;
    BLASLONG ix = 0, iy = 0;

    if (n <= 0) return (0);

    if ((inc_x == 1) && (inc_y == 1)) {

        BLASLONG n1 = n & -8;

        if (n1) { 
            zaxpy_kernel_8(n1, x, y, da_r,da_i);
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


