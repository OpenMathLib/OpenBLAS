/***************************************************************************
Copyright (c) 2013-2017, The OpenBLAS Project
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
#include <math.h>

#if defined(DOUBLE)
#define ABS fabs
#else
#define ABS fabsf
#endif

static FLOAT samax_kernel_64(BLASLONG n, FLOAT *x)
{
    FLOAT amax;

    __asm__ volatile (
        "vl     %%v0,0(%2)                \n\t"
        "vflpsb %%v0,%%v0                 \n\t"
        "srlg  %%r0,%1,6                  \n\t"
        "xgr %%r1,%%r1                    \n\t"
        "0:                               \n\t"
        "pfd 1, 1024(%%r1,%2)             \n\t"

        "vl  %%v16,0(%%r1,%2)             \n\t"
        "vl  %%v17,16(%%r1,%2)            \n\t"
        "vl  %%v18,32(%%r1,%2)            \n\t"
        "vl  %%v19,48(%%r1,%2)            \n\t"
        "vl  %%v20,64(%%r1,%2)            \n\t"
        "vl  %%v21,80(%%r1,%2)            \n\t"
        "vl  %%v22,96(%%r1,%2)            \n\t"
        "vl  %%v23,112(%%r1,%2)           \n\t"
        "vflpsb  %%v16, %%v16             \n\t"
        "vflpsb  %%v17, %%v17             \n\t"
        "vflpsb  %%v18, %%v18             \n\t"
        "vflpsb  %%v19, %%v19             \n\t"
        "vflpsb  %%v20, %%v20             \n\t"
        "vflpsb  %%v21, %%v21             \n\t"
        "vflpsb  %%v22, %%v22             \n\t"
        "vflpsb  %%v23, %%v23             \n\t"
        
        "vfchsb  %%v24,%%v16,%%v17        \n\t"
        "vfchsb  %%v25,%%v18,%%v19        \n\t"
        "vfchsb  %%v26,%%v20,%%v21        \n\t"
        "vfchsb  %%v27,%%v22,%%v23        \n\t"
        "vsel    %%v24,%%v16,%%v17,%%v24  \n\t"
        "vsel    %%v25,%%v18,%%v19,%%v25  \n\t"
        "vsel    %%v26,%%v20,%%v21,%%v26  \n\t"
        "vsel    %%v27,%%v22,%%v23,%%v27  \n\t"

        "vfchsb  %%v28,%%v24,%%v25        \n\t"
        "vfchsb  %%v29,%%v26,%%v27        \n\t"
        "vsel    %%v28,%%v24,%%v25,%%v28  \n\t"
        "vsel    %%v29,%%v26,%%v27,%%v29  \n\t"

        "vfchsb  %%v30,%%v28,%%v29        \n\t"
        "vsel    %%v30,%%v28,%%v29,%%v30  \n\t"

        "vfchsb  %%v31,%%v30,%%v0         \n\t"
        "vsel    %%v0,%%v30,%%v0,%%v31    \n\t"

        "vl  %%v16,128(%%r1,%2)           \n\t"
        "vl  %%v17,144(%%r1,%2)           \n\t"
        "vl  %%v18,160(%%r1,%2)           \n\t"
        "vl  %%v19,176(%%r1,%2)           \n\t"
        "vl  %%v20,192(%%r1,%2)           \n\t"
        "vl  %%v21,208(%%r1,%2)           \n\t"
        "vl  %%v22,224(%%r1,%2)           \n\t"
        "vl  %%v23,240(%%r1,%2)           \n\t"
        "vflpsb  %%v16, %%v16             \n\t"
        "vflpsb  %%v17, %%v17             \n\t"
        "vflpsb  %%v18, %%v18             \n\t"
        "vflpsb  %%v19, %%v19             \n\t"
        "vflpsb  %%v20, %%v20             \n\t"
        "vflpsb  %%v21, %%v21             \n\t"
        "vflpsb  %%v22, %%v22             \n\t"
        "vflpsb  %%v23, %%v23             \n\t"
        
        "vfchsb  %%v24,%%v16,%%v17        \n\t"
        "vfchsb  %%v25,%%v18,%%v19        \n\t"
        "vfchsb  %%v26,%%v20,%%v21        \n\t"
        "vfchsb  %%v27,%%v22,%%v23        \n\t"
        "vsel    %%v24,%%v16,%%v17,%%v24  \n\t"
        "vsel    %%v25,%%v18,%%v19,%%v25  \n\t"
        "vsel    %%v26,%%v20,%%v21,%%v26  \n\t"
        "vsel    %%v27,%%v22,%%v23,%%v27  \n\t"

        "vfchsb  %%v28,%%v24,%%v25        \n\t"
        "vfchsb  %%v29,%%v26,%%v27        \n\t"
        "vsel    %%v28,%%v24,%%v25,%%v28  \n\t"
        "vsel    %%v29,%%v26,%%v27,%%v29  \n\t"

        "vfchsb  %%v30,%%v28,%%v29        \n\t"
        "vsel    %%v30,%%v28,%%v29,%%v30  \n\t"

        "vfchsb  %%v31,%%v30,%%v0         \n\t"
        "vsel    %%v0,%%v30,%%v0,%%v31    \n\t"

        "agfi    %%r1, 256                \n\t"
        "brctg   %%r0, 0b                 \n\t"

        "veslg   %%v16,%%v0,32            \n\t"
        "vfchsb  %%v17,%%v16,%%v0         \n\t"
        "vsel    %%v0,%%v16,%%v0,%%v17    \n\t"

        "vrepf  %%v16,%%v0,2              \n\t"
        "wfchsb %%v17,%%v16,%%v0          \n\t"
        "vsel   %%v0,%%v16,%%v0,%%v17     \n\t"
        "ler    %0,%%f0                       "
        :"=f"(amax)
        :"r"(n),"ZR"((const FLOAT (*)[n])x)
        :"memory","cc","r0","r1","v0","v16","v17","v18","v19","v20","v21","v22","v23","v24","v25","v26","v27","v28","v29","v30","v31"
    );

    return amax;
}
 
FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x) {
    BLASLONG i = 0;
    BLASLONG j = 0;
    FLOAT maxf = 0.0;

    if (n <= 0 || inc_x <= 0) return (maxf);

    if (inc_x == 1) {

        BLASLONG n1 = n & -64;
        if (n1 > 0) {

            maxf = samax_kernel_64(n1, x);

            i = n1;
        }
        else
        {
            maxf=ABS(x[0]);
            i++;
        }

        while (i < n) {
            if (ABS(x[i]) > maxf) {
                maxf = ABS(x[i]);
            }
            i++;
        }
        return (maxf);

    } else {

        maxf=ABS(x[0]);
        i += inc_x;
        j++;

        BLASLONG n1 = (n - 1) & -4;
        while (j < n1) {

            if (ABS(x[i]) > maxf) {
                maxf = ABS(x[i]);
            }
            if (ABS(x[i + inc_x]) > maxf) {
                maxf = ABS(x[i + inc_x]);
            }
            if (ABS(x[i + 2 * inc_x]) > maxf) {
                maxf = ABS(x[i + 2 * inc_x]);
            }
            if (ABS(x[i + 3 * inc_x]) > maxf) {
                maxf = ABS(x[i + 3 * inc_x]);
            }

            i += inc_x * 4;

            j += 4;

        }


        while (j < n) {
            if (ABS(x[i]) > maxf) {
                maxf = ABS(x[i]);
            }
            i += inc_x;
            j++;
        }
        return (maxf);
    }
}
