/***************************************************************************
Copyright (c) 2013-2018, The OpenBLAS Project
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

static FLOAT dasum_kernel_32(BLASLONG n, FLOAT *x)
{
    FLOAT asum;

    __asm__ (
        "vzero   %%v0               \n\t"
        "vzero   %%v1               \n\t"
        "vzero   %%v2               \n\t"
        "vzero   %%v3               \n\t"
        "srlg  %%r0,%1,5            \n\t"
        "xgr %%r1,%%r1              \n\t"
        "0:                         \n\t"
        "pfd  1, 1024(%%r1,%2)      \n\t"
        "vl  %%v16, 0(%%r1,%2)      \n\t"
        "vl  %%v17, 16(%%r1,%2)     \n\t"
        "vl  %%v18, 32(%%r1,%2)     \n\t"
        "vl  %%v19, 48(%%r1,%2)     \n\t"
        "vl  %%v20, 64(%%r1,%2)     \n\t"
        "vl  %%v21, 80(%%r1,%2)     \n\t"
        "vl  %%v22, 96(%%r1,%2)     \n\t"
        "vl  %%v23, 112(%%r1,%2)    \n\t"

        "vflpdb  %%v16, %%v16       \n\t"
        "vflpdb  %%v17, %%v17       \n\t"
        "vflpdb  %%v18, %%v18       \n\t"
        "vflpdb  %%v19, %%v19       \n\t"
        "vflpdb  %%v20, %%v20       \n\t"
        "vflpdb  %%v21, %%v21       \n\t"
        "vflpdb  %%v22, %%v22       \n\t"
        "vflpdb  %%v23, %%v23       \n\t"

        "vfadb   %%v0,%%v0,%%v16    \n\t"
        "vfadb   %%v1,%%v1,%%v17    \n\t"
        "vfadb   %%v2,%%v2,%%v18    \n\t"
        "vfadb   %%v3,%%v3,%%v19    \n\t"
        "vfadb   %%v0,%%v0,%%v20    \n\t"
        "vfadb   %%v1,%%v1,%%v21    \n\t"
        "vfadb   %%v2,%%v2,%%v22    \n\t"
        "vfadb   %%v3,%%v3,%%v23    \n\t"

        "vl  %%v16, 128(%%r1,%2)    \n\t"
        "vl  %%v17, 144(%%r1,%2)    \n\t"
        "vl  %%v18, 160(%%r1,%2)    \n\t"
        "vl  %%v19, 176(%%r1,%2)    \n\t"
        "vl  %%v20, 192(%%r1,%2)    \n\t"
        "vl  %%v21, 208(%%r1,%2)    \n\t"
        "vl  %%v22, 224(%%r1,%2)    \n\t"
        "vl  %%v23, 240(%%r1,%2)    \n\t"

        "vflpdb  %%v16, %%v16       \n\t"
        "vflpdb  %%v17, %%v17       \n\t"
        "vflpdb  %%v18, %%v18       \n\t"
        "vflpdb  %%v19, %%v19       \n\t"
        "vflpdb  %%v20, %%v20       \n\t"
        "vflpdb  %%v21, %%v21       \n\t"
        "vflpdb  %%v22, %%v22       \n\t"
        "vflpdb  %%v23, %%v23       \n\t"

        "vfadb   %%v0,%%v0,%%v16    \n\t"
        "vfadb   %%v1,%%v1,%%v17    \n\t"
        "vfadb   %%v2,%%v2,%%v18    \n\t"
        "vfadb   %%v3,%%v3,%%v19    \n\t"
        "vfadb   %%v0,%%v0,%%v20    \n\t"
        "vfadb   %%v1,%%v1,%%v21    \n\t"
        "vfadb   %%v2,%%v2,%%v22    \n\t"
        "vfadb   %%v3,%%v3,%%v23    \n\t"
        
        "agfi  %%r1,256             \n\t"
        "brctg %%r0,0b              \n\t"
        "vfadb   %%v0,%%v0,%%v1     \n\t"
        "vfadb   %%v0,%%v0,%%v2     \n\t"
        "vfadb   %%v0,%%v0,%%v3     \n\t"
        "vrepg   %%v1,%%v0,1        \n\t"
        "adbr    %%f0,%%f1          \n\t"
        "ldr     %0,%%f0                "
        :"=f"(asum)
        :"r"(n),"ZR"((const FLOAT (*)[n])x)
        :"memory","cc","r0","r1","v0","v1","v2","v3","v16","v17","v18","v19","v20","v21","v22","v23"
    );

    return asum;
}

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x) {
    BLASLONG i = 0;
    BLASLONG j = 0;
    FLOAT sumf = 0.0;
    BLASLONG n1;

    if (n <= 0 || inc_x <= 0) return sumf;

    if (inc_x == 1) {

        n1 = n & -32;
               
        if (n1 > 0) {

            sumf = dasum_kernel_32(n1, x);
            i = n1;
        }

        while (i < n) {
            sumf += ABS(x[i]);
            i++;
        }

    } else {
        BLASLONG n1 = n & -4;
        register FLOAT sum1, sum2;
        sum1 = 0.0;
        sum2 = 0.0;
        while (j < n1) {

            sum1 += ABS(x[i]);
            sum2 += ABS(x[i + inc_x]);
            sum1 += ABS(x[i + 2 * inc_x]);
            sum2 += ABS(x[i + 3 * inc_x]);

            i += inc_x * 4;
            j += 4;

        }
        sumf = sum1 + sum2;
        while (j < n) {

            sumf += ABS(x[i]);
            i += inc_x;
            j++;
        }


    }
    return sumf;
}


