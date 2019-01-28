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
#include <math.h>
 
#if defined(DOUBLE)
#define ABS fabs
#else
#define ABS fabsf
#endif
#define CABS1(x,i) (ABS(x[i]) + ABS(x[i + 1]))

static BLASLONG icamin_kernel_32(BLASLONG n, FLOAT *x, FLOAT *amin)
{ 
    BLASLONG iamin;

    __asm__ volatile (
        "vlef   %%v0,0(%3),0             \n\t"
        "vlef   %%v1,4(%3),0             \n\t"
        "vlef   %%v0,8(%3),1             \n\t"
        "vlef   %%v1,12(%3),1            \n\t"
        "vlef   %%v0,16(%3),2            \n\t"
        "vlef   %%v1,20(%3),2            \n\t"
        "vlef   %%v0,24(%3),3            \n\t"
        "vlef   %%v1,28(%3),3            \n\t"
        "vflpsb %%v0,%%v0                \n\t"
        "vflpsb %%v1,%%v1                \n\t"
        "vfasb  %%v0,%%v0,%%v1           \n\t"
        "vleig  %%v1,0,0                 \n\t"
        "vleig  %%v1,2,1                 \n\t"
        "vleig  %%v2,1,0                 \n\t"
        "vleig  %%v2,3,1                 \n\t"
        "vrepig %%v3,16                  \n\t"
        "vzero  %%v4                     \n\t"
        "vleif  %%v24,0,0                \n\t"
        "vleif  %%v24,1,1                \n\t"
        "vleif  %%v24,2,2                \n\t"
        "vleif  %%v24,3,3                \n\t"
        "vleif  %%v25,4,0                \n\t"
        "vleif  %%v25,5,1                \n\t"
        "vleif  %%v25,6,2                \n\t"
        "vleif  %%v25,7,3                \n\t"
        "vleif  %%v26,8,0                \n\t"
        "vleif  %%v26,9,1                \n\t"
        "vleif  %%v26,10,2               \n\t"
        "vleif  %%v26,11,3               \n\t"
        "vleif  %%v27,12,0               \n\t"
        "vleif  %%v27,13,1               \n\t"
        "vleif  %%v27,14,2               \n\t"
        "vleif  %%v27,15,3               \n\t"
        "srlg  %%r0,%2,5                 \n\t"
        "xgr %%r1,%%r1                   \n\t"
        "0:                              \n\t"
        "pfd 1, 1024(%%r1,%3)            \n\t"

        "vlef  %%v16,0(%%r1,%3),0        \n\t"
        "vlef  %%v17,4(%%r1,%3),0        \n\t"
        "vlef  %%v16,8(%%r1,%3),1        \n\t"
        "vlef  %%v17,12(%%r1,%3),1       \n\t"
        "vlef  %%v16,16(%%r1,%3),2       \n\t"
        "vlef  %%v17,20(%%r1,%3),2       \n\t"
        "vlef  %%v16,24(%%r1,%3),3       \n\t"
        "vlef  %%v17,28(%%r1,%3),3       \n\t"

        "vlef  %%v18,32(%%r1,%3),0       \n\t"
        "vlef  %%v19,36(%%r1,%3),0       \n\t"
        "vlef  %%v18,40(%%r1,%3),1       \n\t"
        "vlef  %%v19,44(%%r1,%3),1       \n\t"
        "vlef  %%v18,48(%%r1,%3),2       \n\t"
        "vlef  %%v19,52(%%r1,%3),2       \n\t"
        "vlef  %%v18,56(%%r1,%3),3       \n\t"
        "vlef  %%v19,30(%%r1,%3),3       \n\t"

        "vlef  %%v20,64(%%r1,%3),0       \n\t"
        "vlef  %%v21,68(%%r1,%3),0       \n\t"
        "vlef  %%v20,72(%%r1,%3),1       \n\t"
        "vlef  %%v21,76(%%r1,%3),1       \n\t"
        "vlef  %%v20,80(%%r1,%3),2       \n\t"
        "vlef  %%v21,84(%%r1,%3),2       \n\t"
        "vlef  %%v20,88(%%r1,%3),3       \n\t"
        "vlef  %%v21,92(%%r1,%3),3       \n\t"

        "vlef  %%v22,96(%%r1,%3),0       \n\t"
        "vlef  %%v23,100(%%r1,%3),0      \n\t"
        "vlef  %%v22,104(%%r1,%3),1      \n\t"
        "vlef  %%v23,108(%%r1,%3),1      \n\t"
        "vlef  %%v22,112(%%r1,%3),2      \n\t"
        "vlef  %%v23,116(%%r1,%3),2      \n\t"
        "vlef  %%v22,120(%%r1,%3),3      \n\t"
        "vlef  %%v23,124(%%r1,%3),3      \n\t"

        "vflpsb  %%v16, %%v16            \n\t"
        "vflpsb  %%v17, %%v17            \n\t"
        "vflpsb  %%v18, %%v18            \n\t"
        "vflpsb  %%v19, %%v19            \n\t"
        "vflpsb  %%v20, %%v20            \n\t"
        "vflpsb  %%v21, %%v21            \n\t"
        "vflpsb  %%v22, %%v22            \n\t"
        "vflpsb  %%v23, %%v23            \n\t"
        "vfasb %%v16,%%v16,%%v17         \n\t"
        "vfasb %%v17,%%v18,%%v19         \n\t"
        "vfasb %%v18,%%v20,%%v21         \n\t"
        "vfasb %%v19,%%v22,%%v23         \n\t"
        
        "vfchesb  %%v5,%%v17,%%v16       \n\t"
        "vfchesb  %%v6,%%v19,%%v18       \n\t"
        "vsel    %%v16,%%v16,%%v17,%%v5  \n\t"
        "vsel    %%v5,%%v24,%%v25,%%v5   \n\t"
        "vsel    %%v17,%%v18,%%v19,%%v6  \n\t"
        "vsel    %%v6,%%v26,%%v27,%%v6   \n\t"

        "vfchesb  %%v18,%%v17,%%v16      \n\t"
        "vsel    %%v16,%%v16,%%v17,%%v18 \n\t"
        "vsel    %%v5,%%v5,%%v6,%%v18    \n\t"
        "vsegf   %%v6,%%v5               \n\t"
        "vesrlg  %%v5,%%v5,32            \n\t"
        "vag     %%v5,%%v5,%%v4          \n\t"
        "vag     %%v6,%%v6,%%v4          \n\t"

        "vfchesb  %%v7,%%v16,%%v0        \n\t"
        "vsel    %%v0,%%v0,%%v16,%%v7    \n\t"
        "vsegf   %%v8,%%v7               \n\t"
        "vesrlg  %%v7,%%v7,32            \n\t"
        "vsegf   %%v7,%%v7               \n\t"
        "vsel    %%v1,%%v1,%%v5,%%v7     \n\t"
        "vsel    %%v2,%%v2,%%v6,%%v8     \n\t"
        "vag     %%v4,%%v4,%%v3          \n\t"

        "vlef  %%v16,128(%%r1,%3),0      \n\t"
        "vlef  %%v17,132(%%r1,%3),0      \n\t"
        "vlef  %%v16,136(%%r1,%3),1      \n\t"
        "vlef  %%v17,140(%%r1,%3),1      \n\t"
        "vlef  %%v16,144(%%r1,%3),2      \n\t"
        "vlef  %%v17,148(%%r1,%3),2      \n\t"
        "vlef  %%v16,152(%%r1,%3),3      \n\t"
        "vlef  %%v17,156(%%r1,%3),3      \n\t"

        "vlef  %%v18,160(%%r1,%3),0      \n\t"
        "vlef  %%v19,164(%%r1,%3),0      \n\t"
        "vlef  %%v18,168(%%r1,%3),1      \n\t"
        "vlef  %%v19,172(%%r1,%3),1      \n\t"
        "vlef  %%v18,176(%%r1,%3),2      \n\t"
        "vlef  %%v19,180(%%r1,%3),2      \n\t"
        "vlef  %%v18,184(%%r1,%3),3      \n\t"
        "vlef  %%v19,188(%%r1,%3),3      \n\t"

        "vlef  %%v20,192(%%r1,%3),0      \n\t"
        "vlef  %%v21,196(%%r1,%3),0      \n\t"
        "vlef  %%v20,200(%%r1,%3),1      \n\t"
        "vlef  %%v21,204(%%r1,%3),1      \n\t"
        "vlef  %%v20,208(%%r1,%3),2      \n\t"
        "vlef  %%v21,212(%%r1,%3),2      \n\t"
        "vlef  %%v20,216(%%r1,%3),3      \n\t"
        "vlef  %%v21,220(%%r1,%3),3      \n\t"

        "vlef  %%v22,224(%%r1,%3),0      \n\t"
        "vlef  %%v23,228(%%r1,%3),0      \n\t"
        "vlef  %%v22,232(%%r1,%3),1      \n\t"
        "vlef  %%v23,236(%%r1,%3),1      \n\t"
        "vlef  %%v22,240(%%r1,%3),2      \n\t"
        "vlef  %%v23,244(%%r1,%3),2      \n\t"
        "vlef  %%v22,248(%%r1,%3),3      \n\t"
        "vlef  %%v23,252(%%r1,%3),3      \n\t"

        "vflpsb  %%v16, %%v16            \n\t"
        "vflpsb  %%v17, %%v17            \n\t"
        "vflpsb  %%v18, %%v18            \n\t"
        "vflpsb  %%v19, %%v19            \n\t"
        "vflpsb  %%v20, %%v20            \n\t"
        "vflpsb  %%v21, %%v21            \n\t"
        "vflpsb  %%v22, %%v22            \n\t"
        "vflpsb  %%v23, %%v23            \n\t"
        "vfasb %%v16,%%v16,%%v17         \n\t"
        "vfasb %%v17,%%v18,%%v19         \n\t"
        "vfasb %%v18,%%v20,%%v21         \n\t"
        "vfasb %%v19,%%v22,%%v23         \n\t"
        
        "vfchesb  %%v5,%%v17,%%v16       \n\t"
        "vfchesb  %%v6,%%v19,%%v18       \n\t"
        "vsel    %%v16,%%v16,%%v17,%%v5  \n\t"
        "vsel    %%v5,%%v24,%%v25,%%v5   \n\t"
        "vsel    %%v17,%%v18,%%v19,%%v6  \n\t"
        "vsel    %%v6,%%v26,%%v27,%%v6   \n\t"

        "vfchesb  %%v18,%%v17,%%v16      \n\t"
        "vsel    %%v16,%%v16,%%v17,%%v18 \n\t"
        "vsel    %%v5,%%v5,%%v6,%%v18    \n\t"
        "vsegf   %%v6,%%v5               \n\t"
        "vesrlg  %%v5,%%v5,32            \n\t"
        "vag     %%v5,%%v5,%%v4          \n\t"
        "vag     %%v6,%%v6,%%v4          \n\t"

        "vfchesb  %%v7,%%v16,%%v0        \n\t"
        "vsel    %%v0,%%v0,%%v16,%%v7    \n\t"
        "vsegf   %%v8,%%v7               \n\t"
        "vesrlg  %%v7,%%v7,32            \n\t"
        "vsegf   %%v7,%%v7               \n\t"
        "vsel    %%v1,%%v1,%%v5,%%v7     \n\t"
        "vsel    %%v2,%%v2,%%v6,%%v8     \n\t"
        "vag     %%v4,%%v4,%%v3          \n\t"

        "agfi    %%r1, 256               \n\t"
        "brctg   %%r0, 0b                \n\t"

        "veslg   %%v3,%%v0,32            \n\t"
        "vfchsb  %%v4,%%v3,%%v0          \n\t"
        "vchlg   %%v5,%%v2,%%v1          \n\t"
        "vfcesb  %%v6,%%v0,%%v3          \n\t"
        "vn      %%v5,%%v5,%%v6          \n\t"
        "vo      %%v4,%%v4,%%v5          \n\t"
        "vsel    %%v0,%%v0,%%v3,%%v4     \n\t"
        "vesrlg  %%v4,%%v4,32            \n\t"
        "vsegf   %%v4,%%v4               \n\t"
        "vsel    %%v1,%%v1,%%v2,%%v4     \n\t"

        "vrepf  %%v2,%%v0,2              \n\t"
        "vrepg  %%v3,%%v1,1              \n\t"
        "wfcsb  %%v2,%%v0                \n\t"
        "jne 1f                          \n\t"
        "vstef  %%v0,%1,0                \n\t"
        "vmnlg  %%v0,%%v1,%%v3           \n\t"
        "vlgvg  %0,%%v0,0                \n\t"
        "j 2f                            \n\t"
        "1:                              \n\t"
        "wfchsb %%v4,%%v0,%%v2           \n\t"
        "vsel   %%v1,%%v3,%%v1,%%v4      \n\t"
        "vsel   %%v0,%%v2,%%v0,%%v4      \n\t"
        "ste    %%f0,%1                  \n\t"
        "vlgvg  %0,%%v1,0                \n\t"
        "2:                              \n\t"
        "nop                                 "
        :"=r"(iamin),"=m"(*amin)
        :"r"(n),"ZR"((const FLOAT (*)[n * 2])x)
        :"memory","cc","r0","r1","v0","v1","v2","v3","v4","v5","v6","v7","v8","v16","v17","v18","v19","v20","v21","v22","v23","v24","v25","v26","v27"
    );

    return iamin;
}

BLASLONG CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
    BLASLONG i = 0;
    BLASLONG ix = 0;
    FLOAT minf = 0;
    BLASLONG min = 0;
    BLASLONG inc_x2;

    if (n <= 0 || inc_x <= 0) return(min);
     
    if (inc_x == 1) {

      BLASLONG n1 = n & -32;
      if (n1 > 0) {

            min = icamin_kernel_32(n1, x, &minf);

            i = n1;
      }
      else
        {
            minf = CABS1(x,0);
            ix += 2;
            i++;
        }

      while(i < n)
    {
        if( CABS1(x,ix) < minf )
        {
            min = i;
            minf = CABS1(x,ix);
        }
        ix += 2;
        i++;
    }
        return (min + 1);

    } else {
 
    min = 0;
    minf = CABS1(x,0);
    inc_x2 = 2 * inc_x;
    ix += inc_x2;
    i++;

    while(i < n)
    {
        if( CABS1(x,ix) < minf )
        {
            min = i;
            minf = CABS1(x,ix);
        }
        ix += inc_x2;
        i++;
    }
        return (min + 1);
    }
}


