/***************************************************************************
Copyright (c) 2019, The OpenBLAS Project
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
#if !defined(__VEC__) || !defined(__ALTIVEC__)
#include "../arm/gemv_t.c"

#else

#include "common.h"

#define NBMAX 2048

#include <altivec.h> 
 
static void sgemv_kernel_4x8(BLASLONG n, BLASLONG lda, FLOAT *ap, FLOAT *x, FLOAT *y, FLOAT alpha) {
    BLASLONG i;  
    FLOAT *a0, *a1, *a2, *a3, *a4, *a5, *a6, *a7;
    __vector float *va0, *va1, *va2, *va3, *va4, *va5, *va6, *va7, *v_x;
    register __vector float temp0 = {0,0,0,0};
    register __vector float temp1 = {0,0,0,0};
    register __vector float temp2 = {0,0,0,0};
    register __vector float temp3 = {0,0,0,0};
    register __vector float temp4 = {0,0,0,0};
    register __vector float temp5 = {0,0,0,0};
    register __vector float temp6 = {0,0,0,0};
    register __vector float temp7 = {0,0,0,0};

    a0 = ap;
    a1 = ap + lda;
    a2 = a1 + lda;
    a3 = a2 + lda;
    a4 = a3 + lda;
    a5 = a4 + lda;
    a6 = a5 + lda;
    a7 = a6 + lda;
    va0 = (__vector float*) a0;
    va1 = (__vector float*) a1;
    va2 = (__vector float*) a2;
    va3 = (__vector float*) a3;
    va4 = (__vector float*) a4;
    va5 = (__vector float*) a5;
    va6 = (__vector float*) a6;
    va7 = (__vector float*) a7;
    v_x = (__vector float*) x;
 
   
        for (i = 0; i < n/4; i ++) {
            temp0 += v_x[i] * va0[i];
            temp1 += v_x[i] * va1[i];
            temp2 += v_x[i] * va2[i];
            temp3 += v_x[i] * va3[i];
            temp4 += v_x[i] * va4[i];
            temp5 += v_x[i] * va5[i];
            temp6 += v_x[i] * va6[i];
            temp7 += v_x[i] * va7[i]; 
        }
    
  
    y[0] += alpha * (temp0[0] + temp0[1]+temp0[2] + temp0[3]);
    y[1] += alpha * (temp1[0] + temp1[1]+temp1[2] + temp1[3]);
    y[2] += alpha * (temp2[0] + temp2[1]+temp2[2] + temp2[3]);
    y[3] += alpha * (temp3[0] + temp3[1]+temp3[2] + temp3[3]);

    y[4] += alpha * (temp4[0] + temp4[1]+temp4[2] + temp4[3]);
    y[5] += alpha * (temp5[0] + temp5[1]+temp5[2] + temp5[3]);
    y[6] += alpha * (temp6[0] + temp6[1]+temp6[2] + temp6[3]);
    y[7] += alpha * (temp7[0] + temp7[1]+temp7[2] + temp7[3]);

}
 

static void sgemv_kernel_4x4(BLASLONG n, BLASLONG lda, FLOAT *ap, FLOAT *x, FLOAT *y, FLOAT alpha) {
    BLASLONG i = 0;
    FLOAT *a0, *a1, *a2, *a3;
    a0 = ap;
    a1 = ap + lda;
    a2 = a1 + lda;
    a3 = a2 + lda;
    __vector float* va0 = (__vector float*) a0;
    __vector float* va1 = (__vector float*) a1;
    __vector float* va2 = (__vector float*) a2;
    __vector float* va3 = (__vector float*) a3;
    __vector float* v_x = (__vector float*) x;
    register __vector float temp0 = {0,0,0,0};
    register __vector float temp1 = {0,0,0,0};
    register __vector float temp2 = {0,0,0,0};
    register __vector float temp3 = {0,0,0,0}; 

    for (i = 0; i < n / 4; i ++) {
        temp0 += v_x[i] * va0[i];
        temp1 += v_x[i] * va1[i];
        temp2 += v_x[i] * va2[i];
        temp3 += v_x[i] * va3[i]; 
    }
 
    y[0] += alpha * (temp0[0] + temp0[1]+temp0[2] + temp0[3]);
    y[1] += alpha * (temp1[0] + temp1[1]+temp1[2] + temp1[3]);
    y[2] += alpha * (temp2[0] + temp2[1]+temp2[2] + temp2[3]);
    y[3] += alpha * (temp3[0] + temp3[1]+temp3[2] + temp3[3]);

}
 

static void sgemv_kernel_4x2(BLASLONG n, BLASLONG lda, FLOAT *ap, FLOAT *x, FLOAT *y, FLOAT alpha, BLASLONG inc_y) {

    BLASLONG i;
    FLOAT *a0, *a1;
    a0 = ap;
    a1 = ap + lda;
    __vector float* va0 = (__vector float*) a0;
    __vector float* va1 = (__vector float*) a1;
    __vector float* v_x = (__vector float*) x;
    __vector float temp0 = {0,0,0,0};
    __vector float temp1 = {0,0,0,0};
    for (i = 0; i < n / 4; i ++) {
        temp0 += v_x[i] * va0[i];
        temp1 += v_x[i] * va1[i];
    }



    y[0] += alpha * (temp0[0] + temp0[1]+temp0[2] + temp0[3]);
    y[inc_y] += alpha * (temp1[0] + temp1[1]+temp1[2] + temp1[3]); 
}

static void sgemv_kernel_4x1(BLASLONG n, FLOAT *ap, FLOAT *x, FLOAT *y, FLOAT alpha) {

    BLASLONG i;
    FLOAT *a0;
    a0 = ap;
    __vector float* va0 = (__vector float*) a0;
    __vector float* v_x = (__vector float*) x;
    __vector float temp0 = {0,0,0,0};
    for (i = 0; i < n / 4; i ++) {
        temp0 += v_x[i] * va0[i] ;
    }

    y[0] += alpha * (temp0[0] + temp0[1]+temp0[2] + temp0[3]);

}

static void copy_x(BLASLONG n, FLOAT *src, FLOAT *dest, BLASLONG inc_src) {
    BLASLONG i;
    for (i = 0; i < n; i++) {
        *dest++ = *src;
        src += inc_src;
    }
}

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *buffer) {
    BLASLONG i;
    BLASLONG j;
    FLOAT *a_ptr;
    FLOAT *x_ptr;
    FLOAT *y_ptr;

    BLASLONG n1;
    BLASLONG m1;
    BLASLONG m2;
    BLASLONG m3;
    BLASLONG n2; 
    FLOAT ybuffer[8] __attribute__((aligned(16)));
    FLOAT *xbuffer; 
    if (m < 1) return (0);
    if (n < 1) return (0);

    xbuffer = buffer;

    n1 = n >> 3;
    n2 = n & 7;

    m3 = m & 3;
    m1 = m - m3;
    m2 = (m & (NBMAX - 1)) - m3;

    BLASLONG NB = NBMAX;

    while (NB == NBMAX) {

        m1 -= NB;
        if (m1 < 0) {
            if (m2 == 0) break;
            NB = m2;
        }

        y_ptr = y;
        a_ptr = a;
        x_ptr = x;

        if (inc_x != 1)
            copy_x(NB, x_ptr, xbuffer, inc_x);
        else
            xbuffer = x_ptr;

        BLASLONG lda8 = lda << 3;


        if (inc_y == 1) {

            for (i = 0; i < n1; i++) {
                 
                sgemv_kernel_4x8(NB, lda, a_ptr, xbuffer, y_ptr, alpha);
 
                y_ptr += 8;
                a_ptr += lda8;
        
            }

        } else {
                   
            for (i = 0; i < n1; i++) {
                ybuffer[0] = 0;
                ybuffer[1] = 0;
                ybuffer[2] = 0;
                ybuffer[3] = 0;
                ybuffer[4] = 0;
                ybuffer[5] = 0;
                ybuffer[6] = 0;
                ybuffer[7] = 0;
                sgemv_kernel_4x8(NB, lda, a_ptr, xbuffer, ybuffer, alpha);

 

                *y_ptr += ybuffer[0];
                y_ptr += inc_y;
                *y_ptr += ybuffer[1];
                y_ptr += inc_y;
                *y_ptr += ybuffer[2];
                y_ptr += inc_y;
                *y_ptr += ybuffer[3];
                y_ptr += inc_y;

                *y_ptr += ybuffer[4];
                y_ptr += inc_y;
                *y_ptr += ybuffer[5];
                y_ptr += inc_y;
                *y_ptr += ybuffer[6];
                y_ptr += inc_y;
                *y_ptr += ybuffer[7];
                y_ptr += inc_y;

                a_ptr += lda8;
            }

        }


        if (n2 & 4) {
            ybuffer[0] = 0;
            ybuffer[1] = 0;
            ybuffer[2] = 0;
            ybuffer[3] = 0;
            sgemv_kernel_4x4(NB, lda, a_ptr, xbuffer, ybuffer, alpha);

            a_ptr += lda<<2;

            *y_ptr += ybuffer[0];
            y_ptr += inc_y;
            *y_ptr += ybuffer[1];
            y_ptr += inc_y;
            *y_ptr += ybuffer[2];
            y_ptr += inc_y;
            *y_ptr += ybuffer[3];
            y_ptr += inc_y;
        }

        if (n2 & 2) {
            sgemv_kernel_4x2(NB, lda, a_ptr, xbuffer, y_ptr, alpha, inc_y);
            a_ptr += lda << 1;
            y_ptr += 2 * inc_y;

        }

        if (n2 & 1) {
            sgemv_kernel_4x1(NB, a_ptr, xbuffer, y_ptr, alpha);
            a_ptr += lda;
            y_ptr += inc_y;

        }

        a += NB;
        x += NB * inc_x;


    }

    if (m3 == 0) return (0);

    x_ptr = x;
    a_ptr = a;
    if (m3 == 3) {
        FLOAT xtemp0 = *x_ptr * alpha;
        x_ptr += inc_x;
        FLOAT xtemp1 = *x_ptr * alpha;
        x_ptr += inc_x;
        FLOAT xtemp2 = *x_ptr * alpha;

        FLOAT *aj = a_ptr;
        y_ptr = y;

        if (lda == 3 && inc_y == 1) {

            for (j = 0; j < (n & -4); j += 4) {

                y_ptr[j] += aj[0] * xtemp0 + aj[1] * xtemp1 + aj[2] * xtemp2;
                y_ptr[j + 1] += aj[3] * xtemp0 + aj[4] * xtemp1 + aj[5] * xtemp2;
                y_ptr[j + 2] += aj[6] * xtemp0 + aj[7] * xtemp1 + aj[8] * xtemp2;
                y_ptr[j + 3] += aj[9] * xtemp0 + aj[10] * xtemp1 + aj[11] * xtemp2;
                aj += 12;
            }

            for (; j < n; j++) {
                y_ptr[j] += aj[0] * xtemp0 + aj[1] * xtemp1 + aj[2] * xtemp2;
                aj += 3;
            }

        } else {

            if (inc_y == 1) {

                BLASLONG register lda2 = lda << 1;
                BLASLONG register lda4 = lda << 2;
                BLASLONG register lda3 = lda2 + lda;

                for (j = 0; j < (n & -4); j += 4) {

                    y_ptr[j] += *aj * xtemp0 + *(aj + 1) * xtemp1 + *(aj + 2) * xtemp2;
                    y_ptr[j + 1] += *(aj + lda) * xtemp0 + *(aj + lda + 1) * xtemp1 + *(aj + lda + 2) * xtemp2;
                    y_ptr[j + 2] += *(aj + lda2) * xtemp0 + *(aj + lda2 + 1) * xtemp1 + *(aj + lda2 + 2) * xtemp2;
                    y_ptr[j + 3] += *(aj + lda3) * xtemp0 + *(aj + lda3 + 1) * xtemp1 + *(aj + lda3 + 2) * xtemp2;
                    aj += lda4;
                }

                for (; j < n; j++) {

                    y_ptr[j] += *aj * xtemp0 + *(aj + 1) * xtemp1 + *(aj + 2) * xtemp2;
                    aj += lda;
                }

            } else {

                for (j = 0; j < n; j++) {
                    *y_ptr += *aj * xtemp0 + *(aj + 1) * xtemp1 + *(aj + 2) * xtemp2;
                    y_ptr += inc_y;
                    aj += lda;
                }

            }

        }
        return (0);
    }

    if (m3 == 2) {
        FLOAT xtemp0 = *x_ptr * alpha;
        x_ptr += inc_x;
        FLOAT xtemp1 = *x_ptr * alpha;

        FLOAT *aj = a_ptr;
        y_ptr = y;

        if (lda == 2 && inc_y == 1) {

            for (j = 0; j < (n & -4); j += 4) {
                y_ptr[j] += aj[0] * xtemp0 + aj[1] * xtemp1;
                y_ptr[j + 1] += aj[2] * xtemp0 + aj[3] * xtemp1;
                y_ptr[j + 2] += aj[4] * xtemp0 + aj[5] * xtemp1;
                y_ptr[j + 3] += aj[6] * xtemp0 + aj[7] * xtemp1;
                aj += 8;

            }

            for (; j < n; j++) {
                y_ptr[j] += aj[0] * xtemp0 + aj[1] * xtemp1;
                aj += 2;
            }

        } else {
            if (inc_y == 1) {

                BLASLONG register lda2 = lda << 1;
                BLASLONG register lda4 = lda << 2;
                BLASLONG register lda3 = lda2 + lda;

                for (j = 0; j < (n & -4); j += 4) {

                    y_ptr[j] += *aj * xtemp0 + *(aj + 1) * xtemp1;
                    y_ptr[j + 1] += *(aj + lda) * xtemp0 + *(aj + lda + 1) * xtemp1;
                    y_ptr[j + 2] += *(aj + lda2) * xtemp0 + *(aj + lda2 + 1) * xtemp1;
                    y_ptr[j + 3] += *(aj + lda3) * xtemp0 + *(aj + lda3 + 1) * xtemp1;
                    aj += lda4;
                }

                for (; j < n; j++) {

                    y_ptr[j] += *aj * xtemp0 + *(aj + 1) * xtemp1;
                    aj += lda;
                }

            } else {
                for (j = 0; j < n; j++) {
                    *y_ptr += *aj * xtemp0 + *(aj + 1) * xtemp1;
                    y_ptr += inc_y;
                    aj += lda;
                }
            }

        }
        return (0);

    }

    FLOAT xtemp = *x_ptr * alpha;
    FLOAT *aj = a_ptr;
    y_ptr = y;
    if (lda == 1 && inc_y == 1) {
        for (j = 0; j < (n & -4); j += 4) {
            y_ptr[j] += aj[j] * xtemp;
            y_ptr[j + 1] += aj[j + 1] * xtemp;
            y_ptr[j + 2] += aj[j + 2] * xtemp;
            y_ptr[j + 3] += aj[j + 3] * xtemp;
        }
        for (; j < n; j++) {
            y_ptr[j] += aj[j] * xtemp;
        }


    } else {
        if (inc_y == 1) {

            BLASLONG register lda2 = lda << 1;
            BLASLONG register lda4 = lda << 2;
            BLASLONG register lda3 = lda2 + lda;
            for (j = 0; j < (n & -4); j += 4) {
                y_ptr[j] += *aj * xtemp;
                y_ptr[j + 1] += *(aj + lda) * xtemp;
                y_ptr[j + 2] += *(aj + lda2) * xtemp;
                y_ptr[j + 3] += *(aj + lda3) * xtemp;
                aj += lda4;
            }

            for (; j < n; j++) {
                y_ptr[j] += *aj * xtemp;
                aj += lda;
            }

        } else {
            for (j = 0; j < n; j++) {
                *y_ptr += *aj * xtemp;
                y_ptr += inc_y;
                aj += lda;
            }

        }
    }

    return (0);

}

#endif
