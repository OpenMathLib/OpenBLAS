/***************************************************************************
Copyright (c) 2020, The OpenBLAS Project
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

#if !defined(DOUBLE)
#define VSETVL(n)               RISCV_RVV(vsetvl_e32m8)(n)
#define VSETVLMAX               RISCV_RVV(vsetvlmax_e32m8)()
#define FLOAT_V_T               vfloat32m8_t
#define VLEV_FLOAT              RISCV_RVV(vle32_v_f32m8)
#define VLSEV_FLOAT             RISCV_RVV(vlse32_v_f32m8)
#define VSEV_FLOAT              RISCV_RVV(vse32_v_f32m8)
#define VSSEV_FLOAT             RISCV_RVV(vsse32_v_f32m8)
#define VFMACCVF_FLOAT          RISCV_RVV(vfmacc_vf_f32m8)
#else
#define VSETVL(n)               RISCV_RVV(vsetvl_e64m8)(n)
#define VSETVLMAX               RISCV_RVV(vsetvlmax_e64m8)()
#define FLOAT_V_T               vfloat64m8_t
#define VLEV_FLOAT              RISCV_RVV(vle64_v_f64m8)
#define VLSEV_FLOAT             RISCV_RVV(vlse64_v_f64m8)
#define VSEV_FLOAT              RISCV_RVV(vse64_v_f64m8)
#define VSSEV_FLOAT             RISCV_RVV(vsse64_v_f64m8)
#define VFMACCVF_FLOAT          RISCV_RVV(vfmacc_vf_f64m8)
#endif

/* Number of consecutive columns whose contribution is accumulated while a
   single y vector slice stays resident in vector registers.  Keeping the y
   slice in registers removes the per-column y load/store round-trip, while
   the fixed VL used by the main loop lets the compiler hoist the vsetvli and
   strength-reduce the pointer increments out of the hot loop. */
#define GEMV_N_BLOCK            4

int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *buffer)
{
    if (n < 0) return(0);

    FLOAT *a_ptr, *y_ptr, *a2_ptr, temp, temp2;
    BLASLONG i, j, vl;
    FLOAT_V_T va, vy, va2;

    if (inc_y == 1) {
        const BLASLONG epr = VSETVLMAX;

        /* Register-blocked main loop.  A single y vector slice is kept
           resident in registers while GEMV_N_BLOCK consecutive columns are
           accumulated into it, so the y load/store round-trip happens once
           per block instead of once per column.  VL is fixed to the hardware
           maximum for the whole main loop, which lets the compiler hoist the
           vsetvli and strength-reduce the pointer increments out of the hot
           loop.  Columns are still consumed in ascending order, therefore the
           FP64 accumulation sequence of every y[i] is bit-identical to the
           original column-at-a-time loop: the intermediate y vector is only
           held in a register rather than stored and reloaded. */
        for (j = 0; j + GEMV_N_BLOCK <= n; j += GEMV_N_BLOCK) {
            FLOAT *ap0 = a + (j + 0) * lda;
            FLOAT *ap1 = a + (j + 1) * lda;
            FLOAT *ap2 = a + (j + 2) * lda;
            FLOAT *ap3 = a + (j + 3) * lda;
            FLOAT t0 = alpha * x[(j + 0) * inc_x];
            FLOAT t1 = alpha * x[(j + 1) * inc_x];
            FLOAT t2 = alpha * x[(j + 2) * inc_x];
            FLOAT t3 = alpha * x[(j + 3) * inc_x];

            for (i = 0; i + epr <= m; i += epr) {
                y_ptr = y + i;
                vy = VLEV_FLOAT(y_ptr, epr);
                vy = VFMACCVF_FLOAT(vy, t0, VLEV_FLOAT(ap0 + i, epr), epr);
                vy = VFMACCVF_FLOAT(vy, t1, VLEV_FLOAT(ap1 + i, epr), epr);
                vy = VFMACCVF_FLOAT(vy, t2, VLEV_FLOAT(ap2 + i, epr), epr);
                vy = VFMACCVF_FLOAT(vy, t3, VLEV_FLOAT(ap3 + i, epr), epr);
                VSEV_FLOAT(y_ptr, vy, epr);
            }
            if (i < m) {
                vl = VSETVL(m - i);
                y_ptr = y + i;
                vy = VLEV_FLOAT(y_ptr, vl);
                vy = VFMACCVF_FLOAT(vy, t0, VLEV_FLOAT(ap0 + i, vl), vl);
                vy = VFMACCVF_FLOAT(vy, t1, VLEV_FLOAT(ap1 + i, vl), vl);
                vy = VFMACCVF_FLOAT(vy, t2, VLEV_FLOAT(ap2 + i, vl), vl);
                vy = VFMACCVF_FLOAT(vy, t3, VLEV_FLOAT(ap3 + i, vl), vl);
                VSEV_FLOAT(y_ptr, vy, vl);
            }
        }

        /* Remaining columns keep the original two-column (n >> 1) structure. */
        for (; j + 2 <= n; j += 2) {
            temp = alpha * x[(j + 0) * inc_x];
            temp2 = alpha * x[(j + 1) * inc_x];
            a_ptr = a + (j + 0) * lda;
            a2_ptr = a + (j + 1) * lda;
            for (i = 0; i + epr <= m; i += epr) {
                y_ptr = y + i;
                vy = VLEV_FLOAT(y_ptr, epr);
                va = VLEV_FLOAT(a_ptr + i, epr);
                va2 = VLEV_FLOAT(a2_ptr + i, epr);
                vy = VFMACCVF_FLOAT(vy, temp, va, epr);
                vy = VFMACCVF_FLOAT(vy, temp2, va2, epr);
                VSEV_FLOAT(y_ptr, vy, epr);
            }
            if (i < m) {
                vl = VSETVL(m - i);
                y_ptr = y + i;
                vy = VLEV_FLOAT(y_ptr, vl);
                va = VLEV_FLOAT(a_ptr + i, vl);
                va2 = VLEV_FLOAT(a2_ptr + i, vl);
                vy = VFMACCVF_FLOAT(vy, temp, va, vl);
                vy = VFMACCVF_FLOAT(vy, temp2, va2, vl);
                VSEV_FLOAT(y_ptr, vy, vl);
            }
        }
        /* Remaining odd column keeps the original single-column structure. */
        if (j < n) {
            temp = alpha * x[j * inc_x];
            a_ptr = a + j * lda;
            for (i = 0; i + epr <= m; i += epr) {
                y_ptr = y + i;
                vy = VLEV_FLOAT(y_ptr, epr);
                va = VLEV_FLOAT(a_ptr + i, epr);
                vy = VFMACCVF_FLOAT(vy, temp, va, epr);
                VSEV_FLOAT(y_ptr, vy, epr);
            }
            if (i < m) {
                vl = VSETVL(m - i);
                y_ptr = y + i;
                vy = VLEV_FLOAT(y_ptr, vl);
                va = VLEV_FLOAT(a_ptr + i, vl);
                vy = VFMACCVF_FLOAT(vy, temp, va, vl);
                VSEV_FLOAT(y_ptr, vy, vl);
            }
        }
    } else {
        BLASLONG stride_y = inc_y * sizeof(FLOAT);
        for (j = 0; j < (n >> 1); j++) {
            temp = alpha * x[0];
            temp2 = alpha * x[inc_x];
            y_ptr = y;
            a_ptr = a;
            a2_ptr = a + lda;
            for (i = m; i > 0; i -= vl) {
                vl = VSETVL(i);
                vy = VLSEV_FLOAT(y_ptr, stride_y, vl);
                va = VLEV_FLOAT(a_ptr, vl);
                va2 = VLEV_FLOAT(a2_ptr, vl);
                vy = VFMACCVF_FLOAT(vy, temp, va, vl);
                vy = VFMACCVF_FLOAT(vy, temp2, va2, vl);
                VSSEV_FLOAT(y_ptr, stride_y, vy, vl);
                y_ptr += vl * inc_y;
                a_ptr += vl;
                a2_ptr += vl;
            }
            x += inc_x * 2;
            a += lda * 2;
        }
        if (n & 1) {
            temp = alpha * x[0];
            y_ptr = y;
            a_ptr = a;
            for (i = m; i > 0; i -= vl) {
                vl = VSETVL(i);
                vy = VLSEV_FLOAT(y_ptr, stride_y, vl);
                va = VLEV_FLOAT(a_ptr, vl);
                vy = VFMACCVF_FLOAT(vy, temp, va, vl);
                VSSEV_FLOAT(y_ptr, stride_y, vy, vl);
                y_ptr += vl * inc_y;
                a_ptr += vl;
            }
        }
    }
    return(0);
}
