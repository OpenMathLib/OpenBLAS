/***************************************************************************
Copyright (c) 2026, The OpenBLAS Project
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
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#include "common.h"
#include <math.h>

/*
 * Two paths, chosen from measurements on the SpacemiT X100 (VLEN 256):
 *  - data up to STREAM_BYTES: one LMUL-8 accumulator and a single final
 *    reduction (several times faster in cache than reducing every strip);
 *  - larger unit-stride data, which streams from DRAM: four LMUL-1
 *    accumulators. On the X100 small vector loads sustain a markedly higher
 *    single-core DRAM read bandwidth than LMUL-8 loads.
 */
#define STREAM_BYTES (2 * 1024 * 1024)

#if !defined(DOUBLE)
#define VSETVL_M8(n)            __riscv_vsetvl_e32m8(n)
#define VSETVL_MAX_M8           __riscv_vsetvlmax_e32m8()
#define FLOAT_V_T_M8            vfloat32m8_t
#define VLEV_M8                 __riscv_vle32_v_f32m8
#define VLSEV_M8                __riscv_vlse32_v_f32m8
#define VFABS_M8                __riscv_vfabs_v_f32m8
#define VFMVVF_M8               __riscv_vfmv_v_f_f32m8
#define VSETVL_M1(n)            __riscv_vsetvl_e32m1(n)
#define VSETVL_MAX_M1           __riscv_vsetvlmax_e32m1()
#define FLOAT_V_T_M1            vfloat32m1_t
#define VLEV_M1                 __riscv_vle32_v_f32m1
#define VFABS_M1                __riscv_vfabs_v_f32m1
#define VFMVVF_M1               __riscv_vfmv_v_f_f32m1
#define VFMVFS_M1               __riscv_vfmv_f_s_f32m1_f32
#define ABS_FLOAT               fabsf
#else
#define VSETVL_M8(n)            __riscv_vsetvl_e64m8(n)
#define VSETVL_MAX_M8           __riscv_vsetvlmax_e64m8()
#define FLOAT_V_T_M8            vfloat64m8_t
#define VLEV_M8                 __riscv_vle64_v_f64m8
#define VLSEV_M8                __riscv_vlse64_v_f64m8
#define VFABS_M8                __riscv_vfabs_v_f64m8
#define VFMVVF_M8               __riscv_vfmv_v_f_f64m8
#define VSETVL_M1(n)            __riscv_vsetvl_e64m1(n)
#define VSETVL_MAX_M1           __riscv_vsetvlmax_e64m1()
#define FLOAT_V_T_M1            vfloat64m1_t
#define VLEV_M1                 __riscv_vle64_v_f64m1
#define VFABS_M1                __riscv_vfabs_v_f64m1
#define VFMVVF_M1               __riscv_vfmv_v_f_f64m1
#define VFMVFS_M1               __riscv_vfmv_f_s_f64m1_f64
#define ABS_FLOAT               fabs
#endif
#if !defined(DOUBLE)
#define VFMAX_M8_TU              __riscv_vfmax_vv_f32m8_tu
#define VFREDMAX_M8              __riscv_vfredmax_vs_f32m8_f32m1
#define VFMAX_M1                 __riscv_vfmax_vv_f32m1
#define VFMAX_M1_TU              __riscv_vfmax_vv_f32m1_tu
#define VFREDMAX_M1              __riscv_vfredmax_vs_f32m1_f32m1
#endif
#if defined(DOUBLE)
#define VFMAX_M8_TU              __riscv_vfmax_vv_f64m8_tu
#define VFREDMAX_M8              __riscv_vfredmax_vs_f64m8_f64m1
#define VFMAX_M1                 __riscv_vfmax_vv_f64m1
#define VFMAX_M1_TU              __riscv_vfmax_vv_f64m1_tu
#define VFREDMAX_M1              __riscv_vfredmax_vs_f64m1_f64m1
#endif

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
    if (n <= 0 || inc_x <= 0) return(0.0);

    // The result starts from 0, as in amax_vector.c.
    FLOAT init = 0.0;
    FLOAT_V_T_M1 v_res = VFMVVF_M1(init, VSETVL_MAX_M1);

    if (inc_x == 1 && (size_t)n * sizeof(FLOAT) >= STREAM_BYTES) {
        size_t vl1 = VSETVL_MAX_M1;
        FLOAT_V_T_M1 a0 = VFMVVF_M1(init, vl1), a1 = a0, a2 = a0, a3 = a0;
        BLASLONG j = 0, step = 4 * (BLASLONG)vl1;
        for (; j + step <= n; j += step) {
            a0 = VFMAX_M1(a0, VFABS_M1(VLEV_M1(&x[j], vl1), vl1), vl1);
            a1 = VFMAX_M1(a1, VFABS_M1(VLEV_M1(&x[j + vl1], vl1), vl1), vl1);
            a2 = VFMAX_M1(a2, VFABS_M1(VLEV_M1(&x[j + 2 * vl1], vl1), vl1), vl1);
            a3 = VFMAX_M1(a3, VFABS_M1(VLEV_M1(&x[j + 3 * vl1], vl1), vl1), vl1);
        }
        for (size_t vl; j < n; j += vl) {
            vl = VSETVL_M1(n - j);
            a0 = VFMAX_M1_TU(a0, a0, VFABS_M1(VLEV_M1(&x[j], vl), vl), vl);
        }
        a0 = VFMAX_M1(VFMAX_M1(a0, a1, vl1), VFMAX_M1(a2, a3, vl1), vl1);
        v_res = VFREDMAX_M1(a0, v_res, vl1);
        return(VFMVFS_M1(v_res));
    }

    size_t vlmax = VSETVL_MAX_M8;
    FLOAT_V_T_M8 acc = VFMVVF_M8(init, vlmax);
    if (inc_x == 1) {
        for (size_t vl; n > 0; n -= vl, x += vl) {
            vl = VSETVL_M8(n);
            acc = VFMAX_M8_TU(acc, acc, VFABS_M8(VLEV_M8(x, vl), vl), vl);
        }
    } else {
        BLASLONG stride_x = inc_x * sizeof(FLOAT);
        for (size_t vl; n > 0; n -= vl, x += vl * inc_x) {
            vl = VSETVL_M8(n);
            acc = VFMAX_M8_TU(acc, acc, VFABS_M8(VLSEV_M8(x, stride_x, vl), vl), vl);
        }
    }
    v_res = VFREDMAX_M8(acc, v_res, vlmax);
    return(VFMVFS_M1(v_res));
}
