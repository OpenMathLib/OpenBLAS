/***************************************************************************
Copyright (c) 2022, The OpenBLAS Project
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
#define VSETVL(n)               __riscv_vsetvl_e32m2(n)
#define FLOAT_V_T               vfloat32m2_t
#define FLOAT_VX2_T             vfloat32m2x2_t
#define VGET_VX2                __riscv_vget_v_f32m2x2_f32m2
#define VSET_VX2                __riscv_vset_v_f32m2_f32m2x2
#define VLEV_FLOAT              __riscv_vle32_v_f32m2
#define VSEV_FLOAT              __riscv_vse32_v_f32m2
#define VSSEG2_FLOAT            __riscv_vsseg2e32_v_f32m2x2
#define VLSEG2_FLOAT            __riscv_vlseg2e32_v_f32m2x2
#define VFMACCVF_FLOAT          __riscv_vfmacc_vf_f32m2
#define VFNMSACVF_FLOAT         __riscv_vfnmsac_vf_f32m2
#define VFMULVF_FLOAT           __riscv_vfmul_vf_f32m2
#else
#define VSETVL(n)               __riscv_vsetvl_e64m2(n)
#define FLOAT_V_T               vfloat64m2_t
#define FLOAT_VX2_T             vfloat64m2x2_t
#define VGET_VX2                __riscv_vget_v_f64m2x2_f64m2
#define VSET_VX2                __riscv_vset_v_f64m2_f64m2x2
#define VLEV_FLOAT              __riscv_vle64_v_f64m2
#define VSEV_FLOAT              __riscv_vse64_v_f64m2
#define VSSEG2_FLOAT            __riscv_vsseg2e64_v_f64m2x2
#define VLSEG2_FLOAT            __riscv_vlseg2e64_v_f64m2x2
#define VFMVVF_FLOAT            __riscv_vfmv_v_f_f64m2
#define VFMACCVF_FLOAT          __riscv_vfmacc_vf_f64m2
#define VFNMSACVF_FLOAT         __riscv_vfnmsac_vf_f64m2
#define VFMULVF_FLOAT           __riscv_vfmul_vf_f64m2
#endif


static FLOAT dm1 = -1.;

#ifdef CONJ
#define GEMM_KERNEL   GEMM_KERNEL_L
#else
#define GEMM_KERNEL   GEMM_KERNEL_N
#endif

#if GEMM_DEFAULT_UNROLL_N == 1
#define GEMM_UNROLL_N_SHIFT 0
#endif

#if GEMM_DEFAULT_UNROLL_N == 2
#define GEMM_UNROLL_N_SHIFT 1
#endif

#if GEMM_DEFAULT_UNROLL_N == 4
#define GEMM_UNROLL_N_SHIFT 2
#endif

#if GEMM_DEFAULT_UNROLL_N == 8
#define GEMM_UNROLL_N_SHIFT 3
#endif

#if GEMM_DEFAULT_UNROLL_N == 16
#define GEMM_UNROLL_N_SHIFT 4
#endif

// Optimizes the implementation in ../arm64/trsm_kernel_LT_sve.c

#ifndef COMPLEX

static inline void solve(BLASLONG m, BLASLONG n, FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc) {

    FLOAT aa;
    FLOAT* pc;

    int i, j, k;

    FLOAT_V_T va, vc;

    size_t vl;

    for (i = 0; i < m; i++) {

        aa = *(a + i);
        for (j = 0; j < n; j++) {
            FLOAT bb;

            pc = c + j * ldc;
            bb = *(pc + i) * aa;
            *(b + j) = bb;
            *(pc + i) = bb;
        }

        for (k = i + 1; k < m; k += vl) {
            vl = VSETVL(m - k);
            va = VLEV_FLOAT(a + k, vl);
            for (j = 0; j < n; j++) {
                pc = c + j * ldc;
                vc = VLEV_FLOAT(pc + k, vl);
                vc = VFNMSACVF_FLOAT(vc, *(b + j), va, vl);
                VSEV_FLOAT(pc + k, vc, vl);
            }
        }

        b += n;
        a += m;
    }
}

#else

static inline void solve(BLASLONG m, BLASLONG n, FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc) {

    FLOAT aa1, aa2;
    FLOAT *pc;
    int i, j, k;

    FLOAT_VX2_T vax2, vcx2;
    FLOAT_V_T va1, va2, vc1, vc2;
    size_t vl;

    ldc *= 2;

    for (i = 0; i < m; i++) {
        aa1 = *(a + i * 2 + 0);
        aa2 = *(a + i * 2 + 1);
        for (j = 0; j < n; j++) {
            FLOAT bb1, bb2, ss1, ss2;

            pc = c + j * ldc;
            bb1 = *(pc + i * 2 + 0);
            bb2 = *(pc + i * 2 + 1);
#ifndef CONJ
            ss1 = aa1 * bb1 - aa2 * bb2;
            ss2 = aa1 * bb2 + aa2 * bb1;
#else
            ss1 = aa1 * bb1 + aa2 * bb2;
            ss2 = aa1 * bb2 - aa2 * bb1;
#endif
            *(b + j * 2 + 0) = ss1;
            *(b + j * 2 + 1) = ss2;
            *(pc + i * 2 + 0) = ss1;
            *(pc + i * 2 + 1) = ss2;
        }

        for (k = i + 1; k < m; k += vl) {
            vl = VSETVL(m - k);
            vax2 = VLSEG2_FLOAT(a + k * 2, vl);
            va1 = VGET_VX2(vax2, 0);
            va2 = VGET_VX2(vax2, 1);
            for (j = 0; j < n; j++) {
                FLOAT ss1 = *(b + j * 2 + 0);
                FLOAT ss2 = *(b + j * 2 + 1);

                pc = c + j * ldc;
                vcx2 = VLSEG2_FLOAT(pc + k * 2, vl);
                vc1 = VGET_VX2(vcx2, 0);
                vc2 = VGET_VX2(vcx2, 1);
#ifndef CONJ
                vc1 =  VFMACCVF_FLOAT(vc1, ss2, va2, vl);
                vc1 = VFNMSACVF_FLOAT(vc1, ss1, va1, vl);
                vc2 = VFNMSACVF_FLOAT(vc2, ss1, va2, vl);
                vc2 = VFNMSACVF_FLOAT(vc2, ss2, va1, vl);
#else
                vc1 = VFNMSACVF_FLOAT(vc1, ss2, va2, vl);
                vc1 = VFNMSACVF_FLOAT(vc1, ss1, va1, vl);
                vc2 =  VFMACCVF_FLOAT(vc2, ss1, va2, vl);
                vc2 = VFNMSACVF_FLOAT(vc2, ss2, va1, vl);
#endif
                vcx2 = VSET_VX2(vcx2, 0, vc1);
                vcx2 = VSET_VX2(vcx2, 1, vc2);
                VSSEG2_FLOAT(pc + k * 2, vcx2, vl);
            }
        }

        b += n * 2;
        a += m * 2;
    }
}

#endif

int CNAME(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT dummy1,
#ifdef COMPLEX
	   FLOAT dummy2,
#endif
	   FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG offset){

  FLOAT *aa, *cc;
  BLASLONG  kk;
  BLASLONG i, j;

#ifndef COMPLEX
#define PROCESS_LT_M_BLOCK(MB, NB) do { \
    if (kk > 0) { \
      GEMM_KERNEL((MB), (NB), kk, dm1, aa, b, cc, ldc); \
    } \
    solve((MB), (NB), \
          aa + kk * (MB) * COMPSIZE, \
          b  + kk * (NB) * COMPSIZE, \
          cc, ldc); \
    aa += (MB) * k * COMPSIZE; \
    cc += (MB)     * COMPSIZE; \
    kk += (MB); \
  } while (0)
#else
#define PROCESS_LT_M_BLOCK(MB, NB) do { \
    if (kk > 0) { \
      GEMM_KERNEL((MB), (NB), kk, dm1, ZERO, aa, b, cc, ldc); \
    } \
    solve((MB), (NB), \
          aa + kk * (MB) * COMPSIZE, \
          b  + kk * (NB) * COMPSIZE, \
          cc, ldc); \
    aa += (MB) * k * COMPSIZE; \
    cc += (MB)     * COMPSIZE; \
    kk += (MB); \
  } while (0)
#endif

  j = (n >> GEMM_UNROLL_N_SHIFT);

  while (j > 0) {

    kk = offset;
    aa = a;
    cc = c;

    i = 0;
    while (i + GEMM_DEFAULT_UNROLL_M <= m) {
      PROCESS_LT_M_BLOCK(GEMM_DEFAULT_UNROLL_M, GEMM_UNROLL_N);
      i += GEMM_DEFAULT_UNROLL_M;
    }

    if (m & (GEMM_DEFAULT_UNROLL_M - 1)) {
      BLASLONG mm = (GEMM_DEFAULT_UNROLL_M >> 1);
      while (mm > 0) {
        if ((m - i) & mm) {
          PROCESS_LT_M_BLOCK(mm, GEMM_UNROLL_N);
          i += mm;
        }
        mm >>= 1;
      }
    }

    b += GEMM_UNROLL_N * k   * COMPSIZE;
    c += GEMM_UNROLL_N * ldc * COMPSIZE;
    j --;
  }

  if (n & (GEMM_UNROLL_N - 1)) {

    j = (GEMM_UNROLL_N >> 1);
    while (j > 0) {
      if (n & j) {

        kk = offset;
        aa = a;
        cc = c;

        i = 0;
        while (i + GEMM_DEFAULT_UNROLL_M <= m) {
          PROCESS_LT_M_BLOCK(GEMM_DEFAULT_UNROLL_M, j);
          i += GEMM_DEFAULT_UNROLL_M;
        }

        if (m & (GEMM_DEFAULT_UNROLL_M - 1)) {
          BLASLONG mm = (GEMM_DEFAULT_UNROLL_M >> 1);
          while (mm > 0) {
            if ((m - i) & mm) {
              PROCESS_LT_M_BLOCK(mm, j);
              i += mm;
            }
            mm >>= 1;
          }
        }

        b += j * k   * COMPSIZE;
        c += j * ldc * COMPSIZE;
      }
      j >>= 1;
    }
  }

  return 0;

#undef PROCESS_LT_M_BLOCK
}
