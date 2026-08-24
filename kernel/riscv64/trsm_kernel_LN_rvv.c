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

// Optimizes the implementation in ../arm64/trsm_kernel_LN_sve.c

#ifndef COMPLEX

static inline void solve(BLASLONG m, BLASLONG n, FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc) {
    FLOAT aa;
    FLOAT* pc;

    int i, j, k;

    FLOAT_V_T va, vc;

    size_t vl;

    a += (m - 1) * m;
    b += (m - 1) * n;

    for (i = m - 1; i >= 0; i--) {

        aa = *(a + i);
        for (j = 0; j < n; j++) {
            FLOAT bb;

            pc = c + j * ldc;
            bb = *(pc + i) * aa;
            *(b + j) = bb;
            *(pc + i) = bb;
        }

        for (k = 0; k < i; k += vl) {
            vl = VSETVL(i - k);
            va = VLEV_FLOAT(a + k, vl);
            for (j = 0; j < n; j++) {
                pc = c + j * ldc;
                vc = VLEV_FLOAT(pc + k, vl);
                vc = VFNMSACVF_FLOAT(vc, *(b + j), va, vl);
                VSEV_FLOAT(pc + k, vc, vl);
            }
        }

        a -= m;
        b += n;
        b -= 2 * n;
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
    BLASLONG ldc2 = ldc * 2;

    a += (m - 1) * m * 2;
    b += (m - 1) * n * 2;

    for (i = m - 1; i >= 0; i--) {

        aa1 = *(a + i * 2 + 0);
        aa2 = *(a + i * 2 + 1);
        for (j = 0; j < n; j++) {
            FLOAT bb1, bb2, ss1, ss2;

            pc = c + j * ldc2;
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

        for (k = 0; k < i; k += vl) {
            vl = VSETVL(i - k);
            vax2 = VLSEG2_FLOAT(a + k * 2, vl);
            va1 = VGET_VX2(vax2, 0);
            va2 = VGET_VX2(vax2, 1);
            for (j = 0; j < n; j++) {
                FLOAT ss1 = *(b + j * 2 + 0);
                FLOAT ss2 = *(b + j * 2 + 1);

                pc = c + j * ldc2;
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

        a -= m * 2;
        b += n * 2;
        b -= 4 * n;
    }
}


#endif

int CNAME(BLASLONG m, BLASLONG n, BLASLONG k,  FLOAT dummy1,
#ifdef COMPLEX
    FLOAT dummy2,
#endif
    FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG offset){

  BLASLONG i, j;
  FLOAT *aa, *cc;
  BLASLONG  kk;

#ifndef COMPLEX
#define PROCESS_LN_M_BLOCK(MB, NB) do { \
    if (k - kk > 0) { \
      GEMM_KERNEL((MB), (NB), k - kk, dm1, \
                  aa + (MB) * kk * COMPSIZE, \
                  b  + (NB) * kk * COMPSIZE, \
                  cc, ldc); \
    } \
    solve((MB), (NB), \
          aa + (kk - (MB)) * (MB) * COMPSIZE, \
          b  + (kk - (MB)) * (NB) * COMPSIZE, \
          cc, ldc); \
  } while (0)
#else
#define PROCESS_LN_M_BLOCK(MB, NB) do { \
    if (k - kk > 0) { \
      GEMM_KERNEL((MB), (NB), k - kk, dm1, ZERO, \
                  aa + (MB) * kk * COMPSIZE, \
                  b  + (NB) * kk * COMPSIZE, \
                  cc, ldc); \
    } \
    solve((MB), (NB), \
          aa + (kk - (MB)) * (MB) * COMPSIZE, \
          b  + (kk - (MB)) * (NB) * COMPSIZE, \
          cc, ldc); \
  } while (0)
#endif

  j = (n >> GEMM_UNROLL_N_SHIFT);

  while (j > 0) {

    kk = m + offset;

    if (m & (GEMM_DEFAULT_UNROLL_M - 1)) {
      for (i = 1; i < GEMM_DEFAULT_UNROLL_M; i <<= 1) {
        if (m & i) {
          aa = a + ((m & ~(i - 1)) - i) * k * COMPSIZE;
          cc = c + ((m & ~(i - 1)) - i)     * COMPSIZE;

          PROCESS_LN_M_BLOCK(i, GEMM_UNROLL_N);
          kk -= i;
        }
      }
    }

    i = (m & ~(GEMM_DEFAULT_UNROLL_M - 1));
    if (i > 0) {
      aa = a + (i - GEMM_DEFAULT_UNROLL_M) * k * COMPSIZE;
      cc = c + (i - GEMM_DEFAULT_UNROLL_M)     * COMPSIZE;

      do {
        PROCESS_LN_M_BLOCK(GEMM_DEFAULT_UNROLL_M, GEMM_UNROLL_N);

        aa -= GEMM_DEFAULT_UNROLL_M * k * COMPSIZE;
        cc -= GEMM_DEFAULT_UNROLL_M     * COMPSIZE;
        kk -= GEMM_DEFAULT_UNROLL_M;
        i -= GEMM_DEFAULT_UNROLL_M;
      } while (i > 0);
    }

    b += GEMM_UNROLL_N * k * COMPSIZE;
    c += GEMM_UNROLL_N * ldc * COMPSIZE;
    j --;
  }

  if (n & (GEMM_UNROLL_N - 1)) {

    j = (GEMM_UNROLL_N >> 1);
    while (j > 0) {
      if (n & j) {

        kk = m + offset;

        if (m & (GEMM_DEFAULT_UNROLL_M - 1)) {
          for (i = 1; i < GEMM_DEFAULT_UNROLL_M; i <<= 1) {
            if (m & i) {
              aa = a + ((m & ~(i - 1)) - i) * k * COMPSIZE;
              cc = c + ((m & ~(i - 1)) - i)     * COMPSIZE;

              PROCESS_LN_M_BLOCK(i, j);
              kk -= i;
            }
          }
        }

        i = (m & ~(GEMM_DEFAULT_UNROLL_M - 1));
        if (i > 0) {
          aa = a + (i - GEMM_DEFAULT_UNROLL_M) * k * COMPSIZE;
          cc = c + (i - GEMM_DEFAULT_UNROLL_M)     * COMPSIZE;

          do {
            PROCESS_LN_M_BLOCK(GEMM_DEFAULT_UNROLL_M, j);

            aa -= GEMM_DEFAULT_UNROLL_M * k * COMPSIZE;
            cc -= GEMM_DEFAULT_UNROLL_M     * COMPSIZE;
            kk -= GEMM_DEFAULT_UNROLL_M;
            i -= GEMM_DEFAULT_UNROLL_M;
          } while (i > 0);
        }

        b += j * k   * COMPSIZE;
        c += j * ldc * COMPSIZE;
      }
      j >>= 1;
    }
  }

  return 0;

#undef PROCESS_LN_M_BLOCK
}
