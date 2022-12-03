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
#define VSETVL(n) vsetvl_e32m2(n)
#define VSETVL_MAX vsetvlmax_e32m2()
#define FLOAT_V_T vfloat32m2_t
#define VLEV_FLOAT vle32_v_f32m2
#define VLSEV_FLOAT vlse32_v_f32m2
#define VLSEG2_FLOAT vlseg2e32_v_f32m2
#define VSEV_FLOAT vse32_v_f32m2
#define VSSEV_FLOAT vsse32_v_f32m2
#define VSSEG2_FLOAT vsseg2e32_v_f32m2
#define VFMACCVF_FLOAT vfmacc_vf_f32m2
#define VFMULVF_FLOAT vfmul_vf_f32m2
#define VFNMSACVF_FLOAT vfnmsac_vf_f32m2
#else
#define VSETVL(n) vsetvl_e64m2(n)
#define VSETVL_MAX vsetvlmax_e64m2()
#define FLOAT_V_T vfloat64m2_t
#define VLEV_FLOAT vle64_v_f64m2
#define VLSEV_FLOAT vlse64_v_f64m2
#define VLSEG2_FLOAT vlseg2e64_v_f64m2
#define VSEV_FLOAT vse64_v_f64m2
#define VSSEV_FLOAT vsse64_v_f64m2
#define VSSEG2_FLOAT vsseg2e64_v_f64m2
#define VFMACCVF_FLOAT vfmacc_vf_f64m2
#define VFMULVF_FLOAT vfmul_vf_f64m2
#define VFNMSACVF_FLOAT vfnmsac_vf_f64m2
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

#if GEMM_DEFAULT_UNROLL_N == 1

static inline void solve(BLASLONG m, BLASLONG n, FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc) {

    FLOAT aa,  bb;
    FLOAT *pa, *pc;

    int i, j, k;
    //fprintf(stderr, "%s , %s, m = %4ld  n = %4ld  offset = %4ld\n", __FILE__, __FUNCTION__, m, n, ldc); // Debug

    size_t vl;
    FLOAT_V_T va, vc;

    a += (m - 1) * m;
    b += (m - 1) * n;

    for (i = m - 1; i >= 0; i--) 
    {
        aa = *(a + i);
        for (j = 0; j < n; j ++) 
        {
            bb = *(c + i + j * ldc);
            bb *= aa;
            *b             = bb;
            *(c + i + j * ldc) = bb;
            b ++;

            pa = a;
            pc = c + j * ldc;
            for (k = i; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                vc = VLEV_FLOAT(pc, vl);
                va = VLEV_FLOAT(pa, vl);
                vc = VFNMSACVF_FLOAT(vc, bb, va, vl);
                VSEV_FLOAT(pc, vc, vl);
                pa += vl;
                pc += vl;
            }
        }
        a -= m;
        b -= 2 * n;
    }

}
#elif GEMM_DEFAULT_UNROLL_N == 2

static inline void solve(BLASLONG m, BLASLONG n, FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc) {

    FLOAT aa,  bb0, bb1;
    FLOAT *pa, *pc, *pc0, *pc1;
    FLOAT *pb0, *pb1;

    int i, j, k;
    fprintf(stderr, "%s , %s, m = %4ld  n = %4ld  offset = %4ld\n", __FILE__, __FUNCTION__, m, n, ldc); // Debug

    size_t vl;
    FLOAT_V_T va, vc0, vc1;

    a += (m - 1) * m;
    b += (m - 1) * n;

    for (i = m - 1; i >= 0; i--) 
    {
        aa = *(a + i);
        pc = c + i;
        for (j = 0; j < n/2; j ++) 
        {
            //bb = *(c + i + j * ldc);
            pb0 = pc + j * ldc * 2;
            pb1 = pb0 + ldc;
            //bb *= aa;
            bb0 = (*pb0) * aa;
            bb1 = (*pb1) * aa;
            //*b             = bb;
            *b      = bb0;
            *(b+1)  = bb1;
            *pb0    = bb0;
            *pb1    = bb1;

            //*(c + i + j * ldc) = bb;
            //b ++;

            b += 2;
            //pa = a + i + 1;
            pc0 = c + j * ldc * 2;
            pc1 = pc0 + ldc;
            pa = a;
            //pc = c + j * ldc;
            for (k = i; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                vc0 = VLEV_FLOAT(pc0, vl);
                vc1 = VLEV_FLOAT(pc1, vl);
                va = VLEV_FLOAT(pa, vl);
                vc0 = VFNMSACVF_FLOAT(vc0, bb0, va, vl);
                vc1 = VFNMSACVF_FLOAT(vc1, bb1, va, vl);
                VSEV_FLOAT(pc0, vc0, vl);
                VSEV_FLOAT(pc1, vc1, vl);

                pa += vl;
                pc0 += vl;
                pc1 += vl;
            }
        }
        pc += ldc * (n/2) * 2;
        if (n & 1)
        {
            pb0 = pc;
            bb0 = (*pb0) * aa;
            *b      = bb0;
            *pb0    = bb0;
            b += 1;

            pc0 = pc - i;
            pa = a;
            for (k = i; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                vc0 = VLEV_FLOAT(pc0, vl);
                va = VLEV_FLOAT(pa, vl);
                vc0 = VFNMSACVF_FLOAT(vc0, bb0, va, vl);
                VSEV_FLOAT(pc0, vc0, vl);

                pa += vl;
                pc0 += vl;
            }
        }

        a -= m;
        b -= 2 * n;
    }

}

#elif GEMM_DEFAULT_UNROLL_N == 4

static inline void solve(BLASLONG m, BLASLONG n, FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc) {

    FLOAT aa,  bb0, bb1, bb2, bb3;
    FLOAT *pa, *pc, *pc0, *pc1, *pc2, *pc3;
    FLOAT *pb0, *pb1, *pb2, *pb3;

    int i, j, k;

    size_t vl;
    FLOAT_V_T va, vc0, vc1, vc2, vc3;

    a += (m - 1) * m;
    b += (m - 1) * n;

    for (i = m - 1; i >= 0; i--) 
    {
        aa = *(a + i);
        pc = c + i;
        for (j = 0; j < n/4; j ++) 
        {
            pb0 = pc + j * ldc * 4;
            pb1 = pb0 + ldc;
            pb2 = pb1 + ldc;
            pb3 = pb2 + ldc;
            
            bb0 = (*pb0) * aa;
            bb1 = (*pb1) * aa;
            bb2 = (*pb2) * aa;
            bb3 = (*pb3) * aa;

            *b      = bb0;
            *(b+1)  = bb1;
            *(b+2)  = bb2;
            *(b+3)  = bb3;

            *pb0    = bb0;
            *pb1    = bb1;
            *pb2    = bb2;
            *pb3    = bb3;

            b += 4;

            pc0 = c + j * ldc * 4;
            pc1 = pc0 + ldc;
            pc2 = pc1 + ldc;
            pc3 = pc2 + ldc;

            pa = a;
            for (k = i; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                vc0 = VLEV_FLOAT(pc0, vl);
                vc1 = VLEV_FLOAT(pc1, vl);
                vc2 = VLEV_FLOAT(pc2, vl);
                vc3 = VLEV_FLOAT(pc3, vl);
                va = VLEV_FLOAT(pa, vl);
                vc0 = VFNMSACVF_FLOAT(vc0, bb0, va, vl);
                vc1 = VFNMSACVF_FLOAT(vc1, bb1, va, vl);
                vc2 = VFNMSACVF_FLOAT(vc2, bb2, va, vl);
                vc3 = VFNMSACVF_FLOAT(vc3, bb3, va, vl);
                VSEV_FLOAT(pc0, vc0, vl);
                VSEV_FLOAT(pc1, vc1, vl);
                VSEV_FLOAT(pc2, vc2, vl);
                VSEV_FLOAT(pc3, vc3, vl);

                pa += vl;
                pc0 += vl;
                pc1 += vl;
                pc2 += vl;
                pc3 += vl;
            }
        }
        pc += ldc * (n/4) * 4;

        if (n & 2)
        {
            pb0 = pc + j * ldc * 2;
            pb1 = pb0 + ldc;
            
            bb0 = (*pb0) * aa;
            bb1 = (*pb1) * aa;

            *b      = bb0;
            *(b+1)  = bb1;

            *pb0    = bb0;
            *pb1    = bb1;

            b += 2;

            pc0 = c + j * ldc * 2;
            pc1 = pc0 + ldc;

            pa = a;
            for (k = i; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                vc0 = VLEV_FLOAT(pc0, vl);
                vc1 = VLEV_FLOAT(pc1, vl);
                va = VLEV_FLOAT(pa, vl);
                vc0 = VFNMSACVF_FLOAT(vc0, bb0, va, vl);
                vc1 = VFNMSACVF_FLOAT(vc1, bb1, va, vl);
                VSEV_FLOAT(pc0, vc0, vl);
                VSEV_FLOAT(pc1, vc1, vl);

                pa += vl;
                pc0 += vl;
                pc1 += vl;
            }
            pc += ldc * 2;
        }

        if (n & 1)
        {
            pb0 = pc;
            bb0 = (*pb0) * aa;
            *b      = bb0;
            *pb0    = bb0;
            b += 1;

            pc0 = pc - i;
            pa = a;
            for (k = i; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                vc0 = VLEV_FLOAT(pc0, vl);
                va = VLEV_FLOAT(pa, vl);
                vc0 = VFNMSACVF_FLOAT(vc0, bb0, va, vl);
                VSEV_FLOAT(pc0, vc0, vl);

                pa += vl;
                pc0 += vl;
            }
        }

        a -= m;
        b -= 2 * n;
    }

}
#elif GEMM_DEFAULT_UNROLL_N == 8

static inline void solve(BLASLONG m, BLASLONG n, FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc) {

    FLOAT aa,  bb0, bb1, bb2, bb3, bb4, bb5, bb6, bb7;
    FLOAT *pa, *pc, *pc0, *pc1, *pc2, *pc3, *pc4, *pc5, *pc6, *pc7;
    FLOAT *pb0, *pb1, *pb2, *pb3, *pb4, *pb5, *pb6, *pb7;

    int i, j, k;

    size_t vl;
    FLOAT_V_T va, vc0, vc1, vc2, vc3, vc4, vc5, vc6, vc7;

    a += (m - 1) * m;
    b += (m - 1) * n;

    for (i = m - 1; i >= 0; i--) 
    {
        aa = *(a + i);
        pc = c + i;
        for (j = 0; j < n/8; j ++) 
        {
            pb0 = pc + j * ldc * 8;
            pb1 = pb0 + ldc;
            pb2 = pb1 + ldc;
            pb3 = pb2 + ldc;
            pb4 = pb3 + ldc;
            pb5 = pb4 + ldc;
            pb6 = pb5 + ldc;
            pb7 = pb6 + ldc;
            
            bb0 = (*pb0) * aa;
            bb1 = (*pb1) * aa;
            bb2 = (*pb2) * aa;
            bb3 = (*pb3) * aa;
            bb4 = (*pb4) * aa;
            bb5 = (*pb5) * aa;
            bb6 = (*pb6) * aa;
            bb7 = (*pb7) * aa;

            *b      = bb0;
            *(b+1)  = bb1;
            *(b+2)  = bb2;
            *(b+3)  = bb3;
            *(b+4)  = bb4;
            *(b+5)  = bb5;
            *(b+6)  = bb6;
            *(b+7)  = bb7;

            *pb0    = bb0;
            *pb1    = bb1;
            *pb2    = bb2;
            *pb3    = bb3;
            *pb4    = bb4;
            *pb5    = bb5;
            *pb6    = bb6;
            *pb7    = bb7;

            b += 8;

            pc0 = c + j * ldc * 8;
            pc1 = pc0 + ldc;
            pc2 = pc1 + ldc;
            pc3 = pc2 + ldc;
            pc4 = pc3 + ldc;
            pc5 = pc4 + ldc;
            pc6 = pc5 + ldc;
            pc7 = pc6 + ldc;

            pa = a;
            for (k = i; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                vc0 = VLEV_FLOAT(pc0, vl);
                vc1 = VLEV_FLOAT(pc1, vl);
                vc2 = VLEV_FLOAT(pc2, vl);
                vc3 = VLEV_FLOAT(pc3, vl);
                vc4 = VLEV_FLOAT(pc4, vl);
                vc5 = VLEV_FLOAT(pc5, vl);
                vc6 = VLEV_FLOAT(pc6, vl);
                vc7 = VLEV_FLOAT(pc7, vl);
                va = VLEV_FLOAT(pa, vl);
                vc0 = VFNMSACVF_FLOAT(vc0, bb0, va, vl);
                vc1 = VFNMSACVF_FLOAT(vc1, bb1, va, vl);
                vc2 = VFNMSACVF_FLOAT(vc2, bb2, va, vl);
                vc3 = VFNMSACVF_FLOAT(vc3, bb3, va, vl);
                vc4 = VFNMSACVF_FLOAT(vc4, bb4, va, vl);
                vc5 = VFNMSACVF_FLOAT(vc5, bb5, va, vl);
                vc6 = VFNMSACVF_FLOAT(vc6, bb6, va, vl);
                vc7 = VFNMSACVF_FLOAT(vc7, bb7, va, vl);
                VSEV_FLOAT(pc0, vc0, vl);
                VSEV_FLOAT(pc1, vc1, vl);
                VSEV_FLOAT(pc2, vc2, vl);
                VSEV_FLOAT(pc3, vc3, vl);
                VSEV_FLOAT(pc4, vc4, vl);
                VSEV_FLOAT(pc5, vc5, vl);
                VSEV_FLOAT(pc6, vc6, vl);
                VSEV_FLOAT(pc7, vc7, vl);

                pa += vl;
                pc0 += vl;
                pc1 += vl;
                pc2 += vl;
                pc3 += vl;
                pc4 += vl;
                pc5 += vl;
                pc6 += vl;
                pc7 += vl;
            }
        }
        pc += ldc * (n/8) * 8;

        if (n & 4)
        {
            pb0 = pc + j * ldc * 4;
            pb1 = pb0 + ldc;
            pb2 = pb1 + ldc;
            pb3 = pb2 + ldc;
            
            bb0 = (*pb0) * aa;
            bb1 = (*pb1) * aa;
            bb2 = (*pb2) * aa;
            bb3 = (*pb3) * aa;

            *b      = bb0;
            *(b+1)  = bb1;
            *(b+2)  = bb2;
            *(b+3)  = bb3;

            *pb0    = bb0;
            *pb1    = bb1;
            *pb2    = bb2;
            *pb3    = bb3;

            b += 4;

            pc0 = c + j * ldc * 4;
            pc1 = pc0 + ldc;
            pc2 = pc1 + ldc;
            pc3 = pc2 + ldc;

            pa = a;
            for (k = i; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                vc0 = VLEV_FLOAT(pc0, vl);
                vc1 = VLEV_FLOAT(pc1, vl);
                vc2 = VLEV_FLOAT(pc2, vl);
                vc3 = VLEV_FLOAT(pc3, vl);
                va = VLEV_FLOAT(pa, vl);
                vc0 = VFNMSACVF_FLOAT(vc0, bb0, va, vl);
                vc1 = VFNMSACVF_FLOAT(vc1, bb1, va, vl);
                vc2 = VFNMSACVF_FLOAT(vc2, bb2, va, vl);
                vc3 = VFNMSACVF_FLOAT(vc3, bb3, va, vl);
                VSEV_FLOAT(pc0, vc0, vl);
                VSEV_FLOAT(pc1, vc1, vl);
                VSEV_FLOAT(pc2, vc2, vl);
                VSEV_FLOAT(pc3, vc3, vl);

                pa += vl;
                pc0 += vl;
                pc1 += vl;
                pc2 += vl;
                pc3 += vl;
            }
            pc += ldc * 4;
        }

        if (n & 2)
        {
            pb0 = pc + j * ldc * 2;
            pb1 = pb0 + ldc;
            
            bb0 = (*pb0) * aa;
            bb1 = (*pb1) * aa;

            *b      = bb0;
            *(b+1)  = bb1;

            *pb0    = bb0;
            *pb1    = bb1;

            b += 2;

            pc0 = c + j * ldc * 2;
            pc1 = pc0 + ldc;

            pa = a;
            for (k = i; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                vc0 = VLEV_FLOAT(pc0, vl);
                vc1 = VLEV_FLOAT(pc1, vl);
                va = VLEV_FLOAT(pa, vl);
                vc0 = VFNMSACVF_FLOAT(vc0, bb0, va, vl);
                vc1 = VFNMSACVF_FLOAT(vc1, bb1, va, vl);
                VSEV_FLOAT(pc0, vc0, vl);
                VSEV_FLOAT(pc1, vc1, vl);

                pa += vl;
                pc0 += vl;
                pc1 += vl;
            }
            pc += ldc * 2;
        }

        if (n & 1)
        {
            pb0 = pc;
            bb0 = (*pb0) * aa;
            *b      = bb0;
            *pb0    = bb0;
            b += 1;

            pc0 = pc - i;
            pa = a;
            for (k = i; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                vc0 = VLEV_FLOAT(pc0, vl);
                va = VLEV_FLOAT(pa, vl);
                vc0 = VFNMSACVF_FLOAT(vc0, bb0, va, vl);
                VSEV_FLOAT(pc0, vc0, vl);

                pa += vl;
                pc0 += vl;
            }
        }

        a -= m;
        b -= 2 * n;
    }

}
#else
static inline void solve_generic(BLASLONG m, BLASLONG n, FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc) {

  FLOAT aa,  bb;

  int i, j, k;

  a += (m - 1) * m;
  b += (m - 1) * n;

  for (i = m - 1; i >= 0; i--) {

    aa = *(a + i);

    for (j = 0; j < n; j ++) {
      bb = *(c + i + j * ldc);
      bb *= aa;
      *b             = bb;
      *(c + i + j * ldc) = bb;
      b ++;

      for (k = 0; k < i; k ++){
        *(c + k + j * ldc) -= bb * *(a + k);
      }

    }
    a -= m;
    b -= 2 * n;
  }

}

#endif

#else

static inline void solve(BLASLONG m, BLASLONG n, FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc) {

  FLOAT aa1, aa2;
  FLOAT bb1, bb2;
  FLOAT cc1, cc2;

  int i, j, k;

  ldc *= 2;
  a += (m - 1) * m * 2;
  b += (m - 1) * n * 2;

  for (i = m - 1; i >= 0; i--) {

    aa1 = *(a + i * 2 + 0);
    aa2 = *(a + i * 2 + 1);

    for (j = 0; j < n; j ++) {
      bb1 = *(c + i * 2 + 0 + j * ldc);
      bb2 = *(c + i * 2 + 1 + j * ldc);

#ifndef CONJ
      cc1 = aa1 * bb1 - aa2 * bb2;
      cc2 = aa1 * bb2 + aa2 * bb1;
#else
      cc1 = aa1 * bb1 + aa2 * bb2;
      cc2 = aa1 * bb2 - aa2 * bb1;
#endif


      *(b + 0) = cc1;
      *(b + 1) = cc2;
      *(c + i * 2 + 0 + j * ldc) = cc1;
      *(c + i * 2 + 1 + j * ldc) = cc2;
      b += 2;

      for (k = 0; k < i; k ++){
#ifndef CONJ
        *(c + k * 2 + 0 + j * ldc) -= cc1 * *(a + k * 2 + 0) - cc2 * *(a + k * 2 + 1);
        *(c + k * 2 + 1 + j * ldc) -= cc1 * *(a + k * 2 + 1) + cc2 * *(a + k * 2 + 0);
#else
        *(c + k * 2 + 0 + j * ldc) -=   cc1 * *(a + k * 2 + 0) + cc2 * *(a + k * 2 + 1);
        *(c + k * 2 + 1 + j * ldc) -= - cc1 * *(a + k * 2 + 1) + cc2 * *(a + k * 2 + 0);
#endif
      }

    }
    a -= m * 2;
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
  
  size_t vl = VSETVL_MAX;

    //fprintf(stderr, "%s , %s, m = %4ld  n = %4ld  k = %4ld offset = %4ld\n", __FILE__, __FUNCTION__, m, n, k, offset); // Debug

  j = (n >> GEMM_UNROLL_N_SHIFT);

  while (j > 0) {

    kk = m + offset;

    i = m % vl;
    if (i) {
      aa = a + (m - i) * k * COMPSIZE;
      cc = c + (m - i)     * COMPSIZE;

      if (k - kk > 0) {
        GEMM_KERNEL(i, GEMM_UNROLL_N, k - kk, dm1,
#ifdef COMPLEX
            ZERO,
#endif
            aa + i             * kk * COMPSIZE,
            b  + GEMM_UNROLL_N * kk * COMPSIZE,
            cc,
            ldc);
      }

      solve(i, GEMM_UNROLL_N,
          aa + (kk - i) * i             * COMPSIZE,
          b  + (kk - i) * GEMM_UNROLL_N * COMPSIZE,
          cc, ldc);

      kk -= i;

    }

    int mod = i;
    i = vl;
    if (i <= m) {
      aa = a + (m - mod - vl) * k * COMPSIZE;
      cc = c + (m - mod - vl)     * COMPSIZE;

      do {
        if (k - kk > 0) {
          GEMM_KERNEL(vl, GEMM_UNROLL_N, k - kk, dm1,
#ifdef COMPLEX
              ZERO,
#endif
              aa + vl * kk * COMPSIZE,
              b +  GEMM_UNROLL_N * kk * COMPSIZE,
              cc,
              ldc);
        }

        solve(vl, GEMM_UNROLL_N,
            aa + (kk - vl) * vl * COMPSIZE,
            b  + (kk - vl) * GEMM_UNROLL_N * COMPSIZE,
            cc, ldc);

        aa -= vl * k * COMPSIZE;
        cc -= vl     * COMPSIZE;
        kk -= vl;

        i += vl;
      } while (i <= m);
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

        i = m % vl;
        if (i) {
          aa = a + (m - i) * k * COMPSIZE;
          cc = c + (m - i)     * COMPSIZE;

          if (k - kk > 0) {
            GEMM_KERNEL(i, j, k - kk, dm1,
#ifdef COMPLEX
                ZERO,
#endif
                aa + i * kk * COMPSIZE,
                b  + j * kk * COMPSIZE,
                cc, ldc);
          }

          solve(i, j,
              aa + (kk - i) * i * COMPSIZE,
              b  + (kk - i) * j * COMPSIZE,
              cc, ldc);

          kk -= i;

        }

        int mod = i;
        i = vl;
        if (i <= m) {
          aa = a + (m - mod - vl) * k * COMPSIZE;
          cc = c + (m - mod - vl)     * COMPSIZE;

          do {
            if (k - kk > 0) {
              GEMM_KERNEL(vl, j, k - kk, dm1,
#ifdef COMPLEX
                  ZERO,
#endif
                  aa + vl * kk * COMPSIZE,
                  b +  j             * kk * COMPSIZE,
                  cc,
                  ldc);
            }

            solve(vl, j,
                aa + (kk - vl) * vl * COMPSIZE,
                b  + (kk - vl) * j             * COMPSIZE,
                cc, ldc);

            aa -= vl * k * COMPSIZE;
            cc -= vl     * COMPSIZE;
            kk -= vl;

            i += vl;
          } while (i <= m);
        }

        b += j * k   * COMPSIZE;
        c += j * ldc * COMPSIZE;
      }
      j >>= 1;
    }
  }

  return 0;
}
