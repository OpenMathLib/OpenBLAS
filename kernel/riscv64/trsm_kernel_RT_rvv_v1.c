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
#define VFNMSACVF_FLOAT vfnmsac_vf_f64m2
#endif


static FLOAT dm1 = -1.;

#ifdef CONJ
#define GEMM_KERNEL   GEMM_KERNEL_R
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

// Optimizes the implementation in ../arm64/trsm_kernel_RT_sve.c

#ifndef COMPLEX

#if GEMM_DEFAULT_UNROLL_N == 1
static inline void solve(BLASLONG m, BLASLONG n, FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc) {

  FLOAT aa,  bb;
  FLOAT *pb, *pc;
  BLASLONG stride_ldc = sizeof(FLOAT) * ldc;

  int i, j, k;
  size_t vl;
  FLOAT_V_T vb, vc;

  a += (n - 1) * m;
  b += (n - 1) * n;

  for (i = n - 1; i >= 0; i--) {

    bb = *(b + i);

    for (j = 0; j < m; j ++) {
      aa = *(c + j + i * ldc);
      aa *= bb;
      *a   = aa;
      *(c + j + i * ldc) = aa;
      a ++;

        pb = b;
        pc = c + j;
        for (k = i; k > 0; k -= vl)
        {
            vl = VSETVL(k);
            vc = VLSEV_FLOAT(pc, stride_ldc, vl);
            vb = VLEV_FLOAT(pb, vl);
            vc = VFNMSACVF_FLOAT(vc, aa, vb, vl);
            VSSEV_FLOAT(pc, stride_ldc, vc, vl);
            pb += vl;
            pc++;
        }
    }
    b -= n;
    a -= 2 * m;
  }

}
#elif GEMM_DEFAULT_UNROLL_N == 2

static inline void solve(BLASLONG m, BLASLONG n, FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc) {

    FLOAT aa0, aa1, bb;
    FLOAT *pb, *pc;
    FLOAT *pa0, *pa1, *pc0, *pc1;
    BLASLONG stride_ldc = sizeof(FLOAT) * ldc;
    int i, j, k;
    size_t vl;
    FLOAT_V_T vb, vc0, vc1;

    a += (n - 1) * m;
    b += (n - 1) * n;

    for (i = n - 1; i >= 0; i--)
    {
        bb = *(b + i);
        pc = c + i * ldc;
        for (j = 0; j < m/2; j ++) 
        {
            pa0 = pc + j * 2;
            pa1 = pc + j * 2 + 1;
            aa0 = *pa0 * bb;
            aa1 = *pa1 * bb;

            *pa0    = aa0;
            *pa1    = aa1;
            *a      = aa0;
            *(a + 1)= aa1;
            a  += 2;

            pb  = b;
            pc0 = c + j * 2;
            pc1 = pc0 + 1;
            for (k = i; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                vc0 = VLSEV_FLOAT(pc0, stride_ldc, vl);
                vc1 = VLSEV_FLOAT(pc1, stride_ldc, vl);
                vb = VLEV_FLOAT(pb, vl);
                vc0 = VFNMSACVF_FLOAT(vc0, aa0, vb, vl);
                vc1 = VFNMSACVF_FLOAT(vc1, aa1, vb, vl);
                VSSEV_FLOAT(pc0, stride_ldc, vc0, vl);
                VSSEV_FLOAT(pc1, stride_ldc, vc1, vl);
                pb += vl;
                pc0++;
                pc1++;
            }
        }
        pc += (m/2)*2;

        if (m & 1)
        {
            pa0 = pc;
            aa0 = *pa0 * bb;
            
            *pa0    = aa0;
            *a      = aa0;
            a  += 1;
           
            pb = b;
            pc0 = pc - i * ldc;
            for (k = i; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                vc0 = VLSEV_FLOAT(pc0, stride_ldc, vl);
                vb = VLEV_FLOAT(pb, vl);
                vc0 = VFNMSACVF_FLOAT(vc0, aa0, vb, vl);
                VSSEV_FLOAT(pc0, stride_ldc, vc0, vl);
                pb += vl;
                pc0++;
            }
        }
        b -= n;
        a -= 2 * m;
    }
}

#elif GEMM_DEFAULT_UNROLL_N == 4

static inline void solve(BLASLONG m, BLASLONG n, FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc) {

    FLOAT aa0, aa1, aa2, aa3;
    FLOAT bb;
    FLOAT *pb, *pc;
    FLOAT *pa0, *pa1, *pa2, *pa3;
    FLOAT *pc0, *pc1, *pc2, *pc3;
    BLASLONG stride_ldc = sizeof(FLOAT) * ldc;
    int i, j, k;
    size_t vl;
    FLOAT_V_T vb, vc0, vc1, vc2, vc3;

    a += (n - 1) * m;
    b += (n - 1) * n;

    for (i = n - 1; i >= 0; i--)
    {
        bb = *(b + i);
        pc = c + i * ldc;
        for (j = 0; j < m/4; j ++) 
        {
            pa0 = pc + j * 4;
            pa1 = pa0 + 1;
            pa2 = pa1 + 1;
            pa3 = pa2 + 1;
            
            aa0 = *pa0 * bb;
            aa1 = *pa1 * bb;
            aa2 = *pa2 * bb;
            aa3 = *pa3 * bb;

            *pa0    = aa0;
            *pa1    = aa1;
            *pa2    = aa2;
            *pa3    = aa3;

            *a      = aa0;
            *(a + 1)= aa1;
            *(a + 2)= aa2;
            *(a + 3)= aa3;
            a  += 4;

            pb  = b;
            pc0 = c + j * 4;
            pc1 = pc0 + 1;
            pc2 = pc1 + 1;
            pc3 = pc2 + 1;
            for (k = i; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                vc0 = VLSEV_FLOAT(pc0, stride_ldc, vl);
                vc1 = VLSEV_FLOAT(pc1, stride_ldc, vl);
                vc2 = VLSEV_FLOAT(pc2, stride_ldc, vl);
                vc3 = VLSEV_FLOAT(pc3, stride_ldc, vl);
                vb = VLEV_FLOAT(pb, vl);

                vc0 = VFNMSACVF_FLOAT(vc0, aa0, vb, vl);
                vc1 = VFNMSACVF_FLOAT(vc1, aa1, vb, vl);
                vc2 = VFNMSACVF_FLOAT(vc2, aa2, vb, vl);
                vc3 = VFNMSACVF_FLOAT(vc3, aa3, vb, vl);

                VSSEV_FLOAT(pc0, stride_ldc, vc0, vl);
                VSSEV_FLOAT(pc1, stride_ldc, vc1, vl);
                VSSEV_FLOAT(pc2, stride_ldc, vc2, vl);
                VSSEV_FLOAT(pc3, stride_ldc, vc3, vl);

                pb += vl;
                pc0++;
                pc1++;
                pc2++;
                pc3++;
            }
        }
        pc += (m/4)*4;

        if (m & 2)
        {
            pa0 = pc + j * 2;
            pa1 = pa0 + 1;
            
            aa0 = *pa0 * bb;
            aa1 = *pa1 * bb;

            *pa0    = aa0;
            *pa1    = aa1;

            *a      = aa0;
            *(a + 1)= aa1;
            a  += 2;

            pb  = b;
            pc0 = c + j * 4;
            pc1 = pc0 + 1;
            for (k = i; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                vc0 = VLSEV_FLOAT(pc0, stride_ldc, vl);
                vc1 = VLSEV_FLOAT(pc1, stride_ldc, vl);
                vb = VLEV_FLOAT(pb, vl);

                vc0 = VFNMSACVF_FLOAT(vc0, aa0, vb, vl);
                vc1 = VFNMSACVF_FLOAT(vc1, aa1, vb, vl);

                VSSEV_FLOAT(pc0, stride_ldc, vc0, vl);
                VSSEV_FLOAT(pc1, stride_ldc, vc1, vl);

                pb += vl;
                pc0++;
                pc1++;
            }
            pc += 2;
        }

        if (m & 1)
        {
            pa0 = pc;
            aa0 = *pa0 * bb;
            
            *pa0    = aa0;
            *a      = aa0;
            a  += 1;
           
            pb = b;
            pc0 = pc - i * ldc;
            for (k = i; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                vc0 = VLSEV_FLOAT(pc0, stride_ldc, vl);
                vb = VLEV_FLOAT(pb, vl);
                vc0 = VFNMSACVF_FLOAT(vc0, aa0, vb, vl);
                VSSEV_FLOAT(pc0, stride_ldc, vc0, vl);
                pb += vl;
                pc0++;
            }
        }
        b -= n;
        a -= 2 * m;
    }
}
#elif GEMM_DEFAULT_UNROLL_N == 8

static inline void solve(BLASLONG m, BLASLONG n, FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc) {

    FLOAT aa0, aa1, aa2, aa3, aa4, aa5, aa6, aa7;
    FLOAT bb;
    FLOAT *pb, *pc;
    FLOAT *pa0, *pa1, *pa2, *pa3, *pa4, *pa5, *pa6, *pa7;
    FLOAT *pc0, *pc1, *pc2, *pc3, *pc4, *pc5, *pc6, *pc7;
    BLASLONG stride_ldc = sizeof(FLOAT) * ldc;
    int i, j, k;
    size_t vl;
    FLOAT_V_T vb, vc0, vc1, vc2, vc3, vc4, vc5, vc6, vc7;

    a += (n - 1) * m;
    b += (n - 1) * n;

    for (i = n - 1; i >= 0; i--)
    {
        bb = *(b + i);
        pc = c + i * ldc;
        for (j = 0; j < m/8; j ++) 
        {
            pa0 = pc + j * 8;
            pa1 = pa0 + 1;
            pa2 = pa1 + 1;
            pa3 = pa2 + 1;
            pa4 = pa3 + 1;
            pa5 = pa4 + 1;
            pa6 = pa5 + 1;
            pa7 = pa6 + 1;
            
            aa0 = *pa0 * bb;
            aa1 = *pa1 * bb;
            aa2 = *pa2 * bb;
            aa3 = *pa3 * bb;
            aa4 = *pa4 * bb;
            aa5 = *pa5 * bb;
            aa6 = *pa6 * bb;
            aa7 = *pa7 * bb;

            *pa0    = aa0;
            *pa1    = aa1;
            *pa2    = aa2;
            *pa3    = aa3;
            *pa4    = aa4;
            *pa5    = aa5;
            *pa6    = aa6;
            *pa7    = aa7;

            *a      = aa0;
            *(a + 1)= aa1;
            *(a + 2)= aa2;
            *(a + 3)= aa3;
            *(a + 4)= aa4;
            *(a + 5)= aa5;
            *(a + 6)= aa6;
            *(a + 7)= aa7;
            a  += 8;

            pb  = b;
            pc0 = c + j * 8;
            pc1 = pc0 + 1;
            pc2 = pc1 + 1;
            pc3 = pc2 + 1;
            pc4 = pc3 + 1;
            pc5 = pc4 + 1;
            pc6 = pc5 + 1;
            pc7 = pc6 + 1;
            for (k = i; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                vc0 = VLSEV_FLOAT(pc0, stride_ldc, vl);
                vc1 = VLSEV_FLOAT(pc1, stride_ldc, vl);
                vc2 = VLSEV_FLOAT(pc2, stride_ldc, vl);
                vc3 = VLSEV_FLOAT(pc3, stride_ldc, vl);
                vc4 = VLSEV_FLOAT(pc4, stride_ldc, vl);
                vc5 = VLSEV_FLOAT(pc5, stride_ldc, vl);
                vc6 = VLSEV_FLOAT(pc6, stride_ldc, vl);
                vc7 = VLSEV_FLOAT(pc7, stride_ldc, vl);
                vb = VLEV_FLOAT(pb, vl);

                vc0 = VFNMSACVF_FLOAT(vc0, aa0, vb, vl);
                vc1 = VFNMSACVF_FLOAT(vc1, aa1, vb, vl);
                vc2 = VFNMSACVF_FLOAT(vc2, aa2, vb, vl);
                vc3 = VFNMSACVF_FLOAT(vc3, aa3, vb, vl);
                vc4 = VFNMSACVF_FLOAT(vc4, aa4, vb, vl);
                vc5 = VFNMSACVF_FLOAT(vc5, aa5, vb, vl);
                vc6 = VFNMSACVF_FLOAT(vc6, aa6, vb, vl);
                vc7 = VFNMSACVF_FLOAT(vc7, aa7, vb, vl);

                VSSEV_FLOAT(pc0, stride_ldc, vc0, vl);
                VSSEV_FLOAT(pc1, stride_ldc, vc1, vl);
                VSSEV_FLOAT(pc2, stride_ldc, vc2, vl);
                VSSEV_FLOAT(pc3, stride_ldc, vc3, vl);
                VSSEV_FLOAT(pc4, stride_ldc, vc4, vl);
                VSSEV_FLOAT(pc5, stride_ldc, vc5, vl);
                VSSEV_FLOAT(pc6, stride_ldc, vc6, vl);
                VSSEV_FLOAT(pc7, stride_ldc, vc7, vl);

                pb += vl;
                pc0++;
                pc1++;
                pc2++;
                pc3++;
                pc4++;
                pc5++;
                pc6++;
                pc7++;
            }
        }
        pc += (m/8)*8;

        if (m & 4)
        {
            pa0 = pc;
            pa1 = pa0 + 1;
            pa2 = pa1 + 1;
            pa3 = pa2 + 1;
            
            aa0 = *pa0 * bb;
            aa1 = *pa1 * bb;
            aa2 = *pa2 * bb;
            aa3 = *pa3 * bb;

            *pa0    = aa0;
            *pa1    = aa1;
            *pa2    = aa2;
            *pa3    = aa3;

            *a      = aa0;
            *(a + 1)= aa1;
            *(a + 2)= aa2;
            *(a + 3)= aa3;
            a  += 4;

            pb  = b;
            pc0 = pc - i * ldc;
            pc1 = pc0 + 1;
            pc2 = pc1 + 1;
            pc3 = pc2 + 1;
            for (k = i; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                vc0 = VLSEV_FLOAT(pc0, stride_ldc, vl);
                vc1 = VLSEV_FLOAT(pc1, stride_ldc, vl);
                vc2 = VLSEV_FLOAT(pc2, stride_ldc, vl);
                vc3 = VLSEV_FLOAT(pc3, stride_ldc, vl);
                vb = VLEV_FLOAT(pb, vl);

                vc0 = VFNMSACVF_FLOAT(vc0, aa0, vb, vl);
                vc1 = VFNMSACVF_FLOAT(vc1, aa1, vb, vl);
                vc2 = VFNMSACVF_FLOAT(vc2, aa2, vb, vl);
                vc3 = VFNMSACVF_FLOAT(vc3, aa3, vb, vl);

                VSSEV_FLOAT(pc0, stride_ldc, vc0, vl);
                VSSEV_FLOAT(pc1, stride_ldc, vc1, vl);
                VSSEV_FLOAT(pc2, stride_ldc, vc2, vl);
                VSSEV_FLOAT(pc3, stride_ldc, vc3, vl);

                pb += vl;
                pc0++;
                pc1++;
                pc2++;
                pc3++;
            }
            pc += 4;
        }

        if (m & 2)
        {
            pa0 = pc;
            pa1 = pa0 + 1;
            
            aa0 = *pa0 * bb;
            aa1 = *pa1 * bb;

            *pa0    = aa0;
            *pa1    = aa1;

            *a      = aa0;
            *(a + 1)= aa1;
            a  += 2;

            pb  = b;
            pc0 = pc - i * ldc;
            pc1 = pc0 + 1;
            for (k = i; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                vc0 = VLSEV_FLOAT(pc0, stride_ldc, vl);
                vc1 = VLSEV_FLOAT(pc1, stride_ldc, vl);
                vb = VLEV_FLOAT(pb, vl);

                vc0 = VFNMSACVF_FLOAT(vc0, aa0, vb, vl);
                vc1 = VFNMSACVF_FLOAT(vc1, aa1, vb, vl);

                VSSEV_FLOAT(pc0, stride_ldc, vc0, vl);
                VSSEV_FLOAT(pc1, stride_ldc, vc1, vl);

                pb += vl;
                pc0++;
                pc1++;
            }
            pc += 2;
        }

        if (m & 1)
        {
            pa0 = pc;
            aa0 = *pa0 * bb;
            
            *pa0    = aa0;
            *a      = aa0;
            a  += 1;
           
            pb = b;
            pc0 = pc - i * ldc;
            for (k = i; k > 0; k -= vl)
            {
                vl = VSETVL(k);
                vc0 = VLSEV_FLOAT(pc0, stride_ldc, vl);
                vb = VLEV_FLOAT(pb, vl);
                vc0 = VFNMSACVF_FLOAT(vc0, aa0, vb, vl);
                VSSEV_FLOAT(pc0, stride_ldc, vc0, vl);
                pb += vl;
                pc0++;
            }
        }
        b -= n;
        a -= 2 * m;
    }
}

#else

static inline void solve(BLASLONG m, BLASLONG n, FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc) {

  FLOAT aa,  bb;

  int i, j, k;

  a += (n - 1) * m;
  b += (n - 1) * n;

  for (i = n - 1; i >= 0; i--) {

    bb = *(b + i);

    for (j = 0; j < m; j ++) {
      aa = *(c + j + i * ldc);
      aa *= bb;
      *a   = aa;
      *(c + j + i * ldc) = aa;
      a ++;

      for (k = 0; k < i; k ++){
	*(c + j + k * ldc) -= aa * *(b + k);
      }

    }
    b -= n;
    a -= 2 * m;
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

  a += (n - 1) * m * 2;
  b += (n - 1) * n * 2;

  for (i = n - 1; i >= 0; i--) {

    bb1 = *(b + i * 2 + 0);
    bb2 = *(b + i * 2 + 1);

    for (j = 0; j < m; j ++) {

      aa1 = *(c + j * 2 + 0 + i * ldc);
      aa2 = *(c + j * 2 + 1 + i * ldc);

#ifndef CONJ
      cc1 = aa1 * bb1 - aa2 * bb2;
      cc2 = aa1 * bb2 + aa2 * bb1;
#else
      cc1 =  aa1 * bb1  + aa2 * bb2;
      cc2 = - aa1 * bb2 + aa2 * bb1;
#endif

      *(a + 0) = cc1;
      *(a + 1) = cc2;

      *(c + j * 2 + 0 + i * ldc) = cc1;
      *(c + j * 2 + 1 + i * ldc) = cc2;
      a += 2;

      for (k = 0; k < i; k ++){
#ifndef CONJ
	*(c + j * 2 + 0 + k * ldc) -= cc1 * *(b + k * 2 + 0) - cc2 * *(b + k * 2 + 1);
	*(c + j * 2 + 1 + k * ldc) -= cc1 * *(b + k * 2 + 1) + cc2 * *(b + k * 2 + 0);
#else
	*(c + j * 2 + 0 + k * ldc) -=   cc1 * *(b + k * 2 + 0) + cc2 * *(b + k * 2 + 1);
	*(c + j * 2 + 1 + k * ldc) -=  -cc1 * *(b + k * 2 + 1) + cc2 * *(b + k * 2 + 0);
#endif
      }

    }
    b -= n * 2;
    a -= 4 * m;
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

  kk = n - offset;
  c += n * ldc * COMPSIZE;
  b += n * k   * COMPSIZE;

  if (n & (GEMM_UNROLL_N - 1)) {

    j = 1;
    while (j < GEMM_UNROLL_N) {
      if (n & j) {

        aa  = a;
        b -= j * k  * COMPSIZE;
        c -= j * ldc* COMPSIZE;
        cc  = c;

        i = vl;
        if (i <= m) {

          do {
            if (k - kk > 0) {
              GEMM_KERNEL(vl, j, k - kk, dm1,
#ifdef COMPLEX
                  ZERO,
#endif
                  aa + vl * kk * COMPSIZE,
                  b  +  j            * kk * COMPSIZE,
                  cc,
                  ldc);
            }

            solve(vl, j,
                aa + (kk - j) * vl * COMPSIZE,
                b  + (kk - j) * j             * COMPSIZE,
                cc, ldc);

            aa += vl * k * COMPSIZE;
            cc += vl     * COMPSIZE;
            i += vl;
          } while (i <= m);
        }

        i = m % vl;
        if (i) {
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
              aa + (kk - j) * i * COMPSIZE,
              b  + (kk - j) * j * COMPSIZE,
              cc, ldc);

          aa += i * k * COMPSIZE;
          cc += i     * COMPSIZE;

        }
        kk -= j;
      }
      j <<= 1;
    }
  }

  j = (n >> GEMM_UNROLL_N_SHIFT);

  if (j > 0) {

    do {
      aa  = a;
      b -= GEMM_UNROLL_N * k   * COMPSIZE;
      c -= GEMM_UNROLL_N * ldc * COMPSIZE;
      cc  = c;

      i = vl;
      if (i <= m) {
	do {
	  if (k - kk > 0) {
	    GEMM_KERNEL(vl, GEMM_UNROLL_N, k - kk, dm1,
#ifdef COMPLEX
			ZERO,
#endif
			aa + vl * kk * COMPSIZE,
			b  + GEMM_UNROLL_N * kk * COMPSIZE,
			cc,
			ldc);
	  }

	  solve(vl, GEMM_UNROLL_N,
		aa + (kk - GEMM_UNROLL_N) * vl * COMPSIZE,
		b  + (kk - GEMM_UNROLL_N) * GEMM_UNROLL_N * COMPSIZE,
		cc, ldc);

	  aa += vl * k * COMPSIZE;
	  cc += vl     * COMPSIZE;
	  i += vl;
	} while (i <= m);
      }

      i = m % vl;
      if (i) {
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
		  aa + (kk - GEMM_UNROLL_N) * i             * COMPSIZE,
		  b  + (kk - GEMM_UNROLL_N) * GEMM_UNROLL_N * COMPSIZE,
		  cc, ldc);

	    aa += i * k * COMPSIZE;
	    cc += i     * COMPSIZE;

      }

      kk -= GEMM_UNROLL_N;
      j --;
    } while (j > 0);
  }

  return 0;
}


