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
#define FLOAT_V_T vfloat32m2_t
#define VLSEG2_FLOAT vlseg2e32_v_f32m2
#define VSSSEG2_FLOAT vssseg2e32_v_f32m2
#define VSSSEG4_FLOAT vssseg4e32_v_f32m2
#else
#define VSETVL(n) vsetvl_e64m2(n)
#define FLOAT_V_T vfloat64m2_t
#define VLSEG2_FLOAT vlseg2e64_v_f64m2
#define VSSSEG2_FLOAT vssseg2e64_v_f64m2
#define VSSSEG4_FLOAT vssseg4e64_v_f64m2
#endif

// Optimizes the implementation in ../generic/gemm_tcopy_4.c

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, FLOAT *b)
{
    BLASLONG i, j;

    FLOAT *a_offset, *a_offset1, *a_offset2, *a_offset3, *a_offset4;
    FLOAT *b_offset, *b_offset1, *b_offset2, *b_offset3;
    FLOAT  ctemp1,  ctemp2,  ctemp3,  ctemp4;
    FLOAT  ctemp5,  ctemp6,  ctemp7,  ctemp8;
    FLOAT  ctemp9, ctemp10, ctemp11, ctemp12;
    FLOAT ctemp13, ctemp14, ctemp15, ctemp16;

    //fprintf(stderr, "gemm_tcopy_4 m=%ld n=%ld lda=%ld\n", m, n, lda);

    a_offset = a;
    b_offset = b;

    b_offset2 = b + m  * (n & ~3);
    b_offset3 = b + m  * (n & ~1);

    for(j = (m >> 2); j > 0; j--) {
        a_offset1 = a_offset;
        a_offset2 = a_offset1 + lda;
        a_offset3 = a_offset2 + lda;
        a_offset4 = a_offset3 + lda;
        a_offset += 4 * lda;

        b_offset1 = b_offset;
        b_offset += 16;

        for(i = (n >> 2); i > 0; i--) {
            v1 = VLEV_FLOAT(a_offset1, 4);
            v2 = VLEV_FLOAT(a_offset2, 4);
            v3 = VLEV_FLOAT(a_offset3, 4);
            v4 = VLEV_FLOAT(a_offset4, 4);

            a_offset1 += 4;
            a_offset2 += 4;
            a_offset3 += 4;
            a_offset4 += 4;

            VSEV_FLOAT(b_offset1, v1, 4);
            VSEV_FLOAT(b_offset2+4, v2, 4);
            VSEV_FLOAT(b_offset2+8, v3, 4);
            VSEV_FLOAT(b_offset2+12, v4, 4);

            b_offset1 += m * 4;
        }

        if (n & 2) {
            v1 = VLEV_FLOAT(a_offset1, 2);
            v2 = VLEV_FLOAT(a_offset2, 2);
            v3 = VLEV_FLOAT(a_offset3, 2);
            v4 = VLEV_FLOAT(a_offset4, 2);

            a_offset1 += 2;
            a_offset2 += 2;
            a_offset3 += 2;
            a_offset4 += 2;

            VSEV_FLOAT(b_offset2, v1, 2);
            VSEV_FLOAT(b_offset2+2, v2, 2);
            VSEV_FLOAT(b_offset2+4, v3, 2);
            VSEV_FLOAT(b_offset2+6, v4, 2);

            b_offset2 += 8;
        }

        if (n & 1) {
            v1 = VLEV_FLOAT(a_offset1, 1);
            v2 = VLEV_FLOAT(a_offset2, 1);
            v3 = VLEV_FLOAT(a_offset3, 1);
            v4 = VLEV_FLOAT(a_offset4, 1);

            VSSEG4_FLOAT(b_offset3, v1, v2, v3, v4, 1);

            b_offset3 += 4;
        }

    }

// TODO cleanup

  if (m & 2){
    a_offset1  = a_offset;
    a_offset2  = a_offset1 + lda;
    a_offset  += 2 * lda;

    b_offset1  = b_offset;
    b_offset  += 8;

    i = (n >> 2);
    if (i > 0){
      do{
	ctemp1  = *(a_offset1 + 0);
	ctemp2  = *(a_offset1 + 1);
	ctemp3  = *(a_offset1 + 2);
	ctemp4  = *(a_offset1 + 3);

	ctemp5  = *(a_offset2 + 0);
	ctemp6  = *(a_offset2 + 1);
	ctemp7  = *(a_offset2 + 2);
	ctemp8  = *(a_offset2 + 3);

	a_offset1 += 4;
	a_offset2 += 4;

	*(b_offset1 +  0) = ctemp1;
	*(b_offset1 +  1) = ctemp2;
	*(b_offset1 +  2) = ctemp3;
	*(b_offset1 +  3) = ctemp4;

	*(b_offset1 +  4) = ctemp5;
	*(b_offset1 +  5) = ctemp6;
	*(b_offset1 +  6) = ctemp7;
	*(b_offset1 +  7) = ctemp8;

	b_offset1 += m * 4;
	i --;
      }while(i > 0);
    }

    if (n & 2) {
      ctemp1  = *(a_offset1 + 0);
      ctemp2  = *(a_offset1 + 1);

      ctemp3  = *(a_offset2 + 0);
      ctemp4  = *(a_offset2 + 1);

      a_offset1 += 2;
      a_offset2 += 2;

      *(b_offset2 +  0) = ctemp1;
      *(b_offset2 +  1) = ctemp2;
      *(b_offset2 +  2) = ctemp3;
      *(b_offset2 +  3) = ctemp4;

      b_offset2 += 4;
    }

    if (n & 1) {
      ctemp1  = *(a_offset1 + 0);
      ctemp2  = *(a_offset2 + 0);

      *(b_offset3 +  0) = ctemp1;
      *(b_offset3 +  1) = ctemp2;
      b_offset3 += 2;
    }
  }

  if (m & 1){
    a_offset1  = a_offset;
    b_offset1  = b_offset;

    i = (n >> 2);
    if (i > 0){
      do{
	ctemp1  = *(a_offset1 + 0);
	ctemp2  = *(a_offset1 + 1);
	ctemp3  = *(a_offset1 + 2);
	ctemp4  = *(a_offset1 + 3);

	a_offset1 += 4;

	*(b_offset1 +  0) = ctemp1;
	*(b_offset1 +  1) = ctemp2;
	*(b_offset1 +  2) = ctemp3;
	*(b_offset1 +  3) = ctemp4;

	b_offset1 += 4 * m;

	i --;
      }while(i > 0);
    }

    if (n & 2) {
      ctemp1  = *(a_offset1 + 0);
      ctemp2  = *(a_offset1 + 1);
      a_offset1 += 2;

      *(b_offset2 +  0) = ctemp1;
      *(b_offset2 +  1) = ctemp2;
    }

    if (n & 1) {
      ctemp1  = *(a_offset1 + 0);
      *(b_offset3 +  0) = ctemp1;
    }
  }

  return 0;
}
