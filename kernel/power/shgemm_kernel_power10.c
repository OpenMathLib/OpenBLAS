/*********************************************************************************
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
**********************************************************************************/
#include "common.h"
#include <altivec.h>
#if defined(HALF) && defined(HALFCONVERSION)
static float
bfloat16tof32 (bfloat16 f16)
{
  float result = 0;
  unsigned short *q = (unsigned short *) (&result);
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
  q[0] = f16;
#else
  q[1] = f16;
#endif
  return result;
}

#define BF16TOF32(x) (bfloat16tof32(x))
#else
#define BF16TOF32(x) x
#endif

typedef unsigned char vec_t __attribute__ ((vector_size (16)));
typedef FLOAT v4sf_t __attribute__ ((vector_size (16)));
typedef FLOAT v2sf_t __attribute__ ((vector_size (8)));

vector char mask =
  { 0x0, 0x1, 0x8, 0x9, 0x2, 0x3, 0xa, 0xb, 0x4, 0x5, 0xc, 0xd, 0x6, 0x7, 0xe,
  0xf
};

/* 
 * BFLOAT16 xvbf16ger2pp instruction needs 4×2 matrix of
 * bfloat16 floating-point values as input. Hence this
 * merging is needed on A and B matrices. 
 */
#define MERGE_ROW(x) vec_perm(x, x, mask)
#define MERGE_HIGH(x, y) (vec_t) vec_mergeh ((vector short)x, (vector short)y)
#define MERGE_LOW(x, y) (vec_t) vec_mergel ((vector short)x, (vector short)y)

#define SAVE_ACC(ACC, J)  \
	  __builtin_mma_disassemble_acc (result, ACC); \
	  rowC = (v4sf_t *) &CO[0* ldc+J]; \
          rowC[0] += result[3] * alpha; \
          rowC = (v4sf_t *) &CO[1*ldc+J]; \
          rowC[0] += result[2] * alpha; \
          rowC = (v4sf_t *) &CO[2*ldc+J]; \
          rowC[0] += result[1] * alpha; \
          rowC = (v4sf_t *) &CO[3*ldc+J]; \
          rowC[0] += result[0] * alpha;
#define SAVE_ACC1(ACC, J)  \
	  __builtin_mma_disassemble_acc (result, ACC); \
	  rowC = (v4sf_t *) &CO[4* ldc+J]; \
          rowC[0] += result[3] * alpha; \
          rowC = (v4sf_t *) &CO[5*ldc+J]; \
          rowC[0] += result[2] * alpha; \
          rowC = (v4sf_t *) &CO[6*ldc+J]; \
          rowC[0] += result[1] * alpha; \
          rowC = (v4sf_t *) &CO[7*ldc+J]; \
          rowC[0] += result[0] * alpha;
#define  SAVE4x2_ACC(ACC, J)  \
	  __builtin_mma_disassemble_acc (result, ACC); \
	  rowC = (v2sf_t *) &CO[0* ldc+J]; \
          rowC[0] += result[6] * alpha; \
	  rowC = (v2sf_t *) &CO[1* ldc+J]; \
          rowC[0] += result[4] * alpha; \
	  rowC = (v2sf_t *) &CO[2* ldc+J]; \
          rowC[0] += result[2] * alpha; \
	  rowC = (v2sf_t *) &CO[3* ldc+J]; \
          rowC[0] += result[0] * alpha;
#define  SAVE4x2_ACC1(ACC, J)  \
	  __builtin_mma_disassemble_acc (result, ACC); \
	  rowC = (v2sf_t *) &CO[4* ldc+J]; \
          rowC[0] += result[6] * alpha; \
	  rowC = (v2sf_t *) &CO[5* ldc+J]; \
          rowC[0] += result[4] * alpha; \
	  rowC = (v2sf_t *) &CO[6* ldc+J]; \
          rowC[0] += result[2] * alpha; \
	  rowC = (v2sf_t *) &CO[7* ldc+J]; \
          rowC[0] += result[0] * alpha;

#define MMA __builtin_mma_xvbf16ger2pp

#define  SAVE2x4_ACC(ACC, J)  \
	  __builtin_mma_disassemble_acc (result, ACC); \
	  rowC = (v4sf_t *) &CO[0* ldc+J]; \
          rowC[0] += result[3] * alpha; \
	  rowC = (v4sf_t *) &CO[1* ldc+J]; \
          rowC[0] += result[2] * alpha;

#define SET_ACC_ZERO4() \
	  __builtin_mma_xxsetaccz (&acc0); \
	  __builtin_mma_xxsetaccz (&acc1); \
	  __builtin_mma_xxsetaccz (&acc2); \
	  __builtin_mma_xxsetaccz (&acc3);

#define SET_ACC_ZERO8() \
	  __builtin_mma_xxsetaccz (&acc0); \
	  __builtin_mma_xxsetaccz (&acc1); \
	  __builtin_mma_xxsetaccz (&acc2); \
	  __builtin_mma_xxsetaccz (&acc3); \
	  __builtin_mma_xxsetaccz (&acc4); \
	  __builtin_mma_xxsetaccz (&acc5); \
	  __builtin_mma_xxsetaccz (&acc6); \
	  __builtin_mma_xxsetaccz (&acc7);

#define PREFETCH1(x, y) asm volatile ("dcbt %0, %1" : : "r" (x), "b" (y) : "memory");
/*************************************************************************************
* SHGEMM Kernel
*************************************************************************************/
int
CNAME (BLASLONG m, BLASLONG n, BLASLONG k, FLOAT alpha, IFLOAT * A,
       IFLOAT * B, FLOAT * C, BLASLONG ldc)
{
  BLASLONG N = n;
  BLASLONG i1;
  v4sf_t valpha = { alpha, alpha, alpha, alpha };
  vector short vzero = { 0, 0, 0, 0, 0, 0, 0, 0 };
  N = n >> 3;
  /* Loop for n >= 8. */
  for (i1 = 0; i1 < N; i1++)
    {
      BLASLONG i, j;
      FLOAT *CO;
      IFLOAT *AO;
      CO = C;
      C += ldc << 3;
      AO = A;
      PREFETCH1 (A, 128);
      PREFETCH1 (A, 256);
      i = m >> 4;
      /* Loop for m >= 16. */
      for (j = 0; j < i; j++)
	{
	  IFLOAT *BO = B;
	  v4sf_t *rowC;
	  v4sf_t result[4];
	  __vector_quad acc0, acc1, acc2, acc3, acc4, acc5, acc6, acc7;
	  SET_ACC_ZERO8 ();
	  BLASLONG l = 0;
	  for (l = 0; l < k / 2; l++)
	    {
	      vec_t *rowA = (vec_t *) & (AO[l << 5]);
	      vec_t *rowB = (vec_t *) & (BO[l << 4]);
	      vec_t rowB_h = MERGE_HIGH (rowB[0], rowB[1]);
	      vec_t rowB_l = MERGE_LOW (rowB[0], rowB[1]);
	      vec_t rowA_h = MERGE_HIGH (rowA[0], rowA[2]);
	      vec_t rowA_l = MERGE_LOW (rowA[0], rowA[2]);
	      vec_t rowA2_h = MERGE_HIGH (rowA[1], rowA[3]);
	      vec_t rowA2_l = MERGE_LOW (rowA[1], rowA[3]);
	      MMA (&acc0, rowB_h, rowA_h);
	      MMA (&acc1, rowB_l, rowA_h);
	      MMA (&acc2, rowB_h, rowA_l);
	      MMA (&acc3, rowB_l, rowA_l);
	      MMA (&acc4, rowB_h, rowA2_h);
	      MMA (&acc5, rowB_l, rowA2_h);
	      MMA (&acc6, rowB_h, rowA2_l);
	      MMA (&acc7, rowB_l, rowA2_l);
	    }
	  if (k % 2 == 1)
	    {
	      if (k > 1)
		l = (k / 2) << 4;
	      vec_t *rowA = (vec_t *) & (AO[l << 1]);
	      vec_t *rowB = (vec_t *) & (BO[l]);
	      vec_t rowB_h = MERGE_HIGH (rowB[0], rowB[1]);
	      vec_t rowB_l = MERGE_LOW (rowB[0], rowB[1]);
	      vec_t rowA_h = MERGE_HIGH (rowA[0], vzero);
	      vec_t rowA_l = MERGE_LOW (rowA[0], vzero);
	      vec_t rowA2_h = MERGE_HIGH (rowA[1], vzero);
	      vec_t rowA2_l = MERGE_LOW (rowA[1], vzero);
	      MMA (&acc0, rowB_h, rowA_h);
	      MMA (&acc1, rowB_l, rowA_h);
	      MMA (&acc2, rowB_h, rowA_l);
	      MMA (&acc3, rowB_l, rowA_l);
	      MMA (&acc4, rowB_h, rowA2_h);
	      MMA (&acc5, rowB_l, rowA2_h);
	      MMA (&acc6, rowB_h, rowA2_l);
	      MMA (&acc7, rowB_l, rowA2_l);
	    }
	  SAVE_ACC (&acc0, 0);
	  SAVE_ACC (&acc2, 4);
	  SAVE_ACC1 (&acc1, 0);
	  SAVE_ACC1 (&acc3, 4);
	  SAVE_ACC (&acc4, 8);
	  SAVE_ACC (&acc6, 12);
	  SAVE_ACC1 (&acc5, 8);
	  SAVE_ACC1 (&acc7, 12);
	  CO += 16;

	  AO += (k << 4);
	  BO += (k << 3);
	}
      i = (m & 15) >> 3;
      /* Loop for m >= 8. */
      for (j = 0; j < i; j++)
	{
	  IFLOAT *BO = B;
	  v4sf_t *rowC;
	  v4sf_t result[4];
	  __vector_quad acc0, acc1, acc2, acc3;
	  SET_ACC_ZERO4 ();
	  BLASLONG l = 0;
	  for (l = 0; l < k / 2; l++)
	    {
	      vec_t *rowA = (vec_t *) & (AO[l << 4]);
	      vec_t *rowB = (vec_t *) & (BO[l << 4]);
	      vec_t rowB_h = MERGE_HIGH (rowB[0], rowB[1]);
	      vec_t rowB_l = MERGE_LOW (rowB[0], rowB[1]);
	      vec_t rowA_h = MERGE_HIGH (rowA[0], rowA[1]);
	      vec_t rowA_l = MERGE_LOW (rowA[0], rowA[1]);
	      MMA (&acc0, rowB_h, rowA_h);
	      MMA (&acc1, rowB_l, rowA_h);
	      MMA (&acc2, rowB_h, rowA_l);
	      MMA (&acc3, rowB_l, rowA_l);
	    }
	  if (k % 2 == 1)
	    {
	      if (k > 1)
		l = (k / 2) << 4;
	      vec_t *rowA = (vec_t *) & (AO[l]);
	      vec_t *rowB = (vec_t *) & (BO[l]);
	      vec_t rowB_h = MERGE_HIGH (rowB[0], rowB[1]);
	      vec_t rowB_l = MERGE_LOW (rowB[0], rowB[1]);
	      vec_t rowA_h = MERGE_HIGH (rowA[0], vzero);
	      vec_t rowA_l = MERGE_LOW (rowA[0], vzero);
	      MMA (&acc0, rowB_h, rowA_h);
	      MMA (&acc1, rowB_l, rowA_h);
	      MMA (&acc2, rowB_h, rowA_l);
	      MMA (&acc3, rowB_l, rowA_l);
	    }
	  SAVE_ACC (&acc0, 0);
	  SAVE_ACC (&acc2, 4);
	  SAVE_ACC1 (&acc1, 0);
	  SAVE_ACC1 (&acc3, 4);
	  CO += 8;
	  AO += (k << 3);
	  BO += (k << 3);
	}
      i = (m & 7) >> 2;
      /* Loop for m >= 4. */
      for (j = 0; j < i; j++)
	{
	  IFLOAT *BO = B;
	  v4sf_t *rowC;
	  v4sf_t result[4];
	  __vector_quad acc0, acc1;
	  __builtin_mma_xxsetaccz (&acc0);
	  __builtin_mma_xxsetaccz (&acc1);
	  BLASLONG l = 0;
	  for (l = 0; l < k / 2; l++)
	    {
	      vec_t *rowA = (vec_t *) & (AO[l << 3]);
	      vec_t *rowB = (vec_t *) & (BO[l << 4]);
	      vec_t rowA_mrg = MERGE_ROW (rowA[0]);
	      MMA (&acc0, MERGE_HIGH (rowB[0], rowB[1]), rowA_mrg);
	      MMA (&acc1, MERGE_LOW (rowB[0], rowB[1]), rowA_mrg);
	    }
	  if (k % 2 == 1)
	    {
	      if (k > 1)
		l = (k / 2) << 3;
	      vector short rowA =
		{ AO[l + 0], 0, AO[l + 1], 0, AO[l + 2], 0, AO[l + 3], 0 };
	      vec_t *rowB = (vec_t *) & (BO[l << 1]);
	      MMA (&acc0, MERGE_HIGH (rowB[0], rowB[1]), (vec_t) rowA);
	      MMA (&acc1, MERGE_LOW (rowB[0], rowB[1]), (vec_t) rowA);
	    }
	  SAVE_ACC (&acc0, 0);
	  SAVE_ACC1 (&acc1, 0);
	  CO += 4;
	  AO += (k << 2);
	  BO += (k << 3);
	}
      i = (m & 3) >> 1;
      /* Loop for m >= 2. */
      for (j = 0; j < i; j++)
	{
	  IFLOAT *BO = B;
	  v2sf_t *rowC;
	  v2sf_t result[8];
	  __vector_quad acc0, acc1;
	  __builtin_mma_xxsetaccz (&acc0);
	  __builtin_mma_xxsetaccz (&acc1);
	  BLASLONG l = 0;
	  for (l = 0; l < k / 2; l++)
	    {
	      vector short rowA =
		{ AO[(l << 2) + 0], AO[(l << 2) + 2], AO[(l << 2) + 1],
		AO[(l << 2) + 3],
		0, 0, 0, 0
	      };
	      vec_t *rowB = (vec_t *) & (BO[l << 4]);
	      MMA (&acc0, MERGE_HIGH (rowB[0], rowB[1]), (vec_t) rowA);
	      MMA (&acc1, MERGE_LOW (rowB[0], rowB[1]), (vec_t) rowA);
	    }
	  if (k % 2 == 1)
	    {
	      if (k > 1)
		l = (k / 2) << 2;
	      vector short rowA = { AO[l + 0], 0, AO[l + 1], 0, 0, 0, 0, 0 };
	      vec_t *rowB = (vec_t *) & (BO[(l << 2)]);
	      MMA (&acc0, MERGE_HIGH (rowB[0], rowB[1]), (vec_t) rowA);
	      MMA (&acc1, MERGE_LOW (rowB[0], rowB[1]), (vec_t) rowA);
	    }
	  SAVE4x2_ACC (&acc0, 0);
	  SAVE4x2_ACC1 (&acc1, 0);
	  CO += 2;
	  AO += (k << 1);
	  BO += (k << 3);
	}
      i = (m & 1) >> 0;
      /* Loop for m = 1. */
      for (j = 0; j < i; j++)
	{
	  IFLOAT *BO = B;
	  BLASLONG l = 0;
	  v4sf_t t = { 0, 0, 0, 0 }
	  , t1 =
	  {
	  0, 0, 0, 0};
	  for (l = 0; l < k; l++)
	    {
	      v4sf_t rowA =
		{ BF16TOF32 (AO[l]), BF16TOF32 (AO[l]), BF16TOF32 (AO[l]),
		BF16TOF32 (AO[l])
	      };
	      v4sf_t rowB =
		{ BF16TOF32 (BO[l << 3]), BF16TOF32 (BO[(l << 3) + 1]),
		BF16TOF32 (BO[(l << 3) + 2]),
		BF16TOF32 (BO[(l << 3) + 3])
	      };
	      v4sf_t rowB1 =
		{ BF16TOF32 (BO[(l << 3) + 4]), BF16TOF32 (BO[(l << 3) + 5]),
		BF16TOF32 (BO[(l << 3) + 6]),
		BF16TOF32 (BO[(l << 3) + 7])
	      };
	      t += rowA * rowB;
	      t1 += rowA * rowB1;
	    }
	  t = t * valpha;
	  t1 = t1 * valpha;
	  CO[0 * ldc] += t[0];
	  CO[1 * ldc] += t[1];
	  CO[2 * ldc] += t[2];
	  CO[3 * ldc] += t[3];
	  CO[4 * ldc] += t1[0];
	  CO[5 * ldc] += t1[1];
	  CO[6 * ldc] += t1[2];
	  CO[7 * ldc] += t1[3];
	  CO += 1;
	  AO += k;
	  BO += (k << 3);
	}
      B += k << 3;
    }
  N = (n & 7) >> 2;
  /* Loop for n >= 4. */
  for (i1 = 0; i1 < N; i1++)
    {
      BLASLONG i, j;
      FLOAT *CO;
      IFLOAT *AO;
      CO = C;
      C += ldc << 2;
      AO = A;
      i = m >> 5;
      /* Loop for m >= 32. */
      for (j = 0; j < i; j++)
	{
	  IFLOAT *BO = B;
	  IFLOAT *A1 = AO + (16 * k);
	  v4sf_t *rowC;
	  v4sf_t result[4];
	  __vector_quad acc0, acc1, acc2, acc3, acc4, acc5, acc6, acc7;
	  SET_ACC_ZERO8 ();
	  BLASLONG l = 0;
	  for (l = 0; l < k / 2; l++)
	    {
	      vec_t *rowA = (vec_t *) & (AO[l << 5]);
	      vec_t *rowA1 = (vec_t *) & (A1[l << 5]);
	      vec_t *rowB = (vec_t *) & (BO[l << 3]);
	      vec_t rowB_mrg = MERGE_ROW (rowB[0]);
	      MMA (&acc0, rowB_mrg, MERGE_HIGH (rowA[0], rowA[2]));
	      MMA (&acc1, rowB_mrg, MERGE_LOW (rowA[0], rowA[2]));
	      MMA (&acc2, rowB_mrg, MERGE_HIGH (rowA[1], rowA[3]));
	      MMA (&acc3, rowB_mrg, MERGE_LOW (rowA[1], rowA[3]));
	      MMA (&acc4, rowB_mrg, MERGE_HIGH (rowA1[0], rowA1[2]));
	      MMA (&acc5, rowB_mrg, MERGE_LOW (rowA1[0], rowA1[2]));
	      MMA (&acc6, rowB_mrg, MERGE_HIGH (rowA1[1], rowA1[3]));
	      MMA (&acc7, rowB_mrg, MERGE_LOW (rowA1[1], rowA1[3]));
	    }
	  if (k % 2 == 1)
	    {
	      if (k > 1)
		l = (k / 2) << 3;
	      vec_t *rowA = (vec_t *) & (AO[(l << 2)]);
	      vec_t *rowA1 = (vec_t *) & (A1[(l << 2)]);
	      vec_t *rowB = (vec_t *) & (BO[l]);
	      vec_t rowB_mrg = MERGE_ROW (rowB[0]);
	      MMA (&acc0, rowB_mrg, MERGE_HIGH (rowA[0], vzero));
	      MMA (&acc1, rowB_mrg, MERGE_LOW (rowA[0], vzero));
	      MMA (&acc2, rowB_mrg, MERGE_HIGH (rowA[1], vzero));
	      MMA (&acc3, rowB_mrg, MERGE_LOW (rowA[1], vzero));
	      MMA (&acc4, rowB_mrg, MERGE_HIGH (rowA1[0], vzero));
	      MMA (&acc5, rowB_mrg, MERGE_LOW (rowA1[0], vzero));
	      MMA (&acc6, rowB_mrg, MERGE_HIGH (rowA1[1], vzero));
	      MMA (&acc7, rowB_mrg, MERGE_LOW (rowA1[1], vzero));
	    }

	  SAVE_ACC (&acc0, 0);
	  SAVE_ACC (&acc1, 4);
	  CO += 8;
	  SAVE_ACC (&acc2, 0);
	  SAVE_ACC (&acc3, 4);
	  CO += 8;
	  SAVE_ACC (&acc4, 0);
	  SAVE_ACC (&acc5, 4);
	  CO += 8;
	  SAVE_ACC (&acc6, 0);
	  SAVE_ACC (&acc7, 4);
	  CO += 8;
	  AO += k << 5;
	  BO += k << 2;
	}
      i = (m & 31) >> 4;
      /* Loop for m >= 16. */
      for (j = 0; j < i; j++)
	{
	  IFLOAT *BO = B;
	  v4sf_t *rowC;
	  v4sf_t result[4];
	  __vector_quad acc0, acc1, acc2, acc3;
	  SET_ACC_ZERO4 ();
	  BLASLONG l = 0;
	  for (l = 0; l < k / 2; l++)
	    {
	      vec_t *rowA = (vec_t *) & (AO[l << 5]);
	      vec_t *rowB = (vec_t *) & (BO[l << 3]);
	      vec_t rowB_mrg = MERGE_ROW (rowB[0]);
	      MMA (&acc0, rowB_mrg, MERGE_HIGH (rowA[0], rowA[2]));
	      MMA (&acc1, rowB_mrg, MERGE_LOW (rowA[0], rowA[2]));
	      MMA (&acc2, rowB_mrg, MERGE_HIGH (rowA[1], rowA[3]));
	      MMA (&acc3, rowB_mrg, MERGE_LOW (rowA[1], rowA[3]));
	    }
	  if (k % 2 == 1)
	    {
	      if (k > 1)
		l = (k / 2) << 3;
	      vec_t *rowA = (vec_t *) & (AO[(l << 2)]);
	      vec_t *rowB = (vec_t *) & (BO[l]);
	      vec_t rowB_mrg = MERGE_ROW (rowB[0]);
	      MMA (&acc0, rowB_mrg, MERGE_HIGH (rowA[0], vzero));
	      MMA (&acc1, rowB_mrg, MERGE_LOW (rowA[0], vzero));
	      MMA (&acc2, rowB_mrg, MERGE_HIGH (rowA[1], vzero));
	      MMA (&acc3, rowB_mrg, MERGE_LOW (rowA[1], vzero));
	    }

	  SAVE_ACC (&acc0, 0);
	  SAVE_ACC (&acc1, 4);
	  CO += 8;
	  SAVE_ACC (&acc2, 0);
	  SAVE_ACC (&acc3, 4);
	  CO += 8;
	  AO += k << 4;
	  BO += k << 2;
	}
      i = (m & 15) >> 3;
      /* Loop for m >= 8. */
      for (j = 0; j < i; j++)
	{
	  IFLOAT *BO = B;
	  v4sf_t *rowC;
	  v4sf_t result[4];
	  __vector_quad acc0, acc1;
	  __builtin_mma_xxsetaccz (&acc0);
	  __builtin_mma_xxsetaccz (&acc1);
	  BLASLONG l = 0;
	  for (l = 0; l < k / 2; l++)
	    {
	      vec_t *rowA = (vec_t *) & (AO[l << 4]);
	      vec_t *rowB = (vec_t *) & (BO[l << 3]);
	      vec_t rowB_mrg = MERGE_ROW (rowB[0]);
	      MMA (&acc0, rowB_mrg, MERGE_HIGH (rowA[0], rowA[1]));
	      MMA (&acc1, rowB_mrg, MERGE_LOW (rowA[0], rowA[1]));
	    }
	  if (k % 2 == 1)
	    {
	      if (k > 1)
		l = (k / 2) << 3;
	      vec_t *rowA = (vec_t *) & (AO[l << 1]);
	      vec_t *rowB = (vec_t *) & (BO[l]);
	      vec_t rowB_mrg = MERGE_ROW (rowB[0]);
	      MMA (&acc0, rowB_mrg, MERGE_HIGH (rowA[0], vzero));
	      MMA (&acc1, rowB_mrg, MERGE_LOW (rowA[0], vzero));
	    }
	  SAVE_ACC (&acc0, 0);
	  SAVE_ACC (&acc1, 4);
	  CO += 8;
	  AO += k << 3;
	  BO += k << 2;
	}
      i = (m & 7) >> 2;
      /* Loop for m >= 4. */
      for (j = 0; j < i; j++)
	{
	  IFLOAT *BO = B;
	  v4sf_t *rowC;
	  __vector_quad acc0;
	  v4sf_t result[4];
	  BLASLONG l = 0;
	  __builtin_mma_xxsetaccz (&acc0);
	  for (l = 0; l < k / 2; l++)
	    {
	      vec_t *rowA = (vec_t *) & (AO[l << 3]);
	      vec_t *rowB = (vec_t *) & (BO[l << 3]);
	      MMA (&acc0, MERGE_ROW (rowB[0]), MERGE_ROW (rowA[0]));
	    }
	  if (k % 2 == 1)
	    {
	      if (k > 1)
		l = (k / 2) << 3;
	      vector short rowA =
		{ AO[l], 0, AO[l + 1], 0, AO[l + 2], 0, AO[l + 3], 0 };
	      vec_t *rowB = (vec_t *) & (BO[l]);
	      MMA (&acc0, MERGE_ROW (rowB[0]), (vec_t) rowA);
	    }
	  SAVE_ACC (&acc0, 0);
	  CO += 4;
	  AO += k << 2;
	  BO += k << 2;
	}
      i = (m & 3) >> 1;
      /* Loop for m >= 2. */
      for (j = 0; j < i; j++)
	{
	  IFLOAT *BO = B;
	  v2sf_t *rowC;
	  v2sf_t result[8];
	  __vector_quad acc0;
	  BLASLONG l = 0;
	  __builtin_mma_xxsetaccz (&acc0);
	  for (l = 0; l < k / 2; l++)
	    {
	      vector short rowA =
		{ AO[(l << 2) + 0], AO[(l << 2) + 2], AO[(l << 2) + 1],
		AO[(l << 2) + 3],
		0, 0, 0, 0
	      };
	      vec_t *rowB = (vec_t *) & (BO[l << 3]);
	      MMA (&acc0, MERGE_ROW (rowB[0]), (vec_t) rowA);
	    }
	  if (k % 2 == 1)
	    {
	      if (k > 1)
		l = (k / 2) << 2;
	      vector short rowA = { AO[l], 0, AO[l + 1], 0, 0, 0, 0, 0 };
	      vec_t *rowB = (vec_t *) & (BO[l << 1]);
	      MMA (&acc0, MERGE_ROW (rowB[0]), (vec_t) rowA);
	    }
	  SAVE4x2_ACC (&acc0, 0);
	  CO += 2;
	  AO += k << 1;
	  BO += k << 2;
	}
      i = (m & 1) >> 0;
      /* Loop for m = 1. */
      for (j = 0; j < i; j++)
	{
	  IFLOAT *BO = B;
	  BLASLONG l = 0;
	  v4sf_t t = { 0, 0, 0, 0 };
	  for (l = 0; l < k; l++)
	    {
	      v4sf_t rowA =
		{ BF16TOF32 (AO[l]), BF16TOF32 (AO[l]), BF16TOF32 (AO[l]),
		BF16TOF32 (AO[l])
	      };
	      v4sf_t rowB =
		{ BF16TOF32 (BO[l << 2]), BF16TOF32 (BO[(l << 2) + 1]),
		BF16TOF32 (BO[(l << 2) + 2]),
		BF16TOF32 (BO[(l << 2) + 3])
	      };
	      t += rowA * rowB;
	    }
	  t = t * valpha;
	  CO[0 * ldc] += t[0];
	  CO[1 * ldc] += t[1];
	  CO[2 * ldc] += t[2];
	  CO[3 * ldc] += t[3];
	  AO += k;
	  BO += (k << 2);
	  CO += 1;
	}

      B += k << 2;
    }
  N = (n & 3) >> 1;
  /* Loop for n >= 2. */
  for (i1 = 0; i1 < N; i1++)
    {
      BLASLONG i, j;
      FLOAT *CO;
      IFLOAT *AO;
      CO = C;
      C += ldc << 1;
      AO = A;
      i = m >> 5;
      /* Loop for m >= 32. */
      for (j = 0; j < i; j++)
	{
	  IFLOAT *BO = B;
	  v4sf_t *rowC;
	  v4sf_t result[4];
	  IFLOAT *A1 = AO + (16 * k);
	  __vector_quad acc0, acc1, acc2, acc3, acc4, acc5, acc6, acc7;
	  SET_ACC_ZERO8 ();
	  BLASLONG l = 0;
	  for (l = 0; l < k / 2; l++)
	    {
	      vector short rowB =
		{ BO[(l << 2) + 0], BO[(l << 2) + 2], BO[(l << 2) + 1],
		BO[(l << 2) + 3],
		0, 0, 0, 0
	      };
	      vec_t *rowA = (vec_t *) & (AO[l << 5]);
	      vec_t *rowA1 = (vec_t *) & (A1[l << 5]);
	      MMA (&acc0, (vec_t) rowB, MERGE_HIGH (rowA[0], rowA[2]));
	      MMA (&acc1, (vec_t) rowB, MERGE_LOW (rowA[0], rowA[2]));
	      MMA (&acc2, (vec_t) rowB, MERGE_HIGH (rowA[1], rowA[3]));
	      MMA (&acc3, (vec_t) rowB, MERGE_LOW (rowA[1], rowA[3]));
	      MMA (&acc4, (vec_t) rowB, MERGE_HIGH (rowA1[0], rowA1[2]));
	      MMA (&acc5, (vec_t) rowB, MERGE_LOW (rowA1[0], rowA1[2]));
	      MMA (&acc6, (vec_t) rowB, MERGE_HIGH (rowA1[1], rowA1[3]));
	      MMA (&acc7, (vec_t) rowB, MERGE_LOW (rowA1[1], rowA1[3]));
	    }
	  if (k % 2 == 1)
	    {
	      if (k > 1)
		l = (k / 2) << 2;
	      vector short rowB = { BO[l + 0], 0, BO[l + 1], 0, 0, 0, 0, 0 };
	      vec_t *rowA = (vec_t *) & (AO[l << 3]);
	      vec_t *rowA1 = (vec_t *) & (A1[l << 3]);
	      MMA (&acc0, (vec_t) rowB, MERGE_HIGH (rowA[0], rowA[2]));
	      MMA (&acc1, (vec_t) rowB, MERGE_LOW (rowA[0], rowA[2]));
	      MMA (&acc2, (vec_t) rowB, MERGE_HIGH (rowA[1], rowA[3]));
	      MMA (&acc3, (vec_t) rowB, MERGE_LOW (rowA[1], rowA[3]));
	      MMA (&acc4, (vec_t) rowB, MERGE_HIGH (rowA1[0], rowA1[2]));
	      MMA (&acc5, (vec_t) rowB, MERGE_LOW (rowA1[0], rowA1[2]));
	      MMA (&acc6, (vec_t) rowB, MERGE_HIGH (rowA1[1], rowA1[3]));
	      MMA (&acc7, (vec_t) rowB, MERGE_LOW (rowA1[1], rowA1[3]));
	    }
	  SAVE2x4_ACC (&acc0, 0);
	  SAVE2x4_ACC (&acc1, 4);
	  SAVE2x4_ACC (&acc2, 8);
	  SAVE2x4_ACC (&acc3, 12);
	  CO += 16;
	  SAVE2x4_ACC (&acc4, 0);
	  SAVE2x4_ACC (&acc5, 4);
	  SAVE2x4_ACC (&acc6, 8);
	  SAVE2x4_ACC (&acc7, 12);
	  CO += 16;
	  AO += k << 5;
	  BO += k << 1;
	}
      i = (m & 31) >> 4;
      /* Loop for m >= 16. */
      for (j = 0; j < i; j++)
	{
	  IFLOAT *BO = B;
	  v4sf_t *rowC;
	  v4sf_t result[4];
	  __vector_quad acc0, acc1, acc2, acc3;
	  SET_ACC_ZERO4 ();
	  BLASLONG l = 0;
	  for (l = 0; l < k / 2; l++)
	    {
	      vector short rowB =
		{ BO[(l << 2) + 0], BO[(l << 2) + 2], BO[(l << 2) + 1],
		BO[(l << 2) + 3],
		0, 0, 0, 0
	      };
	      vec_t *rowA = (vec_t *) & (AO[l << 5]);
	      MMA (&acc0, (vec_t) rowB, MERGE_HIGH (rowA[0], rowA[2]));
	      MMA (&acc1, (vec_t) rowB, MERGE_LOW (rowA[0], rowA[2]));
	      MMA (&acc2, (vec_t) rowB, MERGE_HIGH (rowA[1], rowA[3]));
	      MMA (&acc3, (vec_t) rowB, MERGE_LOW (rowA[1], rowA[3]));
	    }
	  if (k % 2 == 1)
	    {
	      if (k > 1)
		l = (k / 2) << 2;
	      vector short rowB = { BO[l + 0], 0, BO[l + 1], 0, 0, 0, 0, 0 };
	      vec_t *rowA = (vec_t *) & (AO[l << 3]);
	      MMA (&acc0, (vec_t) rowB, MERGE_HIGH (rowA[0], rowA[2]));
	      MMA (&acc1, (vec_t) rowB, MERGE_LOW (rowA[0], rowA[2]));
	      MMA (&acc2, (vec_t) rowB, MERGE_HIGH (rowA[1], rowA[3]));
	      MMA (&acc3, (vec_t) rowB, MERGE_LOW (rowA[1], rowA[3]));
	    }
	  SAVE2x4_ACC (&acc0, 0);
	  SAVE2x4_ACC (&acc1, 4);
	  SAVE2x4_ACC (&acc2, 8);
	  SAVE2x4_ACC (&acc3, 12);
	  CO += 16;
	  AO += k << 4;
	  BO += k << 1;
	}
      i = (m & 15) >> 3;
      /* Loop for m >= 8. */
      for (j = 0; j < i; j++)
	{
	  IFLOAT *BO = B;
	  v4sf_t *rowC;
	  v4sf_t result[4];
	  __vector_quad acc0, acc1;
	  __builtin_mma_xxsetaccz (&acc0);
	  __builtin_mma_xxsetaccz (&acc1);
	  BLASLONG l = 0;
	  for (l = 0; l < k / 2; l++)
	    {
	      vector short rowB =
		{ BO[(l << 2) + 0], BO[(l << 2) + 2], BO[(l << 2) + 1],
		BO[(l << 2) + 3],
		0, 0, 0, 0
	      };
	      vec_t *rowA = (vec_t *) & (AO[l << 4]);
	      MMA (&acc0, (vec_t) rowB, MERGE_HIGH (rowA[0], rowA[1]));
	      MMA (&acc1, (vec_t) rowB, MERGE_LOW (rowA[0], rowA[1]));
	    }
	  if (k % 2 == 1)
	    {
	      if (k > 1)
		l = (k / 2) << 2;
	      vector short rowB = { BO[l + 0], 0, BO[l + 1], 0, 0, 0, 0, 0 };
	      vec_t *rowA = (vec_t *) & (AO[(l << 2)]);
	      MMA (&acc0, (vec_t) rowB, MERGE_HIGH (rowA[0], rowA[1]));
	      MMA (&acc1, (vec_t) rowB, MERGE_LOW (rowA[0], rowA[1]));
	    }
	  SAVE2x4_ACC (&acc0, 0);
	  SAVE2x4_ACC (&acc1, 4);
	  CO += 8;
	  AO += k << 3;
	  BO += k << 1;
	}
      i = (m & 7) >> 2;
      /* Loop for m >= 4. */
      for (j = 0; j < i; j++)
	{
	  IFLOAT *BO = B;
	  v4sf_t *rowC;
	  v4sf_t result[4];
	  __vector_quad acc0;
	  __builtin_mma_xxsetaccz (&acc0);
	  BLASLONG l = 0;
	  for (l = 0; l < k / 2; l++)
	    {
	      vector short rowB =
		{ BO[(l << 2) + 0], BO[(l << 2) + 2], BO[(l << 2) + 1],
		BO[(l << 2) + 3],
		0, 0, 0, 0
	      };
	      vec_t *rowA = (vec_t *) & (AO[l << 3]);
	      MMA (&acc0, (vec_t) rowB, MERGE_ROW (rowA[0]));
	    }
	  if (k % 2 == 1)
	    {
	      if (k > 1)
		l = (k / 2) << 2;
	      vector short rowB = { BO[l + 0], 0, BO[l + 1], 0, 0, 0, 0, 0 };
	      vec_t *rowA = (vec_t *) & (AO[l << 1]);
	      MMA (&acc0, (vec_t) rowB, MERGE_ROW (rowA[0]));
	    }
	  SAVE2x4_ACC (&acc0, 0);
	  CO += 4;
	  AO += k << 2;
	  BO += k << 1;
	}
      i = (m & 3) >> 1;
      /* Loop for m >= 2. */
      for (j = 0; j < i; j++)
	{
	  IFLOAT *BO = B;
	  BLASLONG l = 0;
	  v4sf_t t = { 0, 0, 0, 0 };
	  for (l = 0; l < (k << 1); l += 2)
	    {
	      v4sf_t rowA =
		{ BF16TOF32 (AO[l]), BF16TOF32 (AO[l]), BF16TOF32 (AO[l + 1]),
		BF16TOF32 (AO[l + 1])
	      };
	      v4sf_t rowB =
		{ BF16TOF32 (BO[l]), BF16TOF32 (BO[l + 1]), BF16TOF32 (BO[l]),
		BF16TOF32 (BO[l + 1])
	      };
	      t += rowA * rowB;
	    }
	  t = t * valpha;
	  CO[0 * ldc] += t[0];
	  CO[1 * ldc] += t[1];
	  CO[0 * ldc + 1] += t[2];
	  CO[1 * ldc + 1] += t[3];
	  CO += 2;
	  AO += k << 1;
	  BO += k << 1;
	}
      i = (m & 1) >> 0;
      /* Loop for m = 1. */
      for (j = 0; j < i; j++)
	{
	  IFLOAT *BO = B;
	  BLASLONG l = 0;
	  v4sf_t t = { 0, 0, 0, 0 };
	  for (l = 0; l < k; l++)
	    {
	      v4sf_t rowA = { BF16TOF32 (AO[l]), BF16TOF32 (AO[l]), 0, 0 };
	      v4sf_t rowB =
		{ BF16TOF32 (BO[l << 1]), BF16TOF32 (BO[(l << 1) + 1]), 0,
		0
	      };
	      t += rowA * rowB;
	    }
	  CO[0 * ldc] += t[0] * alpha;
	  CO[1 * ldc] += t[1] * alpha;
	  CO += 1;
	  AO += k;
	  BO += k << 1;
	}
      B += k << 1;
    }
  N = (n & 1) >> 0;
  /* Loop for n = 1. */
  for (i1 = 0; i1 < N; i1++)
    {
      BLASLONG i;
      FLOAT *CO;
      IFLOAT *AO;
      CO = C;
      C += ldc;
      AO = A;
      i = m;
      /* Loop for m >= 16. */
      while (i >= 16)
	{
	  IFLOAT *BO = B;
	  BLASLONG l = 0;
	  v4sf_t t = { 0, 0, 0, 0 };
	  v4sf_t t1 = { 0, 0, 0, 0 };
	  v4sf_t t2 = { 0, 0, 0, 0 };
	  v4sf_t t3 = { 0, 0, 0, 0 };
	  for (l = 0; l < k; l++)
	    {
	      v4sf_t rowB =
		{ BF16TOF32 (BO[l]), BF16TOF32 (BO[l]), BF16TOF32 (BO[l]),
		BF16TOF32 (BO[l])
	      };
	      v4sf_t rowA =
		{ BF16TOF32 (AO[l << 4]), BF16TOF32 (AO[(l << 4) + 1]),
		BF16TOF32 (AO[(l << 4) + 2]),
		BF16TOF32 (AO[(l << 4) + 3])
	      };
	      v4sf_t rowA1 =
		{ BF16TOF32 (AO[(l << 4) + 4]), BF16TOF32 (AO[(l << 4) + 5]),
		BF16TOF32 (AO[(l << 4) + 6]),
		BF16TOF32 (AO[(l << 4) + 7])
	      };
	      v4sf_t rowA2 =
		{ BF16TOF32 (AO[(l << 4) + 8]), BF16TOF32 (AO[(l << 4) + 9]),
		BF16TOF32 (AO[(l << 4) + 10]),
		BF16TOF32 (AO[(l << 4) + 11])
	      };
	      v4sf_t rowA3 = { BF16TOF32 (AO[(l << 4) + 12]),
		BF16TOF32 (AO[(l << 4) + 13]), BF16TOF32 (AO[(l << 4) + 14]),
		BF16TOF32 (AO[(l << 4) + 15])
	      };
	      t += rowA * rowB;
	      t1 += rowA1 * rowB;
	      t2 += rowA2 * rowB;
	      t3 += rowA3 * rowB;
	    }
	  t = t * valpha;
	  t1 = t1 * valpha;
	  t2 = t2 * valpha;
	  t3 = t3 * valpha;
	  CO[0] += t[0];
	  CO[1] += t[1];
	  CO[2] += t[2];
	  CO[3] += t[3];
	  CO[4] += t1[0];
	  CO[5] += t1[1];
	  CO[6] += t1[2];
	  CO[7] += t1[3];
	  CO[8] += t2[0];
	  CO[9] += t2[1];
	  CO[10] += t2[2];
	  CO[11] += t2[3];
	  CO[12] += t3[0];
	  CO[13] += t3[1];
	  CO[14] += t3[2];
	  CO[15] += t3[3];
	  AO += k << 4;
	  BO += k;
	  CO += 16;
	  i -= 16;
	}
      /* Loop for m >= 8. */
      while (i >= 8)
	{
	  IFLOAT *BO = B;
	  BLASLONG l = 0;
	  v4sf_t t = { 0, 0, 0, 0 };
	  v4sf_t t1 = { 0, 0, 0, 0 };
	  for (l = 0; l < k; l++)
	    {
	      v4sf_t rowB =
		{ BF16TOF32 (BO[l]), BF16TOF32 (BO[l]), BF16TOF32 (BO[l]),
		BF16TOF32 (BO[l])
	      };
	      v4sf_t rowA =
		{ BF16TOF32 (AO[l << 3]), BF16TOF32 (AO[(l << 3) + 1]),
		BF16TOF32 (AO[(l << 3) + 2]),
		BF16TOF32 (AO[(l << 3) + 3])
	      };
	      v4sf_t rowA1 =
		{ BF16TOF32 (AO[(l << 3) + 4]), BF16TOF32 (AO[(l << 3) + 5]),
		BF16TOF32 (AO[(l << 3) + 6]),
		BF16TOF32 (AO[(l << 3) + 7])
	      };
	      t += rowA * rowB;
	      t1 += rowA1 * rowB;
	    }
	  t = t * valpha;
	  t1 = t1 * valpha;
	  CO[0] += t[0];
	  CO[1] += t[1];
	  CO[2] += t[2];
	  CO[3] += t[3];
	  CO[4] += t1[0];
	  CO[5] += t1[1];
	  CO[6] += t1[2];
	  CO[7] += t1[3];
	  AO += k << 3;
	  BO += k;
	  CO += 8;
	  i -= 8;
	}
      /* Loop for m >= 4. */
      while (i >= 4)
	{
	  IFLOAT *BO = B;
	  BLASLONG l = 0;
	  v4sf_t t = { 0, 0, 0, 0 };
	  for (l = 0; l < k; l++)
	    {
	      v4sf_t rowB =
		{ BF16TOF32 (BO[l]), BF16TOF32 (BO[l]), BF16TOF32 (BO[l]),
		BF16TOF32 (BO[l])
	      };
	      v4sf_t rowA =
		{ BF16TOF32 (AO[l << 2]), BF16TOF32 (AO[(l << 2) + 1]),
		BF16TOF32 (AO[(l << 2) + 2]),
		BF16TOF32 (AO[(l << 2) + 3])
	      };
	      t += rowA * rowB;
	    }
	  t = t * valpha;
	  CO[0] += t[0];
	  CO[1] += t[1];
	  CO[2] += t[2];
	  CO[3] += t[3];
	  AO += k << 2;
	  BO += k;
	  CO += 4;
	  i -= 4;
	}
      /* Loop for m >= 2. */
      while (i >= 2)
	{
	  IFLOAT *BO = B;
	  BLASLONG l = 0;
	  v4sf_t t = { 0, 0, 0, 0 };
	  for (l = 0; l < k; l++)
	    {
	      v4sf_t rowB = { BF16TOF32 (BO[l]), BF16TOF32 (BO[l]), 0, 0 };
	      v4sf_t rowA =
		{ BF16TOF32 (AO[l << 1]), BF16TOF32 (AO[(l << 1) + 1]), 0,
		0
	      };
	      t += rowA * rowB;
	    }
	  t = t * valpha;
	  CO[0] += t[0];
	  CO[1] += t[1];
	  AO += k << 1;
	  BO += k;
	  CO += 2;
	  i -= 2;
	}
      /* Loop for m = 1. */
      while (i >= 1)
	{
	  IFLOAT *BO = B;
	  BLASLONG l = 0;
	  FLOAT t = 0;
	  for (l = 0; l < k; l++)
	    {
	      t += BF16TOF32 (AO[l]) * BF16TOF32 (BO[l]);
	    }
	  AO += k;
	  BO += k;
	  CO[0] += t * alpha;
	  CO += 1;
	  i -= 1;
	}

      B += k;
    }

  return 0;
}
