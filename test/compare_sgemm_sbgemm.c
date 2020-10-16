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
#include <stdio.h>
#include <stdint.h>
#include "../common.h"
#define SGEMM   BLASFUNC(sgemm)
#define SBGEMM   BLASFUNC(sbgemm)
typedef union
{
  unsigned short v;
  struct
  {
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    unsigned short s:1;
    unsigned short e:8;
    unsigned short m:7;
#else
    unsigned short m:7;
    unsigned short e:8;
    unsigned short s:1;
#endif
  } bits;
} bfloat16_bits;

typedef union
{
  float v;
  struct
  {
    uint32_t m:23;
    uint32_t e:8;
    uint32_t s:1;
  } bits;
} float32_bits;

float
float16to32 (bfloat16_bits f16)
{
  float32_bits f32;
  f32.bits.s = f16.bits.s;
  f32.bits.e = f16.bits.e;
  f32.bits.m = (uint32_t) f16.bits.m << 16;
  return f32.v;
}

int
main (int argc, char *argv[])
{
  int m, n, k;
  int i, j, l;
  int x;
  int ret = 0;
  int loop = 100;
  char transA = 'N', transB = 'N';
  float alpha = 1.0, beta = 0.0;

  for (x = 0; x <= loop; x++)
    {
      m = k = n = x;
      float A[m * k];
      float B[k * n];
      float C[m * n];
      bfloat16_bits AA[m * k], BB[k * n];
      float DD[m * n], CC[m * n];

      for (j = 0; j < m; j++)
	{
	  for (i = 0; i < m; i++)
	    {
	      A[j * k + i] = ((FLOAT) rand () / (FLOAT) RAND_MAX) + 0.5;
	      B[j * k + i] = ((FLOAT) rand () / (FLOAT) RAND_MAX) + 0.5;
	      C[j * k + i] = 0;
	      AA[j * k + i].v = *(uint32_t *) & A[j * k + i] >> 16;
	      BB[j * k + i].v = *(uint32_t *) & B[j * k + i] >> 16;
	      CC[j * k + i] = 0;
	      DD[j * k + i] = 0;
	    }
	}
      SGEMM (&transA, &transB, &m, &n, &k, &alpha, A,
	     &m, B, &k, &beta, C, &m);
      SBGEMM (&transA, &transB, &m, &n, &k, &alpha, AA,
	      &m, BB, &k, &beta, CC, &m);
      for (i = 0; i < n; i++)
	for (j = 0; j < m; j++)
	  for (l = 0; l < k; l++)
	    if (fabs (CC[i * m + j] - C[i * m + j]) > 1.0)
	      ret++;
      if (transA == 'N' && transB == 'N')
	{
	  for (i = 0; i < n; i++)
	    for (j = 0; j < m; j++)
	      for (l = 0; l < k; l++)
		{
		  DD[i * m + j] +=
		    float16to32 (AA[l * m + j]) * float16to32 (BB[l + k * i]);
		}
	  for (i = 0; i < n; i++)
	    for (j = 0; j < m; j++)
	      for (l = 0; l < k; l++)
		if (CC[i * m + j] != DD[i * m + j])
		  ret++;
	}
    }
  if (ret != 0)
    fprintf (stderr, "FATAL ERROR SBGEMM - Return code: %d\n", ret);
  return ret;
}
