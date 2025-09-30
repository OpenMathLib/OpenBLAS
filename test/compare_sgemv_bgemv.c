/***************************************************************************
Copyright (c) 2020,2025 The OpenBLAS Project
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

#include "test_helpers.h"

#define SGEMV BLASFUNC(sgemv)
#define BGEMV BLASFUNC(bgemv)
#define BGEMV_LARGEST 256

int main(int argc, char *argv[])
{
  blasint k;
  int i, j, l;
  blasint x, y;
  blasint one = 1;
  int ret = 0;
  int loop = BGEMV_LARGEST;
  char transA = 'N';
  float alpha = 1.0, beta = 0.0;
  bfloat16 alpha_bf16, beta_bf16;

  for (beta = 0; beta < 3; beta += 1)
  {
    for (alpha = 0; alpha < 3; alpha += 1)
    {
      for (l = 0; l < 2; l++)
      { // l = 1 to test inc_x & inc_y not equal to one.
        for (x = 1; x <= loop; x++)
        {
          k = (x == 0) ? 0 : l + 1;
          float *A = (float *)malloc_safe(x * x * sizeof(FLOAT));
          float *B = (float *)malloc_safe(x * sizeof(FLOAT) << l);
          float *C = (float *)malloc_safe(x * sizeof(FLOAT) << l);
          bfloat16 *AA = (bfloat16 *)malloc_safe(x * x * sizeof(bfloat16));
          bfloat16 *BB = (bfloat16 *)malloc_safe(x * sizeof(bfloat16) << l);
          bfloat16 *CC = (bfloat16 *)malloc_safe(x * sizeof(bfloat16) << l);
          float *DD = (float *)malloc_safe(x * sizeof(FLOAT));
          if ((A == NULL) || (B == NULL) || (C == NULL) || (AA == NULL) || (BB == NULL) ||
              (CC == NULL) || (DD == NULL))
            return 1;

          for (j = 0; j < x; j++)
          {
            for (i = 0; i < x; i++)
            {
              A[j * x + i] = ((FLOAT)rand() / (FLOAT)RAND_MAX) + 0.5;
              sbstobf16_(&one, &A[j * x + i], &one, &AA[j * x + i], &one);
            }
            B[j << l] = ((FLOAT)rand() / (FLOAT)RAND_MAX) + 0.5;
            sbstobf16_(&one, &B[j << l], &one, &BB[j << l], &one);

            C[j << l] = ((FLOAT)rand() / (FLOAT)RAND_MAX) + 0.5;
            sbstobf16_(&one, &B[j << l], &one, &CC[j << l], &one);
          }

          for (y = 0; y < 2; y++)
          {
            if (y == 0)
            {
              transA = 'N';
            }
            else
            {
              transA = 'T';
            }

            memset(C, 0, x * sizeof(FLOAT) << l);
            memset(CC, 0, x * sizeof(bfloat16) << l);
            memset(DD, 0, x * sizeof(FLOAT));

            sbstobf16_(&one, &alpha, &one, &alpha_bf16, &one);
            sbstobf16_(&one, &beta, &one, &beta_bf16, &one);
            SGEMV(&transA, &x, &x, &alpha, A, &x, B, &k, &beta, C, &k);
            BGEMV(&transA, &x, &x, &alpha_bf16, AA, &x, BB, &k, &beta_bf16, CC, &k);

            for (i = 0; i < x; i++)
              DD[i] *= beta;

            for (j = 0; j < x; j++)
              for (i = 0; i < x; i++)
                if (transA == 'N')
                {
                  DD[i] += alpha * float16to32(AA[j * x + i]) * float16to32(BB[j << l]);
                }
                else if (transA == 'T')
                {
                  DD[j] += alpha * float16to32(AA[j * x + i]) * float16to32(BB[i << l]);
                }

            for (j = 0; j < x; j++)
            {
              if (!is_close(float16to32(CC[j << l]), truncate_float32_to_bfloat16(C[j << l]), 0.01, 0.001))
              {
#ifdef DEBUG
                printf("Mismatch at trans=%c, alpha=%.2f, beta=%.2f, i=%d, j=%d, k=%d: CC=%.6f, C=%.6f\n",
                       transA, alpha, beta, i, j, k, float16to32(CC[j << l]), truncate_float32_to_bfloat16(C[j << l]));
#endif
                ret++;
              }
              if (!is_close(float16to32(CC[j << l]), truncate_float32_to_bfloat16(DD[j]), 0.001, 0.0001))
              {
#ifdef DEBUG
                printf("Mismatch at trans=%c, alpha=%.2f, beta=%.2f, i=%d, j=%d, k=%d: CC=%.6f, DD=%.6f\n",
                       transA, alpha, beta, i, j, k, float16to32(CC[j << l]), truncate_float32_to_bfloat16(DD[j]));
#endif
                ret++;
              }
            }
          }

          free(A);
          free(B);
          free(C);
          free(AA);
          free(BB);
          free(CC);
          free(DD);
        } // x
      } // l
    } // alpha
  } // beta

  if (ret != 0)
    fprintf(stderr, "FATAL ERROR BGEMV - Return code: %d\n", ret);
  return ret;
}
