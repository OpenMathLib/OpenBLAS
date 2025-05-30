/***************************************************************************
Copyright (c) 2025 The OpenBLAS Project
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
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/
#include "../common.h"
#include <stdint.h>
#include <stdio.h>

#include <arm_neon.h>

#define SGEMM BLASFUNC(sgemm)
#define BGEMM BLASFUNC(bgemm)
#define BGEMM_LARGEST 256

void *malloc_safe(size_t size) {
  if (size == 0)
    return malloc(1);
  else
    return malloc(size);
}

bfloat16 convert_to_bf16(float x) {
  bfloat16_t src = x;
  bfloat16 dst = 0;
  memcpy(&dst, &src, sizeof(src));
  return dst;
}

int main(int argc, char *argv[]) {
  blasint m, n, k;
  int i, j, l;
  blasint x, y;
  int ret = 0;
  int loop = BGEMM_LARGEST;
  char transA = 'N', transB = 'N';

  float alpha = 1.0, beta = 0.0;

  for (x = 1; x <= loop; x++) {
    if ((x > 100) && (x != BGEMM_LARGEST))
      continue;
    m = k = n = x;

    float *A = (float *)malloc_safe(m * k * sizeof(FLOAT));
    float *B = (float *)malloc_safe(k * n * sizeof(FLOAT));
    float *C = (float *)malloc_safe(m * n * sizeof(FLOAT));

    bfloat16_t *AA = (bfloat16_t *)malloc_safe(m * k * sizeof(bfloat16));
    bfloat16_t *BB = (bfloat16_t *)malloc_safe(k * n * sizeof(bfloat16));
    bfloat16_t *CC = (bfloat16_t *)malloc_safe(m * n * sizeof(bfloat16));

    if ((A == NULL) || (B == NULL) || (C == NULL) || (AA == NULL) ||
        (BB == NULL) || (CC == NULL))
      return 1;

    for (int i = 0; i < m; i++) {
      for (int j = 0; j < k; j++) {
        A[i * k + j] = ((FLOAT) rand () / (FLOAT) RAND_MAX) + 0.5;
        AA[i * k + j] =  A[i * k + j] ;
      }
    }

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < k; j++) {
        // BB[i * k + j] = (i * k + j + 1) % 100;
        B[i * k + j] = ((FLOAT) rand () / (FLOAT) RAND_MAX) + 0.5;
        BB[i * k + j] =  B[i * k + j] ;
      }
    }

    for (y = 0; y < 1; y++) {
      if ((y == 0) || (y == 2)) {
        transA = 'N';
      } else {
        transA = 'T';
      }
      if ((y == 0) || (y == 1)) {
        transB = 'N';
      } else {
        transB = 'T';
      }
      // printf("******** x = %d, y = %d********\n", x, y);
      // printf("Matrix AA (m x k):\n");
      // for (int i = 0; i < m; i++) {
      //   for (int j = 0; j < k; j++) {
      //     printf("%.2f ", (float)AA[i * k + j]);  // or %4.1f if float
      //   }
      //   printf("\n");
      // }

      // printf("Matrix A (copy of AA):\n");
      // for (int i = 0; i < m; i++) {
      //   for (int j = 0; j < k; j++) {
      //     printf("%.2f ", A[i * k + j]);
      //   }
      //   printf("\n");
      // }

      // printf("Matrix BB (n x k):\n");
      // for (int i = 0; i < n; i++) {
      //   for (int j = 0; j < k; j++) {
      //     printf("%.2f ", (float)BB[i * k + j]);
      //   }
      //   printf("\n");
      // }

      // printf("Matrix B (copy of BB):\n");
      // for (int i = 0; i < n; i++) {
      //   for (int j = 0; j < k; j++) {
      //     printf("%.2f ", B[i * k + j]);
      //   }
      //   printf("\n");
      // }

      memset(C, 0, m * n * sizeof(FLOAT));
      memset(CC, 0, m * n * sizeof(bfloat16));
      SGEMM(&transA, &transB, &m, &n, &k, &alpha, A, &m, B, &k, &beta,
            C, &m);
      BGEMM(&transA, &transB, &m, &n, &k, &alpha, (bfloat16 *)AA, &m,
            (bfloat16 *)BB, &k, &beta, (bfloat16 *)CC, &m);
      

      // printf("Matrix CC (n x m):\n");
      // for (int i = 0; i < n; i++) {
      //   for (int j = 0; j < m; j++) {
      //     printf("%.2f ", (float)CC[i * m + j]);
      //   }
      //   printf("\n");
      // }

      // printf("Matrix C :\n");
      // for (int i = 0; i < n; i++) {
      //   for (int j = 0; j < k; j++) {
      //     printf("%.2f ", C[i * k + j]);
      //   }
      //   printf("\n");
      // }

      for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            if (fabs((float)CC[i * m + j] - C[i * m + j]) > 1.0) {
              ret ++;
            }
        }
      }

      printf("x = %d, err = %d\n", x, ret);
      ret = 0;
    }

    free(A);
    free(B);
    free(C);
    free(AA);
    free(BB);
    free(CC);
  }

  if (ret != 0) {
    fprintf(stderr, "FATAL ERROR BGEMM - Return code: %d\n", ret);
    return ret;
  }

  return 0;
}