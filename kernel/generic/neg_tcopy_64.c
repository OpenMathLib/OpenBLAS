/*********************************************************************/
/* Copyright 2025 The OpenBLAS Project.                              */
/* Bitsliced neg_tcopy MR=64 — panels 64,32,16,8,4,2,1.              */
/*********************************************************************/
#include "common.h"

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, FLOAT *b){
  FLOAT *boffset = b;
  BLASLONG widths[] = {64, 32, 16, 8, 4, 2, 1};
  BLASLONG col = 0;
  int wi;

  for (wi = 0; wi < 7; wi++) {
    BLASLONG w = widths[wi];
    while (n >= w) {
      BLASLONG k, r;
      for (k = 0; k < m; k++) {
        for (r = 0; r < w; r++)
          *boffset++ = -a[(col + r) + k * lda];
      }
      col += w;
      n -= w;
    }
  }
  return 0;
}
