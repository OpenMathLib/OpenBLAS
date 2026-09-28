/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* Copyright 2025 The OpenBLAS Project.                              */
/*********************************************************************/

#include "common.h"

#ifndef UNIT
#define INV(a) (ONE / (a))
#else
#define INV(a) (ONE)
#endif

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, BLASLONG offset, FLOAT *b){

  BLASLONG i, ii, j, jj, k;

  FLOAT *a1, *a2, *a3, *a4, *a5, *a6, *a7, *a8, *a9, *a10, *a11, *a12, *a13, *a14, *a15, *a16, *a17, *a18, *a19, *a20, *a21, *a22, *a23, *a24, *a25, *a26, *a27, *a28, *a29, *a30, *a31, *a32;

  jj = offset;

  j = (n >> 5);
  while (j > 0){

    a1  = a +  0 * lda;
    a2  = a +  1 * lda;
    a3  = a +  2 * lda;
    a4  = a +  3 * lda;
    a5  = a +  4 * lda;
    a6  = a +  5 * lda;
    a7  = a +  6 * lda;
    a8  = a +  7 * lda;
    a9  = a +  8 * lda;
    a10  = a +  9 * lda;
    a11  = a +  10 * lda;
    a12  = a +  11 * lda;
    a13  = a +  12 * lda;
    a14  = a +  13 * lda;
    a15  = a +  14 * lda;
    a16  = a +  15 * lda;
    a17  = a +  16 * lda;
    a18  = a +  17 * lda;
    a19  = a +  18 * lda;
    a20  = a +  19 * lda;
    a21  = a +  20 * lda;
    a22  = a +  21 * lda;
    a23  = a +  22 * lda;
    a24  = a +  23 * lda;
    a25  = a +  24 * lda;
    a26  = a +  25 * lda;
    a27  = a +  26 * lda;
    a28  = a +  27 * lda;
    a29  = a +  28 * lda;
    a30  = a +  29 * lda;
    a31  = a +  30 * lda;
    a32  = a +  31 * lda;
    a += 32 * lda;
    ii = 0;

    for (i = 0; i < m; i++) {

      if ((ii >= jj ) && (ii - jj < 32)) {
	*(b +  ii - jj) = INV(*(a1 + (ii - jj) * lda));
	for (k = ii - jj + 1; k < 32; k ++) {
	  *(b +  k) = *(a1 +  k * lda);
	}
      }

      if (ii - jj < 0) {
	*(b +  0) = *(a1  +  0);
	*(b +  1) = *(a2  +  0);
	*(b +  2) = *(a3  +  0);
	*(b +  3) = *(a4  +  0);
	*(b +  4) = *(a5  +  0);
	*(b +  5) = *(a6  +  0);
	*(b +  6) = *(a7  +  0);
	*(b +  7) = *(a8  +  0);
	*(b +  8) = *(a9  +  0);
	*(b +  9) = *(a10  +  0);
	*(b +  10) = *(a11  +  0);
	*(b +  11) = *(a12  +  0);
	*(b +  12) = *(a13  +  0);
	*(b +  13) = *(a14  +  0);
	*(b +  14) = *(a15  +  0);
	*(b +  15) = *(a16  +  0);
	*(b +  16) = *(a17  +  0);
	*(b +  17) = *(a18  +  0);
	*(b +  18) = *(a19  +  0);
	*(b +  19) = *(a20  +  0);
	*(b +  20) = *(a21  +  0);
	*(b +  21) = *(a22  +  0);
	*(b +  22) = *(a23  +  0);
	*(b +  23) = *(a24  +  0);
	*(b +  24) = *(a25  +  0);
	*(b +  25) = *(a26  +  0);
	*(b +  26) = *(a27  +  0);
	*(b +  27) = *(a28  +  0);
	*(b +  28) = *(a29  +  0);
	*(b +  29) = *(a30  +  0);
	*(b +  30) = *(a31  +  0);
	*(b +  31) = *(a32  +  0);
      }

      a1  ++;
      a2  ++;
      a3  ++;
      a4  ++;
      a5  ++;
      a6  ++;
      a7  ++;
      a8  ++;
      a9  ++;
      a10  ++;
      a11  ++;
      a12  ++;
      a13  ++;
      a14  ++;
      a15  ++;
      a16  ++;
      a17  ++;
      a18  ++;
      a19  ++;
      a20  ++;
      a21  ++;
      a22  ++;
      a23  ++;
      a24  ++;
      a25  ++;
      a26  ++;
      a27  ++;
      a28  ++;
      a29  ++;
      a30  ++;
      a31  ++;
      a32  ++;
      b  += 32;
      ii ++;
    }

    jj += 32;
    j --;
  }

  if (n & 16) {

    a1  = a +  0 * lda;
    a2  = a +  1 * lda;
    a3  = a +  2 * lda;
    a4  = a +  3 * lda;
    a5  = a +  4 * lda;
    a6  = a +  5 * lda;
    a7  = a +  6 * lda;
    a8  = a +  7 * lda;
    a9  = a +  8 * lda;
    a10  = a +  9 * lda;
    a11  = a +  10 * lda;
    a12  = a +  11 * lda;
    a13  = a +  12 * lda;
    a14  = a +  13 * lda;
    a15  = a +  14 * lda;
    a16  = a +  15 * lda;
    a += 16 * lda;
    ii = 0;

    for (i = 0; i < m; i++) {

      if ((ii >= jj ) && (ii - jj < 16)) {
	*(b +  ii - jj) = INV(*(a1 + (ii - jj) * lda));
	for (k = ii - jj + 1; k < 16; k ++) {
	  *(b +  k) = *(a1 +  k * lda);
	}
      }

      if (ii - jj < 0) {
	*(b +  0) = *(a1  +  0);
	*(b +  1) = *(a2  +  0);
	*(b +  2) = *(a3  +  0);
	*(b +  3) = *(a4  +  0);
	*(b +  4) = *(a5  +  0);
	*(b +  5) = *(a6  +  0);
	*(b +  6) = *(a7  +  0);
	*(b +  7) = *(a8  +  0);
	*(b +  8) = *(a9  +  0);
	*(b +  9) = *(a10  +  0);
	*(b +  10) = *(a11  +  0);
	*(b +  11) = *(a12  +  0);
	*(b +  12) = *(a13  +  0);
	*(b +  13) = *(a14  +  0);
	*(b +  14) = *(a15  +  0);
	*(b +  15) = *(a16  +  0);
      }

      a1  ++;
      a2  ++;
      a3  ++;
      a4  ++;
      a5  ++;
      a6  ++;
      a7  ++;
      a8  ++;
      a9  ++;
      a10  ++;
      a11  ++;
      a12  ++;
      a13  ++;
      a14  ++;
      a15  ++;
      a16  ++;
      b  += 16;
      ii ++;
    }

    jj += 16;
  }

  if (n & 8) {

    a1  = a +  0 * lda;
    a2  = a +  1 * lda;
    a3  = a +  2 * lda;
    a4  = a +  3 * lda;
    a5  = a +  4 * lda;
    a6  = a +  5 * lda;
    a7  = a +  6 * lda;
    a8  = a +  7 * lda;
    a += 8 * lda;
    ii = 0;

    for (i = 0; i < m; i++) {

      if ((ii >= jj ) && (ii - jj < 8)) {
	*(b +  ii - jj) = INV(*(a1 + (ii - jj) * lda));
	for (k = ii - jj + 1; k < 8; k ++) {
	  *(b +  k) = *(a1 +  k * lda);
	}
      }

      if (ii - jj < 0) {
	*(b +  0) = *(a1  +  0);
	*(b +  1) = *(a2  +  0);
	*(b +  2) = *(a3  +  0);
	*(b +  3) = *(a4  +  0);
	*(b +  4) = *(a5  +  0);
	*(b +  5) = *(a6  +  0);
	*(b +  6) = *(a7  +  0);
	*(b +  7) = *(a8  +  0);
      }

      a1  ++;
      a2  ++;
      a3  ++;
      a4  ++;
      a5  ++;
      a6  ++;
      a7  ++;
      a8  ++;
      b  += 8;
      ii ++;
    }

    jj += 8;
  }

  if (n & 4) {

    a1  = a +  0 * lda;
    a2  = a +  1 * lda;
    a3  = a +  2 * lda;
    a4  = a +  3 * lda;
    a += 4 * lda;
    ii = 0;

    for (i = 0; i < m; i++) {

      if ((ii >= jj ) && (ii - jj < 4)) {
	*(b +  ii - jj) = INV(*(a1 + (ii - jj) * lda));
	for (k = ii - jj + 1; k < 4; k ++) {
	  *(b +  k) = *(a1 +  k * lda);
	}
      }

      if (ii - jj < 0) {
	*(b +  0) = *(a1  +  0);
	*(b +  1) = *(a2  +  0);
	*(b +  2) = *(a3  +  0);
	*(b +  3) = *(a4  +  0);
      }

      a1  ++;
      a2  ++;
      a3  ++;
      a4  ++;
      b  += 4;
      ii ++;
    }

    jj += 4;
  }

  if (n & 2) {

    a1  = a +  0 * lda;
    a2  = a +  1 * lda;
    a += 2 * lda;
    ii = 0;

    for (i = 0; i < m; i++) {

      if ((ii >= jj ) && (ii - jj < 2)) {
	*(b +  ii - jj) = INV(*(a1 + (ii - jj) * lda));
	for (k = ii - jj + 1; k < 2; k ++) {
	  *(b +  k) = *(a1 +  k * lda);
	}
      }

      if (ii - jj < 0) {
	*(b +  0) = *(a1  +  0);
	*(b +  1) = *(a2  +  0);
      }

      a1  ++;
      a2  ++;
      b  += 2;
      ii ++;
    }

    jj += 2;
  }

  if (n & 1) {

    a1  = a +  0 * lda;
    a += 1 * lda;
    ii = 0;

    for (i = 0; i < m; i++) {

      if ((ii >= jj ) && (ii - jj < 1)) {
	*(b +  ii - jj) = INV(*(a1 + (ii - jj) * lda));
	for (k = ii - jj + 1; k < 1; k ++) {
	  *(b +  k) = *(a1 +  k * lda);
	}
      }

      if (ii - jj < 0) {
	*(b +  0) = *(a1  +  0);
      }

      a1  ++;
      b  += 1;
      ii ++;
    }

    jj += 1;
  }

  return 0;
}
