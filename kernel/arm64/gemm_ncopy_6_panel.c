/* True MR=6 no-transpose pack (INCOPY) for OpenBLAS DGEMM 6xN.
 *
 * Used when A is transposed (TN/TT/...): INCOPY(m=Kc, n=Mc).
 * Packs 6 columns at a time into contiguous 6-wide micropanels along m
 * (same layout as loongarch64/gemm_ncopy_6.prefx.c).
 *
 * Stock generic/gemm_ncopy_6.c packs as 4-wide and is wrong for 6x8.
 */
#include "common.h"

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, FLOAT *b)
{
  BLASLONG i, j;
  FLOAT *aoffset, *aoffset1, *aoffset2, *aoffset3, *aoffset4, *aoffset5, *aoffset6;
  FLOAT *boffset;
  FLOAT c1, c2, c3, c4, c5, c6;

  aoffset = a;
  boffset = b;

  j = n / 6;
  if (j > 0) {
    do {
      aoffset1 = aoffset;
      aoffset2 = aoffset1 + lda;
      aoffset3 = aoffset2 + lda;
      aoffset4 = aoffset3 + lda;
      aoffset5 = aoffset4 + lda;
      aoffset6 = aoffset5 + lda;
      aoffset += 6 * lda;

      i = m;
      if (i > 0) {
        do {
          c1 = *(aoffset1);
          c2 = *(aoffset2);
          c3 = *(aoffset3);
          c4 = *(aoffset4);
          c5 = *(aoffset5);
          c6 = *(aoffset6);
          aoffset1++; aoffset2++; aoffset3++;
          aoffset4++; aoffset5++; aoffset6++;
          *(boffset + 0) = c1;
          *(boffset + 1) = c2;
          *(boffset + 2) = c3;
          *(boffset + 3) = c4;
          *(boffset + 4) = c5;
          *(boffset + 5) = c6;
          boffset += 6;
          i--;
        } while (i > 0);
      }
      j--;
    } while (j > 0);
  }

  /* n remainder: 4, 2, 1 column groups (matches ukernel N-edge order for A-side) */
  {
    BLASLONG nr = n - (n / 6) * 6;
    if (nr >= 4) {
      aoffset1 = aoffset;
      aoffset2 = aoffset1 + lda;
      aoffset3 = aoffset2 + lda;
      aoffset4 = aoffset3 + lda;
      aoffset += 4 * lda;
      for (i = 0; i < m; i++) {
        boffset[0] = *aoffset1++;
        boffset[1] = *aoffset2++;
        boffset[2] = *aoffset3++;
        boffset[3] = *aoffset4++;
        boffset += 4;
      }
      nr -= 4;
    }
    if (nr >= 2) {
      aoffset1 = aoffset;
      aoffset2 = aoffset1 + lda;
      aoffset += 2 * lda;
      for (i = 0; i < m; i++) {
        boffset[0] = *aoffset1++;
        boffset[1] = *aoffset2++;
        boffset += 2;
      }
      nr -= 2;
    }
    if (nr >= 1) {
      aoffset1 = aoffset;
      for (i = 0; i < m; i++) {
        *boffset++ = *aoffset1++;
      }
    }
  }
  return 0;
}
