/* True MR=6 transpose pack (ITCOPY) for OpenBLAS DGEMM 6xN.
 *
 * Call signature matches gemm_tcopy_*: (m, n, a, lda, b) where for NN
 * ICOPY uses ITCOPY(m=Kc, n=Mc). Output is Goto-style micropanels:
 *   for each panel of 6 rows: Kc contiguous packs of 6 doubles.
 * Remainders after full panels: 4, then 2, then 1 (kernel cascade).
 *
 * Stock generic/gemm_tcopy_6.c is a clone of tcopy_4 (4+2) and does NOT
 * match a contiguous 6-wide ukernel — do not use it for 6x8.
 */
#include "common.h"

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, FLOAT *b)
{
  BLASLONG j, js;
  FLOAT *b_ptr = b;
  BLASLONG n6 = n / 6;
  BLASLONG nr = n - n6 * 6;

  for (js = 0; js < n6; js++) {
    FLOAT *a_row = a + js * 6;
    for (j = 0; j < m; j++) {
      FLOAT *ap = a_row + j * lda;
      b_ptr[0] = ap[0];
      b_ptr[1] = ap[1];
      b_ptr[2] = ap[2];
      b_ptr[3] = ap[3];
      b_ptr[4] = ap[4];
      b_ptr[5] = ap[5];
      b_ptr += 6;
    }
  }

  FLOAT *a_row = a + n6 * 6;
  if (nr >= 4) {
    for (j = 0; j < m; j++) {
      FLOAT *ap = a_row + j * lda;
      b_ptr[0] = ap[0];
      b_ptr[1] = ap[1];
      b_ptr[2] = ap[2];
      b_ptr[3] = ap[3];
      b_ptr += 4;
    }
    a_row += 4;
    nr -= 4;
  }
  if (nr >= 2) {
    for (j = 0; j < m; j++) {
      FLOAT *ap = a_row + j * lda;
      b_ptr[0] = ap[0];
      b_ptr[1] = ap[1];
      b_ptr += 2;
    }
    a_row += 2;
    nr -= 2;
  }
  if (nr >= 1) {
    for (j = 0; j < m; j++) {
      *b_ptr++ = a_row[j * lda];
    }
  }
  return 0;
}
