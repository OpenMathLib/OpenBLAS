#include "common.h"

/* Correctness-first TRMM microkernel for UNROLL_M=6, UNROLL_N=8.
 * Matches contiguous MR=6 / NR=8 Goto packs (TARGET=CORTEXA72, #5997).
 * Not performance-tuned; replaces the incorrect trmmkernel_2x2 fallback.
 */

int CNAME(BLASLONG bm, BLASLONG bn, BLASLONG bk, FLOAT alpha, FLOAT *ba, FLOAT *bb, FLOAT *C, BLASLONG ldc
#ifdef TRMMKERNEL
	  , BLASLONG offset
#endif
	  )
{
  BLASLONG i, j, k, ii, jj;
  FLOAT *Cptr, *ptrba, *ptrbb;
  FLOAT res[6][8];
  FLOAT a;
  BLASLONG off, temp;
  BLASLONG mr, nr, mrem, nrem;

#if defined(TRMMKERNEL) && !defined(LEFT)
  off = -offset;
#else
  off = 0;
#endif

  for (j = 0; j < bn / 8; j++) {
    Cptr = C;
#if defined(TRMMKERNEL) && defined(LEFT)
    off = offset;
#endif
    ptrba = ba;

    for (i = 0; i < bm / 6; i++) {
#if (defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
      ptrbb = bb;
#else
      ptrba += off * 6;
      ptrbb = bb + off * 8;
#endif
      for (ii = 0; ii < 6; ii++)
	for (jj = 0; jj < 8; jj++)
	  res[ii][jj] = ZERO;

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
      temp = bk - off;
#elif defined(LEFT)
      temp = off + 6;
#else
      temp = off + 8;
#endif
      for (k = 0; k < temp; k++) {
	for (ii = 0; ii < 6; ii++) {
	  a = ptrba[ii];
	  for (jj = 0; jj < 8; jj++)
	    res[ii][jj] += a * ptrbb[jj];
	}
	ptrba += 6;
	ptrbb += 8;
      }
      for (ii = 0; ii < 6; ii++)
	for (jj = 0; jj < 8; jj++)
	  Cptr[ii + jj * ldc] = res[ii][jj] * alpha;

#if (defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
      temp = bk - off;
#ifdef LEFT
      temp -= 6;
#else
      temp -= 8;
#endif
      ptrba += temp * 6;
      ptrbb += temp * 8;
#endif
#ifdef LEFT
      off += 6;
#endif
      Cptr += 6;
    }

    /* M remainder after m/6 panels is 0..5; cascade 4,2,1. */
    mrem = bm % 6;
    for (mr = 4; mr >= 1; mr >>= 1) {
      if (mrem & mr) {
#if (defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	ptrbb = bb;
#else
	ptrba += off * mr;
	ptrbb = bb + off * 8;
#endif
	for (ii = 0; ii < mr; ii++)
	  for (jj = 0; jj < 8; jj++)
	    res[ii][jj] = ZERO;
#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	temp = bk - off;
#elif defined(LEFT)
	temp = off + mr;
#else
	temp = off + 8;
#endif
	for (k = 0; k < temp; k++) {
	  for (ii = 0; ii < mr; ii++) {
	    a = ptrba[ii];
	    for (jj = 0; jj < 8; jj++)
	      res[ii][jj] += a * ptrbb[jj];
	  }
	  ptrba += mr;
	  ptrbb += 8;
	}
	for (ii = 0; ii < mr; ii++)
	  for (jj = 0; jj < 8; jj++)
	    Cptr[ii + jj * ldc] = res[ii][jj] * alpha;
#if (defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	temp = bk - off;
#ifdef LEFT
	temp -= mr;
#else
	temp -= 8;
#endif
	ptrba += temp * mr;
	ptrbb += temp * 8;
#endif
#ifdef LEFT
	off += mr;
#endif
	Cptr += mr;
      }
    }

#if defined(TRMMKERNEL) && !defined(LEFT)
    off += 8;
#endif
    bb += bk * 8;
    C += ldc * 8;
  }

  /* N remainder after n/8 panels is 0..7; cascade 4,2,1. */
  nrem = bn % 8;
  for (nr = 4; nr >= 1; nr >>= 1) {
    if (nrem & nr) {
      Cptr = C;
#if defined(TRMMKERNEL) && defined(LEFT)
      off = offset;
#endif
      ptrba = ba;

      for (i = 0; i < bm / 6; i++) {
#if (defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	ptrbb = bb;
#else
	ptrba += off * 6;
	ptrbb = bb + off * nr;
#endif
	for (ii = 0; ii < 6; ii++)
	  for (jj = 0; jj < nr; jj++)
	    res[ii][jj] = ZERO;
#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	temp = bk - off;
#elif defined(LEFT)
	temp = off + 6;
#else
	temp = off + nr;
#endif
	for (k = 0; k < temp; k++) {
	  for (ii = 0; ii < 6; ii++) {
	    a = ptrba[ii];
	    for (jj = 0; jj < nr; jj++)
	      res[ii][jj] += a * ptrbb[jj];
	  }
	  ptrba += 6;
	  ptrbb += nr;
	}
	for (ii = 0; ii < 6; ii++)
	  for (jj = 0; jj < nr; jj++)
	    Cptr[ii + jj * ldc] = res[ii][jj] * alpha;
#if (defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	temp = bk - off;
#ifdef LEFT
	temp -= 6;
#else
	temp -= nr;
#endif
	ptrba += temp * 6;
	ptrbb += temp * nr;
#endif
#ifdef LEFT
	off += 6;
#endif
	Cptr += 6;
      }

      mrem = bm % 6;
      for (mr = 4; mr >= 1; mr >>= 1) {
	if (mrem & mr) {
#if (defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	  ptrbb = bb;
#else
	  ptrba += off * mr;
	  ptrbb = bb + off * nr;
#endif
	  for (ii = 0; ii < mr; ii++)
	    for (jj = 0; jj < nr; jj++)
	      res[ii][jj] = ZERO;
#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	  temp = bk - off;
#elif defined(LEFT)
	  temp = off + mr;
#else
	  temp = off + nr;
#endif
	  for (k = 0; k < temp; k++) {
	    for (ii = 0; ii < mr; ii++) {
	      a = ptrba[ii];
	      for (jj = 0; jj < nr; jj++)
		res[ii][jj] += a * ptrbb[jj];
	    }
	    ptrba += mr;
	    ptrbb += nr;
	  }
	  for (ii = 0; ii < mr; ii++)
	    for (jj = 0; jj < nr; jj++)
	      Cptr[ii + jj * ldc] = res[ii][jj] * alpha;
#if (defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	  temp = bk - off;
#ifdef LEFT
	  temp -= mr;
#else
	  temp -= nr;
#endif
	  ptrba += temp * mr;
	  ptrbb += temp * nr;
#endif
#ifdef LEFT
	  off += mr;
#endif
	  Cptr += mr;
	}
      }

#if defined(TRMMKERNEL) && !defined(LEFT)
      off += nr;
#endif
      bb += bk * nr;
      C += ldc * nr;
    }
  }
  return 0;
}
