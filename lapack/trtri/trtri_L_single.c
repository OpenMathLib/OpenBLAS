/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/*   1. Redistributions of source code must retain the above         */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer.                                                  */
/*                                                                   */
/*   2. Redistributions in binary form must reproduce the above      */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer in the documentation and/or other materials       */
/*      provided with the distribution.                              */
/*                                                                   */
/*    THIS  SOFTWARE IS PROVIDED  BY THE  UNIVERSITY OF  TEXAS AT    */
/*    AUSTIN  ``AS IS''  AND ANY  EXPRESS OR  IMPLIED WARRANTIES,    */
/*    INCLUDING, BUT  NOT LIMITED  TO, THE IMPLIED  WARRANTIES OF    */
/*    MERCHANTABILITY  AND FITNESS FOR  A PARTICULAR  PURPOSE ARE    */
/*    DISCLAIMED.  IN  NO EVENT SHALL THE UNIVERSITY  OF TEXAS AT    */
/*    AUSTIN OR CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT,    */
/*    INCIDENTAL,  SPECIAL, EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES    */
/*    (INCLUDING, BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE    */
/*    GOODS  OR  SERVICES; LOSS  OF  USE,  DATA,  OR PROFITS;  OR    */
/*    BUSINESS INTERRUPTION) HOWEVER CAUSED  AND ON ANY THEORY OF    */
/*    LIABILITY, WHETHER  IN CONTRACT, STRICT  LIABILITY, OR TORT    */
/*    (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY WAY OUT    */
/*    OF  THE  USE OF  THIS  SOFTWARE,  EVEN  IF ADVISED  OF  THE    */
/*    POSSIBILITY OF SUCH DAMAGE.                                    */
/*                                                                   */
/* The views and conclusions contained in the software and           */
/* documentation are those of the authors and should not be          */
/* interpreted as representing official policies, either expressed   */
/* or implied, of The University of Texas at Austin.                 */
/*********************************************************************/

#include <stdio.h>
#include "common.h"

static FLOAT dp1 =  1.;
static FLOAT dm1 = -1.;

#ifdef UNIT
#define TRTI2	TRTI2_LU
#else
#define TRTI2	TRTI2_LN
#endif

#if 0
#undef GEMM_P
#undef GEMM_Q
#undef GEMM_R

#define GEMM_P 8
#define GEMM_Q 20
#define GEMM_R 64
#endif

#define GEMM_PQ  MAX(GEMM_P, GEMM_Q)
#define REAL_GEMM_R (GEMM_R - 2 * GEMM_PQ)

blasint CNAME(blas_arg_t *args, BLASLONG *range_m, BLASLONG *range_n, FLOAT *sa, FLOAT *sb, BLASLONG myid) {

  BLASLONG  n, lda;
  FLOAT *a;

  BLASLONG i, is, min_i, start_i;
  BLASLONG ls, min_l;
  BLASLONG bk;
  BLASLONG blocking;
  BLASLONG range_N[2];

  FLOAT *sa_trsm = (FLOAT *)((BLASLONG)sb);
  FLOAT *sa_trmm = (FLOAT *)((((BLASLONG)sb
			       + GEMM_PQ * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)
			     + GEMM_OFFSET_A);
  FLOAT *sb_gemm = (FLOAT *)((((BLASLONG)sa_trmm
			       + GEMM_PQ * GEMM_Q * COMPSIZE * SIZE + GEMM_ALIGN) & ~GEMM_ALIGN)
			     + GEMM_OFFSET_B);

  n      = args -> n;
  a      = (FLOAT *)args -> a;
  lda    = args -> lda;

  if (range_n) {
    n      = range_n[1] - range_n[0];
    a     += range_n[0] * (lda + 1) * COMPSIZE;
  }

  if (n <= DTB_ENTRIES) {
    TRTI2(args, NULL, range_n, sa, sb, 0);
    return 0;
  }

  blocking = GEMM_Q;
  if (n <= 4 * GEMM_Q) blocking = (n + 3) / 4;

  start_i = 0;
  while (start_i < n) start_i += blocking;
  start_i -= blocking;

  for (i = start_i; i >= 0; i -= blocking) {
    bk = MIN(blocking, n - i);
    
    if (n - bk - i > 0) TRSM_OLNCOPY(bk, bk, a + (i + i * lda) * COMPSIZE, lda, 0, sa_trsm);
	
    if (!range_n) {
      range_N[0] = i;
      range_N[1] = i + bk;
    } else {
      range_N[0] = range_n[0] + i;
      range_N[1] = range_n[0] + i + bk;
    }

    CNAME(args, NULL, range_N, sa, sa_trmm, 0);

    if (i > 0) {
      TRMM_ILTCOPY(bk, bk, a + (i + i * lda) * COMPSIZE, lda, 0, 0, sa_trmm);

      for (ls = 0; ls < i; ls += REAL_GEMM_R) {
	min_l = i - ls;
	if (min_l > REAL_GEMM_R) min_l = REAL_GEMM_R;
	
	GEMM_ONCOPY (bk, min_l, a + (i + ls * lda) * COMPSIZE, lda, sb_gemm);
	
	if (n - bk - i > 0) {
	  for (is = i + bk; is < n; is += GEMM_P) {
	    min_i = n - is;
	    if (min_i > GEMM_P) min_i = GEMM_P;
	    
	    if (ls == 0) {
	      NEG_TCOPY (bk, min_i, a + (is + i * lda) * COMPSIZE, lda, sa);

	      TRSM_KERNEL_RT(min_i, bk, bk, dm1, 
#ifdef COMPLEX
			     ZERO,
#endif
			     sa, sa_trsm,
			     a + (is + i * lda) * COMPSIZE, lda, 0);
	    } else {
	      GEMM_ITCOPY (bk, min_i, a + (is + i * lda) * COMPSIZE, lda, sa);
	    }

	    GEMM_KERNEL_N(min_i, min_l, bk, dp1, 
#ifdef COMPLEX
			  ZERO,
#endif
			  sa, sb_gemm,
			  a + (is + ls * lda) * COMPSIZE, lda);
	  }
	}
	
	for (is = 0; is < bk; is += GEMM_P) {
	  min_i = bk - is;
	  if (min_i > GEMM_P) min_i = GEMM_P;
	  
	  TRMM_KERNEL_LT(min_i, min_l, bk, dp1, 
#ifdef COMPLEX
			 ZERO,
#endif
			 sa_trmm + is * bk * COMPSIZE, sb_gemm,
			 a + (i + is + ls * lda) * COMPSIZE, lda, is);
	}
      }

    } else {

      if (n - bk - i > 0) {
	for (is = 0; is < n - bk - i; is += GEMM_P) {
	  min_i = n - bk - i - is;
	  if (min_i > GEMM_P) min_i = GEMM_P;
	  
	  NEG_TCOPY (bk, min_i, a + (i + bk + is + i * lda) * COMPSIZE, lda, sa);
	  
	  TRSM_KERNEL_RT(min_i, bk, bk, dm1, 
#ifdef COMPLEX
			 ZERO,
#endif
			 sa, sa_trsm,
			 a + (i + bk + is + i * lda) * COMPSIZE, lda, 0);
	}
      }
      
    }
  }

  return 0;
}
