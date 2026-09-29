/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* Copyright 2025 The OpenBLAS Project.                              */
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

#include <stdio.h>
#include "common.h"

/* Scale the part of C that GEMMT actually references by beta.  This is
   the same operation as syrk_beta() in syrk_k.c. */
static __inline int gemmt_beta(BLASLONG m_from, BLASLONG m_to, BLASLONG n_from, BLASLONG n_to,
			       FLOAT *alpha, FLOAT *c, BLASLONG ldc) {

  BLASLONG i;

#ifndef LOWER
  if (m_from > n_from) n_from = m_from;
  if (m_to   > n_to  ) m_to   = n_to;
#else
  if (m_from < n_from) m_from = n_from;
  if (m_to   < n_to  ) n_to   = m_to;
#endif

  if ((m_to <= m_from) || (n_to <= n_from)) return 0;

  c += (m_from + n_from * ldc) * COMPSIZE;

  m_to -= m_from;
  n_to -= n_from;

  for (i = 0; i < n_to; i++){

#ifndef LOWER

    SCAL_K(MIN(i + n_from - m_from + 1, m_to), 0, 0, alpha[0],
#ifdef COMPLEX
	   alpha[1],
#endif
	   c, 1, NULL, 0, NULL, 0);

    c += ldc * COMPSIZE;

#else

    SCAL_K(MIN(m_to - i + m_from - n_from, m_to), 0, 0, alpha[0],
#ifdef COMPLEX
	 alpha[1],
#endif
	 c, 1, NULL, 0, NULL, 0);

    if (i < m_from - n_from) {
      c += ldc * COMPSIZE;
    } else {
      c += (1 + ldc) * COMPSIZE;
    }
#endif

  }

  return 0;
}

#define GEMMT_BETA(M_FROM, M_TO, N_FROM, N_TO, BETA, C, LDC) \
	gemmt_beta(M_FROM, M_TO, N_FROM, N_TO, BETA, C, LDC)

/* The blocks on the diagonal are computed by the SYRK kernel, which
   clips a block against the diagonal; the conjugation variants needed
   for complex GEMMT are the same file built with CONJA / CONJB. */
#ifndef KERNEL_FUNC
#ifndef LOWER
#if defined(NN) || defined(NT) || defined(TN) || defined(TT)
#define KERNEL_FUNC	SYRK_KERNEL_U
#elif defined(CN) || defined(CT) || defined(RN) || defined(RT)
#define KERNEL_FUNC	GEMMT_KERNEL_UCN
#elif defined(NC) || defined(TC) || defined(NR) || defined(TR)
#define KERNEL_FUNC	GEMMT_KERNEL_UNC
#else
#define KERNEL_FUNC	GEMMT_KERNEL_UCC
#endif
#else
#if defined(NN) || defined(NT) || defined(TN) || defined(TT)
#define KERNEL_FUNC	SYRK_KERNEL_L
#elif defined(CN) || defined(CT) || defined(RN) || defined(RT)
#define KERNEL_FUNC	GEMMT_KERNEL_LCN
#elif defined(NC) || defined(TC) || defined(NR) || defined(TR)
#define KERNEL_FUNC	GEMMT_KERNEL_LNC
#else
#define KERNEL_FUNC	GEMMT_KERNEL_LCC
#endif
#endif
#endif

/* GEMMT is threaded by splitting the triangle into equal-area column
   ranges (syrk_thread) and running this driver on each of them, so there
   is no separate THREADED_LEVEL3 variant. */
#include "level3_gemmt.c"
