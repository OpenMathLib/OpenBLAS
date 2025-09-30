/***************************************************************************
Copyright (c) 2025, The OpenBLAS Project
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
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#include <stdio.h>
#include "common.h"

#define copysign(x,y) ((y) < 0 ? ((x) < 0 ? (x) : -(x)) : ((x) < 0 ? -(x) : (x)))

#if defined(DOUBLE)
#define LAMC3  BLASFUNC(dlamc3)
#define LAED4  BLASFUNC(dlaed4)
#define GEMM   BLASFUNC(dgemm)
#define NRM2   BLASFUNC(dnrm2)
#define COPY   BLASFUNC(dcopy)
#define LACPY  BLASFUNC(dlacpy)
#define LASET  BLASFUNC(dlaset)
#else
#define LAMC3  BLASFUNC(slamc3)
#define LAED4  BLASFUNC(slaed4)
#define GEMM   BLASFUNC(sgemm)
#define NRM2   BLASFUNC(snrm2)
#define COPY   BLASFUNC(scopy)
#define LACPY  BLASFUNC(slacpy)
#define LASET  BLASFUNC(slaset)
#endif

FLOAT LAMC3(FLOAT *, FLOAT *);
void LAED4(blasint *, blasint *, FLOAT *, FLOAT *, FLOAT *, FLOAT *, FLOAT *, blasint *);
void LACPY(char *, blasint *, blasint *, FLOAT *, blasint *, FLOAT *, blasint *);
void LASET(char *, blasint *, blasint *, FLOAT *, FLOAT *, FLOAT *, blasint *);

/* Table of constant values */
static blasint c1 = 1;
static FLOAT c1f = 1.;
static FLOAT c0f = 0.;

/* ===================================================================== */
blasint CNAME(blasint *k, blasint *n, blasint *n1, FLOAT *d, 
        FLOAT *q, blasint *ldq, FLOAT *rho, FLOAT *dlamda,
        FLOAT *q2, blasint *indx, blasint *ctot, FLOAT *w, 
        FLOAT *s, blasint *info)
{
  FLOAT temp;
  blasint kval, qdim;
  blasint i, j, itmp;
  blasint n2, n12, ii, n23, iq2;

  qdim = *ldq;
  kval = *k;

/*     Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can */
/*     be computed with high relative accuracy (barring over/underflow). */

  for (i = 0; i < kval; i++) {
    dlamda[i] = LAMC3(&dlamda[i], &dlamda[i]) - dlamda[i];
  }

  for (j = 1; j <= kval; j++) {
    LAED4(k, &j, dlamda, w, &q[(j - 1) * qdim], rho, &d[j - 1], info);
    if(*info != 0) break;
  }

/*   If the zero finder fails, the computation is terminated. */

  if (*info != 0) {
    return 0;
  }

  if (kval == 2) {
    for (j = 0; j < kval; j++) {
      w[0] = q[j * qdim];
      w[1] = q[j * qdim + 1];
      ii = indx[0] - 1;
      q[j * qdim] = w[ii];
      ii = indx[1] - 1;
      q[j * qdim + 1] = w[ii];
    }
  } else if (kval != 1) {

/*   Compute updated W. */

    COPY(k, w, &c1, s, &c1);

/*   Initialize W(I) = Q(I,I) */

    itmp = qdim + 1;
    COPY(k, q, &itmp, w, &c1);
    for (j = 0; j < kval; j++) {
      for (i = 0; i < j; i++) {
        w[i] *= q[j * qdim + i] / (dlamda[i] - dlamda[j]);
      }
      for (i = j + 1; i < kval; i++) {
        w[i] *= q[j * qdim + i] / (dlamda[i] - dlamda[j]);
      }
    }
    for (i = 0; i < kval; i++) {
      temp = sqrt(-w[i]);
      w[i] = copysign(temp, s[i]);
    }

/*   Compute eigenvectors of the modified rank-1 modification. */

    for (j = 0; j < kval; j++) {
      for (i = 0; i < kval; i++) {
        s[i] = w[i] / q[j * qdim + i];
      }
      temp = NRM2(k, s, &c1);
      for (i = 0; i < kval; i++) {
        ii = indx[i] - 1;
        q[j * qdim + i] = s[ii] / temp;
      }
    }
  }

/*   Compute the updated eigenvectors. */

  n2 = *n - *n1;
  n12 = ctot[0] + ctot[1];
  n23 = ctot[1] + ctot[2];

  LACPY("A", &n23, k, &q[ctot[0]], ldq, s, &n23);
  iq2 = *n1 * n12;
  if (n23 != 0) {
    GEMM("N", "N", &n2, k, &n23, &c1f, &q2[iq2], &n2, s, &n23, &c0f, &q[*n1], ldq);
  } else {
    LASET("A", &n2, k, &c0f, &c0f, &q[*n1], ldq);
  }

  LACPY("A", &n12, k, q, ldq, s, &n12);
  if (n12 != 0) {
    GEMM("N", "N", n1, k, &n12, &c1f, q2, n1, s, &n12, &c0f, q, ldq);
  } else {
    LASET("A", n1, k, &c0f, &c0f, q, ldq);
  }

  return 0;
}
