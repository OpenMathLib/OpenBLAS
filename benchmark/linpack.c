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

#include "bench.h"

double fabs(double);

#undef GETRF
#undef GETRS

#ifndef COMPLEX
#ifdef XDOUBLE
#define GETRF   BLASFUNC(qgetrf)
#define GETRS   BLASFUNC(qgetrs)
#elif defined(DOUBLE)
#define GETRF   BLASFUNC(dgetrf)
#define GETRS   BLASFUNC(dgetrs)
#else
#define GETRF   BLASFUNC(sgetrf)
#define GETRS   BLASFUNC(sgetrs)
#endif
#else
#ifdef XDOUBLE
#define GETRF   BLASFUNC(xgetrf)
#define GETRS   BLASFUNC(xgetrs)
#elif defined(DOUBLE)
#define GETRF   BLASFUNC(zgetrf)
#define GETRS   BLASFUNC(zgetrs)
#else
#define GETRF   BLASFUNC(cgetrf)
#define GETRS   BLASFUNC(cgetrs)
#endif
#endif

int main(int argc, char *argv[]){

  FLOAT *a, *b;
  blasint *ipiv;

  blasint m, i, j, l, info;
  blasint unit =   1;

  int from =   1;
  int to   = 200;
  int step =   1;
  int loops =  1;

  FLOAT maxerr;

  double time1, time2, timeg1,timeg2;

  char *p;
  if ((p = getenv("OPENBLAS_LOOPS"))) loops=atoi(p);
  
  argc--;argv++;

  if (argc > 0) { from     = atol(*argv);		argc--; argv++;}
  if (argc > 0) { to       = MAX(atol(*argv), from);	argc--; argv++;}
  if (argc > 0) { step     = atol(*argv);		argc--; argv++;}

  fprintf(stderr, "From : %3d  To : %3d Step = %3d\n", from, to, step);

  if (( a = (FLOAT *)malloc(sizeof(FLOAT) * to * to * COMPSIZE)) == NULL){
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }

  if (( b = (FLOAT *)malloc(sizeof(FLOAT) * to * COMPSIZE)) == NULL){
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }

  if (( ipiv = (blasint *)malloc(sizeof(blasint) * to * COMPSIZE)) == NULL){
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }

#ifdef __linux
  srandom(getpid());
#endif

  fprintf(stderr, "   SIZE       Residual     Decompose            Solve           Total\n");

  for(m = from; m <= to; m += step){
    timeg1 = timeg2 = 0.;
    fprintf(stderr, " %6d : ", (int)m);
    for (l = 0; l < loops; l++) {
    for(j = 0; j < m; j++){
      for(i = 0; i < m * COMPSIZE; i++){
	a[(long)i + (long)j * (long)m * COMPSIZE] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
      }
    }

    for (i = 0; i < m * COMPSIZE; ++i) b[i] = 0.;

    for (j = 0; j < m; ++j) {
      for (i = 0; i < m * COMPSIZE; ++i) {
	b[i] += a[(long)i + (long)j * (long)m * COMPSIZE];
      }
    }

    begin();

    GETRF (&m, &m, a, &m, ipiv, &info);

    end();

    if (info) {
      fprintf(stderr, "Matrix is not singular .. %d\n", info);
      exit(1);
    }

    timeg1 += getsec();

    begin();

    GETRS("N", &m, &unit, a, &m, ipiv, b, &m, &info);

    end();

    if (info) {
      fprintf(stderr, "Matrix is not singular .. %d\n", info);
      exit(1);
    }

    timeg2 += getsec();
    } //loops
    time1=timeg1/(double)loops;
    time2=timeg2/(double)loops;
    maxerr = 0.;

    for(i = 0; i < m; i++){
#ifndef XDOUBLE
      if (maxerr < fabs(b[i * COMPSIZE] - 1.0)) maxerr = fabs(b[i * COMPSIZE] - 1.0);
#ifdef COMPLEX
      if (maxerr < fabs(b[i * COMPSIZE] + 1)) maxerr = fabs(b[i * COMPSIZE + 1]);
#endif
#else
      if (maxerr < fabsl(b[i * COMPSIZE] - 1.0L)) maxerr = fabsl(b[i * COMPSIZE] - 1.0L);
#ifdef COMPLEX
      if (maxerr < fabsl(b[i * COMPSIZE] + 1)) maxerr = fabsl(b[i * COMPSIZE + 1]);
#endif
#endif
    }

#ifdef XDOUBLE
    fprintf(stderr,"  %Le ", maxerr);
#else
    fprintf(stderr,"  %e ", maxerr);
#endif

    fprintf(stderr,
	    " %10.2f MFlops %10.2f MFlops %10.2f MFlops\n",
	    COMPSIZE * COMPSIZE * 2. / 3. * (double)m * (double)m * (double)m / time1 * 1.e-6,
	    COMPSIZE * COMPSIZE * 2.      * (double)m * (double)m             / time2 * 1.e-6,
	    COMPSIZE * COMPSIZE * (2. / 3. * (double)m * (double)m * (double)m + 2. * (double)m * (double)m) / (time1 + time2) * 1.e-6);

#if 0
    if (
#ifdef DOUBLE
	maxerr > 1.e-8
#else
	maxerr > 1.e-1
#endif
	) {
      fprintf(stderr, "Error is too large.\n");
      exit(1);
    }
#endif

  }

  return 0;
}

// void main(int argc, char *argv[]) __attribute__((weak, alias("MAIN__")));
