/***************************************************************************
Copyright (c) 2014, The OpenBLAS Project
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
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#include "bench.h"

#undef SYR2


#ifdef DOUBLE
#define SYR2   BLASFUNC(dsyr2)
#else
#define SYR2   BLASFUNC(ssyr2)
#endif

int main(int argc, char *argv[]){

  FLOAT *x, *y, *a;
  FLOAT alpha[] = {1.0, 1.0};
  char *p;

  char uplo='U';

  if ((p = getenv("OPENBLAS_UPLO"))) uplo=*p; 

  blasint m, i, j, l;
  blasint inc_x= 1;
  blasint inc_y= 1;
  int from =   1;
  int to   = 200;
  int step =   1;
  int loops =  1;

  if ((p = getenv("OPENBLAS_LOOPS"))) loops=atoi(p);

  double time1,timeg;

  argc--;argv++;

  if (argc > 0) { from     = atol(*argv);		argc--; argv++;}
  if (argc > 0) { to       = MAX(atol(*argv), from);	argc--; argv++;}
  if (argc > 0) { step     = atol(*argv);		argc--; argv++;}

  fprintf(stderr, "From : %3d  To : %3d Step = %3d Uplo = %c Inc_x = %d Inc_y = %d\n", from, to, step,uplo,inc_x,inc_y);

  if (( a = (FLOAT *)malloc(sizeof(FLOAT) * to * to * COMPSIZE)) == NULL){
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }

   if (( x = (FLOAT *)malloc(sizeof(FLOAT) * to * abs(inc_x) * COMPSIZE)) == NULL){
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }
  
   if (( y = (FLOAT *)malloc(sizeof(FLOAT) * to * abs(inc_y) * COMPSIZE)) == NULL){
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }



#ifdef __linux
  srandom(getpid());
#endif

  fprintf(stderr, "   SIZE       Flops\n");

  for(m = from; m <= to; m += step)
  {
    timeg = 0.;	
    fprintf(stderr, " %6d : ", (int)m);
    for (l = 0; l < loops; l++) {
	for(i = 0; i < m * COMPSIZE * abs(inc_x); i++){
	x[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
   	}
	
	for(i = 0; i < m * COMPSIZE * abs(inc_y); i++){
	y[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
   	}

    for(j = 0; j < m; j++){
      for(i = 0; i < m * COMPSIZE; i++){
	a[(long)i + (long)j * (long)m * COMPSIZE] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
      }
    }

    begin();

    SYR2 (&uplo, &m, alpha, x, &inc_x, y, &inc_y, a, &m );

    end();

    timeg += getsec();
    } // loops

    time1 = timeg/(double)loops;
    fprintf(stderr,
	    " %10.2f MFlops\n",
	    COMPSIZE * COMPSIZE * 2. * (double)m * (double)m / time1 * 1.e-6);

  }

  return 0;
}

// void main(int argc, char *argv[]) __attribute__((weak, alias("MAIN__")));
