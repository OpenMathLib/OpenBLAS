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
#include "cblas.h"
#undef GEMM

#ifndef COMPLEX

#ifdef DOUBLE
#define GEMM   cblas_dgemm
#elif defined(BFLOAT16) && defined(BGEMM)
#define GEMM   cblas_bgemm
#elif defined(BFLOAT16)
#define GEMM   cblas_sbgemm
#undef IFLOAT
#define IFLOAT bfloat16
#elif defined(HFLOAT16)
#define GEMM   cblas_shgemm
#undef IFLOAT
#define IFLOAT hfloat16
#else
#define GEMM   cblas_sgemm
#undef IFLOAT
#define IFLOAT float
#endif

#else

#ifdef DOUBLE
#define GEMM   cblas_zgemm
#else
#define GEMM   cblas_cgemm
#endif

#endif

int main(int argc, char *argv[]){

  IFLOAT *a, *b;
  //IFLOAT *aa, *bb;
  FLOAT *c;
  //FLOAT *cc;
#ifdef BGEMM
  blasint one=1;
  blasint two=2;
  float alpha_in[] = {1.0, 0.0};
  float beta_in[] = {0.0, 0.0};
  FLOAT alpha[2], beta[2];
  sbstobf16_(&two, alpha_in, &one, alpha, &one);
  sbstobf16_(&two, beta_in, &one, beta, &one);
#else
#ifdef COMPLEX
  FLOAT alpha[] = {1.0, 0.0};
  FLOAT beta [] = {0.0, 0.0};
#else
  FLOAT alpha = 1.0;
  FLOAT beta = 0.0;
#endif
#endif
  CBLAS_TRANSPOSE transa = CblasNoTrans;
  CBLAS_TRANSPOSE transb = CblasNoTrans;
 char transac, transbc;
  blasint m, n, k, i, j, lda, ldb, ldc;
  int loops = 1;
  int has_param_m = 0;
  int has_param_n = 0;
  int has_param_k = 0;
  int has_param_lda = 0;
  int has_param_ldb = 0;
  char *p;
//blasint sme=0;
  int from =   1;
  int to   = 200;
  int step =   1;

  double time1, timeg;

  argc--;argv++;

  if (argc > 0) { from = atol(*argv);            argc--; argv++; }
  if (argc > 0) { to   = MAX(atol(*argv), from); argc--; argv++; }
  if (argc > 0) { step = atol(*argv);            argc--; argv++; }

  if ((p = getenv("OPENBLAS_TRANS"))) {
  transa=(*p=='N') ? CblasNoTrans : CblasTrans;
  transb=(*p=='N') ? CblasNoTrans : CblasTrans;
  }
  if ((p = getenv("OPENBLAS_TRANSA"))) {
  transa=(*p=='N') ? CblasNoTrans : CblasTrans;
  }
  if ((p = getenv("OPENBLAS_TRANSB"))) {
  transb=(*p=='N') ? CblasNoTrans : CblasTrans;
  }

  transac=(transa==CblasNoTrans) ? 'N' : 'T';
  transbc=(transb==CblasNoTrans) ? 'N' : 'T';
  fprintf(stderr, "From : %3d  To : %3d Step=%d : Transa=%c : Transb=%c\n", from, to, step, transac, transbc);

  p = getenv("OPENBLAS_LOOPS");
  if ( p != NULL ) {
    loops = atoi(p);
  }

  if ((p = getenv("OPENBLAS_PARAM_M"))) {
    m = atoi(p);
    has_param_m=1;
  } else {
    m = to;
  }
  if ((p = getenv("OPENBLAS_PARAM_N"))) {
    n = atoi(p);
    has_param_n=1;
  } else {
    n = to;
  }
  if ((p = getenv("OPENBLAS_PARAM_K"))) {
    k = atoi(p);
    has_param_k=1;
  } else {
    k = to;
  }
  if ((p = getenv("OPENBLAS_PARAM_LDA"))) {
    lda = atoi(p);
    has_param_lda=1;
  } 
  if ((p = getenv("OPENBLAS_PARAM_LDB"))) {
    ldb = atoi(p);
    has_param_ldb=1;
  } 

  if (( a = (IFLOAT *)malloc(sizeof(IFLOAT) * m * k * COMPSIZE)) == NULL) {
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }
  if (( b = (IFLOAT *)malloc(sizeof(IFLOAT) * k * n * COMPSIZE)) == NULL) {
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }
  if (( c = (FLOAT *)malloc(sizeof(FLOAT) * m * n * COMPSIZE)) == NULL) {
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }
  //if (( aa = (IFLOAT *)malloc(sizeof(IFLOAT) * m * k * COMPSIZE)) == NULL) {
  //  fprintf(stderr,"Out of Memory!!\n");exit(1);
  //}
  //if (( bb = (IFLOAT *)malloc(sizeof(IFLOAT) * k * n * COMPSIZE)) == NULL) {
  //  fprintf(stderr,"Out of Memory!!\n");exit(1);
  //}
  //if (( cc = (FLOAT *)malloc(sizeof(FLOAT) * m * n * COMPSIZE)) == NULL) {
  //  fprintf(stderr,"Out of Memory!!\n");exit(1);
  //}

#ifdef __linux
  srandom(getpid());
#endif

  for (i = 0; i < m * k * COMPSIZE; i++) {
    a[i] = ((IFLOAT) rand() / (IFLOAT) RAND_MAX) - 0.5;
  //  aa[i]=a[i];
  }
  for (i = 0; i < k * n * COMPSIZE; i++) {
    b[i] = ((IFLOAT) rand() / (IFLOAT) RAND_MAX) - 0.5;
  //  bb[i]=b[i];
  }
  for (i = 0; i < m * n * COMPSIZE; i++) {
    c[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
  //  cc[i]=c[i];
  }

  fprintf(stderr, "          SIZE                   Flops             Time\n");

  for (i = from; i <= to; i += step) {
    
    timeg=0;

    if (!has_param_m) { m = i; }
    if (!has_param_n) { n = i; }
    if (!has_param_k) { k = i; }
    
    if (!has_param_lda) {
      if (transa == CblasNoTrans) { lda = k; }
      else { lda = m; }
    }
    if (!has_param_ldb) {
      if (transb == CblasNoTrans) { ldb = n; }
      else { ldb = k; }
    }
    ldc = n;

    fprintf(stderr, " M=%4d, N=%4d, K=%4d : ", (int)m, (int)n, (int)k);
    begin();

    for (j=0; j<loops; j++) {
      GEMM (CblasRowMajor,transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }

 // for (ii = 0; ii < m * n * COMPSIZE; ii++) if (fabsf(c[ii]-cc[ii])>1.5e-5){fprintf(stderr,"mismatch %d %f !=  %f: %g\n",ii,c[ii],cc[ii],fabsf(c[ii]-cc[ii]));}
    end();
    time1 = getsec();

    timeg = time1/loops;
    fprintf(stderr,
	    " %10.2f MFlops %10.6f sec\n",
	    COMPSIZE * COMPSIZE * 2. * (double)k * (double)m * (double)n / timeg * 1.e-6, time1);
    
  }

  return 0;
}

// void main(int argc, char *argv[]) __attribute__((weak, alias("MAIN__")));
