
#include "common.h"

// Some compilers will report feature support for SVE without the appropriate
// header available
#ifdef HAVE_SVE
#if defined __has_include
#if __has_include(<arm_sve.h>) && __ARM_FEATURE_SVE
#define USE_SVE
#endif
#endif
#endif

#include "dgemv_kernel_sve.c"
#include "dgemv_kernel_c.c"


int CNAME(BLASLONG m,  BLASLONG n ,  BLASLONG dummy,  FLOAT  alpha,  FLOAT* a, BLASLONG lda ,  FLOAT *x, BLASLONG incx,  FLOAT *y, BLASLONG incy,  FLOAT *buffer){
  
  if (  incx == 1 && incy == 1){
 // if(alpha!=1) for(BLASLONG i=0; i<n; ++i)X[i]=alpha*X[i]; 
  for(BLASLONG i=0; i<n; ++i){
  // Y[i*incy]+= dgemv_kernel_sve(i,A,lda,X,incx,n);
  y[i]+= dgemv_kernel_sve(i,lda,a,m,x,alpha,n);
 }
  }
//  BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *buffer)
  else  dgemv_kernel_c( m,  n, dummy,  alpha, a, lda,  x,  incx, y,  incy, buffer );
	
 return 0;
}

