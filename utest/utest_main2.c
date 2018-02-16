/*****************************************************************************
Copyright (c) 2011-2016, The OpenBLAS Project
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

**********************************************************************************/

#include <stdio.h>
#include <complex.h>

#define CTEST_MAIN
#define CTEST_SEGFAULT
#define CTEST_ADD_TESTS_MANUALLY

#include "openblas_utest.h"

CTEST(amax, samax){
  blasint N=3, inc=1;
  float te_max=0.0, tr_max=0.0;
  float x[]={-1.1, 2.2, -3.3};
  te_max=BLASFUNC(samax)(&N, x, &inc);
  tr_max=3.3;
  
  ASSERT_DBL_NEAR_TOL((double)(tr_max), (double)(te_max), SINGLE_EPS);
}

CTEST (drotmg,rotmg){
	double te_d1, tr_d1;
	double te_d2, tr_d2;
	double te_x1, tr_x1;
	double te_y1, tr_y1;
	double te_param[5];
	double tr_param[5];
	int i=0;
	// original test case for libGoto bug fixed by feb2014 rewrite
	te_d1= 0.21149573940783739;
	te_d2= 0.046892057172954082;
	te_x1= -0.42272687517106533;
	te_y1= 0.42211309121921659;


	for(i=0; i<5; i++){
	  te_param[i]=tr_param[i]=0.0;
	}

	//reference values as calulated by netlib blas

        tr_d1= 0.1732048;
        tr_d2= 0.03840234;
        tr_x1= -0.516180;
        tr_y1= 0.422113;
        tr_d1= 0.17320483687975;
        tr_d2= 0.03840233915037;
        tr_x1= -0.51618034832329;
        tr_y1= 0.42211309121922;

	tr_param[0]= 0.0;
	tr_param[1]= 0.0;
	tr_param[2]= 0.99854803659786; 
	tr_param[3]= -0.22139439665872;
	tr_param[4]= 0.0;

	BLASFUNC(drotmg)(&te_d1, &te_d2, &te_x1, &te_y1, te_param);
	ASSERT_DBL_NEAR_TOL(te_d1, tr_d1, DOUBLE_EPS);
	ASSERT_DBL_NEAR_TOL(te_d2, tr_d2, DOUBLE_EPS);
	ASSERT_DBL_NEAR_TOL(te_x1, tr_x1, DOUBLE_EPS);
	ASSERT_DBL_NEAR_TOL(te_y1, tr_y1, DOUBLE_EPS);

	for(i=0; i<5; i++){
		ASSERT_DBL_NEAR_TOL(te_param[i], tr_param[i], DOUBLE_EPS);
	}
}

CTEST (drotmg,rotmg_issue1452){
	double te_d1, tr_d1;
	double te_d2, tr_d2;
	double te_x1, tr_x1;
	double te_y1, tr_y1;
	double te_param[5];
	double tr_param[5];
	int i=0;

	// from issue #1452, buggy version returned 0.000244 for param[3]
	te_d1 = 5.9e-8;
	te_d2 = 5.960464e-8;
	te_x1 = 1.0;
	te_y1 = 150.0;

	for(i=0; i<5; i++){
	  te_param[i]=tr_param[i]=0.0;
	}

	//reference values as calulated by netlib blas
	tr_d1= 0.99995592822897;
	tr_d2= 0.98981219860583;
	tr_x1= 0.03662270484346;
	tr_y1= 150.000000000000;

	tr_param[0]= -1.0;
	tr_param[1]= 0.00000161109346;
	tr_param[2]= -0.00024414062500;
	tr_param[3]= 1.0;
	tr_param[4]= 0.00000162760417;

	//OpenBLAS
	BLASFUNC(drotmg)(&te_d1, &te_d2, &te_x1, &te_y1, te_param);

	ASSERT_DBL_NEAR_TOL(te_d1, tr_d1, DOUBLE_EPS);
	ASSERT_DBL_NEAR_TOL(te_d2, tr_d2, DOUBLE_EPS);
	ASSERT_DBL_NEAR_TOL(te_x1, tr_x1, DOUBLE_EPS);
	ASSERT_DBL_NEAR_TOL(te_y1, tr_y1, DOUBLE_EPS);

	for(i=0; i<5; i++){
		ASSERT_DBL_NEAR_TOL(te_param[i], tr_param[i], DOUBLE_EPS);
	}

}

CTEST(drotmg, rotmg_D1eqD2_X1eqX2){
	double te_d1, tr_d1;
	double te_d2, tr_d2;
	double te_x1, tr_x1;
	double te_y1, tr_y1;
	double te_param[5];
	double tr_param[5];
	int i=0;
	te_d1= tr_d1=2.;
	te_d2= tr_d2=2.;
	te_x1= tr_x1=8.;
	te_y1= tr_y1=8.;

	for(i=0; i<5; i++){
	  te_param[i]=tr_param[i]=0.0;
	}
	
	//reference values as calulated by netlib blas
        tr_d1= 1.0;
        tr_d2= 1.0;
        tr_x1= 16.0;
        tr_y1= 8.0;

	tr_param[0]=1.0;
	tr_param[1]=1.0;
	tr_param[2]=0.0;
	tr_param[3]=0.0;
	tr_param[4]=1.0;

	//OpenBLAS
	BLASFUNC(drotmg)(&te_d1, &te_d2, &te_x1, &te_y1, te_param);

	ASSERT_DBL_NEAR_TOL(te_d1, tr_d1, DOUBLE_EPS);
	ASSERT_DBL_NEAR_TOL(te_d2, tr_d2, DOUBLE_EPS);
	ASSERT_DBL_NEAR_TOL(te_x1, tr_x1, DOUBLE_EPS);
	ASSERT_DBL_NEAR_TOL(te_y1, tr_y1, DOUBLE_EPS);

	for(i=0; i<5; i++){
		ASSERT_DBL_NEAR_TOL(te_param[i], tr_param[i], DOUBLE_EPS);
	}
}

CTEST(axpy,daxpy_inc_0)
{
	int i;
	int N=8,incX=0,incY=0;
	double a=0.25;
	double x1[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	double y1[]={2.0,4.0,6.0,8.0,2.0,4.0,6.0,8.0};

	double x2[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	double y2[]={4.0,4.0,6.0,8.0,2.0,4.0,6.0,8.0};

	//OpenBLAS
	BLASFUNC(daxpy)(&N,&a,x1,&incX,y1,&incY);

	for(i=0; i<N; i++){
		ASSERT_DBL_NEAR_TOL(x1[i], x2[i], DOUBLE_EPS);
		ASSERT_DBL_NEAR_TOL(y1[i], y2[i], DOUBLE_EPS);
	}
}

CTEST(axpy,zaxpy_inc_0)
{
	int i;
	int N=4,incX=0,incY=0;
	double a[2]={0.25,0.5};
	double x1[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	double y1[]={2.0,4.0,6.0,8.0,2.0,4.0,6.0,8.0};
	double x2[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	double y2[]={-3.0,9.0,6.0,8.0,2.0,4.0,6.0,8.0};

	//OpenBLAS
	BLASFUNC(zaxpy)(&N,a,x1,&incX,y1,&incY);

	for(i=0; i<2*N; i++){
		ASSERT_DBL_NEAR_TOL(x1[i], x2[i], DOUBLE_EPS);
		ASSERT_DBL_NEAR_TOL(y1[i], y2[i], DOUBLE_EPS);
	}
}

CTEST(axpy,saxpy_inc_0)
{
	int i;
	int N=8,incX=0,incY=0;
	float a=0.25;
	float x1[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	float y1[]={2.0,4.0,6.0,8.0,2.0,4.0,6.0,8.0};
	float x2[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	float y2[]={4.0,4.0,6.0,8.0,2.0,4.0,6.0,8.0};

	//OpenBLAS
	BLASFUNC(saxpy)(&N,&a,x1,&incX,y1,&incY);

	for(i=0; i<N; i++){
		ASSERT_DBL_NEAR_TOL(x1[i], x2[i], DOUBLE_EPS);
		ASSERT_DBL_NEAR_TOL(y1[i], y2[i], DOUBLE_EPS);
	}
}

CTEST(axpy,caxpy_inc_0)
{
	int i;
	int N=4,incX=0,incY=0;
	float a[2]={0.25,0.5};
	float x1[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	float y1[]={2.0,4.0,6.0,8.0,2.0,4.0,6.0,8.0};
	float x2[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	float y2[]={-3.0,9.0,6.0,8.0,2.0,4.0,6.0,8.0};

	//OpenBLAS
	BLASFUNC(caxpy)(&N,a,x1,&incX,y1,&incY);

	for(i=0; i<2*N; i++){
		ASSERT_DBL_NEAR_TOL(x1[i], x2[i], DOUBLE_EPS);
		ASSERT_DBL_NEAR_TOL(y1[i], y2[i], DOUBLE_EPS);
	}
}

CTEST( zdotu,zdotu_n_1)
{
	int N=1,incX=1,incY=1;
	double x1[]={1.0,1.0};
	double y1[]={1.0,2.0};
	double x2[]={1.0,1.0};
	double y2[]={1.0,2.0};
	_Complex result1=0.0;
        _Complex result2={-1.0000+3.0000*I};
	//OpenBLAS
	result1=BLASFUNC(zdotu)(&N,x1,&incX,y1,&incY);

	ASSERT_DBL_NEAR_TOL(creal(result1), creal(result2), DOUBLE_EPS);
	ASSERT_DBL_NEAR_TOL(cimag(result1), cimag(result2), DOUBLE_EPS);

}

CTEST(zdotu, zdotu_offset_1)
{
	int N=1,incX=1,incY=1;
	double x1[]={1.0,2.0,3.0,4.0};
	double y1[]={5.0,6.0,7.0,8.0};
	double x2[]={1.0,2.0,3.0,4.0};
	double y2[]={5.0,6.0,7.0,8.0};
	_Complex result1=0.0;
	_Complex result2={-9.0000+32.0000*I};
	//OpenBLAS
	result1=BLASFUNC(zdotu)(&N,x1+1,&incX,y1+1,&incY);

	ASSERT_DBL_NEAR_TOL(creal(result1), creal(result2), DOUBLE_EPS);
	ASSERT_DBL_NEAR_TOL(cimag(result1), cimag(result2), DOUBLE_EPS);

}

CTEST(dsdot,dsdot_n_1)
{
	float x= 0.172555164F;
	float y= -0.0138700781F;
	int incx=1;
	int incy=1;
	int n=1;

	double res1=0.0f, res2=-0.00239335360107;

	res1=BLASFUNC(dsdot)(&n, &x, &incx, &y, &incy);
	ASSERT_DBL_NEAR_TOL(res1, res2, DOUBLE_EPS);

}

CTEST(rot,drot_inc_0)
{
	int i=0;
	int N=4,incX=0,incY=0;
	double c=0.25,s=0.5;
	double x1[]={1.0,3.0,5.0,7.0};
	double y1[]={2.0,4.0,6.0,8.0};
	double x2[]={-0.21484375000000,3.0,5.0,7.0};
	double y2[]={ 0.03906250000000,4.0,6.0,8.0};


	//OpenBLAS
	BLASFUNC(drot)(&N,x1,&incX,y1,&incY,&c,&s);

	for(i=0; i<N; i++){
		ASSERT_DBL_NEAR_TOL(x1[i], x2[i], DOUBLE_EPS);
		ASSERT_DBL_NEAR_TOL(y1[i], y2[i], DOUBLE_EPS);
	}
}

CTEST(rot,zdrot_inc_0)
{
	int i=0;
	int N=4,incX=0,incY=0;
	double c=0.25,s=0.5;
	double x1[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	double y1[]={2.0,4.0,6.0,8.0,2.0,4.0,6.0,8.0};
	double x2[]={-0.21484375000000,-0.45703125000000 ,5.0,7.0,1.0,3.0,5.0,7.0};
	double y2[]={ 0.03906250000000, 0.17187500000000 ,6.0,8.0,2.0,4.0,6.0,8.0};
	

	//OpenBLAS
	BLASFUNC(zdrot)(&N,x1,&incX,y1,&incY,&c,&s);

	for(i=0; i<2*N; i++){
		ASSERT_DBL_NEAR_TOL(x1[i], x2[i], DOUBLE_EPS);
		ASSERT_DBL_NEAR_TOL(y1[i], y2[i], DOUBLE_EPS);
	}
}

CTEST(rot,srot_inc_0)
{
	int i=0;
	int N=4,incX=0,incY=0;
	float c=0.25,s=0.5;
	float x1[]={1.0,3.0,5.0,7.0};
	float y1[]={2.0,4.0,6.0,8.0};
	float x2[]={-0.21484375000000,3.0,5.0,7.0};
	float y2[]={ 0.03906250000000,4.0,6.0,8.0};

	//OpenBLAS
	BLASFUNC(srot)(&N,x1,&incX,y1,&incY,&c,&s);

	for(i=0; i<N; i++){
		ASSERT_DBL_NEAR_TOL(x1[i], x2[i], SINGLE_EPS);
		ASSERT_DBL_NEAR_TOL(y1[i], y2[i], SINGLE_EPS);
	}
}

CTEST(rot, csrot_inc_0)
{
	int i=0;
	int N=4,incX=0,incY=0;
	float c=0.25,s=0.5;
	float x1[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	float y1[]={2.0,4.0,6.0,8.0,2.0,4.0,6.0,8.0};
	float x2[]={-0.21484375000000,-0.45703125000000 ,5.0,7.0,1.0,3.0,5.0,7.0};
	float y2[]={ 0.03906250000000, 0.17187500000000 ,6.0,8.0,2.0,4.0,6.0,8.0};
	
	//OpenBLAS
	BLASFUNC(csrot)(&N,x1,&incX,y1,&incY,&c,&s);

	for(i=0; i<2*N; i++){
		ASSERT_DBL_NEAR_TOL(x1[i], x2[i], SINGLE_EPS);
		ASSERT_DBL_NEAR_TOL(y1[i], y2[i], SINGLE_EPS);
	}
}

CTEST(swap,dswap_inc_0)
{
	int i=0;
	int N=4,incX=0,incY=0;
	double x1[]={1.0,3.0,5.0,7.0};
	double y1[]={2.0,4.0,6.0,8.0};
	double x2[]={1.0,3.0,5.0,7.0};
	double y2[]={2.0,4.0,6.0,8.0};
	
	//OpenBLAS
	BLASFUNC(dswap)(&N,x1,&incX,y1,&incY);

	for(i=0; i<N; i++){
		ASSERT_DBL_NEAR_TOL(x1[i], x2[i], DOUBLE_EPS);
		ASSERT_DBL_NEAR_TOL(y1[i], y2[i], DOUBLE_EPS);
	}
}

CTEST(swap,zswap_inc_0)
{
	int i=0;
	int N=4,incX=0,incY=0;
	double x1[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	double y1[]={2.0,4.0,6.0,8.0,2.0,4.0,6.0,8.0};
	double x2[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	double y2[]={2.0,4.0,6.0,8.0,2.0,4.0,6.0,8.0};

	//OpenBLAS
	BLASFUNC(zswap)(&N,x1,&incX,y1,&incY);

	for(i=0; i<2*N; i++){
		ASSERT_DBL_NEAR_TOL(x1[i], x2[i], DOUBLE_EPS);
		ASSERT_DBL_NEAR_TOL(y1[i], y2[i], DOUBLE_EPS);
	}
}

CTEST(swap,sswap_inc_0)
{
	int i=0;
	int N=4,incX=0,incY=0;
	float x1[]={1.0,3.0,5.0,7.0};
	float y1[]={2.0,4.0,6.0,8.0};
	float x2[]={1.0,3.0,5.0,7.0};
	float y2[]={2.0,4.0,6.0,8.0};

	//OpenBLAS
	BLASFUNC(sswap)(&N,x1,&incX,y1,&incY);

	for(i=0; i<N; i++){
		ASSERT_DBL_NEAR_TOL(x1[i], x2[i], SINGLE_EPS);
		ASSERT_DBL_NEAR_TOL(y1[i], y2[i], SINGLE_EPS);
	}
}

CTEST(swap,cswap_inc_0)
{
	int i=0;
	int N=4,incX=0,incY=0;
	float x1[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	float y1[]={2.0,4.0,6.0,8.0,2.0,4.0,6.0,8.0};
	float x2[]={1.0,3.0,5.0,7.0,1.0,3.0,5.0,7.0};
	float y2[]={2.0,4.0,6.0,8.0,2.0,4.0,6.0,8.0};

	//OpenBLAS
	BLASFUNC(cswap)(&N,x1,&incX,y1,&incY);

	for(i=0; i<2*N; i++){
		ASSERT_DBL_NEAR_TOL(x1[i], x2[i], SINGLE_EPS);
		ASSERT_DBL_NEAR_TOL(y1[i], y2[i], SINGLE_EPS);
	}
}

int main(int argc, const char ** argv){

  CTEST_ADD(amax, samax);
  int num_fail=0;

  num_fail=ctest_main(argc, argv);

  return num_fail;
}

