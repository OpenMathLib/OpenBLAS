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
/* interpreted as representing official policies, either expressed   */
/* or implied, of The University of Texas at Austin.                 */
/*********************************************************************/

#ifndef ASSEMBLER

#ifdef __cplusplus
extern "C" {
	/* Assume C declarations for C++ */
#endif  /* __cplusplus */

int    BLASFUNC(xerbla)(char *, blasint *info, blasint);

OPENBLAS_EXPORT void extern openblas_set_num_threads_(int *);

/*Set the threading backend to a custom callback.*/
typedef void (*openblas_dojob_callback)(int thread_num, void *jobdata, int dojob_data);
typedef void (*openblas_threads_callback)(int sync, openblas_dojob_callback dojob, int numjobs, size_t jobdata_elsize, void *jobdata, int dojob_data);
extern openblas_threads_callback openblas_threads_callback_;

FLOATRET  OPENBLAS_API(sdot)  (blasint *, float  *, blasint *, float  *, blasint *);
FLOATRET  OPENBLAS_API(sdsdot)(blasint *, float  *,        float  *, blasint *, float  *, blasint *);

double OPENBLAS_API(dsdot) (blasint *, float  *, blasint *, float  *, blasint *);
double OPENBLAS_API(ddot)  (blasint *, double *, blasint *, double *, blasint *);
xdouble OPENBLAS_API(qdot)  (blasint *, xdouble *, blasint *, xdouble *, blasint *);

void   OPENBLAS_API(bscal) (blasint *,  bfloat16  *, bfloat16  *, blasint *);
float  OPENBLAS_API(sbdot)     (blasint *, bfloat16 *, blasint *, bfloat16 *, blasint *);
void   OPENBLAS_API(sbstobf16) (blasint *, float *,    blasint *, bfloat16 *, blasint *);
void   OPENBLAS_API(sbdtobf16) (blasint *, double *,   blasint *, bfloat16 *, blasint *);
void   OPENBLAS_API(sbf16tos)  (blasint *, bfloat16 *, blasint *, float *,    blasint *);
void   OPENBLAS_API(dbf16tod)  (blasint *, bfloat16 *, blasint *, double *,   blasint *);

#ifdef RETURN_BY_STRUCT
typedef struct {
  float r, i;
} myccomplex_t;

typedef struct {
  double r, i;
} myzcomplex_t;

typedef struct {
  xdouble r, i;
} myxcomplex_t;

myccomplex_t    OPENBLAS_API(cdotu)  (blasint *, float  *, blasint *, float  *, blasint *);
myccomplex_t    OPENBLAS_API(cdotc)  (blasint *, float  *, blasint *, float  *, blasint *);
myzcomplex_t    OPENBLAS_API(zdotu)  (blasint *, double  *, blasint *, double  *, blasint *);
myzcomplex_t    OPENBLAS_API(zdotc)  (blasint *, double  *, blasint *, double  *, blasint *);
myxcomplex_t    OPENBLAS_API(xdotu)  (blasint *, xdouble  *, blasint *, xdouble  *, blasint *);
myxcomplex_t    OPENBLAS_API(xdotc)  (blasint *, xdouble  *, blasint *, xdouble  *, blasint *);

#elif defined RETURN_BY_STACK
void  OPENBLAS_API(cdotu)  (openblas_complex_float   *,  blasint *, float  * , blasint *, float  *,  blasint *);
void  OPENBLAS_API(cdotc)  (openblas_complex_float   *,  blasint *, float  *,  blasint *, float  *,  blasint *);
void  OPENBLAS_API(zdotu)  (openblas_complex_double  *, blasint *, double  *, blasint *, double  *, blasint *);
void  OPENBLAS_API(zdotc)  (openblas_complex_double  *, blasint *, double  *, blasint *, double  *, blasint *);
void  OPENBLAS_API(xdotu)  (openblas_complex_xdouble *, blasint *, xdouble  *, blasint *, xdouble  *, blasint *);
void  OPENBLAS_API(xdotc)  (openblas_complex_xdouble *, blasint *, xdouble  *, blasint *, xdouble  *, blasint *);
#else
openblas_complex_float   OPENBLAS_API(cdotu)  (blasint *, float  *, blasint *, float  *, blasint *);
openblas_complex_float   OPENBLAS_API(cdotc)  (blasint *, float  *, blasint *, float  *, blasint *);
openblas_complex_double  OPENBLAS_API(zdotu)  (blasint *, double  *, blasint *, double  *, blasint *);
openblas_complex_double  OPENBLAS_API(zdotc)  (blasint *, double  *, blasint *, double  *, blasint *);
openblas_complex_xdouble OPENBLAS_API(xdotu)  (blasint *, xdouble  *, blasint *, xdouble  *, blasint *);
openblas_complex_xdouble OPENBLAS_API(xdotc)  (blasint *, xdouble  *, blasint *, xdouble  *, blasint *);
#endif

void    OPENBLAS_API(saxpy) (blasint *, float  *, float  *, blasint *, float  *, blasint *);
void    OPENBLAS_API(daxpy) (blasint *, double *, double *, blasint *, double *, blasint *);
void    OPENBLAS_API(qaxpy) (blasint *, xdouble *, xdouble *, blasint *, xdouble *, blasint *);
void    OPENBLAS_API(caxpy) (blasint *, float  *, float  *, blasint *, float  *, blasint *);
void    OPENBLAS_API(zaxpy) (blasint *, double *, double *, blasint *, double *, blasint *);
void    OPENBLAS_API(xaxpy) (blasint *, xdouble *, xdouble *, blasint *, xdouble *, blasint *);
void    OPENBLAS_API(caxpyc)(blasint *, float  *, float  *, blasint *, float  *, blasint *);
void    OPENBLAS_API(zaxpyc)(blasint *, double *, double *, blasint *, double *, blasint *);
void    OPENBLAS_API(xaxpyc)(blasint *, xdouble *, xdouble *, blasint *, xdouble *, blasint *);

void    OPENBLAS_API(scopy) (blasint *, float  *, blasint *, float  *, blasint *);
void    OPENBLAS_API(dcopy) (blasint *, double *, blasint *, double *, blasint *);
void    OPENBLAS_API(qcopy) (blasint *, xdouble *, blasint *, xdouble *, blasint *);
void    OPENBLAS_API(ccopy) (blasint *, float  *, blasint *, float  *, blasint *);
void    OPENBLAS_API(zcopy) (blasint *, double *, blasint *, double *, blasint *);
void    OPENBLAS_API(xcopy) (blasint *, xdouble *, blasint *, xdouble *, blasint *);

void    OPENBLAS_API(sswap) (blasint *, float  *, blasint *, float  *, blasint *);
void    OPENBLAS_API(dswap) (blasint *, double *, blasint *, double *, blasint *);
void    OPENBLAS_API(qswap) (blasint *, xdouble *, blasint *, xdouble *, blasint *);
void    OPENBLAS_API(cswap) (blasint *, float  *, blasint *, float  *, blasint *);
void    OPENBLAS_API(zswap) (blasint *, double *, blasint *, double *, blasint *);
void    OPENBLAS_API(xswap) (blasint *, xdouble *, blasint *, xdouble *, blasint *);

FLOATRET  OPENBLAS_API(sasum) (blasint *, float  *, blasint *);
FLOATRET  OPENBLAS_API(scasum)(blasint *, float  *, blasint *);
double OPENBLAS_API(dasum) (blasint *, double *, blasint *);
xdouble OPENBLAS_API(qasum) (blasint *, xdouble *, blasint *);
double OPENBLAS_API(dzasum)(blasint *, double *, blasint *);
xdouble OPENBLAS_API(qxasum)(blasint *, xdouble *, blasint *);

FLOATRET  OPENBLAS_API(ssum) (blasint *, float  *, blasint *);
FLOATRET  OPENBLAS_API(scsum)(blasint *, float  *, blasint *);
double OPENBLAS_API(dsum) (blasint *, double *, blasint *);
xdouble OPENBLAS_API(qsum) (blasint *, xdouble *, blasint *);
double OPENBLAS_API(dzsum)(blasint *, double *, blasint *);
xdouble OPENBLAS_API(qxsum)(blasint *, xdouble *, blasint *);

blasint    OPENBLAS_API(isamax)(blasint *, float  *, blasint *);
blasint    OPENBLAS_API(idamax)(blasint *, double *, blasint *);
blasint    OPENBLAS_API(iqamax)(blasint *, xdouble *, blasint *);
blasint    OPENBLAS_API(icamax)(blasint *, float  *, blasint *);
blasint    OPENBLAS_API(izamax)(blasint *, double *, blasint *);
blasint    OPENBLAS_API(ixamax)(blasint *, xdouble *, blasint *);

blasint    OPENBLAS_API(ismax) (blasint *, float  *, blasint *);
blasint    OPENBLAS_API(idmax) (blasint *, double *, blasint *);
blasint    OPENBLAS_API(iqmax) (blasint *, xdouble *, blasint *);
blasint    OPENBLAS_API(icmax) (blasint *, float  *, blasint *);
blasint    OPENBLAS_API(izmax) (blasint *, double *, blasint *);
blasint    OPENBLAS_API(ixmax) (blasint *, xdouble *, blasint *);

blasint    OPENBLAS_API(isamin)(blasint *, float  *, blasint *);
blasint    OPENBLAS_API(idamin)(blasint *, double *, blasint *);
blasint    OPENBLAS_API(iqamin)(blasint *, xdouble *, blasint *);
blasint    OPENBLAS_API(icamin)(blasint *, float  *, blasint *);
blasint    OPENBLAS_API(izamin)(blasint *, double *, blasint *);
blasint    OPENBLAS_API(ixamin)(blasint *, xdouble *, blasint *);

blasint    OPENBLAS_API(ismin)(blasint *, float  *, blasint *);
blasint    OPENBLAS_API(idmin)(blasint *, double *, blasint *);
blasint    OPENBLAS_API(iqmin)(blasint *, xdouble *, blasint *);
blasint    OPENBLAS_API(icmin)(blasint *, float  *, blasint *);
blasint    OPENBLAS_API(izmin)(blasint *, double *, blasint *);
blasint    OPENBLAS_API(ixmin)(blasint *, xdouble *, blasint *);

FLOATRET  OPENBLAS_API(samax) (blasint *, float  *, blasint *);
double OPENBLAS_API(damax) (blasint *, double *, blasint *);
xdouble OPENBLAS_API(qamax) (blasint *, xdouble *, blasint *);
FLOATRET  OPENBLAS_API(scamax)(blasint *, float  *, blasint *);
double OPENBLAS_API(dzamax)(blasint *, double *, blasint *);
xdouble OPENBLAS_API(qxamax)(blasint *, xdouble *, blasint *);

FLOATRET  OPENBLAS_API(samin) (blasint *, float  *, blasint *);
double OPENBLAS_API(damin) (blasint *, double *, blasint *);
xdouble OPENBLAS_API(qamin) (blasint *, xdouble *, blasint *);
FLOATRET  OPENBLAS_API(scamin)(blasint *, float  *, blasint *);
double OPENBLAS_API(dzamin)(blasint *, double *, blasint *);
xdouble OPENBLAS_API(qxamin)(blasint *, xdouble *, blasint *);

FLOATRET  OPENBLAS_API(smax)  (blasint *, float  *, blasint *);
double OPENBLAS_API(dmax)  (blasint *, double *, blasint *);
xdouble OPENBLAS_API(qmax)  (blasint *, xdouble *, blasint *);
FLOATRET  OPENBLAS_API(scmax) (blasint *, float  *, blasint *);
double OPENBLAS_API(dzmax) (blasint *, double *, blasint *);
xdouble OPENBLAS_API(qxmax) (blasint *, xdouble *, blasint *);

FLOATRET  OPENBLAS_API(smin)  (blasint *, float  *, blasint *);
double OPENBLAS_API(dmin)  (blasint *, double *, blasint *);
xdouble OPENBLAS_API(qmin)  (blasint *, xdouble *, blasint *);
FLOATRET  OPENBLAS_API(scmin) (blasint *, float  *, blasint *);
double OPENBLAS_API(dzmin) (blasint *, double *, blasint *);
xdouble OPENBLAS_API(qxmin) (blasint *, xdouble *, blasint *);

void    OPENBLAS_API(sscal) (blasint *,  float  *, float  *, blasint *);
void    OPENBLAS_API(dscal) (blasint *,  double *, double *, blasint *);
void    OPENBLAS_API(qscal) (blasint *,  xdouble *, xdouble *, blasint *);
void    OPENBLAS_API(cscal) (blasint *,  float  *, float  *, blasint *);
void    OPENBLAS_API(zscal) (blasint *,  double *, double *, blasint *);
void    OPENBLAS_API(xscal) (blasint *,  xdouble *, xdouble *, blasint *);
void    OPENBLAS_API(csscal)(blasint *,  float  *, float  *, blasint *);
void    OPENBLAS_API(zdscal)(blasint *,  double *, double *, blasint *);
void    OPENBLAS_API(xqscal)(blasint *,  xdouble *, xdouble *, blasint *);

FLOATRET  OPENBLAS_API(snrm2) (blasint *, float  *, blasint *);
FLOATRET  OPENBLAS_API(scnrm2)(blasint *, float  *, blasint *);

double OPENBLAS_API(dnrm2) (blasint *, double *, blasint *);
xdouble OPENBLAS_API(qnrm2) (blasint *, xdouble *, blasint *);
double OPENBLAS_API(dznrm2)(blasint *, double *, blasint *);
xdouble OPENBLAS_API(qxnrm2)(blasint *, xdouble *, blasint *);

void  OPENBLAS_API(srot)  (blasint *, float  *, blasint *, float  *, blasint *, float  *, float  *);
void  OPENBLAS_API(drot)  (blasint *, double *, blasint *, double *, blasint *, double *, double *);
void  OPENBLAS_API(qrot)  (blasint *, xdouble *, blasint *, xdouble *, blasint *, xdouble *, xdouble *);
void  OPENBLAS_API(csrot) (blasint *, float  *, blasint *, float  *, blasint *, float  *, float  *);
void  OPENBLAS_API(zdrot) (blasint *, double *, blasint *, double *, blasint *, double *, double *);
void  OPENBLAS_API(xqrot) (blasint *, xdouble *, blasint *, xdouble *, blasint *, xdouble *, xdouble *);

void  OPENBLAS_API(srotg) (float  *, float  *, float  *, float  *);
void  OPENBLAS_API(drotg) (double *, double *, double *, double *);
void  OPENBLAS_API(qrotg) (xdouble *, xdouble *, xdouble *, xdouble *);
void  OPENBLAS_API(crotg) (float  *, float  *, float  *, float  *);
void  OPENBLAS_API(zrotg) (double *, double *, double *, double *);
void  OPENBLAS_API(xrotg) (xdouble *, xdouble *, xdouble *, xdouble *);

void  OPENBLAS_API(srotmg)(float  *, float  *, float  *, float  *, float  *);
void  OPENBLAS_API(drotmg)(double *, double *, double *, double *, double *);

void  OPENBLAS_API(srotm) (blasint *, float  *, blasint *, float  *, blasint *, float  *);
void  OPENBLAS_API(drotm) (blasint *, double *, blasint *, double *, blasint *, double *);
void  OPENBLAS_API(qrotm) (blasint *, xdouble *, blasint *, xdouble *, blasint *, xdouble *);

/* Level 2 routines */

void OPENBLAS_API(sger)(blasint *,    blasint *, float *,  float *, blasint *,
		   float *,  blasint *, float *,  blasint *);
void OPENBLAS_API(dger)(blasint *,    blasint *, double *, double *, blasint *,
		   double *, blasint *, double *, blasint *);
void OPENBLAS_API(qger)(blasint *,    blasint *, xdouble *, xdouble *, blasint *,
		   xdouble *, blasint *, xdouble *, blasint *);
void OPENBLAS_API(cgeru)(blasint *,    blasint *, float *,  float *, blasint *,
		    float *,  blasint *, float *,  blasint *);
void OPENBLAS_API(cgerc)(blasint *,    blasint *, float *,  float *, blasint *,
		    float *,  blasint *, float *,  blasint *);
void OPENBLAS_API(zgeru)(blasint *,    blasint *, double *, double *, blasint *,
		    double *, blasint *, double *, blasint *);
void OPENBLAS_API(zgerc)(blasint *,    blasint *, double *, double *, blasint *,
		    double *, blasint *, double *, blasint *);
void OPENBLAS_API(xgeru)(blasint *,    blasint *, xdouble *, xdouble *, blasint *,
		    xdouble *, blasint *, xdouble *, blasint *);
void OPENBLAS_API(xgerc)(blasint *,    blasint *, xdouble *, xdouble *, blasint *,
		    xdouble *, blasint *, xdouble *, blasint *);

void OPENBLAS_API(bgemv)(char *, blasint *, blasint *, bfloat16  *, bfloat16 *, blasint *,
            bfloat16  *, blasint *, bfloat16  *, bfloat16  *, blasint *);
void OPENBLAS_API(sbgemv)(char *, blasint *, blasint *, float  *, bfloat16 *, blasint *,
            bfloat16  *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(shgemv)(char *, blasint *, blasint *, float  *, hfloat16 *, blasint *,
            hfloat16  *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(sgemv)(char *, blasint *, blasint *, float  *, float  *, blasint *,
		    float  *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(dgemv)(char *, blasint *, blasint *, double *, double *, blasint *,
		    double *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(qgemv)(char *, blasint *, blasint *, xdouble *, xdouble *, blasint *,
		    xdouble *, blasint *, xdouble *, xdouble *, blasint *);
void OPENBLAS_API(cgemv)(char *, blasint *, blasint *, float  *, float  *, blasint *,
		    float  *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(zgemv)(char *, blasint *, blasint *, double *, double *, blasint *,
		    double *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(xgemv)(char *, blasint *, blasint *, xdouble *, xdouble *, blasint *,
		    xdouble *, blasint *, xdouble *, xdouble *, blasint *);

void OPENBLAS_API(strsv) (char *, char *, char *, blasint *, float  *, blasint *,
		     float  *, blasint *);
void OPENBLAS_API(dtrsv) (char *, char *, char *, blasint *, double *, blasint *,
		     double *, blasint *);
void OPENBLAS_API(qtrsv) (char *, char *, char *, blasint *, xdouble *, blasint *,
		     xdouble *, blasint *);
void OPENBLAS_API(ctrsv) (char *, char *, char *, blasint *, float  *, blasint *,
		     float  *, blasint *);
void OPENBLAS_API(ztrsv) (char *, char *, char *, blasint *, double *, blasint *,
		     double *, blasint *);
void OPENBLAS_API(xtrsv) (char *, char *, char *, blasint *, xdouble *, blasint *,
		     xdouble *, blasint *);

void OPENBLAS_API(strmv) (char *, char *, char *, blasint *, float  *, blasint *,
		     float  *, blasint *);
void OPENBLAS_API(dtrmv) (char *, char *, char *, blasint *, double *, blasint *,
		     double *, blasint *);
void OPENBLAS_API(qtrmv) (char *, char *, char *, blasint *, xdouble *, blasint *,
		     xdouble *, blasint *);
void OPENBLAS_API(ctrmv) (char *, char *, char *, blasint *, float  *, blasint *,
		     float  *, blasint *);
void OPENBLAS_API(ztrmv) (char *, char *, char *, blasint *, double *, blasint *,
		     double *, blasint *);
void OPENBLAS_API(xtrmv) (char *, char *, char *, blasint *, xdouble *, blasint *,
		     xdouble *, blasint *);

void OPENBLAS_API(stpsv) (char *, char *, char *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(dtpsv) (char *, char *, char *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(qtpsv) (char *, char *, char *, blasint *, xdouble *, xdouble *, blasint *);
void OPENBLAS_API(ctpsv) (char *, char *, char *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(ztpsv) (char *, char *, char *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(xtpsv) (char *, char *, char *, blasint *, xdouble *, xdouble *, blasint *);

void OPENBLAS_API(stpmv) (char *, char *, char *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(dtpmv) (char *, char *, char *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(qtpmv) (char *, char *, char *, blasint *, xdouble *, xdouble *, blasint *);
void OPENBLAS_API(ctpmv) (char *, char *, char *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(ztpmv) (char *, char *, char *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(xtpmv) (char *, char *, char *, blasint *, xdouble *, xdouble *, blasint *);

void OPENBLAS_API(stbmv) (char *, char *, char *, blasint *, blasint *, float  *, blasint *, float  *, blasint *);
void OPENBLAS_API(dtbmv) (char *, char *, char *, blasint *, blasint *, double *, blasint *, double *, blasint *);
void OPENBLAS_API(qtbmv) (char *, char *, char *, blasint *, blasint *, xdouble *, blasint *, xdouble *, blasint *);
void OPENBLAS_API(ctbmv) (char *, char *, char *, blasint *, blasint *, float  *, blasint *, float  *, blasint *);
void OPENBLAS_API(ztbmv) (char *, char *, char *, blasint *, blasint *, double *, blasint *, double *, blasint *);
void OPENBLAS_API(xtbmv) (char *, char *, char *, blasint *, blasint *, xdouble *, blasint *, xdouble *, blasint *);

void OPENBLAS_API(stbsv) (char *, char *, char *, blasint *, blasint *, float  *, blasint *, float  *, blasint *);
void OPENBLAS_API(dtbsv) (char *, char *, char *, blasint *, blasint *, double *, blasint *, double *, blasint *);
void OPENBLAS_API(qtbsv) (char *, char *, char *, blasint *, blasint *, xdouble *, blasint *, xdouble *, blasint *);
void OPENBLAS_API(ctbsv) (char *, char *, char *, blasint *, blasint *, float  *, blasint *, float  *, blasint *);
void OPENBLAS_API(ztbsv) (char *, char *, char *, blasint *, blasint *, double *, blasint *, double *, blasint *);
void OPENBLAS_API(xtbsv) (char *, char *, char *, blasint *, blasint *, xdouble *, blasint *, xdouble *, blasint *);

void OPENBLAS_API(ssymv) (char *, blasint *, float  *, float *, blasint *,
		     float  *, blasint *, float *, float *, blasint *);
void OPENBLAS_API(dsymv) (char *, blasint *, double  *, double *, blasint *,
		     double  *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(qsymv) (char *, blasint *, xdouble  *, xdouble *, blasint *,
		     xdouble  *, blasint *, xdouble *, xdouble *, blasint *);
void OPENBLAS_API(csymv) (char *, blasint *, float  *, float *, blasint *,
		     float  *, blasint *, float *, float *, blasint *);
void OPENBLAS_API(zsymv) (char *, blasint *, double  *, double *, blasint *,
		     double  *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(xsymv) (char *, blasint *, xdouble  *, xdouble *, blasint *,
		     xdouble  *, blasint *, xdouble *, xdouble *, blasint *);

void OPENBLAS_API(sspmv) (char *, blasint *, float  *, float *,
		     float  *, blasint *, float *, float *, blasint *);
void OPENBLAS_API(dspmv) (char *, blasint *, double  *, double *,
		     double  *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(qspmv) (char *, blasint *, xdouble  *, xdouble *,
		     xdouble  *, blasint *, xdouble *, xdouble *, blasint *);
void OPENBLAS_API(cspmv) (char *, blasint *, float  *, float *,
		     float  *, blasint *, float *, float *, blasint *);
void OPENBLAS_API(zspmv) (char *, blasint *, double  *, double *,
		     double  *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(xspmv) (char *, blasint *, xdouble  *, xdouble *,
		     xdouble  *, blasint *, xdouble *, xdouble *, blasint *);

void OPENBLAS_API(ssyr) (char *, blasint *, float   *, float  *, blasint *,
		    float  *, blasint *);
void OPENBLAS_API(dsyr) (char *, blasint *, double  *, double *, blasint *,
		    double *, blasint *);
void OPENBLAS_API(qsyr) (char *, blasint *, xdouble  *, xdouble *, blasint *,
		    xdouble *, blasint *);
void OPENBLAS_API(csyr) (char *, blasint *, float   *, float  *, blasint *,
		    float  *, blasint *);
void OPENBLAS_API(zsyr) (char *, blasint *, double  *, double *, blasint *,
		    double *, blasint *);
void OPENBLAS_API(xsyr) (char *, blasint *, xdouble  *, xdouble *, blasint *,
		    xdouble *, blasint *);

void OPENBLAS_API(ssyr2) (char *, blasint *, float   *,
		     float  *, blasint *, float  *, blasint *, float  *, blasint *);
void OPENBLAS_API(dsyr2) (char *, blasint *, double  *,
		     double *, blasint *, double *, blasint *, double *, blasint *);
void OPENBLAS_API(qsyr2) (char *, blasint *, xdouble  *,
		     xdouble *, blasint *, xdouble *, blasint *, xdouble *, blasint *);
void OPENBLAS_API(csyr2) (char *, blasint *, float   *,
		     float  *, blasint *, float  *, blasint *, float  *, blasint *);
void OPENBLAS_API(zsyr2) (char *, blasint *, double  *,
		     double *, blasint *, double *, blasint *, double *, blasint *);
void OPENBLAS_API(xsyr2) (char *, blasint *, xdouble  *,
		     xdouble *, blasint *, xdouble *, blasint *, xdouble *, blasint *);

void OPENBLAS_API(sspr) (char *, blasint *, float   *, float  *, blasint *,
		    float  *);
void OPENBLAS_API(dspr) (char *, blasint *, double  *, double *, blasint *,
		    double *);
void OPENBLAS_API(qspr) (char *, blasint *, xdouble  *, xdouble *, blasint *,
		    xdouble *);
void OPENBLAS_API(cspr) (char *, blasint *, float   *, float  *, blasint *,
		    float  *);
void OPENBLAS_API(zspr) (char *, blasint *, double  *, double *, blasint *,
		    double *);
void OPENBLAS_API(xspr) (char *, blasint *, xdouble  *, xdouble *, blasint *,
		    xdouble *);

void OPENBLAS_API(sspr2) (char *, blasint *, float   *,
		     float  *, blasint *, float  *, blasint *, float  *);
void OPENBLAS_API(dspr2) (char *, blasint *, double  *,
		     double *, blasint *, double *, blasint *, double *);
void OPENBLAS_API(qspr2) (char *, blasint *, xdouble  *,
		     xdouble *, blasint *, xdouble *, blasint *, xdouble *);
void OPENBLAS_API(cspr2) (char *, blasint *, float   *,
		     float  *, blasint *, float  *, blasint *, float  *);
void OPENBLAS_API(zspr2) (char *, blasint *, double  *,
		     double *, blasint *, double *, blasint *, double *);
void OPENBLAS_API(xspr2) (char *, blasint *, xdouble  *,
		     xdouble *, blasint *, xdouble *, blasint *, xdouble *);

void OPENBLAS_API(cher) (char *, blasint *, float   *, float  *, blasint *,
		    float  *, blasint *);
void OPENBLAS_API(zher) (char *, blasint *, double  *, double *, blasint *,
		    double *, blasint *);
void OPENBLAS_API(xher) (char *, blasint *, xdouble  *, xdouble *, blasint *,
		    xdouble *, blasint *);

void OPENBLAS_API(chpr) (char *, blasint *, float   *, float  *, blasint *, float  *);
void OPENBLAS_API(zhpr) (char *, blasint *, double  *, double *, blasint *, double *);
void OPENBLAS_API(xhpr) (char *, blasint *, xdouble  *, xdouble *, blasint *, xdouble *);

void OPENBLAS_API(cher2) (char *, blasint *, float   *,
		     float  *, blasint *, float  *, blasint *, float  *, blasint *);
void OPENBLAS_API(zher2) (char *, blasint *, double  *,
		     double *, blasint *, double *, blasint *, double *, blasint *);
void OPENBLAS_API(xher2) (char *, blasint *, xdouble  *,
		     xdouble *, blasint *, xdouble *, blasint *, xdouble *, blasint *);

void OPENBLAS_API(chpr2) (char *, blasint *, float   *,
		     float  *, blasint *, float  *, blasint *, float  *);
void OPENBLAS_API(zhpr2) (char *, blasint *, double  *,
		     double *, blasint *, double *, blasint *, double *);
void OPENBLAS_API(xhpr2) (char *, blasint *, xdouble  *,
		     xdouble *, blasint *, xdouble *, blasint *, xdouble *);

void OPENBLAS_API(chemv) (char *, blasint *, float  *, float *, blasint *,
		     float  *, blasint *, float *, float *, blasint *);
void OPENBLAS_API(zhemv) (char *, blasint *, double  *, double *, blasint *,
		     double  *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(xhemv) (char *, blasint *, xdouble  *, xdouble *, blasint *,
		     xdouble  *, blasint *, xdouble *, xdouble *, blasint *);

void OPENBLAS_API(chpmv) (char *, blasint *, float  *, float *,
		     float  *, blasint *, float *, float *, blasint *);
void OPENBLAS_API(zhpmv) (char *, blasint *, double  *, double *,
		     double  *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(xhpmv) (char *, blasint *, xdouble  *, xdouble *,
		     xdouble  *, blasint *, xdouble *, xdouble *, blasint *);

int OPENBLAS_API(snorm)(char *, blasint *, blasint *, float  *, blasint *);
int OPENBLAS_API(dnorm)(char *, blasint *, blasint *, double *, blasint *);
int OPENBLAS_API(cnorm)(char *, blasint *, blasint *, float  *, blasint *);
int OPENBLAS_API(znorm)(char *, blasint *, blasint *, double *, blasint *);

void OPENBLAS_API(sgbmv)(char *, blasint *, blasint *, blasint *, blasint *, float  *, float  *, blasint *,
		    float  *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(dgbmv)(char *, blasint *, blasint *, blasint *, blasint *, double *, double *, blasint *,
		    double *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(qgbmv)(char *, blasint *, blasint *, blasint *, blasint *, xdouble *, xdouble *, blasint *,
		    xdouble *, blasint *, xdouble *, xdouble *, blasint *);
void OPENBLAS_API(cgbmv)(char *, blasint *, blasint *, blasint *, blasint *, float  *, float  *, blasint *,
		    float  *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(zgbmv)(char *, blasint *, blasint *, blasint *, blasint *, double *, double *, blasint *,
		    double *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(xgbmv)(char *, blasint *, blasint *, blasint *, blasint *, xdouble *, xdouble *, blasint *,
		    xdouble *, blasint *, xdouble *, xdouble *, blasint *);

void OPENBLAS_API(ssbmv)(char *, blasint *, blasint *, float  *, float  *, blasint *,
		    float  *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(dsbmv)(char *, blasint *, blasint *, double *, double *, blasint *,
		    double *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(qsbmv)(char *, blasint *, blasint *, xdouble *, xdouble *, blasint *,
		    xdouble *, blasint *, xdouble *, xdouble *, blasint *);
void OPENBLAS_API(csbmv)(char *, blasint *, blasint *, float  *, float  *, blasint *,
		    float  *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(zsbmv)(char *, blasint *, blasint *, double *, double *, blasint *,
		    double *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(xsbmv)(char *, blasint *, blasint *, xdouble *, xdouble *, blasint *,
		    xdouble *, blasint *, xdouble *, xdouble *, blasint *);

void OPENBLAS_API(chbmv)(char *, blasint *, blasint *, float  *, float  *, blasint *,
		    float  *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(zhbmv)(char *, blasint *, blasint *, double *, double *, blasint *,
		    double *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(xhbmv)(char *, blasint *, blasint *, xdouble *, xdouble *, blasint *,
		    xdouble *, blasint *, xdouble *, xdouble *, blasint *);

/* Level 3 routines */

void OPENBLAS_API(shgemm)(char *, char *, blasint *, blasint *, blasint *, float *,
	   hfloat16  *, blasint *, hfloat16 *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(bgemm)(char *, char *, blasint *, blasint *, blasint *, bfloat16 *,
	   bfloat16 *, blasint *, bfloat16 *, blasint *, bfloat16 *, bfloat16 *, blasint *);
void OPENBLAS_API(sbgemm)(char *, char *, blasint *, blasint *, blasint *, float *,
	   bfloat16 *, blasint *, bfloat16 *, blasint *, float *, float *, blasint *);
void OPENBLAS_API(sgemm)(char *, char *, blasint *, blasint *, blasint *, float *,
	   float  *, blasint *, float  *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(dgemm)(char *, char *, blasint *, blasint *, blasint *, double *,
	   double *, blasint *, double *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(qgemm)(char *, char *, blasint *, blasint *, blasint *, xdouble *,
	   xdouble *, blasint *, xdouble *, blasint *, xdouble *, xdouble *, blasint *);
void OPENBLAS_API(cgemm)(char *, char *, blasint *, blasint *, blasint *, float *,
	   float  *, blasint *, float  *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(zgemm)(char *, char *, blasint *, blasint *, blasint *, double *,
	   double *, blasint *, double *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(xgemm)(char *, char *, blasint *, blasint *, blasint *, xdouble *,
	   xdouble *, blasint *, xdouble *, blasint *, xdouble *, xdouble *, blasint *);

void OPENBLAS_API(cgemm3m)(char *, char *, blasint *, blasint *, blasint *, float *,
	   float  *, blasint *, float  *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(zgemm3m)(char *, char *, blasint *, blasint *, blasint *, double *,
	   double *, blasint *, double *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(xgemm3m)(char *, char *, blasint *, blasint *, blasint *, xdouble *,
	   xdouble *, blasint *, xdouble *, blasint *, xdouble *, xdouble *, blasint *);

void OPENBLAS_API(sgemmt)(char*, char *, char *, blasint *, blasint *, float *,
	   float  *, blasint *, float  *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(dgemmt)(char*, char *, char *, blasint *, blasint *, double *,
	   double *, blasint *, double *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(cgemmt)(char*, char *, char *, blasint *, blasint *, float *,
	   float  *, blasint *, float  *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(zgemmt)(char*, char *, char *, blasint *, blasint *, double *,
	   double *, blasint *, double *, blasint *, double *, double *, blasint *);

int OPENBLAS_API(sge2mm)(char *, char *, char *, blasint *, blasint *,
		     float *, float  *, blasint *, float  *, blasint *,
		     float *, float  *, blasint *);
int OPENBLAS_API(dge2mm)(char *, char *, char *, blasint *, blasint *,
		     double *, double  *, blasint *, double  *, blasint *,
		     double *, double  *, blasint *);
int OPENBLAS_API(cge2mm)(char *, char *, char *, blasint *, blasint *,
		     float *, float  *, blasint *, float  *, blasint *,
		     float *, float  *, blasint *);
int OPENBLAS_API(zge2mm)(char *, char *, char *, blasint *, blasint *,
		     double *, double  *, blasint *, double  *, blasint *,
		     double *, double  *, blasint *);

void OPENBLAS_API(strsm)(char *, char *, char *, char *, blasint *, blasint *,
	   float *,  float *, blasint *, float *, blasint *);
void OPENBLAS_API(dtrsm)(char *, char *, char *, char *, blasint *, blasint *,
	   double *,  double *, blasint *, double *, blasint *);
void OPENBLAS_API(qtrsm)(char *, char *, char *, char *, blasint *, blasint *,
	   xdouble *,  xdouble *, blasint *, xdouble *, blasint *);
void OPENBLAS_API(ctrsm)(char *, char *, char *, char *, blasint *, blasint *,
	   float *,  float *, blasint *, float *, blasint *);
void OPENBLAS_API(ztrsm)(char *, char *, char *, char *, blasint *, blasint *,
	   double *,  double *, blasint *, double *, blasint *);
void OPENBLAS_API(xtrsm)(char *, char *, char *, char *, blasint *, blasint *,
	   xdouble *,  xdouble *, blasint *, xdouble *, blasint *);

void OPENBLAS_API(strmm)(char *, char *, char *, char *, blasint *, blasint *,
	   float *,  float *, blasint *, float *, blasint *);
void OPENBLAS_API(dtrmm)(char *, char *, char *, char *, blasint *, blasint *,
	   double *,  double *, blasint *, double *, blasint *);
void OPENBLAS_API(qtrmm)(char *, char *, char *, char *, blasint *, blasint *,
	   xdouble *,  xdouble *, blasint *, xdouble *, blasint *);
void OPENBLAS_API(ctrmm)(char *, char *, char *, char *, blasint *, blasint *,
	   float *,  float *, blasint *, float *, blasint *);
void OPENBLAS_API(ztrmm)(char *, char *, char *, char *, blasint *, blasint *,
	   double *,  double *, blasint *, double *, blasint *);
void OPENBLAS_API(xtrmm)(char *, char *, char *, char *, blasint *, blasint *,
	   xdouble *,  xdouble *, blasint *, xdouble *, blasint *);

void OPENBLAS_API(ssymm)(char *, char *, blasint *, blasint *, float  *, float  *, blasint *,
	   float  *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(dsymm)(char *, char *, blasint *, blasint *, double *, double *, blasint *,
	   double *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(qsymm)(char *, char *, blasint *, blasint *, xdouble *, xdouble *, blasint *,
	   xdouble *, blasint *, xdouble *, xdouble *, blasint *);
void OPENBLAS_API(csymm)(char *, char *, blasint *, blasint *, float  *, float  *, blasint *,
	   float  *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(zsymm)(char *, char *, blasint *, blasint *, double *, double *, blasint *,
	   double *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(xsymm)(char *, char *, blasint *, blasint *, xdouble *, xdouble *, blasint *,
	   xdouble *, blasint *, xdouble *, xdouble *, blasint *);

void OPENBLAS_API(csymm3m)(char *, char *, blasint *, blasint *, float  *, float  *, blasint *,
	   float  *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(zsymm3m)(char *, char *, blasint *, blasint *, double *, double *, blasint *,
	   double *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(xsymm3m)(char *, char *, blasint *, blasint *, xdouble *, xdouble *, blasint *,
	   xdouble *, blasint *, xdouble *, xdouble *, blasint *);

void OPENBLAS_API(ssyrk)(char *, char *, blasint *, blasint *, float  *, float  *, blasint *,
	   float  *, float  *, blasint *);
void OPENBLAS_API(dsyrk)(char *, char *, blasint *, blasint *, double *, double *, blasint *,
	   double *, double *, blasint *);
void OPENBLAS_API(qsyrk)(char *, char *, blasint *, blasint *, xdouble *, xdouble *, blasint *,
	   xdouble *, xdouble *, blasint *);
void OPENBLAS_API(csyrk)(char *, char *, blasint *, blasint *, float  *, float  *, blasint *,
	   float  *, float  *, blasint *);
void OPENBLAS_API(zsyrk)(char *, char *, blasint *, blasint *, double *, double *, blasint *,
	   double *, double *, blasint *);
void OPENBLAS_API(xsyrk)(char *, char *, blasint *, blasint *, xdouble *, xdouble *, blasint *,
	   xdouble *, xdouble *, blasint *);

void OPENBLAS_API(ssyr2k)(char *, char *, blasint *, blasint *, float  *, float  *, blasint *,
	   float *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(dsyr2k)(char *, char *, blasint *, blasint *, double *, double *, blasint *,
	   double*, blasint *, double *, double *, blasint *);
void OPENBLAS_API(qsyr2k)(char *, char *, blasint *, blasint *, xdouble *, xdouble *, blasint *,
	   xdouble*, blasint *, xdouble *, xdouble *, blasint *);
void OPENBLAS_API(csyr2k)(char *, char *, blasint *, blasint *, float  *, float  *, blasint *,
	   float *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(zsyr2k)(char *, char *, blasint *, blasint *, double *, double *, blasint *,
	   double*, blasint *, double *, double *, blasint *);
void OPENBLAS_API(xsyr2k)(char *, char *, blasint *, blasint *, xdouble *, xdouble *, blasint *,
	   xdouble*, blasint *, xdouble *, xdouble *, blasint *);

void OPENBLAS_API(chemm)(char *, char *, blasint *, blasint *, float  *, float  *, blasint *,
	   float  *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(zhemm)(char *, char *, blasint *, blasint *, double *, double *, blasint *,
	   double *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(xhemm)(char *, char *, blasint *, blasint *, xdouble *, xdouble *, blasint *,
	   xdouble *, blasint *, xdouble *, xdouble *, blasint *);

void OPENBLAS_API(chemm3m)(char *, char *, blasint *, blasint *, float  *, float  *, blasint *,
	   float  *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(zhemm3m)(char *, char *, blasint *, blasint *, double *, double *, blasint *,
	   double *, blasint *, double *, double *, blasint *);
void OPENBLAS_API(xhemm3m)(char *, char *, blasint *, blasint *, xdouble *, xdouble *, blasint *,
	   xdouble *, blasint *, xdouble *, xdouble *, blasint *);

void OPENBLAS_API(cherk)(char *, char *, blasint *, blasint *, float  *, float  *, blasint *,
	   float  *, float  *, blasint *);
void OPENBLAS_API(zherk)(char *, char *, blasint *, blasint *, double *, double *, blasint *,
	   double *, double *, blasint *);
void OPENBLAS_API(xherk)(char *, char *, blasint *, blasint *, xdouble *, xdouble *, blasint *,
	   xdouble *, xdouble *, blasint *);

void OPENBLAS_API(cher2k)(char *, char *, blasint *, blasint *, float  *, float  *, blasint *,
	   float *, blasint *, float  *, float  *, blasint *);
void OPENBLAS_API(zher2k)(char *, char *, blasint *, blasint *, double *, double *, blasint *,
	   double*, blasint *, double *, double *, blasint *);
void OPENBLAS_API(xher2k)(char *, char *, blasint *, blasint *, xdouble *, xdouble *, blasint *,
	   xdouble*, blasint *, xdouble *, xdouble *, blasint *);

int OPENBLAS_API(cher2m)(char *, char *, char *, blasint *, blasint *, float  *, float  *, blasint *,
	   float *, blasint *, float  *, float  *, blasint *);
int OPENBLAS_API(zher2m)(char *, char *, char *, blasint *, blasint *, double *, double *, blasint *,
	   double*, blasint *, double *, double *, blasint *);
int OPENBLAS_API(xher2m)(char *, char *, char *, blasint *, blasint *, xdouble *, xdouble *, blasint *,
	   xdouble*, blasint *, xdouble *, xdouble *, blasint *);

int OPENBLAS_API(sgemt)(char *, blasint *, blasint *, float  *, float  *, blasint *,
		    float  *, blasint *);
int OPENBLAS_API(dgemt)(char *, blasint *, blasint *, double *, double *, blasint *,
		    double *, blasint *);
int OPENBLAS_API(cgemt)(char *, blasint *, blasint *, float  *, float  *, blasint *,
		    float  *, blasint *);
int OPENBLAS_API(zgemt)(char *, blasint *, blasint *, double *, double *, blasint *,
		    double *, blasint *);

int OPENBLAS_API(sgema)(char *, char *, blasint *, blasint *, float  *,
		    float  *, blasint *, float *, float  *, blasint *, float *, blasint *);
int OPENBLAS_API(dgema)(char *, char *, blasint *, blasint *, double *,
		    double *, blasint *, double*, double *, blasint *, double*, blasint *);
int OPENBLAS_API(cgema)(char *, char *, blasint *, blasint *, float  *,
		    float  *, blasint *, float *, float  *, blasint *, float *, blasint *);
int OPENBLAS_API(zgema)(char *, char *, blasint *, blasint *, double *,
		    double *, blasint *, double*, double *, blasint *, double*, blasint *);

int OPENBLAS_API(sgems)(char *, char *, blasint *, blasint *, float  *,
		    float  *, blasint *, float *, float  *, blasint *, float *, blasint *);
int OPENBLAS_API(dgems)(char *, char *, blasint *, blasint *, double *,
		    double *, blasint *, double*, double *, blasint *, double*, blasint *);
int OPENBLAS_API(cgems)(char *, char *, blasint *, blasint *, float  *,
		    float  *, blasint *, float *, float  *, blasint *, float *, blasint *);
int OPENBLAS_API(zgems)(char *, char *, blasint *, blasint *, double *,
		    double *, blasint *, double*, double *, blasint *, double*, blasint *);

int OPENBLAS_API(sgemc)(char *, char *, blasint *, blasint *, blasint *, float *,
	   float  *, blasint *, float  *, blasint *, float  *, blasint *, float  *, float  *, blasint *);
int OPENBLAS_API(dgemc)(char *, char *, blasint *, blasint *, blasint *, double *,
	   double *, blasint *, double *, blasint *, double *, blasint *, double *, double *, blasint *);
int OPENBLAS_API(qgemc)(char *, char *, blasint *, blasint *, blasint *, xdouble *,
	   xdouble *, blasint *, xdouble *, blasint *, xdouble *, blasint *,  xdouble *, xdouble *, blasint *);
int OPENBLAS_API(cgemc)(char *, char *, blasint *, blasint *, blasint *, float *,
	   float  *, blasint *, float  *, blasint *, float  *, blasint *, float  *, float  *, blasint *);
int OPENBLAS_API(zgemc)(char *, char *, blasint *, blasint *, blasint *, double *,
	   double *, blasint *, double *, blasint *, double *, blasint *, double *, double *, blasint *);
int OPENBLAS_API(xgemc)(char *, char *, blasint *, blasint *, blasint *, xdouble *,
	   xdouble *, blasint *, xdouble *, blasint *, xdouble *, blasint *, xdouble *, xdouble *, blasint *);

/* Lapack routines */

int OPENBLAS_API(sgetf2)(blasint *, blasint *, float  *, blasint *, blasint *, blasint *);
int OPENBLAS_API(dgetf2)(blasint *, blasint *, double *, blasint *, blasint *, blasint *);
int OPENBLAS_API(qgetf2)(blasint *, blasint *, xdouble *, blasint *, blasint *, blasint *);
int OPENBLAS_API(cgetf2)(blasint *, blasint *, float  *, blasint *, blasint *, blasint *);
int OPENBLAS_API(zgetf2)(blasint *, blasint *, double *, blasint *, blasint *, blasint *);
int OPENBLAS_API(xgetf2)(blasint *, blasint *, xdouble *, blasint *, blasint *, blasint *);

int OPENBLAS_API(sgetrf)(blasint *, blasint *, float  *, blasint *, blasint *, blasint *);
int OPENBLAS_API(dgetrf)(blasint *, blasint *, double *, blasint *, blasint *, blasint *);
int OPENBLAS_API(qgetrf)(blasint *, blasint *, xdouble *, blasint *, blasint *, blasint *);
int OPENBLAS_API(cgetrf)(blasint *, blasint *, float  *, blasint *, blasint *, blasint *);
int OPENBLAS_API(zgetrf)(blasint *, blasint *, double *, blasint *, blasint *, blasint *);
int OPENBLAS_API(xgetrf)(blasint *, blasint *, xdouble *, blasint *, blasint *, blasint *);

int OPENBLAS_API(slaswp)(blasint *, float  *, blasint *, blasint *, blasint *, blasint *, blasint *);
int OPENBLAS_API(dlaswp)(blasint *, double *, blasint *, blasint *, blasint *, blasint *, blasint *);
int OPENBLAS_API(qlaswp)(blasint *, xdouble *, blasint *, blasint *, blasint *, blasint *, blasint *);
int OPENBLAS_API(claswp)(blasint *, float  *, blasint *, blasint *, blasint *, blasint *, blasint *);
int OPENBLAS_API(zlaswp)(blasint *, double *, blasint *, blasint *, blasint *, blasint *, blasint *);
int OPENBLAS_API(xlaswp)(blasint *, xdouble *, blasint *, blasint *, blasint *, blasint *, blasint *);

int OPENBLAS_API(sgetrs)(char *, blasint *, blasint *, float  *, blasint *, blasint *, float  *, blasint *, blasint *);
int OPENBLAS_API(dgetrs)(char *, blasint *, blasint *, double *, blasint *, blasint *, double *, blasint *, blasint *);
int OPENBLAS_API(qgetrs)(char *, blasint *, blasint *, xdouble *, blasint *, blasint *, xdouble *, blasint *, blasint *);
int OPENBLAS_API(cgetrs)(char *, blasint *, blasint *, float  *, blasint *, blasint *, float  *, blasint *, blasint *);
int OPENBLAS_API(zgetrs)(char *, blasint *, blasint *, double *, blasint *, blasint *, double *, blasint *, blasint *);
int OPENBLAS_API(xgetrs)(char *, blasint *, blasint *, xdouble *, blasint *, blasint *, xdouble *, blasint *, blasint *);

int OPENBLAS_API(sgesv)(blasint *, blasint *, float  *, blasint *, blasint *, float *, blasint *, blasint *);
int OPENBLAS_API(dgesv)(blasint *, blasint *, double *, blasint *, blasint *, double*, blasint *, blasint *);
int OPENBLAS_API(qgesv)(blasint *, blasint *, xdouble *, blasint *, blasint *, xdouble*, blasint *, blasint *);
int OPENBLAS_API(cgesv)(blasint *, blasint *, float  *, blasint *, blasint *, float *, blasint *, blasint *);
int OPENBLAS_API(zgesv)(blasint *, blasint *, double *, blasint *, blasint *, double*, blasint *, blasint *);
int OPENBLAS_API(xgesv)(blasint *, blasint *, xdouble *, blasint *, blasint *, xdouble*, blasint *, blasint *);

int OPENBLAS_API(spotf2)(char *, blasint *, float  *, blasint *, blasint *);
int OPENBLAS_API(dpotf2)(char *, blasint *, double *, blasint *, blasint *);
int OPENBLAS_API(qpotf2)(char *, blasint *, xdouble *, blasint *, blasint *);
int OPENBLAS_API(cpotf2)(char *, blasint *, float  *, blasint *, blasint *);
int OPENBLAS_API(zpotf2)(char *, blasint *, double *, blasint *, blasint *);
int OPENBLAS_API(xpotf2)(char *, blasint *, xdouble *, blasint *, blasint *);

int OPENBLAS_API(spotrf)(char *, blasint *, float  *, blasint *, blasint *);
int OPENBLAS_API(dpotrf)(char *, blasint *, double *, blasint *, blasint *);
int OPENBLAS_API(qpotrf)(char *, blasint *, xdouble *, blasint *, blasint *);
int OPENBLAS_API(cpotrf)(char *, blasint *, float  *, blasint *, blasint *);
int OPENBLAS_API(zpotrf)(char *, blasint *, double *, blasint *, blasint *);
int OPENBLAS_API(xpotrf)(char *, blasint *, xdouble *, blasint *, blasint *);

int OPENBLAS_API(spotri)(char *, blasint *, float  *, blasint *, blasint *);
int OPENBLAS_API(dpotri)(char *, blasint *, double *, blasint *, blasint *);
int OPENBLAS_API(qpotri)(char *, blasint *, xdouble *, blasint *, blasint *);
int OPENBLAS_API(cpotri)(char *, blasint *, float  *, blasint *, blasint *);
int OPENBLAS_API(zpotri)(char *, blasint *, double *, blasint *, blasint *);
int OPENBLAS_API(xpotri)(char *, blasint *, xdouble *, blasint *, blasint *);

int OPENBLAS_API(spotrs)(char *, blasint *, blasint *, float   *, blasint *, float   *, blasint *, blasint *);
int OPENBLAS_API(dpotrs)(char *, blasint *, blasint *, double  *, blasint *, double  *, blasint *, blasint *);
int OPENBLAS_API(qpotrs)(char *, blasint *, blasint *, xdouble *, blasint *, xdouble *, blasint *, blasint *);
int OPENBLAS_API(cpotrs)(char *, blasint *, blasint *, float   *, blasint *, float   *, blasint *, blasint *);
int OPENBLAS_API(zpotrs)(char *, blasint *, blasint *, double  *, blasint *, double  *, blasint *, blasint *);
int OPENBLAS_API(xpotrs)(char *, blasint *, blasint *, xdouble *, blasint *, xdouble *, blasint *, blasint *);

int OPENBLAS_API(slauu2)(char *, blasint *, float  *, blasint *, blasint *);
int OPENBLAS_API(dlauu2)(char *, blasint *, double *, blasint *, blasint *);
int OPENBLAS_API(qlauu2)(char *, blasint *, xdouble *, blasint *, blasint *);
int OPENBLAS_API(clauu2)(char *, blasint *, float  *, blasint *, blasint *);
int OPENBLAS_API(zlauu2)(char *, blasint *, double *, blasint *, blasint *);
int OPENBLAS_API(xlauu2)(char *, blasint *, xdouble *, blasint *, blasint *);

int OPENBLAS_API(slauum)(char *, blasint *, float  *, blasint *, blasint *);
int OPENBLAS_API(dlauum)(char *, blasint *, double *, blasint *, blasint *);
int OPENBLAS_API(qlauum)(char *, blasint *, xdouble *, blasint *, blasint *);
int OPENBLAS_API(clauum)(char *, blasint *, float  *, blasint *, blasint *);
int OPENBLAS_API(zlauum)(char *, blasint *, double *, blasint *, blasint *);
int OPENBLAS_API(xlauum)(char *, blasint *, xdouble *, blasint *, blasint *);

int OPENBLAS_API(strti2)(char *, char *, blasint *, float  *, blasint *, blasint *);
int OPENBLAS_API(dtrti2)(char *, char *, blasint *, double *, blasint *, blasint *);
int OPENBLAS_API(qtrti2)(char *, char *, blasint *, xdouble *, blasint *, blasint *);
int OPENBLAS_API(ctrti2)(char *, char *, blasint *, float  *, blasint *, blasint *);
int OPENBLAS_API(ztrti2)(char *, char *, blasint *, double *, blasint *, blasint *);
int OPENBLAS_API(xtrti2)(char *, char *, blasint *, xdouble *, blasint *, blasint *);

int OPENBLAS_API(strtri)(char *, char *, blasint *, float  *, blasint *, blasint *);
int OPENBLAS_API(dtrtri)(char *, char *, blasint *, double *, blasint *, blasint *);
int OPENBLAS_API(qtrtri)(char *, char *, blasint *, xdouble *, blasint *, blasint *);
int OPENBLAS_API(ctrtri)(char *, char *, blasint *, float  *, blasint *, blasint *);
int OPENBLAS_API(ztrtri)(char *, char *, blasint *, double *, blasint *, blasint *);
int OPENBLAS_API(xtrtri)(char *, char *, blasint *, xdouble *, blasint *, blasint *);


FLOATRET  OPENBLAS_API(slamch)(char *);
double    OPENBLAS_API(dlamch)(char *);
xdouble   OPENBLAS_API(qlamch)(char *);

FLOATRET  OPENBLAS_API(slamc3)(float *, float *);
double    OPENBLAS_API(dlamc3)(double *, double *);
xdouble   OPENBLAS_API(qlamc3)(xdouble *, xdouble *);

/* BLAS extensions */

void    OPENBLAS_API(saxpby) (blasint *, float  *, float  *, blasint *, float *, float  *, blasint *);
void    OPENBLAS_API(daxpby) (blasint *, double  *, double  *, blasint *, double *, double  *, blasint *);
void    OPENBLAS_API(caxpby) (blasint *, void  *, float  *, blasint *, void *, float  *, blasint *);
void    OPENBLAS_API(zaxpby) (blasint *, void  *, double *, blasint *, void *, double  *, blasint *);

void    OPENBLAS_API(somatcopy) (char *, char *, blasint *, blasint *, float  *, float  *, blasint *, float  *, blasint *);
void    OPENBLAS_API(domatcopy) (char *, char *, blasint *, blasint *, double  *, double  *, blasint *, double  *, blasint *);
void    OPENBLAS_API(comatcopy) (char *, char *, blasint *, blasint *, float  *, float  *, blasint *, float  *, blasint *);
void    OPENBLAS_API(zomatcopy) (char *, char *, blasint *, blasint *, double  *, double  *, blasint *, double  *, blasint *);

void    OPENBLAS_API(simatcopy) (char *, char *, blasint *, blasint *, float  *, float  *, blasint *, blasint *);
void    OPENBLAS_API(dimatcopy) (char *, char *, blasint *, blasint *, double  *, double  *, blasint *, blasint *);
void    OPENBLAS_API(cimatcopy) (char *, char *, blasint *, blasint *, float  *, float  *, blasint *, blasint *);
void    OPENBLAS_API(zimatcopy) (char *, char *, blasint *, blasint *, double  *, double  *, blasint *, blasint *);

void    OPENBLAS_API(sgeadd) (blasint *, blasint *, float *, float *, blasint *, float *, float *, blasint*,char*, char*); 
void    OPENBLAS_API(dgeadd) (blasint *, blasint *, double *, double *, blasint *, double *, double *, blasint*,char *, char *); 
void    OPENBLAS_API(cgeadd) (blasint *, blasint *, float *, float *, blasint *, float *, float *, blasint*,char *, char *); 
void    OPENBLAS_API(zgeadd) (blasint *, blasint *, double *, double *, blasint *, double *, double *, blasint*,char *, char *); 


#ifdef __cplusplus
}

#endif  /* __cplusplus */

#endif
