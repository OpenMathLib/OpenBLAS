/* f2c.h  --  Standard Fortran to C header file */

/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."

	- From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */

#ifndef F2C_INCLUDE
#define F2C_INCLUDE

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <complex.h>
#ifdef complex
#undef complex
#endif
#ifdef I
#undef I
#endif

typedef int integer;
typedef unsigned int uinteger;
typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
static inline _Complex float Cf(complex *z) {return z->r + z->i*_Complex_I;}
static inline _Complex double Cd(doublecomplex *z) {return z->r + z->i*_Complex_I;}
static inline _Complex float * _pCf(complex *z) {return (_Complex float*)z;}
static inline _Complex double * _pCd(doublecomplex *z) {return (_Complex double*)z;}
#define pCf(z) (*_pCf(z))
#define pCd(z) (*_pCd(z))
typedef int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;

#define TRUE_ (1)
#define FALSE_ (0)

/* Extern is for use with -E */
#ifndef Extern
#define Extern extern
#endif

/* I/O stuff */

typedef int flag;
typedef int ftnlen;
typedef int ftnint;

/*external read, write*/
typedef struct
{	flag cierr;
	ftnint ciunit;
	flag ciend;
	char *cifmt;
	ftnint cirec;
} cilist;

/*internal read, write*/
typedef struct
{	flag icierr;
	char *iciunit;
	flag iciend;
	char *icifmt;
	ftnint icirlen;
	ftnint icirnum;
} icilist;

/*open*/
typedef struct
{	flag oerr;
	ftnint ounit;
	char *ofnm;
	ftnlen ofnmlen;
	char *osta;
	char *oacc;
	char *ofm;
	ftnint orl;
	char *oblnk;
} olist;

/*close*/
typedef struct
{	flag cerr;
	ftnint cunit;
	char *csta;
} cllist;

/*rewind, backspace, endfile*/
typedef struct
{	flag aerr;
	ftnint aunit;
} alist;

/* inquire */
typedef struct
{	flag inerr;
	ftnint inunit;
	char *infile;
	ftnlen infilen;
	ftnint	*inex;	/*parameters in standard's order*/
	ftnint	*inopen;
	ftnint	*innum;
	ftnint	*innamed;
	char	*inname;
	ftnlen	innamlen;
	char	*inacc;
	ftnlen	inacclen;
	char	*inseq;
	ftnlen	inseqlen;
	char 	*indir;
	ftnlen	indirlen;
	char	*infmt;
	ftnlen	infmtlen;
	char	*inform;
	ftnint	informlen;
	char	*inunf;
	ftnlen	inunflen;
	ftnint	*inrecl;
	ftnint	*innrec;
	char	*inblank;
	ftnlen	inblanklen;
} inlist;

#define VOID void

union Multitype {	/* for multiple entry points */
	integer1 g;
	shortint h;
	integer i;
	/* longint j; */
	real r;
	doublereal d;
	complex c;
	doublecomplex z;
	};

typedef union Multitype Multitype;

struct Vardesc {	/* for Namelist */
	char *name;
	char *addr;
	ftnlen *dims;
	int  type;
	};
typedef struct Vardesc Vardesc;

struct Namelist {
	char *name;
	Vardesc **vars;
	int nvars;
	};
typedef struct Namelist Namelist;

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (fabs(x))
#define f2cmin(a,b) ((a) <= (b) ? (a) : (b))
#define f2cmax(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (f2cmin(a,b))
#define dmax(a,b) (f2cmax(a,b))
#define bit_test(a,b)	((a) >> (b) & 1)
#define bit_clear(a,b)	((a) & ~((uinteger)1 << (b)))
#define bit_set(a,b)	((a) |  ((uinteger)1 << (b)))

#define abort_() { sig_die("Fortran abort routine called", 1); }
#define c_abs(z) (cabsf(Cf(z)))
#define c_cos(R,Z) { pCf(R)=ccos(Cf(Z)); }
#define c_div(c, a, b) {pCf(c) = Cf(a)/Cf(b);}
#define z_div(c, a, b) {pCd(c) = Cd(a)/Cd(b);}
#define c_exp(R, Z) {pCf(R) = cexpf(Cf(Z));}
#define c_log(R, Z) {pCf(R) = clogf(Cf(Z));}
#define c_sin(R, Z) {pCf(R) = csinf(Cf(Z));}
//#define c_sqrt(R, Z) {*(R) = csqrtf(Cf(Z));}
#define c_sqrt(R, Z) {pCf(R) = csqrtf(Cf(Z));}
#define d_abs(x) (fabs(*(x)))
#define d_acos(x) (acos(*(x)))
#define d_asin(x) (asin(*(x)))
#define d_atan(x) (atan(*(x)))
#define d_atn2(x, y) (atan2(*(x),*(y)))
#define d_cnjg(R, Z) { pCd(R) = conj(Cd(Z)); }
#define r_cnjg(R, Z) { pCf(R) = conj(Cf(Z)); }
#define d_cos(x) (cos(*(x)))
#define d_cosh(x) (cosh(*(x)))
#define d_dim(__a, __b) ( *(__a) > *(__b) ? *(__a) - *(__b) : 0.0 )
#define d_exp(x) (exp(*(x)))
#define d_imag(z) (cimag(Cd(z)))
#define r_imag(z) (cimag(Cf(z)))
#define d_int(__x) (*(__x)>0 ? floor(*(__x)) : -floor(- *(__x)))
#define r_int(__x) (*(__x)>0 ? floor(*(__x)) : -floor(- *(__x)))
#define d_lg10(x) ( 0.43429448190325182765 * log(*(x)) )
#define r_lg10(x) ( 0.43429448190325182765 * log(*(x)) )
#define d_log(x) (log(*(x)))
#define d_mod(x, y) (fmod(*(x), *(y)))
#define u_nint(__x) ((__x)>=0 ? floor((__x) + .5) : -floor(.5 - (__x)))
#define d_nint(x) u_nint(*(x))
#define u_sign(__a,__b) ((__b) >= 0 ? ((__a) >= 0 ? (__a) : -(__a)) : -((__a) >= 0 ? (__a) : -(__a)))
#define d_sign(a,b) u_sign(*(a),*(b))
#define r_sign(a,b) u_sign(*(a),*(b))
#define d_sin(x) (sin(*(x)))
#define d_sinh(x) (sinh(*(x)))
#define d_sqrt(x) (sqrt(*(x)))
#define d_tan(x) (tan(*(x)))
#define d_tanh(x) (tanh(*(x)))
#define i_abs(x) abs(*(x))
#define i_dnnt(x) ((integer)u_nint(*(x)))
#define i_len(s, n) (n)
#define i_nint(x) ((integer)u_nint(*(x)))
#define i_sign(a,b) ((integer)u_sign((integer)*(a),(integer)*(b)))
#define pow_dd(ap, bp) ( pow(*(ap), *(bp)))
#define pow_si(B,E) spow_ui(*(B),*(E))
#define pow_ri(B,E) spow_ui(*(B),*(E))
#define pow_di(B,E) dpow_ui(*(B),*(E))
#define pow_zi(p, a, b) {pCd(p) = zpow_ui(Cd(a), *(b));}
#define pow_ci(p, a, b) {pCf(p) = cpow_ui(Cf(a), *(b));}
#define pow_zz(R,A,B) {pCd(R) = cpow(Cd(A),*(B));}
#define s_cat(lpp, rpp, rnp, np, llp) { 	ftnlen i, nc, ll; char *f__rp, *lp; 	ll = (llp); lp = (lpp); 	for(i=0; i < (int)*(np); ++i) {         	nc = ll; 	        if((rnp)[i] < nc) nc = (rnp)[i]; 	        ll -= nc;         	f__rp = (rpp)[i]; 	        while(--nc >= 0) *lp++ = *(f__rp)++;         } 	while(--ll >= 0) *lp++ = ' '; }
#define s_cmp(a,b,c,d) ((integer)strncmp((a),(b),f2cmin((c),(d))))
#define s_copy(A,B,C,D) { int __i,__m; for (__i=0, __m=f2cmin((C),(D)); __i<__m && (B)[__i] != 0; ++__i) (A)[__i] = (B)[__i]; }
#define sig_die(s, kill) { exit(1); }
#define s_stop(s, n) {exit(0);}
static char junk[] = "\n@(#)LIBF77 VERSION 19990503\n";
#define z_abs(z) (cabs(Cd(z)))
#define z_exp(R, Z) {pCd(R) = cexp(Cd(Z));}
#define z_sqrt(R, Z) {pCd(R) = csqrt(Cd(Z));}
#define myexit_() break;
#define mycycle() continue;
#define myceiling(w) {ceil(w)}
#define myhuge(w) {HUGE_VAL}
//#define mymaxloc_(w,s,e,n) {if (sizeof(*(w)) == sizeof(double)) dmaxloc_((w),*(s),*(e),n); else dmaxloc_((w),*(s),*(e),n);}
#define mymaxloc(w,s,e,n) {dmaxloc_(w,*(s),*(e),n)}

/* procedure parameter types for -A and -C++ */

#define F2C_proc_par_types 1
#ifdef __cplusplus
typedef logical (*L_fp)(...);
#else
typedef logical (*L_fp)();
#endif

static float spow_ui(float x, integer n) {
	float pow=1.0; unsigned long int u;
	if(n != 0) {
		if(n < 0) n = -n, x = 1/x;
		for(u = n; ; ) {
			if(u & 01) pow *= x;
			if(u >>= 1) x *= x;
			else break;
		}
	}
	return pow;
}
static double dpow_ui(double x, integer n) {
	double pow=1.0; unsigned long int u;
	if(n != 0) {
		if(n < 0) n = -n, x = 1/x;
		for(u = n; ; ) {
			if(u & 01) pow *= x;
			if(u >>= 1) x *= x;
			else break;
		}
	}
	return pow;
}
static _Complex float cpow_ui(_Complex float x, integer n) {
	_Complex float pow=1.0; unsigned long int u;
	if(n != 0) {
		if(n < 0) n = -n, x = 1/x;
		for(u = n; ; ) {
			if(u & 01) pow *= x;
			if(u >>= 1) x *= x;
			else break;
		}
	}
	return pow;
}
static _Complex double zpow_ui(_Complex double x, integer n) {
	_Complex double pow=1.0; unsigned long int u;
	if(n != 0) {
		if(n < 0) n = -n, x = 1/x;
		for(u = n; ; ) {
			if(u & 01) pow *= x;
			if(u >>= 1) x *= x;
			else break;
		}
	}
	return pow;
}
static integer pow_ii(integer x, integer n) {
	integer pow; unsigned long int u;
	if (n <= 0) {
		if (n == 0 || x == 1) pow = 1;
		else if (x != -1) pow = x == 0 ? 1/x : 0;
		else n = -n;
	}
	if ((n > 0) || !(n == 0 || x == 1 || x != -1)) {
		u = n;
		for(pow = 1; ; ) {
			if(u & 01) pow *= x;
			if(u >>= 1) x *= x;
			else break;
		}
	}
	return pow;
}
static integer dmaxloc_(double *w, integer s, integer e, integer *n)
{
	double m; integer i, mi;
	for(m=w[s-1], mi=s, i=s+1; i<=e; i++)
		if (w[i-1]>m) mi=i ,m=w[i-1];
	return mi-s+1;
}
static integer smaxloc_(float *w, integer s, integer e, integer *n)
{
	float m; integer i, mi;
	for(m=w[s-1], mi=s, i=s+1; i<=e; i++)
		if (w[i-1]>m) mi=i ,m=w[i-1];
	return mi-s+1;
}
static inline void cdotc_(complex *z, integer *n_, complex *x, integer *incx_, complex *y, integer *incy_) {
	integer n = *n_, incx = *incx_, incy = *incy_, i;
	_Complex float zdotc = 0.0;
	if (incx == 1 && incy == 1) {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += conjf(Cf(&x[i])) * Cf(&y[i]);
		}
	} else {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += conjf(Cf(&x[i*incx])) * Cf(&y[i*incy]);
		}
	}
	pCf(z) = zdotc;
}
static inline void zdotc_(doublecomplex *z, integer *n_, doublecomplex *x, integer *incx_, doublecomplex *y, integer *incy_) {
	integer n = *n_, incx = *incx_, incy = *incy_, i;
	_Complex double zdotc = 0.0;
	if (incx == 1 && incy == 1) {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += conj(Cd(&x[i])) * Cd(&y[i]);
		}
	} else {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += conj(Cd(&x[i*incx])) * Cd(&y[i*incy]);
		}
	}
	pCd(z) = zdotc;
}	
static inline void cdotu_(complex *z, integer *n_, complex *x, integer *incx_, complex *y, integer *incy_) {
	integer n = *n_, incx = *incx_, incy = *incy_, i;
	_Complex float zdotc = 0.0;
	if (incx == 1 && incy == 1) {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += Cf(&x[i]) * Cf(&y[i]);
		}
	} else {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += Cf(&x[i*incx]) * Cf(&y[i*incy]);
		}
	}
	pCf(z) = zdotc;
}
static inline void zdotu_(doublecomplex *z, integer *n_, doublecomplex *x, integer *incx_, doublecomplex *y, integer *incy_) {
	integer n = *n_, incx = *incx_, incy = *incy_, i;
	_Complex double zdotc = 0.0;
	if (incx == 1 && incy == 1) {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += Cd(&x[i]) * Cd(&y[i]);
		}
	} else {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += Cd(&x[i*incx]) * Cd(&y[i*incy]);
		}
	}
	pCd(z) = zdotc;
}
#endif
/*  -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/



/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SLANTR returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the ele
ment of largest absolute value of a trapezoidal or triangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLANTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slantr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slantr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slantr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       REAL             FUNCTION SLANTR( NORM, UPLO, DIAG, M, N, A, LDA, */
/*                        WORK ) */

/*       CHARACTER          DIAG, NORM, UPLO */
/*       INTEGER            LDA, M, N */
/*       REAL               A( LDA, * ), WORK( * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLANTR  returns the value of the one norm,  or the Frobenius norm, or */
/* > the  infinity norm,  or the  element of  largest absolute value  of a */
/* > trapezoidal or triangular matrix A. */
/* > \endverbatim */
/* > */
/* > \return SLANTR */
/* > \verbatim */
/* > */
/* >    SLANTR = ( f2cmax(abs(A(i,j))), NORM = 'M' or 'm' */
/* >             ( */
/* >             ( norm1(A),         NORM = '1', 'O' or 'o' */
/* >             ( */
/* >             ( normI(A),         NORM = 'I' or 'i' */
/* >             ( */
/* >             ( normF(A),         NORM = 'F', 'f', 'E' or 'e' */
/* > */
/* > where  norm1  denotes the  one norm of a matrix (maximum column sum), */
/* > normI  denotes the  infinity norm  of a matrix  (maximum row sum) and */
/* > normF  denotes the  Frobenius norm of a matrix (square root of sum of */
/* > squares).  Note that  f2cmax(abs(A(i,j)))  is not a consistent matrix norm. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] NORM */
/* > \verbatim */
/* >          NORM is CHARACTER*1 */
/* >          Specifies the value to be returned in SLANTR as described */
/* >          above. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the matrix A is upper or lower trapezoidal. */
/* >          = 'U':  Upper trapezoidal */
/* >          = 'L':  Lower trapezoidal */
/* >          Note that A is triangular instead of trapezoidal if M = N. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >          Specifies whether or not the matrix A has unit diagonal. */
/* >          = 'N':  Non-unit diagonal */
/* >          = 'U':  Unit diagonal */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0, and if */
/* >          UPLO = 'U', M <= N.  When M = 0, SLANTR is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0, and if */
/* >          UPLO = 'L', N <= M.  When N = 0, SLANTR is set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          The trapezoidal matrix A (A is triangular if M = N). */
/* >          If UPLO = 'U', the leading m by n upper trapezoidal part of */
/* >          the array A contains the upper trapezoidal matrix, and the */
/* >          strictly lower triangular part of A is not referenced. */
/* >          If UPLO = 'L', the leading m by n lower trapezoidal part of */
/* >          the array A contains the lower trapezoidal matrix, and the */
/* >          strictly upper triangular part of A is not referenced.  Note */
/* >          that when DIAG = 'U', the diagonal elements of A are not */
/* >          referenced and are assumed to be one. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= f2cmax(M,1). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (MAX(1,LWORK)), */
/* >          where LWORK >= M when NORM = 'I'; otherwise, WORK is not */
/* >          referenced. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
real slantr_(char *norm, char *uplo, char *diag, integer *m, integer *n, real 
	*a, integer *lda, real *work)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    real ret_val, r__1;

    /* Local variables */
    extern /* Subroutine */ int scombssq_(real *, real *);
    integer i__, j;
    logical udiag;
    extern logical lsame_(char *, char *);
    real value;
    extern logical sisnan_(real *);
    real colssq[2];
    extern /* Subroutine */ int slassq_(integer *, real *, integer *, real *, 
	    real *);
    real sum, ssq[2];


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */


/* ===================================================================== */


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --work;

    /* Function Body */
    if (f2cmin(*m,*n) == 0) {
	value = 0.f;
    } else if (lsame_(norm, "M")) {

/*        Find f2cmax(abs(A(i,j))). */

	if (lsame_(diag, "U")) {
	    value = 1.f;
	    if (lsame_(uplo, "U")) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
		    i__3 = *m, i__4 = j - 1;
		    i__2 = f2cmin(i__3,i__4);
		    for (i__ = 1; i__ <= i__2; ++i__) {
			sum = (r__1 = a[i__ + j * a_dim1], abs(r__1));
			if (value < sum || sisnan_(&sum)) {
			    value = sum;
			}
/* L10: */
		    }
/* L20: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			sum = (r__1 = a[i__ + j * a_dim1], abs(r__1));
			if (value < sum || sisnan_(&sum)) {
			    value = sum;
			}
/* L30: */
		    }
/* L40: */
		}
	    }
	} else {
	    value = 0.f;
	    if (lsame_(uplo, "U")) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = f2cmin(*m,j);
		    for (i__ = 1; i__ <= i__2; ++i__) {
			sum = (r__1 = a[i__ + j * a_dim1], abs(r__1));
			if (value < sum || sisnan_(&sum)) {
			    value = sum;
			}
/* L50: */
		    }
/* L60: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = j; i__ <= i__2; ++i__) {
			sum = (r__1 = a[i__ + j * a_dim1], abs(r__1));
			if (value < sum || sisnan_(&sum)) {
			    value = sum;
			}
/* L70: */
		    }
/* L80: */
		}
	    }
	}
    } else if (lsame_(norm, "O") || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

	value = 0.f;
	udiag = lsame_(diag, "U");
	if (lsame_(uplo, "U")) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (udiag && j <= *m) {
		    sum = 1.f;
		    i__2 = j - 1;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			sum += (r__1 = a[i__ + j * a_dim1], abs(r__1));
/* L90: */
		    }
		} else {
		    sum = 0.f;
		    i__2 = f2cmin(*m,j);
		    for (i__ = 1; i__ <= i__2; ++i__) {
			sum += (r__1 = a[i__ + j * a_dim1], abs(r__1));
/* L100: */
		    }
		}
		if (value < sum || sisnan_(&sum)) {
		    value = sum;
		}
/* L110: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (udiag) {
		    sum = 1.f;
		    i__2 = *m;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			sum += (r__1 = a[i__ + j * a_dim1], abs(r__1));
/* L120: */
		    }
		} else {
		    sum = 0.f;
		    i__2 = *m;
		    for (i__ = j; i__ <= i__2; ++i__) {
			sum += (r__1 = a[i__ + j * a_dim1], abs(r__1));
/* L130: */
		    }
		}
		if (value < sum || sisnan_(&sum)) {
		    value = sum;
		}
/* L140: */
	    }
	}
    } else if (lsame_(norm, "I")) {

/*        Find normI(A). */

	if (lsame_(uplo, "U")) {
	    if (lsame_(diag, "U")) {
		i__1 = *m;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    work[i__] = 1.f;
/* L150: */
		}
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
		    i__3 = *m, i__4 = j - 1;
		    i__2 = f2cmin(i__3,i__4);
		    for (i__ = 1; i__ <= i__2; ++i__) {
			work[i__] += (r__1 = a[i__ + j * a_dim1], abs(r__1));
/* L160: */
		    }
/* L170: */
		}
	    } else {
		i__1 = *m;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    work[i__] = 0.f;
/* L180: */
		}
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = f2cmin(*m,j);
		    for (i__ = 1; i__ <= i__2; ++i__) {
			work[i__] += (r__1 = a[i__ + j * a_dim1], abs(r__1));
/* L190: */
		    }
/* L200: */
		}
	    }
	} else {
	    if (lsame_(diag, "U")) {
		i__1 = f2cmin(*m,*n);
		for (i__ = 1; i__ <= i__1; ++i__) {
		    work[i__] = 1.f;
/* L210: */
		}
		i__1 = *m;
		for (i__ = *n + 1; i__ <= i__1; ++i__) {
		    work[i__] = 0.f;
/* L220: */
		}
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			work[i__] += (r__1 = a[i__ + j * a_dim1], abs(r__1));
/* L230: */
		    }
/* L240: */
		}
	    } else {
		i__1 = *m;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    work[i__] = 0.f;
/* L250: */
		}
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = j; i__ <= i__2; ++i__) {
			work[i__] += (r__1 = a[i__ + j * a_dim1], abs(r__1));
/* L260: */
		    }
/* L270: */
		}
	    }
	}
	value = 0.f;
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sum = work[i__];
	    if (value < sum || sisnan_(&sum)) {
		value = sum;
	    }
/* L280: */
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */
/*        SSQ(1) is scale */
/*        SSQ(2) is sum-of-squares */
/*        For better accuracy, sum each column separately. */

	if (lsame_(uplo, "U")) {
	    if (lsame_(diag, "U")) {
		ssq[0] = 1.f;
		ssq[1] = (real) f2cmin(*m,*n);
		i__1 = *n;
		for (j = 2; j <= i__1; ++j) {
		    colssq[0] = 0.f;
		    colssq[1] = 1.f;
/* Computing MIN */
		    i__3 = *m, i__4 = j - 1;
		    i__2 = f2cmin(i__3,i__4);
		    slassq_(&i__2, &a[j * a_dim1 + 1], &c__1, colssq, &colssq[
			    1]);
		    scombssq_(ssq, colssq);
/* L290: */
		}
	    } else {
		ssq[0] = 0.f;
		ssq[1] = 1.f;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    colssq[0] = 0.f;
		    colssq[1] = 1.f;
		    i__2 = f2cmin(*m,j);
		    slassq_(&i__2, &a[j * a_dim1 + 1], &c__1, colssq, &colssq[
			    1]);
		    scombssq_(ssq, colssq);
/* L300: */
		}
	    }
	} else {
	    if (lsame_(diag, "U")) {
		ssq[0] = 1.f;
		ssq[1] = (real) f2cmin(*m,*n);
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    colssq[0] = 0.f;
		    colssq[1] = 1.f;
		    i__2 = *m - j;
/* Computing MIN */
		    i__3 = *m, i__4 = j + 1;
		    slassq_(&i__2, &a[f2cmin(i__3,i__4) + j * a_dim1], &c__1, 
			    colssq, &colssq[1]);
		    scombssq_(ssq, colssq);
/* L310: */
		}
	    } else {
		ssq[0] = 0.f;
		ssq[1] = 1.f;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    colssq[0] = 0.f;
		    colssq[1] = 1.f;
		    i__2 = *m - j + 1;
		    slassq_(&i__2, &a[j + j * a_dim1], &c__1, colssq, &colssq[
			    1]);
		    scombssq_(ssq, colssq);
/* L320: */
		}
	    }
	}
	value = ssq[0] * sqrt(ssq[1]);
    }

    ret_val = value;
    return ret_val;

/*     End of SLANTR */

} /* slantr_ */

