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

/* > \brief <b> DGBSVX computes the solution to system of linear equations A * X = B for GB matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGBSVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbsvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbsvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbsvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGBSVX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, */
/*                          LDAFB, IPIV, EQUED, R, C, B, LDB, X, LDX, */
/*                          RCOND, FERR, BERR, WORK, IWORK, INFO ) */

/*       CHARACTER          EQUED, FACT, TRANS */
/*       INTEGER            INFO, KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS */
/*       DOUBLE PRECISION   RCOND */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       DOUBLE PRECISION   AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), */
/*      $                   BERR( * ), C( * ), FERR( * ), R( * ), */
/*      $                   WORK( * ), X( LDX, * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGBSVX uses the LU factorization to compute the solution to a real */
/* > system of linear equations A * X = B, A**T * X = B, or A**H * X = B, */
/* > where A is a band matrix of order N with KL subdiagonals and KU */
/* > superdiagonals, and X and B are N-by-NRHS matrices. */
/* > */
/* > Error bounds on the solution and a condition estimate are also */
/* > provided. */
/* > \endverbatim */

/* > \par Description: */
/*  ================= */
/* > */
/* > \verbatim */
/* > */
/* > The following steps are performed by this subroutine: */
/* > */
/* > 1. If FACT = 'E', real scaling factors are computed to equilibrate */
/* >    the system: */
/* >       TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B */
/* >       TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B */
/* >       TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B */
/* >    Whether or not the system will be equilibrated depends on the */
/* >    scaling of the matrix A, but if equilibration is used, A is */
/* >    overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N') */
/* >    or diag(C)*B (if TRANS = 'T' or 'C'). */
/* > */
/* > 2. If FACT = 'N' or 'E', the LU decomposition is used to factor the */
/* >    matrix A (after equilibration if FACT = 'E') as */
/* >       A = L * U, */
/* >    where L is a product of permutation and unit lower triangular */
/* >    matrices with KL subdiagonals, and U is upper triangular with */
/* >    KL+KU superdiagonals. */
/* > */
/* > 3. If some U(i,i)=0, so that U is exactly singular, then the routine */
/* >    returns with INFO = i. Otherwise, the factored form of A is used */
/* >    to estimate the condition number of the matrix A.  If the */
/* >    reciprocal of the condition number is less than machine precision, */
/* >    INFO = N+1 is returned as a warning, but the routine still goes on */
/* >    to solve for X and compute error bounds as described below. */
/* > */
/* > 4. The system of equations is solved for X using the factored form */
/* >    of A. */
/* > */
/* > 5. Iterative refinement is applied to improve the computed solution */
/* >    matrix and calculate error bounds and backward error estimates */
/* >    for it. */
/* > */
/* > 6. If equilibration was used, the matrix X is premultiplied by */
/* >    diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so */
/* >    that it solves the original system before equilibration. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] FACT */
/* > \verbatim */
/* >          FACT is CHARACTER*1 */
/* >          Specifies whether or not the factored form of the matrix A is */
/* >          supplied on entry, and if not, whether the matrix A should be */
/* >          equilibrated before it is factored. */
/* >          = 'F':  On entry, AFB and IPIV contain the factored form of */
/* >                  A.  If EQUED is not 'N', the matrix A has been */
/* >                  equilibrated with scaling factors given by R and C. */
/* >                  AB, AFB, and IPIV are not modified. */
/* >          = 'N':  The matrix A will be copied to AFB and factored. */
/* >          = 'E':  The matrix A will be equilibrated if necessary, then */
/* >                  copied to AFB and factored. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          Specifies the form of the system of equations. */
/* >          = 'N':  A * X = B     (No transpose) */
/* >          = 'T':  A**T * X = B  (Transpose) */
/* >          = 'C':  A**H * X = B  (Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of linear equations, i.e., the order of the */
/* >          matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KL */
/* > \verbatim */
/* >          KL is INTEGER */
/* >          The number of subdiagonals within the band of A.  KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* >          KU is INTEGER */
/* >          The number of superdiagonals within the band of A.  KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrices B and X.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is DOUBLE PRECISION array, dimension (LDAB,N) */
/* >          On entry, the matrix A in band storage, in rows 1 to KL+KU+1. */
/* >          The j-th column of A is stored in the j-th column of the */
/* >          array AB as follows: */
/* >          AB(KU+1+i-j,j) = A(i,j) for f2cmax(1,j-KU)<=i<=f2cmin(N,j+kl) */
/* > */
/* >          If FACT = 'F' and EQUED is not 'N', then A must have been */
/* >          equilibrated by the scaling factors in R and/or C.  AB is not */
/* >          modified if FACT = 'F' or 'N', or if FACT = 'E' and */
/* >          EQUED = 'N' on exit. */
/* > */
/* >          On exit, if EQUED .ne. 'N', A is scaled as follows: */
/* >          EQUED = 'R':  A := diag(R) * A */
/* >          EQUED = 'C':  A := A * diag(C) */
/* >          EQUED = 'B':  A := diag(R) * A * diag(C). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AFB */
/* > \verbatim */
/* >          AFB is DOUBLE PRECISION array, dimension (LDAFB,N) */
/* >          If FACT = 'F', then AFB is an input argument and on entry */
/* >          contains details of the LU factorization of the band matrix */
/* >          A, as computed by DGBTRF.  U is stored as an upper triangular */
/* >          band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, */
/* >          and the multipliers used during the factorization are stored */
/* >          in rows KL+KU+2 to 2*KL+KU+1.  If EQUED .ne. 'N', then AFB is */
/* >          the factored form of the equilibrated matrix A. */
/* > */
/* >          If FACT = 'N', then AFB is an output argument and on exit */
/* >          returns details of the LU factorization of A. */
/* > */
/* >          If FACT = 'E', then AFB is an output argument and on exit */
/* >          returns details of the LU factorization of the equilibrated */
/* >          matrix A (see the description of AB for the form of the */
/* >          equilibrated matrix). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAFB */
/* > \verbatim */
/* >          LDAFB is INTEGER */
/* >          The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          If FACT = 'F', then IPIV is an input argument and on entry */
/* >          contains the pivot indices from the factorization A = L*U */
/* >          as computed by DGBTRF; row i of the matrix was interchanged */
/* >          with row IPIV(i). */
/* > */
/* >          If FACT = 'N', then IPIV is an output argument and on exit */
/* >          contains the pivot indices from the factorization A = L*U */
/* >          of the original matrix A. */
/* > */
/* >          If FACT = 'E', then IPIV is an output argument and on exit */
/* >          contains the pivot indices from the factorization A = L*U */
/* >          of the equilibrated matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] EQUED */
/* > \verbatim */
/* >          EQUED is CHARACTER*1 */
/* >          Specifies the form of equilibration that was done. */
/* >          = 'N':  No equilibration (always true if FACT = 'N'). */
/* >          = 'R':  Row equilibration, i.e., A has been premultiplied by */
/* >                  diag(R). */
/* >          = 'C':  Column equilibration, i.e., A has been postmultiplied */
/* >                  by diag(C). */
/* >          = 'B':  Both row and column equilibration, i.e., A has been */
/* >                  replaced by diag(R) * A * diag(C). */
/* >          EQUED is an input argument if FACT = 'F'; otherwise, it is an */
/* >          output argument. */
/* > \endverbatim */
/* > */
/* > \param[in,out] R */
/* > \verbatim */
/* >          R is DOUBLE PRECISION array, dimension (N) */
/* >          The row scale factors for A.  If EQUED = 'R' or 'B', A is */
/* >          multiplied on the left by diag(R); if EQUED = 'N' or 'C', R */
/* >          is not accessed.  R is an input argument if FACT = 'F'; */
/* >          otherwise, R is an output argument.  If FACT = 'F' and */
/* >          EQUED = 'R' or 'B', each element of R must be positive. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (N) */
/* >          The column scale factors for A.  If EQUED = 'C' or 'B', A is */
/* >          multiplied on the right by diag(C); if EQUED = 'N' or 'R', C */
/* >          is not accessed.  C is an input argument if FACT = 'F'; */
/* >          otherwise, C is an output argument.  If FACT = 'F' and */
/* >          EQUED = 'C' or 'B', each element of C must be positive. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
/* >          On entry, the right hand side matrix B. */
/* >          On exit, */
/* >          if EQUED = 'N', B is not modified; */
/* >          if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by */
/* >          diag(R)*B; */
/* >          if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is */
/* >          overwritten by diag(C)*B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= f2cmax(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension (LDX,NRHS) */
/* >          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X */
/* >          to the original system of equations.  Note that A and B are */
/* >          modified on exit if EQUED .ne. 'N', and the solution to the */
/* >          equilibrated system is inv(diag(C))*X if TRANS = 'N' and */
/* >          EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C' */
/* >          and EQUED = 'R' or 'B'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* >          LDX is INTEGER */
/* >          The leading dimension of the array X.  LDX >= f2cmax(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* >          RCOND is DOUBLE PRECISION */
/* >          The estimate of the reciprocal condition number of the matrix */
/* >          A after equilibration (if done).  If RCOND is less than the */
/* >          machine precision (in particular, if RCOND = 0), the matrix */
/* >          is singular to working precision.  This condition is */
/* >          indicated by a return code of INFO > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] FERR */
/* > \verbatim */
/* >          FERR is DOUBLE PRECISION array, dimension (NRHS) */
/* >          The estimated forward error bound for each solution vector */
/* >          X(j) (the j-th column of the solution matrix X). */
/* >          If XTRUE is the true solution corresponding to X(j), FERR(j) */
/* >          is an estimated upper bound for the magnitude of the largest */
/* >          element in (X(j) - XTRUE) divided by the magnitude of the */
/* >          largest element in X(j).  The estimate is as reliable as */
/* >          the estimate for RCOND, and is almost always a slight */
/* >          overestimate of the true error. */
/* > \endverbatim */
/* > */
/* > \param[out] BERR */
/* > \verbatim */
/* >          BERR is DOUBLE PRECISION array, dimension (NRHS) */
/* >          The componentwise relative backward error of each solution */
/* >          vector X(j) (i.e., the smallest relative change in */
/* >          any element of A or B that makes X(j) an exact solution). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (3*N) */
/* >          On exit, WORK(1) contains the reciprocal pivot growth */
/* >          factor norm(A)/norm(U). The "f2cmax absolute element" norm is */
/* >          used. If WORK(1) is much less than 1, then the stability */
/* >          of the LU factorization of the (equilibrated) matrix A */
/* >          could be poor. This also means that the solution X, condition */
/* >          estimator RCOND, and forward error bound FERR could be */
/* >          unreliable. If factorization fails with 0<INFO<=N, then */
/* >          WORK(1) contains the reciprocal pivot growth factor for the */
/* >          leading INFO columns of A. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, and i is */
/* >                <= N:  U(i,i) is exactly zero.  The factorization */
/* >                       has been completed, but the factor U is exactly */
/* >                       singular, so the solution and error bounds */
/* >                       could not be computed. RCOND = 0 is returned. */
/* >                = N+1: U is nonsingular, but RCOND is less than machine */
/* >                       precision, meaning that the matrix is singular */
/* >                       to working precision.  Nevertheless, the */
/* >                       solution and error bounds are computed because */
/* >                       there are a number of situations where the */
/* >                       computed solution can be more accurate than the */
/* >                       value of RCOND would suggest. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date April 2012 */

/* > \ingroup doubleGBsolve */

/*  ===================================================================== */
/* Subroutine */ int dgbsvx_(char *fact, char *trans, integer *n, integer *kl,
	 integer *ku, integer *nrhs, doublereal *ab, integer *ldab, 
	doublereal *afb, integer *ldafb, integer *ipiv, char *equed, 
	doublereal *r__, doublereal *c__, doublereal *b, integer *ldb, 
	doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr, 
	doublereal *berr, doublereal *work, integer *iwork, integer *info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset, 
	    x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    doublereal amax;
    char norm[1];
    integer i__, j;
    extern logical lsame_(char *, char *);
    doublereal rcmin, rcmax, anorm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    logical equil;
    integer j1, j2;
    extern doublereal dlamch_(char *), dlangb_(char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *);
    extern /* Subroutine */ int dlaqgb_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, char *), 
	    dgbcon_(char *, integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *);
    doublereal colcnd;
    extern doublereal dlantb_(char *, char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *);
    extern /* Subroutine */ int dgbequ_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *), dgbrfs_(
	    char *, integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *), dgbtrf_(integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, integer 
	    *, integer *);
    logical nofact;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *), 
	    xerbla_(char *, integer *, ftnlen);
    doublereal bignum;
    extern /* Subroutine */ int dgbtrs_(char *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *);
    integer infequ;
    logical colequ;
    doublereal rowcnd;
    logical notran;
    doublereal smlnum;
    logical rowequ;
    doublereal rpvgrw;


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     April 2012 */


/*  ===================================================================== */


    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1 * 1;
    ab -= ab_offset;
    afb_dim1 = *ldafb;
    afb_offset = 1 + afb_dim1 * 1;
    afb -= afb_offset;
    --ipiv;
    --r__;
    --c__;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1 * 1;
    x -= x_offset;
    --ferr;
    --berr;
    --work;
    --iwork;

    /* Function Body */
    *info = 0;
    nofact = lsame_(fact, "N");
    equil = lsame_(fact, "E");
    notran = lsame_(trans, "N");
    if (nofact || equil) {
	*(unsigned char *)equed = 'N';
	rowequ = FALSE_;
	colequ = FALSE_;
    } else {
	rowequ = lsame_(equed, "R") || lsame_(equed, 
		"B");
	colequ = lsame_(equed, "C") || lsame_(equed, 
		"B");
	smlnum = dlamch_("Safe minimum");
	bignum = 1. / smlnum;
    }

/*     Test the input parameters. */

    if (! nofact && ! equil && ! lsame_(fact, "F")) {
	*info = -1;
    } else if (! notran && ! lsame_(trans, "T") && ! 
	    lsame_(trans, "C")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*kl < 0) {
	*info = -4;
    } else if (*ku < 0) {
	*info = -5;
    } else if (*nrhs < 0) {
	*info = -6;
    } else if (*ldab < *kl + *ku + 1) {
	*info = -8;
    } else if (*ldafb < (*kl << 1) + *ku + 1) {
	*info = -10;
    } else if (lsame_(fact, "F") && ! (rowequ || colequ 
	    || lsame_(equed, "N"))) {
	*info = -12;
    } else {
	if (rowequ) {
	    rcmin = bignum;
	    rcmax = 0.;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
		d__1 = rcmin, d__2 = r__[j];
		rcmin = f2cmin(d__1,d__2);
/* Computing MAX */
		d__1 = rcmax, d__2 = r__[j];
		rcmax = f2cmax(d__1,d__2);
/* L10: */
	    }
	    if (rcmin <= 0.) {
		*info = -13;
	    } else if (*n > 0) {
		rowcnd = f2cmax(rcmin,smlnum) / f2cmin(rcmax,bignum);
	    } else {
		rowcnd = 1.;
	    }
	}
	if (colequ && *info == 0) {
	    rcmin = bignum;
	    rcmax = 0.;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
		d__1 = rcmin, d__2 = c__[j];
		rcmin = f2cmin(d__1,d__2);
/* Computing MAX */
		d__1 = rcmax, d__2 = c__[j];
		rcmax = f2cmax(d__1,d__2);
/* L20: */
	    }
	    if (rcmin <= 0.) {
		*info = -14;
	    } else if (*n > 0) {
		colcnd = f2cmax(rcmin,smlnum) / f2cmin(rcmax,bignum);
	    } else {
		colcnd = 1.;
	    }
	}
	if (*info == 0) {
	    if (*ldb < f2cmax(1,*n)) {
		*info = -16;
	    } else if (*ldx < f2cmax(1,*n)) {
		*info = -18;
	    }
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGBSVX", &i__1, (ftnlen)6);
	return 0;
    }

    if (equil) {

/*        Compute row and column scalings to equilibrate the matrix A. */

	dgbequ_(n, n, kl, ku, &ab[ab_offset], ldab, &r__[1], &c__[1], &rowcnd,
		 &colcnd, &amax, &infequ);
	if (infequ == 0) {

/*           Equilibrate the matrix. */

	    dlaqgb_(n, n, kl, ku, &ab[ab_offset], ldab, &r__[1], &c__[1], &
		    rowcnd, &colcnd, &amax, equed);
	    rowequ = lsame_(equed, "R") || lsame_(equed,
		     "B");
	    colequ = lsame_(equed, "C") || lsame_(equed,
		     "B");
	}
    }

/*     Scale the right hand side. */

    if (notran) {
	if (rowequ) {
	    i__1 = *nrhs;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    b[i__ + j * b_dim1] = r__[i__] * b[i__ + j * b_dim1];
/* L30: */
		}
/* L40: */
	    }
	}
    } else if (colequ) {
	i__1 = *nrhs;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		b[i__ + j * b_dim1] = c__[i__] * b[i__ + j * b_dim1];
/* L50: */
	    }
/* L60: */
	}
    }

    if (nofact || equil) {

/*        Compute the LU factorization of the band matrix A. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	    i__2 = j - *ku;
	    j1 = f2cmax(i__2,1);
/* Computing MIN */
	    i__2 = j + *kl;
	    j2 = f2cmin(i__2,*n);
	    i__2 = j2 - j1 + 1;
	    dcopy_(&i__2, &ab[*ku + 1 - j + j1 + j * ab_dim1], &c__1, &afb[*
		    kl + *ku + 1 - j + j1 + j * afb_dim1], &c__1);
/* L70: */
	}

	dgbtrf_(n, n, kl, ku, &afb[afb_offset], ldafb, &ipiv[1], info);

/*        Return if INFO is non-zero. */

	if (*info > 0) {

/*           Compute the reciprocal pivot growth factor of the */
/*           leading rank-deficient INFO columns of A. */

	    anorm = 0.;
	    i__1 = *info;
	    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
		i__2 = *ku + 2 - j;
/* Computing MIN */
		i__4 = *n + *ku + 1 - j, i__5 = *kl + *ku + 1;
		i__3 = f2cmin(i__4,i__5);
		for (i__ = f2cmax(i__2,1); i__ <= i__3; ++i__) {
/* Computing MAX */
		    d__2 = anorm, d__3 = (d__1 = ab[i__ + j * ab_dim1], abs(
			    d__1));
		    anorm = f2cmax(d__2,d__3);
/* L80: */
		}
/* L90: */
	    }
/* Computing MIN */
	    i__3 = *info - 1, i__2 = *kl + *ku;
	    i__1 = f2cmin(i__3,i__2);
/* Computing MAX */
	    i__4 = 1, i__5 = *kl + *ku + 2 - *info;
	    rpvgrw = dlantb_("M", "U", "N", info, &i__1, &afb[f2cmax(i__4,i__5) 
		    + afb_dim1], ldafb, &work[1]);
	    if (rpvgrw == 0.) {
		rpvgrw = 1.;
	    } else {
		rpvgrw = anorm / rpvgrw;
	    }
	    work[1] = rpvgrw;
	    *rcond = 0.;
	    return 0;
	}
    }

/*     Compute the norm of the matrix A and the */
/*     reciprocal pivot growth factor RPVGRW. */

    if (notran) {
	*(unsigned char *)norm = '1';
    } else {
	*(unsigned char *)norm = 'I';
    }
    anorm = dlangb_(norm, n, kl, ku, &ab[ab_offset], ldab, &work[1]);
    i__1 = *kl + *ku;
    rpvgrw = dlantb_("M", "U", "N", n, &i__1, &afb[afb_offset], ldafb, &work[
	    1]);
    if (rpvgrw == 0.) {
	rpvgrw = 1.;
    } else {
	rpvgrw = dlangb_("M", n, kl, ku, &ab[ab_offset], ldab, &work[1]) / rpvgrw;
    }

/*     Compute the reciprocal of the condition number of A. */

    dgbcon_(norm, n, kl, ku, &afb[afb_offset], ldafb, &ipiv[1], &anorm, rcond,
	     &work[1], &iwork[1], info);

/*     Compute the solution matrix X. */

    dlacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx);
    dgbtrs_(trans, n, kl, ku, nrhs, &afb[afb_offset], ldafb, &ipiv[1], &x[
	    x_offset], ldx, info);

/*     Use iterative refinement to improve the computed solution and */
/*     compute error bounds and backward error estimates for it. */

    dgbrfs_(trans, n, kl, ku, nrhs, &ab[ab_offset], ldab, &afb[afb_offset], 
	    ldafb, &ipiv[1], &b[b_offset], ldb, &x[x_offset], ldx, &ferr[1], &
	    berr[1], &work[1], &iwork[1], info);

/*     Transform the solution matrix X to a solution of the original */
/*     system. */

    if (notran) {
	if (colequ) {
	    i__1 = *nrhs;
	    for (j = 1; j <= i__1; ++j) {
		i__3 = *n;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    x[i__ + j * x_dim1] = c__[i__] * x[i__ + j * x_dim1];
/* L100: */
		}
/* L110: */
	    }
	    i__1 = *nrhs;
	    for (j = 1; j <= i__1; ++j) {
		ferr[j] /= colcnd;
/* L120: */
	    }
	}
    } else if (rowequ) {
	i__1 = *nrhs;
	for (j = 1; j <= i__1; ++j) {
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		x[i__ + j * x_dim1] = r__[i__] * x[i__ + j * x_dim1];
/* L130: */
	    }
/* L140: */
	}
	i__1 = *nrhs;
	for (j = 1; j <= i__1; ++j) {
	    ferr[j] /= rowcnd;
/* L150: */
	}
    }

/*     Set INFO = N+1 if the matrix is singular to working precision. */

    if (*rcond < dlamch_("Epsilon")) {
	*info = *n + 1;
    }

    work[1] = rpvgrw;
    return 0;

/*     End of DGBSVX */

} /* dgbsvx_ */

