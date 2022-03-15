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

#if defined(OS_WINDOWS) && defined(__64BIT__)
typedef long long BLASLONG;
typedef unsigned long long BLASULONG;
#else
typedef long BLASLONG;
typedef unsigned long BLASULONG;
#endif

#ifdef LAPACK_ILP64
typedef BLASLONG blasint;
#if defined(OS_WINDOWS) && defined(__64BIT__)
#define blasabs(x) llabs(x)
#else
#define blasabs(x) labs(x)
#endif
#else
typedef int blasint;
#define blasabs(x) abs(x)
#endif

typedef blasint integer;

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

static doublecomplex c_b1 = {-1.,0.};
static integer c__1 = 1;

/* > \brief \b ZUNBDB4 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZUNBDB4 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunbdb4
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunbdb4
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunbdb4
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZUNBDB4( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI, */
/*                           TAUP1, TAUP2, TAUQ1, PHANTOM, WORK, LWORK, */
/*                           INFO ) */

/*       INTEGER            INFO, LWORK, M, P, Q, LDX11, LDX21 */
/*       DOUBLE PRECISION   PHI(*), THETA(*) */
/*       COMPLEX*16         PHANTOM(*), TAUP1(*), TAUP2(*), TAUQ1(*), */
/*      $                   WORK(*), X11(LDX11,*), X21(LDX21,*) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* >\verbatim */
/* > */
/* > ZUNBDB4 simultaneously bidiagonalizes the blocks of a tall and skinny */
/* > matrix X with orthonomal columns: */
/* > */
/* >                            [ B11 ] */
/* >      [ X11 ]   [ P1 |    ] [  0  ] */
/* >      [-----] = [---------] [-----] Q1**T . */
/* >      [ X21 ]   [    | P2 ] [ B21 ] */
/* >                            [  0  ] */
/* > */
/* > X11 is P-by-Q, and X21 is (M-P)-by-Q. M-Q must be no larger than P, */
/* > M-P, or Q. Routines ZUNBDB1, ZUNBDB2, and ZUNBDB3 handle cases in */
/* > which M-Q is not the minimum dimension. */
/* > */
/* > The unitary matrices P1, P2, and Q1 are P-by-P, (M-P)-by-(M-P), */
/* > and (M-Q)-by-(M-Q), respectively. They are represented implicitly by */
/* > Householder vectors. */
/* > */
/* > B11 and B12 are (M-Q)-by-(M-Q) bidiagonal matrices represented */
/* > implicitly by angles THETA, PHI. */
/* > */
/* >\endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           The number of rows X11 plus the number of rows in X21. */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* >          P is INTEGER */
/* >           The number of rows in X11. 0 <= P <= M. */
/* > \endverbatim */
/* > */
/* > \param[in] Q */
/* > \verbatim */
/* >          Q is INTEGER */
/* >           The number of columns in X11 and X21. 0 <= Q <= M and */
/* >           M-Q <= f2cmin(P,M-P,Q). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X11 */
/* > \verbatim */
/* >          X11 is COMPLEX*16 array, dimension (LDX11,Q) */
/* >           On entry, the top block of the matrix X to be reduced. On */
/* >           exit, the columns of tril(X11) specify reflectors for P1 and */
/* >           the rows of triu(X11,1) specify reflectors for Q1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX11 */
/* > \verbatim */
/* >          LDX11 is INTEGER */
/* >           The leading dimension of X11. LDX11 >= P. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X21 */
/* > \verbatim */
/* >          X21 is COMPLEX*16 array, dimension (LDX21,Q) */
/* >           On entry, the bottom block of the matrix X to be reduced. On */
/* >           exit, the columns of tril(X21) specify reflectors for P2. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX21 */
/* > \verbatim */
/* >          LDX21 is INTEGER */
/* >           The leading dimension of X21. LDX21 >= M-P. */
/* > \endverbatim */
/* > */
/* > \param[out] THETA */
/* > \verbatim */
/* >          THETA is DOUBLE PRECISION array, dimension (Q) */
/* >           The entries of the bidiagonal blocks B11, B21 are defined by */
/* >           THETA and PHI. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] PHI */
/* > \verbatim */
/* >          PHI is DOUBLE PRECISION array, dimension (Q-1) */
/* >           The entries of the bidiagonal blocks B11, B21 are defined by */
/* >           THETA and PHI. See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP1 */
/* > \verbatim */
/* >          TAUP1 is COMPLEX*16 array, dimension (P) */
/* >           The scalar factors of the elementary reflectors that define */
/* >           P1. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUP2 */
/* > \verbatim */
/* >          TAUP2 is COMPLEX*16 array, dimension (M-P) */
/* >           The scalar factors of the elementary reflectors that define */
/* >           P2. */
/* > \endverbatim */
/* > */
/* > \param[out] TAUQ1 */
/* > \verbatim */
/* >          TAUQ1 is COMPLEX*16 array, dimension (Q) */
/* >           The scalar factors of the elementary reflectors that define */
/* >           Q1. */
/* > \endverbatim */
/* > */
/* > \param[out] PHANTOM */
/* > \verbatim */
/* >          PHANTOM is COMPLEX*16 array, dimension (M) */
/* >           The routine computes an M-by-1 column vector Y that is */
/* >           orthogonal to the columns of [ X11; X21 ]. PHANTOM(1:P) and */
/* >           PHANTOM(P+1:M) contain Householder vectors for Y(1:P) and */
/* >           Y(P+1:M), respectively. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (LWORK) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >           The dimension of the array WORK. LWORK >= M-Q. */
/* > */
/* >           If LWORK = -1, then a workspace query is assumed; the routine */
/* >           only calculates the optimal size of the WORK array, returns */
/* >           this value as the first entry of the WORK array, and no error */
/* >           message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >           = 0:  successful exit. */
/* >           < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date July 2012 */

/* > \ingroup complex16OTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The upper-bidiagonal blocks B11, B21 are represented implicitly by */
/* >  angles THETA(1), ..., THETA(Q) and PHI(1), ..., PHI(Q-1). Every entry */
/* >  in each bidiagonal band is a product of a sine or cosine of a THETA */
/* >  with a sine or cosine of a PHI. See [1] or ZUNCSD for details. */
/* > */
/* >  P1, P2, and Q1 are represented as products of elementary reflectors. */
/* >  See ZUNCSD2BY1 for details on generating P1, P2, and Q1 using ZUNGQR */
/* >  and ZUNGLQ. */
/* > \endverbatim */

/* > \par References: */
/*  ================ */
/* > */
/* >  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer. */
/* >      Algorithms, 50(1):33-65, 2009. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zunbdb4_(integer *m, integer *p, integer *q, 
	doublecomplex *x11, integer *ldx11, doublecomplex *x21, integer *
	ldx21, doublereal *theta, doublereal *phi, doublecomplex *taup1, 
	doublecomplex *taup2, doublecomplex *tauq1, doublecomplex *phantom, 
	doublecomplex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer x11_dim1, x11_offset, x21_dim1, x21_offset, i__1, i__2, i__3, 
	    i__4;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Local variables */
    integer lworkmin, lworkopt;
    doublereal c__;
    integer i__, j;
    doublereal s;
    integer ilarf, llarf;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    integer childinfo;
    extern /* Subroutine */ int zlarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *), zdrot_(integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, doublereal *);
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zlacgv_(
	    integer *, doublecomplex *, integer *);
    logical lquery;
    integer iorbdb5, lorbdb5;
    extern /* Subroutine */ int zunbdb5_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *), zlarfgp_(integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *);


/*  -- LAPACK computational routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     July 2012 */


/*  ==================================================================== */


/*     Test input arguments */

    /* Parameter adjustments */
    x11_dim1 = *ldx11;
    x11_offset = 1 + x11_dim1 * 1;
    x11 -= x11_offset;
    x21_dim1 = *ldx21;
    x21_offset = 1 + x21_dim1 * 1;
    x21 -= x21_offset;
    --theta;
    --phi;
    --taup1;
    --taup2;
    --tauq1;
    --phantom;
    --work;

    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;

    if (*m < 0) {
	*info = -1;
    } else if (*p < *m - *q || *m - *p < *m - *q) {
	*info = -2;
    } else if (*q < *m - *q || *q > *m) {
	*info = -3;
    } else if (*ldx11 < f2cmax(1,*p)) {
	*info = -5;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *m - *p;
	if (*ldx21 < f2cmax(i__1,i__2)) {
	    *info = -7;
	}
    }

/*     Compute workspace */

    if (*info == 0) {
	ilarf = 2;
/* Computing MAX */
	i__1 = *q - 1, i__2 = *p - 1, i__1 = f2cmax(i__1,i__2), i__2 = *m - *p - 
		1;
	llarf = f2cmax(i__1,i__2);
	iorbdb5 = 2;
	lorbdb5 = *q;
	lworkopt = ilarf + llarf - 1;
/* Computing MAX */
	i__1 = lworkopt, i__2 = iorbdb5 + lorbdb5 - 1;
	lworkopt = f2cmax(i__1,i__2);
	lworkmin = lworkopt;
	work[1].r = (doublereal) lworkopt, work[1].i = 0.;
	if (*lwork < lworkmin && ! lquery) {
	    *info = -14;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZUNBDB4", &i__1, (ftnlen)7);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Reduce columns 1, ..., M-Q of X11 and X21 */

    i__1 = *m - *q;
    for (i__ = 1; i__ <= i__1; ++i__) {

	if (i__ == 1) {
	    i__2 = *m;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = j;
		phantom[i__3].r = 0., phantom[i__3].i = 0.;
	    }
	    i__2 = *m - *p;
	    zunbdb5_(p, &i__2, q, &phantom[1], &c__1, &phantom[*p + 1], &c__1,
		     &x11[x11_offset], ldx11, &x21[x21_offset], ldx21, &work[
		    iorbdb5], &lorbdb5, &childinfo);
	    zscal_(p, &c_b1, &phantom[1], &c__1);
	    zlarfgp_(p, &phantom[1], &phantom[2], &c__1, &taup1[1]);
	    i__2 = *m - *p;
	    zlarfgp_(&i__2, &phantom[*p + 1], &phantom[*p + 2], &c__1, &taup2[
		    1]);
	    theta[i__] = atan2((doublereal) phantom[1].r, (doublereal) 
		    phantom[*p + 1].r);
	    c__ = cos(theta[i__]);
	    s = sin(theta[i__]);
	    phantom[1].r = 1., phantom[1].i = 0.;
	    i__2 = *p + 1;
	    phantom[i__2].r = 1., phantom[i__2].i = 0.;
	    d_cnjg(&z__1, &taup1[1]);
	    zlarf_("L", p, q, &phantom[1], &c__1, &z__1, &x11[x11_offset], 
		    ldx11, &work[ilarf]);
	    i__2 = *m - *p;
	    d_cnjg(&z__1, &taup2[1]);
	    zlarf_("L", &i__2, q, &phantom[*p + 1], &c__1, &z__1, &x21[
		    x21_offset], ldx21, &work[ilarf]);
	} else {
	    i__2 = *p - i__ + 1;
	    i__3 = *m - *p - i__ + 1;
	    i__4 = *q - i__ + 1;
	    zunbdb5_(&i__2, &i__3, &i__4, &x11[i__ + (i__ - 1) * x11_dim1], &
		    c__1, &x21[i__ + (i__ - 1) * x21_dim1], &c__1, &x11[i__ + 
		    i__ * x11_dim1], ldx11, &x21[i__ + i__ * x21_dim1], ldx21,
		     &work[iorbdb5], &lorbdb5, &childinfo);
	    i__2 = *p - i__ + 1;
	    zscal_(&i__2, &c_b1, &x11[i__ + (i__ - 1) * x11_dim1], &c__1);
	    i__2 = *p - i__ + 1;
	    zlarfgp_(&i__2, &x11[i__ + (i__ - 1) * x11_dim1], &x11[i__ + 1 + (
		    i__ - 1) * x11_dim1], &c__1, &taup1[i__]);
	    i__2 = *m - *p - i__ + 1;
	    zlarfgp_(&i__2, &x21[i__ + (i__ - 1) * x21_dim1], &x21[i__ + 1 + (
		    i__ - 1) * x21_dim1], &c__1, &taup2[i__]);
	    theta[i__] = atan2((doublereal) x11[i__ + (i__ - 1) * x11_dim1].r,
		     (doublereal) x21[i__ + (i__ - 1) * x21_dim1].r);
	    c__ = cos(theta[i__]);
	    s = sin(theta[i__]);
	    i__2 = i__ + (i__ - 1) * x11_dim1;
	    x11[i__2].r = 1., x11[i__2].i = 0.;
	    i__2 = i__ + (i__ - 1) * x21_dim1;
	    x21[i__2].r = 1., x21[i__2].i = 0.;
	    i__2 = *p - i__ + 1;
	    i__3 = *q - i__ + 1;
	    d_cnjg(&z__1, &taup1[i__]);
	    zlarf_("L", &i__2, &i__3, &x11[i__ + (i__ - 1) * x11_dim1], &c__1,
		     &z__1, &x11[i__ + i__ * x11_dim1], ldx11, &work[ilarf]);
	    i__2 = *m - *p - i__ + 1;
	    i__3 = *q - i__ + 1;
	    d_cnjg(&z__1, &taup2[i__]);
	    zlarf_("L", &i__2, &i__3, &x21[i__ + (i__ - 1) * x21_dim1], &c__1,
		     &z__1, &x21[i__ + i__ * x21_dim1], ldx21, &work[ilarf]);
	}

	i__2 = *q - i__ + 1;
	d__1 = -c__;
	zdrot_(&i__2, &x11[i__ + i__ * x11_dim1], ldx11, &x21[i__ + i__ * 
		x21_dim1], ldx21, &s, &d__1);
	i__2 = *q - i__ + 1;
	zlacgv_(&i__2, &x21[i__ + i__ * x21_dim1], ldx21);
	i__2 = *q - i__ + 1;
	zlarfgp_(&i__2, &x21[i__ + i__ * x21_dim1], &x21[i__ + (i__ + 1) * 
		x21_dim1], ldx21, &tauq1[i__]);
	i__2 = i__ + i__ * x21_dim1;
	c__ = x21[i__2].r;
	i__2 = i__ + i__ * x21_dim1;
	x21[i__2].r = 1., x21[i__2].i = 0.;
	i__2 = *p - i__;
	i__3 = *q - i__ + 1;
	zlarf_("R", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], ldx21, &tauq1[
		i__], &x11[i__ + 1 + i__ * x11_dim1], ldx11, &work[ilarf]);
	i__2 = *m - *p - i__;
	i__3 = *q - i__ + 1;
	zlarf_("R", &i__2, &i__3, &x21[i__ + i__ * x21_dim1], ldx21, &tauq1[
		i__], &x21[i__ + 1 + i__ * x21_dim1], ldx21, &work[ilarf]);
	i__2 = *q - i__ + 1;
	zlacgv_(&i__2, &x21[i__ + i__ * x21_dim1], ldx21);
	if (i__ < *m - *q) {
	    i__2 = *p - i__;
/* Computing 2nd power */
	    d__1 = dznrm2_(&i__2, &x11[i__ + 1 + i__ * x11_dim1], &c__1);
	    i__3 = *m - *p - i__;
/* Computing 2nd power */
	    d__2 = dznrm2_(&i__3, &x21[i__ + 1 + i__ * x21_dim1], &c__1);
	    s = sqrt(d__1 * d__1 + d__2 * d__2);
	    phi[i__] = atan2(s, c__);
	}

    }

/*     Reduce the bottom-right portion of X11 to [ I 0 ] */

    i__1 = *p;
    for (i__ = *m - *q + 1; i__ <= i__1; ++i__) {
	i__2 = *q - i__ + 1;
	zlacgv_(&i__2, &x11[i__ + i__ * x11_dim1], ldx11);
	i__2 = *q - i__ + 1;
	zlarfgp_(&i__2, &x11[i__ + i__ * x11_dim1], &x11[i__ + (i__ + 1) * 
		x11_dim1], ldx11, &tauq1[i__]);
	i__2 = i__ + i__ * x11_dim1;
	x11[i__2].r = 1., x11[i__2].i = 0.;
	i__2 = *p - i__;
	i__3 = *q - i__ + 1;
	zlarf_("R", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], ldx11, &tauq1[
		i__], &x11[i__ + 1 + i__ * x11_dim1], ldx11, &work[ilarf]);
	i__2 = *q - *p;
	i__3 = *q - i__ + 1;
	zlarf_("R", &i__2, &i__3, &x11[i__ + i__ * x11_dim1], ldx11, &tauq1[
		i__], &x21[*m - *q + 1 + i__ * x21_dim1], ldx21, &work[ilarf]);
	i__2 = *q - i__ + 1;
	zlacgv_(&i__2, &x11[i__ + i__ * x11_dim1], ldx11);
    }

/*     Reduce the bottom-right portion of X21 to [ 0 I ] */

    i__1 = *q;
    for (i__ = *p + 1; i__ <= i__1; ++i__) {
	i__2 = *q - i__ + 1;
	zlacgv_(&i__2, &x21[*m - *q + i__ - *p + i__ * x21_dim1], ldx21);
	i__2 = *q - i__ + 1;
	zlarfgp_(&i__2, &x21[*m - *q + i__ - *p + i__ * x21_dim1], &x21[*m - *
		q + i__ - *p + (i__ + 1) * x21_dim1], ldx21, &tauq1[i__]);
	i__2 = *m - *q + i__ - *p + i__ * x21_dim1;
	x21[i__2].r = 1., x21[i__2].i = 0.;
	i__2 = *q - i__;
	i__3 = *q - i__ + 1;
	zlarf_("R", &i__2, &i__3, &x21[*m - *q + i__ - *p + i__ * x21_dim1], 
		ldx21, &tauq1[i__], &x21[*m - *q + i__ - *p + 1 + i__ * 
		x21_dim1], ldx21, &work[ilarf]);
	i__2 = *q - i__ + 1;
	zlacgv_(&i__2, &x21[*m - *q + i__ - *p + i__ * x21_dim1], ldx21);
    }

    return 0;

/*     End of ZUNBDB4 */

} /* zunbdb4_ */

