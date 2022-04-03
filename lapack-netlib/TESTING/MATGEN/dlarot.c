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

#if defined(_WIN64)
typedef long long BLASLONG;
typedef unsigned long long BLASULONG;
#else
typedef long BLASLONG;
typedef unsigned long BLASULONG;
#endif

#ifdef LAPACK_ILP64
typedef BLASLONG blasint;
#if defined(_WIN64)
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
#define mycycle_() continue;
#define myceiling_(w) ceil(w)
#define myhuge_(w) HUGE_VAL
//#define mymaxloc_(w,s,e,n) {if (sizeof(*(w)) == sizeof(double)) dmaxloc_((w),*(s),*(e),n); else dmaxloc_((w),*(s),*(e),n);}
#define mymaxloc_(w,s,e,n) dmaxloc_(w,*(s),*(e),n)

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

static integer c__4 = 4;
static integer c__8 = 8;
static integer c__1 = 1;

/* > \brief \b DLAROT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAROT( LROWS, LLEFT, LRIGHT, NL, C, S, A, LDA, XLEFT, */
/*                          XRIGHT ) */

/*       LOGICAL            LLEFT, LRIGHT, LROWS */
/*       INTEGER            LDA, NL */
/*       DOUBLE PRECISION   C, S, XLEFT, XRIGHT */
/*       DOUBLE PRECISION   A( * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    DLAROT applies a (Givens) rotation to two adjacent rows or */
/* >    columns, where one element of the first and/or last column/row */
/* >    for use on matrices stored in some format other than GE, so */
/* >    that elements of the matrix may be used or modified for which */
/* >    no array element is provided. */
/* > */
/* >    One example is a symmetric matrix in SB format (bandwidth=4), for */
/* >    which UPLO='L':  Two adjacent rows will have the format: */
/* > */
/* >    row j:     C> C> C> C> C> .  .  .  . */
/* >    row j+1:      C> C> C> C> C> .  .  .  . */
/* > */
/* >    '*' indicates elements for which storage is provided, */
/* >    '.' indicates elements for which no storage is provided, but */
/* >    are not necessarily zero; their values are determined by */
/* >    symmetry.  ' ' indicates elements which are necessarily zero, */
/* >     and have no storage provided. */
/* > */
/* >    Those columns which have two '*'s can be handled by DROT. */
/* >    Those columns which have no '*'s can be ignored, since as long */
/* >    as the Givens rotations are carefully applied to preserve */
/* >    symmetry, their values are determined. */
/* >    Those columns which have one '*' have to be handled separately, */
/* >    by using separate variables "p" and "q": */
/* > */
/* >    row j:     C> C> C> C> C> p  .  .  . */
/* >    row j+1:   q  C> C> C> C> C> .  .  .  . */
/* > */
/* >    The element p would have to be set correctly, then that column */
/* >    is rotated, setting p to its new value.  The next call to */
/* >    DLAROT would rotate columns j and j+1, using p, and restore */
/* >    symmetry.  The element q would start out being zero, and be */
/* >    made non-zero by the rotation.  Later, rotations would presumably */
/* >    be chosen to zero q out. */
/* > */
/* >    Typical Calling Sequences: rotating the i-th and (i+1)-st rows. */
/* >    ------- ------- --------- */
/* > */
/* >      General dense matrix: */
/* > */
/* >              CALL DLAROT(.TRUE.,.FALSE.,.FALSE., N, C,S, */
/* >                      A(i,1),LDA, DUMMY, DUMMY) */
/* > */
/* >      General banded matrix in GB format: */
/* > */
/* >              j = MAX(1, i-KL ) */
/* >              NL = MIN( N, i+KU+1 ) + 1-j */
/* >              CALL DLAROT( .TRUE., i-KL.GE.1, i+KU.LT.N, NL, C,S, */
/* >                      A(KU+i+1-j,j),LDA-1, XLEFT, XRIGHT ) */
/* > */
/* >              [ note that i+1-j is just MIN(i,KL+1) ] */
/* > */
/* >      Symmetric banded matrix in SY format, bandwidth K, */
/* >      lower triangle only: */
/* > */
/* >              j = MAX(1, i-K ) */
/* >              NL = MIN( K+1, i ) + 1 */
/* >              CALL DLAROT( .TRUE., i-K.GE.1, .TRUE., NL, C,S, */
/* >                      A(i,j), LDA, XLEFT, XRIGHT ) */
/* > */
/* >      Same, but upper triangle only: */
/* > */
/* >              NL = MIN( K+1, N-i ) + 1 */
/* >              CALL DLAROT( .TRUE., .TRUE., i+K.LT.N, NL, C,S, */
/* >                      A(i,i), LDA, XLEFT, XRIGHT ) */
/* > */
/* >      Symmetric banded matrix in SB format, bandwidth K, */
/* >      lower triangle only: */
/* > */
/* >              [ same as for SY, except:] */
/* >                  . . . . */
/* >                      A(i+1-j,j), LDA-1, XLEFT, XRIGHT ) */
/* > */
/* >              [ note that i+1-j is just MIN(i,K+1) ] */
/* > */
/* >      Same, but upper triangle only: */
/* >                   . . . */
/* >                      A(K+1,i), LDA-1, XLEFT, XRIGHT ) */
/* > */
/* >      Rotating columns is just the transpose of rotating rows, except */
/* >      for GB and SB: (rotating columns i and i+1) */
/* > */
/* >      GB: */
/* >              j = MAX(1, i-KU ) */
/* >              NL = MIN( N, i+KL+1 ) + 1-j */
/* >              CALL DLAROT( .TRUE., i-KU.GE.1, i+KL.LT.N, NL, C,S, */
/* >                      A(KU+j+1-i,i),LDA-1, XTOP, XBOTTM ) */
/* > */
/* >              [note that KU+j+1-i is just MAX(1,KU+2-i)] */
/* > */
/* >      SB: (upper triangle) */
/* > */
/* >                   . . . . . . */
/* >                      A(K+j+1-i,i),LDA-1, XTOP, XBOTTM ) */
/* > */
/* >      SB: (lower triangle) */
/* > */
/* >                   . . . . . . */
/* >                      A(1,i),LDA-1, XTOP, XBOTTM ) */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \verbatim */
/* >  LROWS  - LOGICAL */
/* >           If .TRUE., then DLAROT will rotate two rows.  If .FALSE., */
/* >           then it will rotate two columns. */
/* >           Not modified. */
/* > */
/* >  LLEFT  - LOGICAL */
/* >           If .TRUE., then XLEFT will be used instead of the */
/* >           corresponding element of A for the first element in the */
/* >           second row (if LROWS=.FALSE.) or column (if LROWS=.TRUE.) */
/* >           If .FALSE., then the corresponding element of A will be */
/* >           used. */
/* >           Not modified. */
/* > */
/* >  LRIGHT - LOGICAL */
/* >           If .TRUE., then XRIGHT will be used instead of the */
/* >           corresponding element of A for the last element in the */
/* >           first row (if LROWS=.FALSE.) or column (if LROWS=.TRUE.) If */
/* >           .FALSE., then the corresponding element of A will be used. */
/* >           Not modified. */
/* > */
/* >  NL     - INTEGER */
/* >           The length of the rows (if LROWS=.TRUE.) or columns (if */
/* >           LROWS=.FALSE.) to be rotated.  If XLEFT and/or XRIGHT are */
/* >           used, the columns/rows they are in should be included in */
/* >           NL, e.g., if LLEFT = LRIGHT = .TRUE., then NL must be at */
/* >           least 2.  The number of rows/columns to be rotated */
/* >           exclusive of those involving XLEFT and/or XRIGHT may */
/* >           not be negative, i.e., NL minus how many of LLEFT and */
/* >           LRIGHT are .TRUE. must be at least zero; if not, XERBLA */
/* >           will be called. */
/* >           Not modified. */
/* > */
/* >  C, S   - DOUBLE PRECISION */
/* >           Specify the Givens rotation to be applied.  If LROWS is */
/* >           true, then the matrix ( c  s ) */
/* >                                 (-s  c )  is applied from the left; */
/* >           if false, then the transpose thereof is applied from the */
/* >           right.  For a Givens rotation, C**2 + S**2 should be 1, */
/* >           but this is not checked. */
/* >           Not modified. */
/* > */
/* >  A      - DOUBLE PRECISION array. */
/* >           The array containing the rows/columns to be rotated.  The */
/* >           first element of A should be the upper left element to */
/* >           be rotated. */
/* >           Read and modified. */
/* > */
/* >  LDA    - INTEGER */
/* >           The "effective" leading dimension of A.  If A contains */
/* >           a matrix stored in GE or SY format, then this is just */
/* >           the leading dimension of A as dimensioned in the calling */
/* >           routine.  If A contains a matrix stored in band (GB or SB) */
/* >           format, then this should be *one less* than the leading */
/* >           dimension used in the calling routine.  Thus, if */
/* >           A were dimensioned A(LDA,*) in DLAROT, then A(1,j) would */
/* >           be the j-th element in the first of the two rows */
/* >           to be rotated, and A(2,j) would be the j-th in the second, */
/* >           regardless of how the array may be stored in the calling */
/* >           routine.  [A cannot, however, actually be dimensioned thus, */
/* >           since for band format, the row number may exceed LDA, which */
/* >           is not legal FORTRAN.] */
/* >           If LROWS=.TRUE., then LDA must be at least 1, otherwise */
/* >           it must be at least NL minus the number of .TRUE. values */
/* >           in XLEFT and XRIGHT. */
/* >           Not modified. */
/* > */
/* >  XLEFT  - DOUBLE PRECISION */
/* >           If LLEFT is .TRUE., then XLEFT will be used and modified */
/* >           instead of A(2,1) (if LROWS=.TRUE.) or A(1,2) */
/* >           (if LROWS=.FALSE.). */
/* >           Read and modified. */
/* > */
/* >  XRIGHT - DOUBLE PRECISION */
/* >           If LRIGHT is .TRUE., then XRIGHT will be used and modified */
/* >           instead of A(1,NL) (if LROWS=.TRUE.) or A(NL,1) */
/* >           (if LROWS=.FALSE.). */
/* >           Read and modified. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup double_matgen */

/*  ===================================================================== */
/* Subroutine */ int dlarot_(logical *lrows, logical *lleft, logical *lright, 
	integer *nl, doublereal *c__, doublereal *s, doublereal *a, integer *
	lda, doublereal *xleft, doublereal *xright)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer iinc;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    integer inext, ix, iy, nt;
    doublereal xt[2], yt[2];
    extern /* Subroutine */ int xerbla_(char *, integer *);
    integer iyt;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */


/*  ===================================================================== */


/*     Set up indices, arrays for ends */

    /* Parameter adjustments */
    --a;

    /* Function Body */
    if (*lrows) {
	iinc = *lda;
	inext = 1;
    } else {
	iinc = 1;
	inext = *lda;
    }

    if (*lleft) {
	nt = 1;
	ix = iinc + 1;
	iy = *lda + 2;
	xt[0] = a[1];
	yt[0] = *xleft;
    } else {
	nt = 0;
	ix = 1;
	iy = inext + 1;
    }

    if (*lright) {
	iyt = inext + 1 + (*nl - 1) * iinc;
	++nt;
	xt[nt - 1] = *xright;
	yt[nt - 1] = a[iyt];
    }

/*     Check for errors */

    if (*nl < nt) {
	xerbla_("DLAROT", &c__4);
	return 0;
    }
    if (*lda <= 0 || ! (*lrows) && *lda < *nl - nt) {
	xerbla_("DLAROT", &c__8);
	return 0;
    }

/*     Rotate */

    i__1 = *nl - nt;
    drot_(&i__1, &a[ix], &iinc, &a[iy], &iinc, c__, s);
    drot_(&nt, xt, &c__1, yt, &c__1, c__, s);

/*     Stuff values back into XLEFT, XRIGHT, etc. */

    if (*lleft) {
	a[1] = xt[0];
	*xleft = yt[0];
    }

    if (*lright) {
	*xright = xt[nt - 1];
	a[iyt] = yt[nt - 1];
    }

    return 0;

/*     End of DLAROT */

} /* dlarot_ */

