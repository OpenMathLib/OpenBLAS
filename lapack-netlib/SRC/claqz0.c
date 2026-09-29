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
#ifdef _MSC_VER
static inline _Fcomplex Cf(complex *z) {_Fcomplex zz={z->r , z->i}; return zz;}
static inline _Fcomplex * _pCf(complex *z) {return (_Fcomplex*)z;}
#else
static inline _Complex float Cf(complex *z) {return z->r + z->i*_Complex_I;}
static inline _Complex float * _pCf(complex *z) {return (_Complex float*)z;}
#endif
#define pCf(z) (*_pCf(z))
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
#ifdef _MSC_VER
#define c_div(c, a, b) {Cf(c)._Val[0] = (Cf(a)._Val[0]/Cf(b)._Val[0]); Cf(c)._Val[1]=(Cf(a)._Val[1]/Cf(b)._Val[1]);}
#define z_div(c, a, b) {Cd(c)._Val[0] = (Cd(a)._Val[0]/Cd(b)._Val[0]); Cd(c)._Val[1]=(Cd(a)._Val[1]/Cd(b)._Val[1]);}
#else
#define c_div(c, a, b) {pCf(c) = Cf(a)/Cf(b);}
#define z_div(c, a, b) {pCd(c) = Cd(a)/Cd(b);}
#endif
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
#define r_cnjg(R, Z) { pCf(R) = conjf(Cf(Z)); }
#define d_cos(x) (cos(*(x)))
#define d_cosh(x) (cosh(*(x)))
#define d_dim(__a, __b) ( *(__a) > *(__b) ? *(__a) - *(__b) : 0.0 )
#define d_exp(x) (exp(*(x)))
#define d_imag(z) (cimag(Cd(z)))
#define r_imag(z) (cimagf(Cf(z)))
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
#define z_abs(z) (cabs(Cd(z)))
#define z_exp(R, Z) {pCd(R) = cexp(Cd(Z));}
#define z_sqrt(R, Z) {pCd(R) = csqrt(Cd(Z));}
#define myexit_() break;
#define mycycle_() continue;
#define myceiling_(w) {ceil(w)}
#define myhuge_(w) {HUGE_VAL}
//#define mymaxloc_(w,s,e,n) {if (sizeof(*(w)) == sizeof(double)) dmaxloc_((w),*(s),*(e),n); else dmaxloc_((w),*(s),*(e),n);}
#define mymaxloc_(w,s,e,n) dmaxloc_(w,*(s),*(e),n)

/* procedure parameter types for -A and -C++ */

#define F2C_proc_par_types 1
#ifdef __cplusplus
typedef logical (*L_fp)(...);
#else
typedef logical (*L_fp)();
#endif
#if 0
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
#ifdef _MSC_VER
static _Fcomplex cpow_ui(complex x, integer n) {
	complex pow={1.0,0.0}; unsigned long int u;
		if(n != 0) {
		if(n < 0) n = -n, x.r = 1/x.r, x.i=1/x.i;
		for(u = n; ; ) {
			if(u & 01) pow.r *= x.r, pow.i *= x.i;
			if(u >>= 1) x.r *= x.r, x.i *= x.i;
			else break;
		}
	}
	_Fcomplex p={pow.r, pow.i};
	return p;
}
#else
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
#endif
#ifdef _MSC_VER
static _Dcomplex zpow_ui(_Dcomplex x, integer n) {
	_Dcomplex pow={1.0,0.0}; unsigned long int u;
	if(n != 0) {
		if(n < 0) n = -n, x._Val[0] = 1/x._Val[0], x._Val[1] =1/x._Val[1];
		for(u = n; ; ) {
			if(u & 01) pow._Val[0] *= x._Val[0], pow._Val[1] *= x._Val[1];
			if(u >>= 1) x._Val[0] *= x._Val[0], x._Val[1] *= x._Val[1];
			else break;
		}
	}
	_Dcomplex p = {pow._Val[0], pow._Val[1]};
	return p;
}
#else
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
#endif
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
#ifdef _MSC_VER
	_Fcomplex zdotc = {0.0, 0.0};
	if (incx == 1 && incy == 1) {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc._Val[0] += conjf(Cf(&x[i]))._Val[0] * Cf(&y[i])._Val[0];
			zdotc._Val[1] += conjf(Cf(&x[i]))._Val[1] * Cf(&y[i])._Val[1];
		}
	} else {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc._Val[0] += conjf(Cf(&x[i*incx]))._Val[0] * Cf(&y[i*incy])._Val[0];
			zdotc._Val[1] += conjf(Cf(&x[i*incx]))._Val[1] * Cf(&y[i*incy])._Val[1];
		}
	}
	pCf(z) = zdotc;
}
#else
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
#endif
static inline void zdotc_(doublecomplex *z, integer *n_, doublecomplex *x, integer *incx_, doublecomplex *y, integer *incy_) {
	integer n = *n_, incx = *incx_, incy = *incy_, i;
#ifdef _MSC_VER
	_Dcomplex zdotc = {0.0, 0.0};
	if (incx == 1 && incy == 1) {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc._Val[0] += conj(Cd(&x[i]))._Val[0] * Cd(&y[i])._Val[0];
			zdotc._Val[1] += conj(Cd(&x[i]))._Val[1] * Cd(&y[i])._Val[1];
		}
	} else {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc._Val[0] += conj(Cd(&x[i*incx]))._Val[0] * Cd(&y[i*incy])._Val[0];
			zdotc._Val[1] += conj(Cd(&x[i*incx]))._Val[1] * Cd(&y[i*incy])._Val[1];
		}
	}
	pCd(z) = zdotc;
}
#else
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
#endif	
static inline void cdotu_(complex *z, integer *n_, complex *x, integer *incx_, complex *y, integer *incy_) {
	integer n = *n_, incx = *incx_, incy = *incy_, i;
#ifdef _MSC_VER
	_Fcomplex zdotc = {0.0, 0.0};
	if (incx == 1 && incy == 1) {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc._Val[0] += Cf(&x[i])._Val[0] * Cf(&y[i])._Val[0];
			zdotc._Val[1] += Cf(&x[i])._Val[1] * Cf(&y[i])._Val[1];
		}
	} else {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc._Val[0] += Cf(&x[i*incx])._Val[0] * Cf(&y[i*incy])._Val[0];
			zdotc._Val[1] += Cf(&x[i*incx])._Val[1] * Cf(&y[i*incy])._Val[1];
		}
	}
	pCf(z) = zdotc;
}
#else
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
#endif
static inline void zdotu_(doublecomplex *z, integer *n_, doublecomplex *x, integer *incx_, doublecomplex *y, integer *incy_) {
	integer n = *n_, incx = *incx_, incy = *incy_, i;
#ifdef _MSC_VER
	_Dcomplex zdotc = {0.0, 0.0};
	if (incx == 1 && incy == 1) {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc._Val[0] += Cd(&x[i])._Val[0] * Cd(&y[i])._Val[0];
			zdotc._Val[1] += Cd(&x[i])._Val[1] * Cd(&y[i])._Val[1];
		}
	} else {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc._Val[0] += Cd(&x[i*incx])._Val[0] * Cd(&y[i*incy])._Val[0];
			zdotc._Val[1] += Cd(&x[i*incx])._Val[1] * Cd(&y[i*incy])._Val[1];
		}
	}
	pCd(z) = zdotc;
}
#else
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
#endif
/*  -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/



/*  -- translated by f2c (version 20200916).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/



/* Table of constant values */

static complex c_b1 = {0.f,0.f};
static complex c_b2 = {1.f,0.f};
static integer c__12 = 12;
static integer c__13 = 13;
static integer c__14 = 14;
static integer c__15 = 15;
static integer c__17 = 17;
static integer c_n1 = -1;
static integer c__1 = 1;

/* Subroutine */ void claqz0_(char *wants, char *wantq, char *wantz, integer *
	n, integer *ilo, integer *ihi, complex *a, integer *lda, complex *b, 
	integer *ldb, complex *alpha, complex *beta, complex *q, integer *ldq,
	 complex *z__, integer *ldz, complex *work, integer *lwork, real *
	rwork, integer *rec, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2;
    complex q__1, q__2;

    /* Local variables */
    integer aed_info__, shiftpos, lworkreq, k;
    real c1;
    integer k2;
    complex s1;
    integer norm_info__, ld, ns, n_deflated__, nw, sweep_info__, nbr;
    logical ilq, ilz;
    real ulp;
    integer nsr, nwr;
    real btol;
    integer nmin;
    complex temp;
    extern /* Subroutine */ void crot_(integer *, complex *, integer *, 
	    complex *, integer *, real *, complex *);
    integer n_undeflated__;
    extern logical lsame_(char *, char *);
    integer iiter;
    real bnorm;
    integer maxit, rcost, istop;
    extern /* Subroutine */ void claqz2_(logical *, logical *, logical *, 
	    integer *, integer *, integer *, integer *, complex *, integer *, 
	    complex *, integer *, complex *, integer *, complex *, integer *, 
	    integer *, integer *, complex *, complex *, complex *, integer *, 
	    complex *, integer *, complex *, integer *, real *, integer *, 
	    integer *), claqz3_(logical *, logical *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, complex *, complex *, 
	    complex *, integer *, complex *, integer *, complex *, integer *, 
	    complex *, integer *, complex *, integer *, complex *, integer *, 
	    complex *, integer *, integer *);
    integer itemp1, itemp2, nibble;
    extern real slamch_(char *);
    integer nblock;
    extern real clanhs_(char *, integer *, complex *, integer *, real *);
    extern /* Subroutine */ void claset_(char *, integer *, integer *, complex 
	    *, complex *, complex *, integer *), clartg_(complex *, 
	    complex *, real *, complex *, complex *);
    real safmin;
    extern /* Subroutine */ void xerbla_(char *, integer *, ftnlen);
    real safmax;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ void chgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, complex *, integer *, complex *, integer *, 
	    complex *, complex *, complex *, integer *, complex *, integer *, 
	    complex *, integer *, real *, integer *);
    complex eshift;
    char jbcmpz[3];
    integer iwantq, iwants, istart;
    real smlnum;
    integer istopm, iwantz;
    integer istart2;
    logical ilschur;
    integer nshifts, istartm;

/*     Arguments */
/*     Parameters */
/*     Local scalars */
/*     External Functions */

/*     Decode wantS,wantQ,wantZ */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --alpha;
    --beta;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --rwork;

    /* Function Body */
    if (lsame_(wants, "E")) {
	ilschur = FALSE_;
	iwants = 1;
    } else if (lsame_(wants, "S")) {
	ilschur = TRUE_;
	iwants = 2;
    } else {
	iwants = 0;
    }
    if (lsame_(wantq, "N")) {
	ilq = FALSE_;
	iwantq = 1;
    } else if (lsame_(wantq, "V")) {
	ilq = TRUE_;
	iwantq = 2;
    } else if (lsame_(wantq, "I")) {
	ilq = TRUE_;
	iwantq = 3;
    } else {
	iwantq = 0;
    }
    if (lsame_(wantz, "N")) {
	ilz = FALSE_;
	iwantz = 1;
    } else if (lsame_(wantz, "V")) {
	ilz = TRUE_;
	iwantz = 2;
    } else if (lsame_(wantz, "I")) {
	ilz = TRUE_;
	iwantz = 3;
    } else {
	iwantz = 0;
    }

/*     Check Argument Values */

    *info = 0;
    if (iwants == 0) {
	*info = -1;
    } else if (iwantq == 0) {
	*info = -2;
    } else if (iwantz == 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*ilo < 1) {
	*info = -5;
    } else if (*ihi > *n || *ihi < *ilo - 1) {
	*info = -6;
    } else if (*lda < *n) {
	*info = -8;
    } else if (*ldb < *n) {
	*info = -10;
    } else if (*ldq < 1 || (ilq && *ldq < *n)) {
	*info = -15;
    } else if (*ldz < 1 || (ilz && *ldz < *n)) {
	*info = -17;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CLAQZ0", &i__1, (ftnlen)6);
	return;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	work[1].r = 1.f, work[1].i = 0.f;
	return;
    }

/*     Get the parameters */

    *(unsigned char *)jbcmpz = *(unsigned char *)wants;
    *(unsigned char *)&jbcmpz[1] = *(unsigned char *)wantq;
    *(unsigned char *)&jbcmpz[2] = *(unsigned char *)wantz;
    nmin = ilaenv_(&c__12, "CLAQZ0", jbcmpz, n, ilo, ihi, lwork, (ftnlen)6, (
	    ftnlen)3);
    nwr = ilaenv_(&c__13, "CLAQZ0", jbcmpz, n, ilo, ihi, lwork, (ftnlen)6, (
	    ftnlen)3);
    nwr = f2cmax(2,nwr);
/* Computing MIN */
    i__1 = *ihi - *ilo + 1, i__2 = (*n - 1) / 3, i__1 = f2cmin(i__1,i__2);
    nwr = f2cmin(i__1,nwr);
    nibble = ilaenv_(&c__14, "CLAQZ0", jbcmpz, n, ilo, ihi, lwork, (ftnlen)6, 
	    (ftnlen)3);
    nsr = ilaenv_(&c__15, "CLAQZ0", jbcmpz, n, ilo, ihi, lwork, (ftnlen)6, (
	    ftnlen)3);
/* Computing MIN */
    i__1 = nsr, i__2 = (*n + 6) / 9, i__1 = f2cmin(i__1,i__2), i__2 = *ihi - *
	    ilo;
    nsr = f2cmin(i__1,i__2);
/* Computing MAX */
    i__1 = 2, i__2 = nsr - nsr % 2;
    nsr = f2cmax(i__1,i__2);
    rcost = ilaenv_(&c__17, "CLAQZ0", jbcmpz, n, ilo, ihi, lwork, (ftnlen)6, (
	    ftnlen)3);
    itemp1 = (integer) ((real) nsr / sqrt((real) nsr * 2 / ((real) rcost / 
	    100 * (real) (*n)) + 1));
    itemp1 = ((itemp1 - 1) / 4 << 2) + 4;
    nbr = nsr + itemp1;
    if (*n < nmin || *rec >= 2) {
	chgeqz_(wants, wantq, wantz, n, ilo, ihi, &a[a_offset], lda, &b[
		b_offset], ldb, &alpha[1], &beta[1], &q[q_offset], ldq, &z__[
		z_offset], ldz, &work[1], lwork, &rwork[1], info);
	return;
    }

/*     Find out required workspace */

/*     Workspace query to CLAQZ2 */
    nw = f2cmax(nwr,nmin);
    claqz2_(&ilschur, &ilq, &ilz, n, ilo, ihi, &nw, &a[a_offset], lda, &b[
	    b_offset], ldb, &q[q_offset], ldq, &z__[z_offset], ldz, &
	    n_undeflated__, &n_deflated__, &alpha[1], &beta[1], &work[1], &nw,
	     &work[1], &nw, &work[1], &c_n1, &rwork[1], rec, &aed_info__);
    itemp1 = (integer) work[1].r;
/*     Workspace query to CLAQZ3 */
    claqz3_(&ilschur, &ilq, &ilz, n, ilo, ihi, &nsr, &nbr, &alpha[1], &beta[1]
	    , &a[a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq, &z__[
	    z_offset], ldz, &work[1], &nbr, &work[1], &nbr, &work[1], &c_n1, &
	    sweep_info__);
    itemp2 = (integer) work[1].r;
/* Computing MAX */
/* Computing 2nd power */
    i__3 = nw;
/* Computing 2nd power */
    i__4 = nbr;
    i__1 = itemp1 + (i__3 * i__3 << 1), i__2 = itemp2 + (i__4 * i__4 << 1);
    lworkreq = f2cmax(i__1,i__2);
    if (*lwork == -1) {
	r__1 = (real) lworkreq;
	work[1].r = r__1, work[1].i = 0.f;
	return;
    } else if (*lwork < lworkreq) {
	*info = -18;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CLAQZ0", &i__1, (ftnlen)6);
	return;
    }

/*     Initialize Q and Z */

    if (iwantq == 3) {
	claset_("FULL", n, n, &c_b1, &c_b2, &q[q_offset], ldq);
    }
    if (iwantz == 3) {
	claset_("FULL", n, n, &c_b1, &c_b2, &z__[z_offset], ldz);
    }
/*     Get machine constants */
    safmin = slamch_("SAFE MINIMUM");
    safmax = 1.f / safmin;
    ulp = slamch_("PRECISION");
    smlnum = safmin * ((real) (*n) / ulp);
    i__1 = *ihi - *ilo + 1;
    bnorm = clanhs_("F", &i__1, &b[*ilo + *ilo * b_dim1], ldb, &rwork[1]);
/* Computing MAX */
    r__1 = safmin, r__2 = ulp * bnorm;
    btol = f2cmax(r__1,r__2);
    istart = *ilo;
    istop = *ihi;
    maxit = (*ihi - *ilo + 1) * 30;
    ld = 0;
    i__1 = maxit;
    for (iiter = 1; iiter <= i__1; ++iiter) {
	if (iiter >= maxit) {
	    *info = istop + 1;
	    goto L80;
	}
	if (istart + 1 >= istop) {
	    istop = istart;
	    myexit_();
	}
/*        Check deflations at the end */
/* Computing MAX */
	r__1 = smlnum, r__2 = ulp * (c_abs(&a[istop + istop * a_dim1]) + 
		c_abs(&a[istop - 1 + (istop - 1) * a_dim1]));
	if (c_abs(&a[istop + (istop - 1) * a_dim1]) <= f2cmax(r__1,r__2)) {
	    i__2 = istop + (istop - 1) * a_dim1;
	    a[i__2].r = 0.f, a[i__2].i = 0.f;
	    --istop;
	    ld = 0;
	    eshift.r = 0.f, eshift.i = 0.f;
	}
/*        Check deflations at the start */
/* Computing MAX */
	r__1 = smlnum, r__2 = ulp * (c_abs(&a[istart + istart * a_dim1]) + 
		c_abs(&a[istart + 1 + (istart + 1) * a_dim1]));
	if (c_abs(&a[istart + 1 + istart * a_dim1]) <= f2cmax(r__1,r__2)) {
	    i__2 = istart + 1 + istart * a_dim1;
	    a[i__2].r = 0.f, a[i__2].i = 0.f;
	    ++istart;
	    ld = 0;
	    eshift.r = 0.f, eshift.i = 0.f;
	}
	if (istart + 1 >= istop) {
	    myexit_();
	}
/*        Check interior deflations */
	istart2 = istart;
	i__2 = istart + 1;
	for (k = istop; k >= i__2; --k) {
/* Computing MAX */
	    r__1 = smlnum, r__2 = ulp * (c_abs(&a[k + k * a_dim1]) + c_abs(&a[
		    k - 1 + (k - 1) * a_dim1]));
	    if (c_abs(&a[k + (k - 1) * a_dim1]) <= f2cmax(r__1,r__2)) {
		i__3 = k + (k - 1) * a_dim1;
		a[i__3].r = 0.f, a[i__3].i = 0.f;
		istart2 = k;
		myexit_();
	    }
	}
/*        Get range to apply rotations to */
	if (ilschur) {
	    istartm = 1;
	    istopm = *n;
	} else {
	    istartm = istart2;
	    istopm = istop;
	}
/*        Check infinite eigenvalues, this is done without blocking so might */
/*        slow down the method when many infinite eigenvalues are present */
	k = istop;
	while(k >= istart2) {
	    if (c_abs(&b[k + k * b_dim1]) < btol) {
/*              A diagonal element of B is negligible, move it */
/*              to the top and deflate it */
		i__2 = istart2 + 1;
		for (k2 = k; k2 >= i__2; --k2) {
		    clartg_(&b[k2 - 1 + k2 * b_dim1], &b[k2 - 1 + (k2 - 1) * 
			    b_dim1], &c1, &s1, &temp);
		    i__3 = k2 - 1 + k2 * b_dim1;
		    b[i__3].r = temp.r, b[i__3].i = temp.i;
		    i__3 = k2 - 1 + (k2 - 1) * b_dim1;
		    b[i__3].r = 0.f, b[i__3].i = 0.f;
		    i__3 = k2 - 2 - istartm + 1;
		    crot_(&i__3, &b[istartm + k2 * b_dim1], &c__1, &b[istartm 
			    + (k2 - 1) * b_dim1], &c__1, &c1, &s1);
/* Computing MIN */
		    i__4 = k2 + 1;
		    i__3 = f2cmin(i__4,istop) - istartm + 1;
		    crot_(&i__3, &a[istartm + k2 * a_dim1], &c__1, &a[istartm 
			    + (k2 - 1) * a_dim1], &c__1, &c1, &s1);
		    if (ilz) {
			crot_(n, &z__[k2 * z_dim1 + 1], &c__1, &z__[(k2 - 1) *
				 z_dim1 + 1], &c__1, &c1, &s1);
		    }
		    if (k2 < istop) {
			clartg_(&a[k2 + (k2 - 1) * a_dim1], &a[k2 + 1 + (k2 - 
				1) * a_dim1], &c1, &s1, &temp);
			i__3 = k2 + (k2 - 1) * a_dim1;
			a[i__3].r = temp.r, a[i__3].i = temp.i;
			i__3 = k2 + 1 + (k2 - 1) * a_dim1;
			a[i__3].r = 0.f, a[i__3].i = 0.f;
			i__3 = istopm - k2 + 1;
			crot_(&i__3, &a[k2 + k2 * a_dim1], lda, &a[k2 + 1 + 
				k2 * a_dim1], lda, &c1, &s1);
			i__3 = istopm - k2 + 1;
			crot_(&i__3, &b[k2 + k2 * b_dim1], ldb, &b[k2 + 1 + 
				k2 * b_dim1], ldb, &c1, &s1);
			if (ilq) {
			    r_cnjg(&q__1, &s1);
			    crot_(n, &q[k2 * q_dim1 + 1], &c__1, &q[(k2 + 1) *
				     q_dim1 + 1], &c__1, &c1, &q__1);
			}
		    }
		}
		if (istart2 < istop) {
		    clartg_(&a[istart2 + istart2 * a_dim1], &a[istart2 + 1 + 
			    istart2 * a_dim1], &c1, &s1, &temp);
		    i__2 = istart2 + istart2 * a_dim1;
		    a[i__2].r = temp.r, a[i__2].i = temp.i;
		    i__2 = istart2 + 1 + istart2 * a_dim1;
		    a[i__2].r = 0.f, a[i__2].i = 0.f;
		    i__2 = istopm - (istart2 + 1) + 1;
		    crot_(&i__2, &a[istart2 + (istart2 + 1) * a_dim1], lda, &
			    a[istart2 + 1 + (istart2 + 1) * a_dim1], lda, &c1,
			     &s1);
		    i__2 = istopm - (istart2 + 1) + 1;
		    crot_(&i__2, &b[istart2 + (istart2 + 1) * b_dim1], ldb, &
			    b[istart2 + 1 + (istart2 + 1) * b_dim1], ldb, &c1,
			     &s1);
		    if (ilq) {
			r_cnjg(&q__1, &s1);
			crot_(n, &q[istart2 * q_dim1 + 1], &c__1, &q[(istart2 
				+ 1) * q_dim1 + 1], &c__1, &c1, &q__1);
		    }
		}
		++istart2;
	    }
	    --k;
	}
/*        istart2 now points to the top of the bottom right */
/*        unreduced Hessenberg block */
	if (istart2 >= istop) {
	    istop = istart2 - 1;
	    ld = 0;
	    eshift.r = 0.f, eshift.i = 0.f;
	    mycycle_();
	}
	nw = nwr;
	nshifts = nsr;
	nblock = nbr;
	if (istop - istart2 + 1 < nmin) {
/*           Setting nw to the size of the subblock will make AED deflate */
/*           all the eigenvalues. This is slightly more efficient than just */
/*           using CHGEQZ because the off diagonal part gets updated via BLAS. */
	    if (istop - istart + 1 < nmin) {
		nw = istop - istart + 1;
		istart2 = istart;
	    } else {
		nw = istop - istart2 + 1;
	    }
	}

/*        Time for AED */

/* Computing 2nd power */
	i__2 = nw;
/* Computing 2nd power */
	i__3 = nw;
/* Computing 2nd power */
	i__5 = nw;
	i__4 = *lwork - (i__5 * i__5 << 1);
	claqz2_(&ilschur, &ilq, &ilz, n, &istart2, &istop, &nw, &a[a_offset], 
		lda, &b[b_offset], ldb, &q[q_offset], ldq, &z__[z_offset], 
		ldz, &n_undeflated__, &n_deflated__, &alpha[1], &beta[1], &
		work[1], &nw, &work[i__2 * i__2 + 1], &nw, &work[(i__3 * i__3 
		<< 1) + 1], &i__4, &rwork[1], rec, &aed_info__);
	if (n_deflated__ > 0) {
	    istop -= n_deflated__;
	    ld = 0;
	    eshift.r = 0.f, eshift.i = 0.f;
	}
	if (n_deflated__ * 100 > nibble * (n_deflated__ + n_undeflated__) || 
		istop - istart2 + 1 < nmin) {
/*           AED has uncovered many eigenvalues. Skip a QZ sweep and run */
/*           AED again. */
	    mycycle_();
	}
	++ld;
/* Computing MIN */
	i__2 = nshifts, i__3 = istop - istart2;
	ns = f2cmin(i__2,i__3);
	ns = f2cmin(ns,n_undeflated__);
	shiftpos = istop - n_undeflated__ + 1;
	if (ld % 6 == 0) {

/*           Exceptional shift.  Chosen for no particularly good reason. */

	    if ((real) maxit * safmin * c_abs(&a[istop + (istop - 1) * a_dim1]
		    ) < c_abs(&a[istop - 1 + (istop - 1) * a_dim1])) {
		c_div(&q__1, &a[istop + (istop - 1) * a_dim1], &b[istop - 1 + 
			(istop - 1) * b_dim1]);
		eshift.r = q__1.r, eshift.i = q__1.i;
	    } else {
		r__1 = safmin * (real) maxit;
		q__2.r = 1.f / r__1, q__2.i = 0.f / r__1;
		q__1.r = eshift.r + q__2.r, q__1.i = eshift.i + q__2.i;
		eshift.r = q__1.r, eshift.i = q__1.i;
	    }
	    i__2 = shiftpos;
	    alpha[i__2].r = 1.f, alpha[i__2].i = 0.f;
	    i__2 = shiftpos;
	    beta[i__2].r = eshift.r, beta[i__2].i = eshift.i;
	    ns = 1;
	}

/*        Time for a QZ sweep */

/* Computing 2nd power */
	i__2 = nblock;
/* Computing 2nd power */
	i__3 = nblock;
/* Computing 2nd power */
	i__5 = nblock;
	i__4 = *lwork - (i__5 * i__5 << 1);
	claqz3_(&ilschur, &ilq, &ilz, n, &istart2, &istop, &ns, &nblock, &
		alpha[shiftpos], &beta[shiftpos], &a[a_offset], lda, &b[
		b_offset], ldb, &q[q_offset], ldq, &z__[z_offset], ldz, &work[
		1], &nblock, &work[i__2 * i__2 + 1], &nblock, &work[(i__3 * 
		i__3 << 1) + 1], &i__4, &sweep_info__);
    }

/*     Call CHGEQZ to normalize the eigenvalue blocks and set the eigenvalues */
/*     If all the eigenvalues have been found, CHGEQZ will not do any iterations */
/*     and only normalize the blocks. In case of a rare convergence failure, */
/*     the single shift might perform better. */

L80:
    chgeqz_(wants, wantq, wantz, n, ilo, ihi, &a[a_offset], lda, &b[b_offset],
	     ldb, &alpha[1], &beta[1], &q[q_offset], ldq, &z__[z_offset], ldz,
	     &work[1], lwork, &rwork[1], &norm_info__);
    *info = norm_info__;
    return;
} /* claqz0_ */

