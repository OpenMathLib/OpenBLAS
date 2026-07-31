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
static inline _Dcomplex Cd(doublecomplex *z) {_Dcomplex zz={z->r , z->i};return zz;}
static inline _Fcomplex * _pCf(complex *z) {return (_Fcomplex*)z;}
static inline _Dcomplex * _pCd(doublecomplex *z) {return (_Dcomplex*)z;}
#else
static inline _Complex float Cf(complex *z) {return z->r + z->i*_Complex_I;}
static inline _Complex double Cd(doublecomplex *z) {return z->r + z->i*_Complex_I;}
static inline _Complex float * _pCf(complex *z) {return (_Complex float*)z;}
static inline _Complex double * _pCd(doublecomplex *z) {return (_Complex double*)z;}
#endif
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
static char junk[] = "\n@(#)LIBF77 VERSION 19990503\n";
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

static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b15 = -1.;
static integer c__1 = 1;

/* Subroutine */ int zgecxx_(char *fact, char *usesd, integer *m, integer *n, 
	integer *desel_rows__, integer *sel_desel_cols__, integer *kmaxfree, 
	doublereal *abstol, doublereal *reltol, doublecomplex *a, integer *
	lda, integer *k, doublereal *maxc2nrmk, doublereal *relmaxc2nrmk, 
	doublereal *fnrmk, integer *ipiv, integer *jpiv, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *qrc, integer *ldqrc, 
	doublecomplex *x, integer *ldx, doublecomplex *work, integer *lwork, 
	doublereal *rwork, integer *lrwork, integer *iwork, integer *liwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, qrc_dim1, qrc_offset, x_dim1, 
	    x_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Local variables */
    doublereal maxc2nrm;
    extern /* Subroutine */ int zgeqp3rk_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, integer *,
	     integer *, doublereal *, doublereal *, integer *, doublecomplex *
	    , doublecomplex *, integer *, doublereal *, integer *, integer *);
    doublereal relmaxc2nrmkfree;
    integer i__, j, minmnfree, ip, jp;
    doublereal abstolfree;
    integer kp0;
    doublereal reltolfree;
    logical use_sel_desel_cols__;
    doublereal eps;
    integer nsel, msub, nsub, kfree, mfree, nfree;
    extern logical lsame_(char *, char *);
    doublereal maxc2nrmkfree;
    integer iinfo, itemp, minmn;
    extern /* Subroutine */ int zgels_(char *, integer *, integer *, integer *
	    , doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *), zcopy_(integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), zswap_(
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *)
	    ;
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *), dlamch_(
	    char *);
    integer jdesel, mdesel, ndesel;
    extern logical disnan_(doublereal *);
    doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *);
    integer mresid, nresid;
    extern integer izamax_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int zgeqrf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    );
    integer kmaxls, lwkmin;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    logical usetol, lquery;
    integer lwkopt;
    extern /* Subroutine */ int zunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *);
    logical use_desel_rows__;
    integer liwkmin;
    logical returnc;
    integer lrwkmin, liwkopt, lrwkopt;
    logical returnx;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */


/*  ===================================================================== */


/*     Test the input arguments */

    /* Parameter adjustments */
    --desel_rows__;
    --sel_desel_cols__;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    --jpiv;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    qrc_dim1 = *ldqrc;
    qrc_offset = 1 + qrc_dim1;
    qrc -= qrc_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --work;
    --rwork;
    --iwork;

    /* Function Body */
    *info = 0;
    mdesel = 0;
    nsel = 0;
    ndesel = 0;
    msub = *m;
    nsub = *n;
    mfree = msub;
    nfree = nsub;
    minmn = f2cmin(*m,*n);

    lquery = *lwork == -1 || *lrwork == -1 || *liwork == -1;

    returnx = lsame_(fact, "X");
    returnc = lsame_(fact, "C") || returnx;

    use_desel_rows__ = lsame_(usesd, "R") || lsame_(
	    usesd, "A");
    use_sel_desel_cols__ = lsame_(usesd, "C") || lsame_(
	    usesd, "A");

    if (! (returnc || lsame_(fact, "P"))) {
	*info = -1;
    } else if (! (use_desel_rows__ || use_sel_desel_cols__ || lsame_(usesd, 
	    "N"))) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else {

/*        This is to check that the number of preselected columns NSEL */
/*        cannot be larger than MSUB, which is the number of rows */
/*        without MDESEL deselected rows. When the number of */
/*        preselected columns NSEL is larger than MSUB, */
/*        the factorization of all preselected NSEL columns cannot be */
/*        completed. MSUB also will be used for LDX argument check */
/*        later. */

	if (use_desel_rows__) {

/*           Count the number of free rows MSUB. */

	    i__1 = *m;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (desel_rows__[i__] == -1) {
		    ++mdesel;
		}
	    }
	    msub = *m - mdesel;
	    mfree = msub;
	}

	if (use_sel_desel_cols__) {

/*           Count the number of preselected columns NSEL and the */
/*           number of preselected and free columns NSUB = N - NDESEL. */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (sel_desel_cols__[j] == 1) {
		    ++nsel;
		}
		if (sel_desel_cols__[j] == -1) {
		    ++ndesel;
		}
	    }
	    nsub = *n - ndesel;
	    mfree = msub - nsel;
	    nfree = nsub - nsel;

	}
	minmnfree = f2cmin(mfree,nfree);

	if (nsel > msub) {
	    *info = -6;
	} else if (*kmaxfree < 0) {
	    *info = -7;
	} else if (disnan_(abstol)) {
	    *info = -8;
	} else if (disnan_(reltol)) {
	    *info = -9;
	} else if (*lda < f2cmax(1,*m)) {
	    *info = -11;
/*        This is a check for LDC */
	} else if (returnc && *ldc < f2cmax(1,*m) || ! returnc && *ldc < 1) {
	    *info = -20;
/*        This is a check for LDQRC */
	} else if (returnx && *ldqrc < f2cmax(1,*m) || ! returnx && *ldqrc < 1) {
	    *info = -22;
/*        This is a check for LDX */
	} else if (returnx && *ldx < f2cmax(1,*m) || ! returnx && *ldx < 1) {
	    *info = -24;
	}

    }

/*     ================================================================== */

/*       a) Test the input workspace size LWORK, LRWORK, LIWORK for the */
/*          minimum size requirement LWKMIN, LRWKMIN, LIWKMIN */
/*          respectively. */
/*       b) Determine the optimal workspace sizes LWKOPT, LRWKOPT, */
/*          and LIWKOPT to be returned in */
/*          WORK( 1 ), RWORK( 1 ) and IWORK( 1 ) respectively, */
/*          if INFO >= 0 in cases: */
/*           (1) LQUERY = .TRUE., */
/*           (2) when the routine exits. */
/*     Here, LWKMIN, LRWKMIN and LIWKMIN are the minimum workspaces */
/*     required for unblocked code. */

    if (*info == 0) {
	if (minmn == 0) {
	    lwkmin = 1;
	    lwkopt = 1;
	    lrwkmin = 1;
	    lrwkopt = 1;
	    liwkmin = 1;
	    liwkopt = 1;
	} else {

/*           (Complex_wk_part_1) Complex minimum and optimal workspace */
/*           computation. */

	    lwkmin = 1;
	    lwkopt = lwkmin;

/*           (Real_wk_part_1) Real minimum workspace computation. */
/*           LRWKMIN = MAX(1, NSUB) for column 2-norm computation */

	    lrwkmin = f2cmax(1,nsub);

/*           (Int_wk_part_1) Integer minimum workspace computation. */

	    liwkmin = 1;

/*           Call of ZGEQRF. */

	    if (nsel > 0) {

/*              (Complex_wk_part_2) Complex minimum workspace */
/*              computation. */

		lwkmin = f2cmax(lwkmin,nsel);

/*              Query for optimal workspace size for ZGEQRF. */

		zgeqrf_(&msub, &nsel, &a[a_offset], lda, &tau[1], &work[1], &
			c_n1, &iinfo);
/* Computing MAX */
		i__1 = lwkopt, i__2 = (integer) work[1].r;
		lwkopt = f2cmax(i__1,i__2);

/*              Call of ZUNMQR. */

		if (nfree > 0) {

/*                 (Complex_wk_part_3) Complex minimum workspace */
/*                 computation. */

		    lwkmin = f2cmax(lwkmin,nfree);

/*                 Query for optimal workspace size for ZUNMQR. */

		    zunmqr_("L", "C", &msub, &nfree, &nsel, &a[a_offset], lda,
			     &tau[1], &a[(nsel + 1) * a_dim1 + 1], lda, &work[
			    1], &c_n1, &iinfo);
/* Computing MAX */
		    i__1 = lwkopt, i__2 = (integer) work[1].r;
		    lwkopt = f2cmax(i__1,i__2);
		}

	    }

/*           Call of ZGEQP3RK. */

	    if (minmnfree != 0) {

/*              (Complex_wk_part_4) Complex minimum workspace */
/*              computation. */
/*              LWKMIN = MAX(1, NFREE-1) for the call of ZGEQP3RK. */

/* Computing MAX */
		i__1 = lwkmin, i__2 = nfree - 1;
		lwkmin = f2cmax(i__1,i__2);

/*              Query for optimal workspace size for ZGEQP3RK. */

		zgeqp3rk_(&mfree, &nfree, &c__0, &nfree, &c_b15, &c_b15, &a[
			a_dim1 + 1], lda, &kfree, &maxc2nrmkfree, &
			relmaxc2nrmkfree, &jpiv[1], &tau[1], &work[1], &c_n1, 
			&rwork[1], &iwork[1], &iinfo);
/* Computing MAX */
		i__1 = lwkopt, i__2 = (integer) work[1].r;
		lwkopt = f2cmax(i__1,i__2);

/*              (Real_wk_part_2) Real minimum workspace computation. */
/*              LRWKMIN = MAX(1, 2*NFREE) for the call of ZGEQP3RK. */

/* Computing MAX */
		i__1 = lrwkmin, i__2 = nfree << 1;
		lrwkmin = f2cmax(i__1,i__2);

/*              (Int_wk_part_2) Integer minimum workspace computation. */
/*              LIWKMIN =  NFREE-1 for the call of ZGEQP3RK. */

/* Computing MAX */
		i__1 = liwkmin, i__2 = nfree - 1;
		liwkmin = f2cmax(i__1,i__2);

		if (nsel != 0) {

/*                 (Int_wk_part_3) Integer minimum workspace computation. */
/*                 NFREE is for ZGEQP3RK and NFREE-1 for JPIV adjustment. */

/* Computing MAX */
		    i__1 = liwkmin, i__2 = nfree + nfree - 1;
		    liwkmin = f2cmax(i__1,i__2);
		}

	    }

	    if (returnc) {

/*              Integer minimum workspace computation. */
/*              (Int_wk_part_4) LIWKMIN = 2*N for applying the */
/*              interchanges for the columns in the matrix C. */

/* Computing MAX */
		i__1 = liwkmin, i__2 = *n << 1;
		liwkmin = f2cmax(i__1,i__2);
	    }

/*           Real and Integer optimal workspace computation. */

	    lrwkopt = lrwkmin;
	    liwkopt = liwkmin;

/*           Call of ZGELS. */

	    if (returnx) {

/*              (Complex_wk_part_5) Complex minimum workspace computation. */
/*              LWKMIN = f2cmax( 1, MINMN + f2cmax( MINMN, N ) ) = */
/*                     = f2cmax( 1, MINMN + N ) for the call of ZGELS. */

/* Computing MAX */
		i__1 = lwkmin, i__2 = minmn + *n;
		lwkmin = f2cmax(i__1,i__2);

/*              Query for optimal workspace size for ZGELS. */

		kmaxls = minmn;

		zgels_("N", m, &kmaxls, n, &qrc[qrc_offset], ldqrc, &x[
			x_offset], ldx, &work[1], &c_n1, &iinfo);
/* Computing MAX */
		i__1 = lwkopt, i__2 = (integer) work[1].r;
		lwkopt = f2cmax(i__1,i__2);

	    }

/*           End of ELSE for IF( MINMN.EQ.0 ) */

	}

	if (*lwork < lwkmin && ! lquery) {
	    *info = -26;
	} else if (*lrwork < lrwkmin && ! lquery) {
	    *info = -28;
	} else if (*liwork < liwkmin && ! lquery) {
	    *info = -30;
	}
    }

    if (*info == 0) {
	z__1.r = (doublereal) lwkopt, z__1.i = 0.;
	work[1].r = z__1.r, work[1].i = z__1.i;
	rwork[1] = (doublereal) lrwkopt;
	iwork[1] = liwkopt;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZGECXX", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     ================================================================== */

/*     Quick return if possible for: */
/*     a)  M = 0 or N = 0. There is no matrix A(1:M,1:N). */
/*     b)  MSUB = 0 or NSUB = 0. There is no matrix A_sub(1:MSUB,1:NSUB). */
/*     NOTE: f2cmin( M, N) = 0 implies f2cmin( MSUB, NSUB) = 0. */
/*     We need to return correct values for all scalar output parameters, */
/*     (including WORK(1) and IWORK(1), which are set above). */

    if (f2cmin(msub,nsub) == 0) {
	*k = 0;
	*maxc2nrmk = 0.;
	*relmaxc2nrmk = 0.;
	*fnrmk = 0.;
	return 0;
    }

/*     ================================================================== */

    *k = 0;

/*     If we need to return factor X, copy the original untouched matrix */
/*     A into the array X. */

    if (returnx) {
	zlacpy_("F", m, n, &a[a_offset], lda, &x[x_offset], ldx);
    }

/*     If we need to return the factor C, copy the original matrix A */
/*     into the array C, only if do not return the factor X. In this */
/*     case, we need to choose the columns of the matrix A in the array C */
/*     in place, otherwise we can copy the columns of the matrix A from */
/*     the array X. */

    if (returnc && ! returnx) {
	zlacpy_("F", m, n, &a[a_offset], lda, &c__[c_offset], ldc);
    }

/*     ================================================================== */
/*     Permute the deselected rows to the bottom of the matrix A. */
/*     1) The initial order of included rows in their block is preserved. */
/*     2) The initial order of deselected rows in their block is not */
/*        preserved. */
/*     ================================================================== */

/*     I is an index of DESEL_ROWS array and a row index of */
/*     the matrix A. MSUB is the number of processed included rows, which */
/*     is also an index pointer to the last included row in the matrix A. */
/*     We can think of I as a row source index, and MSUB as a destination */
/*     index for moving an included row in the matrix A. */

/*     ( We start with MSUB = 0. We loop over index I in (1:M), and */
/*     for each position I in DESEL_ROWS  array, we check if the row at */
/*     the position I in the matrix A is an included row (not -1 value). */
/*     If it is an included row, we increment MSUB pointer, otherwise */
/*     we do not change MSUB index pointer. Then, we bring this included */
/*     row from the index I in the matrix A into smaller (or same) */
/*     MSUB index in the matrix A.  If I = MSUB, then the included row */
/*     is already in place. Due to row swap, the deselected row */
/*     at MSUB index will move into I index in the matrix A. In this way, */
/*     we move all the included rows to the top matrix block preserving */
/*     their initial order within the included block. The initial order */
/*     of deselected rows will not be preserved within their block. */

    if (use_desel_rows__) {

	msub = 0;
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Initialize the row pivot array IPIV. */
	    ipiv[i__] = i__;

/*           The row at the index I is an included row and should be */
/*           moved to the top of the matrix A. */

	    if (desel_rows__[i__] != -1) {
		++msub;

/*              This is a check whether the included row is */
/*              on the included place already. */

		if (i__ != msub) {

/*                 Here, we swap A(I,1:N) into A(MSUB,1:N). */

		    zswap_(n, &a[i__ + a_dim1], lda, &a[msub + a_dim1], lda);

/*                 Save the interchange. */

		    ipiv[i__] = ipiv[msub];
		    ipiv[msub] = i__;
		    desel_rows__[msub] = desel_rows__[i__];
		    desel_rows__[i__] = -1;
		}
	    }

	}

    } else {

/*        We do not use the row deselection DESEL_ROWS array. */
/*        Initialize the row pivot array IPIV. */
/*        NOTE: MSUB=M has default value, */
/*        which is set at the beginning of the routine, before argument */
/*        checks. */

	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ipiv[i__] = i__;
	}
    }

/*     ================================================================== */
/*     Permute the preselected columns to the left and deselected */
/*     columns to the right of the matrix A. */
/*     1) The order of preselected columns is preserved. */
/*     2) The order of free columns is not preserved. */
/*     3) The order of deselected columns is not preserved. */
/*     ================================================================== */

/*     J is the index of SEL_DESEL_COLS array and column J */
/*     of the matrix A. */

    if (use_sel_desel_cols__) {

/*        Column selection. */
/*        NSEL is the number of selected columns, also the pointer to */
/*        the last selected column. */

	nsel = 0;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {

/*           Initialize column pivot array JPIV. */
	    jpiv[j] = j;

	    if (sel_desel_cols__[j] == 1) {
		++nsel;

/*              This is the check whether the selected column is */
/*              on the selected place already. */

		if (j != nsel) {

/*                 Here, we swap the column A(1:M,J) into A(1:M,NSEL) */

		    zswap_(m, &a[j * a_dim1 + 1], &c__1, &a[nsel * a_dim1 + 1]
			    , &c__1);
		    jpiv[j] = jpiv[nsel];
		    jpiv[nsel] = j;
		    sel_desel_cols__[j] = sel_desel_cols__[nsel];
		    sel_desel_cols__[nsel] = 1;
		}
	    }
	}

/*        Column deselection. */
/*        JDESEL the pointer to the last */
/*        deselected column counting right-to-left. */

	jdesel = *n + 1;
	i__1 = nsel + 1;
	for (j = *n; j >= i__1; --j) {
	    if (sel_desel_cols__[j] == -1) {
		--jdesel;

/*              This is the check whether the deselected column is */
/*              on the deselected place already. */

		if (j != jdesel) {

/*                 Here, we swap the column A(1:M,J) into A(1:M,JDESEL) */

		    zswap_(m, &a[j * a_dim1 + 1], &c__1, &a[jdesel * a_dim1 + 
			    1], &c__1);
		    itemp = jpiv[j];
		    jpiv[j] = jpiv[jdesel];
		    jpiv[jdesel] = itemp;
		    sel_desel_cols__[j] = sel_desel_cols__[jdesel];
		    sel_desel_cols__[jdesel] = -1;
		}
	    }
	}

	nsub = jdesel - 1;

    } else {

/*        We do not use the column selection deselection */
/*        SEL_DESEL_COLS array. */
/*        Initialize column pivot array JPIV. */
/*        NOTE: NSUB=N has default value, */
/*        which is set at the beginning of the routine, before argument */
/*        checks. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    jpiv[j] = j;
	}

    }

/*     ================================================================== */
/*     Compute the complete column 2-norms of the submatrix */
/*     A_sub = A(1:MSUB, 1:NSUB) and store them in WORK(1:NSUB). */

    i__1 = nsub;
    for (j = 1; j <= i__1; ++j) {
	rwork[j] = dznrm2_(&msub, &a[j * a_dim1 + 1], &c__1);
    }

/*     Compute the column index of the maximum column 2-norm and */
/*     the maximum column 2-norm itself for the submatrix */
/*     A_sub = A(1:MSUB, 1:NSUB). */

    kp0 = izamax_(&nsub, &work[1], &c__1);
    maxc2nrm = rwork[kp0];

/*     ================================================================== */
/*     Process preselected columns */

/*     Compute the QR factorization of NSEL preselected columns (1:NSEL) */
/*     in the submatrix A_sub = A(1:MSUB, 1:NSUB) and update */
/*     remaining NFREE free columns (NSEL+1:NSUB). */
/*     NSUB = NSEL + NFREE */

    if (nsel > 0) {

/*           Case (a): MSUB < NSEL. */

/*              This is handled at the argument check stage in the */
/*              beginning of the routine. When the number of preselected */
/*              columns is larger than MSUB, hence the factorization of */
/*              all NSEL columns cannot be completed. Return from the */
/*              routine with the error of COL_SEL_DESEL parameter. */

/*           Case (b): MSUB = NSEL. */
/*           Case (c-1): MSUB > NSEL and NSEL = NSUB. */

/*              For cases (b) and (c-1), there will be no residual */
/*              submatrix  after factorization of NSEL columns */
/*              at step K = NSEL: */
/*              A_sub_resid(NSEL) = A(NSEL+1:MSUB, NSEL+1:NSUB). */

/*           Case (c-2): MSUB > NSEL and NSEL < NSUB. */

/*              For Case (c-2) is a submatrix residual at step K=NSEL */
/*              A_sub_resid(NSEL) = A(NSEL+1:MSUB, NSEL+1:NSUB) */

	zgeqrf_(&msub, &nsel, &a[a_offset], lda, &tau[1], &work[1], lwork, &
		iinfo);

/*        Apply Q**T from the left to A(NSEL+1:MSUB, NSEL+1:NSUB) */

	if (nfree > 0) {

/*           This is only for case (c-2) ('L' = Left, 'T' = Transpose) */

	    zunmqr_("L", "C", &msub, &nfree, &nsel, &a[a_offset], lda, &tau[1]
		    , &a[(nsel + 1) * a_dim1 + 1], lda, &work[1], lwork, &
		    iinfo);
	}

	*k += nsel;

/*        End of IF(NSEL.GT.0) */

    }

/*     ================================================================== */

    kfree = 0;

    if (minmnfree != 0) {

/*        Factorize NFREE free columns of */
/*        A_free = A_sub_resid(NSEL) = A(NSEL+1:MSUB, NSEL+1:NSUB), */
/*        KFREE is the number of columns that were actually factorized */
/*        among NFREE columns. */

/*     ================================================================== */

	eps = dlamch_("Epsilon");

	usetol = FALSE_;

/*        Adjust ABSTOL only if nonnegative. Negative value means disabled. */
/*        We need to keep negative value for later use in criterion */
/*        check. */

	if (*abstol >= 0.) {
	    safmin = dlamch_("Safe minimum");
/* Computing MAX */
	    d__1 = *abstol, d__2 = safmin * 2.;
	    *abstol = f2cmax(d__1,d__2);
	    usetol = TRUE_;
	}

/*        Adjust RELTOL only if nonnegative. Negative value means disabled. */
/*        We need to keep negative value for later use in criterion */
/*        check. */

	if (*reltol >= 0.) {
	    *reltol = f2cmax(*reltol,eps);
	    usetol = TRUE_;
	}

/*     ================================================================== */

/*        Disable RELTOLFREE when calling ZGEQP3RK for free columns */
/*        factorization, since ZGEQP3RK expects RELTOLFREE with respect */
/*        to the residual matrix A_sub_resid(NSEL), not the whole */
/*        original matrix A. We can use RELTOL criterion by passing it */
/*        to ABSTOLFREE as RELTOL*MAXC2NRM. We need to make sure that */
/*        the negative values of ABSTOL and RELTOL are propagated */
/*        to ABSTOLFREE and RELTOLFREE, since negative values means */
/*        that the criterion is disabled. */

	if (usetol) {
/* Computing MAX */
	    d__1 = *abstol, d__2 = *reltol * maxc2nrm;
	    abstolfree = f2cmax(d__1,d__2);
	} else {
	    abstolfree = -1.;
	}
	reltolfree = -1.;

/*        Save JPIV(NSEL+1:NSUB) into WORK(NFREE+1:2*NFREE-1) */

	if (nsel != 0) {
	    i__1 = nfree;
	    for (j = 1; j <= i__1; ++j) {
		iwork[nfree + j] = jpiv[nsel + j];
	    }
	}

	zgeqp3rk_(&mfree, &nfree, &c__0, kmaxfree, &abstolfree, &reltolfree, &
		a[nsel + 1 + (nsel + 1) * a_dim1], lda, &kfree, &
		maxc2nrmkfree, &relmaxc2nrmkfree, &jpiv[nsel + 1], &tau[nsel 
		+ 1], &work[1], lwork, &rwork[1], &iwork[1], &iinfo);

/*        Adjust JPIV */

	if (nsel != 0) {
	    i__1 = nfree;
	    for (j = 1; j <= i__1; ++j) {
		jpiv[nsel + j] = iwork[nfree + jpiv[nsel + j]];
	    }
	}

/*        1) Adjust the return value for the number of factorized */
/*           columns K for the whole submatrix A_sub. */
/*        2) MAXC2NRMK is returned transparently without change */
/*           as MAXC2NRMKFREE is returned from ZGEQP3RK. */
/*        3) Adjust the return value RELMAXC2NRMK for the whole */
/*           submatrix A_sub. We do not use RELMAXC2NRMKFREE */
/*           returned from ZGEQP3RK. */

	*k += kfree;
	*maxc2nrmk = maxc2nrmkfree;
	*relmaxc2nrmk = *maxc2nrmk / maxc2nrm;

    } else {

/*        Set norms to zero */

	*maxc2nrmk = 0.;
	*relmaxc2nrmk = 0.;

    }

/*     Now, MRESID and NRESID is the number of rows and columns */
/*     respectively in  A_free_resid = A(K+1:MSUB,K+1:NSUB). */

    mresid = mfree - kfree;
    nresid = nfree - kfree;

    if (f2cmin(mresid,nresid) != 0) {
	*fnrmk = zlange_("F", &mresid, &nresid, &a[*k + 1 + (*k + 1) * a_dim1]
		, lda, &work[1]);
    } else {
	*fnrmk = 0.;
    }

/*     ================================================================== */

/*     Return the matrix C. */

    if (returnc && *k > 0) {

	if (returnx) {

/*        Copy the selected K columns of the original matrix A (that was */
/*        saved into the array X) into the array C according to */
/*        the pivot array JPIV. If we return X, then the matrix A is */
/*        saved in the array X, and it is faster to copy into C than */
/*        doing column permutation in place, as it is the ELSE case. */

	    i__1 = *k;
	    for (j = 1; j <= i__1; ++j) {
		zcopy_(m, &x[jpiv[j] * x_dim1 + 1], &c__1, &c__[j * c_dim1 + 
			1], &c__1);
	    }

	} else {

/*        Swap the columns of the original matrix A copied into */
/*        the array C in place. */

/*        The original M-by-N matrix A was copied into the array C at */
/*        the beginning of the routine, if RETURNC = .TRUE.. */
/*        Apply the column permutation matrix P stored in JPIV(1:K) */
/*        to the columns 1:K in the M-by-N array C in place. */
/*        After column interchanges, the first K columns of C should */
/*        be the same as the first K columns of A*P, i.e. */
/*        (A*P)(1:M,1:K) = C(1:M,1:K). The complexity of this algorithm */
/*        is f2cmin(K,N-1). */

/*        Index I is the original column index in the */
/*        array C before interchanges. */
/*        J is the current column index of the original column I at */
/*        each step of interchanges. */

/*        Auxiliary array IWORK(1:N) stores the inverse P_inv(J) */
/*        of the current column permutation matrix P(J) at each */
/*        column interchange step J only for the array */
/*        values >= J:N. */
/*        C_prev  = P_inv(J) * C_next. */
/*        Each IWORK(I) contains JJ corresponding to I */
/*        Initialize IWORK(1:N) as (1:N). */

	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		iwork[i__] = i__;
	    }

/*        Auxiliary array IWORK(N+1:2N) stores the current column */
/*        permutation matrix P_(J) at each column interchange step J */
/*        only for the array index >= J:N. */
/*        C_prev * P_(J) = C_next. */
/*        Each IWORK(N+JJ) contains I corresponding to JJ. */
/*        Initialize IWORK(N+1:2*N) as (1:N). */

	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		iwork[*n + j] = j;
	    }

/*        Loop over the columns J = ( 1:f2cmin( K, N-1 ) ) in C. */

/* Computing MIN */
	    i__2 = *k, i__3 = *n - 1;
	    i__1 = f2cmin(i__2,i__3);
	    for (j = 1; j <= i__1; ++j) {

/*           IP is the original pivot column, i.e. is the original */
/*           column that should be placed in the current column index */
/*           J in the array C. */

		ip = jpiv[j];

/*           I is the original column that is */
/*           currently in the column index J in the array C after */
/*           previous column interchanges. */

		i__ = iwork[*n + j];

		if (i__ != ip) {

/*              JP is the current index of the original pivot */
/*              column IP in the array C after previous column */
/*              interchanges. */

		    jp = iwork[ip];
/*              Swap the original pivot column IP = JPIV( J ), */
/*              at the current pivot index JP = IWORK( IP ) into */
/*              index J. */

		    zswap_(m, &c__[j * c_dim1 + 1], &c__1, &c__[jp * c_dim1 + 
			    1], &c__1);

/*              Update the array IWORK(1:N) for the original column */
/*              I that was swapped with IP. */

		    iwork[i__] = iwork[ip];

/*              Update the array IWORK(N+1:2*N) for the current column */
/*              index JP that was swapped with the current column */
/*              index J. */

		    iwork[*n + jp] = iwork[*n + j];

		}

	    }

/*     End of ELSE( RETURNX ) */

	}

/*     End of IF( RETURNC .AND. K.GT.0 ) */

    }

/*     ================================================================== */

/*     Return the matrix X. */

    if (returnx && *k > 0) {

/*        We need to use C and A to compute X = pseudoinv(C) * A, as */
/*        the linear least squares solution to the overdetermined system */
/*        C*X = A. We use LLS routine that uses the QR factorization. For */
/*        that purpose, we store the matrix C into the array QRC. */
/*        The matrix A was copied into the array X at the beginning */
/*        of the routine. */

	zlacpy_("F", m, k, &c__[c_offset], ldc, &qrc[qrc_offset], ldqrc);

	zgels_("N", m, k, n, &qrc[qrc_offset], ldqrc, &x[x_offset], ldx, &
		work[1], lwork, &iinfo);
	*info = iinfo;

    }

    z__1.r = (doublereal) lwkopt, z__1.i = 0.;
    work[1].r = z__1.r, work[1].i = z__1.i;
    rwork[1] = (doublereal) lrwkopt;
    iwork[1] = liwkopt;

/*     End of ZGECXX */

    return 0;
} /* zgecxx_ */

