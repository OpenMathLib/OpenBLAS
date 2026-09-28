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
static integer c__1 = 1;
static logical c_true = TRUE_;

/* Subroutine */ void claqz3_(logical *ilschur, logical *ilq, logical *ilz, 
	integer *n, integer *ilo, integer *ihi, integer *nshifts, integer *
	nblock_desired__, complex *alpha, complex *beta, complex *a, integer *
	lda, complex *b, integer *ldb, complex *q, integer *ldq, complex *z__,
	 integer *ldz, complex *qc, integer *ldqc, complex *zc, integer *ldzc,
	 complex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, qc_dim1, qc_offset, zc_dim1, zc_offset, i__1, i__2, 
	    i__3, i__4, i__5;
    complex q__1, q__2, q__3;

    /* Local variables */
    real c__;
    integer i__, j, k;
    complex s;
    integer np, ns;
    complex temp;
    extern /* Subroutine */ void crot_(integer *, complex *, integer *, 
	    complex *, integer *, real *, complex *);
    integer npos;
    complex temp2, temp3;
    real scale;
    extern /* Subroutine */ void cgemm_(char *, char *, integer *, integer *, 
	    integer *, complex *, complex *, integer *, complex *, integer *, 
	    complex *, complex *, integer *), claqz1_(logical 
	    *, logical *, integer *, integer *, integer *, integer *, complex 
	    *, integer *, complex *, integer *, integer *, integer *, complex 
	    *, integer *, integer *, integer *, complex *, integer *), 
	    slabad_(real *, real *);
    extern real slamch_(char *);
    integer nblock;
    extern /* Subroutine */ void claset_(char *, integer *, integer *, complex 
	    *, complex *, complex *, integer *), clartg_(complex *, 
	    complex *, real *, complex *, complex *);
    real safmin;
    extern /* Subroutine */ void xerbla_(char *, integer *, ftnlen);
    real safmax;
    extern /* Subroutine */ void clacpy_(char *, integer *, integer *, complex 
	    *, integer *, complex *, integer *);
    integer ishift, istopb, swidth, istopm, sheight, istartb, istartm;

/*     Function arguments */
/*     Parameters */
/*     Local scalars */
/*     External Functions */
    /* Parameter adjustments */
    --alpha;
    --beta;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    qc_dim1 = *ldqc;
    qc_offset = 1 + qc_dim1;
    qc -= qc_offset;
    zc_dim1 = *ldzc;
    zc_offset = 1 + zc_dim1;
    zc -= zc_offset;
    --work;

    /* Function Body */
    *info = 0;
    if (*nblock_desired__ < *nshifts + 1) {
	*info = -8;
    }
    if (*lwork == -1) {
/*        workspace query, quick return */
	i__1 = *n * *nblock_desired__;
	work[1].r = (real) i__1, work[1].i = 0.f;
	return;
    } else if (*lwork < *n * *nblock_desired__) {
	*info = -25;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CLAQZ3", &i__1, (ftnlen)6);
	return;
    }

/*     Executable statements */

/*     Get machine constants */
    safmin = slamch_("SAFE MINIMUM");
    safmax = 1.f / safmin;
    slabad_(&safmin, &safmax);
    if (*ilo >= *ihi) {
	return;
    }
    if (*ilschur) {
	istartm = 1;
	istopm = *n;
    } else {
	istartm = *ilo;
	istopm = *ihi;
    }
    ns = *nshifts;
/* Computing MAX */
    i__1 = *nblock_desired__ - ns;
    npos = f2cmax(i__1,1);
/*     The following block introduces the shifts and chases */
/*     them down one by one just enough to make space for */
/*     the other shifts. The near-the-diagonal block is */
/*     of size (ns+1) x ns. */
    i__1 = ns + 1;
    i__2 = ns + 1;
    claset_("FULL", &i__1, &i__2, &c_b1, &c_b2, &qc[qc_offset], ldqc);
    claset_("FULL", &ns, &ns, &c_b1, &c_b2, &zc[zc_offset], ldzc);
    i__1 = ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*        Introduce the shift */
	scale = sqrt(c_abs(&alpha[i__])) * sqrt(c_abs(&beta[i__]));
	if (scale >= safmin && scale <= safmax) {
	    i__2 = i__;
	    i__3 = i__;
	    q__1.r = alpha[i__3].r / scale, q__1.i = alpha[i__3].i / scale;
	    alpha[i__2].r = q__1.r, alpha[i__2].i = q__1.i;
	    i__2 = i__;
	    i__3 = i__;
	    q__1.r = beta[i__3].r / scale, q__1.i = beta[i__3].i / scale;
	    beta[i__2].r = q__1.r, beta[i__2].i = q__1.i;
	}
	i__2 = i__;
	i__3 = *ilo + *ilo * a_dim1;
	q__2.r = beta[i__2].r * a[i__3].r - beta[i__2].i * a[i__3].i, q__2.i =
		 beta[i__2].r * a[i__3].i + beta[i__2].i * a[i__3].r;
	i__4 = i__;
	i__5 = *ilo + *ilo * b_dim1;
	q__3.r = alpha[i__4].r * b[i__5].r - alpha[i__4].i * b[i__5].i, 
		q__3.i = alpha[i__4].r * b[i__5].i + alpha[i__4].i * b[i__5]
		.r;
	q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	temp2.r = q__1.r, temp2.i = q__1.i;
	i__2 = i__;
	i__3 = *ilo + 1 + *ilo * a_dim1;
	q__1.r = beta[i__2].r * a[i__3].r - beta[i__2].i * a[i__3].i, q__1.i =
		 beta[i__2].r * a[i__3].i + beta[i__2].i * a[i__3].r;
	temp3.r = q__1.r, temp3.i = q__1.i;
	if (c_abs(&temp2) > safmax || c_abs(&temp3) > safmax) {
	    temp2.r = 1.f, temp2.i = 0.f;
	    temp3.r = 0.f, temp3.i = 0.f;
	}
	clartg_(&temp2, &temp3, &c__, &s, &temp);
	crot_(&ns, &a[*ilo + *ilo * a_dim1], lda, &a[*ilo + 1 + *ilo * a_dim1]
		, lda, &c__, &s);
	crot_(&ns, &b[*ilo + *ilo * b_dim1], ldb, &b[*ilo + 1 + *ilo * b_dim1]
		, ldb, &c__, &s);
	i__2 = ns + 1;
	r_cnjg(&q__1, &s);
	crot_(&i__2, &qc[qc_dim1 + 1], &c__1, &qc[(qc_dim1 << 1) + 1], &c__1, 
		&c__, &q__1);
/*        Chase the shift down */
	i__2 = ns - i__;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = *ihi - *ilo + 1;
	    i__4 = ns + 1;
	    claqz1_(&c_true, &c_true, &j, &c__1, &ns, &i__3, &a[*ilo + *ilo * 
		    a_dim1], lda, &b[*ilo + *ilo * b_dim1], ldb, &i__4, &c__1,
		     &qc[qc_offset], ldqc, &ns, &c__1, &zc[zc_offset], ldzc);
	}
    }
/*     Update the rest of the pencil */
/*     Update A(ilo:ilo+ns,ilo+ns:istopm) and B(ilo:ilo+ns,ilo+ns:istopm) */
/*     from the left with Qc(1:ns+1,1:ns+1)' */
    sheight = ns + 1;
    swidth = istopm - (*ilo + ns) + 1;
    if (swidth > 0) {
	cgemm_("C", "N", &sheight, &swidth, &sheight, &c_b2, &qc[qc_offset], 
		ldqc, &a[*ilo + (*ilo + ns) * a_dim1], lda, &c_b1, &work[1], &
		sheight);
	clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[*ilo + (*ilo 
		+ ns) * a_dim1], lda);
	cgemm_("C", "N", &sheight, &swidth, &sheight, &c_b2, &qc[qc_offset], 
		ldqc, &b[*ilo + (*ilo + ns) * b_dim1], ldb, &c_b1, &work[1], &
		sheight);
	clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[*ilo + (*ilo 
		+ ns) * b_dim1], ldb);
    }
    if (*ilq) {
	cgemm_("N", "N", n, &sheight, &sheight, &c_b2, &q[*ilo * q_dim1 + 1], 
		ldq, &qc[qc_offset], ldqc, &c_b1, &work[1], n);
	clacpy_("ALL", n, &sheight, &work[1], n, &q[*ilo * q_dim1 + 1], ldq);
    }
/*     Update A(istartm:ilo-1,ilo:ilo+ns-1) and B(istartm:ilo-1,ilo:ilo+ns-1) */
/*     from the right with Zc(1:ns,1:ns) */
    sheight = *ilo - 1 - istartm + 1;
    swidth = ns;
    if (sheight > 0) {
	cgemm_("N", "N", &sheight, &swidth, &swidth, &c_b2, &a[istartm + *ilo 
		* a_dim1], lda, &zc[zc_offset], ldzc, &c_b1, &work[1], &
		sheight);
	clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[istartm + *
		ilo * a_dim1], lda);
	cgemm_("N", "N", &sheight, &swidth, &swidth, &c_b2, &b[istartm + *ilo 
		* b_dim1], ldb, &zc[zc_offset], ldzc, &c_b1, &work[1], &
		sheight);
	clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[istartm + *
		ilo * b_dim1], ldb);
    }
    if (*ilz) {
	cgemm_("N", "N", n, &swidth, &swidth, &c_b2, &z__[*ilo * z_dim1 + 1], 
		ldz, &zc[zc_offset], ldzc, &c_b1, &work[1], n);
	clacpy_("ALL", n, &swidth, &work[1], n, &z__[*ilo * z_dim1 + 1], ldz);
    }
/*     The following block chases the shifts down to the bottom */
/*     right block. If possible, a shift is moved down npos */
/*     positions at a time */
    k = *ilo;
    while(k < *ihi - ns) {
/* Computing MIN */
	i__1 = *ihi - ns - k;
	np = f2cmin(i__1,npos);
/*        Size of the near-the-diagonal block */
	nblock = ns + np;
/*        istartb points to the first row we will be updating */
	istartb = k + 1;
/*        istopb points to the last column we will be updating */
	istopb = k + nblock - 1;
	i__1 = ns + np;
	i__2 = ns + np;
	claset_("FULL", &i__1, &i__2, &c_b1, &c_b2, &qc[qc_offset], ldqc);
	i__1 = ns + np;
	i__2 = ns + np;
	claset_("FULL", &i__1, &i__2, &c_b1, &c_b2, &zc[zc_offset], ldzc);
/*        Near the diagonal shift chase */
	for (i__ = ns - 1; i__ >= 0; --i__) {
	    i__1 = np - 1;
	    for (j = 0; j <= i__1; ++j) {
/*              Move down the block with index k+i+j, updating */
/*              the (ns+np x ns+np) block: */
/*              (k:k+ns+np,k:k+ns+np-1) */
		i__2 = k + i__ + j;
		i__3 = k + 1;
		claqz1_(&c_true, &c_true, &i__2, &istartb, &istopb, ihi, &a[
			a_offset], lda, &b[b_offset], ldb, &nblock, &i__3, &
			qc[qc_offset], ldqc, &nblock, &k, &zc[zc_offset], 
			ldzc);
	    }
	}
/*        Update rest of the pencil */
/*        Update A(k+1:k+ns+np, k+ns+np:istopm) and */
/*        B(k+1:k+ns+np, k+ns+np:istopm) */
/*        from the left with Qc(1:ns+np,1:ns+np)' */
	sheight = ns + np;
	swidth = istopm - (k + ns + np) + 1;
	if (swidth > 0) {
	    cgemm_("C", "N", &sheight, &swidth, &sheight, &c_b2, &qc[
		    qc_offset], ldqc, &a[k + 1 + (k + ns + np) * a_dim1], lda,
		     &c_b1, &work[1], &sheight);
	    clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[k + 1 + (
		    k + ns + np) * a_dim1], lda);
	    cgemm_("C", "N", &sheight, &swidth, &sheight, &c_b2, &qc[
		    qc_offset], ldqc, &b[k + 1 + (k + ns + np) * b_dim1], ldb,
		     &c_b1, &work[1], &sheight);
	    clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[k + 1 + (
		    k + ns + np) * b_dim1], ldb);
	}
	if (*ilq) {
	    cgemm_("N", "N", n, &nblock, &nblock, &c_b2, &q[(k + 1) * q_dim1 
		    + 1], ldq, &qc[qc_offset], ldqc, &c_b1, &work[1], n);
	    clacpy_("ALL", n, &nblock, &work[1], n, &q[(k + 1) * q_dim1 + 1], 
		    ldq);
	}
/*        Update A(istartm:k,k:k+ns+npos-1) and B(istartm:k,k:k+ns+npos-1) */
/*        from the right with Zc(1:ns+np,1:ns+np) */
	sheight = k - istartm + 1;
	swidth = nblock;
	if (sheight > 0) {
	    cgemm_("N", "N", &sheight, &swidth, &swidth, &c_b2, &a[istartm + 
		    k * a_dim1], lda, &zc[zc_offset], ldzc, &c_b1, &work[1], &
		    sheight);
	    clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[istartm 
		    + k * a_dim1], lda);
	    cgemm_("N", "N", &sheight, &swidth, &swidth, &c_b2, &b[istartm + 
		    k * b_dim1], ldb, &zc[zc_offset], ldzc, &c_b1, &work[1], &
		    sheight);
	    clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[istartm 
		    + k * b_dim1], ldb);
	}
	if (*ilz) {
	    cgemm_("N", "N", n, &nblock, &nblock, &c_b2, &z__[k * z_dim1 + 1],
		     ldz, &zc[zc_offset], ldzc, &c_b1, &work[1], n);
	    clacpy_("ALL", n, &nblock, &work[1], n, &z__[k * z_dim1 + 1], ldz);
	}
	k += np;
    }
/*     The following block removes the shifts from the bottom right corner */
/*     one by one. Updates are initially applied to A(ihi-ns+1:ihi,ihi-ns:ihi). */
    claset_("FULL", &ns, &ns, &c_b1, &c_b2, &qc[qc_offset], ldqc);
    i__1 = ns + 1;
    i__2 = ns + 1;
    claset_("FULL", &i__1, &i__2, &c_b1, &c_b2, &zc[zc_offset], ldzc);
/*     istartb points to the first row we will be updating */
    istartb = *ihi - ns + 1;
/*     istopb points to the last column we will be updating */
    istopb = *ihi;
    i__1 = ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*        Chase the shift down to the bottom right corner */
	i__2 = *ihi - 1;
	for (ishift = *ihi - i__; ishift <= i__2; ++ishift) {
	    i__3 = *ihi - ns + 1;
	    i__4 = ns + 1;
	    i__5 = *ihi - ns;
	    claqz1_(&c_true, &c_true, &ishift, &istartb, &istopb, ihi, &a[
		    a_offset], lda, &b[b_offset], ldb, &ns, &i__3, &qc[
		    qc_offset], ldqc, &i__4, &i__5, &zc[zc_offset], ldzc);
	}
    }
/*     Update rest of the pencil */
/*     Update A(ihi-ns+1:ihi, ihi+1:istopm) */
/*     from the left with Qc(1:ns,1:ns)' */
    sheight = ns;
    swidth = istopm - (*ihi + 1) + 1;
    if (swidth > 0) {
	cgemm_("C", "N", &sheight, &swidth, &sheight, &c_b2, &qc[qc_offset], 
		ldqc, &a[*ihi - ns + 1 + (*ihi + 1) * a_dim1], lda, &c_b1, &
		work[1], &sheight);
	clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[*ihi - ns + 
		1 + (*ihi + 1) * a_dim1], lda);
	cgemm_("C", "N", &sheight, &swidth, &sheight, &c_b2, &qc[qc_offset], 
		ldqc, &b[*ihi - ns + 1 + (*ihi + 1) * b_dim1], ldb, &c_b1, &
		work[1], &sheight);
	clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[*ihi - ns + 
		1 + (*ihi + 1) * b_dim1], ldb);
    }
    if (*ilq) {
	cgemm_("N", "N", n, &ns, &ns, &c_b2, &q[(*ihi - ns + 1) * q_dim1 + 1],
		 ldq, &qc[qc_offset], ldqc, &c_b1, &work[1], n);
	clacpy_("ALL", n, &ns, &work[1], n, &q[(*ihi - ns + 1) * q_dim1 + 1], 
		ldq);
    }
/*     Update A(istartm:ihi-ns,ihi-ns:ihi) */
/*     from the right with Zc(1:ns+1,1:ns+1) */
    sheight = *ihi - ns - istartm + 1;
    swidth = ns + 1;
    if (sheight > 0) {
	cgemm_("N", "N", &sheight, &swidth, &swidth, &c_b2, &a[istartm + (*
		ihi - ns) * a_dim1], lda, &zc[zc_offset], ldzc, &c_b1, &work[
		1], &sheight);
	clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[istartm + (*
		ihi - ns) * a_dim1], lda);
	cgemm_("N", "N", &sheight, &swidth, &swidth, &c_b2, &b[istartm + (*
		ihi - ns) * b_dim1], ldb, &zc[zc_offset], ldzc, &c_b1, &work[
		1], &sheight);
	clacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[istartm + (*
		ihi - ns) * b_dim1], ldb);
    }
    if (*ilz) {
	i__1 = ns + 1;
	i__2 = ns + 1;
	cgemm_("N", "N", n, &i__1, &i__2, &c_b2, &z__[(*ihi - ns) * z_dim1 + 
		1], ldz, &zc[zc_offset], ldzc, &c_b1, &work[1], n);
	i__1 = ns + 1;
	clacpy_("ALL", n, &i__1, &work[1], n, &z__[(*ihi - ns) * z_dim1 + 1], 
		ldz);
    }
    return;
} /* claqz3_ */

