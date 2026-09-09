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
#define d_abs(x) (fabs(*(x)))
#define d_acos(x) (acos(*(x)))
#define d_asin(x) (asin(*(x)))
#define d_atan(x) (atan(*(x)))
#define d_atn2(x, y) (atan2(*(x),*(y)))
#define d_cos(x) (cos(*(x)))
#define d_cosh(x) (cosh(*(x)))
#define d_dim(__a, __b) ( *(__a) > *(__b) ? *(__a) - *(__b) : 0.0 )
#define d_exp(x) (exp(*(x)))
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
#define s_cat(lpp, rpp, rnp, np, llp) { 	ftnlen i, nc, ll; char *f__rp, *lp; 	ll = (llp); lp = (lpp); 	for(i=0; i < (int)*(np); ++i) {         	nc = ll; 	        if((rnp)[i] < nc) nc = (rnp)[i]; 	        ll -= nc;         	f__rp = (rpp)[i]; 	        while(--nc >= 0) *lp++ = *(f__rp)++;         } 	while(--ll >= 0) *lp++ = ' '; }
#define s_cmp(a,b,c,d) ((integer)strncmp((a),(b),f2cmin((c),(d))))
#define s_copy(A,B,C,D) { int __i,__m; for (__i=0, __m=f2cmin((C),(D)); __i<__m && (B)[__i] != 0; ++__i) (A)[__i] = (B)[__i]; }
#define sig_die(s, kill) { exit(1); }
#define s_stop(s, n) {exit(0);}
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

/* Table of constant values */

static logical c_true = TRUE_;
static integer c_n1 = -1;
static integer c__1 = 1;
static real c_b16 = 0.f;
static real c_b17 = 1.f;

/*  ===================================================================== */
/* Subroutine */ void slaqz3_(logical *ilschur, logical *ilq, logical *ilz, 
	integer *n, integer *ilo, integer *ihi, integer *nw, real *a, integer 
	*lda, real *b, integer *ldb, real *q, integer *ldq, real *z__, 
	integer *ldz, integer *ns, integer *nd, real *alphar, real *alphai, 
	real *beta, real *qc, integer *ldqc, real *zc, integer *ldzc, real *
	work, integer *lwork, integer *rec, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, qc_dim1, qc_offset, zc_dim1, zc_offset, i__1, i__2, 
	    i__3, i__4;
    real r__1, r__2, r__3, r__4, r__5, r__6;

    /* Local variables */
    integer lworkreq, k;
    real s, c1;
    integer k2;
    real s1;
    integer jw, imk;
    real ulp;
    integer stgexc_info__, ifst;
    real temp;
    integer ilst;
    extern /* Subroutine */ void srot_(integer *, real *, integer *, real *, 
	    integer *, real *, real *), slag2_(real *, integer *, real *, 
	    integer *, real *, real *, real *, real *, real *, real *);
    logical bulge;
    extern /* Subroutine */ void sgemm_(char *, char *, integer *, integer *, 
	    integer *, real *, real *, integer *, real *, integer *, real *, 
	    real *, integer *);
    integer kwbot;
    real mktmp;
    integer kwtop, qz_small_info__;
    extern /* Subroutine */ void slaqz0_(char *, char *, char *, integer *, 
	    integer *, integer *, real *, integer *, real *, integer *, real *
	    , real *, real *, real *, integer *, real *, integer *, real *, 
	    integer *, integer *, integer *), slaqz2_(
	    logical *, logical *, integer *, integer *, integer *, integer *, 
	    real *, integer *, real *, integer *, integer *, integer *, real *
	    , integer *, integer *, integer *, real *, integer *);
    extern real slamch_(char *);
    real safmin;
    extern /* Subroutine */ void xerbla_(char *, integer *, ftnlen);
    real safmax;
    extern /* Subroutine */ void slacpy_(char *, integer *, integer *, real *, 
	    integer *, real *, integer *), slaset_(char *, integer *, 
	    integer *, real *, real *, real *, integer *), stgexc_(
	    logical *, logical *, integer *, real *, integer *, real *, 
	    integer *, real *, integer *, real *, integer *, integer *, 
	    integer *, real *, integer *, integer *), slartg_(real *, real *, 
	    real *, real *, real *);
    integer istopm;
    real smlnum;
    extern real sroundup_lwork__(integer *);
    integer istartm;

/*     Arguments */
/*     Parameters */
/*     Local Scalars */
/*     External Functions */
    /* Parameter adjustments */
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
    --alphar;
    --alphai;
    --beta;
    qc_dim1 = *ldqc;
    qc_offset = 1 + qc_dim1;
    qc -= qc_offset;
    zc_dim1 = *ldzc;
    zc_offset = 1 + zc_dim1;
    zc -= zc_offset;
    --work;

    /* Function Body */
    *info = 0;
/*     Set up deflation window */
/* Computing MIN */
    i__1 = *nw, i__2 = *ihi - *ilo + 1;
    jw = f2cmin(i__1,i__2);
    kwtop = *ihi - jw + 1;
    if (kwtop == *ilo) {
	s = 0.f;
    } else {
	s = a[kwtop + (kwtop - 1) * a_dim1];
    }
/*     Determine required workspace */
    ifst = 1;
    ilst = jw;
    stgexc_(&c_true, &c_true, &jw, &a[a_offset], lda, &b[b_offset], ldb, &qc[
	    qc_offset], ldqc, &zc[zc_offset], ldzc, &ifst, &ilst, &work[1], &
	    c_n1, &stgexc_info__);
    lworkreq = (integer) work[1];
    i__1 = *rec + 1;
    slaqz0_("S", "V", "V", &jw, &c__1, &jw, &a[kwtop + kwtop * a_dim1], lda, &
	    b[kwtop + kwtop * b_dim1], ldb, &alphar[1], &alphai[1], &beta[1], 
	    &qc[qc_offset], ldqc, &zc[zc_offset], ldzc, &work[1], &c_n1, &
	    i__1, &qz_small_info__);
/* Computing MAX */
/* Computing 2nd power */
    i__3 = jw;
    i__1 = lworkreq, i__2 = (integer) work[1] + (i__3 * i__3 << 1);
    lworkreq = f2cmax(i__1,i__2);
/* Computing MAX */
/* Computing 2nd power */
    i__3 = *nw;
    i__1 = lworkreq, i__2 = *n * *nw, i__1 = f2cmax(i__1,i__2), i__2 = (i__3 * 
	    i__3 << 1) + *n;
    lworkreq = f2cmax(i__1,i__2);
    if (*lwork == -1) {
/*        workspace query, quick return */
	work[1] = sroundup_lwork__(&lworkreq);
	return;
    } else if (*lwork < lworkreq) {
	*info = -26;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SLAQZ3", &i__1, (ftnlen)6);
	return;
    }
/*     Get machine constants */
    safmin = slamch_("SAFE MINIMUM");
    safmax = 1.f / safmin;
    ulp = slamch_("PRECISION");
    smlnum = safmin * ((real) (*n) / ulp);
    if (*ihi == kwtop) {
/*        1 by 1 deflation window, just try a regular deflation */
	alphar[kwtop] = a[kwtop + kwtop * a_dim1];
	alphai[kwtop] = 0.f;
	beta[kwtop] = b[kwtop + kwtop * b_dim1];
	*ns = 1;
	*nd = 0;
/* Computing MAX */
	r__2 = smlnum, r__3 = ulp * (r__1 = a[kwtop + kwtop * a_dim1], abs(
		r__1));
	if (abs(s) <= f2cmax(r__2,r__3)) {
	    *ns = 0;
	    *nd = 1;
	    if (kwtop > *ilo) {
		a[kwtop + (kwtop - 1) * a_dim1] = 0.f;
	    }
	}
    }
/*     Store window in case of convergence failure */
    slacpy_("ALL", &jw, &jw, &a[kwtop + kwtop * a_dim1], lda, &work[1], &jw);
/* Computing 2nd power */
    i__1 = jw;
    slacpy_("ALL", &jw, &jw, &b[kwtop + kwtop * b_dim1], ldb, &work[i__1 * 
	    i__1 + 1], &jw);
/*     Transform window to real schur form */
    slaset_("FULL", &jw, &jw, &c_b16, &c_b17, &qc[qc_offset], ldqc)
	    ;
    slaset_("FULL", &jw, &jw, &c_b16, &c_b17, &zc[zc_offset], ldzc)
	    ;
/* Computing 2nd power */
    i__1 = jw;
/* Computing 2nd power */
    i__3 = jw;
    i__2 = *lwork - (i__3 * i__3 << 1);
    i__4 = *rec + 1;
    slaqz0_("S", "V", "V", &jw, &c__1, &jw, &a[kwtop + kwtop * a_dim1], lda, &
	    b[kwtop + kwtop * b_dim1], ldb, &alphar[1], &alphai[1], &beta[1], 
	    &qc[qc_offset], ldqc, &zc[zc_offset], ldzc, &work[(i__1 * i__1 << 
	    1) + 1], &i__2, &i__4, &qz_small_info__);
    if (qz_small_info__ != 0) {
/*        Convergence failure, restore the window and exit */
	*nd = 0;
	*ns = jw - qz_small_info__;
	slacpy_("ALL", &jw, &jw, &work[1], &jw, &a[kwtop + kwtop * a_dim1], 
		lda);
/* Computing 2nd power */
	i__1 = jw;
	slacpy_("ALL", &jw, &jw, &work[i__1 * i__1 + 1], &jw, &b[kwtop + 
		kwtop * b_dim1], ldb);
	return;
    }
/*     Deflation detection loop */
    if (kwtop == *ilo || s == 0.f) {
	kwbot = kwtop - 1;
    } else {
	kwbot = *ihi;
	k = 1;
	k2 = 1;
	while(k <= jw) {
	    bulge = FALSE_;
	    if (kwbot - kwtop + 1 >= 2) {
		bulge = a[kwbot + (kwbot - 1) * a_dim1] != 0.f;
	    }
	    if (bulge) {
/*              Try to deflate complex conjugate eigenvalue pair */
		temp = (r__3 = a[kwbot + kwbot * a_dim1], abs(r__3)) + sqrt((
			r__1 = a[kwbot + (kwbot - 1) * a_dim1], abs(r__1))) * 
			sqrt((r__2 = a[kwbot - 1 + kwbot * a_dim1], abs(r__2))
			);
		if (temp == 0.f) {
		    temp = abs(s);
		}
/* Computing MAX */
		r__3 = (r__1 = s * qc[(kwbot - kwtop) * qc_dim1 + 1], abs(
			r__1)), r__4 = (r__2 = s * qc[(kwbot - kwtop + 1) * 
			qc_dim1 + 1], abs(r__2));
/* Computing MAX */
		r__5 = smlnum, r__6 = ulp * temp;
		if (f2cmax(r__3,r__4) <= f2cmax(r__5,r__6)) {
/*                 Deflatable */
		    kwbot += -2;
		} else {
/*                 Not deflatable, move out of the way */
		    ifst = kwbot - kwtop + 1;
		    ilst = k2;
		    stgexc_(&c_true, &c_true, &jw, &a[kwtop + kwtop * a_dim1],
			     lda, &b[kwtop + kwtop * b_dim1], ldb, &qc[
			    qc_offset], ldqc, &zc[zc_offset], ldzc, &ifst, &
			    ilst, &work[1], lwork, &stgexc_info__);
		    k2 += 2;
		}
		k += 2;
	    } else {
/*              Try to deflate real eigenvalue */
		temp = (r__1 = a[kwbot + kwbot * a_dim1], abs(r__1));
		if (temp == 0.f) {
		    temp = abs(s);
		}
/* Computing MAX */
		r__2 = ulp * temp;
		if ((r__1 = s * qc[(kwbot - kwtop + 1) * qc_dim1 + 1], abs(
			r__1)) <= f2cmax(r__2,smlnum)) {
/*                 Deflatable */
		    --kwbot;
		} else {
/*                 Not deflatable, move out of the way */
		    ifst = kwbot - kwtop + 1;
		    ilst = k2;
		    stgexc_(&c_true, &c_true, &jw, &a[kwtop + kwtop * a_dim1],
			     lda, &b[kwtop + kwtop * b_dim1], ldb, &qc[
			    qc_offset], ldqc, &zc[zc_offset], ldzc, &ifst, &
			    ilst, &work[1], lwork, &stgexc_info__);
		    ++k2;
		}
		++k;
	    }
	}
    }
/*     Store eigenvalues */
    *nd = *ihi - kwbot;
    *ns = jw - *nd;
    k = kwtop;
    while(k <= *ihi) {
	bulge = FALSE_;
	if (k < *ihi) {
	    if (a[k + 1 + k * a_dim1] != 0.f) {
		bulge = TRUE_;
	    }
	}
	if (bulge) {
/*           2x2 eigenvalue block */
	    slag2_(&a[k + k * a_dim1], lda, &b[k + k * b_dim1], ldb, &safmin, 
		    &beta[k], &beta[k + 1], &alphar[k], &alphar[k + 1], &
		    alphai[k]);
	    alphai[k + 1] = -alphai[k];
	    k += 2;
	} else {
/*           1x1 eigenvalue block */
	    alphar[k] = a[k + k * a_dim1];
	    alphai[k] = 0.f;
	    beta[k] = b[k + k * b_dim1];
	    ++k;
	}
    }
    if (kwtop != *ilo && s != 0.f) {
/*        Reflect spike back, this will create optimally packed bulges */
/*         A( KWTOP:KWBOT, KWTOP-1 ) = A( KWTOP, KWTOP-1 )*QC( 1, */
/*     $      1:JW-ND ) */
	i__1 = jw - *nd;
	for (imk = 1; imk <= i__1; ++imk) {
	    mktmp = a[kwtop + (kwtop - 1) * a_dim1] * qc[imk * qc_dim1 + 1];
	}
	i__1 = kwbot;
	for (imk = kwtop; imk <= i__1; ++imk) {
	    a[imk + (kwtop - 1) * a_dim1] = mktmp;
	}
	i__1 = kwtop;
	for (k = kwbot - 1; k >= i__1; --k) {
	    slartg_(&a[k + (kwtop - 1) * a_dim1], &a[k + 1 + (kwtop - 1) * 
		    a_dim1], &c1, &s1, &temp);
	    a[k + (kwtop - 1) * a_dim1] = temp;
	    a[k + 1 + (kwtop - 1) * a_dim1] = 0.f;
/* Computing MAX */
	    i__2 = kwtop, i__3 = k - 1;
	    k2 = f2cmax(i__2,i__3);
	    i__2 = *ihi - k2 + 1;
	    srot_(&i__2, &a[k + k2 * a_dim1], lda, &a[k + 1 + k2 * a_dim1], 
		    lda, &c1, &s1);
	    i__2 = *ihi - (k - 1) + 1;
	    srot_(&i__2, &b[k + (k - 1) * b_dim1], ldb, &b[k + 1 + (k - 1) * 
		    b_dim1], ldb, &c1, &s1);
	    srot_(&jw, &qc[(k - kwtop + 1) * qc_dim1 + 1], &c__1, &qc[(k + 1 
		    - kwtop + 1) * qc_dim1 + 1], &c__1, &c1, &s1);
	}
/*        Chase bulges down */
	istartm = kwtop;
	istopm = *ihi;
	k = kwbot - 1;
	while(k >= kwtop) {
	    if (k >= kwtop + 1 && a[k + 1 + (k - 1) * a_dim1] != 0.f) {
/*              Move double pole block down and remove it */
		i__1 = kwbot - 2;
		for (k2 = k - 1; k2 <= i__1; ++k2) {
		    i__2 = kwtop + jw - 1;
		    slaqz2_(&c_true, &c_true, &k2, &kwtop, &i__2, &kwbot, &a[
			    a_offset], lda, &b[b_offset], ldb, &jw, &kwtop, &
			    qc[qc_offset], ldqc, &jw, &kwtop, &zc[zc_offset], 
			    ldzc);
		}
		k += -2;
	    } else {
/*              k points to single shift */
		i__1 = kwbot - 2;
		for (k2 = k; k2 <= i__1; ++k2) {
/*                 Move shift down */
		    slartg_(&b[k2 + 1 + (k2 + 1) * b_dim1], &b[k2 + 1 + k2 * 
			    b_dim1], &c1, &s1, &temp);
		    b[k2 + 1 + (k2 + 1) * b_dim1] = temp;
		    b[k2 + 1 + k2 * b_dim1] = 0.f;
		    i__2 = k2 + 2 - istartm + 1;
		    srot_(&i__2, &a[istartm + (k2 + 1) * a_dim1], &c__1, &a[
			    istartm + k2 * a_dim1], &c__1, &c1, &s1);
		    i__2 = k2 - istartm + 1;
		    srot_(&i__2, &b[istartm + (k2 + 1) * b_dim1], &c__1, &b[
			    istartm + k2 * b_dim1], &c__1, &c1, &s1);
		    srot_(&jw, &zc[(k2 + 1 - kwtop + 1) * zc_dim1 + 1], &c__1,
			     &zc[(k2 - kwtop + 1) * zc_dim1 + 1], &c__1, &c1, 
			    &s1);
		    slartg_(&a[k2 + 1 + k2 * a_dim1], &a[k2 + 2 + k2 * a_dim1]
			    , &c1, &s1, &temp);
		    a[k2 + 1 + k2 * a_dim1] = temp;
		    a[k2 + 2 + k2 * a_dim1] = 0.f;
		    i__2 = istopm - k2;
		    srot_(&i__2, &a[k2 + 1 + (k2 + 1) * a_dim1], lda, &a[k2 + 
			    2 + (k2 + 1) * a_dim1], lda, &c1, &s1);
		    i__2 = istopm - k2;
		    srot_(&i__2, &b[k2 + 1 + (k2 + 1) * b_dim1], ldb, &b[k2 + 
			    2 + (k2 + 1) * b_dim1], ldb, &c1, &s1);
		    srot_(&jw, &qc[(k2 + 1 - kwtop + 1) * qc_dim1 + 1], &c__1,
			     &qc[(k2 + 2 - kwtop + 1) * qc_dim1 + 1], &c__1, &
			    c1, &s1);
		}
/*              Remove the shift */
		slartg_(&b[kwbot + kwbot * b_dim1], &b[kwbot + (kwbot - 1) * 
			b_dim1], &c1, &s1, &temp);
		b[kwbot + kwbot * b_dim1] = temp;
		b[kwbot + (kwbot - 1) * b_dim1] = 0.f;
		i__1 = kwbot - istartm;
		srot_(&i__1, &b[istartm + kwbot * b_dim1], &c__1, &b[istartm 
			+ (kwbot - 1) * b_dim1], &c__1, &c1, &s1);
		i__1 = kwbot - istartm + 1;
		srot_(&i__1, &a[istartm + kwbot * a_dim1], &c__1, &a[istartm 
			+ (kwbot - 1) * a_dim1], &c__1, &c1, &s1);
		srot_(&jw, &zc[(kwbot - kwtop + 1) * zc_dim1 + 1], &c__1, &zc[
			(kwbot - 1 - kwtop + 1) * zc_dim1 + 1], &c__1, &c1, &
			s1);
		--k;
	    }
	}
    }
/*     Apply Qc and Zc to rest of the matrix */
    if (*ilschur) {
	istartm = 1;
	istopm = *n;
    } else {
	istartm = *ilo;
	istopm = *ihi;
    }
    if (istopm - *ihi > 0) {
	i__1 = istopm - *ihi;
	sgemm_("T", "N", &jw, &i__1, &jw, &c_b17, &qc[qc_offset], ldqc, &a[
		kwtop + (*ihi + 1) * a_dim1], lda, &c_b16, &work[1], &jw);
	i__1 = istopm - *ihi;
	slacpy_("ALL", &jw, &i__1, &work[1], &jw, &a[kwtop + (*ihi + 1) * 
		a_dim1], lda);
	i__1 = istopm - *ihi;
	sgemm_("T", "N", &jw, &i__1, &jw, &c_b17, &qc[qc_offset], ldqc, &b[
		kwtop + (*ihi + 1) * b_dim1], ldb, &c_b16, &work[1], &jw);
	i__1 = istopm - *ihi;
	slacpy_("ALL", &jw, &i__1, &work[1], &jw, &b[kwtop + (*ihi + 1) * 
		b_dim1], ldb);
    }
    if (*ilq) {
	sgemm_("N", "N", n, &jw, &jw, &c_b17, &q[kwtop * q_dim1 + 1], ldq, &
		qc[qc_offset], ldqc, &c_b16, &work[1], n);
	slacpy_("ALL", n, &jw, &work[1], n, &q[kwtop * q_dim1 + 1], ldq);
    }
    if (kwtop - 1 - istartm + 1 > 0) {
	i__1 = kwtop - istartm;
	i__2 = kwtop - istartm;
	sgemm_("N", "N", &i__1, &jw, &jw, &c_b17, &a[istartm + kwtop * a_dim1]
		, lda, &zc[zc_offset], ldzc, &c_b16, &work[1], &i__2);
	i__1 = kwtop - istartm;
	i__2 = kwtop - istartm;
	slacpy_("ALL", &i__1, &jw, &work[1], &i__2, &a[istartm + kwtop * 
		a_dim1], lda);
	i__1 = kwtop - istartm;
	i__2 = kwtop - istartm;
	sgemm_("N", "N", &i__1, &jw, &jw, &c_b17, &b[istartm + kwtop * b_dim1]
		, ldb, &zc[zc_offset], ldzc, &c_b16, &work[1], &i__2);
	i__1 = kwtop - istartm;
	i__2 = kwtop - istartm;
	slacpy_("ALL", &i__1, &jw, &work[1], &i__2, &b[istartm + kwtop * 
		b_dim1], ldb);
    }
    if (*ilz) {
	sgemm_("N", "N", n, &jw, &jw, &c_b17, &z__[kwtop * z_dim1 + 1], ldz, &
		zc[zc_offset], ldzc, &c_b16, &work[1], n);
	slacpy_("ALL", n, &jw, &work[1], n, &z__[kwtop * z_dim1 + 1], ldz);
    }
    return;
} /* slaqz3_ */

