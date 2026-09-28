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
static integer c_n1 = -1;
static logical c_true = TRUE_;

/* Subroutine */ void claqz2_(logical *ilschur, logical *ilq, logical *ilz, 
	integer *n, integer *ilo, integer *ihi, integer *nw, complex *a, 
	integer *lda, complex *b, integer *ldb, complex *q, integer *ldq, 
	complex *z__, integer *ldz, integer *ns, integer *nd, complex *alpha, 
	complex *beta, complex *qc, integer *ldqc, complex *zc, integer *ldzc,
	 complex *work, integer *lwork, real *rwork, integer *rec, integer *
	info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, qc_dim1, qc_offset, zc_dim1, zc_offset, i__1, i__2, 
	    i__3, i__4;
    real r__1, r__2;
    complex q__1, q__2;

    /* Local variables */
    integer lworkreq, k;
    complex s;
    real c1;
    integer k2;
    complex s1;
    integer jw, imk;
    real ulp;
    integer ctgexc_info__, ifst;
    complex temp;
    extern /* Subroutine */ void crot_(integer *, complex *, integer *, 
	    complex *, integer *, real *, complex *);
    integer ilst;
    extern /* Subroutine */ void cgemm_(char *, char *, integer *, integer *, 
	    integer *, complex *, complex *, integer *, complex *, integer *, 
	    complex *, complex *, integer *);
    integer kwbot;
    real tempr;
    complex mktmp;
    integer kwtop;
    extern /* Subroutine */ void claqz0_(char *, char *, char *, integer *, 
	    integer *, integer *, complex *, integer *, complex *, integer *, 
	    complex *, complex *, complex *, integer *, complex *, integer *, 
	    complex *, integer *, real *, integer *, integer *);
    integer qz_small_info__;
    extern /* Subroutine */ void claqz1_(logical *, logical *, integer *, 
	    integer *, integer *, integer *, complex *, integer *, complex *, 
	    integer *, integer *, integer *, complex *, integer *, integer *, 
	    integer *, complex *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */ void clacpy_(char *, integer *, integer *, complex 
	    *, integer *, complex *, integer *), claset_(char *, 
	    integer *, integer *, complex *, complex *, complex *, integer *);
    real safmin;
    extern /* Subroutine */ void xerbla_(char *, integer *, ftnlen);
    real safmax;
    extern /* Subroutine */ void ctgexc_(logical *, logical *, integer *, 
	    complex *, integer *, complex *, integer *, complex *, integer *, 
	    complex *, integer *, integer *, integer *, integer *), clartg_(
	    complex *, complex *, real *, complex *, complex *);
    integer istopm;
    real smlnum;
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
    --alpha;
    --beta;
    qc_dim1 = *ldqc;
    qc_offset = 1 + qc_dim1;
    qc -= qc_offset;
    zc_dim1 = *ldzc;
    zc_offset = 1 + zc_dim1;
    zc -= zc_offset;
    --work;
    --rwork;

    /* Function Body */
    *info = 0;
/*     Set up deflation window */
/* Computing MIN */
    i__1 = *nw, i__2 = *ihi - *ilo + 1;
    jw = f2cmin(i__1,i__2);
    kwtop = *ihi - jw + 1;
    if (kwtop == *ilo) {
	s.r = 0.f, s.i = 0.f;
    } else {
	i__1 = kwtop + (kwtop - 1) * a_dim1;
	s.r = a[i__1].r, s.i = a[i__1].i;
    }
/*     Determine required workspace */
    ifst = 1;
    ilst = jw;
    i__1 = *rec + 1;
    claqz0_("S", "V", "V", &jw, &c__1, &jw, &a[kwtop + kwtop * a_dim1], lda, &
	    b[kwtop + kwtop * b_dim1], ldb, &alpha[1], &beta[1], &qc[
	    qc_offset], ldqc, &zc[zc_offset], ldzc, &work[1], &c_n1, &rwork[1]
	    , &i__1, &qz_small_info__);
/* Computing 2nd power */
    i__1 = jw;
    lworkreq = (integer) work[1].r + (i__1 * i__1 << 1);
/* Computing MAX */
/* Computing 2nd power */
    i__3 = *nw;
    i__1 = lworkreq, i__2 = *n * *nw, i__1 = f2cmax(i__1,i__2), i__2 = (i__3 * 
	    i__3 << 1) + *n;
    lworkreq = f2cmax(i__1,i__2);
    if (*lwork == -1) {
/*        workspace query, quick return */
	q__1.r = (real) lworkreq, q__1.i = 0.f;
	work[1].r = q__1.r, work[1].i = q__1.i;
	return;
    } else if (*lwork < lworkreq) {
	*info = -25;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CLAQZ2", &i__1, (ftnlen)6);
	return;
    }
/*     Get machine constants */
    safmin = slamch_("SAFE MINIMUM");
    safmax = 1.f / safmin;
    ulp = slamch_("PRECISION");
    smlnum = safmin * ((real) (*n) / ulp);
    if (*ihi == kwtop) {
/*        1 by 1 deflation window, just try a regular deflation */
	i__1 = kwtop;
	i__2 = kwtop + kwtop * a_dim1;
	alpha[i__1].r = a[i__2].r, alpha[i__1].i = a[i__2].i;
	i__1 = kwtop;
	i__2 = kwtop + kwtop * b_dim1;
	beta[i__1].r = b[i__2].r, beta[i__1].i = b[i__2].i;
	*ns = 1;
	*nd = 0;
/* Computing MAX */
	r__1 = smlnum, r__2 = ulp * c_abs(&a[kwtop + kwtop * a_dim1]);
	if (c_abs(&s) <= f2cmax(r__1,r__2)) {
	    *ns = 0;
	    *nd = 1;
	    if (kwtop > *ilo) {
		i__1 = kwtop + (kwtop - 1) * a_dim1;
		a[i__1].r = 0.f, a[i__1].i = 0.f;
	    }
	}
    }
/*     Store window in case of convergence failure */
    clacpy_("ALL", &jw, &jw, &a[kwtop + kwtop * a_dim1], lda, &work[1], &jw);
/* Computing 2nd power */
    i__1 = jw;
    clacpy_("ALL", &jw, &jw, &b[kwtop + kwtop * b_dim1], ldb, &work[i__1 * 
	    i__1 + 1], &jw);
/*     Transform window to real schur form */
    claset_("FULL", &jw, &jw, &c_b1, &c_b2, &qc[qc_offset], ldqc);
    claset_("FULL", &jw, &jw, &c_b1, &c_b2, &zc[zc_offset], ldzc);
/* Computing 2nd power */
    i__1 = jw;
/* Computing 2nd power */
    i__3 = jw;
    i__2 = *lwork - (i__3 * i__3 << 1);
    i__4 = *rec + 1;
    claqz0_("S", "V", "V", &jw, &c__1, &jw, &a[kwtop + kwtop * a_dim1], lda, &
	    b[kwtop + kwtop * b_dim1], ldb, &alpha[1], &beta[1], &qc[
	    qc_offset], ldqc, &zc[zc_offset], ldzc, &work[(i__1 * i__1 << 1) 
	    + 1], &i__2, &rwork[1], &i__4, &qz_small_info__);
    if (qz_small_info__ != 0) {
/*        Convergence failure, restore the window and exit */
	*nd = 0;
	*ns = jw - qz_small_info__;
	clacpy_("ALL", &jw, &jw, &work[1], &jw, &a[kwtop + kwtop * a_dim1], 
		lda);
/* Computing 2nd power */
	i__1 = jw;
	clacpy_("ALL", &jw, &jw, &work[i__1 * i__1 + 1], &jw, &b[kwtop + 
		kwtop * b_dim1], ldb);
	return;
    }
/*     Deflation detection loop */
    if (kwtop == *ilo || (s.r == 0.f && s.i == 0.f)) {
	kwbot = kwtop - 1;
    } else {
	kwbot = *ihi;
	k = 1;
	k2 = 1;
	while(k <= jw) {
/*              Try to deflate eigenvalue */
	    tempr = c_abs(&a[kwbot + kwbot * a_dim1]);
	    if (tempr == 0.f) {
		tempr = c_abs(&s);
	    }
	    i__1 = (kwbot - kwtop + 1) * qc_dim1 + 1;
	    q__1.r = s.r * qc[i__1].r - s.i * qc[i__1].i, q__1.i = s.r * qc[
		    i__1].i + s.i * qc[i__1].r;
/* Computing MAX */
	    r__1 = ulp * tempr;
	    if (c_abs(&q__1) <= f2cmax(r__1,smlnum)) {
/*                 Deflatable */
		--kwbot;
	    } else {
/*                 Not deflatable, move out of the way */
		ifst = kwbot - kwtop + 1;
		ilst = k2;
		ctgexc_(&c_true, &c_true, &jw, &a[kwtop + kwtop * a_dim1], 
			lda, &b[kwtop + kwtop * b_dim1], ldb, &qc[qc_offset], 
			ldqc, &zc[zc_offset], ldzc, &ifst, &ilst, &
			ctgexc_info__);
		++k2;
	    }
	    ++k;
	}
    }
/*     Store eigenvalues */
    *nd = *ihi - kwbot;
    *ns = jw - *nd;
    k = kwtop;
    while(k <= *ihi) {
	i__1 = k;
	i__2 = k + k * a_dim1;
	alpha[i__1].r = a[i__2].r, alpha[i__1].i = a[i__2].i;
	i__1 = k;
	i__2 = k + k * b_dim1;
	beta[i__1].r = b[i__2].r, beta[i__1].i = b[i__2].i;
	++k;
    }
    if (kwtop != *ilo && (s.r != 0.f || s.i != 0.f)) {
/*        Reflect spike back, this will create optimally packed bulges */
/*         A( KWTOP:KWBOT, KWTOP-1 ) = A( KWTOP, KWTOP-1 ) *CONJG( QC( 1, */
/*     $      1:JW-ND ) ) */
	i__1 = jw - *nd;
	for (imk = 1; imk <= i__1; ++imk) {
	    i__2 = kwtop + (kwtop - 1) * a_dim1;
	    r_cnjg(&q__2, &qc[imk * qc_dim1 + 1]);
	    q__1.r = a[i__2].r * q__2.r - a[i__2].i * q__2.i, q__1.i = a[i__2]
		    .r * q__2.i + a[i__2].i * q__2.r;
	    mktmp.r = q__1.r, mktmp.i = q__1.i;
	}
	i__1 = kwbot;
	for (imk = kwtop; imk <= i__1; ++imk) {
	    i__2 = imk + (kwtop - 1) * a_dim1;
	    a[i__2].r = mktmp.r, a[i__2].i = mktmp.i;
	}
	i__1 = kwtop;
	for (k = kwbot - 1; k >= i__1; --k) {
	    clartg_(&a[k + (kwtop - 1) * a_dim1], &a[k + 1 + (kwtop - 1) * 
		    a_dim1], &c1, &s1, &temp);
	    i__2 = k + (kwtop - 1) * a_dim1;
	    a[i__2].r = temp.r, a[i__2].i = temp.i;
	    i__2 = k + 1 + (kwtop - 1) * a_dim1;
	    a[i__2].r = 0.f, a[i__2].i = 0.f;
/* Computing MAX */
	    i__2 = kwtop, i__3 = k - 1;
	    k2 = f2cmax(i__2,i__3);
	    i__2 = *ihi - k2 + 1;
	    crot_(&i__2, &a[k + k2 * a_dim1], lda, &a[k + 1 + k2 * a_dim1], 
		    lda, &c1, &s1);
	    i__2 = *ihi - (k - 1) + 1;
	    crot_(&i__2, &b[k + (k - 1) * b_dim1], ldb, &b[k + 1 + (k - 1) * 
		    b_dim1], ldb, &c1, &s1);
	    r_cnjg(&q__1, &s1);
	    crot_(&jw, &qc[(k - kwtop + 1) * qc_dim1 + 1], &c__1, &qc[(k + 1 
		    - kwtop + 1) * qc_dim1 + 1], &c__1, &c1, &q__1);
	}
/*        Chase bulges down */
	istartm = kwtop;
	istopm = *ihi;
	k = kwbot - 1;
	while(k >= kwtop) {
/*           Move bulge down and remove it */
	    i__1 = kwbot - 1;
	    for (k2 = k; k2 <= i__1; ++k2) {
		i__2 = kwtop + jw - 1;
		claqz1_(&c_true, &c_true, &k2, &kwtop, &i__2, &kwbot, &a[
			a_offset], lda, &b[b_offset], ldb, &jw, &kwtop, &qc[
			qc_offset], ldqc, &jw, &kwtop, &zc[zc_offset], ldzc);
	    }
	    --k;
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
	cgemm_("C", "N", &jw, &i__1, &jw, &c_b2, &qc[qc_offset], ldqc, &a[
		kwtop + (*ihi + 1) * a_dim1], lda, &c_b1, &work[1], &jw);
	i__1 = istopm - *ihi;
	clacpy_("ALL", &jw, &i__1, &work[1], &jw, &a[kwtop + (*ihi + 1) * 
		a_dim1], lda);
	i__1 = istopm - *ihi;
	cgemm_("C", "N", &jw, &i__1, &jw, &c_b2, &qc[qc_offset], ldqc, &b[
		kwtop + (*ihi + 1) * b_dim1], ldb, &c_b1, &work[1], &jw);
	i__1 = istopm - *ihi;
	clacpy_("ALL", &jw, &i__1, &work[1], &jw, &b[kwtop + (*ihi + 1) * 
		b_dim1], ldb);
    }
    if (*ilq) {
	cgemm_("N", "N", n, &jw, &jw, &c_b2, &q[kwtop * q_dim1 + 1], ldq, &qc[
		qc_offset], ldqc, &c_b1, &work[1], n);
	clacpy_("ALL", n, &jw, &work[1], n, &q[kwtop * q_dim1 + 1], ldq);
    }
    if (kwtop - 1 - istartm + 1 > 0) {
	i__1 = kwtop - istartm;
	i__2 = kwtop - istartm;
	cgemm_("N", "N", &i__1, &jw, &jw, &c_b2, &a[istartm + kwtop * a_dim1],
		 lda, &zc[zc_offset], ldzc, &c_b1, &work[1], &i__2);
	i__1 = kwtop - istartm;
	i__2 = kwtop - istartm;
	clacpy_("ALL", &i__1, &jw, &work[1], &i__2, &a[istartm + kwtop * 
		a_dim1], lda);
	i__1 = kwtop - istartm;
	i__2 = kwtop - istartm;
	cgemm_("N", "N", &i__1, &jw, &jw, &c_b2, &b[istartm + kwtop * b_dim1],
		 ldb, &zc[zc_offset], ldzc, &c_b1, &work[1], &i__2);
	i__1 = kwtop - istartm;
	i__2 = kwtop - istartm;
	clacpy_("ALL", &i__1, &jw, &work[1], &i__2, &b[istartm + kwtop * 
		b_dim1], ldb);
    }
    if (*ilz) {
	cgemm_("N", "N", n, &jw, &jw, &c_b2, &z__[kwtop * z_dim1 + 1], ldz, &
		zc[zc_offset], ldzc, &c_b1, &work[1], n)
		;
	clacpy_("ALL", n, &jw, &work[1], n, &z__[kwtop * z_dim1 + 1], ldz);
    }
    return;
} /* claqz2_ */

