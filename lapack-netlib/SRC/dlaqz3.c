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
#define pow_dd(ap, bp) ( pow(*(ap), *(bp)))
#define pow_si(B,E) spow_ui(*(B),*(E))
#define pow_ri(B,E) spow_ui(*(B),*(E))
#define pow_di(B,E) dpow_ui(*(B),*(E))
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

static logical c_true = TRUE_;
static integer c_n1 = -1;
static integer c__1 = 1;
static doublereal c_b16 = 0.;
static doublereal c_b17 = 1.;

/* Subroutine */ void dlaqz3_(logical *ilschur, logical *ilq, logical *ilz, 
	integer *n, integer *ilo, integer *ihi, integer *nw, doublereal *a, 
	integer *lda, doublereal *b, integer *ldb, doublereal *q, integer *
	ldq, doublereal *z__, integer *ldz, integer *ns, integer *nd, 
	doublereal *alphar, doublereal *alphai, doublereal *beta, doublereal *
	qc, integer *ldqc, doublereal *zc, integer *ldzc, doublereal *work, 
	integer *lwork, integer *rec, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, qc_dim1, qc_offset, zc_dim1, zc_offset, i__1, i__2, 
	    i__3, i__4;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    integer lworkreq, k;
    doublereal s, c1;
    integer k2;
    doublereal s1;
    integer jw, imk;
    doublereal ulp;
    integer dtgexc_info__, ifst;
    doublereal temp;
    extern /* Subroutine */ void drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    integer ilst;
    extern /* Subroutine */ void dlag2_(doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *), dgemm_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *);
    logical bulge;
    integer kwbot;
    doublereal mktmp;
    integer kwtop, qz_small_info__;
    extern /* Subroutine */ void dlaqz0_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *), dlaqz2_(logical *, 
	    logical *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *), dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *);
    extern /* Subroutine */ void dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    doublereal safmin;
    extern /* Subroutine */ void xerbla_(char *, integer *);
    doublereal safmax;
    extern /* Subroutine */ void dtgexc_(logical *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *), dlartg_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    integer istopm;
    doublereal smlnum;
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
	s = 0.;
    } else {
	s = a[kwtop + (kwtop - 1) * a_dim1];
    }
/*     Determine required workspace */
    ifst = 1;
    ilst = jw;
    dtgexc_(&c_true, &c_true, &jw, &a[a_offset], lda, &b[b_offset], ldb, &qc[
	    qc_offset], ldqc, &zc[zc_offset], ldzc, &ifst, &ilst, &work[1], &
	    c_n1, &dtgexc_info__);
    lworkreq = (integer) work[1];
    i__1 = *rec + 1;
    dlaqz0_("S", "V", "V", &jw, &c__1, &jw, &a[kwtop + kwtop * a_dim1], lda, &
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
	work[1] = (doublereal) lworkreq;
	return;
    } else if (*lwork < lworkreq) {
	*info = -26;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DLAQZ3", &i__1);
	return;
    }
/*     Get machine constants */
    safmin = dlamch_("SAFE MINIMUM");
    safmax = 1. / safmin;
    dlabad_(&safmin, &safmax);
    ulp = dlamch_("PRECISION");
    smlnum = safmin * ((doublereal) (*n) / ulp);
    if (*ihi == kwtop) {
/*        1 by 1 deflation window, just try a regular deflation */
	alphar[kwtop] = a[kwtop + kwtop * a_dim1];
	alphai[kwtop] = 0.;
	beta[kwtop] = b[kwtop + kwtop * b_dim1];
	*ns = 1;
	*nd = 0;
/* Computing MAX */
	d__2 = smlnum, d__3 = ulp * (d__1 = a[kwtop + kwtop * a_dim1], abs(
		d__1));
	if (abs(s) <= f2cmax(d__2,d__3)) {
	    *ns = 0;
	    *nd = 1;
	    if (kwtop > *ilo) {
		a[kwtop + (kwtop - 1) * a_dim1] = 0.;
	    }
	}
    }
/*     Store window in case of convergence failure */
    dlacpy_("ALL", &jw, &jw, &a[kwtop + kwtop * a_dim1], lda, &work[1], &jw);
/* Computing 2nd power */
    i__1 = jw;
    dlacpy_("ALL", &jw, &jw, &b[kwtop + kwtop * b_dim1], ldb, &work[i__1 * 
	    i__1 + 1], &jw);
/*     Transform window to real schur form */
    dlaset_("FULL", &jw, &jw, &c_b16, &c_b17, &qc[qc_offset], ldqc)
	    ;
    dlaset_("FULL", &jw, &jw, &c_b16, &c_b17, &zc[zc_offset], ldzc)
	    ;
/* Computing 2nd power */
    i__1 = jw;
/* Computing 2nd power */
    i__3 = jw;
    i__2 = *lwork - (i__3 * i__3 << 1);
    i__4 = *rec + 1;
    dlaqz0_("S", "V", "V", &jw, &c__1, &jw, &a[kwtop + kwtop * a_dim1], lda, &
	    b[kwtop + kwtop * b_dim1], ldb, &alphar[1], &alphai[1], &beta[1], 
	    &qc[qc_offset], ldqc, &zc[zc_offset], ldzc, &work[(i__1 * i__1 << 
	    1) + 1], &i__2, &i__4, &qz_small_info__);
    if (qz_small_info__ != 0) {
/*        Convergence failure, restore the window and exit */
	*nd = 0;
	*ns = jw - qz_small_info__;
	dlacpy_("ALL", &jw, &jw, &work[1], &jw, &a[kwtop + kwtop * a_dim1], 
		lda);
/* Computing 2nd power */
	i__1 = jw;
	dlacpy_("ALL", &jw, &jw, &work[i__1 * i__1 + 1], &jw, &b[kwtop + 
		kwtop * b_dim1], ldb);
	return;
    }
/*     Deflation detection loop */
    if (kwtop == *ilo || s == 0.) {
	kwbot = kwtop - 1;
    } else {
	kwbot = *ihi;
	k = 1;
	k2 = 1;
	while(k <= jw) {
	    bulge = FALSE_;
	    if (kwbot - kwtop + 1 >= 2) {
		bulge = a[kwbot + (kwbot - 1) * a_dim1] != 0.;
	    }
	    if (bulge) {
/*              Try to deflate complex conjugate eigenvalue pair */
		temp = (d__3 = a[kwbot + kwbot * a_dim1], abs(d__3)) + sqrt((
			d__1 = a[kwbot + (kwbot - 1) * a_dim1], abs(d__1))) * 
			sqrt((d__2 = a[kwbot - 1 + kwbot * a_dim1], abs(d__2))
			);
		if (temp == 0.) {
		    temp = abs(s);
		}
/* Computing MAX */
		d__3 = (d__1 = s * qc[(kwbot - kwtop) * qc_dim1 + 1], abs(
			d__1)), d__4 = (d__2 = s * qc[(kwbot - kwtop + 1) * 
			qc_dim1 + 1], abs(d__2));
/* Computing MAX */
		d__5 = smlnum, d__6 = ulp * temp;
		if (f2cmax(d__3,d__4) <= f2cmax(d__5,d__6)) {
/*                 Deflatable */
		    kwbot += -2;
		} else {
/*                 Not deflatable, move out of the way */
		    ifst = kwbot - kwtop + 1;
		    ilst = k2;
		    dtgexc_(&c_true, &c_true, &jw, &a[kwtop + kwtop * a_dim1],
			     lda, &b[kwtop + kwtop * b_dim1], ldb, &qc[
			    qc_offset], ldqc, &zc[zc_offset], ldzc, &ifst, &
			    ilst, &work[1], lwork, &dtgexc_info__);
		    k2 += 2;
		}
		k += 2;
	    } else {
/*              Try to deflate real eigenvalue */
		temp = (d__1 = a[kwbot + kwbot * a_dim1], abs(d__1));
		if (temp == 0.) {
		    temp = abs(s);
		}
/* Computing MAX */
		d__2 = ulp * temp;
		if ((d__1 = s * qc[(kwbot - kwtop + 1) * qc_dim1 + 1], abs(
			d__1)) <= f2cmax(d__2,smlnum)) {
/*                 Deflatable */
		    --kwbot;
		} else {
/*                 Not deflatable, move out of the way */
		    ifst = kwbot - kwtop + 1;
		    ilst = k2;
		    dtgexc_(&c_true, &c_true, &jw, &a[kwtop + kwtop * a_dim1],
			     lda, &b[kwtop + kwtop * b_dim1], ldb, &qc[
			    qc_offset], ldqc, &zc[zc_offset], ldzc, &ifst, &
			    ilst, &work[1], lwork, &dtgexc_info__);
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
	    if (a[k + 1 + k * a_dim1] != 0.) {
		bulge = TRUE_;
	    }
	}
	if (bulge) {
/*           2x2 eigenvalue block */
	    dlag2_(&a[k + k * a_dim1], lda, &b[k + k * b_dim1], ldb, &safmin, 
		    &beta[k], &beta[k + 1], &alphar[k], &alphar[k + 1], &
		    alphai[k]);
	    alphai[k + 1] = -alphai[k];
	    k += 2;
	} else {
/*           1x1 eigenvalue block */
	    alphar[k] = a[k + k * a_dim1];
	    alphai[k] = 0.;
	    beta[k] = b[k + k * b_dim1];
	    ++k;
	}
    }
    if (kwtop != *ilo && s != 0.) {
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
	    dlartg_(&a[k + (kwtop - 1) * a_dim1], &a[k + 1 + (kwtop - 1) * 
		    a_dim1], &c1, &s1, &temp);
	    a[k + (kwtop - 1) * a_dim1] = temp;
	    a[k + 1 + (kwtop - 1) * a_dim1] = 0.;
/* Computing MAX */
	    i__2 = kwtop, i__3 = k - 1;
	    k2 = f2cmax(i__2,i__3);
	    i__2 = *ihi - k2 + 1;
	    drot_(&i__2, &a[k + k2 * a_dim1], lda, &a[k + 1 + k2 * a_dim1], 
		    lda, &c1, &s1);
	    i__2 = *ihi - (k - 1) + 1;
	    drot_(&i__2, &b[k + (k - 1) * b_dim1], ldb, &b[k + 1 + (k - 1) * 
		    b_dim1], ldb, &c1, &s1);
	    drot_(&jw, &qc[(k - kwtop + 1) * qc_dim1 + 1], &c__1, &qc[(k + 1 
		    - kwtop + 1) * qc_dim1 + 1], &c__1, &c1, &s1);
	}
/*        Chase bulges down */
	istartm = kwtop;
	istopm = *ihi;
	k = kwbot - 1;
	while(k >= kwtop) {
	    if (k >= kwtop + 1 && a[k + 1 + (k - 1) * a_dim1] != 0.) {
/*              Move double pole block down and remove it */
		i__1 = kwbot - 2;
		for (k2 = k - 1; k2 <= i__1; ++k2) {
		    i__2 = kwtop + jw - 1;
		    dlaqz2_(&c_true, &c_true, &k2, &kwtop, &i__2, &kwbot, &a[
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
		    dlartg_(&b[k2 + 1 + (k2 + 1) * b_dim1], &b[k2 + 1 + k2 * 
			    b_dim1], &c1, &s1, &temp);
		    b[k2 + 1 + (k2 + 1) * b_dim1] = temp;
		    b[k2 + 1 + k2 * b_dim1] = 0.;
		    i__2 = k2 + 2 - istartm + 1;
		    drot_(&i__2, &a[istartm + (k2 + 1) * a_dim1], &c__1, &a[
			    istartm + k2 * a_dim1], &c__1, &c1, &s1);
		    i__2 = k2 - istartm + 1;
		    drot_(&i__2, &b[istartm + (k2 + 1) * b_dim1], &c__1, &b[
			    istartm + k2 * b_dim1], &c__1, &c1, &s1);
		    drot_(&jw, &zc[(k2 + 1 - kwtop + 1) * zc_dim1 + 1], &c__1,
			     &zc[(k2 - kwtop + 1) * zc_dim1 + 1], &c__1, &c1, 
			    &s1);
		    dlartg_(&a[k2 + 1 + k2 * a_dim1], &a[k2 + 2 + k2 * a_dim1]
			    , &c1, &s1, &temp);
		    a[k2 + 1 + k2 * a_dim1] = temp;
		    a[k2 + 2 + k2 * a_dim1] = 0.;
		    i__2 = istopm - k2;
		    drot_(&i__2, &a[k2 + 1 + (k2 + 1) * a_dim1], lda, &a[k2 + 
			    2 + (k2 + 1) * a_dim1], lda, &c1, &s1);
		    i__2 = istopm - k2;
		    drot_(&i__2, &b[k2 + 1 + (k2 + 1) * b_dim1], ldb, &b[k2 + 
			    2 + (k2 + 1) * b_dim1], ldb, &c1, &s1);
		    drot_(&jw, &qc[(k2 + 1 - kwtop + 1) * qc_dim1 + 1], &c__1,
			     &qc[(k2 + 2 - kwtop + 1) * qc_dim1 + 1], &c__1, &
			    c1, &s1);
		}
/*              Remove the shift */
		dlartg_(&b[kwbot + kwbot * b_dim1], &b[kwbot + (kwbot - 1) * 
			b_dim1], &c1, &s1, &temp);
		b[kwbot + kwbot * b_dim1] = temp;
		b[kwbot + (kwbot - 1) * b_dim1] = 0.;
		i__1 = kwbot - istartm;
		drot_(&i__1, &b[istartm + kwbot * b_dim1], &c__1, &b[istartm 
			+ (kwbot - 1) * b_dim1], &c__1, &c1, &s1);
		i__1 = kwbot - istartm + 1;
		drot_(&i__1, &a[istartm + kwbot * a_dim1], &c__1, &a[istartm 
			+ (kwbot - 1) * a_dim1], &c__1, &c1, &s1);
		drot_(&jw, &zc[(kwbot - kwtop + 1) * zc_dim1 + 1], &c__1, &zc[
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
	dgemm_("T", "N", &jw, &i__1, &jw, &c_b17, &qc[qc_offset], ldqc, &a[
		kwtop + (*ihi + 1) * a_dim1], lda, &c_b16, &work[1], &jw);
	i__1 = istopm - *ihi;
	dlacpy_("ALL", &jw, &i__1, &work[1], &jw, &a[kwtop + (*ihi + 1) * 
		a_dim1], lda);
	i__1 = istopm - *ihi;
	dgemm_("T", "N", &jw, &i__1, &jw, &c_b17, &qc[qc_offset], ldqc, &b[
		kwtop + (*ihi + 1) * b_dim1], ldb, &c_b16, &work[1], &jw);
	i__1 = istopm - *ihi;
	dlacpy_("ALL", &jw, &i__1, &work[1], &jw, &b[kwtop + (*ihi + 1) * 
		b_dim1], ldb);
    }
    if (*ilq) {
	dgemm_("N", "N", n, &jw, &jw, &c_b17, &q[kwtop * q_dim1 + 1], ldq, &
		qc[qc_offset], ldqc, &c_b16, &work[1], n);
	dlacpy_("ALL", n, &jw, &work[1], n, &q[kwtop * q_dim1 + 1], ldq);
    }
    if (kwtop - 1 - istartm + 1 > 0) {
	i__1 = kwtop - istartm;
	i__2 = kwtop - istartm;
	dgemm_("N", "N", &i__1, &jw, &jw, &c_b17, &a[istartm + kwtop * a_dim1]
		, lda, &zc[zc_offset], ldzc, &c_b16, &work[1], &i__2);
	i__1 = kwtop - istartm;
	i__2 = kwtop - istartm;
	dlacpy_("ALL", &i__1, &jw, &work[1], &i__2, &a[istartm + kwtop * 
		a_dim1], lda);
	i__1 = kwtop - istartm;
	i__2 = kwtop - istartm;
	dgemm_("N", "N", &i__1, &jw, &jw, &c_b17, &b[istartm + kwtop * b_dim1]
		, ldb, &zc[zc_offset], ldzc, &c_b16, &work[1], &i__2);
	i__1 = kwtop - istartm;
	i__2 = kwtop - istartm;
	dlacpy_("ALL", &i__1, &jw, &work[1], &i__2, &b[istartm + kwtop * 
		b_dim1], ldb);
    }
    if (*ilz) {
	dgemm_("N", "N", n, &jw, &jw, &c_b17, &z__[kwtop * z_dim1 + 1], ldz, &
		zc[zc_offset], ldzc, &c_b16, &work[1], n);
	dlacpy_("ALL", n, &jw, &work[1], n, &z__[kwtop * z_dim1 + 1], ldz);
    }
    return;
} /* dlaqz3_ */

