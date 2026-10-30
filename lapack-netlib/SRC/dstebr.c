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
typedef blasint logical;

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
#define mycycle() continue;
#define myceiling(w) {ceil(w)}
#define myhuge(w) {HUGE_VAL}
//#define mymaxloc_(w,s,e,n) {if (sizeof(*(w)) == sizeof(double)) dmaxloc_((w),*(s),*(e),n); else dmaxloc_((w),*(s),*(e),n);}
#define mymaxloc(w,s,e,n) {dmaxloc_(w,*(s),*(e),n)}

/*  -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


/* Table of constant values */

static integer c__9 = 9;
static integer c__0 = 0;
static doublereal c_b19 = 1.;
static integer c__1 = 1;

/* > \brief \b DSTEBR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSTEBR( N, D, E, WORK, LWORK, IWORK, LIWORK, INFO ) */

/*       INTEGER            INFO, LIWORK, LWORK, N */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   D( * ), E( * ), WORK( * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSTEBR computes all eigenvalues of a symmetric tridiagonal matrix */
/* > using a values-only divide-and-conquer method that propagates only */
/* > the first and last boundary rows of each local eigenvector block */
/* > (STEBR = Symmetric tridiagonal Eigenvalues, Boundary-Row DC). */
/* > */
/* > Unlike DSTEDC with COMPZ = 'N', DSTEBR does not store full secular */
/* > eigenvector blocks across merge levels for later replay.  It is a */
/* > distinct routine and does not replace DSTEDC or DSTERF. */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the n diagonal elements of the tridiagonal matrix. */
/* >          On exit, if INFO = 0, the eigenvalues in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* >          matrix. */
/* >          On exit, E has been destroyed. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* >          If N <= SMLSIZ, where SMLSIZ is returned by ILAENV and is */
/* >          typically about 25, LWORK must be at least 1. */
/* >          If N > SMLSIZ, LWORK must be at least */
/* >          MAX( 16*N, 2*N + SMLSIZ**2 + 2*SMLSIZ - 2 ). */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the optimal sizes of the WORK and */
/* >          IWORK arrays, returns these values as the first entries of */
/* >          the WORK and IWORK arrays, and no error message related to */
/* >          LWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* >          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* >          LIWORK is INTEGER */
/* >          The dimension of the array IWORK. */
/* >          If N <= SMLSIZ, LIWORK must be at least 1. */
/* >          If N > SMLSIZ, LIWORK must be at least 7*N. */
/* > */
/* >          If LIWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the optimal sizes of the WORK and */
/* >          IWORK arrays, returns these values as the first entries of */
/* >          the WORK and IWORK arrays, and no error message related to */
/* >          LWORK or LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  The algorithm failed to compute an eigenvalue while */
/* >                working on a submatrix lying in rows and columns */
/* >                INFO/(N+1) through MOD(INFO,N+1). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup stebr */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Ruiyi Zhan and Shaoshuai Zhang, University of Electronic Science */
/* > and Technology of China */

/*  ===================================================================== */
/* Subroutine */ void dstebr_(integer *n, doublereal *d__, doublereal *e, 
	doublereal *work, integer *lwork, integer *iwork, integer *liwork, 
	integer *info)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    extern /* Subroutine */ void dlaed0_br_(integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *);
    integer m;
    doublereal eps, tiny;
    integer lwmin, start;
    extern doublereal dlamch_(char *);
    extern /* Subroutine */ void dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    integer finish;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *);
    extern /* Subroutine */ void dsterf_(integer *, doublereal *, doublereal *,
	     integer *), dlasrt_(char *, integer *, doublereal *, integer *);
    integer liwmin;
    doublereal orgnrm;
    logical lquery;
    integer smlsiz;


/*  -- LAPACK computational routine -- */


/*  ===================================================================== */


    /* Parameter adjustments */
    --iwork;
    --work;
    --e;
    --d__;

    /* Function Body */
    *info = 0;
    lquery = *lwork == -1 || *liwork == -1;

    if (*n < 0) {
	*info = -1;
    }

    if (*info == 0) {
	smlsiz = ilaenv_(&c__9, "DSTEBR", " ", &c__0, &c__0, &c__0, &c__0, (
		ftnlen)6, (ftnlen)1);
	if (*n <= smlsiz) {
	    lwmin = 1;
	    liwmin = 1;
	} else {

/*           BLO/BHI: 2*N; merge scratch: 14*N (see DLAED7_BR). */
/*           Leaf solves need 2*N + SMLSIZ**2 + 2*SMLSIZ - 2: the */
/*           SMLSIZ-by-SMLSIZ leaf eigenvector block plus the */
/*           DSTEQR work array. */

/* Computing MAX */
	    i__1 = *n << 4, i__2 = (*n << 1) + smlsiz * smlsiz + (smlsiz << 1)
		     - 2;
	    lwmin = f2cmax(i__1,i__2);
/* Computing MAX */
	    i__1 = 1, i__2 = *n * 7;
	    liwmin = f2cmax(i__1,i__2);
	}
	work[1] = (doublereal) lwmin;
	iwork[1] = liwmin;

	if (*lwork < lwmin && ! lquery) {
	    *info = -5;
	} else if (*liwork < liwmin && ! lquery) {
	    *info = -7;
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSTEBR", &i__1, (ftnlen)6);
	return;
    } else if (lquery) {
	return;
    }

    if (*n == 0) {
	return;
    }
    if (*n == 1) {
	return;
    }

    if (*n <= smlsiz) {
	dsterf_(n, &d__[1], &e[1], info);
	goto L50;
    }

    orgnrm = dlanst_("M", n, &d__[1], &e[1]);
    if (orgnrm <= 0.) {
	goto L50;
    }

    eps = dlamch_("Epsilon");
    start = 1;

L10:
    if (start <= *n) {
	finish = start;
L20:
	if (finish < *n) {
	    tiny = eps * sqrt((d__1 = d__[finish], abs(d__1))) * sqrt((d__2 = 
		    d__[finish + 1], abs(d__2)));
	    if ((d__1 = e[finish], abs(d__1)) > tiny) {
		++finish;
		goto L20;
	    }
	}

	m = finish - start + 1;
	if (m == 1) {
	    start = finish + 1;
	    goto L10;
	}

	if (m > smlsiz) {
	    orgnrm = dlanst_("M", &m, &d__[start], &e[start]);
	    dlascl_("G", &c__0, &c__0, &orgnrm, &c_b19, &m, &c__1, &d__[start]
		    , &m, info);
	    i__1 = m - 1;
	    i__2 = m - 1;
	    dlascl_("G", &c__0, &c__0, &orgnrm, &c_b19, &i__1, &c__1, &e[
		    start], &i__2, info);

	    dlaed0_br_(&m, &d__[start], &e[start], &work[1], &iwork[1], info)
		    ;
	    if (*info != 0) {
		*info = (*info / (m + 1) + start - 1) * (*n + 1) + *info % (m 
			+ 1) + start - 1;
		goto L50;
	    }

	    dlascl_("G", &c__0, &c__0, &c_b19, &orgnrm, &m, &c__1, &d__[start]
		    , &m, info);
	} else {
	    dsterf_(&m, &d__[start], &e[start], info);
	    if (*info != 0) {
		*info = start * (*n + 1) + finish;
		goto L50;
	    }
	}

	start = finish + 1;
	goto L10;
    }

    dlasrt_("I", n, &d__[1], info);

L50:
    work[1] = (doublereal) lwmin;
    iwork[1] = liwmin;
    return;

/*     End of DSTEBR */

} /* dstebr_ */

