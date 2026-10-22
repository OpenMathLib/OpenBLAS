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

/* procedure parameter types for -A and -C++ */


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
#endif
/*  -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


/* > \brief \b SLARRB_TGK bisection refinement of eigenvalues of an LDL^T representation (Golub-Kahan variant 
of SLARRB); also SLARRBTGK_STURM4 (four-lane Sturm count) and SLARRBTGK_STURM_SLOW (reference Sturm co
unt). */
/* > */
/* > \ingroup bdsvdmr3 */
/* > */
/* > \par Purpose: */
/* > ============= */
/* > */
/* > \verbatim */
/* > bisection refinement of eigenvalues of an LDL^T representation (Golub-Kahan variant of SLARRB); also SLAR
RBTGK_STURM4 (four-lane Sturm count) and SLARRBTGK_STURM_SLOW (reference Sturm count). */
/* > All routines in this file are auxiliary to SBDSVDMR3 and derive */
/* > from the LAPACK 3.0 (1999) MRRR kernels SLARRV, SLAR1V, SLARRB and */
/* > SLARRF rather than from the current SSTEMR kernels. */
/* > \endverbatim */

/* Subroutine */ void slarrb_tgk_(integer *n, real *d__, real *l, real *ld, 
	real *lld, integer *ifirst, integer *ilast, real *rtol1, real *rtol2, 
	integer *offset, real *w, real *wgap, real *werr, real *work, integer 
	*iwork, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1, r__2;

    /* Local variables */
    integer b, i__, j, k, p;
    real s;
    integer i1, i2, bi[4], nb, ii, kk;
    real fac, gap, mid;
    integer cnt;
    real eps, tmp, bmid[4];
    integer bcnt[4];
    real left;
    integer iter;
    real ovfl;
    integer nint, nexp, prev, next, nleft;
    real right, width;
    integer bnext[4];
    real dplus, error;
    extern /* Subroutine */ void slarrbtgk_sturm_slow_(integer *, real *, 
	    real *, real *, integer *);
    integer bafter;
    extern real slamch_(char *);
    integer bnrght[4], nright, olnint;


/*  -- LAPACK auxiliary routine -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     Derived from the LAPACK 3.0 SLARRB (1999); TGK adaptation 2026 */


/*  Purpose */
/*  ======= */

/*  Given the relatively robust representation(RRR) L D L^T, SLARRB */
/*  does "limited" bisection to refine the eigenvalues of L D L^T, */
/*  W( IFIRST-OFFSET ) thru' W( ILAST-OFFSET ), to more accuracy. Initial */
/*  guesses for these eigenvalues are input in W, the corresponding estimate */
/*  of the error in these guesses and their gaps are input in WERR */
/*  and WGAP, respectively. During bisection, intervals */
/*  [left, right] are maintained by storing their mid-points and */
/*  semi-widths in the arrays W and WERR respectively. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The order of the matrix. */

/*  D       (input) REAL array, dimension (N) */
/*          The N diagonal elements of the diagonal matrix D. */

/*  L       (input) REAL array, dimension (N-1) */
/*          The (N-1) subdiagonal elements of the unit bidiagonal matrix L. */

/*  LD      (input) REAL array, dimension (N-1) */
/*          The (N-1) elements L(i)*D(i). */

/*  LLD     (input) REAL array, dimension (N-1) */
/*          The (N-1) elements L(i)*L(i)*D(i). */

/*  IFIRST  (input) INTEGER */
/*          The index of the first eigenvalue to be computed. */

/*  ILAST   (input) INTEGER */
/*          The index of the last eigenvalue to be computed. */

/*  RTOL1   (input) REAL */
/*          Tolerance for the convergence of the bisection intervals. */
/*          An interval [LEFT,RIGHT] has converged if */
/*          RIGHT-LEFT.LT.MAX( RTOL1*GAP, RTOL2*ABS((LEFT+RIGHT)/2) ), */
/*          where GAP is the (estimated) distance to the nearest */
/*          eigenvalue. */

/*  RTOL2   (input) REAL */
/*          Tolerance for the convergence of the bisection intervals. */
/*          An interval [LEFT,RIGHT] has converged if */
/*          RIGHT-LEFT.LT.MAX( RTOL1*GAP, RTOL2*ABS((LEFT+RIGHT)/2) ), */
/*          where GAP is the (estimated) distance to the nearest */
/*          eigenvalue. */

/*  OFFSET  (input) INTEGER */
/*          Offset for the arrays W, WGAP and WERR, i.e., the IFIRST-OFFSET */
/*          thru' ILAST-OFFSET elements of these arrays are to be used. */

/*  W       (input/output) REAL array, dimension (N) */
/*          On input, W( IFIRST-OFFSET ) thru' W( ILAST-OFFSET ) are */
/*          estimates of the eigenvalues of L D L^T indexed IFIRST thru' ILAST. */
/*          On output, these estimates are refined. */

/*  WGAP    (input/output) REAL array, dimension (N-1) */
/*          On input, the (estimated) gaps between consecutive */
/*          eigenvalues of L D L^T, i.e., WGAP(I-OFFSET) is the gap between */
/*          eigenvalues I and I+1. Note that if IFIRST.EQ.ILAST */
/*          then WGAP(IFIRST-OFFSET) must be set to ZERO. */
/*          On output, these gaps are refined. */

/*  WERR    (input/output) REAL array, dimension (N) */
/*          On input, WERR( IFIRST-OFFSET ) thru' WERR( ILAST-OFFSET ) are */
/*          the errors in the estimates of the corresponding elements in W. */
/*          On output, these errors are refined. */

/*  WORK    (workspace) REAL array, dimension (2*N) */
/*          Workspace. */

/*  IWORK   (workspace) INTEGER array, dimension (2*N) */
/*          Workspace. */

/*  INFO    (output) INTEGER */
/*          = 0: successful exit. */
/*          = 1: a bracket expansion was nonfinite or made no progress. */
/*          = 2: a bracket expansion exceeded MAXEXP iterations. */
/*          = 3: interval refinement exceeded MAXREF passes. */

/*  Further Details */
/*  =============== */

/*  Based on contributions by */
/*     Inderjit Dhillon, University of Texas, Austin, USA */
/*     Osni Marques, LBNL/NERSC, USA */

/*  ===================================================================== */



/*     Check to see if any of the initial eigenvalue */
/*     estimates is acceptable. */

    /* Parameter adjustments */
    --iwork;
    --work;
    --werr;
    --wgap;
    --w;
    --lld;
    --ld;
    --l;
    --d__;

    /* Function Body */
    *info = 0;
    eps = slamch_("Precision");
    ovfl = slamch_("Overflow");
    i__1 = *n << 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iwork[i__] = 0;
/* L5: */
    }
    i1 = *ifirst;
    i2 = *ifirst;
    prev = 0;
    i__1 = *ilast;
    for (i__ = *ifirst; i__ <= i__1; ++i__) {
	ii = i__ - *offset;
	if (i__ == *ifirst) {
	    gap = wgap[ii];
	} else if (i__ == *ilast) {
	    gap = wgap[ii - 1];
	} else {
/* Computing MIN */
	    r__1 = wgap[ii - 1], r__2 = wgap[ii];
	    gap = f2cmin(r__1,r__2);
	}
	error = werr[ii];
	k = i__ << 1;
/*        IF( ERROR.LT.RTOL1*GAP ) THEN */
/*           WORK( K-1 ) = W( II ) - ERROR */
/*           WORK( K ) = W( II ) + ERROR */
/*           IWORK( K-1 ) = -1 */
/*           IF( I1.EQ.I ) THEN */
/*              I1 = I1 + 1 */
/*              PREV = I */
/*           END IF */
/*        ELSE */
	iwork[k - 1] = 1;
	i2 = i__;
/*        END IF */
/* L10: */
    }

/*     Initialize unconverged intervals in [ WORK(2*I-1), WORK(2*I) ]. */
/*     The Sturm Count, Count( WORK(2*I-1) ) is arranged to be I-1, while */
/*     Count( WORK(2*I) ) is stored in IWORK( 2*I ). The integer IWORK( 2*I-1 ) */
/*     for an unconverged interval is set to the index of the next unconverged */
/*     interval, and is -1 or 0 for a converged interval. Thus a linked */
/*     list of unconverged intervals is set up. */

    i__ = i1;
    nint = 0;
L30:
    if (i__ <= i2) {
	ii = i__ - *offset;
	if (iwork[(i__ << 1) - 1] == 1) {
	    fac = 1.f;
	    nexp = 0;
	    left = w[ii] - werr[ii];

/*           Do while( CNT(LEFT).GT.I-1 ) */

L40:
	    if (i__ > i1 && left <= right) {
		left = right;
		cnt = i__ - 1;
	    } else {
		s = -left;
		cnt = 0;
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    dplus = d__[j] + s;
		    s = s * lld[j] / dplus - left;
		    if (dplus < 0.f) {
			++cnt;
		    }
/* L50: */
		}
		dplus = d__[*n] + s;
		if (dplus < 0.f) {
		    ++cnt;
		}
		if (! (s > 0.f || s < 1.f)) {

/*                 Runs a slower version of the above loop if a NaN is detected */

		    cnt = 0;
		    s = -left;
		    i__1 = *n - 1;
		    for (j = 1; j <= i__1; ++j) {
			dplus = d__[j] + s;
			if (dplus < 0.f) {
			    ++cnt;
			}
			tmp = lld[j] / dplus;
			if (tmp == 0.f) {
			    s = lld[j] - left;
			} else {
			    s = s * tmp - left;
			}
/* L55: */
		    }
		    dplus = d__[*n] + s;
		    if (dplus < 0.f) {
			++cnt;
		    }
		}
		if (cnt > i__ - 1) {
		    tmp = left - werr[ii] * fac;
/*                 A zero, NaN, or overflowing step cannot enlarge the */
/*                 bracket.  Return a positive INFO instead of repeating */
/*                 this Sturm count forever. */
		    if (! (tmp < left && abs(tmp) <= ovfl)) {
			*info = 1;
			return;
		    }
		    ++nexp;
		    if (nexp > 256) {
			*info = 2;
			return;
		    }
		    left = tmp;
		    fac *= 2.f;
		    goto L40;
		}
	    }
	    nleft = cnt + 1;
	    i1 = f2cmin(i1,nleft);
	    fac = 1.f;
	    nexp = 0;
	    right = w[ii] + werr[ii];

/*           Do while( CNT(RIGHT).LT.I ) */

L60:
	    s = -right;
	    cnt = 0;
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
		dplus = d__[j] + s;
		s = s * lld[j] / dplus - right;
		if (dplus < 0.f) {
		    ++cnt;
		}
/* L70: */
	    }
	    dplus = d__[*n] + s;
	    if (dplus < 0.f) {
		++cnt;
	    }
	    if (! (s > 0.f || s < 1.f)) {

/*                 Runs a slower version of the above loop if a NaN is detected */

		cnt = 0;
		s = -right;
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
		    dplus = d__[j] + s;
		    if (dplus < 0.f) {
			++cnt;
		    }
		    tmp = lld[j] / dplus;
		    if (tmp == 0.f) {
			s = lld[j] - right;
		    } else {
			s = s * tmp - right;
		    }
/* L75: */
		}
		dplus = d__[*n] + s;
		if (dplus < 0.f) {
		    ++cnt;
		}
	    }
	    if (cnt < i__) {
		tmp = right + werr[ii] * fac;
/*              Apply the same progress/finite checks to the upper bracket. */
		if (! (tmp > right && abs(tmp) <= ovfl)) {
		    *info = 1;
		    return;
		}
		++nexp;
		if (nexp > 256) {
		    *info = 2;
		    return;
		}
		right = tmp;
		fac *= 2.f;
		goto L60;
	    }
	    cnt = f2cmin(cnt,i2);
	    ++nint;
	    k = nleft << 1;
	    work[k - 1] = left;
	    work[k] = right;
	    i__ = cnt + 1;
	    iwork[k - 1] = i__;
	    iwork[k] = cnt;
	    if (prev != nleft - 1) {
		work[k - 2] = left;
	    }
	    prev = nleft;
	} else {
	    right = work[i__ * 2];

/*           Remove converged interval from linked list */

	    ++iwork[k - 1];
	    prev = i__;
	    ++i__;
	}
	goto L30;
    }
    if (i__ <= *n && iwork[(i__ << 1) - 1] != -1) {
	work[(i__ << 1) - 1] = work[prev * 2];
    }

/*     Do while( NINT.GT.0 ).  In exact arithmetic each pass halves every */
/*     active interval; cap the passes so NaN/non-progress cannot spin. */

    iter = 0;
L80:
    ++iter;
    if (iter > 256) {
	*info = 3;
	return;
    }
    prev = i1 - 1;
    olnint = nint;
    i__ = i1;
    p = 1;
L85:

/*     Snapshot up to NSTURM intervals from the pass's original linked */
/*     list, then interleave their independent Sturm recurrences.  The */
/*     interval updates below remain in their original order.  A split can */
/*     only insert an interval inside the current [I,NRIGHT] range, so it */
/*     cannot change the saved bounds or successor of a later original */
/*     interval in this pass. */

/* Computing MIN */
    i__1 = 4, i__2 = olnint - p + 1;
    nb = f2cmin(i__1,i__2);
    bafter = i__;
    i__1 = nb;
    for (b = 1; b <= i__1; ++b) {
	bi[b - 1] = bafter;
	k = bafter << 1;
	bmid[b - 1] = (work[k - 1] + work[k]) * .5f;
	bnext[b - 1] = iwork[k - 1];
	bnrght[b - 1] = iwork[k];
	bafter = bnext[b - 1];
/* L87: */
    }
/*     The four-lane count SLARRBTGK_STURM4 is kept below but not used; */
/*     each midpoint is counted by the reference one-shift routine. */
    i__1 = nb;
    for (b = 1; b <= i__1; ++b) {
	slarrbtgk_sturm_slow_(n, &d__[1], &lld[1], &bmid[b - 1], &bcnt[b - 1]
		);
/* L88: */
    }

    i__1 = nb;
    for (b = 1; b <= i__1; ++b) {
	i__ = bi[b - 1];
	k = i__ << 1;
	left = work[k - 1];
	right = work[k];
	next = bnext[b - 1];
	nright = bnrght[b - 1];
	mid = bmid[b - 1];
	width = right - mid;
/* Computing MAX */
	r__1 = abs(left), r__2 = abs(right);
	tmp = f2cmax(r__1,r__2);

/*        Check for convergence of interval if there is only one */
/*        eigenvalue in the interval. */

	gap = 0.f;
	if (i__ == nright) {
	    if (prev > 0 && next <= *n) {
/* Computing MIN */
		r__1 = left - work[k - 2], r__2 = work[k + 1] - right;
		gap = f2cmin(r__1,r__2);
	    } else if (prev > 0) {
		gap = left - work[k - 2];
	    } else if (next <= *n) {
		gap = work[k + 1] - right;
	    }
	}
/* Computing MAX */
	r__1 = *rtol1 * gap, r__2 = *rtol2 * tmp;
	if (width < f2cmax(r__1,r__2)) {
	    --nint;
	    iwork[k - 1] = 0;
	    kk = k;
	    i__2 = nright;
	    for (j = i__ + 1; j <= i__2; ++j) {
		kk += 2;
		iwork[kk - 1] = 0;
		work[kk - 1] = left;
		work[kk] = right;
		wgap[j - 1 - *offset] = 0.f;
	    }
	    if (i1 == i__) {
		i1 = next;
	    } else {
		iwork[(prev << 1) - 1] = next;
	    }
	    i__ = next;
	    goto L100;
	}
	prev = i__;

/*        Perform one bisection step.  BCNT(B) is exactly the scalar Sturm */
/*        count for MID; only independent recurrences were interleaved. */

	cnt = bcnt[b - 1];
/* Computing MAX */
	i__2 = i__ - 1, i__3 = f2cmin(nright,cnt);
	cnt = f2cmax(i__2,i__3);
	if (cnt == i__ - 1) {
	    work[k - 1] = mid;
	} else if (cnt == nright) {
	    work[k] = mid;
	} else {
	    iwork[k] = cnt;
	    ++cnt;
	    iwork[k - 1] = cnt;
	    kk = cnt << 1;
	    iwork[kk - 1] = next;
	    iwork[kk] = nright;
	    work[k] = mid;
	    work[kk - 1] = mid;
	    work[kk] = right;
	    prev = cnt;
	    if (cnt - 1 > i__) {
		work[kk - 2] = mid;
	    }
	    if (cnt > *ifirst && cnt <= *ilast) {
		++nint;
	    } else if (cnt <= *ifirst) {
		i1 = cnt;
	    }
	}
	i__ = next;
L100:
	;
    }
    i__ = bafter;
    p += nb;
    if (p <= olnint) {
	goto L85;
    }
    if (nint > 0) {
	goto L80;
    }
    i__1 = *ilast;
    for (i__ = *ifirst; i__ <= i__1; ++i__) {
	k = i__ << 1;
	ii = i__ - *offset;
	if (iwork[k - 1] != -1) {
	    w[ii] = (work[k - 1] + work[k]) * .5f;
	    werr[ii] = work[k] - w[ii];
	    if (i__ != *ilast) {
		wgap[ii] = work[k + 1] - work[k];
	    }
	}
/* L110: */
    }

    return;

/*     End of SLARRB */

} /* slarrb_tgk__ */

/* Subroutine */ void slarrbtgk_sturm4_(integer *n, real *d__, real *lld, 
	integer *nb, real *x, integer *cnt)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer j, c1, c2, c3, c4;
    real s1, s2, s3, s4, dp1, dp2, dp3, dp4;
    extern /* Subroutine */ void slarrbtgk_sturm_slow_(integer *, real *, 
	    real *, real *, integer *);


/*     Evaluate one to four independent Sturm counts while preserving the */
/*     scalar operation order of each recurrence.  Interleaving exposes the */
/*     otherwise serial floating-point divisions to an out-of-order CPU. */


    /* Parameter adjustments */
    --cnt;
    --x;
    --lld;
    --d__;

    /* Function Body */
    c1 = 0;
    s1 = -x[1];
    if (*nb == 1) {
	i__1 = *n - 1;
	for (j = 1; j <= i__1; ++j) {
	    dp1 = d__[j] + s1;
	    s1 = s1 * lld[j] / dp1 - x[1];
	    if (dp1 < 0.f) {
		++c1;
	    }
/* L10: */
	}
	dp1 = d__[*n] + s1;
	if (dp1 < 0.f) {
	    ++c1;
	}
	cnt[1] = c1;
    } else if (*nb == 2) {
	c2 = 0;
	s2 = -x[2];
	i__1 = *n - 1;
	for (j = 1; j <= i__1; ++j) {
	    dp1 = d__[j] + s1;
	    dp2 = d__[j] + s2;
	    s1 = s1 * lld[j] / dp1 - x[1];
	    s2 = s2 * lld[j] / dp2 - x[2];
	    if (dp1 < 0.f) {
		++c1;
	    }
	    if (dp2 < 0.f) {
		++c2;
	    }
/* L20: */
	}
	dp1 = d__[*n] + s1;
	dp2 = d__[*n] + s2;
	if (dp1 < 0.f) {
	    ++c1;
	}
	if (dp2 < 0.f) {
	    ++c2;
	}
	cnt[1] = c1;
	cnt[2] = c2;
    } else if (*nb == 3) {
	c2 = 0;
	c3 = 0;
	s2 = -x[2];
	s3 = -x[3];
	i__1 = *n - 1;
	for (j = 1; j <= i__1; ++j) {
	    dp1 = d__[j] + s1;
	    dp2 = d__[j] + s2;
	    dp3 = d__[j] + s3;
	    s1 = s1 * lld[j] / dp1 - x[1];
	    s2 = s2 * lld[j] / dp2 - x[2];
	    s3 = s3 * lld[j] / dp3 - x[3];
	    if (dp1 < 0.f) {
		++c1;
	    }
	    if (dp2 < 0.f) {
		++c2;
	    }
	    if (dp3 < 0.f) {
		++c3;
	    }
/* L30: */
	}
	dp1 = d__[*n] + s1;
	dp2 = d__[*n] + s2;
	dp3 = d__[*n] + s3;
	if (dp1 < 0.f) {
	    ++c1;
	}
	if (dp2 < 0.f) {
	    ++c2;
	}
	if (dp3 < 0.f) {
	    ++c3;
	}
	cnt[1] = c1;
	cnt[2] = c2;
	cnt[3] = c3;
    } else {
	c2 = 0;
	c3 = 0;
	c4 = 0;
	s2 = -x[2];
	s3 = -x[3];
	s4 = -x[4];
	i__1 = *n - 1;
	for (j = 1; j <= i__1; ++j) {
	    dp1 = d__[j] + s1;
	    dp2 = d__[j] + s2;
	    dp3 = d__[j] + s3;
	    dp4 = d__[j] + s4;
	    s1 = s1 * lld[j] / dp1 - x[1];
	    s2 = s2 * lld[j] / dp2 - x[2];
	    s3 = s3 * lld[j] / dp3 - x[3];
	    s4 = s4 * lld[j] / dp4 - x[4];
	    if (dp1 < 0.f) {
		++c1;
	    }
	    if (dp2 < 0.f) {
		++c2;
	    }
	    if (dp3 < 0.f) {
		++c3;
	    }
	    if (dp4 < 0.f) {
		++c4;
	    }
/* L40: */
	}
	dp1 = d__[*n] + s1;
	dp2 = d__[*n] + s2;
	dp3 = d__[*n] + s3;
	dp4 = d__[*n] + s4;
	if (dp1 < 0.f) {
	    ++c1;
	}
	if (dp2 < 0.f) {
	    ++c2;
	}
	if (dp3 < 0.f) {
	    ++c3;
	}
	if (dp4 < 0.f) {
	    ++c4;
	}
	cnt[1] = c1;
	cnt[2] = c2;
	cnt[3] = c3;
	cnt[4] = c4;
    }

/*     Match the original guarded recurrence exactly if a fast recurrence */
/*     produced NaN.  This path is exceptional; the normal path above is */
/*     the performance target. */

    if (! (s1 > 0.f || s1 < 1.f)) {
	slarrbtgk_sturm_slow_(n, &d__[1], &lld[1], &x[1], &cnt[1]);
    }
    if (*nb >= 2) {
	if (! (s2 > 0.f || s2 < 1.f)) {
	    slarrbtgk_sturm_slow_(n, &d__[1], &lld[1], &x[2], &cnt[2]);
	}
    }
    if (*nb >= 3) {
	if (! (s3 > 0.f || s3 < 1.f)) {
	    slarrbtgk_sturm_slow_(n, &d__[1], &lld[1], &x[3], &cnt[3]);
	}
    }
    if (*nb >= 4) {
	if (! (s4 > 0.f || s4 < 1.f)) {
	    slarrbtgk_sturm_slow_(n, &d__[1], &lld[1], &x[4], &cnt[4]);
	}
    }
    return;
} /* slarrbtgk_sturm4__ */

/* Subroutine */ void slarrbtgk_sturm_slow_(integer *n, real *d__, real *lld, 
	real *x, integer *cnt)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer j;
    real s, tmp, dplus;


/*     NaN-safe scalar recurrence, identical to SLARRB's original slow */
/*     path.  It is split out only so each lane of SLARRBTGK_STURM4 can use it. */


    /* Parameter adjustments */
    --lld;
    --d__;

    /* Function Body */
    *cnt = 0;
    s = -(*x);
    i__1 = *n - 1;
    for (j = 1; j <= i__1; ++j) {
	dplus = d__[j] + s;
	if (dplus < 0.f) {
	    ++(*cnt);
	}
	tmp = lld[j] / dplus;
	if (tmp == 0.f) {
	    s = lld[j] - *x;
	} else {
	    s = s * tmp - *x;
	}
/* L10: */
    }
    dplus = d__[*n] + s;
    if (dplus < 0.f) {
	++(*cnt);
    }
    return;
} /* slarrbtgk_sturm_slow__ */

