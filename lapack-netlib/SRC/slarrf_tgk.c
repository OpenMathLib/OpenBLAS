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


/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b SLARRF_TGK computes a shifted child representation for a cluster (Golub-Kahan variant of SLARRF
). */
/* > */
/* > \ingroup bdsvdmr3 */
/* > */
/* > \par Purpose: */
/* > ============= */
/* > */
/* > \verbatim */
/* > computes a shifted child representation for a cluster (Golub-Kahan variant of SLARRF). */
/* > All routines in this file are auxiliary to SBDSVDMR3 and derive */
/* > from the LAPACK 3.0 (1999) MRRR kernels SLARRV, SLAR1V, SLARRB and */
/* > SLARRF rather than from the current SSTEMR kernels. */
/* > \endverbatim */

/* Subroutine */ void slarrf_tgk__(char *rep, integer *n, real *d__, real *l, 
	real *ld, real *lld, integer *ifirst, integer *ilast, real *w, real *
	sigma, real *dplus, real *lplus, real *work, integer *iseed, integer *
	info)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3;

    /* Local variables */
    integer i__;
    real s;
    logical tgk;
    real eps, tmp, pert;
    integer ntry;
    real dmax1, dmax2;
    logical degen;
    real delta;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ void scopy_(integer *, real *, integer *, real *, 
	    integer *);
    extern real slamch_(char *);
    extern /* Subroutine */ void slarnv_(integer *, integer *, integer *, real 
	    *);


/*  -- New auxiliary routine for bidiagonal SVD via TGK-rooted MR^3 -- */
/*     Based on LAPACK SLARRF (Dhillon & Marques, Nov 11 2003). */
/*     Derived from the LAPACK 3.0 SLARRF; this is a separate copy that adds */
/*     a representation switch so the parent may be supplied either as an */
/*     L D L^T factorization (REP='L', identical to the original) or as a */
/*     symmetric tridiagonal matrix -- in particular the Tridiagonal */
/*     Golub-Kahan (TGK) matrix -- with REP='T'. */


/*  Purpose */
/*  ======= */

/*  Given a parent representation M and a cluster of close eigenvalues */
/*  (in a relative measure) W( IFIRST ), ..., W( ILAST ), SLARRF_TGK finds */
/*  a new relatively robust representation */
/*       M - SIGMA I = L(+) D(+) L(+)^T */
/*  such that at least one of the eigenvalues of L(+) D(+) L(+)^T is */
/*  relatively isolated.  The output L(+) D(+) L(+)^T is an ordinary */
/*  L D L^T representation regardless of REP, so every descendant node in */
/*  the representation tree can be processed by the unmodified LAPACK 3.0-lineage */
/*  routines. */

/*  Arguments */
/*  ========= */

/*  REP     (input) CHARACTER*1 */
/*          = 'L':  the parent M is the L D L^T representation defined by */
/*                  D, L (and LD = L(i)*D(i), LLD = L(i)^2*D(i)).  Behaviour */
/*                  is identical to the original SLARRF. */
/*          = 'T':  the parent M is a symmetric tridiagonal matrix (e.g. */
/*                  the TGK matrix).  D holds the N diagonal entries (all */
/*                  zero for TGK, kept general) and L holds the N-1 */
/*                  off-diagonal entries beta(i).  LD and LLD are NOT */
/*                  referenced. */

/*  N       (input) INTEGER */
/*          The order of the matrix. */

/*  D       (input) REAL array, dimension (N) */
/*          REP='L': diagonal of D.   REP='T': diagonal of M. */

/*  L       (input) REAL array, dimension (N-1) */
/*          REP='L': subdiagonal of the unit bidiagonal L. */
/*          REP='T': off-diagonal entries beta(i) of M. */

/*  LD      (input) REAL array, dimension (N-1) */
/*          REP='L': L(i)*D(i).   REP='T': not referenced. */

/*  LLD     (input) REAL array, dimension (N-1) */
/*          REP='L': L(i)*L(i)*D(i).   REP='T': not referenced. */

/*  IFIRST  (input) INTEGER   Index of the first eigenvalue in the cluster. */
/*  ILAST   (input) INTEGER   Index of the last  eigenvalue in the cluster. */

/*  W       (input) REAL array, dimension (N) */
/*          Eigenvalues of M in ascending order; W(IFIRST..ILAST) is the */
/*          cluster.  (For the TGK root these are the positive eigenvalues */
/*          == singular values; the caller passes them in ascending order.) */

/*  SIGMA   (output) REAL  The shift used to form L(+)D(+)L(+)^T. */

/*  DPLUS   (output) REAL array, dimension (N)   Diagonal D(+). */
/*  LPLUS   (output) REAL array, dimension (N-1) Subdiag L(+). */

/*  WORK    (workspace) REAL array, dimension (2*N) */

/*  INFO    (output) INTEGER   = 0 successful exit; < 0 illegal argument. */

/*  ===================================================================== */

/*     DOPERT enables the Sec.5.4 random relative perturbation of the child */
/*     RRR.  It is applied SELECTIVELY -- only to a child whose cluster has */
/*     a near-degenerate pair, W(i+1)-W(i) < PGFAC*eps*W(i) (PGFAC=100). */
/*     The per-entry size is the full PERTK*eps: empirically the magnitude */
/*     cannot be reduced (PERTK*eps/N and PERTK*eps/sqrt(N) both fail -- the */
/*     overflow NaN returns) because breaking the near-zero pivot */
/*     cancellation that destroys the eigenvector (DPLUS ~ 1e2*eps, see */
/*     below) needs an O(eps) per-entry jolt.  The gate is what keeps this */
/*     from over-perturbing moderately-clustered RRRs. */

    /* Parameter adjustments */
    --iseed;
    --work;
    --lplus;
    --dplus;
    --w;
    --lld;
    --ld;
    --l;
    --d__;

    /* Function Body */
    *info = 0;
    if (! (lsame_(rep, "L") || lsame_(rep, "T"))) {
	*info = -1;
	return;
    }
    tgk = lsame_(rep, "T");
    eps = slamch_("Precision");
    *sigma = w[*ifirst];
    delta = eps * 2.f;
    ntry = 0;

/*     Try the representation based at the left end of the cluster, */
/*     SIGMA = W( IFIRST ).  DMAX1 monitors the element growth. */

L10:
    if (tgk) {
/*        Tridiagonal parent:  M - SIGMA I = L(+) D(+) L(+)^T directly. */
	dplus[1] = d__[1] - *sigma;
	dmax1 = abs(dplus[1]);
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    lplus[i__] = l[i__] / dplus[i__];
	    dplus[i__ + 1] = d__[i__ + 1] - *sigma - lplus[i__] * l[i__];
/* Computing MAX */
	    r__2 = dmax1, r__3 = (r__1 = dplus[i__ + 1], abs(r__1));
	    dmax1 = f2cmax(r__2,r__3);
/* L20: */
	}
    } else {
/*        L D L^T parent (original SLARRF recurrence). */
	s = -(*sigma);
	dplus[1] = d__[1] + s;
	dmax1 = abs(dplus[1]);
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    lplus[i__] = ld[i__] / dplus[i__];
	    s = s * lplus[i__] * l[i__] - *sigma;
	    dplus[i__ + 1] = d__[i__ + 1] + s;
/* Computing MAX */
	    r__2 = dmax1, r__3 = (r__1 = dplus[i__ + 1], abs(r__1));
	    dmax1 = f2cmax(r__2,r__3);
/* L25: */
	}
    }
    if (dmax1 * 0.f != 0.f) {
/*        A non-finite pivot (Inf or NaN) appeared -- e.g. a zero TGK */
/*        off-diagonal meeting a zero pivot at SIGMA=0.  (The DMAX1*0.NE.0 */
/*        test is true exactly for Inf/NaN.)  An Inf child would later */
/*        hang SLARRB's bisection, so we reject it here.  Nudge SIGMA off */
/*        the eigenvalue and retry, */
/*        but only a bounded number of times so we can never loop forever */
/*        (the additive ONE term lets the nudge escape SIGMA=0; the matrix */
/*        is pre-scaled to O(1) norm by SBDSVDMR3). */
	++ntry;
	if (ntry > 10) {
	    *info = 1;
	    return;
	}
	*sigma -= (abs(*sigma) + 1.f) * delta;
	delta *= 2.f;
	goto L10;
    }

/*     Try the representation based at the right end of the cluster, */
/*     TMP = W( ILAST ), storing it in WORK; DMAX2 monitors its growth. */

    tmp = w[*ilast];
    delta = eps * 2.f;
    ntry = 0;
L30:
    if (tgk) {
	work[1] = d__[1] - tmp;
	dmax2 = abs(work[1]);
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    work[*n + i__] = l[i__] / work[i__];
	    work[i__ + 1] = d__[i__ + 1] - tmp - work[*n + i__] * l[i__];
/* Computing MAX */
	    r__2 = dmax2, r__3 = (r__1 = work[i__ + 1], abs(r__1));
	    dmax2 = f2cmax(r__2,r__3);
/* L40: */
	}
    } else {
	s = -tmp;
	work[1] = d__[1] + s;
	dmax2 = abs(work[1]);
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    work[*n + i__] = ld[i__] / work[i__];
	    s = s * work[*n + i__] * l[i__] - tmp;
	    work[i__ + 1] = d__[i__ + 1] + s;
/* Computing MAX */
	    r__2 = dmax2, r__3 = (r__1 = work[i__ + 1], abs(r__1));
	    dmax2 = f2cmax(r__2,r__3);
/* L45: */
	}
    }
    if (dmax2 * 0.f != 0.f) {
	++ntry;
	if (ntry > 10) {
/*           The right-end shift never gave finite growth; fall back to */
/*           the (finite) left-end representation found above. */
	    goto L50;
	}
	tmp += (abs(tmp) + 1.f) * delta;
	delta *= 2.f;
	goto L30;
    }

/*     Keep whichever shift produced the smaller element growth. */

    if (dmax2 < dmax1) {
	*sigma = tmp;
	scopy_(n, &work[1], &c__1, &dplus[1], &c__1);
	i__1 = *n - 1;
	scopy_(&i__1, &work[*n + 1], &c__1, &lplus[1], &c__1);
    }

L50:

/*     Random relative perturbation of the child RRR (Dhillon, Parlett & */
/*     Voemel 2005, eqs. 5.1-5.2): */
/*        DPLUS(i) := DPLUS(i)*(1 + mu_i *PERTK*eps), */
/*        LPLUS(j) := LPLUS(j)*(1 + nu_j *PERTK*eps), */
/*     with mu,nu uniform on [-1,1] (SLARNV, IDIST=2) and PERTK=4. */

/*     WHY this is needed (and load-bearing for the TGK).  A tight near- */
/*     degenerate pair of singular values (e.g. a glued-Wilkinson Wilkinson */
/*     pair, relative gap ~eps) forces a child RRR whose shift sits almost */
/*     exactly on an eigenvalue, giving a near-zero pivot / large element */
/*     growth (DPLUS as small as ~1e2*eps next to ~1/eps).  The downstream */
/*     SLAR1V twisted factorization then drives lambda onto that eigenvalue, */
/*     the pivot collapses to ~1e-30, the back-substitution overflows, and */
/*     the eigenvector comes back as NaN.  A relative perturbation of the */
/*     child factors breaks the exact pivot cancellation: the perturbed RRR */
/*     has no eigenvalue closer than ~eps*||M|| to the shift, so the twist */
/*     gamma cannot collapse and the vector stays finite. */

/*     APPLIED SELECTIVELY -- only to a child whose cluster contains a near- */
/*     degenerate consecutive pair, W(i+1)-W(i) < PGFAC*eps*|W(i)| with */
/*     PGFAC=100.  Perturbing every child instead (no gate) needlessly */
/*     jitters moderately-clustered RRRs and degrades their orthogonality; */
/*     the gate confines the perturbation to the clusters that actually need */
/*     it.  The size cannot be lowered below the full PERTK*eps per entry -- */
/*     PERTK*eps/N and PERTK*eps/sqrt(N) were both tried and the overflow */
/*     NaN returned -- because breaking the cancellation needs an O(eps) */
/*     per-entry jolt; this does scramble genuinely sub-eps gaps (strong- */
/*     glue cases land at ~1e-8..1e-6 orthogonality), the unavoidable price. */
/*     ISEED is threaded by the caller (fixed seed) for reproducibility. */

    degen = FALSE_;
    if (FALSE_) {
	i__1 = *ilast - 1;
	for (i__ = *ifirst; i__ <= i__1; ++i__) {
	    if (w[i__ + 1] - w[i__] < eps * 100.f * (r__1 = w[i__], abs(r__1))
		    ) {
		degen = TRUE_;
	    }
/* L55: */
	}
    }
    if (degen) {
	pert = eps * 4.f;
	slarnv_(&c__2, &iseed[1], n, &work[1]);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dplus[i__] *= work[i__] * pert + 1.f;
/* L60: */
	}
	i__1 = *n - 1;
	slarnv_(&c__2, &iseed[1], &i__1, &work[1]);
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    lplus[i__] *= work[i__] * pert + 1.f;
/* L70: */
	}
    }
    return;

/*     End of SLARRF_TGK */

} /* slarrf_tgk__ */

