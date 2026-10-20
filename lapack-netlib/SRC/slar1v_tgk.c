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


/* > \brief \b SLAR1V_TGK computes one eigenvector of a shifted tridiagonal from a twisted factorization, keep
ing the full support (Golub-Kahan variant of SLAR1V); also SLAR1V2_TGK, the two-eigenvalue variant. */
/* > */
/* > \ingroup bdsvdmr3 */
/* > */
/* > \par Purpose: */
/* > ============= */
/* > */
/* > \verbatim */
/* > computes one eigenvector of a shifted tridiagonal from a twisted factorization, keeping the full support 
(Golub-Kahan variant of SLAR1V); also SLAR1V2_TGK, the two-eigenvalue variant. */
/* > All routines in this file are auxiliary to SBDSVDMR3 and derive */
/* > from the LAPACK 3.0 (1999) MRRR kernels SLARRV, SLAR1V, SLARRB and */
/* > SLARRF rather than from the current SSTEMR kernels. */
/* > \endverbatim */

/* Subroutine */ void slar1v_tgk__(char *rep, integer *n, integer *b1, integer 
	*bn, real *sigma, real *d__, real *l, real *ld, real *lld, real *eval,
	 real *gersch, real *z__, real *ztz, real *mingma, integer *r__, 
	integer *isuppz, real *work)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, j;
    real s;
    integer r1, r2;
    logical tgk;
    real eps, tmp;
    integer indp, inds;
    extern logical lsame_(char *, char *);
    real dplus;
    extern real slamch_(char *);
    logical sawnan;
    integer indumn;
    real dminus;


/*  -- New auxiliary routine for bidiagonal SVD via TGK-rooted MR^3 -- */
/*     Based on LAPACK SLAR1V (Dhillon & Marques, Nov 11 2003). */
/*     Derived from the LAPACK 3.0 SLAR1V; this is a separate copy that adds a */
/*     representation switch so the matrix M may be supplied either as an */
/*     L D L^T factorization (REP='L', identical to the original) or as a */
/*     symmetric tridiagonal -- in particular the Tridiagonal Golub-Kahan */
/*     (TGK) matrix -- with REP='T'. */


/*  Purpose */
/*  ======= */

/*  SLAR1V_TGK computes the (scaled) r-th column of the inverse of the */
/*  submatrix in rows B1 through BN of  M - sigma I,  where M is either */
/*  an L D L^T representation (REP='L') or a symmetric tridiagonal matrix */
/*  such as the TGK matrix (REP='T').  When sigma is close to an */
/*  eigenvalue, the computed vector is an accurate eigenvector.  The steps */
/*  are those of SLAR1V: */
/*  (a) stationary qd transform   M - sigma I = L(+) D(+) L(+)^T, */
/*  (b) progressive qd transform  M - sigma I = U(-) D(-) U(-)^T, */
/*  (c) the twist index r = argmin_k |gamma_k|, where gamma_k is the k-th */
/*      diagonal element of the twisted factorization and equals the */
/*      reciprocal of (M - sigma I)^{-1}_{kk}, */
/*  (d) the scaled r-th column of the inverse from the twisted factor. */

/*  For REP='T' the qd transforms run on the tridiagonal entries directly */
/*  (diagonal in D, off-diagonals beta_i in L), so LD and LLD are not */
/*  referenced.  Using the identity  gamma_k = D(+)(k) + D(-)(k) - J_kk */
/*  with J_kk = D(k) - sigma (Dhillon, PhD thesis, UC Berkeley 1997, eq. 3.1.12), the WORK slots are filled */
/*  so that the twist selection and back-substitution code below are */
/*  shared, unchanged, with the REP='L' path: */
/*       WORK(INDS + k-1) = D(+)(k) */
/*       WORK(INDP + k-1) = D(-)(k) - ( D(k) - sigma ) */
/*       WORK(k)          = L(+)(k) = beta_k / D(+)(k) */
/*       WORK(INDUMN + k) = U(-)(k) = beta_k / D(-)(k+1) */
/*  so WORK(INDS+k-1) + WORK(INDP+k-1) = gamma_k.  For the TGK matrix */
/*  (D == 0) this is gamma_k = D(+)(k) + D(-)(k) + sigma. */

/*  Zero pivots are expected for the indefinite zero-diagonal TGK root. */
/*  In the direct recurrence a zero pivot self-heals through IEEE */
/*  arithmetic (D(+) -> -Inf, then L(+) -> 0), so no stationary/progressive */
/*  restart is needed; only the back-substitution uses a slow path (to */
/*  avoid Inf*0 = NaN), with beta-ratios per Algorithm 3.3.1 of Dhillon's thesis (UC Berkeley, 1997). */

/*  Arguments */
/*  ========= */

/*  REP   (input) CHARACTER*1   'L' = L D L^T input; 'T' = tridiagonal/TGK. */
/*  N,B1,BN,SIGMA,EVAL,GERSCH,Z,ZTZ,MINGMA,R,ISUPPZ,WORK : as in SLAR1V. */
/*  D     (input) REAL (N)   REP='L': diag(D).  REP='T': diag(M). */
/*  L     (input) REAL (N-1) REP='L': subdiag(L). REP='T': beta_i. */
/*  LD,LLD(input) REAL (N-1) REP='L': L*D, L*L*D. REP='T': unused. */

/*  ===================================================================== */


    /* Parameter adjustments */
    --work;
    --isuppz;
    --z__;
    --gersch;
    --lld;
    --ld;
    --l;
    --d__;

    /* Function Body */
    tgk = lsame_(rep, "T");
    eps = slamch_("Precision");
    if (*r__ == 0) {

/*        Eliminate the top and bottom indices from among the possible */
/*        values of R, using the Gerschgorin intervals. */

	r1 = *b1;
	r2 = *bn;
	i__1 = *bn;
	for (i__ = *b1; i__ <= i__1; ++i__) {
	    if (*eval >= gersch[(i__ << 1) - 1] && *eval <= gersch[i__ * 2]) {
		r1 = i__;
		goto L20;
	    }
/* L10: */
	}
	goto L40;
L20:
	i__1 = *b1;
	for (i__ = *bn; i__ >= i__1; --i__) {
	    if (*eval >= gersch[(i__ << 1) - 1] && *eval <= gersch[i__ * 2]) {
		r2 = i__;
		goto L40;
	    }
/* L30: */
	}
    } else {
	r1 = *r__;
	r2 = *r__;
    }

L40:
    indumn = *n;
    inds = (*n << 1) + 1;
    indp = *n * 3 + 1;
    sawnan = FALSE_;

/*     ================================================================ */
/*     Stationary transform (differential form) up to the index R2. */
/*     ================================================================ */

    if (tgk) {

/*        Tridiagonal/TGK parent.  Direct recurrence; store D(+)(k). */

	dplus = d__[*b1] - *sigma;
	work[inds + *b1 - 1] = dplus;
	if (dplus == 0.f) {
	    sawnan = TRUE_;
	}
	i__1 = r2 - 1;
	for (i__ = *b1; i__ <= i__1; ++i__) {
	    work[i__] = l[i__] / dplus;
	    dplus = d__[i__ + 1] - *sigma - work[i__] * l[i__];
	    work[inds + i__] = dplus;
	    if (dplus == 0.f) {
		sawnan = TRUE_;
	    }
/* L45: */
	}
    } else {

/*        L D L^T parent (original SLAR1V stationary loop). */

	if (*b1 == 1) {
	    work[inds] = 0.f;
	} else {
	    work[inds] = lld[*b1 - 1];
	}
	s = work[inds] - *sigma;
	i__1 = r2 - 1;
	for (i__ = *b1; i__ <= i__1; ++i__) {
	    dplus = d__[i__] + s;
	    work[i__] = ld[i__] / dplus;
	    work[inds + i__] = s * work[i__] * l[i__];
	    s = work[inds + i__] - *sigma;
/* L50: */
	}

	if (! (s > 0.f || s < 1.f)) {
/*           Slower version of the above loop if a NaN is detected. */
	    sawnan = TRUE_;
	    j = *b1 + 1;
L60:
	    if (work[inds + j] > 0.f || work[inds + j] < 1.f) {
		++j;
		goto L60;
	    }
	    work[inds + j] = lld[j];
	    s = work[inds + j] - *sigma;
	    i__1 = r2 - 1;
	    for (i__ = j + 1; i__ <= i__1; ++i__) {
		dplus = d__[i__] + s;
		work[i__] = ld[i__] / dplus;
		if (work[i__] == 0.f) {
		    work[inds + i__] = lld[i__];
		} else {
		    work[inds + i__] = s * work[i__] * l[i__];
		}
		s = work[inds + i__] - *sigma;
/* L70: */
	    }
	}
    }

/*     ================================================================ */
/*     Progressive transform (differential form) down to the index R1. */
/*     ================================================================ */

    if (tgk) {

/*        Tridiagonal/TGK parent.  Direct recurrence.  Store */
/*        WORK(INDP+i-1) = D(-)(i) - (D(i)-sigma)  and  U(-)(i). */

	dminus = d__[*bn] - *sigma;
	work[indp + *bn - 1] = 0.f;
	if (dminus == 0.f) {
	    sawnan = TRUE_;
	}
	i__1 = r1;
	for (i__ = *bn - 1; i__ >= i__1; --i__) {
	    work[indumn + i__] = l[i__] / dminus;
	    dminus = d__[i__] - *sigma - work[indumn + i__] * l[i__];
	    work[indp + i__ - 1] = dminus - (d__[i__] - *sigma);
	    if (dminus == 0.f) {
		sawnan = TRUE_;
	    }
/* L75: */
	}
    } else {

/*        L D L^T parent (original SLAR1V progressive loop). */

	work[indp + *bn - 1] = d__[*bn] - *sigma;
	i__1 = r1;
	for (i__ = *bn - 1; i__ >= i__1; --i__) {
	    dminus = lld[i__] + work[indp + i__];
	    tmp = d__[i__] / dminus;
	    work[indumn + i__] = l[i__] * tmp;
	    work[indp + i__ - 1] = work[indp + i__] * tmp - *sigma;
/* L80: */
	}
	tmp = work[indp + r1 - 1];
	if (! (tmp > 0.f || tmp < 1.f)) {
/*           Slower version of the above loop if a NaN is detected. */
	    sawnan = TRUE_;
	    j = *bn - 3;
L90:
	    if (work[indp + j] > 0.f || work[indp + j] < 1.f) {
		--j;
		goto L90;
	    }
	    work[indp + j] = d__[j + 1] - *sigma;
	    i__1 = r1;
	    for (i__ = j; i__ >= i__1; --i__) {
		dminus = lld[i__] + work[indp + i__];
		tmp = d__[i__] / dminus;
		work[indumn + i__] = l[i__] * tmp;
		if (tmp == 0.f) {
		    work[indp + i__ - 1] = d__[i__] - *sigma;
		} else {
		    work[indp + i__ - 1] = work[indp + i__] * tmp - *sigma;
		}
/* L100: */
	    }
	}
    }

/*     ================================================================ */
/*     Twist index: r = argmin_k |gamma_k|, gamma_k stored as */
/*     WORK(INDS+k-1) + WORK(INDP+k-1) for both representations. */
/*     ================================================================ */

    *mingma = work[inds + r1 - 1] + work[indp + r1 - 1];
    if (*mingma == 0.f) {
	*mingma = eps * work[inds + r1 - 1];
    }
    *r__ = r1;
    i__1 = r2 - 1;
    for (i__ = r1; i__ <= i__1; ++i__) {
	tmp = work[inds + i__] + work[indp + i__];
	if (tmp == 0.f) {
	    tmp = eps * work[inds + i__];
	}
	if (abs(tmp) < abs(*mingma)) {
	    *mingma = tmp;
	    *r__ = i__ + 1;
	}
/* L110: */
    }

/*     ================================================================ */
/*     Compute the (scaled) r-th column of the inverse.  WORK(k)=L(+)(k) */
/*     and WORK(INDUMN+k)=U(-)(k) for both representations, so the fast */
/*     path is shared.  The slow (SAWNAN) path uses LD-ratios for REP='L' */
/*     and beta-ratios (L) for REP='T'. */
/*     ================================================================ */

    isuppz[1] = *b1;
    isuppz[2] = *bn;
    z__[*r__] = 1.f;
    *ztz = 1.f;
    if (! sawnan) {
/*        No support truncation.  A TGK eigenvector can be BIMODAL -- the */
/*        left singular vector u occupies one end of z and the right vector */
/*        v the other, with a near-zero trough between -- so the usual */
/*        "two consecutive tiny entries => the eigenvector has died out, */
/*        stop" early-out would wrongly amputate the far hump.  This is */
/*        true for EVERY vector this routine computes, whether the */
/*        representation is the root TGK (REP='T') or a deeper L D L^T */
/*        child of it (REP='L'), since the eigenvectors are TGK */
/*        eigenvectors in both cases.  Hence we always form the full */
/*        vector and leave ISUPPZ = [B1,BN]. */
	i__1 = *b1;
	for (i__ = *r__ - 1; i__ >= i__1; --i__) {
	    z__[i__] = -(work[i__] * z__[i__ + 1]);
	    *ztz += z__[i__] * z__[i__];
/* L130: */
	}
	i__1 = *bn;
	for (i__ = *r__ + 1; i__ <= i__1; ++i__) {
	    z__[i__] = -(work[indumn + i__ - 1] * z__[i__ - 1]);
	    *ztz += z__[i__] * z__[i__];
/* L150: */
	}
    } else {

/*        Slow path: some pivot vanished, so L(+)/U(-) may be +-Inf. */
/*        Use the recurrence of Algorithm 3.3.1 in Dhillon's thesis, switching the */
/*        coupling ratio between LD (REP='L') and beta=L (REP='T'). */

	i__1 = *b1;
	for (i__ = *r__ - 1; i__ >= i__1; --i__) {
	    if (z__[i__ + 1] == 0.f) {
		if (tgk) {
		    z__[i__] = -(l[i__ + 1] / l[i__]) * z__[i__ + 2];
		} else {
		    z__[i__] = -(ld[i__ + 1] / ld[i__]) * z__[i__ + 2];
		}
	    } else {
		z__[i__] = -(work[i__] * z__[i__ + 1]);
	    }
	    *ztz += z__[i__] * z__[i__];
/* L160: */
	}
/* L170: */
	i__1 = *bn - 1;
	for (i__ = *r__; i__ <= i__1; ++i__) {
	    if (z__[i__] == 0.f) {
		if (tgk) {
		    z__[i__ + 1] = -(l[i__ - 1] / l[i__]) * z__[i__ - 1];
		} else {
		    z__[i__ + 1] = -(ld[i__ - 1] / ld[i__]) * z__[i__ - 1];
		}
	    } else {
		z__[i__ + 1] = -(work[indumn + i__] * z__[i__]);
	    }
	    *ztz += z__[i__ + 1] * z__[i__ + 1];
/* L180: */
	}
/* L190: */
    }

    return;

/*     End of SLAR1V_TGK */

} /* slar1v_tgk__ */

/* Subroutine */ void slar1v2_tgk__(char *rep, integer *n, integer *b1, 
	integer *bn, real *sigma1, real *sigma2, real *d__, real *l, real *ld,
	 real *lld, real *eval1, real *eval2, real *gersch, real *z1, real *
	z2, real *ztz1, real *ztz2, real *mingma1, real *mingma2, integer *
	rone, integer *rtwo, integer *isup1, integer *isup2, real *work1, 
	real *work2)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;
    real sa, sb;
    extern /* Subroutine */ void slar1v_tgk__(char *, integer *, integer *, 
	    integer *, real *, real *, real *, real *, real *, real *, real *,
	     real *, real *, real *, integer *, integer *, real *);
    integer r1a, r1b, r2a, r2b;
    real dma, dmb, dpa, dpb;
    logical tgk;
    real eps;
    logical bad1, bad2;
    integer indp, inds;
    real tmpa, tmpb;
    extern logical lsame_(char *, char *);
    extern real slamch_(char *);
    integer indumn;


/*     Two independent SLAR1V_TGK fast paths interleaved to expose the */
/*     serial qd divisions to the CPU.  Every operation within either lane */
/*     has the same order as the scalar routine.  Exceptional zero/NaN */
/*     pivots fall back to two scalar calls, preserving their exact path. */


    /* Parameter adjustments */
    --work2;
    --work1;
    --isup2;
    --isup1;
    --z2;
    --z1;
    --gersch;
    --lld;
    --ld;
    --l;
    --d__;

    /* Function Body */
    tgk = lsame_(rep, "T");
    eps = slamch_("Precision");
    if (*rone == 0) {
	r1a = *b1;
	r2a = *bn;
	i__1 = *bn;
	for (i__ = *b1; i__ <= i__1; ++i__) {
	    if (*eval1 >= gersch[(i__ << 1) - 1] && *eval1 <= gersch[i__ * 2])
		     {
		r1a = i__;
		goto L20;
	    }
/* L10: */
	}
	goto L40;
L20:
	i__1 = *b1;
	for (i__ = *bn; i__ >= i__1; --i__) {
	    if (*eval1 >= gersch[(i__ << 1) - 1] && *eval1 <= gersch[i__ * 2])
		     {
		r2a = i__;
		goto L40;
	    }
/* L30: */
	}
L40:
	;
    } else {
	r1a = *rone;
	r2a = *rone;
    }
    if (*rtwo == 0) {
	r1b = *b1;
	r2b = *bn;
	i__1 = *bn;
	for (i__ = *b1; i__ <= i__1; ++i__) {
	    if (*eval2 >= gersch[(i__ << 1) - 1] && *eval2 <= gersch[i__ * 2])
		     {
		r1b = i__;
		goto L60;
	    }
/* L50: */
	}
	goto L80;
L60:
	i__1 = *b1;
	for (i__ = *bn; i__ >= i__1; --i__) {
	    if (*eval2 >= gersch[(i__ << 1) - 1] && *eval2 <= gersch[i__ * 2])
		     {
		r2b = i__;
		goto L80;
	    }
/* L70: */
	}
L80:
	;
    } else {
	r1b = *rtwo;
	r2b = *rtwo;
    }

    indumn = *n;
    inds = (*n << 1) + 1;
    indp = *n * 3 + 1;
    bad1 = FALSE_;
    bad2 = FALSE_;
    if (tgk) {
	dpa = d__[*b1] - *sigma1;
	work1[inds + *b1 - 1] = dpa;
	if (dpa == 0.f) {
	    bad1 = TRUE_;
	}
	dpb = d__[*b1] - *sigma2;
	work2[inds + *b1 - 1] = dpb;
	if (dpb == 0.f) {
	    bad2 = TRUE_;
	}
	i__1 = f2cmax(r2a,r2b) - 1;
	for (i__ = *b1; i__ <= i__1; ++i__) {
	    if (i__ < r2a) {
		work1[i__] = l[i__] / dpa;
		dpa = d__[i__ + 1] - *sigma1 - work1[i__] * l[i__];
		work1[inds + i__] = dpa;
		if (dpa == 0.f) {
		    bad1 = TRUE_;
		}
	    }
	    if (i__ < r2b) {
		work2[i__] = l[i__] / dpb;
		dpb = d__[i__ + 1] - *sigma2 - work2[i__] * l[i__];
		work2[inds + i__] = dpb;
		if (dpb == 0.f) {
		    bad2 = TRUE_;
		}
	    }
/* L90: */
	}
    } else {
	if (*b1 == 1) {
	    work1[inds] = 0.f;
	    work2[inds] = 0.f;
	} else {
	    work1[inds] = lld[*b1 - 1];
	    work2[inds] = lld[*b1 - 1];
	}
	sa = work1[inds] - *sigma1;
	sb = work2[inds] - *sigma2;
	i__1 = f2cmax(r2a,r2b) - 1;
	for (i__ = *b1; i__ <= i__1; ++i__) {
	    if (i__ < r2a) {
		dpa = d__[i__] + sa;
		work1[i__] = ld[i__] / dpa;
		work1[inds + i__] = sa * work1[i__] * l[i__];
		sa = work1[inds + i__] - *sigma1;
	    }
	    if (i__ < r2b) {
		dpb = d__[i__] + sb;
		work2[i__] = ld[i__] / dpb;
		work2[inds + i__] = sb * work2[i__] * l[i__];
		sb = work2[inds + i__] - *sigma2;
	    }
/* L100: */
	}
	if (! (sa > 0.f || sa < 1.f)) {
	    bad1 = TRUE_;
	}
	if (! (sb > 0.f || sb < 1.f)) {
	    bad2 = TRUE_;
	}
    }

    if (tgk) {
	dma = d__[*bn] - *sigma1;
	work1[indp + *bn - 1] = 0.f;
	if (dma == 0.f) {
	    bad1 = TRUE_;
	}
	dmb = d__[*bn] - *sigma2;
	work2[indp + *bn - 1] = 0.f;
	if (dmb == 0.f) {
	    bad2 = TRUE_;
	}
	i__1 = f2cmin(r1a,r1b);
	for (i__ = *bn - 1; i__ >= i__1; --i__) {
	    if (i__ >= r1a) {
		work1[indumn + i__] = l[i__] / dma;
		dma = d__[i__] - *sigma1 - work1[indumn + i__] * l[i__];
		work1[indp + i__ - 1] = dma - (d__[i__] - *sigma1);
		if (dma == 0.f) {
		    bad1 = TRUE_;
		}
	    }
	    if (i__ >= r1b) {
		work2[indumn + i__] = l[i__] / dmb;
		dmb = d__[i__] - *sigma2 - work2[indumn + i__] * l[i__];
		work2[indp + i__ - 1] = dmb - (d__[i__] - *sigma2);
		if (dmb == 0.f) {
		    bad2 = TRUE_;
		}
	    }
/* L110: */
	}
    } else {
	work1[indp + *bn - 1] = d__[*bn] - *sigma1;
	work2[indp + *bn - 1] = d__[*bn] - *sigma2;
	i__1 = f2cmin(r1a,r1b);
	for (i__ = *bn - 1; i__ >= i__1; --i__) {
	    if (i__ >= r1a) {
		dma = lld[i__] + work1[indp + i__];
		tmpa = d__[i__] / dma;
		work1[indumn + i__] = l[i__] * tmpa;
		work1[indp + i__ - 1] = work1[indp + i__] * tmpa - *sigma1;
	    }
	    if (i__ >= r1b) {
		dmb = lld[i__] + work2[indp + i__];
		tmpb = d__[i__] / dmb;
		work2[indumn + i__] = l[i__] * tmpb;
		work2[indp + i__ - 1] = work2[indp + i__] * tmpb - *sigma2;
	    }
/* L120: */
	}
	tmpa = work1[indp + r1a - 1];
	tmpb = work2[indp + r1b - 1];
	if (! (tmpa > 0.f || tmpa < 1.f)) {
	    bad1 = TRUE_;
	}
	if (! (tmpb > 0.f || tmpb < 1.f)) {
	    bad2 = TRUE_;
	}
    }
    if (bad1 || bad2) {
	goto L300;
    }

    *mingma1 = work1[inds + r1a - 1] + work1[indp + r1a - 1];
    if (*mingma1 == 0.f) {
	*mingma1 = eps * work1[inds + r1a - 1];
    }
    *rone = r1a;
    i__1 = r2a - 1;
    for (i__ = r1a; i__ <= i__1; ++i__) {
	tmpa = work1[inds + i__] + work1[indp + i__];
	if (tmpa == 0.f) {
	    tmpa = eps * work1[inds + i__];
	}
	if (abs(tmpa) < abs(*mingma1)) {
	    *mingma1 = tmpa;
	    *rone = i__ + 1;
	}
/* L130: */
    }
    *mingma2 = work2[inds + r1b - 1] + work2[indp + r1b - 1];
    if (*mingma2 == 0.f) {
	*mingma2 = eps * work2[inds + r1b - 1];
    }
    *rtwo = r1b;
    i__1 = r2b - 1;
    for (i__ = r1b; i__ <= i__1; ++i__) {
	tmpb = work2[inds + i__] + work2[indp + i__];
	if (tmpb == 0.f) {
	    tmpb = eps * work2[inds + i__];
	}
	if (abs(tmpb) < abs(*mingma2)) {
	    *mingma2 = tmpb;
	    *rtwo = i__ + 1;
	}
/* L140: */
    }

    isup1[1] = *b1;
    isup1[2] = *bn;
    isup2[1] = *b1;
    isup2[2] = *bn;
    z1[*rone] = 1.f;
    z2[*rtwo] = 1.f;
    *ztz1 = 1.f;
    *ztz2 = 1.f;
    i__1 = *b1;
    for (i__ = f2cmax(*rone,*rtwo) - 1; i__ >= i__1; --i__) {
	if (i__ < *rone) {
	    z1[i__] = -(work1[i__] * z1[i__ + 1]);
	    *ztz1 += z1[i__] * z1[i__];
	}
	if (i__ < *rtwo) {
	    z2[i__] = -(work2[i__] * z2[i__ + 1]);
	    *ztz2 += z2[i__] * z2[i__];
	}
/* L150: */
    }
    i__1 = *bn;
    for (i__ = f2cmin(*rone,*rtwo) + 1; i__ <= i__1; ++i__) {
	if (i__ > *rone) {
	    z1[i__] = -(work1[indumn + i__ - 1] * z1[i__ - 1]);
	    *ztz1 += z1[i__] * z1[i__];
	}
	if (i__ > *rtwo) {
	    z2[i__] = -(work2[indumn + i__ - 1] * z2[i__ - 1]);
	    *ztz2 += z2[i__] * z2[i__];
	}
/* L160: */
    }
    return;

L300:
    slar1v_tgk__(rep, n, b1, bn, sigma1, &d__[1], &l[1], &ld[1], &lld[1], 
	    eval1, &gersch[1], &z1[1], ztz1, mingma1, rone, &isup1[1], &work1[
	    1]);
    slar1v_tgk__(rep, n, b1, bn, sigma2, &d__[1], &l[1], &ld[1], &lld[1], 
	    eval2, &gersch[1], &z2[1], ztz2, mingma2, rtwo, &isup2[1], &work2[
	    1]);
    return;
} /* slar1v2_tgk__ */

