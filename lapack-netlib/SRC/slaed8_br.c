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

/* Table of constant values */

static real c_b3 = -1.f;
static integer c__1 = 1;

/* > \brief \b SLAED8_BR selected-row variant of SLAED8. */

/*  SLAED8_BR mirrors SLAED8 deflation and bookkeeping but applies the */
/*  column rotations/permutations to NBR selected rows rather than to a */
/*  full eigenvector matrix.  It is an internal helper for SLAED7_BR. */

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

/* Subroutine */ void slaed8_br_(integer *nbr, integer *k, integer *n, real *
	d__, real *q, integer *ldq, integer *indxq, real *rho, integer *
	cutpnt, real *z__, real *dlambda, real *q2, integer *ldq2, real *w, 
	integer *perm, integer *givptr, integer *givcol, real *givnum, 
	integer *indxp, integer *indx, integer *info)
{
    /* System generated locals */
    integer q_dim1, q_offset, q2_dim1, q2_offset, i__1;
    real r__1;

    /* Local variables */
    real c__;
    integer i__, j;
    real s, t;
    integer k2, n1, n2, jp, n1p1;
    real eps, tau, tol;
    integer jlam, imax, jmax;
    extern /* Subroutine */ void srot_(integer *, real *, integer *, real *, 
	    integer *, real *, real *), sscal_(integer *, real *, real *, 
	    integer *), scopy_(integer *, real *, integer *, real *, integer *
	    );
    extern real slapy2_(real *, real *), slamch_(char *);
    extern /* Subroutine */ void xerbla_(char *, integer *, ftnlen);
    extern integer isamax_(integer *, real *, integer *);
    extern /* Subroutine */ void slamrg_(integer *, integer *, real *, integer 
	    *, integer *, integer *), slacpy_(char *, integer *, integer *, 
	    real *, integer *, real *, integer *);



/*  ===================================================================== */


    /* Parameter adjustments */
    --d__;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --indxq;
    --z__;
    --dlambda;
    q2_dim1 = *ldq2;
    q2_offset = 1 + q2_dim1;
    q2 -= q2_offset;
    --w;
    --perm;
    givcol -= 3;
    givnum -= 3;
    --indxp;
    --indx;

    /* Function Body */
    *info = 0;
    if (*nbr < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ldq < f2cmax(1,*nbr)) {
	*info = -6;
    } else if (*cutpnt < f2cmin(1,*n) || *cutpnt > *n) {
	*info = -9;
    } else if (*ldq2 < f2cmax(1,*nbr)) {
	*info = -13;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SLAED8_BR", &i__1, (ftnlen)9);
	return;
    }

    *givptr = 0;
    if (*n == 0) {
	return;
    }

    n1 = *cutpnt;
    n2 = *n - n1;
    n1p1 = n1 + 1;

    if (*rho < 0.f) {
	sscal_(&n2, &c_b3, &z__[n1p1], &c__1);
    }

    t = 1.f / sqrt(2.f);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	indx[j] = j;
/* L80: */
    }
    sscal_(n, &t, &z__[1], &c__1);
    *rho = (r__1 = *rho * 2.f, abs(r__1));

/*     Sort the eigenvalues into increasing order. */

    i__1 = *n;
    for (i__ = *cutpnt + 1; i__ <= i__1; ++i__) {
	indxq[i__] += *cutpnt;
/* L90: */
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dlambda[i__] = d__[indxq[i__]];
	w[i__] = z__[indxq[i__]];
/* L100: */
    }
    slamrg_(&n1, &n2, &dlambda[1], &c__1, &c__1, &indx[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__[i__] = dlambda[indx[i__]];
	z__[i__] = w[indx[i__]];
/* L110: */
    }

    imax = isamax_(n, &z__[1], &c__1);
    jmax = isamax_(n, &d__[1], &c__1);
    eps = slamch_("Epsilon");
    tol = eps * 8.f * (r__1 = d__[jmax], abs(r__1));

    if (*rho * (r__1 = z__[imax], abs(r__1)) <= tol) {
	*k = 0;
	if (*nbr == 0) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		perm[j] = indxq[indx[j]];
/* L115: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		perm[j] = indxq[indx[j]];
		scopy_(nbr, &q[perm[j] * q_dim1 + 1], &c__1, &q2[j * q2_dim1 
			+ 1], &c__1);
/* L120: */
	    }
	    slacpy_("A", nbr, n, &q2[q2_offset], ldq2, &q[q_offset], ldq);
	}
	return;
    }

    *k = 0;
    k2 = *n + 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (*rho * (r__1 = z__[j], abs(r__1)) <= tol) {
	    --k2;
	    indxp[k2] = j;
	    if (j == *n) {
		goto L170;
	    }
	} else {
	    jlam = j;
	    goto L140;
	}
/* L130: */
    }
L140:
    ++j;
    if (j > *n) {
	goto L160;
    }
    if (*rho * (r__1 = z__[j], abs(r__1)) <= tol) {
	--k2;
	indxp[k2] = j;
    } else {
	s = z__[jlam];
	c__ = z__[j];
	tau = slapy2_(&c__, &s);
	t = d__[j] - d__[jlam];
	c__ /= tau;
	s = -s / tau;
	if ((r__1 = t * c__ * s, abs(r__1)) <= tol) {
	    z__[j] = tau;
	    z__[jlam] = 0.f;

	    ++(*givptr);
	    givcol[(*givptr << 1) + 1] = indxq[indx[jlam]];
	    givcol[(*givptr << 1) + 2] = indxq[indx[j]];
	    givnum[(*givptr << 1) + 1] = c__;
	    givnum[(*givptr << 1) + 2] = s;
	    if (*nbr != 0) {
		srot_(nbr, &q[indxq[indx[jlam]] * q_dim1 + 1], &c__1, &q[
			indxq[indx[j]] * q_dim1 + 1], &c__1, &c__, &s);
	    }
	    t = d__[jlam] * c__ * c__ + d__[j] * s * s;
	    d__[j] = d__[jlam] * s * s + d__[j] * c__ * c__;
	    d__[jlam] = t;
	    --k2;
	    i__ = 1;
L150:
	    if (k2 + i__ <= *n) {
		if (d__[jlam] < d__[indxp[k2 + i__]]) {
		    indxp[k2 + i__ - 1] = indxp[k2 + i__];
		    indxp[k2 + i__] = jlam;
		    ++i__;
		    goto L150;
		} else {
		    indxp[k2 + i__ - 1] = jlam;
		}
	    } else {
		indxp[k2 + i__ - 1] = jlam;
	    }
	    jlam = j;
	} else {
	    ++(*k);
	    w[*k] = z__[jlam];
	    dlambda[*k] = d__[jlam];
	    indxp[*k] = jlam;
	    jlam = j;
	}
    }
    goto L140;
L160:

    ++(*k);
    w[*k] = z__[jlam];
    dlambda[*k] = d__[jlam];
    indxp[*k] = jlam;

L170:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jp = indxp[j];
	dlambda[j] = d__[jp];
	perm[j] = indxq[indx[jp]];
	if (*nbr != 0) {
	    scopy_(nbr, &q[perm[j] * q_dim1 + 1], &c__1, &q2[j * q2_dim1 + 1],
		     &c__1);
	}
/* L180: */
    }

    if (*k < *n) {
	i__1 = *n - *k;
	scopy_(&i__1, &dlambda[*k + 1], &c__1, &d__[*k + 1], &c__1);
	if (*nbr != 0) {
	    i__1 = *n - *k;
	    slacpy_("A", nbr, &i__1, &q2[(*k + 1) * q2_dim1 + 1], ldq2, &q[(*
		    k + 1) * q_dim1 + 1], ldq);
	}
    }

    return;

/*     End of SLAED8_BR */

} /* slaed8_br__ */

