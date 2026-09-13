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




/* > \brief \b DPTTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPTTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpttrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpttrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpttrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPTTRF( N, D, E, INFO ) */

/*       INTEGER            INFO, N */
/*       DOUBLE PRECISION   D( * ), E( * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPTTRF computes the L*D*L**T factorization of a real symmetric */
/* > positive definite tridiagonal matrix A.  The factorization may also */
/* > be regarded as having the form A = U**T*D*U. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          On entry, the n diagonal elements of the tridiagonal matrix */
/* >          A.  On exit, the n diagonal elements of the diagonal matrix */
/* >          D from the L*D*L**T factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* >          matrix A.  On exit, the (n-1) subdiagonal elements of the */
/* >          unit bidiagonal factor L from the L*D*L**T factorization of A. */
/* >          E can also be regarded as the superdiagonal of the unit */
/* >          bidiagonal factor U from the U**T*D*U factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -k, the k-th argument had an illegal value */
/* >          > 0: if INFO = k, the leading minor of order k is not */
/* >               positive definite; if k < N, the factorization could not */
/* >               be completed, while if k = N, the factorization was */
/* >               completed, but D(N) <= 0. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doublePTcomputational */

/*  ===================================================================== */
/* Subroutine */ void dpttrf_(integer *n, doublereal *d__, doublereal *e, 
	integer *info)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, i4;
    doublereal ei;
    extern /* Subroutine */ void xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */


/*  ===================================================================== */


/*     Test the input parameters. */

    /* Parameter adjustments */
    --e;
    --d__;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
	i__1 = -(*info);
	xerbla_("DPTTRF", &i__1, (ftnlen)6);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/*     Compute the L*D*L**T (or U**T*D*U) factorization of A. */

    i4 = (*n - 1) % 4;
    i__1 = i4;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (d__[i__] <= 0.) {
	    *info = i__;
	    goto L30;
	}
	ei = e[i__];
	e[i__] = ei / d__[i__];
	d__[i__ + 1] -= e[i__] * ei;
/* L10: */
    }

    i__1 = *n - 4;
    for (i__ = i4 + 1; i__ <= i__1; i__ += 4) {

/*        Drop out of the loop if d(i) <= 0: the matrix is not positive */
/*        definite. */

	if (d__[i__] <= 0.) {
	    *info = i__;
	    goto L30;
	}

/*        Solve for e(i) and d(i+1). */

	ei = e[i__];
	e[i__] = ei / d__[i__];
	d__[i__ + 1] -= e[i__] * ei;

	if (d__[i__ + 1] <= 0.) {
	    *info = i__ + 1;
	    goto L30;
	}

/*        Solve for e(i+1) and d(i+2). */

	ei = e[i__ + 1];
	e[i__ + 1] = ei / d__[i__ + 1];
	d__[i__ + 2] -= e[i__ + 1] * ei;

	if (d__[i__ + 2] <= 0.) {
	    *info = i__ + 2;
	    goto L30;
	}

/*        Solve for e(i+2) and d(i+3). */

	ei = e[i__ + 2];
	e[i__ + 2] = ei / d__[i__ + 2];
	d__[i__ + 3] -= e[i__ + 2] * ei;

	if (d__[i__ + 3] <= 0.) {
	    *info = i__ + 3;
	    goto L30;
	}

/*        Solve for e(i+3) and d(i+4). */

	ei = e[i__ + 3];
	e[i__ + 3] = ei / d__[i__ + 3];
	d__[i__ + 4] -= e[i__ + 3] * ei;
/* L20: */
    }

/*     Check d(n) for positive definiteness. */

    if (d__[*n] <= 0.) {
	*info = *n;
    }

L30:
    return;

/*     End of DPTTRF */

} /* dpttrf_ */

