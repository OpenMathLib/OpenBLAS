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




/* > \brief \b DLASQ5 computes one dqds transform in ping-pong form. Used by sbdsqr and sstegr. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASQ5 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq5.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasq5.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasq5.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASQ5( I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, DN, */
/*                          DNM1, DNM2, IEEE, EPS ) */

/*       LOGICAL            IEEE */
/*       INTEGER            I0, N0, PP */
/*       DOUBLE PRECISION   DMIN, DMIN1, DMIN2, DN, DNM1, DNM2, TAU, SIGMA, EPS */
/*       DOUBLE PRECISION   Z( * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASQ5 computes one dqds transform in ping-pong form, one */
/* > version for IEEE machines another for non IEEE machines. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] I0 */
/* > \verbatim */
/* >          I0 is INTEGER */
/* >        First index. */
/* > \endverbatim */
/* > */
/* > \param[in] N0 */
/* > \verbatim */
/* >          N0 is INTEGER */
/* >        Last index. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension ( 4*N ) */
/* >        Z holds the qd array. EMIN is stored in Z(4*N0) to avoid */
/* >        an extra argument. */
/* > \endverbatim */
/* > */
/* > \param[in] PP */
/* > \verbatim */
/* >          PP is INTEGER */
/* >        PP=0 for ping, PP=1 for pong. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is DOUBLE PRECISION */
/* >        This is the shift. */
/* > \endverbatim */
/* > */
/* > \param[in] SIGMA */
/* > \verbatim */
/* >          SIGMA is DOUBLE PRECISION */
/* >        This is the accumulated shift up to this step. */
/* > \endverbatim */
/* > */
/* > \param[out] DMIN */
/* > \verbatim */
/* >          DMIN is DOUBLE PRECISION */
/* >        Minimum value of d. */
/* > \endverbatim */
/* > */
/* > \param[out] DMIN1 */
/* > \verbatim */
/* >          DMIN1 is DOUBLE PRECISION */
/* >        Minimum value of d, excluding D( N0 ). */
/* > \endverbatim */
/* > */
/* > \param[out] DMIN2 */
/* > \verbatim */
/* >          DMIN2 is DOUBLE PRECISION */
/* >        Minimum value of d, excluding D( N0 ) and D( N0-1 ). */
/* > \endverbatim */
/* > */
/* > \param[out] DN */
/* > \verbatim */
/* >          DN is DOUBLE PRECISION */
/* >        d(N0), the last value of d. */
/* > \endverbatim */
/* > */
/* > \param[out] DNM1 */
/* > \verbatim */
/* >          DNM1 is DOUBLE PRECISION */
/* >        d(N0-1). */
/* > \endverbatim */
/* > */
/* > \param[out] DNM2 */
/* > \verbatim */
/* >          DNM2 is DOUBLE PRECISION */
/* >        d(N0-2). */
/* > \endverbatim */
/* > */
/* > \param[in] IEEE */
/* > \verbatim */
/* >          IEEE is LOGICAL */
/* >        Flag for IEEE or non IEEE arithmetic. */
/* > \endverbatim */
/* > */
/* > \param[in] EPS */
/* > \verbatim */
/* >          EPS is DOUBLE PRECISION */
/* >        This is the value of epsilon used. */
/* > \endverbatim */
/* > */
/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2017 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ void dlasq5_(integer *i0, integer *n0, doublereal *z__, 
	integer *pp, doublereal *tau, doublereal *sigma, doublereal *dmin__, 
	doublereal *dmin1, doublereal *dmin2, doublereal *dn, doublereal *
	dnm1, doublereal *dnm2, logical *ieee, doublereal *eps)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    doublereal emin, temp, d__;
    integer j4, j4p2;
    doublereal dthresh;


/*  -- LAPACK computational routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */


/*  ===================================================================== */


    /* Parameter adjustments */
    --z__;

    /* Function Body */
    if (*n0 - *i0 - 1 <= 0) {
	return;
    }

    dthresh = *eps * (*sigma + *tau);
    if (*tau < dthresh * .5) {
	*tau = 0.;
    }
    if (*tau != 0.) {
	j4 = (*i0 << 2) + *pp - 3;
	emin = z__[j4 + 4];
	d__ = z__[j4] - *tau;
	*dmin__ = d__;
	*dmin1 = -z__[j4];

	if (*ieee) {

/*        Code for IEEE arithmetic. */

	    if (*pp == 0) {
		i__1 = (*n0 - 3) << 2;
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		    z__[j4 - 2] = d__ + z__[j4 - 1];
		    temp = z__[j4 + 1] / z__[j4 - 2];
		    d__ = d__ * temp - *tau;
		    *dmin__ = f2cmin(*dmin__,d__);
		    z__[j4] = z__[j4 - 1] * temp;
/* Computing MIN */
		    d__1 = z__[j4];
		    emin = f2cmin(d__1,emin);
/* L10: */
		}
	    } else {
		i__1 = (*n0 - 3) << 2;
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		    z__[j4 - 3] = d__ + z__[j4];
		    temp = z__[j4 + 2] / z__[j4 - 3];
		    d__ = d__ * temp - *tau;
		    *dmin__ = f2cmin(*dmin__,d__);
		    z__[j4 - 1] = z__[j4] * temp;
/* Computing MIN */
		    d__1 = z__[j4 - 1];
		    emin = f2cmin(d__1,emin);
/* L20: */
		}
	    }

/*        Unroll last two steps. */

	    *dnm2 = d__;
	    *dmin2 = *dmin__;
	    j4 = ((*n0 - 2) << 2) - *pp;
	    j4p2 = j4 + (*pp << 1) - 1;
	    z__[j4 - 2] = *dnm2 + z__[j4p2];
	    z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	    *dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
	    *dmin__ = f2cmin(*dmin__,*dnm1);

	    *dmin1 = *dmin__;
	    j4 += 4;
	    j4p2 = j4 + (*pp << 1) - 1;
	    z__[j4 - 2] = *dnm1 + z__[j4p2];
	    z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	    *dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
	    *dmin__ = f2cmin(*dmin__,*dn);

	} else {

/*        Code for non IEEE arithmetic. */

	    if (*pp == 0) {
		i__1 = (*n0 - 3) << 2;
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		    z__[j4 - 2] = d__ + z__[j4 - 1];
		    if (d__ < 0.) {
			return;
		    } else {
			z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
			d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]) - *tau;
		    }
		    *dmin__ = f2cmin(*dmin__,d__);
/* Computing MIN */
		    d__1 = emin, d__2 = z__[j4];
		    emin = f2cmin(d__1,d__2);
/* L30: */
		}
	    } else {
		i__1 = (*n0 - 3) << 2;
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		    z__[j4 - 3] = d__ + z__[j4];
		    if (d__ < 0.) {
			return;
		    } else {
			z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
			d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]) - *tau;
		    }
		    *dmin__ = f2cmin(*dmin__,d__);
/* Computing MIN */
		    d__1 = emin, d__2 = z__[j4 - 1];
		    emin = f2cmin(d__1,d__2);
/* L40: */
		}
	    }

/*        Unroll last two steps. */

	    *dnm2 = d__;
	    *dmin2 = *dmin__;
	    j4 = ((*n0 - 2) << 2) - *pp;
	    j4p2 = j4 + (*pp << 1) - 1;
	    z__[j4 - 2] = *dnm2 + z__[j4p2];
	    if (*dnm2 < 0.) {
		return;
	    } else {
		z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
		*dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
	    }
	    *dmin__ = f2cmin(*dmin__,*dnm1);

	    *dmin1 = *dmin__;
	    j4 += 4;
	    j4p2 = j4 + (*pp << 1) - 1;
	    z__[j4 - 2] = *dnm1 + z__[j4p2];
	    if (*dnm1 < 0.) {
		return;
	    } else {
		z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
		*dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
	    }
	    *dmin__ = f2cmin(*dmin__,*dn);

	}
    } else {
/*     This is the version that sets d's to zero if they are small enough */
	j4 = (*i0 << 2) + *pp - 3;
	emin = z__[j4 + 4];
	d__ = z__[j4] - *tau;
	*dmin__ = d__;
	*dmin1 = -z__[j4];
	if (*ieee) {

/*     Code for IEEE arithmetic. */

	    if (*pp == 0) {
		i__1 = (*n0 - 3) << 2;
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		    z__[j4 - 2] = d__ + z__[j4 - 1];
		    temp = z__[j4 + 1] / z__[j4 - 2];
		    d__ = d__ * temp - *tau;
		    if (d__ < dthresh) {
			d__ = 0.;
		    }
		    *dmin__ = f2cmin(*dmin__,d__);
		    z__[j4] = z__[j4 - 1] * temp;
/* Computing MIN */
		    d__1 = z__[j4];
		    emin = f2cmin(d__1,emin);
/* L50: */
		}
	    } else {
		i__1 = (*n0 - 3) << 2;
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		    z__[j4 - 3] = d__ + z__[j4];
		    temp = z__[j4 + 2] / z__[j4 - 3];
		    d__ = d__ * temp - *tau;
		    if (d__ < dthresh) {
			d__ = 0.;
		    }
		    *dmin__ = f2cmin(*dmin__,d__);
		    z__[j4 - 1] = z__[j4] * temp;
/* Computing MIN */
		    d__1 = z__[j4 - 1];
		    emin = f2cmin(d__1,emin);
/* L60: */
		}
	    }

/*     Unroll last two steps. */

	    *dnm2 = d__;
	    *dmin2 = *dmin__;
	    j4 = ((*n0 - 2) << 2) - *pp;
	    j4p2 = j4 + (*pp << 1) - 1;
	    z__[j4 - 2] = *dnm2 + z__[j4p2];
	    z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	    *dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
	    *dmin__ = f2cmin(*dmin__,*dnm1);

	    *dmin1 = *dmin__;
	    j4 += 4;
	    j4p2 = j4 + (*pp << 1) - 1;
	    z__[j4 - 2] = *dnm1 + z__[j4p2];
	    z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
	    *dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
	    *dmin__ = f2cmin(*dmin__,*dn);

	} else {

/*     Code for non IEEE arithmetic. */

	    if (*pp == 0) {
		i__1 = (*n0 - 3) << 2;
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		    z__[j4 - 2] = d__ + z__[j4 - 1];
		    if (d__ < 0.) {
			return;
		    } else {
			z__[j4] = z__[j4 + 1] * (z__[j4 - 1] / z__[j4 - 2]);
			d__ = z__[j4 + 1] * (d__ / z__[j4 - 2]) - *tau;
		    }
		    if (d__ < dthresh) {
			d__ = 0.;
		    }
		    *dmin__ = f2cmin(*dmin__,d__);
/* Computing MIN */
		    d__1 = emin, d__2 = z__[j4];
		    emin = f2cmin(d__1,d__2);
/* L70: */
		}
	    } else {
		i__1 = (*n0 - 3) << 2;
		for (j4 = *i0 << 2; j4 <= i__1; j4 += 4) {
		    z__[j4 - 3] = d__ + z__[j4];
		    if (d__ < 0.) {
			return;
		    } else {
			z__[j4 - 1] = z__[j4 + 2] * (z__[j4] / z__[j4 - 3]);
			d__ = z__[j4 + 2] * (d__ / z__[j4 - 3]) - *tau;
		    }
		    if (d__ < dthresh) {
			d__ = 0.;
		    }
		    *dmin__ = f2cmin(*dmin__,d__);
/* Computing MIN */
		    d__1 = emin, d__2 = z__[j4 - 1];
		    emin = f2cmin(d__1,d__2);
/* L80: */
		}
	    }

/*     Unroll last two steps. */

	    *dnm2 = d__;
	    *dmin2 = *dmin__;
	    j4 = ((*n0 - 2) << 2) - *pp;
	    j4p2 = j4 + (*pp << 1) - 1;
	    z__[j4 - 2] = *dnm2 + z__[j4p2];
	    if (*dnm2 < 0.) {
		return;
	    } else {
		z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
		*dnm1 = z__[j4p2 + 2] * (*dnm2 / z__[j4 - 2]) - *tau;
	    }
	    *dmin__ = f2cmin(*dmin__,*dnm1);

	    *dmin1 = *dmin__;
	    j4 += 4;
	    j4p2 = j4 + (*pp << 1) - 1;
	    z__[j4 - 2] = *dnm1 + z__[j4p2];
	    if (*dnm1 < 0.) {
		return;
	    } else {
		z__[j4] = z__[j4p2 + 2] * (z__[j4p2] / z__[j4 - 2]);
		*dn = z__[j4p2 + 2] * (*dnm1 / z__[j4 - 2]) - *tau;
	    }
	    *dmin__ = f2cmin(*dmin__,*dn);

	}
    }

    z__[j4 + 2] = *dn;
    z__[(*n0 << 2) - *pp] = emin;
    return;

/*     End of DLASQ5 */

} /* dlasq5_ */

