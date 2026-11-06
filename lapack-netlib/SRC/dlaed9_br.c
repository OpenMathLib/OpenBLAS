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

/* > \brief \b DLAED9_BR used by DSTEBR. */

/*  DLAED9_BR finds the roots of the secular equation for one */
/*  boundary-row DC merge and directly updates a small set of selected */
/*  rows.  It mirrors DLAED9's root solve and eigenvector normalization, */
/*  but avoids materializing the full K-by-K secular eigenvector block. */
/*  DLAED4_BR returns high-relative-accuracy DELTA values for each root; */
/*  consume those DELTA values immediately to form WPROD, following the */
/*  same numerical idea as LAPACK's SVD divide-and-conquer helpers. */

/*  Serial Reference-LAPACK form: no OpenMP or thread-private */
/*  accumulators. */

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

/* Subroutine */ void dlaed9_br_(integer *k, integer *nbr, doublereal *d__, 
	doublereal *qin, integer *ldqin, doublereal *qout, integer *ldqout, 
	doublereal *rho, doublereal *dlambda, doublereal *w, doublereal *
	delta, doublereal *wprod, doublereal *wsgn, doublereal *tau, 
	doublereal *near__, integer *info)
{
    /* System generated locals */
    integer qin_dim1, qin_offset, qout_dim1, qout_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    extern /* Subroutine */ void dlaed4_br_(integer *, integer *, doublereal *
	    , doublereal *, doublereal *, doublereal *, doublereal *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *);
    integer i__, j, r__, n1, n2;
    doublereal dlt, acc1, acc2;
    integer iorg;
    doublereal temp;
    extern /* Subroutine */ void dlaed4_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    , xerbla_(char *, integer *, ftnlen);
    doublereal rowsum;



/*  ===================================================================== */


    /* Parameter adjustments */
    --d__;
    qin_dim1 = *ldqin;
    qin_offset = 1 + qin_dim1;
    qin -= qin_offset;
    qout_dim1 = *ldqout;
    qout_offset = 1 + qout_dim1;
    qout -= qout_offset;
    --dlambda;
    --w;
    --delta;
    --wprod;
    --wsgn;
    --tau;
    near__ -= 3;

    /* Function Body */
    *info = 0;
    if (*k < 0) {
	*info = -1;
    } else if (*nbr < 0) {
	*info = -2;
    } else if (*ldqin < f2cmax(1,*nbr)) {
	*info = -5;
    } else if (*ldqout < f2cmax(1,*nbr)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DLAED9_BR", &i__1, (ftnlen)9);
	return;
    }

    if (*k == 0) {
	return;
    }

/*     Root-only final merge: no parent boundary rows will be observed by */
/*     a higher merge level, so updated W and secular vectors are not */
/*     needed. */

    if (*nbr == 0) {
	i__1 = *k;
	for (j = 1; j <= i__1; ++j) {
	    dlaed4_(k, &j, &dlambda[1], &w[1], &delta[1], rho, &d__[j], info);
	    if (*info != 0) {
		goto L170;
	    }
/* L5: */
	}
	goto L170;
    }

/*     For K = 1 or 2, DLAED4/DLAED5 returns normalized eigenvectors */
/*     directly in DELTA.  Keep this path identical to DLAED9's */
/*     special-case copy to S. */

    if (*k == 1 || *k == 2) {
	i__1 = *k;
	for (j = 1; j <= i__1; ++j) {
	    dlaed4_(k, &j, &dlambda[1], &w[1], &delta[1], rho, &d__[j], info);
	    if (*info != 0) {
		goto L170;
	    }
	    if (*nbr == 2) {
		acc1 = 0.;
		acc2 = 0.;
		i__2 = *k;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    acc1 += qin[i__ * qin_dim1 + 1] * delta[i__];
		    acc2 += qin[i__ * qin_dim1 + 2] * delta[i__];
/* L10: */
		}
		qout[j * qout_dim1 + 1] = acc1;
		qout[j * qout_dim1 + 2] = acc2;
	    } else {
		i__2 = *nbr;
		for (r__ = 1; r__ <= i__2; ++r__) {
		    rowsum = 0.;
		    i__3 = *k;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			rowsum += qin[r__ + i__ * qin_dim1] * delta[i__];
/* L15: */
		    }
		    qout[r__ + j * qout_dim1] = rowsum;
/* L20: */
		}
	    }
/* L30: */
	}
	goto L170;
    }

/*     Preserve the update vector, solve all secular roots, then */
/*     accumulate the product needed to reconstruct DLAED9's updated */
/*     W.  The high-relative-accuracy DELTA values produced by */
/*     DLAED4_BR are consumed immediately rather than reconstructed */
/*     later from the compact origin/tau representation. */

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wsgn[i__] = w[i__];
/* L40: */
    }

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wprod[i__] = 1.;
/* L41: */
    }

    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	dlaed4_br_(k, &j, &dlambda[1], &wsgn[1], &delta[1], rho, &d__[j], &
		iorg, &tau[j], &near__[(j << 1) + 1], &near__[(j << 1) + 2], 
		info);
	if (*info != 0) {
	    goto L170;
	}

	wprod[j] *= delta[j];
	i__2 = j - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    wprod[i__] *= delta[i__] / (dlambda[i__] - dlambda[j]);
/* L42: */
	}
	i__2 = *k;
	for (i__ = j + 1; i__ <= i__2; ++i__) {
	    wprod[i__] *= delta[i__] / (dlambda[i__] - dlambda[j]);
/* L43: */
	}
/* L44: */
    }

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__1 = sqrt(-wprod[i__]);
	w[i__] = d_sign(&d__1, &wsgn[i__]);
/* L80: */
    }

/*     Stream each secular eigenvector column through the selected rows. */
/*     Reconstruct the DELTA column from the compact origin/tau state, */
/*     using cached pole-adjacent values where cancellation is largest. */

    if (*nbr == 2) {
	i__1 = *k;
	for (j = 1; j <= i__1; ++j) {
	    if (j == *k) {
		iorg = *k;
		n1 = *k - 1;
		n2 = *k;
	    } else if (tau[j] >= 0.) {
		iorg = j;
		n1 = j;
		n2 = j + 1;
	    } else {
		iorg = j + 1;
		n1 = j;
		n2 = j + 1;
	    }
	    temp = 0.;
	    acc1 = 0.;
	    acc2 = 0.;

	    i__2 = n1 - 1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		dlt = w[i__] / (dlambda[i__] - dlambda[iorg] - tau[j]);
		temp += dlt * dlt;
		acc1 += qin[i__ * qin_dim1 + 1] * dlt;
		acc2 += qin[i__ * qin_dim1 + 2] * dlt;
/* L85: */
	    }

	    dlt = w[n1] / near__[(j << 1) + 1];
	    temp += dlt * dlt;
	    acc1 += qin[n1 * qin_dim1 + 1] * dlt;
	    acc2 += qin[n1 * qin_dim1 + 2] * dlt;

	    i__2 = n2 - 1;
	    for (i__ = n1 + 1; i__ <= i__2; ++i__) {
		dlt = w[i__] / (dlambda[i__] - dlambda[iorg] - tau[j]);
		temp += dlt * dlt;
		acc1 += qin[i__ * qin_dim1 + 1] * dlt;
		acc2 += qin[i__ * qin_dim1 + 2] * dlt;
/* L86: */
	    }

	    dlt = w[n2] / near__[(j << 1) + 2];
	    temp += dlt * dlt;
	    acc1 += qin[n2 * qin_dim1 + 1] * dlt;
	    acc2 += qin[n2 * qin_dim1 + 2] * dlt;

	    i__2 = *k;
	    for (i__ = n2 + 1; i__ <= i__2; ++i__) {
		dlt = w[i__] / (dlambda[i__] - dlambda[iorg] - tau[j]);
		temp += dlt * dlt;
		acc1 += qin[i__ * qin_dim1 + 1] * dlt;
		acc2 += qin[i__ * qin_dim1 + 2] * dlt;
/* L87: */
	    }

	    temp = sqrt(temp);
	    qout[j * qout_dim1 + 1] = acc1 / temp;
	    qout[j * qout_dim1 + 2] = acc2 / temp;
/* L115: */
	}
	goto L170;
    }

    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	if (j == *k) {
	    iorg = *k;
	} else if (tau[j] >= 0.) {
	    iorg = j;
	} else {
	    iorg = j + 1;
	}
	temp = 0.;
	i__2 = *nbr;
	for (r__ = 1; r__ <= i__2; ++r__) {
	    qout[r__ + j * qout_dim1] = 0.;
/* L90: */
	}
	i__2 = *k;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dlt = dlambda[i__] - dlambda[iorg] - tau[j];
	    if (j == *k) {
		if (i__ == *k - 1) {
		    dlt = near__[(j << 1) + 1];
		}
		if (i__ == *k) {
		    dlt = near__[(j << 1) + 2];
		}
	    } else {
		if (i__ == j) {
		    dlt = near__[(j << 1) + 1];
		}
		if (i__ == j + 1) {
		    dlt = near__[(j << 1) + 2];
		}
	    }
	    dlt = w[i__] / dlt;
	    temp += dlt * dlt;
	    i__3 = *nbr;
	    for (r__ = 1; r__ <= i__3; ++r__) {
		qout[r__ + j * qout_dim1] += qin[r__ + i__ * qin_dim1] * dlt;
/* L95: */
	    }
/* L100: */
	}
	temp = sqrt(temp);
	i__2 = *nbr;
	for (r__ = 1; r__ <= i__2; ++r__) {
	    qout[r__ + j * qout_dim1] /= temp;
/* L120: */
	}
/* L130: */
    }

L170:
    return;

/*     End of DLAED9_BR */

} /* dlaed9_br__ */

