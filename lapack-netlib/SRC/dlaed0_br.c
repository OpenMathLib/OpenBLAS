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
static integer c__14 = 14;
static integer c__1 = 1;

/* > \brief \b DLAED0_BR used by DSTEBR. */

/*  DLAED0_BR computes eigenvalues of a symmetric tridiagonal matrix with */
/*  the boundary-row divide-and-conquer method used by DSTEBR.  Unlike */
/*  DLAED0(ICOMPQ=0), it keeps only first/last boundary rows for each */
/*  active block across merge levels instead of saving complete */
/*  eigenvector blocks. */

/*  Serial Reference-LAPACK form: no OpenMP, environment-variable */
/*  scheduling, or thread-scaled workspace. */

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

/* Subroutine */ void dlaed0_br_(integer *n, doublereal *d__, doublereal *e, 
	doublereal *work, integer *iwork, integer *info)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    extern /* Subroutine */ void dlaed0_merge_pair_(integer *, integer *, 
	    integer *, logical *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *);
    integer i__, j, k;
    extern /* Subroutine */ void dsteqr_boundary_br_(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *);
    integer msd2, smm1, spm1, ibhi, iblo, iscr, kdefl, qbase, merge, wbase, 
	    infom;
    extern /* Subroutine */ void dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    integer indxq, iiwrk;
    logical wantq;
    integer iwbase, nmerge;
    extern /* Subroutine */ void xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    integer submat, subpbs, matsiz, smlsiz, wrkbase;



/*  ===================================================================== */


    /* Parameter adjustments */
    --iwork;
    --work;
    --e;
    --d__;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DLAED0_BR", &i__1, (ftnlen)9);
	return;
    }

    if (*n == 0) {
	return;
    }

    smlsiz = ilaenv_(&c__9, "DLAED0", " ", &c__0, &c__0, &c__0, &c__0, (
	    ftnlen)6, (ftnlen)1);

/*     Determine the size and placement of leaf submatrices. */

    iwork[1] = *n;
    subpbs = 1;
L10:
    if (iwork[subpbs] > smlsiz) {
	for (j = subpbs; j >= 1; --j) {
	    iwork[j * 2] = (iwork[j] + 1) / 2;
	    iwork[(j << 1) - 1] = iwork[j] / 2;
/* L20: */
	}
	subpbs <<= 1;
	goto L10;
    }
    i__1 = subpbs;
    for (j = 2; j <= i__1; ++j) {
	iwork[j] += iwork[j - 1];
/* L30: */
    }

/*     Apply rank-1 cuts. */

    spm1 = subpbs - 1;
    i__1 = spm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	submat = iwork[i__] + 1;
	smm1 = submat - 1;
	d__[smm1] -= (d__1 = e[smm1], abs(d__1));
	d__[submat] -= (d__1 = e[smm1], abs(d__1));
/* L40: */
    }

/*     Workspace layout: */
/*       WORK(IBLO:IBLO+N-1)       first boundary rows */
/*       WORK(IBHI:IBHI+N-1)       last boundary rows */
/*       WORK(ISCR:)               leaf scratch and per-merge scratch */

/*     Integer workspace layout: */
/*       IWORK(1:N)                split endpoints */
/*       IWORK(INDXQ:INDXQ+N-1)   sorting permutations */
/*       IWORK(IIWRK:)             DLAED7_BR scratch */

    iblo = 1;
    ibhi = iblo + *n;
    iscr = ibhi + *n;
    indxq = *n + 1;
    iiwrk = indxq + *n;

/*     Solve leaf eigenproblems and initialize boundary row invariants. */

    i__1 = subpbs;
    for (i__ = 1; i__ <= i__1; ++i__) {
	qbase = iscr;
	wbase = qbase + smlsiz * smlsiz;
	if (i__ == 1) {
	    submat = 1;
	    matsiz = iwork[1];
	} else {
	    submat = iwork[i__ - 1] + 1;
	    matsiz = iwork[i__] - iwork[i__ - 1];
	}

	dsteqr_boundary_br_(&matsiz, &d__[submat], &e[submat], &work[iblo + 
		submat - 1], &work[ibhi + submat - 1], &work[qbase], &smlsiz, 
		&work[wbase], info);
	if (*info != 0) {
	    goto L130;
	}
	k = 1;
	i__2 = iwork[i__];
	for (j = submat; j <= i__2; ++j) {
	    iwork[indxq + j - 1] = k;
	    ++k;
/* L65: */
	}
/* L70: */
    }

/*     Merge adjacent eigensystems.  Boundary rows are updated in place */
/*     for the merged physical D order. */

L80:
    if (subpbs > 1) {
	wantq = subpbs > 2;
	nmerge = subpbs / 2;
	i__1 = nmerge;
	for (merge = 1; merge <= i__1; ++merge) {
	    dlaed0_merge_pair_(&merge, &nmerge, &iwork[1], &wantq, &d__[1], &
		    e[1], &work[1], &iblo, &ibhi, &iscr, &iiwrk, &indxq, &
		    submat, &matsiz, &msd2, &wrkbase, &iwbase, &kdefl, &infom,
		     &c__14);
	    if (infom != 0) {
		*info = infom;
		goto L130;
	    }
/* L90: */
	}

	i__1 = nmerge;
	for (merge = 1; merge <= i__1; ++merge) {
	    iwork[merge] = iwork[merge * 2];
/* L95: */
	}
	subpbs /= 2;
	goto L80;
    }

/*     Reorder eigenvalues from the final physical merge order. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iwork[indxq + i__ - 1];
	work[iscr + i__ - 1] = d__[j];
/* L100: */
    }
    dcopy_(n, &work[iscr], &c__1, &d__[1], &c__1);
    goto L140;

L130:
    *info = submat * (*n + 1) + submat + matsiz - 1;

L140:
    return;

/*     End of DLAED0_BR */

} /* dlaed0_br__ */

