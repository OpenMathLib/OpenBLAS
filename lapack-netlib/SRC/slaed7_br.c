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

static integer c__2 = 2;
static integer c__1 = 1;
static integer c_n1 = -1;

/* > \brief \b SLAED7_BR used by SSTEBR. */

/*  SLAED7_BR performs one divide-and-conquer merge while propagating */
/*  only the first and last rows of the current local eigenvector block. */
/*  It is called by SLAED0_BR and is not itself a public LAPACK API. */

/*  Serial Reference-LAPACK form: no OpenMP or thread-scaled workspace. */

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

/* Subroutine */ void slaed7_br_(logical *wantq, integer *n, integer *cutpnt, 
	real *d__, integer *indxq, real *rho, real *blo, real *bhi, real *
	work, integer *iwork, integer *kdefl, integer *info)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ void slaed8_br_(integer *, integer *, integer *, 
	    real *, real *, integer *, integer *, real *, integer *, real *, 
	    real *, real *, integer *, real *, integer *, integer *, integer *
	    , real *, integer *, integer *, integer *), slaed9_br_(integer *,
	     integer *, real *, real *, integer *, real *, integer *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, integer *)
	    ;
    integer i__, j, k, n1, n2, iq, iw, iz, iq2, nbr, iscr, itau, inear, iindx,
	     iperm, iwsgn, idlmda, igivcl;
    extern /* Subroutine */ void xerbla_(char *, integer *, ftnlen), slamrg_(
	    integer *, integer *, real *, integer *, integer *, integer *);
    integer igivnm, iindxp, givptr;



/*  ===================================================================== */


    /* Parameter adjustments */
    --iwork;
    --work;
    --bhi;
    --blo;
    --indxq;
    --d__;

    /* Function Body */
    *info = 0;
    *kdefl = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*cutpnt < f2cmin(1,*n) || *cutpnt > *n) {
	*info = -2;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SLAED7_BR", &i__1, (ftnlen)9);
	return;
    }

    if (*n == 0) {
	return;
    }

/*     Workspace layout (14*N doubles): */

/*       Z, DLAMBDA, W                  : N each */
/*       QBR, Q2BR, GIVNUM/reused WPROD : 2*N each */
/*       DELTA, WSGN, TAU, NEAR          : N, N, N, 2*N */

    iz = 1;
    idlmda = iz + *n;
    iw = idlmda + *n;
    iq = iw + *n;
    iq2 = iq + (*n << 1);
    igivnm = iq2 + (*n << 1);
    iscr = igivnm + (*n << 1);
    iwsgn = iscr + *n;
    itau = iwsgn + *n;
    inear = itau + *n;

/*     Integer workspace layout. */

    iindxp = 1;
    iindx = iindxp + *n;
    iperm = iindx + *n;
    igivcl = iperm + *n;

    if (*wantq) {
	nbr = 2;
    } else {
	nbr = 0;
    }

/*     Z is the rank-one update vector.  QBR has two selected rows of */
/*     blockdiag(Q_left,Q_right): parent first boundary and parent last */
/*     boundary, represented in the child eigenbases.  The root merge has */
/*     no parent, so it only needs Z for the secular eigenvalues. */

    if (*wantq) {
	i__1 = *cutpnt;
	for (j = 1; j <= i__1; ++j) {
	    work[iz + j - 1] = bhi[j];
	    work[iq + (j - 1 << 1)] = blo[j];
	    work[iq + (j - 1 << 1) + 1] = 0.f;
/* L20: */
	}
	i__1 = *n;
	for (j = *cutpnt + 1; j <= i__1; ++j) {
	    work[iz + j - 1] = blo[j];
	    work[iq + (j - 1 << 1)] = 0.f;
	    work[iq + (j - 1 << 1) + 1] = bhi[j];
/* L30: */
	}
    } else {
	i__1 = *cutpnt;
	for (j = 1; j <= i__1; ++j) {
	    work[iz + j - 1] = bhi[j];
/* L40: */
	}
	i__1 = *n;
	for (j = *cutpnt + 1; j <= i__1; ++j) {
	    work[iz + j - 1] = blo[j];
/* L50: */
	}
    }

/*     Sort and deflate, applying the same Givens/permutation actions to */
/*     the selected boundary rows when a parent merge will need them. */


    slaed8_br_(&nbr, &k, n, &d__[1], &work[iq], &c__2, &indxq[1], rho, 
	    cutpnt, &work[iz], &work[idlmda], &work[iq2], &c__2, &work[iw], &
	    iwork[iperm], &givptr, &iwork[igivcl], &work[igivnm], &iwork[
	    iindxp], &iwork[iindx], info);
    if (*info != 0) {
	goto L70;
    }
    *kdefl = k;

/*     Solve the secular equation and directly update the selected rows. */
/*     WORK(ISCR) holds one SLAED4 delta vector at a time; no K-by-K */
/*     delta or eigenvector block is materialized. */

    if (k != 0) {
	slaed9_br_(&k, &nbr, &d__[1], &work[iq2], &c__2, &work[iq], &c__2, 
		rho, &work[idlmda], &work[iw], &work[iscr], &work[igivnm], &
		work[iwsgn], &work[itau], &work[inear], info);
	if (*info != 0) {
	    goto L70;
	}

/*        Prepare INDXQ sorting permutation. */

	n1 = k;
	n2 = *n - k;
	slamrg_(&n1, &n2, &d__[1], &c__1, &c_n1, &indxq[1]);
    } else {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    indxq[i__] = i__;
/* L60: */
	}
    }

/*     Persist the parent boundary rows in the current physical D order. */

    if (*wantq) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    blo[j] = work[iq + (j - 1 << 1)];
	    bhi[j] = work[iq + (j - 1 << 1) + 1];
/* L65: */
	}
    }

L70:
    return;

/*     End of SLAED7_BR */

} /* slaed7_br__ */

