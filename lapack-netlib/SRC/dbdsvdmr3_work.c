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

static doublereal c_b7 = 1.;
static integer c__1 = 1;
static logical c_false = FALSE_;
static logical c_true = TRUE_;

/* > \brief \b DBDSVDMR3_WORK computational core of DBDSVDMR3. */
/* > */
/* > \ingroup bdsvdmr3 */
/* > */
/* > \par Purpose: */
/* > ============= */
/* > */
/* > \verbatim */
/* > computational core of DBDSVDMR3: preprocessing, dqds singular values, MRRR eigenvectors of the Golub-Kaha
n matrix, extraction of singular vectors. */
/* > All routines in this file are auxiliary to DBDSVDMR3 and derive */
/* > from the LAPACK 3.0 (1999) MRRR kernels DLARRV, DLAR1V, DLARRB and */
/* > DLARRF rather than from the current DSTEMR kernels. */
/* > \endverbatim */

/* Subroutine */ void dbdsvdmr3_work__(char *jobz, char *uplo, integer *n, 
	doublereal *d__, doublereal *e, doublereal *s, doublereal *u, integer 
	*ldu, doublereal *vt, integer *ldvt, integer *m, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info)
{
    /* System generated locals */
    integer u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    integer i__, j, k, l, kb, iq;
    doublereal ze, zg;
    extern /* Subroutine */ void dlarrv_tgk__(integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *);
    doublereal rt2, scl, eps, tol;
    integer kzu, kzv, idbl, iebl, bdim, tbeg, jblk, tend;
    doublereal bnrm;
    integer mpos, npos;
    doublereal znrm;
    integer nzsv, inde2;
    logical diagb;
    extern /* Subroutine */ void dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    integer indbl;
    extern logical lsame_(char *, char *);
    integer iinfo, indwk;
    doublereal sfmin;
    integer lwmin;
    logical upper, wantz;
    extern /* Subroutine */ void dlasq1_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *);
    extern doublereal dlamch_(char *);
    extern /* Subroutine */ void dbdtgk_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    integer indgsc;
    extern /* Subroutine */ void xerbla_(char *, integer *, ftnlen);
    integer indtgk, indiwk, indisp;
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *);
    doublereal thresh;
    extern /* Subroutine */ void dbsort_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, logical *);
    integer liwmin, indixw, nsplit;
    logical lquery;
    integer indisup;


/*  -- New driver for the bidiagonal SVD via TGK-rooted MR^3 -- */
/*     Uses the LAPACK 3.0-lineage MRRR tree code below the root via */
/*     DLARRV_TGK.  The bundled DLARRB has a bounded-progress guard so an */
/*     invalid bracket radius returns INFO instead of looping forever. */
/*     Standard LAPACK dependencies include DLASQ1, DLAMCH, DLANST, */
/*     DSCAL, and DSWAP. */

/*  Purpose */
/*  ======= */

/*  DBDSVDMR3 computes the SVD of a real N-by-N bidiagonal matrix */
/*  B = U * S * VT, with S diagonal (the singular values, ascending) and */
/*  U, VT orthonormal.  Method: */
/*    1. Build the Tridiagonal Golub-Kahan (TGK) matrix of order 2N; its */
/*       off-diagonals beta interleave the bidiagonal entries (d,e). */
/*    2. Split the TGK wherever any beta is negligible -- this covers both */
/*       the bidiagonal off-diagonals e_i (even positions) and zero */
/*       diagonal entries d_i (odd positions), each of which makes the TGK */
/*       reducible. */
/*    3. Per TGK block, compute the positive eigenvalues (= singular */
/*       values) with DLASQ1 (dqds); an odd-dimension block is a */
/*       rectangular bidiagonal handled by appending a zero diagonal. */
/*    4. Compute the positive eigenvectors by MR^3 rooted at the TGK */
/*       (DLARRV_TGK) and de-interleave each into its left/right singular */
/*       vectors. */
/*    5. Each odd-dimension block also carries one structural ZERO */
/*       eigenvalue; its one-sided null vector (closed form) supplies a */
/*       singular vector for a zero singular value. */

/*  Arguments  (see below; UPLO='U'/'L'; JOBZ='N'/'V'; RANGE = all) */
/*  ========= */

/*  JOBZ   (input) CHARACTER*1   'N' values only; 'V' values + vectors. */
/*  UPLO   (input) CHARACTER*1   'U' upper bidiagonal; 'L' lower. */
/*  N      (input) INTEGER       Order of B.  N >= 0. */
/*  D      (in/out) DOUBLE (N)   Diagonal of B (overwritten). */
/*  E      (in/out) DOUBLE (N)   Off-diagonal in 1..N-1 (overwritten). */
/*  S      (output) DOUBLE (N)   Singular values, ascending. */
/*  U      (output) DOUBLE (LDU,N)   Left singular vectors (col i <-> S(i)). */
/*  LDU    (input) INTEGER */
/*  VT     (output) DOUBLE (LDVT,N)  Right singular vectors, by ROWS. */
/*  LDVT   (input) INTEGER */
/*  M      (output) INTEGER      Number of singular values found (= N). */
/*  WORK   (workspace) DOUBLE (LWORK).  LWORK >= f2cmax(1, 2*N*N + 34*N). */
/*  LWORK  (input) INTEGER       -1 => workspace query. */
/*  IWORK  (workspace) INTEGER (LIWORK).  LIWORK >= f2cmax(1, 20*N). */
/*  LIWORK (input) INTEGER       -1 => workspace query. */
/*  INFO   (output) INTEGER  0 ok; <0 illegal arg; 1 DLASQ1; 2 DLARRV_TGK. */

/*  ===================================================================== */


    /* Parameter adjustments */
    --d__;
    --e;
    --s;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --work;
    --iwork;

    /* Function Body */
    wantz = lsame_(jobz, "V");
    upper = lsame_(uplo, "U");
    lquery = *lwork == -1 || *liwork == -1;
/* Computing MAX */
    i__1 = 1, i__2 = (*n << 1) * *n + *n * 34;
    lwmin = f2cmax(i__1,i__2);
/* Computing MAX */
    i__1 = 1, i__2 = *n * 24;
    liwmin = f2cmax(i__1,i__2);

    *info = 0;
    if (! (wantz || lsame_(jobz, "N"))) {
	*info = -1;
    } else if (! (upper || lsame_(uplo, "L"))) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ldu < 1 || wantz && *ldu < *n) {
	*info = -8;
    } else if (*ldvt < 1 || wantz && *ldvt < *n) {
	*info = -10;
    } else if (*lwork < lwmin && ! lquery) {
	*info = -13;
    } else if (*liwork < liwmin && ! lquery) {
	*info = -15;
    }
    if (*info == 0) {
	work[1] = (doublereal) lwmin;
	iwork[1] = liwmin;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DBDSVDMR3_WORK", &i__1, (ftnlen)14);
	return;
    } else if (lquery) {
	return;
    }

    *m = 0;
    if (*n == 0) {
	return;
    }
    if (*n == 1) {
	*m = 1;
	s[1] = abs(d__[1]);
	if (wantz) {
	    u[u_dim1 + 1] = d_sign(&c_b7, &d__[1]);
	    vt[vt_dim1 + 1] = 1.;
	}
	return;
    }

/*     An exactly diagonal bidiagonal has an exact canonical SVD.  Taking */
/*     it through the TGK path would normalize 2-by-2 eigenvectors by */
/*     1/sqrt(2) and then multiply them by sqrt(2), needlessly turning */
/*     exact identity/permutation vectors into one-ULP approximations. */
/*     Besides removing that roundoff, this path avoids the MR^3 setup for */
/*     the diagonal cases that occur in the LAPACK test suite. */

    diagb = TRUE_;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (e[i__] != 0.) {
	    diagb = FALSE_;
	}
/* L5: */
    }
    if (diagb) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s[i__] = (d__1 = d__[i__], abs(d__1));
/* L8: */
	}
	if (wantz) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    u[i__ + j * u_dim1] = 0.;
		    vt[i__ + j * vt_dim1] = 0.;
/* L10: */
		}
		u[j + j * u_dim1] = d_sign(&c_b7, &d__[j]);
		vt[j + j * vt_dim1] = 1.;
/* L12: */
	    }
	}
	dbsort_(n, &s[1], &u[u_offset], ldu, &vt[vt_offset], ldvt, &wantz);
	*m = *n;
	return;
    }

    eps = dlamch_("Precision");
    sfmin = dlamch_("Safe minimum");
    rt2 = sqrt(2.);

/*     Workspace map (DOUBLE):  INDTGK 2N (TGK diag=0); INDE2 2N (TGK */
/*     off-diagonals beta); INDGSC 4N (Gerschgorin); INDWK = scratch / */
/*     eigenvector array Z (LDZ=2N) / DLARRV_TGK WORK. */
/*     INTEGER:  INDISP 2N (ISPLIT); INDBL 2N (IBLOCK); INDIXW 2N (INDEXW); */
/*     INDISUP 2N (ISUPPZ scratch); INDIWK (DLARRV_TGK IWORK). */

    indtgk = 1;
    inde2 = indtgk + (*n << 1);
    indgsc = inde2 + (*n << 1);
    indwk = indgsc + (*n << 2);

    indisp = 1;
    indbl = indisp + (*n << 1);
    indixw = indbl + (*n << 1);
    indisup = indixw + (*n << 1);
    indiwk = indisup + (*n << 1);

/*     Scale B so that its norm is in a safe range. */

    bnrm = dlanst_("M", n, &d__[1], &e[1]);
    scl = 1.;
    if (bnrm > 0. && bnrm < sqrt(sfmin)) {
	scl = sqrt(sfmin) / bnrm;
    } else if (bnrm > 1. / sqrt(sfmin)) {
	scl = 1. / sqrt(sfmin) / bnrm;
    }
    if (scl != 1.) {
	dscal_(n, &scl, &d__[1], &c__1);
	i__1 = *n - 1;
	dscal_(&i__1, &scl, &e[1], &c__1);
	bnrm *= scl;
    }

/*     Random relative perturbation (Dhillon, Parlett & Voemel 2005, */
/*     sec. 5.4) is NOT applied blindly to the root here.  A blanket */
/*     O(n*eps) perturbation jolts every singular value by ~n*eps, which */
/*     OVERWRITES the genuine fine structure of clusters whose true */
/*     internal gaps are smaller than n*eps (e.g. the Lipshitz matrices, */
/*     whose ~800-fold cluster at sigma~1 has real relative spacing */
/*     ~7e-15 < n*eps): there the perturbation destroyed orthogonality. */
/*     Instead the perturbation is applied SELECTIVELY by DLARRF_TGK, and */
/*     only to a child RRR whose cluster contains EXACTLY-equal */
/*     eigenvalues (true ties, e.g. from identical glued blocks) that */
/*     recursive shifting alone can never separate; the fixed seed is */
/*     set inside DLARRV_TGK. */


/*     Build the TGK matrix of order 2N from the (perturbed) bidiagonal. */
/*     WORK(INDE2 .. INDE2+2N-2) holds beta_1..beta_{2N-1}, where */
/*     beta_{2k-1}=D(k), beta_{2k}=E(k). */

    dbdtgk_(uplo, n, &d__[1], &e[1], &work[indtgk], &work[inde2], &work[
	    indgsc]);

/*     Split the TGK wherever an off-diagonal beta_j is negligible.  This */
/*     handles both the bidiagonal off-diagonals e_i (even j) and zero */
/*     diagonal entries d_i (odd j).  ISPLIT lists the block-end TGK rows. */

    thresh = eps * bnrm;
    nsplit = 0;
    i__1 = (*n << 1) - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if ((d__1 = work[inde2 + i__ - 1], abs(d__1)) <= thresh) {
	    work[inde2 + i__ - 1] = 0.;
	    ++nsplit;
	    iwork[indisp + nsplit - 1] = i__;
	}
/* L20: */
    }
    ++nsplit;
    iwork[indisp + nsplit - 1] = *n << 1;

/*     For each TGK block compute the NPOS positive eigenvalues */
/*     (= singular values) with DLASQ1.  Block of dimension BDIM has its */
/*     beta's at WORK(INDE2 + TBEG-1 .. INDE2 + TEND-2); the sub-bidiagonal */
/*     is d~(l)=beta(local 2l-1), e~(l)=beta(local 2l). */
/*        even BDIM=2k : square k-by-k bidiagonal -> k positive sv; */
/*        odd  BDIM=2k+1: rectangular k-by-(k+1); append a zero diagonal to */
/*                        make it (k+1)-square -> k positive sv + 1 zero. */
/*     The NPOS positives sit in the upper part of the block's BDIM-point */
/*     spectrum, so INDEXW(j-th positive) = (BDIM-NPOS)+j. */

    idbl = indwk;
    iebl = indwk + *n + 1;
    iq = indwk + (*n << 1) + 2;
    mpos = 0;
    tbeg = 1;
    i__1 = nsplit;
    for (jblk = 1; jblk <= i__1; ++jblk) {
	tend = iwork[indisp + jblk - 1];
	bdim = tend - tbeg + 1;
	npos = bdim / 2;
	kb = npos;
/*        diagonals d~(1..KB) at local odd positions 1,3,... */
	i__2 = kb;
	for (l = 1; l <= i__2; ++l) {
	    work[idbl + l - 1] = (d__1 = work[inde2 + tbeg - 2 + (l << 1) - 1]
		    , abs(d__1));
/* L30: */
	}
	if (bdim % 2 == 0) {
/*           even block: square KB-by-KB, KB-1 off-diagonals */
	    i__2 = kb - 1;
	    for (l = 1; l <= i__2; ++l) {
		work[iebl + l - 1] = (d__1 = work[inde2 + tbeg - 2 + (l << 1)]
			, abs(d__1));
/* L32: */
	    }
	    dlasq1_(&kb, &work[idbl], &work[iebl], &work[iq], &iinfo);
	} else {
/*           odd block: append a zero diagonal -> (KB+1)-square, KB offdiag */
	    work[idbl + kb] = 0.;
	    i__2 = kb;
	    for (l = 1; l <= i__2; ++l) {
		work[iebl + l - 1] = (d__1 = work[inde2 + tbeg - 2 + (l << 1)]
			, abs(d__1));
/* L34: */
	    }
	    i__2 = kb + 1;
	    dlasq1_(&i__2, &work[idbl], &work[iebl], &work[iq], &iinfo);
	}
	if (iinfo != 0) {
	    *info = 1;
	    return;
	}
/*        DLASQ1 returns the singular values in DECREASING order in IDBL; */
/*        the appended structural zero (odd block) is the smallest entry, */
/*        so the NPOS positives are IDBL(0..NPOS-1).  Store ascending. */
	i__2 = npos;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s[mpos + i__] = work[idbl + npos - i__];
	    iwork[indbl + mpos + i__ - 1] = jblk;
	    iwork[indixw + mpos + i__ - 1] = bdim - npos + i__;
/* L40: */
	}
	mpos += npos;
	tbeg = tend + 1;
/* L60: */
    }

    nzsv = *n - mpos;

/*     Eigenvectors are not wanted: append the NZSV zero singular values, */
/*     sort ascending, rescale, and return. */

    if (! wantz) {
	i__1 = nzsv;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s[mpos + i__] = 0.;
/* L65: */
	}
	dbsort_(n, &s[1], &u[u_offset], ldu, &vt[vt_offset], ldvt, &c_false);
	if (scl != 1.) {
	    d__1 = 1. / scl;
	    dscal_(n, &d__1, &s[1], &c__1);
	}
	*m = *n;
	work[1] = (doublereal) lwmin;
	iwork[1] = liwmin;
	return;
    }

/*     NOTE: this must run BEFORE DLARRV_TGK, which overwrites the */
/*     off-diagonal array WORK(INDE2) in place with child representations. */
/*     Zero singular values.  Each odd-dimension TGK block has one */
/*     structural zero eigenvalue whose null vector lives on the ODD local */
/*     positions of the block (the even local entries are zero): */
/*        z(1)=1,  z(2l+1) = -beta(2l-1)/beta(2l) * z(2l-1). */
/*     De-interleaving by global parity, the vector is one-sided -- it is */
/*     a left singular vector (u) if the block starts at a u-row, else a */
/*     right singular vector (v) -- supplying one half of a zero singular */
/*     triplet (the other half comes from the paired block).  Because the */
/*     singular value is zero, any orthonormal completion is valid, so we */
/*     simply collect the u's into the trailing U-columns and the v's into */
/*     the trailing VT-rows.  Build the 2N null vector in WORK(INDWK) (the */
/*     eigenvector array, not yet in use). */

    kzu = 0;
    kzv = 0;
    tbeg = 1;
    i__1 = nsplit;
    for (jblk = 1; jblk <= i__1; ++jblk) {
	tend = iwork[indisp + jblk - 1];
	bdim = tend - tbeg + 1;
	if (bdim % 2 == 1) {
/*           build the local null vector spread over global rows TBEG..TEND */
	    i__2 = *n << 1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		work[indwk + i__ - 1] = 0.;
/* L110: */
	    }
	    work[indwk + tbeg - 1] = 1.;
	    znrm = 1.;
	    zg = 1.;
	    i__2 = (bdim - 1) / 2;
	    for (l = 1; l <= i__2; ++l) {
		ze = -(work[inde2 + tbeg - 2 + (l << 1) - 1] / work[inde2 + 
			tbeg - 2 + (l << 1)]) * zg;
		work[indwk + tbeg - 1 + (l << 1)] = ze;
		znrm += ze * ze;
		zg = ze;
/* L120: */
	    }
	    znrm = sqrt(znrm);
/*           The (one-sided) null vector lives on global-ODD rows when the */
/*           block starts at an odd row (MOD(TBEG,2)=1), else on global- */
/*           EVEN rows.  For LOWER B the odd half is u and the even half is */
/*           v; for UPPER the two roles swap.  Route the nonzero half to a */
/*           u-column or a v-row accordingly, normalizing by its own norm. */
	    if (tbeg % 2 == 1) {
		if (upper) {
		    ++kzv;
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			vt[mpos + kzv + i__ * vt_dim1] = work[indwk + (i__ << 
				1) - 2] / znrm;
/* L140: */
		    }
		} else {
		    ++kzu;
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			u[i__ + (mpos + kzu) * u_dim1] = work[indwk + (i__ << 
				1) - 2] / znrm;
/* L142: */
		    }
		}
	    } else {
		if (upper) {
		    ++kzu;
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			u[i__ + (mpos + kzu) * u_dim1] = work[indwk + (i__ << 
				1) - 1] / znrm;
/* L144: */
		    }
		} else {
		    ++kzv;
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			vt[mpos + kzv + i__ * vt_dim1] = work[indwk + (i__ << 
				1) - 1] / znrm;
/* L146: */
		    }
		}
	    }
	}
	tbeg = tend + 1;
/* L160: */
    }
/*     Compute the eigenvectors of the MPOS positive eigenvalues by MR^3 */
/*     rooted at the TGK matrix.  The 2N-by-MPOS eigenvector array is held */
/*     in WORK(INDWK) with leading dimension 2N. */

    if (mpos > 0) {
/* Computing MAX */
	d__1 = (doublereal) (*n << 1) * eps;
	tol = f2cmax(d__1,0.);
	i__1 = *n << 1;
	i__2 = *n << 1;
	dlarrv_tgk__(&i__1, &work[indtgk], &work[inde2], &iwork[indisp], &
		mpos, &s[1], &iwork[indbl], &iwork[indixw], &work[indgsc], &
		tol, &work[indwk], &i__2, &iwork[indisup], &work[indwk + (*n 
		<< 1) * mpos], &iwork[indiwk], &iinfo);
	if (iinfo != 0) {
	    *info = 2;
	    return;
	}

/*        De-interleave each 2N-eigenvector into its left/right singular */
/*        vectors by GLOBAL position parity (this stays correct across */
/*        split blocks).  beta_{2k-1}=d_k, beta_{2k}=e_k, so global-odd */
/*        rows carry one half and global-even rows the other.  For LOWER B */
/*        the odd half is u and the even half is v; for UPPER (B^T is the */
/*        lower bidiagonal with the same d,e) the two roles swap. */

	i__1 = mpos;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (upper) {
		    vt[k + i__ * vt_dim1] = rt2 * work[indwk + (k - 1 << 1) * 
			    *n + ((i__ << 1) - 2)];
		    u[i__ + k * u_dim1] = rt2 * work[indwk + (k - 1 << 1) * *
			    n + ((i__ << 1) - 1)];
		} else {
		    u[i__ + k * u_dim1] = rt2 * work[indwk + (k - 1 << 1) * *
			    n + ((i__ << 1) - 2)];
		    vt[k + i__ * vt_dim1] = rt2 * work[indwk + (k - 1 << 1) * 
			    *n + ((i__ << 1) - 1)];
		}
/* L90: */
	    }
/* L100: */
	}
    }


/*     Set the zero singular values and (defensively) any unmatched */
/*     trailing vectors, then sort everything ascending and rescale. */

    i__1 = nzsv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s[mpos + i__] = 0.;
/* L170: */
    }
    *m = *n;
    dbsort_(n, &s[1], &u[u_offset], ldu, &vt[vt_offset], ldvt, &c_true);
    if (scl != 1.) {
	d__1 = 1. / scl;
	dscal_(n, &d__1, &s[1], &c__1);
    }

    work[1] = (doublereal) lwmin;
    iwork[1] = liwmin;
    return;

/*     End of DBDSVDMR3 */

} /* dbdsvdmr3_work__ */


/*     ================================================================= */

/* Subroutine */ void dbsort_(integer *n, doublereal *s, doublereal *u, 
	integer *ldu, doublereal *vt, integer *ldvt, logical *wantz)
{
    /* System generated locals */
    integer u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;

    /* Local variables */
    integer i__, j, k;
    doublereal p;
    extern /* Subroutine */ void dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);


/*     Sort S(1:N) into ascending order by straight selection, applying */
/*     the same permutation to the columns of U and the rows of VT when */
/*     WANTZ is .TRUE. */

    /* Parameter adjustments */
    --s;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;

    /* Function Body */
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = i__;
	p = s[i__];
	i__2 = *n;
	for (j = i__ + 1; j <= i__2; ++j) {
	    if (s[j] < p) {
		k = j;
		p = s[j];
	    }
/* L10: */
	}
	if (k != i__) {
	    s[k] = s[i__];
	    s[i__] = p;
	    if (*wantz) {
		dswap_(n, &u[i__ * u_dim1 + 1], &c__1, &u[k * u_dim1 + 1], &
			c__1);
		dswap_(n, &vt[i__ + vt_dim1], ldvt, &vt[k + vt_dim1], ldvt);
	    }
	}
/* L20: */
    }
    return;
} /* dbsort_ */

