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

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;

/* Subroutine */ void cheevx_2stage_(char *jobz, char *range, char *uplo, 
	integer *n, complex *a, integer *lda, real *vl, real *vu, integer *il,
	 integer *iu, real *abstol, integer *m, real *w, complex *z__, 
	integer *ldz, complex *work, integer *lwork, real *rwork, integer *
	iwork, integer *ifail, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    integer iscale, i__, j, ib, kd, jj;
    real eps, vll, vuu, tmp1;
    integer indd, inde;
    extern integer ilaenv2stage_(integer *, char *, char *, integer *, 
	    integer *, integer *, integer *);
    real anrm, rmin, rmax;
    logical test;
    integer itmp1;
    extern /* Subroutine */ void chetrd_2stage_(char *, char *, integer *, 
	    complex *, integer *, real *, real *, complex *, complex *, 
	    integer *, complex *, integer *, integer *);
    integer indee;
    real sigma;
    extern logical lsame_(char *, char *);
    integer iinfo;
    extern /* Subroutine */ void sscal_(integer *, real *, real *, integer *);
    char order[1];
    integer lhtrd;
    extern /* Subroutine */ void cswap_(integer *, complex *, integer *, 
	    complex *, integer *);
    integer lwmin;
    logical lower;
    integer lwtrd;
    extern /* Subroutine */ void scopy_(integer *, real *, integer *, real *, 
	    integer *);
    logical wantz;
    extern real clanhe_(char *, char *, integer *, complex *, integer *, real 
	    *);
    logical alleig, indeig;
    integer indibl;
    logical valeig;
    extern real slamch_(char *);
    extern /* Subroutine */ void csscal_(integer *, real *, complex *, integer 
	    *), clacpy_(char *, integer *, integer *, complex *, integer *, 
	    complex *, integer *);
    real safmin;
    extern /* Subroutine */ void xerbla_(char *, integer *, ftnlen);
    real abstll, bignum;
    integer indiwk, indisp, indtau;
    extern /* Subroutine */ void cstein_(integer *, real *, real *, integer *, 
	    real *, integer *, integer *, complex *, integer *, real *, 
	    integer *, integer *, integer *);
    integer indrwk, indwrk;
    extern /* Subroutine */ void csteqr_(char *, integer *, real *, real *, 
	    complex *, integer *, real *, integer *), cungtr_(char *, 
	    integer *, complex *, integer *, complex *, complex *, integer *, 
	    integer *), ssterf_(integer *, real *, real *, integer *),
	     cunmtr_(char *, char *, char *, integer *, integer *, complex *, 
	    integer *, complex *, complex *, integer *, complex *, integer *, 
	    integer *);
    integer nsplit, llwork;
    real smlnum;
    extern /* Subroutine */ void sstebz_(char *, char *, integer *, real *, 
	    real *, integer *, integer *, real *, real *, real *, integer *, 
	    integer *, real *, integer *, integer *, real *, integer *, 
	    integer *);
    logical lquery;
    extern real sroundup_lwork__(integer *);
    integer indhous;



/*  -- LAPACK driver routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */


/*  ===================================================================== */


/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --rwork;
    --iwork;
    --ifail;

    /* Function Body */
    lower = lsame_(uplo, "L");
    wantz = lsame_(jobz, "V");
    alleig = lsame_(range, "A");
    valeig = lsame_(range, "V");
    indeig = lsame_(range, "I");
    lquery = *lwork == -1;

    *info = 0;
    if (! lsame_(jobz, "N")) {
	*info = -1;
    } else if (! (alleig || valeig || indeig)) {
	*info = -2;
    } else if (! (lower || lsame_(uplo, "U"))) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*lda < f2cmax(1,*n)) {
	*info = -6;
    } else {
	if (valeig) {
	    if (*n > 0 && *vu <= *vl) {
		*info = -8;
	    }
	} else if (indeig) {
	    if (*il < 1 || *il > f2cmax(1,*n)) {
		*info = -9;
	    } else if (*iu < f2cmin(*n,*il) || *iu > *n) {
		*info = -10;
	    }
	}
    }
    if (*info == 0) {
	if (*ldz < 1 || (wantz && *ldz < *n)) {
	    *info = -15;
	}
    }

    if (*info == 0) {
	if (*n <= 1) {
	    lwmin = 1;
	    r__1 = sroundup_lwork__(&lwmin);
	    work[1].r = r__1, work[1].i = 0.f;
	} else {
	    kd = ilaenv2stage_(&c__1, "CHETRD_2STAGE", jobz, n, &c_n1, &c_n1, 
		    &c_n1);
	    ib = ilaenv2stage_(&c__2, "CHETRD_2STAGE", jobz, n, &kd, &c_n1, &
		    c_n1);
	    lhtrd = ilaenv2stage_(&c__3, "CHETRD_2STAGE", jobz, n, &kd, &ib, &
		    c_n1);
	    lwtrd = ilaenv2stage_(&c__4, "CHETRD_2STAGE", jobz, n, &kd, &ib, &
		    c_n1);
	    lwmin = *n + lhtrd + lwtrd;
	    r__1 = sroundup_lwork__(&lwmin);
	    work[1].r = r__1, work[1].i = 0.f;
	}

	if (*lwork < lwmin && ! lquery) {
	    *info = -17;
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CHEEVX_2STAGE", &i__1, (ftnlen)13);
	return;
    } else if (lquery) {
	return;
    }

/*     Quick return if possible */

    *m = 0;
    if (*n == 0) {
	return;
    }

    if (*n == 1) {
	if (alleig || indeig) {
	    *m = 1;
	    i__1 = a_dim1 + 1;
	    w[1] = a[i__1].r;
	} else if (valeig) {
	    i__1 = a_dim1 + 1;
	    i__2 = a_dim1 + 1;
	    if (*vl < a[i__1].r && *vu >= a[i__2].r) {
		*m = 1;
		i__1 = a_dim1 + 1;
		w[1] = a[i__1].r;
	    }
	}
	if (wantz) {
	    i__1 = z_dim1 + 1;
	    z__[i__1].r = 1.f, z__[i__1].i = 0.f;
	}
	return;
    }

/*     Get machine constants. */

    safmin = slamch_("Safe minimum");
    eps = slamch_("Precision");
    smlnum = safmin / eps;
    bignum = 1.f / smlnum;
    rmin = sqrt(smlnum);
/* Computing MIN */
    r__1 = sqrt(bignum), r__2 = 1.f / sqrt(sqrt(safmin));
    rmax = f2cmin(r__1,r__2);

/*     Scale matrix to allowable range, if necessary. */

    iscale = 0;
    abstll = *abstol;
    if (valeig) {
	vll = *vl;
	vuu = *vu;
    }
    anrm = clanhe_("M", uplo, n, &a[a_offset], lda, &rwork[1]);
    if (anrm > 0.f && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
/*        An infinite norm gives SIGMA = 0; leave such a matrix */
/*        unscaled so that VL and VU remain a valid interval. */
	sigma = rmax / anrm;
	if (sigma > 0.f) {
	    iscale = 1;
	}
    }
    if (iscale == 1) {
	if (lower) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n - j + 1;
		csscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
/* L10: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		csscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
/* L20: */
	    }
	}
	if (*abstol > 0.f) {
	    abstll = *abstol * sigma;
	}
	if (valeig) {
	    vll = *vl * sigma;
	    vuu = *vu * sigma;
	}
    }

/*     Call CHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form. */

    indd = 1;
    inde = indd + *n;
    indrwk = inde + *n;
    indtau = 1;
    indhous = indtau + *n;
    indwrk = indhous + lhtrd;
    llwork = *lwork - indwrk + 1;

    chetrd_2stage_(jobz, uplo, n, &a[a_offset], lda, &rwork[indd], &rwork[
	    inde], &work[indtau], &work[indhous], &lhtrd, &work[indwrk], &
	    llwork, &iinfo);

/*     If all eigenvalues are desired and ABSTOL is less than or equal to */
/*     zero, then call SSTERF or CUNGTR and CSTEQR.  If this fails for */
/*     some eigenvalue, then try SSTEBZ. */

    test = FALSE_;
    if (indeig) {
	if (*il == 1 && *iu == *n) {
	    test = TRUE_;
	}
    }
    if ((alleig || test) && *abstol <= 0.f) {
	scopy_(n, &rwork[indd], &c__1, &w[1], &c__1);
	indee = indrwk + (*n << 1);
	if (! wantz) {
	    i__1 = *n - 1;
	    scopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
	    ssterf_(n, &w[1], &rwork[indee], info);
	} else {
	    clacpy_("A", n, n, &a[a_offset], lda, &z__[z_offset], ldz);
	    cungtr_(uplo, n, &z__[z_offset], ldz, &work[indtau], &work[indwrk]
		    , &llwork, &iinfo);
	    i__1 = *n - 1;
	    scopy_(&i__1, &rwork[inde], &c__1, &rwork[indee], &c__1);
	    csteqr_(jobz, n, &w[1], &rwork[indee], &z__[z_offset], ldz, &
		    rwork[indrwk], info);
	    if (*info == 0) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    ifail[i__] = 0;
/* L30: */
		}
	    }
	}
	if (*info == 0) {
	    *m = *n;
	    goto L40;
	}
	*info = 0;
    }

/*     Otherwise, call SSTEBZ and, if eigenvectors are desired, CSTEIN. */

    if (wantz) {
	*(unsigned char *)order = 'B';
    } else {
	*(unsigned char *)order = 'E';
    }
    indibl = 1;
    indisp = indibl + *n;
    indiwk = indisp + *n;
    sstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &rwork[indd], &
	    rwork[inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &
	    rwork[indrwk], &iwork[indiwk], &iinfo);
    if (iinfo != 0) {
	*info = *n + iinfo;
	if (iinfo != 1) {
	    goto L40;
	}
    }

    if (wantz) {
	cstein_(n, &rwork[indd], &rwork[inde], m, &w[1], &iwork[indibl], &
		iwork[indisp], &z__[z_offset], ldz, &rwork[indrwk], &iwork[
		indiwk], &ifail[1], &iinfo);
	if (iinfo != 0 && *info == 0) {
	    *info = iinfo;
	}

/*        Apply unitary matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by CSTEIN. */

	cunmtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau], &z__[
		z_offset], ldz, &work[indwrk], &llwork, &iinfo);
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

L40:
    if (iscale == 1) {
	r__1 = 1.f / sigma;
	sscal_(m, &r__1, &w[1], &c__1);
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors. */

    if (wantz) {
	i__1 = *m - 1;
	for (j = 1; j <= i__1; ++j) {
	    i__ = 0;
	    tmp1 = w[j];
	    i__2 = *m;
	    for (jj = j + 1; jj <= i__2; ++jj) {
		if (w[jj] < tmp1) {
		    i__ = jj;
		    tmp1 = w[jj];
		}
/* L50: */
	    }

	    if (i__ != 0) {
		itmp1 = iwork[indibl + i__ - 1];
		w[i__] = w[j];
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
		w[j] = tmp1;
		iwork[indibl + j - 1] = itmp1;
		cswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
		if (*info != 0) {
		    itmp1 = ifail[i__];
		    ifail[i__] = ifail[j];
		    ifail[j] = itmp1;
		}
	    }
/* L60: */
	}
    }

/*     Set WORK(1) to optimal complex workspace size. */

    r__1 = sroundup_lwork__(&lwmin);
    work[1].r = r__1, work[1].i = 0.f;

    return;

/*     End of CHEEVX_2STAGE */

} /* cheevx_2stage__ */

