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

static integer c__10 = 10;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c_n1 = -1;

/* Subroutine */ void ssyevr_2stage_(char *jobz, char *range, char *uplo, 
	integer *n, real *a, integer *lda, real *vl, real *vu, integer *il, 
	integer *iu, real *abstol, integer *m, real *w, real *z__, integer *
	ldz, integer *isuppz, real *work, integer *lwork, integer *iwork, 
	integer *liwork, integer *info)
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
    integer inddd, indee;
    real sigma;
    extern logical lsame_(char *, char *);
    integer iinfo;
    extern /* Subroutine */ void sscal_(integer *, real *, real *, integer *);
    char order[1];
    integer indwk, lhtrd, lwmin;
    logical lower;
    integer lwtrd;
    extern /* Subroutine */ void scopy_(integer *, real *, integer *, real *, 
	    integer *), sswap_(integer *, real *, integer *, real *, integer *
	    ), ssytrd_2stage_(char *, char *, integer *, real *, integer *, 
	    real *, real *, real *, real *, integer *, real *, integer *, 
	    integer *);
    logical wantz, alleig, indeig;
    integer ieeeok, indibl, indifl;
    logical valeig;
    extern real slamch_(char *);
    real safmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ void xerbla_(char *, integer *, ftnlen);
    real abstll, bignum;
    integer indtau, indisp, indiwo, indwkn, liwmin;
    logical tryrac;
    extern /* Subroutine */ void sstein_(integer *, real *, real *, integer *, 
	    real *, integer *, integer *, real *, integer *, real *, integer *
	    , integer *, integer *), ssterf_(integer *, real *, real *, 
	    integer *);
    integer llwrkn, llwork, nsplit;
    real smlnum;
    extern real slansy_(char *, char *, integer *, real *, integer *, real *);
    extern /* Subroutine */ void sstebz_(char *, char *, integer *, real *, 
	    real *, integer *, integer *, real *, real *, real *, integer *, 
	    integer *, real *, integer *, integer *, real *, integer *, 
	    integer *), sstemr_(char *, char *, integer *, 
	    real *, real *, real *, real *, integer *, integer *, integer *, 
	    real *, real *, integer *, integer *, integer *, logical *, real *
	    , integer *, integer *, integer *, integer *);
    logical lquery;
    extern /* Subroutine */ void sormtr_(char *, char *, char *, integer *, 
	    integer *, real *, integer *, real *, real *, integer *, real *, 
	    integer *, integer *);
    extern real sroundup_lwork__(integer *);
    integer indhous;



/*  -- LAPACK driver routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */


/* ===================================================================== */


/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --work;
    --iwork;

    /* Function Body */
    ieeeok = ilaenv_(&c__10, "SSYEVR", "N", &c__1, &c__2, &c__3, &c__4, (
	    ftnlen)6, (ftnlen)1);

    lower = lsame_(uplo, "L");
    wantz = lsame_(jobz, "V");
    alleig = lsame_(range, "A");
    valeig = lsame_(range, "V");
    indeig = lsame_(range, "I");

    lquery = *lwork == -1 || *liwork == -1;

    kd = ilaenv2stage_(&c__1, "SSYTRD_2STAGE", jobz, n, &c_n1, &c_n1, &c_n1);
    ib = ilaenv2stage_(&c__2, "SSYTRD_2STAGE", jobz, n, &kd, &c_n1, &c_n1);
    lhtrd = ilaenv2stage_(&c__3, "SSYTRD_2STAGE", jobz, n, &kd, &ib, &c_n1);
    lwtrd = ilaenv2stage_(&c__4, "SSYTRD_2STAGE", jobz, n, &kd, &ib, &c_n1);

    if (*n <= 1) {
	lwmin = 1;
	liwmin = 1;
    } else {
/* Computing MAX */
	i__1 = *n * 26, i__2 = *n * 5 + lhtrd + lwtrd;
	lwmin = f2cmax(i__1,i__2);
	liwmin = *n * 10;
    }

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
	} else if (*lwork < lwmin && ! lquery) {
	    *info = -18;
	} else if (*liwork < liwmin && ! lquery) {
	    *info = -20;
	}
    }

    if (*info == 0) {
/*         NB = ILAENV( 1, 'SSYTRD', UPLO, N, -1, -1, -1 ) */
/*         NB = MAX( NB, ILAENV( 1, 'SORMTR', UPLO, N, -1, -1, -1 ) ) */
/*         LWKOPT = MAX( ( NB+1 )*N, LWMIN ) */
	work[1] = sroundup_lwork__(&lwmin);
	iwork[1] = liwmin;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SSYEVR_2STAGE", &i__1, (ftnlen)13);
	return;
    } else if (lquery) {
	return;
    }

/*     Quick return if possible */

    *m = 0;
    if (*n == 0) {
	work[1] = 1.f;
	return;
    }

    if (*n == 1) {
	work[1] = 1.f;
	if (alleig || indeig) {
	    *m = 1;
	    w[1] = a[a_dim1 + 1];
	} else {
	    if (*vl < a[a_dim1 + 1] && *vu >= a[a_dim1 + 1]) {
		*m = 1;
		w[1] = a[a_dim1 + 1];
	    }
	}
	if (wantz) {
	    z__[z_dim1 + 1] = 1.f;
	    isuppz[1] = 1;
	    isuppz[2] = 1;
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
    anrm = slansy_("M", uplo, n, &a[a_offset], lda, &work[1]);
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
		sscal_(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
/* L10: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		sscal_(&j, &sigma, &a[j * a_dim1 + 1], &c__1);
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
/*     Initialize indices into workspaces.  Note: The IWORK indices are */
/*     used only if SSTERF or SSTEMR fail. */
/*     WORK(INDTAU:INDTAU+N-1) stores the scalar factors of the */
/*     elementary reflectors used in SSYTRD. */
    indtau = 1;
/*     WORK(INDD:INDD+N-1) stores the tridiagonal's diagonal entries. */
    indd = indtau + *n;
/*     WORK(INDE:INDE+N-1) stores the off-diagonal entries of the */
/*     tridiagonal matrix from SSYTRD. */
    inde = indd + *n;
/*     WORK(INDDD:INDDD+N-1) is a copy of the diagonal entries over */
/*     -written by SSTEMR (the SSTERF path copies the diagonal to W). */
    inddd = inde + *n;
/*     WORK(INDEE:INDEE+N-1) is a copy of the off-diagonal entries over */
/*     -written while computing the eigenvalues in SSTERF and SSTEMR. */
    indee = inddd + *n;
/*     INDHOUS is the starting offset Householder storage of stage 2 */
    indhous = indee + *n;
/*     INDWK is the starting offset of the left-over workspace, and */
/*     LLWORK is the remaining workspace size. */
    indwk = indhous + lhtrd;
    llwork = *lwork - indwk + 1;
/*     IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in SSTEBZ and */
/*     stores the block indices of each of the M<=N eigenvalues. */
    indibl = 1;
/*     IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in SSTEBZ and */
/*     stores the starting and finishing indices of each block. */
    indisp = indibl + *n;
/*     IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors */
/*     that corresponding to eigenvectors that fail to converge in */
/*     SSTEIN.  This information is discarded; if any fail, the driver */
/*     returns INFO > 0. */
    indifl = indisp + *n;
/*     INDIWO is the offset of the remaining integer workspace. */
    indiwo = indifl + *n;

/*     Call SSYTRD_2STAGE to reduce symmetric matrix to tridiagonal form. */


    ssytrd_2stage_(jobz, uplo, n, &a[a_offset], lda, &work[indd], &work[inde]
	    , &work[indtau], &work[indhous], &lhtrd, &work[indwk], &llwork, &
	    iinfo);

/*     If all eigenvalues are desired */
/*     then call SSTERF or SSTEMR and SORMTR. */

    test = FALSE_;
    if (indeig) {
	if (*il == 1 && *iu == *n) {
	    test = TRUE_;
	}
    }
    if ((alleig || test) && ieeeok == 1) {
	if (! wantz) {
	    scopy_(n, &work[indd], &c__1, &w[1], &c__1);
	    i__1 = *n - 1;
	    scopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
	    ssterf_(n, &w[1], &work[indee], info);
	} else {
	    i__1 = *n - 1;
	    scopy_(&i__1, &work[inde], &c__1, &work[indee], &c__1);
	    scopy_(n, &work[indd], &c__1, &work[inddd], &c__1);

	    if (*abstol <= (real) (*n) * 2.f * eps) {
		tryrac = TRUE_;
	    } else {
		tryrac = FALSE_;
	    }
	    sstemr_(jobz, "A", n, &work[inddd], &work[indee], vl, vu, il, iu, 
		    m, &w[1], &z__[z_offset], ldz, n, &isuppz[1], &tryrac, &
		    work[indwk], lwork, &iwork[1], liwork, info);



/*        Apply orthogonal matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by SSTEMR. */

	    if (wantz && *info == 0) {
		indwkn = inde;
		llwrkn = *lwork - indwkn + 1;
		sormtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau]
			, &z__[z_offset], ldz, &work[indwkn], &llwrkn, &iinfo);
	    }
	}


	if (*info == 0) {
/*           Everything worked.  Skip SSTEBZ/SSTEIN.  IWORK(:) are */
/*           undefined. */
	    *m = *n;
	    goto L30;
	}
	*info = 0;
    }

/*     Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN. */
/*     Also call SSTEBZ and SSTEIN if SSTEMR fails. */

    if (wantz) {
	*(unsigned char *)order = 'B';
    } else {
	*(unsigned char *)order = 'E';
    }
    sstebz_(range, order, n, &vll, &vuu, il, iu, &abstll, &work[indd], &work[
	    inde], m, &nsplit, &w[1], &iwork[indibl], &iwork[indisp], &work[
	    indwk], &iwork[indiwo], info);

    if (wantz) {
	sstein_(n, &work[indd], &work[inde], m, &w[1], &iwork[indibl], &iwork[
		indisp], &z__[z_offset], ldz, &work[indwk], &iwork[indiwo], &
		iwork[indifl], info);

/*        Apply orthogonal matrix used in reduction to tridiagonal */
/*        form to eigenvectors returned by SSTEIN. */

	indwkn = inde;
	llwrkn = *lwork - indwkn + 1;
	sormtr_("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau], &z__[
		z_offset], ldz, &work[indwkn], &llwrkn, &iinfo);
    }

/*     If matrix was scaled, then rescale eigenvalues appropriately. */

/*  Jump here if SSTEMR/SSTEIN succeeded. */
L30:
    if (iscale == 1) {
	r__1 = 1.f / sigma;
	sscal_(m, &r__1, &w[1], &c__1);
    }

/*     If eigenvalues are not in order, then sort them, along with */
/*     eigenvectors.  Note: We do not sort the IFAIL portion of IWORK. */
/*     It may not be initialized (if SSTEMR/SSTEIN succeeded), and we do */
/*     not return this detailed information to the user. */

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
/* L40: */
	    }

	    if (i__ != 0) {
		w[i__] = w[j];
		w[j] = tmp1;
		sswap_(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
	    }
/* L50: */
	}
    }

/*     Set WORK(1) to optimal workspace size. */

    work[1] = sroundup_lwork__(&lwmin);
    iwork[1] = liwmin;

    return;

/*     End of SSYEVR_2STAGE */

} /* ssyevr_2stage__ */

