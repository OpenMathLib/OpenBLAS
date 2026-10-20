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

static real c_b5 = 0.f;
static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b SLARRV_TGK MRRR representation-tree traversal computing eigenvectors of the Golub-Kahan matrix 
(Golub-Kahan variant of the LAPACK 3.0 SLARRV). */
/* > */
/* > \ingroup bdsvdmr3 */
/* > */
/* > \par Purpose: */
/* > ============= */
/* > */
/* > \verbatim */
/* > MRRR representation-tree traversal computing eigenvectors of the Golub-Kahan matrix (Golub-Kahan variant 
of the LAPACK 3.0 SLARRV). */
/* > All routines in this file are auxiliary to SBDSVDMR3 and derive */
/* > from the LAPACK 3.0 (1999) MRRR kernels SLARRV, SLAR1V, SLARRB and */
/* > SLARRF rather than from the current SSTEMR kernels. */
/* > \endverbatim */

/* Subroutine */ void slarrv_tgk__(integer *n, real *d__, real *l, integer *
	isplit, integer *m, real *w, integer *iblock, integer *indexw, real *
	gersch, real *tol, real *z__, integer *ldz, integer *isuppz, real *
	work, integer *iwork, integer *info)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    real r__1, r__2;

    /* Local variables */
    integer i__, j, k, p, q, k2, im, in;
    extern /* Subroutine */ void slar1v_tgk__(char *, integer *, integer *, 
	    integer *, real *, real *, real *, real *, real *, real *, real *,
	     real *, real *, real *, integer *, integer *, real *), 
	    slarrb_tgk__(integer *, real *, real *, real *, real *, integer *,
	     integer *, real *, real *, integer *, real *, real *, real *, 
	    real *, integer *, integer *), slarrf_tgk__(char *, integer *, 
	    real *, real *, real *, real *, integer *, integer *, real *, 
	    real *, real *, real *, real *, integer *, integer *);
    real gap;
    char rep[1];
    real eps, tmp;
    integer zto;
    real ztz, gap2, tmp2, ztz2;
    integer iend, jblk, wend, iter, temp[1];
    extern real sdot_(integer *, real *, integer *, real *, integer *);
    integer ktot;
    logical done1, done2;
    integer iter2, itmp1, itmp2;
    extern real snrm2_(integer *, real *, integer *);
    integer ktot2, iseed[4], indld;
    real sigma;
    integer ndone, iinfo, iindr;
    real resid;
    extern /* Subroutine */ void sscal_(integer *, real *, real *, integer *);
    logical nomgs;
    integer nclus, zfrom;
    extern /* Subroutine */ void scopy_(integer *, real *, integer *, real *, 
	    integer *), saxpy_(integer *, real *, real *, integer *, real *, 
	    integer *);
    integer iindc1, iindc2;
    real resid2, lambda;
    integer ibegin, indgap, indlld;
    extern real slamch_(char *);
    real mingma;
    integer oldien, oldncl, wbegin;
    real relgap;
    integer oldcls, ndepth, inderr, iindwk, newcls, oldfst;
    real minrgp;
    integer indwrk, oldlst;
    real lambda2, reltol;
    integer newfrs, newftt, parity;
    real mgstol, nrminv, rqcorr;
    integer newlst;
    extern /* Subroutine */ void slaset_(char *, integer *, integer *, real *, 
	    real *, real *, integer *);
    integer newsiz;
    extern /* Subroutine */ void sstein_(integer *, real *, real *, integer *, 
	    real *, integer *, integer *, real *, integer *, real *, integer *
	    , integer *, integer *);
    real mingma2, rqcorr2, nrminv2;
    logical ktot2ok;


/*  -- New auxiliary routine for bidiagonal SVD via TGK-rooted MR^3 -- */
/*     Based on LAPACK SLARRV (Dhillon & Marques, Nov 11 2003). */
/*     Derived from the LAPACK 3.0 SLARRV.  This is a separate copy with one */
/*     structural change: the ROOT of the representation tree (tree depth */
/*     NDEPTH = 0) is a symmetric tridiagonal -- the Tridiagonal Golub- */
/*     Kahan (TGK) matrix -- rather than an L D L^T factorization.  At the */
/*     root the eigenvectors / child RRRs are obtained with the TGK-input */
/*     paths of SLAR1V_TGK / SLARRF_TGK (REP = 'T'); every node at depth */
/*     >= 1 is an ordinary L D L^T child and uses REP = 'L'.  SLARRB */
/*     retains that algorithm with an added bounded-progress guard; */
/*     its positive INFO is propagated to the driver for safe fallback. */

/*  Usage / calling convention (set up by the driver SBDSVDMR3) */
/*  ---------------------------------------------------------- */
/*  On entry, for each split block the arrays D and L hold the TGK */
/*  representation of that block: D = 0 (the TGK diagonal) and L = the */
/*  TGK off-diagonals beta_i.  The block has TGK order IN = 2*nb where nb */
/*  is the bidiagonal block size; W holds the nb positive eigenvalues */
/*  (= singular values) of the block in ascending order, and */
/*      INDEXW(i) = nb + (position of W(i) among the block's positives) */
/*  because, after shifting by ~+sigma, the positive cluster occupies the */
/*  upper half of the block's 2*nb-point spectrum (all shifted negatives */
/*  lie below all shifted positives).  GERSCH holds the 2N TGK Gerschgorin */
/*  intervals from SBDTGK.  Otherwise the interface matches SLARRV. */

/*  ===================================================================== */

/*     Fixed seed for the random relative perturbations that SLARRF_TGK */
/*     applies to each child RRR (Sec. 5.4).  Threaded (in/out) through all */
/*     SLARRF_TGK calls; fixed initial value makes the result reproducible. */
    /* Parameter adjustments */
    --d__;
    --l;
    --isplit;
    --w;
    --iblock;
    --indexw;
    --gersch;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --work;
    --iwork;

    /* Function Body */
    iseed[0] = 1;
    iseed[1] = 13;
    iseed[2] = 257;
    iseed[3] = 2047;

    inderr = *n;
    indld = *n << 1;
    indlld = *n * 3;
    indgap = *n << 2;
    indwrk = *n * 5 + 1;

    iindr = *n;
    iindc1 = *n << 1;
    iindc2 = *n * 3;
    iindwk = (*n << 2) + 1;

    eps = slamch_("Precision");

    i__1 = *n << 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iwork[i__] = 0;
/* L10: */
    }
/*     Zero the N-by-M eigenvector array (here N is the TGK order 2*nb and */
/*     M is the number of selected eigenvalues; Z has only M columns, unlike */
/*     the original SLARRV where the matrix order and column count coincide). */
    slaset_("Full", n, m, &c_b5, &c_b5, &z__[z_offset], ldz);
    mgstol = eps * 100.f;

    ibegin = 1;
    wbegin = 1;
    i__1 = iblock[*m];
    for (jblk = 1; jblk <= i__1; ++jblk) {
	iend = isplit[jblk];

/*        Find the eigenvectors of the submatrix indexed IBEGIN */
/*        through IEND. */

	wend = wbegin - 1;
L171:
	if (wend < *m) {
	    if (iblock[wend + 1] == jblk) {
		++wend;
		goto L171;
	    }
	}
	if (wend < wbegin) {
	    ibegin = iend + 1;
	    goto L170;
	}

	if (ibegin == iend) {
	    z__[ibegin + wbegin * z_dim1] = 1.f;
	    isuppz[(wbegin << 1) - 1] = ibegin;
	    isuppz[wbegin * 2] = ibegin;
	    ibegin = iend + 1;
	    wbegin = wend + 1;
	    goto L170;
	}
	oldien = ibegin - 1;
	in = iend - oldien;
/* Computing MIN */
	r__1 = .003f, r__2 = 1.f / (real) in;
	reltol = f2cmin(r__1,r__2);
	im = wend - wbegin + 1;
	scopy_(&im, &w[wbegin], &c__1, &work[1], &c__1);
	i__2 = im - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[inderr + i__] = eps * (r__1 = work[i__], abs(r__1));
	    work[indgap + i__] = work[i__ + 1] - work[i__];
/* L30: */
	}
	work[inderr + im] = eps * (r__1 = work[im], abs(r__1));
/* Computing MAX */
	r__2 = (r__1 = work[im], abs(r__1));
	work[indgap + im] = f2cmax(r__2,eps);
	ndone = 0;

	ndepth = 0;
	parity = 1;
	nclus = 1;
	iwork[iindc1 + 1] = 1;
	iwork[iindc1 + 2] = im;

/*        While( NDONE.LT.IM ) do */

L40:
	if (ndone < im && ndepth <= *n << 2) {

/*           REP selects the representation of the current tree level: */
/*           'T' (TGK tridiagonal) at the root, 'L' (L D L^T) deeper down. */

	    if (ndepth == 0) {
		*(unsigned char *)rep = 'T';
	    } else {
		*(unsigned char *)rep = 'L';
	    }
	    oldncl = nclus;
	    nclus = 0;
	    parity = 1 - parity;
	    if (parity == 0) {
		oldcls = iindc1;
		newcls = iindc2;
	    } else {
		oldcls = iindc2;
		newcls = iindc1;
	    }
	    i__2 = oldncl;
	    for (i__ = 1; i__ <= i__2; ++i__) {

/*              If NDEPTH > 1, retrieve the relatively robust */
/*              representation (RRR) and perform limited bisection */
/*              (if necessary) to get approximate eigenvalues. */

		j = oldcls + (i__ << 1);
		oldfst = iwork[j - 1];
		oldlst = iwork[j];
		if (ndepth > 0) {
		    j = wbegin + oldfst - 1;
		    scopy_(&in, &z__[ibegin + j * z_dim1], &c__1, &d__[ibegin]
			    , &c__1);
		    i__3 = in - 1;
		    scopy_(&i__3, &z__[ibegin + (j + 1) * z_dim1], &c__1, &l[
			    ibegin], &c__1);
		    slaset_("Full", &in, &c__2, &c_b5, &c_b5, &z__[ibegin + j 
			    * z_dim1], ldz);
		}

/*              Compute LD and LLD of the current representation.  For the */
/*              TGK root (REP='T') the diagonal is zero so these are zero */
/*              and unused (SLARRB is skipped, and the TGK paths of */
/*              SLARRF_TGK/SLAR1V_TGK ignore LD,LLD); they are needed for */
/*              every L D L^T child (REP='L'). */

		k = ibegin;
		i__3 = in - 1;
		for (j = 1; j <= i__3; ++j) {
		    tmp = d__[k] * l[k];
		    work[indld + j] = tmp;
		    work[indlld + j] = tmp * l[k];
		    ++k;
/* L50: */
		}
		if (ndepth > 0) {

/*  The way OLDFST and OLDLST are used in the following is a problem */
/*  when eigenvalues in a range are desired! */

		    p = indexw[wbegin - 1 + oldfst];
		    q = indexw[wbegin - 1 + oldlst];
		    r__1 = eps * 4.f;
		    i__3 = p - oldfst;
		    slarrb_tgk__(&in, &d__[ibegin], &l[ibegin], &work[indld + 
			    1], &work[indlld + 1], &p, &q, &reltol, &r__1, &
			    i__3, &work[1], &work[indgap + 1], &work[inderr + 
			    1], &work[indwrk + in], &iwork[iindwk], &iinfo);
		    if (iinfo != 0) {
			*info = 3;
			return;
		    }
		}

/*              Classify eigenvalues of the current representation (RRR) */
/*              as (i) isolated, or (ii) clustered. */

		newfrs = oldfst;
		i__3 = oldlst;
		for (j = oldfst; j <= i__3; ++j) {
		    if (j == oldlst || work[indgap + j] >= reltol * (r__1 = 
			    work[j], abs(r__1))) {
			newlst = j;
		    } else {

/*                    continue (to the next loop) */

			relgap = work[indgap + j] / (r__1 = work[j], abs(r__1)
				);
			if (j == newfrs) {
			    minrgp = relgap;
			} else {
			    minrgp = f2cmin(minrgp,relgap);
			}
			goto L140;
		    }
		    newsiz = newlst - newfrs + 1;
		    newftt = wbegin + newfrs - 1;
		    nomgs = newsiz == 1 || newsiz > 1 || minrgp < mgstol;
		    if (newsiz > 1 && nomgs) {

/*                    Find a new L D L^T representation if the cluster */
/*                    size is larger than MGSSIZ or the minimum */
/*                    relative gap within the cluster is too small. */
/*                    REP='T' at the root forms the child from the TGK */
/*                    matrix; REP='L' forms it from the parent L D L^T. */

			slarrf_tgk__(rep, &in, &d__[ibegin], &l[ibegin], &
				work[indld + 1], &work[indlld + 1], &newfrs, &
				newlst, &work[1], &sigma, &z__[ibegin + 
				newftt * z_dim1], &z__[ibegin + (newftt + 1) *
				 z_dim1], &work[indwrk], iseed, info);
			if (*info == 0) {
			    tmp = eps * abs(sigma);
			    i__4 = newlst;
			    for (k = newfrs; k <= i__4; ++k) {
				work[k] -= sigma;
/* Computing MAX */
				r__1 = work[indgap + k];
				work[indgap + k] = f2cmax(r__1,tmp);
				work[inderr + k] += tmp;
/* L51: */
			    }
			    ++nclus;
			    k = newcls + (nclus << 1);
			    iwork[k - 1] = newfrs;
			    iwork[k] = newlst;
			} else {
			    *info = 0;
			    if (minrgp >= mgstol) {
/*                          No relatively robust child representation */
/*                          could be found for this cluster; report it */
/*                          instead of leaving the cluster unprocessed */
/*                          (the caller falls back to SBDSDC). */
				*info = 3;
				return;
			    } else if (*(unsigned char *)rep == 'T') {
/*                          At the root the representation is the */
/*                          Golub-Kahan matrix by its entries, so the */
/*                          LDL^T-based inverse-iteration fallback does */
/*                          not apply. */
				*info = 2;
				return;
			    } else {

/*                          Call SSTEIN to process this tight cluster. */
/*                          This happens only if MINRGP <= MGSTOL */
/*                          and SLARRF_TGK returns INFO = 1. The latter */
/*                          means that a new RRR to "break" the */
/*                          cluster could not be found. */

				work[indwrk] = d__[ibegin];
				i__4 = in - 1;
				for (k = 1; k <= i__4; ++k) {
				    work[indwrk + k] = d__[ibegin + k] + work[
					    indlld + k];
/* L60: */
				}
				i__4 = newsiz;
				for (k = 1; k <= i__4; ++k) {
				    iwork[iindwk + k - 1] = 1;
/* L70: */
				}
				i__4 = newlst;
				for (k = newfrs; k <= i__4; ++k) {
				    isuppz[(oldien + k << 1) - 1] = 1;
				    isuppz[(oldien + k) * 2] = in;
/* L80: */
				}
				temp[0] = in;
				sstein_(&in, &work[indwrk], &work[indld + 1], 
					&newsiz, &work[newfrs], &iwork[iindwk]
					, temp, &z__[ibegin + newftt * z_dim1]
					, ldz, &work[indwrk + in], &iwork[
					iindwk + in], &iwork[iindwk + (in << 
					1)], &iinfo);
				if (iinfo != 0) {
				    *info = 2;
				    return;
				}
				ndone += newsiz;
			    }
			}
		    } else if (newsiz > 1) {
			ktot = newftt;
			i__4 = newlst;
			for (k = newfrs; k <= i__4; ++k) {
			    iter = 0;
L90:
			    lambda = work[k];

/*                       Given LAMBDA, compute the eigenvector.  REP='T' */
/*                       uses the TGK matrix as input at the root. */

			    slar1v_tgk__(rep, &in, &c__1, &in, &lambda, &d__[
				    ibegin], &l[ibegin], &work[indld + 1], &
				    work[indlld + 1], &w[wbegin + k - 1], &
				    gersch[(oldien << 1) + 1], &z__[ibegin + 
				    ktot * z_dim1], &ztz, &mingma, &iwork[
				    iindr + ktot], &isuppz[(ktot << 1) - 1], &
				    work[indwrk]);
			    tmp = 1.f / ztz;
			    nrminv = sqrt(tmp);
			    resid = abs(mingma) * nrminv;
			    rqcorr = mingma * tmp;
			    if (k == in) {
				gap = work[indgap + k - 1];
			    } else if (k == 1) {
				gap = work[indgap + k];
			    } else {
/* Computing MIN */
				r__1 = work[indgap + k - 1], r__2 = work[
					indgap + k];
				gap = f2cmin(r__1,r__2);
			    }
			    ++iter;
			    if (resid > *tol * gap && abs(rqcorr) > eps * 4.f 
				    * abs(lambda)) {
				work[k] = lambda + rqcorr;
				if (iter < 8) {
				    goto L90;
				}
			    }
			    iwork[ktot] = 1;
			    if (newsiz == 1) {
				++ndone;
			    }
			    zfrom = isuppz[(ktot << 1) - 1];
			    zto = isuppz[ktot * 2];
			    i__5 = zto - zfrom + 1;
			    sscal_(&i__5, &nrminv, &z__[ibegin + zfrom - 1 + 
				    ktot * z_dim1], &c__1);
			    ++ktot;
/* L100: */
			}
			if (newsiz > 1) {
			    itmp1 = isuppz[(newftt << 1) - 1];
			    itmp2 = isuppz[newftt * 2];
			    ktot = oldien + newlst;
			    i__4 = ktot;
			    for (p = newftt + 1; p <= i__4; ++p) {
				i__5 = p - 1;
				for (q = newftt; q <= i__5; ++q) {
				    tmp = -sdot_(&in, &z__[ibegin + p * 
					    z_dim1], &c__1, &z__[ibegin + q * 
					    z_dim1], &c__1);
				    saxpy_(&in, &tmp, &z__[ibegin + q * 
					    z_dim1], &c__1, &z__[ibegin + p * 
					    z_dim1], &c__1);
/* L110: */
				}
				tmp = 1.f / snrm2_(&in, &z__[ibegin + p * 
					z_dim1], &c__1);
				sscal_(&in, &tmp, &z__[ibegin + p * z_dim1], &
					c__1);
/* Computing MIN */
				i__5 = itmp1, i__6 = isuppz[(p << 1) - 1];
				itmp1 = f2cmin(i__5,i__6);
/* Computing MAX */
				i__5 = itmp2, i__6 = isuppz[p * 2];
				itmp2 = f2cmax(i__5,i__6);
/* L120: */
			    }
			    i__4 = ktot;
			    for (p = newftt; p <= i__4; ++p) {
				isuppz[(p << 1) - 1] = itmp1;
				isuppz[p * 2] = itmp2;
/* L130: */
			    }
			    ndone += newsiz;
			}
		    } else {
/*                    Defer isolated eigenvectors until this complete parent */
/*                    RRR has been classified.  Consecutive marked vectors */
/*                    can then run through two independent SLAR1V recurrences */
/*                    together without changing either lane's FP order. */
			iwork[newftt] = -1;
		    }
		    newfrs = j + 1;
L140:
		    ;
		}

/*              Process the isolated eigenvectors belonging to this parent */
/*              representation.  Pair adjacent lanes when possible; an odd */
/*              tail follows the original scalar path verbatim.  The second */
/*              4*IN scratch lane ends at WORK(13*N), exactly within the */
/*              private SLARRV_TGK workspace supplied by SBDSVDMR3. */

		k = oldfst;
L141:
		if (k <= oldlst) {
		    ktot = wbegin + k - 1;
		    if (iwork[ktot] != -1) {
			++k;
			goto L141;
		    }
		    k2 = k + 1;
		    if (k2 <= oldlst) {
			ktot2 = wbegin + k2 - 1;
		    } else {
			ktot2 = 0;
		    }
		    ktot2ok = FALSE_;
		    if (ktot2 > 0) {
			if (iwork[ktot2] == -1) {
			    ktot2ok = TRUE_;
			}
		    }
		    if (ktot2ok) {
			iter = 0;
			iter2 = 0;
			done1 = FALSE_;
			done2 = FALSE_;
L142:
			lambda = work[k];
			lambda2 = work[k2];
/*                    The paired two-lane variant SLAR1V2_TGK is kept in */
/*                    slar1v_tgk.f but not used: the two eigenvectors are */
/*                    computed by separate one-lane calls (reference code). */
			if (! done1) {
			    slar1v_tgk__(rep, &in, &c__1, &in, &lambda, &d__[
				    ibegin], &l[ibegin], &work[indld + 1], &
				    work[indlld + 1], &w[wbegin + k - 1], &
				    gersch[(oldien << 1) + 1], &z__[ibegin + 
				    ktot * z_dim1], &ztz, &mingma, &iwork[
				    iindr + ktot], &isuppz[(ktot << 1) - 1], &
				    work[indwrk]);
			}
			if (! done2) {
			    slar1v_tgk__(rep, &in, &c__1, &in, &lambda2, &d__[
				    ibegin], &l[ibegin], &work[indld + 1], &
				    work[indlld + 1], &w[wbegin + k2 - 1], &
				    gersch[(oldien << 1) + 1], &z__[ibegin + 
				    ktot2 * z_dim1], &ztz2, &mingma2, &iwork[
				    iindr + ktot2], &isuppz[(ktot2 << 1) - 1],
				     &work[indwrk + (in << 2)]);
			}
			if (! done1) {
			    tmp = 1.f / ztz;
			    nrminv = sqrt(tmp);
			    resid = abs(mingma) * nrminv;
			    rqcorr = mingma * tmp;
			    if (k == in) {
				gap = work[indgap + k - 1];
			    } else if (k == 1) {
				gap = work[indgap + k];
			    } else {
/* Computing MIN */
				r__1 = work[indgap + k - 1], r__2 = work[
					indgap + k];
				gap = f2cmin(r__1,r__2);
			    }
			    ++iter;
			    done1 = TRUE_;
			    if (resid > *tol * gap && abs(rqcorr) > eps * 4.f 
				    * abs(lambda)) {
				work[k] = lambda + rqcorr;
				if (iter < 8) {
				    done1 = FALSE_;
				}
			    }
			}
			if (! done2) {
			    tmp2 = 1.f / ztz2;
			    nrminv2 = sqrt(tmp2);
			    resid2 = abs(mingma2) * nrminv2;
			    rqcorr2 = mingma2 * tmp2;
			    if (k2 == in) {
				gap2 = work[indgap + k2 - 1];
			    } else if (k2 == 1) {
				gap2 = work[indgap + k2];
			    } else {
/* Computing MIN */
				r__1 = work[indgap + k2 - 1], r__2 = work[
					indgap + k2];
				gap2 = f2cmin(r__1,r__2);
			    }
			    ++iter2;
			    done2 = TRUE_;
			    if (resid2 > *tol * gap2 && abs(rqcorr2) > eps * 
				    4.f * abs(lambda2)) {
				work[k2] = lambda2 + rqcorr2;
				if (iter2 < 8) {
				    done2 = FALSE_;
				}
			    }
			}
			if (! done1 || ! done2) {
			    goto L142;
			}
			iwork[ktot] = 1;
			++ndone;
			zfrom = isuppz[(ktot << 1) - 1];
			zto = isuppz[ktot * 2];
			i__3 = zto - zfrom + 1;
			sscal_(&i__3, &nrminv, &z__[ibegin + zfrom - 1 + ktot 
				* z_dim1], &c__1);
			iwork[ktot2] = 1;
			++ndone;
			zfrom = isuppz[(ktot2 << 1) - 1];
			zto = isuppz[ktot2 * 2];
			i__3 = zto - zfrom + 1;
			sscal_(&i__3, &nrminv2, &z__[ibegin + zfrom - 1 + 
				ktot2 * z_dim1], &c__1);
			k += 2;
		    } else {
			iter = 0;
L143:
			lambda = work[k];
			slar1v_tgk__(rep, &in, &c__1, &in, &lambda, &d__[
				ibegin], &l[ibegin], &work[indld + 1], &work[
				indlld + 1], &w[wbegin + k - 1], &gersch[(
				oldien << 1) + 1], &z__[ibegin + ktot * 
				z_dim1], &ztz, &mingma, &iwork[iindr + ktot], 
				&isuppz[(ktot << 1) - 1], &work[indwrk]);
			tmp = 1.f / ztz;
			nrminv = sqrt(tmp);
			resid = abs(mingma) * nrminv;
			rqcorr = mingma * tmp;
			if (k == in) {
			    gap = work[indgap + k - 1];
			} else if (k == 1) {
			    gap = work[indgap + k];
			} else {
/* Computing MIN */
			    r__1 = work[indgap + k - 1], r__2 = work[indgap + 
				    k];
			    gap = f2cmin(r__1,r__2);
			}
			++iter;
			if (resid > *tol * gap && abs(rqcorr) > eps * 4.f * 
				abs(lambda)) {
			    work[k] = lambda + rqcorr;
			    if (iter < 8) {
				goto L143;
			    }
			}
			iwork[ktot] = 1;
			++ndone;
			zfrom = isuppz[(ktot << 1) - 1];
			zto = isuppz[ktot * 2];
			i__3 = zto - zfrom + 1;
			sscal_(&i__3, &nrminv, &z__[ibegin + zfrom - 1 + ktot 
				* z_dim1], &c__1);
			++k;
		    }
		    goto L141;
		}
/* L150: */
	    }
	    ++ndepth;
	    goto L40;
	}
	j = wbegin << 1;
	i__2 = wend;
	for (i__ = wbegin; i__ <= i__2; ++i__) {
	    isuppz[j - 1] += oldien;
	    isuppz[j] += oldien;
	    j += 2;
/* L160: */
	}
	ibegin = iend + 1;
	wbegin = wend + 1;
L170:
	;
    }

    return;

/*     End of SLARRV_TGK */

} /* slarrv_tgk__ */

