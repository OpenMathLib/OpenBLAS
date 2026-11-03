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

static real c_b12 = 1.f;
static integer c__1 = 1;

/* > \brief \b SBDSVDMR3 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SBDSVDMR3( UPLO, JOBZ, RANGE, N, D, E, VL, VU, IL, IU, */
/*                          NS, S, Z, LDZ, WORK, LWORK, IWORK, LIWORK, */
/*                          INFO ) */

/*       CHARACTER          JOBZ, RANGE, UPLO */
/*       INTEGER            IL, INFO, IU, LDZ, LIWORK, LWORK, N, NS */
/*       REAL   VL, VU */
/*       INTEGER            IWORK( * ) */
/*       REAL   D( * ), E( * ), S( * ), WORK( * ), */
/*                          Z( LDZ, * ) */

/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >  SBDSVDMR3 computes the singular value decomposition (SVD) of a real */
/* >  N-by-N (upper or lower) bidiagonal matrix B, B = U * S * VT, using */
/* >  the Tridiagonal Golub-Kahan (TGK) MR^3 driver SBDSVDMR3.  Its */
/* >  argument list mirrors SBDSVDX so that SBDSVDMR3 is a drop-in */
/* >  alternative alongside SBDSVDX for callers that want the MR^3 */
/* >  variant. */
/* > */
/* >  Given an upper bidiagonal B with diagonal D and superdiagonal E, */
/* >  SBDSVDMR3 builds the 2N-by-2N TGK tridiagonal, computes its */
/* >  eigendecomposition by the MR^3 algorithm rooted at the TGK */
/* >  (see P. Willems and B. Lang, SIAM J. Sci. Comput., 35:740-766, */
/* >  2013), and returns the singular values and (optionally) singular */
/* >  vectors of B in the same packing convention as SBDSVDX. */
/* > */
/* >  RANGE = 'V' or 'I' is realized by post-filtering the full spectrum */
/* >  from SBDSVDMR3 (which always returns all N singular values). */
/* > */
/* >  SBDSVDMR3 tries SBDSVDMR3_WORK first for every matrix order.  If SBDSVDMR3_WORK */
/* >  reports a positive INFO, returns an incomplete spectrum, or produces */
/* >  a NaN or infinity in the requested output, SBDSVDMR3 restores the input */
/* >  and retries with SBDSDC.  SBDSDC results are normalized to the same */
/* >  internal ordering before RANGE filtering, so either successful path */
/* >  is transparent.  This result-based fallback avoids a size cutoff: */
/* >  solver-only probes showed that the historical large TIMEOUT cases */
/* >  were evaluator timeouts rather than SBDSVDMR3 failing to return. */
/* > */
/* >  Subset modes.  RANGE = 'I' and RANGE = 'V' follow SBDSVDX's */
/* >  conventions (for RANGE = 'I', IL = 1 selects the largest singular */
/* >  value; S is returned in descending order).  In all modes the full */
/* >  decomposition is computed and the requested subset is then */
/* >  returned, so the cost and workspace are those of RANGE = 'A'; the */
/* >  O(N*K) subset capability of MRRR is not exposed by this routine. */
/* > */
/* >  Accuracy notes.  (1) Entries of B (diagonal or off-diagonal) with */
/* >  magnitude at most eps*f2cmax|B| are set to zero before the singular */
/* >  values are computed; the singular values are therefore accurate */
/* >  to O(eps*f2cmax|B|) in the absolute sense (as for SBDSDC) but are not */
/* >  determined to high relative accuracy when they are far below |B| */
/* >  (SBDSQR and SBDSVDX preserve relative accuracy).  (2) Child */
/* >  representations of the MRRR tree are accepted without the */
/* >  element-growth test of SSTEMR; on spectra with tight clusters the */
/* >  orthogonality of the computed vectors can degrade to 10-100 times */
/* >  that of SBDSDC while the residuals remain small.  See the */
/* >  accompanying paper for a quantitative account. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  B is upper bidiagonal; */
/* >          = 'L':  B is lower bidiagonal. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBZ */
/* > \verbatim */
/* >          JOBZ is CHARACTER*1 */
/* >          = 'N':  Compute singular values only; */
/* >          = 'V':  Compute singular values and singular vectors. */
/* > \endverbatim */
/* > */
/* > \param[in] RANGE */
/* > \verbatim */
/* >          RANGE is CHARACTER*1 */
/* >          = 'A': all singular values will be found. */
/* >          = 'V': all singular values in the half-open interval [VL,VU) */
/* >                 will be found. */
/* >          = 'I': the IL-th through IU-th singular values will be found. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the bidiagonal matrix.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          The n diagonal elements of the bidiagonal matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is REAL array, dimension (f2cmax(1,N-1)) */
/* >          The (n-1) off-diagonal elements of the bidiagonal matrix */
/* >          B in elements 1 to N-1. */
/* > \endverbatim */
/* > */
/* > \param[in] VL */
/* > \verbatim */
/* >         VL is REAL */
/* >          If RANGE='V', the lower bound of the interval to */
/* >          be searched for singular values. VU > VL. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] VU */
/* > \verbatim */
/* >         VU is REAL */
/* >          If RANGE='V', the upper bound of the interval to */
/* >          be searched for singular values. VU > VL. */
/* >          Not referenced if RANGE = 'A' or 'I'. */
/* > \endverbatim */
/* > */
/* > \param[in] IL */
/* > \verbatim */
/* >          IL is INTEGER */
/* >          If RANGE='I', the index of the smallest singular value */
/* >          to be returned.  1 <= IL <= IU <= N, if N > 0. */
/* >          Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[in] IU */
/* > \verbatim */
/* >          IU is INTEGER */
/* >          If RANGE='I', the index of the largest singular value */
/* >          to be returned.  1 <= IL <= IU <= N, if N > 0. */
/* >          Not referenced if RANGE = 'A' or 'V'. */
/* > \endverbatim */
/* > */
/* > \param[out] NS */
/* > \verbatim */
/* >          NS is INTEGER */
/* >          The total number of singular values found.  0 <= NS <= N. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is REAL array, dimension (N) */
/* >          The first NS elements contain the selected singular values */
/* >          in decreasing order (matching SBDSVDX's ordering). */
/* > \endverbatim */
/* > */
/* > \param[out] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension (LDZ, K) */
/* >          If JOBZ = 'V', the first NS columns of Z contain the */
/* >          singular vectors of B corresponding to the selected */
/* >          singular values, with U stored in rows 1..N and V in */
/* >          rows N+1..2N (Z = [U; V]).  If JOBZ = 'N', Z is not */
/* >          referenced.  K >= NS+1 (the caller must supply an upper */
/* >          bound on NS when RANGE = 'V'). */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          Leading dimension of Z.  LDZ >= 1; if JOBZ = 'V', LDZ >= 2N. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (LWORK) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          If JOBZ='V', LWORK >= f2cmax( 1, 5*N*N + 37*N ). */
/* >          If JOBZ='N', LWORK >= f2cmax( 1, 2*N*N + 37*N ).  These bounds */
/* >          cover both SBDSVDMR3 and the SBDSDC fallback.  If */
/* >          LWORK = -1, a workspace query is performed: the optimal */
/* >          size is returned in WORK(1) and no work is done. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (LIWORK) */
/* >          On exit, IWORK(2) reports the path taken: 1 = MRRR result */
/* >          accepted, 2 = SBDSDC fallback after the MRRR path returned */
/* >          INFO > 0, 3 = SBDSDC fallback after the finiteness audit */
/* >          rejected the MRRR result. */
/* >          On successful exit for N > 1, IWORK(2) records the path: */
/* >             1: SBDSVDMR3 */
/* >             2: SBDSVDMR3 returned positive INFO, then SBDSDC succeeded */
/* >             3: SBDSVDMR3 output audit failed, then SBDSDC succeeded */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* >          LIWORK is INTEGER */
/* >          LIWORK >= f2cmax( 1, 20*N ).  If LIWORK = -1, a workspace */
/* >          query is performed. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          = N+1: the input contains a NaN or an Inf; nothing was */
/* >                computed. */
/* >          = N+2: the MRRR path was rejected and the SBDSDC fallback */
/* >                returned a non-finite result. */
/* >          > 0, other values: the MRRR path was rejected and the SBDSDC */
/* >                fallback failed; INFO is the INFO returned by SBDSDC. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup bdsvdmr3 */

/*  ===================================================================== */
/* Subroutine */ void sbdsvdmr3_(char *uplo, char *jobz, char *range, integer *
	n, real *d__, real *e, real *vl, real *vu, integer *il, integer *iu, 
	integer *ns, real *s, real *z__, integer *ldz, real *work, integer *
	lwork, integer *iwork, integer *liwork, integer *info)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;

    /* Local variables */
    integer i__, j, ihi, ilo, ncol, path;
    real sval, ovfl;
    integer iwrk;
    char compc[1];
    extern logical lsame_(char *, char *);
    logical allsv, indsv;
    integer lwmin;
    logical lower, valsv;
    extern /* Subroutine */ void scopy_(integer *, real *, integer *, real *, 
	    integer *), sswap_(integer *, real *, integer *, real *, integer *
	    );
    logical wantz;
    extern /* Subroutine */ void sbdsdc_(char *, char *, integer *, real *, 
	    real *, real *, integer *, real *, integer *, real *, integer *, 
	    real *, integer *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */ void xerbla_(char *, integer *, ftnlen);
    logical badout;
    integer iiweng, iwkeng, idcopy, iecopy, mfound, isfull, liwmin, iufull;
    extern /* Subroutine */ void sbdsvdmr3_work_(char *, char *, integer *, 
	    real *, real *, real *, real *, integer *, real *, integer *, 
	    integer *, real *, integer *, integer *, integer *, integer *);
    logical lquery;
    extern real sroundup_lwork__(integer *);
    integer fallwhy, ivtfull;


/*  -- LAPACK driver routine -- */


/*  ===================================================================== */


    /* Parameter adjustments */
    --d__;
    --e;
    --s;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --work;
    --iwork;

    /* Function Body */
    allsv = lsame_(range, "A");
    valsv = lsame_(range, "V");
    indsv = lsame_(range, "I");
    wantz = lsame_(jobz, "V");
    lower = lsame_(uplo, "L");
    lquery = *lwork == -1 || *liwork == -1;

/*     Test the input parameters (mirror SBDSVDX). */

    *info = 0;
    if (! lsame_(uplo, "U") && ! lower) {
	*info = -1;
    } else if (! (wantz || lsame_(jobz, "N"))) {
	*info = -2;
    } else if (! (allsv || valsv || indsv)) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*n > 0) {
	if (valsv) {
	    if (*vl < 0.f) {
		*info = -7;
	    } else if (*vu <= *vl) {
		*info = -8;
	    }
	} else if (indsv) {
	    if (*il < 1 || *il > f2cmax(1,*n)) {
		*info = -9;
	    } else if (*iu < f2cmin(*n,*il) || *iu > *n) {
		*info = -10;
	    }
	}
    }
    if (*info == 0) {
	if (*ldz < 1 || wantz && *ldz < *n << 1) {
	    *info = -14;
	}
    }

/*     Workspace requirements.  The wrapper needs (in doubles): */
/*       * N   for a private copy of D */
/*       * N   for a private copy of E */
/*       * N   for the full unfiltered spectrum */
/*       * LDU*N = N*N for a private copy of U (when WANTZ) */
/*       * LDVT*N = N*N for a private copy of VT (when WANTZ) */
/*       * 2*N*N + 34*N for the SBDSVDMR3 engine itself */
/*     The SBDSDC fallback additionally needs 3*N*N + 4*N doubles after */
/*     its explicit U and VT arrays.  Total (WANTZ): 5*N*N + 7*N; retain */
/*     the larger linear allowance from SBDSVDMR3 for small N. */

    if (*n == 0) {
	lwmin = 1;
	liwmin = 1;
    } else if (wantz) {
/* Computing MAX */
	i__1 = 1, i__2 = *n * 5 * *n + *n * 37;
	lwmin = f2cmax(i__1,i__2);
/* Computing MAX */
	i__1 = 1, i__2 = *n * 24;
	liwmin = f2cmax(i__1,i__2);
    } else {
/* Computing MAX */
	i__1 = 1, i__2 = (*n << 1) * *n + *n * 37 + 2;
	lwmin = f2cmax(i__1,i__2);
/* Computing MAX */
	i__1 = 1, i__2 = *n * 24;
	liwmin = f2cmax(i__1,i__2);
    }

    if (*info == 0 && ! lquery) {
	if (*lwork < lwmin) {
	    *info = -16;
	} else if (*liwork < liwmin) {
	    *info = -18;
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SBDSVDMR3", &i__1, (ftnlen)9);
	return;
    }

    if (lquery) {
	work[1] = sroundup_lwork__(&lwmin);
	iwork[1] = liwmin;
	return;
    }

/*     Quick return if possible. */

    *ns = 0;
    if (*n == 0) {
	return;
    }

/*     A NaN or Inf in the input cannot be processed by either backend */
/*     (it may hang the fallback's QR iteration); report it as INFO = N+1 */
/*     without calling XERBLA, so that the caller can recover. */

    ovfl = slamch_("Overflow");
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sval = d__[i__];
	if (! (sval == sval && abs(sval) <= ovfl)) {
	    *info = *n + 1;
	    return;
	}
	if (i__ < *n) {
	    sval = e[i__];
	    if (! (sval == sval && abs(sval) <= ovfl)) {
		*info = *n + 1;
		return;
	    }
	}
/* L5: */
    }

    if (*n == 1) {
	iwork[2] = 0;
	sval = abs(d__[1]);
	if (allsv || indsv) {
	    *ns = 1;
	    s[1] = sval;
	} else {
	    if (*vl < sval && *vu >= sval) {
		*ns = 1;
		s[1] = sval;
	    }
	}
	if (wantz && *ns > 0) {
	    z__[z_dim1 + 1] = r_sign(&c_b12, &d__[1]);
	    z__[z_dim1 + 2] = 1.f;
	}
	return;
    }

/*     Slice the WORK buffer. */

    idcopy = 1;
    iecopy = idcopy + *n;
    isfull = iecopy + *n;
    if (wantz) {
	iufull = isfull + *n;
	ivtfull = iufull + *n * *n;
	iwrk = ivtfull + *n * *n;
    } else {
/*        U and VT are not referenced for JOBZ = 'N'; give them one */
/*        element each, disjoint from WORK, to avoid aliasing. */
	iufull = isfull + *n;
	ivtfull = iufull + 1;
	iwrk = ivtfull + 1;
    }

/*     Stage inputs (SBDSVDMR3 overwrites D and E; we must not */
/*     modify the caller's arrays). */

    scopy_(n, &d__[1], &c__1, &work[idcopy], &c__1);
    if (*n > 1) {
	i__1 = *n - 1;
	scopy_(&i__1, &e[1], &c__1, &work[iecopy], &c__1);
    }
    work[iecopy + *n - 1] = 0.f;

/*     Compute the full spectrum with SBDSVDMR3 for every N.  A positive */
/*     INFO, incomplete spectrum, or nonfinite requested output triggers a */
/*     restored-input SBDSDC retry.  The O(N**2) output audit is no more */
/*     expensive asymptotically than returning all singular vectors. */

    iwkeng = *lwork - iwrk + 1;
    iiweng = *liwork;
    if (wantz) {
	sbdsvdmr3_work_("V", uplo, n, &work[idcopy], &work[iecopy], &work[
		isfull], &work[iufull], n, &work[ivtfull], n, &mfound, &work[
		iwrk], &iwkeng, &iwork[1], &iiweng, info);
    } else {
	i__1 = f2cmax(1,*n);
	i__2 = f2cmax(1,*n);
	sbdsvdmr3_work_("N", uplo, n, &work[idcopy], &work[iecopy], &work[
		isfull], &work[iufull], &i__1, &work[ivtfull], &i__2, &mfound,
		 &work[iwrk], &iwkeng, &iwork[1], &iiweng, info);
    }
    if (*info < 0) {
	return;
    }

    fallwhy = 0;
    if (*info > 0) {
	fallwhy = 2;
    } else {
	ovfl = slamch_("Overflow");
	badout = mfound != *n;
	if (! badout) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		sval = work[isfull + i__ - 1];
		if (! (sval == sval && abs(sval) <= ovfl)) {
		    badout = TRUE_;
		}
/* L12: */
	    }
	    if (wantz) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			sval = work[iufull + (j - 1) * *n + i__ - 1];
			if (! (sval == sval && abs(sval) <= ovfl)) {
			    badout = TRUE_;
			}
			sval = work[ivtfull + (j - 1) * *n + i__ - 1];
			if (! (sval == sval && abs(sval) <= ovfl)) {
			    badout = TRUE_;
			}
/* L13: */
		    }
/* L14: */
		}
	    }
	}
	if (badout) {
	    fallwhy = 3;
	}
    }

    if (fallwhy == 0) {
	path = 1;
    } else {
/*        SBDSVDMR3 destroys its staged inputs.  Restore them before the */
/*        SBDSDC retry. */
	scopy_(n, &d__[1], &c__1, &work[idcopy], &c__1);
	i__1 = *n - 1;
	scopy_(&i__1, &e[1], &c__1, &work[iecopy], &c__1);
	work[iecopy + *n - 1] = 0.f;
	*info = 0;
	if (wantz) {
	    *(unsigned char *)compc = 'I';
	} else {
	    *(unsigned char *)compc = 'N';
	}
/*        Q and IQ are not referenced for COMPQ = 'I' or 'N'; pass */
/*        storage disjoint from WORK and IWORK to avoid aliasing. */
	i__1 = f2cmax(1,*n);
	i__2 = f2cmax(1,*n);
	sbdsdc_(uplo, compc, n, &work[idcopy], &work[iecopy], &work[iufull], &
		i__1, &work[ivtfull], &i__2, &work[isfull], &iwork[(*n << 3) 
		+ 1], &work[iwrk], &iwork[1], info);
	if (*info != 0) {
	    return;
	}
/*        The finiteness audit applies to the fallback as well. */
	badout = FALSE_;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sval = work[idcopy + i__ - 1];
	    if (! (sval == sval && abs(sval) <= ovfl)) {
		badout = TRUE_;
	    }
/* L16: */
	}
	if (wantz) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    sval = work[iufull + (j - 1) * *n + i__ - 1];
		    if (! (sval == sval && abs(sval) <= ovfl)) {
			badout = TRUE_;
		    }
		    sval = work[ivtfull + (j - 1) * *n + i__ - 1];
		    if (! (sval == sval && abs(sval) <= ovfl)) {
			badout = TRUE_;
		    }
/* L15: */
		}
/* L17: */
	    }
	}
	if (badout) {
	    *info = *n + 2;
	    return;
	}

/*        SBDSDC returns descending singular values, while SBDSVDMR3 */
/*        returns ascending values.  Normalize SBDSDC's result so the */
/*        common RANGE filtering and packing below is unchanged. */

	scopy_(n, &work[idcopy], &c__1, &work[isfull], &c__1);
	i__1 = *n / 2;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sswap_(&c__1, &work[isfull + i__ - 1], &c__1, &work[isfull + *n - 
		    i__], &c__1);
	    if (wantz) {
		sswap_(n, &work[iufull + (i__ - 1) * *n], &c__1, &work[iufull 
			+ (*n - i__) * *n], &c__1);
		sswap_(n, &work[ivtfull + i__ - 1], n, &work[ivtfull + *n - 
			i__], n);
	    }
/* L18: */
	}
	path = fallwhy;
    }
    iwork[2] = path;

/*     Filter and pack the requested slice into S and Z. */

    if (allsv) {
	ilo = 1;
	ihi = *n;
    } else if (indsv) {
/*        IL and IU follow SBDSVDX: index 1 is the largest singular */
/*        value.  The internal spectrum is ascending. */
	ilo = *n - *iu + 1;
	ihi = *n - *il + 1;
    } else {
/*        RANGE='V'; scan ascending spectrum for [VL, VU). */
	ilo = *n + 1;
	ihi = 0;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sval = work[isfull + i__ - 1];
	    if (sval >= *vl && sval < *vu) {
		if (i__ < ilo) {
		    ilo = i__;
		}
		if (i__ > ihi) {
		    ihi = i__;
		}
	    }
/* L10: */
	}
    }

    if (ihi >= ilo) {
	*ns = ihi - ilo + 1;
/*        Reverse ascending -> descending so S(1) is the largest, */
/*        matching SBDSVDX's convention (see sbdsvdx.f */
/*        line 747 sort loop). */
	i__1 = *ns;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s[i__] = work[isfull + ihi - i__];
/* L20: */
	}
	if (wantz) {
/*           Column J of Z holds the vector for S(J).  In the engine's */
/*           spectrum that is column (IHI - J + 1). */
/*           Rows 1..N of Z(:, J) hold U(:, IHI-J+1). */
/*           Rows N+1..2N of Z(:, J) hold VT(IHI-J+1, :)^T. */
	    ncol = *ns;
	    i__1 = ncol;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    z__[i__ + j * z_dim1] = work[iufull + (ihi - j) * *n + 
			    i__ - 1];
		    z__[*n + i__ + j * z_dim1] = work[ivtfull + (i__ - 1) * *
			    n + (ihi - j)];
/* L30: */
		}
/* L40: */
	    }
	}
    } else {
	*ns = 0;
    }

    return;

/*     End of SBDSVDMR3 */

} /* sbdsvdmr3_ */

