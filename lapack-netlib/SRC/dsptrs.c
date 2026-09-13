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

static doublereal c_b7 = -1.;
static integer c__1 = 1;
static doublereal c_b19 = 1.;

/* > \brief \b DSPTRS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSPTRS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsptrs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsptrs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsptrs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO ) */

/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDB, N, NRHS */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   AP( * ), B( LDB, * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSPTRS solves a system of linear equations A*X = B with a real */
/* > symmetric matrix A stored in packed format using the factorization */
/* > A = U*D*U**T or A = L*D*L**T computed by DSPTRF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the details of the factorization are stored */
/* >          as an upper or lower triangular matrix. */
/* >          = 'U':  Upper triangular, form is A = U*D*U**T; */
/* >          = 'L':  Lower triangular, form is A = L*D*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* >          The block diagonal matrix D and the multipliers used to */
/* >          obtain the factor U or L as computed by DSPTRF, stored as a */
/* >          packed triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by DSPTRF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
/* >          On entry, the right hand side matrix B. */
/* >          On exit, the solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= f2cmax(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ void dsptrs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *ap, integer *ipiv, doublereal *b, integer *ldb, integer *
	info)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1;
    doublereal d__1;

    /* Local variables */
    extern /* Subroutine */ void dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    doublereal akm1k;
    integer j, k;
    extern /* Subroutine */ void dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *);
    doublereal denom;
    extern /* Subroutine */ void dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *), dswap_(integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    logical upper;
    doublereal ak, bk;
    integer kc, kp;
    extern /* Subroutine */ void xerbla_(char *, integer *, ftnlen);
    doublereal akm1, bkm1;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */


/*  ===================================================================== */


    /* Parameter adjustments */
    --ap;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*ldb < f2cmax(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSPTRS", &i__1, (ftnlen)6);
	return;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	return;
    }

    if (upper) {

/*        Solve A*X = B, where A = U*D*U**T. */

/*        First solve U*D*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

	k = *n;
	kc = *n * (*n + 1) / 2 + 1;
L10:

/*        If K < 1, exit from loop. */

	if (k < 1) {
	    goto L30;
	}

	kc -= k;
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

	    kp = ipiv[k];
	    if (kp != k) {
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in column K of A. */

	    i__1 = k - 1;
	    dger_(&i__1, nrhs, &c_b7, &ap[kc], &c__1, &b[k + b_dim1], ldb, &b[
		    b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

	    d__1 = 1. / ap[kc + k - 1];
	    dscal_(nrhs, &d__1, &b[k + b_dim1], ldb);
	    --k;
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K-1 and -IPIV(K). */

	    kp = -ipiv[k];
	    if (kp != k - 1) {
		dswap_(nrhs, &b[k - 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
	    }

/*           Multiply by inv(U(K)), where U(K) is the transformation */
/*           stored in columns K-1 and K of A. */

	    i__1 = k - 2;
	    dger_(&i__1, nrhs, &c_b7, &ap[kc], &c__1, &b[k + b_dim1], ldb, &b[
		    b_dim1 + 1], ldb);
	    i__1 = k - 2;
	    dger_(&i__1, nrhs, &c_b7, &ap[kc - (k - 1)], &c__1, &b[k - 1 + 
		    b_dim1], ldb, &b[b_dim1 + 1], ldb);

/*           Multiply by the inverse of the diagonal block. */

	    akm1k = ap[kc + k - 2];
	    akm1 = ap[kc - 1] / akm1k;
	    ak = ap[kc + k - 1] / akm1k;
	    denom = akm1 * ak - 1.;
	    i__1 = *nrhs;
	    for (j = 1; j <= i__1; ++j) {
		bkm1 = b[k - 1 + j * b_dim1] / akm1k;
		bk = b[k + j * b_dim1] / akm1k;
		b[k - 1 + j * b_dim1] = (ak * bkm1 - bk) / denom;
		b[k + j * b_dim1] = (akm1 * bk - bkm1) / denom;
/* L20: */
	    }
	    kc = kc - k + 1;
	    k += -2;
	}

	goto L10;
L30:

/*        Next solve U**T*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

	k = 1;
	kc = 1;
L40:

/*        If K > N, exit from loop. */

	if (k > *n) {
	    goto L50;
	}

	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(U**T(K)), where U(K) is the transformation */
/*           stored in column K of A. */

	    i__1 = k - 1;
	    dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &ap[kc]
		    , &c__1, &c_b19, &b[k + b_dim1], ldb);

/*           Interchange rows K and IPIV(K). */

	    kp = ipiv[k];
	    if (kp != k) {
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
	    }
	    kc += k;
	    ++k;
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(U**T(K+1)), where U(K+1) is the transformation */
/*           stored in columns K and K+1 of A. */

	    i__1 = k - 1;
	    dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &ap[kc]
		    , &c__1, &c_b19, &b[k + b_dim1], ldb);
	    i__1 = k - 1;
	    dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[b_offset], ldb, &ap[kc 
		    + k], &c__1, &c_b19, &b[k + 1 + b_dim1], ldb);

/*           Interchange rows K and -IPIV(K). */

	    kp = -ipiv[k];
	    if (kp != k) {
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
	    }
	    kc = kc + (k << 1) + 1;
	    k += 2;
	}

	goto L40;
L50:

	;
    } else {

/*        Solve A*X = B, where A = L*D*L**T. */

/*        First solve L*D*X = B, overwriting B with X. */

/*        K is the main loop index, increasing from 1 to N in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

	k = 1;
	kc = 1;
L60:

/*        If K > N, exit from loop. */

	if (k > *n) {
	    goto L80;
	}

	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Interchange rows K and IPIV(K). */

	    kp = ipiv[k];
	    if (kp != k) {
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in column K of A. */

	    if (k < *n) {
		i__1 = *n - k;
		dger_(&i__1, nrhs, &c_b7, &ap[kc + 1], &c__1, &b[k + b_dim1], 
			ldb, &b[k + 1 + b_dim1], ldb);
	    }

/*           Multiply by the inverse of the diagonal block. */

	    d__1 = 1. / ap[kc];
	    dscal_(nrhs, &d__1, &b[k + b_dim1], ldb);
	    kc = kc + *n - k + 1;
	    ++k;
	} else {

/*           2 x 2 diagonal block */

/*           Interchange rows K+1 and -IPIV(K). */

	    kp = -ipiv[k];
	    if (kp != k + 1) {
		dswap_(nrhs, &b[k + 1 + b_dim1], ldb, &b[kp + b_dim1], ldb);
	    }

/*           Multiply by inv(L(K)), where L(K) is the transformation */
/*           stored in columns K and K+1 of A. */

	    if (k < *n - 1) {
		i__1 = *n - k - 1;
		dger_(&i__1, nrhs, &c_b7, &ap[kc + 2], &c__1, &b[k + b_dim1], 
			ldb, &b[k + 2 + b_dim1], ldb);
		i__1 = *n - k - 1;
		dger_(&i__1, nrhs, &c_b7, &ap[kc + *n - k + 2], &c__1, &b[k + 
			1 + b_dim1], ldb, &b[k + 2 + b_dim1], ldb);
	    }

/*           Multiply by the inverse of the diagonal block. */

	    akm1k = ap[kc + 1];
	    akm1 = ap[kc] / akm1k;
	    ak = ap[kc + *n - k + 1] / akm1k;
	    denom = akm1 * ak - 1.;
	    i__1 = *nrhs;
	    for (j = 1; j <= i__1; ++j) {
		bkm1 = b[k + j * b_dim1] / akm1k;
		bk = b[k + 1 + j * b_dim1] / akm1k;
		b[k + j * b_dim1] = (ak * bkm1 - bk) / denom;
		b[k + 1 + j * b_dim1] = (akm1 * bk - bkm1) / denom;
/* L70: */
	    }
	    kc = kc + ((*n - k) << 1) + 1;
	    k += 2;
	}

	goto L60;
L80:

/*        Next solve L**T*X = B, overwriting B with X. */

/*        K is the main loop index, decreasing from N to 1 in steps of */
/*        1 or 2, depending on the size of the diagonal blocks. */

	k = *n;
	kc = *n * (*n + 1) / 2 + 1;
L90:

/*        If K < 1, exit from loop. */

	if (k < 1) {
	    goto L100;
	}

	kc -= *n - k + 1;
	if (ipiv[k] > 0) {

/*           1 x 1 diagonal block */

/*           Multiply by inv(L**T(K)), where L(K) is the transformation */
/*           stored in column K of A. */

	    if (k < *n) {
		i__1 = *n - k;
		dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], 
			ldb, &ap[kc + 1], &c__1, &c_b19, &b[k + b_dim1], ldb);
	    }

/*           Interchange rows K and IPIV(K). */

	    kp = ipiv[k];
	    if (kp != k) {
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
	    }
	    --k;
	} else {

/*           2 x 2 diagonal block */

/*           Multiply by inv(L**T(K-1)), where L(K-1) is the transformation */
/*           stored in columns K-1 and K of A. */

	    if (k < *n) {
		i__1 = *n - k;
		dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], 
			ldb, &ap[kc + 1], &c__1, &c_b19, &b[k + b_dim1], ldb);
		i__1 = *n - k;
		dgemv_("Transpose", &i__1, nrhs, &c_b7, &b[k + 1 + b_dim1], 
			ldb, &ap[kc - (*n - k)], &c__1, &c_b19, &b[k - 1 + 
			b_dim1], ldb);
	    }

/*           Interchange rows K and -IPIV(K). */

	    kp = -ipiv[k];
	    if (kp != k) {
		dswap_(nrhs, &b[k + b_dim1], ldb, &b[kp + b_dim1], ldb);
	    }
	    kc -= *n - k + 2;
	    k += -2;
	}

	goto L90;
L100:
	;
    }

    return;

/*     End of DSPTRS */

} /* dsptrs_ */

