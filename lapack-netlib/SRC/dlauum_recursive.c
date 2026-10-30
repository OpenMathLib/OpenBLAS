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

static integer c__3 = 3;
static integer c_n1 = -1;
static doublereal c_b12 = 1.;

/* > \brief \b DLAUUM_ computes the product UUH or LHL, where U and L are upper or lower triangular matrices (
recursive algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAUUM_( UPLO, N, A, LDA, INFO ) */

/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, N */
/*       DOUBLE PRECISION   A( LDA, * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAUUM_ computes the product U * U**T or L**T * L, where the triangular */
/* > factor U or L is stored in the upper or lower triangular part of */
/* > the array A. */
/* > */
/* > If UPLO = 'U' or 'u' then the upper triangle of the result is stored, */
/* > overwriting the factor U in A, and the strictly lower triangular part */
/* > of A is not referenced. */
/* > If UPLO = 'L' or 'l' then the lower triangle of the result is stored, */
/* > overwriting the factor L in A, and the strictly upper triangular part */
/* > of A is not referenced. */
/* > */
/* > This is the recursive version of the algorithm */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the triangular factor stored in the array A */
/* >          is upper or lower triangular: */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the triangular factor U or L.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the triangular factor U or L. */
/* >          On exit, if UPLO = 'U', the upper triangle of A is */
/* >          overwritten with the upper triangle of the product U * U**T, */
/* >          and the strictly lower triangular part of A is not referenced. */
/* >          If UPLO = 'L', the lower triangle of A is overwritten with */
/* >          the lower triangle of the product L**T * L, and the */
/* >          strictly upper triangular part of A is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= f2cmax(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -k, the k-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup lauum */

/*  ===================================================================== */
/* Subroutine */ void dlauum_recursive_(char *uplo, integer *n, doublereal *a, integer 
	*lda, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    integer k, nx;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ void dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *);
    logical upper;
    extern /* Subroutine */ void dsyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *), dlauu2_(char *, integer *, 
	    doublereal *, integer *, integer *), xerbla_(char *, 
	    integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);


/*  -- LAPACK auxiliary routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */


/*  ===================================================================== */


/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < f2cmax(1,*n)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DLAUUM_RECURSIVE", &i__1, (ftnlen)16);
	return;
    }

/*     Early termination criteria */

    if (*n == 0) {
	return;
    }

/*     Base Case */

    if (*n == 1) {
	a[a_dim1 + 1] *= a[a_dim1 + 1];
	return;
    }

/*     Determine crossover point for when to bail to level2 */

    nx = ilaenv_(&c__3, "DLAUUM_RECURSIVE", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)16, (
	    ftnlen)1);
    if (*n < nx) {
	dlauu2_(uplo, n, &a[a_offset], lda, info);
	return;
    }

/*     Beginning of executable statements for the recursive case */

    k = *n / 2;
    if (upper) {

/*        We are computing A = ut(U*U**T). */

/*        Break apart U as follows */
/*              |-----------------| */
/*        U =   | U_{11} U_{12}   | */
/*              | 0      U_{22}   | */
/*              |-----------------| */

/*        Where */
/*           U_{11}\in\R^{k\times k} U_{12}\in\R^{  k\times n-k} */
/*                                   U_{22}\in\R^{n-k\times n-k} */

/*        and U_{11},U_{22} are upper triangular and U_{12} is rectangular */

/*        This gives us our operations as */
/*                          |--------------------| |------------------------| */
/*        ut(U * U**T)   =  |  U_{11}   U_{12}   | |  U_{11}**T  0          | */
/*                          |  0        U_{22}   | |  U_{12}**T  U_{22}**T  | */
/*                          |--------------------| |------------------------| */

/*        Thus we get */

/*        U_{11} = U_{11}U_{11}**T + U_{12}U_{12}**T */
/*        U_{12} = U_{12}U_{22}**T */

/*        U_{22} = U_{22}U_{22}**T */

/*        We break these operations apart as follows */

/*        U_{11} = U_{11}U_{11}**T            (This subroutine) */
/*        U_{11} = U_{12}U_{12}**T + U_{11}   (SYRK) */

/*        U_{12} = U_{12}U_{22}**T            (TRMM) */

/*        U_{22} = U_{22}U_{22}**T            (This subroutine) */


/*        Compute U_{11} */

	dlauum_recursive_(uplo, &k, &a[a_offset], lda, info);
	i__1 = *n - k;
	dsyrk_("Upper", "No Transpose", &k, &i__1, &c_b12, &a[(k + 1) * 
		a_dim1 + 1], lda, &c_b12, &a[a_offset], lda);

/*        Compute U_{12} */

	i__1 = *n - k;
	dtrmm_("Right", "Upper", "Transpose", "Non-unit", &k, &i__1, &c_b12, &
		a[k + 1 + (k + 1) * a_dim1], lda, &a[(k + 1) * a_dim1 + 1], 
		lda);

/*        Compute U_{22} */

	i__1 = *n - k;
	dlauum_recursive_(uplo, &i__1, &a[k + 1 + (k + 1) * a_dim1], lda, info);
    } else {

/*        We are computing A = lt(L**T*L). */

/*        Break apart L as follows */
/*              |-----------------| */
/*        L =   | L_{11} 0        | */
/*              | L_{21} L_{22}   | */
/*              |-----------------| */

/*        Where */
/*           L_{11}\in\R^{  k\times k} */
/*           L_{21}\in\R^{n-k\times k} l_{22}\in\R^{n-k\times n-k} */

/*        and L_{11},L_{22} are lower triangular and L_{21} is rectangular */

/*        This gives us our operations as */
/*                          |--------------------------| |-----------------| */
/*        lt(L**T * L)   =  |  L_{11}**T   L_{21}**T   | | L_{11} 0        | */
/*                          |  0           L_{22}**T   | | L_{21} L_{22}   | */
/*                          |--------------------------| |-----------------| */

/*        Thus we get */

/*        L_{11} = L_{11}**T L_{11} + L_{21}**T L_{21} */
/*        L_{21} = L_{22}**T L_{21} */

/*        L_{22} = L_{22}**T L_{22} */

/*        We break these operations apart as follows */

/*        L_{11} = L_{11}**T L_{11}           (This subroutine) */
/*        L_{11} = L_{21}**T L_{21} + L_{11}  (SYRK) */

/*        L_{21} = L_{22}**T L_{21}           (TRMM) */

/*        L_{22} = L_{22}**T L_{22}           (This subroutine) */


/*        Compute L_{11} */

	dlauum_recursive_(uplo, &k, &a[a_offset], lda, info);
	i__1 = *n - k;
	dsyrk_("Lower", "Transpose", &k, &i__1, &c_b12, &a[k + 1 + a_dim1], 
		lda, &c_b12, &a[a_offset], lda);

/*        Compute L_{21} */

	i__1 = *n - k;
	dtrmm_("Left", "Lower", "Transpose", "Non-Unit", &i__1, &k, &c_b12, &
		a[k + 1 + (k + 1) * a_dim1], lda, &a[k + 1 + a_dim1], lda);

/*        Compute L_{22} */

	i__1 = *n - k;
	dlauum_recursive_(uplo, &i__1, &a[k + 1 + (k + 1) * a_dim1], lda, info);
    }
    return;
} /* dlauum___ */

