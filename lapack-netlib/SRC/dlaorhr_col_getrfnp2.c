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

static doublereal c_b3 = 1.;
static integer c__1 = 1;
static doublereal c_b19 = -1.;

/* > \brief \b DLAORHR_COL_GETRFNP2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAORHR_GETRF2NP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaorhr
_col_getrfnp2.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaorhr
_col_getrfnp2.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaorhr
_col_getrfnp2.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*        SUBROUTINE DLAORHR_COL_GETRFNP2( M, N, A, LDA, D, INFO ) */

/*       INTEGER            INFO, LDA, M, N */
/*       DOUBLE PRECISION   A( LDA, * ), D( * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAORHR_COL_GETRFNP2 computes the modified LU factorization without */
/* > pivoting of a real general M-by-N matrix A. The factorization has */
/* > the form: */
/* > */
/* >     A - S = L * U, */
/* > */
/* > where: */
/* >    S is a m-by-n diagonal sign matrix with the diagonal D, so that */
/* >    D(i) = S(i,i), 1 <= i <= f2cmin(M,N). The diagonal D is constructed */
/* >    as D(i)=-SIGN(A(i,i)), where A(i,i) is the value after performing */
/* >    i-1 steps of Gaussian elimination. This means that the diagonal */
/* >    element at each step of "modified" Gaussian elimination is at */
/* >    least one in absolute value (so that division-by-zero not */
/* >    possible during the division by the diagonal element); */
/* > */
/* >    L is a M-by-N lower triangular matrix with unit diagonal elements */
/* >    (lower trapezoidal if M > N); */
/* > */
/* >    and U is a M-by-N upper triangular matrix */
/* >    (upper trapezoidal if M < N). */
/* > */
/* > This routine is an auxiliary routine used in the Householder */
/* > reconstruction routine DORHR_COL. In DORHR_COL, this routine is */
/* > applied to an M-by-N matrix A with orthonormal columns, where each */
/* > element is bounded by one in absolute value. With the choice of */
/* > the matrix S above, one can show that the diagonal element at each */
/* > step of Gaussian elimination is the largest (in absolute value) in */
/* > the column on or below the diagonal, so that no pivoting is required */
/* > for numerical stability [1]. */
/* > */
/* > For more details on the Householder reconstruction algorithm, */
/* > including the modified LU factorization, see [1]. */
/* > */
/* > This is the recursive version of the LU factorization algorithm. */
/* > Denote A - S by B. The algorithm divides the matrix B into four */
/* > submatrices: */
/* > */
/* >        [  B11 | B12  ]  where B11 is n1 by n1, */
/* >    B = [ -----|----- ]        B21 is (m-n1) by n1, */
/* >        [  B21 | B22  ]        B12 is n1 by n2, */
/* >                               B22 is (m-n1) by n2, */
/* >                               with n1 = f2cmin(m,n)/2, n2 = n-n1. */
/* > */
/* > */
/* > The subroutine calls itself to factor B11, solves for B21, */
/* > solves for B12, updates B22, then calls itself to factor B22. */
/* > */
/* > For more details on the recursive LU algorithm, see [2]. */
/* > */
/* > DLAORHR_COL_GETRFNP2 is called to factorize a block by the blocked */
/* > routine DLAORHR_COL_GETRFNP, which uses blocked code calling */
/* . Level 3 BLAS to update the submatrix. However, DLAORHR_COL_GETRFNP2 */
/* > is self-sufficient and can be used without DLAORHR_COL_GETRFNP. */
/* > */
/* > [1] "Reconstructing Householder vectors from tall-skinny QR", */
/* >     G. Ballard, J. Demmel, L. Grigori, M. Jacquelin, H.D. Nguyen, */
/* >     E. Solomonik, J. Parallel Distrib. Comput., */
/* >     vol. 85, pp. 3-31, 2015. */
/* > */
/* > [2] "Recursion leads to automatic variable blocking for dense linear */
/* >     algebra algorithms", F. Gustavson, IBM J. of Res. and Dev., */
/* >     vol. 41, no. 6, pp. 737-755, 1997. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix to be factored. */
/* >          On exit, the factors L and U from the factorization */
/* >          A-S=L*U; the unit diagonal elements of L are not stored. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= f2cmax(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension f2cmin(M,N) */
/* >          The diagonal elements of the diagonal M-by-N sign matrix S, */
/* >          D(i) = S(i,i), where 1 <= i <= f2cmin(M,N). The elements can */
/* >          be only plus or minus one. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* > */
/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2019 */

/* > \ingroup doubleGEcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* > November 2019, Igor Kozachenko, */
/* >                Computer Science Division, */
/* >                University of California, Berkeley */
/* > */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ void dlaorhr_col_getrfnp2_(integer *m, integer *n, 
	doublereal *a, integer *lda, doublereal *d__, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1;

    /* Local variables */
    integer i__;
    extern /* Subroutine */ void dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *);
    integer iinfo;
    doublereal sfmin;
    extern /* Subroutine */ void dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *);
    integer n1, n2;
    extern doublereal dlamch_(char *);
    extern /* Subroutine */ void xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.9.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2019 */


/*  ===================================================================== */


/*     Test the input parameters */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --d__;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < f2cmax(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DLAORHR_COL_GETRFNP2", &i__1, (ftnlen)20);
	return;
    }

/*     Quick return if possible */

    if (f2cmin(*m,*n) == 0) {
	return;
    }
    if (*m == 1) {

/*        One row case, (also recursion termination case), */
/*        use unblocked code */

/*        Transfer the sign */

	d__[1] = -d_sign(&c_b3, &a[a_dim1 + 1]);

/*        Construct the row of U */

	a[a_dim1 + 1] -= d__[1];

    } else if (*n == 1) {

/*        One column case, (also recursion termination case), */
/*        use unblocked code */

/*        Transfer the sign */

	d__[1] = -d_sign(&c_b3, &a[a_dim1 + 1]);

/*        Construct the row of U */

	a[a_dim1 + 1] -= d__[1];

/*        Scale the elements 2:M of the column */

/*        Determine machine safe minimum */

	sfmin = dlamch_("S");

/*        Construct the subdiagonal elements of L */

	if ((d__1 = a[a_dim1 + 1], abs(d__1)) >= sfmin) {
	    i__1 = *m - 1;
	    d__1 = 1. / a[a_dim1 + 1];
	    dscal_(&i__1, &d__1, &a[a_dim1 + 2], &c__1);
	} else {
	    i__1 = *m;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		a[i__ + a_dim1] /= a[a_dim1 + 1];
	    }
	}

    } else {

/*        Divide the matrix B into four submatrices */

	n1 = f2cmin(*m,*n) / 2;
	n2 = *n - n1;

/*        Factor B11, recursive call */

	dlaorhr_col_getrfnp2_(&n1, &n1, &a[a_offset], lda, &d__[1], &iinfo);

/*        Solve for B21 */

	i__1 = *m - n1;
	dtrsm_("R", "U", "N", "N", &i__1, &n1, &c_b3, &a[a_offset], lda, &a[
		n1 + 1 + a_dim1], lda);

/*        Solve for B12 */

	dtrsm_("L", "L", "N", "U", &n1, &n2, &c_b3, &a[a_offset], lda, &a[(n1 
		+ 1) * a_dim1 + 1], lda);

/*        Update B22, i.e. compute the Schur complement */
/*        B22 := B22 - B21*B12 */

	i__1 = *m - n1;
	dgemm_("N", "N", &i__1, &n2, &n1, &c_b19, &a[n1 + 1 + a_dim1], lda, &
		a[(n1 + 1) * a_dim1 + 1], lda, &c_b3, &a[n1 + 1 + (n1 + 1) * 
		a_dim1], lda);

/*        Factor B22, recursive call */

	i__1 = *m - n1;
	dlaorhr_col_getrfnp2_(&i__1, &n2, &a[n1 + 1 + (n1 + 1) * a_dim1], 
		lda, &d__[n1 + 1], &iinfo);

    }
    return;

/*     End of DLAORHR_COL_GETRFNP2 */

} /* dlaorhr_col_getrfnp2__ */

