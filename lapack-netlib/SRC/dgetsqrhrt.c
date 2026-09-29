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

/* > \brief \b DGETSQRHRT */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGETSQRHRT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetsqr
hrt.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetsqr
hrt.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetsqr
hrt.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGETSQRHRT( M, N, MB1, NB1, NB2, A, LDA, T, LDT, WORK, */
/*      $                       LWORK, INFO ) */
/*       IMPLICIT NONE */

/*       INTEGER           INFO, LDA, LDT, LWORK, M, N, NB1, NB2, MB1 */
/*       DOUBLE PRECISION  A( LDA, * ), T( LDT, * ), WORK( * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGETSQRHRT computes a NB2-sized column blocked QR-factorization */
/* > of a real M-by-N matrix A with M >= N, */
/* > */
/* >    A = Q * R. */
/* > */
/* > The routine uses internally a NB1-sized column blocked and MB1-sized */
/* > row blocked TSQR-factorization and perfors the reconstruction */
/* > of the Householder vectors from the TSQR output. The routine also */
/* > converts the R_tsqr factor from the TSQR-factorization output into */
/* > the R factor that corresponds to the Householder QR-factorization, */
/* > */
/* >    A = Q_tsqr * R_tsqr = Q * R. */
/* > */
/* > The output Q and R factors are stored in the same format as in DGEQRT */
/* > (Q is in blocked compact WY-representation). See the documentation */
/* > of DGEQRT for more details on the format. */
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
/* >          The number of columns of the matrix A. M >= N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] MB1 */
/* > \verbatim */
/* >          MB1 is INTEGER */
/* >          The row block size to be used in the blocked TSQR. */
/* >          MB1 > N. */
/* > \endverbatim */
/* > */
/* > \param[in] NB1 */
/* > \verbatim */
/* >          NB1 is INTEGER */
/* >          The column block size to be used in the blocked TSQR. */
/* >          N >= NB1 >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in] NB2 */
/* > \verbatim */
/* >          NB2 is INTEGER */
/* >          The block size to be used in the blocked QR that is */
/* >          output. NB2 >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* > */
/* >          On entry: an M-by-N matrix A. */
/* > */
/* >          On exit: */
/* >           a) the elements on and above the diagonal */
/* >              of the array contain the N-by-N upper-triangular */
/* >              matrix R corresponding to the Householder QR; */
/* >           b) the elements below the diagonal represent Q by */
/* >              the columns of blocked V (compact WY-representation). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= f2cmax(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          T is DOUBLE PRECISION array, dimension (LDT,N)) */
/* >          The upper triangular block reflectors stored in compact form */
/* >          as a sequence of upper triangular blocks. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T.  LDT >= NB2. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          The dimension of the array WORK. */
/* >          LWORK >= MAX( LWT + LW1, MAX( LWT+N*N+LW2, LWT+N*N+N ) ), */
/* >          where */
/* >             NUM_ALL_ROW_BLOCKS = CEIL((M-N)/(MB1-N)), */
/* >             NB1LOCAL = MIN(NB1,N). */
/* >             LWT = NUM_ALL_ROW_BLOCKS * N * NB1LOCAL, */
/* >             LW1 = NB1LOCAL * N, */
/* >             LW2 = NB1LOCAL * MAX( NB1LOCAL, ( N - NB1LOCAL ) ), */
/* >          If LWORK = -1, then a workspace query is assumed. */
/* >          The routine only calculates the optimal size of the WORK */
/* >          array, returns this value as the first entry of the WORK */
/* >          array, and no error message related to LWORK is issued */
/* >          by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup doubleOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* > November 2020, Igor Kozachenko, */
/* >                Computer Science Division, */
/* >                University of California, Berkeley */
/* > */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ void dgetsqrhrt_(integer *m, integer *n, integer *mb1, 
	integer *nb1, integer *nb2, doublereal *a, integer *lda, doublereal *
	t, integer *ldt, doublereal *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    integer ldwt, lworkopt, i__, j;
    extern /* Subroutine */ void dorgtsqr_row_(integer *, integer *, integer *
	    , integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *);
    integer iinfo;
    extern /* Subroutine */ void dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dorhr_col_(integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    extern void xerbla_(char *, integer *, ftnlen);
    logical lquery;
    integer lw1, lw2, num_all_row_blocks__, lwt;
    extern /* Subroutine */ void dlatsqr_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *);
    integer nb1local, nb2local;


/*  -- LAPACK computational routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */


/*  ===================================================================== */


/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1 * 1;
    t -= t_offset;
    --work;

    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0 || *m < *n) {
	*info = -2;
    } else if (*mb1 <= *n) {
	*info = -3;
    } else if (*nb1 < 1) {
	*info = -4;
    } else if (*nb2 < 1) {
	*info = -5;
    } else if (*lda < f2cmax(1,*m)) {
	*info = -7;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = f2cmin(*nb2,*n);
	if (*ldt < f2cmax(i__1,i__2)) {
	    *info = -9;
	} else {

/*        Test the input LWORK for the dimension of the array WORK. */
/*        This workspace is used to store array: */
/*        a) Matrix T and WORK for DLATSQR; */
/*        b) N-by-N upper-triangular factor R_tsqr; */
/*        c) Matrix T and array WORK for DORGTSQR_ROW; */
/*        d) Diagonal D for DORHR_COL. */

	    if (*lwork < *n * *n + 1 && ! lquery) {
		*info = -11;
	    } else {

/*           Set block size for column blocks */

		nb1local = f2cmin(*nb1,*n);

/* Computing MAX */
		d__3 = (doublereal) (*m - *n) / (doublereal) (*mb1 - *n) + 
			.5f;
		d__1 = 1., d__2 = d_int(&d__3);
		num_all_row_blocks__ = (integer) f2cmax(d__1,d__2);

/*           Length and leading dimension of WORK array to place */
/*           T array in TSQR. */

		lwt = num_all_row_blocks__ * *n * nb1local;
		ldwt = nb1local;

/*           Length of TSQR work array */

		lw1 = nb1local * *n;

/*           Length of DORGTSQR_ROW work array. */

/* Computing MAX */
		i__1 = nb1local, i__2 = *n - nb1local;
		lw2 = nb1local * f2cmax(i__1,i__2);

/* Computing MAX */
/* Computing MAX */
		i__3 = lwt + *n * *n + lw2, i__4 = lwt + *n * *n + *n;
		i__1 = lwt + lw1, i__2 = f2cmax(i__3,i__4);
		lworkopt = f2cmax(i__1,i__2);

		if (*lwork < f2cmax(1,lworkopt) && ! lquery) {
		    *info = -11;
		}

	    }
	}
    }

/*     Handle error in the input parameters and return workspace query. */

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGETSQRHRT", &i__1, (ftnlen)10);
	return;
    } else if (lquery) {
	work[1] = (doublereal) lworkopt;
	return;
    }

/*     Quick return if possible */

    if (f2cmin(*m,*n) == 0) {
	work[1] = (doublereal) lworkopt;
	return;
    }

    nb2local = f2cmin(*nb2,*n);


/*     (1) Perform TSQR-factorization of the M-by-N matrix A. */

    dlatsqr_(m, n, mb1, &nb1local, &a[a_offset], lda, &work[1], &ldwt, &work[
	    lwt + 1], &lw1, &iinfo);

/*     (2) Copy the factor R_tsqr stored in the upper-triangular part */
/*         of A into the square matrix in the work array */
/*         WORK(LWT+1:LWT+N*N) column-by-column. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	dcopy_(&j, &a[j * a_dim1 + 1], &c__1, &work[lwt + *n * (j - 1) + 1], &
		c__1);
    }

/*     (3) Generate a M-by-N matrix Q with orthonormal columns from */
/*     the result stored below the diagonal in the array A in place. */

    dorgtsqr_row_(m, n, mb1, &nb1local, &a[a_offset], lda, &work[1], &ldwt, &
	    work[lwt + *n * *n + 1], &lw2, &iinfo);

/*     (4) Perform the reconstruction of Householder vectors from */
/*     the matrix Q (stored in A) in place. */

    dorhr_col_(m, n, &nb2local, &a[a_offset], lda, &t[t_offset], ldt, &work[
	    lwt + *n * *n + 1], &iinfo);

/*     (5) Copy the factor R_tsqr stored in the square matrix in the */
/*     work array WORK(LWT+1:LWT+N*N) into the upper-triangular */
/*     part of A. */

/*     (6) Compute from R_tsqr the factor R_hr corresponding to */
/*     the reconstructed Householder vectors, i.e. R_hr = S * R_tsqr. */
/*     This multiplication by the sign matrix S on the left means */
/*     changing the sign of I-th row of the matrix R_tsqr according */
/*     to sign of the I-th diagonal element DIAG(I) of the matrix S. */
/*     DIAG is stored in WORK( LWT+N*N+1 ) from the DORHR_COL output. */

/*     (5) and (6) can be combined in a single loop, so the rows in A */
/*     are accessed only once. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (work[lwt + *n * *n + i__] == -1.) {
	    i__2 = *n;
	    for (j = i__; j <= i__2; ++j) {
		a[i__ + j * a_dim1] = work[lwt + *n * (j - 1) + i__] * -1.;
	    }
	} else {
	    i__2 = *n - i__ + 1;
	    dcopy_(&i__2, &work[lwt + *n * (i__ - 1) + i__], n, &a[i__ + i__ *
		     a_dim1], lda);
	}
    }

    work[1] = (doublereal) lworkopt;
    return;

/*     End of DGETSQRHRT */

} /* dgetsqrhrt_ */

