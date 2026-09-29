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
#ifdef _MSC_VER
static inline _Fcomplex Cf(complex *z) {_Fcomplex zz={z->r , z->i}; return zz;}
static inline _Fcomplex * _pCf(complex *z) {return (_Fcomplex*)z;}
#else
static inline _Complex float Cf(complex *z) {return z->r + z->i*_Complex_I;}
static inline _Complex float * _pCf(complex *z) {return (_Complex float*)z;}
#endif
#define pCf(z) (*_pCf(z))

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
#define r_cnjg(R, Z) { pCf(R) = conjf(Cf(Z)); }
#define r_imag(z) (cimagf(Cf(z)))
#ifdef _MSC_VER
#define c_div(c, a, b) {float nenn=crealf(_FCmulcc(Cf(b),conjf(Cf(b)))); _Fcomplex zaehl=_FCmulcc(Cf(a),conjf(Cf(b))); pCf(c)=_FCbuild(crealf(zaehl)/nenn,cimagf(zaehl)/nenn);}
#else
#define c_div(c, a, b) {pCf(c) = Cf(a)/Cf(b);}
#endif
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




/* > \brief \b CGTTRF */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGTTRF + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgttrf.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgttrf.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgttrf.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGTTRF( N, DL, D, DU, DU2, IPIV, INFO ) */

/*       INTEGER            INFO, N */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            D( * ), DL( * ), DU( * ), DU2( * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGTTRF computes an LU factorization of a complex tridiagonal matrix A */
/* > using elimination with partial pivoting and row interchanges. */
/* > */
/* > The factorization has the form */
/* >    A = L * U */
/* > where L is a product of permutation and unit lower bidiagonal */
/* > matrices and U is upper triangular with nonzeros in only the main */
/* > diagonal and first two superdiagonals. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DL */
/* > \verbatim */
/* >          DL is COMPLEX array, dimension (N-1) */
/* >          On entry, DL must contain the (n-1) sub-diagonal elements of */
/* >          A. */
/* > */
/* >          On exit, DL is overwritten by the (n-1) multipliers that */
/* >          define the matrix L from the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is COMPLEX array, dimension (N) */
/* >          On entry, D must contain the diagonal elements of A. */
/* > */
/* >          On exit, D is overwritten by the n diagonal elements of the */
/* >          upper triangular matrix U from the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DU */
/* > \verbatim */
/* >          DU is COMPLEX array, dimension (N-1) */
/* >          On entry, DU must contain the (n-1) super-diagonal elements */
/* >          of A. */
/* > */
/* >          On exit, DU is overwritten by the (n-1) elements of the first */
/* >          super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[out] DU2 */
/* > \verbatim */
/* >          DU2 is COMPLEX array, dimension (N-2) */
/* >          On exit, DU2 is overwritten by the (n-2) elements of the */
/* >          second super-diagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          The pivot indices; for 1 <= i <= n, row i of the matrix was */
/* >          interchanged with row IPIV(i).  IPIV(i) will always be either */
/* >          i or i+1; IPIV(i) = i indicates a row interchange was not */
/* >          required. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -k, the k-th argument had an illegal value */
/* >          > 0:  if INFO = k, U(k,k) is exactly zero. The factorization */
/* >                has been completed, but the factor U is exactly */
/* >                singular, and division by zero will occur if it is used */
/* >                to solve a system of equations. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complexGTcomputational */

/*  ===================================================================== */
/* Subroutine */ void cgttrf_(integer *n, complex *dl, complex *d__, complex *
	du, complex *du2, integer *ipiv, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    real r__1=0., r__2=0., r__3=0., r__4=0.;
    complex q__1={0.,0.}, q__2={0.,0.};

    /* Local variables */
    complex fact, temp;
    integer i__;
    extern /* Subroutine */ void xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */


/*  ===================================================================== */


    /* Parameter adjustments */
    --ipiv;
    --du2;
    --du;
    --d__;
    --dl;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
	i__1 = -(*info);
	xerbla_("CGTTRF", &i__1, (ftnlen)6);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/*     Initialize IPIV(i) = i and DU2(i) = 0 */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ipiv[i__] = i__;
/* L10: */
    }
    i__1 = *n - 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	du2[i__2].r = 0.f, du2[i__2].i = 0.f;
/* L20: */
    }

    i__1 = *n - 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	i__3 = i__;
	if ((r__1 = d__[i__2].r, abs(r__1)) + (r__2 = r_imag(&d__[i__]), abs(
		r__2)) >= (r__3 = dl[i__3].r, abs(r__3)) + (r__4 = r_imag(&dl[
		i__]), abs(r__4))) {

/*           No row interchange required, eliminate DL(I) */

	    i__2 = i__;
	    if ((r__1 = d__[i__2].r, abs(r__1)) + (r__2 = r_imag(&d__[i__]), 
		    abs(r__2)) != 0.f) {
		c_div(&q__1, &dl[i__], &d__[i__]);
		fact.r = q__1.r, fact.i = q__1.i;
		i__2 = i__;
		dl[i__2].r = fact.r, dl[i__2].i = fact.i;
		i__2 = i__ + 1;
		i__3 = i__ + 1;
		i__4 = i__;
		q__2.r = fact.r * du[i__4].r - fact.i * du[i__4].i, q__2.i = 
			fact.r * du[i__4].i + fact.i * du[i__4].r;
		q__1.r = d__[i__3].r - q__2.r, q__1.i = d__[i__3].i - q__2.i;
		d__[i__2].r = q__1.r, d__[i__2].i = q__1.i;
	    }
	} else {

/*           Interchange rows I and I+1, eliminate DL(I) */

	    c_div(&q__1, &d__[i__], &dl[i__]);
	    fact.r = q__1.r, fact.i = q__1.i;
	    i__2 = i__;
	    i__3 = i__;
	    d__[i__2].r = dl[i__3].r, d__[i__2].i = dl[i__3].i;
	    i__2 = i__;
	    dl[i__2].r = fact.r, dl[i__2].i = fact.i;
	    i__2 = i__;
	    temp.r = du[i__2].r, temp.i = du[i__2].i;
	    i__2 = i__;
	    i__3 = i__ + 1;
	    du[i__2].r = d__[i__3].r, du[i__2].i = d__[i__3].i;
	    i__2 = i__ + 1;
	    i__3 = i__ + 1;
	    q__2.r = fact.r * d__[i__3].r - fact.i * d__[i__3].i, q__2.i = 
		    fact.r * d__[i__3].i + fact.i * d__[i__3].r;
	    q__1.r = temp.r - q__2.r, q__1.i = temp.i - q__2.i;
	    d__[i__2].r = q__1.r, d__[i__2].i = q__1.i;
	    i__2 = i__;
	    i__3 = i__ + 1;
	    du2[i__2].r = du[i__3].r, du2[i__2].i = du[i__3].i;
	    i__2 = i__ + 1;
	    q__2.r = -fact.r, q__2.i = -fact.i;
	    i__3 = i__ + 1;
	    q__1.r = q__2.r * du[i__3].r - q__2.i * du[i__3].i, q__1.i = 
		    q__2.r * du[i__3].i + q__2.i * du[i__3].r;
	    du[i__2].r = q__1.r, du[i__2].i = q__1.i;
	    ipiv[i__] = i__ + 1;
	}
/* L30: */
    }
    if (*n > 1) {
	i__ = *n - 1;
	i__1 = i__;
	i__2 = i__;
	if ((r__1 = d__[i__1].r, abs(r__1)) + (r__2 = r_imag(&d__[i__]), abs(
		r__2)) >= (r__3 = dl[i__2].r, abs(r__3)) + (r__4 = r_imag(&dl[
		i__]), abs(r__4))) {
	    i__1 = i__;
	    if ((r__1 = d__[i__1].r, abs(r__1)) + (r__2 = r_imag(&d__[i__]), 
		    abs(r__2)) != 0.f) {
		c_div(&q__1, &dl[i__], &d__[i__]);
		fact.r = q__1.r, fact.i = q__1.i;
		i__1 = i__;
		dl[i__1].r = fact.r, dl[i__1].i = fact.i;
		i__1 = i__ + 1;
		i__2 = i__ + 1;
		i__3 = i__;
		q__2.r = fact.r * du[i__3].r - fact.i * du[i__3].i, q__2.i = 
			fact.r * du[i__3].i + fact.i * du[i__3].r;
		q__1.r = d__[i__2].r - q__2.r, q__1.i = d__[i__2].i - q__2.i;
		d__[i__1].r = q__1.r, d__[i__1].i = q__1.i;
	    }
	} else {
	    c_div(&q__1, &d__[i__], &dl[i__]);
	    fact.r = q__1.r, fact.i = q__1.i;
	    i__1 = i__;
	    i__2 = i__;
	    d__[i__1].r = dl[i__2].r, d__[i__1].i = dl[i__2].i;
	    i__1 = i__;
	    dl[i__1].r = fact.r, dl[i__1].i = fact.i;
	    i__1 = i__;
	    temp.r = du[i__1].r, temp.i = du[i__1].i;
	    i__1 = i__;
	    i__2 = i__ + 1;
	    du[i__1].r = d__[i__2].r, du[i__1].i = d__[i__2].i;
	    i__1 = i__ + 1;
	    i__2 = i__ + 1;
	    q__2.r = fact.r * d__[i__2].r - fact.i * d__[i__2].i, q__2.i = 
		    fact.r * d__[i__2].i + fact.i * d__[i__2].r;
	    q__1.r = temp.r - q__2.r, q__1.i = temp.i - q__2.i;
	    d__[i__1].r = q__1.r, d__[i__1].i = q__1.i;
	    ipiv[i__] = i__ + 1;
	}
    }

/*     Check for a zero on the diagonal of U. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	if ((r__1 = d__[i__2].r, abs(r__1)) + (r__2 = r_imag(&d__[i__]), abs(
		r__2)) == 0.f) {
	    *info = i__;
	    goto L50;
	}
/* L40: */
    }
L50:

    return;

/*     End of CGTTRF */

} /* cgttrf_ */

