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




/* > \brief \b DLA_SYAMV computes a matrix-vector product using a symmetric indefinite matrix to calculate err
or bounds. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLA_SYAMV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dla_sya
mv.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dla_sya
mv.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dla_sya
mv.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLA_SYAMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, */
/*                             INCY ) */

/*       DOUBLE PRECISION   ALPHA, BETA */
/*       INTEGER            INCX, INCY, LDA, N, UPLO */
/*       DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLA_SYAMV  performs the matrix-vector operation */
/* > */
/* >         y := alpha*abs(A)*abs(x) + beta*abs(y), */
/* > */
/* > where alpha and beta are scalars, x and y are vectors and A is an */
/* > n by n symmetric matrix. */
/* > */
/* > This function is primarily used in calculating error bounds. */
/* > To protect against underflow during evaluation, components in */
/* > the resulting vector are perturbed away from zero by (N+1) */
/* > times the underflow threshold.  To prevent unnecessarily large */
/* > errors for block-structure embedded in general matrices, */
/* > "symbolically" zero components are not perturbed.  A zero */
/* > entry is considered "symbolic" if all multiplications involved */
/* > in computing that entry have at least one zero multiplicand. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is INTEGER */
/* >           On entry, UPLO specifies whether the upper or lower */
/* >           triangular part of the array A is to be referenced as */
/* >           follows: */
/* > */
/* >              UPLO = BLAS_UPPER   Only the upper triangular part of A */
/* >                                  is to be referenced. */
/* > */
/* >              UPLO = BLAS_LOWER   Only the lower triangular part of A */
/* >                                  is to be referenced. */
/* > */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the number of columns of the matrix A. */
/* >           N must be at least zero. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is DOUBLE PRECISION . */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension ( LDA, n ). */
/* >           Before entry, the leading m by n part of the array A must */
/* >           contain the matrix of coefficients. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           f2cmax( 1, n ). */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ) */
/* >           Before entry, the incremented array X must contain the */
/* >           vector x. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           X. INCX must not be zero. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is DOUBLE PRECISION . */
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is DOUBLE PRECISION array, dimension */
/* >           ( 1 + ( n - 1 )*abs( INCY ) ) */
/* >           Before entry with BETA non-zero, the incremented array Y */
/* >           must contain the vector y. On exit, Y is overwritten by the */
/* >           updated vector y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >           On entry, INCY specifies the increment for the elements of */
/* >           Y. INCY must not be zero. */
/* >           Unchanged on exit. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2017 */

/* > \ingroup doubleSYcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Level 2 Blas routine. */
/* > */
/* >  -- Written on 22-October-1986. */
/* >     Jack Dongarra, Argonne National Lab. */
/* >     Jeremy Du Croz, Nag Central Office. */
/* >     Sven Hammarling, Nag Central Office. */
/* >     Richard Hanson, Sandia National Labs. */
/* >  -- Modified for the absolute-value product, April 2006 */
/* >     Jason Riedy, UC Berkeley */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ void dla_syamv_(integer *uplo, integer *n, doublereal *alpha,
	 doublereal *a, integer *lda, doublereal *x, integer *incx, 
	doublereal *beta, doublereal *y, integer *incy)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    integer info;
    doublereal temp, safe1;
    integer i__, j;
    logical symb_zero__;
    extern doublereal dlamch_(char *);
    integer iy, jx, kx, ky;
    extern /* Subroutine */ void xerbla_(char *, integer *, ftnlen);
    extern integer ilauplo_(char *);


/*  -- LAPACK computational routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */


/*  ===================================================================== */


/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    info = 0;
    if (*uplo != ilauplo_("U") && *uplo != ilauplo_("L")
	    ) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*lda < f2cmax(1,*n)) {
	info = 5;
    } else if (*incx == 0) {
	info = 7;
    } else if (*incy == 0) {
	info = 10;
    }
    if (info != 0) {
	xerbla_("DLA_SYAMV", &info, (ftnlen)9);
	return;
    }

/*     Quick return if possible. */

    if (*n == 0 || (*alpha == 0. && *beta == 1.)) {
	return;
    }

/*     Set up the start points in  X  and  Y. */

    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (*n - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (*n - 1) * *incy;
    }

/*     Set SAFE1 essentially to be the underflow threshold times the */
/*     number of additions in each row. */

    safe1 = dlamch_("Safe minimum");
    safe1 = (*n + 1) * safe1;

/*     Form  y := alpha*abs(A)*abs(x) + beta*abs(y). */

/*     The O(N^2) SYMB_ZERO tests could be replaced by O(N) queries to */
/*     the inexact flag.  Still doesn't help change the iteration order */
/*     to per-column. */

    iy = ky;
    if (*incx == 1) {
	if (*uplo == ilauplo_("U")) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (*beta == 0.) {
		    symb_zero__ = TRUE_;
		    y[iy] = 0.;
		} else if (y[iy] == 0.) {
		    symb_zero__ = TRUE_;
		} else {
		    symb_zero__ = FALSE_;
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
		}
		if (*alpha != 0.) {
		    i__2 = i__;
		    for (j = 1; j <= i__2; ++j) {
			temp = (d__1 = a[j + i__ * a_dim1], abs(d__1));
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
			y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
		    }
		    i__2 = *n;
		    for (j = i__ + 1; j <= i__2; ++j) {
			temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
			y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
		    }
		}
		if (! symb_zero__) {
		    y[iy] += d_sign(&safe1, &y[iy]);
		}
		iy += *incy;
	    }
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (*beta == 0.) {
		    symb_zero__ = TRUE_;
		    y[iy] = 0.;
		} else if (y[iy] == 0.) {
		    symb_zero__ = TRUE_;
		} else {
		    symb_zero__ = FALSE_;
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
		}
		if (*alpha != 0.) {
		    i__2 = i__;
		    for (j = 1; j <= i__2; ++j) {
			temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
			y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
		    }
		    i__2 = *n;
		    for (j = i__ + 1; j <= i__2; ++j) {
			temp = (d__1 = a[j + i__ * a_dim1], abs(d__1));
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
			y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
		    }
		}
		if (! symb_zero__) {
		    y[iy] += d_sign(&safe1, &y[iy]);
		}
		iy += *incy;
	    }
	}
    } else {
	if (*uplo == ilauplo_("U")) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (*beta == 0.) {
		    symb_zero__ = TRUE_;
		    y[iy] = 0.;
		} else if (y[iy] == 0.) {
		    symb_zero__ = TRUE_;
		} else {
		    symb_zero__ = FALSE_;
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
		}
		jx = kx;
		if (*alpha != 0.) {
		    i__2 = i__;
		    for (j = 1; j <= i__2; ++j) {
			temp = (d__1 = a[j + i__ * a_dim1], abs(d__1));
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
			y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
			jx += *incx;
		    }
		    i__2 = *n;
		    for (j = i__ + 1; j <= i__2; ++j) {
			temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
			y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
			jx += *incx;
		    }
		}
		if (! symb_zero__) {
		    y[iy] += d_sign(&safe1, &y[iy]);
		}
		iy += *incy;
	    }
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (*beta == 0.) {
		    symb_zero__ = TRUE_;
		    y[iy] = 0.;
		} else if (y[iy] == 0.) {
		    symb_zero__ = TRUE_;
		} else {
		    symb_zero__ = FALSE_;
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
		}
		jx = kx;
		if (*alpha != 0.) {
		    i__2 = i__;
		    for (j = 1; j <= i__2; ++j) {
			temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
			y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
			jx += *incx;
		    }
		    i__2 = *n;
		    for (j = i__ + 1; j <= i__2; ++j) {
			temp = (d__1 = a[j + i__ * a_dim1], abs(d__1));
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
			y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
			jx += *incx;
		    }
		}
		if (! symb_zero__) {
		    y[iy] += d_sign(&safe1, &y[iy]);
		}
		iy += *incy;
	    }
	}
    }

    return;

/*     End of DLA_SYAMV */

} /* dla_syamv__ */

