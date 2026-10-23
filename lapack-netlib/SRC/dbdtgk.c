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


/* > \brief \b DBDTGK forms the Golub-Kahan tridiagonal matrix of a bidiagonal matrix. */
/* > */
/* > \ingroup bdsvdmr3 */
/* > */
/* > \par Purpose: */
/* > ============= */
/* > */
/* > \verbatim */
/* > forms the Golub-Kahan tridiagonal matrix of a bidiagonal matrix. */
/* > All routines in this file are auxiliary to DBDSVDMR3 and derive */
/* > from the LAPACK 3.0 (1999) MRRR kernels DLARRV, DLAR1V, DLARRB and */
/* > DLARRF rather than from the current DSTEMR kernels. */
/* > \endverbatim */

/* Subroutine */ void dbdtgk_(char *uplo, integer *n, doublereal *d__, 
	doublereal *e, doublereal *dtgk, doublereal *etgk, doublereal *gersch)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    integer i__, m;
    doublereal rad;


/*  -- New auxiliary routine for bidiagonal SVD via TGK-rooted MR^3 -- */
/*     Builds the Tridiagonal Golub-Kahan (TGK) matrix of a bidiagonal B. */
/*     Standalone; does not modify any other routine. */

/*  Purpose */
/*  ======= */

/*  Given a bidiagonal matrix B of order N (diagonal D, off-diagonal E), */
/*  DBDTGK forms the symmetric tridiagonal Golub-Kahan matrix T_GK of */
/*  order 2N.  T_GK has a ZERO diagonal and off-diagonals that interleave */
/*  the entries of B: */
/*       offdiag(T_GK) = [ D(1), E(1), D(2), E(2), ..., E(N-1), D(N) ]. */
/*  Its eigenvalues are +/- the singular values of B; the eigenvector for */
/*  +sigma_i interleaves the left and right singular vectors: */
/*       z(2k-1) ~ u_i(k),  z(2k) ~ v_i(k)   (each scaled by 1/sqrt(2)). */
/*  The Gerschgorin intervals of T_GK are also returned (used by */
/*  DLAR1V_TGK to restrict the twist-index search); since the diagonal is */
/*  zero, interval j is [ -rad_j, +rad_j ] with rad_j the row sum of the */
/*  off-diagonal magnitudes. */

/*  Arguments */
/*  ========= */

/*  UPLO   (input) CHARACTER*1 */
/*         = 'U':  B is upper bidiagonal (E is the superdiagonal); */
/*         = 'L':  B is lower bidiagonal (E is the subdiagonal). */
/*         The TGK off-diagonal ordering is identical for both; UPLO is */
/*         retained for interface clarity and future use (it only affects */
/*         how the recovered z-entries map to u vs v, handled by caller). */

/*  N      (input) INTEGER   The order of B (B is N-by-N). */

/*  D      (input) DOUBLE PRECISION array, dimension (N) */
/*         The diagonal entries of B. */

/*  E      (input) DOUBLE PRECISION array, dimension (N-1) */
/*         The off-diagonal entries of B. */

/*  DTGK   (output) DOUBLE PRECISION array, dimension (2*N) */
/*         The diagonal of T_GK (all zeros). */

/*  ETGK   (output) DOUBLE PRECISION array, dimension (2*N-1) */
/*         The off-diagonal beta_i of T_GK (interleaved D and E). */

/*  GERSCH (output) DOUBLE PRECISION array, dimension (2*(2*N)) */
/*         The 2N Gerschgorin intervals; interval j is */
/*         (GERSCH(2*j-1), GERSCH(2*j)) = (-rad_j, +rad_j). */

/*  ===================================================================== */


    /* Parameter adjustments */
    --gersch;
    --etgk;
    --dtgk;
    --e;
    --d__;

    /* Function Body */
    m = *n << 1;

/*     Diagonal of T_GK is identically zero. */

    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtgk[i__] = 0.;
/* L10: */
    }

/*     Interleave D and E into the off-diagonal of T_GK: */
/*        ETGK(2k-1) = D(k),  ETGK(2k) = E(k). */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	etgk[(i__ << 1) - 1] = d__[i__];
/* L20: */
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	etgk[i__ * 2] = e[i__];
/* L30: */
    }

/*     Gerschgorin intervals: zero diagonal => centered at 0, radius is the */
/*     sum of the magnitudes of the (at most two) adjacent off-diagonals. */

    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rad = 0.;
	if (i__ > 1) {
	    rad += (d__1 = etgk[i__ - 1], abs(d__1));
	}
	if (i__ < m) {
	    rad += (d__1 = etgk[i__], abs(d__1));
	}
	gersch[(i__ << 1) - 1] = -rad;
	gersch[i__ * 2] = rad;
/* L40: */
    }

    return;

/*     End of DBDTGK */

} /* dbdtgk_ */

