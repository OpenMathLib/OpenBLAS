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
static inline _Dcomplex Cd(doublecomplex *z) {_Dcomplex zz={z->r , z->i};return zz;}
static inline _Dcomplex * _pCd(doublecomplex *z) {return (_Dcomplex*)z;}
#else
static inline _Complex double Cd(doublecomplex *z) {return z->r + z->i*_Complex_I;}
static inline _Complex double * _pCd(doublecomplex *z) {return (_Complex double*)z;}
#endif
#define pCd(z) (*_pCd(z))
typedef int logical;
typedef short int shortlogical;
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
#define c_abs(z) (cabsf(Cf(z)))
#define c_cos(R,Z) { pCf(R)=ccos(Cf(Z)); }
#ifdef _MSC_VER
#define c_div(c, a, b) {Cf(c)._Val[0] = (Cf(a)._Val[0]/Cf(b)._Val[0]); Cf(c)._Val[1]=(Cf(a)._Val[1]/Cf(b)._Val[1]);}
#define z_div(c, a, b) {Cd(c)._Val[0] = (Cd(a)._Val[0]/Cd(b)._Val[0]); Cd(c)._Val[1]=(Cd(a)._Val[1]/Cd(b)._Val[1]);}
#else
#define c_div(c, a, b) {pCf(c) = Cf(a)/Cf(b);}
#define z_div(c, a, b) {pCd(c) = Cd(a)/Cd(b);}
#endif
#define c_exp(R, Z) {pCf(R) = cexpf(Cf(Z));}
#define c_log(R, Z) {pCf(R) = clogf(Cf(Z));}
#define c_sin(R, Z) {pCf(R) = csinf(Cf(Z));}
//#define c_sqrt(R, Z) {*(R) = csqrtf(Cf(Z));}
#define c_sqrt(R, Z) {pCf(R) = csqrtf(Cf(Z));}
#define d_abs(x) (fabs(*(x)))
#define d_acos(x) (acos(*(x)))
#define d_asin(x) (asin(*(x)))
#define d_atan(x) (atan(*(x)))
#define d_atn2(x, y) (atan2(*(x),*(y)))
#define d_cnjg(R, Z) { pCd(R) = conj(Cd(Z)); }
#define r_cnjg(R, Z) { pCf(R) = conjf(Cf(Z)); }
#define d_cos(x) (cos(*(x)))
#define d_cosh(x) (cosh(*(x)))
#define d_dim(__a, __b) ( *(__a) > *(__b) ? *(__a) - *(__b) : 0.0 )
#define d_exp(x) (exp(*(x)))
#define d_imag(z) (cimag(Cd(z)))
#define r_imag(z) (cimagf(Cf(z)))
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
#define s_cat(lpp, rpp, rnp, np, llp) { 	ftnlen i, nc, ll; char *f__rp, *lp; 	ll = (llp); lp = (lpp); 	for(i=0; i < (int)*(np); ++i) {         	nc = ll; 	        if((rnp)[i] < nc) nc = (rnp)[i]; 	        ll -= nc;         	f__rp = (rpp)[i]; 	        while(--nc >= 0) *lp++ = *(f__rp)++;         } 	while(--ll >= 0) *lp++ = ' '; }
#define s_cmp(a,b,c,d) ((integer)strncmp((a),(b),f2cmin((c),(d))))
#define s_copy(A,B,C,D) { int __i,__m; for (__i=0, __m=f2cmin((C),(D)); __i<__m && (B)[__i] != 0; ++__i) (A)[__i] = (B)[__i]; }
#define sig_die(s, kill) { exit(1); }
#define s_stop(s, n) {exit(0);}
#define z_abs(z) (cabs(Cd(z)))
#define z_exp(R, Z) {pCd(R) = cexp(Cd(Z));}
#define z_sqrt(R, Z) {pCd(R) = csqrt(Cd(Z));}

/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};
static integer c__1 = 1;

/* > \brief \b ZLARF1F applies an elementary reflector to a general rectangular */
/*              matrix assuming v(1) = 1. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > Download ZLARF1F + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarf1f
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarf1f
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarf1f
.f"> */
/* > [TXT]</a> */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLARF1F( SIDE, M, N, V, INCV, TAU, C, LDC, WORK ) */

/*       CHARACTER          SIDE */
/*       INTEGER            INCV, LDC, M, N */
/*       COMPLEX*16         TAU */
/*       COMPLEX*16         C( LDC, * ), V( * ), WORK( * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLARF1F applies a complex elementary reflector H to a real m by n matrix */
/* > C, from either the left or the right. H is represented in the form */
/* > */
/* >       H = I - tau * v * v**H */
/* > */
/* > where tau is a complex scalar and v is a complex vector. */
/* > */
/* > If tau = 0, then H is taken to be the unit matrix. */
/* > */
/* > To apply H**H, supply conjg(tau) instead */
/* > tau. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': form  H * C */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* >          V is COMPLEX*16 array, dimension */
/* >                     (1 + (M-1)*abs(INCV)) if SIDE = 'L' */
/* >                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R' */
/* >          The vector v in the representation of H. V is not used if */
/* >          TAU = 0. V(1) is not referenced or modified. */
/* > \endverbatim */
/* > */
/* > \param[in] INCV */
/* > \verbatim */
/* >          INCV is INTEGER */
/* >          The increment between elements of v. INCV <> 0. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX*16 */
/* >          The value tau in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX*16 array, dimension (LDC,N) */
/* >          On entry, the m by n matrix C. */
/* >          On exit, C is overwritten by the matrix H * C if SIDE = 'L', */
/* >          or C * H if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C. LDC >= f2cmax(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension */
/* >                         (N) if SIDE = 'L' */
/* >                      or (M) if SIDE = 'R' */
/* > \endverbatim */
/*  To take advantage of the fact that v(1) = 1, we do the following */
/*     v = [ 1 v_2 ]**T */
/*     If SIDE='L' */
/*           |-----| */
/*           | C_1 | */
/*        C =| C_2 | */
/*           |-----| */
/*        C_1\in\mathbb{C}^{1\times n}, C_2\in\mathbb{C}^{m-1\times n} */
/*        So we compute: */
/*        C = HC   = (I - \tau vv**T)C */
/*                 = C - \tau vv**T C */
/*        w = C**T v  = [ C_1**T C_2**T ] [ 1 v_2 ]**T */
/*                    = C_1**T + C_2**T v ( ZGEMM then ZAXPYC-like ) */
/*        C  = C - \tau vv**T C */
/*           = C - \tau vw**T */
/*        Giving us   C_1 = C_1 - \tau w**T ( ZAXPYC-like ) */
/*                 and */
/*                    C_2 = C_2 - \tau v_2w**T ( ZGERC ) */
/*     If SIDE='R' */

/*        C = [ C_1 C_2 ] */
/*        C_1\in\mathbb{C}^{m\times 1}, C_2\in\mathbb{C}^{m\times n-1} */
/*        So we compute: */
/*        C = CH   = C(I - \tau vv**T) */
/*                 = C - \tau Cvv**T */

/*        w = Cv   = [ C_1 C_2 ] [ 1 v_2 ]**T */
/*                 = C_1 + C_2v_2 ( ZGEMM then ZAXPYC-like ) */
/*        C  = C - \tau Cvv**T */
/*           = C - \tau wv**T */
/*        Giving us   C_1 = C_1 - \tau w ( ZAXPYC-like ) */
/*                 and */
/*                    C_2 = C_2 - \tau wv_2**T ( ZGERC ) */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup larf */

/*  ===================================================================== */
/* Subroutine */ void zlarf1f_(char *side, integer *m, integer *n, 
	doublecomplex *v, integer *incv, doublecomplex *tau, doublecomplex *
	c__, integer *ldc, doublecomplex *work)
{
    /* System generated locals */
    integer c_dim1, c_offset, i__1, i__2, i__3;
    doublecomplex z__1, z__2, z__3;

    /* Local variables */
    integer i__;
    logical applyleft;
    extern logical lsame_(char *, char *);
    integer lastc;
    extern /* Subroutine */ void zgerc_(integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *);
    integer lastv;
    extern /* Subroutine */ void zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern integer ilazlc_(integer *, integer *, doublecomplex *, integer *), 
	    ilazlr_(integer *, integer *, doublecomplex *, integer *);


/*  -- LAPACK auxiliary routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */


/*  ===================================================================== */


    /* Parameter adjustments */
    --v;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;

    /* Function Body */
    applyleft = lsame_(side, "L");
    lastv = 1;
    lastc = 0;
    if (tau->r != 0. || tau->i != 0.) {
/*     Set up variables for scanning V.  LASTV begins pointing to the end */
/*     of V. */
	if (applyleft) {
	    lastv = *m;
	} else {
	    lastv = *n;
	}
	if (*incv > 0) {
	    i__ = (lastv - 1) * *incv + 1;
	} else {
	    i__ = 1;
	}
/*     Look for the last non-zero row in V. */
/*        Since we are assuming that V(1) = 1, and it is not stored, so we */
/*        shouldn't access it. */
	for(;;) { /* while(complicated condition) */
	    i__1 = i__;
	    if (!(lastv > 1 && (v[i__1].r == 0. && v[i__1].i == 0.)))
	    	break;
	    --lastv;
	    i__ -= *incv;
	}
	if (applyleft) {
/*     Scan for the last non-zero column in C(1:lastv,:). */
	    lastc = ilazlc_(&lastv, n, &c__[c_offset], ldc);
	} else {
/*     Scan for the last non-zero row in C(:,1:lastv). */
	    lastc = ilazlr_(m, &lastv, &c__[c_offset], ldc);
	}
    }
    if (lastc == 0) {
	return;
    }
    if (applyleft) {

/*        Form  H * C */

/* Check if m = 1. This means v = 1, So we just need to compu */
/* C := HC = (1-\tau)C. */
	if (lastv == 1) {
	    z__1.r = 1. - tau->r, z__1.i = 0. - tau->i;
	    zscal_(&lastc, &z__1, &c__[c_offset], ldc);
	} else {

/*              w(1:lastc,1) := C(1:lastv,1:lastc)**H * v(1:lastv,1) */

/* (I - tvv**H)C = C - tvv**H C */
/* First compute w**H = v**H c -> w = C**H v */
/* C = [ C_1 C_2 ]**T, v = [1 v_2]**T */
/* w = C_1**H + C_2**Hv_2 */
/* w = C_2**Hv_2 */
	    i__1 = lastv - 1;
	    zgemv_("Conjugate transpose", &i__1, &lastc, &c_b1, &c__[c_dim1 + 
		    2], ldc, &v[*incv + 1], incv, &c_b2, &work[1], &c__1);

/*              w(1:lastc,1) += v(1,1) * C(1,1:lastc)**H */

	    i__1 = lastc;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = i__;
		i__3 = i__;
		d_cnjg(&z__2, &c__[i__ * c_dim1 + 1]);
		z__1.r = work[i__3].r + z__2.r, z__1.i = work[i__3].i + 
			z__2.i;
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
	    }

/*           C(1:lastv,1:lastc) := C(...) - tau * v(1:lastv,1) * w(1:lastc,1)**H */

/* C(1, 1:lastc)   := C(...) - tau * v(1,1) * w(1:lastc,1)**H */
/*                  = C(...) - tau * Conj(w(1:lastc,1)) */
/* This is essentially a zaxpyc */
	    i__1 = lastc;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = i__ * c_dim1 + 1;
		i__3 = i__ * c_dim1 + 1;
		d_cnjg(&z__3, &work[i__]);
		z__2.r = tau->r * z__3.r - tau->i * z__3.i, z__2.i = tau->r * 
			z__3.i + tau->i * z__3.r;
		z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
		c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
	    }

/*        C(2:lastv,1:lastc) += - tau * v(2:lastv,1) * w(1:lastc,1)**H */

	    i__1 = lastv - 1;
	    z__1.r = -tau->r, z__1.i = -tau->i;
	    zgerc_(&i__1, &lastc, &z__1, &v[*incv + 1], incv, &work[1], &c__1,
		     &c__[c_dim1 + 2], ldc);
	}
    } else {

/*        Form  C * H */

/* Check if n = 1. This means v = 1, so we just need to compu */
/* C := CH = C(1-\tau). */
	if (lastv == 1) {
	    z__1.r = 1. - tau->r, z__1.i = 0. - tau->i;
	    zscal_(&lastc, &z__1, &c__[c_offset], &c__1);
	} else {

/*              w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1) */

/* w(1:lastc,1) := C(1:lastc,2:lastv) * v(2:lastv,1) */
	    i__1 = lastv - 1;
	    zgemv_("No transpose", &lastc, &i__1, &c_b1, &c__[(c_dim1 << 1) + 
		    1], ldc, &v[*incv + 1], incv, &c_b2, &work[1], &c__1);
/* w(1:lastc,1) += C(1:lastc,1) v(1,1) = C(1:lastc,1) */
	    zaxpy_(&lastc, &c_b1, &c__[c_offset], &c__1, &work[1], &c__1);

/*              C(1:lastc,1:lastv) := C(...) - tau * w(1:lastc,1) * v(1:lastv,1)**T */

/* C(1:lastc,1)     := C(...) - tau * w(1:lastc,1) * v(1,1 */
/*                   = C(...) - tau * w(1:lastc,1) */
	    z__1.r = -tau->r, z__1.i = -tau->i;
	    zaxpy_(&lastc, &z__1, &work[1], &c__1, &c__[c_offset], &c__1);
/* C(1:lastc,2:lastv) := C(...) - tau * w(1:lastc,1) * v(2 */
	    i__1 = lastv - 1;
	    z__1.r = -tau->r, z__1.i = -tau->i;
	    zgerc_(&lastc, &i__1, &z__1, &work[1], &c__1, &v[*incv + 1], incv,
		     &c__[(c_dim1 << 1) + 1], ldc);
	}
    }
    return;

/*     End of ZLARF1F */

} /* zlarf1f_ */

