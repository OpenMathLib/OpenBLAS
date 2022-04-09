/* f2c.h  --  Standard Fortran to C header file */

/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."

	- From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */

#ifndef F2C_INCLUDE
#define F2C_INCLUDE

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
static inline _Complex float Cf(complex *z) {return z->r + z->i*_Complex_I;}
static inline _Complex double Cd(doublecomplex *z) {return z->r + z->i*_Complex_I;}
static inline _Complex float * _pCf(complex *z) {return (_Complex float*)z;}
static inline _Complex double * _pCd(doublecomplex *z) {return (_Complex double*)z;}
#define pCf(z) (*_pCf(z))
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
#define c_div(c, a, b) {pCf(c) = Cf(a)/Cf(b);}
#define z_div(c, a, b) {pCd(c) = Cd(a)/Cd(b);}
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
#define r_cnjg(R, Z) { pCf(R) = conj(Cf(Z)); }
#define d_cos(x) (cos(*(x)))
#define d_cosh(x) (cosh(*(x)))
#define d_dim(__a, __b) ( *(__a) > *(__b) ? *(__a) - *(__b) : 0.0 )
#define d_exp(x) (exp(*(x)))
#define d_imag(z) (cimag(Cd(z)))
#define r_imag(z) (cimag(Cf(z)))
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
#define pow_zi(p, a, b) {pCd(p) = zpow_ui(Cd(a), *(b));}
#define pow_ci(p, a, b) {pCf(p) = cpow_ui(Cf(a), *(b));}
#define pow_zz(R,A,B) {pCd(R) = cpow(Cd(A),*(B));}
#define s_cat(lpp, rpp, rnp, np, llp) { 	ftnlen i, nc, ll; char *f__rp, *lp; 	ll = (llp); lp = (lpp); 	for(i=0; i < (int)*(np); ++i) {         	nc = ll; 	        if((rnp)[i] < nc) nc = (rnp)[i]; 	        ll -= nc;         	f__rp = (rpp)[i]; 	        while(--nc >= 0) *lp++ = *(f__rp)++;         } 	while(--ll >= 0) *lp++ = ' '; }
#define s_cmp(a,b,c,d) ((integer)strncmp((a),(b),f2cmin((c),(d))))
#define s_copy(A,B,C,D) { int __i,__m; for (__i=0, __m=f2cmin((C),(D)); __i<__m && (B)[__i] != 0; ++__i) (A)[__i] = (B)[__i]; }
#define sig_die(s, kill) { exit(1); }
#define s_stop(s, n) {exit(0);}
static char junk[] = "\n@(#)LIBF77 VERSION 19990503\n";
#define z_abs(z) (cabs(Cd(z)))
#define z_exp(R, Z) {pCd(R) = cexp(Cd(Z));}
#define z_sqrt(R, Z) {pCd(R) = csqrt(Cd(Z));}
#define myexit_() break;
#define mycycle_() continue;
#define myceiling_(w) ceil(w)
#define myhuge_(w) HUGE_VAL
//#define mymaxloc_(w,s,e,n) {if (sizeof(*(w)) == sizeof(double)) dmaxloc_((w),*(s),*(e),n); else dmaxloc_((w),*(s),*(e),n);}
#define mymaxloc_(w,s,e,n) dmaxloc_(w,*(s),*(e),n)

/* procedure parameter types for -A and -C++ */

#define F2C_proc_par_types 1
#ifdef __cplusplus
typedef logical (*L_fp)(...);
#else
typedef logical (*L_fp)();
#endif

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
static _Complex float cpow_ui(_Complex float x, integer n) {
	_Complex float pow=1.0; unsigned long int u;
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
static _Complex double zpow_ui(_Complex double x, integer n) {
	_Complex double pow=1.0; unsigned long int u;
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
static inline void cdotc_(complex *z, integer *n_, complex *x, integer *incx_, complex *y, integer *incy_) {
	integer n = *n_, incx = *incx_, incy = *incy_, i;
	_Complex float zdotc = 0.0;
	if (incx == 1 && incy == 1) {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += conjf(Cf(&x[i])) * Cf(&y[i]);
		}
	} else {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += conjf(Cf(&x[i*incx])) * Cf(&y[i*incy]);
		}
	}
	pCf(z) = zdotc;
}
static inline void zdotc_(doublecomplex *z, integer *n_, doublecomplex *x, integer *incx_, doublecomplex *y, integer *incy_) {
	integer n = *n_, incx = *incx_, incy = *incy_, i;
	_Complex double zdotc = 0.0;
	if (incx == 1 && incy == 1) {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += conj(Cd(&x[i])) * Cd(&y[i]);
		}
	} else {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += conj(Cd(&x[i*incx])) * Cd(&y[i*incy]);
		}
	}
	pCd(z) = zdotc;
}	
static inline void cdotu_(complex *z, integer *n_, complex *x, integer *incx_, complex *y, integer *incy_) {
	integer n = *n_, incx = *incx_, incy = *incy_, i;
	_Complex float zdotc = 0.0;
	if (incx == 1 && incy == 1) {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += Cf(&x[i]) * Cf(&y[i]);
		}
	} else {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += Cf(&x[i*incx]) * Cf(&y[i*incy]);
		}
	}
	pCf(z) = zdotc;
}
static inline void zdotu_(doublecomplex *z, integer *n_, doublecomplex *x, integer *incx_, doublecomplex *y, integer *incy_) {
	integer n = *n_, incx = *incx_, incy = *incy_, i;
	_Complex double zdotc = 0.0;
	if (incx == 1 && incy == 1) {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += Cd(&x[i]) * Cd(&y[i]);
		}
	} else {
		for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
			zdotc += Cd(&x[i*incx]) * Cd(&y[i*incy]);
		}
	}
	pCd(z) = zdotc;
}
#endif
/*  -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/



/* Table of constant values */

static integer c__2 = 2;
static complex c_b6 = {0.f,0.f};

/* > \brief \b CLAHILB */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLAHILB( N, NRHS, A, LDA, X, LDX, B, LDB, WORK, */
/*            INFO, PATH) */

/*       INTEGER N, NRHS, LDA, LDX, LDB, INFO */
/*       REAL WORK(N) */
/*       COMPLEX A(LDA,N), X(LDX, NRHS), B(LDB, NRHS) */
/*       CHARACTER*3        PATH */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAHILB generates an N by N scaled Hilbert matrix in A along with */
/* > NRHS right-hand sides in B and solutions in X such that A*X=B. */
/* > */
/* > The Hilbert matrix is scaled by M = LCM(1, 2, ..., 2*N-1) so that all */
/* > entries are integers.  The right-hand sides are the first NRHS */
/* > columns of M * the identity matrix, and the solutions are the */
/* > first NRHS columns of the inverse Hilbert matrix. */
/* > */
/* > The condition number of the Hilbert matrix grows exponentially with */
/* > its size, roughly as O(e ** (3.5*N)).  Additionally, the inverse */
/* > Hilbert matrices beyond a relatively small dimension cannot be */
/* > generated exactly without extra precision.  Precision is exhausted */
/* > when the largest entry in the inverse Hilbert matrix is greater than */
/* > 2 to the power of the number of bits in the fraction of the data type */
/* > used plus one, which is 24 for single precision. */
/* > */
/* > In single, the generated solution is exact for N <= 6 and has */
/* > small componentwise error for 7 <= N <= 11. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The dimension of the matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The requested number of right-hand sides. */
/* > \endverbatim */
/* > */
/* > \param[out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA, N) */
/* >          The generated scaled Hilbert matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* >          X is COMPLEX array, dimension (LDX, NRHS) */
/* >          The generated exact solutions.  Currently, the first NRHS */
/* >          columns of the inverse Hilbert matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* >          LDX is INTEGER */
/* >          The leading dimension of the array X.  LDX >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB, NRHS) */
/* >          The generated right-hand sides.  Currently, the first NRHS */
/* >          columns of LCM(1, 2, ..., 2*N-1) * the identity matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          = 1: N is too large; the data is still generated but may not */
/* >               be not exact. */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* > */
/* > \param[in] PATH */
/* > \verbatim */
/* >          PATH is CHARACTER*3 */
/* >          The LAPACK path name. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup complex_matgen */

/*  ===================================================================== */
/* Subroutine */ int clahilb_(integer *n, integer *nrhs, complex *a, integer *
	lda, complex *x, integer *ldx, complex *b, integer *ldb, real *work, 
	integer *info, char *path)
{
    /* Initialized data */

    static complex d1[8] = { {-1.f,0.f},{0.f,1.f},{-1.f,-1.f},{0.f,-1.f},{1.f,
	    0.f},{-1.f,1.f},{1.f,1.f},{1.f,-1.f} };
    static complex d2[8] = { {-1.f,0.f},{0.f,-1.f},{-1.f,1.f},{0.f,1.f},{1.f,
	    0.f},{-1.f,-1.f},{1.f,-1.f},{1.f,1.f} };
    static complex invd1[8] = { {-1.f,0.f},{0.f,-1.f},{-.5f,.5f},{0.f,1.f},{
	    1.f,0.f},{-.5f,-.5f},{.5f,-.5f},{.5f,.5f} };
    static complex invd2[8] = { {-1.f,0.f},{0.f,1.f},{-.5f,-.5f},{0.f,-1.f},{
	    1.f,0.f},{-.5f,.5f},{.5f,.5f},{.5f,-.5f} };

    /* System generated locals */
    integer a_dim1, a_offset, x_dim1, x_offset, b_dim1, b_offset, i__1, i__2, 
	    i__3, i__4, i__5;
    real r__1;
    complex q__1, q__2;

    /* Local variables */
    integer i__, j, m, r__;
    char c2[2];
    integer ti, tm;
    extern /* Subroutine */ int claset_(char *, integer *, integer *, complex 
	    *, complex *, complex *, integer *), xerbla_(char *, 
	    integer *);
    extern logical lsamen_(integer *, char *, char *);
    complex tmp;


/*  -- LAPACK test routine (version 3.8.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */


/*  ===================================================================== */
/*     NMAX_EXACT   the largest dimension where the generated data is */
/*                  exact. */
/*     NMAX_APPROX  the largest dimension where the generated data has */
/*                  a small componentwise relative error. */
/*     ??? complex uses how many bits ??? */

/*     d's are generated from random permutation of those eight elements. */
    /* Parameter adjustments */
    --work;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1 * 1;
    x -= x_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;

    /* Function Body */
    s_copy(c2, path + 1, (ftnlen)2, (ftnlen)2);

/*     Test the input arguments */

    *info = 0;
    if (*n < 0 || *n > 11) {
	*info = -1;
    } else if (*nrhs < 0) {
	*info = -2;
    } else if (*lda < *n) {
	*info = -4;
    } else if (*ldx < *n) {
	*info = -6;
    } else if (*ldb < *n) {
	*info = -8;
    }
    if (*info < 0) {
	i__1 = -(*info);
	xerbla_("CLAHILB", &i__1);
	return 0;
    }
    if (*n > 6) {
	*info = 1;
    }

/*     Compute M = the LCM of the integers [1, 2*N-1].  The largest */
/*     reasonable N is small enough that integers suffice (up to N = 11). */
    m = 1;
    i__1 = (*n << 1) - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	tm = m;
	ti = i__;
	r__ = tm % ti;
	while(r__ != 0) {
	    tm = ti;
	    ti = r__;
	    r__ = tm % ti;
	}
	m = m / ti * i__;
    }

/*     Generate the scaled Hilbert matrix in A */
/*     If we are testing SY routines, take */
/*          D1_i = D2_i, else, D1_i = D2_i* */
    if (lsamen_(&c__2, c2, "SY")) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__3 = i__ + j * a_dim1;
		i__4 = j % 8;
		r__1 = (real) m / (i__ + j - 1);
		q__2.r = r__1 * d1[i__4].r, q__2.i = r__1 * d1[i__4].i;
		i__5 = i__ % 8;
		q__1.r = q__2.r * d1[i__5].r - q__2.i * d1[i__5].i, q__1.i = 
			q__2.r * d1[i__5].i + q__2.i * d1[i__5].r;
		a[i__3].r = q__1.r, a[i__3].i = q__1.i;
	    }
	}
    } else {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__3 = i__ + j * a_dim1;
		i__4 = j % 8;
		r__1 = (real) m / (i__ + j - 1);
		q__2.r = r__1 * d1[i__4].r, q__2.i = r__1 * d1[i__4].i;
		i__5 = i__ % 8;
		q__1.r = q__2.r * d2[i__5].r - q__2.i * d2[i__5].i, q__1.i = 
			q__2.r * d2[i__5].i + q__2.i * d2[i__5].r;
		a[i__3].r = q__1.r, a[i__3].i = q__1.i;
	    }
	}
    }

/*     Generate matrix B as simply the first NRHS columns of M * the */
/*     identity. */
    r__1 = (real) m;
    tmp.r = r__1, tmp.i = 0.f;
    claset_("Full", n, nrhs, &c_b6, &tmp, &b[b_offset], ldb);

/*     Generate the true solutions in X.  Because B = the first NRHS */
/*     columns of M*I, the true solutions are just the first NRHS columns */
/*     of the inverse Hilbert matrix. */
    work[1] = (real) (*n);
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	work[j] = work[j - 1] / (j - 1) * (j - 1 - *n) / (j - 1) * (*n + j - 
		1);
    }
/*     If we are testing SY routines, */
/*            take D1_i = D2_i, else, D1_i = D2_i* */
    if (lsamen_(&c__2, c2, "SY")) {
	i__1 = *nrhs;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__3 = i__ + j * x_dim1;
		i__4 = j % 8;
		r__1 = work[i__] * work[j] / (i__ + j - 1);
		q__2.r = r__1 * invd1[i__4].r, q__2.i = r__1 * invd1[i__4].i;
		i__5 = i__ % 8;
		q__1.r = q__2.r * invd1[i__5].r - q__2.i * invd1[i__5].i, 
			q__1.i = q__2.r * invd1[i__5].i + q__2.i * invd1[i__5]
			.r;
		x[i__3].r = q__1.r, x[i__3].i = q__1.i;
	    }
	}
    } else {
	i__1 = *nrhs;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__3 = i__ + j * x_dim1;
		i__4 = j % 8;
		r__1 = work[i__] * work[j] / (i__ + j - 1);
		q__2.r = r__1 * invd2[i__4].r, q__2.i = r__1 * invd2[i__4].i;
		i__5 = i__ % 8;
		q__1.r = q__2.r * invd1[i__5].r - q__2.i * invd1[i__5].i, 
			q__1.i = q__2.r * invd1[i__5].i + q__2.i * invd1[i__5]
			.r;
		x[i__3].r = q__1.r, x[i__3].i = q__1.i;
	    }
	}
    }
    return 0;
} /* clahilb_ */

