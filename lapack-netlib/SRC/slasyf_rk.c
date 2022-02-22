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

typedef int integer;
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
#define mycycle() continue;
#define myceiling(w) {ceil(w)}
#define myhuge(w) {HUGE_VAL}
//#define mymaxloc_(w,s,e,n) {if (sizeof(*(w)) == sizeof(double)) dmaxloc_((w),*(s),*(e),n); else dmaxloc_((w),*(s),*(e),n);}
#define mymaxloc(w,s,e,n) {dmaxloc_(w,*(s),*(e),n)}

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

static integer c__1 = 1;
static real c_b9 = -1.f;
static real c_b10 = 1.f;

/* > \brief \b SLASYF_RK computes a partial factorization of a real symmetric indefinite matrix using bounded 
Bunch-Kaufman (rook) diagonal pivoting method. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLASYF_RK + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasyf_
rk.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasyf_
rk.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasyf_
rk.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLASYF_RK( UPLO, N, NB, KB, A, LDA, E, IPIV, W, LDW, */
/*                             INFO ) */

/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KB, LDA, LDW, N, NB */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ), E( * ), W( LDW, * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > SLASYF_RK computes a partial factorization of a real symmetric */
/* > matrix A using the bounded Bunch-Kaufman (rook) diagonal */
/* > pivoting method. The partial factorization has the form: */
/* > */
/* > A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or: */
/* >       ( 0  U22 ) (  0   D  ) ( U12**T U22**T ) */
/* > */
/* > A  =  ( L11  0 ) (  D   0  ) ( L11**T L21**T )  if UPLO = 'L', */
/* >       ( L21  I ) (  0  A22 ) (  0       I    ) */
/* > */
/* > where the order of D is at most NB. The actual order is returned in */
/* > the argument KB, and is either NB or NB-1, or N if N <= NB. */
/* > */
/* > SLASYF_RK is an auxiliary routine called by SSYTRF_RK. It uses */
/* > blocked code (calling Level 3 BLAS) to update the submatrix */
/* > A11 (if UPLO = 'U') or A22 (if UPLO = 'L'). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          Specifies whether the upper or lower triangular part of the */
/* >          symmetric matrix A is stored: */
/* >          = 'U':  Upper triangular */
/* >          = 'L':  Lower triangular */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* >          NB is INTEGER */
/* >          The maximum number of columns of the matrix A that should be */
/* >          factored.  NB should be at least 2 to allow for 2-by-2 pivot */
/* >          blocks. */
/* > \endverbatim */
/* > */
/* > \param[out] KB */
/* > \verbatim */
/* >          KB is INTEGER */
/* >          The number of columns of A that were actually factored. */
/* >          KB is either NB-1 or NB, or N if N <= NB. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the symmetric matrix A. */
/* >            If UPLO = 'U': the leading N-by-N upper triangular part */
/* >            of A contains the upper triangular part of the matrix A, */
/* >            and the strictly lower triangular part of A is not */
/* >            referenced. */
/* > */
/* >            If UPLO = 'L': the leading N-by-N lower triangular part */
/* >            of A contains the lower triangular part of the matrix A, */
/* >            and the strictly upper triangular part of A is not */
/* >            referenced. */
/* > */
/* >          On exit, contains: */
/* >            a) ONLY diagonal elements of the symmetric block diagonal */
/* >               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k); */
/* >               (superdiagonal (or subdiagonal) elements of D */
/* >                are stored on exit in array E), and */
/* >            b) If UPLO = 'U': factor U in the superdiagonal part of A. */
/* >               If UPLO = 'L': factor L in the subdiagonal part of A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= f2cmax(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] E */
/* > \verbatim */
/* >          E is REAL array, dimension (N) */
/* >          On exit, contains the superdiagonal (or subdiagonal) */
/* >          elements of the symmetric block diagonal matrix D */
/* >          with 1-by-1 or 2-by-2 diagonal blocks, where */
/* >          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0; */
/* >          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0. */
/* > */
/* >          NOTE: For 1-by-1 diagonal block D(k), where */
/* >          1 <= k <= N, the element E(k) is set to 0 in both */
/* >          UPLO = 'U' or UPLO = 'L' cases. */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          IPIV describes the permutation matrix P in the factorization */
/* >          of matrix A as follows. The absolute value of IPIV(k) */
/* >          represents the index of row and column that were */
/* >          interchanged with the k-th row and column. The value of UPLO */
/* >          describes the order in which the interchanges were applied. */
/* >          Also, the sign of IPIV represents the block structure of */
/* >          the symmetric block diagonal matrix D with 1-by-1 or 2-by-2 */
/* >          diagonal blocks which correspond to 1 or 2 interchanges */
/* >          at each factorization step. */
/* > */
/* >          If UPLO = 'U', */
/* >          ( in factorization order, k decreases from N to 1 ): */
/* >            a) A single positive entry IPIV(k) > 0 means: */
/* >               D(k,k) is a 1-by-1 diagonal block. */
/* >               If IPIV(k) != k, rows and columns k and IPIV(k) were */
/* >               interchanged in the submatrix A(1:N,N-KB+1:N); */
/* >               If IPIV(k) = k, no interchange occurred. */
/* > */
/* > */
/* >            b) A pair of consecutive negative entries */
/* >               IPIV(k) < 0 and IPIV(k-1) < 0 means: */
/* >               D(k-1:k,k-1:k) is a 2-by-2 diagonal block. */
/* >               (NOTE: negative entries in IPIV appear ONLY in pairs). */
/* >               1) If -IPIV(k) != k, rows and columns */
/* >                  k and -IPIV(k) were interchanged */
/* >                  in the matrix A(1:N,N-KB+1:N). */
/* >                  If -IPIV(k) = k, no interchange occurred. */
/* >               2) If -IPIV(k-1) != k-1, rows and columns */
/* >                  k-1 and -IPIV(k-1) were interchanged */
/* >                  in the submatrix A(1:N,N-KB+1:N). */
/* >                  If -IPIV(k-1) = k-1, no interchange occurred. */
/* > */
/* >            c) In both cases a) and b) is always ABS( IPIV(k) ) <= k. */
/* > */
/* >            d) NOTE: Any entry IPIV(k) is always NONZERO on output. */
/* > */
/* >          If UPLO = 'L', */
/* >          ( in factorization order, k increases from 1 to N ): */
/* >            a) A single positive entry IPIV(k) > 0 means: */
/* >               D(k,k) is a 1-by-1 diagonal block. */
/* >               If IPIV(k) != k, rows and columns k and IPIV(k) were */
/* >               interchanged in the submatrix A(1:N,1:KB). */
/* >               If IPIV(k) = k, no interchange occurred. */
/* > */
/* >            b) A pair of consecutive negative entries */
/* >               IPIV(k) < 0 and IPIV(k+1) < 0 means: */
/* >               D(k:k+1,k:k+1) is a 2-by-2 diagonal block. */
/* >               (NOTE: negative entries in IPIV appear ONLY in pairs). */
/* >               1) If -IPIV(k) != k, rows and columns */
/* >                  k and -IPIV(k) were interchanged */
/* >                  in the submatrix A(1:N,1:KB). */
/* >                  If -IPIV(k) = k, no interchange occurred. */
/* >               2) If -IPIV(k+1) != k+1, rows and columns */
/* >                  k-1 and -IPIV(k-1) were interchanged */
/* >                  in the submatrix A(1:N,1:KB). */
/* >                  If -IPIV(k+1) = k+1, no interchange occurred. */
/* > */
/* >            c) In both cases a) and b) is always ABS( IPIV(k) ) >= k. */
/* > */
/* >            d) NOTE: Any entry IPIV(k) is always NONZERO on output. */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is REAL array, dimension (LDW,NB) */
/* > \endverbatim */
/* > */
/* > \param[in] LDW */
/* > \verbatim */
/* >          LDW is INTEGER */
/* >          The leading dimension of the array W.  LDW >= f2cmax(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* > */
/* >          < 0: If INFO = -k, the k-th argument had an illegal value */
/* > */
/* >          > 0: If INFO = k, the matrix A is singular, because: */
/* >                 If UPLO = 'U': column k in the upper */
/* >                 triangular part of A contains all zeros. */
/* >                 If UPLO = 'L': column k in the lower */
/* >                 triangular part of A contains all zeros. */
/* > */
/* >               Therefore D(k,k) is exactly zero, and superdiagonal */
/* >               elements of column k of U (or subdiagonal elements of */
/* >               column k of L ) are all zeros. The factorization has */
/* >               been completed, but the block diagonal matrix D is */
/* >               exactly singular, and division by zero will occur if */
/* >               it is used to solve a system of equations. */
/* > */
/* >               NOTE: INFO only stores the first occurrence of */
/* >               a singularity, any subsequent occurrence of singularity */
/* >               is not stored in INFO even though the factorization */
/* >               always completes. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup singleSYcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > \verbatim */
/* > */
/* >  December 2016,  Igor Kozachenko, */
/* >                  Computer Science Division, */
/* >                  University of California, Berkeley */
/* > */
/* >  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* >                  School of Mathematics, */
/* >                  University of Manchester */
/* > */
/* > \endverbatim */

/*  ===================================================================== */
/* Subroutine */ int slasyf_rk_(char *uplo, integer *n, integer *nb, integer 
	*kb, real *a, integer *lda, real *e, integer *ipiv, real *w, integer *
	ldw, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, w_dim1, w_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1;

    /* Local variables */
    logical done;
    integer imax, jmax, j, k, p;
    real t, alpha;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *), 
	    sgemm_(char *, char *, integer *, integer *, integer *, real *, 
	    real *, integer *, real *, integer *, real *, real *, integer *);
    real sfmin;
    integer itemp;
    extern /* Subroutine */ int sgemv_(char *, integer *, integer *, real *, 
	    real *, integer *, real *, integer *, real *, real *, integer *);
    integer kstep;
    real stemp;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *), sswap_(integer *, real *, integer *, real *, integer *
	    );
    real r1, d11, d12, d21, d22;
    integer jb, ii, jj, kk, kp;
    real absakk;
    integer kw;
    extern real slamch_(char *);
    extern integer isamax_(integer *, real *, integer *);
    real colmax, rowmax;
    integer kkw;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */


/*  ===================================================================== */


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --e;
    --ipiv;
    w_dim1 = *ldw;
    w_offset = 1 + w_dim1 * 1;
    w -= w_offset;

    /* Function Body */
    *info = 0;

/*     Initialize ALPHA for use in choosing pivot block size. */

    alpha = (sqrt(17.f) + 1.f) / 8.f;

/*     Compute machine safe minimum */

    sfmin = slamch_("S");

    if (lsame_(uplo, "U")) {

/*        Factorize the trailing columns of A using the upper triangle */
/*        of A and working backwards, and compute the matrix W = U12*D */
/*        for use in updating A11 */

/*        Initialize the first entry of array E, where superdiagonal */
/*        elements of D are stored */

	e[1] = 0.f;

/*        K is the main loop index, decreasing from N in steps of 1 or 2 */

	k = *n;
L10:

/*        KW is the column of W which corresponds to column K of A */

	kw = *nb + k - *n;

/*        Exit from loop */

	if (k <= *n - *nb + 1 && *nb < *n || k < 1) {
	    goto L30;
	}

	kstep = 1;
	p = k;

/*        Copy column K of A to column KW of W and update it */

	scopy_(&k, &a[k * a_dim1 + 1], &c__1, &w[kw * w_dim1 + 1], &c__1);
	if (k < *n) {
	    i__1 = *n - k;
	    sgemv_("No transpose", &k, &i__1, &c_b9, &a[(k + 1) * a_dim1 + 1],
		     lda, &w[k + (kw + 1) * w_dim1], ldw, &c_b10, &w[kw * 
		    w_dim1 + 1], &c__1);
	}

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

	absakk = (r__1 = w[k + kw * w_dim1], abs(r__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

	if (k > 1) {
	    i__1 = k - 1;
	    imax = isamax_(&i__1, &w[kw * w_dim1 + 1], &c__1);
	    colmax = (r__1 = w[imax + kw * w_dim1], abs(r__1));
	} else {
	    colmax = 0.f;
	}

	if (f2cmax(absakk,colmax) == 0.f) {

/*           Column K is zero or underflow: set INFO and continue */

	    if (*info == 0) {
		*info = k;
	    }
	    kp = k;
	    scopy_(&k, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);

/*           Set E( K ) to zero */

	    if (k > 1) {
		e[k] = 0.f;
	    }

	} else {

/*           ============================================================ */

/*           Test for interchange */

/*           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
/*           (used to handle NaN and Inf) */

	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

		kp = k;

	    } else {

		done = FALSE_;

/*              Loop until pivot found */

L12:

/*                 Begin pivot search loop body */


/*                 Copy column IMAX to column KW-1 of W and update it */

		scopy_(&imax, &a[imax * a_dim1 + 1], &c__1, &w[(kw - 1) * 
			w_dim1 + 1], &c__1);
		i__1 = k - imax;
		scopy_(&i__1, &a[imax + (imax + 1) * a_dim1], lda, &w[imax + 
			1 + (kw - 1) * w_dim1], &c__1);

		if (k < *n) {
		    i__1 = *n - k;
		    sgemv_("No transpose", &k, &i__1, &c_b9, &a[(k + 1) * 
			    a_dim1 + 1], lda, &w[imax + (kw + 1) * w_dim1], 
			    ldw, &c_b10, &w[(kw - 1) * w_dim1 + 1], &c__1);
		}

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

		if (imax != k) {
		    i__1 = k - imax;
		    jmax = imax + isamax_(&i__1, &w[imax + 1 + (kw - 1) * 
			    w_dim1], &c__1);
		    rowmax = (r__1 = w[jmax + (kw - 1) * w_dim1], abs(r__1));
		} else {
		    rowmax = 0.f;
		}

		if (imax > 1) {
		    i__1 = imax - 1;
		    itemp = isamax_(&i__1, &w[(kw - 1) * w_dim1 + 1], &c__1);
		    stemp = (r__1 = w[itemp + (kw - 1) * w_dim1], abs(r__1));
		    if (stemp > rowmax) {
			rowmax = stemp;
			jmax = itemp;
		    }
		}

/*                 Equivalent to testing for */
/*                 ABS( W( IMAX, KW-1 ) ).GE.ALPHA*ROWMAX */
/*                 (used to handle NaN and Inf) */

		if (! ((r__1 = w[imax + (kw - 1) * w_dim1], abs(r__1)) < 
			alpha * rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

		    kp = imax;

/*                    copy column KW-1 of W to column KW of W */

		    scopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * 
			    w_dim1 + 1], &c__1);

		    done = TRUE_;

/*                 Equivalent to testing for ROWMAX.EQ.COLMAX, */
/*                 (used to handle NaN and Inf) */

		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K-1 and IMAX, */
/*                    use 2-by-2 pivot block */

		    kp = imax;
		    kstep = 2;
		    done = TRUE_;
		} else {

/*                    Pivot not found: set params and repeat */

		    p = imax;
		    colmax = rowmax;
		    imax = jmax;

/*                    Copy updated JMAXth (next IMAXth) column to Kth of W */

		    scopy_(&k, &w[(kw - 1) * w_dim1 + 1], &c__1, &w[kw * 
			    w_dim1 + 1], &c__1);

		}

/*                 End pivot search loop body */

		if (! done) {
		    goto L12;
		}

	    }

/*           ============================================================ */

	    kk = k - kstep + 1;

/*           KKW is the column of W which corresponds to column KK of A */

	    kkw = *nb + kk - *n;

	    if (kstep == 2 && p != k) {

/*              Copy non-updated column K to column P */

		i__1 = k - p;
		scopy_(&i__1, &a[p + 1 + k * a_dim1], &c__1, &a[p + (p + 1) * 
			a_dim1], lda);
		scopy_(&p, &a[k * a_dim1 + 1], &c__1, &a[p * a_dim1 + 1], &
			c__1);

/*              Interchange rows K and P in last N-K+1 columns of A */
/*              and last N-K+2 columns of W */

		i__1 = *n - k + 1;
		sswap_(&i__1, &a[k + k * a_dim1], lda, &a[p + k * a_dim1], 
			lda);
		i__1 = *n - kk + 1;
		sswap_(&i__1, &w[k + kkw * w_dim1], ldw, &w[p + kkw * w_dim1],
			 ldw);
	    }

/*           Updated column KP is already stored in column KKW of W */

	    if (kp != kk) {

/*              Copy non-updated column KK to column KP */

		a[kp + k * a_dim1] = a[kk + k * a_dim1];
		i__1 = k - 1 - kp;
		scopy_(&i__1, &a[kp + 1 + kk * a_dim1], &c__1, &a[kp + (kp + 
			1) * a_dim1], lda);
		scopy_(&kp, &a[kk * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &
			c__1);

/*              Interchange rows KK and KP in last N-KK+1 columns */
/*              of A and W */

		i__1 = *n - kk + 1;
		sswap_(&i__1, &a[kk + kk * a_dim1], lda, &a[kp + kk * a_dim1],
			 lda);
		i__1 = *n - kk + 1;
		sswap_(&i__1, &w[kk + kkw * w_dim1], ldw, &w[kp + kkw * 
			w_dim1], ldw);
	    }

	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column KW of W now holds */

/*              W(k) = U(k)*D(k) */

/*              where U(k) is the k-th column of U */

/*              Store U(k) in column k of A */

		scopy_(&k, &w[kw * w_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
			c__1);
		if (k > 1) {
		    if ((r__1 = a[k + k * a_dim1], abs(r__1)) >= sfmin) {
			r1 = 1.f / a[k + k * a_dim1];
			i__1 = k - 1;
			sscal_(&i__1, &r1, &a[k * a_dim1 + 1], &c__1);
		    } else if (a[k + k * a_dim1] != 0.f) {
			i__1 = k - 1;
			for (ii = 1; ii <= i__1; ++ii) {
			    a[ii + k * a_dim1] /= a[k + k * a_dim1];
/* L14: */
			}
		    }

/*                 Store the superdiagonal element of D in array E */

		    e[k] = 0.f;

		}

	    } else {

/*              2-by-2 pivot block D(k): columns KW and KW-1 of W now */
/*              hold */

/*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */

/*              where U(k) and U(k-1) are the k-th and (k-1)-th columns */
/*              of U */

		if (k > 2) {

/*                 Store U(k) and U(k-1) in columns k and k-1 of A */

		    d12 = w[k - 1 + kw * w_dim1];
		    d11 = w[k + kw * w_dim1] / d12;
		    d22 = w[k - 1 + (kw - 1) * w_dim1] / d12;
		    t = 1.f / (d11 * d22 - 1.f);
		    i__1 = k - 2;
		    for (j = 1; j <= i__1; ++j) {
			a[j + (k - 1) * a_dim1] = t * ((d11 * w[j + (kw - 1) *
				 w_dim1] - w[j + kw * w_dim1]) / d12);
			a[j + k * a_dim1] = t * ((d22 * w[j + kw * w_dim1] - 
				w[j + (kw - 1) * w_dim1]) / d12);
/* L20: */
		    }
		}

/*              Copy diagonal elements of D(K) to A, */
/*              copy superdiagonal element of D(K) to E(K) and */
/*              ZERO out superdiagonal entry of A */

		a[k - 1 + (k - 1) * a_dim1] = w[k - 1 + (kw - 1) * w_dim1];
		a[k - 1 + k * a_dim1] = 0.f;
		a[k + k * a_dim1] = w[k + kw * w_dim1];
		e[k] = w[k - 1 + kw * w_dim1];
		e[k - 1] = 0.f;

	    }

/*           End column K is nonsingular */

	}

/*        Store details of the interchanges in IPIV */

	if (kstep == 1) {
	    ipiv[k] = kp;
	} else {
	    ipiv[k] = -p;
	    ipiv[k - 1] = -kp;
	}

/*        Decrease K and return to the start of the main loop */

	k -= kstep;
	goto L10;

L30:

/*        Update the upper triangle of A11 (= A(1:k,1:k)) as */

/*        A11 := A11 - U12*D*U12**T = A11 - U12*W**T */

/*        computing blocks of NB columns at a time */

	i__1 = -(*nb);
	for (j = (k - 1) / *nb * *nb + 1; i__1 < 0 ? j >= 1 : j <= 1; j += 
		i__1) {
/* Computing MIN */
	    i__2 = *nb, i__3 = k - j + 1;
	    jb = f2cmin(i__2,i__3);

/*           Update the upper triangle of the diagonal block */

	    i__2 = j + jb - 1;
	    for (jj = j; jj <= i__2; ++jj) {
		i__3 = jj - j + 1;
		i__4 = *n - k;
		sgemv_("No transpose", &i__3, &i__4, &c_b9, &a[j + (k + 1) * 
			a_dim1], lda, &w[jj + (kw + 1) * w_dim1], ldw, &c_b10,
			 &a[j + jj * a_dim1], &c__1);
/* L40: */
	    }

/*           Update the rectangular superdiagonal block */

	    if (j >= 2) {
		i__2 = j - 1;
		i__3 = *n - k;
		sgemm_("No transpose", "Transpose", &i__2, &jb, &i__3, &c_b9, 
			&a[(k + 1) * a_dim1 + 1], lda, &w[j + (kw + 1) * 
			w_dim1], ldw, &c_b10, &a[j * a_dim1 + 1], lda);
	    }
/* L50: */
	}

/*        Set KB to the number of columns factorized */

	*kb = *n - k;

    } else {

/*        Factorize the leading columns of A using the lower triangle */
/*        of A and working forwards, and compute the matrix W = L21*D */
/*        for use in updating A22 */

/*        Initialize the unused last entry of the subdiagonal array E. */

	e[*n] = 0.f;

/*        K is the main loop index, increasing from 1 in steps of 1 or 2 */

	k = 1;
L70:

/*        Exit from loop */

	if (k >= *nb && *nb < *n || k > *n) {
	    goto L90;
	}

	kstep = 1;
	p = k;

/*        Copy column K of A to column K of W and update it */

	i__1 = *n - k + 1;
	scopy_(&i__1, &a[k + k * a_dim1], &c__1, &w[k + k * w_dim1], &c__1);
	if (k > 1) {
	    i__1 = *n - k + 1;
	    i__2 = k - 1;
	    sgemv_("No transpose", &i__1, &i__2, &c_b9, &a[k + a_dim1], lda, &
		    w[k + w_dim1], ldw, &c_b10, &w[k + k * w_dim1], &c__1);
	}

/*        Determine rows and columns to be interchanged and whether */
/*        a 1-by-1 or 2-by-2 pivot block will be used */

	absakk = (r__1 = w[k + k * w_dim1], abs(r__1));

/*        IMAX is the row-index of the largest off-diagonal element in */
/*        column K, and COLMAX is its absolute value. */
/*        Determine both COLMAX and IMAX. */

	if (k < *n) {
	    i__1 = *n - k;
	    imax = k + isamax_(&i__1, &w[k + 1 + k * w_dim1], &c__1);
	    colmax = (r__1 = w[imax + k * w_dim1], abs(r__1));
	} else {
	    colmax = 0.f;
	}

	if (f2cmax(absakk,colmax) == 0.f) {

/*           Column K is zero or underflow: set INFO and continue */

	    if (*info == 0) {
		*info = k;
	    }
	    kp = k;
	    i__1 = *n - k + 1;
	    scopy_(&i__1, &w[k + k * w_dim1], &c__1, &a[k + k * a_dim1], &
		    c__1);

/*           Set E( K ) to zero */

	    if (k < *n) {
		e[k] = 0.f;
	    }

	} else {

/*           ============================================================ */

/*           Test for interchange */

/*           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX */
/*           (used to handle NaN and Inf) */

	    if (! (absakk < alpha * colmax)) {

/*              no interchange, use 1-by-1 pivot block */

		kp = k;

	    } else {

		done = FALSE_;

/*              Loop until pivot found */

L72:

/*                 Begin pivot search loop body */


/*                 Copy column IMAX to column K+1 of W and update it */

		i__1 = imax - k;
		scopy_(&i__1, &a[imax + k * a_dim1], lda, &w[k + (k + 1) * 
			w_dim1], &c__1);
		i__1 = *n - imax + 1;
		scopy_(&i__1, &a[imax + imax * a_dim1], &c__1, &w[imax + (k + 
			1) * w_dim1], &c__1);
		if (k > 1) {
		    i__1 = *n - k + 1;
		    i__2 = k - 1;
		    sgemv_("No transpose", &i__1, &i__2, &c_b9, &a[k + a_dim1]
			    , lda, &w[imax + w_dim1], ldw, &c_b10, &w[k + (k 
			    + 1) * w_dim1], &c__1);
		}

/*                 JMAX is the column-index of the largest off-diagonal */
/*                 element in row IMAX, and ROWMAX is its absolute value. */
/*                 Determine both ROWMAX and JMAX. */

		if (imax != k) {
		    i__1 = imax - k;
		    jmax = k - 1 + isamax_(&i__1, &w[k + (k + 1) * w_dim1], &
			    c__1);
		    rowmax = (r__1 = w[jmax + (k + 1) * w_dim1], abs(r__1));
		} else {
		    rowmax = 0.f;
		}

		if (imax < *n) {
		    i__1 = *n - imax;
		    itemp = imax + isamax_(&i__1, &w[imax + 1 + (k + 1) * 
			    w_dim1], &c__1);
		    stemp = (r__1 = w[itemp + (k + 1) * w_dim1], abs(r__1));
		    if (stemp > rowmax) {
			rowmax = stemp;
			jmax = itemp;
		    }
		}

/*                 Equivalent to testing for */
/*                 ABS( W( IMAX, K+1 ) ).GE.ALPHA*ROWMAX */
/*                 (used to handle NaN and Inf) */

		if (! ((r__1 = w[imax + (k + 1) * w_dim1], abs(r__1)) < alpha 
			* rowmax)) {

/*                    interchange rows and columns K and IMAX, */
/*                    use 1-by-1 pivot block */

		    kp = imax;

/*                    copy column K+1 of W to column K of W */

		    i__1 = *n - k + 1;
		    scopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + k * 
			    w_dim1], &c__1);

		    done = TRUE_;

/*                 Equivalent to testing for ROWMAX.EQ.COLMAX, */
/*                 (used to handle NaN and Inf) */

		} else if (p == jmax || rowmax <= colmax) {

/*                    interchange rows and columns K+1 and IMAX, */
/*                    use 2-by-2 pivot block */

		    kp = imax;
		    kstep = 2;
		    done = TRUE_;
		} else {

/*                    Pivot not found: set params and repeat */

		    p = imax;
		    colmax = rowmax;
		    imax = jmax;

/*                    Copy updated JMAXth (next IMAXth) column to Kth of W */

		    i__1 = *n - k + 1;
		    scopy_(&i__1, &w[k + (k + 1) * w_dim1], &c__1, &w[k + k * 
			    w_dim1], &c__1);

		}

/*                 End pivot search loop body */

		if (! done) {
		    goto L72;
		}

	    }

/*           ============================================================ */

	    kk = k + kstep - 1;

	    if (kstep == 2 && p != k) {

/*              Copy non-updated column K to column P */

		i__1 = p - k;
		scopy_(&i__1, &a[k + k * a_dim1], &c__1, &a[p + k * a_dim1], 
			lda);
		i__1 = *n - p + 1;
		scopy_(&i__1, &a[p + k * a_dim1], &c__1, &a[p + p * a_dim1], &
			c__1);

/*              Interchange rows K and P in first K columns of A */
/*              and first K+1 columns of W */

		sswap_(&k, &a[k + a_dim1], lda, &a[p + a_dim1], lda);
		sswap_(&kk, &w[k + w_dim1], ldw, &w[p + w_dim1], ldw);
	    }

/*           Updated column KP is already stored in column KK of W */

	    if (kp != kk) {

/*              Copy non-updated column KK to column KP */

		a[kp + k * a_dim1] = a[kk + k * a_dim1];
		i__1 = kp - k - 1;
		scopy_(&i__1, &a[k + 1 + kk * a_dim1], &c__1, &a[kp + (k + 1) 
			* a_dim1], lda);
		i__1 = *n - kp + 1;
		scopy_(&i__1, &a[kp + kk * a_dim1], &c__1, &a[kp + kp * 
			a_dim1], &c__1);

/*              Interchange rows KK and KP in first KK columns of A and W */

		sswap_(&kk, &a[kk + a_dim1], lda, &a[kp + a_dim1], lda);
		sswap_(&kk, &w[kk + w_dim1], ldw, &w[kp + w_dim1], ldw);
	    }

	    if (kstep == 1) {

/*              1-by-1 pivot block D(k): column k of W now holds */

/*              W(k) = L(k)*D(k) */

/*              where L(k) is the k-th column of L */

/*              Store L(k) in column k of A */

		i__1 = *n - k + 1;
		scopy_(&i__1, &w[k + k * w_dim1], &c__1, &a[k + k * a_dim1], &
			c__1);
		if (k < *n) {
		    if ((r__1 = a[k + k * a_dim1], abs(r__1)) >= sfmin) {
			r1 = 1.f / a[k + k * a_dim1];
			i__1 = *n - k;
			sscal_(&i__1, &r1, &a[k + 1 + k * a_dim1], &c__1);
		    } else if (a[k + k * a_dim1] != 0.f) {
			i__1 = *n;
			for (ii = k + 1; ii <= i__1; ++ii) {
			    a[ii + k * a_dim1] /= a[k + k * a_dim1];
/* L74: */
			}
		    }

/*                 Store the subdiagonal element of D in array E */

		    e[k] = 0.f;

		}

	    } else {

/*              2-by-2 pivot block D(k): columns k and k+1 of W now hold */

/*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */

/*              where L(k) and L(k+1) are the k-th and (k+1)-th columns */
/*              of L */

		if (k < *n - 1) {

/*                 Store L(k) and L(k+1) in columns k and k+1 of A */

		    d21 = w[k + 1 + k * w_dim1];
		    d11 = w[k + 1 + (k + 1) * w_dim1] / d21;
		    d22 = w[k + k * w_dim1] / d21;
		    t = 1.f / (d11 * d22 - 1.f);
		    i__1 = *n;
		    for (j = k + 2; j <= i__1; ++j) {
			a[j + k * a_dim1] = t * ((d11 * w[j + k * w_dim1] - w[
				j + (k + 1) * w_dim1]) / d21);
			a[j + (k + 1) * a_dim1] = t * ((d22 * w[j + (k + 1) * 
				w_dim1] - w[j + k * w_dim1]) / d21);
/* L80: */
		    }
		}

/*              Copy diagonal elements of D(K) to A, */
/*              copy subdiagonal element of D(K) to E(K) and */
/*              ZERO out subdiagonal entry of A */

		a[k + k * a_dim1] = w[k + k * w_dim1];
		a[k + 1 + k * a_dim1] = 0.f;
		a[k + 1 + (k + 1) * a_dim1] = w[k + 1 + (k + 1) * w_dim1];
		e[k] = w[k + 1 + k * w_dim1];
		e[k + 1] = 0.f;

	    }

/*           End column K is nonsingular */

	}

/*        Store details of the interchanges in IPIV */

	if (kstep == 1) {
	    ipiv[k] = kp;
	} else {
	    ipiv[k] = -p;
	    ipiv[k + 1] = -kp;
	}

/*        Increase K and return to the start of the main loop */

	k += kstep;
	goto L70;

L90:

/*        Update the lower triangle of A22 (= A(k:n,k:n)) as */

/*        A22 := A22 - L21*D*L21**T = A22 - L21*W**T */

/*        computing blocks of NB columns at a time */

	i__1 = *n;
	i__2 = *nb;
	for (j = k; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
	    i__3 = *nb, i__4 = *n - j + 1;
	    jb = f2cmin(i__3,i__4);

/*           Update the lower triangle of the diagonal block */

	    i__3 = j + jb - 1;
	    for (jj = j; jj <= i__3; ++jj) {
		i__4 = j + jb - jj;
		i__5 = k - 1;
		sgemv_("No transpose", &i__4, &i__5, &c_b9, &a[jj + a_dim1], 
			lda, &w[jj + w_dim1], ldw, &c_b10, &a[jj + jj * 
			a_dim1], &c__1);
/* L100: */
	    }

/*           Update the rectangular subdiagonal block */

	    if (j + jb <= *n) {
		i__3 = *n - j - jb + 1;
		i__4 = k - 1;
		sgemm_("No transpose", "Transpose", &i__3, &jb, &i__4, &c_b9, 
			&a[j + jb + a_dim1], lda, &w[j + w_dim1], ldw, &c_b10,
			 &a[j + jb + j * a_dim1], lda);
	    }
/* L110: */
	}

/*        Set KB to the number of columns factorized */

	*kb = k - 1;

    }

    return 0;

/*     End of SLASYF_RK */

} /* slasyf_rk__ */

