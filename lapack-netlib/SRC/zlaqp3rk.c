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
#define c_abs(z) (cabsf(Cf(z)))
#define c_cos(R,Z) { pCf(R)=ccos(Cf(Z)); }
#ifdef _MSC_VER
#define c_div(c, a, b) {float nenn=crealf(_FCmulcc(Cf(b),conjf(Cf(b)))); _Fcomplex zaehl=_FCmulcc(Cf(a),conjf(Cf(b))); pCf(c)=_FCbuild(crealf(zaehl)/nenn,cimagf(zaehl)/nenn);}
#define z_div(c, a, b) {double nenn=creal(_Cmulcc(Cd(b),conj(Cd(b)))); _Dcomplex zaehl=_Cmulcc(Cd(a),conj(Cd(b))); pCd(c)=_Cbuild(creal(zaehl)/nenn,cimag(zaehl)/nenn);}
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
#define z_abs(z) (cabs(Cd(z)))
#define z_exp(R, Z) {pCd(R) = cexp(Cd(Z));}
#define z_sqrt(R, Z) {pCd(R) = csqrt(Cd(Z));}
#define myexit_() break;
#define mycycle_() continue;
#define myceiling_(w) {ceil(w)}
#define myhuge_(w) {HUGE_VAL}
//#define mymaxloc_(w,s,e,n) {if (sizeof(*(w)) == sizeof(double)) dmaxloc_((w),*(s),*(e),n); else dmaxloc_((w),*(s),*(e),n);}
#define mymaxloc_(w,s,e,n) dmaxloc_(w,*(s),*(e),n)

/*  -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/



/*  -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/



/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* Subroutine */ void zlaqp3rk_(integer *m, integer *n, integer *nrhs, integer 
	*ioffset, integer *nb, doublereal *abstol, doublereal *reltol, 
	integer *kp1, doublereal *maxc2nrm, doublecomplex *a, integer *lda, 
	logical *done, integer *kb, doublereal *maxc2nrmk, doublereal *
	relmaxc2nrmk, integer *jpiv, doublecomplex *tau, doublereal *vn1, 
	doublereal *vn2, doublecomplex *auxv, doublecomplex *f, integer *ldf, 
	integer *iwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, f_dim1, f_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Local variables */
    doublereal temp, temp2;
    integer i__, j, k;
    doublereal tol3z;
    integer itemp;
    extern /* Subroutine */ void zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *), zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *);
    integer minmnfact;
    extern /* Subroutine */ void zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    doublereal myhugeval;
    integer minmnupdt;
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *);
    integer if__;
    extern doublereal dlamch_(char *);
    integer kp;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern logical disnan_(doublereal *);
    integer lsticc;
    extern /* Subroutine */ void zlarfg_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *);
    doublereal taunan;
    doublecomplex aik;


/*  -- LAPACK auxiliary routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */


/*  ===================================================================== */


/*     Initialize INFO */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --jpiv;
    --tau;
    --vn1;
    --vn2;
    --auxv;
    f_dim1 = *ldf;
    f_offset = 1 + f_dim1 * 1;
    f -= f_offset;
    --iwork;

    /* Function Body */
    *info = 0;

/*     MINMNFACT in the smallest dimension of the submatrix */
/*     A(IOFFSET+1:M,1:N) to be factorized. */

/* Computing MIN */
    i__1 = *m - *ioffset;
    minmnfact = f2cmin(i__1,*n);
/* Computing MIN */
    i__1 = *m - *ioffset, i__2 = *n + *nrhs;
    minmnupdt = f2cmin(i__1,i__2);
    *nb = f2cmin(*nb,minmnfact);
    tol3z = sqrt(dlamch_("Epsilon"));
    myhugeval = dlamch_("Overflow");

/*     Compute factorization in a while loop over NB columns, */
/*     K is the column index in the block A(1:M,1:N). */

    k = 0;
    lsticc = 0;
    *done = FALSE_;

    while(k < *nb && lsticc == 0) {
	++k;
	i__ = *ioffset + k;

	if (i__ == 1) {

/*           We are at the first column of the original whole matrix A_orig, */
/*           therefore we use the computed KP1 and MAXC2NRM from the */
/*           main routine. */

	    kp = *kp1;

	} else {

/*           Determine the pivot column in K-th step, i.e. the index */
/*           of the column with the maximum 2-norm in the */
/*           submatrix A(I:M,K:N). */

	    i__1 = *n - k + 1;
	    kp = k - 1 + idamax_(&i__1, &vn1[k], &c__1);

/*           Determine the maximum column 2-norm and the relative maximum */
/*           column 2-norm of the submatrix A(I:M,K:N) in step K. */

	    *maxc2nrmk = vn1[kp];

/*           ============================================================ */

/*           Check if the submatrix A(I:M,K:N) contains NaN, set */
/*           INFO parameter to the column number, where the first NaN */
/*           is found and return from the routine. */
/*           We need to check the condition only if the */
/*           column index (same as row index) of the original whole */
/*           matrix is larger than 1, since the condition for whole */
/*           original matrix is checked in the main routine. */

	    if (disnan_(maxc2nrmk)) {

		*done = TRUE_;

/*              Set KB, the number of factorized partial columns */
/*                      that are non-zero in each step in the block, */
/*                      i.e. the rank of the factor R. */
/*              Set IF, the number of processed rows in the block, which */
/*                      is the same as the number of processed rows in */
/*                      the original whole matrix A_orig. */

		*kb = k - 1;
		if__ = i__ - 1;
		*info = *kb + kp;

/*              Set RELMAXC2NRMK to NaN. */

		*relmaxc2nrmk = *maxc2nrmk;

/*              There is no need to apply the block reflector to the */
/*              residual of the matrix A stored in A(KB+1:M,KB+1:N), */
/*              since the submatrix contains NaN and we stop */
/*              the computation. */
/*              But, we need to apply the block reflector to the residual */
/*              right hand sides stored in A(KB+1:M,N+1:N+NRHS), if the */
/*              residual right hand sides exist.  This occurs */
/*              when ( NRHS != 0 AND KB <= (M-IOFFSET) ): */

/*              A(I+1:M,N+1:N+NRHS) := A(I+1:M,N+1:N+NRHS) - */
/*                               A(I+1:M,1:KB) * F(N+1:N+NRHS,1:KB)**H. */
		if (*nrhs > 0 && *kb < *m - *ioffset) {
		    i__1 = *m - if__;
		    z__1.r = -1., z__1.i = 0.;
		    zgemm_("No transpose", "Conjugate transpose", &i__1, nrhs,
			     kb, &z__1, &a[if__ + 1 + a_dim1], lda, &f[*n + 1 
			    + f_dim1], ldf, &c_b2, &a[if__ + 1 + (*n + 1) * 
			    a_dim1], lda);
		}

/*              There is no need to recompute the 2-norm of the */
/*              difficult columns, since we stop the factorization. */

/*              Array TAU(KF+1:MINMNFACT) is not set and contains */
/*              undefined elements. */

/*              Return from the routine. */

		return;
	    }

/*           Quick return, if the submatrix A(I:M,K:N) is */
/*           a zero matrix. We need to check it only if the column index */
/*           (same as row index) is larger than 1, since the condition */
/*           for the whole original matrix A_orig is checked in the main */
/*           routine. */

	    if (*maxc2nrmk == 0.) {

		*done = TRUE_;

/*              Set KB, the number of factorized partial columns */
/*                      that are non-zero in each step in the block, */
/*                      i.e. the rank of the factor R. */
/*              Set IF, the number of processed rows in the block, which */
/*                      is the same as the number of processed rows in */
/*                      the original whole matrix A_orig. */

		*kb = k - 1;
		if__ = i__ - 1;
		*relmaxc2nrmk = 0.;

/*              There is no need to apply the block reflector to the */
/*              residual of the matrix A stored in A(KB+1:M,KB+1:N), */
/*              since the submatrix is zero and we stop the computation. */
/*              But, we need to apply the block reflector to the residual */
/*              right hand sides stored in A(KB+1:M,N+1:N+NRHS), if the */
/*              residual right hand sides exist.  This occurs */
/*              when ( NRHS != 0 AND KB <= (M-IOFFSET) ): */

/*              A(I+1:M,N+1:N+NRHS) := A(I+1:M,N+1:N+NRHS) - */
/*                               A(I+1:M,1:KB) * F(N+1:N+NRHS,1:KB)**H. */

		if (*nrhs > 0 && *kb < *m - *ioffset) {
		    i__1 = *m - if__;
		    z__1.r = -1., z__1.i = 0.;
		    zgemm_("No transpose", "Conjugate transpose", &i__1, nrhs,
			     kb, &z__1, &a[if__ + 1 + a_dim1], lda, &f[*n + 1 
			    + f_dim1], ldf, &c_b2, &a[if__ + 1 + (*n + 1) * 
			    a_dim1], lda);
		}

/*              There is no need to recompute the 2-norm of the */
/*              difficult columns, since we stop the factorization. */

/*              Set TAUs corresponding to the columns that were not */
/*              factorized to ZERO, i.e. set TAU(KB+1:MINMNFACT) = CZERO, */
/*              which is equivalent to seting TAU(K:MINMNFACT) = CZERO. */

		i__1 = minmnfact;
		for (j = k; j <= i__1; ++j) {
		    i__2 = j;
		    tau[i__2].r = 0., tau[i__2].i = 0.;
		}

/*              Return from the routine. */

		return;

	    }

/*           ============================================================ */

/*           Check if the submatrix A(I:M,K:N) contains Inf, */
/*           set INFO parameter to the column number, where */
/*           the first Inf is found plus N, and continue */
/*           the computation. */
/*           We need to check the condition only if the */
/*           column index (same as row index) of the original whole */
/*           matrix is larger than 1, since the condition for whole */
/*           original matrix is checked in the main routine. */

	    if (*info == 0 && *maxc2nrmk > myhugeval) {
		*info = *n + k - 1 + kp;
	    }

/*           ============================================================ */

/*           Test for the second and third tolerance stopping criteria. */
/*           NOTE: There is no need to test for ABSTOL.GE.ZERO, since */
/*           MAXC2NRMK is non-negative. Similarly, there is no need */
/*           to test for RELTOL.GE.ZERO, since RELMAXC2NRMK is */
/*           non-negative. */
/*           We need to check the condition only if the */
/*           column index (same as row index) of the original whole */
/*           matrix is larger than 1, since the condition for whole */
/*           original matrix is checked in the main routine. */

	    *relmaxc2nrmk = *maxc2nrmk / *maxc2nrm;

	    if (*maxc2nrmk <= *abstol || *relmaxc2nrmk <= *reltol) {

		*done = TRUE_;

/*              Set KB, the number of factorized partial columns */
/*                      that are non-zero in each step in the block, */
/*                      i.e. the rank of the factor R. */
/*              Set IF, the number of processed rows in the block, which */
/*                      is the same as the number of processed rows in */
/*                      the original whole matrix A_orig; */

		*kb = k - 1;
		if__ = i__ - 1;

/*              Apply the block reflector to the residual of the */
/*              matrix A and the residual of the right hand sides B, if */
/*              the residual matrix and and/or the residual of the right */
/*              hand sides exist,  i.e. if the submatrix */
/*              A(I+1:M,KB+1:N+NRHS) exists.  This occurs when */
/*                 KB < MINMNUPDT = f2cmin( M-IOFFSET, N+NRHS ): */

/*              A(IF+1:M,K+1:N+NRHS) := A(IF+1:M,KB+1:N+NRHS) - */
/*                             A(IF+1:M,1:KB) * F(KB+1:N+NRHS,1:KB)**H. */

		if (*kb < minmnupdt) {
		    i__1 = *m - if__;
		    i__2 = *n + *nrhs - *kb;
		    z__1.r = -1., z__1.i = 0.;
		    zgemm_("No transpose", "Conjugate transpose", &i__1, &
			    i__2, kb, &z__1, &a[if__ + 1 + a_dim1], lda, &f[*
			    kb + 1 + f_dim1], ldf, &c_b2, &a[if__ + 1 + (*kb 
			    + 1) * a_dim1], lda);
		}

/*              There is no need to recompute the 2-norm of the */
/*              difficult columns, since we stop the factorization. */

/*              Set TAUs corresponding to the columns that were not */
/*              factorized to ZERO, i.e. set TAU(KB+1:MINMNFACT) = CZERO, */
/*              which is equivalent to seting TAU(K:MINMNFACT) = CZERO. */

		i__1 = minmnfact;
		for (j = k; j <= i__1; ++j) {
		    i__2 = j;
		    tau[i__2].r = 0., tau[i__2].i = 0.;
		}

/*              Return from the routine. */

		return;

	    }

/*           ============================================================ */

/*           End ELSE of IF(I.EQ.1) */

	}

/*        =============================================================== */

/*        If the pivot column is not the first column of the */
/*        subblock A(1:M,K:N): */
/*        1) swap the K-th column and the KP-th pivot column */
/*           in A(1:M,1:N); */
/*        2) swap the K-th row and the KP-th row in F(1:N,1:K-1) */
/*        3) copy the K-th element into the KP-th element of the partial */
/*           and exact 2-norm vectors VN1 and VN2. (Swap is not needed */
/*           for VN1 and VN2 since we use the element with the index */
/*           larger than K in the next loop step.) */
/*        4) Save the pivot interchange with the indices relative to the */
/*           the original matrix A_orig, not the block A(1:M,1:N). */

	if (kp != k) {
	    zswap_(m, &a[kp * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
	    i__1 = k - 1;
	    zswap_(&i__1, &f[kp + f_dim1], ldf, &f[k + f_dim1], ldf);
	    vn1[kp] = vn1[k];
	    vn2[kp] = vn2[k];
	    itemp = jpiv[kp];
	    jpiv[kp] = jpiv[k];
	    jpiv[k] = itemp;
	}

/*        Apply previous Householder reflectors to column K: */
/*        A(I:M,K) := A(I:M,K) - A(I:M,1:K-1)*F(K,1:K-1)**H. */

	if (k > 1) {
	    i__1 = k - 1;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = k + j * f_dim1;
		d_cnjg(&z__1, &f[k + j * f_dim1]);
		f[i__2].r = z__1.r, f[i__2].i = z__1.i;
	    }
	    i__1 = *m - i__ + 1;
	    i__2 = k - 1;
	    z__1.r = -1., z__1.i = 0.;
	    zgemv_("No transpose", &i__1, &i__2, &z__1, &a[i__ + a_dim1], lda,
		     &f[k + f_dim1], ldf, &c_b2, &a[i__ + k * a_dim1], &c__1);
	    i__1 = k - 1;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = k + j * f_dim1;
		d_cnjg(&z__1, &f[k + j * f_dim1]);
		f[i__2].r = z__1.r, f[i__2].i = z__1.i;
	    }
	}

/*        Generate elementary reflector H(k) using the column A(I:M,K). */

	if (i__ < *m) {
	    i__1 = *m - i__ + 1;
	    zlarfg_(&i__1, &a[i__ + k * a_dim1], &a[i__ + 1 + k * a_dim1], &
		    c__1, &tau[k]);
	} else {
	    i__1 = k;
	    tau[i__1].r = 0., tau[i__1].i = 0.;
	}

/*        Check if TAU(K) contains NaN, set INFO parameter */
/*        to the column number where NaN is found and return from */
/*        the routine. */
/*        NOTE: There is no need to check TAU(K) for Inf, */
/*        since ZLARFG cannot produce TAU(KK) or Householder vector */
/*        below the diagonal containing Inf. Only BETA on the diagonal, */
/*        returned by ZLARFG can contain Inf, which requires */
/*        TAU(K) to contain NaN. Therefore, this case of generating Inf */
/*        by ZLARFG is covered by checking TAU(K) for NaN. */

	i__1 = k;
	d__1 = tau[i__1].r;
	if (disnan_(&d__1)) {
	    i__1 = k;
	    taunan = tau[i__1].r;
	} else /* if(complicated condition) */ {
	    d__1 = d_imag(&tau[k]);
	    if (disnan_(&d__1)) {
		taunan = d_imag(&tau[k]);
	    } else {
		taunan = 0.;
	    }
	}

	if (disnan_(&taunan)) {

	    *done = TRUE_;

/*           Set KB, the number of factorized partial columns */
/*                   that are non-zero in each step in the block, */
/*                   i.e. the rank of the factor R. */
/*           Set IF, the number of processed rows in the block, which */
/*                   is the same as the number of processed rows in */
/*                   the original whole matrix A_orig. */

	    *kb = k - 1;
	    if__ = i__ - 1;
	    *info = k;

/*           Set MAXC2NRMK and  RELMAXC2NRMK to NaN. */

	    *maxc2nrmk = taunan;
	    *relmaxc2nrmk = taunan;

/*           There is no need to apply the block reflector to the */
/*           residual of the matrix A stored in A(KB+1:M,KB+1:N), */
/*           since the submatrix contains NaN and we stop */
/*           the computation. */
/*           But, we need to apply the block reflector to the residual */
/*           right hand sides stored in A(KB+1:M,N+1:N+NRHS), if the */
/*           residual right hand sides exist.  This occurs */
/*           when ( NRHS != 0 AND KB <= (M-IOFFSET) ): */

/*           A(I+1:M,N+1:N+NRHS) := A(I+1:M,N+1:N+NRHS) - */
/*                            A(I+1:M,1:KB) * F(N+1:N+NRHS,1:KB)**H. */

	    if (*nrhs > 0 && *kb < *m - *ioffset) {
		i__1 = *m - if__;
		z__1.r = -1., z__1.i = 0.;
		zgemm_("No transpose", "Conjugate transpose", &i__1, nrhs, kb,
			 &z__1, &a[if__ + 1 + a_dim1], lda, &f[*n + 1 + 
			f_dim1], ldf, &c_b2, &a[if__ + 1 + (*n + 1) * a_dim1],
			 lda);
	    }

/*           There is no need to recompute the 2-norm of the */
/*           difficult columns, since we stop the factorization. */

/*           Array TAU(KF+1:MINMNFACT) is not set and contains */
/*           undefined elements. */

/*           Return from the routine. */

	    return;
	}

/*        =============================================================== */

	i__1 = i__ + k * a_dim1;
	aik.r = a[i__1].r, aik.i = a[i__1].i;
	i__1 = i__ + k * a_dim1;
	a[i__1].r = 1., a[i__1].i = 0.;

/*        =============================================================== */

/*        Compute the current K-th column of F: */
/*          1) F(K+1:N,K) := tau(K) * A(I:M,K+1:N)**H * A(I:M,K). */

	if (k < *n + *nrhs) {
	    i__1 = *m - i__ + 1;
	    i__2 = *n + *nrhs - k;
	    zgemv_("Conjugate transpose", &i__1, &i__2, &tau[k], &a[i__ + (k 
		    + 1) * a_dim1], lda, &a[i__ + k * a_dim1], &c__1, &c_b1, &
		    f[k + 1 + k * f_dim1], &c__1);
	}

/*           2) Zero out elements above and on the diagonal of the */
/*              column K in matrix F, i.e elements F(1:K,K). */

	i__1 = k;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = j + k * f_dim1;
	    f[i__2].r = 0., f[i__2].i = 0.;
	}

/*         3) Incremental updating of the K-th column of F: */
/*        F(1:N,K) := F(1:N,K) - tau(K) * F(1:N,1:K-1) * A(I:M,1:K-1)**H */
/*                    * A(I:M,K). */

	if (k > 1) {
	    i__1 = *m - i__ + 1;
	    i__2 = k - 1;
	    i__3 = k;
	    z__1.r = -tau[i__3].r, z__1.i = -tau[i__3].i;
	    zgemv_("Conjugate Transpose", &i__1, &i__2, &z__1, &a[i__ + 
		    a_dim1], lda, &a[i__ + k * a_dim1], &c__1, &c_b1, &auxv[1]
		    , &c__1);

	    i__1 = *n + *nrhs;
	    i__2 = k - 1;
	    zgemv_("No transpose", &i__1, &i__2, &c_b2, &f[f_dim1 + 1], ldf, &
		    auxv[1], &c__1, &c_b2, &f[k * f_dim1 + 1], &c__1);
	}

/*        =============================================================== */

/*        Update the current I-th row of A: */
/*        A(I,K+1:N+NRHS) := A(I,K+1:N+NRHS) */
/*                         - A(I,1:K)*F(K+1:N+NRHS,1:K)**H. */

	if (k < *n + *nrhs) {
	    i__1 = *n + *nrhs - k;
	    z__1.r = -1., z__1.i = 0.;
	    zgemm_("No transpose", "Conjugate transpose", &c__1, &i__1, &k, &
		    z__1, &a[i__ + a_dim1], lda, &f[k + 1 + f_dim1], ldf, &
		    c_b2, &a[i__ + (k + 1) * a_dim1], lda);
	}

	i__1 = i__ + k * a_dim1;
	a[i__1].r = aik.r, a[i__1].i = aik.i;

/*        Update the partial column 2-norms for the residual matrix, */
/*        only if the residual matrix A(I+1:M,K+1:N) exists, i.e. */
/*        when K < MINMNFACT = f2cmin( M-IOFFSET, N ). */

	if (k < minmnfact) {

	    i__1 = *n;
	    for (j = k + 1; j <= i__1; ++j) {
		if (vn1[j] != 0.) {

/*                 NOTE: The following lines follow from the analysis in */
/*                 Lapack Working Note 176. */

		    temp = z_abs(&a[i__ + j * a_dim1]) / vn1[j];
/* Computing MAX */
		    d__1 = 0., d__2 = (temp + 1.) * (1. - temp);
		    temp = f2cmax(d__1,d__2);
/* Computing 2nd power */
		    d__1 = vn1[j] / vn2[j];
		    temp2 = temp * (d__1 * d__1);
		    if (temp2 <= tol3z) {

/*                    At J-index, we have a difficult column for the */
/*                    update of the 2-norm. Save the index of the previous */
/*                    difficult column in IWORK(J-1). */
/*                    NOTE: ILSTCC > 1, threfore we can use IWORK only */
/*                    with N-1 elements, where the elements are */
/*                    shifted by 1 to the left. */

			iwork[j - 1] = lsticc;

/*                    Set the index of the last difficult column LSTICC. */

			lsticc = j;

		    } else {
			vn1[j] *= sqrt(temp);
		    }
		}
	    }

	}

/*        End of while loop. */

    }

/*     Now, afler the loop: */
/*        Set KB, the number of factorized columns in the block; */
/*        Set IF, the number of processed rows in the block, which */
/*                is the same as the number of processed rows in */
/*                the original whole matrix A_orig, IF = IOFFSET + KB. */

    *kb = k;
    if__ = i__;

/*     Apply the block reflector to the residual of the matrix A */
/*     and the residual of the right hand sides B, if the residual */
/*     matrix and and/or the residual of the right hand sides */
/*     exist,  i.e. if the submatrix A(I+1:M,KB+1:N+NRHS) exists. */
/*     This occurs when KB < MINMNUPDT = f2cmin( M-IOFFSET, N+NRHS ): */

/*     A(IF+1:M,K+1:N+NRHS) := A(IF+1:M,KB+1:N+NRHS) - */
/*                         A(IF+1:M,1:KB) * F(KB+1:N+NRHS,1:KB)**H. */

    if (*kb < minmnupdt) {
	i__1 = *m - if__;
	i__2 = *n + *nrhs - *kb;
	z__1.r = -1., z__1.i = 0.;
	zgemm_("No transpose", "Conjugate transpose", &i__1, &i__2, kb, &z__1,
		 &a[if__ + 1 + a_dim1], lda, &f[*kb + 1 + f_dim1], ldf, &c_b2,
		 &a[if__ + 1 + (*kb + 1) * a_dim1], lda);
    }

/*     Recompute the 2-norm of the difficult columns. */
/*     Loop over the index of the difficult columns from the largest */
/*     to the smallest index. */

    while(lsticc > 0) {

/*        LSTICC is the index of the last difficult column is greater */
/*        than 1. */
/*        ITEMP is the index of the previous difficult column. */

	itemp = iwork[lsticc - 1];

/*        Compute the 2-norm explicilty for the last difficult column and */
/*        save it in the partial and exact 2-norm vectors VN1 and VN2. */

/*        NOTE: The computation of VN1( LSTICC ) relies on the fact that */
/*        DZNRM2 does not fail on vectors with norm below the value of */
/*        SQRT(DLAMCH('S')) */

	i__1 = *m - if__;
	vn1[lsticc] = dznrm2_(&i__1, &a[if__ + 1 + lsticc * a_dim1], &c__1);
	vn2[lsticc] = vn1[lsticc];

/*        Downdate the index of the last difficult column to */
/*        the index of the previous difficult column. */

	lsticc = itemp;

    }

    return;

/*     End of ZLAQP3RK */

} /* zlaqp3rk_ */

