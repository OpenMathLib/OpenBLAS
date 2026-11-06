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

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZLARFT_LVL2 forms the triangular factor T of a block reflector H = I - vtvH */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLARFT_LVL2( DIRECT, STOREV, N, K, V, LDV, TAU, */
/*                    T, LDT ) */

/*       CHARACTER          DIRECT, STOREV */
/*       INTEGER            K, LDT, LDV, N */
/*       COMPLEX*16         T( LDT, * ), TAU( * ), V( LDV, * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLARFT_LVL2 forms the triangular factor T of a complex block reflector H */
/* > of order n, which is defined as a product of k elementary reflectors. */
/* > */
/* > If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular; */
/* > */
/* > If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular. */
/* > */
/* > If STOREV = 'C', the vector which defines the elementary reflector */
/* > H(i) is stored in the i-th column of the array V, and */
/* > */
/* >    H  =  I - V * T * V**H */
/* > */
/* > If STOREV = 'R', the vector which defines the elementary reflector */
/* > H(i) is stored in the i-th row of the array V, and */
/* > */
/* >    H  =  I - V**H * T * V */
/* > */
/* > If DIRECT or STOREV = 'T', see Further Details for the shape of T */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] DIRECT */
/* > \verbatim */
/* >          DIRECT is CHARACTER*1 */
/* >          Specifies the order in which the elementary reflectors are */
/* >          multiplied to form the block reflector: */
/* >          = 'F': H = H(1) H(2) . . . H(k) (Forward) */
/* >          = 'B': H = H(k) . . . H(2) H(1) (Backward) */
/* >          = 'T' (with STOREV='R'): H = H(k) . . . H(2) H(1) (backward) */
/* >                but we return the T matrix that is already (conjugate) transposed */
/* > \endverbatim */
/* > */
/* > \param[in] STOREV */
/* > \verbatim */
/* >          STOREV is CHARACTER*1 */
/* >          Specifies how the vectors which define the elementary */
/* >          reflectors are stored (see also Further Details): */
/* >          = 'C': column-wise */
/* >          = 'R': row-wise */
/* >          = 'T': (With DIRECT='F') Row-wise, but we return the T */
/* >                matrix that is already (conjugate) transposed. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the block reflector H. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The order of the triangular factor T (= the number of */
/* >          elementary reflectors). K >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* >          V is COMPLEX*16 array, dimension */
/* >                               (LDV,K) if STOREV = 'C' */
/* >                               (LDV,N) if STOREV = 'R' */
/* >          The matrix V. See further details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the array V. */
/* >          If STOREV = 'C', LDV >= f2cmax(1,N); if STOREV = 'R', LDV >= K. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX*16 array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          T is COMPLEX*16 array, dimension (LDT,K) */
/* >          The k by k triangular factor T of the block reflector. */
/* >          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is */
/* >          lower triangular. The rest of the array is not used. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T. LDT >= K. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \ingroup larft */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The shape of the matrix V and the storage of the vectors which define */
/* >  the H(i) is best illustrated by the following example with n = 5 and */
/* >  k = 3. The elements equal to 1 are not stored. */
/* > */
/* >  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R': */
/* > */
/* >               V = (  1       )                 V = (  1 v1 v1 v1 v1 ) */
/* >                   ( v1  1    )                     (     1 v2 v2 v2 ) */
/* >                   ( v1 v2  1 )                     (        1 v3 v3 ) */
/* >                   ( v1 v2 v3 ) */
/* >                   ( v1 v2 v3 ) */
/* > */
/* >  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R': */
/* > */
/* >               V = ( v1 v2 v3 )                 V = ( v1 v1  1       ) */
/* >                   ( v1 v2 v3 )                     ( v2 v2 v2  1    ) */
/* >                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 ) */
/* >                   (     1 v3 ) */
/* >                   (        1 ) */
/* > */
/* >  In addition, the shape of T is determined by these same flags as */
/* >   in the below table. */
/* >    'U' denotes upper triangular */
/* >    'L' denotes lower triangular */
/* >    'X' denotes no current implementation */
/* >    We also provide the logical variable that represents the case */
/* >    in the code if it is implemented */
/* > */
/* >    |-----------------------------------------------------------| */
/* >    |              | DIRECT = 'F' | DIRECT = 'B' | DIRECT = 'T' | */
/* >    |--------------+--------------+--------------+--------------| */
/* >    | STOREV = 'C' | U (QR)       | L (QL)       | X            | */
/* >    | STOREV = 'R' | U (LQ)       | L (RQ)       | U (RQT)      | */
/* >    | STOREV = 'T' | L (LQT)      | X            | X            | */
/* >    |-----------------------------------------------------------| */
/* > */
/* >    Finally, the relationship between the (conjugate) transposed T matrices */
/* >    are as follows: (Note that T_{FC} denotes the T associated with calling */
/* >    this routine with DIRECT = 'F' and STOREV = 'C') */
/* > */
/* >    T_{TR} = (T_{BR})**H */
/* >    T_{FT} = (T_{FR})**H */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ void zlarft_lvl2__(char *direct, char *storev, integer *n, 
	integer *k, doublecomplex *v, integer *ldv, doublecomplex *tau, 
	doublecomplex *t, integer *ldt)
{
    /* System generated locals */
    integer t_dim1, t_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3;

    /* Local variables */
    integer i__, j;
    logical lq, ql, qr, rq;
    integer kmi, nmi;
    logical lqt, rqt, dirb, dirf;
    integer info;
    logical dirt;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ void zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *), zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *), 
	    ztrmv_(char *, char *, char *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), 
	    xerbla_(char *, integer *, ftnlen);
    logical storec, storer, storet;


/*  -- LAPACK auxiliary routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */













/*     Convert our character flags to logical flags for later */

    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --tau;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;

    /* Function Body */
    dirf = lsame_(direct, "F");
    dirb = lsame_(direct, "B");
    dirt = lsame_(direct, "T");
    storec = lsame_(storev, "C");
    storer = lsame_(storev, "R");
    storet = lsame_(storev, "T");

/*     Error handling for our character flags */

    info = 0;
    if (! (dirf || dirb || dirt)) {

/*        DIRECT holds an illegal value */

	info = 1;
    } else if (! (storec || storer || storet)) {

/*        STOREV holds an illegal value */

	info = 2;
    } else if (dirb && storet) {

/*        This case is purposefully not implemented, but any other value for */
/*        STOREV is valid, so we report STOREV as the invalid input */

	info = 2;
    } else if (dirt && storec) {

/*        This case is purposefully not implemented, but any other value for */
/*        DIRECT is valid, so we report DIRECT as the invalid input */

	info = 1;
    } else if (dirt && storet) {

/*        This case is purposefully not implemented, and is ambiguous what */
/*        the user wants to do, so we arbitrarily say DIRECT is the incorrect */
/*        character flag. */

	info = 1;
    }
    if (info != 0) {
	xerbla_("ZLARFT_LVL2", &info, (ftnlen)11);
	return;
    }

/*     Quick return if possible */

    if (*n == 0 || *k == 0) {
	return;
    }

/*     Now we determine what factorization our flags are associated with */

/*     QR happens when we have forward direction in column storage */

    qr = dirf && storec;

/*     LQT happens when we have forward direction in row storage and want to compute the transpose of */
/*     the T we would normally compute */

    lqt = dirf && storet;

/*     LQ happens when we have forward direction in row storage and want to compute the T we would */
/*     normally compute */

    lq = dirf && storer;

/*     QL happens when we have backward direction in column storage */

    ql = dirb && storec;

/*     RQT happens when we have backward direction in row storage and want to compute the transpose */
/*     of the T we would normally compute */

    rqt = dirt && storer;

/*     RQ happens when we have backward direction in row storage and want to compute the T that we */
/*     would normally compute */

    rq = dirb && storer;
    if (*n == 1 || *k == 1) {
	if (lqt || rqt) {
	    i__1 = t_dim1 + 1;
	    d_cnjg(&z__1, &tau[1]);
	    t[i__1].r = z__1.r, t[i__1].i = z__1.i;
	} else {
	    i__1 = t_dim1 + 1;
	    t[i__1].r = tau[1].r, t[i__1].i = tau[1].i;
	}
	return;
    }
    if (qr) {

/*        Break V into 9 components */

/*        V = |-----------------------| */
/*            |V_{1,1} 0       0      | i-1 */
/*            |V_{2,1} V_{2,2} 0      | 1 */
/*            |V_{3,1} V_{3,2} V_{3,3}| n-i */
/*            |-----------------------| */
/*             i-1     1       k-i */

/*        V_{1,1}, V_{2,2} and V_{3,3} are unit lower triangular */

/*        This is how we are going to view the matrix V at each step */
/*        i=2,\dots,k, then we grow into V_{3,3} and repeat until we */
/*        reach the end. On each iteration V_{3,3} is not referenced */

/*        We will construct T one column at a time from left to right */
/*        after initializing T(1,1) = TAU(1) */

/*        T = |-------------------------| */
/*            | T_{1,1} T_{1,2} T_{1,3} | i-1 */
/*            | 0       T_{2,2} T_{2,3} | 1 */
/*            | 0       0       T_{3,3} | k-i */
/*            |-------------------------| */
/*              i-1     1       k-i */

/*        T_{1,1}, T_{2,2}, and T_{3,3} are non-unit lower triangular */

/*        Similarly as above, we will construct T_{1,2} and T_{2,2} at */
/*        each iteration i = 2, \dots k, and then grow into T_{1:3,3}. On */
/*        each iteration, T_{1:3,3} are not referenced. See zlarft.f */
/*        for details on how these formulae were constructed. */

/*        We now get */

/*        T_{1,2} = -T_{1,1}[V_{1,1}\\V_{2,1}\\V_{3,1}]' */
/*                       [0\\V_{2,2}\\V_{3,2}]T_{2,2} */

/*        T_{1,2} = -T_{1,1}(V_{2,1}' + V_{3,1}'V_{3,2})T_{2,2} */

/*        This means we will do the following */

/*        T_{1,2} = -V_{2,1}'T_{2,2} = -\tau_{i}V_{2,1}' */
/*        T_{1,2} = -\tau_{i}V_{3,2}' V_{3,1} + T_{1,2} */
/*        T_{1,2} = T_{1,1}T_{1,2} */
/*        T_{2,2} = \tau{i} */

	i__1 = t_dim1 + 1;
	t[i__1].r = tau[1].r, t[i__1].i = tau[1].i;
	i__1 = *k;
	for (i__ = 2; i__ <= i__1; ++i__) {

/*           T_{1,2} = -V_{2,1}'V_{2,2}T_{2,2} = -\tau_i V_{2,1}' */
/*           We must do this at copy time as otherwise gemv will do nothing */
/*           on the last column when n=k, but we neet to make sure we are */
/*           scaled by this value */

	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = j + i__ * t_dim1;
		i__4 = i__;
		z__2.r = -tau[i__4].r, z__2.i = -tau[i__4].i;
		d_cnjg(&z__3, &v[i__ + j * v_dim1]);
		z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = z__2.r * 
			z__3.i + z__2.i * z__3.r;
		t[i__3].r = z__1.r, t[i__3].i = z__1.i;
	    }

/*           T_{1,2} = -V_{3,1}'V_{3,2}T_{2,2} + T_{1,2} */
/*                   = -\tau{i} V_{3,2}'V_{3,1} + T_{1,2} */

	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    i__4 = i__;
	    z__1.r = -tau[i__4].r, z__1.i = -tau[i__4].i;
	    zgemv_("Conjugate Transpose", &i__2, &i__3, &z__1, &v[i__ + 1 + 
		    v_dim1], ldv, &v[i__ + 1 + i__ * v_dim1], &c__1, &c_b1, &
		    t[i__ * t_dim1 + 1], &c__1);

/*           T_{1,2} = T_{1,1}T_{1,2} */

	    i__2 = i__ - 1;
	    ztrmv_("Upper", "No Transpose", "Non-unit", &i__2, &t[t_offset], 
		    ldt, &t[i__ * t_dim1 + 1], &c__1);

/*           T_{2,2} = \tau{i} */

	    i__2 = i__ + i__ * t_dim1;
	    i__3 = i__;
	    t[i__2].r = tau[i__3].r, t[i__2].i = tau[i__3].i;
	}
    } else if (lq) {

/*        Break V into 9 components */

/*        V = |-------------------------| */
/*            | V_{1,1} V_{1,2} V_{1,3} | i-1 */
/*            | 0       V_{2,2} V_{2,3} | 1 */
/*            | 0       0       V_{3,3} | k-i */
/*            |-------------------------| */
/*              i-1     1       n-i */

/*        V_{1,1}, V_{2,2} and V_{3,3} are unit upper triangular */

/*        This is how we are going to view the matrix V at each step */
/*        i=2,\dots,k, then we grow into V_{3,3} and repeat until we */
/*        reach the end. On each iteration V_{3,3} is not referenced */

/*        We will construct T one column at a time from left to right */
/*        after initializing T(1,1) = TAU(1) */

/*        T = |-------------------------| */
/*            | T_{1,1} T_{1,2} T_{1,3} | i-1 */
/*            | 0       T_{2,2} T_{2,3} | 1 */
/*            | 0       0       T_{3,3} | k-i */
/*            |-------------------------| */
/*              i-1     1       k-i */

/*        Similarly as above, we will construct T_{1,2} and T_{2,2} at */
/*        each iteration i = 2, \dots k, and then grow into T_{1:3,3}. On */
/*        each iteration, T_{1:3,3} are not referenced. See zlarft.f */
/*        for details on how these formulae were constructed. */

/*        We now get */

/*        T_{1,2} = -T_{1,1}[V_{1,1} V_{1,2} V_{1,3}][ 0 V_{2,2} V_{2,3} ]'T_{2,2} */

/*        T_{1,2} = -T_{1,1}(V_{1,2} + V_{1,3}V_{2,3}')T_{2,2} */

/*        This means we will do the following */

/*        T_{1,2} = -V_{1,2}T_{2,2} = -\tau_{i}V_{1,2} */
/*        T_{1,2} = -\tau_{i}V_{1,3}V_{2,3}' + T_{1,2} */
/*        T_{1,2} = T_{1,1}T_{1,2} */
/*        T_{2,2} = \tau{i} */

	i__1 = t_dim1 + 1;
	t[i__1].r = tau[1].r, t[i__1].i = tau[1].i;
	i__1 = *k;
	for (i__ = 2; i__ <= i__1; ++i__) {

/*           T_{1,2} = -\tau_{i}V_{1,2} */

	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = j + i__ * t_dim1;
		i__4 = i__;
		z__2.r = -tau[i__4].r, z__2.i = -tau[i__4].i;
		i__5 = j + i__ * v_dim1;
		z__1.r = z__2.r * v[i__5].r - z__2.i * v[i__5].i, z__1.i = 
			z__2.r * v[i__5].i + z__2.i * v[i__5].r;
		t[i__3].r = z__1.r, t[i__3].i = z__1.i;
	    }

/*           T_{1,2} = -\tau_{i}V_{1,3}V_{2,3}' + T_{1,2} */

	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    i__4 = i__;
	    z__1.r = -tau[i__4].r, z__1.i = -tau[i__4].i;
	    zgemm_("No Transpose", "Conjugate Transpose", &i__2, &c__1, &i__3,
		     &z__1, &v[(i__ + 1) * v_dim1 + 1], ldv, &v[i__ + (i__ + 
		    1) * v_dim1], ldv, &c_b1, &t[i__ * t_dim1 + 1], ldt);

/*           T_{1,2} = T_{1,1}T_{1,2} */

	    i__2 = i__ - 1;
	    ztrmv_("Upper", "No Transpose", "Non-unit", &i__2, &t[t_offset], 
		    ldt, &t[i__ * t_dim1 + 1], &c__1);

/*           T_{2,2} = \tau{i} */

	    i__2 = i__ + i__ * t_dim1;
	    i__3 = i__;
	    t[i__2].r = tau[i__3].r, t[i__2].i = tau[i__3].i;
	}
    } else if (lqt) {

/*        Break V into 9 components */

/*        V = |-------------------------| */
/*            | V_{1,1} V_{1,2} V_{1,3} | i-1 */
/*            | 0       V_{2,2} V_{2,3} | 1 */
/*            | 0       0       V_{3,3} | k-i */
/*            |-------------------------| */
/*              i-1     1       n-i */

/*        V_{1,1}, V_{2,2} and V_{3,3} are unit upper triangular */

/*        This is how we are going to view the matrix V at each step */
/*        i=2,\dots,k, then we grow into V_{3,3} and repeat until we */
/*        reach the end. On each iteration V_{3,3} is not referenced */

/*        We will construct T one column at a time from left to right */
/*        after initializing T(1,1) = TAU(1) */

/*        T = |-------------------------| */
/*            | T_{1,1} 0       0       | i-1 */
/*            | T_{2,1} T_{2,2} 0       | 1 */
/*            | T_{3,1} T_{3,2} T_{3,3} | k-i */
/*            |-------------------------| */
/*              i-1     1       k-i */

/*        Similarly as above, we will construct T_{2,1} and T_{2,2} at */
/*        each iteration i = 2, \dots k, and then grow into T_{3,1:3}. On */
/*        each iteration, T_{3,1:3} are not referenced. See zlarft.f */
/*        for details on how these formulae were constructed. */

/*        We now get */

/*        T_{2,1} = -T_{2,2}[0 V_{2,2} V_{2,3}][V_{1,1} V_{1,2} V_{1,3}]'T_{1,1} */

/*        T_{2,1} = -T_{2,2}(V_{1,2}' + V_{2,3}V_{1,3}')T_{1,1} */

/*        This means we will do the following */

/*        T_{2,1} = -T_{2,2}V_{1,2}' = -\tau_{i}V_{1,2}' */
/*        T_{2,1} = -\tau_{i}V_{1,3}V_{2,3}' + T_{2,1} */
/*        T_{2,1} = T_{1,1}'T_{2,1} */
/*        T_{2,2} = \tau{i} */

	i__1 = t_dim1 + 1;
	d_cnjg(&z__1, &tau[1]);
	t[i__1].r = z__1.r, t[i__1].i = z__1.i;
	i__1 = *k;
	for (i__ = 2; i__ <= i__1; ++i__) {

/*           T_{2,1} = -\tau_{i}V_{1,2}' */

	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = i__ + j * t_dim1;
		i__4 = i__;
		i__5 = j + i__ * v_dim1;
		z__3.r = tau[i__4].r * v[i__5].r - tau[i__4].i * v[i__5].i, 
			z__3.i = tau[i__4].r * v[i__5].i + tau[i__4].i * v[
			i__5].r;
		d_cnjg(&z__2, &z__3);
		z__1.r = -z__2.r, z__1.i = -z__2.i;
		t[i__3].r = z__1.r, t[i__3].i = z__1.i;
	    }

/*           T_{2,1} = -\tau_{i}V_{2,3}V_{1,3}' + T_{2,1} */

	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    d_cnjg(&z__2, &tau[i__]);
	    z__1.r = -z__2.r, z__1.i = -z__2.i;
	    zgemm_("No Transpose", "Conjugate Transpose", &c__1, &i__2, &i__3,
		     &z__1, &v[i__ + (i__ + 1) * v_dim1], ldv, &v[(i__ + 1) * 
		    v_dim1 + 1], ldv, &c_b1, &t[i__ + t_dim1], ldt);

/*           T_{2,1} = T_{1,1}'T_{2,1} */

	    i__2 = i__ - 1;
	    ztrmv_("Lower", "Transpose", "Non-unit", &i__2, &t[t_offset], ldt,
		     &t[i__ + t_dim1], ldt);
	    i__2 = i__ + i__ * t_dim1;
	    d_cnjg(&z__1, &tau[i__]);
	    t[i__2].r = z__1.r, t[i__2].i = z__1.i;
	}
    } else if (ql) {

/*     Break V into 9 components */

/*     V = |-------------------------| */
/*         | V_{1,1} V_{1,2} V_{1,3} | n-i */
/*         | 0       V_{2,2} V_{2,3} | 1 */
/*         | 0       0       V_{3,3} | i-1 */
/*         |-------------------------| */
/*           k-i     1       i-1 */

/*        V_{1,1}, V_{2,2} and V_{3,3} are unit upper triangular */

/*        This is how we are going to view the matrix V at each step */
/*        i=2,\dots,k, then we grow into V_{1,1} and repeat until we */
/*        reach the end. On each iteration V_{1,1} is not referenced */

/*        We will construct T one column at a time from right to left */
/*        after initializing T(K,K) = TAU(K) */

/*     T = |-------------------------| */
/*         | T_{1,1} 0       0       | k-i */
/*         | T_{2,1} T_{2,2} 0       | 1 */
/*         | T_{3,1} T_{3,2} T_{3,3} | i-1 */
/*         |-------------------------| */
/*           k-i     1       i-1 */

/*        T_{1,1}, T_{2,2}, and T_{3,3} are non-unit lower triangular */

/*        Similarly as above, we will construct T_{2,2} and T_{3,2} at */
/*        each iteration i = 2, \dots k, and then grow into T_{1:3,1}. On */
/*        each iteration, T_{1:3,1} are not referenced. See zlarft.f */
/*        for details on how these formulae were constructed. */

/*        We get that */

/*        T_{3,2} = -T_{3,3}[V_{1,3}\\V_{2,3}\\V_{3,3}]' */
/*        [V_{1,2}\\V_{2,2}\\0]T_{2,2} */

/*        T_{3,2} = -T_{3,3}(V_{1,3}'V_{1,2} + V_{2,3}')T_{2,2} */

/*        Thus, we will compute */

/*        T_{2,2} = \tau_{k-i+1} */
/*        T_{3,2} = -T_{3,3}V_{3,2}' = -\tau_{k-i+1}V_{3,2}' */
/*        T_{3,2} = -\tau_{k-i+1}V_{1,3}'V_{1,2} + T_{3,2} */
/*        T_{3,2} = T_{3,3}T_{3,2} */

	i__1 = *k + *k * t_dim1;
	i__2 = *k;
	t[i__1].r = tau[i__2].r, t[i__1].i = tau[i__2].i;
	i__1 = *k;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    kmi = *k - i__ + 1;
	    nmi = *n - i__ + 1;

/*             T_{2,2} = \tau(k-i+1) */

	    i__2 = kmi + kmi * t_dim1;
	    i__3 = kmi;
	    t[i__2].r = tau[i__3].r, t[i__2].i = tau[i__3].i;

/*             T_{3,2} = -\tau(k-i+1)V_{2,3}' */

	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = kmi + j + kmi * t_dim1;
		i__4 = kmi;
		z__2.r = -tau[i__4].r, z__2.i = -tau[i__4].i;
		d_cnjg(&z__3, &v[nmi + (kmi + j) * v_dim1]);
		z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = z__2.r * 
			z__3.i + z__2.i * z__3.r;
		t[i__3].r = z__1.r, t[i__3].i = z__1.i;
	    }

/*             T_{3,2} = -\tau(k-i+1)V_{1,3}'V_{1,2} + T_{3,2} */

	    i__2 = *n - i__;
	    i__3 = i__ - 1;
	    i__4 = kmi;
	    z__1.r = -tau[i__4].r, z__1.i = -tau[i__4].i;
	    zgemv_("Conjugate Transpose", &i__2, &i__3, &z__1, &v[(kmi + 1) * 
		    v_dim1 + 1], ldv, &v[kmi * v_dim1 + 1], &c__1, &c_b1, &t[
		    kmi + 1 + kmi * t_dim1], &c__1);

/*             T_{3,2} = T_{3,3}T_{3,2} */

	    i__2 = i__ - 1;
	    ztrmv_("Lower", "No Transpose", "Non-unit", &i__2, &t[kmi + 1 + (
		    kmi + 1) * t_dim1], ldt, &t[kmi + 1 + kmi * t_dim1], &
		    c__1);
	}
    } else if (rq) {

/*     Break V into 9 components */

/*     V = |-------------------------| */
/*         | V_{1,1} 0       0       | k-i */
/*         | V_{2,1} V_{2,2} 0       | 1 */
/*         | V_{3,1} V_{3,2} V_{3,3} | i-1 */
/*         |-------------------------| */
/*           n-i     1       i-1 */

/*        V_{1,1}, V_{2,2} and V_{3,3} are unit lower triangular */

/*        This is how we are going to view the matrix V at each step */
/*        i=2,\dots,k, then we grow into V_{1,1} and repeat until we */
/*        reach the end. On each iteration V_{1,1} is not referenced */

/*        We will construct T one column at a time from right to left */
/*        after initializing T(K,K) = TAU(K) */

/*     T = |-------------------------| */
/*         | T_{1,1} 0       0       | k-i */
/*         | T_{2,1} T_{2,2} 0       | 1 */
/*         | T_{3,1} T_{3,2} T_{3,3} | i-1 */
/*         |-------------------------| */
/*           k-i     1       i-1 */

/*        T_{1,1}, T_{2,2}, and T_{3,3} are non-unit lower triangular */

/*        Similarly as above, we will construct T_{2,2} and T_{3,2} at */
/*        each iteration i = 2, \dots k, and then grow into T_{1:3,1}. On */
/*        each iteration, T_{1:3,1} are not referenced. See zlarft.f */
/*        for details on how these formulae were constructed. */

/*        We get that */

/*        T_{3,2} = -T_{3,3}[V_{3,1} V_{3,2} V_{3,3}][V_{2,1} V_{2,2} 0]'T_{2,2} */

/*        T_{3,2} = -T_{3,3}(V_{3,1}V_{2,1}' + V_{3,2})T_{2,2} */

/*        Thus, we will compute */

/*        T_{2,2} = \tau_{k-i+1} */
/*        T_{3,2} = -\tau_{k-i+1}V_{3,2} */
/*        T_{3,2} = -\tau_{k-i+1}V_{3,1}V_{2,1}' + T_{3,2} */
/*        T_{3,2} = T_{3,3}T_{3,2} */

	i__1 = *k + *k * t_dim1;
	i__2 = *k;
	t[i__1].r = tau[i__2].r, t[i__1].i = tau[i__2].i;
	i__1 = *k;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    kmi = *k - i__ + 1;
	    nmi = *n - i__ + 1;

/*           T_{2,2} = \tau_{k-i+1} */

	    i__2 = kmi + kmi * t_dim1;
	    i__3 = kmi;
	    t[i__2].r = tau[i__3].r, t[i__2].i = tau[i__3].i;

/*           T_{3,2} = -\tau_{k-i+1}V_{3,2} */

	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = kmi + j + kmi * t_dim1;
		i__4 = kmi;
		z__2.r = -tau[i__4].r, z__2.i = -tau[i__4].i;
		i__5 = kmi + j + nmi * v_dim1;
		z__1.r = z__2.r * v[i__5].r - z__2.i * v[i__5].i, z__1.i = 
			z__2.r * v[i__5].i + z__2.i * v[i__5].r;
		t[i__3].r = z__1.r, t[i__3].i = z__1.i;
	    }

/*           T_{3,2} = -\tau_{k-i+1}V_{3,1}V_{2,1}' + T_{3,2} */

	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    i__4 = kmi;
	    z__1.r = -tau[i__4].r, z__1.i = -tau[i__4].i;
	    zgemm_("No Transpose", "Conjugate Transpose", &i__2, &c__1, &i__3,
		     &z__1, &v[kmi + 1 + v_dim1], ldv, &v[kmi + v_dim1], ldv, 
		    &c_b1, &t[kmi + 1 + kmi * t_dim1], ldt);

/*           T_{3,2} = T_{3,3}T_{3,2} */

	    i__2 = i__ - 1;
	    ztrmv_("Lower", "No Transpose", "Non-unit", &i__2, &t[kmi + 1 + (
		    kmi + 1) * t_dim1], ldt, &t[kmi + 1 + kmi * t_dim1], &
		    c__1);
	}
    } else if (rqt) {

/*     Break V into 9 components */

/*     V = |-------------------------| */
/*         | V_{1,1} 0       0       | k-i */
/*         | V_{2,1} V_{2,2} 0       | 1 */
/*         | V_{3,1} V_{3,2} V_{3,3} | i-1 */
/*         |-------------------------| */
/*           n-i     1       i-1 */

/*        V_{1,1}, V_{2,2} and V_{3,3} are unit lower triangular */

/*        This is how we are going to view the matrix V at each step */
/*        i=2,\dots,k, then we grow into V_{1,1} and repeat until we */
/*        reach the end. On each iteration V_{1,1} is not referenced */

/*        We will construct T one column at a time from right to left */
/*        after initializing T(K,K) = TAU(K) */

/*     T = |-------------------------| */
/*         | T_{1,1} T_{1,2} T_{1,3} | k-i */
/*         | 0       T_{2,2} T_{2,3} | 1 */
/*         | 0       0       T_{3,3} | i-1 */
/*         |-------------------------| */
/*           k-i     1       i-1 */

/*        T_{1,1}, T_{2,2}, and T_{3,3} are non-unit lower triangular */

/*        Similarly as above, we will construct T_{2,2} and T_{2,3} at */
/*        each iteration i = 2, \dots k, and then grow into T_{1,1:3}. On */
/*        each iteration, T_{1,1:3} are not referenced. See zlarft.f */
/*        for details on how these formulae were constructed. */

/*        We get that */

/*        T_{2,3} = -T_{2,2}[V_{2,1} V_{2,2} 0][V_{3,1} V_{3,2} V_{3,3}]'T_{3,3} */

/*        T_{3,2} = -T_{2,2}(V_{2,1}V_{3,1}' + V_{3,2}')T_{3,3} */

/*        Thus, we will compute */

/*        T_{2,2} = \tau_{k-i+1} */
/*        T_{2,3} = -\tau_{k-i+1}V_{3,2}' */
/*        T_{2,3} = -\tau_{k-i+1}V_{2,1}V_{3,1}' + T_{2,3} */
/*        T_{2,3} = T_{3,3}'T_{2,3} */

	i__1 = *k + *k * t_dim1;
	d_cnjg(&z__1, &tau[*k]);
	t[i__1].r = z__1.r, t[i__1].i = z__1.i;
	i__1 = *k;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    kmi = *k - i__ + 1;
	    nmi = *n - i__ + 1;

/*           T_{2,2} = \tau_{k-i+1} */

	    i__2 = kmi + kmi * t_dim1;
	    d_cnjg(&z__1, &tau[kmi]);
	    t[i__2].r = z__1.r, t[i__2].i = z__1.i;

/*           T_{2,3} = -\tau_{k-i+1}V_{3,2}' */

	    i__2 = i__ - 1;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = kmi + (kmi + j) * t_dim1;
		i__4 = kmi;
		i__5 = kmi + j + nmi * v_dim1;
		z__3.r = tau[i__4].r * v[i__5].r - tau[i__4].i * v[i__5].i, 
			z__3.i = tau[i__4].r * v[i__5].i + tau[i__4].i * v[
			i__5].r;
		d_cnjg(&z__2, &z__3);
		z__1.r = -z__2.r, z__1.i = -z__2.i;
		t[i__3].r = z__1.r, t[i__3].i = z__1.i;
	    }

/*           T_{2,3} = -\tau_{k-i+1}V_{2,1}V_{3,1}' + T_{2,3} */

	    i__2 = i__ - 1;
	    i__3 = *n - i__;
	    d_cnjg(&z__2, &tau[kmi]);
	    z__1.r = -z__2.r, z__1.i = -z__2.i;
	    zgemm_("No Transpose", "Conjugate Transpose", &c__1, &i__2, &i__3,
		     &z__1, &v[kmi + v_dim1], ldv, &v[kmi + 1 + v_dim1], ldv, 
		    &c_b1, &t[kmi + (kmi + 1) * t_dim1], ldt);

/*           T_{2,3} = T_{3,3}'T_{2,3} */

	    i__2 = i__ - 1;
	    ztrmv_("Upper", "Transpose", "Non-unit", &i__2, &t[kmi + 1 + (kmi 
		    + 1) * t_dim1], ldt, &t[kmi + (kmi + 1) * t_dim1], ldt);
	}
    }
    return;
} /* zlarft_lvl2__ */

