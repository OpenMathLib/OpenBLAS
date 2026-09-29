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

static complex c_b1 = {1.f,0.f};
static complex c_b3 = {-1.f,0.f};
static integer c__3 = 3;
static integer c_n1 = -1;
static integer c__2 = 2;

/* > \brief \b CLARFT forms the triangular factor T of a block reflector H = I - vtvH */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > Download CLARFT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarft.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarft.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarft.
f"> */
/* > [TXT]</a> */

/*  Definition: */
/*  =========== */

/*        SUBROUTINE CLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT ) */

/*       CHARACTER          DIRECT, STOREV */
/*       INTEGER            K, LDT, LDV, N */
/*       COMPLEX            T( LDT, * ), TAU( * ), V( LDV, * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARFT forms the triangular factor T of a complex block reflector H */
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
/* > \endverbatim */
/* > */
/* > \param[in] STOREV */
/* > \verbatim */
/* >          STOREV is CHARACTER*1 */
/* >          Specifies how the vectors which define the elementary */
/* >          reflectors are stored (see also Further Details): */
/* >          = 'C': columnwise */
/* >          = 'R': rowwise */
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
/* >          V is COMPLEX array, dimension */
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
/* >          TAU is COMPLEX array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          T is COMPLEX array, dimension (LDT,K) */
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
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ void clarft_(char *direct, char *storev, integer *n, integer *
	k, complex *v, integer *ldv, complex *tau, complex *t, integer *ldt)
{
    /* System generated locals */
    address a__1[2];
    integer t_dim1, t_offset, v_dim1, v_offset, i__1, i__2[2], i__3, i__4;
    complex q__1;
    char ch__1[2];

    /* Local variables */
    integer i__, j, l;
    logical lq, ql, qr;
    integer nx;
    extern /* Subroutine */ void clarft_lvl2__(char *, char *, integer *, 
	    integer *, complex *, integer *, complex *, complex *, integer *);
    logical dirf, colv;
    extern /* Subroutine */ void cgemm_(char *, char *, integer *, integer *, 
	    integer *, complex *, complex *, integer *, complex *, integer *, 
	    complex *, complex *, integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ void ctrmm_(char *, char *, char *, char *, 
	    integer *, integer *, complex *, complex *, integer *, complex *, 
	    integer *), clacpy_(char *, 
	    integer *, integer *, complex *, integer *, complex *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);


/*  -- LAPACK auxiliary routine -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */














/*     The general scheme used is inspired by the approach inside DGEQRT3 */
/*     which was (at the time of writing this code): */
/*     Based on the algorithm of Elmroth and Gustavson, */
/*     IBM J. Res. Develop. Vol 44 No. 4 July 2000. */

/*     Quick return if possible */

    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --tau;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;

    /* Function Body */
    if (*n == 0 || *k == 0) {
	return;
    }

/*     Base case */

    if (*n == 1 || *k == 1) {
	i__1 = t_dim1 + 1;
	t[i__1].r = tau[1].r, t[i__1].i = tau[1].i;
	return;
    }

/*     Determine when to cross over into the level 2 based implementation */

/* Writing concatenation */
    i__2[0] = 1, a__1[0] = direct;
    i__2[1] = 1, a__1[1] = storev;
    s_cat(ch__1, a__1, i__2, &c__2, (ftnlen)2);
    nx = ilaenv_(&c__3, "CLARFT", ch__1, n, k, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)2);
    if (*k < nx) {
	clarft_lvl2__(direct, storev, n, k, &v[v_offset], ldv, &tau[1], &t[
		t_offset], ldt);
	return;
    }

/*     Beginning of executable statements */

    l = *k / 2;

/*     Determine what kind of Q we need to compute */
/*     We assume that if the user doesn't provide 'F' for DIRECT, */
/*     then they meant to provide 'B' and if they don't provide */
/*     'C' for STOREV, then they meant to provide 'R' */

    dirf = lsame_(direct, "F");
    colv = lsame_(storev, "C");

/*     QR happens when we have forward direction in column storage */

    qr = dirf && colv;

/*     LQ happens when we have forward direction in row storage */

    lq = dirf && ! colv;

/*     QL happens when we have backward direction in column storage */

    ql = ! dirf && colv;

/*     The last case is RQ. Due to how we structured this, if the */
/*     above 3 are false, then RQ must be true, so we never store */
/*     this */
/*     RQ happens when we have backward direction in row storage */
/*     RQ = (.NOT.DIRF).AND.(.NOT.COLV) */

    if (qr) {

/*        Break V apart into 6 components */

/*        V = |---------------| */
/*            |V_{1,1} 0      | */
/*            |V_{2,1} V_{2,2}| */
/*            |V_{3,1} V_{3,2}| */
/*            |---------------| */

/*        V_{1,1}\in\C^{l,l}      unit lower triangular */
/*        V_{2,1}\in\C^{k-l,l}    rectangular */
/*        V_{3,1}\in\C^{n-k,l}    rectangular */

/*        V_{2,2}\in\C^{k-l,k-l}  unit lower triangular */
/*        V_{3,2}\in\C^{n-k,k-l}  rectangular */

/*        We will construct the T matrix */
/*        T = |---------------| */
/*            |T_{1,1} T_{1,2}| */
/*            |0       T_{2,2}| */
/*            |---------------| */

/*        T is the triangular factor obtained from block reflectors. */
/*        To motivate the structure, assume we have already computed T_{1,1} */
/*        and T_{2,2}. Then collect the associated reflectors in V_1 and V_2 */

/*        T_{1,1}\in\C^{l, l}     upper triangular */
/*        T_{2,2}\in\C^{k-l, k-l} upper triangular */
/*        T_{1,2}\in\C^{l, k-l}   rectangular */

/*        Where l = floor(k/2) */

/*        Then, consider the product: */

/*        (I - V_1*T_{1,1}*V_1')*(I - V_2*T_{2,2}*V_2') */
/*        = I - V_1*T_{1,1}*V_1' - V_2*T_{2,2}*V_2' + V_1*T_{1,1}*V_1'*V_2*T_{2,2}*V_2' */

/*        Define T{1,2} = -T_{1,1}*V_1'*V_2*T_{2,2} */

/*        Then, we can define the matrix V as */
/*        V = |-------| */
/*            |V_1 V_2| */
/*            |-------| */

/*        So, our product is equivalent to the matrix product */
/*        I - V*T*V' */
/*        This means, we can compute T_{1,1} and T_{2,2}, then use this information */
/*        to compute T_{1,2} */

/*        Compute T_{1,1} recursively */

	clarft_(direct, storev, n, &l, &v[v_offset], ldv, &tau[1], &t[
		t_offset], ldt);

/*        Compute T_{2,2} recursively */

	i__1 = *n - l;
	i__3 = *k - l;
	clarft_(direct, storev, &i__1, &i__3, &v[l + 1 + (l + 1) * v_dim1], 
		ldv, &tau[l + 1], &t[l + 1 + (l + 1) * t_dim1], ldt);

/*        Compute T_{1,2} */
/*        T_{1,2} = V_{2,1}' */

	i__1 = l;
	for (j = 1; j <= i__1; ++j) {
	    i__3 = *k - l;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__4 = j + (l + i__) * t_dim1;
		r_cnjg(&q__1, &v[l + i__ + j * v_dim1]);
		t[i__4].r = q__1.r, t[i__4].i = q__1.i;
	    }
	}

/*        T_{1,2} = T_{1,2}*V_{2,2} */

	i__1 = *k - l;
	ctrmm_("Right", "Lower", "No transpose", "Unit", &l, &i__1, &c_b1, &v[
		l + 1 + (l + 1) * v_dim1], ldv, &t[(l + 1) * t_dim1 + 1], ldt);

/*        T_{1,2} = V_{3,1}'*V_{3,2} + T_{1,2} */
/*        Note: We assume K <= N, and GEMM will do nothing if N=K */

	i__1 = *k - l;
	i__3 = *n - *k;
	cgemm_("Conjugate", "No transpose", &l, &i__1, &i__3, &c_b1, &v[*k + 
		1 + v_dim1], ldv, &v[*k + 1 + (l + 1) * v_dim1], ldv, &c_b1, &
		t[(l + 1) * t_dim1 + 1], ldt);

/*        At this point, we have that T_{1,2} = V_1'*V_2 */
/*        All that is left is to pre and post multiply by -T_{1,1} and T_{2,2} */
/*        respectively. */

/*        T_{1,2} = -T_{1,1}*T_{1,2} */

	i__1 = *k - l;
	ctrmm_("Left", "Upper", "No transpose", "Non-unit", &l, &i__1, &c_b3, 
		&t[t_offset], ldt, &t[(l + 1) * t_dim1 + 1], ldt);

/*        T_{1,2} = T_{1,2}*T_{2,2} */

	i__1 = *k - l;
	ctrmm_("Right", "Upper", "No transpose", "Non-unit", &l, &i__1, &c_b1,
		 &t[l + 1 + (l + 1) * t_dim1], ldt, &t[(l + 1) * t_dim1 + 1], 
		ldt);
    } else if (lq) {

/*        Break V apart into 6 components */

/*        V = |----------------------| */
/*            |V_{1,1} V_{1,2} V{1,3}| */
/*            |0       V_{2,2} V{2,3}| */
/*            |----------------------| */

/*        V_{1,1}\in\C^{l,l}      unit upper triangular */
/*        V_{1,2}\in\C^{l,k-l}    rectangular */
/*        V_{1,3}\in\C^{l,n-k}    rectangular */

/*        V_{2,2}\in\C^{k-l,k-l}  unit upper triangular */
/*        V_{2,3}\in\C^{k-l,n-k}  rectangular */

/*        Where l = floor(k/2) */

/*        We will construct the T matrix */
/*        T = |---------------| */
/*            |T_{1,1} T_{1,2}| */
/*            |0       T_{2,2}| */
/*            |---------------| */

/*        T is the triangular factor obtained from block reflectors. */
/*        To motivate the structure, assume we have already computed T_{1,1} */
/*        and T_{2,2}. Then collect the associated reflectors in V_1 and V_2 */

/*        T_{1,1}\in\C^{l, l}     upper triangular */
/*        T_{2,2}\in\C^{k-l, k-l} upper triangular */
/*        T_{1,2}\in\C^{l, k-l}   rectangular */

/*        Then, consider the product: */

/*        (I - V_1'*T_{1,1}*V_1)*(I - V_2'*T_{2,2}*V_2) */
/*        = I - V_1'*T_{1,1}*V_1 - V_2'*T_{2,2}*V_2 + V_1'*T_{1,1}*V_1*V_2'*T_{2,2}*V_2 */

/*        Define T_{1,2} = -T_{1,1}*V_1*V_2'*T_{2,2} */

/*        Then, we can define the matrix V as */
/*        V = |---| */
/*            |V_1| */
/*            |V_2| */
/*            |---| */

/*        So, our product is equivalent to the matrix product */
/*        I - V'*T*V */
/*        This means, we can compute T_{1,1} and T_{2,2}, then use this information */
/*        to compute T_{1,2} */

/*        Compute T_{1,1} recursively */

	clarft_(direct, storev, n, &l, &v[v_offset], ldv, &tau[1], &t[
		t_offset], ldt);

/*        Compute T_{2,2} recursively */

	i__1 = *n - l;
	i__3 = *k - l;
	clarft_(direct, storev, &i__1, &i__3, &v[l + 1 + (l + 1) * v_dim1], 
		ldv, &tau[l + 1], &t[l + 1 + (l + 1) * t_dim1], ldt);

/*        Compute T_{1,2} */
/*        T_{1,2} = V_{1,2} */

	i__1 = *k - l;
	clacpy_("All", &l, &i__1, &v[(l + 1) * v_dim1 + 1], ldv, &t[(l + 1) * 
		t_dim1 + 1], ldt);

/*        T_{1,2} = T_{1,2}*V_{2,2}' */

	i__1 = *k - l;
	ctrmm_("Right", "Upper", "Conjugate", "Unit", &l, &i__1, &c_b1, &v[l 
		+ 1 + (l + 1) * v_dim1], ldv, &t[(l + 1) * t_dim1 + 1], ldt);

/*        T_{1,2} = V_{1,3}*V_{2,3}' + T_{1,2} */
/*        Note: We assume K <= N, and GEMM will do nothing if N=K */

	i__1 = *k - l;
	i__3 = *n - *k;
	cgemm_("No transpose", "Conjugate", &l, &i__1, &i__3, &c_b1, &v[(*k + 
		1) * v_dim1 + 1], ldv, &v[l + 1 + (*k + 1) * v_dim1], ldv, &
		c_b1, &t[(l + 1) * t_dim1 + 1], ldt);

/*        At this point, we have that T_{1,2} = V_1*V_2' */
/*        All that is left is to pre and post multiply by -T_{1,1} and T_{2,2} */
/*        respectively. */

/*        T_{1,2} = -T_{1,1}*T_{1,2} */

	i__1 = *k - l;
	ctrmm_("Left", "Upper", "No transpose", "Non-unit", &l, &i__1, &c_b3, 
		&t[t_offset], ldt, &t[(l + 1) * t_dim1 + 1], ldt);

/*        T_{1,2} = T_{1,2}*T_{2,2} */

	i__1 = *k - l;
	ctrmm_("Right", "Upper", "No transpose", "Non-unit", &l, &i__1, &c_b1,
		 &t[l + 1 + (l + 1) * t_dim1], ldt, &t[(l + 1) * t_dim1 + 1], 
		ldt);
    } else if (ql) {

/*        Break V apart into 6 components */

/*        V = |---------------| */
/*            |V_{1,1} V_{1,2}| */
/*            |V_{2,1} V_{2,2}| */
/*            |0       V_{3,2}| */
/*            |---------------| */

/*        V_{1,1}\in\C^{n-k,k-l}  rectangular */
/*        V_{2,1}\in\C^{k-l,k-l}  unit upper triangular */

/*        V_{1,2}\in\C^{n-k,l}    rectangular */
/*        V_{2,2}\in\C^{k-l,l}    rectangular */
/*        V_{3,2}\in\C^{l,l}      unit upper triangular */

/*        We will construct the T matrix */
/*        T = |---------------| */
/*            |T_{1,1} 0      | */
/*            |T_{2,1} T_{2,2}| */
/*            |---------------| */

/*        T is the triangular factor obtained from block reflectors. */
/*        To motivate the structure, assume we have already computed T_{1,1} */
/*        and T_{2,2}. Then collect the associated reflectors in V_1 and V_2 */

/*        T_{1,1}\in\C^{k-l, k-l} non-unit lower triangular */
/*        T_{2,2}\in\C^{l, l}     non-unit lower triangular */
/*        T_{2,1}\in\C^{k-l, l}   rectangular */

/*        Where l = floor(k/2) */

/*        Then, consider the product: */

/*        (I - V_2*T_{2,2}*V_2')*(I - V_1*T_{1,1}*V_1') */
/*        = I - V_2*T_{2,2}*V_2' - V_1*T_{1,1}*V_1' + V_2*T_{2,2}*V_2'*V_1*T_{1,1}*V_1' */

/*        Define T_{2,1} = -T_{2,2}*V_2'*V_1*T_{1,1} */

/*        Then, we can define the matrix V as */
/*        V = |-------| */
/*            |V_1 V_2| */
/*            |-------| */

/*        So, our product is equivalent to the matrix product */
/*        I - V*T*V' */
/*        This means, we can compute T_{1,1} and T_{2,2}, then use this information */
/*        to compute T_{2,1} */

/*        Compute T_{1,1} recursively */

	i__1 = *n - l;
	i__3 = *k - l;
	clarft_(direct, storev, &i__1, &i__3, &v[v_offset], ldv, &tau[1], &t[
		t_offset], ldt);

/*        Compute T_{2,2} recursively */

	clarft_(direct, storev, n, &l, &v[(*k - l + 1) * v_dim1 + 1], ldv, &
		tau[*k - l + 1], &t[*k - l + 1 + (*k - l + 1) * t_dim1], ldt);

/*        Compute T_{2,1} */
/*        T_{2,1} = V_{2,2}' */

	i__1 = *k - l;
	for (j = 1; j <= i__1; ++j) {
	    i__3 = l;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__4 = *k - l + i__ + j * t_dim1;
		r_cnjg(&q__1, &v[*n - *k + j + (*k - l + i__) * v_dim1]);
		t[i__4].r = q__1.r, t[i__4].i = q__1.i;
	    }
	}

/*        T_{2,1} = T_{2,1}*V_{2,1} */

	i__1 = *k - l;
	ctrmm_("Right", "Upper", "No transpose", "Unit", &l, &i__1, &c_b1, &v[
		*n - *k + 1 + v_dim1], ldv, &t[*k - l + 1 + t_dim1], ldt);

/*        T_{2,1} = V_{2,2}'*V_{2,1} + T_{2,1} */
/*        Note: We assume K <= N, and GEMM will do nothing if N=K */

	i__1 = *k - l;
	i__3 = *n - *k;
	cgemm_("Conjugate", "No transpose", &l, &i__1, &i__3, &c_b1, &v[(*k - 
		l + 1) * v_dim1 + 1], ldv, &v[v_offset], ldv, &c_b1, &t[*k - 
		l + 1 + t_dim1], ldt);

/*        At this point, we have that T_{2,1} = V_2'*V_1 */
/*        All that is left is to pre and post multiply by -T_{2,2} and T_{1,1} */
/*        respectively. */

/*        T_{2,1} = -T_{2,2}*T_{2,1} */

	i__1 = *k - l;
	ctrmm_("Left", "Lower", "No transpose", "Non-unit", &l, &i__1, &c_b3, 
		&t[*k - l + 1 + (*k - l + 1) * t_dim1], ldt, &t[*k - l + 1 + 
		t_dim1], ldt);

/*        T_{2,1} = T_{2,1}*T_{1,1} */

	i__1 = *k - l;
	ctrmm_("Right", "Lower", "No transpose", "Non-unit", &l, &i__1, &c_b1,
		 &t[t_offset], ldt, &t[*k - l + 1 + t_dim1], ldt);
    } else {

/*        Else means RQ case */

/*        Break V apart into 6 components */

/*        V = |-----------------------| */
/*            |V_{1,1} V_{1,2} 0      | */
/*            |V_{2,1} V_{2,2} V_{2,3}| */
/*            |-----------------------| */

/*        V_{1,1}\in\C^{k-l,n-k}  rectangular */
/*        V_{1,2}\in\C^{k-l,k-l}  unit lower triangular */

/*        V_{2,1}\in\C^{l,n-k}    rectangular */
/*        V_{2,2}\in\C^{l,k-l}    rectangular */
/*        V_{2,3}\in\C^{l,l}      unit lower triangular */

/*        We will construct the T matrix */
/*        T = |---------------| */
/*            |T_{1,1} 0      | */
/*            |T_{2,1} T_{2,2}| */
/*            |---------------| */

/*        T is the triangular factor obtained from block reflectors. */
/*        To motivate the structure, assume we have already computed T_{1,1} */
/*        and T_{2,2}. Then collect the associated reflectors in V_1 and V_2 */

/*        T_{1,1}\in\C^{k-l, k-l} non-unit lower triangular */
/*        T_{2,2}\in\C^{l, l}     non-unit lower triangular */
/*        T_{2,1}\in\C^{k-l, l}   rectangular */

/*        Where l = floor(k/2) */

/*        Then, consider the product: */

/*        (I - V_2'*T_{2,2}*V_2)*(I - V_1'*T_{1,1}*V_1) */
/*        = I - V_2'*T_{2,2}*V_2 - V_1'*T_{1,1}*V_1 + V_2'*T_{2,2}*V_2*V_1'*T_{1,1}*V_1 */

/*        Define T_{2,1} = -T_{2,2}*V_2*V_1'*T_{1,1} */

/*        Then, we can define the matrix V as */
/*        V = |---| */
/*            |V_1| */
/*            |V_2| */
/*            |---| */

/*        So, our product is equivalent to the matrix product */
/*        I - V'*T*V */
/*        This means, we can compute T_{1,1} and T_{2,2}, then use this information */
/*        to compute T_{2,1} */

/*        Compute T_{1,1} recursively */

	i__1 = *n - l;
	i__3 = *k - l;
	clarft_(direct, storev, &i__1, &i__3, &v[v_offset], ldv, &tau[1], &t[
		t_offset], ldt);

/*        Compute T_{2,2} recursively */

	clarft_(direct, storev, n, &l, &v[*k - l + 1 + v_dim1], ldv, &tau[*k 
		- l + 1], &t[*k - l + 1 + (*k - l + 1) * t_dim1], ldt);

/*        Compute T_{2,1} */
/*        T_{2,1} = V_{2,2} */

	i__1 = *k - l;
	clacpy_("All", &l, &i__1, &v[*k - l + 1 + (*n - *k + 1) * v_dim1], 
		ldv, &t[*k - l + 1 + t_dim1], ldt);

/*        T_{2,1} = T_{2,1}*V_{1,2}' */

	i__1 = *k - l;
	ctrmm_("Right", "Lower", "Conjugate", "Unit", &l, &i__1, &c_b1, &v[(*
		n - *k + 1) * v_dim1 + 1], ldv, &t[*k - l + 1 + t_dim1], ldt);

/*        T_{2,1} = V_{2,1}*V_{1,1}' + T_{2,1} */
/*        Note: We assume K <= N, and GEMM will do nothing if N=K */

	i__1 = *k - l;
	i__3 = *n - *k;
	cgemm_("No transpose", "Conjugate", &l, &i__1, &i__3, &c_b1, &v[*k - 
		l + 1 + v_dim1], ldv, &v[v_offset], ldv, &c_b1, &t[*k - l + 1 
		+ t_dim1], ldt);

/*        At this point, we have that T_{2,1} = V_2*V_1' */
/*        All that is left is to pre and post multiply by -T_{2,2} and T_{1,1} */
/*        respectively. */

/*        T_{2,1} = -T_{2,2}*T_{2,1} */

	i__1 = *k - l;
	ctrmm_("Left", "Lower", "No tranpose", "Non-unit", &l, &i__1, &c_b3, &
		t[*k - l + 1 + (*k - l + 1) * t_dim1], ldt, &t[*k - l + 1 + 
		t_dim1], ldt);

/*        T_{2,1} = T_{2,1}*T_{1,1} */

	i__1 = *k - l;
	ctrmm_("Right", "Lower", "No tranpose", "Non-unit", &l, &i__1, &c_b1, 
		&t[t_offset], ldt, &t[*k - l + 1 + t_dim1], ldt);
    }
    return;
} /* clarft_ */

