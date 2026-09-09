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

/* procedure parameter types for -A and -C++ */


#ifdef __cplusplus
typedef logical (*L_fp)(...);
#else
typedef logical (*L_fp)();
#endif
#if 0
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
#endif
/*  -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


/* Table of constant values */

static integer c__3 = 3;
static integer c_n1 = -1;
static integer c__2 = 2;
static real c_b13 = 1.f;
static real c_b22 = -1.f;

/* > \brief \b SLARFT forms the triangular factor T of a block reflector H = I - vtvH */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > Download SLARFT + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarft.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarft.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarft.
f"> */
/* > [TXT]</a> */

/*  Definition: */
/*  =========== */

/*        SUBROUTINE SLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT ) */

/*       CHARACTER          DIRECT, STOREV */
/*       INTEGER            K, LDT, LDV, N */
/*       REAL               T( LDT, * ), TAU( * ), V( LDV, * ) */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLARFT forms the triangular factor T of a real block reflector H */
/* > of order n, which is defined as a product of k elementary reflectors. */
/* > */
/* > If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular; */
/* > */
/* > If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular. */
/* > */
/* > If STOREV = 'C', the vector which defines the elementary reflector */
/* > H(i) is stored in the i-th column of the array V, and */
/* > */
/* >    H  =  I - V * T * V**T */
/* > */
/* > If STOREV = 'R', the vector which defines the elementary reflector */
/* > H(i) is stored in the i-th row of the array V, and */
/* > */
/* >    H  =  I - V**T * T * V */
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
/* >          V is REAL array, dimension */
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
/* >          TAU is REAL array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i). */
/* > \endverbatim */
/* > */
/* > \param[out] T */
/* > \verbatim */
/* >          T is REAL array, dimension (LDT,K) */
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
/* > \author Johnathan Rhyne, Univ. of Colorado Denver (original author, 2024) */
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
/* Subroutine */ void slarft_(char *direct, char *storev, integer *n, integer *
	k, real *v, integer *ldv, real *tau, real *t, integer *ldt)
{
    /* System generated locals */
    address a__1[2];
    integer t_dim1, t_offset, v_dim1, v_offset, i__1[2], i__2, i__3;
    char ch__1[2];

    /* Local variables */
    integer i__, j, l;
    logical lq, ql, qr;
    integer nx;
    extern /* Subroutine */ void slarft_lvl2__(char *, char *, integer *, 
	    integer *, real *, integer *, real *, real *, integer *);
    logical dirf, colv;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ void sgemm_(char *, char *, integer *, integer *, 
	    integer *, real *, real *, integer *, real *, integer *, real *, 
	    real *, integer *), strmm_(char *, char *, char *,
	     char *, integer *, integer *, real *, real *, integer *, real *, 
	    integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ void slacpy_(char *, integer *, integer *, real *, 
	    integer *, real *, integer *);


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
	t[t_dim1 + 1] = tau[1];
	return;
    }

/*     Determine when to cross over into the level 2 based implementation */

/* Writing concatenation */
    i__1[0] = 1, a__1[0] = direct;
    i__1[1] = 1, a__1[1] = storev;
    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
    nx = ilaenv_(&c__3, "SLARFT", ch__1, n, k, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)2);
    if (*k < nx) {
	slarft_lvl2__(direct, storev, n, k, &v[v_offset], ldv, &tau[1], &t[
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

/*        V_{1,1}\in\R^{l,l}      unit lower triangular */
/*        V_{2,1}\in\R^{k-l,l}    rectangular */
/*        V_{3,1}\in\R^{n-k,l}    rectangular */

/*        V_{2,2}\in\R^{k-l,k-l}  unit lower triangular */
/*        V_{3,2}\in\R^{n-k,k-l}  rectangular */

/*        We will construct the T matrix */
/*        T = |---------------| */
/*            |T_{1,1} T_{1,2}| */
/*            |0       T_{2,2}| */
/*            |---------------| */

/*        T is the triangular factor obtained from block reflectors. */
/*        To motivate the structure, assume we have already computed T_{1,1} */
/*        and T_{2,2}. Then collect the associated reflectors in V_1 and V_2 */

/*        T_{1,1}\in\R^{l, l}     upper triangular */
/*        T_{2,2}\in\R^{k-l, k-l} upper triangular */
/*        T_{1,2}\in\R^{l, k-l}   rectangular */

/*        Where l = floor(k/2) */

/*        Then, consider the product: */

/*        (I - V_1*T_{1,1}*V_1')*(I - V_2*T_{2,2}*V_2') */
/*        = I - V_1*T_{1,1}*V_1' - V_2*T_{2,2}*V_2' + V_1*T_{1,1}*V_1'*V_2*T_{2,2}*V_2' */

/*        Define T_{1,2} = -T_{1,1}*V_1'*V_2*T_{2,2} */

/*        Then, we can define the matrix V as */
/*        V = |-------| */
/*            |V_1 V_2| */
/*            |-------| */

/*        So, our product is equivalent to the matrix product */
/*        I - V*T*V' */
/*        This means, we can compute T_{1,1} and T_{2,2}, then use this information */
/*        to compute T_{1,2} */

/*        Compute T_{1,1} recursively */

	slarft_(direct, storev, n, &l, &v[v_offset], ldv, &tau[1], &t[
		t_offset], ldt);

/*        Compute T_{2,2} recursively */

	i__2 = *n - l;
	i__3 = *k - l;
	slarft_(direct, storev, &i__2, &i__3, &v[l + 1 + (l + 1) * v_dim1], 
		ldv, &tau[l + 1], &t[l + 1 + (l + 1) * t_dim1], ldt);

/*        Compute T_{1,2} */
/*        T_{1,2} = V_{2,1}' */

	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = *k - l;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		t[j + (l + i__) * t_dim1] = v[l + i__ + j * v_dim1];
	    }
	}

/*        T_{1,2} = T_{1,2}*V_{2,2} */

	i__2 = *k - l;
	strmm_("Right", "Lower", "No transpose", "Unit", &l, &i__2, &c_b13, &
		v[l + 1 + (l + 1) * v_dim1], ldv, &t[(l + 1) * t_dim1 + 1], 
		ldt);

/*        T_{1,2} = V_{3,1}'*V_{3,2} + T_{1,2} */
/*        Note: We assume K <= N, and GEMM will do nothing if N=K */

	i__2 = *k - l;
	i__3 = *n - *k;
	sgemm_("Transpose", "No transpose", &l, &i__2, &i__3, &c_b13, &v[*k + 
		1 + v_dim1], ldv, &v[*k + 1 + (l + 1) * v_dim1], ldv, &c_b13, 
		&t[(l + 1) * t_dim1 + 1], ldt);

/*        At this point, we have that T_{1,2} = V_1'*V_2 */
/*        All that is left is to pre and post multiply by -T_{1,1} and T_{2,2} */
/*        respectively. */

/*        T_{1,2} = -T_{1,1}*T_{1,2} */

	i__2 = *k - l;
	strmm_("Left", "Upper", "No transpose", "Non-unit", &l, &i__2, &c_b22,
		 &t[t_offset], ldt, &t[(l + 1) * t_dim1 + 1], ldt);

/*        T_{1,2} = T_{1,2}*T_{2,2} */

	i__2 = *k - l;
	strmm_("Right", "Upper", "No transpose", "Non-unit", &l, &i__2, &
		c_b13, &t[l + 1 + (l + 1) * t_dim1], ldt, &t[(l + 1) * t_dim1 
		+ 1], ldt);
    } else if (lq) {

/*        Break V apart into 6 components */

/*        V = |----------------------| */
/*            |V_{1,1} V_{1,2} V{1,3}| */
/*            |0       V_{2,2} V{2,3}| */
/*            |----------------------| */

/*        V_{1,1}\in\R^{l,l}      unit upper triangular */
/*        V_{1,2}\in\R^{l,k-l}    rectangular */
/*        V_{1,3}\in\R^{l,n-k}    rectangular */

/*        V_{2,2}\in\R^{k-l,k-l}  unit upper triangular */
/*        V_{2,3}\in\R^{k-l,n-k}  rectangular */

/*        Where l = floor(k/2) */

/*        We will construct the T matrix */
/*        T = |---------------| */
/*            |T_{1,1} T_{1,2}| */
/*            |0       T_{2,2}| */
/*            |---------------| */

/*        T is the triangular factor obtained from block reflectors. */
/*        To motivate the structure, assume we have already computed T_{1,1} */
/*        and T_{2,2}. Then collect the associated reflectors in V_1 and V_2 */

/*        T_{1,1}\in\R^{l, l}     upper triangular */
/*        T_{2,2}\in\R^{k-l, k-l} upper triangular */
/*        T_{1,2}\in\R^{l, k-l}   rectangular */

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

	slarft_(direct, storev, n, &l, &v[v_offset], ldv, &tau[1], &t[
		t_offset], ldt);

/*        Compute T_{2,2} recursively */

	i__2 = *n - l;
	i__3 = *k - l;
	slarft_(direct, storev, &i__2, &i__3, &v[l + 1 + (l + 1) * v_dim1], 
		ldv, &tau[l + 1], &t[l + 1 + (l + 1) * t_dim1], ldt);

/*        Compute T_{1,2} */
/*        T_{1,2} = V_{1,2} */

	i__2 = *k - l;
	slacpy_("All", &l, &i__2, &v[(l + 1) * v_dim1 + 1], ldv, &t[(l + 1) * 
		t_dim1 + 1], ldt);

/*        T_{1,2} = T_{1,2}*V_{2,2}' */

	i__2 = *k - l;
	strmm_("Right", "Upper", "Transpose", "Unit", &l, &i__2, &c_b13, &v[l 
		+ 1 + (l + 1) * v_dim1], ldv, &t[(l + 1) * t_dim1 + 1], ldt);

/*        T_{1,2} = V_{1,3}*V_{2,3}' + T_{1,2} */
/*        Note: We assume K <= N, and GEMM will do nothing if N=K */

	i__2 = *k - l;
	i__3 = *n - *k;
	sgemm_("No transpose", "Transpose", &l, &i__2, &i__3, &c_b13, &v[(*k 
		+ 1) * v_dim1 + 1], ldv, &v[l + 1 + (*k + 1) * v_dim1], ldv, &
		c_b13, &t[(l + 1) * t_dim1 + 1], ldt);

/*        At this point, we have that T_{1,2} = V_1*V_2' */
/*        All that is left is to pre and post multiply by -T_{1,1} and T_{2,2} */
/*        respectively. */

/*        T_{1,2} = -T_{1,1}*T_{1,2} */

	i__2 = *k - l;
	strmm_("Left", "Upper", "No transpose", "Non-unit", &l, &i__2, &c_b22,
		 &t[t_offset], ldt, &t[(l + 1) * t_dim1 + 1], ldt);

/*        T_{1,2} = T_{1,2}*T_{2,2} */

	i__2 = *k - l;
	strmm_("Right", "Upper", "No transpose", "Non-unit", &l, &i__2, &
		c_b13, &t[l + 1 + (l + 1) * t_dim1], ldt, &t[(l + 1) * t_dim1 
		+ 1], ldt);
    } else if (ql) {

/*        Break V apart into 6 components */

/*        V = |---------------| */
/*            |V_{1,1} V_{1,2}| */
/*            |V_{2,1} V_{2,2}| */
/*            |0       V_{3,2}| */
/*            |---------------| */

/*        V_{1,1}\in\R^{n-k,k-l}  rectangular */
/*        V_{2,1}\in\R^{k-l,k-l}  unit upper triangular */

/*        V_{1,2}\in\R^{n-k,l}    rectangular */
/*        V_{2,2}\in\R^{k-l,l}    rectangular */
/*        V_{3,2}\in\R^{l,l}      unit upper triangular */

/*        We will construct the T matrix */
/*        T = |---------------| */
/*            |T_{1,1} 0      | */
/*            |T_{2,1} T_{2,2}| */
/*            |---------------| */

/*        T is the triangular factor obtained from block reflectors. */
/*        To motivate the structure, assume we have already computed T_{1,1} */
/*        and T_{2,2}. Then collect the associated reflectors in V_1 and V_2 */

/*        T_{1,1}\in\R^{k-l, k-l} non-unit lower triangular */
/*        T_{2,2}\in\R^{l, l}     non-unit lower triangular */
/*        T_{2,1}\in\R^{k-l, l}   rectangular */

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

	i__2 = *n - l;
	i__3 = *k - l;
	slarft_(direct, storev, &i__2, &i__3, &v[v_offset], ldv, &tau[1], &t[
		t_offset], ldt);

/*        Compute T_{2,2} recursively */

	slarft_(direct, storev, n, &l, &v[(*k - l + 1) * v_dim1 + 1], ldv, &
		tau[*k - l + 1], &t[*k - l + 1 + (*k - l + 1) * t_dim1], ldt);

/*        Compute T_{2,1} */
/*        T_{2,1} = V_{2,2}' */

	i__2 = *k - l;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = l;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		t[*k - l + i__ + j * t_dim1] = v[*n - *k + j + (*k - l + i__) 
			* v_dim1];
	    }
	}

/*        T_{2,1} = T_{2,1}*V_{2,1} */

	i__2 = *k - l;
	strmm_("Right", "Upper", "No transpose", "Unit", &l, &i__2, &c_b13, &
		v[*n - *k + 1 + v_dim1], ldv, &t[*k - l + 1 + t_dim1], ldt);

/*        T_{2,1} = V_{2,2}'*V_{2,1} + T_{2,1} */
/*        Note: We assume K <= N, and GEMM will do nothing if N=K */

	i__2 = *k - l;
	i__3 = *n - *k;
	sgemm_("Transpose", "No transpose", &l, &i__2, &i__3, &c_b13, &v[(*k 
		- l + 1) * v_dim1 + 1], ldv, &v[v_offset], ldv, &c_b13, &t[*k 
		- l + 1 + t_dim1], ldt);

/*        At this point, we have that T_{2,1} = V_2'*V_1 */
/*        All that is left is to pre and post multiply by -T_{2,2} and T_{1,1} */
/*        respectively. */

/*        T_{2,1} = -T_{2,2}*T_{2,1} */

	i__2 = *k - l;
	strmm_("Left", "Lower", "No transpose", "Non-unit", &l, &i__2, &c_b22,
		 &t[*k - l + 1 + (*k - l + 1) * t_dim1], ldt, &t[*k - l + 1 + 
		t_dim1], ldt);

/*        T_{2,1} = T_{2,1}*T_{1,1} */

	i__2 = *k - l;
	strmm_("Right", "Lower", "No transpose", "Non-unit", &l, &i__2, &
		c_b13, &t[t_offset], ldt, &t[*k - l + 1 + t_dim1], ldt);
    } else {

/*        Else means RQ case */

/*        Break V apart into 6 components */

/*        V = |-----------------------| */
/*            |V_{1,1} V_{1,2} 0      | */
/*            |V_{2,1} V_{2,2} V_{2,3}| */
/*            |-----------------------| */

/*        V_{1,1}\in\R^{k-l,n-k}  rectangular */
/*        V_{1,2}\in\R^{k-l,k-l}  unit lower triangular */

/*        V_{2,1}\in\R^{l,n-k}    rectangular */
/*        V_{2,2}\in\R^{l,k-l}    rectangular */
/*        V_{2,3}\in\R^{l,l}      unit lower triangular */

/*        We will construct the T matrix */
/*        T = |---------------| */
/*            |T_{1,1} 0      | */
/*            |T_{2,1} T_{2,2}| */
/*            |---------------| */

/*        T is the triangular factor obtained from block reflectors. */
/*        To motivate the structure, assume we have already computed T_{1,1} */
/*        and T_{2,2}. Then collect the associated reflectors in V_1 and V_2 */

/*        T_{1,1}\in\R^{k-l, k-l} non-unit lower triangular */
/*        T_{2,2}\in\R^{l, l}     non-unit lower triangular */
/*        T_{2,1}\in\R^{k-l, l}   rectangular */

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
/*        I - V'TV */
/*        This means, we can compute T_{1,1} and T_{2,2}, then use this information */
/*        to compute T_{2,1} */

/*        Compute T_{1,1} recursively */

	i__2 = *n - l;
	i__3 = *k - l;
	slarft_(direct, storev, &i__2, &i__3, &v[v_offset], ldv, &tau[1], &t[
		t_offset], ldt);

/*        Compute T_{2,2} recursively */

	slarft_(direct, storev, n, &l, &v[*k - l + 1 + v_dim1], ldv, &tau[*k 
		- l + 1], &t[*k - l + 1 + (*k - l + 1) * t_dim1], ldt);

/*        Compute T_{2,1} */
/*        T_{2,1} = V_{2,2} */

	i__2 = *k - l;
	slacpy_("All", &l, &i__2, &v[*k - l + 1 + (*n - *k + 1) * v_dim1], 
		ldv, &t[*k - l + 1 + t_dim1], ldt);

/*        T_{2,1} = T_{2,1}*V_{1,2}' */

	i__2 = *k - l;
	strmm_("Right", "Lower", "Transpose", "Unit", &l, &i__2, &c_b13, &v[(*
		n - *k + 1) * v_dim1 + 1], ldv, &t[*k - l + 1 + t_dim1], ldt);

/*        T_{2,1} = V_{2,1}*V_{1,1}' + T_{2,1} */
/*        Note: We assume K <= N, and GEMM will do nothing if N=K */

	i__2 = *k - l;
	i__3 = *n - *k;
	sgemm_("No transpose", "Transpose", &l, &i__2, &i__3, &c_b13, &v[*k - 
		l + 1 + v_dim1], ldv, &v[v_offset], ldv, &c_b13, &t[*k - l + 
		1 + t_dim1], ldt);

/*        At this point, we have that T_{2,1} = V_2*V_1' */
/*        All that is left is to pre and post multiply by -T_{2,2} and T_{1,1} */
/*        respectively. */

/*        T_{2,1} = -T_{2,2}*T_{2,1} */

	i__2 = *k - l;
	strmm_("Left", "Lower", "No tranpose", "Non-unit", &l, &i__2, &c_b22, 
		&t[*k - l + 1 + (*k - l + 1) * t_dim1], ldt, &t[*k - l + 1 + 
		t_dim1], ldt);

/*        T_{2,1} = T_{2,1}*T_{1,1} */

	i__2 = *k - l;
	strmm_("Right", "Lower", "No tranpose", "Non-unit", &l, &i__2, &c_b13,
		 &t[t_offset], ldt, &t[*k - l + 1 + t_dim1], ldt);
    }
    return;
} /* slarft_ */

