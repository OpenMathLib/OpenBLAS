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
#define s_cat(lpp, rpp, rnp, np, llp) { 	ftnlen i, nc, ll; char *f__rp, *lp; 	ll = (llp); lp = (lpp); 	for(i=0; i < (int)*(np); ++i) {         	nc = ll; 	        if((rnp)[i] < nc) nc = (rnp)[i]; 	        ll -= nc;         	f__rp = (rpp)[i]; 	        while(--nc >= 0) *lp++ = *(f__rp)++;         } 	while(--ll >= 0) *lp++ = ' '; }
#define s_cmp(a,b,c,d) ((integer)strncmp((a),(b),f2cmin((c),(d))))
#define s_copy(A,B,C,D) { int __i,__m; for (__i=0, __m=f2cmin((C),(D)); __i<__m && (B)[__i] != 0; ++__i) (A)[__i] = (B)[__i]; }
#define sig_die(s, kill) { exit(1); }
#define s_stop(s, n) {exit(0);}
#define myexit_() break;
#define mycycle_() continue;
#define myceiling_(w) {ceil(w)}
#define myhuge_(w) {HUGE_VAL}
//#define mymaxloc_(w,s,e,n) {if (sizeof(*(w)) == sizeof(double)) dmaxloc_((w),*(s),*(e),n); else dmaxloc_((w),*(s),*(e),n);}
#define mymaxloc_(w,s,e,n) dmaxloc_(w,*(s),*(e),n)

/* procedure parameter types for -A and -C++ */

#define F2C_proc_par_types 1
#ifdef __cplusplus
typedef logical (*L_fp)(...);
#else
typedef logical (*L_fp)();
#endif


/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;

/* Subroutine */ void slaqz2_(logical *ilq, logical *ilz, integer *k, integer *
	istartm, integer *istopm, integer *ihi, real *a, integer *lda, real *
	b, integer *ldb, integer *nq, integer *qstart, real *q, integer *ldq, 
	integer *nz, integer *zstart, real *z__, integer *ldz)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1;

    /* Local variables */
    static real h__[6]	/* was [2][3] */, c1, c2, s1, s2, temp;
    extern /* Subroutine */ void srot_(integer *, real *, integer *, real *, 
	    integer *, real *, real *), slartg_(real *, real *, real *, real *
	    , real *);


/*     Arguments */

/*     Parameters */

/*     Local variables */

/*     External functions */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;

    /* Function Body */
    if (*k + 2 == *ihi) {
/*        Shift is located on the edge of the matrix, remove it */
/*         H = B( IHI-1:IHI, IHI-2:IHI ) */
	h__[0] = b[*ihi - 1 + (*ihi - 2) * b_dim1];
	h__[1] = b[*ihi + (*ihi - 2) * b_dim1];
	h__[2] = b[*ihi - 1 + (*ihi - 1) * b_dim1];
	h__[4] = b[*ihi - 1 + *ihi * b_dim1];
	h__[3] = b[*ihi + (*ihi - 1) * b_dim1];
	h__[5] = b[*ihi + *ihi * b_dim1];
/*        Make H upper triangular */
	slartg_(h__, &h__[1], &c1, &s1, &temp);
	h__[1] = 0.f;
	h__[0] = temp;
	srot_(&c__2, &h__[2], &c__2, &h__[3], &c__2, &c1, &s1);

	slartg_(&h__[5], &h__[3], &c1, &s1, &temp);
	srot_(&c__1, &h__[4], &c__1, &h__[2], &c__1, &c1, &s1);
	slartg_(&h__[2], h__, &c2, &s2, &temp);

	i__1 = *ihi - *istartm + 1;
	srot_(&i__1, &b[*istartm + *ihi * b_dim1], &c__1, &b[*istartm + (*ihi 
		- 1) * b_dim1], &c__1, &c1, &s1);
	i__1 = *ihi - *istartm + 1;
	srot_(&i__1, &b[*istartm + (*ihi - 1) * b_dim1], &c__1, &b[*istartm + 
		(*ihi - 2) * b_dim1], &c__1, &c2, &s2);
	b[*ihi - 1 + (*ihi - 2) * b_dim1] = 0.f;
	b[*ihi + (*ihi - 2) * b_dim1] = 0.f;
	i__1 = *ihi - *istartm + 1;
	srot_(&i__1, &a[*istartm + *ihi * a_dim1], &c__1, &a[*istartm + (*ihi 
		- 1) * a_dim1], &c__1, &c1, &s1);
	i__1 = *ihi - *istartm + 1;
	srot_(&i__1, &a[*istartm + (*ihi - 1) * a_dim1], &c__1, &a[*istartm + 
		(*ihi - 2) * a_dim1], &c__1, &c2, &s2);
	if (*ilz) {
	    srot_(nz, &z__[(*ihi - *zstart + 1) * z_dim1 + 1], &c__1, &z__[(*
		    ihi - 1 - *zstart + 1) * z_dim1 + 1], &c__1, &c1, &s1);
	    srot_(nz, &z__[(*ihi - 1 - *zstart + 1) * z_dim1 + 1], &c__1, &
		    z__[(*ihi - 2 - *zstart + 1) * z_dim1 + 1], &c__1, &c2, &
		    s2);
	}

	slartg_(&a[*ihi - 1 + (*ihi - 2) * a_dim1], &a[*ihi + (*ihi - 2) * 
		a_dim1], &c1, &s1, &temp);
	a[*ihi - 1 + (*ihi - 2) * a_dim1] = temp;
	a[*ihi + (*ihi - 2) * a_dim1] = 0.f;
	i__1 = *istopm - *ihi + 2;
	srot_(&i__1, &a[*ihi - 1 + (*ihi - 1) * a_dim1], lda, &a[*ihi + (*ihi 
		- 1) * a_dim1], lda, &c1, &s1);
	i__1 = *istopm - *ihi + 2;
	srot_(&i__1, &b[*ihi - 1 + (*ihi - 1) * b_dim1], ldb, &b[*ihi + (*ihi 
		- 1) * b_dim1], ldb, &c1, &s1);
	if (*ilq) {
	    srot_(nq, &q[(*ihi - 1 - *qstart + 1) * q_dim1 + 1], &c__1, &q[(*
		    ihi - *qstart + 1) * q_dim1 + 1], &c__1, &c1, &s1);
	}

	slartg_(&b[*ihi + *ihi * b_dim1], &b[*ihi + (*ihi - 1) * b_dim1], &c1,
		 &s1, &temp);
	b[*ihi + *ihi * b_dim1] = temp;
	b[*ihi + (*ihi - 1) * b_dim1] = 0.f;
	i__1 = *ihi - *istartm;
	srot_(&i__1, &b[*istartm + *ihi * b_dim1], &c__1, &b[*istartm + (*ihi 
		- 1) * b_dim1], &c__1, &c1, &s1);
	i__1 = *ihi - *istartm + 1;
	srot_(&i__1, &a[*istartm + *ihi * a_dim1], &c__1, &a[*istartm + (*ihi 
		- 1) * a_dim1], &c__1, &c1, &s1);
	if (*ilz) {
	    srot_(nz, &z__[(*ihi - *zstart + 1) * z_dim1 + 1], &c__1, &z__[(*
		    ihi - 1 - *zstart + 1) * z_dim1 + 1], &c__1, &c1, &s1);
	}

    } else {

/*        Normal operation, move bulge down */

/*         H = B( K+1:K+2, K:K+2 ) */
	h__[0] = b[*k + 1 + *k * b_dim1];
	h__[1] = b[*k + 2 + *k * b_dim1];
	h__[2] = b[*k + 1 + (*k + 1) * b_dim1];
	h__[4] = b[*k + 1 + (*k + 2) * b_dim1];
	h__[3] = b[*k + 2 + (*k + 1) * b_dim1];
	h__[5] = b[*k + 2 + (*k + 2) * b_dim1];

/*        Make H upper triangular */

	slartg_(h__, &h__[1], &c1, &s1, &temp);
	h__[1] = 0.f;
	h__[0] = temp;
	srot_(&c__2, &h__[2], &c__2, &h__[3], &c__2, &c1, &s1);

/*        Calculate Z1 and Z2 */

	slartg_(&h__[5], &h__[3], &c1, &s1, &temp);
	srot_(&c__1, &h__[4], &c__1, &h__[2], &c__1, &c1, &s1);
	slartg_(&h__[2], h__, &c2, &s2, &temp);

/*        Apply transformations from the right */

	i__1 = *k + 3 - *istartm + 1;
	srot_(&i__1, &a[*istartm + (*k + 2) * a_dim1], &c__1, &a[*istartm + (*
		k + 1) * a_dim1], &c__1, &c1, &s1);
	i__1 = *k + 3 - *istartm + 1;
	srot_(&i__1, &a[*istartm + (*k + 1) * a_dim1], &c__1, &a[*istartm + *
		k * a_dim1], &c__1, &c2, &s2);
	i__1 = *k + 2 - *istartm + 1;
	srot_(&i__1, &b[*istartm + (*k + 2) * b_dim1], &c__1, &b[*istartm + (*
		k + 1) * b_dim1], &c__1, &c1, &s1);
	i__1 = *k + 2 - *istartm + 1;
	srot_(&i__1, &b[*istartm + (*k + 1) * b_dim1], &c__1, &b[*istartm + *
		k * b_dim1], &c__1, &c2, &s2);
	if (*ilz) {
	    srot_(nz, &z__[(*k + 2 - *zstart + 1) * z_dim1 + 1], &c__1, &z__[(
		    *k + 1 - *zstart + 1) * z_dim1 + 1], &c__1, &c1, &s1);
	    srot_(nz, &z__[(*k + 1 - *zstart + 1) * z_dim1 + 1], &c__1, &z__[(
		    *k - *zstart + 1) * z_dim1 + 1], &c__1, &c2, &s2);
	}
	b[*k + 1 + *k * b_dim1] = 0.f;
	b[*k + 2 + *k * b_dim1] = 0.f;

/*        Calculate Q1 and Q2 */

	slartg_(&a[*k + 2 + *k * a_dim1], &a[*k + 3 + *k * a_dim1], &c1, &s1, 
		&temp);
	a[*k + 2 + *k * a_dim1] = temp;
	a[*k + 3 + *k * a_dim1] = 0.f;
	slartg_(&a[*k + 1 + *k * a_dim1], &a[*k + 2 + *k * a_dim1], &c2, &s2, 
		&temp);
	a[*k + 1 + *k * a_dim1] = temp;
	a[*k + 2 + *k * a_dim1] = 0.f;

/*     Apply transformations from the left */

	i__1 = *istopm - *k;
	srot_(&i__1, &a[*k + 2 + (*k + 1) * a_dim1], lda, &a[*k + 3 + (*k + 1)
		 * a_dim1], lda, &c1, &s1);
	i__1 = *istopm - *k;
	srot_(&i__1, &a[*k + 1 + (*k + 1) * a_dim1], lda, &a[*k + 2 + (*k + 1)
		 * a_dim1], lda, &c2, &s2);

	i__1 = *istopm - *k;
	srot_(&i__1, &b[*k + 2 + (*k + 1) * b_dim1], ldb, &b[*k + 3 + (*k + 1)
		 * b_dim1], ldb, &c1, &s1);
	i__1 = *istopm - *k;
	srot_(&i__1, &b[*k + 1 + (*k + 1) * b_dim1], ldb, &b[*k + 2 + (*k + 1)
		 * b_dim1], ldb, &c2, &s2);
	if (*ilq) {
	    srot_(nq, &q[(*k + 2 - *qstart + 1) * q_dim1 + 1], &c__1, &q[(*k 
		    + 3 - *qstart + 1) * q_dim1 + 1], &c__1, &c1, &s1);
	    srot_(nq, &q[(*k + 1 - *qstart + 1) * q_dim1 + 1], &c__1, &q[(*k 
		    + 2 - *qstart + 1) * q_dim1 + 1], &c__1, &c2, &s2);
	}

    }

/*     End of SLAQZ2 */

    return;
} /* slaqz2_ */

