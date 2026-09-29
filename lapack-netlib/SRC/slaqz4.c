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

static real c_b4 = 0.f;
static real c_b5 = 1.f;
static integer c__1 = 1;
static logical c_true = TRUE_;

/* Subroutine */ void slaqz4_(logical *ilschur, logical *ilq, logical *ilz, 
	integer *n, integer *ilo, integer *ihi, integer *nshifts, integer *
	nblock_desired__, real *sr, real *si, real *ss, real *a, integer *lda,
	 real *b, integer *ldb, real *q, integer *ldq, real *z__, integer *
	ldz, real *qc, integer *ldqc, real *zc, integer *ldzc, real *work, 
	integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, qc_dim1, qc_offset, zc_dim1, zc_offset, i__1, i__2, 
	    i__3, i__4, i__5;

    /* Local variables */
    static integer i__, j, k;
    static real v[3], c1, c2, s1, s2;
    static integer np, ns;
    static real temp, swap;
    static integer npos;
    extern /* Subroutine */ void srot_(integer *, real *, integer *, real *, 
	    integer *, real *, real *), sgemm_(char *, char *, integer *, 
	    integer *, integer *, real *, real *, integer *, real *, integer *
	    , real *, real *, integer * /*, ftnlen, ftnlen*/), slaqz1_(real *, 
	    integer *, real *, integer *, real *, real *, real *, real *, 
	    real *, real *), slaqz2_(logical *, logical *, integer *, integer 
	    *, integer *, integer *, real *, integer *, real *, integer *, 
	    integer *, integer *, real *, integer *, integer *, integer *, 
	    real *, integer *);
    static integer nblock;
    extern /* Subroutine */ void xerbla_(char *, integer *, ftnlen);
    static integer ishift;
    extern /* Subroutine */ void slaset_(char *, integer *, integer *, real *, 
	    real *, real *, integer * /*, ftnlen*/), slartg_(real *, real *, real *
	    , real *, real *), slacpy_(char *, integer *, integer *, real *, 
	    integer *, real *, integer * /*, ftnlen*/);
    static integer istopb, swidth, istopm, sheight;
    extern real sroundup_lwork__(integer *);
    static integer istartb, istartm;

/*     Function arguments */
/*     Parameters */
/*     Local scalars */

/*     External functions */
    /* Parameter adjustments */
    --sr;
    --si;
    --ss;
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
    qc_dim1 = *ldqc;
    qc_offset = 1 + qc_dim1;
    qc -= qc_offset;
    zc_dim1 = *ldzc;
    zc_offset = 1 + zc_dim1;
    zc -= zc_offset;
    --work;

    /* Function Body */
    *info = 0;
    if (*nblock_desired__ < *nshifts + 1) {
	*info = -8;
    }
    if (*lwork == -1) {
/*        workspace query, quick return */
	i__1 = *n * *nblock_desired__;
	work[1] = sroundup_lwork__(&i__1);
	return;
    } else if (*lwork < *n * *nblock_desired__) {
	*info = -25;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SLAQZ4", &i__1 , (ftnlen)6);
	return;
    }
/*     Executable statements */
    if (*nshifts < 2) {
	return;
    }
    if (*ilo >= *ihi) {
	return;
    }
    if (*ilschur) {
	istartm = 1;
	istopm = *n;
    } else {
	istartm = *ilo;
	istopm = *ihi;
    }
/*     Shuffle shifts into pairs of real shifts and pairs */
/*     of complex conjugate shifts assuming complex */
/*     conjugate shifts are already adjacent to one */
/*     another */
    i__1 = *nshifts - 2;
    for (i__ = 1; i__ <= i__1; i__ += 2) {
	if (si[i__] != -si[i__ + 1]) {

	    swap = sr[i__];
	    sr[i__] = sr[i__ + 1];
	    sr[i__ + 1] = sr[i__ + 2];
	    sr[i__ + 2] = swap;
	    swap = si[i__];
	    si[i__] = si[i__ + 1];
	    si[i__ + 1] = si[i__ + 2];
	    si[i__ + 2] = swap;
	    swap = ss[i__];
	    ss[i__] = ss[i__ + 1];
	    ss[i__ + 1] = ss[i__ + 2];
	    ss[i__ + 2] = swap;
	}
    }
/*     NSHFTS is supposed to be even, but if it is odd, */
/*     then simply reduce it by one.  The shuffle above */
/*     ensures that the dropped shift is real and that */
/*     the remaining shifts are paired. */
    ns = *nshifts - *nshifts % 2;
/* Computing MAX */
    i__1 = *nblock_desired__ - ns;
    npos = f2cmax(i__1,1);
/*     The following block introduces the shifts and chases */
/*     them down one by one just enough to make space for */
/*     the other shifts. The near-the-diagonal block is */
/*     of size (ns+1) x ns. */
    i__1 = ns + 1;
    i__2 = ns + 1;
    slaset_("FULL", &i__1, &i__2, &c_b4, &c_b5, &qc[qc_offset], ldqc /*, (ftnlen)4*/);
    slaset_("FULL", &ns, &ns, &c_b4, &c_b5, &zc[zc_offset], ldzc /*, (ftnlen)4*/);
    i__1 = ns;
    for (i__ = 1; i__ <= i__1; i__ += 2) {
/*        Introduce the shift */
	slaqz1_(&a[*ilo + *ilo * a_dim1], lda, &b[*ilo + *ilo * b_dim1], ldb, 
		&sr[i__], &sr[i__ + 1], &si[i__], &ss[i__], &ss[i__ + 1], v);
	temp = v[1];
	slartg_(&temp, &v[2], &c1, &s1, &v[1]);
	slartg_(v, &v[1], &c2, &s2, &temp);
	srot_(&ns, &a[*ilo + 1 + *ilo * a_dim1], lda, &a[*ilo + 2 + *ilo * 
		a_dim1], lda, &c1, &s1);
	srot_(&ns, &a[*ilo + *ilo * a_dim1], lda, &a[*ilo + 1 + *ilo * a_dim1]
		, lda, &c2, &s2);
	srot_(&ns, &b[*ilo + 1 + *ilo * b_dim1], ldb, &b[*ilo + 2 + *ilo * 
		b_dim1], ldb, &c1, &s1);
	srot_(&ns, &b[*ilo + *ilo * b_dim1], ldb, &b[*ilo + 1 + *ilo * b_dim1]
		, ldb, &c2, &s2);
	i__2 = ns + 1;
	srot_(&i__2, &qc[(qc_dim1 << 1) + 1], &c__1, &qc[qc_dim1 * 3 + 1], &
		c__1, &c1, &s1);
	i__2 = ns + 1;
	srot_(&i__2, &qc[qc_dim1 + 1], &c__1, &qc[(qc_dim1 << 1) + 1], &c__1, 
		&c2, &s2);
/*        Chase the shift down */
	i__2 = ns - 1 - i__;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = *ihi - *ilo + 1;
	    i__4 = ns + 1;
	    slaqz2_(&c_true, &c_true, &j, &c__1, &ns, &i__3, &a[*ilo + *ilo * 
		    a_dim1], lda, &b[*ilo + *ilo * b_dim1], ldb, &i__4, &c__1,
		     &qc[qc_offset], ldqc, &ns, &c__1, &zc[zc_offset], ldzc);
	}
    }
/*     Update the rest of the pencil */
/*     Update A(ilo:ilo+ns,ilo+ns:istopm) and B(ilo:ilo+ns,ilo+ns:istopm) */
/*     from the left with Qc(1:ns+1,1:ns+1)' */
    sheight = ns + 1;
    swidth = istopm - (*ilo + ns) + 1;
    if (swidth > 0) {
	sgemm_("T", "N", &sheight, &swidth, &sheight, &c_b5, &qc[qc_offset], 
		ldqc, &a[*ilo + (*ilo + ns) * a_dim1], lda, &c_b4, &work[1], &
		sheight/*, (ftnlen)1, (ftnlen)1*/);
	slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[*ilo + (*ilo 
		+ ns) * a_dim1], lda/*, (ftnlen)3*/);
	sgemm_("T", "N", &sheight, &swidth, &sheight, &c_b5, &qc[qc_offset], 
		ldqc, &b[*ilo + (*ilo + ns) * b_dim1], ldb, &c_b4, &work[1], &
		sheight/*, (ftnlen)1, (ftnlen)1*/);
	slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[*ilo + (*ilo 
		+ ns) * b_dim1], ldb/*, (ftnlen)3*/);
    }
    if (*ilq) {
	sgemm_("N", "N", n, &sheight, &sheight, &c_b5, &q[*ilo * q_dim1 + 1], 
		ldq, &qc[qc_offset], ldqc, &c_b4, &work[1], n/*, (ftnlen)1, (
		ftnlen)1*/);
	slacpy_("ALL", n, &sheight, &work[1], n, &q[*ilo * q_dim1 + 1], ldq/*, (
		ftnlen)3*/);
    }
/*     Update A(istartm:ilo-1,ilo:ilo+ns-1) and B(istartm:ilo-1,ilo:ilo+ns-1) */
/*     from the right with Zc(1:ns,1:ns) */
    sheight = *ilo - 1 - istartm + 1;
    swidth = ns;
    if (sheight > 0) {
	sgemm_("N", "N", &sheight, &swidth, &swidth, &c_b5, &a[istartm + *ilo 
		* a_dim1], lda, &zc[zc_offset], ldzc, &c_b4, &work[1], &
		sheight/*, (ftnlen)1, (ftnlen)1*/);
	slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[istartm + *
		ilo * a_dim1], lda/*, (ftnlen)3*/);
	sgemm_("N", "N", &sheight, &swidth, &swidth, &c_b5, &b[istartm + *ilo 
		* b_dim1], ldb, &zc[zc_offset], ldzc, &c_b4, &work[1], &
		sheight/*, (ftnlen)1, (ftnlen)1*/);
	slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[istartm + *
		ilo * b_dim1], ldb/*, (ftnlen)3*/);
    }
    if (*ilz) {
	sgemm_("N", "N", n, &swidth, &swidth, &c_b5, &z__[*ilo * z_dim1 + 1], 
		ldz, &zc[zc_offset], ldzc, &c_b4, &work[1], n/*, (ftnlen)1, (
		ftnlen)1*/);
	slacpy_("ALL", n, &swidth, &work[1], n, &z__[*ilo * z_dim1 + 1], ldz
			/*, (ftnlen)3*/);
    }
/*     The following block chases the shifts down to the bottom */
/*     right block. If possible, a shift is moved down npos */
/*     positions at a time */
    k = *ilo;
    while(k < *ihi - ns) {
/* Computing MIN */
	i__1 = *ihi - ns - k;
	np = f2cmin(i__1,npos);
/*        Size of the near-the-diagonal block */
	nblock = ns + np;
/*        istartb points to the first row we will be updating */
	istartb = k + 1;
/*        istopb points to the last column we will be updating */
	istopb = k + nblock - 1;
	i__1 = ns + np;
	i__2 = ns + np;
	slaset_("FULL", &i__1, &i__2, &c_b4, &c_b5, &qc[qc_offset], ldqc
			/*, (ftnlen)4*/);
	i__1 = ns + np;
	i__2 = ns + np;
	slaset_("FULL", &i__1, &i__2, &c_b4, &c_b5, &zc[zc_offset], ldzc
			/*, (ftnlen)4*/);
/*        Near the diagonal shift chase */
	for (i__ = ns - 1; i__ >= 0; i__ += -2) {
	    i__1 = np - 1;
	    for (j = 0; j <= i__1; ++j) {
/*              Move down the block with index k+i+j-1, updating */
/*              the (ns+np x ns+np) block: */
/*              (k:k+ns+np,k:k+ns+np-1) */
		i__2 = k + i__ + j - 1;
		i__3 = k + 1;
		slaqz2_(&c_true, &c_true, &i__2, &istartb, &istopb, ihi, &a[
			a_offset], lda, &b[b_offset], ldb, &nblock, &i__3, &
			qc[qc_offset], ldqc, &nblock, &k, &zc[zc_offset], 
			ldzc);
	    }
	}
/*        Update rest of the pencil */
/*        Update A(k+1:k+ns+np, k+ns+np:istopm) and */
/*        B(k+1:k+ns+np, k+ns+np:istopm) */
/*        from the left with Qc(1:ns+np,1:ns+np)' */
	sheight = ns + np;
	swidth = istopm - (k + ns + np) + 1;
	if (swidth > 0) {
	    sgemm_("T", "N", &sheight, &swidth, &sheight, &c_b5, &qc[
		    qc_offset], ldqc, &a[k + 1 + (k + ns + np) * a_dim1], lda,
		     &c_b4, &work[1], &sheight /*, (ftnlen)1, (ftnlen)1*/);
	    slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[k + 1 + (
		    k + ns + np) * a_dim1], lda /*, (ftnlen)3*/);
	    sgemm_("T", "N", &sheight, &swidth, &sheight, &c_b5, &qc[
		    qc_offset], ldqc, &b[k + 1 + (k + ns + np) * b_dim1], ldb,
		     &c_b4, &work[1], &sheight/*, (ftnlen)1, (ftnlen)1*/);
	    slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[k + 1 + (
		    k + ns + np) * b_dim1], ldb/*, (ftnlen)3*/);
	}
	if (*ilq) {
	    sgemm_("N", "N", n, &nblock, &nblock, &c_b5, &q[(k + 1) * q_dim1 
		    + 1], ldq, &qc[qc_offset], ldqc, &c_b4, &work[1], n
		    /*, (ftnlen)1, (ftnlen)1*/);
	    slacpy_("ALL", n, &nblock, &work[1], n, &q[(k + 1) * q_dim1 + 1], 
		    ldq /*, (ftnlen)3*/ );
	}
/*        Update A(istartm:k,k:k+ns+npos-1) and B(istartm:k,k:k+ns+npos-1) */
/*        from the right with Zc(1:ns+np,1:ns+np) */
	sheight = k - istartm + 1;
	swidth = nblock;
	if (sheight > 0) {
	    sgemm_("N", "N", &sheight, &swidth, &swidth, &c_b5, &a[istartm + 
		    k * a_dim1], lda, &zc[zc_offset], ldzc, &c_b4, &work[1], &
		    sheight /*, (ftnlen)1, (ftnlen)1*/ );
	    slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[istartm 
		    + k * a_dim1], lda /*, (ftnlen)3 */);
	    sgemm_("N", "N", &sheight, &swidth, &swidth, &c_b5, &b[istartm + 
		    k * b_dim1], ldb, &zc[zc_offset], ldzc, &c_b4, &work[1], &
		    sheight /*, (ftnlen)1, (ftnlen)1*/ );
	    slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[istartm 
		    + k * b_dim1], ldb /*, (ftnlen)3*/ );
	}
	if (*ilz) {
	    sgemm_("N", "N", n, &nblock, &nblock, &c_b5, &z__[k * z_dim1 + 1],
		     ldz, &zc[zc_offset], ldzc, &c_b4, &work[1], n
		     /* , (ftnlen)1,(ftnlen)1*/ );
	    slacpy_("ALL", n, &nblock, &work[1], n, &z__[k * z_dim1 + 1], ldz
			    /* ,(ftnlen)3*/ );
	}
	k += np;
    }
/*     The following block removes the shifts from the bottom right corner */
/*     one by one. Updates are initially applied to A(ihi-ns+1:ihi,ihi-ns:ihi). */
    slaset_("FULL", &ns, &ns, &c_b4, &c_b5, &qc[qc_offset], ldqc /*, (ftnlen)4*/ );
    i__1 = ns + 1;
    i__2 = ns + 1;
    slaset_("FULL", &i__1, &i__2, &c_b4, &c_b5, &zc[zc_offset], ldzc /*, (ftnlen)4*/ );
/*     istartb points to the first row we will be updating */
    istartb = *ihi - ns + 1;
/*     istopb points to the last column we will be updating */
    istopb = *ihi;
    i__1 = ns;
    for (i__ = 1; i__ <= i__1; i__ += 2) {
/*        Chase the shift down to the bottom right corner */
	i__2 = *ihi - 2;
	for (ishift = *ihi - i__ - 1; ishift <= i__2; ++ishift) {
	    i__3 = *ihi - ns + 1;
	    i__4 = ns + 1;
	    i__5 = *ihi - ns;
	    slaqz2_(&c_true, &c_true, &ishift, &istartb, &istopb, ihi, &a[
		    a_offset], lda, &b[b_offset], ldb, &ns, &i__3, &qc[
		    qc_offset], ldqc, &i__4, &i__5, &zc[zc_offset], ldzc);
	}
    }
/*     Update rest of the pencil */
/*     Update A(ihi-ns+1:ihi, ihi+1:istopm) */
/*     from the left with Qc(1:ns,1:ns)' */
    sheight = ns;
    swidth = istopm - (*ihi + 1) + 1;
    if (swidth > 0) {
	sgemm_("T", "N", &sheight, &swidth, &sheight, &c_b5, &qc[qc_offset], 
		ldqc, &a[*ihi - ns + 1 + (*ihi + 1) * a_dim1], lda, &c_b4, &
		work[1], &sheight /*, (ftnlen)1, (ftnlen)1*/ );
	slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[*ihi - ns + 
		1 + (*ihi + 1) * a_dim1], lda /*, (ftnlen)3*/);
	sgemm_("T", "N", &sheight, &swidth, &sheight, &c_b5, &qc[qc_offset], 
		ldqc, &b[*ihi - ns + 1 + (*ihi + 1) * b_dim1], ldb, &c_b4, &
		work[1], &sheight /*, (ftnlen)1, (ftnlen)1*/);
	slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[*ihi - ns + 
		1 + (*ihi + 1) * b_dim1], ldb /*, (ftnlen)3*/);
    }
    if (*ilq) {
	sgemm_("N", "N", n, &ns, &ns, &c_b5, &q[(*ihi - ns + 1) * q_dim1 + 1],
		 ldq, &qc[qc_offset], ldqc, &c_b4, &work[1], n
		 /* , (ftnlen)1, (ftnlen)1 */ );
	slacpy_("ALL", n, &ns, &work[1], n, &q[(*ihi - ns + 1) * q_dim1 + 1], 
		ldq /*, (ftnlen)3*/ );
    }
/*     Update A(istartm:ihi-ns,ihi-ns:ihi) */
/*     from the right with Zc(1:ns+1,1:ns+1) */
    sheight = *ihi - ns - istartm + 1;
    swidth = ns + 1;
    if (sheight > 0) {
	sgemm_("N", "N", &sheight, &swidth, &swidth, &c_b5, &a[istartm + (*
		ihi - ns) * a_dim1], lda, &zc[zc_offset], ldzc, &c_b4, &work[
		1], &sheight /*, (ftnlen)1, (ftnlen)1*/ );
	slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &a[istartm + (*
		ihi - ns) * a_dim1], lda /*, (ftnlen)3*/);
	sgemm_("N", "N", &sheight, &swidth, &swidth, &c_b5, &b[istartm + (*
		ihi - ns) * b_dim1], ldb, &zc[zc_offset], ldzc, &c_b4, &work[
		1], &sheight /*, (ftnlen)1, (ftnlen)1*/);
	slacpy_("ALL", &sheight, &swidth, &work[1], &sheight, &b[istartm + (*
		ihi - ns) * b_dim1], ldb /*, (ftnlen)3*/);
    }
    if (*ilz) {
	i__1 = ns + 1;
	i__2 = ns + 1;
	sgemm_("N", "N", n, &i__1, &i__2, &c_b5, &z__[(*ihi - ns) * z_dim1 + 
		1], ldz, &zc[zc_offset], ldzc, &c_b4, &work[1], n
		/* , (ftnlen)1, (ftnlen)1*/ );
	i__1 = ns + 1;
	slacpy_("ALL", n, &i__1, &work[1], n, &z__[(*ihi - ns) * z_dim1 + 1], 
		ldz /* , (ftnlen)3*/ );
    }
    return;
} /* slaqz4_ */

