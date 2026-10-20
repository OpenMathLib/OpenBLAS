

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* > \brief \b DLARRF_TGK computes a shifted child representation for a cluster (Golub-Kahan variant of DLARRF
). */
/* > */
/* > \ingroup bdsvdmr3 */
/* > */
/* > \par Purpose: */
/* > ============= */
/* > */
/* > \verbatim */
/* > computes a shifted child representation for a cluster (Golub-Kahan variant of DLARRF). */
/* > All routines in this file are auxiliary to DBDSVDMR3 and derive */
/* > from the LAPACK 3.0 (1999) MRRR kernels DLARRV, DLAR1V, DLARRB and */
/* > DLARRF rather than from the current DSTEMR kernels. */
/* > \endverbatim */

/* Subroutine */ void dlarrf_tgk__(char *rep, integer *n, doublereal *d__, 
	doublereal *l, doublereal *ld, doublereal *lld, integer *ifirst, 
	integer *ilast, doublereal *w, doublereal *sigma, doublereal *dplus, 
	doublereal *lplus, doublereal *work, integer *iseed, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    integer i__;
    doublereal s;
    logical tgk;
    doublereal eps, tmp, pert;
    integer ntry;
    doublereal dmax1, dmax2;
    logical degen;
    doublereal delta;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ void dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    extern doublereal dlamch_(char *);
    extern /* Subroutine */ void dlarnv_(integer *, integer *, integer *, 
	    doublereal *);


/*  -- New auxiliary routine for bidiagonal SVD via TGK-rooted MR^3 -- */
/*     Based on LAPACK DLARRF (Dhillon & Marques, Nov 11 2003). */
/*     Derived from the LAPACK 3.0 DLARRF; this is a separate copy that adds */
/*     a representation switch so the parent may be supplied either as an */
/*     L D L^T factorization (REP='L', identical to the original) or as a */
/*     symmetric tridiagonal matrix -- in particular the Tridiagonal */
/*     Golub-Kahan (TGK) matrix -- with REP='T'. */


/*  Purpose */
/*  ======= */

/*  Given a parent representation M and a cluster of close eigenvalues */
/*  (in a relative measure) W( IFIRST ), ..., W( ILAST ), DLARRF_TGK finds */
/*  a new relatively robust representation */
/*       M - SIGMA I = L(+) D(+) L(+)^T */
/*  such that at least one of the eigenvalues of L(+) D(+) L(+)^T is */
/*  relatively isolated.  The output L(+) D(+) L(+)^T is an ordinary */
/*  L D L^T representation regardless of REP, so every descendant node in */
/*  the representation tree can be processed by the unmodified LAPACK 3.0-lineage */
/*  routines. */

/*  Arguments */
/*  ========= */

/*  REP     (input) CHARACTER*1 */
/*          = 'L':  the parent M is the L D L^T representation defined by */
/*                  D, L (and LD = L(i)*D(i), LLD = L(i)^2*D(i)).  Behaviour */
/*                  is identical to the original DLARRF. */
/*          = 'T':  the parent M is a symmetric tridiagonal matrix (e.g. */
/*                  the TGK matrix).  D holds the N diagonal entries (all */
/*                  zero for TGK, kept general) and L holds the N-1 */
/*                  off-diagonal entries beta(i).  LD and LLD are NOT */
/*                  referenced. */

/*  N       (input) INTEGER */
/*          The order of the matrix. */

/*  D       (input) DOUBLE PRECISION array, dimension (N) */
/*          REP='L': diagonal of D.   REP='T': diagonal of M. */

/*  L       (input) DOUBLE PRECISION array, dimension (N-1) */
/*          REP='L': subdiagonal of the unit bidiagonal L. */
/*          REP='T': off-diagonal entries beta(i) of M. */

/*  LD      (input) DOUBLE PRECISION array, dimension (N-1) */
/*          REP='L': L(i)*D(i).   REP='T': not referenced. */

/*  LLD     (input) DOUBLE PRECISION array, dimension (N-1) */
/*          REP='L': L(i)*L(i)*D(i).   REP='T': not referenced. */

/*  IFIRST  (input) INTEGER   Index of the first eigenvalue in the cluster. */
/*  ILAST   (input) INTEGER   Index of the last  eigenvalue in the cluster. */

/*  W       (input) DOUBLE PRECISION array, dimension (N) */
/*          Eigenvalues of M in ascending order; W(IFIRST..ILAST) is the */
/*          cluster.  (For the TGK root these are the positive eigenvalues */
/*          == singular values; the caller passes them in ascending order.) */

/*  SIGMA   (output) DOUBLE PRECISION  The shift used to form L(+)D(+)L(+)^T. */

/*  DPLUS   (output) DOUBLE PRECISION array, dimension (N)   Diagonal D(+). */
/*  LPLUS   (output) DOUBLE PRECISION array, dimension (N-1) Subdiag L(+). */

/*  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N) */

/*  INFO    (output) INTEGER   = 0 successful exit; < 0 illegal argument. */

/*  ===================================================================== */

/*     DOPERT enables the Sec.5.4 random relative perturbation of the child */
/*     RRR.  It is applied SELECTIVELY -- only to a child whose cluster has */
/*     a near-degenerate pair, W(i+1)-W(i) < PGFAC*eps*W(i) (PGFAC=100). */
/*     The per-entry size is the full PERTK*eps: empirically the magnitude */
/*     cannot be reduced (PERTK*eps/N and PERTK*eps/sqrt(N) both fail -- the */
/*     overflow NaN returns) because breaking the near-zero pivot */
/*     cancellation that destroys the eigenvector (DPLUS ~ 1e2*eps, see */
/*     below) needs an O(eps) per-entry jolt.  The gate is what keeps this */
/*     from over-perturbing moderately-clustered RRRs. */

    /* Parameter adjustments */
    --iseed;
    --work;
    --lplus;
    --dplus;
    --w;
    --lld;
    --ld;
    --l;
    --d__;

    /* Function Body */
    *info = 0;
    if (! (lsame_(rep, "L") || lsame_(rep, "T"))) {
	*info = -1;
	return;
    }
    tgk = lsame_(rep, "T");
    eps = dlamch_("Precision");
    *sigma = w[*ifirst];
    delta = eps * 2.;
    ntry = 0;

/*     Try the representation based at the left end of the cluster, */
/*     SIGMA = W( IFIRST ).  DMAX1 monitors the element growth. */

L10:
    if (tgk) {
/*        Tridiagonal parent:  M - SIGMA I = L(+) D(+) L(+)^T directly. */
	dplus[1] = d__[1] - *sigma;
	dmax1 = abs(dplus[1]);
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    lplus[i__] = l[i__] / dplus[i__];
	    dplus[i__ + 1] = d__[i__ + 1] - *sigma - lplus[i__] * l[i__];
/* Computing MAX */
	    d__2 = dmax1, d__3 = (d__1 = dplus[i__ + 1], abs(d__1));
	    dmax1 = f2cmax(d__2,d__3);
/* L20: */
	}
    } else {
/*        L D L^T parent (original DLARRF recurrence). */
	s = -(*sigma);
	dplus[1] = d__[1] + s;
	dmax1 = abs(dplus[1]);
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    lplus[i__] = ld[i__] / dplus[i__];
	    s = s * lplus[i__] * l[i__] - *sigma;
	    dplus[i__ + 1] = d__[i__ + 1] + s;
/* Computing MAX */
	    d__2 = dmax1, d__3 = (d__1 = dplus[i__ + 1], abs(d__1));
	    dmax1 = f2cmax(d__2,d__3);
/* L25: */
	}
    }
    if (dmax1 * 0. != 0.) {
/*        A non-finite pivot (Inf or NaN) appeared -- e.g. a zero TGK */
/*        off-diagonal meeting a zero pivot at SIGMA=0.  (The DMAX1*0.NE.0 */
/*        test is true exactly for Inf/NaN.)  An Inf child would later */
/*        hang DLARRB's bisection, so we reject it here.  Nudge SIGMA off */
/*        the eigenvalue and retry, */
/*        but only a bounded number of times so we can never loop forever */
/*        (the additive ONE term lets the nudge escape SIGMA=0; the matrix */
/*        is pre-scaled to O(1) norm by DBDSVDMR3). */
	++ntry;
	if (ntry > 10) {
	    *info = 1;
	    return;
	}
	*sigma -= (abs(*sigma) + 1.) * delta;
	delta *= 2.;
	goto L10;
    }

/*     Try the representation based at the right end of the cluster, */
/*     TMP = W( ILAST ), storing it in WORK; DMAX2 monitors its growth. */

    tmp = w[*ilast];
    delta = eps * 2.;
    ntry = 0;
L30:
    if (tgk) {
	work[1] = d__[1] - tmp;
	dmax2 = abs(work[1]);
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    work[*n + i__] = l[i__] / work[i__];
	    work[i__ + 1] = d__[i__ + 1] - tmp - work[*n + i__] * l[i__];
/* Computing MAX */
	    d__2 = dmax2, d__3 = (d__1 = work[i__ + 1], abs(d__1));
	    dmax2 = f2cmax(d__2,d__3);
/* L40: */
	}
    } else {
	s = -tmp;
	work[1] = d__[1] + s;
	dmax2 = abs(work[1]);
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    work[*n + i__] = ld[i__] / work[i__];
	    s = s * work[*n + i__] * l[i__] - tmp;
	    work[i__ + 1] = d__[i__ + 1] + s;
/* Computing MAX */
	    d__2 = dmax2, d__3 = (d__1 = work[i__ + 1], abs(d__1));
	    dmax2 = f2cmax(d__2,d__3);
/* L45: */
	}
    }
    if (dmax2 * 0. != 0.) {
	++ntry;
	if (ntry > 10) {
/*           The right-end shift never gave finite growth; fall back to */
/*           the (finite) left-end representation found above. */
	    goto L50;
	}
	tmp += (abs(tmp) + 1.) * delta;
	delta *= 2.;
	goto L30;
    }

/*     Keep whichever shift produced the smaller element growth. */

    if (dmax2 < dmax1) {
	*sigma = tmp;
	dcopy_(n, &work[1], &c__1, &dplus[1], &c__1);
	i__1 = *n - 1;
	dcopy_(&i__1, &work[*n + 1], &c__1, &lplus[1], &c__1);
    }

L50:

/*     Random relative perturbation of the child RRR (Dhillon, Parlett & */
/*     Voemel 2005, eqs. 5.1-5.2): */
/*        DPLUS(i) := DPLUS(i)*(1 + mu_i *PERTK*eps), */
/*        LPLUS(j) := LPLUS(j)*(1 + nu_j *PERTK*eps), */
/*     with mu,nu uniform on [-1,1] (DLARNV, IDIST=2) and PERTK=4. */

/*     WHY this is needed (and load-bearing for the TGK).  A tight near- */
/*     degenerate pair of singular values (e.g. a glued-Wilkinson Wilkinson */
/*     pair, relative gap ~eps) forces a child RRR whose shift sits almost */
/*     exactly on an eigenvalue, giving a near-zero pivot / large element */
/*     growth (DPLUS as small as ~1e2*eps next to ~1/eps).  The downstream */
/*     DLAR1V twisted factorization then drives lambda onto that eigenvalue, */
/*     the pivot collapses to ~1e-30, the back-substitution overflows, and */
/*     the eigenvector comes back as NaN.  A relative perturbation of the */
/*     child factors breaks the exact pivot cancellation: the perturbed RRR */
/*     has no eigenvalue closer than ~eps*||M|| to the shift, so the twist */
/*     gamma cannot collapse and the vector stays finite. */

/*     APPLIED SELECTIVELY -- only to a child whose cluster contains a near- */
/*     degenerate consecutive pair, W(i+1)-W(i) < PGFAC*eps*|W(i)| with */
/*     PGFAC=100.  Perturbing every child instead (no gate) needlessly */
/*     jitters moderately-clustered RRRs and degrades their orthogonality; */
/*     the gate confines the perturbation to the clusters that actually need */
/*     it.  The size cannot be lowered below the full PERTK*eps per entry -- */
/*     PERTK*eps/N and PERTK*eps/sqrt(N) were both tried and the overflow */
/*     NaN returned -- because breaking the cancellation needs an O(eps) */
/*     per-entry jolt; this does scramble genuinely sub-eps gaps (strong- */
/*     glue cases land at ~1e-8..1e-6 orthogonality), the unavoidable price. */
/*     ISEED is threaded by the caller (fixed seed) for reproducibility. */

    degen = FALSE_;
    if (FALSE_) {
	i__1 = *ilast - 1;
	for (i__ = *ifirst; i__ <= i__1; ++i__) {
	    if (w[i__ + 1] - w[i__] < eps * 100. * (d__1 = w[i__], abs(d__1)))
		     {
		degen = TRUE_;
	    }
/* L55: */
	}
    }
    if (degen) {
	pert = eps * 4.;
	dlarnv_(&c__2, &iseed[1], n, &work[1]);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dplus[i__] *= work[i__] * pert + 1.;
/* L60: */
	}
	i__1 = *n - 1;
	dlarnv_(&c__2, &iseed[1], &i__1, &work[1]);
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    lplus[i__] *= work[i__] * pert + 1.;
/* L70: */
	}
    }
    return;

/*     End of DLARRF_TGK */

} /* dlarrf_tgk__ */

