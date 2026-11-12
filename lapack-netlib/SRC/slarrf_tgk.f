*> \brief \b SLARRF_TGK computes a shifted child representation for a cluster (Golub-Kahan variant of SLARRF).
*>
*> \ingroup bdsvdmr3
*>
*> \par Purpose:
*> =============
*>
*> \verbatim
*> computes a shifted child representation for a cluster (Golub-Kahan variant of SLARRF).
*> All routines in this file are auxiliary to SBDSVDMR3 and derive
*> from the LAPACK 3.0 (1999) MRRR kernels SLARRV, SLAR1V, SLARRB and
*> SLARRF rather than from the current SSTEMR kernels.
*> \endverbatim
*
      SUBROUTINE SLARRF_TGK( REP, N, D, L, LD, LLD, IFIRST, ILAST, W,
     $                   SIGMA, DPLUS, LPLUS, WORK, ISEED, INFO )
*
*  -- New auxiliary routine for bidiagonal SVD via TGK-rooted MR^3 --
*     Based on LAPACK SLARRF (Dhillon & Marques, Nov 11 2003).
*     Derived from the LAPACK 3.0 SLARRF; this is a separate copy that adds
*     a representation switch so the parent may be supplied either as an
*     L D L^T factorization (REP='L', identical to the original) or as a
*     symmetric tridiagonal matrix -- in particular the Tridiagonal
*     Golub-Kahan (TGK) matrix -- with REP='T'.
*
*     .. Scalar Arguments ..
      CHARACTER          REP
      INTEGER            IFIRST, ILAST, INFO, N
      REAL   SIGMA
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      REAL   D( * ), DPLUS( * ), L( * ), LD( * ), LLD( * ),
     $                   LPLUS( * ), W( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  Given a parent representation M and a cluster of close eigenvalues
*  (in a relative measure) W( IFIRST ), ..., W( ILAST ), SLARRF_TGK finds
*  a new relatively robust representation
*       M - SIGMA I = L(+) D(+) L(+)^T
*  such that at least one of the eigenvalues of L(+) D(+) L(+)^T is
*  relatively isolated.  The output L(+) D(+) L(+)^T is an ordinary
*  L D L^T representation regardless of REP, so every descendant node in
*  the representation tree can be processed by the unmodified LAPACK 3.0-lineage
*  routines.
*
*  Arguments
*  =========
*
*  REP     (input) CHARACTER*1
*          = 'L':  the parent M is the L D L^T representation defined by
*                  D, L (and LD = L(i)*D(i), LLD = L(i)^2*D(i)).  Behaviour
*                  is identical to the original SLARRF.
*          = 'T':  the parent M is a symmetric tridiagonal matrix (e.g.
*                  the TGK matrix).  D holds the N diagonal entries (all
*                  zero for TGK, kept general) and L holds the N-1
*                  off-diagonal entries beta(i).  LD and LLD are NOT
*                  referenced.
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  D       (input) REAL array, dimension (N)
*          REP='L': diagonal of D.   REP='T': diagonal of M.
*
*  L       (input) REAL array, dimension (N-1)
*          REP='L': subdiagonal of the unit bidiagonal L.
*          REP='T': off-diagonal entries beta(i) of M.
*
*  LD      (input) REAL array, dimension (N-1)
*          REP='L': L(i)*D(i).   REP='T': not referenced.
*
*  LLD     (input) REAL array, dimension (N-1)
*          REP='L': L(i)*L(i)*D(i).   REP='T': not referenced.
*
*  IFIRST  (input) INTEGER   Index of the first eigenvalue in the cluster.
*  ILAST   (input) INTEGER   Index of the last  eigenvalue in the cluster.
*
*  W       (input) REAL array, dimension (N)
*          Eigenvalues of M in ascending order; W(IFIRST..ILAST) is the
*          cluster.  (For the TGK root these are the positive eigenvalues
*          == singular values; the caller passes them in ascending order.)
*
*  SIGMA   (output) REAL  The shift used to form L(+)D(+)L(+)^T.
*
*  DPLUS   (output) REAL array, dimension (N)   Diagonal D(+).
*  LPLUS   (output) REAL array, dimension (N-1) Subdiag L(+).
*
*  WORK    (workspace) REAL array, dimension (2*N)
*
*  INFO    (output) INTEGER   = 0 successful exit; < 0 illegal argument.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL   ZERO, ONE, TWO, PERTK
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0,
     $                   PERTK = 4.0E0 )
      INTEGER            MAXTRY
      PARAMETER          ( MAXTRY = 10 )
*     DOPERT enables the Sec.5.4 random relative perturbation of the child
*     RRR.  It is applied SELECTIVELY -- only to a child whose cluster has
*     a near-degenerate pair, W(i+1)-W(i) < PGFAC*eps*W(i) (PGFAC=100).
*     The per-entry size is the full PERTK*eps: empirically the magnitude
*     cannot be reduced (PERTK*eps/N and PERTK*eps/sqrt(N) both fail -- the
*     overflow NaN returns) because breaking the near-zero pivot
*     cancellation that destroys the eigenvector (DPLUS ~ 1e2*eps, see
*     below) needs an O(eps) per-entry jolt.  The gate is what keeps this
*     from over-perturbing moderately-clustered RRRs.
      LOGICAL            DOPERT
      REAL   PGFAC
      PARAMETER          ( DOPERT = .FALSE., PGFAC = 1.0E2 )
*     ..
*     .. Local Scalars ..
      LOGICAL            TGK, DEGEN
      INTEGER            I, NTRY
      REAL   DELTA, DMAX1, DMAX2, EPS, S, TMP, PERT
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL   SLAMCH
      EXTERNAL           LSAME, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SLARNV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      IF( .NOT.( LSAME( REP, 'L' ) .OR. LSAME( REP, 'T' ) ) ) THEN
         INFO = -1
         RETURN
      END IF
      TGK = LSAME( REP, 'T' )
      EPS = SLAMCH( 'Precision' )
      SIGMA = W( IFIRST )
      DELTA = TWO * EPS
      NTRY = 0
*
*     Try the representation based at the left end of the cluster,
*     SIGMA = W( IFIRST ).  DMAX1 monitors the element growth.
*
   10 CONTINUE
      IF( TGK ) THEN
*        Tridiagonal parent:  M - SIGMA I = L(+) D(+) L(+)^T directly.
         DPLUS( 1 ) = D( 1 ) - SIGMA
         DMAX1 = ABS( DPLUS( 1 ) )
         DO 20 I = 1, N - 1
            LPLUS( I ) = L( I ) / DPLUS( I )
            DPLUS( I+1 ) = D( I+1 ) - SIGMA - LPLUS( I )*L( I )
            DMAX1 = MAX( DMAX1, ABS( DPLUS( I+1 ) ) )
   20    CONTINUE
      ELSE
*        L D L^T parent (original SLARRF recurrence).
         S = -SIGMA
         DPLUS( 1 ) = D( 1 ) + S
         DMAX1 = ABS( DPLUS( 1 ) )
         DO 25 I = 1, N - 1
            LPLUS( I ) = LD( I ) / DPLUS( I )
            S = S*LPLUS( I )*L( I ) - SIGMA
            DPLUS( I+1 ) = D( I+1 ) + S
            DMAX1 = MAX( DMAX1, ABS( DPLUS( I+1 ) ) )
   25    CONTINUE
      END IF
      IF( DMAX1*ZERO.NE.ZERO ) THEN
*        A non-finite pivot (Inf or NaN) appeared -- e.g. a zero TGK
*        off-diagonal meeting a zero pivot at SIGMA=0.  (The DMAX1*0.NE.0
*        test is true exactly for Inf/NaN.)  An Inf child would later
*        hang SLARRB's bisection, so we reject it here.  Nudge SIGMA off
*        the eigenvalue and retry,
*        but only a bounded number of times so we can never loop forever
*        (the additive ONE term lets the nudge escape SIGMA=0; the matrix
*        is pre-scaled to O(1) norm by SBDSVDMR3).
         NTRY = NTRY + 1
         IF( NTRY.GT.MAXTRY ) THEN
            INFO = 1
            RETURN
         END IF
         SIGMA = SIGMA - ( ABS( SIGMA ) + ONE )*DELTA
         DELTA = TWO*DELTA
         GO TO 10
      END IF
*
*     Try the representation based at the right end of the cluster,
*     TMP = W( ILAST ), storing it in WORK; DMAX2 monitors its growth.
*
      TMP = W( ILAST )
      DELTA = TWO * EPS
      NTRY = 0
   30 CONTINUE
      IF( TGK ) THEN
         WORK( 1 ) = D( 1 ) - TMP
         DMAX2 = ABS( WORK( 1 ) )
         DO 40 I = 1, N - 1
            WORK( N+I ) = L( I ) / WORK( I )
            WORK( I+1 ) = D( I+1 ) - TMP - WORK( N+I )*L( I )
            DMAX2 = MAX( DMAX2, ABS( WORK( I+1 ) ) )
   40    CONTINUE
      ELSE
         S = -TMP
         WORK( 1 ) = D( 1 ) + S
         DMAX2 = ABS( WORK( 1 ) )
         DO 45 I = 1, N - 1
            WORK( N+I ) = LD( I ) / WORK( I )
            S = S*WORK( N+I )*L( I ) - TMP
            WORK( I+1 ) = D( I+1 ) + S
            DMAX2 = MAX( DMAX2, ABS( WORK( I+1 ) ) )
   45    CONTINUE
      END IF
      IF( DMAX2*ZERO.NE.ZERO ) THEN
         NTRY = NTRY + 1
         IF( NTRY.GT.MAXTRY ) THEN
*           The right-end shift never gave finite growth; fall back to
*           the (finite) left-end representation found above.
            GO TO 50
         END IF
         TMP = TMP + ( ABS( TMP ) + ONE )*DELTA
         DELTA = TWO*DELTA
         GO TO 30
      END IF
*
*     Keep whichever shift produced the smaller element growth.
*
      IF( DMAX2.LT.DMAX1 ) THEN
         SIGMA = TMP
         CALL SCOPY( N, WORK, 1, DPLUS, 1 )
         CALL SCOPY( N-1, WORK(N+1), 1, LPLUS, 1 )
      END IF
*
   50 CONTINUE
*
*     Random relative perturbation of the child RRR (Dhillon, Parlett &
*     Voemel 2005, eqs. 5.1-5.2):
*        DPLUS(i) := DPLUS(i)*(1 + mu_i *PERTK*eps),
*        LPLUS(j) := LPLUS(j)*(1 + nu_j *PERTK*eps),
*     with mu,nu uniform on [-1,1] (SLARNV, IDIST=2) and PERTK=4.
*
*     WHY this is needed (and load-bearing for the TGK).  A tight near-
*     degenerate pair of singular values (e.g. a glued-Wilkinson Wilkinson
*     pair, relative gap ~eps) forces a child RRR whose shift sits almost
*     exactly on an eigenvalue, giving a near-zero pivot / large element
*     growth (DPLUS as small as ~1e2*eps next to ~1/eps).  The downstream
*     SLAR1V twisted factorization then drives lambda onto that eigenvalue,
*     the pivot collapses to ~1e-30, the back-substitution overflows, and
*     the eigenvector comes back as NaN.  A relative perturbation of the
*     child factors breaks the exact pivot cancellation: the perturbed RRR
*     has no eigenvalue closer than ~eps*||M|| to the shift, so the twist
*     gamma cannot collapse and the vector stays finite.
*
*     APPLIED SELECTIVELY -- only to a child whose cluster contains a near-
*     degenerate consecutive pair, W(i+1)-W(i) < PGFAC*eps*|W(i)| with
*     PGFAC=100.  Perturbing every child instead (no gate) needlessly
*     jitters moderately-clustered RRRs and degrades their orthogonality;
*     the gate confines the perturbation to the clusters that actually need
*     it.  The size cannot be lowered below the full PERTK*eps per entry --
*     PERTK*eps/N and PERTK*eps/sqrt(N) were both tried and the overflow
*     NaN returned -- because breaking the cancellation needs an O(eps)
*     per-entry jolt; this does scramble genuinely sub-eps gaps (strong-
*     glue cases land at ~1e-8..1e-6 orthogonality), the unavoidable price.
*     ISEED is threaded by the caller (fixed seed) for reproducibility.
*
      DEGEN = .FALSE.
      IF( DOPERT ) THEN
         DO 55 I = IFIRST, ILAST - 1
            IF( W( I+1 )-W( I ).LT.PGFAC*EPS*ABS( W( I ) ) )
     $         DEGEN = .TRUE.
   55    CONTINUE
      END IF
      IF( DEGEN ) THEN
         PERT = PERTK*EPS
         CALL SLARNV( 2, ISEED, N, WORK )
         DO 60 I = 1, N
            DPLUS( I ) = DPLUS( I )*( ONE + WORK( I )*PERT )
   60    CONTINUE
         CALL SLARNV( 2, ISEED, N-1, WORK )
         DO 70 I = 1, N - 1
            LPLUS( I ) = LPLUS( I )*( ONE + WORK( I )*PERT )
   70    CONTINUE
      END IF
      RETURN
*
*     End of SLARRF_TGK
*
      END
