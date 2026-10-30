*> \brief \b DLARRB_TGK bisection refinement of eigenvalues of an LDL^T representation (Golub-Kahan variant of DLARRB); also DLARRBTGK_STURM4 (four-lane Sturm count) and DLARRBTGK_STURM_SLOW (reference Sturm count).
*>
*> \ingroup bdsvdmr3
*>
*> \par Purpose:
*> =============
*>
*> \verbatim
*> bisection refinement of eigenvalues of an LDL^T representation (Golub-Kahan variant of DLARRB); also DLARRBTGK_STURM4 (four-lane Sturm count) and DLARRBTGK_STURM_SLOW (reference Sturm count).
*> All routines in this file are auxiliary to DBDSVDMR3 and derive
*> from the LAPACK 3.0 (1999) MRRR kernels DLARRV, DLAR1V, DLARRB and
*> DLARRF rather than from the current DSTEMR kernels.
*> \endverbatim
*
      SUBROUTINE DLARRB_TGK( N, D, L, LD, LLD, IFIRST, ILAST, RTOL1,
     $                   RTOL2, OFFSET, W, WGAP, WERR, WORK, IWORK,
     $                   INFO )
*
*  -- LAPACK auxiliary routine --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     Derived from the LAPACK 3.0 DLARRB (1999); TGK adaptation 2026
*
*     .. Scalar Arguments ..
      INTEGER            IFIRST, ILAST, INFO, N, OFFSET
      DOUBLE PRECISION   RTOL1, RTOL2
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   D( * ), L( * ), LD( * ), LLD( * ), W( * ),
     $                   WERR( * ), WGAP( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  Given the relatively robust representation(RRR) L D L^T, DLARRB
*  does "limited" bisection to refine the eigenvalues of L D L^T,
*  W( IFIRST-OFFSET ) thru' W( ILAST-OFFSET ), to more accuracy. Initial
*  guesses for these eigenvalues are input in W, the corresponding estimate
*  of the error in these guesses and their gaps are input in WERR
*  and WGAP, respectively. During bisection, intervals
*  [left, right] are maintained by storing their mid-points and
*  semi-widths in the arrays W and WERR respectively.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          The N diagonal elements of the diagonal matrix D.
*
*  L       (input) DOUBLE PRECISION array, dimension (N-1)
*          The (N-1) subdiagonal elements of the unit bidiagonal matrix L.
*
*  LD      (input) DOUBLE PRECISION array, dimension (N-1)
*          The (N-1) elements L(i)*D(i).
*
*  LLD     (input) DOUBLE PRECISION array, dimension (N-1)
*          The (N-1) elements L(i)*L(i)*D(i).
*
*  IFIRST  (input) INTEGER
*          The index of the first eigenvalue to be computed.
*
*  ILAST   (input) INTEGER
*          The index of the last eigenvalue to be computed.
*
*  RTOL1   (input) DOUBLE PRECISION
*          Tolerance for the convergence of the bisection intervals.
*          An interval [LEFT,RIGHT] has converged if
*          RIGHT-LEFT.LT.MAX( RTOL1*GAP, RTOL2*ABS((LEFT+RIGHT)/2) ),
*          where GAP is the (estimated) distance to the nearest
*          eigenvalue.
*
*  RTOL2   (input) DOUBLE PRECISION
*          Tolerance for the convergence of the bisection intervals.
*          An interval [LEFT,RIGHT] has converged if
*          RIGHT-LEFT.LT.MAX( RTOL1*GAP, RTOL2*ABS((LEFT+RIGHT)/2) ),
*          where GAP is the (estimated) distance to the nearest
*          eigenvalue.
*
*  OFFSET  (input) INTEGER
*          Offset for the arrays W, WGAP and WERR, i.e., the IFIRST-OFFSET
*          thru' ILAST-OFFSET elements of these arrays are to be used.
*
*  W       (input/output) DOUBLE PRECISION array, dimension (N)
*          On input, W( IFIRST-OFFSET ) thru' W( ILAST-OFFSET ) are
*          estimates of the eigenvalues of L D L^T indexed IFIRST thru' ILAST.
*          On output, these estimates are refined.
*
*  WGAP    (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On input, the (estimated) gaps between consecutive
*          eigenvalues of L D L^T, i.e., WGAP(I-OFFSET) is the gap between
*          eigenvalues I and I+1. Note that if IFIRST.EQ.ILAST
*          then WGAP(IFIRST-OFFSET) must be set to ZERO.
*          On output, these gaps are refined.
*
*  WERR    (input/output) DOUBLE PRECISION array, dimension (N)
*          On input, WERR( IFIRST-OFFSET ) thru' WERR( ILAST-OFFSET ) are
*          the errors in the estimates of the corresponding elements in W.
*          On output, these errors are refined.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)
*          Workspace.
*
*  IWORK   (workspace) INTEGER array, dimension (2*N)
*          Workspace.
*
*  INFO    (output) INTEGER
*          = 0: successful exit.
*          = 1: a bracket expansion was nonfinite or made no progress.
*          = 2: a bracket expansion exceeded MAXEXP iterations.
*          = 3: interval refinement exceeded MAXREF passes.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Inderjit Dhillon, University of Texas, Austin, USA
*     Osni Marques, LBNL/NERSC, USA
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            MAXEXP, MAXREF, NSTURM
      PARAMETER          ( MAXEXP = 256, MAXREF = 256, NSTURM = 4 )
      DOUBLE PRECISION   ZERO, ONE, TWO, HALF, TEN
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   HALF = 0.5D0, TEN = 10.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            B, BAFTER, CNT, I, II, I1, I2, ITER, J, K,
     $                   KK, NB, NEXT, NEXP, NLEFT, NINT, NRIGHT,
     $                   OLNINT, P, PREV
      DOUBLE PRECISION   DPLUS, EPS, ERROR, FAC, GAP, LEFT, MID, RIGHT,
     $                   OVFL, S, TMP, WIDTH
      INTEGER            BCNT( NSTURM ), BI( NSTURM ),
     $                   BNEXT( NSTURM ), BNRGHT( NSTURM )
      DOUBLE PRECISION   BMID( NSTURM )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*
*     Check to see if any of the initial eigenvalue
*     estimates is acceptable.
*
      INFO = 0
      EPS = DLAMCH( 'Precision' )
      OVFL = DLAMCH( 'Overflow' )
      DO 5 I = 1, 2*N
         IWORK( I ) = 0
   5  CONTINUE
      I1 = IFIRST
      I2 = IFIRST
      PREV = 0
      DO 10 I = IFIRST, ILAST
         II = I - OFFSET
         IF( I.EQ.IFIRST ) THEN
            GAP = WGAP( II )
         ELSE IF( I.EQ.ILAST ) THEN
            GAP = WGAP( II-1 )
         ELSE
            GAP = MIN( WGAP( II-1 ), WGAP( II ) )
         END IF
         ERROR = WERR( II )
         K = 2*I
*        IF( ERROR.LT.RTOL1*GAP ) THEN
*           WORK( K-1 ) = W( II ) - ERROR
*           WORK( K ) = W( II ) + ERROR
*           IWORK( K-1 ) = -1
*           IF( I1.EQ.I ) THEN
*              I1 = I1 + 1
*              PREV = I
*           END IF
*        ELSE
            IWORK( K-1 ) = 1
            I2 = I
*        END IF
   10 CONTINUE
*
*     Initialize unconverged intervals in [ WORK(2*I-1), WORK(2*I) ].
*     The Sturm Count, Count( WORK(2*I-1) ) is arranged to be I-1, while
*     Count( WORK(2*I) ) is stored in IWORK( 2*I ). The integer IWORK( 2*I-1 )
*     for an unconverged interval is set to the index of the next unconverged
*     interval, and is -1 or 0 for a converged interval. Thus a linked
*     list of unconverged intervals is set up.
*
      I = I1
      NINT = 0
   30 CONTINUE
      IF( I.LE.I2 ) THEN
         II = I - OFFSET
         IF( IWORK( 2*I-1 ).EQ.1 ) THEN
            FAC = ONE
            NEXP = 0
            LEFT = W( II ) - WERR( II )
*
*           Do while( CNT(LEFT).GT.I-1 )
*
   40       CONTINUE
            IF( I.GT.I1 .AND. LEFT.LE.RIGHT ) THEN
               LEFT = RIGHT
               CNT = I - 1
            ELSE
               S = -LEFT
               CNT = 0
               DO 50 J = 1, N - 1
                  DPLUS = D( J ) + S
                  S = S*LLD( J ) / DPLUS - LEFT
                  IF( DPLUS.LT.ZERO )
     $               CNT = CNT + 1
   50          CONTINUE
               DPLUS = D( N ) + S
               IF( DPLUS.LT.ZERO )
     $            CNT = CNT + 1
               IF( .NOT.( S.GT.ZERO .OR. S.LT.ONE ) ) THEN
*
*                 Runs a slower version of the above loop if a NaN is detected
*
                 CNT = 0
                 S = -LEFT
                 DO 55 J = 1, N - 1
                    DPLUS = D( J ) + S
                    IF( DPLUS.LT.ZERO )
     $                 CNT = CNT + 1
                    TMP = LLD( J ) / DPLUS
                    IF( TMP.EQ.ZERO ) THEN
                       S = LLD( J ) - LEFT
                    ELSE
                       S = S*TMP - LEFT
                    END IF
   55            CONTINUE
                 DPLUS = D( N ) + S
                 IF( DPLUS.LT.ZERO )
     $             CNT = CNT + 1
               END IF
               IF( CNT.GT.I-1 ) THEN
                  TMP = LEFT - WERR( II )*FAC
*                 A zero, NaN, or overflowing step cannot enlarge the
*                 bracket.  Return a positive INFO instead of repeating
*                 this Sturm count forever.
                  IF( .NOT.( TMP.LT.LEFT .AND.
     $                ABS( TMP ).LE.OVFL ) ) THEN
                     INFO = 1
                     RETURN
                  END IF
                  NEXP = NEXP + 1
                  IF( NEXP.GT.MAXEXP ) THEN
                     INFO = 2
                     RETURN
                  END IF
                  LEFT = TMP
                  FAC = TWO*FAC
                  GO TO 40
               END IF
            END IF
            NLEFT = CNT + 1
            I1 = MIN( I1, NLEFT )
            FAC = ONE
            NEXP = 0
            RIGHT = W( II ) + WERR( II )
*
*           Do while( CNT(RIGHT).LT.I )
*
   60       CONTINUE
            S = -RIGHT
            CNT = 0
            DO 70 J = 1, N - 1
               DPLUS = D( J ) + S
               S = S*LLD( J ) / DPLUS - RIGHT
               IF( DPLUS.LT.ZERO )
     $            CNT = CNT + 1
   70       CONTINUE
            DPLUS = D( N ) + S
            IF( DPLUS.LT.ZERO )
     $         CNT = CNT + 1
               IF( .NOT.( S.GT.ZERO .OR. S.LT.ONE ) ) THEN
*
*                 Runs a slower version of the above loop if a NaN is detected
*
                 CNT = 0
                 S = -RIGHT
                 DO 75 J = 1, N - 1
                    DPLUS = D( J ) + S
                    IF( DPLUS.LT.ZERO )
     $                 CNT = CNT + 1
                    TMP = LLD( J ) / DPLUS
                    IF( TMP.EQ.ZERO ) THEN
                       S = LLD( J ) - RIGHT
                    ELSE
                       S = S*TMP - RIGHT
                    END IF
   75            CONTINUE
                 DPLUS = D( N ) + S
                 IF( DPLUS.LT.ZERO )
     $             CNT = CNT + 1
               END IF
            IF( CNT.LT.I ) THEN
               TMP = RIGHT + WERR( II )*FAC
*              Apply the same progress/finite checks to the upper bracket.
               IF( .NOT.( TMP.GT.RIGHT .AND.
     $             ABS( TMP ).LE.OVFL ) ) THEN
                  INFO = 1
                  RETURN
               END IF
               NEXP = NEXP + 1
               IF( NEXP.GT.MAXEXP ) THEN
                  INFO = 2
                  RETURN
               END IF
               RIGHT = TMP
               FAC = TWO*FAC
               GO TO 60
            END IF
            CNT = MIN( CNT, I2 )
            NINT = NINT + 1
            K = 2*NLEFT
            WORK( K-1 ) = LEFT
            WORK( K ) = RIGHT
            I = CNT + 1
            IWORK( K-1 ) = I
            IWORK( K ) = CNT
            IF( PREV.NE.NLEFT-1 ) THEN
               WORK( K-2 ) = LEFT
            END IF
            PREV = NLEFT
         ELSE
            RIGHT = WORK( 2*I )
*
*           Remove converged interval from linked list
*
            IWORK( K-1 ) = IWORK( K-1 ) + 1
            PREV = I
            I = I + 1
         END IF
         GO TO 30
      END IF
      IF( I.LE.N .AND. IWORK( 2*I-1 ).NE.-1 )
     $   WORK( 2*I-1 ) = WORK( 2*PREV )
*
*     Do while( NINT.GT.0 ).  In exact arithmetic each pass halves every
*     active interval; cap the passes so NaN/non-progress cannot spin.
*
      ITER = 0
   80 CONTINUE
      ITER = ITER + 1
      IF( ITER.GT.MAXREF ) THEN
         INFO = 3
         RETURN
      END IF
      PREV = I1 - 1
      OLNINT = NINT
      I = I1
      P = 1
   85 CONTINUE
*
*     Snapshot up to NSTURM intervals from the pass's original linked
*     list, then interleave their independent Sturm recurrences.  The
*     interval updates below remain in their original order.  A split can
*     only insert an interval inside the current [I,NRIGHT] range, so it
*     cannot change the saved bounds or successor of a later original
*     interval in this pass.
*
      NB = MIN( NSTURM, OLNINT-P+1 )
      BAFTER = I
      DO 87 B = 1, NB
         BI( B ) = BAFTER
         K = 2*BAFTER
         BMID( B ) = HALF*( WORK( K-1 ) + WORK( K ) )
         BNEXT( B ) = IWORK( K-1 )
         BNRGHT( B ) = IWORK( K )
         BAFTER = BNEXT( B )
   87 CONTINUE
*     The four-lane count DLARRBTGK_STURM4 is kept below but not used;
*     each midpoint is counted by the reference one-shift routine.
      DO 88 B = 1, NB
         CALL DLARRBTGK_STURM_SLOW( N, D, LLD, BMID( B ), BCNT( B ) )
   88 CONTINUE
*
      DO 100 B = 1, NB
         I = BI( B )
         K = 2*I
         LEFT = WORK( K-1 )
         RIGHT = WORK( K )
         NEXT = BNEXT( B )
         NRIGHT = BNRGHT( B )
         MID = BMID( B )
         WIDTH = RIGHT - MID
         TMP = MAX( ABS( LEFT ), ABS( RIGHT ) )
*
*        Check for convergence of interval if there is only one
*        eigenvalue in the interval.
*
         GAP = ZERO
         IF( I.EQ.NRIGHT ) THEN
            IF( PREV.GT.0 .AND. NEXT.LE.N ) THEN
               GAP = MIN( LEFT-WORK( K-2 ), WORK( K+1 )-RIGHT )
            ELSE IF( PREV.GT.0 ) THEN
               GAP = LEFT - WORK( K-2 )
            ELSE IF( NEXT.LE.N ) THEN
               GAP = WORK( K+1 ) - RIGHT
            END IF
         END IF
         IF( WIDTH.LT.MAX( RTOL1*GAP, RTOL2*TMP ) ) THEN
            NINT = NINT - 1
            IWORK( K-1 ) = 0
            KK = K
            DO J = I+1, NRIGHT
               KK = KK+2
               IWORK( KK-1 ) = 0
               WORK( KK-1 ) = LEFT
               WORK( KK ) = RIGHT
               WGAP( J-1-OFFSET ) = ZERO
            END DO
            IF( I1.EQ.I ) THEN
               I1 = NEXT
            ELSE
               IWORK( 2*PREV-1 ) = NEXT
            END IF
            I = NEXT
            GO TO 100
         END IF
         PREV = I
*
*        Perform one bisection step.  BCNT(B) is exactly the scalar Sturm
*        count for MID; only independent recurrences were interleaved.
*
         CNT = BCNT( B )
         CNT = MAX( I-1, MIN( NRIGHT, CNT ) )
         IF( CNT.EQ.I-1 ) THEN
            WORK( K-1 ) = MID
         ELSE IF( CNT.EQ.NRIGHT ) THEN
            WORK( K ) = MID
         ELSE
            IWORK( K ) = CNT
            CNT = CNT + 1
            IWORK( K-1 ) = CNT
            KK = 2*CNT
            IWORK( KK-1 ) = NEXT
            IWORK( KK ) = NRIGHT
            WORK( K ) = MID
            WORK( KK-1 ) = MID
            WORK( KK ) = RIGHT
            PREV = CNT
            IF( CNT-1.GT.I ) THEN
               WORK( KK-2 ) = MID
            END IF
            IF( CNT.GT.IFIRST .AND. CNT.LE.ILAST ) THEN
               NINT = NINT + 1
            ELSE IF( CNT.LE.IFIRST ) THEN
               I1 = CNT
            END IF
         END IF
         I = NEXT
  100 CONTINUE
      I = BAFTER
      P = P + NB
      IF( P.LE.OLNINT ) GO TO 85
      IF( NINT.GT.0 )
     $   GO TO 80
      DO 110 I = IFIRST, ILAST
         K = 2*I
         II = I - OFFSET
         IF( IWORK( K-1 ).NE.-1 ) THEN
            W( II ) = HALF*( WORK( K-1 )+WORK( K ) )
            WERR( II ) = WORK( K ) - W( II )
            IF( I.NE.ILAST ) THEN
               WGAP( II ) = WORK( K+1 ) - WORK( K )
            END IF
         END IF
  110 CONTINUE
*
      RETURN
*
*     End of DLARRB
*
      END

      SUBROUTINE DLARRBTGK_STURM4( N, D, LLD, NB, X, CNT )
*
*     Evaluate one to four independent Sturm counts while preserving the
*     scalar operation order of each recurrence.  Interleaving exposes the
*     otherwise serial floating-point divisions to an out-of-order CPU.
*
      INTEGER            N, NB, CNT( * )
      DOUBLE PRECISION   D( * ), LLD( * ), X( * )
      INTEGER            C1, C2, C3, C4, J
      DOUBLE PRECISION   DP1, DP2, DP3, DP4, S1, S2, S3, S4
*
      C1 = 0
      S1 = -X( 1 )
      IF( NB.EQ.1 ) THEN
         DO 10 J = 1, N - 1
            DP1 = D( J ) + S1
            S1 = S1*LLD( J ) / DP1 - X( 1 )
            IF( DP1.LT.0.0D0 ) C1 = C1 + 1
   10    CONTINUE
         DP1 = D( N ) + S1
         IF( DP1.LT.0.0D0 ) C1 = C1 + 1
         CNT( 1 ) = C1
      ELSE IF( NB.EQ.2 ) THEN
         C2 = 0
         S2 = -X( 2 )
         DO 20 J = 1, N - 1
            DP1 = D( J ) + S1
            DP2 = D( J ) + S2
            S1 = S1*LLD( J ) / DP1 - X( 1 )
            S2 = S2*LLD( J ) / DP2 - X( 2 )
            IF( DP1.LT.0.0D0 ) C1 = C1 + 1
            IF( DP2.LT.0.0D0 ) C2 = C2 + 1
   20    CONTINUE
         DP1 = D( N ) + S1
         DP2 = D( N ) + S2
         IF( DP1.LT.0.0D0 ) C1 = C1 + 1
         IF( DP2.LT.0.0D0 ) C2 = C2 + 1
         CNT( 1 ) = C1
         CNT( 2 ) = C2
      ELSE IF( NB.EQ.3 ) THEN
         C2 = 0
         C3 = 0
         S2 = -X( 2 )
         S3 = -X( 3 )
         DO 30 J = 1, N - 1
            DP1 = D( J ) + S1
            DP2 = D( J ) + S2
            DP3 = D( J ) + S3
            S1 = S1*LLD( J ) / DP1 - X( 1 )
            S2 = S2*LLD( J ) / DP2 - X( 2 )
            S3 = S3*LLD( J ) / DP3 - X( 3 )
            IF( DP1.LT.0.0D0 ) C1 = C1 + 1
            IF( DP2.LT.0.0D0 ) C2 = C2 + 1
            IF( DP3.LT.0.0D0 ) C3 = C3 + 1
   30    CONTINUE
         DP1 = D( N ) + S1
         DP2 = D( N ) + S2
         DP3 = D( N ) + S3
         IF( DP1.LT.0.0D0 ) C1 = C1 + 1
         IF( DP2.LT.0.0D0 ) C2 = C2 + 1
         IF( DP3.LT.0.0D0 ) C3 = C3 + 1
         CNT( 1 ) = C1
         CNT( 2 ) = C2
         CNT( 3 ) = C3
      ELSE
         C2 = 0
         C3 = 0
         C4 = 0
         S2 = -X( 2 )
         S3 = -X( 3 )
         S4 = -X( 4 )
         DO 40 J = 1, N - 1
            DP1 = D( J ) + S1
            DP2 = D( J ) + S2
            DP3 = D( J ) + S3
            DP4 = D( J ) + S4
            S1 = S1*LLD( J ) / DP1 - X( 1 )
            S2 = S2*LLD( J ) / DP2 - X( 2 )
            S3 = S3*LLD( J ) / DP3 - X( 3 )
            S4 = S4*LLD( J ) / DP4 - X( 4 )
            IF( DP1.LT.0.0D0 ) C1 = C1 + 1
            IF( DP2.LT.0.0D0 ) C2 = C2 + 1
            IF( DP3.LT.0.0D0 ) C3 = C3 + 1
            IF( DP4.LT.0.0D0 ) C4 = C4 + 1
   40    CONTINUE
         DP1 = D( N ) + S1
         DP2 = D( N ) + S2
         DP3 = D( N ) + S3
         DP4 = D( N ) + S4
         IF( DP1.LT.0.0D0 ) C1 = C1 + 1
         IF( DP2.LT.0.0D0 ) C2 = C2 + 1
         IF( DP3.LT.0.0D0 ) C3 = C3 + 1
         IF( DP4.LT.0.0D0 ) C4 = C4 + 1
         CNT( 1 ) = C1
         CNT( 2 ) = C2
         CNT( 3 ) = C3
         CNT( 4 ) = C4
      END IF
*
*     Match the original guarded recurrence exactly if a fast recurrence
*     produced NaN.  This path is exceptional; the normal path above is
*     the performance target.
*
      IF( .NOT.( S1.GT.0.0D0 .OR. S1.LT.1.0D0 ) )
     $   CALL DLARRBTGK_STURM_SLOW( N, D, LLD, X( 1 ), CNT( 1 ) )
      IF( NB.GE.2 ) THEN
         IF( .NOT.( S2.GT.0.0D0 .OR. S2.LT.1.0D0 ) )
     $      CALL DLARRBTGK_STURM_SLOW( N, D, LLD, X( 2 ), CNT( 2 ) )
      END IF
      IF( NB.GE.3 ) THEN
         IF( .NOT.( S3.GT.0.0D0 .OR. S3.LT.1.0D0 ) )
     $      CALL DLARRBTGK_STURM_SLOW( N, D, LLD, X( 3 ), CNT( 3 ) )
      END IF
      IF( NB.GE.4 ) THEN
         IF( .NOT.( S4.GT.0.0D0 .OR. S4.LT.1.0D0 ) )
     $      CALL DLARRBTGK_STURM_SLOW( N, D, LLD, X( 4 ), CNT( 4 ) )
      END IF
      RETURN
      END

      SUBROUTINE DLARRBTGK_STURM_SLOW( N, D, LLD, X, CNT )
*
*     NaN-safe scalar recurrence, identical to DLARRB's original slow
*     path.  It is split out only so each lane of DLARRBTGK_STURM4 can use it.
*
      INTEGER            N, CNT
      DOUBLE PRECISION   D( * ), LLD( * ), X
      INTEGER            J
      DOUBLE PRECISION   DPLUS, S, TMP
*
      CNT = 0
      S = -X
      DO 10 J = 1, N - 1
         DPLUS = D( J ) + S
         IF( DPLUS.LT.0.0D0 ) CNT = CNT + 1
         TMP = LLD( J ) / DPLUS
         IF( TMP.EQ.0.0D0 ) THEN
            S = LLD( J ) - X
         ELSE
            S = S*TMP - X
         END IF
   10 CONTINUE
      DPLUS = D( N ) + S
      IF( DPLUS.LT.0.0D0 ) CNT = CNT + 1
      RETURN
      END
