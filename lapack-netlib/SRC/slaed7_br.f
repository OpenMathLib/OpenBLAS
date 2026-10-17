*> \brief \b SLAED7_BR used by SSTEBR.
*
*  SLAED7_BR performs one divide-and-conquer merge while propagating
*  only the first and last rows of the current local eigenvector block.
*  It is called by SLAED0_BR and is not itself a public LAPACK API.
*
*  Serial Reference-LAPACK form: no OpenMP or thread-scaled workspace.
*
      SUBROUTINE SLAED7_BR( WANTQ, N, CUTPNT, D, INDXQ, RHO, BLO, BHI,
     $                       WORK, IWORK, KDEFL, INFO )
*
*     .. Scalar Arguments ..
      LOGICAL            WANTQ
      INTEGER            CUTPNT, INFO, KDEFL, N
      REAL               RHO
*     ..
*     .. Array Arguments ..
      INTEGER            INDXQ( * ), IWORK( * )
      REAL               BLO( * ), BHI( * ), D( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            GIVPTR, IDLMDA, IGIVCL, IGIVNM, I, IINDX,
     $                   IINDXP, INEAR, IPERM, IQ, IQ2, ISCR, ITAU, IW,
     $                   IWSGN, IZ, J, K, N1, N2, NBR
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLAED8_BR, SLAED9_BR, SLAMRG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      KDEFL = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( CUTPNT.LT.MIN( 1, N ) .OR. CUTPNT.GT.N ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLAED7_BR', -INFO )
         RETURN
      END IF
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Workspace layout (14*N doubles):
*
*       Z, DLAMBDA, W                  : N each
*       QBR, Q2BR, GIVNUM/reused WPROD : 2*N each
*       DELTA, WSGN, TAU, NEAR          : N, N, N, 2*N
*
      IZ = 1
      IDLMDA = IZ + N
      IW = IDLMDA + N
      IQ = IW + N
      IQ2 = IQ + 2*N
      IGIVNM = IQ2 + 2*N
      ISCR = IGIVNM + 2*N
      IWSGN = ISCR + N
      ITAU = IWSGN + N
      INEAR = ITAU + N
*
*     Integer workspace layout.
*
      IINDXP = 1
      IINDX = IINDXP + N
      IPERM = IINDX + N
      IGIVCL = IPERM + N
*
      IF( WANTQ ) THEN
         NBR = 2
      ELSE
         NBR = 0
      END IF
*
*     Z is the rank-one update vector.  QBR has two selected rows of
*     blockdiag(Q_left,Q_right): parent first boundary and parent last
*     boundary, represented in the child eigenbases.  The root merge has
*     no parent, so it only needs Z for the secular eigenvalues.
*
      IF( WANTQ ) THEN
         DO 20 J = 1, CUTPNT
            WORK( IZ+J-1 ) = BHI( J )
            WORK( IQ+2*( J-1 ) ) = BLO( J )
            WORK( IQ+2*( J-1 )+1 ) = ZERO
   20    CONTINUE
         DO 30 J = CUTPNT + 1, N
            WORK( IZ+J-1 ) = BLO( J )
            WORK( IQ+2*( J-1 ) ) = ZERO
            WORK( IQ+2*( J-1 )+1 ) = BHI( J )
   30    CONTINUE
      ELSE
         DO 40 J = 1, CUTPNT
            WORK( IZ+J-1 ) = BHI( J )
   40    CONTINUE
         DO 50 J = CUTPNT + 1, N
            WORK( IZ+J-1 ) = BLO( J )
   50    CONTINUE
      END IF
*
*     Sort and deflate, applying the same Givens/permutation actions to
*     the selected boundary rows when a parent merge will need them.
*
*
      CALL SLAED8_BR( NBR, K, N, D, WORK( IQ ), 2, INDXQ, RHO,
     $                 CUTPNT, WORK( IZ ), WORK( IDLMDA ),
     $                 WORK( IQ2 ), 2, WORK( IW ), IWORK( IPERM ),
     $                 GIVPTR, IWORK( IGIVCL ), WORK( IGIVNM ),
     $                 IWORK( IINDXP ), IWORK( IINDX ), INFO )
      IF( INFO.NE.0 )
     $   GO TO 70
      KDEFL = K
*
*     Solve the secular equation and directly update the selected rows.
*     WORK(ISCR) holds one SLAED4 delta vector at a time; no K-by-K
*     delta or eigenvector block is materialized.
*
      IF( K.NE.0 ) THEN
         CALL SLAED9_BR( K, NBR, D, WORK( IQ2 ), 2, WORK( IQ ), 2,
     $                   RHO, WORK( IDLMDA ), WORK( IW ),
     $                   WORK( ISCR ), WORK( IGIVNM ),
     $                   WORK( IWSGN ), WORK( ITAU ),
     $                   WORK( INEAR ), INFO )
         IF( INFO.NE.0 )
     $      GO TO 70
*
*        Prepare INDXQ sorting permutation.
*
         N1 = K
         N2 = N - K
         CALL SLAMRG( N1, N2, D, 1, -1, INDXQ )
      ELSE
         DO 60 I = 1, N
            INDXQ( I ) = I
   60    CONTINUE
      END IF
*
*     Persist the parent boundary rows in the current physical D order.
*
      IF( WANTQ ) THEN
         DO 65 J = 1, N
            BLO( J ) = WORK( IQ+2*( J-1 ) )
            BHI( J ) = WORK( IQ+2*( J-1 )+1 )
   65    CONTINUE
      END IF
*
   70 CONTINUE
      RETURN
*
*     End of SLAED7_BR
*
      END
*
*> \brief \b SLAED8_BR selected-row variant of SLAED8.
*
*  SLAED8_BR mirrors SLAED8 deflation and bookkeeping but applies the
*  column rotations/permutations to NBR selected rows rather than to a
*  full eigenvector matrix.  It is an internal helper for SLAED7_BR.
*
      SUBROUTINE SLAED8_BR( NBR, K, N, D, Q, LDQ, INDXQ, RHO,
     $                      CUTPNT, Z, DLAMBDA, Q2, LDQ2, W, PERM,
     $                      GIVPTR, GIVCOL, GIVNUM, INDXP, INDX,
     $                      INFO )
*
*     .. Scalar Arguments ..
      INTEGER            CUTPNT, GIVPTR, INFO, K, LDQ, LDQ2, N, NBR
      REAL               RHO
*     ..
*     .. Array Arguments ..
      INTEGER            GIVCOL( 2, * ), INDX( * ), INDXP( * ),
     $                   INDXQ( * ), PERM( * )
      REAL               D( * ), DLAMBDA( * ), GIVNUM( 2, * ),
     $                   Q( LDQ, * ), Q2( LDQ2, * ), W( * ), Z( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               MONE, ZERO, ONE, TWO, EIGHT
      PARAMETER          ( MONE = -1.0E0, ZERO = 0.0E0,
     $                   ONE = 1.0E0, TWO = 2.0E0, EIGHT = 8.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IMAX, J, JLAM, JMAX, JP, K2, N1, N1P1, N2
      REAL               C, EPS, S, T, TAU, TOL
*     ..
*     .. External Functions ..
      INTEGER            ISAMAX
      REAL               SLAMCH, SLAPY2
      EXTERNAL           ISAMAX, SLAMCH, SLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SLACPY, SLAMRG, SROT, SSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      IF( NBR.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDQ.LT.MAX( 1, NBR ) ) THEN
         INFO = -6
      ELSE IF( CUTPNT.LT.MIN( 1, N ) .OR. CUTPNT.GT.N ) THEN
         INFO = -9
      ELSE IF( LDQ2.LT.MAX( 1, NBR ) ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLAED8_BR', -INFO )
         RETURN
      END IF
*
      GIVPTR = 0
      IF( N.EQ.0 )
     $   RETURN
*
      N1 = CUTPNT
      N2 = N - N1
      N1P1 = N1 + 1
*
      IF( RHO.LT.ZERO ) THEN
         CALL SSCAL( N2, MONE, Z( N1P1 ), 1 )
      END IF
*
      T = ONE / SQRT( TWO )
      DO 80 J = 1, N
         INDX( J ) = J
   80 CONTINUE
      CALL SSCAL( N, T, Z, 1 )
      RHO = ABS( TWO*RHO )
*
*     Sort the eigenvalues into increasing order.
*
      DO 90 I = CUTPNT + 1, N
         INDXQ( I ) = INDXQ( I ) + CUTPNT
   90 CONTINUE
      DO 100 I = 1, N
         DLAMBDA( I ) = D( INDXQ( I ) )
         W( I ) = Z( INDXQ( I ) )
  100 CONTINUE
      CALL SLAMRG( N1, N2, DLAMBDA, 1, 1, INDX )
      DO 110 I = 1, N
         D( I ) = DLAMBDA( INDX( I ) )
         Z( I ) = W( INDX( I ) )
  110 CONTINUE
*
      IMAX = ISAMAX( N, Z, 1 )
      JMAX = ISAMAX( N, D, 1 )
      EPS = SLAMCH( 'Epsilon' )
      TOL = EIGHT*EPS*ABS( D( JMAX ) )
*
      IF( RHO*ABS( Z( IMAX ) ).LE.TOL ) THEN
         K = 0
         IF( NBR.EQ.0 ) THEN
            DO 115 J = 1, N
               PERM( J ) = INDXQ( INDX( J ) )
  115       CONTINUE
         ELSE
            DO 120 J = 1, N
               PERM( J ) = INDXQ( INDX( J ) )
               CALL SCOPY( NBR, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 )
  120       CONTINUE
            CALL SLACPY( 'A', NBR, N, Q2, LDQ2, Q, LDQ )
         END IF
         RETURN
      END IF
*
      K = 0
      K2 = N + 1
      DO 130 J = 1, N
         IF( RHO*ABS( Z( J ) ).LE.TOL ) THEN
            K2 = K2 - 1
            INDXP( K2 ) = J
            IF( J.EQ.N )
     $         GO TO 170
         ELSE
            JLAM = J
            GO TO 140
         END IF
  130 CONTINUE
  140 CONTINUE
      J = J + 1
      IF( J.GT.N )
     $   GO TO 160
      IF( RHO*ABS( Z( J ) ).LE.TOL ) THEN
         K2 = K2 - 1
         INDXP( K2 ) = J
      ELSE
         S = Z( JLAM )
         C = Z( J )
         TAU = SLAPY2( C, S )
         T = D( J ) - D( JLAM )
         C = C / TAU
         S = -S / TAU
         IF( ABS( T*C*S ).LE.TOL ) THEN
            Z( J ) = TAU
            Z( JLAM ) = ZERO
*
            GIVPTR = GIVPTR + 1
            GIVCOL( 1, GIVPTR ) = INDXQ( INDX( JLAM ) )
            GIVCOL( 2, GIVPTR ) = INDXQ( INDX( J ) )
            GIVNUM( 1, GIVPTR ) = C
            GIVNUM( 2, GIVPTR ) = S
            IF( NBR.NE.0 )
     $         CALL SROT( NBR, Q( 1, INDXQ( INDX( JLAM ) ) ), 1,
     $                    Q( 1, INDXQ( INDX( J ) ) ), 1, C, S )
            T = D( JLAM )*C*C + D( J )*S*S
            D( J ) = D( JLAM )*S*S + D( J )*C*C
            D( JLAM ) = T
            K2 = K2 - 1
            I = 1
  150       CONTINUE
            IF( K2+I.LE.N ) THEN
               IF( D( JLAM ).LT.D( INDXP( K2+I ) ) ) THEN
                  INDXP( K2+I-1 ) = INDXP( K2+I )
                  INDXP( K2+I ) = JLAM
                  I = I + 1
                  GO TO 150
               ELSE
                  INDXP( K2+I-1 ) = JLAM
               END IF
            ELSE
               INDXP( K2+I-1 ) = JLAM
            END IF
            JLAM = J
         ELSE
            K = K + 1
            W( K ) = Z( JLAM )
            DLAMBDA( K ) = D( JLAM )
            INDXP( K ) = JLAM
            JLAM = J
         END IF
      END IF
      GO TO 140
  160 CONTINUE
*
      K = K + 1
      W( K ) = Z( JLAM )
      DLAMBDA( K ) = D( JLAM )
      INDXP( K ) = JLAM
*
  170 CONTINUE
      DO 180 J = 1, N
         JP = INDXP( J )
         DLAMBDA( J ) = D( JP )
         PERM( J ) = INDXQ( INDX( JP ) )
         IF( NBR.NE.0 )
     $      CALL SCOPY( NBR, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 )
  180 CONTINUE
*
      IF( K.LT.N ) THEN
         CALL SCOPY( N-K, DLAMBDA( K+1 ), 1, D( K+1 ), 1 )
         IF( NBR.NE.0 )
     $      CALL SLACPY( 'A', NBR, N-K, Q2( 1, K+1 ), LDQ2,
     $                   Q( 1, K+1 ), LDQ )
      END IF
*
      RETURN
*
*     End of SLAED8_BR
*
      END
