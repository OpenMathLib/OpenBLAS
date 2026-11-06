*> \brief \b DLAED8_BR selected-row variant of DLAED8.
*
*  DLAED8_BR mirrors DLAED8 deflation and bookkeeping but applies the
*  column rotations/permutations to NBR selected rows rather than to a
*  full eigenvector matrix.  It is an internal helper for DLAED7_BR.
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup stebr
*
*> \par Contributors:
*  ==================
*>
*> Ruiyi Zhan and Shaoshuai Zhang, University of Electronic Science
*> and Technology of China
*
      SUBROUTINE DLAED8_BR( NBR, K, N, D, Q, LDQ, INDXQ, RHO,
     $                      CUTPNT, Z, DLAMBDA, Q2, LDQ2, W, PERM,
     $                      GIVPTR, GIVCOL, GIVNUM, INDXP, INDX,
     $                      INFO )
*
*     .. Scalar Arguments ..
      INTEGER            CUTPNT, GIVPTR, INFO, K, LDQ, LDQ2, N, NBR
      DOUBLE PRECISION   RHO
*     ..
*     .. Array Arguments ..
      INTEGER            GIVCOL( 2, * ), INDX( * ), INDXP( * ),
     $                   INDXQ( * ), PERM( * )
      DOUBLE PRECISION   D( * ), DLAMBDA( * ), GIVNUM( 2, * ),
     $                   Q( LDQ, * ), Q2( LDQ2, * ), W( * ), Z( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   MONE, ZERO, ONE, TWO, EIGHT
      PARAMETER          ( MONE = -1.0D0, ZERO = 0.0D0,
     $                   ONE = 1.0D0, TWO = 2.0D0, EIGHT = 8.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IMAX, J, JLAM, JMAX, JP, K2, N1, N1P1, N2
      DOUBLE PRECISION   C, EPS, S, T, TAU, TOL
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH, DLAPY2
      EXTERNAL           IDAMAX, DLAMCH, DLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLACPY, DLAMRG, DROT, DSCAL, XERBLA
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
         CALL XERBLA( 'DLAED8_BR', -INFO )
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
         CALL DSCAL( N2, MONE, Z( N1P1 ), 1 )
      END IF
*
      T = ONE / SQRT( TWO )
      DO 80 J = 1, N
         INDX( J ) = J
   80 CONTINUE
      CALL DSCAL( N, T, Z, 1 )
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
      CALL DLAMRG( N1, N2, DLAMBDA, 1, 1, INDX )
      DO 110 I = 1, N
         D( I ) = DLAMBDA( INDX( I ) )
         Z( I ) = W( INDX( I ) )
  110 CONTINUE
*
      IMAX = IDAMAX( N, Z, 1 )
      JMAX = IDAMAX( N, D, 1 )
      EPS = DLAMCH( 'Epsilon' )
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
               CALL DCOPY( NBR, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 )
  120       CONTINUE
            CALL DLACPY( 'A', NBR, N, Q2, LDQ2, Q, LDQ )
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
         TAU = DLAPY2( C, S )
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
     $         CALL DROT( NBR, Q( 1, INDXQ( INDX( JLAM ) ) ), 1,
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
     $      CALL DCOPY( NBR, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 )
  180 CONTINUE
*
      IF( K.LT.N ) THEN
         CALL DCOPY( N-K, DLAMBDA( K+1 ), 1, D( K+1 ), 1 )
         IF( NBR.NE.0 )
     $      CALL DLACPY( 'A', NBR, N-K, Q2( 1, K+1 ), LDQ2,
     $                   Q( 1, K+1 ), LDQ )
      END IF
*
      RETURN
*
*     End of DLAED8_BR
*
      END
