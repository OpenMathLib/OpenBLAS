*> \brief \b DLAED9_BR used by DSTEBR.
*
*  DLAED9_BR finds the roots of the secular equation for one
*  boundary-row DC merge and directly updates a small set of selected
*  rows.  It mirrors DLAED9's root solve and eigenvector normalization,
*  but avoids materializing the full K-by-K secular eigenvector block.
*  DLAED4_BR returns high-relative-accuracy DELTA values for each root;
*  consume those DELTA values immediately to form WPROD, following the
*  same numerical idea as LAPACK's SVD divide-and-conquer helpers.
*
*  Serial Reference-LAPACK form: no OpenMP or thread-private
*  accumulators.
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
      SUBROUTINE DLAED9_BR( K, NBR, D, QIN, LDQIN, QOUT, LDQOUT,
     $                       RHO, DLAMBDA, W, DELTA, WPROD, WSGN,
     $                       TAU, NEAR, INFO )
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDQIN, LDQOUT, NBR
      DOUBLE PRECISION   RHO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), DELTA( * ), DLAMBDA( * ),
     $                   NEAR( 2, * ), QIN( LDQIN, * ),
     $                   QOUT( LDQOUT, * ), TAU( * ), W( * ),
     $                   WPROD( * ), WSGN( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IORG, J, N1, N2, R
      DOUBLE PRECISION   ACC1, ACC2, DLT, ROWSUM, TEMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAED4, DLAED4_BR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      IF( K.LT.0 ) THEN
         INFO = -1
      ELSE IF( NBR.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDQIN.LT.MAX( 1, NBR ) ) THEN
         INFO = -5
      ELSE IF( LDQOUT.LT.MAX( 1, NBR ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED9_BR', -INFO )
         RETURN
      END IF
*
      IF( K.EQ.0 )
     $   RETURN
*
*     Root-only final merge: no parent boundary rows will be observed by
*     a higher merge level, so updated W and secular vectors are not
*     needed.
*
      IF( NBR.EQ.0 ) THEN
         DO 5 J = 1, K
            CALL DLAED4( K, J, DLAMBDA, W, DELTA, RHO, D( J ), INFO )
            IF( INFO.NE.0 )
     $         GO TO 170
    5    CONTINUE
         GO TO 170
      END IF
*
*     For K = 1 or 2, DLAED4/DLAED5 returns normalized eigenvectors
*     directly in DELTA.  Keep this path identical to DLAED9's
*     special-case copy to S.
*
      IF( K.EQ.1 .OR. K.EQ.2 ) THEN
         DO 30 J = 1, K
            CALL DLAED4( K, J, DLAMBDA, W, DELTA, RHO, D( J ),
     $                   INFO )
            IF( INFO.NE.0 )
     $         GO TO 170
            IF( NBR.EQ.2 ) THEN
               ACC1 = ZERO
               ACC2 = ZERO
               DO 10 I = 1, K
                  ACC1 = ACC1 + QIN( 1, I )*DELTA( I )
                  ACC2 = ACC2 + QIN( 2, I )*DELTA( I )
   10          CONTINUE
               QOUT( 1, J ) = ACC1
               QOUT( 2, J ) = ACC2
            ELSE
               DO 20 R = 1, NBR
                  ROWSUM = ZERO
                  DO 15 I = 1, K
                     ROWSUM = ROWSUM + QIN( R, I )*DELTA( I )
   15             CONTINUE
                  QOUT( R, J ) = ROWSUM
   20          CONTINUE
            END IF
   30    CONTINUE
         GO TO 170
      END IF
*
*     Preserve the update vector, solve all secular roots, then
*     accumulate the product needed to reconstruct DLAED9's updated
*     W.  The high-relative-accuracy DELTA values produced by
*     DLAED4_BR are consumed immediately rather than reconstructed
*     later from the compact origin/tau representation.
*
      DO 40 I = 1, K
         WSGN( I ) = W( I )
   40 CONTINUE
*
      DO 41 I = 1, K
         WPROD( I ) = ONE
   41 CONTINUE
*
      DO 44 J = 1, K
         CALL DLAED4_BR( K, J, DLAMBDA, WSGN, DELTA, RHO, D( J ),
     $                   IORG, TAU( J ), NEAR( 1, J ),
     $                   NEAR( 2, J ), INFO )
         IF( INFO.NE.0 )
     $      GO TO 170
*
         WPROD( J ) = WPROD( J )*DELTA( J )
         DO 42 I = 1, J - 1
            WPROD( I ) = WPROD( I )*
     $                   ( DELTA( I ) /
     $                     ( DLAMBDA( I )-DLAMBDA( J ) ) )
   42    CONTINUE
         DO 43 I = J + 1, K
            WPROD( I ) = WPROD( I )*
     $                   ( DELTA( I ) /
     $                     ( DLAMBDA( I )-DLAMBDA( J ) ) )
   43    CONTINUE
   44 CONTINUE
*
      DO 80 I = 1, K
         W( I ) = SIGN( SQRT( -WPROD( I ) ), WSGN( I ) )
   80 CONTINUE
*
*     Stream each secular eigenvector column through the selected rows.
*     Reconstruct the DELTA column from the compact origin/tau state,
*     using cached pole-adjacent values where cancellation is largest.
*
      IF( NBR.EQ.2 ) THEN
         DO 115 J = 1, K
            IF( J.EQ.K ) THEN
               IORG = K
               N1 = K - 1
               N2 = K
            ELSE IF( TAU( J ).GE.ZERO ) THEN
               IORG = J
               N1 = J
               N2 = J + 1
            ELSE
               IORG = J + 1
               N1 = J
               N2 = J + 1
            END IF
            TEMP = ZERO
            ACC1 = ZERO
            ACC2 = ZERO
*
            DO 85 I = 1, N1 - 1
               DLT = W( I ) /
     $               ( ( DLAMBDA( I )-DLAMBDA( IORG ) ) - TAU( J ) )
               TEMP = TEMP + DLT*DLT
               ACC1 = ACC1 + QIN( 1, I )*DLT
               ACC2 = ACC2 + QIN( 2, I )*DLT
   85       CONTINUE
*
            DLT = W( N1 ) / NEAR( 1, J )
            TEMP = TEMP + DLT*DLT
            ACC1 = ACC1 + QIN( 1, N1 )*DLT
            ACC2 = ACC2 + QIN( 2, N1 )*DLT
*
            DO 86 I = N1 + 1, N2 - 1
               DLT = W( I ) /
     $               ( ( DLAMBDA( I )-DLAMBDA( IORG ) ) - TAU( J ) )
               TEMP = TEMP + DLT*DLT
               ACC1 = ACC1 + QIN( 1, I )*DLT
               ACC2 = ACC2 + QIN( 2, I )*DLT
   86       CONTINUE
*
            DLT = W( N2 ) / NEAR( 2, J )
            TEMP = TEMP + DLT*DLT
            ACC1 = ACC1 + QIN( 1, N2 )*DLT
            ACC2 = ACC2 + QIN( 2, N2 )*DLT
*
            DO 87 I = N2 + 1, K
               DLT = W( I ) /
     $               ( ( DLAMBDA( I )-DLAMBDA( IORG ) ) - TAU( J ) )
               TEMP = TEMP + DLT*DLT
               ACC1 = ACC1 + QIN( 1, I )*DLT
               ACC2 = ACC2 + QIN( 2, I )*DLT
   87       CONTINUE
*
            TEMP = SQRT( TEMP )
            QOUT( 1, J ) = ACC1 / TEMP
            QOUT( 2, J ) = ACC2 / TEMP
  115    CONTINUE
         GO TO 170
      END IF
*
      DO 130 J = 1, K
         IF( J.EQ.K ) THEN
            IORG = K
         ELSE IF( TAU( J ).GE.ZERO ) THEN
            IORG = J
         ELSE
            IORG = J + 1
         END IF
         TEMP = ZERO
         DO 90 R = 1, NBR
            QOUT( R, J ) = ZERO
   90    CONTINUE
         DO 100 I = 1, K
            DLT = ( DLAMBDA( I )-DLAMBDA( IORG ) ) - TAU( J )
            IF( J.EQ.K ) THEN
               IF( I.EQ.K-1 )
     $            DLT = NEAR( 1, J )
               IF( I.EQ.K )
     $            DLT = NEAR( 2, J )
            ELSE
               IF( I.EQ.J )
     $            DLT = NEAR( 1, J )
               IF( I.EQ.J+1 )
     $            DLT = NEAR( 2, J )
            END IF
            DLT = W( I ) / DLT
            TEMP = TEMP + DLT*DLT
            DO 95 R = 1, NBR
               QOUT( R, J ) = QOUT( R, J ) + QIN( R, I )*DLT
   95       CONTINUE
  100    CONTINUE
         TEMP = SQRT( TEMP )
         DO 120 R = 1, NBR
            QOUT( R, J ) = QOUT( R, J ) / TEMP
  120    CONTINUE
  130 CONTINUE
*
  170 CONTINUE
      RETURN
*
*     End of DLAED9_BR
*
      END
