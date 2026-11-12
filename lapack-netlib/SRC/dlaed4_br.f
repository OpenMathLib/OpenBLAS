*> \brief \b DLAED4_BR used by DSTEBR.
*
*  DLAED4_BR computes one secular root with DLAED4 and exposes a compact
*  representation of the returned high-relative-accuracy DELTA vector.
*  For N > 2, DLAED4's DELTA entries are all relative to either D(I)
*  or D(I+1), except for the final root where the origin is D(N).
*  ORG records that origin index and TAU satisfies
*
*     lambda = D(ORG) + TAU,
*     DELTA(j) = D(j) - D(ORG) - TAU.
*
*  NEAR1/NEAR2 cache the two pole-adjacent DELTA entries that are most
*  sensitive to cancellation when reconstructed from D and TAU.
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
      SUBROUTINE DLAED4_BR( N, I, D, Z, DELTA, RHO, DLAM, ORG,
     $                       TAU, NEAR1, NEAR2, INFO )
*
*     .. Scalar Arguments ..
      INTEGER            I, INFO, N, ORG
      DOUBLE PRECISION   DLAM, NEAR1, NEAR2, RHO, TAU
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), DELTA( * ), Z( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   HALF, ZERO
      PARAMETER          ( HALF = 0.5D0, ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   MID
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAED4, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      ORG = 1
      TAU = ZERO
      NEAR1 = ZERO
      NEAR2 = ZERO
*
      IF( N.LT.1 ) THEN
         INFO = -1
      ELSE IF( I.LT.1 .OR. I.GT.N ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED4_BR', -INFO )
         RETURN
      END IF
*
      CALL DLAED4( N, I, D, Z, DELTA, RHO, DLAM, INFO )
      IF( INFO.NE.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         ORG = 1
         TAU = DLAM - D( 1 )
         NEAR1 = DELTA( 1 )
         NEAR2 = ZERO
         RETURN
      END IF
*
      IF( N.EQ.2 ) THEN
         ORG = 1
         TAU = DLAM - D( 1 )
         NEAR1 = DELTA( 1 )
         NEAR2 = DELTA( 2 )
         RETURN
      END IF
*
      IF( I.EQ.N ) THEN
         ORG = N
      ELSE
         MID = HALF*( D( I )+D( I+1 ) )
         IF( DLAM.LT.MID ) THEN
            ORG = I
         ELSE
            ORG = I + 1
         END IF
      END IF
*
*     Avoid reconstructing TAU from DLAM-D(ORG).  DLAED4 already
*     returned DELTA(ORG) in its shifted coordinate, so using
*     TAU=-DELTA(ORG) preserves that final rounded denominator and
*     avoids an extra add/subtract rounding through DLAM.
*
      TAU = -DELTA( ORG )
*
      IF( I.EQ.N ) THEN
         NEAR1 = DELTA( N-1 )
         NEAR2 = DELTA( N )
      ELSE
         NEAR1 = DELTA( I )
         NEAR2 = DELTA( I+1 )
      END IF
*
      RETURN
*
*     End of DLAED4_BR
*
      END
