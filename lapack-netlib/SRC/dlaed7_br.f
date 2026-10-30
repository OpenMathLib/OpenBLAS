*> \brief \b DLAED7_BR used by DSTEBR.
*
*  DLAED7_BR performs one divide-and-conquer merge while propagating
*  only the first and last rows of the current local eigenvector block.
*  It is called by DLAED0_BR and is not itself a public LAPACK API.
*
*  Serial Reference-LAPACK form: no OpenMP or thread-scaled workspace.
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
      SUBROUTINE DLAED7_BR( WANTQ, N, CUTPNT, D, INDXQ, RHO, BLO, BHI,
     $                       WORK, IWORK, KDEFL, INFO )
*
*     .. Scalar Arguments ..
      LOGICAL            WANTQ
      INTEGER            CUTPNT, INFO, KDEFL, N
      DOUBLE PRECISION   RHO
*     ..
*     .. Array Arguments ..
      INTEGER            INDXQ( * ), IWORK( * )
      DOUBLE PRECISION   BLO( * ), BHI( * ), D( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            GIVPTR, IDLMDA, IGIVCL, IGIVNM, I, IINDX,
     $                   IINDXP, INEAR, IPERM, IQ, IQ2, ISCR, ITAU, IW,
     $                   IWSGN, IZ, J, K, N1, N2, NBR
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAED8_BR, DLAED9_BR, DLAMRG, XERBLA
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
         CALL XERBLA( 'DLAED7_BR', -INFO )
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
      CALL DLAED8_BR( NBR, K, N, D, WORK( IQ ), 2, INDXQ, RHO,
     $                 CUTPNT, WORK( IZ ), WORK( IDLMDA ),
     $                 WORK( IQ2 ), 2, WORK( IW ), IWORK( IPERM ),
     $                 GIVPTR, IWORK( IGIVCL ), WORK( IGIVNM ),
     $                 IWORK( IINDXP ), IWORK( IINDX ), INFO )
      IF( INFO.NE.0 )
     $   GO TO 70
      KDEFL = K
*
*     Solve the secular equation and directly update the selected rows.
*     WORK(ISCR) holds one DLAED4 delta vector at a time; no K-by-K
*     delta or eigenvector block is materialized.
*
      IF( K.NE.0 ) THEN
         CALL DLAED9_BR( K, NBR, D, WORK( IQ2 ), 2, WORK( IQ ), 2,
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
         CALL DLAMRG( N1, N2, D, 1, -1, INDXQ )
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
*     End of DLAED7_BR
*
      END
