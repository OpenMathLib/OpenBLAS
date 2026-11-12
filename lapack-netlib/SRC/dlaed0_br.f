*> \brief \b DLAED0_BR used by DSTEBR.
*
*  DLAED0_BR computes eigenvalues of a symmetric tridiagonal matrix with
*  the boundary-row divide-and-conquer method used by DSTEBR.  Unlike
*  DLAED0(ICOMPQ=0), it keeps only first/last boundary rows for each
*  active block across merge levels instead of saving complete
*  eigenvector blocks.
*
*  Serial Reference-LAPACK form: no OpenMP, environment-variable
*  scheduling, or thread-scaled workspace.
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
      SUBROUTINE DLAED0_BR( N, D, E, WORK, IWORK, INFO )
*
*     .. Scalar Arguments ..
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            WRKSTR
      PARAMETER          ( WRKSTR = 14 )
*     ..
*     .. Local Scalars ..
      LOGICAL            WANTQ
      INTEGER            I, IBHI, IBLO, IIWRK, INDXQ, INFOM, ISCR,
     $                   IWBASE, J, K, KDEFL, MATSIZ, MERGE, MSD2,
     $                   NMERGE, QBASE, SMLSIZ, SMM1, SPM1, SUBMAT,
     $                   SUBPBS, WBASE, WRKBASE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLAED0_MERGE_PAIR,
     $                   DSTEQR_BOUNDARY_BR, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED0_BR', -INFO )
         RETURN
      END IF
*
      IF( N.EQ.0 )
     $   RETURN
*
      SMLSIZ = ILAENV( 9, 'DLAED0', ' ', 0, 0, 0, 0 )
*
*     Determine the size and placement of leaf submatrices.
*
      IWORK( 1 ) = N
      SUBPBS = 1
   10 CONTINUE
      IF( IWORK( SUBPBS ).GT.SMLSIZ ) THEN
         DO 20 J = SUBPBS, 1, -1
            IWORK( 2*J ) = ( IWORK( J )+1 ) / 2
            IWORK( 2*J-1 ) = IWORK( J ) / 2
   20    CONTINUE
         SUBPBS = 2*SUBPBS
         GO TO 10
      END IF
      DO 30 J = 2, SUBPBS
         IWORK( J ) = IWORK( J ) + IWORK( J-1 )
   30 CONTINUE
*
*     Apply rank-1 cuts.
*
      SPM1 = SUBPBS - 1
      DO 40 I = 1, SPM1
         SUBMAT = IWORK( I ) + 1
         SMM1 = SUBMAT - 1
         D( SMM1 ) = D( SMM1 ) - ABS( E( SMM1 ) )
         D( SUBMAT ) = D( SUBMAT ) - ABS( E( SMM1 ) )
   40 CONTINUE
*
*     Workspace layout:
*       WORK(IBLO:IBLO+N-1)       first boundary rows
*       WORK(IBHI:IBHI+N-1)       last boundary rows
*       WORK(ISCR:)               leaf scratch and per-merge scratch
*
*     Integer workspace layout:
*       IWORK(1:N)                split endpoints
*       IWORK(INDXQ:INDXQ+N-1)   sorting permutations
*       IWORK(IIWRK:)             DLAED7_BR scratch
*
      IBLO = 1
      IBHI = IBLO + N
      ISCR = IBHI + N
      INDXQ = N + 1
      IIWRK = INDXQ + N
*
*     Solve leaf eigenproblems and initialize boundary row invariants.
*
      DO 70 I = 1, SUBPBS
         QBASE = ISCR
         WBASE = QBASE + SMLSIZ*SMLSIZ
         IF( I.EQ.1 ) THEN
            SUBMAT = 1
            MATSIZ = IWORK( 1 )
         ELSE
            SUBMAT = IWORK( I-1 ) + 1
            MATSIZ = IWORK( I ) - IWORK( I-1 )
         END IF
*
         CALL DSTEQR_BOUNDARY_BR( MATSIZ, D( SUBMAT ), E( SUBMAT ),
     $                            WORK( IBLO+SUBMAT-1 ),
     $                            WORK( IBHI+SUBMAT-1 ),
     $                            WORK( QBASE ), SMLSIZ, WORK( WBASE ),
     $                            INFO )
         IF( INFO.NE.0 )
     $      GO TO 130
         K = 1
         DO 65 J = SUBMAT, IWORK( I )
            IWORK( INDXQ+J-1 ) = K
            K = K + 1
   65    CONTINUE
   70 CONTINUE
*
*     Merge adjacent eigensystems.  Boundary rows are updated in place
*     for the merged physical D order.
*
   80 CONTINUE
      IF( SUBPBS.GT.1 ) THEN
         WANTQ = SUBPBS.GT.2
         NMERGE = SUBPBS / 2
         DO 90 MERGE = 1, NMERGE
            CALL DLAED0_MERGE_PAIR( MERGE, NMERGE, IWORK, WANTQ, D, E,
     $           WORK, IBLO, IBHI, ISCR, IIWRK, INDXQ, SUBMAT, MATSIZ,
     $           MSD2, WRKBASE, IWBASE, KDEFL, INFOM, WRKSTR )
            IF( INFOM.NE.0 ) THEN
               INFO = INFOM
               GO TO 130
            END IF
   90    CONTINUE
*
         DO 95 MERGE = 1, NMERGE
            IWORK( MERGE ) = IWORK( 2*MERGE )
   95    CONTINUE
         SUBPBS = SUBPBS / 2
         GO TO 80
      END IF
*
*     Reorder eigenvalues from the final physical merge order.
*
      DO 100 I = 1, N
         J = IWORK( INDXQ+I-1 )
         WORK( ISCR+I-1 ) = D( J )
  100 CONTINUE
      CALL DCOPY( N, WORK( ISCR ), 1, D, 1 )
      GO TO 140
*
  130 CONTINUE
      INFO = SUBMAT*( N+1 ) + SUBMAT + MATSIZ - 1
*
  140 CONTINUE
      RETURN
*
*     End of DLAED0_BR
*
      END
