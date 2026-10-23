*> \brief \b DSTEQR_BOUNDARY_BR solves one leaf and keeps boundary rows.
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
      SUBROUTINE DSTEQR_BOUNDARY_BR( N, D, E, BLO, BHI, Q, LDQ, WORK,
     $                               INFO )
*
      INTEGER            INFO, J, LDQ, N
      DOUBLE PRECISION   BLO( * ), BHI( * ), D( * ), E( * ), Q( * ),
     $                   WORK( * )
      EXTERNAL           DSTEQR
*
      CALL DSTEQR( 'I', N, D, E, Q, LDQ, WORK, INFO )
      IF( INFO.NE.0 )
     $   RETURN
*
      DO 10 J = 1, N
         BLO( J ) = Q( 1 + ( J - 1 )*LDQ )
         BHI( J ) = Q( N + ( J - 1 )*LDQ )
   10 CONTINUE
      RETURN
      END
