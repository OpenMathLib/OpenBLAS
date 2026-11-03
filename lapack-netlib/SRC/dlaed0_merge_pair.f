*> \brief \b DLAED0_MERGE_PAIR executes one merge in the DC tree.
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
      SUBROUTINE DLAED0_MERGE_PAIR( MERGE, NMERGE, IWORK, WANTQ, D, E,
     $                               WORK, IBLO, IBHI, ISCR, IIWRK,
     $                               INDXQ, SUBMAT, MATSIZ, MSD2,
     $                               WRKBASE, IWBASE, KDEFL, INFO,
     $                               WRKSTR )
*
      INTEGER            IBHI, IBLO, IIWRK, INDXQ, INFO, ISCR, IWBASE,
     $                   KDEFL, MATSIZ, MERGE, MSD2, NMERGE, RIGHT,
     $                   SUBMAT, WRKBASE, WRKSTR
      LOGICAL            WANTQ
      INTEGER            IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), WORK( * )
      EXTERNAL           DLAED7_BR
*
      INTEGER            LEFT
*
      LEFT = 2*MERGE - 1
      RIGHT = LEFT + 1
      IF( MERGE.EQ.1 ) THEN
         SUBMAT = 1
         MATSIZ = IWORK( RIGHT )
         MSD2 = IWORK( 1 )
      ELSE
         SUBMAT = IWORK( LEFT-1 ) + 1
         MATSIZ = IWORK( RIGHT ) - IWORK( LEFT-1 )
         MSD2 = IWORK( LEFT ) - IWORK( LEFT-1 )
      END IF
*
      WRKBASE = ISCR + WRKSTR*( SUBMAT-1 )
      IWBASE = IIWRK + 5*( SUBMAT-1 )
*
      CALL DLAED7_BR( WANTQ, MATSIZ, MSD2, D( SUBMAT ),
     $                IWORK( INDXQ+SUBMAT-1 ), E( SUBMAT+MSD2-1 ),
     $                WORK( IBLO+SUBMAT-1 ), WORK( IBHI+SUBMAT-1 ),
     $                WORK( WRKBASE ), IWORK( IWBASE ), KDEFL, INFO )
      RETURN
      END
