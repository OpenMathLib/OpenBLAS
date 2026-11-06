*> \brief \b SLSETS
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SLSETS( M, P, N, A, AF, LDA, B, BF, LDB, C, CF,
*                          D, DF, X, WORK, LWORK, RWORK, RESULT )
*
*       .. Scalar Arguments ..
*       INTEGER            LDA, LDB, LWORK, M, P, N
*       ..
*       .. Array Arguments ..
*       REAL               A( LDA, * ), AF( LDA, * ), B( LDB, * ),
*      $                   BF( LDB, * ), RESULT( 2 ), RWORK( * ),
*      $                   C( * ), D( * ), CF( * ), DF( * ),
*      $                   WORK( LWORK ), X( * )
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SLSETS tests SGGLSE - a subroutine for solving linear equality
*> constrained least square problem (LSE), including scaled inputs.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          P is INTEGER
*>          The number of rows of the matrix B.  P >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrices A and B.  N >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          The M-by-N matrix A.
*> \endverbatim
*>
*> \param[out] AF
*> \verbatim
*>          AF is REAL array, dimension (LDA,N)
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the arrays A, AF, Q and R.
*>          LDA >= max(M,N).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is REAL array, dimension (LDB,N)
*>          The P-by-N matrix A.
*> \endverbatim
*>
*> \param[out] BF
*> \verbatim
*>          BF is REAL array, dimension (LDB,N)
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the arrays B, BF, V and S.
*>          LDB >= max(P,N).
*> \endverbatim
*>
*> \param[in] C
*> \verbatim
*>          C is REAL array, dimension( M )
*>          the vector C in the LSE problem.
*> \endverbatim
*>
*> \param[out] CF
*> \verbatim
*>          CF is REAL array, dimension( M )
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is REAL array, dimension( P )
*>          the vector D in the LSE problem.
*> \endverbatim
*>
*> \param[out] DF
*> \verbatim
*>          DF is REAL array, dimension( P )
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is REAL array, dimension( N )
*>          solution vector X in the LSE problem.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (LWORK)
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL array, dimension (M)
*> \endverbatim
*>
*> \param[out] RESULT
*> \verbatim
*>          RESULT is REAL array, dimension (2)
*>          The test ratios:
*>            RESULT(1) = norm( A*x - c )/ norm(A)*norm(X)*EPS
*>            RESULT(2) = norm( B*x - d )/ norm(B)*norm(X)*EPS
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup single_eig
*
*  =====================================================================
      SUBROUTINE SLSETS( M, P, N, A, AF, LDA, B, BF, LDB, C, CF,
     $                   D, DF, X, WORK, LWORK, RWORK, RESULT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LWORK, M, P, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), AF( LDA, * ), B( LDB, * ),
     $                   BF( LDB, * ), RESULT( 2 ), RWORK( * ),
     $                   C( * ), D( * ), CF( * ), DF( * ),
     $                   WORK( LWORK ), X( * )
*
*  ====================================================================
*
*     ..
*     .. Local Scalars ..
      INTEGER            INFO, ISCALE, J
      REAL              ASCL, BSCL, CSCL, DSCL, RESID, SCL, TNRM
*     ..
*     .. Parameters ..
      REAL              ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      INTEGER            MAXEXP
      PARAMETER          ( MAXEXP = MAXEXPONENT( ZERO ) - 2 )
*     ..
*     .. External Functions ..
      LOGICAL            SISNAN
      REAL              SLANGE, SLAMCH
      EXTERNAL           SISNAN, SLANGE, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SGGLSE, SLACPY, SGET02, SSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          EXPONENT, MAX, MAXEXPONENT, SCALE
*     ..
*     .. Executable Statements ..
*
*     Copy the matrices A and B to the arrays AF and BF,
*     and the vectors C and D to the arrays CF and DF,
*
      CALL SLACPY( 'Full', M, N, A, LDA, AF, LDA )
      CALL SLACPY( 'Full', P, N, B, LDB, BF, LDB )
      CALL SCOPY( M, C, 1, CF, 1 )
      CALL SCOPY( P, D, 1, DF, 1 )
*
*     Solve LSE problem
*
      CALL SGGLSE( M, N, P, AF, LDA, BF, LDB, CF, DF, X,
     $             WORK, LWORK, INFO )
*
*     Test the residual for the solution of LSE
*
*     Compute RESULT(1) = norm( A*x - c ) / norm(A)*norm(X)*EPS
*
      CALL SCOPY( M, C, 1, CF, 1 )
      CALL SCOPY( P, D, 1, DF, 1 )
      CALL SGET02( 'No transpose', M, N, 1, A, LDA, X, N, CF, M,
     $             RWORK, RESULT( 1 ) )
*
*     Compute result(2) = norm( B*x - d ) / norm(B)*norm(X)*EPS
*
      CALL SGET02( 'No transpose', P, N, 1, B, LDB, X, N, DF, P,
     $             RWORK, RESULT( 2 ) )
*
*     Scaling (A,c) and (B,d) independently leaves x unchanged.
*     Scaling only (c,d) scales x by the same factor.  Exercise
*     both ends of the range and check residuals at the input scale.
*
      TNRM = MAX( SLANGE( 'M', M, N, A, LDA, RWORK ),
     $           SLANGE( 'M', P, N, B, LDB, RWORK ),
     $           SLANGE( 'M', M, 1, C, M, RWORK ),
     $           SLANGE( 'M', P, 1, D, P, RWORK ) )
      IF( TNRM.GT.ZERO .AND. TNRM.LE.SLAMCH( 'Overflow' ) ) THEN
*
*        Cases: common large, common tiny, tiny (c,d), tiny (A,c),
*        and tiny (B,d).
*
         DO 30 ISCALE = 1, 5
            IF( ISCALE.EQ.1 ) THEN
               SCL = SCALE( ONE, MAXEXP-EXPONENT( TNRM ) )
            ELSE
               SCL = SLAMCH( 'Safe minimum' ) /
     $               SLAMCH( 'Precision' )
               SCL = SCALE( ONE, EXPONENT( SCL )-4-
     $                      EXPONENT( TNRM ) )
            END IF
            ASCL = SCL
            BSCL = SCL
            CSCL = SCL
            DSCL = SCL
            IF( ISCALE.EQ.3 .OR. ISCALE.EQ.5 ) ASCL = ONE
            IF( ISCALE.EQ.3 .OR. ISCALE.EQ.4 ) BSCL = ONE
            IF( ISCALE.EQ.4 ) DSCL = ONE
            IF( ISCALE.EQ.5 ) CSCL = ONE
            CALL SLACPY( 'Full', M, N, A, LDA, AF, LDA )
            CALL SLACPY( 'Full', P, N, B, LDB, BF, LDB )
            CALL SCOPY( M, C, 1, CF, 1 )
            CALL SCOPY( P, D, 1, DF, 1 )
            DO 10 J = 1, N
               CALL SSCAL( M, ASCL, AF( 1, J ), 1 )
               CALL SSCAL( P, BSCL, BF( 1, J ), 1 )
   10       CONTINUE
            CALL SSCAL( M, CSCL, CF, 1 )
            CALL SSCAL( P, DSCL, DF, 1 )
*
            CALL SGGLSE( M, N, P, AF, LDA, BF, LDB, CF, DF, X, WORK,
     $                   LWORK, INFO )
            IF( INFO.NE.0 ) THEN
               RESULT( 1 ) = ONE / SLAMCH( 'Epsilon' )
               RESULT( 2 ) = RESULT( 1 )
               GO TO 30
            END IF
            IF( ISCALE.EQ.3 )
     $         CALL SSCAL( N, ONE / SCL, X, 1 )
*
            CALL SCOPY( M, C, 1, CF, 1 )
            CALL SCOPY( P, D, 1, DF, 1 )
            CALL SGET02( 'No transpose', M, N, 1, A, LDA, X, N, CF, M,
     $                   RWORK, RESID )
            IF( SISNAN( RESID ) ) THEN
               RESULT( 1 ) = ONE / SLAMCH( 'Epsilon' )
            ELSE
               RESULT( 1 ) = MAX( RESULT( 1 ), RESID )
            END IF
            CALL SGET02( 'No transpose', P, N, 1, B, LDB, X, N, DF, P,
     $                   RWORK, RESID )
            IF( SISNAN( RESID ) ) THEN
               RESULT( 2 ) = ONE / SLAMCH( 'Epsilon' )
            ELSE
               RESULT( 2 ) = MAX( RESULT( 2 ), RESID )
            END IF
   30    CONTINUE
      END IF
*
      RETURN
*
*     End of SLSETS
*
      END
