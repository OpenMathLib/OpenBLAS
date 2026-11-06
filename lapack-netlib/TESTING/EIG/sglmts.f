*> \brief \b SGLMTS
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGLMTS( N, M, P, A, AF, LDA, B, BF, LDB, D, DF,
*                          X, U, WORK, LWORK, RWORK, RESULT )
*
*       .. Scalar Arguments ..
*       INTEGER            LDA, LDB, LWORK, M, P, N
*       REAL               RESULT
*       ..
*       .. Array Arguments ..
*       REAL               A( LDA, * ), AF( LDA, * ), B( LDB, * ),
*      $                   BF( LDB, * ), RWORK( * ), D( * ), DF( * ),
*      $                   U( * ), WORK( LWORK ), X( * )
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGLMTS tests SGGGLM - a subroutine for solving the generalized
*> linear model problem, including independent scaling of its inputs.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of rows of the matrices A and B.  N >= 0.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of columns of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          P is INTEGER
*>          The number of columns of the matrix B.  P >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is REAL array, dimension (LDA,M)
*>          The N-by-M matrix A.
*> \endverbatim
*>
*> \param[out] AF
*> \verbatim
*>          AF is REAL array, dimension (LDA,M)
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the arrays A, AF. LDA >= max(M,N).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is REAL array, dimension (LDB,P)
*>          The N-by-P matrix A.
*> \endverbatim
*>
*> \param[out] BF
*> \verbatim
*>          BF is REAL array, dimension (LDB,P)
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the arrays B, BF. LDB >= max(P,N).
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is REAL array, dimension( N )
*>          On input, the left hand side of the GLM.
*> \endverbatim
*>
*> \param[out] DF
*> \verbatim
*>          DF is REAL array, dimension( N )
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is REAL array, dimension( M )
*>          solution vector X in the GLM problem.
*> \endverbatim
*>
*> \param[out] U
*> \verbatim
*>          U is REAL array, dimension( P )
*>          solution vector U in the GLM problem.
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
*>          RESULT is REAL
*>          The test ratio:
*>                           norm( d - A*x - B*u )
*>            RESULT = -----------------------------------------
*>                     (norm(A)+norm(B))*(norm(x)+norm(u))*EPS
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
      SUBROUTINE SGLMTS( N, M, P, A, AF, LDA, B, BF, LDB, D, DF,
     $                   X, U, WORK, LWORK, RWORK, RESULT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LWORK, M, P, N
      REAL               RESULT
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), AF( LDA, * ), B( LDB, * ),
     $                   BF( LDB, * ), RWORK( * ), D( * ), DF( * ),
     $                   U( * ), WORK( LWORK ), X( * )
*
*  ====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      INTEGER            MAXEXP
      PARAMETER          ( MAXEXP = MAXEXPONENT( ZERO ) - 2 )
*     ..
*     .. Local Scalars ..
      INTEGER            INFO, ISCALE, J
      REAL               ANORM, BNORM, EPS, XNORM, YNORM, DNORM, UNFL
      REAL               ASCL, BSCL, DSCL, SCL, TNRM
*     ..
*     .. External Functions ..
      REAL               SASUM, SLAMCH, SLANGE
      LOGICAL            SISNAN
      EXTERNAL           SISNAN, SASUM, SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SGEMV, SGGGLM, SLACPY, SSCAL
*
*     .. Intrinsic Functions ..
      INTRINSIC          EXPONENT, MAX, MAXEXPONENT, SCALE
*     ..
*     .. Executable Statements ..
*
      EPS = SLAMCH( 'Epsilon' )
      UNFL = SLAMCH( 'Safe minimum' )
      ANORM = MAX( SLANGE( '1', N, M, A, LDA, RWORK ), UNFL )
      BNORM = MAX( SLANGE( '1', N, P, B, LDB, RWORK ), UNFL )
*
*     Copy the matrices A and B to the arrays AF and BF,
*     and the vector D the array DF.
*
      CALL SLACPY( 'Full', N, M, A, LDA, AF, LDA )
      CALL SLACPY( 'Full', N, P, B, LDB, BF, LDB )
      CALL SCOPY( N, D, 1, DF, 1 )
*
*     Solve GLM problem
*
      CALL SGGGLM( N, M, P, AF, LDA, BF, LDB, DF, X, U, WORK, LWORK,
     $             INFO )
*
*     Test the residual for the solution of LSE
*
*                       norm( d - A*x - B*u )
*       RESULT = -----------------------------------------
*                (norm(A)+norm(B))*(norm(x)+norm(u))*EPS
*
      CALL SCOPY( N, D, 1, DF, 1 )
      CALL SGEMV( 'No transpose', N, M, -ONE, A, LDA, X, 1,
     $             ONE, DF, 1 )
*
      CALL SGEMV( 'No transpose', N, P, -ONE, B, LDB, U, 1,
     $             ONE, DF, 1 )
*
      DNORM = SASUM( N, DF, 1 )
      XNORM = SASUM( M, X, 1 ) + SASUM( P, U, 1 )
      YNORM = ANORM + BNORM
*
      IF( XNORM.LE.ZERO ) THEN
         RESULT = ZERO
      ELSE
         RESULT =  ( ( DNORM / YNORM ) / XNORM ) /EPS
      END IF
*
*     A, B and d may be scaled independently: for factors a, b, d,
*     the solutions become (d/a)*x and (d/b)*u.  Test large and
*     tiny inputs, then undo solution scaling before the residual.
*
      TNRM = MAX( SLANGE( 'M', N, M, A, LDA, RWORK ),
     $           SLANGE( 'M', N, P, B, LDB, RWORK ),
     $           SLANGE( 'M', N, 1, D, N, RWORK ) )
      IF( TNRM.GT.ZERO .AND. TNRM.LE.SLAMCH( 'Overflow' ) ) THEN
*
*        Cases: common large, common tiny, tiny d, tiny (A,d),
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
            DSCL = SCL
            IF( ISCALE.EQ.3 .OR. ISCALE.EQ.5 ) ASCL = ONE
            IF( ISCALE.EQ.3 .OR. ISCALE.EQ.4 ) BSCL = ONE
            CALL SLACPY( 'Full', N, M, A, LDA, AF, LDA )
            CALL SLACPY( 'Full', N, P, B, LDB, BF, LDB )
            CALL SCOPY( N, D, 1, DF, 1 )
            DO 10 J = 1, M
               CALL SSCAL( N, ASCL, AF( 1, J ), 1 )
   10       CONTINUE
            DO 20 J = 1, P
               CALL SSCAL( N, BSCL, BF( 1, J ), 1 )
   20       CONTINUE
            CALL SSCAL( N, DSCL, DF, 1 )
*
            CALL SGGGLM( N, M, P, AF, LDA, BF, LDB, DF, X, U, WORK,
     $                   LWORK, INFO )
            IF( INFO.NE.0 ) THEN
               RESULT = ONE / EPS
               GO TO 30
            END IF
            CALL SSCAL( M, ASCL / DSCL, X, 1 )
            CALL SSCAL( P, BSCL / DSCL, U, 1 )
*
            CALL SCOPY( N, D, 1, DF, 1 )
            CALL SGEMV( 'No transpose', N, M, -ONE, A, LDA, X, 1, ONE,
     $                  DF, 1 )
            CALL SGEMV( 'No transpose', N, P, -ONE, B, LDB, U, 1, ONE,
     $                  DF, 1 )
            DNORM = SASUM( N, DF, 1 )
            XNORM = SASUM( M, X, 1 ) + SASUM( P, U, 1 )
            IF( SISNAN( DNORM ) .OR. SISNAN( XNORM ) ) THEN
               RESULT = ONE / EPS
            ELSE IF( XNORM.GT.ZERO ) THEN
               RESULT = MAX( RESULT,
     $                       ( ( DNORM / YNORM ) / XNORM ) / EPS )
            END IF
   30    CONTINUE
      END IF
*
      RETURN
*
*     End of SGLMTS
*
      END
