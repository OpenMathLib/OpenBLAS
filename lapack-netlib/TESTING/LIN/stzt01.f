*> \brief \b STZT01
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       REAL             FUNCTION STZT01( M, N, A, AF, LDA, TAU, WORK,
*                        LWORK )
* 
*       .. Scalar Arguments ..
*       INTEGER            LDA, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       REAL               A( LDA, * ), AF( LDA, * ), TAU( * ),
*      $                   WORK( LWORK )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> STZT01 returns
*>      || A - R*Q || / ( M * eps * ||A|| )
*> for an upper trapezoidal A that was factored with STZRQF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrices A and AF.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrices A and AF.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          The original upper trapezoidal M by N matrix A.
*> \endverbatim
*>
*> \param[in] AF
*> \verbatim
*>          AF is REAL array, dimension (LDA,N)
*>          The output of STZRQF for input matrix A.
*>          The lower triangle is not referenced.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the arrays A and AF.
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is REAL array, dimension (M)
*>          Details of the  Householder transformations as returned by
*>          STZRQF.
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
*>          The length of the array WORK.  LWORK >= m*n + m.
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
*> \date November 2011
*
*> \ingroup single_lin
*
*  =====================================================================
      REAL             FUNCTION STZT01( M, N, A, AF, LDA, TAU, WORK,
     $                 LWORK )
*
*  -- LAPACK test routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), AF( LDA, * ), TAU( * ),
     $                   WORK( LWORK )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      REAL               NORMA
*     ..
*     .. Local Arrays ..
      REAL               RWORK( 1 )
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLANGE
      EXTERNAL           SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SAXPY, SLATZM, SLASET, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, REAL
*     ..
*     .. Executable Statements ..
*
      STZT01 = ZERO
*
      IF( LWORK.LT.M*N+M ) THEN
         CALL XERBLA( 'STZT01', 8 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      NORMA = SLANGE( 'One-norm', M, N, A, LDA, RWORK )
*
*     Copy upper triangle R
*
      CALL SLASET( 'Full', M, N, ZERO, ZERO, WORK, M )
      DO 20 J = 1, M
         DO 10 I = 1, J
            WORK( ( J-1 )*M+I ) = AF( I, J )
   10    CONTINUE
   20 CONTINUE
*
*     R = R * P(1) * ... *P(m)
*
      DO 30 I = 1, M
         CALL SLATZM( 'Right', I, N-M+1, AF( I, M+1 ), LDA, TAU( I ),
     $                WORK( ( I-1 )*M+1 ), WORK( M*M+1 ), M,
     $                WORK( M*N+1 ) )
   30 CONTINUE
*
*     R = R - A
*
      DO 40 I = 1, N
         CALL SAXPY( M, -ONE, A( 1, I ), 1, WORK( ( I-1 )*M+1 ), 1 )
   40 CONTINUE
*
      STZT01 = SLANGE( 'One-norm', M, N, WORK, M, RWORK )
*
      STZT01 = STZT01 / ( SLAMCH( 'Epsilon' )*REAL( MAX( M, N ) ) )
      IF( NORMA.NE.ZERO )
     $   STZT01 = STZT01 / NORMA
*
      RETURN
*
*     End of STZT01
*
      END
