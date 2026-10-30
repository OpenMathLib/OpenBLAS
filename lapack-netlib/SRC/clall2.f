*> \brief \b CLALL2 computes the product LLH or UHU, where U and L are upper or lower triangular matrices (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CLALL2( UPLO, N, A, LDA, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       COMPLEX   A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CLALL2 computes the product L * L**H or U**H * U, where the triangular
*> factor U or L is stored in the upper or lower triangular part of
*> the array A.
*>
*> If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
*> overwriting the factor U in A, and the strictly lower triangular part
*> of A is not referenced.
*> If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
*> overwriting the factor L in A, and the strictly upper triangular part
*> of A is not referenced.
*>
*> This is the unblocked form of the algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the triangular factor stored in the array A
*>          is upper or lower triangular:
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the triangular factor U or L.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>          On entry, the triangular factor U or L.
*>          On exit, if UPLO = 'U', the upper triangle of A is
*>          overwritten with the upper triangle of the product U**h * U,
*>          and the strictly lower triangular part of A is not referenced.
*>          If UPLO = 'L', the lower triangle of A is overwritten with
*>          the lower triangle of the product L * L**H, and the strictly
*>          upper triangular part of A is not referenced. Regardless,
*>          the diagonal elements are explicitly set to be real by casting
*>          to REAL before exit.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -k, the k-th argument had an illegal value
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
*> \ingroup lall2
*
*  =====================================================================
      SUBROUTINE CLALL2( UPLO, N, A, LDA, INFO )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = (1.0E+0,0.0E+0) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I,J
      COMPLEX            AII,AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, CSCAL, CGEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, CONJG
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLALL2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
      IF( UPPER ) THEN
*
*        Compute and store ut(U**H U)
*
         DO J = N, 1, -1
*
*           Complex will want to store the conjugate
*
            AJJ = CONJG(A(J,J))
            IF( J.EQ.1 ) THEN
               CALL CSCAL(N, AJJ, A, LDA)
            ELSE
*
*              complex version will need {c,z}gemm since the vector x is
*              conjugated (or a {c,z}gemcv-like statement)
*
               CALL CGEMM('Conjugate Transpose', 'No Transpose',
     $            1, N-J+1, J-1, ONE, A(1,J), LDA, A(1,J), LDA,
     $            AJJ, A(J,J), LDA)
            END IF
            A(J,J) = REAL(A(J,J))
         END DO
      ELSE
*
*        Compute and store lt(LL**H)
*
         DO I = N, 1, -1
*
*           Complex will want to store the conjugate
*
            AII = CONJG(A(I,I))
            IF( I.EQ.1 ) THEN
               CALL CSCAL(N, AII, A, 1)
            ELSE
*
*              complex version will need {c,z}gemm since the vector x is
*              conjugated (or a {c,z}gemcv-like statement)
*
               CALL CGEMM('No Transpose', 'Conjugate Transpose',
     $            N-I+1, 1, I-1, ONE, A(I,1), LDA, A(I,1), LDA,
     $            AII, A(I,I), LDA)
            END IF
            A(I,I) = REAL(A(I,I))
         END DO
      END IF
      END SUBROUTINE
