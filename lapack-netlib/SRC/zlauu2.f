*> \brief \b ZLAUU2 computes the product UUH or LHL, where U and L are upper or lower triangular matrices (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download ZLAUU2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlauu2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlauu2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlauu2.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLAUU2( UPLO, N, A, LDA, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLAUU2 computes the product U * U**H or L**H * L, where the triangular
*> factor U or L is stored in the upper or lower triangular part of
*> the array A.
*>
*> If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
*> overwriting the factor U in A.
*> If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
*> overwriting the factor L in A.
*>
*> This is the unblocked form of the algorithm, calling Level 2 BLAS.
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
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the triangular factor U or L.
*>          On exit, if UPLO = 'U', the upper triangle of A is
*>          overwritten with the upper triangle of the product U * U**H;
*>          if UPLO = 'L', the lower triangle of A is overwritten with
*>          the lower triangle of the product L**H * L. Regardless,
*>          the diagonal elements are explicitly set to be real by casting
*.          to DBLE before exit.
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
*> \ingroup lauu2
*
*  =====================================================================
      SUBROUTINE ZLAUU2( UPLO, N, A, LDA, INFO )
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
      COMPLEX*16         A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, J
      COMPLEX*16         AII, AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZSCAL, ZGEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCONJG, MAX
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
         CALL XERBLA( 'ZLAUU2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Compute the product U * U**H.
*
         DO J = 1, N
            AJJ = DCONJG(A(J,J))
            IF( J.EQ.N ) THEN
*
*              We have this case separate just to ensure we never
*              access A out of bounds when we do A(J,J+1) = A(N,N+1)
*
               CALL ZSCAL( N, AJJ, A(1,N), 1 )
            ELSE
               CALL ZGEMM('No Transpose', 'Conjugate Transpose',
     $            J, 1, N-J, ONE, A(1,J+1), LDA, A(J,J+1), LDA,
     $            AJJ, A(1,J), LDA)
            END IF
*
*           Since in exact arithmetic the diagonal will be real on output
*           we explicitly cast the diagonal to real (zero out the
*           imaginary component)
*
            A(J,J) = DBLE(A(J,J))
         END DO
*
      ELSE
*
*        Compute the product L**H * L.
*
         DO I = 1, N
            AII = DCONJG(A(I,I))
            IF( I.EQ.N ) THEN
*
*              We have this case separate just to ensure we never
*              access A out of bounds when we do A(I+1,I) = A(N+1,N)
*
               CALL ZSCAL(N, AII, A(N,1), LDA)
            ELSE
               CALL ZGEMM('Conjugate Transpose', 'No Transpose',
     $            1, I, N-I, ONE, A(I+1,I), LDA, A(I+1,1), LDA,
     $            AII, A(I,1), LDA)
            END IF
*
*           Since in exact arithmetic the diagonal will be real on output
*           we explicitly cast the diagonal to real (zero out the
*           imaginary component)
*
            A(I,I) = DBLE(A(I,I))
         END DO
      END IF
*
      RETURN
*
*     End of ZLAUU2
*
      END
