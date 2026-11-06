*> \brief \b SLALLM_BLOCKED computes the product LLT or UTU, where U and L are upper or lower triangular matrices (blocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SLALLM_BLOCKED( UPLO, N, A, LDA, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       REAL               A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SLALLM_BLOCKED computes the product L * L**T or U**T * U, where the triangular
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
*> This is the blocked form of the algorithm, calling Level 3 BLAS.
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
*>          A is REAL array, dimension (LDA,N)
*>          On entry, the triangular factor U or L.
*>          On exit, if UPLO = 'U', the upper triangle of A is
*>          overwritten with the upper triangle of the product U**T * U,
*>          and the strictly lower triangular part of A is not referenced.
*>          If UPLO = 'L', the lower triangle of A is overwritten with
*>          the lower triangle of the product L * L**T, and the strictly
*>          upper triangular part of A is not referenced.
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
*> \ingroup lallm
*
*  =====================================================================
      SUBROUTINE SLALLM_BLOCKED( UPLO, N, A, LDA, INFO )
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
      REAL               A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL   ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IB, NB, K, KI
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, SLALL2, STRMM, SGEMM, SSYRK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
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
         CALL XERBLA( 'SLALLM_BLOCKED', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'SLALLM_BLOCKED', UPLO, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code
*
         CALL SLALL2( UPLO, N, A, LDA, INFO )
         RETURN
      END IF
*
*     Determine the size of the last block
*
      K = N / NB
*
*     The last block starts at KI
*
      KI = NB*K + 1
*
*     IB will hold the size of the block we are currently interested in
*
      IB = N - KI + 1
      IF( UPPER ) THEN
*
*        We go bottom to top computing one row at a time
*
         DO I = KI, 1, -NB
*
*           If not the first iteration
*
            IF( I+IB.LE.N ) THEN
*
*              Use the current diagonal block to update the right of it
*              We guard this to prevent potentially referencing invalid memory
*
               CALL STRMM('Left', 'Upper', 'Transpose', 'Non-Unit',
     $            IB, N-(I+IB)+1, ONE, A(I,I), LDA, A(I,I+IB), LDA)
*
*              Update the right of the current diagonal block with the
*              part of A above the current diagonal block
*              We don't do this on the last (first in the array)
*              However since we will not be accessing invalid
*              memory and sgemm promises to do nothing when
*              any input dimension is 0, we don't have to guard here
*
               CALL SGEMM('Transpose', 'No Transpose',
     $            IB, N-(I+IB)+1, I-1,
     $            ONE, A(1,I), LDA, A(1,I+IB), LDA,
     $            ONE, A(I,I+IB), LDA)
            END IF
*
*           Now, we update the diagonal block of interest
*
            CALL SLALL2('Upper', IB, A(I,I), LDA, INFO)
*
*           Update the diagonal block with the part of A above the diagonal
*
            CALL SSYRK('Upper', 'Transpose', IB, I-1,
     $         ONE, A(1,I), LDA, ONE, A(I,I), LDA)
*
*           Only the first iteration will have IB be different from NB
*
            IB = NB
         END DO
      ELSE
         DO I = KI, 1, -NB
*
*           If not the first iteration
*
            IF( I+IB.LE.N ) THEN
*
*              Use the current diagonal block to update below it
*              We guard this to prevent potentially referencing invalid memory
*
               CALL STRMM('Right', 'Lower', 'Transpose', 'Non-Unit',
     $            N-(I+IB)+1, IB, ONE, A(I,I), LDA, A(I+IB,I), LDA)
*
*              Update the right of the current diagonal block with the
*              part of A above the current diagonal block
*              We don't do this on the last (first in the array)
*              However since we will not be accessing invalid
*              memory and sgemm promises to do nothing when
*              any input dimension is 0, we don't have to guard here
*
               CALL SGEMM('No Transpose', 'Transpose',
     $            N-(I+IB)+1, IB, I-1,
     $            ONE, A(I+IB,1), LDA, A(I,1), LDA,
     $            ONE, A(I+IB,I), LDA)
            END IF
*
*           Now, we update the diagonal block of interest
*
            CALL SLALL2('Lower', IB, A(I,I), LDA, INFO)
*
*           Update the diagonal block with the part of A above the diagonal
*
            CALL SSYRK('Lower', 'No Transpose', IB, I-1,
     $         ONE, A(I,1), LDA, ONE, A(I,I), LDA)
*
*           Only the first iteration will have IB be different from NB
*
            IB = NB
         END DO
      END IF
      END SUBROUTINE
