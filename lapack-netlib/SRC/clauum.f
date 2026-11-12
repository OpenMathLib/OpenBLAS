*> \brief \b CLAUUM computes the product UUH or LHL, where U and L are upper or lower triangular matrices (driver algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CLAUUM + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clauum.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clauum.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clauum.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CLAUUM( UPLO, N, A, LDA, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       COMPLEX            A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CLAUUM computes the product U * U**H or L**H * L, where the triangular
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
*> This is the driver that dispatches to either blocked or recursive
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
*>          overwritten with the upper triangle of the product U * U**H,
*>          and the strictly lower triangular part of A is not referenced.
*>          If UPLO = 'L', the lower triangle of A is overwritten with
*>          the lower triangle of the product L**H * L, and the strictly upper
*>          triangular part of A is not referenced. Regardless,
*>          the diagonal elements are explicitly set to be real by casting
*.          to REAL before exit.
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
*> \ingroup lauum
*
*  =====================================================================
      SUBROUTINE CLAUUM( UPLO, N, A, LDA, INFO )
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
      COMPLEX            A( LDA, * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, CLAUUM_RECURSIVE,
     $                   CLAUUM_BLOCKED
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
         CALL XERBLA( 'CLAUUM', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Here we dispatch to whatever is more efficient in a particular environment
*     We are defaulting to recursive, but if you want to use the blocked variant
*     Comment out the line starting with `CALL CLAUUM_RECURSIVE...`
*     and uncomment the line starting with `CALL CLAUUM_BLOCKED...`
*
      CALL CLAUUM_RECURSIVE(UPLO, N, A, LDA, INFO)
*      CALL CLAUUM_BLOCKED(UPLO, N, A, LDA, INFO)
      RETURN
*
*     End of CLAUUM
*
      END
