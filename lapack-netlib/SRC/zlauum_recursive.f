*> \brief \b ZLAUUM_RECURSIVE computes the product UUH or LHL, where U and L are upper or lower triangular matrices (recursive algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLAUUM_RECURSIVE( UPLO, N, A, LDA, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16   A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLAUUM_RECURSIVE computes the product U * U**H or L**H * L, where the triangular
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
*> This is the recursive version of the algorithm
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
*>          overwritten with the upper triangle of the product U * U**H,
*>          and the strictly lower triangular part of A is not referenced.
*>          If UPLO = 'L', the lower triangle of A is overwritten with
*>          the lower triangle of the product L**H * L, and the
*>          strictly upper triangular part of A is not referenced. Regardless,
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
*> \ingroup lauum
*
*  =====================================================================
      RECURSIVE SUBROUTINE ZLAUUM_RECURSIVE( UPLO, N, A, LDA, INFO )
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
      COMPLEX*16   A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZONE
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0, ZONE = (1.0D+0, 0.0D+0) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            K, NX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLAUU2, ZHERK, ZTRMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, DCONJG, DBLE
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
         CALL XERBLA( 'ZLAUUM_RECURSIVE', -INFO )
         RETURN
      END IF
*
*     Early termination criteria
*
      IF( N.EQ.0 ) THEN
         RETURN
      END IF
*
*     Base Case
*
      IF( N.EQ.1 ) THEN
         A(1,1) = DBLE(A(1,1) * DCONJG(A(1,1)))
         RETURN
      END IF
*
*     Determine crossover point for when to bail to level2
*
      NX = ILAENV(3, "ZLAUUM_RECURSIVE", UPLO, N, -1, -1, -1)
      IF( N.LT.NX ) THEN
         CALL ZLAUU2(UPLO, N, A, LDA, INFO)
         RETURN
      END IF
*
*     Beginning of executable statements for the recursive case
*
      K = N/2
      IF( UPPER ) THEN
*
*        We are computing A = ut(U*U**H).
*
*        Break apart U as follows
*              |-----------------|
*        U =   | U_{11} U_{12}   |
*              | 0      U_{22}   |
*              |-----------------|
*
*        Where
*           U_{11}\in\C^{k\times k} U_{12}\in\C^{  k\times n-k}
*                                   U_{22}\in\C^{n-k\times n-k}
*
*        and U_{11},U_{22} are upper triangular and U_{12} is rectangular
*
*        This gives us our operations as
*                          |--------------------| |------------------------|
*        ut(U * U**H)   =  |  U_{11}   U_{12}   | |  U_{11}**H  0          |
*                          |  0        U_{22}   | |  U_{12}**H  U_{22}**H  |
*                          |--------------------| |------------------------|
*
*        Thus we get
*
*        U_{11} = U_{11}U_{11}**H + U_{12}U_{12}**H
*        U_{12} = U_{12}U_{22}**H
*
*        U_{22} = U_{22}U_{22}**H
*
*        We break these operations apart as follows
*
*        U_{11} = U_{11}U_{11}**H            (This subroutine)
*        U_{11} = U_{12}U_{12}**H + U_{11}   (HERK)
*
*        U_{12} = U_{12}U_{22}**H            (TRMM)
*
*        U_{22} = U_{22}U_{22}**H            (This subroutine)
*
*
*        Compute U_{11}
*
         CALL ZLAUUM_RECURSIVE(UPLO, K, A, LDA, INFO)
         CALL ZHERK('Upper', 'No Transpose', K, N-K,
     $      ONE, A(1,K+1), LDA, ONE, A, LDA)
*
*        Compute U_{12}
*
         CALL ZTRMM('Right', 'Upper', 'Conjugate', 'Non-unit',
     $      K, N-K, ZONE, A(K+1,K+1), LDA, A(1,K+1), LDA)
*
*        Compute U_{22}
*
         CALL ZLAUUM_RECURSIVE(UPLO, N-K, A(K+1,K+1), LDA, INFO)
      ELSE
*
*        We are computing A = lt(L**H*L).
*
*        Break apart L as follows
*              |-----------------|
*        L =   | L_{11} 0        |
*              | L_{21} L_{22}   |
*              |-----------------|
*
*        Where
*           L_{11}\in\C^{  k\times k}
*           L_{21}\in\C^{n-k\times k} l_{22}\in\C^{n-k\times n-k}
*
*        and L_{11},L_{22} are lower triangular and L_{21} is rectangular
*
*        This gives us our operations as
*                          |--------------------------| |-----------------|
*        lt(L**H * L)   =  |  L_{11}**H   L_{21}**H   | | L_{11} 0        |
*                          |  0           L_{22}**H   | | L_{21} L_{22}   |
*                          |--------------------------| |-----------------|
*
*        Thus we get
*
*        L_{11} = L_{11}**H L_{11} + L_{21}**H L_{21}
*        L_{21} = L_{22}**H L_{21}
*
*        L_{22} = L_{22}**H L_{22}
*
*        We break these operations apart as follows
*
*        L_{11} = L_{11}**H L_{11}           (This subroutine)
*        L_{11} = L_{21}**H L_{21} + L_{11}  (HERK)
*
*        L_{21} = L_{22}**H L_{21}           (TRMM)
*
*        L_{22} = L_{22}**H L_{22}           (This subroutine)
*
*
*        Compute L_{11}
*
         CALL ZLAUUM_RECURSIVE(UPLO, K, A, LDA, INFO)
         CALL ZHERK('Lower', 'Conjugate', K, N-K,
     $      ONE, A(K+1,1), LDA, ONE, A, LDA)
*
*        Compute L_{21}
*
         CALL ZTRMM('Left', 'Lower', 'Conjugate', 'Non-Unit',
     $      N-K, K, ZONE, A(K+1,K+1), LDA, A(K+1,1), LDA)
*
*        Compute L_{22}
*
         CALL ZLAUUM_RECURSIVE(UPLO, N-K, A(K+1,K+1), LDA, INFO)
      END IF
      END SUBROUTINE
