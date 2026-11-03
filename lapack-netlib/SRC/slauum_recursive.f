*> \brief \b SLAUUM_RECURSIVE computes the product UUH or LHL, where U and L are upper or lower triangular matrices (recursive algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SLAUUM_RECURSIVE( UPLO, N, A, LDA, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       REAL   A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SLAUUM_RECURSIVE computes the product U * U**T or L**T * L, where the triangular
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
*>          A is REAL array, dimension (LDA,N)
*>          On entry, the triangular factor U or L.
*>          On exit, if UPLO = 'U', the upper triangle of A is
*>          overwritten with the upper triangle of the product U * U**T,
*>          and the strictly lower triangular part of A is not referenced.
*>          If UPLO = 'L', the lower triangle of A is overwritten with
*>          the lower triangle of the product L**T * L, and the
*>          strictly upper triangular part of A is not referenced.
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
      RECURSIVE SUBROUTINE SLAUUM_RECURSIVE( UPLO, N, A, LDA, INFO )
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
      REAL   A( LDA, * )
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
      INTEGER            K, NX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, SLAUU2, SSYRK, STRMM
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
         CALL XERBLA( 'SLAUUM_RECURSIVE', -INFO )
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
         A(1,1) = A(1,1) * A(1,1)
         RETURN
      END IF
*
*     Determine crossover point for when to bail to level2
*
      NX = ILAENV(3, "SLAUUM_RECURSIVE", UPLO, N, -1, -1, -1)
      IF( N.LT.NX ) THEN
         CALL SLAUU2(UPLO, N, A, LDA, INFO)
         RETURN
      END IF
*
*     Beginning of executable statements for the recursive case
*
      K = N/2
      IF( UPPER ) THEN
*
*        We are computing A = ut(U*U**T).
*
*        Break apart U as follows
*              |-----------------|
*        U =   | U_{11} U_{12}   |
*              | 0      U_{22}   |
*              |-----------------|
*
*        Where
*           U_{11}\in\R^{k\times k} U_{12}\in\R^{  k\times n-k}
*                                   U_{22}\in\R^{n-k\times n-k}
*
*        and U_{11},U_{22} are upper triangular and U_{12} is rectangular
*
*        This gives us our operations as
*                          |--------------------| |------------------------|
*        ut(U * U**T)   =  |  U_{11}   U_{12}   | |  U_{11}**T  0          |
*                          |  0        U_{22}   | |  U_{12}**T  U_{22}**T  |
*                          |--------------------| |------------------------|
*
*        Thus we get
*
*        U_{11} = U_{11}U_{11}**T + U_{12}U_{12}**T
*        U_{12} = U_{12}U_{22}**T
*
*        U_{22} = U_{22}U_{22}**T
*
*        We break these operations apart as follows
*
*        U_{11} = U_{11}U_{11}**T            (This subroutine)
*        U_{11} = U_{12}U_{12}**T + U_{11}   (SYRK)
*
*        U_{12} = U_{12}U_{22}**T            (TRMM)
*
*        U_{22} = U_{22}U_{22}**T            (This subroutine)
*
*
*        Compute U_{11}
*
         CALL SLAUUM_RECURSIVE(UPLO, K, A, LDA, INFO)
         CALL SSYRK('Upper', 'No Transpose', K, N-K,
     $      ONE, A(1,K+1), LDA, ONE, A, LDA)
*
*        Compute U_{12}
*
         CALL STRMM('Right', 'Upper', 'Transpose', 'Non-unit',
     $      K, N-K, ONE, A(K+1,K+1), LDA, A(1,K+1), LDA)
*
*        Compute U_{22}
*
         CALL SLAUUM_RECURSIVE(UPLO, N-K, A(K+1,K+1), LDA, INFO)
      ELSE
*
*        We are computing A = lt(L**T*L).
*
*        Break apart L as follows
*              |-----------------|
*        L =   | L_{11} 0        |
*              | L_{21} L_{22}   |
*              |-----------------|
*
*        Where
*           L_{11}\in\R^{  k\times k}
*           L_{21}\in\R^{n-k\times k} l_{22}\in\R^{n-k\times n-k}
*
*        and L_{11},L_{22} are lower triangular and L_{21} is rectangular
*
*        This gives us our operations as
*                          |--------------------------| |-----------------|
*        lt(L**T * L)   =  |  L_{11}**T   L_{21}**T   | | L_{11} 0        |
*                          |  0           L_{22}**T   | | L_{21} L_{22}   |
*                          |--------------------------| |-----------------|
*
*        Thus we get
*
*        L_{11} = L_{11}**T L_{11} + L_{21}**T L_{21}
*        L_{21} = L_{22}**T L_{21}
*
*        L_{22} = L_{22}**T L_{22}
*
*        We break these operations apart as follows
*
*        L_{11} = L_{11}**T L_{11}           (This subroutine)
*        L_{11} = L_{21}**T L_{21} + L_{11}  (SYRK)
*
*        L_{21} = L_{22}**T L_{21}           (TRMM)
*
*        L_{22} = L_{22}**T L_{22}           (This subroutine)
*
*
*        Compute L_{11}
*
         CALL SLAUUM_RECURSIVE(UPLO, K, A, LDA, INFO)
         CALL SSYRK('Lower', 'Transpose', K, N-K,
     $      ONE, A(K+1,1), LDA, ONE, A, LDA)
*
*        Compute L_{21}
*
         CALL STRMM('Left', 'Lower', 'Transpose', 'Non-Unit',
     $      N-K, K, ONE, A(K+1,K+1), LDA, A(K+1,1), LDA)
*
*        Compute L_{22}
*
         CALL SLAUUM_RECURSIVE(UPLO, N-K, A(K+1,K+1), LDA, INFO)
      END IF
      END SUBROUTINE
