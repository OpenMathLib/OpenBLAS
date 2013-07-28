*> \brief \b ZLAVHE
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLAVHE( UPLO, TRANS, DIAG, N, NRHS, A, LDA, IPIV, B,
*                          LDB, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          DIAG, TRANS, UPLO
*       INTEGER            INFO, LDA, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX*16         A( LDA, * ), B( LDB, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    ZLAVHE  performs one of the matrix-vector operations
*>       x := A*x  or  x := A^H*x,
*>    where x is an N element vector and  A is one of the factors
*>    from the symmetric factorization computed by ZHETRF.
*>    ZHETRF produces a factorization of the form
*>         U * D * U^H     or     L * D * L^H,
*>    where U (or L) is a product of permutation and unit upper (lower)
*>    triangular matrices, U^H (or L^H) is the conjugate transpose of
*>    U (or L), and D is Hermitian and block diagonal with 1 x 1 and
*>    2 x 2 diagonal blocks.  The multipliers for the transformations
*>    and the upper or lower triangular parts of the diagonal blocks
*>    are stored in the leading upper or lower triangle of the 2-D
*>    array A.
*>
*>    If TRANS = 'N' or 'n', ZLAVHE multiplies either by U or U * D
*>    (or L or L * D).
*>    If TRANS = 'C' or 'c', ZLAVHE multiplies either by U^H or D * U^H
*>    (or L^H or D * L^H ).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \verbatim
*>  UPLO   - CHARACTER*1
*>           On entry, UPLO specifies whether the triangular matrix
*>           stored in A is upper or lower triangular.
*>              UPLO = 'U' or 'u'   The matrix is upper triangular.
*>              UPLO = 'L' or 'l'   The matrix is lower triangular.
*>           Unchanged on exit.
*>
*>  TRANS  - CHARACTER*1
*>           On entry, TRANS specifies the operation to be performed as
*>           follows:
*>              TRANS = 'N' or 'n'   x := A*x.
*>              TRANS = 'C' or 'c'   x := A^H*x.
*>           Unchanged on exit.
*>
*>  DIAG   - CHARACTER*1
*>           On entry, DIAG specifies whether the diagonal blocks are
*>           assumed to be unit matrices:
*>              DIAG = 'U' or 'u'   Diagonal blocks are unit matrices.
*>              DIAG = 'N' or 'n'   Diagonal blocks are non-unit.
*>           Unchanged on exit.
*>
*>  N      - INTEGER
*>           On entry, N specifies the order of the matrix A.
*>           N must be at least zero.
*>           Unchanged on exit.
*>
*>  NRHS   - INTEGER
*>           On entry, NRHS specifies the number of right hand sides,
*>           i.e., the number of vectors x to be multiplied by A.
*>           NRHS must be at least zero.
*>           Unchanged on exit.
*>
*>  A      - COMPLEX*16 array, dimension( LDA, N )
*>           On entry, A contains a block diagonal matrix and the
*>           multipliers of the transformations used to obtain it,
*>           stored as a 2-D triangular matrix.
*>           Unchanged on exit.
*>
*>  LDA    - INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling ( sub ) program. LDA must be at least
*>           max( 1, N ).
*>           Unchanged on exit.
*>
*>  IPIV   - INTEGER array, dimension( N )
*>           On entry, IPIV contains the vector of pivot indices as
*>           determined by ZSYTRF or ZHETRF.
*>           If IPIV( K ) = K, no interchange was done.
*>           If IPIV( K ) <> K but IPIV( K ) > 0, then row K was inter-
*>           changed with row IPIV( K ) and a 1 x 1 pivot block was used.
*>           If IPIV( K ) < 0 and UPLO = 'U', then row K-1 was exchanged
*>           with row | IPIV( K ) | and a 2 x 2 pivot block was used.
*>           If IPIV( K ) < 0 and UPLO = 'L', then row K+1 was exchanged
*>           with row | IPIV( K ) | and a 2 x 2 pivot block was used.
*>
*>  B      - COMPLEX*16 array, dimension( LDB, NRHS )
*>           On entry, B contains NRHS vectors of length N.
*>           On exit, B is overwritten with the product A * B.
*>
*>  LDB    - INTEGER
*>           On entry, LDB contains the leading dimension of B as
*>           declared in the calling program.  LDB must be at least
*>           max( 1, N ).
*>           Unchanged on exit.
*>
*>  INFO   - INTEGER
*>           INFO is the error flag.
*>           On exit, a value of 0 indicates a successful exit.
*>           A negative value, say -K, indicates that the K-th argument
*>           has an illegal value.
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
*> \ingroup complex16_lin
*
*  =====================================================================
      SUBROUTINE ZLAVHE( UPLO, TRANS, DIAG, N, NRHS, A, LDA, IPIV, B,
     $                   LDB, INFO )
*
*  -- LAPACK test routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT
      INTEGER            J, K, KP
      COMPLEX*16         D11, D12, D21, D22, T1, T2
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMV, ZGERU, ZLACGV, ZSCAL, ZSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DCONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.LSAME( TRANS, 'C' ) )
     $          THEN
         INFO = -2
      ELSE IF( .NOT.LSAME( DIAG, 'U' ) .AND. .NOT.LSAME( DIAG, 'N' ) )
     $          THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLAVHE ', -INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
      NOUNIT = LSAME( DIAG, 'N' )
*------------------------------------------
*
*     Compute  B := A * B  (No transpose)
*
*------------------------------------------
      IF( LSAME( TRANS, 'N' ) ) THEN
*
*        Compute  B := U*B
*        where U = P(m)*inv(U(m))* ... *P(1)*inv(U(1))
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Loop forward applying the transformations.
*
            K = 1
   10       CONTINUE
            IF( K.GT.N )
     $         GO TO 30
            IF( IPIV( K ).GT.0 ) THEN
*
*              1 x 1 pivot block
*
*              Multiply by the diagonal element if forming U * D.
*
               IF( NOUNIT )
     $            CALL ZSCAL( NRHS, A( K, K ), B( K, 1 ), LDB )
*
*              Multiply by  P(K) * inv(U(K))  if K > 1.
*
               IF( K.GT.1 ) THEN
*
*                 Apply the transformation.
*
                  CALL ZGERU( K-1, NRHS, ONE, A( 1, K ), 1, B( K, 1 ),
     $                        LDB, B( 1, 1 ), LDB )
*
*                 Interchange if P(K) != I.
*
                  KP = IPIV( K )
                  IF( KP.NE.K )
     $               CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               END IF
               K = K + 1
            ELSE
*
*              2 x 2 pivot block
*
*              Multiply by the diagonal block if forming U * D.
*
               IF( NOUNIT ) THEN
                  D11 = A( K, K )
                  D22 = A( K+1, K+1 )
                  D12 = A( K, K+1 )
                  D21 = DCONJG( D12 )
                  DO 20 J = 1, NRHS
                     T1 = B( K, J )
                     T2 = B( K+1, J )
                     B( K, J ) = D11*T1 + D12*T2
                     B( K+1, J ) = D21*T1 + D22*T2
   20             CONTINUE
               END IF
*
*              Multiply by  P(K) * inv(U(K))  if K > 1.
*
               IF( K.GT.1 ) THEN
*
*                 Apply the transformations.
*
                  CALL ZGERU( K-1, NRHS, ONE, A( 1, K ), 1, B( K, 1 ),
     $                        LDB, B( 1, 1 ), LDB )
                  CALL ZGERU( K-1, NRHS, ONE, A( 1, K+1 ), 1,
     $                        B( K+1, 1 ), LDB, B( 1, 1 ), LDB )
*
*                 Interchange if P(K) != I.
*
                  KP = ABS( IPIV( K ) )
                  IF( KP.NE.K )
     $               CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               END IF
               K = K + 2
            END IF
            GO TO 10
   30       CONTINUE
*
*        Compute  B := L*B
*        where L = P(1)*inv(L(1))* ... *P(m)*inv(L(m)) .
*
         ELSE
*
*           Loop backward applying the transformations to B.
*
            K = N
   40       CONTINUE
            IF( K.LT.1 )
     $         GO TO 60
*
*           Test the pivot index.  If greater than zero, a 1 x 1
*           pivot was used, otherwise a 2 x 2 pivot was used.
*
            IF( IPIV( K ).GT.0 ) THEN
*
*              1 x 1 pivot block:
*
*              Multiply by the diagonal element if forming L * D.
*
               IF( NOUNIT )
     $            CALL ZSCAL( NRHS, A( K, K ), B( K, 1 ), LDB )
*
*              Multiply by  P(K) * inv(L(K))  if K < N.
*
               IF( K.NE.N ) THEN
                  KP = IPIV( K )
*
*                 Apply the transformation.
*
                  CALL ZGERU( N-K, NRHS, ONE, A( K+1, K ), 1, B( K, 1 ),
     $                        LDB, B( K+1, 1 ), LDB )
*
*                 Interchange if a permutation was applied at the
*                 K-th step of the factorization.
*
                  IF( KP.NE.K )
     $               CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               END IF
               K = K - 1
*
            ELSE
*
*              2 x 2 pivot block:
*
*              Multiply by the diagonal block if forming L * D.
*
               IF( NOUNIT ) THEN
                  D11 = A( K-1, K-1 )
                  D22 = A( K, K )
                  D21 = A( K, K-1 )
                  D12 = DCONJG( D21 )
                  DO 50 J = 1, NRHS
                     T1 = B( K-1, J )
                     T2 = B( K, J )
                     B( K-1, J ) = D11*T1 + D12*T2
                     B( K, J ) = D21*T1 + D22*T2
   50             CONTINUE
               END IF
*
*              Multiply by  P(K) * inv(L(K))  if K < N.
*
               IF( K.NE.N ) THEN
*
*                 Apply the transformation.
*
                  CALL ZGERU( N-K, NRHS, ONE, A( K+1, K ), 1, B( K, 1 ),
     $                        LDB, B( K+1, 1 ), LDB )
                  CALL ZGERU( N-K, NRHS, ONE, A( K+1, K-1 ), 1,
     $                        B( K-1, 1 ), LDB, B( K+1, 1 ), LDB )
*
*                 Interchange if a permutation was applied at the
*                 K-th step of the factorization.
*
                  KP = ABS( IPIV( K ) )
                  IF( KP.NE.K )
     $               CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               END IF
               K = K - 2
            END IF
            GO TO 40
   60       CONTINUE
         END IF
*--------------------------------------------------
*
*     Compute  B := A^H * B  (conjugate transpose)
*
*--------------------------------------------------
      ELSE
*
*        Form  B := U^H*B
*        where U  = P(m)*inv(U(m))* ... *P(1)*inv(U(1))
*        and   U^H = inv(U^H(1))*P(1)* ... *inv(U^H(m))*P(m)
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           Loop backward applying the transformations.
*
            K = N
   70       CONTINUE
            IF( K.LT.1 )
     $         GO TO 90
*
*           1 x 1 pivot block.
*
            IF( IPIV( K ).GT.0 ) THEN
               IF( K.GT.1 ) THEN
*
*                 Interchange if P(K) != I.
*
                  KP = IPIV( K )
                  IF( KP.NE.K )
     $               CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
*
*                 Apply the transformation
*                    y = y - B' conjg(x),
*                 where x is a column of A and y is a row of B.
*
                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
                  CALL ZGEMV( 'Conjugate', K-1, NRHS, ONE, B, LDB,
     $                        A( 1, K ), 1, ONE, B( K, 1 ), LDB )
                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
               END IF
               IF( NOUNIT )
     $            CALL ZSCAL( NRHS, A( K, K ), B( K, 1 ), LDB )
               K = K - 1
*
*           2 x 2 pivot block.
*
            ELSE
               IF( K.GT.2 ) THEN
*
*                 Interchange if P(K) != I.
*
                  KP = ABS( IPIV( K ) )
                  IF( KP.NE.K-1 )
     $               CALL ZSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ),
     $                           LDB )
*
*                 Apply the transformations
*                    y = y - B' conjg(x),
*                 where x is a block column of A and y is a block
*                 row of B.
*
                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
                  CALL ZGEMV( 'Conjugate', K-2, NRHS, ONE, B, LDB,
     $                        A( 1, K ), 1, ONE, B( K, 1 ), LDB )
                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
*
                  CALL ZLACGV( NRHS, B( K-1, 1 ), LDB )
                  CALL ZGEMV( 'Conjugate', K-2, NRHS, ONE, B, LDB,
     $                        A( 1, K-1 ), 1, ONE, B( K-1, 1 ), LDB )
                  CALL ZLACGV( NRHS, B( K-1, 1 ), LDB )
               END IF
*
*              Multiply by the diagonal block if non-unit.
*
               IF( NOUNIT ) THEN
                  D11 = A( K-1, K-1 )
                  D22 = A( K, K )
                  D12 = A( K-1, K )
                  D21 = DCONJG( D12 )
                  DO 80 J = 1, NRHS
                     T1 = B( K-1, J )
                     T2 = B( K, J )
                     B( K-1, J ) = D11*T1 + D12*T2
                     B( K, J ) = D21*T1 + D22*T2
   80             CONTINUE
               END IF
               K = K - 2
            END IF
            GO TO 70
   90       CONTINUE
*
*        Form  B := L^H*B
*        where L  = P(1)*inv(L(1))* ... *P(m)*inv(L(m))
*        and   L^H = inv(L^H(m))*P(m)* ... *inv(L^H(1))*P(1)
*
         ELSE
*
*           Loop forward applying the L-transformations.
*
            K = 1
  100       CONTINUE
            IF( K.GT.N )
     $         GO TO 120
*
*           1 x 1 pivot block
*
            IF( IPIV( K ).GT.0 ) THEN
               IF( K.LT.N ) THEN
*
*                 Interchange if P(K) != I.
*
                  KP = IPIV( K )
                  IF( KP.NE.K )
     $               CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
*
*                 Apply the transformation
*
                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
                  CALL ZGEMV( 'Conjugate', N-K, NRHS, ONE, B( K+1, 1 ),
     $                        LDB, A( K+1, K ), 1, ONE, B( K, 1 ), LDB )
                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
               END IF
               IF( NOUNIT )
     $            CALL ZSCAL( NRHS, A( K, K ), B( K, 1 ), LDB )
               K = K + 1
*
*           2 x 2 pivot block.
*
            ELSE
               IF( K.LT.N-1 ) THEN
*
*              Interchange if P(K) != I.
*
                  KP = ABS( IPIV( K ) )
                  IF( KP.NE.K+1 )
     $               CALL ZSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ),
     $                           LDB )
*
*                 Apply the transformation
*
                  CALL ZLACGV( NRHS, B( K+1, 1 ), LDB )
                  CALL ZGEMV( 'Conjugate', N-K-1, NRHS, ONE,
     $                        B( K+2, 1 ), LDB, A( K+2, K+1 ), 1, ONE,
     $                        B( K+1, 1 ), LDB )
                  CALL ZLACGV( NRHS, B( K+1, 1 ), LDB )
*
                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
                  CALL ZGEMV( 'Conjugate', N-K-1, NRHS, ONE,
     $                        B( K+2, 1 ), LDB, A( K+2, K ), 1, ONE,
     $                        B( K, 1 ), LDB )
                  CALL ZLACGV( NRHS, B( K, 1 ), LDB )
               END IF
*
*              Multiply by the diagonal block if non-unit.
*
               IF( NOUNIT ) THEN
                  D11 = A( K, K )
                  D22 = A( K+1, K+1 )
                  D21 = A( K+1, K )
                  D12 = DCONJG( D21 )
                  DO 110 J = 1, NRHS
                     T1 = B( K, J )
                     T2 = B( K+1, J )
                     B( K, J ) = D11*T1 + D12*T2
                     B( K+1, J ) = D21*T1 + D22*T2
  110             CONTINUE
               END IF
               K = K + 2
            END IF
            GO TO 100
  120       CONTINUE
         END IF
*
      END IF
      RETURN
*
*     End of ZLAVHE
*
      END
