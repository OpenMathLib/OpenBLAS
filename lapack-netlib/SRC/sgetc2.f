*> \brief \b SGETC2 computes the LU factorization with complete pivoting of the general n-by-n matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SGETC2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgetc2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgetc2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgetc2.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGETC2( N, A, LDA, IPIV, JPIV, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * ), JPIV( * )
*       REAL               A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGETC2 computes an LU factorization with complete pivoting of the
*> n-by-n matrix A. The factorization has the form A = P * L * U * Q,
*> where P and Q are permutation matrices, L is lower triangular with
*> unit diagonal elements and U is upper triangular.
*>
*> This is the Level 2 BLAS algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A. N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA, N)
*>          On entry, the n-by-n matrix A to be factored.
*>          On exit, the factors L and U from the factorization
*>          A = P*L*U*Q; the unit diagonal elements of L are not stored.
*>          If U(k, k) appears to be less than SMIN, U(k, k) is given the
*>          value of SMIN, i.e., giving a nonsingular perturbed system.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension(N).
*>          The pivot indices; for 1 <= i <= N, row i of the
*>          matrix has been interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[out] JPIV
*> \verbatim
*>          JPIV is INTEGER array, dimension(N).
*>          The pivot indices; for 1 <= j <= N, column j of the
*>          matrix has been interchanged with column JPIV(j).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>           = 0: successful exit
*>           > 0: if INFO = k, U(k, k) is likely to produce owerflow if
*>                we try to solve for x in Ax = b. So U is perturbed to
*>                avoid the overflow.
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
*> \date November 2013
*
*> \ingroup realGEauxiliary
*
*> \par Contributors:
*  ==================
*>
*>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
*>     Umea University, S-901 87 Umea, Sweden.
*
*  =====================================================================
      SUBROUTINE SGETC2( N, A, LDA, IPIV, JPIV, INFO )
*
*  -- LAPACK auxiliary routine (version 3.5.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2013
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), JPIV( * )
      REAL               A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IP, IPV, J, JP, JPV
      REAL               BIGNUM, EPS, SMIN, SMLNUM, XMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGER, SLABAD, SSWAP
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Set constants to control overflow
*
      INFO = 0
      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
      CALL SLABAD( SMLNUM, BIGNUM )
*
*     Factorize A using complete pivoting.
*     Set pivots less than SMIN to SMIN.
*
      DO 40 I = 1, N - 1
*
*        Find max element in matrix A
*
         XMAX = ZERO
         DO 20 IP = I, N
            DO 10 JP = I, N
               IF( ABS( A( IP, JP ) ).GE.XMAX ) THEN
                  XMAX = ABS( A( IP, JP ) )
                  IPV = IP
                  JPV = JP
               END IF
   10       CONTINUE
   20    CONTINUE
         IF( I.EQ.1 )
     $      SMIN = MAX( EPS*XMAX, SMLNUM )
*
*        Swap rows
*
         IF( IPV.NE.I )
     $      CALL SSWAP( N, A( IPV, 1 ), LDA, A( I, 1 ), LDA )
         IPIV( I ) = IPV
*
*        Swap columns
*
         IF( JPV.NE.I )
     $      CALL SSWAP( N, A( 1, JPV ), 1, A( 1, I ), 1 )
         JPIV( I ) = JPV
*
*        Check for singularity
*
         IF( ABS( A( I, I ) ).LT.SMIN ) THEN
            INFO = I
            A( I, I ) = SMIN
         END IF
         DO 30 J = I + 1, N
            A( J, I ) = A( J, I ) / A( I, I )
   30    CONTINUE
         CALL SGER( N-I, N-I, -ONE, A( I+1, I ), 1, A( I, I+1 ), LDA,
     $              A( I+1, I+1 ), LDA )
   40 CONTINUE
*
      IF( ABS( A( N, N ) ).LT.SMIN ) THEN
         INFO = N
         A( N, N ) = SMIN
      END IF
*
*     Set last pivots to N
*
      IPIV( N ) = N
      JPIV( N ) = N
*
      RETURN
*
*     End of SGETC2
*
      END
