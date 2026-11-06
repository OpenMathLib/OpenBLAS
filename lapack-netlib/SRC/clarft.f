*> \brief \b CLARFT forms the triangular factor T of a block reflector H = I - vtvH
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CLARFT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarft.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarft.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarft.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       RECURSIVE SUBROUTINE CLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
*
*       .. Scalar Arguments ..
*       CHARACTER          DIRECT, STOREV
*       INTEGER            K, LDT, LDV, N
*       ..
*       .. Array Arguments ..
*       COMPLEX            T( LDT, * ), TAU( * ), V( LDV, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CLARFT forms the triangular factor T of a complex block reflector H
*> of order n, which is defined as a product of k elementary reflectors.
*>
*> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
*>
*> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
*>
*> If STOREV = 'C', the vector which defines the elementary reflector
*> H(i) is stored in the i-th column of the array V, and
*>
*>    H  =  I - V * T * V**H
*>
*> If STOREV = 'R', the vector which defines the elementary reflector
*> H(i) is stored in the i-th row of the array V, and
*>
*>    H  =  I - V**H * T * V
*>
*> If DIRECT or STOREV = 'T', see Further Details for the shape of T
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] DIRECT
*> \verbatim
*>          DIRECT is CHARACTER*1
*>          Specifies the order in which the elementary reflectors are
*>          multiplied to form the block reflector:
*>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*>          = 'T' (with STOREV='R'): H = H(k) . . . H(2) H(1) (backward)
*>                but we return the T matrix that is already (conjugate) transposed
*> \endverbatim
*>
*> \param[in] STOREV
*> \verbatim
*>          STOREV is CHARACTER*1
*>          Specifies how the vectors which define the elementary
*>          reflectors are stored (see also Further Details):
*>          = 'C': column-wise
*>          = 'R': row-wise
*>          = 'T': (With DIRECT='F') Row-wise, but we return the T
*>                matrix that is already (conjugate) transposed.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the block reflector H. N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The order of the triangular factor T (= the number of
*>          elementary reflectors). K >= 1.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is COMPLEX array, dimension
*>                               (LDV,K) if STOREV = 'C'
*>                               (LDV,N) if STOREV = 'R'
*>          The matrix V. See further details.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V.
*>          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is COMPLEX array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i).
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          T is COMPLEX array, dimension (LDT,K)
*>          The k by k triangular factor T of the block reflector.
*>          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
*>          lower triangular. The rest of the array is not used.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T. LDT >= K.
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
*> \ingroup larft
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The shape of the matrix V and the storage of the vectors which define
*>  the H(i) is best illustrated by the following example with n = 5 and
*>  k = 3. The elements equal to 1 are not stored.
*>
*>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and (STOREV = 'R' or STOREV = 'T'):
*>
*>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
*>                   ( v1  1    )                     (     1 v2 v2 v2 )
*>                   ( v1 v2  1 )                     (        1 v3 v3 )
*>                   ( v1 v2 v3 )
*>                   ( v1 v2 v3 )
*>
*>  DIRECT = 'B' and STOREV = 'C':         (DIRECT = 'B' or DIRECT = 'T') and STOREV = 'R':
*>
*>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
*>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
*>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
*>                   (     1 v3 )
*>                   (        1 )
*>
*>  In addition, the shape of T is determined by these same flags as
*>   in the below table.
*>    'U' denotes upper triangular
*>    'L' denotes lower triangular
*>    'X' denotes no current implementation
*>    We also provide the logical variable that represents the case
*>    in the code if it is implemented
*>
*>    |-----------------------------------------------------------|
*>    |              | DIRECT = 'F' | DIRECT = 'B' | DIRECT = 'T' |
*>    |--------------+--------------+--------------+--------------|
*>    | STOREV = 'C' | U (QR)       | L (QL)       | X            |
*>    | STOREV = 'R' | U (LQ)       | L (RQ)       | U (RQT)      |
*>    | STOREV = 'T' | L (LQT)      | X            | X            |
*>    |-----------------------------------------------------------|
*>
*>    Finally, the relationship between the (conjugate) transposed T matrices
*>    are as follows: (Note that T_{FC} denotes the T associated with calling
*>    this routine with DIRECT = 'F' and STOREV = 'C')
*>
*>    T_{TR} = (T_{BR})**H
*>    T_{FT} = (T_{FR})**H
*> \endverbatim
*>
*  =====================================================================
      RECURSIVE SUBROUTINE CLARFT( DIRECT, STOREV, N, K, V, LDV,
     $                             TAU, T, LDT )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*        .. Scalar Arguments
*
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
*     ..
*     .. Array Arguments ..
*
      COMPLEX            T( LDT, * ), TAU( * ), V( LDV, * )
*     ..
*
*     .. Parameters ..
*
      COMPLEX            ONE, NEG_ONE
      PARAMETER(ONE=(1.0E+0,0.0E+0),
     $          NEG_ONE=(-1.0E+0,0.0E+0))
*
*     .. Local Scalars ..
*
      INTEGER           I,J,L,NX,INFO
      LOGICAL           QR, LQ, QL, RQ, LQT, RQT,
     $                  DIRF, DIRB, DIRT,
     $                  STOREC, STORER, STORET
*
*     .. External Subroutines ..
*
      EXTERNAL          CLARFT_LVL2, CTRMM, CGEMM, CLACPY, XERBLA
*
*     .. External Functions..
*
      LOGICAL           LSAME
      INTEGER           ILAENV
      EXTERNAL          LSAME, ILAENV
*
*     .. Intrinsic Functions..
*
      INTRINSIC         CONJG
*
*     The general scheme used is inspired by the approach inside DGEQRT3
*     which was (at the time of writing this code):
*     Based on the algorithm of Elmroth and Gustavson,
*     IBM J. Res. Develop. Vol 44 No. 4 July 2000.
*     ..
*     .. Executable Statements ..
*
*     Convert our character flags to logical flags for later
*
      DIRF = LSAME(DIRECT,'F')
      DIRB = LSAME(DIRECT,'B')
      DIRT = LSAME(DIRECT,'T')
      STOREC = LSAME(STOREV,'C')
      STORER = LSAME(STOREV,'R')
      STORET = LSAME(STOREV,'T')
*
*     Error handling for our character flags
*
      INFO = 0
      IF( .NOT.(DIRF.OR.DIRB.OR.DIRT) ) THEN
*
*        DIRECT holds an illegal value
*
         INFO = 1
      ELSE IF( .NOT.(STOREC.OR.STORER.OR.STORET) ) THEN
*
*        STOREV holds an illegal value
*
         INFO = 2
      ELSE IF( DIRB.AND.STORET ) THEN
*
*        This case is purposefully not implemented, but any other value for
*        STOREV is valid, so we report STOREV as the invalid input
*
         INFO = 2
      ELSE IF( DIRT.AND.STOREC ) THEN
*
*        This case is purposefully not implemented, but any other value for
*        DIRECT is valid, so we report DIRECT as the invalid input
*
         INFO = 1
      ELSE IF( DIRT.AND.STORET ) THEN
*
*        This case is purposefully not implemented, and is ambiguous what
*        the user wants to do, so we arbitrarily say DIRECT is the incorrect
*        character flag.
*
         INFO = 1
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA('CLARFT', INFO)
         RETURN
      END IF
*
*     Quick return if possible
*
      IF(N.EQ.0.OR.K.EQ.0) THEN
         RETURN
      END IF
*
*     Now we determine what factorization our flags are associated with
*
*     QR happens when we have forward direction in column storage
*
      QR = DIRF.AND.STOREC
*
*     LQT happens when we have forward direction in row storage and want to compute the transpose of
*     the T we would normally compute
*
      LQT = DIRF.AND.STORET
*
*     LQ happens when we have forward direction in row storage and want to compute the T we would
*     normally compute
*
      LQ = DIRF.AND.STORER
*
*     QL happens when we have backward direction in column storage
*
      QL = DIRB.AND.STOREC
*
*     RQT happens when we have backward direction in row storage and want to compute the transpose
*     of the T we would normally compute
*
      RQT = DIRT.AND.STORER
*
*     RQ happens when we have backward direction in row storage and want to compute the T that we
*     would normally compute
*
      RQ = DIRB.AND.STORER
*
*     Base case
*
      IF(N.EQ.1.OR.K.EQ.1) THEN
         IF( LQT.OR.RQT ) THEN
            T(1,1) = CONJG(TAU(1))
         ELSE
            T(1,1) = TAU(1)
         END IF
         RETURN
      END IF
*
*     Determine crossover point from level 2 to level 3 BLAS implementation
*
      NX = ILAENV(3, "CLARFT", DIRECT // STOREV, N, K, -1, -1)
      IF(K.LT.NX) THEN
*
*        Finish this component with a level 2 BLAS implementation
*
         CALL CLARFT_LVL2(DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT)
         RETURN
      END IF
*
*     Beginning of executable statements for the recursive case
*
      L = K / 2
      IF(QR) THEN
*
*        Break V apart into 6 components
*
*        V = |---------------|
*            |V_{1,1} 0      |
*            |V_{2,1} V_{2,2}|
*            |V_{3,1} V_{3,2}|
*            |---------------|
*
*        V_{1,1}\in\C^{l,l}      unit lower triangular
*        V_{2,1}\in\C^{k-l,l}    rectangular
*        V_{3,1}\in\C^{n-k,l}    rectangular
*
*        V_{2,2}\in\C^{k-l,k-l}  unit lower triangular
*        V_{3,2}\in\C^{n-k,k-l}  rectangular
*
*        We will construct the T matrix
*        T = |---------------|
*            |T_{1,1} T_{1,2}|
*            |0       T_{2,2}|
*            |---------------|
*
*        T is the triangular factor obtained from block reflectors.
*        To motivate the structure, assume we have already computed T_{1,1}
*        and T_{2,2}. Then collect the associated reflectors in V_1 and V_2
*
*        T_{1,1}\in\C^{l, l}     upper triangular
*        T_{2,2}\in\C^{k-l, k-l} upper triangular
*        T_{1,2}\in\C^{l, k-l}   rectangular
*
*        Where l = floor(k/2)
*
*        Then, consider the product:
*
*        (I - V_1*T_{1,1}*V_1')*(I - V_2*T_{2,2}*V_2')
*        = I - V_1*T_{1,1}*V_1' - V_2*T_{2,2}*V_2' + V_1*T_{1,1}*V_1'*V_2*T_{2,2}*V_2'
*
*        Define T_{1,2} = -T_{1,1}*V_1'*V_2*T_{2,2}
*
*        Then, we can define the matrix V as
*        V = |-------|
*            |V_1 V_2|
*            |-------|
*
*        So, our product is equivalent to the matrix product
*        I - V*T*V'
*        This means, we can compute T_{1,1} and T_{2,2}, then use this information
*        to compute T_{1,2}
*
*        Compute T_{1,1} recursively
*
         CALL CLARFT(DIRECT, STOREV, N, L, V, LDV, TAU, T, LDT)
*
*        Compute T_{2,2} recursively
*
         CALL CLARFT(DIRECT, STOREV, N-L, K-L, V(L+1,L+1), LDV,
     $      TAU(L+1), T(L+1,L+1), LDT)
*
*        Compute T_{1,2}
*        T_{1,2} = V_{2,1}'
*
         DO J = 1, L
            DO I = 1, K-L
               T(J,L+I) = CONJG(V(L+I,J))
            END DO
         END DO
*
*        T_{1,2} = T_{1,2}*V_{2,2}
*
         CALL CTRMM('Right', 'Lower', 'No transpose', 'Unit', L,
     $         K-L, ONE, V(L+1, L+1), LDV, T(1, L+1), LDT)

*
*        T_{1,2} = V_{3,1}'*V_{3,2} + T_{1,2}
*        Note: We assume K <= N, and GEMM will do nothing if N=K
*
         CALL CGEMM('Conjugate', 'No transpose', L, K-L, N-K, ONE,
     $         V(K+1, 1), LDV, V(K+1,L+1), LDV, ONE, T(1, L+1), LDT)
*
*        At this point, we have that T_{1,2} = V_1'*V_2
*        All that is left is to pre and post multiply by -T_{1,1} and T_{2,2}
*        respectively.
*
*        T_{1,2} = -T_{1,1}*T_{1,2}
*
         CALL CTRMM('Left', 'Upper', 'No transpose', 'Non-unit', L,
     $         K-L, NEG_ONE, T, LDT, T(1, L+1), LDT)
*
*        T_{1,2} = T_{1,2}*T_{2,2}
*
         CALL CTRMM('Right', 'Upper', 'No transpose', 'Non-unit', L,
     $         K-L, ONE, T(L+1,L+1), LDT, T(1, L+1), LDT)

      ELSE IF(LQ) THEN
*
*        Break V apart into 6 components
*
*        V = |----------------------|
*            |V_{1,1} V_{1,2} V{1,3}|
*            |0       V_{2,2} V{2,3}|
*            |----------------------|
*
*        V_{1,1}\in\C^{l,l}      unit upper triangular
*        V_{1,2}\in\C^{l,k-l}    rectangular
*        V_{1,3}\in\C^{l,n-k}    rectangular
*
*        V_{2,2}\in\C^{k-l,k-l}  unit upper triangular
*        V_{2,3}\in\C^{k-l,n-k}  rectangular
*
*        Where l = floor(k/2)
*
*        We will construct the T matrix
*        T = |---------------|
*            |T_{1,1} T_{1,2}|
*            |0       T_{2,2}|
*            |---------------|
*
*        T is the triangular factor obtained from block reflectors.
*        To motivate the structure, assume we have already computed T_{1,1}
*        and T_{2,2}. Then collect the associated reflectors in V_1 and V_2
*
*        T_{1,1}\in\C^{l, l}     upper triangular
*        T_{2,2}\in\C^{k-l, k-l} upper triangular
*        T_{1,2}\in\C^{l, k-l}   rectangular
*
*        Then, consider the product:
*
*        (I - V_1'*T_{1,1}'*V_1)*(I - V_2'*T_{2,2}'*V_2)
*        = I - V_1'*T_{1,1}'*V_1 - V_2'*T_{2,2}'*V_2 + V_1'*T_{1,1}'*V_1*V_2'*T_{2,2}'*V_2
*
*        Define T_{1,2} = -T_{1,1}*V_1*V_2'*T_{2,2}
*
*        Then, we can define the matrix V as
*        V = |---|
*            |V_1|
*            |V_2|
*            |---|
*
*        So, our product is equivalent to the matrix product
*        I - V'*T'*V
*        This means, we can compute T_{1,1} and T_{2,2}, then use this information
*        to compute T_{1,2}
*
*        Compute T_{1,1} recursively
*
         CALL CLARFT(DIRECT, STOREV, N, L, V, LDV, TAU, T, LDT)
*
*        Compute T_{2,2} recursively
*
         CALL CLARFT(DIRECT, STOREV, N-L, K-L, V(L+1,L+1), LDV,
     $      TAU(L+1), T(L+1,L+1), LDT)

*
*        Compute T_{1,2}
*        T_{1,2} = V_{1,2}
*
         CALL CLACPY('All', L, K - L, V(1,L+1), LDV, T(1, L+1), LDT)
*
*        T_{1,2} = T_{1,2}*V_{2,2}'
*
         CALL CTRMM('Right', 'Upper', 'Conjugate', 'Unit', L, K-L,
     $      ONE, V(L+1, L+1), LDV, T(1, L+1), LDT)

*
*        T_{1,2} = V_{1,3}*V_{2,3}' + T_{1,2}
*        Note: We assume K <= N, and GEMM will do nothing if N=K
*
         CALL CGEMM('No transpose', 'Conjugate', L, K-L, N-K, ONE,
     $      V(1, K+1), LDV, V(L+1, K+1), LDV, ONE, T(1, L+1), LDT)
*
*        At this point, we have that T_{1,2} = V_1*V_2'
*        All that is left is to pre and post multiply by -T_{1,1} and T_{2,2}
*        respectively.
*
*        T_{1,2} = -T_{1,1}*T_{1,2}
*
         CALL CTRMM('Left', 'Upper', 'No transpose', 'Non-unit', L,
     $      K-L, NEG_ONE, T, LDT, T(1, L+1), LDT)

*
*        T_{1,2} = T_{1,2}*T_{2,2}
*
         CALL CTRMM('Right', 'Upper', 'No transpose', 'Non-unit', L,
     $      K-L, ONE, T(L+1,L+1), LDT, T(1, L+1), LDT)
      ELSE IF(LQT) THEN
*
*        Break V apart into 6 components
*
*        V = |----------------------|
*            |V_{1,1} V_{1,2} V{1,3}|
*            |0       V_{2,2} V{2,3}|
*            |----------------------|
*
*        V_{1,1}\in\C^{l,l}      unit upper triangular
*        V_{1,2}\in\C^{l,k-l}    rectangular
*        V_{1,3}\in\C^{l,n-k}    rectangular
*
*        V_{2,2}\in\C^{k-l,k-l}  unit upper triangular
*        V_{2,3}\in\C^{k-l,n-k}  rectangular
*
*        Where l = floor(k/2)
*
*        We will construct the T matrix
*        T = |---------------|
*            |T_{1,1} 0      |
*            |T_{2,1} T_{2,2}|
*            |---------------|
*
*        T is the triangular factor obtained from block reflectors.
*        To motivate the structure, assume we have already computed T_{1,1}
*        and T_{2,2}. Then collect the associated reflectors in V_1 and V_2
*
*        T_{1,1}\in\C^{l, l}     lower triangular
*        T_{2,2}\in\C^{k-l, k-l} lower triangular
*        T_{2,1}\in\C^{k-l, l}   rectangular
*
*        Then, consider the product:
*
*        (I - V_1'*T_{1,1}*V_1)*(I - V_2'*T_{2,2}*V_2)
*        = I - V_1'*T_{1,1}*V_1 - V_2'*T_{2,2}*V_2 + V_1'*T_{1,1}*V_1*V_2'*T_{2,2}*V_2
*
*        Define T_{2,1} = -T_{2,2}*V_2*V_1'*T_{1,1}
*
*        Then, we can define the matrix V as
*        V = |---|
*            |V_1|
*            |V_2|
*            |---|
*
*        So, our product is equivalent to the matrix product
*        I - V'*T*V
*        This means, we can compute T_{1,1} and T_{2,2}, then use this information
*        to compute T_{2,1}
*
*        Compute T_{1,1} recursively
*
         CALL CLARFT(DIRECT, STOREV, N, L, V, LDV, TAU, T, LDT)
*
*        Compute T_{2,2} recursively
*
         CALL CLARFT(DIRECT, STOREV, N-L, K-L, V(L+1,L+1), LDV,
     $      TAU(L+1), T(L+1,L+1), LDT)

*
*        Compute T_{2,1}
*        T_{2,1} = V_{1,2}'
*
         DO I = 1, K-L
            DO J = 1, L
               T(L+I,J) = CONJG(V(J,L+I))
            END DO
         END DO
*
*        T_{2,1} = V_{2,2}*T_{2,1}
*
         CALL CTRMM('Left', 'Upper', 'No Transpose', 'Unit', K-L, L,
     $         ONE, V(L+1,L+1), LDV, T(L+1,1), LDT)
*
*        T_{2,1} = V_{2,3}*V_{1,3}' + T_{2,1}
*        Note: We assume K <= N, and GEMM will do nothing if N=K
*
         CALL CGEMM('No Transpose', 'Conjugate', K-L, L, N-K, ONE,
     $         V(L+1,K+1), LDV, V(1, K+1), LDV, ONE, T(L+1,1), LDT)
*
*        At this point, we have that T_{2,1} = V_2*V_1'
*        All that is left is to pre and post multiply by -T_{2,2} and T_{1,1}
*        respectively.
*
*        T_{2,1} = -T_{2,2}*T_{2,1}
*
         CALL CTRMM('Left', 'Lower', 'No Transpose', 'Non-unit',
     $         K-L, L, NEG_ONE, T(L+1,L+1), LDT, T(L+1,1), LDT)
*
*        T_{2,1} = T_{2,1}*T_{1,1}
*
         CALL CTRMM('Right', 'Lower', 'No Transpose', 'Non-unit',
     $         K-L, L, ONE, T, LDT, T(L+1,1), LDT)
      ELSE IF(QL) THEN
*
*        Break V apart into 6 components
*
*        V = |---------------|
*            |V_{1,1} V_{1,2}|
*            |V_{2,1} V_{2,2}|
*            |0       V_{3,2}|
*            |---------------|
*
*        V_{1,1}\in\C^{n-k,k-l}  rectangular
*        V_{2,1}\in\C^{k-l,k-l}  unit upper triangular
*
*        V_{1,2}\in\C^{n-k,l}    rectangular
*        V_{2,2}\in\C^{k-l,l}    rectangular
*        V_{3,2}\in\C^{l,l}      unit upper triangular
*
*        We will construct the T matrix
*        T = |---------------|
*            |T_{1,1} 0      |
*            |T_{2,1} T_{2,2}|
*            |---------------|
*
*        T is the triangular factor obtained from block reflectors.
*        To motivate the structure, assume we have already computed T_{1,1}
*        and T_{2,2}. Then collect the associated reflectors in V_1 and V_2
*
*        T_{1,1}\in\C^{k-l, k-l} non-unit lower triangular
*        T_{2,2}\in\C^{l, l}     non-unit lower triangular
*        T_{2,1}\in\C^{k-l, l}   rectangular
*
*        Where l = floor(k/2)
*
*        Then, consider the product:
*
*        (I - V_2*T_{2,2}*V_2')*(I - V_1*T_{1,1}*V_1')
*        = I - V_2*T_{2,2}*V_2' - V_1*T_{1,1}*V_1' + V_2*T_{2,2}*V_2'*V_1*T_{1,1}*V_1'
*
*        Define T_{2,1} = -T_{2,2}*V_2'*V_1*T_{1,1}
*
*        Then, we can define the matrix V as
*        V = |-------|
*            |V_1 V_2|
*            |-------|
*
*        So, our product is equivalent to the matrix product
*        I - V*T*V'
*        This means, we can compute T_{1,1} and T_{2,2}, then use this information
*        to compute T_{2,1}
*
*        Compute T_{1,1} recursively
*
         CALL CLARFT(DIRECT, STOREV, N-L, K-L, V, LDV, TAU, T,
     $         LDT)
*
*        Compute T_{2,2} recursively
*
         CALL CLARFT(DIRECT, STOREV, N, L, V(1, K-L+1), LDV,
     $      TAU(K-L+1), T(K-L+1,K-L+1), LDT)
*
*        Compute T_{2,1}
*        T_{2,1} = V_{2,2}'
*
         DO J = 1, K-L
            DO I = 1, L
               T(K-L+I,J) = CONJG(V(N-K+J, K-L+I))
            END DO
         END DO
*
*        T_{2,1} = T_{2,1}*V_{2,1}
*
         CALL CTRMM('Right', 'Upper', 'No transpose', 'Unit', L,
     $      K-L, ONE, V(N-K+1,1), LDV, T(K-L+1,1), LDT)

*
*        T_{2,1} = V_{2,2}'*V_{2,1} + T_{2,1}
*        Note: We assume K <= N, and GEMM will do nothing if N=K
*
         CALL CGEMM('Conjugate', 'No transpose', L, K-L, N-K, ONE,
     $      V(1,K-L+1), LDV, V, LDV, ONE, T(K-L+1,1), LDT)
*
*        At this point, we have that T_{2,1} = V_2'*V_1
*        All that is left is to pre and post multiply by -T_{2,2} and T_{1,1}
*        respectively.
*
*        T_{2,1} = -T_{2,2}*T_{2,1}
*
         CALL CTRMM('Left', 'Lower', 'No transpose', 'Non-unit', L,
     $      K-L, NEG_ONE, T(K-L+1,K-L+1), LDT, T(K-L+1,1), LDT)
*
*        T_{2,1} = T_{2,1}*T_{1,1}
*
         CALL CTRMM('Right', 'Lower', 'No transpose', 'Non-unit', L,
     $      K-L, ONE, T, LDT, T(K-L+1,1), LDT)
      ELSE IF(RQ) THEN
*        Break V apart into 6 components
*
*        V = |-----------------------|
*            |V_{1,1} V_{1,2} 0      |
*            |V_{2,1} V_{2,2} V_{2,3}|
*            |-----------------------|
*
*        V_{1,1}\in\C^{k-l,n-k}  rectangular
*        V_{1,2}\in\C^{k-l,k-l}  unit lower triangular
*
*        V_{2,1}\in\C^{l,n-k}    rectangular
*        V_{2,2}\in\C^{l,k-l}    rectangular
*        V_{2,3}\in\C^{l,l}      unit lower triangular
*
*        We will construct the T matrix
*        T = |---------------|
*            |T_{1,1} 0      |
*            |T_{2,1} T_{2,2}|
*            |---------------|
*
*        T is the triangular factor obtained from block reflectors.
*        To motivate the structure, assume we have already computed T_{1,1}
*        and T_{2,2}. Then collect the associated reflectors in V_1 and V_2
*
*        T_{1,1}\in\C^{k-l, k-l} non-unit lower triangular
*        T_{2,2}\in\C^{l, l}     non-unit lower triangular
*        T_{2,1}\in\C^{k-l, l}   rectangular
*
*        Where l = floor(k/2)
*
*        Then, consider the product:
*
*        (I - V_2'*T_{2,2}'*V_2)*(I - V_1'*T_{1,1}'*V_1)
*        = I - V_2'*T_{2,2}'*V_2 - V_1'*T_{1,1}'*V_1 + V_2'*T_{2,2}'*V_2*V_1'*T_{1,1}'*V_1
*
*        Define T_{2,1} = -T_{2,2}*V_2*V_1'*T_{1,1}
*
*        Then, we can define the matrix V as
*        V = |---|
*            |V_1|
*            |V_2|
*            |---|
*
*        So, our product is equivalent to the matrix product
*        I - V'*T'*V
*        This means, we can compute T_{1,1} and T_{2,2}, then use this information
*        to compute T_{2,1}
*
*        Compute T_{1,1} recursively
*
         CALL CLARFT(DIRECT, STOREV, N-L, K-L, V, LDV, TAU, T,
     $         LDT)
*
*        Compute T_{2,2} recursively
*
         CALL CLARFT(DIRECT, STOREV, N, L, V(K-L+1,1), LDV,
     $      TAU(K-L+1), T(K-L+1,K-L+1), LDT)
*
*        Compute T_{2,1}
*        T_{2,1} = V_{2,2}
*
         CALL CLACPY('All', L, K-L, V(K-L+1,N-K+1), LDV, T(K-L+1,1),
     $      LDT)

*
*        T_{2,1} = T_{2,1}*V_{1,2}'
*
         CALL CTRMM('Right', 'Lower', 'Conjugate', 'Unit', L, K-L,
     $      ONE, V(1, N-K+1), LDV, T(K-L+1,1), LDT)

*
*        T_{2,1} = V_{2,1}*V_{1,1}' + T_{2,1}
*        Note: We assume K <= N, and GEMM will do nothing if N=K
*
         CALL CGEMM('No transpose', 'Conjugate', L, K-L, N-K, ONE,
     $      V(K-L+1,1), LDV, V, LDV, ONE, T(K-L+1,1), LDT)

*
*        At this point, we have that T_{2,1} = V_2*V_1'
*        All that is left is to pre and post multiply by -T_{2,2} and T_{1,1}
*        respectively.
*
*        T_{2,1} = -T_{2,2}*T_{2,1}
*
         CALL CTRMM('Left', 'Lower', 'No transpose', 'Non-unit', L,
     $      K-L, NEG_ONE, T(K-L+1,K-L+1), LDT, T(K-L+1,1), LDT)

*
*        T_{2,1} = T_{2,1}*T_{1,1}
*
         CALL CTRMM('Right', 'Lower', 'No transpose', 'Non-unit', L,
     $      K-L, ONE, T, LDT, T(K-L+1,1), LDT)
      ELSE IF(RQT) THEN
*        Break V apart into 6 components
*
*        V = |-----------------------|
*            |V_{1,1} V_{1,2} 0      |
*            |V_{2,1} V_{2,2} V_{2,3}|
*            |-----------------------|
*
*        V_{1,1}\in\C^{k-l,n-k}  rectangular
*        V_{1,2}\in\C^{k-l,k-l}  unit lower triangular
*
*        V_{2,1}\in\C^{l,n-k}    rectangular
*        V_{2,2}\in\C^{l,k-l}    rectangular
*        V_{2,3}\in\C^{l,l}      unit lower triangular
*
*        We will construct the T matrix
*        T = |---------------|
*            |T_{1,1} T_{1,2}|
*            | 0      T_{2,2}|
*            |---------------|
*
*        T is the triangular factor obtained from block reflectors.
*        To motivate the structure, assume we have already computed T_{1,1}
*        and T_{2,2}. Then collect the associated reflectors in V_1 and V_2
*
*        T_{1,1}\in\C^{k-l, k-l} non-unit upper triangular
*        T_{2,2}\in\C^{l, l}     non-unit upper triangular
*        T_{1,2}\in\C^{k-l, l}   rectangular
*
*        Where l = floor(k/2)
*
*        Then, consider the product:
*
*        (I - V_2'*T_{2,2}*V_2)*(I - V_1'*T_{1,1}*V_1)
*        = I - V_2'*T_{2,2}*V_2 - V_1'*T_{1,1}*V_1 + V_2'*T_{2,2}*V_2*V_1'*T_{1,1}*V_1
*
*        Define T_{1,2} = -T_{1,1}*V_1*V_2'*T_{2,2}
*
*        Then, we can define the matrix V as
*        V = |---|
*            |V_1|
*            |V_2|
*            |---|
*
*        So, our product is equivalent to the matrix product
*        I - V'*T*V
*        This means, we can compute T_{1,1} and T_{2,2}, then use this information
*        to compute T_{1,2}
*
*        Compute T_{1,1} recursively
*
         CALL CLARFT(DIRECT, STOREV, N-L, K-L, V, LDV, TAU, T,
     $         LDT)
*
*        Compute T_{2,2} recursively
*
         CALL CLARFT(DIRECT, STOREV, N, L, V(K-L+1,1), LDV,
     $      TAU(K-L+1), T(K-L+1,K-L+1), LDT)
*
*        Compute T_{1,2}
*        T_{1,2} = V_{2,2}'
*
         DO I = 1, K-L
            DO J = 1, L
               T(I,K-L+J) = CONJG(V(K-L+J, N-K+I))
            END DO
         END DO
*
*        T_{1,2} = V_{1,2}T_{1,2}
*
         CALL CTRMM('Left', 'Lower', 'No Transpose', 'Unit', K-L, L,
     $         ONE, V(1,N-K+1), LDV, T(1,K-L+1), LDT)
*
*        T_{1,2} = V_{1,1}V_{2,1}' + T_{1,2}
*
         CALL CGEMM('No Transpose', 'Conjugate', K-L, L, N-K, ONE, V,
     $         LDV, V(K-L+1,1), LDV, ONE, T(1, K-L+1), LDT)
*
*        At this point, we have that T_{1,2} = V_1*V_2'
*        All that is left is to pre and post multiply by -T_{1,1} and T_{2,2}
*        respectively.
*
*        T_{1,2} = -T_{1,1}*T_{1,2}
*
         CALL CTRMM('Left', 'Upper', 'No Transpose', 'Non-Unit',
     $         K-L, L, NEG_ONE, T, LDT, T(1, K-L+1), LDT)
*
*        T_{1,2} = T_{1,2}*T_{2,2}
*
         CALL CTRMM('Right', 'Upper', 'No Transpose', 'Non-Unit',
     $         K-L, L, ONE, T(K-L+1,K-L+1), LDT, T(1, K-L+1), LDT)
      END IF
      END SUBROUTINE
