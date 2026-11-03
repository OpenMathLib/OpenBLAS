*> \brief \b CLARFT_LVL2 forms the triangular factor T of a block reflector H = I - vtvH
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CLARFT_LVL2( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
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
*> CLARFT_LVL2 forms the triangular factor T of a complex block reflector H
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
*>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
*>
*>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
*>                   ( v1  1    )                     (     1 v2 v2 v2 )
*>                   ( v1 v2  1 )                     (        1 v3 v3 )
*>                   ( v1 v2 v3 )
*>                   ( v1 v2 v3 )
*>
*>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
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
      SUBROUTINE CLARFT_LVL2( DIRECT, STOREV, N, K, V, LDV, TAU,
     $            T, LDT )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
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
      COMPLEX            ONE
      PARAMETER(ONE=(1.0E+0,0.0E+0))
*
*     .. Local Scalars ..
*
      INTEGER           I,J,KMI,NMI,INFO
      LOGICAL           QR, LQ, QL, RQ, LQT, RQT,
     $                  DIRF, DIRB, DIRT,
     $                  STOREC, STORER, STORET
*
*     .. External Subroutines ..
*
      EXTERNAL          CTRMV,CGEMV,CGEMM,XERBLA
*
*     .. External Functions..
*
      LOGICAL           LSAME
      EXTERNAL          LSAME
*
*     .. Intrinsic Functions..
*
      INTRINSIC         CONJG
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
         CALL XERBLA('CLARFT_LVL2', INFO)
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
      IF( N.EQ.1.OR.K.EQ.1) THEN
         IF( LQT.OR.RQT ) THEN
            T(1,1) = CONJG(TAU(1))
         ELSE
            T(1,1) = TAU(1)
         END IF
         RETURN
      END IF
      IF (QR) THEN
*
*        Break V into 9 components
*
*        V = |-----------------------|
*            |V_{1,1} 0       0      | i-1
*            |V_{2,1} V_{2,2} 0      | 1
*            |V_{3,1} V_{3,2} V_{3,3}| n-i
*            |-----------------------|
*             i-1     1       k-i
*
*        V_{1,1}, V_{2,2} and V_{3,3} are unit lower triangular
*
*        This is how we are going to view the matrix V at each step
*        i=2,\dots,k, then we grow into V_{3,3} and repeat until we
*        reach the end. On each iteration V_{3,3} is not referenced
*
*        We will construct T one column at a time from left to right
*        after initializing T(1,1) = TAU(1)
*
*        T = |-------------------------|
*            | T_{1,1} T_{1,2} T_{1,3} | i-1
*            | 0       T_{2,2} T_{2,3} | 1
*            | 0       0       T_{3,3} | k-i
*            |-------------------------|
*              i-1     1       k-i
*
*        T_{1,1}, T_{2,2}, and T_{3,3} are non-unit lower triangular
*
*        Similarly as above, we will construct T_{1,2} and T_{2,2} at
*        each iteration i = 2, \dots k, and then grow into T_{1:3,3}. On
*        each iteration, T_{1:3,3} are not referenced. See clarft.f
*        for details on how these formulae were constructed.
*
*        We now get
*
*        T_{1,2} = -T_{1,1}[V_{1,1}\\V_{2,1}\\V_{3,1}]'
*                       [0\\V_{2,2}\\V_{3,2}]T_{2,2}
*
*        T_{1,2} = -T_{1,1}(V_{2,1}' + V_{3,1}'V_{3,2})T_{2,2}
*
*        This means we will do the following
*
*        T_{1,2} = -V_{2,1}'T_{2,2} = -\tau_{i}V_{2,1}'
*        T_{1,2} = -\tau_{i}V_{3,2}' V_{3,1} + T_{1,2}
*        T_{1,2} = T_{1,1}T_{1,2}
*        T_{2,2} = \tau{i}
*
         T(1,1) = TAU(1)

         DO I = 2, K
*
*           T_{1,2} = -V_{2,1}'V_{2,2}T_{2,2} = -\tau_i V_{2,1}'
*           We must do this at copy time as otherwise gemv will do nothing
*           on the last column when n=k, but we neet to make sure we are
*           scaled by this value
*
            DO J = 1, I-1
               T(J,I) = -TAU(I)*CONJG(V(I,J))
            END DO

*
*           T_{1,2} = -V_{3,1}'V_{3,2}T_{2,2} + T_{1,2}
*                   = -\tau{i} V_{3,2}'V_{3,1} + T_{1,2}
*
            CALL CGEMV('Conjugate Transpose', N-I, I-1, -TAU(I),
     $            V(I+1,1), LDV, V(I+1,I), 1, ONE, T(1, I), 1)


*
*           T_{1,2} = T_{1,1}T_{1,2}
*
            CALL CTRMV('Upper', 'No Transpose', 'Non-unit', I-1,
     $            T, LDT, T(1,I), 1)

*
*           T_{2,2} = \tau{i}
*
            T(I,I) = TAU(I)
         END DO
      ELSE IF (LQ) THEN
*
*        Break V into 9 components
*
*        V = |-------------------------|
*            | V_{1,1} V_{1,2} V_{1,3} | i-1
*            | 0       V_{2,2} V_{2,3} | 1
*            | 0       0       V_{3,3} | k-i
*            |-------------------------|
*              i-1     1       n-i
*
*        V_{1,1}, V_{2,2} and V_{3,3} are unit upper triangular
*
*        This is how we are going to view the matrix V at each step
*        i=2,\dots,k, then we grow into V_{3,3} and repeat until we
*        reach the end. On each iteration V_{3,3} is not referenced
*
*        We will construct T one column at a time from left to right
*        after initializing T(1,1) = TAU(1)
*
*        T = |-------------------------|
*            | T_{1,1} T_{1,2} T_{1,3} | i-1
*            | 0       T_{2,2} T_{2,3} | 1
*            | 0       0       T_{3,3} | k-i
*            |-------------------------|
*              i-1     1       k-i
*
*        Similarly as above, we will construct T_{1,2} and T_{2,2} at
*        each iteration i = 2, \dots k, and then grow into T_{1:3,3}. On
*        each iteration, T_{1:3,3} are not referenced. See clarft.f
*        for details on how these formulae were constructed.
*
*        We now get
*
*        T_{1,2} = -T_{1,1}[V_{1,1} V_{1,2} V_{1,3}][ 0 V_{2,2} V_{2,3} ]'T_{2,2}
*
*        T_{1,2} = -T_{1,1}(V_{1,2} + V_{1,3}V_{2,3}')T_{2,2}
*
*        This means we will do the following
*
*        T_{1,2} = -V_{1,2}T_{2,2} = -\tau_{i}V_{1,2}
*        T_{1,2} = -\tau_{i}V_{1,3}V_{2,3}' + T_{1,2}
*        T_{1,2} = T_{1,1}T_{1,2}
*        T_{2,2} = \tau{i}
*
         T(1,1) = TAU(1)

         DO I = 2, K
*
*           T_{1,2} = -\tau_{i}V_{1,2}
*
            DO J = 1, I-1
               T(J, I) = -TAU(I)*V(J, I)
            END DO

*
*           T_{1,2} = -\tau_{i}V_{1,3}V_{2,3}' + T_{1,2}
*
            CALL CGEMM('No Transpose', 'Conjugate Transpose', I-1,
     $            1, N-I, -TAU(I), V(1,I+1), LDV, V(I, I+1), LDV, ONE,
     $            T(1, I), LDT)
*
*           T_{1,2} = T_{1,1}T_{1,2}
*
            CALL CTRMV('Upper', 'No Transpose', 'Non-unit', I-1,
     $            T, LDT, T(1,I), 1)

*
*           T_{2,2} = \tau{i}
*
            T(I,I) = TAU(I)
         END DO
      ELSE IF (LQT) THEN
*
*        Break V into 9 components
*
*        V = |-------------------------|
*            | V_{1,1} V_{1,2} V_{1,3} | i-1
*            | 0       V_{2,2} V_{2,3} | 1
*            | 0       0       V_{3,3} | k-i
*            |-------------------------|
*              i-1     1       n-i
*
*        V_{1,1}, V_{2,2} and V_{3,3} are unit upper triangular
*
*        This is how we are going to view the matrix V at each step
*        i=2,\dots,k, then we grow into V_{3,3} and repeat until we
*        reach the end. On each iteration V_{3,3} is not referenced
*
*        We will construct T one column at a time from left to right
*        after initializing T(1,1) = TAU(1)
*
*        T = |-------------------------|
*            | T_{1,1} 0       0       | i-1
*            | T_{2,1} T_{2,2} 0       | 1
*            | T_{3,1} T_{3,2} T_{3,3} | k-i
*            |-------------------------|
*              i-1     1       k-i
*
*        Similarly as above, we will construct T_{2,1} and T_{2,2} at
*        each iteration i = 2, \dots k, and then grow into T_{3,1:3}. On
*        each iteration, T_{3,1:3} are not referenced. See clarft.f
*        for details on how these formulae were constructed.
*
*        We now get
*
*        T_{2,1} = -T_{2,2}[0 V_{2,2} V_{2,3}][V_{1,1} V_{1,2} V_{1,3}]'T_{1,1}
*
*        T_{2,1} = -T_{2,2}(V_{1,2}' + V_{2,3}V_{1,3}')T_{1,1}
*
*        This means we will do the following
*
*        T_{2,1} = -T_{2,2}V_{1,2}' = -\tau_{i}V_{1,2}'
*        T_{2,1} = -\tau_{i}V_{1,3}V_{2,3}' + T_{2,1}
*        T_{2,1} = T_{1,1}'T_{2,1}
*        T_{2,2} = \tau{i}
*
         T(1,1) = CONJG(TAU(1))

         DO I = 2, K
*
*           T_{2,1} = -\tau_{i}V_{1,2}'
*
            DO J = 1, I-1
               T(I,J) = -CONJG(TAU(I)*V(J,I))
            END DO
*
*           T_{2,1} = -\tau_{i}V_{2,3}V_{1,3}' + T_{2,1}
*
            CALL CGEMM('No Transpose', 'Conjugate Transpose', 1,
     $            I-1, N-I, -CONJG(TAU(I)), V(I,I+1), LDV, V(1, I+1),
     $            LDV, ONE, T(I, 1), LDT)
*
*           T_{2,1} = T_{1,1}'T_{2,1}
*
            CALL CTRMV('Lower', 'Transpose', 'Non-unit',
     $            I-1, T, LDT, T(I,1), LDT)

            T(I,I) = CONJG(TAU(I))
         END DO
      ELSE IF (QL) THEN
*
*     Break V into 9 components
*
*     V = |-------------------------|
*         | V_{1,1} V_{1,2} V_{1,3} | n-i
*         | 0       V_{2,2} V_{2,3} | 1
*         | 0       0       V_{3,3} | i-1
*         |-------------------------|
*           k-i     1       i-1
*
*        V_{1,1}, V_{2,2} and V_{3,3} are unit upper triangular
*
*        This is how we are going to view the matrix V at each step
*        i=2,\dots,k, then we grow into V_{1,1} and repeat until we
*        reach the end. On each iteration V_{1,1} is not referenced
*
*        We will construct T one column at a time from right to left
*        after initializing T(K,K) = TAU(K)
*
*     T = |-------------------------|
*         | T_{1,1} 0       0       | k-i
*         | T_{2,1} T_{2,2} 0       | 1
*         | T_{3,1} T_{3,2} T_{3,3} | i-1
*         |-------------------------|
*           k-i     1       i-1
*
*        T_{1,1}, T_{2,2}, and T_{3,3} are non-unit lower triangular
*
*        Similarly as above, we will construct T_{2,2} and T_{3,2} at
*        each iteration i = 2, \dots k, and then grow into T_{1:3,1}. On
*        each iteration, T_{1:3,1} are not referenced. See clarft.f
*        for details on how these formulae were constructed.
*
*        We get that
*
*        T_{3,2} = -T_{3,3}[V_{1,3}\\V_{2,3}\\V_{3,3}]'
*        [V_{1,2}\\V_{2,2}\\0]T_{2,2}
*
*        T_{3,2} = -T_{3,3}(V_{1,3}'V_{1,2} + V_{2,3}')T_{2,2}
*
*        Thus, we will compute
*
*        T_{2,2} = \tau_{k-i+1}
*        T_{3,2} = -T_{3,3}V_{3,2}' = -\tau_{k-i+1}V_{3,2}'
*        T_{3,2} = -\tau_{k-i+1}V_{1,3}'V_{1,2} + T_{3,2}
*        T_{3,2} = T_{3,3}T_{3,2}
*
         T(K,K) = TAU(K)
         DO I = 2, K
            KMI = K-I+1
            NMI = N-I+1
*
*             T_{2,2} = \tau(k-i+1)
*
            T(KMI,KMI) = TAU(KMI)
*
*             T_{3,2} = -\tau(k-i+1)V_{2,3}'
*
            DO J = 1, I-1
               T(KMI + J, KMI) = -TAU(KMI)*CONJG(V(NMI, KMI + J))
            END DO
*
*             T_{3,2} = -\tau(k-i+1)V_{1,3}'V_{1,2} + T_{3,2}
*
            CALL CGEMV('Conjugate Transpose', N-I, I-1, -TAU(KMI),
     $            V(1, KMI + 1), LDV, V(1, KMI), 1, ONE,
     $            T(KMI+1, KMI), 1)
*
*             T_{3,2} = T_{3,3}T_{3,2}
*
            CALL CTRMV('Lower', 'No Transpose', 'Non-unit', I-1,
     $            T(KMI + 1, KMI + 1), LDT, T(KMI + 1, KMI), 1)

         END DO
      ELSE IF (RQ) THEN
*
*     Break V into 9 components
*
*     V = |-------------------------|
*         | V_{1,1} 0       0       | k-i
*         | V_{2,1} V_{2,2} 0       | 1
*         | V_{3,1} V_{3,2} V_{3,3} | i-1
*         |-------------------------|
*           n-i     1       i-1
*
*        V_{1,1}, V_{2,2} and V_{3,3} are unit lower triangular
*
*        This is how we are going to view the matrix V at each step
*        i=2,\dots,k, then we grow into V_{1,1} and repeat until we
*        reach the end. On each iteration V_{1,1} is not referenced
*
*        We will construct T one column at a time from right to left
*        after initializing T(K,K) = TAU(K)
*
*     T = |-------------------------|
*         | T_{1,1} 0       0       | k-i
*         | T_{2,1} T_{2,2} 0       | 1
*         | T_{3,1} T_{3,2} T_{3,3} | i-1
*         |-------------------------|
*           k-i     1       i-1
*
*        T_{1,1}, T_{2,2}, and T_{3,3} are non-unit lower triangular
*
*        Similarly as above, we will construct T_{2,2} and T_{3,2} at
*        each iteration i = 2, \dots k, and then grow into T_{1:3,1}. On
*        each iteration, T_{1:3,1} are not referenced. See clarft.f
*        for details on how these formulae were constructed.
*
*        We get that
*
*        T_{3,2} = -T_{3,3}[V_{3,1} V_{3,2} V_{3,3}][V_{2,1} V_{2,2} 0]'T_{2,2}
*
*        T_{3,2} = -T_{3,3}(V_{3,1}V_{2,1}' + V_{3,2})T_{2,2}
*
*        Thus, we will compute
*
*        T_{2,2} = \tau_{k-i+1}
*        T_{3,2} = -\tau_{k-i+1}V_{3,2}
*        T_{3,2} = -\tau_{k-i+1}V_{3,1}V_{2,1}' + T_{3,2}
*        T_{3,2} = T_{3,3}T_{3,2}
*
         T(K,K) = TAU(K)
         DO I = 2, K
            KMI = K-I+1
            NMI = N-I+1
*
*           T_{2,2} = \tau_{k-i+1}
*
            T(KMI,KMI) = TAU(KMI)
*
*           T_{3,2} = -\tau_{k-i+1}V_{3,2}
*
            DO J = 1, I-1
               T(KMI + J, KMI) = -TAU(KMI)*V(KMI + J, NMI)
            END DO
*
*           T_{3,2} = -\tau_{k-i+1}V_{3,1}V_{2,1}' + T_{3,2}
*
            CALL CGEMM('No Transpose', 'Conjugate Transpose', I-1,
     $            1, N-I, -TAU(KMI), V(KMI+1, 1), LDV, V(KMI, 1), LDV,
     $            ONE, T(KMI+1, KMI), LDT)
*
*           T_{3,2} = T_{3,3}T_{3,2}
*
            CALL CTRMV('Lower', 'No Transpose', 'Non-unit', I-1,
     $            T(KMI+1, KMI+1), LDT, T(KMI+1, KMI), 1)
         END DO
      ELSE IF (RQT) THEN
*
*     Break V into 9 components
*
*     V = |-------------------------|
*         | V_{1,1} 0       0       | k-i
*         | V_{2,1} V_{2,2} 0       | 1
*         | V_{3,1} V_{3,2} V_{3,3} | i-1
*         |-------------------------|
*           n-i     1       i-1
*
*        V_{1,1}, V_{2,2} and V_{3,3} are unit lower triangular
*
*        This is how we are going to view the matrix V at each step
*        i=2,\dots,k, then we grow into V_{1,1} and repeat until we
*        reach the end. On each iteration V_{1,1} is not referenced
*
*        We will construct T one column at a time from right to left
*        after initializing T(K,K) = TAU(K)
*
*     T = |-------------------------|
*         | T_{1,1} T_{1,2} T_{1,3} | k-i
*         | 0       T_{2,2} T_{2,3} | 1
*         | 0       0       T_{3,3} | i-1
*         |-------------------------|
*           k-i     1       i-1
*
*        T_{1,1}, T_{2,2}, and T_{3,3} are non-unit lower triangular
*
*        Similarly as above, we will construct T_{2,2} and T_{2,3} at
*        each iteration i = 2, \dots k, and then grow into T_{1,1:3}. On
*        each iteration, T_{1,1:3} are not referenced. See clarft.f
*        for details on how these formulae were constructed.
*
*        We get that
*
*        T_{2,3} = -T_{2,2}[V_{2,1} V_{2,2} 0][V_{3,1} V_{3,2} V_{3,3}]'T_{3,3}
*
*        T_{3,2} = -T_{2,2}(V_{2,1}V_{3,1}' + V_{3,2}')T_{3,3}
*
*        Thus, we will compute
*
*        T_{2,2} = \tau_{k-i+1}
*        T_{2,3} = -\tau_{k-i+1}V_{3,2}'
*        T_{2,3} = -\tau_{k-i+1}V_{2,1}V_{3,1}' + T_{2,3}
*        T_{2,3} = T_{3,3}'T_{2,3}
*
         T(K,K) = CONJG(TAU(K))
         DO I = 2, K
            KMI = K-I+1
            NMI = N-I+1
*
*           T_{2,2} = \tau_{k-i+1}
*
            T(KMI,KMI) = CONJG(TAU(KMI))
*
*           T_{2,3} = -\tau_{k-i+1}V_{3,2}'
*
            DO J = 1, I-1
               T(KMI, KMI + J) = -CONJG(TAU(KMI)*V(KMI + J, NMI))
            END DO
*
*           T_{2,3} = -\tau_{k-i+1}V_{2,1}V_{3,1}' + T_{2,3}
*
            CALL CGEMM('No Transpose', 'Conjugate Transpose', 1,
     $            I-1, N-I, -CONJG(TAU(KMI)), V(KMI, 1), LDV,
     $            V(KMI+1,1), LDV, ONE, T(KMI, KMI+1), LDT)
*
*           T_{2,3} = T_{3,3}'T_{2,3}
*
            CALL CTRMV('Upper', 'Transpose', 'Non-unit', I-1,
     $            T(KMI+1, KMI+1), LDT, T(KMI, KMI+1), LDT)
         END DO
      END IF
      END SUBROUTINE
