*> \brief \b DSTEBR
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSTEBR( N, D, E, WORK, LWORK, IWORK, LIWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LIWORK, LWORK, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       DOUBLE PRECISION   D( * ), E( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSTEBR computes all eigenvalues of a symmetric tridiagonal matrix
*> using a values-only divide-and-conquer method that propagates only
*> the first and last boundary rows of each local eigenvector block
*> (STEBR = Symmetric tridiagonal Eigenvalues, Boundary-Row DC).
*>
*> Unlike DSTEDC with COMPZ = 'N', DSTEBR does not store full secular
*> eigenvector blocks across merge levels for later replay.  It is a
*> distinct routine and does not replace DSTEDC or DSTERF.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] D
*> \verbatim
*>          D is DOUBLE PRECISION array, dimension (N)
*>          On entry, the n diagonal elements of the tridiagonal matrix.
*>          On exit, if INFO = 0, the eigenvalues in ascending order.
*> \endverbatim
*>
*> \param[in,out] E
*> \verbatim
*>          E is DOUBLE PRECISION array, dimension (N-1)
*>          On entry, the (n-1) subdiagonal elements of the tridiagonal
*>          matrix.
*>          On exit, E has been destroyed.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>          If N <= SMLSIZ, where SMLSIZ is returned by ILAENV and is
*>          typically about 25, LWORK must be at least 1.
*>          If N > SMLSIZ, LWORK must be at least
*>          MAX( 16*N, 2*N + SMLSIZ**2 + 2*SMLSIZ - 2 ).
*>
*>          If LWORK = -1, then a workspace query is assumed; the
*>          routine only calculates the optimal sizes of the WORK and
*>          IWORK arrays, returns these values as the first entries of
*>          the WORK and IWORK arrays, and no error message related to
*>          LWORK or LIWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
*>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*> \endverbatim
*>
*> \param[in] LIWORK
*> \verbatim
*>          LIWORK is INTEGER
*>          The dimension of the array IWORK.
*>          If N <= SMLSIZ, LIWORK must be at least 1.
*>          If N > SMLSIZ, LIWORK must be at least 7*N.
*>
*>          If LIWORK = -1, then a workspace query is assumed; the
*>          routine only calculates the optimal sizes of the WORK and
*>          IWORK arrays, returns these values as the first entries of
*>          the WORK and IWORK arrays, and no error message related to
*>          LWORK or LIWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          > 0:  The algorithm failed to compute an eigenvalue while
*>                working on a submatrix lying in rows and columns
*>                INFO/(N+1) through MOD(INFO,N+1).
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
*> \ingroup stebr
*
*> \par Contributors:
*  ==================
*>
*> Ruiyi Zhan and Shaoshuai Zhang, University of Electronic Science
*> and Technology of China
*
*  =====================================================================
      SUBROUTINE DSTEBR( N, D, E, WORK, LWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK computational routine --
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LIWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            FINISH, LIWMIN, LWMIN, M, SMLSIZ, START
      DOUBLE PRECISION   EPS, ORGNRM, TINY
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           ILAENV, DLAMCH, DLANST
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAED0_BR, DLASCL, DLASRT, DSTERF, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MOD, SQRT
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
      IF( N.LT.0 ) THEN
         INFO = -1
      END IF
*
      IF( INFO.EQ.0 ) THEN
         SMLSIZ = ILAENV( 9, 'DSTEBR', ' ', 0, 0, 0, 0 )
         IF( N.LE.SMLSIZ ) THEN
            LWMIN = 1
            LIWMIN = 1
         ELSE
*
*           BLO/BHI: 2*N; merge scratch: 14*N (see DLAED7_BR).
*           Leaf solves need 2*N + SMLSIZ**2 + 2*SMLSIZ - 2: the
*           SMLSIZ-by-SMLSIZ leaf eigenvector block plus the
*           DSTEQR work array.
*
            LWMIN = MAX( 16*N, 2*N + SMLSIZ*SMLSIZ + 2*SMLSIZ - 2 )
            LIWMIN = MAX( 1, 7*N )
         END IF
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
*
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -5
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -7
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEBR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      IF( N.EQ.0 )
     $   RETURN
      IF( N.EQ.1 )
     $   RETURN
*
      IF( N.LE.SMLSIZ ) THEN
         CALL DSTERF( N, D, E, INFO )
         GO TO 50
      END IF
*
      ORGNRM = DLANST( 'M', N, D, E )
      IF( ORGNRM.LE.ZERO )
     $   GO TO 50
*
      EPS = DLAMCH( 'Epsilon' )
      START = 1
*
   10 CONTINUE
      IF( START.LE.N ) THEN
         FINISH = START
   20    CONTINUE
         IF( FINISH.LT.N ) THEN
            TINY = EPS*SQRT( ABS( D( FINISH ) ) )*
     $             SQRT( ABS( D( FINISH+1 ) ) )
            IF( ABS( E( FINISH ) ).GT.TINY ) THEN
               FINISH = FINISH + 1
               GO TO 20
            END IF
         END IF
*
         M = FINISH - START + 1
         IF( M.EQ.1 ) THEN
            START = FINISH + 1
            GO TO 10
         END IF
*
         IF( M.GT.SMLSIZ ) THEN
            ORGNRM = DLANST( 'M', M, D( START ), E( START ) )
            CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, M, 1, D( START ),
     $                   M, INFO )
            CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, M-1, 1,
     $                   E( START ), M-1, INFO )
*
            CALL DLAED0_BR( M, D( START ), E( START ), WORK, IWORK,
     $                      INFO )
            IF( INFO.NE.0 ) THEN
               INFO = ( INFO / ( M+1 )+START-1 )*( N+1 ) +
     $                MOD( INFO, ( M+1 ) ) + START - 1
               GO TO 50
            END IF
*
            CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, M, 1, D( START ),
     $                   M, INFO )
         ELSE
            CALL DSTERF( M, D( START ), E( START ), INFO )
            IF( INFO.NE.0 ) THEN
               INFO = START*( N+1 ) + FINISH
               GO TO 50
            END IF
         END IF
*
         START = FINISH + 1
         GO TO 10
      END IF
*
      CALL DLASRT( 'I', N, D, INFO )
*
   50 CONTINUE
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
      RETURN
*
*     End of DSTEBR
*
      END
