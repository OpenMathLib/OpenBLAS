*> \brief \b ZDRVMR3
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZDRVMR3( NSIZES, MM, NN, NTYPES, DOTYPE, ISEED, THRESH,
*                          A, LDA, U, LDU, VT, LDVT, ASAV, USAV, VTSAV, S,
*                          SSAV, E, WORK, LWORK, RWORK, IWORK, NOUNIT,
*                          INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, LDU, LDVT, LWORK, NOUNIT, NSIZES,
*      $                   NTYPES
*       DOUBLE PRECISION   THRESH
*       ..
*       .. Array Arguments ..
*       LOGICAL            DOTYPE( * )
*       INTEGER            ISEED( 4 ), IWORK( * ), MM( * ), NN( * )
*       DOUBLE PRECISION   E( * ), RWORK( * ), S( * ), SSAV( * )
*       COMPLEX*16         A( LDA, * ), ASAV( LDA, * ), U( LDU, * ),
*      $                   USAV( LDU, * ), VT( LDVT, * ),
*      $                   VTSAV( LDVT, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZDRVMR3 checks the singular value decomposition (SVD) driver ZGESVD,
*> ZGESDD, ZGESVJ, ZGEJSV, ZGESVDMR3, and ZGESVDQ.
*>
*> ZGESVD and ZGESDD factors A = U diag(S) VT, where U and VT are
*> unitary and diag(S) is diagonal with the entries of the array S on
*> its diagonal. The entries of S are the singular values, nonnegative
*> and stored in decreasing order.  U and VT can be optionally not
*> computed, overwritten on A, or computed partially.
*>
*> A is M by N. Let MNMIN = min( M, N ). S has dimension MNMIN.
*> U can be M by M or M by MNMIN. VT can be N by N or MNMIN by N.
*>
*> When ZDRVMR3 is called, a number of matrix "sizes" (M's and N's)
*> and a number of matrix "types" are specified.  For each size (M,N)
*> and each type of matrix, and for the minimal workspace as well as
*> workspace adequate to permit blocking, an  M x N  matrix "A" will be
*> generated and used to test the SVD routines.  For each matrix, A will
*> be factored as A = U diag(S) VT and the following 12 tests computed:
*>
*> Test for ZGESVD:
*>
*> (1)   | A - U diag(S) VT | / ( |A| max(M,N) ulp )
*>
*> (2)   | I - U'U | / ( M ulp )
*>
*> (3)   | I - VT VT' | / ( N ulp )
*>
*> (4)   S contains MNMIN nonnegative values in decreasing order.
*>       (Return 0 if true, 1/ULP if false.)
*>
*> (5)   | U - Upartial | / ( M ulp ) where Upartial is a partially
*>       computed U.
*>
*> (6)   | VT - VTpartial | / ( N ulp ) where VTpartial is a partially
*>       computed VT.
*>
*> (7)   | S - Spartial | / ( MNMIN ulp |S| ) where Spartial is the
*>       vector of singular values from the partial SVD
*>
*> Test for ZGESDD:
*>
*> (8)   | A - U diag(S) VT | / ( |A| max(M,N) ulp )
*>
*> (9)   | I - U'U | / ( M ulp )
*>
*> (10)  | I - VT VT' | / ( N ulp )
*>
*> (11)  S contains MNMIN nonnegative values in decreasing order.
*>       (Return 0 if true, 1/ULP if false.)
*>
*> (12)  | U - Upartial | / ( M ulp ) where Upartial is a partially
*>       computed U.
*>
*> (13)  | VT - VTpartial | / ( N ulp ) where VTpartial is a partially
*>       computed VT.
*>
*> (14)  | S - Spartial | / ( MNMIN ulp |S| ) where Spartial is the
*>       vector of singular values from the partial SVD
*>
*> Test for ZGESVDQ:
*>
*> (36)  | A - U diag(S) VT | / ( |A| max(M,N) ulp )
*>
*> (37)  | I - U'U | / ( M ulp )
*>
*> (38)  | I - VT VT' | / ( N ulp )
*>
*> (39)  S contains MNMIN nonnegative values in decreasing order.
*>       (Return 0 if true, 1/ULP if false.)
*>
*> Test for ZGESVJ:
*>
*> (15)  | A - U diag(S) VT | / ( |A| max(M,N) ulp )
*>
*> (16)  | I - U'U | / ( M ulp )
*>
*> (17)  | I - VT VT' | / ( N ulp )
*>
*> (18)  S contains MNMIN nonnegative values in decreasing order.
*>       (Return 0 if true, 1/ULP if false.)
*>
*> Test for ZGEJSV:
*>
*> (19)  | A - U diag(S) VT | / ( |A| max(M,N) ulp )
*>
*> (20)  | I - U'U | / ( M ulp )
*>
*> (21)  | I - VT VT' | / ( N ulp )
*>
*> (22)  S contains MNMIN nonnegative values in decreasing order.
*>        (Return 0 if true, 1/ULP if false.)
*>
*> Test for ZGESVDMR3( 'V', 'V', 'A' )/ZGESVDMR3( 'N', 'N', 'A' )
*>
*> (23)  | A - U diag(S) VT | / ( |A| max(M,N) ulp )
*>
*> (24)  | I - U'U | / ( M ulp )
*>
*> (25)  | I - VT VT' | / ( N ulp )
*>
*> (26)  S contains MNMIN nonnegative values in decreasing order.
*>       (Return 0 if true, 1/ULP if false.)
*>
*> (27)  | U - Upartial | / ( M ulp ) where Upartial is a partially
*>       computed U.
*>
*> (28)  | VT - VTpartial | / ( N ulp ) where VTpartial is a partially
*>       computed VT.
*>
*> (29)  | S - Spartial | / ( MNMIN ulp |S| ) where Spartial is the
*>       vector of singular values from the partial SVD
*>
*> Test for ZGESVDMR3( 'V', 'V', 'I' )
*>
*> (30)  | U' A VT''' - diag(S) | / ( |A| max(M,N) ulp )
*>
*> (31)  | I - U'U | / ( M ulp )
*>
*> (32)  | I - VT VT' | / ( N ulp )
*>
*> Test for ZGESVDMR3( 'V', 'V', 'V' )
*>
*> (33)   | U' A VT''' - diag(S) | / ( |A| max(M,N) ulp )
*>
*> (34)   | I - U'U | / ( M ulp )
*>
*> (35)   | I - VT VT' | / ( N ulp )
*>
*> The "sizes" are specified by the arrays MM(1:NSIZES) and
*> NN(1:NSIZES); the value of each element pair (MM(j),NN(j))
*> specifies one size.  The "types" are specified by a logical array
*> DOTYPE( 1:NTYPES ); if DOTYPE(j) is .TRUE., then matrix type "j"
*> will be generated.
*> Currently, the list of possible types is:
*>
*> (1)  The zero matrix.
*> (2)  The identity matrix.
*> (3)  A matrix of the form  U D V, where U and V are orthogonal and
*>      D has an arithmetic spectrum with condition number sqrt(1/ULP).
*> (4)  As (3), with a geometric spectrum, condition number sqrt(1/ULP).
*> (5)  As (3), with a log-uniform spectrum, condition number sqrt(1/ULP).
*> (6)  As (3), with an arithmetic spectrum, condition number 1/ULP.
*> (7)  As (3), with a geometric spectrum, condition number 1/ULP.
*> (8)  As (3), with one singular value equal to 1 and the rest equal to
*>      sqrt(ULP) (a clustered spectrum).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] NSIZES
*> \verbatim
*>          NSIZES is INTEGER
*>          The number of sizes of matrices to use.  If it is zero,
*>          ZDRVMR3 does nothing.  It must be at least zero.
*> \endverbatim
*>
*> \param[in] MM
*> \verbatim
*>          MM is INTEGER array, dimension (NSIZES)
*>          An array containing the matrix "heights" to be used.  For
*>          each j=1,...,NSIZES, if MM(j) is zero, then MM(j) and NN(j)
*>          will be ignored.  The MM(j) values must be at least zero.
*> \endverbatim
*>
*> \param[in] NN
*> \verbatim
*>          NN is INTEGER array, dimension (NSIZES)
*>          An array containing the matrix "widths" to be used.  For
*>          each j=1,...,NSIZES, if NN(j) is zero, then MM(j) and NN(j)
*>          will be ignored.  The NN(j) values must be at least zero.
*> \endverbatim
*>
*> \param[in] NTYPES
*> \verbatim
*>          NTYPES is INTEGER
*>          The number of elements in DOTYPE.   If it is zero, ZDRVMR3
*>          does nothing.  It must be at least zero.  If it is MAXTYP+1
*>          and NSIZES is 1, then an additional type, MAXTYP+1 is
*>          defined, which is to use whatever matrices are in A and B.
*>          This is only useful if DOTYPE(1:MAXTYP) is .FALSE. and
*>          DOTYPE(MAXTYP+1) is .TRUE. .
*> \endverbatim
*>
*> \param[in] DOTYPE
*> \verbatim
*>          DOTYPE is LOGICAL array, dimension (NTYPES)
*>          If DOTYPE(j) is .TRUE., then for each size (m,n), a matrix
*>          of type j will be generated.  If NTYPES is smaller than the
*>          maximum number of types defined (PARAMETER MAXTYP), then
*>          types NTYPES+1 through MAXTYP will not be generated.  If
*>          NTYPES is larger than MAXTYP, DOTYPE(MAXTYP+1) through
*>          DOTYPE(NTYPES) will be ignored.
*> \endverbatim
*>
*> \param[in,out] ISEED
*> \verbatim
*>          ISEED is INTEGER array, dimension (4)
*>          On entry ISEED specifies the seed of the random number
*>          generator. The array elements should be between 0 and 4095;
*>          if not they will be reduced mod 4096.  Also, ISEED(4) must
*>          be odd.  The random number generator uses a linear
*>          congruential sequence limited to small integers, and so
*>          should produce machine independent random numbers. The
*>          values of ISEED are changed on exit, and can be used in the
*>          next call to ZDRVMR3 to continue the same random number
*>          sequence.
*> \endverbatim
*>
*> \param[in] THRESH
*> \verbatim
*>          THRESH is DOUBLE PRECISION
*>          A test will count as "failed" if the "error", computed as
*>          described above, exceeds THRESH.  Note that the error
*>          is scaled to be O(1), so THRESH should be a reasonably
*>          small multiple of 1, e.g., 10 or 100.  In particular,
*>          it should not depend on the precision (single vs. double)
*>          or the size of the matrix.  It must be at least zero.
*> \endverbatim
*>
*> \param[out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,max(NN))
*>          Used to hold the matrix whose singular values are to be
*>          computed.  On exit, A contains the last matrix actually
*>          used.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of A.  It must be at
*>          least 1 and at least max( MM ).
*> \endverbatim
*>
*> \param[out] U
*> \verbatim
*>          U is COMPLEX*16 array, dimension (LDU,max(MM))
*>          Used to hold the computed matrix of right singular vectors.
*>          On exit, U contains the last such vectors actually computed.
*> \endverbatim
*>
*> \param[in] LDU
*> \verbatim
*>          LDU is INTEGER
*>          The leading dimension of U.  It must be at
*>          least 1 and at least max( MM ).
*> \endverbatim
*>
*> \param[out] VT
*> \verbatim
*>          VT is COMPLEX*16 array, dimension (LDVT,max(NN))
*>          Used to hold the computed matrix of left singular vectors.
*>          On exit, VT contains the last such vectors actually computed.
*> \endverbatim
*>
*> \param[in] LDVT
*> \verbatim
*>          LDVT is INTEGER
*>          The leading dimension of VT.  It must be at
*>          least 1 and at least max( NN ).
*> \endverbatim
*>
*> \param[out] ASAV
*> \verbatim
*>          ASAV is COMPLEX*16 array, dimension (LDA,max(NN))
*>          Used to hold a different copy of the matrix whose singular
*>          values are to be computed.  On exit, A contains the last
*>          matrix actually used.
*> \endverbatim
*>
*> \param[out] USAV
*> \verbatim
*>          USAV is COMPLEX*16 array, dimension (LDU,max(MM))
*>          Used to hold a different copy of the computed matrix of
*>          right singular vectors. On exit, USAV contains the last such
*>          vectors actually computed.
*> \endverbatim
*>
*> \param[out] VTSAV
*> \verbatim
*>          VTSAV is COMPLEX*16 array, dimension (LDVT,max(NN))
*>          Used to hold a different copy of the computed matrix of
*>          left singular vectors. On exit, VTSAV contains the last such
*>          vectors actually computed.
*> \endverbatim
*>
*> \param[out] S
*> \verbatim
*>          S is DOUBLE PRECISION array, dimension (max(min(MM,NN)))
*>          Contains the computed singular values.
*> \endverbatim
*>
*> \param[out] SSAV
*> \verbatim
*>          SSAV is DOUBLE PRECISION array, dimension (max(min(MM,NN)))
*>          Contains another copy of the computed singular values.
*> \endverbatim
*>
*> \param[out] E
*> \verbatim
*>          E is DOUBLE PRECISION array, dimension (max(min(MM,NN)))
*>          Workspace for ZGESVD.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (LWORK)
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The number of entries in WORK.  This must be at least
*>          MAX(3*MIN(M,N)+MAX(M,N)**2,5*MIN(M,N),3*MAX(M,N)) for all
*>          pairs  (M,N)=(MM(j),NN(j))
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array,
*>                      dimension ( 5*max(max(MM,NN)) )
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension at least 8*min(M,N)
*> \endverbatim
*>
*> \param[in] NOUNIT
*> \verbatim
*>          NOUNIT is INTEGER
*>          The FORTRAN unit number for printing out error messages
*>          (e.g., if a routine returns IINFO not equal to 0.)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          If 0, then everything ran OK.
*>           -1: NSIZES < 0
*>           -2: Some MM(j) < 0
*>           -3: Some NN(j) < 0
*>           -4: NTYPES < 0
*>           -7: THRESH < 0
*>          -10: LDA < 1 or LDA < MMAX, where MMAX is max( MM(j) ).
*>          -12: LDU < 1 or LDU < MMAX.
*>          -14: LDVT < 1 or LDVT < NMAX, where NMAX is max( NN(j) ).
*>          -21: LWORK too small.
*>          If  ZLATMS, or ZGESVD returns an error code, the
*>              absolute value of it is returned.
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
*> \ingroup complex16_eig
*
*  =====================================================================
      SUBROUTINE ZDRVMR3( NSIZES, MM, NN, NTYPES, DOTYPE, ISEED, THRESH,
     $                   A, LDA, U, LDU, VT, LDVT, ASAV, USAV, VTSAV, S,
     $                   SSAV, E, WORK, LWORK, RWORK, IWORK, NOUNIT,
     $                   INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDU, LDVT, LWORK, NOUNIT, NSIZES,
     $                   NTYPES
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            ISEED( 4 ), IWORK( * ), MM( * ), NN( * )
      DOUBLE PRECISION   E( * ), RWORK( * ), S( * ), SSAV( * )
      COMPLEX*16         A( LDA, * ), ASAV( LDA, * ), U( LDU, * ),
     $                   USAV( LDU, * ), VT( LDVT, * ),
     $                   VTSAV( LDVT, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, HALF
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   HALF = 0.5D0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
      INTEGER            MAXTYP
      PARAMETER          ( MAXTYP = 8 )
*     ..
*     .. Local Scalars ..
      LOGICAL            BADMM, BADNN
      CHARACTER          JOBQ, JOBU, JOBVT, RANGE
      INTEGER            I, IINFO, IJQ, IJU, IJVT, IL, IU, ITEMP,
     $                   IWSPC, IWTMP, J, JSIZE, JTYPE, LSWORK, M,
     $                   MINWRK, MMAX, MNMAX, MNMIN, MTYPES, N,
     $                   NERRS, NFAIL, NMAX, NS, NSI, NSV, NTEST,
     $                   NTESTF, NTESTT, LRWORK
      DOUBLE PRECISION   ANORM, DIF, DIV, OVFL, RTUNFL, ULP, ULPINV,
     $                   UNFL, VL, VU
*     ..
*     .. Local Scalars for ZGESVDQ ..
      INTEGER            LIWORK, NUMRANK
*     ..
*     .. Local Arrays ..
      CHARACTER          CJOB( 4 ), CJOBR( 3 ), CJOBV( 2 )
      INTEGER            IOLDSD( 4 ), ISEED2( 4 )
      DOUBLE PRECISION   RESULT( 13 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLARND, ZLANGE
      EXTERNAL           DLAMCH, DLARND, ZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALASVM, XERBLA, ZBDT01, ZBDT05, ZGESDD,
     $                   ZGESVD, ZGESVDQ, ZGESVJ, ZGEJSV, ZGESVDMR3,
     $                   ZLACPY, ZLASET, ZLATMS, ZUNT01, ZUNT03
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN
*     ..
*     .. Scalars in Common ..
      CHARACTER*32       SRNAMT
*     ..
*     .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Data statements ..
      DATA               CJOB / 'N', 'O', 'S', 'A' /
      DATA               CJOBR / 'A', 'V', 'I' /
      DATA               CJOBV / 'N', 'V' /
*     ..
*     .. Executable Statements ..
*
*     Check for errors
*
      INFO = 0
*
*     Important constants
*
      NERRS = 0
      NTESTT = 0
      NTESTF = 0
      BADMM = .FALSE.
      BADNN = .FALSE.
      MMAX = 1
      NMAX = 1
      MNMAX = 1
      MINWRK = 1
      DO 10 J = 1, NSIZES
         MMAX = MAX( MMAX, MM( J ) )
         IF( MM( J ).LT.0 )
     $      BADMM = .TRUE.
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ).LT.0 )
     $      BADNN = .TRUE.
         MNMAX = MAX( MNMAX, MIN( MM( J ), NN( J ) ) )
         MINWRK = MAX( MINWRK, MAX( 3*MIN( MM( J ),
     $            NN( J ) )+MAX( MM( J ), NN( J ) )**2, 5*MIN( MM( J ),
     $            NN( J ) ), 3*MAX( MM( J ), NN( J ) ) ) )
   10 CONTINUE
*
*     Check for errors
*
      IF( NSIZES.LT.0 ) THEN
         INFO = -1
      ELSE IF( BADMM ) THEN
         INFO = -2
      ELSE IF( BADNN ) THEN
         INFO = -3
      ELSE IF( NTYPES.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, MMAX ) ) THEN
         INFO = -10
      ELSE IF( LDU.LT.MAX( 1, MMAX ) ) THEN
         INFO = -12
      ELSE IF( LDVT.LT.MAX( 1, NMAX ) ) THEN
         INFO = -14
      ELSE IF( MINWRK.GT.LWORK ) THEN
         INFO = -21
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZDRVMR3', -INFO )
         RETURN
      END IF
*
*     Quick return if nothing to do
*
      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 )
     $   RETURN
*
*     More Important constants
*
      UNFL = DLAMCH( 'S' )
      OVFL = ONE / UNFL
      ULP = DLAMCH( 'E' )
      ULPINV = ONE / ULP
      RTUNFL = SQRT( UNFL )
*
*     Loop over sizes, types
*
      NERRS = 0
*
      DO 230 JSIZE = 1, NSIZES
         M = MM( JSIZE )
         N = NN( JSIZE )
         MNMIN = MIN( M, N )
*
         IF( NSIZES.NE.1 ) THEN
            MTYPES = MIN( MAXTYP, NTYPES )
         ELSE
            MTYPES = MIN( MAXTYP+1, NTYPES )
         END IF
*
         DO 220 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) )
     $         GO TO 220
            NTEST = 0
*
            DO 20 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   20       CONTINUE
*
*           Compute "A"
*
            IF( MTYPES.GT.MAXTYP )
     $         GO TO 50
*
            IINFO = 0
            ANORM = ONE
            IF( JTYPE.EQ.1 ) THEN
               CALL ZLASET( 'Full', M, N, CZERO, CZERO, A, LDA )
            ELSE IF( JTYPE.EQ.2 ) THEN
               CALL ZLASET( 'Full', M, N, CZERO, CONE, A, LDA )
            ELSE IF( JTYPE.EQ.3 ) THEN
*              arithmetic spectrum, condition sqrt(1/ulp)
               CALL ZLATMS( M, N, 'U', ISEED, 'N', S, 4,
     $                      SQRT( ULPINV ), ONE, M-1, N-1, 'N', A, LDA,
     $                      WORK, IINFO )
            ELSE IF( JTYPE.EQ.4 ) THEN
*              geometric spectrum, condition sqrt(1/ulp)
               CALL ZLATMS( M, N, 'U', ISEED, 'N', S, 3,
     $                      SQRT( ULPINV ), ONE, M-1, N-1, 'N', A, LDA,
     $                      WORK, IINFO )
            ELSE IF( JTYPE.EQ.5 ) THEN
*              log-uniform spectrum, condition sqrt(1/ulp)
               CALL ZLATMS( M, N, 'U', ISEED, 'N', S, 5,
     $                      SQRT( ULPINV ), ONE, M-1, N-1, 'N', A, LDA,
     $                      WORK, IINFO )
            ELSE IF( JTYPE.EQ.6 ) THEN
*              arithmetic spectrum, condition 1/ulp
               CALL ZLATMS( M, N, 'U', ISEED, 'N', S, 4,
     $                      ULPINV, ONE, M-1, N-1, 'N', A, LDA,
     $                      WORK, IINFO )
            ELSE IF( JTYPE.EQ.7 ) THEN
*              geometric spectrum, condition 1/ulp
               CALL ZLATMS( M, N, 'U', ISEED, 'N', S, 3,
     $                      ULPINV, ONE, M-1, N-1, 'N', A, LDA,
     $                      WORK, IINFO )
            ELSE
*              clustered spectrum: one singular value 1, the rest sqrt(ulp)
               CALL ZLATMS( M, N, 'U', ISEED, 'N', S, 1,
     $                      SQRT( ULPINV ), ONE, M-1, N-1, 'N', A, LDA,
     $                      WORK, IINFO )
            END IF
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9996 )'Generator', IINFO, M, N,
     $            JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
*
   50       CONTINUE
            CALL ZLACPY( 'F', M, N, A, LDA, ASAV, LDA )
*
*           Do for minimal and adequate (for blocking) workspace
*
            DO 210 IWSPC = 1, 4
*
*              Test for ZGESVD
*
               IWTMP = 2*MIN( M, N )+MAX( M, N )
               LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
               LSWORK = MIN( LSWORK, LWORK )
               LSWORK = MAX( LSWORK, 1 )
               IF( IWSPC.EQ.4 )
     $            LSWORK = LWORK
*
               DO 60 J = 1, 13
                  RESULT( J ) = -ONE
   60          CONTINUE
*
*              Test ZGESVDMR3
*
*              Factorize A
*
               CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               SRNAMT = 'ZGESVDMR3'
               CALL ZGESVDMR3( 'V', 'V', 'A', M, N, A, LDA,
     $                       VL, VU, IL, IU, NS, SSAV, USAV, LDU,
     $                       VTSAV, LDVT, WORK, LWORK, RWORK,
     $                       IWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9995 )'GESVDMR3', IINFO, M, N,
     $               JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF
*
*              Do tests -21---18
*
               RESULT( 1 ) = ZERO
               RESULT( 2 ) = ZERO
               RESULT( 3 ) = ZERO
               CALL ZBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E,
     $                      VTSAV, LDVT, WORK, RWORK, RESULT( 1 ) )
               IF( M.NE.0 .AND. N.NE.0 ) THEN
                  CALL ZUNT01( 'Columns', M, MNMIN, USAV, LDU, WORK,
     $                         LWORK, RWORK, RESULT( 2 ) )
                  CALL ZUNT01( 'Rows', MNMIN, N, VTSAV, LDVT, WORK,
     $                         LWORK, RWORK, RESULT( 3 ) )
               END IF
               RESULT( 4 ) = ZERO
               DO 140 I = 1, MNMIN - 1
                  IF( SSAV( I ).LT.SSAV( I+1 ) )
     $               RESULT( 4 ) = ULPINV
                  IF( SSAV( I ).LT.ZERO )
     $               RESULT( 4 ) = ULPINV
  140          CONTINUE
               IF( MNMIN.GE.1 ) THEN
                  IF( SSAV( MNMIN ).LT.ZERO )
     $               RESULT( 4 ) = ULPINV
               END IF
*
*              Do partial SVDs, comparing to SSAV, USAV, and VTSAV
*
               RESULT( 5 ) = ZERO
               RESULT( 6 ) = ZERO
               RESULT( 7 ) = ZERO
               DO 170 IJU = 0, 1
                  DO 160 IJVT = 0, 1
                     IF( ( IJU.EQ.0 .AND. IJVT.EQ.0 ) .OR.
     $                   ( IJU.EQ.1 .AND. IJVT.EQ.1 ) ) GO TO 160
                     JOBU = CJOBV( IJU+1 )
                     JOBVT = CJOBV( IJVT+1 )
                     RANGE = CJOBR( 1 )
                     CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
                     SRNAMT = 'ZGESVDMR3'
                     CALL ZGESVDMR3( JOBU, JOBVT, 'A', M, N, A, LDA,
     $                            VL, VU, IL, IU, NS, S, U, LDU,
     $                            VT, LDVT, WORK, LWORK, RWORK,
     $                            IWORK, IINFO )
*
*                    Compare U
*
                     DIF = ZERO
                     IF( M.GT.0 .AND. N.GT.0 ) THEN
                        IF( IJU.EQ.1 ) THEN
                           CALL ZUNT03( 'C', M, MNMIN, M, MNMIN, USAV,
     $                                  LDU, U, LDU, WORK, LWORK, RWORK,
     $                                  DIF, IINFO )
                        END IF
                     END IF
                     RESULT( 5 ) = MAX( RESULT( 5 ), DIF )
*
*                    Compare VT
*
                     DIF = ZERO
                     IF( M.GT.0 .AND. N.GT.0 ) THEN
                        IF( IJVT.EQ.1 ) THEN
                           CALL ZUNT03( 'R', N, MNMIN, N, MNMIN, VTSAV,
     $                                  LDVT, VT, LDVT, WORK, LWORK,
     $                                  RWORK, DIF, IINFO )
                        END IF
                     END IF
                     RESULT( 6 ) = MAX( RESULT( 6 ), DIF )
*
*                    Compare S
*
                     DIF = ZERO
                     DIV = MAX( DBLE( MNMIN )*ULP*S( 1 ),
     $                     DLAMCH( 'Safe minimum' ) )
                     DO 150 I = 1, MNMIN - 1
                        IF( SSAV( I ).LT.SSAV( I+1 ) )
     $                     DIF = ULPINV
                        IF( SSAV( I ).LT.ZERO )
     $                     DIF = ULPINV
                        DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
  150                CONTINUE
                     RESULT( 7 ) = MAX( RESULT( 7 ), DIF )
  160             CONTINUE
  170          CONTINUE
*
*              Do tests -14---12
*
               DO 180 I = 1, 4
                  ISEED2( I ) = ISEED( I )
  180          CONTINUE
               IF( MNMIN.LE.1 ) THEN
                  IL = 1
                  IU = MAX( 1, MNMIN )
               ELSE
                  IL = 1 + INT( MNMIN*DLARND( 1, ISEED2 ) )
                  IU = 1 + INT( MNMIN*DLARND( 1, ISEED2 ) )
                  IF( IU.LT.IL ) THEN
                     ITEMP = IU
                     IU = IL
                     IL = ITEMP
                  END IF
               END IF
               CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               SRNAMT = 'ZGESVDMR3'
               CALL ZGESVDMR3( 'V', 'V', 'I', M, N, A, LDA,
     $                       VL, VU, IL, IU, NSI, S, U, LDU,
     $                       VT, LDVT, WORK, LWORK, RWORK,
     $                       IWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9995 )'GESVDMR3', IINFO, M, N,
     $               JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF
*
               RESULT( 8 ) = ZERO
               RESULT( 9 ) = ZERO
               RESULT( 10 ) = ZERO
               CALL ZBDT05( M, N, ASAV, LDA, S, NSI, U, LDU,
     $                      VT, LDVT, WORK, RESULT( 8 ) )
*
*              ZBDT05 normalizes by max |A(i,j)|; rescale to the 1-norm
*              used by ZBDT01 in test (1), so that tests (1), (8) and
*              (11) measure the same quantity.
*
               ANORM = ZLANGE( '1', M, N, ASAV, LDA, RWORK )
               IF( ANORM.GT.ZERO )
     $            RESULT( 8 ) = RESULT( 8 )*
     $                 ( ZLANGE( 'M', M, N, ASAV, LDA, RWORK ) / ANORM )
               IF( M.NE.0 .AND. N.NE.0 ) THEN
                  CALL ZUNT01( 'Columns', M, NSI, U, LDU, WORK,
     $                         LWORK, RWORK, RESULT( 9 ) )
                  CALL ZUNT01( 'Rows', NSI, N, VT, LDVT, WORK,
     $                         LWORK, RWORK, RESULT( 10 ) )
               END IF
*
*              Do tests -11---9
*
               IF( MNMIN.GT.0 .AND. NSI.GT.1 ) THEN
                  IF( IL.NE.1 ) THEN
                     VU = SSAV( IL ) +
     $                    MAX( HALF*ABS( SSAV( IL )-SSAV( IL-1 ) ),
     $                    ULP*ANORM, TWO*RTUNFL )
                  ELSE
                     VU = SSAV( 1 ) +
     $                    MAX( HALF*ABS( SSAV( NS )-SSAV( 1 ) ),
     $                    ULP*ANORM, TWO*RTUNFL )
                  END IF
                  IF( IU.NE.NS ) THEN
                     VL = SSAV( IU ) - MAX( ULP*ANORM, TWO*RTUNFL,
     $                    HALF*ABS( SSAV( IU+1 )-SSAV( IU ) ) )
                  ELSE
                     VL = SSAV( NS ) - MAX( ULP*ANORM, TWO*RTUNFL,
     $                    HALF*ABS( SSAV( NS )-SSAV( 1 ) ) )
                  END IF
                  VL = MAX( VL,ZERO )
                  VU = MAX( VU,ZERO )
                  IF( VL.GE.VU ) VU = MAX( VU*2, VU+VL+HALF )
               ELSE
                  VL = ZERO
                  VU = ONE
               END IF
               CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               SRNAMT = 'ZGESVDMR3'
               CALL ZGESVDMR3( 'V', 'V', 'V', M, N, A, LDA,
     $                       VL, VU, IL, IU, NSV, S, U, LDU,
     $                       VT, LDVT, WORK, LWORK, RWORK,
     $                       IWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9995 )'GESVDMR3', IINFO, M, N,
     $               JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF
*
               RESULT( 11 ) = ZERO
               RESULT( 12 ) = ZERO
               RESULT( 13 ) = ZERO
               CALL ZBDT05( M, N, ASAV, LDA, S, NSV, U, LDU,
     $                      VT, LDVT, WORK, RESULT( 11 ) )
*
*              ZBDT05 normalizes by max |A(i,j)|; rescale to the 1-norm
*              used by ZBDT01 in test (1), so that tests (1), (8) and
*              (11) measure the same quantity.
*
               ANORM = ZLANGE( '1', M, N, ASAV, LDA, RWORK )
               IF( ANORM.GT.ZERO )
     $            RESULT( 11 ) = RESULT( 11 )*
     $                 ( ZLANGE( 'M', M, N, ASAV, LDA, RWORK ) / ANORM )
               IF( M.NE.0 .AND. N.NE.0 ) THEN
                  CALL ZUNT01( 'Columns', M, NSV, U, LDU, WORK,
     $                         LWORK, RWORK, RESULT( 12 ) )
                  CALL ZUNT01( 'Rows', NSV, N, VT, LDVT, WORK,
     $                         LWORK, RWORK, RESULT( 13 ) )
               END IF
*
*              End of Loop -- Check for RESULT(j) > THRESH
*
               NTEST = 0
               NFAIL = 0
               DO 190 J = 1, 13
                  IF( RESULT( J ).GE.ZERO )
     $               NTEST = NTEST + 1
                  IF( RESULT( J ).GE.THRESH )
     $               NFAIL = NFAIL + 1
  190          CONTINUE
*
               IF( NFAIL.GT.0 )
     $            NTESTF = NTESTF + 1
               IF( NTESTF.EQ.1 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )
                  WRITE( NOUNIT, FMT = 9998 )THRESH
                  NTESTF = 2
               END IF
*
               DO 200 J = 1, 13
                  IF( RESULT( J ).GE.THRESH ) THEN
                     WRITE( NOUNIT, FMT = 9997 )M, N, JTYPE, IWSPC,
     $                  IOLDSD, J, RESULT( J )
                  END IF
  200          CONTINUE
*
               NERRS = NERRS + NFAIL
               NTESTT = NTESTT + NTEST
*
  210       CONTINUE
*
  220    CONTINUE
  230 CONTINUE
*
*     Summary
*
      CALL ALASVM( 'ZMR', NOUNIT, NERRS, NTESTT, 0 )
*
 9999 FORMAT( ' SVD -- ', A, ' Singular Value Decomposition Driver ',
     $      / ' ZGESVDMR3 (MRRR dense driver)',
     $      / ' Matrix types (see ZDRVMR3 for details):',
     $      / ' 1 = Zero matrix', / ' 2 = Identity matrix',
     $      / ' 3 = Arithmetic singular values, condition sqrt(1/ulp)',
     $      / ' 4 = Geometric singular values, condition sqrt(1/ulp)',
     $      / ' 5 = Log-uniform singular values, condition sqrt(1/ulp)',
     $      / ' 6 = Arithmetic singular values, condition 1/ulp',
     $      / ' 7 = Geometric singular values, condition 1/ulp',
     $      / ' 8 = Clustered: one singular value 1, the rest',
     $      'sqrt(ulp)',
     $      / ' Tests performed: ( A is dense, U and VT',
     $      'unitary/orthogonal )',
     $      / '  1 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ),',
     $      '(V,V,A)',
     $      / '  2 = | I - U**H U | / ( M ulp ) ',
     $      / '  3 = | I - VT VT**H | / ( N ulp ) ',
     $      / '  4 = 0 if S contains min(M,N) nonnegative values in',
     $      ' decreasing order, else 1/ulp',
     $      / '  5 = | U - Upartial | / ( M ulp )',
     $      / '  6 = | VT - VTpartial | / ( N ulp )',
     $      / '  7 = | S - Spartial | / ( min(M,N) ulp |S| )',
     $      / '  8 = | U**H A VT**H - diag(S) | / ( |A| max(M,N)',
     $      'ulp ), (V,V,I)',
     $      / '  9 = | I - U**H U | / ( M ulp ) ',
     $      / ' 10 = | I - VT VT**H | / ( N ulp ) ',
     $      / ' 11 = | U**H A VT**H - diag(S) | / ( |A| max(M,N)',
     $      'ulp ), (V,V,V)',
     $      / ' 12 = | I - U**H U | / ( M ulp ) ',
     $      / ' 13 = | I - VT VT**H | / ( N ulp ) ',
     $      / / )
 9998 FORMAT( ' Tests performed with Test Threshold = ', F8.2 )
 9997 FORMAT( ' M=', I5, ', N=', I5, ', type ', I1, ', IWS=', I1,
     $      ', seed=', 4( I4, ',' ), ' test(', I2, ')=', G11.4 )
 9996 FORMAT( ' ZDRVMR3: ', A, ' returned INFO=', I6, '.', / 9X, 'M=',
     $      I6, ', N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ),
     $      I5, ')' )
 9995 FORMAT( ' ZDRVMR3: ', A, ' returned INFO=', I6, '.', / 9X, 'M=',
     $      I6, ', N=', I6, ', JTYPE=', I6, ', LSWORK=', I6, / 9X,
     $      'ISEED=(', 3( I5, ',' ), I5, ')' )
*
      RETURN
*
*     End of ZDRVMR3
*
      END
