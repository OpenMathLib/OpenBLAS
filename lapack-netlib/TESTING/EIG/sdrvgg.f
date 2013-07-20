*> \brief \b SDRVGG
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE SDRVGG( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
*                          THRSHN, NOUNIT, A, LDA, B, S, T, S2, T2, Q,
*                          LDQ, Z, ALPHR1, ALPHI1, BETA1, ALPHR2, ALPHI2,
*                          BETA2, VL, VR, WORK, LWORK, RESULT, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, LDQ, LWORK, NOUNIT, NSIZES, NTYPES
*       REAL               THRESH, THRSHN
*       ..
*       .. Array Arguments ..
*       LOGICAL            DOTYPE( * )
*       INTEGER            ISEED( 4 ), NN( * )
*       REAL               A( LDA, * ), ALPHI1( * ), ALPHI2( * ),
*      $                   ALPHR1( * ), ALPHR2( * ), B( LDA, * ),
*      $                   BETA1( * ), BETA2( * ), Q( LDQ, * ),
*      $                   RESULT( * ), S( LDA, * ), S2( LDA, * ),
*      $                   T( LDA, * ), T2( LDA, * ), VL( LDQ, * ),
*      $                   VR( LDQ, * ), WORK( * ), Z( LDQ, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SDRVGG  checks the nonsymmetric generalized eigenvalue driver
*> routines.
*>                               T          T        T
*> SGEGS factors A and B as Q S Z  and Q T Z , where   means
*> transpose, T is upper triangular, S is in generalized Schur form
*> (block upper triangular, with 1x1 and 2x2 blocks on the diagonal,
*> the 2x2 blocks corresponding to complex conjugate pairs of
*> generalized eigenvalues), and Q and Z are orthogonal.  It also
*> computes the generalized eigenvalues (alpha(1),beta(1)), ...,
*> (alpha(n),beta(n)), where alpha(j)=S(j,j) and beta(j)=P(j,j) --
*> thus, w(j) = alpha(j)/beta(j) is a root of the generalized
*> eigenvalue problem
*>
*>     det( A - w(j) B ) = 0
*>
*> and m(j) = beta(j)/alpha(j) is a root of the essentially equivalent
*> problem
*>
*>     det( m(j) A - B ) = 0
*>
*> SGEGV computes the generalized eigenvalues (alpha(1),beta(1)), ...,
*> (alpha(n),beta(n)), the matrix L whose columns contain the
*> generalized left eigenvectors l, and the matrix R whose columns
*> contain the generalized right eigenvectors r for the pair (A,B).
*>
*> When SDRVGG is called, a number of matrix "sizes" ("n's") and a
*> number of matrix "types" are specified.  For each size ("n")
*> and each type of matrix, one matrix will be generated and used
*> to test the nonsymmetric eigenroutines.  For each matrix, 7
*> tests will be performed and compared with the threshhold THRESH:
*>
*> Results from SGEGS:
*>
*>                  T
*> (1)   | A - Q S Z  | / ( |A| n ulp )
*>
*>                  T
*> (2)   | B - Q T Z  | / ( |B| n ulp )
*>
*>               T
*> (3)   | I - QQ  | / ( n ulp )
*>
*>               T
*> (4)   | I - ZZ  | / ( n ulp )
*>
*> (5)   maximum over j of D(j)  where:
*>
*> if alpha(j) is real:
*>                     |alpha(j) - S(j,j)|        |beta(j) - T(j,j)|
*>           D(j) = ------------------------ + -----------------------
*>                  max(|alpha(j)|,|S(j,j)|)   max(|beta(j)|,|T(j,j)|)
*>
*> if alpha(j) is complex:
*>                                 | det( s S - w T ) |
*>           D(j) = ---------------------------------------------------
*>                  ulp max( s norm(S), |w| norm(T) )*norm( s S - w T )
*>
*>           and S and T are here the 2 x 2 diagonal blocks of S and T
*>           corresponding to the j-th eigenvalue.
*>
*> Results from SGEGV:
*>
*> (6)   max over all left eigenvalue/-vector pairs (beta/alpha,l) of
*>
*>    | l**H * (beta A - alpha B) | / ( ulp max( |beta A|, |alpha B| ) )
*>
*>       where l**H is the conjugate tranpose of l.
*>
*> (7)   max over all right eigenvalue/-vector pairs (beta/alpha,r) of
*>
*>       | (beta A - alpha B) r | / ( ulp max( |beta A|, |alpha B| ) )
*>
*> Test Matrices
*> ---- --------
*>
*> The sizes of the test matrices are specified by an array
*> NN(1:NSIZES); the value of each element NN(j) specifies one size.
*> The "types" are specified by a logical array DOTYPE( 1:NTYPES ); if
*> DOTYPE(j) is .TRUE., then matrix type "j" will be generated.
*> Currently, the list of possible types is:
*>
*> (1)  ( 0, 0 )         (a pair of zero matrices)
*>
*> (2)  ( I, 0 )         (an identity and a zero matrix)
*>
*> (3)  ( 0, I )         (an identity and a zero matrix)
*>
*> (4)  ( I, I )         (a pair of identity matrices)
*>
*>         t   t
*> (5)  ( J , J  )       (a pair of transposed Jordan blocks)
*>
*>                                     t                ( I   0  )
*> (6)  ( X, Y )         where  X = ( J   0  )  and Y = (      t )
*>                                  ( 0   I  )          ( 0   J  )
*>                       and I is a k x k identity and J a (k+1)x(k+1)
*>                       Jordan block; k=(N-1)/2
*>
*> (7)  ( D, I )         where D is diag( 0, 1,..., N-1 ) (a diagonal
*>                       matrix with those diagonal entries.)
*> (8)  ( I, D )
*>
*> (9)  ( big*D, small*I ) where "big" is near overflow and small=1/big
*>
*> (10) ( small*D, big*I )
*>
*> (11) ( big*I, small*D )
*>
*> (12) ( small*I, big*D )
*>
*> (13) ( big*D, big*I )
*>
*> (14) ( small*D, small*I )
*>
*> (15) ( D1, D2 )        where D1 is diag( 0, 0, 1, ..., N-3, 0 ) and
*>                        D2 is diag( 0, N-3, N-4,..., 1, 0, 0 )
*>           t   t
*> (16) Q ( J , J ) Z     where Q and Z are random orthogonal matrices.
*>
*> (17) Q ( T1, T2 ) Z    where T1 and T2 are upper triangular matrices
*>                        with random O(1) entries above the diagonal
*>                        and diagonal entries diag(T1) =
*>                        ( 0, 0, 1, ..., N-3, 0 ) and diag(T2) =
*>                        ( 0, N-3, N-4,..., 1, 0, 0 )
*>
*> (18) Q ( T1, T2 ) Z    diag(T1) = ( 0, 0, 1, 1, s, ..., s, 0 )
*>                        diag(T2) = ( 0, 1, 0, 1,..., 1, 0 )
*>                        s = machine precision.
*>
*> (19) Q ( T1, T2 ) Z    diag(T1)=( 0,0,1,1, 1-d, ..., 1-(N-5)*d=s, 0 )
*>                        diag(T2) = ( 0, 1, 0, 1, ..., 1, 0 )
*>
*>                                                        N-5
*> (20) Q ( T1, T2 ) Z    diag(T1)=( 0, 0, 1, 1, a, ..., a   =s, 0 )
*>                        diag(T2) = ( 0, 1, 0, 1, ..., 1, 0, 0 )
*>
*> (21) Q ( T1, T2 ) Z    diag(T1)=( 0, 0, 1, r1, r2, ..., r(N-4), 0 )
*>                        diag(T2) = ( 0, 1, 0, 1, ..., 1, 0, 0 )
*>                        where r1,..., r(N-4) are random.
*>
*> (22) Q ( big*T1, small*T2 ) Z    diag(T1) = ( 0, 0, 1, ..., N-3, 0 )
*>                                  diag(T2) = ( 0, 1, ..., 1, 0, 0 )
*>
*> (23) Q ( small*T1, big*T2 ) Z    diag(T1) = ( 0, 0, 1, ..., N-3, 0 )
*>                                  diag(T2) = ( 0, 1, ..., 1, 0, 0 )
*>
*> (24) Q ( small*T1, small*T2 ) Z  diag(T1) = ( 0, 0, 1, ..., N-3, 0 )
*>                                  diag(T2) = ( 0, 1, ..., 1, 0, 0 )
*>
*> (25) Q ( big*T1, big*T2 ) Z      diag(T1) = ( 0, 0, 1, ..., N-3, 0 )
*>                                  diag(T2) = ( 0, 1, ..., 1, 0, 0 )
*>
*> (26) Q ( T1, T2 ) Z     where T1 and T2 are random upper-triangular
*>                         matrices.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] NSIZES
*> \verbatim
*>          NSIZES is INTEGER
*>          The number of sizes of matrices to use.  If it is zero,
*>          SDRVGG does nothing.  It must be at least zero.
*> \endverbatim
*>
*> \param[in] NN
*> \verbatim
*>          NN is INTEGER array, dimension (NSIZES)
*>          An array containing the sizes to be used for the matrices.
*>          Zero values will be skipped.  The values must be at least
*>          zero.
*> \endverbatim
*>
*> \param[in] NTYPES
*> \verbatim
*>          NTYPES is INTEGER
*>          The number of elements in DOTYPE.   If it is zero, SDRVGG
*>          does nothing.  It must be at least zero.  If it is MAXTYP+1
*>          and NSIZES is 1, then an additional type, MAXTYP+1 is
*>          defined, which is to use whatever matrix is in A.  This
*>          is only useful if DOTYPE(1:MAXTYP) is .FALSE. and
*>          DOTYPE(MAXTYP+1) is .TRUE. .
*> \endverbatim
*>
*> \param[in] DOTYPE
*> \verbatim
*>          DOTYPE is LOGICAL array, dimension (NTYPES)
*>          If DOTYPE(j) is .TRUE., then for each size in NN a
*>          matrix of that size and of type j will be generated.
*>          If NTYPES is smaller than the maximum number of types
*>          defined (PARAMETER MAXTYP), then types NTYPES+1 through
*>          MAXTYP will not be generated.  If NTYPES is larger
*>          than MAXTYP, DOTYPE(MAXTYP+1) through DOTYPE(NTYPES)
*>          will be ignored.
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
*>          next call to SDRVGG to continue the same random number
*>          sequence.
*> \endverbatim
*>
*> \param[in] THRESH
*> \verbatim
*>          THRESH is REAL
*>          A test will count as "failed" if the "error", computed as
*>          described above, exceeds THRESH.  Note that the error is
*>          scaled to be O(1), so THRESH should be a reasonably small
*>          multiple of 1, e.g., 10 or 100.  In particular, it should
*>          not depend on the precision (single vs. double) or the size
*>          of the matrix.  It must be at least zero.
*> \endverbatim
*>
*> \param[in] THRSHN
*> \verbatim
*>          THRSHN is REAL
*>          Threshhold for reporting eigenvector normalization error.
*>          If the normalization of any eigenvector differs from 1 by
*>          more than THRSHN*ulp, then a special error message will be
*>          printed.  (This is handled separately from the other tests,
*>          since only a compiler or programming error should cause an
*>          error message, at least if THRSHN is at least 5--10.)
*> \endverbatim
*>
*> \param[in] NOUNIT
*> \verbatim
*>          NOUNIT is INTEGER
*>          The FORTRAN unit number for printing out error messages
*>          (e.g., if a routine returns IINFO not equal to 0.)
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension
*>                            (LDA, max(NN))
*>          Used to hold the original A matrix.  Used as input only
*>          if NTYPES=MAXTYP+1, DOTYPE(1:MAXTYP)=.FALSE., and
*>          DOTYPE(MAXTYP+1)=.TRUE.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of A, B, S, T, S2, and T2.
*>          It must be at least 1 and at least max( NN ).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is REAL array, dimension
*>                            (LDA, max(NN))
*>          Used to hold the original B matrix.  Used as input only
*>          if NTYPES=MAXTYP+1, DOTYPE(1:MAXTYP)=.FALSE., and
*>          DOTYPE(MAXTYP+1)=.TRUE.
*> \endverbatim
*>
*> \param[out] S
*> \verbatim
*>          S is REAL array, dimension (LDA, max(NN))
*>          The Schur form matrix computed from A by SGEGS.  On exit, S
*>          contains the Schur form matrix corresponding to the matrix
*>          in A.
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          T is REAL array, dimension (LDA, max(NN))
*>          The upper triangular matrix computed from B by SGEGS.
*> \endverbatim
*>
*> \param[out] S2
*> \verbatim
*>          S2 is REAL array, dimension (LDA, max(NN))
*>          The matrix computed from A by SGEGV.  This will be the
*>          Schur form of some matrix related to A, but will not, in
*>          general, be the same as S.
*> \endverbatim
*>
*> \param[out] T2
*> \verbatim
*>          T2 is REAL array, dimension (LDA, max(NN))
*>          The matrix computed from B by SGEGV.  This will be the
*>          Schur form of some matrix related to B, but will not, in
*>          general, be the same as T.
*> \endverbatim
*>
*> \param[out] Q
*> \verbatim
*>          Q is REAL array, dimension (LDQ, max(NN))
*>          The (left) orthogonal matrix computed by SGEGS.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>          The leading dimension of Q, Z, VL, and VR.  It must
*>          be at least 1 and at least max( NN ).
*> \endverbatim
*>
*> \param[out] Z
*> \verbatim
*>          Z is REAL array of
*>                             dimension( LDQ, max(NN) )
*>          The (right) orthogonal matrix computed by SGEGS.
*> \endverbatim
*>
*> \param[out] ALPHR1
*> \verbatim
*>          ALPHR1 is REAL array, dimension (max(NN))
*> \endverbatim
*>
*> \param[out] ALPHI1
*> \verbatim
*>          ALPHI1 is REAL array, dimension (max(NN))
*> \endverbatim
*>
*> \param[out] BETA1
*> \verbatim
*>          BETA1 is REAL array, dimension (max(NN))
*>
*>          The generalized eigenvalues of (A,B) computed by SGEGS.
*>          ( ALPHR1(k)+ALPHI1(k)*i ) / BETA1(k) is the k-th
*>          generalized eigenvalue of the matrices in A and B.
*> \endverbatim
*>
*> \param[out] ALPHR2
*> \verbatim
*>          ALPHR2 is REAL array, dimension (max(NN))
*> \endverbatim
*>
*> \param[out] ALPHI2
*> \verbatim
*>          ALPHI2 is REAL array, dimension (max(NN))
*> \endverbatim
*>
*> \param[out] BETA2
*> \verbatim
*>          BETA2 is REAL array, dimension (max(NN))
*>
*>          The generalized eigenvalues of (A,B) computed by SGEGV.
*>          ( ALPHR2(k)+ALPHI2(k)*i ) / BETA2(k) is the k-th
*>          generalized eigenvalue of the matrices in A and B.
*> \endverbatim
*>
*> \param[out] VL
*> \verbatim
*>          VL is REAL array, dimension (LDQ, max(NN))
*>          The (block lower triangular) left eigenvector matrix for
*>          the matrices in A and B.  (See STGEVC for the format.)
*> \endverbatim
*>
*> \param[out] VR
*> \verbatim
*>          VR is REAL array, dimension (LDQ, max(NN))
*>          The (block upper triangular) right eigenvector matrix for
*>          the matrices in A and B.  (See STGEVC for the format.)
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (LWORK)
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The number of entries in WORK.  This must be at least
*>          2*N + MAX( 6*N, N*(NB+1), (k+1)*(2*k+N+1) ), where
*>          "k" is the sum of the blocksize and number-of-shifts for
*>          SHGEQZ, and NB is the greatest of the blocksizes for
*>          SGEQRF, SORMQR, and SORGQR.  (The blocksizes and the
*>          number-of-shifts are retrieved through calls to ILAENV.)
*> \endverbatim
*>
*> \param[out] RESULT
*> \verbatim
*>          RESULT is REAL array, dimension (15)
*>          The values computed by the tests described above.
*>          The values are currently limited to 1/ulp, to avoid
*>          overflow.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          > 0:  A routine returned an error code.  INFO is the
*>                absolute value of the INFO value returned.
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
*> \ingroup single_eig
*
*  =====================================================================
      SUBROUTINE SDRVGG( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH,
     $                   THRSHN, NOUNIT, A, LDA, B, S, T, S2, T2, Q,
     $                   LDQ, Z, ALPHR1, ALPHI1, BETA1, ALPHR2, ALPHI2,
     $                   BETA2, VL, VR, WORK, LWORK, RESULT, INFO )
*
*  -- LAPACK test routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDQ, LWORK, NOUNIT, NSIZES, NTYPES
      REAL               THRESH, THRSHN
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            ISEED( 4 ), NN( * )
      REAL               A( LDA, * ), ALPHI1( * ), ALPHI2( * ),
     $                   ALPHR1( * ), ALPHR2( * ), B( LDA, * ),
     $                   BETA1( * ), BETA2( * ), Q( LDQ, * ),
     $                   RESULT( * ), S( LDA, * ), S2( LDA, * ),
     $                   T( LDA, * ), T2( LDA, * ), VL( LDQ, * ),
     $                   VR( LDQ, * ), WORK( * ), Z( LDQ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0, ONE = 1.0 )
      INTEGER            MAXTYP
      PARAMETER          ( MAXTYP = 26 )
*     ..
*     .. Local Scalars ..
      LOGICAL            BADNN, ILABAD
      INTEGER            I1, IADD, IINFO, IN, J, JC, JR, JSIZE, JTYPE,
     $                   LWKOPT, MTYPES, N, N1, NB, NBZ, NERRS, NMATS,
     $                   NMAX, NS, NTEST, NTESTT
      REAL               SAFMAX, SAFMIN, TEMP1, TEMP2, ULP, ULPINV
*     ..
*     .. Local Arrays ..
      INTEGER            IASIGN( MAXTYP ), IBSIGN( MAXTYP ),
     $                   IOLDSD( 4 ), KADD( 6 ), KAMAGN( MAXTYP ),
     $                   KATYPE( MAXTYP ), KAZERO( MAXTYP ),
     $                   KBMAGN( MAXTYP ), KBTYPE( MAXTYP ),
     $                   KBZERO( MAXTYP ), KCLASS( MAXTYP ),
     $                   KTRIAN( MAXTYP ), KZ1( 6 ), KZ2( 6 )
      REAL               DUMMA( 4 ), RMAGN( 0: 3 )
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      REAL               SLAMCH, SLARND
      EXTERNAL           ILAENV, SLAMCH, SLARND
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALASVM, SGEGS, SGEGV, SGET51, SGET52, SGET53,
     $                   SLABAD, SLACPY, SLARFG, SLASET, SLATM4, SORM2R,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, REAL, SIGN
*     ..
*     .. Data statements ..
      DATA               KCLASS / 15*1, 10*2, 1*3 /
      DATA               KZ1 / 0, 1, 2, 1, 3, 3 /
      DATA               KZ2 / 0, 0, 1, 2, 1, 1 /
      DATA               KADD / 0, 0, 0, 0, 3, 2 /
      DATA               KATYPE / 0, 1, 0, 1, 2, 3, 4, 1, 4, 4, 1, 1, 4,
     $                   4, 4, 2, 4, 5, 8, 7, 9, 4*4, 0 /
      DATA               KBTYPE / 0, 0, 1, 1, 2, -3, 1, 4, 1, 1, 4, 4,
     $                   1, 1, -4, 2, -4, 8*8, 0 /
      DATA               KAZERO / 6*1, 2, 1, 2*2, 2*1, 2*2, 3, 1, 3,
     $                   4*5, 4*3, 1 /
      DATA               KBZERO / 6*1, 1, 2, 2*1, 2*2, 2*1, 4, 1, 4,
     $                   4*6, 4*4, 1 /
      DATA               KAMAGN / 8*1, 2, 3, 2, 3, 2, 3, 7*1, 2, 3, 3,
     $                   2, 1 /
      DATA               KBMAGN / 8*1, 3, 2, 3, 2, 2, 3, 7*1, 3, 2, 3,
     $                   2, 1 /
      DATA               KTRIAN / 16*0, 10*1 /
      DATA               IASIGN / 6*0, 2, 0, 2*2, 2*0, 3*2, 0, 2, 3*0,
     $                   5*2, 0 /
      DATA               IBSIGN / 7*0, 2, 2*0, 2*2, 2*0, 2, 0, 2, 9*0 /
*     ..
*     .. Executable Statements ..
*
*     Check for errors
*
      INFO = 0
*
      BADNN = .FALSE.
      NMAX = 1
      DO 10 J = 1, NSIZES
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ).LT.0 )
     $      BADNN = .TRUE.
   10 CONTINUE
*
*     Maximum blocksize and shift -- we assume that blocksize and number
*     of shifts are monotone increasing functions of N.
*
      NB = MAX( 1, ILAENV( 1, 'SGEQRF', ' ', NMAX, NMAX, -1, -1 ),
     $     ILAENV( 1, 'SORMQR', 'LT', NMAX, NMAX, NMAX, -1 ),
     $     ILAENV( 1, 'SORGQR', ' ', NMAX, NMAX, NMAX, -1 ) )
      NBZ = ILAENV( 1, 'SHGEQZ', 'SII', NMAX, 1, NMAX, 0 )
      NS = ILAENV( 4, 'SHGEQZ', 'SII', NMAX, 1, NMAX, 0 )
      I1 = NBZ + NS
      LWKOPT = 2*NMAX + MAX( 6*NMAX, NMAX*( NB+1 ),
     $         ( 2*I1+NMAX+1 )*( I1+1 ) )
*
*     Check for errors
*
      IF( NSIZES.LT.0 ) THEN
         INFO = -1
      ELSE IF( BADNN ) THEN
         INFO = -2
      ELSE IF( NTYPES.LT.0 ) THEN
         INFO = -3
      ELSE IF( THRESH.LT.ZERO ) THEN
         INFO = -6
      ELSE IF( LDA.LE.1 .OR. LDA.LT.NMAX ) THEN
         INFO = -10
      ELSE IF( LDQ.LE.1 .OR. LDQ.LT.NMAX ) THEN
         INFO = -19
      ELSE IF( LWKOPT.GT.LWORK ) THEN
         INFO = -30
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SDRVGG', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 )
     $   RETURN
*
      SAFMIN = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
      SAFMIN = SAFMIN / ULP
      SAFMAX = ONE / SAFMIN
      CALL SLABAD( SAFMIN, SAFMAX )
      ULPINV = ONE / ULP
*
*     The values RMAGN(2:3) depend on N, see below.
*
      RMAGN( 0 ) = ZERO
      RMAGN( 1 ) = ONE
*
*     Loop over sizes, types
*
      NTESTT = 0
      NERRS = 0
      NMATS = 0
*
      DO 170 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         N1 = MAX( 1, N )
         RMAGN( 2 ) = SAFMAX*ULP / REAL( N1 )
         RMAGN( 3 ) = SAFMIN*ULPINV*N1
*
         IF( NSIZES.NE.1 ) THEN
            MTYPES = MIN( MAXTYP, NTYPES )
         ELSE
            MTYPES = MIN( MAXTYP+1, NTYPES )
         END IF
*
         DO 160 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) )
     $         GO TO 160
            NMATS = NMATS + 1
            NTEST = 0
*
*           Save ISEED in case of an error.
*
            DO 20 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   20       CONTINUE
*
*           Initialize RESULT
*
            DO 30 J = 1, 15
               RESULT( J ) = ZERO
   30       CONTINUE
*
*           Compute A and B
*
*           Description of control parameters:
*
*           KCLASS: =1 means w/o rotation, =2 means w/ rotation,
*                   =3 means random.
*           KATYPE: the "type" to be passed to SLATM4 for computing A.
*           KAZERO: the pattern of zeros on the diagonal for A:
*                   =1: ( xxx ), =2: (0, xxx ) =3: ( 0, 0, xxx, 0 ),
*                   =4: ( 0, xxx, 0, 0 ), =5: ( 0, 0, 1, xxx, 0 ),
*                   =6: ( 0, 1, 0, xxx, 0 ).  (xxx means a string of
*                   non-zero entries.)
*           KAMAGN: the magnitude of the matrix: =0: zero, =1: O(1),
*                   =2: large, =3: small.
*           IASIGN: 1 if the diagonal elements of A are to be
*                   multiplied by a random magnitude 1 number, =2 if
*                   randomly chosen diagonal blocks are to be rotated
*                   to form 2x2 blocks.
*           KBTYPE, KBZERO, KBMAGN, IBSIGN: the same, but for B.
*           KTRIAN: =0: don't fill in the upper triangle, =1: do.
*           KZ1, KZ2, KADD: used to implement KAZERO and KBZERO.
*           RMAGN: used to implement KAMAGN and KBMAGN.
*
            IF( MTYPES.GT.MAXTYP )
     $         GO TO 110
            IINFO = 0
            IF( KCLASS( JTYPE ).LT.3 ) THEN
*
*              Generate A (w/o rotation)
*
               IF( ABS( KATYPE( JTYPE ) ).EQ.3 ) THEN
                  IN = 2*( ( N-1 ) / 2 ) + 1
                  IF( IN.NE.N )
     $               CALL SLASET( 'Full', N, N, ZERO, ZERO, A, LDA )
               ELSE
                  IN = N
               END IF
               CALL SLATM4( KATYPE( JTYPE ), IN, KZ1( KAZERO( JTYPE ) ),
     $                      KZ2( KAZERO( JTYPE ) ), IASIGN( JTYPE ),
     $                      RMAGN( KAMAGN( JTYPE ) ), ULP,
     $                      RMAGN( KTRIAN( JTYPE )*KAMAGN( JTYPE ) ), 2,
     $                      ISEED, A, LDA )
               IADD = KADD( KAZERO( JTYPE ) )
               IF( IADD.GT.0 .AND. IADD.LE.N )
     $            A( IADD, IADD ) = ONE
*
*              Generate B (w/o rotation)
*
               IF( ABS( KBTYPE( JTYPE ) ).EQ.3 ) THEN
                  IN = 2*( ( N-1 ) / 2 ) + 1
                  IF( IN.NE.N )
     $               CALL SLASET( 'Full', N, N, ZERO, ZERO, B, LDA )
               ELSE
                  IN = N
               END IF
               CALL SLATM4( KBTYPE( JTYPE ), IN, KZ1( KBZERO( JTYPE ) ),
     $                      KZ2( KBZERO( JTYPE ) ), IBSIGN( JTYPE ),
     $                      RMAGN( KBMAGN( JTYPE ) ), ONE,
     $                      RMAGN( KTRIAN( JTYPE )*KBMAGN( JTYPE ) ), 2,
     $                      ISEED, B, LDA )
               IADD = KADD( KBZERO( JTYPE ) )
               IF( IADD.NE.0 .AND. IADD.LE.N )
     $            B( IADD, IADD ) = ONE
*
               IF( KCLASS( JTYPE ).EQ.2 .AND. N.GT.0 ) THEN
*
*                 Include rotations
*
*                 Generate Q, Z as Householder transformations times
*                 a diagonal matrix.
*
                  DO 50 JC = 1, N - 1
                     DO 40 JR = JC, N
                        Q( JR, JC ) = SLARND( 3, ISEED )
                        Z( JR, JC ) = SLARND( 3, ISEED )
   40                CONTINUE
                     CALL SLARFG( N+1-JC, Q( JC, JC ), Q( JC+1, JC ), 1,
     $                            WORK( JC ) )
                     WORK( 2*N+JC ) = SIGN( ONE, Q( JC, JC ) )
                     Q( JC, JC ) = ONE
                     CALL SLARFG( N+1-JC, Z( JC, JC ), Z( JC+1, JC ), 1,
     $                            WORK( N+JC ) )
                     WORK( 3*N+JC ) = SIGN( ONE, Z( JC, JC ) )
                     Z( JC, JC ) = ONE
   50             CONTINUE
                  Q( N, N ) = ONE
                  WORK( N ) = ZERO
                  WORK( 3*N ) = SIGN( ONE, SLARND( 2, ISEED ) )
                  Z( N, N ) = ONE
                  WORK( 2*N ) = ZERO
                  WORK( 4*N ) = SIGN( ONE, SLARND( 2, ISEED ) )
*
*                 Apply the diagonal matrices
*
                  DO 70 JC = 1, N
                     DO 60 JR = 1, N
                        A( JR, JC ) = WORK( 2*N+JR )*WORK( 3*N+JC )*
     $                                A( JR, JC )
                        B( JR, JC ) = WORK( 2*N+JR )*WORK( 3*N+JC )*
     $                                B( JR, JC )
   60                CONTINUE
   70             CONTINUE
                  CALL SORM2R( 'L', 'N', N, N, N-1, Q, LDQ, WORK, A,
     $                         LDA, WORK( 2*N+1 ), IINFO )
                  IF( IINFO.NE.0 )
     $               GO TO 100
                  CALL SORM2R( 'R', 'T', N, N, N-1, Z, LDQ, WORK( N+1 ),
     $                         A, LDA, WORK( 2*N+1 ), IINFO )
                  IF( IINFO.NE.0 )
     $               GO TO 100
                  CALL SORM2R( 'L', 'N', N, N, N-1, Q, LDQ, WORK, B,
     $                         LDA, WORK( 2*N+1 ), IINFO )
                  IF( IINFO.NE.0 )
     $               GO TO 100
                  CALL SORM2R( 'R', 'T', N, N, N-1, Z, LDQ, WORK( N+1 ),
     $                         B, LDA, WORK( 2*N+1 ), IINFO )
                  IF( IINFO.NE.0 )
     $               GO TO 100
               END IF
            ELSE
*
*              Random matrices
*
               DO 90 JC = 1, N
                  DO 80 JR = 1, N
                     A( JR, JC ) = RMAGN( KAMAGN( JTYPE ) )*
     $                             SLARND( 2, ISEED )
                     B( JR, JC ) = RMAGN( KBMAGN( JTYPE ) )*
     $                             SLARND( 2, ISEED )
   80             CONTINUE
   90          CONTINUE
            END IF
*
  100       CONTINUE
*
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE,
     $            IOLDSD
               INFO = ABS( IINFO )
               RETURN
            END IF
*
  110       CONTINUE
*
*           Call SGEGS to compute H, T, Q, Z, alpha, and beta.
*
            CALL SLACPY( ' ', N, N, A, LDA, S, LDA )
            CALL SLACPY( ' ', N, N, B, LDA, T, LDA )
            NTEST = 1
            RESULT( 1 ) = ULPINV
*
            CALL SGEGS( 'V', 'V', N, S, LDA, T, LDA, ALPHR1, ALPHI1,
     $                  BETA1, Q, LDQ, Z, LDQ, WORK, LWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SGEGS', IINFO, N, JTYPE,
     $            IOLDSD
               INFO = ABS( IINFO )
               GO TO 140
            END IF
*
            NTEST = 4
*
*           Do tests 1--4
*
            CALL SGET51( 1, N, A, LDA, S, LDA, Q, LDQ, Z, LDQ, WORK,
     $                   RESULT( 1 ) )
            CALL SGET51( 1, N, B, LDA, T, LDA, Q, LDQ, Z, LDQ, WORK,
     $                   RESULT( 2 ) )
            CALL SGET51( 3, N, B, LDA, T, LDA, Q, LDQ, Q, LDQ, WORK,
     $                   RESULT( 3 ) )
            CALL SGET51( 3, N, B, LDA, T, LDA, Z, LDQ, Z, LDQ, WORK,
     $                   RESULT( 4 ) )
*
*           Do test 5: compare eigenvalues with diagonals.
*           Also check Schur form of A.
*
            TEMP1 = ZERO
*
            DO 120 J = 1, N
               ILABAD = .FALSE.
               IF( ALPHI1( J ).EQ.ZERO ) THEN
                  TEMP2 = ( ABS( ALPHR1( J )-S( J, J ) ) /
     $                    MAX( SAFMIN, ABS( ALPHR1( J ) ), ABS( S( J,
     $                    J ) ) )+ABS( BETA1( J )-T( J, J ) ) /
     $                    MAX( SAFMIN, ABS( BETA1( J ) ), ABS( T( J,
     $                    J ) ) ) ) / ULP
                  IF( J.LT.N ) THEN
                     IF( S( J+1, J ).NE.ZERO )
     $                  ILABAD = .TRUE.
                  END IF
                  IF( J.GT.1 ) THEN
                     IF( S( J, J-1 ).NE.ZERO )
     $                  ILABAD = .TRUE.
                  END IF
               ELSE
                  IF( ALPHI1( J ).GT.ZERO ) THEN
                     I1 = J
                  ELSE
                     I1 = J - 1
                  END IF
                  IF( I1.LE.0 .OR. I1.GE.N ) THEN
                     ILABAD = .TRUE.
                  ELSE IF( I1.LT.N-1 ) THEN
                     IF( S( I1+2, I1+1 ).NE.ZERO )
     $                  ILABAD = .TRUE.
                  ELSE IF( I1.GT.1 ) THEN
                     IF( S( I1, I1-1 ).NE.ZERO )
     $                  ILABAD = .TRUE.
                  END IF
                  IF( .NOT.ILABAD ) THEN
                     CALL SGET53( S( I1, I1 ), LDA, T( I1, I1 ), LDA,
     $                            BETA1( J ), ALPHR1( J ), ALPHI1( J ),
     $                            TEMP2, IINFO )
                     IF( IINFO.GE.3 ) THEN
                        WRITE( NOUNIT, FMT = 9997 )IINFO, J, N, JTYPE,
     $                     IOLDSD
                        INFO = ABS( IINFO )
                     END IF
                  ELSE
                     TEMP2 = ULPINV
                  END IF
               END IF
               TEMP1 = MAX( TEMP1, TEMP2 )
               IF( ILABAD ) THEN
                  WRITE( NOUNIT, FMT = 9996 )J, N, JTYPE, IOLDSD
               END IF
  120       CONTINUE
            RESULT( 5 ) = TEMP1
*
*           Call SGEGV to compute S2, T2, VL, and VR, do tests.
*
*           Eigenvalues and Eigenvectors
*
            CALL SLACPY( ' ', N, N, A, LDA, S2, LDA )
            CALL SLACPY( ' ', N, N, B, LDA, T2, LDA )
            NTEST = 6
            RESULT( 6 ) = ULPINV
*
            CALL SGEGV( 'V', 'V', N, S2, LDA, T2, LDA, ALPHR2, ALPHI2,
     $                  BETA2, VL, LDQ, VR, LDQ, WORK, LWORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUNIT, FMT = 9999 )'SGEGV', IINFO, N, JTYPE,
     $            IOLDSD
               INFO = ABS( IINFO )
               GO TO 140
            END IF
*
            NTEST = 7
*
*           Do Tests 6 and 7
*
            CALL SGET52( .TRUE., N, A, LDA, B, LDA, VL, LDQ, ALPHR2,
     $                   ALPHI2, BETA2, WORK, DUMMA( 1 ) )
            RESULT( 6 ) = DUMMA( 1 )
            IF( DUMMA( 2 ).GT.THRSHN ) THEN
               WRITE( NOUNIT, FMT = 9998 )'Left', 'SGEGV', DUMMA( 2 ),
     $            N, JTYPE, IOLDSD
            END IF
*
            CALL SGET52( .FALSE., N, A, LDA, B, LDA, VR, LDQ, ALPHR2,
     $                   ALPHI2, BETA2, WORK, DUMMA( 1 ) )
            RESULT( 7 ) = DUMMA( 1 )
            IF( DUMMA( 2 ).GT.THRESH ) THEN
               WRITE( NOUNIT, FMT = 9998 )'Right', 'SGEGV', DUMMA( 2 ),
     $            N, JTYPE, IOLDSD
            END IF
*
*           Check form of Complex eigenvalues.
*
            DO 130 J = 1, N
               ILABAD = .FALSE.
               IF( ALPHI2( J ).GT.ZERO ) THEN
                  IF( J.EQ.N ) THEN
                     ILABAD = .TRUE.
                  ELSE IF( ALPHI2( J+1 ).GE.ZERO ) THEN
                     ILABAD = .TRUE.
                  END IF
               ELSE IF( ALPHI2( J ).LT.ZERO ) THEN
                  IF( J.EQ.1 ) THEN
                     ILABAD = .TRUE.
                  ELSE IF( ALPHI2( J-1 ).LE.ZERO ) THEN
                     ILABAD = .TRUE.
                  END IF
               END IF
               IF( ILABAD ) THEN
                  WRITE( NOUNIT, FMT = 9996 )J, N, JTYPE, IOLDSD
               END IF
  130       CONTINUE
*
*           End of Loop -- Check for RESULT(j) > THRESH
*
  140       CONTINUE
*
            NTESTT = NTESTT + NTEST
*
*           Print out tests which fail.
*
            DO 150 JR = 1, NTEST
               IF( RESULT( JR ).GE.THRESH ) THEN
*
*                 If this is the first test to fail,
*                 print a header to the data file.
*
                  IF( NERRS.EQ.0 ) THEN
                     WRITE( NOUNIT, FMT = 9995 )'SGG'
*
*                    Matrix types
*
                     WRITE( NOUNIT, FMT = 9994 )
                     WRITE( NOUNIT, FMT = 9993 )
                     WRITE( NOUNIT, FMT = 9992 )'Orthogonal'
*
*                    Tests performed
*
                     WRITE( NOUNIT, FMT = 9991 )'orthogonal', '''',
     $                  'transpose', ( '''', J = 1, 5 )
*
                  END IF
                  NERRS = NERRS + 1
                  IF( RESULT( JR ).LT.10000.0 ) THEN
                     WRITE( NOUNIT, FMT = 9990 )N, JTYPE, IOLDSD, JR,
     $                  RESULT( JR )
                  ELSE
                     WRITE( NOUNIT, FMT = 9989 )N, JTYPE, IOLDSD, JR,
     $                  RESULT( JR )
                  END IF
               END IF
  150       CONTINUE
*
  160    CONTINUE
  170 CONTINUE
*
*     Summary
*
      CALL ALASVM( 'SGG', NOUNIT, NERRS, NTESTT, 0 )
      RETURN
*
 9999 FORMAT( ' SDRVGG: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
*
 9998 FORMAT( ' SDRVGG: ', A, ' Eigenvectors from ', A, ' incorrectly ',
     $      'normalized.', / ' Bits of error=', 0P, G10.3, ',', 9X,
     $      'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5,
     $      ')' )
*
 9997 FORMAT( ' SDRVGG: SGET53 returned INFO=', I1, ' for eigenvalue ',
     $      I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(',
     $      3( I5, ',' ), I5, ')' )
*
 9996 FORMAT( ' SDRVGG: S not in Schur form at eigenvalue ', I6, '.',
     $      / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ),
     $      I5, ')' )
*
 9995 FORMAT( / 1X, A3, ' -- Real Generalized eigenvalue problem driver'
     $       )
*
 9994 FORMAT( ' Matrix types (see SDRVGG for details): ' )
*
 9993 FORMAT( ' Special Matrices:', 23X,
     $      '(J''=transposed Jordan block)',
     $      / '   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J'',J'')  ',
     $      '6=(diag(J'',I), diag(I,J''))', / ' Diagonal Matrices:  ( ',
     $      'D=diag(0,1,2,...) )', / '   7=(D,I)   9=(large*D, small*I',
     $      ')  11=(large*I, small*D)  13=(large*D, large*I)', /
     $      '   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D) ',
     $      ' 14=(small*D, small*I)', / '  15=(D, reversed D)' )
 9992 FORMAT( ' Matrices Rotated by Random ', A, ' Matrices U, V:',
     $      / '  16=Transposed Jordan Blocks             19=geometric ',
     $      'alpha, beta=0,1', / '  17=arithm. alpha&beta             ',
     $      '      20=arithmetic alpha, beta=0,1', / '  18=clustered ',
     $      'alpha, beta=0,1            21=random alpha, beta=0,1',
     $      / ' Large & Small Matrices:', / '  22=(large, small)   ',
     $      '23=(small,large)    24=(small,small)    25=(large,large)',
     $      / '  26=random O(1) matrices.' )
*
 9991 FORMAT( / ' Tests performed:  (S is Schur, T is triangular, ',
     $      'Q and Z are ', A, ',', / 20X,
     $      'l and r are the appropriate left and right', / 19X,
     $      'eigenvectors, resp., a is alpha, b is beta, and', / 19X, A,
     $      ' means ', A, '.)', / ' 1 = | A - Q S Z', A,
     $      ' | / ( |A| n ulp )      2 = | B - Q T Z', A,
     $      ' | / ( |B| n ulp )', / ' 3 = | I - QQ', A,
     $      ' | / ( n ulp )             4 = | I - ZZ', A,
     $      ' | / ( n ulp )', /
     $      ' 5 = difference between (alpha,beta) and diagonals of',
     $      ' (S,T)', / ' 6 = max | ( b A - a B )', A,
     $      ' l | / const.   7 = max | ( b A - a B ) r | / const.',
     $      / 1X )
 9990 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=',
     $      4( I4, ',' ), ' result ', I3, ' is', 0P, F8.2 )
 9989 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=',
     $      4( I4, ',' ), ' result ', I3, ' is', 1P, E10.3 )
*
*     End of SDRVGG
*
      END
