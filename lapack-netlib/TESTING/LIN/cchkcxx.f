*> \brief \b CCHKCXX
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*      SUBROUTINE CCHKCXX( DOTYPE, NM, MVAL, NN, NVAL,
*     $                    NNB, NBVAL, NXVAL, THRESH, TSTERR,
*     $                    A, COPYA,
*     $                    C, COPYC, QRC, COPYQRC, X, COPYX, S, TAU,
*     $                    DESEL_ROWS, COPY_DESEL_ROWS,
*     $                    SEL_DESEL_COLS, COPY_SEL_DESEL_COLS,
*     $                    IPIV, COPY_IPIV, JPIV, COPY_JPIV,
*     $                    WORK, RWORK, IWORK, NOUT )
*      IMPLICIT NONE
*
*      .. Scalar Arguments ..
*      LOGICAL            TSTERR
*      INTEGER            NM, NN, NNB, NOUT
*      REAL               THRESH
*      ..
*      .. Array Arguments ..
*      LOGICAL            DOTYPE( * )
*      INTEGER            IWORK( * ), NBVAL( * ), MVAL( * ), NVAL( * ),
*     $                   NXVAL( * ),
*     $                   DESEL_ROWS( * ), COPY_DESEL_ROWS( * ),
*     $                   SEL_DESEL_COLS( * ), COPY_SEL_DESEL_COLS( * ),
*     $                   IPIV( * ), COPY_IPIV( * ),
*     $                   JPIV( * ), COPY_JPIV( * )
*      REAL               RWORK( * ), S( * )
*      COMPLEX            A( * ), COPYA( * ), C( * ), COPYC( * ),
*     $                   QRC( * ), COPYQRC( * ), X( * ), COPYX( * ),
*     $                   TAU( * ), WORK( * )
*     ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CCHKCXX tests CGECXX.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] DOTYPE
*> \verbatim
*>          DOTYPE is LOGICAL array, dimension (NTYPES)
*>          The matrix types to be used for testing.  Matrices of type j
*>          (for 1 <= j <= NTYPES) are used for testing if DOTYPE(j) =
*>          .TRUE.; if DOTYPE(j) = .FALSE., then type j is not used.
*> \endverbatim
*>
*> \param[in] NM
*> \verbatim
*>          NM is INTEGER
*>          The number of values of M contained in the vector MVAL.
*> \endverbatim
*>
*> \param[in] MVAL
*> \verbatim
*>          MVAL is INTEGER array, dimension (NM)
*>          The values of the matrix row dimension M.
*> \endverbatim
*>
*> \param[in] NN
*> \verbatim
*>          NN is INTEGER
*>          The number of values of N contained in the vector NVAL.
*> \endverbatim
*>
*> \param[in] NVAL
*> \verbatim
*>          NVAL is INTEGER array, dimension (NN)
*>          The values of the matrix column dimension N.
*> \endverbatim
*>
*> \param[in] NNB
*> \verbatim
*>          NNB is INTEGER
*>          The number of values of NB and NX contained in the
*>          vectors NBVAL and NXVAL.  The blocking parameters are used
*>          in pairs (NB,NX).
*> \endverbatim
*>
*> \param[in] NBVAL
*> \verbatim
*>          NBVAL is INTEGER array, dimension (NNB)
*>          The values of the blocksize NB.
*> \endverbatim
*>
*> \param[in] NXVAL
*> \verbatim
*>          NXVAL is INTEGER array, dimension (NNB)
*>          The values of the crossover point NX.
*> \endverbatim
*>
*> \param[in] THRESH
*> \verbatim
*>          THRESH is REAL
*>          The threshold value for the test ratios.  A result is
*>          included in the output file if RESULT >= THRESH.  To have
*>          every test ratio printed, use THRESH = 0.
*> \endverbatim
*>
*> \param[in] TSTERR
*> \verbatim
*>          TSTERR is LOGICAL
*>          Flag that indicates whether error exits are to be tested.
*> \endverbatim
*>
*> \param[out] A
*> \verbatim
*>          A is COMPLEX    array, dimension (MMAX*NMAX)
*>          where MMAX is the maximum value of M in MVAL and NMAX is the
*>          maximum value of N in NVAL.
*> \endverbatim
*>
*> \param[out] COPYA
*> \verbatim
*>          COPYA is COMPLEX array, dimension (MMAX*NMAX)
*>          where MMAX is the maximum value of M in MVAL and NMAX is the
*>          maximum value of N in NVAL.
*> \endverbatim
*>
*> \param[out] C
*> \verbatim
*>          C is COMPLEX array, dimension (MMAX*NMAX)
*>          where MMAX is the maximum value of M in MVAL and NMAX is the
*>          maximum value of N in NVAL.
*> \endverbatim
*>
*> \param[out] COPYC
*> \verbatim
*>          COPYC is COMPLEX array, dimension (MMAX*NMAX)
*>          where MMAX is the maximum value of M in MVAL and NMAX is the
*>          maximum value of N in NVAL.
*> \endverbatim
*>
*> \param[out] QRC
*> \verbatim
*>          QRC is COMPLEX array, dimension (MMAX*NMAX)
*>          where MMAX is the maximum value of M in MVAL and NMAX is the
*>          maximum value of N in NVAL.
*> \endverbatim
*>
*> \param[out] COPYQRC
*> \verbatim
*>          COPYQRC is COMPLEX array, dimension (MMAX*NMAX)
*>          where MMAX is the maximum value of M in MVAL and NMAX is the
*>          maximum value of N in NVAL.
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is COMPLEX array, dimension (NMAX*NMAX)
*>          NMAX is the maximum value of N in NVAL.
*> \endverbatim
*>
*> \param[out] COPYX
*> \verbatim
*>          COPYX is COMPLEX array, dimension (NMAX*NMAX)
*>          NMAX is the maximum value of N in NVAL.
*> \endverbatim
*>
*> \param[out] S
*> \verbatim
*>          S is REAL array, dimension
*>                      (min(MMAX,NMAX))
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is COMPLEX array, dimension (MMAX)
*> \endverbatim
*>
*> \param[out] DESEL_ROWS
*> \verbatim
*>          DESEL_ROWS is INTEGER array, dimension (MMAX)
*> \endverbatim
*>
*> \param[out] COPY_DESEL_ROWS
*> \verbatim
*>          COPY_DESEL_ROWS is INTEGER array, dimension (MMAX)
*> \endverbatim
*>
*> \param[out] SEL_DESEL_COLS
*> \verbatim
*>          SEL_DESEL_COLS is INTEGER array, dimension (NMAX)
*> \endverbatim
*>
*> \param[out] COPY_SEL_DESEL_COLS
*> \verbatim
*>          COPY_SEL_DESEL_COLS is INTEGER array, dimension (NMAX)
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (MMAX)
*> \endverbatim
*>
*> \param[out] COPY_IPIV
*> \verbatim
*>          COPY_IPIV is INTEGER array, dimension (MMAX)
*> \endverbatim
*>
*> \param[out] JPIV
*> \verbatim
*>          JPIV is INTEGER array, dimension (NMAX)
*> \endverbatim
*>
*> \param[out] COPY_JPIV
*> \verbatim
*>          COPY_JPIV is INTEGER array, dimension (NMAX)
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array.
*>          Dimension is the maximum of the following two expressions:
*>          (1) Optimal complex workspace dimension for matrix generation
*>              and test routines.
*>              (MMAX + 3) * max(MMAX,NMAX)
*>              This is an upper bound for:
*>                a) CLATMS: 3*max(M,N)
*>                b) CQRT12: M*N + 2*min(M,N) + max(M,N)
*>                c) CQPT01: M*N + N
*>                d) CQRT11: M*M + M
*>
*>          (2) Optimal complex workspace dimension for CGECXX.
*>              max( NMAX*NBMAX,                      \\ for CGEQRF inside
*>                   NMAX*min(NBMAX_UNMQR,NBMAX)      \\ for CUNMQR inside
*>                   + (NBMAX_UNMQR+1)*NBMAX_UNMQR ),
*>                   NBMAX*( NMAX + 1 ),              \\ for CGEQP3RK inside
*>                   min(MMAX,NMAX) + NMAX*NBMAX )    \\ for CGELS inside
*>              where NBMAX_UNMQR=64 is hardwired in CUNMQR.
*>
*>          Assuming MMAX = NMAX, and NBMAX = NMAX, the expressions become:
*>          (1) NMAX*NMAX + 3*NMAX
*>          (2) NMAX * min(64,NMAX) + 4160
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL array, dimension (2*NMAX)
*>
*>          Dimension is the maximum of the following two expressions:
*>          (1) Optimal real workspace dimension for matrix generation and test routines.
*>              2*min(MMAX,NMAX)
*>              This is an upper bound for CQRT12 routine 2*min(M,N).
*>
*>          (2) Optimal real workspace dimension for CGECXX.
*>              2*NMAX
*>              This is an upper bound for CGEQP3RK routine 2*N.
*>
*>          Assuming MMAX = NMAX, the expressions (1) anf (2) become 2*NMAX.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (2*NMAX)
*>                   for CGECXX optimal IWORK size.
*> \endverbatim
*>
*> \param[in] NOUT
*> \verbatim
*>          NOUT is INTEGER
*>          The unit number for output.
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
*> \ingroup complex_lin
*
*  =====================================================================
      SUBROUTINE CCHKCXX( DOTYPE, NM, MVAL, NN, NVAL,
     $                    NNB, NBVAL, NXVAL, THRESH, TSTERR,
     $                    A, COPYA,
     $                    C, COPYC, QRC, COPYQRC, X, COPYX, S, TAU,
     $                    DESEL_ROWS, COPY_DESEL_ROWS,
     $                    SEL_DESEL_COLS, COPY_SEL_DESEL_COLS,
     $                    IPIV, COPY_IPIV, JPIV, COPY_JPIV,
     $                    WORK, RWORK, IWORK, NOUT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      INTEGER            NM, NN, NNB, NOUT
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            IWORK( * ), NBVAL( * ), MVAL( * ), NVAL( * ),
     $                   NXVAL( * ),
     $                   DESEL_ROWS( * ), COPY_DESEL_ROWS( * ),
     $                   SEL_DESEL_COLS( * ), COPY_SEL_DESEL_COLS( * ),
     $                   IPIV( * ), COPY_IPIV( * ),
     $                   JPIV( * ), COPY_JPIV( * )
      REAL               RWORK( * ), S( * )
      COMPLEX            A( * ), COPYA( * ), C( * ), COPYC( * ),
     $                   QRC( * ), COPYQRC( * ), X( * ), COPYX( * ),
     $                   TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NTYPES
      PARAMETER          ( NTYPES = 19 )
      INTEGER            NTESTS
      PARAMETER          ( NTESTS = 5 )
      REAL               ONE, ZERO, BIGNUM
      COMPLEX            CZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0,
     $                     CZERO = ( 0.0E+0, 0.0E+0 ),
     $                     BIGNUM = 1.0E+38 )
*     ..
*     .. Local Scalars ..
      CHARACTER          DIST, TYPE, FACT, USESD
      CHARACTER*3        PATH
      INTEGER            I, IM, IMAT, IN, INB, IND_OFFSET_GEN,
     $                   IND_IN, IND_OUT, INFO, J, J_INC, J_FIRST_NZ,
     $                   JB_ZERO, K, KL, KMAXFREE, KU, LDA, LDC,
     $                   LDQRC, LDX, LIWORK, LRWORK, LWORK, LWKTST,
     $                   M, MINMN, MINMNB_GEN, MODE, N,
     $                   NB, NBMAX_UNMQR, NB_ZERO, NERRS, NFAIL,
     $                   NB_GEN, NRUN, NX, T
      REAL               ANORM, CNDNUM, EPS, ABSTOL, RELTOL,
     $                   DTEMP, MAXC2NRMK, RELMAXC2NRMK, FNRMK
*     ..
*     .. Local Arrays ..
      INTEGER            ISEED( 4 ), ISEEDY( 4 )
      REAL               RESULT( NTESTS )
*     ..
*     .. External Functions ..
      REAL               SLAMCH, CQPT01, CQRT11, CQRT12
      EXTERNAL           SLAMCH, CQPT01, CQRT11, CQRT12
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAERH, ALAHD, ALASUM, CERRCXX,
     $                   CGECXX, CLACPY, SLAORD, CLASET,
     $                   CLATB4, CLATMS, CSWAP, ICOPY, XLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, MOD
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*32       SRNAMT
      INTEGER            INFOT, IOUNIT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, IOUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'CX'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO I = 1, 4
         ISEED( I ) = ISEEDY( I )
      END DO
      EPS = SLAMCH( 'Epsilon' )
*
*     Test the error exits
*
      IF( TSTERR )
     $   CALL CERRCXX( PATH, NOUT )
*
      INFOT = 0
*
      DO IM = 1, NM
*
*        Do for each value of M in MVAL.
*
         M = MVAL( IM )
         LDA = MAX( 1, M )
         LDC = MAX( 1, M )
         LDQRC = MAX( 1, M )
*
         DO IN = 1, NN
*
*           Do for each value of N in NVAL.
*
            N = NVAL( IN )
            MINMN = MIN( M, N )
            LDX = MAX( 1, N )
*
*           1) NOTE: for matrix generation routine CLATMS, the workspace length
*           LWKTMS = 3*MAX( M, N ).  LWKTMS not used in the code.
*
*           2) Set workspace length for testing routines.
*             a) for CQRT12, real LRWKTST = 2*MIN(M,N). LRWKTST not used in the code.
*                for CQRT12, complex:
*
            LWKTST = MAX( 1, M*N + 2*MINMN + MAX( M, N ) )
*
*             b) for CQPT01, complex:
*
            LWKTST = MAX( LWKTST, M*N + N )
*
*             c) for CQRT11, complex:
*
            LWKTST = MAX( LWKTST, M*M + M )
*
            DO IMAT = 1, NTYPES
*
*              Do for each value of IMAT in NTYPES.
*
*              Do the tests only if DOTYPE( IMAT ) is true.
*
               IF( .NOT.DOTYPE( IMAT ) )
     $            CYCLE
*
*              The type of distribution used to generate the random
*              eigen-/singular values:
*              ( 'S' for symmetric distribution ) => UNIFORM( -1, 1 )
*
*           Do for each type of NON-SYMMETRIC matrix:                               CNDNUM                          NORM                                     MODE
*            1. Zero matrix                                                         CNDNUM = Inf                    0                                        N/A
*            2. Random, Diagonal                                                    CNDNUM = 2                      1                                        3 ( geometric distribution of singular values )
*            3. Random, Upper triangular                                            CNDNUM = 2                      1                                        3 ( geometric distribution of singular values )
*            4. Random, Lower triangular                                            CNDNUM = 2                      1                                        3 ( geometric distribution of singular values )
*            5. Random, First column is zero                                        CNDNUM = 2                      1                                        3 ( geometric distribution of singular values )
*            6. Random, Last MINMN column is zero                                   CNDNUM = 2                      1                                        3 ( geometric distribution of singular values )
*            7. Random, Last N column is zero                                       CNDNUM = 2                      1                                        3 ( geometric distribution of singular values )
*            8. Random, Middle column in MINMN is zero                              CNDNUM = 2                      1                                        3 ( geometric distribution of singular values )
*            9. Random, First half of MINMN columns are zero,
*                       zero block size MINMN/2                                     CNDNUM = 2                      1                                        3 ( geometric distribution of singular values )
*           10. Random, Last columns are zero starting from MINMN/2+1 column,
*                       zero block size N - MINMN/2                                 CNDNUM = 2                      1                                        3 ( geometric distribution of singular values )
*           11. Random, Half of MINMN columns in the middle are zero starting
*                       from  MINMN/2-(MINMN/2)/2+1 column,
*                       zero block size MINMN/2                                     CNDNUM = 2                      1                                        3 ( geometric distribution of singular values )
*           12. Random, Odd columns are ZERO                                        CNDNUM = 2                      1                                        3 ( geometric distribution of singular values )
*           13. Random, Even columns are ZERO                                       CNDNUM = 2                      1                                        3 ( geometric distribution of singular values )
*           14. Random, CNDNUM = 2                                                  CNDNUM = 2                      1                                        3 ( geometric distribution of singular values )
*           15. Random, CNDNUM = sqrt(0.1/EPS)                                      CNDNUM = BADC1 = sqrt(0.1/EPS)  1                                        3 ( geometric distribution of singular values )
*           16. Random, CNDNUM = 0.1/EPS                                            CNDNUM = BADC2 = 0.1/EPS        1                                        3 ( geometric distribution of singular values )
*           17. Random, CNDNUM = 0.1/EPS, one small singular value S(N)=1/CNDNUM    CNDNUM = BADC2 = 0.1/EPS        1                                        2 ( one small singular value, S(N)=1/CNDNUM )
*           18. Random, CNDNUM = 2, scaled near underflow                           CNDNUM = 2                      SMALL = SAFMIN                           3 ( geometric distribution of singular values )
*           19. Random, CNDNUM = 2, scaled near overflow                            CNDNUM = 2                      LARGE = 1.0/( 0.25 * ( SAFMIN / EPS ) )  3 ( geometric distribution of singular values )
*
*              Generate matrices.
*
               IF( IMAT.EQ.1 ) THEN
*
*                 Matrix 1 (Zero matrix).
*
                  CALL CLASET( 'Full', M, N, CZERO, CZERO, COPYA, LDA )
*
*                 Array S(1:min(M,N)) should contain svd(A), the sigular
*                 values of the generated matrix A in decreasing absolute
*                 value order. S in this format will be used later in the test.
*                 We set the array S explicitly here, since we are not using
*                 CLATMS (which sets the array S) to generate zero matrix.
*
                  DO I = 1, MINMN
                     S( I ) = ZERO
                  END DO
*
               ELSE IF( ( IMAT.EQ.2 .OR. IMAT.EQ.3 .OR. IMAT.EQ.4 )
     $            .OR.  ( IMAT.GE.14 .AND. IMAT.LE.19 ) ) THEN
*
*                 Matrix 2 (Diagonal),
*                 Matrix 3 (Upper triangular),
*                 Matrix 4 (Lower triangular),
*                 Matrices 14-19 (Various rectangular random matrices
*                          without zero columns).
*
*                 Set up parameters with CLATB4 and generate a test
*                 matrix with CLATMS.
*
                  CALL CLATB4( PATH, IMAT, M, N, TYPE, KL, KU, ANORM,
     $                         MODE, CNDNUM, DIST )
*
                  SRNAMT = 'CLATMS'
                  CALL CLATMS( M, N, DIST, ISEED, TYPE, S, MODE,
     $                         CNDNUM, ANORM, KL, KU, 'No packing',
     $                         COPYA, LDA, WORK, INFO )
*
*                 Check error code from CLATMS.
*
                  IF( INFO.NE.0 ) THEN
                     CALL ALAERH( PATH, 'CLATMS', INFO, 0, ' ', M, N,
     $                            -1, -1, -1, IMAT, NFAIL, NERRS,
     $                            NOUT )
                     CYCLE
                  END IF
*
*                 Array S(1:min(M,N)) should contain svd(A), the sigular
*                 values of the generated matrix A in decreasing absolute
*                 value order. S in this format will be used later in
*                 the test. Unordered singular values are returned by
*                 CLATMS in S. We need to order singular values in S.
*
                  CALL SLAORD( 'Decreasing', MINMN, S, 1 )
*
               ELSE IF( MINMN.GE.2
     $                  .AND. IMAT.GE.5 .AND. IMAT.LE.13 ) THEN
*
*                 Matrices 5-13 (Rectangular random matrices that
*                 contain zero columns). Only for matrices MINMN >= 2.
*
*                 JB_ZERO is the column index of ZERO block.
*                 NB_ZERO is the column block size of ZERO block.
*                 NB_GEN is the column blcok size of the
*                 generated block.
*                 J_INC in the non_zero column index increment
*                 to generate matrix 12 and 13.
*                 J_FIRS_NZ is the index of the first non-zero
*                 column to generate matrix 12 and 13.
*
                  IF( IMAT.EQ.5 ) THEN
*
*                    Matrix 5. First column is zero.
*
                     JB_ZERO = 1
                     NB_ZERO = 1
                     NB_GEN = N - NB_ZERO
*
                  ELSE IF( IMAT.EQ.6 ) THEN
*
*                    Matrix 6. Last column MINMN is zero.
*
                     JB_ZERO = MINMN
                     NB_ZERO = 1
                     NB_GEN = N - NB_ZERO
*
                  ELSE IF( IMAT.EQ.7 ) THEN
*
*                    Matrix 7. Last column N is zero.
*
                     JB_ZERO = N
                     NB_ZERO = 1
                     NB_GEN = N - NB_ZERO
*
                  ELSE IF( IMAT.EQ.8 ) THEN
*
*                    MAtrix 8. Middle column in MINMN is zero.
*
                     JB_ZERO = MINMN / 2 + 1
                     NB_ZERO = 1
                     NB_GEN = N - NB_ZERO
*
                  ELSE IF( IMAT.EQ.9 ) THEN
*
*                    Matrix 9. First half of MINMN columns is zero, zero block size MINMN/2.
*
                     JB_ZERO = 1
                     NB_ZERO = MINMN / 2
                     NB_GEN = N - NB_ZERO
*
                  ELSE IF( IMAT.EQ.10 ) THEN
*
*                    Matrix 10. Last columns are zero columns,
*                    starting from (MINMN / 2 + 1) column,zero block size N - MINMN/2
*
                     JB_ZERO = MINMN / 2 + 1
                     NB_ZERO = N - MINMN / 2
                     NB_GEN = N - NB_ZERO
*
                  ELSE IF( IMAT.EQ.11 ) THEN
*
*                    Matrix 11. Half of the columns in the middle of first MINMN
*                    columns is zero, starting from MINMN/2 - (MINMN/2)/2 + 1 column,
*                    zero block size MINMN/2.
*
                     JB_ZERO = MINMN / 2 - (MINMN / 2) / 2 + 1
                     NB_ZERO = MINMN / 2
                     NB_GEN = N - NB_ZERO
*
                  ELSE IF( IMAT.EQ.12 ) THEN
*
*                    Matrix 12. Odd-numbered columns are zero,
*
                     NB_GEN = N / 2
                     NB_ZERO = N - NB_GEN
                     J_INC = 2
                     J_FIRST_NZ = 2
*
                  ELSE IF( IMAT.EQ.13 ) THEN
*
*                    Matrix 13. Even-numbered columns are zero.
*
                     NB_ZERO = N / 2
                     NB_GEN = N - NB_ZERO
                     J_INC = 2
                     J_FIRST_NZ = 1
*
                  END IF
*
*
*                 1) Set the first NB_ZERO columns in COPYA(1:M,1:N)
*                    to zero.
*
                  CALL CLASET( 'Full', M, NB_ZERO, CZERO, CZERO,
     $                         COPYA, LDA )
*
*                 2) Generate an M-by-(N-NB_ZERO) matrix with the
*                    chosen singular value distribution
*                    in COPYA(1:M,NB_ZERO+1:N).
*
                  CALL CLATB4( PATH, IMAT, M, NB_GEN, TYPE, KL, KU,
     $                         ANORM, MODE, CNDNUM, DIST )
*
                  SRNAMT = 'CLATMS'
*
                  IND_OFFSET_GEN = NB_ZERO * LDA
*
                  CALL CLATMS( M, NB_GEN, DIST, ISEED, TYPE, S, MODE,
     $                        CNDNUM, ANORM, KL, KU, 'No packing',
     $                        COPYA( IND_OFFSET_GEN + 1 ), LDA,
     $                        WORK, INFO )
*
*                 Check error code from CLATMS.
*
                  IF( INFO.NE.0 ) THEN
                     CALL ALAERH( PATH, 'CLATMS', INFO, 0, ' ', M,
     $                            NB_GEN, -1, -1, -1, IMAT, NFAIL,
     $                            NERRS, NOUT )
                     CYCLE
                  END IF
*
*                 3) Swap the gererated colums from the right side
*                 NB_GEN-size block in COPYA into correct column
*                 positions.
*
                  IF( IMAT.EQ.6
     $                    .OR. IMAT.EQ.7
     $                    .OR. IMAT.EQ.8
     $                    .OR. IMAT.EQ.10
     $                    .OR. IMAT.EQ.11 ) THEN
*
*                    Move by swapping the generated columns
*                    from the right NB_GEN-size block from
*                    (NB_ZERO+1:NB_ZERO+JB_ZERO)
*                    into columns (1:JB_ZERO-1).
*
                     DO J = 1, JB_ZERO-1, 1
                        CALL CSWAP( M,
     $                        COPYA( ( NB_ZERO+J-1)*LDA+1), 1,
     $                        COPYA( (J-1)*LDA + 1 ), 1 )
                     END DO
*
                  ELSE IF( IMAT.EQ.12 .OR. IMAT.EQ.13 ) THEN
*
*                    ( IMAT = 12, Odd-numbered ZERO columns. )
*                    Swap the generated columns from the right
*                    NB_GEN-size block into the even zero colums in the
*                    left NB_ZERO-size block.
*
*                    ( IMAT = 13, Even-numbered ZERO columns. )
*                    Swap the generated columns from the right
*                    NB_GEN-size block into the odd zero colums in the
*                    left NB_ZERO-size block.
*
                     DO J = 1, NB_GEN, 1
                        IND_OUT = ( NB_ZERO+J-1 )*LDA + 1
                        IND_IN = ( J_INC*(J-1)+(J_FIRST_NZ-1) )*LDA
     $                            + 1
                        CALL CSWAP( M,
     $                              COPYA( IND_OUT ), 1,
     $                              COPYA( IND_IN ), 1 )
                        END DO
*
                  END IF
*
*                 5) Order the singular values generated by
*                    CLAMTS in decreasing absolute value order and
*                    add trailing zeros that correspond to zero columns.
*                    The total number of singular values is MINMN.
*
                  MINMNB_GEN = MIN( M, NB_GEN )
                  CALL SLAORD( 'Decreasing', MINMNB_GEN, S, 1 )
*
                  DO I = MINMNB_GEN+1, MINMN
                     S( I ) = ZERO
                  END DO
*
               ELSE
*
*                 IF( MINMN.LT.2 .AND. ( IMAT.GE.5 .AND. IMAT.LE.13 ) )
*                 skip this size for this matrix type.
*
                  CYCLE
               END IF
*
*              End generate COPYA matrix.
*
*              Initialize COPYC matrix with zeros.
*
               CALL CLASET( 'Full', M, N, CZERO, CZERO,
     $                      COPYC, LDC )
*
*              Initialize COPYQRC matrix with zeros.
*
               CALL CLASET( 'Full', M, N, CZERO, CZERO,
     $                      COPYQRC, LDQRC )
*
*              Initialize COPYX matrix with zeros.
*
               CALL CLASET( 'Full', MINMN, N, CZERO, CZERO,
     $                      COPYX, LDX )
*
*              Initialize a copy array for pivot IPIV for CGECXX.
*
               DO I = 1, M
                  COPY_IPIV( I ) = 0
               END DO
*
*              Initialize a copy array for pivot JPIV for CGECXX.
*
               DO J = 1, N
                  COPY_JPIV( J ) = 0
               END DO
*
*              Initialize a copy array COPY_DESEL_ROWS for CGECXX.
*
               DO I = 1, M
                  COPY_DESEL_ROWS( I ) = 0
               END DO
*
*              Initialize a copy array COPY_SEL_DESEL_COLS for CGECXX.
*
               DO J = 1, N
                  COPY_SEL_DESEL_COLS( J ) = 0
               END DO
*
               DO INB = 1, NNB
*
*                 Do for each pair of values (NB,NX) in NBVAL and NXVAL.
*
                  NB = NBVAL( INB )
                  CALL XLAENV( 1, NB )
                  NX = NXVAL( INB )
                  CALL XLAENV( 3, NX )
*
*                 We do MIN(M,N)+1 because we need a test for KMAX > N,
*                 when KMAX is larger than MIN(M,N), KMAX should be
*                 KMAX = MIN(M,N)
*
                  DO KMAXFREE = 0, MIN(M,N)+1
*
*                 Get a working copy of COPYA into A( 1:M,1:N ).
*                 Get a working copy of COPYC into C( 1:M,1:N ).
*                 Get a working copy of COPYQRC into QRC( 1:M,1:N ).
*                 Get a working copy of COPYX into X( 1:N,1:N ).
*                 Get a working copy of COPY_IPIV(1:M) into IPIV(1:M).
*                 Get a working copy of COPY_JPIV(1:N) into JPIV(1:N).
*                 Get a working copy of COPY_DESEL_ROWS(1:M) into DESEL_ROWS(1:M).
*                 Get a working copy of COPY_SEL_DESEL_COLS(1:N) into SEL_DESEL_COLS(1:N).
*
                  CALL CLACPY( 'All', M, N, COPYA, LDA, A, LDA )
                  CALL CLACPY( 'All', M, N, COPYC, LDC, C, LDC )
                  CALL CLACPY( 'All', M, N, COPYQRC, LDQRC, QRC, LDQRC )
                  CALL CLACPY( 'All', MINMN, N, COPYX, LDX, X, LDX )
                  CALL ICOPY( M, COPY_IPIV, 1, IPIV, 1 )
                  CALL ICOPY( N, COPY_JPIV, 1, JPIV, 1 )
                  CALL ICOPY( M, COPY_DESEL_ROWS, 1, DESEL_ROWS, 1 )
                  CALL ICOPY( N, COPY_SEL_DESEL_COLS, 1,
     $                           SEL_DESEL_COLS, 1 )
*
*                 Set test ratios for all tests to zero.
*
                  DO I = 1, NTESTS
                    RESULT( I ) = ZERO
                  END DO
*
*                 We are not testing with ABSTOL and RELTOL stopping criteria.
*                 Disable them.
*
                  FACT = 'C'
                  USESD = 'N'
                  ABSTOL = -ONE
                  RELTOL = -ONE
*
*                 Compute the QR factorization with pivoting of A
*
*                 Determine LWORK
*
*                 NBMAX_UNMQR is hardwired in CUNMQR as NBMAX = 64.
*
                  NBMAX_UNMQR = 64
*
*                 a) For CGEQRF inside CGECXX, complex
*
                  LWORK = MAX( 1, N*NB )
*
*                 b) For CUNMQR inside CGECXX, complex
*
                  LWORK = MAX( LWORK,
     $             N*MIN(NBMAX_UNMQR,NB)+(NBMAX_UNMQR+1)*NBMAX_UNMQR )
*
*                 c1) For CGEQP3RK inside CGECXX, complex
*
                  LWORK = MAX( LWORK, NB*( N + 1 ) )
*
*                 c2) For CGEQP3RK inside CGECXX, real
*
                  LRWORK = MAX( 1, 2*N )
*
*                 d) For CGELS inside CGECXX, complex
*
                  LWORK = MAX( LWORK, MIN(M,N) + N*NB )
*
*                 Determine LIWORK
*
                  LIWORK = MAX( 1, 2*N )
*
*                 Compute CGECXX factorization of A.
*
                  SRNAMT = 'CGECXX'
                  CALL CGECXX( FACT, USESD, M, N,
     $                         DESEL_ROWS, SEL_DESEL_COLS,
     $                         KMAXFREE, ABSTOL, RELTOL, A, LDA,
     $                         K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $                         IPIV, JPIV, TAU, C, LDC, QRC, LDQRC,
     $                         X, LDX, WORK, LWORK, RWORK, LRWORK,
     $                         IWORK, LIWORK, INFO )
*
*                 Check an error code from CGECXX.
*
                  IF( INFO.LT.0 )
     $               CALL ALAERH( PATH, 'CGECXX', INFO, 0, ' ',
     $                            M, N, NX, -1, NB, IMAT,
     $                            NFAIL, NERRS, NOUT )
*
*                 Compute test 1:
*
*                 This test in only for the full rank factorization of
*                 the matrix A.
*
*                 Array S(1:min(M,N)) contains svd(A) the sigular values
*                 of the original matrix A in decreasing absolute value
*                 order. The test computes svd(R), the vector sigular
*                 values of the upper trapezoid of A(1:M,1:N) that
*                 contains the factor R, in decreasing order. The test
*                 returns the ratio:
*
*                 2-norm(svd(R) - svd(A)) / ( max(M,N) * 2-norm(svd(A)) * EPS )
*
                  IF( K.EQ.MINMN ) THEN
*
                     RESULT( 1 ) = CQRT12( M, N, A, LDA, S, WORK,
     $                                     LWKTST, RWORK )
*
                     NRUN = NRUN + 1
*
*                   End test 1
*
                  END IF
*
*
*                 Compute test 2:
*
*                 The test returns the ratio:
*
*                 1-norm( A*P - Q*R ) / ( max(M,N) * 1-norm(A) * EPS )
*
                  RESULT( 2 ) = CQPT01( M, N, K, COPYA, A, LDA, TAU,
     $                                  JPIV, WORK, LWKTST )
*
*                 Compute test 3:
*
*                 The test returns the ratio:
*
*                 1-norm( Q**T * Q - I ) / ( M * EPS )
*
                  RESULT( 3 ) = CQRT11( M, K, A, LDA, TAU, WORK,
     $                                  LWKTST )
*
                  NRUN = NRUN + 2
*
*                 Compute test 4:
*
*                 This test is only for the factorizations with the
*                 rank greater then 1.
*                 The elements on the diagonal of R should be non-
*                 increasing.
*
*                 The test returns the ratio:
*
*                 Returns 1.0E+38 if abs(R(j+1,j+1)) > abs(R(j,j)),
*                 j=1:K-1
*
                  IF( MIN(K, MINMN).GT.1 ) THEN
*
                     DO J = 1, K-1, 1

                        DTEMP = (( ABS( A( (J-1)*LDA+J ) ) -
     $                          ABS( A( (J)*LDA+J+1 ) ) ) /
     $                          ABS( A(1) ) )
*
                        IF( DTEMP.LT.ZERO ) THEN
                           RESULT( 4 ) = BIGNUM
                        END IF
*
                     END DO
*
                     NRUN = NRUN + 1
*
*                    End test 4.
*
                  END IF
*
*                 ===============
*                 Compute test 5:
*                 ===============
*                 This test is only for the factorizations with the
*                 rank greater than 0.
*                 For J=1:K, the J-th column of C should be elementwise
*                 equal (including NaN and Inf)
*                 to the JPIV(J)-th column of A.
*
                  RESULT( 5 ) = ZERO
*                 Disable for now, incomplete test.
                  IF(.FALSE.) THEN
                  DO J = 1, K, 1
                     DO I = 1, M, 1
                        IF( .NOT. (C( (J-1)*LDC+I )
     $                     .EQ. A( (JPIV( J )-1)*LDA+I ) ) ) THEN
                           RESULT( 5 ) = BIGNUM
                        END IF
                     END DO
                  END DO
                  END IF
*
*
*                 Print information about the tests that did not
*                 pass the threshold.
*
                  DO T = 1, NTESTS
                     IF( RESULT( T ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                     CALL ALAHD( NOUT, PATH )
                        WRITE( NOUT, FMT = 9999 ) 'CGECXX', M, N,
     $                     FACT, USESD, KMAXFREE, ABSTOL, RELTOL,
     $                     NB, NX, IMAT, T, RESULT( T )
                        NFAIL = NFAIL + 1
                     END IF
                  END DO
*
*                 END DO KMAX = 1, MIN(M,N)+1
*
                  END DO
*
*                 END DO for INB = 1, NNB
*
               END DO
*
*              END DO  for IMAT = 1, NTYPES
*
            END DO
*
*           END DO for IN = 1, NN
*
         END DO
*
*        END DO for IM = 1, NM
*
      END DO
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( 1X, A, ' M =', I5, ', N =', I5,
     $        ', FACT = ''', A1, ''', USESD = ''', A1,
     $        ''', KMAXFREE =', I5, ', ABSTOL =', G12.5,
     $        ', RELTOL =', G12.5, ', NB =', I4, ', NX =', I4,
     $        ', type ', I2, ', test ', I2, ', ratio =', G12.5 )
*
*     End of CCHKCXX
*
      END
