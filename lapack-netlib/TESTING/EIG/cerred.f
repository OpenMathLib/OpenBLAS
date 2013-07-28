*> \brief \b CERRED
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE CERRED( PATH, NUNIT )
* 
*       .. Scalar Arguments ..
*       CHARACTER*3        PATH
*       INTEGER            NUNIT
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CERRED tests the error exits for the eigenvalue driver routines for
*> REAL matrices:
*>
*> PATH  driver   description
*> ----  ------   -----------
*> CEV   CGEEV    find eigenvalues/eigenvectors for nonsymmetric A
*> CES   CGEES    find eigenvalues/Schur form for nonsymmetric A
*> CVX   CGEEVX   CGEEV + balancing and condition estimation
*> CSX   CGEESX   CGEES + balancing and condition estimation
*> CBD   CGESVD   compute SVD of an M-by-N matrix A
*>       CGESDD   compute SVD of an M-by-N matrix A(by divide and
*>                conquer)
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] PATH
*> \verbatim
*>          PATH is CHARACTER*3
*>          The LAPACK path name for the routines to be tested.
*> \endverbatim
*>
*> \param[in] NUNIT
*> \verbatim
*>          NUNIT is INTEGER
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
*> \date November 2011
*
*> \ingroup complex_eig
*
*  =====================================================================
      SUBROUTINE CERRED( PATH, NUNIT )
*
*  -- LAPACK test routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER*3        PATH
      INTEGER            NUNIT
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NMAX, LW
      PARAMETER          ( NMAX = 4, LW = 5*NMAX )
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E0, ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      CHARACTER*2        C2
      INTEGER            I, IHI, ILO, INFO, J, NT, SDIM
      REAL               ABNRM
*     ..
*     .. Local Arrays ..
      LOGICAL            B( NMAX )
      INTEGER            IW( 4*NMAX )
      REAL               R1( NMAX ), R2( NMAX ), RW( LW ), S( NMAX )
      COMPLEX            A( NMAX, NMAX ), U( NMAX, NMAX ),
     $                   VL( NMAX, NMAX ), VR( NMAX, NMAX ),
     $                   VT( NMAX, NMAX ), W( 4*NMAX ), X( NMAX )
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEES, CGEESX, CGEEV, CGEEVX, CGESDD, CGESVD,
     $                   CHKXER
*     ..
*     .. External Functions ..
      LOGICAL            CSLECT, LSAMEN
      EXTERNAL           CSLECT, LSAMEN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LEN_TRIM
*     ..
*     .. Arrays in Common ..
      LOGICAL            SELVAL( 20 )
      REAL               SELWI( 20 ), SELWR( 20 )
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*32       SRNAMT
      INTEGER            INFOT, NOUT, SELDIM, SELOPT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      COMMON             / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )
*
*     Initialize A
*
      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            A( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
      DO 30 I = 1, NMAX
         A( I, I ) = ONE
   30 CONTINUE
      OK = .TRUE.
      NT = 0
*
      IF( LSAMEN( 2, C2, 'EV' ) ) THEN
*
*        Test CGEEV
*
         SRNAMT = 'CGEEV '
         INFOT = 1
         CALL CGEEV( 'X', 'N', 0, A, 1, X, VL, 1, VR, 1, W, 1, RW,
     $               INFO )
         CALL CHKXER( 'CGEEV ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGEEV( 'N', 'X', 0, A, 1, X, VL, 1, VR, 1, W, 1, RW,
     $               INFO )
         CALL CHKXER( 'CGEEV ', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGEEV( 'N', 'N', -1, A, 1, X, VL, 1, VR, 1, W, 1, RW,
     $               INFO )
         CALL CHKXER( 'CGEEV ', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGEEV( 'N', 'N', 2, A, 1, X, VL, 1, VR, 1, W, 4, RW,
     $               INFO )
         CALL CHKXER( 'CGEEV ', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CGEEV( 'V', 'N', 2, A, 2, X, VL, 1, VR, 1, W, 4, RW,
     $               INFO )
         CALL CHKXER( 'CGEEV ', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL CGEEV( 'N', 'V', 2, A, 2, X, VL, 1, VR, 1, W, 4, RW,
     $               INFO )
         CALL CHKXER( 'CGEEV ', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL CGEEV( 'V', 'V', 1, A, 1, X, VL, 1, VR, 1, W, 1, RW,
     $               INFO )
         CALL CHKXER( 'CGEEV ', INFOT, NOUT, LERR, OK )
         NT = NT + 7
*
      ELSE IF( LSAMEN( 2, C2, 'ES' ) ) THEN
*
*        Test CGEES
*
         SRNAMT = 'CGEES '
         INFOT = 1
         CALL CGEES( 'X', 'N', CSLECT, 0, A, 1, SDIM, X, VL, 1, W, 1,
     $               RW, B, INFO )
         CALL CHKXER( 'CGEES ', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGEES( 'N', 'X', CSLECT, 0, A, 1, SDIM, X, VL, 1, W, 1,
     $               RW, B, INFO )
         CALL CHKXER( 'CGEES ', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CGEES( 'N', 'S', CSLECT, -1, A, 1, SDIM, X, VL, 1, W, 1,
     $               RW, B, INFO )
         CALL CHKXER( 'CGEES ', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CGEES( 'N', 'S', CSLECT, 2, A, 1, SDIM, X, VL, 1, W, 4,
     $               RW, B, INFO )
         CALL CHKXER( 'CGEES ', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL CGEES( 'V', 'S', CSLECT, 2, A, 2, SDIM, X, VL, 1, W, 4,
     $               RW, B, INFO )
         CALL CHKXER( 'CGEES ', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL CGEES( 'N', 'S', CSLECT, 1, A, 1, SDIM, X, VL, 1, W, 1,
     $               RW, B, INFO )
         CALL CHKXER( 'CGEES ', INFOT, NOUT, LERR, OK )
         NT = NT + 6
*
      ELSE IF( LSAMEN( 2, C2, 'VX' ) ) THEN
*
*        Test CGEEVX
*
         SRNAMT = 'CGEEVX'
         INFOT = 1
         CALL CGEEVX( 'X', 'N', 'N', 'N', 0, A, 1, X, VL, 1, VR, 1, ILO,
     $                IHI, S, ABNRM, R1, R2, W, 1, RW, INFO )
         CALL CHKXER( 'CGEEVX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGEEVX( 'N', 'X', 'N', 'N', 0, A, 1, X, VL, 1, VR, 1, ILO,
     $                IHI, S, ABNRM, R1, R2, W, 1, RW, INFO )
         CALL CHKXER( 'CGEEVX', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGEEVX( 'N', 'N', 'X', 'N', 0, A, 1, X, VL, 1, VR, 1, ILO,
     $                IHI, S, ABNRM, R1, R2, W, 1, RW, INFO )
         CALL CHKXER( 'CGEEVX', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CGEEVX( 'N', 'N', 'N', 'X', 0, A, 1, X, VL, 1, VR, 1, ILO,
     $                IHI, S, ABNRM, R1, R2, W, 1, RW, INFO )
         CALL CHKXER( 'CGEEVX', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGEEVX( 'N', 'N', 'N', 'N', -1, A, 1, X, VL, 1, VR, 1,
     $                ILO, IHI, S, ABNRM, R1, R2, W, 1, RW, INFO )
         CALL CHKXER( 'CGEEVX', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CGEEVX( 'N', 'N', 'N', 'N', 2, A, 1, X, VL, 1, VR, 1, ILO,
     $                IHI, S, ABNRM, R1, R2, W, 4, RW, INFO )
         CALL CHKXER( 'CGEEVX', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL CGEEVX( 'N', 'V', 'N', 'N', 2, A, 2, X, VL, 1, VR, 1, ILO,
     $                IHI, S, ABNRM, R1, R2, W, 4, RW, INFO )
         CALL CHKXER( 'CGEEVX', INFOT, NOUT, LERR, OK )
         INFOT = 12
         CALL CGEEVX( 'N', 'N', 'V', 'N', 2, A, 2, X, VL, 1, VR, 1, ILO,
     $                IHI, S, ABNRM, R1, R2, W, 4, RW, INFO )
         CALL CHKXER( 'CGEEVX', INFOT, NOUT, LERR, OK )
         INFOT = 20
         CALL CGEEVX( 'N', 'N', 'N', 'N', 1, A, 1, X, VL, 1, VR, 1, ILO,
     $                IHI, S, ABNRM, R1, R2, W, 1, RW, INFO )
         CALL CHKXER( 'CGEEVX', INFOT, NOUT, LERR, OK )
         INFOT = 20
         CALL CGEEVX( 'N', 'N', 'V', 'V', 1, A, 1, X, VL, 1, VR, 1, ILO,
     $                IHI, S, ABNRM, R1, R2, W, 2, RW, INFO )
         CALL CHKXER( 'CGEEVX', INFOT, NOUT, LERR, OK )
         NT = NT + 10
*
      ELSE IF( LSAMEN( 2, C2, 'SX' ) ) THEN
*
*        Test CGEESX
*
         SRNAMT = 'CGEESX'
         INFOT = 1
         CALL CGEESX( 'X', 'N', CSLECT, 'N', 0, A, 1, SDIM, X, VL, 1,
     $                R1( 1 ), R2( 1 ), W, 1, RW, B, INFO )
         CALL CHKXER( 'CGEESX', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGEESX( 'N', 'X', CSLECT, 'N', 0, A, 1, SDIM, X, VL, 1,
     $                R1( 1 ), R2( 1 ), W, 1, RW, B, INFO )
         CALL CHKXER( 'CGEESX', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CGEESX( 'N', 'N', CSLECT, 'X', 0, A, 1, SDIM, X, VL, 1,
     $                R1( 1 ), R2( 1 ), W, 1, RW, B, INFO )
         CALL CHKXER( 'CGEESX', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGEESX( 'N', 'N', CSLECT, 'N', -1, A, 1, SDIM, X, VL, 1,
     $                R1( 1 ), R2( 1 ), W, 1, RW, B, INFO )
         CALL CHKXER( 'CGEESX', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL CGEESX( 'N', 'N', CSLECT, 'N', 2, A, 1, SDIM, X, VL, 1,
     $                R1( 1 ), R2( 1 ), W, 4, RW, B, INFO )
         CALL CHKXER( 'CGEESX', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL CGEESX( 'V', 'N', CSLECT, 'N', 2, A, 2, SDIM, X, VL, 1,
     $                R1( 1 ), R2( 1 ), W, 4, RW, B, INFO )
         CALL CHKXER( 'CGEESX', INFOT, NOUT, LERR, OK )
         INFOT = 15
         CALL CGEESX( 'N', 'N', CSLECT, 'N', 1, A, 1, SDIM, X, VL, 1,
     $                R1( 1 ), R2( 1 ), W, 1, RW, B, INFO )
         CALL CHKXER( 'CGEESX', INFOT, NOUT, LERR, OK )
         NT = NT + 7
*
      ELSE IF( LSAMEN( 2, C2, 'BD' ) ) THEN
*
*        Test CGESVD
*
         SRNAMT = 'CGESVD'
         INFOT = 1
         CALL CGESVD( 'X', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW,
     $                INFO )
         CALL CHKXER( 'CGESVD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGESVD( 'N', 'X', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW,
     $                INFO )
         CALL CHKXER( 'CGESVD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGESVD( 'O', 'O', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW,
     $                INFO )
         CALL CHKXER( 'CGESVD', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGESVD( 'N', 'N', -1, 0, A, 1, S, U, 1, VT, 1, W, 1, RW,
     $                INFO )
         CALL CHKXER( 'CGESVD', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL CGESVD( 'N', 'N', 0, -1, A, 1, S, U, 1, VT, 1, W, 1, RW,
     $                INFO )
         CALL CHKXER( 'CGESVD', INFOT, NOUT, LERR, OK )
         INFOT = 6
         CALL CGESVD( 'N', 'N', 2, 1, A, 1, S, U, 1, VT, 1, W, 5, RW,
     $                INFO )
         CALL CHKXER( 'CGESVD', INFOT, NOUT, LERR, OK )
         INFOT = 9
         CALL CGESVD( 'A', 'N', 2, 1, A, 2, S, U, 1, VT, 1, W, 5, RW,
     $                INFO )
         CALL CHKXER( 'CGESVD', INFOT, NOUT, LERR, OK )
         INFOT = 11
         CALL CGESVD( 'N', 'A', 1, 2, A, 1, S, U, 1, VT, 1, W, 5, RW,
     $                INFO )
         CALL CHKXER( 'CGESVD', INFOT, NOUT, LERR, OK )
         NT = NT + 8
         IF( OK ) THEN
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ),
     $           NT
         ELSE
            WRITE( NOUT, FMT = 9998 )
         END IF
*
*        Test CGESDD
*
         SRNAMT = 'CGESDD'
         INFOT = 1
         CALL CGESDD( 'X', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, IW,
     $                INFO )
         CALL CHKXER( 'CGESDD', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL CGESDD( 'N', -1, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, IW,
     $                INFO )
         CALL CHKXER( 'CGESDD', INFOT, NOUT, LERR, OK )
         INFOT = 3
         CALL CGESDD( 'N', 0, -1, A, 1, S, U, 1, VT, 1, W, 1, RW, IW,
     $                INFO )
         CALL CHKXER( 'CGESDD', INFOT, NOUT, LERR, OK )
         INFOT = 5
         CALL CGESDD( 'N', 2, 1, A, 1, S, U, 1, VT, 1, W, 5, RW, IW,
     $                INFO )
         CALL CHKXER( 'CGESDD', INFOT, NOUT, LERR, OK )
         INFOT = 8
         CALL CGESDD( 'A', 2, 1, A, 2, S, U, 1, VT, 1, W, 5, RW, IW,
     $                INFO )
         CALL CHKXER( 'CGESDD', INFOT, NOUT, LERR, OK )
         INFOT = 10
         CALL CGESDD( 'A', 1, 2, A, 1, S, U, 1, VT, 1, W, 5, RW, IW,
     $                INFO )
         CALL CHKXER( 'CGESDD', INFOT, NOUT, LERR, OK )
         NT = NT - 2
         IF( OK ) THEN
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ),
     $           NT
         ELSE
            WRITE( NOUT, FMT = 9998 )
         END IF
      END IF
*
*     Print a summary line.
*
      IF( .NOT.LSAMEN( 2, C2, 'BD' ) ) THEN
         IF( OK ) THEN
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ),
     $           NT
         ELSE
            WRITE( NOUT, FMT = 9998 )
         END IF
      END IF
*
 9999 FORMAT( 1X, A, ' passed the tests of the error exits (', I3,
     $      ' tests done)' )
 9998 FORMAT( ' *** ', A, ' failed the tests of the error exits ***' )
      RETURN
*
*     End of CERRED
*
      END
