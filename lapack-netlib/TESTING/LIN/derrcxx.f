*> \brief \b DERRCXX
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DERRCXX( PATH, NUNIT )
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
*> DERRCXX tests the error exits for DERRCXX that does
*> CX decomposition.
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
*> \ingroup double_lin
*
*  =====================================================================
      SUBROUTINE DERRCXX( PATH, NUNIT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER(LEN=3)   PATH
      INTEGER            NUNIT
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NMAX
      PARAMETER          ( NMAX = 5 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J, K
      DOUBLE PRECISION   MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $                   NAN, ONE, ZERO
*     ..
*     .. Local Arrays ..
      INTEGER            DESEL_ROWS( NMAX ), SEL_DESEL_COLS( NMAX ),
     $                   IPIV( NMAX ), JPIV( NMAX ), IW( NMAX )
      DOUBLE PRECISION   A( NMAX, NMAX ), C( NMAX, NMAX ),
     $                   QRC( NMAX, NMAX ), X( NMAX, NMAX ),
     $                   TAU( NMAX ), W( NMAX )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, DGECXX
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER(LEN=32)  SRNAMT
      INTEGER            INFOT, NOUT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, SQRT
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
*
*     Set the variables to innocuous values.
*
      DO J = 1, NMAX
         DESEL_ROWS( J ) = 0
         SEL_DESEL_COLS( J ) = 0
         IPIV( J ) = 0
         JPIV( J ) = 0
         TAU( J ) = 1.D+0 / DBLE( J )
         W( J ) = 1.D+0 / DBLE( J ) 
         IW( J ) = -J     
         DO I = 1, NMAX
            A( I, J ) = 1.D+0 / DBLE( I+J )
            C( I, J ) = 1.D+0 / DBLE( I+J )
            QRC( I, J ) = 1.D+0 / DBLE( I+J )
            X( I, J ) = 1.D+0 / DBLE( I+J )
         END DO
      END DO
*
*     Create a NaN
*
      ONE = 1.0D+0 
      ZERO = 0.0D+0
      NAN = SQRT( -ONE )
*
      OK = .TRUE.
*
*     Error exits for CX decomposition
*
*     DGECXX
*
      SRNAMT = 'DGECXX'
*      
*     ====================== 
*     Test parameter FACT
*     ====================== 
      INFOT = 1
      CALL DGECXX( '/', 'A', 0, 0,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 1,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 1, QRC, 1,
     $             X, 1, W, 1, IW, 1, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*      
*     ====================== 
*     Test parameter USESD
*     ======================
*
      INFOT = 2
*      
      CALL DGECXX( 'P', '/', 0, 0,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 1,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 1, QRC, 1,
     $             X, 1, W, 1, IW, 1, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )      
*
*     ====================== 
*     Test parameter M
*     ======================
* 
      INFOT = 3
*      
      CALL DGECXX( 'P', 'A', -1, 0,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 1,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 1, QRC, 1,
     $             X, 1, W, 1, IW, 1, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     ======================= 
*     Test parameter N
*     ======================= 
*
      INFOT = 4
*      
      CALL DGECXX( 'P', 'A', 0, -1,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 1,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 1, QRC, 1,
     $             X, 1, W, 1, IW, 1, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     =======================      
*     Test parameter SEL_DESEL_COLS
*     =======================
*
*     NSEL (the number of preselected columns in SEL_DESEL_COLS
*     (element value = 1)) cannot be greater then MSUB.
*
      INFOT = 6
*      
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = 0
      SEL_DESEL_COLS( 1 ) = 1
      SEL_DESEL_COLS( 2 ) = 1
      CALL DGECXX( 'P', 'A', 1, 2,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 1,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 1, QRC, 1,
     $             X, 1, W, 1, IW, 1, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     ======================= 
*     Test parameter KMAXFREE
*     ======================= 
*
      INFOT = 7
*      
      CALL DGECXX( 'P', 'A', 0, 0,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             -1, ONE, ONE, A, 1,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 1, QRC, 1,
     $             X, 1, W, 1, IW, 1, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     =======================      
*     Test parameter ABSTOL
*     =======================
*
      INFOT = 8
*      
      CALL DGECXX( 'P', 'A', 0, 0,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, NAN, ONE, A, 1,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 1, QRC, 1,
     $             X, 1, W, 1, IW, 1, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )

*
*     =======================
*     Test parameter RELTOL
*     =======================
*
      INFOT = 9
*      
      CALL DGECXX( 'P', 'A', 0, 0,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, NAN, A, 1,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 1, QRC, 1,
     $             X, 1, W, 1, IW, 1, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     =======================
*     Test parameter LDA
*     =======================
*
      INFOT = 11
*
*     min(M,N) = 0
*
      CALL DGECXX( 'P', 'A', 0, 0,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 0,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 1, QRC, 1,
     $             X, 1, W, 1, IW, 1, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )   
*
      CALL DGECXX( 'P', 'A', 2, 2,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 1,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 2, QRC, 2,
     $             X, 2, W, 20, IW, 20, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )      
*
*     =======================
*     Test parameter LDC
*     =======================
*
      INFOT = 20
*
*     min(M,N) = 0
*
      CALL DGECXX( 'P', 'A', 0, 0,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 1,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 0, QRC, 1,
     $             X, 1, W, 20, IW, 20, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     FACT = 'P'
*
      CALL DGECXX( 'P', 'A', 2, 2,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 2,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 0, QRC, 2,
     $             X, 2, W, 20, IW, 20, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     FACT = 'C'
*               
      CALL DGECXX( 'C', 'A', 2, 2,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 2,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 1, QRC, 2,
     $             X, 2, W, 20, IW, 20, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     FACT = 'X'
*
      CALL DGECXX( 'X', 'A', 2, 2,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 2,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 1, QRC, 2,
     $             X, 2, W, 20, IW, 20, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     =======================
*     Test parameter LDQRC
*     =======================
*
*     QRC is used only when the matrix X is returned.
*
      INFOT = 22
*
*     min(M,N) = 0
*
      CALL DGECXX( 'P', 'A', 0, 0,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 1,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 1, QRC, 0,
     $             X, 1, W, 20, IW, 20, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )      
*
*     FACT = 'P'
*
      CALL DGECXX( 'P', 'A', 2, 2,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 2,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 2, QRC, 0,
     $             X, 1, W, 20, IW, 20, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     FACT = 'C'
*            
      CALL DGECXX( 'C', 'A', 2, 2,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 2,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 2, QRC, 0,
     $             X, 2, W, 20, IW, 20, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     FACT = 'X'
*
      CALL DGECXX( 'X', 'A', 2, 2,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 2,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 2, QRC, 1,
     $             X, 2, W, 20, IW, 20, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*      
*     =======================
*     Test parameter LDX
*     =======================
*
      INFOT = 24
*
*     min(M,N) = 0
*
      CALL DGECXX( 'P', 'A', 0, 0,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 1,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 1, QRC, 1,
     $             X, 0, W, 20, IW, 20, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     FACT = 'P'
*     
      CALL DGECXX( 'P', 'A', 2, 2,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 2,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 2, QRC, 2,
     $             X, 0, W, 20, IW, 20, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )      
*
*     FACT = 'C'
*            
      CALL DGECXX( 'C', 'A', 2, 2,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 2,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 2, QRC, 2,
     $             X, 0, W, 20, IW, 20, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     FACT = 'X'
*
      CALL DGECXX( 'X', 'A', 4, 2,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 4,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 4, QRC, 4,
     $             X, 3, W, 20, IW, 20, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     =======================
*     Test parameter LWORK
*     =======================
*
      INFOT = 26
*      
*     Test group 1. LWORK test for MIN(M,N) = 0, then LWKMIN  => 1
*     ==========================================
*   
      CALL DGECXX( 'X', 'A', 0, 0,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 1,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 1, QRC, 1,
     $             X, 1, W, 0, IW, 1, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     Test group 2. LWORK tests for USESD = 'N'.
*     ==========================================
*     if FACT = 'P',  LWKMIN = MAX(1, 3*N - 1)
*     if FACT = 'C',  LWKMIN = MAX(1, 3*N - 1)
*     if FACT = 'X',  LWKMIN = MAX(1, 3*N - 1, MINMN + N) = MAX(1, 3*N - 1)
*
      CALL DGECXX( 'P', 'N', 2, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 2,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 2, QRC, 2,
     $             X, 4, W, 10, IW, 20, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )  
*          
      CALL DGECXX( 'C', 'N', 2, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 2,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 2, QRC, 2,
     $             X, 4, W, 10, IW, 20, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*     
      CALL DGECXX( 'X', 'N', 2, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 2,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 2, QRC, 2,
     $             X, 4, W, 10, IW, 20, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )



*
*     Test group 3. LWORK tests for USESD = 'R'.
*     ==========================================
*     if FACT = 'P',  LWKMIN = MAX(1, 3*N - 1)
*     if FACT = 'C',  LWKMIN = MAX(1, 3*N - 1)
*     if FACT = 'X',  LWKMIN = MAX(1, 3*N - 1, MINMN + N) = MAX(1, 3*N - 1) 
*  
      DESEL_ROWS( 1 ) = -1
      DESEL_ROWS( 2 ) = -1
      DESEL_ROWS( 3 ) = 0
      DESEL_ROWS( 4 ) = 0
      DESEL_ROWS( 5 ) = -1 
*          
      CALL DGECXX( 'P', 'R', 5, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 10, IW, 20, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
      DESEL_ROWS( 1 ) = -1
      DESEL_ROWS( 2 ) = 0
      DESEL_ROWS( 3 ) = -1
      DESEL_ROWS( 4 ) = 0
      DESEL_ROWS( 5 ) = -1 
*               
      CALL DGECXX( 'C', 'R', 5, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 10, IW, 20, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*  
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = -1
      DESEL_ROWS( 3 ) = 0
      DESEL_ROWS( 4 ) = -1
      DESEL_ROWS( 5 ) = -1 
*          
      CALL DGECXX( 'X', 'R', 5, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 10, IW, 20, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     Test group 4. LWORK tests for USESD = 'C'.
*     ==========================================
*     (a) if FACT = 'P',  LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) )
*     (b) if FACT = 'C',  LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) )
*     (c) if FACT = 'X',  LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1), min(M,N)+N )
*
*     Parameter LWORK.
*     Case g4(a). USESD = 'C', if FACT = 'P', then LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) ). 
*     Test g4(a1). Set min(1,MINMNFREE == 1 ( i.e. enable (3*N_free-1) ) AND (3*N_free-1) is the largest component.
*          M = 4, N = 4, 
*          M_sub = M = 4, N_sub = N = 4,
*          N_sel = 0, M_free = M_sub - N_sel = 4,  N_free = N_sub - N_sel = 4,
*          MINMNFREE = min( M_free, N_free ) = min( 4, 4 ) = 4,
*          (3*N_free - 1) = 11
*          LWKMIN = (3*N_free - 1) = 11
*
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = 0
      SEL_DESEL_COLS( 3 ) = 0
      SEL_DESEL_COLS( 4 ) = 0      
*      
      CALL DGECXX( 'P', 'C', 4, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 4,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 4, QRC, 4,
     $             X, 4, W, 10, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     Parameter LWORK.
*     Case g4(a). USESD = 'C', if FACT = 'P', then LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) ). 
*     Test g4(a2). Set min(1,MINMNFREE) = 1 ( i.e. enable (3*N_free-1) ) AND N_sub is the largest component.
*          M = 4, N = 4, 
*          M_sub = M = 4, N_sub = N = 4,
*          N_sel = 3, M_free = M_sub - N_sel = 1,  N_free = N_sub - N_sel = 1,
*          MINMNFREE = min( M_free, N_free ) = min( 1, 1 ) = 1,
*          3*N_free - 1 = 2
*          LWKMIN = N_sub = 4
*
      SEL_DESEL_COLS( 1 ) = 1
      SEL_DESEL_COLS( 2 ) = 1
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 0       
*      
      CALL DGECXX( 'P', 'C', 4, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 4,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 4, QRC, 4,
     $             X, 4, W, 3, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*   
*     Parameter LWORK.  
*     Case g4(a). USESD = 'C', if FACT = 'P', then LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) ).
*     Test g4(a3). Set min(1,MINMNFREE) = 0 ( i.e. disable (3*N_free-1) ) AND (3*N_free-1) is the largest component.
*          M = 2, N = 4, 
*          M_sub = M = 2, N_sub = N = 4,
*          N_sel = 2, M_free = M_sub - N_sel = 0,  N_free = N_sub - N_sel = 2,
*          MINMNFREE = min( M_free, N_free ) = min( 0, 2 ) = 0,
*          (3*N_free - 1) = 5
*          LWKMIN = N_sub = 4
*
      SEL_DESEL_COLS( 1 ) = 1
      SEL_DESEL_COLS( 2 ) = 1
      SEL_DESEL_COLS( 3 ) = 0
      SEL_DESEL_COLS( 4 ) = 0       
*      
      CALL DGECXX( 'P', 'C', 2, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 2,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 2, QRC, 2,
     $             X, 4, W, 3, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     Parameter LWORK.     
*     Case g4(a). USESD = 'C', if FACT = 'P', then LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) ).
*     Test g4(a4). Set min(1,MINMNFREE) = 0 ( i.e. disable (3*N_free-1) ) AND N_sub is the largest component.
*          M = 3, N = 4, 
*          M_sub = M = 3, N_sub = N = 4,
*          N_sel = 3, M_free = M_sub - N_sel = 0,  N_free = N_sub - N_sel = 1,
*          MINMNFREE = min( M_free, N_free ) = min( 0, 1 ) = 0,
*          (3*N_free - 1) = 2
*          LWKMIN = N_sub = 4
*
      SEL_DESEL_COLS( 1 ) = 1
      SEL_DESEL_COLS( 2 ) = 1
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 0       
*      
      CALL DGECXX( 'P', 'C', 3, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 3,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 3, QRC, 3,
     $             X, 4, W, 3, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     Parameter LWORK.
*     Case g4(b). USESD = 'C', FACT = 'C', then LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) ).
*     Test g4(b1). min(1,MINMNFREE) = 1 ( i.e. enable (3*N_free-1) ) AND (3*N_free-1) is the largest component.
*          M = 4, N = 4, 
*          M_sub = M = 4, N_sub = N = 4,
*          N_sel = 0, M_free = M_sub - N_sel = 4,  N_free = N_sub - N_sel = 4,
*          MINMNFREE = min( M_free, N_free ) = min( 4, 4 ) = 4,
*          (3*N_free - 1) = 11
*          LWKMIN = (3*N_free - 1) = 11
*
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = 0
      SEL_DESEL_COLS( 3 ) = 0
      SEL_DESEL_COLS( 4 ) = 0      
*      
      CALL DGECXX( 'C', 'C', 4, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 4,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 4, QRC, 4,
     $             X, 4, W, 10, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )     
*     
*     Parameter LWORK.
*     Case g4(b). USESD = 'C', FACT = 'C', then LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) ).
*     Test g4(b2). Set min(1,MINMNFREE) = 1 ( i.e. enable (3*N_free-1) ) AND N_sub is the largest component.
*          M = 4, N = 4, 
*          M_sub = M = 4, N_sub = N = 4,
*          N_sel = 3, M_free = M_sub - N_sel = 1,  N_free = N_sub - N_sel = 1,
*          MINMNFREE = min( M_free, N_free ) = min( 1, 1 ) = 1,
*          3*N_free - 1 = 2
*          LWKMIN = N_sub = 4
*
      SEL_DESEL_COLS( 1 ) = 1
      SEL_DESEL_COLS( 2 ) = 1
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 0       
*      
      CALL DGECXX( 'C', 'C', 4, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 4,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 4, QRC, 4,
     $             X, 4, W, 3, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     Parameter LWORK.     
*     Case g4(b). USESD = 'C', FACT = 'C', then LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) ).
*     Test g4(b3). Set min(1,MINMNFREE) = 0 ( i.e. disable (3*N_free-1) ) AND (3*N_free-1) is the largest component.
*          M = 2, N = 4, 
*          M_sub = M = 2, N_sub = N = 4,
*          N_sel = 2, M_free = M_sub - N_sel = 0,  N_free = N_sub - N_sel = 2,
*          MINMNFREE = min( M_free, N_free ) = min( 0, 2 ) = 0,
*          (3*N_free - 1) = 5
*          LWKMIN = N_sub = 4
*
      SEL_DESEL_COLS( 1 ) = 1
      SEL_DESEL_COLS( 2 ) = 1
      SEL_DESEL_COLS( 3 ) = 0
      SEL_DESEL_COLS( 4 ) = 0       
*      
      CALL DGECXX( 'C', 'C', 2, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 2,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 2, QRC, 2,
     $             X, 4, W, 3, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*  
*     Parameter LWORK.   
*     Case g4(b). USESD = 'C', FACT = 'C', then LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) ).
*     Test g4(b4). Set min(1,MINMNFREE) = 0 ( i.e. disable (3*N_free-1) ) AND N_sub is the largest component.
*          M = 3, N = 4, 
*          M_sub = M = 3, N_sub = N = 4,
*          N_sel = 3, M_free = M_sub - N_sel = 0,  N_free = N_sub - N_sel = 1,
*          MINMNFREE = min( M_free, N_free ) = min( 0, 1 ) = 0,
*          (3*N_free - 1) = 2
*          LWKMIN = N_sub = 4
*
      SEL_DESEL_COLS( 1 ) = 1
      SEL_DESEL_COLS( 2 ) = 1
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 0       
*      
      CALL DGECXX( 'C', 'C', 3, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 3,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 3, QRC, 3,
     $             X, 4, W, 3, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*      
*     Parameter LWORK.      
*     Case g4(c). USESD = 'C', FACT = 'X', then LWKMIN = max( 1, min(M,N)+N, min(1,MINMNFREE)*(3*N_free-1) ).
*     Test g4(c1). Set min(1,MINMNFREE) = 1 ( i.e. enable (3*N_free-1) ) AND (3*N_free-1) is the largest component.
*          M = 4, N = 4, 
*          M_sub = M = 4, N_sub = N = 4,
*          N_sel = 0, M_free = M_sub - N_sel = 4,  N_free = N_sub - N_sel = 4,
*          MINMNFREE = min( M_free, N_free ) = min( 4, 4 ) = 4,
*          (3*N_free - 1) = 11
*          (min(M,N)+N) = 4 + 4 = 8
*          LWKMIN = (3*N_free - 1) = 11
*
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = 0
      SEL_DESEL_COLS( 3 ) = 0
      SEL_DESEL_COLS( 4 ) = 0      
*      
      CALL DGECXX( 'X', 'C', 4, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 4,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 4, QRC, 4,
     $             X, 4, W, 10, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
* 
*     Parameter LWORK.    
*     Case g4(c). USESD = 'C', FACT = 'X', then LWKMIN = max( 1, min(M,N)+N, min(1,MINMNFREE)*(3*N_free-1) ).
*     Test g4(c2). Set min(1,MINMNFREE) = 1 ( i.e. enable (3*N_free-1) ) AND (min(M,N)+N) is the largest component.
*          M = 4, N = 4, 
*          M_sub = M = 4, N_sub = N = 4,
*          N_sel = 3, M_free = M_sub - N_sel = 1,  N_free = N_sub - N_sel = 1,
*          MINMNFREE = min( M_free, N_free ) = min( 1, 1 ) = 1,
*          (3*N_free - 1) = 2
*          (min(M,N)+N) = 4 + 4 = 8
*          LWKMIN = (min(M,N)+N) = 8
*
      SEL_DESEL_COLS( 1 ) = 1
      SEL_DESEL_COLS( 2 ) = 1
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 0       
*      
      CALL DGECXX( 'X', 'C', 4, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 4,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 4, QRC, 4,
     $             X, 4, W, 7, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     Parameter LWORK.     
*     Case g4(c). USESD = 'C', FACT = 'X', then LWKMIN = max( 1, min(M,N)+N, min(1,MINMNFREE)*(3*N_free-1) ).
*     Test g4(c3). Set min(1,MINMNFREE) = 0 ( i.e. disable (3*N_free-1) ) AND (3*N_free-1) is the largest component.
*          M = 2, N = 5, 
*          M_sub = M = 2, N_sub = N = 5,
*          N_sel = 2, M_free = M_sub - N_sel = 0,  N_free = N_sub - N_sel = 3,
*          MINMNFREE = min( M_free, N_free ) = min( 0, 3 ) = 0,
*          (3*N_free - 1) = 8
*          (min(M,N)+N) = 2 + 5 = 7
*          LWKMIN = (min(M,N)+N) = 7
*
      SEL_DESEL_COLS( 1 ) = 1
      SEL_DESEL_COLS( 2 ) = 1
      SEL_DESEL_COLS( 3 ) = 0
      SEL_DESEL_COLS( 4 ) = 0
      SEL_DESEL_COLS( 5 ) = 0         
*
      CALL DGECXX( 'X', 'C', 2, 5,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 2,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 2, QRC, 2,
     $             X, 2, W, 6, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     Parameter LWORK.     
*     Case g4(c). USESD = 'C', FACT = 'X', then LWKMIN = max( 1, min(M,N)+N, min(1,MINMNFREE)*(3*N_free-1) ).
*     Test g4(c4). Set min(1,MINMNFREE) = 0 ( i.e. disable (3*N_free-1) ) AND (min(M,N)+N) is the largest component.
*          M = 3, N = 4, 
*          M_sub = M = 3, N_sub = N = 4,
*          N_sel = 3, M_free = M_sub - N_sel = 0,  N_free = N_sub - N_sel = 1,
*          MINMNFREE = min( M_free, N_free ) = min( 0, 1 ) = 0,
*          (3*N_free - 1) = 2
*          (min(M,N)+N) = 3+4 = 7
*          LWKMIN = (min(M,N)+N) = 7
*
      SEL_DESEL_COLS( 1 ) = 1
      SEL_DESEL_COLS( 2 ) = 1
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 0       
*      
      CALL DGECXX( 'X', 'C', 3, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 3,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 3, QRC, 3,
     $             X, 4, W, 6, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     Test group 5. LWORK tests for USESD = 'A'.
*     ==========================================
*     (a) if FACT = 'P',  LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) )
*     (b) if FACT = 'C',  LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) )
*     (c) if FACT = 'X',  LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1), min(M,N)+N )
*
*     Parameter LWORK.
*     Case g4(a). USESD = 'A', if FACT = 'P', then LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) ). 
*     Test g4(a1). Set min(1,MINMNFREE) = 1 ( i.e. enable (3*N_free-1) ) AND (3*N_free-1) is the largest component.
*          M = 5, N = 4, 
*          M_sub = 4, N_sub = N = 4,
*          N_sel = 0, M_free = M_sub - N_sel = 4,  N_free = N_sub - N_sel = 4,
*          MINMNFREE = min( M_free, N_free ) = min( 4, 4 ) = 4,
*          (3*N_free - 1) = 11
*          LWKMIN = (3*N_free - 1) = 11
*
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = -1
      DESEL_ROWS( 3 ) = 0
      DESEL_ROWS( 4 ) = 0
      DESEL_ROWS( 5 ) = 0
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = 0
      SEL_DESEL_COLS( 3 ) = 0
      SEL_DESEL_COLS( 4 ) = 0      
*      
      CALL DGECXX( 'P', 'A', 5, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 10, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     Parameter LWORK.
*     Case g4(a). USESD = 'A', if FACT = 'P', then LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) ). 
*     Test g4(a2). Set min(1,MINMNFREE) = 1 ( i.e. enable (3*N_free-1) ) AND N_sub is the largest component.
*          M = 5, N = 4, 
*          M_sub = 4, N_sub = N = 4,
*          N_sel = 3, M_free = M_sub - N_sel = 1,  N_free = N_sub - N_sel = 1,
*          MINMNFREE = min( M_free, N_free ) = min( 1, 1 ) = 1,
*          3*N_free - 1 = 2
*          LWKMIN = N_sub = 4
*
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = 0
      DESEL_ROWS( 3 ) = 0
      DESEL_ROWS( 4 ) = 0
      DESEL_ROWS( 5 ) = -1
      SEL_DESEL_COLS( 1 ) = 1
      SEL_DESEL_COLS( 2 ) = 1
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 0       
*      
      CALL DGECXX( 'P', 'A', 5, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 3, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*   
*     Parameter LWORK.  
*     Case g4(a). USESD = 'A', if FACT = 'P', then LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) ).
*     Test g4(a3). Set min(1,MINMNFREE) = 0 ( i.e. disable (3*N_free-1) ) AND (3*N_free-1) is the largest component.
*          M = 5, N = 4, 
*          M_sub = 2, N_sub = N = 4,
*          N_sel = 2, M_free = M_sub - N_sel = 0,  N_free = N_sub - N_sel = 2,
*          MINMNFREE = min( M_free, N_free ) = min( 0, 2 ) = 0,
*          (3*N_free - 1) = 5
*          LWKMIN = N_sub = 4

      DESEL_ROWS( 1 ) = -1
      DESEL_ROWS( 2 ) = 0
      DESEL_ROWS( 3 ) = -1
      DESEL_ROWS( 4 ) = 0
      DESEL_ROWS( 5 ) = -1      
      SEL_DESEL_COLS( 1 ) = 1
      SEL_DESEL_COLS( 2 ) = 1
      SEL_DESEL_COLS( 3 ) = 0
      SEL_DESEL_COLS( 4 ) = 0       
*      
      CALL DGECXX( 'P', 'A', 5, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 3, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     Parameter LWORK.     
*     Case g4(a). USESD = 'A', if FACT = 'P', then LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) ).
*     Test g4(a4). Set min(1,MINMNFREE) = 0 ( i.e. disable (3*N_free-1) ) AND N_sub is the largest component.
*          M = 5, N = 4, 
*          M_sub = M = 3, N_sub = N = 4,
*          N_sel = 3, M_free = M_sub - N_sel = 0,  N_free = N_sub - N_sel = 1,
*          MINMNFREE = min( M_free, N_free ) = min( 0, 1 ) = 0,
*          (3*N_free - 1) = 2
*          LWKMIN = N_sub = 4
*
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = 0
      DESEL_ROWS( 3 ) = -1
      DESEL_ROWS( 4 ) = -1
      DESEL_ROWS( 5 ) = 0 
      SEL_DESEL_COLS( 1 ) = 1
      SEL_DESEL_COLS( 2 ) = 1
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 0       
*      
      CALL DGECXX( 'P', 'A', 5, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 3, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     Parameter LWORK.
*     Case g4(b). USESD = 'A', FACT = 'C', then LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) ).
*     Test g4(b1). min(1,MINMNFREE) = 1 ( i.e. enable (3*N_free-1) ) AND (3*N_free-1) is the largest component.
*          M = 5, N = 4, 
*          M_sub = 4, N_sub = N = 4,
*          N_sel = 0, M_free = M_sub - N_sel = 4,  N_free = N_sub - N_sel = 4,
*          MINMNFREE = min( M_free, N_free ) = min( 4, 4 ) = 4,
*          (3*N_free - 1) = 11
*          LWKMIN = (3*N_free - 1) = 11
*
      DESEL_ROWS( 1 ) = -1
      DESEL_ROWS( 2 ) = 0
      DESEL_ROWS( 3 ) = 0
      DESEL_ROWS( 4 ) = 0
      DESEL_ROWS( 5 ) = 0
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = 0
      SEL_DESEL_COLS( 3 ) = 0
      SEL_DESEL_COLS( 4 ) = 0      
*      
      CALL DGECXX( 'C', 'A', 5, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 10, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )     
*     
*     Parameter LWORK.
*     Case g4(b). USESD = 'A', FACT = 'C', then LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) ).
*     Test g4(b2). Set min(1,MINMNFREE) = 1 ( i.e. enable (3*N_free-1) ) AND N_sub is the largest component.
*          M = 5, N = 4, 
*          M_sub = M = 4, N_sub = N = 4,
*          N_sel = 3, M_free = M_sub - N_sel = 1,  N_free = N_sub - N_sel = 1,
*          MINMNFREE = min( M_free, N_free ) = min( 1, 1 ) = 1,
*          3*N_free - 1 = 2
*          LWKMIN = N_sub = 4
*
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = 0
      DESEL_ROWS( 3 ) = -1
      DESEL_ROWS( 4 ) = 0
      DESEL_ROWS( 5 ) = 0      
      SEL_DESEL_COLS( 1 ) = 1
      SEL_DESEL_COLS( 2 ) = 1
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 0       
*      
      CALL DGECXX( 'C', 'A', 5, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 3, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     Parameter LWORK.     
*     Case g4(b). USESD = 'A', FACT = 'C', then LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) ).
*     Test g4(b3). Set min(1,MINMNFREE) = 0 ( i.e. disable (3*N_free-1) ) AND (3*N_free-1) is the largest component.
*          M = 5, N = 4, 
*          M_sub = 2, N_sub = N = 4,
*          N_sel = 2, M_free = M_sub - N_sel = 0,  N_free = N_sub - N_sel = 2,
*          MINMNFREE = min( M_free, N_free ) = min( 0, 2 ) = 0,
*          (3*N_free - 1) = 5
*          LWKMIN = N_sub = 4
*
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = -1
      DESEL_ROWS( 3 ) = -1
      DESEL_ROWS( 4 ) = -1
      DESEL_ROWS( 5 ) = 0
      SEL_DESEL_COLS( 1 ) = 1
      SEL_DESEL_COLS( 2 ) = 1
      SEL_DESEL_COLS( 3 ) = 0
      SEL_DESEL_COLS( 4 ) = 0       
*      
      CALL DGECXX( 'C', 'A', 5, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 3, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*  
*     Parameter LWORK.   
*     Case g4(b). USESD = 'A', FACT = 'C', then LWKMIN = max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) ).
*     Test g4(b4). Set min(1,MINMNFREE) = 0 ( i.e. disable (3*N_free-1) ) AND N_sub is the largest component.
*          M = 5, N = 4, 
*          M_sub = 3, N_sub = N = 4,
*          N_sel = 3, M_free = M_sub - N_sel = 0,  N_free = N_sub - N_sel = 1,
*          MINMNFREE = min( M_free, N_free ) = min( 0, 1 ) = 0,
*          (3*N_free - 1) = 2
*          LWKMIN = N_sub = 4
*
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = -1
      DESEL_ROWS( 3 ) = -1
      DESEL_ROWS( 4 ) = 0
      DESEL_ROWS( 4 ) = 0
      SEL_DESEL_COLS( 1 ) = 1
      SEL_DESEL_COLS( 2 ) = 1
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 0       
*      
      CALL DGECXX( 'C', 'A', 5, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 3, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*      
*     Parameter LWORK.      
*     Case g4(c). USESD = 'A', FACT = 'X', then LWKMIN = max( 1, min(M,N)+N, min(1,MINMNFREE)*(3*N_free-1) ).
*     Test g4(c1). Set min(1,MINMNFREE) = 1 ( i.e. enable (3*N_free-1) ) AND (3*N_free-1) is the largest component.
*          M = 5, N = 4, 
*          M_sub = 4, N_sub = N = 4,
*          N_sel = 0, M_free = M_sub - N_sel = 4,  N_free = N_sub - N_sel = 4,
*          MINMNFREE = min( M_free, N_free ) = min( 4, 4 ) = 4,
*          (3*N_free - 1) = 11
*          (min(M,N)+N) = 4 + 4 = 8
*          LWKMIN = (3*N_free - 1) = 11
*
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = 0
      DESEL_ROWS( 3 ) = 0
      DESEL_ROWS( 4 ) = -1
      DESEL_ROWS( 5 ) = 0
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = 0
      SEL_DESEL_COLS( 3 ) = 0
      SEL_DESEL_COLS( 4 ) = 0      
*      
      CALL DGECXX( 'X', 'A', 5, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 10, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
* 
*     Parameter LWORK.    
*     Case g4(c). USESD = 'A', FACT = 'X', then LWKMIN = max( 1, min(M,N)+N, min(1,MINMNFREE)*(3*N_free-1) ).
*     Test g4(c2). Set min(1,MINMNFREE) = 1 ( i.e. enable (3*N_free-1) ) AND (min(M,N)+N) is the largest component.
*          M = 5, N = 4, 
*          M_sub = 4, N_sub = N = 4,
*          N_sel = 3, M_free = M_sub - N_sel = 1,  N_free = N_sub - N_sel = 1,
*          MINMNFREE = min( M_free, N_free ) = min( 1, 1 ) = 1,
*          (3*N_free - 1) = 2
*          (min(M,N)+N) = 4 + 4 = 8
*          LWKMIN = (min(M,N)+N) = 8
*
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = 0
      DESEL_ROWS( 3 ) = 0
      DESEL_ROWS( 4 ) = -1
      DESEL_ROWS( 5 ) = 0
      SEL_DESEL_COLS( 1 ) = 1
      SEL_DESEL_COLS( 2 ) = 1
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 0       
*      
      CALL DGECXX( 'X', 'A', 5, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 7, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     Parameter LWORK.     
*     Case g4(c). USESD = 'A', FACT = 'X', then LWKMIN = max( 1, min(M,N)+N, min(1,MINMNFREE)*(3*N_free-1) ).
*     Test g4(c3). Set min(1,MINMNFREE) = 0 ( i.e. disable (3*N_free-1) ) AND (3*N_free-1) is the largest component.
*          M = 5, N = 5, 
*          M_sub = 2, N_sub = N = 5,
*          N_sel = 2, M_free = M_sub - N_sel = 0,  N_free = N_sub - N_sel = 3,
*          MINMNFREE = min( M_free, N_free ) = min( 0, 3 ) = 0,
*          (3*N_free - 1) = 8
*          (min(M,N)+N) = 2 + 5 = 7
*          LWKMIN = (min(M,N)+N) = 7
*
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = -1
      DESEL_ROWS( 3 ) = -1
      DESEL_ROWS( 4 ) = -1
      DESEL_ROWS( 5 ) = 0
      SEL_DESEL_COLS( 1 ) = 1
      SEL_DESEL_COLS( 2 ) = 1
      SEL_DESEL_COLS( 3 ) = 0
      SEL_DESEL_COLS( 4 ) = 0
      SEL_DESEL_COLS( 5 ) = 0         
*
      CALL DGECXX( 'X', 'A', 5, 5,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 6, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     Parameter LWORK.     
*     Case g4(c). USESD = 'A', FACT = 'X', then LWKMIN = max( 1, min(M,N)+N, min(1,MINMNFREE)*(3*N_free-1) ).
*     Test g4(c4). Set min(1,MINMNFREE) = 0 ( i.e. disable (3*N_free-1) ) AND (min(M,N)+N) is the largest component.
*          M = 5, N = 4, 
*          M_sub = 3, N_sub = N = 4,
*          N_sel = 3, M_free = M_sub - N_sel = 0,  N_free = N_sub - N_sel = 1,
*          MINMNFREE = min( M_free, N_free ) = min( 0, 1 ) = 0,
*          (3*N_free - 1) = 2
*          (min(M,N)+N) = 3+4 = 7
*          LWKMIN = (min(M,N)+N) = 7
*
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = -1
      DESEL_ROWS( 3 ) = -1
      DESEL_ROWS( 4 ) = 0
      DESEL_ROWS( 5 ) = 0
      SEL_DESEL_COLS( 1 ) = 1
      SEL_DESEL_COLS( 2 ) = 1
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 0       
*      
      CALL DGECXX( 'X', 'A', 5, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 6, IW, 10, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     =======================
*     Test parameter LIWORK
*     =======================
*
      INFOT = 28
*      
*     Test group 1. LIWORK test for MIN(M,N) = 0, then LWKMIN  => 1
*     ==========================================
* 
      CALL DGECXX( 'X', 'A', 0, 0,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 1,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 1, QRC, 1,
     $             X, 1, W, 1, IW, 0, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     Test group 2. LIWORK tests for USESD = 'N'
*     ==========================================
*     if FACT = 'P',  LIWKMIN = MAX(1, N-1)
*     if FACT = 'C',  LIWKMIN = MAX(1, 2*N)
*     if FACT = 'X',  LIWKMIN = MAX(1, 2*N)
*
      CALL DGECXX( 'P', 'N', 2, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 2,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 2, QRC, 2,
     $             X, 4, W, 11, IW, 2, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
      CALL DGECXX( 'C', 'N', 2, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 2,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 2, QRC, 2,
     $             X, 4, W, 11, IW, 7, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
      CALL DGECXX( 'X', 'N', 2, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 2,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 2, QRC, 2,
     $             X, 4, W, 11, IW, 7, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     Test group 3. LIWORK tests for USESD = 'R'
*     ==========================================
*     if FACT = 'P',  LIWKMIN = MAX(1, N-1)
*     if FACT = 'C',  LIWKMIN = MAX(1, 2*N)
*     if FACT = 'X',  LIWKMIN = MAX(1, 2*N)
*
      DESEL_ROWS( 1 ) = -1
      DESEL_ROWS( 2 ) = 0
      DESEL_ROWS( 3 ) = -1
      DESEL_ROWS( 4 ) = 0
      DESEL_ROWS( 5 ) = -1  
*      
      CALL DGECXX( 'P', 'R', 5, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 11, IW, 2, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = -1
      DESEL_ROWS( 3 ) = 0
      DESEL_ROWS( 4 ) = -1
      DESEL_ROWS( 5 ) = -1  
*       
      CALL DGECXX( 'C', 'R', 5, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 11, IW, 7, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = -1
      DESEL_ROWS( 3 ) = 0
      DESEL_ROWS( 4 ) = -1
      DESEL_ROWS( 5 ) = -1  
*       
      CALL DGECXX( 'X', 'R', 5, 4,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 11, IW, 7, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     Test group 4. LIWORK tests for USESD = 'C'.
*     ==========================================
*     (a) if FACT = 'P',  LIWKMIN = max( 1, (N_free-1) + min(1,N_sel)*N_free )
*     (b) if FACT = 'C',  LIWKMIN = max( 1, 2*N )
*     (c) if FACT = 'X',  LIWKMIN = max( 1, 2*N )
*
*     Parameter LIWORK.
*     Case g4(a). USESD = 'C', if FACT = 'P', then LIWKMIN = max( 1, (N_free-1) + min(1,N_sel)*N_free ). 
*     Test g4(a1). Set min(1,N_sel) = 0 (i.e. disable N_free term).
*          M = 5, N = 5, 
*          M_sub = M = 5, N_sub = 4,
*          N_sel = 0, M_free = M_sub - N_sel = 5,  N_free = N_sub - N_sel = 4,
*          min(1,N_sel) = 0
*          LIWKMIN = (N_free-1) = 4 - 1 = 3
*
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = -1
      SEL_DESEL_COLS( 3 ) = 0
      SEL_DESEL_COLS( 4 ) = 0
      SEL_DESEL_COLS( 5 ) = 0       
*      
      CALL DGECXX( 'P', 'C', 5, 5,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 11, IW, 2, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*      
*     Parameter LIWORK.
*     Case g4(a). USESD = 'C', if FACT = 'P', then LIWKMIN = max( 1, (N_free-1) + min(1,N_sel)*N_free ). 
*     Test g4(a2). Set min(1,N_sel) = 1 (i.e. enable N_free term).
*          M = 5, N = 5, 
*          M_sub = M = 5, N_sub = 4,
*          N_sel = 1, M_free = M_sub - N_sel = 4,  N_free = N_sub - N_sel = 3,
*          min(1,N_sel) = 1
*          LIWKMIN = (N_free-1) + N_free = 3 - 1 + 3 = 5
*
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = -1
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 0
      SEL_DESEL_COLS( 5 ) = 0       
*      
      CALL DGECXX( 'P', 'C', 5, 5,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 11, IW, 4, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*      
*     Parameter LIWORK.
*     Case g4(b). USESD = 'C', if FACT = 'C', then LIWKMIN = max( 1, 2*N )
*     Test g4(b1). (N_free-1) + min(1,N_sel)*N_free. 
*          Set min(1,N_sel) = 0 (i.e. disable N_free term).
*          M = 5, N = 5, 
*          M_sub = M = 5, N_sub = 4,
*          N_sel = 0, M_free = M_sub - N_sel = 5,  N_free = N_sub - N_sel = 4,
*          min(1,N_sel) = 0
*          LIWKMIN = N = 5
*
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = -1
      SEL_DESEL_COLS( 3 ) = 0
      SEL_DESEL_COLS( 4 ) = 0
      SEL_DESEL_COLS( 5 ) = 0       
*      
      CALL DGECXX( 'C', 'C', 5, 5,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 11, IW, 9, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*      
*     Parameter LIWORK.
*     Case g4(b). USESD = 'C', if FACT = 'C', then LIWKMIN = max( 1, 2*N )
*     Test g4(b2). (N_free-1) + min(1,N_sel)*N_free
*          Set min(1,N_sel) = 1 (i.e. enable N_free term) AND N is the largest component.
*          M = 5, N = 5, 
*          M_sub = M = 5, N_sub = 4,
*          N_sel = 2, M_free = M_sub - N_sel = 3,  N_free = N_sub - N_sel = 2,
*          min(1,N_sel) = 1
*          (N_free - 1) + N_free = (2-1) + 2 = 3
*          LIWKMIN = N = 5
*
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = -1
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 1
      SEL_DESEL_COLS( 5 ) = 0       
*      
      CALL DGECXX( 'C', 'C', 5, 5,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 11, IW, 9, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )      
*      
*     Parameter LIWORK.
*     Case g4(b). USESD = 'C', if FACT = 'C', then LIWKMIN = max( 1, 2*N )
*     Test g4(b3). (N_free-1) + min(1,N_sel)*N_free
*          Set min(1,N_sel) = 1 (i.e. enable N_free term) AND ((N_free - 1) + N_free) is the largest component.
*          M = 5, N = 5, 
*          M_sub = M = 5, N_sub = N = 5,
*          N_sel = 1, M_free = M_sub - N_sel = 4,  N_free = N_sub - N_sel = 4,
*          min(1,N_sel) = 1
*          (N_free - 1) + N_free = (4-1) + 4 = 7
*          LIWKMIN = ((N_free - 1) + N_free) = 7
*
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = 0
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 0
      SEL_DESEL_COLS( 5 ) = 0       
*      
      CALL DGECXX( 'C', 'C', 5, 5,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 11, IW, 9, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*      
*     Parameter LIWORK.
*     Case g4(c). USESD = 'C', if FACT = 'X', then LIWKMIN = max( 1, 2*N )
*     Test g4(c1). (N_free-1) + min(1,N_sel)*N_free
*          Set min(1,N_sel) = 0 (i.e. disable N_free term).
*          M = 5, N = 5, 
*          M_sub = M = 5, N_sub = 4,
*          N_sel = 0, M_free = M_sub - N_sel = 5,  N_free = N_sub - N_sel = 4,
*          min(1,N_sel) = 0
*          LIWKMIN = N = 5
*
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = -1
      SEL_DESEL_COLS( 3 ) = 0
      SEL_DESEL_COLS( 4 ) = 0
      SEL_DESEL_COLS( 5 ) = 0       
*      
      CALL DGECXX( 'X', 'C', 5, 5,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 11, IW, 9, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*      
*     Parameter LIWORK.
*     Case g4(c). USESD = 'C', if FACT = 'X', then LIWKMIN = max( 1, 2*N )
*     Test g4(c2). (N_free-1) + min(1,N_sel)*N_free.
*        Set min(1,N_sel) = 1 (i.e. enable N_free term) AND N is the largest component.
*          M = 5, N = 5, 
*          M_sub = M = 5, N_sub = 4,
*          N_sel = 2, M_free = M_sub - N_sel = 3,  N_free = N_sub - N_sel = 2,
*          min(1,N_sel) = 1
*          (N_free - 1) + N_free = (2-1) + 2 = 3
*          LIWKMIN = N = 5`
*
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = -1
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 1
      SEL_DESEL_COLS( 5 ) = 0       
*      
      CALL DGECXX( 'X', 'C', 5, 5,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 11, IW, 9, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )      
*      
*     Parameter LIWORK.
*     Case g4(c). USESD = 'C', if FACT = 'X', then LIWKMIN = max( 1, 2*N )
*     Test g4(c3). (N_free-1) + min(1,N_sel)*N_free
*         Set min(1,N_sel) = 1 (i.e. enable N_free term) AND ((N_free - 1) + N_free) is the largest component.
*          M = 5, N = 5, 
*          M_sub = M = 5, N_sub = N = 5,
*          N_sel = 1, M_free = M_sub - N_sel = 4,  N_free = N_sub - N_sel = 4,
*          min(1,N_sel) = 1
*          (N_free - 1) + N_free = (4-1) + 4 = 7
*          LIWKMIN = ((N_free - 1) + N_free) = 7
*
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = 0
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 0
      SEL_DESEL_COLS( 5 ) = 0       
*      
      CALL DGECXX( 'X', 'C', 5, 5,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 11, IW, 9, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     Test group 5. LIWORK tests for USESD = 'A'.
*     ==========================================
*     (a) if FACT = 'P',  LIWKMIN = max( 1, (N_free-1) + min(1,N_sel)*N_free )
*     (b) if FACT = 'C',  LIWKMIN = max( 1, (N_free-1) + min(1,N_sel)*N_free, N )
*     (c) if FACT = 'X',  LIWKMIN = max( 1, (N_free-1) + min(1,N_sel)*N_free, N )
*
*     Parameter LIWORK.
*     Case g5(a). USESD = 'A', if FACT = 'P', then LIWKMIN = max( 1, (N_free-1) + min(1,N_sel)*N_free ). 
*     Test g5(a1). Set min(1,N_sel) = 0 (i.e. disable N_free term).
*          M = 5, N = 5, 
*          M_sub = 4, N_sub = 4,
*          N_sel = 0, M_free = M_sub - N_sel = 4,  N_free = N_sub - N_sel = 4,
*          min(1,N_sel) = 0
*          LIWKMIN = (N_free-1) = 4 - 1 = 3
*
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = -1
      DESEL_ROWS( 3 ) = 0
      DESEL_ROWS( 4 ) = 0
      DESEL_ROWS( 5 ) = 0
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = -1
      SEL_DESEL_COLS( 3 ) = 0
      SEL_DESEL_COLS( 4 ) = 0
      SEL_DESEL_COLS( 5 ) = 0       
*      
      CALL DGECXX( 'P', 'A', 5, 5,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 11, IW, 2, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*      
*     Parameter LIWORK.
*     Case g5(a). USESD = 'A', if FACT = 'P', then LIWKMIN = max( 1, (N_free-1) + min(1,N_sel)*N_free ). 
*     Test g5(a2). Set min(1,N_sel) = 1 (i.e. enable N_free term).
*          M = 5, N = 5, 
*          M_sub = 4, N_sub = 4,
*          N_sel = 1, M_free = M_sub - N_sel = 3,  N_free = N_sub - N_sel = 3,
*          min(1,N_sel) = 1
*          LIWKMIN = (N_free-1) + N_free = 3 - 1 + 3 = 5
*
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = -1
      DESEL_ROWS( 3 ) = 0
      DESEL_ROWS( 4 ) = 0
      DESEL_ROWS( 5 ) = 0
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = -1
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 0
      SEL_DESEL_COLS( 5 ) = 0       
*      
      CALL DGECXX( 'P', 'A', 5, 5,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 11, IW, 4, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*      
*     Parameter LIWORK.
*     Case g5(b). USESD = 'A', if FACT = 'C', then LIWKMIN = max( 1, 2*N )
*     Test g5(b1). (N_free-1) + min(1,N_sel)*N_free
*          Set min(1,N_sel) = 0 (i.e. disable N_free term).
*          M = 5, N = 5, 
*          M_sub = 4, N_sub = 4,
*          N_sel = 0, M_free = M_sub - N_sel = 4,  N_free = N_sub - N_sel = 4,
*          min(1,N_sel) = 0
*          LIWKMIN = N = 5
*
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = -1
      DESEL_ROWS( 3 ) = 0
      DESEL_ROWS( 4 ) = 0
      DESEL_ROWS( 5 ) = 0
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = -1
      SEL_DESEL_COLS( 3 ) = 0
      SEL_DESEL_COLS( 4 ) = 0
      SEL_DESEL_COLS( 5 ) = 0       
*      
      CALL DGECXX( 'C', 'A', 5, 5,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 11, IW, 9, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*      
*     Parameter LIWORK.
*     Case g5(b). USESD = 'A', if FACT = 'C', then LIWKMIN = max( 1, 2*N )
*     Test g5(b2). (N_free-1) + min(1,N_sel)*N_free
*          Set min(1,N_sel) = 1 (i.e. enable N_free term) AND N is the largest component.
*          M = 5, N = 5, 
*          M_sub = 4, N_sub = 4,
*          N_sel = 2, M_free = M_sub - N_sel = 2,  N_free = N_sub - N_sel = 2,
*          min(1,N_sel) = 1
*          (N_free - 1) + N_free = (2-1) + 2 = 3
*          LIWKMIN = N = 5
*
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = -1
      DESEL_ROWS( 3 ) = 0
      DESEL_ROWS( 4 ) = 0
      DESEL_ROWS( 5 ) = 0
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = -1
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 1
      SEL_DESEL_COLS( 5 ) = 0       
*      
      CALL DGECXX( 'C', 'A', 5, 5,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 11, IW, 9, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )      
*      
*     Parameter LIWORK.
*     Case g5(b). USESD = 'A', if FACT = 'C', then LIWKMIN = max( 1, 2*N )
*     Test g5(b3). (N_free-1) + min(1,N_sel)*N_free
*          Set min(1,N_sel) = 1 (i.e. enable N_free term) AND ((N_free - 1) + N_free) is the largest component.
*          M = 5, N = 5, 
*          M_sub = 4, N_sub = N = 5,
*          N_sel = 1, M_free = M_sub - N_sel = 3,  N_free = N_sub - N_sel = 4,
*          min(1,N_sel) = 1
*          (N_free - 1) + N_free = (4-1) + 4 = 7
*          LIWKMIN = ((N_free - 1) + N_free) = 7
*
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = -1
      DESEL_ROWS( 3 ) = 0
      DESEL_ROWS( 4 ) = 0
      DESEL_ROWS( 5 ) = 0
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = 0
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 0
      SEL_DESEL_COLS( 5 ) = 0       
*      
      CALL DGECXX( 'C', 'A', 5, 5,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 11, IW, 9, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*      
*     Parameter LIWORK.
*     Case g5(c). USESD = 'A', if FACT = 'X', then LIWKMIN = max( 1, 2*N )
*     Test g5(c1). (N_free-1) + min(1,N_sel)*N_free
*         Set min(1,N_sel) = 0 (i.e. disable N_free term).
*          M = 5, N = 5, 
*          M_sub = 4, N_sub = 4,
*          N_sel = 0, M_free = M_sub - N_sel = 4,  N_free = N_sub - N_sel = 4,
*          min(1,N_sel) = 0
*          LIWKMIN = N = 5
*
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = -1
      DESEL_ROWS( 3 ) = 0
      DESEL_ROWS( 4 ) = 0
      DESEL_ROWS( 5 ) = 0
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = -1
      SEL_DESEL_COLS( 3 ) = 0
      SEL_DESEL_COLS( 4 ) = 0
      SEL_DESEL_COLS( 5 ) = 0       
*      
      CALL DGECXX( 'X', 'A', 5, 5,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 11, IW, 4, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*      
*     Parameter LIWORK.
*     Case g5(c). USESD = 'A', if FACT = 'X', then LIWKMIN = max( 1, 2*N )
*     Test g5(c2). (N_free-1) + min(1,N_sel)*N_free 
*         Set min(1,N_sel) = 1 (i.e. enable N_free term) AND N is the largest component.
*          M = 5, N = 5, 
*          M_sub = 4, N_sub = 4,
*          N_sel = 2, M_free = M_sub - N_sel = 2,  N_free = N_sub - N_sel = 2,
*          min(1,N_sel) = 1
*          (N_free - 1) + N_free = (2-1) + 2 = 3
*          LIWKMIN = N = 5
*
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = -1
      DESEL_ROWS( 3 ) = 0
      DESEL_ROWS( 4 ) = 0
      DESEL_ROWS( 5 ) = 0
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = -1
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 1
      SEL_DESEL_COLS( 5 ) = 0       
*      
      CALL DGECXX( 'X', 'A', 5, 5,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 11, IW, 4, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )      
*      
*     Parameter LIWORK.
*     Case g5(c). USESD = 'A', if FACT = 'X', then LIWKMIN = max( 1, 2*N )
*     Test g5(c3).  (N_free-1) + min(1,N_sel)*N_free
*          Set min(1,N_sel) = 1 (i.e. enable N_free term) AND ((N_free - 1) + N_free) is the largest component.
*          M = 5, N = 5, 
*          M_sub = 4, N_sub = N = 5,
*          N_sel = 1, M_free = M_sub - N_sel = 3,  N_free = N_sub - N_sel = 4,
*          min(1,N_sel) = 1
*          (N_free - 1) + N_free = (4-1) + 4 = 7
*          LIWKMIN = ((N_free - 1) + N_free) = 7
*
      DESEL_ROWS( 1 ) = 0
      DESEL_ROWS( 2 ) = 0
      DESEL_ROWS( 3 ) = 0
      DESEL_ROWS( 4 ) = 0
      DESEL_ROWS( 5 ) = 0
      SEL_DESEL_COLS( 1 ) = 0
      SEL_DESEL_COLS( 2 ) = 0
      SEL_DESEL_COLS( 3 ) = 1
      SEL_DESEL_COLS( 4 ) = 0
      SEL_DESEL_COLS( 5 ) = 0       
*      
      CALL DGECXX( 'X', 'A', 5, 5,
     $             DESEL_ROWS, SEL_DESEL_COLS,
     $             0, ONE, ONE, A, 5,
     $             K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $             IPIV, JPIV, TAU, C, 5, QRC, 5,
     $             X, 5, W, 11, IW, 9, INFO )
      CALL CHKXER( 'DGECXX', INFOT, NOUT, LERR, OK )
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of DERRCXX
*
      END
