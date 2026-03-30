*> \brief \b CLARF1F applies an elementary reflector to a general rectangular
*              matrix assuming v(1) = 1.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CLARF1F + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarf1f.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarf1f.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarf1f.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CLARF1F( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*       .. Scalar Arguments ..
*       CHARACTER          SIDE
*       INTEGER            INCV, LDC, M, N
*       COMPLEX            TAU
*       ..
*       .. Array Arguments ..
*       COMPLEX            C( LDC, * ), V( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CLARF1F applies a complex elementary reflector H to a complex m by n matrix
*> C, from either the left or the right. H is represented in the form
*>
*>       H = I - tau * v * v**H
*>
*> where tau is a complex scalar and v is a complex vector assuming v(1) = 1.
*>
*> If tau = 0, then H is taken to be the unit matrix.
*>
*> To apply H**H (the conjugate transpose of H), supply conjg(tau) instead
*> tau.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': form  H * C
*>          = 'R': form  C * H
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix C.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is COMPLEX array, dimension
*>                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
*>                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
*>          The vector v in the representation of H. V is not used if
*>          TAU = 0.
*> \endverbatim
*>
*> \param[in] INCV
*> \verbatim
*>          INCV is INTEGER
*>          The increment between elements of v. INCV <> 0.
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is COMPLEX
*>          The value tau in the representation of H.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX array, dimension (LDC,N)
*>          On entry, the m by n matrix C.
*>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
*>          or C * H if SIDE = 'R'.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension
*>                         (N) if SIDE = 'L'
*>                      or (M) if SIDE = 'R'
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
*> \ingroup larf1f
*
*  =====================================================================
      SUBROUTINE CLARF1F( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      COMPLEX            TAU
*     ..
*     .. Array Arguments ..
      COMPLEX            C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),
     $                   ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, LASTV, LASTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           CAXPY, CGEMV, CGER, CSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILACLR, ILACLC
      EXTERNAL           LSAME, ILACLR, ILACLC
*     ..
*     .. Executable Statements ..
*
      APPLYLEFT = LSAME( SIDE, 'L' )
      LASTV = 1
      LASTC = 0
      IF( TAU.NE.ZERO ) THEN
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V up to V(1).
         IF( APPLYLEFT ) THEN
            LASTV = M
         ELSE
            LASTV = N
         END IF
         IF( INCV.GT.0 ) THEN
            I = 1 + (LASTV-1) * INCV
         ELSE
            I = 1
         END IF
!     Look for the last non-zero row in V.
         DO WHILE( LASTV.GT.1 .AND. V( I ).EQ.ZERO )
            LASTV = LASTV - 1
            I = I - INCV
         END DO
         IF( APPLYLEFT ) THEN
!     Scan for the last non-zero column in C(1:lastv,:).
            LASTC = ILACLC(LASTV, N, C, LDC)
         ELSE
!     Scan for the last non-zero row in C(:,1:lastv).
            LASTC = ILACLR(M, LASTV, C, LDC)
         END IF
      END IF
      IF( LASTC.EQ.0 ) THEN
         RETURN
      END IF
      IF( APPLYLEFT ) THEN
*
*        Form  H * C
*
         IF( LASTV.EQ.1 ) THEN
*
*           C(1,1:lastc) := ( 1 - tau ) * C(1,1:lastc)
*
            CALL CSCAL( LASTC, ONE - TAU, C, LDC )
         ELSE
*
*        w(1:lastc,1) := C(2:lastv,1:lastc)**H * v(2:lastv,1)
*
         CALL CGEMV( 'Conjugate transpose', LASTV - 1, LASTC, ONE,
     $               C( 2, 1 ), LDC, V( 1 + INCV ), INCV, ZERO,
     $               WORK, 1 )
*
*        w(1:lastc,1) += v(1,1) * C(1,1:lastc)**H
*
         DO I = 1, LASTC
            WORK( I ) = WORK( I ) + CONJG( C( 1, I ) )
         END DO
*
*        C(1, 1:lastc) += - tau * v(1,1) * w(1:lastc,1)**H
*
         DO I = 1, LASTC
            C( 1, I ) = C( 1, I ) - TAU * CONJG( WORK( I ) )
         END DO
*
*        C(2:lastv,1:lastc) += - tau * v(2:lastv,1) * w(1:lastc,1)**H
*
         CALL CGERC( LASTV - 1, LASTC, -TAU, V( 1 + INCV ), INCV,
     $               WORK, 1, C( 2, 1 ), LDC )
            END IF
      ELSE
*
*        Form  C * H
*
         IF( LASTV.EQ.1 ) THEN
*
*           C(1:lastc,1) := ( 1 - tau ) * C(1:lastc,1)
*
            CALL CSCAL( LASTC, ONE - TAU, C, 1 )
         ELSE
*
*           w(1:lastc,1) := C(1:lastc,2:lastv) * v(2:lastv,1)
*
            CALL CGEMV( 'No transpose', LASTC, LASTV - 1, ONE,
     $                  C( 1, 2 ), LDC, V( 1 + INCV ), INCV, ZERO,
     $                  WORK, 1 )
*
*           w(1:lastc,1) += v(1,1) * C(1:lastc,1)
*
            CALL CAXPY( LASTC, ONE, C, 1, WORK, 1 )
*
*           C(1:lastc,1) += - tau * v(1,1) * w(1:lastc,1)
*
            CALL CAXPY( LASTC, -TAU, WORK, 1, C, 1 )
*
*           C(1:lastc,2:lastv) += - tau * w(1:lastc,1) * v(2:lastv)**H
*
            CALL CGERC( LASTC, LASTV - 1, -TAU, WORK, 1,
     $                  V( 1 + INCV ), INCV, C( 1, 2 ), LDC )
         END IF
      END IF
      RETURN
*
*     End of CLARF1F
*
      END
