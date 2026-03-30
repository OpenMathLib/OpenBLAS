*> \brief \b ZLARF1L applies an elementary reflector to a general rectangular
*              matrix assuming v(lastv) = 1, where lastv is the last non-zero
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download ZLARF1L + dependencies
*> <a
*href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarf1l.f">
*> [TGZ]</a>
*> <a
*href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarf1l.f">
*> [ZIP]</a>
*> <a
*href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarf1l.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLARF1L( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*       .. Scalar Arguments ..
*       CHARACTER          SIDE
*       INTEGER            INCV, LDC, M, N
*       COMPLEX*16         TAU
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         C( LDC, * ), V( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLARF1L applies a complex elementary reflector H to a complex m by n matrix
*> C, from either the left or the right. H is represented in the form
*>
*>       H = I - tau * v * v**H
*>
*> where tau is a real scalar and v is a real vector assuming v(lastv) = 1,
*> where lastv is the last non-zero element.
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
*>          V is COMPLEX*16 array, dimension
*>                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
*>                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
*>          The vector v in the representation of H. V is not used if
*>          TAU = 0.
*> \endverbatim
*>
*> \param[in] INCV
*> \verbatim
*>          INCV is INTEGER
*>          The increment between elements of v. INCV > 0.
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is COMPLEX*16
*>          The value tau in the representation of H.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX*16 array, dimension (LDC,N)
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
*>          WORK is COMPLEX*16 array, dimension
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
      SUBROUTINE ZLARF1L( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      COMPLEX*16         TAU
*     ..
*     .. Array Arguments ..
      COMPLEX*16         C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, J, LASTV, LASTC, FIRSTV
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZAXPY, ZGEMV, ZGERC, ZSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAZLR, ILAZLC
      EXTERNAL           LSAME, ILAZLR, ILAZLC
*     ..
*     .. Executable Statements ..
*
      APPLYLEFT = LSAME( SIDE, 'L' )
      FIRSTV = 1
      LASTC = 0
      IF( TAU.NE.ZERO ) THEN
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V up to V(1).
         IF( APPLYLEFT ) THEN
            LASTV = M
         ELSE
            LASTV = N
         END IF
         I = 1
!     Look for the last non-zero row in V.
         DO WHILE( LASTV.GT.FIRSTV .AND. V( I ).EQ.ZERO )
            FIRSTV = FIRSTV + 1
            I = I + INCV
         END DO
         IF( APPLYLEFT ) THEN
!     Scan for the last non-zero column in C(1:lastv,:).
            LASTC = ILAZLC(LASTV, N, C, LDC)
         ELSE
!     Scan for the last non-zero row in C(:,1:lastv).
            LASTC = ILAZLR(M, LASTV, C, LDC)
         END IF
      END IF
      IF( LASTC.EQ.0 ) THEN
         RETURN
      END IF
      IF( APPLYLEFT ) THEN
*
*        Form  H * C
*
         IF( LASTV.EQ.FIRSTV ) THEN        
*
*           C(lastv,1:lastc) := ( 1 - tau ) * C(lastv,1:lastc)
*
            CALL ZSCAL( LASTC, ONE - TAU, C( LASTV, 1 ), LDC )
         ELSE
*
*           w(1:lastc,1) := C(firstv:lastv-1,1:lastc)**T * v(firstv:lastv-1,1)
*
            CALL ZGEMV( 'Conjugate transpose', LASTV - FIRSTV, LASTC,
     $                  ONE, C( FIRSTV, 1 ), LDC, V( I ), INCV, ZERO,
     $                  WORK, 1 )
*
*           w(1:lastc,1) += C(lastv,1:lastc)**H * v(lastv,1)
*
            DO J = 1, LASTC
               WORK( J ) = WORK( J ) + CONJG( C( LASTV, J ) )
            END DO
*
*           C(lastv,1:lastc) += - tau * v(lastv,1) * w(1:lastc,1)**H
*
            DO J = 1, LASTC
               C( LASTV, J ) = C( LASTV, J )
     $                         - TAU * CONJG( WORK( J ) )
            END DO
*
*           C(firstv:lastv-1,1:lastc) += - tau * v(firstv:lastv-1,1) * w(1:lastc,1)**H
*
            CALL ZGERC( LASTV - FIRSTV, LASTC, -TAU, V( I ), INCV,
     $                  WORK, 1, C( FIRSTV, 1 ), LDC)
         END IF
      ELSE
*
*        Form  C * H
*
         IF( LASTV.EQ.FIRSTV ) THEN
*
*           C(1:lastc,lastv) := ( 1 - tau ) * C(1:lastc,lastv)
*
            CALL ZSCAL( LASTC, ONE - TAU, C( 1, LASTV ), 1 )
         ELSE
*
*           w(1:lastc,1) := C(1:lastc,firstv:lastv-1) * v(firstv:lastv-1,1)
*
            CALL ZGEMV( 'No transpose', LASTC, LASTV - FIRSTV, ONE,
     $                  C( 1, FIRSTV ), LDC, V( I ), INCV, ZERO,
     $                  WORK, 1 )
*
*           w(1:lastc,1) += C(1:lastc,lastv) * v(lastv,1)
*
            CALL ZAXPY( LASTC, ONE, C( 1, LASTV ), 1, WORK, 1 )
*
*           C(1:lastc,lastv) += - tau * v(lastv,1) * w(1:lastc,1)
*
            CALL ZAXPY( LASTC, -TAU, WORK, 1, C( 1, LASTV ), 1 )
*
*           C(1:lastc,firstv:lastv-1) += - tau * w(1:lastc,1) * v(firstv:lastv-1)**H
*
            CALL ZGERC( LASTC, LASTV - FIRSTV, -TAU, WORK, 1, V( I ),
     $                  INCV, C( 1, FIRSTV ), LDC )
         END IF
      END IF
      RETURN
*
*     End of ZLARF1L
*
      END
