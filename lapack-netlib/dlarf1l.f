*> \brief \b DLARF1L applies an elementary reflector to a general rectangular
*              matrix assuming v(lastv) = 1 where lastv is the last non-zero
*              element
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DLARF1L + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarf1l.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarf1l.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarf1l.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLARF1L( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*       .. Scalar Arguments ..
*       CHARACTER          SIDE
*       INTEGER            INCV, LDC, M, N
*       DOUBLE PRECISION   TAU
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLARF1L applies a real elementary reflector H to a real m by n matrix
*> C, from either the left or the right. H is represented in the form
*>
*>       H = I - tau * v * v**T
*>
*> where tau is a real scalar and v is a real vector.
*>
*> If tau = 0, then H is taken to be the unit matrix.
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
*>          V is DOUBLE PRECISION array, dimension
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
*>          TAU is DOUBLE PRECISION
*>          The value tau in the representation of H.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array, dimension (LDC,N)
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
*>          WORK is DOUBLE PRECISION array, dimension
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
*> \ingroup larf
*
*  =====================================================================
      SUBROUTINE DLARF1L( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      DOUBLE PRECISION   TAU
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, FIRSTV, LASTV, LASTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DGEMV, DGER, DSCAL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILADLR, ILADLC
      EXTERNAL           LSAME, ILADLR, ILADLC
*     ..
*     .. Executable Statements ..
*
      APPLYLEFT = LSAME( SIDE, 'L' )
      FIRSTV = 1
      LASTC = 0
      IF( TAU.NE.ZERO ) THEN
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V.
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
            LASTC = ILADLC(LASTV, N, C, LDC)
         ELSE
!     Scan for the last non-zero row in C(:,1:lastv).
            LASTC = ILADLR(M, LASTV, C, LDC)
         END IF
      END IF
      IF( LASTC.EQ.0 ) THEN
         RETURN
      END IF
      IF( APPLYLEFT ) THEN
*
*        Form  H * C
*
         IF( LASTV.GT.0 ) THEN
            ! Check if m = 1. This means v = 1, So we just need to compute
            ! C := HC = (1-\tau)C.
            IF( LASTV.EQ.FIRSTV ) THEN
               CALL DSCAL(LASTC, ONE - TAU, C( FIRSTV, 1), LDC)
            ELSE
*
*              w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
*
               ! w(1:lastc,1) := C(1:lastv-1,1:lastc)**T * v(1:lastv-1,1)
               CALL DGEMV( 'Transpose', LASTV-FIRSTV, LASTC, ONE,
     $                     C(FIRSTV,1), LDC, V(I), INCV, ZERO, 
     $                     WORK, 1)
               ! w(1:lastc,1) += C(lastv,1:lastc)**T * v(lastv,1) = C(lastv,1:lastc)**T
               CALL DAXPY(LASTC, ONE, C(LASTV,1), LDC, WORK, 1)
*
*              C(1:lastv,1:lastc) := C(...) - tau * v(1:lastv,1) * w(1:lastc,1)**T
*
               ! C(lastv, 1:lastc)   := C(...) - tau * v(lastv,1) * w(1:lastc,1)**T
               !                      = C(...) - tau * w(1:lastc,1)**T
               CALL DAXPY(LASTC, -TAU, WORK, 1, C(LASTV,1), LDC)
               ! C(1:lastv-1,1:lastc) := C(...) - tau * v(1:lastv-1,1)*w(1:lastc,1)**T
               CALL DGER(LASTV-FIRSTV, LASTC, -TAU, V(I), INCV,
     $                   WORK, 1, C(FIRSTV,1), LDC)
            END IF
         END IF
      ELSE
*
*        Form  C * H
*
         IF( LASTV.GT.0 ) THEN
            ! Check if n = 1. This means v = 1, so we just need to compute
            ! C := CH = C(1-\tau).
            IF( LASTV.EQ.FIRSTV ) THEN
               CALL DSCAL(LASTC, ONE - TAU, C, 1)
            ELSE
*
*              w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
*
               ! w(1:lastc,1) := C(1:lastc,1:lastv-1) * v(1:lastv-1,1)
               CALL DGEMV( 'No transpose', LASTC, LASTV-FIRSTV,
     $            ONE, C(1,FIRSTV), LDC, V(I), INCV, ZERO, WORK, 1 )
               ! w(1:lastc,1) += C(1:lastc,lastv) * v(lastv,1) = C(1:lastc,lastv)
               CALL DAXPY(LASTC, ONE, C(1,LASTV), 1, WORK, 1)
*
*              C(1:lastc,1:lastv) := C(...) - tau * w(1:lastc,1) * v(1:lastv,1)**T
*
               ! C(1:lastc,lastv)     := C(...) - tau * w(1:lastc,1) * v(lastv,1)**T
               !                       = C(...) - tau * w(1:lastc,1)
               CALL DAXPY(LASTC, -TAU, WORK, 1, C(1,LASTV), 1)
               ! C(1:lastc,1:lastv-1) := C(...) - tau * w(1:lastc,1) * v(1:lastv-1)**T
               CALL DGER( LASTC, LASTV-FIRSTV, -TAU, WORK, 1, V(I), 
     $                     INCV, C(1,FIRSTV), LDC )
            END IF
         END IF
      END IF
      RETURN
*
*     End of DLARF1L
*
      END
