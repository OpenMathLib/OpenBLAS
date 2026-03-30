*> \brief \b ZLARF1F applies an elementary reflector to a general rectangular
*              matrix assuming v(1) = 1.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download ZLARF1F + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarf1f.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarf1f.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarf1f.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLARF1F( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
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
*> ZLARF1F applies a complex elementary reflector H to a real m by n matrix
*> C, from either the left or the right. H is represented in the form
*>
*>       H = I - tau * v * v**H
*>
*> where tau is a complex scalar and v is a complex vector.
*>
*> If tau = 0, then H is taken to be the unit matrix.
*>
*> To apply H**H, supply conjg(tau) instead
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
*>          TAU = 0. V(1) is not referenced or modified.
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
*  To take advantage of the fact that v(1) = 1, we do the following
*     v = [ 1 v_2 ]**T
*     If SIDE='L'
*           |-----|
*           | C_1 |
*        C =| C_2 |
*           |-----|
*        C_1\in\mathbb{C}^{1\times n}, C_2\in\mathbb{C}^{m-1\times n}
*        So we compute:
*        C = HC   = (I - \tau vv**T)C
*                 = C - \tau vv**T C
*        w = C**T v  = [ C_1**T C_2**T ] [ 1 v_2 ]**T
*                    = C_1**T + C_2**T v ( ZGEMM then ZAXPYC-like )
*        C  = C - \tau vv**T C
*           = C - \tau vw**T
*        Giving us   C_1 = C_1 - \tau w**T ( ZAXPYC-like )
*                 and
*                    C_2 = C_2 - \tau v_2w**T ( ZGERC )
*     If SIDE='R'
*
*        C = [ C_1 C_2 ]
*        C_1\in\mathbb{C}^{m\times 1}, C_2\in\mathbb{C}^{m\times n-1}
*        So we compute: 
*        C = CH   = C(I - \tau vv**T)
*                 = C - \tau Cvv**T
*
*        w = Cv   = [ C_1 C_2 ] [ 1 v_2 ]**T
*                 = C_1 + C_2v_2 ( ZGEMM then ZAXPYC-like )
*        C  = C - \tau Cvv**T
*           = C - \tau wv**T
*        Giving us   C_1 = C_1 - \tau w ( ZAXPYC-like )
*                 and
*                    C_2 = C_2 - \tau wv_2**T ( ZGERC )
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
      SUBROUTINE ZLARF1F( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
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
      INTEGER            I, LASTV, LASTC, J
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZAXPY, ZGEMV, ZGERC, ZSCAL
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
      LASTV = 1
      LASTC = 0
      IF( TAU.NE.ZERO ) THEN
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V.
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
!        Since we are assuming that V(1) = 1, and it is not stored, so we
!        shouldn't access it.
         DO WHILE( LASTV.GT.1 .AND. V( I ).EQ.ZERO )
            LASTV = LASTV - 1
            I = I - INCV
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
            ! Check if m = 1. This means v = 1, So we just need to compute
            ! C := HC = (1-\tau)C.
            IF( LASTV.EQ.1 ) THEN
               CALL ZSCAL(LASTC, ONE - TAU, C, LDC)
            ELSE
*
*              w(1:lastc,1) := C(1:lastv,1:lastc)**H * v(1:lastv,1)
*
               ! (I - tvv**H)C = C - tvv**H C
               ! First compute w**H = v**H c -> w = C**H v
               ! C = [ C_1 C_2 ]**T, v = [1 v_2]**T
               ! w = C_1**H + C_2**Hv_2
               ! w = C_2**Hv_2
               CALL ZGEMV( 'Conjugate transpose', LASTV - 1,
     $               LASTC, ONE, C( 1+1, 1 ), LDC, V( 1 + INCV ),
     $               INCV, ZERO, WORK, 1 )
*
*              w(1:lastc,1) += v(1,1) * C(1,1:lastc)**H
*
               DO I = 1, LASTC
                  WORK( I ) = WORK( I ) + DCONJG( C( 1, I ) )
               END DO
*
*           C(1:lastv,1:lastc) := C(...) - tau * v(1:lastv,1) * w(1:lastc,1)**H
*
            ! C(1, 1:lastc)   := C(...) - tau * v(1,1) * w(1:lastc,1)**H
            !                  = C(...) - tau * Conj(w(1:lastc,1))
            ! This is essentially a zaxpyc
               DO I = 1, LASTC
                  C( 1, I ) = C( 1, I ) - TAU * DCONJG( WORK( I ) )
               END DO
*
*        C(2:lastv,1:lastc) += - tau * v(2:lastv,1) * w(1:lastc,1)**H
*
               CALL ZGERC( LASTV - 1, LASTC, -TAU, V( 1 + INCV ),
     $               INCV, WORK, 1, C( 1+1, 1 ), LDC )
            END IF
      ELSE
*
*        Form  C * H
*
            ! Check if n = 1. This means v = 1, so we just need to compute
            ! C := CH = C(1-\tau).
            IF( LASTV.EQ.1 ) THEN
               CALL ZSCAL(LASTC, ONE - TAU, C, 1)
            ELSE
*
*              w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
*
               ! w(1:lastc,1) := C(1:lastc,2:lastv) * v(2:lastv,1)
               CALL ZGEMV( 'No transpose', LASTC, LASTV-1, ONE, 
     $            C(1,1+1), LDC, V(1+INCV), INCV, ZERO, WORK, 1 )
               ! w(1:lastc,1) += C(1:lastc,1) v(1,1) = C(1:lastc,1)
               CALL ZAXPY(LASTC, ONE, C, 1, WORK, 1)
*
*              C(1:lastc,1:lastv) := C(...) - tau * w(1:lastc,1) * v(1:lastv,1)**T
*
               ! C(1:lastc,1)     := C(...) - tau * w(1:lastc,1) * v(1,1)**T
               !                   = C(...) - tau * w(1:lastc,1)
               CALL ZAXPY(LASTC, -TAU, WORK, 1, C, 1)
               ! C(1:lastc,2:lastv) := C(...) - tau * w(1:lastc,1) * v(2:lastv)**T
               CALL ZGERC( LASTC, LASTV-1, -TAU, WORK, 1, V(1+INCV),
     $                     INCV, C(1,1+1), LDC )
            END IF
      END IF
      RETURN
*
*     End of ZLARF1F
*
      END
