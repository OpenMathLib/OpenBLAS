*> \brief \b SLAR1V_TGK computes one eigenvector of a shifted tridiagonal from a twisted factorization, keeping the full support (Golub-Kahan variant of SLAR1V); also SLAR1V2_TGK, the two-eigenvalue variant.
*>
*> \ingroup bdsvdmr3
*>
*> \par Purpose:
*> =============
*>
*> \verbatim
*> computes one eigenvector of a shifted tridiagonal from a twisted factorization, keeping the full support (Golub-Kahan variant of SLAR1V); also SLAR1V2_TGK, the two-eigenvalue variant.
*> All routines in this file are auxiliary to SBDSVDMR3 and derive
*> from the LAPACK 3.0 (1999) MRRR kernels SLARRV, SLAR1V, SLARRB and
*> SLARRF rather than from the current SSTEMR kernels.
*> \endverbatim
*
      SUBROUTINE SLAR1V_TGK( REP, N, B1, BN, SIGMA, D, L, LD, LLD, EVAL,
     $                   GERSCH, Z, ZTZ, MINGMA, R, ISUPPZ, WORK )
*
*  -- New auxiliary routine for bidiagonal SVD via TGK-rooted MR^3 --
*     Based on LAPACK SLAR1V (Dhillon & Marques, Nov 11 2003).
*     Derived from the LAPACK 3.0 SLAR1V; this is a separate copy that adds a
*     representation switch so the matrix M may be supplied either as an
*     L D L^T factorization (REP='L', identical to the original) or as a
*     symmetric tridiagonal -- in particular the Tridiagonal Golub-Kahan
*     (TGK) matrix -- with REP='T'.
*
*     .. Scalar Arguments ..
      CHARACTER          REP
      INTEGER            B1, BN, N, R
      REAL   EVAL, MINGMA, SIGMA, ZTZ
*     ..
*     .. Array Arguments ..
      INTEGER            ISUPPZ( * )
      REAL   D( * ), GERSCH( * ), L( * ), LD( * ), LLD( * ),
     $                   WORK( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  SLAR1V_TGK computes the (scaled) r-th column of the inverse of the
*  submatrix in rows B1 through BN of  M - sigma I,  where M is either
*  an L D L^T representation (REP='L') or a symmetric tridiagonal matrix
*  such as the TGK matrix (REP='T').  When sigma is close to an
*  eigenvalue, the computed vector is an accurate eigenvector.  The steps
*  are those of SLAR1V:
*  (a) stationary qd transform   M - sigma I = L(+) D(+) L(+)^T,
*  (b) progressive qd transform  M - sigma I = U(-) D(-) U(-)^T,
*  (c) the twist index r = argmin_k |gamma_k|, where gamma_k is the k-th
*      diagonal element of the twisted factorization and equals the
*      reciprocal of (M - sigma I)^{-1}_{kk},
*  (d) the scaled r-th column of the inverse from the twisted factor.
*
*  For REP='T' the qd transforms run on the tridiagonal entries directly
*  (diagonal in D, off-diagonals beta_i in L), so LD and LLD are not
*  referenced.  Using the identity  gamma_k = D(+)(k) + D(-)(k) - J_kk
*  with J_kk = D(k) - sigma (Dhillon, PhD thesis, UC Berkeley 1997, eq. 3.1.12), the WORK slots are filled
*  so that the twist selection and back-substitution code below are
*  shared, unchanged, with the REP='L' path:
*       WORK(INDS + k-1) = D(+)(k)
*       WORK(INDP + k-1) = D(-)(k) - ( D(k) - sigma )
*       WORK(k)          = L(+)(k) = beta_k / D(+)(k)
*       WORK(INDUMN + k) = U(-)(k) = beta_k / D(-)(k+1)
*  so WORK(INDS+k-1) + WORK(INDP+k-1) = gamma_k.  For the TGK matrix
*  (D == 0) this is gamma_k = D(+)(k) + D(-)(k) + sigma.
*
*  Zero pivots are expected for the indefinite zero-diagonal TGK root.
*  In the direct recurrence a zero pivot self-heals through IEEE
*  arithmetic (D(+) -> -Inf, then L(+) -> 0), so no stationary/progressive
*  restart is needed; only the back-substitution uses a slow path (to
*  avoid Inf*0 = NaN), with beta-ratios per Algorithm 3.3.1 of Dhillon's thesis (UC Berkeley, 1997).
*
*  Arguments
*  =========
*
*  REP   (input) CHARACTER*1   'L' = L D L^T input; 'T' = tridiagonal/TGK.
*  N,B1,BN,SIGMA,EVAL,GERSCH,Z,ZTZ,MINGMA,R,ISUPPZ,WORK : as in SLAR1V.
*  D     (input) REAL (N)   REP='L': diag(D).  REP='T': diag(M).
*  L     (input) REAL (N-1) REP='L': subdiag(L). REP='T': beta_i.
*  LD,LLD(input) REAL (N-1) REP='L': L*D, L*L*D. REP='T': unused.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLKSIZ
      PARAMETER          ( BLKSIZ = 32 )
      REAL   ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            SAWNAN, TGK
      INTEGER            FROM, I, INDP, INDS, INDUMN, J, R1, R2, TO
      REAL   DMINUS, DPLUS, EPS, S, TMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL   SLAMCH
      EXTERNAL           LSAME, SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      TGK = LSAME( REP, 'T' )
      EPS = SLAMCH( 'Precision' )
      IF( R.EQ.0 ) THEN
*
*        Eliminate the top and bottom indices from among the possible
*        values of R, using the Gerschgorin intervals.
*
         R1 = B1
         R2 = BN
         DO 10 I = B1, BN
            IF( EVAL.GE.GERSCH( 2*I-1 ) .AND. EVAL.LE.GERSCH( 2*I ) )
     $           THEN
               R1 = I
               GO TO 20
            END IF
   10    CONTINUE
         GO TO 40
   20    CONTINUE
         DO 30 I = BN, B1, -1
            IF( EVAL.GE.GERSCH( 2*I-1 ) .AND. EVAL.LE.GERSCH( 2*I ) )
     $           THEN
               R2 = I
               GO TO 40
            END IF
   30    CONTINUE
      ELSE
         R1 = R
         R2 = R
      END IF
*
   40 CONTINUE
      INDUMN = N
      INDS = 2*N + 1
      INDP = 3*N + 1
      SAWNAN = .FALSE.
*
*     ================================================================
*     Stationary transform (differential form) up to the index R2.
*     ================================================================
*
      IF( TGK ) THEN
*
*        Tridiagonal/TGK parent.  Direct recurrence; store D(+)(k).
*
         DPLUS = D( B1 ) - SIGMA
         WORK( INDS+B1-1 ) = DPLUS
         IF( DPLUS.EQ.ZERO ) SAWNAN = .TRUE.
         DO 45 I = B1, R2 - 1
            WORK( I ) = L( I ) / DPLUS
            DPLUS = D( I+1 ) - SIGMA - WORK( I )*L( I )
            WORK( INDS+I ) = DPLUS
            IF( DPLUS.EQ.ZERO ) SAWNAN = .TRUE.
   45    CONTINUE
      ELSE
*
*        L D L^T parent (original SLAR1V stationary loop).
*
         IF( B1.EQ.1 ) THEN
            WORK( INDS ) = ZERO
         ELSE
            WORK( INDS ) = LLD( B1-1 )
         END IF
         S = WORK( INDS ) - SIGMA
         DO 50 I = B1, R2 - 1
            DPLUS = D( I ) + S
            WORK( I ) = LD( I ) / DPLUS
            WORK( INDS+I ) = S*WORK( I )*L( I )
            S = WORK( INDS+I ) - SIGMA
   50    CONTINUE
*
         IF( .NOT.( S.GT.ZERO .OR. S.LT.ONE ) ) THEN
*           Slower version of the above loop if a NaN is detected.
            SAWNAN = .TRUE.
            J = B1 + 1
   60       CONTINUE
            IF( WORK( INDS+J ).GT.ZERO .OR. WORK( INDS+J ).LT.ONE ) THEN
               J = J + 1
               GO TO 60
            END IF
            WORK( INDS+J ) = LLD( J )
            S = WORK( INDS+J ) - SIGMA
            DO 70 I = J + 1, R2 - 1
               DPLUS = D( I ) + S
               WORK( I ) = LD( I ) / DPLUS
               IF( WORK( I ).EQ.ZERO ) THEN
                  WORK( INDS+I ) = LLD( I )
               ELSE
                  WORK( INDS+I ) = S*WORK( I )*L( I )
               END IF
               S = WORK( INDS+I ) - SIGMA
   70       CONTINUE
         END IF
      END IF
*
*     ================================================================
*     Progressive transform (differential form) down to the index R1.
*     ================================================================
*
      IF( TGK ) THEN
*
*        Tridiagonal/TGK parent.  Direct recurrence.  Store
*        WORK(INDP+i-1) = D(-)(i) - (D(i)-sigma)  and  U(-)(i).
*
         DMINUS = D( BN ) - SIGMA
         WORK( INDP+BN-1 ) = ZERO
         IF( DMINUS.EQ.ZERO ) SAWNAN = .TRUE.
         DO 75 I = BN - 1, R1, -1
            WORK( INDUMN+I ) = L( I ) / DMINUS
            DMINUS = D( I ) - SIGMA - WORK( INDUMN+I )*L( I )
            WORK( INDP+I-1 ) = DMINUS - ( D( I ) - SIGMA )
            IF( DMINUS.EQ.ZERO ) SAWNAN = .TRUE.
   75    CONTINUE
      ELSE
*
*        L D L^T parent (original SLAR1V progressive loop).
*
         WORK( INDP+BN-1 ) = D( BN ) - SIGMA
         DO 80 I = BN - 1, R1, -1
            DMINUS = LLD( I ) + WORK( INDP+I )
            TMP = D( I ) / DMINUS
            WORK( INDUMN+I ) = L( I )*TMP
            WORK( INDP+I-1 ) = WORK( INDP+I )*TMP - SIGMA
   80    CONTINUE
         TMP = WORK( INDP+R1-1 )
         IF( .NOT.( TMP.GT.ZERO .OR. TMP.LT.ONE ) ) THEN
*           Slower version of the above loop if a NaN is detected.
            SAWNAN = .TRUE.
            J = BN - 3
   90       CONTINUE
            IF( WORK( INDP+J ).GT.ZERO .OR. WORK( INDP+J ).LT.ONE ) THEN
               J = J - 1
               GO TO 90
            END IF
            WORK( INDP+J ) = D( J+1 ) - SIGMA
            DO 100 I = J, R1, -1
               DMINUS = LLD( I ) + WORK( INDP+I )
               TMP = D( I ) / DMINUS
               WORK( INDUMN+I ) = L( I )*TMP
               IF( TMP.EQ.ZERO ) THEN
                  WORK( INDP+I-1 ) = D( I ) - SIGMA
               ELSE
                  WORK( INDP+I-1 ) = WORK( INDP+I )*TMP - SIGMA
               END IF
  100       CONTINUE
         END IF
      END IF
*
*     ================================================================
*     Twist index: r = argmin_k |gamma_k|, gamma_k stored as
*     WORK(INDS+k-1) + WORK(INDP+k-1) for both representations.
*     ================================================================
*
      MINGMA = WORK( INDS+R1-1 ) + WORK( INDP+R1-1 )
      IF( MINGMA.EQ.ZERO )
     $   MINGMA = EPS*WORK( INDS+R1-1 )
      R = R1
      DO 110 I = R1, R2 - 1
         TMP = WORK( INDS+I ) + WORK( INDP+I )
         IF( TMP.EQ.ZERO )
     $      TMP = EPS*WORK( INDS+I )
         IF( ABS( TMP ).LT.ABS( MINGMA ) ) THEN
            MINGMA = TMP
            R = I + 1
         END IF
  110 CONTINUE
*
*     ================================================================
*     Compute the (scaled) r-th column of the inverse.  WORK(k)=L(+)(k)
*     and WORK(INDUMN+k)=U(-)(k) for both representations, so the fast
*     path is shared.  The slow (SAWNAN) path uses LD-ratios for REP='L'
*     and beta-ratios (L) for REP='T'.
*     ================================================================
*
      ISUPPZ( 1 ) = B1
      ISUPPZ( 2 ) = BN
      Z( R ) = ONE
      ZTZ = ONE
      IF( .NOT.SAWNAN ) THEN
*        No support truncation.  A TGK eigenvector can be BIMODAL -- the
*        left singular vector u occupies one end of z and the right vector
*        v the other, with a near-zero trough between -- so the usual
*        "two consecutive tiny entries => the eigenvector has died out,
*        stop" early-out would wrongly amputate the far hump.  This is
*        true for EVERY vector this routine computes, whether the
*        representation is the root TGK (REP='T') or a deeper L D L^T
*        child of it (REP='L'), since the eigenvectors are TGK
*        eigenvectors in both cases.  Hence we always form the full
*        vector and leave ISUPPZ = [B1,BN].
         DO 130 I = R - 1, B1, -1
            Z( I ) = -( WORK( I )*Z( I+1 ) )
            ZTZ = ZTZ + Z( I )*Z( I )
  130    CONTINUE
         DO 150 I = R + 1, BN
            Z( I ) = -( WORK( INDUMN+I-1 )*Z( I-1 ) )
            ZTZ = ZTZ + Z( I )*Z( I )
  150    CONTINUE
      ELSE
*
*        Slow path: some pivot vanished, so L(+)/U(-) may be +-Inf.
*        Use the recurrence of Algorithm 3.3.1 in Dhillon's thesis, switching the
*        coupling ratio between LD (REP='L') and beta=L (REP='T').
*
         DO 160 I = R - 1, B1, -1
            IF( Z( I+1 ).EQ.ZERO ) THEN
               IF( TGK ) THEN
                  Z( I ) = -( L( I+1 ) / L( I ) )*Z( I+2 )
               ELSE
                  Z( I ) = -( LD( I+1 ) / LD( I ) )*Z( I+2 )
               END IF
            ELSE
               Z( I ) = -( WORK( I )*Z( I+1 ) )
            END IF
            ZTZ = ZTZ + Z( I )*Z( I )
  160    CONTINUE
  170    CONTINUE
         DO 180 I = R, BN - 1
            IF( Z( I ).EQ.ZERO ) THEN
               IF( TGK ) THEN
                  Z( I+1 ) = -( L( I-1 ) / L( I ) )*Z( I-1 )
               ELSE
                  Z( I+1 ) = -( LD( I-1 ) / LD( I ) )*Z( I-1 )
               END IF
            ELSE
               Z( I+1 ) = -( WORK( INDUMN+I )*Z( I ) )
            END IF
            ZTZ = ZTZ + Z( I+1 )*Z( I+1 )
  180    CONTINUE
  190    CONTINUE
      END IF
*
      RETURN
*
*     End of SLAR1V_TGK
*
      END

      SUBROUTINE SLAR1V2_TGK( REP, N, B1, BN, SIGMA1, SIGMA2, D, L,
     $                   LD, LLD, EVAL1, EVAL2, GERSCH, Z1, Z2,
     $                   ZTZ1, ZTZ2, MINGMA1, MINGMA2, RONE, RTWO,
     $                   ISUP1, ISUP2, WORK1, WORK2 )
*
*     Two independent SLAR1V_TGK fast paths interleaved to expose the
*     serial qd divisions to the CPU.  Every operation within either lane
*     has the same order as the scalar routine.  Exceptional zero/NaN
*     pivots fall back to two scalar calls, preserving their exact path.
*
      CHARACTER          REP
      INTEGER            B1, BN, N, RONE, RTWO
      INTEGER            ISUP1( * ), ISUP2( * )
      REAL   D( * ), EVAL1, EVAL2, GERSCH( * ), L( * ),
     $                   LD( * ), LLD( * ), MINGMA1, MINGMA2,
     $                   SIGMA1, SIGMA2, WORK1( * ), WORK2( * ),
     $                   Z1( * ), Z2( * ), ZTZ1, ZTZ2
      REAL   ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
      LOGICAL            BAD1, BAD2, TGK
      INTEGER            I, INDP, INDS, INDUMN, R1A, R1B, R2A, R2B
      REAL   DMA, DMB, DPA, DPB, EPS, SA, SB, TMPA, TMPB
      LOGICAL            LSAME
      REAL   SLAMCH
      EXTERNAL           SLAMCH, LSAME
      INTRINSIC          ABS, MAX, MIN
*
      TGK = LSAME( REP, 'T' )
      EPS = SLAMCH( 'Precision' )
      IF( RONE.EQ.0 ) THEN
         R1A = B1
         R2A = BN
         DO 10 I = B1, BN
            IF( EVAL1.GE.GERSCH( 2*I-1 ) .AND.
     $          EVAL1.LE.GERSCH( 2*I ) ) THEN
               R1A = I
               GO TO 20
            END IF
   10    CONTINUE
         GO TO 40
   20    CONTINUE
         DO 30 I = BN, B1, -1
            IF( EVAL1.GE.GERSCH( 2*I-1 ) .AND.
     $          EVAL1.LE.GERSCH( 2*I ) ) THEN
               R2A = I
               GO TO 40
            END IF
   30    CONTINUE
   40    CONTINUE
      ELSE
         R1A = RONE
         R2A = RONE
      END IF
      IF( RTWO.EQ.0 ) THEN
         R1B = B1
         R2B = BN
         DO 50 I = B1, BN
            IF( EVAL2.GE.GERSCH( 2*I-1 ) .AND.
     $          EVAL2.LE.GERSCH( 2*I ) ) THEN
               R1B = I
               GO TO 60
            END IF
   50    CONTINUE
         GO TO 80
   60    CONTINUE
         DO 70 I = BN, B1, -1
            IF( EVAL2.GE.GERSCH( 2*I-1 ) .AND.
     $          EVAL2.LE.GERSCH( 2*I ) ) THEN
               R2B = I
               GO TO 80
            END IF
   70    CONTINUE
   80    CONTINUE
      ELSE
         R1B = RTWO
         R2B = RTWO
      END IF
*
      INDUMN = N
      INDS = 2*N + 1
      INDP = 3*N + 1
      BAD1 = .FALSE.
      BAD2 = .FALSE.
      IF( TGK ) THEN
         DPA = D( B1 ) - SIGMA1
         WORK1( INDS+B1-1 ) = DPA
         IF( DPA.EQ.ZERO ) BAD1 = .TRUE.
         DPB = D( B1 ) - SIGMA2
         WORK2( INDS+B1-1 ) = DPB
         IF( DPB.EQ.ZERO ) BAD2 = .TRUE.
         DO 90 I = B1, MAX( R2A, R2B ) - 1
            IF( I.LT.R2A ) THEN
               WORK1( I ) = L( I ) / DPA
               DPA = D( I+1 ) - SIGMA1 - WORK1( I )*L( I )
               WORK1( INDS+I ) = DPA
               IF( DPA.EQ.ZERO ) BAD1 = .TRUE.
            END IF
            IF( I.LT.R2B ) THEN
               WORK2( I ) = L( I ) / DPB
               DPB = D( I+1 ) - SIGMA2 - WORK2( I )*L( I )
               WORK2( INDS+I ) = DPB
               IF( DPB.EQ.ZERO ) BAD2 = .TRUE.
            END IF
   90    CONTINUE
      ELSE
         IF( B1.EQ.1 ) THEN
            WORK1( INDS ) = ZERO
            WORK2( INDS ) = ZERO
         ELSE
            WORK1( INDS ) = LLD( B1-1 )
            WORK2( INDS ) = LLD( B1-1 )
         END IF
         SA = WORK1( INDS ) - SIGMA1
         SB = WORK2( INDS ) - SIGMA2
         DO 100 I = B1, MAX( R2A, R2B ) - 1
            IF( I.LT.R2A ) THEN
               DPA = D( I ) + SA
               WORK1( I ) = LD( I ) / DPA
               WORK1( INDS+I ) = SA*WORK1( I )*L( I )
               SA = WORK1( INDS+I ) - SIGMA1
            END IF
            IF( I.LT.R2B ) THEN
               DPB = D( I ) + SB
               WORK2( I ) = LD( I ) / DPB
               WORK2( INDS+I ) = SB*WORK2( I )*L( I )
               SB = WORK2( INDS+I ) - SIGMA2
            END IF
  100    CONTINUE
         IF( .NOT.( SA.GT.ZERO .OR. SA.LT.ONE ) ) BAD1 = .TRUE.
         IF( .NOT.( SB.GT.ZERO .OR. SB.LT.ONE ) ) BAD2 = .TRUE.
      END IF
*
      IF( TGK ) THEN
         DMA = D( BN ) - SIGMA1
         WORK1( INDP+BN-1 ) = ZERO
         IF( DMA.EQ.ZERO ) BAD1 = .TRUE.
         DMB = D( BN ) - SIGMA2
         WORK2( INDP+BN-1 ) = ZERO
         IF( DMB.EQ.ZERO ) BAD2 = .TRUE.
         DO 110 I = BN - 1, MIN( R1A, R1B ), -1
            IF( I.GE.R1A ) THEN
               WORK1( INDUMN+I ) = L( I ) / DMA
               DMA = D( I ) - SIGMA1 - WORK1( INDUMN+I )*L( I )
               WORK1( INDP+I-1 ) = DMA - ( D( I ) - SIGMA1 )
               IF( DMA.EQ.ZERO ) BAD1 = .TRUE.
            END IF
            IF( I.GE.R1B ) THEN
               WORK2( INDUMN+I ) = L( I ) / DMB
               DMB = D( I ) - SIGMA2 - WORK2( INDUMN+I )*L( I )
               WORK2( INDP+I-1 ) = DMB - ( D( I ) - SIGMA2 )
               IF( DMB.EQ.ZERO ) BAD2 = .TRUE.
            END IF
  110    CONTINUE
      ELSE
         WORK1( INDP+BN-1 ) = D( BN ) - SIGMA1
         WORK2( INDP+BN-1 ) = D( BN ) - SIGMA2
         DO 120 I = BN - 1, MIN( R1A, R1B ), -1
            IF( I.GE.R1A ) THEN
               DMA = LLD( I ) + WORK1( INDP+I )
               TMPA = D( I ) / DMA
               WORK1( INDUMN+I ) = L( I )*TMPA
               WORK1( INDP+I-1 ) = WORK1( INDP+I )*TMPA - SIGMA1
            END IF
            IF( I.GE.R1B ) THEN
               DMB = LLD( I ) + WORK2( INDP+I )
               TMPB = D( I ) / DMB
               WORK2( INDUMN+I ) = L( I )*TMPB
               WORK2( INDP+I-1 ) = WORK2( INDP+I )*TMPB - SIGMA2
            END IF
  120    CONTINUE
         TMPA = WORK1( INDP+R1A-1 )
         TMPB = WORK2( INDP+R1B-1 )
         IF( .NOT.( TMPA.GT.ZERO .OR. TMPA.LT.ONE ) ) BAD1 = .TRUE.
         IF( .NOT.( TMPB.GT.ZERO .OR. TMPB.LT.ONE ) ) BAD2 = .TRUE.
      END IF
      IF( BAD1 .OR. BAD2 ) GO TO 300
*
      MINGMA1 = WORK1( INDS+R1A-1 ) + WORK1( INDP+R1A-1 )
      IF( MINGMA1.EQ.ZERO ) MINGMA1 = EPS*WORK1( INDS+R1A-1 )
      RONE = R1A
      DO 130 I = R1A, R2A - 1
         TMPA = WORK1( INDS+I ) + WORK1( INDP+I )
         IF( TMPA.EQ.ZERO ) TMPA = EPS*WORK1( INDS+I )
         IF( ABS( TMPA ).LT.ABS( MINGMA1 ) ) THEN
            MINGMA1 = TMPA
            RONE = I + 1
         END IF
  130 CONTINUE
      MINGMA2 = WORK2( INDS+R1B-1 ) + WORK2( INDP+R1B-1 )
      IF( MINGMA2.EQ.ZERO ) MINGMA2 = EPS*WORK2( INDS+R1B-1 )
      RTWO = R1B
      DO 140 I = R1B, R2B - 1
         TMPB = WORK2( INDS+I ) + WORK2( INDP+I )
         IF( TMPB.EQ.ZERO ) TMPB = EPS*WORK2( INDS+I )
         IF( ABS( TMPB ).LT.ABS( MINGMA2 ) ) THEN
            MINGMA2 = TMPB
            RTWO = I + 1
         END IF
  140 CONTINUE
*
      ISUP1( 1 ) = B1
      ISUP1( 2 ) = BN
      ISUP2( 1 ) = B1
      ISUP2( 2 ) = BN
      Z1( RONE ) = ONE
      Z2( RTWO ) = ONE
      ZTZ1 = ONE
      ZTZ2 = ONE
      DO 150 I = MAX( RONE, RTWO ) - 1, B1, -1
         IF( I.LT.RONE ) THEN
            Z1( I ) = -( WORK1( I )*Z1( I+1 ) )
            ZTZ1 = ZTZ1 + Z1( I )*Z1( I )
         END IF
         IF( I.LT.RTWO ) THEN
            Z2( I ) = -( WORK2( I )*Z2( I+1 ) )
            ZTZ2 = ZTZ2 + Z2( I )*Z2( I )
         END IF
  150 CONTINUE
      DO 160 I = MIN( RONE, RTWO ) + 1, BN
         IF( I.GT.RONE ) THEN
            Z1( I ) = -( WORK1( INDUMN+I-1 )*Z1( I-1 ) )
            ZTZ1 = ZTZ1 + Z1( I )*Z1( I )
         END IF
         IF( I.GT.RTWO ) THEN
            Z2( I ) = -( WORK2( INDUMN+I-1 )*Z2( I-1 ) )
            ZTZ2 = ZTZ2 + Z2( I )*Z2( I )
         END IF
  160 CONTINUE
      RETURN
*
  300 CONTINUE
      CALL SLAR1V_TGK( REP, N, B1, BN, SIGMA1, D, L, LD, LLD,
     $                 EVAL1, GERSCH, Z1, ZTZ1, MINGMA1, RONE,
     $                 ISUP1, WORK1 )
      CALL SLAR1V_TGK( REP, N, B1, BN, SIGMA2, D, L, LD, LLD,
     $                 EVAL2, GERSCH, Z2, ZTZ2, MINGMA2, RTWO,
     $                 ISUP2, WORK2 )
      RETURN
      END
