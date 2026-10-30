*> \brief \b SLARRV_TGK MRRR representation-tree traversal computing eigenvectors of the Golub-Kahan matrix (Golub-Kahan variant of the LAPACK 3.0 SLARRV).
*>
*> \ingroup bdsvdmr3
*>
*> \par Purpose:
*> =============
*>
*> \verbatim
*> MRRR representation-tree traversal computing eigenvectors of the Golub-Kahan matrix (Golub-Kahan variant of the LAPACK 3.0 SLARRV).
*> All routines in this file are auxiliary to SBDSVDMR3 and derive
*> from the LAPACK 3.0 (1999) MRRR kernels SLARRV, SLAR1V, SLARRB and
*> SLARRF rather than from the current SSTEMR kernels.
*> \endverbatim
*
      SUBROUTINE SLARRV_TGK( N, D, L, ISPLIT, M, W, IBLOCK, INDEXW,
     $                   GERSCH, TOL, Z, LDZ, ISUPPZ, WORK, IWORK,
     $                   INFO )
*
*  -- New auxiliary routine for bidiagonal SVD via TGK-rooted MR^3 --
*     Based on LAPACK SLARRV (Dhillon & Marques, Nov 11 2003).
*     Derived from the LAPACK 3.0 SLARRV.  This is a separate copy with one
*     structural change: the ROOT of the representation tree (tree depth
*     NDEPTH = 0) is a symmetric tridiagonal -- the Tridiagonal Golub-
*     Kahan (TGK) matrix -- rather than an L D L^T factorization.  At the
*     root the eigenvectors / child RRRs are obtained with the TGK-input
*     paths of SLAR1V_TGK / SLARRF_TGK (REP = 'T'); every node at depth
*     >= 1 is an ordinary L D L^T child and uses REP = 'L'.  SLARRB
*     retains that algorithm with an added bounded-progress guard;
*     its positive INFO is propagated to the driver for safe fallback.
*
*  Usage / calling convention (set up by the driver SBDSVDMR3)
*  ----------------------------------------------------------
*  On entry, for each split block the arrays D and L hold the TGK
*  representation of that block: D = 0 (the TGK diagonal) and L = the
*  TGK off-diagonals beta_i.  The block has TGK order IN = 2*nb where nb
*  is the bidiagonal block size; W holds the nb positive eigenvalues
*  (= singular values) of the block in ascending order, and
*      INDEXW(i) = nb + (position of W(i) among the block's positives)
*  because, after shifting by ~+sigma, the positive cluster occupies the
*  upper half of the block's 2*nb-point spectrum (all shifted negatives
*  lie below all shifted positives).  GERSCH holds the 2N TGK Gerschgorin
*  intervals from SBDTGK.  Otherwise the interface matches SLARRV.
*
*  =====================================================================
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDZ, M, N
      REAL   TOL
*     ..
*     .. Array Arguments ..
      INTEGER            IBLOCK( * ), INDEXW( * ), ISPLIT( * ),
     $                   ISUPPZ( * ), IWORK( * )
      REAL   D( * ), GERSCH( * ), L( * ), W( * ), WORK( * ),
     $                   Z( LDZ, * )
*     ..
*     .. Parameters ..
      INTEGER            MAXITR, MGSSIZ
      PARAMETER          ( MAXITR = 8, MGSSIZ = 1 )
      REAL   ZERO, ONE, FOUR, RELTHR
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, FOUR = 4.0E0,
     $                     RELTHR = 3.0E-3 )
*     ..
*     .. Local Scalars ..
      CHARACTER          REP
      LOGICAL            DONE1, DONE2, NOMGS, KTOT2OK
      INTEGER            I, IBEGIN, IEND, IINDC1, IINDC2, IINDR, IINDWK,
     $                   IINFO, IM, IN, INDERR, INDGAP, INDLD, INDLLD,
     $                   INDWRK, ITER, ITER2, ITMP1, ITMP2, J, JBLK, K,
     $                   K2, KTOT, KTOT2,
     $                   NCLUS, NDEPTH, NDONE, NEWCLS, NEWFRS, NEWFTT,
     $                   NEWLST, NEWSIZ, OLDCLS, OLDFST, OLDIEN, OLDLST,
     $                   OLDNCL, P, PARITY, Q, WBEGIN, WEND, ZFROM, ZTO
      REAL   EPS, GAP, GAP2, LAMBDA, LAMBDA2, MGSTOL,
     $                   MINGMA, MINGMA2, MINRGP, NRMINV, NRMINV2,
     $                   RELGAP, RELTOL, RESID, RESID2, RQCORR,
     $                   RQCORR2, SIGMA, TMP, TMP2, ZTZ, ZTZ2
*     ..
*     .. External Functions ..
      REAL   SDOT, SLAMCH, SNRM2
      EXTERNAL           SDOT, SLAMCH, SNRM2
*     ..
*     .. External Subroutines ..
      EXTERNAL           SAXPY, SCOPY, SLAR1V_TGK, SLAR1V2_TGK,
     $                   SLARRB_TGK, SLARRF_TGK, SLASET, SSCAL, SSTEIN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, REAL, MAX, MIN, SQRT
*     ..
*     .. Local Arrays ..
      INTEGER            TEMP( 1 ), ISEED( 4 )
*     ..
*     .. Executable Statements ..
*     ..
*     Fixed seed for the random relative perturbations that SLARRF_TGK
*     applies to each child RRR (Sec. 5.4).  Threaded (in/out) through all
*     SLARRF_TGK calls; fixed initial value makes the result reproducible.
      ISEED( 1 ) = 1
      ISEED( 2 ) = 13
      ISEED( 3 ) = 257
      ISEED( 4 ) = 2047
*
      INDERR = N
      INDLD = 2*N
      INDLLD = 3*N
      INDGAP = 4*N
      INDWRK = 5*N + 1
*
      IINDR = N
      IINDC1 = 2*N
      IINDC2 = 3*N
      IINDWK = 4*N + 1
*
      EPS = SLAMCH( 'Precision' )
*
      DO 10 I = 1, 2*N
         IWORK( I ) = 0
   10 CONTINUE
*     Zero the N-by-M eigenvector array (here N is the TGK order 2*nb and
*     M is the number of selected eigenvalues; Z has only M columns, unlike
*     the original SLARRV where the matrix order and column count coincide).
      CALL SLASET( 'Full', N, M, ZERO, ZERO, Z, LDZ )
      MGSTOL = 100.0E0*EPS
*
      IBEGIN = 1
      WBEGIN = 1
      DO 170 JBLK = 1, IBLOCK( M )
         IEND = ISPLIT( JBLK )
*
*        Find the eigenvectors of the submatrix indexed IBEGIN
*        through IEND.
*
         WEND = WBEGIN - 1
 171     CONTINUE
         IF( WEND.LT.M ) THEN
            IF( IBLOCK( WEND+1 ).EQ.JBLK ) THEN
               WEND = WEND + 1
               GO TO 171
            END IF
         END IF
         IF( WEND.LT.WBEGIN ) THEN
            IBEGIN = IEND + 1
            GO TO 170
         END IF
*
         IF( IBEGIN.EQ.IEND ) THEN
            Z( IBEGIN, WBEGIN ) = ONE
            ISUPPZ( 2*WBEGIN-1 ) = IBEGIN
            ISUPPZ( 2*WBEGIN ) = IBEGIN
            IBEGIN = IEND + 1
            WBEGIN = WEND + 1
            GO TO 170
         END IF
         OLDIEN = IBEGIN - 1
         IN = IEND - OLDIEN
         RELTOL = MIN( RELTHR, ONE / REAL( IN ) )
         IM = WEND - WBEGIN + 1
         CALL SCOPY( IM, W( WBEGIN ), 1, WORK, 1 )
         DO 30 I = 1, IM - 1
            WORK( INDERR+I ) = EPS*ABS( WORK( I ) )
            WORK( INDGAP+I ) = WORK( I+1 ) - WORK( I )
   30    CONTINUE
         WORK( INDERR+IM ) = EPS*ABS( WORK( IM ) )
         WORK( INDGAP+IM ) = MAX( ABS( WORK( IM ) ), EPS )
         NDONE = 0
*
         NDEPTH = 0
         PARITY = 1
         NCLUS = 1
         IWORK( IINDC1+1 ) = 1
         IWORK( IINDC1+2 ) = IM
*
*        While( NDONE.LT.IM ) do
*
   40    CONTINUE
         IF( NDONE.LT.IM .AND. NDEPTH.LE.4*N ) THEN
*
*           REP selects the representation of the current tree level:
*           'T' (TGK tridiagonal) at the root, 'L' (L D L^T) deeper down.
*
            IF( NDEPTH.EQ.0 ) THEN
               REP = 'T'
            ELSE
               REP = 'L'
            END IF
            OLDNCL = NCLUS
            NCLUS = 0
            PARITY = 1 - PARITY
            IF( PARITY.EQ.0 ) THEN
               OLDCLS = IINDC1
               NEWCLS = IINDC2
            ELSE
               OLDCLS = IINDC2
               NEWCLS = IINDC1
            END IF
            DO 150 I = 1, OLDNCL
*
*              If NDEPTH > 1, retrieve the relatively robust
*              representation (RRR) and perform limited bisection
*              (if necessary) to get approximate eigenvalues.
*
               J = OLDCLS + 2*I
               OLDFST = IWORK( J-1 )
               OLDLST = IWORK( J )
               IF( NDEPTH.GT.0 ) THEN
                  J = WBEGIN + OLDFST - 1
                  CALL SCOPY( IN, Z( IBEGIN, J ), 1, D( IBEGIN ), 1 )
                  CALL SCOPY( IN-1, Z( IBEGIN, J+1 ), 1, L( IBEGIN ),
     $               1 )
                  CALL SLASET( 'Full', IN, 2, ZERO, ZERO,
     $                               Z( IBEGIN, J ), LDZ )
               END IF
*
*              Compute LD and LLD of the current representation.  For the
*              TGK root (REP='T') the diagonal is zero so these are zero
*              and unused (SLARRB is skipped, and the TGK paths of
*              SLARRF_TGK/SLAR1V_TGK ignore LD,LLD); they are needed for
*              every L D L^T child (REP='L').
*
               K = IBEGIN
               DO 50 J = 1, IN - 1
                  TMP = D( K )*L( K )
                  WORK( INDLD+J ) = TMP
                  WORK( INDLLD+J ) = TMP*L( K )
                  K = K + 1
   50          CONTINUE
               IF( NDEPTH.GT.0 ) THEN
*
*  The way OLDFST and OLDLST are used in the following is a problem
*  when eigenvalues in a range are desired!
*
                  P = INDEXW( WBEGIN-1+OLDFST )
                  Q = INDEXW( WBEGIN-1+OLDLST )
                  CALL SLARRB_TGK( IN, D( IBEGIN ), L( IBEGIN ),
     $                         WORK( INDLD+1 ), WORK( INDLLD+1 ),
     $                         P, Q, RELTOL, FOUR*EPS, P-OLDFST, WORK,
     $                         WORK( INDGAP+1 ), WORK( INDERR+1 ),
     $                         WORK( INDWRK+IN ), IWORK( IINDWK ),
     $                         IINFO )
                  IF( IINFO.NE.0 ) THEN
                     INFO = 3
                     RETURN
                  END IF
               END IF
*
*              Classify eigenvalues of the current representation (RRR)
*              as (i) isolated, or (ii) clustered.
*
               NEWFRS = OLDFST
               DO 140 J = OLDFST, OLDLST
                  IF( J.EQ.OLDLST .OR. WORK( INDGAP+J ).GE.RELTOL*
     $                ABS( WORK( J ) ) ) THEN
                     NEWLST = J
                  ELSE
*
*                    continue (to the next loop)
*
                     RELGAP = WORK( INDGAP+J ) / ABS( WORK( J ) )
                     IF( J.EQ.NEWFRS ) THEN
                        MINRGP = RELGAP
                     ELSE
                        MINRGP = MIN( MINRGP, RELGAP )
                     END IF
                     GO TO 140
                  END IF
                  NEWSIZ = NEWLST - NEWFRS + 1
                  NEWFTT = WBEGIN + NEWFRS - 1
                  NOMGS = NEWSIZ.EQ.1 .OR. NEWSIZ.GT.MGSSIZ .OR.
     $                    MINRGP.LT.MGSTOL
                  IF( NEWSIZ.GT.1 .AND. NOMGS ) THEN
*
*                    Find a new L D L^T representation if the cluster
*                    size is larger than MGSSIZ or the minimum
*                    relative gap within the cluster is too small.
*                    REP='T' at the root forms the child from the TGK
*                    matrix; REP='L' forms it from the parent L D L^T.
*
                     CALL SLARRF_TGK( REP, IN, D( IBEGIN ), L( IBEGIN ),
     $                            WORK( INDLD+1 ), WORK( INDLLD+1 ),
     $                            NEWFRS, NEWLST, WORK, SIGMA,
     $                            Z( IBEGIN, NEWFTT ),
     $                            Z( IBEGIN, NEWFTT+1 ),
     $                            WORK( INDWRK ), ISEED, INFO )
                     IF( INFO.EQ.0 ) THEN
                        TMP = EPS*( ABS( SIGMA ) )
                        DO 51 K = NEWFRS, NEWLST
                           WORK( K ) = WORK( K ) - SIGMA
                           WORK( INDGAP+K ) = MAX( WORK( INDGAP+K ),
     $                                        TMP )
                           WORK( INDERR+K ) = WORK( INDERR+K ) + TMP
   51                   CONTINUE
                        NCLUS = NCLUS + 1
                        K = NEWCLS + 2*NCLUS
                        IWORK( K-1 ) = NEWFRS
                        IWORK( K ) = NEWLST
                     ELSE
                        INFO = 0
                        IF( MINRGP.GE.MGSTOL ) THEN
*                          No relatively robust child representation
*                          could be found for this cluster; report it
*                          instead of leaving the cluster unprocessed
*                          (the caller falls back to SBDSDC).
                           INFO = 3
                           RETURN
                        ELSE IF( REP.EQ.'T' ) THEN
*                          At the root the representation is the
*                          Golub-Kahan matrix by its entries, so the
*                          LDL^T-based inverse-iteration fallback does
*                          not apply.
                           INFO = 2
                           RETURN
                        ELSE
*
*                          Call SSTEIN to process this tight cluster.
*                          This happens only if MINRGP <= MGSTOL
*                          and SLARRF_TGK returns INFO = 1. The latter
*                          means that a new RRR to "break" the
*                          cluster could not be found.
*
                           WORK( INDWRK ) = D( IBEGIN )
                           DO 60 K = 1, IN - 1
                              WORK( INDWRK+K ) = D( IBEGIN+K ) +
     $                                           WORK( INDLLD+K )
   60                      CONTINUE
                           DO 70 K = 1, NEWSIZ
                              IWORK( IINDWK+K-1 ) = 1
   70                      CONTINUE
                           DO 80 K = NEWFRS, NEWLST
                              ISUPPZ( 2*( OLDIEN+K )-1 ) = 1
                              ISUPPZ( 2*( OLDIEN+K ) ) = IN
   80                      CONTINUE
                           TEMP( 1 ) = IN
                           CALL SSTEIN( IN, WORK( INDWRK ),
     $                                  WORK( INDLD+1 ), NEWSIZ,
     $                                  WORK( NEWFRS ),
     $                                  IWORK( IINDWK ), TEMP( 1 ),
     $                                  Z( IBEGIN, NEWFTT ), LDZ,
     $                                  WORK( INDWRK+IN ),
     $                                  IWORK( IINDWK+IN ),
     $                                  IWORK( IINDWK+2*IN ), IINFO )
                           IF( IINFO.NE.0 ) THEN
                              INFO = 2
                              RETURN
                           END IF
                           NDONE = NDONE + NEWSIZ
                        END IF
                     END IF
                  ELSE IF( NEWSIZ.GT.1 ) THEN
                     KTOT = NEWFTT
                     DO 100 K = NEWFRS, NEWLST
                        ITER = 0
   90                   CONTINUE
                        LAMBDA = WORK( K )
*
*                       Given LAMBDA, compute the eigenvector.  REP='T'
*                       uses the TGK matrix as input at the root.
*
                        CALL SLAR1V_TGK( REP, IN, 1, IN, LAMBDA,
     $                               D( IBEGIN ),
     $                               L( IBEGIN ), WORK( INDLD+1 ),
     $                               WORK( INDLLD+1 ), W( WBEGIN+K-1 ),
     $                               GERSCH( 2*OLDIEN+1 ),
     $                               Z( IBEGIN, KTOT ), ZTZ, MINGMA,
     $                               IWORK( IINDR+KTOT ),
     $                               ISUPPZ( 2*KTOT-1 ),
     $                               WORK( INDWRK ) )
                        TMP = ONE / ZTZ
                        NRMINV = SQRT( TMP )
                        RESID = ABS( MINGMA )*NRMINV
                        RQCORR = MINGMA*TMP
                        IF( K.EQ.IN ) THEN
                           GAP = WORK( INDGAP+K-1 )
                        ELSE IF( K.EQ.1 ) THEN
                           GAP = WORK( INDGAP+K )
                        ELSE
                           GAP = MIN( WORK( INDGAP+K-1 ),
     $                           WORK( INDGAP+K ) )
                        END IF
                        ITER = ITER + 1
                        IF( RESID.GT.TOL*GAP .AND. ABS( RQCORR ).GT.
     $                      FOUR*EPS*ABS( LAMBDA ) ) THEN
                           WORK( K ) = LAMBDA + RQCORR
                           IF( ITER.LT.MAXITR ) THEN
                              GO TO 90
                           END IF
                        END IF
                        IWORK( KTOT ) = 1
                        IF( NEWSIZ.EQ.1 )
     $                     NDONE = NDONE + 1
                        ZFROM = ISUPPZ( 2*KTOT-1 )
                        ZTO = ISUPPZ( 2*KTOT )
                        CALL SSCAL( ZTO-ZFROM+1, NRMINV,
     $                          Z( IBEGIN+ZFROM-1, KTOT ), 1 )
                        KTOT = KTOT + 1
  100                CONTINUE
                     IF( NEWSIZ.GT.1 ) THEN
                        ITMP1 = ISUPPZ( 2*NEWFTT-1 )
                        ITMP2 = ISUPPZ( 2*NEWFTT )
                        KTOT = OLDIEN + NEWLST
                        DO 120 P = NEWFTT + 1, KTOT
                           DO 110 Q = NEWFTT, P - 1
                              TMP = -SDOT( IN, Z( IBEGIN, P ), 1,
     $                               Z( IBEGIN, Q ), 1 )
                              CALL SAXPY( IN, TMP, Z( IBEGIN, Q ), 1,
     $                                    Z( IBEGIN, P ), 1 )
  110                      CONTINUE
                           TMP = ONE / SNRM2( IN, Z( IBEGIN, P ), 1 )
                           CALL SSCAL( IN, TMP, Z( IBEGIN, P ), 1 )
                           ITMP1 = MIN( ITMP1, ISUPPZ( 2*P-1 ) )
                           ITMP2 = MAX( ITMP2, ISUPPZ( 2*P ) )
  120                   CONTINUE
                        DO 130 P = NEWFTT, KTOT
                           ISUPPZ( 2*P-1 ) = ITMP1
                           ISUPPZ( 2*P ) = ITMP2
  130                   CONTINUE
                        NDONE = NDONE + NEWSIZ
                     END IF
                  ELSE
*                    Defer isolated eigenvectors until this complete parent
*                    RRR has been classified.  Consecutive marked vectors
*                    can then run through two independent SLAR1V recurrences
*                    together without changing either lane's FP order.
                     IWORK( NEWFTT ) = -1
                  END IF
                  NEWFRS = J + 1
  140          CONTINUE
*
*              Process the isolated eigenvectors belonging to this parent
*              representation.  Pair adjacent lanes when possible; an odd
*              tail follows the original scalar path verbatim.  The second
*              4*IN scratch lane ends at WORK(13*N), exactly within the
*              private SLARRV_TGK workspace supplied by SBDSVDMR3.
*
               K = OLDFST
  141          CONTINUE
               IF( K.LE.OLDLST ) THEN
                  KTOT = WBEGIN + K - 1
                  IF( IWORK( KTOT ).NE.-1 ) THEN
                     K = K + 1
                     GO TO 141
                  END IF
                  K2 = K + 1
                  IF( K2.LE.OLDLST ) THEN
                     KTOT2 = WBEGIN + K2 - 1
                  ELSE
                     KTOT2 = 0
                  END IF
                  KTOT2OK = .FALSE.
                  IF( KTOT2.GT.0 ) THEN
                     IF( IWORK( KTOT2 ).EQ.-1 ) KTOT2OK = .TRUE.
                  END IF
                  IF( KTOT2OK ) THEN
                     ITER = 0
                     ITER2 = 0
                     DONE1 = .FALSE.
                     DONE2 = .FALSE.
  142                CONTINUE
                     LAMBDA = WORK( K )
                     LAMBDA2 = WORK( K2 )
*                    The paired two-lane variant SLAR1V2_TGK is kept in
*                    slar1v_tgk.f but not used: the two eigenvectors are
*                    computed by separate one-lane calls (reference code).
                     IF( .NOT.DONE1 ) THEN
                        CALL SLAR1V_TGK( REP, IN, 1, IN, LAMBDA,
     $                       D( IBEGIN ), L( IBEGIN ), WORK( INDLD+1 ),
     $                       WORK( INDLLD+1 ), W( WBEGIN+K-1 ),
     $                       GERSCH( 2*OLDIEN+1 ), Z( IBEGIN, KTOT ),
     $                       ZTZ, MINGMA, IWORK( IINDR+KTOT ),
     $                       ISUPPZ( 2*KTOT-1 ), WORK( INDWRK ) )
                     END IF
                     IF( .NOT.DONE2 ) THEN
                        CALL SLAR1V_TGK( REP, IN, 1, IN, LAMBDA2,
     $                       D( IBEGIN ), L( IBEGIN ), WORK( INDLD+1 ),
     $                       WORK( INDLLD+1 ), W( WBEGIN+K2-1 ),
     $                       GERSCH( 2*OLDIEN+1 ), Z( IBEGIN, KTOT2 ),
     $                       ZTZ2, MINGMA2, IWORK( IINDR+KTOT2 ),
     $                       ISUPPZ( 2*KTOT2-1 ),
     $                       WORK( INDWRK+4*IN ) )
                     END IF
                     IF( .NOT.DONE1 ) THEN
                        TMP = ONE / ZTZ
                        NRMINV = SQRT( TMP )
                        RESID = ABS( MINGMA )*NRMINV
                        RQCORR = MINGMA*TMP
                        IF( K.EQ.IN ) THEN
                           GAP = WORK( INDGAP+K-1 )
                        ELSE IF( K.EQ.1 ) THEN
                           GAP = WORK( INDGAP+K )
                        ELSE
                           GAP = MIN( WORK( INDGAP+K-1 ),
     $                           WORK( INDGAP+K ) )
                        END IF
                        ITER = ITER + 1
                        DONE1 = .TRUE.
                        IF( RESID.GT.TOL*GAP .AND. ABS( RQCORR ).GT.
     $                      FOUR*EPS*ABS( LAMBDA ) ) THEN
                           WORK( K ) = LAMBDA + RQCORR
                           IF( ITER.LT.MAXITR ) DONE1 = .FALSE.
                        END IF
                     END IF
                     IF( .NOT.DONE2 ) THEN
                        TMP2 = ONE / ZTZ2
                        NRMINV2 = SQRT( TMP2 )
                        RESID2 = ABS( MINGMA2 )*NRMINV2
                        RQCORR2 = MINGMA2*TMP2
                        IF( K2.EQ.IN ) THEN
                           GAP2 = WORK( INDGAP+K2-1 )
                        ELSE IF( K2.EQ.1 ) THEN
                           GAP2 = WORK( INDGAP+K2 )
                        ELSE
                           GAP2 = MIN( WORK( INDGAP+K2-1 ),
     $                            WORK( INDGAP+K2 ) )
                        END IF
                        ITER2 = ITER2 + 1
                        DONE2 = .TRUE.
                        IF( RESID2.GT.TOL*GAP2 .AND. ABS( RQCORR2 ).GT.
     $                      FOUR*EPS*ABS( LAMBDA2 ) ) THEN
                           WORK( K2 ) = LAMBDA2 + RQCORR2
                           IF( ITER2.LT.MAXITR ) DONE2 = .FALSE.
                        END IF
                     END IF
                     IF( .NOT.DONE1 .OR. .NOT.DONE2 ) GO TO 142
                     IWORK( KTOT ) = 1
                     NDONE = NDONE + 1
                     ZFROM = ISUPPZ( 2*KTOT-1 )
                     ZTO = ISUPPZ( 2*KTOT )
                     CALL SSCAL( ZTO-ZFROM+1, NRMINV,
     $                    Z( IBEGIN+ZFROM-1, KTOT ), 1 )
                     IWORK( KTOT2 ) = 1
                     NDONE = NDONE + 1
                     ZFROM = ISUPPZ( 2*KTOT2-1 )
                     ZTO = ISUPPZ( 2*KTOT2 )
                     CALL SSCAL( ZTO-ZFROM+1, NRMINV2,
     $                    Z( IBEGIN+ZFROM-1, KTOT2 ), 1 )
                     K = K + 2
                  ELSE
                     ITER = 0
  143                CONTINUE
                     LAMBDA = WORK( K )
                     CALL SLAR1V_TGK( REP, IN, 1, IN, LAMBDA,
     $                    D( IBEGIN ), L( IBEGIN ), WORK( INDLD+1 ),
     $                    WORK( INDLLD+1 ), W( WBEGIN+K-1 ),
     $                    GERSCH( 2*OLDIEN+1 ), Z( IBEGIN, KTOT ),
     $                    ZTZ, MINGMA, IWORK( IINDR+KTOT ),
     $                    ISUPPZ( 2*KTOT-1 ), WORK( INDWRK ) )
                     TMP = ONE / ZTZ
                     NRMINV = SQRT( TMP )
                     RESID = ABS( MINGMA )*NRMINV
                     RQCORR = MINGMA*TMP
                     IF( K.EQ.IN ) THEN
                        GAP = WORK( INDGAP+K-1 )
                     ELSE IF( K.EQ.1 ) THEN
                        GAP = WORK( INDGAP+K )
                     ELSE
                        GAP = MIN( WORK( INDGAP+K-1 ),
     $                        WORK( INDGAP+K ) )
                     END IF
                     ITER = ITER + 1
                     IF( RESID.GT.TOL*GAP .AND. ABS( RQCORR ).GT.
     $                   FOUR*EPS*ABS( LAMBDA ) ) THEN
                        WORK( K ) = LAMBDA + RQCORR
                        IF( ITER.LT.MAXITR ) GO TO 143
                     END IF
                     IWORK( KTOT ) = 1
                     NDONE = NDONE + 1
                     ZFROM = ISUPPZ( 2*KTOT-1 )
                     ZTO = ISUPPZ( 2*KTOT )
                     CALL SSCAL( ZTO-ZFROM+1, NRMINV,
     $                    Z( IBEGIN+ZFROM-1, KTOT ), 1 )
                     K = K + 1
                  END IF
                  GO TO 141
               END IF
  150       CONTINUE
            NDEPTH = NDEPTH + 1
            GO TO 40
         END IF
         J = 2*WBEGIN
         DO 160 I = WBEGIN, WEND
            ISUPPZ( J-1 ) = ISUPPZ( J-1 ) + OLDIEN
            ISUPPZ( J ) = ISUPPZ( J ) + OLDIEN
            J = J + 2
  160    CONTINUE
         IBEGIN = IEND + 1
         WBEGIN = WEND + 1
  170 CONTINUE
*
      RETURN
*
*     End of SLARRV_TGK
*
      END
