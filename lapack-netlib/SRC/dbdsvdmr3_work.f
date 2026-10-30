*> \brief \b DBDSVDMR3_WORK computational core of DBDSVDMR3.
*>
*> \ingroup bdsvdmr3
*>
*> \par Purpose:
*> =============
*>
*> \verbatim
*> computational core of DBDSVDMR3: preprocessing, dqds singular values, MRRR eigenvectors of the Golub-Kahan matrix, extraction of singular vectors.
*> All routines in this file are auxiliary to DBDSVDMR3 and derive
*> from the LAPACK 3.0 (1999) MRRR kernels DLARRV, DLAR1V, DLARRB and
*> DLARRF rather than from the current DSTEMR kernels.
*> \endverbatim
*
      SUBROUTINE DBDSVDMR3_WORK( JOBZ, UPLO, N, D, E, S, U, LDU, VT,
     $              LDVT,
     $                   M, WORK, LWORK, IWORK, LIWORK, INFO )
*
*  -- New driver for the bidiagonal SVD via TGK-rooted MR^3 --
*     Uses the LAPACK 3.0-lineage MRRR tree code below the root via
*     DLARRV_TGK.  The bundled DLARRB has a bounded-progress guard so an
*     invalid bracket radius returns INFO instead of looping forever.
*     Standard LAPACK dependencies include DLASQ1, DLAMCH, DLANST,
*     DSCAL, and DSWAP.
*
*  Purpose
*  =======
*
*  DBDSVDMR3 computes the SVD of a real N-by-N bidiagonal matrix
*  B = U * S * VT, with S diagonal (the singular values, ascending) and
*  U, VT orthonormal.  Method:
*    1. Build the Tridiagonal Golub-Kahan (TGK) matrix of order 2N; its
*       off-diagonals beta interleave the bidiagonal entries (d,e).
*    2. Split the TGK wherever any beta is negligible -- this covers both
*       the bidiagonal off-diagonals e_i (even positions) and zero
*       diagonal entries d_i (odd positions), each of which makes the TGK
*       reducible.
*    3. Per TGK block, compute the positive eigenvalues (= singular
*       values) with DLASQ1 (dqds); an odd-dimension block is a
*       rectangular bidiagonal handled by appending a zero diagonal.
*    4. Compute the positive eigenvectors by MR^3 rooted at the TGK
*       (DLARRV_TGK) and de-interleave each into its left/right singular
*       vectors.
*    5. Each odd-dimension block also carries one structural ZERO
*       eigenvalue; its one-sided null vector (closed form) supplies a
*       singular vector for a zero singular value.
*
*  Arguments  (see below; UPLO='U'/'L'; JOBZ='N'/'V'; RANGE = all)
*  =========
*
*  JOBZ   (input) CHARACTER*1   'N' values only; 'V' values + vectors.
*  UPLO   (input) CHARACTER*1   'U' upper bidiagonal; 'L' lower.
*  N      (input) INTEGER       Order of B.  N >= 0.
*  D      (in/out) DOUBLE (N)   Diagonal of B (overwritten).
*  E      (in/out) DOUBLE (N)   Off-diagonal in 1..N-1 (overwritten).
*  S      (output) DOUBLE (N)   Singular values, ascending.
*  U      (output) DOUBLE (LDU,N)   Left singular vectors (col i <-> S(i)).
*  LDU    (input) INTEGER
*  VT     (output) DOUBLE (LDVT,N)  Right singular vectors, by ROWS.
*  LDVT   (input) INTEGER
*  M      (output) INTEGER      Number of singular values found (= N).
*  WORK   (workspace) DOUBLE (LWORK).  LWORK >= max(1, 2*N*N + 34*N).
*  LWORK  (input) INTEGER       -1 => workspace query.
*  IWORK  (workspace) INTEGER (LIWORK).  LIWORK >= max(1, 20*N).
*  LIWORK (input) INTEGER       -1 => workspace query.
*  INFO   (output) INTEGER  0 ok; <0 illegal arg; 1 DLASQ1; 2 DLARRV_TGK.
*
*  =====================================================================
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDU, LDVT, LIWORK, LWORK, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), S( * ), U( LDU, * ),
     $                   VT( LDVT, * ), WORK( * )
*     ..
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
*     ..
*     .. Local Arrays ..
*     ..
*     .. Local Scalars ..
      LOGICAL            DIAGB, WANTZ, LQUERY, UPPER
      INTEGER            I, IDBL, IEBL, IINFO, INDBL, INDE2, INDGSC,
     $                   INDISP, INDISUP, INDIWK, INDIXW, INDTGK, INDWK,
     $                   IQ, J, JBLK, K, KB, KZU, KZV, L, LIWMIN, LWMIN,
     $                   MPOS, NPOS, NSPLIT, NZSV, TBEG, TEND, BDIM
      DOUBLE PRECISION   BNRM, EPS, RT2, SCL, SFMIN, THRESH, TOL,
     $                   ZNRM, ZE, ZG
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           LSAME, DLAMCH, DLANST
*     ..
*     .. External Subroutines ..
      EXTERNAL           DBDTGK, DLARRV_TGK, DLASQ1, DSCAL,
     $                   DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MOD, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
      LWMIN = MAX( 1, 2*N*N + 34*N )
      LIWMIN = MAX( 1, 24*N )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDU.LT.1 .OR. ( WANTZ .AND. LDU.LT.N ) ) THEN
         INFO = -8
      ELSE IF( LDVT.LT.1 .OR. ( WANTZ .AND. LDVT.LT.N ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -13
      ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -15
      END IF
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DBDSVDMR3_WORK', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      M = 0
      IF( N.EQ.0 ) RETURN
      IF( N.EQ.1 ) THEN
         M = 1
         S( 1 ) = ABS( D( 1 ) )
         IF( WANTZ ) THEN
            U( 1, 1 ) = SIGN( ONE, D( 1 ) )
            VT( 1, 1 ) = ONE
         END IF
         RETURN
      END IF
*
*     An exactly diagonal bidiagonal has an exact canonical SVD.  Taking
*     it through the TGK path would normalize 2-by-2 eigenvectors by
*     1/sqrt(2) and then multiply them by sqrt(2), needlessly turning
*     exact identity/permutation vectors into one-ULP approximations.
*     Besides removing that roundoff, this path avoids the MR^3 setup for
*     the diagonal cases that occur in the LAPACK test suite.
*
      DIAGB = .TRUE.
      DO 5 I = 1, N - 1
         IF( E( I ).NE.ZERO ) DIAGB = .FALSE.
    5 CONTINUE
      IF( DIAGB ) THEN
         DO 8 I = 1, N
            S( I ) = ABS( D( I ) )
    8    CONTINUE
         IF( WANTZ ) THEN
            DO 12 J = 1, N
               DO 10 I = 1, N
                  U( I, J ) = ZERO
                  VT( I, J ) = ZERO
   10          CONTINUE
               U( J, J ) = SIGN( ONE, D( J ) )
               VT( J, J ) = ONE
   12       CONTINUE
         END IF
         CALL DBSORT( N, S, U, LDU, VT, LDVT, WANTZ )
         M = N
         RETURN
      END IF
*
      EPS = DLAMCH( 'Precision' )
      SFMIN = DLAMCH( 'Safe minimum' )
      RT2 = SQRT( TWO )
*
*     Workspace map (DOUBLE):  INDTGK 2N (TGK diag=0); INDE2 2N (TGK
*     off-diagonals beta); INDGSC 4N (Gerschgorin); INDWK = scratch /
*     eigenvector array Z (LDZ=2N) / DLARRV_TGK WORK.
*     INTEGER:  INDISP 2N (ISPLIT); INDBL 2N (IBLOCK); INDIXW 2N (INDEXW);
*     INDISUP 2N (ISUPPZ scratch); INDIWK (DLARRV_TGK IWORK).
*
      INDTGK = 1
      INDE2 = INDTGK + 2*N
      INDGSC = INDE2 + 2*N
      INDWK = INDGSC + 4*N
*
      INDISP = 1
      INDBL = INDISP + 2*N
      INDIXW = INDBL + 2*N
      INDISUP = INDIXW + 2*N
      INDIWK = INDISUP + 2*N
*
*     Scale B so that its norm is in a safe range.
*
      BNRM = DLANST( 'M', N, D, E )
      SCL = ONE
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SQRT( SFMIN ) ) THEN
         SCL = SQRT( SFMIN ) / BNRM
      ELSE IF( BNRM.GT.ONE / SQRT( SFMIN ) ) THEN
         SCL = ( ONE / SQRT( SFMIN ) ) / BNRM
      END IF
      IF( SCL.NE.ONE ) THEN
         CALL DSCAL( N, SCL, D, 1 )
         CALL DSCAL( N-1, SCL, E, 1 )
         BNRM = BNRM*SCL
      END IF
*
*     Random relative perturbation (Dhillon, Parlett & Voemel 2005,
*     sec. 5.4) is NOT applied blindly to the root here.  A blanket
*     O(n*eps) perturbation jolts every singular value by ~n*eps, which
*     OVERWRITES the genuine fine structure of clusters whose true
*     internal gaps are smaller than n*eps (e.g. the Lipshitz matrices,
*     whose ~800-fold cluster at sigma~1 has real relative spacing
*     ~7e-15 < n*eps): there the perturbation destroyed orthogonality.
*     Instead the perturbation is applied SELECTIVELY by DLARRF_TGK, and
*     only to a child RRR whose cluster contains EXACTLY-equal
*     eigenvalues (true ties, e.g. from identical glued blocks) that
*     recursive shifting alone can never separate; the fixed seed is
*     set inside DLARRV_TGK.
*
*
*     Build the TGK matrix of order 2N from the (perturbed) bidiagonal.
*     WORK(INDE2 .. INDE2+2N-2) holds beta_1..beta_{2N-1}, where
*     beta_{2k-1}=D(k), beta_{2k}=E(k).
*
      CALL DBDTGK( UPLO, N, D, E, WORK( INDTGK ), WORK( INDE2 ),
     $             WORK( INDGSC ) )
*
*     Split the TGK wherever an off-diagonal beta_j is negligible.  This
*     handles both the bidiagonal off-diagonals e_i (even j) and zero
*     diagonal entries d_i (odd j).  ISPLIT lists the block-end TGK rows.
*
      THRESH = EPS*BNRM
      NSPLIT = 0
      DO 20 I = 1, 2*N - 1
         IF( ABS( WORK( INDE2+I-1 ) ).LE.THRESH ) THEN
            WORK( INDE2+I-1 ) = ZERO
            NSPLIT = NSPLIT + 1
            IWORK( INDISP+NSPLIT-1 ) = I
         END IF
   20 CONTINUE
      NSPLIT = NSPLIT + 1
      IWORK( INDISP+NSPLIT-1 ) = 2*N
*
*     For each TGK block compute the NPOS positive eigenvalues
*     (= singular values) with DLASQ1.  Block of dimension BDIM has its
*     beta's at WORK(INDE2 + TBEG-1 .. INDE2 + TEND-2); the sub-bidiagonal
*     is d~(l)=beta(local 2l-1), e~(l)=beta(local 2l).
*        even BDIM=2k : square k-by-k bidiagonal -> k positive sv;
*        odd  BDIM=2k+1: rectangular k-by-(k+1); append a zero diagonal to
*                        make it (k+1)-square -> k positive sv + 1 zero.
*     The NPOS positives sit in the upper part of the block's BDIM-point
*     spectrum, so INDEXW(j-th positive) = (BDIM-NPOS)+j.
*
      IDBL = INDWK
      IEBL = INDWK + N + 1
      IQ = INDWK + 2*N + 2
      MPOS = 0
      TBEG = 1
      DO 60 JBLK = 1, NSPLIT
         TEND = IWORK( INDISP+JBLK-1 )
         BDIM = TEND - TBEG + 1
         NPOS = BDIM / 2
         KB = NPOS
*        diagonals d~(1..KB) at local odd positions 1,3,...
         DO 30 L = 1, KB
            WORK( IDBL+L-1 ) = ABS( WORK( INDE2+TBEG-2+2*L-1 ) )
   30    CONTINUE
         IF( MOD( BDIM, 2 ).EQ.0 ) THEN
*           even block: square KB-by-KB, KB-1 off-diagonals
            DO 32 L = 1, KB - 1
               WORK( IEBL+L-1 ) = ABS( WORK( INDE2+TBEG-2+2*L ) )
   32       CONTINUE
            CALL DLASQ1( KB, WORK( IDBL ), WORK( IEBL ), WORK( IQ ),
     $                   IINFO )
         ELSE
*           odd block: append a zero diagonal -> (KB+1)-square, KB offdiag
            WORK( IDBL+KB ) = ZERO
            DO 34 L = 1, KB
               WORK( IEBL+L-1 ) = ABS( WORK( INDE2+TBEG-2+2*L ) )
   34       CONTINUE
            CALL DLASQ1( KB+1, WORK( IDBL ), WORK( IEBL ), WORK( IQ ),
     $                   IINFO )
         END IF
         IF( IINFO.NE.0 ) THEN
            INFO = 1
            RETURN
         END IF
*        DLASQ1 returns the singular values in DECREASING order in IDBL;
*        the appended structural zero (odd block) is the smallest entry,
*        so the NPOS positives are IDBL(0..NPOS-1).  Store ascending.
         DO 40 I = 1, NPOS
            S( MPOS+I ) = WORK( IDBL+NPOS-I )
            IWORK( INDBL+MPOS+I-1 ) = JBLK
            IWORK( INDIXW+MPOS+I-1 ) = ( BDIM-NPOS ) + I
   40    CONTINUE
         MPOS = MPOS + NPOS
         TBEG = TEND + 1
   60 CONTINUE
*
      NZSV = N - MPOS
*
*     Eigenvectors are not wanted: append the NZSV zero singular values,
*     sort ascending, rescale, and return.
*
      IF( .NOT.WANTZ ) THEN
         DO 65 I = 1, NZSV
            S( MPOS+I ) = ZERO
   65    CONTINUE
         CALL DBSORT( N, S, U, LDU, VT, LDVT, .FALSE. )
         IF( SCL.NE.ONE )
     $      CALL DSCAL( N, ONE / SCL, S, 1 )
         M = N
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
         RETURN
      END IF
*
*     NOTE: this must run BEFORE DLARRV_TGK, which overwrites the
*     off-diagonal array WORK(INDE2) in place with child representations.
*     Zero singular values.  Each odd-dimension TGK block has one
*     structural zero eigenvalue whose null vector lives on the ODD local
*     positions of the block (the even local entries are zero):
*        z(1)=1,  z(2l+1) = -beta(2l-1)/beta(2l) * z(2l-1).
*     De-interleaving by global parity, the vector is one-sided -- it is
*     a left singular vector (u) if the block starts at a u-row, else a
*     right singular vector (v) -- supplying one half of a zero singular
*     triplet (the other half comes from the paired block).  Because the
*     singular value is zero, any orthonormal completion is valid, so we
*     simply collect the u's into the trailing U-columns and the v's into
*     the trailing VT-rows.  Build the 2N null vector in WORK(INDWK) (the
*     eigenvector array, not yet in use).
*
      KZU = 0
      KZV = 0
      TBEG = 1
      DO 160 JBLK = 1, NSPLIT
         TEND = IWORK( INDISP+JBLK-1 )
         BDIM = TEND - TBEG + 1
         IF( MOD( BDIM, 2 ).EQ.1 ) THEN
*           build the local null vector spread over global rows TBEG..TEND
            DO 110 I = 1, 2*N
               WORK( INDWK+I-1 ) = ZERO
  110       CONTINUE
            WORK( INDWK+TBEG-1 ) = ONE
            ZNRM = ONE
            ZG = ONE
            DO 120 L = 1, ( BDIM-1 ) / 2
               ZE = -( WORK( INDE2+TBEG-2+2*L-1 ) /
     $                 WORK( INDE2+TBEG-2+2*L ) )*ZG
               WORK( INDWK+TBEG-1+2*L ) = ZE
               ZNRM = ZNRM + ZE*ZE
               ZG = ZE
  120       CONTINUE
            ZNRM = SQRT( ZNRM )
*           The (one-sided) null vector lives on global-ODD rows when the
*           block starts at an odd row (MOD(TBEG,2)=1), else on global-
*           EVEN rows.  For LOWER B the odd half is u and the even half is
*           v; for UPPER the two roles swap.  Route the nonzero half to a
*           u-column or a v-row accordingly, normalizing by its own norm.
            IF( MOD( TBEG, 2 ).EQ.1 ) THEN
               IF( UPPER ) THEN
                  KZV = KZV + 1
                  DO 140 I = 1, N
                     VT( MPOS+KZV, I ) = WORK( INDWK+2*I-2 ) / ZNRM
  140             CONTINUE
               ELSE
                  KZU = KZU + 1
                  DO 142 I = 1, N
                     U( I, MPOS+KZU ) = WORK( INDWK+2*I-2 ) / ZNRM
  142             CONTINUE
               END IF
            ELSE
               IF( UPPER ) THEN
                  KZU = KZU + 1
                  DO 144 I = 1, N
                     U( I, MPOS+KZU ) = WORK( INDWK+2*I-1 ) / ZNRM
  144             CONTINUE
               ELSE
                  KZV = KZV + 1
                  DO 146 I = 1, N
                     VT( MPOS+KZV, I ) = WORK( INDWK+2*I-1 ) / ZNRM
  146             CONTINUE
               END IF
            END IF
         END IF
         TBEG = TEND + 1
  160 CONTINUE
*     Compute the eigenvectors of the MPOS positive eigenvalues by MR^3
*     rooted at the TGK matrix.  The 2N-by-MPOS eigenvector array is held
*     in WORK(INDWK) with leading dimension 2N.
*
      IF( MPOS.GT.0 ) THEN
         TOL = MAX( DBLE( 2*N )*EPS, ZERO )
         CALL DLARRV_TGK( 2*N, WORK( INDTGK ), WORK( INDE2 ),
     $                IWORK( INDISP ), MPOS, S, IWORK( INDBL ),
     $                IWORK( INDIXW ), WORK( INDGSC ), TOL,
     $                WORK( INDWK ), 2*N, IWORK( INDISUP ),
     $                WORK( INDWK+2*N*MPOS ), IWORK( INDIWK ), IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = 2
            RETURN
         END IF
*
*        De-interleave each 2N-eigenvector into its left/right singular
*        vectors by GLOBAL position parity (this stays correct across
*        split blocks).  beta_{2k-1}=d_k, beta_{2k}=e_k, so global-odd
*        rows carry one half and global-even rows the other.  For LOWER B
*        the odd half is u and the even half is v; for UPPER (B^T is the
*        lower bidiagonal with the same d,e) the two roles swap.
*
         DO 100 K = 1, MPOS
            DO 90 I = 1, N
               IF( UPPER ) THEN
                  VT( K, I ) = RT2*WORK( INDWK+(K-1)*2*N+(2*I-2) )
                  U( I, K )  = RT2*WORK( INDWK+(K-1)*2*N+(2*I-1) )
               ELSE
                  U( I, K )  = RT2*WORK( INDWK+(K-1)*2*N+(2*I-2) )
                  VT( K, I ) = RT2*WORK( INDWK+(K-1)*2*N+(2*I-1) )
               END IF
   90       CONTINUE
  100    CONTINUE
      END IF
*
*
*     Set the zero singular values and (defensively) any unmatched
*     trailing vectors, then sort everything ascending and rescale.
*
      DO 170 I = 1, NZSV
         S( MPOS+I ) = ZERO
  170 CONTINUE
      M = N
      CALL DBSORT( N, S, U, LDU, VT, LDVT, .TRUE. )
      IF( SCL.NE.ONE )
     $   CALL DSCAL( N, ONE / SCL, S, 1 )
*
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
      RETURN
*
*     End of DBDSVDMR3
*
      END
*
*     =================================================================
*
      SUBROUTINE DBSORT( N, S, U, LDU, VT, LDVT, WANTZ )
*
*     Sort S(1:N) into ascending order by straight selection, applying
*     the same permutation to the columns of U and the rows of VT when
*     WANTZ is .TRUE.
*
      LOGICAL            WANTZ
      INTEGER            N, LDU, LDVT
      DOUBLE PRECISION   S( * ), U( LDU, * ), VT( LDVT, * )
      INTEGER            I, J, K
      DOUBLE PRECISION   P
      EXTERNAL           DSWAP
      INTRINSIC          ABS
      DO 20 I = 1, N - 1
         K = I
         P = S( I )
         DO 10 J = I + 1, N
            IF( S( J ).LT.P ) THEN
               K = J
               P = S( J )
            END IF
   10    CONTINUE
         IF( K.NE.I ) THEN
            S( K ) = S( I )
            S( I ) = P
            IF( WANTZ ) THEN
               CALL DSWAP( N, U( 1, I ), 1, U( 1, K ), 1 )
               CALL DSWAP( N, VT( I, 1 ), LDVT, VT( K, 1 ), LDVT )
            END IF
         END IF
   20 CONTINUE
      RETURN
      END
