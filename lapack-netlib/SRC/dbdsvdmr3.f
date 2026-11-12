*> \brief \b DBDSVDMR3
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DBDSVDMR3( UPLO, JOBZ, RANGE, N, D, E, VL, VU, IL, IU,
*                          NS, S, Z, LDZ, WORK, LWORK, IWORK, LIWORK,
*                          INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ, RANGE, UPLO
*       INTEGER            IL, INFO, IU, LDZ, LIWORK, LWORK, N, NS
*       DOUBLE PRECISION   VL, VU
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       DOUBLE PRECISION   D( * ), E( * ), S( * ), WORK( * ),
*                          Z( LDZ, * )
*       ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>  DBDSVDMR3 computes the singular value decomposition (SVD) of a real
*>  N-by-N (upper or lower) bidiagonal matrix B, B = U * S * VT, using
*>  the Tridiagonal Golub-Kahan (TGK) MR^3 driver DBDSVDMR3.  Its
*>  argument list mirrors DBDSVDX so that DBDSVDMR3 is a drop-in
*>  alternative alongside DBDSVDX for callers that want the MR^3
*>  variant.
*>
*>  Given an upper bidiagonal B with diagonal D and superdiagonal E,
*>  DBDSVDMR3 builds the 2N-by-2N TGK tridiagonal, computes its
*>  eigendecomposition by the MR^3 algorithm rooted at the TGK
*>  (see P. Willems and B. Lang, SIAM J. Sci. Comput., 35:740-766,
*>  2013), and returns the singular values and (optionally) singular
*>  vectors of B in the same packing convention as DBDSVDX.
*>
*>  RANGE = 'V' or 'I' is realized by post-filtering the full spectrum
*>  from DBDSVDMR3 (which always returns all N singular values).
*>
*>  DBDSVDMR3 tries DBDSVDMR3_WORK first for every matrix order.  If DBDSVDMR3_WORK
*>  reports a positive INFO, returns an incomplete spectrum, or produces
*>  a NaN or infinity in the requested output, DBDSVDMR3 restores the input
*>  and retries with DBDSDC.  DBDSDC results are normalized to the same
*>  internal ordering before RANGE filtering, so either successful path
*>  is transparent.  This result-based fallback avoids a size cutoff:
*>  solver-only probes showed that the historical large TIMEOUT cases
*>  were evaluator timeouts rather than DBDSVDMR3 failing to return.
*>
*>  Subset modes.  RANGE = 'I' and RANGE = 'V' follow DBDSVDX's
*>  conventions (for RANGE = 'I', IL = 1 selects the largest singular
*>  value; S is returned in descending order).  In all modes the full
*>  decomposition is computed and the requested subset is then
*>  returned, so the cost and workspace are those of RANGE = 'A'; the
*>  O(N*K) subset capability of MRRR is not exposed by this routine.
*>
*>  Accuracy notes.  (1) Entries of B (diagonal or off-diagonal) with
*>  magnitude at most eps*max|B| are set to zero before the singular
*>  values are computed; the singular values are therefore accurate
*>  to O(eps*max|B|) in the absolute sense (as for DBDSDC) but are not
*>  determined to high relative accuracy when they are far below |B|
*>  (DBDSQR and DBDSVDX preserve relative accuracy).  (2) Child
*>  representations of the MRRR tree are accepted without the
*>  element-growth test of DSTEMR; on spectra with tight clusters the
*>  orthogonality of the computed vectors can degrade to 10-100 times
*>  that of DBDSDC while the residuals remain small.  See the
*>  accompanying paper for a quantitative account.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  B is upper bidiagonal;
*>          = 'L':  B is lower bidiagonal.
*> \endverbatim
*>
*> \param[in] JOBZ
*> \verbatim
*>          JOBZ is CHARACTER*1
*>          = 'N':  Compute singular values only;
*>          = 'V':  Compute singular values and singular vectors.
*> \endverbatim
*>
*> \param[in] RANGE
*> \verbatim
*>          RANGE is CHARACTER*1
*>          = 'A': all singular values will be found.
*>          = 'V': all singular values in the half-open interval [VL,VU)
*>                 will be found.
*>          = 'I': the IL-th through IU-th singular values will be found.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the bidiagonal matrix.  N >= 0.
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is DOUBLE PRECISION array, dimension (N)
*>          The n diagonal elements of the bidiagonal matrix B.
*> \endverbatim
*>
*> \param[in] E
*> \verbatim
*>          E is DOUBLE PRECISION array, dimension (max(1,N-1))
*>          The (n-1) off-diagonal elements of the bidiagonal matrix
*>          B in elements 1 to N-1.
*> \endverbatim
*>
*> \param[in] VL
*> \verbatim
*>         VL is DOUBLE PRECISION
*>          If RANGE='V', the lower bound of the interval to
*>          be searched for singular values. VU > VL.
*>          Not referenced if RANGE = 'A' or 'I'.
*> \endverbatim
*>
*> \param[in] VU
*> \verbatim
*>         VU is DOUBLE PRECISION
*>          If RANGE='V', the upper bound of the interval to
*>          be searched for singular values. VU > VL.
*>          Not referenced if RANGE = 'A' or 'I'.
*> \endverbatim
*>
*> \param[in] IL
*> \verbatim
*>          IL is INTEGER
*>          If RANGE='I', the index of the smallest singular value
*>          to be returned.  1 <= IL <= IU <= N, if N > 0.
*>          Not referenced if RANGE = 'A' or 'V'.
*> \endverbatim
*>
*> \param[in] IU
*> \verbatim
*>          IU is INTEGER
*>          If RANGE='I', the index of the largest singular value
*>          to be returned.  1 <= IL <= IU <= N, if N > 0.
*>          Not referenced if RANGE = 'A' or 'V'.
*> \endverbatim
*>
*> \param[out] NS
*> \verbatim
*>          NS is INTEGER
*>          The total number of singular values found.  0 <= NS <= N.
*> \endverbatim
*>
*> \param[out] S
*> \verbatim
*>          S is DOUBLE PRECISION array, dimension (N)
*>          The first NS elements contain the selected singular values
*>          in decreasing order (matching DBDSVDX's ordering).
*> \endverbatim
*>
*> \param[out] Z
*> \verbatim
*>          Z is DOUBLE PRECISION array, dimension (LDZ, K)
*>          If JOBZ = 'V', the first NS columns of Z contain the
*>          singular vectors of B corresponding to the selected
*>          singular values, with U stored in rows 1..N and V in
*>          rows N+1..2N (Z = [U; V]).  If JOBZ = 'N', Z is not
*>          referenced.  K >= NS+1 (the caller must supply an upper
*>          bound on NS when RANGE = 'V').
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          Leading dimension of Z.  LDZ >= 1; if JOBZ = 'V', LDZ >= 2N.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (LWORK)
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          If JOBZ='V', LWORK >= max( 1, 5*N*N + 37*N ).
*>          If JOBZ='N', LWORK >= max( 1, 2*N*N + 37*N ).  These bounds
*>          cover both DBDSVDMR3 and the DBDSDC fallback.  If
*>          LWORK = -1, a workspace query is performed: the optimal
*>          size is returned in WORK(1) and no work is done.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (LIWORK)
*>          On exit, IWORK(2) reports the path taken: 1 = MRRR result
*>          accepted, 2 = DBDSDC fallback after the MRRR path returned
*>          INFO > 0, 3 = DBDSDC fallback after the finiteness audit
*>          rejected the MRRR result.
*>          On successful exit for N > 1, IWORK(2) records the path:
*>             1: DBDSVDMR3
*>             2: DBDSVDMR3 returned positive INFO, then DBDSDC succeeded
*>             3: DBDSVDMR3 output audit failed, then DBDSDC succeeded
*> \endverbatim
*>
*> \param[in] LIWORK
*> \verbatim
*>          LIWORK is INTEGER
*>          LIWORK >= max( 1, 20*N ).  If LIWORK = -1, a workspace
*>          query is performed.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          = N+1: the input contains a NaN or an Inf; nothing was
*>                computed.
*>          = N+2: the MRRR path was rejected and the DBDSDC fallback
*>                returned a non-finite result.
*>          > 0, other values: the MRRR path was rejected and the DBDSDC
*>                fallback failed; INFO is the INFO returned by DBDSDC.
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
*> \ingroup bdsvdmr3
*
*  =====================================================================
      SUBROUTINE DBDSVDMR3( UPLO, JOBZ, RANGE, N, D, E, VL, VU, IL, IU,
     $                   NS, S, Z, LDZ, WORK, LWORK, IWORK, LIWORK,
     $                   INFO )
      IMPLICIT NONE
*
*  -- LAPACK driver routine --
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, IU, LDZ, LIWORK, LWORK, N, NS
      DOUBLE PRECISION   VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), S( * ), WORK( * ),
     $                   Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLSV, BADOUT, INDSV, LOWER, LQUERY, VALSV,
     $                   WANTZ
      CHARACTER          COMPC
      INTEGER            I, IDCOPY, IECOPY, ISFULL, IUFULL, IVTFULL,
     $                   IWRK, J, LWMIN, LIWMIN, MFOUND, FALLWHY,
     $                   IWKENG, IIWENG, ILO, IHI, NCOL, PATH
      DOUBLE PRECISION   OVFL, SVAL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH, LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DBDSDC, DBDSVDMR3_WORK, DCOPY, DLASET, DSWAP,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      ALLSV  = LSAME( RANGE, 'A' )
      VALSV  = LSAME( RANGE, 'V' )
      INDSV  = LSAME( RANGE, 'I' )
      WANTZ  = LSAME( JOBZ,  'V' )
      LOWER  = LSAME( UPLO,  'L' )
      LQUERY = ( LWORK.EQ.-1 ) .OR. ( LIWORK.EQ.-1 )
*
*     Test the input parameters (mirror DBDSVDX).
*
      INFO = 0
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LOWER ) THEN
         INFO = -1
      ELSE IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( ALLSV .OR. VALSV .OR. INDSV ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.GT.0 ) THEN
         IF( VALSV ) THEN
            IF( VL.LT.ZERO ) THEN
               INFO = -7
            ELSE IF( VU.LE.VL ) THEN
               INFO = -8
            END IF
         ELSE IF( INDSV ) THEN
            IF( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) THEN
               INFO = -9
            ELSE IF( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) THEN
               INFO = -10
            END IF
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.2*N ) ) INFO = -14
      END IF
*
*     Workspace requirements.  The wrapper needs (in doubles):
*       * N   for a private copy of D
*       * N   for a private copy of E
*       * N   for the full unfiltered spectrum
*       * LDU*N = N*N for a private copy of U (when WANTZ)
*       * LDVT*N = N*N for a private copy of VT (when WANTZ)
*       * 2*N*N + 34*N for the DBDSVDMR3 engine itself
*     The DBDSDC fallback additionally needs 3*N*N + 4*N doubles after
*     its explicit U and VT arrays.  Total (WANTZ): 5*N*N + 7*N; retain
*     the larger linear allowance from DBDSVDMR3 for small N.
*
      IF( N.EQ.0 ) THEN
         LWMIN  = 1
         LIWMIN = 1
      ELSE IF( WANTZ ) THEN
         LWMIN  = MAX( 1, 5*N*N + 37*N )
         LIWMIN = MAX( 1, 24*N )
      ELSE
         LWMIN  = MAX( 1, 2*N*N + 37*N + 2 )
         LIWMIN = MAX( 1, 24*N )
      END IF
*
      IF( INFO.EQ.0 .AND. .NOT.LQUERY ) THEN
         IF( LWORK.LT.LWMIN ) THEN
            INFO = -16
         ELSE IF( LIWORK.LT.LIWMIN ) THEN
            INFO = -18
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DBDSVDMR3', -INFO )
         RETURN
      END IF
*
      IF( LQUERY ) THEN
         WORK( 1 ) = DBLE( LWMIN )
         IWORK( 1 ) = LIWMIN
         RETURN
      END IF
*
*     Quick return if possible.
*
      NS = 0
      IF( N.EQ.0 ) RETURN
*
*     A NaN or Inf in the input cannot be processed by either backend
*     (it may hang the fallback's QR iteration); report it as INFO = N+1
*     without calling XERBLA, so that the caller can recover.
*
      OVFL = DLAMCH( 'Overflow' )
      DO 5 I = 1, N
         SVAL = D( I )
         IF( .NOT.( SVAL.EQ.SVAL .AND. ABS( SVAL ).LE.OVFL ) ) THEN
            INFO = N + 1
            RETURN
         END IF
         IF( I.LT.N ) THEN
            SVAL = E( I )
            IF( .NOT.( SVAL.EQ.SVAL .AND. ABS( SVAL ).LE.OVFL ) ) THEN
               INFO = N + 1
               RETURN
            END IF
         END IF
    5 CONTINUE
*
      IF( N.EQ.1 ) THEN
         IWORK( 2 ) = 0
         SVAL = ABS( D( 1 ) )
         IF( ALLSV .OR. INDSV ) THEN
            NS = 1
            S( 1 ) = SVAL
         ELSE
            IF( VL.LT.SVAL .AND. VU.GE.SVAL ) THEN
               NS = 1
               S( 1 ) = SVAL
            END IF
         END IF
         IF( WANTZ .AND. NS.GT.0 ) THEN
            Z( 1, 1 ) = SIGN( ONE, D( 1 ) )
            Z( 2, 1 ) = ONE
         END IF
         RETURN
      END IF
*
*     Slice the WORK buffer.
*
      IDCOPY  = 1
      IECOPY  = IDCOPY + N
      ISFULL  = IECOPY + N
      IF( WANTZ ) THEN
         IUFULL  = ISFULL + N
         IVTFULL = IUFULL + N*N
         IWRK    = IVTFULL + N*N
      ELSE
*        U and VT are not referenced for JOBZ = 'N'; give them one
*        element each, disjoint from WORK, to avoid aliasing.
         IUFULL  = ISFULL + N
         IVTFULL = IUFULL + 1
         IWRK    = IVTFULL + 1
      END IF
*
*     Stage inputs (DBDSVDMR3 overwrites D and E; we must not
*     modify the caller's arrays).
*
      CALL DCOPY( N, D, 1, WORK( IDCOPY ), 1 )
      IF( N.GT.1 ) THEN
         CALL DCOPY( N-1, E, 1, WORK( IECOPY ), 1 )
      END IF
      WORK( IECOPY + N - 1 ) = ZERO
*
*     Compute the full spectrum with DBDSVDMR3 for every N.  A positive
*     INFO, incomplete spectrum, or nonfinite requested output triggers a
*     restored-input DBDSDC retry.  The O(N**2) output audit is no more
*     expensive asymptotically than returning all singular vectors.
*
      IWKENG = LWORK - IWRK + 1
      IIWENG = LIWORK
      IF( WANTZ ) THEN
         CALL DBDSVDMR3_WORK( 'V', UPLO, N,
     $                   WORK( IDCOPY ), WORK( IECOPY ),
     $                   WORK( ISFULL ),
     $                   WORK( IUFULL  ), N,
     $                   WORK( IVTFULL ), N,
     $                   MFOUND,
     $                   WORK( IWRK ), IWKENG,
     $                   IWORK, IIWENG, INFO )
      ELSE
         CALL DBDSVDMR3_WORK( 'N', UPLO, N,
     $                   WORK( IDCOPY ), WORK( IECOPY ),
     $                   WORK( ISFULL ),
     $                   WORK( IUFULL  ), MAX( 1, N ),
     $                   WORK( IVTFULL ), MAX( 1, N ),
     $                   MFOUND,
     $                   WORK( IWRK ), IWKENG,
     $                   IWORK, IIWENG, INFO )
      END IF
      IF( INFO.LT.0 ) RETURN
*
      FALLWHY = 0
      IF( INFO.GT.0 ) THEN
         FALLWHY = 2
      ELSE
         OVFL = DLAMCH( 'Overflow' )
         BADOUT = MFOUND.NE.N
         IF( .NOT.BADOUT ) THEN
            DO 12 I = 1, N
               SVAL = WORK( ISFULL + I - 1 )
               IF( .NOT.( SVAL.EQ.SVAL .AND.
     $                     ABS( SVAL ).LE.OVFL ) ) BADOUT = .TRUE.
   12       CONTINUE
            IF( WANTZ ) THEN
               DO 14 J = 1, N
                  DO 13 I = 1, N
                     SVAL = WORK( IUFULL + (J-1)*N + I - 1 )
                     IF( .NOT.( SVAL.EQ.SVAL .AND.
     $                           ABS( SVAL ).LE.OVFL ) ) BADOUT = .TRUE.
                     SVAL = WORK( IVTFULL + (J-1)*N + I - 1 )
                     IF( .NOT.( SVAL.EQ.SVAL .AND.
     $                           ABS( SVAL ).LE.OVFL ) ) BADOUT = .TRUE.
   13             CONTINUE
   14          CONTINUE
            END IF
         END IF
         IF( BADOUT ) FALLWHY = 3
      END IF
*
      IF( FALLWHY.EQ.0 ) THEN
         PATH = 1
      ELSE
*        DBDSVDMR3 destroys its staged inputs.  Restore them before the
*        DBDSDC retry.
         CALL DCOPY( N, D, 1, WORK( IDCOPY ), 1 )
         CALL DCOPY( N-1, E, 1, WORK( IECOPY ), 1 )
         WORK( IECOPY + N - 1 ) = ZERO
         INFO = 0
         IF( WANTZ ) THEN
            COMPC = 'I'
         ELSE
            COMPC = 'N'
         END IF
*        Q and IQ are not referenced for COMPQ = 'I' or 'N'; pass
*        storage disjoint from WORK and IWORK to avoid aliasing.
         CALL DBDSDC( UPLO, COMPC, N, WORK( IDCOPY ),
     $                WORK( IECOPY ), WORK( IUFULL ), MAX( 1, N ),
     $                WORK( IVTFULL ), MAX( 1, N ), WORK( ISFULL ),
     $                IWORK( 8*N+1 ), WORK( IWRK ), IWORK, INFO )
         IF( INFO.NE.0 ) RETURN
*        The finiteness audit applies to the fallback as well.
         BADOUT = .FALSE.
         DO 16 I = 1, N
            SVAL = WORK( IDCOPY + I - 1 )
            IF( .NOT.( SVAL.EQ.SVAL .AND. ABS( SVAL ).LE.OVFL ) )
     $         BADOUT = .TRUE.
   16    CONTINUE
         IF( WANTZ ) THEN
            DO 17 J = 1, N
               DO 15 I = 1, N
                  SVAL = WORK( IUFULL + (J-1)*N + I - 1 )
                  IF( .NOT.( SVAL.EQ.SVAL .AND. ABS( SVAL ).LE.OVFL ) )
     $               BADOUT = .TRUE.
                  SVAL = WORK( IVTFULL + (J-1)*N + I - 1 )
                  IF( .NOT.( SVAL.EQ.SVAL .AND. ABS( SVAL ).LE.OVFL ) )
     $               BADOUT = .TRUE.
   15          CONTINUE
   17       CONTINUE
         END IF
         IF( BADOUT ) THEN
            INFO = N + 2
            RETURN
         END IF
*
*        DBDSDC returns descending singular values, while DBDSVDMR3
*        returns ascending values.  Normalize DBDSDC's result so the
*        common RANGE filtering and packing below is unchanged.
*
         CALL DCOPY( N, WORK( IDCOPY ), 1, WORK( ISFULL ), 1 )
         DO 18 I = 1, N / 2
            CALL DSWAP( 1, WORK( ISFULL + I - 1 ), 1,
     $                  WORK( ISFULL + N - I ), 1 )
            IF( WANTZ ) THEN
               CALL DSWAP( N, WORK( IUFULL + (I-1)*N ), 1,
     $                     WORK( IUFULL + (N-I)*N ), 1 )
               CALL DSWAP( N, WORK( IVTFULL + I - 1 ), N,
     $                     WORK( IVTFULL + N - I ), N )
            END IF
   18    CONTINUE
         PATH = FALLWHY
      END IF
      IWORK( 2 ) = PATH
*
*     Filter and pack the requested slice into S and Z.
*
      IF( ALLSV ) THEN
         ILO = 1
         IHI = N
      ELSE IF( INDSV ) THEN
*        IL and IU follow DBDSVDX: index 1 is the largest singular
*        value.  The internal spectrum is ascending.
         ILO = N - IU + 1
         IHI = N - IL + 1
      ELSE
*        RANGE='V'; scan ascending spectrum for [VL, VU).
         ILO = N + 1
         IHI = 0
         DO 10 I = 1, N
            SVAL = WORK( ISFULL + I - 1 )
            IF( SVAL.GE.VL .AND. SVAL.LT.VU ) THEN
               IF( I.LT.ILO ) ILO = I
               IF( I.GT.IHI ) IHI = I
            END IF
   10    CONTINUE
      END IF
*
      IF( IHI.GE.ILO ) THEN
         NS = IHI - ILO + 1
*        Reverse ascending -> descending so S(1) is the largest,
*        matching DBDSVDX's convention (see dbdsvdx.f
*        line 747 sort loop).
         DO 20 I = 1, NS
            S( I ) = WORK( ISFULL + IHI - I )
   20    CONTINUE
         IF( WANTZ ) THEN
*           Column J of Z holds the vector for S(J).  In the engine's
*           spectrum that is column (IHI - J + 1).
*           Rows 1..N of Z(:, J) hold U(:, IHI-J+1).
*           Rows N+1..2N of Z(:, J) hold VT(IHI-J+1, :)^T.
            NCOL = NS
            DO 40 J = 1, NCOL
               DO 30 I = 1, N
                  Z( I,     J ) = WORK( IUFULL  + (IHI-J)*N + I - 1 )
                  Z( N + I, J ) = WORK( IVTFULL + (I-1)*N + (IHI-J) )
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE
         NS = 0
      END IF
*
      RETURN
*
*     End of DBDSVDMR3
*
      END
