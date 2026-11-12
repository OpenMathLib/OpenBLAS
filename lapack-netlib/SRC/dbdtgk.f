*> \brief \b DBDTGK forms the Golub-Kahan tridiagonal matrix of a bidiagonal matrix.
*>
*> \ingroup bdsvdmr3
*>
*> \par Purpose:
*> =============
*>
*> \verbatim
*> forms the Golub-Kahan tridiagonal matrix of a bidiagonal matrix.
*> All routines in this file are auxiliary to DBDSVDMR3 and derive
*> from the LAPACK 3.0 (1999) MRRR kernels DLARRV, DLAR1V, DLARRB and
*> DLARRF rather than from the current DSTEMR kernels.
*> \endverbatim
*
      SUBROUTINE DBDTGK( UPLO, N, D, E, DTGK, ETGK, GERSCH )
*
*  -- New auxiliary routine for bidiagonal SVD via TGK-rooted MR^3 --
*     Builds the Tridiagonal Golub-Kahan (TGK) matrix of a bidiagonal B.
*     Standalone; does not modify any other routine.
*
*  Purpose
*  =======
*
*  Given a bidiagonal matrix B of order N (diagonal D, off-diagonal E),
*  DBDTGK forms the symmetric tridiagonal Golub-Kahan matrix T_GK of
*  order 2N.  T_GK has a ZERO diagonal and off-diagonals that interleave
*  the entries of B:
*       offdiag(T_GK) = [ D(1), E(1), D(2), E(2), ..., E(N-1), D(N) ].
*  Its eigenvalues are +/- the singular values of B; the eigenvector for
*  +sigma_i interleaves the left and right singular vectors:
*       z(2k-1) ~ u_i(k),  z(2k) ~ v_i(k)   (each scaled by 1/sqrt(2)).
*  The Gerschgorin intervals of T_GK are also returned (used by
*  DLAR1V_TGK to restrict the twist-index search); since the diagonal is
*  zero, interval j is [ -rad_j, +rad_j ] with rad_j the row sum of the
*  off-diagonal magnitudes.
*
*  Arguments
*  =========
*
*  UPLO   (input) CHARACTER*1
*         = 'U':  B is upper bidiagonal (E is the superdiagonal);
*         = 'L':  B is lower bidiagonal (E is the subdiagonal).
*         The TGK off-diagonal ordering is identical for both; UPLO is
*         retained for interface clarity and future use (it only affects
*         how the recovered z-entries map to u vs v, handled by caller).
*
*  N      (input) INTEGER   The order of B (B is N-by-N).
*
*  D      (input) DOUBLE PRECISION array, dimension (N)
*         The diagonal entries of B.
*
*  E      (input) DOUBLE PRECISION array, dimension (N-1)
*         The off-diagonal entries of B.
*
*  DTGK   (output) DOUBLE PRECISION array, dimension (2*N)
*         The diagonal of T_GK (all zeros).
*
*  ETGK   (output) DOUBLE PRECISION array, dimension (2*N-1)
*         The off-diagonal beta_i of T_GK (interleaved D and E).
*
*  GERSCH (output) DOUBLE PRECISION array, dimension (2*(2*N))
*         The 2N Gerschgorin intervals; interval j is
*         (GERSCH(2*j-1), GERSCH(2*j)) = (-rad_j, +rad_j).
*
*  =====================================================================
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), DTGK( * ), ETGK( * ),
     $                   GERSCH( * )
*     ..
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, M
      DOUBLE PRECISION   RAD
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      M = 2*N
*
*     Diagonal of T_GK is identically zero.
*
      DO 10 I = 1, M
         DTGK( I ) = ZERO
   10 CONTINUE
*
*     Interleave D and E into the off-diagonal of T_GK:
*        ETGK(2k-1) = D(k),  ETGK(2k) = E(k).
*
      DO 20 I = 1, N
         ETGK( 2*I-1 ) = D( I )
   20 CONTINUE
      DO 30 I = 1, N - 1
         ETGK( 2*I ) = E( I )
   30 CONTINUE
*
*     Gerschgorin intervals: zero diagonal => centered at 0, radius is the
*     sum of the magnitudes of the (at most two) adjacent off-diagonals.
*
      DO 40 I = 1, M
         RAD = ZERO
         IF( I.GT.1 ) RAD = RAD + ABS( ETGK( I-1 ) )
         IF( I.LT.M ) RAD = RAD + ABS( ETGK( I ) )
         GERSCH( 2*I-1 ) = -RAD
         GERSCH( 2*I ) = RAD
   40 CONTINUE
*
      RETURN
*
*     End of DBDTGK
*
      END
