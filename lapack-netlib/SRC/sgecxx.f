*> \brief \b SGECXX computes a CX factorization of a real M-by-N matrix A using a truncated (rank k) Householder QR factorization with column pivoting.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SGECXX + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgecxx.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgecxx.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgecxx.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGECXX( FACT, USESD, M, N,
*      $                   DESEL_ROWS, SEL_DESEL_COLS,
*      $                   KMAXFREE, ABSTOL, RELTOL, A, LDA,
*      $                   K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
*      $                   IPIV, JPIV, TAU, C, LDC, QRC, LDQRC,
*      $                   X, LDX, WORK, LWORK, IWORK, LIWORK, INFO )
*       IMPLICIT NONE
*
*      .. Scalar Arguments ..
*       CHARACTER           FACT, USESD
*       INTEGER             INFO, K, KMAXFREE, LDA, LDC, LDQRC,
*      $                    LDX, LIWORK, LWORK, M, N
*       REAL                ABSTOL, FNRMK, MAXC2NRMK,
*      $                    RELMAXC2NRMK, RELTOL
*      ..
*      .. Array Arguments ..
*       INTEGER             DESEL_ROWS( * ), IPIV( * ), IWORK( * ),
*      $                    JPIV( * ), SEL_DESEL_COLS( * )
*       REAL                A( LDA, * ), C( LDC, * ), QRC( LDQRC, * ),
*      $                    TAU( * ), WORK( * ), X( LDX, *)
*      ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGECXX computes a CX factorization of a real M-by-N matrix A using
*> a truncated rank-K Householder QR factorization with a column
*> pivoting algorithm, which is implemented in the SGEQP3RK routine.
*>
*>   A * P = C*X + A_resid, where
*>
*>   C is an M-by-K matrix consisting of K columns selected
*>     from the original matrix A,
*>
*>   X is a K-by-N matrix that minimizes the Frobenius norm of the
*>     residual matrix A_resid, X = pseudoinv(C) * A,
*>
*>   P is an N-by-N permutation matrix chosen so that the first
*>     K columns of A*P equal C,
*>
*>   A_resid is an M-by-N residual matrix.
*>
*> The column selection for the matrix C has two stages.
*>
*> Column preselection stage 1 (optional).
*> =======================================
*>
*> The user can select N_sel columns and deselect N_desel columns
*> of the matrix A that MUST be included and excluded respectively
*> from the matrix C a priori, before running the column selection
*> algorithm. This is controlled by flags in the array
*> SEL_DESEL_COLS. The deselected columns are permuted to the right
*> side of the matrix A and selected columns are permuted to the left
*> side of the matrix A. The details of the column permutation
*> (i.e. the column permutation matrix P) are stored in the
*> array JPIV. This feature can be used when the goal is to approximate
*> the deselected columns by linear combinations of K selected columns,
*> where the K columns MUST include the N_sel preselected columns.
*>
*> Column selection stage 2.
*> =========================
*>
*> The routine runs a column selection algorithm that can
*> be controlled by three stopping criteria described below.
*> For column selection, the routine uses a truncated (rank-K)
*> Householder QR factorization with column pivoting algorithm using
*> the routine SGEQP3RK.
*>
*> Optionally, before running the column selection
*> algorithm, the user can deselect M_desel rows of the matrix A that
*> should NOT be considered by the column selection algorithm (i.e.
*> during the factorization). This is controlled by flags in
*> the array DESEL_ROWS. The deselected rows are permuted to the
*> bottom of the matrix A. The details of the row permutation (i.e. the
*> row permutation matrix) are stored in the array IPIV. This feature
*> can be used when the goal is to use the deselected rows as test data,
*> and the selected rows as training data.
*>
*> This means that the column selection factorization algorithm is
*> effectively running on the submatrix A_sub = A(1:M_sub,1:N_sub) of
*> the matrix A after the permutations described above. Here M_sub is
*> the number of rows of the matrix A minus the number of deselected
*> rows M_desel, i.e. M_sub = M - M_desel, and N_sub is the number
*> of columns of the matrix A minus the number of deselected columns
*> N_desel, i.e. N_sub = N - N_desel.
*>
*> The reported column selection error metrics MAXC2NRMK, RELMAXC2NRMK
*> and FNRMK described below are computed using only A_sub.
*>
*> Column selection criteria.
*> ==========================
*>
*> The column selection criteria (i.e. when to stop the factorization)
*> can be any of the following:
*>
*>   1) KMAXFREE: This input parameter specifies the maximum number of
*>      columns to factorize in addition to the N_sel preselected
*>      columns. The factorization rank is limited to N_sel + KMAXFREE.
*>      If N_sel + KMAXFREE >= min(M_sub, N_sub), this criterion
*>      is not used.
*>
*>   2) ABSTOL: This input parameter specifies the absolute tolerance
*>      for the maximum column 2-norm of the submatrix residual
*>      A_sub_resid(K) = A_sub(K)(K+1:M_sub, K+1:N_sub), where
*>      A_sub(K) denotes the contents of the array
*>      A_sub = A(1:M_sub, 1:N_sub) after K columns were factorized.
*>      This means that the factorization stops if this norm is less
*>      than or equal to ABSTOL. If ABSTOL < 0.0, this criterion is
*>      not used.
*>
*>   3) RELTOL: This input parameter specifies the tolerance for
*>      the maximum column 2-norm of the submatrix residual
*>      A_sub_resid(K) = A_sub(K)(K+1:M_sub, K+1:N_sub) divided
*>      by the maximum column 2-norm of the submatrix
*>      A_sub = A(1:M_sub, 1:N_sub), where A_sub(K) denotes the contents
*>      of the array A_sub after K columns were factorized.
*>      This means that the factorization stops when the ratio of the
*>      maximum column 2-norm of A_sub_resid(K) to the maximum column
*>      2-norm of A_sub is less than or equal to RELTOL.
*>      If RELTOL < 0.0, this criterion is not used.
*>
*> The algorithm stops when any of these conditions is first
*> satisfied, otherwise the entire submatrix A_sub is factorized.
*>
*> To perform a full-rank factorization of the matrix A_sub, use
*> selection criteria that satisfy N_sel + KMAXFREE >= min(M_sub,N_sub)
*> and ABSTOL < 0.0 and RELTOL < 0.0.
*>
*> If the user wishes to verify that the columns of the matrix C are
*> sufficiently linearly independent for their intended use, the user
*> can compute the condition number of its R factor by calling DTRCON
*> on the upper-triangular part of QRC(1:K,1:K) in the output
*> array QRC.
*>
*> How N_sel affects the column selection algorithm.
*> =================================================
*>
*> As mentioned above, the N_sel preselected columns are permuted to the
*> left side of the matrix A, and will be included in the column
*> selection. Then the routine factorizes that block A(1:M_sub,1:N_sel),
*> and if any of the three stopping criteria is met immediately after
*> factoring the first N_sel columns the routine exits
*> (i.e. if the user does not want to select KMAXFREE > 0 extra columns,
*> or if the absolute or relative tolerance of the maximum column 2-norm
*> of the residual is satisfied). In this case, the number
*> of selected columns would be K = N_sel. Otherwise, the factorization
*> routine finds a new column to select with the maximum column 2-norm
*> in the residual A(N_sel+1:M_sub,N_sel+1:N_sub), and swaps that
*> column with the first column of A(1:M,N_sel+1:N_sub). Then the
*> routine checks if the stopping criteria are met in the next residual
*> A(N_sel+2:M_sub,N_sel+2:N_sub), and so on.
*>
*> Computation of the matrix factors.
*> ==================================
*>
*> When the columns are selected for the factor C, and:
*>  (a) If the flag FACT = 'P', the routine returns only the indices of
*>      the selected columns from the original matrix A, which are
*>      stored in the first K elements of the JPIV array.
*>  (b) If the flag FACT = 'C', then in addition to (a), the routine
*>      explicitly returns the matrix C in the array C.
*>  (c) If the flag FACT = 'X', then in addition to (a) and (b),
*>      the routine explicitly computes and returns the factor
*>      X = pseudoinv(C) * A in the array X, and it also returns
*>      the factor R alongside the Householder vectors
*>      of the QR factorization of the matrix C in the array QRC.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] FACT
*> \verbatim
*>          FACT is CHARACTER*1
*>          The flag specifies how the factors of a CX factorization
*>          are returned.
*>
*>          = 'P': the routine returns:
*>                 (1) only the column permutation matrix P in
*>                     the array JPIV.
*>                     (The first K elements of the array JPIV
*>                     contain indices of the columns that were
*>                     selected from the matrix A to form the
*>                     factor C.)
*>                 (fastest option, smallest memory space)
*>
*>          = 'C': the routine returns:
*>                 (1) the column permutation matrix P
*>                     in the array JPIV. (The first K elements are
*>                     indices of the selected columns from
*>                     the matrix A.)
*>                 (2) the M-by-K factor C explicitly in the array C.
*>                 (slower option, more memory space)
*>
*>          = 'X': the routine returns:
*>                 (1) the column permutation matrix P in
*>                     the array JPIV. (The first K elements are
*>                     indices of the selected columns from
*>                     the matrix A.)
*>                 (2) the M-by-K factor C explicitly in the array C.
*>                 (3) the K-by-N factor X explicitly in the array X.
*>                 (4) the K-by-K upper triangular factor R and
*>                     the Householder vectors of the QR factorization
*>                     of the factor C in the array QRC.
*>                     ( The factor R may be useful for checking
*>                       the factor C for singularity, in which case
*>                       R will have a zero on the diagonal, and
*>                       the factor X cannot be computed. )
*>                 (slowest option, largest memory space)
*> \endverbatim
*>
*> \param[in] USESD
*> \verbatim
*>          USESD is CHARACTER*1
*>          The flag specifies whether the row deselection and column
*>          preselection-deselection functionality is turned ON or OFF.
*>
*>          = 'N': Both row deselection and column
*>                 preselection-deselection are OFF.
*>                 Both arrays DESEL_ROWS and SEL_DESEL_COLS
*>                 are not used.
*>
*>          = 'R': Only row deselection is ON.
*>                 Column preselection-deselection is OFF.
*>                 The array SEL_DESEL_COLS is not used.
*>
*>          = 'C': Only column preselection-deselection is ON.
*>                 Row deselection is OFF.
*>                 The array DESEL_ROWS is not used.
*>
*>          = 'A': Means "All". Both row deselection and column
*>                 preselection-deselection are ON.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] DESEL_ROWS
*> \verbatim
*>          DESEL_ROWS is INTEGER array, dimension (M)
*>          DESEL_ROWS is only accessed if USESD = 'R' or 'A'.
*>          This is a row deselection mask array that separates
*>          the rows of matrix A into 2 sets.
*>
*>          On entry:
*>          a) If DESEL_ROWS(i) = -1, the i-th row of the matrix A is
*>             deselected by the user, i.e. chosen to be excluded from
*>             the column selection algorithm (in both preselection and
*>             selection stages) and will be permuted to the bottom
*>             of the matrix A.
*>             The number of deselected rows is denoted by M_desel.
*>
*>          b) If DESEL_ROWS(i) is not equal -1,
*>             the i-th row of A will be used in the column selection
*>             algorithm (in both preselection and selection stages).
*>             This defines a set of M_sub = M - M_desel rows that
*>             the algorithm will use to select columns.
*>             After the permutation, this set will be at the top
*>             of the matrix A.
*>
*>           On exit:
*>             DESEL_ROWS will be permuted according to IPIV(i),
*>             so that, if IPIV(i) = k, then the entry i of DESEL_ROWS
*>             on exit was the entry k of DESEL_ROWS on entry.
*>
*> \endverbatim
*>
*> \param[in,out] SEL_DESEL_COLS
*> \verbatim
*>          SEL_DESEL_COLS is INTEGER array, dimension (N)
*>          SEL_DESEL_COLS is only accessed if USESD = 'C' or 'A'.
*>          This is a column preselection-deselection mask array that
*>          separates the columns of matrix A into 3 sets.
*>
*>          On entry:
*>          a) If SEL_DESEL_COLS(j) = +1, the j-th column of the matrix
*>             A is preselected by the user to be included
*>             in the factor C and will be permuted to the left side
*>             of the array A. The number of selected columns is
*>             denoted by N_sel.
*>
*>          b) If SEL_DESEL_COLS(j) = -1, the j-th column of the matrix
*>             A is deselected by the user, i.e. chosen to be excluded
*>             from the factor C and will be permuted to the right side
*>             of the array A. The number of deselected columns is
*>             denoted by N_desel.
*>
*>          c) If SEL_DESEL_COLS(j) is not equal to 1 and not equal
*>             to -1, the j-th column of A is a free column and will be
*>             used by the column selection algorithm to determine if
*>             this column will be selected. This defines a set of
*>             columns of size N_free = N - N_sel - N_desel.
*>
*>           On exit:
*>             SEL_DESEL_COLS will be permuted according to JPIV(j),
*>             so that, if JPIV(j) = k, then the entry j
*>             of SEL_DESEL_COLS on exit was the entry k
*>             of SEL_DESEL_COLS on entry.
*>
*>          NOTE: An error returned as INFO = -6 means that the number
*>          of preselected N_sel columns is larger than M_sub.
*>          Therefore, the QR factorization of all N_sel preselected
*>          columns cannot be completed.
*> \endverbatim
*>
*> \param[in] KMAXFREE
*> \verbatim
*>          KMAXFREE is INTEGER, KMAXFREE >= 0.
*>
*>          The first column selection stopping criterion from
*>          the N_free columns (N_sel+1:N_sub) of the submatrix
*>          A_sub = A(1:M_sub, 1:N_sub) in the column selection stage 2.
*>
*>          KMAXFREE is the maximum number of columns of the matrix
*>          A_free = A(N_sel+1:M_sub, N_sel+1:N_sub) to select
*>          during the column selection stage 2.
*>
*>          KMAXFREE does not include the preselected N_sel columns.
*>          N_sel + KMAXFREE is the maximum factorization rank of
*>          the matrix A_sub.
*>
*>          a) If N_sel + KMAXFREE >= min(M_sub, N_sub), then this
*>             stopping criterion is not used, i.e. columns are
*>             selected in the factorization stage 2 depending
*>             on ABSTOL and RELTOL.
*>
*>          b) If KMAXFREE = 0, then this stopping criterion is
*>             satisfied on input and the routine exits without
*>             performing column selection stage 2
*>             on the submatrix A_sub. This means that the matrix
*>             A_free = A(N_sel+1:M_sub, N_sel+1:N_sub) is not modified
*>             in the column selection stage 2
*>             and A_free is itself the residual for the factorization.
*> \endverbatim
*>
*> \param[in] ABSTOL
*> \verbatim
*>          ABSTOL is REAL, cannot be NaN.
*>
*>          The second column selection stopping criterion from
*>          the N_free columns (N_sel+1:N_sub) of the submatrix
*>          A_sub = A(1:M_sub, 1:N_sub) in the column selection stage 2.
*>
*>          ABSTOL is the absolute tolerance (stopping threshold)
*>          for maxcol2norm(A_sub_resid(K)), where K >= N_sel.
*>
*>          maxcol2norm(A_sub_resid(K)) is the maximum column 2-norm
*>          of the residual matrix
*>          A_sub_resid(K) = A_sub(K)(K+1:M_sub, K+1:N_sub)
*>          when K columns have been factorized.
*>          The column selection algorithm converges
*>          (stops the factorization) when
*>          maxcol2norm(A_sub_resid(K)) <= ABSTOL, where K >= N_sel.
*>
*>          In the following,
*>                SAFMIN = SLAMCH('S'),
*>                A_free = A(N_sel+1:M_sub, N_sel+1:N_sub),
*>                maxcol2norm(A_free) is the maximum column 2-norm
*>                of the matrix A_free.
*>
*>          a) If ABSTOL is NaN, then no computation is performed
*>                and an error message ( INFO = -8 ) is issued
*>                by XERBLA.
*>
*>          b) If ABSTOL < 0.0, then this stopping criterion is not
*>                used, and the column selection algorithm stops
*>                the factorization of A_free depending
*>                on KMAXFREE and RELTOL.
*>                This includes the case where ABSTOL = -Inf.
*>
*>          c) If 0.0 <= ABSTOL < 2*SAFMIN, then ABSTOL = 2*SAFMIN
*>                is used. This includes the case where ABSTOL = -0.0.
*>
*>          d) If 2*SAFMIN <= ABSTOL then the input value
*>                of ABSTOL is used.
*>
*>          If ABSTOL chosen above is >= maxcol2norm(A_free), then
*>          this stopping criterion is satisfied on input, and
*>          the routine only preselects K = N_sel columns. The leftmost
*>          preselected N_sel columns in the submatrix
*>          A_sub = A(1:M_sub, 1:N_sub) are factorized. The routine
*>          then computes maxcol2norm(A_free) and returns it
*>          in MAXC2NORMK, computes and returns RELMAXC2NORMK of A_free,
*>          and exits immediately.
*>          This means that the factorization residual
*>          A_sub_resid(N_sel) = A_free = A(N_sel+1:M_sub,N_sel+1:N_sub)
*>          is not modified in the column selection stage 2.
*>          This includes the case where ABSTOL = +Inf.
*> \endverbatim
*>
*> \param[in] RELTOL
*> \verbatim
*>          RELTOL is REAL, cannot be NaN.
*>
*>          The third column selection stopping criterion from
*>          the N_free columns (N_sel+1:N_sub) of the submatrix
*>          A_sub = A(1:M_sub, 1:N_sub) in the column selection stage 2.
*>
*>          RELTOL is the tolerance (stopping threshold) for the ratio
*>          relmaxcol2norm(A_sub_resid(K)) =
*>                 = maxcol2norm(A_sub_resid(K))/maxcol2norm(A_sub),
*>          where K >= N_sel.
*>
*>          maxcol2norm(A_sub_resid(K)) is the maximum column 2-norm
*>          of the residual matrix
*>          A_sub_resid(K) = A_sub(K)(K+1:M_sub, K+1:N_sub)
*>          when K columns have been factorized.
*>          maxcol2norm(A_sub) is the maximum column 2-norm
*>          of the original submatrix A_sub = A(1:M_sub, 1:N_sub).
*>          The column selection algorithm converges
*>          (stops the factorization) when the ratio
*>          relmaxcol2norm(A_sub_resid(K)) <= RELTOL, where K >= N_sel.
*>
*>          In the following,
*>                EPS = SLAMCH('E'),
*>                A_free = A(N_sel+1:M_sub, N_sel+1:N_sub).
*>
*>          a) If RELTOL is NaN, then no computation is performed
*>                and an error message ( INFO = -9 ) is issued
*>                by XERBLA.
*>
*>          b) If RELTOL < 0.0, then this stopping criterion is not
*>                used and the column selection algorithm stops
*>                the factorization of A_free depending
*>                on KMAXFREE and ABSTOL.
*>                This includes the case RELTOL = -Inf.
*>
*>          c) If 0.0 <= RELTOL < EPS, then RELTOL = EPS is used.
*>                This includes the case RELTOL = -0.0.
*>
*>          d) If EPS <= RELTOL then the input value of RELTOL
*>                is used.
*>
*>          If RELTOL chosen above is >= 1.0, then this stopping
*>          criterion is satisfied on input, and the routine
*>          only preselects K = N_sel columns. The leftmost
*>          preselected N_sel columns in the submatrix
*>          A_sub = A(1:M_sub, 1:N_sub) are factorized.
*>          The routine then computes maxcol2norm(A_free) and returns
*>          it in MAXC2NORMK, returns RELMAXC2NORMK as 1.0, and exits
*>          immediately.
*>          This means that the factorization residual
*>          A_sub_resid(N_sel) = A_free = A(N_sel+1:M_sub,N_sel+1:N_sub)
*>          is not modified.
*>          This includes the case RELTOL = +Inf.
*>
*>          NOTE: We recommend RELTOL to satisfy
*>                min(max(M_sub,N_sub)*EPS, sqrt(EPS)) <= RELTOL
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>
*>          On entry:
*>            the M-by-N matrix A.
*>
*>          On exit:
*>
*>            NOTE:
*>            The output parameter K, the number of selected
*>            columns, is described later.
*>            A_sub = A(1:M_sub, 1:N_sub).
*>
*>            1) If K = 0, A(1:M,1:N) contains the original matrix A.
*>
*>            2) If K > 0, A(1:M,1:N) contains the following parts:
*>
*>            (a) If M_sub < M (which is the same as M_desel > 0),
*>                the subarray A(M_sub+1:M,1:N) contains the deselected
*>                rows.
*>
*>            (b) If N_sub < N ( which is the same as N_desel > 0 ),
*>                the subarray A(1:M,N_sub+1:N) contains the
*>                deselected columns.
*>
*>            (c) If N_sel > 0,
*>                the union of the subarray A(1:M_sub, 1:N_sel)
*>                and the subarray A(1:N_sel, 1:N_sub) contains parts
*>                of the factors obtained by computing Householder QR
*>                factorization WITHOUT column pivoting of N_sel
*>                preselected columns using the routine SGEQRF.
*>
*>            (d) The subarray A(N_sel+1:M_sub, N_sel+1:N_sub)
*>                contains parts of the factors obtained by computing
*>                a truncated (rank K) Householder QR factorization with
*>                column pivoting using the routine SGEQP3RK on
*>                the matrix A_free = A(N_sel+1:M_sub, N_sel+1:N_sub),
*>                which is the result of applying selection and
*>                deselection of columns, applying deselection of rows
*>                to the original matrix A, and applying orthogonal
*>                transformation from the factorization of the first
*>                N_sel columns as described in part (c).
*>
*>             1. The elements below the diagonal of the subarray
*>                A_sub(1:M_sub,1:K) together with TAU(1:K)
*>                represent the orthogonal matrix Q(K) as a
*>                product of K Householder elementary reflectors.
*>
*>             2. The elements on and above the diagonal of
*>                the subarray A_sub(1:K,1:N_sub) contain the
*>                K-by-N_sub upper-trapezoidal matrix
*>                R_sub_approx(K) = ( R_sub11(K), R_sub12(K) ).
*>                NOTE: If K = min(M_sub,N_sub), i.e. full rank
*>                      factorization, then R_sub_approx(K) is the
*>                      full factor R which is upper-trapezoidal.
*>                      If, in addition, M_sub >= N_sub, then R is
*>                      upper-triangular.
*>
*>             3. The subarray A_sub(K+1:M_sub,K+1:N_sub) contains
*>                the (M_sub-K)-by-(N_sub-K) rectangular matrix
*>                A_sub_resid(K) = A_sub(K)(K+1:M_sub, K+1:N_sub).
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] K
*> \verbatim
*>          K is INTEGER
*>          The number of columns that were selected
*>          (K is the factorization rank).
*>          0 <= K <= min( M_sub, N_sel+KMAXFREE, N_sub ).
*>
*>          NOTE: If K = 0, a) the arrays A is not, modified.
*>                          b) the array TAU(1,min(M_sub,N_sub))
*>                             is set to ZERO.
*> \endverbatim
*>
*> \param[out] MAXC2NRMK
*> \verbatim
*>          MAXC2NRMK is REAL
*>          The maximum column 2-norm of the residual matrix
*>          A_sub_resid(K) = A_sub(K)(K+1:M_sub, K+1:N_sub),
*>          when factorization stopped at rank K. MAXC2NRMK >= 0.
*>
*>          a) If K = 0, i.e. the factorization was not performed, so
*>             the matrix A_sub = A(1:M_sub, 1:N_sub) was not modified
*>             and is itself a residual matrix, then MAXC2NRMK equals
*>             the maximum column 2-norm of the original matrix A_sub.
*>
*>          b) If 0 < K < min(M_sub, N_sub), then MAXC2NRMK is returned.
*>
*>          c) If K = min(M_sub, N_sub), i.e. the whole matrix A_sub was
*>             factorized and there is no residual matrix,
*>             then MAXC2NRMK = 0.0.
*>
*>          NOTE: MAXC2NRMK at the factorization step K is equal
*>                to the diagonal element R_sub(K+1,K+1) of the factor
*>                R_sub in the next factorization step K+1.
*> \endverbatim
*>
*> \param[out] RELMAXC2NRMK
*> \verbatim
*>          RELMAXC2NRMK is REAL
*>          The ratio MAXC2NRMK / MAXC2NRM
*>          of the maximum column 2-norm MAXC2NRMK of the residual
*>          matrix A_sub_resid(K) = A_sub(K+1:M_sub, K+1:N_sub) (when
*>          factorization stopped at rank K) and maximum column 2-norm
*>          MAXC2NRM of the matrix A_sub = A(1:M_sub, 1:N_sub).
*>          RELMAXC2NRMK >= 0.
*>
*>          a) If K = 0, i.e. the factorization was not performed,
*>             the matrix  A_sub was not modified
*>             and is itself a residual matrix,
*>             then RELMAXC2NRMK = 1.0.
*>
*>          b) If 0 < K < min(M_sub,N_sub), then
*>                RELMAXC2NRMK = MAXC2NRMK / MAXC2NRM is returned.
*>
*>          c) If K = min(M_sub,N_sub), i.e. the whole matrix A_sub was
*>             factorized and there is no residual matrix
*>             A_sub_resid(K), then RELMAXC2NRMK = 0.0.
*>
*>         NOTE: RELMAXC2NRMK at the factorization step K would equal
*>               abs(R_sub(K+1,K+1))/MAXC2NRM in the next
*>               factorization step K+1, where R_sub(K+1,K+1) is the
*>               diagonal element of the factor R_sub in the next
*>               factorization step K+1.
*> \endverbatim
*>
*> \param[out] FNRMK
*> \verbatim
*>          FNRMK is REAL
*>          Frobenius norm of the residual matrix
*>          A_sub_resid(K) = A_sub(K+1:M_sub, K+1:N_sub).
*>          FNRMK >= 0.0
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (M)
*>          Row permutation indices due to row deselection,
*>          for 1 <= i <= M.
*>          If IPIV(i) = k, then the row i of A was
*>          the row k of A.
*> \endverbatim
*>
*> \param[out] JPIV
*> \verbatim
*>          JPIV is INTEGER array, dimension (N)
*>          Column permutation indices, for 1 <= j <= N.
*>          If JPIV(j)= k, then the column j of A*P was
*>          the column k of A.
*>
*>          The first K elements of the array JPIV contain
*>          indices of the columns of the factor C that were selected
*>          from the matrix A.
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is REAL array, dimension (min(M_sub,N_sub))
*>          The scalar factors of the elementary reflectors.
*>
*>          If K = 0, all elements TAU(1:min(M_sub,N_sub)) are set
*>             to zero.
*>          If 0 < K <= min(M_sub,N_sub):
*>             only the elements TAU(1:K) may be modified,
*>             the elements TAU(K+1:min(M_sub,N_sub)) are set to zero.
*> \endverbatim
*>
*> \param[out] C
*> \verbatim
*>          C is REAL array.
*>
*>          If FACT = 'P':
*>             the array is not used, the array dimension >= (1,1).
*>
*>          If FACT = 'C':
*>             the array dimension is (LDC,N).
*>             If K = 0:
*>                the M-by-N array C contains a copy of
*>                the original M-by-N matrix A.
*>             If K > 0:
*>              a) columns (1:K) of the array C contain
*>                 the M-by-K factor C (the selected columns
*>                 from the original matrix A).
*>              b) columns (K+1:N) of the array C contain
*>                 the deselected columns from the original
*>                 matrix A.
*>
*>          If FACT = 'X':
*>             the array dimension is (LDC,N).
*>             If K = 0:
*>                the M-by-N array C is not used.
*>             If K > 0:
*>              a) columns (1:K) of the array C contain
*>                 the M-by-K factor C (the selected columns
*>                 from the original matrix A).
*>              b) columns (K+1:N) of the array C are
*>                 not used.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C.
*>          If FACT = 'P', LDC >= 1.
*>          If FACT = 'C' or 'X', LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] QRC
*> \verbatim
*>          QRC is REAL array.
*>
*>          If FACT = 'P' or 'C': The array is not used,
*>                         the array dimension is >= (1,1).
*>
*>          If FACT = 'X': the array dimension is (LDQRC,min(M,N)).
*>
*>             If K = 0, the array is not used.
*>             If K > 0, QRC(1:M,1:K) stores two components from
*>             the QR factorization of the factor C. The K-by-K
*>             factor R is stored in the upper triangle.
*>             The Householder vectors are stored in the lower
*>             trapezoid below the diagonal.
*> \endverbatim
*>
*> \param[in] LDQRC
*> \verbatim
*>          LDQRC is INTEGER
*>          The leading dimension of the array QRC.
*>          If FACT = 'P' or 'C', LDQRC >= 1.
*>          If FACT = 'X', LDQRC >= max(1,M).
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is REAL array.
*>          If FACT = 'P' or 'C': The array is not used,
*>                       the array dimension is >= (1,1).
*>
*>          If FACT = 'X': The array dimension is (LDX,N).
*>            1) If K = 0:
*>                the M-by-N array X contains a copy of
*>                the original M-by-N matrix A.
*>            2) If K > 0:
*>                a) rows (1:K) of the M-by-N array X contain
*>                   the K-by-N factor X, where K <= N.
*>                b) rows (K+1:M) of the M-by-N array X.
*>                   Each column of these rows contains the elements
*>                   whose sum of squares is the residual sum of
*>                   squares for the solution in each column of
*>                   the least squares problem.
*>                   min|| A - C*X ||_F for the unknown X.
*> \endverbatim
*>
*> \param[in] LDX
*> \verbatim
*>          LDX is INTEGER
*>          The leading dimension of the array X.
*>          If FACT = 'P' or 'C', LDX >= 1.
*>          If FACT = 'X', LDX >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (max(1,LWORK)).
*>
*>          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>
*>          Minimal LWORK workspace general requirement.
*>          LWORK >= max( 1, 3*N - 1 ) would be sufficient for all values
*>          of FACT and USESD flags.
*>
*>          For good performance, LWORK should generally be larger, and
*>          the user should query the routine for the optimal LWORK.
*>
*>          If LWORK = -1 or LIWORK =-1 then a workspace query is assumed.
*>          The routine only calculates the optimal size of the WORK and
*>          IWORK arrays, returns these values as the first entry of
*>          the WORK and IWORK  arrays respectively, and no error message
*>          related to LWORK is issued by XERBLA.
*>
*>          Exact minimal workspace requirements.
*>          For USESD = 'N' or 'R' and for all FACT:
*>              LWORK >= max( 1, 3*N - 1 )
*>          For USESD = 'C' or 'A':
*>            a) If FACT = 'P' or 'C':
*>              LWORK >= max( 1, N_sub, min(1,MINMNFREE)*(3*N_free-1) )
*>            b) If FACT = 'X':
*>              LWORK >= max( 1, min(M,N)+N,
*>                               min(1,MINMNFREE)*(3*N_free-1) )
*>          where MINMNFREE = min( M_free, N_free ).
*>
*>          NOTE: The decision, whether the routine uses unblocked
*>          BLAS 2 or blocked BLAS 3 code is based not only on the
*>          dimension LWORK of the available workspace WORK, but
*>          also on:
*>           1a) column preselection stage using SGEQRF:
*>               the optimal block size NB, the crossover point NX
*>               returned by ILAENV for the routine SGEQRF
*>               in comparison to N_sel. (For N_sel <= NX
*>               or N_sel <= NB, unblocked code is used in SGEQRF.)
*>           1b) column preselection stage using SORMQR:
*>               the optimal block size NB returned by ILAENV for
*>               the routine SORMQR in comparison to N_sel. (For
*>               N_sel <= NB, unblocked code is used in SORMQR.)
*>            2) column selection stage via criteria using SGEQRP3RK:
*>               the optimal block size NB, the crossover point NX
*>               returned by ILAENV for the routine SGEQRP3RK
*>               in comparison to min(M,N_sel). (For
*>               min(M_sub, N_free, KMAXFREE) <= NX
*>               or min(M_sub, N_free, KMAXFREE) <= NB, unblocked code
*>               is used in SGEQRP3RK.)
*>           3a) computation of the factor X using SGEQRF in SGELS:
*>               the optimal block size NB, the crossover point NX
*>               returned by ILAENV for the routine SGEQRF
*>               in comparison to K. (For K <= NX or K <= NB,
*>               unblocked code is used in SGEQRF inside SGELS.)
*>           3b) computation of the factor X using SORMQR in SGELS:
*>               the optimal block size NB returned by ILAENV for
*>               the routine SORMQR in comparison to N. (For
*>               N <= NB, unblocked code is used in SORMQR
*>               inside SGELS.)
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (max(1,LIWORK)).
*>
*>          On exit, if INFO >= 0, IWORK(1) returns the optimal LIWORK.
*> \endverbatim
*>
*> \param[in] LIWORK
*> \verbatim
*>          LIWORK is INTEGER
*>          The dimension of the array IWORK.
*>
*>          Minimal LIWORK workspace general requirement.
*>          LIWORK >= max( 1, 2*N ) would be sufficient for all values
*>          of FACT and USESD flags.
*>
*>          The optimal LIWORK is the same as the minimal LIWORK.
*>          The user can still query the routine for the optimal LIWORK.
*>
*>          If LWORK = -1 or LIWORK =-1 then a workspace query is
*>          assumed. The routine only calculates the optimal size of
*>          the WORK and IWORK arrays, returns these values as the first
*>          entry of the WORK and IWORK arrays respectively, and no
*>          error message related to LIWORK is issued by XERBLA.
*>
*>          Exact minimal workspace requirements.
*>          For USESD = 'N' or 'R':
*>           a) If FACT = 'P':
*>            LIWORK >= max( 1, N-1 )
*>           b) If FACT = 'C' or 'X':
*>            LIWORK >= max( 1, 2*N )
*>          For USESD = 'C' or 'A':
*>           a) If FACT = 'P':
*>            LIWORK >= max( 1, (N_free-1) + min(1,N_sel)*N_free )
*>           b) If FACT = 'C' or 'X':
*>            LIWORK >= max( 1, 2*N )
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit.
*>          < 0: if INFO = -i, the i-th argument had an illegal value.
*>          > 0: if INFO =  i, the i-th diagonal element of the
*>               triangular R factor of the QR factorization of
*>               the matrix C is zero. Consequently, C does not have
*>               full rank, and X cannot be computed as the least
*>               squares solution to the overdetermined system C*X = A.
*>               (R is stored in the array QRC.)
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
*> \par Contributors:
*  ==================
*>
*> \verbatim
*>
*>  April 2026, Igor Kozachenko, James Demmel,
*>              EECS Department,
*>              University of California, Berkeley, USA.
*> \endverbatim
*
*> \ingroup gecxx
*
*  =====================================================================
      SUBROUTINE SGECXX( FACT, USESD, M, N,
     $                   DESEL_ROWS, SEL_DESEL_COLS,
     $                   KMAXFREE, ABSTOL, RELTOL, A, LDA,
     $                   K, MAXC2NRMK, RELMAXC2NRMK, FNRMK,
     $                   IPIV, JPIV, TAU, C, LDC, QRC, LDQRC,
     $                   X, LDX, WORK, LWORK, IWORK, LIWORK, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER           FACT, USESD
      INTEGER             INFO, K, KMAXFREE, LDA, LDC, LDQRC,
     $                    LDX, LIWORK, LWORK, M, N
      REAL                ABSTOL, FNRMK, MAXC2NRMK,
     $                    RELMAXC2NRMK, RELTOL
*     ..
*     .. Array Arguments ..
      INTEGER             DESEL_ROWS( * ), IPIV( * ), IWORK( * ),
     $                    JPIV( * ), SEL_DESEL_COLS( * )
      REAL                A( LDA, * ), C( LDC, * ), QRC( LDQRC, * ),
     $                    TAU( * ), WORK( * ), X( LDX, *)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, TWO, MINUSONE
      PARAMETER          ( ZERO = 0.0E+0, TWO = 2.0E+0,
     $                   MINUSONE = -1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, RETURNC, RETURNX,
     $                   USE_DESEL_ROWS, USE_SEL_DESEL_COLS, USETOL
      INTEGER            I, IP, IINFO, ITEMP, J, JDESEL, JP, KFREE,
     $                   KMAXLS, KP0, LIWKMIN, LIWKOPT, LWKMIN,
     $                   LWKOPT, MFREE, MDESEL, MINMN, MINMNFREE,
     $                   MRESID, MSUB, NFREE, NDESEL, NRESID, NSEL,
     $                   NSUB
      REAL               ABSTOLFREE, EPS, MAXC2NRM, MAXC2NRMKFREE,
     $                   RELTOLFREE, RELMAXC2NRMKFREE, SAFMIN

*     .. External Subroutines ..
      EXTERNAL           SCOPY, SGELS, SGEQP3RK, SGEQRF, SLACPY,
     $                   SORMQR, SSWAP, XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            SISNAN, LSAME
      INTEGER            ISAMAX, ILAENV
      REAL               SLAMCH, SLANGE, SNRM2
      EXTERNAL           SISNAN, SLAMCH, SLANGE, SNRM2, ISAMAX,
     $                   ILAENV, LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      MDESEL = 0
      NSEL = 0
      NDESEL = 0
      MSUB = M
      NSUB = N
      MFREE = MSUB
      NFREE = NSUB
      MINMN = MIN( M, N )
*
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
      RETURNX = LSAME( FACT, 'X' )
      RETURNC = LSAME( FACT, 'C' ) .OR. RETURNX
*
      USE_DESEL_ROWS = LSAME( USESD, 'R' )
     $                 .OR. LSAME( USESD, 'A' )
      USE_SEL_DESEL_COLS = LSAME( USESD, 'C' )
     $                 .OR. LSAME( USESD, 'A' )
*
      IF( .NOT.( RETURNC .OR. LSAME( FACT, 'P') ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( USE_DESEL_ROWS .OR. USE_SEL_DESEL_COLS
     $       .OR. LSAME( USESD, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE
*
*        This is to check that the number of preselected columns NSEL
*        cannot be larger than MSUB, which is the number of rows
*        without MDESEL deselected rows. When the number of
*        preselected columns NSEL is larger than MSUB,
*        the factorization of all preselected NSEL columns cannot be
*        completed. MSUB also will be used for LDX argument check
*        later.
*
         IF( USE_DESEL_ROWS ) THEN
*
*           Count the number of free rows MSUB.
*
            DO I = 1, M
               IF( DESEL_ROWS( I ).EQ.-1 ) MDESEL = MDESEL + 1
            END DO
            MSUB = M - MDESEL
            MFREE = MSUB
         END IF
*
         IF( USE_SEL_DESEL_COLS ) THEN
*
*           Count the number of preselected columns NSEL and the
*           number of preselected and free columns NSUB = N - NDESEL.
*
            DO J = 1, N
               IF( SEL_DESEL_COLS( J ).EQ.1 ) NSEL = NSEL + 1
               IF( SEL_DESEL_COLS( J ).EQ.-1 ) NDESEL = NDESEL + 1
            END DO
            NSUB = N - NDESEL
            MFREE = MSUB - NSEL
            NFREE = NSUB - NSEL
*
         END IF
         MINMNFREE = MIN( MFREE, NFREE )
*
         IF( NSEL.GT.MSUB ) THEN
            INFO = -6
         ELSE IF( KMAXFREE.LT.0 ) THEN
            INFO = -7
         ELSE IF( SISNAN( ABSTOL ) ) THEN
            INFO = -8
         ELSE IF( SISNAN( RELTOL ) ) THEN
            INFO = -9
         ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
            INFO = -11
*        This is a check for LDC
         ELSE IF( ( RETURNC .AND. LDC.LT.MAX( 1, M ) )
     $      .OR. ( .NOT.RETURNC .AND. LDC.LT.1 ) ) THEN
            INFO = -20
*        This is a check for LDQRC
         ELSE IF( ( RETURNX .AND. LDQRC.LT.MAX( 1, M ) )
     $      .OR. ( .NOT.RETURNX .AND. LDQRC.LT.1 ) ) THEN
            INFO = -22
*        This is a check for LDX
         ELSE IF( ( RETURNX .AND. LDX.LT.MAX( 1, M ) )
     $      .OR. ( .NOT.RETURNX .AND. LDX.LT.1 ) ) THEN
            INFO = -24
         END IF
*
      END IF
*
*     ==================================================================
*
*       a) Test the input workspace size LWORK and LIWORK for the
*          minimum size requirement LWKMIN and LIWKMIN respectively.
*       b) Determine the optimal workspace sizes LWKOPT and LIWKOPT to
*          be returned in WORK( 1 ) and IWORK( 1 ) respectively,
*          if INFO >= 0 in cases:
*           (1) LQUERY = .TRUE.,
*           (2) when the routine exits.
*     Here, LWKMIN and LIWKMIN are the minimum workspaces required for
*     unblocked code.
*
      IF( INFO.EQ.0 ) THEN
         IF( MINMN.EQ.0 ) THEN
            LWKMIN = 1
            LWKOPT = 1
            LIWKMIN = 1
            LIWKOPT = 1
         ELSE
*
*           (Real_wk_part_1) Real minimum and optimal workspace
*           computation.
*           LWKMIN = MAX(1, NSUB) for column 2-norm computation
*
            LWKMIN = MAX( 1, NSUB )
            LWKOPT = LWKMIN
*
*           (Int_wk_part_1) Integer minimum workspace computation.
*
            LIWKMIN = 1
*
*           Call of SGEQRF.
*
            IF( NSEL.GT.0 ) THEN
*
*              (Real_wk_part_2) Real minimum workspace computation.
*              LWKMIN = MAX(1, NSEL) for the call of SGEQRF.
*              We can skip counting this workspace as
*              LWKMIN = MAX( LWKMIN, NSEL ), since NSEL <= NSUB.
*
*              Query for optimal workspace size for SGEQRF.
*
               CALL SGEQRF( MSUB, NSEL, A, LDA, TAU, WORK,
     $                         -1, IINFO )
               LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
*
*              Call of SORMQR.
*
               IF( NFREE.GT.0 ) THEN
*
*                 (Real_wk_part_3) Real minimum workspace computation.
*                 NOTE: minimum workspace requirement for DORMQR
*                 LWKMIN = MAX(1, NFREE) is smaller than NSUB
*                 and it is smaller than LWKMIN = 3*NFREE-1 for
*                 DGEQP3RK. We can skip counting this workspace as
*                 as LWKMIN = MAX( LWKMIN, NFREE ).).
*
*                 Query for optimal workspace size for SORMQR.
*
                  CALL SORMQR( 'L', 'T', MSUB, NFREE,
     $               NSEL, A, LDA, TAU, A( 1, NSEL+1 ), LDA, WORK,
     $               -1, IINFO )
                  LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
               END IF
*
            END IF
*
*           Call of SGEQP3RK.
*
            IF ( MINMNFREE.NE.0 ) THEN
*
*              (Real_wk_part_4) Real minimum workspace computation.
*              LWKMIN = MAX(1, 3*NFREE-1) for the call of SGEQP3RK.
*
               LWKMIN = MAX( LWKMIN, 3*NFREE - 1 )
*
*              Query for optimal workspace size for SGEQP3RK.
*
               CALL SGEQP3RK( MFREE, NFREE, 0, NFREE,
     $                        MINUSONE, MINUSONE,
     $                        A( 1, 1 ), LDA, KFREE, MAXC2NRMKFREE,
     $                        RELMAXC2NRMKFREE, JPIV( 1 ), TAU( 1 ),
     $                        WORK, -1, IWORK, IINFO )
               LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
*
*              (Int_wk_part_2) Integer minimum workspace computation.
*              LIWKMIN =  NFREE-1 for the call of SGEQP3RK.
*
               LIWKMIN = MAX( LIWKMIN, NFREE-1 )
*
               IF( NSEL.NE.0 ) THEN
*
*                 (Int_wk_part_3) Integer minimum workspace computation.
*                 NFREE is for SGEQP3RK and NFREE-1 for JPIV adjustment.
*
                  LIWKMIN = MAX( LIWKMIN, NFREE + NFREE-1 )
               END IF
*
            END IF
*
            IF( RETURNC ) THEN
*
*              Integer minimum workspace computation.
*              (Int_wk_part_4) LIWKMIN = 2*N for applying the
*              interchanges for the columns in the matrix C.
*
               LIWKMIN = MAX( LIWKMIN, 2*N )
            END IF
*
*           Integer optimal workspace computation.
*
            LIWKOPT = LIWKMIN
*
*           Call of SGELS.
*
            IF( RETURNX ) THEN
*
*              (Real_wk_part_5) Real minimum workspace computation.
*              LWKMIN = max( 1, MINMN + max( MINMN, N ) ) =
*                     = max( 1, MINMN + N ) for the call of SGELS.
*
               LWKMIN = MAX( LWKMIN, MINMN + N )
*
*              Query for optimal workspace size for SGELS.
*
               KMAXLS = MINMN
*
               CALL SGELS( 'N', M, KMAXLS, N, QRC, LDQRC, X, LDX,
     $                     WORK, -1, IINFO )
               LWKOPT = MAX( LWKOPT, INT( WORK(1) ) )
*
            END IF
*
*           End of ELSE for IF( MINMN.EQ.0 )
*
         END IF
*
         IF( ( LWORK.LT.LWKMIN ) .AND. .NOT.LQUERY ) THEN
            INFO = -26
         ELSE IF( ( LIWORK.LT.LIWKMIN ) .AND. .NOT.LQUERY ) THEN
            INFO = -28
         END IF
      END IF
*
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = REAL( LWKOPT )
         IWORK( 1 ) = LIWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGECXX', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     ==================================================================
*
*     Quick return if possible for:
*     a)  M = 0 or N = 0. There is no matrix A(1:M,1:N).
*     b)  MSUB = 0 or NSUB = 0. There is no matrix A_sub(1:MSUB,1:NSUB).
*     NOTE: min( M, N) = 0 implies min( MSUB, NSUB) = 0.
*     We need to return correct values for all scalar output parameters,
*     (including WORK(1) and IWORK(1), which are set above).
*
      IF( MIN( MSUB, NSUB ).EQ.0 ) THEN
         K = 0
         MAXC2NRMK = ZERO
         RELMAXC2NRMK = ZERO
         FNRMK = ZERO
         RETURN
      END IF
*
*     ==================================================================
*
      K = 0
*
*     If we need to return factor X, copy the original untouched matrix
*     A into the array X.
*
      IF( RETURNX ) THEN
         CALL SLACPY( 'F', M, N, A, LDA, X, LDX )
      END IF
*
*     If we need to return the factor C, copy the original matrix A
*     into the array C, only if do not return the factor X. In this
*     case, we need to choose the columns of the matrix A in the array C
*     in place, otherwise we can copy the columns of the matrix A from
*     the array X.
*
      IF( RETURNC .AND. .NOT. RETURNX ) THEN
         CALL SLACPY( 'F', M, N, A, LDA, C, LDC )
      END IF
*
*     ==================================================================
*     Permute the deselected rows to the bottom of the matrix A.
*     1) The initial order of included rows in their block is preserved.
*     2) The initial order of deselected rows in their block is not
*        preserved.
*     ==================================================================
*
*     I is an index of DESEL_ROWS array and a row index of
*     the matrix A. MSUB is the number of processed included rows, which
*     is also an index pointer to the last included row in the matrix A.
*     We can think of I as a row source index, and MSUB as a destination
*     index for moving an included row in the matrix A.
*
*     ( We start with MSUB = 0. We loop over index I in (1:M), and
*     for each position I in DESEL_ROWS  array, we check if the row at
*     the position I in the matrix A is an included row (not -1 value).
*     If it is an included row, we increment MSUB pointer, otherwise
*     we do not change MSUB index pointer. Then, we bring this included
*     row from the index I in the matrix A into smaller (or same)
*     MSUB index in the matrix A.  If I = MSUB, then the included row
*     is already in place. Due to row swap, the deselected row
*     at MSUB index will move into I index in the matrix A. In this way,
*     we move all the included rows to the top matrix block preserving
*     their initial order within the included block. The initial order
*     of deselected rows will not be preserved within their block.
*
      IF( USE_DESEL_ROWS ) THEN
*
         MSUB = 0
         DO I = 1, M, 1
*
*           Initialize the row pivot array IPIV.
            IPIV( I ) = I
*
*           The row at the index I is an included row and should be
*           moved to the top of the matrix A.
*
            IF( DESEL_ROWS( I ).NE.-1 ) THEN
               MSUB = MSUB + 1
*
*              This is a check whether the included row is
*              on the included place already.
*
               IF( I.NE.MSUB ) THEN
*
*                 Here, we swap A(I,1:N) into A(MSUB,1:N).
*
                  CALL SSWAP( N, A( I, 1 ), LDA, A( MSUB, 1 ), LDA )
*
*                 Save the interchange.
*
                  IPIV( I ) = IPIV( MSUB )
                  IPIV( MSUB ) = I
                  DESEL_ROWS( MSUB ) = DESEL_ROWS( I )
                  DESEL_ROWS( I ) = -1
               END IF
            END IF
*
         END DO
*
      ELSE
*
*        We do not use the row deselection DESEL_ROWS array.
*        Initialize the row pivot array IPIV.
*        NOTE: MSUB=M has default value,
*        which is set at the beginning of the routine, before argument
*        checks.
*
         DO I = 1, M, 1
            IPIV( I ) = I
         END DO
      END IF
*
*     ==================================================================
*     Permute the preselected columns to the left and deselected
*     columns to the right of the matrix A.
*     1) The order of preselected columns is preserved.
*     2) The order of free columns is not preserved.
*     3) The order of deselected columns is not preserved.
*     ==================================================================
*
*     J is the index of SEL_DESEL_COLS array and column J
*     of the matrix A.
*
      IF( USE_SEL_DESEL_COLS ) THEN
*
*        Column selection.
*        NSEL is the number of selected columns, also the pointer to
*        the last selected column.
*
         NSEL = 0
         DO J = 1, N, 1
*
*           Initialize column pivot array JPIV.
            JPIV( J ) = J
*
            IF( SEL_DESEL_COLS( J ).EQ.1 ) THEN
               NSEL = NSEL + 1
*
*              This is the check whether the selected column is
*              on the selected place already.
*
               IF( J.NE.NSEL ) THEN
*
*                 Here, we swap the column A(1:M,J) into A(1:M,NSEL)
*
                  CALL SSWAP( M, A( 1, J ), 1, A( 1, NSEL ), 1 )
                  JPIV( J ) = JPIV( NSEL )
                  JPIV( NSEL ) = J
                  SEL_DESEL_COLS( J ) = SEL_DESEL_COLS( NSEL )
                  SEL_DESEL_COLS( NSEL ) = 1
               END IF
            END IF
         END DO
*
*        Column deselection.
*        JDESEL the pointer to the last
*        deselected column counting right-to-left.
*
         JDESEL = N+1
         DO J = N, NSEL+1, -1
            IF( SEL_DESEL_COLS( J ).EQ.-1 ) THEN
               JDESEL = JDESEL - 1
*
*              This is the check whether the deselected column is
*              on the deselected place already.
*
               IF( J.NE.JDESEL ) THEN
*
*                 Here, we swap the column A(1:M,J) into A(1:M,JDESEL)
*
                  CALL SSWAP( M, A( 1, J ), 1, A( 1, JDESEL ), 1 )
                  ITEMP = JPIV( J )
                  JPIV( J ) = JPIV( JDESEL )
                  JPIV( JDESEL ) = ITEMP
                  SEL_DESEL_COLS( J ) = SEL_DESEL_COLS( JDESEL )
                  SEL_DESEL_COLS( JDESEL ) = -1
               END IF
            END IF
         END DO
*
         NSUB = JDESEL - 1
*
      ELSE
*
*        We do not use the column selection deselection
*        SEL_DESEL_COLS array.
*        Initialize column pivot array JPIV.
*        NOTE: NSUB=N has default value,
*        which is set at the beginning of the routine, before argument
*        checks.
*
         DO J = 1, N, 1
            JPIV( J ) = J
         END DO
*
      END IF
*
*     ==================================================================
*     Compute the complete column 2-norms of the submatrix
*     A_sub = A(1:MSUB, 1:NSUB) and store them in WORK(1:NSUB).
*
      DO J = 1, NSUB
         WORK( J ) = SNRM2( MSUB, A( 1, J ), 1 )
      END DO
*
*     Compute the column index of the maximum column 2-norm and
*     the maximum column 2-norm itself for the submatrix
*     A_sub = A(1:MSUB, 1:NSUB).
*
      KP0 = ISAMAX( NSUB, WORK( 1 ), 1 )
      MAXC2NRM = WORK( KP0 )
*
*     ==================================================================
*     Process preselected columns
*
*     Compute the QR factorization of NSEL preselected columns (1:NSEL)
*     in the submatrix A_sub = A(1:MSUB, 1:NSUB) and update
*     remaining NFREE free columns (NSEL+1:NSUB).
*     NSUB = NSEL + NFREE
*
      IF( NSEL.GT.0 ) THEN
*
*           Case (a): MSUB < NSEL.
*
*              This is handled at the argument check stage in the
*              beginning of the routine. When the number of preselected
*              columns is larger than MSUB, hence the factorization of
*              all NSEL columns cannot be completed. Return from the
*              routine with the error of COL_SEL_DESEL parameter.
*
*           Case (b): MSUB = NSEL.
*           Case (c-1): MSUB > NSEL and NSEL = NSUB.
*
*              For cases (b) and (c-1), there will be no residual
*              submatrix  after factorization of NSEL columns
*              at step K = NSEL:
*              A_sub_resid(NSEL) = A(NSEL+1:MSUB, NSEL+1:NSUB).
*
*           Case (c-2): MSUB > NSEL and NSEL < NSUB.
*
*              For Case (c-2) is a submatrix residual at step K=NSEL
*              A_sub_resid(NSEL) = A(NSEL+1:MSUB, NSEL+1:NSUB)
*
         CALL SGEQRF( MSUB, NSEL, A, LDA, TAU, WORK, LWORK, IINFO )
*
*        Apply Q**T from the left to A(NSEL+1:MSUB, NSEL+1:NSUB)
*
         IF( NFREE.GT.0 ) THEN
*
*           This is only for case (c-2) ('L' = Left, 'T' = Transpose)
*
            CALL SORMQR( 'L', 'T', MSUB, NFREE, NSEL,
     $                   A, LDA, TAU, A( 1, NSEL+1 ), LDA, WORK,
     $                   LWORK, IINFO )
         END IF
*
         K = K + NSEL
*
*        End of IF(NSEL.GT.0)
*
      END IF
*
*     ==================================================================
*
      KFREE = 0
*
      IF( MINMNFREE.NE.0 ) THEN
*
*        Factorize NFREE free columns of
*        A_free = A_sub_resid(NSEL) = A(NSEL+1:MSUB, NSEL+1:NSUB),
*        KFREE is the number of columns that were actually factorized
*        among NFREE columns.
*
*     ==================================================================
*
         EPS = SLAMCH('Epsilon')
*
         USETOL = .FALSE.
*
*        Adjust ABSTOL only if nonnegative. Negative value means disabled.
*        We need to keep negative value for later use in criterion
*        check.
*
         IF( ABSTOL.GE.ZERO ) THEN
            SAFMIN = SLAMCH('Safe minimum')
            ABSTOL = MAX( ABSTOL, TWO*SAFMIN )
            USETOL = .TRUE.
         END IF
*
*        Adjust RELTOL only if nonnegative. Negative value means disabled.
*        We need to keep negative value for later use in criterion
*        check.
*
         IF( RELTOL.GE.ZERO ) THEN
            RELTOL = MAX( RELTOL, EPS )
            USETOL = .TRUE.
         END IF
*
*     ==================================================================
*
*        Disable RELTOLFREE when calling SGEQP3RK for free columns
*        factorization, since SGEQP3RK expects RELTOLFREE with respect
*        to the residual matrix A_sub_resid(NSEL), not the whole
*        original matrix A. We can use RELTOL criterion by passing it
*        to ABSTOLFREE as RELTOL*MAXC2NRM. We need to make sure that
*        the negative values of ABSTOL and RELTOL are propagated
*        to ABSTOLFREE and RELTOLFREE, since negative values means
*        that the criterion is disabled.
*
         IF( USETOL ) THEN
            ABSTOLFREE = MAX( ABSTOL, RELTOL * MAXC2NRM )
         ELSE
            ABSTOLFREE = MINUSONE
         END IF
         RELTOLFREE = MINUSONE
*
*        Save JPIV(NSEL+1:NSUB) into WORK(NFREE+1:2*NFREE-1)
*
         IF( NSEL.NE.0 ) THEN
            DO J = 1, NFREE, 1
               IWORK( NFREE + J ) = JPIV( NSEL+J )
            END DO
         END IF
*
         CALL SGEQP3RK( MFREE, NFREE, 0, KMAXFREE,
     $                  ABSTOLFREE, RELTOLFREE,
     $                  A( NSEL+1, NSEL+1 ), LDA, KFREE, MAXC2NRMKFREE,
     $                  RELMAXC2NRMKFREE, JPIV( NSEL+1 ),
     $                  TAU( NSEL+1 ), WORK, LWORK, IWORK, IINFO )
*
*        Adjust JPIV
*
         IF( NSEL.NE.0 ) THEN
            DO J = 1, NFREE, 1
               JPIV( NSEL+J ) = IWORK( NFREE + JPIV( NSEL+J ) )
            END DO
         END IF
*
*        1) Adjust the return value for the number of factorized
*           columns K for the whole submatrix A_sub.
*        2) MAXC2NRMK is returned transparently without change
*           as MAXC2NRMKFREE is returned from SGEQP3RK.
*        3) Adjust the return value RELMAXC2NRMK for the whole
*           submatrix A_sub. We do not use RELMAXC2NRMKFREE
*           returned from SGEQP3RK.
*
         K = K + KFREE
         MAXC2NRMK =  MAXC2NRMKFREE
         RELMAXC2NRMK =  MAXC2NRMK / MAXC2NRM
*
      ELSE
*
*        Set norms to zero
*
         MAXC2NRMK = ZERO
         RELMAXC2NRMK = ZERO
*
      END IF
*
*     Now, MRESID and NRESID is the number of rows and columns
*     respectively in  A_free_resid = A(K+1:MSUB,K+1:NSUB).
*
      MRESID = MFREE-KFREE
      NRESID = NFREE-KFREE
*
      IF( MIN( MRESID, NRESID ).NE.0 ) THEN
         FNRMK = SLANGE( 'F',  MRESID, NRESID, A( K+1, K+1 ),
     $                   LDA, WORK )
      ELSE
         FNRMK = ZERO
      END IF
*
*     ==================================================================
*
*     Return the matrix C.
*
      IF( RETURNC .AND. K.GT.0 ) THEN
*
      IF( RETURNX ) THEN
*
*        Copy the selected K columns of the original matrix A (that was
*        saved into the array X) into the array C according to
*        the pivot array JPIV. If we return X, then the matrix A is
*        saved in the array X, and it is faster to copy into C than
*        doing column permutation in place, as it is the ELSE case.
*
         DO J = 1, K, 1
            CALL SCOPY( M, X( 1, JPIV( J ) ), 1, C( 1, J ), 1 )
         END DO
*
      ELSE
*
*        Swap the columns of the original matrix A copied into
*        the array C in place.
*
*        The original M-by-N matrix A was copied into the array C at
*        the beginning of the routine, if RETURNC = .TRUE..

*        Apply the column permutation matrix P stored in JPIV(1:K)
*        to the columns 1:K in the M-by-N array C in place.
*        After column interchanges, the first K columns of C should
*        be the same as the first K columns of A*P, i.e.
*        (A*P)(1:M,1:K) = C(1:M,1:K). The complexity of this algorithm
*        is min(K,N-1).
*
*        Index I is the original column index in the
*        array C before interchanges.
*        J is the current column index of the original column I at
*        each step of interchanges.
*
*        Auxiliary array IWORK(1:N) stores the inverse P_inv(J)
*        of the current column permutation matrix P(J) at each
*        column interchange step J only for the array
*        values >= J:N.
*        C_prev  = P_inv(J) * C_next.
*        Each IWORK(I) contains JJ corresponding to I
*        Initialize IWORK(1:N) as (1:N).
*
         DO I = 1, N, 1
            IWORK( I ) = I
         END DO
*
*        Auxiliary array IWORK(N+1:2N) stores the current column
*        permutation matrix P_(J) at each column interchange step J
*        only for the array index >= J:N.
*        C_prev * P_(J) = C_next.
*        Each IWORK(N+JJ) contains I corresponding to JJ.
*        Initialize IWORK(N+1:2*N) as (1:N).
*
         DO J = 1, N, 1
            IWORK( N + J ) = J
         END DO
*
*        Loop over the columns J = ( 1:min( K, N-1 ) ) in C.
*
         DO J = 1, MIN( K, N-1 ), 1
*
*           IP is the original pivot column, i.e. is the original
*           column that should be placed in the current column index
*           J in the array C.
*
            IP = JPIV( J )
*
*           I is the original column that is
*           currently in the column index J in the array C after
*           previous column interchanges.
*
            I = IWORK( N+J )
*
            IF( I.NE.IP ) THEN
*
*              JP is the current index of the original pivot
*              column IP in the array C after previous column
*              interchanges.
*
               JP = IWORK( IP )

*              Swap the original pivot column IP = JPIV( J ),
*              at the current pivot index JP = IWORK( IP ) into
*              index J.
*
               CALL SSWAP( M, C( 1, J ), 1, C( 1, JP ), 1 )
*
*              Update the array IWORK(1:N) for the original column
*              I that was swapped with IP.
*
               IWORK( I ) = IWORK( IP )
*
*              Update the array IWORK(N+1:2*N) for the current column
*              index JP that was swapped with the current column
*              index J.
*
               IWORK( N + JP ) = IWORK( N + J )
*
            END IF
*
         END DO
*
*     End of ELSE( RETURNX )
*
      END IF
*
*     End of IF( RETURNC .AND. K.GT.0 )
*
      END IF
*
*     ==================================================================
*
*     Return the matrix X.
*
      IF( RETURNX .AND. K.GT.0 ) THEN
*
*        We need to use C and A to compute X = pseudoinv(C) * A, as
*        the linear least squares solution to the overdetermined system
*        C*X = A. We use LLS routine that uses the QR factorization. For
*        that purpose, we store the matrix C into the array QRC.
*        The matrix A was copied into the array X at the beginning
*        of the routine.
*
         CALL SLACPY( 'F', M, K, C, LDC, QRC, LDQRC )
*
         CALL SGELS( 'N', M, K, N, QRC, LDQRC, X, LDX,
     $               WORK, LWORK, IINFO )
         INFO = IINFO
*
      END IF
*
      WORK( 1 ) = REAL( LWKOPT )
      IWORK( 1 ) = LIWKOPT
*
*     End of SGECXX
*
      END
