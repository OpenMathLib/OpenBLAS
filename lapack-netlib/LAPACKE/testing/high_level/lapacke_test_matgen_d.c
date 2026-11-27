/*****************************************************************************
 * Test matrix generation for the dedicated high-level LAPACKE tests
 * (double precision).
 * See lapacke_test.h.
 *****************************************************************************/

#include <stdio.h>

#include "lapacke_test.h"

/**
 * \brief Generate a random, well-conditioned rows-by-cols matrix.
 *
 * Uses ?LATMS (the generator the LAPACK test suite uses) to generate the
 * matrix in column-major order and scatters it into the target layout;
 * padding positions get a finite sentinel. The seed is fixed per fill type,
 * so every call produces the same matrix and test failures are
 * reproducible. A ?LATMS failure is recorded as a failed check.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  rows   Number of matrix rows (at most LAPACKE_TEST_LD).
 * \param[in]  cols   Number of matrix columns (at most LAPACKE_TEST_LD).
 * \param[out] a      Buffer of lapacke_test_alloc_len(...) doubles.
 * \param[in]  ld     Leading dimension of a in the given layout.
 * \param[in]  sym    ?LATMS SYM parameter: 'N' for a general matrix, 'P'
 *                    for symmetric positive definite, 'S' for symmetric
 *                    with random-sign eigenvalues (indefinite).
 * \param[in]  kl     Lower bandwidth (0 gives an upper triangular matrix).
 * \param[in]  ku     Upper bandwidth (0 gives a lower triangular matrix).
 * \param[in]  seed   The ?LATMS ISEED to start from (not modified; the
 *                    last entry must be odd).
 */
static void lapacke_test_dfill_latms(int layout, lapack_int rows,
                                     lapack_int cols, double *a, lapack_int ld,
                                     char sym, lapack_int kl, lapack_int ku,
                                     const lapack_int seed[4])
{
    double tmp[LAPACKE_TEST_LD * LAPACKE_TEST_LD];
    double d[LAPACKE_TEST_LD];
    double work[3 * LAPACKE_TEST_LD];
    lapack_int iseed[4];
    const char dist = 'U';
    const char pack = 'N';
    const lapack_int mode = 3;
    const double cond = 10.0;
    const double dmax = 1.0;
    const lapack_int ld_tmp = LAPACKE_TEST_LD;
    lapack_int info = 0;
    const size_t len = lapacke_test_alloc_len(layout, rows, cols, ld);
    size_t p;
    lapack_int i, j;

    iseed[0] = seed[0];
    iseed[1] = seed[1];
    iseed[2] = seed[2];
    iseed[3] = seed[3];
    LAPACK_dlatms(&rows, &cols, &dist, iseed, &sym, d, &mode, &cond, &dmax, &kl,
                  &ku, &pack, tmp, &ld_tmp, work, &info);
    if (info != 0) {
        lapacke_test_checks++;
        lapacke_test_failures++;
        printf("FAIL matrix generation: info = %d\n", (int)info);
    }

    for (p = 0; p < len; p++) {
        if (lapacke_test_map_position(layout, rows, cols, ld, p, &i, &j)) {
            a[p] = tmp[i + j * LAPACKE_TEST_LD];
        } else {
            a[p] = 0.25; /* padding, never referenced */
        }
    }
}

/**
 * \brief Fill a with a random general rows-by-cols matrix, double precision.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  rows   Number of matrix rows (at most LAPACKE_TEST_LD).
 * \param[in]  cols   Number of matrix columns (at most LAPACKE_TEST_LD).
 * \param[out] a      Buffer of lapacke_test_alloc_len(...) doubles.
 * \param[in]  ld     Leading dimension of a in the given layout.
 */
void lapacke_test_dfill(int layout, lapack_int rows, lapack_int cols, double *a,
                        lapack_int ld)
{
    static const lapack_int seed[4] = {1988, 1989, 1990, 1991};
    lapacke_test_dfill_latms(layout, rows, cols, a, ld, 'N', rows - 1, cols - 1,
                             seed);
}

/**
 * \brief Fill a with a random symmetric positive definite n-by-n matrix, double
 * precision.
 *
 * The eigenvalues lie in [0.1, 1], so the matrix is safely positive
 * definite for the Cholesky-based tests.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  n      Matrix order (at most LAPACKE_TEST_LD).
 * \param[out] a      Buffer of lapacke_test_alloc_len(...) doubles.
 * \param[in]  ld     Leading dimension of a in the given layout.
 */
void lapacke_test_dfill_spd(int layout, lapack_int n, double *a, lapack_int ld)
{
    static const lapack_int seed[4] = {1990, 1991, 1992, 1993};
    lapacke_test_dfill_latms(layout, n, n, a, ld, 'P', n - 1, n - 1, seed);
}

/**
 * \brief Fill a with a random symmetric indefinite n-by-n matrix, double
 * precision.
 *
 * The eigenvalues have magnitudes in [0.1, 1] and random (seed-fixed)
 * signs, so the matrix is well conditioned but not positive definite --
 * suitable for the symmetric indefinite (sy) factorizations.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  n      Matrix order (at most LAPACKE_TEST_LD).
 * \param[out] a      Buffer of lapacke_test_alloc_len(...) doubles.
 * \param[in]  ld     Leading dimension of a in the given layout.
 */
void lapacke_test_dfill_sym(int layout, lapack_int n, double *a, lapack_int ld)
{
    static const lapack_int seed[4] = {2010, 2011, 2012, 2013};
    lapacke_test_dfill_latms(layout, n, n, a, ld, 'S', n - 1, n - 1, seed);
}

/**
 * \brief Fill a with a random nonsingular triangular n-by-n matrix, double
 * precision.
 *
 * Generated by restricting the bandwidth of a general ?LATMS matrix, so
 * the singular values stay in [0.1, 1]: the matrix is well conditioned
 * and nonsingular, and the entries beyond the uplo triangle are exactly
 * zero. Routines with a unit-diagonal option do not read the diagonal;
 * that is a question for their access region, not for the fill.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  uplo   'U' for upper, 'L' for lower triangular.
 * \param[in]  n      Matrix order (at most LAPACKE_TEST_LD).
 * \param[out] a      Buffer of lapacke_test_alloc_len(...) doubles.
 * \param[in]  ld     Leading dimension of a in the given layout.
 */
void lapacke_test_dfill_tri(int layout, char uplo, lapack_int n, double *a,
                            lapack_int ld)
{
    static const lapack_int seed[4] = {2020, 2021, 2022, 2023};
    const int upper = (uplo == 'U' || uplo == 'u');
    lapacke_test_dfill_latms(layout, n, n, a, ld, 'N', upper ? 0 : n - 1,
                             upper ? n - 1 : 0, seed);
}

/**
 * \brief Fill b with a random right-hand side, double precision.
 *
 * Same generator as lapacke_test_dfill with a different seed, so right-hand
 * sides differ from system matrices.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  rows   Number of matrix rows (at most LAPACKE_TEST_LD).
 * \param[in]  cols   Number of matrix columns (at most LAPACKE_TEST_LD).
 * \param[out] b      Buffer of lapacke_test_alloc_len(...) doubles.
 * \param[in]  ld     Leading dimension of b in the given layout.
 */
void lapacke_test_dfill_rhs(int layout, lapack_int rows, lapack_int cols,
                            double *b, lapack_int ld)
{
    static const lapack_int seed[4] = {2000, 2001, 2002, 2003};
    lapacke_test_dfill_latms(layout, rows, cols, b, ld, 'N', rows - 1, cols - 1,
                             seed);
}

/**
 * \brief Fill every allocated position of a matrix buffer with NaN, double
 * precision.
 *
 * Used to verify that a routine with NaN checking disabled does not reject
 * NaNs anywhere in the buffer.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  rows   Number of matrix rows.
 * \param[in]  cols   Number of matrix columns.
 * \param[out] a      Buffer of lapacke_test_alloc_len(...) doubles.
 * \param[in]  ld     Leading dimension of a in the given layout.
 */
void lapacke_test_dfill_nan(int layout, lapack_int rows, lapack_int cols,
                            double *a, lapack_int ld)
{
    const size_t len = lapacke_test_alloc_len(layout, rows, cols, ld);
    for (size_t p = 0; p < len; p++) {
        a[p] = lapacke_create_nan();
    }
}

/**
 * \brief A quiet NaN, double precision.
 *
 * \return A quiet NaN.
 */
double lapacke_test_dnan(void)
{
    return lapacke_create_nan();
}
