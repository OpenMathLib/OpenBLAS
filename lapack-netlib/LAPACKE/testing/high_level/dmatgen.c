#include <string.h>

#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* Self-test for the matrix generation helpers: the routine tests rely on
 * the fills being deterministic and on the structural properties promised
 * by each fill type (also for the fills no routine test uses yet). */

/* Buffer index of logical element (i, j). */
static size_t idx(int layout, lapack_int i, lapack_int j)
{
    return layout == LAPACK_COL_MAJOR ? i + j * (size_t)LD : i * (size_t)LD + j;
}

/* The fills are deterministic, so identical runs must agree bit for bit. */
static int equal_buffers(const double *a, const double *b, size_t len)
{
    return memcmp(a, b, len * sizeof(*a)) == 0;
}

static int is_zero(double v)
{
    return v == 0.0;
}

/* a[i][j] == conj(a[j][i]); plain equality in the real precisions. */
static int conj_equal(double x, double y)
{
    return x == y;
}

static int is_triangular(const double *a, int layout, char uplo)
{
    lapack_int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            const double v = a[idx(layout, i, j)];
            const int outside = (uplo == 'U') ? (i > j) : (i < j);
            if (outside && !is_zero(v)) return 0;
            if (i == j && is_zero(v)) return 0;
        }
    }
    return 1;
}

LAPACKE_TEST(dmatgen)
{
    double a[LD * LD], b[LD * LD];

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];
        const char *lname = lapacke_test_layout_names[l];
        const size_t len_sq = lapacke_test_alloc_len(layout, N, N, LD);
        lapack_int i, j;
        int ok;

        /* The seeds are fixed, so repeated fills of every type must
         * produce bit-identical buffers. */
        ok = 1;
        lapacke_test_dfill(layout, M, N, a, LD);
        lapacke_test_dfill(layout, M, N, b, LD);
        ok = ok &&
             equal_buffers(a, b, lapacke_test_alloc_len(layout, M, N, LD));
        lapacke_test_dfill_rhs(layout, M, NRHS, a, LD);
        lapacke_test_dfill_rhs(layout, M, NRHS, b, LD);
        ok = ok &&
             equal_buffers(a, b, lapacke_test_alloc_len(layout, M, NRHS, LD));
        lapacke_test_dfill_spd(layout, N, a, LD);
        lapacke_test_dfill_spd(layout, N, b, LD);
        ok = ok && equal_buffers(a, b, len_sq);
        lapacke_test_dfill_sym(layout, N, a, LD);
        lapacke_test_dfill_sym(layout, N, b, LD);
        ok = ok && equal_buffers(a, b, len_sq);
        lapacke_test_dfill_tri(layout, 'U', N, a, LD);
        lapacke_test_dfill_tri(layout, 'U', N, b, LD);
        ok = ok && equal_buffers(a, b, len_sq);
        lapacke_test_check("dmatgen fills deterministic", lname, ok ? 0 : -999,
                           0);

        lapacke_test_dfill_sym(layout, N, a, LD);
        ok = 1;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                ok = ok &&
                     conj_equal(a[idx(layout, i, j)], a[idx(layout, j, i)]);
            }
        }
        lapacke_test_check("dmatgen sym is symmetric", lname, ok ? 0 : -999, 0);

        lapacke_test_dfill_tri(layout, 'U', N, a, LD);
        ok = is_triangular(a, layout, 'U');
        lapacke_test_dfill_tri(layout, 'L', N, a, LD);
        ok = ok && is_triangular(a, layout, 'L');
        lapacke_test_check("dmatgen tri is triangular", lname, ok ? 0 : -999,
                           0);
    }

    /* The SPD fill must Cholesky-factorize; the indefinite fill must be
     * rejected with a positive info (not positive definite). */
    lapacke_test_dfill_spd(LAPACK_COL_MAJOR, N, a, LD);
    lapacke_test_check(
        "dmatgen spd is positive definite", NULL,
        API_SUFFIX(LAPACKE_dpotrf)(LAPACK_COL_MAJOR, 'U', N, a, LD), 0);

    lapacke_test_dfill_sym(LAPACK_COL_MAJOR, N, a, LD);
    lapacke_test_check(
        "dmatgen sym is indefinite", NULL,
        API_SUFFIX(LAPACKE_dpotrf)(LAPACK_COL_MAJOR, 'U', N, a, LD) > 0 ? 0
                                                                        : -999,
        0);
}
