#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* One allocation test: refill the inputs, schedule the malloc failure
 * countdown and check the info result of a dgels call in the indexed
 * layout. */
#define LAPACKE_DGELS_ALLOC_TEST(layout_index, countdown, name, expected)      \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_dfill(layout, M, N, a, LD);                               \
        lapacke_test_dfill_rhs(layout, M, NRHS, b, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_dgels)(layout, 'N', M, N, NRHS, a, LD, b, LD),  \
            expected);                                                         \
    } while (0)

LAPACKE_TEST(dgels)
{
    double a[LD * LD], b[LD * LD];

    /* B has max(M, N) rows; with M > N and trans = 'N' all of them are
     * documented inputs. (For M < N only the first M rows would be inputs
     * while LAPACKE checks all max(M, N) rows -- a known divergence kept
     * out of the spike.) */
    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_DNAN_SWEEP(
            "dgels a", l, M, N, a, LD, lapacke_test_region_full, -6,
            (lapacke_test_dfill(layout, M, N, a, LD),
             lapacke_test_dfill_rhs(layout, M, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dgels)(layout, 'N', M, N, NRHS, a, LD, b, LD));

        LAPACKE_TEST_DNAN_SWEEP(
            "dgels b", l, M, NRHS, b, LD, lapacke_test_region_full, -8,
            (lapacke_test_dfill(layout, M, N, a, LD),
             lapacke_test_dfill_rhs(layout, M, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_dgels)(layout, 'N', M, N, NRHS, a, LD, b, LD));

        /* With NaN checking disabled even all-NaN input must go through to
         * the Fortran routine (valid arguments, so info must not be
         * negative). */
        LAPACKE_set_nancheck(0);
        lapacke_test_dfill_nan(layout, M, N, a, LD);
        lapacke_test_dfill_nan(layout, M, NRHS, b, LD);
        lapacke_test_check("dgels NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_dgels)(layout, 'N', M, N, NRHS, a,
                                                     LD, b, LD) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* Column-major: the internal workspace is the only allocation. */
    LAPACKE_DGELS_ALLOC_TEST(0, 0, "dgels work alloc failure",
                             LAPACK_WORK_MEMORY_ERROR);

    /* Scheduled one past the last column-major allocation: the failure must
     * not fire, so the call must succeed and the count must match exactly. */
    LAPACKE_DGELS_ALLOC_TEST(0, 1, "dgels allocation count", 0);
    lapacke_test_check_alloc_count("dgels col-major allocation count");

    /* Row-major: allocation order is workspace, then the transposed copies
     * of A and B inside the middle level. */
    LAPACKE_DGELS_ALLOC_TEST(1, 0, "dgels work alloc failure",
                             LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_DGELS_ALLOC_TEST(1, 1, "dgels transpose alloc failure (a)",
                             LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_DGELS_ALLOC_TEST(1, 2, "dgels transpose alloc failure (b)",
                             LAPACK_TRANSPOSE_MEMORY_ERROR);

    /* Scheduled one past the last row-major allocation: fires if the call
     * allocates more than expected. */
    LAPACKE_DGELS_ALLOC_TEST(1, 3, "dgels allocation count", 0);
    lapacke_test_check_alloc_count("dgels row-major allocation count");

    /* An invalid matrix_layout must be rejected as an error in argument 1,
     * before any allocation: the scheduled failure must not fire. */
    LAPACKE_DGELS_ALLOC_TEST(2, 0, "dgels invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("dgels invalid layout allocation count");
}
