#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* One allocation test: refill the input, schedule the malloc failure
 * countdown and check the info result of a zgetrf call in the indexed
 * layout. */
#define LAPACKE_ZGETRF_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill(layout, M, N, a, LD);                               \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(                                                    \
            name, lapacke_test_layout_names[layout_index],                     \
            API_SUFFIX(LAPACKE_zgetrf)(layout, M, N, a, LD, ipiv), expected);  \
    } while (0)

LAPACKE_TEST(zgetrf)
{
    lapack_complex_double a[LD * LD];
    lapack_int ipiv[M];

    /* The whole M-by-N matrix is a documented input. */
    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];
        LAPACKE_TEST_ZNAN_SWEEP(
            "zgetrf a", l, M, N, a, LD, lapacke_test_region_full, -4,
            (lapacke_test_zfill(layout, M, N, a, LD)),
            API_SUFFIX(LAPACKE_zgetrf)(layout, M, N, a, LD, ipiv));

        /* With NaN checking disabled even all-NaN input must go through to
         * the Fortran routine (valid arguments, so info must not be
         * negative). */
        LAPACKE_set_nancheck(0);
        lapacke_test_zfill_nan(layout, M, N, a, LD);
        lapacke_test_check(
            "zgetrf NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zgetrf)(layout, M, N, a, LD, ipiv) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* Column-major neither transposes nor allocates a workspace: the
     * scheduled failure must not fire at all. */
    LAPACKE_ZGETRF_ALLOC_TEST(0, 0, "zgetrf allocation count", 0);
    lapacke_test_check_alloc_count("zgetrf col-major allocation count");

    /* Row-major allocates the transposed copy of A (no workspace). */
    LAPACKE_ZGETRF_ALLOC_TEST(1, 0, "zgetrf transpose alloc failure (a)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);

    /* Scheduled one past the last row-major allocation: fires if the call
     * allocates more than expected. */
    LAPACKE_ZGETRF_ALLOC_TEST(1, 1, "zgetrf allocation count", 0);
    lapacke_test_check_alloc_count("zgetrf row-major allocation count");

    /* An invalid matrix_layout must be rejected as an error in argument 1,
     * before any allocation: the scheduled failure must not fire. */
    LAPACKE_ZGETRF_ALLOC_TEST(2, 0, "zgetrf invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zgetrf invalid layout allocation count");
}
