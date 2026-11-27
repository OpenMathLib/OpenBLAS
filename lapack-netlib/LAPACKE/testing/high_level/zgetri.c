#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* One allocation test: refill the input, schedule the malloc failure
 * countdown and check the info result of a zgetri call in the indexed
 * layout. */
#define LAPACKE_ZGETRI_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_zfill(layout, N, N, a, LD);                               \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_zgetri)(layout, N, a, LD, ipiv), \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(zgetri)
{
    lapack_complex_double a[LD * LD];
    const lapack_int ipiv[N] = {1, 2, 3};

    /* A holds both the L and the U factor: fully read. */
    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];
        LAPACKE_TEST_ZNAN_SWEEP(
            "zgetri a", l, N, N, a, LD, lapacke_test_region_full, -3,
            (lapacke_test_zfill(layout, N, N, a, LD)),
            API_SUFFIX(LAPACKE_zgetri)(layout, N, a, LD, ipiv));

        /* With NaN checking disabled even all-NaN input must go through to
         * the Fortran routine (valid arguments, so info must not be
         * negative). */
        LAPACKE_set_nancheck(0);
        lapacke_test_zfill_nan(layout, N, N, a, LD);
        lapacke_test_check(
            "zgetri NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_zgetri)(layout, N, a, LD, ipiv) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* The high level allocates the workspace (column-major: the only
     * allocation). */
    LAPACKE_ZGETRI_ALLOC_TEST(0, 0, "zgetri work alloc failure",
                              LAPACK_WORK_MEMORY_ERROR);

    /* Scheduled one past the last column-major allocation: the failure must
     * not fire, so the call must succeed and the count must match exactly. */
    LAPACKE_ZGETRI_ALLOC_TEST(0, 1, "zgetri allocation count", 0);
    lapacke_test_check_alloc_count("zgetri col-major allocation count");

    /* Row-major allocates the workspace, then the transposed copy of A. */
    LAPACKE_ZGETRI_ALLOC_TEST(1, 0, "zgetri work alloc failure",
                              LAPACK_WORK_MEMORY_ERROR);
    LAPACKE_ZGETRI_ALLOC_TEST(1, 1, "zgetri transpose alloc failure (a)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);

    /* Scheduled one past the last row-major allocation: fires if the call
     * allocates more than expected. */
    LAPACKE_ZGETRI_ALLOC_TEST(1, 2, "zgetri allocation count", 0);
    lapacke_test_check_alloc_count("zgetri row-major allocation count");

    /* An invalid matrix_layout must be rejected as an error in argument 1,
     * before any allocation: the scheduled failure must not fire. */
    LAPACKE_ZGETRI_ALLOC_TEST(2, 0, "zgetri invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("zgetri invalid layout allocation count");
}
