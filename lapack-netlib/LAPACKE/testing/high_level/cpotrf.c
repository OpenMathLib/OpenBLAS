#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* One allocation test: refill the input, schedule the malloc failure
 * countdown and check the info result of a cpotrf call in the indexed
 * layout. */
#define LAPACKE_CPOTRF_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill_spd(layout, N, a, LD);                              \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_cpotrf)(layout, 'U', N, a, LD),  \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(cpotrf)
{
    lapack_complex_float a[LD * LD];

    /* Only the uplo triangle is a documented input; a NaN in the other
     * (never referenced) triangle must be ignored. */
    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        LAPACKE_TEST_CNAN_SWEEP(
            "cpotrf a uplo=U", l, N, N, a, LD, lapacke_test_region_upper, -4,
            (lapacke_test_cfill_spd(layout, N, a, LD)),
            API_SUFFIX(LAPACKE_cpotrf)(layout, 'U', N, a, LD));

        LAPACKE_TEST_CNAN_SWEEP(
            "cpotrf a uplo=L", l, N, N, a, LD, lapacke_test_region_lower, -4,
            (lapacke_test_cfill_spd(layout, N, a, LD)),
            API_SUFFIX(LAPACKE_cpotrf)(layout, 'L', N, a, LD));

        /* With NaN checking disabled the NaN must go through to the Fortran
         * routine, which reports it as a not-positive-definite leading
         * minor (info > 0), never as a LAPACKE argument error. */
        LAPACKE_set_nancheck(0);
        lapacke_test_cfill_nan(layout, N, N, a, LD);
        lapacke_test_check(
            "cpotrf NaN with nancheck off", lapacke_test_layout_names[l],
            API_SUFFIX(LAPACKE_cpotrf)(layout, 'U', N, a, LD) < 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* Column-major neither transposes nor allocates a workspace: the
     * scheduled failure must not fire at all. */
    LAPACKE_CPOTRF_ALLOC_TEST(0, 0, "cpotrf allocation count", 0);
    lapacke_test_check_alloc_count("cpotrf col-major allocation count");

    /* Row-major allocates the transposed copy of A (no workspace). */
    LAPACKE_CPOTRF_ALLOC_TEST(1, 0, "cpotrf transpose alloc failure (a)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);

    /* Scheduled one past the last row-major allocation: fires if the call
     * allocates more than expected. */
    LAPACKE_CPOTRF_ALLOC_TEST(1, 1, "cpotrf allocation count", 0);
    lapacke_test_check_alloc_count("cpotrf row-major allocation count");

    /* An invalid matrix_layout must be rejected as an error in argument 1,
     * before any allocation: the scheduled failure must not fire. */
    LAPACKE_CPOTRF_ALLOC_TEST(2, 0, "cpotrf invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("cpotrf invalid layout allocation count");
}
