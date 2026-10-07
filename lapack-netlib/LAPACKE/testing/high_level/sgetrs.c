#include "lapacke_test.h"

#define N LAPACKE_TEST_N
#define NRHS LAPACKE_TEST_NRHS
#define LD LAPACKE_TEST_LD

/* One allocation test: refill the inputs, schedule the malloc failure
 * countdown and check the info result of a sgetrs call in the indexed
 * layout. */
#define LAPACKE_SGETRS_ALLOC_TEST(layout_index, countdown, name, expected)     \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_sfill(layout, N, N, a, LD);                               \
        lapacke_test_sfill_rhs(layout, N, NRHS, b, LD);                        \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           API_SUFFIX(LAPACKE_sgetrs)(layout, 'N', N, NRHS, a, \
                                                      LD, ipiv, b, LD),        \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(sgetrs)
{
    float a[LD * LD], b[LD * LD];
    const lapack_int ipiv[N] = {1, 2, 3};

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];

        /* A holds both the L and the U factor: fully read. */
        LAPACKE_TEST_SNAN_SWEEP(
            "sgetrs a", l, N, N, a, LD, lapacke_test_region_full, -5,
            (lapacke_test_sfill(layout, N, N, a, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_sgetrs)(layout, 'N', N, NRHS, a, LD, ipiv, b,
                                       LD));

        LAPACKE_TEST_SNAN_SWEEP(
            "sgetrs b", l, N, NRHS, b, LD, lapacke_test_region_full, -8,
            (lapacke_test_sfill(layout, N, N, a, LD),
             lapacke_test_sfill_rhs(layout, N, NRHS, b, LD)),
            API_SUFFIX(LAPACKE_sgetrs)(layout, 'N', N, NRHS, a, LD, ipiv, b,
                                       LD));

        /* With NaN checking disabled even all-NaN input must go through to
         * the Fortran routine (valid arguments, so info must not be
         * negative). */
        LAPACKE_set_nancheck(0);
        lapacke_test_sfill_nan(layout, N, N, a, LD);
        lapacke_test_sfill_nan(layout, N, NRHS, b, LD);
        lapacke_test_check("sgetrs NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           API_SUFFIX(LAPACKE_sgetrs)(layout, 'N', N, NRHS, a,
                                                      LD, ipiv, b, LD) < 0,
                           0);
        LAPACKE_set_nancheck(1);
    }

    /* Column-major neither transposes nor allocates a workspace: the
     * scheduled failure must not fire at all. */
    LAPACKE_SGETRS_ALLOC_TEST(0, 0, "sgetrs allocation count", 0);
    lapacke_test_check_alloc_count("sgetrs col-major allocation count");

    /* Row-major allocates the transposed copies of A, then B. */
    LAPACKE_SGETRS_ALLOC_TEST(1, 0, "sgetrs transpose alloc failure (a)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);
    LAPACKE_SGETRS_ALLOC_TEST(1, 1, "sgetrs transpose alloc failure (b)",
                              LAPACK_TRANSPOSE_MEMORY_ERROR);

    /* Scheduled one past the last row-major allocation: fires if the call
     * allocates more than expected. */
    LAPACKE_SGETRS_ALLOC_TEST(1, 2, "sgetrs allocation count", 0);
    lapacke_test_check_alloc_count("sgetrs row-major allocation count");

    /* An invalid matrix_layout must be rejected as an error in argument 1,
     * before any allocation: the scheduled failure must not fire. */
    LAPACKE_SGETRS_ALLOC_TEST(2, 0, "sgetrs invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("sgetrs invalid layout allocation count");
}
