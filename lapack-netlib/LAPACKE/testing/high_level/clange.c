#include "lapacke_test.h"

#define M LAPACKE_TEST_M
#define N LAPACKE_TEST_N
#define LD LAPACKE_TEST_LD

/* clange is a value-returning function; on NaN detection it returns the
 * (negative) info code as its value. Map a plausible norm (finite,
 * nonnegative) to 0 so the sweep can compare info codes. A NaN result
 * means the routine consumed the swept NaN instead of rejecting it; map
 * it to a sentinel instead of casting (casting NaN to an integer is
 * undefined behavior and yields 0 on AArch64, masking the failure). */
static lapack_int clange_info(float value)
{
    if (value != value) {
        return -999;
    }
    if (value >= 0.0f) {
        return 0;
    }
    return (lapack_int)value;
}

/* One allocation test: refill the input, schedule the malloc failure
 * countdown and check the mapped result of a clange call with the given
 * norm in the indexed layout. */
#define LAPACKE_CLANGE_ALLOC_TEST(layout_index, countdown, norm, name,         \
                                  expected)                                    \
    do {                                                                       \
        const int layout = lapacke_test_layouts[layout_index];                 \
        lapacke_test_cfill(layout, M, N, a, LD);                               \
        lapacke_test_schedule_malloc_failure(countdown);                       \
        lapacke_test_check(name, lapacke_test_layout_names[layout_index],      \
                           clange_info(API_SUFFIX(LAPACKE_clange)(             \
                               layout, norm, M, N, a, LD)),                    \
                           expected);                                          \
    } while (0)

LAPACKE_TEST(clange)
{
    lapack_complex_float a[LD * LD];
    float res;

    for (size_t l = 0; l < 2; l++) {
        const int layout = lapacke_test_layouts[l];
        LAPACKE_TEST_CNAN_SWEEP(
            "clange a", l, M, N, a, LD, lapacke_test_region_full, -5,
            (lapacke_test_cfill(layout, M, N, a, LD)),
            clange_info(API_SUFFIX(LAPACKE_clange)(layout, '1', M, N, a, LD)));

        /* With NaN checking disabled the NaN must not be rejected with
         * -5. */
        LAPACKE_set_nancheck(0);
        lapacke_test_cfill_nan(layout, M, N, a, LD);
        res = API_SUFFIX(LAPACKE_clange)(layout, '1', M, N, a, LD);
        lapacke_test_check("clange NaN with nancheck off",
                           lapacke_test_layout_names[l],
                           clange_info(res) == -5 ? -5 : 0, 0);
        LAPACKE_set_nancheck(1);
    }

    /* The infinity norm allocates a work array in the high level. On
     * allocation failure LAPACKE_clange returns 0.0 -- indistinguishable
     * from a genuine zero norm; checked exactly here (the macro's mapping
     * would accept any plausible norm). Pinned as the current behavior;
     * questionable enough to raise upstream. */
    lapacke_test_cfill(LAPACK_COL_MAJOR, M, N, a, LD);
    lapacke_test_schedule_malloc_failure(0);
    res = API_SUFFIX(LAPACKE_clange)(LAPACK_COL_MAJOR, 'I', M, N, a, LD);
    lapacke_test_check("clange work alloc failure returns 0.0", "col-major",
                       res == 0.0f ? 0 : -999, 0);

    /* Scheduled one past the last column-major allocation: the failure must
     * not fire, so the call must succeed and the count must match exactly. */
    LAPACKE_CLANGE_ALLOC_TEST(0, 1, 'I', "clange allocation count", 0);
    lapacke_test_check_alloc_count("clange col-major allocation count");

    /* An invalid matrix_layout must be rejected as an error in argument 1,
     * before any allocation: the scheduled failure must not fire. */
    LAPACKE_CLANGE_ALLOC_TEST(2, 0, '1', "clange invalid matrix_layout", -1);
    lapacke_test_check_alloc_count("clange invalid layout allocation count");
}
