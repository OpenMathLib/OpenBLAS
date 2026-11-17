/*****************************************************************************
 * Check reporting for the dedicated high-level LAPACKE tests.
 * See lapacke_test.h.
 *****************************************************************************/

#include <stdio.h>

#include "lapacke_test.h"

int lapacke_test_checks = 0;
int lapacke_test_failures = 0;

const int lapacke_test_layouts[3] = {LAPACK_COL_MAJOR, LAPACK_ROW_MAJOR,
                                     LAPACKE_TEST_INVALID_LAYOUT};
const char *lapacke_test_layout_names[3] = {"col-major", "row-major",
                                            "invalid-layout"};

/**
 * \brief A quiet NaN produced at run time (not folded by the compiler).
 *
 * The zero is read from a volatile object so that the division cannot be
 * folded (and warned about) at compile time.
 *
 * \return A quiet NaN.
 */
double lapacke_create_nan(void)
{
    static volatile double zero = 0.0;
    return 0.0 / zero;
}

/**
 * \brief Record one check of an info result against its expected value.
 *
 * Prints one report line, counts the check and counts the failure if info
 * and expected differ.
 *
 * \param[in] name     Name of the check, printed on the report line.
 * \param[in] variant  Optional variant tag (e.g. the layout name); may be
 *                     NULL.
 * \param[in] info     Value returned by the call under test.
 * \param[in] expected Expected value.
 */
void lapacke_test_check(const char *name, const char *variant, lapack_int info,
                        lapack_int expected)
{
    const char *status = (info == expected) ? "ok  " : "FAIL";
    lapacke_test_checks++;
    if (info != expected) {
        lapacke_test_failures++;
    }
    if (variant != NULL) {
        printf("%s %s (%s): info = %6d, expected %6d\n", status, name, variant,
               (int)info, (int)expected);
    } else {
        printf("%s %s: info = %6d, expected %6d\n", status, name, (int)info,
               (int)expected);
    }
}
