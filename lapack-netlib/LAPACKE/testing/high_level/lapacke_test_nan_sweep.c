/*****************************************************************************
 * NaN position sweeps for the dedicated high-level LAPACKE tests.
 * See the LAPACKE_TEST_NAN_SWEEP macro in lapacke_test.h.
 *****************************************************************************/

#include <stdio.h>

#include "lapacke_test.h"

/**
 * \brief Number of buffer positions a rows-by-cols matrix with leading
 * dimension ld occupies in the given layout.
 *
 * \param[in] layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in] rows   Number of matrix rows.
 * \param[in] cols   Number of matrix columns.
 * \param[in] ld     Leading dimension in the given layout.
 * \return The number of buffer positions the rows-by-cols matrix occupies,
 *         including the leading-dimension padding.
 */
size_t lapacke_test_alloc_len(int layout, lapack_int rows, lapack_int cols,
                              lapack_int ld)
{
    if (layout == LAPACK_COL_MAJOR) {
        return (size_t)ld * cols;
    }
    return (size_t)rows * ld;
}

/**
 * \brief Map a buffer position to logical matrix coordinates.
 *
 * \param[in]  layout LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR.
 * \param[in]  rows   Number of matrix rows.
 * \param[in]  cols   Number of matrix columns.
 * \param[in]  ld     Leading dimension in the given layout.
 * \param[in]  p      Buffer position to map.
 * \param[out] i      Row of the position.
 * \param[out] j      Column of the position.
 * \return Nonzero if the position lies inside the logical rows-by-cols
 *         matrix, zero for leading-dimension padding.
 */
int lapacke_test_map_position(int layout, lapack_int rows, lapack_int cols,
                              lapack_int ld, size_t p, lapack_int *i,
                              lapack_int *j)
{
    if (layout == LAPACK_COL_MAJOR) {
        *i = (lapack_int)(p % (size_t)ld);
        *j = (lapack_int)(p / (size_t)ld);
        return *i < rows;
    }

    *j = (lapack_int)(p % (size_t)ld);
    *i = (lapack_int)(p / (size_t)ld);
    return *j < cols;
}

/**
 * \brief Print the detail line for one misreported NaN sweep position.
 *
 * Called by LAPACKE_TEST_NAN_SWEEP; does not count as a check, the sweep
 * records one aggregated check via lapacke_test_sweep_result.
 *
 * \param[in] name     Name of the sweep.
 * \param[in] variant  Variant tag (the layout name).
 * \param[in] i        Row of the misreported position.
 * \param[in] j        Column of the misreported position.
 * \param[in] info     Value returned by the call under test.
 * \param[in] expected Expected value at this position.
 */
void lapacke_test_sweep_mismatch(const char *name, const char *variant,
                                 lapack_int i, lapack_int j, lapack_int info,
                                 lapack_int expected)
{
    printf("     %s (%s): NaN at (%d, %d): info = %6d, expected %6d\n", name,
           variant, (int)i, (int)j, (int)info, (int)expected);
}

/**
 * \brief Record the one aggregated check of a finished NaN sweep.
 *
 * \param[in] name       Name of the sweep.
 * \param[in] variant    Variant tag (the layout name).
 * \param[in] positions  Number of buffer positions swept.
 * \param[in] mismatches Number of misreported positions; zero passes.
 */
void lapacke_test_sweep_result(const char *name, const char *variant,
                               size_t positions, int mismatches)
{
    lapacke_test_checks++;
    if (mismatches != 0) {
        lapacke_test_failures++;
        printf("FAIL %s (%s): %d of %zu NaN positions misreported\n", name,
               variant, mismatches, positions);
    } else {
        printf("ok   %s (%s): swept %zu NaN positions\n", name, variant,
               positions);
    }
}
