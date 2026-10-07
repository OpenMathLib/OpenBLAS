#include "openblas_utest.h"
#include <cblas.h>

#ifndef NO_CBLAS

#define SYMM_PACKING(CBLAS, TYPE)                                             \
    do {                                                                     \
        TYPE a[64], b[32], c[32] = {0};                                       \
        int i;                                                               \
        for (i = 0; i < 64; ++i) {                                            \
            a[i] = i / 8 == i % 8 ? 2 : 1;                                   \
        }                                                                    \
        for (i = 0; i < 32; ++i) {                                            \
            b[i] = i + 1;                                                    \
        }                                                                    \
        CBLAS(CblasRowMajor, CblasLeft, CblasLower, 8, 4, 1, a, 8, b, 4,       \
              0, c, 4);                                                     \
        for (i = 0; i < 32; ++i) {                                            \
            ASSERT_TRUE(c[i] == i + 1 + 120 + 8 * (i % 4));                   \
        }                                                                    \
    } while (0)

#define TRMM_PACKING(CBLAS, TYPE)                                             \
    do {                                                                     \
        TYPE a[64], b[32], expected[32] = {0};                                \
        int i, j, k;                                                         \
        for (i = 0; i < 64; ++i) {                                            \
            a[i] = i / 8 < i % 8 ? 0 : i / 8 == i % 8 ? 2 : 1;               \
        }                                                                    \
        for (i = 0; i < 32; ++i) {                                            \
            b[i] = i + 1;                                                    \
        }                                                                    \
        for (i = 0; i < 8; ++i) {                                             \
            for (j = 0; j < 4; ++j) {                                        \
                for (k = 0; k <= i; ++k) {                                   \
                    expected[i * 4 + j] += a[i * 8 + k] * b[k * 4 + j];       \
                }                                                            \
            }                                                                \
        }                                                                    \
        CBLAS(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, \
              8, 4, 1, a, 8, b, 4);                                         \
        for (i = 0; i < 32; ++i) {                                            \
            ASSERT_TRUE(b[i] == expected[i]);                                 \
        }                                                                    \
    } while (0)

#ifdef BUILD_DOUBLE
CTEST(matrix_packing, dsymm)
{
    SYMM_PACKING(cblas_dsymm, double);
}

CTEST(matrix_packing, dtrmm)
{
    TRMM_PACKING(cblas_dtrmm, double);
}
#endif

#ifdef BUILD_SINGLE
CTEST(matrix_packing, ssymm)
{
    SYMM_PACKING(cblas_ssymm, float);
}

CTEST(matrix_packing, strmm)
{
    TRMM_PACKING(cblas_strmm, float);
}
#endif

#endif
