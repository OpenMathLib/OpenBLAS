#include "openblas_utest.h"
#include <cblas.h>

#ifndef NAN
#define NAN 0.0/0.0
#endif
#ifndef INFINITY
#define INFINITY 1.0/0.0
#endif

/* With alpha == 0 the reference GEMM never references A or B, so
 * non-finite values there must not reach C. */

#define GEMM_ALPHA0(CBLAS, TYPE, TA, TB, M, N, K, BETA)             \
    do {                                                                    \
        static TYPE a[64 * 64], b[64 * 64], c[64 * 64];                     \
        int i;                                                              \
        for (i = 0; i < 64 * 64; i++) {                                     \
            a[i] = (i & 1) ? NAN : INFINITY;                                \
            b[i] = (i & 1) ? INFINITY : NAN;                                \
            c[i] = 3.0;                                                     \
        }                                                                   \
        CBLAS(CblasColMajor, TA, TB, M, N, K, 0.0, a, 64, b, 64, BETA,      \
              c, 64);                                                       \
        for (i = 0; i < N; i++) {                                           \
            int j;                                                          \
            for (j = 0; j < M; j++)                                         \
                ASSERT_TRUE(c[i * 64 + j] == (BETA) * 3.0);                 \
        }                                                                   \
    } while (0)

#ifdef BUILD_DOUBLE

CTEST(dgemm, alpha_zero_nonfinite_ab_nn)
{
    GEMM_ALPHA0(cblas_dgemm, double, CblasNoTrans, CblasNoTrans, 17, 19, 23, 1.0);
}

CTEST(dgemm, alpha_zero_nonfinite_ab_nt)
{
    GEMM_ALPHA0(cblas_dgemm, double, CblasNoTrans, CblasTrans, 17, 19, 23, 1.0);
}

CTEST(dgemm, alpha_zero_nonfinite_ab_tn)
{
    GEMM_ALPHA0(cblas_dgemm, double, CblasTrans, CblasNoTrans, 8, 9, 40, 1.0);
}

CTEST(dgemm, alpha_zero_nonfinite_ab_tt)
{
    GEMM_ALPHA0(cblas_dgemm, double, CblasTrans, CblasTrans, 17, 19, 23, 1.0);
}

CTEST(dgemm, alpha_zero_nonfinite_ab_beta2)
{
    GEMM_ALPHA0(cblas_dgemm, double, CblasNoTrans, CblasNoTrans, 17, 19, 23, 2.0);
}

CTEST(dgemm, alpha_zero_nonfinite_ab_beta0)
{
    GEMM_ALPHA0(cblas_dgemm, double, CblasNoTrans, CblasNoTrans, 17, 19, 23, 0.0);
}

#endif

#ifdef BUILD_SINGLE

CTEST(sgemm, alpha_zero_nonfinite_ab_nn)
{
    GEMM_ALPHA0(cblas_sgemm, float, CblasNoTrans, CblasNoTrans, 17, 19, 23, 1.0f);
}

CTEST(sgemm, alpha_zero_nonfinite_ab_nt)
{
    GEMM_ALPHA0(cblas_sgemm, float, CblasNoTrans, CblasTrans, 17, 19, 23, 1.0f);
}

CTEST(sgemm, alpha_zero_nonfinite_ab_tn)
{
    GEMM_ALPHA0(cblas_sgemm, float, CblasTrans, CblasNoTrans, 8, 9, 40, 1.0f);
}

CTEST(sgemm, alpha_zero_nonfinite_ab_tt)
{
    GEMM_ALPHA0(cblas_sgemm, float, CblasTrans, CblasTrans, 17, 19, 23, 1.0f);
}

CTEST(sgemm, alpha_zero_nonfinite_ab_beta2)
{
    GEMM_ALPHA0(cblas_sgemm, float, CblasNoTrans, CblasNoTrans, 17, 19, 23, 2.0f);
}

CTEST(sgemm, alpha_zero_nonfinite_ab_beta0)
{
    GEMM_ALPHA0(cblas_sgemm, float, CblasNoTrans, CblasNoTrans, 17, 19, 23, 0.0f);
}

#endif
