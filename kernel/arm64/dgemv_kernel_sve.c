#include "common.h"

#include <arm_sve.h>

#ifdef DOUBLE
#define SVE_TYPE svfloat64_t
#define SVE_ZERO svdup_f64(0.0)
#define SVE_WHILELT svwhilelt_b64
#define SVE_ALL svptrue_b64()
#define SVE_WIDTH svcntd()
#else
#define SVE_TYPE svfloat32_t
#define SVE_ZERO svdup_f32(0.0)
#define SVE_WHILELT svwhilelt_b32
#define SVE_ALL svptrue_b32()
#define SVE_WIDTH svcntw()
#endif

static  FLOAT dgemv_kernel_sve(BLASLONG i, FLOAT  *x, BLASLONG lda, FLOAT  *y, BLASLONG incx, BLASLONG  n){
        SVE_TYPE acc_a = SVE_ZERO;
        SVE_TYPE acc_b = SVE_ZERO;

        BLASLONG  sve_width = SVE_WIDTH;

        for (BLASLONG j = 0; j < n; j += sve_width * 2) {
                svbool_t pg_a = SVE_WHILELT(j, n);
                svbool_t pg_b = SVE_WHILELT(j + sve_width, n);

                SVE_TYPE x_vec_a = svld1(pg_a, &x[i*lda+j]);
                SVE_TYPE y_vec_a = svld1(pg_a, &y[j*incx]);
                SVE_TYPE x_vec_b = svld1(pg_b, &x[i*lda+j + sve_width]);
                SVE_TYPE y_vec_b = svld1(pg_b, &y[j*incx + sve_width]);

                acc_a = svmla_m(pg_a, acc_a, x_vec_a, y_vec_a);
                acc_b = svmla_m(pg_b, acc_b, x_vec_b, y_vec_b);
        }

        return svaddv(SVE_ALL, acc_a) + svaddv(SVE_ALL, acc_b);

}

