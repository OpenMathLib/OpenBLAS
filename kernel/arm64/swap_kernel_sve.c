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

static int  swap_kernel_sve(BLASLONG n, FLOAT *x,BLASLONG inc_x, FLOAT *y, BLASLONG inc_y) {
          BLASLONG sve_width = SVE_WIDTH;

        for (BLASLONG i = 0; i < n; i += sve_width * 2) {
                svbool_t pg_a = SVE_WHILELT(i, n);
                svbool_t pg_b = SVE_WHILELT((i + sve_width), n);
                SVE_TYPE x_vec_a = svld1(pg_a, &x[i]);
                SVE_TYPE y_vec_a = svld1(pg_a, &y[i]);
                SVE_TYPE x_vec_b = svld1(pg_b, &x[i + sve_width]);
                SVE_TYPE y_vec_b = svld1(pg_b, &y[i + sve_width]);

        svst1(pg_a, &x[i], y_vec_a);
        svst1(pg_a, &y[i], x_vec_a);
	svst1(pg_b, &x[i+sve_width], y_vec_b);
        svst1(pg_b, &y[i+sve_width], x_vec_b);
             }
	return 0;
}

