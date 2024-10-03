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
static int scal_kernel_sve(int n, FLOAT *x, FLOAT da)
{
      for (int i = 0; i < n; i += SVE_WIDTH){
        svbool_t pg = SVE_WHILELT(i, n);
        SVE_TYPE  x_vec = svld1(pg, &x[i]);
        SVE_TYPE  result= svmul_z(pg,x_vec,da);
        svst1(pg,&x[i],result);
    }
   return (0);
}

