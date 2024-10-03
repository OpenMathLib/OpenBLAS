#include "common.h"


#ifdef HAVE_SVE
#if defined __has_include
#if __has_include(<arm_sve.h>) && __ARM_FEATURE_SVE
#define USE_SVE
#endif
#endif
#endif

#include "swap_kernel_sve.c"

//(BLASLONG,  BLASLONG,  BLASLONG,  float,  float *, BLASLONG,  float *, BLASLONG,  float *, BLASLONG)
//int CNAME(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT dummy3, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *dummy, BLASLONG dummy2)
int CNAME(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT  dummy3, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT* dummy, BLASLONG dummy2)
{
      swap_kernel_sve(n, x,inc_x, y, inc_y);
	return 0;

}

