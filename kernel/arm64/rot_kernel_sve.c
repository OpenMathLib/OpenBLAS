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

static void  rot_kernel_sve(BLASLONG n, FLOAT *x, FLOAT *y, FLOAT c, FLOAT s){
       
       for(int i=0; i<n; i+=SVE_WIDTH){     
       	svbool_t pg =SVE_WHILELT((uint32_t)i,(uint32_t) n);
   
        SVE_TYPE  x_vec = svld1(pg, &x[i]);
        SVE_TYPE  y_vec = svld1(pg, &y[i]);
        SVE_TYPE  cx_vec=svmul_z(pg,x_vec,c);
        SVE_TYPE  sy_vec=svmul_z(pg,y_vec,s);

        SVE_TYPE  sx_vec=svmul_z(pg,x_vec,s);
        SVE_TYPE  cy_vec=svmul_z(pg,y_vec,c);

       svst1(pg,&x[i],svadd_z(pg,cx_vec,sy_vec));
       svst1(pg,&y[i],svsub_z(pg,cy_vec,sx_vec));
	
        }
}
    

