#include "common.h"

#ifndef ABS_K
#define ABS_K(a) ((a) > 0 ? (a) : (-(a)))
#endif

#if defined(SKYLAKEX)
#include "dasum_microk_skylakex-2.c"
#elif defined(HASWELL)
#include "dasum_microk_haswell-2.c"
#endif

#ifndef HAVE_DASUM_KERNEL
static FLOAT dasum_kernel(BLASLONG n, FLOAT *x1)
{

    BLASLONG i=0;
    BLASLONG n_8 = n & -8;
    FLOAT *x = x1;
    FLOAT temp0, temp1, temp2, temp3;
    FLOAT temp4, temp5, temp6, temp7;
    FLOAT sum0 = 0.0;
    FLOAT sum1 = 0.0;
    FLOAT sum2 = 0.0;
    FLOAT sum3 = 0.0;
    FLOAT sum4 = 0.0;
    
    while (i < n_8) {
        temp0 = ABS_K(x[0]);
        temp1 = ABS_K(x[1]);
        temp2 = ABS_K(x[2]);
        temp3 = ABS_K(x[3]);
        temp4 = ABS_K(x[4]);
        temp5 = ABS_K(x[5]);
        temp6 = ABS_K(x[6]);
        temp7 = ABS_K(x[7]);
        
        sum0 += temp0;
        sum1 += temp1;
        sum2 += temp2;
        sum3 += temp3;
        
        sum0 += temp4;
        sum1 += temp5;
        sum2 += temp6;
        sum3 += temp7;
        
        x+=8;
        i+=8;
    }

     while (i < n) {
        sum4 += ABS_K(x1[i]);
        i++;
     }

    return sum0+sum1+sum2+sum3+sum4;
}

#endif

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
    BLASLONG i=0;
    FLOAT sumf = 0.0;

    if (n <= 0 || inc_x <= 0) return(sumf);

    if ( inc_x == 1 ) {
        sumf = dasum_kernel(n, x);
    } 
    else {
        n *= inc_x;
       
        while(i < n) {
            sumf += ABS_K(x[i]);
            i += inc_x;
        }
    }
    return(sumf);
}

