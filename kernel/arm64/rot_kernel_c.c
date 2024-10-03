#include "common.h"

static int rot_kernel_c(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT c, FLOAT s)
{
        BLASLONG i=0;
        BLASLONG ix=0,iy=0;
        FLOAT temp;

        if ( n <= 0     )  return(0);

        while(i < n)
        {
                temp   = c*x[ix] + s*y[iy] ;
                y[iy]  = c*y[iy] - s*x[ix] ;
                x[ix]  = temp ;

                ix += inc_x ;
                iy += inc_y ;
                i++ ;

        }
        return(0);

}

