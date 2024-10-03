
#include "common.h"

static int scal_kernel_c(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT da, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *dummy, BLASLONG dummy2)
{
         BLASLONG i=0,j=0;
     
        while(j < n)
        {

                if ( da == 0.0 )
                        x[i]=0.0;
                else
                        x[i] = da * x[i] ;

                i += inc_x ;
                j++;

        }
        return 0;

}

