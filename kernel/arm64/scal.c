#include "common.h"

#include "scal_kernel_sve.c"
#include "scal_kernel_c.c"

int CNAME(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT da, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *dummy, BLASLONG dummy2)
{
        
        if ( (n <= 0) || (inc_x <= 0))
                return(0);
         if (inc_x == 1)
		 scal_kernel_sve( n, x, da);
	 else
		 scal_kernel_c(n,dummy0,dummy1,da,x,inc_x,y,inc_y,dummy,dummy2);

        
        return 0;

}

