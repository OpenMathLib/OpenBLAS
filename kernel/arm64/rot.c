#include "common.h"

#include "rot_kernel_sve.c"
#include "rot_kernel_c.c"

int CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT c, FLOAT s)
{
  
       
        if ( n <= 0 )  return(0);
        if ( inc_x == 1 && inc_y==1)
	   rot_kernel_sve( n, x, y, c, s);
        else 
	   rot_kernel_c ( n, x, inc_x, y, inc_y, c, s);	
       
       	return(0);

}

