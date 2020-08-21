#include "common.h"
#include <math.h>

#if defined(DOUBLE)

#error supports float only

#else

#define ABS fabsf

#endif

#if defined(SKYLAKEX)
#include "sasum_microk_skylakex-2.c"
#elif defined(HASWELL)
#include "sasum_microk_haswell-2.c"
#endif

#ifndef HAVE_KERNEL_32

static FLOAT sasum_kernel_32(BLASLONG n, FLOAT *x1)
{

	BLASLONG i=0;
	FLOAT *x = x1;
	FLOAT temp0, temp1, temp2, temp3;
	FLOAT temp4, temp5, temp6, temp7;
	FLOAT sum0 = 0.0;
	FLOAT sum1 = 0.0;
	FLOAT sum2 = 0.0;
	FLOAT sum3 = 0.0;

	while ( i< n )
	{

		temp0 = ABS(x[0]);
		temp1 = ABS(x[1]);
		temp2 = ABS(x[2]);
		temp3 = ABS(x[3]);
		temp4 = ABS(x[4]);
		temp5 = ABS(x[5]);
		temp6 = ABS(x[6]);
		temp7 = ABS(x[7]);

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

	return sum0+sum1+sum2+sum3;
}

#endif

FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
	BLASLONG i=0;
	FLOAT sumf = 0.0;
	BLASLONG n1;

	if (n <= 0 || inc_x <= 0) return(sumf);

	if ( inc_x == 1 )
	{

		n1 = n & -32;
		if ( n1 > 0 )
		{

			sumf = sasum_kernel_32(n1, x);
			i=n1;
		}

		while(i < n)
		{
			sumf += ABS(x[i]);
			i++;
		}

	}
	else
	{

		n *= inc_x;
		while(i < n)
		{
			sumf += ABS(x[i]);
			i += inc_x;
		}

	}
	return(sumf);
}
