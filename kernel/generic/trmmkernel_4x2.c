#include "common.h"
#include <stdbool.h>

int CNAME(BLASLONG bm,BLASLONG bn,BLASLONG bk,FLOAT alpha,FLOAT* ba,FLOAT* bb,FLOAT* C,BLASLONG ldc ,BLASLONG offset)
{

   BLASLONG i,j,k;
   FLOAT *C0,*C1,*ptrba,*ptrbb;

   FLOAT res0_0;
   FLOAT res0_1;
   FLOAT res0_2;
   FLOAT res0_3;

   FLOAT res1_0;
   FLOAT res1_1;
   FLOAT res1_2;
   FLOAT res1_3;

   FLOAT a0;
   FLOAT a1;

   FLOAT b0;
   FLOAT b1;

   BLASLONG off, temp;

   bool left;
   bool transposed;
   bool backwards;

#ifdef LEFT
   left = true;
#else
   left = false;
#endif

#ifdef TRANSA
   transposed = true;
#else
   transposed = false;
#endif

   backwards = left != transposed;

   if (!left) {
      off = -offset;
   }

   for (j=0; j<(bn/2); j+=2) // do the Mx2 loops 
   {
        C0 = C;
        C1 = C0+ldc;

#if defined(TRMMKERNEL) && defined(LEFT)
		off = offset;
#endif


        ptrba = ba;

        for (i=0; i<bm/4; i+=1) // do blocks of 4x2
	{

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		ptrbb = bb;
#else
		ptrba += off*4;
		ptrbb = bb + off*2;
#endif

		res0_0 = 0;
		res0_1 = 0;
		res0_2 = 0;
		res0_3 = 0;

		res1_0 = 0;
		res1_1 = 0;
		res1_2 = 0;
		res1_3 = 0;


#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
		temp = bk-off;
#elif defined(LEFT)
		temp = off+4;	// number of values in A
#else
		temp = off+2;	// number of values in B
#endif

		for (k=0; k<temp; k++)
                {
			b0 = ptrbb[0];
			b1 = ptrbb[1];

			a0 = ptrba[0];
			res0_0 += a0*b0;
			res1_0 += a0*b1;

			a1 = ptrba[1];
			res0_1 += a1*b0;
			res1_1 += a1*b1;

			a0 = ptrba[2];
			res0_2 += a0*b0;
			res1_2 += a0*b1;

			a1 = ptrba[3];
			res0_3 += a1*b0;
			res1_3 += a1*b1;

			ptrba = ptrba+4;
			ptrbb = ptrbb+2;
                }

		res0_0 *= alpha;
		res0_1 *= alpha;
		res0_2 *= alpha;
		res0_3 *= alpha;

		res1_0 *= alpha;
		res1_1 *= alpha;
		res1_2 *= alpha;
		res1_3 *= alpha;

		C0[0] = res0_0;
		C0[1] = res0_1;
		C0[2] = res0_2;
		C0[3] = res0_3;

		C1[0] = res1_0;
		C1[1] = res1_1;
		C1[2] = res1_2;
		C1[3] = res1_3;


#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		temp = bk - off;
#ifdef LEFT
		temp -= 4; // number of values in A
#else
		temp -= 2; // number of values in B
#endif
		ptrba += temp*4;
		ptrbb += temp*2;
#endif

#ifdef LEFT
		off += 4; // number of values in A
#endif

		C0 = C0+4;
		C1 = C1+4;

	}

	if ( bm & 2 ) // do any 2x2 loop
	{

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		ptrbb = bb;
#else
		ptrba += off*2;
		ptrbb = bb + off*2;
#endif

		res0_0 = 0;
		res0_1 = 0;

		res1_0 = 0;
		res1_1 = 0;


#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
		temp = bk-off;
#elif defined(LEFT)
		temp = off+2;	// number of values in A
#else
		temp = off+2;	// number of values in B
#endif

		for (k=0; k<temp; k++)
                {
			b0 = ptrbb[0];
			b1 = ptrbb[1];

			a0 = ptrba[0];
			res0_0 += a0*b0;
			res1_0 += a0*b1;

			a1 = ptrba[1];
			res0_1 += a1*b0;
			res1_1 += a1*b1;

			ptrba = ptrba+2;
			ptrbb = ptrbb+2;
                }

		res0_0 *= alpha;
		res0_1 *= alpha;

		res1_0 *= alpha;
		res1_1 *= alpha;

		C0[0] = res0_0;
		C0[1] = res0_1;

		C1[0] = res1_0;
		C1[1] = res1_1;


#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		temp = bk - off;
#ifdef LEFT
		temp -= 2; // number of values in A
#else
		temp -= 2; // number of values in B
#endif
		ptrba += temp*2;
		ptrbb += temp*2;
#endif

#ifdef LEFT
		off += 2; // number of values in A
#endif

		C0 = C0+2;
		C1 = C1+2;

	}

	if ( bm & 1 ) // do any 1x2 loop
	{

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		ptrbb = bb;
#else
		ptrba += off*1;
		ptrbb = bb + off*2;
#endif

		res0_0 = 0;

		res1_0 = 0;


#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
		temp = bk-off;
#elif defined(LEFT)
		temp = off+1;	// number of values in A
#else
		temp = off+2;	// number of values in B
#endif

		for (k=0; k<temp; k++)
                {
			b0 = ptrbb[0];
			b1 = ptrbb[1];

			a0 = ptrba[0];
			res0_0 += a0*b0;
			res1_0 += a0*b1;

			ptrba = ptrba+1;
			ptrbb = ptrbb+2;
                }

		res0_0 *= alpha;

		res1_0 *= alpha;

		C0[0] = res0_0;

		C1[0] = res1_0;


#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		temp = bk - off;
#ifdef LEFT
		temp -= 1; // number of values in A
#else
		temp -= 2; // number of values in B
#endif
		ptrba += temp*1;
		ptrbb += temp*2;
#endif

#ifdef LEFT
		off += 1; // number of values in A
#endif

		C0 = C0+1;
		C1 = C1+1;

	}


#if defined(TRMMKERNEL) && !defined(LEFT)
		off += 2;
#endif

        k = (bk<<1);
        bb = bb+k;
        i = (ldc<<1);
        C = C+i;
    }







   for (j=0; j<(bn&1); j+=1) // do the Mx1 loops
   {
        C0 = C;

#if defined(TRMMKERNEL) &&  defined(LEFT)
	off = offset;
#endif

        ptrba = ba;

        for (i=0; i<bm/4; i+=1) // do blocks of 4x1 loops
	{

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		ptrbb = bb;
#else
		ptrba += off*4;
		ptrbb = bb + off*1;
#endif

		res0_0 = 0;
		res0_1 = 0;
		res0_2 = 0;
		res0_3 = 0;


#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
		temp = bk-off;
#elif defined(LEFT)
		temp = off+4;	// number of values in A
#else
		temp = off+1;	// number of values in B
#endif

		for (k=0; k<temp; k++)
                {
			b0 = ptrbb[0];

			a0 = ptrba[0];
			res0_0 += a0*b0;

			a1 = ptrba[1];
			res0_1 += a1*b0;

			a0 = ptrba[2];
			res0_2 += a0*b0;

			a1 = ptrba[3];
			res0_3 += a1*b0;

			ptrba = ptrba+4;
			ptrbb = ptrbb+1;
                }

		res0_0 *= alpha;
		res0_1 *= alpha;
		res0_2 *= alpha;
		res0_3 *= alpha;

		C0[0] = res0_0;
		C0[1] = res0_1;
		C0[2] = res0_2;
		C0[3] = res0_3;


#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		temp = bk - off;
#ifdef LEFT
		temp -= 4; // number of values in A
#else
		temp -= 1; // number of values in B
#endif
		ptrba += temp*4;
		ptrbb += temp*1;
#endif

#ifdef LEFT
		off += 4; // number of values in A
#endif

		C0 = C0+4;

	}

	if ( bm & 2 ) // do any 2x1 loop
	{

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		ptrbb = bb;
#else
		ptrba += off*2;
		ptrbb = bb + off*1;
#endif

		res0_0 = 0;
		res0_1 = 0;



#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
		temp = bk-off;
#elif defined(LEFT)
		temp = off+2;	// number of values in A
#else
		temp = off+1;	// number of values in B
#endif

		for (k=0; k<temp; k++)
                {
			b0 = ptrbb[0];

			a0 = ptrba[0];
			res0_0 += a0*b0;

			a1 = ptrba[1];
			res0_1 += a1*b0;

			ptrba = ptrba+2;
			ptrbb = ptrbb+1;
                }

		res0_0 *= alpha;
		res0_1 *= alpha;

		C0[0] = res0_0;
		C0[1] = res0_1;


#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		temp = bk - off;
#ifdef LEFT
		temp -= 2; // number of values in A
#else
		temp -= 1; // number of values in B
#endif
		ptrba += temp*2;
		ptrbb += temp*1;
#endif

#ifdef LEFT
		off += 2; // number of values in A
#endif

		C0 = C0+2;

	}

	if ( bm & 1 ) // do any 1x1 loop
	{

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		ptrbb = bb;
#else
		ptrba += off*1;
		ptrbb = bb + off*1;
#endif

		res0_0 = 0;


#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
		temp = bk-off;
#elif defined(LEFT)
		temp = off+1;	// number of values in A
#else
		temp = off+1;	// number of values in B
#endif

		for (k=0; k<temp; k++)
                {
			b0 = ptrbb[0];

			a0 = ptrba[0];
			res0_0 += a0*b0;

			ptrba = ptrba+1;
			ptrbb = ptrbb+1;
                }

		res0_0 *= alpha;

		C0[0] = res0_0;


#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		temp = bk - off;
#ifdef LEFT
		temp -= 1; // number of values in A
#else
		temp -= 1; // number of values in B
#endif
		ptrba += temp*1;
		ptrbb += temp*1;
#endif

#ifdef LEFT
		off += 1; // number of values in A
#endif

		C0 = C0+1;

	}



#if defined(TRMMKERNEL) && !defined(LEFT)
		off += 1;
#endif

        k = (bk<<0);
        bb = bb+k;
        C = C+ldc;
   }
   return 0;
}
