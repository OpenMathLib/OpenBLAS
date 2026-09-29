#include "common.h"
#include <stdbool.h>
#if defined(__wasm_simd128__)
#include <wasm_simd128.h>
#if defined(__wasm_relaxed_simd__)
#define MADD_F32(a, b, c) wasm_f32x4_relaxed_madd((a), (b), (c))
#define MADD_F64(a, b, c) wasm_f64x2_relaxed_madd((a), (b), (c))
#else
#define MADD_F32(a, b, c) wasm_f32x4_add((c), wasm_f32x4_mul((a), (b)))
#define MADD_F64(a, b, c) wasm_f64x2_add((c), wasm_f64x2_mul((a), (b)))
#endif
#endif

int CNAME(BLASLONG bm,BLASLONG bn,BLASLONG bk,FLOAT alpha,FLOAT* ba,FLOAT* bb,FLOAT* C,BLASLONG ldc ,BLASLONG offset)
{

   BLASLONG i,j,k;
   FLOAT *C0,*C1,*C2,*C3,*ptrba,*ptrbb;

   FLOAT res0_0;
   FLOAT res0_1;
   FLOAT res0_2;
   FLOAT res0_3;

   FLOAT res1_0;
   FLOAT res1_1;
   FLOAT res1_2;
   FLOAT res1_3;

   FLOAT res2_0;
   FLOAT res2_1;
   FLOAT res2_2;
   FLOAT res2_3;

   FLOAT res3_0;
   FLOAT res3_1;
   FLOAT res3_2;
   FLOAT res3_3;

   FLOAT a0;
   FLOAT a1;

   FLOAT b0;
   FLOAT b1;
   FLOAT b2;
   FLOAT b3;

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


   for (j=0; j<bn/4; j+=1) // do blocks of the Mx4 loops 
   {
        C0 = C;
        C1 = C0+ldc;
        C2 = C1+ldc;
        C3 = C2+ldc;


        if (left) {
            off = offset;
        }

        ptrba = ba;

        for (i=0; i<bm/4; i+=1) // do blocks of 4x4
	{

		ptrbb = bb;
                if (backwards)
                {
		   ptrba += off*4; // number of values in A
		   ptrbb += off*4; // number of values in B
                }

		res0_0 = 0;
		res0_1 = 0;
		res0_2 = 0;
		res0_3 = 0;

		res1_0 = 0;
		res1_1 = 0;
		res1_2 = 0;
		res1_3 = 0;

		res2_0 = 0;
		res2_1 = 0;
		res2_2 = 0;
		res2_3 = 0;

		res3_0 = 0;
		res3_1 = 0;
		res3_2 = 0;
		res3_3 = 0;

                temp = backwards ? bk-off :
                             left ? off + 4 : // number of values in A
                                    off + 4;  // number of values in B

#if defined(__wasm_simd128__) && !defined(DOUBLE)
		{
		  v128_t acc0 = wasm_f32x4_splat(0.0f);
		  v128_t acc1 = wasm_f32x4_splat(0.0f);
		  v128_t acc2 = wasm_f32x4_splat(0.0f);
		  v128_t acc3 = wasm_f32x4_splat(0.0f);
		  for (k = 0; k < temp; k++) {
		    v128_t va = wasm_v128_load(ptrba);
		    v128_t vb = wasm_v128_load(ptrbb);
		    acc0 = MADD_F32(va, wasm_i32x4_shuffle(vb, vb, 0, 0, 0, 0), acc0);
		    acc1 = MADD_F32(va, wasm_i32x4_shuffle(vb, vb, 1, 1, 1, 1), acc1);
		    acc2 = MADD_F32(va, wasm_i32x4_shuffle(vb, vb, 2, 2, 2, 2), acc2);
		    acc3 = MADD_F32(va, wasm_i32x4_shuffle(vb, vb, 3, 3, 3, 3), acc3);
		    ptrba += 4;
		    ptrbb += 4;
		  }
		  v128_t valpha = wasm_f32x4_splat(alpha);
		  wasm_v128_store(C0, wasm_f32x4_mul(acc0, valpha));
		  wasm_v128_store(C1, wasm_f32x4_mul(acc1, valpha));
		  wasm_v128_store(C2, wasm_f32x4_mul(acc2, valpha));
		  wasm_v128_store(C3, wasm_f32x4_mul(acc3, valpha));
		}
#elif defined(__wasm_simd128__) && defined(DOUBLE)
		{
		  v128_t a0l = wasm_f64x2_splat(0.0), a0h = wasm_f64x2_splat(0.0);
		  v128_t a1l = wasm_f64x2_splat(0.0), a1h = wasm_f64x2_splat(0.0);
		  v128_t a2l = wasm_f64x2_splat(0.0), a2h = wasm_f64x2_splat(0.0);
		  v128_t a3l = wasm_f64x2_splat(0.0), a3h = wasm_f64x2_splat(0.0);
		  for (k = 0; k < temp; k++) {
		    v128_t va01 = wasm_v128_load(ptrba);
		    v128_t va23 = wasm_v128_load(ptrba + 2);
		    v128_t vb01 = wasm_v128_load(ptrbb);
		    v128_t vb23 = wasm_v128_load(ptrbb + 2);
		    v128_t b0v = wasm_i64x2_shuffle(vb01, vb01, 0, 0);
		    v128_t b1v = wasm_i64x2_shuffle(vb01, vb01, 1, 1);
		    v128_t b2v = wasm_i64x2_shuffle(vb23, vb23, 0, 0);
		    v128_t b3v = wasm_i64x2_shuffle(vb23, vb23, 1, 1);
		    a0l = MADD_F64(va01, b0v, a0l);
		    a0h = MADD_F64(va23, b0v, a0h);
		    a1l = MADD_F64(va01, b1v, a1l);
		    a1h = MADD_F64(va23, b1v, a1h);
		    a2l = MADD_F64(va01, b2v, a2l);
		    a2h = MADD_F64(va23, b2v, a2h);
		    a3l = MADD_F64(va01, b3v, a3l);
		    a3h = MADD_F64(va23, b3v, a3h);
		    ptrba += 4;
		    ptrbb += 4;
		  }
		  v128_t valpha = wasm_f64x2_splat(alpha);
		  wasm_v128_store(C0, wasm_f64x2_mul(a0l, valpha));
		  wasm_v128_store(C0 + 2, wasm_f64x2_mul(a0h, valpha));
		  wasm_v128_store(C1, wasm_f64x2_mul(a1l, valpha));
		  wasm_v128_store(C1 + 2, wasm_f64x2_mul(a1h, valpha));
		  wasm_v128_store(C2, wasm_f64x2_mul(a2l, valpha));
		  wasm_v128_store(C2 + 2, wasm_f64x2_mul(a2h, valpha));
		  wasm_v128_store(C3, wasm_f64x2_mul(a3l, valpha));
		  wasm_v128_store(C3 + 2, wasm_f64x2_mul(a3h, valpha));
		}
#else
		for (k=0; k<temp; k++)
                {
			b0 = ptrbb[0];
			b1 = ptrbb[1];
			b2 = ptrbb[2];
			b3 = ptrbb[3];

			a0 = ptrba[0];
			res0_0 += a0*b0;
			res1_0 += a0*b1;
			res2_0 += a0*b2;
			res3_0 += a0*b3;

			a1 = ptrba[1];
			res0_1 += a1*b0;
			res1_1 += a1*b1;
			res2_1 += a1*b2;
			res3_1 += a1*b3;

			a0 = ptrba[2];
			res0_2 += a0*b0;
			res1_2 += a0*b1;
			res2_2 += a0*b2;
			res3_2 += a0*b3;

			a1 = ptrba[3];
			res0_3 += a1*b0;
			res1_3 += a1*b1;
			res2_3 += a1*b2;
			res3_3 += a1*b3;

			ptrba = ptrba+4;
			ptrbb = ptrbb+4;
                }

		res0_0 *= alpha;
		res0_1 *= alpha;
		res0_2 *= alpha;
		res0_3 *= alpha;

		res1_0 *= alpha;
		res1_1 *= alpha;
		res1_2 *= alpha;
		res1_3 *= alpha;

		res2_0 *= alpha;
		res2_1 *= alpha;
		res2_2 *= alpha;
		res2_3 *= alpha;

		res3_0 *= alpha;
		res3_1 *= alpha;
		res3_2 *= alpha;
		res3_3 *= alpha;

		C0[0] = res0_0;
		C0[1] = res0_1;
		C0[2] = res0_2;
		C0[3] = res0_3;

		C1[0] = res1_0;
		C1[1] = res1_1;
		C1[2] = res1_2;
		C1[3] = res1_3;

		C2[0] = res2_0;
		C2[1] = res2_1;
		C2[2] = res2_2;
		C2[3] = res2_3;

		C3[0] = res3_0;
		C3[1] = res3_1;
		C3[2] = res3_2;
		C3[3] = res3_3;
#endif

		if (!backwards) {
                    temp = bk-off;
                    temp = left ? temp - 4 : // number of values in A
                                  temp - 4;  // number of values in B

                    ptrba += temp*4; // number of values in A
		    ptrbb += temp*4; // number of values in B
                }
#ifdef LEFT
		off += 4; // number of values in A
#endif

		C0 = C0+4;
		C1 = C1+4;
		C2 = C2+4;
		C3 = C3+4;

	}

	if ( bm & 2 ) // do any 2x4 loop
	{

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		ptrbb = bb;
#else
		ptrba += off*2;
		ptrbb = bb + off*4;
#endif

		res0_0 = 0;
		res0_1 = 0;

		res1_0 = 0;
		res1_1 = 0;

		res2_0 = 0;
		res2_1 = 0;

		res3_0 = 0;
		res3_1 = 0;


#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
		temp = bk-off;
#elif defined(LEFT)
		temp = off+2;	// number of values in A
#else
		temp = off+4;	// number of values in B
#endif

		for (k=0; k<temp; k++)
                {
			b0 = ptrbb[0];
			b1 = ptrbb[1];
			b2 = ptrbb[2];
			b3 = ptrbb[3];

			a0 = ptrba[0];
			res0_0 += a0*b0;
			res1_0 += a0*b1;
			res2_0 += a0*b2;
			res3_0 += a0*b3;

			a1 = ptrba[1];
			res0_1 += a1*b0;
			res1_1 += a1*b1;
			res2_1 += a1*b2;
			res3_1 += a1*b3;

			ptrba = ptrba+2;
			ptrbb = ptrbb+4;
                }

		res0_0 *= alpha;
		res0_1 *= alpha;

		res1_0 *= alpha;
		res1_1 *= alpha;

		res2_0 *= alpha;
		res2_1 *= alpha;

		res3_0 *= alpha;
		res3_1 *= alpha;

		C0[0] = res0_0;
		C0[1] = res0_1;

		C1[0] = res1_0;
		C1[1] = res1_1;

		C2[0] = res2_0;
		C2[1] = res2_1;

		C3[0] = res3_0;
		C3[1] = res3_1;


#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		temp = bk - off;
#ifdef LEFT
		temp -= 2; // number of values in A
#else
		temp -= 4; // number of values in B
#endif
		ptrba += temp*2;
		ptrbb += temp*4;
#endif

#ifdef LEFT
		off += 2; // number of values in A
#endif

		C0 = C0+2;
		C1 = C1+2;
		C2 = C2+2;
		C3 = C3+2;

	}

	if ( bm & 1 ) // do any 1x4 loop
	{

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		ptrbb = bb;
#else
		ptrba += off*1;
		ptrbb = bb + off*4;
#endif

		res0_0 = 0;
		res1_0 = 0;
		res2_0 = 0;
		res3_0 = 0;


#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
		temp = bk-off;
#elif defined(LEFT)
		temp = off+1;	// number of values in A
#else
		temp = off+4;	// number of values in B
#endif

		for (k=0; k<temp; k++)
                {
			b0 = ptrbb[0];
			b1 = ptrbb[1];
			b2 = ptrbb[2];
			b3 = ptrbb[3];

			a0 = ptrba[0];
			res0_0 += a0*b0;
			res1_0 += a0*b1;
			res2_0 += a0*b2;
			res3_0 += a0*b3;

			ptrba = ptrba+1;
			ptrbb = ptrbb+4;
                }

		res0_0 *= alpha;

		res1_0 *= alpha;

		res2_0 *= alpha;

		res3_0 *= alpha;

		C0[0] = res0_0;

		C1[0] = res1_0;

		C2[0] = res2_0;

		C3[0] = res3_0;


#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		temp = bk - off;
#ifdef LEFT
		temp -= 1; // number of values in A
#else
		temp -= 4; // number of values in B
#endif
		ptrba += temp*1;
		ptrbb += temp*4;
#endif

#ifdef LEFT
		off += 1; // number of values in A
#endif

		C0 = C0+1;
		C1 = C1+1;
		C2 = C2+1;
		C3 = C3+1;

	}


#if defined(TRMMKERNEL) && !defined(LEFT)
		off += 4;
#endif

        k = (bk<<2);
        bb = bb+k;
        i = (ldc<<2);
        C = C+i;
    }

   for (j=0; j<(bn&2); j+=2) // do the Mx2 loops 
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
