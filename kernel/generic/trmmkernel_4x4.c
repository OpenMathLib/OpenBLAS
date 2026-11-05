#include "common.h"

/* 4x4 TRMM micro-kernel. Structure mirrors trmmkernel_2x2.c (proven) expanded
 * to the 4/2/1 tiling of gemmkernel_4x4.c so it matches ncopy_4/tcopy_4 and
 * SGEMM/DGEMM UNROLL 4. TRMM overwrites C (no accumulation into existing C). */
int CNAME(BLASLONG bm,BLASLONG bn,BLASLONG bk,FLOAT alpha,FLOAT* ba,FLOAT* bb,FLOAT* C,BLASLONG ldc
#ifdef TRMMKERNEL
		,BLASLONG offset
#endif
		)
{
   BLASLONG i,j,k;
   FLOAT *C0,*C1,*C2,*C3,*ptrba,*ptrbb;
   FLOAT r0c0,r1c0,r2c0,r3c0,r0c1,r1c1,r2c1,r3c1,r0c2,r1c2,r2c2,r3c2,r0c3,r1c3,r2c3,r3c3;
   FLOAT a0,a1,a2,a3,b0,b1,b2,b3;
   BLASLONG off, temp;
#if defined(TRMMKERNEL) && !defined(LEFT)
   off = -offset;
#else
   off = 0;
#endif

#define BKWD ((defined(LEFT)&&defined(TRANSA))||(!defined(LEFT)&&!defined(TRANSA)))

   /* ==================== N panels of 4 ==================== */
   for (j=0; j<bn/4; j+=1)
     {
        C0=C; C1=C0+ldc; C2=C1+ldc; C3=C2+ldc;
#if defined(TRMMKERNEL) && defined(LEFT)
	off = offset;
#endif
        ptrba = ba;
        for (i=0; i<bm/4; i+=1) /* 4x4, mblk=4 nblk=4 */
          {
#if BKWD
             ptrbb = bb;
#else
	     ptrba += off*4;
	     ptrbb = bb + off*4;
#endif
             r0c0=r1c0=r2c0=r3c0=0; r0c1=r1c1=r2c1=r3c1=0; r0c2=r1c2=r2c2=r3c2=0; r0c3=r1c3=r2c3=r3c3=0;
#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	     temp = bk-off;
#elif defined(LEFT)
	     temp = off+4;
#else
	     temp = off+4;
#endif
             for (k=0; k<temp; k+=1)
               {
                  b0=ptrbb[0]; b1=ptrbb[1]; b2=ptrbb[2]; b3=ptrbb[3];
                  a0=ptrba[0]; a1=ptrba[1]; a2=ptrba[2]; a3=ptrba[3];
                  r0c0+=a0*b0; r1c0+=a1*b0; r2c0+=a2*b0; r3c0+=a3*b0;
                  r0c1+=a0*b1; r1c1+=a1*b1; r2c1+=a2*b1; r3c1+=a3*b1;
                  r0c2+=a0*b2; r1c2+=a1*b2; r2c2+=a2*b2; r3c2+=a3*b2;
                  r0c3+=a0*b3; r1c3+=a1*b3; r2c3+=a2*b3; r3c3+=a3*b3;
                  ptrba+=4; ptrbb+=4;
               }
             C0[0]=r0c0*alpha; C0[1]=r1c0*alpha; C0[2]=r2c0*alpha; C0[3]=r3c0*alpha;
             C1[0]=r0c1*alpha; C1[1]=r1c1*alpha; C1[2]=r2c1*alpha; C1[3]=r3c1*alpha;
             C2[0]=r0c2*alpha; C2[1]=r1c2*alpha; C2[2]=r2c2*alpha; C2[3]=r3c2*alpha;
             C3[0]=r0c3*alpha; C3[1]=r1c3*alpha; C3[2]=r2c3*alpha; C3[3]=r3c3*alpha;
#if BKWD
	     temp = bk - off;
#ifdef LEFT
	     temp -= 4;
#else
	     temp -= 4;
#endif
	     ptrba += temp*4;
	     ptrbb += temp*4;
#endif
#ifdef LEFT
	     off += 4;
#endif
             C0+=4; C1+=4; C2+=4; C3+=4;
          }
        if (bm & 2) /* 2x4, mblk=2 nblk=4 */
          {
#if BKWD
             ptrbb = bb;
#else
	     ptrba += off*2;
	     ptrbb = bb + off*4;
#endif
             r0c0=r1c0=0; r0c1=r1c1=0; r0c2=r1c2=0; r0c3=r1c3=0;
#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	     temp = bk-off;
#elif defined(LEFT)
	     temp = off+2;
#else
	     temp = off+4;
#endif
             for (k=0; k<temp; k+=1)
               {
                  b0=ptrbb[0]; b1=ptrbb[1]; b2=ptrbb[2]; b3=ptrbb[3];
                  a0=ptrba[0]; a1=ptrba[1];
                  r0c0+=a0*b0; r1c0+=a1*b0;
                  r0c1+=a0*b1; r1c1+=a1*b1;
                  r0c2+=a0*b2; r1c2+=a1*b2;
                  r0c3+=a0*b3; r1c3+=a1*b3;
                  ptrba+=2; ptrbb+=4;
               }
             C0[0]=r0c0*alpha; C0[1]=r1c0*alpha;
             C1[0]=r0c1*alpha; C1[1]=r1c1*alpha;
             C2[0]=r0c2*alpha; C2[1]=r1c2*alpha;
             C3[0]=r0c3*alpha; C3[1]=r1c3*alpha;
#if BKWD
	     temp = bk-off;
#ifdef LEFT
	     temp -= 2;
#else
	     temp -= 4;
#endif
	     ptrba += temp*2;
	     ptrbb += temp*4;
#endif
#ifdef LEFT
	     off += 2;
#endif
             C0+=2; C1+=2; C2+=2; C3+=2;
          }
        if (bm & 1) /* 1x4, mblk=1 nblk=4 */
          {
#if BKWD
             ptrbb = bb;
#else
	     ptrba += off*1;
	     ptrbb = bb + off*4;
#endif
             r0c0=0; r0c1=0; r0c2=0; r0c3=0;
#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	     temp = bk-off;
#elif defined(LEFT)
	     temp = off+1;
#else
	     temp = off+4;
#endif
             for (k=0; k<temp; k+=1)
               {
                  b0=ptrbb[0]; b1=ptrbb[1]; b2=ptrbb[2]; b3=ptrbb[3]; a0=ptrba[0];
                  r0c0+=a0*b0; r0c1+=a0*b1; r0c2+=a0*b2; r0c3+=a0*b3;
                  ptrba+=1; ptrbb+=4;
               }
             C0[0]=r0c0*alpha; C1[0]=r0c1*alpha; C2[0]=r0c2*alpha; C3[0]=r0c3*alpha;
#if BKWD
	     temp = bk-off;
#ifdef LEFT
	     temp -= 1;
#else
	     temp -= 4;
#endif
	     ptrba += temp*1;
	     ptrbb += temp*4;
#endif
#ifdef LEFT
	     off += 1;
#endif
             C0+=1; C1+=1; C2+=1; C3+=1;
          }
#if defined(TRMMKERNEL) && !defined(LEFT)
	off += 4;
#endif
        bb = bb + bk*4;
        C  = C  + ldc*4;
     }

   /* ==================== N panel of 2 ==================== */
   if (bn & 2)
     {
        C0=C; C1=C0+ldc;
#if defined(TRMMKERNEL) && defined(LEFT)
	off = offset;
#endif
        ptrba = ba;
        for (i=0; i<bm/4; i+=1) /* 4x2, mblk=4 nblk=2 */
          {
#if BKWD
             ptrbb = bb;
#else
	     ptrba += off*4;
	     ptrbb = bb + off*2;
#endif
             r0c0=r1c0=r2c0=r3c0=0; r0c1=r1c1=r2c1=r3c1=0;
#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	     temp = bk-off;
#elif defined(LEFT)
	     temp = off+4;
#else
	     temp = off+2;
#endif
             for (k=0; k<temp; k+=1)
               {
                  b0=ptrbb[0]; b1=ptrbb[1]; a0=ptrba[0]; a1=ptrba[1]; a2=ptrba[2]; a3=ptrba[3];
                  r0c0+=a0*b0; r1c0+=a1*b0; r2c0+=a2*b0; r3c0+=a3*b0;
                  r0c1+=a0*b1; r1c1+=a1*b1; r2c1+=a2*b1; r3c1+=a3*b1;
                  ptrba+=4; ptrbb+=2;
               }
             C0[0]=r0c0*alpha; C0[1]=r1c0*alpha; C0[2]=r2c0*alpha; C0[3]=r3c0*alpha;
             C1[0]=r0c1*alpha; C1[1]=r1c1*alpha; C1[2]=r2c1*alpha; C1[3]=r3c1*alpha;
#if BKWD
	     temp = bk-off;
#ifdef LEFT
	     temp -= 4;
#else
	     temp -= 2;
#endif
	     ptrba += temp*4;
	     ptrbb += temp*2;
#endif
#ifdef LEFT
	     off += 4;
#endif
             C0+=4; C1+=4;
          }
        if (bm & 2) /* 2x2 */
          {
#if BKWD
             ptrbb = bb;
#else
	     ptrba += off*2;
	     ptrbb = bb + off*2;
#endif
             r0c0=r1c0=0; r0c1=r1c1=0;
#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	     temp = bk-off;
#elif defined(LEFT)
	     temp = off+2;
#else
	     temp = off+2;
#endif
             for (k=0; k<temp; k+=1)
               {
                  b0=ptrbb[0]; b1=ptrbb[1]; a0=ptrba[0]; a1=ptrba[1];
                  r0c0+=a0*b0; r1c0+=a1*b0; r0c1+=a0*b1; r1c1+=a1*b1;
                  ptrba+=2; ptrbb+=2;
               }
             C0[0]=r0c0*alpha; C0[1]=r1c0*alpha; C1[0]=r0c1*alpha; C1[1]=r1c1*alpha;
#if BKWD
	     temp = bk-off;
#ifdef LEFT
	     temp -= 2;
#else
	     temp -= 2;
#endif
	     ptrba += temp*2;
	     ptrbb += temp*2;
#endif
#ifdef LEFT
	     off += 2;
#endif
             C0+=2; C1+=2;
          }
        if (bm & 1) /* 1x2, mblk=1 nblk=2 */
          {
#if BKWD
             ptrbb = bb;
#else
	     ptrba += off*1;
	     ptrbb = bb + off*2;
#endif
             r0c0=0; r0c1=0;
#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	     temp = bk-off;
#elif defined(LEFT)
	     temp = off+1;
#else
	     temp = off+2;
#endif
             for (k=0; k<temp; k+=1)
               {
                  b0=ptrbb[0]; b1=ptrbb[1]; a0=ptrba[0];
                  r0c0+=a0*b0; r0c1+=a0*b1; ptrba+=1; ptrbb+=2;
               }
             C0[0]=r0c0*alpha; C1[0]=r0c1*alpha;
#if BKWD
	     temp = bk-off;
#ifdef LEFT
	     temp -= 1;
#else
	     temp -= 2;
#endif
	     ptrba += temp*1;
	     ptrbb += temp*2;
#endif
#ifdef LEFT
	     off += 1;
#endif
             C0+=1; C1+=1;
          }
#if defined(TRMMKERNEL) && !defined(LEFT)
	off += 2;
#endif
        bb = bb + bk*2;
        C  = C  + ldc*2;
     }

   /* ==================== N panel of 1 ==================== */

   if (bn & 1)
     {
        C0=C;
#if defined(TRMMKERNEL) && defined(LEFT)
	off = offset;
#endif
        ptrba = ba;
        for (i=0; i<bm/4; i+=1) /* 4x1, mblk=4 nblk=1 */
          {
#if BKWD
             ptrbb = bb;
#else
	     ptrba += off*4;
	     ptrbb = bb + off*1;
#endif
             r0c0=r1c0=r2c0=r3c0=0;
#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	     temp = bk-off;
#elif defined(LEFT)
	     temp = off+4;
#else
	     temp = off+1;
#endif
             for (k=0; k<temp; k+=1)
               {
                  b0=ptrbb[0]; a0=ptrba[0]; a1=ptrba[1]; a2=ptrba[2]; a3=ptrba[3];
                  r0c0+=a0*b0; r1c0+=a1*b0; r2c0+=a2*b0; r3c0+=a3*b0;
                  ptrba+=4; ptrbb+=1;
               }
             C0[0]=r0c0*alpha; C0[1]=r1c0*alpha; C0[2]=r2c0*alpha; C0[3]=r3c0*alpha;
#if BKWD
	     temp = bk-off;
#ifdef LEFT
	     temp -= 4;
#else
	     temp -= 1;
#endif
	     ptrba += temp*4;
	     ptrbb += temp*1;
#endif
#ifdef LEFT
	     off += 4;
#endif
             C0+=4;
          }
        if (bm & 2) /* 2x1 */
          {
#if BKWD
             ptrbb = bb;
#else
	     ptrba += off*2;
	     ptrbb = bb + off*1;
#endif
             r0c0=r1c0=0;
#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	     temp = bk-off;
#elif defined(LEFT)
	     temp = off+2;
#else
	     temp = off+1;
#endif
             for (k=0; k<temp; k+=1)
               {
                  b0=ptrbb[0]; a0=ptrba[0]; a1=ptrba[1];
                  r0c0+=a0*b0; r1c0+=a1*b0; ptrba+=2; ptrbb+=1;
               }
             C0[0]=r0c0*alpha; C0[1]=r1c0*alpha;
#if BKWD
	     temp = bk-off;
#ifdef LEFT
	     temp -= 2;
#else
	     temp -= 1;
#endif
	     ptrba += temp*2;
	     ptrbb += temp*1;
#endif
#ifdef LEFT
	     off += 2;
#endif
             C0+=2;
          }
        if (bm & 1) /* 1x1 */
          {
#if BKWD
             ptrbb = bb;
#else
	     ptrba += off*1;
	     ptrbb = bb + off*1;
#endif
             r0c0=0;
#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	     temp = bk-off;
#elif defined(LEFT)
	     temp = off+1;
#else
	     temp = off+1;
#endif
             for (k=0; k<temp; k+=1)
               {
                  r0c0+=ptrba[0]*ptrbb[0]; ptrba+=1; ptrbb+=1;
               }
             C0[0]=r0c0*alpha;
#if BKWD
	     temp = bk-off;
#ifdef LEFT
	     temp -= 1;
#else
	     temp -= 1;
#endif
	     ptrba += temp*1;
	     ptrbb += temp*1;
#endif
#ifdef LEFT
	     off += 1;
#endif
             C0+=1;
          }
#if defined(TRMMKERNEL) && !defined(LEFT)
	off += 1;
#endif
        bb = bb + bk*1;
        C  = C  + ldc;
     }

#undef BKWD
   return 0;
}
