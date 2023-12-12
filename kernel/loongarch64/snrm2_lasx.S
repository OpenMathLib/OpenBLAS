#define ASSEMBLER

#include "common.h"

#define N      $r4
#define X      $r5
#define INCX   $r6
#define I      $r17
#define TEMP   $r18
#define t1     $r12
#define t2     $r13
#define t3     $r14
#define t4     $r15
#define VX0    $xr15
#define VX1    $xr16
#define VX2    $xr17
#define VX3    $xr18
#define res1   $xr19
#define res2   $xr20

   PROLOGUE

#ifdef F_INTERFACE
   LDINT   N,     0(N)
   LDINT   INCX,  0(INCX)
#endif

   xvxor.v res1, res1, res1
   xvxor.v res2, res2, res2
   bge $r0, N, .L999
   beq $r0, INCX, .L999
   li.d  TEMP, SIZE
   slli.d INCX, INCX, BASE_SHIFT
   srai.d I, N, 3
   bne INCX, TEMP, .L20
   bge $r0,  I, .L997
   .align 3

.L10:
   xvld VX0, X, 0 * SIZE
   xvld VX1, X, 0 * SIZE
   xvfcvtl.d.s VX0, VX0
   xvfcvth.d.s VX1, VX1
   xvfmadd.d res1, VX0, VX0, res1
   xvfmadd.d res2, VX1, VX1, res2
   addi.d I, I, -1
   addi.d X, X, 8 * SIZE 
   blt $r0, I, .L10
   .align 3
   b .L996

.L20:
   bge $r0, I, .L997
   .align 3

.L21:
   ld.w t1, X, 0 * SIZE
   add.d X, X, INCX
   ld.w t2, X, 0 * SIZE
   add.d X, X, INCX
   ld.w t3, X, 0 * SIZE
   add.d X, X, INCX
   ld.w t4, X, 0 * SIZE
   add.d X, X, INCX
   xvinsgr2vr.w VX0, t1, 0
   xvinsgr2vr.w VX0, t2, 1
   xvinsgr2vr.w VX0, t3, 2
   xvinsgr2vr.w VX0, t4, 3
   ld.w t1, X, 0 * SIZE
   add.d X, X, INCX
   ld.w t2, X, 0 * SIZE
   add.d X, X, INCX
   ld.w t3, X, 0 * SIZE
   add.d X, X, INCX
   ld.w t4, X, 0 * SIZE
   add.d X, X, INCX
   xvinsgr2vr.w VX0, t1, 4
   xvinsgr2vr.w VX0, t2, 5
   xvinsgr2vr.w VX0, t3, 6
   xvinsgr2vr.w VX0, t4, 7
   ld.w t1, X, 0 * SIZE
   add.d X, X, INCX
   ld.w t2, X, 0 * SIZE
   add.d X, X, INCX
   ld.w t3, X, 0 * SIZE
   add.d X, X, INCX
   ld.w t4, X, 0 * SIZE
   add.d X, X, INCX
   xvinsgr2vr.w VX1, t1, 0
   xvinsgr2vr.w VX1, t2, 1
   xvinsgr2vr.w VX1, t3, 2
   xvinsgr2vr.w VX1, t4, 3
   ld.w t1, X, 0 * SIZE
   add.d X, X, INCX
   ld.w t2, X, 0 * SIZE
   add.d X, X, INCX
   ld.w t3, X, 0 * SIZE
   add.d X, X, INCX
   ld.w t4, X, 0 * SIZE
   add.d X, X, INCX
   xvinsgr2vr.w VX1, t1, 4
   xvinsgr2vr.w VX1, t2, 5
   xvinsgr2vr.w VX1, t3, 6
   xvinsgr2vr.w VX1, t4, 7
   xvfcvtl.d.s VX0, VX0
   xvfcvth.d.s VX1, VX1
   xvfmadd.d res1, VX0, VX0, res1
   xvfmadd.d res2, VX1, VX1, res2
   addi.d  I, I, -1
   blt $r0, I, .L21
   b .L996

.L996:
   xvfadd.d res1, res1, res2
   xvpickve.w VX1, res1, 1
   xvpickve.w VX2, res1, 2
   xvpickve.w VX3, res1, 3
   xvfadd.s res1, VX1, res1
   xvfadd.s res1, VX2, res1
   xvfadd.s res1, VX3, res1
   .align 3

.L997:
   andi I, N, 7
   bge $r0, I, .L999
   .align 3

.L998:
   fld.s $f15, X, 0 * SIZE
   addi.d I, I, -1
   fcvt.d.s $f15, $f15
   fmadd.d $f19, $f15, $f15, $f19
   add.d X, X, INCX
   blt $r0, I, .L998
   .align 3

.L999:
   fsqrt.d $f19, $f19
   move $r4, $r17
   fcvt.s.d $f0, $f19
   jirl $r0, $r1, 0x0

   EPILOGUE
