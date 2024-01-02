/*****************************************************************************
Copyright (c) 2023, The OpenBLAS Project
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.
   3. Neither the name of the OpenBLAS project nor the names of
      its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********************************************************************************/

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

/* Don't change following FR unless you know the effects. */
#define VX0    $xr15
#define VX1    $xr16
#define VX2    $xr17
#define VX3    $xr18
#define VX4    $xr21
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
   xvld VX0, X, 0
   xvfcvtl.d.s VX1, VX0
   xvfcvth.d.s VX2, VX0
   xvfmadd.d res1, VX1, VX1, res1
   xvfmadd.d res2, VX2, VX2, res2
   addi.d I, I, -1
   addi.d X, X, 8 * SIZE
   blt $r0, I, .L10
   .align 3
   b .L996

.L20:
   bge $r0, I, .L997
   .align 3

.L21:
   ld.w t1, X, 0
   add.d X, X, INCX
   ld.w t2, X, 0
   add.d X, X, INCX
   ld.w t3, X, 0
   add.d X, X, INCX
   ld.w t4, X, 0
   add.d X, X, INCX
   xvinsgr2vr.w VX0, t1, 0
   xvinsgr2vr.w VX0, t2, 1
   xvinsgr2vr.w VX0, t3, 2
   xvinsgr2vr.w VX0, t4, 3
   ld.w t1, X, 0
   add.d X, X, INCX
   ld.w t2, X, 0
   add.d X, X, INCX
   ld.w t3, X, 0
   add.d X, X, INCX
   ld.w t4, X, 0
   xvinsgr2vr.w VX0, t1, 4
   xvinsgr2vr.w VX0, t2, 5
   xvinsgr2vr.w VX0, t3, 6
   xvinsgr2vr.w VX0, t4, 7
   add.d X, X, INCX
   xvfcvtl.d.s VX1, VX0
   xvfcvth.d.s VX2, VX0
   xvfmadd.d res1, VX1, VX1, res1
   xvfmadd.d res2, VX2, VX2, res2
   addi.d  I, I, -1
   blt $r0, I, .L21
   b .L996

.L996:
   xvfadd.d res1, res1, res2
   xvpickve.d VX1, res1, 1
   xvpickve.d VX2, res1, 2
   xvpickve.d VX3, res1, 3
   fadd.d   $f19, $f19, $f16
   fadd.d   $f19, $f19, $f17
   fadd.d   $f19, $f19, $f18
   .align 3

.L997:
   andi I, N, 7
   bge $r0, I, .L999
   .align 3

.L998:
   fld.s $f15, X, 0
   add.d    X, X, INCX
   addi.d   I, I, -1
   fcvt.d.s $f15, $f15
   fmadd.d  $f19, $f15, $f15, $f19
   blt $r0, I, .L998
   .align 3

.L999:
   fsqrt.d $f19, $f19
   move $r4, $r17
   fcvt.s.d $f0, $f19
   jirl $r0, $r1, 0x0

   EPILOGUE
