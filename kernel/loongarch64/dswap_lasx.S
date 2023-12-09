#define ASSEMBLER

#include "common.h"
#define N      $r4
#define X      $r7
#define INCX   $r8
#define Y      $r9
#define INCY   $r10

#define I      $r17
#define TEMP   $r18
#define XX     $r5
#define YY     $r6
#define t1     $r14
#define t2     $r15
#define t3     $r16
#define t4     $r19
#define a1     $f12
#define a2     $f13
#define a3     $f14
#define a4     $f15
#define b1     $f16
#define b2     $f17
#define b3     $f18
#define b4     $f19
#define VX0    $xr12
#define VX1    $xr13
#define VX2    $xr14
#define VX3    $xr15


    PROLOGUE
    bge $r0, N, .L999
    li.d TEMP, 1
    slli.d  TEMP, TEMP, BASE_SHIFT
    slli.d  INCX, INCX, BASE_SHIFT
    slli.d  INCY, INCY, BASE_SHIFT
    srai.d I, N, 3
    bne INCX, TEMP, .L20
    bne INCY, TEMP, .L12 // INCX==1 and INCY!=1
    b .L11  // INCX==1 and INCY==1
.L20:
    bne INCY, TEMP, .L22 // INCX!=1 and INCY!=1
    b .L21 // INCX!=1 and INCY==1

.L11:
    bge $r0, I, .L112
    .align 3

.L111:
    xvld VX0, X, 0 * SIZE
    xvld VX1, X, 4 * SIZE
    xvld VX2, Y, 0 * SIZE
    xvld VX3, Y, 4 * SIZE
    addi.d  I, I, -1
    xvst VX2, X, 0 * SIZE
    xvst VX3, X, 4 * SIZE
    xvst VX0, Y, 0 * SIZE
    xvst VX1, Y, 4 * SIZE
    addi.d X, X, 8 * SIZE
    addi.d Y, Y, 8 * SIZE
    blt $r0, I, .L111
    .align 3

.L112:
    andi I, N, 7
    bge $r0, I, .L999
    .align 3

.L113:
    fld.d $f12, X, 0 * SIZE
    fld.d $f14, Y, 0 * SIZE
    addi.d I, I, -1
    fst.d $f12, Y, 0 * SIZE
    fst.d $f14, X, 0 * SIZE
    addi.d  X, X, SIZE
    addi.d  Y, Y, SIZE
    blt $r0, I, .L113
    b .L999
    .align 3

.L12: // INCX==1 and INCY!=1
    bge $r0, I, .L122
    .align 3

.L121:
    xvld VX0, X, 0 * SIZE
    ld.d t1, Y, 0 * SIZE
    xvstelm.d VX0, Y, 0, 0
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    xvstelm.d VX0, Y, 0, 1
    add.d Y, Y, INCY
    ld.d t3, Y, 0 * SIZE
    xvstelm.d VX0, Y, 0, 2
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    xvstelm.d VX0, Y, 0, 3
    xvinsgr2vr.d VX2, t1, 0
    xvinsgr2vr.d VX2, t2, 1
    xvinsgr2vr.d VX2, t3, 2
    xvinsgr2vr.d VX2, t4, 3
    add.d Y, Y, INCY
    xvst VX2, X, 0 * SIZE
    xvld VX1, X, 4 * SIZE
    ld.d t1, Y, 0 * SIZE
    xvstelm.d VX1, Y, 0, 0
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    xvstelm.d VX1, Y, 0, 1
    add.d Y, Y, INCY
    ld.d t3, Y, 0 * SIZE
    xvstelm.d VX1, Y, 0, 2
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    xvstelm.d VX1, Y, 0, 3
    xvinsgr2vr.d VX3, t1, 0
    xvinsgr2vr.d VX3, t2, 1
    xvinsgr2vr.d VX3, t3, 2
    xvinsgr2vr.d VX3, t4, 3
    add.d Y, Y, INCY
    xvst VX3, X, 4 * SIZE
    addi.d X, X, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L121
    .align 3

.L122:
    andi I, N, 7
    bge $r0, I, .L999
    .align 3

.L123:
    fld.d $f12, X, 0 * SIZE
    fld.d $f14, Y, 0 * SIZE
    addi.d I, I, -1
    fst.d $f12, Y, 0 * SIZE
    fst.d $f14, X, 0 * SIZE
    addi.d  X, X, SIZE
    add.d  Y, Y, INCY
    blt $r0, I, .L123
    b .L999
    .align 3

.L21:
    bge $r0, I, .L212
    .align 3

.L211:
    xvld VX2, Y, 0 * SIZE
    ld.d t1, X, 0 * SIZE
    xvstelm.d VX2, X, 0, 0
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    xvstelm.d VX2, X, 0, 1
    add.d X, X, INCX
    ld.d t3, X, 0 * SIZE
    xvstelm.d VX2, X, 0, 2
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    xvstelm.d VX2, X, 0, 3
    xvinsgr2vr.d VX0, t1, 0
    xvinsgr2vr.d VX0, t2, 1
    xvinsgr2vr.d VX0, t3, 2
    xvinsgr2vr.d VX0, t4, 3
    add.d X, X, INCX
    xvst VX0, Y, 0 * SIZE
    xvld VX3, Y, 4 * SIZE
    ld.d t1, X, 0 * SIZE
    xvstelm.d VX3, X, 0, 0
    add.d X, X, INCY
    ld.d t2, X, 0 * SIZE
    xvstelm.d VX3, X, 0, 1
    add.d X, X, INCX
    ld.d t3, X, 0 * SIZE
    xvstelm.d VX3, X, 0, 2
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    xvstelm.d VX3, X, 0, 3
    xvinsgr2vr.d VX1, t1, 0
    xvinsgr2vr.d VX1, t2, 1
    xvinsgr2vr.d VX1, t3, 2
    xvinsgr2vr.d VX1, t4, 3
    add.d X, X, INCX
    xvst VX1, Y, 0 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L211
    .align 3

.L212:
    andi I, N, 7
    bge $r0, I, .L999
    .align 3

.L213:
    fld.d $f12, X, 0 * SIZE
    fld.d $f14, Y, 0 * SIZE
    addi.d I, I, -1
    fst.d $f12, Y, 0 * SIZE
    fst.d $f14, X, 0 * SIZE
    add.d  X, X, INCX
    addi.d  Y, Y, SIZE
    blt $r0, I, .L213
    b .L999
    .align 3

.L22:
    bgez INCX, .L220
    //addi.d TEMP, N, -1
    //mul.d TEMP, TEMP, INCX
    //sub.d X, X, TEMP
    .align 3

.L220:
    bge $r0, I, .L223
    .align 3
    move XX, X

.L222:
    fld.d a1, X, 0 * SIZE
    add.d X, X, INCX
    fld.d a2, X, 0 * SIZE
    add.d X, X, INCX
    fld.d a3, X, 0 * SIZE
    add.d X, X, INCX
    fld.d a4, X, 0 * SIZE
    add.d X, X, INCX
    fld.d b1, Y, 0 * SIZE
    fst.d a1, Y, 0 * SIZE
    add.d Y, Y, INCY
    fld.d b2, Y, 0 * SIZE
    fst.d a2, Y, 0 * SIZE
    add.d Y, Y, INCY
    fld.d b3, Y, 0 * SIZE
    fst.d a3, Y, 0 * SIZE
    add.d Y, Y, INCY
    fld.d b4, Y, 0 * SIZE
    fst.d a4, Y, 0 * SIZE
    add.d Y, Y, INCY
    fld.d a1, X, 0 * SIZE
    add.d X, X, INCX
    fst.d b1, XX, 0 * SIZE
    add.d XX, XX, INCX
    fld.d b1, Y, 0 * SIZE
    fst.d a1, Y, 0 * SIZE
    add.d Y, Y, INCY
    fld.d a2, X, 0 * SIZE
    add.d X, X, INCX
    fst.d b2, XX, 0 * SIZE
    add.d XX, XX, INCX
    fld.d b2, Y, 0 * SIZE
    fst.d a2, Y, 0 * SIZE
    add.d Y, Y, INCY
    fld.d a3, X, 0 * SIZE
    add.d X, X, INCX
    fst.d b3, XX, 0 * SIZE
    add.d XX, XX, INCX
    fld.d b3, Y, 0 * SIZE
    fst.d a3, Y, 0 * SIZE
    fld.d a4, X, 0 * SIZE
    add.d X, X, INCX
    fst.d b4, XX, 0 * SIZE
    add.d XX, XX, INCX
    fld.d b4, Y, 0 * SIZE
    fst.d a4, Y, 0 * SIZE
    add.d Y, Y, INCY
    fst.d b1, XX, 0 * SIZE
    add.d XX, XX, INCX
    fst.d b2, XX, 0 * SIZE
    add.d XX, XX, INCX
    fst.d b3, XX, 0 * SIZE
    add.d XX, XX, INCX
    fst.d b4, XX, 0 * SIZE
    add.d XX, XX, INCX
    addi.d  I, I, -1
    blt $r0, I, .L222
    .align 3

.L223:
    andi I, N, 7
    bge $r0, I, .L999
    .align 3

.L224:
    fld.d $f12, X, 0 * SIZE
    fld.d $f14, Y, 0 * SIZE
    addi.d I, I, -1
    fst.d $f12, Y, 0 * SIZE
    fst.d $f14, X, 0 * SIZE
    add.d  X, X, INCX
    add.d  Y, Y, INCY
    blt $r0, I, .L224
    .align 3

.L999:
    move $r4, $r12
    jirl $r0, $r1, 0x0
    .align 3

    EPILOGUE
