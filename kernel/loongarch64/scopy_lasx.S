#define ASSEMBLER

#include "common.h"
#define N      $r4
#define X      $r5
#define INCX   $r6
#define Y      $r7
#define INCY   $r8
#define I      $r17
#define TEMP   $r18
#define t1     $r14
#define t2     $r15
#define t3     $r16
#define t4     $r19
#define a1     $f12
#define a2     $f13
#define a3     $f14
#define a4     $f15
#define VX0    $xr12
#define VX1    $xr13

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
    addi.d  I, I, -1
    xvst VX0, Y, 0 * SIZE
    addi.d X, X, 8 * SIZE
    addi.d Y, Y, 8 * SIZE
    blt $r0, I, .L111
    .align 3

.L112:
    andi I, N, 7
    bge $r0, I, .L999
    .align 3

.L113:
    fld.s $f12, X, 0 * SIZE
    addi.d I, I, -1
    addi.d  X, X, SIZE
    fst.s $f12, Y, 0 * SIZE
    addi.d  Y, Y, SIZE
    blt $r0, I, .L113
    b .L999
    .align 3

.L12:
    bge $r0, I, .L122
    .align 3

.L121:
    xvld VX0, X, 0 * SIZE
    xvstelm.w VX0, Y, 0, 0
    add.d Y, Y, INCY
    xvstelm.w VX0, Y, 0, 1
    add.d Y, Y, INCY
    xvstelm.w VX0, Y, 0, 2
    add.d Y, Y, INCY
    xvstelm.w VX0, Y, 0, 3
    add.d Y, Y, INCY
    xvstelm.w VX0, Y, 0, 4
    add.d Y, Y, INCY
    xvstelm.w VX0, Y, 0, 5
    add.d Y, Y, INCY
    xvstelm.w VX0, Y, 0, 6
    add.d Y, Y, INCY
    xvstelm.w VX0, Y, 0, 7
    add.d Y, Y, INCY
    addi.d X, X, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L121
    .align 3

.L122:
    andi I, N, 7
    bge $r0, I, .L999
    .align 3

.L123:
    fld.s $f12, X, 0 * SIZE
    addi.d I, I, -1
    addi.d  X, X, SIZE
    fst.s $f12, Y, 0 * SIZE
    add.d  Y, Y, INCY
    blt $r0, I, .L123
    b .L999
    .align 3

.L21:
    bge $r0, I, .L212
    .align 3

.L211:
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
    xvst VX0, Y, 0 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L211
    .align 3

.L212:
    andi I, N, 7
    bge $r0, I, .L999
    .align 3

.L213:
    fld.s $f12, X, 0 * SIZE
    addi.d I, I, -1
    fst.s $f12, Y, 0 * SIZE
    add.d  X, X, INCX
    addi.d  Y, Y, SIZE
    blt $r0, I, .L213
    b .L999
    .align 3

.L22:
    bge $r0, I, .L223
    .align 3

.L222:
    fld.s a1, X, 0 * SIZE
    add.d X, X, INCX
    fld.s a2, X, 0 * SIZE
    add.d X, X, INCX
    fld.s a3, X, 0 * SIZE
    add.d X, X, INCX
    fld.s a4, X, 0 * SIZE
    add.d X, X, INCX
    fst.s a1, Y, 0 * SIZE
    add.d Y, Y, INCY
    fst.s a2, Y, 0 * SIZE
    add.d Y, Y, INCY
    fst.s a3, X, 0 * SIZE
    add.d Y, Y, INCY
    fst.s a4, X, 0 * SIZE
    add.d Y, Y, INCY
    fld.s a1, X, 0 * SIZE
    add.d X, X, INCX
    fld.s a2, X, 0 * SIZE
    add.d X, X, INCX
    fld.s a3, X, 0 * SIZE
    add.d X, X, INCX
    fld.s a4, X, 0 * SIZE
    add.d X, X, INCX
    fst.s a1, Y, 0 * SIZE
    add.d Y, Y, INCY
    fst.s a2, Y, 0 * SIZE
    add.d Y, Y, INCY
    fst.s a3, X, 0 * SIZE
    add.d Y, Y, INCY
    fst.s a4, X, 0 * SIZE
    add.d Y, Y, INCY
    addi.d  I, I, -1
    blt $r0, I, .L222
    .align 3

.L223:
    andi I, N, 7
    bge $r0, I, .L999
    .align 3

.L224:
    fld.s $f12, X, 0 * SIZE
    addi.d I, I, -1
    fst.s $f12, Y, 0 * SIZE
    add.d  X, X, INCX
    add.d  Y, Y, INCY
    blt $r0, I, .L224
    .align 3

.L999:
    move $r4, $r12
    jirl $r0, $r1, 0x0
    .align 3

    EPILOGUE
