#define ASSEMBLER

#include "common.h"
#define N      $r4
#define X      $r5
#define INCX   $r6
#define Y      $r7
#define INCY   $r8
#define C      $f0
#define S      $f1

#define I      $r12
#define TEMP   $r13
#define t1     $r14
#define t2     $r16
#define t3     $r15
#define t4     $r17
#define XX     $r18
#define YY     $r19
#define a1     $f12
#define VX0    $xr8
#define VX1    $xr20
#define VX2    $xr21
#define VX3    $xr22
#define VT0    $xr10
#define VT1    $xr18
#define VXC    $xr23
#define VXS    $xr9
#define VXZ    $xr19

    PROLOGUE

    bge $r0, N, .L999
    li.d TEMP, 1
    movgr2fr.d a1, $r0
    ffint.d.l a1, a1
    slli.d  TEMP, TEMP, BASE_SHIFT
    slli.d  INCX, INCX, BASE_SHIFT
    slli.d  INCY, INCY, BASE_SHIFT
    movfr2gr.d t1, C
    xvreplgr2vr.d VXC, t1
    movfr2gr.d t2, S
    xvreplgr2vr.d VXS, t2
    movfr2gr.d t3, a1
    xvreplgr2vr.d VXZ, t3
    srai.d I, N, 3
    bne INCX, TEMP, .L20
    bne INCY, TEMP, .L12 // INCX==1 and INCY!=1
    b .L11  // INCX==1 and INCY==1
.L20:
    bne INCY, TEMP, .L22 // INCX!=1 and INCY!=1
    b .L21 // INCX!=1 and INCY==1

.L11:
    bge $r0, I, .L997
    fcmp.ceq.d $fcc0, C, a1
    bcnez $fcc0, .L110
    fcmp.ceq.d $fcc0, S, a1
    bcnez $fcc0, .L112 // C!=0 S==0
    b .L111 // C!=0 S!=0
    .align 3

.L110:
    fcmp.ceq.d $fcc0, S, a1
    bcnez $fcc0, .L114 // C==0 S==0
    b .L113 // C==0 S!=0
    .align 3

.L111: // C!=0 S!=0
    xvld VX0, X, 0 * SIZE
    xvld VX2, Y, 0 * SIZE
    xvld VX1, X, 4 * SIZE
    xvld VX3, Y, 4 * SIZE
    xvfmul.d VT0, VX0, VXC
    xvfmadd.d VT0, VX2, VXS, VT0
    xvfmul.d VT1, VX0, VXS
    xvfmsub.d VT1, VX2, VXC, VT1
    xvst VT0, X, 0 * SIZE
    xvst VT1, Y, 0 * SIZE
    xvfmul.d VT0, VX1, VXC
    xvfmadd.d VT0, VX3, VXS, VT0
    xvfmul.d VT1, VX1, VXS
    xvfmsub.d VT1, VX3, VXC, VT1
    xvst VT0, X, 4 * SIZE
    xvst VT1, Y, 4 * SIZE
    addi.d X, X, 8 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L111
    b .L997
    .align 3

.L112: // C!=0 S==0
    xvld VX0, X, 0 * SIZE
    xvld VX2, Y, 0 * SIZE
    xvld VX1, X, 4 * SIZE
    xvld VX3, Y, 4 * SIZE
    xvfmul.d VT0, VX0, VXC
    xvfmul.d VT1, VX2, VXC
    xvst VT0, X, 0 * SIZE
    xvst VT1, Y, 0 * SIZE
    xvfmul.d VT0, VX1, VXC
    xvfmul.d VT1, VX3, VXC
    xvst VT0, X, 4 * SIZE
    xvst VT1, Y, 4 * SIZE
    addi.d X, X, 8 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L112
    b .L997
    .align 3

.L113: // C==0 S!=0
    xvld VX0, X, 0 * SIZE
    xvld VX2, Y, 0 * SIZE
    xvld VX1, X, 4 * SIZE
    xvld VX3, Y, 4 * SIZE
    xvfmul.d VT0, VX2, VXS
    xvfmul.d VT1, VX0, VXS
    xvfsub.d VT1, VXZ, VT1
    xvst VT0, X, 0 * SIZE
    xvst VT1, Y, 0 * SIZE
    xvfmul.d VT0, VX3, VXS
    xvfmul.d VT1, VX1, VXS
    xvfsub.d VT1, VXZ, VT1
    xvst VT0, X, 4 * SIZE
    xvst VT1, Y, 4 * SIZE
    addi.d X, X, 8 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L113
    b .L997
    .align 3

.L114: // C==0 S==0
    xvst VXZ, X, 0 * SIZE
    xvst VXZ, Y, 0 * SIZE
    xvst VXZ, X, 4 * SIZE
    xvst VXZ, Y, 4 * SIZE
    addi.d X, X, 8 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L114
    b .L997
    .align 3

.L12: // INCX==1 and INCY!=1
    bge $r0, I, .L997
    move YY, Y
    move XX, X
    fcmp.ceq.d $fcc0, C, a1
    bcnez $fcc0, .L120
    fcmp.ceq.d $fcc0, S, a1
    bcnez $fcc0, .L122 // C!=0 S==0
    b .L121 // C!=0 S!=0
    .align 3

.L120:
    fcmp.ceq.d $fcc0, S, a1
    bcnez $fcc0, .L124 // C==0 S==0
    b .L123 // C==0 S!=0
    .align 3

.L121: // C!=0 S!=0
    xvld VX0, X, 0 * SIZE
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    xvinsgr2vr.d VX2, t1, 0
    xvinsgr2vr.d VX2, t2, 1
    xvinsgr2vr.d VX2, t3, 2
    xvinsgr2vr.d VX2, t4, 3
    add.d Y, Y, INCY
    xvfmul.d VT0, VX0, VXC
    xvfmadd.d VT0, VX2, VXS, VT0
    xvfmul.d VT1, VX0, VXS
    xvfmsub.d VT1, VX2, VXC, VT1
    xvld VX1, X, 4 * SIZE
    xvst VT0, X, 0 * SIZE
    xvstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 3
    add.d YY, YY, INCY
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    xvinsgr2vr.d VX3, t1, 0
    xvinsgr2vr.d VX3, t2, 1
    xvinsgr2vr.d VX3, t3, 2
    xvinsgr2vr.d VX3, t4, 3
    add.d Y, Y, INCY
    xvfmul.d VT0, VX1, VXC
    xvfmadd.d VT0, VX3, VXS, VT0
    xvfmul.d VT1, VX1, VXS
    xvfmsub.d VT1, VX3, VXC, VT1
    addi.d  I, I, -1
    xvst VT0, X, 4 * SIZE
    xvstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 3
    add.d YY, YY, INCY
    addi.d X, X, 8 * SIZE
    blt $r0, I, .L121
    b .L997
    .align 3

.L122: // C!=0 S==0
    xvld VX0, X, 0 * SIZE
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    xvinsgr2vr.d VX2, t1, 0
    xvinsgr2vr.d VX2, t2, 1
    xvinsgr2vr.d VX2, t3, 2
    xvinsgr2vr.d VX2, t4, 3
    add.d Y, Y, INCY
    xvfmul.d VT0, VX0, VXC
    xvfmul.d VT1, VX2, VXC
    xvld VX1, X, 4 * SIZE
    xvst VT0, X, 0 * SIZE
    xvstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 3
    add.d YY, YY, INCY
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    xvinsgr2vr.d VX3, t1, 0
    xvinsgr2vr.d VX3, t2, 1
    xvinsgr2vr.d VX3, t3, 2
    xvinsgr2vr.d VX3, t4, 3
    add.d Y, Y, INCY
    xvfmul.d VT0, VX1, VXC
    xvfmul.d VT1, VX3, VXC
    addi.d  I, I, -1
    xvst VT0, X, 4 * SIZE
    xvstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 3
    add.d YY, YY, INCY
    addi.d X, X, 8 * SIZE
    blt $r0, I, .L122
    b .L997
    .align 3

.L123: // C==0 S!=0
    xvld VX0, X, 0 * SIZE
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    xvinsgr2vr.d VX2, t1, 0
    xvinsgr2vr.d VX2, t2, 1
    xvinsgr2vr.d VX2, t3, 2
    xvinsgr2vr.d VX2, t4, 3
    add.d Y, Y, INCY
    xvfmul.d VT0, VX2, VXS
    xvfmul.d VT1, VX0, VXS
    xvfsub.d VT1, VXZ, VT1
    xvst VT0, X, 0 * SIZE
    xvstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 3
    add.d YY, YY, INCY
    xvld VX1, X, 4 * SIZE
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    xvinsgr2vr.d VX3, t1, 0
    xvinsgr2vr.d VX3, t2, 1
    xvinsgr2vr.d VX3, t3, 2
    xvinsgr2vr.d VX3, t4, 3
    add.d Y, Y, INCY
    xvfmul.d VT0, VX3, VXS
    xvfmul.d VT1, VX1, VXS
    xvfsub.d VT1, VXZ, VT1
    addi.d  I, I, -1
    xvst VT0, X, 4 * SIZE
    xvstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 3
    add.d YY, YY, INCY
    addi.d X, X, 8 * SIZE
    blt $r0, I, .L123
    b .L997
    .align 3

.L124: // C==0 S==0
    xvst VXZ, X, 0 * SIZE
    xvst VXZ, X, 4 * SIZE
    xvstelm.d VXZ, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VXZ, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VXZ, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VXZ, YY, 0, 3
    add.d YY, YY, INCY
    xvstelm.d VXZ, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VXZ, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VXZ, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VXZ, YY, 0, 3
    add.d YY, YY, INCY
    addi.d  I, I, -1
    blt $r0, I, .L124
    b .L997
    .align 3

.L21:// INCX!=1 and INCY==1
    bge $r0, I, .L997
    move XX, X
    fcmp.ceq.d $fcc0, C, a1
    bcnez $fcc0, .L210
    fcmp.ceq.d $fcc0, S, a1
    bcnez $fcc0, .L212 // C!=0 S==0
    b .L211 // C!=0 S!=0
    .align 3

.L210:
    fcmp.ceq.d $fcc0, S, a1
    bcnez $fcc0, .L214 // C==0 S==0
    b .L213 // C==0 S!=0
    .align 3

.L211: // C!=0 S!=0
    xvld VX2, Y, 0 * SIZE
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    xvinsgr2vr.d VX0, t1, 0
    xvinsgr2vr.d VX0, t2, 1
    xvinsgr2vr.d VX0, t3, 2
    xvinsgr2vr.d VX0, t4, 3
    add.d X, X, INCX
    xvfmul.d VT0, VXC, VX0
    xvfmadd.d VT0, VX2, VXS, VT0
    xvfmul.d VT1, VXS, VX0
    xvfmsub.d VT1, VX2, VXC, VT1
    xvstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 3
    add.d XX, XX, INCX
    xvst VT1, Y, 0 * SIZE
    xvld VX3, Y, 4 * SIZE
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    xvinsgr2vr.d VX1, t1, 0
    xvinsgr2vr.d VX1, t2, 1
    xvinsgr2vr.d VX1, t3, 2
    xvinsgr2vr.d VX1, t4, 3
    add.d X, X, INCX
    xvfmul.d VT0, VX1, VXC
    xvfmadd.d VT0, VX3, VXS, VT0
    xvfmul.d VT1, VX1, VXS
    xvfmsub.d VT1, VX3, VXC, VT1
    xvstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 3
    add.d XX, XX, INCX
    xvst VT1, Y, 4 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L211
    b .L997
    .align 3

.L212: // C!=0 S==0
    xvld VX2, Y, 0 * SIZE
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    xvinsgr2vr.d VX0, t1, 0
    xvinsgr2vr.d VX0, t2, 1
    xvinsgr2vr.d VX0, t3, 2
    xvinsgr2vr.d VX0, t4, 3
    add.d X, X, INCX
    xvfmul.d VT0, VXC, VX0
    xvstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 3
    add.d XX, XX, INCX
    xvfmul.d VT1, VX2, VXC
    xvst VT1, Y, 0 * SIZE
    xvld VX3, Y, 4 * SIZE
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    xvinsgr2vr.d VX1, t1, 0
    xvinsgr2vr.d VX1, t2, 1
    xvinsgr2vr.d VX1, t3, 2
    xvinsgr2vr.d VX1, t4, 3
    add.d X, X, INCX
    xvfmul.d VT0, VX1, VXC
    xvstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 3
    add.d XX, XX, INCX
    xvfmul.d VT1, VX3, VXS
    xvst VT1, Y, 4 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L212
    b .L997
    .align 3

.L213: // C==0 S!=0
    xvld VX2, Y, 0 * SIZE
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    xvinsgr2vr.d VX0, t1, 0
    xvinsgr2vr.d VX0, t2, 1
    xvinsgr2vr.d VX0, t3, 2
    xvinsgr2vr.d VX0, t4, 3
    add.d X, X, INCX
    xvfmul.d VT0, VXS, VX2
    xvfmul.d VT1, VXS, VX0
    xvfsub.d VT1, VXZ, VT1
    xvstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 3
    add.d XX, XX, INCX
    xvst VT1, Y, 0 * SIZE
    xvld VX3, Y, 4 * SIZE
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    xvinsgr2vr.d VX1, t1, 0
    xvinsgr2vr.d VX1, t2, 1
    xvinsgr2vr.d VX1, t3, 2
    xvinsgr2vr.d VX1, t4, 3
    add.d X, X, INCX
    xvfmul.d VT0, VX3, VXS
    xvfmul.d VT1, VX1, VXS
    xvfsub.d VT1, VXZ, VT1
    xvstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 3
    add.d XX, XX, INCX
    xvst VT1, Y, 4 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L213
    b .L997
    .align 3

.L214: // C==0 S==0
    xvstelm.d VXZ, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.d VXZ, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.d VXZ, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.d VXZ, XX, 0, 3
    add.d XX, XX, INCX
    xvst VT1, Y, 0 * SIZE
    xvstelm.d VXZ, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.d VXZ, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.d VXZ, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.d VXZ, XX, 0, 3
    add.d XX, XX, INCX
    xvst VT1, Y, 4 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L211
    b .L997
    .align 3

.L22:
    bge $r0, I, .L997
    move YY, Y
    move XX, X
    fcmp.ceq.d $fcc0, C, a1
    bcnez $fcc0, .L220
    fcmp.ceq.d $fcc0, S, a1
    bcnez $fcc0, .L222 // C!=0 S==0
    b .L221 // C!=0 S!=0
    .align 3

.L220:
    fcmp.ceq.d $fcc0, S, a1
    bcnez $fcc0, .L224 // C==0 S==0
    b .L223 // C==0 S!=0
    .align 3

.L221: // C!=0 S!=0
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    add.d X, X, INCX
    xvinsgr2vr.d VX0, t1, 0
    xvinsgr2vr.d VX0, t2, 1
    xvinsgr2vr.d VX0, t3, 2
    xvinsgr2vr.d VX0, t4, 3
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    xvinsgr2vr.d VX2, t1, 0
    xvinsgr2vr.d VX2, t2, 1
    xvinsgr2vr.d VX2, t3, 2
    xvinsgr2vr.d VX2, t4, 3
    add.d Y, Y, INCY
    xvfmul.d VT0, VX0, VXC
    xvfmadd.d VT0, VX2, VXS, VT0
    xvfmul.d VT1, VX0, VXS
    xvfmsub.d VT1, VX2, VXC, VT1
    xvstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 3
    add.d XX, XX, INCX
    xvstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 3
    add.d YY, YY, INCY
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    add.d X, X, INCX
    xvinsgr2vr.d VX1, t1, 0
    xvinsgr2vr.d VX1, t2, 1
    xvinsgr2vr.d VX1, t3, 2
    xvinsgr2vr.d VX1, t4, 3
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    xvinsgr2vr.d VX3, t1, 0
    xvinsgr2vr.d VX3, t2, 1
    xvinsgr2vr.d VX3, t3, 2
    xvinsgr2vr.d VX3, t4, 3
    add.d Y, Y, INCY
    xvfmul.d VT0, VX1, VXC
    xvfmadd.d VT0, VX3, VXS, VT0
    xvfmul.d VT1, VX0, VXS
    xvfmsub.d VT1, VX3, VXC, VT1
    xvstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 3
    add.d XX, XX, INCX
    xvstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 3
    add.d YY, YY, INCY
    addi.d  I, I, -1
    blt $r0, I, .L221
    b .L997
    .align 3

.L222: // C!=0 S==0
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    add.d X, X, INCX
    xvinsgr2vr.d VX0, t1, 0
    xvinsgr2vr.d VX0, t2, 1
    xvinsgr2vr.d VX0, t3, 2
    xvinsgr2vr.d VX0, t4, 3
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    xvinsgr2vr.d VX2, t1, 0
    xvinsgr2vr.d VX2, t2, 1
    xvinsgr2vr.d VX2, t3, 2
    xvinsgr2vr.d VX2, t4, 3
    add.d Y, Y, INCY
    xvfmul.d VT0, VX0, VXC
    xvfmul.d VT1, VX2, VXC
    xvstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 3
    add.d XX, XX, INCX
    xvstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 3
    add.d YY, YY, INCY
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    add.d X, X, INCX
    xvinsgr2vr.d VX1, t1, 0
    xvinsgr2vr.d VX1, t2, 1
    xvinsgr2vr.d VX1, t3, 2
    xvinsgr2vr.d VX1, t4, 3
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    xvinsgr2vr.d VX3, t1, 0
    xvinsgr2vr.d VX3, t2, 1
    xvinsgr2vr.d VX3, t3, 2
    xvinsgr2vr.d VX3, t4, 3
    add.d Y, Y, INCY
    xvfmul.d VT0, VX1, VXC
    xvfmul.d VT1, VX3, VXC
    xvstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 3
    add.d XX, XX, INCX
    xvstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 3
    add.d YY, YY, INCY
    addi.d  I, I, -1
    blt $r0, I, .L222
    b .L997
    .align 3

.L223: // C==0 S!=0
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    add.d X, X, INCX
    xvinsgr2vr.d VX0, t1, 0
    xvinsgr2vr.d VX0, t2, 1
    xvinsgr2vr.d VX0, t3, 2
    xvinsgr2vr.d VX0, t4, 3
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    xvinsgr2vr.d VX2, t1, 0
    xvinsgr2vr.d VX2, t2, 1
    xvinsgr2vr.d VX2, t3, 2
    xvinsgr2vr.d VX2, t4, 3
    add.d Y, Y, INCY
    xvfmul.d VT0, VX2, VXS
    xvfmul.d VT1, VX0, VXS
    xvfsub.d VT1, VXZ, VT1
    xvstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 3
    add.d XX, XX, INCX
    xvstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 3
    add.d YY, YY, INCY
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    add.d X, X, INCX
    xvinsgr2vr.d VX1, t1, 0
    xvinsgr2vr.d VX1, t2, 1
    xvinsgr2vr.d VX1, t3, 2
    xvinsgr2vr.d VX1, t4, 3
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    xvinsgr2vr.d VX3, t1, 0
    xvinsgr2vr.d VX3, t2, 1
    xvinsgr2vr.d VX3, t3, 2
    xvinsgr2vr.d VX3, t4, 3
    add.d Y, Y, INCY
    xvfmul.d VT0, VX3, VXS
    xvfmul.d VT1, VX0, VXS
    xvfsub.d VT1, VXZ, VT1
    xvstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.d VT0, XX, 0, 3
    add.d XX, XX, INCX
    xvstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VT1, YY, 0, 3
    add.d YY, YY, INCY
    addi.d  I, I, -1
    blt $r0, I, .L223
    b .L997
    .align 3

.L224: // C==0 S==0
    xvstelm.d VXZ, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.d VXZ, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.d VXZ, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.d VXZ, XX, 0, 3
    add.d XX, XX, INCX
    xvstelm.d VXZ, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VXZ, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VXZ, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VXZ, YY, 0, 3
    add.d YY, YY, INCY
    xvstelm.d VXZ, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.d VXZ, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.d VXZ, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.d VXZ, XX, 0, 3
    add.d XX, XX, INCX
    xvstelm.d VXZ, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VXZ, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VXZ, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VXZ, YY, 0, 3
    add.d YY, YY, INCY
    addi.d  I, I, -1
    blt $r0, I, .L224
    b .L997
    .align 3

.L997:
    andi I, N, 7
    bge $r0, I, .L999
    .align 3

.L998:
    fld.d $f12, X, 0 * SIZE
    fld.d $f13, Y, 0 * SIZE
    fmul.d $f10, $f12, C
    fmadd.d $f10, $f13, S, $f10
    fst.d $f10, X, 0 * SIZE
    addi.d I, I, -1
    fmul.d $f20, $f12, S
    fmsub.d $f20, $f13, C, $f20
    fst.d $f20, Y, 0 * SIZE
    add.d  X, X, INCX
    add.d  Y, Y, INCY
    blt $r0, I, .L998
    .align 3

.L999:
    move $r4, $r12
    jirl $r0, $r1, 0x0
    .align 3

    EPILOGUE