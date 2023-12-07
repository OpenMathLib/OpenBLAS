#define ASSEMBLER

#include "common.h"
#define N      $r4
#define ALPHA  $f0
#define X      $r5
#define INCX   $r6
#define BETA   $f1
#define Y      $r7
#define INCY   $r8

#define I      $r12
#define TEMP   $r13
#define t1     $r14
#define t2     $r16
#define t3     $r15
#define t4     $r17
#define XX     $r18
#define YY     $r19
#define a1     $f12
#define a2     $f13
#define VX0    $xr8
#define VX1    $xr20
#define VX2    $xr21
#define VX3    $xr22
#define VXA    $xr23
#define VXB    $xr9
#define VXZ    $xr19

    PROLOGUE

    bge $r0, N, .L999
    li.d TEMP, 1
    movgr2fr.d a1, $r0
    ffint.d.l a1, a1
    slli.d  TEMP, TEMP, BASE_SHIFT
    slli.d  INCX, INCX, BASE_SHIFT
    slli.d  INCY, INCY, BASE_SHIFT
    movfr2gr.d t1, ALPHA
    xvreplgr2vr.d VXA, t1
    movfr2gr.d t2, BETA
    xvreplgr2vr.d VXB, t2
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
    fcmp.ceq.d $fcc0, ALPHA, a1
    bcnez $fcc0, .L110
    fcmp.ceq.d $fcc0, BETA, a1
    bcnez $fcc0, .L112 // ALPHA!=0 BETA==0
    b .L111 // ALPHA!=0 BETA!=0
    .align 3

.L110:
    fcmp.ceq.d $fcc0, BETA, a1
    bcnez $fcc0, .L114 // ALPHA==0 BETA==0
    b .L113 // ALPHA==0 BETA!=0
    .align 3

.L111: // ALPHA!=0 BETA!=0
    xvld VX0, X, 0 * SIZE
    xvld VX2, Y, 0 * SIZE
    xvld VX1, X, 4 * SIZE
    xvld VX3, Y, 4 * SIZE
    xvfmul.d VX0, VX0, VXA
    xvfmul.d VX1, VX1, VXA
    xvfmadd.d VX2, VX2, VXB, VX0
    xvfmadd.d VX3, VX3, VXB, VX1
    addi.d  I, I, -1
    xvst VX2, Y, 0 * SIZE
    xvst VX3, Y, 4 * SIZE
    addi.d X, X, 8 * SIZE
    addi.d Y, Y, 8 * SIZE
    blt $r0, I, .L111
    b .L997
    .align 3

.L112: // ALPHA!=0 BETA==0
    xvld VX0, X, 0 * SIZE
    xvld VX1, X, 4 * SIZE
    xvfmul.d VX0, VX0, VXA
    xvfmul.d VX1, VX1, VXA
    xvst VX0, Y, 0 * SIZE
    xvst VX1, Y, 4 * SIZE
    addi.d X, X, 8 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L112
    b .L997
    .align 3

.L113: // ALPHA==0 BETA!=0
    xvld VX2, Y, 0 * SIZE
    xvld VX3, Y, 4 * SIZE
    xvfmul.d VX2, VX2, VXB
    xvfmul.d VX3, VX3, VXB
    xvst VX2, Y, 0 * SIZE
    xvst VX3, Y, 4 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L113
    b .L997
    .align 3

.L114: // ALPHA==0 BETA==0
    xvst VXZ, Y, 0 * SIZE
    xvst VXZ, Y, 4 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L114
    b .L997
    .align 3

.L12: // INCX==1 and INCY!=1
    bge $r0, I, .L997
    move YY, Y
    fcmp.ceq.d $fcc0, ALPHA, a1
    bcnez $fcc0, .L120
    fcmp.ceq.d $fcc0, BETA, a1
    bcnez $fcc0, .L122 // ALPHA!=0 BETA==0
    b .L121 // ALPHA!=0 BETA!=0
    .align 3

.L120:
    fcmp.ceq.d $fcc0, BETA, a1
    bcnez $fcc0, .L124 // ALPHA==0 BETA==0
    b .L123 // ALPHA==0 BETA!=0
    .align 3

.L121: // ALPHA!=0 BETA!=0
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
    xvfmul.d VX0, VX0, VXA
    xvld VX1, X, 4 * SIZE
    xvfmadd.d VX2, VX2, VXB, VX0
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    add.d Y, Y, INCY
    xvinsgr2vr.d VX3, t1, 0
    xvinsgr2vr.d VX3, t2, 1
    xvinsgr2vr.d VX3, t3, 2
    xvinsgr2vr.d VX3, t4, 3
    xvstelm.d VX2, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VX2, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VX2, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VX2, YY, 0, 3
    add.d YY, YY, INCY
    xvfmul.d VX1, VX1, VXA
    xvfmadd.d VX3, VX3, VXB, VX1
    addi.d  I, I, -1
    xvstelm.d VX3, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VX3, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VX3, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VX3, YY, 0, 3
    add.d YY, YY, INCY
    addi.d X, X, 8 * SIZE
    blt $r0, I, .L121
    b .L997
    .align 3

.L122: // ALPHA!=0 BETA==0
    xvld VX0, X, 0 * SIZE
    xvld VX1, X, 4 * SIZE
    xvfmul.d VX0, VX0, VXA
    xvfmul.d VX1, VX1, VXA
    xvstelm.d VX0, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VX0, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VX0, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VX0, YY, 0, 3
    add.d YY, YY, INCY
    xvstelm.d VX1, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VX1, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VX1, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VX1, YY, 0, 3
    add.d YY, YY, INCY
    addi.d X, X, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L122
    b .L997
    .align 3

.L123: // ALPHA==0 BETA!=0
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
    xvfmul.d VX2, VX2, VXB
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    add.d Y, Y, INCY
    xvinsgr2vr.d VX3, t1, 0
    xvinsgr2vr.d VX3, t2, 1
    xvinsgr2vr.d VX3, t3, 2
    xvinsgr2vr.d VX3, t4, 3
    xvstelm.d VX2, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VX2, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VX2, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VX2, YY, 0, 3
    add.d YY, YY, INCY
    xvfmul.d VX3, VX3, VXB
    addi.d  I, I, -1
    xvstelm.d VX3, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VX3, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VX3, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VX3, YY, 0, 3
    add.d YY, YY, INCY
    blt $r0, I, .L123
    b .L997
    .align 3

.L124: // ALPHA==0 BETA==0
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
    fcmp.ceq.d $fcc0, ALPHA, a1
    bcnez $fcc0, .L210
    fcmp.ceq.d $fcc0, BETA, a1
    bcnez $fcc0, .L212 // ALPHA!=0 BETA==0
    b .L211 // ALPHA!=0 BETA!=0
    .align 3

.L210:
    fcmp.ceq.d $fcc0, BETA, a1
    bcnez $fcc0, .L214 // ALPHA==0 BETA==0
    b .L213 // ALPHA==0 BETA!=0
    .align 3

.L211: // ALPHA!=0 BETA!=0
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
    xvfmul.d VX0, VXA, VX0
    xvfmadd.d VX2, VX2, VXB, VX0
    xvld VX3, Y, 4 * SIZE
    xvst VX2, Y, 0 * SIZE
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
    xvfmul.d VX1, VX1, VXA
    xvfmadd.d VX3, VX3, VXB, VX1
    addi.d  I, I, -1
    xvst VX3, Y, 4 * SIZE
    addi.d Y, Y, 8 * SIZE
    blt $r0, I, .L211
    b .L997
    .align 3

.L212: // ALPHA!=0 BETA==0
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
    xvfmul.d VX0, VXA, VX0
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
    xvst VX0, Y, 0 * SIZE
    xvfmul.d VX1, VX1, VXA
    addi.d  I, I, -1
    xvst VX1, Y, 4 * SIZE
    addi.d Y, Y, 8 * SIZE
    blt $r0, I, .L212
    b .L997
    .align 3

.L213: // ALPHA==0 BETA!=0
    xvld VX2, Y, 0 * SIZE
    xvld VX3, Y, 4 * SIZE
    xvfmul.d VX2, VX2, VXB
    xvfmul.d VX3, VX3, VXB
    addi.d  I, I, -1
    xvst VX2, Y, 0 * SIZE
    xvst VX3, Y, 4 * SIZE
    addi.d Y, Y, 8 * SIZE
    blt $r0, I, .L213
    b .L997
    .align 3

.L214: // ALPHA==0 BETA==0
    xvst VXZ, Y, 0 * SIZE
    xvst VXZ, Y, 4 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L214
    b .L997
    .align 3

.L22:
    bge $r0, I, .L997
    move YY, Y
    fcmp.ceq.d $fcc0, ALPHA, a1
    bcnez $fcc0, .L220
    fcmp.ceq.d $fcc0, BETA, a1
    bcnez $fcc0, .L222 // ALPHA!=0 BETA==0
    b .L221 // ALPHA!=0 BETA!=0
    .align 3

.L220:
    fcmp.ceq.d $fcc0, BETA, a1
    bcnez $fcc0, .L224 // ALPHA==0 BETA==0
    b .L223 // ALPHA==0 BETA!=0
    .align 3

.L221: // ALPHA!=0 BETA!=0
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
    xvfmul.d VX0, VX0, VXA
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    xvfmadd.d VX2, VX2, VXB, VX0
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
    xvstelm.d VX2, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VX2, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VX2, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VX2, YY, 0, 3
    add.d YY, YY, INCY
    ld.d t1, Y, 0 * SIZE
    xvinsgr2vr.d VX3, t1, 0
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    xvinsgr2vr.d VX3, t2, 1
    xvinsgr2vr.d VX3, t3, 2
    xvinsgr2vr.d VX3, t4, 3
    add.d Y, Y, INCY
    xvfmul.d VX1, VX1, VXA
    xvfmadd.d VX3, VX3, VXB, VX1
    addi.d  I, I, -1
    xvstelm.d VX3, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VX3, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VX3, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VX3, YY, 0, 3
    add.d YY, YY, INCY
    blt $r0, I, .L221
    b .L997
    .align 3

.L222: // ALPHA!=0 BETA==0
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
    xvfmul.d VX0, VX0, VXA
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
    xvstelm.d VX0, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VX0, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VX0, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VX0, YY, 0, 3
    add.d YY, YY, INCY
    xvfmul.d VX1, VX1, VXA
    addi.d  I, I, -1
    xvstelm.d VX1, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VX1, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VX1, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VX1, YY, 0, 3
    add.d YY, YY, INCY
    blt $r0, I, .L222
    b .L997
    .align 3

.L223: // ALPHA==0 BETA!=0
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
    xvfmul.d VX2, VX2, VXB
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    add.d Y, Y, INCY
    xvinsgr2vr.d VX3, t1, 0
    xvinsgr2vr.d VX3, t2, 1
    xvinsgr2vr.d VX3, t3, 2
    xvinsgr2vr.d VX3, t4, 3
    xvstelm.d VX2, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VX2, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VX2, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VX2, YY, 0, 3
    add.d YY, YY, INCY
    xvfmul.d VX3, VX3, VXB
    addi.d  I, I, -1
    xvstelm.d VX3, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.d VX3, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.d VX3, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.d VX3, YY, 0, 3
    add.d YY, YY, INCY
    blt $r0, I, .L223
    b .L997
    .align 3

.L224: // ALPHA==0 BETA==0
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
    addi.d I, I, -1
    fmul.d $f12, $f12, ALPHA
    fmadd.d $f13, $f13, BETA, $f12
    fst.d $f13, Y, 0 * SIZE
    add.d  X, X, INCX
    add.d  Y, Y, INCY
    blt $r0, I, .L998
    .align 3

.L999:
    move $r4, $r12
    jirl $r0, $r1, 0x0
    .align 3

    EPILOGUE
