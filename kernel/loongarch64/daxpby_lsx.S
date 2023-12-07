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
#define VX0    $vr8
#define VX1    $vr20
#define VX2    $vr21
#define VX3    $vr22
#define VXA    $vr23
#define VXB    $vr9
#define VXZ    $vr19

    PROLOGUE

    bge $r0, N, .L999
    li.d TEMP, 1
    movgr2fr.d a1, $r0
    ffint.d.l a1, a1
    slli.d  TEMP, TEMP, BASE_SHIFT
    slli.d  INCX, INCX, BASE_SHIFT
    slli.d  INCY, INCY, BASE_SHIFT
    movfr2gr.d t1, ALPHA
    vreplgr2vr.d VXA, t1
    movfr2gr.d t2, BETA
    vreplgr2vr.d VXB, t2
    movfr2gr.d t3, a1
    vreplgr2vr.d VXZ, t3
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
    vld VX0, X, 0 * SIZE
    vld VX2, Y, 0 * SIZE
    vld VX1, X, 2 * SIZE
    vld VX3, Y, 2 * SIZE
    vfmul.d VX0, VX0, VXA
    vfmul.d VX1, VX1, VXA
    vfmadd.d VX2, VX2, VXB, VX0
    vfmadd.d VX3, VX3, VXB, VX1
    vst VX2, Y, 0 * SIZE
    vst VX3, Y, 2 * SIZE
    vld VX0, X, 4 * SIZE
    vld VX2, Y, 4 * SIZE
    vld VX1, X, 6 * SIZE
    vld VX3, Y, 6 * SIZE
    vfmul.d VX0, VX0, VXA
    vfmul.d VX1, VX1, VXA
    vfmadd.d VX2, VX2, VXB, VX0
    vfmadd.d VX3, VX3, VXB, VX1
    vst VX2, Y, 4 * SIZE
    vst VX3, Y, 6 * SIZE
    addi.d X, X, 8 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L111
    b .L997
    .align 3

.L112: // ALPHA!=0 BETA==0
    vld VX0, X, 0 * SIZE
    vld VX1, X, 2 * SIZE
    vfmul.d VX0, VX0, VXA
    vfmul.d VX1, VX1, VXA
    vst VX0, Y, 0 * SIZE
    vst VX1, Y, 2 * SIZE
    vld VX2, X, 4 * SIZE
    vld VX3, X, 6 * SIZE
    vfmul.d VX2, VX2, VXA
    vfmul.d VX3, VX3, VXA
    vst VX2, Y, 4 * SIZE
    vst VX3, Y, 6 * SIZE
    addi.d X, X, 8 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L112
    b .L997
    .align 3

.L113: // ALPHA==0 BETA!=0\
    vld VX0, Y, 0 * SIZE
    vld VX1, Y, 2 * SIZE
    vfmul.d VX0, VX0, VXB
    vfmul.d VX1, VX1, VXB
    vst VX0, Y, 0 * SIZE
    vst VX1, Y, 2 * SIZE
    vld VX2, Y, 4 * SIZE
    vld VX3, Y, 6 * SIZE
    vfmul.d VX2, VX2, VXB
    vfmul.d VX3, VX3, VXB
    vst VX2, Y, 4 * SIZE
    vst VX3, Y, 6 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L113
    b .L997
    .align 3

.L114: // ALPHA==0 BETA==0
    vst VXZ, Y, 0 * SIZE
    vst VXZ, Y, 2 * SIZE
    vst VXZ, Y, 4 * SIZE
    vst VXZ, Y, 6 * SIZE
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
    vld VX0, X, 0 * SIZE
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    vinsgr2vr.d VX2, t1, 0
    vinsgr2vr.d VX2, t2, 1
    vfmul.d VX0, VX0, VXA
    vld VX1, X, 2 * SIZE
    vfmadd.d VX2, VX2, VXB, VX0
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    add.d Y, Y, INCY
    vinsgr2vr.d VX3, t3, 0
    vinsgr2vr.d VX3, t4, 1
    vstelm.d VX2, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX2, YY, 0, 1
    add.d YY, YY, INCY
    vfmul.d VX1, VX1, VXA
    vld VX0, X, 4 * SIZE
    vfmadd.d VX3, VX3, VXB, VX1
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    vinsgr2vr.d VX2, t1, 0
    vinsgr2vr.d VX2, t2, 1
    vstelm.d VX3, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX3, YY, 0, 1
    add.d YY, YY, INCY
    vfmul.d VX0, VX0, VXA
    vld VX1, X, 6 * SIZE
    vfmadd.d VX2, VX2, VXB, VX0
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    add.d Y, Y, INCY
    vinsgr2vr.d VX3, t3, 0
    vinsgr2vr.d VX3, t4, 1
    vstelm.d VX2, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX2, YY, 0, 1
    add.d YY, YY, INCY
    vfmul.d VX1, VX1, VXA
    vfmadd.d VX3, VX3, VXB, VX1
    addi.d  I, I, -1
    vstelm.d VX3, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX3, YY, 0, 1
    add.d YY, YY, INCY
    addi.d X, X, 8 * SIZE
    blt $r0, I, .L121
    b .L997
    .align 3

.L122: // ALPHA!=0 BETA==0
    vld VX0, X, 0 * SIZE
    vld VX1, X, 2 * SIZE
    vfmul.d VX0, VX0, VXA
    vfmul.d VX1, VX1, VXA
    vstelm.d VX0, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX0, YY, 0, 1
    add.d YY, YY, INCY
    vstelm.d VX1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX1, YY, 0, 1
    add.d YY, YY, INCY
    vld VX0, X, 4 * SIZE
    vld VX1, X, 6 * SIZE
    vfmul.d VX0, VX0, VXA
    vfmul.d VX1, VX1, VXA
    vstelm.d VX0, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX0, YY, 0, 1
    add.d YY, YY, INCY
    vstelm.d VX1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX1, YY, 0, 1
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
    vinsgr2vr.d VX2, t1, 0
    vinsgr2vr.d VX2, t2, 1
    add.d Y, Y, INCY
    vfmul.d VX2, VX2, VXB
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    add.d Y, Y, INCY
    vinsgr2vr.d VX3, t3, 0
    vinsgr2vr.d VX3, t4, 1
    vstelm.d VX2, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX2, YY, 0, 1
    add.d YY, YY, INCY
    vfmul.d VX3, VX3, VXB
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    vinsgr2vr.d VX2, t1, 0
    vinsgr2vr.d VX2, t2, 1
    vstelm.d VX3, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX3, YY, 0, 1
    add.d YY, YY, INCY
    vfmul.d VX2, VX2, VXB
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    add.d Y, Y, INCY
    vinsgr2vr.d VX3, t3, 0
    vinsgr2vr.d VX3, t4, 1
    vstelm.d VX2, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX2, YY, 0, 1
    add.d YY, YY, INCY
    vfmul.d VX3, VX3, VXB
    addi.d  I, I, -1
    vstelm.d VX3, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX3, YY, 0, 1
    add.d YY, YY, INCY
    blt $r0, I, .L123
    b .L997
    .align 3

.L124: // ALPHA==0 BETA==0
    vstelm.d VXZ, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VXZ, YY, 0, 1
    add.d YY, YY, INCY
    vstelm.d VXZ, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VXZ, YY, 0, 1
    add.d YY, YY, INCY
    vstelm.d VXZ, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VXZ, YY, 0, 1
    add.d YY, YY, INCY
    vstelm.d VXZ, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VXZ, YY, 0, 1
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
    vld VX2, Y, 0 * SIZE
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    vinsgr2vr.d VX0, t1, 0
    vinsgr2vr.d VX0, t2, 1
    add.d X, X, INCX
    vfmul.d VX0, VXA, VX0
    vld VX3, Y, 2 * SIZE
    vfmadd.d VX2, VX2, VXB, VX0
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX1, t3, 0
    vinsgr2vr.d VX1, t4, 1
    vst VX2, Y, 0 * SIZE
    vfmul.d VX1, VXA, VX1
    vld VX2, Y, 4 * SIZE
    vfmadd.d VX3, VX3, VXB, VX1
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX0, t1, 0
    vinsgr2vr.d VX0, t2, 1
    vst VX3, Y, 2 * SIZE
    vfmul.d VX0, VX0, VXA
    vld VX3, Y, 6 * SIZE
    vfmadd.d VX2, VX2, VXB, VX0
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX1, t3, 0
    vinsgr2vr.d VX1, t4, 1
    vst VX2, Y, 4 * SIZE
    vfmul.d VX1, VX1, VXA
    vfmadd.d VX3, VX3, VXB, VX1
    addi.d  I, I, -1
    vst VX3, Y, 6 * SIZE
    addi.d Y, Y, 8 * SIZE
    blt $r0, I, .L211
    b .L997
    .align 3

.L212: // ALPHA!=0 BETA==0
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    vinsgr2vr.d VX0, t1, 0
    vinsgr2vr.d VX0, t2, 1
    add.d X, X, INCX
    vfmul.d VX0, VXA, VX0
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX1, t3, 0
    vinsgr2vr.d VX1, t4, 1
    vst VX0, Y, 0 * SIZE
    vfmul.d VX1, VXA, VX1
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX0, t1, 0
    vinsgr2vr.d VX0, t2, 1
    vst VX1, Y, 2 * SIZE
    vfmul.d VX0, VX0, VXA
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX1, t3, 0
    vinsgr2vr.d VX1, t4, 1
    vst VX0, Y, 4 * SIZE
    vfmul.d VX1, VX1, VXA
    addi.d  I, I, -1
    vst VX1, Y, 6 * SIZE
    addi.d Y, Y, 8 * SIZE
    blt $r0, I, .L212
    b .L997
    .align 3

.L213: // ALPHA==0 BETA!=0
    vld VX2, Y, 0 * SIZE
    vld VX3, Y, 2 * SIZE
    vfmul.d VX2, VX2, VXB
    vfmul.d VX3, VX3, VXB
    vst VX2, Y, 0 * SIZE
    vst VX3, Y, 2 * SIZE
    vld VX2, Y, 4 * SIZE
    vld VX3, Y, 6 * SIZE
    vfmul.d VX2, VX2, VXB
    vfmul.d VX3, VX3, VXB
    addi.d  I, I, -1
    vst VX2, Y, 4 * SIZE
    vst VX3, Y, 6 * SIZE
    addi.d Y, Y, 8 * SIZE
    blt $r0, I, .L213
    b .L997
    .align 3

.L214: // ALPHA==0 BETA==0
    vst VXZ, Y, 0 * SIZE
    vst VXZ, Y, 2 * SIZE
    vst VXZ, Y, 4 * SIZE
    vst VXZ, Y, 6 * SIZE
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
    vinsgr2vr.d VX0, t1, 0
    vinsgr2vr.d VX0, t2, 1
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    vinsgr2vr.d VX2, t3, 0
    vinsgr2vr.d VX2, t4, 1
    add.d Y, Y, INCY
    vfmul.d VX0, VX0, VXA
    vfmadd.d VX2, VX2, VXB, VX0
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX1, t1, 0
    vinsgr2vr.d VX1, t2, 1
    vstelm.d VX2, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX2, YY, 0, 1
    add.d YY, YY, INCY
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    vinsgr2vr.d VX3, t3, 0
    vinsgr2vr.d VX3, t4, 1
    add.d Y, Y, INCY
    vfmul.d VX1, VX1, VXA
    vfmadd.d VX3, VX3, VXB, VX1
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX0, t1, 0
    vinsgr2vr.d VX0, t2, 1
    vstelm.d VX3, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX3, YY, 0, 1
    add.d YY, YY, INCY
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    vinsgr2vr.d VX2, t3, 0
    vinsgr2vr.d VX2, t4, 1
    add.d Y, Y, INCY
    vfmul.d VX0, VX0, VXA
    vfmadd.d VX2, VX2, VXB, VX0
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX1, t3, 0
    vinsgr2vr.d VX1, t4, 1
    vstelm.d VX2, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX2, YY, 0, 1
    add.d YY, YY, INCY
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    vinsgr2vr.d VX3, t1, 0
    vinsgr2vr.d VX3, t2, 1
    add.d Y, Y, INCY
    vfmul.d VX1, VX1, VXA
    vfmadd.d VX3, VX3, VXB, VX1
    addi.d  I, I, -1
    vstelm.d VX3, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX3, YY, 0, 1
    add.d YY, YY, INCY
    blt $r0, I, .L221
    b .L997
    .align 3

.L222: // ALPHA!=0 BETA==0
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    vinsgr2vr.d VX0, t1, 0
    vinsgr2vr.d VX0, t2, 1
    add.d X, X, INCX
    vfmul.d VX0, VX0, VXA
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX1, t3, 0
    vinsgr2vr.d VX1, t4, 1
    vstelm.d VX0, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX0, YY, 0, 1
    add.d YY, YY, INCY
    vfmul.d VX1, VX1, VXA
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX0, t1, 0
    vinsgr2vr.d VX0, t2, 1
    vstelm.d VX1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX1, YY, 0, 1
    add.d YY, YY, INCY
    vfmul.d VX0, VX0, VXA
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX1, t3, 0
    vinsgr2vr.d VX1, t4, 1
    vstelm.d VX0, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX0, YY, 0, 1
    add.d YY, YY, INCY
    vfmul.d VX1, VX1, VXA
    addi.d  I, I, -1
    vstelm.d VX1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX1, YY, 0, 1
    add.d YY, YY, INCY
    blt $r0, I, .L222
    b .L997
    .align 3

.L223: // ALPHA==0 BETA!=0
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    vinsgr2vr.d VX2, t1, 0
    vinsgr2vr.d VX2, t2, 1
    add.d Y, Y, INCY
    vfmul.d VX2, VX2, VXB
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    add.d Y, Y, INCY
    vinsgr2vr.d VX3, t3, 0
    vinsgr2vr.d VX3, t4, 1
    vstelm.d VX2, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX2, YY, 0, 1
    add.d YY, YY, INCY
    vfmul.d VX3, VX3, VXB
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    vinsgr2vr.d VX2, t1, 0
    vinsgr2vr.d VX2, t2, 1
    vstelm.d VX3, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX3, YY, 0, 1
    add.d YY, YY, INCY
    vfmul.d VX2, VX2, VXB
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    add.d Y, Y, INCY
    vinsgr2vr.d VX3, t3, 0
    vinsgr2vr.d VX3, t4, 1
    vstelm.d VX2, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX2, YY, 0, 1
    add.d YY, YY, INCY
    vfmul.d VX3, VX3, VXB
    addi.d  I, I, -1
    vstelm.d VX3, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VX3, YY, 0, 1
    add.d YY, YY, INCY
    blt $r0, I, .L223
    b .L997
    .align 3

.L224: // ALPHA==0 BETA==0
    vstelm.d VXZ, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VXZ, YY, 0, 1
    add.d YY, YY, INCY
    vstelm.d VXZ, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VXZ, YY, 0, 1
    add.d YY, YY, INCY
    vstelm.d VXZ, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VXZ, YY, 0, 1
    add.d YY, YY, INCY
    vstelm.d VXZ, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VXZ, YY, 0, 1
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
