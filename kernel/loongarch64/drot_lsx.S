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
#define VX0    $vr8
#define VX1    $vr20
#define VX2    $vr21
#define VX3    $vr22
#define VT0    $vr10
#define VT1    $vr18
#define VXC    $vr23
#define VXS    $vr9
#define VXZ    $vr19

    PROLOGUE

    bge $r0, N, .L999
    li.d TEMP, 1
    movgr2fr.d a1, $r0
    ffint.d.l a1, a1
    slli.d  TEMP, TEMP, BASE_SHIFT
    slli.d  INCX, INCX, BASE_SHIFT
    slli.d  INCY, INCY, BASE_SHIFT
    movfr2gr.d t1, C
    vreplgr2vr.d VXC, t1
    movfr2gr.d t2, S
    vreplgr2vr.d VXS, t2
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
    vld VX0, X, 0 * SIZE
    vld VX2, Y, 0 * SIZE
    vld VX1, X, 2 * SIZE
    vld VX3, Y, 2 * SIZE
    vfmul.d VT0, VX0, VXC
    vfmadd.d VT0, VX2, VXS, VT0
    vfmul.d VT1, VX0, VXS
    vfmsub.d VT1, VX2, VXC, VT1
    vst VT0, X, 0 * SIZE
    vst VT1, Y, 0 * SIZE
    vfmul.d VT0, VX1, VXC
    vfmadd.d VT0, VX3, VXS, VT0
    vfmul.d VT1, VX1, VXS
    vfmsub.d VT1, VX3, VXC, VT1
    vst VT0, X, 2 * SIZE
    vst VT1, Y, 2 * SIZE
    vld VX0, X, 4 * SIZE
    vld VX2, Y, 4 * SIZE
    vld VX1, X, 6 * SIZE
    vld VX3, Y, 6 * SIZE
    vfmul.d VT0, VX0, VXC
    vfmadd.d VT0, VX2, VXS, VT0
    vfmul.d VT1, VX0, VXS
    vfmsub.d VT1, VX2, VXC, VT1
    vst VT0, X, 4 * SIZE
    vst VT1, Y, 4 * SIZE
    vfmul.d VT0, VX1, VXC
    vfmadd.d VT0, VX3, VXS, VT0
    vfmul.d VT1, VX1, VXS
    vfmsub.d VT1, VX3, VXC, VT1
    vst VT0, X, 6 * SIZE
    vst VT1, Y, 6 * SIZE
    addi.d X, X, 8 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L111
    b .L997
    .align 3

.L112: // C!=0 S==0
    vld VX0, X, 0 * SIZE
    vld VX2, Y, 0 * SIZE
    vld VX1, X, 2 * SIZE
    vld VX3, Y, 2 * SIZE
    vfmul.d VT0, VX0, VXC
    vfmul.d VT1, VX2, VXC
    vst VT0, X, 0 * SIZE
    vst VT1, Y, 0 * SIZE
    vfmul.d VT0, VX1, VXC
    vfmul.d VT1, VX3, VXC
    vst VT0, X, 2 * SIZE
    vst VT1, Y, 2 * SIZE
    vld VX0, X, 4 * SIZE
    vld VX2, Y, 4 * SIZE
    vld VX1, X, 6 * SIZE
    vld VX3, Y, 6 * SIZE
    vfmul.d VT0, VX0, VXC
    vfmul.d VT1, VX2, VXC
    vst VT0, X, 4 * SIZE
    vst VT1, Y, 4 * SIZE
    vfmul.d VT0, VX1, VXC
    vfmul.d VT1, VX3, VXC
    vst VT0, X, 6 * SIZE
    vst VT1, Y, 6 * SIZE
    addi.d X, X, 8 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L112
    b .L997
    .align 3

.L113: // C==0 S!=0
    vld VX0, X, 0 * SIZE
    vld VX2, Y, 0 * SIZE
    vld VX1, X, 2 * SIZE
    vld VX3, Y, 2 * SIZE
    vfmul.d VT0, VX2, VXS
    vfmul.d VT1, VX0, VXS
    vfsub.d VT1, VXZ, VT1
    vst VT0, X, 0 * SIZE
    vst VT1, Y, 0 * SIZE
    vfmul.d VT0, VX3, VXS
    vfmul.d VT1, VX1, VXS
    vfsub.d VT1, VXZ, VT1
    vst VT0, X, 2 * SIZE
    vst VT1, Y, 2 * SIZE
    vld VX0, X, 4 * SIZE
    vld VX2, Y, 4 * SIZE
    vld VX1, X, 6 * SIZE
    vld VX3, Y, 6 * SIZE
    vfmul.d VT0, VX2, VXS
    vfmul.d VT1, VX0, VXS
    vfsub.d VT1, VXZ, VT1
    vst VT0, X, 4 * SIZE
    vst VT1, Y, 4 * SIZE
    vfmul.d VT0, VX3, VXS
    vfmul.d VT1, VX1, VXS
    vfsub.d VT1, VXZ, VT1
    vst VT0, X, 6 * SIZE
    vst VT1, Y, 6 * SIZE
    addi.d X, X, 8 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L113
    b .L997
    .align 3

.L114: // C==0 S==0
    vst VXZ, X, 0 * SIZE
    vst VXZ, Y, 0 * SIZE
    vst VXZ, X, 2 * SIZE
    vst VXZ, Y, 2 * SIZE
    vst VXZ, X, 4 * SIZE
    vst VXZ, Y, 4 * SIZE
    vst VXZ, X, 6 * SIZE
    vst VXZ, Y, 6 * SIZE
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
    vld VX0, X, 0 * SIZE
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    vinsgr2vr.d VX2, t1, 0
    vinsgr2vr.d VX2, t2, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX0, VXC
    vfmadd.d VT0, VX2, VXS, VT0
    vfmul.d VT1, VX0, VXS
    vfmsub.d VT1, VX2, VXC, VT1
    vst VT0, X, 0 * SIZE
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    vld VX0, X, 2 * SIZE
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    vinsgr2vr.d VX2, t3, 0
    vinsgr2vr.d VX2, t4, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX0, VXC
    vfmadd.d VT0, VX2, VXS, VT0
    vfmul.d VT1, VX0, VXS
    vfmsub.d VT1, VX2, VXC, VT1
    vst VT0, X, 2 * SIZE
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    vld VX1, X, 4 * SIZE
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    vinsgr2vr.d VX3, t1, 0
    vinsgr2vr.d VX3, t2, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX1, VXC
    vfmadd.d VT0, VX3, VXS, VT0
    vfmul.d VT1, VX1, VXS
    vfmsub.d VT1, VX3, VXC, VT1
    vst VT0, X, 4 * SIZE
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    vld VX1, X, 6 * SIZE
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    vinsgr2vr.d VX3, t3, 0
    vinsgr2vr.d VX3, t4, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX1, VXC
    vfmadd.d VT0, VX3, VXS, VT0
    vfmul.d VT1, VX1, VXS
    vfmsub.d VT1, VX3, VXC, VT1
    vst VT0, X, 6 * SIZE
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    addi.d X, X, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L121
    b .L997
    .align 3

.L122: // C!=0 S==0
    vld VX0, X, 0 * SIZE
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    vinsgr2vr.d VX2, t1, 0
    vinsgr2vr.d VX2, t2, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX0, VXC
    vfmul.d VT1, VX2, VXC
    vst VT0, X, 0 * SIZE
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    vld VX0, X, 2 * SIZE
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    vinsgr2vr.d VX2, t3, 0
    vinsgr2vr.d VX2, t4, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX0, VXC
    vfmul.d VT1, VX2, VXC
    vst VT0, X, 2 * SIZE
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    vld VX1, X, 4 * SIZE
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    vinsgr2vr.d VX3, t1, 0
    vinsgr2vr.d VX3, t2, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX1, VXC
    vfmul.d VT1, VX3, VXC
    vst VT0, X, 4 * SIZE
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    vld VX1, X, 6 * SIZE
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    vinsgr2vr.d VX3, t3, 0
    vinsgr2vr.d VX3, t4, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX1, VXC
    vfmul.d VT1, VX3, VXC
    vst VT0, X, 6 * SIZE
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    addi.d X, X, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L122
    b .L997
    .align 3

.L123: // C==0 S!=0
    vld VX0, X, 0 * SIZE
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    vinsgr2vr.d VX2, t1, 0
    vinsgr2vr.d VX2, t2, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX2, VXS
    vfmul.d VT1, VX0, VXS
    vfsub.d VT1, VXZ, VT1
    vst VT0, X, 0 * SIZE
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    vld VX0, X, 2 * SIZE
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    vinsgr2vr.d VX2, t3, 0
    vinsgr2vr.d VX2, t4, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX2, VXS
    vfmul.d VT1, VX0, VXS
    vfsub.d VT1, VXZ, VT1
    vst VT0, X, 2 * SIZE
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    vld VX1, X, 4 * SIZE
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    vinsgr2vr.d VX3, t1, 0
    vinsgr2vr.d VX3, t2, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX3, VXS
    vfmul.d VT1, VX1, VXS
    vfsub.d VT1, VXZ, VT1
    vst VT0, X, 4 * SIZE
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    vld VX1, X, 6 * SIZE
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    vinsgr2vr.d VX3, t3, 0
    vinsgr2vr.d VX3, t4, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX3, VXS
    vfmul.d VT1, VX1, VXS
    vfsub.d VT1, VXZ, VT1
    vst VT0, X, 6 * SIZE
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    addi.d X, X, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L123
    b .L997
    .align 3

.L124: // C==0 S==0
    vst VXZ, X, 0 * SIZE
    vst VXZ, X, 4 * SIZE
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
    vld VX2, Y, 0 * SIZE
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    vinsgr2vr.d VX0, t1, 0
    vinsgr2vr.d VX0, t2, 1
    add.d X, X, INCX
    vfmul.d VT0, VXC, VX0
    vfmadd.d VT0, VX2, VXS, VT0
    vfmul.d VT1, VXS, VX0
    vfmsub.d VT1, VX2, VXC, VT1
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vst VT1, Y, 0 * SIZE
    vld VX2, Y, 2 * SIZE
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    vinsgr2vr.d VX0, t3, 0
    vinsgr2vr.d VX0, t4, 1
    add.d X, X, INCX
    vfmul.d VT0, VXC, VX0
    vfmadd.d VT0, VX2, VXS, VT0
    vfmul.d VT1, VXS, VX0
    vfmsub.d VT1, VX2, VXC, VT1
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vst VT1, Y, 2 * SIZE
    vld VX3, Y, 4 * SIZE
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    vinsgr2vr.d VX1, t1, 0
    vinsgr2vr.d VX1, t2, 1
    add.d X, X, INCX
    vfmul.d VT0, VX1, VXC
    vfmadd.d VT0, VX3, VXS, VT0
    vfmul.d VT1, VX1, VXS
    vfmsub.d VT1, VX3, VXC, VT1
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vst VT1, Y, 4 * SIZE
    vld VX3, Y, 6 * SIZE
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    vinsgr2vr.d VX1, t3, 0
    vinsgr2vr.d VX1, t4, 1
    add.d X, X, INCX
    vfmul.d VT0, VX1, VXC
    vfmadd.d VT0, VX3, VXS, VT0
    vfmul.d VT1, VX1, VXS
    vfmsub.d VT1, VX3, VXC, VT1
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vst VT1, Y, 6 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L211
    b .L997
    .align 3

.L212: // C!=0 S==0
    vld VX2, Y, 0 * SIZE
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    vinsgr2vr.d VX0, t1, 0
    vinsgr2vr.d VX0, t2, 1
    add.d X, X, INCX
    vfmul.d VT0, VXC, VX0
    vfmul.d VT1, VX2, VXC
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vst VT1, Y, 0 * SIZE
    vld VX2, Y, 2 * SIZE
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    vinsgr2vr.d VX0, t3, 0
    vinsgr2vr.d VX0, t4, 1
    add.d X, X, INCX
    vfmul.d VT0, VXC, VX0
    vfmul.d VT1, VX2, VXC
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vst VT1, Y, 2 * SIZE
    vld VX3, Y, 4 * SIZE
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    vinsgr2vr.d VX1, t1, 0
    vinsgr2vr.d VX1, t2, 1
    add.d X, X, INCX
    vfmul.d VT0, VX1, VXC
    vfmul.d VT1, VX3, VXS
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vst VT1, Y, 4 * SIZE
    vld VX3, Y, 6 * SIZE
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    vinsgr2vr.d VX1, t3, 0
    vinsgr2vr.d VX1, t4, 1
    add.d X, X, INCX
    vfmul.d VT0, VX1, VXC
    vfmul.d VT1, VX3, VXS
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    vst VT1, Y, 6 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L212
    b .L997
    .align 3

.L213: // C==0 S!=0
    vld VX2, Y, 0 * SIZE
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    vinsgr2vr.d VX0, t1, 0
    vinsgr2vr.d VX0, t2, 1
    add.d X, X, INCX
    vfmul.d VT0, VXS, VX2
    vfmul.d VT1, VXS, VX0
    vfsub.d VT1, VXZ, VT1
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vst VT1, Y, 0 * SIZE
    vld VX2, Y, 2 * SIZE
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    vinsgr2vr.d VX0, t3, 0
    vinsgr2vr.d VX0, t4, 1
    add.d X, X, INCX
    vfmul.d VT0, VXS, VX2
    vfmul.d VT1, VXS, VX0
    vfsub.d VT1, VXZ, VT1
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vst VT1, Y, 2 * SIZE
    vld VX3, Y, 4 * SIZE
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    vinsgr2vr.d VX1, t1, 0
    vinsgr2vr.d VX1, t2, 1
    add.d X, X, INCX
    vfmul.d VT0, VX3, VXS
    vfmul.d VT1, VX1, VXS
    vfsub.d VT1, VXZ, VT1
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vst VT1, Y, 4 * SIZE
    vld VX3, Y, 6 * SIZE
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    vinsgr2vr.d VX1, t3, 0
    vinsgr2vr.d VX1, t4, 1
    add.d X, X, INCX
    vfmul.d VT0, VX3, VXS
    vfmul.d VT1, VX1, VXS
    vfsub.d VT1, VXZ, VT1
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vst VT1, Y, 6 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L213
    b .L997
    .align 3

.L214: // C==0 S==0
    vstelm.d VXZ, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VXZ, XX, 0, 1
    add.d XX, XX, INCX
    vstelm.d VXZ, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VXZ, XX, 0, 1
    add.d XX, XX, INCX
    vst VT1, Y, 0 * SIZE
    vstelm.d VXZ, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VXZ, XX, 0, 1
    add.d XX, XX, INCX
    vstelm.d VXZ, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VXZ, XX, 0, 1
    add.d XX, XX, INCX
    vst VT1, Y, 4 * SIZE
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
    vinsgr2vr.d VX0, t1, 0
    vinsgr2vr.d VX0, t2, 1
    add.d X, X, INCX
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    vinsgr2vr.d VX2, t1, 0
    vinsgr2vr.d VX2, t2, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX0, VXC
    vfmadd.d VT0, VX2, VXS, VT0
    vfmul.d VT1, VX0, VXS
    vfmsub.d VT1, VX2, VXC, VT1
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX0, t3, 0
    vinsgr2vr.d VX0, t4, 1
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    vinsgr2vr.d VX2, t3, 0
    vinsgr2vr.d VX2, t4, 1
    add.d Y, Y, INCY
    vfmul.d  VT0, VX0, VXC
    vfmadd.d VT0, VX2, VXS, VT0
    vfmul.d  VT1, VX0, VXS
    vfmsub.d VT1, VX2, VXC, VT1
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX1, t1, 0
    vinsgr2vr.d VX1, t2, 1
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    vinsgr2vr.d VX3, t1, 0
    vinsgr2vr.d VX3, t2, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX1, VXC
    vfmadd.d VT0, VX3, VXS, VT0
    vfmul.d VT1, VX0, VXS
    vfmsub.d VT1, VX3, VXC, VT1
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    vinsgr2vr.d VX1, t3, 0
    vinsgr2vr.d VX1, t4, 1
    add.d X, X, INCX
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    vinsgr2vr.d VX3, t3, 0
    vinsgr2vr.d VX3, t4, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX1, VXC
    vfmadd.d VT0, VX3, VXS, VT0
    vfmul.d VT1, VX0, VXS
    vfmsub.d VT1, VX3, VXC, VT1
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
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
    vinsgr2vr.d VX0, t1, 0
    vinsgr2vr.d VX0, t2, 1
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    vinsgr2vr.d VX2, t1, 0
    vinsgr2vr.d VX2, t2, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX0, VXC
    vfmul.d VT1, VX2, VXC
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX0, t3, 0
    vinsgr2vr.d VX0, t4, 1
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    vinsgr2vr.d VX2, t3, 0
    vinsgr2vr.d VX2, t4, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX0, VXC
    vfmul.d VT1, VX2, VXC
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX1, t1, 0
    vinsgr2vr.d VX1, t2, 1
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    vinsgr2vr.d VX3, t1, 0
    vinsgr2vr.d VX3, t2, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX1, VXC
    vfmul.d VT1, VX3, VXC
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX1, t3, 0
    vinsgr2vr.d VX1, t4, 1
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    vinsgr2vr.d VX3, t3, 0
    vinsgr2vr.d VX3, t4, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX1, VXC
    vfmul.d VT1, VX3, VXC
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
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
    vinsgr2vr.d VX0, t1, 0
    vinsgr2vr.d VX0, t2, 1
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    vinsgr2vr.d VX2, t1, 0
    vinsgr2vr.d VX2, t2, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX2, VXS
    vfmul.d VT1, VX0, VXS
    vfsub.d VT1, VXZ, VT1
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX0, t3, 0
    vinsgr2vr.d VX0, t4, 1
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    vinsgr2vr.d VX2, t3, 0
    vinsgr2vr.d VX2, t4, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX2, VXS
    vfmul.d VT1, VX0, VXS
    vfsub.d VT1, VXZ, VT1
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX1, t1, 0
    vinsgr2vr.d VX1, t2, 1
    ld.d t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t2, Y, 0 * SIZE
    vinsgr2vr.d VX3, t1, 0
    vinsgr2vr.d VX3, t2, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX3, VXS
    vfmul.d VT1, VX0, VXS
    vfsub.d VT1, VXZ, VT1
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX1, t3, 0
    vinsgr2vr.d VX1, t4, 1
    ld.d t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.d t4, Y, 0 * SIZE
    vinsgr2vr.d VX3, t3, 0
    vinsgr2vr.d VX3, t4, 1
    add.d Y, Y, INCY
    vfmul.d VT0, VX3, VXS
    vfmul.d VT1, VX0, VXS
    vfsub.d VT1, VXZ, VT1
    vstelm.d VT0, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VT0, XX, 0, 1
    add.d XX, XX, INCX
    vstelm.d VT1, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VT1, YY, 0, 1
    add.d YY, YY, INCY
    addi.d  I, I, -1
    blt $r0, I, .L223
    b .L997
    .align 3

.L224: // C==0 S==0
    vstelm.d VXZ, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VXZ, XX, 0, 1
    add.d XX, XX, INCX
    vstelm.d VXZ, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VXZ, XX, 0, 1
    add.d XX, XX, INCX
    vstelm.d VXZ, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VXZ, YY, 0, 1
    add.d YY, YY, INCY
    vstelm.d VXZ, YY, 0, 0
    add.d YY, YY, INCY
    vstelm.d VXZ, YY, 0, 1
    add.d YY, YY, INCY
    vstelm.d VXZ, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VXZ, XX, 0, 1
    add.d XX, XX, INCX
    vstelm.d VXZ, XX, 0, 0
    add.d XX, XX, INCX
    vstelm.d VXZ, XX, 0, 1
    add.d XX, XX, INCX
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