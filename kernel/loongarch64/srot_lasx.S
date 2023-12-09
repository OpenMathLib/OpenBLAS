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
    ffint.s.l a1, a1
    slli.d  TEMP, TEMP, BASE_SHIFT
    slli.d  INCX, INCX, BASE_SHIFT
    slli.d  INCY, INCY, BASE_SHIFT
    movfr2gr.s t1, C
    xvreplgr2vr.w VXC, t1
    movfr2gr.s t2, S
    xvreplgr2vr.w VXS, t2
    movfr2gr.s t3, a1
    xvreplgr2vr.w VXZ, t3
    srai.d I, N, 3
    bne INCX, TEMP, .L20
    bne INCY, TEMP, .L12 // INCX==1 and INCY!=1
    b .L11  // INCX==1 and INCY==1
.L20:
    bne INCY, TEMP, .L22 // INCX!=1 and INCY!=1
    b .L21 // INCX!=1 and INCY==1

.L11:
    bge $r0, I, .L997
    fcmp.ceq.s $fcc0, C, a1
    bcnez $fcc0, .L110
    fcmp.ceq.s $fcc0, S, a1
    bcnez $fcc0, .L112 // C!=0 S==0
    b .L111 // C!=0 S!=0
    .align 3

.L110:
    fcmp.ceq.s $fcc0, S, a1
    bcnez $fcc0, .L114 // C==0 S==0
    b .L113 // C==0 S!=0
    .align 3

.L111: // C!=0 S!=0
    xvld VX0, X, 0 * SIZE
    xvld VX2, Y, 0 * SIZE
    xvfmul.s VT0, VX0, VXC
    xvfmadd.s VT0, VX2, VXS, VT0
    xvfmul.s VT1, VX0, VXS
    xvfmsub.s VT1, VX2, VXC, VT1
    xvst VT0, X, 0 * SIZE
    xvst VT1, Y, 0 * SIZE
    addi.d X, X, 8 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L111
    b .L997
    .align 3

.L112: // C!=0 S==0
    xvld VX0, X, 0 * SIZE
    xvld VX2, Y, 0 * SIZE
    xvfmul.s VT0, VX0, VXC
    xvfmul.s VT1, VX2, VXC
    xvst VT0, X, 0 * SIZE
    xvst VT1, Y, 0 * SIZE
    addi.d X, X, 8 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L112
    b .L997
    .align 3

.L113: // C==0 S!=0
    xvld VX0, X, 0 * SIZE
    xvld VX2, Y, 0 * SIZE
    xvfmul.s VT0, VX2, VXS
    xvfmul.s VT1, VX0, VXS
    xvfsub.s VT1, VXZ, VT1
    xvst VT0, X, 0 * SIZE
    xvst VT1, Y, 0 * SIZE
    addi.d X, X, 8 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L113
    b .L997
    .align 3

.L114: // C==0 S==0
    xvst VXZ, X, 0 * SIZE
    xvst VXZ, Y, 0 * SIZE
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
    fcmp.ceq.s $fcc0, C, a1
    bcnez $fcc0, .L120
    fcmp.ceq.s $fcc0, S, a1
    bcnez $fcc0, .L122 // C!=0 S==0
    b .L121 // C!=0 S!=0
    .align 3

.L120:
    fcmp.ceq.s $fcc0, S, a1
    bcnez $fcc0, .L124 // C==0 S==0
    b .L123 // C==0 S!=0
    .align 3

.L121: // C!=0 S!=0
    xvld VX0, X, 0 * SIZE
    ld.w t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t4, Y, 0 * SIZE
    add.d Y, Y, INCY
    xvinsgr2vr.w VX2, t1, 0
    xvinsgr2vr.w VX2, t2, 1
    xvinsgr2vr.w VX2, t3, 2
    xvinsgr2vr.w VX2, t4, 3
    ld.w t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t4, Y, 0 * SIZE
    xvinsgr2vr.w VX2, t1, 4
    xvinsgr2vr.w VX2, t2, 5
    xvinsgr2vr.w VX2, t3, 6
    xvinsgr2vr.w VX2, t4, 7
    add.d Y, Y, INCY
    xvfmul.s VT0, VX0, VXC
    xvfmadd.s VT0, VX2, VXS, VT0
    xvfmul.s VT1, VX0, VXS
    xvfmsub.s VT1, VX2, VXC, VT1
    xvst VT0, X, 0 * SIZE
    xvstelm.w VT1, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 3
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 4
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 5
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 6
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 7
    add.d YY, YY, INCY
    addi.d X, X, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L121
    b .L997
    .align 3

.L122: // C!=0 S==0
    xvld VX0, X, 0 * SIZE
    ld.w t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t4, Y, 0 * SIZE
    add.d Y, Y, INCY
    xvinsgr2vr.w VX2, t1, 0
    xvinsgr2vr.w VX2, t2, 1
    xvinsgr2vr.w VX2, t3, 2
    xvinsgr2vr.w VX2, t4, 3
    ld.w t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t4, Y, 0 * SIZE
    xvinsgr2vr.w VX2, t1, 4
    xvinsgr2vr.w VX2, t2, 5
    xvinsgr2vr.w VX2, t3, 6
    xvinsgr2vr.w VX2, t4, 7
    add.d Y, Y, INCY
    xvfmul.s VT0, VX0, VXC
    xvfmul.s VT1, VX2, VXC
    xvst VT0, X, 0 * SIZE
    xvstelm.w VT1, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 3
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 4
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 5
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 6
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 7
    add.d YY, YY, INCY
    addi.d X, X, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L122
    b .L997
    .align 3

.L123: // C==0 S!=0
    xvld VX0, X, 0 * SIZE
    ld.w t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t4, Y, 0 * SIZE
    add.d Y, Y, INCY
    xvinsgr2vr.w VX2, t1, 0
    xvinsgr2vr.w VX2, t2, 1
    xvinsgr2vr.w VX2, t3, 2
    xvinsgr2vr.w VX2, t4, 3
    ld.w t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t4, Y, 0 * SIZE
    xvinsgr2vr.w VX2, t1, 4
    xvinsgr2vr.w VX2, t2, 5
    xvinsgr2vr.w VX2, t3, 6
    xvinsgr2vr.w VX2, t4, 7
    add.d Y, Y, INCY
    xvfmul.s VT0, VX2, VXS
    xvfmul.s VT1, VX0, VXS
    xvfsub.s VT1, VXZ, VT1
    xvst VT0, X, 0 * SIZE
    xvstelm.w VT1, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 3
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 4
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 5
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 6
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 7
    add.d YY, YY, INCY
    addi.d X, X, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L123
    b .L997
    .align 3

.L124: // C==0 S==0
    xvst VXZ, X, 0 * SIZE
    xvstelm.w VXZ, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.w VXZ, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.w VXZ, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.w VXZ, YY, 0, 3
    add.d YY, YY, INCY
    xvstelm.w VXZ, YY, 0, 4
    add.d YY, YY, INCY
    xvstelm.w VXZ, YY, 0, 5
    add.d YY, YY, INCY
    xvstelm.w VXZ, YY, 0, 6
    add.d YY, YY, INCY
    xvstelm.w VXZ, YY, 0, 7
    add.d YY, YY, INCY
    addi.d  I, I, -1
    blt $r0, I, .L124
    b .L997
    .align 3

.L21:// INCX!=1 and INCY==1
    bge $r0, I, .L997
    move XX, X
    fcmp.ceq.s $fcc0, C, a1
    bcnez $fcc0, .L210
    fcmp.ceq.s $fcc0, S, a1
    bcnez $fcc0, .L212 // C!=0 S==0
    b .L211 // C!=0 S!=0
    .align 3

.L210:
    fcmp.ceq.s $fcc0, S, a1
    bcnez $fcc0, .L214 // C==0 S==0
    b .L213 // C==0 S!=0
    .align 3

.L211: // C!=0 S!=0
    xvld VX2, Y, 0 * SIZE
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
    xvinsgr2vr.w VX0, t1, 4
    xvinsgr2vr.w VX0, t2, 5
    xvinsgr2vr.w VX0, t3, 6
    xvinsgr2vr.w VX0, t4, 7
    add.d X, X, INCX
    xvfmul.s VT0, VXC, VX0
    xvfmadd.s VT0, VX2, VXS, VT0
    xvfmul.s VT1, VX0, VXS
    xvfmsub.s VT1, VX2, VXC, VT1
    xvstelm.w VT0, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 3
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 4
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 5
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 6
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 7
    add.d XX, XX, INCX
    xvst VT1, Y, 0 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L211
    b .L997
    .align 3

.L212: // C!=0 S==0
    xvld VX2, Y, 0 * SIZE
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
    xvinsgr2vr.w VX0, t1, 4
    xvinsgr2vr.w VX0, t2, 5
    xvinsgr2vr.w VX0, t3, 6
    xvinsgr2vr.w VX0, t4, 7
    add.d X, X, INCX
    xvfmul.s VT0, VXC, VX0
    xvfmul.s VT1, VX2, VXC
    xvstelm.w VT0, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 3
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 4
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 5
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 6
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 7
    add.d XX, XX, INCX
    xvst VT1, Y, 0 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L212
    b .L997
    .align 3

.L213: // C==0 S!=0
    xvld VX2, Y, 0 * SIZE
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
    xvinsgr2vr.w VX0, t1, 4
    xvinsgr2vr.w VX0, t2, 5
    xvinsgr2vr.w VX0, t3, 6
    xvinsgr2vr.w VX0, t4, 7
    add.d X, X, INCX
    xvfmul.s VT0, VXS, VX2
    xvfmul.s VT1, VXS, VX0
    xvfsub.s VT1, VXZ, VT1
    xvstelm.w VT0, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 3
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 4
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 5
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 6
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 7
    add.d XX, XX, INCX
    xvst VT1, Y, 0 * SIZE
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L213
    b .L997
    .align 3

.L214: // C==0 S==0
    xvstelm.w VXZ, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.w VXZ, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.w VXZ, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.w VXZ, XX, 0, 3
    add.d XX, XX, INCX
    xvst VT1, Y, 0 * SIZE
    xvstelm.w VXZ, XX, 0, 4
    add.d XX, XX, INCX
    xvstelm.w VXZ, XX, 0, 5
    add.d XX, XX, INCX
    xvstelm.w VXZ, XX, 0, 6
    add.d XX, XX, INCX
    xvstelm.w VXZ, XX, 0, 7
    add.d XX, XX, INCX
    addi.d Y, Y, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L211
    b .L997
    .align 3

.L22:
    bge $r0, I, .L997
    move YY, Y
    move XX, X
    fcmp.ceq.s $fcc0, C, a1
    bcnez $fcc0, .L220
    fcmp.ceq.s $fcc0, S, a1
    bcnez $fcc0, .L222 // C!=0 S==0
    b .L221 // C!=0 S!=0
    .align 3

.L220:
    fcmp.ceq.s $fcc0, S, a1
    bcnez $fcc0, .L224 // C==0 S==0
    b .L223 // C==0 S!=0
    .align 3

.L221: // C!=0 S!=0
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
    ld.w t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t4, Y, 0 * SIZE
    add.d Y, Y, INCY
    xvinsgr2vr.w VX2, t1, 0
    xvinsgr2vr.w VX2, t2, 1
    xvinsgr2vr.w VX2, t3, 2
    xvinsgr2vr.w VX2, t4, 3
    ld.w t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t4, Y, 0 * SIZE
    xvinsgr2vr.w VX2, t1, 4
    xvinsgr2vr.w VX2, t2, 5
    xvinsgr2vr.w VX2, t3, 6
    xvinsgr2vr.w VX2, t4, 7
    add.d Y, Y, INCY
    xvfmul.s  VT0, VX0, VXC
    xvfmadd.s VT0, VX2, VXS, VT0
    xvfmul.s  VT1, VX0, VXS
    xvfmsub.s VT1, VX2, VXC, VT1
    xvstelm.w VT0, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 3
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 4
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 5
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 6
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 7
    add.d XX, XX, INCX
    xvstelm.w VT1, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 3
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 4
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 5
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 6
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 7
    add.d YY, YY, INCY
    addi.d  I, I, -1
    blt $r0, I, .L221
    b .L997
    .align 3

.L222: // C!=0 S==0
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
    ld.w t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t4, Y, 0 * SIZE
    add.d Y, Y, INCY
    xvinsgr2vr.w VX2, t1, 0
    xvinsgr2vr.w VX2, t2, 1
    xvinsgr2vr.w VX2, t3, 2
    xvinsgr2vr.w VX2, t4, 3
    ld.w t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t4, Y, 0 * SIZE
    xvinsgr2vr.w VX2, t1, 4
    xvinsgr2vr.w VX2, t2, 5
    xvinsgr2vr.w VX2, t3, 6
    xvinsgr2vr.w VX2, t4, 7
    add.d Y, Y, INCY
    xvfmul.s VT0, VX0, VXC
    xvfmul.s VT1, VX2, VXC
    xvstelm.w VT0, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 3
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 4
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 5
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 6
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 7
    add.d XX, XX, INCX
    xvstelm.w VT1, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 3
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 4
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 5
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 6
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 7
    add.d YY, YY, INCY
    addi.d  I, I, -1
    blt $r0, I, .L222
    b .L997
    .align 3

.L223: // C==0 S!=0
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
    ld.w t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t4, Y, 0 * SIZE
    add.d Y, Y, INCY
    xvinsgr2vr.w VX2, t1, 0
    xvinsgr2vr.w VX2, t2, 1
    xvinsgr2vr.w VX2, t3, 2
    xvinsgr2vr.w VX2, t4, 3
    ld.w t1, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t2, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t3, Y, 0 * SIZE
    add.d Y, Y, INCY
    ld.w t4, Y, 0 * SIZE
    xvinsgr2vr.w VX2, t1, 4
    xvinsgr2vr.w VX2, t2, 5
    xvinsgr2vr.w VX2, t3, 6
    xvinsgr2vr.w VX2, t4, 7
    add.d Y, Y, INCY
    xvfmul.s VT0, VX2, VXS
    xvfmul.s VT1, VX0, VXS
    xvfsub.s VT1, VXZ, VT1
    xvstelm.w VT0, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 3
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 4
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 5
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 6
    add.d XX, XX, INCX
    xvstelm.w VT0, XX, 0, 7
    add.d XX, XX, INCX
    xvstelm.w VT1, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 3
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 4
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 5
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 6
    add.d YY, YY, INCY
    xvstelm.w VT1, YY, 0, 7
    add.d YY, YY, INCY
    addi.d  I, I, -1
    blt $r0, I, .L223
    b .L997
    .align 3

.L224: // C==0 S==0
    xvstelm.w VXZ, XX, 0, 0
    add.d XX, XX, INCX
    xvstelm.w VXZ, XX, 0, 1
    add.d XX, XX, INCX
    xvstelm.w VXZ, XX, 0, 2
    add.d XX, XX, INCX
    xvstelm.w VXZ, XX, 0, 3
    add.d XX, XX, INCX
    xvstelm.w VXZ, YY, 0, 0
    add.d YY, YY, INCY
    xvstelm.w VXZ, YY, 0, 1
    add.d YY, YY, INCY
    xvstelm.w VXZ, YY, 0, 2
    add.d YY, YY, INCY
    xvstelm.w VXZ, YY, 0, 3
    add.d YY, YY, INCY
    xvstelm.w VXZ, XX, 0, 4
    add.d XX, XX, INCX
    xvstelm.w VXZ, XX, 0, 5
    add.d XX, XX, INCX
    xvstelm.w VXZ, XX, 0, 6
    add.d XX, XX, INCX
    xvstelm.w VXZ, XX, 0, 7
    add.d XX, XX, INCX
    xvstelm.w VXZ, YY, 0, 4
    add.d YY, YY, INCY
    xvstelm.w VXZ, YY, 0, 5
    add.d YY, YY, INCY
    xvstelm.w VXZ, YY, 0, 6
    add.d YY, YY, INCY
    xvstelm.w VXZ, YY, 0, 7
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
    fld.s $f12, X, 0 * SIZE
    fld.s $f13, Y, 0 * SIZE
    fmul.s $f10, $f12, C
    fmadd.s $f10, $f13, S, $f10
    fst.s $f10, X, 0 * SIZE
    addi.d I, I, -1
    fmul.s $f20, $f12, S
    fmsub.s $f20, $f13, C, $f20
    fst.s $f20, Y, 0 * SIZE
    add.d  X, X, INCX
    add.d  Y, Y, INCY
    blt $r0, I, .L998
    .align 3

.L999:
    move $r4, $r12
    jirl $r0, $r1, 0x0
    .align 3

    EPILOGUE