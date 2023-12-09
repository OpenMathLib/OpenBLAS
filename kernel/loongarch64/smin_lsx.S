#define ASSEMBLER

#include "common.h"

#define N $r4
#define X $r5
#define INCX $r6
#define I $r12
#define J $r13
#define t1 $r14
#define t2 $r18
#define t3 $r15
#define t4 $r17
#define TEMP $r16
#define m0 $vr8
#define x1 $vr9
#define x2 $vr10
#define x3 $vr11
#define x4 $vr12
#define VX0 $vr20
#define VX1 $vr21
#define VM0 $vr22
#define VM1 $vr23

    PROLOGUE

    bge $r0, N, .L999
    bge $r0, INCX, .L999
    li.d TEMP, 1
    slli.d TEMP, TEMP, BASE_SHIFT
    slli.d INCX, INCX, BASE_SHIFT
    bne INCX, TEMP, .L20
    vld VM0, X, 0
    srai.d I, N, 3
    bge $r0, I, .L12
    .align 3

.L10:
    vld VX0, X, 0 * SIZE
    vld VX1, X, 4 * SIZE
    addi.d I, I, -1
    vfmin.s VM1, VX0, VX1
    addi.d X, X, 8 * SIZE
    vfmin.s VM0, VM0, VM1
    blt $r0, I, .L10
    .align 3

.L11:
    vreplvei.w x1, VM0, 0
    vreplvei.w x2, VM0, 1
    vreplvei.w x3, VM0, 2
    vreplvei.w x4, VM0, 3
    vfmin.s VM1, x1, x2
    vfmin.s VM0, x3, x4
    vfmin.s VM0, VM0, VM1
    .align 3

.L12: //INCX==1 and N<8
    andi I, N, 7
    li.d J, 4
    bge J, I, .L13 // 4<N<8
    vld VX0, X, 0
    slli.d J, J, 1 // 8
    sub.d I, J, I
    slli.d I, I, BASE_SHIFT
    vldx VX1, X, I
    vfmin.s m0, VX0, VX1 //patial repeat read
    vreplvei.w x1, m0, 0
    vreplvei.w x2, m0, 1
    vreplvei.w x3, m0, 2
    vreplvei.w x4, m0, 3
    vfmin.s m0, x1, x2
    vfmin.s VM1, x3, x4
    vfmin.s m0, m0, VM1
    vfmin.s VM0, m0, VM0
    fmov.s $f0, $f22
    jirl $r0, $r1, 0x0
    .align 3

.L13: //INCX==1 and 0<=N<=4
    bge $r0, I, .L15
    .align 3

.L14:
    vld x1, X, 0
    addi.d I, I, -1
    vfmin.s VM0, VM0, x1
    addi.d  X, X, SIZE
    blt $r0, I, .L14
    .align 3

.L15:
    fmov.s $f0, $f22
    jirl $r0, $r1, 0x0
    .align 3

.L20: // INCX!=1
    move TEMP, X
    ld.w t1, TEMP, 0 * SIZE
    add.d TEMP, TEMP, INCX
    vinsgr2vr.w VM0, t1, 0
    srai.d I, N, 3
    bge $r0, I, .L23
    ld.w t2, TEMP, 0 * SIZE
    add.d TEMP, TEMP, INCX
    ld.w t3, TEMP, 0 * SIZE
    add.d TEMP, TEMP, INCX
    ld.w t4, TEMP, 0 * SIZE
    add.d TEMP, TEMP, INCX
    vinsgr2vr.w VM0, t2, 1
    vinsgr2vr.w VM0, t3, 2
    vinsgr2vr.w VM0, t4, 3
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
    vinsgr2vr.w VX0, t1, 0
    vinsgr2vr.w VX0, t2, 1
    vinsgr2vr.w VX0, t3, 2
    vinsgr2vr.w VX0, t4, 3
    ld.w t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.w t2, X, 0 * SIZE
    add.d X, X, INCX
    ld.w t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.w t4, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.w VX1, t1, 0
    vinsgr2vr.w VX1, t2, 1
    vinsgr2vr.w VX1, t3, 2
    vinsgr2vr.w VX1, t4, 3
    vfmin.s VM1, VX0, VX1
    addi.d I, I, -1
    vfmin.s VM0, VM0, VM1
    blt $r0, I, .L21
    .align 3

.L22:
    vreplvei.w x1, VM0, 0
    vreplvei.w x2, VM0, 1
    vreplvei.w x3, VM0, 2
    vreplvei.w x4, VM0, 3
    vfmin.s VM1, x1, x2
    vfmin.s VM0, x3, x4
    vfmin.s VM0, VM0, VM1
    .align 3

.L23: //INCX!=1 and N<8
    andi I, N, 7
    bge $r0, I, .L999
    .align 3

.L24:
    vld x1, X, 0
    vfmin.s VM0, VM0, x1
    addi.d I, I, -1
    add.d  X, X, INCX
    blt $r0, I, .L24
    .align 3

.L999:
    fmov.s $f0, $f22
    jirl $r0, $r1, 0x0
    .align 3

    EPILOGUE
