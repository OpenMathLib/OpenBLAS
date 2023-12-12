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
#define m0 $xr8
#define x1 $xr9
#define x2 $xr10
#define x3 $xr11
#define x4 $xr12
#define VX0 $xr20
#define VX1 $xr21
#define VM0 $xr22
#define VM1 $xr23
#define VM2 $xr19

    PROLOGUE

    bge $r0, N, .L999
    bge $r0, INCX, .L999
    li.d TEMP, 1
    slli.d TEMP, TEMP, BASE_SHIFT
    slli.d INCX, INCX, BASE_SHIFT
    bne INCX, TEMP, .L20
    xvld VM0, X, 0
    srai.d I, N, 3
    bge $r0, I, .L12
    .align 3

.L10:
    xvld VX0, X, 0 * SIZE
    addi.d I, I, -1
    xvld VX1, X, 4 * SIZE
    xvfmina.d VM1, VX1, VX0
    addi.d X, X, 8 * SIZE
    xvfmina.d VM0, VM0, VM1
    blt $r0, I, .L10
    .align 3

.L11:
    xvpickve.d x1, VM0, 0
    xvpickve.d x2, VM0, 1
    xvpickve.d x3, VM0, 2
    xvpickve.d x4, VM0, 3
    xvfmina.d VM1, x1, x2
    xvfmina.d VM2, x3, x4
    xvfmina.d VM0, VM1, VM2
    .align 3

.L12: //INCX==1 and N<8
    andi I, N, 7
    li.d J, 4
    bge J, I, .L13 // 4<N<8
    xvld VX0, X, 0
    slli.d J, J, 1 // 8
    sub.d I, J, I
    slli.d I, I, BASE_SHIFT
    xvldx VX1, X, I
    xvfmina.d m0, VX0, VX1 //patial repeat read
    xvpickve.d x1, m0, 0
    xvpickve.d x2, m0, 1
    xvpickve.d x3, m0, 2
    xvpickve.d x4, m0, 3
    xvfmina.d VM1, x1, x2
    xvfmina.d m0, x3, x4
    xvfmina.d m0, m0, VM1
    xvfmina.d VM0, m0, VM0
    fabs.d $f22, $f22
    fmov.d $f0, $f22
    jirl $r0, $r1, 0x0
    .align 3

.L13: //INCX==1 and 0<=N<=4
    bge $r0, I, .L15
    .align 3

.L14:
    xvld x1, X, 0
    addi.d I, I, -1
    xvfmina.d VM0, VM0, x1
    addi.d  X, X, SIZE
    blt $r0, I, .L14
    .align 3

.L15:
    fabs.d $f22, $f22
    fmov.d $f0, $f22
    jirl $r0, $r1, 0x0
    .align 3

.L20: // INCX!=1
    move TEMP, X // initialize the mina value
    ld.d t1, TEMP, 0 * SIZE
    add.d TEMP, TEMP, INCX
    xvinsgr2vr.d VM0, t1, 0
    srai.d I, N, 3
    bge $r0, I, .L23
    ld.d t2, TEMP, 0 * SIZE
    add.d TEMP, TEMP, INCX
    ld.d t3, TEMP, 0 * SIZE
    add.d TEMP, TEMP, INCX
    ld.d t4, TEMP, 0 * SIZE
    add.d TEMP, TEMP, INCX
    xvinsgr2vr.d VM0, t2, 1
    xvinsgr2vr.d VM0, t3, 2
    xvinsgr2vr.d VM0, t4, 3
    .align 3

.L21:
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
    addi.d I, I, -1
    xvfmina.d VM1, VX1, VX0
    xvfmina.d VM0, VM1, VM0
    blt $r0, I, .L21
    .align 3

.L22:
    xvpickve.d x1, VM0, 0
    xvpickve.d x2, VM0, 1
    xvpickve.d x3, VM0, 2
    xvpickve.d x4, VM0, 3
    xvfmina.d VM1, x1, x2
    xvfmina.d VM2, x3, x4
    xvfmina.d VM0, VM1, VM2
    .align 3

.L23: //INCX!=1 and N<8
    andi I, N, 7
    bge $r0, I, .L999
    .align 3

.L24:
    xvld x1, X, 0
    addi.d I, I, -1
    xvfmina.d VM0, VM0, x1
    add.d  X, X, INCX
    blt $r0, I, .L24
    .align 3

.L999:
    fabs.d $f22, $f22
    fmov.d $f0, $f22
    jirl $r0, $r1, 0x0
    .align 3

    EPILOGUE