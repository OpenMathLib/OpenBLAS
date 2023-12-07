#define ASSEMBLER

#include "common.h"

#define N $r4
#define X $r5
#define INCX $r6
#define I $r12
#define t1 $r14
#define t2 $r18
#define t3 $r15
#define t4 $r17
#define TEMP $r16
#define x1 $vr9
#define x2 $vr10
#define VX0 $vr20
#define VX1 $vr21
#define VM0 $vr22
#define VM1 $vr23
#define VM2 $vr18
#define VM3 $vr19

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
    vld VX1, X, 2 * SIZE
    addi.d I, I, -1
    vfmaxa.d VM1, VX1, VX0
    vld VX0, X, 4 * SIZE
    vld VX1, X, 6 * SIZE
    vfmaxa.d VM2, VX1, VX0
    vfmaxa.d VM3, VM1, VM2
    addi.d X, X, 8 * SIZE
    vfmaxa.d VM0, VM0, VM3
    blt $r0, I, .L10
    .align 3

.L11:
    vreplvei.d x1, VM0, 0
    vreplvei.d x2, VM0, 1
    vfmaxa.d VM0, x1, x2
    .align 3

.L12: //INCX==1 and N<8
    andi I, N, 7
    bge $r0, I, .L14
    .align 3

.L13:
    vld x1, X, 0
    addi.d I, I, -1
    vfmaxa.d VM0, VM0, x1
    addi.d  X, X, SIZE
    blt $r0, I, .L13
    .align 3

.L14:
    fabs.d $f22, $f22
    fmov.d $f0, $f22
    jirl $r0, $r1, 0x0
    .align 3

.L20: // INCX!=1
    move TEMP, X // initialize the maxa value
    ld.d t1, TEMP, 0 * SIZE
    add.d TEMP, TEMP, INCX
    vinsgr2vr.d VM0, t1, 0
    srai.d I, N, 3
    bge $r0, I, .L23
    ld.d t2, TEMP, 0 * SIZE
    add.d TEMP, TEMP, INCX
    vinsgr2vr.d VM0, t2, 1
    .align 3

.L21:
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX0, t1, 0
    vinsgr2vr.d VX0, t2, 1
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX1, t3, 0
    vinsgr2vr.d VX1, t4, 1
    vfmaxa.d VM1, VX0, VX1
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX0, t1, 0
    vinsgr2vr.d VX0, t2, 1
    ld.d t3, X, 0 * SIZE
    add.d X, X, INCX
    ld.d t4, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX1, t3, 0
    vinsgr2vr.d VX1, t4, 1
    vfmaxa.d VM2, VX0, VX1
    vfmaxa.d VM3, VM1, VM2
    vfmaxa.d VM0, VM0, VM3
    addi.d I, I, -1
    blt $r0, I, .L21
    .align 3

.L22:
    vreplvei.d x1, VM0, 0
    vreplvei.d x2, VM0, 1
    vfmaxa.d VM0, x1, x2
    .align 3

.L23: //INCX!=1 and N<8
    andi I, N, 7
    bge $r0, I, .L999
    .align 3

.L24:
    vld x1, X, 0
    addi.d I, I, -1
    vfmaxa.d VM0, VM0, x1
    add.d  X, X, INCX
    blt $r0, I, .L24
    .align 3

.L999:
    fabs.d $f22, $f22
    fmov.d $f0, $f22
    jirl $r0, $r1, 0x0
    .align 3

    EPILOGUE
