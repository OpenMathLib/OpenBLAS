#define ASSEMBLER
#include "common.h"
#define N      $r4
#define X      $r5
#define INCX   $r6
#define I      $r17
#define TEMP   $r18
#define t1     $r15
#define t2     $r12
#define t3     $r13
#define t4     $r14
#define VX0    $vr12
#define VX1    $vr13
#define VX2    $vr14
#define VX3    $vr15
#define VT0    $vr23
#define VT1    $vr22
#define res1   $vr16
#define res2   $vr17
#define res0   $vr18
#define neg1   $vr19

    PROLOGUE
    vxor.v res1, res1, res1
    vxor.v res2, res2, res2
    vxor.v res0, res0, res0
    bge $r0, N, .L999
    bge $r0, INCX, .L999
    li.w t1, -1
    vreplgr2vr.w neg1, t1
    vffint.s.w neg1, neg1
    li.d  TEMP, SIZE
    slli.d  INCX, INCX, BASE_SHIFT
    srai.d I, N, 3
    bne INCX, TEMP, .L20
    bge $r0, I, .L13
    .align 3

.L11:
    vld VX0, X, 0 * SIZE
    vld VX1, X, 4 * SIZE
    vfmul.s VX2, neg1, VX0
    vfmul.s VX3, neg1, VX1
    vfcmp.clt.s VT0, VX0, res0
    vfcmp.clt.s VT1, VX1, res0
    vbitsel.v VX0, VX0, VX2, VT0
    vbitsel.v VX1, VX1, VX3, VT1
    vfadd.s res2, VX0, VX1
    vfadd.s res1, res1, res2
    addi.d X, X, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L11
    .align 3

.L12:
    vreplvei.w VX1, res1, 1
    vreplvei.w VX2, res1, 2
    vreplvei.w VX3, res1, 3
    vfadd.s res1, VX1, res1
    vfadd.s res1, VX2, res1
    vfadd.s res1, VX3, res1
    .align 3

.L13:
    andi I, N, 7
    bge $r0, I, .L999
    .align 3

.L14:
    fld.s $f12, X, 0 * SIZE
    fabs.s $f12, $f12
    fadd.s $f16, $f12, $f16
    addi.d I, I, -1
    addi.d  X, X, SIZE
    blt $r0, I, .L14
    b .L999
    .align 3

.L20:
    bge $r0, I, .L23
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
    vfmul.s VX2, neg1, VX0
    vfmul.s VX3, neg1, VX1
    vfcmp.clt.s VT0, VX0, res0
    vfcmp.clt.s VT1, VX1, res0
    vbitsel.v VX0, VX0, VX2, VT0
    vbitsel.v VX1, VX1, VX3, VT1
    vfadd.s res2, VX0, VX1
    vfadd.s res1, res1, res2
    addi.d  I, I, -1
    blt $r0, I, .L21
    .align 3

.L22:
    vreplvei.w VX1, res1, 1
    vreplvei.w VX2, res1, 2
    vreplvei.w VX3, res1, 3
    vfadd.s res1, VX1, res1
    vfadd.s res1, VX2, res1
    vfadd.s res1, VX3, res1
    .align 3

.L23:
    andi I, N, 7
    bge $r0, I, .L999
    .align 3

.L24:
    fld.s $f12, X, 0 * SIZE
    fabs.s $f12, $f12
    fadd.s $f16, $f12, $f16
    addi.d I, I, -1
    add.d  X, X, INCX
    blt $r0, I, .L24
    .align 3

.L999:
    fmov.s $f0, $f16
    jirl $r0, $r1, 0x0
    .align 3

    EPILOGUE