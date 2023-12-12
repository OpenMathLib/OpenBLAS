#define ASSEMBLER
#include "common.h"

#define N $r4
#define ALPHA $f0
#define X  $r7
#define INCX $r8
#define I $r12
#define TEMP $r13
#define t1 $r14
#define t2 $r18
#define t3 $r15
#define t4 $r17
#define XX $r16
#define VX0 $xr12
#define VX1 $xr13
#define VT0 $xr14
#define VT1 $xr15
#define VALPHA $xr19
#define a1 $f8
#define a2 $f23

    PROLOGUE

    bge $r0, N, .L999
    bge $r0, INCX, .L999
    li.d TEMP, 1
    movgr2fr.d a1, $r0
    ffint.s.l a1, a1
    movgr2fr.d a2, TEMP
    ffint.s.l a2, a2
    slli.d TEMP, TEMP, BASE_SHIFT
    slli.d INCX, INCX, BASE_SHIFT
    fcmp.ceq.s $fcc0, ALPHA, a1
    bcnez $fcc0, .L20   //ALPHA==0
    fcmp.ceq.s $fcc0, ALPHA, a2
    bcnez $fcc0, .L999  //ALPHA==1 return
    srai.d I, N, 3
    beq INCX, TEMP, .L30    //ALPHA！=0|1 and INCX==1
    movfr2gr.s TEMP, ALPHA
    xvreplgr2vr.w VALPHA, TEMP
    move XX, X

.L10:   //ALPHA！=0|1 and INCX!=1
    bge $r0, I, .L32
    .align 3
.L11:
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
    xvfmul.s VT0, VX0, VALPHA
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
    addi.d  I, I, -1
    blt $r0, I, .L11
    b .L32
    .align 3

.L20:
    srai.d I, N, 3
    beq INCX, TEMP, .L24
    bge $r0, I, .L22
    .align 3

.L21:
    fst.s a1, X, 0
    add.d X, X, INCX
    fst.s a1, X, 0
    add.d X, X, INCX
    fst.s a1, X, 0
    add.d X, X, INCX
    fst.s a1, X, 0
    add.d X, X, INCX
    fst.s a1, X, 0
    add.d X, X, INCX
    fst.s a1, X, 0
    add.d X, X, INCX
    fst.s a1, X, 0
    add.d X, X, INCX
    fst.s a1, X, 0
    add.d X, X, INCX
    addi.d  I, I, -1
    blt $r0, I, .L21
    .align 3

.L22:
    andi I, N, 7
    bge $r0, I, .L999
    .align 3
.L23:
    fst.s a1, X, 0 * SIZE
    addi.d I, I, -1
    add.d X, X, INCX
    blt $r0, I, .L23
    jirl $r0, $r1, 0
    .align 3

.L24:
    bge $r0, I, .L26 /*N<8 INCX==1*/
    .align 3
.L25:
    xvxor.v VX0, VX0, VX0
    xvst VX0, X, 0 * SIZE
    addi.d I, I, -1
    addi.d X, X, 8 * SIZE
    blt $r0, I, .L25
    .align 3

.L26:
    andi I, N, 7
    bge $r0, I, .L999
    .align 3
.L27:
    fst.s a1, X, 0 * SIZE
    addi.d I, I, -1
    addi.d X, X, SIZE
    blt $r0, I, .L27
    jirl $r0, $r1, 0
    .align 3

.L30:
    bge $r0, I, .L32/*N<8 INCX==1*/
    movfr2gr.s TEMP, ALPHA
    xvreplgr2vr.w VALPHA , TEMP
    .align 3

.L31:
    xvld VX0, X, 0 * SIZE
    addi.d I, I, -1
    xvfmul.s VT0, VX0, VALPHA
    xvst VT0, X, 0 * SIZE
    addi.d X, X, 8 * SIZE
    blt $r0, I, .L31
    .align 3

.L32:
    andi I, N, 7
    bge $r0, I, .L999
    .align 3
.L33:
    fld.s a1, X, 0 * SIZE
    addi.d I, I, -1
    fmul.s a1, ALPHA, a1
    fst.s a1, X, 0 * SIZE
    add.d X, X, INCX
    blt $r0, I, .L33
    jirl $r0, $r1, 0
    .align 3

.L999:
    jirl    $r0, $r1, 0x0

    EPILOGUE
