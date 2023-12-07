#define ASSEMBLER

#include "common.h"

#define N $r4
#define X $r5
#define INCX $r6
#define I $r12
#define t1 $r13
#define t2 $r15
#define t3 $r18
#define t4 $r16
#define i0 $r17
#define i1 $r14
#define TEMP $r19
#define x1 $xr9
#define x2 $xr10
#define x3 $xr11
#define x4 $xr12
#define VX0 $xr13
#define VX1 $xr14
#define VM0 $xr15
#define VM1 $xr16
#define VINC4 $xr17
#define VINC8 $xr18
#define VI0 $xr20
#define VI1 $xr21
#define VI2 $xr22
#define VI3 $xr8
#define VI4 $xr19
#define VT0 $xr23

    PROLOGUE
    li.d  i0, 0
    bge $r0, N, .L999
    bge $r0, INCX, .L999
    li.d TEMP, 1
    slli.d TEMP, TEMP, BASE_SHIFT
    slli.d INCX, INCX, BASE_SHIFT
    bne INCX, TEMP, .L20
    xvld VM0, X, 0
    addi.d i0, i0, 1
    srai.d I, N, 3
    bge $r0, I, .L21
    slli.d i0, i0, 2 //4
    xvreplgr2vr.d VINC4, i0
    slli.d i0, i0, 1 //8
    xvreplgr2vr.d VINC8, i0
    addi.d i0, i0, -15
    xvinsgr2vr.d VI1, i0, 0 //initialize the index value for vectorization
    addi.d i0, i0, 1
    xvinsgr2vr.d VI1, i0, 1
    addi.d i0, i0, 1
    xvinsgr2vr.d VI1, i0, 2
    addi.d i0, i0, 1
    xvinsgr2vr.d VI1, i0, 3
    addi.d i0, i0, 5
    xvinsgr2vr.d VI0, i0, 0 //1
    addi.d i0, i0, 1
    xvinsgr2vr.d VI0, i0, 1 //2
    addi.d i0, i0, 1
    xvinsgr2vr.d VI0, i0, 2 //3
    addi.d i0, i0, 1
    xvinsgr2vr.d VI0, i0, 3 //4
    .align 3

.L10:
    xvld VX0, X, 0 * SIZE
    xvadd.d VI1, VI1, VINC8
    xvld VX1, X, 4 * SIZE
    xvadd.d VI2, VI1, VINC4
    xvfmina.d VM1, VX0, VX1
    xvfcmp.ceq.d VT0, VX0, VM1
    addi.d I, I, -1
    xvbitsel.v VI2, VI2, VI1, VT0
    xvfmina.d VM1, VM0, VM1
    xvfcmp.ceq.d VT0, VM0, VM1
    addi.d X, X, 8 * SIZE
    xvbitsel.v VM0, VM1, VM0, VT0
    xvbitsel.v VI0, VI2, VI0, VT0
    blt $r0, I, .L10
    .align 3

.L15:
    xvpickve.d VI1, VI0, 0
    xvpickve.d VI2, VI0, 1
    xvpickve.d VI3, VI0, 2
    xvpickve.d VI4, VI0, 3
    xvpickve.d x1, VM0, 0
    xvpickve.d x2, VM0, 1
    xvpickve.d x3, VM0, 2
    xvpickve.d x4, VM0, 3
    xvfmina.d VM1, x1, x2
    xvfcmp.ceq.d VT0, x1, VM1
    xvbitsel.v VINC4, VI2, VI1, VT0
    xvfmina.d VM0, x4, x3
    xvfcmp.ceq.d VT0, x3, VM0
    xvbitsel.v VINC8, VI4, VI3, VT0
    xvfmina.d VM0, VM0, VM1
    xvfcmp.ceq.d VT0, VM0, VM1
    xvbitsel.v VI0, VINC8, VINC4, VT0
    li.d TEMP, 1 //处理尾数相等时取最小序号
    movgr2fr.d $f17, TEMP
    ffint.d.l $f17, $f17
    xvfcmp.ceq.d VT0, VM0, x1
    fcmp.ceq.d $fcc0, $f23, $f17
    bceqz $fcc0, .L26
    xvfcmp.clt.d VT0, VI1, VI0
    xvbitsel.v VI0, VI0, VI1, VT0
    b .L26
    .align 3

.L20: // INCX!=1
    move TEMP, X
    addi.d i0, i0, 1
    ld.d t1, TEMP, 0 * SIZE
    add.d TEMP, TEMP, INCX
    xvinsgr2vr.d VM0, t1, 0
    srai.d I, N, 3
    bge $r0, I, .L21
    ld.d t2, TEMP, 0 * SIZE
    add.d TEMP, TEMP, INCX
    ld.d t3, TEMP, 0 * SIZE
    add.d TEMP, TEMP, INCX
    ld.d t4, TEMP, 0 * SIZE
    add.d TEMP, TEMP, INCX
    xvinsgr2vr.d VM0, t2, 1
    xvinsgr2vr.d VM0, t3, 2
    xvinsgr2vr.d VM0, t4, 3
    slli.d i0, i0, 2 //4
    xvreplgr2vr.d VINC4, i0
    slli.d i0, i0, 1 //8
    xvreplgr2vr.d VINC8, i0
    addi.d i0, i0, -15
    xvinsgr2vr.d VI1, i0, 0 //initialize the index value for vectorization
    addi.d i0, i0, 1
    xvinsgr2vr.d VI1, i0, 1
    addi.d i0, i0, 1
    xvinsgr2vr.d VI1, i0, 2
    addi.d i0, i0, 1
    xvinsgr2vr.d VI1, i0, 3
    addi.d i0, i0, 5
    xvinsgr2vr.d VI0, i0, 0 //1
    addi.d i0, i0, 1
    xvinsgr2vr.d VI0, i0, 1 //2
    addi.d i0, i0, 1
    xvinsgr2vr.d VI0, i0, 2 //3
    addi.d i0, i0, 1
    xvinsgr2vr.d VI0, i0, 3 //4
    .align 3

.L24:
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
    xvadd.d VI1, VI1, VINC8
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
    xvadd.d VI2, VI1, VINC4
    xvfmina.d VM1, VX0, VX1
    xvfcmp.ceq.d VT0, VX0, VM1
    xvbitsel.v VI2, VI2, VI1, VT0
    xvfmina.d VM1, VM0, VM1
    xvfcmp.ceq.d VT0, VM0, VM1
    addi.d I, I, -1
    xvbitsel.v VM0, VM1, VM0, VT0
    xvbitsel.v VI0, VI2, VI0, VT0
    blt $r0, I, .L24
    .align 3

.L25:
    xvpickve.d VI1, VI0, 0
    xvpickve.d VI2, VI0, 1
    xvpickve.d VI3, VI0, 2
    xvpickve.d VI4, VI0, 3
    xvpickve.d x1, VM0, 0
    xvpickve.d x2, VM0, 1
    xvpickve.d x3, VM0, 2
    xvpickve.d x4, VM0, 3
    xvfmina.d VM1, x1, x2
    xvfcmp.ceq.d VT0, x1, VM1
    xvbitsel.v VINC4, VI2, VI1, VT0
    xvfmina.d VM0, x4, x3
    xvfcmp.ceq.d VT0, x3, VM0
    xvbitsel.v VINC8, VI4, VI3, VT0
    xvfmina.d VM0, VM0, VM1
    xvfcmp.ceq.d VT0, VM0, VM1
    xvbitsel.v VI0, VINC8, VINC4, VT0
    li.d TEMP, 1 //处理尾数相等时取最小序号
    movgr2fr.d $f17, TEMP
    ffint.d.l $f17, $f17
    xvfcmp.ceq.d VT0, VM0, x1
    fcmp.ceq.d $fcc0, $f23, $f17
    bceqz $fcc0, .L26
    xvfcmp.clt.d VT0, VI1, VI0
    xvbitsel.v VI0, VI0, VI1, VT0
    .align 3

.L26:
    xvfcmp.ceq.d VT0, VM0, x2
    fcmp.ceq.d $fcc0, $f23, $f17
    bceqz $fcc0, .L27
    xvfcmp.clt.d VT0, VI2, VI0
    xvbitsel.v VI0, VI0, VI2, VT0
    .align 3

.L27:
    xvfcmp.ceq.d VT0, VM0, x3
    fcmp.ceq.d $fcc0, $f23, $f17
    bceqz $fcc0, .L28
    xvfcmp.clt.d VT0, VI3, VI0
    xvbitsel.v VI0, VI0, VI3, VT0
    .align 3

.L28:
    xvfcmp.ceq.d VT0, VM0, x4
    fcmp.ceq.d $fcc0, $f23, $f17
    bceqz $fcc0, .L29
    xvfcmp.clt.d VT0, VI4, VI0
    xvbitsel.v VI0, VI0, VI4, VT0
    .align 3

.L29:
    movfr2gr.d i0, $f20
    .align 3

.L21: // N<8
    andi I, N, 7
    bge $r0, I, .L999
    srai.d i1, N, 3
    slli.d i1, i1, 3
    addi.d i1, i1, 1 //current index
    movgr2fr.d $f21, i1
    movgr2fr.d $f20, i0
    .align 3

.L22:
    fld.d $f9, X, 0
    addi.d I, I, -1
    xvfmina.d VM1, x1, VM0
    xvfcmp.ceq.d VT0, VM0, VM1
    add.d  X, X, INCX
    xvbitsel.v VM0, VM1, VM0, VT0
    xvbitsel.v VI0, VI1, VI0, VT0
    addi.d i1, i1, 1
    movgr2fr.d $f21, i1
    blt $r0, I, .L22
    movfr2gr.d i0, $f20
    .align 3

.L999:
    move $r4, $r17
    jirl $r0, $r1, 0x0
    .align 3

    EPILOGUE
