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
#define x1 $vr9
#define x2 $vr10
#define x3 $vr11
#define x4 $vr12
#define VX0 $vr13
#define VX1 $vr14
#define VM0 $vr15
#define VM1 $vr16
#define VINC2 $vr17
#define VINC4 $vr18
#define VI0 $vr20
#define VI1 $vr21
#define VI2 $vr22
#define VI3 $vr8
#define VI4 $vr19
#define VT0 $vr23

    PROLOGUE
    li.d  i0, 0
    bge $r0, N, .L999
    bge $r0, INCX, .L999
    li.d TEMP, 1
    slli.d TEMP, TEMP, BASE_SHIFT
    slli.d INCX, INCX, BASE_SHIFT
    bne INCX, TEMP, .L20
    vld VM0, X, 0
    addi.d i0, i0, 1
    srai.d I, N, 3
    bge $r0, I, .L11
    slli.d i0, i0, 1 //2
    vreplgr2vr.d VINC2, i0
    slli.d i0, i0, 1 //4
    vreplgr2vr.d VINC4, i0
    addi.d i0, i0, -7
    vinsgr2vr.d VI1, i0, 0 //initialize the index value for vectorization
    addi.d i0, i0, 1
    vinsgr2vr.d VI1, i0, 1
    addi.d i0, i0, 3
    vinsgr2vr.d VI0, i0, 0 //1
    addi.d i0, i0, 1
    vinsgr2vr.d VI0, i0, 1 //2
    .align 3

.L10:
    vld VX0, X, 0 * SIZE
    vadd.d VI1, VI1, VINC4
    vld VX1, X, 2 * SIZE
    vadd.d VI2, VI1, VINC2
    vfmaxa.d x1, VX0, VX1
    vfcmp.ceq.d VT0, VX0, x1
    vbitsel.v x2, VI2, VI1, VT0
    vld VX0, X, 4 * SIZE
    vadd.d VI1, VI2, VINC2
    vld VX1, X, 6 * SIZE
    vadd.d VI2, VI1, VINC2
    vfmaxa.d x3, VX0, VX1
    vfcmp.ceq.d VT0, VX0, x3
    vbitsel.v x4, VI2, VI1, VT0
    vfmaxa.d x3, x1, x3
    vfcmp.ceq.d VT0, x1, x3
    vbitsel.v x2, x4, x2, VT0
    vfmaxa.d VM1, VM0, x3
    vfcmp.ceq.d VT0, VM0, VM1
    vbitsel.v VM0, VM1, VM0, VT0
    vbitsel.v VI0, x2, VI0, VT0
    addi.d I, I, -1
    addi.d X, X, 8 * SIZE
    blt $r0, I, .L10
    .align 3

.L15:
    vreplvei.d VI1, VI0, 0
    vreplvei.d VI2, VI0, 1
    vreplvei.d x1, VM0, 0
    vreplvei.d x2, VM0, 1
    li.d TEMP, 1 //处理尾数相等时取最小序号
    movgr2fr.d $f17, TEMP
    ffint.d.l $f17, $f17
    vfcmp.ceq.d VT0, x2, x1
    fcmp.ceq.d $fcc0, $f23, $f17
    bceqz $fcc0, .L16
    vfcmp.clt.d VT0, VI1, VI0
    vbitsel.v VI0, VI0, VI1, VT0
    b .L17
    .align 3

.L16:
    vfmaxa.d VM0, x1, x2
    vfcmp.ceq.d VT0, x1, VM0
    vbitsel.v VI0, VI2, VI1, VT0
    .align 3

.L17:
    movfr2gr.d i0, $f20
    .align 3

.L11: //INCX==1 and N<8
    andi I, N, 7
    bge $r0, I, .L14
    srai.d i1, N, 3
    slli.d i1, i1, 3
    addi.d i1, i1, 1 //current index
    movgr2fr.d $f21, i1
    movgr2fr.d $f20, i0
    .align 3

.L13:
    fld.d $f9, X, 0
    vfmaxa.d VM1, x1, VM0
    vfcmp.ceq.d VT0, VM0, VM1
    vbitsel.v VM0, VM1, VM0, VT0
    vbitsel.v VI0, VI1, VI0, VT0
    addi.d I, I, -1
    addi.d i1, i1, 1
    addi.d  X, X, SIZE
    movgr2fr.d $f21, i1
    blt $r0, I, .L13
    movfr2gr.d i0, $f20
    .align 3

.L14:
    move $r4, $r17
    jirl $r0, $r1, 0x0
    .align 3

.L20: // INCX!=1
    move TEMP, X
    addi.d i0, i0, 1
    ld.d t1, TEMP, 0 * SIZE
    add.d TEMP, TEMP, INCX
    vinsgr2vr.d VM0, t1, 0
    srai.d I, N, 3
    bge $r0, I, .L21
    ld.d t2, TEMP, 0 * SIZE
    add.d TEMP, TEMP, INCX
    vinsgr2vr.d VM0, t2, 1
    slli.d i0, i0, 1 //2
    vreplgr2vr.d VINC2, i0
    slli.d i0, i0, 1 //4
    vreplgr2vr.d VINC4, i0
    addi.d i0, i0, -7
    vinsgr2vr.d VI1, i0, 0 //initialize the index value for vectorization
    addi.d i0, i0, 1
    vinsgr2vr.d VI1, i0, 1
    addi.d i0, i0, 3
    vinsgr2vr.d VI0, i0, 0 //1
    addi.d i0, i0, 1
    vinsgr2vr.d VI0, i0, 1 //2
    .align 3

.L24:
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX0, t1, 0
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX0, t2, 1
    vadd.d VI1, VI1, VINC4
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX1, t1, 0
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX1, t2, 1
    vadd.d VI2, VI1, VINC2
    vfmaxa.d x1, VX0, VX1
    vfcmp.ceq.d VT0, VX0, x1
    vbitsel.v x2, VI2, VI1, VT0
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX0, t1, 0
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX0, t2, 1
    vadd.d VI1, VI2, VINC2
    ld.d t1, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX1, t1, 0
    ld.d t2, X, 0 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d VX1, t2, 1
    vadd.d VI2, VI1, VINC2
    vfmaxa.d x3, VX0, VX1
    vfcmp.ceq.d VT0, VX0, x3
    vbitsel.v x4, VI2, VI1, VT0
    vfmaxa.d x3, x1, x3
    vfcmp.ceq.d VT0, x1, x3
    vbitsel.v x2, x4, x2, VT0
    vfmaxa.d VM1, VM0, x3
    vbitsel.v VM0, VM1, VM0, VT0
    vfcmp.ceq.d VT0, VM0, VM1
    vbitsel.v VI0, x2, VI0, VT0
    addi.d I, I, -1
    blt $r0, I, .L24
    .align 3

.L25:
    vreplvei.d VI1, VI0, 0
    vreplvei.d VI2, VI0, 1
    vreplvei.d x1, VM0, 0
    vreplvei.d x2, VM0, 1
    li.d TEMP, 1 //处理尾数相等时取最小序号
    movgr2fr.d $f17, TEMP
    ffint.d.l $f17, $f17
    vfcmp.ceq.d VT0, x2, x1
    fcmp.ceq.d $fcc0, $f23, $f17
    bceqz $fcc0, .L26
    vfcmp.clt.d VT0, VI1, VI0
    vbitsel.v VI0, VI0, VI1, VT0
    b .L27
    .align 3

.L26:
    vfmaxa.d VM0, x1, x2
    vfcmp.ceq.d VT0, x1, VM0
    vbitsel.v VI0, VI2, VI1, VT0
    .align 3

.L27:
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
    vfmaxa.d VM1, x1, VM0
    vfcmp.ceq.d VT0, VM0, VM1
    vbitsel.v VM0, VM1, VM0, VT0
    vbitsel.v VI0, VI1, VI0, VT0
    addi.d I, I, -1
    addi.d i1, i1, 1
    add.d  X, X, INCX
    movgr2fr.d $f21, i1
    blt $r0, I, .L22
    movfr2gr.d i0, $f20
    .align 3

.L999:
    move $r4, $r17
    jirl $r0, $r1, 0x0
    .align 3

    EPILOGUE
