/***************************************************************************
Copyright (c) 2023, The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#define ASSEMBLER

#include "common.h"

#define N    $r4
#define X    $r5
#define INCX $r6
#define I    $r12
#define t1   $f14
#define t2   $f18
#define t3   $f15
#define t4   $f17
#define s1   $f22
#define s2   $f9
#define s3   $f10
#define s4   $f11
#define TEMP $r16
#define a0   $f20
#define a1   $f21
#define x1   $vr9
#define x2   $vr10
#define x3   $vr11
#define x4   $vr12
#define VT0  $vr13
#define VT1  $vr14
#define res0 $vr18
#define neg1 $vr19
#define VX0  $vr20
#define VX1  $vr21
#define VM0  $vr22
#define VM1  $vr23

    PROLOGUE
    vxor.v VM0, VM0, VM0
    vxor.v res0, res0, res0
    bge $r0, N, .L999
    bge $r0, INCX, .L999
    li.d TEMP, 1
    li.w I, -1
    slli.d TEMP, TEMP, ZBASE_SHIFT
    slli.d INCX, INCX, ZBASE_SHIFT
    vreplgr2vr.w neg1, I
    vffint.s.w neg1, neg1
    srai.d I, N, 3
    bne INCX, TEMP, .L20
    bge $r0, I, .L23
    .align 3

.L10:
    vld VX0, X, 0 * SIZE
    vld VX1, X, 4 * SIZE
    addi.d I, I, -1
    vpickev.w x1, VX1, VX0
    vpickod.w x2, VX1, VX0
    vfmul.s x3, neg1, x1
    vfmul.s x4, neg1, x2
    vfcmp.clt.s VT0, x1, res0
    vfcmp.clt.s VT1, x2, res0
    vld VX0, X, 8 * SIZE
    vbitsel.v x1, x1, x3, VT0
    vbitsel.v x2, x2, x4, VT1
    vld VX1, X, 12 * SIZE
    vfadd.s VM1, x1, x2
    vpickev.w x1, VX1, VX0
    vpickod.w x2, VX1, VX0
    vfmul.s x3, neg1, x1
    vfmul.s x4, neg1, x2
    vfcmp.clt.s VT0, x1, res0
    vfcmp.clt.s VT1, x2, res0
    addi.d X, X, 16 * SIZE
    vbitsel.v x1, x1, x3, VT0
    vbitsel.v x2, x2, x4, VT1
    vfadd.s x1, x1, x2
    vfmax.s VM1, x1, VM1
    vfmax.s VM0, VM0, VM1
    blt $r0, I, .L10
    .align 3

.L11:
    vreplvei.w x1, VM0, 0
    vreplvei.w x2, VM0, 1
    vreplvei.w x3, VM0, 2
    vreplvei.w x4, VM0, 3
    vfmax.s VM1, x1, x2
    vfmax.s VM0, x3, x4
    vfmax.s VM0, VM0, VM1
    b .L23
    .align 3

.L20: // INCX!=1
    bge $r0, I, .L23
    .align 3

.L21:
    fld.s t1, X, 0 * SIZE
    fld.s t2, X, 1 * SIZE
    add.d X, X, INCX
    fld.s t3, X, 0 * SIZE
    fld.s t4, X, 1 * SIZE
    add.d X, X, INCX
    fabs.s t1, t1
    fabs.s t2, t2
    fabs.s t3, t3
    fabs.s t4, t4
    fadd.s t1, t1, t2
    fadd.s t3, t3, t4
    fmax.s s1, t1, t3
    fld.s t1, X, 0 * SIZE
    fld.s t2, X, 1 * SIZE
    add.d X, X, INCX
    fld.s t3, X, 0 * SIZE
    fld.s t4, X, 1 * SIZE
    add.d X, X, INCX
    fabs.s t1, t1
    fabs.s t2, t2
    fabs.s t3, t3
    fabs.s t4, t4
    fadd.s t1, t1, t2
    fadd.s t3, t3, t4
    fmax.s s1, t1, t3
    fld.s t1, X, 0 * SIZE
    fld.s t2, X, 1 * SIZE
    add.d X, X, INCX
    fld.s t3, X, 0 * SIZE
    fld.s t4, X, 1 * SIZE
    add.d X, X, INCX
    fabs.s t1, t1
    fabs.s t2, t2
    fabs.s t3, t3
    fabs.s t4, t4
    addi.d I, I, -1
    fadd.s t1, t1, t2
    fadd.s t3, t3, t4
    fmax.s s3, t1, t3
    fld.s t1, X, 0 * SIZE
    fld.s t2, X, 1 * SIZE
    add.d X, X, INCX
    fld.s t3, X, 0 * SIZE
    fld.s t4, X, 1 * SIZE
    add.d X, X, INCX
    fabs.s t1, t1
    fabs.s t2, t2
    fabs.s t3, t3
    fabs.s t4, t4
    fadd.s t1, t1, t2
    fadd.s t3, t3, t4
    fmax.s s4, t1, t3
    blt $r0, I, .L21
    .align 3

.L22:
    fmax.s s1, s1, s2
    fmax.s s3, s3, s4
    fmax.s s1, s1, s3
    .align 3

.L23: //N<8
    andi I, N, 7
    bge $r0, I, .L999
    .align 3

.L24:
    fld.s a0, X, 0 * SIZE
    fld.s a1, X, 1 * SIZE
    addi.d I, I, -1
    fabs.s a0, a0
    fabs.s a1, a1
    fadd.s a0, a0, a1
    add.d  X, X, INCX
    fmax.s s1, a0, s1
    blt $r0, I, .L24
    .align 3

.L999:
    fmov.s $f0, $f22
    jirl $r0, $r1, 0x0
    .align 3

    EPILOGUE
