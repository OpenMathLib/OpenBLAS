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

#define N      $r4
#define ALPHAR $f0
#define ALPHAI $f1
#define X      $r7
#define INCX   $r8

#define I      $r12
#define TEMP   $r13
#define t1     $r14
#define t2     $r16
#define t3     $r15
#define t4     $r17
#define XX     $r18
#define a1     $f12
#define a2     $f13
#define a3     $f14
#define a4     $f15
#define s1     $f16
#define s2     $f17
#define s3     $f18
#define s4     $f19
#define VX0    $vr8
#define VX1    $vr20
#define VX2    $vr21
#define VX3    $vr22
#define VXAR   $vr23
#define VXAI   $vr19
#define VXZ    $vr12
#define x1     $vr18
#define x2     $vr17
#define x3     $vr16
#define x4     $vr15

    PROLOGUE

    bge $r0, N, .L999
    bge $r0, INCX, .L999
    li.d TEMP, 1
    movgr2fr.d a1, $r0
    FFINT    a1, a1
    slli.d  TEMP, TEMP, ZBASE_SHIFT
    slli.d  INCX, INCX, ZBASE_SHIFT
    MTG  t1, ALPHAR
#ifdef DOUBLE
    vreplgr2vr.d VXAR, t1
    movfr2gr.d t2, ALPHAI
    vreplgr2vr.d VXAI, t2
#else
    vreplgr2vr.w VXAR, t1
    movfr2gr.s t2, ALPHAI
    vreplgr2vr.w VXAI, t2
#endif
    vxor.v VXZ, VXZ, VXZ
    srai.d I, N, 2
    bne INCX, TEMP, .L22

.L11:
    bge $r0, I, .L997
    CMPEQ  $fcc0, ALPHAR, a1
    CMPEQ  $fcc1, ALPHAI, a1
    bceqz $fcc0, .L13
    b .L14
    .align 3

.L13:
    bceqz $fcc1, .L114 //alpha_r != 0.0 && alpha_i != 0.0
    b .L113 //alpha_r != 0.0 && alpha_i == 0.0

.L14:
    bceqz $fcc1, .L114  //alpha_r == 0.0 && alpha_i != 0.0
    b .L111 //alpha_r == 0.0 && alpha_i == 0.0
    .align 3

.L111:  //alpha_r == 0.0 && alpha_i == 0.0
    vst VXZ, X, 0 * SIZE
#ifdef DOUBLE
    vst VXZ, X, 2 * SIZE
    vst VXZ, X, 4 * SIZE
    vst VXZ, X, 6 * SIZE
#else
    vst VXZ, X, 4 * SIZE
#endif
    addi.d X, X, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L111
    b .L997
    .align 3

.L113: //alpha_r != 0.0 && alpha_i == 0.0
    vld VX0, X, 0 * SIZE
#ifdef DOUBLE
    vld VX1, X, 2 * SIZE
    vpickev.d x1, VX1, VX0
    vpickod.d x2, VX1, VX0
    vfmul.d x3, VXAR, x1
    vfmul.d x4, VXAR, x2
    vilvl.d VX2, x4 ,x3
    vilvh.d VX3, x4, x3
    vst VX2, X, 0 * SIZE
    vst VX3, X, 2 * SIZE
    vld VX0, X, 4 * SIZE
    vld VX1, X, 6 * SIZE
    vpickev.d x1, VX1, VX0
    vpickod.d x2, VX1, VX0
    vfmul.d x3, VXAR, x1
    vfmul.d x4, VXAR, x2
    vilvl.d VX2, x4 ,x3
    vilvh.d VX3, x4, x3
    vst VX2, X, 4 * SIZE
    vst VX3, X, 6 * SIZE
#else
    vld VX1, X, 4 * SIZE
    vpickev.w x1, VX1, VX0
    vpickod.w x2, VX1, VX0
    vfmul.s x3, VXAR, x1
    vfmul.s x4, VXAR, x2
    vilvl.w VX2, x4 ,x3
    vilvh.w VX3, x4, x3
    vst VX2, X, 0 * SIZE
    vst VX3, X, 4 * SIZE
#endif
    addi.d X, X, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L113
    b .L997
    .align 3

.L114:  //alpha_r != 0.0 && alpha_i != 0.0
    vld VX0, X, 0 * SIZE
#ifdef DOUBLE
    vld VX1, X, 2 * SIZE
    vpickev.d x1, VX1, VX0
    vpickod.d x2, VX1, VX0
    vfmul.d VX0, VXAI, x2
    vfmsub.d x3, VXAR, x1, VX0
    vfmul.d VX1, VXAI, x1
    vfmadd.d x4, VXAR, x2, VX1
    vilvl.d VX2, x4 ,x3
    vilvh.d VX3, x4, x3
    vst VX2, X, 0 * SIZE
    vst VX3, X, 2 * SIZE
    vld VX0, X, 4 * SIZE
    vld VX1, X, 6 * SIZE
    vpickev.d x1, VX1, VX0
    vpickod.d x2, VX1, VX0
    vfmul.d VX0, VXAI, x2
    vfmsub.d x3, VXAR, x1, VX0
    vfmul.d VX1, VXAI, x1
    vfmadd.d x4, VXAR, x2, VX1
    vilvl.d VX2, x4 ,x3
    vilvh.d VX3, x4, x3
    vst VX2, X, 4 * SIZE
    vst VX3, X, 6 * SIZE
#else
    vld VX1, X, 4 * SIZE
    vpickev.w x1, VX1, VX0
    vpickod.w x2, VX1, VX0
    vfmul.s VX0, VXAI, x2
    vfmsub.s x3, VXAR, x1, VX0
    vfmul.s VX1, VXAI, x1
    vfmadd.s x4, VXAR, x2, VX1
    vilvl.w VX2, x4 ,x3
    vilvh.w VX3, x4, x3
    vst VX2, X, 0 * SIZE
    vst VX3, X, 4 * SIZE
#endif
    addi.d X, X, 8 * SIZE
    addi.d  I, I, -1
    blt $r0, I, .L114
    b .L997
    .align 3

.L22:
    bge $r0, I, .L997
    move XX, X
    CMPEQ  $fcc0, ALPHAR, a1
    CMPEQ  $fcc1, ALPHAI, a1
    bceqz $fcc0, .L23
    b .L24
    .align 3

.L23:
    bceqz $fcc1, .L224 //alpha_r != 0.0 && alpha_i != 0.0
    b .L223 //alpha_r != 0.0 && alpha_i == 0.0

.L24:
    bceqz $fcc1, .L224  //alpha_r == 0.0 && alpha_i != 0.0
    b .L221 //alpha_r == 0.0 && alpha_i == 0.0
    .align 3

.L221:  //alpha_r == 0.0 && alpha_i == 0.0
#ifdef DOUBLE
    vstelm.d VXZ, X, 0, 0
    vstelm.d VXZ, X, 1 * SIZE, 0
    add.d X, X, INCX
    vstelm.d VXZ, X, 0, 0
    vstelm.d VXZ, X, 1 * SIZE, 0
    add.d X, X, INCX
    vstelm.d VXZ, X, 0, 0
    vstelm.d VXZ, X, 1 * SIZE, 0
    add.d X, X, INCX
    vstelm.d VXZ, X, 0, 0
    vstelm.d VXZ, X, 1 * SIZE, 0
#else
    vstelm.w VXZ, X, 0, 0
    vstelm.w VXZ, X, 1 * SIZE, 0
    add.d X, X, INCX
    vstelm.w VXZ, X, 0, 0
    vstelm.w VXZ, X, 1 * SIZE, 0
    add.d X, X, INCX
    vstelm.w VXZ, X, 0, 0
    vstelm.w VXZ, X, 1 * SIZE, 0
    add.d X, X, INCX
    vstelm.w VXZ, X, 0, 0
    vstelm.w VXZ, X, 1 * SIZE, 0
#endif
    add.d X, X, INCX
    addi.d  I, I, -1
    blt $r0, I, .L221
    b .L997
    .align 3

.L223: //alpha_r != 0.0 && alpha_i == 0.0
#ifdef DOUBLE
    ld.d t1, X, 0 * SIZE
    ld.d t2, X, 1 * SIZE
    add.d X, X, INCX
    ld.d t3, X, 0 * SIZE
    ld.d t4, X, 1 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d x1, t1, 0
    vinsgr2vr.d x2, t2, 0
    vinsgr2vr.d x1, t3, 1
    vinsgr2vr.d x2, t4, 1
    vfmul.d x3, VXAR, x1
    vfmul.d x4, VXAR, x2
    vstelm.d x3, XX, 0 * SIZE, 0
    vstelm.d x4, XX, 1 * SIZE, 0
    add.d XX, XX, INCX
    vstelm.d x3, XX, 0 * SIZE, 1
    vstelm.d x4, XX, 1 * SIZE, 1
    add.d XX, XX, INCX

    ld.d t1, X, 0 * SIZE
    ld.d t2, X, 1 * SIZE
    add.d X, X, INCX
    ld.d t3, X, 0 * SIZE
    ld.d t4, X, 1 * SIZE
    vinsgr2vr.d x1, t1, 0
    vinsgr2vr.d x2, t2, 0
    vinsgr2vr.d x1, t3, 1
    vinsgr2vr.d x2, t4, 1
    add.d X, X, INCX
    vfmul.d x3, VXAR, x1
    vfmul.d x4, VXAR, x2
    addi.d  I, I, -1
    vstelm.d x3, XX, 0 * SIZE, 0
    vstelm.d x4, XX, 1 * SIZE, 0
    add.d XX, XX, INCX
    vstelm.d x3, XX, 0 * SIZE, 1
    vstelm.d x4, XX, 1 * SIZE, 1
#else
    ld.w t1, X, 0 * SIZE
    ld.w t2, X, 1 * SIZE
    add.d X, X, INCX
    ld.w t3, X, 0 * SIZE
    ld.w t4, X, 1 * SIZE
    add.d X, X, INCX
    vinsgr2vr.w x1, t1, 0
    vinsgr2vr.w x2, t2, 0
    vinsgr2vr.w x1, t3, 1
    vinsgr2vr.w x2, t4, 1
    ld.w t1, X, 0 * SIZE
    ld.w t2, X, 1 * SIZE
    add.d X, X, INCX
    ld.w t3, X, 0 * SIZE
    ld.w t4, X, 1 * SIZE
    vinsgr2vr.w x1, t1, 2
    vinsgr2vr.w x2, t2, 2
    vinsgr2vr.w x1, t3, 3
    vinsgr2vr.w x2, t4, 3
    add.d X, X, INCX

    vfmul.s x3, VXAR, x1
    vfmul.s x4, VXAR, x2
    addi.d  I, I, -1
    vstelm.w x3, XX, 0 * SIZE, 0
    vstelm.w x4, XX, 1 * SIZE, 0
    add.d XX, XX, INCX
    vstelm.w x3, XX, 0 * SIZE, 1
    vstelm.w x4, XX, 1 * SIZE, 1
    add.d XX, XX, INCX
    vstelm.w x3, XX, 0 * SIZE, 2
    vstelm.w x4, XX, 1 * SIZE, 2
    add.d XX, XX, INCX
    vstelm.w x3, XX, 0 * SIZE, 3
    vstelm.w x4, XX, 1 * SIZE, 3
#endif
    add.d XX, XX, INCX
    blt $r0, I, .L223
    b .L997
    .align 3

.L224:  //alpha_r != 0.0 && alpha_i != 0.0
#ifdef DOUBLE
    ld.d t1, X, 0 * SIZE
    ld.d t2, X, 1 * SIZE
    add.d X, X, INCX
    ld.d t3, X, 0 * SIZE
    ld.d t4, X, 1 * SIZE
    add.d X, X, INCX
    vinsgr2vr.d x1, t1, 0
    vinsgr2vr.d x2, t2, 0
    vinsgr2vr.d x1, t3, 1
    vinsgr2vr.d x2, t4, 1
    vfmul.d VX0, VXAI, x2
    vfmsub.d x3, VXAR, x1, VX0
    vfmul.d VX1, VXAI, x1
    vfmadd.d x4, VXAR, x2, VX1
    vstelm.d x3, XX, 0 * SIZE, 0
    vstelm.d x4, XX, 1 * SIZE, 0
    add.d XX, XX, INCX
    vstelm.d x3, XX, 0 * SIZE, 1
    vstelm.d x4, XX, 1 * SIZE, 1
    add.d XX, XX, INCX

    ld.d t1, X, 0 * SIZE
    ld.d t2, X, 1 * SIZE
    add.d X, X, INCX
    ld.d t3, X, 0 * SIZE
    ld.d t4, X, 1 * SIZE
    vinsgr2vr.d x1, t1, 0
    vinsgr2vr.d x2, t2, 0
    vinsgr2vr.d x1, t3, 1
    vinsgr2vr.d x2, t4, 1
    add.d X, X, INCX
    vfmul.d VX0, VXAI, x2
    vfmsub.d x3, VXAR, x1, VX0
    vfmul.d VX1, VXAI, x1
    vfmadd.d x4, VXAR, x2, VX1
    addi.d  I, I, -1
    vstelm.d x3, XX, 0 * SIZE, 0
    vstelm.d x4, XX, 1 * SIZE, 0
    add.d XX, XX, INCX
    vstelm.d x3, XX, 0 * SIZE, 1
    vstelm.d x4, XX, 1 * SIZE, 1
#else
    ld.w t1, X, 0 * SIZE
    ld.w t2, X, 1 * SIZE
    add.d X, X, INCX
    ld.w t3, X, 0 * SIZE
    ld.w t4, X, 1 * SIZE
    add.d X, X, INCX
    vinsgr2vr.w x1, t1, 0
    vinsgr2vr.w x2, t2, 0
    vinsgr2vr.w x1, t3, 1
    vinsgr2vr.w x2, t4, 1
    ld.w t1, X, 0 * SIZE
    ld.w t2, X, 1 * SIZE
    add.d X, X, INCX
    ld.w t3, X, 0 * SIZE
    ld.w t4, X, 1 * SIZE
    vinsgr2vr.w x1, t1, 2
    vinsgr2vr.w x2, t2, 2
    vinsgr2vr.w x1, t3, 3
    vinsgr2vr.w x2, t4, 3
    add.d X, X, INCX

    vfmul.s VX0, VXAI, x2
    vfmsub.s x3, VXAR, x1, VX0
    vfmul.s VX1, VXAI, x1
    vfmadd.s x4, VXAR, x2, VX1
    addi.d  I, I, -1
    vstelm.w x3, XX, 0 * SIZE, 0
    vstelm.w x4, XX, 1 * SIZE, 0
    add.d XX, XX, INCX
    vstelm.w x3, XX, 0 * SIZE, 1
    vstelm.w x4, XX, 1 * SIZE, 1
    add.d XX, XX, INCX
    vstelm.w x3, XX, 0 * SIZE, 2
    vstelm.w x4, XX, 1 * SIZE, 2
    add.d XX, XX, INCX
    vstelm.w x3, XX, 0 * SIZE, 3
    vstelm.w x4, XX, 1 * SIZE, 3
#endif
    add.d XX, XX, INCX
    blt $r0, I, .L224
    b .L997
    .align 3

.L997:
    andi I, N, 3
    bge $r0, I, .L999
    .align 3

.L998:
    LD   a1, X, 0 * SIZE
    LD   a2, X, 1 * SIZE
    addi.d I, I, -1
    MUL  s1, ALPHAI, a2
    MUL  s2, ALPHAI, a1
    MSUB s1, ALPHAR, a1, s1
    MADD s2, ALPHAR, a2, s2
    ST   s1, X, 0 * SIZE
    ST   s2, X, 1 * SIZE
    add.d X, X, INCX
    blt $r0, I, .L998
    .align 3

.L999:
    move $r4, $r12
    jirl $r0, $r1, 0x0
    .align 3

    EPILOGUE
