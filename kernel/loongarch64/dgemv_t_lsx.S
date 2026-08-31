/*******************************************************************************
Copyright (c) 2024, The OpenBLAS Project
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
*******************************************************************************/

#define ASSEMBLER
#include "common.h"

#define M      $r4
#define N      $r5
#define A      $r7
#define LDA    $r8
#define X      $r9
#define INCX   $r10
#define Y      $r11
#define INCY   $r6
#define BUFFER $r16
#define ALPHA  $f0

#define I      $r14
#define J      $r15
#define AO     $r23
#define XORG   $r24
#define IX     $r25
#define IY     $r26
#define II     $r27
#define T0     $r19
#define T1     $r20
#define T2     $r28
#define T3     $r29
#define T4     $r30

#define ACC    $vr11
#define VX     $vr12
#define VA     $vr13
#define VT     $vr14
#define FAC    $f17
#define FAX    $f18
#define FAY    $f19

    PROLOGUE

    LDARG  INCY,    $sp,  0
    LDARG  BUFFER,  $sp,  8

    addi.d  $sp, $sp, -80
    SDARG  $r23,  $sp,    0
    SDARG  $r24,  $sp,    8
    SDARG  $r25,  $sp,    16
    SDARG  $r26,  $sp,    32
    SDARG  $r27,  $sp,    40
    SDARG  $r28,  $sp,    48
    SDARG  $r29,  $sp,    56
    SDARG  $r30,  $sp,    64
    ST     ALPHA, $sp,    72

    slli.d     LDA,     LDA,  BASE_SHIFT
    slli.d     INCX,    INCX, BASE_SHIFT
    slli.d     INCY,    INCY, BASE_SHIFT

    bge    $r0,    M,    .L999
    bge    $r0,    N,    .L999

    move   J,      $r0
    move   IY,     $r0
    move   XORG,   X
    move   AO,     A
    li.d   T1,     8

.L01:  /* for each column J */
    vxor.v    ACC,    ACC,   ACC
    MTC       FAC,    $r0
    move      X,      XORG
    move      II,     $r0
    move      I,      $r0
    srli.d    T2,     M,     1
    slli.d    T2,     T2,    1       /* T2 = M & ~1 */

    bne       INCX,   T1,    .L_SCALAR

.L02:
    beq       I,      T2,    .L03
    vld       VX,     X,     0
    vld       VA,     AO,    0
    vfmadd.d  ACC,    VA,    VX,   ACC
    addi.d    II,     II,    16
    addi.d    I,      I,     2
    addi.d    X,      X,     16
    addi.d    AO,     AO,    16
    b         .L02

.L03:
    andi      T3,     M,     1
    beq       $r0,    T3,    .L05
    fld.d     FAY,    AO,    0
    fld.d     FAX,    X,     0
    fmadd.d   FAC,    FAY,   FAX,  FAC
.L05:
    /* dot = (lane0+lane1) + tail */
    vextrins.d VT,    ACC,   0x01
    vfadd.d   VT,     ACC,    VT
    vpickve2gr.d T0,  VT,    0
    movgr2fr.d  FAX,  T0
    fadd.d      FAC,  FAX,   FAC
    /* y += alpha*dot */
    fldx.d    FAY,    Y,     IY
    fmadd.d   FAY,    ALPHA, FAC,  FAY
    fstx.d    FAY,    Y,     IY
    /* next column: AO = A + (J+1)*LDA */
    addi.d    T3,     J,     1
    mul.d     T3,     T3,    LDA
    add.d     AO,     A,     T3
    add.d     IY,     IY,    INCY
    addi.d    J,      J,     1
    blt       J,      N,     .L01
    b         .L999

.L_SCALAR:
    move      IX,     $r0
    move      I,      $r0
    move      II,     $r0
.L06:
    beq       I,      M,     .L07
    fldx.d    FAX,    X,     IX
    fldx.d    FAY,    AO,    II
    fmadd.d   FAC,    FAY,   FAX,  FAC
    add.d     IX,     IX,    INCX
    addi.d    II,     II,    8
    addi.d    I,      I,     1
    b         .L06
.L07:
    fldx.d    FAY,    Y,     IY
    fmadd.d   FAY,    ALPHA, FAC,  FAY
    fstx.d    FAY,    Y,     IY
    addi.d    T3,     J,     1
    mul.d     T3,     T3,    LDA
    add.d     AO,     A,     T3
    add.d     IY,     IY,    INCY
    addi.d    J,      J,     1
    blt       J,      N,     .L01
    b         .L999

.L999:
    LDARG  $r23,  $sp,    0
    LDARG  $r24,  $sp,    8
    LDARG  $r25,  $sp,    16
    LDARG  $r26,  $sp,    32
    LDARG  $r27,  $sp,    40
    LDARG  $r28,  $sp,    48
    LDARG  $r29,  $sp,    56
    LDARG  $r30,  $sp,    64
    LD     ALPHA, $sp,    72
    addi.d  $sp, $sp, 80
    jirl    $r0, $r1, 0x0

    EPILOGUE
