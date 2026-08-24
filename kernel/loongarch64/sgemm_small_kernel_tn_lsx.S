/***************************************************************************
Copyright (c) 2026, The OpenBLAS Project
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
#include "loongarch64_asm.S"

#define M       $a0
#define N       $a1
#define K       $a2
#define A       $a3
#define LDA     $a4
#define ALPHA   $f0
#define B       $a5
#define LDB     $a6
#define C       $a7
#define LDC     $t0
#ifdef B0
#define BETA    $f1
#endif
#undef ZERO
#define ZERO    $r0

#define M16     $t1
#define M8      $t1
#define M4      $t1
#define M2      $t1
#define M1      $t1
#define N4      $t2
#define N2      $t2
#define N1      $t2
#define K4      $t3
#define A0      $t4
#define X0      $t5
#define B1      $t6
#define B2      $t7
#define B3      $t8
#define C0      $s0
#define C1      $s1
#define C2      $s2
#define C3      $s3
#define K1      $s4

#define A_cur   $s5
#define TEMP    $s6

#define VALPHA  $vr0
#ifndef B0
#define VBETA   $vr1
#endif
#define D0     $vr2
#define D1     $vr3
#define D2     $vr4
#define D3     $vr5
#define D4     $vr6
#define D5     $vr7
#define D6     $vr8
#define D7     $vr9
#define D8     $vr10
#define D9     $vr11
#define D10    $vr12
#define D11    $vr13
#define D12    $vr14
#define D13    $vr15
#define D14    $vr16
#define D15    $vr17
#define S0     $vr18
#define S1     $vr19
#define S2     $vr20
#define S3     $vr21
#define Z0     $vr22
#define Z1     $vr23
#define Z2     $vr24
#define Z3     $vr25

#define V0     $vr2
#define V1     $vr3
#define V2     $vr4
#define V3     $vr5

#define F0     $f2
#define F1     $f3
#define F2     $f4
#define F3     $f5
#define FS0    $f18

.macro SGEMM_SMALL_KERNEL_TN_TAIL M
    PTR_SRAI    N4,     N,      2
    move        A0,     A
    move        X0,     B
    PTR_ADD     B1,     X0,     LDB
    PTR_ADD     B2,     B1,     LDB
    PTR_ADD     B3,     B2,     LDB
    move        C0,     C
    PTR_ADD     C1,     C0,     LDC
    PTR_ADD     C2,     C1,     LDC
    PTR_ADD     C3,     C2,     LDC
    beqz        N4,     .L_TN_M\M\()_N3
.L_TN_M\M\()_N4:
    vxor.v      D0, D0, D0
    vxor.v      D1, D1, D1
    vxor.v      D2, D2, D2
    vxor.v      D3, D3, D3
    move        K1,     K
    PTR_ADDI    N4,     N4,     -1
    bge         ZERO,   K,      .L_TN_M\M\()_N4_END
.L_TN_M\M\()_N4_K1:
    PTR_ADDI    K1,     K1,    -1
    move        A_cur,  A0
    vxor.v      S0, S0, S0
.if \M >= 1
    ld.w        TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 0
.endif
.if \M >= 2
    PTR_ADD     A_cur, A_cur, LDA; ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 1
.endif
.if \M >= 4
    PTR_ADD     A_cur, A_cur, LDA; ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 2
    PTR_ADD     A_cur, A_cur, LDA; ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 3
.endif
    vldrepl.w   Z0, X0, 0x00
    vldrepl.w   Z1, B1, 0x00
    vldrepl.w   Z2, B2, 0x00
    vldrepl.w   Z3, B3, 0x00
    vfmadd.s    D0, S0, Z0, D0
    vfmadd.s    D1, S0, Z1, D1
    vfmadd.s    D2, S0, Z2, D2
    vfmadd.s    D3, S0, Z3, D3
    PTR_ADDI    X0,     X0,     0x04
    PTR_ADDI    B1,     B1,     0x04
    PTR_ADDI    B2,     B2,     0x04
    PTR_ADDI    B3,     B3,     0x04
    PTR_ADDI    A0,     A0,     0x04
    bnez        K1,     .L_TN_M\M\()_N4_K1
.L_TN_M\M\()_N4_END:
    vfmul.s     D0, D0, VALPHA
    vfmul.s     D1, D1, VALPHA
    vfmul.s     D2, D2, VALPHA
    vfmul.s     D3, D3, VALPHA
#ifndef B0
.if \M == 4
    vld         S0, C0, 0x00; vfmadd.s D0, S0, VBETA, D0
    vld         S0, C1, 0x00; vfmadd.s D1, S0, VBETA, D1
    vld         S0, C2, 0x00; vfmadd.s D2, S0, VBETA, D2
    vld         S0, C3, 0x00; vfmadd.s D3, S0, VBETA, D3
.elseif \M == 2
    fld.d       FS0, C0, 0x00; vfmadd.s D0, S0, VBETA, D0
    fld.d       FS0, C1, 0x00; vfmadd.s D1, S0, VBETA, D1
    fld.d       FS0, C2, 0x00; vfmadd.s D2, S0, VBETA, D2
    fld.d       FS0, C3, 0x00; vfmadd.s D3, S0, VBETA, D3
.elseif \M == 1
    fld.s       FS0, C0, 0x00; vfmadd.s D0, S0, VBETA, D0
    fld.s       FS0, C1, 0x00; vfmadd.s D1, S0, VBETA, D1
    fld.s       FS0, C2, 0x00; vfmadd.s D2, S0, VBETA, D2
    fld.s       FS0, C3, 0x00; vfmadd.s D3, S0, VBETA, D3
.endif
#endif
.if \M == 4
    vst         V0, C0, 0x00
    vst         V1, C1, 0x00
    vst         V2, C2, 0x00
    vst         V3, C3, 0x00
.elseif \M == 2
    fst.d       F0, C0, 0x00
    fst.d       F1, C1, 0x00
    fst.d       F2, C2, 0x00
    fst.d       F3, C3, 0x00
.elseif \M == 1
    fst.s       F0, C0, 0x00
    fst.s       F1, C1, 0x00
    fst.s       F2, C2, 0x00
    fst.s       F3, C3, 0x00
.endif
    PTR_ALSL    C0,     LDC,    C0,     2
    PTR_ALSL    C1,     LDC,    C1,     2
    PTR_ALSL    C2,     LDC,    C2,     2
    PTR_ALSL    C3,     LDC,    C3,     2
    PTR_SUB     X0,     X0,    K4
    PTR_SUB     B1,     B1,    K4
    PTR_SUB     B2,     B2,    K4
    PTR_SUB     B3,     B3,    K4
    PTR_ALSL    X0,     LDB,    X0,     2
    PTR_ALSL    B1,     LDB,    B1,     2
    PTR_ALSL    B2,     LDB,    B2,     2
    PTR_ALSL    B3,     LDB,    B3,     2
    move        A0,     A
    bnez        N4,     .L_TN_M\M\()_N4
.L_TN_M\M\()_N3:
    andi        N2,     N,      0x02
    beqz        N2,     .L_TN_M\M\()_N1
.L_TN_M\M\()_N2:
    vxor.v      D0, D0, D0
    vxor.v      D1, D1, D1
    move        K1,     K
    bge         ZERO,   K,      .L_TN_M\M\()_N2_END
.L_TN_M\M\()_N2_K1:
    PTR_ADDI    K1,     K1,    -1
    move        A_cur,  A0
    vxor.v      S0, S0, S0
.if \M >= 1
    ld.w        TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 0
.endif
.if \M >= 2
    PTR_ADD     A_cur, A_cur, LDA; ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 1
.endif
.if \M >= 4
    PTR_ADD     A_cur, A_cur, LDA; ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 2
    PTR_ADD     A_cur, A_cur, LDA; ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 3
.endif
    vldrepl.w   Z0, X0, 0x00
    vldrepl.w   Z1, B1, 0x00
    vfmadd.s    D0, S0, Z0, D0
    vfmadd.s    D1, S0, Z1, D1
    PTR_ADDI    X0,     X0,     0x04
    PTR_ADDI    B1,     B1,     0x04
    PTR_ADDI    A0,     A0,     0x04
    bnez        K1,     .L_TN_M\M\()_N2_K1
.L_TN_M\M\()_N2_END:
    vfmul.s     D0,  D0, VALPHA
    vfmul.s     D1,  D1, VALPHA
#ifndef B0
.if \M == 4
    vld         S0, C0, 0x00; vfmadd.s D0, S0, VBETA, D0
    vld         S0, C1, 0x00; vfmadd.s D1, S0, VBETA, D1
.elseif \M == 2
    fld.d       FS0, C0, 0x00; vfmadd.s D0, S0, VBETA, D0
    fld.d       FS0, C1, 0x00; vfmadd.s D1, S0, VBETA, D1
.elseif \M == 1
    fld.s       FS0, C0, 0x00; vfmadd.s D0, S0, VBETA, D0
    fld.s       FS0, C1, 0x00; vfmadd.s D1, S0, VBETA, D1
.endif
#endif
.if \M == 4
    vst         V0, C0, 0x00
    vst         V1, C1, 0x00
.elseif \M == 2
    fst.d       F0, C0, 0x00
    fst.d       F1, C1, 0x00
.elseif \M == 1
    fst.s       F0, C0, 0x00
    fst.s       F1, C1, 0x00
.endif
    PTR_ALSL    C0,     LDC,    C0,     1
    PTR_ALSL    C1,     LDC,    C1,     1
    PTR_SUB     X0,     X0,    K4
    PTR_SUB     B1,     B1,    K4
    PTR_ALSL    X0,     LDB,   X0,      1
    PTR_ALSL    B1,     LDB,   B1,      1
    move        A0,     A
.L_TN_M\M\()_N1:
    andi        N1,     N,      0x01
    beqz        N1,     .L_TN_M\M\()_END
    vxor.v      D0, D0, D0
    move        K1,     K
    bge         ZERO,   K,      .L_TN_M\M\()_N1_END
.L_TN_M\M\()_N1_K1:
    PTR_ADDI    K1,     K1,    -1
    move        A_cur,  A0
    vxor.v      S0, S0, S0
.if \M >= 1
    ld.w        TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 0
.endif
.if \M >= 2
    PTR_ADD     A_cur, A_cur, LDA; ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 1
.endif
.if \M >= 4
    PTR_ADD     A_cur, A_cur, LDA; ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 2
    PTR_ADD     A_cur, A_cur, LDA; ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 3
.endif
    vldrepl.w   Z0, X0, 0x00
    vfmadd.s    D0, S0, Z0, D0
    PTR_ADDI    X0,     X0,     0x04
    PTR_ADDI    A0,     A0,     0x04
    bnez        K1,     .L_TN_M\M\()_N1_K1
.L_TN_M\M\()_N1_END:
    vfmul.s     D0, D0, VALPHA
#ifndef B0
.if \M == 4
    vld         S0, C0, 0x00; vfmadd.s D0, S0, VBETA, D0
.elseif \M == 2
    fld.d       FS0, C0, 0x00; vfmadd.s D0, S0, VBETA, D0
.elseif \M == 1
    fld.s       FS0, C0, 0x00; vfmadd.s D0, S0, VBETA, D0
.endif
#endif
.if \M == 4
    vst         V0, C0, 0x00
.elseif \M == 2
    fst.d       F0, C0, 0x00
.elseif \M == 1
    fst.s       F0, C0, 0x00
.endif
.L_TN_M\M\()_END:
.if \M == 4
    PTR_SLLI    TEMP,   LDA,    2
    PTR_ADD     A,      A,      TEMP
    PTR_ADDI    C,      C,     0x10
.elseif \M == 2
    PTR_SLLI    TEMP,   LDA,    1
    PTR_ADD     A,      A,      TEMP
    PTR_ADDI    C,      C,     0x08
.elseif \M == 1
    PTR_ADD     A,      A,      LDA
    PTR_ADDI    C,      C,     0x04
.endif
.endm

    PROLOGUE
    PTR_LD          LDC,    $sp,    0
    push_if_used    7, 2
    vreplvei.w  VALPHA, VALPHA, 0
#ifndef B0
    vreplvei.w  VBETA,  VBETA,  0
#endif
    PTR_SLLI    LDA,    LDA,    2
    PTR_SLLI    LDB,    LDB,    2
    PTR_SLLI    LDC,    LDC,    2
    PTR_SLLI    K4,     K,      2
    PTR_SRAI    M16,    M,      4
    beqz        M16,    .L_M15
.L_M16:
    PTR_SRAI    N4,     N,      2
    move        A0,     A
    move        X0,     B
    PTR_ADD     B1,     X0,     LDB
    PTR_ADD     B2,     B1,     LDB
    PTR_ADD     B3,     B2,     LDB
    move        C0,     C
    PTR_ADD     C1,     C0,     LDC
    PTR_ADD     C2,     C1,     LDC
    PTR_ADD     C3,     C2,     LDC
    beqz        N4,     .L_M16_N3
.L_M16_N4:
    vxor.v      D0, D0, D0; vxor.v D1, D1, D1; vxor.v D2, D2, D2; vxor.v D3, D3, D3
    vxor.v      D4, D4, D4; vxor.v D5, D5, D5; vxor.v D6, D6, D6; vxor.v D7, D7, D7
    vxor.v      D8, D8, D8; vxor.v D9, D9, D9; vxor.v D10, D10, D10; vxor.v D11, D11, D11
    vxor.v      D12, D12, D12; vxor.v D13, D13, D13; vxor.v D14, D14, D14; vxor.v D15, D15, D15
    move        K1,     K
    PTR_ADDI    N4,     N4,     -1
    bge         ZERO,   K,      .L_M16_N4_END
.L_M16_N4_K1:
    PTR_ADDI    K1,     K1,    -1
    move        A_cur,  A0
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 0; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 1; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 2; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 3; PTR_ADD A_cur, A_cur, LDA

    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 0; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 1; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 2; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 3; PTR_ADD A_cur, A_cur, LDA

    ld.w TEMP, A_cur, 0; vinsgr2vr.w S2, TEMP, 0; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S2, TEMP, 1; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S2, TEMP, 2; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S2, TEMP, 3; PTR_ADD A_cur, A_cur, LDA

    ld.w TEMP, A_cur, 0; vinsgr2vr.w S3, TEMP, 0; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S3, TEMP, 1; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S3, TEMP, 2; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S3, TEMP, 3

    vldrepl.w   Z0, X0, 0x00
    vldrepl.w   Z1, B1, 0x00
    vldrepl.w   Z2, B2, 0x00
    vldrepl.w   Z3, B3, 0x00

    vfmadd.s    D0,  S0, Z0, D0; vfmadd.s D1,  S1, Z0, D1; vfmadd.s D2,  S2, Z0, D2; vfmadd.s D3,  S3, Z0, D3
    vfmadd.s    D4,  S0, Z1, D4; vfmadd.s D5,  S1, Z1, D5; vfmadd.s D6,  S2, Z1, D6; vfmadd.s D7,  S3, Z1, D7
    vfmadd.s    D8,  S0, Z2, D8; vfmadd.s D9,  S1, Z2, D9; vfmadd.s D10, S2, Z2, D10; vfmadd.s D11, S3, Z2, D11
    vfmadd.s    D12, S0, Z3, D12; vfmadd.s D13, S1, Z3, D13; vfmadd.s D14, S2, Z3, D14; vfmadd.s D15, S3, Z3, D15

    PTR_ADDI    X0,     X0,     0x04
    PTR_ADDI    B1,     B1,     0x04
    PTR_ADDI    B2,     B2,     0x04
    PTR_ADDI    B3,     B3,     0x04
    PTR_ADDI    A0,     A0,     0x04
    bnez        K1,     .L_M16_N4_K1
.L_M16_N4_END:
    vfmul.s     D0,  D0,  VALPHA; vfmul.s D1,  D1,  VALPHA; vfmul.s D2,  D2,  VALPHA; vfmul.s D3,  D3,  VALPHA
    vfmul.s     D4,  D4,  VALPHA; vfmul.s D5,  D5,  VALPHA; vfmul.s D6,  D6,  VALPHA; vfmul.s D7,  D7,  VALPHA
    vfmul.s     D8,  D8,  VALPHA; vfmul.s D9,  D9,  VALPHA; vfmul.s D10, D10, VALPHA; vfmul.s D11, D11, VALPHA
    vfmul.s     D12, D12, VALPHA; vfmul.s D13, D13, VALPHA; vfmul.s D14, D14, VALPHA; vfmul.s D15, D15, VALPHA
#ifndef B0
    vld         S0, C0, 0x00; vld S1, C0, 0x10; vld S2, C0, 0x20; vld S3, C0, 0x30
    vfmadd.s    D0, S0, VBETA, D0; vfmadd.s D1, S1, VBETA, D1; vfmadd.s D2, S2, VBETA, D2; vfmadd.s D3, S3, VBETA, D3

    vld         S0, C1, 0x00; vld S1, C1, 0x10; vld S2, C1, 0x20; vld S3, C1, 0x30
    vfmadd.s    D4, S0, VBETA, D4; vfmadd.s D5, S1, VBETA, D5; vfmadd.s D6, S2, VBETA, D6; vfmadd.s D7, S3, VBETA, D7

    vld         S0, C2, 0x00; vld S1, C2, 0x10; vld S2, C2, 0x20; vld S3, C2, 0x30
    vfmadd.s    D8, S0, VBETA, D8; vfmadd.s D9, S1, VBETA, D9; vfmadd.s D10, S2, VBETA, D10; vfmadd.s D11, S3, VBETA, D11

    vld         S0, C3, 0x00; vld S1, C3, 0x10; vld S2, C3, 0x20; vld S3, C3, 0x30
    vfmadd.s    D12, S0, VBETA, D12; vfmadd.s D13, S1, VBETA, D13; vfmadd.s D14, S2, VBETA, D14; vfmadd.s D15, S3, VBETA, D15
#endif
    vst         D0,  C0, 0x00; vst D1,  C0, 0x10; vst D2,  C0, 0x20; vst D3,  C0, 0x30
    vst         D4,  C1, 0x00; vst D5,  C1, 0x10; vst D6,  C1, 0x20; vst D7,  C1, 0x30
    vst         D8,  C2, 0x00; vst D9,  C2, 0x10; vst D10, C2, 0x20; vst D11, C2, 0x30
    vst         D12, C3, 0x00; vst D13, C3, 0x10; vst D14, C3, 0x20; vst D15, C3, 0x30
    PTR_ALSL    C0,     LDC,    C0,     2
    PTR_ALSL    C1,     LDC,    C1,     2
    PTR_ALSL    C2,     LDC,    C2,     2
    PTR_ALSL    C3,     LDC,    C3,     2
    PTR_SUB     X0,     X0,    K4
    PTR_SUB     B1,     B1,    K4
    PTR_SUB     B2,     B2,    K4
    PTR_SUB     B3,     B3,    K4
    PTR_ALSL    X0,     LDB,    X0,     2
    PTR_ALSL    B1,     LDB,    B1,     2
    PTR_ALSL    B2,     LDB,    B2,     2
    PTR_ALSL    B3,     LDB,    B3,     2
    move        A0,     A
    bnez        N4,     .L_M16_N4
.L_M16_N3:
    andi        N2,     N,      0x02
    beqz        N2,     .L_M16_N1
.L_M16_N2:
    vxor.v      D0, D0, D0; vxor.v D1, D1, D1; vxor.v D2, D2, D2; vxor.v D3, D3, D3
    vxor.v      D4, D4, D4; vxor.v D5, D5, D5; vxor.v D6, D6, D6; vxor.v D7, D7, D7
    move        K1,     K
    bge         ZERO,   K,      .L_M16_N2_END
.L_M16_N2_K1:
    PTR_ADDI    K1,     K1,    -1
    move        A_cur,  A0
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 0; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 1; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 2; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 3; PTR_ADD A_cur, A_cur, LDA

    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 0; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 1; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 2; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 3; PTR_ADD A_cur, A_cur, LDA

    ld.w TEMP, A_cur, 0; vinsgr2vr.w S2, TEMP, 0; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S2, TEMP, 1; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S2, TEMP, 2; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S2, TEMP, 3; PTR_ADD A_cur, A_cur, LDA

    ld.w TEMP, A_cur, 0; vinsgr2vr.w S3, TEMP, 0; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S3, TEMP, 1; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S3, TEMP, 2; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S3, TEMP, 3

    vldrepl.w   Z0, X0, 0x00
    vldrepl.w   Z1, B1, 0x00

    vfmadd.s    D0,  S0, Z0, D0; vfmadd.s D1,  S1, Z0, D1; vfmadd.s D2,  S2, Z0, D2; vfmadd.s D3,  S3, Z0, D3
    vfmadd.s    D4,  S0, Z1, D4; vfmadd.s D5,  S1, Z1, D5; vfmadd.s D6,  S2, Z1, D6; vfmadd.s D7,  S3, Z1, D7

    PTR_ADDI    X0,     X0,     0x04
    PTR_ADDI    B1,     B1,     0x04
    PTR_ADDI    A0,     A0,     0x04
    bnez        K1,     .L_M16_N2_K1
.L_M16_N2_END:
    vfmul.s     D0,  D0,  VALPHA; vfmul.s D1,  D1,  VALPHA; vfmul.s D2,  D2,  VALPHA; vfmul.s D3,  D3,  VALPHA
    vfmul.s     D4,  D4,  VALPHA; vfmul.s D5,  D5,  VALPHA; vfmul.s D6,  D6,  VALPHA; vfmul.s D7,  D7,  VALPHA
#ifndef B0
    vld         S0, C0, 0x00; vld S1, C0, 0x10; vld S2, C0, 0x20; vld S3, C0, 0x30
    vfmadd.s    D0, S0, VBETA, D0; vfmadd.s D1, S1, VBETA, D1; vfmadd.s D2, S2, VBETA, D2; vfmadd.s D3, S3, VBETA, D3

    vld         S0, C1, 0x00; vld S1, C1, 0x10; vld S2, C1, 0x20; vld S3, C1, 0x30
    vfmadd.s    D4, S0, VBETA, D4; vfmadd.s D5, S1, VBETA, D5; vfmadd.s D6, S2, VBETA, D6; vfmadd.s D7, S3, VBETA, D7
#endif
    vst         D0,  C0, 0x00; vst D1,  C0, 0x10; vst D2,  C0, 0x20; vst D3,  C0, 0x30
    vst         D4,  C1, 0x00; vst D5,  C1, 0x10; vst D6,  C1, 0x20; vst D7,  C1, 0x30
    PTR_ALSL    C0,     LDC,    C0,     1
    PTR_ALSL    C1,     LDC,    C1,     1
    PTR_SUB     X0,     X0,    K4
    PTR_SUB     B1,     B1,    K4
    PTR_ALSL    X0,     LDB,    X0,     1
    PTR_ALSL    B1,     LDB,    B1,     1
    move        A0,     A
.L_M16_N1:
    andi        N1,     N,      0x01
    beqz        N1,     .L_M16_END
    vxor.v      D0, D0, D0; vxor.v D1, D1, D1; vxor.v D2, D2, D2; vxor.v D3, D3, D3
    move        K1,     K
    bge         ZERO,   K,      .L_M16_N1_END
.L_M16_N1_K1:
    PTR_ADDI    K1,     K1,    -1
    move        A_cur,  A0
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 0; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 1; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 2; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 3; PTR_ADD A_cur, A_cur, LDA

    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 0; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 1; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 2; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 3; PTR_ADD A_cur, A_cur, LDA

    ld.w TEMP, A_cur, 0; vinsgr2vr.w S2, TEMP, 0; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S2, TEMP, 1; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S2, TEMP, 2; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S2, TEMP, 3; PTR_ADD A_cur, A_cur, LDA

    ld.w TEMP, A_cur, 0; vinsgr2vr.w S3, TEMP, 0; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S3, TEMP, 1; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S3, TEMP, 2; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S3, TEMP, 3

    vldrepl.w   Z0, X0, 0x00
    vfmadd.s    D0,  S0, Z0, D0; vfmadd.s D1,  S1, Z0, D1; vfmadd.s D2,  S2, Z0, D2; vfmadd.s D3,  S3, Z0, D3
    PTR_ADDI    X0,     X0,     0x04
    PTR_ADDI    A0,     A0,     0x04
    bnez        K1,     .L_M16_N1_K1
.L_M16_N1_END:
    vfmul.s     D0,  D0,  VALPHA; vfmul.s D1,  D1,  VALPHA; vfmul.s D2,  D2,  VALPHA; vfmul.s D3,  D3,  VALPHA
#ifndef B0
    vld         S0, C0, 0x00; vld S1, C0, 0x10; vld S2, C0, 0x20; vld S3, C0, 0x30
    vfmadd.s    D0, S0, VBETA, D0; vfmadd.s D1, S1, VBETA, D1; vfmadd.s D2, S2, VBETA, D2; vfmadd.s D3, S3, VBETA, D3
#endif
    vst         D0,  C0, 0x00; vst D1,  C0, 0x10; vst D2,  C0, 0x20; vst D3,  C0, 0x30
    PTR_ALSL    C0,     LDC,    C0,     2
    PTR_SUB     X0,     X0,    K4
    PTR_ALSL    X0,     LDB,   X0,     2
    move        A0,     A
.L_M16_END:
    PTR_ADDI    M16,    M16,   -1
    PTR_SLLI    TEMP,   LDA,    4
    PTR_ADD     A,      A,      TEMP
    PTR_ADDI    C,      C,     0x40
    bnez        M16,    .L_M16
.L_M15:
    andi        M8,     M,      0x08
    beqz        M8,     .L_M7
.L_M8:
    PTR_SRAI    N4,     N,      2
    move        A0,     A
    move        X0,     B
    PTR_ADD     B1,     X0,     LDB
    PTR_ADD     B2,     B1,     LDB
    PTR_ADD     B3,     B2,     LDB
    move        C0,     C
    PTR_ADD     C1,     C0,     LDC
    PTR_ADD     C2,     C1,     LDC
    PTR_ADD     C3,     C2,     LDC
    beqz        N4,     .L_M8_N3
.L_M8_N4:
    vxor.v      D0, D0, D0; vxor.v D1, D1, D1; vxor.v D2, D2, D2; vxor.v D3, D3, D3
    vxor.v      D4, D4, D4; vxor.v D5, D5, D5; vxor.v D6, D6, D6; vxor.v D7, D7, D7
    move        K1,     K
    PTR_ADDI    N4,     N4,     -1
    bge         ZERO,   K,      .L_M8_N4_END
.L_M8_N4_K1:
    PTR_ADDI    K1,     K1,    -1
    move        A_cur,  A0
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 0; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 1; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 2; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 3; PTR_ADD A_cur, A_cur, LDA

    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 0; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 1; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 2; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 3

    vldrepl.w   Z0, X0, 0x00
    vldrepl.w   Z1, B1, 0x00
    vldrepl.w   Z2, B2, 0x00
    vldrepl.w   Z3, B3, 0x00

    vfmadd.s    D0, S0, Z0, D0; vfmadd.s D1, S1, Z0, D1
    vfmadd.s    D2, S0, Z1, D2; vfmadd.s D3, S1, Z1, D3
    vfmadd.s    D4, S0, Z2, D4; vfmadd.s D5, S1, Z2, D5
    vfmadd.s    D6, S0, Z3, D6; vfmadd.s D7, S1, Z3, D7

    PTR_ADDI    X0,     X0,     0x04
    PTR_ADDI    B1,     B1,     0x04
    PTR_ADDI    B2,     B2,     0x04
    PTR_ADDI    B3,     B3,     0x04
    PTR_ADDI    A0,     A0,     0x04
    bnez        K1,     .L_M8_N4_K1
.L_M8_N4_END:
    vfmul.s     D0,  D0,  VALPHA; vfmul.s D1,  D1,  VALPHA; vfmul.s D2,  D2,  VALPHA; vfmul.s D3,  D3,  VALPHA
    vfmul.s     D4,  D4,  VALPHA; vfmul.s D5,  D5,  VALPHA; vfmul.s D6,  D6,  VALPHA; vfmul.s D7,  D7,  VALPHA
#ifndef B0
    vld         S0, C0, 0x00; vld S1, C0, 0x10
    vfmadd.s    D0, S0, VBETA, D0; vfmadd.s D1, S1, VBETA, D1

    vld         S0, C1, 0x00; vld S1, C1, 0x10
    vfmadd.s    D2, S0, VBETA, D2; vfmadd.s D3, S1, VBETA, D3

    vld         S0, C2, 0x00; vld S1, C2, 0x10
    vfmadd.s    D4, S0, VBETA, D4; vfmadd.s D5, S1, VBETA, D5

    vld         S0, C3, 0x00; vld S1, C3, 0x10
    vfmadd.s    D6, S0, VBETA, D6; vfmadd.s D7, S1, VBETA, D7
#endif
    vst         D0,  C0, 0x00; vst D1,  C0, 0x10
    vst         D2,  C1, 0x00; vst D3,  C1, 0x10
    vst         D4,  C2, 0x00; vst D5,  C2, 0x10
    vst         D6,  C3, 0x00; vst D7,  C3, 0x10
    PTR_ALSL    C0,     LDC,    C0,     2
    PTR_ALSL    C1,     LDC,    C1,     2
    PTR_ALSL    C2,     LDC,    C2,     2
    PTR_ALSL    C3,     LDC,    C3,     2
    PTR_SUB     X0,     X0,    K4
    PTR_SUB     B1,     B1,    K4
    PTR_SUB     B2,     B2,    K4
    PTR_SUB     B3,     B3,    K4
    PTR_ALSL    X0,     LDB,    X0,     2
    PTR_ALSL    B1,     LDB,    B1,     2
    PTR_ALSL    B2,     LDB,    B2,     2
    PTR_ALSL    B3,     LDB,    B3,     2
    move        A0,     A
    bnez        N4,     .L_M8_N4
.L_M8_N3:
    andi        N2,     N,      0x02
    beqz        N2,     .L_M8_N1
.L_M8_N2:
    vxor.v      D0, D0, D0; vxor.v D1, D1, D1; vxor.v D2, D2, D2; vxor.v D3, D3, D3
    move        K1,     K
    bge         ZERO,   K,      .L_M8_N2_END
.L_M8_N2_K1:
    PTR_ADDI    K1,     K1,    -1
    move        A_cur,  A0
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 0; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 1; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 2; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 3; PTR_ADD A_cur, A_cur, LDA

    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 0; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 1; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 2; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 3

    vldrepl.w   Z0, X0, 0x00; vldrepl.w Z1, B1, 0x00
    vfmadd.s    D0, S0, Z0, D0; vfmadd.s D1, S1, Z0, D1
    vfmadd.s    D2, S0, Z1, D2; vfmadd.s D3, S1, Z1, D3
    PTR_ADDI    X0,     X0,     0x04
    PTR_ADDI    B1,     B1,     0x04
    PTR_ADDI    A0,     A0,     0x04
    bnez        K1,     .L_M8_N2_K1
.L_M8_N2_END:
    vfmul.s     D0,  D0,  VALPHA; vfmul.s D1,  D1,  VALPHA; vfmul.s D2,  D2,  VALPHA; vfmul.s D3,  D3,  VALPHA
#ifndef B0
    vld         S0, C0, 0x00; vld S1, C0, 0x10
    vfmadd.s    D0, S0, VBETA, D0; vfmadd.s D1, S1, VBETA, D1

    vld         S0, C1, 0x00; vld S1, C1, 0x10
    vfmadd.s    D2, S0, VBETA, D2; vfmadd.s D3, S1, VBETA, D3
#endif
    vst         D0,  C0, 0x00; vst D1,  C0, 0x10
    vst         D2,  C1, 0x00; vst D3,  C1, 0x10
    PTR_ALSL    C0,     LDC,    C0,     1
    PTR_ALSL    C1,     LDC,    C1,     1
    PTR_SUB     X0,     X0,    K4
    PTR_SUB     B1,     B1,    K4
    PTR_ALSL    X0,     LDB,    X0,     1
    PTR_ALSL    B1,     LDB,    B1,     1
    move        A0,     A
.L_M8_N1:
    andi        N1,     N,      0x01
    beqz        N1,     .L_M8_END
    vxor.v      D0, D0, D0; vxor.v D1, D1, D1
    move        K1,     K
    bge         ZERO,   K,      .L_M8_N1_END
.L_M8_N1_K1:
    PTR_ADDI    K1,     K1,    -1
    move        A_cur,  A0
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 0; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 1; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 2; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S0, TEMP, 3; PTR_ADD A_cur, A_cur, LDA

    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 0; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 1; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 2; PTR_ADD A_cur, A_cur, LDA
    ld.w TEMP, A_cur, 0; vinsgr2vr.w S1, TEMP, 3

    vldrepl.w   Z0, X0, 0x00
    vfmadd.s    D0, S0, Z0, D0; vfmadd.s D1, S1, Z0, D1
    PTR_ADDI    X0,     X0,     0x04
    PTR_ADDI    A0,     A0,     0x04
    bnez        K1,     .L_M8_N1_K1
.L_M8_N1_END:
    vfmul.s     D0,  D0,  VALPHA; vfmul.s D1,  D1,  VALPHA
#ifndef B0
    vld         S0, C0, 0x00; vld S1, C0, 0x10
    vfmadd.s    D0, S0, VBETA, D0; vfmadd.s D1, S1, VBETA, D1
#endif
    vst         D0, C0, 0x00; vst D1, C0, 0x10
.L_M8_END:
    PTR_SLLI    TEMP,   LDA,    3
    PTR_ADD     A,      A,      TEMP
    PTR_ADDI    C,      C,     0x20
.L_M7:
    andi        M4,     M,     0x04
    beqz        M4,     .L_M3
.L_M4:
    SGEMM_SMALL_KERNEL_TN_TAIL 4
.L_M3:
    andi        M2,     M,     0x02
    beqz        M2,     .L_M1
.L_M2:
    SGEMM_SMALL_KERNEL_TN_TAIL 2
.L_M1:
    andi        M1,     M,     0x01
    beqz        M1,     .L_M0
    SGEMM_SMALL_KERNEL_TN_TAIL 1
.L_M0:
    pop_if_used  7, 2
    jirl    $r0, $r1, 0x0
    EPILOGUE