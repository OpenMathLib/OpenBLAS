/*******************************************************************************
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
*******************************************************************************/

#if __loongarch_grlen == 64
#define LA_REG    int64_t
#define REG_SIZE  8
#define REG_LOG   3
#define PTR_ADDI  addi.d
#define PTR_ADD   add.d
#define PTR_SUB   sub.d
#define PTR_LD    ld.d
#define PTR_ST    st.d
#define PTR_SLLI  slli.d
#define PTR_SRLI  srli.d
#define PTR_SRAI  srai.d
#define PTR_MUL   mul.d
#define PTR_ALSL  alsl.d
#elif __loongarch_grlen == 32
#define LA_REG    int32_t
#define REG_SIZE  4
#define REG_LOG   2
#define PTR_ADDI  addi.w
#define PTR_ADD   add.w
#define PTR_SUB   sub.w
#define PTR_LD    ld.w
#define PTR_ST    st.w
#define PTR_SLLI  slli.w
#define PTR_SRLI  srli.w
#define PTR_SRAI  srai.w
#define PTR_MUL   mul.w
#define PTR_ALSL  alsl.w
#else
// If neither of the above two conditions is supported, it means this is an early
// internal toolchain. To ensure maximum compatibility, the following approach is taken:
#define LA_REG    int64_t
#define REG_SIZE  8
#define REG_LOG   3
#define PTR_ADDI  addi.d
#define PTR_ADD   add.d
#define PTR_SUB   sub.d
#define PTR_LD    ld.d
#define PTR_ST    st.d
#define PTR_SLLI  slli.d
#define PTR_SRLI  srli.d
#define PTR_SRAI  srai.d
#define PTR_MUL   mul.d
#define PTR_ALSL  alsl.d
#endif

#if __loongarch_frlen == 64
#define FREG_SIZE 8
#define FREG_LOG  3
#define PTR_FLD   fld.d
#define PTR_FST   fst.d
#elif __loongarch_frlen == 32
#define FREG_SIZE 4
#define FREG_LOG  2
#define PTR_FLD   fld.s
#define PTR_FST   fst.s
#else
// If neither of the above two conditions is supported, it means this is an early
// internal toolchain. To ensure maximum compatibility, the following approach is taken:
#define FREG_SIZE 8
#define FREG_LOG  3
#define PTR_FLD   fld.d
#define PTR_FST   fst.d
#endif

// The max registers available to the user which
// do not need to be preserved across calls.
// Ref: https://loongson.github.io/LoongArch-Documentation/LoongArch-ELF-ABI-CN.html
#define MAX_INT_CALLER_SAVED 17
#define MAX_FP_CALLER_SAVED  24

.altmacro // Enable alternate macro mode

.macro push_if_used regs, fregs
.if \regs > MAX_INT_CALLER_SAVED
    PTR_ADDI      $sp,    $sp,    -((\regs - MAX_INT_CALLER_SAVED) << REG_LOG)
    push_regs 0, \regs - MAX_INT_CALLER_SAVED - 1
.endif
.if \fregs > MAX_FP_CALLER_SAVED
    PTR_ADDI      $sp,    $sp,    -((\fregs - MAX_FP_CALLER_SAVED) << FREG_LOG)
    push_fregs 0, \fregs - MAX_FP_CALLER_SAVED - 1
.endif
.endm // End push_if_used
.macro pop_if_used regs, fregs
.if \fregs > MAX_FP_CALLER_SAVED
    pop_fregs 0, \fregs - MAX_FP_CALLER_SAVED - 1
    PTR_ADDI      $sp,    $sp,    (\fregs - MAX_FP_CALLER_SAVED) << FREG_LOG
.endif
.if \regs > MAX_INT_CALLER_SAVED
    pop_regs 0, \regs - MAX_INT_CALLER_SAVED - 1
    PTR_ADDI      $sp,    $sp,    (\regs - MAX_INT_CALLER_SAVED) << REG_LOG
.endif
.endm // End pop_if_used
.macro push_regs from, to
    PTR_ST    $s\()\from,     $sp,    \from << REG_LOG
.if \to - \from
    push_regs %from + 1, \to
.endif
.endm // End push_regs
.macro pop_regs from, to
    PTR_LD    $s\()\from,     $sp,    \from << REG_LOG
.if \to - \from
    pop_regs %from + 1, \to
.endif
.endm // End pop_regs
.macro push_fregs from, to
    PTR_FST   $fs\()\from,    $sp,    \from << FREG_LOG
.if \to - \from
    push_fregs %from + 1, \to
.endif
.endm // End push_fregs
.macro pop_fregs from, to
    PTR_FLD   $fs\()\from,    $sp,    \from << FREG_LOG
.if \to - \from
    pop_fregs %from + 1, \to
.endif
.endm // End pop_fregs

//
// Instruction Related Macros
//
// GLD
//
.macro GLD pre_op:req, suf_op=0, out:req, src:req, offset:req/* imm */, more:vararg
.ifeqs "\suf_op", "0"
    \pre_op\()ld    \out,   \src,   \offset
.else
    \pre_op\()ld.\suf_op    \out,   \src,   \offset
.endif
.ifnb \more
    GLD \pre_op, \suf_op, \more
.endif
.endm

//
// GLD_INC
//
.macro GLD_INC pre_op:req, suf_op=0, inc:req, out:req, src:req, offset:req/* imm */, more:vararg
.ifeqs "\suf_op", "0"
    \pre_op\()ld    \out,   \src,   \offset
.else
    \pre_op\()ld.\suf_op    \out,   \src,   \offset
.endif
    PTR_ADDI  \src,   \src,   \inc
.ifnb \more
    GLD_INC \pre_op, \suf_op, \inc, \more
.endif
.endm
//
// GLDX is same as GLD except the stride is a register
//
.macro GLDX pre_op:req, suf_op=0, out:req, src:req, offset:req/* reg */, more:vararg
.ifeqs "\suf_op", "0"
    \pre_op\()ldx    \out,   \src,   \offset
.else
    \pre_op\()ldx.\suf_op    \out,   \src,   \offset
.endif
.ifnb \more
    GLDX \pre_op, \suf_op, \more
.endif
.endm
//
// GLDREPL
//
.macro GLDREPL pre_op:req, suf_op:req, out:req, src:req, offset:req/* imm */, more:vararg
    \pre_op\()ldrepl.\suf_op    \out,   \src,   \offset
.ifnb \more
    GLDREPL  \pre_op, \suf_op, \more
.endif
.endm
//
// GST
//
.macro GST pre_op:req, suf_op=0, src:req, dst:req, offset:req/* imm */, more:vararg
.ifeqs "\suf_op", "0"
    \pre_op\()st    \src,   \dst,   \offset
.else
    \pre_op\()st.\suf_op \src,  \dst, \offset
.endif
.ifnb \more
    GST \pre_op, \suf_op, \more
.endif
.endm
//
// GMUL
//
.macro GMUL pre_op, suf_op:req, out:req, in0:req, in1:req, more:vararg
    \pre_op\()mul.\suf_op   \out,   \in0,   \in1
.ifnb \more
    GMUL \pre_op, \suf_op, \more
.endif
.endm
//
// GMADD
//
.macro GMADD pre_op, suf_op:req, out:req, in0:req, in1:req, in2:req, more:vararg
    \pre_op\()madd.\suf_op \out, \in0, \in1, \in2
.ifnb \more
    GMADD \pre_op, \suf_op, \more
.endif
.endm
//
// GADD
//
.macro GADD pre_op, suf_op:req, out:req, in0:req, in1:req, more:vararg
    \pre_op\()add.\suf_op \out, \in0, \in1
.ifnb \more
    GADD \pre_op, \suf_op, \more
.endif
.endm
//
// GADDI
//
.macro GADDI pre_op, suf_op:req, out:req, in0:req, in1:req, more:vararg
    \pre_op\()addi.\suf_op  \out,   \in0,   \in1
.ifnb \more
    GADDI \pre_op, \suf_op, \more
.endif
.endm
//
// GSUB
//
.macro GSUB pre_op, suf_op:req, out:req, in0:req, in1:req, more:vararg
    \pre_op\()sub.\suf_op \out, \in0, \in1
.ifnb \more
    GSUB \pre_op, \suf_op, \more
.endif
.endm
//
// GSLLI
//
.macro GSLLI pre_op, suf_op:req, out:req, in0:req, in1:req, more:vararg
    \pre_op\()slli.\suf_op  \out,   \in0,   \in1
.ifnb \more
    GSLLI \pre_op, \suf_op, \more
.endif
.endm
//
// GINSVE0
//
.macro GINSVE0 pre_op:req, suf_op:req, out:req, in0:req, in1:req, more:vararg
    \pre_op\()insve0.\suf_op    \out,   \in0,   \in1
.ifnb \more
    GINSVE0 \pre_op, \suf_op, \more
.endif
.endm
//
// GXOR
//
.macro GXOR pre_op:req, suf_op:req, out:req, in0:req, in1:req, more:vararg
    \pre_op\()xor.\suf_op    \out,   \in0,   \in1
.ifnb \more
    GXOR \pre_op, \suf_op, \more
.endif
.endm
//
// GPERMI
//
.macro GPERMI pre_op:req, suf_op:req, out:req, in0:req, in1:req, more:vararg
    \pre_op\()permi.\suf_op     \out,   \in0,   \in1
.ifnb \more
    GPERMI \pre_op, \suf_op, \more
.endif
.endm
//
// GNMSUB
//
.macro GNMSUB pre_op:req, suf_op:req, out:req, in0:req, in1:req, in2:req, more:vararg
    \pre_op\()nmsub.\suf_op     \out,   \in0,   \in1,   \in2
.ifnb \more
    GNMSUB \pre_op, \suf_op, \more
.endif
.endm
//
// GPRELD
//
.macro GPRELD in0:req, in1:req, in2:req, more:vararg
    preld \in0, \in1, \in2
.ifnb \more
    GPRELD \more
.endif
.endm

//
// Compound instructions
//
// GACC: Accumulate the values of vector registers
//
.macro GACC pre_op:req, suf_op:req, out:req, in:req, more:vararg
.ifeqs "\pre_op", "xvf"
    xvpermi.q              \out,   \in,    0x01
    \pre_op\()add.\suf_op  \in,    \out,   \in
    xvpackod.d             \out,   \in,    \in
    \pre_op\()add.\suf_op  \out,   \out,   \in
.ifeqs "\suf_op", "s"
    xvpackod.w             \in,    \out,   \out
    \pre_op\()add.\suf_op  \out,   \out,   \in
.endif
.endif

.ifeqs "\pre_op", "vf"
    vpackod.d              \out,   \in,    \in
    \pre_op\()add.\suf_op  \out,   \out,   \in
.ifeqs "\suf_op", "s"
    vpackod.w              \in,    \out,   \out
    \pre_op\()add.\suf_op  \out,   \out,   \in
.endif
.endif

.ifeqs "\pre_op", "xv"
    xvpermi.q              \out,   \in,    0x01
    \pre_op\()add.\suf_op  \in,    \out,   \in
    xvpackod.d             \out,   \in,    \in
    \pre_op\()add.\suf_op  \out,   \out,   \in
.ifnc "\suf_op", "d"
    xvpackod.w             \in,    \out,   \out
    \pre_op\()add.\suf_op  \out,   \out,   \in
.ifnc "\suf_op", "w"
    xvpackod.h             \in,    \out,   \out
    \pre_op\()add.\suf_op  \out,   \out,   \in
.ifnc "\suf_op", "h"
    xvpackod.b             \in,    \out,   \out
    \pre_op\()add.\suf_op  \out,   \out,   \in
.endif
.endif
.endif
.endif

.ifeqs "\pre_op", "v"
    vpackod.d              \out,   \in,    \in
    \pre_op\()add.\suf_op  \out,   \out,   \in
.ifnc "\suf_op", "d"
    vpackod.w              \in,    \out,   \out
    \pre_op\()add.\suf_op  \out,   \out,   \in
.ifnc "\suf_op", "w"
    vpackod.h              \in,    \out,   \out
    \pre_op\()add.\suf_op  \out,   \out,   \in
.ifnc "\suf_op", "h"
    vpackod.b              \in,    \out,   \out
    \pre_op\()add.\suf_op  \out,   \out,   \in
.endif
.endif
.endif
.endif

.ifnb \more
    GACC \pre_op, \suf_op, \more
.endif
.endm
//
// GMOV
//
.macro GMOV pre_op:req, out:req, in:req, more:vararg
    \pre_op\()or.v  \out, \in, \in
.ifnb \more
    GMOV \pre_op, \more
.endif
.endm

//
// Media Related Macros
//
.macro GSBUTTERFLY pre_op, suf_op, out0, out1, in0, in1
    \pre_op\()ilvl.\suf_op \out0, \in0, \in1
    \pre_op\()ilvh.\suf_op \out1, \in0, \in1
.endm
.macro GINTERLACE pre_op, suf_op, out0, out1, in0, in1
    \pre_op\()pickev.\suf_op \out0, \in0, \in1
    \pre_op\()pickod.\suf_op \out1, \in0, \in1
.endm

//
// TRANSPOSE4x4_D: Transpose 4x4 block with double-word elements in vectors,
// has no pre_op param. 128-bit vector instructions are not supported.
//
.macro GTRANSPOSE4x4_D in0, in1, in2, in3, out0, out1, out2, out3, \
                       vt0, vt1
    GSBUTTERFLY xv, d, \vt0, \out1, \in1, \in0
    GSBUTTERFLY xv, d, \vt1, \out3, \in3, \in2
    GMOV xv, \out0, \vt0, \out2, \vt1, \vt1, \out3
    GPERMI xv, q, \out0, \out2, 0x02, \out2, \vt0, 0x31, \out3, \out1, 0x31, \out1, \vt1, 0x02
.endm

.macro GTRANSPOSE8x8_W out0, out1, out2, out3, out4, out5, out6, out7, \
                       in0, in1, in2, in3, in4, in5, in6, in7, \
                       tmp0, tmp1, tmp2, tmp3
    GSBUTTERFLY xv, w, \tmp0, \tmp2, \in2, \in0
    GSBUTTERFLY xv, w, \tmp1, \tmp3, \in3, \in1
    GSBUTTERFLY xv, w, \out0, \out1, \tmp1, \tmp0
    GSBUTTERFLY xv, w, \out2, \out3, \tmp3, \tmp2

    GSBUTTERFLY xv, w, \tmp0, \tmp2, \in6, \in4
    GSBUTTERFLY xv, w, \tmp1, \tmp3, \in7, \in5
    GSBUTTERFLY xv, w, \out4, \out5, \tmp1, \tmp0
    GSBUTTERFLY xv, w, \out6, \out7, \tmp3, \tmp2

    GMOV xv, \tmp0, \out0, \tmp1, \out1, \tmp2, \out2, \tmp3, \out3

    GPERMI xv, q, \out0, \out4, 0x02, \out1, \out5, 0x02, \
                  \out2, \out6, 0x02, \out3, \out7, 0x02, \
                  \out4, \tmp0, 0x31, \out5, \tmp1, 0x31, \
                  \out6, \tmp2, 0x31, \out7, \tmp3, 0x31
.endm
