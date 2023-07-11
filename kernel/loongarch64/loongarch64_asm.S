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
#define PTR_ALSL  alsl.d
#else
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
#define PTR_ALSL  alsl.w
#endif

#if __loongarch_frlen == 64
#define FREG_SIZE 8
#define FREG_LOG  3
#define PTR_FLD   fld.d
#define PTR_FST   fst.d
#else
#define FREG_SIZE 4
#define FREG_LOG  2
#define PTR_FLD   fld.s
#define PTR_FST   fst.s
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
