/*****************************************************************************
Copyright (c) 2011-2020, The OpenBLAS Project
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
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********************************************************************************/

/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/*   1. Redistributions of source code must retain the above         */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer.                                                  */
/*                                                                   */
/*   2. Redistributions in binary form must reproduce the above      */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer in the documentation and/or other materials       */
/*      provided with the distribution.                              */
/*                                                                   */
/*    THIS  SOFTWARE IS PROVIDED  BY THE  UNIVERSITY OF  TEXAS AT    */
/*    AUSTIN  ``AS IS''  AND ANY  EXPRESS OR  IMPLIED WARRANTIES,    */
/*    INCLUDING, BUT  NOT LIMITED  TO, THE IMPLIED  WARRANTIES OF    */
/*    MERCHANTABILITY  AND FITNESS FOR  A PARTICULAR  PURPOSE ARE    */
/*    DISCLAIMED.  IN  NO EVENT SHALL THE UNIVERSITY  OF TEXAS AT    */
/*    AUSTIN OR CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT,    */
/*    INCIDENTAL,  SPECIAL, EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES    */
/*    (INCLUDING, BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE    */
/*    GOODS  OR  SERVICES; LOSS  OF  USE,  DATA,  OR PROFITS;  OR    */
/*    BUSINESS INTERRUPTION) HOWEVER CAUSED  AND ON ANY THEORY OF    */
/*    LIABILITY, WHETHER  IN CONTRACT, STRICT  LIABILITY, OR TORT    */
/*    (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY WAY OUT    */
/*    OF  THE  USE OF  THIS  SOFTWARE,  EVEN  IF ADVISED  OF  THE    */
/*    POSSIBILITY OF SUCH DAMAGE.                                    */
/*                                                                   */
/* The views and conclusions contained in the software and           */
/* documentation are those of the authors and should not be          */
/* interpreted as representing official policies, either expressed   */
/* or implied, of The University of Texas at Austin.                 */
/*********************************************************************/

#ifndef COMMON_LOONGARCH64
#define COMMON_LOONGARCH64

#define MB  __sync_synchronize()
#define WMB __sync_synchronize()
#define RMB __sync_synchronize()

#ifndef ASSEMBLER

static inline int blas_quickdivide(blasint x, blasint y){
  return x / y;
}

#ifndef NO_AFFINITY
static inline int WhereAmI(void){
  int ret = 0, counter = 0;
  __asm__ volatile (
    "rdtimel.w  %[counter],   %[id]"
    : [id]"=r"(ret), [counter]"=r"(counter)
    :
    : "memory"
  );
  return ret;
}
#endif

#ifdef DOUBLE
#define GET_IMAGE(res)  __asm__ __volatile__("fmov.d %0, $f2" : "=f"(res)  : : "memory")
#else
#define GET_IMAGE(res)  __asm__ __volatile__("fmov.s %0, $f2" : "=f"(res)  : : "memory")
#endif

#define GET_IMAGE_CANCEL

#else

#ifdef DOUBLE
#define LD      fld.d
#define ST      fst.d
#define MADD    fmadd.d
#define NMADD   fnmadd.d
#define MSUB    fmsub.d
#define NMSUB   fnmsub.d
#define ADD     fadd.d
#define SUB     fsub.d
#define MUL     fmul.d
#define MOV     fmov.d
#define CMOVT   fsel
#define MTC     movgr2fr.d
#define FABS    fabs.d
#define CMPEQ   fcmp.ceq.d
#define CMPLE   fcmp.cle.d
#define CMPLT   fcmp.clt.d
#define NEG     fneg.d
#else
#define LD      fld.s
#define ST      fst.s
#define MADD    fmadd.s
#define NMADD   fnmadd.s
#define MSUB    fmsub.s
#define NMSUB   fnmsub.s
#define ADD     fadd.s
#define SUB     fsub.s
#define MUL     fmul.s
#define MOV     fmov.s
#define CMOVT   fsel
#define MTC     movgr2fr.w
#define FABS    fabs.s
#define CMPEQ   fcmp.ceq.s
#define CMPLE   fcmp.cle.s
#define CMPLT   fcmp.clt.s
#define NEG     fneg.s
#endif /* defined(DOUBLE) */

#if defined(__64BIT__) && defined(USE64BITINT)
#define LDINT   ld.d
#define LDARG   ld.d
#define SDARG   st.d
#elif defined(__64BIT__) && !defined(USE64BITINT)
#define LDINT   ld.w
#define LDARG   ld.d
#define SDARG   st.d
#else
#define LDINT   ld.w
#define LDARG   ld.w
#define SDARG   st.w
#endif


#ifndef F_INTERFACE
#define REALNAME ASMNAME
#else
#define REALNAME ASMFNAME
#endif /* defined(F_INTERFACE) */

#if defined(ASSEMBLER) && !defined(NEEDPARAM)

#define PROLOGUE \
    .text ;\
    .align 5 ;\
    .globl  REALNAME ;\
    .type   REALNAME, @function ;\
REALNAME: ;\

#if defined(__linux__) && defined(__ELF__)
#define GNUSTACK .section .note.GNU-stack,"",@progbits
#else
#define GNUSTACK
#endif /* defined(__linux__) && defined(__ELF__) */

#define EPILOGUE      \
    .end    REALNAME ;\
    GNUSTACK

#define PROFCODE

#define MOVT(dst, src, cc)  \
    bceqz  cc,   1f;        \
    add.d dst,  src,  $r0;  \
    1:

#endif /* defined(ASSEMBLER) && !defined(NEEDPARAM) */

#endif /* defined(ASSEMBLER) */

#define SEEK_ADDRESS

#define BUFFER_SIZE     ( 32 << 20)

#define PAGESIZE        (16UL << 10)
#define FIXED_PAGESIZE  (16UL << 10)
#define HUGE_PAGESIZE   ( 2 << 20)

#define BASE_ADDRESS (START_ADDRESS - BUFFER_SIZE * MAX_CPU_NUMBER)

#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif

#endif
