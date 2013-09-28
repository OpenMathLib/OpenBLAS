/*****************************************************************************
Copyright (c) 2011, Lab of Parallel Software and Computational Science,ICSAS
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
   3. Neither the name of the ISCAS nor the names of its contributors may 
      be used to endorse or promote products derived from this software 
      without specific prior written permission.

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

#ifndef COMMON_ARM
#define COMMON_ARM

#define MB
#define WMB

#define INLINE inline

#define RETURN_BY_COMPLEX

#ifndef ASSEMBLER

static void INLINE blas_lock(volatile unsigned long *address){

//  long int ret, val = 1;
/*
  do {
    while (*address) {YIELDING;};

    __asm__ __volatile__(
			 "1:	ll	%0, %3\n"
			 "	ori	%2, %0, 1\n"
			 "	sc	%2, %1\n"
			 "	beqz	%2, 1b\n"
			 "	 andi	%2, %0, 1\n"
			 "	sync\n"
			 : "=&r" (val), "=m" (address), "=&r" (ret)
			 : "m" (address)
			 : "memory");

  } while (ret);
*/
}

static inline unsigned int rpcc(void){
  unsigned long ret=0;

  return ret;
}

static inline int blas_quickdivide(blasint x, blasint y){
  return x / y;
}

#if defined(DOUBLE)
#define GET_IMAGE(res)  __asm__ __volatile__("vstr.f64 d1, %0" : "=m"(res) : : "memory")
#else
#define GET_IMAGE(res)  __asm__ __volatile__("vstr.f32 s1, %0" : "=m"(res) : : "memory")
#endif

#define GET_IMAGE_CANCEL

#endif


#ifndef F_INTERFACE
#define REALNAME ASMNAME
#else
#define REALNAME ASMFNAME
#endif

#if defined(ASSEMBLER) && !defined(NEEDPARAM)

#define PROLOGUE \
	.arm		 ;\
	.global	REALNAME ;\
	.func	REALNAME  ;\
REALNAME:

#define EPILOGUE 

#define PROFCODE

#endif


#define SEEK_ADDRESS

#ifndef PAGESIZE
#define PAGESIZE        ( 4 << 10)
#endif
#define HUGE_PAGESIZE   ( 4 << 20)

#define BUFFER_SIZE     (16 << 20)


#define BASE_ADDRESS (START_ADDRESS - BUFFER_SIZE * MAX_CPU_NUMBER)

#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif

#endif
