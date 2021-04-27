/***************************************************************************
Copyright (c) 2020, The OpenBLAS Project
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

#define HAVE_KERNEL 1

static void copy_kernel (BLASLONG n, FLOAT *x, FLOAT *y)
{
  __asm__
    (
       "lxvp		32, 0(%2)	\n\t"
       "lxvp		34, 32(%2)	\n\t"
       "lxvp		36, 64(%2)	\n\t"
       "lxvp		38, 96(%2)	\n\t"
       "lxvp		40, 128(%2)	\n\t"
       "lxvp		42, 160(%2)	\n\t"
       "lxvp		44, 192(%2)	\n\t"
       "lxvp		46, 224(%2)	\n\t"

       "addi		%2, %2, 256	\n\t"
       "addic.		%1, %1, -32	\n\t"
       "ble		two%=		\n\t"

       ".align	5		\n"
     "one%=:				\n\t"

       "stxv		33, 0(%3)	\n\t"
       "stxv		32, 16(%3)	\n\t"
       "stxv		35, 32(%3)	\n\t"
       "stxv		34, 48(%3)	\n\t"
       "stxv		37, 64(%3)	\n\t"
       "stxv		36, 80(%3)	\n\t"
       "stxv		39, 96(%3)	\n\t"
       "stxv		38, 112(%3)	\n\t"
       "lxvp		32, 0(%2)	\n\t"
       "lxvp		34, 32(%2)	\n\t"
       "lxvp		36, 64(%2)	\n\t"
       "lxvp		38, 96(%2)	\n\t"

       "stxv		41, 128(%3)	\n\t"
       "stxv		40, 144(%3)	\n\t"
       "stxv		43, 160(%3)	\n\t"
       "stxv		42, 176(%3)	\n\t"
       "stxv		45, 192(%3)	\n\t"
       "stxv		44, 208(%3)	\n\t"
       "stxv		47, 224(%3)	\n\t"
       "stxv		46, 240(%3)	\n\t"
       "lxvp		40, 128(%2)	\n\t"
       "lxvp		42, 160(%2)	\n\t"
       "lxvp		44, 192(%2)	\n\t"
       "lxvp		46, 224(%2)	\n\t"


       "addi		%3, %3, 256	\n\t"
       "addi		%2, %2, 256	\n\t"

       "addic.		%1, %1, -32	\n\t"
       "bgt		one%=		\n"

     "two%=:				\n\t"

       "stxv		33, 0(%3)	\n\t"
       "stxv		32, 16(%3)	\n\t"
       "stxv		35, 32(%3)	\n\t"
       "stxv		34, 48(%3)	\n\t"
       "stxv		37, 64(%3)	\n\t"
       "stxv		36, 80(%3)	\n\t"
       "stxv		39, 96(%3)	\n\t"
       "stxv		38, 112(%3)	\n\t"
       "stxv		41, 128(%3)	\n\t"
       "stxv		40, 144(%3)	\n\t"
       "stxv		43, 160(%3)	\n\t"
       "stxv		42, 176(%3)	\n\t"
       "stxv		45, 192(%3)	\n\t"
       "stxv		44, 208(%3)	\n\t"
       "stxv		47, 224(%3)	\n\t"
       "stxv		46, 240(%3)	\n\t"

     "#n=%1 x=%4=%2 y=%0=%3"
     :
       "=m" (*y),
       "+r" (n),	// 1
       "+b" (x),	// 2
       "+b" (y) 	// 3
     :
       "m" (*x)
     :
       "cr0",
       "vs32","vs33","vs34","vs35","vs36","vs37","vs38","vs39",
       "vs40","vs41","vs42","vs43","vs44","vs45","vs46","vs47"
     );
}
