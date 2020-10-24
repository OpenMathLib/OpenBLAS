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

       "lxvp		48, 256(%2)	\n\t"
       "lxvp		50, 288(%2)	\n\t"
       "lxvp		52, 320(%2)	\n\t"
       "lxvp		54, 352(%2)	\n\t"
       "lxvp		56, 384(%2)	\n\t"
       "lxvp		58, 416(%2)	\n\t"
       "lxvp		60, 448(%2)	\n\t"
       "lxvp		62, 480(%2)	\n\t"
       "addi		%2, %2, 512	\n\t"
#if !defined(COMPLEX) && !defined(DOUBLE)
       "addic.		%1, %1, -128	\n\t"
#elif defined(COMPLEX) && defined(DOUBLE)
       "addic.		%1, %1, -32	\n\t"
#else
       "addic.		%1, %1, -64	\n\t"
#endif
       "ble		two%=		\n\t"

       ".align	5		\n"
     "one%=:				\n\t"

       "stxvp		32, 0(%3)	\n\t"
       "lxvp		32, 0(%2)	\n\t"
       "stxvp		34, 32(%3)	\n\t"
       "lxvp		34, 32(%2)	\n\t"
       "stxvp		36, 64(%3)	\n\t"
       "lxvp		36, 64(%2)	\n\t"
       "stxvp		38, 96(%3)	\n\t"
       "lxvp		38, 96(%2)	\n\t"

       "stxvp		40, 128(%3)	\n\t"
       "lxvp		40, 128(%2)	\n\t"
       "stxvp		42, 160(%3)	\n\t"
       "lxvp		42, 160(%2)	\n\t"
       "stxvp		44, 192(%3)	\n\t"
       "lxvp		44, 192(%2)	\n\t"
       "stxvp		46, 224(%3)	\n\t"
       "lxvp		46, 224(%2)	\n\t"

       "stxvp		48, 256(%3)	\n\t"
       "lxvp		48, 256(%2)	\n\t"
       "stxvp		50, 288(%3)	\n\t"
       "lxvp		50, 288(%2)	\n\t"
       "stxvp		52, 320(%3)	\n\t"
       "lxvp		52, 320(%2)	\n\t"
       "stxvp		54, 352(%3)	\n\t"
       "lxvp		54, 352(%2)	\n\t"
       "stxvp		56, 384(%3)	\n\t"
       "lxvp		56, 384(%2)	\n\t"
       "stxvp		58, 416(%3)	\n\t"
       "lxvp		58, 416(%2)	\n\t"
       "stxvp		60, 448(%3)	\n\t"
       "lxvp		60, 448(%2)	\n\t"
       "stxvp		62, 480(%3)	\n\t"
       "lxvp		62, 480(%2)	\n\t"

       "addi		%3, %3, 512	\n\t"
       "addi		%2, %2, 512	\n\t"

#if !defined(COMPLEX) && !defined(DOUBLE)
       "addic.		%1, %1, -128	\n\t"
#elif defined(COMPLEX) && defined(DOUBLE)
       "addic.		%1, %1, -32	\n\t"
#else
       "addic.		%1, %1, -64	\n\t"
#endif
       "bgt		one%=		\n"

     "two%=:				\n\t"

       "stxvp		32, 0(%3)	\n\t"
       "stxvp		34, 32(%3)	\n\t"
       "stxvp		36, 64(%3)	\n\t"
       "stxvp		38, 96(%3)	\n\t"
       "stxvp		40, 128(%3)	\n\t"
       "stxvp		42, 160(%3)	\n\t"
       "stxvp		44, 192(%3)	\n\t"
       "stxvp		46, 224(%3)	\n\t"
       "stxvp		48, 256(%3)	\n\t"
       "stxvp		50, 288(%3)	\n\t"
       "stxvp		52, 320(%3)	\n\t"
       "stxvp		54, 352(%3)	\n\t"
       "stxvp		56, 384(%3)	\n\t"
       "stxvp		58, 416(%3)	\n\t"
       "stxvp		60, 448(%3)	\n\t"
       "stxvp		62, 480(%3)	\n\t"

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
       "vs40","vs41","vs42","vs43","vs44","vs45","vs46","vs47",
       "vs48","vs49","vs50","vs51","vs52","vs53","vs54","vs55",
       "vs56","vs57","vs58","vs59","vs60","vs61","vs62","vs63"
     );
}
