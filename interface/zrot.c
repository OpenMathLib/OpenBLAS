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

#include <stdio.h>
#include "common.h"
#ifdef FUNCTION_PROFILE
#include "functable.h"
#endif

#if defined(ARCH_X86_64) && !defined(XDOUBLE)
/* The real c and s act on the real and imaginary parts separately, so with unit
   increments the 2n reals can be rotated as one real vector. For long vectors the
   srot/drot kernels of Haswell and newer (AVX2/AVX-512 FMA) are about twice as fast
   as the complex SSE kernels; for short ones they cost more per call, and on the older
   cores the real SSE kernels are not faster at all. So this is only done for
   RROT_MIN <= n <= RROT_MAX: the upper limit keeps 2n within the range of n that the
   real kernels get from srot/drot. The test is marked unlikely so that the short
   vectors keep a straight jump to the complex kernel. */
#ifdef DOUBLE
#define RROT_K DROT_K
#else
#define RROT_K SROT_K
#endif
#define RROT_MIN 4096
#define RROT_MAX ((BLASLONG)1 << 29)
#if defined(__GNUC__)
#define RROT_USE(n) __builtin_expect((n) >= RROT_MIN && (n) <= RROT_MAX, 0)
#else
#define RROT_USE(n) ((n) >= RROT_MIN && (n) <= RROT_MAX)
#endif
#endif

OPENBLAS_EXPORT

#ifndef CBLAS
void NAME(blasint *N, FLOAT *x, blasint *INCX, FLOAT *y, blasint *INCY, FLOAT *C, FLOAT *S){
  BLASLONG n    = *N;
  BLASLONG incx = *INCX;
  BLASLONG incy = *INCY;
  FLOAT c = *C;
  FLOAT s = *S;

#else
void CNAME(blasint n, void *VX, blasint incx, void *VY, blasint incy, FLOAT c, FLOAT s) {
    FLOAT *x = (FLOAT*) VX;
    FLOAT *y = (FLOAT*) VY;
#endif /* CBLAS */

  PRINT_DEBUG_NAME;

  if (n <= 0) return;

  IDEBUG_START;

  FUNCTION_PROFILE_START();

  if (incx < 0) x -= (n - 1) * 2 * incx;
  if (incy < 0) y -= (n - 1) * 2 * incy;

#ifdef RROT_K
  if (RROT_USE(n) && incx == 1 && incy == 1)
    RROT_K(2 * n, x, 1, y, 1, c, s);
  else
    ROT_K(n, x, incx, y, incy, c, s);
#else
  ROT_K(n, x, incx, y, incy, c, s);
#endif

  FUNCTION_PROFILE_END(4, n, n);

  IDEBUG_END;

  return;

}
