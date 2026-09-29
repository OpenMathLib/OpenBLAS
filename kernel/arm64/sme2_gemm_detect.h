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
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

/* Runtime check for the SME2 GEMM of sme2_gemm_impl.h; also used to keep the SME1 direct sgemm out of its way. */
#ifndef SME2_GEMM_DETECT_H
#define SME2_GEMM_DETECT_H

#include <arm_sme.h>
#if defined(__APPLE__)
#include <sys/sysctl.h>
#elif defined(__linux__)
#include <sys/auxv.h>
#endif

/* SME2 with a 512-bit streaming vector length (the kernels assume 16 fp32 lanes). */
static inline int s2_usable(void) {
  static int ok = -1;
  if (ok < 0) {
    int sme2 = 0;
#if defined(__APPLE__)
    int v = 0;
    size_t len = sizeof(v);
    sme2 = sysctlbyname("hw.optional.arm.FEAT_SME2", &v, &len, NULL, 0) == 0 && v;
#elif defined(__linux__)
    sme2 = (getauxval(AT_HWCAP2) & (1UL << 37)) != 0; /* HWCAP2_SME2 */
#endif
    ok = sme2 && svcntsw() == 16;
  }
  return ok;
}

/* The kernel keeps M, N and K in int; larger problems (INTERFACE64) stay on the other kernels. */
#define S2_FITS(m, n, k) ((m) < (1L << 30) && (n) < (1L << 30) && (k) < (1L << 30))

#endif
