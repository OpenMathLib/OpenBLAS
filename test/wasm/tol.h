/*
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
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef TEST_WASM_TOL_H
#define TEST_WASM_TOL_H

/* WASM_RELAXED_SIMD builds may use f32x4.relaxed_madd / f64x2.relaxed_madd
 * in GEMM/TRMM-class kernels. Allow a modestly larger budget then. */

#ifdef TEST_WASM_RELAXED
#define TOL_S_SCALE 8e-4f
#define TOL_D_SCALE 8e-12
#else
#define TOL_S_SCALE 2e-4f
#define TOL_D_SCALE 2e-12
#endif

/* L1 (AXPY) stays IEEE mul+add even with -mrelaxed-simd in current kernels. */
#define TOL_S_L1  4e-6f
#define TOL_D_L1  4e-14

static inline float tol_s_l3(int n) { return TOL_S_SCALE * (float)(n > 0 ? n : 1); }
static inline double tol_d_l3(int n) { return TOL_D_SCALE * (double)(n > 0 ? n : 1); }
static inline float tol_s_l2(int n) { return TOL_S_SCALE * (float)(n > 0 ? n : 1); }
static inline double tol_d_l2(int n) { return TOL_D_SCALE * (double)(n > 0 ? n : 1); }

#endif /* TEST_WASM_TOL_H */
