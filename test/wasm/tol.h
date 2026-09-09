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

/*
 * Relative pass condition (see expect_close_* in common.h):
 *
 *   max_i |got[i] - ref[i]|  <=  tol * (1 + max_i |ref[i]|)
 *
 * So `tol` is a mixed absolute/relative budget: near zero it behaves like an
 * absolute tolerance; for O(1) results it is relative.
 *
 * Choice of numbers
 * -----------------
 * The scalar oracle uses only IEEE `*` / `+` (no FMA). OpenBLAS WASM kernels
 * may reorder sums and, under WASM_RELAXED_SIMD, use f32x4.relaxed_madd /
 * f64x2.relaxed_madd, which need not match IEEE fused multiply-add.
 *
 * Machine epsilons are ~1.2e-7 (f32) and ~2.2e-16 (f64). The budgets below
 * are intentionally larger than a few ulps so the suite gates kernel bugs
 * (wrong tiles, strides, remainders) rather than benign rounding differences.
 *
 * L1 (AXPY and friends): current WASM L1 paths stay ordinary mul+add even
 * with -mrelaxed-simd, so a fixed, tight budget is enough (~30 ulp at |ref|~1
 * for f32). No dependence on n: depth of accumulation is small.
 *
 * L2 / L3: error can grow with the reduction length (inner dimension ~ n for
 * the square cases we run). Budgets are therefore `SCALE * max(n, 1)`.
 * IEEE SCALE is ~1.7e3 ulp (f32) / ~9e3 ulp (f64) per unit of n — loose
 * enough for blocked GEMM/TRMM association, still far below what a wrong
 * microkernel typically produces. Relaxed SCALE is 4× IEEE to cover
 * relaxed_madd without hiding clear functional failures.
 */

/* L2/L3 scale factors; multiplied by problem size in tol_*_l2 / tol_*_l3. */
#ifdef TEST_WASM_RELAXED
#define TOL_S_SCALE 8e-4f
#define TOL_D_SCALE 8e-12
#else
#define TOL_S_SCALE 2e-4f
#define TOL_D_SCALE 2e-12
#endif

/* Fixed L1 budgets (not scaled by n). */
#define TOL_S_L1 4e-6f
#define TOL_D_L1 4e-14

static inline float tol_s_l3(int n) {
  return TOL_S_SCALE * (float)(n > 0 ? n : 1);
}
static inline double tol_d_l3(int n) {
  return TOL_D_SCALE * (double)(n > 0 ? n : 1);
}
static inline float tol_s_l2(int n) {
  return TOL_S_SCALE * (float)(n > 0 ? n : 1);
}
static inline double tol_d_l2(int n) {
  return TOL_D_SCALE * (double)(n > 0 ? n : 1);
}

#endif /* TEST_WASM_TOL_H */
