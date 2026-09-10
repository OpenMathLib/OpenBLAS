# WASM numerical validation suite

Deep CBLAS correctness checks for `TARGET=WASM128_GENERIC`, run under Node / Emscripten.

## Oracle

Results from OpenBLAS (public CBLAS API) are compared to hand-written scalar
C references in `ref_l1.c`, `ref_l2.c`, and `ref_l3.c` (IEEE `*` / `+` only —
no SIMD or FMA intrinsics). This is not Netlib BLAS and not a second OpenBLAS
build.

Tolerances live in `tol.h`. Builds with `WASM_RELAXED_SIMD=1` use a slightly larger L2/L3 budget (`TEST_WASM_RELAXED`).

## Run

```bash
JOBS=20 ./test/wasm/run.sh
```

This:

1. Builds OpenBLAS wasm with `WASM_RELAXED_SIMD=0`, links and runs the suite (IEEE tolerances).
2. Rebuilds with `WASM_RELAXED_SIMD=1`, links and runs again (relaxed tolerances).

Requires `emcc` on `PATH`, or an emscripten-forge prefix via `OPENBLAS_EM_PREFIX` / auto-discovery used by `benchmark/wasm/build.sh`.

CI runs the same script via `.github/workflows/wasm.yml` (Emscripten + Node on `ubuntu-latest`) on changes under `test/wasm/` and `kernel/wasm/`.

## Coverage

The suite covers the complete standard CBLAS Level 1/2/3 families:

- Level 1: rotations, swap, scaling, copy, axpy, dot products, norms, absolute
  sums, and maximum-index operations for all applicable S/D/C/Z types.
- Level 2 dense: general, symmetric, Hermitian, triangular, and rank-update
  operations.
- Level 2 banded and packed: general, symmetric/Hermitian, triangular, solve,
  and rank-update operations.
- Level 3: GEMM, SYMM/HEMM, SYRK/HERK, SYR2K/HER2K, and TRMM/TRSM.

Dense vectors and matrices are filled by `fill_vec_*` / `fill_mat_*` in
`common.h` (existing `fill_f32` / `fill_f64` / `fill_c32` / `fill_c64` wrap
those). Drivers cycle six `FillSpec` cases across the size grid so each
remainder hits a different sign domain and magnitude spread:

| Domain | Spread | f32 magnitudes | f64 magnitudes |
| --- | --- | --- | --- |
| R⁺ (strictly positive) | near 0 | log-uniform in [1e-4, 1] | [1e-8, 1] |
| R⁺ | far from 0 | [1e2, 1e4] | [1e4, 1e8] |
| R⁻ (strictly negative) | near 0 / far | same magnitudes, negated | same |
| R \ {0} (mixed signs) | near 0 / far | same magnitudes, random sign | same |

Values are never exactly 0. Near-0 lower bounds stay large enough that
`tol * (1 + maxv)` still flags a wrong kernel; far-from-0 upper bounds stay
small enough that L3 GEMM at n≈129 does not overflow f32. Failure messages
include the active spec (e.g. `R+ far from 0`).

Triangular solve / multiply fixtures (`make_tri_*`) and band builders stay
O(1) and diagonally dominant so those problems remain well-conditioned.
Symmetric and Hermitian inputs are mirrored explicitly; band and packed
layouts include their off-diagonals. Triangular solves are checked by
constructing a right-hand side with the matching scalar matrix product and
recovering the original input.

This scope is standard BLAS only. OpenBLAS extensions such as `axpby`, `gemmt`,
`imatcopy`, and bfloat16 routines are intentionally excluded.

Size grids emphasize tile remainders around 4×4 / 8×4 / 2×2 (see `cases.h`).

Netlib `ctest` / `utest` remain a separate light gate (`benchmark/wasm/test.sh`).
