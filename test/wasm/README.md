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

Matrices used by the full checks include off-diagonal values. Symmetric and
Hermitian inputs are mirrored explicitly; triangular solve inputs are
diagonally dominant; band and packed layouts include their off-diagonals.
Triangular solves are checked by constructing a right-hand side with the
matching scalar matrix product and recovering the original input.

This scope is standard BLAS only. OpenBLAS extensions such as `axpby`, `gemmt`,
`imatcopy`, and bfloat16 routines are intentionally excluded.

Size grids emphasize tile remainders around 4×4 / 8×4 / 2×2 (see `cases.h`).

Netlib `ctest` / `utest` remain a separate light gate (`benchmark/wasm/test.sh`).
