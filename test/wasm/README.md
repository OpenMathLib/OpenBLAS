# WASM numerical validation suite

Deep CBLAS correctness checks for `TARGET=WASM128_GENERIC`, run under Node / Emscripten.

## Oracle

Results from OpenBLAS (public CBLAS API) are compared to a **hand-written scalar C reference** in `ref.c` (IEEE `*` / `+` only — no SIMD, no FMA intrinsics). This is not Netlib BLAS and not a second OpenBLAS build.

Tolerances live in `tol.h`. Builds with `WASM_RELAXED_SIMD=1` use a slightly larger L2/L3 budget (`TEST_WASM_RELAXED`).

## Run

```bash
JOBS=20 ./test/wasm/run.sh
```

This:

1. Builds OpenBLAS wasm with `WASM_RELAXED_SIMD=0`, links and runs the suite (IEEE tolerances).
2. Rebuilds with `WASM_RELAXED_SIMD=1`, links and runs again (relaxed tolerances).

Requires `emcc` on `PATH`, or an emscripten-forge prefix via `OPENBLAS_EM_PREFIX` / auto-discovery used by `benchmark/wasm/build.sh`.

## Coverage (MVP)

- L1: `saxpy` / `daxpy` (unit stride, non-unit, `inc==0`)
- L2: `sgemv` / `dgemv` (N/T, square and rectangular, some non-unit strides)
- L3: `sgemm` / `dgemm` / `cgemm` / `zgemm`, `ssyrk` / `dsyrk`, `strmm` / `dtrmm`, `strsm` / `dtrsm`

Size grids emphasize tile remainders around 4×4 / 8×4 / 2×2 (see `cases.h`).

Netlib `ctest` / `utest` remain a separate light gate (`benchmark/wasm/test.sh`).
