# RDNA2 (gfx1030 / gfx1031) and RDNA4 (gfx1200 / gfx1201) int8 kernels

Standalone HIP kernels for RDNA2 and RDNA4, written against the GFX10.3 and
GFX12 instruction sets. This is **not** part of the OpenBLAS CPU build and does
not replace any OpenBLAS routine. It is a separate, opt-in GPU component,
written directly in HIP C++: no CUDA-to-ROCm translation layer and no rocBLAS.

Two tile shapes are compiled from the same sources: the 64x64 / 4x4-8x4 shape
this tree started with, and a 128x128 / 8x8 shape aimed at RDNA4 (RX 9070 XT and
9070 are gfx1201, RX 9060 XT is gfx1200). Both run on both architectures. On the
same 128x128 tile the int8 GEMM can be driven either by 8x8 dp4a or by one
16x16x16 matrix instruction per k step, so there are three int8 shapes in all.
Which one the arch selectors pick is decided by measurement, not by generation:
dp4a is faster on the 64x64 tile than on the 128x128 one at every size tried,
and f32 is faster on the 64x64 tile until the two 128-row operand panels fall
out of cache, so RDNA2 keeps 64x64. RDNA4 takes the 128x128 tile only because
the matrix instruction available there beats 64x64 dp4a by 3.1-3.3x. Nothing is
exclusive to one generation: every shape is still compiled, still tested, and
can be requested explicitly on any device.

## Files

- `rdna2_kernels.hip`: the 17 device kernels and their host launchers.
- `rdna2_kernels_math.h`: the per-thread arithmetic of the GEMM/GEMV kernels, with no
  HIP types. The device kernels and the host test include this same header, so the
  exact tiling, bounds and accumulation code also runs on a CPU-only machine.
- `rdna2_gemm_ref.h`: plain C++ references (dp4a, gemv, gemm, sgemm, dot2, dpp
  prefix), used as the CPU fallback and as the test oracle.
- `test_rdna2_host.cpp`: host test, 141 checks against a device and 132 on the
  CPU fallback path (the kernel bodies and the board shapes are exercised on the
  host too, so it is useful without a GPU).
- `Makefile`: `make`, `make test`, `make multi`, `make clean`.

## Kernels

| Kernel | Notes |
|---|---|
| `rdna2_dp4a_gemv_i8` | int8 GEMV; on gfx1031 four packed int32 per 16-byte vector load when `k4 % 4 == 0`, scalar `v_dot4c` loads otherwise (a row of a non-multiple-of-4 `k4` is not 16-byte aligned). Any `k4` is valid |
| `rdna2_dp4a_gemv_i8_wave` | one wave per row; lanes stride the row so a wave reads 32 consecutive words (one 128-byte line) per step, then `__shfl_xor` trees the per-lane dp4a partials (lane 0 stores). Launcher uses it from `k4 >= 8` |
| `rdna2_sgemv_f32` | f32 GEMV, one thread per row; the row loop is `__builtin_fmaf` so the FMA contraction does not depend on `-ffast-math`. Launcher uses it below `K = 256` |
| `rdna2_sgemv_f32_wave` | f32 GEMV, one wave per row; lanes stride the row (consecutive floats per wave) and the partials are `__shfl_xor`-tree'd. Launcher uses it from `K >= RDNA2_SGEMV_WAVE_MIN_K` (256), because a short row does not give every lane enough work to pay for the reduction |
| `rdna2_dp4a_gemm_i8_naive` | one output element per thread, `v_dot4c_i32_i8` accumulation |
| `rdna2_dp4a_gemm_i8_tile` | 4x4 outputs per thread, 16x16 block (256 threads), 64x64 output tile |
| `rdna2_dp4a_gemm_i8_lds` | same tile and 4x4 blocking, K operands staged in LDS (`ds_read`/`ds_write`) instead of re-read from global |
| `rdna2_dp4a_gemm_i8_lds8` | same LDS staging, 8x4 blocking: 128 threads (16x8), 12 LDS reads per 32 dp4a instead of 8 per 16. **Fastest dp4a int8 shape measured, and the shape both RDNA2 arch selectors pick** |
| `rdna2_dp4a_gemm_i8_lds128` | 128x128 output tile with 8x8 blocking, 256 threads (16x16); 16 LDS reads per 64 dp4a, i.e. 4 multiplies per LDS read instead of 2.67 (8x4) or 2 (4x4). Less LDS traffic, but measured slower than `_lds8` on a 9070 XT |
| `rdna2_dp4a_gemm_i8_wmma` | the same 128x128 tile and 16x16 block, 8x8 dp4a replaced by one `v_wmma_i32_16x16x16_iu8` per k step: eight `int32x8_t` accumulators per lane, operands staged in the same 17408 B of LDS. **Default int8 shape on gfx1200/gfx1201.** Outside gfx12 the instruction does not exist, so the same kernel body runs a shuffle-based emulation there: correct, never selected by a dispatcher, and no timing claim |
| `rdna2_i8_row_sums` | per-row sum of the packed int8 activations: one signed `v_dot4_i32_i8` per word against an all-ones operand. Launched only when the format has a nonzero zero point; symmetric weights skip it |
| `rdna2_gemm_i8_epilogue` | applies the quantized format's scale and zero point to the finished int32 accumulators, `C = (acc - zp_b*rowsum_a[m]) * sa[m] * sb[n]`, with the correction in 64-bit. Converts one finished dot product per output element and nothing else |
| `rdna2_sgemm_f32` | f32 GEMM, same 4x4 / 16x16 / 64x64 blocking as the int8 tiled kernel |
| `rdna2_sgemm_f32_lds` | f32 GEMM with LDS-staged K chunks; accumulates through `__builtin_fmaf`, so the `v_fmac_f32` contraction does not depend on `-ffast-math` |
| `rdna2_sgemm_f32_lds128` | f32 counterpart of the 128x128 tile: same 8x8 blocking and 8-float K chunk as `rdna2_sgemm_f32_lds`, four times the tile area |
| `rdna2_sgemm_f32_lds16` | 64x64 tile with a 16-deep K chunk (`RDNA2_LDS_F32B_KSTEP`, stride 17, four floats staged per thread) instead of 8. **Default f32 shape on all four arch codes**: measured faster than the 8-deep chunk at every shape tried, the margin shrinking as the problem grows |
| `rdna2_dot2_f16` | f16 dot product via `__builtin_amdgcn_fdot2`, emitted as `v_dot2_f32_f16` |
| `rdna2_dp4a_dot_u8` | unsigned dot product via `v_dot4_u32_u8` |
| `rdna2_dpp_prefix` | intra-wave inclusive prefix sum using `v_add_nc_u32_dpp` + cross-16-lane fixup |

The naive and tiled int8 GEMMs are both kept: on a device they cross-check each
other, and on the host the test compares both against the reference. The same
holds for the three LDS variants, which differ only in tile shape and register
blocking, and for the matrix-instruction kernel, which runs the same tile as
`_lds128` through a different pipe. Nothing selects a shape from the kernel
side: the launcher's `tiled` / `variant` argument decides, and
`rdna2_gemm_i8_host_tuned` / `rdna2_sgemm_f32_host_tuned` pass the shape that
matches the device that is actually present. `tiled == 5` is the one case that
the host launcher always launches and the *kernel* decides on: `RDNA2_WMMA_IU8`
is a device-compilation constant, so it is 0 in host code, and gating the launch
on it there would silently select the dp4a kernel instead.

## Quantized data contract

The int8 entry points take packed int8 (`const int*` A and Bt, four int8 per word)
and return `int` C. Nothing in that path converts a quantized value to a floating
one, so weights loaded from a quantized file are used byte-for-byte as they sit in
memory:

- No stored quantized value is ever converted. All ten int8 kernels - the nine
  GEMM/GEMV shapes plus `rdna2_i8_row_sums` - contain **zero `v_cvt` of any
  kind** and issue only `v_dot4_*` (dp4a) or `v_wmma_i32_16x16x16_iu8` with
  integer accumulators; none contains `v_fmac_f32` or `v_mul_f32`. The single
  conversion in the family is one `v_cvt_f32_i32` in `rdna2_gemm_i8_epilogue`,
  which converts the finished int32 dot product once per output element.
- Semantics are signed int8. On gfx12 the signed form is expressed through the
  `neg_lo:[1,1,0]` modifier on `v_wmma_i32_16x16x16_iu8` (`v_dot4_i32_iu8` there is
  the unsigned form; see *dp4a on gfx12*).
- Exactness is tested, not assumed: every int8 variant (tiled 0-5) is compared
  against a 64-bit integer reference on inputs packed with `-128`/`+127`, with
  `|sum|` up to 1.33e8 - eight times past 2^24, where any f32 stage would have to
  round - and matches bit-for-bit. The host input buffers are also memcmp'd after
  the call to confirm the launcher never rewrites the caller's quants.
- No arch selector can route int8 onto a floating path: `rdna2_tiled_for_arch`
  returns 1, 3 or 5, all integer, and an unknown GPU gets 1.

The library applies the format's scale and zero point itself, so a loaded
quantized tensor goes straight into `rdna2_gemm_i8_scale_host`:
`out[m][n] = (acc[m][n] - zp_b * sum_a[m]) * sa[m] * sb[n]`. The correction
`-zp_b * sum_a[m]` is the exact integer form of B's asymmetric storage (a stored
value `q` represents `q - zp_b`), is computed in 64-bit on the int32 accumulator,
and skips the row-sum pass when `zp_b == 0`. This is not dequantization of the
stored data: the weights stay int8 and are never rewritten - the test memcmp's
the caller's buffers after the call - and the only float produced per output is
one finished dot product.

4-bit formats (Q4_0 / Q4_K and similar) are still not handled: dp4a consumes whole
int8 bytes, so they would need an in-kernel nibble unpack, again integer-only,
before this path could take them.

## OpenBLAS entry points

The int8 kernels are reachable from the library through two CBLAS-style
functions declared in `cblas.h` and implemented in `interface/gemm_i8.c`:

```c
int cblas_gemm_i8(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA,
                  enum CBLAS_TRANSPOSE TransB, blasint M, blasint N,
                  blasint K4, const int *A, blasint lda, const int *B,
                  blasint ldb, int *C, blasint ldc);

int cblas_gemm_i8_scale(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA,
                        enum CBLAS_TRANSPOSE TransB, blasint M, blasint N,
                        blasint K4, const int *A, blasint lda, const int *B,
                        blasint ldb, float *C, blasint ldc, const float *sa,
                        const float *sb, int zp_b);
```

Both compute `C := op(A) * op(B)` and both take the operands in the packed form
above, so a quantized weight buffer goes in byte for byte. `K4` is the number of
words per row (`K = 4*K4`); a partial trailing word is not accepted. `lda`,
`ldb` and `ldc` are leading dimensions in words (or C elements), with the usual
BLAS meaning per order: the stored row count in a column-major call and the row
stride in a row-major one. `CblasRowMajor` is supported and returns the same
numbers as the transposed column-major call.

The first form returns raw int32 accumulators and has no alpha/beta: the sum is
exact 32-bit two's complement, so it matches the naive integer reference bit for
bit and no tolerance applies. The second form is the quantized epilogue of
`rdna2_gemm_i8_scale_host`:

```
C[i][j] = (acc[i][j] - zp_b * rowsum_j) * sa[i] * sb[j]
```

where `rowsum_j` is the int8 sum of column `j` of `op(B)` (`rdna2_i8_row_sums` on
the device, skipped when `zp_b == 0`), the correction is taken in 64 bits, `sa`
has `M` entries for the rows of `op(A)`, `sb` has `N` for the columns of
`op(B)`, and either scale may be `NULL` for 1.0. Only the finished dot product
is converted to float, so no stored quant is rewritten.

Both return 0 on success, -1 on an invalid argument or layout and -2 if a host
workspace allocation failed, leaving `C` untouched. `M == 0` or `N == 0` is a
quick return, and `K4 == 0` zeroes `C`.

Dispatch: on a target built with `USE_HIP=1` (which defines `OPENBLAS_USE_HIP`
and an RDNA2 arch code) the call packs the operands into the kernels' row-major
layout and runs the tuned int8 shape for the device that is present. Everything
else - a CPU-only build, a machine with no device, or a launch that fails -
runs the portable C implementation in the same file, which computes the same
integer product and therefore the same bits. `OPENBLAS_RDNA2_GPU=0` forces the
CPU path on a GPU build. Neither path needs rocBLAS or a CUDA translation layer.

`test/test_gemm_i8.c` covers both entry points and runs from `make -C test
level3`: five shapes from 1x1 to 131x37 over four transpose combinations, both
orders, padded and unpadded leading dimensions, an accumulation past 2^24
(`K4 = 40000`), `K4 = 0`, argument validation, and a `512x512xK4=128` timing
pass with the GPU path forced on and off. On an RX 9070 XT (gfx1201,
`-DGFX1201`) all 185 checks pass with either path and no HIP device is needed
for the CPU build.

## LDS staging

The thread-per-cell int8 GEMM issues 8 packed global loads per k step per
thread; over a 64x64 tile that is 2048 words of traffic for 128 words of data.
The LDS kernels stage each `RDNA2_LDS_I8_KSTEP` (16 words) chunk cooperatively
and then reuse it from LDS:

- int8: two tiles of 64 x 17 words (stride 17 keeps rows off the same banks) =
  **8704 bytes** per block. Staging: 256 threads x 4 words (16-byte loads where
  `ncols % 4 == 0`), or 128 threads x 8 words for the 8x4 variant.
- f32: two tiles of 64 x 9 floats = **4608 bytes**. Staging: 256 threads x 2
  floats.
- f32, 16-deep K chunk (`rdna2_sgemm_f32_lds16`): two tiles of 64 x 17 floats =
  **8704 bytes**. 256 threads stage four floats each, one 16-byte vector load
  per thread per chunk.
- int8, 128x128 tile: two tiles of 128 x 17 words = **17408 bytes**. The same
  256 threads now stage 8 words each per chunk; the stride stays 17, so the
  bank-conflict argument is unchanged. `rdna2_dp4a_gemm_i8_wmma` reuses this
exact staging unchanged: its `RDNA2_WMMA_I8_KSTEP_WORDS` (4 words per lane per
16-byte K step) plus the `2 * (lane >> 4)` column-group offset is what makes one
LDS tile readable as the 16x16x16 A and B fragments, so the two 128-tile
kernels differ only in what consumes the staged data.
- f32, 128x128 tile: two tiles of 128 x 9 floats = **9216 bytes**.

Measured resource use on gfx1201 (`--cuda-device-only -S`, `-O2`) is 59 VGPRs /
8704 B for the 64x64 8x4 kernel, 93 VGPRs / 17408 B for the 128x128 dp4a
kernel and 96 VGPRs / 17408 B for the matrix-instruction kernel, with no
`scratch_` spills and no private segment in any of them. That register/LDS step
is what decides dp4a throughput: with 64 KiB of LDS per CU the 128x128 tile
admits 3 blocks where the 64x64 tile admits 7, and dp4a on this part is already
LDS-throughput bound, so the extra arithmetic intensity of the 8x8 blocking does
not pay for the lost residency (see *Performance*). The matrix instruction is
the exception: it does the same 3-block residency but needs 16x fewer issued
instructions for the same tile, which is why the same 17408 B of LDS is no
longer the limit there.

Rows and columns past the matrix edge, and words past `K4`, are staged as zero,
so the compute loop needs no bounds checks; `__syncthreads()` before the first
read and after the last read is what keeps the single-buffered tile legal.

## Launchers

All launchers return `0` on success, `1` on a HIP error, `2` on a bad argument,
and free every allocation on the error path.

- `rdna2_gemv_i8_host(y, A, x, rows, k4)`: `y[rows] = A[rows x k4] . x[k4]`, packed
  int8. Any `k4` (including `k4 % 4 != 0`) and any host buffer alignment are
  accepted; the kernel copies into an aligned `hipMalloc` staging buffer and the
  remainder is handled by the scalar tail.
- `rdna2_gemv_i8_host_variant(y, A, x, rows, k4, variant)`: forces `variant` 0
  (one thread per row) or 1 (one wave per row). `rdna2_gemv_i8_host` picks 1 from
  `k4 >= RDNA2_GEMV_WAVE_MIN_K4` (8) and 0 below it; `variant` 1 also falls back
  to the thread-per-row kernel if the wave size cannot be queried, because the
  wave kernel derives rows-per-block from `blockDim/warpSize` and a host guess
  would mis-size the grid.
- `rdna2_gemm_i8_host(C, A, Bt, M, N, K4, tiled)`: `C[M x N] = A[M x K4] * Bt[N x K4]^T`,
  packed int8. `tiled` selects 0 = naive, 1 = 4x4 tiled on global loads, 2 =
  4x4 with LDS staging, 3 = 8x4 with LDS staging, 4 = 128x128 tile with 8x8
  dp4a blocking, 5 = the same 128x128 tile driven by the 16x16x16 matrix
  instruction. 4 and 5 share a grid and a launch configuration; on a device
  without the instruction 5 runs the emulation described above (correct, slow),
  so the choice is safe but only worth making on gfx12. `M` and `N` need no
  multiple of the tile: edge threads read
  `0` and only store in-range elements, and the last K chunk is zero-padded.
- `rdna2_gemm_i8_scale_host(C, A, Bt, M, N, K4, tiled, sa, sb, zp_b)` and
  `rdna2_gemm_i8_scale_host_tuned(...)`: the same int8 GEMM followed by the
  format's scale and zero point, so a loaded quantized tensor goes in exactly as
  stored and the caller never expands it. `sa` scales rows of `A`, `sb` scales
  columns of `Bt` (either may be `NULL` for 1.0), and `zp_b` is B's zero point
  (`0` for symmetric weights, which skips the row-sum pass entirely). The result
  is f32 - `out[m][n] = (acc[m][n] - zp_b*sum_a[m]) * sa[m] * sb[n]` - because a
  GEMM result feeds the next layer's float path; no stored weight is converted.
  `tiled` selects the same shapes and grids as `rdna2_gemm_i8_host`.
- `rdna2_sgemv_f32_host(y, A, x, rows, K)` and
  `rdna2_sgemv_f32_host_variant(..., variant)`: f32 GEMV, `0` = one thread per
  row, `1` = one wave per row. `rdna2_sgemv_f32_host` picks 1 from
  `K >= RDNA2_SGEMV_WAVE_MIN_K` (256) and 0 below it; `variant` 1 also falls back
  to thread-per-row when the wave size cannot be queried, for the same grid-size
  reason as the int8 GEMV.
- `rdna2_sgemm_f32_host(C, A, Bt, M, N, K)`: f32 version of the tiled GEMM.
  `rdna2_sgemm_f32_host_v(..., variant)` selects 0 = global-load 4x4 tiling,
  1 = 64x64 LDS staging, 2 = 128x128 LDS staging, 3 = 64x64 LDS staging with a
  16-deep K chunk (the default).
- `rdna2_gemm_i8_host_tuned` / `rdna2_sgemm_f32_host_tuned`: same arguments,
  except that the shape is picked from the device that is present
  (`rdna2_tiled_for_arch` / `rdna2_sgemm_variant_for_arch`): gfx1030 and gfx1031
  get 3 / 3, gfx1200 and gfx1201 get 5 / 3, anything else gets the generic 1 / 1
  path.
  Measured on a 9070 XT, the 128x128 dp4a int8 tile never won, and the 128x128
  f32 tile only wins while its operand panels stay cache-resident, so `variant 2`
  is opt-in rather than an arch default. The 16-deep K chunk (variant 3) won
  everywhere it was tried, so it is the f32 default on all four arch codes.
  Without a device both entry points fall through to the CPU reference.
- `rdna2_dot_u8_host`, `rdna2_dot_f16_host`, `rdna2_prefix_host`.
- `rdna2_device_info(name, len, warp_size, arch_ok)`: reports the device name and
  wave size, plus a target code: `RDNA2_ARCH_GFX1031` (1), `RDNA2_ARCH_GFX1030`
  (2), `RDNA2_ARCH_GFX1201` (3), `RDNA2_ARCH_GFX1200` (4), or
  `RDNA2_ARCH_UNKNOWN` (0) for any other AMD GPU, which still runs.

When HIP reports no device, the GEMV/GEMM/SGEMM/dot launchers compute the same
result with `rdna2_gemm_ref.h`, so callers work on machines without a GPU.

## Build

    make test                                   # build and run the host test
    make HIP_ARCH=gfx1030                       # build the object for gfx1030
    make HIP_ARCH=gfx1201 test                  # RX 9070 / 9070 XT (RDNA4)
    make HIP_ARCH=gfx1200 test                  # RX 9060 XT (RDNA4)
    make multi                                  # fat binary: gfx1030 gfx1031 gfx1200 gfx1201

or directly:

    clang++ -x hip --offload-arch=gfx1031 -O2 -c rdna2_kernels.hip -o rdna2_kernels.o
    clang++ -x hip --offload-arch=gfx1201 --offload-arch=gfx1200 -O2 -c rdna2_kernels.hip -o rdna2_kernels.o

Use the AMD ROCm LLVM clang (or `amdclang++`) with the HIP headers and
`libamdhip64`. Pass `--offload-arch=` once per target to build a fat binary;
gfx1200 and gfx1201 are separate targets and a binary built for one will not be
picked up by the other, so a 9070 XT needs gfx1201 (or both) in the offload set.
`HIP_ARCH`/`ARCHS` changes do not invalidate existing objects, so `make clean`
between architectures or use `make multi`.

### Building the whole library with the GPU side

The knobs mirror the CUDA ones and live in `Makefile.rule`/`Makefile.system`:
`USE_HIP` turns the feature on, `GPU_ARCH` replaces `CUDA_ARCH`, `HIPCC` replaces
`NVCC`, and `HIP_PATH`/`HIP_LIB_PATH` replace `CUDADIR`.

    make USE_HIP=1 libs                          # CPU core from TARGET + gfx1201
    make USE_HIP=1 GPU_ARCH=gfx1031 libs         # same, for RDNA2
    make USE_HIP=1 GPU_ARCH="gfx1031 gfx1201" libs   # one fat library
    make USE_HIP=1 HIPCC=/opt/rocm/bin/hipcc libs

The CPU side is still whatever `TARGET` selects, so `TARGET=zen USE_HIP=1`
yields an ordinary Zen library in which large level-3 calls additionally run on
the device. `GPU_ARCH` defaults to the architecture of a `GFX*` TARGET, else
`gfx1201`. Switching `GPU_ARCH` rebuilds the device object (a stamp records the
arch list, since make cannot see flag changes); switching `USE_HIP` on or off
changes the compile flags of every object, so run `make clean` first.

The resulting library needs the HIP runtime at link time: `make` adds
`-L$(HIP_LIB_PATH) -lamdhip64 -lstdc++` for its own tests, and a program linking
the static archive must add the same.

On Windows with the ROCm HIP SDK three extra steps are environment, not code:

- `hipcc` cannot find the device libraries unless it is told where they are.
  clang looks under `lib/amdgcn`, while HIP 7.x installs them in
  `lib/llvm/amdgcn/bitcode`, so `rdna2_kernels.hip` fails with *cannot find
  ROCm device library*:

      make USE_HIP=1 GPU_ARCH=gfx1031 HIPFLAGS="-O2 -std=c++17 --rocm-path=C:/Rocm10/lib/llvm" libs

  `--rocm-device-lib-path` looks like the right flag for this and is not: clang
  23 accepts it and still fails to resolve the bitcode.

- The SDK ships only the MSVC import library (`amdhip64.lib`), so a mingw link
  needs one generated from the DLL:

      gendef C:/Rocm10/bin/amdhip64_7.dll
      dlltool -d amdhip64_7.def -l libamdhip64.a

- Which `amdhip64_*.dll` the process loads decides whether it sees a device at
  all. The SDK installs a versioned runtime next to `hipcc` and the display
  driver installs another one into `System32`, which Windows resolves first
  unless the executable sits in the SDK's `bin`. Against the driver's copy
  `hipGetDeviceCount` returns `hipErrorNoDevice` (100) even though `hipInfo.exe`
  lists the card, and every launcher silently takes its CPU reference path (see
  *CPU fallback is silent*). Run the test binaries from the SDK `bin`
  directory, or copy the SDK's `amdhip64_7.dll` beside them, before trusting any
  GPU number.

### Level-3 bridge test

    make level3-test                 # build the library, then run the test
    make -C gpu/rdna2 oblas-test     # run it against an existing library

The test checks every dispatched level-3 entry point -- SGEMM, SYMM, SYRK,
SYR2K, TRMM and TRSM -- against its own double-precision reference (an LU solve
for TRSM, so the reference does not share the bridge's triangle logic). Each op
runs over its full layout range: both store orders, side, uplo, trans and diag.
The contracts the bridges must keep are checked too: the unstored triangle and
the `ldb`/`ldc` padding are left untouched, and the unused triangle of a
symmetric or triangular factor plus the stored diagonal of a unit-diagonal one
are never read (they hold a value that no correct result can produce), which is
what fails loudly if a bridge reads the wrong half.

All six ops are seen at the default size gate as well, so the real dispatch
decision is covered end to end, not only the forced-open small shapes.

The scalars get their own pass on top of the layouts: `alpha` and `beta` at
zero, negative and large values (0, -0.5, -2, 1000). A zero scalar has to drop
its term exactly -- `alpha == 0` leaves `beta*C` and nothing else, `beta == 0`
drops the old `C`, both zero means exact zeros -- and the large values keep the
accumulation honest. TRMM and TRSM have no `beta`, so they sweep `alpha` only.

The TRSM factor is built row-wise diagonally dominant, each stored row summing
to a quarter of its diagonal. A triangular solve over a merely "small"
off-diagonal grows by orders of magnitude across a thousand rows, and the answer
then depends on the order the terms happened to be summed in rather than on the
triangle: the loosely conditioned 1300-row version reached `||A^-1|| ~ 2e8` and
flipped between builds of equivalent code. Dominance keeps the reference usable
at the gate size. The above-gate TRSM case is a dispatch check rather than an
accuracy one: at 1300 rows a float solve on this path is only good to about a
percent however it is compiled, so the layout and scalar cases carry the tight
tolerances.

It builds and passes on a CPU-only library, so it runs in ordinary CI (see `.github/workflows/rdna2_level3.yml`);
against a `USE_HIP` library it exercises the real dispatch, which
`OPENBLAS_RDNA2_GPU_DEBUG=1` makes visible as one `handled=1` line per call.
The link line adds the HIP runtime by itself when the archive contains the
device objects, so no ROCm is needed on machines that build without it.

### dp4a on gfx12

gfx12 dropped the `dot1-insts` feature that `__builtin_amdgcn_sdot4` needs, so on
gfx1200/gfx1201 the builtin does not compile: clang 23 fails with
`'__builtin_amdgcn_sdot4' needs target feature dot1-insts`. The hardware
instruction is still there under
the older mnemonic, so `rdna2_dot4_i8_impl` emits `v_dot4_i32_i8` directly for
those two targets and keeps the builtin for gfx10.3. Two things are worth
knowing, both checked on a real gfx1201 card:

- `v_dot4_i32_i8` on gfx12 computes **signed** results, byte for byte identical
  to the scalar reference on six adversarial vectors including `0xffffffff` and
  `0x00000080`.
- `v_dot4_i32_iu8`, the renamed dot10 form, assembles on gfx12 but computes
  **unsigned x unsigned**, so it must not be substituted for the signed mnemonic
  in a signed GEMM. It is only correct for the u8 dot kernel, which is what
  `rdna2_dp4a_dot_u8` uses.

`-DRDNA2_NO_DOT1=1` drops the asm and takes the other branch in
`rdna2_gemm_ref.h`: each int8 operand is sign-extended to int32 and LLVM forms
the signed dot itself (`v_dot4_i32_iu8 ... neg_lo:[1,1,0]` on gfx12, `v_dot4c_i32_i8`
on gfx10.3). Same results, more instructions; useful only if a toolchain rejects
the mnemonics. That build passed the whole suite as it stood when it was run,
before the matrix-instruction tests existed.

### Matrix instruction on gfx12 (`tiled == 5`)

gfx1200/gfx1201 have `v_wmma_i32_16x16x16_iu8`: a 16x16x16 integer matrix
multiply-accumulate per instruction, which for int8 is 4096 MACs where a dp4a is
4. `rdna2_dp4a_gemm_i8_wmma` is the 128x128 tile with that instruction, and it
is the int8 default for those two parts. Three things matter here; the first two
were measured on the card rather than taken from the ISA description.

- The operands are **signed**, so `RDNA2_WMMA_IU8` sets `negA` and `negB`. With
  those flags the dp4a `lds128` kernel and this one agree element for element on
  1024^3, 2048x2049x2048 and 4096^3 (16.7M outputs at the largest size).
- The lane-to-fragment mapping is fixed by the hardware and is not the one a
  natural reading of the ISA gives. On this part, lane `t` byte `e` of A carries
  `A[t % 16][(t / 16) * 8 + e]`, of B `B[(t / 16) * 8 + e][t % 16]`, and
  accumulator slot `l` of lane `t` carries `D[(t / 16) * 8 + l][t % 16]`. The
  host test checks exactly this: it runs one instruction's worth of
  fragments through `sim_wmma_i8_wave_step` and compares all 256 output slots
  against a per-slot dot product, plus a nonzero-accumulator pass to catch a
  kernel that overwrites instead of accumulating.
- The helper that performs one step, `rdna2_wmma_iu8_step`, has three branches.
  On gfx12 device passes (`RDNA2_WMMA_IU8`) it is the instruction. On any other
  *device* pass, i.e. gfx1030/gfx1031 and gfx11, there is no matrix unit, so the
  same result is built from explicit byte dots, and each slot's operand bytes are
  gathered with `__shfl` out of the two lanes that hold that row's fragment (the
  low half of K from lead lane `d`, the high half from `d + 16`). That branch is
  correct and is never selected by a dispatcher; no timing is claimed for it.
  In the plain C++ pass (the CPU test harness, and the host half of every HIP
  build) no kernel can run at all, so the body is a stub that returns `acc`
  unchanged: it exists only to keep the kernel defined so the unconditional
  launch site links.

That last point is the one that bit this tree once. `RDNA2_WMMA_IU8` is 0 in the
host pass, so a host-side `if (RDNA2_WMMA_IU8)` around the launch compiled out the
real kernel and quietly ran `_lds128` instead - a passing test that never touched
the matrix instruction. The launch is therefore unconditional and only the
*instruction* is arch-gated, which is also why `tiled == 5` on a non-gfx12 device
is safe but pointless.

## Performance

Measured on the RX 9070 XT in this machine (gfx1201, HIP 7.16, clang 23,
`-O2`). Times are kernel-only: allocation and transfers are hoisted out of the
loop, so they are not the launcher's end-to-end cost. Best of 3 runs x 10-20
iterations each.

int8, TFLOP/s equivalent (2 flops per MAC), `K = 1024` unless noted:

| M = N | naive | tiled 2 (4x4 LDS) | **tiled 3 (8x4 LDS)** | tiled 4 (128x128 dp4a) | **tiled 5 (128x128 matrix)** |
|---|---|---|---|---|---|
| 1024 | 0.74 | 33.6 | 38.2 | 30.1 | **91.3** |
| 2048 | 0.63 | 34.0 | 46.1 | 34.7 | **134.0** |
| 4096 | 0.59 | 34.7 | 46.4 | 41.6 | **125.0** |
| 2048, K = 4096 | 0.42 | 35.5 | **44.0** | 40.9 | - |

The three tiled-5 figures were timed in the same harness as tiled 4 (same grid,
launch and buffers, only the instruction differs): 29.4 / 40.6 / 40.1 against
91.3 / 134.0 / 125.0, i.e. 3.1x at 1024^3, 3.3x at 2048^2 and 3.1x at 4096^3
over the 128x128 dp4a kernel. Against the 64x64 8x4 dp4a kernel, which is the
previous int8 default, the gain is 2.4x, 2.9x and 2.7x.

f32, TFLOP/s, `K = 1024` unless noted:

| M = N | v0 (global 4x4) | v1 (64x64 LDS) | v2 (128x128 LDS) |
|---|---|---|---|
| 1024 | 2.94 | **8.11** | 7.51 |
| 2048 | - | 8.53 | **8.85** |
| 4096 | - | 8.01 | **9.47** |
| 6144 | - | 7.84 | **9.48** |
| 8192, K = 512 | - | 7.81 | **9.55** |
| 2048, K = 2048 | - | 6.39 | **8.14** |
| 4096, K = 2048 | - | 8.41 | **9.37** |
| 2048, K = 4096 | - | **7.34** | 3.86 |
| 2048, K = 8192 | - | **7.25** | 4.14 |
| 4096, K = 4096 | - | **7.60** | 6.62 |
| 1024, K = 4096 | - | 4.84 | **6.46** |

The v0/v1/v2 columns are the shapes this tree started with. Variant 3, the
64x64 tile with a 16-deep K chunk, was added later and measured against v1 at
1024^3, 1000^3, 2048^2 x K=2048/3072/4096/8192, 2304^2, 3584^2 and 4096^3: it
won at every one, by 3-42%, the margin shrinking as the problem grows. It is the
f32 default on all four arch codes for that reason. The v3 rows above come from
that separate run; the same-harness comparison of all four shapes is below.

All four f32 shapes in one run (same grid, launch and buffers, only the kernel
body differs; TMAC/s, with TFLOP/s in brackets):

| M = N = K | v0 (global 4x4) | v1 (64x64, KSTEP 8) | v2 (128x128) | **v3 (64x64, KSTEP 16)** |
|---|---|---|---|---|
| 1024 | 1.48 (2.96) | 4.06 (8.12) | 3.92 (7.84) | **4.29 (8.59)** |
| 2048 | 0.48 (0.96) | 3.05 (6.10) | 4.10 (8.19) | **4.23 (8.46)** |
| 4096 | 0.10 (0.20) | 4.29 (8.58) | 3.55 (7.10) | **4.46 (8.93)** |

What that says:

- Among the dp4a shapes, the 64x64 8x4 LDS kernel is the right int8 default on
  this part: it beats the 128x128 dp4a tile at every size measured, by 8% at
  2048^2 x K=4096, 11% at 4096^2 x K=1024 and 27% at 1024^3. The 128 tile still
  has the better arithmetic intensity on paper; in practice its 93 VGPRs and
  17408 B of LDS per block cost more residency than the extra reuse buys.
- The matrix instruction changes that conclusion without changing the tile: same
  128x128 tile, same 96 VGPR / 17408 B residency, 3.1-3.3x the throughput, so it
  is the int8 default on gfx12. It is also the one int8 number that does not sit
  under the dp4a ceiling below - the matrix unit is a separate pipe, so 134
  TFLOP/s equivalent is not 134 TFLOP/s of dp4a issue.
- The f32 128x128 tile is genuinely faster while the two operand panels fit in
  cache (up to +22% at 8192^2 x K=512, +21% at 6144^2 x K=1024, +18% at
  4096^2 x K=1024) and falls off a cliff once they do not: 2048^2 with K = 4096 or
  8192 drops to 3.9-4.1 TFLOP/s against 7.3 for the 64x64 shape, and 1024^3 also
  loses. That is why `variant 2` is opt-in rather than an arch default.
- For scale, the ceilings below are measured on this card by a probe kernel
  that does nothing but issue the instruction: dp4a 43.86 TMAC/s (87.7 int8
  TOP/s) and f32 fma 17.76 TMAC/s (35.5 TFLOP/s). `hipGetDeviceProperties`
  reports 32 CUs at 2.40 GHz, which would put the f32 issue rate at about 19.7
  TFLOPS and the dp4a rate at 39.3 TMAC/s; the probe beats both, so the reported
  CU count is half the physical one and the derived figures in this paragraph
  used to understate every ratio by 2x. Against the measured ceilings the best
  dp4a shape above (46.4 TFLOP/s equivalent, i.e. 23.2 TMAC/s) is 53% of the
  dp4a issue rate and the best f32 shape (9.55 TFLOP/s, i.e. 4.78 TMAC/s) is 27%
  of the 35.5 TFLOP/s f32 probe. The tiled-5 figures are not comparable to either ceiling: the
  matrix instruction replaces 4096 of the MACs a dp4a does 4 at a time.

### Where the time goes

A `trace`-style decomposition, all kernel-only and all in one harness, to check
whether the gaps above are memory, staging or inner-loop limits. Each GEMM row
below is the same kernel compiled three ways: `full`, `staging+sync` (the same
global loads and `__syncthreads`, math removed) and `math-only` (LDS reads and
fma kept, staging replaced by zeros). Dimensions 4096^3 unless noted.

| Kernel | full us | staging+sync us | math-only us | What it says |
|---|---|---|---|---|
| int8 tiled 5 (matrix), 4096^3 | 1219.0 | 1232.7 | 753.6 | staging >= full: the staged traffic alone costs as much as the whole kernel |
| int8 tiled 4 (128x128 dp4a) | 3436.2 | 1186.0 | 3164.0 | staging hidden, inner loop is 92% of full: this kernel is math-bound |
| f32 v3, 1024^3 | 255.1 | 83.7 | 233.8 | inner loop 92% of full |
| f32 v3, 2048^3 | 2075.9 | 620.2 | 1702.5 | inner loop 82% of full |
| f32 v3, 4096^3 | 15745.4 | 6903.4 | 13429.2 | inner loop 85% of full |

The f32 64x64 tile consumes 8.59 GB of panel reads at 4096^3, which at 15745 us
is 546 GB/s - close enough to the 558 GB/s copy ceiling to look DRAM bound. The
staging-only row disproves that: the same traffic alone takes 6903 us (1244 GB/s
effective, because consecutive tiles in a block row reuse the panels in L2), so
most of those bytes never reach DRAM and the kernel is not transfer limited.

The inner loops, isolated from staging entirely (fixed LDS tile, no global
access, 20000 iterations, 64/256/1024 blocks):

| Inner loop | 64 blocks | 256 blocks | 1024 blocks |
|---|---|---|---|
| f32 4x4 block, 16-deep K chunk (v3) | 4.77 TMAC/s | 4.93 | 4.66 |
| f32 4x4 block, unroll-4 variant | 4.81 TMAC/s | 4.91 | 4.64 |
| f32 8x8 block (the 128x128 tile, v2) | 5.71 TMAC/s | 5.79 | 5.58 |

Against the register-only fma probe of 17.76 TMAC/s that is 27% for the 4x4
block and 32% for the 8x8 block. Doubling the fma-per-LDS-read ratio from 2 to 4
buys 20%, and unrolling by 4 buys nothing, so this is not LDS bandwidth, not
LDS byte traffic and not instruction level parallelism that limits the f32 inner
loop. The same numbers explain the v2 anomaly from the other direction: its
inner loop is the fastest of the three (5.71) yet the whole kernel is the
slowest at 4096^3 (3.55 TMAC/s), so what the 128x128 tile loses is residency -
94 VGPR and 9216 B of LDS per block leave too few blocks in flight to cover the
staging it does not overlap.

Ceilings measured on the same card, for the numbers above to be read against:
copy 558.6 GB/s, HBM read 619.7 GB/s, LDS 8 streams 7644.1 GB/s, L2/Infinity
Cache read 7873.5 GB/s, dp4a issue 43.86 TMAC/s, f32 fma issue 17.76 TMAC/s.
The only int8 kernel near its ceiling is tiled 5, and it is near the *staging*
ceiling, not the ALU one; tiled 3 at 22.5 TMAC/s sits at 51% of the dp4a issue
rate, which the ISA explains (8 batched `ds_load_2addr_b32`, then counted
`v_dot4_i32_i8` with the LDS latency hidden behind them).

Tried and rejected on measurement, one variable at a time, all reverted:

- **Unroll the f32 K loop by 4.** 4.81 vs 4.77 TMAC/s isolated, 4.33 vs 4.46 at
  4096^3. No gain: the LDS latency was already covered, so more in-flight loads
  only add registers.
- **`hipHostRegister` the operand panels.** Pinning costs 484 us once per
  buffer, then copies 4 MB at 17 GB/s instead of 12.9 GB/s. End-to-end the
  registration cost is never repaid at any size the host wrappers take.
- **A 128x128 tile for f32 (variant 2) as the default.** Faster inner loop,
  slower kernel past L2-sized problems; see above.
- **A 64x128 tile with 8x4 per-thread blocking (32 accumulators, 61 VGPR,
  13056 B of LDS against v3's 38 / 8704).** Bit-identical to v3 at 1024^3 (0
  mismatches in 1048576 entries) and 3-5% faster at 2048^3 (1936 vs 1998 us), but
  ~2% slower at both 1024^3 (257 vs 250) and 4096^3 (15774 vs 15420) in two
  consecutive runs, so it cannot replace v3 as one default for every size. Its
  only advantage is issue density in the inner loop - 6 `ds_load_2addr_b32` per
  16 `v_dual_fmac` against v3's 4 per 8 - and the larger LDS footprint and
  61-VGPR register bill cancel it. A 17-shape interleaved A/B (128^2 to
  8192x2048, K from 128 to 8192) then settled the dispatch question: the tile
  never wins outside the measured run-to-run spread, and it loses 75-115% at
  128^2/256^2 - its 128-wide tile halves the block count, and the GPU already has
  too few blocks to fill - plus 6% at 512^2 and 21% at 768^2. A size-dependent
  dispatch is therefore not worth building; variant 3 stays the only f32 path.

The host wrappers are dominated by transfer and not by launch or allocation:
`hipMalloc`+`hipFree` of 1 MB is 1.6-2.1 us, a full 1024^3 int8 trio allocates in
5.2 us, but one 4 MB host-to-device copy is 311.7 us (12.9 GB/s) and a
4 MB H2D+D2H round trip 820 us. That floor is in the platform, not in this code,
and it is why `rdna2_gemv_i8_host` takes 558 us end to end for a 4096x1024 GEMV
whose kernel is 4.4 us, `rdna2_sgemv_f32_host` 1615 us for a 12.2 us kernel,
`rdna2_sgemm_f32_host` 1571 us at 1024^3 for a 251 us kernel, and
`rdna2_gemm_i8_host` 1035 us at 1024^3 tiled 5 for a 19.5 us kernel.

## Verification status

| Check | Status |
|---|---|
| Compiles for gfx1031, gfx1030, gfx1200 and gfx1201, `-Wall -Wextra` clean (0 warnings) | Run, passed (clang 23 / HIP 7.16) |
| Host test, 141 checks on a device / 132 on the CPU fallback path: GEMV int8 k4=8, k4=5, misaligned k4=5, wave variants; GEMV f32 K=128, 5 and 133 for both variants; GEMM 70x67xK4=10 naive/tiled/LDS/LDS8 plus grid simulations; LDS paths also at 35x33xK4=40 so a partial last chunk and edge tiles are staged; 135x130xK4=40 so both dimensions span two 128 tiles and the last chunk is 8 of 16, run for the 128x128 dp4a kernel and for the matrix-instruction kernel; a single-instruction fragment check (one step of `rdna2_wmma_i8_wave_step`, all 256 output slots against a per-slot dot product, plus a nonzero-accumulator pass); SGEMM 67x70xK=19 for both RDNA2 variants plus simulation, and 135x130xK=19 for variant 2; u8 dot n=1000; f16 dot n=500; DPP prefix; arch selectors; the `*_tuned` entry points; the quantized epilogue at 37x45xK4=12 (odd sizes, so partial 64 and 128 tiles) for tiled 0, 3, 4 and 5 with unit, row-only, column-only and both scales, plus `zp_b=8`, an exactness check on the zero-point correction, a memcmp of the caller's quants after the call, and the tuned entry point | Run, passed: 132 checks, 0 failures, `path=CPU fallback` (plain C++ incl. the file's host half) |
| Same suite on the **GPU** (`rdna2_kernels.hip` + host test linked with hipcc, device launches) | Run, passed on an RX 9070 XT (gfx1201, warp 32): `PASSED: 141 check(s), 0 failure(s), path=GPU` and `PASSED: 132 check(s), 0 failure(s), path=CPU fallback` with the device hidden (`HIP_VISIBLE_DEVICES=""`). The 9 extra device checks are the `rdna2_sgemv_f32_host` default-variant ones (three per K, times K = 128, 5 and 133): they only exist when a device is present, because on the fallback path the default launcher is the reference itself |
| Matrix-instruction GEMM (`tiled == 5`) against dp4a in ISA and in numbers | Run, passed: `rdna2_dp4a_gemm_i8_wmma` contains 8 `v_wmma_i32_16x16x16_iu8` and no dp4a; the `_lds128` kernel contains zero `v_wmma_*`. Cross-checked against `_lds128` at 1024^3, 2048x2049x2048 and 4096^3: **0 mismatches**, 16.7M outputs at the largest size |
| ISA of all 19 kernels for gfx1201 (`--cuda-device-only -S -O2`) | Run, passed: `ds_load_2addr_b32` / `ds_store_2addr_b32` in every LDS kernel, signed dp4a (`v_dot4_i32_i8`, 137 occurrences across the int8 kernels), `v_dot2_f32_f16` (1, from `__builtin_amdgcn_fdot2`), `v_wmma_i32_16x16x16_iu8` (8, all in the matrix kernel), `v_add_nc_u32_dpp` (4), **no `scratch_` spills and no private segment anywhere**. VGPR / LDS per kernel: gemv_i8 16/0, gemv_i8_wave 12/0, sgemv_f32 14/0, sgemv_f32_wave 11/0, gemm_i8_naive 7/0, gemm_i8_tile 43/0, gemm_i8_lds 43/8704, gemm_i8_lds8 59/8704, gemm_i8_lds128 93/17408, gemm_i8_wmma 96/17408, sgemm_f32 43/0, sgemm_f32_lds 35/4608, sgemm_f32_lds16 38/8704, sgemm_f32_lds128 94/9216, dot2_f16 8/0, dp4a_dot_u8 8/0, dpp_prefix 4/0, i8_row_sums 7/0, gemm_i8_epilogue 8/0 |
| No float conversion of stored quants | Run, passed: all ten int8 kernels contain **0 `v_cvt` of any kind** and **0 f32 arithmetic** (only `v_dot4_*` / `v_wmma_i32_16x16x16_iu8`); the only conversion in the family is one `v_cvt_f32_i32` in `rdna2_gemm_i8_epilogue`, on the finished accumulator. Adversarial inputs packed with `-128`/`+127` give `|sum|` up to 1.33e8 - eight times past 2^24, where any f32 stage would have to round - and all six `tiled` variants match a 64-bit integer reference bit-for-bit, with the caller's packed buffers byte-identical after the call |
| Signed dp4a semantics on gfx12 | Run, passed on the card: `v_dot4_i32_i8` matches the scalar reference on 6 vectors, `v_dot4_i32_iu8` does not (it is unsigned). gfx12 has no `sdot4` builtin; see *dp4a on gfx12*. The `-DRDNA2_NO_DOT1=1` build (no explicit asm) also passed the whole suite at the time it was run, which was before the matrix-instruction tests existed (78 checks then) |
| Device-vs-CPU correctness of every kernel and shape | Run, passed: the same test binary that the table above calls, all 141 checks on the device, including 135x130 tiles that straddle two 128-wide tiles |
| Performance of the five int8 shapes and four f32 shapes | Run, measured, see *Performance*; kernel-only, single card, one input distribution. The int8 shapes and the four f32 shapes were each timed in one harness, so both comparisons are like for like. Re-measured with staging / math phase splits and isolated inner-loop probes: see *Where the time goes* || Correctness of the timed shapes (1024^3 int8, tiled 3, 4 and 5, full 1024x1024 output vs `rdna2_gemm_i8_ref`, plus tiled 4 vs tiled 5 at 2048x2049x2048 and 4096^3) | Run, passed: exact match everywhere, 0 mismatches |

## CPU vs GPU (measured through the built library)

End-to-end (including allocation and PCIe transfers) against OpenBLAS shipped
with numpy 2.5.4, 32 threads on the host:

| size | GPU f32 GFLOPS | OpenBLAS CPU f32 GFLOPS | ratio |
|---|---|---|---|
| 4096³ | 3896 | 960 | **4.06x** |
| 2048³ | 2619 (standalone) / 81 (via lib, tiled transpose) | 1004 | **2.61x** (standalone) |
| 1024³ | 538 | 982 | 0.55x |

Correctness through the library (MATCH = max |diff| < 1e-3):
- `test_oblas_gpu.exe 512` → **MATCH** (0.000e+00 diff, 2.00 ms both paths, 134 GFLOPS)
- `test_oblas_gpu.exe 2048` → **MATCH** (1.717e‑04 diff, 216 ms GPU / 95 ms CPU, rms_rel 8.393e‑07)
- `probe_cpu.exe` → MATCH at 64, 200, 1024 (naive vs cblas_sgemm)
- int8 path exact (`max |diff| = 0`) at 256–2048, 12.49 TOPS end-to-end at 2048³;
  quant epilogue max_rel 1.19e-07.
- Library transpose is now tiled 64×64 (cache-resident); the naive copy was the
  ~200 ms strided-store bottleneck at k=2048.

Coverage limits: the timing tables in *Performance* are one RDNA4 card
(gfx1201) and are not re-timed on RDNA2. The gfx1030/gfx1031 path has since
been executed through the library on an RX 6700 XT - see *RDNA2 hardware* -
so the dispatch and the numerics are covered there, while the kernel-level
shapes remain RDNA4-only measurements.

The shuffle-based emulation of the matrix instruction is compiled (its
branch is visible in a fresh gfx1031 dump: 0 `v_wmma_*`, 167 `v_dot4c_i32_i8`,
83 `ds_bpermute_b32` of which 64 are the emulation's lane gathers, 0 warnings)
and never selected by a dispatcher, but it is not covered by a runtime device
test. Wave-level behaviour is covered on the device; the 141
checks are functional, not a stress or fuzz suite.

## RDNA2 hardware (measured on an RX 6700 XT, gfx1031)

The gfx1031 path has now been executed on a real RDNA2 card: RX 6700 XT (20 CUs
as `hipGetDeviceProperties` reports them, 12.0 GB), HIP 7.16, clang 23,
Windows 11, 12 host threads, `-O2`, with a library built
`USE_HIP=1 GPU_ARCH=gfx1031` as described above. gfx1030 is the same ISA family
and is still compiled and host-simulated only. Through that library, end to end
and including the transfers:

| size | CPU f32 GFLOPS (`OPENBLAS_RDNA2_GPU=0`) | GPU f32 GFLOPS | ratio |
|---|---|---|---|
| 1024^3 | 620.6 | 574.2 (below the 2e9 gate, so CPU) | 1.00x |
| 2048^3 | 589.8 | **693.3** | **1.18x** |
| 4096^3 | 654.5 | **1031.2** | **1.58x** |

On RDNA2 the bridge is a real but modest win - about 1.2x at 2048 and 1.6x at
4096 against a 12-thread host - and the default gate keeps everything below
roughly 1260^3 on the CPU, where the device would lose. That is well short of
the 2.6-4.1x the same dispatch gives on RDNA4 (above), which is what the kernel
shapes predict: gfx12 has the 16x16x16 matrix instruction for int8 and the
128x128 f32 tile, while gfx1031 has neither and runs the 64x64 dp4a (int8) and
the 64x64 KSTEP-16 (f32) shapes.

The level-3 bridge suite (`test_oblas_gpu.exe`) with the device visible:
**517 cases, 0 failures**, every op dispatch-confirmed, with the above-gate 1300
cases timed at ssymm 10 ms, ssyrk 11 ms, ssyr2k 20 ms, strmm 13 ms, trsm 110 ms,
and `sgemm 2048` 201 ms on the first call (25 ms on the same call warmed up; the
first one carries the code-object load).

int8 through the C API (`test_gemm_i8.c`): **185 checks, 0 failures**, both paths
exact, and 512x512xK4=128 in 3.0 ms on the device against 87.0 ms on the CPU
path - a 29x win.

### Where the 1.2x goes

Per shape, end to end through `rdna2_sgemm_f32_host_v` - allocation, upload and
download included, no bridge repacking - same card, `C = A * Bt^T`, ms / GFLOPS:

| M = N = K | v0 (global 4x4) | v1 (64x64, KSTEP 8) | v2 (128x128) | **v3 (64x64, KSTEP 16)** |
|---|---|---|---|---|
| 1024 | 21.0 / 102 | 2.42 / 886 | 2.46 / 874 | **1.85 / 1161** |
| 2048 | 107.1 / 160 | 8.73 / 1968 | 9.31 / 1846 | **8.28 / 2074** |
| 4096 | 779.0 / 176 | 53.2 / 2584 | 49.2 / 2796 | **47.0 / 2925** |

`rdna2_sgemm_f32_host_tuned` is v3 and reproduces that column. v3 is already the
arch default, and v0 - the shape the first version of this tree shipped - is
11-17x slower, so the kernel choice costs nothing here; the loss is elsewhere.

The same problem through `cblas_sgemm` costs 25.0 ms at 2048 and 133 ms at 4096,
against 8.28 ms and 47.0 ms for the call above. The difference, 16.7 ms and
86.0 ms, is `rdna2_gemm_cm_core`'s unconditional host work: it allocates an
`m x k` buffer and runs `rdna2_cm_to_rm` over A on every call. Nothing else is
in it - B passes through untouched when `ldb == k`, and with `alpha = 1`,
`beta = 0`, `ldc == m` the kernel writes directly into C, so there is no output
pass and no temporary. Those copies move 16.8 MB at 2048 and 67.1 MB at 4096 at
an effective 1.0 and 0.78 GB/s: the 64x64 tiling does not make a column-major
read plus a row-strided store pay for itself, and the destination is a fresh
allocation, first-touched, on every call.

| size | CPU (12 threads) | GPU call | bridge transpose | library total | ratio | without it |
|---|---|---|---|---|---|---|
| 2048^3 | 29.1 ms | 8.28 ms | 16.7 ms | 25.0 ms | 1.18x | 3.5x |
| 4096^3 | 210 ms | 47.0 ms | 86.0 ms | 133 ms | 1.58x | 4.5x |

So the kernel and the transfers are not what limits RDNA2: the device does the
4096^3 product in 47 ms end to end, 4.5x the 12-thread CPU, and 2.9 TFLOP/s. A
single per-call column-major to row-major pass over A gives most of that back.
Removing it - a per-thread scratch buffer, transposing once when the caller
reuses A, or teaching the tiled staging loop to read column-major A, which is a
stride swap - is worth about 3x at the gate sizes, far more than further kernel
tuning. The same pass is the most likely reason the RDNA4 library numbers above
sit below the kernel-only tables there.

### CPU fallback is silent

This matters before reading any number above as a GPU result. When
`rdna2_device_count() == 0`, the host layers call the file's own `rdna2_*_ref`
CPU implementation, return `RDNA2_OK`, and the bridge logs `handled=1` exactly
as it does after a real launch. Reached through the `System32` `amdhip64_7.dll`
above, the first run of the level-3 suite reported 517 passes and 518
`handled=1` lines **without touching the GPU**, at 3.1 GFLOPS - 200x slower than
the CPU - for 2048^3 and 4096^3 alike. The flat throughput and the absent size
dependence are the only symptoms. The same binary, run where the device is
visible, gives the 693 / 1031 GFLOPS and the 201 ms first-call timing in the
table above. Until the fallback returns a status distinct from a launch,
`handled=1` means "the call was answered", not "the device ran".

## Not included

- Full OpenBLAS (BLAS levels 1-3, LAPACK) on the GPU. OpenBLAS has no GPU backend;
  this directory provides GEMM/GEMV/dot kernels only.
