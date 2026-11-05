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
| Performance of the five int8 shapes and four f32 shapes | Run, measured, see *Performance*; kernel-only, single card, one input distribution. The int8 shapes and the four f32 shapes were each timed in one harness, so both comparisons are like for like. Re-measured with staging / math phase splits and isolated inner-loop probes: see *Where the time goes* |
| Correctness of the timed shapes (1024^3 int8, tiled 3, 4 and 5, full 1024x1024 output vs `rdna2_gemm_i8_ref`, plus tiled 4 vs tiled 5 at 2048x2049x2048 and 4096^3) | Run, passed: exact match everywhere, 0 mismatches |

Coverage limits: everything above is one RDNA4 card (gfx1201). The gfx1030 and
gfx1031 paths are compiled and host-simulated but were not executed on RDNA2
hardware, so the dp4a-in-LDS timing on RDNA2 is the same code measured only on
RDNA4. The shuffle-based emulation of the matrix instruction is compiled (its
branch is visible in a fresh gfx1031 dump: 0 `v_wmma_*`, 167 `v_dot4c_i32_i8`,
83 `ds_bpermute_b32` of which 64 are the emulation's lane gathers, 0 warnings)
and never selected by a dispatcher, but it is not covered by a runtime device
test. Wave-level behaviour is covered on the device; the 141
checks are functional, not a stress or fuzz suite.

## Not included

- Full OpenBLAS (BLAS levels 1-3, LAPACK) on the GPU. OpenBLAS has no GPU backend;
  this directory provides GEMM/GEMV/dot kernels only.
