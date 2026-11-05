// Kernel arithmetic shared by the device kernels and the host test.
//
// These functions contain no HIP types and no wave intrinsics, so the exact
// same source compiles into the GPU kernels (rdna2_kernels.hip) and runs on the
// host in test_rdna2_host.cpp. That is how the tiling, bounds handling and
// accumulation order of the GEMM/GEMV kernels get executed on a machine with no
// GPU: the test reproduces the kernel's grid and calls the same body.
//
// What this does NOT cover: wave-level behaviour (shuffles, DPP, atomics) and
// anything the compiler does differently for the device. Those need hardware.
#ifndef RDNA2_KERNELS_MATH_H
#define RDNA2_KERNELS_MATH_H

#include <cstdint>
#include <cstddef>

#include "rdna2_gemm_ref.h"

// LDS tiling shape shared by the devices kernels and the host simulation: a
// block of 16x16 threads owns a 64x64 output tile and stages a K chunk of both
// operands in LDS. The row stride is padded by one element so that the lanes of
// a wave, which read rows 4 apart, do not all land on the same LDS banks.
#define RDNA2_TILE_MN 64
#define RDNA2_LDS_I8_KSTEP 16  // packed words staged per chunk
#define RDNA2_LDS_I8_STRIDE (RDNA2_LDS_I8_KSTEP + 1)
#define RDNA2_LDS_F32_KSTEP 8
#define RDNA2_LDS_F32_STRIDE (RDNA2_LDS_F32_KSTEP + 1)
// Wider f32 K chunk: 16 floats per tile row, so a chunk touches a full 64-byte
// line of each row instead of half of one, and the K loop pays half as many
// __syncthreads. Two 64x17 tiles are 8704 bytes, the same LDS as the int8 8x4
// kernel. Select with variant 3 (see README for the measurement).
#define RDNA2_LDS_F32B_KSTEP 16
#define RDNA2_LDS_F32B_STRIDE (RDNA2_LDS_F32B_KSTEP + 1)

// RDNA4 (gfx1200 / gfx1201, e.g. RX 9070 XT) shape: same 16x16 block of 256
// threads and the same K chunk, but a 128x128 output tile and 8x8 outputs per
// thread. gfx12 has 1536 VGPRs per SIMD against 1024 on gfx10/11, so 64 live
// accumulators plus 16 staged operands still leave room for several waves per
// SIMD. Sixteen LDS reads per k step now feed 64 dp4a instead of 8 feeding 16.
#define RDNA2_T4_TILE_MN 128
#define RDNA2_T4_LDS_I8_KSTEP 16
#define RDNA2_T4_LDS_I8_STRIDE (RDNA2_T4_LDS_I8_KSTEP + 1)
#define RDNA2_T4_LDS_F32_KSTEP 8
#define RDNA2_T4_LDS_F32_STRIDE (RDNA2_T4_LDS_F32_KSTEP + 1)
// Column groups of a 128-wide tile: one per 16x16x16 matrix instruction, so a
// wave carries eight accumulators in the tiled = 5 kernel.
#define RDNA2_T4_WMMA_NTILE (RDNA2_T4_TILE_MN / 16)

// RDNA4 matrix-instruction fragments (tiled = 5).
// V_WMMA_I32_16X16X16_IU8 spreads a 16x16x16 signed product over a wave32 like
// this (measured on gfx1201, 256/256 outputs matched a scalar reference):
//   A: lane t, byte e -> A[M = t%16][K = (t/16)*8 + e]
//   B: lane t, byte e -> B[K = (t/16)*8 + e][N = t%16]
//   D: lane t, slot l -> D[M = (t/16)*8 + l][N = t%16]
// The LDS tile holds packed words (four int8 along K per word), so a lane reads
// two consecutive words of row t%16, i.e. bytes 8*(t/16)..+7 of a 16-byte K
// step. These helpers keep the kernel and the host model on one set of indices;
// the instruction semantics themselves are the hardware's, not this file's.
#define RDNA2_WMMA_I8_KSTEP_WORDS 4  // staged packed words per 16-byte K step
// Staged row a lane takes its fragment from: A's M index and B's N index are
// both lane%16, and both operands stage rows of one tile, so the caller adds the
// wave's 16-row band or the column group it is working on.
RDNA2_HD static inline int rdna2_wmma_i8_row(int lane) { return lane & 15; }
// First of the two staged words a lane reads for K step `step` of a chunk.
RDNA2_HD static inline int rdna2_wmma_i8_frag_word(int lane, int step) {
  return RDNA2_WMMA_I8_KSTEP_WORDS * step + 2 * (lane >> 4);
}
// Output row (slot l) and column a lane's accumulators hold in a 16x16 result.
RDNA2_HD static inline int rdna2_wmma_i8_acc_row(int lane, int slot) {
  return 8 * (lane >> 4) + slot;
}
RDNA2_HD static inline int rdna2_wmma_i8_acc_col(int lane) { return lane & 15; }

// Target codes reported by rdna2_device_info(..., arch_ok) and consumed by the
// rdna2_*_for_arch selectors. 0 means "an AMD GPU from neither family"; the
// launchers still run on it with the generic 64x64 / 4x4 shapes.
#define RDNA2_ARCH_UNKNOWN 0
#define RDNA2_ARCH_GFX1031 1  // RDNA2, Navi 22 (RX 6700 XT / 6750 XT)
#define RDNA2_ARCH_GFX1030 2  // RDNA2, Navi 21 (RX 6800 / 6900)
#define RDNA2_ARCH_GFX1201 3  // RDNA4, Navi 48 (RX 9070 XT / 9070)
#define RDNA2_ARCH_GFX1200 4  // RDNA4, Navi 44 (RX 9060 XT)

// One row of the dp4a int8 GEMV: y[row] = sum over k4 of dot4(A_row[k], x[k]).
RDNA2_HD static inline int32_t rdna2_gemv_row_i8(const int32_t* a_row, const int32_t* x,
                                                int k4) {
  int32_t acc = 0;
  for (int j = 0; j < k4; ++j) acc = rdna2_dot4_i8_impl(a_row[j], x[j], acc);
  return acc;
}

// One output element of the int8 GEMM.
RDNA2_HD static inline int32_t rdna2_gemm_cell_i8(const int32_t* a_row, const int32_t* b_row,
                                                  int K4) {
  int32_t acc = 0;
  for (int k = 0; k < K4; ++k) acc = rdna2_dot4_i8_impl(a_row[k], b_row[k], acc);
  return acc;
}

// One thread's 4x4 output block of the tiled int8 GEMM, starting at (m0, n0).
// Rows/columns past M/N read 0 so the accumulation stays exact at the edges;
// the caller only stores the in-range elements.
RDNA2_HD static inline void rdna2_gemm_tile_i8(int32_t acc[4][4], const int32_t* A,
                                               const int32_t* Bt, int M, int N, int K4,
                                               int m0, int n0) {
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) acc[i][j] = 0;
  for (int k = 0; k < K4; ++k) {
    int32_t a[4];
    int32_t b[4];
    for (int i = 0; i < 4; ++i) a[i] = (m0 + i < M) ? A[(size_t)(m0 + i) * K4 + k] : 0;
    for (int j = 0; j < 4; ++j) b[j] = (n0 + j < N) ? Bt[(size_t)(n0 + j) * K4 + k] : 0;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) acc[i][j] = rdna2_dot4_i8_impl(a[i], b[j], acc[i][j]);
  }
}

// One thread's 4x4 output block of the f32 GEMM, same blocking as the int8 one.
RDNA2_HD static inline void rdna2_sgemm_tile_f32(float acc[4][4], const float* A,
                                                 const float* Bt, int M, int N, int K,
                                                 int m0, int n0) {
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) acc[i][j] = 0.0f;
  for (int k = 0; k < K; ++k) {
    float a[4];
    float b[4];
    for (int i = 0; i < 4; ++i) a[i] = (m0 + i < M) ? A[(size_t)(m0 + i) * K + k] : 0.0f;
    for (int j = 0; j < 4; ++j) b[j] = (n0 + j < N) ? Bt[(size_t)(n0 + j) * K + k] : 0.0f;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) acc[i][j] = __builtin_fmaf(a[i], b[j], acc[i][j]);
  }
}

// --- LDS-staged variants ---------------------------------------------------
// The tile pointer types are template parameters on purpose: on the device they
// deduce to __shared__-space pointers (so the compiler keeps the accesses as
// ds_read/ds_write instead of flat ones), on the host they are plain pointers.
// The tile is row-major with `ld` words per row holding `kwords` K columns;
// rows past the matrix edge are already zero, so no bounds checks are needed.
// acc is accumulated into, not zeroed, so a K loop can stream tile chunks in.
template <typename AT, typename BT>
RDNA2_HD static inline void rdna2_gemm_lds_block_i8(int32_t acc[4][4], AT a_tile, BT b_tile,
                                                    int ld, int lm, int ln, int kwords) {
  for (int kk = 0; kk < kwords; ++kk) {
    int32_t a[4];
    int32_t b[4];
    for (int i = 0; i < 4; ++i) a[i] = a_tile[(lm + i) * ld + kk];
    for (int j = 0; j < 4; ++j) b[j] = b_tile[(ln + j) * ld + kk];
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) acc[i][j] = rdna2_dot4_i8_impl(a[i], b[j], acc[i][j]);
  }
}

// 8x4 register blocking: 12 LDS reads feed 32 dp4a instead of 8 feeding 16, so
// the LDS traffic per multiply-accumulate drops by a quarter at the same math
// count. Register cost is 32 accumulators plus the 12 staged operands.
template <typename AT, typename BT>
RDNA2_HD static inline void rdna2_gemm_lds_block8_i8(int32_t acc[8][4], AT a_tile, BT b_tile,
                                                     int ld, int lm, int ln, int kwords) {
  for (int kk = 0; kk < kwords; ++kk) {
    int32_t a[8];
    int32_t b[4];
    for (int i = 0; i < 8; ++i) a[i] = a_tile[(lm + i) * ld + kk];
    for (int j = 0; j < 4; ++j) b[j] = b_tile[(ln + j) * ld + kk];
    for (int i = 0; i < 8; ++i)
      for (int j = 0; j < 4; ++j) acc[i][j] = rdna2_dot4_i8_impl(a[i], b[j], acc[i][j]);
  }
}

// 8x8 register blocking for the 128x128 tile. Sixteen LDS reads per k step feed
// 64 dp4a, i.e. 4 multiplies per LDS read against 2.67 for 8x4 and 2 for 4x4.
template <typename AT, typename BT>
RDNA2_HD static inline void rdna2_gemm_lds_block8x8_i8(int32_t acc[8][8], AT a_tile, BT b_tile,
                                                       int ld, int lm, int ln, int kwords) {
  for (int kk = 0; kk < kwords; ++kk) {
    int32_t a[8];
    int32_t b[8];
    for (int i = 0; i < 8; ++i) a[i] = a_tile[(lm + i) * ld + kk];
    for (int j = 0; j < 8; ++j) b[j] = b_tile[(ln + j) * ld + kk];
    for (int i = 0; i < 8; ++i)
      for (int j = 0; j < 8; ++j) acc[i][j] = rdna2_dot4_i8_impl(a[i], b[j], acc[i][j]);
  }
}

template <typename AT, typename BT>
RDNA2_HD static inline void rdna2_sgemm_lds_block8x8_f32(float acc[8][8], AT a_tile, BT b_tile,
                                                        int ld, int lm, int ln, int kwords) {
  for (int kk = 0; kk < kwords; ++kk) {
    float a[8];
    float b[8];
    for (int i = 0; i < 8; ++i) a[i] = a_tile[(lm + i) * ld + kk];
    for (int j = 0; j < 8; ++j) b[j] = b_tile[(ln + j) * ld + kk];
    for (int i = 0; i < 8; ++i)
      for (int j = 0; j < 8; ++j) acc[i][j] = __builtin_fmaf(a[i], b[j], acc[i][j]);
  }
}

template <typename AT, typename BT>
RDNA2_HD static inline void rdna2_sgemm_lds_block_f32(float acc[4][4], AT a_tile, BT b_tile,
                                                      int ld, int lm, int ln, int kwords) {
  for (int kk = 0; kk < kwords; ++kk) {
    float a[4];
    float b[4];
    for (int i = 0; i < 4; ++i) a[i] = a_tile[(lm + i) * ld + kk];
    for (int j = 0; j < 4; ++j) b[j] = b_tile[(ln + j) * ld + kk];
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) acc[i][j] = __builtin_fmaf(a[i], b[j], acc[i][j]);
  }
}

// One lane's contribution to a row when a whole wave cooperates on that row
// (lanes stride the row, so consecutive lanes read consecutive words). The
// caller reduces the result across the wave.
RDNA2_HD static inline int32_t rdna2_gemv_lane_partial_i8(const int32_t* a_row, const int32_t* x,
                                                          int k4, int lane, int wave_size) {
  int32_t acc = 0;
  for (int k = lane; k < k4; k += wave_size) acc = rdna2_dot4_i8_impl(a_row[k], x[k], acc);
  return acc;
}

// One row of the f32 GEMV: y[row] = sum over k of A_row[k] * x[k]. Four
// independent accumulators, because a single fma chain would leave the fma
// pipes waiting on the previous result; the final sum fixes the evaluation
// order, so the kernel and the host simulation agree exactly.
RDNA2_HD static inline float rdna2_sgemv_row_f32(const float* a_row, const float* x, int K) {
  float acc0 = 0.0f, acc1 = 0.0f, acc2 = 0.0f, acc3 = 0.0f;
  int k = 0;
  for (; k + 4 <= K; k += 4) {
    acc0 = __builtin_fmaf(a_row[k + 0], x[k + 0], acc0);
    acc1 = __builtin_fmaf(a_row[k + 1], x[k + 1], acc1);
    acc2 = __builtin_fmaf(a_row[k + 2], x[k + 2], acc2);
    acc3 = __builtin_fmaf(a_row[k + 3], x[k + 3], acc3);
  }
  for (; k < K; ++k) acc0 = __builtin_fmaf(a_row[k], x[k], acc0);
  return (acc0 + acc1) + (acc2 + acc3);
}

// One lane's contribution to a row when a wave cooperates on that row.
RDNA2_HD static inline float rdna2_sgemv_lane_partial_f32(const float* a_row, const float* x,
                                                          int K, int lane, int wave_size) {
  float acc = 0.0f;
  for (int k = lane; k < K; k += wave_size) acc = __builtin_fmaf(a_row[k], x[k], acc);
  return acc;
}

#endif  // RDNA2_KERNELS_MATH_H
