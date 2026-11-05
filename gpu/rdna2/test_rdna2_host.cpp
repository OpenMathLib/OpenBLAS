// Correctness test for the RDNA2 kernel launchers.
//
// On a machine with a RDNA2 GPU, every check runs on the device (the launchers
// pick the device path automatically) and is compared against the plain C++
// reference. With no GPU the same checks exercise the CPU fallback, so the
// binary is useful either way; the header line says which path ran.
//
// Build and run: see Makefile (make test).
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>

#include "rdna2_gemm_ref.h"
#include "rdna2_kernels_math.h"

// Launchers from rdna2_kernels.hip.
int rdna2_gemv_i8_host(int* y, const int* A, const int* x, int rows, int k4);
int rdna2_gemv_i8_host_variant(int* y, const int* A, const int* x, int rows, int k4, int variant);
int rdna2_gemm_i8_host(int* C, const int* A, const int* Bt, int M, int N, int K4, int tiled);
int rdna2_gemm_i8_scale_host(float* C, const int* A, const int* Bt, int M, int N, int K4,
                             int tiled, const float* sa, const float* sb, int zp_b);
int rdna2_gemm_i8_scale_host_tuned(float* C, const int* A, const int* Bt, int M, int N, int K4,
                                   const float* sa, const float* sb, int zp_b);
int rdna2_sgemm_f32_host(float* C, const float* A, const float* Bt, int M, int N, int K);
int rdna2_sgemm_f32_host_v(float* C, const float* A, const float* Bt, int M, int N, int K,
                           int variant);
int rdna2_sgemv_f32_host(float* y, const float* A, const float* x, int rows, int K);
int rdna2_sgemv_f32_host_variant(float* y, const float* A, const float* x, int rows, int K,
                                 int variant);
int rdna2_dot_u8_host(unsigned* out, const unsigned* a, const unsigned* b, int n);
int rdna2_dot_f16_host(float* out, const unsigned* a, const unsigned* b, int n);
int rdna2_prefix_host(int* out, const int* in, int n);
int rdna2_device_info(char* name, size_t len, int* warp_size, int* arch_ok);
int rdna2_gemm_i8_host_tuned(int* C, const int* A, const int* Bt, int M, int N, int K4);
int rdna2_sgemm_f32_host_tuned(float* C, const float* A, const float* Bt, int M, int N, int K);
int rdna2_tiled_for_arch(int arch);
int rdna2_sgemm_variant_for_arch(int arch);

static int failures = 0;
static int checks = 0;

// Written nowhere if a kernel misses an output element.
static const int32_t SENTINEL = 0x7eadbeef;

#define CHECK(cond, msg)                                                     \
  do {                                                                       \
    ++checks;                                                                \
    if (!(cond)) {                                                           \
      std::printf("FAIL: %s (line %d)\n", msg, __LINE__);                    \
      ++failures;                                                            \
    }                                                                        \
  } while (0)

// Deterministic pseudo-random int32 (full 32-bit range, so every byte pattern
// including negative int8 values appears).
static uint32_t rng_state = 12345u;
static uint32_t rng(void) {
  rng_state ^= rng_state << 13;
  rng_state ^= rng_state >> 17;
  rng_state ^= rng_state << 5;
  return rng_state;
}

// 16-byte aligned view of a vector, optionally offset by `misalign` bytes.
static int32_t* align16(std::vector<int32_t>& store, size_t n, size_t misalign = 0) {
  store.assign(n + 8, 0);
  uintptr_t p = reinterpret_cast<uintptr_t>(store.data());
  uintptr_t base = (p + 15) & ~static_cast<uintptr_t>(15);
  return reinterpret_cast<int32_t*>(base + misalign);
}

// ---------------------------------------------------------------------------
// Host simulation of the device kernels.
//
// These reproduce the grid, block size and per-thread indexing of the kernels
// in rdna2_kernels.hip and call the same per-thread body from
// rdna2_kernels_math.h. On a machine with no GPU this is the only way to run
// the tiling/bounds logic at all; it does not model waves, shuffles or atomics.
// Out-of-range outputs are left untouched, so the caller can seed the result
// buffer with a sentinel to detect stray writes.
// ---------------------------------------------------------------------------
static void sim_gemm_tile_i8(int32_t* C, const int32_t* A, const int32_t* Bt, int M, int N,
                             int K4) {
  const int gx = (N + 63) / 64, gy = (M + 63) / 64;  // block(16,16), 4x4 per thread
  for (int by = 0; by < gy; ++by)
    for (int bx = 0; bx < gx; ++bx)
      for (int ty = 0; ty < 16; ++ty)
        for (int tx = 0; tx < 16; ++tx) {
          const int m0 = (by * 16 + ty) * 4;
          const int n0 = (bx * 16 + tx) * 4;
          if (m0 >= M || n0 >= N) continue;
          int32_t acc[4][4];
          rdna2_gemm_tile_i8(acc, A, Bt, M, N, K4, m0, n0);
          for (int i = 0; i < 4 && m0 + i < M; ++i)
            for (int j = 0; j < 4 && n0 + j < N; ++j)
              C[(size_t)(m0 + i) * N + (n0 + j)] = acc[i][j];
        }
}

static void sim_sgemm_tile_f32(float* C, const float* A, const float* Bt, int M, int N, int K) {
  const int gx = (N + 63) / 64, gy = (M + 63) / 64;
  for (int by = 0; by < gy; ++by)
    for (int bx = 0; bx < gx; ++bx)
      for (int ty = 0; ty < 16; ++ty)
        for (int tx = 0; tx < 16; ++tx) {
          const int m0 = (by * 16 + ty) * 4;
          const int n0 = (bx * 16 + tx) * 4;
          if (m0 >= M || n0 >= N) continue;
          float acc[4][4];
          rdna2_sgemm_tile_f32(acc, A, Bt, M, N, K, m0, n0);
          for (int i = 0; i < 4 && m0 + i < M; ++i)
            for (int j = 0; j < 4 && n0 + j < N; ++j)
              C[(size_t)(m0 + i) * N + (n0 + j)] = acc[i][j];
        }
}

static void sim_gemv_i8(int32_t* y, const int32_t* A, const int32_t* x, int rows, int k4) {
  for (int r = 0; r < rows; ++r) y[r] = rdna2_gemv_row_i8(A + (size_t)r * k4, x, k4);
}

// Wave-per-row GEMV: 128 threads per block, so blockDim/warp_size independent
// rows per block, lanes striding the row. The shuffle reduction is modelled as
// the exact sum of the lane partials (integer adds, so order cannot matter).
static void sim_gemv_wave_i8(int32_t* y, const int32_t* A, const int32_t* x, int rows, int k4,
                             int warp_size) {
  const int threads = 128;
  const int waves = threads / warp_size;
  const int blocks = (rows + waves - 1) / waves;
  for (int b = 0; b < blocks; ++b)
    for (int w = 0; w < waves; ++w) {
      const int row = b * waves + w;
      if (row >= rows) continue;
      int64_t sum = 0;
      for (int lane = 0; lane < warp_size; ++lane)
        sum += rdna2_gemv_lane_partial_i8(A + (size_t)row * k4, x, k4, lane, warp_size);
      y[row] = (int32_t)sum;
    }
}

static void sim_sgemv_f32(float* y, const float* A, const float* x, int rows, int K) {
  const int threads = 64;
  const int blocks = (rows + threads - 1) / threads;
  for (int b = 0; b < blocks; ++b)
    for (int t = 0; t < threads; ++t) {
      const int row = b * threads + t;
      if (row >= rows) continue;
      y[row] = rdna2_sgemv_row_f32(A + (size_t)row * K, x, K);
    }
}

// The device reduces with __shfl_xor, so its grouping differs from this serial
// sum of the lane partials; the test compares with a tolerance, not bit-for-bit.
static void sim_sgemv_wave_f32(float* y, const float* A, const float* x, int rows, int K,
                               int warp_size) {
  const int threads = 128;
  const int waves = threads / warp_size;
  const int blocks = (rows + waves - 1) / waves;
  for (int b = 0; b < blocks; ++b)
    for (int w = 0; w < waves; ++w) {
      const int row = b * waves + w;
      if (row >= rows) continue;
      float sum = 0.0f;
      for (int lane = 0; lane < warp_size; ++lane)
        sum += rdna2_sgemv_lane_partial_f32(A + (size_t)row * K, x, K, lane, warp_size);
      y[row] = sum;
    }
}

// Staging helper of rdna2_stage_i8 / rdna2_stage_f32: one thread copies 4 (i8)
// or 2 (f32) consecutive words of one tile row, writing zeros past the edges.
static void sim_stage_i8(int32_t* lds, const int32_t* src, int nrows, int ncols, int row0, int k0,
                         int ld_row, int ld_col) {
  const int gm = row0 + ld_row;
  int32_t v[4] = {0, 0, 0, 0};
  if (gm < nrows) {
    const int32_t* p = src + (size_t)gm * ncols + k0 + ld_col;
    for (int c = 0; c < 4; ++c)
      if (k0 + ld_col + c < ncols) v[c] = p[c];
  }
  for (int c = 0; c < 4; ++c) lds[(size_t)ld_row * RDNA2_LDS_I8_STRIDE + ld_col + c] = v[c];
}

static void sim_stage_f32(float* lds, const float* src, int nrows, int ncols, int row0, int k0,
                          int ld_row, int ld_col) {
  const int gm = row0 + ld_row;
  float v[2] = {0.0f, 0.0f};
  if (gm < nrows) {
    const float* p = src + (size_t)gm * ncols + k0 + ld_col;
    for (int c = 0; c < 2; ++c)
      if (k0 + ld_col + c < ncols) v[c] = p[c];
  }
  for (int c = 0; c < 2; ++c) lds[(size_t)ld_row * RDNA2_LDS_F32_STRIDE + ld_col + c] = v[c];
}

static void sim_stage_f32_4(float* lds, const float* src, int nrows, int ncols, int row0, int k0,
                            int ld_row, int ld_col) {
  const int gm = row0 + ld_row;
  float v[4] = {0.0f, 0.0f, 0.0f, 0.0f};
  if (gm < nrows) {
    const float* p = src + (size_t)gm * ncols + k0 + ld_col;
    for (int c = 0; c < 4; ++c)
      if (k0 + ld_col + c < ncols) v[c] = p[c];
  }
  for (int c = 0; c < 4; ++c) lds[(size_t)ld_row * RDNA2_LDS_F32B_STRIDE + ld_col + c] = v[c];
}

// LDS-staged int8 GEMM: same output tile and register blocking as the plain
// tiled kernel, but each K chunk is staged for the whole block first (modelling
// the two __syncthreads) and then read back per thread.
static void sim_gemm_lds_i8(int32_t* C, const int32_t* A, const int32_t* Bt, int M, int N,
                            int K4) {
  std::vector<int32_t> at((size_t)RDNA2_TILE_MN * RDNA2_LDS_I8_STRIDE);
  std::vector<int32_t> bt(at.size());
  std::vector<int32_t> accs((size_t)256 * 16);
  const int gx = (N + RDNA2_TILE_MN - 1) / RDNA2_TILE_MN;
  const int gy = (M + RDNA2_TILE_MN - 1) / RDNA2_TILE_MN;
  for (int by = 0; by < gy; ++by)
    for (int bx = 0; bx < gx; ++bx) {
      const int m_blk = by * RDNA2_TILE_MN, n_blk = bx * RDNA2_TILE_MN;
      for (auto& v : accs) v = 0;
      for (int k0 = 0; k0 < K4; k0 += RDNA2_LDS_I8_KSTEP) {
        for (int tid = 0; tid < 256; ++tid) {
          sim_stage_i8(at.data(), A, M, K4, m_blk, k0, tid >> 2, (tid & 3) << 2);
          sim_stage_i8(bt.data(), Bt, N, K4, n_blk, k0, tid >> 2, (tid & 3) << 2);
        }
        int kend = K4 - k0;
        if (kend > RDNA2_LDS_I8_KSTEP) kend = RDNA2_LDS_I8_KSTEP;
        for (int tid = 0; tid < 256; ++tid) {
          int32_t(*acc)[4] = reinterpret_cast<int32_t(*)[4]>(&accs[(size_t)tid * 16]);
          rdna2_gemm_lds_block_i8(acc, at.data(), bt.data(), RDNA2_LDS_I8_STRIDE, (tid >> 4) * 4,
                                  (tid & 15) * 4, kend);
        }
      }
      for (int tid = 0; tid < 256; ++tid) {
        const int m0 = m_blk + (tid >> 4) * 4, n0 = n_blk + (tid & 15) * 4;
        const int32_t(*acc)[4] = reinterpret_cast<const int32_t(*)[4]>(&accs[(size_t)tid * 16]);
        for (int i = 0; i < 4 && m0 + i < M; ++i)
          for (int j = 0; j < 4 && n0 + j < N; ++j) C[(size_t)(m0 + i) * N + (n0 + j)] = acc[i][j];
      }
    }
}

// Same LDS staging, 8x4 register blocking: 128 threads per 64x64 tile and two
// staged words per thread per half-row group.
static void sim_gemm_lds8_i8(int32_t* C, const int32_t* A, const int32_t* Bt, int M, int N,
                             int K4) {
  std::vector<int32_t> at((size_t)RDNA2_TILE_MN * RDNA2_LDS_I8_STRIDE);
  std::vector<int32_t> bt(at.size());
  std::vector<int32_t> accs((size_t)128 * 32);
  const int gx = (N + RDNA2_TILE_MN - 1) / RDNA2_TILE_MN;
  const int gy = (M + RDNA2_TILE_MN - 1) / RDNA2_TILE_MN;
  for (int by = 0; by < gy; ++by)
    for (int bx = 0; bx < gx; ++bx) {
      const int m_blk = by * RDNA2_TILE_MN, n_blk = bx * RDNA2_TILE_MN;
      for (auto& v : accs) v = 0;
      for (int k0 = 0; k0 < K4; k0 += RDNA2_LDS_I8_KSTEP) {
        for (int tid = 0; tid < 128; ++tid) {
          sim_stage_i8(at.data(), A, M, K4, m_blk, k0, tid >> 1, (tid & 1) << 3);
          sim_stage_i8(at.data(), A, M, K4, m_blk, k0, tid >> 1, ((tid & 1) << 3) + 4);
          sim_stage_i8(bt.data(), Bt, N, K4, n_blk, k0, tid >> 1, (tid & 1) << 3);
          sim_stage_i8(bt.data(), Bt, N, K4, n_blk, k0, tid >> 1, ((tid & 1) << 3) + 4);
        }
        int kend = K4 - k0;
        if (kend > RDNA2_LDS_I8_KSTEP) kend = RDNA2_LDS_I8_KSTEP;
        for (int tid = 0; tid < 128; ++tid) {
          int32_t(*acc)[4] = reinterpret_cast<int32_t(*)[4]>(&accs[(size_t)tid * 32]);
          rdna2_gemm_lds_block8_i8(acc, at.data(), bt.data(), RDNA2_LDS_I8_STRIDE, (tid >> 4) * 8,
                                   (tid & 15) * 4, kend);
        }
      }
      for (int tid = 0; tid < 128; ++tid) {
        const int m0 = m_blk + (tid >> 4) * 8, n0 = n_blk + (tid & 15) * 4;
        const int32_t(*acc)[4] = reinterpret_cast<const int32_t(*)[4]>(&accs[(size_t)tid * 32]);
        for (int i = 0; i < 8 && m0 + i < M; ++i)
          for (int j = 0; j < 4 && n0 + j < N; ++j) C[(size_t)(m0 + i) * N + (n0 + j)] = acc[i][j];
      }
    }
}

// RDNA4 shape: 128x128 tile, 256 threads, 16 staged words per row per K step,
// 8x8 accumulators each. The staging helper is the same one the 64x64 kernels
// use (the padded stride is 17 words on both shapes), only the row count grows.
static void sim_gemm_lds8x8_i8(int32_t* C, const int32_t* A, const int32_t* Bt, int M, int N,
                               int K4) {
  std::vector<int32_t> at((size_t)RDNA2_T4_TILE_MN * RDNA2_T4_LDS_I8_STRIDE);
  std::vector<int32_t> bt(at.size());
  std::vector<int32_t> accs((size_t)256 * 64);
  const int gx = (N + RDNA2_T4_TILE_MN - 1) / RDNA2_T4_TILE_MN;
  const int gy = (M + RDNA2_T4_TILE_MN - 1) / RDNA2_T4_TILE_MN;
  for (int by = 0; by < gy; ++by)
    for (int bx = 0; bx < gx; ++bx) {
      const int m_blk = by * RDNA2_T4_TILE_MN, n_blk = bx * RDNA2_T4_TILE_MN;
      for (auto& v : accs) v = 0;
      for (int k0 = 0; k0 < K4; k0 += RDNA2_T4_LDS_I8_KSTEP) {
        for (int tid = 0; tid < 256; ++tid) {
          sim_stage_i8(at.data(), A, M, K4, m_blk, k0, tid >> 1, (tid & 1) << 3);
          sim_stage_i8(at.data(), A, M, K4, m_blk, k0, tid >> 1, ((tid & 1) << 3) + 4);
          sim_stage_i8(bt.data(), Bt, N, K4, n_blk, k0, tid >> 1, (tid & 1) << 3);
          sim_stage_i8(bt.data(), Bt, N, K4, n_blk, k0, tid >> 1, ((tid & 1) << 3) + 4);
        }
        int kend = K4 - k0;
        if (kend > RDNA2_T4_LDS_I8_KSTEP) kend = RDNA2_T4_LDS_I8_KSTEP;
        for (int tid = 0; tid < 256; ++tid) {
          int32_t(*acc)[8] = reinterpret_cast<int32_t(*)[8]>(&accs[(size_t)tid * 64]);
          rdna2_gemm_lds_block8x8_i8(acc, at.data(), bt.data(), RDNA2_T4_LDS_I8_STRIDE,
                                     (tid >> 4) * 8, (tid & 15) * 8, kend);
        }
      }
      for (int tid = 0; tid < 256; ++tid) {
        const int m0 = m_blk + (tid >> 4) * 8, n0 = n_blk + (tid & 15) * 8;
        const int32_t(*acc)[8] = reinterpret_cast<const int32_t(*)[8]>(&accs[(size_t)tid * 64]);
        for (int i = 0; i < 8 && m0 + i < M; ++i)
          for (int j = 0; j < 8 && n0 + j < N; ++j) C[(size_t)(m0 + i) * N + (n0 + j)] = acc[i][j];
      }
    }
}

// Wave-level model of V_WMMA_I32_16X16X16_IU8 (signed x signed, no clamp): one
// call is one wave executing one matrix instruction, acc is [lane][slot]. What
// is modeled is the operand/result placement contract the kernel is written
// against, stated once in rdna2_kernels_math.h and measured on a gfx1201 card
// (256/256 outputs of a signed product). The multiplication itself is trivially
// a 16-deep int8 dot; the hardware is trusted for that, and the GPU run of the
// tiled = 5 test is what confirms this model against silicon.
static void sim_wmma_i8_wave_step(int32_t acc[32][8], const int32_t a_frag[32][2],
                                  const int32_t b_frag[32][2]) {
  for (int t = 0; t < 32; ++t) {
    const int dn = rdna2_wmma_i8_acc_col(t);
    for (int l = 0; l < 8; ++l) {
      // D row this lane's slot l holds. Byte k of that row sits in lane
      // row + 16*(k/8), at byte k%8 of that lane's eight-byte fragment.
      const int dm = rdna2_wmma_i8_acc_row(t, l);
      int32_t d = 0;
      for (int k = 0; k < 16; ++k) {
        const int8_t av = (int8_t)((a_frag[dm + (k >> 3) * 16][(k & 7) >> 2] >> (8 * (k & 3))) & 0xff);
        const int8_t bv = (int8_t)((b_frag[dn + (k >> 3) * 16][(k & 7) >> 2] >> (8 * (k & 3))) & 0xff);
        d += (int32_t)av * (int32_t)bv;
      }
      acc[t][l] += d;
    }
  }
}

// The tiled = 5 kernel on the host: 128x128 tile, 256 threads, eight waves each
// owning a 16-row band and eight 16-column groups, staged exactly as the lds128
// kernel stages. Vector arithmetic is unrolled to whole-wave arrays here.
static void sim_gemm_wmma_i8(int32_t* C, const int32_t* A, const int32_t* Bt, int M, int N, int K4) {
  const int S = RDNA2_T4_LDS_I8_STRIDE;
  std::vector<int32_t> at((size_t)RDNA2_T4_TILE_MN * S);
  std::vector<int32_t> bt(at.size());
  // [wave][column group][lane][slot]
  std::vector<int32_t> accs((size_t)8 * RDNA2_T4_WMMA_NTILE * 32 * 8);
  int32_t af[32][2], bf[32][2];
  const int gx = (N + RDNA2_T4_TILE_MN - 1) / RDNA2_T4_TILE_MN;
  const int gy = (M + RDNA2_T4_TILE_MN - 1) / RDNA2_T4_TILE_MN;
  for (int by = 0; by < gy; ++by)
    for (int bx = 0; bx < gx; ++bx) {
      const int m_blk = by * RDNA2_T4_TILE_MN, n_blk = bx * RDNA2_T4_TILE_MN;
      for (auto& v : accs) v = 0;
      for (int k0 = 0; k0 < K4; k0 += RDNA2_T4_LDS_I8_KSTEP) {
        for (int tid = 0; tid < 256; ++tid) {
          sim_stage_i8(at.data(), A, M, K4, m_blk, k0, tid >> 1, (tid & 1) << 3);
          sim_stage_i8(at.data(), A, M, K4, m_blk, k0, tid >> 1, ((tid & 1) << 3) + 4);
          sim_stage_i8(bt.data(), Bt, N, K4, n_blk, k0, tid >> 1, (tid & 1) << 3);
          sim_stage_i8(bt.data(), Bt, N, K4, n_blk, k0, tid >> 1, ((tid & 1) << 3) + 4);
        }
        int kend = K4 - k0;
        if (kend > RDNA2_T4_LDS_I8_KSTEP) kend = RDNA2_T4_LDS_I8_KSTEP;
        const int ksteps = (kend + RDNA2_WMMA_I8_KSTEP_WORDS - 1) / RDNA2_WMMA_I8_KSTEP_WORDS;
        for (int wave = 0; wave < 8; ++wave)
          for (int ks = 0; ks < ksteps; ++ks) {
            for (int t = 0; t < 32; ++t) {
              const int w = rdna2_wmma_i8_frag_word(t, ks);
              const int r = wave * 16 + rdna2_wmma_i8_row(t);
              af[t][0] = at[(size_t)r * S + w];
              af[t][1] = at[(size_t)r * S + w + 1];
            }
            for (int j = 0; j < RDNA2_T4_WMMA_NTILE; ++j) {
              for (int t = 0; t < 32; ++t) {
                const int w = rdna2_wmma_i8_frag_word(t, ks);
                const int r = j * 16 + rdna2_wmma_i8_row(t);
                bf[t][0] = bt[(size_t)r * S + w];
                bf[t][1] = bt[(size_t)r * S + w + 1];
              }
              int32_t(*acc)[8] = reinterpret_cast<int32_t(*)[8]>(
                  &accs[((size_t)(wave * RDNA2_T4_WMMA_NTILE + j) * 32) * 8]);
              sim_wmma_i8_wave_step(acc, af, bf);
            }
          }
      }
      for (int wave = 0; wave < 8; ++wave)
        for (int j = 0; j < RDNA2_T4_WMMA_NTILE; ++j) {
          const int32_t(*acc)[8] = reinterpret_cast<const int32_t(*)[8]>(
              &accs[((size_t)(wave * RDNA2_T4_WMMA_NTILE + j) * 32) * 8]);
          for (int t = 0; t < 32; ++t) {
            const int n = n_blk + j * 16 + rdna2_wmma_i8_acc_col(t);
            if (n >= N) continue;
            for (int l = 0; l < 8; ++l) {
              const int m = m_blk + wave * 16 + rdna2_wmma_i8_acc_row(t, l);
              if (m >= M) break;  // rows 8*(t/16)+l are contiguous, break is safe
              C[(size_t)m * N + n] = acc[t][l];
            }
          }
        }
    }
}

static void sim_sgemm_lds_f32(float* C, const float* A, const float* Bt, int M, int N, int K) {
  std::vector<float> at((size_t)RDNA2_TILE_MN * RDNA2_LDS_F32_STRIDE);
  std::vector<float> bt(at.size());
  std::vector<float> accs((size_t)256 * 16);
  const int gx = (N + RDNA2_TILE_MN - 1) / RDNA2_TILE_MN;
  const int gy = (M + RDNA2_TILE_MN - 1) / RDNA2_TILE_MN;
  for (int by = 0; by < gy; ++by)
    for (int bx = 0; bx < gx; ++bx) {
      const int m_blk = by * RDNA2_TILE_MN, n_blk = bx * RDNA2_TILE_MN;
      for (auto& v : accs) v = 0.0f;
      for (int k0 = 0; k0 < K; k0 += RDNA2_LDS_F32_KSTEP) {
        for (int tid = 0; tid < 256; ++tid) {
          sim_stage_f32(at.data(), A, M, K, m_blk, k0, tid >> 2, (tid & 3) << 1);
          sim_stage_f32(bt.data(), Bt, N, K, n_blk, k0, tid >> 2, (tid & 3) << 1);
        }
        int kend = K - k0;
        if (kend > RDNA2_LDS_F32_KSTEP) kend = RDNA2_LDS_F32_KSTEP;
        for (int tid = 0; tid < 256; ++tid) {
          float(*acc)[4] = reinterpret_cast<float(*)[4]>(&accs[(size_t)tid * 16]);
          rdna2_sgemm_lds_block_f32(acc, at.data(), bt.data(), RDNA2_LDS_F32_STRIDE, (tid >> 4) * 4,
                                    (tid & 15) * 4, kend);
        }
      }
      for (int tid = 0; tid < 256; ++tid) {
        const int m0 = m_blk + (tid >> 4) * 4, n0 = n_blk + (tid & 15) * 4;
        const float(*acc)[4] = reinterpret_cast<const float(*)[4]>(&accs[(size_t)tid * 16]);
        for (int i = 0; i < 4 && m0 + i < M; ++i)
          for (int j = 0; j < 4 && n0 + j < N; ++j) C[(size_t)(m0 + i) * N + (n0 + j)] = acc[i][j];
      }
    }
}

// Variant 3: same 64x64 tile and 4x4 blocking as the kernel above, but 16
// floats staged per row per K step (four per thread) instead of two.
static void sim_sgemm_lds16_f32(float* C, const float* A, const float* Bt, int M, int N, int K) {
  std::vector<float> at((size_t)RDNA2_TILE_MN * RDNA2_LDS_F32B_STRIDE);
  std::vector<float> bt(at.size());
  std::vector<float> accs((size_t)256 * 16);
  const int gx = (N + RDNA2_TILE_MN - 1) / RDNA2_TILE_MN;
  const int gy = (M + RDNA2_TILE_MN - 1) / RDNA2_TILE_MN;
  for (int by = 0; by < gy; ++by)
    for (int bx = 0; bx < gx; ++bx) {
      const int m_blk = by * RDNA2_TILE_MN, n_blk = bx * RDNA2_TILE_MN;
      for (auto& v : accs) v = 0.0f;
      for (int k0 = 0; k0 < K; k0 += RDNA2_LDS_F32B_KSTEP) {
        for (int tid = 0; tid < 256; ++tid) {
          sim_stage_f32_4(at.data(), A, M, K, m_blk, k0, tid >> 2, (tid & 3) << 2);
          sim_stage_f32_4(bt.data(), Bt, N, K, n_blk, k0, tid >> 2, (tid & 3) << 2);
        }
        int kend = K - k0;
        if (kend > RDNA2_LDS_F32B_KSTEP) kend = RDNA2_LDS_F32B_KSTEP;
        for (int tid = 0; tid < 256; ++tid) {
          float(*acc)[4] = reinterpret_cast<float(*)[4]>(&accs[(size_t)tid * 16]);
          rdna2_sgemm_lds_block_f32(acc, at.data(), bt.data(), RDNA2_LDS_F32B_STRIDE, (tid >> 4) * 4,
                                    (tid & 15) * 4, kend);
        }
      }
      for (int tid = 0; tid < 256; ++tid) {
        const int m0 = m_blk + (tid >> 4) * 4, n0 = n_blk + (tid & 15) * 4;
        const float(*acc)[4] = reinterpret_cast<const float(*)[4]>(&accs[(size_t)tid * 16]);
        for (int i = 0; i < 4 && m0 + i < M; ++i)
          for (int j = 0; j < 4 && n0 + j < N; ++j) C[(size_t)(m0 + i) * N + (n0 + j)] = acc[i][j];
      }
    }
}

// f32 counterpart of the 128x128 tile: two staged floats per thread per call,
// four calls per K step, 8x8 accumulators.
static void sim_sgemm_lds8x8_f32(float* C, const float* A, const float* Bt, int M, int N, int K) {
  std::vector<float> at((size_t)RDNA2_T4_TILE_MN * RDNA2_T4_LDS_F32_STRIDE);
  std::vector<float> bt(at.size());
  std::vector<float> accs((size_t)256 * 64);
  const int gx = (N + RDNA2_T4_TILE_MN - 1) / RDNA2_T4_TILE_MN;
  const int gy = (M + RDNA2_T4_TILE_MN - 1) / RDNA2_T4_TILE_MN;
  for (int by = 0; by < gy; ++by)
    for (int bx = 0; bx < gx; ++bx) {
      const int m_blk = by * RDNA2_T4_TILE_MN, n_blk = bx * RDNA2_T4_TILE_MN;
      for (auto& v : accs) v = 0.0f;
      for (int k0 = 0; k0 < K; k0 += RDNA2_T4_LDS_F32_KSTEP) {
        for (int tid = 0; tid < 256; ++tid) {
          sim_stage_f32(at.data(), A, M, K, m_blk, k0, tid >> 1, (tid & 1) << 2);
          sim_stage_f32(at.data(), A, M, K, m_blk, k0, tid >> 1, ((tid & 1) << 2) + 2);
          sim_stage_f32(bt.data(), Bt, N, K, n_blk, k0, tid >> 1, (tid & 1) << 2);
          sim_stage_f32(bt.data(), Bt, N, K, n_blk, k0, tid >> 1, ((tid & 1) << 2) + 2);
        }
        int kend = K - k0;
        if (kend > RDNA2_T4_LDS_F32_KSTEP) kend = RDNA2_T4_LDS_F32_KSTEP;
        for (int tid = 0; tid < 256; ++tid) {
          float(*acc)[8] = reinterpret_cast<float(*)[8]>(&accs[(size_t)tid * 64]);
          rdna2_sgemm_lds_block8x8_f32(acc, at.data(), bt.data(), RDNA2_T4_LDS_F32_STRIDE,
                                       (tid >> 4) * 8, (tid & 15) * 8, kend);
        }
      }
      for (int tid = 0; tid < 256; ++tid) {
        const int m0 = m_blk + (tid >> 4) * 8, n0 = n_blk + (tid & 15) * 8;
        const float(*acc)[8] = reinterpret_cast<const float(*)[8]>(&accs[(size_t)tid * 64]);
        for (int i = 0; i < 8 && m0 + i < M; ++i)
          for (int j = 0; j < 8 && n0 + j < N; ++j) C[(size_t)(m0 + i) * N + (n0 + j)] = acc[i][j];
      }
    }
}

int main() {
  char devname[256];
  int warp_size = 0, arch_ok = 0;
  if (rdna2_device_info(devname, sizeof(devname), &warp_size, &arch_ok) != 0) {
    std::printf("device query failed\n");
    return 1;
  }
  const int has_device = (devname[0] != '\0');
  if (has_device) {
    const char* arch = arch_ok == RDNA2_ARCH_GFX1031  ? "gfx1031 (RDNA2 target)"
                       : arch_ok == RDNA2_ARCH_GFX1030 ? "gfx1030 (RDNA2 target)"
                       : arch_ok == RDNA2_ARCH_GFX1201 ? "gfx1201 (RDNA4, 9070 XT)"
                       : arch_ok == RDNA2_ARCH_GFX1200 ? "gfx1200 (RDNA4, 9060 XT)"
                                                      : "unknown (generic 64x64 tile)";
    std::printf("device: %s  warp=%d  %s\n", devname, warp_size, arch);
    std::printf("arch: code=%d  gemm i8 tiled=%d  sgemm variant=%d\n", arch_ok,
                rdna2_tiled_for_arch(arch_ok), rdna2_sgemm_variant_for_arch(arch_ok));
  } else {
    std::printf("device: none - running CPU fallback paths only\n");
    warp_size = 32;
  }
  std::printf("path: %s\n", has_device ? "GPU (launchers use the device)" : "CPU fallback");

  // --- 1. Argument validation (must hold on every path) -------------------
  {
    std::vector<int32_t> sa, sx, sy;
    int32_t* A = align16(sa, 256);
    int32_t* x = align16(sx, 16);
    int32_t* y = align16(sy, 8);
    CHECK(rdna2_gemv_i8_host(nullptr, A, x, 4, 4) == 2, "null y rejected");
    CHECK(rdna2_gemv_i8_host(y, A, x, 0, 4) == 2, "rows=0 rejected");
    CHECK(rdna2_gemv_i8_host(y, A, x, 4, 0) == 2, "k4=0 rejected");
    CHECK(rdna2_gemv_i8_host(nullptr, nullptr, nullptr, 4, 4) == 2, "gemv null rejected");
    // k4 not a multiple of 4 and odd host alignment are valid input: the kernel
    // has a scalar tail and the launcher copies into an aligned device buffer.
    CHECK(rdna2_gemv_i8_host(y, A, x, 4, 3) == 0, "k4=3 accepted (scalar tail)");
    CHECK(rdna2_gemv_i8_host(y, align16(sa, 256, 4), x, 4, 4) == 0, "misaligned A accepted");
    CHECK(rdna2_gemv_i8_host(y, A, align16(sx, 16, 8), 4, 4) == 0, "misaligned x accepted");
    CHECK(rdna2_gemv_i8_host(y, A, x, 4, 4) == 0, "valid gemv accepted");

    std::vector<float> fA(256), fX(16), fY(16);
    CHECK(rdna2_sgemv_f32_host(nullptr, fA.data(), fX.data(), 4, 4) == 2,
          "sgemv null y rejected");
    CHECK(rdna2_sgemv_f32_host(fY.data(), fA.data(), fX.data(), 0, 4) == 2, "sgemv rows=0 rejected");
    CHECK(rdna2_sgemv_f32_host(fY.data(), fA.data(), fX.data(), 4, 0) == 2, "sgemv K=0 rejected");
    CHECK(rdna2_sgemv_f32_host(fY.data(), fA.data(), fX.data(), 4, 4) == 0, "valid sgemv accepted");
    CHECK(rdna2_sgemv_f32_host_variant(fY.data(), fA.data(), fX.data(), 4, 4, 1) == 0,
          "valid sgemv wave accepted");

    int32_t c = 0;
    CHECK(rdna2_gemm_i8_host(&c, A, A, 0, 1, 1, 1) == 2, "gemm M=0 rejected");
    CHECK(rdna2_sgemm_f32_host(nullptr, nullptr, nullptr, 1, 1, 1) == 2, "sgemm null rejected");
    CHECK(rdna2_dot_u8_host(nullptr, nullptr, nullptr, 4) == 2, "dot_u8 null rejected");
    CHECK(rdna2_prefix_host(y, x, 0) == 2, "prefix n=0 rejected");
    CHECK(rdna2_prefix_host(y, x, warp_size + 1) == 2, "prefix n>warp rejected");
  }

  // --- 1b. Architecture selectors -----------------------------------------
  // Both RDNA4 parts take the 128x128 matrix-instruction shape for int8 and the
  // 16-deep f32 K step; both RDNA2 parts take the 64x64 dp4a tile (no matrix
  // unit there), and anything else the 4x4 tiled kernel with global loads. The
  // 128-tile dp4a kernel measured slower than the 64x64 one on a 9070 XT (see
  // README), so tiled == 4 is exercised explicitly rather than by selector.
  {
    CHECK(rdna2_tiled_for_arch(RDNA2_ARCH_GFX1201) == 5, "9070 XT takes the wmma tile");
    CHECK(rdna2_tiled_for_arch(RDNA2_ARCH_GFX1200) == 5, "9060 XT takes the wmma tile");
    CHECK(rdna2_tiled_for_arch(RDNA2_ARCH_GFX1031) == 3, "gfx1031 keeps the 64x64 8x4 tile");
    CHECK(rdna2_tiled_for_arch(RDNA2_ARCH_GFX1030) == 3, "gfx1030 keeps the 64x64 8x4 tile");
    CHECK(rdna2_tiled_for_arch(RDNA2_ARCH_UNKNOWN) == 1, "unknown GPU falls back to 4x4 tiled");
    CHECK(rdna2_sgemm_variant_for_arch(RDNA2_ARCH_GFX1201) == 3, "9070 XT takes the 16-deep f32 K step");
    CHECK(rdna2_sgemm_variant_for_arch(RDNA2_ARCH_GFX1200) == 3, "9060 XT takes the 16-deep f32 K step");
    CHECK(rdna2_sgemm_variant_for_arch(RDNA2_ARCH_GFX1031) == 3, "gfx1031 takes the 16-deep f32 K step");
    CHECK(rdna2_sgemm_variant_for_arch(RDNA2_ARCH_GFX1030) == 3, "gfx1030 takes the 16-deep f32 K step");
    CHECK(rdna2_sgemm_variant_for_arch(RDNA2_ARCH_UNKNOWN) == 1, "unknown GPU keeps the 8-deep f32 K step");
    // The tuned entry points query the running device; with no GPU they take the
    // CPU reference path, which must still agree with the plain entry points.
    const int M = 17, N = 9, K = 6, K4 = 3;
    std::vector<int32_t> sA(M * K4), sB(N * K4), rI(M * N), gI(M * N);
    std::vector<float> fA(M * K), fB(N * K), rF(M * N), gF(M * N);
    for (int i = 0; i < M * K4; ++i) sA[i] = (int32_t)rng();
    for (int i = 0; i < N * K4; ++i) sB[i] = (int32_t)rng();
    for (int i = 0; i < M * K; ++i) fA[i] = (float)((int)(rng() % 2001) - 1000) / 64.0f;
    for (int i = 0; i < N * K; ++i) fB[i] = (float)((int)(rng() % 2001) - 1000) / 64.0f;
    rdna2_gemm_i8_ref(rI.data(), sA.data(), sB.data(), M, N, K4);
    rdna2_sgemm_f32_ref(rF.data(), fA.data(), fB.data(), M, N, K);
    CHECK(rdna2_gemm_i8_host_tuned(gI.data(), sA.data(), sB.data(), M, N, K4) == 0,
          "tuned gemm returned error");
    CHECK(rdna2_sgemm_f32_host_tuned(gF.data(), fA.data(), fB.data(), M, N, K) == 0,
          "tuned sgemm returned error");
    int badI = 0, badF = 0;
    for (int i = 0; i < M * N; ++i) {
      badI += (gI[i] != rI[i]);
      badF += (std::fabs((double)gF[i] - (double)rF[i]) > 1e-4 * (std::fabs((double)rF[i]) + 1.0));
    }
    CHECK(badI == 0, "tuned gemm matches reference");
    CHECK(badF == 0, "tuned sgemm matches reference");
  }

  // --- 2. GEMV, dp4a int8 (8 words = 32 int8 per row) ---------------------
  // Three cases: k4 a multiple of 4 (vector path only), k4=5 (vector + tail),
  // and a deliberately misaligned host buffer (launcher copies, so still exact).
  for (int pass = 0; pass < 3; ++pass) {
    const int rows = 37;
    const int k4 = (pass == 0) ? 8 : 5;
    const size_t misalign = (pass == 2) ? 4 : 0;
    std::vector<int32_t> sa, sx, sy;
    int32_t* A = align16(sa, (size_t)rows * k4, misalign);
    int32_t* x = align16(sx, k4, (pass == 2) ? 12 : 0);
    int32_t* y = align16(sy, rows);
    for (int i = 0; i < rows * k4; ++i) A[i] = (int32_t)rng();
    for (int i = 0; i < k4; ++i) x[i] = (int32_t)rng();
    std::vector<int32_t> ref(rows);
    rdna2_gemv_i8_ref(ref.data(), A, x, rows, k4);
    CHECK(rdna2_gemv_i8_host(y, A, x, rows, k4) == 0, "gemv returned error");
    int bad = 0;
    for (int r = 0; r < rows; ++r) bad += (y[r] != ref[r]);
    CHECK(bad == 0, "gemv matches reference");
    // Same body the kernel runs, driven by a host replica of its grid.
    std::vector<int32_t> sim(rows, SENTINEL);
    sim_gemv_i8(sim.data(), A, x, rows, k4);
    int sim_bad = 0;
    for (int r = 0; r < rows; ++r) sim_bad += (sim[r] != ref[r]);
    CHECK(sim_bad == 0, "gemv kernel body matches reference");
    std::printf("gemv   i8  rows=%-4d k4=%-4d misalign=%-2d mismatch=%d sim=%d\n", rows, k4,
                (int)misalign, bad, sim_bad);
  }

  // --- 2b. GEMV, wave-per-row variant (what k4 >= 8 selects by default) ---
  // Narrows down to two cases here because the point is the cooperative row:
  // a short row (k4=5, several lanes idle) and a row that is a multiple of the
  // wave width (k4=8 gives the vector loads plenty of work).
  for (int pass = 0; pass < 2; ++pass) {
    const int rows = 37;
    const int k4 = (pass == 0) ? 8 : 5;
    std::vector<int32_t> sa, sx, sy;
    int32_t* A = align16(sa, (size_t)rows * k4);
    int32_t* x = align16(sx, k4);
    int32_t* y = align16(sy, rows);
    for (int i = 0; i < rows * k4; ++i) A[i] = (int32_t)rng();
    for (int i = 0; i < k4; ++i) x[i] = (int32_t)rng();
    std::vector<int32_t> ref(rows);
    rdna2_gemv_i8_ref(ref.data(), A, x, rows, k4);
    CHECK(rdna2_gemv_i8_host_variant(y, A, x, rows, k4, 1) == 0, "gemv wave returned error");
    int bad = 0;
    for (int r = 0; r < rows; ++r) bad += (y[r] != ref[r]);
    CHECK(bad == 0, "gemv wave variant matches reference");
    std::vector<int32_t> sim(rows, SENTINEL);
    sim_gemv_wave_i8(sim.data(), A, x, rows, k4, warp_size);
    int sim_bad = 0;
    for (int r = 0; r < rows; ++r) sim_bad += (sim[r] != ref[r]);
    CHECK(sim_bad == 0, "gemv wave body matches reference");
    std::printf("gemv   i8w rows=%-4d k4=%-4d wave=%d mismatch=%d sim=%d\n", rows, k4, warp_size, bad,
                sim_bad);
  }

  // --- 2c. GEMV, f32: both variants plus their host replicas ---------------
  // pass 0: K a multiple of 4 and >= the wave threshold (wave variant), pass 1:
  // short row (thread variant), pass 2: odd K with an unaligned A view.
  for (int pass = 0; pass < 3; ++pass) {
    const int rows = 37;
    const int K = (pass == 0) ? 128 : (pass == 1 ? 5 : 133);
    const size_t shift = (pass == 2) ? 1 : 0;
    std::vector<float> A((size_t)rows * K + shift), x((size_t)K), ref((size_t)rows);
    for (auto& v : A) v = (float)((int)(rng() % 2001) - 1000) / 512.0f;
    for (auto& v : x) v = (float)((int)(rng() % 2001) - 1000) / 512.0f;
    const float* Av = A.data() + shift;
    rdna2_sgemv_f32_ref(ref.data(), Av, x.data(), rows, K);
    const int want = K >= 128 ? 1 : 0;
    for (int variant = 0; variant < 2; ++variant) {
      std::vector<float> y((size_t)rows, 1e30f);
      CHECK(rdna2_sgemv_f32_host_variant(y.data(), Av, x.data(), rows, K, variant) == 0,
            "sgemv returned error");
      double worst = 0.0;
      for (int r = 0; r < rows; ++r) {
        const double s = std::fabs((double)ref[r]) > 1e-6 ? std::fabs((double)ref[r]) : 1.0;
        worst = std::fmax(worst, std::fabs((double)y[r] - (double)ref[r]) / s);
      }
      CHECK(worst < 1e-6, "sgemv variant matches reference");
      std::vector<float> sim((size_t)rows, 1e30f);
      if (variant == 1)
        sim_sgemv_wave_f32(sim.data(), Av, x.data(), rows, K, warp_size);
      else
        sim_sgemv_f32(sim.data(), Av, x.data(), rows, K);
      double sim_worst = 0.0;
      for (int r = 0; r < rows; ++r) {
        const double s = std::fabs((double)ref[r]) > 1e-6 ? std::fabs((double)ref[r]) : 1.0;
        sim_worst = std::fmax(sim_worst, std::fabs((double)sim[r] - (double)ref[r]) / s);
      }
      CHECK(sim_worst < 1e-6, "sgemv kernel body matches reference");
      std::printf("sgemv f32 rows=%-4d K=%-4d variant=%d mismatch=%.1e sim=%.1e\n", rows, K, variant,
                  worst, sim_worst);
    }
    if (has_device) {
      std::vector<float> y((size_t)rows, 1e30f);
      CHECK(rdna2_sgemv_f32_host(y.data(), Av, x.data(), rows, K) == 0, "sgemv default path failed");
      double worst = 0.0;
      for (int r = 0; r < rows; ++r) {
        const double s = std::fabs((double)ref[r]) > 1e-6 ? std::fabs((double)ref[r]) : 1.0;
        worst = std::fmax(worst, std::fabs((double)y[r] - (double)ref[r]) / s);
      }
      CHECK(worst < 1e-6, "sgemv default variant matches reference");
      CHECK((K >= 128 ? 1 : 0) == want, "sgemv threshold selects the wave for long rows");
    }
  }

  // --- 3. GEMM int8: naive, tiled and kernel-body sim, crossing the 64x64
  //        tile boundary with partial edge tiles -------------------------
  {
    const int M = 70, N = 67, K4 = 10;
    std::vector<int32_t> sA(M * K4), sB(N * K4), ref(M * N), cN(M * N), cT(M * N);
    for (auto& v : sA) v = (int32_t)rng();
    for (auto& v : sB) v = (int32_t)rng();
    rdna2_gemm_i8_ref(ref.data(), sA.data(), sB.data(), M, N, K4);
    cN.assign(M * N, SENTINEL);
    cT.assign(M * N, SENTINEL);
    CHECK(rdna2_gemm_i8_host(cN.data(), sA.data(), sB.data(), M, N, K4, 0) == 0,
          "gemm naive returned error");
    CHECK(rdna2_gemm_i8_host(cT.data(), sA.data(), sB.data(), M, N, K4, 1) == 0,
          "gemm tiled returned error");
    std::vector<int32_t> cL(M * N, SENTINEL);
    CHECK(rdna2_gemm_i8_host(cL.data(), sA.data(), sB.data(), M, N, K4, 2) == 0,
          "gemm lds returned error");
    int badN = 0, badT = 0, badL = 0;
    for (int i = 0; i < M * N; ++i) {
      badN += (cN[i] != ref[i]);
      badT += (cT[i] != ref[i]);
      badL += (cL[i] != ref[i]);
    }
    CHECK(badN == 0, "gemm naive matches reference");
    CHECK(badT == 0, "gemm tiled matches reference");
    CHECK(badL == 0, "gemm lds matches reference");
    // K4=10 stages one short chunk; K4=40 stages 16+16+8, so restaging and the
    // partial last chunk are exercised on odd matrix sizes too.
    std::vector<int32_t> lds_sim(M * N, SENTINEL);
    sim_gemm_lds_i8(lds_sim.data(), sA.data(), sB.data(), M, N, K4);
    int badLS = 0;
    for (int i = 0; i < M * N; ++i) badLS += (lds_sim[i] != ref[i]);
    CHECK(badLS == 0, "gemm lds simulation covers every element correctly");
    std::vector<int32_t> cL8(M * N, SENTINEL), lds8_sim(M * N, SENTINEL);
    CHECK(rdna2_gemm_i8_host(cL8.data(), sA.data(), sB.data(), M, N, K4, 3) == 0,
          "gemm lds8 returned error");
    int badL8 = 0;
    for (int i = 0; i < M * N; ++i) badL8 += (cL8[i] != ref[i]);
    CHECK(badL8 == 0, "gemm lds8 matches reference");
    sim_gemm_lds8_i8(lds8_sim.data(), sA.data(), sB.data(), M, N, K4);
    int badL8S = 0;
    for (int i = 0; i < M * N; ++i) badL8S += (lds8_sim[i] != ref[i]);
    CHECK(badL8S == 0, "gemm lds8 simulation covers every element correctly");
    // RDNA4 shape: big enough for two 128x128 tiles per dimension (partial edge
    // tiles) and K4=40, which stages 16+16+8.
    {
      const int M4 = 135, N4 = 130, K4b = 40;
      std::vector<int32_t> A4((size_t)M4 * K4b), B4((size_t)N4 * K4b), ref4((size_t)M4 * N4),
          got4((size_t)M4 * N4, SENTINEL), sim4((size_t)M4 * N4, SENTINEL);
      for (auto& v : A4) v = (int32_t)rng();
      for (auto& v : B4) v = (int32_t)rng();
      rdna2_gemm_i8_ref(ref4.data(), A4.data(), B4.data(), M4, N4, K4b);
      CHECK(rdna2_gemm_i8_host(got4.data(), A4.data(), B4.data(), M4, N4, K4b, 4) == 0,
            "gemm lds128 returned error");
      int b4 = 0;
      for (int i = 0; i < M4 * N4; ++i) b4 += (got4[i] != ref4[i]);
      CHECK(b4 == 0, "gemm lds128 matches reference");
      sim_gemm_lds8x8_i8(sim4.data(), A4.data(), B4.data(), M4, N4, K4b);
      int b4s = 0;
      for (int i = 0; i < M4 * N4; ++i) b4s += (sim4[i] != ref4[i]);
      CHECK(b4s == 0, "gemm lds128 simulation covers every element correctly");
      std::printf("gemm i8ld128 M=%-4d N=%-4d K4=%-3d mismatch=%d sim=%d (128x128, 8x8)\n", M4, N4,
                  K4b, b4, b4s);
    }
    // RDNA4 matrix path: the same 128x128 tile driven by 16x16x16 WMMA. M4/N4
    // are not multiples of 128 and K4b=40 stages 16+16+8, so edge tiles and a
    // partial last chunk are covered here too.
    {
      const int M4 = 135, N4 = 130, K4b = 40;
      std::vector<int32_t> A4((size_t)M4 * K4b), B4((size_t)N4 * K4b), ref4((size_t)M4 * N4),
          got5((size_t)M4 * N4, SENTINEL), sim5((size_t)M4 * N4, SENTINEL);
      for (auto& v : A4) v = (int32_t)rng();
      for (auto& v : B4) v = (int32_t)rng();
      rdna2_gemm_i8_ref(ref4.data(), A4.data(), B4.data(), M4, N4, K4b);
      CHECK(rdna2_gemm_i8_host(got5.data(), A4.data(), B4.data(), M4, N4, K4b, 5) == 0,
            "gemm wmma returned error");
      int b5 = 0;
      for (int i = 0; i < M4 * N4; ++i) b5 += (got5[i] != ref4[i]);
      CHECK(b5 == 0, "gemm wmma matches reference");
      sim_gemm_wmma_i8(sim5.data(), A4.data(), B4.data(), M4, N4, K4b);
      int b5s = 0;
      for (int i = 0; i < M4 * N4; ++i) b5s += (sim5[i] != ref4[i]);
      CHECK(b5s == 0, "gemm wmma simulation covers every element correctly");
      std::printf("gemm i8wmma M=%-4d N=%-4d K4=%-3d mismatch=%d sim=%d (128x128, wmma)\n", M4, N4,
                  K4b, b5, b5s);
    }
    // Fragment placement on its own: a single 16x16x16 step with distinct byte
    // values must produce the scalar dot product in every one of the 256 output
    // slots, and must add into a nonzero accumulator rather than overwrite it.
    {
      int8_t ta[16][16], tb[16][16];
      int32_t acc[32][8], af[32][2], bf[32][2];
      for (int m = 0; m < 16; ++m)
        for (int k = 0; k < 16; ++k) {
          ta[m][k] = (int8_t)((int)(m * 7 + k * 3) % 251 - 125);
          tb[k][m] = (int8_t)((int)(m * 5 + k * 11) % 241 - 120);
        }
      for (int t = 0; t < 32; ++t) {
        for (int e = 0; e < 2; ++e) {
          uint32_t wa = 0, wb = 0;
          for (int b = 0; b < 4; ++b) {
            const int k = (t / 16) * 8 + e * 4 + b;
            wa |= (uint32_t)(uint8_t)ta[rdna2_wmma_i8_row(t)][k] << (8 * b);
            wb |= (uint32_t)(uint8_t)tb[k][rdna2_wmma_i8_row(t)] << (8 * b);
          }
          af[t][e] = (int32_t)wa;  // Bt rows are K-major, so tb[k][n] as built
          bf[t][e] = (int32_t)wb;
        }
        for (int l = 0; l < 8; ++l) acc[t][l] = 1000;  // nonzero, must survive
      }
      sim_wmma_i8_wave_step(acc, af, bf);
      int badf = 0;
      for (int t = 0; t < 32; ++t)
        for (int l = 0; l < 8; ++l) {
          const int m = rdna2_wmma_i8_acc_row(t, l), n = rdna2_wmma_i8_acc_col(t);
          int32_t want = 1000;
          for (int k = 0; k < 16; ++k) want += (int32_t)ta[m][k] * (int32_t)tb[k][n];
          badf += (acc[t][l] != want);
        }
      CHECK(badf == 0, "wmma fragment placement maps every result slot");
      std::printf("wmma i8 fragment 16x16x16 mismatch=%d (per-slot dot + accumulate)\n", badf);
    }
    {
      const int M2 = 35, N2 = 33, K2 = 40;
      std::vector<int32_t> A2(M2 * K2), B2(N2 * K2), ref2(M2 * N2), got2(M2 * N2, SENTINEL),
          sim2(M2 * N2, SENTINEL);
      for (auto& v : A2) v = (int32_t)rng();
      for (auto& v : B2) v = (int32_t)rng();
      rdna2_gemm_i8_ref(ref2.data(), A2.data(), B2.data(), M2, N2, K2);
      CHECK(rdna2_gemm_i8_host(got2.data(), A2.data(), B2.data(), M2, N2, K2, 2) == 0,
            "gemm lds multi-chunk returned error");
      int b2 = 0;
      for (int i = 0; i < M2 * N2; ++i) b2 += (got2[i] != ref2[i]);
      CHECK(b2 == 0, "gemm lds multi-chunk matches reference");
      sim_gemm_lds_i8(sim2.data(), A2.data(), B2.data(), M2, N2, K2);
      int b2s = 0;
      for (int i = 0; i < M2 * N2; ++i) b2s += (sim2[i] != ref2[i]);
      CHECK(b2s == 0, "gemm lds multi-chunk simulation matches reference");
      CHECK(badL == 0 && b2s == 0 && b2 == 0, "gemm lds chunk loop is correct");
      // The 8x4 variant stages two words per thread, so it has its own coverage
      // of the partial chunk and the edge tiles.
      std::vector<int32_t> got8(M2 * N2, SENTINEL), sim8(M2 * N2, SENTINEL);
      CHECK(rdna2_gemm_i8_host(got8.data(), A2.data(), B2.data(), M2, N2, K2, 3) == 0,
            "gemm lds8 multi-chunk returned error");
      int b8 = 0;
      for (int i = 0; i < M2 * N2; ++i) b8 += (got8[i] != ref2[i]);
      sim_gemm_lds8_i8(sim8.data(), A2.data(), B2.data(), M2, N2, K2);
      int b8s = 0;
      for (int i = 0; i < M2 * N2; ++i) b8s += (sim8[i] != ref2[i]);
      CHECK(b8 == 0 && b8s == 0, "gemm lds8 multi-chunk matches reference");
      std::printf("gemm i8lds M=%-4d N=%-4d K4=%-4d mismatch=%d sim=%d (partial chunk)\n", M2, N2,
                  K2, b2, b2s);
      std::printf("gemm i8ld8 M=%-4d N=%-4d K4=%-4d mismatch=%d sim=%d (8x4 blocking)\n", M2, N2,
                  K2, b8, b8s);
    }
    // Kernel body on the host, seeded so a missing or stray write is caught.
    std::vector<int32_t> cS(M * N, SENTINEL);
    sim_gemm_tile_i8(cS.data(), sA.data(), sB.data(), M, N, K4);
    int badS = 0;
    for (int i = 0; i < M * N; ++i) badS += (cS[i] != ref[i]);
    CHECK(badS == 0, "gemm tiled kernel body covers every element correctly");
    std::printf("gemm   i8  M=%-4d N=%-4d K4=%-4d naive=%d tiled=%d lds=%d lds8=%d sim=%d lds_sim=%d lds8_sim=%d\n",
                M, N, K4, badN, badT, badL, badL8, badS, badLS, badL8S);
  }

  // --- 3b. Quantized epilogue: scale and zero point applied to the int32
  //         accumulator, with the stored quants left untouched ----------------
  {
    const int M = 37, N = 45, K4 = 12;  // odd: partial 64x64 and 128x128 tiles
    std::vector<int32_t> sA((size_t)M * K4), sB((size_t)N * K4);
    for (auto& v : sA) v = (int32_t)rng();
    for (auto& v : sB) v = (int32_t)rng();
    const std::vector<int32_t> A0 = sA, B0 = sB;
    std::vector<float> sa(M), sb(N), ref((size_t)M * N), got((size_t)M * N);
    for (int i = 0; i < M; ++i) sa[i] = 0.001f * (float)(1 + i % 7);
    for (int i = 0; i < N; ++i) sb[i] = 0.002f * (float)(1 + i % 5);

    auto worst_rel = [&]() {
      double worst = 0.0;
      for (size_t i = 0; i < ref.size(); ++i) {
        const double d = std::fabs((double)got[i] - (double)ref[i]);
        const double s = std::fabs((double)ref[i]) > 1e-6 ? std::fabs((double)ref[i]) : 1e-6;
        if (d / s > worst) worst = d / s;
      }
      return worst;
    };

    struct ScaleCase {
      const char* name;
      const float* a;
      const float* b;
      int zp;
      int tiled;
    };
    const ScaleCase cases[] = {
        {"unit scales", nullptr, nullptr, 0, 3},
        {"row scale", sa.data(), nullptr, 0, 3},
        {"column scale", nullptr, sb.data(), 0, 3},
        {"both scales", sa.data(), sb.data(), 0, 3},
        {"zero point", sa.data(), sb.data(), 8, 3},
        {"zero point, 128 tile", sa.data(), sb.data(), 8, 4},
        {"zero point, wmma", sa.data(), sb.data(), 8, 5},
        {"zero point, naive", sa.data(), sb.data(), 8, 0},
    };
    for (const ScaleCase& c : cases) {
      rdna2_gemm_i8_scale_ref(ref.data(), sA.data(), sB.data(), M, N, K4, c.a, c.b, c.zp);
      got.assign(got.size(), -1.0f);
      CHECK(rdna2_gemm_i8_scale_host(got.data(), sA.data(), sB.data(), M, N, K4, c.tiled, c.a, c.b,
                                     c.zp) == 0,
            "scaled gemm returned error");
      const double worst = worst_rel();
      CHECK(worst < 1e-5, "scaled gemm matches reference");
      std::printf("gemm i8+scale %-21s zp=%d tiled=%d worst_rel=%.2e\n", c.name, c.zp, c.tiled,
                  worst);
    }

    // The correction must be exact, not merely close: with unit scales nothing
    // rounds, so (acc - zp*rowsum) has to come back as that exact integer.
    {
      rdna2_gemm_i8_scale_ref(ref.data(), sA.data(), sB.data(), M, N, K4, nullptr, nullptr, 8);
      CHECK(rdna2_gemm_i8_scale_host(got.data(), sA.data(), sB.data(), M, N, K4, 3, nullptr, nullptr,
                                     8) == 0,
            "zero point exactness returned error");
      int bad = 0;
      for (size_t i = 0; i < ref.size(); ++i) bad += (got[i] != ref[i]);
      CHECK(bad == 0, "zero point correction is bit-exact");
    }

    // Loaded quants are read-only. The paths that dequantise in place, or that
    // expand to a float copy, would show up here as a changed input buffer.
    CHECK(std::memcmp(sA.data(), A0.data(), sA.size() * sizeof(int32_t)) == 0 &&
              std::memcmp(sB.data(), B0.data(), sB.size() * sizeof(int32_t)) == 0,
          "scaled gemm leaves the stored quants untouched");

    CHECK(rdna2_gemm_i8_scale_host(nullptr, sA.data(), sB.data(), M, N, K4, 3, nullptr, nullptr, 0) ==
              2,
          "scaled gemm rejects null C");
    CHECK(rdna2_gemm_i8_scale_host(got.data(), sA.data(), sB.data(), 0, N, K4, 3, nullptr, nullptr,
                                   0) == 2,
          "scaled gemm rejects M=0");

    // The arch-aware entry point must agree with the explicit one on this GPU.
    {
      std::vector<float> got2((size_t)M * N);
      rdna2_gemm_i8_scale_ref(ref.data(), sA.data(), sB.data(), M, N, K4, sa.data(), sb.data(), 8);
      CHECK(rdna2_gemm_i8_scale_host_tuned(got2.data(), sA.data(), sB.data(), M, N, K4, sa.data(),
                                           sb.data(), 8) == 0,
            "tuned scaled gemm returned error");
      int bad = 0;
      for (size_t i = 0; i < ref.size(); ++i) bad += (got2[i] != ref[i]);
      CHECK(bad == 0, "tuned scaled gemm matches reference");
    }
  }

  // --- 4. SGEMM f32, also across the tile boundary ------------------------
  {
    const int M = 67, N = 70, K = 19;
    std::vector<float> A(M * K), B(N * K), ref(M * N), c(M * N, 1e30f);
    for (auto& v : A) v = (float)((int)(rng() % 2001) - 1000) / 512.0f;
    for (auto& v : B) v = (float)((int)(rng() % 2001) - 1000) / 512.0f;
    rdna2_sgemm_f32_ref(ref.data(), A.data(), B.data(), M, N, K);
    CHECK(rdna2_sgemm_f32_host(c.data(), A.data(), B.data(), M, N, K) == 0,
          "sgemm returned error");
    double max_rel = 0.0;
    for (int i = 0; i < M * N; ++i) {
      double d = std::fabs((double)c[i] - (double)ref[i]);
      double s = std::fabs((double)ref[i]) > 1e-6 ? std::fabs((double)ref[i]) : 1.0;
      max_rel = std::fmax(max_rel, d / s);
    }
    // fma vs separate multiply/add can differ in the last bit.
    CHECK(max_rel < 1e-5, "sgemm matches reference within tolerance");
    std::vector<float> cs(M * N, 1e30f);
    sim_sgemm_tile_f32(cs.data(), A.data(), B.data(), M, N, K);
    double max_sim = 0.0;
    for (int i = 0; i < M * N; ++i) {
      double d = std::fabs((double)cs[i] - (double)ref[i]);
      double s = std::fabs((double)ref[i]) > 1e-6 ? std::fabs((double)ref[i]) : 1.0;
      max_sim = std::fmax(max_sim, d / s);
    }
    CHECK(max_sim < 1e-5, "sgemm kernel body matches reference within tolerance");
    // Both device variants explicitly: 0 = global loads, 1 = LDS staged.
    std::vector<float> c0(M * N, 1e30f), c1(M * N, 1e30f);
    CHECK(rdna2_sgemm_f32_host_v(c0.data(), A.data(), B.data(), M, N, K, 0) == 0,
          "sgemm variant 0 returned error");
    CHECK(rdna2_sgemm_f32_host_v(c1.data(), A.data(), B.data(), M, N, K, 1) == 0,
          "sgemm variant 1 returned error");
    double rel0 = 0.0, rel1 = 0.0;
    for (int i = 0; i < M * N; ++i) {
      double s = std::fabs((double)ref[i]) > 1e-6 ? std::fabs((double)ref[i]) : 1.0;
      rel0 = std::fmax(rel0, std::fabs((double)c0[i] - (double)ref[i]) / s);
      rel1 = std::fmax(rel1, std::fabs((double)c1[i] - (double)ref[i]) / s);
    }
    CHECK(rel0 < 1e-5, "sgemm variant 0 matches reference");
    CHECK(rel1 < 1e-5, "sgemm variant 1 matches reference");
    std::vector<float> cls(M * N, 1e30f);
    sim_sgemm_lds_f32(cls.data(), A.data(), B.data(), M, N, K);
    double max_ls = 0.0;
    for (int i = 0; i < M * N; ++i) {
      double d = std::fabs((double)cls[i] - (double)ref[i]);
      double s = std::fabs((double)ref[i]) > 1e-6 ? std::fabs((double)ref[i]) : 1.0;
      max_ls = std::fmax(max_ls, d / s);
    }
    CHECK(max_ls < 1e-5, "sgemm lds simulation matches reference within tolerance");
    std::printf("sgemm  f32 M=%-4d N=%-4d K=%-5d max_rel=%.2e sim_rel=%.2e v0=%.2e lds_sim=%.2e\n",
                M, N, K, max_rel, max_sim, rel0 > rel1 ? rel0 : rel1, max_ls);

    // Variant 3: 16-float K chunks, edges included (M and N cross a 64 tile and
    // K=19 leaves a 3-wide tail chunk).
    const int M3 = 70, N3 = 66, K3 = 19;
    std::vector<float> A3((size_t)M3 * K3), B3((size_t)N3 * K3), ref3((size_t)M3 * N3),
        got3((size_t)M3 * N3, 1e30f), sim3((size_t)M3 * N3, 1e30f);
    for (auto& v : A3) v = (float)((int)(rng() % 2001) - 1000) / 512.0f;
    for (auto& v : B3) v = (float)((int)(rng() % 2001) - 1000) / 512.0f;
    rdna2_sgemm_f32_ref(ref3.data(), A3.data(), B3.data(), M3, N3, K3);
    CHECK(rdna2_sgemm_f32_host_v(got3.data(), A3.data(), B3.data(), M3, N3, K3, 3) == 0,
          "sgemm variant 3 returned error");
    sim_sgemm_lds16_f32(sim3.data(), A3.data(), B3.data(), M3, N3, K3);
    double rel3 = 0.0, rel3s = 0.0;
    for (int i = 0; i < M3 * N3; ++i) {
      const double s = std::fabs((double)ref3[i]) > 1e-6 ? std::fabs((double)ref3[i]) : 1.0;
      rel3 = std::fmax(rel3, std::fabs((double)got3[i] - (double)ref3[i]) / s);
      rel3s = std::fmax(rel3s, std::fabs((double)sim3[i] - (double)ref3[i]) / s);
    }
    CHECK(rel3 < 1e-5, "sgemm variant 3 matches reference");
    CHECK(rel3s < 1e-5, "sgemm variant 3 simulation matches reference");
    std::printf("sgemm f16k M=%-4d N=%-4d K=%-5d max_rel=%.2e sim_rel=%.2e (64x64, KSTEP=16)\n",
                M3, N3, K3, rel3, rel3s);

    // RDNA4 variant: two 128x128 tiles per dimension, K=19 -> 8+8+3 staged.
    const int M4 = 135, N4 = 130, K4 = 19;
    std::vector<float> A4((size_t)M4 * K4), B4((size_t)N4 * K4), ref4((size_t)M4 * N4),
        got4((size_t)M4 * N4, 1e30f), sim4((size_t)M4 * N4, 1e30f);
    for (auto& v : A4) v = (float)((int)(rng() % 2001) - 1000) / 512.0f;
    for (auto& v : B4) v = (float)((int)(rng() % 2001) - 1000) / 512.0f;
    rdna2_sgemm_f32_ref(ref4.data(), A4.data(), B4.data(), M4, N4, K4);
    CHECK(rdna2_sgemm_f32_host_v(got4.data(), A4.data(), B4.data(), M4, N4, K4, 2) == 0,
          "sgemm variant 2 returned error");
    sim_sgemm_lds8x8_f32(sim4.data(), A4.data(), B4.data(), M4, N4, K4);
    double rel2 = 0.0, rel2s = 0.0;
    for (int i = 0; i < M4 * N4; ++i) {
      const double s = std::fabs((double)ref4[i]) > 1e-6 ? std::fabs((double)ref4[i]) : 1.0;
      rel2 = std::fmax(rel2, std::fabs((double)got4[i] - (double)ref4[i]) / s);
      rel2s = std::fmax(rel2s, std::fabs((double)sim4[i] - (double)ref4[i]) / s);
    }
    CHECK(rel2 < 1e-5, "sgemm variant 2 matches reference");
    CHECK(rel2s < 1e-5, "sgemm variant 2 simulation matches reference");
    std::printf("sgemm f128 M=%-4d N=%-4d K=%-5d max_rel=%.2e sim_rel=%.2e (128x128, 8x8)\n", M4, N4,
                K4, rel2, rel2s);
  }

  // --- 5. dp4a u8 dot -----------------------------------------------------
  {
    const int n = 1000;
    std::vector<unsigned> a(n), b(n);
    for (auto& v : a) v = rng();
    for (auto& v : b) v = rng();
    unsigned got = 0;
    const unsigned ref = rdna2_dot_u8_ref(a.data(), b.data(), n);
    CHECK(rdna2_dot_u8_host(&got, a.data(), b.data(), n) == 0, "dot_u8 returned error");
    CHECK(got == ref, "dot_u8 matches reference");
    std::printf("dot    u8  n=%-5d got=%u ref=%u\n", n, got, ref);
  }

  // --- 6. packed f16 dot --------------------------------------------------
  {
    const int n = 500;
    std::vector<unsigned> a(n), b(n);
    std::vector<float> av(n * 2), bv(n * 2);
    for (int i = 0; i < n * 2; ++i) {
      av[i] = (float)((int)(rng() % 4001) - 2000) / 256.0f;
      bv[i] = (float)((int)(rng() % 4001) - 2000) / 256.0f;
    }
#if defined(__FLT16_MAX__)
    for (int i = 0; i < n; ++i) {
      a[i] = rdna2_pack_half2((_Float16)av[2 * i], (_Float16)av[2 * i + 1]);
      b[i] = rdna2_pack_half2((_Float16)bv[2 * i], (_Float16)bv[2 * i + 1]);
    }
    float got = 0.0f;
    const float ref = rdna2_dot_f16_ref(a.data(), b.data(), n);
    CHECK(rdna2_dot_f16_host(&got, a.data(), b.data(), n) == 0, "dot_f16 returned error");
    double rel = std::fabs((double)got - (double)ref) / (std::fabs((double)ref) + 1e-9);
    // Atomic accumulation order is not fixed, so allow f32 rounding drift.
    CHECK(rel < 1e-4, "dot_f16 matches reference within tolerance");
    std::printf("dot    f16 n=%-5d got=%.6f ref=%.6f rel_err=%.2e\n", n, got, ref, rel);
#else
    std::printf("dot    f16 skipped (host compiler lacks _Float16)\n");
#endif
  }

  // --- 7. DPP prefix sum over one wave ------------------------------------
  {
    const int n = warp_size;
    std::vector<int32_t> in(n), ref(n), out(n);
    for (int i = 0; i < n; ++i) in[i] = (int32_t)(rng() % 2001) - 1000;
    rdna2_prefix_ref(ref.data(), in.data(), n);
    CHECK(rdna2_prefix_host(out.data(), in.data(), n) == 0, "prefix returned error");
    int bad = 0;
    for (int i = 0; i < n; ++i) bad += (out[i] != ref[i]);
    CHECK(bad == 0, "prefix matches reference");
    std::printf("prefix i32 lanes=%-4d mismatches=%d\n", n, bad);
  }

  std::printf("\n%s: %d check(s), %d failure(s), path=%s\n", failures ? "FAILED" : "PASSED",
              checks, failures, has_device ? "GPU" : "CPU fallback");
  return failures ? 1 : 0;
}
