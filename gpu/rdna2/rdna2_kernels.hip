// RDNA2 (gfx1030 / gfx1031) kernels, written directly against the GFX10.3 ISA,
// plus an RDNA4 (gfx1200 / gfx1201, RX 9060 XT / 9070 XT) tile shape that spends
// gfx12's larger register file on a 128x128 tile.
// Not a CUDA-to-HIP port: each kernel is shaped around one RDNA2 instruction
// group, and every one of them is checked against a plain-C++ reference
// (rdna2_gemm_ref.h) by test_rdna2_host.cpp.
//
// Instruction groups used:
//   dot4    : v_dot4_i32_i8 / v_dot4c_i32_i8 / v_dot4_u32_u8   dp4a, int8 x4
//   dot2    : v_dot2_f32_f16                                   packed f16 dot
//   dpp     : v_add_nc_u32_dpp row_shr                        16-lane wave scan
//   ds_*    : ds_read_b32 / ds_write_b32                         LDS staging
//
// Build: clang++ -x hip --offload-arch=gfx1031 -O2 -c rdna2_kernels.hip
//        clang++ -x hip --offload-arch=gfx1201 -O2 -c rdna2_kernels.hip  (9070 XT)
//
// Return codes shared by all launchers:
//   0 = success, 1 = HIP/runtime error, 2 = invalid argument.

#include <hip/hip_runtime.h>
#include <cstdint>
#include <cstddef>

#include "rdna2_gemm_ref.h"
#include "rdna2_kernels_math.h"

#define RDNA2_OK 0
#define RDNA2_ERR_HIP 1
#define RDNA2_ERR_ARG 2

// The RDNA2_ARCH_* target codes live in rdna2_kernels_math.h, next to the tile
// shapes they select, so the host test uses the same numbers.

// Four packed int32 = 16 bytes: one global_load_dwordx4 on the gfx1031 path.
typedef int int4_t __attribute__((ext_vector_type(4)));
// Two packed f16: the operand of v_dot2_f32_f16.
typedef _Float16 half2_t __attribute__((ext_vector_type(2)));
// WMMA operands and accumulators: V_WMMA_I32_16X16X16_IU8 takes two pairs of
// packed int32 (16 operand bytes per lane) and eight int32 accumulators.
typedef int int32x2_t __attribute__((ext_vector_type(2)));
typedef int int32x8_t __attribute__((ext_vector_type(8)));

// Whether this compilation pass can emit V_WMMA_I32_16X16X16_IU8: only the
// device passes of gfx1200 / gfx1201. Note the host pass of the same build never
// has it, so this must not decide which kernel gets launched -- the launch site
// is unconditional and only the instruction inside the kernel is arch-gated.
#if defined(__HIP_DEVICE_COMPILE__) && (defined(__gfx1200__) || defined(__gfx1201__))
#define RDNA2_WMMA_IU8 1
#else
#define RDNA2_WMMA_IU8 0
#endif


// The LDS tiling shape (RDNA2_TILE_MN, RDNA2_LDS_*_KSTEP/STRIDE) lives in
// rdna2_kernels_math.h so the host simulation uses the identical numbers.

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------
static inline int rdna2_device_count(void) {
  int n = 0;
  if (hipGetDeviceCount(&n) != hipSuccess) return 0;
  return n;
}

// 0 when the query fails or reports nonsense. Callers must then avoid shapes
// that depend on the wave width rather than guess one: the wave-per-row GEMV
// derives rows-per-block from blockDim/warpSize on the device, so a host guess
// of 32 with an actual wave64 would leave half the rows unprocessed.
static inline int rdna2_device_warp_size(void) {
  hipDeviceProp_t prop{};
  if (hipGetDeviceProperties(&prop, 0) != hipSuccess || prop.warpSize <= 0) return 0;
  const int ws = prop.warpSize;
  if (ws > 64 || (ws & (ws - 1)) != 0) return 0;  // not a power of two <= 64
  return ws;
}

static inline int rdna2_warp_size(void) {
  const int ws = rdna2_device_warp_size();
  return ws > 0 ? ws : 32;
}

// ---------------------------------------------------------------------------
// Group 1: dp4a int8 GEMV.  y[rows] = A[rows x k4] . x[k4]
// gfx1031 reads four packed words per load (16-byte vector access); any k4 is
// accepted because the remainder is handled by a scalar tail. Alignment is not
// a caller requirement: the launcher copies into a device buffer, which HIP
// allocates with far stronger alignment than the 16 bytes the vectors need.
// ---------------------------------------------------------------------------
__global__ void rdna2_dp4a_gemv_i8(int* __restrict__ y,
                                   const int* __restrict__ A,  // rows x k4 packed
                                   const int* __restrict__ x,  // k4 packed
                                   int rows, int k4) {
  int row = blockIdx.x * blockDim.x + threadIdx.x;
  if (row >= rows) return;
  const int* a_s = A + (size_t)row * k4;
  // Four independent accumulators: v_dot4 has a multi-cycle result latency, so
  // one dependent chain would leave the dot4 pipes idle. Integer addition is
  // exact, so summing the four at the end is bit-identical to one chain.
  int acc0 = 0, acc1 = 0, acc2 = 0, acc3 = 0;
#if defined(__gfx1031__) || defined(__gfx1200__) || defined(__gfx1201__)
  int j = 0;
  // The 16-byte vector load is only legal when this row also starts on a
  // 16-byte boundary. A row starts at (row * k4) packed words, so that holds
  // for every row only when k4 is a multiple of 4. For any other k4 the row
  // stride is not 16-byte aligned and v_dot4 is fed from scalar loads instead
  // (misaligned dwordx4 loads are not handled by the hardware).
  if ((k4 & 3) == 0) {
    const int4_t* a = reinterpret_cast<const int4_t*>(a_s);
    const int4_t* xv = reinterpret_cast<const int4_t*>(x);
    for (; j + 4 <= k4; j += 4) {
      int4_t av = a[j / 4];
      int4_t xw = xv[j / 4];
      acc0 = rdna2_dot4_i8_impl(av.x, xw.x, acc0);
      acc1 = rdna2_dot4_i8_impl(av.y, xw.y, acc1);
      acc2 = rdna2_dot4_i8_impl(av.z, xw.z, acc2);
      acc3 = rdna2_dot4_i8_impl(av.w, xw.w, acc3);
    }
  }
  for (; j < k4; ++j) acc0 = rdna2_dot4_i8_impl(a_s[j], x[j], acc0);
#else
  for (int j = 0; j < k4; ++j) {
    acc0 = rdna2_dot4_i8_impl(a_s[j], x[j], acc0);
  }
#endif
  y[row] = (acc0 + acc1) + (acc2 + acc3);
}

// One wave per row. Lanes stride the row, so a wave reads 32 consecutive words
// (one 128-byte line) per step instead of the 32 scattered row starts the
// thread-per-row kernel touches, and the row is read exactly once. Better than
// one thread per row once a row holds enough words to keep a wave busy; the
// launcher switches on k4 (see RDNA2_GEMV_WAVE_MIN_K4).
__global__ void rdna2_dp4a_gemv_i8_wave(int* __restrict__ y,
                                        const int* __restrict__ A,
                                        const int* __restrict__ x, int rows, int k4) {
  const int ws = (int)warpSize;
  const int lane = (int)(threadIdx.x % warpSize);
  const int wave = (int)(threadIdx.x / warpSize);
  const int waves = (int)(blockDim.x / warpSize);
  const int row = (int)blockIdx.x * waves + wave;
  if (row >= rows) return;  // wave-uniform: a wave never splits across rows
  const int* a_s = A + (size_t)row * k4;
  int acc = 0;
  for (int k = lane; k < k4; k += ws) acc = rdna2_dot4_i8_impl(a_s[k], x[k], acc);
  for (int off = ws >> 1; off > 0; off >>= 1) acc += __shfl_xor(acc, off);
  if (lane == 0) y[row] = acc;
}

// ---------------------------------------------------------------------------
// Group 1c: f32 GEMV, y[M] = A[M x K] . x[K]. The same two shapes as the int8
// GEMV above: one thread per row (short rows) and one wave per row (long rows,
// lanes stride the row so a wave touches consecutive floats).
// ---------------------------------------------------------------------------
__global__ void rdna2_sgemv_f32(float* __restrict__ y, const float* __restrict__ A,
                                const float* __restrict__ x, int rows, int K) {
  const int row = (int)blockIdx.x * blockDim.x + (int)threadIdx.x;
  if (row >= rows) return;
  const float* a_s = A + (size_t)row * K;
  y[row] = rdna2_sgemv_row_f32(a_s, x, K);
}

__global__ void rdna2_sgemv_f32_wave(float* __restrict__ y, const float* __restrict__ A,
                                     const float* __restrict__ x, int rows, int K) {
  const int ws = (int)warpSize;
  const int lane = (int)(threadIdx.x % warpSize);
  const int wave = (int)(threadIdx.x / warpSize);
  const int waves = (int)(blockDim.x / warpSize);
  const int row = (int)blockIdx.x * waves + wave;
  if (row >= rows) return;  // wave-uniform: a wave never splits across rows
  const float* a_s = A + (size_t)row * K;
  float acc = rdna2_sgemv_lane_partial_f32(a_s, x, K, lane, ws);
  for (int off = ws >> 1; off > 0; off >>= 1) acc += __shfl_xor(acc, off);
  if (lane == 0) y[row] = acc;
}

// ---------------------------------------------------------------------------
// Group 1b: dp4a int8 GEMM.
//   C[M x N] = A[M x K4] * Bt[N x K4]^T   (B stored transposed, N rows of K4)
// Two implementations: a naive one (one thread per element) and a 4x4
// register-tiled one. The device test compares both against the CPU reference.
// ---------------------------------------------------------------------------
__global__ void rdna2_dp4a_gemm_i8_naive(int* __restrict__ C,
                                         const int* __restrict__ A,
                                         const int* __restrict__ Bt,
                                         int M, int N, int K4) {
  int n = blockIdx.x * blockDim.x + threadIdx.x;
  int m = blockIdx.y;
  if (m >= M || n >= N) return;
  const int* a = A + (size_t)m * K4;
  const int* b = Bt + (size_t)n * K4;
  C[(size_t)m * N + n] = rdna2_gemm_cell_i8(a, b, K4);
}

// 4x4 outputs per thread; block (16,16) covers a 64x64 output tile.
__global__ void rdna2_dp4a_gemm_i8_tile(int* __restrict__ C,
                                        const int* __restrict__ A,
                                        const int* __restrict__ Bt,
                                        int M, int N, int K4) {
  const int m0 = (blockIdx.y * blockDim.y + threadIdx.y) * 4;
  const int n0 = (blockIdx.x * blockDim.x + threadIdx.x) * 4;
  if (m0 >= M || n0 >= N) return;

  int32_t acc[4][4];
  rdna2_gemm_tile_i8(acc, A, Bt, M, N, K4, m0, n0);
  for (int i = 0; i < 4; ++i) {
    if (m0 + i >= M) break;
    for (int j = 0; j < 4; ++j) {
      if (n0 + j >= N) break;
      C[(size_t)(m0 + i) * N + (n0 + j)] = acc[i][j];
    }
  }
}

// ---------------------------------------------------------------------------
// Group 1c: LDS-staged tiled int8 GEMM, same 64x64 output tile and 4x4 register
// blocking as above, but the K loop reads the operands from LDS instead of
// global memory. The thread-per-cell version issues 8 packed loads per k step
// per thread; over a 64x64 tile that is 2048 words of traffic for 128 words of
// actual data (16x). Staging each K chunk cooperatively costs 1024 loads spread
// over 256 threads (4 each) and then every reuse is LDS.
// ---------------------------------------------------------------------------
// Stage one K chunk of a matrix into LDS. 256 threads cover 64 rows x 16 words;
// each thread owns 4 consecutive words of one row. Rows past the edge and words
// past K4 are written as zero, so the compute loop needs no bounds checks.
// lds is a template parameter so the __shared__ address space survives: a fixed
// int (*)[STRIDE] parameter would genericize it to flat stores.
template <typename LDP>
__device__ static inline void rdna2_stage_i8(LDP lds, const int* __restrict__ src, int nrows,
                                            int ncols, int row0, int k0, int ld_row, int ld_col) {
  const int gm = row0 + ld_row;
  int4_t v;
  v.x = v.y = v.z = v.w = 0;
  if (gm < nrows) {
    const int* p = src + (size_t)gm * ncols + k0 + ld_col;
    if ((ncols & 3) == 0 && k0 + ld_col + 4 <= ncols) {
      v = *reinterpret_cast<const int4_t*>(p);  // 16-byte access, aligned by ncols%4==0
    } else {
      for (int c = 0; c < 4; ++c)
        if (k0 + ld_col + c < ncols) v[c] = p[c];
    }
  }
  lds[ld_row][ld_col + 0] = v.x;
  lds[ld_row][ld_col + 1] = v.y;
  lds[ld_row][ld_col + 2] = v.z;
  lds[ld_row][ld_col + 3] = v.w;
}

__global__ void rdna2_dp4a_gemm_i8_lds(int* __restrict__ C,
                                       const int* __restrict__ A,
                                       const int* __restrict__ Bt,
                                       int M, int N, int K4) {
  __shared__ int a_lds[RDNA2_TILE_MN][RDNA2_LDS_I8_STRIDE];
  __shared__ int b_lds[RDNA2_TILE_MN][RDNA2_LDS_I8_STRIDE];
  const int tid = threadIdx.y * blockDim.x + threadIdx.x;
  const int m_blk = blockIdx.y * RDNA2_TILE_MN;
  const int n_blk = blockIdx.x * RDNA2_TILE_MN;
  const int lm = threadIdx.y * 4;  // this thread's first tile row
  const int ln = threadIdx.x * 4;  // this thread's first tile column
  const int ld_row = tid >> 2;
  const int ld_col = (tid & 3) << 2;

  int32_t acc[4][4];
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) acc[i][j] = 0;

  for (int k0 = 0; k0 < K4; k0 += RDNA2_LDS_I8_KSTEP) {
    rdna2_stage_i8(a_lds, A, M, K4, m_blk, k0, ld_row, ld_col);
    rdna2_stage_i8(b_lds, Bt, N, K4, n_blk, k0, ld_row, ld_col);
    __syncthreads();
    int kend = K4 - k0;
    if (kend > RDNA2_LDS_I8_KSTEP) kend = RDNA2_LDS_I8_KSTEP;
    rdna2_gemm_lds_block_i8(acc, &a_lds[0][0], &b_lds[0][0], RDNA2_LDS_I8_STRIDE, lm, ln, kend);
    __syncthreads();  // nobody may restage a tile that is still being read
  }

  for (int i = 0; i < 4; ++i) {
    if (m_blk + lm + i >= M) break;
    for (int j = 0; j < 4; ++j) {
      if (n_blk + ln + j >= N) break;
      C[(size_t)(m_blk + lm + i) * N + (n_blk + ln + j)] = acc[i][j];
    }
  }
}

// Same LDS tile and K loop, 8x4 register blocking: 128 threads (16,8) cover the
// same 64x64 output tile, each owning 8 rows x 4 columns. Twelve LDS reads per k
// step feed 32 dp4a instead of eight feeding 16, i.e. a quarter less LDS traffic
// for identical math. Staging spreads the same 64x16 words over 128 threads, so
// each thread writes 8 words (two 16-byte loads) instead of 4.
__global__ void rdna2_dp4a_gemm_i8_lds8(int* __restrict__ C,
                                        const int* __restrict__ A,
                                        const int* __restrict__ Bt,
                                        int M, int N, int K4) {
  __shared__ int a_lds[RDNA2_TILE_MN][RDNA2_LDS_I8_STRIDE];
  __shared__ int b_lds[RDNA2_TILE_MN][RDNA2_LDS_I8_STRIDE];
  const int tid = threadIdx.y * blockDim.x + threadIdx.x;
  const int m_blk = blockIdx.y * RDNA2_TILE_MN;
  const int n_blk = blockIdx.x * RDNA2_TILE_MN;
  const int lm = threadIdx.y * 8;  // this thread's first tile row
  const int ln = threadIdx.x * 4;  // this thread's first tile column
  const int ld_row = tid >> 1;     // 0..63, two threads per staged row
  const int ld_col = (tid & 1) << 3;

  int32_t acc[8][4];
  for (int i = 0; i < 8; ++i)
    for (int j = 0; j < 4; ++j) acc[i][j] = 0;

  for (int k0 = 0; k0 < K4; k0 += RDNA2_LDS_I8_KSTEP) {
    rdna2_stage_i8(a_lds, A, M, K4, m_blk, k0, ld_row, ld_col);
    rdna2_stage_i8(a_lds, A, M, K4, m_blk, k0, ld_row, ld_col + 4);
    rdna2_stage_i8(b_lds, Bt, N, K4, n_blk, k0, ld_row, ld_col);
    rdna2_stage_i8(b_lds, Bt, N, K4, n_blk, k0, ld_row, ld_col + 4);
    __syncthreads();
    int kend = K4 - k0;
    if (kend > RDNA2_LDS_I8_KSTEP) kend = RDNA2_LDS_I8_KSTEP;
    rdna2_gemm_lds_block8_i8(acc, &a_lds[0][0], &b_lds[0][0], RDNA2_LDS_I8_STRIDE, lm, ln, kend);
    __syncthreads();  // nobody may restage a tile that is still being read
  }

  for (int i = 0; i < 8; ++i) {
    if (m_blk + lm + i >= M) break;
    for (int j = 0; j < 4; ++j) {
      if (n_blk + ln + j >= N) break;
      C[(size_t)(m_blk + lm + i) * N + (n_blk + ln + j)] = acc[i][j];
    }
  }
}

// ---------------------------------------------------------------------------
// Group 1e: RDNA4 (gfx1200 / gfx1201) path. Same 16x16 block of 256 threads and
// the same 16-word K chunk as groups 1c/1d, but a 128x128 output tile with 8x8
// outputs per thread: sixteen LDS reads per k step feed 64 dp4a instead of eight
// feeding 16, and every staged word pays for four times as many results. gfx12
// has 1536 VGPRs per SIMD against 1024 on gfx10/11, so 64 accumulators plus 16
// staged operands still leave several waves resident. Two staged tiles of
// 128 x 17 words = 17408 bytes, inside the 64 KiB a gfx12 workgroup addresses.
// Select with tiled = 4 (see rdna2_tiled_for_arch).
// ---------------------------------------------------------------------------
__global__ void rdna2_dp4a_gemm_i8_lds128(int* __restrict__ C,
                                          const int* __restrict__ A,
                                          const int* __restrict__ Bt,
                                          int M, int N, int K4) {
  __shared__ int a_lds[RDNA2_T4_TILE_MN][RDNA2_T4_LDS_I8_STRIDE];
  __shared__ int b_lds[RDNA2_T4_TILE_MN][RDNA2_T4_LDS_I8_STRIDE];
  const int tid = threadIdx.y * blockDim.x + threadIdx.x;
  const int m_blk = blockIdx.y * RDNA2_T4_TILE_MN;
  const int n_blk = blockIdx.x * RDNA2_T4_TILE_MN;
  const int lm = threadIdx.y * 8;     // this thread's first tile row
  const int ln = threadIdx.x * 8;     // this thread's first tile column
  const int ld_row = tid >> 1;        // 0..127: two threads per staged row
  const int ld_col = (tid & 1) << 3;  // 0 or 8, eight staged words each

  int32_t acc[8][8];
  for (int i = 0; i < 8; ++i)
    for (int j = 0; j < 8; ++j) acc[i][j] = 0;

  for (int k0 = 0; k0 < K4; k0 += RDNA2_T4_LDS_I8_KSTEP) {
    rdna2_stage_i8(a_lds, A, M, K4, m_blk, k0, ld_row, ld_col);
    rdna2_stage_i8(a_lds, A, M, K4, m_blk, k0, ld_row, ld_col + 4);
    rdna2_stage_i8(b_lds, Bt, N, K4, n_blk, k0, ld_row, ld_col);
    rdna2_stage_i8(b_lds, Bt, N, K4, n_blk, k0, ld_row, ld_col + 4);
    __syncthreads();
    int kend = K4 - k0;
    if (kend > RDNA2_T4_LDS_I8_KSTEP) kend = RDNA2_T4_LDS_I8_KSTEP;
    rdna2_gemm_lds_block8x8_i8(acc, &a_lds[0][0], &b_lds[0][0], RDNA2_T4_LDS_I8_STRIDE, lm,
                               ln, kend);
    __syncthreads();  // nobody may restage a tile that is still being read
  }

  for (int i = 0; i < 8; ++i) {
    if (m_blk + lm + i >= M) break;
    for (int j = 0; j < 8; ++j) {
      if (n_blk + ln + j >= N) break;
      C[(size_t)(m_blk + lm + i) * N + (n_blk + ln + j)] = acc[i][j];
    }
  }
}

// ---------------------------------------------------------------------------
// Group 1f: RDNA4 (gfx1200 / gfx1201) matrix path. Same 128x128 output tile,
// same 16-word K chunk and same LDS staging as group 1e, but the arithmetic is
// the gfx12 matrix instruction V_WMMA_I32_16X16X16_IU8 instead of 64 dp4a per
// thread: one instruction covers a 16x16x16 signed product that the 8x8 dp4a
// kernel needs 64 instructions and 16 LDS reads per k step to reach.
//
// The eight waves of the block each own a 16-row band of the tile and sweep the
// eight 16-column groups, keeping 8 x int32x8 accumulators (64 VGPRs, the same
// budget as group 1e's 8x8 block). Per K chunk a wave issues 32 WMMA -- 8 per
// 16-byte K step -- and pays only 2 LDS words for A and 2 per B fragment, i.e.
// 36 word reads for 32 matrix ops instead of 256 reads for 1024 dp4a.
//
// Fragment layout, measured on a gfx1201 card (256/256 outputs of a signed
// 16x16x16 product matched the scalar reference; the layout and the negA/negB
// semantics were both derived experimentally, see README):
//   A: lane t, byte e -> A[M = t%16][K = (t/16)*8 + e]
//   B: lane t, byte e -> B[K = (t/16)*8 + e][N = t%16]
//   D: lane t, slot l -> D[M = (t/16)*8 + l][N = t%16]
// Both operands are signed, so negA and negB are true; __builtin_amdgcn_wmma_*
// needs those flags as compile-time constants, hence the literals.
// Select with tiled = 5 (see rdna2_tiled_for_arch).
// ---------------------------------------------------------------------------
// One 16x16x16 signed step for this lane: acc holds this lane's eight output
// slots, a and b are its eight-byte A and B fragments.
//   gfx12          -- the matrix instruction.
//   other AMD GPU  -- no matrix unit on gfx10/11, so the same result is built
//                     from explicit byte dots, gathering each slot's operands
//                     out of lanes row and row+16 with shuffles. Correct, but
//                     with a shuffle per operand byte it cannot beat the dp4a
//                     kernel, and rdna2_tiled_for_arch never picks tiled = 5 on
//                     those chips. Not measured: no RDNA2 hardware here.
//   plain C++ pass -- the file also compiles as host code (CPU test harness and
//                     the host half of a HIP build), where no kernel can run; the
//                     body is only there to keep the kernel defined so the
//                     unconditional launch site links.
__device__ static inline int32x8_t rdna2_wmma_iu8_step(int32x2_t a, int32x2_t b, int32x8_t acc,
                                                       int lane) {
  (void)lane;  // the matrix instruction addresses the other lanes in hardware
#if RDNA2_WMMA_IU8
  return __builtin_amdgcn_wmma_i32_16x16x16_iu8_w32_gfx12(true, a, true, b, acc, false);
#elif defined(__HIP_DEVICE_COMPILE__)
  for (int l = 0; l < 8; ++l) {
    const int dm = rdna2_wmma_i8_acc_row(lane, l);
    const int dn = rdna2_wmma_i8_acc_col(lane);
    int32_t d = 0;
    for (int k = 0; k < 16; ++k) {
      const int shift = 8 * (k & 3), word = (k & 7) >> 2;
      const int av = (int8_t)(__shfl((int)a[word], dm + (k >> 3) * 16) >> shift);
      const int bv = (int8_t)(__shfl((int)b[word], dn + (k >> 3) * 16) >> shift);
      d += av * bv;
    }
    acc[l] += d;
  }
  return acc;
#else
  (void)a;
  (void)b;
  return acc;
#endif
}

__global__ void rdna2_dp4a_gemm_i8_wmma(int* __restrict__ C,
                                        const int* __restrict__ A,
                                        const int* __restrict__ Bt,
                                        int M, int N, int K4) {
  __shared__ int a_lds[RDNA2_T4_TILE_MN][RDNA2_T4_LDS_I8_STRIDE];
  __shared__ int b_lds[RDNA2_T4_TILE_MN][RDNA2_T4_LDS_I8_STRIDE];
  const int tid = threadIdx.y * blockDim.x + threadIdx.x;
  const int lane = tid & 31;
  const int wave = tid >> 5;  // 0..7: this wave's 16-row band of the tile
  const int m_blk = blockIdx.y * RDNA2_T4_TILE_MN;
  const int n_blk = blockIdx.x * RDNA2_T4_TILE_MN;
  const int m_row = wave * 16 + rdna2_wmma_i8_row(lane);  // staged A row this lane feeds
  const int n_col = rdna2_wmma_i8_row(lane);              // column group offset inside a 16-wide group
  const int ld_row = tid >> 1;                            // 0..127: two threads per staged row
  const int ld_col = (tid & 1) << 3;                      // 0 or 8, eight staged words each

  int32x8_t acc[RDNA2_T4_WMMA_NTILE];
  for (int j = 0; j < RDNA2_T4_WMMA_NTILE; ++j) {
    for (int l = 0; l < 8; ++l) acc[j][l] = 0;
  }

  for (int k0 = 0; k0 < K4; k0 += RDNA2_T4_LDS_I8_KSTEP) {
    rdna2_stage_i8(a_lds, A, M, K4, m_blk, k0, ld_row, ld_col);
    rdna2_stage_i8(a_lds, A, M, K4, m_blk, k0, ld_row, ld_col + 4);
    rdna2_stage_i8(b_lds, Bt, N, K4, n_blk, k0, ld_row, ld_col);
    rdna2_stage_i8(b_lds, Bt, N, K4, n_blk, k0, ld_row, ld_col + 4);
    __syncthreads();
    int kend = K4 - k0;
    if (kend > RDNA2_T4_LDS_I8_KSTEP) kend = RDNA2_T4_LDS_I8_KSTEP;
    // One WMMA covers 4 staged words (16 K bytes); a short last chunk runs the
    // steps it has and relies on the staging zeros for the rest of the tile.
    const int ksteps = (kend + 3) / 4;
    for (int ks = 0; ks < ksteps; ++ks) {
      const int w = rdna2_wmma_i8_frag_word(lane, ks);
      int32x2_t a;
      a[0] = a_lds[m_row][w];
      a[1] = a_lds[m_row][w + 1];
      for (int j = 0; j < RDNA2_T4_WMMA_NTILE; ++j) {
        int32x2_t b;
        b[0] = b_lds[j * 16 + n_col][w];
        b[1] = b_lds[j * 16 + n_col][w + 1];
        acc[j] = rdna2_wmma_iu8_step(a, b, acc[j], lane);
      }
    }
    __syncthreads();  // nobody may restage a tile that is still being read
  }

  // D map: slot l of lane t lands at row 8*(t/16) + l, column t%16 of the
  // 16x16 tile, so a lane writes one column of its wave's 16-row band.
  for (int j = 0; j < RDNA2_T4_WMMA_NTILE; ++j) {
    const int n = n_blk + j * 16 + rdna2_wmma_i8_acc_col(lane);
    if (n >= N) continue;
    for (int l = 0; l < 8; ++l) {
      const int m = m_blk + wave * 16 + rdna2_wmma_i8_acc_row(lane, l);
      if (m >= M) break;  // rows 8*(t/16)+l are contiguous, so a break is safe
      C[(size_t)m * N + n] = acc[j][l];
    }
  }
}

// ---------------------------------------------------------------------------
// Group 2: f32 GEMM, 4x4 register tiled.  C[M x N] = A[M x K] * Bt[N x K]^T
// ---------------------------------------------------------------------------
__global__ void rdna2_sgemm_f32(float* __restrict__ C,
                                const float* __restrict__ A,
                                const float* __restrict__ Bt,
                                int M, int N, int K) {
  const int m0 = (blockIdx.y * blockDim.y + threadIdx.y) * 4;
  const int n0 = (blockIdx.x * blockDim.x + threadIdx.x) * 4;
  if (m0 >= M || n0 >= N) return;

  float acc[4][4];
  rdna2_sgemm_tile_f32(acc, A, Bt, M, N, K, m0, n0);
  for (int i = 0; i < 4; ++i) {
    if (m0 + i >= M) break;
    for (int j = 0; j < 4; ++j) {
      if (n0 + j >= N) break;
      C[(size_t)(m0 + i) * N + (n0 + j)] = acc[i][j];
    }
  }
}

// ---------------------------------------------------------------------------
// Group 2b: LDS-staged f32 GEMM, same 64x64 tile and 4x4 register blocking as
// group 1c: the K loop reads A and Bt from LDS instead of re-reading global
// memory in every thread.
// ---------------------------------------------------------------------------
// 256 threads cover 64 rows x 8 floats: two consecutive floats per thread.
template <typename LDP>
__device__ static inline void rdna2_stage_f32(LDP lds, const float* __restrict__ src, int nrows,
                                              int ncols, int row0, int k0, int ld_row, int ld_col) {
  const int gm = row0 + ld_row;
  float v0 = 0.0f, v1 = 0.0f;
  if (gm < nrows) {
    const float* p = src + (size_t)gm * ncols + k0 + ld_col;
    if (k0 + ld_col + 0 < ncols) v0 = p[0];
    if (k0 + ld_col + 1 < ncols) v1 = p[1];
  }
  lds[ld_row][ld_col + 0] = v0;
  lds[ld_row][ld_col + 1] = v1;
}

__global__ void rdna2_sgemm_f32_lds(float* __restrict__ C,
                                    const float* __restrict__ A,
                                    const float* __restrict__ Bt,
                                    int M, int N, int K) {
  __shared__ float a_lds[RDNA2_TILE_MN][RDNA2_LDS_F32_STRIDE];
  __shared__ float b_lds[RDNA2_TILE_MN][RDNA2_LDS_F32_STRIDE];
  const int tid = threadIdx.y * blockDim.x + threadIdx.x;
  const int m_blk = blockIdx.y * RDNA2_TILE_MN;
  const int n_blk = blockIdx.x * RDNA2_TILE_MN;
  const int lm = threadIdx.y * 4;
  const int ln = threadIdx.x * 4;
  const int ld_row = tid >> 2;
  const int ld_col = (tid & 3) << 1;

  float acc[4][4];
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) acc[i][j] = 0.0f;

  for (int k0 = 0; k0 < K; k0 += RDNA2_LDS_F32_KSTEP) {
    rdna2_stage_f32(a_lds, A, M, K, m_blk, k0, ld_row, ld_col);
    rdna2_stage_f32(b_lds, Bt, N, K, n_blk, k0, ld_row, ld_col);
    __syncthreads();
    int kend = K - k0;
    if (kend > RDNA2_LDS_F32_KSTEP) kend = RDNA2_LDS_F32_KSTEP;
    rdna2_sgemm_lds_block_f32(acc, &a_lds[0][0], &b_lds[0][0], RDNA2_LDS_F32_STRIDE, lm, ln, kend);
    __syncthreads();
  }

  for (int i = 0; i < 4; ++i) {
    if (m_blk + lm + i >= M) break;
    for (int j = 0; j < 4; ++j) {
      if (n_blk + ln + j >= N) break;
      C[(size_t)(m_blk + lm + i) * N + (n_blk + ln + j)] = acc[i][j];
    }
  }
}

// ---------------------------------------------------------------------------
// Group 2b: 64x64 f32 tile with a 16-float K chunk (variant 3). Same block
// shape and 4x4 per-thread blocking as group 2, but 256 threads stage 1024
// floats (4 each, RDNA2_LDS_F32B_STRIDE = 17 per row) instead of 512 floats
// (2 each) per iteration. That makes each staged chunk cover a whole 64-byte
// cache line of every row, and halves the number of __syncthreads over K.
// ---------------------------------------------------------------------------
template <typename LDP>
__device__ static inline void rdna2_stage_f32_4(LDP lds, const float* __restrict__ src, int nrows,
                                                int ncols, int row0, int k0, int ld_row,
                                                int ld_col) {
  const int gm = row0 + ld_row;
  float v0 = 0.0f, v1 = 0.0f, v2 = 0.0f, v3 = 0.0f;
  if (gm < nrows) {
    const float* p = src + (size_t)gm * ncols + k0 + ld_col;
    const int rem = ncols - (k0 + ld_col);
    if (rem > 0) v0 = p[0];
    if (rem > 1) v1 = p[1];
    if (rem > 2) v2 = p[2];
    if (rem > 3) v3 = p[3];
  }
  lds[ld_row][ld_col + 0] = v0;
  lds[ld_row][ld_col + 1] = v1;
  lds[ld_row][ld_col + 2] = v2;
  lds[ld_row][ld_col + 3] = v3;
}

__global__ void rdna2_sgemm_f32_lds16(float* __restrict__ C,
                                      const float* __restrict__ A,
                                      const float* __restrict__ Bt,
                                      int M, int N, int K) {
  __shared__ float a_lds[RDNA2_TILE_MN][RDNA2_LDS_F32B_STRIDE];
  __shared__ float b_lds[RDNA2_TILE_MN][RDNA2_LDS_F32B_STRIDE];
  const int tid = threadIdx.y * blockDim.x + threadIdx.x;
  const int m_blk = blockIdx.y * RDNA2_TILE_MN;
  const int n_blk = blockIdx.x * RDNA2_TILE_MN;
  const int lm = threadIdx.y * 4;
  const int ln = threadIdx.x * 4;
  const int ld_row = tid >> 2;
  const int ld_col = (tid & 3) << 2;

  float acc[4][4];
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) acc[i][j] = 0.0f;

  for (int k0 = 0; k0 < K; k0 += RDNA2_LDS_F32B_KSTEP) {
    rdna2_stage_f32_4(a_lds, A, M, K, m_blk, k0, ld_row, ld_col);
    rdna2_stage_f32_4(b_lds, Bt, N, K, n_blk, k0, ld_row, ld_col);
    __syncthreads();
    int kend = K - k0;
    if (kend > RDNA2_LDS_F32B_KSTEP) kend = RDNA2_LDS_F32B_KSTEP;
    rdna2_sgemm_lds_block_f32(acc, &a_lds[0][0], &b_lds[0][0], RDNA2_LDS_F32B_STRIDE, lm, ln, kend);
    __syncthreads();
  }

  for (int i = 0; i < 4; ++i) {
    if (m_blk + lm + i >= M) break;
    for (int j = 0; j < 4; ++j) {
      if (n_blk + ln + j >= N) break;
      C[(size_t)(m_blk + lm + i) * N + (n_blk + ln + j)] = acc[i][j];
    }
  }
}

// ---------------------------------------------------------------------------
// Group 2c: RDNA4 (gfx1200 / gfx1201) f32 GEMM on the same 128x128 tile and 8x8
// per-thread blocking as group 1e. Two staged tiles of 128 x 9 floats = 9216
// bytes; sixteen LDS reads per k step feed 64 fma. Select with variant 2.
// ---------------------------------------------------------------------------
__global__ void rdna2_sgemm_f32_lds128(float* __restrict__ C,
                                       const float* __restrict__ A,
                                       const float* __restrict__ Bt,
                                       int M, int N, int K) {
  __shared__ float a_lds[RDNA2_T4_TILE_MN][RDNA2_T4_LDS_F32_STRIDE];
  __shared__ float b_lds[RDNA2_T4_TILE_MN][RDNA2_T4_LDS_F32_STRIDE];
  const int tid = threadIdx.y * blockDim.x + threadIdx.x;
  const int m_blk = blockIdx.y * RDNA2_T4_TILE_MN;
  const int n_blk = blockIdx.x * RDNA2_T4_TILE_MN;
  const int lm = threadIdx.y * 8;
  const int ln = threadIdx.x * 8;
  const int ld_row = tid >> 1;        // 0..127: two threads per staged row
  const int ld_col = (tid & 1) << 2;  // 0 or 4, four staged floats each

  float acc[8][8];
  for (int i = 0; i < 8; ++i)
    for (int j = 0; j < 8; ++j) acc[i][j] = 0.0f;

  for (int k0 = 0; k0 < K; k0 += RDNA2_T4_LDS_F32_KSTEP) {
    rdna2_stage_f32(a_lds, A, M, K, m_blk, k0, ld_row, ld_col);
    rdna2_stage_f32(a_lds, A, M, K, m_blk, k0, ld_row, ld_col + 2);
    rdna2_stage_f32(b_lds, Bt, N, K, n_blk, k0, ld_row, ld_col);
    rdna2_stage_f32(b_lds, Bt, N, K, n_blk, k0, ld_row, ld_col + 2);
    __syncthreads();
    int kend = K - k0;
    if (kend > RDNA2_T4_LDS_F32_KSTEP) kend = RDNA2_T4_LDS_F32_KSTEP;
    rdna2_sgemm_lds_block8x8_f32(acc, &a_lds[0][0], &b_lds[0][0], RDNA2_T4_LDS_F32_STRIDE, lm,
                                 ln, kend);
    __syncthreads();
  }

  for (int i = 0; i < 8; ++i) {
    if (m_blk + lm + i >= M) break;
    for (int j = 0; j < 8; ++j) {
      if (n_blk + ln + j >= N) break;
      C[(size_t)(m_blk + lm + i) * N + (n_blk + ln + j)] = acc[i][j];
    }
  }
}

// ---------------------------------------------------------------------------
// Group 3: packed f16 dot (v_dot2_f32_f16; clang emits the accumulate form
// from __builtin_amdgcn_fdot2, and the ISA mnemonic has no 'c').
// out is accumulated with atomics and must be zeroed by the launcher.
// Correct for both wave32 and wave64: the shuffle reduction covers one
// 16-lane row, then one atomic per row.
// ---------------------------------------------------------------------------
__global__ void rdna2_dot2_f16(float* __restrict__ out,
                               const half2_t* __restrict__ a,
                               const half2_t* __restrict__ b,
                               int n) {
  float acc = 0.0f;
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n;
       i += gridDim.x * blockDim.x) {
    acc = __builtin_amdgcn_fdot2(a[i], b[i], acc, false);
  }
  acc += __shfl_xor(acc, 1);
  acc += __shfl_xor(acc, 2);
  acc += __shfl_xor(acc, 4);
  acc += __shfl_xor(acc, 8);
  if ((threadIdx.x & 15) == 0) atomicAdd(out, acc);
}

// ---------------------------------------------------------------------------
// Group 4: dp4a u8 dot (v_dot4_u32_u8).
// out is accumulated with atomics and must be zeroed by the launcher.
// ---------------------------------------------------------------------------
__global__ void rdna2_dp4a_dot_u8(unsigned* __restrict__ out,
                                  const unsigned* __restrict__ a,
                                  const unsigned* __restrict__ b,
                                  int n) {
  unsigned acc = 0;
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n;
       i += gridDim.x * blockDim.x) {
    acc = __builtin_amdgcn_udot4(a[i], b[i], acc, false);
  }
  acc += __shfl_xor(acc, 1);
  acc += __shfl_xor(acc, 2);
  acc += __shfl_xor(acc, 4);
  acc += __shfl_xor(acc, 8);
  if ((threadIdx.x & 15) == 0) atomicAdd(out, acc);
}

// ---------------------------------------------------------------------------
// Group 5: inclusive prefix sum over one wave.
// DPP row_shr builds a 16-lane Hillis-Steele scan with no LDS traffic; the
// cross-half total then arrives through a shuffle (v_readlane/v_permlane).
// out[tid] = in[0] + ... + in[tid] for every lane tid < n <= warpSize.
// ---------------------------------------------------------------------------
__global__ void rdna2_dpp_prefix(int* __restrict__ out,
                                 const int* __restrict__ in,
                                 int n) {
  const int lane = threadIdx.x & (warpSize - 1);
  int v = (lane < n) ? in[lane] : 0;

  // row_shr:1,2,4,8 with bound_ctrl off: lanes shifted out of the row get 0.
  v += __builtin_amdgcn_update_dpp(0, v, 0x111, 0xf, 0xf, false);
  v += __builtin_amdgcn_update_dpp(0, v, 0x112, 0xf, 0xf, false);
  v += __builtin_amdgcn_update_dpp(0, v, 0x114, 0xf, 0xf, false);
  v += __builtin_amdgcn_update_dpp(0, v, 0x118, 0xf, 0xf, false);

  // Lane 15 holds the total of lanes 0-15; add it to the upper half.
  int lo = __shfl(v, 15);
  if (lane >= 16) v += lo;
  // On wave64, lanes 32-63 also need the total of lanes 0-31.
  if (warpSize == 64) {
    int mid = __shfl(v, 31);
    if (lane >= 32) v += mid;
  }
  out[lane] = v;
}

// ---------------------------------------------------------------------------
// Group 6: quantized epilogue (scale and zero point)
// ---------------------------------------------------------------------------
// Per-row sum of the int8 activations. dp4a against an all-ones operand sums the
// four signed bytes of a packed word in one instruction, so this is one v_dot4
// per word. Only needed for asymmetric weights (a nonzero zero point); a
// symmetric format passes zp_b = 0 and skips it entirely.
__global__ void rdna2_i8_row_sums(int* __restrict__ sums, const int* __restrict__ A, int rows,
                                  int k4) {
  const int row = blockIdx.x * blockDim.x + threadIdx.x;
  if (row >= rows) return;
  const int* a = A + (size_t)row * k4;
  const int ones = 0x01010101;
  int acc = 0;
  for (int k = 0; k < k4; ++k) acc = rdna2_dot4_i8_impl(a[k], ones, acc);
  sums[row] = acc;
}

// Applies the quantized format's scale and zero point to the finished integer
// accumulators:
//   C[m][n] = (acc[m][n] - zp_b * sums[m]) * sa[m] * sb[n]
// The correction is the exact integer form of B's asymmetric storage: B's stored
// value at k represents (q - zp_b), so sum_k a*(q - zp_b) = acc - zp_b*sum_k a.
// No stored weight is converted here - only the completed dot product is, once,
// after the correction has been applied in 64-bit so it cannot overflow. sa is
// NULL for unit row scale, sb NULL for unit column scale.
__global__ void rdna2_gemm_i8_epilogue(float* __restrict__ C, const int* __restrict__ acc,
                                       const int* __restrict__ sums, int M, int N,
                                       const float* __restrict__ sa, const float* __restrict__ sb,
                                       int zp_b) {
  const int n = blockIdx.x * blockDim.x + threadIdx.x;
  if (n >= N) return;
  // Rows stride over gridDim.y so the row count is not capped at the grid's
  // y limit; each iteration is an independent output element.
  for (int m = blockIdx.y; m < M; m += gridDim.y) {
    long long v = (long long)acc[(size_t)m * N + n];
    if (zp_b != 0 && sums) v -= (long long)zp_b * (long long)sums[m];
    float f = (float)v;
    if (sa) f *= sa[m];
    if (sb) f *= sb[n];
    C[(size_t)m * N + n] = f;
  }
}

// ---------------------------------------------------------------------------
// Host launchers
// ---------------------------------------------------------------------------
// One thread per row is the better shape for short rows (no reduction to pay
// for); from this many packed words per row on, a whole wave per row wins on
// coalescing.
#define RDNA2_GEMV_WAVE_MIN_K4 8
// f32 rows need a real body before a wave pays for its reduction. Measured on a
// 9070 XT at 2048-16384 rows: below K=256 the five shuffle steps cost more than
// the coalescing gains (thread-per-row is 1.2-3.8x ahead at K=64..160), and from
// K=256 up the wave-per-row kernel wins every time (1.15x at K=256, ~2x by
// K=1024, 3.5x at K=4096 with 16384 rows).
#define RDNA2_SGEMV_WAVE_MIN_K 256

static int rdna2_gemv_i8_device(int* y, const int* A, const int* x, int rows, int k4,
                                int variant) {
  int *dA = nullptr, *dx = nullptr, *dy = nullptr;
  const size_t bytesA = (size_t)rows * k4 * sizeof(int);
  const size_t bytesx = (size_t)k4 * sizeof(int);
  const size_t bytesy = (size_t)rows * sizeof(int);
  int rc = RDNA2_OK;
  int ws = 0;  // wave width; stays 0 unless the wave-per-row variant can size its grid
  if (hipMalloc(&dA, bytesA) != hipSuccess || hipMalloc(&dx, bytesx) != hipSuccess ||
      hipMalloc(&dy, bytesy) != hipSuccess) {
    rc = RDNA2_ERR_HIP;
    goto cleanup;
  }
  if (hipMemcpy(dA, A, bytesA, hipMemcpyHostToDevice) != hipSuccess ||
      hipMemcpy(dx, x, bytesx, hipMemcpyHostToDevice) != hipSuccess) {
    rc = RDNA2_ERR_HIP;
    goto cleanup;
  }
  if (variant == 1) ws = rdna2_device_warp_size();
  if (ws > 0) {
    const int threads = 128;  // 4 waves of 32 or 2 waves of 64
    const int waves = threads / ws;
    const int blocks = (rows + waves - 1) / waves;
    hipLaunchKernelGGL(rdna2_dp4a_gemv_i8_wave, dim3(blocks), dim3(threads), 0, 0,
                       dy, dA, dx, rows, k4);
  } else {
    const int threads = 64;
    const int blocks = (rows + threads - 1) / threads;
    hipLaunchKernelGGL(rdna2_dp4a_gemv_i8, dim3(blocks), dim3(threads), 0, 0,
                       dy, dA, dx, rows, k4);
  }
  if (hipGetLastError() != hipSuccess ||
      hipMemcpy(y, dy, bytesy, hipMemcpyDeviceToHost) != hipSuccess) {
    rc = RDNA2_ERR_HIP;
  }
cleanup:
  if (dA) (void)hipFree(dA);
  if (dx) (void)hipFree(dx);
  if (dy) (void)hipFree(dy);
  return rc;
}

// variant 0 = one thread per row, 1 = one wave per row.
int rdna2_gemv_i8_host_variant(int* y, const int* A, const int* x, int rows, int k4, int variant) {
  if (!y || !A || !x || rows <= 0 || k4 <= 0) return RDNA2_ERR_ARG;
  if (rdna2_device_count() == 0) {
    rdna2_gemv_i8_ref(y, A, x, rows, k4);
    return RDNA2_OK;
  }
  return rdna2_gemv_i8_device(y, A, x, rows, k4, variant);
}

int rdna2_gemv_i8_host(int* y, const int* A, const int* x, int rows, int k4) {
  return rdna2_gemv_i8_host_variant(y, A, x, rows, k4, k4 >= RDNA2_GEMV_WAVE_MIN_K4 ? 1 : 0);
}

static int rdna2_sgemv_f32_device(float* y, const float* A, const float* x, int rows, int K,
                                  int variant) {
  float *dA = nullptr, *dx = nullptr, *dy = nullptr;
  const size_t bytesA = (size_t)rows * K * sizeof(float);
  const size_t bytesx = (size_t)K * sizeof(float);
  const size_t bytesy = (size_t)rows * sizeof(float);
  int rc = RDNA2_OK;
  int ws = 0;  // wave width; stays 0 unless the wave-per-row variant can size its grid
  if (hipMalloc(&dA, bytesA) != hipSuccess || hipMalloc(&dx, bytesx) != hipSuccess ||
      hipMalloc(&dy, bytesy) != hipSuccess) {
    rc = RDNA2_ERR_HIP;
    goto cleanup;
  }
  if (hipMemcpy(dA, A, bytesA, hipMemcpyHostToDevice) != hipSuccess ||
      hipMemcpy(dx, x, bytesx, hipMemcpyHostToDevice) != hipSuccess) {
    rc = RDNA2_ERR_HIP;
    goto cleanup;
  }
  if (variant == 1) ws = rdna2_device_warp_size();
  if (ws > 0) {
    const int threads = 128;  // 4 waves of 32 or 2 waves of 64
    const int waves = threads / ws;
    const int blocks = (rows + waves - 1) / waves;
    hipLaunchKernelGGL(rdna2_sgemv_f32_wave, dim3(blocks), dim3(threads), 0, 0,
                       dy, dA, dx, rows, K);
  } else {
    const int threads = 64;
    const int blocks = (rows + threads - 1) / threads;
    hipLaunchKernelGGL(rdna2_sgemv_f32, dim3(blocks), dim3(threads), 0, 0,
                       dy, dA, dx, rows, K);
  }
  if (hipGetLastError() != hipSuccess ||
      hipMemcpy(y, dy, bytesy, hipMemcpyDeviceToHost) != hipSuccess) {
    rc = RDNA2_ERR_HIP;
  }
cleanup:
  if (dA) (void)hipFree(dA);
  if (dx) (void)hipFree(dx);
  if (dy) (void)hipFree(dy);
  return rc;
}

// variant 0 = one thread per row, 1 = one wave per row.
int rdna2_sgemv_f32_host_variant(float* y, const float* A, const float* x, int rows, int K,
                                 int variant) {
  if (!y || !A || !x || rows <= 0 || K <= 0) return RDNA2_ERR_ARG;
  if (rdna2_device_count() == 0) {
    rdna2_sgemv_f32_ref(y, A, x, rows, K);
    return RDNA2_OK;
  }
  return rdna2_sgemv_f32_device(y, A, x, rows, K, variant);
}

int rdna2_sgemv_f32_host(float* y, const float* A, const float* x, int rows, int K) {
  return rdna2_sgemv_f32_host_variant(y, A, x, rows, K,
                                      K >= RDNA2_SGEMV_WAVE_MIN_K ? 1 : 0);
}

// Shapes the int8 GEMM onto a grid and launches it. Shared by the plain entry
// point and the quantized-epilogue one so both measure the same kernels.
static void rdna2_launch_gemm_i8(int* dC, const int* dA, const int* dB, int M, int N, int K4,
                                 int tiled) {
  if (tiled == 4 || tiled == 5) {
    // RDNA4 shape: 128x128 tile, 16x16 threads. 17408 B of LDS per block for
    // int8, well inside the 64 KiB workgroup limit. tiled == 5 is the same tile
    // driven by the matrix instruction instead of dp4a, so it shares the grid;
    // which instruction runs is decided inside the kernel (see
    // rdna2_wmma_iu8_step), not here -- this function is host code and cannot
    // see the offload arch.
    const dim3 grid((N + RDNA2_T4_TILE_MN - 1) / RDNA2_T4_TILE_MN,
                    (M + RDNA2_T4_TILE_MN - 1) / RDNA2_T4_TILE_MN, 1);
    const dim3 block(16, 16, 1);
    if (tiled == 5)
      hipLaunchKernelGGL(rdna2_dp4a_gemm_i8_wmma, grid, block, 0, 0, dC, dA, dB, M, N, K4);
    else
      hipLaunchKernelGGL(rdna2_dp4a_gemm_i8_lds128, grid, block, 0, 0, dC, dA, dB, M, N, K4);
  } else if (tiled) {
    const dim3 grid((N + RDNA2_TILE_MN - 1) / RDNA2_TILE_MN,
                    (M + RDNA2_TILE_MN - 1) / RDNA2_TILE_MN, 1);
    if (tiled == 3) {
      const dim3 block(16, 8, 1);  // 128 threads, 8x4 outputs each
      hipLaunchKernelGGL(rdna2_dp4a_gemm_i8_lds8, grid, block, 0, 0, dC, dA, dB, M, N, K4);
    } else if (tiled == 2) {
      const dim3 block(16, 16, 1);
      hipLaunchKernelGGL(rdna2_dp4a_gemm_i8_lds, grid, block, 0, 0, dC, dA, dB, M, N, K4);
    } else {
      const dim3 block(16, 16, 1);
      hipLaunchKernelGGL(rdna2_dp4a_gemm_i8_tile, grid, block, 0, 0, dC, dA, dB, M, N, K4);
    }
  } else {
    const dim3 block(64, 1, 1);
    const dim3 grid((N + 63) / 64, M, 1);
    hipLaunchKernelGGL(rdna2_dp4a_gemm_i8_naive, grid, block, 0, 0, dC, dA, dB, M, N, K4);
  }
}

static int rdna2_gemm_i8_device(int* C, const int* A, const int* Bt, int M, int N, int K4,
                                int tiled) {
  int *dA = nullptr, *dB = nullptr, *dC = nullptr;
  const size_t bytesA = (size_t)M * K4 * sizeof(int);
  const size_t bytesB = (size_t)N * K4 * sizeof(int);
  const size_t bytesC = (size_t)M * N * sizeof(int);
  int rc = RDNA2_OK;
  if (hipMalloc(&dA, bytesA) != hipSuccess || hipMalloc(&dB, bytesB) != hipSuccess ||
      hipMalloc(&dC, bytesC) != hipSuccess) {
    rc = RDNA2_ERR_HIP;
    goto cleanup;
  }
  if (hipMemcpy(dA, A, bytesA, hipMemcpyHostToDevice) != hipSuccess ||
      hipMemcpy(dB, Bt, bytesB, hipMemcpyHostToDevice) != hipSuccess) {
    rc = RDNA2_ERR_HIP;
    goto cleanup;
  }
  rdna2_launch_gemm_i8(dC, dA, dB, M, N, K4, tiled);
  if (hipGetLastError() != hipSuccess ||
      hipMemcpy(C, dC, bytesC, hipMemcpyDeviceToHost) != hipSuccess) {
    rc = RDNA2_ERR_HIP;
  }
cleanup:
  if (dA) (void)hipFree(dA);
  if (dB) (void)hipFree(dB);
  if (dC) (void)hipFree(dC);
  return rc;
}

// tiled: 0 = naive, 1 = 4x4 register tiled (global loads), 2 = same tiling with
// LDS-staged K chunks, 3 = LDS staging with 8x4 register blocking (quarter less
// LDS traffic for the same math), 4 = RDNA4 128x128 tile with 8x8 blocking.
// Any non-zero value keeps the tiled shape; 4 picks its own grid.
int rdna2_gemm_i8_host(int* C, const int* A, const int* Bt, int M, int N, int K4, int tiled) {
  if (!C || !A || !Bt || M <= 0 || N <= 0 || K4 <= 0) return RDNA2_ERR_ARG;
  if (rdna2_device_count() == 0) {
    rdna2_gemm_i8_ref(C, A, Bt, M, N, K4);
    return RDNA2_OK;
  }
  return rdna2_gemm_i8_device(C, A, Bt, M, N, K4, tiled);
}

// Quantized entry point: the same int8 GEMM followed by the format's scale and
// zero point. A and Bt are consumed exactly as stored - no weight is expanded to
// float, and the only float produced is the finished dot product, once. C is f32
// because a GEMM result feeds the next layer's float activation path; sa scales
// rows of A and sb columns of Bt (either may be NULL for 1.0), zp_b is B's zero
// point (0 for symmetric weights, in which case the row-sum pass is skipped).
static int rdna2_gemm_i8_scale_device(float* C, const int* A, const int* Bt, int M, int N, int K4,
                                      int tiled, const float* sa, const float* sb, int zp_b) {
  int *dA = nullptr, *dB = nullptr, *dAcc = nullptr, *dSums = nullptr;
  float *dSa = nullptr, *dSb = nullptr, *dC = nullptr;
  const size_t bytesA = (size_t)M * K4 * sizeof(int);
  const size_t bytesB = (size_t)N * K4 * sizeof(int);
  const size_t nOut = (size_t)M * N;
  int rc = RDNA2_OK;
  if (hipMalloc(&dA, bytesA) != hipSuccess || hipMalloc(&dB, bytesB) != hipSuccess ||
      hipMalloc(&dAcc, nOut * sizeof(int)) != hipSuccess ||
      hipMalloc(&dC, nOut * sizeof(float)) != hipSuccess) {
    rc = RDNA2_ERR_HIP;
    goto cleanup;
  }
  if ((zp_b != 0 && hipMalloc(&dSums, (size_t)M * sizeof(int)) != hipSuccess) ||
      (sa && hipMalloc(&dSa, (size_t)M * sizeof(float)) != hipSuccess) ||
      (sb && hipMalloc(&dSb, (size_t)N * sizeof(float)) != hipSuccess)) {
    rc = RDNA2_ERR_HIP;
    goto cleanup;
  }
  if (hipMemcpy(dA, A, bytesA, hipMemcpyHostToDevice) != hipSuccess ||
      hipMemcpy(dB, Bt, bytesB, hipMemcpyHostToDevice) != hipSuccess ||
      (sa && hipMemcpy(dSa, sa, (size_t)M * sizeof(float), hipMemcpyHostToDevice) != hipSuccess) ||
      (sb && hipMemcpy(dSb, sb, (size_t)N * sizeof(float), hipMemcpyHostToDevice) != hipSuccess)) {
    rc = RDNA2_ERR_HIP;
    goto cleanup;
  }
  rdna2_launch_gemm_i8(dAcc, dA, dB, M, N, K4, tiled);
  if (zp_b != 0) {
    const int threads = 128;
    hipLaunchKernelGGL(rdna2_i8_row_sums, dim3((M + threads - 1) / threads), dim3(threads), 0, 0,
                       dSums, dA, M, K4);
  }
  {
    const dim3 block(256, 1, 1);
    const int rows_per_grid = 4096;  // rows stride, so grid.y stays inside limits
    const dim3 grid((N + 255) / 256, M < rows_per_grid ? M : rows_per_grid, 1);
    hipLaunchKernelGGL(rdna2_gemm_i8_epilogue, grid, block, 0, 0, dC, dAcc, dSums, M, N, dSa, dSb,
                       zp_b);
  }
  if (hipGetLastError() != hipSuccess ||
      hipMemcpy(C, dC, nOut * sizeof(float), hipMemcpyDeviceToHost) != hipSuccess) {
    rc = RDNA2_ERR_HIP;
  }
cleanup:
  if (dA) (void)hipFree(dA);
  if (dB) (void)hipFree(dB);
  if (dAcc) (void)hipFree(dAcc);
  if (dSums) (void)hipFree(dSums);
  if (dSa) (void)hipFree(dSa);
  if (dSb) (void)hipFree(dSb);
  if (dC) (void)hipFree(dC);
  return rc;
}

int rdna2_gemm_i8_scale_host(float* C, const int* A, const int* Bt, int M, int N, int K4,
                             int tiled, const float* sa, const float* sb, int zp_b) {
  if (!C || !A || !Bt || M <= 0 || N <= 0 || K4 <= 0) return RDNA2_ERR_ARG;
  if (rdna2_device_count() == 0) {
    rdna2_gemm_i8_scale_ref(C, A, Bt, M, N, K4, sa, sb, zp_b);
    return RDNA2_OK;
  }
  return rdna2_gemm_i8_scale_device(C, A, Bt, M, N, K4, tiled, sa, sb, zp_b);
}

// variant 0 = 4x4 tiled reading global memory, 1 = the same tiling with LDS
// staging (fewer global reads per output), 2 = RDNA4 128x128 tile with 8x8
// blocking.
int rdna2_sgemm_f32_host_v(float* C, const float* A, const float* Bt, int M, int N, int K,
                           int variant) {
  if (!C || !A || !Bt || M <= 0 || N <= 0 || K <= 0) return RDNA2_ERR_ARG;
  if (rdna2_device_count() == 0) {
    rdna2_sgemm_f32_ref(C, A, Bt, M, N, K);
    return RDNA2_OK;
  }
  float *dA = nullptr, *dB = nullptr, *dC = nullptr;
  const size_t bytesA = (size_t)M * K * sizeof(float);
  const size_t bytesB = (size_t)N * K * sizeof(float);
  const size_t bytesC = (size_t)M * N * sizeof(float);
  int rc = RDNA2_OK;
  if (hipMalloc(&dA, bytesA) != hipSuccess || hipMalloc(&dB, bytesB) != hipSuccess ||
      hipMalloc(&dC, bytesC) != hipSuccess) {
    rc = RDNA2_ERR_HIP;
    goto cleanup;
  }
  if (hipMemcpy(dA, A, bytesA, hipMemcpyHostToDevice) != hipSuccess ||
      hipMemcpy(dB, Bt, bytesB, hipMemcpyHostToDevice) != hipSuccess) {
    rc = RDNA2_ERR_HIP;
    goto cleanup;
  }
  if (variant == 2) {
    // RDNA4 shape: 128x128 f32 tile, 8x8 outputs per thread, 9216 B of LDS.
    const dim3 block(16, 16, 1);
    const dim3 grid((N + RDNA2_T4_TILE_MN - 1) / RDNA2_T4_TILE_MN,
                    (M + RDNA2_T4_TILE_MN - 1) / RDNA2_T4_TILE_MN, 1);
    hipLaunchKernelGGL(rdna2_sgemm_f32_lds128, grid, block, 0, 0, dC, dA, dB, M, N, K);
    if (hipGetLastError() != hipSuccess ||
        hipMemcpy(C, dC, bytesC, hipMemcpyDeviceToHost) != hipSuccess) {
      rc = RDNA2_ERR_HIP;
    }
  } else {
    const dim3 block(16, 16, 1);
    const dim3 grid((N + 63) / 64, (M + 63) / 64, 1);
    if (variant == 1) {
      hipLaunchKernelGGL(rdna2_sgemm_f32_lds, grid, block, 0, 0, dC, dA, dB, M, N, K);
    } else if (variant == 3) {
      hipLaunchKernelGGL(rdna2_sgemm_f32_lds16, grid, block, 0, 0, dC, dA, dB, M, N, K);
    } else {
      hipLaunchKernelGGL(rdna2_sgemm_f32, grid, block, 0, 0, dC, dA, dB, M, N, K);
    }
    if (hipGetLastError() != hipSuccess ||
        hipMemcpy(C, dC, bytesC, hipMemcpyDeviceToHost) != hipSuccess) {
      rc = RDNA2_ERR_HIP;
    }
  }
cleanup:
  if (dA) (void)hipFree(dA);
  if (dB) (void)hipFree(dB);
  if (dC) (void)hipFree(dC);
  return rc;
}

int rdna2_sgemm_f32_host(float* C, const float* A, const float* Bt, int M, int N, int K) {
  return rdna2_sgemm_f32_host_v(C, A, Bt, M, N, K, 1);
}

int rdna2_dot_u8_host(unsigned* out, const unsigned* a, const unsigned* b, int n) {
  if (!out || !a || !b || n <= 0) return RDNA2_ERR_ARG;
  if (rdna2_device_count() == 0) {
    *out = rdna2_dot_u8_ref(a, b, n);
    return RDNA2_OK;
  }
  unsigned* da = nullptr;
  unsigned* db = nullptr;
  unsigned* dout = nullptr;
  int rc = RDNA2_OK;
  const size_t bytes = (size_t)n * sizeof(unsigned);
  if (hipMalloc(&da, bytes) != hipSuccess || hipMalloc(&db, bytes) != hipSuccess ||
      hipMalloc(&dout, sizeof(unsigned)) != hipSuccess) {
    rc = RDNA2_ERR_HIP;
    goto cleanup;
  }
  if (hipMemcpy(da, a, bytes, hipMemcpyHostToDevice) != hipSuccess ||
      hipMemcpy(db, b, bytes, hipMemcpyHostToDevice) != hipSuccess ||
      hipMemset(dout, 0, sizeof(unsigned)) != hipSuccess) {
    rc = RDNA2_ERR_HIP;
    goto cleanup;
  }
  {
    const int threads = 256;
    int blocks = (n + threads - 1) / threads;
    if (blocks > 1024) blocks = 1024;
    hipLaunchKernelGGL(rdna2_dp4a_dot_u8, dim3(blocks), dim3(threads), 0, 0, dout, da, db, n);
    if (hipGetLastError() != hipSuccess ||
        hipMemcpy(out, dout, sizeof(unsigned), hipMemcpyDeviceToHost) != hipSuccess) {
      rc = RDNA2_ERR_HIP;
    }
  }
cleanup:
  if (da) (void)hipFree(da);
  if (db) (void)hipFree(db);
  if (dout) (void)hipFree(dout);
  return rc;
}

// a and b hold n packed pairs of _Float16 (two per uint32).
int rdna2_dot_f16_host(float* out, const unsigned* a, const unsigned* b, int n) {
  if (!out || !a || !b || n <= 0) return RDNA2_ERR_ARG;
  if (rdna2_device_count() == 0) {
#if defined(__FLT16_MAX__)
    *out = rdna2_dot_f16_ref(a, b, n);
    return RDNA2_OK;
#else
    return RDNA2_ERR_HIP;  // host has no _Float16: reference unavailable
#endif
  }
  unsigned* da = nullptr;
  unsigned* db = nullptr;
  float* dout = nullptr;
  int rc = RDNA2_OK;
  const size_t bytes = (size_t)n * sizeof(unsigned);
  if (hipMalloc(&da, bytes) != hipSuccess || hipMalloc(&db, bytes) != hipSuccess ||
      hipMalloc(&dout, sizeof(float)) != hipSuccess) {
    rc = RDNA2_ERR_HIP;
    goto cleanup;
  }
  if (hipMemcpy(da, a, bytes, hipMemcpyHostToDevice) != hipSuccess ||
      hipMemcpy(db, b, bytes, hipMemcpyHostToDevice) != hipSuccess ||
      hipMemset(dout, 0, sizeof(float)) != hipSuccess) {
    rc = RDNA2_ERR_HIP;
    goto cleanup;
  }
  {
    const int threads = 256;
    int blocks = (n + threads - 1) / threads;
    if (blocks > 1024) blocks = 1024;
    hipLaunchKernelGGL(rdna2_dot2_f16, dim3(blocks), dim3(threads), 0, 0,
                       dout, reinterpret_cast<const half2_t*>(da),
                       reinterpret_cast<const half2_t*>(db), n);
    if (hipGetLastError() != hipSuccess ||
        hipMemcpy(out, dout, sizeof(float), hipMemcpyDeviceToHost) != hipSuccess) {
      rc = RDNA2_ERR_HIP;
    }
  }
cleanup:
  if (da) (void)hipFree(da);
  if (db) (void)hipFree(db);
  if (dout) (void)hipFree(dout);
  return rc;
}

// Inclusive prefix sum over a single wave; n must satisfy 1 <= n <= warpSize.
int rdna2_prefix_host(int* out, const int* in, int n) {
  if (!out || !in || n <= 0) return RDNA2_ERR_ARG;
  const int ws = rdna2_warp_size();
  if (n > ws) return RDNA2_ERR_ARG;
  if (rdna2_device_count() == 0) {
    rdna2_prefix_ref(out, in, n);
    return RDNA2_OK;
  }
  int* din = nullptr;
  int* dout = nullptr;
  int rc = RDNA2_OK;
  // dout is warpSize wide: the kernel writes one value per lane.
  if (hipMalloc(&din, (size_t)n * sizeof(int)) != hipSuccess ||
      hipMalloc(&dout, (size_t)ws * sizeof(int)) != hipSuccess) {
    rc = RDNA2_ERR_HIP;
    goto cleanup;
  }
  if (hipMemcpy(din, in, (size_t)n * sizeof(int), hipMemcpyHostToDevice) != hipSuccess) {
    rc = RDNA2_ERR_HIP;
    goto cleanup;
  }
  hipLaunchKernelGGL(rdna2_dpp_prefix, dim3(1), dim3(ws), 0, 0, dout, din, n);
  if (hipGetLastError() != hipSuccess ||
      hipMemcpy(out, dout, (size_t)n * sizeof(int), hipMemcpyDeviceToHost) != hipSuccess) {
    rc = RDNA2_ERR_HIP;
  }
cleanup:
  if (din) (void)hipFree(din);
  if (dout) (void)hipFree(dout);
  return rc;
}

// ---------------------------------------------------------------------------
// Device information, for the test's report line. arch_ok is set to one of the
// RDNA2_ARCH_* codes (1 = gfx1031, 2 = gfx1030, 3 = gfx1201, 4 = gfx1200), or 0
// when the running device is none of those.
// ---------------------------------------------------------------------------
int rdna2_device_info(char* name, size_t len, int* warp_size, int* arch_ok) {
  if (warp_size) *warp_size = 0;
  if (arch_ok) *arch_ok = 0;
  if (rdna2_device_count() == 0) {
    if (name && len) name[0] = '\0';
    return RDNA2_OK;
  }
  hipDeviceProp_t prop{};
  if (hipGetDeviceProperties(&prop, 0) != hipSuccess) return RDNA2_ERR_HIP;
  if (name && len) {
    size_t i = 0;
    for (; i + 1 < len && prop.name[i]; ++i) name[i] = prop.name[i];
    name[i] = '\0';
  }
  if (warp_size) *warp_size = prop.warpSize;
  if (arch_ok) {
    // "gfx1031", "gfx1201:xnack-" and friends: match the six target digits only.
    const char* gcn = prop.gcnArchName;
    if (gcn && gcn[0] == 'g' && gcn[1] == 'f' && gcn[2] == 'x') {
      if (gcn[3] == '1' && gcn[4] == '0' && gcn[5] == '3') {
        if (gcn[6] == '1') *arch_ok = RDNA2_ARCH_GFX1031;
        else if (gcn[6] == '0') *arch_ok = RDNA2_ARCH_GFX1030;
      } else if (gcn[3] == '1' && gcn[4] == '2' && gcn[5] == '0') {
        if (gcn[6] == '1') *arch_ok = RDNA2_ARCH_GFX1201;
        else if (gcn[6] == '0') *arch_ok = RDNA2_ARCH_GFX1200;
      }
    }
  }
  return RDNA2_OK;
}

// ---------------------------------------------------------------------------
// Architecture selection
// ---------------------------------------------------------------------------
// The 128x128 / 8x8 shapes (tiled == 4, sgemm variant == 2) were written for
// RDNA4 on the theory that the bigger register file wants a bigger tile. The
// 9070 XT disagrees: measured against the 64x64 8x4 kernel at 1024^3, 2048^2,
// 4096^2 and 2048^2 x K=4096, the 128 tile is 8-33% slower for int8 (82 VGPRs
// per thread plus 17408 B of LDS per block caps residency at 3 blocks, and the
// dp4a issue rate is already limited by LDS bandwidth), and for f32 it only
// wins while the operand working set stays inside the cache. So both RDNA4
// parts keep the 64x64 8x4 int8 shape and the 64x64 LDS f32 shape;
// tiled == 4 / variant == 2 stay selectable and are still what the test suite
// exercises on every device. Numbers are in README.md.
//
// RDNA4 then overrides that for int8: the matrix instruction (tiled == 5) on the
// same 128x128 tile beat the 64x64 8x4 dp4a kernel by 3.1x at 1024^3, 3.3x at
// 2048^2 and 3.1x at 4096^3 on the 9070 XT (29.4/40.6/40.1 TOPS ->
// 91.3/134.0/125.0), so gfx1200/gfx1201 take tiled == 5 and RDNA2 keeps dp4a.
int rdna2_tiled_for_arch(int arch) {
  if (arch == RDNA2_ARCH_GFX1200 || arch == RDNA2_ARCH_GFX1201) {
    return 5;  // 128x128 tile, one 16x16x16 matrix instruction per k step
  }
  if (arch == RDNA2_ARCH_GFX1030 || arch == RDNA2_ARCH_GFX1031) {
    return 3;  // 64x64 tile, 8x4 outputs per thread, LDS staged K chunks
  }
  return 1;  // unknown AMD GPU: 4x4 tiled, global loads
}

// variant 2 is measurably faster on RDNA4 only while the two 128-row operand
// panels stay cache-resident (K <= 2048 at 2048^2); past that it loses up to 2x
// to the 64x64 shape, so it is opt-in rather than the arch default.
// variant 3 is the 64x64 tile with a 16-deep K step (KSTEP=16, stride 17, four
// floats staged per thread). On the 9070 XT it beat variant 1 at every shape
// tried (1024^3, 1000^3, 2048^2 x K=2048/3072/4096/8192, 2304^2, 3584^2, 4096^3)
// by 3-42%, the margin shrinking as the problem grows, so it is the default on
// the RDNA2 and RDNA4 parts and the generic 64x64 shape is kept for anything
// unrecognised.
int rdna2_sgemm_variant_for_arch(int arch) {
  if (arch == RDNA2_ARCH_GFX1030 || arch == RDNA2_ARCH_GFX1031 ||
      arch == RDNA2_ARCH_GFX1200 || arch == RDNA2_ARCH_GFX1201) {
    return 3;  // 64x64 tile, 16-deep K step, four floats staged per thread
  }
  return 1;  // unknown AMD GPU: 64x64 LDS staging, 8-deep K step
}

// The arch-aware entry points: same contract as rdna2_gemm_i8_host /
// rdna2_sgemm_f32_host, but the tile shape follows the device that is actually
// present instead of a fixed default. On a machine with no GPU both fall back to
// the CPU reference exactly like the plain entry points do.
int rdna2_gemm_i8_host_tuned(int* C, const int* A, const int* Bt, int M, int N, int K4) {
  int arch = RDNA2_ARCH_UNKNOWN;
  (void)rdna2_device_info(nullptr, 0, nullptr, &arch);
  return rdna2_gemm_i8_host(C, A, Bt, M, N, K4, rdna2_tiled_for_arch(arch));
}

int rdna2_gemm_i8_scale_host_tuned(float* C, const int* A, const int* Bt, int M, int N, int K4,
                                   const float* sa, const float* sb, int zp_b) {
  int arch = RDNA2_ARCH_UNKNOWN;
  (void)rdna2_device_info(nullptr, 0, nullptr, &arch);
  return rdna2_gemm_i8_scale_host(C, A, Bt, M, N, K4, rdna2_tiled_for_arch(arch), sa, sb, zp_b);
}

int rdna2_sgemm_f32_host_tuned(float* C, const float* A, const float* Bt, int M, int N, int K) {
  int arch = RDNA2_ARCH_UNKNOWN;
  (void)rdna2_device_info(nullptr, 0, nullptr, &arch);
  return rdna2_sgemm_f32_host_v(C, A, Bt, M, N, K, rdna2_sgemm_variant_for_arch(arch));
}

