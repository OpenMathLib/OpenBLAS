// CPU references for the RDNA2 kernels. Plain C++ (no HIP), so every
// reference can be built and tested on a machine without a GPU. The launchers
// in rdna2_kernels.hip use these as the fallback when no device is present,
// and the test uses them as the correctness oracle.
//
// Packed int8 convention (matches v_dot4_*): four int8 values per int32,
// value i at bits [8*i, 8*i+8), i.e. little-endian byte order.
#ifndef RDNA2_GEMM_REF_H
#define RDNA2_GEMM_REF_H

#include <cstdint>
#include <cstddef>
#include <cstring>

// Ref functions that device code calls must be __host__ __device__. Under a
// plain host compiler the macro expands to nothing.
#ifndef RDNA2_HD
#if defined(__HIP__) || defined(__HIPCC__) || defined(__CUDACC__)
#define RDNA2_HD __host__ __device__
#else
#define RDNA2_HD
#endif
#endif

// ---------------------------------------------------------------------------
// dp4a int8 / u8
// ---------------------------------------------------------------------------
RDNA2_HD static inline int32_t rdna2_dot4_i8_ref(int32_t a, int32_t b, int32_t acc) {
  for (int i = 0; i < 4; ++i) {
    int8_t x = static_cast<int8_t>((static_cast<uint32_t>(a) >> (8 * i)) & 0xffu);
    int8_t y = static_cast<int8_t>((static_cast<uint32_t>(b) >> (8 * i)) & 0xffu);
    uint32_t p = static_cast<uint32_t>(static_cast<int32_t>(x) * static_cast<int32_t>(y));
    acc = static_cast<int32_t>(static_cast<uint32_t>(acc) + p);
  }
  return acc;
}

static inline uint32_t rdna2_dot4_u8_ref(uint32_t a, uint32_t b, uint32_t acc) {
  for (int i = 0; i < 4; ++i) {
    uint32_t x = (a >> (8 * i)) & 0xffu;
    uint32_t y = (b >> (8 * i)) & 0xffu;
    acc += x * y;
  }
  return acc;
}

// The dot4 primitive used by the kernel bodies. On the device this is the real
// v_dot4c_i32_i8 (dp4a) instruction, independent of what the optimiser can
// pattern-match; on the host it is the bit-exact scalar reference above, so the
// same kernel body still runs on a CPU-only machine.
RDNA2_HD static inline int32_t rdna2_dot4_i8_impl(int32_t a, int32_t b, int32_t acc) {
#if defined(__HIP_DEVICE_COMPILE__) && !defined(RDNA2_NO_DOT1)
#if defined(__gfx1200__) || defined(__gfx1201__)
  // gfx12 (RDNA4) has no dot1-insts, so `__builtin_amdgcn_sdot4` is rejected
  // ("needs target feature dot1-insts"). The dp4a hardware is still there: the
  // legacy v_dot4_i32_i8 mnemonic assembles and computes exact signed results
  // (checked against the scalar reference on a gfx1201 card). Note that
  // v_dot4_i32_iu8 is *not* interchangeable: on gfx12 it is the unsigned form.
  asm volatile("v_dot4_i32_i8 %0, %1, %2, %0" : "+v"(acc) : "v"(a), "v"(b));
  return acc;
#else
  return __builtin_amdgcn_sdot4(a, b, acc, false);
#endif
#else
  // Either a host build or -DRDNA2_NO_DOT1=1 (a target whose dot1-insts feature
  // is missing, e.g. an older gfx12 toolchain): same arithmetic, four scalar
  // multiplies. Both branches are compared in the host test.
  return rdna2_dot4_i8_ref(a, b, acc);
#endif
}

// ---------------------------------------------------------------------------
// GEMV: y[M] = A[M x K4] . x[K4], packed int8
// ---------------------------------------------------------------------------
static inline void rdna2_gemv_i8_ref(int32_t* y, const int32_t* A, const int32_t* x,
                                     int rows, int K4) {
  for (int r = 0; r < rows; ++r) {
    int32_t acc = 0;
    for (int k = 0; k < K4; ++k) acc = rdna2_dot4_i8_ref(A[(size_t)r * K4 + k], x[k], acc);
    y[r] = acc;
  }
}

// ---------------------------------------------------------------------------
// GEMV: y[M] = A[M x K] . x[K], float
// ---------------------------------------------------------------------------
static inline void rdna2_sgemv_f32_ref(float* y, const float* A, const float* x,
                                      int rows, int K) {
  for (int r = 0; r < rows; ++r) {
    // The reference accumulates in double and rounds once, so it does not
    // share the kernel's four-way partial-sum split; the test compares with a
    // tolerance rather than bit-for-bit.
    double acc = 0.0;
    for (int k = 0; k < K; ++k) acc += (double)A[(size_t)r * K + k] * (double)x[k];
    y[r] = (float)acc;
  }
}

// ---------------------------------------------------------------------------
// GEMM: C[M x N] = A[M x K4] * Bt[N x K4]^T, packed int8
// ---------------------------------------------------------------------------
static inline void rdna2_gemm_i8_ref(int32_t* C, const int32_t* A, const int32_t* Bt,
                                     int M, int N, int K4) {
  for (int m = 0; m < M; ++m) {
    for (int n = 0; n < N; ++n) {
      const int32_t* a = A + static_cast<size_t>(m) * K4;
      const int32_t* b = Bt + static_cast<size_t>(n) * K4;
      int32_t acc = 0;
      for (int k = 0; k < K4; ++k) acc = rdna2_dot4_i8_ref(a[k], b[k], acc);
      C[static_cast<size_t>(m) * N + n] = acc;
    }
  }
}

// ---------------------------------------------------------------------------
// GEMM with the quantized epilogue applied on the CPU:
//   C[m][n] = (acc[m][n] - zp_b * rowsum_a[m]) * sa[m] * sb[n]
// rowsum_a[m] is the sum of the int8 activations in row m. That is the exact
// integer form of B's asymmetric storage - B's stored value at k represents
// (q - zp_b), so sum_k a*(q - zp_b) = acc - zp_b*sum_k a - and it never converts
// a stored weight: only the finished dot product becomes float. sa / sb may be
// NULL (unit scale). The correction and the product are computed in 64-bit so
// the reference cannot overflow before its single float conversion.
// ---------------------------------------------------------------------------
static inline void rdna2_gemm_i8_scale_ref(float* C, const int32_t* A, const int32_t* Bt,
                                          int M, int N, int K4, const float* sa,
                                          const float* sb, int zp_b) {
  for (int m = 0; m < M; ++m) {
    const int32_t* a = A + static_cast<size_t>(m) * K4;
    int64_t rowsum = 0;
    if (zp_b != 0) {
      for (int k = 0; k < K4; ++k) {
        for (int b = 0; b < 4; ++b) {
          rowsum += static_cast<int8_t>((static_cast<uint32_t>(a[k]) >> (8 * b)) & 0xffu);
        }
      }
    }
    const float s_a = sa ? sa[m] : 1.0f;
    for (int n = 0; n < N; ++n) {
      const int32_t* b = Bt + static_cast<size_t>(n) * K4;
      int32_t acc = 0;
      for (int k = 0; k < K4; ++k) acc = rdna2_dot4_i8_ref(a[k], b[k], acc);
      const int64_t corr = static_cast<int64_t>(acc) - static_cast<int64_t>(zp_b) * rowsum;
      C[static_cast<size_t>(m) * N + n] =
          static_cast<float>(corr) * s_a * (sb ? sb[n] : 1.0f);
    }
  }
}

// ---------------------------------------------------------------------------
// SGEMM: C[M x N] = A[M x K] * Bt[N x K]^T, float
// ---------------------------------------------------------------------------
static inline void rdna2_sgemm_f32_ref(float* C, const float* A, const float* Bt,
                                       int M, int N, int K) {
  for (int m = 0; m < M; ++m) {
    for (int n = 0; n < N; ++n) {
      const float* a = A + static_cast<size_t>(m) * K;
      const float* b = Bt + static_cast<size_t>(n) * K;
      float acc = 0.0f;
      for (int k = 0; k < K; ++k) acc += a[k] * b[k];
      C[static_cast<size_t>(m) * N + n] = acc;
    }
  }
}

// ---------------------------------------------------------------------------
// f16 dot (v_dot2c_f32_f16): acc = sum over i of a.xy * b.xy, f32 accumulate.
// Pairs are packed two _Float16 per int32, low half first.
// ---------------------------------------------------------------------------
#if defined(__FLT16_MAX__)
static inline void rdna2_unpack_half2(uint32_t v, _Float16* lo, _Float16* hi) {
  uint16_t l = static_cast<uint16_t>(v & 0xffffu);
  uint16_t h = static_cast<uint16_t>(v >> 16);
  std::memcpy(lo, &l, sizeof(uint16_t));
  std::memcpy(hi, &h, sizeof(uint16_t));
}

static inline uint32_t rdna2_pack_half2(_Float16 lo, _Float16 hi) {
  uint16_t l, h;
  std::memcpy(&l, &lo, sizeof(uint16_t));
  std::memcpy(&h, &hi, sizeof(uint16_t));
  return static_cast<uint32_t>(l) | (static_cast<uint32_t>(h) << 16);
}

// Reference matches the hardware order: acc += a.x*b.x then acc += a.y*b.y.
static inline float rdna2_dot_f16_ref(const uint32_t* a, const uint32_t* b, int n) {
  float acc = 0.0f;
  for (int i = 0; i < n; ++i) {
    _Float16 ax, ay, bx, by;
    rdna2_unpack_half2(a[i], &ax, &ay);
    rdna2_unpack_half2(b[i], &bx, &by);
    acc += static_cast<float>(ax) * static_cast<float>(bx);
    acc += static_cast<float>(ay) * static_cast<float>(by);
  }
  return acc;
}
#endif

// ---------------------------------------------------------------------------
// u8 dot: out = sum over i of dot4_u8(a[i], b[i])
// ---------------------------------------------------------------------------
static inline uint32_t rdna2_dot_u8_ref(const uint32_t* a, const uint32_t* b, int n) {
  uint32_t acc = 0;
  for (int i = 0; i < n; ++i) acc = rdna2_dot4_u8_ref(a[i], b[i], acc);
  return acc;
}

// ---------------------------------------------------------------------------
// 32-lane inclusive prefix sum (matches rdna2_dpp_prefix on wave32).
// ---------------------------------------------------------------------------
static inline void rdna2_prefix_ref(int32_t* out, const int32_t* in, int n) {
  int32_t run = 0;
  for (int i = 0; i < n; ++i) {
    run += in[i];
    out[i] = run;
  }
}

#endif  // RDNA2_GEMM_REF_H
