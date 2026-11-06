// CPU vs GPU comparison for the RDNA2 / RDNA4 HIP GEMM kernels.
//
// Everything below runs the *same* operands through both sides and reports
// measured time, GFLOPS and the GPU's numerical error against the CPU result,
// so a speedup is only ever claimed when this program printed one.
//
// GPU timings are end-to-end: the public host entry points allocate, copy in,
// launch, and copy out, exactly like a full cblas_sgemm call is end-to-end on
// the CPU. Scaling matrices to a device that cannot read host memory is not a
// transfer we can hide, so it is included rather than quietly excluded; the
// copy-only time is printed separately so the split is visible.
//
// The GPU columns are only produced when the library itself can see a device.
// Its host entry points run a scalar CPU reference when they cannot, and they
// return success while doing it, so no timing distinguishes that from a slow
// GPU: this program asks the library first and exits instead of printing it.
//
// CPU baselines:
//   cpu1t  - this file's blocked -O3 single-thread f32 gemm (always available)
//   openblas - the OpenBLAS DLL shipped inside numpy, multithreaded, when the
//              path is given or found. This is a tuned, threaded BLAS: the fair
//              floor for the GPU to beat.
//
// Usage:
//   bench_rdna2 [openblas.dll] [max_size]
//
#include <hip/hip_runtime.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <dlfcn.h>
#endif

extern "C" int rdna2_sgemm_f32_host_tuned(float* C, const float* A, const float* Bt, int M, int N,
                                          int K);
extern "C" int rdna2_gemm_i8_host_tuned(int* C, const int* A, const int* Bt, int M, int N, int K4);
// sa is the M-element row scale (NULL for unit), sb the N-element column scale
// (NULL for unit); the zero point is B's, and the correction is subtracted from
// the finished int32 dot, never from a stored weight.
extern "C" int rdna2_gemm_i8_scale_host_tuned(float* C, const int* A, const int* Bt, int M, int N,
                                              int K4, const float* sa, const float* sb, int zp_b);
// The library's own view of the device: the same hipGetDeviceCount() the host
// entry points use, so this is what decides whether they compute or fall back.
extern "C" int rdna2_device_info(char* name, size_t len, int* warp_size, int* arch_ok);

// ---------------------------------------------------------------------------
// CPU reference: blocked f32 GEMM, C[MxN] = A[MxK] * Bt[NxK]^T
//
// Bt is packed row-major N x K (one row per output column), which is the layout
// the GPU kernels take. A tuned BLAS also repacks operands before its main
// loop, so transposing Bt into a K x N panel once, outside the timed region,
// mirrors that and keeps the loop below a plain k-outer accumulation.
// ---------------------------------------------------------------------------
static void cpu_sgemm_blocked(float* C, const float* A, const float* Bt_kn, int M, int N, int K) {
  const int BS = 64;
  for (int i = 0; i < M * N; ++i) C[i] = 0.0f;
  for (int ii = 0; ii < M; ii += BS) {
    const int iend = std::min(ii + BS, M);
    for (int kk = 0; kk < K; kk += BS) {
      const int kend = std::min(kk + BS, K);
      for (int jj = 0; jj < N; jj += BS) {
        const int jend = std::min(jj + BS, N);
        for (int i = ii; i < iend; ++i) {
          float* c = C + (size_t)i * N;
          const float* a = A + (size_t)i * K;
          for (int k = kk; k < kend; ++k) {
            const float av = a[k];
            const float* b = Bt_kn + (size_t)k * N;
            for (int j = jj; j < jend; ++j) c[j] += av * b[j];
          }
        }
      }
    }
  }
}

// CPU reference for the int8 path: int32 accumulation of 4 dot4 groups, exactly
// the arithmetic the dp4a kernels perform, in scalar form.
static void cpu_gemm_i8(int* C, const int* A, const int* Bt, int M, int N, int K4) {
  for (int i = 0; i < M; ++i) {
    const int* a = A + (size_t)i * K4;
    for (int j = 0; j < N; ++j) {
      const int* b = Bt + (size_t)j * K4;
      int acc = 0;
      for (int k = 0; k < K4; ++k) {
        int p = a[k], q = b[k];
        for (int t = 0; t < 4; ++t) {
          int8_t x = (int8_t)((p >> (8 * t)) & 0xFF);
          int8_t y = (int8_t)((q >> (8 * t)) & 0xFF);
          acc += (int)x * (int)y;
        }
      }
      C[(size_t)i * N + j] = acc;
    }
  }
}

// ---------------------------------------------------------------------------
// Optional OpenBLAS baseline, loaded at run time so the bench still builds and
// runs where no such library exists.
// ---------------------------------------------------------------------------
typedef void (*cblas_sgemm64_fn)(int order, int transA, int transB, int64_t M, int64_t N,
                                 int64_t K, float alpha, const float* A, int64_t lda,
                                 const float* B, int64_t ldb, float beta, float* C, int64_t ldc);

struct OpenBlas {
  cblas_sgemm64_fn sgemm = nullptr;
  void* lib = nullptr;
  std::string path;

  bool open(const char* p) {
    if (!p || !*p) return false;
#ifdef _WIN32
    lib = (void*)LoadLibraryA(p);
    if (!lib) return false;
    sgemm = (cblas_sgemm64_fn)GetProcAddress((HMODULE)lib, "scipy_cblas_sgemm64_");
    if (!sgemm) sgemm = (cblas_sgemm64_fn)GetProcAddress((HMODULE)lib, "cblas_sgemm");
#else
    lib = dlopen(p, RTLD_NOW);
    if (!lib) return false;
    sgemm = (cblas_sgemm64_fn)dlsym(lib, "scipy_cblas_sgemm64_");
    if (!sgemm) sgemm = (cblas_sgemm64_fn)dlsym(lib, "cblas_sgemm");
#endif
    path = p;
    return sgemm != nullptr;
  }
  // C[MxN] = A[MxK] * Bt[NxK]^T
  void run(float* C, const float* A, const float* Bt, int M, int N, int K) {
    // CblasRowMajor=101, CblasNoTrans=111, CblasTrans=112
    sgemm(101, 111, 112, M, N, K, 1.0f, A, K, Bt, K, 0.0f, C, N);
  }
};

// ---------------------------------------------------------------------------
static double now_s() {
  using clock = std::chrono::steady_clock;
  return std::chrono::duration<double>(clock::now().time_since_epoch()).count();
}

template <typename F>
static double best_of(int reps, F&& f) {
  double best = 1e30;
  for (int r = 0; r < reps; ++r) {
    const double t0 = now_s();
    f();
    const double dt = now_s() - t0;
    best = std::min(best, dt);
  }
  return best;
}

static double gflops(double secs, int M, int N, int K) {
  return (2.0 * (double)M * (double)N * (double)K) / (secs * 1e9);
}

// best_of() for the launchers. They return 0 on success and non-zero when the
// call could not be made on the device, and a failed call returns immediately,
// so timing one would print as an excellent result. Stop at the first non-zero
// return code instead.
template <typename F>
static double best_of_gpu(int reps, const char* op, int n, F&& f) {
  double best = 1e30;
  for (int r = 0; r < reps; ++r) {
    const double t0 = now_s();
    const int rc = f();
    const double dt = now_s() - t0;
    if (rc != 0) {
      printf("bench: %s at n=%d returned %d - not a GPU run, aborting\n", op, n, rc);
      exit(2);
    }
    best = std::min(best, dt);
  }
  return best;
}

// Without a device the launchers run the scalar CPU reference and report
// success, so no timing can tell that apart from a slow GPU. Ask the library
// what it sees and refuse to run at all when the answer is "no device".
static void require_library_device(void) {
  char name[256] = {0};
  int warp = 0, arch = 0;
  const int rc = rdna2_device_info(name, sizeof(name), &warp, &arch);
  printf("bench: library sees rc=%d name='%s' warp=%d arch=%d\n", rc, name, warp, arch);
  if (rc != 0 || arch == 0 || name[0] == '\0') {
    printf(
        "bench: the library cannot see a device, so every GPU column below would be\n"
        "bench: the CPU reference. Refusing to report that as a GPU result. Run this\n"
        "bench: binary where the HIP runtime finds the card - on the ROCm Windows\n"
        "bench: build, from the bin directory that holds amdhip64_*.dll.\n");
    exit(1);
  }
}

static void fill_rand(float* p, size_t n, unsigned seed) {
  unsigned s = seed;
  for (size_t i = 0; i < n; ++i) {
    s = s * 1664525u + 1013904223u;
    p[i] = (float)((s >> 8) & 0xFFFF) / 65535.0f - 0.5f;
  }
}

static void fill_rand_i8_packed(int* p, size_t n, unsigned seed) {
  unsigned s = seed;
  for (size_t i = 0; i < n; ++i) {
    unsigned w = 0;
    for (int t = 0; t < 4; ++t) {
      s = s * 1664525u + 1013904223u;
      w |= ((s >> 16) & 0xFF) << (8 * t);
    }
    p[i] = (int)w;
  }
}

static void report_error(const char* label, const float* got, const float* want, size_t n) {
  double maxabs = 0.0, num = 0.0, den = 0.0, maxrel = 0.0;
  for (size_t i = 0; i < n; ++i) {
    const double d = (double)got[i] - (double)want[i];
    const double ad = d < 0 ? -d : d;
    if (ad > maxabs) maxabs = ad;
    num += d * d;
    den += (double)want[i] * (double)want[i];
    const double ref = want[i] < 0 ? -(double)want[i] : (double)want[i];
    if (ref > 1e-3) {
      const double rel = ad / ref;
      if (rel > maxrel) maxrel = rel;
    }
  }
  const double rms = den > 0 ? std::sqrt(num / den) : 0.0;
  printf("  %-10s vs %-10s max_abs=%.3e  max_rel=%.3e  rms_rel=%.3e\n", label, "reference",
         maxabs, maxrel, rms);
}

int main(int argc, char** argv) {
  setvbuf(stdout, nullptr, _IONBF, 0);
  const char* ob_path = (argc > 1) ? argv[1] : nullptr;
  int max_size = (argc > 2) ? atoi(argv[2]) : 4096;
  if (max_size <= 0) max_size = 4096;

  printf("bench: hipsgemm cpu-vs-gpu\n");
  char devname[256] = {0};
  int ndev = 0;
  const hipError_t dc = hipGetDeviceCount(&ndev);
  printf("bench: hipGetDeviceCount -> %d, ndev=%d\n", (int)dc, ndev);
  if (dc == hipSuccess && ndev > 0) {
    hipDeviceProp_t prop{};
    if (hipGetDeviceProperties(&prop, 0) == hipSuccess) {
      snprintf(devname, sizeof(devname), "%s", prop.name);
      printf("GPU  : %s  (gcnArchName=%s, %d CUs)\n", prop.name, prop.gcnArchName,
             prop.multiProcessorCount);
    }
  } else {
    printf("GPU  : none detected by the HIP runtime in this process\n");
  }
  require_library_device();

  OpenBlas ob;
  if (!ob.open(ob_path)) {
    printf("CPU  : openblas baseline NOT loaded (pass its DLL path as argv[1])\n");
  } else {
    printf("CPU  : openblas baseline %s\n", ob.path.c_str());
  }

  printf("\nf32 SGEMM  C[MxN] = A[MxK] * B[NxK]^T   (GPU column is end-to-end, transfers included)\n");
  printf("(cpu1t, openblas and gpu are GFLOPS; the last column is the same rate in TF/s)\n");
  printf("%-16s %10s %10s %10s %10s %10s %10s\n", "M=N=K", "cpu1t", "openblas", "gpu",
         "gpu/cpu1t", "gpu/oblas", "gpu_TF/s");

  const int sizes[] = {256, 512, 1024, 2048, 4096};
  for (int si = 0; si < 5; ++si) {
    const int n = sizes[si];
    if (n > max_size) break;
    const size_t na = (size_t)n * n;
    std::vector<float> A(na), Bt(na), Bt_kn(na), C_cpu(na), C_ob(na), C_gpu(na);
    fill_rand(A.data(), na, 12345u + n);
    fill_rand(Bt.data(), na, 999u + n);
    for (int j = 0; j < n; ++j)
      for (int k = 0; k < n; ++k) Bt_kn[(size_t)k * n + j] = Bt[(size_t)j * n + k];

    const int reps = n >= 2048 ? 3 : 5;
    const double t_cpu =
        best_of(reps, [&] { cpu_sgemm_blocked(C_cpu.data(), A.data(), Bt_kn.data(), n, n, n); });
    double t_ob = 0.0;
    if (ob.sgemm) {
      ob.run(C_ob.data(), A.data(), Bt.data(), n, n, n);  // warm up threads
      t_ob = best_of(reps, [&] { ob.run(C_ob.data(), A.data(), Bt.data(), n, n, n); });
    }
    const double t_gpu = best_of_gpu(reps, "sgemm_f32", n, [&] {
      return rdna2_sgemm_f32_host_tuned(C_gpu.data(), A.data(), Bt.data(), n, n, n);
    });

    const double gf_cpu = gflops(t_cpu, n, n, n), gf_gpu = gflops(t_gpu, n, n, n);
    const double gf_ob = t_ob > 0 ? gflops(t_ob, n, n, n) : 0.0;

    char label[32];
    snprintf(label, sizeof(label), "%d", n);
    printf("%-16s %10.1f %10.1f %10.1f %9.2fx %9s %10.3f\n", label, gf_cpu, gf_ob, gf_gpu,
           gf_gpu / gf_cpu, t_ob > 0 ? (std::to_string(gf_gpu / gf_ob) + "x").c_str() : "-",
           gf_gpu / 1e3);

    if (n == 1024 || n == max_size) {
      report_error("gpu", C_gpu.data(), C_cpu.data(), na);
      if (ob.sgemm) report_error("openblas", C_ob.data(), C_cpu.data(), na);
    }
  }

  // The host entry points allocate and copy on every call, so this is the part
  // of the GPU column that is not math. Printed so the split is visible instead
  // of being argued about.
  printf("\nGPU-side fixed cost per call (hipMalloc x3 + copy in + copy out + free)\n");
  printf("%-16s %10s %10s\n", "M=N=K", "ms", "GB/s eff");
  for (int si = 0; si < 5; ++si) {
    const int n = sizes[si];
    if (n > max_size) break;
    const size_t bytes = (size_t)n * n * sizeof(float);
    std::vector<float> A(bytes / sizeof(float)), C(bytes / sizeof(float));
    const double t = best_of(3, [&] {
      float *dA = nullptr, *dB = nullptr, *dC = nullptr;
      if (hipMalloc(&dA, bytes) != hipSuccess || hipMalloc(&dB, bytes) != hipSuccess ||
          hipMalloc(&dC, bytes) != hipSuccess ||
          hipMemcpy(dA, A.data(), bytes, hipMemcpyHostToDevice) != hipSuccess ||
          hipMemcpy(dB, A.data(), bytes, hipMemcpyHostToDevice) != hipSuccess ||
          hipMemcpy(C.data(), dC, bytes, hipMemcpyDeviceToHost) != hipSuccess) {
        printf("bench: hipMalloc/hipMemcpy failed at n=%d - no copy-only number\n", n);
        exit(3);
      }
      (void)hipFree(dA);
      (void)hipFree(dB);
      (void)hipFree(dC);
    });
    const double gbs = (3.0 * (double)bytes) / (t * 1e9);
    char label[32];
    snprintf(label, sizeof(label), "%d", n);
    printf("%-16s %10.3f %10.1f\n", label, t * 1e3, gbs);
  }

  // ---------------------------------------------------------------------
  // int8: packed int8 weights stay packed; only the accumulator leaves the
  // dot in int32/float. GPU vs scalar CPU dp4a-equivalent arithmetic.
  // ---------------------------------------------------------------------
  printf("\nint8 GEMM  int32 accumulate, 4 int8 per word (dp4a) - K4 words = 4*K4 columns\n");
  printf("%-16s %10s %10s %10s %12s\n", "M=N=K", "cpu1t", "gpu", "gpu/cpu1t", "gpu_Tops/s");
  const int isizes[] = {256, 512, 1024, 2048};
  for (int si = 0; si < 4; ++si) {
    const int n = isizes[si];
    if (n > max_size) break;
    const size_t na = (size_t)n * n;
    std::vector<int> A(na), Bt(na), C_cpu(na), C_gpu(na);
    fill_rand_i8_packed(A.data(), na, 777u + n);
    fill_rand_i8_packed(Bt.data(), na, 555u + n);

    const int reps = n >= 2048 ? 3 : 5;
    const double t_cpu =
        best_of(reps, [&] { cpu_gemm_i8(C_cpu.data(), A.data(), Bt.data(), n, n, n); });
    const double t_gpu = best_of_gpu(reps, "gemm_i8", n, [&] {
      return rdna2_gemm_i8_host_tuned(C_gpu.data(), A.data(), Bt.data(), n, n, n);
    });

    const double ops = 2.0 * 4.0 * (double)n * (double)n * (double)n;  // 2 flops per int8 MAC
    const double tops_cpu = ops / (t_cpu * 1e12), tops_gpu = ops / (t_gpu * 1e12);
    char label[32];
    snprintf(label, sizeof(label), "%d", n);
    printf("%-16s %10.2f %10.2f %9.2fx %12.2f\n", label, tops_cpu, tops_gpu, tops_gpu / tops_cpu,
           tops_gpu);

    long long maxdiff = 0;
    for (size_t i = 0; i < na; ++i) {
      long long d = (long long)C_gpu[i] - (long long)C_cpu[i];
      if (d < 0) d = -d;
      if (d > maxdiff) maxdiff = d;
    }
    printf("  int8 gpu vs cpu: max |diff| = %lld (exact match: %s)\n", maxdiff,
           maxdiff == 0 ? "yes" : "NO");
  }

  // Quantized path: int8 weights stay packed through the whole GEMM, and the
  // scale / zero point is applied once to each finished int32 dot.
  {
    const int n = 1024, K4 = n;
    const size_t na = (size_t)n * K4;
    std::vector<int> A(na), Bt(na), Ci((size_t)n * n);
    std::vector<float> C_gpu((size_t)n * n), C_ref((size_t)n * n);
    std::vector<float> sa(n), sb(n);
    fill_rand_i8_packed(A.data(), na, 4242u);
    fill_rand_i8_packed(Bt.data(), na, 2424u);
    for (int i = 0; i < n; ++i) {
      sa[i] = 1.0f / 127.0f;
      sb[i] = 1.0f / 63.0f;
    }
    const int zp_b = 3;

    const double t_gpu = best_of_gpu(3, "gemm_i8_scale", n, [&] {
      return rdna2_gemm_i8_scale_host_tuned(C_gpu.data(), A.data(), Bt.data(), n, n, K4,
                                            sa.data(), sb.data(), zp_b);
    });

    // Reference: C[m][n] = (acc - zp_b * rowsum_a[m]) * sa[m] * sb[n], exactly
    // what the epilogue does, with the int8 values never dequantised.
    cpu_gemm_i8(Ci.data(), A.data(), Bt.data(), n, n, K4);
    std::vector<long long> rowsum(n, 0);
    for (int m = 0; m < n; ++m) {
      const int* a = A.data() + (size_t)m * K4;
      long long s = 0;
      for (int k = 0; k < K4; ++k)
        for (int t = 0; t < 4; ++t) s += (int)(int8_t)((a[k] >> (8 * t)) & 0xFF);
      rowsum[m] = s;
    }
    for (int m = 0; m < n; ++m) {
      for (int j = 0; j < n; ++j) {
        long long v = (long long)Ci[(size_t)m * n + j] - (long long)zp_b * rowsum[m];
        C_ref[(size_t)m * n + j] = (float)((double)v * (double)sa[m] * (double)sb[j]);
      }
    }
    report_error("quant epi", C_gpu.data(), C_ref.data(), C_gpu.size());
    printf("  quant epilogue 1024^3: %.3f ms (%.1f GFLOPS)\n", t_gpu * 1e3,
           gflops(t_gpu, n, n, n));
  }

  return 0;
}
