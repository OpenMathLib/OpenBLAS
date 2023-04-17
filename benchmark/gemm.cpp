/***************************************************************************
Copyright (c) 2014, 2023. The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#include <algorithm>
#include <iostream>
#include <random>

#include <common.h>

#include "nanobench.h"

#undef GEMM

#ifndef COMPLEX
#ifdef DOUBLE
#define GEMM BLASFUNC(dgemm)
#elif defined(HALF)
#define GEMM BLASFUNC(sbgemm)
#else
#define GEMM BLASFUNC(sgemm)
#endif
#else
#ifdef DOUBLE
#define GEMM BLASFUNC(zgemm)
#else
#define GEMM BLASFUNC(cgemm)
#endif
#endif

template <typename T> static void fill_vector(std::vector<T> vec) {
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_real_distribution<T> distribution(std::numeric_limits<T>::min(),
                                                 std::numeric_limits<T>::max());

  std::generate(vec.begin(), vec.end(),
                [&]() { return distribution(generator); });
}

static std::pair<bool, std::string>
env_param(const std::string &name, const std::string &default_value) {
  const char *value = getenv(name.c_str());
  return {value == nullptr, value ? value : default_value};
}

static std::string env_value(const std::string &name,
                             const std::string &default_value) {
  return env_param(name, default_value).second;
}

int main(int argc, char *argv[]) {
  int from = (argc > 1) ? atol(argv[1]) : 1;
  int to = (argc > 2) ? MAX(atol(argv[2]), from) : 200;
  int step = (argc > 3) ? atol(argv[3]) : 1;

  FLOAT alpha[] = {1.0, 0.0};
  FLOAT beta[] = {0.0, 0.0};

  int epochs = atoi(env_value("OPENBLAS_EPOCHS", "1").c_str());
  bool json_output = env_value("OPENBLAS_OUTPUT_JSON", "0").front() == '1';

  std::pair<bool, std::string> param_m = env_param("OPENBLAS_PARAM_M", "100");
  std::pair<bool, std::string> param_n = env_param("OPENBLAS_PARAM_N", "100");
  std::pair<bool, std::string> param_k = env_param("OPENBLAS_PARAM_K", "100");
  blasint m = param_m.first ? atoi(param_m.second.c_str()) : to;
  blasint n = param_n.first ? atoi(param_n.second.c_str()) : to;
  blasint k = param_k.first ? atoi(param_k.second.c_str()) : to;

  char transpose = toupper(env_value("OPENBLAS_TRANS", "N").front());
  char transpose_a = toupper(env_value("OPENBLAS_TRANSA", "N").front());
  char transpose_b = toupper(env_value("OPENBLAS_TRANSB", "N").front());

  bool is_specific_size = param_m.first && param_n.first && param_k.first;
  if (is_specific_size) {
    from = 1;
    to = 1;
    step = 1;
  }

  std::vector<IFLOAT> a(m * k);
  std::vector<IFLOAT> b(n * k);
  std::vector<FLOAT> c(m * n);
  fill_vector(a);
  fill_vector(b);
  fill_vector(c);

  if (!is_specific_size) {
    std::cout << "From: " << std::to_string(from) << " To: " << std::to_string(to)
              << " Step: " << std::to_string(step) << " TransA: " << transpose_a
              << " TransB: " << transpose_b << "\n";
  } else {
    std::cout << "M: " << std::to_string(m) << " N: " << std::to_string(n)
              << " K: " << std::to_string(k) << " TransA: " << transpose_a
              << " TransB: " << transpose_b << "\n";
  }

  for (int i = from; i <= to; i += step) {
    if (!param_m.first) {
      m = i;
    }
    if (!param_n.first) {
      n = i;
    }
    if (!param_k.first) {
      k = i;
    }

    blasint lda = transpose == 'N' && transpose_a == 'N' ? m : k;
    blasint ldb = transpose == 'N' && transpose_b == 'N' ? k : n;
    blasint ldc = m;

    ankerl::nanobench::Bench bench;
    if (json_output) {
      bench.output(nullptr);
    }

    std::string bench_name = "M=" + std::to_string(m) +
                             " N=" + std::to_string(n) +
                             " K=" + std::to_string(k);
    bench.minEpochIterations(epochs).run(bench_name, [&]() {
      GEMM(&transpose_a, &transpose_b, &m, &n, &k, alpha, a.data(), &lda,
           b.data(), &ldb, beta, c.data(), &ldc);
    });
    if (json_output) {
      bench.render(ankerl::nanobench::templates::json(), std::cout);
    }
  }
}
