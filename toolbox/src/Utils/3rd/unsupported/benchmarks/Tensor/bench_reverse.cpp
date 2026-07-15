// Benchmarks for Eigen TensorReverse.
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#include <benchmark/benchmark.h>
#include <unsupported/Eigen/Tensor>

using namespace Eigen;

typedef float Scalar;

// --- Reverse only the inner-most (contiguous) dimension. SIMD preverse case. ---
static void BM_Reverse_Inner(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);

  Tensor<Scalar, 2> A(M, N);
  A.setRandom();

  array<bool, 2> dim_rev = {true, false};

  for (auto _ : state) {
    Tensor<Scalar, 2> B = A.reverse(dim_rev);
    benchmark::DoNotOptimize(B.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * static_cast<int64_t>(M) * N * sizeof(Scalar));
}

// --- Reverse only an outer dimension. Inner dim stays contiguous. ---
static void BM_Reverse_Outer(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);

  Tensor<Scalar, 2> A(M, N);
  A.setRandom();

  array<bool, 2> dim_rev = {false, true};

  for (auto _ : state) {
    Tensor<Scalar, 2> B = A.reverse(dim_rev);
    benchmark::DoNotOptimize(B.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * static_cast<int64_t>(M) * N * sizeof(Scalar));
}

// --- Reverse every dimension. ---
static void BM_Reverse_All(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);

  Tensor<Scalar, 2> A(M, N);
  A.setRandom();

  array<bool, 2> dim_rev = {true, true};

  for (auto _ : state) {
    Tensor<Scalar, 2> B = A.reverse(dim_rev);
    benchmark::DoNotOptimize(B.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * static_cast<int64_t>(M) * N * sizeof(Scalar));
}

// --- 3D reverse with the inner dim reversed (typical CNN-style layout). ---
static void BM_Reverse_3D_Inner(benchmark::State& state) {
  const int D0 = state.range(0);
  const int D1 = state.range(1);
  const int D2 = state.range(2);

  Tensor<Scalar, 3> A(D0, D1, D2);
  A.setRandom();

  array<bool, 3> dim_rev = {true, false, false};

  for (auto _ : state) {
    Tensor<Scalar, 3> B = A.reverse(dim_rev);
    benchmark::DoNotOptimize(B.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * static_cast<int64_t>(D0) * D1 * D2 * sizeof(Scalar));
}

// Sweep sizes that span L1 (~32 KB), L2 (~256 KB), and LLC (~MBs) for float
// tensors. Bytes per element = 4, so per-side sizes:
//   64x64    = 16 KB (L1)
//   256x256  = 256 KB (L2)
//   1024x1024 = 4 MB (LLC / DRAM)
// clang-format off
#define REVERSE_SIZES \
  ->Args({64, 64})->Args({256, 256})->Args({1024, 1024})

// 128 KB / 1 MB / 8 MB
#define REVERSE_3D_SIZES \
  ->Args({32, 32, 32})->Args({64, 64, 64})->Args({128, 128, 128})
// clang-format on

BENCHMARK(BM_Reverse_Inner) REVERSE_SIZES;
BENCHMARK(BM_Reverse_Outer) REVERSE_SIZES;
BENCHMARK(BM_Reverse_All) REVERSE_SIZES;
BENCHMARK(BM_Reverse_3D_Inner) REVERSE_3D_SIZES;
