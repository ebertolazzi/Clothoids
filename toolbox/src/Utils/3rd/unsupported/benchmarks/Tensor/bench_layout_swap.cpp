// Benchmarks for Eigen TensorLayoutSwap.
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#include <benchmark/benchmark.h>
#include <unsupported/Eigen/Tensor>

using namespace Eigen;

typedef float Scalar;

static void BM_LayoutSwap_2D(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);

  Tensor<Scalar, 2, ColMajor> A(M, N);
  A.setRandom();

  for (auto _ : state) {
    Tensor<Scalar, 2, RowMajor> B = A.swap_layout();
    benchmark::DoNotOptimize(B.data());
    benchmark::ClobberMemory();
  }
  // 1 read (A) + 1 write (B).
  state.SetBytesProcessed(state.iterations() * 2ll * static_cast<int64_t>(M) * N * sizeof(Scalar));
}

static void BM_LayoutSwap_3D(benchmark::State& state) {
  const int D0 = state.range(0);
  const int D1 = state.range(1);
  const int D2 = state.range(2);

  Tensor<Scalar, 3, ColMajor> A(D0, D1, D2);
  A.setRandom();

  for (auto _ : state) {
    Tensor<Scalar, 3, RowMajor> B = A.swap_layout();
    benchmark::DoNotOptimize(B.data());
    benchmark::ClobberMemory();
  }
  // 1 read (A) + 1 write (B).
  state.SetBytesProcessed(state.iterations() * 2ll * static_cast<int64_t>(D0) * D1 * D2 * sizeof(Scalar));
}

// Composing swap_layout with a coefficient-wise op forces evaluation through
// the executor and exercises any subsequent block consumers.
static void BM_LayoutSwap_Composed(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);

  Tensor<Scalar, 2, ColMajor> A(M, N);
  Tensor<Scalar, 2, ColMajor> B(M, N);
  A.setRandom();
  B.setRandom();

  for (auto _ : state) {
    Tensor<Scalar, 2, RowMajor> C = (A + B).swap_layout();
    benchmark::DoNotOptimize(C.data());
    benchmark::ClobberMemory();
  }
  // 2 reads (A, B) + 1 write (C).
  state.SetBytesProcessed(state.iterations() * 3ll * static_cast<int64_t>(M) * N * sizeof(Scalar));
}

// {n, n} and {n, n, n}: explicit because dims are repeated.
#define LAYOUT_SWAP_SIZES ->Args({64, 64})->Args({256, 256})->Args({1024, 1024})
#define LAYOUT_SWAP_3D_SIZES ->Args({32, 32, 32})->Args({64, 64, 64})->Args({128, 128, 128})

BENCHMARK(BM_LayoutSwap_2D) LAYOUT_SWAP_SIZES;
BENCHMARK(BM_LayoutSwap_3D) LAYOUT_SWAP_3D_SIZES;
BENCHMARK(BM_LayoutSwap_Composed) LAYOUT_SWAP_SIZES;
