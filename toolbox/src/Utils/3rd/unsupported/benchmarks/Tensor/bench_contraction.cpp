// Benchmarks for Eigen Tensor contraction (generalized GEMM).
// Tests single-threaded (DefaultDevice) and multi-threaded (ThreadPoolDevice) variants.
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#define EIGEN_USE_THREADS

#include <benchmark/benchmark.h>
#include <unsupported/Eigen/Tensor>
#include <unsupported/Eigen/ThreadPool>

using namespace Eigen;

#ifndef SCALAR
#define SCALAR float
#endif

typedef SCALAR Scalar;

// --- DefaultDevice contraction (rank-2, equivalent to matrix multiply) ---
static void BM_Contraction(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);
  const int K = state.range(2);

  Tensor<Scalar, 2> A(M, K);
  Tensor<Scalar, 2> B(K, N);
  Tensor<Scalar, 2> C(M, N);
  A.setRandom();
  B.setRandom();

  using ContractDims = Tensor<Scalar, 2>::DimensionPair;
  Eigen::array<ContractDims, 1> contract_dims = {ContractDims(1, 0)};

  for (auto _ : state) {
    C = A.contract(B, contract_dims);
    benchmark::DoNotOptimize(C.data());
    benchmark::ClobberMemory();
  }
  state.counters["GFLOPS"] =
      benchmark::Counter(2.0 * M * N * K, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
}

// --- ThreadPoolDevice contraction ---
static void BM_Contraction_ThreadPool(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);
  const int K = state.range(2);
  const int threads = state.range(3);

  Tensor<Scalar, 2> A(M, K);
  Tensor<Scalar, 2> B(K, N);
  Tensor<Scalar, 2> C(M, N);
  A.setRandom();
  B.setRandom();

  ThreadPool tp(threads);
  ThreadPoolDevice dev(&tp, threads);

  using ContractDims = Tensor<Scalar, 2>::DimensionPair;
  Eigen::array<ContractDims, 1> contract_dims = {ContractDims(1, 0)};

  for (auto _ : state) {
    C.device(dev) = A.contract(B, contract_dims);
    benchmark::DoNotOptimize(C.data());
    benchmark::ClobberMemory();
  }
  state.counters["GFLOPS"] =
      benchmark::Counter(2.0 * M * N * K, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
  state.counters["threads"] = threads;
}

// --- Rank-3 batch contraction ---
// Contracts A(batch, M, K) with B(batch, K, N) over batch dim (0<->0)
// and K dim (2<->1), producing C(M, N). This sums over both the batch
// and inner dimensions: C(m, n) = sum_b sum_k A(b, m, k) * B(b, k, n).
static void BM_BatchContraction(benchmark::State& state) {
  const int batch = state.range(0);
  const int M = state.range(1);
  const int N = state.range(2);
  const int K = state.range(3);

  Tensor<Scalar, 3> A(batch, M, K);
  Tensor<Scalar, 3> B(batch, K, N);
  Tensor<Scalar, 2> C(M, N);
  A.setRandom();
  B.setRandom();

  using ContractDims = Tensor<Scalar, 3>::DimensionPair;
  Eigen::array<ContractDims, 2> contract_dims = {ContractDims(0, 0), ContractDims(2, 1)};

  for (auto _ : state) {
    C = A.contract(B, contract_dims);
    benchmark::DoNotOptimize(C.data());
    benchmark::ClobberMemory();
  }
  state.counters["GFLOPS"] = benchmark::Counter(2.0 * batch * M * N * K, benchmark::Counter::kIsIterationInvariantRate,
                                                benchmark::Counter::kIs1000);
}

// --- RowMajor contraction ---
static void BM_Contraction_RowMajor(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);
  const int K = state.range(2);

  Tensor<Scalar, 2, RowMajor> A(M, K);
  Tensor<Scalar, 2, RowMajor> B(K, N);
  Tensor<Scalar, 2, RowMajor> C(M, N);
  A.setRandom();
  B.setRandom();

  using ContractDims = Tensor<Scalar, 2, RowMajor>::DimensionPair;
  Eigen::array<ContractDims, 1> contract_dims = {ContractDims(1, 0)};

  for (auto _ : state) {
    C = A.contract(B, contract_dims);
    benchmark::DoNotOptimize(C.data());
    benchmark::ClobberMemory();
  }
  state.counters["GFLOPS"] =
      benchmark::Counter(2.0 * M * N * K, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
}

// clang-format off
#define CONTRACTION_SIZES \
  ->Args({32, 32, 32})->Args({64, 64, 64})->Args({128, 128, 128}) \
  ->Args({256, 256, 256})->Args({512, 512, 512})->Args({1024, 1024, 1024}) \
  ->Args({256, 256, 1024})->Args({1024, 64, 64})

#define CONTRACTION_THREADPOOL_SIZES \
  ->Args({64, 64, 64, 1})->Args({64, 64, 64, 2})->Args({64, 64, 64, 4}) \
  ->Args({64, 64, 64, 8})->Args({64, 64, 64, 16}) \
  ->Args({256, 256, 256, 1})->Args({256, 256, 256, 2})->Args({256, 256, 256, 4}) \
  ->Args({256, 256, 256, 8})->Args({256, 256, 256, 16}) \
  ->Args({512, 512, 512, 1})->Args({512, 512, 512, 2})->Args({512, 512, 512, 4}) \
  ->Args({512, 512, 512, 8})->Args({512, 512, 512, 16}) \
  ->Args({1024, 1024, 1024, 1})->Args({1024, 1024, 1024, 2})->Args({1024, 1024, 1024, 4}) \
  ->Args({1024, 1024, 1024, 8})->Args({1024, 1024, 1024, 16})

#define BATCH_SIZES \
  ->Args({1, 64, 64, 64})->Args({1, 256, 256, 256}) \
  ->Args({8, 64, 64, 64})->Args({8, 256, 256, 256}) \
  ->Args({32, 64, 64, 64})->Args({32, 256, 256, 256})
// clang-format on

BENCHMARK(BM_Contraction) CONTRACTION_SIZES;
BENCHMARK(BM_Contraction_RowMajor) CONTRACTION_SIZES;
BENCHMARK(BM_Contraction_ThreadPool) CONTRACTION_THREADPOOL_SIZES;
BENCHMARK(BM_BatchContraction) BATCH_SIZES;
