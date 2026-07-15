// Benchmarks for chained tensor expressions with ThreadPool.
// Tests realistic compound expressions spanning memory-bound to compute-bound.
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#define EIGEN_USE_THREADS

#include <benchmark/benchmark.h>
#include <unsupported/Eigen/Tensor>
#include <unsupported/Eigen/ThreadPool>

using namespace Eigen;

typedef float Scalar;

// --- Pure memory-bound baseline (dst = src) ---
static void BM_Copy_ThreadPool(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);
  const int threads = state.range(2);

  Tensor<Scalar, 2> src(M, N);
  Tensor<Scalar, 2> dst(M, N);
  src.setRandom();

  ThreadPool tp(threads);
  ThreadPoolDevice dev(&tp, threads);

  for (auto _ : state) {
    dst.device(dev) = src;
    benchmark::DoNotOptimize(dst.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * N * sizeof(Scalar) * 2);
  state.counters["threads"] = threads;
}

// --- Near-memory-bound: bias + ReLU ---
// Pattern: (x + bias.broadcast()).cwiseMax(0)
static void BM_BiasReLU_ThreadPool(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);
  const int threads = state.range(2);

  Tensor<Scalar, 2> x(M, N);
  Tensor<Scalar, 2> bias(1, N);
  Tensor<Scalar, 2> result(M, N);
  x.setRandom();
  bias.setRandom();

  ThreadPool tp(threads);
  ThreadPoolDevice dev(&tp, threads);

  Eigen::array<int, 2> bcast = {M, 1};

  for (auto _ : state) {
    result.device(dev) = (x + bias.broadcast(bcast)).cwiseMax(Scalar(0));
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * N * sizeof(Scalar) * 2);
  state.counters["threads"] = threads;
}

// --- Compute-bound: Horner polynomial ((a*x+b)*x+c)*x+d ---
static void BM_Polynomial_ThreadPool(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);
  const int threads = state.range(2);

  Tensor<Scalar, 2> x(M, N);
  Tensor<Scalar, 2> result(M, N);
  x.setRandom();

  ThreadPool tp(threads);
  ThreadPoolDevice dev(&tp, threads);

  const Scalar a = 0.5f, b = 1.2f, c = -0.3f, d = 0.7f;

  for (auto _ : state) {
    result.device(dev) = ((x.constant(a) * x + x.constant(b)) * x + x.constant(c)) * x + x.constant(d);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * N * sizeof(Scalar) * 2);
  state.counters["threads"] = threads;
}

// --- Compute-bound: exp (expensive transcendental) ---
static void BM_ExpNormalize_ThreadPool(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);
  const int threads = state.range(2);

  Tensor<Scalar, 2> x(M, N);
  Tensor<Scalar, 2> result(M, N);
  x.setRandom();

  ThreadPool tp(threads);
  ThreadPoolDevice dev(&tp, threads);

  for (auto _ : state) {
    result.device(dev) = x.exp();
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * N * sizeof(Scalar) * 2);
  state.counters["threads"] = threads;
}

// --- Batch normalization: gamma * (x - mean) / sqrt(var + eps) + beta ---
static void BM_BatchNorm_ThreadPool(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);
  const int threads = state.range(2);

  Tensor<Scalar, 2> x(M, N);
  Tensor<Scalar, 2> result(M, N);
  Tensor<Scalar, 2> gamma(1, N);
  Tensor<Scalar, 2> beta(1, N);
  Tensor<Scalar, 2> mean(1, N);
  Tensor<Scalar, 2> var(1, N);
  x.setRandom();
  gamma.setRandom();
  beta.setRandom();
  mean.setRandom();
  var.setRandom();
  var = var.abs() + var.constant(Scalar(0.1));

  ThreadPool tp(threads);
  ThreadPoolDevice dev(&tp, threads);

  Eigen::array<int, 2> bcast = {M, 1};
  const Scalar eps = 1e-5f;

  for (auto _ : state) {
    result.device(dev) =
        gamma.broadcast(bcast) * (x - mean.broadcast(bcast)) * (var.broadcast(bcast) + x.constant(eps)).rsqrt() +
        beta.broadcast(bcast);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * N * sizeof(Scalar) * 2);
  state.counters["threads"] = threads;
}

// clang-format off
#define CHAINED_SIZES \
  ->Args({256, 256, 1})->Args({256, 256, 2})->Args({256, 256, 4}) \
  ->Args({256, 256, 8})->Args({256, 256, 12})->Args({256, 256, 16}) \
  ->Args({1024, 1024, 1})->Args({1024, 1024, 2})->Args({1024, 1024, 4}) \
  ->Args({1024, 1024, 8})->Args({1024, 1024, 12})->Args({1024, 1024, 16}) \
  ->Args({4096, 4096, 1})->Args({4096, 4096, 2})->Args({4096, 4096, 4}) \
  ->Args({4096, 4096, 8})->Args({4096, 4096, 12})->Args({4096, 4096, 16})
// clang-format on

BENCHMARK(BM_Copy_ThreadPool) CHAINED_SIZES->UseRealTime();
BENCHMARK(BM_BiasReLU_ThreadPool) CHAINED_SIZES->UseRealTime();
BENCHMARK(BM_Polynomial_ThreadPool) CHAINED_SIZES->UseRealTime();
BENCHMARK(BM_ExpNormalize_ThreadPool) CHAINED_SIZES->UseRealTime();
BENCHMARK(BM_BatchNorm_ThreadPool) CHAINED_SIZES->UseRealTime();
