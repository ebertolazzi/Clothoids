// Benchmarks for Eigen Tensor coefficient-wise operations.
// Covers activation functions, normalization, and element-wise arithmetic.
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#define EIGEN_USE_THREADS

#include <benchmark/benchmark.h>
#include <unsupported/Eigen/Tensor>
#include <unsupported/Eigen/ThreadPool>

using namespace Eigen;

typedef float Scalar;

// Macro to define a benchmark for a unary tensor operation.
#define BENCH_TENSOR_UNARY(NAME, EXPR)                                        \
  static void BM_##NAME(benchmark::State& state) {                            \
    const int M = state.range(0);                                             \
    const int N = state.range(1);                                             \
    Tensor<Scalar, 2> a(M, N);                                                \
    a.setRandom();                                                            \
    Tensor<Scalar, 2> b(M, N);                                                \
    for (auto _ : state) {                                                    \
      b = EXPR;                                                               \
      benchmark::DoNotOptimize(b.data());                                     \
      benchmark::ClobberMemory();                                             \
    }                                                                         \
    state.SetBytesProcessed(state.iterations() * M * N * sizeof(Scalar) * 2); \
  }

// Macro for ThreadPool variant of a unary tensor operation.
#define BENCH_TENSOR_UNARY_THREADPOOL(NAME, EXPR)                             \
  static void BM_##NAME##_ThreadPool(benchmark::State& state) {               \
    const int M = state.range(0);                                             \
    const int N = state.range(1);                                             \
    const int threads = state.range(2);                                       \
    Tensor<Scalar, 2> a(M, N);                                                \
    a.setRandom();                                                            \
    Tensor<Scalar, 2> b(M, N);                                                \
    ThreadPool tp(threads);                                                   \
    ThreadPoolDevice dev(&tp, threads);                                       \
    for (auto _ : state) {                                                    \
      b.device(dev) = EXPR;                                                   \
      benchmark::DoNotOptimize(b.data());                                     \
      benchmark::ClobberMemory();                                             \
    }                                                                         \
    state.SetBytesProcessed(state.iterations() * M * N * sizeof(Scalar) * 2); \
    state.counters["threads"] = threads;                                      \
  }

BENCH_TENSOR_UNARY(Exp, a.exp())
BENCH_TENSOR_UNARY(Log, a.abs().log())
BENCH_TENSOR_UNARY(Tanh, a.tanh())
BENCH_TENSOR_UNARY(Sigmoid, a.sigmoid())
BENCH_TENSOR_UNARY(ReLU, a.cwiseMax(Scalar(0)))
BENCH_TENSOR_UNARY(Sqrt, a.abs().sqrt())

BENCH_TENSOR_UNARY_THREADPOOL(Exp, a.exp())
BENCH_TENSOR_UNARY_THREADPOOL(Tanh, a.tanh())
BENCH_TENSOR_UNARY_THREADPOOL(ReLU, a.cwiseMax(Scalar(0)))

// --- Element-wise binary operations ---
static void BM_Add(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);

  Tensor<Scalar, 2> a(M, N);
  Tensor<Scalar, 2> b(M, N);
  Tensor<Scalar, 2> c(M, N);
  a.setRandom();
  b.setRandom();

  for (auto _ : state) {
    c = a + b;
    benchmark::DoNotOptimize(c.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * N * sizeof(Scalar) * 3);
}

static void BM_Mul(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);

  Tensor<Scalar, 2> a(M, N);
  Tensor<Scalar, 2> b(M, N);
  Tensor<Scalar, 2> c(M, N);
  a.setRandom();
  b.setRandom();

  for (auto _ : state) {
    c = a * b;
    benchmark::DoNotOptimize(c.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * N * sizeof(Scalar) * 3);
}

// --- Fused multiply-add ---
static void BM_FMA(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);

  Tensor<Scalar, 2> a(M, N);
  Tensor<Scalar, 2> b(M, N);
  Tensor<Scalar, 2> c(M, N);
  Tensor<Scalar, 2> d(M, N);
  a.setRandom();
  b.setRandom();
  c.setRandom();

  for (auto _ : state) {
    d = a * b + c;
    benchmark::DoNotOptimize(d.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * N * sizeof(Scalar) * 4);
}

// --- ThreadPool binary operations ---
static void BM_Add_ThreadPool(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);
  const int threads = state.range(2);

  Tensor<Scalar, 2> a(M, N);
  Tensor<Scalar, 2> b(M, N);
  Tensor<Scalar, 2> c(M, N);
  a.setRandom();
  b.setRandom();

  ThreadPool tp(threads);
  ThreadPoolDevice dev(&tp, threads);

  for (auto _ : state) {
    c.device(dev) = a + b;
    benchmark::DoNotOptimize(c.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * N * sizeof(Scalar) * 3);
  state.counters["threads"] = threads;
}

static void BM_Mul_ThreadPool(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);
  const int threads = state.range(2);

  Tensor<Scalar, 2> a(M, N);
  Tensor<Scalar, 2> b(M, N);
  Tensor<Scalar, 2> c(M, N);
  a.setRandom();
  b.setRandom();

  ThreadPool tp(threads);
  ThreadPoolDevice dev(&tp, threads);

  for (auto _ : state) {
    c.device(dev) = a * b;
    benchmark::DoNotOptimize(c.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * N * sizeof(Scalar) * 3);
  state.counters["threads"] = threads;
}

static void BM_FMA_ThreadPool(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);
  const int threads = state.range(2);

  Tensor<Scalar, 2> a(M, N);
  Tensor<Scalar, 2> b(M, N);
  Tensor<Scalar, 2> c(M, N);
  Tensor<Scalar, 2> d(M, N);
  a.setRandom();
  b.setRandom();
  c.setRandom();

  ThreadPool tp(threads);
  ThreadPoolDevice dev(&tp, threads);

  for (auto _ : state) {
    d.device(dev) = a * b + c;
    benchmark::DoNotOptimize(d.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * N * sizeof(Scalar) * 4);
  state.counters["threads"] = threads;
}

// --- Rank-4 coefficient-wise (CNN feature maps) ---
static void BM_ReLU_Rank4(benchmark::State& state) {
  const int batch = state.range(0);
  const int C = state.range(1);
  const int H = state.range(2);

  Tensor<Scalar, 4> a(batch, C, H, H);
  Tensor<Scalar, 4> b(batch, C, H, H);
  a.setRandom();

  for (auto _ : state) {
    b = a.cwiseMax(Scalar(0));
    benchmark::DoNotOptimize(b.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * batch * C * H * H * sizeof(Scalar) * 2);
}

// clang-format off
#define CWISE_SIZES \
  ->Args({256, 256})->Args({1024, 1024})

#define CWISE_THREADPOOL_SIZES \
  ->Args({256, 256, 1})->Args({256, 256, 2})->Args({256, 256, 4}) \
  ->Args({256, 256, 8})->Args({256, 256, 12})->Args({256, 256, 16}) \
  ->Args({1024, 1024, 1})->Args({1024, 1024, 2})->Args({1024, 1024, 4}) \
  ->Args({1024, 1024, 8})->Args({1024, 1024, 12})->Args({1024, 1024, 16})

#define RANK4_SIZES \
  ->Args({32, 64, 16})->Args({8, 128, 32})->Args({1, 256, 64})
// clang-format on

BENCHMARK(BM_Exp) CWISE_SIZES;
BENCHMARK(BM_Log) CWISE_SIZES;
BENCHMARK(BM_Tanh) CWISE_SIZES;
BENCHMARK(BM_Sigmoid) CWISE_SIZES;
BENCHMARK(BM_ReLU) CWISE_SIZES;
BENCHMARK(BM_Sqrt) CWISE_SIZES;
BENCHMARK(BM_Add) CWISE_SIZES;
BENCHMARK(BM_Mul) CWISE_SIZES;
BENCHMARK(BM_FMA) CWISE_SIZES;
BENCHMARK(BM_ReLU_Rank4) RANK4_SIZES;
BENCHMARK(BM_Add_ThreadPool) CWISE_THREADPOOL_SIZES->UseRealTime();
BENCHMARK(BM_Mul_ThreadPool) CWISE_THREADPOOL_SIZES->UseRealTime();
BENCHMARK(BM_FMA_ThreadPool) CWISE_THREADPOOL_SIZES->UseRealTime();
BENCHMARK(BM_Exp_ThreadPool) CWISE_THREADPOOL_SIZES->UseRealTime();
BENCHMARK(BM_Tanh_ThreadPool) CWISE_THREADPOOL_SIZES->UseRealTime();
BENCHMARK(BM_ReLU_ThreadPool) CWISE_THREADPOOL_SIZES->UseRealTime();
