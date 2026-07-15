// Benchmarks for Eigen Tensor broadcasting.
// Tests broadcasting along various dimensions and ranks.
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#define EIGEN_USE_THREADS

#include <benchmark/benchmark.h>
#include <unsupported/Eigen/Tensor>
#include <unsupported/Eigen/ThreadPool>

using namespace Eigen;

typedef float Scalar;

// --- Broadcast row vector {1,N} -> {M,N} ---
static void BM_BroadcastRow(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);

  Tensor<Scalar, 2> row(1, N);
  Tensor<Scalar, 2> result(M, N);
  row.setRandom();

  Eigen::array<int, 2> bcast = {M, 1};

  for (auto _ : state) {
    result = row.broadcast(bcast);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * N * sizeof(Scalar));
}

// --- Broadcast col vector {M,1} -> {M,N} ---
static void BM_BroadcastCol(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);

  Tensor<Scalar, 2> col(M, 1);
  Tensor<Scalar, 2> result(M, N);
  col.setRandom();

  Eigen::array<int, 2> bcast = {1, N};

  for (auto _ : state) {
    result = col.broadcast(bcast);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * N * sizeof(Scalar));
}

// --- Broadcast + element-wise add (bias addition pattern) ---
static void BM_BroadcastAdd(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);

  Tensor<Scalar, 2> mat(M, N);
  Tensor<Scalar, 2> bias(1, N);
  Tensor<Scalar, 2> result(M, N);
  mat.setRandom();
  bias.setRandom();

  Eigen::array<int, 2> bcast = {M, 1};

  for (auto _ : state) {
    result = mat + bias.broadcast(bcast);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * N * sizeof(Scalar) * 2);
}

// --- Rank-4 broadcast (batch x channels x 1 x 1) -> (batch x channels x H x W) ---
static void BM_BroadcastRank4(benchmark::State& state) {
  const int batch = state.range(0);
  const int C = state.range(1);
  const int H = state.range(2);

  Tensor<Scalar, 4> bias(batch, C, 1, 1);
  Tensor<Scalar, 4> result(batch, C, H, H);
  bias.setRandom();

  Eigen::array<int, 4> bcast = {1, 1, H, H};

  for (auto _ : state) {
    result = bias.broadcast(bcast);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * batch * C * H * H * sizeof(Scalar));
}

// --- ThreadPool variants ---

static void BM_BroadcastRow_ThreadPool(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);
  const int threads = state.range(2);

  Tensor<Scalar, 2> row(1, N);
  Tensor<Scalar, 2> result(M, N);
  row.setRandom();

  ThreadPool tp(threads);
  ThreadPoolDevice dev(&tp, threads);

  Eigen::array<int, 2> bcast = {M, 1};

  for (auto _ : state) {
    result.device(dev) = row.broadcast(bcast);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * N * sizeof(Scalar));
  state.counters["threads"] = threads;
}

static void BM_BroadcastAdd_ThreadPool(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);
  const int threads = state.range(2);

  Tensor<Scalar, 2> mat(M, N);
  Tensor<Scalar, 2> bias(1, N);
  Tensor<Scalar, 2> result(M, N);
  mat.setRandom();
  bias.setRandom();

  ThreadPool tp(threads);
  ThreadPoolDevice dev(&tp, threads);

  Eigen::array<int, 2> bcast = {M, 1};

  for (auto _ : state) {
    result.device(dev) = mat + bias.broadcast(bcast);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * N * sizeof(Scalar) * 2);
  state.counters["threads"] = threads;
}

// {m, n} and {batch, c, h}: pure Cartesian products.
#define BROADCAST_SIZES ->ArgsProduct({{64, 256, 1024}, {64, 256, 1024}})
#define BROADCAST_RANK4_SIZES ->ArgsProduct({{1, 8}, {64, 256}, {16, 32}})

// {size, size, threads}: explicit because size is repeated.
// clang-format off
#define BROADCAST_THREADPOOL_SIZES \
  ->Args({256, 256, 1})->Args({256, 256, 2})->Args({256, 256, 4}) \
  ->Args({256, 256, 8})->Args({256, 256, 12})->Args({256, 256, 16}) \
  ->Args({1024, 1024, 1})->Args({1024, 1024, 2})->Args({1024, 1024, 4}) \
  ->Args({1024, 1024, 8})->Args({1024, 1024, 12})->Args({1024, 1024, 16})
// clang-format on

BENCHMARK(BM_BroadcastRow) BROADCAST_SIZES;
BENCHMARK(BM_BroadcastCol) BROADCAST_SIZES;
BENCHMARK(BM_BroadcastAdd) BROADCAST_SIZES;
BENCHMARK(BM_BroadcastRank4) BROADCAST_RANK4_SIZES;
BENCHMARK(BM_BroadcastRow_ThreadPool) BROADCAST_THREADPOOL_SIZES->UseRealTime();
BENCHMARK(BM_BroadcastAdd_ThreadPool) BROADCAST_THREADPOOL_SIZES->UseRealTime();
