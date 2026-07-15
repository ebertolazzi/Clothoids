// Benchmarks for Eigen Tensor morphing operations: reshape, slice, chip, pad, stride.
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#define EIGEN_USE_THREADS

#include <benchmark/benchmark.h>
#include <unsupported/Eigen/Tensor>
#include <unsupported/Eigen/ThreadPool>

using namespace Eigen;

typedef float Scalar;

// --- Reshape (zero-cost if no evaluation needed; force eval via assignment) ---
static void BM_Reshape(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);

  Tensor<Scalar, 2> A(M, N);
  A.setRandom();

  Eigen::array<Index, 1> new_shape = {M * N};

  for (auto _ : state) {
    Tensor<Scalar, 1> B = A.reshape(new_shape);
    benchmark::DoNotOptimize(B.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * N * sizeof(Scalar));
}

// --- Slice ---
static void BM_Slice(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);

  Tensor<Scalar, 2> A(M, N);
  A.setRandom();

  int sliceM = M / 2;
  int sliceN = N / 2;
  Eigen::array<Index, 2> offsets = {0, 0};
  Eigen::array<Index, 2> extents = {sliceM, sliceN};

  for (auto _ : state) {
    Tensor<Scalar, 2> B = A.slice(offsets, extents);
    benchmark::DoNotOptimize(B.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * sliceM * sliceN * sizeof(Scalar));
}

// --- Chip (extract a sub-tensor along one dimension) ---
static void BM_Chip(benchmark::State& state) {
  const int D0 = state.range(0);
  const int D1 = state.range(1);
  const int D2 = state.range(2);

  Tensor<Scalar, 3> A(D0, D1, D2);
  A.setRandom();

  for (auto _ : state) {
    Tensor<Scalar, 2> B = A.chip(0, 0);
    benchmark::DoNotOptimize(B.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * D1 * D2 * sizeof(Scalar));
}

// --- Pad ---
static void BM_Pad(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);
  const int padSize = state.range(2);

  Tensor<Scalar, 2> A(M, N);
  A.setRandom();

  Eigen::array<std::pair<int, int>, 2> paddings;
  paddings[0] = {padSize, padSize};
  paddings[1] = {padSize, padSize};

  for (auto _ : state) {
    Tensor<Scalar, 2> B = A.pad(paddings);
    benchmark::DoNotOptimize(B.data());
    benchmark::ClobberMemory();
  }
  int outM = M + 2 * padSize;
  int outN = N + 2 * padSize;
  state.SetBytesProcessed(state.iterations() * outM * outN * sizeof(Scalar));
}

// --- Stride ---
static void BM_Stride(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);
  const int stride = state.range(2);

  Tensor<Scalar, 2> A(M, N);
  A.setRandom();

  Eigen::array<Index, 2> strides_arr = {stride, stride};

  for (auto _ : state) {
    Tensor<Scalar, 2> B = A.stride(strides_arr);
    benchmark::DoNotOptimize(B.data());
    benchmark::ClobberMemory();
  }
  int outM = (M + stride - 1) / stride;
  int outN = (N + stride - 1) / stride;
  state.SetBytesProcessed(state.iterations() * outM * outN * sizeof(Scalar));
}

// --- ThreadPool variants ---

static void BM_Slice_ThreadPool(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);
  const int threads = state.range(2);

  Tensor<Scalar, 2> A(M, N);
  A.setRandom();

  int sliceM = M / 2;
  int sliceN = N / 2;
  Eigen::array<Index, 2> offsets = {0, 0};
  Eigen::array<Index, 2> extents = {sliceM, sliceN};

  ThreadPool tp(threads);
  ThreadPoolDevice dev(&tp, threads);

  Tensor<Scalar, 2> B(sliceM, sliceN);

  for (auto _ : state) {
    B.device(dev) = A.slice(offsets, extents);
    benchmark::DoNotOptimize(B.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * sliceM * sliceN * sizeof(Scalar));
  state.counters["threads"] = threads;
}

static void BM_Pad_ThreadPool(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);
  const int threads = state.range(2);

  Tensor<Scalar, 2> A(M, N);
  A.setRandom();

  Eigen::array<std::pair<int, int>, 2> paddings;
  paddings[0] = {4, 4};
  paddings[1] = {4, 4};

  int outM = M + 8;
  int outN = N + 8;

  ThreadPool tp(threads);
  ThreadPoolDevice dev(&tp, threads);

  Tensor<Scalar, 2> B(outM, outN);

  for (auto _ : state) {
    B.device(dev) = A.pad(paddings);
    benchmark::DoNotOptimize(B.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * outM * outN * sizeof(Scalar));
  state.counters["threads"] = threads;
}

// clang-format off
#define MORPH_SIZES \
  ->Args({256, 256})->Args({1024, 1024})

#define CHIP_SIZES \
  ->Args({32, 256, 256})->Args({64, 128, 128})->Args({8, 512, 512})

#define PAD_SIZES \
  ->Args({256, 256, 1})->Args({256, 256, 4})->Args({256, 256, 16}) \
  ->Args({1024, 1024, 1})->Args({1024, 1024, 4})->Args({1024, 1024, 16})

#define STRIDE_SIZES \
  ->Args({256, 256, 2})->Args({256, 256, 4}) \
  ->Args({1024, 1024, 2})->Args({1024, 1024, 4})

#define MORPH_THREADPOOL_SIZES \
  ->Args({256, 256, 1})->Args({256, 256, 2})->Args({256, 256, 4}) \
  ->Args({256, 256, 8})->Args({256, 256, 12})->Args({256, 256, 16}) \
  ->Args({1024, 1024, 1})->Args({1024, 1024, 2})->Args({1024, 1024, 4}) \
  ->Args({1024, 1024, 8})->Args({1024, 1024, 12})->Args({1024, 1024, 16})
// clang-format on

BENCHMARK(BM_Reshape) MORPH_SIZES;
BENCHMARK(BM_Slice) MORPH_SIZES;
BENCHMARK(BM_Chip) CHIP_SIZES;
BENCHMARK(BM_Pad) PAD_SIZES;
BENCHMARK(BM_Stride) STRIDE_SIZES;
BENCHMARK(BM_Slice_ThreadPool) MORPH_THREADPOOL_SIZES->UseRealTime();
BENCHMARK(BM_Pad_ThreadPool) MORPH_THREADPOOL_SIZES->UseRealTime();
