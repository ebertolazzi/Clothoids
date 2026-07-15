// Benchmarks for Eigen Tensor convolution (1D and 2D).
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#define EIGEN_USE_THREADS

#include <benchmark/benchmark.h>
#include <unsupported/Eigen/Tensor>
#include <unsupported/Eigen/ThreadPool>

using namespace Eigen;

typedef float Scalar;

// --- 1D convolution ---
static void BM_Convolve1D(benchmark::State& state) {
  const int input_size = state.range(0);
  const int kernel_size = state.range(1);

  Tensor<Scalar, 1> input(input_size);
  Tensor<Scalar, 1> kernel(kernel_size);
  input.setRandom();
  kernel.setRandom();

  Eigen::array<int, 1> dims = {0};

  for (auto _ : state) {
    Tensor<Scalar, 1> result = input.convolve(kernel, dims);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  double flops = 2.0 * (input_size - kernel_size + 1) * kernel_size;
  state.counters["GFLOPS"] =
      benchmark::Counter(flops, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
}

// --- 2D convolution ---
static void BM_Convolve2D(benchmark::State& state) {
  const int H = state.range(0);
  const int W = state.range(1);
  const int kH = state.range(2);
  const int kW = state.range(3);

  Tensor<Scalar, 2> input(H, W);
  Tensor<Scalar, 2> kernel(kH, kW);
  input.setRandom();
  kernel.setRandom();

  Eigen::array<int, 2> dims = {0, 1};

  for (auto _ : state) {
    Tensor<Scalar, 2> result = input.convolve(kernel, dims);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  double flops = 2.0 * (H - kH + 1) * (W - kW + 1) * kH * kW;
  state.counters["GFLOPS"] =
      benchmark::Counter(flops, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
}

// --- 2D convolution with channels (rank-3: C x H x W, convolve on H,W) ---
static void BM_Convolve2D_Channels(benchmark::State& state) {
  const int C = state.range(0);
  const int H = state.range(1);
  const int kH = state.range(2);

  Tensor<Scalar, 3> input(C, H, H);
  Tensor<Scalar, 2> kernel(kH, kH);
  input.setRandom();
  kernel.setRandom();

  Eigen::array<int, 2> dims = {1, 2};

  for (auto _ : state) {
    Tensor<Scalar, 3> result = input.convolve(kernel, dims);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  int outH = H - kH + 1;
  double flops = 2.0 * C * outH * outH * kH * kH;
  state.counters["GFLOPS"] =
      benchmark::Counter(flops, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
}

// --- 2D convolution with ThreadPool ---
static void BM_Convolve2D_ThreadPool(benchmark::State& state) {
  const int H = state.range(0);
  const int kH = state.range(1);
  const int threads = state.range(2);

  Tensor<Scalar, 2> input(H, H);
  Tensor<Scalar, 2> kernel(kH, kH);
  Tensor<Scalar, 2> result(H - kH + 1, H - kH + 1);
  input.setRandom();
  kernel.setRandom();

  ThreadPool tp(threads);
  ThreadPoolDevice dev(&tp, threads);

  Eigen::array<int, 2> dims = {0, 1};

  for (auto _ : state) {
    result.device(dev) = input.convolve(kernel, dims);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  int outH = H - kH + 1;
  double flops = 2.0 * outH * outH * kH * kH;
  state.counters["GFLOPS"] =
      benchmark::Counter(flops, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
  state.counters["threads"] = threads;
}

// {input, kernel}, {channels, hw, k}, {hw, k, threads}: pure Cartesian products.
#define CONV1D_SIZES ->ArgsProduct({{128, 512, 2048}, {3, 5, 11}})
#define CONV2D_CHANNEL_SIZES ->ArgsProduct({{3, 64, 128}, {16, 32, 56}, {3, 5}})
#define CONV2D_THREADPOOL_SIZES ->ArgsProduct({{64, 128, 224}, {3, 5}, {2, 4, 8}})

// {hw, hw, k, k}: explicit because hw and k are repeated.
// clang-format off
#define CONV2D_SIZES \
  ->Args({32, 32, 3, 3})->Args({32, 32, 5, 5})->Args({32, 32, 7, 7}) \
  ->Args({64, 64, 3, 3})->Args({64, 64, 5, 5})->Args({64, 64, 7, 7}) \
  ->Args({128, 128, 3, 3})->Args({128, 128, 5, 5})->Args({128, 128, 7, 7}) \
  ->Args({224, 224, 3, 3})->Args({224, 224, 5, 5})->Args({224, 224, 7, 7})
// clang-format on

BENCHMARK(BM_Convolve1D) CONV1D_SIZES;
BENCHMARK(BM_Convolve2D) CONV2D_SIZES;
BENCHMARK(BM_Convolve2D_Channels) CONV2D_CHANNEL_SIZES;
BENCHMARK(BM_Convolve2D_ThreadPool) CONV2D_THREADPOOL_SIZES;
