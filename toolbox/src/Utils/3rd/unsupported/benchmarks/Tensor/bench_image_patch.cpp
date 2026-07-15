// Benchmarks for Eigen TensorImagePatch extraction.
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#define EIGEN_USE_THREADS

#include <benchmark/benchmark.h>
#include <unsupported/Eigen/Tensor>
#include <unsupported/Eigen/ThreadPool>

using namespace Eigen;

typedef float Scalar;

// --- Basic image patch extraction with PADDING_VALID ---
static void BM_ImagePatch_Valid(benchmark::State& state) {
  const int C = state.range(0);
  const int H = state.range(1);
  const int W = state.range(2);
  const int kH = state.range(3);
  const int kW = state.range(4);

  Tensor<Scalar, 4> input(C, H, W, 1);
  input.setRandom();
  const int outH = H - kH + 1;
  const int outW = W - kW + 1;
  Tensor<Scalar, 5> result(C, kH, kW, outH * outW, 1);

  for (auto _ : state) {
    result = input.extract_image_patches(kH, kW, 1, 1, 1, 1, PADDING_VALID);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }

  const double bytes = static_cast<double>(C) * outH * outW * kH * kW * sizeof(Scalar);
  state.SetBytesProcessed(state.iterations() * bytes);
}

// --- Basic image patch extraction with PADDING_SAME ---
static void BM_ImagePatch_Same(benchmark::State& state) {
  const int C = state.range(0);
  const int H = state.range(1);
  const int W = state.range(2);
  const int kH = state.range(3);
  const int kW = state.range(4);

  Tensor<Scalar, 4> input(C, H, W, 1);
  input.setRandom();
  const int outH = H;
  const int outW = W;
  Tensor<Scalar, 5> result(C, kH, kW, outH * outW, 1);

  for (auto _ : state) {
    result = input.extract_image_patches(kH, kW, 1, 1, 1, 1, PADDING_SAME);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }

  const double bytes = static_cast<double>(C) * H * W * kH * kW * sizeof(Scalar);
  state.SetBytesProcessed(state.iterations() * bytes);
}

// --- Image patch with strides (simulates strided convolution) ---
static void BM_ImagePatch_Strided(benchmark::State& state) {
  const int C = state.range(0);
  const int H = state.range(1);
  const int kH = state.range(2);
  const int stride = state.range(3);

  Tensor<Scalar, 4> input(C, H, H, 1);
  input.setRandom();
  const int outH = (H + stride - 1) / stride;
  Tensor<Scalar, 5> result(C, kH, kH, outH * outH, 1);

  for (auto _ : state) {
    result = input.extract_image_patches(kH, kH, stride, stride, 1, 1, PADDING_SAME);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }

  const double bytes = static_cast<double>(C) * outH * outH * kH * kH * sizeof(Scalar);
  state.SetBytesProcessed(state.iterations() * bytes);
}

// --- Image patch with dilation (atrous/dilated convolution) ---
static void BM_ImagePatch_Dilated(benchmark::State& state) {
  const int C = state.range(0);
  const int H = state.range(1);
  const int kH = state.range(2);
  const int dilation = state.range(3);

  Tensor<Scalar, 4> input(C, H, H, 1);
  input.setRandom();
  const int outH = H;
  Tensor<Scalar, 5> result(C, kH, kH, outH * outH, 1);

  for (auto _ : state) {
    result = input.extract_image_patches(kH, kH, 1, 1, dilation, dilation, PADDING_SAME);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }

  const double bytes = static_cast<double>(C) * H * H * kH * kH * sizeof(Scalar);
  state.SetBytesProcessed(state.iterations() * bytes);
}

// --- Image patch with explicit padding ---
static void BM_ImagePatch_ExplicitPadding(benchmark::State& state) {
  const int C = state.range(0);
  const int H = state.range(1);
  const int W = state.range(2);
  const int kH = state.range(3);

  const int pad = kH / 2;

  Tensor<Scalar, 4> input(C, H, W, 1);
  input.setRandom();
  const int outH = H + 2 * pad - kH + 1;
  const int outW = W + 2 * pad - kH + 1;
  Tensor<Scalar, 5> result(C, kH, kH, outH * outW, 1);

  for (auto _ : state) {
    result = input.extract_image_patches(kH, kH, 1, 1, 1, 1, 1, 1, pad, pad, pad, pad, Scalar(0));
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }

  const double bytes = static_cast<double>(C) * outH * outW * kH * kH * sizeof(Scalar);
  state.SetBytesProcessed(state.iterations() * bytes);
}

// --- Batched image patch (multiple images) ---
static void BM_ImagePatch_Batched(benchmark::State& state) {
  const int C = state.range(0);
  const int H = state.range(1);
  const int kH = state.range(2);
  const int batch = state.range(3);

  Tensor<Scalar, 4> input(C, H, H, batch);
  input.setRandom();
  const int outH = H;
  Tensor<Scalar, 5> result(C, kH, kH, outH * outH, batch);

  for (auto _ : state) {
    result = input.extract_image_patches(kH, kH, 1, 1, 1, 1, PADDING_SAME);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }

  const double bytes = static_cast<double>(batch) * C * H * H * kH * kH * sizeof(Scalar);
  state.SetBytesProcessed(state.iterations() * bytes);
}

// --- ImageNet-style configurations (realistic CNN layer sizes) ---
static void BM_ImagePatch_ImageNet(benchmark::State& state) {
  const int C = state.range(0);
  const int H = state.range(1);
  const int kH = state.range(2);
  const int stride = state.range(3);

  Tensor<Scalar, 4> input(C, H, H, 1);
  input.setRandom();
  const int outH = (H + stride - 1) / stride;
  Tensor<Scalar, 5> result(C, kH, kH, outH * outH, 1);

  for (auto _ : state) {
    result = input.extract_image_patches(kH, kH, stride, stride, 1, 1, PADDING_SAME);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }

  const double bytes = static_cast<double>(C) * outH * outH * kH * kH * sizeof(Scalar);
  state.SetBytesProcessed(state.iterations() * bytes);
}

// --- ThreadPool variant ---
static void BM_ImagePatch_ThreadPool(benchmark::State& state) {
  const int C = state.range(0);
  const int H = state.range(1);
  const int kH = state.range(2);
  const int threads = state.range(3);

  Tensor<Scalar, 4> input(C, H, H, 1);
  input.setRandom();

  ThreadPool tp(threads);
  ThreadPoolDevice dev(&tp, threads);

  const int outH = H;
  Tensor<Scalar, 5> result(C, kH, kH, outH * outH, 1);

  for (auto _ : state) {
    result.device(dev) = input.extract_image_patches(kH, kH, 1, 1, 1, 1, PADDING_SAME);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }

  const double bytes = static_cast<double>(C) * H * H * kH * kH * sizeof(Scalar);
  state.SetBytesProcessed(state.iterations() * bytes);
  state.counters["threads"] = threads;
}

// --- Size configurations ---

// channels, H, W, kH, kW (H==W and kH==kW); explicit because of duplicated dims.
// clang-format off
#define PATCH_SIZES \
  ->Args({3, 32, 32, 3, 3})->Args({3, 32, 32, 5, 5})->Args({3, 32, 32, 7, 7}) \
  ->Args({3, 64, 64, 3, 3})->Args({3, 64, 64, 5, 5})->Args({3, 64, 64, 7, 7}) \
  ->Args({3, 128, 128, 3, 3})->Args({3, 128, 128, 5, 5})->Args({3, 128, 128, 7, 7}) \
  ->Args({32, 32, 32, 3, 3})->Args({32, 32, 32, 5, 5})->Args({32, 32, 32, 7, 7}) \
  ->Args({32, 64, 64, 3, 3})->Args({32, 64, 64, 5, 5})->Args({32, 64, 64, 7, 7}) \
  ->Args({32, 128, 128, 3, 3})->Args({32, 128, 128, 5, 5})->Args({32, 128, 128, 7, 7}) \
  ->Args({64, 32, 32, 3, 3})->Args({64, 32, 32, 5, 5})->Args({64, 32, 32, 7, 7}) \
  ->Args({64, 64, 64, 3, 3})->Args({64, 64, 64, 5, 5})->Args({64, 64, 64, 7, 7}) \
  ->Args({64, 128, 128, 3, 3})->Args({64, 128, 128, 5, 5})->Args({64, 128, 128, 7, 7})

// channels, H, W, kH (H==W); explicit because of duplicated H/W dim.
#define EXPLICIT_PADDING_SIZES \
  ->Args({3, 32, 32, 3})->Args({3, 32, 32, 5})->Args({3, 64, 64, 3})->Args({3, 64, 64, 5}) \
  ->Args({3, 128, 128, 3})->Args({3, 128, 128, 5})->Args({64, 32, 32, 3})->Args({64, 32, 32, 5}) \
  ->Args({64, 64, 64, 3})->Args({64, 64, 64, 5})->Args({64, 128, 128, 3})->Args({64, 128, 128, 5})

// {channels, spatial, kernel, stride/dilation/threads/batch}: pure Cartesian products.
#define STRIDED_SIZES ->ArgsProduct({{3, 64}, {56, 112, 224}, {3, 5}, {1, 2}})
#define DILATED_SIZES ->ArgsProduct({{3, 64}, {32, 64}, {3, 5}, {2, 4}})
#define BATCHED_SIZES ->ArgsProduct({{3, 64}, {32, 56}, {3, 5}, {4, 16, 32}})
#define THREAD_POOL_SIZES ->ArgsProduct({{64, 128}, {56, 112}, {3, 5}, {2, 4, 8}})

// Realistic CNN layer configurations: channels, spatial_size, kernel, stride.
// AlexNet conv1; VGG, VGG deeper x2; ResNet, ResNet downsample, ResNet deeper x2;
// MobileNet depthwise; Inception 1x1 (degenerate patch).
#define IMAGENET_SIZES \
  ->Args({3, 227, 11, 4}) \
  ->Args({64, 224, 3, 1})->Args({128, 112, 3, 1})->Args({256, 56, 3, 1}) \
  ->Args({64, 56, 3, 1})->Args({128, 56, 3, 2})->Args({256, 28, 3, 1})->Args({512, 14, 3, 1}) \
  ->Args({32, 112, 3, 1})->Args({192, 28, 1, 1})
// clang-format on

BENCHMARK(BM_ImagePatch_Valid) PATCH_SIZES;
BENCHMARK(BM_ImagePatch_Same) PATCH_SIZES;
BENCHMARK(BM_ImagePatch_Strided) STRIDED_SIZES;
BENCHMARK(BM_ImagePatch_Dilated) DILATED_SIZES;
BENCHMARK(BM_ImagePatch_ExplicitPadding) EXPLICIT_PADDING_SIZES;
BENCHMARK(BM_ImagePatch_Batched) BATCHED_SIZES;
BENCHMARK(BM_ImagePatch_ImageNet) IMAGENET_SIZES;
BENCHMARK(BM_ImagePatch_ThreadPool) THREAD_POOL_SIZES;
