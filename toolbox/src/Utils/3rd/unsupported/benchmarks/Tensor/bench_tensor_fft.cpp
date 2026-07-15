// Benchmarks for Eigen Tensor FFT.
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#include <benchmark/benchmark.h>
#include <unsupported/Eigen/Tensor>

using namespace Eigen;

#ifndef SCALAR
#define SCALAR float
#endif

typedef SCALAR Scalar;

// 5*N*log2(N) is the conventional flop count for a length-N complex FFT.
static double FFTFlops(double n) { return 5.0 * n * std::log2(n); }

// --- 1D FFT (real input, complex output) ---
static void BM_TensorFFT_1D(benchmark::State& state) {
  const int N = state.range(0);
  Tensor<Scalar, 1> input(N);
  input.setRandom();
  Tensor<std::complex<Scalar>, 1> result(N);
  Eigen::array<int, 1> fft_dims = {0};
  for (auto _ : state) {
    result = input.template fft<BothParts, FFT_FORWARD>(fft_dims);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  state.counters["MFLOPS"] =
      benchmark::Counter(FFTFlops(N), benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
}

// --- 1D FFT (complex input, complex output) ---
static void BM_TensorFFT_1D_C2C(benchmark::State& state) {
  const int N = state.range(0);
  Tensor<std::complex<Scalar>, 1> input(N);
  input.setRandom();
  Tensor<std::complex<Scalar>, 1> result(N);
  Eigen::array<int, 1> fft_dims = {0};
  for (auto _ : state) {
    result = input.template fft<BothParts, FFT_FORWARD>(fft_dims);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  state.counters["MFLOPS"] =
      benchmark::Counter(FFTFlops(N), benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
}

// --- 1D inverse FFT ---
static void BM_TensorIFFT_1D(benchmark::State& state) {
  const int N = state.range(0);
  Tensor<std::complex<Scalar>, 1> input(N);
  input.setRandom();
  Tensor<std::complex<Scalar>, 1> result(N);
  Eigen::array<int, 1> fft_dims = {0};
  for (auto _ : state) {
    result = input.template fft<BothParts, FFT_REVERSE>(fft_dims);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  state.counters["MFLOPS"] =
      benchmark::Counter(FFTFlops(N), benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
}

// --- 2D FFT (square) ---
static void BM_TensorFFT_2D(benchmark::State& state) {
  const int N = state.range(0);
  Tensor<Scalar, 2> input(N, N);
  input.setRandom();
  Tensor<std::complex<Scalar>, 2> result(N, N);
  Eigen::array<int, 2> fft_dims = {0, 1};
  for (auto _ : state) {
    result = input.template fft<BothParts, FFT_FORWARD>(fft_dims);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  // A length-N 2D FFT is 2N length-N 1D FFTs.
  state.counters["MFLOPS"] = benchmark::Counter(2.0 * N * FFTFlops(N), benchmark::Counter::kIsIterationInvariantRate,
                                                benchmark::Counter::kIs1000);
}

// --- Batched 1D FFT along contiguous (innermost) dim ---
// Many parallel transforms; tests per-line setup overhead amortization.
static void BM_TensorFFT_Batched_Inner(benchmark::State& state) {
  const int N = state.range(0);
  const int batch = state.range(1);
  // ColMajor: dim 0 is innermost (stride 1).
  Tensor<std::complex<Scalar>, 2> input(N, batch);
  input.setRandom();
  Tensor<std::complex<Scalar>, 2> result(N, batch);
  Eigen::array<int, 1> fft_dims = {0};
  for (auto _ : state) {
    result = input.template fft<BothParts, FFT_FORWARD>(fft_dims);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  state.counters["MFLOPS"] = benchmark::Counter(batch * FFTFlops(N), benchmark::Counter::kIsIterationInvariantRate,
                                                benchmark::Counter::kIs1000);
}

// --- Batched 1D FFT along strided (outer) dim ---
// Stride != 1; exercises the gather/scatter copy path in evalToBuf.
static void BM_TensorFFT_Batched_Outer(benchmark::State& state) {
  const int N = state.range(0);
  const int batch = state.range(1);
  // ColMajor: dim 1 is outermost (stride = batch_inner).
  Tensor<std::complex<Scalar>, 2> input(batch, N);
  input.setRandom();
  Tensor<std::complex<Scalar>, 2> result(batch, N);
  Eigen::array<int, 1> fft_dims = {1};
  for (auto _ : state) {
    result = input.template fft<BothParts, FFT_FORWARD>(fft_dims);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  state.counters["MFLOPS"] = benchmark::Counter(batch * FFTFlops(N), benchmark::Counter::kIsIterationInvariantRate,
                                                benchmark::Counter::kIs1000);
}

// --- 1D Bluestein (non-power-of-2) ---
static void BM_TensorFFT_1D_Bluestein(benchmark::State& state) {
  const int N = state.range(0);
  Tensor<std::complex<Scalar>, 1> input(N);
  input.setRandom();
  Tensor<std::complex<Scalar>, 1> result(N);
  Eigen::array<int, 1> fft_dims = {0};
  for (auto _ : state) {
    result = input.template fft<BothParts, FFT_FORWARD>(fft_dims);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  state.counters["MFLOPS"] =
      benchmark::Counter(FFTFlops(N), benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
}

// --- Batched Bluestein (non-power-of-2, contiguous inner dim) ---
// Same b-sequence FFT shared by every line; tests the per-dim caching path.
static void BM_TensorFFT_Bluestein_Batched(benchmark::State& state) {
  const int N = state.range(0);
  const int batch = state.range(1);
  Tensor<std::complex<Scalar>, 2> input(N, batch);
  input.setRandom();
  Tensor<std::complex<Scalar>, 2> result(N, batch);
  Eigen::array<int, 1> fft_dims = {0};
  for (auto _ : state) {
    result = input.template fft<BothParts, FFT_FORWARD>(fft_dims);
    benchmark::DoNotOptimize(result.data());
    benchmark::ClobberMemory();
  }
  state.counters["MFLOPS"] = benchmark::Counter(batch * FFTFlops(N), benchmark::Counter::kIsIterationInvariantRate,
                                                benchmark::Counter::kIs1000);
}

// clang-format off
#define POW2_SIZES      ->Arg(64)->Arg(256)->Arg(1024)->Arg(4096)->Arg(16384)->Arg(65536)
#define POW2_2D         ->Arg(64)->Arg(256)->Arg(1024)
#define BLUESTEIN       ->Arg(100)->Arg(1000)->Arg(4099)
#define BATCH_SIZES     ->Args({64,64})->Args({256,64})->Args({1024,64})->Args({4096,64})
#define BLUESTEIN_BATCH ->Args({100,64})->Args({1000,64})->Args({4099,32})
// clang-format on

BENCHMARK(BM_TensorFFT_1D) POW2_SIZES;
BENCHMARK(BM_TensorFFT_1D_C2C) POW2_SIZES;
BENCHMARK(BM_TensorIFFT_1D) POW2_SIZES;
BENCHMARK(BM_TensorFFT_2D) POW2_2D;
BENCHMARK(BM_TensorFFT_Batched_Inner) BATCH_SIZES;
BENCHMARK(BM_TensorFFT_Batched_Outer) BATCH_SIZES;
BENCHMARK(BM_TensorFFT_1D_Bluestein) BLUESTEIN;
BENCHMARK(BM_TensorFFT_Bluestein_Batched) BLUESTEIN_BATCH;
