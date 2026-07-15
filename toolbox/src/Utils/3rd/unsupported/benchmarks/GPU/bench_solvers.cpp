// GPU solver benchmarks: gpu::LLT and gpu::LU compute + solve throughput.
//
// Measures factorization and solve performance for the host-matrix and
// DeviceMatrix code paths across a range of matrix sizes.
//
// For Nsight Systems profiling:
//   nsys profile --trace=cuda,nvtx ./bench_solvers
//
// For Nsight Compute kernel analysis:
//   ncu --set full -o profile ./bench_solvers --benchmark_filter=BM_GpuLLT_Compute/4096
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#include <benchmark/benchmark.h>

#include <Eigen/Cholesky>
#include <unsupported/Eigen/GPU>
#include <Eigen/LU>

using namespace Eigen;

#ifndef SCALAR
#define SCALAR double
#endif

using Scalar = SCALAR;
using Mat = Matrix<Scalar, Dynamic, Dynamic>;

// --------------------------------------------------------------------------
// Helpers
// --------------------------------------------------------------------------

static Mat make_spd(Index n) {
  Mat M = Mat::Random(n, n);
  return M.adjoint() * M + Mat::Identity(n, n) * static_cast<Scalar>(n);
}

// CUDA warm-up: ensure the GPU is initialized before timing.
static void cuda_warmup() {
  static bool done = false;
  if (!done) {
    void* p;
    EIGEN_CUDA_RUNTIME_CHECK(cudaMalloc(&p, 1));
    EIGEN_CUDA_RUNTIME_CHECK(cudaFree(p));
    done = true;
  }
}

// --------------------------------------------------------------------------
// GpuLLT benchmarks
// --------------------------------------------------------------------------

// Factorize from host matrix (includes H2D upload).
static void BM_GpuLLT_Compute_Host(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  Mat A = make_spd(n);
  gpu::LLT<Scalar> llt;

  for (auto _ : state) {
    llt.compute(A);
    if (llt.info() != Success) state.SkipWithError("factorization failed");
  }

  double flops = static_cast<double>(n) * static_cast<double>(n) * static_cast<double>(n) / 3.0;
  state.counters["GFLOPS"] =
      benchmark::Counter(flops, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
  state.counters["n"] = n;
}

// Factorize from DeviceMatrix (D2D copy path).
static void BM_GpuLLT_Compute_Device(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  Mat A = make_spd(n);
  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  gpu::LLT<Scalar> llt;

  for (auto _ : state) {
    llt.compute(d_A);
    if (llt.info() != Success) state.SkipWithError("factorization failed");
  }

  double flops = static_cast<double>(n) * static_cast<double>(n) * static_cast<double>(n) / 3.0;
  state.counters["GFLOPS"] =
      benchmark::Counter(flops, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
  state.counters["n"] = n;
}

// Factorize from DeviceMatrix (move path, no copy).
static void BM_GpuLLT_Compute_DeviceMove(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  Mat A = make_spd(n);
  gpu::LLT<Scalar> llt;

  for (auto _ : state) {
    auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
    llt.compute(std::move(d_A));
    if (llt.info() != Success) state.SkipWithError("factorization failed");
  }

  double flops = static_cast<double>(n) * static_cast<double>(n) * static_cast<double>(n) / 3.0;
  state.counters["GFLOPS"] =
      benchmark::Counter(flops, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
  state.counters["n"] = n;
}

// Solve from host matrix (H2D + potrs + D2H).
static void BM_GpuLLT_Solve_Host(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  const Index nrhs = state.range(1);
  Mat A = make_spd(n);
  Mat B = Mat::Random(n, nrhs);
  gpu::LLT<Scalar> llt(A);

  for (auto _ : state) {
    Mat X = llt.solve(B);
    benchmark::DoNotOptimize(X.data());
  }

  state.counters["n"] = n;
  state.counters["nrhs"] = nrhs;
}

// Solve from DeviceMatrix (D2D + potrs, async, toHost at end).
static void BM_GpuLLT_Solve_Device(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  const Index nrhs = state.range(1);
  Mat A = make_spd(n);
  Mat B = Mat::Random(n, nrhs);
  gpu::LLT<Scalar> llt(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  for (auto _ : state) {
    gpu::DeviceMatrix<Scalar> d_X = llt.solve(d_B);
    Mat X = d_X.toHost();
    benchmark::DoNotOptimize(X.data());
  }

  state.counters["n"] = n;
  state.counters["nrhs"] = nrhs;
}

// Solve staying entirely on device (no toHost — measures pure GPU time).
static void BM_GpuLLT_Solve_DeviceOnly(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  const Index nrhs = state.range(1);
  Mat A = make_spd(n);
  Mat B = Mat::Random(n, nrhs);
  gpu::LLT<Scalar> llt(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  for (auto _ : state) {
    gpu::DeviceMatrix<Scalar> d_X = llt.solve(d_B);
    // Force completion without D2H transfer.
    EIGEN_CUDA_RUNTIME_CHECK(cudaStreamSynchronize(llt.stream()));
    benchmark::DoNotOptimize(d_X.data());
  }

  state.counters["n"] = n;
  state.counters["nrhs"] = nrhs;
}

// --------------------------------------------------------------------------
// GpuLU benchmarks
// --------------------------------------------------------------------------

static void BM_GpuLU_Compute_Host(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  Mat A = Mat::Random(n, n);
  gpu::LU<Scalar> lu;

  for (auto _ : state) {
    lu.compute(A);
    if (lu.info() != Success) state.SkipWithError("factorization failed");
  }

  double flops = 2.0 / 3.0 * static_cast<double>(n) * static_cast<double>(n) * static_cast<double>(n);
  state.counters["GFLOPS"] =
      benchmark::Counter(flops, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
  state.counters["n"] = n;
}

static void BM_GpuLU_Compute_Device(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  Mat A = Mat::Random(n, n);
  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  gpu::LU<Scalar> lu;

  for (auto _ : state) {
    lu.compute(d_A);
    if (lu.info() != Success) state.SkipWithError("factorization failed");
  }

  double flops = 2.0 / 3.0 * static_cast<double>(n) * static_cast<double>(n) * static_cast<double>(n);
  state.counters["GFLOPS"] =
      benchmark::Counter(flops, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
  state.counters["n"] = n;
}

static void BM_GpuLU_Solve_Host(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  const Index nrhs = state.range(1);
  Mat A = Mat::Random(n, n);
  Mat B = Mat::Random(n, nrhs);
  gpu::LU<Scalar> lu(A);

  for (auto _ : state) {
    Mat X = lu.solve(B);
    benchmark::DoNotOptimize(X.data());
  }

  state.counters["n"] = n;
  state.counters["nrhs"] = nrhs;
}

static void BM_GpuLU_Solve_Device(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  const Index nrhs = state.range(1);
  Mat A = Mat::Random(n, n);
  Mat B = Mat::Random(n, nrhs);
  gpu::LU<Scalar> lu(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  for (auto _ : state) {
    gpu::DeviceMatrix<Scalar> d_X = lu.solve(d_B);
    Mat X = d_X.toHost();
    benchmark::DoNotOptimize(X.data());
  }

  state.counters["n"] = n;
  state.counters["nrhs"] = nrhs;
}

// --------------------------------------------------------------------------
// CPU baselines for comparison
// --------------------------------------------------------------------------

static void BM_CpuLLT_Compute(benchmark::State& state) {
  const Index n = state.range(0);
  Mat A = make_spd(n);
  LLT<Mat> llt;

  for (auto _ : state) {
    llt.compute(A);
    benchmark::DoNotOptimize(llt.matrixLLT().data());
  }

  double flops = static_cast<double>(n) * static_cast<double>(n) * static_cast<double>(n) / 3.0;
  state.counters["GFLOPS"] =
      benchmark::Counter(flops, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
  state.counters["n"] = n;
}

static void BM_CpuLU_Compute(benchmark::State& state) {
  const Index n = state.range(0);
  Mat A = Mat::Random(n, n);
  PartialPivLU<Mat> lu;

  for (auto _ : state) {
    lu.compute(A);
    benchmark::DoNotOptimize(lu.matrixLU().data());
  }

  double flops = 2.0 / 3.0 * static_cast<double>(n) * static_cast<double>(n) * static_cast<double>(n);
  state.counters["GFLOPS"] =
      benchmark::Counter(flops, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
  state.counters["n"] = n;
}

// --------------------------------------------------------------------------
// Registration
// --------------------------------------------------------------------------

// clang-format off
BENCHMARK(BM_GpuLLT_Compute_Host)->ArgsProduct({{64, 128, 256, 512, 1024, 2048, 4096}})->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_GpuLLT_Compute_Device)->ArgsProduct({{64, 128, 256, 512, 1024, 2048, 4096}})->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_GpuLLT_Compute_DeviceMove)->ArgsProduct({{64, 128, 256, 512, 1024, 2048, 4096}})->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_GpuLLT_Solve_Host)->ArgsProduct({{64, 256, 1024, 4096}, {1, 16}})->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_GpuLLT_Solve_Device)->ArgsProduct({{64, 256, 1024, 4096}, {1, 16}})->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_GpuLLT_Solve_DeviceOnly)->ArgsProduct({{64, 256, 1024, 4096}, {1, 16}})->Unit(benchmark::kMicrosecond);

BENCHMARK(BM_GpuLU_Compute_Host)->ArgsProduct({{64, 128, 256, 512, 1024, 2048, 4096}})->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_GpuLU_Compute_Device)->ArgsProduct({{64, 128, 256, 512, 1024, 2048, 4096}})->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_GpuLU_Solve_Host)->ArgsProduct({{64, 256, 1024, 4096}, {1, 16}})->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_GpuLU_Solve_Device)->ArgsProduct({{64, 256, 1024, 4096}, {1, 16}})->Unit(benchmark::kMicrosecond);

BENCHMARK(BM_CpuLLT_Compute)->ArgsProduct({{64, 128, 256, 512, 1024, 2048, 4096}})->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_CpuLU_Compute)->ArgsProduct({{64, 128, 256, 512, 1024, 2048, 4096}})->Unit(benchmark::kMicrosecond);
// clang-format on
