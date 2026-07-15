// GEMM dispatch benchmark: measure DeviceMatrix GEMM throughput across sizes.
//
// Compares cublasLtMatmul with plan cache against a raw cublasGemmEx call
// (no descriptor overhead) to verify no regression.
//
// Usage:
//   cmake --build build-bench-gpu --target bench_gemm
//   ./build-bench-gpu/bench_gemm
//
// Profiling:
//   nsys profile --trace=cuda ./build-bench-gpu/bench_gemm
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#include <benchmark/benchmark.h>

// EIGEN_USE_GPU is set by the CMake target (eigen_add_gpu_benchmark).
#include <unsupported/Eigen/GPU>

using namespace Eigen;

#ifndef SCALAR
#define SCALAR double
#endif

using Scalar = SCALAR;
using Mat = Matrix<Scalar, Dynamic, Dynamic>;

static void cuda_warmup() {
  static bool done = false;
  if (!done) {
    void* p;
    cudaMalloc(&p, 256);
    cudaFree(p);
    // Force context creation and JIT.
    gpu::Context ctx;
    Mat A = Mat::Random(64, 64);
    Mat B = Mat::Random(64, 64);
    auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A, ctx.stream());
    auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B, ctx.stream());
    gpu::DeviceMatrix<Scalar> d_C;
    d_C.device(ctx) = d_A * d_B;
    if (cudaDeviceSynchronize() != cudaSuccess) abort();
    done = true;
  }
}

// --------------------------------------------------------------------------
// DeviceMatrix GEMM (uses cublasLtMatmul with plan cache)
// --------------------------------------------------------------------------

static void BM_DeviceMatrix_Gemm(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);

  gpu::Context ctx;
  Mat hostA = Mat::Random(n, n);
  Mat hostB = Mat::Random(n, n);
  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(hostA, ctx.stream());
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(hostB, ctx.stream());
  gpu::DeviceMatrix<Scalar> d_C;

  // Warmup: run a few GEMMs to stabilize clocks and populate plan cache.
  for (int i = 0; i < 5; ++i) {
    d_C.device(ctx) = d_A * d_B;
  }
  if (cudaDeviceSynchronize() != cudaSuccess) abort();

  for (auto _ : state) {
    d_C.device(ctx) = d_A * d_B;
    if (cudaDeviceSynchronize() != cudaSuccess) abort();
  }

  double flops = 2.0 * n * n * n;
  state.counters["GFLOPS"] =
      benchmark::Counter(flops, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
  state.counters["n"] = n;
}

// --------------------------------------------------------------------------
// Raw cublasGemmEx (direct call, no descriptor overhead)
// --------------------------------------------------------------------------

static void BM_Raw_CublasGemmEx(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);

  gpu::Context ctx;
  Mat hostA = Mat::Random(n, n);
  Mat hostB = Mat::Random(n, n);
  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(hostA, ctx.stream());
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(hostB, ctx.stream());
  gpu::DeviceMatrix<Scalar> d_C(n, n);

  constexpr cudaDataType_t dtype = gpu::internal::cuda_data_type<Scalar>::value;
  constexpr cublasComputeType_t compute = gpu::internal::cuda_compute_type<Scalar>::value;
  Scalar alpha = Scalar(1);
  Scalar beta = Scalar(0);
  const int ni = static_cast<int>(n);

  constexpr cublasGemmAlgo_t algo = gpu::internal::cuda_gemm_algo();

  // Warmup.
  for (int i = 0; i < 5; ++i) {
    cublasStatus_t s = cublasGemmEx(ctx.cublasHandle(), CUBLAS_OP_N, CUBLAS_OP_N, ni, ni, ni, &alpha, d_A.data(), dtype,
                                    ni, d_B.data(), dtype, ni, &beta, d_C.data(), dtype, ni, compute, algo);
    if (s != CUBLAS_STATUS_SUCCESS) {
      state.SkipWithError("cublasGemmEx failed");
      return;
    }
  }
  if (cudaDeviceSynchronize() != cudaSuccess) abort();

  for (auto _ : state) {
    cublasStatus_t s = cublasGemmEx(ctx.cublasHandle(), CUBLAS_OP_N, CUBLAS_OP_N, ni, ni, ni, &alpha, d_A.data(), dtype,
                                    ni, d_B.data(), dtype, ni, &beta, d_C.data(), dtype, ni, compute, algo);
    if (s != CUBLAS_STATUS_SUCCESS) {
      state.SkipWithError("cublasGemmEx failed");
      return;
    }
    if (cudaDeviceSynchronize() != cudaSuccess) abort();
  }

  double flops = 2.0 * n * n * n;
  state.counters["GFLOPS"] =
      benchmark::Counter(flops, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
  state.counters["n"] = n;
}

// --------------------------------------------------------------------------
// DeviceMatrix GEMM with transpose: C = A^T * B
// --------------------------------------------------------------------------

static void BM_DeviceMatrix_Gemm_TransA(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);

  gpu::Context ctx;
  Mat hostA = Mat::Random(n, n);
  Mat hostB = Mat::Random(n, n);
  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(hostA, ctx.stream());
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(hostB, ctx.stream());
  gpu::DeviceMatrix<Scalar> d_C;

  for (int i = 0; i < 5; ++i) {
    d_C.device(ctx) = d_A.transpose() * d_B;
  }
  if (cudaDeviceSynchronize() != cudaSuccess) abort();

  for (auto _ : state) {
    d_C.device(ctx) = d_A.transpose() * d_B;
    if (cudaDeviceSynchronize() != cudaSuccess) abort();
  }

  double flops = 2.0 * n * n * n;
  state.counters["GFLOPS"] =
      benchmark::Counter(flops, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::kIs1000);
  state.counters["n"] = n;
}

// Square GEMM: range of sizes from small (where descriptor overhead matters)
// to large (where compute dominates).
BENCHMARK(BM_DeviceMatrix_Gemm)
    ->ArgsProduct({{16, 32, 64, 128, 256, 512, 1024, 2048, 4096}})
    ->Unit(benchmark::kMicrosecond);

BENCHMARK(BM_Raw_CublasGemmEx)
    ->ArgsProduct({{16, 32, 64, 128, 256, 512, 1024, 2048, 4096}})
    ->Unit(benchmark::kMicrosecond);

BENCHMARK(BM_DeviceMatrix_Gemm_TransA)
    ->ArgsProduct({{16, 32, 64, 128, 256, 512, 1024, 2048, 4096}})
    ->Unit(benchmark::kMicrosecond);
