// GPU batching benchmarks: multi-stream concurrency for many small solves.
//
// Each gpu::LLT/gpu::LU owns its own CUDA stream. This benchmark measures how
// well multiple solver instances overlap on the GPU, which is critical for
// workloads like robotics (many small systems) and SLAM (batched poses).
//
// Compares:
//   1. Sequential: one solver handles all systems one by one
//   2. Batched: N solvers on N streams, all launched before any sync
//   3. CPU baseline: Eigen LLT on host
//
// For Nsight Systems: batched mode should show overlapping kernels on
// different streams in the timeline view.
//
//   nsys profile --trace=cuda ./bench_batching
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#include <benchmark/benchmark.h>

#include <Eigen/Cholesky>
#include <unsupported/Eigen/GPU>

#include <memory>
#include <vector>

using namespace Eigen;

#ifndef SCALAR
#define SCALAR double
#endif

using Scalar = SCALAR;
using Mat = Matrix<Scalar, Dynamic, Dynamic>;

static Mat make_spd(Index n) {
  Mat M = Mat::Random(n, n);
  return M.adjoint() * M + Mat::Identity(n, n) * static_cast<Scalar>(n);
}

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
// Sequential: one solver, N systems solved one after another
// --------------------------------------------------------------------------

static void BM_Batch_Sequential(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  const int batch_size = static_cast<int>(state.range(1));

  // Pre-generate all SPD matrices and RHS vectors.
  std::vector<Mat> As(batch_size);
  std::vector<Mat> Bs(batch_size);
  for (int i = 0; i < batch_size; ++i) {
    As[i] = make_spd(n);
    Bs[i] = Mat::Random(n, 1);
  }

  gpu::LLT<Scalar> llt;

  for (auto _ : state) {
    std::vector<Mat> results(batch_size);
    for (int i = 0; i < batch_size; ++i) {
      llt.compute(As[i]);
      results[i] = llt.solve(Bs[i]);
    }
    benchmark::DoNotOptimize(results.back().data());
  }

  state.counters["n"] = n;
  state.counters["batch"] = batch_size;
  state.counters["total_solves"] = batch_size;
}

// --------------------------------------------------------------------------
// Sequential with DeviceMatrix (avoid re-upload of A each iteration)
// --------------------------------------------------------------------------

static void BM_Batch_Sequential_Device(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  const int batch_size = static_cast<int>(state.range(1));

  std::vector<Mat> As(batch_size);
  std::vector<Mat> Bs(batch_size);
  std::vector<gpu::DeviceMatrix<Scalar>> d_As(batch_size);
  std::vector<gpu::DeviceMatrix<Scalar>> d_Bs(batch_size);
  for (int i = 0; i < batch_size; ++i) {
    As[i] = make_spd(n);
    Bs[i] = Mat::Random(n, 1);
    d_As[i] = gpu::DeviceMatrix<Scalar>::fromHost(As[i]);
    d_Bs[i] = gpu::DeviceMatrix<Scalar>::fromHost(Bs[i]);
  }

  gpu::LLT<Scalar> llt;

  for (auto _ : state) {
    std::vector<Mat> results(batch_size);
    for (int i = 0; i < batch_size; ++i) {
      llt.compute(d_As[i]);
      gpu::DeviceMatrix<Scalar> d_X = llt.solve(d_Bs[i]);
      results[i] = d_X.toHost();
    }
    benchmark::DoNotOptimize(results.back().data());
  }

  state.counters["n"] = n;
  state.counters["batch"] = batch_size;
  state.counters["total_solves"] = batch_size;
}

// --------------------------------------------------------------------------
// Batched: N solvers on N streams, overlapping execution
// --------------------------------------------------------------------------

static void BM_Batch_MultiStream(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  const int batch_size = static_cast<int>(state.range(1));

  std::vector<Mat> As(batch_size);
  std::vector<Mat> Bs(batch_size);
  std::vector<gpu::DeviceMatrix<Scalar>> d_As(batch_size);
  std::vector<gpu::DeviceMatrix<Scalar>> d_Bs(batch_size);
  for (int i = 0; i < batch_size; ++i) {
    As[i] = make_spd(n);
    Bs[i] = Mat::Random(n, 1);
    d_As[i] = gpu::DeviceMatrix<Scalar>::fromHost(As[i]);
    d_Bs[i] = gpu::DeviceMatrix<Scalar>::fromHost(Bs[i]);
  }

  // N solvers = N independent CUDA streams.
  std::vector<std::unique_ptr<gpu::LLT<Scalar>>> solvers(batch_size);
  for (int i = 0; i < batch_size; ++i) {
    solvers[i] = std::make_unique<gpu::LLT<Scalar>>();
  }

  for (auto _ : state) {
    // Phase 1: launch all factorizations (async, different streams).
    for (int i = 0; i < batch_size; ++i) {
      solvers[i]->compute(d_As[i]);
    }

    // Phase 2: launch all solves (async, different streams).
    std::vector<gpu::DeviceMatrix<Scalar>> d_Xs(batch_size);
    for (int i = 0; i < batch_size; ++i) {
      d_Xs[i] = solvers[i]->solve(d_Bs[i]);
    }

    // Phase 3: download all results.
    std::vector<Mat> results(batch_size);
    for (int i = 0; i < batch_size; ++i) {
      results[i] = d_Xs[i].toHost();
    }
    benchmark::DoNotOptimize(results.back().data());
  }

  state.counters["n"] = n;
  state.counters["batch"] = batch_size;
  state.counters["streams"] = batch_size;
  state.counters["total_solves"] = batch_size;
}

// --------------------------------------------------------------------------
// Batched with async download (overlap D2H with computation)
// --------------------------------------------------------------------------

static void BM_Batch_MultiStream_AsyncDownload(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  const int batch_size = static_cast<int>(state.range(1));

  std::vector<Mat> As(batch_size);
  std::vector<Mat> Bs(batch_size);
  std::vector<gpu::DeviceMatrix<Scalar>> d_As(batch_size);
  std::vector<gpu::DeviceMatrix<Scalar>> d_Bs(batch_size);
  for (int i = 0; i < batch_size; ++i) {
    As[i] = make_spd(n);
    Bs[i] = Mat::Random(n, 1);
    d_As[i] = gpu::DeviceMatrix<Scalar>::fromHost(As[i]);
    d_Bs[i] = gpu::DeviceMatrix<Scalar>::fromHost(Bs[i]);
  }

  std::vector<std::unique_ptr<gpu::LLT<Scalar>>> solvers(batch_size);
  for (int i = 0; i < batch_size; ++i) {
    solvers[i] = std::make_unique<gpu::LLT<Scalar>>();
  }

  for (auto _ : state) {
    // Launch all compute + solve.
    std::vector<gpu::DeviceMatrix<Scalar>> d_Xs(batch_size);
    for (int i = 0; i < batch_size; ++i) {
      solvers[i]->compute(d_As[i]);
      d_Xs[i] = solvers[i]->solve(d_Bs[i]);
    }

    // Enqueue all async downloads.
    std::vector<gpu::HostTransfer<Scalar>> transfers;
    transfers.reserve(batch_size);
    for (int i = 0; i < batch_size; ++i) {
      transfers.push_back(d_Xs[i].toHostAsync());
    }

    // Collect all results.
    for (int i = 0; i < batch_size; ++i) {
      benchmark::DoNotOptimize(transfers[i].get().data());
    }
  }

  state.counters["n"] = n;
  state.counters["batch"] = batch_size;
  state.counters["streams"] = batch_size;
  state.counters["total_solves"] = batch_size;
}

// --------------------------------------------------------------------------
// CPU baseline: Eigen LLT on host, sequential
// --------------------------------------------------------------------------

static void BM_Batch_CPU(benchmark::State& state) {
  const Index n = state.range(0);
  const int batch_size = static_cast<int>(state.range(1));

  std::vector<Mat> As(batch_size);
  std::vector<Mat> Bs(batch_size);
  for (int i = 0; i < batch_size; ++i) {
    As[i] = make_spd(n);
    Bs[i] = Mat::Random(n, 1);
  }

  for (auto _ : state) {
    std::vector<Mat> results(batch_size);
    for (int i = 0; i < batch_size; ++i) {
      LLT<Mat> llt(As[i]);
      results[i] = llt.solve(Bs[i]);
    }
    benchmark::DoNotOptimize(results.back().data());
  }

  state.counters["n"] = n;
  state.counters["batch"] = batch_size;
  state.counters["total_solves"] = batch_size;
}

// --------------------------------------------------------------------------
// Registration
// --------------------------------------------------------------------------

// clang-format off
// Args: {matrix_size, batch_size}
// Small matrices with large batches are the interesting case for multi-stream.
BENCHMARK(BM_Batch_Sequential)->ArgsProduct({{16, 32, 64, 128, 256, 512}, {1, 4, 16, 64}})->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_Batch_Sequential_Device)->ArgsProduct({{16, 32, 64, 128, 256, 512}, {1, 4, 16, 64}})->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_Batch_MultiStream)->ArgsProduct({{16, 32, 64, 128, 256, 512}, {1, 4, 16, 64}})->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_Batch_MultiStream_AsyncDownload)->ArgsProduct({{16, 32, 64, 128, 256, 512}, {1, 4, 16, 64}})->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_Batch_CPU)->ArgsProduct({{16, 32, 64, 128, 256, 512}, {1, 4, 16, 64}})->Unit(benchmark::kMicrosecond);

// Also run larger sizes with moderate batching.
BENCHMARK(BM_Batch_MultiStream)->ArgsProduct({{512, 1024, 2048}, {1, 4, 8}})->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_Batch_MultiStream_AsyncDownload)->ArgsProduct({{512, 1024, 2048}, {1, 4, 8}})->Unit(benchmark::kMicrosecond);
// clang-format on
