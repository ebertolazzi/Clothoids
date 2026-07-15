// GPU chaining benchmarks: measure async pipeline efficiency.
//
// Compares:
//   1. Host round-trip per solve (baseline)
//   2. DeviceMatrix chaining (no host round-trip between solves)
//   3. Varying chain lengths (1, 2, 4, 8 consecutive solves)
//
// For Nsight Systems: look for gaps between kernel launches in the timeline.
// Host round-trip creates visible idle gaps; chaining should show back-to-back kernels.
//
//   nsys profile --trace=cuda,nvtx ./bench_chaining
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#include <benchmark/benchmark.h>

#include <Eigen/Cholesky>
#include <unsupported/Eigen/GPU>

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
// Baseline: host round-trip between every solve
// --------------------------------------------------------------------------

static void BM_Chain_HostRoundtrip(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  const int chain_len = static_cast<int>(state.range(1));

  Mat A = make_spd(n);
  Mat B = Mat::Random(n, 1);
  gpu::LLT<Scalar> llt(A);

  for (auto _ : state) {
    Mat X = B;
    for (int i = 0; i < chain_len; ++i) {
      X = llt.solve(X);  // host → device → host each time
    }
    benchmark::DoNotOptimize(X.data());
  }

  state.counters["n"] = n;
  state.counters["chain"] = chain_len;
  state.counters["solves/iter"] = chain_len;
}

// --------------------------------------------------------------------------
// DeviceMatrix chaining: no host round-trip between solves
// --------------------------------------------------------------------------

static void BM_Chain_Device(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  const int chain_len = static_cast<int>(state.range(1));

  Mat A = make_spd(n);
  Mat B = Mat::Random(n, 1);
  gpu::LLT<Scalar> llt(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  for (auto _ : state) {
    gpu::DeviceMatrix<Scalar> d_X = llt.solve(d_B);
    for (int i = 1; i < chain_len; ++i) {
      d_X = llt.solve(d_X);  // device → device, fully async
    }
    Mat X = d_X.toHost();  // single sync at end
    benchmark::DoNotOptimize(X.data());
  }

  state.counters["n"] = n;
  state.counters["chain"] = chain_len;
  state.counters["solves/iter"] = chain_len;
}

// --------------------------------------------------------------------------
// DeviceMatrix chaining with async download (overlap D2H with next iteration)
//
// Double-buffered: each loop body issues iteration N+1's chain *before*
// draining iteration N's D2H, so the host wait overlaps with on-device
// compute instead of stalling the pipeline.
// --------------------------------------------------------------------------

static void BM_Chain_DeviceAsync(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  const int chain_len = static_cast<int>(state.range(1));

  Mat A = make_spd(n);
  Mat B = Mat::Random(n, 1);
  gpu::LLT<Scalar> llt(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  auto run_chain = [&]() {
    gpu::DeviceMatrix<Scalar> d_X = llt.solve(d_B);
    for (int i = 1; i < chain_len; ++i) {
      d_X = llt.solve(d_X);
    }
    return d_X.toHostAsync();
  };

  // Prime the pipeline so each timed iteration overlaps a fresh chain
  // with the previous iteration's D2H.
  auto prev = run_chain();

  for (auto _ : state) {
    auto next = run_chain();  // kick off N+1 while D2H of N is in flight
    Mat X = prev.get();       // drain N (overlaps with `next` compute)
    benchmark::DoNotOptimize(X.data());
    prev = std::move(next);
  }

  // Drain the trailing transfer outside the timed region.
  prev.get();

  state.counters["n"] = n;
  state.counters["chain"] = chain_len;
  state.counters["solves/iter"] = chain_len;
}

// --------------------------------------------------------------------------
// Pure GPU chain (no download — measures kernel-only throughput)
// --------------------------------------------------------------------------

static void BM_Chain_DeviceNoDownload(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  const int chain_len = static_cast<int>(state.range(1));

  Mat A = make_spd(n);
  Mat B = Mat::Random(n, 1);
  gpu::LLT<Scalar> llt(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  for (auto _ : state) {
    gpu::DeviceMatrix<Scalar> d_X = llt.solve(d_B);
    for (int i = 1; i < chain_len; ++i) {
      d_X = llt.solve(d_X);
    }
    EIGEN_CUDA_RUNTIME_CHECK(cudaStreamSynchronize(llt.stream()));
    benchmark::DoNotOptimize(d_X.data());
  }

  state.counters["n"] = n;
  state.counters["chain"] = chain_len;
  state.counters["solves/iter"] = chain_len;
}

// --------------------------------------------------------------------------
// Compute + solve chain (full pipeline: factorize, then chain solves)
// --------------------------------------------------------------------------

static void BM_FullPipeline_Host(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  const int chain_len = static_cast<int>(state.range(1));

  Mat A = make_spd(n);
  Mat B = Mat::Random(n, 1);

  for (auto _ : state) {
    gpu::LLT<Scalar> llt(A);
    Mat X = B;
    for (int i = 0; i < chain_len; ++i) {
      X = llt.solve(X);
    }
    benchmark::DoNotOptimize(X.data());
  }

  state.counters["n"] = n;
  state.counters["chain"] = chain_len;
}

static void BM_FullPipeline_Device(benchmark::State& state) {
  cuda_warmup();
  const Index n = state.range(0);
  const int chain_len = static_cast<int>(state.range(1));

  Mat A = make_spd(n);
  Mat B = Mat::Random(n, 1);

  for (auto _ : state) {
    auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
    auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);
    gpu::LLT<Scalar> llt;
    llt.compute(d_A);
    gpu::DeviceMatrix<Scalar> d_X = llt.solve(d_B);
    for (int i = 1; i < chain_len; ++i) {
      d_X = llt.solve(d_X);
    }
    Mat X = d_X.toHost();
    benchmark::DoNotOptimize(X.data());
  }

  state.counters["n"] = n;
  state.counters["chain"] = chain_len;
}

// --------------------------------------------------------------------------
// Registration
// --------------------------------------------------------------------------

// clang-format off
// Args: {matrix_size, chain_length}
BENCHMARK(BM_Chain_HostRoundtrip)->ArgsProduct({{64, 256, 1024, 4096}, {1, 2, 4, 8}})->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_Chain_Device)->ArgsProduct({{64, 256, 1024, 4096}, {1, 2, 4, 8}})->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_Chain_DeviceAsync)->ArgsProduct({{64, 256, 1024, 4096}, {1, 2, 4, 8}})->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_Chain_DeviceNoDownload)->ArgsProduct({{64, 256, 1024, 4096}, {1, 2, 4, 8}})->Unit(benchmark::kMicrosecond);

BENCHMARK(BM_FullPipeline_Host)->ArgsProduct({{256, 1024, 4096}, {1, 4}})->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_FullPipeline_Device)->ArgsProduct({{256, 1024, 4096}, {1, 4}})->Unit(benchmark::kMicrosecond);
// clang-format on
