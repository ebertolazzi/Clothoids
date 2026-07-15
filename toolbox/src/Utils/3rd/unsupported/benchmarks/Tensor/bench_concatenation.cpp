// Benchmarks for Eigen Tensor concatenation: standalone assigns and
// compound expressions that read from a concatenation. Tests both the
// outer-axis path (each operand slab is a single contiguous run, so
// TensorBlockIO collapses to one memcpy per slab) and the inner-axis
// path (the slab walk per row is short).
//
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#define EIGEN_USE_THREADS

#include <benchmark/benchmark.h>
#include <unsupported/Eigen/Tensor>
#include <unsupported/Eigen/ThreadPool>

using namespace Eigen;

typedef float Scalar;

// --- Standalone concat along axis 0 (innermost for ColMajor) ---
static void BM_Concat2D_Axis0(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);

  Tensor<Scalar, 2> A(M, N);
  Tensor<Scalar, 2> B(M, N);
  Tensor<Scalar, 2> dst(2 * M, N);
  A.setRandom();
  B.setRandom();

  for (auto _ : state) {
    dst = A.concatenate(B, 0);
    benchmark::DoNotOptimize(dst.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * 2 * M * N * sizeof(Scalar) * 2);
}

// --- Standalone concat along axis 1 (outermost for ColMajor 2D) ---
static void BM_Concat2D_Axis1(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);

  Tensor<Scalar, 2> A(M, N);
  Tensor<Scalar, 2> B(M, N);
  Tensor<Scalar, 2> dst(M, 2 * N);
  A.setRandom();
  B.setRandom();

  for (auto _ : state) {
    dst = A.concatenate(B, 1);
    benchmark::DoNotOptimize(dst.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * 2 * N * sizeof(Scalar) * 2);
}

// --- Standalone concat 3D along axis 2 (outermost for ColMajor) ---
static void BM_Concat3D_OuterAxis(benchmark::State& state) {
  const int D0 = state.range(0);
  const int D1 = state.range(1);
  const int D2 = state.range(2);

  Tensor<Scalar, 3> A(D0, D1, D2);
  Tensor<Scalar, 3> B(D0, D1, D2);
  Tensor<Scalar, 3> dst(D0, D1, 2 * D2);
  A.setRandom();
  B.setRandom();

  for (auto _ : state) {
    dst = A.concatenate(B, 2);
    benchmark::DoNotOptimize(dst.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * D0 * D1 * 2 * D2 * sizeof(Scalar) * 2);
}

// --- Compound: (concat(A, B) + C) — concat feeds cwise binary ---
static void BM_Concat2D_Plus(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);

  Tensor<Scalar, 2> A(M, N);
  Tensor<Scalar, 2> B(M, N);
  Tensor<Scalar, 2> C(2 * M, N);
  Tensor<Scalar, 2> dst(2 * M, N);
  A.setRandom();
  B.setRandom();
  C.setRandom();

  for (auto _ : state) {
    dst = A.concatenate(B, 0) + C;
    benchmark::DoNotOptimize(dst.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * 2 * M * N * sizeof(Scalar) * 3);
}

// --- Compound: ((concat(A, B) + C) * D).sqrt() — concat in longer chain ---
static void BM_Concat2D_Chain(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);

  Tensor<Scalar, 2> A(M, N);
  Tensor<Scalar, 2> B(M, N);
  Tensor<Scalar, 2> C(2 * M, N);
  Tensor<Scalar, 2> D(2 * M, N);
  Tensor<Scalar, 2> dst(2 * M, N);
  A.setRandom();
  B.setRandom();
  C.setRandom();
  C = C.abs();
  D.setRandom();
  D = D.abs();

  for (auto _ : state) {
    dst = ((A.concatenate(B, 0) + C) * D).sqrt();
    benchmark::DoNotOptimize(dst.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * 2 * M * N * sizeof(Scalar) * 4);
}

// --- Compound: (concat(A, B) along axis 1 + C) — outer-axis concat feeds a
// cwise binary. Each operand slab is a contiguous run, so block evaluation
// returns it as a zero-copy view and the chain stays fully streaming
// (no scratch round-trip), unlike the inner-axis BM_Concat2D_Plus. ---
static void BM_Concat2D_Axis1_Plus(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);

  Tensor<Scalar, 2> A(M, N);
  Tensor<Scalar, 2> B(M, N);
  Tensor<Scalar, 2> C(M, 2 * N);
  Tensor<Scalar, 2> dst(M, 2 * N);
  A.setRandom();
  B.setRandom();
  C.setRandom();

  for (auto _ : state) {
    dst = A.concatenate(B, 1) + C;
    benchmark::DoNotOptimize(dst.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * 2 * N * sizeof(Scalar) * 3);
}

// --- Compound: ((concat(A, B) along axis 1 + C) * D).sqrt() — outer-axis
// analogue of BM_Concat2D_Chain. The zero-copy view path keeps this case off
// the materialize-then-reread trade-off that BM_Concat2D_Chain pays. ---
static void BM_Concat2D_Axis1_Chain(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);

  Tensor<Scalar, 2> A(M, N);
  Tensor<Scalar, 2> B(M, N);
  Tensor<Scalar, 2> C(M, 2 * N);
  Tensor<Scalar, 2> D(M, 2 * N);
  Tensor<Scalar, 2> dst(M, 2 * N);
  A.setRandom();
  B.setRandom();
  C.setRandom();
  C = C.abs();
  D.setRandom();
  D = D.abs();

  for (auto _ : state) {
    dst = ((A.concatenate(B, 1) + C) * D).sqrt();
    benchmark::DoNotOptimize(dst.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * 2 * N * sizeof(Scalar) * 4);
}

// --- Compound: concat of two raw tensors with an outer reduction-ish op ---
// .square() forces another cwise pass after the concat.
static void BM_Concat3D_OuterAxis_Square(benchmark::State& state) {
  const int D0 = state.range(0);
  const int D1 = state.range(1);
  const int D2 = state.range(2);

  Tensor<Scalar, 3> A(D0, D1, D2);
  Tensor<Scalar, 3> B(D0, D1, D2);
  Tensor<Scalar, 3> dst(D0, D1, 2 * D2);
  A.setRandom();
  B.setRandom();

  for (auto _ : state) {
    dst = A.concatenate(B, 2).square();
    benchmark::DoNotOptimize(dst.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * D0 * D1 * 2 * D2 * sizeof(Scalar) * 2);
}

// --- Threaded standalone concat ---
static void BM_Concat2D_Axis1_ThreadPool(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);
  const int threads = state.range(2);

  Tensor<Scalar, 2> A(M, N);
  Tensor<Scalar, 2> B(M, N);
  Tensor<Scalar, 2> dst(M, 2 * N);
  A.setRandom();
  B.setRandom();

  ThreadPool tp(threads);
  ThreadPoolDevice dev(&tp, threads);

  for (auto _ : state) {
    dst.device(dev) = A.concatenate(B, 1);
    benchmark::DoNotOptimize(dst.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * 2 * N * sizeof(Scalar) * 2);
  state.counters["threads"] = threads;
}

// --- Threaded compound: (concat(A, B) + C) * D ---
static void BM_Concat2D_Chain_ThreadPool(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);
  const int threads = state.range(2);

  Tensor<Scalar, 2> A(M, N);
  Tensor<Scalar, 2> B(M, N);
  Tensor<Scalar, 2> C(2 * M, N);
  Tensor<Scalar, 2> D(2 * M, N);
  Tensor<Scalar, 2> dst(2 * M, N);
  A.setRandom();
  B.setRandom();
  C.setRandom();
  D.setRandom();

  ThreadPool tp(threads);
  ThreadPoolDevice dev(&tp, threads);

  for (auto _ : state) {
    dst.device(dev) = (A.concatenate(B, 0) + C) * D;
    benchmark::DoNotOptimize(dst.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * 2 * M * N * sizeof(Scalar) * 4);
  state.counters["threads"] = threads;
}

// --- Threaded outer-axis compound: (concat(A, B) along axis 1 + C) * D ---
static void BM_Concat2D_Axis1_Chain_ThreadPool(benchmark::State& state) {
  const int M = state.range(0);
  const int N = state.range(1);
  const int threads = state.range(2);

  Tensor<Scalar, 2> A(M, N);
  Tensor<Scalar, 2> B(M, N);
  Tensor<Scalar, 2> C(M, 2 * N);
  Tensor<Scalar, 2> D(M, 2 * N);
  Tensor<Scalar, 2> dst(M, 2 * N);
  A.setRandom();
  B.setRandom();
  C.setRandom();
  D.setRandom();

  ThreadPool tp(threads);
  ThreadPoolDevice dev(&tp, threads);

  for (auto _ : state) {
    dst.device(dev) = (A.concatenate(B, 1) + C) * D;
    benchmark::DoNotOptimize(dst.data());
    benchmark::ClobberMemory();
  }
  state.SetBytesProcessed(state.iterations() * M * 2 * N * sizeof(Scalar) * 4);
  state.counters["threads"] = threads;
}

// clang-format off
#define CONCAT2D_SIZES \
  ->Args({256, 256})->Args({1024, 1024})->Args({4096, 4096})
#define CONCAT3D_SIZES \
  ->Args({64, 64, 64})->Args({128, 128, 128})->Args({256, 256, 64})
#define CONCAT2D_THREAD_SIZES \
  ->Args({1024, 1024, 1})->Args({1024, 1024, 4})->Args({1024, 1024, 8}) \
  ->Args({4096, 4096, 1})->Args({4096, 4096, 4})->Args({4096, 4096, 8})
// clang-format on

BENCHMARK(BM_Concat2D_Axis0) CONCAT2D_SIZES->UseRealTime();
BENCHMARK(BM_Concat2D_Axis1) CONCAT2D_SIZES->UseRealTime();
BENCHMARK(BM_Concat3D_OuterAxis) CONCAT3D_SIZES->UseRealTime();
BENCHMARK(BM_Concat2D_Plus) CONCAT2D_SIZES->UseRealTime();
BENCHMARK(BM_Concat2D_Chain) CONCAT2D_SIZES->UseRealTime();
BENCHMARK(BM_Concat2D_Axis1_Plus) CONCAT2D_SIZES->UseRealTime();
BENCHMARK(BM_Concat2D_Axis1_Chain) CONCAT2D_SIZES->UseRealTime();
BENCHMARK(BM_Concat3D_OuterAxis_Square) CONCAT3D_SIZES->UseRealTime();
BENCHMARK(BM_Concat2D_Axis1_ThreadPool) CONCAT2D_THREAD_SIZES->UseRealTime();
BENCHMARK(BM_Concat2D_Chain_ThreadPool) CONCAT2D_THREAD_SIZES->UseRealTime();
BENCHMARK(BM_Concat2D_Axis1_Chain_ThreadPool) CONCAT2D_THREAD_SIZES->UseRealTime();

BENCHMARK_MAIN();
