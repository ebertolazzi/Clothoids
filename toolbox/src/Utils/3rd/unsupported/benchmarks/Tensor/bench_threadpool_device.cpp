// Microbenchmarks for ThreadPoolDevice scheduling overhead.
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#define EIGEN_USE_THREADS

#include <atomic>
#include <map>
#include <memory>
#include <mutex>

#include <benchmark/benchmark.h>
#include <unsupported/Eigen/CXX11/ThreadPool>
#include <unsupported/Eigen/CXX11/Tensor>

using Eigen::Index;
using Eigen::TensorOpCost;
using Eigen::ThreadPool;
using Eigen::ThreadPoolDevice;

// Returns a per-thread-count singleton device, lazily created on first use.
// Pools are reused across benchmark iterations so thread-creation cost does
// not enter the measurement, and varying thread counts across benchmarks is
// supported without later calls silently colliding with the first one.
static ThreadPoolDevice& device(int threads) {
  static std::mutex mu;
  static std::map<int, std::unique_ptr<ThreadPool>> pools;
  static std::map<int, std::unique_ptr<ThreadPoolDevice>> devs;
  std::lock_guard<std::mutex> lock(mu);
  auto it = devs.find(threads);
  if (it != devs.end()) return *it->second;
  auto pool = std::unique_ptr<ThreadPool>(new ThreadPool(threads));
  auto* dev = new ThreadPoolDevice(pool.get(), threads);
  pools.emplace(threads, std::move(pool));
  devs.emplace(threads, std::unique_ptr<ThreadPoolDevice>(dev));
  return *dev;
}

// ---- enqueue overhead: empty lambda, wait via shared counter -------------
static void BM_EnqueueEmpty(benchmark::State& state) {
  const int threads = state.range(0);
  auto& dev = device(threads);
  const int batch = 1024;
  for (auto _ : state) {
    std::atomic<int> remaining{batch};
    for (int i = 0; i < batch; ++i) {
      dev.enqueue([&remaining] { remaining.fetch_sub(1, std::memory_order_release); });
    }
    while (remaining.load(std::memory_order_acquire) != 0) {
    }
  }
  state.SetItemsProcessed(state.iterations() * batch);
  state.counters["threads"] = threads;
}
BENCHMARK(BM_EnqueueEmpty)->Arg(8)->UseRealTime();

// ---- parallelFor with tiny per-element cost: dispatch overhead dominates --
static void BM_ParallelForTiny(benchmark::State& state) {
  const Index n = state.range(0);
  const int threads = state.range(1);
  auto& dev = device(threads);
  // Cost just large enough that the cost model decides to parallelize.
  const TensorOpCost cost(/*bytes_loaded=*/4, /*bytes_stored=*/4, /*compute_cycles=*/4);
  std::atomic<Index> sink{0};
  for (auto _ : state) {
    dev.parallelFor(n, cost,
                    [&sink](Index first, Index last) { sink.fetch_add(last - first, std::memory_order_relaxed); });
  }
  benchmark::DoNotOptimize(sink.load());
  state.SetItemsProcessed(state.iterations() * n);
  state.counters["threads"] = threads;
}
BENCHMARK(BM_ParallelForTiny)->Args({256, 8})->Args({4096, 8})->Args({65536, 8})->Args({1048576, 8})->UseRealTime();

// ---- parallelFor with realistic per-element compute work ------------------
static void BM_ParallelForCompute(benchmark::State& state) {
  const Index n = state.range(0);
  const int threads = state.range(1);
  auto& dev = device(threads);
  const TensorOpCost cost(/*bytes_loaded=*/8, /*bytes_stored=*/8, /*compute_cycles=*/32);
  std::vector<double> data(n, 1.0);
  for (auto _ : state) {
    dev.parallelFor(n, cost, [&data](Index first, Index last) {
      for (Index i = first; i < last; ++i) {
        data[i] = data[i] * 1.0000001 + 1e-12;
      }
    });
    benchmark::DoNotOptimize(data.data());
  }
  state.SetItemsProcessed(state.iterations() * n);
  state.counters["threads"] = threads;
}
BENCHMARK(BM_ParallelForCompute)->Args({4096, 8})->Args({65536, 8})->Args({1048576, 8})->UseRealTime();

// ---- parallelForAsync: measures async dispatch + heap context overhead ----
static void BM_ParallelForAsync(benchmark::State& state) {
  const Index n = state.range(0);
  const int threads = state.range(1);
  auto& dev = device(threads);
  const TensorOpCost cost(8, 8, 4);
  for (auto _ : state) {
    std::atomic<int> done_flag{0};
    dev.parallelForAsync(
        n, cost, [](Index, Index) {}, [&done_flag] { done_flag.store(1, std::memory_order_release); });
    while (done_flag.load(std::memory_order_acquire) == 0) {
    }
  }
  state.SetItemsProcessed(state.iterations() * n);
  state.counters["threads"] = threads;
}
BENCHMARK(BM_ParallelForAsync)->Args({4096, 8})->Args({65536, 8})->Args({1048576, 8})->UseRealTime();

// ---- ThreadPoolDevice::memcpy at various sizes ----------------------------
static void BM_DeviceMemcpy(benchmark::State& state) {
  const size_t bytes = state.range(0);
  const int threads = state.range(1);
  auto& dev = device(threads);
  std::vector<char> src(bytes), dst(bytes);
  for (auto _ : state) {
    dev.memcpy(dst.data(), src.data(), bytes);
    benchmark::DoNotOptimize(dst.data());
  }
  state.SetBytesProcessed(state.iterations() * bytes);
  state.counters["threads"] = threads;
}
BENCHMARK(BM_DeviceMemcpy)
    ->Args({16 << 10, 8})
    ->Args({256 << 10, 8})
    ->Args({4 << 20, 8})
    ->Args({64 << 20, 8})
    ->UseRealTime();

BENCHMARK_MAIN();
