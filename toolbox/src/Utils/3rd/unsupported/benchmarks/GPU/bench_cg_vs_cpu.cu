// Benchmark: GPU CG vs CPU CG on realistic sparse systems.
//
// Tests 2D Laplacian (5-point stencil) and 3D Laplacian (7-point stencil)
// in both float and double precision.
//
// Usage:
//   cmake --build build-bench-gpu --target bench_gpu_cg_vs_cpu
//   ./build-bench-gpu/bench_gpu_cg_vs_cpu
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#include <benchmark/benchmark.h>

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/GPU>

using namespace Eigen;

// ---- Sparse matrix generators -----------------------------------------------

template <typename Scalar>
SparseMatrix<Scalar, ColMajor, int> make_laplacian_2d(int grid_n) {
  using SpMat = SparseMatrix<Scalar, ColMajor, int>;
  const int n = grid_n * grid_n;
  SpMat A(n, n);
  A.reserve(VectorXi::Constant(n, 5));
  for (int i = 0; i < grid_n; ++i) {
    for (int j = 0; j < grid_n; ++j) {
      int idx = i * grid_n + j;
      A.insert(idx, idx) = Scalar(4);
      if (i > 0) A.insert(idx, idx - grid_n) = Scalar(-1);
      if (i < grid_n - 1) A.insert(idx, idx + grid_n) = Scalar(-1);
      if (j > 0) A.insert(idx, idx - 1) = Scalar(-1);
      if (j < grid_n - 1) A.insert(idx, idx + 1) = Scalar(-1);
    }
  }
  A.makeCompressed();
  return A;
}

template <typename Scalar>
SparseMatrix<Scalar, ColMajor, int> make_laplacian_3d(int grid_n) {
  using SpMat = SparseMatrix<Scalar, ColMajor, int>;
  const int n = grid_n * grid_n * grid_n;
  const int n2 = grid_n * grid_n;
  SpMat A(n, n);
  A.reserve(VectorXi::Constant(n, 7));
  for (int i = 0; i < grid_n; ++i) {
    for (int j = 0; j < grid_n; ++j) {
      for (int k = 0; k < grid_n; ++k) {
        int idx = i * n2 + j * grid_n + k;
        A.insert(idx, idx) = Scalar(6);
        if (i > 0) A.insert(idx, idx - n2) = Scalar(-1);
        if (i < grid_n - 1) A.insert(idx, idx + n2) = Scalar(-1);
        if (j > 0) A.insert(idx, idx - grid_n) = Scalar(-1);
        if (j < grid_n - 1) A.insert(idx, idx + grid_n) = Scalar(-1);
        if (k > 0) A.insert(idx, idx - 1) = Scalar(-1);
        if (k < grid_n - 1) A.insert(idx, idx + 1) = Scalar(-1);
      }
    }
  }
  A.makeCompressed();
  return A;
}

static void cuda_warmup() {
  static bool done = false;
  if (!done) {
    void* p;
    cudaMalloc(&p, 1);
    cudaFree(p);
    done = true;
  }
}

// ---- CPU CG -----------------------------------------------------------------

template <typename Scalar, typename MatGen>
void run_cpu_cg(benchmark::State& state, MatGen make_matrix) {
  using SpMat = SparseMatrix<Scalar, ColMajor, int>;
  using Vec = Matrix<Scalar, Dynamic, 1>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  const int grid_n = state.range(0);
  SpMat A = make_matrix(grid_n);
  Vec b = Vec::Random(A.rows());

  ConjugateGradient<SpMat, Lower | Upper> cg;
  cg.setMaxIterations(10000);
  cg.setTolerance(RealScalar(1e-8));
  cg.compute(A);

  int last_iters = 0;
  for (auto _ : state) {
    Vec x = cg.solve(b);
    benchmark::DoNotOptimize(x.data());
    last_iters = cg.iterations();
  }
  state.counters["n"] = A.rows();
  state.counters["nnz"] = A.nonZeros();
  state.counters["iters"] = last_iters;
  state.counters["error"] = cg.error();
}

// ---- GPU CG -----------------------------------------------------------------

template <typename Scalar, typename MatGen>
void run_gpu_cg(benchmark::State& state, MatGen make_matrix) {
  using SpMat = SparseMatrix<Scalar, ColMajor, int>;
  using Vec = Matrix<Scalar, Dynamic, 1>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  cuda_warmup();
  const int grid_n = state.range(0);
  SpMat A = make_matrix(grid_n);
  const Index n = A.rows();
  Vec b = Vec::Random(n);

  // Extract inverse diagonal.
  Vec invdiag(n);
  for (Index j = 0; j < A.outerSize(); ++j) {
    typename SpMat::InnerIterator it(A, j);
    while (it && it.index() != j) ++it;
    if (it && it.index() == j && it.value() != Scalar(0))
      invdiag(j) = Scalar(1) / it.value();
    else
      invdiag(j) = Scalar(1);
  }

  gpu::Context ctx;
  gpu::Context::setThreadLocal(&ctx);
  gpu::SparseContext<Scalar> spmv_ctx(ctx);
  auto mat = spmv_ctx.deviceView(A);
  auto d_invdiag = gpu::DeviceMatrix<Scalar>::fromHost(invdiag, ctx.stream());
  auto d_b = gpu::DeviceMatrix<Scalar>::fromHost(b, ctx.stream());

  int last_iters = 0;
  RealScalar last_error = 0;

  for (auto _ : state) {
    gpu::DeviceMatrix<Scalar> d_x(n, 1);
    d_x.setZero(ctx);
    gpu::DeviceMatrix<Scalar> residual(n, 1);
    residual.copyFrom(ctx, d_b);

    RealScalar rhsNorm2 = d_b.squaredNorm(ctx);
    RealScalar tol = RealScalar(1e-8);
    RealScalar threshold = tol * tol * rhsNorm2;
    RealScalar residualNorm2 = residual.squaredNorm(ctx);

    gpu::DeviceMatrix<Scalar> p = d_invdiag.cwiseProduct(ctx, residual);
    gpu::DeviceMatrix<Scalar> z(n, 1), tmp(n, 1);

    auto absNew = residual.dot(ctx, p);
    Index i = 0;
    Index maxIters = 10000;
    while (i < maxIters) {
      tmp.noalias() = mat * p;
      auto alpha = absNew / p.dot(ctx, tmp);
      d_x += alpha * p;
      residual -= alpha * tmp;

      residualNorm2 = residual.squaredNorm(ctx);
      if (residualNorm2 < threshold) break;

      z.cwiseProduct(ctx, d_invdiag, residual);
      auto absOld = std::move(absNew);
      absNew = residual.dot(ctx, z);
      auto beta = absNew / absOld;

      p *= beta;
      p += z;
      i++;
    }
    benchmark::DoNotOptimize(d_x.data());
    last_iters = i;
    last_error = numext::sqrt(residualNorm2 / rhsNorm2);
  }

  gpu::Context::setThreadLocal(nullptr);
  state.counters["n"] = n;
  state.counters["nnz"] = A.nonZeros();
  state.counters["iters"] = last_iters;
  state.counters["error"] = last_error;
}

// ---- 2D Laplacian, double ---------------------------------------------------

static void BM_CG_CPU_2D_double(benchmark::State& state) { run_cpu_cg<double>(state, make_laplacian_2d<double>); }
static void BM_CG_GPU_2D_double(benchmark::State& state) { run_gpu_cg<double>(state, make_laplacian_2d<double>); }

BENCHMARK(BM_CG_CPU_2D_double)->ArgsProduct({{32, 64, 128, 256, 512}});
BENCHMARK(BM_CG_GPU_2D_double)->ArgsProduct({{32, 64, 128, 256, 512}});

// ---- 2D Laplacian, float ----------------------------------------------------

static void BM_CG_CPU_2D_float(benchmark::State& state) { run_cpu_cg<float>(state, make_laplacian_2d<float>); }
static void BM_CG_GPU_2D_float(benchmark::State& state) { run_gpu_cg<float>(state, make_laplacian_2d<float>); }

BENCHMARK(BM_CG_CPU_2D_float)->ArgsProduct({{32, 64, 128, 256, 512}});
BENCHMARK(BM_CG_GPU_2D_float)->ArgsProduct({{32, 64, 128, 256, 512}});

// ---- 3D Laplacian, double ---------------------------------------------------

static void BM_CG_CPU_3D_double(benchmark::State& state) { run_cpu_cg<double>(state, make_laplacian_3d<double>); }
static void BM_CG_GPU_3D_double(benchmark::State& state) { run_gpu_cg<double>(state, make_laplacian_3d<double>); }

BENCHMARK(BM_CG_CPU_3D_double)->ArgsProduct({{16, 32, 48, 64}});
BENCHMARK(BM_CG_GPU_3D_double)->ArgsProduct({{16, 32, 48, 64}});

// ---- 3D Laplacian, float ----------------------------------------------------

static void BM_CG_CPU_3D_float(benchmark::State& state) { run_cpu_cg<float>(state, make_laplacian_3d<float>); }
static void BM_CG_GPU_3D_float(benchmark::State& state) { run_gpu_cg<float>(state, make_laplacian_3d<float>); }

BENCHMARK(BM_CG_CPU_3D_float)->ArgsProduct({{16, 32, 48, 64}});
BENCHMARK(BM_CG_GPU_3D_float)->ArgsProduct({{16, 32, 48, 64}});
