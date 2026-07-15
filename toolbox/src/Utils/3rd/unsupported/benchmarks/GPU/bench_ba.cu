// Bundle Adjustment benchmark: GPU CG vs CPU CG on real BAL datasets.
//
// Tests Eigen's GPU CG pipeline (gpu::DeviceMatrix + gpu::SparseContext + DeviceScalar)
// on the normal equations (J^T*J) arising from bundle adjustment problems.
//
// Reads a BAL (Bundle Adjustment in the Large) format file, computes the
// Jacobian and residual, forms the normal equations H = J^T*J + lambda*I,
// then solves H*dx = -J^T*r with both CPU and GPU conjugate gradients.
//
// BAL format: http://grail.cs.washington.edu/projects/bal/
//
// Usage:
//   cmake --build build-bench-gpu --target bench_gpu_ba
//
//   # Download a BAL dataset (bz2-compressed):
//   wget http://grail.cs.washington.edu/projects/bal/data/ladybug/problem-49-7776-pre.txt.bz2
//   bunzip2 problem-49-7776-pre.txt.bz2
//
//   # Run on a specific problem:
//   BAL_FILE=problem-49-7776-pre.txt ./build-bench-gpu/bench_gpu_ba
//
//   # Append results to the log:
//   BAL_FILE=problem-49-7776-pre.txt ./build-bench-gpu/bench_gpu_ba \
//     --benchmark_format=console 2>&1 | tee -a benchmarks/GPU/ba_results.log
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#include <benchmark/benchmark.h>

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/GPU>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using namespace Eigen;

// ============================================================================
// BAL problem data
// ============================================================================

struct BALProblem {
  int num_cameras = 0;
  int num_points = 0;
  int num_observations = 0;

  // Observations: (camera_idx, point_idx, observed_x, observed_y).
  std::vector<int> camera_index;
  std::vector<int> point_index;
  std::vector<double> observations_x;
  std::vector<double> observations_y;

  // Camera parameters: 9 per camera (Rodrigues r[3], translation t[3], f, k1, k2).
  std::vector<double> cameras;  // [num_cameras * 9]

  // 3D points: 3 per point.
  std::vector<double> points;  // [num_points * 3]

  const double* camera(int i) const { return &cameras[i * 9]; }
  const double* point(int i) const { return &points[i * 3]; }

  bool load(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) {
      fprintf(stderr, "ERROR: Cannot open BAL file: %s\n", filename.c_str());
      return false;
    }

    in >> num_cameras >> num_points >> num_observations;
    if (!in || num_cameras <= 0 || num_points <= 0 || num_observations <= 0) {
      fprintf(stderr, "ERROR: Invalid BAL header in %s\n", filename.c_str());
      return false;
    }

    camera_index.resize(num_observations);
    point_index.resize(num_observations);
    observations_x.resize(num_observations);
    observations_y.resize(num_observations);
    for (int i = 0; i < num_observations; ++i) {
      in >> camera_index[i] >> point_index[i] >> observations_x[i] >> observations_y[i];
    }

    cameras.resize(num_cameras * 9);
    for (int i = 0; i < num_cameras * 9; ++i) {
      in >> cameras[i];
    }

    points.resize(num_points * 3);
    for (int i = 0; i < num_points * 3; ++i) {
      in >> points[i];
    }

    if (!in) {
      fprintf(stderr, "ERROR: Truncated BAL file: %s\n", filename.c_str());
      return false;
    }

    fprintf(stderr, "Loaded BAL: %d cameras, %d points, %d observations\n", num_cameras, num_points, num_observations);
    return true;
  }
};

// ============================================================================
// Camera projection model (BAL convention)
// ============================================================================

// Rodrigues rotation: rotate point X by axis-angle vector omega.
static void rodrigues_rotate(const double* omega, const double* X, double* result) {
  double theta2 = omega[0] * omega[0] + omega[1] * omega[1] + omega[2] * omega[2];
  if (theta2 > 1e-30) {
    double theta = std::sqrt(theta2);
    double costh = std::cos(theta);
    double sinth = std::sin(theta);
    double k = (1.0 - costh) / theta2;

    // Cross product omega x X.
    double wx = omega[1] * X[2] - omega[2] * X[1];
    double wy = omega[2] * X[0] - omega[0] * X[2];
    double wz = omega[0] * X[1] - omega[1] * X[0];

    // Dot product omega . X.
    double dot = omega[0] * X[0] + omega[1] * X[1] + omega[2] * X[2];

    result[0] = X[0] * costh + wx * (sinth / theta) + omega[0] * dot * k;
    result[1] = X[1] * costh + wy * (sinth / theta) + omega[1] * dot * k;
    result[2] = X[2] * costh + wz * (sinth / theta) + omega[2] * dot * k;
  } else {
    // Small angle: R ≈ I + [omega]×.
    result[0] = X[0] + omega[1] * X[2] - omega[2] * X[1];
    result[1] = X[1] + omega[2] * X[0] - omega[0] * X[2];
    result[2] = X[2] + omega[0] * X[1] - omega[1] * X[0];
  }
}

// Project a 3D point through a camera, returning the 2D residual.
// camera: [r0,r1,r2, t0,t1,t2, f, k1, k2]
// point:  [X, Y, Z]
// observed: [ox, oy]
// residual: [rx, ry] = projected - observed
static void project(const double* camera, const double* point, const double* observed, double* residual) {
  // Rotate.
  double P[3];
  rodrigues_rotate(camera, point, P);

  // Translate.
  P[0] += camera[3];
  P[1] += camera[4];
  P[2] += camera[5];

  // Normalize (BAL convention: negative z).
  double xp = -P[0] / P[2];
  double yp = -P[1] / P[2];

  // Radial distortion.
  double r2 = xp * xp + yp * yp;
  double distortion = 1.0 + camera[7] * r2 + camera[8] * r2 * r2;

  // Apply focal length.
  double predicted_x = camera[6] * distortion * xp;
  double predicted_y = camera[6] * distortion * yp;

  residual[0] = predicted_x - observed[0];
  residual[1] = predicted_y - observed[1];
}

// ============================================================================
// Jacobian computation (numerical differentiation)
// ============================================================================

// Compute the 2x9 Jacobian block w.r.t. camera params and 2x3 block w.r.t.
// point coords for a single observation, using central finite differences.
static void compute_jacobian_block(const double* camera, const double* point, const double* observed,
                                   double* J_cam,    // 2x9, row-major
                                   double* J_point)  // 2x3, row-major
{
  constexpr double eps = 1e-8;

  // Camera parameters (9).
  double cam_pert[9];
  std::copy(camera, camera + 9, cam_pert);
  for (int j = 0; j < 9; ++j) {
    double orig = cam_pert[j];
    double rp[2], rm[2];

    cam_pert[j] = orig + eps;
    project(cam_pert, point, observed, rp);
    cam_pert[j] = orig - eps;
    project(cam_pert, point, observed, rm);
    cam_pert[j] = orig;

    J_cam[0 * 9 + j] = (rp[0] - rm[0]) / (2.0 * eps);
    J_cam[1 * 9 + j] = (rp[1] - rm[1]) / (2.0 * eps);
  }

  // Point coordinates (3).
  double pt_pert[3];
  std::copy(point, point + 3, pt_pert);
  for (int j = 0; j < 3; ++j) {
    double orig = pt_pert[j];
    double rp[2], rm[2];

    pt_pert[j] = orig + eps;
    project(camera, pt_pert, observed, rp);
    pt_pert[j] = orig - eps;
    project(camera, pt_pert, observed, rm);
    pt_pert[j] = orig;

    J_point[0 * 3 + j] = (rp[0] - rm[0]) / (2.0 * eps);
    J_point[1 * 3 + j] = (rp[1] - rm[1]) / (2.0 * eps);
  }
}

// ============================================================================
// Build normal equations: H = J^T*J + lambda*I, g = -J^T*r
// ============================================================================

struct NormalEquations {
  SparseMatrix<double, ColMajor, int> H;
  VectorXd g;
  VectorXd residual;
  double residual_norm;
  int jacobian_rows;
  int jacobian_cols;
  long jacobian_nnz;
};

static NormalEquations build_normal_equations(const BALProblem& problem, double lambda = 1.0) {
  const int num_cam_params = problem.num_cameras * 9;
  const int num_pt_params = problem.num_points * 3;
  const int num_params = num_cam_params + num_pt_params;
  const int num_residuals = problem.num_observations * 2;

  fprintf(stderr, "Building Jacobian: %d x %d, %ld nonzeros\n", num_residuals, num_params,
          (long)problem.num_observations * 24);

  // Build J as a triplet list.
  using Triplet = Eigen::Triplet<double>;
  std::vector<Triplet> triplets;
  triplets.reserve(problem.num_observations * 24);  // 2 rows × 12 nonzeros = 24 entries per obs

  VectorXd residual(num_residuals);

  for (int obs = 0; obs < problem.num_observations; ++obs) {
    int ci = problem.camera_index[obs];
    int pi = problem.point_index[obs];
    double observed[2] = {problem.observations_x[obs], problem.observations_y[obs]};

    // Compute residual.
    double r[2];
    project(problem.camera(ci), problem.point(pi), observed, r);
    residual[obs * 2 + 0] = r[0];
    residual[obs * 2 + 1] = r[1];

    // Compute Jacobian blocks.
    double J_cam[18], J_pt[6];  // 2x9 and 2x3
    compute_jacobian_block(problem.camera(ci), problem.point(pi), observed, J_cam, J_pt);

    // Insert camera block: rows [2*obs, 2*obs+1], cols [9*ci, 9*ci+8].
    for (int row = 0; row < 2; ++row) {
      for (int col = 0; col < 9; ++col) {
        double val = J_cam[row * 9 + col];
        if (val != 0.0) {
          triplets.emplace_back(obs * 2 + row, ci * 9 + col, val);
        }
      }
    }

    // Insert point block: rows [2*obs, 2*obs+1], cols [num_cam_params + 3*pi, ...].
    for (int row = 0; row < 2; ++row) {
      for (int col = 0; col < 3; ++col) {
        double val = J_pt[row * 3 + col];
        if (val != 0.0) {
          triplets.emplace_back(obs * 2 + row, num_cam_params + pi * 3 + col, val);
        }
      }
    }
  }

  // Build sparse Jacobian.
  SparseMatrix<double, ColMajor, int> J(num_residuals, num_params);
  J.setFromTriplets(triplets.begin(), triplets.end());

  fprintf(stderr, "Jacobian: %dx%d, nnz=%ld\n", (int)J.rows(), (int)J.cols(), (long)J.nonZeros());

  // Form normal equations: H = J^T*J + lambda*I.
  SparseMatrix<double, ColMajor, int> H = (J.transpose() * J).pruned();

  // Add Levenberg-Marquardt damping.
  for (int i = 0; i < num_params; ++i) {
    H.coeffRef(i, i) += lambda;
  }
  H.makeCompressed();

  // Gradient: g = -J^T * r.
  VectorXd g = -(J.transpose() * residual);

  double rnorm = residual.norm();
  fprintf(stderr, "Normal equations: H is %dx%d, nnz=%ld, |r|=%.6e\n", (int)H.rows(), (int)H.cols(), (long)H.nonZeros(),
          rnorm);

  return {std::move(H), std::move(g), std::move(residual), rnorm, num_residuals, num_params, (long)J.nonZeros()};
}

// ============================================================================
// Global problem state (loaded once before benchmarks run)
// ============================================================================

static BALProblem g_problem;
static NormalEquations g_neq;
static bool g_loaded = false;

static void ensure_loaded() {
  if (g_loaded) return;

  const char* bal_file = std::getenv("BAL_FILE");
  if (!bal_file) {
    fprintf(stderr,
            "ERROR: Set BAL_FILE environment variable to a BAL problem file.\n"
            "  Download from: http://grail.cs.washington.edu/projects/bal/\n"
            "  Example:\n"
            "    wget http://grail.cs.washington.edu/projects/bal/data/ladybug/"
            "problem-49-7776-pre.txt.bz2\n"
            "    bunzip2 problem-49-7776-pre.txt.bz2\n"
            "    BAL_FILE=problem-49-7776-pre.txt ./build-bench-gpu/bench_gpu_ba\n");
    std::exit(1);
  }

  if (!g_problem.load(bal_file)) {
    std::exit(1);
  }

  g_neq = build_normal_equations(g_problem);
  g_loaded = true;
}

// ============================================================================
// CPU CG benchmark
// ============================================================================

static void BM_BA_CPU_CG(benchmark::State& state) {
  ensure_loaded();
  const auto& H = g_neq.H;
  const auto& g = g_neq.g;

  ConjugateGradient<SparseMatrix<double, ColMajor, int>, Lower | Upper> cg;
  cg.setMaxIterations(10000);
  cg.setTolerance(1e-8);
  cg.compute(H);

  int last_iters = 0;
  double last_error = 0;
  for (auto _ : state) {
    VectorXd dx = cg.solve(g);
    benchmark::DoNotOptimize(dx.data());
    last_iters = cg.iterations();
    last_error = cg.error();
  }

  state.counters["n"] = H.rows();
  state.counters["nnz"] = H.nonZeros();
  state.counters["iters"] = last_iters;
  state.counters["error"] = last_error;
  state.counters["cameras"] = g_problem.num_cameras;
  state.counters["points"] = g_problem.num_points;
  state.counters["observations"] = g_problem.num_observations;
}

// ============================================================================
// GPU CG benchmark (with Jacobi preconditioner)
// ============================================================================

static void cuda_warmup() {
  static bool done = false;
  if (!done) {
    void* p;
    cudaMalloc(&p, 1);
    cudaFree(p);
    done = true;
  }
}

static void BM_BA_GPU_CG(benchmark::State& state) {
  ensure_loaded();
  cuda_warmup();

  const auto& H = g_neq.H;
  const auto& g = g_neq.g;
  const Index n = H.rows();

  // Extract inverse diagonal (Jacobi preconditioner).
  using SpMat = SparseMatrix<double, ColMajor, int>;
  VectorXd invdiag(n);
  for (Index j = 0; j < H.outerSize(); ++j) {
    SpMat::InnerIterator it(H, j);
    while (it && it.index() != j) ++it;
    if (it && it.index() == j && it.value() != 0.0)
      invdiag(j) = 1.0 / it.value();
    else
      invdiag(j) = 1.0;
  }

  // Set up GPU context and upload data.
  gpu::Context ctx;
  gpu::Context::setThreadLocal(&ctx);
  gpu::SparseContext<double> spmv_ctx(ctx);
  auto mat = spmv_ctx.deviceView(H);
  auto d_invdiag = gpu::DeviceMatrix<double>::fromHost(invdiag, ctx.stream());
  auto d_g = gpu::DeviceMatrix<double>::fromHost(g, ctx.stream());

  int last_iters = 0;
  double last_error = 0;

  for (auto _ : state) {
    gpu::DeviceMatrix<double> d_x(n, 1);
    d_x.setZero(ctx);
    gpu::DeviceMatrix<double> residual(n, 1);
    residual.copyFrom(ctx, d_g);

    double rhsNorm2 = d_g.squaredNorm(ctx);
    double threshold = 1e-8 * 1e-8 * rhsNorm2;
    double residualNorm2 = residual.squaredNorm(ctx);

    gpu::DeviceMatrix<double> p = d_invdiag.cwiseProduct(ctx, residual);
    gpu::DeviceMatrix<double> z(n, 1), tmp(n, 1);

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

      z.cwiseProduct(ctx, d_invdiag, residual);  // in-place, no allocation
      auto absOld = std::move(absNew);
      absNew = residual.dot(ctx, z);
      auto beta = absNew / absOld;

      p *= beta;  // device-pointer scal, no host sync
      p += z;
      i++;
    }
    benchmark::DoNotOptimize(d_x.data());
    last_iters = i;
    last_error = std::sqrt(residualNorm2 / rhsNorm2);
  }

  gpu::Context::setThreadLocal(nullptr);

  state.counters["n"] = n;
  state.counters["nnz"] = H.nonZeros();
  state.counters["iters"] = last_iters;
  state.counters["error"] = last_error;
  state.counters["cameras"] = g_problem.num_cameras;
  state.counters["points"] = g_problem.num_points;
  state.counters["observations"] = g_problem.num_observations;
}

// ============================================================================
// CPU CG with Jacobi preconditioner (apples-to-apples comparison)
// ============================================================================

static void BM_BA_CPU_CG_Jacobi(benchmark::State& state) {
  ensure_loaded();
  const auto& H = g_neq.H;
  const auto& g = g_neq.g;

  // Eigen's DiagonalPreconditioner is effectively Jacobi.
  ConjugateGradient<SparseMatrix<double, ColMajor, int>, Lower | Upper> cg;
  cg.setMaxIterations(10000);
  cg.setTolerance(1e-8);
  cg.compute(H);

  int last_iters = 0;
  double last_error = 0;
  for (auto _ : state) {
    VectorXd dx = cg.solve(g);
    benchmark::DoNotOptimize(dx.data());
    last_iters = cg.iterations();
    last_error = cg.error();
  }

  state.counters["n"] = H.rows();
  state.counters["nnz"] = H.nonZeros();
  state.counters["iters"] = last_iters;
  state.counters["error"] = last_error;
}

// ============================================================================
// Register benchmarks
// ============================================================================

BENCHMARK(BM_BA_CPU_CG)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_BA_CPU_CG_Jacobi)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_BA_GPU_CG)->Unit(benchmark::kMillisecond);

// ============================================================================
// Custom main: print summary after benchmarks
// ============================================================================

int main(int argc, char** argv) {
  benchmark::Initialize(&argc, argv);

  // Print problem info before benchmarks.
  const char* bal_file = std::getenv("BAL_FILE");
  if (bal_file) {
    ensure_loaded();
    fprintf(stderr,
            "\n"
            "=== Bundle Adjustment GPU CG Benchmark ===\n"
            "BAL file:      %s\n"
            "Cameras:       %d\n"
            "Points:        %d\n"
            "Observations:  %d\n"
            "J size:        %d x %d, nnz=%ld\n"
            "H size:        %d x %d, nnz=%ld\n"
            "|residual|:    %.6e\n"
            "==========================================\n\n",
            bal_file, g_problem.num_cameras, g_problem.num_points, g_problem.num_observations, g_neq.jacobian_rows,
            g_neq.jacobian_cols, g_neq.jacobian_nnz, (int)g_neq.H.rows(), (int)g_neq.H.cols(), (long)g_neq.H.nonZeros(),
            g_neq.residual_norm);
  }

  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return 0;
}
