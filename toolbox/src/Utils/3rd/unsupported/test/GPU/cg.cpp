// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Rasmus Munk Larsen <rmlarsen@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

// End-to-end test: CG algorithm running on GPU via gpu::DeviceMatrix.
//
// Uses DeviceSparseView for SpMV, gpu::DeviceMatrix for vectors, DeviceScalar
// for deferred reductions. Verifies correctness against CPU ConjugateGradient.

#define EIGEN_USE_GPU
#include "main.h"
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/GPU>
#include "gpu_test_helpers.h"

using namespace Eigen;

// ---- Helper: build a sparse SPD matrix --------------------------------------

template <typename Scalar>
SparseMatrix<Scalar, ColMajor, int> make_spd(Index n, double density = 0.1) {
  using SpMat = SparseMatrix<Scalar, ColMajor, int>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  SpMat R(n, n);
  R.reserve(VectorXi::Constant(n, static_cast<int>(n * density) + 1));
  for (Index j = 0; j < n; ++j) {
    for (Index i = 0; i < n; ++i) {
      if (i == j || (std::rand() / double(RAND_MAX)) < density) {
        R.insert(i, j) = Scalar(std::rand() / double(RAND_MAX) - 0.5);
      }
    }
  }
  R.makeCompressed();
  SpMat A = R.adjoint() * R;
  for (Index i = 0; i < n; ++i) A.coeffRef(i, i) += Scalar(RealScalar(n));
  A.makeCompressed();
  return A;
}

// ---- GPU CG without preconditioner ------------------------------------------

template <typename Scalar>
void test_gpu_cg(Index n) {
  using SpMat = SparseMatrix<Scalar, ColMajor, int>;
  using Vec = Matrix<Scalar, Dynamic, 1>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  SpMat A = make_spd<Scalar>(n);
  Vec b = Vec::Random(n);

  // CPU reference (identity preconditioner to match GPU).
  ConjugateGradient<SpMat, Lower | Upper, IdentityPreconditioner> cpu_cg;
  cpu_cg.setMaxIterations(1000);
  cpu_cg.setTolerance(RealScalar(1e-8));
  cpu_cg.compute(A);
  Vec x_cpu = cpu_cg.solve(b);
  VERIFY_IS_EQUAL(cpu_cg.info(), Success);

  // GPU CG: mirrors Eigen's conjugate_gradient() using gpu::DeviceMatrix ops.
  gpu::Context ctx;
  gpu::Context::setThreadLocal(&ctx);
  gpu::SparseContext<Scalar> spmv_ctx(ctx);
  auto mat = spmv_ctx.deviceView(A);

  auto d_b = gpu::DeviceMatrix<Scalar>::fromHost(b, ctx.stream());
  gpu::DeviceMatrix<Scalar> d_x(n, 1);
  d_x.setZero(ctx);

  // r = b (since x=0)
  gpu::DeviceMatrix<Scalar> residual(n, 1);
  residual.copyFrom(ctx, d_b);

  RealScalar rhsNorm2 = d_b.squaredNorm(ctx);
  RealScalar tol = RealScalar(1e-8);
  RealScalar threshold = tol * tol * rhsNorm2;
  RealScalar residualNorm2 = residual.squaredNorm(ctx);

  // p = r (no preconditioner)
  gpu::DeviceMatrix<Scalar> p(n, 1);
  p.copyFrom(ctx, residual);
  gpu::DeviceMatrix<Scalar> z(n, 1), tmp(n, 1);

  auto absNew = residual.dot(ctx, p);
  Index maxIters = 1000;
  Index i = 0;
  while (i < maxIters) {
    tmp.noalias() = mat * p;

    auto alpha = absNew / p.dot(ctx, tmp);
    d_x += alpha * p;
    residual -= alpha * tmp;

    residualNorm2 = residual.squaredNorm(ctx);
    if (residualNorm2 < threshold) break;

    // z = r (no preconditioner)
    z.copyFrom(ctx, residual);

    auto absOld = std::move(absNew);
    absNew = residual.dot(ctx, z);
    auto beta = absNew / absOld;

    p *= beta;
    p += z;
    i++;
  }

  gpu::Context::setThreadLocal(nullptr);

  Vec x_gpu = d_x.toHost(ctx.stream());

  // Verify residual.
  Vec r = A * x_gpu - b;
  RealScalar relres = r.norm() / b.norm();
  VERIFY(relres < RealScalar(1e-6));

  // Compare with CPU.
  RealScalar sol_tol = RealScalar(100) * RealScalar(n) * NumTraits<Scalar>::epsilon();
  VERIFY((x_gpu - x_cpu).norm() / (x_cpu.norm() + RealScalar(1)) < sol_tol);
}

// ---- GPU CG with Jacobi preconditioner --------------------------------------

template <typename Scalar>
void test_gpu_cg_jacobi(Index n) {
  using SpMat = SparseMatrix<Scalar, ColMajor, int>;
  using Vec = Matrix<Scalar, Dynamic, 1>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  SpMat A = make_spd<Scalar>(n);
  Vec b = Vec::Random(n);

  // CPU reference.
  ConjugateGradient<SpMat, Lower | Upper> cpu_cg;
  cpu_cg.setMaxIterations(1000);
  cpu_cg.setTolerance(RealScalar(1e-8));
  cpu_cg.compute(A);
  Vec x_cpu = cpu_cg.solve(b);

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

  // GPU CG with Jacobi preconditioner.
  gpu::Context ctx;
  gpu::Context::setThreadLocal(&ctx);
  gpu::SparseContext<Scalar> spmv_ctx(ctx);
  auto mat = spmv_ctx.deviceView(A);
  auto d_invdiag = gpu::DeviceMatrix<Scalar>::fromHost(invdiag, ctx.stream());

  auto d_b = gpu::DeviceMatrix<Scalar>::fromHost(b, ctx.stream());
  gpu::DeviceMatrix<Scalar> d_x(n, 1);
  d_x.setZero(ctx);

  gpu::DeviceMatrix<Scalar> residual(n, 1);
  residual.copyFrom(ctx, d_b);

  RealScalar rhsNorm2 = d_b.squaredNorm(ctx);
  RealScalar tol = RealScalar(1e-8);
  RealScalar threshold = tol * tol * rhsNorm2;
  RealScalar residualNorm2 = residual.squaredNorm(ctx);

  // p = precond.solve(r) = invdiag .* r
  gpu::DeviceMatrix<Scalar> p = d_invdiag.cwiseProduct(ctx, residual);
  gpu::DeviceMatrix<Scalar> z(n, 1), tmp(n, 1);

  auto absNew = residual.dot(ctx, p);
  Index maxIters = 1000;
  Index i = 0;
  while (i < maxIters) {
    tmp.noalias() = mat * p;

    auto alpha = absNew / p.dot(ctx, tmp);
    d_x += alpha * p;
    residual -= alpha * tmp;

    residualNorm2 = residual.squaredNorm(ctx);
    if (residualNorm2 < threshold) break;

    // z = precond.solve(r) = invdiag .* r
    z.cwiseProduct(ctx, d_invdiag, residual);

    auto absOld = std::move(absNew);
    absNew = residual.dot(ctx, z);
    auto beta = absNew / absOld;

    p *= beta;
    p += z;
    i++;
  }

  gpu::Context::setThreadLocal(nullptr);

  Vec x_gpu = d_x.toHost(ctx.stream());

  Vec r = A * x_gpu - b;
  RealScalar relres = r.norm() / b.norm();
  VERIFY(relres < RealScalar(1e-6));

  RealScalar sol_tol = RealScalar(100) * RealScalar(n) * NumTraits<Scalar>::epsilon();
  VERIFY((x_gpu - x_cpu).norm() / (x_cpu.norm() + RealScalar(1)) < sol_tol);
}

EIGEN_DECLARE_TEST(gpu_cg) {
  gpu_test::require_cusparse_context();

  // Split by scalar so each part compiles in parallel.
  CALL_SUBTEST_1(test_gpu_cg<double>(64));
  CALL_SUBTEST_1(test_gpu_cg<double>(256));
  CALL_SUBTEST_1(test_gpu_cg_jacobi<double>(64));
  CALL_SUBTEST_1(test_gpu_cg_jacobi<double>(256));
  CALL_SUBTEST_2(test_gpu_cg<float>(64));
  CALL_SUBTEST_2(test_gpu_cg<float>(256));
  CALL_SUBTEST_2(test_gpu_cg_jacobi<float>(64));
  CALL_SUBTEST_2(test_gpu_cg_jacobi<float>(256));
}
