// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Rasmus Munk Larsen <rmlarsen@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

// Tests for GpuSparseLLT: GPU sparse Cholesky via cuDSS.

#define EIGEN_USE_GPU
#include "main.h"
#include <Eigen/Sparse>
#include <unsupported/Eigen/GPU>
#include "gpu_test_helpers.h"

using namespace Eigen;

// ---- Helper: build a random sparse SPD matrix -------------------------------

template <typename Scalar>
SparseMatrix<Scalar, ColMajor, int> make_spd(Index n, double density = 0.1) {
  using SpMat = SparseMatrix<Scalar, ColMajor, int>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  // Off-diagonal entries carry nonzero imaginary parts for complex Scalar so
  // the cuDSS HPD path (CUDSS_MTYPE_HPD) is genuinely exercised; A = R^H * R
  // is Hermitian regardless of R, and the +nI shift keeps the diagonal real.
  SpMat R(n, n);
  R.reserve(VectorXi::Constant(n, static_cast<int>(n * density) + 1));
  for (Index j = 0; j < n; ++j) {
    for (Index i = 0; i < n; ++i) {
      if (i == j || (std::rand() / double(RAND_MAX)) < density) {
        const RealScalar re = RealScalar(std::rand() / double(RAND_MAX) - 0.5);
        const RealScalar im = RealScalar(std::rand() / double(RAND_MAX) - 0.5);
        R.insert(i, j) = gpu_test::make_test_value<Scalar>(re, im);
      }
    }
  }
  R.makeCompressed();

  // A = R^H * R + n * I  (guaranteed SPD).
  SpMat A = R.adjoint() * R;
  for (Index i = 0; i < n; ++i) A.coeffRef(i, i) += Scalar(RealScalar(n));
  A.makeCompressed();
  return A;
}

// ---- Solve and check residual -----------------------------------------------

template <typename Scalar>
void test_solve(Index n) {
  using SpMat = SparseMatrix<Scalar, ColMajor, int>;
  using Vec = Matrix<Scalar, Dynamic, 1>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  SpMat A = make_spd<Scalar>(n);
  Vec b = Vec::Random(n);

  gpu::SparseLLT<Scalar> llt(A);
  VERIFY_IS_EQUAL(llt.info(), Success);

  Vec x = llt.solve(b);
  VERIFY_IS_EQUAL(x.rows(), n);

  // Check residual: ||Ax - b|| / ||b||.
  Vec r = A * x - b;
  RealScalar tol = RealScalar(100) * RealScalar(n) * NumTraits<Scalar>::epsilon();
  VERIFY(r.norm() / b.norm() < tol);
}

// ---- Compare with CPU SimplicialLLT -----------------------------------------

template <typename Scalar>
void test_vs_cpu(Index n) {
  using SpMat = SparseMatrix<Scalar, ColMajor, int>;
  using Vec = Matrix<Scalar, Dynamic, 1>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  SpMat A = make_spd<Scalar>(n);
  Vec b = Vec::Random(n);

  gpu::SparseLLT<Scalar> gpu_llt(A);
  VERIFY_IS_EQUAL(gpu_llt.info(), Success);
  Vec x_gpu = gpu_llt.solve(b);

  SimplicialLLT<SpMat> cpu_llt(A);
  VERIFY_IS_EQUAL(cpu_llt.info(), Success);
  Vec x_cpu = cpu_llt.solve(b);

  RealScalar tol = RealScalar(100) * RealScalar(n) * NumTraits<Scalar>::epsilon();
  VERIFY((x_gpu - x_cpu).norm() / x_cpu.norm() < tol);
}

// ---- Multiple RHS -----------------------------------------------------------

template <typename Scalar>
void test_multiple_rhs(Index n, Index nrhs) {
  using SpMat = SparseMatrix<Scalar, ColMajor, int>;
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  SpMat A = make_spd<Scalar>(n);
  Mat B = Mat::Random(n, nrhs);

  gpu::SparseLLT<Scalar> llt(A);
  VERIFY_IS_EQUAL(llt.info(), Success);

  Mat X = llt.solve(B);
  VERIFY_IS_EQUAL(X.rows(), n);
  VERIFY_IS_EQUAL(X.cols(), nrhs);

  Mat R = A * X - B;
  RealScalar tol = RealScalar(100) * RealScalar(n) * NumTraits<Scalar>::epsilon();
  VERIFY(R.norm() / B.norm() < tol);
}

// ---- Separate analyze + factorize (refactorization) -------------------------

template <typename Scalar>
void test_refactorize(Index n) {
  using SpMat = SparseMatrix<Scalar, ColMajor, int>;
  using Vec = Matrix<Scalar, Dynamic, 1>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  SpMat A = make_spd<Scalar>(n);
  Vec b = Vec::Random(n);

  gpu::SparseLLT<Scalar> llt;
  llt.analyzePattern(A);
  VERIFY_IS_EQUAL(llt.info(), Success);

  // First factorize + solve.
  llt.factorize(A);
  VERIFY_IS_EQUAL(llt.info(), Success);
  Vec x1 = llt.solve(b);

  // Modify values (keep same pattern): scale diagonal.
  SpMat A2 = A;
  for (Index i = 0; i < n; ++i) A2.coeffRef(i, i) *= Scalar(RealScalar(2));

  // Refactorize with same pattern.
  llt.factorize(A2);
  VERIFY_IS_EQUAL(llt.info(), Success);
  Vec x2 = llt.solve(b);

  // Both solutions should satisfy their respective systems.
  RealScalar tol = RealScalar(100) * RealScalar(n) * NumTraits<Scalar>::epsilon();
  VERIFY((A * x1 - b).norm() / b.norm() < tol);
  VERIFY((A2 * x2 - b).norm() / b.norm() < tol);

  // Solutions should differ meaningfully: diagonal was scaled 2x, so x1 vs x2
  // should differ by a substantial fraction of x1's magnitude. A near-epsilon
  // bound here would pass even if refactorize silently reused the stale factor.
  VERIFY((x1 - x2).norm() > RealScalar(0.01) * x1.norm());
}

// ---- Empty matrix -----------------------------------------------------------

template <typename Scalar>
void test_empty() {
  using SpMat = SparseMatrix<Scalar, ColMajor, int>;
  SpMat A(0, 0);
  A.makeCompressed();
  gpu::SparseLLT<Scalar> llt(A);
  VERIFY_IS_EQUAL(llt.info(), Success);
  VERIFY_IS_EQUAL(llt.rows(), 0);
  VERIFY_IS_EQUAL(llt.cols(), 0);
}

// ---- Upper triangle ---------------------------------------------------------

template <typename Scalar>
void test_upper(Index n) {
  using SpMat = SparseMatrix<Scalar, ColMajor, int>;
  using Vec = Matrix<Scalar, Dynamic, 1>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  SpMat A = make_spd<Scalar>(n);
  Vec b = Vec::Random(n);

  gpu::SparseLLT<Scalar, Upper> llt(A);
  VERIFY_IS_EQUAL(llt.info(), Success);

  Vec x = llt.solve(b);
  Vec r = A * x - b;
  RealScalar tol = RealScalar(100) * RealScalar(n) * NumTraits<Scalar>::epsilon();
  VERIFY(r.norm() / b.norm() < tol);
}

// ---- Per-scalar driver ------------------------------------------------------

template <typename Scalar>
void test_scalar() {
  CALL_SUBTEST(test_solve<Scalar>(64));
  CALL_SUBTEST(test_solve<Scalar>(256));
  CALL_SUBTEST(test_vs_cpu<Scalar>(64));
  CALL_SUBTEST(test_multiple_rhs<Scalar>(64, 4));
  CALL_SUBTEST(test_refactorize<Scalar>(64));
  CALL_SUBTEST(test_upper<Scalar>(64));
}

EIGEN_DECLARE_TEST(gpu_cudss_llt) {
  gpu_test::require_cudss_context();
  // Split by scalar so each part compiles in parallel.
  CALL_SUBTEST_1(test_scalar<float>());
  CALL_SUBTEST_2(test_scalar<double>());
  CALL_SUBTEST_3(test_scalar<std::complex<float>>());
  CALL_SUBTEST_4(test_scalar<std::complex<double>>());
  CALL_SUBTEST_5(test_empty<float>());
  CALL_SUBTEST_5(test_empty<double>());
  CALL_SUBTEST_5(test_empty<std::complex<float>>());
  CALL_SUBTEST_5(test_empty<std::complex<double>>());
}
