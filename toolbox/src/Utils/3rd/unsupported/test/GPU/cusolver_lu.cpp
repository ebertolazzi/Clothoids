// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Eigen Authors
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

// Tests for GpuLU: GPU partial-pivoting LU decomposition via cuSOLVER.
// Covers cusolverDnXgetrf (factorization) and cusolverDnXgetrs (solve)
// for float, double, complex<float>, complex<double>.
//
#define EIGEN_USE_GPU
#include "main.h"
#include <unsupported/Eigen/GPU>

using namespace Eigen;

// ---- Test factorization + NoTrans solve: residual ||A*X - B|| / ||B|| -------

template <typename Scalar>
void test_getrf(Index n) {
  using MatrixType = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  MatrixType A = MatrixType::Random(n, n);
  MatrixType B = MatrixType::Random(n, 4);

  gpu::LU<Scalar> lu(A);
  VERIFY_IS_EQUAL(lu.info(), Success);

  MatrixType X = lu.solve(B);
  // Backward error bound for LU: ||A*X - B|| <= O(n*u) * ||A|| * ||X||.
  // Normalize by ||A||*||X|| rather than ||B|| to be condition-number agnostic.
  RealScalar residual = (A * X - B).norm() / (A.norm() * X.norm());
  VERIFY(residual < RealScalar(10) * RealScalar(n) * NumTraits<Scalar>::epsilon());
}

// ---- Test solve: A^T*X = B and A^H*X = B ------------------------------------

template <typename Scalar>
void test_getrs_trans(Index n) {
  using MatrixType = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  MatrixType A = MatrixType::Random(n, n);
  MatrixType B = MatrixType::Random(n, 3);
  RealScalar tol = RealScalar(10) * RealScalar(n) * NumTraits<Scalar>::epsilon();

  gpu::LU<Scalar> lu(A);
  VERIFY_IS_EQUAL(lu.info(), Success);

  MatrixType Xt = lu.solve(B, gpu::GpuOp::Trans);
  VERIFY((A.transpose() * Xt - B).norm() / (A.norm() * Xt.norm()) < tol);

  MatrixType Xc = lu.solve(B, gpu::GpuOp::ConjTrans);
  VERIFY((A.adjoint() * Xc - B).norm() / (A.norm() * Xc.norm()) < tol);
}

// ---- Test multiple solves reuse the device-resident LU ----------------------

template <typename Scalar>
void test_multiple_solves(Index n) {
  using MatrixType = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  MatrixType A = MatrixType::Random(n, n);
  gpu::LU<Scalar> lu(A);
  VERIFY_IS_EQUAL(lu.info(), Success);

  RealScalar tol = RealScalar(10) * RealScalar(n) * NumTraits<Scalar>::epsilon();
  for (int k = 0; k < 5; ++k) {
    MatrixType B = MatrixType::Random(n, 3);
    MatrixType X = lu.solve(B);
    VERIFY((A * X - B).norm() / (A.norm() * X.norm()) < tol);
  }
}

// ---- Residual check for host solve ------------------------------------------

template <typename Scalar>
void test_vs_cpu(Index n) {
  using MatrixType = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  MatrixType A = MatrixType::Random(n, n);
  MatrixType B = MatrixType::Random(n, 5);

  gpu::LU<Scalar> gpu_lu(A);
  VERIFY_IS_EQUAL(gpu_lu.info(), Success);

  MatrixType X_gpu = gpu_lu.solve(B);

  RealScalar residual = (A * X_gpu - B).norm() / (A.norm() * X_gpu.norm());
  VERIFY(residual < RealScalar(10) * RealScalar(n) * NumTraits<Scalar>::epsilon());
}

// ---- Singular matrix detection ----------------------------------------------

void test_singular() {
  MatrixXd A = MatrixXd::Zero(8, 8);
  gpu::LU<double> lu(A);
  VERIFY_IS_EQUAL(lu.info(), NumericalIssue);
}

// solve(DeviceMatrix) must not silently return garbage when the factorization
// failed: it must sync the info word and assert just like solve(MatrixBase).
void test_singular_device_solve_asserts() {
  MatrixXd A = MatrixXd::Zero(8, 8);
  MatrixXd B = MatrixXd::Random(8, 4);
  gpu::LU<double> lu(A);
  VERIFY_IS_EQUAL(lu.info(), NumericalIssue);
  auto d_B = gpu::DeviceMatrix<double>::fromHost(B);
  VERIFY_RAISES_ASSERT(lu.solve(d_B));
}

// ---- DeviceMatrix integration tests -----------------------------------------

template <typename Scalar>
void test_device_matrix_solve(Index n) {
  using MatrixType = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  MatrixType A = MatrixType::Random(n, n);
  MatrixType B = MatrixType::Random(n, 4);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::LU<Scalar> lu;
  lu.compute(d_A);
  VERIFY_IS_EQUAL(lu.info(), Success);

  gpu::DeviceMatrix<Scalar> d_X = lu.solve(d_B);
  MatrixType X = d_X.toHost();

  RealScalar residual = (A * X - B).norm() / (A.norm() * X.norm());
  VERIFY(residual < RealScalar(10) * RealScalar(n) * NumTraits<Scalar>::epsilon());
}

template <typename Scalar>
void test_device_matrix_move_compute(Index n) {
  using MatrixType = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  MatrixType A = MatrixType::Random(n, n);
  MatrixType B = MatrixType::Random(n, 1);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  gpu::LU<Scalar> lu;
  lu.compute(std::move(d_A));
  VERIFY_IS_EQUAL(lu.info(), Success);
  VERIFY(d_A.empty());

  MatrixType X = lu.solve(B);
  RealScalar residual = (A * X - B).norm() / (A.norm() * X.norm());
  VERIFY(residual < RealScalar(10) * RealScalar(n) * NumTraits<Scalar>::epsilon());
}

template <typename Scalar>
void test_chaining(Index n) {
  using MatrixType = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  MatrixType A = MatrixType::Random(n, n);
  MatrixType B = MatrixType::Random(n, 3);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::LU<Scalar> lu;
  lu.compute(d_A);
  VERIFY_IS_EQUAL(lu.info(), Success);

  // Chain: solve → use result as RHS
  gpu::DeviceMatrix<Scalar> d_X = lu.solve(d_B);
  gpu::DeviceMatrix<Scalar> d_Y = lu.solve(d_X);
  MatrixType X = d_X.toHost();
  MatrixType Y = d_Y.toHost();

  RealScalar tol = RealScalar(10) * RealScalar(n) * NumTraits<Scalar>::epsilon();
  VERIFY((A * X - B).norm() / (A.norm() * X.norm()) < tol);
  VERIFY((A * Y - X).norm() / (A.norm() * Y.norm()) < tol);
}

// ---- Per-scalar driver -------------------------------------------------------

template <typename Scalar>
void test_scalar() {
  CALL_SUBTEST(test_getrf<Scalar>(1));
  CALL_SUBTEST(test_getrf<Scalar>(64));
  CALL_SUBTEST(test_getrf<Scalar>(256));

  CALL_SUBTEST(test_getrs_trans<Scalar>(64));
  CALL_SUBTEST(test_getrs_trans<Scalar>(128));

  CALL_SUBTEST(test_multiple_solves<Scalar>(128));

  CALL_SUBTEST(test_vs_cpu<Scalar>(64));
  CALL_SUBTEST(test_vs_cpu<Scalar>(256));

  CALL_SUBTEST(test_device_matrix_solve<Scalar>(64));
  CALL_SUBTEST(test_device_matrix_move_compute<Scalar>(64));
  CALL_SUBTEST(test_chaining<Scalar>(64));
}

EIGEN_DECLARE_TEST(gpu_cusolver_lu) {
  // Split by scalar so each part compiles in parallel.
  CALL_SUBTEST_1(test_scalar<float>());
  CALL_SUBTEST_2(test_scalar<double>());
  CALL_SUBTEST_3(test_scalar<std::complex<float>>());
  CALL_SUBTEST_4(test_scalar<std::complex<double>>());
  CALL_SUBTEST_5(test_singular());
  CALL_SUBTEST_5(test_singular_device_solve_asserts());
}
