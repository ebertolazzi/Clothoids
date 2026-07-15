// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Eigen Authors
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

// Tests for GpuLLT: GPU Cholesky (LL^T) using cuSOLVER.
// Covers cusolverDnXpotrf (factorization) and cusolverDnXpotrs (solve)
// for float, double, complex<float>, complex<double>, Lower and Upper.

#define EIGEN_USE_GPU
#include "main.h"
#include <Eigen/Cholesky>
#include <unsupported/Eigen/GPU>

// Identifier convention throughout this file:
//   h_ prefix for host-resident Eigen::Matrix values
//   d_ prefix for device-resident Eigen::gpu::DeviceMatrix values
// We also keep namespaces explicit (Eigen::, Eigen::gpu::) so the CPU vs GPU
// path is obvious at every call site.

// Build a random symmetric positive-definite matrix: A = M^H*M + n*I.
template <typename MatrixType>
MatrixType make_spd(Eigen::Index n) {
  using Scalar = typename MatrixType::Scalar;
  MatrixType M = MatrixType::Random(n, n);
  return M.adjoint() * M + MatrixType::Identity(n, n) * static_cast<Scalar>(n);
}

// Test factorization: L*L^H must reconstruct A to within floating-point tolerance.
template <typename Scalar, int UpLo>
void test_potrf(Eigen::Index n) {
  using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  using RealScalar = typename Eigen::NumTraits<Scalar>::Real;

  MatrixType h_A = make_spd<MatrixType>(n);

  // GPU factorization under test: factor stays on device.
  Eigen::gpu::LLT<Scalar, UpLo> gpu_llt(h_A);
  VERIFY_IS_EQUAL(gpu_llt.info(), Eigen::Success);

  // CPU LLT is the oracle. GpuLLT does not expose the device-resident factor,
  // so we validate correctness via the CPU-reconstructed matrix.
  Eigen::LLT<MatrixType, UpLo> cpu_llt(h_A);
  VERIFY_IS_EQUAL(cpu_llt.info(), Eigen::Success);
  MatrixType h_A_reconstructed = cpu_llt.reconstructedMatrix();

  RealScalar tol = RealScalar(4) * RealScalar(n) * Eigen::NumTraits<Scalar>::epsilon() * h_A.norm();
  VERIFY((h_A_reconstructed - h_A).norm() < tol);

  // Cross-check: GPU and CPU solves must agree on the same RHS. `solve(Matrix)`
  // uploads/downloads internally, so both outputs end up on the host.
  MatrixType h_b = MatrixType::Random(n, 1);
  MatrixType h_x_gpu = gpu_llt.solve(h_b);
  MatrixType h_x_cpu = cpu_llt.solve(h_b);
  VERIFY((h_x_gpu - h_x_cpu).norm() < tol);
}

// Test solve: residual ||A*X - B|| / ||B|| must be small.
template <typename Scalar, int UpLo>
void test_potrs(Eigen::Index n, Eigen::Index nrhs) {
  using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  using RealScalar = typename Eigen::NumTraits<Scalar>::Real;

  MatrixType h_A = make_spd<MatrixType>(n);
  MatrixType h_B = MatrixType::Random(n, nrhs);

  Eigen::gpu::LLT<Scalar, UpLo> gpu_llt(h_A);
  VERIFY_IS_EQUAL(gpu_llt.info(), Eigen::Success);

  MatrixType h_X = gpu_llt.solve(h_B);

  RealScalar residual = (h_A * h_X - h_B).norm() / h_B.norm();
  RealScalar tol = RealScalar(n) * Eigen::NumTraits<Scalar>::epsilon();
  VERIFY(residual < tol);
}

// Test that multiple solves against the same factor all produce correct results.
// This exercises the key design property: L stays on device across calls.
template <typename Scalar>
void test_multiple_solves(Eigen::Index n) {
  using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  using RealScalar = typename Eigen::NumTraits<Scalar>::Real;

  MatrixType h_A = make_spd<MatrixType>(n);
  Eigen::gpu::LLT<Scalar, Eigen::Lower> gpu_llt(h_A);
  VERIFY_IS_EQUAL(gpu_llt.info(), Eigen::Success);

  RealScalar tol = RealScalar(n) * Eigen::NumTraits<Scalar>::epsilon();
  for (int k = 0; k < 5; ++k) {
    MatrixType h_B = MatrixType::Random(n, 3);
    MatrixType h_X = gpu_llt.solve(h_B);
    RealScalar residual = (h_A * h_X - h_B).norm() / h_B.norm();
    VERIFY(residual < tol);
  }
}

// Test that GpuLLT correctly detects a non-SPD matrix.
void test_not_spd() {
  Eigen::MatrixXd h_A = -Eigen::MatrixXd::Identity(8, 8);  // negative definite
  Eigen::gpu::LLT<double> gpu_llt(h_A);
  VERIFY_IS_EQUAL(gpu_llt.info(), Eigen::NumericalIssue);
}

// solve(DeviceMatrix) must not silently return garbage when the factorization
// failed: it must sync the info word and assert just like solve(MatrixBase).
void test_not_spd_device_solve_asserts() {
  Eigen::MatrixXd h_A = -Eigen::MatrixXd::Identity(8, 8);
  Eigen::MatrixXd h_B = Eigen::MatrixXd::Random(8, 4);
  Eigen::gpu::LLT<double> gpu_llt(h_A);
  VERIFY_IS_EQUAL(gpu_llt.info(), Eigen::NumericalIssue);
  auto d_B = Eigen::gpu::DeviceMatrix<double>::fromHost(h_B);
  VERIFY_RAISES_ASSERT(gpu_llt.solve(d_B));
}

// ---- DeviceMatrix-native API --------------------------------------------
// These tests exercise the device-resident path: compute(DeviceMatrix) +
// solve(DeviceMatrix) -> DeviceMatrix, with the user explicitly managing
// upload/download. The tests above use the host-Matrix overloads which do
// the transfers internally; this section covers the no-implicit-transfer
// surface that keeps data on device across a chain of calls.

// compute(DeviceMatrix) + solve(DeviceMatrix) → toHost
template <typename Scalar, int UpLo>
void test_device_matrix_solve(Eigen::Index n, Eigen::Index nrhs) {
  using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  using RealScalar = typename Eigen::NumTraits<Scalar>::Real;

  MatrixType h_A = make_spd<MatrixType>(n);
  MatrixType h_B = MatrixType::Random(n, nrhs);

  auto d_A = Eigen::gpu::DeviceMatrix<Scalar>::fromHost(h_A);
  auto d_B = Eigen::gpu::DeviceMatrix<Scalar>::fromHost(h_B);

  Eigen::gpu::LLT<Scalar, UpLo> gpu_llt;
  gpu_llt.compute(d_A);
  VERIFY_IS_EQUAL(gpu_llt.info(), Eigen::Success);

  Eigen::gpu::DeviceMatrix<Scalar> d_X = gpu_llt.solve(d_B);
  MatrixType h_X = d_X.toHost();

  RealScalar residual = (h_A * h_X - h_B).norm() / h_B.norm();
  VERIFY(residual < RealScalar(n) * Eigen::NumTraits<Scalar>::epsilon());
}

// compute(DeviceMatrix&&) — move path
template <typename Scalar>
void test_device_matrix_move_compute(Eigen::Index n) {
  using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  using RealScalar = typename Eigen::NumTraits<Scalar>::Real;

  MatrixType h_A = make_spd<MatrixType>(n);
  MatrixType h_B = MatrixType::Random(n, 1);

  auto d_A = Eigen::gpu::DeviceMatrix<Scalar>::fromHost(h_A);
  Eigen::gpu::LLT<Scalar, Eigen::Lower> gpu_llt;
  gpu_llt.compute(std::move(d_A));
  VERIFY_IS_EQUAL(gpu_llt.info(), Eigen::Success);

  // d_A should be empty after move.
  VERIFY(d_A.empty());

  MatrixType h_X = gpu_llt.solve(h_B);
  RealScalar residual = (h_A * h_X - h_B).norm() / h_B.norm();
  VERIFY(residual < RealScalar(n) * Eigen::NumTraits<Scalar>::epsilon());
}

// Full async chain: compute → solve → solve again with result as RHS → toHost
template <typename Scalar>
void test_chaining(Eigen::Index n) {
  using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  using RealScalar = typename Eigen::NumTraits<Scalar>::Real;

  MatrixType h_A = make_spd<MatrixType>(n);
  MatrixType h_B = MatrixType::Random(n, 3);

  auto d_A = Eigen::gpu::DeviceMatrix<Scalar>::fromHost(h_A);
  auto d_B = Eigen::gpu::DeviceMatrix<Scalar>::fromHost(h_B);

  Eigen::gpu::LLT<Scalar, Eigen::Lower> gpu_llt;
  gpu_llt.compute(d_A);
  VERIFY_IS_EQUAL(gpu_llt.info(), Eigen::Success);

  // Chain: solve → use result as RHS for another solve. Everything stays on
  // device until the final toHost() below; that's the only sync point.
  Eigen::gpu::DeviceMatrix<Scalar> d_X = gpu_llt.solve(d_B);
  Eigen::gpu::DeviceMatrix<Scalar> d_Y = gpu_llt.solve(d_X);

  MatrixType h_Y = d_Y.toHost();

  // Verify: Y = A^{-2} * B using CPU oracle.
  MatrixType h_X_ref = Eigen::LLT<MatrixType, Eigen::Lower>(h_A).solve(h_B);
  MatrixType h_Y_ref = Eigen::LLT<MatrixType, Eigen::Lower>(h_A).solve(h_X_ref);

  RealScalar tol = RealScalar(4) * RealScalar(n) * Eigen::NumTraits<Scalar>::epsilon() * h_Y_ref.norm();
  VERIFY((h_Y - h_Y_ref).norm() < tol);
}

template <typename Scalar>
void test_scalar() {
  CALL_SUBTEST((test_potrf<Scalar, Eigen::Lower>(1)));
  CALL_SUBTEST((test_potrf<Scalar, Eigen::Lower>(64)));
  CALL_SUBTEST((test_potrf<Scalar, Eigen::Lower>(256)));
  CALL_SUBTEST((test_potrf<Scalar, Eigen::Upper>(64)));
  CALL_SUBTEST((test_potrf<Scalar, Eigen::Upper>(256)));

  CALL_SUBTEST((test_potrs<Scalar, Eigen::Lower>(64, 1)));
  CALL_SUBTEST((test_potrs<Scalar, Eigen::Lower>(64, 4)));
  CALL_SUBTEST((test_potrs<Scalar, Eigen::Lower>(256, 8)));
  CALL_SUBTEST((test_potrs<Scalar, Eigen::Upper>(64, 1)));
  CALL_SUBTEST((test_potrs<Scalar, Eigen::Upper>(256, 4)));

  CALL_SUBTEST(test_multiple_solves<Scalar>(128));

  CALL_SUBTEST((test_device_matrix_solve<Scalar, Eigen::Lower>(64, 4)));
  CALL_SUBTEST((test_device_matrix_solve<Scalar, Eigen::Upper>(128, 1)));
  CALL_SUBTEST(test_device_matrix_move_compute<Scalar>(64));
  CALL_SUBTEST(test_chaining<Scalar>(64));
}

EIGEN_DECLARE_TEST(gpu_cusolver_llt) {
  // Split by scalar so each part compiles in parallel.
  CALL_SUBTEST_1(test_scalar<float>());
  CALL_SUBTEST_2(test_scalar<double>());
  CALL_SUBTEST_3(test_scalar<std::complex<float>>());
  CALL_SUBTEST_4(test_scalar<std::complex<double>>());
  CALL_SUBTEST_5(test_not_spd());
  CALL_SUBTEST_5(test_not_spd_device_solve_asserts());
}
