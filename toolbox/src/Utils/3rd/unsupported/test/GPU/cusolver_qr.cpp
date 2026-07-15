// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Rasmus Munk Larsen <rmlarsen@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

// Tests for GpuQR: GPU QR decomposition via cuSOLVER.

#define EIGEN_USE_GPU
#include "main.h"
#include <Eigen/QR>
#include <Eigen/SVD>
#include <unsupported/Eigen/GPU>

using namespace Eigen;

// ---- Solve square system: A * X = B -----------------------------------------

template <typename Scalar>
void test_qr_solve_square(Index n, Index nrhs) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(n, n);
  Mat B = Mat::Random(n, nrhs);

  gpu::QR<Scalar> qr(A);
  VERIFY_IS_EQUAL(qr.info(), Success);

  Mat X = qr.solve(B);
  RealScalar residual = (A * X - B).norm() / (A.norm() * X.norm());
  VERIFY(residual < RealScalar(10) * RealScalar(n) * NumTraits<Scalar>::epsilon());
}

// ---- Solve overdetermined system: m > n (least-squares) ---------------------

template <typename Scalar>
void test_qr_solve_overdetermined(Index m, Index n, Index nrhs) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  eigen_assert(m >= n);
  Mat A = Mat::Random(m, n);
  Mat B = Mat::Random(m, nrhs);

  gpu::QR<Scalar> qr(A);
  VERIFY_IS_EQUAL(qr.info(), Success);

  Mat X = qr.solve(B);
  VERIFY_IS_EQUAL(X.rows(), n);
  VERIFY_IS_EQUAL(X.cols(), nrhs);

  // For an overdetermined system the residual r = A X - B is generally large
  // (O(||B||) when B is not in col(A)). What backward-stable least-squares
  // makes small is the gradient of the LS objective, A^H r, which is zero at
  // the optimum. Higham's bound for QR-based LS gives ||A^H r|| <= O(m * eps)
  // * ||A|| * (||A|| ||X|| + ||B||), regardless of kappa(A).
  Mat X_cpu = HouseholderQR<Mat>(A).solve(B);
  RealScalar tol = RealScalar(10) * RealScalar(m) * NumTraits<Scalar>::epsilon();
  RealScalar A_norm = A.norm();
  RealScalar denom_gpu = A_norm * (A_norm * X.norm() + B.norm());
  RealScalar denom_cpu = A_norm * (A_norm * X_cpu.norm() + B.norm());
  VERIFY((A.adjoint() * (A * X - B)).norm() / denom_gpu < tol);
  VERIFY((A.adjoint() * (A * X_cpu - B)).norm() / denom_cpu < tol);
}

// ---- Solve with DeviceMatrix input ------------------------------------------

template <typename Scalar>
void test_qr_solve_device(Index n, Index nrhs) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(n, n);
  Mat B = Mat::Random(n, nrhs);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::QR<Scalar> qr;
  qr.compute(d_A);
  VERIFY_IS_EQUAL(qr.info(), Success);

  gpu::DeviceMatrix<Scalar> d_X = qr.solve(d_B);
  Mat X = d_X.toHost();

  RealScalar residual = (A * X - B).norm() / (A.norm() * X.norm());
  VERIFY(residual < RealScalar(10) * RealScalar(n) * NumTraits<Scalar>::epsilon());
}

// ---- Solve overdetermined via device path -----------------------------------

template <typename Scalar>
void test_qr_solve_overdetermined_device(Index m, Index n, Index nrhs) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  eigen_assert(m >= n);
  Mat A = Mat::Random(m, n);
  Mat B = Mat::Random(m, nrhs);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::QR<Scalar> qr;
  qr.compute(d_A);
  VERIFY_IS_EQUAL(qr.info(), Success);

  gpu::DeviceMatrix<Scalar> d_X = qr.solve(d_B);
  VERIFY_IS_EQUAL(d_X.rows(), n);
  VERIFY_IS_EQUAL(d_X.cols(), nrhs);

  Mat X = d_X.toHost();
  Mat X_cpu = HouseholderQR<Mat>(A).solve(B);
  // See test_qr_solve_overdetermined: backward error for LS is bounded on
  // the gradient A^H r, not on r itself.
  RealScalar tol = RealScalar(10) * RealScalar(m) * NumTraits<Scalar>::epsilon();
  RealScalar A_norm = A.norm();
  RealScalar denom_gpu = A_norm * (A_norm * X.norm() + B.norm());
  RealScalar denom_cpu = A_norm * (A_norm * X_cpu.norm() + B.norm());
  VERIFY((A.adjoint() * (A * X - B)).norm() / denom_gpu < tol);
  VERIFY((A.adjoint() * (A * X_cpu - B)).norm() / denom_cpu < tol);
}

// ---- Solve underdetermined system: m < n (minimum norm) ---------------------

template <typename Scalar>
void test_qr_solve_underdetermined(Index m, Index n, Index nrhs) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  eigen_assert(m < n);
  Mat A = Mat::Random(m, n);
  Mat B = Mat::Random(m, nrhs);

  gpu::QR<Scalar> qr(A);
  VERIFY_IS_EQUAL(qr.info(), Success);

  Mat X = qr.solve(B);
  VERIFY_IS_EQUAL(X.rows(), n);
  VERIFY_IS_EQUAL(X.cols(), nrhs);

  // The minimum-norm solution exactly satisfies A X = B (m equations, n > m
  // unknowns). Backward error is O(m * eps) * ||A|| * ||X||.
  RealScalar tol = RealScalar(20) * RealScalar(n) * NumTraits<Scalar>::epsilon();
  VERIFY((A * X - B).norm() / (A.norm() * X.norm() + B.norm()) < tol);

  // Min-norm property: X ⟂ null(A), so X ∈ row space(A) = col space(A^H).
  // ||X|| <= ||X_any|| for any other solution. Compare to the SVD min-norm reference.
  Mat X_svd = BDCSVD<Mat, ComputeThinU | ComputeThinV>(A).solve(B);
  // Both should have the same norm (up to rounding).
  VERIFY(numext::abs(X.norm() - X_svd.norm()) / X_svd.norm() <
         RealScalar(50) * RealScalar(n) * NumTraits<Scalar>::epsilon());
}

template <typename Scalar>
void test_qr_solve_underdetermined_device(Index m, Index n, Index nrhs) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  eigen_assert(m < n);
  Mat A = Mat::Random(m, n);
  Mat B = Mat::Random(m, nrhs);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::QR<Scalar> qr;
  qr.compute(d_A);
  VERIFY_IS_EQUAL(qr.info(), Success);

  gpu::DeviceMatrix<Scalar> d_X = qr.solve(d_B);
  VERIFY_IS_EQUAL(d_X.rows(), n);
  VERIFY_IS_EQUAL(d_X.cols(), nrhs);

  Mat X = d_X.toHost();
  RealScalar tol = RealScalar(20) * RealScalar(n) * NumTraits<Scalar>::epsilon();
  VERIFY((A * X - B).norm() / (A.norm() * X.norm() + B.norm()) < tol);
}

// ---- matrixR() returns the upper-triangular factor --------------------------

template <typename Scalar>
void test_qr_matrixR(Index m, Index n) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  eigen_assert(m >= n);
  Mat A = Mat::Random(m, n);
  gpu::QR<Scalar> qr(A);
  VERIFY_IS_EQUAL(qr.info(), Success);

  Mat R = qr.matrixR();
  VERIFY_IS_EQUAL(R.rows(), n);
  VERIFY_IS_EQUAL(R.cols(), n);

  // Verify R is upper-triangular (lower triangle exactly zero).
  for (Index j = 0; j < n; ++j) {
    for (Index i = j + 1; i < n; ++i) {
      VERIFY_IS_EQUAL(R(i, j), Scalar(0));
    }
  }

  // Verify that the diagonal of R matches the absolute values of CPU R diagonal
  // up to sign (Householder QR has free diagonal sign).
  Mat R_cpu = HouseholderQR<Mat>(A).matrixQR().topRows(n).template triangularView<Upper>();
  RealScalar tol = RealScalar(20) * RealScalar(m) * NumTraits<Scalar>::epsilon() * A.norm();
  for (Index i = 0; i < n; ++i) {
    VERIFY(numext::abs(numext::abs(R(i, i)) - numext::abs(R_cpu(i, i))) < tol);
  }
}

// ---- Multiple solves reuse the factorization --------------------------------

template <typename Scalar>
void test_qr_multiple_solves(Index n) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(n, n);
  gpu::QR<Scalar> qr(A);
  VERIFY_IS_EQUAL(qr.info(), Success);

  RealScalar tol = RealScalar(10) * RealScalar(n) * NumTraits<Scalar>::epsilon();
  for (int k = 0; k < 5; ++k) {
    Mat B = Mat::Random(n, 3);
    Mat X = qr.solve(B);
    RealScalar residual = (A * X - B).norm() / (A.norm() * X.norm());
    VERIFY(residual < tol);
  }
}

// ---- Agreement with CPU HouseholderQR ---------------------------------------

template <typename Scalar>
void test_qr_vs_cpu(Index n, Index nrhs) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(n, n);
  Mat B = Mat::Random(n, nrhs);

  gpu::QR<Scalar> gpu_qr(A);
  VERIFY_IS_EQUAL(gpu_qr.info(), Success);

  Mat X_gpu = gpu_qr.solve(B);
  Mat X_cpu = HouseholderQR<Mat>(A).solve(B);

  // Compare via residual rather than directly between X_gpu and X_cpu: for an
  // ill-conditioned A (a random N(0,1) n*n matrix easily reaches kappa(A) ~ n),
  // the forward error of each correct solve is bounded by O(kappa * eps), so
  // ||X_gpu - X_cpu|| / ||X_cpu|| can legitimately exceed n*eps even though
  // both are correct. The relative residual ||A X - B|| / (||A||*||X|| + ||B||)
  // is bounded by Higham's backward-stable solve result at O(n * eps),
  // independent of kappa.
  RealScalar tol = RealScalar(10) * RealScalar(n) * NumTraits<Scalar>::epsilon();
  RealScalar A_norm = A.norm();
  RealScalar B_norm = B.norm();
  RealScalar denom_gpu = A_norm * X_gpu.norm() + B_norm;
  RealScalar denom_cpu = A_norm * X_cpu.norm() + B_norm;
  VERIFY((A * X_gpu - B).norm() / denom_gpu < tol);
  VERIFY((A * X_cpu - B).norm() / denom_cpu < tol);
}

// ---- Per-scalar driver ------------------------------------------------------

template <typename Scalar>
void test_scalar() {
  CALL_SUBTEST(test_qr_solve_square<Scalar>(1, 1));
  CALL_SUBTEST(test_qr_solve_square<Scalar>(64, 1));
  CALL_SUBTEST(test_qr_solve_square<Scalar>(64, 4));
  CALL_SUBTEST(test_qr_solve_square<Scalar>(256, 8));

  CALL_SUBTEST(test_qr_solve_overdetermined<Scalar>(128, 64, 4));
  CALL_SUBTEST(test_qr_solve_overdetermined<Scalar>(256, 128, 1));

  CALL_SUBTEST(test_qr_solve_underdetermined<Scalar>(64, 128, 4));
  CALL_SUBTEST(test_qr_solve_underdetermined<Scalar>(64, 128, 1));
  CALL_SUBTEST(test_qr_solve_underdetermined_device<Scalar>(64, 128, 4));

  CALL_SUBTEST(test_qr_matrixR<Scalar>(64, 64));
  CALL_SUBTEST(test_qr_matrixR<Scalar>(128, 64));

  CALL_SUBTEST(test_qr_solve_device<Scalar>(64, 4));
  CALL_SUBTEST(test_qr_solve_overdetermined_device<Scalar>(128, 64, 4));
  CALL_SUBTEST(test_qr_multiple_solves<Scalar>(64));
  CALL_SUBTEST(test_qr_vs_cpu<Scalar>(64, 4));
  CALL_SUBTEST(test_qr_vs_cpu<Scalar>(256, 8));
}

void test_qr_empty() {
  gpu::QR<double> qr(MatrixXd(0, 0));
  VERIFY_IS_EQUAL(qr.info(), Success);
  VERIFY_IS_EQUAL(qr.rows(), 0);
  VERIFY_IS_EQUAL(qr.cols(), 0);
}

EIGEN_DECLARE_TEST(gpu_cusolver_qr) {
  // Split by scalar so each part compiles in parallel.
  CALL_SUBTEST_1(test_scalar<float>());
  CALL_SUBTEST_2(test_scalar<double>());
  CALL_SUBTEST_3(test_scalar<std::complex<float>>());
  CALL_SUBTEST_4(test_scalar<std::complex<double>>());
  CALL_SUBTEST_5(test_qr_empty());
}
