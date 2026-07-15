// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Rasmus Munk Larsen <rmlarsen@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

// Tests for GpuSelfAdjointEigenSolver: GPU symmetric/Hermitian eigenvalue
// decomposition via cuSOLVER.

#define EIGEN_USE_GPU
#include "main.h"
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/GPU>

using namespace Eigen;

// ---- Reconstruction: V * diag(W) * V^H ≈ A ---------------------------------

template <typename Scalar>
void test_eigen_reconstruction(Index n) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  // Build a symmetric/Hermitian matrix.
  Mat R = Mat::Random(n, n);
  Mat A = R + R.adjoint();

  gpu::SelfAdjointEigenSolver<Scalar> es(A);
  VERIFY_IS_EQUAL(es.info(), Success);

  auto W = es.eigenvalues();
  Mat V = es.eigenvectors();

  VERIFY_IS_EQUAL(W.size(), n);
  VERIFY_IS_EQUAL(V.rows(), n);
  VERIFY_IS_EQUAL(V.cols(), n);

  // Reconstruct: A_hat = V * diag(W) * V^H.
  Mat A_hat = V * W.asDiagonal() * V.adjoint();
  RealScalar tol = RealScalar(8) * static_cast<RealScalar>(n) * NumTraits<Scalar>::epsilon() * A.norm();
  VERIFY((A_hat - A).norm() < tol);

  // Orthogonality: V^H * V ≈ I.
  Mat VhV = V.adjoint() * V;
  Mat eye = Mat::Identity(n, n);
  VERIFY((VhV - eye).norm() < tol);
}

// ---- Eigenvalues match CPU SelfAdjointEigenSolver ---------------------------

template <typename Scalar>
void test_eigen_values(Index n) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat R = Mat::Random(n, n);
  Mat A = R + R.adjoint();

  gpu::SelfAdjointEigenSolver<Scalar> gpu_es(A);
  VERIFY_IS_EQUAL(gpu_es.info(), Success);
  auto W_gpu = gpu_es.eigenvalues();

  SelfAdjointEigenSolver<Mat> cpu_es(A);
  auto W_cpu = cpu_es.eigenvalues();

  RealScalar tol =
      RealScalar(2) * static_cast<RealScalar>(n) * NumTraits<Scalar>::epsilon() * W_cpu.cwiseAbs().maxCoeff();
  VERIFY((W_gpu - W_cpu).norm() < tol);
}

// ---- Eigenvalues-only mode --------------------------------------------------

template <typename Scalar>
void test_eigen_values_only(Index n) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat R = Mat::Random(n, n);
  Mat A = R + R.adjoint();

  gpu::SelfAdjointEigenSolver<Scalar> gpu_es(A, EigenvaluesOnly);
  VERIFY_IS_EQUAL(gpu_es.info(), Success);
  auto W_gpu = gpu_es.eigenvalues();

  SelfAdjointEigenSolver<Mat> cpu_es(A, EigenvaluesOnly);
  auto W_cpu = cpu_es.eigenvalues();

  RealScalar tol =
      RealScalar(2) * static_cast<RealScalar>(n) * NumTraits<Scalar>::epsilon() * W_cpu.cwiseAbs().maxCoeff();
  VERIFY((W_gpu - W_cpu).norm() < tol);
}

// ---- DeviceMatrix input path ------------------------------------------------

template <typename Scalar>
void test_eigen_device_matrix(Index n) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat R = Mat::Random(n, n);
  Mat A = R + R.adjoint();

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  gpu::SelfAdjointEigenSolver<Scalar> es;
  es.compute(d_A);
  VERIFY_IS_EQUAL(es.info(), Success);

  auto W_gpu = es.eigenvalues();
  Mat V = es.eigenvectors();

  // Verify reconstruction.
  Mat A_hat = V * W_gpu.asDiagonal() * V.adjoint();
  RealScalar tol = RealScalar(8) * static_cast<RealScalar>(n) * NumTraits<Scalar>::epsilon() * A.norm();
  VERIFY((A_hat - A).norm() < tol);
}

// ---- Recompute (reuse solver object) ----------------------------------------

template <typename Scalar>
void test_eigen_recompute(Index n) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  gpu::SelfAdjointEigenSolver<Scalar> es;

  for (int trial = 0; trial < 3; ++trial) {
    Mat R = Mat::Random(n, n);
    Mat A = R + R.adjoint();
    es.compute(A);
    VERIFY_IS_EQUAL(es.info(), Success);

    auto W = es.eigenvalues();
    Mat V = es.eigenvectors();
    Mat A_hat = V * W.asDiagonal() * V.adjoint();
    RealScalar tol = RealScalar(8) * static_cast<RealScalar>(n) * NumTraits<Scalar>::epsilon() * A.norm();
    VERIFY((A_hat - A).norm() < tol);
  }
}

// ---- Repeated eigenvalues ---------------------------------------------------
//
// Builds A with a degenerate eigenvalue cluster. Eigenvectors within the cluster
// are not uniquely determined, but eigenvalues, reconstruction, and orthogonality
// must still hold.

template <typename Scalar>
void test_eigen_repeated_values(Index n) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using Vec = Matrix<typename NumTraits<Scalar>::Real, Dynamic, 1>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  // Build a unitary V via QR of a random matrix.
  Mat seed = Mat::Random(n, n);
  Mat V = HouseholderQR<Mat>(seed).householderQ();

  // Spectrum: a 4-fold cluster at value 1, then increasing distinct values.
  Vec eigs(n);
  for (Index i = 0; i < n; ++i) {
    eigs(i) = (i < 4) ? RealScalar(1) : RealScalar(i);
  }
  Mat A = V * eigs.asDiagonal() * V.adjoint();
  Mat A_sym = (A + A.adjoint()) * RealScalar(0.5);  // symmetrize numerically

  gpu::SelfAdjointEigenSolver<Scalar> es(A_sym);
  VERIFY_IS_EQUAL(es.info(), Success);

  // Eigenvalues must match the constructed spectrum (cuSOLVER returns ascending).
  Vec W_gpu = es.eigenvalues();
  Vec W_expected = eigs;
  std::sort(W_expected.data(), W_expected.data() + n);

  RealScalar tol_w =
      RealScalar(2) * static_cast<RealScalar>(n) * NumTraits<Scalar>::epsilon() * W_expected.cwiseAbs().maxCoeff();
  VERIFY((W_gpu - W_expected).norm() < tol_w);

  // Reconstruction must hold even with a degenerate cluster (eigenvectors form a
  // valid orthonormal basis for the cluster's invariant subspace).
  Mat V_gpu = es.eigenvectors();
  Mat A_hat = V_gpu * W_gpu.asDiagonal() * V_gpu.adjoint();
  RealScalar tol = RealScalar(20) * static_cast<RealScalar>(n) * NumTraits<Scalar>::epsilon() * A.norm();
  VERIFY((A_hat - A_sym).norm() < tol);

  // Orthogonality of computed eigenvectors must hold.
  VERIFY((V_gpu.adjoint() * V_gpu - Mat::Identity(n, n)).norm() < tol);
}

// ---- Device-side accessors --------------------------------------------------
//
// d_eigenvalues / d_eigenvectors return non-owning views over the solver's
// internal d_W_ / d_A_ buffers. Verify they agree with the host accessors and
// that repeated invocations leave the solver state intact.

template <typename Scalar>
void test_eigen_device_accessors(Index n) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;

  Mat R = Mat::Random(n, n);
  Mat A = R + R.adjoint();

  gpu::SelfAdjointEigenSolver<Scalar> es(A);
  VERIFY_IS_EQUAL(es.info(), Success);

  auto d_W = es.d_eigenvalues();
  auto d_V = es.d_eigenvectors();

  VERIFY_IS_APPROX(d_W.toHost(), es.eigenvalues());
  VERIFY_IS_APPROX(d_V.toHost(), es.eigenvectors());

  // Re-fetch must keep working (view destruction must not free the underlying buffers).
  (void)es.d_eigenvalues();
  (void)es.d_eigenvectors();
  VERIFY_IS_APPROX(es.eigenvalues(), d_W.toHost());
  VERIFY_IS_APPROX(es.eigenvectors(), d_V.toHost());
}

// ---- Chain device views into a downstream cuBLAS GEMM (no D2D copy) ---------

template <typename Scalar>
void test_eigen_chain_orthogonality(Index n) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat R = Mat::Random(n, n);
  Mat A = R + R.adjoint();

  gpu::SelfAdjointEigenSolver<Scalar> es(A);
  VERIFY_IS_EQUAL(es.info(), Success);

  // V^H * V = I — V is the device view of d_A_; the GEMM consumes the view
  // without an intervening D2D copy.
  gpu::Context ctx;
  gpu::DeviceMatrix<Scalar> d_VtV;
  {
    auto d_V = es.d_eigenvectors();
    d_VtV.device(ctx) = d_V.adjoint() * d_V;
  }
  Mat VtV = d_VtV.toHost();
  RealScalar tol = RealScalar(20) * static_cast<RealScalar>(n) * NumTraits<Scalar>::epsilon();
  VERIFY((VtV - Mat::Identity(n, n)).norm() < tol);

  // After view destruction, the solver's state must remain valid.
  Mat V_again = es.eigenvectors();
  VERIFY_IS_EQUAL(V_again.rows(), n);
}

// ---- Move support -------------------------------------------------------------

template <typename Scalar>
void test_eigen_move(Index n) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat R = Mat::Random(n, n);
  Mat A = R + R.adjoint();

  gpu::SelfAdjointEigenSolver<Scalar> es(A);
  VERIFY_IS_EQUAL(es.info(), Success);

  gpu::SelfAdjointEigenSolver<Scalar> moved(std::move(es));
  VERIFY_IS_EQUAL(moved.info(), Success);
  Mat V = moved.eigenvectors();
  auto W = moved.eigenvalues();
  RealScalar tol = RealScalar(8) * static_cast<RealScalar>(n) * NumTraits<Scalar>::epsilon() * A.norm();
  VERIFY((V * W.asDiagonal() * V.adjoint() - A).norm() < tol);

  gpu::SelfAdjointEigenSolver<Scalar> assigned;
  assigned = std::move(moved);
  VERIFY_IS_EQUAL(assigned.info(), Success);
  V = assigned.eigenvectors();
  W = assigned.eigenvalues();
  VERIFY((V * W.asDiagonal() * V.adjoint() - A).norm() < tol);
}

// ---- Empty matrix -----------------------------------------------------------

void test_eigen_empty() {
  gpu::SelfAdjointEigenSolver<double> es(MatrixXd(0, 0));
  VERIFY_IS_EQUAL(es.info(), Success);
  VERIFY_IS_EQUAL(es.rows(), 0);
  VERIFY_IS_EQUAL(es.cols(), 0);
}

// ---- Per-scalar driver ------------------------------------------------------

template <typename Scalar>
void test_scalar() {
  // Reconstruction + orthogonality.
  CALL_SUBTEST(test_eigen_reconstruction<Scalar>(64));
  CALL_SUBTEST(test_eigen_reconstruction<Scalar>(128));

  // Eigenvalues match CPU.
  CALL_SUBTEST(test_eigen_values<Scalar>(64));
  CALL_SUBTEST(test_eigen_values<Scalar>(128));

  // Values-only mode.
  CALL_SUBTEST(test_eigen_values_only<Scalar>(64));

  // DeviceMatrix input.
  CALL_SUBTEST(test_eigen_device_matrix<Scalar>(64));

  // Recompute.
  CALL_SUBTEST(test_eigen_recompute<Scalar>(32));

  // Repeated eigenvalues (degenerate cluster).
  CALL_SUBTEST(test_eigen_repeated_values<Scalar>(64));

  // Device accessors.
  CALL_SUBTEST(test_eigen_device_accessors<Scalar>(64));

  // Chain device view into a downstream GEMM.
  CALL_SUBTEST(test_eigen_chain_orthogonality<Scalar>(64));

  // Move constructor/assignment.
  CALL_SUBTEST(test_eigen_move<Scalar>(32));
}

EIGEN_DECLARE_TEST(gpu_cusolver_eigen) {
  // Split by scalar so each part compiles in parallel.
  CALL_SUBTEST_1(test_scalar<float>());
  CALL_SUBTEST_2(test_scalar<double>());
  CALL_SUBTEST_3(test_scalar<std::complex<float>>());
  CALL_SUBTEST_4(test_scalar<std::complex<double>>());
  CALL_SUBTEST_5(test_eigen_empty());
}
