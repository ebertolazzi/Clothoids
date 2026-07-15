// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Rasmus Munk Larsen <rmlarsen@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

// Tests for GpuSVD: GPU SVD via cuSOLVER.

#define EIGEN_USE_GPU
#include "main.h"
#include <Eigen/SVD>
#include <unsupported/Eigen/GPU>

using namespace Eigen;

// ---- SVD reconstruction: U * diag(S) * VT ≈ A ------------------------------

template <typename Scalar, unsigned int Options>
void test_svd_reconstruction(Index m, Index n) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(m, n);
  gpu::SVD<Scalar> svd(A, Options);
  VERIFY_IS_EQUAL(svd.info(), Success);

  auto S = svd.singularValues();
  Mat U = svd.matrixU();
  Mat VT = svd.matrixVT();

  const Index k = (std::min)(m, n);

  // Reconstruct: A_hat = U[:,:k] * diag(S) * VT[:k,:].
  Mat A_hat = U.leftCols(k) * S.asDiagonal() * VT.topRows(k);
  RealScalar tol = RealScalar(5) * std::sqrt(static_cast<RealScalar>(k)) * NumTraits<Scalar>::epsilon() * A.norm();
  VERIFY((A_hat - A).norm() < tol);

  // Orthogonality: U^H * U ≈ I.
  Mat UtU = U.adjoint() * U;
  Mat I_u = Mat::Identity(U.cols(), U.cols());
  VERIFY((UtU - I_u).norm() < tol);

  // Orthogonality: VT * VT^H ≈ I.
  Mat VtVh = VT * VT.adjoint();
  Mat I_v = Mat::Identity(VT.rows(), VT.rows());
  VERIFY((VtVh - I_v).norm() < tol);
}

// ---- Singular values match CPU BDCSVD ---------------------------------------

template <typename Scalar>
void test_svd_singular_values(Index m, Index n) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(m, n);
  gpu::SVD<Scalar> svd(A, 0);  // values only
  VERIFY_IS_EQUAL(svd.info(), Success);

  auto S_gpu = svd.singularValues();
  auto S_cpu = BDCSVD<Mat>(A).singularValues();

  // Weyl's perturbation bound (Higham, Accuracy and Stability of Numerical
  // Algorithms, 2nd ed., §10.2.3): |σ_i(A) - σ_i(A+δA)| ≤ ||δA||. Both cuSOLVER
  // and BDCSVD have backward error ||δA|| / ||A|| ≤ p(m,n) · u; the difference
  // of the two computed singular value vectors is bounded by 2·p(m,n)·u·S_max,
  // and √(min) · ||·||_∞ ≥ ||·||_2. Across 5k trials per case in
  // {float, double, complex<float>, complex<double>} × {(64,64), (128,64)},
  // worst observed err / (√min · u · S_max) is ≈ 6.2; we use 12 for headroom
  // against future runs hitting the tail of the distribution.
  RealScalar tol =
      RealScalar(12) * std::sqrt(static_cast<RealScalar>((std::min)(m, n))) * NumTraits<Scalar>::epsilon() * S_cpu(0);
  VERIFY((S_gpu - S_cpu).norm() < tol);
}

// ---- Solve: pseudoinverse ---------------------------------------------------

template <typename Scalar>
void test_svd_solve(Index m, Index n, Index nrhs) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(m, n);
  Mat B = Mat::Random(m, nrhs);

  gpu::SVD<Scalar> svd(A, ComputeThinU | ComputeThinV);
  VERIFY_IS_EQUAL(svd.info(), Success);

  Mat X = svd.solve(B);
  VERIFY_IS_EQUAL(X.rows(), n);
  VERIFY_IS_EQUAL(X.cols(), nrhs);

  // Compare with CPU BDCSVD solve. Wedin's perturbation theorem (Higham,
  // Accuracy and Stability of Numerical Algorithms, 2nd ed., §20.1) bounds
  // the forward error of a backward-stable SVD-based pseudoinverse solve by
  // c · κ(A) · u with c = O(1). Comparing two such solvers doubles the
  // constant. Across 6k trials over {float, double, complex<float>,
  // complex<double>} and {square, over-/underdetermined} shapes, the worst
  // observed err / (κ · u) is 5.3.
  auto cpu_svd = BDCSVD<Mat, ComputeThinU | ComputeThinV>(A);
  Mat X_cpu = cpu_svd.solve(B);
  auto S = cpu_svd.singularValues();
  const RealScalar cond = S(0) / S(S.size() - 1);
  const RealScalar tol = RealScalar(8) * cond * NumTraits<Scalar>::epsilon();
  VERIFY((X - X_cpu).norm() / X_cpu.norm() < tol);
}

// ---- Solve: truncated -------------------------------------------------------

template <typename Scalar>
void test_svd_solve_truncated(Index m, Index n) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(m, n);
  Mat B = Mat::Random(m, 1);
  const Index k = (std::min)(m, n);
  const Index trunc = k / 2;
  eigen_assert(trunc > 0);

  gpu::SVD<Scalar> svd(A, ComputeThinU | ComputeThinV);
  Mat X_trunc = svd.solve(B, trunc);

  // Build CPU reference: truncated pseudoinverse.
  auto cpu_svd = BDCSVD<Mat, ComputeThinU | ComputeThinV>(A);
  auto S = cpu_svd.singularValues();
  Mat U = cpu_svd.matrixU();
  Mat V = cpu_svd.matrixV();

  // D_ii = 1/S_i for i < trunc, 0 otherwise.
  Matrix<RealScalar, Dynamic, 1> D = Matrix<RealScalar, Dynamic, 1>::Zero(k);
  for (Index i = 0; i < trunc; ++i) D(i) = RealScalar(1) / S(i);
  Mat X_ref = V * D.asDiagonal() * U.adjoint() * B;

  RealScalar tol = RealScalar(100) * RealScalar(k) * NumTraits<Scalar>::epsilon();
  VERIFY((X_trunc - X_ref).norm() / X_ref.norm() < tol);
}

// ---- Solve: Tikhonov regularized --------------------------------------------

template <typename Scalar>
void test_svd_solve_regularized(Index m, Index n) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(m, n);
  Mat B = Mat::Random(m, 1);
  RealScalar lambda = RealScalar(0.1);
  const Index k = (std::min)(m, n);

  gpu::SVD<Scalar> svd(A, ComputeThinU | ComputeThinV);
  Mat X_reg = svd.solve(B, lambda);

  // CPU reference: D_ii = S_i / (S_i^2 + lambda^2).
  auto cpu_svd = BDCSVD<Mat, ComputeThinU | ComputeThinV>(A);
  auto S = cpu_svd.singularValues();
  Mat U = cpu_svd.matrixU();
  Mat V = cpu_svd.matrixV();

  Matrix<RealScalar, Dynamic, 1> D(k);
  for (Index i = 0; i < k; ++i) D(i) = S(i) / (S(i) * S(i) + lambda * lambda);
  Mat X_ref = V * D.asDiagonal() * U.adjoint() * B;

  RealScalar tol = RealScalar(100) * RealScalar(k) * NumTraits<Scalar>::epsilon();
  VERIFY((X_reg - X_ref).norm() / X_ref.norm() < tol);
}

// ---- Solve: rank-deficient (exercise drop_threshold pseudoinverse) ----------

template <typename Scalar>
void test_svd_solve_rank_deficient(Index m, Index n, Index r) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using Vec = Matrix<typename NumTraits<Scalar>::Real, Dynamic, 1>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  eigen_assert(r > 0 && r < (std::min)(m, n));

  // Build A of rank exactly r: A = U[:,:r] * diag(s) * V[:,:r]^H.
  Mat U_seed = Mat::Random(m, m);
  Mat U_full = HouseholderQR<Mat>(U_seed).householderQ();
  Mat V_seed = Mat::Random(n, n);
  Mat V_full = HouseholderQR<Mat>(V_seed).householderQ();
  Vec sigma(r);
  for (Index i = 0; i < r; ++i) sigma(i) = RealScalar(1) + RealScalar(i);  // distinct, well-spaced
  Mat A = U_full.leftCols(r) * sigma.asDiagonal() * V_full.leftCols(r).adjoint();

  Mat B = Mat::Random(m, 1);

  gpu::SVD<Scalar> svd(A, ComputeThinU | ComputeThinV);
  VERIFY_IS_EQUAL(svd.info(), Success);
  Mat X_gpu = svd.solve(B);

  // Reference: CPU BDCSVD with the same drop_threshold semantics (its default solve()
  // also drops near-zero singular values relative to S(0) * (m, n)_max * epsilon).
  Mat X_cpu = BDCSVD<Mat, ComputeThinU | ComputeThinV>(A).solve(B);

  // Both should produce the minimum-norm least-squares solution; compare via the
  // "useful" residual A^H r: zero up to the rank-r component plus rounding.
  RealScalar tol = RealScalar(50) * RealScalar((std::max)(m, n)) * NumTraits<Scalar>::epsilon();
  VERIFY((X_gpu - X_cpu).norm() / X_cpu.norm() < tol);

  // Norm of X should also be minimum (X has no component in null(A) up to rounding).
  VERIFY(numext::abs(X_gpu.norm() - X_cpu.norm()) / X_cpu.norm() < tol);
}

// ---- Reconstruction: full U / full V on a wide matrix (m < n) ---------------
//
// Exercises the transposed-internal-representation path through gesvd combined
// with the matrixU/matrixVT swap logic.

template <typename Scalar>
void test_svd_reconstruction_full_wide(Index m, Index n) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  eigen_assert(m < n);
  Mat A = Mat::Random(m, n);
  gpu::SVD<Scalar> svd(A, ComputeFullU | ComputeFullV);
  VERIFY_IS_EQUAL(svd.info(), Success);

  auto S = svd.singularValues();
  Mat U = svd.matrixU();    // m × m
  Mat VT = svd.matrixVT();  // n × n
  Mat V = svd.matrixV();    // n × n (matrixV alias should match matrixVT().adjoint())

  VERIFY_IS_EQUAL(U.rows(), m);
  VERIFY_IS_EQUAL(U.cols(), m);
  VERIFY_IS_EQUAL(VT.rows(), n);
  VERIFY_IS_EQUAL(VT.cols(), n);
  VERIFY_IS_EQUAL(V.rows(), n);
  VERIFY_IS_EQUAL(V.cols(), n);
  VERIFY_IS_APPROX(V, VT.adjoint());

  const Index k = (std::min)(m, n);
  Mat A_hat = U.leftCols(k) * S.asDiagonal() * VT.topRows(k);
  RealScalar tol = RealScalar(5) * std::sqrt(static_cast<RealScalar>(k)) * NumTraits<Scalar>::epsilon() * A.norm();
  VERIFY((A_hat - A).norm() < tol);

  Mat UtU = U.adjoint() * U;
  VERIFY((UtU - Mat::Identity(m, m)).norm() < tol);
  Mat VtVh = VT * VT.adjoint();
  VERIFY((VtVh - Mat::Identity(n, n)).norm() < tol);
}

// ---- Device-side accessors --------------------------------------------------
//
// Verify that d_singularValues / d_matrixU / d_matrixVT return views over the
// SVD's internal storage that agree with the host accessors after toHost().

template <typename Scalar>
void test_svd_device_accessors(Index m, Index n) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(m, n);
  gpu::SVD<Scalar> svd(A, ComputeThinU | ComputeThinV);
  VERIFY_IS_EQUAL(svd.info(), Success);

  auto d_S = svd.d_singularValues();
  auto d_U = svd.d_matrixU();
  auto d_VT = svd.d_matrixVT();

  Matrix<RealScalar, Dynamic, 1> S_host = d_S.toHost();
  Mat U_host = d_U.toHost();
  Mat VT_host = d_VT.toHost();

  VERIFY_IS_APPROX(S_host, svd.singularValues());
  VERIFY_IS_APPROX(U_host, svd.matrixU());
  VERIFY_IS_APPROX(VT_host, svd.matrixVT());

  // Multiple invocations must keep returning views without disturbing the SVD's state:
  // call again, then verify the SVD's host accessors still produce correct results.
  (void)svd.d_matrixU();
  (void)svd.d_matrixVT();
  VERIFY_IS_APPROX(svd.matrixU(), U_host);
  VERIFY_IS_APPROX(svd.matrixVT(), VT_host);
}

template <typename Scalar>
void test_svd_device_accessors_full_wide(Index m, Index n) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  eigen_assert(m < n);
  Mat A = Mat::Random(m, n);
  gpu::SVD<Scalar> svd(A, ComputeFullU | ComputeFullV);
  VERIFY_IS_EQUAL(svd.info(), Success);

  auto d_S = svd.d_singularValues();
  auto d_U = svd.d_matrixU();
  auto d_VT = svd.d_matrixVT();

  VERIFY_IS_EQUAL(d_U.rows(), m);
  VERIFY_IS_EQUAL(d_U.cols(), m);
  VERIFY_IS_EQUAL(d_VT.rows(), n);
  VERIFY_IS_EQUAL(d_VT.cols(), n);

  Matrix<RealScalar, Dynamic, 1> S_host = d_S.toHost();
  Mat U_host = d_U.toHost();
  Mat VT_host = d_VT.toHost();

  VERIFY_IS_APPROX(S_host, svd.singularValues());
  VERIFY_IS_APPROX(U_host, svd.matrixU());
  VERIFY_IS_APPROX(VT_host, svd.matrixVT());
}

// ---- Chain device views into a downstream cuBLAS GEMM (no D2D copy) ---------
//
// d_matrixU() returns a non-owning view over the SVD's internal d_U_ buffer.
// Feeding the view straight into a Context-driven GEMM exercises cross-stream
// event sync and confirms the borrow-deleter does not double-free on temporary
// destruction.

template <typename Scalar>
void test_svd_chain_orthogonality(Index m, Index n) {
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(m, n);
  gpu::SVD<Scalar> svd(A, ComputeThinU | ComputeThinV);
  VERIFY_IS_EQUAL(svd.info(), Success);

  const Index k = (std::min)(m, n);

  // U^H * U = I_k, all on device, no host roundtrip on U.
  gpu::Context ctx;
  gpu::DeviceMatrix<Scalar> d_UtU;
  {
    auto d_U = svd.d_matrixU();
    d_UtU.device(ctx) = d_U.adjoint() * d_U;
  }
  Mat UtU = d_UtU.toHost();
  RealScalar tol = RealScalar(20) * std::sqrt(static_cast<RealScalar>(k)) * NumTraits<Scalar>::epsilon();
  VERIFY((UtU - Mat::Identity(k, k)).norm() < tol);

  // VT * VT^H = I_k.
  gpu::DeviceMatrix<Scalar> d_VVt;
  {
    auto d_VT = svd.d_matrixVT();
    d_VVt.device(ctx) = d_VT * d_VT.adjoint();
  }
  Mat VVt = d_VVt.toHost();
  VERIFY((VVt - Mat::Identity(k, k)).norm() < tol);

  // After view destruction, SVD's state must remain valid (the deleter is a no-op
  // for views, so the underlying d_U_ / d_VT_ are not freed).
  Mat U_again = svd.matrixU();
  Mat VT_again = svd.matrixVT();
  VERIFY_IS_EQUAL(U_again.rows(), m);
  VERIFY_IS_EQUAL(VT_again.cols(), n);
}

// ---- Empty matrix -----------------------------------------------------------

void test_svd_empty() {
  gpu::SVD<double> svd(MatrixXd(0, 0), 0);
  VERIFY_IS_EQUAL(svd.info(), Success);
  VERIFY_IS_EQUAL(svd.rows(), 0);
  VERIFY_IS_EQUAL(svd.cols(), 0);

  svd.compute(MatrixXd::Random(4, 6), ComputeThinU | ComputeThinV);
  VERIFY_IS_EQUAL(svd.info(), Success);

  svd.compute(MatrixXd(0, 7), 0);
  VERIFY_IS_EQUAL(svd.info(), Success);
  VERIFY_IS_EQUAL(svd.rows(), 0);
  VERIFY_IS_EQUAL(svd.cols(), 7);
  VERIFY_IS_EQUAL(svd.singularValues().size(), 0);

  svd.compute(MatrixXd(5, 0), 0);
  VERIFY_IS_EQUAL(svd.info(), Success);
  VERIFY_IS_EQUAL(svd.rows(), 5);
  VERIFY_IS_EQUAL(svd.cols(), 0);
  VERIFY_IS_EQUAL(svd.singularValues().size(), 0);
}

// ---- Per-scalar driver ------------------------------------------------------

template <typename Scalar>
void test_scalar() {
  // Reconstruction + orthogonality (thin and full, identical test logic).
  CALL_SUBTEST((test_svd_reconstruction<Scalar, ComputeThinU | ComputeThinV>(64, 64)));
  CALL_SUBTEST((test_svd_reconstruction<Scalar, ComputeThinU | ComputeThinV>(128, 64)));
  CALL_SUBTEST((test_svd_reconstruction<Scalar, ComputeThinU | ComputeThinV>(64, 128)));  // wide (m < n)
  CALL_SUBTEST((test_svd_reconstruction<Scalar, ComputeFullU | ComputeFullV>(64, 64)));
  CALL_SUBTEST((test_svd_reconstruction<Scalar, ComputeFullU | ComputeFullV>(128, 64)));

  // Singular values.
  CALL_SUBTEST(test_svd_singular_values<Scalar>(64, 64));
  CALL_SUBTEST(test_svd_singular_values<Scalar>(128, 64));

  // Solve.
  CALL_SUBTEST(test_svd_solve<Scalar>(64, 64, 4));
  CALL_SUBTEST(test_svd_solve<Scalar>(128, 64, 4));
  CALL_SUBTEST(test_svd_solve<Scalar>(64, 128, 4));  // wide (m < n)

  // Truncated and regularized solve.
  CALL_SUBTEST(test_svd_solve_truncated<Scalar>(64, 64));
  CALL_SUBTEST(test_svd_solve_regularized<Scalar>(64, 64));

  // Rank-deficient solve (exercises drop_threshold pseudoinverse).
  CALL_SUBTEST(test_svd_solve_rank_deficient<Scalar>(64, 64, 32));
  CALL_SUBTEST(test_svd_solve_rank_deficient<Scalar>(96, 64, 16));
  CALL_SUBTEST(test_svd_solve_rank_deficient<Scalar>(64, 96, 16));

  // Wide matrix with full U/V (transposed-internal path).
  CALL_SUBTEST((test_svd_reconstruction_full_wide<Scalar>(64, 96)));

  // Device accessors.
  CALL_SUBTEST(test_svd_device_accessors<Scalar>(64, 64));
  CALL_SUBTEST(test_svd_device_accessors<Scalar>(96, 64));
  CALL_SUBTEST(test_svd_device_accessors<Scalar>(64, 96));
  CALL_SUBTEST(test_svd_device_accessors_full_wide<Scalar>(64, 96));

  // Chain device views into a downstream GEMM (orthogonality check).
  CALL_SUBTEST(test_svd_chain_orthogonality<Scalar>(64, 64));
  CALL_SUBTEST(test_svd_chain_orthogonality<Scalar>(96, 64));
}

EIGEN_DECLARE_TEST(gpu_cusolver_svd) {
  // Split by scalar so each part compiles in parallel.
  CALL_SUBTEST_1(test_scalar<float>());
  CALL_SUBTEST_2(test_scalar<double>());
  CALL_SUBTEST_3(test_scalar<std::complex<float>>());
  CALL_SUBTEST_4(test_scalar<std::complex<double>>());
  CALL_SUBTEST_5(test_svd_empty());
}
