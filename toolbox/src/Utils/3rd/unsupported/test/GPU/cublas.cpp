// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Rasmus Munk Larsen <rmlarsen@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

// Tests for cuBLAS GEMM dispatch via DeviceMatrix expression syntax.
// Covers: d_C = d_A * d_B, adjoint, transpose, scaled, +=, .device(ctx).

#define EIGEN_USE_GPU
#include "main.h"
#include <unsupported/Eigen/GPU>

using namespace Eigen;

// ---- Basic GEMM: C = A * B -------------------------------------------------

template <typename Scalar>
void test_gemm_basic(Index m, Index n, Index k) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(m, k);
  Mat B = Mat::Random(k, n);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  // Expression: d_C = d_A * d_B
  gpu::DeviceMatrix<Scalar> d_C;
  d_C = d_A * d_B;

  Mat C = d_C.toHost();
  Mat C_ref = A * B;

  RealScalar tol = RealScalar(k) * NumTraits<Scalar>::epsilon() * C_ref.norm();
  VERIFY((C - C_ref).norm() < tol);
}

// ---- GEMM with adjoint: C = A^H * B ----------------------------------------

template <typename Scalar>
void test_gemm_adjoint_lhs(Index m, Index n, Index k) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(k, m);  // A is k×m, A^H is m×k
  Mat B = Mat::Random(k, n);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::DeviceMatrix<Scalar> d_C;
  d_C = d_A.adjoint() * d_B;

  Mat C = d_C.toHost();
  Mat C_ref = A.adjoint() * B;

  RealScalar tol = RealScalar(k) * NumTraits<Scalar>::epsilon() * C_ref.norm();
  VERIFY((C - C_ref).norm() < tol);
}

// ---- GEMM with transpose: C = A * B^T --------------------------------------

template <typename Scalar>
void test_gemm_transpose_rhs(Index m, Index n, Index k) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(m, k);
  Mat B = Mat::Random(n, k);  // B is n×k, B^T is k×n

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::DeviceMatrix<Scalar> d_C;
  d_C = d_A * d_B.transpose();

  Mat C = d_C.toHost();
  Mat C_ref = A * B.transpose();

  RealScalar tol = RealScalar(k) * NumTraits<Scalar>::epsilon() * C_ref.norm();
  VERIFY((C - C_ref).norm() < tol);
}

// ---- GEMM with view/view operands: C = alpha * A^T * B^H -------------------

template <typename Scalar>
void test_gemm_view_pairs(Index m, Index n, Index k) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(k, m);  // A^T is m x k
  Mat B = Mat::Random(n, k);  // B^H is k x n
  Scalar alpha = Scalar(2);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::DeviceMatrix<Scalar> d_C;
  d_C = (alpha * d_A.transpose()) * d_B.adjoint();

  Mat C = d_C.toHost();
  Mat C_ref = (alpha * A.transpose()) * B.adjoint();

  RealScalar tol = RealScalar(k) * NumTraits<Scalar>::epsilon() * C_ref.norm();
  VERIFY((C - C_ref).norm() < tol);
}

// ---- GEMM with scaled: C = alpha * A * B ------------------------------------

template <typename Scalar>
void test_gemm_scaled(Index m, Index n, Index k) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(m, k);
  Mat B = Mat::Random(k, n);
  Scalar alpha = Scalar(2.5);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::DeviceMatrix<Scalar> d_C;
  d_C = alpha * d_A * d_B;

  Mat C = d_C.toHost();
  Mat C_ref = alpha * A * B;

  RealScalar tol = RealScalar(k) * NumTraits<Scalar>::epsilon() * C_ref.norm();
  VERIFY((C - C_ref).norm() < tol);
}

// ---- GEMM accumulate: C += A * B (beta=1) -----------------------------------

template <typename Scalar>
void test_gemm_accumulate(Index m, Index n, Index k) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(m, k);
  Mat B = Mat::Random(k, n);
  Mat C_init = Mat::Random(m, n);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);
  auto d_C = gpu::DeviceMatrix<Scalar>::fromHost(C_init);

  d_C += d_A * d_B;

  Mat C = d_C.toHost();
  Mat C_ref = C_init + A * B;

  RealScalar tol = RealScalar(k) * NumTraits<Scalar>::epsilon() * C_ref.norm();
  VERIFY((C - C_ref).norm() < tol);
}

// ---- GEMM accumulate into empty destination ---------------------------------

template <typename Scalar>
void test_gemm_accumulate_empty(Index m, Index n, Index k) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(m, k);
  Mat B = Mat::Random(k, n);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);
  gpu::DeviceMatrix<Scalar> d_C;

  d_C += d_A * d_B;

  Mat C = d_C.toHost();
  Mat C_ref = A * B;

  RealScalar tol = RealScalar(k) * NumTraits<Scalar>::epsilon() * C_ref.norm();
  VERIFY((C - C_ref).norm() < tol);
}

// ---- GEMM subtract: C -= A * B (beta=1, alpha=-1) --------------------------

template <typename Scalar>
void test_gemm_subtract(Index m, Index n, Index k) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(m, k);
  Mat B = Mat::Random(k, n);
  Mat C_init = Mat::Random(m, n);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);
  auto d_C = gpu::DeviceMatrix<Scalar>::fromHost(C_init);

  gpu::Context ctx;
  d_C.device(ctx) -= d_A * d_B;

  Mat C = d_C.toHost();
  Mat C_ref = C_init - A * B;

  RealScalar tol = RealScalar(k) * NumTraits<Scalar>::epsilon() * C_ref.norm();
  VERIFY((C - C_ref).norm() < tol);
}

// ---- GEMM subtract from empty destination -----------------------------------

template <typename Scalar>
void test_gemm_subtract_empty(Index m, Index n, Index k) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(m, k);
  Mat B = Mat::Random(k, n);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::Context ctx;
  gpu::DeviceMatrix<Scalar> d_C;
  d_C.device(ctx) -= d_A * d_B;

  Mat C = d_C.toHost();
  Mat C_ref = -(A * B);

  RealScalar tol = RealScalar(k) * NumTraits<Scalar>::epsilon() * C_ref.norm();
  VERIFY((C - C_ref).norm() < tol);
}

// ---- GEMM with scaled RHS: C = A * (alpha * B) -----------------------------

template <typename Scalar>
void test_gemm_scaled_rhs(Index m, Index n, Index k) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(m, k);
  Mat B = Mat::Random(k, n);
  Scalar alpha = Scalar(3.0);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::DeviceMatrix<Scalar> d_C;
  d_C = d_A * (alpha * d_B);

  Mat C = d_C.toHost();
  Mat C_ref = A * (alpha * B);

  RealScalar tol = RealScalar(k) * NumTraits<Scalar>::epsilon() * C_ref.norm();
  VERIFY((C - C_ref).norm() < tol);
}

// ---- GEMM dimension mismatch must assert ------------------------------------
// Note: we do NOT use VERIFY_RAISES_ASSERT here because it relies on
// setjmp/longjmp which skips RAII destructors for DeviceMatrix (GPU memory)
// and cuBLAS/cuSOLVER handles, corrupting state for subsequent tests.

// ---- GEMM with explicit gpu::Context ------------------------------------------

template <typename Scalar>
void test_gemm_explicit_context(Index m, Index n, Index k) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(m, k);
  Mat B = Mat::Random(k, n);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::Context ctx;
  gpu::DeviceMatrix<Scalar> d_C;
  d_C.device(ctx) = d_A * d_B;

  Mat C = d_C.toHost();
  Mat C_ref = A * B;

  RealScalar tol = RealScalar(k) * NumTraits<Scalar>::epsilon() * C_ref.norm();
  VERIFY((C - C_ref).norm() < tol);
}

// ---- GEMM cross-context reuse of the same destination -----------------------

template <typename Scalar>
void test_gemm_cross_context_reuse(Index n) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(n, n);
  Mat B = Mat::Random(n, n);
  Mat D = Mat::Random(n, n);
  Mat E = Mat::Random(n, n);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);
  auto d_D = gpu::DeviceMatrix<Scalar>::fromHost(D);
  auto d_E = gpu::DeviceMatrix<Scalar>::fromHost(E);

  gpu::Context ctx1;
  gpu::Context ctx2;
  gpu::DeviceMatrix<Scalar> d_C;
  d_C.device(ctx1) = d_A * d_B;
  d_C.device(ctx2) += d_D * d_E;

  Mat C = d_C.toHost();
  Mat C_ref = A * B + D * E;

  RealScalar tol = RealScalar(2) * RealScalar(n) * NumTraits<Scalar>::epsilon() * C_ref.norm();
  VERIFY((C - C_ref).norm() < tol);
}

// ---- GEMM cross-context resize of the destination ---------------------------

template <typename Scalar>
void test_gemm_cross_context_resize() {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(64, 64);
  Mat B = Mat::Random(64, 64);
  Mat D = Mat::Random(32, 16);
  Mat E = Mat::Random(16, 8);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);
  auto d_D = gpu::DeviceMatrix<Scalar>::fromHost(D);
  auto d_E = gpu::DeviceMatrix<Scalar>::fromHost(E);

  gpu::Context ctx1;
  gpu::Context ctx2;
  gpu::DeviceMatrix<Scalar> d_C;
  d_C.device(ctx1) = d_A * d_B;
  d_C.device(ctx2) = d_D * d_E;

  Mat C = d_C.toHost();
  Mat C_ref = D * E;

  RealScalar tol = RealScalar(16) * NumTraits<Scalar>::epsilon() * C_ref.norm();
  VERIFY((C - C_ref).norm() < tol);
}

// ---- GEMM chaining: C = (A * B) then D = C * E -----------------------------

template <typename Scalar>
void test_gemm_chain(Index n) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(n, n);
  Mat B = Mat::Random(n, n);
  Mat E = Mat::Random(n, n);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);
  auto d_E = gpu::DeviceMatrix<Scalar>::fromHost(E);

  gpu::DeviceMatrix<Scalar> d_C;
  d_C = d_A * d_B;
  gpu::DeviceMatrix<Scalar> d_D;
  d_D = d_C * d_E;

  Mat D = d_D.toHost();
  Mat D_ref = (A * B) * E;

  RealScalar tol = RealScalar(2) * RealScalar(n) * NumTraits<Scalar>::epsilon() * D_ref.norm();
  VERIFY((D - D_ref).norm() < tol);
}

// ---- Square identity check: A * I = A ---------------------------------------

template <typename Scalar>
void test_gemm_identity(Index n) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;

  Mat A = Mat::Random(n, n);
  Mat eye = Mat::Identity(n, n);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_I = gpu::DeviceMatrix<Scalar>::fromHost(eye);

  gpu::DeviceMatrix<Scalar> d_C;
  d_C = d_A * d_I;

  Mat C = d_C.toHost();
  VERIFY_IS_APPROX(C, A);
}

// ---- LLT solve expression: d_X = d_A.llt().solve(d_B) ----------------------

template <typename MatrixType>
MatrixType make_spd(Index n) {
  using Scalar = typename MatrixType::Scalar;
  MatrixType M = MatrixType::Random(n, n);
  return M.adjoint() * M + MatrixType::Identity(n, n) * static_cast<Scalar>(n);
}

template <typename Scalar>
void test_llt_solve_expr(Index n, Index nrhs) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = make_spd<Mat>(n);
  Mat B = Mat::Random(n, nrhs);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::DeviceMatrix<Scalar> d_X;
  d_X = d_A.llt().solve(d_B);

  Mat X = d_X.toHost();
  RealScalar residual = (A * X - B).norm() / B.norm();
  VERIFY(residual < RealScalar(n) * NumTraits<Scalar>::epsilon());
}

// ---- LLT solve with explicit context ----------------------------------------

template <typename Scalar>
void test_llt_solve_expr_context(Index n, Index nrhs) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = make_spd<Mat>(n);
  Mat B = Mat::Random(n, nrhs);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::Context ctx;
  gpu::DeviceMatrix<Scalar> d_X;
  d_X.device(ctx) = d_A.llt().solve(d_B);

  Mat X = d_X.toHost();
  RealScalar residual = (A * X - B).norm() / B.norm();
  VERIFY(residual < RealScalar(n) * NumTraits<Scalar>::epsilon());
}

// ---- LU solve expression: d_X = d_A.lu().solve(d_B) ------------------------

template <typename Scalar>
void test_lu_solve_expr(Index n, Index nrhs) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(n, n);
  Mat B = Mat::Random(n, nrhs);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::DeviceMatrix<Scalar> d_X;
  d_X = d_A.lu().solve(d_B);

  Mat X = d_X.toHost();
  RealScalar residual = (A * X - B).norm() / (A.norm() * X.norm());
  VERIFY(residual < RealScalar(10) * RealScalar(n) * NumTraits<Scalar>::epsilon());
}

// ---- GEMM + solver chain: C = A * B, X = C.llt().solve(D) ------------------

template <typename Scalar>
void test_gemm_then_solve(Index n) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(n, n);
  Mat D = Mat::Random(n, 1);

  // Make SPD: C = A^H * A + n*I
  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  gpu::DeviceMatrix<Scalar> d_C;
  d_C = d_A.adjoint() * d_A;

  // Add n*I on host (no element-wise ops on DeviceMatrix yet).
  Mat C = d_C.toHost();
  C += Mat::Identity(n, n) * static_cast<Scalar>(n);
  d_C = gpu::DeviceMatrix<Scalar>::fromHost(C);

  auto d_D = gpu::DeviceMatrix<Scalar>::fromHost(D);

  gpu::DeviceMatrix<Scalar> d_X;
  d_X = d_C.llt().solve(d_D);

  Mat X = d_X.toHost();
  RealScalar residual = (C * X - D).norm() / D.norm();
  VERIFY(residual < RealScalar(n) * NumTraits<Scalar>::epsilon());
}

// ---- LLT solve with Upper triangle -----------------------------------------

template <typename Scalar>
void test_llt_solve_upper(Index n, Index nrhs) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = make_spd<Mat>(n);
  Mat B = Mat::Random(n, nrhs);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::DeviceMatrix<Scalar> d_X;
  d_X = d_A.template llt<Upper>().solve(d_B);

  Mat X = d_X.toHost();
  RealScalar residual = (A * X - B).norm() / B.norm();
  VERIFY(residual < RealScalar(n) * NumTraits<Scalar>::epsilon());
}

// ---- LU solve with explicit context -----------------------------------------

template <typename Scalar>
void test_lu_solve_expr_context(Index n, Index nrhs) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(n, n);
  Mat B = Mat::Random(n, nrhs);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::Context ctx;
  gpu::DeviceMatrix<Scalar> d_X;
  d_X.device(ctx) = d_A.lu().solve(d_B);

  Mat X = d_X.toHost();
  RealScalar residual = (A * X - B).norm() / (A.norm() * X.norm());
  VERIFY(residual < RealScalar(10) * RealScalar(n) * NumTraits<Scalar>::epsilon());
}

// ---- Zero-nrhs solver expressions ------------------------------------------

template <typename Scalar>
void test_llt_solve_zero_nrhs(Index n) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;

  Mat A = make_spd<Mat>(n);
  Mat B = Mat::Random(n, 0);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::DeviceMatrix<Scalar> d_X;
  d_X = d_A.llt().solve(d_B);

  VERIFY_IS_EQUAL(d_X.rows(), n);
  VERIFY_IS_EQUAL(d_X.cols(), 0);
}

template <typename Scalar>
void test_lu_solve_zero_nrhs(Index n) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;

  Mat A = Mat::Random(n, n);
  Mat B = Mat::Random(n, 0);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::DeviceMatrix<Scalar> d_X;
  d_X = d_A.lu().solve(d_B);

  VERIFY_IS_EQUAL(d_X.rows(), n);
  VERIFY_IS_EQUAL(d_X.cols(), 0);
}

// ---- TRSM: triangularView<UpLo>().solve(B) ----------------------------------

template <typename Scalar, int UpLo>
void test_trsm(Index n, Index nrhs) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  // Build a well-conditioned triangular matrix.
  Mat A = Mat::Random(n, n);
  A.diagonal().array() += static_cast<Scalar>(n);  // ensure non-singular
  if (UpLo == Lower)
    A = A.template triangularView<Lower>();
  else
    A = A.template triangularView<Upper>();

  Mat B = Mat::Random(n, nrhs);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::DeviceMatrix<Scalar> d_X;
  d_X = d_A.template triangularView<UpLo>().solve(d_B);

  Mat X = d_X.toHost();
  RealScalar residual = (A * X - B).norm() / B.norm();
  VERIFY(residual < RealScalar(n) * NumTraits<Scalar>::epsilon());
}

// ---- SYMM/HEMM: selfadjointView<UpLo>() * B --------------------------------

template <typename Scalar, int UpLo>
void test_symm(Index n, Index nrhs) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = make_spd<Mat>(n);  // SPD is also self-adjoint
  Mat B = Mat::Random(n, nrhs);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);
  auto d_B = gpu::DeviceMatrix<Scalar>::fromHost(B);

  gpu::DeviceMatrix<Scalar> d_C;
  d_C = d_A.template selfadjointView<UpLo>() * d_B;

  Mat C = d_C.toHost();
  Mat C_ref = A * B;  // A is symmetric, so full multiply == symm

  RealScalar tol = RealScalar(n) * NumTraits<Scalar>::epsilon() * C_ref.norm();
  VERIFY((C - C_ref).norm() < tol);
}

// ---- SYRK/HERK: rankUpdate(A) → C = A * A^H --------------------------------

template <typename Scalar>
void test_syrk(Index n, Index k) {
  using Mat = Eigen::Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  Mat A = Mat::Random(n, k);

  auto d_A = gpu::DeviceMatrix<Scalar>::fromHost(A);

  gpu::DeviceMatrix<Scalar> d_C;
  d_C.template selfadjointView<Lower>().rankUpdate(d_A);

  Mat C = d_C.toHost();
  // Only lower triangle is meaningful for SYRK. Compare lower triangle.
  Mat C_ref = A * A.adjoint();

  // Extract lower triangle for comparison.
  Mat C_lower = C.template triangularView<Lower>();
  Mat C_ref_lower = C_ref.template triangularView<Lower>();

  RealScalar tol = RealScalar(k) * NumTraits<Scalar>::epsilon() * C_ref.norm();
  VERIFY((C_lower - C_ref_lower).norm() < tol);
}

// ---- Per-scalar driver ------------------------------------------------------

template <typename Scalar>
void test_scalar() {
  CALL_SUBTEST(test_gemm_basic<Scalar>(64, 64, 64));
  CALL_SUBTEST(test_gemm_basic<Scalar>(128, 64, 32));
  CALL_SUBTEST(test_gemm_basic<Scalar>(1, 1, 1));
  CALL_SUBTEST(test_gemm_basic<Scalar>(256, 256, 256));

  CALL_SUBTEST(test_gemm_adjoint_lhs<Scalar>(64, 64, 64));
  CALL_SUBTEST(test_gemm_adjoint_lhs<Scalar>(128, 32, 64));

  CALL_SUBTEST(test_gemm_transpose_rhs<Scalar>(64, 64, 64));
  CALL_SUBTEST(test_gemm_transpose_rhs<Scalar>(128, 32, 64));
  CALL_SUBTEST(test_gemm_view_pairs<Scalar>(64, 32, 48));

  CALL_SUBTEST(test_gemm_scaled<Scalar>(64, 64, 64));
  CALL_SUBTEST(test_gemm_scaled_rhs<Scalar>(64, 64, 64));
  CALL_SUBTEST(test_gemm_accumulate<Scalar>(64, 64, 64));
  CALL_SUBTEST(test_gemm_accumulate_empty<Scalar>(64, 64, 64));
  CALL_SUBTEST(test_gemm_subtract<Scalar>(64, 64, 64));
  CALL_SUBTEST(test_gemm_subtract_empty<Scalar>(64, 64, 64));
  CALL_SUBTEST(test_gemm_explicit_context<Scalar>(64, 64, 64));
  CALL_SUBTEST(test_gemm_cross_context_reuse<Scalar>(64));
  CALL_SUBTEST(test_gemm_cross_context_resize<Scalar>());
  CALL_SUBTEST(test_gemm_chain<Scalar>(64));
  CALL_SUBTEST(test_gemm_identity<Scalar>(64));

  // Solver expressions — zero-size edge cases (use dedicated tests, not residual-based)

  // Solver expressions
  CALL_SUBTEST(test_llt_solve_expr<Scalar>(64, 1));
  CALL_SUBTEST(test_llt_solve_expr<Scalar>(64, 4));
  CALL_SUBTEST(test_llt_solve_expr<Scalar>(256, 8));
  CALL_SUBTEST(test_llt_solve_expr_context<Scalar>(64, 4));
  CALL_SUBTEST(test_llt_solve_upper<Scalar>(64, 4));
  CALL_SUBTEST(test_lu_solve_expr<Scalar>(64, 1));
  CALL_SUBTEST(test_lu_solve_expr<Scalar>(64, 4));
  CALL_SUBTEST(test_lu_solve_expr<Scalar>(256, 8));
  CALL_SUBTEST(test_lu_solve_expr_context<Scalar>(64, 4));
  CALL_SUBTEST(test_llt_solve_zero_nrhs<Scalar>(64));
  CALL_SUBTEST(test_llt_solve_zero_nrhs<Scalar>(0));
  CALL_SUBTEST(test_lu_solve_zero_nrhs<Scalar>(64));
  CALL_SUBTEST(test_lu_solve_zero_nrhs<Scalar>(0));
  CALL_SUBTEST(test_gemm_then_solve<Scalar>(64));

  // TRSM
  CALL_SUBTEST((test_trsm<Scalar, Lower>(64, 1)));
  CALL_SUBTEST((test_trsm<Scalar, Lower>(64, 4)));
  CALL_SUBTEST((test_trsm<Scalar, Upper>(64, 4)));
  CALL_SUBTEST((test_trsm<Scalar, Lower>(256, 8)));

  // SYMM/HEMM
  CALL_SUBTEST((test_symm<Scalar, Lower>(64, 4)));
  CALL_SUBTEST((test_symm<Scalar, Upper>(64, 4)));
  CALL_SUBTEST((test_symm<Scalar, Lower>(128, 8)));

  // SYRK/HERK
  CALL_SUBTEST(test_syrk<Scalar>(64, 64));
  CALL_SUBTEST(test_syrk<Scalar>(64, 32));
  CALL_SUBTEST(test_syrk<Scalar>(128, 64));
}

// ---- Solver failure mode tests (not templated on Scalar) --------------------
// Use the cached GpuLLT/GpuLU API which reports failure via info() rather than
// the expression API which asserts inside dispatch (incompatible with
// VERIFY_RAISES_ASSERT due to longjmp skipping RAII destructors).

void test_llt_not_spd() {
  // Negative definite matrix — LLT factorization must fail.
  MatrixXd A = -MatrixXd::Identity(8, 8);
  gpu::LLT<double> llt(A);
  VERIFY_IS_EQUAL(llt.info(), NumericalIssue);
}

void test_lu_singular() {
  // Zero matrix — LU factorization must detect singularity.
  MatrixXd A = MatrixXd::Zero(8, 8);
  gpu::LU<double> lu(A);
  VERIFY_IS_EQUAL(lu.info(), NumericalIssue);
}

EIGEN_DECLARE_TEST(gpu_cublas) {
  // Split by scalar so each part compiles in parallel.
  CALL_SUBTEST_1(test_scalar<float>());
  CALL_SUBTEST_2(test_scalar<double>());
  CALL_SUBTEST_3(test_scalar<std::complex<float>>());
  CALL_SUBTEST_4(test_scalar<std::complex<double>>());
  CALL_SUBTEST_5(test_llt_not_spd());
  CALL_SUBTEST_5(test_lu_singular());
}
