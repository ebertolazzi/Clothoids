// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Rasmus Munk Larsen <rmlarsen@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

// Tests for GpuSparseLDLT: GPU sparse LDL^T via cuDSS.

#define EIGEN_USE_GPU
#include "main.h"
#include <Eigen/Sparse>
#include <unsupported/Eigen/GPU>
#include "gpu_test_helpers.h"

using namespace Eigen;

// ---- Helper: build a random sparse symmetric indefinite matrix ---------------

template <typename Scalar>
SparseMatrix<Scalar, ColMajor, int> make_symmetric_indefinite(Index n, double density = 0.1) {
  using SpMat = SparseMatrix<Scalar, ColMajor, int>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  // Diagonal has mixed signs for indefiniteness. For complex Scalar:
  // off-diagonals carry nonzero imaginary parts so the Hermitian path is
  // genuinely exercised; diagonal stays real so A = R + R.adjoint() has a
  // real diagonal explicitly (R + R.adjoint() at i,i = 2*Re(R_ii) regardless).
  SpMat R(n, n);
  R.reserve(VectorXi::Constant(n, static_cast<int>(n * density) + 1));
  for (Index j = 0; j < n; ++j) {
    for (Index i = 0; i < n; ++i) {
      if (i == j || (std::rand() / double(RAND_MAX)) < density) {
        const RealScalar re = RealScalar(std::rand() / double(RAND_MAX) - 0.5);
        const RealScalar im = (i == j) ? RealScalar(0) : RealScalar(std::rand() / double(RAND_MAX) - 0.5);
        R.insert(i, j) = gpu_test::make_test_value<Scalar>(re, im);
      }
    }
  }
  R.makeCompressed();

  // A = R + R^H (symmetric), then add diagonal with alternating signs for indefiniteness.
  SpMat A = R + SparseMatrix<Scalar, ColMajor, int>(R.adjoint());
  for (Index i = 0; i < n; ++i) {
    Scalar diag_val = Scalar((i % 2 == 0) ? n : -n);
    A.coeffRef(i, i) += diag_val;
  }
  A.makeCompressed();
  return A;
}

// ---- Solve and check residual -----------------------------------------------

template <typename Scalar>
void test_solve(Index n) {
  using SpMat = SparseMatrix<Scalar, ColMajor, int>;
  using Vec = Matrix<Scalar, Dynamic, 1>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  SpMat A = make_symmetric_indefinite<Scalar>(n);
  Vec b = Vec::Random(n);

  gpu::SparseLDLT<Scalar> ldlt(A);
  VERIFY_IS_EQUAL(ldlt.info(), Success);

  Vec x = ldlt.solve(b);
  VERIFY_IS_EQUAL(x.rows(), n);

  Vec r = A * x - b;
  RealScalar tol = RealScalar(100) * RealScalar(n) * NumTraits<Scalar>::epsilon();
  VERIFY(r.norm() / b.norm() < tol);
}

// ---- Multiple RHS -----------------------------------------------------------

template <typename Scalar>
void test_multiple_rhs(Index n, Index nrhs) {
  using SpMat = SparseMatrix<Scalar, ColMajor, int>;
  using Mat = Matrix<Scalar, Dynamic, Dynamic>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  SpMat A = make_symmetric_indefinite<Scalar>(n);
  Mat B = Mat::Random(n, nrhs);

  gpu::SparseLDLT<Scalar> ldlt(A);
  VERIFY_IS_EQUAL(ldlt.info(), Success);

  Mat X = ldlt.solve(B);
  VERIFY_IS_EQUAL(X.rows(), n);
  VERIFY_IS_EQUAL(X.cols(), nrhs);

  Mat R = A * X - B;
  RealScalar tol = RealScalar(100) * RealScalar(n) * NumTraits<Scalar>::epsilon();
  VERIFY(R.norm() / B.norm() < tol);
}

// ---- Refactorize ------------------------------------------------------------

template <typename Scalar>
void test_refactorize(Index n) {
  using SpMat = SparseMatrix<Scalar, ColMajor, int>;
  using Vec = Matrix<Scalar, Dynamic, 1>;
  using RealScalar = typename NumTraits<Scalar>::Real;

  SpMat A = make_symmetric_indefinite<Scalar>(n);
  Vec b = Vec::Random(n);

  gpu::SparseLDLT<Scalar> ldlt;
  ldlt.analyzePattern(A);
  VERIFY_IS_EQUAL(ldlt.info(), Success);

  ldlt.factorize(A);
  VERIFY_IS_EQUAL(ldlt.info(), Success);
  Vec x1 = ldlt.solve(b);

  // Modify values, keep pattern.
  SpMat A2 = A;
  for (Index i = 0; i < n; ++i) A2.coeffRef(i, i) *= Scalar(RealScalar(2));

  ldlt.factorize(A2);
  VERIFY_IS_EQUAL(ldlt.info(), Success);
  Vec x2 = ldlt.solve(b);

  RealScalar tol = RealScalar(100) * RealScalar(n) * NumTraits<Scalar>::epsilon();
  VERIFY((A * x1 - b).norm() / b.norm() < tol);
  VERIFY((A2 * x2 - b).norm() / b.norm() < tol);
  // Diagonal scaled 2x; x1 and x2 must differ by a substantial fraction.
  VERIFY((x1 - x2).norm() > RealScalar(0.01) * x1.norm());
}

// ---- Empty ------------------------------------------------------------------

template <typename Scalar>
void test_empty() {
  using SpMat = SparseMatrix<Scalar, ColMajor, int>;
  SpMat A(0, 0);
  A.makeCompressed();
  gpu::SparseLDLT<Scalar> ldlt(A);
  VERIFY_IS_EQUAL(ldlt.info(), Success);
  VERIFY_IS_EQUAL(ldlt.rows(), 0);
  VERIFY_IS_EQUAL(ldlt.cols(), 0);
}

// ---- Per-scalar driver ------------------------------------------------------

template <typename Scalar>
void test_scalar() {
  CALL_SUBTEST(test_solve<Scalar>(64));
  CALL_SUBTEST(test_solve<Scalar>(256));
  CALL_SUBTEST(test_multiple_rhs<Scalar>(64, 4));
  CALL_SUBTEST(test_refactorize<Scalar>(64));
}

EIGEN_DECLARE_TEST(gpu_cudss_ldlt) {
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
