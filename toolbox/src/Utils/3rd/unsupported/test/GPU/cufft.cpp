// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Rasmus Munk Larsen <rmlarsen@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

// Tests for GpuFFT: GPU FFT via cuFFT.

#define EIGEN_USE_GPU
#include "main.h"
#include <unsupported/Eigen/GPU>
#include "gpu_test_helpers.h"

using namespace Eigen;

// ---- 1D C2C roundtrip: inv(fwd(x)) ≈ x -------------------------------------

template <typename Scalar>
void test_c2c_roundtrip(Index n) {
  using Complex = std::complex<Scalar>;
  using Vec = Matrix<Complex, Dynamic, 1>;
  using RealScalar = Scalar;

  Vec x = Vec::Random(n);

  gpu::FFT<Scalar> fft;
  Vec X = fft.fwd(x);
  VERIFY_IS_EQUAL(X.size(), n);

  Vec y = fft.inv(X);
  VERIFY_IS_EQUAL(y.size(), n);

  RealScalar tol = RealScalar(10) * RealScalar(n) * NumTraits<Scalar>::epsilon();
  VERIFY((y - x).norm() / x.norm() < tol);
}

// ---- 1D C2C known signal: FFT of constant = delta --------------------------

template <typename Scalar>
void test_c2c_constant() {
  using Complex = std::complex<Scalar>;
  using Vec = Matrix<Complex, Dynamic, 1>;
  using RealScalar = Scalar;

  const int n = 64;
  Vec x = Vec::Constant(n, Complex(3.0, 0.0));

  gpu::FFT<Scalar> fft;
  Vec X = fft.fwd(x);

  // FFT of constant c: X[0] = c*n, X[k] = 0 for k > 0.
  RealScalar tol = RealScalar(10) * NumTraits<Scalar>::epsilon() * RealScalar(n);
  VERIFY(std::abs(X(0) - Complex(3.0 * n, 0.0)) < tol);
  for (int k = 1; k < n; ++k) {
    VERIFY(std::abs(X(k)) < tol);
  }
}

// ---- 1D R2C/C2R roundtrip: invReal(fwd(r), n) ≈ r --------------------------

template <typename Scalar>
void test_r2c_roundtrip(Index n) {
  using Complex = std::complex<Scalar>;
  using CVec = Matrix<Complex, Dynamic, 1>;
  using RVec = Matrix<Scalar, Dynamic, 1>;
  using RealScalar = Scalar;

  RVec r = RVec::Random(n);

  gpu::FFT<Scalar> fft;
  CVec R = fft.fwd(r);

  // R2C returns n/2+1 complex values.
  VERIFY_IS_EQUAL(R.size(), n / 2 + 1);

  RVec s = fft.invReal(R, n);
  VERIFY_IS_EQUAL(s.size(), n);

  RealScalar tol = RealScalar(10) * RealScalar(n) * NumTraits<Scalar>::epsilon();
  VERIFY((s - r).norm() / r.norm() < tol);
}

// ---- 2D C2C roundtrip: inv2(fwd2(A)) ≈ A ---------------------------------

template <typename Scalar>
void test_2d_roundtrip(Index rows, Index cols) {
  using Complex = std::complex<Scalar>;
  using Mat = Matrix<Complex, Dynamic, Dynamic>;
  using RealScalar = Scalar;

  Mat A = Mat::Random(rows, cols);

  gpu::FFT<Scalar> fft;
  Mat B = fft.fwd2(A);
  VERIFY_IS_EQUAL(B.rows(), rows);
  VERIFY_IS_EQUAL(B.cols(), cols);

  Mat C = fft.inv2(B);
  VERIFY_IS_EQUAL(C.rows(), rows);
  VERIFY_IS_EQUAL(C.cols(), cols);

  RealScalar tol = RealScalar(10) * RealScalar(rows * cols) * NumTraits<Scalar>::epsilon();
  VERIFY((C - A).norm() / A.norm() < tol);
}

// ---- 2D C2C known signal: constant matrix -----------------------------------

template <typename Scalar>
void test_2d_constant() {
  using Complex = std::complex<Scalar>;
  using Mat = Matrix<Complex, Dynamic, Dynamic>;
  using RealScalar = Scalar;

  const int rows = 16, cols = 32;
  Mat A = Mat::Constant(rows, cols, Complex(2.0, 0.0));

  gpu::FFT<Scalar> fft;
  Mat B = fft.fwd2(A);

  // 2D FFT of constant c: B(0,0) = c*rows*cols, all others = 0.
  RealScalar tol = RealScalar(10) * NumTraits<Scalar>::epsilon() * RealScalar(rows * cols);
  VERIFY(std::abs(B(0, 0) - Complex(2.0 * rows * cols, 0.0)) < tol);
  for (int j = 0; j < cols; ++j) {
    for (int i = 0; i < rows; ++i) {
      if (i == 0 && j == 0) continue;
      VERIFY(std::abs(B(i, j)) < tol);
    }
  }
}

// ---- Plan reuse: repeated calls should work ---------------------------------

template <typename Scalar>
void test_plan_reuse() {
  using Complex = std::complex<Scalar>;
  using Vec = Matrix<Complex, Dynamic, 1>;
  using RealScalar = Scalar;

  gpu::FFT<Scalar> fft;
  for (int trial = 0; trial < 5; ++trial) {
    Vec x = Vec::Random(128);
    Vec X = fft.fwd(x);
    Vec y = fft.inv(X);
    RealScalar tol = RealScalar(10) * RealScalar(128) * NumTraits<Scalar>::epsilon();
    VERIFY((y - x).norm() / x.norm() < tol);
  }
}

// ---- LRU eviction + user-configurable capacity ------------------------------
// Constructs an FFT with an explicit small plan-cache capacity (also exercising
// the user-configurable constructor), pushes 2 * capacity distinct shapes
// through it so the cache is forced to evict and recreate plans, verifies the
// cache stays clamped at capacity, then revisits the first (now-evicted) shape.
// Two correctness conditions: every transform stays numerically correct, and
// cufftDestroy on eviction leaves no dangling state on the stream.

template <typename Scalar>
void test_plan_cache_eviction() {
  using Complex = std::complex<Scalar>;
  using Vec = Matrix<Complex, Dynamic, 1>;
  using RealScalar = Scalar;

  constexpr std::size_t capacity = 4;
  gpu::FFT<Scalar> fft(capacity);
  VERIFY_IS_EQUAL(fft.plan_cache_capacity(), capacity);
  VERIFY_IS_EQUAL(fft.plan_cache_size(), std::size_t(0));

  const int num_shapes = static_cast<int>(2 * capacity);
  for (int k = 0; k < num_shapes; ++k) {
    const Index n = 32 + 16 * k;  // 32, 48, 64, ... all distinct.
    Vec x = Vec::Random(n);
    Vec X = fft.fwd(x);
    Vec y = fft.inv(X);
    RealScalar tol = RealScalar(10) * RealScalar(n) * NumTraits<Scalar>::epsilon();
    VERIFY((y - x).norm() / x.norm() < tol);
  }
  // Cache should be clamped at capacity, never exceeding it.
  VERIFY_IS_EQUAL(fft.plan_cache_size(), capacity);

  // Smallest shape has been evicted by now; re-creating it must succeed.
  Vec x = Vec::Random(32);
  Vec X = fft.fwd(x);
  Vec y = fft.inv(X);
  RealScalar tol = RealScalar(10) * RealScalar(32) * NumTraits<Scalar>::epsilon();
  VERIFY((y - x).norm() / x.norm() < tol);
  VERIFY_IS_EQUAL(fft.plan_cache_size(), capacity);
}

// ---- Empty ------------------------------------------------------------------

template <typename Scalar>
void test_empty() {
  using Complex = std::complex<Scalar>;
  using Vec = Matrix<Complex, Dynamic, 1>;

  gpu::FFT<Scalar> fft;
  Vec x(0);
  Vec X = fft.fwd(x);
  VERIFY_IS_EQUAL(X.size(), 0);
  Vec y = fft.inv(X);
  VERIFY_IS_EQUAL(y.size(), 0);
}

// ---- Explicit-Context ctor: stream and cuBLAS plumbing ----------------------
// Verifies that FFT(Context&) actually binds to the supplied Context's stream
// and cuBLAS handle (rather than falling back to Context::threadLocal()), and
// that fwd/inv exercise both: fwd schedules the cuFFT plan on ctx.stream(),
// inv additionally scales via ctx.cublasHandle().

template <typename Scalar>
void test_explicit_context(Index n) {
  using Complex = std::complex<Scalar>;
  using Vec = Matrix<Complex, Dynamic, 1>;
  using RealScalar = Scalar;

  // A fresh Context owns a distinct stream from the thread-local default.
  gpu::Context ctx;
  VERIFY(ctx.stream() != gpu::Context::threadLocal().stream());

  gpu::FFT<Scalar> fft(ctx);
  VERIFY(fft.stream() == ctx.stream());
  VERIFY(&fft.context() == &ctx);

  Vec x = Vec::Random(n);
  Vec X = fft.fwd(x);
  Vec y = fft.inv(X);
  RealScalar tol = RealScalar(10) * RealScalar(n) * NumTraits<Scalar>::epsilon();
  VERIFY((y - x).norm() / x.norm() < tol);
}

// ---- Per-scalar driver ------------------------------------------------------

template <typename Scalar>
void test_scalar() {
  CALL_SUBTEST(test_c2c_roundtrip<Scalar>(64));
  CALL_SUBTEST(test_c2c_roundtrip<Scalar>(256));
  CALL_SUBTEST(test_c2c_roundtrip<Scalar>(1000));  // non-power-of-2
  CALL_SUBTEST(test_c2c_constant<Scalar>());
  CALL_SUBTEST(test_r2c_roundtrip<Scalar>(64));
  CALL_SUBTEST(test_r2c_roundtrip<Scalar>(256));
  CALL_SUBTEST(test_2d_roundtrip<Scalar>(32, 32));
  CALL_SUBTEST(test_2d_roundtrip<Scalar>(16, 64));  // non-square
  CALL_SUBTEST(test_2d_constant<Scalar>());
  CALL_SUBTEST(test_plan_reuse<Scalar>());
  CALL_SUBTEST(test_plan_cache_eviction<Scalar>());
  CALL_SUBTEST(test_empty<Scalar>());
  CALL_SUBTEST(test_explicit_context<Scalar>(64));
}

EIGEN_DECLARE_TEST(gpu_cufft) {
  gpu_test::require_cufft_context();
  // Split by scalar so each part compiles in parallel.
  CALL_SUBTEST_1(test_scalar<float>());
  CALL_SUBTEST_2(test_scalar<double>());
}
