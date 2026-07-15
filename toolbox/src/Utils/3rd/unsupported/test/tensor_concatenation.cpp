// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2014 Benoit Steiner <benoit.steiner.goog@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

#include "main.h"

#include <complex>

#include <Eigen/Tensor>

using Eigen::Tensor;

template <int DataLayout>
static void test_dimension_failures() {
  Tensor<int, 3, DataLayout> left(2, 3, 1);
  Tensor<int, 3, DataLayout> right(3, 3, 1);
  left.setRandom();
  right.setRandom();

  // Okay; other dimensions are equal.
  Tensor<int, 3, DataLayout> concatenation = left.concatenate(right, 0);

  // Dimension mismatches.
  VERIFY_RAISES_ASSERT(concatenation = left.concatenate(right, 1));
  VERIFY_RAISES_ASSERT(concatenation = left.concatenate(right, 2));

  // Axis > NumDims or < 0.
  VERIFY_RAISES_ASSERT(concatenation = left.concatenate(right, 3));
  VERIFY_RAISES_ASSERT(concatenation = left.concatenate(right, -1));
}

template <int DataLayout>
static void test_static_dimension_failure() {
  Tensor<int, 2, DataLayout> left(2, 3);
  Tensor<int, 3, DataLayout> right(2, 3, 1);
  left.setRandom();
  right.setRandom();

  // TensorConcatenationOp requires both operands to have the same static rank.
  // To join tensors of different ranks, reshape one of the operands at the
  // call site; both directions are exercised here.
  Tensor<int, 3, DataLayout> concatenation = left.reshape(Tensor<int, 3>::Dimensions(2, 3, 1)).concatenate(right, 0);
  Tensor<int, 2, DataLayout> alternative = left.concatenate(right.reshape(Tensor<int, 2>::Dimensions(2, 3)), 0);
}

template <int DataLayout>
static void test_simple_concatenation() {
  Tensor<int, 3, DataLayout> left(2, 3, 1);
  Tensor<int, 3, DataLayout> right(2, 3, 1);
  left.setRandom();
  right.setRandom();

  Tensor<int, 3, DataLayout> concatenation = left.concatenate(right, 0);
  VERIFY_IS_EQUAL(concatenation.dimension(0), 4);
  VERIFY_IS_EQUAL(concatenation.dimension(1), 3);
  VERIFY_IS_EQUAL(concatenation.dimension(2), 1);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 2; ++i) {
      VERIFY_IS_EQUAL(concatenation(i, j, 0), left(i, j, 0));
    }
    for (int i = 2; i < 4; ++i) {
      VERIFY_IS_EQUAL(concatenation(i, j, 0), right(i - 2, j, 0));
    }
  }

  concatenation = left.concatenate(right, 1);
  VERIFY_IS_EQUAL(concatenation.dimension(0), 2);
  VERIFY_IS_EQUAL(concatenation.dimension(1), 6);
  VERIFY_IS_EQUAL(concatenation.dimension(2), 1);
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      VERIFY_IS_EQUAL(concatenation(i, j, 0), left(i, j, 0));
    }
    for (int j = 3; j < 6; ++j) {
      VERIFY_IS_EQUAL(concatenation(i, j, 0), right(i, j - 3, 0));
    }
  }

  concatenation = left.concatenate(right, 2);
  VERIFY_IS_EQUAL(concatenation.dimension(0), 2);
  VERIFY_IS_EQUAL(concatenation.dimension(1), 3);
  VERIFY_IS_EQUAL(concatenation.dimension(2), 2);
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      VERIFY_IS_EQUAL(concatenation(i, j, 0), left(i, j, 0));
      VERIFY_IS_EQUAL(concatenation(i, j, 1), right(i, j, 0));
    }
  }
}

// Exercise the packet() fast path when the concat axis is not the innermost
// dim and the inner dim is small enough that a packet load spans multiple
// rows -- including rows that fall on the right side of the boundary. The
// guard in packet() must reject this case and fall back to scalars.
template <int DataLayout>
static void test_concatenation_packet_axis_not_innermost() {
  // Output shape (8, 6, 1) with concat along axis 1: each packet load whose
  // first/last linear indices land on the left side will sweep through right
  // rows in between unless the fast path is correctly guarded.
  Tensor<float, 3, DataLayout> left(8, 3, 1);
  Tensor<float, 3, DataLayout> right(8, 3, 1);
  left.setRandom();
  right.setRandom();

  Tensor<float, 3, DataLayout> concatenation = left.concatenate(right, 1);
  VERIFY_IS_EQUAL(concatenation.dimension(0), 8);
  VERIFY_IS_EQUAL(concatenation.dimension(1), 6);
  VERIFY_IS_EQUAL(concatenation.dimension(2), 1);
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 3; ++j) {
      VERIFY_IS_EQUAL(concatenation(i, j, 0), left(i, j, 0));
      VERIFY_IS_EQUAL(concatenation(i, j + 3, 0), right(i, j, 0));
    }
  }

  // Force evaluation through the packet path with a coefficient-wise op so
  // the executor will request packets aligned to the output strides.
  Tensor<float, 3, DataLayout> doubled = concatenation * concatenation.constant(2.0f);
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 3; ++j) {
      VERIFY_IS_APPROX(doubled(i, j, 0), 2.0f * left(i, j, 0));
      VERIFY_IS_APPROX(doubled(i, j + 3, 0), 2.0f * right(i, j, 0));
    }
  }
}

static void test_concatenation_as_lvalue() {
  Tensor<int, 2> t1(2, 3);
  Tensor<int, 2> t2(2, 3);
  t1.setRandom();
  t2.setRandom();

  Tensor<int, 2> result(4, 3);
  result.setRandom();
  t1.concatenate(t2, 0) = result;

  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      VERIFY_IS_EQUAL(t1(i, j), result(i, j));
      VERIFY_IS_EQUAL(t2(i, j), result(i + 2, j));
    }
  }
}

// Regression tests: when a scalar-changing consumer (TensorCwiseUnaryOp,
// TensorConversionOp) sits above concat in an assign, the assign forwards a
// destination buffer sized for its *output* scalar (e.g. float in
// `abs(complex)`, double in `int.cast<double>()`). Before the producer-side
// fix in TensorCwiseUnaryOp::block / TensorConversionOp::block, concat's
// prepareStorage would reuse that buffer as its own (different) scalar,
// asserting in debug and corrupting output in release. These cases exercise
// the block path with the consumer dropping the buffer.

template <int DataLayout>
static void test_complex_concatenation_through_abs() {
  Tensor<std::complex<float>, 2, DataLayout> a(2, 3);
  Tensor<std::complex<float>, 2, DataLayout> b(2, 3);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 2; ++i) {
      a(i, j) = std::complex<float>(static_cast<float>(i + 1), static_cast<float>(j + 1));
      b(i, j) = std::complex<float>(static_cast<float>(i + 5), static_cast<float>(j + 2));
    }
  }

  Tensor<float, 2, DataLayout> out(4, 3);
  out = a.concatenate(b, 0).abs();

  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 2; ++i) {
      VERIFY_IS_APPROX(out(i, j), std::abs(a(i, j)));
      VERIFY_IS_APPROX(out(i + 2, j), std::abs(b(i, j)));
    }
  }
}

template <int DataLayout>
static void test_concatenation_through_cast() {
  Tensor<int, 2, DataLayout> a(2, 3);
  Tensor<int, 2, DataLayout> b(2, 3);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 2; ++i) {
      a(i, j) = i + 1 + 10 * j;
      b(i, j) = i + 5 + 10 * j;
    }
  }

  Tensor<double, 2, DataLayout> out(4, 3);
  out = a.concatenate(b, 0).template cast<double>();

  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 2; ++i) {
      VERIFY_IS_APPROX(out(i, j), static_cast<double>(a(i, j)));
      VERIFY_IS_APPROX(out(i + 2, j), static_cast<double>(b(i, j)));
    }
  }
}

EIGEN_DECLARE_TEST(tensor_concatenation) {
  CALL_SUBTEST(test_dimension_failures<ColMajor>());
  CALL_SUBTEST(test_dimension_failures<RowMajor>());
  CALL_SUBTEST(test_static_dimension_failure<ColMajor>());
  CALL_SUBTEST(test_static_dimension_failure<RowMajor>());
  CALL_SUBTEST(test_simple_concatenation<ColMajor>());
  CALL_SUBTEST(test_simple_concatenation<RowMajor>());
  CALL_SUBTEST(test_concatenation_packet_axis_not_innermost<ColMajor>());
  CALL_SUBTEST(test_concatenation_packet_axis_not_innermost<RowMajor>());
  CALL_SUBTEST(test_concatenation_as_lvalue());
  CALL_SUBTEST(test_complex_concatenation_through_abs<ColMajor>());
  CALL_SUBTEST(test_complex_concatenation_through_abs<RowMajor>());
  CALL_SUBTEST(test_concatenation_through_cast<ColMajor>());
  CALL_SUBTEST(test_concatenation_through_cast<RowMajor>());
}
