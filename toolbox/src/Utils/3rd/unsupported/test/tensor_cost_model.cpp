// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 The Eigen Authors.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

#include "main.h"

#include <Eigen/Tensor>

using Eigen::DefaultDevice;
using Eigen::Tensor;
using Eigen::TensorEvaluator;

// Regression: TensorCwiseNullaryOp must report bytes_loaded == 0.
// NullaryOps (constants, Zero, Identity, Random, sequence generators)
// produce values from registers or minimal state without loading from
// memory. If they reported nonzero bytes_loaded, expressions dominated
// by constants (e.g. Horner-form polynomials) would be misclassified as
// memory-bound and the threadpool cost model would over-restrict
// parallelism. See TensorEvaluator.h, the CwiseNullaryOp specialization
// of costPerCoeff().
template <typename Scalar>
static void test_nullary_zero_bytes_loaded() {
  Tensor<Scalar, 1> shape(/*size=*/16);
  auto zeros = shape.constant(Scalar(0));
  auto sevens = shape.constant(Scalar(7));

  using ZeroEval = TensorEvaluator<const decltype(zeros), DefaultDevice>;
  using ConstEval = TensorEvaluator<const decltype(sevens), DefaultDevice>;

  DefaultDevice device;
  ZeroEval zero_eval(zeros, device);
  ConstEval const_eval(sevens, device);

  for (bool vectorized : {false, true}) {
    const auto zero_cost = zero_eval.costPerCoeff(vectorized);
    const auto const_cost = const_eval.costPerCoeff(vectorized);
    VERIFY_IS_EQUAL(zero_cost.bytes_loaded(), 0.0);
    VERIFY_IS_EQUAL(zero_cost.bytes_stored(), 0.0);
    VERIFY_IS_EQUAL(const_cost.bytes_loaded(), 0.0);
    VERIFY_IS_EQUAL(const_cost.bytes_stored(), 0.0);
  }
}

EIGEN_DECLARE_TEST(tensor_cost_model) {
  CALL_SUBTEST(test_nullary_zero_bytes_loaded<float>());
  CALL_SUBTEST(test_nullary_zero_bytes_loaded<double>());
  CALL_SUBTEST(test_nullary_zero_bytes_loaded<int>());
}
