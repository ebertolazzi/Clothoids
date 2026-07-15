// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#include "main.h"

#include <Eigen/Tensor>

using DynamicDims = DSizes<Index, 3>;
using StaticDims = Sizes<2, 3, 4>;
using IndexTuple = internal::IndexTuple<Index, type2index<3>>;
using DynamicIndexList = IndexList<Index, type2index<3>>;

static_assert(std::is_trivially_destructible<TensorOpCost>::value, "TensorOpCost should have a trivial destructor");
static_assert(std::is_move_constructible<TensorOpCost>::value, "TensorOpCost should be move constructible");
static_assert(std::is_move_constructible<DynamicDims>::value, "DSizes should be move constructible");
static_assert(std::is_trivially_destructible<DynamicDims>::value, "DSizes should have a trivial destructor");
static_assert(std::is_move_constructible<StaticDims>::value, "Sizes should be move constructible");
static_assert(std::is_trivially_destructible<StaticDims>::value, "Sizes should have a trivial destructor");
static_assert(std::is_move_constructible<Pair<int, int>>::value, "Pair should be move constructible");
static_assert(std::is_trivially_destructible<Pair<int, int>>::value, "Pair should have a trivial destructor");
static_assert(std::is_move_constructible<IndexPair<Index>>::value, "IndexPair should be move constructible");
static_assert(std::is_trivially_destructible<IndexPair<Index>>::value, "IndexPair should have a trivial destructor");
static_assert(std::is_move_constructible<IndexTuple>::value, "IndexTuple should be move constructible");
static_assert(std::is_trivially_destructible<IndexTuple>::value, "IndexTuple should have a trivial destructor");
static_assert(std::is_move_constructible<DynamicIndexList>::value, "IndexList should be move constructible");
static_assert(std::is_trivially_destructible<DynamicIndexList>::value, "IndexList should have a trivial destructor");

constexpr TensorOpCost default_cost;
static_assert(default_cost.bytes_loaded() == 0.0, "TensorOpCost default bytes_loaded should be zero");
static_assert(default_cost.bytes_stored() == 0.0, "TensorOpCost default bytes_stored should be zero");
static_assert(default_cost.compute_cycles() == 0.0, "TensorOpCost default compute_cycles should be zero");

EIGEN_DECLARE_TEST(tensor_special_members) {}
