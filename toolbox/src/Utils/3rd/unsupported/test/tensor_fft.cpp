// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2014 Jianwei Cui <thucjw@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

#include "main.h"
#include <Eigen/Tensor>

using Eigen::Tensor;

template <int DataLayout>
static void test_fft_2D_golden() {
  Tensor<float, 2, DataLayout> input(2, 3);
  input(0, 0) = 1;
  input(0, 1) = 2;
  input(0, 2) = 3;
  input(1, 0) = 4;
  input(1, 1) = 5;
  input(1, 2) = 6;

  array<ptrdiff_t, 2> fft;
  fft[0] = 0;
  fft[1] = 1;

  Tensor<std::complex<float>, 2, DataLayout> output = input.template fft<Eigen::BothParts, Eigen::FFT_FORWARD>(fft);

  std::complex<float> output_golden[6];  // in ColMajor order
  output_golden[0] = std::complex<float>(21, 0);
  output_golden[1] = std::complex<float>(-9, 0);
  output_golden[2] = std::complex<float>(-3, 1.73205);
  output_golden[3] = std::complex<float>(0, 0);
  output_golden[4] = std::complex<float>(-3, -1.73205);
  output_golden[5] = std::complex<float>(0, 0);

  std::complex<float> c_offset = std::complex<float>(1.0, 1.0);

  if (DataLayout == ColMajor) {
    VERIFY_IS_APPROX(output(0) + c_offset, output_golden[0] + c_offset);
    VERIFY_IS_APPROX(output(1) + c_offset, output_golden[1] + c_offset);
    VERIFY_IS_APPROX(output(2) + c_offset, output_golden[2] + c_offset);
    VERIFY_IS_APPROX(output(3) + c_offset, output_golden[3] + c_offset);
    VERIFY_IS_APPROX(output(4) + c_offset, output_golden[4] + c_offset);
    VERIFY_IS_APPROX(output(5) + c_offset, output_golden[5] + c_offset);
  } else {
    VERIFY_IS_APPROX(output(0) + c_offset, output_golden[0] + c_offset);
    VERIFY_IS_APPROX(output(1) + c_offset, output_golden[2] + c_offset);
    VERIFY_IS_APPROX(output(2) + c_offset, output_golden[4] + c_offset);
    VERIFY_IS_APPROX(output(3) + c_offset, output_golden[1] + c_offset);
    VERIFY_IS_APPROX(output(4) + c_offset, output_golden[3] + c_offset);
    VERIFY_IS_APPROX(output(5) + c_offset, output_golden[5] + c_offset);
  }
}

static void test_fft_complex_input_golden() {
  Tensor<std::complex<float>, 1, ColMajor> input(5);
  input(0) = std::complex<float>(1, 1);
  input(1) = std::complex<float>(2, 2);
  input(2) = std::complex<float>(3, 3);
  input(3) = std::complex<float>(4, 4);
  input(4) = std::complex<float>(5, 5);

  array<ptrdiff_t, 1> fft;
  fft[0] = 0;

  Tensor<std::complex<float>, 1, ColMajor> forward_output_both_parts = input.fft<BothParts, FFT_FORWARD>(fft);
  Tensor<std::complex<float>, 1, ColMajor> reverse_output_both_parts = input.fft<BothParts, FFT_REVERSE>(fft);

  Tensor<float, 1, ColMajor> forward_output_real_part = input.fft<RealPart, FFT_FORWARD>(fft);
  Tensor<float, 1, ColMajor> reverse_output_real_part = input.fft<RealPart, FFT_REVERSE>(fft);

  Tensor<float, 1, ColMajor> forward_output_imag_part = input.fft<ImagPart, FFT_FORWARD>(fft);
  Tensor<float, 1, ColMajor> reverse_output_imag_part = input.fft<ImagPart, FFT_REVERSE>(fft);

  VERIFY_IS_EQUAL(forward_output_both_parts.dimension(0), input.dimension(0));
  VERIFY_IS_EQUAL(reverse_output_both_parts.dimension(0), input.dimension(0));

  VERIFY_IS_EQUAL(forward_output_real_part.dimension(0), input.dimension(0));
  VERIFY_IS_EQUAL(reverse_output_real_part.dimension(0), input.dimension(0));

  VERIFY_IS_EQUAL(forward_output_imag_part.dimension(0), input.dimension(0));
  VERIFY_IS_EQUAL(reverse_output_imag_part.dimension(0), input.dimension(0));

  std::complex<float> forward_golden_result[5];
  std::complex<float> reverse_golden_result[5];

  forward_golden_result[0] = std::complex<float>(15.000000000000000, +15.000000000000000);
  forward_golden_result[1] = std::complex<float>(-5.940954801177935, +0.940954801177934);
  forward_golden_result[2] = std::complex<float>(-3.312299240582266, -1.687700759417735);
  forward_golden_result[3] = std::complex<float>(-1.687700759417735, -3.312299240582266);
  forward_golden_result[4] = std::complex<float>(0.940954801177934, -5.940954801177935);

  reverse_golden_result[0] = std::complex<float>(3.000000000000000, +3.000000000000000);
  reverse_golden_result[1] = std::complex<float>(0.188190960235587, -1.188190960235587);
  reverse_golden_result[2] = std::complex<float>(-0.337540151883547, -0.662459848116453);
  reverse_golden_result[3] = std::complex<float>(-0.662459848116453, -0.337540151883547);
  reverse_golden_result[4] = std::complex<float>(-1.188190960235587, +0.188190960235587);

  for (int i = 0; i < 5; ++i) {
    VERIFY_IS_APPROX(forward_output_both_parts(i), forward_golden_result[i]);
    VERIFY_IS_APPROX(forward_output_real_part(i), forward_golden_result[i].real());
    VERIFY_IS_APPROX(forward_output_imag_part(i), forward_golden_result[i].imag());
  }

  for (int i = 0; i < 5; ++i) {
    VERIFY_IS_APPROX(reverse_output_both_parts(i), reverse_golden_result[i]);
    VERIFY_IS_APPROX(reverse_output_real_part(i), reverse_golden_result[i].real());
    VERIFY_IS_APPROX(reverse_output_imag_part(i), reverse_golden_result[i].imag());
  }
}

static void test_fft_real_input_golden() {
  Tensor<float, 1, ColMajor> input(5);
  input(0) = 1.0;
  input(1) = 2.0;
  input(2) = 3.0;
  input(3) = 4.0;
  input(4) = 5.0;

  array<ptrdiff_t, 1> fft;
  fft[0] = 0;

  Tensor<std::complex<float>, 1, ColMajor> forward_output_both_parts = input.fft<BothParts, FFT_FORWARD>(fft);
  Tensor<std::complex<float>, 1, ColMajor> reverse_output_both_parts = input.fft<BothParts, FFT_REVERSE>(fft);

  Tensor<float, 1, ColMajor> forward_output_real_part = input.fft<RealPart, FFT_FORWARD>(fft);
  Tensor<float, 1, ColMajor> reverse_output_real_part = input.fft<RealPart, FFT_REVERSE>(fft);

  Tensor<float, 1, ColMajor> forward_output_imag_part = input.fft<ImagPart, FFT_FORWARD>(fft);
  Tensor<float, 1, ColMajor> reverse_output_imag_part = input.fft<ImagPart, FFT_REVERSE>(fft);

  VERIFY_IS_EQUAL(forward_output_both_parts.dimension(0), input.dimension(0));
  VERIFY_IS_EQUAL(reverse_output_both_parts.dimension(0), input.dimension(0));

  VERIFY_IS_EQUAL(forward_output_real_part.dimension(0), input.dimension(0));
  VERIFY_IS_EQUAL(reverse_output_real_part.dimension(0), input.dimension(0));

  VERIFY_IS_EQUAL(forward_output_imag_part.dimension(0), input.dimension(0));
  VERIFY_IS_EQUAL(reverse_output_imag_part.dimension(0), input.dimension(0));

  std::complex<float> forward_golden_result[5];
  std::complex<float> reverse_golden_result[5];

  forward_golden_result[0] = std::complex<float>(15, 0);
  forward_golden_result[1] = std::complex<float>(-2.5, +3.44095480117793);
  forward_golden_result[2] = std::complex<float>(-2.5, +0.81229924058227);
  forward_golden_result[3] = std::complex<float>(-2.5, -0.81229924058227);
  forward_golden_result[4] = std::complex<float>(-2.5, -3.44095480117793);

  reverse_golden_result[0] = std::complex<float>(3.0, 0);
  reverse_golden_result[1] = std::complex<float>(-0.5, -0.688190960235587);
  reverse_golden_result[2] = std::complex<float>(-0.5, -0.162459848116453);
  reverse_golden_result[3] = std::complex<float>(-0.5, +0.162459848116453);
  reverse_golden_result[4] = std::complex<float>(-0.5, +0.688190960235587);

  std::complex<float> c_offset(1.0, 1.0);
  float r_offset = 1.0;

  for (int i = 0; i < 5; ++i) {
    VERIFY_IS_APPROX(forward_output_both_parts(i) + c_offset, forward_golden_result[i] + c_offset);
    VERIFY_IS_APPROX(forward_output_real_part(i) + r_offset, forward_golden_result[i].real() + r_offset);
    VERIFY_IS_APPROX(forward_output_imag_part(i) + r_offset, forward_golden_result[i].imag() + r_offset);
  }

  for (int i = 0; i < 5; ++i) {
    VERIFY_IS_APPROX(reverse_output_both_parts(i) + c_offset, reverse_golden_result[i] + c_offset);
    VERIFY_IS_APPROX(reverse_output_real_part(i) + r_offset, reverse_golden_result[i].real() + r_offset);
    VERIFY_IS_APPROX(reverse_output_imag_part(i) + r_offset, reverse_golden_result[i].imag() + r_offset);
  }
}

template <int DataLayout, typename RealScalar, bool isComplexInput, int FFTResultType, int FFTDirection, int TensorRank>
static void test_fft_real_input_energy() {
  Eigen::DSizes<ptrdiff_t, TensorRank> dimensions;
  ptrdiff_t total_size = 1;
  for (int i = 0; i < TensorRank; ++i) {
    dimensions[i] = rand() % 20 + 1;
    total_size *= dimensions[i];
  }
  const DSizes<ptrdiff_t, TensorRank> arr = dimensions;

  typedef std::conditional_t<isComplexInput == true, std::complex<RealScalar>, RealScalar> InputScalar;

  Tensor<InputScalar, TensorRank, DataLayout> input;
  input.resize(arr);
  input.setRandom();

  array<ptrdiff_t, TensorRank> fft;
  for (int i = 0; i < TensorRank; ++i) {
    fft[i] = i;
  }

  typedef std::conditional_t<FFTResultType == Eigen::BothParts, std::complex<RealScalar>, RealScalar> OutputScalar;
  Tensor<OutputScalar, TensorRank, DataLayout> output;
  output = input.template fft<FFTResultType, FFTDirection>(fft);

  for (int i = 0; i < TensorRank; ++i) {
    VERIFY_IS_EQUAL(output.dimension(i), input.dimension(i));
  }

  RealScalar energy_original = 0.0;
  RealScalar energy_after_fft = 0.0;

  for (int i = 0; i < total_size; ++i) {
    energy_original += numext::abs2(input(i));
  }

  for (int i = 0; i < total_size; ++i) {
    energy_after_fft += numext::abs2(output(i));
  }

  if (FFTDirection == FFT_FORWARD) {
    VERIFY_IS_APPROX(energy_original, energy_after_fft / total_size);
  } else {
    VERIFY_IS_APPROX(energy_original, energy_after_fft * total_size);
  }
}

template <typename RealScalar>
static void test_fft_non_power_of_2_round_trip(int exponent) {
  int n = (1 << exponent) + 1;

  Eigen::DSizes<ptrdiff_t, 1> dimensions;
  dimensions[0] = n;
  const DSizes<ptrdiff_t, 1> arr = dimensions;
  Tensor<RealScalar, 1, ColMajor, ptrdiff_t> input;

  input.resize(arr);
  input.setRandom();

  array<int, 1> fft;
  fft[0] = 0;

  Tensor<std::complex<RealScalar>, 1, ColMajor> forward = input.template fft<BothParts, FFT_FORWARD>(fft);

  Tensor<RealScalar, 1, ColMajor, ptrdiff_t> output = forward.template fft<RealPart, FFT_REVERSE>(fft);

  for (int i = 0; i < n; ++i) {
    RealScalar tol = test_precision<RealScalar>() * (std::abs(input[i]) + std::abs(output[i]) + 1);
    VERIFY_IS_APPROX_OR_LESS_THAN(std::abs(input[i] - output[i]), tol);
  }
}

// Documented precondition: finite inputs must produce finite outputs. The hot-
// path complex multiply uses the naive form (no NaN/inf disambiguation), so
// this verifies the optimization didn't introduce spurious NaN/inf on legal
// inputs. Covers both Cooley-Tukey (power-of-2) and Bluestein (non-pow-2).
template <typename RealScalar>
static void test_fft_finite_input_stays_finite(int n) {
  using Complex = std::complex<RealScalar>;
  Tensor<Complex, 1, ColMajor> input(n);
  // Large but finite magnitudes — close to overflow but not at it — to catch
  // any cmul that produces inf - inf = NaN under unfavorable scheduling.
  const RealScalar big = RealScalar(1e18);
  for (int i = 0; i < n; ++i) {
    input(i) = Complex(big * std::cos(RealScalar(i)), big * std::sin(RealScalar(0.5) * RealScalar(i)));
  }
  array<int, 1> fft = {0};
  Tensor<Complex, 1, ColMajor> fwd = input.template fft<Eigen::BothParts, FFT_FORWARD>(fft);
  for (int i = 0; i < n; ++i) {
    VERIFY((numext::isfinite)(fwd(i).real()));
    VERIFY((numext::isfinite)(fwd(i).imag()));
  }
  Tensor<Complex, 1, ColMajor> rev = fwd.template fft<Eigen::BothParts, FFT_REVERSE>(fft);
  for (int i = 0; i < n; ++i) {
    VERIFY((numext::isfinite)(rev(i).real()));
    VERIFY((numext::isfinite)(rev(i).imag()));
  }
}

// `x = x.fft(...)` aliases the FFT input and output onto the same storage.
// Regression test for the memcpy fast-path's self-copy guard in evalToBuf.
template <typename RealScalar>
static void test_fft_in_place_assign(int n) {
  using Complex = std::complex<RealScalar>;
  Tensor<Complex, 1, ColMajor> reference(n);
  reference.setRandom();
  Tensor<Complex, 1, ColMajor> x = reference;
  array<int, 1> fft = {0};
  x = x.template fft<Eigen::BothParts, FFT_FORWARD>(fft);
  x = x.template fft<Eigen::BothParts, FFT_REVERSE>(fft);
  for (int i = 0; i < n; ++i) {
    RealScalar tol = test_precision<RealScalar>() * (std::abs(reference(i)) + std::abs(x(i)) + RealScalar(1));
    VERIFY_IS_APPROX_OR_LESS_THAN(std::abs(reference(i) - x(i)) / RealScalar(n), tol);
  }
}

EIGEN_DECLARE_TEST(tensor_fft) {
  test_fft_complex_input_golden();
  test_fft_real_input_golden();

  test_fft_2D_golden<ColMajor>();
  test_fft_2D_golden<RowMajor>();

  test_fft_real_input_energy<ColMajor, float, true, Eigen::BothParts, FFT_FORWARD, 1>();
  test_fft_real_input_energy<ColMajor, double, true, Eigen::BothParts, FFT_FORWARD, 1>();
  test_fft_real_input_energy<ColMajor, float, false, Eigen::BothParts, FFT_FORWARD, 1>();
  test_fft_real_input_energy<ColMajor, double, false, Eigen::BothParts, FFT_FORWARD, 1>();

  test_fft_real_input_energy<ColMajor, float, true, Eigen::BothParts, FFT_FORWARD, 2>();
  test_fft_real_input_energy<ColMajor, double, true, Eigen::BothParts, FFT_FORWARD, 2>();
  test_fft_real_input_energy<ColMajor, float, false, Eigen::BothParts, FFT_FORWARD, 2>();
  test_fft_real_input_energy<ColMajor, double, false, Eigen::BothParts, FFT_FORWARD, 2>();

  test_fft_real_input_energy<ColMajor, float, true, Eigen::BothParts, FFT_FORWARD, 3>();
  test_fft_real_input_energy<ColMajor, double, true, Eigen::BothParts, FFT_FORWARD, 3>();
  test_fft_real_input_energy<ColMajor, float, false, Eigen::BothParts, FFT_FORWARD, 3>();
  test_fft_real_input_energy<ColMajor, double, false, Eigen::BothParts, FFT_FORWARD, 3>();

  test_fft_real_input_energy<ColMajor, float, true, Eigen::BothParts, FFT_FORWARD, 4>();
  test_fft_real_input_energy<ColMajor, double, true, Eigen::BothParts, FFT_FORWARD, 4>();
  test_fft_real_input_energy<ColMajor, float, false, Eigen::BothParts, FFT_FORWARD, 4>();
  test_fft_real_input_energy<ColMajor, double, false, Eigen::BothParts, FFT_FORWARD, 4>();

  test_fft_real_input_energy<RowMajor, float, true, Eigen::BothParts, FFT_FORWARD, 1>();
  test_fft_real_input_energy<RowMajor, double, true, Eigen::BothParts, FFT_FORWARD, 1>();
  test_fft_real_input_energy<RowMajor, float, false, Eigen::BothParts, FFT_FORWARD, 1>();
  test_fft_real_input_energy<RowMajor, double, false, Eigen::BothParts, FFT_FORWARD, 1>();

  test_fft_real_input_energy<RowMajor, float, true, Eigen::BothParts, FFT_FORWARD, 2>();
  test_fft_real_input_energy<RowMajor, double, true, Eigen::BothParts, FFT_FORWARD, 2>();
  test_fft_real_input_energy<RowMajor, float, false, Eigen::BothParts, FFT_FORWARD, 2>();
  test_fft_real_input_energy<RowMajor, double, false, Eigen::BothParts, FFT_FORWARD, 2>();

  test_fft_real_input_energy<RowMajor, float, true, Eigen::BothParts, FFT_FORWARD, 3>();
  test_fft_real_input_energy<RowMajor, double, true, Eigen::BothParts, FFT_FORWARD, 3>();
  test_fft_real_input_energy<RowMajor, float, false, Eigen::BothParts, FFT_FORWARD, 3>();
  test_fft_real_input_energy<RowMajor, double, false, Eigen::BothParts, FFT_FORWARD, 3>();

  test_fft_real_input_energy<RowMajor, float, true, Eigen::BothParts, FFT_FORWARD, 4>();
  test_fft_real_input_energy<RowMajor, double, true, Eigen::BothParts, FFT_FORWARD, 4>();
  test_fft_real_input_energy<RowMajor, float, false, Eigen::BothParts, FFT_FORWARD, 4>();
  test_fft_real_input_energy<RowMajor, double, false, Eigen::BothParts, FFT_FORWARD, 4>();

  test_fft_non_power_of_2_round_trip<float>(7);
  test_fft_non_power_of_2_round_trip<double>(7);

  // Pow-2 (Cooley-Tukey) and non-pow-2 (Bluestein) finite-input precondition.
  test_fft_finite_input_stays_finite<float>(64);
  test_fft_finite_input_stays_finite<float>(100);
  test_fft_finite_input_stays_finite<double>(64);
  test_fft_finite_input_stays_finite<double>(100);

  // x = x.fft(...) — exercises the self-copy guard in evalToBuf.
  test_fft_in_place_assign<float>(64);
  test_fft_in_place_assign<float>(100);
  test_fft_in_place_assign<double>(64);
  test_fft_in_place_assign<double>(100);
}
