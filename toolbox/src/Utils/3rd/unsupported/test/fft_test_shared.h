// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2009 Mark Borgerding mark a borgerding net
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

#ifndef EIGEN_UNSUPPORTED_TEST_FFT_TEST_SHARED_H
#define EIGEN_UNSUPPORTED_TEST_FFT_TEST_SHARED_H

// Enable runtime malloc tracking so test_inplace_complex<>() can assert the
// in-place scratch comes from the stack.  Allocation defaults to allowed; the
// tracking only fires inside explicit set_is_malloc_allowed(false) windows.
#define EIGEN_RUNTIME_NO_MALLOC

#include "main.h"
#include <unsupported/Eigen/FFT>

template <typename T>
inline std::complex<T> RandomCpx() {
  return std::complex<T>((T)(rand() / (T)RAND_MAX - .5), (T)(rand() / (T)RAND_MAX - .5));
}

using namespace std;
using namespace Eigen;

template <typename T>
inline complex<long double> promote(complex<T> x) {
  return complex<long double>((long double)x.real(), (long double)x.imag());
}

inline complex<long double> promote(float x) { return complex<long double>((long double)x); }
inline complex<long double> promote(double x) { return complex<long double>((long double)x); }
inline complex<long double> promote(long double x) { return complex<long double>((long double)x); }

template <typename VT1, typename VT2>
long double fft_rmse(const VT1& fftbuf, const VT2& timebuf) {
  long double totalpower = 0;
  long double difpower = 0;
  long double pi = acos((long double)-1);
  for (size_t k0 = 0; k0 < (size_t)fftbuf.size(); ++k0) {
    complex<long double> acc = 0;
    long double phinc = (long double)(-2.) * k0 * pi / timebuf.size();
    for (size_t k1 = 0; k1 < (size_t)timebuf.size(); ++k1) {
      acc += promote(timebuf[k1]) * exp(complex<long double>(0, k1 * phinc));
    }
    totalpower += numext::abs2(acc);
    complex<long double> x = promote(fftbuf[k0]);
    complex<long double> dif = acc - x;
    difpower += numext::abs2(dif);
    // cerr << k0 << "\t" << acc << "\t" <<  x << "\t" << sqrt(numext::abs2(dif)) << endl;
  }
  // cerr << "rmse:" << sqrt(difpower/totalpower) << endl;
  return sqrt(difpower / totalpower);
}

template <typename VT1, typename VT2>
long double dif_rmse(const VT1 buf1, const VT2 buf2) {
  long double totalpower = 0;
  long double difpower = 0;
  size_t n = (min)(buf1.size(), buf2.size());
  for (size_t k = 0; k < n; ++k) {
    totalpower += (long double)((numext::abs2(buf1[k]) + numext::abs2(buf2[k])) / 2);
    difpower += (long double)(numext::abs2(buf1[k] - buf2[k]));
  }
  return sqrt(difpower / totalpower);
}

enum {
  StdVectorContainer,
  EigenVectorContainer,
  EigenRowVectorContainer,
  EigenArrayXContainer,
  EigenRowArrayXContainer
};

template <int Container, typename Scalar>
struct VectorType;

template <typename Scalar>
struct VectorType<StdVectorContainer, Scalar> {
  typedef vector<Scalar> type;
};

template <typename Scalar>
struct VectorType<EigenVectorContainer, Scalar> {
  typedef Matrix<Scalar, Dynamic, 1> type;
};

template <typename Scalar>
struct VectorType<EigenRowVectorContainer, Scalar> {
  typedef Matrix<Scalar, 1, Dynamic> type;
};

template <typename Scalar>
struct VectorType<EigenArrayXContainer, Scalar> {
  typedef Array<Scalar, Dynamic, 1> type;
};

template <typename Scalar>
struct VectorType<EigenRowArrayXContainer, Scalar> {
  typedef Array<Scalar, 1, Dynamic> type;
};

template <int ScalarContainer, int ComplexContainer, typename T>
void test_scalar_generic(int nfft) {
  typedef typename FFT<T>::Complex Complex;
  typedef typename FFT<T>::Scalar Scalar;
  typedef typename VectorType<ScalarContainer, Scalar>::type ScalarVector;
  typedef typename VectorType<ComplexContainer, Complex>::type ComplexVector;

  FFT<T> fft;
  ScalarVector tbuf(nfft);
  ComplexVector freqBuf;
  for (int k = 0; k < nfft; ++k) tbuf[k] = (T)(rand() / (double)RAND_MAX - .5);

  // make sure it DOESN'T give the right full spectrum answer
  // if we've asked for half-spectrum
  fft.SetFlag(fft.HalfSpectrum);
  fft.fwd(freqBuf, tbuf);
  VERIFY((size_t)freqBuf.size() == (size_t)((nfft >> 1) + 1));
  VERIFY(T(fft_rmse(freqBuf, tbuf)) < test_precision<T>());  // gross check

  fft.ClearFlag(fft.HalfSpectrum);
  fft.fwd(freqBuf, tbuf);
  VERIFY((size_t)freqBuf.size() == (size_t)nfft);
  VERIFY(T(fft_rmse(freqBuf, tbuf)) < test_precision<T>());  // gross check

  if (nfft & 1) return;  // odd FFTs get the wrong size inverse FFT

  ScalarVector tbuf2;
  fft.inv(tbuf2, freqBuf);
  VERIFY(T(dif_rmse(tbuf, tbuf2)) < test_precision<T>());  // gross check

  // verify that the Unscaled flag takes effect
  ScalarVector tbuf3;
  fft.SetFlag(fft.Unscaled);

  fft.inv(tbuf3, freqBuf);

  for (int k = 0; k < nfft; ++k) tbuf3[k] *= T(1. / nfft);

  // for (size_t i=0;i<(size_t) tbuf.size();++i)
  //     cout << "freqBuf=" << freqBuf[i] << " in2=" << tbuf3[i] << " -  in=" << tbuf[i] << " => " << (tbuf3[i] -
  //     tbuf[i] ) <<  endl;

  VERIFY(T(dif_rmse(tbuf, tbuf3)) < test_precision<T>());  // gross check

  // verify that ClearFlag works
  fft.ClearFlag(fft.Unscaled);
  fft.inv(tbuf2, freqBuf);
  VERIFY(T(dif_rmse(tbuf, tbuf2)) < test_precision<T>());  // gross check
}

template <typename T>
void test_scalar(int nfft) {
  // std:vector is a special case that does not interact with DenseBase types
  test_scalar_generic<StdVectorContainer, StdVectorContainer, T>(nfft);

  // All Dense types as dst and src in various combinations
  test_scalar_generic<EigenVectorContainer, EigenVectorContainer, T>(nfft);
  test_scalar_generic<EigenVectorContainer, EigenRowVectorContainer, T>(nfft);
  test_scalar_generic<EigenVectorContainer, EigenArrayXContainer, T>(nfft);
  test_scalar_generic<EigenVectorContainer, EigenRowArrayXContainer, T>(nfft);

  test_scalar_generic<EigenRowVectorContainer, EigenVectorContainer, T>(nfft);
  test_scalar_generic<EigenRowVectorContainer, EigenRowVectorContainer, T>(nfft);
  test_scalar_generic<EigenRowVectorContainer, EigenArrayXContainer, T>(nfft);
  test_scalar_generic<EigenRowVectorContainer, EigenRowArrayXContainer, T>(nfft);

  test_scalar_generic<EigenArrayXContainer, EigenVectorContainer, T>(nfft);
  test_scalar_generic<EigenArrayXContainer, EigenRowVectorContainer, T>(nfft);
  test_scalar_generic<EigenArrayXContainer, EigenArrayXContainer, T>(nfft);
  test_scalar_generic<EigenArrayXContainer, EigenRowArrayXContainer, T>(nfft);

  test_scalar_generic<EigenRowArrayXContainer, EigenVectorContainer, T>(nfft);
  test_scalar_generic<EigenRowArrayXContainer, EigenRowVectorContainer, T>(nfft);
  test_scalar_generic<EigenRowArrayXContainer, EigenArrayXContainer, T>(nfft);
  test_scalar_generic<EigenRowArrayXContainer, EigenRowArrayXContainer, T>(nfft);
}

template <int ContainerA, int ContainerB, typename T>
void test_complex_generic(int nfft) {
  typedef typename FFT<T>::Complex Complex;
  typedef typename VectorType<ContainerA, Complex>::type ComplexVectorA;
  typedef typename VectorType<ContainerB, Complex>::type ComplexVectorB;

  FFT<T> fft;

  ComplexVectorB inbuf(nfft);
  ComplexVectorA outbuf;
  ComplexVectorB buf3;
  for (int k = 0; k < nfft; ++k)
    inbuf[k] = Complex((T)(rand() / (double)RAND_MAX - .5), (T)(rand() / (double)RAND_MAX - .5));
  fft.fwd(outbuf, inbuf);

  VERIFY(T(fft_rmse(outbuf, inbuf)) < test_precision<T>());  // gross check
  fft.inv(buf3, outbuf);

  VERIFY(T(dif_rmse(inbuf, buf3)) < test_precision<T>());  // gross check

  // verify that the Unscaled flag takes effect
  ComplexVectorA buf4;
  fft.SetFlag(fft.Unscaled);
  fft.inv(buf4, outbuf);
  for (int k = 0; k < nfft; ++k) buf4[k] *= T(1. / nfft);
  VERIFY(T(dif_rmse(inbuf, buf4)) < test_precision<T>());  // gross check

  // verify that ClearFlag works
  fft.ClearFlag(fft.Unscaled);
  fft.inv(buf3, outbuf);
  VERIFY(T(dif_rmse(inbuf, buf3)) < test_precision<T>());  // gross check
}

template <typename T>
void test_complex_strided(int nfft) {
  typedef typename FFT<T>::Complex Complex;
  typedef typename Eigen::Vector<Complex, Dynamic> ComplexVector;
  constexpr int kInputStride = 3;
  constexpr int kOutputStride = 7;
  constexpr int kInvOutputStride = 13;

  FFT<T> fft;

  ComplexVector inbuf(nfft * kInputStride);
  inbuf.setRandom();
  ComplexVector outbuf(nfft * kOutputStride);
  outbuf.setRandom();
  ComplexVector invoutbuf(nfft * kInvOutputStride);
  invoutbuf.setRandom();

  using StridedComplexVector = Map<ComplexVector, /*MapOptions=*/0, InnerStride<Dynamic>>;
  StridedComplexVector input(inbuf.data(), nfft, InnerStride<Dynamic>(kInputStride));
  StridedComplexVector output(outbuf.data(), nfft, InnerStride<Dynamic>(kOutputStride));
  StridedComplexVector inv_output(invoutbuf.data(), nfft, InnerStride<Dynamic>(kInvOutputStride));

  for (int k = 0; k < nfft; ++k)
    input[k] = Complex((T)(rand() / (double)RAND_MAX - .5), (T)(rand() / (double)RAND_MAX - .5));
  fft.fwd(output, input);

  VERIFY(T(fft_rmse(output, input)) < test_precision<T>());  // gross check
  fft.inv(inv_output, output);
  VERIFY(T(dif_rmse(inv_output, input)) < test_precision<T>());  // gross check
}

template <typename T>
void test_complex(int nfft) {
  // std:vector is a special case that does not interact with DenseBase types
  test_complex_generic<StdVectorContainer, StdVectorContainer, T>(nfft);

  // All Dense types as dst and src in various combinations
  test_complex_generic<EigenVectorContainer, EigenVectorContainer, T>(nfft);
  test_complex_generic<EigenVectorContainer, EigenRowVectorContainer, T>(nfft);
  test_complex_generic<EigenVectorContainer, EigenArrayXContainer, T>(nfft);
  test_complex_generic<EigenVectorContainer, EigenRowArrayXContainer, T>(nfft);

  test_complex_generic<EigenRowVectorContainer, EigenVectorContainer, T>(nfft);
  test_complex_generic<EigenRowVectorContainer, EigenRowVectorContainer, T>(nfft);
  test_complex_generic<EigenRowVectorContainer, EigenArrayXContainer, T>(nfft);
  test_complex_generic<EigenRowVectorContainer, EigenRowArrayXContainer, T>(nfft);

  test_complex_generic<EigenArrayXContainer, EigenVectorContainer, T>(nfft);
  test_complex_generic<EigenArrayXContainer, EigenRowVectorContainer, T>(nfft);
  test_complex_generic<EigenArrayXContainer, EigenArrayXContainer, T>(nfft);
  test_complex_generic<EigenArrayXContainer, EigenRowArrayXContainer, T>(nfft);

  test_complex_generic<EigenRowArrayXContainer, EigenVectorContainer, T>(nfft);
  test_complex_generic<EigenRowArrayXContainer, EigenRowVectorContainer, T>(nfft);
  test_complex_generic<EigenRowArrayXContainer, EigenArrayXContainer, T>(nfft);
  test_complex_generic<EigenRowArrayXContainer, EigenRowArrayXContainer, T>(nfft);

  test_complex_strided<T>(nfft);
}

template <typename T, int nrows, int ncols>
void test_complex2d() {
  typedef typename Eigen::FFT<T>::Complex Complex;
  FFT<T> fft;
  Eigen::Matrix<Complex, nrows, ncols> src, src2, dst, dst2;

  src = Eigen::Matrix<Complex, nrows, ncols>::Random();
  // src =  Eigen::Matrix<Complex,nrows,ncols>::Identity();

  for (int k = 0; k < ncols; k++) {
    Eigen::Matrix<Complex, nrows, 1> tmpOut;
    fft.fwd(tmpOut, src.col(k));
    dst2.col(k) = tmpOut;
  }

  for (int k = 0; k < nrows; k++) {
    Eigen::Matrix<Complex, 1, ncols> tmpOut;
    fft.fwd(tmpOut, dst2.row(k));
    dst2.row(k) = tmpOut;
  }

  fft.fwd2(dst.data(), src.data(), ncols, nrows);
  fft.inv2(src2.data(), dst.data(), ncols, nrows);
  VERIFY((src - src2).norm() < test_precision<T>());
  VERIFY((dst - dst2).norm() < test_precision<T>());
}

// Regression for issue #868: fft.fwd(buf, buf) / fft.inv(buf, buf) with the
// same buffer as input and output must produce the out-of-place result.
// Also pins down that the in-place scratch comes from the stack for typical
// sizes so EIGEN_RUNTIME_NO_MALLOC users aren't forced to heap-allocate.
template <typename T>
void test_inplace_complex(int nfft) {
  typedef typename FFT<T>::Complex Complex;
  typedef Matrix<Complex, Dynamic, 1> ComplexVector;

  ComplexVector in(nfft);
  for (int k = 0; k < nfft; ++k)
    in[k] = Complex((T)(rand() / (double)RAND_MAX - .5), (T)(rand() / (double)RAND_MAX - .5));

  FFT<T> fft;
  ComplexVector out_ref;
  fft.fwd(out_ref, in);
  ComplexVector inv_ref;
  fft.inv(inv_ref, out_ref);

  ComplexVector inout_fwd = in;
  ComplexVector inout_inv = out_ref;

  Eigen::internal::set_is_malloc_allowed(false);
  fft.fwd(inout_fwd, inout_fwd);
  fft.inv(inout_inv, inout_inv);
  Eigen::internal::set_is_malloc_allowed(true);

  VERIFY((out_ref - inout_fwd).cwiseAbs().maxCoeff() < test_precision<T>());
  VERIFY((inv_ref - inout_inv).cwiseAbs().maxCoeff() < test_precision<T>());
}

// Regression for issue #675: zero-padding a fixed-size vector must not assume
// a column block in the temporary contiguous FFT input.
template <typename T>
void test_fwd_padding(int nfft) {
  typedef typename FFT<T>::Complex Complex;
  typedef Matrix<T, 10, 1> FixedRealColumn;
  typedef Matrix<T, 1, 10> FixedRealRow;
  typedef Matrix<Complex, 10, 1> FixedComplexColumn;
  typedef Matrix<Complex, 1, 10> FixedComplexRow;
  typedef Matrix<T, Dynamic, 1> RealVector;
  typedef Matrix<Complex, Dynamic, 1> ComplexVector;

  FixedRealColumn real_column;
  FixedComplexColumn complex_column;
  for (int k = 0; k < real_column.size(); ++k) {
    real_column[k] = (T)(rand() / (double)RAND_MAX - .5);
    complex_column[k] = Complex((T)(rand() / (double)RAND_MAX - .5), (T)(rand() / (double)RAND_MAX - .5));
  }
  FixedRealRow real_row = real_column.transpose();
  FixedComplexRow complex_row = complex_column.transpose();

  RealVector real_padded = RealVector::Zero(nfft);
  real_padded.head(real_column.size()) = real_column;
  ComplexVector complex_padded = ComplexVector::Zero(nfft);
  complex_padded.head(complex_column.size()) = complex_column;

  FFT<T> fft;
  ComplexVector expected;
  ComplexVector actual;

  fft.fwd(expected, real_padded);
  fft.fwd(actual, real_column, nfft);
  VERIFY_IS_APPROX(actual, expected);
  fft.fwd(actual, real_row, nfft);
  VERIFY_IS_APPROX(actual, expected);

  fft.fwd(expected, complex_padded);
  fft.fwd(actual, complex_column, nfft);
  VERIFY_IS_APPROX(actual, expected);
  fft.fwd(actual, complex_row, nfft);
  VERIFY_IS_APPROX(actual, expected);

  fft.SetFlag(fft.HalfSpectrum);
  fft.fwd(expected, real_padded);
  fft.fwd(actual, real_column, nfft);
  VERIFY_IS_APPROX(actual, expected);
  fft.fwd(actual, real_row, nfft);
  VERIFY_IS_APPROX(actual, expected);
}

inline void test_return_by_value(int len) {
  VectorXf in;
  VectorXf in1;
  in.setRandom(len);
  VectorXcf out1, out2;
  FFT<float> fft;

  fft.SetFlag(fft.HalfSpectrum);

  fft.fwd(out1, in);
  out2 = fft.fwd(in);
  VERIFY((out1 - out2).norm() < test_precision<float>());
  in1 = fft.inv(out1);
  VERIFY((in1 - in).norm() < test_precision<float>());
}

// Regression for issue #1537: reusing the same FFT object across real-input
// and complex-input transforms of the same size must produce correct results.
// Before the fix, the FFTW backend's plan cache keyed only on (nfft, inverse,
// inplace, aligned), so an r2c plan could be returned for a c2c call (or
// vice versa).
template <typename T>
void test_reuse_real_and_complex(int nfft) {
  typedef typename FFT<T>::Complex Complex;
  typedef Matrix<T, Dynamic, 1> ScalarVector;
  typedef Matrix<Complex, Dynamic, 1> ComplexVector;

  ScalarVector real_in(nfft);
  ComplexVector complex_in(nfft);
  for (int k = 0; k < nfft; ++k) {
    real_in[k] = (T)(rand() / (double)RAND_MAX - .5);
    complex_in[k] = Complex((T)(rand() / (double)RAND_MAX - .5), (T)(rand() / (double)RAND_MAX - .5));
  }

  FFT<T> fft;
  ComplexVector r2c_out;
  ComplexVector c2c_out;

  fft.fwd(r2c_out, real_in);
  fft.fwd(c2c_out, complex_in);
  VERIFY(T(fft_rmse(r2c_out, real_in)) < test_precision<T>());
  VERIFY(T(fft_rmse(c2c_out, complex_in)) < test_precision<T>());

  // Repeat with the reverse first-call ordering so the cache miss happens on
  // the opposite transform kind; this catches the symmetric c2c-then-r2c case.
  fft.fwd(c2c_out, complex_in);
  fft.fwd(r2c_out, real_in);
  VERIFY(T(fft_rmse(r2c_out, real_in)) < test_precision<T>());
  VERIFY(T(fft_rmse(c2c_out, complex_in)) < test_precision<T>());

  // Round-trip fwd->inv on the same shared FFT object exercises the c2r and
  // c2c inverse plans alongside the forward plans cached above.
  ScalarVector real_round;
  ComplexVector complex_round;
  fft.inv(real_round, r2c_out);
  VERIFY(T(dif_rmse(real_in, real_round)) < test_precision<T>());
  fft.inv(complex_round, c2c_out);
  VERIFY(T(dif_rmse(complex_in, complex_round)) < test_precision<T>());
}

EIGEN_DECLARE_TEST(FFTW) {
  CALL_SUBTEST(test_return_by_value(32));
  // Regression test for #1537 -- reuse one FFT object for both real and
  // complex inputs of the same size.
  CALL_SUBTEST(test_reuse_real_and_complex<float>(32));
  CALL_SUBTEST(test_reuse_real_and_complex<double>(32));
  CALL_SUBTEST(test_reuse_real_and_complex<float>(256));
  CALL_SUBTEST(test_reuse_real_and_complex<double>(256));
  CALL_SUBTEST(test_inplace_complex<float>(32));
  CALL_SUBTEST(test_inplace_complex<double>(32));
  CALL_SUBTEST(test_inplace_complex<float>(256));
  CALL_SUBTEST(test_inplace_complex<double>(256));
  CALL_SUBTEST(test_fwd_padding<float>(16));
  CALL_SUBTEST(test_fwd_padding<double>(16));
  CALL_SUBTEST(test_complex<float>(32));
  CALL_SUBTEST(test_complex<double>(32));
  CALL_SUBTEST(test_complex<float>(256));
  CALL_SUBTEST(test_complex<double>(256));
  CALL_SUBTEST(test_complex<float>(3 * 8));
  CALL_SUBTEST(test_complex<double>(3 * 8));
  CALL_SUBTEST(test_complex<float>(5 * 32));
  CALL_SUBTEST(test_complex<double>(5 * 32));
  CALL_SUBTEST(test_complex<float>(2 * 3 * 4));
  CALL_SUBTEST(test_complex<double>(2 * 3 * 4));
  CALL_SUBTEST(test_complex<float>(2 * 3 * 4 * 5));
  CALL_SUBTEST(test_complex<double>(2 * 3 * 4 * 5));
  CALL_SUBTEST(test_complex<float>(2 * 3 * 4 * 5 * 7));
  CALL_SUBTEST(test_complex<double>(2 * 3 * 4 * 5 * 7));

  CALL_SUBTEST(test_scalar<float>(32));
  CALL_SUBTEST(test_scalar<double>(32));
  CALL_SUBTEST(test_scalar<float>(45));
  CALL_SUBTEST(test_scalar<double>(45));
  CALL_SUBTEST(test_scalar<float>(50));
  CALL_SUBTEST(test_scalar<double>(50));
  CALL_SUBTEST(test_scalar<float>(256));
  CALL_SUBTEST(test_scalar<double>(256));
  CALL_SUBTEST(test_scalar<float>(2 * 3 * 4 * 5 * 7));
  CALL_SUBTEST(test_scalar<double>(2 * 3 * 4 * 5 * 7));

#if defined EIGEN_HAS_FFTWL || defined EIGEN_POCKETFFT_DEFAULT || defined EIGEN_DUCCFFT_DEFAULT
  CALL_SUBTEST(test_complex<long double>(32));
  CALL_SUBTEST(test_complex<long double>(256));
  CALL_SUBTEST(test_complex<long double>(3 * 8));
  CALL_SUBTEST(test_complex<long double>(5 * 32));
  CALL_SUBTEST(test_complex<long double>(2 * 3 * 4));
  CALL_SUBTEST(test_complex<long double>(2 * 3 * 4 * 5));
  CALL_SUBTEST(test_complex<long double>(2 * 3 * 4 * 5 * 7));

  CALL_SUBTEST(test_scalar<long double>(32));
  CALL_SUBTEST(test_scalar<long double>(45));
  CALL_SUBTEST(test_scalar<long double>(50));
  CALL_SUBTEST(test_scalar<long double>(256));
  CALL_SUBTEST(test_scalar<long double>(2 * 3 * 4 * 5 * 7));

  CALL_SUBTEST((test_complex2d<long double, 2 * 3 * 4, 2 * 3 * 4>()));
  CALL_SUBTEST((test_complex2d<long double, 3 * 4 * 5, 3 * 4 * 5>()));
  CALL_SUBTEST((test_complex2d<long double, 24, 60>()));
  CALL_SUBTEST((test_complex2d<long double, 60, 24>()));
// fail to build since Eigen limit the stack allocation size,too big here.
// CALL_SUBTEST( ( test_complex2d<long double, 256, 256> () ) );
#endif
#if defined EIGEN_FFTW_DEFAULT || defined EIGEN_POCKETFFT_DEFAULT || defined EIGEN_DUCCFFT_DEFAULT || \
    defined EIGEN_MKL_DEFAULT
  CALL_SUBTEST((test_complex2d<float, 24, 24>()));
  CALL_SUBTEST((test_complex2d<float, 60, 60>()));
  CALL_SUBTEST((test_complex2d<float, 24, 60>()));
  CALL_SUBTEST((test_complex2d<float, 60, 24>()));
#endif
#if defined EIGEN_FFTW_DEFAULT || defined EIGEN_POCKETFFT_DEFAULT || defined EIGEN_DUCCFFT_DEFAULT || \
    defined EIGEN_MKL_DEFAULT
  CALL_SUBTEST((test_complex2d<double, 24, 24>()));
  CALL_SUBTEST((test_complex2d<double, 60, 60>()));
  CALL_SUBTEST((test_complex2d<double, 24, 60>()));
  CALL_SUBTEST((test_complex2d<double, 60, 24>()));
#endif
}

#endif  // EIGEN_UNSUPPORTED_TEST_FFT_TEST_SHARED_H
