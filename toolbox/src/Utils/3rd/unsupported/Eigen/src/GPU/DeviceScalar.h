// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Rasmus Munk Larsen <rmlarsen@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

// Device-resident scalar for deferred host synchronization.
//
// DeviceScalar<Scalar> wraps a single value in device memory. Reductions
// (dot, nrm2) write results directly to device memory via
// CUBLAS_POINTER_MODE_DEVICE, deferring host sync until the value is read.
//
// Implicit conversion to Scalar triggers cudaStreamSynchronize + download.
// In CG, this reduces 3 syncs/iter to effectively 1: the first conversion
// syncs the stream, subsequent conversions in the same expression just
// download (the stream is already flushed).
//
// Usage:
//   auto dot_val = d_x.dot(d_y);       // DeviceScalar, no sync
//   auto norm_val = d_r.squaredNorm();  // DeviceScalar, no sync
//   Scalar alpha = absNew / dot_val;    // sync here (both values downloaded)
//   d_x += alpha * d_p;                 // host-scalar axpy (as before)

#ifndef EIGEN_GPU_DEVICE_SCALAR_H
#define EIGEN_GPU_DEVICE_SCALAR_H

// IWYU pragma: private
#include "./InternalHeaderCheck.h"

#include "./GpuSupport.h"
#include "./DeviceScalarOps.h"

namespace Eigen {
namespace gpu {

template <typename Scalar_>
class DeviceScalar {
 public:
  using Scalar = Scalar_;

  /** Allocate uninitialized device scalar. Contents are undefined until written
   * (e.g., by cuBLAS dot/nrm2 with POINTER_MODE_DEVICE). Consistent with
   * DeviceMatrix(rows, cols) which also does not zero-initialize. */
  explicit DeviceScalar(cudaStream_t stream = nullptr) : d_val_(sizeof(Scalar)), stream_(stream) {}

  DeviceScalar(Scalar host_val, cudaStream_t stream) : d_val_(sizeof(Scalar)), stream_(stream) {
    EIGEN_CUDA_RUNTIME_CHECK(cudaMemcpyAsync(d_val_.get(), &host_val, sizeof(Scalar), cudaMemcpyHostToDevice, stream_));
  }

  DeviceScalar(DeviceScalar&& o) noexcept : d_val_(std::move(o.d_val_)), stream_(o.stream_) { o.stream_ = nullptr; }

  DeviceScalar& operator=(DeviceScalar&& o) noexcept {
    if (this != &o) {
      d_val_ = std::move(o.d_val_);
      stream_ = o.stream_;
      o.stream_ = nullptr;
    }
    return *this;
  }

  DeviceScalar(const DeviceScalar&) = delete;
  DeviceScalar& operator=(const DeviceScalar&) = delete;

  /** Download from device. Synchronizes the stream on first call;
   * subsequent calls in the same expression are cheap (stream already flushed). */
  Scalar get() const {
    Scalar result;
    EIGEN_CUDA_RUNTIME_CHECK(cudaMemcpyAsync(&result, d_val_.get(), sizeof(Scalar), cudaMemcpyDeviceToHost, stream_));
    EIGEN_CUDA_RUNTIME_CHECK(cudaStreamSynchronize(stream_));
    return result;
  }

  /** Implicit conversion — allows `Scalar alpha = deviceScalar` and
   * `if (deviceScalar < threshold)`. Triggers sync. */
  operator Scalar() const { return get(); }

  Scalar* devicePtr() { return static_cast<Scalar*>(d_val_.get()); }
  const Scalar* devicePtr() const { return static_cast<const Scalar*>(d_val_.get()); }
  cudaStream_t stream() const { return stream_; }

  // ---- Device-side arithmetic (no host sync) ---------------------------------
  // Uses NPP from DeviceScalarOps.h. All results stay on device.
  // Currently supports real types only (float, double). Complex types
  // fall back to implicit conversion (host sync) for division.
  //
  // Note: DeviceScalar has no cross-stream readiness tracking. All
  // operations must be on the same CUDA stream. This is the natural
  // pattern in iterative solvers where one GpuContext owns all work.

  friend DeviceScalar operator/(const DeviceScalar& a, const DeviceScalar& b) {
    eigen_assert(a.stream_ == b.stream_ && "DeviceScalar operator/: operands must share the same stream");
    DeviceScalar result(a.stream_);
    gpu::internal::device_scalar_div(a.devicePtr(), b.devicePtr(), result.devicePtr(), a.stream_);
    return result;
  }

  friend DeviceScalar operator/(Scalar a, const DeviceScalar& b) {
    DeviceScalar d_a(a, b.stream_);
    return d_a / b;
  }

  friend DeviceScalar operator/(const DeviceScalar& a, Scalar b) {
    DeviceScalar d_b(b, a.stream_);
    return a / d_b;
  }

  DeviceScalar operator-() const {
    DeviceScalar result(stream_);
    gpu::internal::device_scalar_neg(devicePtr(), result.devicePtr(), stream_);
    return result;
  }

 private:
  internal::DeviceBuffer d_val_;
  cudaStream_t stream_ = nullptr;
};

}  // namespace gpu
}  // namespace Eigen

#endif  // EIGEN_GPU_DEVICE_SCALAR_H
