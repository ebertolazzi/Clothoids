// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Rasmus Munk Larsen <rmlarsen@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

// Generic CUDA runtime support shared across all GPU library integrations
// (cuSOLVER and cuBLAS):
//   - Error-checking macros
//   - RAII device buffer
//
// Only depends on <cuda_runtime.h>. No NVIDIA library headers.

#ifndef EIGEN_GPU_SUPPORT_H
#define EIGEN_GPU_SUPPORT_H

// IWYU pragma: private
#include "./InternalHeaderCheck.h"

#include <cuda_runtime.h>
#include <vector>

#include <limits>
#include <memory>

namespace Eigen {
namespace gpu {

// ---- Generic operation flag -------------------------------------------------
// Public flag for transpose/adjoint in BLAS-, solver-, and sparse-style calls.
// Each library's support header maps this to its own enum (cublasOperation_t,
// cusparseOperation_t, etc.) via a small to_<lib>_op() helper.

enum class GpuOp { NoTrans, Trans, ConjTrans };

namespace internal {

// ---- Error-checking macros --------------------------------------------------
// These abort (via eigen_assert) on failure. Not for use in destructors.

#define EIGEN_CUDA_RUNTIME_CHECK(expr)                             \
  do {                                                             \
    cudaError_t _e = (expr);                                       \
    eigen_assert(_e == cudaSuccess && "CUDA runtime call failed"); \
  } while (0)

// ---- Bounds-checked narrowing for cuBLAS/cuSOLVER int parameters ------------
// cuBLAS and the legacy cuSOLVER APIs take dimensions and leading dimensions as
// `int` (32-bit signed). Modern GPUs can host allocations whose dimensions
// exceed INT_MAX, and Eigen's Index is 64-bit by default. Use this helper at
// every narrowing call site so an out-of-range value triggers an assert
// instead of silently overflowing the BLAS argument.

inline int to_blas_int(int64_t v) {
  eigen_assert(v >= 0 && v <= static_cast<int64_t>((std::numeric_limits<int>::max)()) &&
               "dimension exceeds the int range supported by cuBLAS / cuSOLVER");
  return static_cast<int>(v);
}

// ---- Custom deleters for CUDA-allocated memory ------------------------------
// Used with std::unique_ptr to give CUDA allocations RAII semantics with no
// hand-rolled move/dtor boilerplate.

struct CudaFreeDeleter {
  // When `borrow == true`, the unique_ptr does not free the pointer. Used by
  // DeviceMatrix::view() to wrap a non-owning device pointer with the same
  // smart-pointer machinery as owning storage, without changing the type.
  bool borrow = false;
  void operator()(void* p) const noexcept {
    if (p && !borrow) (void)cudaFree(p);
  }
};

struct CudaFreeHostDeleter {
  void operator()(void* p) const noexcept {
    if (p) (void)cudaFreeHost(p);
  }
};

// ---- Thread-local pool of small device buffers ------------------------------
// Recycles allocations up to kSmallBufferThreshold bytes (e.g., DeviceScalar)
// to avoid cudaMalloc/cudaFree overhead. Larger allocations bypass the pool.

template <size_t SmallBufferThreshold = 256, size_t MaxPoolSize = 64>
struct DeviceBufferPool {
  static constexpr size_t kSmallBufferThreshold = SmallBufferThreshold;
  static constexpr size_t kMaxPoolSize = MaxPoolSize;

  struct Entry {
    void* ptr;
    size_t bytes;
  };

  ~DeviceBufferPool() {
    for (auto& e : free_list_) (void)cudaFree(e.ptr);
  }

  void* allocate(size_t bytes) {
    for (size_t i = 0; i < free_list_.size(); ++i) {
      if (free_list_[i].bytes >= bytes) {
        void* p = free_list_[i].ptr;
        free_list_[i] = free_list_.back();
        free_list_.pop_back();
        return p;
      }
    }
    void* p = nullptr;
    EIGEN_CUDA_RUNTIME_CHECK(cudaMalloc(&p, bytes));
    return p;
  }

  void deallocate(void* p, size_t bytes) {
    if (free_list_.size() < kMaxPoolSize) {
      free_list_.push_back({p, bytes});
    } else {
      (void)cudaFree(p);
    }
  }

  static DeviceBufferPool& threadLocal() {
    thread_local DeviceBufferPool pool;
    return pool;
  }

 private:
  std::vector<Entry> free_list_;
};

// Stateful deleter that returns small buffers to the thread-local pool and
// cudaFree's larger ones. size==0 means "always cudaFree" (for adopted ptrs).
struct PooledCudaFreeDeleter {
  size_t size = 0;

  void operator()(void* p) const noexcept {
    if (!p) return;
    if (size > 0 && size <= DeviceBufferPool<>::kSmallBufferThreshold) {
      DeviceBufferPool<>::threadLocal().deallocate(p, size);
    } else {
      (void)cudaFree(p);
    }
  }
};

// ---- RAII: device buffer ----------------------------------------------------

class DeviceBuffer {
 public:
  DeviceBuffer() = default;

  explicit DeviceBuffer(size_t bytes) {
    if (bytes > 0) {
      void* p = nullptr;
      if (bytes <= DeviceBufferPool<>::kSmallBufferThreshold) {
        p = DeviceBufferPool<>::threadLocal().allocate(bytes);
      } else {
        EIGEN_CUDA_RUNTIME_CHECK(cudaMalloc(&p, bytes));
      }
      ptr_ = std::unique_ptr<void, PooledCudaFreeDeleter>(p, PooledCudaFreeDeleter{bytes});
    }
  }

  void* get() const noexcept { return ptr_.get(); }
  void* release() noexcept { return ptr_.release(); }
  explicit operator bool() const noexcept { return static_cast<bool>(ptr_); }

  size_t size() const noexcept { return ptr_.get_deleter().size; }

  // Adopt an existing device pointer. Caller relinquishes ownership.
  // Adopted buffers bypass the pool on destruction (deleter size == 0).
  static DeviceBuffer adopt(void* p) noexcept {
    DeviceBuffer b;
    b.ptr_ = std::unique_ptr<void, PooledCudaFreeDeleter>(p, PooledCudaFreeDeleter{});
    return b;
  }

 private:
  std::unique_ptr<void, PooledCudaFreeDeleter> ptr_;
};

// ---- RAII: pinned host buffer -----------------------------------------------
// For async D2H copies (cudaMemcpyAsync requires pinned host memory for true
// asynchrony and to avoid compute-sanitizer warnings).

class PinnedHostBuffer {
 public:
  PinnedHostBuffer() = default;

  explicit PinnedHostBuffer(size_t bytes) {
    if (bytes > 0) {
      void* p = nullptr;
      EIGEN_CUDA_RUNTIME_CHECK(cudaMallocHost(&p, bytes));
      ptr_.reset(p);
    }
  }

  void* get() const noexcept { return ptr_.get(); }
  explicit operator bool() const noexcept { return static_cast<bool>(ptr_); }

 private:
  std::unique_ptr<void, CudaFreeHostDeleter> ptr_;
};

// ---- Scalar → cudaDataType_t ------------------------------------------------
// Shared by cuBLAS and cuSOLVER. cudaDataType_t is defined in library_types.h
// which is included transitively by cuda_runtime.h.

template <typename Scalar>
struct cuda_data_type;

template <>
struct cuda_data_type<float> {
  static constexpr cudaDataType_t value = CUDA_R_32F;
};
template <>
struct cuda_data_type<double> {
  static constexpr cudaDataType_t value = CUDA_R_64F;
};
template <>
struct cuda_data_type<std::complex<float>> {
  static constexpr cudaDataType_t value = CUDA_C_32F;
};
template <>
struct cuda_data_type<std::complex<double>> {
  static constexpr cudaDataType_t value = CUDA_C_64F;
};

}  // namespace internal
}  // namespace gpu
}  // namespace Eigen

#endif  // EIGEN_GPU_SUPPORT_H
