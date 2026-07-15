// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Rasmus Munk Larsen <rmlarsen@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

// Unified GPU execution context.
//
// gpu::Context owns a CUDA stream and NVIDIA library handles (cuBLAS
// eagerly, cuSOLVER / cuBLASLt / cuSPARSE lazily on first use). It is the
// entry point for all GPU operations on gpu::DeviceMatrix.
//
// The cuSOLVER handle is created on the first call to cusolverHandle()
// so that translation units which only use cuFFT or cuBLAS paths (e.g.
// the cufft test) do not pull cusolverDn* symbols into the link.
//
// Usage:
//   gpu::Context ctx;                        // explicit context
//   d_C.device(ctx) = d_A * d_B;            // GEMM on ctx's stream
//
//   d_C = d_A * d_B;                        // thread-local default context
//   gpu::Context& ctx = gpu::Context::threadLocal();

#ifndef EIGEN_GPU_CONTEXT_H
#define EIGEN_GPU_CONTEXT_H

// IWYU pragma: private
#include "./InternalHeaderCheck.h"

#include "./CuBlasSupport.h"
#include "./CuSolverSupport.h"
#include <cusparse.h>

namespace Eigen {
namespace gpu {

/** \ingroup GPU_Module
 * \class Context
 * \brief Unified GPU execution context owning a CUDA stream and library handles.
 *
 * Each Context instance creates a dedicated CUDA stream and a cuBLAS handle
 * bound to that stream. The cuSOLVER handle is created on first use via
 * cusolverHandle(); translation units that never call it do not require
 * cuSOLVER at link time. cuBLASLt and cuSPARSE handles are similarly lazy.
 * Multiple contexts enable concurrent execution on independent streams.
 *
 * A lazily-created thread-local default is available via threadLocal() for
 * simple single-stream usage. A single Context is not thread-safe — use one
 * per thread, or external synchronization, since cuBLAS / cuSOLVER handles
 * are not thread-safe per handle and lazy-init of secondary handles is racy.
 */
class Context {
 public:
  /** Create a new context with a dedicated CUDA stream. */
  Context() {
    EIGEN_CUDA_RUNTIME_CHECK(cudaStreamCreate(&stream_));
    init_cublas();
  }

  /** Create a context on an existing stream (e.g., stream 0 = nullptr).
   * The caller retains ownership of the stream — this context will not destroy it. */
  explicit Context(cudaStream_t stream) : stream_(stream), owns_stream_(false) { init_cublas(); }

  ~Context() {
    // Indirect calls keep cusolverDnDestroy / cusparseDestroy out of TUs that
    // never call cusolverHandle() / cusparseHandle() (e.g. the cufft test).
    if (cusparse_destroyer_) (void)cusparse_destroyer_(cusparse_);
    if (cusolver_destroyer_) (void)cusolver_destroyer_(cusolver_);
    // Release plan-cache descriptors before tearing down the cuBLASLt handle.
    gemm_plan_cache_.clear();
    if (cublas_lt_) (void)cublasLtDestroy(cublas_lt_);
    if (cublas_) (void)cublasDestroy(cublas_);
    if (owns_stream_ && stream_) (void)cudaStreamDestroy(stream_);
  }

  // Non-copyable, non-movable (owns library handles).
  Context(const Context&) = delete;
  Context& operator=(const Context&) = delete;
  Context(Context&&) = delete;
  Context& operator=(Context&&) = delete;

  /** Get the thread-local default context.
   *
   * If setThreadLocal() has been called, returns that context.
   * Otherwise lazily creates a new context with a dedicated stream.
   *
   * \note The thread-local instance is destroyed when the thread exits (or at
   * static destruction time for the main thread). On some CUDA driver
   * configurations this may print "CUDA_ERROR_DEINITIALIZED" to stderr if the
   * CUDA context has already been torn down. These errors are harmless and are
   * suppressed in the destructor, but they can produce noise in test output.
   * To avoid this, call cudaDeviceReset() only after all Context instances
   * (including thread-local ones) have been destroyed. */
  static Context& threadLocal() {
    Context* override = tl_override_ptr();
    if (override) return *override;
    thread_local Context ctx;
    return ctx;
  }

  /** Override the thread-local default context for this thread.
   * The caller retains ownership of \p ctx — it must outlive all uses.
   * Pass nullptr to restore the lazily-created default. */
  static void setThreadLocal(Context* ctx) { tl_override_ptr() = ctx; }

  cudaStream_t stream() const { return stream_; }
  cublasHandle_t cublasHandle() const { return cublas_; }

  /** Returns the cuSOLVER handle, creating it on first call. */
  cusolverDnHandle_t cusolverHandle() {
    if (!cusolver_) {
      EIGEN_CUSOLVER_CHECK(cusolverDnCreate(&cusolver_));
      EIGEN_CUSOLVER_CHECK(cusolverDnSetStream(cusolver_, stream_));
      cusolver_destroyer_ = &destroyCusolver;
    }
    return cusolver_;
  }

  /** cuBLASLt handle (lazy-initialized on first GEMM call). */
  cublasLtHandle_t cublasLtHandle() const {
    if (!cublas_lt_) {
      EIGEN_CUBLAS_CHECK(cublasLtCreate(&cublas_lt_));
    }
    return cublas_lt_;
  }

  /** Workspace buffer for cublasLtMatmul (grown lazily by cublaslt_gemm).
   * Not thread-safe — all GEMM calls must be on this context's stream. */
  internal::DeviceBuffer* gemmWorkspace() const { return &gemm_workspace_; }

  /** Plan cache for cublasLtMatmul (caches descriptors and selected algorithm
   * by shape to avoid per-call overhead). Same thread-safety as workspace. */
  internal::CublasLtPlanCache* gemmPlanCache() const { return &gemm_plan_cache_; }

  /** Workspace ceiling passed to the cublasLtMatmul heuristic at plan-creation time.
   * Defaults to internal::kCublasLtMaxWorkspaceBytes (compile-time configurable via
   * EIGEN_CUDA_CUBLASLT_MAX_WORKSPACE_BYTES). */
  std::size_t cublasLtMaxWorkspaceBytes() const { return cublaslt_max_workspace_bytes_; }

  /** Override the workspace ceiling for future plan-cache misses on this context.
   * The cap is consulted at plan-creation time only; pre-existing cached plans
   * keep the cap they were built with. Call gemmPlanCache()->clear() to force
   * re-selection under the new cap. */
  void setCublasLtMaxWorkspaceBytes(std::size_t bytes) { cublaslt_max_workspace_bytes_ = bytes; }

  /** cuSPARSE handle (lazy-initialized on first call). */
  cusparseHandle_t cusparseHandle() const {
    if (!cusparse_) {
      cusparseStatus_t s1 = cusparseCreate(&cusparse_);
      eigen_assert(s1 == CUSPARSE_STATUS_SUCCESS && "cusparseCreate failed");
      EIGEN_UNUSED_VARIABLE(s1);
      cusparseStatus_t s2 = cusparseSetStream(cusparse_, stream_);
      eigen_assert(s2 == CUSPARSE_STATUS_SUCCESS && "cusparseSetStream failed");
      EIGEN_UNUSED_VARIABLE(s2);
      cusparse_destroyer_ = &destroyCusparse;
    }
    return cusparse_;
  }

 private:
  static cusolverStatus_t destroyCusolver(cusolverDnHandle_t h) { return cusolverDnDestroy(h); }
  static cusparseStatus_t destroyCusparse(cusparseHandle_t h) { return cusparseDestroy(h); }

  cudaStream_t stream_ = nullptr;
  cublasHandle_t cublas_ = nullptr;
  cusolverDnHandle_t cusolver_ = nullptr;
  cusolverStatus_t (*cusolver_destroyer_)(cusolverDnHandle_t) = nullptr;
  mutable cublasLtHandle_t cublas_lt_ = nullptr;  // lazy
  mutable cusparseHandle_t cusparse_ = nullptr;   // lazy
  mutable cusparseStatus_t (*cusparse_destroyer_)(cusparseHandle_t) = nullptr;
  mutable internal::DeviceBuffer gemm_workspace_;  // lazy
  mutable internal::CublasLtPlanCache gemm_plan_cache_{internal::kCublasLtPlanCacheCapacity};
  std::size_t cublaslt_max_workspace_bytes_ = internal::kCublasLtMaxWorkspaceBytes;
  bool owns_stream_ = true;

  static Context*& tl_override_ptr() {
    thread_local Context* ptr = nullptr;
    return ptr;
  }

  void init_cublas() {
    EIGEN_CUBLAS_CHECK(cublasCreate(&cublas_));
    EIGEN_CUBLAS_CHECK(cublasSetStream(cublas_, stream_));
  }
};

}  // namespace gpu
}  // namespace Eigen

#endif  // EIGEN_GPU_CONTEXT_H
