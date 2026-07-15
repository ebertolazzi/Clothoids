// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Rasmus Munk Larsen <rmlarsen@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

// GPU sparse matrix-vector multiply (SpMV) and sparse matrix-dense matrix
// multiply (SpMM) via cuSPARSE.
//
// SparseContext manages cuSPARSE descriptors and device buffers. It accepts
// Eigen SparseMatrix<Scalar, ColMajor> (CSC) and performs SpMV/SpMM on the GPU.
// RowMajor input is implicitly converted to ColMajor.
//
// Can borrow a Context for same-stream execution with BLAS-1 ops (zero
// event overhead in iterative solvers like CG).
//
// Thread safety: not thread-safe. Concurrent multiply* calls on a single
// SparseContext race on the cuSPARSE handle, the bound stream, and the
// cached device buffers. Use one SparseContext per thread.
//
// Usage:
//   // Standalone (own stream):
//   gpu::SparseContext<double> ctx;
//   VectorXd y = ctx.multiply(A, x);                  // y = A * x
//   ctx.multiply(A, x, y, 2.0, 1.0);                  // y = 2*A*x + y
//   ctx.multiply(A, x, y, 1.0, 0.0, gpu::GpuOp::ConjTrans);  // y = A^H * x
//   VectorXd z = ctx.multiplyT(A, x);                 // z = A^T * x
//   VectorXcd w = ctx.multiplyAdjoint(A, x);          // w = A^H * x (complex)
//   MatrixXd Y = ctx.multiplyMat(A, X);               // Y = A * X (multiple RHS)
//
//   // Shared context (same stream as BLAS-1 ops):
//   gpu::Context gpu_ctx;
//   gpu::SparseContext<double> sparse_ctx(gpu_ctx);
//   VectorXd y = sparse_ctx.multiply(A, x);
//
//   // Device-resident (no host roundtrip):
//   sparse_ctx.multiply(A, d_x, d_y);                 // DeviceMatrix in/out

#ifndef EIGEN_GPU_SPARSE_CONTEXT_H
#define EIGEN_GPU_SPARSE_CONTEXT_H

// IWYU pragma: private
#include "./InternalHeaderCheck.h"

#include "./CuSparseSupport.h"

namespace Eigen {
namespace gpu {

// Forward declarations.
template <typename Scalar_>
class SparseContext;
template <typename Scalar_>
class DeviceSparseView;

/** SpMV expression: DeviceSparseView * DeviceMatrix → SpMVExpr.
 * Evaluated by DeviceMatrix::operator=(SpMVExpr). */
template <typename Scalar_>
class SpMVExpr {
 public:
  using Scalar = Scalar_;
  SpMVExpr(const DeviceSparseView<Scalar>& view, const DeviceMatrix<Scalar>& x) : view_(view), x_(x) {}
  const DeviceSparseView<Scalar>& view() const { return view_; }
  const DeviceMatrix<Scalar>& x() const { return x_; }

 private:
  const DeviceSparseView<Scalar>& view_;
  const DeviceMatrix<Scalar>& x_;
};

/** Device-resident sparse matrix view. Returned by SparseContext::deviceView().
 * Lightweight handle referencing the context's cached device data.
 *
 * \warning One SparseContext caches one sparse matrix at a time.
 * Creating a second deviceView on the same context overwrites the first.
 * For multiple simultaneous sparse matrices, use separate SparseContext
 * instances (they can share a Context for same-stream execution).
 *
 * Supports `d_y = d_A * d_x` via SpMVExpr. */
template <typename Scalar_>
class DeviceSparseView {
 public:
  using Scalar = Scalar_;
  using SpMat = SparseMatrix<Scalar, ColMajor, int>;

  DeviceSparseView(SparseContext<Scalar>& ctx, Index rows, Index cols) : ctx_(ctx), rows_(rows), cols_(cols) {}

  /** SpMV expression: d_A * d_x. Evaluated by DeviceMatrix::operator=. */
  SpMVExpr<Scalar> operator*(const DeviceMatrix<Scalar>& x) const { return SpMVExpr<Scalar>(*this, x); }

  Index rows() const { return rows_; }
  Index cols() const { return cols_; }
  const SparseContext<Scalar>& context() const { return ctx_; }

 private:
  SparseContext<Scalar>& ctx_;
  Index rows_;
  Index cols_;
};

template <typename Scalar_>
class SparseContext {
 public:
  using Scalar = Scalar_;
  using RealScalar = typename NumTraits<Scalar>::Real;
  using StorageIndex = int;
  using SpMat = SparseMatrix<Scalar, ColMajor, StorageIndex>;
  using DenseVector = Matrix<Scalar, Dynamic, 1>;
  using DenseMatrix = Matrix<Scalar, Dynamic, Dynamic, ColMajor>;

  /** Standalone: creates own stream and cuSPARSE handle. */
  SparseContext() : owns_handle_(true) {
    EIGEN_CUDA_RUNTIME_CHECK(cudaStreamCreate(&stream_));
    owns_stream_ = true;
    EIGEN_CUSPARSE_CHECK(cusparseCreate(&handle_));
    EIGEN_CUSPARSE_CHECK(cusparseSetStream(handle_, stream_));
  }

  /** Borrow a Context: shares stream and cuSPARSE handle.
   * The Context must outlive this SparseContext. */
  explicit SparseContext(Context& ctx)
      : stream_(ctx.stream()), handle_(ctx.cusparseHandle()), owns_stream_(false), owns_handle_(false) {}

  ~SparseContext() {
    destroy_descriptors_unchecked();
    if (owns_handle_ && handle_) (void)cusparseDestroy(handle_);
    if (owns_stream_ && stream_) (void)cudaStreamDestroy(stream_);
  }

  SparseContext(const SparseContext&) = delete;
  SparseContext& operator=(const SparseContext&) = delete;

  // ---- Device sparse view (for expression syntax: d_y = d_A * d_x) ----------

  /** Upload a sparse matrix to device and return a lightweight view.
   * The sparse data is uploaded immediately and cached in this context.
   * The returned view can be used for repeated SpMV without re-uploading.
   * If the matrix values change, call deviceView() again to re-upload.
   *
   * \warning One context caches one matrix. Calling deviceView() again
   * overwrites the previous upload. For multiple simultaneous matrices,
   * use separate SparseContext instances sharing the same Context.
   *
   * Supports `d_y = d_A * d_x` expression syntax. */
  DeviceSparseView<Scalar> deviceView(const SpMat& A) {
    eigen_assert(A.isCompressed());
    upload_sparse(A);
    return DeviceSparseView<Scalar>(*this, A.rows(), A.cols());
  }

  // ---- SpMV: y = A * x (host vectors) --------------------------------------

  /** Compute y = A * x. Returns y as a new dense vector. */
  template <typename InputType, typename Rhs>
  DenseVector multiply(const SparseMatrixBase<InputType>& A, const MatrixBase<Rhs>& x) {
    const InputType& input = A.derived();
    check_storage_index_bounds(input.rows(), input.cols(), input.nonZeros());
    const SpMat mat(input);
    DenseVector y(mat.rows());
    y.setZero();
    multiply_host_impl(mat, x.derived(), y, Scalar(1), Scalar(0), CUSPARSE_OPERATION_NON_TRANSPOSE);
    return y;
  }

  /** Compute y = alpha * op(A) * x + beta * y (in-place, host vectors). */
  template <typename InputType, typename Rhs, typename Dest>
  void multiply(const SparseMatrixBase<InputType>& A, const MatrixBase<Rhs>& x, MatrixBase<Dest>& y,
                Scalar alpha = Scalar(1), Scalar beta = Scalar(0), GpuOp op = GpuOp::NoTrans) {
    const InputType& input = A.derived();
    check_storage_index_bounds(input.rows(), input.cols(), input.nonZeros());
    const SpMat mat(input);
    multiply_host_impl(mat, x.derived(), y.derived(), alpha, beta, internal::to_cusparse_op_for_scalar<Scalar>(op));
  }

  // ---- SpMV: y = A * x (DeviceMatrix, no host roundtrip) -------------------

  /** Compute d_y = A * d_x. Device-resident, no host transfer.
   * Sparse matrix A is uploaded to device (cached). Dense vectors stay on device. */
  template <typename InputType>
  void multiply(const SparseMatrixBase<InputType>& A, const DeviceMatrix<Scalar>& d_x, DeviceMatrix<Scalar>& d_y) {
    const SpMat mat(A.derived());
    multiply_device_impl(mat, d_x, d_y, Scalar(1), Scalar(0), CUSPARSE_OPERATION_NON_TRANSPOSE);
  }

  /** Compute d_y = alpha * op(A) * d_x + beta * d_y (DeviceMatrix, in-place). */
  template <typename InputType>
  void multiply(const SparseMatrixBase<InputType>& A, const DeviceMatrix<Scalar>& d_x, DeviceMatrix<Scalar>& d_y,
                Scalar alpha, Scalar beta, cusparseOperation_t op = CUSPARSE_OPERATION_NON_TRANSPOSE) {
    const SpMat mat(A.derived());
    multiply_device_impl(mat, d_x, d_y, alpha, beta, op);
  }

  // ---- SpMV transpose -------------------------------------------------------

  /** Compute y = A^T * x (host vectors). */
  template <typename InputType, typename Rhs>
  DenseVector multiplyT(const SparseMatrixBase<InputType>& A, const MatrixBase<Rhs>& x) {
    const InputType& input = A.derived();
    check_storage_index_bounds(input.rows(), input.cols(), input.nonZeros());
    const SpMat mat(input);
    DenseVector y(mat.cols());
    y.setZero();
    multiply_host_impl(mat, x.derived(), y, Scalar(1), Scalar(0), CUSPARSE_OPERATION_TRANSPOSE);
    return y;
  }

  // ---- SpMV adjoint: y = A^H * x -------------------------------------------

  /** Compute y = A^H * x (conjugate transpose). For real Scalar this is equivalent to multiplyT. */
  template <typename InputType, typename Rhs>
  DenseVector multiplyAdjoint(const SparseMatrixBase<InputType>& A, const MatrixBase<Rhs>& x) {
    const InputType& input = A.derived();
    check_storage_index_bounds(input.rows(), input.cols(), input.nonZeros());
    const SpMat mat(input);
    DenseVector y(mat.cols());
    y.setZero();
    multiply_host_impl(mat, x.derived(), y, Scalar(1), Scalar(0),
                       internal::to_cusparse_op_for_scalar<Scalar>(GpuOp::ConjTrans));
    return y;
  }

  // ---- SpMM: Y = op(A) * X (multiple RHS) ----------------------------------

  /** Compute Y = op(A) * X where X is a dense matrix (multiple RHS). Returns Y. */
  template <typename InputType, typename Rhs>
  DenseMatrix multiplyMat(const SparseMatrixBase<InputType>& A, const MatrixBase<Rhs>& X, GpuOp op = GpuOp::NoTrans) {
    const InputType& input = A.derived();
    check_storage_index_bounds(input.rows(), input.cols(), input.nonZeros());
    const SpMat mat(input);
    const DenseMatrix rhs(X.derived());

    const cusparseOperation_t cu_op = internal::to_cusparse_op_for_scalar<Scalar>(op);
    const Index m = (op == GpuOp::NoTrans) ? mat.rows() : mat.cols();
    const Index k = (op == GpuOp::NoTrans) ? mat.cols() : mat.rows();
    eigen_assert(k == rhs.rows());

    const Index n = rhs.cols();
    if (m == 0 || n == 0 || mat.nonZeros() == 0) return DenseMatrix::Zero(m, n);

    DenseMatrix Y = DenseMatrix::Zero(m, n);
    spmm_impl(mat, rhs, Y, Scalar(1), Scalar(0), cu_op);
    return Y;
  }

  // ---- Accessors ------------------------------------------------------------

  cudaStream_t stream() const { return stream_; }

 private:
  cudaStream_t stream_ = nullptr;
  cusparseHandle_t handle_ = nullptr;
  bool owns_stream_ = false;
  bool owns_handle_ = false;

  // Cached device buffers for sparse matrix (grow-only).
  internal::DeviceBuffer d_outerPtr_;
  internal::DeviceBuffer d_innerIdx_;
  internal::DeviceBuffer d_values_;
  size_t d_outerPtr_size_ = 0;
  size_t d_innerIdx_size_ = 0;
  size_t d_values_size_ = 0;

  // Cached device buffers for host-API dense vectors (grow-only).
  internal::DeviceBuffer d_x_;
  internal::DeviceBuffer d_y_;
  size_t d_x_size_ = 0;
  size_t d_y_size_ = 0;

  mutable internal::DeviceBuffer d_workspace_;
  mutable size_t d_workspace_size_ = 0;

  // Cached cuSPARSE sparse matrix descriptor.
  cusparseSpMatDescr_t spmat_desc_ = nullptr;
  Index cached_rows_ = -1;
  Index cached_cols_ = -1;
  Index cached_nnz_ = -1;

  // ---- SpMV with host vectors (upload/download per call) --------------------

  template <typename RhsDerived, typename DestDerived>
  void multiply_host_impl(const SpMat& A, const RhsDerived& x, DestDerived& y, Scalar alpha, Scalar beta,
                          cusparseOperation_t op) {
    eigen_assert(A.isCompressed());

    const Index m = A.rows();
    const Index n = A.cols();
    const Index nnz = A.nonZeros();
    const Index x_size = (op == CUSPARSE_OPERATION_NON_TRANSPOSE) ? n : m;
    const Index y_size = (op == CUSPARSE_OPERATION_NON_TRANSPOSE) ? m : n;

    eigen_assert(x.size() == x_size);
    eigen_assert(y.size() == y_size);

    if (m == 0 || n == 0 || nnz == 0) {
      if (beta == Scalar(0))
        y.setZero();
      else
        y *= beta;
      return;
    }

    upload_sparse(A);

    ensure_buffer(d_x_, d_x_size_, static_cast<size_t>(x_size) * sizeof(Scalar));
    const DenseVector x_tmp(x);
    EIGEN_CUDA_RUNTIME_CHECK(
        cudaMemcpyAsync(d_x_.get(), x_tmp.data(), x_size * sizeof(Scalar), cudaMemcpyHostToDevice, stream_));

    ensure_buffer(d_y_, d_y_size_, static_cast<size_t>(y_size) * sizeof(Scalar));
    if (beta != Scalar(0)) {
      const DenseVector y_tmp(y);
      EIGEN_CUDA_RUNTIME_CHECK(
          cudaMemcpyAsync(d_y_.get(), y_tmp.data(), y_size * sizeof(Scalar), cudaMemcpyHostToDevice, stream_));
    }

    exec_spmv(x_size, y_size, d_x_.get(), d_y_.get(), alpha, beta, op);

    EIGEN_CUDA_RUNTIME_CHECK(
        cudaMemcpyAsync(y.data(), d_y_.get(), y_size * sizeof(Scalar), cudaMemcpyDeviceToHost, stream_));
    EIGEN_CUDA_RUNTIME_CHECK(cudaStreamSynchronize(stream_));
  }

  // ---- SpMV with DeviceMatrix (no host transfer) ----------------------------

  // Called by public multiply(A, d_x, d_y) — always re-uploads A.
  void multiply_device_impl(const SpMat& A, const DeviceMatrix<Scalar>& d_x, DeviceMatrix<Scalar>& d_y, Scalar alpha,
                            Scalar beta, cusparseOperation_t op) {
    upload_sparse(A);
    spmv_device_exec(d_x, d_y, alpha, beta, op);
  }

 public:
  /** Execute SpMV using the already-uploaded sparse matrix (no re-upload).
   * Used by SpMVExpr (d_y = d_A * d_x) for cached deviceView() paths.
   * The sparse matrix must have been uploaded via deviceView() or multiply(). */
  void spmv_device_exec(const DeviceMatrix<Scalar>& d_x, DeviceMatrix<Scalar>& d_y, Scalar alpha = Scalar(1),
                        Scalar beta = Scalar(0), cusparseOperation_t op = CUSPARSE_OPERATION_NON_TRANSPOSE) const {
    eigen_assert(spmat_desc_ && "sparse matrix not uploaded — call deviceView() or multiply() first");
    // cuSPARSE SpMV: y must not alias x (undefined behavior).
    eigen_assert(d_x.data() != d_y.data() && "SpMV: output aliases input vector");

    const Index m = cached_rows_;
    const Index n = cached_cols_;
    const Index x_size = (op == CUSPARSE_OPERATION_NON_TRANSPOSE) ? n : m;
    const Index y_size = (op == CUSPARSE_OPERATION_NON_TRANSPOSE) ? m : n;

    eigen_assert(d_x.rows() * d_x.cols() == x_size);

    if (m == 0 || n == 0 || cached_nnz_ == 0) {
      // Empty A reduces SpMV to y <- beta*y; SparseContext owns no cuBLAS
      // handle for the scale, so the beta != 0 case must be handled by the caller.
      eigen_assert(beta == Scalar(0) && "SpMV with empty A and beta != 0 is unsupported; scale d_y externally");
      if (d_y.rows() * d_y.cols() != y_size) d_y.resize(y_size, 1);
      d_y.setZero(stream_);
      return;
    }

    // Ensure d_y is allocated.
    if (d_y.rows() * d_y.cols() != y_size) {
      d_y.resize(y_size, 1);
    }

    // Wait for input data to be ready on this stream.
    d_x.waitReady(stream_);
    d_y.waitReady(stream_);

    exec_spmv(x_size, y_size, const_cast<void*>(static_cast<const void*>(d_x.data())), static_cast<void*>(d_y.data()),
              alpha, beta, op);

    d_y.recordReady(stream_);
  }

 private:
  // cuSPARSE 11.x's cusparseSpMM rejects CSC for matA (CSC support landed in
  // CUDA 12.0). On 11.x we register the same buffers as CSR-of-A^T (dims
  // swapped) and invert the user-facing op before each cuSPARSE call. On 12+
  // we keep the natural CSC path so users pay no extra cost.
#if !defined(CUSPARSE_VERSION) || CUSPARSE_VERSION < 12000
  static constexpr bool kUseCsrOfTranspose = true;
  static constexpr cusparseSpMMAlg_t kSpMMAlg = CUSPARSE_SPMM_CSR_ALG2;
#else
  static constexpr bool kUseCsrOfTranspose = false;
  static constexpr cusparseSpMMAlg_t kSpMMAlg = CUSPARSE_SPMM_ALG_DEFAULT;
#endif

  // Map a user-facing op on A to the cuSPARSE op on the cached descriptor.
  // Identity on cuSPARSE 12+ (descriptor is CSC of A); inverted on 11.x
  // (descriptor is CSR of A^T).
  static cusparseOperation_t descriptor_op(cusparseOperation_t user_op) {
    if (!kUseCsrOfTranspose) return user_op;
    switch (user_op) {
      case CUSPARSE_OPERATION_NON_TRANSPOSE:
        return CUSPARSE_OPERATION_TRANSPOSE;
      case CUSPARSE_OPERATION_TRANSPOSE:
        return CUSPARSE_OPERATION_NON_TRANSPOSE;
      default:
        // CONJUGATE_TRANSPOSE on the CSR-of-A^T descriptor would compute
        // conj(A) * x, not A^H * x — not supported via this representation.
        eigen_assert(false && "CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE not supported on cuSPARSE < 12.0");
        return user_op;
    }
  }

  // ---- Shared SpMV execution ------------------------------------------------

  void exec_spmv(Index x_size, Index y_size, void* d_x_ptr, void* d_y_ptr, Scalar alpha, Scalar beta,
                 cusparseOperation_t op) const {
    constexpr cudaDataType_t dtype = internal::cuda_data_type<Scalar>::value;
    const cusparseOperation_t cu_op = descriptor_op(op);
    cusparseDnVecDescr_t x_desc = nullptr, y_desc = nullptr;
    EIGEN_CUSPARSE_CHECK(cusparseCreateDnVec(&x_desc, x_size, d_x_ptr, dtype));
    EIGEN_CUSPARSE_CHECK(cusparseCreateDnVec(&y_desc, y_size, d_y_ptr, dtype));

    size_t ws_size = 0;
    EIGEN_CUSPARSE_CHECK(cusparseSpMV_bufferSize(handle_, cu_op, &alpha, spmat_desc_, x_desc, &beta, y_desc, dtype,
                                                 CUSPARSE_SPMV_ALG_DEFAULT, &ws_size));
    ensure_buffer(d_workspace_, d_workspace_size_, ws_size);

    EIGEN_CUSPARSE_CHECK(cusparseSpMV(handle_, cu_op, &alpha, spmat_desc_, x_desc, &beta, y_desc, dtype,
                                      CUSPARSE_SPMV_ALG_DEFAULT, d_workspace_.get()));

    EIGEN_CUSPARSE_CHECK(cusparseDestroyDnVec(x_desc));
    EIGEN_CUSPARSE_CHECK(cusparseDestroyDnVec(y_desc));
  }

  // ---- SpMM implementation --------------------------------------------------

  void spmm_impl(const SpMat& A, const DenseMatrix& X, DenseMatrix& Y, Scalar alpha, Scalar beta,
                 cusparseOperation_t op) {
    eigen_assert(A.isCompressed());

    // For op != NON_TRANSPOSE, Y = op(A) * X. The dense X / Y descriptors must
    // describe the *post-op* shapes: X has k_op rows (= input dim of op(A)),
    // Y has m_op rows (= output dim of op(A)).
    const bool transposed = (op != CUSPARSE_OPERATION_NON_TRANSPOSE);
    const Index m_op = transposed ? A.cols() : A.rows();
    const Index k_op = transposed ? A.rows() : A.cols();
    const Index n = X.cols();
    const Index nnz = A.nonZeros();

    if (m_op == 0 || n == 0 || k_op == 0 || nnz == 0) {
      if (beta == Scalar(0))
        Y.setZero();
      else
        Y *= beta;
      return;
    }

    upload_sparse(A);

    // Upload X to device. X is k_op x n, Y is m_op x n (column-major).
    const size_t x_bytes = static_cast<size_t>(k_op) * static_cast<size_t>(n) * sizeof(Scalar);
    const size_t y_bytes = static_cast<size_t>(m_op) * static_cast<size_t>(n) * sizeof(Scalar);
    ensure_buffer(d_x_, d_x_size_, x_bytes);
    ensure_buffer(d_y_, d_y_size_, y_bytes);
    EIGEN_CUDA_RUNTIME_CHECK(cudaMemcpyAsync(d_x_.get(), X.data(), x_bytes, cudaMemcpyHostToDevice, stream_));
    if (beta != Scalar(0)) {
      EIGEN_CUDA_RUNTIME_CHECK(cudaMemcpyAsync(d_y_.get(), Y.data(), y_bytes, cudaMemcpyHostToDevice, stream_));
    }

    constexpr cudaDataType_t dtype = internal::cuda_data_type<Scalar>::value;
    const cusparseOperation_t cu_op = descriptor_op(op);
    cusparseDnMatDescr_t x_desc = nullptr, y_desc = nullptr;
    // Eigen is column-major, so ld = rows.
    EIGEN_CUSPARSE_CHECK(cusparseCreateDnMat(&x_desc, k_op, n, k_op, d_x_.get(), dtype, CUSPARSE_ORDER_COL));
    EIGEN_CUSPARSE_CHECK(cusparseCreateDnMat(&y_desc, m_op, n, m_op, d_y_.get(), dtype, CUSPARSE_ORDER_COL));

    size_t ws_size = 0;
    EIGEN_CUSPARSE_CHECK(cusparseSpMM_bufferSize(handle_, cu_op, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, spmat_desc_,
                                                 x_desc, &beta, y_desc, dtype, kSpMMAlg, &ws_size));
    ensure_buffer(d_workspace_, d_workspace_size_, ws_size);

    EIGEN_CUSPARSE_CHECK(cusparseSpMM(handle_, cu_op, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, spmat_desc_, x_desc,
                                      &beta, y_desc, dtype, kSpMMAlg, d_workspace_.get()));

    EIGEN_CUDA_RUNTIME_CHECK(cudaMemcpyAsync(Y.data(), d_y_.get(), y_bytes, cudaMemcpyDeviceToHost, stream_));
    EIGEN_CUDA_RUNTIME_CHECK(cudaStreamSynchronize(stream_));

    EIGEN_CUSPARSE_CHECK(cusparseDestroyDnMat(x_desc));
    EIGEN_CUSPARSE_CHECK(cusparseDestroyDnMat(y_desc));
  }

  // ---- Helpers --------------------------------------------------------------

  static void check_storage_index_bounds(Index rows, Index cols, Index nnz) {
    const Index max_storage_index = static_cast<Index>((std::numeric_limits<StorageIndex>::max)());
    eigen_assert(rows <= max_storage_index && cols <= max_storage_index && nnz <= max_storage_index &&
                 "gpu::SparseContext currently uses int StorageIndex; matrix dimensions or nonzeros exceed int range");
    EIGEN_UNUSED_VARIABLE(rows);
    EIGEN_UNUSED_VARIABLE(cols);
    EIGEN_UNUSED_VARIABLE(nnz);
    EIGEN_UNUSED_VARIABLE(max_storage_index);
  }

  void upload_sparse(const SpMat& A) {
    // cuSPARSE 12.0+ accepts CSC directly. On cuSPARSE 11.x, cusparseSpMM
    // rejects CSC and CONJUGATE_TRANSPOSE on CSC+complex SpMV silently
    // demotes to TRANSPOSE. We register the same CSC buffers as CSR-of-A^T
    // (dims swapped) on 11.x and invert the op at exec time via
    // descriptor_op() — no transpose-copy required.
    upload_compressed_arrays(A.rows(), A.cols(), A.nonZeros(),
                             /*outer_count=*/A.cols() + 1, A.outerIndexPtr(), A.innerIndexPtr(), A.valuePtr());
  }

  void upload_compressed_arrays(Index m, Index n, Index nnz, Index outer_count, const StorageIndex* host_outer,
                                const StorageIndex* host_inner, const Scalar* host_values) {
    const size_t outer_bytes = static_cast<size_t>(outer_count) * sizeof(StorageIndex);
    const size_t inner_bytes = static_cast<size_t>(nnz) * sizeof(StorageIndex);
    const size_t val_bytes = static_cast<size_t>(nnz) * sizeof(Scalar);

    ensure_buffer(d_outerPtr_, d_outerPtr_size_, outer_bytes);
    ensure_buffer(d_innerIdx_, d_innerIdx_size_, inner_bytes);
    ensure_buffer(d_values_, d_values_size_, val_bytes);

    EIGEN_CUDA_RUNTIME_CHECK(
        cudaMemcpyAsync(d_outerPtr_.get(), host_outer, outer_bytes, cudaMemcpyHostToDevice, stream_));
    EIGEN_CUDA_RUNTIME_CHECK(
        cudaMemcpyAsync(d_innerIdx_.get(), host_inner, inner_bytes, cudaMemcpyHostToDevice, stream_));
    EIGEN_CUDA_RUNTIME_CHECK(cudaMemcpyAsync(d_values_.get(), host_values, val_bytes, cudaMemcpyHostToDevice, stream_));

    if (m != cached_rows_ || n != cached_cols_ || nnz != cached_nnz_) {
      destroy_descriptors_checked();

      constexpr cusparseIndexType_t idx_type = (sizeof(StorageIndex) == 4) ? CUSPARSE_INDEX_32I : CUSPARSE_INDEX_64I;
      constexpr cudaDataType_t val_type = internal::cuda_data_type<Scalar>::value;

      if (kUseCsrOfTranspose) {
        // cuSPARSE 11.x: cusparseSpMM rejects CSC for matA. CSC of A and CSR of
        // A^T share the same buffers, so register the data as CSR-of-A^T (dims
        // swapped) and invert the op in exec_spmv / spmm_impl via descriptor_op.
        EIGEN_CUSPARSE_CHECK(cusparseCreateCsr(&spmat_desc_, n, m, nnz, d_outerPtr_.get(), d_innerIdx_.get(),
                                               d_values_.get(), idx_type, idx_type, CUSPARSE_INDEX_BASE_ZERO,
                                               val_type));
      } else {
        EIGEN_CUSPARSE_CHECK(cusparseCreateCsc(&spmat_desc_, m, n, nnz, d_outerPtr_.get(), d_innerIdx_.get(),
                                               d_values_.get(), idx_type, idx_type, CUSPARSE_INDEX_BASE_ZERO,
                                               val_type));
      }
      cached_rows_ = m;
      cached_cols_ = n;
      cached_nnz_ = nnz;
    } else if (kUseCsrOfTranspose) {
      EIGEN_CUSPARSE_CHECK(cusparseCsrSetPointers(spmat_desc_, d_outerPtr_.get(), d_innerIdx_.get(), d_values_.get()));
    } else {
      EIGEN_CUSPARSE_CHECK(cusparseCscSetPointers(spmat_desc_, d_outerPtr_.get(), d_innerIdx_.get(), d_values_.get()));
    }
  }

  // Destructor-only cleanup: there is no useful recovery path for failures.
  void destroy_descriptors_unchecked() {
    if (spmat_desc_) {
      (void)cusparseDestroySpMat(spmat_desc_);
      spmat_desc_ = nullptr;
    }
    cached_rows_ = -1;
    cached_cols_ = -1;
    cached_nnz_ = -1;
  }

  void destroy_descriptors_checked() {
    if (spmat_desc_) {
      EIGEN_CUSPARSE_CHECK(cusparseDestroySpMat(spmat_desc_));
      spmat_desc_ = nullptr;
    }
    cached_rows_ = -1;
    cached_cols_ = -1;
    cached_nnz_ = -1;
  }

  void ensure_buffer(internal::DeviceBuffer& buf, size_t& current_size, size_t needed) const {
    if (needed > current_size) {
      if (buf) EIGEN_CUDA_RUNTIME_CHECK(cudaStreamSynchronize(stream_));
      buf = internal::DeviceBuffer(needed);
      current_size = needed;
    }
  }
};

// ---- DeviceMatrix::operator=(SpMVExpr) out-of-line definition ----------------
// Defined here because it needs the full SparseContext definition.

template <typename Scalar_>
DeviceMatrix<Scalar_>& DeviceMatrix<Scalar_>::operator=(const SpMVExpr<Scalar_>& expr) {
  // Use spmv_device_exec — the sparse matrix was already uploaded by deviceView().
  // No re-upload on repeated SpMV with the same view.
  expr.view().context().spmv_device_exec(expr.x(), *this, Scalar_(1), Scalar_(0), CUSPARSE_OPERATION_NON_TRANSPOSE);
  return *this;
}

}  // namespace gpu
}  // namespace Eigen

#endif  // EIGEN_GPU_SPARSE_CONTEXT_H
