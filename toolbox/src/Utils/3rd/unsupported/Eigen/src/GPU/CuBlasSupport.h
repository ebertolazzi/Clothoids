// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Rasmus Munk Larsen <rmlarsen@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

// cuBLAS-specific support types:
//   - Error-checking macro
//   - Operation enum and mapping to cublasOperation_t
//
// Generic CUDA runtime utilities (DeviceBuffer, cuda_data_type) are in GpuSupport.h.

#ifndef EIGEN_GPU_CUBLAS_SUPPORT_H
#define EIGEN_GPU_CUBLAS_SUPPORT_H

// IWYU pragma: private
#include "./InternalHeaderCheck.h"

#include "./GpuSupport.h"
#include <cublas_v2.h>
#include <cublasLt.h>
#include <climits>
#include <cstring>
#include <utility>

namespace Eigen {
namespace gpu {
namespace internal {

// ---- Error-checking macro ---------------------------------------------------

#define EIGEN_CUBLAS_CHECK(expr)                                       \
  do {                                                                 \
    cublasStatus_t _s = (expr);                                        \
    eigen_assert(_s == CUBLAS_STATUS_SUCCESS && "cuBLAS call failed"); \
  } while (0)

constexpr cublasOperation_t to_cublas_op(GpuOp op) {
  switch (op) {
    case GpuOp::Trans:
      return CUBLAS_OP_T;
    case GpuOp::ConjTrans:
      return CUBLAS_OP_C;
    default:
      return CUBLAS_OP_N;
  }
}

// ---- Scalar → cublasComputeType_t -------------------------------------------
// cublasLtMatmul requires a compute type (separate from the data type).
//
// Precision policy:
//   - Default: tensor core algorithms enabled via cublasLtMatmul heuristics.
//     For double, cuBLAS may use Ozaki emulation on sm_80+ tensor cores.
//   - EIGEN_CUDA_TF32: opt-in to TF32 for float (~2x faster, 10-bit mantissa).
//   - EIGEN_NO_CUDA_TENSOR_OPS: disables all tensor core usage. Uses pedantic
//     compute types. For bit-exact reproducibility.

// Single-precision (real or complex) and double-precision (real or complex) each
// pick their compute type from one set of preprocessor switches; specializations
// just dispatch to the right precision tag.
namespace cuda_compute_type_detail {
#if defined(EIGEN_NO_CUDA_TENSOR_OPS)
constexpr cublasComputeType_t kFloat = CUBLAS_COMPUTE_32F_PEDANTIC;
constexpr cublasComputeType_t kDouble = CUBLAS_COMPUTE_64F_PEDANTIC;
#elif defined(EIGEN_CUDA_TF32)
constexpr cublasComputeType_t kFloat = CUBLAS_COMPUTE_32F_FAST_TF32;
constexpr cublasComputeType_t kDouble = CUBLAS_COMPUTE_64F;
#else
constexpr cublasComputeType_t kFloat = CUBLAS_COMPUTE_32F;
constexpr cublasComputeType_t kDouble = CUBLAS_COMPUTE_64F;
#endif
}  // namespace cuda_compute_type_detail

template <typename Scalar>
struct cuda_compute_type;

template <>
struct cuda_compute_type<float> {
  static constexpr cublasComputeType_t value = cuda_compute_type_detail::kFloat;
};
template <>
struct cuda_compute_type<double> {
  static constexpr cublasComputeType_t value = cuda_compute_type_detail::kDouble;
};
template <>
struct cuda_compute_type<std::complex<float>> {
  static constexpr cublasComputeType_t value = cuda_compute_type_detail::kFloat;
};
template <>
struct cuda_compute_type<std::complex<double>> {
  static constexpr cublasComputeType_t value = cuda_compute_type_detail::kDouble;
};

// ---- cublasLt GEMM dispatch -------------------------------------------------
// Wraps cublasLtMatmul with descriptor setup, heuristic algorithm selection,
// and a small per-context plan cache. Supports 64-bit dimensions natively.
//
// The workspace buffer (DeviceBuffer*) is grown lazily to match the selected
// algorithm's actual requirement. Growth is monotonic.
//
// EIGEN_NO_CUDA_TENSOR_OPS: pedantic compute types prevent cublasLt from
// selecting tensor core algorithms, matching the previous cublasGemmEx behavior.
//
// Thread safety: the workspace buffer and plan cache are not thread-safe. All
// GEMM calls sharing them must be on the same CUDA stream (guaranteed by
// gpu::Context's single-stream design and by each solver owning its own stream).

#define EIGEN_CUBLASLT_CHECK(expr)                                       \
  do {                                                                   \
    cublasStatus_t _s = (expr);                                          \
    eigen_assert(_s == CUBLAS_STATUS_SUCCESS && "cuBLASLt call failed"); \
  } while (0)

// Maximum workspace the heuristic is allowed to consider. This is a preference
// ceiling, not an allocation — actual allocation matches the selected algorithm.
// Override at compile time via EIGEN_CUDA_CUBLASLT_MAX_WORKSPACE_BYTES.
#ifndef EIGEN_CUDA_CUBLASLT_MAX_WORKSPACE_BYTES
#define EIGEN_CUDA_CUBLASLT_MAX_WORKSPACE_BYTES (32 * 1024 * 1024)  // 32 MB
#endif
static constexpr size_t kCublasLtMaxWorkspaceBytes = EIGEN_CUDA_CUBLASLT_MAX_WORKSPACE_BYTES;

// cublasGemmEx fallback algorithm hint (used when cublasLt heuristic returns no results).
constexpr cublasGemmAlgo_t cuda_gemm_algo() {
#ifdef EIGEN_NO_CUDA_TENSOR_OPS
  return CUBLAS_GEMM_DEFAULT;
#else
  return CUBLAS_GEMM_DEFAULT_TENSOR_OP;
#endif
}

// ---- cublasLt plan cache ----------------------------------------------------
// Caches matmul descriptors, matrix layouts, and the selected algorithm for
// repeated GEMM calls with the same (m, n, k, dtype, transA, transB) shape.
// Eliminates per-call descriptor creation and heuristic lookup overhead, which
// can be 5-35% of total time for small/medium matrices.
//
// Backed by Eigen::internal::LruCache. GEMM shapes in typical workloads (CG
// iteration, chained solves) have very low cardinality (usually 1-3 distinct
// shapes), so the cache is small.

static constexpr std::size_t kCublasLtPlanCacheCapacity = 8;

struct CublasLtPlanKey {
  int64_t m, n, k;
  int64_t lda, ldb, ldc;
  cudaDataType_t dtype;
  cublasOperation_t transA, transB;

  bool operator==(const CublasLtPlanKey& o) const {
    return m == o.m && n == o.n && k == o.k && lda == o.lda && ldb == o.ldb && ldc == o.ldc && dtype == o.dtype &&
           transA == o.transA && transB == o.transB;
  }
};

struct CublasLtPlanKeyHash {
  std::size_t operator()(const CublasLtPlanKey& k) const noexcept {
    // boost-style hash_combine: mix each field into the rolling hash.
    auto mix = [](std::size_t a, std::size_t b) { return a ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2)); };
    std::size_t r = std::hash<int64_t>{}(k.m);
    r = mix(r, std::hash<int64_t>{}(k.n));
    r = mix(r, std::hash<int64_t>{}(k.k));
    r = mix(r, std::hash<int64_t>{}(k.lda));
    r = mix(r, std::hash<int64_t>{}(k.ldb));
    r = mix(r, std::hash<int64_t>{}(k.ldc));
    r = mix(r, std::hash<int>{}(static_cast<int>(k.dtype)));
    r = mix(r, std::hash<int>{}(static_cast<int>(k.transA)));
    r = mix(r, std::hash<int>{}(static_cast<int>(k.transB)));
    return r;
  }
};

// Move-only RAII wrapper for a cached cuBLASLt matmul plan: the descriptor,
// three matrix layouts, and the heuristic-selected algorithm. Destruction
// destroys all cuBLASLt handles, so the cache can manage entry lifetime via
// LruCache's eviction.
//
// Move-only because each instance uniquely owns four cuBLASLt handles
// (matmul_desc and three matrix layouts); copying would alias the handles
// and cause double-destroy.
class CublasLtPlanEntry {
 public:
  // Build descriptors and run the heuristic. If the heuristic returns no
  // usable algorithm, use_cublaslt stays false and the caller takes the
  // cublasGemmEx fallback path. `max_workspace_bytes` is the heuristic's
  // workspace ceiling — see gpu::Context::setCublasLtMaxWorkspaceBytes().
  CublasLtPlanEntry(cublasLtHandle_t lt_handle, const CublasLtPlanKey& key, cublasComputeType_t compute,
                    cudaDataType_t alpha_type, std::size_t max_workspace_bytes) {
    EIGEN_CUBLASLT_CHECK(cublasLtMatmulDescCreate(&matmul_desc, compute, alpha_type));
    EIGEN_CUBLASLT_CHECK(
        cublasLtMatmulDescSetAttribute(matmul_desc, CUBLASLT_MATMUL_DESC_TRANSA, &key.transA, sizeof(key.transA)));
    EIGEN_CUBLASLT_CHECK(
        cublasLtMatmulDescSetAttribute(matmul_desc, CUBLASLT_MATMUL_DESC_TRANSB, &key.transB, sizeof(key.transB)));

    // Layout dimensions are the physical (rows, cols) of the column-major operand;
    // the leading dimension is the actual stride between columns (lda/ldb/ldc),
    // which may exceed the active row count (e.g., a thin view of a wider buffer).
    const int64_t a_rows = (key.transA == CUBLAS_OP_N) ? key.m : key.k;
    const int64_t a_cols = (key.transA == CUBLAS_OP_N) ? key.k : key.m;
    const int64_t b_rows = (key.transB == CUBLAS_OP_N) ? key.k : key.n;
    const int64_t b_cols = (key.transB == CUBLAS_OP_N) ? key.n : key.k;
    EIGEN_CUBLASLT_CHECK(cublasLtMatrixLayoutCreate(&layout_A, key.dtype, a_rows, a_cols, key.lda));
    EIGEN_CUBLASLT_CHECK(cublasLtMatrixLayoutCreate(&layout_B, key.dtype, b_rows, b_cols, key.ldb));
    EIGEN_CUBLASLT_CHECK(cublasLtMatrixLayoutCreate(&layout_C, key.dtype, key.m, key.n, key.ldc));

    cublasLtMatmulPreference_t preference = nullptr;
    EIGEN_CUBLASLT_CHECK(cublasLtMatmulPreferenceCreate(&preference));
    EIGEN_CUBLASLT_CHECK(cublasLtMatmulPreferenceSetAttribute(preference, CUBLASLT_MATMUL_PREF_MAX_WORKSPACE_BYTES,
                                                              &max_workspace_bytes, sizeof(max_workspace_bytes)));

    cublasLtMatmulHeuristicResult_t result;
    int returned_results = 0;
    cublasStatus_t heuristic_status = cublasLtMatmulAlgoGetHeuristic(
        lt_handle, matmul_desc, layout_A, layout_B, layout_C, layout_C, preference, 1, &result, &returned_results);

    EIGEN_CUBLASLT_CHECK(cublasLtMatmulPreferenceDestroy(preference));

    // cublasLtMatmulAlgoGetHeuristic can return CUBLAS_STATUS_SUCCESS overall while
    // marking individual results NOT_SUPPORTED via result.state, so gate on both.
    if (heuristic_status == CUBLAS_STATUS_SUCCESS && returned_results > 0 && result.state == CUBLAS_STATUS_SUCCESS) {
      algo = result.algo;
      workspace_size = result.workspaceSize;
      use_cublaslt = true;
    }
  }

  ~CublasLtPlanEntry() { destroy(); }

  CublasLtPlanEntry(const CublasLtPlanEntry&) = delete;
  CublasLtPlanEntry& operator=(const CublasLtPlanEntry&) = delete;

  CublasLtPlanEntry(CublasLtPlanEntry&& o) noexcept
      : matmul_desc(o.matmul_desc),
        layout_A(o.layout_A),
        layout_B(o.layout_B),
        layout_C(o.layout_C),
        algo(o.algo),
        workspace_size(o.workspace_size),
        use_cublaslt(o.use_cublaslt) {
    o.matmul_desc = nullptr;
    o.layout_A = o.layout_B = o.layout_C = nullptr;
    o.use_cublaslt = false;
  }

  CublasLtPlanEntry& operator=(CublasLtPlanEntry&& o) noexcept {
    if (this != &o) {
      destroy();
      matmul_desc = o.matmul_desc;
      layout_A = o.layout_A;
      layout_B = o.layout_B;
      layout_C = o.layout_C;
      algo = o.algo;
      workspace_size = o.workspace_size;
      use_cublaslt = o.use_cublaslt;
      o.matmul_desc = nullptr;
      o.layout_A = o.layout_B = o.layout_C = nullptr;
      o.use_cublaslt = false;
    }
    return *this;
  }

  // Public read-side state for cublaslt_gemm().
  cublasLtMatmulDesc_t matmul_desc = nullptr;
  cublasLtMatrixLayout_t layout_A = nullptr;
  cublasLtMatrixLayout_t layout_B = nullptr;
  cublasLtMatrixLayout_t layout_C = nullptr;
  cublasLtMatmulAlgo_t algo{};
  std::size_t workspace_size = 0;
  bool use_cublaslt = false;

 private:
  void destroy() noexcept {
    if (layout_C) cublasLtMatrixLayoutDestroy(layout_C);
    if (layout_B) cublasLtMatrixLayoutDestroy(layout_B);
    if (layout_A) cublasLtMatrixLayoutDestroy(layout_A);
    if (matmul_desc) cublasLtMatmulDescDestroy(matmul_desc);
  }
};

using CublasLtPlanCache = Eigen::internal::LruCache<CublasLtPlanKey, CublasLtPlanEntry, CublasLtPlanKeyHash>;

// cublasLtMatmul GEMM with shape-keyed plan cache and lazy workspace.
//
// Falls back to cublasGemmEx for shapes/types where the cublasLt heuristic
// returns no results.
template <typename Scalar>
void cublaslt_gemm(cublasLtHandle_t lt_handle, cublasHandle_t cublas_handle, cublasOperation_t transA,
                   cublasOperation_t transB, int64_t m, int64_t n, int64_t k, const Scalar* alpha, const Scalar* A,
                   int64_t lda, const Scalar* B, int64_t ldb, const Scalar* beta, Scalar* C, int64_t ldc,
                   DeviceBuffer* workspace, CublasLtPlanCache* plan_cache, std::size_t max_workspace_bytes,
                   cudaStream_t stream) {
  constexpr cudaDataType_t dtype = cuda_data_type<Scalar>::value;
  constexpr cublasComputeType_t compute = cuda_compute_type<Scalar>::value;
  constexpr cudaDataType_t alpha_type = cuda_data_type<Scalar>::value;

  // Look up or create a cached plan for this shape (key includes leading dims so
  // strided views — e.g. SVD's thin VT/U slices — get distinct cache entries).
  const CublasLtPlanKey key{m, n, k, lda, ldb, ldc, dtype, transA, transB};
  CublasLtPlanEntry* entry = plan_cache->find(key);
  if (!entry) {
    entry = plan_cache->insert(key, CublasLtPlanEntry(lt_handle, key, compute, alpha_type, max_workspace_bytes));
  }

  if (entry->use_cublaslt) {
    const size_t needed = entry->workspace_size;
    if (needed > workspace->size()) {
      // Sync only when freeing an existing buffer that may be in use.
      if (workspace->get()) EIGEN_CUDA_RUNTIME_CHECK(cudaStreamSynchronize(stream));
      *workspace = DeviceBuffer(needed);
    }

    EIGEN_CUBLASLT_CHECK(cublasLtMatmul(lt_handle, entry->matmul_desc, alpha, A, entry->layout_A, B, entry->layout_B,
                                        beta, C, entry->layout_C, C, entry->layout_C, &entry->algo, workspace->get(),
                                        needed, stream));
  } else {
    // Fallback: cublasGemmEx for shapes/types that cublasLt cannot handle.
    // cublasGemmEx takes int dimensions; cublaslt_gemm itself supports int64_t.
    eigen_assert(m <= INT_MAX && n <= INT_MAX && k <= INT_MAX && lda <= INT_MAX && ldb <= INT_MAX && ldc <= INT_MAX &&
                 "cublasGemmEx fallback: dimensions exceed int range");
    EIGEN_CUBLAS_CHECK(cublasGemmEx(cublas_handle, transA, transB, static_cast<int>(m), static_cast<int>(n),
                                    static_cast<int>(k), alpha, A, dtype, static_cast<int>(lda), B, dtype,
                                    static_cast<int>(ldb), beta, C, dtype, static_cast<int>(ldc), compute,
                                    cuda_gemm_algo()));
  }
}

// ---- Type-specific cuBLAS wrappers ------------------------------------------
// cuBLAS uses separate functions per type (Sgemm, Dgemm, etc.).
// These overloaded wrappers allow calling cublasXgemm/cublasXtrsm/etc.
// with any supported scalar type.

// GEMM wrappers
inline cublasStatus_t cublasXgemm(cublasHandle_t h, cublasOperation_t transA, cublasOperation_t transB, int m, int n,
                                  int k, const float* alpha, const float* A, int lda, const float* B, int ldb,
                                  const float* beta, float* C, int ldc) {
  return cublasSgemm(h, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
inline cublasStatus_t cublasXgemm(cublasHandle_t h, cublasOperation_t transA, cublasOperation_t transB, int m, int n,
                                  int k, const double* alpha, const double* A, int lda, const double* B, int ldb,
                                  const double* beta, double* C, int ldc) {
  return cublasDgemm(h, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
static_assert(sizeof(cuComplex) == sizeof(std::complex<float>), "cuComplex and std::complex<float> layout mismatch");
static_assert(sizeof(cuDoubleComplex) == sizeof(std::complex<double>),
              "cuDoubleComplex and std::complex<double> layout mismatch");

// Complex scalar args (alpha, beta) are type-punned from std::complex<T>*
// to cuComplex*/cuDoubleComplex*.  A reinterpret_cast violates strict
// aliasing: when inlined, clang/MSVC can elide the caller's store (the
// compiler no longer sees a read through the original type), causing
// segfaults.  We use memcpy — the standard-blessed type-pun — for scalars.
// Device array pointers (A, B, C) are opaque to the host compiler, so
// reinterpret_cast is safe there.
inline cublasStatus_t cublasXgemm(cublasHandle_t h, cublasOperation_t transA, cublasOperation_t transB, int m, int n,
                                  int k, const std::complex<float>* alpha, const std::complex<float>* A, int lda,
                                  const std::complex<float>* B, int ldb, const std::complex<float>* beta,
                                  std::complex<float>* C, int ldc) {
  cuComplex a, b;
  std::memcpy(&a, alpha, sizeof(a));
  std::memcpy(&b, beta, sizeof(b));
  return cublasCgemm(h, transA, transB, m, n, k, &a, reinterpret_cast<const cuComplex*>(A), lda,
                     reinterpret_cast<const cuComplex*>(B), ldb, &b, reinterpret_cast<cuComplex*>(C), ldc);
}
inline cublasStatus_t cublasXgemm(cublasHandle_t h, cublasOperation_t transA, cublasOperation_t transB, int m, int n,
                                  int k, const std::complex<double>* alpha, const std::complex<double>* A, int lda,
                                  const std::complex<double>* B, int ldb, const std::complex<double>* beta,
                                  std::complex<double>* C, int ldc) {
  cuDoubleComplex a, b;
  std::memcpy(&a, alpha, sizeof(a));
  std::memcpy(&b, beta, sizeof(b));
  return cublasZgemm(h, transA, transB, m, n, k, &a, reinterpret_cast<const cuDoubleComplex*>(A), lda,
                     reinterpret_cast<const cuDoubleComplex*>(B), ldb, &b, reinterpret_cast<cuDoubleComplex*>(C), ldc);
}

// TRSM wrappers
inline cublasStatus_t cublasXtrsm(cublasHandle_t h, cublasSideMode_t side, cublasFillMode_t uplo,
                                  cublasOperation_t trans, cublasDiagType_t diag, int m, int n, const float* alpha,
                                  const float* A, int lda, float* B, int ldb) {
  return cublasStrsm(h, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb);
}
inline cublasStatus_t cublasXtrsm(cublasHandle_t h, cublasSideMode_t side, cublasFillMode_t uplo,
                                  cublasOperation_t trans, cublasDiagType_t diag, int m, int n, const double* alpha,
                                  const double* A, int lda, double* B, int ldb) {
  return cublasDtrsm(h, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb);
}
inline cublasStatus_t cublasXtrsm(cublasHandle_t h, cublasSideMode_t side, cublasFillMode_t uplo,
                                  cublasOperation_t trans, cublasDiagType_t diag, int m, int n,
                                  const std::complex<float>* alpha, const std::complex<float>* A, int lda,
                                  std::complex<float>* B, int ldb) {
  cuComplex a;
  std::memcpy(&a, alpha, sizeof(a));
  return cublasCtrsm(h, side, uplo, trans, diag, m, n, &a, reinterpret_cast<const cuComplex*>(A), lda,
                     reinterpret_cast<cuComplex*>(B), ldb);
}
inline cublasStatus_t cublasXtrsm(cublasHandle_t h, cublasSideMode_t side, cublasFillMode_t uplo,
                                  cublasOperation_t trans, cublasDiagType_t diag, int m, int n,
                                  const std::complex<double>* alpha, const std::complex<double>* A, int lda,
                                  std::complex<double>* B, int ldb) {
  cuDoubleComplex a;
  std::memcpy(&a, alpha, sizeof(a));
  return cublasZtrsm(h, side, uplo, trans, diag, m, n, &a, reinterpret_cast<const cuDoubleComplex*>(A), lda,
                     reinterpret_cast<cuDoubleComplex*>(B), ldb);
}

// SYMM wrappers (real → symm, complex → hemm)
inline cublasStatus_t cublasXsymm(cublasHandle_t h, cublasSideMode_t side, cublasFillMode_t uplo, int m, int n,
                                  const float* alpha, const float* A, int lda, const float* B, int ldb,
                                  const float* beta, float* C, int ldc) {
  return cublasSsymm(h, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}
inline cublasStatus_t cublasXsymm(cublasHandle_t h, cublasSideMode_t side, cublasFillMode_t uplo, int m, int n,
                                  const double* alpha, const double* A, int lda, const double* B, int ldb,
                                  const double* beta, double* C, int ldc) {
  return cublasDsymm(h, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}
inline cublasStatus_t cublasXsymm(cublasHandle_t h, cublasSideMode_t side, cublasFillMode_t uplo, int m, int n,
                                  const std::complex<float>* alpha, const std::complex<float>* A, int lda,
                                  const std::complex<float>* B, int ldb, const std::complex<float>* beta,
                                  std::complex<float>* C, int ldc) {
  cuComplex a, b;
  std::memcpy(&a, alpha, sizeof(a));
  std::memcpy(&b, beta, sizeof(b));
  return cublasChemm(h, side, uplo, m, n, &a, reinterpret_cast<const cuComplex*>(A), lda,
                     reinterpret_cast<const cuComplex*>(B), ldb, &b, reinterpret_cast<cuComplex*>(C), ldc);
}
inline cublasStatus_t cublasXsymm(cublasHandle_t h, cublasSideMode_t side, cublasFillMode_t uplo, int m, int n,
                                  const std::complex<double>* alpha, const std::complex<double>* A, int lda,
                                  const std::complex<double>* B, int ldb, const std::complex<double>* beta,
                                  std::complex<double>* C, int ldc) {
  cuDoubleComplex a, b;
  std::memcpy(&a, alpha, sizeof(a));
  std::memcpy(&b, beta, sizeof(b));
  return cublasZhemm(h, side, uplo, m, n, &a, reinterpret_cast<const cuDoubleComplex*>(A), lda,
                     reinterpret_cast<const cuDoubleComplex*>(B), ldb, &b, reinterpret_cast<cuDoubleComplex*>(C), ldc);
}

// GEAM wrappers: C = alpha * op(A) + beta * op(B)
inline cublasStatus_t cublasXgeam(cublasHandle_t h, cublasOperation_t transA, cublasOperation_t transB, int m, int n,
                                  const float* alpha, const float* A, int lda, const float* beta, const float* B,
                                  int ldb, float* C, int ldc) {
  return cublasSgeam(h, transA, transB, m, n, alpha, A, lda, beta, B, ldb, C, ldc);
}
inline cublasStatus_t cublasXgeam(cublasHandle_t h, cublasOperation_t transA, cublasOperation_t transB, int m, int n,
                                  const double* alpha, const double* A, int lda, const double* beta, const double* B,
                                  int ldb, double* C, int ldc) {
  return cublasDgeam(h, transA, transB, m, n, alpha, A, lda, beta, B, ldb, C, ldc);
}
inline cublasStatus_t cublasXgeam(cublasHandle_t h, cublasOperation_t transA, cublasOperation_t transB, int m, int n,
                                  const std::complex<float>* alpha, const std::complex<float>* A, int lda,
                                  const std::complex<float>* beta, const std::complex<float>* B, int ldb,
                                  std::complex<float>* C, int ldc) {
  cuComplex a, b;
  std::memcpy(&a, alpha, sizeof(a));
  std::memcpy(&b, beta, sizeof(b));
  return cublasCgeam(h, transA, transB, m, n, &a, reinterpret_cast<const cuComplex*>(A), lda, &b,
                     reinterpret_cast<const cuComplex*>(B), ldb, reinterpret_cast<cuComplex*>(C), ldc);
}
inline cublasStatus_t cublasXgeam(cublasHandle_t h, cublasOperation_t transA, cublasOperation_t transB, int m, int n,
                                  const std::complex<double>* alpha, const std::complex<double>* A, int lda,
                                  const std::complex<double>* beta, const std::complex<double>* B, int ldb,
                                  std::complex<double>* C, int ldc) {
  cuDoubleComplex a, b;
  std::memcpy(&a, alpha, sizeof(a));
  std::memcpy(&b, beta, sizeof(b));
  return cublasZgeam(h, transA, transB, m, n, &a, reinterpret_cast<const cuDoubleComplex*>(A), lda, &b,
                     reinterpret_cast<const cuDoubleComplex*>(B), ldb, reinterpret_cast<cuDoubleComplex*>(C), ldc);
}

// SYRK wrappers (real → syrk, complex → herk)
inline cublasStatus_t cublasXsyrk(cublasHandle_t h, cublasFillMode_t uplo, cublasOperation_t trans, int n, int k,
                                  const float* alpha, const float* A, int lda, const float* beta, float* C, int ldc) {
  return cublasSsyrk(h, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}
inline cublasStatus_t cublasXsyrk(cublasHandle_t h, cublasFillMode_t uplo, cublasOperation_t trans, int n, int k,
                                  const double* alpha, const double* A, int lda, const double* beta, double* C,
                                  int ldc) {
  return cublasDsyrk(h, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}
inline cublasStatus_t cublasXsyrk(cublasHandle_t h, cublasFillMode_t uplo, cublasOperation_t trans, int n, int k,
                                  const float* alpha, const std::complex<float>* A, int lda, const float* beta,
                                  std::complex<float>* C, int ldc) {
  return cublasCherk(h, uplo, trans, n, k, alpha, reinterpret_cast<const cuComplex*>(A), lda, beta,
                     reinterpret_cast<cuComplex*>(C), ldc);
}
inline cublasStatus_t cublasXsyrk(cublasHandle_t h, cublasFillMode_t uplo, cublasOperation_t trans, int n, int k,
                                  const double* alpha, const std::complex<double>* A, int lda, const double* beta,
                                  std::complex<double>* C, int ldc) {
  return cublasZherk(h, uplo, trans, n, k, alpha, reinterpret_cast<const cuDoubleComplex*>(A), lda, beta,
                     reinterpret_cast<cuDoubleComplex*>(C), ldc);
}

// SCAL wrappers: x = alpha * x.
// For complex x, alpha is real-valued (Csscal/Zdscal) — this matches the
// 1/n inverse-FFT scaling pattern, where the scale is intrinsically real.
inline cublasStatus_t cublasXscal(cublasHandle_t h, int n, const float* alpha, float* x, int incx) {
  return cublasSscal(h, n, alpha, x, incx);
}
inline cublasStatus_t cublasXscal(cublasHandle_t h, int n, const double* alpha, double* x, int incx) {
  return cublasDscal(h, n, alpha, x, incx);
}
inline cublasStatus_t cublasXscal(cublasHandle_t h, int n, const float* alpha, std::complex<float>* x, int incx) {
  return cublasCsscal(h, n, alpha, reinterpret_cast<cuComplex*>(x), incx);
}
inline cublasStatus_t cublasXscal(cublasHandle_t h, int n, const double* alpha, std::complex<double>* x, int incx) {
  return cublasZdscal(h, n, alpha, reinterpret_cast<cuDoubleComplex*>(x), incx);
}

// By-value alpha overloads: convenience for callers that hold the scale as a
// scalar rather than a host pointer (e.g. inverse-FFT 1/n normalization).
inline cublasStatus_t cublasXscal(cublasHandle_t h, int n, float alpha, float* x, int incx) {
  return cublasSscal(h, n, &alpha, x, incx);
}
inline cublasStatus_t cublasXscal(cublasHandle_t h, int n, double alpha, double* x, int incx) {
  return cublasDscal(h, n, &alpha, x, incx);
}
inline cublasStatus_t cublasXscal(cublasHandle_t h, int n, float alpha, std::complex<float>* x, int incx) {
  return cublasCsscal(h, n, &alpha, reinterpret_cast<cuComplex*>(x), incx);
}
inline cublasStatus_t cublasXscal(cublasHandle_t h, int n, double alpha, std::complex<double>* x, int incx) {
  return cublasZdscal(h, n, &alpha, reinterpret_cast<cuDoubleComplex*>(x), incx);
}

// DGMM wrappers: C = A * diag(x)  (side=RIGHT) or C = diag(x) * A  (side=LEFT).
// Useful for applying a diagonal scaling without materialising diag(x) as a
// dense matrix. cuBLAS docs guarantee in-place is safe when C == A.
inline cublasStatus_t cublasXdgmm(cublasHandle_t h, cublasSideMode_t side, int m, int n, const float* A, int lda,
                                  const float* x, int incx, float* C, int ldc) {
  return cublasSdgmm(h, side, m, n, A, lda, x, incx, C, ldc);
}
inline cublasStatus_t cublasXdgmm(cublasHandle_t h, cublasSideMode_t side, int m, int n, const double* A, int lda,
                                  const double* x, int incx, double* C, int ldc) {
  return cublasDdgmm(h, side, m, n, A, lda, x, incx, C, ldc);
}
inline cublasStatus_t cublasXdgmm(cublasHandle_t h, cublasSideMode_t side, int m, int n, const std::complex<float>* A,
                                  int lda, const std::complex<float>* x, int incx, std::complex<float>* C, int ldc) {
  return cublasCdgmm(h, side, m, n, reinterpret_cast<const cuComplex*>(A), lda, reinterpret_cast<const cuComplex*>(x),
                     incx, reinterpret_cast<cuComplex*>(C), ldc);
}
inline cublasStatus_t cublasXdgmm(cublasHandle_t h, cublasSideMode_t side, int m, int n, const std::complex<double>* A,
                                  int lda, const std::complex<double>* x, int incx, std::complex<double>* C, int ldc) {
  return cublasZdgmm(h, side, m, n, reinterpret_cast<const cuDoubleComplex*>(A), lda,
                     reinterpret_cast<const cuDoubleComplex*>(x), incx, reinterpret_cast<cuDoubleComplex*>(C), ldc);
}

// ---- cuBLAS Level-1 wrappers ------------------------------------------------
// Type-dispatched wrappers for BLAS-1 vector operations: dot, axpy, nrm2, scal, copy.
// These work with CUBLAS_POINTER_MODE_HOST or CUBLAS_POINTER_MODE_DEVICE depending
// on the caller's configuration. For device pointer mode, scalar result pointers
// (dot, nrm2) must point to device memory.

// dot: result = x^T * y (real) or x^H * y (complex conjugate dot)
inline cublasStatus_t cublasXdot(cublasHandle_t h, int n, const float* x, int incx, const float* y, int incy,
                                 float* result) {
  return cublasSdot(h, n, x, incx, y, incy, result);
}
inline cublasStatus_t cublasXdot(cublasHandle_t h, int n, const double* x, int incx, const double* y, int incy,
                                 double* result) {
  return cublasDdot(h, n, x, incx, y, incy, result);
}
inline cublasStatus_t cublasXdot(cublasHandle_t h, int n, const std::complex<float>* x, int incx,
                                 const std::complex<float>* y, int incy, std::complex<float>* result) {
  return cublasCdotc(h, n, reinterpret_cast<const cuComplex*>(x), incx, reinterpret_cast<const cuComplex*>(y), incy,
                     reinterpret_cast<cuComplex*>(result));
}
inline cublasStatus_t cublasXdot(cublasHandle_t h, int n, const std::complex<double>* x, int incx,
                                 const std::complex<double>* y, int incy, std::complex<double>* result) {
  return cublasZdotc(h, n, reinterpret_cast<const cuDoubleComplex*>(x), incx,
                     reinterpret_cast<const cuDoubleComplex*>(y), incy, reinterpret_cast<cuDoubleComplex*>(result));
}

// nrm2: result = ||x||_2 (always returns real)
inline cublasStatus_t cublasXnrm2(cublasHandle_t h, int n, const float* x, int incx, float* result) {
  return cublasSnrm2(h, n, x, incx, result);
}
inline cublasStatus_t cublasXnrm2(cublasHandle_t h, int n, const double* x, int incx, double* result) {
  return cublasDnrm2(h, n, x, incx, result);
}
inline cublasStatus_t cublasXnrm2(cublasHandle_t h, int n, const std::complex<float>* x, int incx, float* result) {
  return cublasScnrm2(h, n, reinterpret_cast<const cuComplex*>(x), incx, result);
}
inline cublasStatus_t cublasXnrm2(cublasHandle_t h, int n, const std::complex<double>* x, int incx, double* result) {
  return cublasDznrm2(h, n, reinterpret_cast<const cuDoubleComplex*>(x), incx, result);
}

// axpy: y += alpha * x
inline cublasStatus_t cublasXaxpy(cublasHandle_t h, int n, const float* alpha, const float* x, int incx, float* y,
                                  int incy) {
  return cublasSaxpy(h, n, alpha, x, incx, y, incy);
}
inline cublasStatus_t cublasXaxpy(cublasHandle_t h, int n, const double* alpha, const double* x, int incx, double* y,
                                  int incy) {
  return cublasDaxpy(h, n, alpha, x, incx, y, incy);
}
inline cublasStatus_t cublasXaxpy(cublasHandle_t h, int n, const std::complex<float>* alpha,
                                  const std::complex<float>* x, int incx, std::complex<float>* y, int incy) {
  cuComplex a;
  std::memcpy(&a, alpha, sizeof(a));
  return cublasCaxpy(h, n, &a, reinterpret_cast<const cuComplex*>(x), incx, reinterpret_cast<cuComplex*>(y), incy);
}
inline cublasStatus_t cublasXaxpy(cublasHandle_t h, int n, const std::complex<double>* alpha,
                                  const std::complex<double>* x, int incx, std::complex<double>* y, int incy) {
  cuDoubleComplex a;
  std::memcpy(&a, alpha, sizeof(a));
  return cublasZaxpy(h, n, &a, reinterpret_cast<const cuDoubleComplex*>(x), incx, reinterpret_cast<cuDoubleComplex*>(y),
                     incy);
}

// SCAL with complex alpha on complex vectors (Cscal/Zscal). The real-alpha
// overloads (Sscal/Dscal/Csscal/Zdscal) live above with the FFT-scaling forms.
inline cublasStatus_t cublasXscal(cublasHandle_t h, int n, const std::complex<float>* alpha, std::complex<float>* x,
                                  int incx) {
  cuComplex a;
  std::memcpy(&a, alpha, sizeof(a));
  return cublasCscal(h, n, &a, reinterpret_cast<cuComplex*>(x), incx);
}
inline cublasStatus_t cublasXscal(cublasHandle_t h, int n, const std::complex<double>* alpha, std::complex<double>* x,
                                  int incx) {
  cuDoubleComplex a;
  std::memcpy(&a, alpha, sizeof(a));
  return cublasZscal(h, n, &a, reinterpret_cast<cuDoubleComplex*>(x), incx);
}

// copy: y = x
inline cublasStatus_t cublasXcopy(cublasHandle_t h, int n, const float* x, int incx, float* y, int incy) {
  return cublasScopy(h, n, x, incx, y, incy);
}
inline cublasStatus_t cublasXcopy(cublasHandle_t h, int n, const double* x, int incx, double* y, int incy) {
  return cublasDcopy(h, n, x, incx, y, incy);
}
inline cublasStatus_t cublasXcopy(cublasHandle_t h, int n, const std::complex<float>* x, int incx,
                                  std::complex<float>* y, int incy) {
  return cublasCcopy(h, n, reinterpret_cast<const cuComplex*>(x), incx, reinterpret_cast<cuComplex*>(y), incy);
}
inline cublasStatus_t cublasXcopy(cublasHandle_t h, int n, const std::complex<double>* x, int incx,
                                  std::complex<double>* y, int incy) {
  return cublasZcopy(h, n, reinterpret_cast<const cuDoubleComplex*>(x), incx, reinterpret_cast<cuDoubleComplex*>(y),
                     incy);
}

}  // namespace internal
}  // namespace gpu
}  // namespace Eigen

#endif  // EIGEN_GPU_CUBLAS_SUPPORT_H
