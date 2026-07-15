// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Rasmus Munk Larsen <rmlarsen@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

#ifndef EIGEN_THREADED_SPARSE_PRODUCT_H
#define EIGEN_THREADED_SPARSE_PRODUCT_H

// IWYU pragma: private
#include "./InternalHeaderCheck.h"

namespace Eigen {

namespace internal {

inline ThreadPool& default_threaded_sparse_pool() {
  static ThreadPool pool(numext::maxi<unsigned>(1u, std::thread::hardware_concurrency()));
  return pool;
}

// nnz-balanced partition of an outer range [0, outerSize) into numChunks
// contiguous chunks. boundaries[t] = first outer index owned by partition t;
// boundaries[numChunks] = outerSize.
//
// The split uses std::lower_bound on the outer-index array. Targets are
// monotonically increasing in t, so each search starts from the previous
// boundary; total work is bounded by O(numChunks + log outerSize) rather than
// numChunks * log(outerSize). Each chunk's nnz count differs from the ideal by
// at most max_nnz_per_outer.
template <typename StorageIndex>
inline void compute_nnz_balanced_partition(const StorageIndex* outer, Index outerSize, Index totalNnz, int numChunks,
                                           std::vector<Index>& boundaries) {
  boundaries.assign(numChunks + 1, 0);
  boundaries[numChunks] = outerSize;
  if (numChunks <= 1 || outerSize == 0 || totalNnz == 0) return;
  const StorageIndex* const last = outer + outerSize + 1;
  const StorageIndex* lo = outer;
  for (int t = 1; t < numChunks; ++t) {
    Index target = (static_cast<Index>(t) * totalNnz) / numChunks;
    lo = std::lower_bound(lo, last, static_cast<StorageIndex>(target));
    boundaries[t] = lo - outer;
  }
}

// Single-row dot-product kernel. Used as the body of a per-row OpenMP
// `parallel for` (which can do its own dynamic scheduling) and from within
// run_dot_chunk for the ThreadPool dispatch path. Marked ALWAYS_INLINE so
// the loop body is visible to the OMP iteration scheduler for vectorization.
template <bool Conjugate, bool Overwrite, typename Scalar, typename StorageIndex, typename XScalar, typename YScalar,
          typename AlphaT>
EIGEN_ALWAYS_INLINE void run_dot_row(const Scalar* EIGEN_RESTRICT vals, const StorageIndex* EIGEN_RESTRICT inner,
                                     const StorageIndex* EIGEN_RESTRICT outer,
                                     const StorageIndex* EIGEN_RESTRICT innerNnz, const XScalar* EIGEN_RESTRICT x,
                                     YScalar* EIGEN_RESTRICT y, Index i, const AlphaT& alpha) {
  const Index k0 = outer[i];
  const Index end = innerNnz ? (k0 + innerNnz[i]) : outer[i + 1];
  Scalar s0(0), s1(0);
  Index k = k0;
  const conj_if<Conjugate> cj{};
  // Loop structure mirrors the existing kernel in SparseDenseProduct.h so
  // the compiler vectorizes identically.
  for (; k < end; ++k) {
    s0 += cj(vals[k]) * x[inner[k]];
    ++k;
    if (k < end) s1 += cj(vals[k]) * x[inner[k]];
  }
  const Scalar s = s0 + s1;
  EIGEN_IF_CONSTEXPR (Overwrite) {
    y[i] = alpha * s;
  } else {
    y[i] += alpha * s;
  }
}

// Dot-product-per-outer kernel. Used by both forward and adjoint paths; the
// adjoint path sets Conjugate=true (no-op for real Scalar).
//
// Processes outer indices [lo, hi). For each outer i:
//   if Overwrite: y[i]  = alpha * sum_k op(val[k]) * x[inner[k]]
//   else:         y[i] += alpha * sum_k op(val[k]) * x[inner[k]]
// where op(.) is conj(.) iff Conjugate.
//
// Writes are independent across threads (each thread owns a contiguous output
// range), so no synchronization is required.
template <bool Conjugate, bool Overwrite, typename Scalar, typename StorageIndex, typename XScalar, typename YScalar,
          typename AlphaT>
EIGEN_STRONG_INLINE void run_dot_chunk(const Scalar* EIGEN_RESTRICT vals, const StorageIndex* EIGEN_RESTRICT inner,
                                       const StorageIndex* EIGEN_RESTRICT outer,
                                       const StorageIndex* EIGEN_RESTRICT innerNnz, const XScalar* EIGEN_RESTRICT x,
                                       YScalar* EIGEN_RESTRICT y, Index lo, Index hi, const AlphaT& alpha) {
  for (Index i = lo; i < hi; ++i) {
    run_dot_row<Conjugate, Overwrite, Scalar, StorageIndex>(vals, inner, outer, innerNnz, x, y, i, alpha);
  }
}

}  // namespace internal

/** \class ThreadedSparseProduct
 * \ingroup SparseCore_Module
 *
 * \brief Cached, thread-parallel sparse matrix * dense vector product.
 *
 * Designed for iterative solvers (CG, BiCGSTAB, GMRES, LSCG) that multiply
 * many times by the same sparse matrix and, in the case of LSCG, by its
 * adjoint. The constructor (or analyzePattern()) computes an nnz-balanced
 * row partition once and reuses it across every apply()/applyAdjoint() call.
 *
 * The forward direction (apply) always runs a dot-product-per-outer kernel,
 * which is embarrassingly parallel (no write conflicts). For matrices whose
 * native storage order doesn't match the forward direction (ColMajor for
 * forward, or RowMajor for adjoint), a transposed mirror is built lazily on
 * first use. The mirror doubles the memory footprint of the matrix but
 * permits conflict-free parallel writes for both directions.
 *
 * The matrix is held by const reference (no copy). Caller is responsible for
 * keeping it alive across apply calls, and for invalidating the cache when
 * the matrix changes:
 *   - sparsity-pattern change: call analyzePattern()/compute() to rebuild
 *     the partition and drop the mirror.
 *   - coefficient-only change (same pattern): call refreshValues() to
 *     drop the (now-stale) mirror; the next mirror-using direction will
 *     rebuild it lazily. Without this call, mirror-backed directions
 *     would read frozen coefficients from the cached copy.
 *
 * Aliasing: apply()/applyAdjoint() do NOT insert temporaries on overlap
 * between x and y; callers must use distinct storage (the kernel writes
 * y[i] from a sum over x[inner[k]] reads, which would mis-compute if
 * y aliases x). Asserted in debug builds.
 *
 * Thread-safety: the const apply()/applyAdjoint()/applyAddTo()/applyAdjointAddTo()
 * methods are safe to call concurrently on the same operator -- the lazy mirror
 * is published through an atomic with double-checked locking. The non-const
 * methods (analyzePattern()/compute()/refreshValues()) and destruction are NOT;
 * they mutate or free the cached mirror, so they must not overlap any in-flight
 * apply() (the same rule as calling any method concurrently with the destructor).
 *
 * Layout: x and y are taken as Ref<const DenseVector> / Ref<DenseVector>.
 * Plain vectors and unit-inner-stride views (e.g. a column of a ColMajor
 * matrix, an unaligned Map of a contiguous buffer) bind without copying;
 * a strided const input is copy-evaluated by Ref; a strided mutable output
 * triggers a Ref runtime error.
 *
 * \tparam SparseMatrixType A compressed Eigen::SparseMatrix<Scalar, Order, StorageIndex>.
 */
template <typename SparseMatrixType_>
class ThreadedSparseProduct {
 public:
  typedef SparseMatrixType_ SparseMatrixType;
  typedef typename SparseMatrixType::Scalar Scalar;
  typedef typename SparseMatrixType::RealScalar RealScalar;
  typedef typename SparseMatrixType::StorageIndex StorageIndex;

  enum { IsRowMajor = static_cast<int>(SparseMatrixType::IsRowMajor) };

  // Opposite-storage-order mirror used by the conflict-free path for whichever
  // direction doesn't match A's native order.
  typedef SparseMatrix<Scalar, IsRowMajor ? ColMajor : RowMajor, StorageIndex> MirrorType;

  // The dot-product kernel walks x/y via raw `Scalar*` with unit inner
  // stride; constraining the public API to `Ref` forces compatible layout
  // (binding for plain matrices, vectors, and unit-stride Maps/Blocks;
  // copy-evaluation for non-conforming const inputs; runtime error for
  // non-conforming mutable outputs).
  typedef Matrix<Scalar, Dynamic, 1> DenseVector;
  typedef Ref<const DenseVector> ConstVectorRef;
  typedef Ref<DenseVector> MutableVectorRef;

  ThreadedSparseProduct() = default;

  explicit ThreadedSparseProduct(const SparseMatrixType& mat, ThreadPool* pool = nullptr) : m_pool(pool) {
    analyzePattern(mat);
  }

  ~ThreadedSparseProduct() { delete m_mirror.load(std::memory_order_acquire); }

  ThreadedSparseProduct(const ThreadedSparseProduct&) = delete;
  ThreadedSparseProduct& operator=(const ThreadedSparseProduct&) = delete;

  /** Bind to \a mat and (re)compute the partition. Drops any previously built
   * adjoint mirror. */
  ThreadedSparseProduct& analyzePattern(const SparseMatrixType& mat) {
    eigen_assert(mat.isCompressed() && "ThreadedSparseProduct requires a compressed SparseMatrix");
    m_mat = &mat;
    delete m_mirror.exchange(nullptr, std::memory_order_acq_rel);
    m_adj_partition.clear();

    build_partition(mat.outerIndexPtr(), mat.outerSize(), mat.nonZeros(), m_native_partition);
    return *this;
  }

  ThreadedSparseProduct& compute(const SparseMatrixType& mat) { return analyzePattern(mat); }

  /** Drops the cached adjoint mirror without recomputing the partition. Use
   * after coefficient-only updates to the bound matrix so the next
   * applyAdjoint() (or forward apply on a ColMajor A) rebuilds the mirror
   * with the new values. */
  ThreadedSparseProduct& refreshValues() {
    delete m_mirror.exchange(nullptr, std::memory_order_acq_rel);
    return *this;
  }

  Index rows() const { return m_mat ? m_mat->rows() : Index(0); }
  Index cols() const { return m_mat ? m_mat->cols() : Index(0); }

  /// Overwriting forward apply: y = A * x.
  void apply(const ConstVectorRef& x, MutableVectorRef y) const {
    apply_impl<false, /*Overwrite=*/true>(x, y, Scalar(1));
  }
  /// Overwriting adjoint apply: y = A^H * x.
  void applyAdjoint(const ConstVectorRef& x, MutableVectorRef y) const {
    apply_impl<true, /*Overwrite=*/true>(x, y, Scalar(1));
  }
  /// Accumulating forward apply: y += alpha * A * x.
  void applyAddTo(const ConstVectorRef& x, MutableVectorRef y, const Scalar& alpha) const {
    apply_impl<false, /*Overwrite=*/false>(x, y, alpha);
  }
  /// Accumulating adjoint apply: y += alpha * A^H * x.
  void applyAdjointAddTo(const ConstVectorRef& x, MutableVectorRef y, const Scalar& alpha) const {
    apply_impl<true, /*Overwrite=*/false>(x, y, alpha);
  }

  /// Returns the thread pool used by this operator.
  ThreadPool* pool() const { return m_pool ? m_pool : &internal::default_threaded_sparse_pool(); }

  /// True iff the lazy adjoint mirror has been materialized.
  bool hasMirror() const { return m_mirror.load(std::memory_order_acquire) != nullptr; }

 private:
  template <bool Adjoint, bool Overwrite>
  void apply_impl(const ConstVectorRef& x, MutableVectorRef& y, const Scalar& alpha) const {
    eigen_assert(m_mat && "ThreadedSparseProduct: matrix not set; call analyzePattern() first");
    if (Adjoint) {
      eigen_assert(x.size() == m_mat->rows());
      eigen_assert(y.size() == m_mat->cols());
    } else {
      eigen_assert(x.size() == m_mat->cols());
      eigen_assert(y.size() == m_mat->rows());
    }
    // The kernel reads x and writes y concurrently across threads; aliasing
    // (x and y overlap) would mis-compute because some x[k] reads see y
    // writes from the same SpMV call. Cheap address-range check; cast to
    // uintptr_t because comparing pointers from unrelated allocations is
    // technically UB in C++.
    eigen_assert((x.size() == 0 || y.size() == 0 || std::uintptr_t(x.data() + x.size()) <= std::uintptr_t(y.data()) ||
                  std::uintptr_t(y.data() + y.size()) <= std::uintptr_t(x.data())) &&
                 "ThreadedSparseProduct: x and y must not overlap");
    // Decide which storage to read from.
    //   Forward kernel iterates a RowMajor view of A; adjoint kernel iterates
    //   a ColMajor view of A. If A already has the required order, read it
    //   directly; otherwise build/use the mirror.
    constexpr bool NeedRowMajorView = !Adjoint;
    constexpr bool NativeIsRowMajor = IsRowMajor;
    constexpr bool UseMirror = NeedRowMajorView != NativeIsRowMajor;

    if (UseMirror) {
      const MirrorType& m = ensure_mirror();
      run<Adjoint, Overwrite>(m.valuePtr(), m.innerIndexPtr(), m.outerIndexPtr(), m.innerNonZeroPtr(),
                              /*outerSize=*/m.outerSize(), m_adj_partition, x, y, alpha);
    } else {
      run<Adjoint, Overwrite>(m_mat->valuePtr(), m_mat->innerIndexPtr(), m_mat->outerIndexPtr(),
                              m_mat->innerNonZeroPtr(),
                              /*outerSize=*/m_mat->outerSize(), m_native_partition, x, y, alpha);
    }
  }

  // Serial-fallback threshold; matches the same constant in SparseDenseProduct.h.
  static constexpr Index kThreadingThreshold = 20000;

  template <bool Conjugate, bool Overwrite>
  void run(const Scalar* vals, const StorageIndex* inner, const StorageIndex* outer, const StorageIndex* innerNnz,
           Index outerSize, const std::vector<Index>& part, const ConstVectorRef& x, MutableVectorRef& y,
           const Scalar& alpha) const {
    // OpenMP path doesn't use the cached partition (dynamic scheduling balances
    // itself), so derive T fresh from the current Eigen::nbThreads() /
    // OMP_NUM_THREADS at apply time -- otherwise `setNbThreads()` after
    // analyzePattern() would be silently ignored. The ThreadPool path is
    // bound to the partition built for a specific T at analyzePattern() time.
#ifdef EIGEN_HAS_OPENMP
    const int T = target_thread_count();
#else
    const int T = static_cast<int>(part.size()) - 1;
#endif
    const Index total_nnz = m_mat->nonZeros();
    // Ref construction already enforced unit inner stride for x/y.
    const Scalar* xp = x.data();
    Scalar* yp = y.data();
    if (T <= 1 || total_nnz < kThreadingThreshold) {
      internal::run_dot_chunk<Conjugate, Overwrite, Scalar, StorageIndex>(vals, inner, outer, innerNnz, xp, yp, 0,
                                                                          outerSize, alpha);
      return;
    }

#ifdef EIGEN_HAS_OPENMP
    // Prefer OpenMP for dispatch when available. Use dynamic scheduling
    // over rows with chunks sized so the OMP runtime gets ~4*T chunks to
    // distribute -- enough granularity to absorb residual nnz imbalance
    // without blowing dispatch overhead.
    const Index chunk = numext::maxi<Index>(Index(1), (outerSize + Index(T) * 4 - 1) / (Index(T) * 4));
#pragma omp parallel for schedule(dynamic, chunk) num_threads(T)
    for (Index i = 0; i < outerSize; ++i) {
      internal::run_dot_row<Conjugate, Overwrite, Scalar, StorageIndex>(vals, inner, outer, innerNnz, xp, yp, i, alpha);
    }
#else
    // ThreadPool path: enqueue T-1 worker tasks, run partition 0 on this
    // thread, then wait on the barrier. Avoids ForkJoin's log(T)-hop critical
    // path before the last leaf starts.
    Barrier barrier(static_cast<unsigned>(T));
    ThreadPool* p = pool();
    for (int t = 1; t < T; ++t) {
      const Index lo = part[t], hi = part[t + 1];
      if (lo == hi) {
        barrier.Notify();
        continue;
      }
      p->Schedule([=, &barrier]() {
        internal::run_dot_chunk<Conjugate, Overwrite, Scalar, StorageIndex>(vals, inner, outer, innerNnz, xp, yp, lo,
                                                                            hi, alpha);
        barrier.Notify();
      });
    }
    internal::run_dot_chunk<Conjugate, Overwrite, Scalar, StorageIndex>(vals, inner, outer, innerNnz, xp, yp, part[0],
                                                                        part[1], alpha);
    barrier.Notify();
    barrier.Wait();
#endif
  }

  // Lazy mirror construction. Double-checked atomic init avoids per-call
  // std::call_once cost (which historically varies across libstdc++ versions).
  const MirrorType& ensure_mirror() const {
    MirrorType* m = m_mirror.load(std::memory_order_acquire);
    if (m) return *m;
    std::lock_guard<std::mutex> lock(m_mirror_init_mu);
    m = m_mirror.load(std::memory_order_relaxed);
    if (!m) {
      // Same logical matrix, opposite storage order. Eigen detects the
      // storage-order mismatch in assignment and reorganizes the data
      // (a "structural transpose"); the logical matrix is preserved, no
      // conjugation. The kernel applies conj at use when needed.
      // Hold the mirror in a unique_ptr while building: the assignment,
      // makeCompressed(), and build_partition() can all throw (bad_alloc), and a
      // raw owning pointer would leak. Hand ownership to the atomic only once the
      // mirror is fully built; on throw the buffer is freed and m_mirror stays
      // null so a later call rebuilds cleanly.
      std::unique_ptr<MirrorType> built(new MirrorType(m_mat->rows(), m_mat->cols()));
      *built = *m_mat;
      built->makeCompressed();
      build_partition(built->outerIndexPtr(), built->outerSize(), built->nonZeros(), m_adj_partition);
      m = built.get();
      m_mirror.store(built.release(), std::memory_order_release);
    }
    return *m;
  }

  // Thread count target.  Under OpenMP, respect Eigen::setNbThreads() /
  // OMP_NUM_THREADS instead of the lazy default ThreadPool's size (which
  // also avoids constructing that pool when OMP does the dispatch).
  int target_thread_count() const {
#ifdef EIGEN_HAS_OPENMP
    return numext::maxi(1, Eigen::nbThreads());
#else
    return pool()->NumThreads();
#endif
  }

  // nnz-balanced partition of A's outer dim into T contiguous chunks, with
  // a guard against pathological skew: if more than half the partitions end
  // up empty (one hub outer holding most of the nnz), fall back to a single
  // serial chunk so the inner kernel doesn't pay parallel dispatch for no
  // parallel work.
  void build_partition(const StorageIndex* outer, Index outerSize, Index nnz, std::vector<Index>& part) const {
    const int T = target_thread_count();
    internal::compute_nnz_balanced_partition(outer, outerSize, nnz, T, part);
    int non_empty = 0;
    for (std::size_t t = 0; t + 1 < part.size(); ++t)
      if (part[t + 1] > part[t]) ++non_empty;
    if (non_empty * 2 < T) {
      part.assign(2, 0);
      part[1] = outerSize;
    }
  }

 private:
  const SparseMatrixType* m_mat = nullptr;
  ThreadPool* m_pool = nullptr;

  // Partition of A's native outer dim by nnz balance. Used by the direction
  // whose kernel matches A's storage order.
  std::vector<Index> m_native_partition;

  // Lazy adjoint-direction mirror in the opposite storage order, plus its
  // own nnz-balanced partition. Constructed on first use.
  mutable std::atomic<MirrorType*> m_mirror{nullptr};
  mutable std::mutex m_mirror_init_mu;
  mutable std::vector<Index> m_adj_partition;
};

}  // namespace Eigen

#endif  // EIGEN_THREADED_SPARSE_PRODUCT_H
