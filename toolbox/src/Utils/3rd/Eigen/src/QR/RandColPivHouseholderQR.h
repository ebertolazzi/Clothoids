// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#ifndef EIGEN_RANDCOLPIVOTINGHOUSEHOLDERQR_H
#define EIGEN_RANDCOLPIVOTINGHOUSEHOLDERQR_H

#include <random>

// IWYU pragma: private
#include "./InternalHeaderCheck.h"

namespace Eigen {

namespace internal {

template <typename MatrixType_, typename PermutationIndex_>
struct traits<RandColPivHouseholderQR<MatrixType_, PermutationIndex_>> : traits<MatrixType_> {
  typedef MatrixXpr XprKind;
  typedef SolverStorage StorageKind;
  typedef PermutationIndex_ PermutationIndex;
  enum { Flags = 0 };
};

// Fill `mat` with iid samples drawn from a real standard normal distribution.
template <typename Derived, typename Engine>
EIGEN_STRONG_INLINE std::enable_if_t<!NumTraits<typename Derived::Scalar>::IsComplex> fill_gaussian(
    MatrixBase<Derived>& mat, Engine& engine) {
  typedef typename Derived::Scalar Scalar;
  std::normal_distribution<Scalar> dist(Scalar(0), Scalar(1));
  for (Index j = 0; j < mat.cols(); ++j)
    for (Index i = 0; i < mat.rows(); ++i) mat.coeffRef(i, j) = dist(engine);
}

// Fill `mat` with iid samples drawn from a complex standard normal
// distribution (real and imaginary parts sampled independently).
template <typename Derived, typename Engine>
EIGEN_STRONG_INLINE std::enable_if_t<NumTraits<typename Derived::Scalar>::IsComplex> fill_gaussian(
    MatrixBase<Derived>& mat, Engine& engine) {
  typedef typename Derived::Scalar Scalar;
  typedef typename NumTraits<Scalar>::Real RealScalar;
  std::normal_distribution<RealScalar> dist(RealScalar(0), RealScalar(1));
  for (Index j = 0; j < mat.cols(); ++j)
    for (Index i = 0; i < mat.rows(); ++i) mat.coeffRef(i, j) = Scalar(dist(engine), dist(engine));
}

// Stable LAWN-176 column-norm downdate after a Householder reflector has
// been applied. `pivot_entry` is the entry of the just-pivoted row above
// column `j`; `tail_norm_fn` recomputes the trailing-column norm from
// scratch when the running estimate has lost too much precision (see
// http://www.netlib.org/lapack/lawnspdf/lawn176.pdf).
template <typename RealScalar, typename Scalar, typename RecomputeFn>
EIGEN_STRONG_INLINE void lawn176_norm_downdate(RealScalar& norm_updated, RealScalar& norm_direct, Scalar pivot_entry,
                                               RealScalar downdate_threshold, RecomputeFn&& tail_norm_fn) {
  using std::abs;
  if (numext::is_exactly_zero(norm_updated)) return;
  RealScalar t = abs(pivot_entry) / norm_updated;
  t = (RealScalar(1) + t) * (RealScalar(1) - t);
  if (t < RealScalar(0)) t = RealScalar(0);
  RealScalar t2 = t * numext::abs2<RealScalar>(norm_updated / norm_direct);
  if (t2 <= downdate_threshold) {
    norm_direct = tail_norm_fn();
    norm_updated = norm_direct;
  } else {
    norm_updated *= numext::sqrt(t);
  }
}

}  // end namespace internal

/** \ingroup QR_Module
 *
 * \class RandColPivHouseholderQR
 *
 * \brief Randomized blocked Householder rank-revealing QR with column pivoting
 *
 * \tparam MatrixType_ the type of the matrix being decomposed.
 * \tparam PermutationIndex_ the type of the permutation indices.
 *
 * Computes \f$ \mathbf{A} \mathbf{P} = \mathbf{Q} \mathbf{R} \f$ using the
 * BQRRP framework introduced by Melnichenko, Murray, Killian, Demmel,
 * Mahoney, Luszczek, and Gates, *Anatomy of High-Performance Column-Pivoted
 * QR Decomposition*, arXiv:2507.00976 (2025). BQRRP is itself a refinement
 * of the earlier randomized blocked QRCP schemes of Duersch and Gu (RQRCP,
 * SIAM J. Sci. Comput. 39(4):C263–C291, 2017, arXiv:1509.06820) and
 * Martinsson, Quintana-Ortí, Heavner, and van de Geijn (HQRRP, SIAM J. Sci.
 * Comput. 39(2):C96–C115, 2017, arXiv:1512.02671).
 *
 * The classical Businger–Golub pivot scan is replaced by selecting blocks
 * of \c b pivots from a small Gaussian sketch \f$ \mathbf{Y} = \mathbf{G}
 * \mathbf{A} \f$ with \f$ \mathbf{G} \in \mathbb{R}^{b \times m} \f$
 * having i.i.d.\ standard normal entries. Following the BQRRP paper, pivot
 * decisions on the sketch are produced by a partial-pivoted LU on the
 * transposed sketch (cheap and robust), the panel itself is factored with
 * unpivoted blocked Householder QR, the trailing block is updated through
 * the compact-WY apply, and the sketch is downdated by the closed-form
 * Duersch–Gu update. After each block step the asymptotic flop count
 * matches that of an unpivoted blocked Householder QR
 * (\f$ 2mn^2 - \tfrac{2}{3}n^3 \f$ for \f$ m \ge n \f$); almost all
 * computation is BLAS-3, which lifts the BLAS-2 ceiling that limits
 * ColPivHouseholderQR.
 *
 * The pivot-quality tradeoff is empirically minor: on the matrices in the
 * cited papers, the rank-revealing behavior is comparable to classical
 * column pivoting (LAPACK \c geqp3).
 *
 * The block size \c b can be configured via setBlockSize(); a value of
 * \c 0 (the default) leaves the algorithm free to pick a size that scales
 * with the input. The seed for the internal RNG can be fixed with
 * setSeed() for reproducible output.
 *
 * \note The \c setOversampling() / \c oversampling() accessors are
 * retained as documented no-ops for source compatibility with earlier
 * HQRRP-style versions of this class. The BQRRP framework with
 * partial-pivoted-LU pivot selection makes oversampling unnecessary
 * (see Section 2.1 of arXiv:2507.00976).
 *
 * This class supports the \link InplaceDecomposition inplace decomposition
 * \endlink mechanism. The same public API as ColPivHouseholderQR is
 * provided.
 *
 * \sa ColPivHouseholderQR, MatrixBase::randColPivHouseholderQr()
 */
template <typename MatrixType_, typename PermutationIndex_>
class RandColPivHouseholderQR : public SolverBase<RandColPivHouseholderQR<MatrixType_, PermutationIndex_>>,
                                public RankRevealingBase<RandColPivHouseholderQR<MatrixType_, PermutationIndex_>> {
 public:
  typedef MatrixType_ MatrixType;
  typedef SolverBase<RandColPivHouseholderQR> Base;
  typedef RankRevealingBase<RandColPivHouseholderQR> RankRevealingBase_;
  friend class SolverBase<RandColPivHouseholderQR>;
  friend class RankRevealingBase<RandColPivHouseholderQR>;
  using RankRevealingBase_::dimensionOfKernel;
  using RankRevealingBase_::isInjective;
  using RankRevealingBase_::isInvertible;
  using RankRevealingBase_::isSurjective;
  using RankRevealingBase_::maxPivot;
  using RankRevealingBase_::nonzeroPivots;
  using RankRevealingBase_::rank;
  using RankRevealingBase_::setThreshold;
  using RankRevealingBase_::threshold;
  typedef PermutationIndex_ PermutationIndex;
  EIGEN_GENERIC_PUBLIC_INTERFACE(RandColPivHouseholderQR)

  enum {
    MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
  };
  typedef typename internal::plain_diag_type<MatrixType>::type HCoeffsType;
  typedef PermutationMatrix<ColsAtCompileTime, MaxColsAtCompileTime, PermutationIndex> PermutationType;
  typedef typename internal::plain_row_type<MatrixType, PermutationIndex>::type IntRowVectorType;
  typedef typename internal::plain_row_type<MatrixType>::type RowVectorType;
  typedef typename internal::plain_row_type<MatrixType, RealScalar>::type RealRowVectorType;
  typedef HouseholderSequence<MatrixType, internal::remove_all_t<typename HCoeffsType::ConjugateReturnType>>
      HouseholderSequenceType;
  typedef typename MatrixType::PlainObject PlainObject;

 private:
  // Default `m_blockSize == 0` means: let computeInPlace pick a size
  // tuned to the input dimensions and the host cache hierarchy.
  static constexpr Index kAutoBlockSize = 0;

  // SIMD-kernel efficiency floor; smaller blocks degrade the panel-QR
  // and trailing-GEMM kernels without speeding up the sketch.
  static constexpr Index kKernelFloor = 48;
  // Hard upper bound on the auto block size.
  static constexpr Index kBlockCeiling = 1024;

  // Minimum min(rows, cols) at which auto-block enters the blocked path.
  // Below this, the per-compute sketch overhead (G fill, G*A, per-panel
  // LU on the sketch) is not amortized by the trailing-GEMM savings vs.
  // classical Businger-Golub pivoting. Empirically calibrated (Apple M4,
  // double): n ~= 128 loses by ~1.5x, n >= 256 wins. Set above the
  // kernel-efficiency floor so that `can_block` (size >= 2*b) cannot
  // fire below it. A user who calls `setBlockSize(b > 0)` explicitly
  // bypasses this gate.
  static constexpr Index kAutoBlockedPathMinSize = 192;

  // L2 cache size in bytes, from Eigen's shared cache-sizes singleton
  // (which handles cross-platform detection plus user overrides via
  // setCpuCacheSizes).
  static Index defaultL2Bytes() {
    std::ptrdiff_t l1, l2, l3;
    internal::manage_caching_sizes(GetAction, &l1, &l2, &l3);
    return Index(l2);
  }

  // Roofline-derived auto block size for `m_blockSize == kAutoBlockSize`.
  // `size = min(rows, cols)` is passed in to avoid recomputing it.
  //
  // The dominant per-iteration cost is the trailing-update GEMM
  // (compact-WY apply: C := C - V*T*V^T*C). With panel V resident in
  // cache, the arithmetic intensity is
  //    AI = 2*b / sizeof(Scalar)     [FLOPS / byte]
  // which crosses the roofline ceiling AI_crit = peak_FLOPS / cache_BW
  // at b roughly 2-3 on current hardware. That bound is dominated by the
  // SIMD-kernel efficiency floor (~48); the roofline therefore gives no
  // useful lower bound on b.
  //
  // The binding constraint is the cache-fit invariant: V must stay
  // resident across the trailing sweep, otherwise the streamed re-reads
  // drop us off the compute roof onto the DRAM-bandwidth side. That
  // yields the upper bound
  //    b * rows * sizeof(Scalar) <= L2 / 2     (half of L2; remainder
  //                                              for the streamed C tile)
  // BQRRP's recommended b ~= n/32 (paper Section 3.2) sits inside this
  // window for matrices up to roughly 16 * L2 / (rows * sizeof(Scalar));
  // beyond that, the cache bound binds and clamps b.
  static Index computeAutoBlockSize(Index rows, Index size) {
    // Degenerate empty inputs (0 rows or 0 cols): there is no work, and
    // the cache-bound calculation below would divide by `rows`. Return 0
    // so the caller's can_block check (b > 0) routes to the unblocked
    // path, which handles empty matrices correctly.
    if (rows <= 0 || size <= 0) return Index(0);
    const Index scalar_bytes = Index(sizeof(Scalar));
    const Index b_cache = numext::maxi(Index(1), (defaultL2Bytes() / Index(2)) / (rows * scalar_bytes));
    const Index b_bqrrp = size / Index(32);

    Index b = numext::maxi(Index(kKernelFloor), b_bqrrp);
    b = numext::mini(b, numext::maxi(Index(kKernelFloor), b_cache));
    b = numext::mini(b, Index(kBlockCeiling));
    return b;
  }

  void init(Index rows, Index cols) {
    Index diag = numext::mini(rows, cols);
    m_hCoeffs.resize(diag);
    m_colsPermutation.resize(cols);
    m_temp.resize(cols);
    m_isInitialized = false;
  }

 public:
  /** \brief Default constructor. */
  RandColPivHouseholderQR() = default;

  /** \brief Constructor with memory preallocation. */
  RandColPivHouseholderQR(Index rows, Index cols) : m_qr(rows, cols) { init(rows, cols); }

  /** \brief Constructs and computes a QR factorization from \a matrix. */
  template <typename InputType>
  explicit RandColPivHouseholderQR(const EigenBase<InputType>& matrix) : m_qr(matrix.rows(), matrix.cols()) {
    init(matrix.rows(), matrix.cols());
    compute(matrix.derived());
  }

  /** \brief Inplace constructor: takes a Ref and decomposes in place. */
  template <typename InputType>
  explicit RandColPivHouseholderQR(EigenBase<InputType>& matrix) : m_qr(matrix.derived()) {
    init(matrix.rows(), matrix.cols());
    computeInPlace();
  }

#ifdef EIGEN_PARSED_BY_DOXYGEN
  template <typename Rhs>
  inline Solve<RandColPivHouseholderQR, Rhs> solve(const MatrixBase<Rhs>& b) const;
#endif

  HouseholderSequenceType householderQ() const;
  HouseholderSequenceType matrixQ() const { return householderQ(); }

  const MatrixType& matrixQR() const {
    eigen_assert(m_isInitialized && "RandColPivHouseholderQR is not initialized.");
    return m_qr;
  }

  const MatrixType& matrixR() const {
    eigen_assert(m_isInitialized && "RandColPivHouseholderQR is not initialized.");
    return m_qr;
  }

  template <typename InputType>
  RandColPivHouseholderQR& compute(const EigenBase<InputType>& matrix);

  const PermutationType& colsPermutation() const {
    eigen_assert(m_isInitialized && "RandColPivHouseholderQR is not initialized.");
    return m_colsPermutation;
  }

  typename MatrixType::Scalar determinant() const;
  typename MatrixType::RealScalar absDeterminant() const;
  typename MatrixType::RealScalar logAbsDeterminant() const;
  typename MatrixType::Scalar signDeterminant() const;

  RealScalar pivotCoeff(Index i) const {
    using std::abs;
    return abs(m_qr.coeff(i, i));
  }

  inline Inverse<RandColPivHouseholderQR> inverse() const {
    eigen_assert(m_isInitialized && "RandColPivHouseholderQR is not initialized.");
    return Inverse<RandColPivHouseholderQR>(*this);
  }

  inline Index rows() const { return m_qr.rows(); }
  inline Index cols() const { return m_qr.cols(); }

  const HCoeffsType& hCoeffs() const { return m_hCoeffs; }

  ComputationInfo info() const {
    eigen_assert(m_isInitialized && "Decomposition is not initialized.");
    return Success;
  }

  /** \brief Sets the panel block size \c b.
   *
   * Larger \c b increases the proportion of work performed in level-3
   * BLAS but raises the per-iteration overhead. Pass \c b = 0 (the
   * default) to let the algorithm pick a size that scales with the
   * input dimensions; this matches the BQRRP paper's recommendation of
   * roughly \c n/32 for large square inputs.
   *
   * For small matrices (where \c min(m,n) < 2b) this class transparently
   * falls back to the unblocked column-pivoted scalar loop.
   */
  RandColPivHouseholderQR& setBlockSize(Index b) {
    eigen_assert(b >= 0 && "Block size must be non-negative.");
    m_blockSize = b;
    return *this;
  }

  /** \brief Sets the oversampling parameter \c p.
   *
   * \deprecated With the BQRRP-style LU-on-transposed-sketch pivot
   * selection used by this class, oversampling provides no benefit (see
   * Section 2.1 of arXiv:2507.00976: the first \c b pivots are
   * independent of the oversampling factor). Retained for source
   * compatibility; the value is ignored.
   */
  RandColPivHouseholderQR& setOversampling(Index /*p*/) { return *this; }

  /** \brief Fixes the seed of the internal RNG for reproducible factorization.
   *
   * If never called, each call to compute() draws a fresh seed from
   * \c std::random_device.
   */
  RandColPivHouseholderQR& setSeed(uint64_t seed) {
    m_seed = seed;
    m_seedSet = true;
    return *this;
  }

  /** \brief Returns the user-set panel block size, or \c 0 if the
   * algorithm should pick automatically. The actual block size used
   * during compute() is not exposed.
   */
  Index blockSize() const { return m_blockSize; }

  /** \deprecated Returns 0; oversampling is no longer used. \sa setOversampling() */
  Index oversampling() const { return Index(0); }

#ifndef EIGEN_PARSED_BY_DOXYGEN
  template <typename RhsType, typename DstType>
  void _solve_impl(const RhsType& rhs, DstType& dst) const;

  template <bool Conjugate, typename RhsType, typename DstType>
  void _solve_impl_transposed(const RhsType& rhs, DstType& dst) const;
#endif

 protected:
  friend class internal::CompleteOrthogonalDecompositionImpl<MatrixType, PermutationIndex, RandColPivHouseholderQR>;

  EIGEN_STATIC_ASSERT_NON_INTEGER(Scalar)

  void computeInPlace();

  // Unblocked column-pivoted Householder QR on rows [row0, m) and columns
  // [col0, col0 + ncols), using full-column swaps. Updates m_hCoeffs in
  // [col0, col0 + ncols), m_colsPermutation, and the maxpivot tracker.
  // Returns the number of column transpositions applied. Used both as
  // the small-matrix fallback path and as the rank-deficient tail of the
  // blocked path.
  Index unblocked_pivoted_qr(Index row0, Index col0, Index ncols, RealScalar threshold_helper);

  // Scan the diagonal of m_qr and set m_nonzero_pivots / m_maxpivot using
  // the LAWN-176 threshold (same scale as ColPivHouseholderQR).
  void finalize_rank(RealScalar threshold_helper);

  MatrixType m_qr;
  HCoeffsType m_hCoeffs;
  PermutationType m_colsPermutation;
  RowVectorType m_temp;
  Index m_blockSize = kAutoBlockSize;
  uint64_t m_seed = 0;
  bool m_seedSet = false;
  bool m_isInitialized = false;
  Index m_det_p = 1;
};

template <typename MatrixType, typename PermutationIndex>
typename MatrixType::Scalar RandColPivHouseholderQR<MatrixType, PermutationIndex>::determinant() const {
  eigen_assert(m_isInitialized && "RandColPivHouseholderQR is not initialized.");
  eigen_assert(m_qr.rows() == m_qr.cols() && "You can't take the determinant of a non-square matrix!");
  Scalar detQ;
  internal::householder_determinant<HCoeffsType, Scalar, NumTraits<Scalar>::IsComplex>::run(m_hCoeffs, detQ);
  return isInjective() ? (detQ * Scalar(m_det_p)) * m_qr.diagonal().prod() : Scalar(0);
}

template <typename MatrixType, typename PermutationIndex>
typename MatrixType::RealScalar RandColPivHouseholderQR<MatrixType, PermutationIndex>::absDeterminant() const {
  using std::abs;
  eigen_assert(m_isInitialized && "RandColPivHouseholderQR is not initialized.");
  eigen_assert(m_qr.rows() == m_qr.cols() && "You can't take the determinant of a non-square matrix!");
  return isInjective() ? abs(m_qr.diagonal().prod()) : RealScalar(0);
}

template <typename MatrixType, typename PermutationIndex>
typename MatrixType::RealScalar RandColPivHouseholderQR<MatrixType, PermutationIndex>::logAbsDeterminant() const {
  eigen_assert(m_isInitialized && "RandColPivHouseholderQR is not initialized.");
  eigen_assert(m_qr.rows() == m_qr.cols() && "You can't take the determinant of a non-square matrix!");
  return isInjective() ? m_qr.diagonal().cwiseAbs().array().log().sum() : -NumTraits<RealScalar>::infinity();
}

template <typename MatrixType, typename PermutationIndex>
typename MatrixType::Scalar RandColPivHouseholderQR<MatrixType, PermutationIndex>::signDeterminant() const {
  eigen_assert(m_isInitialized && "RandColPivHouseholderQR is not initialized.");
  eigen_assert(m_qr.rows() == m_qr.cols() && "You can't take the determinant of a non-square matrix!");
  Scalar detQ;
  internal::householder_determinant<HCoeffsType, Scalar, NumTraits<Scalar>::IsComplex>::run(m_hCoeffs, detQ);
  return isInjective() ? (detQ * Scalar(m_det_p)) * m_qr.diagonal().array().sign().prod() : Scalar(0);
}

template <typename MatrixType, typename PermutationIndex>
template <typename InputType>
RandColPivHouseholderQR<MatrixType, PermutationIndex>& RandColPivHouseholderQR<MatrixType, PermutationIndex>::compute(
    const EigenBase<InputType>& matrix) {
  m_qr = matrix.derived();
  computeInPlace();
  return *this;
}

template <typename MatrixType, typename PermutationIndex>
Index RandColPivHouseholderQR<MatrixType, PermutationIndex>::unblocked_pivoted_qr(Index row0, Index col0, Index ncols,
                                                                                  RealScalar threshold_helper) {
  using std::abs;
  const Index rows = m_qr.rows();
  const Index col_end = col0 + ncols;
  const Index sub_rows = rows - row0;
  Index num_transpositions = 0;

  Matrix<RealScalar, 1, Dynamic> norms_direct = m_qr.block(row0, col0, sub_rows, ncols).colwise().norm();
  Matrix<RealScalar, 1, Dynamic> norms_updated = norms_direct;
  const RealScalar downdate_threshold = numext::sqrt(NumTraits<RealScalar>::epsilon());

  const Index size = (std::min)(sub_rows, ncols);
  for (Index k = 0; k < size; ++k) {
    Index biggest;
    RealScalar biggest_sq = numext::abs2(norms_updated.tail(ncols - k).maxCoeff(&biggest));
    biggest += k;

    if (this->m_nonzero_pivots == m_qr.diagonalSize() && biggest_sq < threshold_helper * RealScalar(rows - row0 - k))
      this->m_nonzero_pivots = col0 + k;

    if (k != biggest) {
      m_qr.col(col0 + k).swap(m_qr.col(col0 + biggest));
      std::swap(norms_updated.coeffRef(k), norms_updated.coeffRef(biggest));
      std::swap(norms_direct.coeffRef(k), norms_direct.coeffRef(biggest));
      m_colsPermutation.applyTranspositionOnTheRight(col0 + k, col0 + biggest);
      ++num_transpositions;
    }

    RealScalar beta;
    m_qr.col(col0 + k).tail(sub_rows - k).makeHouseholderInPlace(m_hCoeffs.coeffRef(col0 + k), beta);
    m_qr.coeffRef(row0 + k, col0 + k) = beta;
    if (abs(beta) > this->m_maxpivot) this->m_maxpivot = abs(beta);

    // Apply the reflector only to the remaining panel columns. Columns at
    // indices >= col_end are owned by the outer caller and updated via a
    // BLAS-3 trailing update; in the tail/fallback path col_end equals
    // m_qr.cols(), so this naturally degrades to a full apply.
    const Index trail_cols = col_end - col0 - k - 1;
    if (trail_cols > 0) {
      m_qr.block(row0 + k, col0 + k + 1, sub_rows - k, trail_cols)
          .applyHouseholderOnTheLeft(m_qr.col(col0 + k).tail(sub_rows - k - 1), m_hCoeffs.coeff(col0 + k),
                                     &m_temp.coeffRef(col0 + k + 1));
    }

    for (Index j = k + 1; j < ncols; ++j) {
      internal::lawn176_norm_downdate(norms_updated.coeffRef(j), norms_direct.coeffRef(j),
                                      m_qr.coeff(row0 + k, col0 + j), downdate_threshold,
                                      [&] { return m_qr.col(col0 + j).tail(sub_rows - k - 1).norm(); });
    }
  }
  return num_transpositions;
}

template <typename MatrixType, typename PermutationIndex>
void RandColPivHouseholderQR<MatrixType, PermutationIndex>::finalize_rank(RealScalar threshold_helper) {
  using std::abs;
  const Index rows = m_qr.rows();
  const Index size = (std::min)(rows, m_qr.cols());
  // m_maxpivot is already up-to-date: the blocked path tracks it per
  // panel, and unblocked_pivoted_qr tracks it per step.
  // If neither the blocked path nor the unblocked tail tightened the
  // rank cap, fall back to the LAWN-176 first-below-threshold scan.
  // This handles the full-rank-blocked-path case where every panel
  // was full rank and the only "small" diagonal entries are in the
  // unblocked tail.
  if (this->m_nonzero_pivots == size) {
    for (Index i = 0; i < size; ++i) {
      RealScalar a = abs(m_qr.coeff(i, i));
      if (numext::abs2(a) < threshold_helper * RealScalar(rows - i)) {
        this->m_nonzero_pivots = i;
        break;
      }
    }
  }
}

template <typename MatrixType, typename PermutationIndex>
void RandColPivHouseholderQR<MatrixType, PermutationIndex>::computeInPlace() {
  eigen_assert(m_qr.cols() <= NumTraits<PermutationIndex>::highest());

  const Index rows = m_qr.rows();
  const Index cols = m_qr.cols();
  const Index size = (std::min)(rows, cols);

  m_hCoeffs.resize(size);
  m_temp.resize(cols);
  m_colsPermutation.resize(cols);
  m_colsPermutation.setIdentity();
  this->m_nonzero_pivots = size;
  this->m_maxpivot = RealScalar(0);
  Index num_transpositions = 0;

  // Global rank-detection scale, derived from the original column norms
  // exactly as ColPivHouseholderQR does (LAWN-176). All panel calls share
  // this threshold so rank revelation is consistent across blocks.
  RealScalar max_initial_norm = RealScalar(0);
  for (Index j = 0; j < cols; ++j) max_initial_norm = numext::maxi(max_initial_norm, m_qr.col(j).norm());
  const RealScalar threshold_helper =
      cols == 0 ? RealScalar(0)
                : numext::abs2<RealScalar>(max_initial_norm * NumTraits<RealScalar>::epsilon()) / RealScalar(rows);

  const bool auto_block = (m_blockSize == kAutoBlockSize);
  const Index requested_b = auto_block ? computeAutoBlockSize(rows, size) : m_blockSize;
  const Index b = (std::min)(requested_b, size);

  // Gates for entering the blocked path:
  //   - `b > 0` and `rows > b` for the trailing update to have work,
  //   - `size >= 2*b` so at least two full blocks fit (one block + tail
  //     is handled more efficiently by the unblocked path),
  //   - in the auto-block path, an additional size threshold so the
  //     sketch overhead is amortized by trailing-GEMM savings (see
  //     kAutoBlockedPathMinSize). Users who pin a block size via
  //     setBlockSize() bypass this — useful for tests that need to
  //     exercise the blocked path on small inputs.
  const bool can_block = (b > 0) && (size >= 2 * b) && (rows > b) && (!auto_block || size >= kAutoBlockedPathMinSize);
  if (!can_block) {
    num_transpositions += unblocked_pivoted_qr(0, 0, cols, threshold_helper);
    m_det_p = (num_transpositions % 2) ? -1 : 1;
    finalize_rank(threshold_helper);
    m_isInitialized = true;
    return;
  }

  typedef Matrix<Scalar, Dynamic, Dynamic, ColMajor> WorkMatrix;
  typedef Matrix<Scalar, Dynamic, 1> WorkVector;
  // Dynamic-sized Ref types so that a fixed-size MatrixType (e.g.
  // Matrix<float, 8, 10>) and its run-time-sized blocks can both be
  // wrapped without tripping Ref's compile-time-size check.
  typedef Ref<WorkMatrix, 0, OuterStride<>> WorkMatrixRef;
  typedef Ref<WorkVector> HCoeffsRef;
  typedef Transpositions<Dynamic, Dynamic, PermutationIndex> IpivType;

  // Hoisted workspaces — allocated once per compute() call. Worst-case
  // sizes are computed from the first iteration (n_remain = cols,
  // trail_cols = cols - b); subsequent iterations use shrinking
  // sub-blocks via leftCols/topRows.
  WorkMatrix Y(b, cols);
  WorkMatrix YT(cols, b);
  WorkMatrix sketch_R(b, b);
  WorkMatrix tmp(b, cols - b);
  IpivType ipiv(b);
  WorkVector y_hcoeffs(b);
  WorkVector y_temp(cols);

  // Initialize the sketch Y = G * A. We don't keep G around; the
  // Duersch-Gu update (step 24) maintains Y deterministically from
  // here on.
  {
    WorkMatrix G(b, rows);
    uint64_t seed = m_seed;
    if (!m_seedSet) {
      std::random_device rd;
      seed = (uint64_t(rd()) << 32) | uint64_t(rd());
    }
    std::mt19937_64 engine(seed);
    internal::fill_gaussian(G, engine);
    Y.noalias() = G * m_qr;
  }

  Index k = 0;
  bool blocked_terminated_early = false;
  while (k + b <= size && cols - k > b) {
    const Index n_remain = cols - k;        // columns from k to end
    const Index trail_cols = n_remain - b;  // columns to the right of the panel
    const Index sub_rows = rows - k;        // rows from k to end

    // === Step 7 (Algorithm 2): "crazy man's QRCP" on the sketch.
    // Form Y_T = Y(:, k:)^T and run partial-pivoted LU on it. The IPIV
    // vector tells us which b columns of Y(:, k:) (and thus which
    // columns of m_qr.middleCols(k, n_remain)) to bring to the front.
    auto Y_curr = Y.middleCols(k, n_remain);
    auto Y_T = YT.topRows(n_remain);
    Y_T.noalias() = Y_curr.transpose();

    typename IpivType::StorageIndex nb_lu_transp = 0;
    internal::partial_lu_inplace(Y_T, ipiv, nb_lu_transp);

    // === Steps 9-11: apply IPIV to columns of m_qr and Y starting at
    // absolute column k (LAPACK-style serial swaps; equivalent to
    // applying Algorithm 4 of the BQRRP paper to the converted pivot
    // vector).
    for (Index i = 0; i < b; ++i) {
      Index dst = static_cast<Index>(ipiv.coeff(i));
      if (dst != i) {
        m_qr.col(k + i).swap(m_qr.col(k + dst));
        Y.col(k + i).swap(Y.col(k + dst));
        m_colsPermutation.applyTranspositionOnTheRight(k + i, k + dst);
        ++num_transpositions;
      }
    }

    // Materialize R_sk in the upper triangle of Y(:, k:) via unpivoted
    // blocked Householder QR. We discard the resulting Householder
    // vectors and tau; only R_sk feeds the step-24 sketch update.
    // Wrap in `Ref` so householder_qr_inplace_blocked sees a flat
    // MatrixQR — its internal `Block<MatrixQR, ...>` typedef does not
    // compose with Block-of-Block expressions.
    {
      WorkMatrixRef Y_curr_ref(Y_curr);
      Ref<WorkVector> y_hc_ref(y_hcoeffs);
      internal::householder_qr_inplace_blocked<WorkMatrixRef, Ref<WorkVector>>::run(
          Y_curr_ref, y_hc_ref, /*maxBlockSize=*/(std::min)(b, Index(48)), y_temp.data());
    }

    // Save R_sk_11 = upper-triangular portion of Y(:, k:k+b) before the
    // panel QR overwrites the corresponding columns of m_qr (via the
    // trailing update); we need it for step 24.
    sketch_R = Y.block(0, k, b, b);

    // === Step 12: tall unpivoted Householder QR on the panel.
    auto panel = m_qr.block(k, k, sub_rows, b);
    auto hCoeffsSegment = m_hCoeffs.segment(k, b);
    {
      WorkMatrixRef panel_ref(panel);
      HCoeffsRef hc_ref(hCoeffsSegment);
      internal::householder_qr_inplace_blocked<WorkMatrixRef, HCoeffsRef>::run(
          panel_ref, hc_ref, /*maxBlockSize=*/(std::min)(b, Index(48)), m_temp.data());
    }

    this->m_maxpivot = (std::max)(this->m_maxpivot, m_qr.diagonal().segment(k, b).cwiseAbs().maxCoeff());

    // Detect rank deficiency by scanning the panel's diagonal against
    // a relative threshold (max pivot seen so far times the default
    // RankRevealingBase tolerance). The LAWN-176 threshold used by
    // ColPivHouseholderQR assumes a strictly monotonic diagonal —
    // BQRRP's unpivoted panel QR breaks that assumption for matrices
    // with closely-spaced singular values (e.g. partial isometries),
    // so we use the same relative tolerance that RankRevealingBase::
    // rank() applies to the final diagonal.
    Index panel_rank = b;
    {
      using std::abs;
      const RealScalar relative_cutoff = this->m_maxpivot * NumTraits<RealScalar>::epsilon() * RealScalar(4 * size);
      for (Index i = 0; i < b; ++i) {
        RealScalar a = abs(m_qr.coeff(k + i, k + i));
        if (a <= relative_cutoff) {
          panel_rank = i;
          break;
        }
      }
    }

    // === Step 17: apply Q^H from the panel to the trailing block.
    // The loop guard `cols - k > b` guarantees `trail_cols > 0`. We
    // pass the un-conjugated coeffs because apply_block_householder_on_the_left
    // conjugates internally when forward=false (matching how HouseholderQR
    // drives this same routine).
    auto trailing = m_qr.block(k, k + b, sub_rows, trail_cols);
    internal::apply_block_householder_on_the_left(trailing, panel, hCoeffsSegment,
                                                  /*forward=*/false);

    if (panel_rank < b) {
      // Pin the rank cap and break out; the unblocked tail will
      // process the remaining columns. The tail's m_nonzero_pivots
      // write is gated on it still being unset (== diagonalSize),
      // so a tighter cap from here survives.
      if (this->m_nonzero_pivots == size) {
        this->m_nonzero_pivots = k + panel_rank;
      }
      k += b;
      blocked_terminated_early = true;
      break;
    }

    // === Step 24: Duersch-Gu sketch update.
    //   Y(:, k+b:) := R_sk_12 - R_sk_11 * R_11^{-1} * R_12
    // R_sk_12 currently sits in Y.middleCols(k+b, trail_cols) (the QR of
    // the sketch above placed it there). R_sk_11 we saved into sketch_R.
    // R_11 and R_12 are the corresponding blocks of m_qr after the
    // panel QR + trailing update.
    {
      auto tmp_block = tmp.leftCols(trail_cols);
      tmp_block = m_qr.block(k, k + b, b, trail_cols);
      m_qr.block(k, k, b, b).template triangularView<Upper>().solveInPlace(tmp_block);

      Y.middleCols(k + b, trail_cols).noalias() -= sketch_R.template triangularView<Upper>() * tmp_block;
    }

    k += b;
  }

  // Tail loop: any remaining columns (last partial block, or the rank-
  // deficient remainder after early termination) handled by the
  // unblocked pivoted QR. This routine also updates m_nonzero_pivots
  // along its diagonal scan, which finalize_rank below cross-checks.
  if (k < cols && (k < size || blocked_terminated_early)) {
    num_transpositions += unblocked_pivoted_qr(k, k, cols - k, threshold_helper);
  }

  m_det_p = (num_transpositions % 2) ? -1 : 1;

  // Final rank determination: scan the entire diagonal with the LAWN-176
  // threshold. This refines whatever the unblocked tail set, and covers
  // the full-rank blocked path where m_nonzero_pivots was never touched.
  finalize_rank(threshold_helper);

  m_isInitialized = true;
}

#ifndef EIGEN_PARSED_BY_DOXYGEN
template <typename MatrixType_, typename PermutationIndex_>
template <typename RhsType, typename DstType>
void RandColPivHouseholderQR<MatrixType_, PermutationIndex_>::_solve_impl(const RhsType& rhs, DstType& dst) const {
  const Index nonzero_pivots = nonzeroPivots();

  if (nonzero_pivots == 0) {
    dst.setZero();
    return;
  }

  typename RhsType::PlainObject c(rhs);

  c.applyOnTheLeft(householderQ().setLength(nonzero_pivots).adjoint());

  m_qr.topLeftCorner(nonzero_pivots, nonzero_pivots)
      .template triangularView<Upper>()
      .solveInPlace(c.topRows(nonzero_pivots));

  for (Index i = 0; i < nonzero_pivots; ++i) dst.row(m_colsPermutation.indices().coeff(i)) = c.row(i);
  for (Index i = nonzero_pivots; i < cols(); ++i) dst.row(m_colsPermutation.indices().coeff(i)).setZero();
}

template <typename MatrixType_, typename PermutationIndex_>
template <bool Conjugate, typename RhsType, typename DstType>
void RandColPivHouseholderQR<MatrixType_, PermutationIndex_>::_solve_impl_transposed(const RhsType& rhs,
                                                                                     DstType& dst) const {
  const Index nonzero_pivots = nonzeroPivots();

  if (nonzero_pivots == 0) {
    dst.setZero();
    return;
  }

  typename RhsType::PlainObject c(m_colsPermutation.transpose() * rhs);

  m_qr.topLeftCorner(nonzero_pivots, nonzero_pivots)
      .template triangularView<Upper>()
      .transpose()
      .template conjugateIf<Conjugate>()
      .solveInPlace(c.topRows(nonzero_pivots));

  dst.topRows(nonzero_pivots) = c.topRows(nonzero_pivots);
  dst.bottomRows(rows() - nonzero_pivots).setZero();

  dst.applyOnTheLeft(householderQ().setLength(nonzero_pivots).template conjugateIf<!Conjugate>());
}
#endif

namespace internal {

template <typename DstXprType, typename MatrixType, typename PermutationIndex>
struct Assignment<DstXprType, Inverse<RandColPivHouseholderQR<MatrixType, PermutationIndex>>,
                  internal::assign_op<typename DstXprType::Scalar,
                                      typename RandColPivHouseholderQR<MatrixType, PermutationIndex>::Scalar>,
                  Dense2Dense> {
  typedef RandColPivHouseholderQR<MatrixType, PermutationIndex> QrType;
  typedef Inverse<QrType> SrcXprType;
  static void run(DstXprType& dst, const SrcXprType& src,
                  const internal::assign_op<typename DstXprType::Scalar, typename QrType::Scalar>&) {
    dst = src.nestedExpression().solve(MatrixType::Identity(src.rows(), src.cols()));
  }
};

}  // end namespace internal

template <typename MatrixType, typename PermutationIndex>
typename RandColPivHouseholderQR<MatrixType, PermutationIndex>::HouseholderSequenceType
RandColPivHouseholderQR<MatrixType, PermutationIndex>::householderQ() const {
  eigen_assert(m_isInitialized && "RandColPivHouseholderQR is not initialized.");
  return HouseholderSequenceType(m_qr, m_hCoeffs.conjugate());
}

/** \return the randomized column-pivoted Householder QR decomposition of \c *this.
 *
 * \sa class RandColPivHouseholderQR
 */
template <typename Derived>
template <typename PermutationIndexType>
RandColPivHouseholderQR<typename MatrixBase<Derived>::PlainObject, PermutationIndexType>
MatrixBase<Derived>::randColPivHouseholderQr() const {
  return RandColPivHouseholderQR<PlainObject, PermutationIndexType>(eval());
}

}  // end namespace Eigen

#endif  // EIGEN_RANDCOLPIVOTINGHOUSEHOLDERQR_H
