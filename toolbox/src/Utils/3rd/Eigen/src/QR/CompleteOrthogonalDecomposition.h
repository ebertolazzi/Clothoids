// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2016 Rasmus Munk Larsen <rmlarsen@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

#ifndef EIGEN_COMPLETEORTHOGONALDECOMPOSITION_H
#define EIGEN_COMPLETEORTHOGONALDECOMPOSITION_H

// IWYU pragma: private
#include "./InternalHeaderCheck.h"

namespace Eigen {

namespace internal {

template <typename MatrixType_, typename PermutationIndex_, template <typename, typename> class RankRevealingQR_>
class CompleteOrthogonalDecompositionImpl;

template <typename MatrixType_, typename PermutationIndex_, template <typename, typename> class RankRevealingQR_>
struct traits<CompleteOrthogonalDecompositionImpl<MatrixType_, PermutationIndex_, RankRevealingQR_>>
    : traits<MatrixType_> {
  typedef MatrixXpr XprKind;
  typedef SolverStorage StorageKind;
  typedef PermutationIndex_ PermutationIndex;
  enum { Flags = 0 };
};

template <typename MatrixType_, typename PermutationIndex_>
struct traits<CompleteOrthogonalDecomposition<MatrixType_, PermutationIndex_>> : traits<MatrixType_> {
  typedef MatrixXpr XprKind;
  typedef SolverStorage StorageKind;
  typedef PermutationIndex_ PermutationIndex;
  enum { Flags = 0 };
};

template <typename MatrixType_, typename PermutationIndex_>
struct traits<RandCompleteOrthogonalDecomposition<MatrixType_, PermutationIndex_>> : traits<MatrixType_> {
  typedef MatrixXpr XprKind;
  typedef SolverStorage StorageKind;
  typedef PermutationIndex_ PermutationIndex;
  enum { Flags = 0 };
};

}  // end namespace internal

namespace internal {

/** \internal
 *
 * \class CompleteOrthogonalDecompositionImpl
 *
 * \brief Implementation backbone for CompleteOrthogonalDecomposition and
 * RandCompleteOrthogonalDecomposition.
 *
 * Parameterized over the rank-revealing QR engine \c RankRevealingQR_. The
 * public, source-compatible class \c CompleteOrthogonalDecomposition is a
 * thin shim around \c CompleteOrthogonalDecompositionImpl with
 * \c ColPivHouseholderQR; \c RandCompleteOrthogonalDecomposition is the
 * parallel shim with \c RandColPivHouseholderQR.
 */
template <typename MatrixType_, typename PermutationIndex_, template <typename, typename> class RankRevealingQR_>
class CompleteOrthogonalDecompositionImpl
    : public SolverBase<CompleteOrthogonalDecompositionImpl<MatrixType_, PermutationIndex_, RankRevealingQR_>> {
 public:
  typedef MatrixType_ MatrixType;
  typedef SolverBase<CompleteOrthogonalDecompositionImpl> Base;

  template <typename Derived>
  friend struct internal::solve_assertion;
  typedef PermutationIndex_ PermutationIndex;
  EIGEN_GENERIC_PUBLIC_INTERFACE(CompleteOrthogonalDecompositionImpl)
  enum {
    MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
  };
  typedef typename internal::plain_diag_type<MatrixType>::type HCoeffsType;
  typedef PermutationMatrix<ColsAtCompileTime, MaxColsAtCompileTime, PermutationIndex> PermutationType;
  typedef typename internal::plain_row_type<MatrixType, Index>::type IntRowVectorType;
  typedef typename internal::plain_row_type<MatrixType>::type RowVectorType;
  typedef typename internal::plain_row_type<MatrixType, RealScalar>::type RealRowVectorType;
  typedef HouseholderSequence<MatrixType, internal::remove_all_t<typename HCoeffsType::ConjugateReturnType>>
      HouseholderSequenceType;
  typedef typename MatrixType::PlainObject PlainObject;
  typedef RankRevealingQR_<MatrixType, PermutationIndex> RankRevealingQRType;

 public:
  CompleteOrthogonalDecompositionImpl() : m_cpqr(), m_zCoeffs(), m_temp() {}

  CompleteOrthogonalDecompositionImpl(Index rows, Index cols)
      : m_cpqr(rows, cols), m_zCoeffs((std::min)(rows, cols)), m_temp(cols) {}

  template <typename InputType>
  explicit CompleteOrthogonalDecompositionImpl(const EigenBase<InputType>& matrix)
      : m_cpqr(matrix.rows(), matrix.cols()),
        m_zCoeffs((std::min)(matrix.rows(), matrix.cols())),
        m_temp(matrix.cols()) {
    compute(matrix.derived());
  }

  template <typename InputType>
  explicit CompleteOrthogonalDecompositionImpl(EigenBase<InputType>& matrix)
      : m_cpqr(matrix.derived()), m_zCoeffs((std::min)(matrix.rows(), matrix.cols())), m_temp(matrix.cols()) {
    computeInPlace();
  }

  HouseholderSequenceType householderQ() const;
  HouseholderSequenceType matrixQ() const { return m_cpqr.householderQ(); }

  MatrixType matrixZ() const {
    MatrixType Z = MatrixType::Identity(m_cpqr.cols(), m_cpqr.cols());
    applyZOnTheLeftInPlace<false>(Z);
    return Z;
  }

  const MatrixType& matrixQTZ() const { return m_cpqr.matrixQR(); }
  const MatrixType& matrixT() const { return m_cpqr.matrixQR(); }

  template <typename InputType>
  CompleteOrthogonalDecompositionImpl& compute(const EigenBase<InputType>& matrix) {
    m_cpqr.compute(matrix);
    computeInPlace();
    return *this;
  }

  const PermutationType& colsPermutation() const { return m_cpqr.colsPermutation(); }

  typename MatrixType::Scalar determinant() const;
  typename MatrixType::RealScalar absDeterminant() const;
  typename MatrixType::RealScalar logAbsDeterminant() const;
  typename MatrixType::Scalar signDeterminant() const;

  inline Index rank() const { return m_cpqr.rank(); }
  inline Index dimensionOfKernel() const { return m_cpqr.dimensionOfKernel(); }
  inline bool isInjective() const { return m_cpqr.isInjective(); }
  inline bool isSurjective() const { return m_cpqr.isSurjective(); }
  inline bool isInvertible() const { return m_cpqr.isInvertible(); }

  inline Index rows() const { return m_cpqr.rows(); }
  inline Index cols() const { return m_cpqr.cols(); }

  inline const HCoeffsType& hCoeffs() const { return m_cpqr.hCoeffs(); }
  const HCoeffsType& zCoeffs() const { return m_zCoeffs; }

  CompleteOrthogonalDecompositionImpl& setThreshold(const RealScalar& threshold) {
    m_cpqr.setThreshold(threshold);
    return *this;
  }

  CompleteOrthogonalDecompositionImpl& setThreshold(Default_t) {
    m_cpqr.setThreshold(Default);
    return *this;
  }

  RealScalar threshold() const { return m_cpqr.threshold(); }

  inline Index nonzeroPivots() const { return m_cpqr.nonzeroPivots(); }
  inline RealScalar maxPivot() const { return m_cpqr.maxPivot(); }

  ComputationInfo info() const {
    eigen_assert(m_cpqr.m_isInitialized && "Decomposition is not initialized.");
    return Success;
  }

  /** \internal Asserts that compute() has been called. */
  void check_initialized() const {
    eigen_assert(m_cpqr.m_isInitialized && "CompleteOrthogonalDecomposition is not initialized.");
  }

#ifndef EIGEN_PARSED_BY_DOXYGEN
  template <typename RhsType, typename DstType>
  void _solve_impl(const RhsType& rhs, DstType& dst) const;

  template <bool Conjugate, typename RhsType, typename DstType>
  void _solve_impl_transposed(const RhsType& rhs, DstType& dst) const;
#endif

 protected:
  EIGEN_STATIC_ASSERT_NON_INTEGER(Scalar)

  template <bool Transpose_, typename Rhs>
  void _check_solve_assertion(const Rhs& b) const {
    EIGEN_ONLY_USED_FOR_DEBUG(b);
    eigen_assert(m_cpqr.m_isInitialized && "CompleteOrthogonalDecomposition is not initialized.");
    eigen_assert((Transpose_ ? this->cols() : this->rows()) == b.rows() &&
                 "CompleteOrthogonalDecomposition::solve(): invalid number of rows of the right hand side matrix b");
  }

  void computeInPlace();

  template <bool Conjugate, typename Rhs>
  void applyZOnTheLeftInPlace(Rhs& rhs) const;

  template <typename Rhs>
  void applyZAdjointOnTheLeftInPlace(Rhs& rhs) const;

  RankRevealingQRType m_cpqr;
  HCoeffsType m_zCoeffs;
  RowVectorType m_temp;
};

template <typename MatrixType, typename PermutationIndex, template <typename, typename> class RankRevealingQR_>
typename MatrixType::Scalar
CompleteOrthogonalDecompositionImpl<MatrixType, PermutationIndex, RankRevealingQR_>::determinant() const {
  return m_cpqr.determinant();
}

template <typename MatrixType, typename PermutationIndex, template <typename, typename> class RankRevealingQR_>
typename MatrixType::RealScalar
CompleteOrthogonalDecompositionImpl<MatrixType, PermutationIndex, RankRevealingQR_>::absDeterminant() const {
  return m_cpqr.absDeterminant();
}

template <typename MatrixType, typename PermutationIndex, template <typename, typename> class RankRevealingQR_>
typename MatrixType::RealScalar
CompleteOrthogonalDecompositionImpl<MatrixType, PermutationIndex, RankRevealingQR_>::logAbsDeterminant() const {
  return m_cpqr.logAbsDeterminant();
}

template <typename MatrixType, typename PermutationIndex, template <typename, typename> class RankRevealingQR_>
typename MatrixType::Scalar
CompleteOrthogonalDecompositionImpl<MatrixType, PermutationIndex, RankRevealingQR_>::signDeterminant() const {
  return m_cpqr.signDeterminant();
}

template <typename MatrixType, typename PermutationIndex, template <typename, typename> class RankRevealingQR_>
void CompleteOrthogonalDecompositionImpl<MatrixType, PermutationIndex, RankRevealingQR_>::computeInPlace() {
  eigen_assert(m_cpqr.cols() <= NumTraits<PermutationIndex>::highest());

  const Index rank = m_cpqr.rank();
  const Index cols = m_cpqr.cols();
  const Index rows = m_cpqr.rows();
  m_zCoeffs.resize((std::min)(rows, cols));
  m_temp.resize(cols);

  if (rank < cols) {
    // We have reduced the (permuted) matrix to the form
    //   [R11 R12]
    //   [ 0  R22]
    // where R11 is r-by-r (r = rank) upper triangular, R12 is
    // r-by-(n-r), and R22 is empty or the norm of R22 is negligible.
    // We now compute the complete orthogonal decomposition by applying
    // Householder transformations from the right to the upper trapezoidal
    // matrix X = [R11 R12] to zero out R12 and obtain the factorization
    // [R11 R12] = [T11 0] * Z, where T11 is r-by-r upper triangular and
    // Z = Z(0) * Z(1) ... Z(r-1) is an n-by-n orthogonal matrix.
    // We store the data representing Z in R12 and m_zCoeffs.
    for (Index k = rank - 1; k >= 0; --k) {
      if (k != rank - 1) {
        // Given the API for Householder reflectors, it is more convenient if
        // we swap the leading parts of columns k and r-1 (zero-based) to form
        // the matrix X_k = [X(0:k, k), X(0:k, r:n)]
        m_cpqr.m_qr.col(k).head(k + 1).swap(m_cpqr.m_qr.col(rank - 1).head(k + 1));
      }
      // Construct Householder reflector Z(k) to zero out the last row of X_k,
      // i.e. choose Z(k) such that
      // [X(k, k), X(k, r:n)] * Z(k) = [beta, 0, .., 0].
      RealScalar beta;
      m_cpqr.m_qr.row(k).tail(cols - rank + 1).makeHouseholderInPlace(m_zCoeffs(k), beta);
      m_cpqr.m_qr(k, rank - 1) = beta;
      if (k > 0) {
        // Apply Z(k) to the first k rows of X_k
        m_cpqr.m_qr.topRightCorner(k, cols - rank + 1)
            .applyHouseholderOnTheRight(m_cpqr.m_qr.row(k).tail(cols - rank).adjoint(), m_zCoeffs(k), &m_temp(0));
      }
      if (k != rank - 1) {
        // Swap X(0:k,k) back to its proper location.
        m_cpqr.m_qr.col(k).head(k + 1).swap(m_cpqr.m_qr.col(rank - 1).head(k + 1));
      }
    }
  }
}

template <typename MatrixType, typename PermutationIndex, template <typename, typename> class RankRevealingQR_>
template <bool Conjugate, typename Rhs>
void CompleteOrthogonalDecompositionImpl<MatrixType, PermutationIndex, RankRevealingQR_>::applyZOnTheLeftInPlace(
    Rhs& rhs) const {
  const Index cols = this->cols();
  const Index nrhs = rhs.cols();
  const Index rank = this->rank();
  Matrix<typename Rhs::Scalar, Dynamic, 1> temp((std::max)(cols, nrhs));
  for (Index k = rank - 1; k >= 0; --k) {
    if (k != rank - 1) {
      rhs.row(k).swap(rhs.row(rank - 1));
    }
    rhs.middleRows(rank - 1, cols - rank + 1)
        .applyHouseholderOnTheLeft(matrixQTZ().row(k).tail(cols - rank).transpose().template conjugateIf<!Conjugate>(),
                                   zCoeffs().template conjugateIf<Conjugate>()(k), &temp(0));
    if (k != rank - 1) {
      rhs.row(k).swap(rhs.row(rank - 1));
    }
  }
}

template <typename MatrixType, typename PermutationIndex, template <typename, typename> class RankRevealingQR_>
template <typename Rhs>
void CompleteOrthogonalDecompositionImpl<MatrixType, PermutationIndex, RankRevealingQR_>::applyZAdjointOnTheLeftInPlace(
    Rhs& rhs) const {
  const Index cols = this->cols();
  const Index nrhs = rhs.cols();
  const Index rank = this->rank();
  Matrix<typename Rhs::Scalar, Dynamic, 1> temp((std::max)(cols, nrhs));
  for (Index k = 0; k < rank; ++k) {
    if (k != rank - 1) {
      rhs.row(k).swap(rhs.row(rank - 1));
    }
    rhs.middleRows(rank - 1, cols - rank + 1)
        .applyHouseholderOnTheLeft(matrixQTZ().row(k).tail(cols - rank).adjoint(), zCoeffs()(k), &temp(0));
    if (k != rank - 1) {
      rhs.row(k).swap(rhs.row(rank - 1));
    }
  }
}

#ifndef EIGEN_PARSED_BY_DOXYGEN
template <typename MatrixType, typename PermutationIndex, template <typename, typename> class RankRevealingQR_>
template <typename RhsType, typename DstType>
void CompleteOrthogonalDecompositionImpl<MatrixType, PermutationIndex, RankRevealingQR_>::_solve_impl(
    const RhsType& rhs, DstType& dst) const {
  const Index rank = this->rank();
  if (rank == 0) {
    dst.setZero();
    return;
  }

  // Compute c = Q^* * rhs
  typename RhsType::PlainObject c(rhs);
  c.applyOnTheLeft(matrixQ().setLength(rank).adjoint());

  // Solve T z = c(1:rank, :)
  dst.topRows(rank) = matrixT().topLeftCorner(rank, rank).template triangularView<Upper>().solve(c.topRows(rank));

  const Index cols = this->cols();
  if (rank < cols) {
    // Compute y = Z^* * [ z ]
    //                   [ 0 ]
    dst.bottomRows(cols - rank).setZero();
    applyZAdjointOnTheLeftInPlace(dst);
  }

  // Undo permutation to get x = P^{-1} * y.
  dst = colsPermutation() * dst;
}

template <typename MatrixType, typename PermutationIndex, template <typename, typename> class RankRevealingQR_>
template <bool Conjugate, typename RhsType, typename DstType>
void CompleteOrthogonalDecompositionImpl<MatrixType, PermutationIndex, RankRevealingQR_>::_solve_impl_transposed(
    const RhsType& rhs, DstType& dst) const {
  const Index rank = this->rank();

  if (rank == 0) {
    dst.setZero();
    return;
  }

  typename RhsType::PlainObject c(colsPermutation().transpose() * rhs);

  if (rank < cols()) {
    applyZOnTheLeftInPlace<!Conjugate>(c);
  }

  matrixT()
      .topLeftCorner(rank, rank)
      .template triangularView<Upper>()
      .transpose()
      .template conjugateIf<Conjugate>()
      .solveInPlace(c.topRows(rank));

  dst.topRows(rank) = c.topRows(rank);
  dst.bottomRows(rows() - rank).setZero();

  dst.applyOnTheLeft(householderQ().setLength(rank).template conjugateIf<!Conjugate>());
}
#endif

template <typename MatrixType, typename PermutationIndex, template <typename, typename> class RankRevealingQR_>
typename CompleteOrthogonalDecompositionImpl<MatrixType, PermutationIndex, RankRevealingQR_>::HouseholderSequenceType
CompleteOrthogonalDecompositionImpl<MatrixType, PermutationIndex, RankRevealingQR_>::householderQ() const {
  return m_cpqr.householderQ();
}

}  // end namespace internal

/** \ingroup QR_Module
 *
 * \class CompleteOrthogonalDecomposition
 *
 * \brief Complete orthogonal decomposition (COD) of a matrix.
 *
 * \tparam MatrixType_ the type of the matrix of which we are computing the COD.
 *
 * This class performs a rank-revealing complete orthogonal decomposition of a
 * matrix  \b A into matrices \b P, \b Q, \b T, and \b Z such that
 * \f[
 *  \mathbf{A} \, \mathbf{P} = \mathbf{Q} \,
 *                     \begin{bmatrix} \mathbf{T} &  \mathbf{0} \\
 *                                     \mathbf{0} & \mathbf{0} \end{bmatrix} \, \mathbf{Z}
 * \f]
 * by using Householder transformations. Here, \b P is a permutation matrix,
 * \b Q and \b Z are unitary matrices and \b T an upper triangular matrix of
 * size rank-by-rank. \b A may be rank deficient.
 *
 * Internally backed by ColPivHouseholderQR. For large matrices, see
 * RandCompleteOrthogonalDecomposition, which uses a randomized blocked
 * column-pivoted QR and is faster on level-3-BLAS-friendly inputs.
 *
 * This class supports the \link InplaceDecomposition inplace decomposition \endlink mechanism.
 *
 * \sa MatrixBase::completeOrthogonalDecomposition(), RandCompleteOrthogonalDecomposition
 */
template <typename MatrixType_, typename PermutationIndex_>
class CompleteOrthogonalDecomposition
    : public internal::CompleteOrthogonalDecompositionImpl<MatrixType_, PermutationIndex_, ColPivHouseholderQR> {
 public:
  typedef internal::CompleteOrthogonalDecompositionImpl<MatrixType_, PermutationIndex_, ColPivHouseholderQR> Base;
  using typename Base::RealScalar;

  CompleteOrthogonalDecomposition() : Base() {}
  CompleteOrthogonalDecomposition(Index rows, Index cols) : Base(rows, cols) {}

  template <typename InputType>
  explicit CompleteOrthogonalDecomposition(const EigenBase<InputType>& matrix) : Base(matrix.derived()) {}

  template <typename InputType>
  explicit CompleteOrthogonalDecomposition(EigenBase<InputType>& matrix) : Base(matrix.derived()) {}

  /** \brief Computes the COD of \a matrix. \sa class CompleteOrthogonalDecomposition */
  template <typename InputType>
  CompleteOrthogonalDecomposition& compute(const EigenBase<InputType>& matrix) {
    Base::compute(matrix);
    return *this;
  }

  CompleteOrthogonalDecomposition& setThreshold(const RealScalar& threshold) {
    Base::setThreshold(threshold);
    return *this;
  }

  CompleteOrthogonalDecomposition& setThreshold(Default_t) {
    Base::setThreshold(Default);
    return *this;
  }

  /** \returns the pseudo-inverse of the matrix of which *this is the complete
   * orthogonal decomposition.
   * \warning Do not compute \c this->pseudoInverse()*rhs to solve a linear system.
   * It is more efficient and numerically stable to call \c this->solve(rhs).
   */
  inline Inverse<CompleteOrthogonalDecomposition> pseudoInverse() const {
    this->check_initialized();
    return Inverse<CompleteOrthogonalDecomposition>(*this);
  }
};

/** \ingroup QR_Module
 *
 * \class RandCompleteOrthogonalDecomposition
 *
 * \brief Complete orthogonal decomposition (COD) of a matrix, backed by
 * RandColPivHouseholderQR.
 *
 * Same factorization as CompleteOrthogonalDecomposition, but the rank-revealing
 * QR step uses the randomized blocked column-pivoted Householder QR
 * (BQRRP framework). On large matrices this is faster than the classical
 * Businger-Golub pivoting in CompleteOrthogonalDecomposition because almost
 * all work is performed in level-3 BLAS. The pivot-quality difference is
 * empirically minor; see RandColPivHouseholderQR for references and details.
 *
 * The block size and RNG seed of the underlying QR can be configured via
 * setBlockSize() and setSeed() before calling compute().
 *
 * For small or fixed-size matrices prefer CompleteOrthogonalDecomposition:
 * the sketch overhead dominates below roughly a few hundred rows/columns.
 *
 * \sa MatrixBase::randCompleteOrthogonalDecomposition(),
 *     CompleteOrthogonalDecomposition, RandColPivHouseholderQR
 */
template <typename MatrixType_, typename PermutationIndex_>
class RandCompleteOrthogonalDecomposition
    : public internal::CompleteOrthogonalDecompositionImpl<MatrixType_, PermutationIndex_, RandColPivHouseholderQR> {
 public:
  typedef internal::CompleteOrthogonalDecompositionImpl<MatrixType_, PermutationIndex_, RandColPivHouseholderQR> Base;
  using typename Base::RealScalar;

  RandCompleteOrthogonalDecomposition() : Base() {}
  RandCompleteOrthogonalDecomposition(Index rows, Index cols) : Base(rows, cols) {}

  template <typename InputType>
  explicit RandCompleteOrthogonalDecomposition(const EigenBase<InputType>& matrix) : Base(matrix.derived()) {}

  template <typename InputType>
  explicit RandCompleteOrthogonalDecomposition(EigenBase<InputType>& matrix) : Base(matrix.derived()) {}

  template <typename InputType>
  RandCompleteOrthogonalDecomposition& compute(const EigenBase<InputType>& matrix) {
    Base::compute(matrix);
    return *this;
  }

  RandCompleteOrthogonalDecomposition& setThreshold(const RealScalar& threshold) {
    Base::setThreshold(threshold);
    return *this;
  }

  RandCompleteOrthogonalDecomposition& setThreshold(Default_t) {
    Base::setThreshold(Default);
    return *this;
  }

  /** \brief Sets the panel block size of the underlying randomized QR.
   * \sa RandColPivHouseholderQR::setBlockSize
   */
  RandCompleteOrthogonalDecomposition& setBlockSize(Index b) {
    this->m_cpqr.setBlockSize(b);
    return *this;
  }

  /** \brief Fixes the RNG seed of the underlying randomized QR.
   * \sa RandColPivHouseholderQR::setSeed
   */
  RandCompleteOrthogonalDecomposition& setSeed(uint64_t seed) {
    this->m_cpqr.setSeed(seed);
    return *this;
  }

  inline Inverse<RandCompleteOrthogonalDecomposition> pseudoInverse() const {
    this->check_initialized();
    return Inverse<RandCompleteOrthogonalDecomposition>(*this);
  }
};

namespace internal {

template <typename MatrixType, typename PermutationIndex>
struct traits<Inverse<CompleteOrthogonalDecomposition<MatrixType, PermutationIndex>>>
    : traits<typename Transpose<typename MatrixType::PlainObject>::PlainObject> {
  enum { Flags = 0 };
};

template <typename DstXprType, typename MatrixType, typename PermutationIndex>
struct Assignment<DstXprType, Inverse<CompleteOrthogonalDecomposition<MatrixType, PermutationIndex>>,
                  internal::assign_op<typename DstXprType::Scalar,
                                      typename CompleteOrthogonalDecomposition<MatrixType, PermutationIndex>::Scalar>,
                  Dense2Dense> {
  typedef CompleteOrthogonalDecomposition<MatrixType, PermutationIndex> CodType;
  typedef Inverse<CodType> SrcXprType;
  static void run(DstXprType& dst, const SrcXprType& src,
                  const internal::assign_op<typename DstXprType::Scalar, typename CodType::Scalar>&) {
    typedef Matrix<typename CodType::Scalar, CodType::RowsAtCompileTime, CodType::RowsAtCompileTime, 0,
                   CodType::MaxRowsAtCompileTime, CodType::MaxRowsAtCompileTime>
        IdentityMatrixType;
    dst = src.nestedExpression().solve(IdentityMatrixType::Identity(src.cols(), src.cols()));
  }
};

template <typename MatrixType, typename PermutationIndex>
struct traits<Inverse<RandCompleteOrthogonalDecomposition<MatrixType, PermutationIndex>>>
    : traits<typename Transpose<typename MatrixType::PlainObject>::PlainObject> {
  enum { Flags = 0 };
};

template <typename DstXprType, typename MatrixType, typename PermutationIndex>
struct Assignment<DstXprType, Inverse<RandCompleteOrthogonalDecomposition<MatrixType, PermutationIndex>>,
                  internal::assign_op<typename DstXprType::Scalar, typename RandCompleteOrthogonalDecomposition<
                                                                       MatrixType, PermutationIndex>::Scalar>,
                  Dense2Dense> {
  typedef RandCompleteOrthogonalDecomposition<MatrixType, PermutationIndex> CodType;
  typedef Inverse<CodType> SrcXprType;
  static void run(DstXprType& dst, const SrcXprType& src,
                  const internal::assign_op<typename DstXprType::Scalar, typename CodType::Scalar>&) {
    typedef Matrix<typename CodType::Scalar, CodType::RowsAtCompileTime, CodType::RowsAtCompileTime, 0,
                   CodType::MaxRowsAtCompileTime, CodType::MaxRowsAtCompileTime>
        IdentityMatrixType;
    dst = src.nestedExpression().solve(IdentityMatrixType::Identity(src.cols(), src.cols()));
  }
};

}  // end namespace internal

/** \return the complete orthogonal decomposition of \c *this.
 *
 * \sa class CompleteOrthogonalDecomposition
 */
template <typename Derived>
template <typename PermutationIndex>
CompleteOrthogonalDecomposition<typename MatrixBase<Derived>::PlainObject, PermutationIndex>
MatrixBase<Derived>::completeOrthogonalDecomposition() const {
  return CompleteOrthogonalDecomposition<PlainObject, PermutationIndex>(eval());
}

/** \return the randomized complete orthogonal decomposition of \c *this.
 *
 * \sa class RandCompleteOrthogonalDecomposition
 */
template <typename Derived>
template <typename PermutationIndex>
RandCompleteOrthogonalDecomposition<typename MatrixBase<Derived>::PlainObject, PermutationIndex>
MatrixBase<Derived>::randCompleteOrthogonalDecomposition() const {
  return RandCompleteOrthogonalDecomposition<PlainObject, PermutationIndex>(eval());
}

}  // end namespace Eigen

#endif  // EIGEN_COMPLETEORTHOGONALDECOMPOSITION_H
