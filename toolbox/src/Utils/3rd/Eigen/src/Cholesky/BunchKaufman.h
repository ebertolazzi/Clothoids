// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Rasmus Munk Larsen <rmlarsen@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

#ifndef EIGEN_BUNCHKAUFMAN_H
#define EIGEN_BUNCHKAUFMAN_H

// IWYU pragma: private
#include "./InternalHeaderCheck.h"

namespace Eigen {

namespace internal {
template <typename MatrixType_, int UpLo_>
struct traits<BunchKaufman<MatrixType_, UpLo_> > : traits<MatrixType_> {
  typedef MatrixXpr XprKind;
  typedef SolverStorage StorageKind;
  typedef int StorageIndex;
  enum { Flags = 0 };
};

template <typename MatrixType, int UpLo>
struct BunchKaufman_Traits;

// Panel width for the blocked factorization (defined below); forward-declared so the size
// constructor can pre-allocate the panel workspace.
template <typename Scalar>
inline Index bunch_kaufman_blocksize();
}  // namespace internal

/** \ingroup Cholesky_Module
 *
 * \class BunchKaufman
 *
 * \brief Bunch-Kaufman factorization of a symmetric / Hermitian indefinite matrix
 *
 * \tparam MatrixType_ the type of the matrix of which to compute the Bunch-Kaufman factorization
 * \tparam UpLo_ the triangular part that will be used for the decomposition: Lower (default) or Upper.
 *             The other triangular part won't be read.
 *
 * Perform a symmetric-indefinite factorization of a real symmetric or complex Hermitian matrix
 * \f$ A \f$ using the Bunch-Kaufman diagonal pivoting method, such that
 * \f$ P A P^T = L D L^* \f$ (or \f$ P A P^T = U^* D U \f$ for the upper variant), where \f$ P \f$ is a
 * permutation matrix, \f$ L \f$ is unit lower triangular and \f$ D \f$ is block diagonal with
 * 1x1 and 2x2 diagonal blocks.
 *
 * Unlike LLT and LDLT, this factorization is numerically stable (backward stable, with a bounded
 * element-growth factor) for \em indefinite symmetric / Hermitian matrices. It is the dense analog of
 * LAPACK's xSYTRF / xHETRF routines. For positive- or negative-(semi)definite matrices LLT / LDLT are
 * cheaper and should be preferred.
 *
 * The diagonal-pivoting threshold \f$ \alpha = (1+\sqrt{17})/8 \f$ and the pivot-selection logic follow:
 *  - J. R. Bunch and L. Kaufman, "Some stable methods for calculating inertia and solving symmetric
 *    linear systems", Math. Comp. 31 (1977), pp. 163-179.
 *  - G. H. Golub and C. F. Van Loan, "Matrix Computations", 4th ed., section 4.4.
 *  - N. J. Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed., chapter 11.
 *
 * The blocked (level-3 BLAS) variant follows the panel/trailing-update structure of LAPACK's
 * xSYTRF / xLASYF.
 *
 * This factorization also reveals the \em inertia of \f$ A \f$ (the number of positive, negative and
 * zero eigenvalues) by Sylvester's law of inertia, accessible through isPositive() / isNegative().
 *
 * This class supports the \link InplaceDecomposition inplace decomposition \endlink mechanism.
 *
 * \sa MatrixBase::bunchKaufman(), SelfAdjointView::bunchKaufman(), class LDLT, class LLT
 */
template <typename MatrixType_, int UpLo_>
class BunchKaufman : public SolverBase<BunchKaufman<MatrixType_, UpLo_> > {
 public:
  typedef MatrixType_ MatrixType;
  typedef SolverBase<BunchKaufman> Base;
  friend class SolverBase<BunchKaufman>;

  EIGEN_GENERIC_PUBLIC_INTERFACE(BunchKaufman)
  enum {
    MaxRowsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime,
    UpLo = UpLo_
  };

  typedef Matrix<Scalar, RowsAtCompileTime, 1, 0, MaxRowsAtCompileTime, 1> TmpVectorType;
  // Panel workspace for the blocked algorithm (only allocated for large dynamic-sized problems).
  typedef Matrix<Scalar, Dynamic, Dynamic> WorkspaceType;
  typedef Transpositions<RowsAtCompileTime, MaxRowsAtCompileTime> TranspositionType;
  typedef PermutationMatrix<RowsAtCompileTime, MaxRowsAtCompileTime> PermutationType;

  typedef internal::BunchKaufman_Traits<MatrixType, UpLo> Traits;

  /** \brief Default Constructor.
   *
   * The default constructor is useful in cases in which the user intends to
   * perform decompositions via BunchKaufman::compute(const MatrixType&).
   */
  BunchKaufman()
      : m_matrix(),
        m_l1_norm(0),
        m_transpositions(),
        m_subdiag(),
        m_n_pos(0),
        m_n_neg(0),
        m_n_zero(0),
        m_isInitialized(false),
        m_info(InvalidInput) {}

  /** \brief Default Constructor with memory preallocation
   *
   * Like the default constructor but with preallocation of the internal data
   * according to the specified problem \a size.
   * \sa BunchKaufman()
   */
  explicit BunchKaufman(Index size)
      : m_matrix(size, size),
        m_l1_norm(0),
        m_transpositions(size),
        m_subdiag(size),
        // Pre-allocate the panel workspace to the exact shape compute() needs, so that a subsequent
        // compute() on a problem of this size performs no heap allocation (the blocked factorization
        // resizes it to n x (nb+1); resizing to the same shape is a no-op).
        m_workspace(size, internal::bunch_kaufman_blocksize<Scalar>() + 1),
        m_n_pos(0),
        m_n_neg(0),
        m_n_zero(0),
        m_isInitialized(false),
        m_info(InvalidInput) {}

  /** \brief Constructor with decomposition
   *
   * This calculates the decomposition for the input \a matrix.
   *
   * \sa BunchKaufman(Index size)
   */
  template <typename InputType>
  explicit BunchKaufman(const EigenBase<InputType>& matrix)
      : m_matrix(matrix.rows(), matrix.cols()),
        m_l1_norm(0),
        m_transpositions(matrix.rows()),
        m_subdiag(matrix.rows()),
        m_n_pos(0),
        m_n_neg(0),
        m_n_zero(0),
        m_isInitialized(false),
        m_info(InvalidInput) {
    compute(matrix.derived());
  }

  /** \brief Constructs a Bunch-Kaufman factorization from a given matrix
   *
   * This overloaded constructor is provided for \link InplaceDecomposition inplace decomposition \endlink when \c
   * MatrixType is a Eigen::Ref.
   *
   * \sa BunchKaufman(const EigenBase&)
   */
  template <typename InputType>
  explicit BunchKaufman(EigenBase<InputType>& matrix)
      : m_matrix(matrix.derived()),
        m_l1_norm(0),
        m_transpositions(matrix.rows()),
        m_subdiag(matrix.rows()),
        m_n_pos(0),
        m_n_neg(0),
        m_n_zero(0),
        m_isInitialized(false),
        m_info(InvalidInput) {
    compute(matrix.derived());
  }

  /** \returns a view of the unit upper triangular matrix U */
  inline typename Traits::MatrixU matrixU() const {
    eigen_assert(m_isInitialized && "BunchKaufman is not initialized.");
    return Traits::getU(m_matrix);
  }

  /** \returns a view of the unit lower triangular matrix L */
  inline typename Traits::MatrixL matrixL() const {
    eigen_assert(m_isInitialized && "BunchKaufman is not initialized.");
    return Traits::getL(m_matrix);
  }

  /** \returns the permutation matrix P as a transposition sequence.
   */
  inline const TranspositionType& transpositionsP() const {
    eigen_assert(m_isInitialized && "BunchKaufman is not initialized.");
    return m_transpositions;
  }

  /** \returns the main diagonal of the block-diagonal matrix D.
   *
   * The block-diagonal matrix D is fully described by vectorD() (its main diagonal) together with
   * subDiagonal() (the sub-diagonal entries of its 2x2 blocks).
   */
  inline Diagonal<const MatrixType> vectorD() const {
    eigen_assert(m_isInitialized && "BunchKaufman is not initialized.");
    return m_matrix.diagonal();
  }

  /** \returns the sub-diagonal of the block-diagonal matrix D.
   *
   * Entry \c k is non-zero if and only if a 2x2 diagonal block of D occupies rows/columns \c k and
   * \c k+1; it then holds \f$ D_{k+1,k} \f$. All other entries are zero. \sa vectorD()
   */
  inline const TmpVectorType& subDiagonal() const {
    eigen_assert(m_isInitialized && "BunchKaufman is not initialized.");
    return m_subdiag;
  }

  /** \returns true if the matrix is positive semidefinite (has no strictly negative eigenvalues). */
  inline bool isPositive() const {
    eigen_assert(m_isInitialized && "BunchKaufman is not initialized.");
    return m_n_neg == 0;
  }

  /** \returns true if the matrix is negative semidefinite (has no strictly positive eigenvalues). */
  inline bool isNegative() const {
    eigen_assert(m_isInitialized && "BunchKaufman is not initialized.");
    return m_n_pos == 0;
  }

#ifdef EIGEN_PARSED_BY_DOXYGEN
  /** \returns a solution x of \f$ A x = b \f$ using the current decomposition of A.
   *
   * This function also supports in-place solves using the syntax <tt>x = decompositionObject.solve(x)</tt> .
   *
   * \note_about_checking_solutions
   *
   * \sa MatrixBase::bunchKaufman(), SelfAdjointView::bunchKaufman()
   */
  template <typename Rhs>
  inline Solve<BunchKaufman, Rhs> solve(const MatrixBase<Rhs>& b) const;
#endif

  template <typename Derived>
  bool solveInPlace(MatrixBase<Derived>& bAndX) const;

  template <typename InputType>
  BunchKaufman& compute(const EigenBase<InputType>& matrix);

  /** \returns an estimate of the reciprocal condition number of the matrix of
   *  which \c *this is the Bunch-Kaufman decomposition.
   */
  RealScalar rcond() const {
    eigen_assert(m_isInitialized && "BunchKaufman is not initialized.");
    return internal::rcond_estimate_helper(m_l1_norm, *this);
  }

  /** \returns the internal Bunch-Kaufman decomposition matrix
   *
   * The strictly lower (resp. upper) triangular part holds the unit triangular factor L (resp. U); the
   * diagonal holds the main diagonal of D. The sub/super-diagonal entries of the 2x2 blocks of D are
   * stored separately, in subDiagonal().
   */
  inline const MatrixType& matrixLDLT() const {
    eigen_assert(m_isInitialized && "BunchKaufman is not initialized.");
    return m_matrix;
  }

  MatrixType reconstructedMatrix() const;

  /** \returns the adjoint of \c *this, that is, a const reference to the decomposition itself as the underlying matrix
   * is self-adjoint.
   *
   * This method is provided for compatibility with other matrix decompositions, thus enabling generic code such as:
   * \code x = decomposition.adjoint().solve(b) \endcode
   */
  const BunchKaufman& adjoint() const { return *this; }

  EIGEN_DEVICE_FUNC constexpr Index rows() const noexcept { return m_matrix.rows(); }
  EIGEN_DEVICE_FUNC constexpr Index cols() const noexcept { return m_matrix.cols(); }

  /** \brief Reports whether previous computation was successful.
   *
   * \returns \c Success if computation was successful,
   *          \c NumericalIssue if the factorization failed because of a zero pivot (the matrix is singular).
   */
  ComputationInfo info() const {
    eigen_assert(m_isInitialized && "BunchKaufman is not initialized.");
    return m_info;
  }

#ifndef EIGEN_PARSED_BY_DOXYGEN
  template <typename RhsType, typename DstType>
  void _solve_impl(const RhsType& rhs, DstType& dst) const;

  template <bool Conjugate, typename RhsType, typename DstType>
  void _solve_impl_transposed(const RhsType& rhs, DstType& dst) const;
#endif

 protected:
  EIGEN_STATIC_ASSERT_NON_INTEGER(Scalar)

  /** \internal Apply the block-diagonal matrix D in place: \c x <- D x. */
  template <typename Derived>
  void applyD(MatrixBase<Derived>& x) const;

  /** \internal Apply the inverse of the block-diagonal matrix D in place: \c x <- D^{-1} x.
   * For \c Conjugate == false the transpose of D is used instead (relevant for the complex case). */
  template <bool Conjugate, typename Derived>
  void solveInPlaceD(MatrixBase<Derived>& x) const;

  /** \internal Compute the inertia (counts of positive / negative / zero eigenvalues) from D. */
  void computeInertia();

  MatrixType m_matrix;
  RealScalar m_l1_norm;
  TranspositionType m_transpositions;
  TmpVectorType m_subdiag;
  WorkspaceType m_workspace;
  Index m_n_pos;
  Index m_n_neg;
  Index m_n_zero;
  bool m_isInitialized;
  ComputationInfo m_info;
};

namespace internal {

/** \internal The Bunch-Kaufman diagonal-pivoting threshold \f$ \alpha = (1+\sqrt{17})/8 \approx 0.6404 \f$.
 * It minimizes the bound on element growth across a 1x1 followed by a 1x1 step versus a single 2x2 step
 * (Bunch & Kaufman, 1977). */
template <typename RealScalar>
EIGEN_DEVICE_FUNC inline RealScalar bunch_kaufman_alpha() {
  using std::sqrt;
  return (RealScalar(1) + sqrt(RealScalar(17))) / RealScalar(8);
}

/** \internal Panel width for the blocked (level-3) Bunch-Kaufman factorization. Can be overridden (mainly
 * for testing the panel logic at small sizes) by defining EIGEN_BUNCHKAUFMAN_BLOCKSIZE. */
template <typename Scalar>
inline Index bunch_kaufman_blocksize() {
#ifdef EIGEN_BUNCHKAUFMAN_BLOCKSIZE
  return Index(EIGEN_BUNCHKAUFMAN_BLOCKSIZE);
#else
  return 64;
#endif
}

template <int UpLo>
struct bunch_kaufman;

template <>
struct bunch_kaufman<Lower> {
  // Interchange rows and columns kk and kp (kp >= kk) in the lower triangle of the Hermitian matrix
  // `mat`, including the already-computed factor columns to the left of column `kfirst` (so that the
  // stored unit triangular factor stays consistent with a single up-front permutation P). `kfirst` is
  // the first column of the pivot block (== kk for a 1x1 pivot, == kk-1 for a 2x2 pivot).
  template <typename MatrixType>
  static void apply_symmetric_pivot(MatrixType& mat, Index kfirst, Index kk, Index kp, Index kstep) {
    typedef typename MatrixType::Scalar Scalar;
    const Index n = mat.rows();
    const Index s = n - kp - 1;
    if (s > 0) mat.col(kk).tail(s).swap(mat.col(kp).tail(s));
    for (Index i = kk + 1; i < kp; ++i) {
      Scalar tmp = mat.coeff(i, kk);
      mat.coeffRef(i, kk) = numext::conj(mat.coeff(kp, i));
      mat.coeffRef(kp, i) = numext::conj(tmp);
    }
    numext::swap(mat.coeffRef(kk, kk), mat.coeffRef(kp, kp));
    EIGEN_IF_CONSTEXPR (NumTraits<Scalar>::IsComplex) {
      mat.coeffRef(kp, kk) = numext::conj(mat.coeff(kp, kk));
    }
    if (kfirst > 0) mat.row(kk).head(kfirst).swap(mat.row(kp).head(kfirst));
    if (kstep == 2) {
      numext::swap(mat.coeffRef(kfirst + 1, kfirst), mat.coeffRef(kp, kfirst));
    }
  }

  // Unblocked (level-2 BLAS) Bunch-Kaufman factorization of the lower triangle of `mat`, in place.
  // Columns [k0, n) are factorized. On output:
  //   - the strictly-lower triangle holds the unit lower factor L (with explicit zeros at the
  //     sub-diagonal positions of 2x2 blocks),
  //   - the diagonal holds the main diagonal of D,
  //   - subdiag(k) holds D(k+1,k) for each 2x2 block starting at column k (and 0 elsewhere),
  //   - transpositions encodes the symmetric permutation P (applied in increasing index order).
  // Returns 0 on success, or the (1-based) index of the first exactly-zero pivot encountered.
  template <typename MatrixType, typename TranspositionType, typename SubDiagType>
  static Index unblocked(MatrixType& mat, TranspositionType& transpositions, SubDiagType& subdiag, Index k0 = 0) {
    using numext::abs;
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::RealScalar RealScalar;
    typedef typename TranspositionType::StorageIndex StorageIndex;
    const Index n = mat.rows();
    const RealScalar alpha = bunch_kaufman_alpha<RealScalar>();
    Index info = 0;

    Index k = k0;
    while (k < n) {
      Index kstep = 1;
      Index kp = k;
      const RealScalar absakk = abs(numext::real(mat.coeff(k, k)));

      // colmax = max_{i>k} |mat(i,k)|, attained at row imax.
      Index imax = k;
      RealScalar colmax(0);
      if (k + 1 < n) {
        Index rel = 0;
        colmax = mat.col(k).tail(n - k - 1).cwiseAbs().maxCoeff(&rel);
        imax = k + 1 + rel;
      }

      if (numext::is_exactly_zero((numext::maxi)(absakk, colmax)) || (numext::isnan)(absakk)) {
        // The whole remaining column is zero: 1x1 zero pivot, matrix is singular.
        kp = k;
        if (info == 0) info = k + 1;
        mat.coeffRef(k, k) = Scalar(numext::real(mat.coeff(k, k)));
      } else if (absakk >= alpha * colmax) {
        // 1x1 pivot at k, no interchange.
        kp = k;
      } else {
        // rowmax = largest off-diagonal magnitude in row/column imax.
        RealScalar rowmax = mat.row(imax).segment(k, imax - k).cwiseAbs().maxCoeff();
        if (imax + 1 < n) {
          RealScalar rowmax2 = mat.col(imax).tail(n - imax - 1).cwiseAbs().maxCoeff();
          rowmax = (numext::maxi)(rowmax, rowmax2);
        }
        if (absakk >= alpha * colmax * (colmax / rowmax)) {
          // 1x1 pivot at k.
          kp = k;
        } else if (abs(numext::real(mat.coeff(imax, imax))) >= alpha * rowmax) {
          // 1x1 pivot at imax: interchange rows/columns k and imax.
          kp = imax;
        } else {
          // 2x2 pivot at (k, imax): interchange rows/columns k+1 and imax.
          kp = imax;
          kstep = 2;
        }
      }

      const Index kk = k + kstep - 1;  // column of the pivot block to interchange with kp
      if (kp != kk) apply_symmetric_pivot(mat, k, kk, kp, kstep);

      if (kstep == 1) {
        transpositions.coeffRef(k) = StorageIndex(kp);
        subdiag.coeffRef(k) = Scalar(0);

        const RealScalar dkk = numext::real(mat.coeff(k, k));
        mat.coeffRef(k, k) = Scalar(dkk);
        const Index rs = n - k - 1;
        if (rs > 0) {
          if (!numext::is_exactly_zero(dkk)) {
            // A22 <- A22 - (1/dkk) w w^*, with w = mat(k+1:n, k); then L column = w / dkk.
            auto w = mat.col(k).tail(rs);
            mat.block(k + 1, k + 1, rs, rs).template selfadjointView<Lower>().rankUpdate(w, RealScalar(-1) / dkk);
            w /= dkk;
          } else if (info == 0) {
            info = k + 1;
          }
        }
        k += 1;
      } else {
        transpositions.coeffRef(k) = StorageIndex(k);
        transpositions.coeffRef(k + 1) = StorageIndex(kp);

        const RealScalar d11 = numext::real(mat.coeff(k, k));
        const RealScalar d22 = numext::real(mat.coeff(k + 1, k + 1));
        const Scalar d21 = mat.coeff(k + 1, k);
        mat.coeffRef(k, k) = Scalar(d11);
        mat.coeffRef(k + 1, k + 1) = Scalar(d22);
        // Scaled 2x2 inverse (LAPACK xSYTF2/xHETF2 strategy). NEVER form det = d11*d22 - |d21|^2 or
        // abs2(d21) directly: those over/underflow for well-conditioned but extreme-scaled blocks
        // (e.g. [[0,s],[s,0]], s=1e200, where det = -s^2 overflows/underflows). Instead divide
        // through by the off-diagonal d21, so the scaled determinant
        //   denom = real(ak*akm1) - 1 = det / |d21|^2     (with ak = d22/d21, akm1 = d11/conj(d21))
        // stays O(1). The reciprocals MUST use Eigen's overflow-safe Scalar division.
        const Scalar id = Scalar(1) / d21;
        const Scalar icjd = numext::conj(id);  // 1 / conj(d21)
        const Scalar ak = d22 * id;
        const Scalar akm1 = d11 * icjd;
        const RealScalar denom = numext::real(ak * akm1) - RealScalar(1);
        // A non-finite 2x2 block (e.g. a NaN pulled in from a candidate row/column) is a numerical
        // failure; flag it so it is reported rather than silently propagated.
        if (info == 0 && (numext::isnan)(denom)) info = k + 1;

        const Index rs = n - k - 2;
        if (rs > 0) {
          const RealScalar t = RealScalar(1) / denom;
          auto A22 = mat.block(k + 2, k + 2, rs, rs);
          auto c0 = mat.col(k).tail(rs);
          auto c1 = mat.col(k + 1).tail(rs);
          // Store the unit lower factor columns L = U D^{-1} (U = [c0 c1]) in scaled form:
          //   L_k = t*(ak*u0 - u1)*icjd,  L_{k+1} = t*(akm1*u1 - u0)*id  (= the unscaled rows of U D^{-1}).
          for (Index i = 0; i < rs; ++i) {
            const Scalar u0 = c0.coeff(i);
            const Scalar u1 = c1.coeff(i);
            c0.coeffRef(i) = t * (ak * u0 - u1) * icjd;
            c1.coeffRef(i) = t * (akm1 * u1 - u0) * id;
          }
          // Trailing update A22 <- A22 - L D L^*  (== A22 - U D^{-1} U^*, since L = U D^{-1}), using the
          // just-stored scaled L columns and the ORIGINAL block entries d11,d22,d21 as coefficients --
          // so no 1/det factor (which would over/underflow) appears. Two rank-1 (syr) and one rank-2
          // (syr2) self-adjoint updates.
          A22.template selfadjointView<Lower>().rankUpdate(c0, -d11);
          A22.template selfadjointView<Lower>().rankUpdate(c1, -d22);
          A22.template selfadjointView<Lower>().rankUpdate(c0, c1, -numext::conj(d21));
        }
        // Move the 2x2 off-diagonal of D out of the L storage.
        subdiag.coeffRef(k) = d21;
        subdiag.coeffRef(k + 1) = Scalar(0);
        mat.coeffRef(k + 1, k) = Scalar(0);
        k += 2;
      }
    }
    return info;
  }

  // Partial factorization of a panel of at most `nb` columns of the lower triangle of `mat`, starting
  // at column k0, using the Bunch-Kaufman method (level-2 within the panel). Follows the panel/workspace
  // structure of LAPACK's xLASYF. The trailing sub-matrix mat(k0+kb:n, k0+kb:n) is left untouched; the
  // workspace `W` (n x nb) returns, in its first kb columns, the product L21*D restricted to the panel
  // columns, so that the caller can apply the deferred trailing update A22 <- A22 - L21 * (L21*D)^* with
  // a single level-3 (triangular) matrix product. Returns the number kb (<= nb) of columns factorized.
  template <typename MatrixType, typename WorkspaceType, typename TranspositionType, typename SubDiagType>
  static Index partial_factor(MatrixType& mat, Index k0, Index nb, WorkspaceType& W, TranspositionType& transpositions,
                              SubDiagType& subdiag, Index& info) {
    using numext::abs;
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::RealScalar RealScalar;
    typedef typename TranspositionType::StorageIndex StorageIndex;
    const Index n = mat.rows();
    const RealScalar alpha = bunch_kaufman_alpha<RealScalar>();
    const bool is_complex = NumTraits<Scalar>::IsComplex;

    Index j = 0;
    while (j < nb) {
      const Index jc = k0 + j;
      const Index h = n - jc;

      // W(jc:n, j) <- updated column jc = (column jc of A) minus the contributions of the panel columns
      // [k0, jc) already factorized (those are stored as L in mat, and L*D in W).
      W.col(j).segment(jc, h) = mat.col(jc).segment(jc, h);
      if (is_complex) W.coeffRef(jc, j) = Scalar(numext::real(W.coeff(jc, j)));
      if (j > 0) {
        W.col(j).segment(jc, h).noalias() -= mat.block(jc, k0, h, j) * W.row(jc).head(j).adjoint();
        if (is_complex) W.coeffRef(jc, j) = Scalar(numext::real(W.coeff(jc, j)));
      }

      Index kstep = 1;
      Index kp = jc;
      const RealScalar absakk = abs(numext::real(W.coeff(jc, j)));
      Index imax = jc;
      RealScalar colmax(0);
      if (jc + 1 < n) {
        Index rel = 0;
        colmax = W.col(j).segment(jc + 1, n - jc - 1).cwiseAbs().maxCoeff(&rel);
        imax = jc + 1 + rel;
      }

      if (numext::is_exactly_zero((numext::maxi)(absakk, colmax)) || (numext::isnan)(absakk)) {
        kp = jc;
        if (info == 0) info = jc + 1;
      } else if (absakk >= alpha * colmax) {
        kp = jc;
      } else {
        // W(jc:n, j+1) <- updated column imax (its leading row part is read from row imax of mat and
        // conjugated, the trailing part from column imax).
        const Index hr = imax - jc;
        if (hr > 0) W.col(j + 1).segment(jc, hr) = mat.row(imax).segment(jc, hr).adjoint();
        W.col(j + 1).segment(imax, n - imax) = mat.col(imax).segment(imax, n - imax);
        if (is_complex) W.coeffRef(imax, j + 1) = Scalar(numext::real(W.coeff(imax, j + 1)));
        if (j > 0) {
          W.col(j + 1).segment(jc, n - jc).noalias() -= mat.block(jc, k0, n - jc, j) * W.row(imax).head(j).adjoint();
          if (is_complex) W.coeffRef(imax, j + 1) = Scalar(numext::real(W.coeff(imax, j + 1)));
        }

        RealScalar rowmax(0);
        if (imax > jc) rowmax = W.col(j + 1).segment(jc, imax - jc).cwiseAbs().maxCoeff();
        if (imax + 1 < n) {
          RealScalar rowmax2 = W.col(j + 1).segment(imax + 1, n - imax - 1).cwiseAbs().maxCoeff();
          rowmax = (numext::maxi)(rowmax, rowmax2);
        }
        if (absakk >= alpha * colmax * (colmax / rowmax)) {
          kp = jc;
        } else if (abs(numext::real(W.coeff(imax, j + 1))) >= alpha * rowmax) {
          kp = imax;
          // imax becomes the 1x1 pivot column: its updated column (in W(:,j+1)) replaces W(:,j).
          W.col(j).segment(jc, n - jc) = W.col(j + 1).segment(jc, n - jc);
        } else {
          kp = imax;
          kstep = 2;
        }
      }

      // A 2x2 block must fit in the panel; otherwise defer this column to the next panel.
      if (kstep == 2 && j + 1 >= nb) break;

      const Index kk = jc + kstep - 1;
      if (kp != kk) {
        apply_symmetric_pivot(mat, jc, kk, kp, kstep);
        const Index nwc = j + kstep;  // number of populated W columns
        W.row(kk).head(nwc).swap(W.row(kp).head(nwc));
      }

      if (kstep == 1) {
        transpositions.coeffRef(jc) = StorageIndex(kp);
        subdiag.coeffRef(jc) = Scalar(0);
        const RealScalar dval = numext::real(W.coeff(jc, j));
        mat.coeffRef(jc, jc) = Scalar(dval);
        const Index rs = n - jc - 1;
        if (rs > 0) {
          if (!numext::is_exactly_zero(dval)) {
            mat.col(jc).tail(rs) = W.col(j).segment(jc + 1, rs) / dval;
          } else {
            // Zero pivot (singular): the Schur column is zero too. Store an explicit zero L column so
            // that matrixL() is well formed, matching the unblocked path (which leaves the already
            // zeroed Schur column in place).
            mat.col(jc).tail(rs).setZero();
            if (info == 0) info = jc + 1;
          }
        }
        j += 1;
      } else {
        transpositions.coeffRef(jc) = StorageIndex(jc);
        transpositions.coeffRef(jc + 1) = StorageIndex(kp);
        const RealScalar d11 = numext::real(W.coeff(jc, j));
        const Scalar d21 = W.coeff(jc + 1, j);
        const RealScalar d22 = numext::real(W.coeff(jc + 1, j + 1));
        mat.coeffRef(jc, jc) = Scalar(d11);
        mat.coeffRef(jc + 1, jc + 1) = Scalar(d22);
        // Scaled 2x2 inverse (see unblocked()): divide through by d21 so the scaled determinant
        // denom = det/|d21|^2 stays O(1); det = d11*d22 - |d21|^2 and abs2(d21) are never formed (they
        // over/underflow on extreme-scaled blocks). The deferred level-3 trailing update below uses W
        // (= L*D, original scale), so it carries no 1/det factor either.
        const Scalar id = Scalar(1) / d21;
        const Scalar icjd = numext::conj(id);  // 1 / conj(d21)
        const Scalar ak = d22 * id;
        const Scalar akm1 = d11 * icjd;
        const RealScalar denom = numext::real(ak * akm1) - RealScalar(1);
        if (info == 0 && (numext::isnan)(denom)) info = jc + 1;
        const Index rs = n - jc - 2;
        if (rs > 0) {
          // L(jc+2:n, jc:jc+1) = W(jc+2:n, j:j+1) * D^{-1}, as vectorized column expressions:
          //   L_k = (t*icjd)*(ak*w0 - w1),  L_{k+1} = (t*id)*(akm1*w1 - w0).
          const RealScalar t = RealScalar(1) / denom;
          const Scalar tic = t * icjd;
          const Scalar tid = t * id;
          auto w0 = W.col(j).segment(jc + 2, rs);
          auto w1 = W.col(j + 1).segment(jc + 2, rs);
          mat.col(jc).tail(rs) = tic * (ak * w0 - w1);
          mat.col(jc + 1).tail(rs) = tid * (akm1 * w1 - w0);
        }
        subdiag.coeffRef(jc) = d21;
        subdiag.coeffRef(jc + 1) = Scalar(0);
        mat.coeffRef(jc + 1, jc) = Scalar(0);
        j += 2;
      }
    }
    return j;
  }

  // Blocked (level-3 BLAS) Bunch-Kaufman factorization of the lower triangle of `mat`, in place.
  template <typename MatrixType, typename TranspositionType, typename SubDiagType, typename WorkspaceType>
  static Index blocked(MatrixType& mat, TranspositionType& transpositions, SubDiagType& subdiag,
                       WorkspaceType& workspace) {
    const Index n = mat.rows();
    const Index nb = bunch_kaufman_blocksize<typename MatrixType::Scalar>();
    if (nb < 2 || n <= nb) return unblocked(mat, transpositions, subdiag, 0);

    // One extra workspace column holds the candidate ("imax") column examined during pivot selection.
    workspace.resize(n, nb + 1);
    Index info = 0;
    Index k = 0;
    while (k < n) {
      if (n - k > nb) {
        const Index kb = partial_factor(mat, k, nb, workspace, transpositions, subdiag, info);
        const Index rs = n - k - kb;
        if (rs > 0) {
          // Deferred trailing update of the lower triangle: A22 <- A22 - L21 * (L21*D)^*.
          mat.block(k + kb, k + kb, rs, rs).template triangularView<Lower>() -=
              mat.block(k + kb, k, rs, kb) * workspace.block(k + kb, 0, rs, kb).adjoint();
        }
        if (kb == 0) {  // defensive: avoid an infinite loop (should not happen for nb >= 2)
          info = (info == 0) ? unblocked(mat, transpositions, subdiag, k) : info;
          break;
        }
        k += kb;
      } else {
        const Index info2 = unblocked(mat, transpositions, subdiag, k);
        if (info == 0) info = info2;
        break;
      }
    }
    return info;
  }
};

template <>
struct bunch_kaufman<Upper> {
  template <typename MatrixType, typename TranspositionType, typename SubDiagType>
  static EIGEN_STRONG_INLINE Index unblocked(MatrixType& mat, TranspositionType& transpositions, SubDiagType& subdiag,
                                             Index k0 = 0) {
    Transpose<MatrixType> matt(mat);
    return bunch_kaufman<Lower>::unblocked(matt, transpositions, subdiag, k0);
  }

  template <typename MatrixType, typename TranspositionType, typename SubDiagType, typename WorkspaceType>
  static EIGEN_STRONG_INLINE Index blocked(MatrixType& mat, TranspositionType& transpositions, SubDiagType& subdiag,
                                           WorkspaceType& workspace) {
    Transpose<MatrixType> matt(mat);
    return bunch_kaufman<Lower>::blocked(matt, transpositions, subdiag, workspace);
  }
};

template <typename MatrixType>
struct BunchKaufman_Traits<MatrixType, Lower> {
  typedef const TriangularView<const MatrixType, UnitLower> MatrixL;
  typedef const TriangularView<const typename MatrixType::AdjointReturnType, UnitUpper> MatrixU;
  static inline MatrixL getL(const MatrixType& m) { return MatrixL(m); }
  static inline MatrixU getU(const MatrixType& m) { return MatrixU(m.adjoint()); }
};

template <typename MatrixType>
struct BunchKaufman_Traits<MatrixType, Upper> {
  typedef const TriangularView<const typename MatrixType::AdjointReturnType, UnitLower> MatrixL;
  typedef const TriangularView<const MatrixType, UnitUpper> MatrixU;
  static inline MatrixL getL(const MatrixType& m) { return MatrixL(m.adjoint()); }
  static inline MatrixU getU(const MatrixType& m) { return MatrixU(m); }
};

}  // end namespace internal

template <typename MatrixType, int UpLo_>
template <typename Derived>
void BunchKaufman<MatrixType, UpLo_>::applyD(MatrixBase<Derived>& x) const {
  const Index n = m_matrix.rows();
  Index k = 0;
  while (k < n) {
    if (k + 1 < n && !numext::is_exactly_zero(m_subdiag.coeff(k))) {
      const RealScalar d11 = numext::real(m_matrix.coeff(k, k));
      const RealScalar d22 = numext::real(m_matrix.coeff(k + 1, k + 1));
      const Scalar d21 = m_subdiag.coeff(k);
      for (Index j = 0; j < x.cols(); ++j) {
        const Scalar x0 = x.coeff(k, j);
        const Scalar x1 = x.coeff(k + 1, j);
        x.coeffRef(k, j) = d11 * x0 + numext::conj(d21) * x1;
        x.coeffRef(k + 1, j) = d21 * x0 + d22 * x1;
      }
      k += 2;
    } else {
      x.row(k) *= numext::real(m_matrix.coeff(k, k));
      k += 1;
    }
  }
}

template <typename MatrixType, int UpLo_>
template <bool Conjugate, typename Derived>
void BunchKaufman<MatrixType, UpLo_>::solveInPlaceD(MatrixBase<Derived>& x) const {
  using numext::abs;
  const Index n = m_matrix.rows();
  // Use the pseudo-inverse of singular 1x1 blocks (see Eigen bug 241 for LDLT). The 2x2 blocks
  // produced by Bunch-Kaufman pivoting are non-singular by construction.
  const RealScalar tol = (std::numeric_limits<RealScalar>::min)();
  Index k = 0;
  while (k < n) {
    if (k + 1 < n && !numext::is_exactly_zero(m_subdiag.coeff(k))) {
      const RealScalar d11 = numext::real(m_matrix.coeff(k, k));
      const RealScalar d22 = numext::real(m_matrix.coeff(k + 1, k + 1));
      // D = [ d11  conj(d21) ; d21  d22 ]; for the transpose solve use conj(d21) instead of d21.
      const Scalar d21 = Conjugate ? m_subdiag.coeff(k) : numext::conj(m_subdiag.coeff(k));
      // Scaled 2x2 solve (LAPACK xSYTRS/xHETRS): divide through by d21 so the scaled determinant
      // denom = det/|d21|^2 is O(1); det = d11*d22 - |d21|^2 is never formed (it over/underflows on
      // extreme-scaled blocks, e.g. [[0,s],[s,0]], s=1e+-200). Reciprocals use overflow-safe division.
      const Scalar id = Scalar(1) / d21;
      const Scalar icjd = numext::conj(id);  // 1 / conj(d21)
      const Scalar ak = d22 * id;
      const Scalar akm1 = d11 * icjd;
      const RealScalar t = RealScalar(1) / (numext::real(ak * akm1) - RealScalar(1));
      for (Index j = 0; j < x.cols(); ++j) {
        const Scalar x0 = x.coeff(k, j);
        const Scalar x1 = x.coeff(k + 1, j);
        const Scalar bk = x1 * id;
        const Scalar bkm1 = x0 * icjd;
        x.coeffRef(k, j) = t * (ak * bkm1 - bk);
        x.coeffRef(k + 1, j) = t * (akm1 * bk - bkm1);
      }
      k += 2;
    } else {
      const RealScalar dk = numext::real(m_matrix.coeff(k, k));
      if (abs(dk) > tol)
        x.row(k) /= dk;
      else
        x.row(k).setZero();
      k += 1;
    }
  }
}

template <typename MatrixType, int UpLo_>
void BunchKaufman<MatrixType, UpLo_>::computeInertia() {
  const Index n = m_matrix.rows();
  m_n_pos = m_n_neg = m_n_zero = 0;
  Index k = 0;
  while (k < n) {
    if (k + 1 < n && !numext::is_exactly_zero(m_subdiag.coeff(k))) {
      const RealScalar d11 = numext::real(m_matrix.coeff(k, k));
      const RealScalar d22 = numext::real(m_matrix.coeff(k + 1, k + 1));
      const Scalar d21 = m_subdiag.coeff(k);
      // Scaled determinant denom = det/|d21|^2 (|d21|^2 > 0 for a 2x2 block), so sign(denom) == sign(det);
      // avoids forming det = d11*d22 - |d21|^2, which over/underflows on extreme-scaled 2x2 blocks.
      const Scalar id = Scalar(1) / d21;
      const RealScalar denom = numext::real((d22 * id) * (d11 * numext::conj(id))) - RealScalar(1);
      if (denom < RealScalar(0)) {
        // Indefinite 2x2 block: one positive and one negative eigenvalue.
        ++m_n_pos;
        ++m_n_neg;
      } else if (numext::is_exactly_zero(denom)) {
        const RealScalar tr = d11 + d22;
        if (tr > RealScalar(0))
          ++m_n_pos;
        else if (tr < RealScalar(0))
          ++m_n_neg;
        else
          ++m_n_zero;
        ++m_n_zero;
      } else {
        // denom > 0: both eigenvalues share the sign of the trace.
        if (d11 + d22 > RealScalar(0))
          m_n_pos += 2;
        else
          m_n_neg += 2;
      }
      k += 2;
    } else {
      const RealScalar dk = numext::real(m_matrix.coeff(k, k));
      if (dk > RealScalar(0))
        ++m_n_pos;
      else if (dk < RealScalar(0))
        ++m_n_neg;
      else
        ++m_n_zero;
      k += 1;
    }
  }
}

/** Compute / recompute the Bunch-Kaufman factorization \f$ P A P^T = L D L^* \f$ (or \f$ U^* D U \f$ for the
 * upper variant) of \a a. */
template <typename MatrixType, int UpLo_>
template <typename InputType>
BunchKaufman<MatrixType, UpLo_>& BunchKaufman<MatrixType, UpLo_>::compute(const EigenBase<InputType>& a) {
  eigen_assert(a.rows() == a.cols());
  const Index size = a.rows();

  m_matrix = a.derived();

  // L1 norm of the implicit self-adjoint matrix, for rcond().
  m_l1_norm = m_matrix.template selfadjointView<UpLo_>().l1Norm();

  m_transpositions.resize(size);
  m_subdiag.resize(size);
  m_isInitialized = false;

  Index info = internal::bunch_kaufman<UpLo>::blocked(m_matrix, m_transpositions, m_subdiag, m_workspace);
  m_info = (info == 0) ? Success : NumericalIssue;

  // The upper variant factorizes the transpose of the (Hermitian) matrix, i.e. its complex conjugate,
  // so the recorded 2x2 off-diagonals of D come out conjugated; undo that so that m_subdiag holds the
  // true sub-diagonal D(k+1,k) consumed by applyD()/solveInPlaceD()/computeInertia(). (No-op when real.)
  EIGEN_IF_CONSTEXPR (int(UpLo) == int(Upper) && NumTraits<Scalar>::IsComplex) {
    m_subdiag = m_subdiag.conjugate();
  }

  computeInertia();

  m_isInitialized = true;
  return *this;
}

#ifndef EIGEN_PARSED_BY_DOXYGEN
template <typename MatrixType_, int UpLo_>
template <typename RhsType, typename DstType>
void BunchKaufman<MatrixType_, UpLo_>::_solve_impl(const RhsType& rhs, DstType& dst) const {
  _solve_impl_transposed<true>(rhs, dst);
}

template <typename MatrixType_, int UpLo_>
template <bool Conjugate, typename RhsType, typename DstType>
void BunchKaufman<MatrixType_, UpLo_>::_solve_impl_transposed(const RhsType& rhs, DstType& dst) const {
  // A^{-1} b = P^T L^{-*} D^{-1} L^{-1} P b   (and the conjugated variants for transpose / adjoint solves).
  dst = m_transpositions * rhs;
  matrixL().template conjugateIf<!Conjugate>().solveInPlace(dst);
  solveInPlaceD<Conjugate>(dst);
  matrixL().transpose().template conjugateIf<Conjugate>().solveInPlace(dst);
  dst = m_transpositions.transpose() * dst;
}
#endif

/** \internal use x = bunchKaufman_object.solve(x);
 *
 * This is the \em in-place version of solve().
 *
 * \returns true always.
 */
template <typename MatrixType, int UpLo_>
template <typename Derived>
bool BunchKaufman<MatrixType, UpLo_>::solveInPlace(MatrixBase<Derived>& bAndX) const {
  eigen_assert(m_isInitialized && "BunchKaufman is not initialized.");
  eigen_assert(m_matrix.rows() == bAndX.rows());
  bAndX = this->solve(bAndX);
  return true;
}

/** \returns the matrix represented by the decomposition, i.e., the product \f$ P^T L D L^* P \f$.
 * This function is provided for debug purposes. */
template <typename MatrixType, int UpLo_>
MatrixType BunchKaufman<MatrixType, UpLo_>::reconstructedMatrix() const {
  eigen_assert(m_isInitialized && "BunchKaufman is not initialized.");
  const Index size = m_matrix.rows();
  MatrixType res(size, size);

  res.setIdentity();
  res = transpositionsP() * res;              // P
  res = matrixU() * res;                      // U P = L^* P
  applyD(res);                                // D L^* P
  res = matrixL() * res;                      // L D L^* P
  res = transpositionsP().transpose() * res;  // P^T L D L^* P

  return res;
}

/** \cholesky_module
 * \returns the Bunch-Kaufman factorization of \c *this
 * \sa SelfAdjointView::bunchKaufman()
 */
template <typename MatrixType, unsigned int UpLo>
inline BunchKaufman<typename SelfAdjointView<MatrixType, UpLo>::PlainObject, UpLo>
SelfAdjointView<MatrixType, UpLo>::bunchKaufman() const {
  return BunchKaufman<PlainObject, UpLo>(m_matrix);
}

/** \cholesky_module
 * \returns the Bunch-Kaufman factorization of \c *this
 * \sa SelfAdjointView::bunchKaufman()
 */
template <typename Derived>
inline BunchKaufman<typename MatrixBase<Derived>::PlainObject> MatrixBase<Derived>::bunchKaufman() const {
  return BunchKaufman<PlainObject>(derived());
}

}  // end namespace Eigen

#endif  // EIGEN_BUNCHKAUFMAN_H
