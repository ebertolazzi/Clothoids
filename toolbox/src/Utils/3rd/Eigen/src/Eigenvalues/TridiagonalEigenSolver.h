// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Rasmus Munk Larsen <rmlarsen@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

#ifndef EIGEN_TRIDIAGONAL_EIGENSOLVER_H
#define EIGEN_TRIDIAGONAL_EIGENSOLVER_H

#include "./TridiagonalBisection.h"
#include "./TridiagonalInverseIteration.h"

// IWYU pragma: private
#include "./InternalHeaderCheck.h"

namespace Eigen {

/** \eigenvalues_module \ingroup Eigenvalues_Module
 *
 * \class TridiagonalEigenSolver
 *
 * \brief Computes eigenvalues and eigenvectors of a real symmetric tridiagonal matrix
 *
 * \tparam Scalar_ the (real floating-point) scalar type of the matrix, e.g. \c float or \c double.
 *
 * This solver computes the eigenvalues of a real symmetric tridiagonal matrix \f$ T \f$, given by
 * its diagonal and sub-diagonal, using SIMD-accelerated, multi-threaded Sturm-sequence spectral
 * bisection (cf. LAPACK's \c xSTEBZ), and the corresponding eigenvectors by inverse iteration
 * (cf. LAPACK's \c xSTEIN). Unlike the implicit-QR path of
 * SelfAdjointEigenSolver::computeFromTridiagonal(), it can compute an arbitrary contiguous subset
 * of the spectrum -- by index range or by value range -- selected with an EigenvalueRange, and it
 * never pays for eigenvectors that are not requested:
 *
 * \code
 * TridiagonalEigenSolver<double> es;
 * es.compute(diag, subdiag);                            // full spectrum, eigenvalues + eigenvectors
 * es.compute(diag, subdiag, EigenvaluesOnly);           // full spectrum, eigenvalues only
 * es.compute(diag, subdiag, ComputeEigenvectors,
 *            EigenvalueRange::indices(0, 10));          // the 10 smallest eigenpairs
 *
 * // Staged: eigenvalues now, eigenvectors later (only if they turn out to be needed).
 * es.computeEigenvalues(diag, subdiag, EigenvalueRange::values(vl, vu));
 * es.computeEigenvectors();
 *
 * // Direct: eigenvectors for already-known eigenvalues.
 * es.computeEigenvectors(diag, subdiag, eigenvalues);
 * \endcode
 *
 * In every mode the computed eigenvalues are returned by eigenvalues() in non-decreasing order,
 * one entry per selected eigenvalue, and eigenvectors() is the \c n x \c m matrix whose column
 * \c j is a unit-norm eigenvector for \c eigenvalues()(j).
 *
 * Besides subset selection, bisection is typically more accurate than the QR algorithm: it
 * resolves each eigenvalue to about one unit in the last place of \f$ \|T\| \f$ independently of
 * the others, whereas the QR forward error accumulates through its O(n) sweeps and grows with the
 * matrix size. (Both are backward stable; this is absolute accuracy relative to \f$ \|T\| \f$, not
 * high relative accuracy for eigenvalues much smaller than \f$ \|T\| \f$.)
 *
 * The eigenvectors are those of the tridiagonal \f$ T \f$ itself; recovering eigenvectors of a
 * dense matrix additionally requires the Householder back-transform performed by
 * SelfAdjointEigenSolver::compute().
 *
 * \note The eigenvalues of a complex Hermitian tridiagonal matrix depend only on the moduli of its
 * off-diagonal entries: it is unitarily similar (via a diagonal phase matrix) to the real symmetric
 * tridiagonal with the same diagonal and off-diagonals \f$ |\beta_k| \f$. So to compute them, pass
 * \c subdiag.cwiseAbs() as the real sub-diagonal. The eigenvectors, by contrast, differ from the
 * real ones by that diagonal phase and cannot be recovered this way.
 *
 * \sa EigenvalueRange, SelfAdjointEigenSolver::computeFromTridiagonal(), class Tridiagonalization
 */
template <typename Scalar_>
class TridiagonalEigenSolver {
 public:
  /** \brief Scalar type of the matrix; must be real. */
  typedef Scalar_ Scalar;
  typedef Scalar RealScalar;
  static_assert(NumTraits<Scalar>::IsComplex == 0 && NumTraits<Scalar>::IsInteger == 0,
                "TridiagonalEigenSolver requires a real floating-point scalar type; for the eigenvalues of a complex "
                "Hermitian tridiagonal matrix, pass subdiag.cwiseAbs() as the real sub-diagonal (see the class "
                "documentation).");

  /** \brief Type for the eigenvalues and the input diagonals: a dynamic-size column vector. */
  typedef Matrix<Scalar, Dynamic, 1> VectorType;
  /** \brief Type for the eigenvector matrix: dynamic-size, one column per selected eigenvalue. */
  typedef Matrix<Scalar, Dynamic, Dynamic> MatrixType;

  /** \brief Default constructor. Call compute() before querying any result. */
  TridiagonalEigenSolver() = default;

  /** \brief Constructor pre-allocating room for the full spectrum of a matrix of dimension \a size. */
  explicit TridiagonalEigenSolver(Index size)
      : m_eivalues(size), m_eivec(size, size), m_diag(size), m_subdiag(size > 1 ? size - 1 : 0) {}

  /** \brief Constructor; computes the eigendecomposition of the given tridiagonal matrix.
   *
   * Equivalent to default construction followed by compute(\a diag, \a subdiag, \a options, \a range).
   */
  template <typename DiagType, typename SubdiagType>
  TridiagonalEigenSolver(const MatrixBase<DiagType>& diag, const MatrixBase<SubdiagType>& subdiag,
                         int options = ComputeEigenvectors, const EigenvalueRange& range = EigenvalueRange::all())
      : TridiagonalEigenSolver() {
    compute(diag, subdiag, options, range);
  }

  /** \brief Computes the selected eigenvalues, and optionally eigenvectors, of a real symmetric
   * tridiagonal matrix.
   *
   * \param[in] diag    The diagonal of the matrix \f$ T \f$ (length \c n).
   * \param[in] subdiag The sub-diagonal of \f$ T \f$ (length \c n-1).
   * \param[in] options Either #ComputeEigenvectors (the default) or #EigenvaluesOnly.
   * \param[in] range   Which eigenvalues to compute (see EigenvalueRange); defaults to the whole
   *            spectrum.
   * \returns Reference to \c *this
   *
   * Equivalent to computeEigenvalues(\a diag, \a subdiag, \a range), followed -- when \a options
   * is #ComputeEigenvectors and the eigenvalues were computed successfully -- by
   * computeEigenvectors().
   */
  template <typename DiagType, typename SubdiagType>
  TridiagonalEigenSolver& compute(const MatrixBase<DiagType>& diag, const MatrixBase<SubdiagType>& subdiag,
                                  int options = ComputeEigenvectors,
                                  const EigenvalueRange& range = EigenvalueRange::all()) {
    eigen_assert((options & ~EigVecMask) == 0 && (options & EigVecMask) != EigVecMask && "invalid option parameter");
    computeEigenvalues(diag, subdiag, range);
    if (m_info == Success && (options & ComputeEigenvectors) == ComputeEigenvectors) computeEigenvectors();
    return *this;
  }

  /** \brief Computes the selected eigenvalues (only) of a real symmetric tridiagonal matrix.
   *
   * \param[in] diag    The diagonal of the matrix \f$ T \f$ (length \c n).
   * \param[in] subdiag The sub-diagonal of \f$ T \f$ (length \c n-1).
   * \param[in] range   Which eigenvalues to compute (see EigenvalueRange); defaults to the whole
   *            spectrum.
   * \returns Reference to \c *this
   *
   * After the call, eigenvalues() returns the selected eigenvalues in non-decreasing order (one
   * entry per selected eigenvalue), and info() reports \c Success, or \c NoConvergence if the
   * input contains a non-finite entry. The tridiagonal is retained, so the eigenvectors can be
   * obtained afterwards with computeEigenvectors() -- and need never be computed if they turn out
   * not to be required.
   */
  template <typename DiagType, typename SubdiagType>
  TridiagonalEigenSolver& computeEigenvalues(const MatrixBase<DiagType>& diag, const MatrixBase<SubdiagType>& subdiag,
                                             const EigenvalueRange& range = EigenvalueRange::all());

  /** \brief Computes eigenvectors for the eigenvalues of the preceding compute step, by inverse
   * iteration.
   *
   * \returns Reference to \c *this
   *
   * Completes a staged solve: call it after computeEigenvalues() to compute the eigenvectors for
   * the eigenvalues just found, using the retained tridiagonal. See
   * computeEigenvectors(const MatrixBase<DiagType>&, const MatrixBase<SubdiagType>&, const
   * MatrixBase<EivalsType>&) for the algorithm and its guarantees.
   *
   * \pre A successful computeEigenvalues() (or compute()) call.
   */
  TridiagonalEigenSolver& computeEigenvectors() {
    eigen_assert(m_isInitialized && "TridiagonalEigenSolver is not initialized.");
    eigen_assert(m_info == Success && "computeEigenvectors() requires a preceding successful eigenvalue computation");
    computeEigenvectorsImpl();
    return *this;
  }

  /** \brief Computes eigenvectors of a tridiagonal matrix by inverse iteration from known eigenvalues.
   *
   * \param[in] diag        The diagonal of the symmetric tridiagonal matrix \f$ T \f$ (length \c n).
   * \param[in] subdiag     The sub-diagonal of \f$ T \f$ (length \c n-1).
   * \param[in] eigenvalues The eigenvalues whose eigenvectors are wanted, in non-decreasing order
   *            (e.g. the output of computeEigenvalues()). Length \c m \f$\le\f$ \c n.
   * \returns Reference to \c *this
   *
   * Computes an eigenvector of \f$ T \f$ for each supplied eigenvalue by inverse iteration (LAPACK's
   * xSTEIN, built on the xLAGTF / xLAGTS factorization and overflow-safe solve of the deliberately
   * near-singular \f$ T - \lambda I \f$), reorthogonalizing within tight clusters so that a degenerate
   * cluster yields an orthonormal basis. After the call, eigenvectors() returns the \c n x \c m matrix
   * whose column \c j is a unit-norm eigenvector for \c eigenvalues[j], and eigenvalues() returns the
   * supplied eigenvalues.
   *
   * \b Subsets. This is the natural companion to the subset-selecting eigenvalue path (see
   * EigenvalueRange): pass that subset as \a eigenvalues and you get back exactly those \c m
   * eigenvectors -- the result has one column per supplied eigenvalue and costs one
   * inverse-iteration solve each, so you never pay for vectors you did not ask for.
   *
   * \note The eigenvalues must be sorted in non-decreasing order, as the cluster
   * reorthogonalization relies on it.
   *
   * \note If any eigenvector fails to converge within the internal inverse-iteration step limit, the
   * vectors are still returned but info() reports \c NoConvergence (cf. LAPACK xSTEIN).
   *
   * \note The inverse iteration is bit-identical for any number of threads. The eigenvectors of a
   * numerically degenerate cluster may nonetheless differ in their last bits between a single- and a
   * multi-threaded run, because the cluster refinement uses parallel (non-associative) matrix products;
   * each result is an equally valid orthonormal basis of the same invariant subspace.
   *
   * \warning Cluster reorthogonalization is performed only \e among the supplied eigenvalues. The
   * returned eigenvectors are always mutually orthonormal, but they are not orthogonalized against
   * eigenvectors of a degenerate cluster that lie outside \a eigenvalues (those are never computed). If
   * a requested subset slices through the middle of a numerically degenerate cluster, supply the whole
   * cluster (widen the range) to obtain the correct invariant subspace. For well-separated eigenvalues,
   * or when the subset covers each cluster entirely, this never matters.
   *
   * \sa computeEigenvalues(), eigenvectors()
   */
  template <typename DiagType, typename SubdiagType, typename EivalsType>
  TridiagonalEigenSolver& computeEigenvectors(const MatrixBase<DiagType>& diag, const MatrixBase<SubdiagType>& subdiag,
                                              const MatrixBase<EivalsType>& eigenvalues);

  /** \brief Returns the computed eigenvalues, in non-decreasing order.
   *
   * \pre compute(), computeEigenvalues() or computeEigenvectors() has been called.
   *
   * The returned vector has one entry per \e selected eigenvalue, so a subset request yields a
   * vector shorter than the matrix dimension. Eigenvalues are repeated according to their
   * algebraic multiplicity.
   */
  const VectorType& eigenvalues() const {
    eigen_assert(m_isInitialized && "TridiagonalEigenSolver is not initialized.");
    return m_eivalues;
  }

  /** \brief Returns the computed eigenvectors, one column per eigenvalue.
   *
   * \pre The eigenvectors have been computed (compute() with #ComputeEigenvectors, or one of the
   * computeEigenvectors() overloads).
   *
   * Column \c j of the returned \c n x \c m matrix is a unit-norm eigenvector of \f$ T \f$ for
   * eigenvalues()(j); the columns are mutually orthonormal.
   */
  const MatrixType& eigenvectors() const {
    eigen_assert(m_isInitialized && "TridiagonalEigenSolver is not initialized.");
    eigen_assert(m_eigenvectorsOk && "The eigenvectors have not been computed.");
    return m_eivec;
  }

  /** \brief Reports whether the computation was successful.
   *
   * \returns \c Success if the computation was successful, \c NoConvergence otherwise.
   */
  ComputationInfo info() const {
    eigen_assert(m_isInitialized && "TridiagonalEigenSolver is not initialized.");
    return m_info;
  }

 protected:
  void computeEigenvectorsImpl();

  VectorType m_eivalues;
  MatrixType m_eivec;
  // The tridiagonal retained by the last compute step, so the staged computeEigenvectors() can
  // re-form T - lambda*I without the caller having to keep diag/subdiag alive.
  VectorType m_diag;
  VectorType m_subdiag;
  ComputationInfo m_info = InvalidInput;
  bool m_isInitialized = false;
  bool m_eigenvectorsOk = false;
};

template <typename Scalar_>
template <typename DiagType, typename SubdiagType>
TridiagonalEigenSolver<Scalar_>& TridiagonalEigenSolver<Scalar_>::computeEigenvalues(
    const MatrixBase<DiagType>& diag, const MatrixBase<SubdiagType>& subdiag, const EigenvalueRange& range) {
  static_assert(internal::is_same<typename DiagType::Scalar, Scalar>::value &&
                    internal::is_same<typename SubdiagType::Scalar, Scalar>::value,
                "diag and subdiag must have the solver's scalar type");
  const Index n = diag.size();
  EIGEN_UNUSED_VARIABLE(n);
  eigen_assert(subdiag.size() == (n > 0 ? n - 1 : 0) && "sub-diagonal must have one fewer entry than the diagonal");

  m_diag = diag;
  m_subdiag = subdiag;
  m_eigenvectorsOk = false;

  // Reject non-finite input up front (cf. SelfAdjointEigenSolver::computeFromTridiagonal()): the
  // Gershgorin bracketing below would otherwise iterate on NaN brackets and return garbage with
  // info() == Success. allFinite() rather than a max-reduction, whose default PropagateFast
  // semantics may not surface a NaN.
  if (!(m_diag.allFinite() && m_subdiag.allFinite())) {
    m_eivalues.resize(0);
    m_info = NoConvergence;
    m_isInitialized = true;
    return *this;
  }

  internal::tridiagonal_bisection(m_diag, m_subdiag, range, RealScalar(0), m_eivalues);
  m_info = Success;
  m_isInitialized = true;
  return *this;
}

template <typename Scalar_>
template <typename DiagType, typename SubdiagType, typename EivalsType>
TridiagonalEigenSolver<Scalar_>& TridiagonalEigenSolver<Scalar_>::computeEigenvectors(
    const MatrixBase<DiagType>& diag, const MatrixBase<SubdiagType>& subdiag,
    const MatrixBase<EivalsType>& eigenvalues) {
  static_assert(internal::is_same<typename DiagType::Scalar, Scalar>::value &&
                    internal::is_same<typename SubdiagType::Scalar, Scalar>::value &&
                    internal::is_same<typename EivalsType::Scalar, Scalar>::value,
                "diag, subdiag and eigenvalues must have the solver's scalar type");
  const Index n = diag.size();
  EIGEN_UNUSED_VARIABLE(n);
  eigen_assert(subdiag.size() == (n > 0 ? n - 1 : 0) && "sub-diagonal must have one fewer entry than the diagonal");
  eigen_assert(eigenvalues.size() <= n && "cannot request more eigenvectors than the size of the matrix");

  m_diag = diag;
  m_subdiag = subdiag;
  m_eivalues = eigenvalues;
  computeEigenvectorsImpl();
  return *this;
}

template <typename Scalar_>
void TridiagonalEigenSolver<Scalar_>::computeEigenvectorsImpl() {
  const Index n = m_diag.size();
  const Index m = m_eivalues.size();
  m_eivec.resize(n, m);
  const Index nonconv = internal::tridiagonal_inverse_iteration(m_diag, m_subdiag, m_eivalues, m_eivec);
  // Refine the eigenvectors of any genuinely degenerate cluster (Rayleigh-Ritz). The eigenvalues are
  // left unchanged, and a non-degenerate spectrum is untouched.
  internal::tridiagonal_rayleigh_ritz_refine(m_diag, m_subdiag, m_eivalues, m_eivec);

  // Like LAPACK xSTEIN, report NoConvergence if any eigenvector failed the inverse-iteration growth
  // test (the columns are still returned, best effort, so eigenvectors() remains usable).
  m_info = (nonconv == 0) ? Success : NoConvergence;
  m_isInitialized = true;
  m_eigenvectorsOk = true;
}

}  // namespace Eigen

#endif  // EIGEN_TRIDIAGONAL_EIGENSOLVER_H
