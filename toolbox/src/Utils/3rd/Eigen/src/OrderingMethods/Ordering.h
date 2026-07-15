// SPDX-License-Identifier: MPL-2.0

// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2012  Désiré Nuentsa-Wakam <desire.nuentsa_wakam@inria.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_ORDERING_H
#define EIGEN_ORDERING_H

// IWYU pragma: private
#include "./InternalHeaderCheck.h"
#include "Eigen_Colamd.h"

namespace Eigen {

/** \ingroup OrderingMethods_Module
 * \class AMDOrdering
 *
 * Functor computing the \em approximate \em minimum \em degree ordering
 * If the matrix is not structurally symmetric, an ordering of A^T+A is computed.
 * Only the sparsity pattern of the input is read — scalar values are not.
 * \tparam  StorageIndex The type of indices of the matrix
 * \sa COLAMDOrdering
 */
template <typename StorageIndex>
class AMDOrdering {
 public:
  typedef PermutationMatrix<Dynamic, Dynamic, StorageIndex> PermutationType;

  /** Compute the permutation vector from a sparse matrix.
   * Only the sparsity pattern of \a mat is read; scalar values are not.
   * This routine is much faster if the input matrix is column-major.
   */
  template <typename MatrixType>
  void operator()(const MatrixType& mat, PermutationType& perm) const {
    // AMD only reads the sparsity pattern. Build a column-major view of mat,
    // then materialize \c pattern(mat + mat^T) directly into a
    // SparseMatrix<signed char> (1-byte placeholder values), bypassing
    // Eigen's generic transpose + sparse-sum evaluators.
    Matrix<StorageIndex, Dynamic, 1> outer_buf;
    Matrix<StorageIndex, Dynamic, 1> inner_buf;
    internal::SparsityPatternRef<StorageIndex> pat = internal::make_col_major_pattern_ref(mat, outer_buf, inner_buf);
    SparseMatrix<signed char, ColMajor, StorageIndex> symm;
    internal::materialize_at_plus_a_pattern(pat, symm);
    internal::minimum_degree_ordering(symm, perm);
  }

  /** Compute the permutation with a selfadjoint matrix.
   * Only the sparsity pattern is used; scalar values are not.
   */
  template <typename SrcType, unsigned int SrcUpLo>
  void operator()(const SparseSelfAdjointView<SrcType, SrcUpLo>& mat, PermutationType& perm) const {
    // Build a column-major pattern view of the underlying matrix and expand
    // its UpLo triangle to the full symmetric pattern in one pass, bypassing
    // Eigen's generic selfadjointView assignment evaluator.
    Matrix<StorageIndex, Dynamic, 1> outer_buf;
    Matrix<StorageIndex, Dynamic, 1> inner_buf;
    internal::SparsityPatternRef<StorageIndex> pat =
        internal::make_col_major_pattern_ref(mat.matrix(), outer_buf, inner_buf);
    SparseMatrix<signed char, ColMajor, StorageIndex> symm;
    internal::materialize_selfadjoint_pattern<SrcUpLo>(pat, symm);
    internal::minimum_degree_ordering(symm, perm);
  }
};

/** \ingroup OrderingMethods_Module
 * \class NaturalOrdering
 *
 * Functor computing the natural ordering (identity)
 *
 * \note Returns an empty permutation matrix
 * \tparam  StorageIndex The type of indices of the matrix
 */
template <typename StorageIndex>
class NaturalOrdering {
 public:
  typedef PermutationMatrix<Dynamic, Dynamic, StorageIndex> PermutationType;

  /** Compute the permutation vector from a column-major sparse matrix */
  template <typename MatrixType>
  void operator()(const MatrixType& /*mat*/, PermutationType& perm) const {
    perm.resize(0);
  }
};

/** \ingroup OrderingMethods_Module
 * \class COLAMDOrdering
 *
 * \tparam  StorageIndex The type of indices of the matrix
 *
 * Functor computing the \em column \em approximate \em minimum \em degree ordering.
 * Only the sparsity pattern of the input is read — scalar values are not.
 */
template <typename StorageIndex>
class COLAMDOrdering {
 public:
  typedef PermutationMatrix<Dynamic, Dynamic, StorageIndex> PermutationType;
  typedef Matrix<StorageIndex, Dynamic, 1> IndexVector;

  /** Compute the permutation vector \a perm from the sparse matrix \a mat. */
  template <typename MatrixType>
  void operator()(const MatrixType& mat, PermutationType& perm) const {
    typedef typename MatrixType::StorageIndex MatrixStorageIndex;
    Matrix<MatrixStorageIndex, Dynamic, 1> outer_buf, inner_buf;
    internal::SparsityPatternRef<MatrixStorageIndex> pat =
        internal::make_col_major_pattern_ref(mat, outer_buf, inner_buf);
    const StorageIndex m = internal::convert_index<StorageIndex>(pat.innerSize);
    const StorageIndex n = internal::convert_index<StorageIndex>(pat.outerSize);
    // Accumulate in Index — Eigen's contract is that any valid nnz fits there
    // (mat.nonZeros() returns Index), so the sum can't overflow. One
    // bounds-checked narrow to StorageIndex at the end catches the only real
    // overflow case (total > StorageIndex range).
    Index total_nnz = 0;
    for (Index j = 0; j < pat.outerSize; ++j) total_nnz += pat.nonZeros(j);
    const StorageIndex nnz = internal::convert_index<StorageIndex>(total_nnz);

    StorageIndex Alen = internal::Colamd::recommended(nnz, m, n);
    double knobs[internal::Colamd::NKnobs];
    StorageIndex stats[internal::Colamd::NStats];
    internal::Colamd::set_defaults(knobs);

    // Colamd writes into A[] in place and needs a contiguous CSC layout, so
    // always compact per column — handles both compressed and uncompressed
    // sources uniformly via SparsityPatternRef::nonZeros(j).
    IndexVector p(n + 1), A(Alen);
    p(0) = 0;
    for (StorageIndex j = 0; j < n; ++j) {
      const Index nz = pat.nonZeros(j);
      const MatrixStorageIndex* src = pat.inner + pat.outer[j];
      copy_colamd_indices(src, nz, A.data() + p(j), std::is_same<MatrixStorageIndex, StorageIndex>());
      p(j + 1) = p(j) + static_cast<StorageIndex>(nz);
    }

    StorageIndex info = internal::Colamd::compute_ordering(m, n, Alen, A.data(), p.data(), knobs, stats);
    EIGEN_UNUSED_VARIABLE(info);
    eigen_assert(info && "COLAMD failed");

    perm.resize(n);
    for (StorageIndex i = 0; i < n; i++) perm.indices()(p(i)) = i;
  }

 private:
  template <typename SrcStorageIndex>
  static void copy_colamd_indices(const SrcStorageIndex* src, Index nz, StorageIndex* dst, std::true_type) {
    std::copy_n(src, nz, dst);
  }

  template <typename SrcStorageIndex>
  static void copy_colamd_indices(const SrcStorageIndex* src, Index nz, StorageIndex* dst, std::false_type) {
    for (Index k = 0; k < nz; ++k) dst[k] = internal::convert_index<StorageIndex>(src[k]);
  }
};

}  // end namespace Eigen

#endif
