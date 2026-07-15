// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2010 Vincent Lejeune
// Copyright (C) 2010 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

#ifndef EIGEN_BLOCK_HOUSEHOLDER_H
#define EIGEN_BLOCK_HOUSEHOLDER_H

// This file contains some helper functions to deal with block householder reflectors

// IWYU pragma: private
#include "./InternalHeaderCheck.h"

namespace Eigen {

namespace internal {

/** \internal */
// This variant avoids modifications in vectors
template <typename TriangularFactorType, typename VectorsType, typename CoeffsType>
void make_block_householder_triangular_factor(TriangularFactorType& triFactor, const VectorsType& vectors,
                                              const CoeffsType& hCoeffs) {
  const Index nbVecs = vectors.cols();
  eigen_assert(triFactor.rows() == nbVecs && triFactor.cols() == nbVecs && vectors.rows() >= nbVecs);

  for (Index i = nbVecs - 1; i >= 0; --i) {
    Index rs = vectors.rows() - i - 1;
    Index rt = nbVecs - i - 1;

    if (rt > 0) {
      triFactor.row(i).tail(rt).noalias() = -hCoeffs(i) * vectors.col(i).tail(rs).adjoint() *
                                            vectors.bottomRightCorner(rs, rt).template triangularView<UnitLower>();

      triFactor.row(i).tail(rt) =
          (triFactor.row(i).tail(rt) * triFactor.bottomRightCorner(rt, rt).template triangularView<Upper>()).eval();
    }
    triFactor(i, i) = hCoeffs(i);
  }
}

/** \internal
 * if forward then perform   mat = H0 * H1 * H2 * mat
 * otherwise perform         mat = H2 * H1 * H0 * mat
 *
 * Implementation note: V (the householder vectors) is unit lower trapezoidal of shape
 * nbRows x nbVecs. Wrapping the *whole* V as TriangularView<UnitLower> would send both
 * V^* * mat and V * tmp through Eigen's triangular_matrix_matrix_product kernel, which
 * is single-threaded (it bypasses parallelize_gemm). We split V into its UnitLower top
 * nbVecs x nbVecs block and the general (nbRows - nbVecs) x nbVecs bottom block; the
 * bottom block — which is the bulk of the work for tall panels — then flows through
 * general_matrix_matrix_product, which parallelizes under OpenMP / EIGEN_GEMM_THREADPOOL.
 */
template <typename MatrixType, typename VectorsType, typename CoeffsType>
void apply_block_householder_on_the_left(MatrixType& mat, const VectorsType& vectors, const CoeffsType& hCoeffs,
                                         bool forward) {
  enum { TFactorSize = VectorsType::ColsAtCompileTime };
  const Index nbVecs = vectors.cols();
  const Index nbBelow = vectors.rows() - nbVecs;
  Matrix<typename MatrixType::Scalar, TFactorSize, TFactorSize, RowMajor> T(nbVecs, nbVecs);

  if (forward)
    make_block_householder_triangular_factor(T, vectors, hCoeffs);
  else
    make_block_householder_triangular_factor(T, vectors, hCoeffs.conjugate());

  const auto V_top = vectors.topRows(nbVecs);

  // tmp = V^* * mat, computed as V_top^* * mat.topRows(nbVecs) + V_bot^* * mat.bottomRows(nbBelow).
  Matrix<typename MatrixType::Scalar, VectorsType::ColsAtCompileTime, MatrixType::ColsAtCompileTime,
         (VectorsType::MaxColsAtCompileTime == 1 && MatrixType::MaxColsAtCompileTime != 1) ? RowMajor : ColMajor,
         VectorsType::MaxColsAtCompileTime, MatrixType::MaxColsAtCompileTime>
      tmp(nbVecs, mat.cols());
  tmp.noalias() = V_top.template triangularView<UnitLower>().adjoint() * mat.topRows(nbVecs);
  if (nbBelow > 0) {
    tmp.noalias() += vectors.bottomRows(nbBelow).adjoint() * mat.bottomRows(nbBelow);
  }

  if (forward)
    tmp = (T.template triangularView<Upper>() * tmp).eval();
  else
    tmp = (T.template triangularView<Upper>().adjoint() * tmp).eval();

  // mat -= V * tmp, split along the same top/bottom partition.
  mat.topRows(nbVecs).noalias() -= V_top.template triangularView<UnitLower>() * tmp;
  if (nbBelow > 0) {
    mat.bottomRows(nbBelow).noalias() -= vectors.bottomRows(nbBelow) * tmp;
  }
}

/** \internal
 * if forward then perform   mat = mat * H0 * H1 * H2
 * otherwise perform         mat = mat * H2 * H1 * H0
 */
template <typename MatrixType, typename VectorsType, typename CoeffsType>
void apply_block_householder_on_the_right(MatrixType& mat, const VectorsType& vectors, const CoeffsType& hCoeffs,
                                          bool forward) {
  enum { TFactorSize = VectorsType::ColsAtCompileTime };
  const Index nbVecs = vectors.cols();
  const Index nbBelow = vectors.rows() - nbVecs;
  Matrix<typename MatrixType::Scalar, TFactorSize, TFactorSize, RowMajor> T(nbVecs, nbVecs);

  if (forward)
    make_block_householder_triangular_factor(T, vectors, hCoeffs);
  else
    make_block_householder_triangular_factor(T, vectors, hCoeffs.conjugate());

  const auto V_top = vectors.topRows(nbVecs);

  // tmp = mat * V, split along V's top/bottom partition (see the left-apply for context).
  Matrix<typename MatrixType::Scalar, MatrixType::RowsAtCompileTime, VectorsType::ColsAtCompileTime,
         (MatrixType::MaxRowsAtCompileTime == 1 && VectorsType::MaxColsAtCompileTime != 1) ? ColMajor : RowMajor,
         MatrixType::MaxRowsAtCompileTime, VectorsType::MaxColsAtCompileTime>
      tmp(mat.rows(), nbVecs);
  tmp.noalias() = mat.leftCols(nbVecs) * V_top.template triangularView<UnitLower>();
  if (nbBelow > 0) {
    tmp.noalias() += mat.rightCols(nbBelow) * vectors.bottomRows(nbBelow);
  }

  if (forward)
    tmp = (tmp * T.template triangularView<Upper>()).eval();
  else
    tmp = (tmp * T.template triangularView<Upper>().adjoint()).eval();

  // mat -= tmp * V^*, split along the same partition.
  mat.leftCols(nbVecs).noalias() -= tmp * V_top.template triangularView<UnitLower>().adjoint();
  if (nbBelow > 0) {
    mat.rightCols(nbBelow).noalias() -= tmp * vectors.bottomRows(nbBelow).adjoint();
  }
}

}  // end namespace internal

}  // end namespace Eigen

#endif  // EIGEN_BLOCK_HOUSEHOLDER_H
