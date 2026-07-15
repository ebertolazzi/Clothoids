// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-FileCopyrightText: The Eigen Authors
// SPDX-License-Identifier: MPL-2.0

#ifndef EIGEN_BLOCKSPARSEMATRIX_H
#define EIGEN_BLOCKSPARSEMATRIX_H

// IWYU pragma: private
#include "./InternalHeaderCheck.h"

#include <algorithm>
#include <numeric>
#include <utility>

namespace Eigen {

// Forward declarations
template <typename, int, bool>
class BlockSparseTriangularView;
template <typename, int, bool>
class BlockSparseSelfAdjointView;
template <typename Scalar_, int Options_, int BlockRows_, int BlockCols_, typename StorageIndex_>
class BlockSparseMatrix;

/** Storage-kind tag for BlockSparseMatrix. */
struct BlockSparse {};

/** Evaluator shape tag for BlockSparseMatrix product dispatch. */
struct BlockSparseShape {
  static std::string debugName() { return "BlockSparseShape"; }
};

namespace internal {
// Returns m.adjoint() when Conj==true, m.transpose() otherwise.
// SFINAE overloads keep the return type concrete under C++14 (no if constexpr).
template <bool Conj, typename T>
std::enable_if_t<Conj, decltype(std::declval<const T&>().adjoint())> adjoint_if(const T& m) {
  return m.adjoint();
}

template <bool Conj, typename T>
std::enable_if_t<!Conj, decltype(std::declval<const T&>().transpose())> adjoint_if(const T& m) {
  return m.transpose();
}
template <>
struct storage_kind_to_evaluator_kind<BlockSparse> {
  using Kind = IndexBased;
};

template <>
struct storage_kind_to_shape<BlockSparse> {
  using Shape = BlockSparseShape;
};

template <typename Scalar_, int Options_, int BlockRows_, int BlockCols_, typename StorageIndex_>
struct traits<BlockSparseMatrix<Scalar_, Options_, BlockRows_, BlockCols_, StorageIndex_>> {
  using Scalar = Scalar_;
  using StorageIndex = StorageIndex_;
  using StorageKind = BlockSparse;
  using XprKind = MatrixXpr;

  static constexpr Index RowsAtCompileTime = Dynamic;
  static constexpr Index ColsAtCompileTime = Dynamic;
  static constexpr Index MaxRowsAtCompileTime = Dynamic;
  static constexpr Index MaxColsAtCompileTime = Dynamic;
  static constexpr int Options = Options_;
  static constexpr unsigned int Flags = Options_ | NestByRefBit | LvalueBit;
};

}  // namespace internal

/** \class BlockTriplet
 * \ingroup SparseCore_Module
 * \brief A (blockRow, blockCol, blockValue) triplet for assembling a BlockSparseMatrix.
 *
 * Coordinates are in \em block space, not element space.
 *
 * \tparam Scalar_       Numeric scalar type.
 * \tparam BlockRows_    Number of rows in each block.
 * \tparam BlockCols_    Number of columns in each block.
 * \tparam StorageIndex_ Signed integer index type (default: int).
 */
template <typename Scalar_, int BlockRows_, int BlockCols_, int Options_ = ColMajor, typename StorageIndex_ = int>
class BlockTriplet {
 public:
  using Scalar = Scalar_;
  using StorageIndex = StorageIndex_;
  using BlockType = Matrix<Scalar, BlockRows_, BlockCols_, Options_>;
  using BlockMapType = Map<BlockType, Unaligned>;
  using ConstBlockMapType = Map<const BlockType, Unaligned>;

  static constexpr int BlockSize = BlockRows_ * BlockCols_;

  BlockTriplet() = default;

  BlockTriplet(StorageIndex blockRow, StorageIndex blockCol, const BlockType& block)
      : m_row(blockRow), m_col(blockCol) {
    BlockMapType{m_value} = block;
  }

  StorageIndex row() const { return m_row; }
  StorageIndex col() const { return m_col; }
  // Implicitly usable wherever a MatrixBase expression is expected.
  ConstBlockMapType value() const { return ConstBlockMapType(m_value); }

 private:
  StorageIndex m_row = 0;
  StorageIndex m_col = 0;
  // Flat array avoids the alignment padding that a Matrix<> member would incur.
  Scalar m_value[BlockSize];
};

/** \class BlockSparseMatrix
 * \ingroup SparseCore_Module
 * \brief A sparse matrix whose stored nonzeros are fixed-size dense blocks.
 *
 * Each nonzero entry is a \c BlockRows x \c BlockCols dense matrix.  The
 * block sparsity pattern is stored in block-level compressed-column
 * (ColMajor) or compressed-row (RowMajor) format.
 *
 * \tparam Scalar_       Numeric scalar type.
 * \tparam Options_      ColMajor (0) or RowMajor.  Controls both the outer
 *                       iteration direction over blocks and the storage
 *                       layout within each block.  Note: vector-shaped blocks
 *                       must use a compatible order — \c BlockCols_=1 requires
 *                       ColMajor, \c BlockRows_=1 requires RowMajor.
 * \tparam BlockRows_    Rows per block; must be a fixed positive integer.
 * \tparam BlockCols_    Columns per block; must be a fixed positive integer.
 * \tparam StorageIndex_ Signed integer type for internal index arrays
 *                       (default: int).
 *
 * ### Assembly
 * Populate the matrix via setFromTriplets(), passing an iterator range of
 * BlockTriplet objects in block coordinates.  Triplets with the same
 * (blockRow, blockCol) are summed.
 *
 * ### Arithmetic
 * Addition, subtraction, and matrix product are supported between compatible
 * BlockSparseMatrix instances.  Scalar multiplication is also available.
 *
 * The matrix product \c C = A * B requires
 * \c A.BlockCols == B.BlockRows (enforced at compile time by the template
 * constraint) and \c A.blockCols() == B.blockRows() (checked at runtime).
 * The result has block type \c Matrix<Scalar,A.BlockRows,B.BlockCols>.
 *
 * ### Conversion
 * toSparse() converts to a standard SparseMatrix with element-level
 * sparsity.  fromSparse() reconstructs a BlockSparseMatrix from an
 * element-level SparseMatrix whose dimensions are divisible by BlockRows
 * and BlockCols.  An implicit conversion operator to SparseMatrix is
 * also provided.
 */
template <typename Scalar_, int Options_, int BlockRows_, int BlockCols_, typename StorageIndex_ = int>
class BlockSparseMatrix
    : public EigenBase<BlockSparseMatrix<Scalar_, Options_, BlockRows_, BlockCols_, StorageIndex_>> {
  EIGEN_STATIC_ASSERT(BlockRows_ >= 1, BLOCKROWS_MUST_BE_A_POSITIVE_COMPILE_TIME_SIZE)
  EIGEN_STATIC_ASSERT(BlockCols_ >= 1, BLOCKCOLS_MUST_BE_A_POSITIVE_COMPILE_TIME_SIZE)
  EIGEN_STATIC_ASSERT(std::is_integral<StorageIndex_>::value&& std::is_signed<StorageIndex_>::value,
                      STORAGEINDEX_MUST_BE_A_SIGNED_INTEGRAL_TYPE)
  // Eigen's Matrix<> requires: a column vector (Cols==1, Rows>1) must be ColMajor;
  // a row vector (Rows==1, Cols>1) must be RowMajor.  Guard those cases here so
  // the error fires at BlockSparseMatrix instantiation rather than inside BlockType.
  EIGEN_STATIC_ASSERT(!(BlockCols_ == 1 && BlockRows_ != 1 && bool(Options_ & RowMajorBit)),
                      INVALID_MATRIX_TEMPLATE_PARAMETERS)
  EIGEN_STATIC_ASSERT(!(BlockRows_ == 1 && BlockCols_ != 1 && !bool(Options_ & RowMajorBit)),
                      INVALID_MATRIX_TEMPLATE_PARAMETERS)

 public:
  // -------------------------------------------------------------------------
  // Type aliases & compile-time constants
  // -------------------------------------------------------------------------
  using Scalar = Scalar_;
  using StorageIndex = StorageIndex_;
  using BlockType = Matrix<Scalar, BlockRows_, BlockCols_, Options_>;
  using TripletType = BlockTriplet<Scalar, BlockRows_, BlockCols_, Options_, StorageIndex>;

  static constexpr int Options = Options_;
  static constexpr Index BlockRows = BlockRows_;
  static constexpr Index BlockCols = BlockCols_;
  static constexpr bool IsRowMajor = Options_ & RowMajorBit;
  static constexpr Index BlockSize = BlockRows_ * BlockCols_;

  // If one block occupies a power-of-two number of bytes, and the values array
  // is Eigen-allocated (guaranteed aligned to EIGEN_MAX_ALIGN_BYTES), then every
  // block pointer is aligned to min(BlockBytes, EIGEN_MAX_ALIGN_BYTES).
  static constexpr std::size_t BlockBytes = std::size_t(BlockSize) * sizeof(Scalar);
  static constexpr int BlockMapAlignment = ((BlockBytes & (BlockBytes - 1)) == 0 && BlockBytes >= 8)
                                               ? int(numext::mini(BlockBytes, std::size_t(EIGEN_MAX_ALIGN_BYTES)))
                                               : 0;

  using BlockMap = Map<BlockType, BlockMapAlignment>;
  using ConstBlockMap = Map<const BlockType, BlockMapAlignment>;

  // -------------------------------------------------------------------------
  // Constructors / copy / move
  // -------------------------------------------------------------------------

  /** Default constructor; creates a 0×0 matrix. */
  BlockSparseMatrix() = default;

  /** Construct a zero matrix with the given number of block-rows and block-columns. */
  BlockSparseMatrix(Index blockRows, Index blockCols)
      : m_blockOuterSize(IsRowMajor ? blockRows : blockCols), m_blockInnerSize(IsRowMajor ? blockCols : blockRows) {}

  // -------------------------------------------------------------------------
  // Dimensions
  // -------------------------------------------------------------------------

  /** Total number of element rows. */
  Index rows() const noexcept { return (IsRowMajor ? m_blockOuterSize : m_blockInnerSize) * BlockRows_; }
  /** Total number of element columns. */
  Index cols() const noexcept { return (IsRowMajor ? m_blockInnerSize : m_blockOuterSize) * BlockCols_; }

  /** Number of block-rows. */
  Index blockRows() const { return IsRowMajor ? m_blockOuterSize : m_blockInnerSize; }
  /** Number of block-columns. */
  Index blockCols() const { return IsRowMajor ? m_blockInnerSize : m_blockOuterSize; }

  /** Outer block dimension (block-cols for ColMajor, block-rows for RowMajor). */
  Index blockOuterSize() const { return m_blockOuterSize; }
  /** Inner block dimension (block-rows for ColMajor, block-cols for RowMajor). */
  Index blockInnerSize() const { return m_blockInnerSize; }

  /** Element-level outer size: cols() for ColMajor, rows() for RowMajor. */
  Index outerSize() const { return IsRowMajor ? rows() : cols(); }
  /** Element-level inner size: rows() for ColMajor, cols() for RowMajor. */
  Index innerSize() const { return IsRowMajor ? cols() : rows(); }

  /** Number of stored (structurally non-zero) blocks. */
  Index nonZeroBlocks() const { return m_outerIndex(m_blockOuterSize); }
  /** Total number of stored scalar coefficients (= nonZeroBlocks() * BlockRows * BlockCols). */
  Index nonZeros() const { return nonZeroBlocks() * BlockSize; }
  /** Number of blocks for which storage is currently allocated (capacity). */
  Index allocatedBlocks() const { return Index(m_innerIndex.size()); }

  // -------------------------------------------------------------------------
  // Raw pointer access (for interoperability)
  // -------------------------------------------------------------------------
  const StorageIndex* outerIndexPtr() const { return m_outerIndex.data(); }
  StorageIndex* outerIndexPtr() { return m_outerIndex.data(); }
  const StorageIndex* innerIndexPtr() const { return m_innerIndex.data(); }
  StorageIndex* innerIndexPtr() { return m_innerIndex.data(); }
  const Scalar* valuePtr() const { return m_values.data(); }
  Scalar* valuePtr() { return m_values.data(); }

  // -------------------------------------------------------------------------
  // Block access by sequential nonzero index
  // -------------------------------------------------------------------------

  /** Read-only Map to the \a k-th stored block (block storage follows \c Options_). */
  ConstBlockMap blockRef(Index k) const { return ConstBlockMap(m_values.data() + k * BlockSize); }
  /** Mutable Map to the \a k-th stored block. */
  BlockMap blockRef(Index k) { return BlockMap(m_values.data() + k * BlockSize); }

  // -------------------------------------------------------------------------
  // Inner iterator over blocks within one outer vector
  // -------------------------------------------------------------------------

  /** \brief Iterates over stored blocks in outer vector \a outer.
   *
   * Usage mirrors SparseMatrix::InnerIterator but value() returns a
   * Map to a BlockRows×BlockCols matrix, not a scalar.
   */
  class InnerIterator {
   public:
    EIGEN_STRONG_INLINE InnerIterator(const BlockSparseMatrix& mat, Index outer)
        : m_mat(mat), m_id(mat.m_outerIndex(outer)), m_end(mat.m_outerIndex(outer + 1)), m_outer(outer) {}

    EIGEN_STRONG_INLINE operator bool() const { return m_id < m_end; }
    EIGEN_STRONG_INLINE InnerIterator& operator++() {
      ++m_id;
      return *this;
    }

    /** Current block outer index (block-col for ColMajor, block-row for RowMajor). */
    EIGEN_STRONG_INLINE Index outer() const { return m_outer; }
    /** Current block inner index (block-row for ColMajor, block-col for RowMajor). */
    EIGEN_STRONG_INLINE Index index() const { return m_mat.m_innerIndex(m_id); }
    /** Block-row of the current block. */
    EIGEN_STRONG_INLINE Index blockRow() const { return IsRowMajor ? m_outer : index(); }
    /** Block-column of the current block. */
    EIGEN_STRONG_INLINE Index blockCol() const { return IsRowMajor ? index() : m_outer; }

    /** Read-only Map to the current block value. */
    EIGEN_STRONG_INLINE ConstBlockMap value() const { return m_mat.blockRef(m_id); }
    /** Mutable Map to the current block value. */
    EIGEN_STRONG_INLINE BlockMap valueRef() {
      return BlockMap(const_cast<Scalar*>(m_mat.m_values.data()) + m_id * BlockSize);
    }

   private:
    const BlockSparseMatrix& m_mat;
    Index m_id;
    Index m_end;
    Index m_outer;
  };

  // -------------------------------------------------------------------------
  // Resize / clear
  // -------------------------------------------------------------------------

  /** Resize to \a blockRows × \a blockCols blocks and set the logical nnz to zero.
   *  Allocated block storage is retained; call squeeze() to release it. */
  void resize(Index blockRows, Index blockCols) {
    m_blockOuterSize = IsRowMajor ? blockRows : blockCols;
    m_blockInnerSize = IsRowMajor ? blockCols : blockRows;
    m_outerIndex.resize(m_blockOuterSize + 1);
    m_outerIndex.setZero();
  }

  /** Clear all stored blocks (logical nnz → 0) while keeping dimensions and allocated storage. */
  void setZero() {
    m_outerIndex.resize(m_blockOuterSize + 1);
    m_outerIndex.setZero();
  }

  /** Pre-allocate storage for at least \a n blocks without changing the logical sparsity pattern.
   *  Existing block data is preserved up to min(n, nonZeroBlocks()). */
  void reserve(Index n) {
    if (n > Index(m_innerIndex.size())) conservativeResizeBlockStorage_(n);
  }

  /** Release any excess allocated block storage so that allocatedBlocks() == nonZeroBlocks(). */
  void squeeze() {
    Index nnz = nonZeroBlocks();
    if (nnz < Index(m_innerIndex.size())) conservativeResizeBlockStorage_(nnz);
  }

  /** Fill the matrix with the block identity: the min(blockRows,blockCols) diagonal blocks
   *  are set to the B×B identity; all other blocks are absent.
   *
   *  \pre BlockRows == BlockCols (square blocks).
   */
  void setIdentity() {
    EIGEN_STATIC_ASSERT(BlockRows_ == BlockCols_, THIS_METHOD_IS_ONLY_FOR_SQUARE_BLOCK_MATRICES)
    Index n = (std::min)(m_blockOuterSize, m_blockInnerSize);
    m_outerIndex.resize(m_blockOuterSize + 1);
    resizeBlockStorage_(n);
    for (Index i = 0; i <= m_blockOuterSize; ++i) m_outerIndex(i) = StorageIndex((std::min)(i, n));
    for (StorageIndex i = 0; i < n; ++i) {
      m_innerIndex(i) = i;
      blockRef(i).setIdentity();
    }
  }

  /** Initialize the block structure directly from compressed outer/inner index arrays,
   *  zero-initializing all block values.
   *
   *  \p outerPtr  has size blockCols+1 (ColMajor) or blockRows+1 (RowMajor).
   *  \p innerPtr  has size nnzBlocks.
   */
  void setFromOuterInner(Index blockRows, Index blockCols, Index nnzBlocks, const StorageIndex_* outerPtr,
                         const StorageIndex_* innerPtr) {
    m_blockOuterSize = IsRowMajor ? blockRows : blockCols;
    m_blockInnerSize = IsRowMajor ? blockCols : blockRows;
    m_outerIndex = Map<const decltype(m_outerIndex)>(outerPtr, m_blockOuterSize + 1);
    m_innerIndex = Map<const decltype(m_innerIndex)>(innerPtr, nnzBlocks);
    resizeBlockStorage_(nnzBlocks);
    m_values.setZero();
  }

  // -------------------------------------------------------------------------
  // Assembly
  // -------------------------------------------------------------------------

  /** Fill the matrix from an iterator range of BlockTriplet objects.
   *
   * Triplet coordinates are in block space.  Triplets with the same
   * (blockRow, blockCol) pair are summed (their block values are added).
   * The input range may be in any order.
   *
   * \tparam InputIterator  Must dereference to a type with \c row(),
   *                        \c col(), and \c value() members, matching
   *                        BlockTriplet's interface.
   */
  template <typename InputIterator>
  void setFromTriplets(InputIterator begin, InputIterator end);

  // -------------------------------------------------------------------------
  // Conversion to / from SparseMatrix
  // -------------------------------------------------------------------------

  /** Convert to a SparseMatrix with scalar nonzeros.
   *
   * Each stored block of size BlockRows×BlockCols expands into up to
   * BlockRows*BlockCols scalar nonzeros.  The resulting SparseMatrix has
   * the same storage order as \c *this.
   */
  SparseMatrix<Scalar, Options_, StorageIndex_> toSparse() const;

  /** Construct a BlockSparseMatrix from an element-level SparseMatrix.
   *
   * \pre  \c sp.rows() % BlockRows == 0 and \c sp.cols() % BlockCols == 0.
   *
   * Each scalar entry \c sp(i,j) is placed into position
   * \c (i%BlockRows, j%BlockCols) of block \c (i/BlockRows, j/BlockCols).
   * Multiple entries mapping to the same element of the same block are
   * accumulated with \c +=.
   */
  static BlockSparseMatrix fromSparse(const SparseMatrix<Scalar_, Options_, StorageIndex_>& sp);

  /** Implicit conversion to SparseMatrix. */
  operator SparseMatrix<Scalar_, Options_, StorageIndex_>() const { return toSparse(); }

  // -------------------------------------------------------------------------
  // Element access
  // -------------------------------------------------------------------------

  /** Read element \c (row, col); returns 0 if no block covers that position. */
  Scalar coeff(Index row, Index col) const {
    eigen_assert(row >= 0 && row < rows() && col >= 0 && col < cols());
    Index bOuter = IsRowMajor ? (row / BlockRows_) : (col / BlockCols_);
    Index bInner = IsRowMajor ? (col / BlockCols_) : (row / BlockRows_);
    Index localRow = row % BlockRows_;
    Index localCol = col % BlockCols_;
    const StorageIndex* beg = m_innerIndex.data() + m_outerIndex(bOuter);
    const StorageIndex* fin = m_innerIndex.data() + m_outerIndex(bOuter + 1);
    const StorageIndex* it = std::lower_bound(beg, fin, StorageIndex(bInner));
    if (it == fin || *it != bInner) return Scalar(0);
    return blockRef(static_cast<Index>(it - m_innerIndex.data()))(localRow, localCol);
  }

  /** Extract the main scalar diagonal as a dense vector.
   *
   * Iterates outer slices once and binary-searches for the diagonal block in each
   * slice, then copies the relevant entries from that block — one search per outer
   * slice (square blocks) or per unique inner-block boundary (non-square blocks),
   * vs. one search per scalar element for the coeff-by-coeff approach.
   */
  Matrix<Scalar, Dynamic, 1> diagonal() const {
    constexpr Index OuterB = IsRowMajor ? BlockRows_ : BlockCols_;
    constexpr Index InnerB = IsRowMajor ? BlockCols_ : BlockRows_;
    const Index diagSize = numext::mini(rows(), cols());
    Matrix<Scalar, Dynamic, 1> diag = Matrix<Scalar, Dynamic, 1>::Zero(diagSize);

    for (Index out = 0; out < m_blockOuterSize; ++out) {
      const Index scalarOuterBegin = out * OuterB;
      if (scalarOuterBegin >= diagSize) break;
      const Index scalarOuterEnd = numext::mini(scalarOuterBegin + OuterB, diagSize);

      // Group consecutive scalar positions that share the same inner block, then
      // binary-search once per group rather than once per scalar element.
      // For square blocks this loop runs exactly once per outer slice.
      Index i = scalarOuterBegin;
      while (i < scalarOuterEnd) {
        const Index bInner = i / InnerB;
        const Index groupEnd = numext::mini((bInner + 1) * InnerB, scalarOuterEnd);

        const StorageIndex* beg = m_innerIndex.data() + m_outerIndex(out);
        const StorageIndex* fin = m_innerIndex.data() + m_outerIndex(out + 1);
        const StorageIndex* it = std::lower_bound(beg, fin, StorageIndex(bInner));
        if (it != fin && *it == StorageIndex(bInner)) {
          const ConstBlockMap blk = blockRef(static_cast<Index>(it - m_innerIndex.data()));
          for (Index j = i; j < groupEnd; ++j) {
            const Index localRow = IsRowMajor ? (j % OuterB) : (j % InnerB);
            const Index localCol = IsRowMajor ? (j % InnerB) : (j % OuterB);
            diag(j) = blk(localRow, localCol);
          }
        }
        i = groupEnd;
      }
    }
    return diag;
  }

  // -------------------------------------------------------------------------
  // Arithmetic
  // -------------------------------------------------------------------------

  /** Element-wise addition.  Both matrices must have the same block dimensions. */
  BlockSparseMatrix operator+(const BlockSparseMatrix& other) const { return disjunctionWith_(other, AddOp_{}); }

  /** Element-wise subtraction. */
  BlockSparseMatrix operator-(const BlockSparseMatrix& other) const { return disjunctionWith_(other, SubOp_{}); }

  /** Element-wise product.  Only blocks present in \em both operands contribute to the result. */
  BlockSparseMatrix cwiseProduct(const BlockSparseMatrix& other) const {
    return conjunctionWith_(other, CwiseMulOp_{});
  }

  /** Applies a scalar unary functor to every stored nonzero, preserving the sparsity pattern. */
  template <typename ScalarFunc>
  BlockSparseMatrix unaryExpr(ScalarFunc func) const {
    return withValues_([&func](const auto& v) { return v.unaryExpr(func); });
  }

  /** Applies a scalar binary functor with union sparsity.
   *
   * \p func must provide three members:
   * \code
   *   Scalar func(Scalar a, Scalar b)  // both present
   *   Scalar func.lhs(Scalar a)        // only lhs present; rhs is implicitly zero
   *   Scalar func.rhs(Scalar b)        // only rhs present; lhs is implicitly zero
   * \endcode
   */
  template <typename ScalarFunc>
  BlockSparseMatrix disjunctionExpr(const BlockSparseMatrix& other, ScalarFunc func) const {
    return disjunctionWith_(other, DisjExprAdapter_<ScalarFunc>{func});
  }

  /** Applies a scalar binary functor with intersection sparsity: only block positions present
   * in \em both matrices contribute; \p func is called as \c func(Scalar a, Scalar b). */
  template <typename ScalarFunc>
  BlockSparseMatrix conjunctionExpr(const BlockSparseMatrix& other, ScalarFunc func) const {
    return conjunctionWith_(other, [&func](const auto& a, const auto& b) { return a.binaryExpr(b, func); });
  }

  /** Unary negation. */
  BlockSparseMatrix operator-() const {
    return withValues_([](const auto& v) { return -v; });
  }

  BlockSparseMatrix& operator+=(const BlockSparseMatrix& other) { return *this = *this + other; }
  BlockSparseMatrix& operator-=(const BlockSparseMatrix& other) { return *this = *this - other; }

  /** Scalar multiplication (returns a new matrix). */
  BlockSparseMatrix operator*(const Scalar& s) const {
    return withValues_([&s](const auto& v) { return v * s; });
  }
  BlockSparseMatrix& operator*=(const Scalar& s) {
    m_values.head(nonZeros()) *= s;
    return *this;
  }
  BlockSparseMatrix operator/(const Scalar& s) const {
    return withValues_([&s](const auto& v) { return v / s; });
  }
  BlockSparseMatrix& operator/=(const Scalar& s) { return *this *= (Scalar(1) / s); }

  /** Scalar-on-left multiplication. */
  friend BlockSparseMatrix operator*(const Scalar& s, const BlockSparseMatrix& m) { return m * s; }

  /** Block-sparse times dense matrix (or vector) product.
   *
   * Returns a lazy \c Product<> expression evaluated via \c generic_product_impl.
   * This enables fused accumulation:
   *   \code
   *   b.noalias()  = A * x;       // no temporary
   *   b.noalias() += A * x;       // fused add
   *   b.noalias() += alpha * (A * x);  // scale then add via evaluator
   *   \endcode
   *
   * \pre  \c this->cols() == rhs.rows().
   * \pre  Scalar types must match.
   * \warning The result uses \c AliasFreeProduct, so assignment goes directly
   *          through \c generic_product_impl::evalTo with no aliasing temporary.
   *          \c x = A * x silently corrupts; use an explicit temporary if needed.
   */
  template <typename OtherDerived>
  Product<BlockSparseMatrix, OtherDerived, AliasFreeProduct> operator*(const MatrixBase<OtherDerived>& rhs) const {
    EIGEN_STATIC_ASSERT(
        (std::is_same<Scalar, typename OtherDerived::Scalar>::value),
        YOU_MIXED_DIFFERENT_NUMERIC_TYPES__YOU_NEED_TO_USE_THE_CAST_METHOD_OF_MATRIXBASE_TO_CAST_NUMERIC_TYPES_EXPLICITLY)
    return Product<BlockSparseMatrix, OtherDerived, AliasFreeProduct>(*this, rhs.derived());
  }

  /** Dense matrix (or vector) times block-sparse product (hidden friend).
   *
   * Returns a lazy \c Product<> expression; evaluated via \c generic_product_impl.
   *
   * \pre  \c lhs.cols() == bsm.rows().
   * \pre  Scalar types must match.
   * \warning The result uses \c AliasFreeProduct, so assignment goes directly
   *          through \c generic_product_impl::evalTo with no aliasing temporary.
   *          \c x = x * A silently corrupts; use an explicit temporary if needed.
   */
  template <typename OtherDerived>
  friend Product<OtherDerived, BlockSparseMatrix, AliasFreeProduct> operator*(const MatrixBase<OtherDerived>& lhs,
                                                                              const BlockSparseMatrix& bsm) {
    EIGEN_STATIC_ASSERT(
        (std::is_same<Scalar_, typename OtherDerived::Scalar>::value),
        YOU_MIXED_DIFFERENT_NUMERIC_TYPES__YOU_NEED_TO_USE_THE_CAST_METHOD_OF_MATRIXBASE_TO_CAST_NUMERIC_TYPES_EXPLICITLY)
    return Product<OtherDerived, BlockSparseMatrix, AliasFreeProduct>(lhs.derived(), bsm);
  }

  /** Block-sparse matrix product.
   *
   * \tparam RhsBlockCols  Block-column count of the right-hand side.
   *
   * The left-hand side has block type \c Matrix<Scalar,BlockRows,BlockCols>
   * and the right-hand side must have block type
   * \c Matrix<Scalar,BlockCols,RhsBlockCols> (enforced by the template
   * parameter).  The result has block type
   * \c Matrix<Scalar,BlockRows,RhsBlockCols>.
   *
   * \pre  \c this->blockCols() == rhs.blockRows() (checked at runtime).
   */
  template <int RhsBlockCols>
  BlockSparseMatrix<Scalar_, Options_, BlockRows_, RhsBlockCols, StorageIndex_> operator*(
      const BlockSparseMatrix<Scalar_, Options_, BlockCols_, RhsBlockCols, StorageIndex_>& rhs) const;

  // -------------------------------------------------------------------------
  // Approximate equality (useful for testing)
  // -------------------------------------------------------------------------
  bool isApprox(const BlockSparseMatrix& other,
                const typename NumTraits<Scalar>::Real& prec = NumTraits<Scalar>::dummy_precision()) const {
    using RealScalar = typename NumTraits<Scalar>::Real;
    // Frobenius-norm comparison, matching SparseMatrixBase::isApprox semantics
    // but computed directly from the block values — no scalar-level SparseMatrix
    // materialization of either operand.  Explicit zero blocks contribute 0 to
    // every norm, so the result is independent of structural differences.
    RealScalar n2a = m_values.head(nonZeros()).matrix().squaredNorm();
    RealScalar n2b = other.m_values.head(other.nonZeros()).matrix().squaredNorm();
    BlockSparseMatrix diff = *this - other;
    RealScalar d2 = diff.m_values.head(diff.nonZeros()).matrix().squaredNorm();
    return d2 <= prec * prec * numext::mini(n2a, n2b);
  }

  // -------------------------------------------------------------------------
  // Transpose / adjoint
  // -------------------------------------------------------------------------

  /** Returns a new BSM whose block dimensions are swapped (BlockRows ↔ BlockCols)
   *  and each stored block is transposed.  Storage order is preserved. */
  BlockSparseMatrix<Scalar_, Options_, BlockCols_, BlockRows_, StorageIndex_> transpose() const;

  /** Returns a new BSM whose block dimensions are swapped and each stored block
   *  is conjugate-transposed (adjoint).  Identical to transpose() for real scalars. */
  BlockSparseMatrix<Scalar_, Options_, BlockCols_, BlockRows_, StorageIndex_> adjoint() const;

  // -------------------------------------------------------------------------
  // View factories
  // -------------------------------------------------------------------------

  /** Returns a block-level triangular view.
   *
   * \tparam Mode             \c Eigen::Upper or \c Eigen::Lower.
   * \tparam DiagIsTriangular When \c false (default), diagonal blocks are
   *   treated as triangular regardless of what is stored in the unused
   *   triangle — products use \c triangularView on the block and \c eval()
   *   explicitly zeros the unused triangle.  When \c true, the caller
   *   guarantees that the unused triangle of every diagonal block is already
   *   zero; products then use a full vectorised GEMV (faster, no zeroing).
   */
  template <int Mode, bool DiagIsTriangular = false>
  BlockSparseTriangularView<BlockSparseMatrix, Mode, DiagIsTriangular> triangularView() const {
    return BlockSparseTriangularView<BlockSparseMatrix, Mode, DiagIsTriangular>(*this);
  }

  /** Returns a block-level self-adjoint view.
   *
   * Only the triangle selected by \p UpLo is read; the opposite triangle is
   * reconstructed on-the-fly as the adjoint of each stored off-diagonal block.
   *
   * \tparam UpLo              \c Eigen::Upper or \c Eigen::Lower.
   * \tparam DiagIsSelfAdjoint Set to \c true when every diagonal block is
   *         itself Hermitian \em and both triangles are explicitly stored.
   *         The dense sub-product then uses a plain product rather than
   *         \c selfadjointView, which is more efficient for the small
   *         fixed-size blocks typical here.  When \c false (the default),
   *         only the \p UpLo triangle of each diagonal block is assumed
   *         valid; \c selfadjointView<UpLo>() reconstructs the full block.
   *
   * \pre BlockRows == BlockCols (diagonal blocks must be square).
   */
  template <int UpLo, bool DiagIsSelfAdjoint = false>
  BlockSparseSelfAdjointView<BlockSparseMatrix, UpLo, DiagIsSelfAdjoint> selfadjointView() const {
    EIGEN_STATIC_ASSERT(BlockRows_ == BlockCols_, THIS_METHOD_IS_ONLY_FOR_SQUARE_BLOCK_MATRICES)
    return BlockSparseSelfAdjointView<BlockSparseMatrix, UpLo, DiagIsSelfAdjoint>(*this);
  }

 private:
  struct AddOp_ {
    template <typename A, typename B>
    BlockType operator()(const A& a, const B& b) const {
      return a + b;
    }
    template <typename A>
    BlockType lhs(const A& a) const {
      return a;
    }
    template <typename B>
    BlockType rhs(const B& b) const {
      return b;
    }
  };

  struct SubOp_ {
    template <typename A, typename B>
    BlockType operator()(const A& a, const B& b) const {
      return a - b;
    }
    template <typename A>
    BlockType lhs(const A& a) const {
      return a;
    }
    template <typename B>
    BlockType rhs(const B& b) const {
      return -b;
    }
  };

  // Scalar-to-block adapter for disjunctionExpr: translates scalar functor with lhs/rhs methods
  // to the block level.
  template <typename ScalarFunc>
  struct DisjExprAdapter_ {
    ScalarFunc func_;
    template <typename A, typename B>
    BlockType operator()(const A& a, const B& b) const {
      return a.binaryExpr(b, func_);
    }
    template <typename A>
    BlockType lhs(const A& a) const {
      return a.unaryExpr([this](const Scalar& x) { return func_.lhs(x); });
    }
    template <typename B>
    BlockType rhs(const B& b) const {
      return b.unaryExpr([this](const Scalar& x) { return func_.rhs(x); });
    }
  };

  struct CwiseMulOp_ {
    template <typename A, typename B>
    BlockType operator()(const A& a, const B& b) const {
      return a.cwiseProduct(b);
    }
  };

  // Returns a copy with the same sparsity structure but m_values replaced by f(m_values).
  // f receives the flat Eigen Array of the logical (non-zero) coefficients and returns any
  // compatible expression.  Only the structure is copied — the source values are never
  // duplicated, and unused tail capacity is neither copied nor evaluated.
  template <typename F>
  BlockSparseMatrix withValues_(F f) const {
    Index nnz = nonZeroBlocks();
    BlockSparseMatrix result(blockRows(), blockCols());
    result.m_outerIndex = m_outerIndex;
    result.m_innerIndex = m_innerIndex.head(nnz);
    result.m_values = f(m_values.head(nnz * BlockSize));
    return result;
  }

  // Disjunction (union-pattern): result has a block wherever *this OR other has one.
  //   lhs-only:  block copied from *this unchanged
  //   rhs-only:  op(b)      — unary overload of op
  //   both:      op(a, b)   — binary overload of op
  template <typename Op>
  BlockSparseMatrix disjunctionWith_(const BlockSparseMatrix& other, Op op) const {
    eigen_assert(blockRows() == other.blockRows() && blockCols() == other.blockCols() &&
                 "BlockSparseMatrix size mismatch");
    BlockSparseMatrix result(blockRows(), blockCols());
    result.resizeBlockStorage_(nonZeroBlocks() + other.nonZeroBlocks());
    Index nnz = 0;
    for (Index j = 0; j < m_blockOuterSize; ++j) {
      result.m_outerIndex(j) = StorageIndex_(nnz);
      Index aId = m_outerIndex(j);
      Index aEnd = m_outerIndex(j + 1);
      Index bId = other.m_outerIndex(j);
      Index bEnd = other.m_outerIndex(j + 1);
      while (aId < aEnd || bId < bEnd) {
        bool hasA = aId < aEnd;
        bool hasB = bId < bEnd;
        Index aInner = hasA ? Index(m_innerIndex(aId)) : -1;
        Index bInner = hasB ? Index(other.m_innerIndex(bId)) : -1;
        // Write the result block in place; the op overloads return a BlockType,
        // which is assigned straight into the result's storage Map.
        if (hasA && (!hasB || aInner < bInner)) {
          result.m_innerIndex(nnz) = StorageIndex_(aInner);
          result.blockRef(nnz) = op.lhs(blockRef(aId++));
        } else if (hasB && (!hasA || bInner < aInner)) {
          result.m_innerIndex(nnz) = StorageIndex_(bInner);
          result.blockRef(nnz) = op.rhs(other.blockRef(bId++));
        } else {
          result.m_innerIndex(nnz) = StorageIndex_(aInner);
          result.blockRef(nnz) = op(blockRef(aId++), other.blockRef(bId++));
        }
        ++nnz;
      }
    }
    result.m_outerIndex(m_blockOuterSize) = StorageIndex_(nnz);
    result.conservativeResizeBlockStorage_(nnz);
    return result;
  }

  // Conjunction (intersection-pattern): result has a block only where *this AND other both have one.
  template <typename BinaryOp>
  BlockSparseMatrix conjunctionWith_(const BlockSparseMatrix& other, BinaryOp func) const {
    eigen_assert(blockRows() == other.blockRows() && blockCols() == other.blockCols() &&
                 "BlockSparseMatrix size mismatch");
    BlockSparseMatrix result(blockRows(), blockCols());
    result.resizeBlockStorage_((std::min)(nonZeroBlocks(), other.nonZeroBlocks()));
    Index nnz = 0;
    for (Index j = 0; j < m_blockOuterSize; ++j) {
      result.m_outerIndex(j) = StorageIndex_(nnz);
      Index aId = m_outerIndex(j);
      Index aEnd = m_outerIndex(j + 1);
      Index bId = other.m_outerIndex(j);
      Index bEnd = other.m_outerIndex(j + 1);
      while (aId < aEnd && bId < bEnd) {
        Index aInner = m_innerIndex(aId);
        Index bInner = other.m_innerIndex(bId);
        if (aInner < bInner) {
          ++aId;
        } else if (bInner < aInner) {
          ++bId;
        } else {
          result.m_innerIndex(nnz) = StorageIndex_(aInner);
          result.blockRef(nnz) = func(blockRef(aId++), other.blockRef(bId++));
          ++nnz;
        }
      }
    }
    result.m_outerIndex(m_blockOuterSize) = StorageIndex_(nnz);
    result.conservativeResizeBlockStorage_(nnz);
    return result;
  }

  template <bool Conjugate>
  BlockSparseMatrix<Scalar_, Options_, BlockCols_, BlockRows_, StorageIndex_> transposeImpl() const;

  // Resize both block storage arrays non-conservatively and update the capacity counter.
  void resizeBlockStorage_(Index n) {
    m_innerIndex.resize(n);
    m_values.resize(n * BlockSize);
  }

  void conservativeResizeBlockStorage_(Index n) {
    m_innerIndex.conservativeResize(n);
    m_values.conservativeResize(n * BlockSize);
  }

  // -------------------------------------------------------------------------
  // Storage
  // -------------------------------------------------------------------------
  Index m_blockOuterSize = 0;  // block-cols (ColMajor) or block-rows (RowMajor)
  Index m_blockInnerSize = 0;  // block-rows (ColMajor) or block-cols (RowMajor)

  Array<StorageIndex, Dynamic, 1> m_outerIndex =  // size: m_blockOuterSize + 1
      decltype(m_outerIndex)::Zero(m_blockOuterSize + 1);
  Array<StorageIndex, Dynamic, 1> m_innerIndex;
  // Block values stored consecutively; block k occupies
  // m_values[k*BlockSize .. (k+1)*BlockSize - 1].
  Array<Scalar, Dynamic, 1> m_values;

  // -------------------------------------------------------------------------
  // MultiInnerIterator
  //
  // Simultaneously walks BlockOuterSize_ consecutive outer vectors of a
  // compressed SparseMatrix, yielding scalar entries one at a time in
  // non-decreasing inner-index order.  On each step the sub-iterator with
  // the smallest current inner index is the "active" one.
  //
  // BlockOuterSize_ is BlockCols_ for ColMajor (one block-column at a time)
  // and BlockRows_ for RowMajor (one block-row at a time).
  // -------------------------------------------------------------------------
  template <typename SparseMatrixType>
  class MultiInnerIterator {
    using StorageIndex = typename SparseMatrixType::StorageIndex;
    static constexpr int BlockOuterSize_ = IsRowMajor ? BlockRows_ : BlockCols_;

   public:
    MultiInnerIterator(const SparseMatrixType& mat, Index outerBase)
        : m_outerPtr(mat.outerIndexPtr()),
          m_innerPtr(mat.innerIndexPtr()),
          m_valuePtr(mat.valuePtr()),
          m_outerBase(outerBase) {
      for (int k = 0; k < BlockOuterSize_; ++k) m_pos[k] = m_outerPtr[outerBase + k];
      advance();
    }

    EIGEN_STRONG_INLINE operator bool() const { return m_valid; }

    EIGEN_STRONG_INLINE MultiInnerIterator& operator++() {
      ++m_pos[m_active];
      advance();
      return *this;
    }

    // Absolute outer index of the current entry.
    EIGEN_STRONG_INLINE Index outer() const { return m_outerBase + m_active; }
    // Inner index of the current entry.
    EIGEN_STRONG_INLINE StorageIndex index() const { return m_innerPtr[m_pos[m_active]]; }
    // Scalar value of the current entry.
    EIGEN_STRONG_INLINE Scalar value() const { return m_valuePtr[m_pos[m_active]]; }

   private:
    // Find the sub-iterator with the smallest inner index.
    void advance() {
      m_valid = false;
      for (int k = 0; k < BlockOuterSize_; ++k) {
        if (m_pos[k] < m_outerPtr[m_outerBase + k + 1]) {
          if (!m_valid || m_innerPtr[m_pos[k]] < m_innerPtr[m_pos[m_active]]) {
            m_active = k;
            m_valid = true;
          }
        }
      }
    }

    const StorageIndex* m_outerPtr;
    const StorageIndex* m_innerPtr;
    const Scalar* m_valuePtr;
    Index m_outerBase;
    StorageIndex m_pos[BlockOuterSize_];
    int m_active = 0;
    bool m_valid = false;
  };

  // Allow other instantiations of BlockSparseMatrix to access private members
  // (needed by operator*).
  template <typename, int, int, int, typename>
  friend class BlockSparseMatrix;
  template <typename, int, bool>
  friend class BlockSparseTriangularView;
  template <typename, int, bool>
  friend class BlockSparseSelfAdjointView;
};

// =============================================================================
// Out-of-line method definitions
// =============================================================================

// -----------------------------------------------------------------------------
// setFromTriplets
// -----------------------------------------------------------------------------

template <typename Scalar_, int Options_, int BlockRows_, int BlockCols_, typename StorageIndex_>
template <typename InputIterator>
void BlockSparseMatrix<Scalar_, Options_, BlockRows_, BlockCols_, StorageIndex_>::setFromTriplets(InputIterator begin,
                                                                                                  InputIterator end) {
  Index n = static_cast<Index>(std::distance(begin, end));

  // Copy triplet coordinates and block values into Eigen arrays.
  Array<StorageIndex, Dynamic, 1> tOuter(n), tInner(n);
  Array<Scalar, Dynamic, 1> tValues(n * BlockSize);

  Index k = 0;
  for (InputIterator it = begin; it != end; ++it, ++k) {
    eigen_assert(it->row() >= 0 && it->row() < blockRows() && "setFromTriplets: block row out of range");
    eigen_assert(it->col() >= 0 && it->col() < blockCols() && "setFromTriplets: block col out of range");
    tOuter(k) = IsRowMajor ? StorageIndex(it->row()) : StorageIndex(it->col());
    tInner(k) = IsRowMajor ? StorageIndex(it->col()) : StorageIndex(it->row());
    BlockMap(tValues.data() + k * BlockSize) = it->value();
  }

  // Order the triplet indices by (outer, inner) using a stable LSD radix sort:
  // pass 1 buckets by inner, pass 2 by outer.  This runs in
  // O(n + blockInnerSize + blockOuterSize) with sequential writes, avoiding the
  // cache-unfriendly indirect comparison sort.
  Array<Index, Dynamic, 1> order(n), scratch(n);
  {
    Array<Index, Dynamic, 1> count = Array<Index, Dynamic, 1>::Zero(m_blockInnerSize + 1);
    for (Index i = 0; i < n; ++i) count(tInner(i) + 1)++;
    for (Index i = 0; i < m_blockInnerSize; ++i) count(i + 1) += count(i);
    for (Index i = 0; i < n; ++i) order(count(tInner(i))++) = i;
  }
  {
    Array<Index, Dynamic, 1> count = Array<Index, Dynamic, 1>::Zero(m_blockOuterSize + 1);
    for (Index i = 0; i < n; ++i) count(tOuter(i) + 1)++;
    for (Index i = 0; i < m_blockOuterSize; ++i) count(i + 1) += count(i);
    for (Index i = 0; i < n; ++i) {
      Index idx = order(i);
      scratch(count(tOuter(idx))++) = idx;
    }
    order.swap(scratch);
  }

  // Reset and pre-allocate (worst case: all n triplets are distinct blocks).
  m_outerIndex.resize(m_blockOuterSize + 1);
  m_outerIndex.setZero();
  resizeBlockStorage_(n);

  Index nnz = 0;
  k = 0;
  while (k < n) {
    Index pi = order(k);
    StorageIndex outer = tOuter(pi);
    StorageIndex inner = tInner(pi);

    BlockType block = ConstBlockMap(tValues.data() + pi * BlockSize);
    ++k;

    // Accumulate duplicate entries at the same (outer, inner) position.
    while (k < n) {
      Index pk = order(k);
      if (tOuter(pk) != outer || tInner(pk) != inner) break;
      block += ConstBlockMap(tValues.data() + pk * BlockSize);
      ++k;
    }

    m_innerIndex(nnz) = inner;
    blockRef(nnz) = block;
    m_outerIndex(outer + 1)++;
    ++nnz;
  }

  // Trim to actual number of unique blocks.
  conservativeResizeBlockStorage_(nnz);

  // Convert per-outer block counts to prefix sums.
  for (Index j = 0; j < m_blockOuterSize; ++j) {
    m_outerIndex(j + 1) += m_outerIndex(j);
  }
}

// -----------------------------------------------------------------------------
// toSparse
// -----------------------------------------------------------------------------

template <typename Scalar_, int Options_, int BlockRows_, int BlockCols_, typename StorageIndex_>
SparseMatrix<Scalar_, Options_, StorageIndex_>
BlockSparseMatrix<Scalar_, Options_, BlockRows_, BlockCols_, StorageIndex_>::toSparse() const {
  SparseMatrix<Scalar_, Options_, StorageIndex_> result(rows(), cols());
  result.reserve(nonZeroBlocks() * BlockSize);

  if (!IsRowMajor) {
    // ColMajor: outer = block-column j.  Emit scalar columns j*BlockCols+c
    // in order c = 0..BlockCols-1.  Within each scalar column, blocks are
    // sorted by bi (block-row), so scalar rows bi*BlockRows+r are increasing.
    for (Index j = 0; j < m_blockOuterSize; ++j) {
      for (Index c = 0; c < BlockCols_; ++c) {
        result.startVec(j * BlockCols_ + c);
        for (Index id = m_outerIndex(j); id < m_outerIndex(j + 1); ++id) {
          Index bi = m_innerIndex(id);
          ConstBlockMap blk = blockRef(id);
          for (Index r = 0; r < BlockRows_; ++r) {
            result.insertBack(bi * BlockRows_ + r, j * BlockCols_ + c) = blk(r, c);
          }
        }
      }
    }
  } else {
    // RowMajor: outer = block-row bi.  Emit scalar rows bi*BlockRows+r
    // in order r = 0..BlockRows-1.  Within each scalar row, blocks are
    // sorted by j (block-col), so scalar cols j*BlockCols+c are increasing.
    for (Index bi = 0; bi < m_blockOuterSize; ++bi) {
      for (Index r = 0; r < BlockRows_; ++r) {
        result.startVec(bi * BlockRows_ + r);
        for (Index id = m_outerIndex(bi); id < m_outerIndex(bi + 1); ++id) {
          Index j = m_innerIndex(id);
          ConstBlockMap blk = blockRef(id);
          for (Index c = 0; c < BlockCols_; ++c) {
            result.insertBack(bi * BlockRows_ + r, j * BlockCols_ + c) = blk(r, c);
          }
        }
      }
    }
  }

  result.finalize();
  return result;
}

// -----------------------------------------------------------------------------
// fromSparse
// -----------------------------------------------------------------------------

template <typename Scalar_, int Options_, int BlockRows_, int BlockCols_, typename StorageIndex_>
BlockSparseMatrix<Scalar_, Options_, BlockRows_, BlockCols_, StorageIndex_>
BlockSparseMatrix<Scalar_, Options_, BlockRows_, BlockCols_, StorageIndex_>::fromSparse(
    const SparseMatrix<Scalar_, Options_, StorageIndex_>& sp) {
  eigen_assert(sp.rows() % BlockRows_ == 0 && "matrix rows not divisible by BlockRows");
  eigen_assert(sp.cols() % BlockCols_ == 0 && "matrix cols not divisible by BlockCols");
  eigen_assert(sp.isCompressed() && "fromSparse requires a compressed SparseMatrix");

  Index bRows = sp.rows() / BlockRows_;
  Index bCols = sp.cols() / BlockCols_;

  // BlockOuterSize: how many consecutive outer vectors form one block-outer strip.
  // BlockInnerSize: the inner dimension of each block.
  constexpr Index BlockOuterSize = IsRowMajor ? BlockRows_ : BlockCols_;
  constexpr Index BlockInnerSize = IsRowMajor ? BlockCols_ : BlockRows_;
  constexpr StorageIndex_ kEmptyIndex = -1;

  using SpMat = SparseMatrix<Scalar_, Options_, StorageIndex_>;

  BlockSparseMatrix result(bRows, bCols);

  // Pass 1: count the number of unique block-inner indices per block-outer,
  // by scanning each group of BlockOuterSize consecutive outer vectors together.
  for (Index outerBlock = 0; outerBlock < result.m_blockOuterSize; ++outerBlock) {
    StorageIndex_ prevInnerBlock = kEmptyIndex;
    for (MultiInnerIterator<SpMat> it(sp, outerBlock * BlockOuterSize); it; ++it) {
      StorageIndex_ innerBlock = it.index() / StorageIndex_(BlockInnerSize);
      if (innerBlock != prevInnerBlock) {
        result.m_outerIndex(outerBlock + 1)++;
        prevInnerBlock = innerBlock;
      }
    }
  }

  // Prefix sum → result.m_outerIndex becomes the standard CSC/CSR outer pointer.
  for (Index j = 0; j < result.m_blockOuterSize; ++j) result.m_outerIndex(j + 1) += result.m_outerIndex(j);

  Index nBlocks = result.m_outerIndex(result.m_blockOuterSize);
  result.resizeBlockStorage_(nBlocks);
  result.m_values.setZero();

  // Pass 2: scatter each scalar entry directly into its position within the
  // pre-zeroed block value array.
  for (Index outerBlock = 0; outerBlock < result.m_blockOuterSize; ++outerBlock) {
    Index blockId = result.m_outerIndex(outerBlock) - 1;  // incremented on first new block
    StorageIndex_ prevInnerBlock = kEmptyIndex;

    for (MultiInnerIterator<SpMat> it(sp, outerBlock * BlockOuterSize); it; ++it) {
      Index absOuter = it.outer();          // absolute outer index in sp
      StorageIndex_ innerIdx = it.index();  // inner index in sp
      StorageIndex_ innerBlock = innerIdx / StorageIndex_(BlockInnerSize);

      if (innerBlock != prevInnerBlock) {
        ++blockId;
        result.m_innerIndex(blockId) = innerBlock;
        prevInnerBlock = innerBlock;
      }

      // Scatter into block storage (layout matches the BSM's Options_).
      // ColMajor blocks: col * BlockRows_ + row = localOuter * BlockInnerSize + localInner
      // RowMajor blocks: row * BlockCols_ + col = localOuter * BlockInnerSize + localInner
      Index localOuter = absOuter % BlockOuterSize;
      Index localInner = innerIdx % BlockInnerSize;
      Index offset = localOuter * BlockInnerSize + localInner;

      result.m_values(blockId * BlockSize + offset) = it.value();
    }
  }

  return result;
}

// -----------------------------------------------------------------------------
// operator* (block-sparse product)
// -----------------------------------------------------------------------------

template <typename Scalar_, int Options_, int BlockRows_, int BlockCols_, typename StorageIndex_>
template <int RhsBlockCols>
BlockSparseMatrix<Scalar_, Options_, BlockRows_, RhsBlockCols, StorageIndex_>
BlockSparseMatrix<Scalar_, Options_, BlockRows_, BlockCols_, StorageIndex_>::operator*(
    const BlockSparseMatrix<Scalar_, Options_, BlockCols_, RhsBlockCols, StorageIndex_>& rhs) const {
  using RhsMatrix = BlockSparseMatrix<Scalar_, Options_, BlockCols_, RhsBlockCols, StorageIndex_>;
  using ResultMatrix = BlockSparseMatrix<Scalar_, Options_, BlockRows_, RhsBlockCols, StorageIndex_>;
  using ResultBlock = Matrix<Scalar_, BlockRows_, RhsBlockCols, Options_>;
  constexpr int ResultBlockSize = BlockRows_ * RhsBlockCols;

  eigen_assert(blockCols() == rhs.blockRows() && "BlockSparseMatrix product: lhs.blockCols() != rhs.blockRows()");

  Index cBlockRows = blockRows();
  Index cBlockCols = rhs.blockCols();
  ResultMatrix result(cBlockRows, cBlockCols);

  // For ColMajor: mask / accum indexed by block-row (size = cBlockRows).
  // For RowMajor: mask / accum indexed by block-col (size = cBlockCols).
  Index maskSize = IsRowMajor ? cBlockCols : cBlockRows;
  Array<uint8_t, Dynamic, 1> mask = Array<uint8_t, Dynamic, 1>::Zero(maskSize);
  Array<Scalar_, Dynamic, 1> accumData(maskSize * ResultBlockSize);
  Array<Index, Dynamic, 1> indices(maskSize);
  Index nIndices = 0;

  // Grow result storage geometrically rather than pre-allocating the dense
  // worst case (cBlockRows*cBlockCols blocks): the product is typically far
  // sparser than that, so the dense bound would blow up peak memory.
  // capacity is always kept <= maxResultNnz (the true upper bound), and since
  // any single outer emits at most maskSize blocks, the initial estimate of
  // maskSize guarantees the first outer fits before the first grow check.
  Index cOuterSize = result.m_blockOuterSize;
  Index maxResultNnz = cBlockRows * cBlockCols;
  Index capacity = numext::mini(maxResultNnz, numext::maxi(maskSize, nonZeroBlocks() + rhs.nonZeroBlocks()));
  result.resizeBlockStorage_(capacity);
  Index nnz = 0;

  for (Index out = 0; out < cOuterSize; ++out) {
    result.m_outerIndex(out) = StorageIndex_(nnz);

    if (!IsRowMajor) {
      // ColMajor: out is block-column j of the result.
      // For each block B(k,j) and each block A(bi,k): C(bi,j) += A(bi,k)*B(k,j).
      Index j = out;
      for (Index rhsId = rhs.m_outerIndex(j); rhsId < rhs.m_outerIndex(j + 1); ++rhsId) {
        Index k = rhs.m_innerIndex(rhsId);
        typename RhsMatrix::ConstBlockMap Bkj = rhs.blockRef(rhsId);
        for (Index lhsId = m_outerIndex(k); lhsId < m_outerIndex(k + 1); ++lhsId) {
          Index bi = m_innerIndex(lhsId);
          if (!mask(bi)) {
            mask(bi) = 1;
            Map<ResultBlock>(accumData.data() + bi * ResultBlockSize).noalias() = blockRef(lhsId) * Bkj;
            indices(nIndices++) = bi;
          } else {
            Map<ResultBlock>(accumData.data() + bi * ResultBlockSize).noalias() += blockRef(lhsId) * Bkj;
          }
        }
      }
    } else {
      // RowMajor: out is block-row bi of the result.
      // For each block A(bi,k) and each block B(k,j): C(bi,j) += A(bi,k)*B(k,j).
      Index bi = out;
      for (Index lhsId = m_outerIndex(bi); lhsId < m_outerIndex(bi + 1); ++lhsId) {
        Index k = m_innerIndex(lhsId);
        ConstBlockMap Aik = blockRef(lhsId);
        for (Index rhsId = rhs.m_outerIndex(k); rhsId < rhs.m_outerIndex(k + 1); ++rhsId) {
          Index j = rhs.m_innerIndex(rhsId);
          if (!mask(j)) {
            mask(j) = 1;
            Map<ResultBlock>(accumData.data() + j * ResultBlockSize).noalias() = Aik * rhs.blockRef(rhsId);
            indices(nIndices++) = j;
          } else {
            Map<ResultBlock>(accumData.data() + j * ResultBlockSize).noalias() += Aik * rhs.blockRef(rhsId);
          }
        }
      }
    }

    // Sort the accumulated indices so the result's inner index array is sorted.
    std::sort(indices.data(), indices.data() + nIndices);
    if (nnz + nIndices > capacity) {
      capacity = numext::mini(maxResultNnz, numext::maxi(2 * capacity, nnz + nIndices));
      result.conservativeResizeBlockStorage_(capacity);
    }
    for (Index ki = 0; ki < nIndices; ++ki) {
      Index idx = indices(ki);
      result.m_innerIndex(nnz) = StorageIndex_(idx);
      result.blockRef(nnz) = Map<ResultBlock>(accumData.data() + idx * ResultBlockSize);
      mask(idx) = 0;
      ++nnz;
    }
    nIndices = 0;
  }
  result.m_outerIndex(cOuterSize) = StorageIndex_(nnz);

  // Trim to actual number of result blocks.
  result.conservativeResizeBlockStorage_(nnz);

  return result;
}

// =============================================================================
// BlockSparseTriangularView
// =============================================================================

/** \class BlockSparseTriangularView
 * \ingroup SparseCore_Module
 * \brief Lazy block-level triangular view of a BlockSparseMatrix.
 *
 * Obtained via \c BSM::triangularView<Mode>() or
 * \c BSM::triangularView<Mode, true>() (the latter asserts that the unused
 * triangle of every diagonal block is already zero, enabling faster vectorised
 * products without an explicit triangularView on the block).
 *
 * By default (\p DiagIsTriangular = \c false) diagonal blocks are treated as
 * triangular regardless of what is stored in the unused half: products call
 * \c block.triangularView<DiagMode>() and \c eval() zeroes the unused
 * triangle.  Set \p DiagIsTriangular = \c true to skip that overhead when you
 * can guarantee the unused triangle is already zero.
 */
template <typename BSM, int Mode, bool DiagIsTriangular = false>
class BlockSparseTriangularView {
 public:
  using Scalar = typename BSM::Scalar;
  using StorageIndex = typename BSM::StorageIndex;
  using BlockType = typename BSM::BlockType;
  using BlockMap = typename BSM::BlockMap;
  using ConstBlockMap = typename BSM::ConstBlockMap;
  static constexpr int BlockRows = BSM::BlockRows;
  static constexpr int BlockCols = BSM::BlockCols;
  static constexpr int BlockSize = BSM::BlockSize;
  static constexpr bool IsRowMajor = BSM::IsRowMajor;
  static constexpr bool IsUpper = (Mode & Upper) != 0;

  explicit BlockSparseTriangularView(const BSM& m) : m_matrix(m) {}

  Index rows() const { return m_matrix.rows(); }
  Index cols() const { return m_matrix.cols(); }

  // ---- Materialize ---------------------------------------------------------

  /** Copy the triangular blocks into a new BSM; off-triangle blocks are dropped.
   *  When DiagIsTriangular is false the unused triangle of each diagonal block
   *  is explicitly zeroed in the output. */
  BSM eval() const {
    constexpr int ZeroMode = IsUpper ? StrictlyLower : StrictlyUpper;
    const BSM& m = m_matrix;
    BSM result(m.blockRows(), m.blockCols());
    result.resizeBlockStorage_(m.nonZeroBlocks());
    Index nnz = 0;

    for (Index out = 0; out < m.m_blockOuterSize; ++out) {
      result.m_outerIndex(out) = StorageIndex(nnz);
      for (Index id = m.m_outerIndex(out); id < m.m_outerIndex(out + 1); ++id) {
        Index inner = m.m_innerIndex(id);
        Index bi = IsRowMajor ? out : inner;
        Index bj = IsRowMajor ? inner : out;
        if (IsUpper ? (bj < bi) : (bj > bi)) continue;
        result.m_innerIndex(nnz) = StorageIndex(inner);
        result.m_values.template segment<BlockSize>(nnz * BlockSize) =
            m.m_values.template segment<BlockSize>(id * BlockSize);
        EIGEN_IF_CONSTEXPR (!DiagIsTriangular) {
          if (bi == bj)
            BlockMap(result.m_values.data() + nnz * BlockSize).template triangularView<ZeroMode>().setZero();
        }
        ++nnz;
      }
    }
    result.m_outerIndex(m.m_blockOuterSize) = StorageIndex(nnz);
    result.conservativeResizeBlockStorage_(nnz);
    return result;
  }

  /** Convert to a scalar-level SparseMatrix (off-triangle blocks zeroed). */
  SparseMatrix<Scalar, BSM::Options, StorageIndex> toSparse() const { return eval().toSparse(); }

  // ---- Arithmetic ----------------------------------------------------------

  BSM operator+(const BlockSparseTriangularView& other) const { return eval() + other.eval(); }
  BSM operator-(const BlockSparseTriangularView& other) const { return eval() - other.eval(); }

  /** Tri × BSM product (materialises this view then delegates). */
  template <int RhsBlockCols>
  BlockSparseMatrix<Scalar, BSM::Options, BlockRows, RhsBlockCols, StorageIndex> operator*(
      const BlockSparseMatrix<Scalar, BSM::Options, BlockCols, RhsBlockCols, StorageIndex>& rhs) const {
    return eval() * rhs;
  }

  // ---- Dense products (no intermediate materialisation) --------------------

  template <typename OtherDerived>
  Matrix<Scalar, Dynamic, OtherDerived::ColsAtCompileTime> operator*(const MatrixBase<OtherDerived>& rhs) const {
    EIGEN_STATIC_ASSERT(
        (std::is_same<Scalar, typename OtherDerived::Scalar>::value),
        YOU_MIXED_DIFFERENT_NUMERIC_TYPES__YOU_NEED_TO_USE_THE_CAST_METHOD_OF_MATRIXBASE_TO_CAST_NUMERIC_TYPES_EXPLICITLY)
    eigen_assert(m_matrix.cols() == rhs.rows() && "BlockSparseTriangularView * Dense: dimension mismatch");
    using ResultType = Matrix<Scalar, Dynamic, OtherDerived::ColsAtCompileTime>;
    ResultType result = ResultType::Zero(m_matrix.rows(), rhs.cols());
    for (Index out = 0; out < m_matrix.m_blockOuterSize; ++out) {
      for (Index id = m_matrix.m_outerIndex(out); id < m_matrix.m_outerIndex(out + 1); ++id) {
        Index inner = m_matrix.m_innerIndex(id);
        Index bi = IsRowMajor ? out : inner;
        Index bj = IsRowMajor ? inner : out;
        if (IsUpper ? (bj < bi) : (bj > bi)) continue;
        constexpr int DiagMode = IsUpper ? Upper : Lower;
        if (!DiagIsTriangular && bi == bj)
          result.template middleRows<BlockRows>(bi * BlockRows).noalias() +=
              m_matrix.blockRef(id).template triangularView<DiagMode>() *
              rhs.template middleRows<BlockCols>(bj * BlockCols);
        else
          result.template middleRows<BlockRows>(bi * BlockRows).noalias() +=
              m_matrix.blockRef(id) * rhs.template middleRows<BlockCols>(bj * BlockCols);
      }
    }
    return result;
  }

  template <typename OtherDerived>
  friend Matrix<Scalar, OtherDerived::RowsAtCompileTime, Dynamic> operator*(const MatrixBase<OtherDerived>& lhs,
                                                                            const BlockSparseTriangularView& tri) {
    EIGEN_STATIC_ASSERT(
        (std::is_same<Scalar, typename OtherDerived::Scalar>::value),
        YOU_MIXED_DIFFERENT_NUMERIC_TYPES__YOU_NEED_TO_USE_THE_CAST_METHOD_OF_MATRIXBASE_TO_CAST_NUMERIC_TYPES_EXPLICITLY)
    eigen_assert(lhs.cols() == tri.m_matrix.rows() && "Dense * BlockSparseTriangularView: dimension mismatch");
    constexpr bool isRM = BSM::IsRowMajor;
    using ResultType = Matrix<Scalar, OtherDerived::RowsAtCompileTime, Dynamic>;
    ResultType result = ResultType::Zero(lhs.rows(), tri.m_matrix.cols());
    for (Index out = 0; out < tri.m_matrix.m_blockOuterSize; ++out) {
      for (Index id = tri.m_matrix.m_outerIndex(out); id < tri.m_matrix.m_outerIndex(out + 1); ++id) {
        Index inner = tri.m_matrix.m_innerIndex(id);
        Index bi = isRM ? out : inner;
        Index bj = isRM ? inner : out;
        if (IsUpper ? (bj < bi) : (bj > bi)) continue;
        constexpr int DiagMode = IsUpper ? Upper : Lower;
        if (!DiagIsTriangular && bi == bj)
          result.template middleCols<BlockCols>(bj * BlockCols).noalias() +=
              lhs.template middleCols<BlockRows>(bi * BlockRows) *
              tri.m_matrix.blockRef(id).template triangularView<DiagMode>();
        else
          result.template middleCols<BlockCols>(bj * BlockCols).noalias() +=
              lhs.template middleCols<BlockRows>(bi * BlockRows) * tri.m_matrix.blockRef(id);
      }
    }
    return result;
  }

  // ---- Triangular solve -------------------------------------------------------

  /** Solve T * x = rhs in-place. Requires square blocks.
   *  ColMajor Lower: forward sub, diagonal first per column.
   *  ColMajor Upper: backward sub, diagonal last per column.
   *  RowMajor Lower: forward sub, diagonal last per row.
   *  RowMajor Upper: backward sub, diagonal first per row.
   */
  template <typename Derived>
  void solveInPlace(MatrixBase<Derived>& x) const {
    doSolveImpl<false, false>(x.derived());
  }

  /** Proxy returned by transpose(): solveInPlace solves T^T x = b. */
  struct TransposeReturnType {
    const BlockSparseTriangularView& m_tri;
    template <typename Derived>
    void solveInPlace(MatrixBase<Derived>& x) const {
      m_tri.template doSolveImpl<true, false>(x.derived());
    }
  };

  /** Proxy returned by adjoint(): solveInPlace solves T^H x = b. */
  struct AdjointReturnType {
    const BlockSparseTriangularView& m_tri;
    template <typename Derived>
    void solveInPlace(MatrixBase<Derived>& x) const {
      m_tri.template doSolveImpl<true, true>(x.derived());
    }
  };

  TransposeReturnType transpose() const { return {*this}; }
  AdjointReturnType adjoint() const { return {*this}; }

 private:
  const BSM& m_matrix;

  // Non-transposed solve for both storage orders.
  //
  // diagFirst = (IsUpper == IsRowMajor): ColMajor Lower→first, ColMajor Upper→last,
  //                                      RowMajor Lower→last,  RowMajor Upper→first.
  // Loop direction: forward for Lower, backward for Upper (same for both storage orders).
  // ColMajor: solve diagonal first, then scatter  x[inner] -= blk * x[k].
  // RowMajor: gather x[k] -= blk * x[inner] first, then solve diagonal.
  template <typename Derived>
  void doSolveDirect(Derived& x) const {
    EIGEN_STATIC_ASSERT(BlockRows == BlockCols, THIS_METHOD_IS_ONLY_FOR_SQUARE_BLOCK_MATRICES)
    constexpr int DiagMode = IsUpper ? Upper : Lower;
    constexpr bool diagFirst = (IsUpper == BSM::IsRowMajor);
    Index nb = m_matrix.blockCols();  // == blockRows() for square matrices
    eigen_assert(x.rows() == m_matrix.rows() && "solveInPlace: size mismatch");

    const StorageIndex* innerPtr = m_matrix.innerIndexPtr();
    const StorageIndex* outerPtr = m_matrix.outerIndexPtr();

    Index outerStart = IsUpper ? nb - 1 : 0;
    Index outerEnd = IsUpper ? -1 : nb;
    constexpr Index kStep = IsUpper ? -1 : 1;

    for (Index k = outerStart; k != outerEnd; k += kStep) {
      const StorageIndex* beg = innerPtr + outerPtr[k];
      const StorageIndex* end = innerPtr + outerPtr[k + 1];
      if (beg == end) continue;
      const StorageIndex* diag_ptr = diagFirst ? beg : end - 1;
      const StorageIndex* off_beg = diagFirst ? beg + 1 : beg;
      const StorageIndex* off_end = diagFirst ? end : end - 1;
      eigen_assert(*diag_ptr == k);
      EIGEN_IF_CONSTEXPR (!BSM::IsRowMajor) {
        m_matrix.blockRef(diag_ptr - innerPtr)
            .template triangularView<DiagMode>()
            .solveInPlace(x.template middleRows<BlockRows>(k * BlockRows));
        for (const StorageIndex* it = off_beg; it != off_end; ++it)
          x.template middleRows<BlockRows>(*it * BlockRows).noalias() -=
              m_matrix.blockRef(it - innerPtr) * x.template middleRows<BlockRows>(k * BlockRows);
      } else {
        for (const StorageIndex* it = off_beg; it != off_end; ++it)
          x.template middleRows<BlockRows>(k * BlockRows).noalias() -=
              m_matrix.blockRef(it - innerPtr) * x.template middleRows<BlockRows>(*it * BlockRows);
        m_matrix.blockRef(diag_ptr - innerPtr)
            .template triangularView<DiagMode>()
            .solveInPlace(x.template middleRows<BlockRows>(k * BlockRows));
      }
    }
  }

  // Transposed/adjoint solve for both storage orders.
  //
  // Loop direction: forward for Upper, backward for Lower (same for both storage orders).
  // ColMajor: gather x[k] -= adj(blk) * x[inner] first, then solve adj(diagonal).
  // RowMajor: solve adj(diagonal) first, then scatter x[inner] -= adj(blk) * x[k].
  template <bool Conjugate, typename Derived>
  void doSolveTransposed(Derived& x) const {
    EIGEN_STATIC_ASSERT(BlockRows == BlockCols, THIS_METHOD_IS_ONLY_FOR_SQUARE_BLOCK_MATRICES)
    constexpr int DiagMode = IsUpper ? Upper : Lower;
    constexpr bool diagFirst = (IsUpper == BSM::IsRowMajor);
    Index nb = m_matrix.blockCols();  // == blockRows() for square matrices
    eigen_assert(x.rows() == m_matrix.rows() && "solveInPlace: size mismatch");

    const StorageIndex* innerPtr = m_matrix.innerIndexPtr();
    const StorageIndex* outerPtr = m_matrix.outerIndexPtr();

    Index outerStart = IsUpper ? 0 : nb - 1;
    Index outerEnd = IsUpper ? nb : -1;
    constexpr Index kStep = IsUpper ? 1 : -1;

    for (Index k = outerStart; k != outerEnd; k += kStep) {
      const StorageIndex* beg = innerPtr + outerPtr[k];
      const StorageIndex* end = innerPtr + outerPtr[k + 1];
      if (beg == end) continue;
      const StorageIndex* diag_ptr = diagFirst ? beg : end - 1;
      const StorageIndex* off_beg = diagFirst ? beg + 1 : beg;
      const StorageIndex* off_end = diagFirst ? end : end - 1;
      eigen_assert(*diag_ptr == k);
      EIGEN_IF_CONSTEXPR (!BSM::IsRowMajor) {
        for (const StorageIndex* it = off_beg; it != off_end; ++it)
          x.template middleRows<BlockRows>(k * BlockRows).noalias() -=
              internal::adjoint_if<Conjugate>(m_matrix.blockRef(it - innerPtr)) *
              x.template middleRows<BlockRows>(*it * BlockRows);
        internal::adjoint_if<Conjugate>(m_matrix.blockRef(diag_ptr - innerPtr).template triangularView<DiagMode>())
            .solveInPlace(x.template middleRows<BlockRows>(k * BlockRows));
      } else {
        internal::adjoint_if<Conjugate>(m_matrix.blockRef(diag_ptr - innerPtr).template triangularView<DiagMode>())
            .solveInPlace(x.template middleRows<BlockRows>(k * BlockRows));
        for (const StorageIndex* it = off_beg; it != off_end; ++it)
          x.template middleRows<BlockRows>(*it * BlockRows).noalias() -=
              internal::adjoint_if<Conjugate>(m_matrix.blockRef(it - innerPtr)) *
              x.template middleRows<BlockRows>(k * BlockRows);
      }
    }
  }

  template <bool Transposed, bool Conjugate, typename Derived>
  void doSolveImpl(Derived& x) const {
    EIGEN_IF_CONSTEXPR (!Transposed)
      doSolveDirect(x);
    else
      doSolveTransposed<Conjugate>(x);
  }
};

// =============================================================================
// BlockSparseSelfAdjointView
// =============================================================================

/** \class BlockSparseSelfAdjointView
 * \ingroup SparseCore_Module
 * \brief Lazy block-level self-adjoint (Hermitian) view of a BlockSparseMatrix.
 *
 * Obtained via \c BSM::selfadjointView<UpLo>() or
 * \c BSM::selfadjointView<UpLo, true>() (the latter signals that every
 * diagonal block is itself Hermitian, enabling DSYMM/ZHEMM on those blocks).
 *
 * Only the triangle selected by \p UpLo is read; the other triangle is
 * reconstructed on-the-fly as the adjoint of each stored off-diagonal block.
 *
 * \pre BSM::BlockRows == BSM::BlockCols.
 */
template <typename BSM, int UpLo, bool DiagIsSelfAdjoint>
class BlockSparseSelfAdjointView {
 public:
  using Scalar = typename BSM::Scalar;
  using StorageIndex = typename BSM::StorageIndex;
  using BlockType = typename BSM::BlockType;
  using BlockMap = typename BSM::BlockMap;
  using ConstBlockMap = typename BSM::ConstBlockMap;
  static constexpr int BlockRows = BSM::BlockRows;  // == BlockCols
  static constexpr int BlockCols = BSM::BlockCols;
  static constexpr int BlockSize = BSM::BlockSize;
  static constexpr bool IsRowMajor = BSM::IsRowMajor;
  static constexpr bool IsUpper = (UpLo & Upper) != 0;
  // UpLo passed to Eigen's dense selfadjointView on diagonal blocks:
  static constexpr int DiagUpLo = IsUpper ? Upper : Lower;

  explicit BlockSparseSelfAdjointView(const BSM& m) : m_matrix(m) {}

  Index rows() const { return m_matrix.rows(); }
  Index cols() const { return m_matrix.cols(); }

  // ---- Materialize ---------------------------------------------------------

  /** Build a full symmetric BSM: stored triangle + adjoint mirror of each
   *  off-diagonal block.  Diagonal blocks: when DiagIsSelfAdjoint is false,
   *  both triangles are filled from the stored triangle via selfadjointView;
   *  when true the block is already fully populated and is copied as-is. */
  BSM eval() const {
    const BSM& m = m_matrix;

    Index nDiag = 0, nOff = 0;
    for (Index out = 0; out < m.m_blockOuterSize; ++out)
      for (Index id = m.m_outerIndex(out); id < m.m_outerIndex(out + 1); ++id) {
        Index inner = m.m_innerIndex(id);
        Index bi = IsRowMajor ? out : inner;
        Index bj = IsRowMajor ? inner : out;
        if (IsUpper ? (bj < bi) : (bj > bi)) continue;
        if (bi == bj)
          ++nDiag;
        else
          ++nOff;
      }

    Index nTotal = nDiag + 2 * nOff;

    Array<StorageIndex, Dynamic, 1> brows(nTotal), bcols(nTotal);
    Array<Scalar, Dynamic, 1> bvals(nTotal * BlockSize);

    Index k = 0;
    for (Index out = 0; out < m.m_blockOuterSize; ++out)
      for (Index id = m.m_outerIndex(out); id < m.m_outerIndex(out + 1); ++id) {
        Index inner = m.m_innerIndex(id);
        Index bi = IsRowMajor ? out : inner;
        Index bj = IsRowMajor ? inner : out;
        if (IsUpper ? (bj < bi) : (bj > bi)) continue;

        brows(k) = StorageIndex(bi);
        bcols(k) = StorageIndex(bj);
        if (!DiagIsSelfAdjoint && bi == bj)
          BlockMap(bvals.data() + k * BlockSize) = m.blockRef(id).template selfadjointView<DiagUpLo>();
        else
          BlockMap(bvals.data() + k * BlockSize) = m.blockRef(id);
        ++k;

        if (bi != bj) {
          brows(k) = StorageIndex(bj);
          bcols(k) = StorageIndex(bi);
          BlockMap(bvals.data() + k * BlockSize) = m.blockRef(id).adjoint();
          ++k;
        }
      }

    // Sort by (outer, inner) then build BSM directly (no duplicates by construction).
    Array<Index, Dynamic, 1> perm(nTotal);
    std::iota(perm.data(), perm.data() + nTotal, Index(0));
    std::sort(perm.data(), perm.data() + nTotal, [&](Index a, Index b) {
      StorageIndex ao = IsRowMajor ? brows(a) : bcols(a);
      StorageIndex bo = IsRowMajor ? brows(b) : bcols(b);
      if (ao != bo) return ao < bo;
      return (IsRowMajor ? bcols(a) : brows(a)) < (IsRowMajor ? bcols(b) : brows(b));
    });

    BSM result(m.blockRows(), m.blockCols());
    result.resizeBlockStorage_(nTotal);

    for (Index ki = 0; ki < nTotal; ++ki) {
      Index pi = perm(ki);
      StorageIndex outer = IsRowMajor ? brows(pi) : bcols(pi);
      StorageIndex inner = IsRowMajor ? bcols(pi) : brows(pi);
      result.m_outerIndex(outer + 1)++;
      result.m_innerIndex(ki) = inner;
      result.m_values.template segment<BlockSize>(ki * BlockSize) = bvals.template segment<BlockSize>(pi * BlockSize);
    }
    for (Index j = 0; j < result.m_blockOuterSize; ++j) result.m_outerIndex(j + 1) += result.m_outerIndex(j);

    return result;
  }

  /** Convert to a symmetrised scalar-level SparseMatrix. */
  SparseMatrix<Scalar, BSM::Options, StorageIndex> toSparse() const { return eval().toSparse(); }

  // ---- Arithmetic ----------------------------------------------------------

  BSM operator+(const BlockSparseSelfAdjointView& other) const { return eval() + other.eval(); }
  BSM operator-(const BlockSparseSelfAdjointView& other) const { return eval() - other.eval(); }

  /** SelfAdj × BSM: materialises the view then uses the general SpGEMM. */
  template <int RhsBlockCols>
  BlockSparseMatrix<Scalar, BSM::Options, BlockRows, RhsBlockCols, StorageIndex> operator*(
      const BlockSparseMatrix<Scalar, BSM::Options, BlockCols, RhsBlockCols, StorageIndex>& rhs) const {
    return eval() * rhs;
  }

  // ---- Dense products (no materialisation; exploits both triangles) ---------

  /** SelfAdj × Dense.
   *
   *  Off-diagonal stored block A(bi,bj) contributes:
   *    result(bi) += A(bi,bj) * rhs(bj)        [stored triangle]
   *    result(bj) += A(bi,bj)^H * rhs(bi)      [implicit mirror]
   *
   *  When DiagIsSelfAdjoint is true, both triangles of each diagonal block
   *  are valid; a plain product is used (faster for small fixed-size blocks).
   *  Otherwise only the UpLo triangle is assumed valid and selfadjointView
   *  reconstructs the full diagonal-block product.
   */
  template <typename OtherDerived>
  Matrix<Scalar, Dynamic, OtherDerived::ColsAtCompileTime> operator*(const MatrixBase<OtherDerived>& rhs) const {
    EIGEN_STATIC_ASSERT(
        (std::is_same<Scalar, typename OtherDerived::Scalar>::value),
        YOU_MIXED_DIFFERENT_NUMERIC_TYPES__YOU_NEED_TO_USE_THE_CAST_METHOD_OF_MATRIXBASE_TO_CAST_NUMERIC_TYPES_EXPLICITLY)
    eigen_assert(m_matrix.cols() == rhs.rows() && "BlockSparseSelfAdjointView * Dense: dimension mismatch");
    using ResultType = Matrix<Scalar, Dynamic, OtherDerived::ColsAtCompileTime>;
    ResultType result = ResultType::Zero(m_matrix.rows(), rhs.cols());

    for (Index out = 0; out < m_matrix.m_blockOuterSize; ++out)
      for (Index id = m_matrix.m_outerIndex(out); id < m_matrix.m_outerIndex(out + 1); ++id) {
        Index inner = m_matrix.m_innerIndex(id);
        Index bi = IsRowMajor ? out : inner;
        Index bj = IsRowMajor ? inner : out;
        if (IsUpper ? (bj < bi) : (bj > bi)) continue;

        if (bi == bj) {
          EIGEN_IF_CONSTEXPR (DiagIsSelfAdjoint) {
            result.template middleRows<BlockRows>(bi * BlockRows).noalias() +=
                m_matrix.blockRef(id) * rhs.template middleRows<BlockCols>(bj * BlockCols);
          } else {
            // Materialize the tiny diagonal block as a fixed-size Hermitian matrix, then use the
            // ordinary (coeff-based for a vector rhs) product.  This avoids the runtime-sized,
            // EIGEN_DONT_INLINE selfadjoint_matrix_vector_product kernel, which is tuned for large
            // matrices and is pure overhead for a 2-4 row block.
            BlockType diag = m_matrix.blockRef(id).template selfadjointView<DiagUpLo>();
            result.template middleRows<BlockRows>(bi * BlockRows).noalias() +=
                diag * rhs.template middleRows<BlockCols>(bj * BlockCols);
          }
        } else {
          result.template middleRows<BlockRows>(bi * BlockRows).noalias() +=
              m_matrix.blockRef(id) * rhs.template middleRows<BlockCols>(bj * BlockCols);
          result.template middleRows<BlockRows>(bj * BlockRows).noalias() +=
              m_matrix.blockRef(id).adjoint() * rhs.template middleRows<BlockRows>(bi * BlockRows);
        }
      }
    return result;
  }

  /** Dense × SelfAdj:  lhs * A == (A^H * lhs^H)^H == (A * lhs^H)^H for Hermitian A. */
  template <typename OtherDerived>
  friend Matrix<Scalar, OtherDerived::RowsAtCompileTime, Dynamic> operator*(const MatrixBase<OtherDerived>& lhs,
                                                                            const BlockSparseSelfAdjointView& view) {
    EIGEN_STATIC_ASSERT(
        (std::is_same<Scalar, typename OtherDerived::Scalar>::value),
        YOU_MIXED_DIFFERENT_NUMERIC_TYPES__YOU_NEED_TO_USE_THE_CAST_METHOD_OF_MATRIXBASE_TO_CAST_NUMERIC_TYPES_EXPLICITLY)
    return (view * lhs.adjoint()).adjoint();
  }

 private:
  const BSM& m_matrix;
};

// =============================================================================
// BlockSparseMatrix::transposeImpl / transpose / adjoint (out-of-line)
// =============================================================================

template <typename Scalar_, int Options_, int BlockRows_, int BlockCols_, typename StorageIndex_>
template <bool Conjugate>
BlockSparseMatrix<Scalar_, Options_, BlockCols_, BlockRows_, StorageIndex_>
BlockSparseMatrix<Scalar_, Options_, BlockRows_, BlockCols_, StorageIndex_>::transposeImpl() const {
  using ResultType = BlockSparseMatrix<Scalar_, Options_, BlockCols_, BlockRows_, StorageIndex_>;
  ResultType result(blockCols(), blockRows());

  // Count entries per new outer (= old inner).
  for (Index id = 0; id < nonZeroBlocks(); ++id) result.m_outerIndex(m_innerIndex(id) + 1)++;

  // Prefix sum.
  for (Index j = 0; j < result.m_blockOuterSize; ++j) result.m_outerIndex(j + 1) += result.m_outerIndex(j);

  Index nnz = nonZeroBlocks();
  result.resizeBlockStorage_(nnz);

  // One insertion cursor per new outer; start at the prefix-sum boundary.
  // Because we iterate oldOuter in increasing order, for each newOuter = oldInner
  // the emitted newInner = oldOuter values are automatically sorted.
  Array<StorageIndex_, Dynamic, 1> pos = result.m_outerIndex.head(result.m_blockOuterSize);

  for (Index oldOuter = 0; oldOuter < m_blockOuterSize; ++oldOuter) {
    for (Index id = m_outerIndex(oldOuter); id < m_outerIndex(oldOuter + 1); ++id) {
      Index newOuter = m_innerIndex(id);
      Index insertAt = pos(newOuter)++;
      result.m_innerIndex(insertAt) = StorageIndex_(oldOuter);
      result.blockRef(insertAt) = internal::adjoint_if<Conjugate>(blockRef(id));
    }
  }
  return result;
}

template <typename Scalar_, int Options_, int BlockRows_, int BlockCols_, typename StorageIndex_>
BlockSparseMatrix<Scalar_, Options_, BlockCols_, BlockRows_, StorageIndex_>
BlockSparseMatrix<Scalar_, Options_, BlockRows_, BlockCols_, StorageIndex_>::transpose() const {
  return transposeImpl<false>();
}

template <typename Scalar_, int Options_, int BlockRows_, int BlockCols_, typename StorageIndex_>
BlockSparseMatrix<Scalar_, Options_, BlockCols_, BlockRows_, StorageIndex_>
BlockSparseMatrix<Scalar_, Options_, BlockRows_, BlockCols_, StorageIndex_>::adjoint() const {
  return transposeImpl<true>();
}

namespace internal {

// ---------------------------------------------------------------------------
// generic_product_impl: BlockSparse × Dense → Dense
// Provides evalTo / addTo / subTo / scaleAndAddTo via generic_product_impl_base.
// ---------------------------------------------------------------------------
template <typename Lhs, typename Rhs, int ProductType>
struct generic_product_impl<Lhs, Rhs, BlockSparseShape, DenseShape, ProductType>
    : generic_product_impl_base<Lhs, Rhs, generic_product_impl<Lhs, Rhs, BlockSparseShape, DenseShape, ProductType>> {
  using Scalar = typename Product<Lhs, Rhs>::Scalar;

  template <typename Dst>
  static void scaleAndAddTo(Dst& dst, const Lhs& lhs, const Rhs& rhs, const Scalar& alpha) {
    constexpr bool IsRM = (Lhs::Options & RowMajorBit) != 0;
    constexpr int BR = Lhs::BlockRows;
    constexpr int BC = Lhs::BlockCols;
    const typename Lhs::StorageIndex* outerPtr = lhs.outerIndexPtr();
    const typename Lhs::StorageIndex* innerPtr = lhs.innerIndexPtr();
    // Branch on alpha before the loop: alpha==1 and alpha==-1 avoid creating a
    // CwiseUnaryOp<scalar_multiple, B×B_block>, which defeats SIMD for complex scalars.
    bool a1 = (alpha == Scalar(1));
    bool am1 = (alpha == Scalar(-1));
    for (Eigen::Index out = 0; out < lhs.blockOuterSize(); ++out) {
      for (Eigen::Index id = outerPtr[out]; id < outerPtr[out + 1]; ++id) {
        Eigen::Index inner = innerPtr[id];
        Eigen::Index bi = IsRM ? out : inner;
        Eigen::Index bj = IsRM ? inner : out;
        auto dst_seg = dst.template middleRows<BR>(bi * BR);
        auto rhs_seg = rhs.template middleRows<BC>(bj * BC);
        if (EIGEN_PREDICT_TRUE(a1))
          dst_seg.noalias() += lhs.blockRef(id) * rhs_seg;
        else if (am1)
          dst_seg.noalias() -= lhs.blockRef(id) * rhs_seg;
        else {
          // Materialize block×rhs_seg into a small fixed-size stack buffer, then
          // scale by alpha.  Keeps the B×B block as a plain Map for vectorization.
          typedef Matrix<Scalar, BR, Rhs::ColsAtCompileTime> TmpType;
          TmpType tmp(BR, rhs.cols());
          tmp.noalias() = lhs.blockRef(id) * rhs_seg;
          dst_seg.noalias() += alpha * tmp;
        }
      }
    }
  }
};

// ---------------------------------------------------------------------------
// generic_product_impl: Dense × BlockSparse → Dense
// ---------------------------------------------------------------------------
template <typename Lhs, typename Rhs, int ProductType>
struct generic_product_impl<Lhs, Rhs, DenseShape, BlockSparseShape, ProductType>
    : generic_product_impl_base<Lhs, Rhs, generic_product_impl<Lhs, Rhs, DenseShape, BlockSparseShape, ProductType>> {
  using Scalar = typename Product<Lhs, Rhs>::Scalar;

  template <typename Dst>
  static void scaleAndAddTo(Dst& dst, const Lhs& lhs, const Rhs& rhs, const Scalar& alpha) {
    constexpr bool IsRM = (Rhs::Options & RowMajorBit) != 0;
    constexpr int BR = Rhs::BlockRows;
    constexpr int BC = Rhs::BlockCols;
    const typename Rhs::StorageIndex* outerPtr = rhs.outerIndexPtr();
    const typename Rhs::StorageIndex* innerPtr = rhs.innerIndexPtr();
    bool a1 = (alpha == Scalar(1));
    bool am1 = (alpha == Scalar(-1));
    for (Eigen::Index out = 0; out < rhs.blockOuterSize(); ++out) {
      for (Eigen::Index id = outerPtr[out]; id < outerPtr[out + 1]; ++id) {
        Eigen::Index inner = innerPtr[id];
        Eigen::Index bi = IsRM ? out : inner;
        Eigen::Index bj = IsRM ? inner : out;
        auto dst_seg = dst.template middleCols<BC>(bj * BC);
        auto lhs_seg = lhs.template middleCols<BR>(bi * BR);
        if (EIGEN_PREDICT_TRUE(a1))
          dst_seg.noalias() += lhs_seg * rhs.blockRef(id);
        else if (am1)
          dst_seg.noalias() -= lhs_seg * rhs.blockRef(id);
        else {
          typedef Matrix<Scalar, Lhs::RowsAtCompileTime, BC> TmpType;
          TmpType tmp(lhs.rows(), BC);
          tmp.noalias() = lhs_seg * rhs.blockRef(id);
          dst_seg.noalias() += alpha * tmp;
        }
      }
    }
  }
};

}  // namespace internal

}  // end namespace Eigen

#endif  // EIGEN_BLOCKSPARSEMATRIX_H
