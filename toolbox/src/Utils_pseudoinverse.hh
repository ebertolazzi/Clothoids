/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2025                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// file: Utils_pseudoinverse.hh
//
// Implementation of Tikhonov solvers for regularized least squares problems:
// \f[
//   \min_{x} \|Ax - b\|^2 + \lambda \|x\|^2
// \f]
//
// Includes implementations for both dense and sparse matrices,
// with QR and KKT approaches.
//

#pragma once

#ifndef UTILS_PSEUDOINVERSE_dot_HH
#define UTILS_PSEUDOINVERSE_dot_HH

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wswitch-enum"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#pragma clang diagnostic ignored "-Wweak-vtables"
#pragma clang diagnostic ignored "-Wshorten-64-to-32"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wdouble-promotion"
#pragma clang diagnostic ignored "-Wsigned-enum-bitfield"
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wweak-vtables"
#pragma clang diagnostic ignored "-Wunused-template"
#pragma clang diagnostic ignored "-Wnon-virtual-dtor"
#pragma clang diagnostic ignored "-Wpadded"
#pragma clang diagnostic ignored "-Wmissing-noreturn"
#endif

#include "Utils.hh"
#include "Utils_eigen.hh"
#include "Utils_fmt.hh"

namespace Utils
{

  using std::abs;
  using std::max;
  using std::min;
  using std::pow;
  using std::sqrt;

  //====================================================================
  // TikhonovSolver: Dense matrices
  //====================================================================

  /**
   * @brief Tikhonov solver for dense matrices using augmented QR approach
   *
   * Solves the regularized least squares problem:
   * \f[
   *   \min_{x} \|Ax - b\|^2 + \lambda^2 \|x\|^2
   * \f]
   *
   * Uses a QR factorization on the augmented matrix:
   * \f[
   *   \begin{bmatrix} A \\ \lambda I \end{bmatrix}
   * \f]
   *
   * Complexity: \f$O((m+n)n^2)\f$ for factorization, \f$O(n^2)\f$ for solve.
   *
   * @tparam Scalar Scalar type (double, float)
   */
  template <typename Scalar>
  class TikhonovSolver
  {
  private:
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    Matrix                             m_A;                        ///< Original m×n matrix
    Eigen::Index                       m_m, m_n;                   ///< Dimensions m (rows), n (cols)
    Eigen::ColPivHouseholderQR<Matrix> m_qr;                       ///< Column-pivoting QR factorization
    bool                               m_lambda_gt_zero{ false };  ///< Flag for λ > 0
    Scalar                             m_lambda;                   ///< Regularization parameter √|λ|

    /**
     * @brief Solve \f$(A^T A + \lambda^2 I) x = c\f$ without forming \f$A^T
     * A\f$
     *
     * Uses the QR factorization: if \f$A = QR\f$, then:
     * \f[
     *   (A^T A + \lambda^2 I) = R^T R + \lambda^2 I
     * \f]
     * Solution is obtained by solving two triangular systems.
     *
     * @param c Right-hand side vector
     * @return Solution x
     */
    Vector
    qr_solve_transpose( Vector const & c ) const
    {
      auto         R   = m_qr.matrixR().topRows( m_n );
      auto const & P   = m_qr.colsPermutation();
      Vector       rhs = P.transpose() * c;
      Vector       y   = R.transpose().template triangularView<Eigen::Lower>().solve( rhs );
      Vector       z   = R.template triangularView<Eigen::Upper>().solve( y );
      return P * z;
    }

  public:
    /**
     * @brief Constructor
     *
     * @param A Dense m×n matrix
     * @param lambda Regularization parameter (≥ 0)
     */
    TikhonovSolver( Matrix const & A, Scalar lambda ) : m_A( A ), m_m( A.rows() ), m_n( A.cols() )
    {
      UTILS_ASSERT( lambda >= 0, "TikhonovSolver( λ={} ) λ must be >= 0", lambda );
      m_lambda         = std::sqrt( lambda );
      m_lambda_gt_zero = m_lambda > Scalar( 0 );
      if ( m_lambda_gt_zero )
      {
        // Build augmented matrix [A; λI] of size (m+n)×n
        Matrix A_aug( m_m + m_n, m_n );
        A_aug.topRows( m_m )    = m_A;
        A_aug.bottomRows( m_n ) = m_lambda * Matrix::Identity( m_n, m_n );
        m_qr.compute( A_aug );
      }
      else
      {
        m_qr.compute( m_A );  // No regularization, use standard QR
      }
    }

    /**
     * @brief Solve the Tikhonov system
     *
     * @param b Right-hand side vector of size m
     * @return Solution x minimizing \f$\|Ax - b\|^2 + \lambda^2 \|x\|^2\f$
     */
    Vector
    solve( Vector const & b ) const
    {
      if ( m_lambda_gt_zero )
      {
        Vector b_aug( m_m + m_n );
        b_aug.head( m_m ) = b;
        b_aug.tail( m_n ).setZero();
        return m_qr.solve( b_aug );
      }
      else
      {
        return m_qr.solve( b );
      }
    }

    /**
     * @brief Compute \f$y = A (A^T A + \lambda^2 I)^{-1} c\f$
     *
     * Useful for dual problems or for computing matrix-vector products
     * with the regularized pseudoinverse.
     *
     * @param c Input vector of size n
     * @return \f$y = A (A^T A + \lambda^2 I)^{-1} c\f$ of size m
     */
    Vector
    solve_transpose( Vector const & c ) const
    {
      return m_A * qr_solve_transpose( c );
    }

    /// @brief Return the regularization parameter λ
    Scalar
    lambda() const
    {
      return m_lambda;
    }

    /// @brief Return number of rows of matrix A
    Eigen::Index
    rows() const
    {
      return m_m;
    }

    /// @brief Return number of columns of matrix A
    Eigen::Index
    cols() const
    {
      return m_n;
    }
  };

  //====================================================================
  // TikhonovSolver2: Dense matrices with KKT formulation
  //====================================================================

  /**
   * @brief Tikhonov solver for dense matrices using KKT formulation
   *
   * Solves the regularized least squares problem:
   * \f[
   *   \min_{x} \|Ax - b\|^2 + \lambda \|x\|^2
   * \f]
   *
   * Via the Karush-Kuhn-Tucker (KKT) system:
   * \f[
   *   \begin{bmatrix}
   *     I_m & A \\
   *     A^T & -\lambda I_n
   *   \end{bmatrix}
   *   \begin{bmatrix}
   *     y \\
   *     x
   *   \end{bmatrix}
   *   =
   *   \begin{bmatrix}
   *     b \\
   *     0
   *   \end{bmatrix}
   * \f]
   *
   * where y is an auxiliary variable. After factorizing the KKT matrix,
   * the linear system is solved and the solution x is extracted.
   *
   * Complexity: \f$O((m+n)^3)\f$ for LU factorization, \f$O((m+n)^2)\f$ for
   * solve.
   *
   * @tparam Scalar Scalar type (double, float)
   */
  template <typename Scalar>
  class TikhonovSolver2
  {
  private:
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    Matrix                      m_A;                  ///< Original m×n matrix
    Eigen::Index                m_m, m_n;             ///< Dimensions m (rows), n (cols)
    Matrix                      m_KKT;                ///< KKT matrix of size (m+n)×(m+n)
    Eigen::PartialPivLU<Matrix> m_lu;                 ///< Partial pivoting LU factorization
    Scalar                      m_lambda;             ///< Regularization parameter (≥ 0)
    bool                        m_reg;                ///< Flag for λ > 0
    bool                        m_use_diag{ false };  ///< Use diagonal matrix D instead of λI
    Vector                      m_D;                  ///< Diagonal regularization matrix (if m_use_diag)

    /**
     * @brief Build and factorize the KKT matrix
     *
     * Block structure:
     * \f[
     *   K = \begin{bmatrix}
     *         I_m   & A     \\
     *         A^T   & -\lambda I_n
     *       \end{bmatrix}
     * \f]
     *
     * For λ = 0, the bottom-right block is zero.
     */
    void
    build_and_factorize()
    {
      const Eigen::Index M = m_m;
      const Eigen::Index N = m_n;
      const Eigen::Index K = M + N;  // Size of KKT matrix

      m_KKT.resize( K, K );

      // Top-left block: I_m
      m_KKT.topLeftCorner( M, M ).setIdentity();

      // Top-right block: A (rows 0..M-1, cols M..M+N-1)
      m_KKT.topRightCorner( M, N ) = m_A;

      // Bottom-left block: A^T (rows M..M+N-1, cols 0..M-1)
      m_KKT.bottomLeftCorner( N, M ) = m_A.transpose();

      // Bottom-right block: -D or -λ I_n
      if ( m_use_diag )
      {
        // Set diagonal elements individually
        for ( Eigen::Index j = 0; j < N; ++j ) { m_KKT( M + j, M + j ) = -m_D( j ); }
      }
      else if ( m_reg && m_lambda != Scalar( 0 ) )
      {
        m_KKT.bottomRightCorner( N, N ) = -m_lambda * Matrix::Identity( N, N );
      }
      else
      {
        m_KKT.bottomRightCorner( N, N ).setZero();
      }

      // Factorize with LU (PartialPivLU is efficient for dense matrices)
      m_lu.compute( m_KKT );

      UTILS_ASSERT( m_lu.info() == Eigen::Success, "TikhonovSolver2: KKT factorization failed" );
    }

  public:
    /**
     * @brief Construct dense Tikhonov solver via KKT factorization
     *
     * @param A Dense m×n matrix
     * @param lambda Regularization parameter (≥ 0)
     */
    TikhonovSolver2( Matrix const & A, Scalar lambda )
      : m_A( A ), m_m( A.rows() ), m_n( A.cols() ), m_reg( lambda > Scalar( 0 ) ), m_use_diag( false )
    {
      UTILS_ASSERT( lambda >= 0, "TikhonovSolver2( λ={} ) λ must be >= 0", lambda );
      m_lambda = lambda;
      build_and_factorize();
    }

    /**
     * @brief Construct dense Tikhonov solver via KKT factorization with diagonal regularization
     *
     * @param A Dense m×n matrix
     * @param D Diagonal regularization vector (≥ 0 element-wise)
     */
    TikhonovSolver2( Matrix const & A, Vector const & D )
      : m_A( A )
      , m_m( A.rows() )
      , m_n( A.cols() )
      , m_reg( D.size() > 0 && D.minCoeff() > 0 )
      , m_use_diag( true )
      , m_D( D )
    {
      UTILS_ASSERT( D.size() == m_n, "TikhonovSolver2: D size {} must match A cols {}", D.size(), m_n );
      UTILS_ASSERT( D.minCoeff() >= 0, "TikhonovSolver2: D elements must be >= 0" );
      m_lambda = 0;
      build_and_factorize();
    }

    /**
     * @brief Solve regularized least squares
     *
     * Solves \f$ \min_{x} \|Ax - b\|^2 + \lambda \|x\|^2 \f$
     *
     * @param b Right-hand side vector of length m
     * @return Solution x of length n
     */
    Vector
    solve( Vector const & b ) const
    {
      UTILS_ASSERT( b.size() == m_m, "TikhonovSolver2::solve: b size mismatch" );

      // RHS for KKT: [b; 0]
      Vector rhsKKT      = Vector::Zero( m_m + m_n );
      rhsKKT.head( m_m ) = b;

      Vector sol = m_lu.solve( rhsKKT );
      UTILS_ASSERT( m_lu.info() == Eigen::Success, "TikhonovSolver2::solve: solve failed" );

      // Extract x (tail of the solution)
      return sol.segment( m_m, m_n );
    }

    /**
     * @brief Compute \f$ y = A (A^T A + \lambda I)^{-1} c \f$
     *
     * Solves KKT * [y; x] = [0; -c], which yields
     * \f$ x = (A^T A + \lambda I)^{-1} c \f$ and returns A*x.
     *
     * @param c Input vector of length n
     * @return \f$ y = A (A^T A + \lambda I)^{-1} c \f$ of length m
     */
    Vector
    solve_transpose( Vector const & c ) const
    {
      UTILS_ASSERT( c.size() == m_n, "TikhonovSolver2::solve_transpose: c size mismatch" );

      Vector rhsKKT = Vector::Zero( m_m + m_n );
      // RHS = [0; -c] so that the resulting x satisfies (A^T A + λI) x = c
      rhsKKT.segment( m_m, m_n ) = -c;

      Vector sol = m_lu.solve( rhsKKT );
      UTILS_ASSERT( m_lu.info() == Eigen::Success, "TikhonovSolver2::solve_transpose: solve failed" );

      Vector x = sol.segment( m_m, m_n );
      return m_A * x;
    }

    /// @brief Return regularization parameter λ (only valid if constructed with scalar)
    Scalar
    lambda() const
    {
      return m_lambda;
    }

    /// @brief Return diagonal regularization vector (only valid if constructed with diagonal)
    Vector const &
    diagonal() const
    {
      return m_D;
    }

    /// @brief Check if using diagonal regularization
    bool
    uses_diagonal() const
    {
      return m_use_diag;
    }

    /// @brief Return number of rows of matrix A
    Eigen::Index
    rows() const
    {
      return m_m;
    }

    /// @brief Return number of columns of matrix A
    Eigen::Index
    cols() const
    {
      return m_n;
    }

    /**
     * @brief Get the KKT matrix (mainly for debugging/testing)
     * @return Constant reference to the KKT matrix
     */
    Matrix const &
    get_KKT() const
    {
      return m_KKT;
    }

    /**
     * @brief Get the LU factorization object (mainly for debugging/testing)
     * @return Constant reference to the LU factorization
     */
    Eigen::PartialPivLU<Matrix> const &
    get_lu() const
    {
      return m_lu;
    }
  };

  //====================================================================
  // SP_TikhonovSolver: Sparse matrices
  //====================================================================

  /**
   * @brief Tikhonov solver for sparse matrices using augmented QR approach
   *
   * Solves the regularized least squares problem:
   * \f[
   *   \min_{x} \|Ax - b\|^2 + \lambda^2 \|x\|^2
   * \f]
   *
   * Uses a sparse QR factorization on the augmented matrix:
   * \f[
   *   \begin{bmatrix} A \\ \lambda I \end{bmatrix}
   * \f]
   *
   * Complexity depends on sparsity pattern. Suitable for large sparse problems.
   *
   * @tparam Scalar Scalar type (double, float)
   */
  template <typename Scalar>
  class SP_TikhonovSolver
  {
  private:
    using SparseMatrix = Eigen::SparseMatrix<Scalar>;
    using Vector       = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using SparseQR     = Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>>;

    SparseMatrix m_A;                        ///< Original sparse m×n matrix
    Eigen::Index m_m, m_n;                   ///< Dimensions m (rows), n (cols)
    SparseQR     m_qr;                       ///< Sparse QR factorization with column ordering
    Scalar       m_lambda;                   ///< Regularization parameter √|λ|
    bool         m_lambda_gt_zero{ false };  ///< Flag for λ > 0

    /**
     * @brief Solve \f$(A^T A + \lambda^2 I) x = c\f$ using sparse QR factors
     *
     * @param c Right-hand side vector
     * @return Solution x
     */
    Vector
    qr_solve_transpose( Vector const & c ) const
    {
      auto const & R = m_qr.matrixR();
      auto const & P = m_qr.colsPermutation();
      Eigen::Index n = m_n;

      SparseMatrix R_top = R.topRows( n );
      Vector       rhs   = P.transpose() * c;
      Vector       y     = R_top.transpose().template triangularView<Eigen::Lower>().solve( rhs );
      Vector       x     = R_top.template triangularView<Eigen::Upper>().solve( y );
      return P * x;
    }

  public:
    /**
     * @brief Constructor for sparse Tikhonov solver
     *
     * @param A Sparse m×n matrix
     * @param lambda Regularization parameter (≥ 0)
     */
    SP_TikhonovSolver( SparseMatrix const & A, Scalar lambda ) : m_A( A ), m_m( A.rows() ), m_n( A.cols() )
    {
      UTILS_ASSERT( lambda >= 0, "SP_TikhonovSolver( λ={} ) λ must be >= 0", lambda );
      m_lambda         = std::sqrt( lambda );
      m_lambda_gt_zero = lambda > Scalar( 0 );

      if ( m_lambda_gt_zero )
      {
        // Build augmented sparse matrix [A; λI] using triplets
        SparseMatrix                        A_aug( m_m + m_n, m_n );
        std::vector<Eigen::Triplet<Scalar>> triplets;
        triplets.reserve( A.nonZeros() + m_n );

        // Copy original matrix entries
        for ( int k = 0; k < A.outerSize(); ++k )
          for ( typename SparseMatrix::InnerIterator it( A, k ); it; ++it )
            triplets.emplace_back( it.row(), it.col(), it.value() );

        // Add λI entries on the bottom
        for ( Eigen::Index i = 0; i < m_n; ++i ) triplets.emplace_back( m_m + i, i, m_lambda );

        A_aug.setFromTriplets( triplets.begin(), triplets.end() );
        m_qr.compute( A_aug );
      }
      else
      {
        m_qr.compute( m_A );
      }

      UTILS_ASSERT( m_qr.info() == Eigen::Success, "QR decomposition failed" );
    }

    /**
     * @brief Solve the Tikhonov system for sparse matrices
     *
     * @param b Right-hand side vector of size m
     * @return Solution x minimizing \f$\|Ax - b\|^2 + \lambda^2 \|x\|^2\f$
     */
    Vector
    solve( Vector const & b ) const
    {
      if ( m_lambda_gt_zero )
      {
        Vector b_aug      = Vector::Zero( m_m + m_n );
        b_aug.head( m_m ) = b;
        return m_qr.solve( b_aug );
      }
      else
      {
        return m_qr.solve( b );
      }
    }

    /**
     * @brief Compute \f$y = A (A^T A + \lambda^2 I)^{-1} c\f$ for sparse
     * matrices
     *
     * @param c Input vector of size n
     * @return \f$y = A (A^T A + \lambda^2 I)^{-1} c\f$ of size m
     */
    Vector
    solve_transpose( Vector const & c ) const
    {
      return m_A * qr_solve_transpose( c );
    }

    /// @brief Return regularization parameter λ
    Scalar
    lambda() const
    {
      return m_lambda;
    }

    /// @brief Return number of rows of matrix A
    Eigen::Index
    rows() const
    {
      return m_m;
    }

    /// @brief Return number of columns of matrix A
    Eigen::Index
    cols() const
    {
      return m_n;
    }
  };

  //====================================================================
  // SP_TikhonovSolver2: Sparse matrices with KKT formulation
  //====================================================================

  /**
   * @brief Tikhonov solver for sparse matrices using KKT formulation
   *
   * Solves the regularized least squares problem:
   * \f[
   *   \min_{x} \|Ax - b\|^2 + \lambda \|x\|^2
   * \f]
   *
   * Via the Karush-Kuhn-Tucker (KKT) system:
   * \f[
   *   \begin{bmatrix}
   *     I_m & A \\
   *     A^T & -\lambda I_n
   *   \end{bmatrix}
   *   \begin{bmatrix}
   *     y \\
   *     x
   *   \end{bmatrix}
   *   =
   *   \begin{bmatrix}
   *     b \\
   *     0
   *   \end{bmatrix}
   * \f]
   *
   * This implementation uses a hybrid factorization approach:
   * 1. First attempts symmetric LDLT factorization (faster for symmetric
   * matrices)
   * 2. Falls back to general LU factorization if LDLT fails due to numerical
   * issues
   *
   * Complexity depends on sparsity pattern. Suitable for large sparse
   * symmetric indefinite systems.
   *
   * @tparam Scalar Scalar type (double, float)
   */
  template <typename Scalar>
  class SP_TikhonovSolver2
  {
  private:
    using SparseMatrix = Eigen::SparseMatrix<Scalar>;
    using Vector       = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Triplet      = Eigen::Triplet<Scalar>;

    SparseMatrix m_A;       ///< Original sparse m×n matrix
    Eigen::Index m_m, m_n;  ///< Dimensions m (rows), n (cols)
    SparseMatrix m_KKT;     ///< Sparse KKT matrix of size (m+n)×(m+n)

    // Sparse symmetric LDLT factorization with AMD ordering or LU
    Eigen::SimplicialLDLT<SparseMatrix, Eigen::Lower, Eigen::AMDOrdering<int>> m_ldlt;
    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>>                  m_lu;

    bool   m_use_lflt;           ///< use symmetric LDLT factorization or back to LU
    Scalar m_lambda;             ///< Regularization parameter (≥ 0)
    bool   m_reg;                ///< Flag for λ > 0
    bool   m_use_diag{ false };  ///< Use diagonal matrix D instead of λI
    Vector m_D;                  ///< Diagonal regularization matrix (if m_use_diag)

    /**
     * @brief Build and factorize the sparse KKT matrix
     *
     * Constructs the symmetric indefinite KKT matrix using triplets
     * and performs sparse LDLT factorization with AMD ordering.
     */
    void
    build_and_factorize()
    {
      const Eigen::Index M = m_m;
      const Eigen::Index N = m_n;
      const Eigen::Index K = M + N;

      std::vector<Triplet> trips;
      // Reserve space: 2*A.nonZeros() + max(M,N) for diagonal
      trips.reserve( m_A.nonZeros() * 2 + std::max( M, N ) );

      // I_m block (top-left)
      for ( Eigen::Index i = 0; i < M; ++i ) trips.emplace_back( i, i, Scalar( 1 ) );

      // A block (top-right)
      // A^T block (bottom-left) - maintaining symmetry
      for ( int k = 0; k < m_A.outerSize(); ++k )
      {
        for ( typename SparseMatrix::InnerIterator it( m_A, k ); it; ++it )
        {
          trips.emplace_back( it.row(), M + it.col(), it.value() );
          trips.emplace_back( M + it.col(), it.row(), it.value() );
        }
      }

      // -D or -λ I_n block (bottom-right)
      if ( m_use_diag )
      {
        for ( Eigen::Index j = 0; j < N; ++j ) trips.emplace_back( M + j, M + j, -m_D( j ) );
      }
      else if ( m_reg )
      {
        for ( Eigen::Index j = 0; j < N; ++j ) trips.emplace_back( M + j, M + j, -m_lambda );
      }

      m_KKT.resize( K, K );
      m_KKT.setFromTriplets( trips.begin(), trips.end() );
      m_KKT.makeCompressed();

      // Try LDLT factorization first (faster for symmetric matrices)
      m_ldlt.analyzePattern( m_KKT );
      m_use_lflt = m_ldlt.info() == Eigen::Success;
      if ( m_use_lflt )
      {
        m_ldlt.factorize( m_KKT );
        m_use_lflt = m_ldlt.info() == Eigen::Success;
      }

      // If LDLT failed, try LU factorization
      if ( !m_use_lflt )
      {
        m_lu.analyzePattern( m_KKT );
        m_lu.factorize( m_KKT );
        UTILS_ASSERT( m_lu.info() == Eigen::Success, "SP_TikhonovSolver2: LU factorization failed" );
      }
    }

  public:
    /**
     * @brief Constructor for sparse KKT Tikhonov solver
     *
     * @param A Sparse m×n matrix
     * @param lambda Regularization parameter (≥ 0)
     */
    SP_TikhonovSolver2( SparseMatrix const & A, Scalar lambda ) : m_A( A ), m_m( A.rows() ), m_n( A.cols() )
    {
      UTILS_ASSERT( lambda >= 0, "SP_TikhonovSolver2( λ={} ) λ must be >= 0", lambda );
      m_lambda   = lambda;
      m_reg      = lambda > Scalar( 0 );
      m_use_diag = false;
      build_and_factorize();
    }

    /**
     * @brief Constructor for sparse KKT Tikhonov solver with diagonal regularization
     *
     * @param A Sparse m×n matrix
     * @param D Diagonal regularization vector (≥ 0 element-wise)
     */
    SP_TikhonovSolver2( SparseMatrix const & A, Vector const & D ) : m_A( A ), m_m( A.rows() ), m_n( A.cols() )
    {
      UTILS_ASSERT( D.size() == m_n, "SP_TikhonovSolver2: D size {} must match A cols {}", D.size(), m_n );
      UTILS_ASSERT( D.minCoeff() >= 0, "SP_TikhonovSolver2: D elements must be >= 0" );

      m_lambda   = 0;
      m_reg      = D.array().abs().maxCoeff() > 0;
      m_use_diag = true;
      m_D        = D;
      build_and_factorize();
    }

    /**
     * @brief Solve regularized least squares for sparse matrices
     *
     * @param b Right-hand side vector of length m
     * @return Solution x of length n
     */
    Vector
    solve( Vector const & b ) const
    {
      Vector rhs      = Vector::Zero( m_m + m_n );
      rhs.head( m_m ) = b;

      if ( m_use_lflt )
      {
        Vector sol = m_ldlt.solve( rhs );
        UTILS_ASSERT( m_ldlt.info() == Eigen::Success, "SP_TikhonovSolver2::solve (LDLT): solve failed" );
        return sol.segment( m_m, m_n );
      }
      else
      {
        Vector sol = m_lu.solve( rhs );
        UTILS_ASSERT( m_lu.info() == Eigen::Success, "SP_TikhonovSolver2::solve (LU): solve failed" );
        return sol.segment( m_m, m_n );
      }
    }

    /**
     * @brief Compute \f$ y = A (A^T A + \lambda I)^{-1} c \f$ for sparse
     * matrices
     *
     * @param c Input vector of length n
     * @return \f$ y = A (A^T A + \lambda I)^{-1} c \f$ of length m
     */
    Vector
    solve_transpose( Vector const & c ) const
    {
      Vector rhs              = Vector::Zero( m_m + m_n );
      rhs.segment( m_m, m_n ) = -c;

      if ( m_use_lflt )
      {
        Vector sol = m_ldlt.solve( rhs );
        UTILS_ASSERT( m_ldlt.info() == Eigen::Success, "SP_TikhonovSolver2::solve_transpose (LDLT): solve failed" );
        return m_A * sol.segment( m_m, m_n );
      }
      else
      {
        Vector sol = m_lu.solve( rhs );
        UTILS_ASSERT( m_lu.info() == Eigen::Success, "SP_TikhonovSolver2::solve_transpose (LU): solve failed" );
        return m_A * sol.segment( m_m, m_n );
      }
    }

    /**
     * @brief Get the KKT matrix (mainly for debugging/testing)
     * @return Constant reference to the KKT matrix
     */
    SparseMatrix const &
    get_KKT() const
    {
      return m_KKT;
    }

    /// @brief Return regularization parameter λ (only valid if constructed with scalar)
    Scalar
    lambda() const
    {
      return m_lambda;
    }

    /// @brief Return diagonal regularization vector (only valid if constructed with diagonal)
    Vector const &
    diagonal() const
    {
      return m_D;
    }

    /// @brief Check if using diagonal regularization
    bool
    uses_diagonal() const
    {
      return m_use_diag;
    }

    /// @brief Return number of rows of matrix A
    Eigen::Index
    rows() const
    {
      return m_m;
    }

    /// @brief Return number of columns of matrix A
    Eigen::Index
    cols() const
    {
      return m_n;
    }
  };

}  // namespace Utils

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif

//
// eof: Utils_pseudoinverse.hh
//
