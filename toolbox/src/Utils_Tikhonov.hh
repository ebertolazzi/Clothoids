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
// file: Utils_Tikhonov.hh
//
// Implementation of Tikhonov solvers for regularized least squares problems:
// \f[
//   \min_{x} \|Ax - b\|^2 + \lambda \|D(x-c)\|^2
// \f]
//
// Includes implementations for both dense and sparse matrices,
// with QR and KKT approaches.
//

#pragma once

#ifndef UTILS_TIKHONOV_dot_HH
#define UTILS_TIKHONOV_dot_HH

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
   * @class TikhonovSolver
   * @brief Tikhonov/ridge regression solver for dense matrices using augmented QR approach
   *
   * Solves the regularized least squares (Tikhonov regularization) problem:
   * \f[
   *   \min_{x \in \mathbb{R}^n} \|Ax - b\|^2 + \lambda \|D(x-c)\|^2
   * \f]
   * where:
   * - \f$ \bm{A} \in \mathbb{R}^{m \times n} \f$ is the design matrix
   * - \f$ \bm{b} \in \mathbb{R}^m \f$ is the observation vector
   * - \f$ \bm{D} \in \mathbb{R}^{n \times n} \f$ is a diagonal regularization matrix (defaults to identity)
   * - \f$ \bm{c} \in \mathbb{R}^n \f$ is the prior/center vector (defaults to zero)
   * - \f$ \lambda \ge 0 \f$ is the regularization parameter
   *
   * Uses an augmented QR factorization approach by constructing:
   * \f[
   *   \begin{bmatrix} \bm{A} \\ \sqrt{\lambda} \bm{D} \end{bmatrix} \bm{x} \approx
   *   \begin{bmatrix} \bm{b} \\ \sqrt{\lambda} \bm{D} \bm{c} \end{bmatrix}
   * \f]
   * and solving via column-pivoting QR decomposition.
   *
   * ## Mathematical Derivation
   *
   * The objective function expands as:
   * \f[
   *   J(\bm{x}) = \|\bm{A}\bm{x} - \bm{b}\|^2 + \lambda \|\bm{D}(\bm{x}-\bm{c})\|^2
   *             = (\bm{A}\bm{x} - \bm{b})^\top (\bm{A}\bm{x} - \bm{b}) +
   *               \lambda (\bm{x}-\bm{c})^\top \bm{D}^\top \bm{D} (\bm{x}-\bm{c})
   * \f]
   * Setting the gradient \f$ \nabla_x J(\bm{x}) = \bm{0} \f$ gives the normal equations:
   * \f[
   *   (\bm{A}^\top \bm{A} + \lambda \bm{D}^\top \bm{D}) \bm{x}
   *   = \bm{A}^\top \bm{b} + \lambda \bm{D}^\top \bm{D} \bm{c}
   * \f]
   * This implementation solves the equivalent augmented least squares problem,
   * which is numerically more stable than directly solving the normal equations.
   *
   * ## Special Cases
   *
   * 1. **Ordinary Least Squares (OLS)**: \f$ \lambda = 0 \f$ or \f$ \bm{D} = 0 \f$
   *    - Reduces to \f$ \min_x \|\bm{A}\bm{x} - \bm{b}\|^2 \f$
   *    - Solved via QR decomposition of A
   *
   * 2. **Standard Ridge Regression**: \f$ \bm{D} = \bm{I} \f$, \f$ \bm{c} = 0 \f$
   *    - \f$ \min_x \|\bm{A}\bm{x} - \bm{b}\|^2 + \lambda \|\bm{x}\|^2 \f$
   *    - Augmented system: \f$ [\bm{A}; \sqrt{\lambda} \bm{I}] \bm{x} \approx [\bm{b}; \bm{0}] \f$
   *
   * 3. **Generalized Tikhonov**: \f$ \bm{D} \neq \bm{I} \f$, \f$ \bm{c} \neq 0 \f$
   *    - Allows different regularization strengths per variable
   *    - Centers regularization around prior estimate \bm{c}
   *
   * ## Numerical Properties
   *
   * - Uses column-pivoting QR (`Eigen::ColPivHouseholderQR`) for rank deficiency handling
   * - Regularization parameter threshold: \f$ \lambda < 10\epsilon \f$ treated as zero
   * - Computes \f$ \sqrt{\lambda} \f$ once and caches it
   * - Time complexity: \f$ O((m+n)n^2) \f$ factorization, \f$ O(n^2) \f$ solve
   * - Memory: \f$ O((m+n)n) \f$
   *
   * ## References
   *
   * 1. Golub, G. H., & Van Loan, C. F. (2013). *Matrix Computations* (4th ed.).
   *    Johns Hopkins University Press. §5.3.2 (Tikhonov regularization)
   * 2. Hansen, P. C. (1998). *Rank-Deficient and Discrete Ill-Posed Problems*.
   *    SIAM. Chapter 4 (Regularization by truncated SVD/Tikhonov)
   * 3. Björck, Å. (1996). *Numerical Methods for Least Squares Problems*. SIAM.
   *    §5.1 (Regularized least squares)
   *
   * @tparam Scalar Floating-point type (double, float)
   *
   * @see TikhonovSolver2 for KKT formulation
   * @see SP_TikhonovSolver for sparse version
   */
  template <typename Scalar> class TikhonovSolver
  {
  private:
    using integer = Eigen::Index;
    using Vector  = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Matrix  = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using QR      = Eigen::ColPivHouseholderQR<Matrix>;

    Matrix  m_A;                 ///< Original m×n matrix \f$ \bm{A} \f$
    integer m_m, m_n;            ///< Dimensions: m (rows), n (cols)
    QR      m_QR;                ///< Column-pivoting QR factorization object
    Scalar  m_lambda;            ///< Regularization parameter \f$ \lambda \f$
    Scalar  m_sqrt_lambda;       ///< Cached \f$ \sqrt{\lambda} \f$ for efficiency
    Vector  m_D;                 ///< Diagonal of \f$ \bm{D} \f$ matrix (empty → identity)
    Vector  m_c;                 ///< Center vector \f$ \bm{c} \f$ (empty → zero)
    bool    m_has_reg  = false;  ///< Flag: \f$ \lambda > \epsilon \f$ or \f$ \bm{D} \neq \bm{0} \f$
    bool    m_use_diag = false;  ///< Flag: using non-identity diagonal \f$ \bm{D} \f$
    bool    m_has_c    = false;  ///< Flag: using non-zero center vector \f$ \bm{c} \f$

    /// Numerical threshold for treating λ as zero: \f$ \epsilon = 10 \times \text{machine epsilon} \f$
    static constexpr Scalar epsilon() { return std::numeric_limits<Scalar>::epsilon() * Scalar( 10 ); }

    /**
     * @brief Build the augmented matrix for QR factorization
     *
     * Constructs the extended system matrix:
     * \f[
     *   \bm{A}_{\text{aug}} = \begin{bmatrix} \bm{A} \\ \sqrt{\lambda} \bm{D} \end{bmatrix}
     *   \in \mathbb{R}^{(m+n) \times n}
     * \f]
     *
     * ## Implementation Details
     *
     * Two cases:
     * 1. **Diagonal D**: \f$ \sqrt{\lambda} \bm{D} \f$ added as \f$ n \times n \f$ diagonal block
     * 2. **Identity**: \f$ \sqrt{\lambda} \bm{I} \f$ added as scaled identity
     *
     * If \f$ \lambda \approx 0 \f$ and \f$ \bm{D} \f$ empty, uses only \f$ \bm{A} \f$.
     *
     * @note Uses cached \f$ \sqrt{\lambda} \f$ to avoid repeated computation
     * @warning Avoids computing \f$ \sqrt{0} \f$ by checking threshold first
     */
    void build_augmented_matrix()
    {
      if ( m_has_reg )
      {
        Matrix A_aug( m_m + m_n, m_n );
        A_aug.topRows( m_m ) = m_A;  // Copy A block

        if ( m_use_diag )
        {
          // Build diagonal matrix: \f$ \sqrt{\lambda} \cdot \text{diag}(\bm{D}) \f$
          A_aug.bottomRows( m_n ) = m_sqrt_lambda * m_D.asDiagonal();
        }
        else
        {
          // Use scaled identity: \f$ \sqrt{\lambda} \bm{I} \f$
          A_aug.bottomRows( m_n ) = m_sqrt_lambda * Matrix::Identity( m_n, m_n );
        }

        m_QR.compute( A_aug );  // QR factorization with column pivoting
      }
      else
      {
        // No regularization → ordinary least squares
        m_QR.compute( m_A );
      }
    }

  public:
    /**
     * @brief Construct a Tikhonov solver for dense matrices
     *
     * Initializes solver for problem:
     * \f[
     *   \min_x \|\bm{A}\bm{x} - \bm{b}\|^2 + \lambda \|\bm{D}(\bm{x}-\bm{c})\|^2
     * \f]
     *
     * @param A Dense matrix \f$ \bm{A} \in \mathbb{R}^{m \times n} \f$
     * @param lambda Regularization parameter \f$ \lambda \ge 0 \f$ (default: 0)
     * @param D Diagonal entries of \f$ \bm{D} \f$ matrix, \f$ D_i \ge 0 \f$ (default: empty → \f$ \bm{D} = \bm{I} \f$)
     * @param c Center vector \f$ \bm{c} \in \mathbb{R}^n \f$ (default: empty → \f$ \bm{c} = \bm{0} \f$)
     *
     * @throws UTILS_ASSERT if \f$ \lambda < 0 \f$
     * @throws UTILS_ASSERT if \f$ \text{size}(\bm{D}) \neq n \f$
     * @throws UTILS_ASSERT if any \f$ D_i < 0 \f$
     * @throws UTILS_ASSERT if \f$ \text{size}(\bm{c}) \neq n \f$
     *
     * @note If \f$ \bm{D} \f$ is empty, uses identity matrix \f$ \bm{I} \f$
     * @note If \f$ \bm{c} \f$ is empty, uses zero vector
     * @note Factorization is performed during construction
     */
    explicit TikhonovSolver(
      Matrix const & A,
      Scalar         lambda = Scalar( 0 ),
      Vector const & D      = Vector(),
      Vector const & c      = Vector() )
      : m_A( A ), m_m( A.rows() ), m_n( A.cols() ), m_lambda( lambda ), m_D( D ), m_c( c )
    {
      UTILS_ASSERT( lambda >= 0, "TikhonovSolver( λ={} ) λ must be >= 0", lambda );

      // Check if lambda is numerically significant
      m_has_reg = lambda > epsilon();

      // Cache sqrt(lambda) to avoid redundant computation
      m_sqrt_lambda = m_has_reg ? std::sqrt( lambda ) : Scalar( 0 );

      // Process diagonal regularization matrix D
      if ( D.size() > 0 )
      {
        UTILS_ASSERT( D.size() == m_n, "TikhonovSolver: D size {} must match A cols {}", D.size(), m_n );
        UTILS_ASSERT( D.minCoeff() >= 0, "TikhonovSolver: D elements must be >= 0" );
        m_use_diag = true;
      }

      // Process center vector c
      if ( c.size() > 0 )
      {
        UTILS_ASSERT( c.size() == m_n, "TikhonovSolver: c size {} must match A cols {}", c.size(), m_n );
        m_has_c = true;
      }

      build_augmented_matrix();
    }

    /**
     * @brief Convenience constructor with diagonal regularization D and λ = 1
     *
     * Equivalent to `TikhonovSolver(A, 1.0, D)`
     *
     * @param A Design matrix \f$ \bm{A} \f$
     * @param D Diagonal regularization matrix entries
     */
    explicit TikhonovSolver( Matrix const & A, Vector const & D ) : TikhonovSolver( A, Scalar( 1 ), D ) {}

    /**
     * @brief Solve the Tikhonov-regularized least squares problem
     *
     * Computes the solution:
     * \f[
     *   \bm{x}^* = \arg\min_{\bm{x}} \|\bm{A}\bm{x} - \bm{b}\|^2 + \lambda \|\bm{D}(\bm{x}-\bm{c})\|^2
     * \f]
     *
     * ## Mathematical Steps
     *
     * 1. Construct augmented right-hand side:
     *    \f[
     *      \bm{b}_{\text{aug}} = \begin{bmatrix} \bm{b} \\ \sqrt{\lambda} \bm{D} \bm{c} \end{bmatrix}
     *    \f]
     *    If \f$ \bm{c} = \bm{0} \f$, bottom block is zero.
     *
     * 2. Solve augmented least squares via pre-computed QR:
     *    \f[
     *      \min_x \|\bm{A}_{\text{aug}} \bm{x} - \bm{b}_{\text{aug}}\|^2
     *    \f]
     *
     * 3. Return solution vector \f$ \bm{x}^* \f$.
     *
     * @param b Observation vector \f$ \bm{b} \in \mathbb{R}^m \f$
     * @return Solution vector \f$ \bm{x}^* \in \mathbb{R}^n \f$
     *
     * @throws UTILS_ASSERT if \f$ \text{size}(\bm{b}) \neq m \f$
     */
    Vector solve( Vector const & b ) const
    {
      UTILS_ASSERT( b.size() == m_m, "TikhonovSolver::solve: b size {} must match A rows {}", b.size(), m_m );

      if ( m_has_reg )
      {
        Vector b_aug( m_m + m_n );
        b_aug.head( m_m ) = b;  // Top block: original RHS b

        if ( m_has_c )
        {
          // Bottom block: \f$ \sqrt{\lambda} D c \f$ (element-wise product)
          if ( m_use_diag )
            b_aug.tail( m_n ) = m_sqrt_lambda * m_D.cwiseProduct( m_c );
          else
            b_aug.tail( m_n ) = m_sqrt_lambda * m_c;
        }
        else
        {
          b_aug.tail( m_n ).setZero();  // Zero prior → zero regularization target
        }

        return m_QR.solve( b_aug );  // Solve via QR decomposition
      }
      else
      {
        // No regularization → standard least squares
        return m_QR.solve( b );
      }
    }

    /// @brief Get regularization parameter \f$ \lambda \f$
    Scalar lambda() const { return m_lambda; }

    /// @brief Get diagonal entries of regularization matrix \f$ D \f$
    Vector const & diagonal() const { return m_D; }

    /// @brief Get center vector \f$ c \f$
    Vector const & center() const { return m_c; }

    /// @brief Check if using non-identity diagonal regularization
    bool uses_diagonal() const { return m_use_diag; }

    /// @brief Check if using non-zero center vector
    bool has_center() const { return m_has_c; }

    /// @brief Check if any regularization is active (\f$ \lambda > 0 \f$ or \f$ D \neq 0 \f$)
    bool has_regularization() const { return m_has_reg; }

    /// @brief Get number of rows \f$ m \f$ of matrix \f$ A \f$
    integer rows() const { return m_m; }

    /// @brief Get number of columns \f$ n \f$ of matrix \f$ A \f$
    integer cols() const { return m_n; }

    /// @brief Get numerical rank of factorized matrix (for diagnostics)
    integer rank() const { return m_QR.rank(); }
  };

  //====================================================================
  // TikhonovSolver2: Dense matrices with KKT formulation
  //====================================================================

  /**
   * @class TikhonovSolver2
   * @brief Tikhonov solver for dense matrices using KKT (Karush-Kuhn-Tucker) formulation
   *
   * Solves the same regularized least squares problem as TikhonovSolver:
   * \f[
   *   \min_{\bm{x}} \|\bm{A}\bm{x} - \bm{b}\|^2 + \lambda \|\bm{D}(\bm{x}-\bm{c})\|^2
   * \f]
   * but uses the KKT (saddle-point) formulation from optimization theory.
   *
   * ## Mathematical Formulation
   *
   * Introducing Lagrange multipliers \f$ \bm{y} \in \mathbb{R}^m \f$ for the equality
   * constraint \f$ \bm{y} = \bm{A}\bm{x} - \bm{b} \f$, we form the Lagrangian:
   * \f[
   *   \mathcal{L}(\bm{x},\bm{y}) = \|\bm{y}\|^2 + \lambda \|\bm{D}(\bm{x}-\bm{c})\|^2 + 2\bm{y}^\top(\bm{A}\bm{x} -
   * \bm{b})
   * \f]
   *
   * The KKT optimality conditions \f$ \nabla_{\bm{x},\bm{y}} \mathcal{L} = 0 \f$ yield:
   * \f[
   *   \begin{cases}
   *     \bm{y} + \bm{A} \bm{x} = \bm{b} \\
   *     \bm{A}^\top \bm{y} - \lambda \bm{D}^\top \bm{D} (\bm{x} - \bm{c}) = 0
   *   \end{cases}
   * \f]
   *
   * In matrix form:
   * \f[
   *   \underbrace{\begin{bmatrix}
   *     \bm{I}      & \bm{A} \\
   *     \bm{A}^\top & -\lambda \bm{D}^\top \bm{D}
   *   \end{bmatrix}}_{\text{KKT matrix}}
   *   \begin{bmatrix} \bm{y} \\ \bm{x} \end{bmatrix}
   *   =
   *   \begin{bmatrix} \bm{b} \\ -\lambda \bm{D}^\top \bm{D} \bm{c} \end{bmatrix}
   * \f]
   *
   * The solution \f$ \bm{x}^* \f$ is extracted from the combined vector \f$ [\bm{y}; \bm{x}] \f$.
   *
   * ## Equivalence to Augmented QR
   *
   * This formulation is mathematically equivalent to TikhonovSolver but:
   * - **Advantages**: Better for certain problem structures, allows access to dual variable \f$ \bm{y} \f$
   * - **Disadvantages**: Larger system \f$ (m+n) \times (m+n) \f$ vs \f$ (m+n) \times n \f$
   *
   * The KKT matrix is symmetric indefinite (has both positive and negative eigenvalues).
   *
   * ## Numerical Properties
   *
   * - Uses partial pivoting LU (`Eigen::PartialPivLU`) for stability
   * - Time complexity: \f$ O((m+n)^3) \f$ factorization, \f$ O((m+n)^2) \f$ solve
   * - Memory: \f$ O((m+n)^2) \f$
   * - For \f$ m \gg n \f$, less efficient than augmented QR approach
   *
   * ## References
   *
   * 1. Nocedal, J., & Wright, S. J. (2006). *Numerical Optimization* (2nd ed.).
   *    Springer. §16.2 (KKT conditions)
   * 2. Boyd, S., & Vandenberghe, L. (2004). *Convex Optimization*.
   *    Cambridge University Press. §10.1 (Equality constrained minimization)
   *
   * @tparam Scalar Floating-point type (double, float)
   *
   * @see TikhonovSolver for augmented QR approach
   */
  template <typename Scalar> class TikhonovSolver2
  {
  private:
    using integer = Eigen::Index;
    using Matrix  = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector  = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using LU      = Eigen::PartialPivLU<Matrix>;

    Matrix  m_A;       ///< Design matrix \f$ \bm{A} \f$
    integer m_m, m_n;  ///< Dimensions: m rows, n cols
    Matrix  m_KKT;  ///< KKT matrix \f$ \begin{bmatrix} \bm{I} & A \\ \bm{A}^\top & -\lambda \bm{D}^2 \end{bmatrix} \f$
    LU      m_LU;   ///< LU factorization of KKT matrix
    Scalar  m_lambda;             ///< Regularization parameter \f$ \lambda \f$
    Vector  m_D2;                 ///< Diagonal of \f$ \bm{D}^2 \f$ (empty → identity)
    Vector  m_c;                  ///< Center vector \f$ \bm{c} \f$ (empty → zero)
    bool    m_use_diag{ false };  ///< Flag: using non-identity diagonal D
    bool    m_has_reg{ false };   ///< Flag: \f$ \lambda > 0 \f$ or \f$ D \neq 0 \f$
    bool    m_has_c{ false };     ///< Flag: using non-zero center vector c

    /// Numerical threshold: \f$ \epsilon = 10 \times \text{machine epsilon} \f$
    static constexpr Scalar epsilon() { return std::numeric_limits<Scalar>::epsilon() * Scalar( 10 ); }

    /**
     * @brief Construct and factorize the KKT matrix
     *
     * Builds the symmetric indefinite KKT matrix:
     * \f[
     *   \bm{K} = \begin{bmatrix}
     *     \bm{I}      & \bm{A} \\
     *     \bm{A}^\top & -\lambda \bm{D}^2
     *   \end{bmatrix}
     *   \in \mathbb{R}^{(m+n) \times (m+n)}
     * \f]
     *
     * ## Implementation Details
     *
     * 1. **Top-left**: \f$ m \times m \f$ identity matrix \f$ \bm{I} \f$
     * 2. **Top-right**: \f$ m \times n \f$ matrix \f$ \bm{A} \f$
     * 3. **Bottom-left**: \f$ n \times m \f$ transpose \f$ \bm{A}^\top \f$
     * 4. **Bottom-right**: \f$ n \times n \f$ diagonal matrix \f$ -\lambda \bm{D}^2 \f$
     *
     * Since \f$ \bm{D} \f$ is diagonal, \f$ \bm{D}^2 = \text{diag}(D_1^2, \dots, D_n^2) \f$.
     * Uses cached \f$ D_i^2 \f$ values for efficiency.
     *
     * @note The KKT matrix is symmetric but not positive definite
     * @warning For \f$ \lambda = 0 \f$, bottom-right block is zero (singular if \f$ A \f$ rank-deficient)
     */
    void build_and_factorize()
    {
      integer const M = m_m;
      integer const N = m_n;
      integer const K = M + N;  // Size of KKT matrix

      m_KKT.resize( K, K );
      m_KKT.setZero();

      // Top-left block: \f$ \bm{I} \f$ (m×m identity)
      m_KKT.topLeftCorner( M, M ).setIdentity();

      // Top-right block: \f$ \bm{A} \f$ (m×n)
      m_KKT.topRightCorner( M, N ) = m_A;

      // Bottom-left block: \f$ \bm{A}^\top \f$ (n×m)
      m_KKT.bottomLeftCorner( N, M ) = m_A.transpose();

      // Bottom-right block: \f$ -\lambda \bm{D}^2 \f$ (n×n diagonal)
      if ( m_has_reg )
      {
        if ( m_use_diag )
        {
          // \f$ \bm{D}^2 = \text{diag}(D_1^2, \dots, D_n^2) \f$
          m_KKT.bottomRightCorner( N, N ).diagonal() = -m_lambda * m_D2;
        }
        else
        {
          // \f$ \bm{D} = \bm{I} \f$ → \f$ \bm{D} = \bm{I} \f$
          m_KKT.bottomRightCorner( N, N ).diagonal().setConstant( -m_lambda );
        }
      }
      // else: bottom-right remains zero (λ = 0 or D = 0)

      // Factorize with partial pivoting LU
      m_LU.compute( m_KKT );
      UTILS_ASSERT( m_LU.info() == Eigen::Success, "TikhonovSolver2: KKT factorization failed" );
    }

  public:
    /**
     * @brief Construct KKT-based Tikhonov solver for dense matrices
     *
     * @param A Design matrix \f$ \bm{A} \in \mathbb{R}^{m \times n} \f$
     * @param lambda Regularization parameter \f$ \lambda \ge 0 \f$ (default: 0)
     * @param D Diagonal entries of \f$ \bm{D} \f$, \f$ D_i \ge 0 \f$ (default: empty → \f$ \bm{D} = \bm{I} \f$)
     * @param c Center vector \f$ \bm{c} \in \mathbb{R}^n \f$ (default: empty → \f$ \bm{c} = 0 \f$)
     *
     * @throws UTILS_ASSERT if \f$ \lambda < 0 \f$
     * @throws UTILS_ASSERT if size mismatches occur
     */
    explicit TikhonovSolver2(
      Matrix const & A,
      Scalar         lambda = Scalar( 0 ),
      Vector const & D      = Vector(),
      Vector const & c      = Vector() )
      : m_A( A ), m_m( A.rows() ), m_n( A.cols() ), m_lambda( lambda ), m_D2( D ), m_c( c )
    {
      UTILS_ASSERT( lambda >= 0, "TikhonovSolver2( λ={} ) λ must be >= 0", lambda );

      m_has_reg = lambda > epsilon();

      if ( D.size() > 0 )
      {
        UTILS_ASSERT( D.size() == m_n, "TikhonovSolver2: D size {} must match A cols {}", D.size(), m_n );
        UTILS_ASSERT( D.minCoeff() >= 0, "TikhonovSolver2: D elements must be >= 0" );
        m_use_diag = true;

        // Cache squared diagonal entries: \f$ D_i^2 \f$
        m_D2.array() *= m_D2.array();
      }

      if ( c.size() > 0 )
      {
        UTILS_ASSERT( c.size() == m_n, "TikhonovSolver2: c size {} must match A cols {}", c.size(), m_n );
        m_has_c = true;
      }

      build_and_factorize();
    }

    /**
     * @brief Convenience constructor with diagonal regularization and λ = 1
     *
     * Equivalent to `TikhonovSolver2(A, 1.0, D)`
     */
    explicit TikhonovSolver2( Matrix const & A, Vector const & D ) : TikhonovSolver2( A, Scalar( 1 ), D ) {}

    /**
     * @brief Solve the Tikhonov problem via KKT system
     *
     * Solves the KKT linear system:
     * \f[
     *   \begin{bmatrix} I_m & A \\ A^\top & -\lambda \bm{D}^2 \end{bmatrix}
     *   \begin{bmatrix} y \\ x \end{bmatrix}
     *   =
     *   \begin{bmatrix} b \\ \lambda \bm{D}^2 \bm{c} \end{bmatrix}
     * \f]
     *
     * and returns the solution \f$ x^* \f$.
     *
     * @param b Observation vector \f$ b \in \mathbb{R}^m \f$
     * @return Solution \f$ x^* \in \mathbb{R}^n \f$
     *
     * @throws UTILS_ASSERT if size mismatch
     */
    Vector solve( Vector const & b ) const
    {
      UTILS_ASSERT( b.size() == m_m, "TikhonovSolver2::solve: b size {} must match A rows {}", b.size(), m_m );

      // Right-hand side for KKT system
      Vector rhsKKT( m_m + m_n );
      rhsKKT.head( m_m ) = b;  // Top block: b

      if ( m_has_c && m_has_reg )
      {
        if ( m_use_diag )
        {
          // Bottom block: \f$ \lambda \bm{D}^2 c = \lambda (D_i^2 c_i) \f$
          rhsKKT.tail( m_n ) = -m_lambda * m_D2.cwiseProduct( m_c );
        }
        else
        {
          // \f$ D = I \f$ → bottom block: \f$ \lambda c \f$
          rhsKKT.tail( m_n ) = -m_lambda * m_c;
        }
      }
      else
      {
        rhsKKT.tail( m_n ).setZero();  // Zero prior or no regularization
      }

      // Solve KKT system
      Vector sol = m_LU.solve( rhsKKT );
      UTILS_ASSERT( m_LU.info() == Eigen::Success, "TikhonovSolver2::solve: solve failed" );

      // Extract x from [y; x] solution vector
      return sol.segment( m_m, m_n );
    }

    /// @brief Get regularization parameter \f$ \lambda \f$
    Scalar lambda() const { return m_lambda; }

    /// @brief Get diagonal entries of \f$ D \f$
    Vector const & diagonal2() const { return m_D2; }

    /// @brief Get center vector \f$ c \f$
    Vector const & center() const { return m_c; }

    /// @brief Check if using non-identity diagonal D
    bool uses_diagonal() const { return m_use_diag; }

    /// @brief Check if using non-zero center vector
    bool has_center() const { return m_has_c; }

    /// @brief Check if any regularization active
    bool has_regularization() const { return m_has_reg; }

    /// @brief Get number of rows \f$ m \f$ of A
    integer rows() const { return m_m; }

    /// @brief Get number of columns \f$ n \f$ of A
    integer cols() const { return m_n; }

    /**
     * @brief Get the KKT matrix (for debugging/testing)
     * @return Const reference to the KKT matrix
     */
    Matrix const & get_KKT() const { return m_KKT; }

    /**
     * @brief Get the LU factorization object (for debugging)
     * @return Const reference to LU factorization
     */
    Eigen::PartialPivLU<Matrix> const & get_LU() const { return m_LU; }
  };

  //====================================================================
  // SP_TikhonovSolver: Sparse matrices
  //====================================================================

  /**
   * @class SP_TikhonovSolver
   * @brief Tikhonov solver for sparse matrices using augmented QR approach
   *
   * Sparse version of TikhonovSolver. Solves:
   * \f[
   *   \min_{\bm{x}} \|\bm{A}\bm{x} - \bm{b}\|^2 + \lambda \|\bm{D}(\bm{x}-\bm{c})\|^2
   * \f]
   * where \f$ \bm{A} \f$ is sparse.
   *
   * Uses sparse QR decomposition (`Eigen::SparseQR`) with COLAMD column ordering
   * to minimize fill-in during factorization.
   *
   * ## Implementation Details
   *
   * Constructs augmented sparse matrix:
   * \f[
   *   \bm{A}_{\text{aug}} = \begin{bmatrix} \bm{A} \\ \sqrt{\lambda} \bm{D} \end{bmatrix}
   * \f]
   *
   * The regularization terms are added as diagonal entries at the bottom,
   * preserving sparsity structure.
   *
   * ## Special Considerations for Sparse Matrices
   *
   * 1. **Triplet format**: Uses `Eigen::Triplet` for efficient sparse matrix construction
   * 2. **Memory allocation**: Pre-reserves space for original nonzeros + regularization diagonal
   * 3. **Fill-in reduction**: COLAMD ordering minimizes fill-in during QR factorization
   * 4. **Compressed format**: Converts to compressed column storage after construction
   *
   * ## Complexity
   *
   * - Depends on sparsity pattern and fill-in during QR
   * - Generally \f$ O(\text{nnz} \cdot \text{rank}) \f$ for sparse QR
   * - Much faster than dense methods for large sparse problems with low fill-in
   *
   * @tparam Scalar Floating-point type (double, float)
   *
   * @see TikhonovSolver for dense version
   * @see SP_TikhonovSolver2 for sparse KKT approach
   */
  template <typename Scalar> class SP_TikhonovSolver
  {
  private:
    using integer      = Eigen::Index;
    using SparseMatrix = Eigen::SparseMatrix<Scalar>;
    using Vector       = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using SparseQR     = Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>>;

    SparseMatrix m_A;            ///< Sparse design matrix \f$ \bm{A} \f$
    integer      m_m, m_n;       ///< Dimensions
    SparseQR     m_QR;           ///< Sparse QR factorization with COLAMD ordering
    Scalar       m_lambda;       ///< Regularization parameter
    Scalar       m_sqrt_lambda;  ///< Cached \f$ \sqrt{\lambda} \f$
    Vector       m_D;            ///< Diagonal of \f$ \bm{D} \f$ (empty → identity)
    Vector       m_c;            ///< Center vector (empty → zero)
    bool         m_use_diag{ false };
    bool         m_has_reg{ false };
    bool         m_has_c{ false };

    /// Numerical threshold: \f$ \epsilon = 10 \times \text{machine epsilon} \f$
    static constexpr Scalar epsilon() { return std::numeric_limits<Scalar>::epsilon() * Scalar( 10 ); }

    /**
     * @brief Build augmented sparse matrix for QR factorization
     *
     * Constructs sparse augmented matrix using triplet format.
     * Preserves sparsity pattern while adding regularization diagonal.
     *
     * ## Steps:
     * 1. Reserve space for original nonzeros + regularization diagonal entries
     * 2. Copy original matrix entries via sparse iterators
     * 3. Add regularization diagonal terms at bottom rows
     * 4. Convert to compressed format and compute QR
     */
    void build_augmented_matrix()
    {
      std::vector<Eigen::Triplet<Scalar>> triplets;

      // Reserve: original nonzeros + regularization diagonal (if any)
      integer reg_entries = m_has_reg ? m_n : 0;
      triplets.reserve( m_A.nonZeros() + reg_entries );

      // Copy original sparse matrix entries
      for ( int k = 0; k < m_A.outerSize(); ++k )
        for ( typename SparseMatrix::InnerIterator it( m_A, k ); it; ++it )
          triplets.emplace_back( it.row(), it.col(), it.value() );

      // Add regularization diagonal terms
      if ( m_has_reg )
      {
        if ( m_use_diag )
        {
          // Add \f$ \sqrt{\lambda} D_i \f$ on diagonal at rows m+i, columns i
          for ( integer i = 0; i < m_n; ++i ) triplets.emplace_back( m_m + i, i, m_sqrt_lambda * m_D( i ) );
        }
        else
        {
          // Add \f$ \sqrt{\lambda} \f$ on diagonal (scaled identity)
          for ( integer i = 0; i < m_n; ++i ) triplets.emplace_back( m_m + i, i, m_sqrt_lambda );
        }
      }

      integer      aug_rows = m_m + ( m_has_reg ? m_n : 0 );
      SparseMatrix A_aug( aug_rows, m_n );
      A_aug.setFromTriplets( triplets.begin(), triplets.end() );
      A_aug.makeCompressed();  // Convert to compressed column storage

      m_QR.compute( A_aug );
      UTILS_ASSERT( m_QR.info() == Eigen::Success, "SP_TikhonovSolver: QR decomposition failed" );
    }

  public:
    /**
     * @brief Construct sparse Tikhonov solver
     *
     * @param A Sparse matrix \f$ \bm{A} \in \mathbb{R}^{m \times n} \f$
     * @param lambda Regularization parameter \f$ \lambda \ge 0 \f$
     * @param D Diagonal entries of \f$ \bm{D} \f$ (default: empty → identity)
     * @param c Center vector (default: empty → zero)
     *
     * @throws UTILS_ASSERT on invalid parameters
     */
    explicit SP_TikhonovSolver(
      SparseMatrix const & A,
      Scalar               lambda = Scalar( 0 ),
      Vector const &       D      = Vector(),
      Vector const &       c      = Vector() )
      : m_A( A ), m_m( A.rows() ), m_n( A.cols() ), m_lambda( lambda ), m_D( D ), m_c( c )
    {
      UTILS_ASSERT( lambda >= 0, "SP_TikhonovSolver( λ={} ) λ must be >= 0", lambda );

      m_has_reg     = ( lambda > epsilon() );
      m_sqrt_lambda = m_has_reg ? std::sqrt( lambda ) : Scalar( 0 );

      if ( D.size() > 0 )
      {
        UTILS_ASSERT( D.size() == m_n, "SP_TikhonovSolver: D size {} must match A cols {}", D.size(), m_n );
        UTILS_ASSERT( D.minCoeff() >= 0, "SP_TikhonovSolver: D elements must be >= 0" );
        m_use_diag = true;
      }

      if ( c.size() > 0 )
      {
        UTILS_ASSERT( c.size() == m_n, "SP_TikhonovSolver: c size {} must match A cols {}", c.size(), m_n );
        m_has_c = true;
      }

      build_augmented_matrix();
    }

    /// @brief Convenience constructor with D and λ = 1
    explicit SP_TikhonovSolver( SparseMatrix const & A, Vector const & D ) : SP_TikhonovSolver( A, Scalar( 1 ), D ) {}

    /**
     * @brief Solve sparse Tikhonov problem
     *
     * @param b Observation vector \f$ b \in \mathbb{R}^m \f$
     * @return Solution \f$ x^* \in \mathbb{R}^n \f$
     */
    Vector solve( Vector const & b ) const
    {
      UTILS_ASSERT( b.size() == m_m, "SP_TikhonovSolver::solve: b size {} must match A rows {}", b.size(), m_m );

      if ( m_has_reg )
      {
        Vector b_aug( m_m + m_n );
        b_aug.head( m_m ) = b;

        if ( m_has_c )
        {
          if ( m_use_diag )
            b_aug.tail( m_n ) = m_sqrt_lambda * m_D.cwiseProduct( m_c );
          else
            b_aug.tail( m_n ) = m_sqrt_lambda * m_c;
        }
        else
        {
          b_aug.tail( m_n ).setZero();
        }

        return m_QR.solve( b_aug );
      }
      else
      {
        return m_QR.solve( b );
      }
    }

    /// @name Accessors
    /// @{
    Scalar         lambda() const { return m_lambda; }
    Vector const & diagonal() const { return m_D; }
    Vector const & center() const { return m_c; }
    bool           uses_diagonal() const { return m_use_diag; }
    bool           has_center() const { return m_has_c; }
    bool           has_regularization() const { return m_has_reg; }
    integer        rows() const { return m_m; }
    integer        cols() const { return m_n; }
    integer        rank() const { return m_QR.rank(); }
    /// @}
  };

  //====================================================================
  // SP_TikhonovSolver2: Sparse matrices with KKT formulation
  //====================================================================

  /**
   * @class SP_TikhonovSolver2
   * @brief Tikhonov solver for sparse matrices using KKT formulation
   *
   * Sparse version of TikhonovSolver2. Solves:
   * \f[
   *   \min_{\bm{x}} \|\bm{A}\bm{x} - \bm{b}\|^2 + \lambda \|\bm{D}(\bm{x}-\bm{c})\|^2
   * \f]
   * via the sparse KKT system.
   *
   * ## Factorization Strategy
   *
   * Uses adaptive approach for numerical stability:
   * 1. **Primary**: Symmetric indefinite LDLT factorization (`Eigen::SimplicialLDLT`)
   *    - Exploits symmetry of KKT matrix
   *    - Uses AMD ordering for fill-in reduction
   *    - Faster and uses less memory when stable
   *
   * 2. **Fallback**: General sparse LU factorization (`Eigen::SparseLU`)
   *    - Used if LDLT fails numerically
   *    - More robust but slower
   *    - Uses COLAMD ordering
   *
   * ## Matrix Structure
   *
   * The KKT matrix has block structure:
   * \f[
   *   \begin{bmatrix}
   *     \bm{I} & \bm{A} \\
   *     \bm{A}^\top & -\lambda \bm{D}^2
   *   \end{bmatrix}
   * \f]
   *
   * Constructed using triplets while maintaining symmetry.
   *
   * @tparam Scalar Floating-point type (double, float)
   *
   * @see SP_TikhonovSolver for sparse augmented QR approach
   */
  template <typename Scalar> class SP_TikhonovSolver2
  {
  private:
    using integer      = Eigen::Index;
    using SparseMatrix = Eigen::SparseMatrix<Scalar>;
    using Vector       = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Triplet      = Eigen::Triplet<Scalar>;
    using LDLT         = Eigen::SimplicialLDLT<SparseMatrix, Eigen::Lower, Eigen::AMDOrdering<int>>;
    using SparseLU     = Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>>;

    SparseMatrix m_A;       ///< Sparse design matrix
    integer      m_m, m_n;  ///< Dimensions
    SparseMatrix m_KKT;     ///< Sparse KKT matrix

    // Dual factorization strategy
    LDLT     m_LDLT;  ///< Symmetric indefinite factorization (preferred)
    SparseLU m_LU;    ///< General LU factorization (fallback)

    bool   m_use_LDLT{ true };  ///< True if using LDLT, false if using LU
    Scalar m_lambda;            ///< Regularization parameter
    Vector m_D2;                ///< Diagonal of \f$ \bm{D}^2 \f$
    Vector m_c;                 ///< Center vector
    bool   m_use_diag{ false };
    bool   m_has_reg{ false };
    bool   m_has_c{ false };

    /// Numerical threshold
    static constexpr Scalar epsilon() { return std::numeric_limits<Scalar>::epsilon() * Scalar( 10 ); }

    /**
     * @brief Build and factorize sparse KKT matrix
     *
     * Constructs symmetric sparse KKT matrix using triplets,
     * then attempts LDLT factorization with fallback to LU.
     */
    void build_and_factorize()
    {
      integer const M = m_m;
      integer const N = m_n;

      std::vector<Triplet> trips;
      // Reserve: identity + 2*A.nonZeros (symmetry) + regularization diagonal
      trips.reserve( M + 2 * m_A.nonZeros() + ( m_has_reg ? N : 0 ) );

      // I_m block (top-left diagonal)
      for ( integer i = 0; i < M; ++i ) trips.emplace_back( i, i, Scalar( 1 ) );

      // A block (top-right) and A^T block (bottom-left) - maintain symmetry
      for ( int k = 0; k < m_A.outerSize(); ++k )
      {
        for ( typename SparseMatrix::InnerIterator it( m_A, k ); it; ++it )
        {
          trips.emplace_back( it.row(), M + it.col(), it.value() );  // A
          trips.emplace_back( M + it.col(), it.row(), it.value() );  // A^T
        }
      }

      // -λ D^2 block (bottom-right diagonal)
      if ( m_has_reg )
      {
        if ( m_use_diag )
        {
          // \f$ \bm{D}^2 = \text{diag}(D_i^2) \f$
          for ( integer j = 0; j < N; ++j ) trips.emplace_back( M + j, M + j, -m_lambda * m_D2( j ) );
        }
        else
        {
          // \f$ D = I \f$
          for ( integer j = 0; j < N; ++j ) trips.emplace_back( M + j, M + j, -m_lambda );
        }
      }

      m_KKT.resize( M + N, M + N );
      m_KKT.setFromTriplets( trips.begin(), trips.end() );
      m_KKT.makeCompressed();

      // Try LDLT first (exploits symmetry)
      m_LDLT.analyzePattern( m_KKT );
      m_use_LDLT = ( m_LDLT.info() == Eigen::Success );

      if ( m_use_LDLT )
      {
        m_LDLT.factorize( m_KKT );
        m_use_LDLT = ( m_LDLT.info() == Eigen::Success );
      }

      // Fallback to LU if LDLT failed
      if ( !m_use_LDLT )
      {
        m_LU.analyzePattern( m_KKT );
        m_LU.factorize( m_KKT );
        UTILS_ASSERT( m_LU.info() == Eigen::Success, "SP_TikhonovSolver2: Both LDLT and LU factorizations failed" );
      }
    }

  public:
    /**
     * @brief Construct sparse KKT Tikhonov solver
     *
     * @param A Sparse design matrix
     * @param lambda Regularization parameter
     * @param D Diagonal of regularization matrix
     * @param c Center vector
     */
    explicit SP_TikhonovSolver2(
      SparseMatrix const & A,
      Scalar               lambda = Scalar( 0 ),
      Vector const &       D      = Vector(),
      Vector const &       c      = Vector() )
      : m_A( A ), m_m( A.rows() ), m_n( A.cols() ), m_lambda( lambda ), m_D2( D ), m_c( c )
    {
      UTILS_ASSERT( lambda >= 0, "SP_TikhonovSolver2( λ={} ) λ must be >= 0", lambda );

      m_has_reg = lambda > epsilon();

      if ( D.size() > 0 )
      {
        UTILS_ASSERT( D.size() == m_n, "SP_TikhonovSolver2: D size {} must match A cols {}", D.size(), m_n );
        UTILS_ASSERT( D.minCoeff() >= 0, "SP_TikhonovSolver2: D elements must be >= 0" );
        m_use_diag = true;
        m_D2.array() *= m_D2.array();
      }

      if ( c.size() > 0 )
      {
        UTILS_ASSERT( c.size() == m_n, "SP_TikhonovSolver2: c size {} must match A cols {}", c.size(), m_n );
        m_has_c = true;
      }

      build_and_factorize();
    }

    /// @brief Convenience constructor with D and λ = 1
    explicit SP_TikhonovSolver2( SparseMatrix const & A, Vector const & D ) : SP_TikhonovSolver2( A, Scalar( 1 ), D ) {}

    /**
     * @brief Solve sparse Tikhonov problem via KKT system
     *
     * @param b Observation vector
     * @return Solution vector
     */
    Vector solve( Vector const & b ) const
    {
      UTILS_ASSERT( b.size() == m_m, "SP_TikhonovSolver2::solve: b size {} must match A rows {}", b.size(), m_m );

      Vector rhs( m_m + m_n );
      rhs.head( m_m ) = b;

      if ( m_has_c && m_has_reg )
      {
        if ( m_use_diag )
          rhs.tail( m_n ) = m_lambda * m_D2.cwiseProduct( m_c );
        else
          rhs.tail( m_n ) = m_lambda * m_c;
      }
      else
      {
        rhs.tail( m_n ).setZero();
      }

      Vector sol;
      if ( m_use_LDLT )
      {
        sol = m_LDLT.solve( rhs );
        UTILS_ASSERT( m_LDLT.info() == Eigen::Success, "SP_TikhonovSolver2::solve (LDLT): solve failed" );
      }
      else
      {
        sol = m_LU.solve( rhs );
        UTILS_ASSERT( m_LU.info() == Eigen::Success, "SP_TikhonovSolver2::solve (LU): solve failed" );
      }

      return sol.segment( m_m, m_n );
    }

    /// @name Accessors and Debug Methods
    /// @{
    SparseMatrix const & get_KKT() const { return m_KKT; }
    Scalar               lambda() const { return m_lambda; }
    Vector const &       diagonal2() const { return m_D2; }
    Vector const &       center() const { return m_c; }
    bool                 uses_diagonal() const { return m_use_diag; }
    bool                 has_center() const { return m_has_c; }
    bool                 has_regularization() const { return m_has_reg; }
    bool                 uses_LDLT() const { return m_use_LDLT; }
    integer              rows() const { return m_m; }
    integer              cols() const { return m_n; }
    /// @}
  };

}  // namespace Utils

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif

//
// eof: Utils_Tikhonov.hh
//
