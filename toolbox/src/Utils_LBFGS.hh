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
// file: Utils_LBFGS.hh
//

#pragma once

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifndef UTILS_LBFGS_dot_HH
#define UTILS_LBFGS_dot_HH

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
#include "Utils_fmt.hh"
#include "Utils_eigen.hh"

#include <optional>
#include <cmath>
#include <algorithm>
#include <limits>
#include <utility>

/**
 * @file Utils_LBFGS.hh
 * @brief Complete header-only implementation of L-BFGS optimization algorithms
 *
 * This file provides a comprehensive implementation of the Limited-memory 
 * Broyden-Fletcher-Goldfarb-Shanno (L-BFGS) quasi-Newton optimization method,
 * including multiple line search strategies and support for box constraints.
 *
 * ## Main Components
 *
 * - **LBFGS**: Core L-BFGS storage and two-loop recursion algorithm
 * - **Line Search Policies**: Multiple strategies (Armijo, Wolfe, Strong Wolfe, 
 *   Goldstein, Hager-Zhang, More-Thuente)
 * - **LBFGS_minimizer**: High-level minimizer with automatic line search and 
 *   projected gradient methods for box-constrained optimization
 *
 * ## Theoretical Background
 *
 * The L-BFGS method approximates the inverse Hessian matrix using a limited number
 * of recent gradient differences. This makes it suitable for large-scale optimization
 * where storing the full Hessian is impractical.
 *
 * ### Algorithm Overview
 *
 * At iteration k, L-BFGS computes a search direction \f$p_k\f$ by:
 * \f[
 *   p_k = -H_k \nabla f(x_k)
 * \f]
 * where \f$H_k\f$ is an approximation to \f$[\nabla^2 f(x_k)]^{-1}\f$ constructed
 * from the m most recent correction pairs \f$(s_i, y_i)\f$:
 * \f[
 *   s_i = x_{i+1} - x_i, \quad y_i = \nabla f(x_{i+1}) - \nabla f(x_i)
 * \f]
 *
 * The key advantage is that \f$H_k\f$ is never formed explicitly; instead,
 * the product \f$H_k g\f$ is computed efficiently via two-loop recursion in O(mn) time.
 *
 * ## References
 *
 * -# J. Nocedal (1980). "Updating Quasi-Newton Matrices with Limited Storage".
 *    Mathematics of Computation, 35(151), 773-782.
 *    DOI: 10.1090/S0025-5718-1980-0572855-7
 *
 * -# D.C. Liu and J. Nocedal (1989). "On the Limited Memory BFGS Method for 
 *    Large Scale Optimization". Mathematical Programming, 45(1-3), 503-528.
 *    DOI: 10.1007/BF01589116
 *
 * -# R.H. Byrd, P. Lu, J. Nocedal, and C. Zhu (1995). "A Limited Memory Algorithm 
 *    for Bound Constrained Optimization". SIAM Journal on Scientific Computing,
 *    16(5), 1190-1208. DOI: 10.1137/0916069
 *
 * -# J. Nocedal and S.J. Wright (2006). "Numerical Optimization", 2nd Edition,
 *    Springer. ISBN: 978-0-387-30303-1
 *
 * -# J.J. Moré and D.J. Thuente (1994). "Line Search Algorithms with Guaranteed
 *    Sufficient Decrease". ACM Transactions on Mathematical Software, 20(3), 286-307.
 *    DOI: 10.1145/192115.192132
 *
 * -# W.W. Hager and H. Zhang (2006). "A New Conjugate Gradient Method with 
 *    Guaranteed Descent and an Efficient Line Search". SIAM Journal on Optimization,
 *    16(1), 170-192. DOI: 10.1137/030601880
 *
 * ## Implementation Features
 *
 * - **Robustness**: Includes Powell damping, curvature checks, and fallback strategies
 * - **Flexibility**: Template-based design supports different scalar types (float/double)
 * - **Efficiency**: Uses Eigen for vectorized operations, minimal memory allocations
 * - **Box Constraints**: Implements projected gradient methods for bound-constrained problems
 * - **Multiple Line Searches**: Choose the best strategy for your problem characteristics
 *
 * ## Usage Example
 *
 * @code{.cpp}
 * using namespace Utils;
 * 
 * // Define objective function
 * auto rosenbrock = [](Vector const& x, Vector* g) -> double {
 *   double f = 100*(x(1)-x(0)*x(0))*(x(1)-x(0)*x(0)) + (1-x(0))*(1-x(0));
 *   if (g) {
 *     (*g)(0) = -400*(x(1)-x(0)*x(0))*x(0) - 2*(1-x(0));
 *     (*g)(1) = 200*(x(1)-x(0)*x(0));
 *   }
 *   return f;
 * };
 * 
 * // Setup minimizer
 * LBFGS_minimizer<double>::Options opts;
 * opts.max_iter = 1000;
 * opts.g_tol = 1e-6;
 * opts.verbose = true;
 * 
 * LBFGS_minimizer<double> minimizer(opts);
 * 
 * // Initial point
 * Vector x0(2);
 * x0 << -1.2, 1.0;
 * 
 * // Minimize
 * auto [status, x_opt, f_opt, data] = minimizer.minimize(
 *   x0, rosenbrock, StrongWolfeLineSearch<double>()
 * );
 * 
 * std::cout << "Solution: " << x_opt.transpose() << std::endl;
 * std::cout << "Minimum: " << f_opt << std::endl;
 * @endcode
 *
 * @author Enrico Bertolazzi
 * @date 2025
 */

namespace Utils {

  using std::abs;
  using std::min;
  using std::max;

  // ===========================================================================
  // LBFGS: Core Two-Loop Recursion Implementation
  // ===========================================================================

  /**
   * @class LBFGS
   * @brief Limited-memory BFGS storage and two-loop recursion implementation
   *
   * This class implements the core L-BFGS algorithm for computing approximate
   * inverse Hessian-vector products without explicitly storing the matrix.
   *
   * ## Algorithm Description
   *
   * The L-BFGS method maintains a history of m correction pairs:
   * \f[
   *   \{(s_i, y_i)\}_{i=k-m}^{k-1}
   * \f]
   * where:
   * - \f$s_i = x_{i+1} - x_i\f$ (displacement vector)
   * - \f$y_i = \nabla f(x_{i+1}) - \nabla f(x_i)\f$ (gradient difference)
   *
   * ### Two-Loop Recursion
   *
   * Given gradient g, compute \f$H_k g\f$ via:
   *
   * **First loop** (backward, from newest to oldest):
   * \f{align*}{
   *   q &\leftarrow g \\
   *   \text{for } i &= k-1, k-2, \ldots, k-m: \\
   *     \alpha_i &\leftarrow \rho_i s_i^T q \\
   *     q &\leftarrow q - \alpha_i y_i
   * \f}
   *
   * **Scaling**:
   * \f[
   *   r \leftarrow H_0 q
   * \f]
   * where \f$H_0 = \gamma I\f$ with \f$\gamma = \frac{s_{k-1}^T y_{k-1}}{y_{k-1}^T y_{k-1}}\f$
   *
   * **Second loop** (forward, from oldest to newest):
   * \f{align*}{
   *   \text{for } i &= k-m, k-m+1, \ldots, k-1: \\
   *     \beta_i &\leftarrow \rho_i y_i^T r \\
   *     r &\leftarrow r + s_i(\alpha_i - \beta_i)
   * \f}
   *
   * **Output**: \f$r = H_k g\f$
   *
   * ### Computational Complexity
   *
   * - Time: O(mn) per iteration, where n is dimension and m is memory size
   * - Space: O(mn) to store correction pairs
   *
   * ### Curvature Condition
   *
   * For numerical stability, pairs are only accepted if:
   * \f[
   *   s_i^T y_i > \epsilon \|s_i\| \|y_i\|
   * \f]
   * This ensures the BFGS update maintains positive definiteness.
   *
   * ## Implementation Details
   *
   * - Uses circular buffer for efficient O(1) insertion of new pairs
   * - Automatic oldest pair removal when capacity is exceeded
   * - Robust curvature checking with both relative and absolute thresholds
   * - Optional Powell damping for improved robustness
   * - Clamping of initial Hessian approximation to avoid extreme scaling
   *
   * @tparam Scalar Floating-point type (typically double or float)
   *
   * @note This implementation is thread-safe for read operations if no
   *       concurrent modifications occur.
   *
   * @see Liu, D.C. and Nocedal, J. (1989) for the original algorithm
   * @see Nocedal and Wright (2006), Chapter 7, for detailed theory
   */

  template <typename Scalar = double>
  class LBFGS {
  public:
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

  private:
    size_t m_capacity;         ///< Maximum number of correction pairs (m)
    size_t m_dimension;        ///< Problem dimension (n)
    Matrix m_S;                ///< Storage for s vectors (n × m matrix, each column is s_i)
    Matrix m_Y;                ///< Storage for y vectors (n × m matrix, each column is y_i)
    Vector m_rho;              ///< Storage for ρ_i = 1/(y_i^T s_i) values
    size_t m_current_size{0};  ///< Current number of stored pairs
    size_t m_oldest_index{0};  ///< Circular buffer index of oldest pair
    size_t m_newest_index{0};  ///< Circular buffer index where next pair will be stored

    /* Robustness parameters */
    bool   m_enable_damping{ true };        ///< Enable Powell-style damping by default
    Scalar m_h0_min{ Scalar(1e-6) };        ///< Minimum allowed initial H0 scaling
    Scalar m_h0_max{ Scalar(1e6) };         ///< Maximum allowed initial H0 scaling

    /* Temporary workspace vectors (mutable for const methods) */
    mutable Vector m_alpha;  ///< Workspace for α_i values in two-loop recursion
    mutable Vector m_q;      ///< Workspace for q vector in two-loop recursion
    mutable Vector m_r;      ///< Workspace for r vector in two-loop recursion

  public:
    /**
     * @brief Construct L-BFGS storage
     *
     * Allocates memory for storing up to maxCorrections correction pairs.
     * The circular buffer strategy allows O(1) insertion of new pairs.
     *
     * @param maxCorrections Maximum number of stored correction pairs (typical: 5-20).
     *                       Larger values provide better Hessian approximation but
     *                       increase memory usage and computational cost.
     * @param problem_dimension Dimension of the optimization problem (n)
     *
     * @note A typical choice is m=10 for most problems. Use smaller m (3-7) for
     *       very large problems or when memory is constrained.
     */
    explicit
    LBFGS( size_t maxCorrections = 10, size_t problem_dimension = 0 )
    : m_capacity( maxCorrections )
    , m_dimension( problem_dimension )
    , m_S( problem_dimension, maxCorrections )
    , m_Y( problem_dimension, maxCorrections )
    , m_rho( maxCorrections )
    , m_alpha( maxCorrections )
    , m_q( problem_dimension )
    , m_r( problem_dimension )
    {
      assert( maxCorrections > 0 );
      clear();
    }

    /**
     * @brief Clear all stored correction pairs
     *
     * Resets the L-BFGS memory to its initial empty state. This is useful when:
     * - Starting optimization of a new problem
     * - Restarting after convergence issues
     * - Periodic refresh to avoid accumulation of numerical errors
     *
     * @post size() == 0 and all internal storage is zeroed
     */
    void
    clear() {
      m_current_size = 0;
      m_oldest_index = 0;
      m_newest_index = 0;
      m_S.setZero();
      m_Y.setZero();
      m_rho.setZero();
      m_alpha.setZero();
      m_q.setZero();
      m_r.setZero();
    }

    /**
     * @brief Return current number of stored correction pairs
     * @return Number of pairs currently in memory (0 ≤ size() ≤ capacity())
     */
    size_t size() const { return m_current_size; }

    /**
     * @brief Return maximum capacity for correction pairs
     * @return Maximum number of pairs (m parameter)
     */
    size_t capacity() const { return m_capacity; }

    /**
     * @brief Return problem dimension
     * @return Dimension of vectors in optimization problem (n)
     */
    size_t dimension() const { return m_dimension; }

    /**
     * @brief Resize for new problem dimension
     *
     * Reallocates storage for a different problem size. This clears all
     * existing correction pairs.
     *
     * @param new_dimension New problem dimension
     *
     * @post size() == 0 and dimension() == new_dimension
     */
    void
    resize( size_t const new_dimension ) {
      m_dimension = new_dimension;
      m_S.resize( new_dimension, m_capacity );
      m_Y.resize( new_dimension, m_capacity );
      m_S.setZero();
      m_Y.setZero();
      clear();
    }

    /**
     * @brief Add a correction pair (s, y) to L-BFGS memory
     *
     * Attempts to add a new correction pair after validating the curvature condition.
     * The pair is accepted only if:
     * \f[
     *   s^T y > \max\{\epsilon_{\text{rel}} \|s\| \|y\|, \epsilon_{\text{abs}} \|s\|^2\}
     * \f]
     *
     * This ensures:
     * - Positive definiteness of the Hessian approximation
     * - Numerical stability (avoids division by very small numbers)
     * - Well-conditioned updates
     *
     * ### Update Strategy
     *
     * When the buffer is full (size() == capacity()), the oldest pair is
     * automatically discarded to make room for the new one.
     *
     * @param s Displacement vector \f$s = x_{k+1} - x_k\f$
     * @param y Gradient difference \f$y = \nabla f(x_{k+1}) - \nabla f(x_k)\f$
     * @param min_curvature_ratio Minimum ratio \f$\frac{s^T y}{\|s\| \|y\|}\f$ 
     *                            for acceptance (default: 1e-8)
     *
     * @return true if pair was accepted and stored, false if rejected due to
     *         insufficient curvature
     *
     * @pre s.size() == y.size()
     * @post If accepted: size() increments by 1 (unless at capacity)
     *
     * @note Rejection typically occurs when:
     *       - Line search is too inaccurate
     *       - Gradient evaluations have numerical errors
     *       - Problem has very small or zero curvature in search direction
     *
     * @see Nocedal & Wright (2006), Section 7.2 for curvature condition theory
     */
    bool
    add_correction( Vector const & s, Vector const & y, Scalar const min_curvature_ratio = 1e-8 ) {
      assert( s.size() == y.size() );
      if ( static_cast<size_t>(s.size()) != m_dimension ) resize( s.size() );

      Scalar const sty  { s.dot(y) };
      Scalar const snrm { s.norm() };
      Scalar const ynrm { y.norm() };

      // basic positivity check
      if ( !(sty > 0) ) return false;

      // relative and absolute thresholds for curvature test (robust)
      Scalar const eps        = std::numeric_limits<Scalar>::epsilon();
      Scalar const rel_thresh = min_curvature_ratio * snrm * ynrm;
      Scalar const abs_thresh = eps * snrm * snrm;
      Scalar const thresh     = max(rel_thresh, abs_thresh);

      if ( sty <= thresh ) return false;

      // Store the new pair
      m_S.col ( m_newest_index ) = s;
      m_Y.col ( m_newest_index ) = y;
      m_rho   ( m_newest_index ) = Scalar(1.0) / sty;

      // Update indices and size
      if ( m_current_size < m_capacity ) {
        // Buffer not full yet
        ++m_current_size;
        m_newest_index = (m_newest_index + 1) % m_capacity;
      } else {
        // Buffer full - overwrite oldest
        m_oldest_index = (m_oldest_index + 1) % m_capacity;
        m_newest_index = (m_newest_index + 1) % m_capacity;
      }

      return true;
    }

    /**
     * @brief Compute H*g using two-loop recursion algorithm
     *
     * This is the core L-BFGS operation that computes the product of the
     * approximate inverse Hessian H with a gradient vector g, without ever
     * forming H explicitly.
     *
     * ### Algorithm Steps
     *
     * 1. **First Loop** (backward through history):
     *    - Start with q = g
     *    - For each pair (s_i, y_i) from newest to oldest:
     *      - Compute α_i = ρ_i (s_i · q)
     *      - Update q ← q - α_i y_i
     *
     * 2. **Initial Hessian Application**:
     *    - Scale: r = h0 * q
     *    - This represents H_0 = h0 * I (scaled identity matrix)
     *
     * 3. **Second Loop** (forward through history):
     *    - For each pair (s_i, y_i) from oldest to newest:
     *      - Compute β_i = ρ_i (y_i · r)
     *      - Update r ← r + s_i (α_i - β_i)
     *
     * 4. **Result**: Return r = H * g
     *
     * ### Computational Details
     *
     * - **Time Complexity**: O(m*n), where m = size() and n = dimension()
     * - **Space Complexity**: O(m + n) for temporary vectors
     * - **Memory Access**: Efficient cache usage due to contiguous storage
     *
     * @param g Gradient vector (typically ∇f(x))
     * @param h0 Initial diagonal Hessian approximation scalar. Common choices:
     *           - h0 = 1.0 (simple scaling)
     *           - h0 = (s^T y)/(y^T y) using most recent pair (recommended)
     *           - Adaptive scaling based on problem characteristics
     *
     * @return Vector r = H * g, representing the quasi-Newton search direction
     *
     * @pre g.size() == dimension()
     * @post result.size() == g.size()
     *
     * @note If no pairs are stored (size() == 0), returns h0 * g (scaled steepest descent)
     *
     * @see Nocedal (1980) for original two-loop recursion algorithm
     * @see Liu & Nocedal (1989) for practical implementation details
     */
    Vector
    two_loop_recursion( Vector const & g, Scalar h0 ) const {
      if ( m_current_size == 0 ) return h0 * g;

      m_q = g; // Start with q = g

      // =====================================================================
      // First loop: Process pairs from newest to oldest
      // =====================================================================
      size_t idx = (m_newest_index == 0) ? m_capacity - 1 : m_newest_index - 1;
        
      size_t i{m_current_size};
      while ( i > 0 ) {
        --i;
        m_alpha(i) = m_rho(idx) * m_S.col(idx).dot(m_q);
        m_q       -= m_alpha(i) * m_Y.col(idx);
            
        // Move to previous index in circular buffer
        idx = (idx == 0) ? m_capacity - 1 : idx - 1;
      }

      // =====================================================================
      // Apply initial Hessian approximation H0 = h0 * I
      // =====================================================================
      m_r = h0 * m_q;

      // =====================================================================
      // Second loop: Process pairs from oldest to newest
      // =====================================================================
      idx = m_oldest_index;
      for ( size_t i{0}; i < m_current_size; ++i ) {
        Scalar beta = m_rho(idx) * m_Y.col(idx).dot(m_r);
        m_r += m_S.col(idx) * (m_alpha(i) - beta);

        // Move to next index in circular buffer
        idx = (idx + 1) % m_capacity;
      }
      
      return m_r;
    }

    /**
     * @brief Compute recommended initial Hessian scaling from latest pair
     *
     * Computes the scalar h0 for the initial Hessian approximation H0 = h0*I
     * using the most recent correction pair:
     * \f[
     *   h_0 = \frac{s_{k-1}^T y_{k-1}}{y_{k-1}^T y_{k-1}}
     * \f]
     *
     * This choice is motivated by:
     * - Self-scaling property: approximates the curvature along the most recent step
     * - Theoretical foundation: relates to Barzilai-Borwein spectral steplength
     * - Empirical success: works well in practice for many problems
     *
     * The result is clamped to [m_h0_min, m_h0_max] to prevent extreme scaling
     * that could lead to numerical instability.
     *
     * @param default_value Value to return if no pairs are stored (default: 1.0)
     *
     * @return Recommended h0 value, or default_value if size() == 0
     *
     * @note The clamping bounds are:
     *       - Lower: 1e-6 (prevents division issues for nearly-flat regions)
     *       - Upper: 1e6 (prevents overflow in steep regions)
     *
     * @see Nocedal & Wright (2006), Section 7.2 for scaling strategies
     */
    Scalar
    compute_initial_h0( Scalar default_value = Scalar(1.0) ) const {
      if ( m_current_size == 0 ) return default_value;
        
      // Get the newest pair (just before newest_index in circular buffer)
      size_t latest_idx = (m_newest_index == 0) ? m_capacity - 1 : m_newest_index - 1;
      auto const & s = m_S.col( latest_idx );
      auto const & y = m_Y.col( latest_idx );
        
      Scalar sty { s.dot(y) };
      Scalar yty { y.dot(y) };
      Scalar h0 = default_value;
      if ( yty > 0 ) h0 = sty / yty;
      
      // Clamp to safe range
      h0 = max(m_h0_min, min(m_h0_max, h0));
      return h0;
    }

    /**
     * @brief Add correction with Powell damping for improved robustness
     *
     * Powell damping (also called modified BFGS) adjusts the gradient difference
     * y to ensure the curvature condition is satisfied, even when the standard
     * pair (s,y) would be rejected.
     *
     * ### Algorithm
     *
     * When s^T y is too small (less than threshold), compute:
     * \f[
     *   \theta = \frac{\text{threshold} - s^T B_0 s}{s^T y - s^T B_0 s}
     * \f]
     * where B_0 = H_0^{-1} = (1/h0)*I, then use damped gradient:
     * \f[
     *   \hat{y} = \theta y + (1-\theta) B_0 s
     * \f]
     *
     * This creates a convex combination that guarantees s^T ŷ ≥ threshold.
     *
     * ### Benefits
     *
     * - Maintains positive definiteness even with inaccurate line searches
     * - Allows continued progress in difficult regions
     * - Reduces need for memory resets
     *
     * @param lb Reference to LBFGS object to update
     * @param s Displacement vector
     * @param y Gradient difference vector
     * @param min_curvature_ratio Minimum curvature ratio (default: 1e-8)
     *
     * @return true if correction was successfully added (possibly damped)
     *
     * @note Falls back to undamped update if damping is unnecessary or fails
     *
     * @see Powell, M.J.D. (1978). "A Fast Algorithm for Nonlinearly Constrained
     *      Optimization Calculations". Numerical Analysis, Lecture Notes in
     *      Mathematics, Springer.
     * @see Nocedal & Wright (2006), Section 18.3 for modified BFGS
     */
    bool
    add_correction_with_damping(
      LBFGS<Scalar> & lb,
      Vector const  & s,
      Vector const  & y,
      Scalar const   min_curvature_ratio = 1e-8
    ) {
      using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
      Scalar const sty        = s.dot(y);
      Scalar const snrm       = s.norm();
      Scalar const ynrm       = y.norm();
      Scalar const eps        = std::numeric_limits<Scalar>::epsilon();
      Scalar const rel_thresh = min_curvature_ratio * snrm * ynrm;
      Scalar const abs_thresh = eps * snrm * snrm;
      Scalar const thresh     = max(rel_thresh, abs_thresh);

      if ( sty > thresh ) return lb.add_correction(s,y,min_curvature_ratio);

      // compute h0 and fallback if invalid
      Scalar h0 = lb.compute_initial_h0(1);
      if ( !(h0 > 0) ) return false;

      Scalar sBs = s.squaredNorm() / h0; // s^T (B0 s) with B0 = H0^{-1} = (1/h0) I

      // If denominator would be zero or very small, don't attempt damping
      Scalar denom = sty - sBs;
      Vector y_hat(y);
      if ( abs(denom) < std::numeric_limits<Scalar>::epsilon() ) {
        // ambiguous, don't damp
        return lb.add_correction(s,y,min_curvature_ratio);
      } else {
        // choose theta so that s^T y_hat = thresh (solve theta*(sty - sBs) + sBs = thresh)
        Scalar theta = (thresh - sBs) / denom;
        // clamp to [0,1]
        if      ( theta < 0 ) theta = 0;
        else if ( theta > 1 ) theta = 1;
        y_hat = theta * y + (1 - theta) * (s / h0);
        Scalar sty_hat = s.dot(y_hat);
        if ( !(sty_hat > 0) ) return false;
        return lb.add_correction(s, y_hat, min_curvature_ratio);
      }
    }

  };

  // ===========================================================================
  // Line Search Helper Functions
  // ===========================================================================

  namespace detail {
    /**
     * @brief Robust cubic interpolation for line-search (Hager–Zhang friendly).
     *
     * This version provides:
     *  - Full finite checks (NaN/inf)
     *  - Degeneracy detection (a ≈ b, fa ≈ fb, etc.)
     *  - Hager–Zhang safeguard window: [a_lo+δΔ, a_hi−σΔ]
     *  - Fallback to bisection in all doubtful cases
     */
    template<class Scalar>
    inline
    Scalar
    cubic_minimizer(
      Scalar a,  Scalar fa,  Scalar fpa,
      Scalar b,  Scalar fb,  Scalar fpb,
      Scalar delta = Scalar(0.1),   // Hager–Zhang recommended
      Scalar sigma = Scalar(0.9)    // Hager–Zhang recommended
    ) {
      // reorder interval
      Scalar lo = min(a,b), hi = max(a,b);
      Scalar width = hi - lo;
    
      // fallback if degenerate interval
      if (width <= 10 * std::numeric_limits<Scalar>::epsilon() * (1+abs(lo)))
          return (a + b) / 2;
    
      // -----------------------------------------------------------------------
      // Standard More–Thuente cubic setup
      // -----------------------------------------------------------------------
      Scalar denom_ab = (a - b);
      if (abs(denom_ab) < std::numeric_limits<Scalar>::epsilon())
        return (a + b) / 2;
    
      Scalar d1 = fpa + fpb - 3 * ((fa - fb) / denom_ab);
      Scalar discr = d1*d1 - fpa*fpb;
    
      // discriminant negative → no real minimizer → bisection
      if ( discr < 0 ) return (a + b) / 2;
    
      Scalar sqrt_disc = sqrt(max(discr, Scalar(0)));
    
      Scalar denom = fpb - fpa + 2 * sqrt_disc;
      if (!std::isfinite(denom) || abs(denom) < std::numeric_limits<Scalar>::epsilon())
        return (a + b) / 2;
    
      // cubic minimizer
      Scalar t = b - (b - a) * ((fpb + sqrt_disc - d1) / denom);
    
      // -----------------------------------------------------------------------
      // Hager–Zhang safeguard: force t inside:
      //   [ a_lo + δΔ , a_hi − σΔ ]
      // -----------------------------------------------------------------------
      Scalar left  = lo + delta * width;
      Scalar right = lo + sigma * width;
    
      if (!std::isfinite(t) || t <= left || t >= right)
        return (a + b) / 2;
    
      return t;
    }

    /**
     * @brief Compute next trial step for More-Thuente line search
     *
     * This function implements the step selection logic from More-Thuente,
     * choosing between cubic interpolation, quadratic interpolation, and
     * bisection based on the current state of the line search.
     *
     * ### Strategy
     *
     * 1. **Primary method**: Cubic interpolation using interval endpoints
     * 2. **Safeguarding**: Switch to bisection or extrapolation if cubic fails
     * 3. **Phase-dependent logic**:
     *    - Bracketing phase: Use extrapolation if cubic invalid
     *    - Zoom phase: Use bisection if cubic invalid
     *
     * @param a_lo Lower bound of current interval with function/derivative info
     * @param phi_lo Function value at a_lo
     * @param der_lo Derivative at a_lo
     * @param a_hi Upper bound of current interval
     * @param phi_hi Function value at a_hi
     * @param der_hi Derivative at a_hi
     * @param a_prev Previous trial point
     * @param phi_prev Function value at previous point
     * @param der_prev Derivative at previous point
     * @param alpha_max Maximum allowable step length
     * @param step_max Maximum step to take in this iteration
     * @param is_bracketing true if in bracketing phase, false if in zoom phase
     *
     * @return Next trial step length
     *
     * @note This is a simplified but robust version of MT step computation
     */
    template <typename Scalar>
    Scalar
    compute_step(
      Scalar a_lo,   Scalar phi_lo,   Scalar der_lo,
      Scalar a_hi,   Scalar phi_hi,   Scalar der_hi,
      Scalar a_prev, Scalar phi_prev, Scalar der_prev,
      Scalar alpha_max, Scalar step_max,
      bool is_bracketing
    ) {
      Scalar a_j = cubic_minimizer(a_lo, phi_lo, der_lo, a_hi, phi_hi, der_hi);
      if (a_j <= a_lo || a_j >= a_hi) {
        a_j = is_bracketing ? (a_hi + 2*(a_hi - a_prev)) : ((a_lo + a_hi)/2);
        if(a_j > alpha_max) a_j = alpha_max;
      }
      Scalar tol = 1e-6;
      a_j = std::clamp(a_j, a_lo+tol*(a_hi-a_lo), a_hi-tol*(a_hi-a_lo));
      return min(a_j, step_max);
    }
     
  } // namespace detail


  // ===========================================================================
  // Line Search Algorithms
  // ===========================================================================

  /**
   * @class ArmijoLineSearch
   * @brief Simple Armijo (sufficient decrease) line search
   *
   * Implements the classical Armijo backtracking line search, which finds a
   * step length satisfying the Armijo condition:
   * \f[
   *   f(x + \alpha p) \leq f(x) + c_1 \alpha \nabla f(x)^T p
   * \f]
   *
   * ### Algorithm
   *
   * 1. Start with initial step α₀ (typically 1.0)
   * 2. Evaluate f(x + α p)
   * 3. If Armijo condition satisfied, accept α
   * 4. Otherwise, reduce α ← ρ α (ρ ≈ 0.5)
   * 5. Repeat until convergence or max iterations
   *
   * ### Characteristics
   *
   * **Advantages:**
   * - Very simple and fast per iteration
   * - Robust: always terminates if gradient is descent direction
   * - No derivative evaluations at trial points needed
   * - Works well for Newton-like methods with good curvature
   *
   * **Disadvantages:**
   * - May accept steps that are too short (inefficient progress)
   * - No guarantee of adequate step length for superlinear convergence
   * - Can lead to many small steps in poorly scaled problems
   *
   * ### When to Use
   *
   * Armijo is recommended for:
   * - Newton or quasi-Newton methods with good Hessian approximation
   * - Problems where function evaluations are cheap
   * - When simplicity and reliability are more important than efficiency
   * - Initial iterations of optimization (before switching to stronger conditions)
   *
   * @tparam Scalar Floating-point type (double or float)
   *
   * @note This is the simplest practical line search. For better efficiency,
   *       consider Wolfe or Strong Wolfe conditions.
   *
   * @see Armijo, L. (1966). "Minimization of functions having Lipschitz 
   *      continuous first partial derivatives". Pacific Journal of Mathematics.
   * @see Nocedal & Wright (2006), Algorithm 3.1
   */
  template <typename Scalar>
  class ArmijoLineSearch {
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    Scalar m_c1          = 1e-4;   ///< Armijo constant (typical: 1e-4)
    Scalar m_step_reduce = 0.5;    ///< Step reduction factor when Armijo fails
    Scalar m_step_expand = 1.2;    ///< Step expansion factor (currently unused)
    size_t m_max_iters   = 50;     ///< Maximum backtracking iterations
    Scalar m_epsi        = 1e-15;  ///< Minimum step size threshold

    mutable Vector m_x_new;        ///< Workspace for trial point

  public:
    /**
     * @brief Perform Armijo backtracking line search
     *
     * @tparam Callback Function type with signature: Scalar(Vector const&, Vector*)
     *
     * @param f0 Function value at current point x
     * @param Df0 Directional derivative ∇f(x)^T d (must be negative)
     * @param x Current point
     * @param d Search direction
     * @param callback Function that evaluates f (and optionally ∇f)
     * @param step0 Initial trial step length (default: 1.0)
     *
     * @return Step length if successful, std::nullopt if failed
     *
     * @note If line search fails repeatedly, returns best step found even if
     *       it doesn't strictly satisfy Armijo
     */
    template <typename Callback>
    std::optional<Scalar>
    operator()(
      Scalar           f0,
      Scalar           Df0,
      Vector   const & x,
      Vector   const & d,
      Callback const & callback,
      Scalar           step0 = 1
    ) const {
      Scalar step = step0;
      Scalar c1_Df0 = m_c1 * Df0;
      auto const n{ x.size() };
      m_x_new.resize(n);
      
      Scalar best_step    = step0;
      Scalar best_f       = std::numeric_limits<Scalar>::max();
      size_t shrink_count = 0;
      
      for (size_t k = 0; k < m_max_iters; ++k) {
        m_x_new.noalias() = x + step * d;
        Scalar f_new = callback(m_x_new, nullptr);
        
        // Track best step found
        if (f_new < best_f) {
          best_f    = f_new;
          best_step = step;
        }
        
        if ( f_new <= f0 + step * c1_Df0 ) return step; // Success
        
        // Adaptive reduction: more aggressive after multiple failures
        Scalar reduction = m_step_reduce;
        if ( shrink_count > 2 ) reduction = 0.1;
            
        step *= reduction;
        shrink_count++;
            
        // Check for excessively small steps
        if ( step < m_epsi ) {
          // Return best step found even if it doesn't satisfy Armijo
          if ( best_f < f0 ) return best_step;
          break;
        }
      }
      
      // Fallback: return best step if it improved objective
      if ( best_f < f0 ) return best_step;

      return std::nullopt;
    }
  };

  /**
   * @class WolfeLineSearch
   * @brief Wolfe conditions line search (sufficient decrease + curvature)
   *
   * Implements a line search satisfying the Wolfe conditions:
   * \f{align*}{
   *   f(x + \alpha p) &\leq f(x) + c_1 \alpha \nabla f(x)^T p \quad &\text{(Armijo)}\\
   *   \nabla f(x + \alpha p)^T p &\geq c_2 \nabla f(x)^T p \quad &\text{(Curvature)}
   * \f}
   * with \f$0 < c_1 < c_2 < 1\f$.
   *
   * ### Algorithm Structure
   *
   * **Phase 1 - Bracketing**: Find interval [α_lo, α_hi] containing acceptable step
   * - Start with α = α₀
   * - Increase α until Armijo violated or curvature satisfied
   * - If Armijo violated, bracket found
   * - If derivative becomes positive, bracket found
   *
   * **Phase 2 - Zoom**: Refine within bracket using cubic interpolation
   * - Repeatedly subdivide interval
   * - Choose new point via cubic interpolation
   * - Update bracket endpoints
   * - Terminate when curvature condition satisfied
   *
   * ### Characteristics
   *
   * **Advantages:**
   * - Guarantees superlinear convergence for quasi-Newton methods
   * - Prevents steps that are too short (via curvature condition)
   * - Well-suited for L-BFGS and conjugate gradient methods
   * - Theoretically sound with convergence guarantees
   *
   * **Disadvantages:**
   * - Requires gradient evaluations at trial points
   * - More complex than Armijo
   * - May require several iterations to converge
   *
   * ### Parameters
   *
   * Typical values:
   * - c₁ = 1e-4 (Armijo parameter)
   * - c₂ = 0.9 for Newton/quasi-Newton, 0.1 for nonlinear CG
   *
   * @tparam Scalar Floating-point type (double or float)
   *
   * @note For L-BFGS, use c₂ = 0.9. For conjugate gradient, use c₂ = 0.1-0.4.
   *
   * @see Nocedal & Wright (2006), Algorithm 3.5 and 3.6
   * @see Wolfe, P. (1969). "Convergence Conditions for Ascent Methods"
   */
  template <typename Scalar>
  class WolfeLineSearch {
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    Scalar m_c1        = 1e-4;  ///< Armijo constant
    Scalar m_c2        = 0.9;   ///< Curvature constant
    Scalar m_alpha_max = 10;    ///< Maximum step length
    Scalar m_epsi      = 1e-12; ///< Convergence tolerance for interval
    size_t m_max_iters = 50;    ///< Maximum total iterations
    
    mutable Vector m_g_new, m_x_new;  ///< Workspace vectors
    
  public:

    /**
     * @brief Perform Wolfe conditions line search
     *
     * @tparam Callback Function type: Scalar(Vector const&, Vector*)
     *                  Must compute both f and ∇f when gradient pointer is non-null
     *
     * @param f0 Function value at x
     * @param Df0 Directional derivative ∇f(x)^T d (must be < 0)
     * @param x Current point
     * @param d Search direction
     * @param callback Function/gradient evaluator
     * @param alpha0 Initial step length (default: 1.0)
     *
     * @return Step length satisfying Wolfe conditions, or std::nullopt if failed
     *
     * @pre Df0 < 0 (d must be descent direction)
     */
    template <typename Callback>
    std::optional<Scalar>
    operator()(
      Scalar           f0,
      Scalar           Df0,
      Vector   const & x,
      Vector   const & d,
      Callback const & callback,
      Scalar           alpha0 = 1
    ) const {

      if ( !(Df0 < 0) ) return std::nullopt;
      m_g_new.resize( x.size() );
      m_x_new.resize( x.size() );
      size_t evals = 0;

      auto eval = [&]( Scalar a, Scalar & f, Scalar & df ) {
        ++evals;
        m_x_new.noalias() = x + a * d;
        f  = callback( m_x_new, &m_g_new );
        df = m_g_new.dot(d);
      };
      
      Scalar c1_Df0{ m_c1 * Df0 };
      Scalar c2_Df0{ m_c2 * Df0 };

      // Initial trial
      Scalar alpha = alpha0;
      Scalar phi, der;
      eval( alpha, phi, der );

      // Check if immediately acceptable
      if ( phi <= f0 + alpha * c1_Df0 && der >= c2_Df0 ) return alpha;

      Scalar alpha_lo   = 0;
      Scalar phi_lo     = f0;
      Scalar der_lo     = Df0;

      Scalar alpha_hi   = 0;
      Scalar phi_hi     = 0;
      Scalar der_hi     = 0;

      bool   bracketed  = false;
      Scalar alpha_prev = 0;
      Scalar phi_prev   = f0;
      Scalar der_prev   = Df0;

      // === Bracketing Phase ===
      while ( evals < m_max_iters ) {
        if ( (phi > f0 + alpha*c1_Df0) || (evals > 1 && phi >= phi_prev) ) {
          alpha_lo  = alpha_prev;
          phi_lo    = phi_prev;
          der_lo    = der_prev;

          alpha_hi  = alpha;
          phi_hi    = phi;
          der_hi    = der;

          bracketed = true;
          break;
        }

        if ( der >= c2_Df0 ) return alpha;

        if ( der >= 0 ) {
          alpha_lo  = alpha_prev;
          phi_lo    = phi_prev;
          der_lo    = der_prev;

          alpha_hi  = alpha;
          phi_hi    = phi;
          der_hi    = der;

          bracketed = true;
          break;
        }

        Scalar new_alpha = min( 2 * alpha, m_alpha_max);
        alpha_prev = alpha;
        phi_prev   = phi;
        der_prev   = der;
        alpha      = new_alpha;
        eval(alpha, phi, der );

        if ( phi <= f0 + alpha * c1_Df0 && der >= c2_Df0 ) return alpha;
      }

      if (!bracketed) return std::nullopt;

      // === Zoom Phase ===
      while ( evals < m_max_iters ) {
        Scalar a_j  = abs(alpha_hi - alpha_lo) > m_epsi ?
                      detail::cubic_minimizer<Scalar>( alpha_lo, phi_lo, der_lo, alpha_hi, phi_hi, der_hi ) :
                      (alpha_lo + alpha_hi) / 2;

        Scalar phi_j, der_j;
        eval(a_j,phi_j, der_j);

        if ( (phi_j > f0 + a_j * c1_Df0) || (phi_j >= phi_lo) ) {
          alpha_hi = a_j;
          phi_hi   = phi_j;
          der_hi   = der_j;
        } else {
          if ( der_j >= c2_Df0 ) return a_j;

          if ( der_j * (alpha_hi - alpha_lo) >= 0 ) {
            alpha_hi = alpha_lo;
            phi_hi   = phi_lo;
            der_hi   = der_lo;
          }

          alpha_lo = a_j;
          phi_lo   = phi_j;
          der_lo   = der_j;
        }

        if ( abs(alpha_hi - alpha_lo) < m_epsi ) return a_j;
      }

      return std::nullopt;
    }
  };

  /**
   * @class StrongWolfeLineSearch
   * @brief Strong Wolfe conditions line search
   *
   * Implements line search satisfying the **strong Wolfe conditions**:
   * \f{align*}{
   *   f(x + \alpha p) &\leq f(x) + c_1 \alpha \nabla f(x)^T p \quad &\text{(Armijo)}\\
   *   |\nabla f(x + \alpha p)^T p| &\leq c_2 |\nabla f(x)^T p| \quad &\text{(Strong Curvature)}
   * \f}
   *
   * The key difference from regular Wolfe is the absolute value in the curvature
   * condition, which prevents the algorithm from accepting steps where the gradient
   * is still strongly negative (indicating the function is still decreasing rapidly).
   *
   * ### Why Strong Wolfe?
   *
   * The strong curvature condition ensures:
   * 1. Step is not too short (function still decreasing rapidly)
   * 2. Step is not too long (function starting to increase)
   * 3. Better guarantee of superlinear convergence
   * 4. More stable behavior near local minima
   *
   * ### Comparison with Regular Wolfe
   *
   * |  | Regular Wolfe | Strong Wolfe |
   * |--|--------------|--------------|
   * | Curvature | φ'(α) ≥ c₂ φ'(0) | \|φ'(α)\| ≤ c₂ \|φ'(0)\| |
   * | Accepts | α with φ'(α) still negative | Only near critical points |
   * | Convergence | Superlinear | Superlinear (better guaranteed) |
   * | Difficulty | Easier to satisfy | Slightly harder |
   *
   * ### When to Use
   *
   * Strong Wolfe is preferred for:
   * - L-BFGS optimization (most common choice)
   * - When high-quality steps are needed for fast convergence
   * - Problems with well-scaled variables
   * - When gradient evaluations are not too expensive
   *
   * @tparam Scalar Floating-point type
   *
   * @see Nocedal & Wright (2006), Section 3.1
   */
  template <typename Scalar>
  class StrongWolfeLineSearch {
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    Scalar m_c1        = 1e-4;   ///< Armijo parameter
    Scalar m_c2        = 0.9;    ///< Strong curvature parameter
    Scalar m_alpha_max = 10;     ///< Maximum step length
    Scalar m_epsi      = 1e-12;  ///< Interval convergence tolerance
    size_t m_max_iters = 50;     ///< Maximum iterations
    
    mutable Vector m_g_new, m_x_new;
    
  public:

    /**
     * @brief Perform Strong Wolfe line search
     *
     * Similar to regular Wolfe but with stricter curvature condition.
     *
     * @param f0 Function value at current point
     * @param Df0 Directional derivative (must be < 0)
     * @param x Current point
     * @param d Search direction
     * @param callback Function/gradient evaluator
     * @param alpha0 Initial step (default: 1.0)
     *
     * @return Step satisfying strong Wolfe, or std::nullopt if failed
     */
    template <typename Callback>
    std::optional<Scalar>
    operator()(
      Scalar           f0,
      Scalar           Df0,
      Vector   const & x,
      Vector   const & d,
      Callback const & callback,
      Scalar           alpha0 = 1
    ) const {

      if ( !(Df0 < 0) ) return std::nullopt;
      m_g_new.resize( x.size() );
      m_x_new.resize( x.size() );
      size_t evals = 0;

      auto eval = [&]( Scalar a, Scalar & f, Scalar & df ) {
        ++evals;
        m_x_new.noalias() = x + a * d;
        f  = callback( m_x_new, &m_g_new );
        df = m_g_new.dot(d);
      };
      
      Scalar c1_Df0{ m_c1 * Df0 };
      Scalar c2_Df0{ m_c2 * Df0 };

      // Initial trial
      Scalar alpha = alpha0;
      Scalar phi, der;
      eval( alpha, phi, der );

      // Check if acceptable (note: abs(der) for strong Wolfe)
      if ( phi <= f0 + alpha * c1_Df0 && abs(der) <= -c2_Df0 ) return alpha;

      Scalar alpha_lo   = 0;
      Scalar phi_lo     = f0;
      Scalar der_lo     = Df0;

      Scalar alpha_hi   = 0;
      Scalar phi_hi     = 0;
      Scalar der_hi     = 0;

      bool   bracketed  = false;
      Scalar alpha_prev = 0;
      Scalar phi_prev   = f0;
      Scalar der_prev   = Df0;

      // === Bracketing ===
      while ( evals < m_max_iters ) {

        if ( (phi > f0 + alpha*c1_Df0) || (evals > 1 && phi >= phi_prev) ) {
          alpha_lo  = alpha_prev;
          phi_lo    = phi_prev;
          der_lo    = der_prev;

          alpha_hi  = alpha;
          phi_hi    = phi;
          der_hi    = der;

          bracketed = true;
          break;
        }

        if ( abs(der) <= -c2_Df0 ) return alpha;

        if ( der >= 0 ) {
          alpha_lo  = alpha_prev;
          phi_lo    = phi_prev;
          der_lo    = der_prev;

          alpha_hi  = alpha;
          phi_hi    = phi;
          der_hi    = der;

          bracketed = true;
          break;
        }

        Scalar new_alpha = min( 2 * alpha, m_alpha_max );
        alpha_prev = alpha;
        phi_prev   = phi;
        der_prev   = der;
        alpha      = new_alpha;
        eval(alpha, phi, der );

        if ( phi <= f0 + alpha * c1_Df0 && abs(der) <= -c2_Df0 ) return alpha;
      }

      if (!bracketed) return std::nullopt;

      // === Zoom ===
      while ( evals < m_max_iters ) {
        Scalar a_j  = abs(alpha_hi - alpha_lo) > m_epsi ?
                      detail::cubic_minimizer<Scalar>( alpha_lo, phi_lo, der_lo, alpha_hi, phi_hi, der_hi ) :
                      (alpha_lo + alpha_hi) / 2;

        Scalar phi_j, der_j;
        eval(a_j,phi_j, der_j);

        if ( (phi_j > f0 + a_j * c1_Df0) || (phi_j >= phi_lo) ) {
          alpha_hi = a_j;
          phi_hi   = phi_j;
          der_hi   = der_j;
        } else {
          if ( abs(der_j) <= -c2_Df0 ) return a_j;

          if ( der_j * (alpha_hi - alpha_lo) >= 0 ) {
            alpha_hi = alpha_lo;
            phi_hi   = phi_lo;
            der_hi   = der_lo;
          }

          alpha_lo = a_j;
          phi_lo   = phi_j;
          der_lo   = der_j;
        }

        if ( abs(alpha_hi - alpha_lo) < m_epsi ) return a_j;
      }

      return std::nullopt;
    }
  };

  /**
   * @class GoldsteinLineSearch
   * @brief Goldstein conditions line search
   *
   * Implements Goldstein conditions, which require the step length to satisfy:
   * \f[
   *   f(x) + (1-c) \alpha \nabla f(x)^T p \leq f(x+\alpha p) \leq f(x) + c \alpha \nabla f(x)^T p
   * \f]
   * where \f$0 < c < 0.5\f$ (typically c = 0.1 to 0.3).
   *
   * ### Interpretation
   *
   * The Goldstein conditions bound the function decrease from both sides:
   * - **Upper bound** (Armijo-like): Ensures sufficient decrease
   * - **Lower bound**: Prevents steps that are too short
   *
   * This creates a "band" of acceptable function values at each step length.
   *
   * ### Comparison with Wolfe
   *
   * | Feature | Goldstein | Wolfe |
   * |---------|-----------|-------|
   * | Conditions | Two-sided function bounds | Function + derivative |
   * | Derivatives | Only at current point | At current and trial |
   * | Cost | Cheaper (no trial gradients) | More expensive |
   * | Theory | Less common in modern software | Standard choice |
   * | Convergence | Guarantees for convex functions | Stronger guarantees |
   *
   * ### Algorithm Strategy
   *
   * 1. Start with α = α₀
   * 2. Check if both Goldstein conditions are satisfied
   * 3. If f too large: reduce α (upper bound violated)
   * 4. If f too small: increase α (lower bound violated)
   * 5. Repeat until both conditions satisfied
   *
   * ### When to Use
   *
   * Goldstein is useful when:
   * - Gradient evaluations are expensive
   * - Problem is moderately well-conditioned
   * - Simple implementation is preferred
   * - Historical compatibility is needed
   *
   * However, Wolfe conditions are generally preferred in modern implementations
   * due to better theoretical properties and broader applicability.
   *
   * @tparam Scalar Floating-point type
   *
   * @note This is mainly included for completeness. For production use,
   *       prefer Strong Wolfe or More-Thuente line searches.
   *
   * @see Goldstein, A.A. (1965). "On Steepest Descent". SIAM Journal on Control.
   * @see Nocedal & Wright (2006), Section 3.1 for comparison with Wolfe
   */
  template <typename Scalar>
  class GoldsteinLineSearch {
  using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    Scalar m_c1          = 0.1;    // Parametro per Armijo (tipicamente 0.1-0.3)
    Scalar m_step_reduce = 0.5;    // Fattore di riduzione del passo
    Scalar m_step_expand = 1.2;    // Fattore di espansione del passo (più conservativo)
    size_t m_max_iters   = 50;

    mutable Vector m_x_new;

  public:

    template <typename Callback>
    std::optional<Scalar>
    operator()(
      Scalar           f0,
      Scalar           Df0,
      Vector   const & x,
      Vector   const & d,
      Callback const & callback,
      Scalar           alpha0 = 1
    ) const {
      // Controllo direzione di discesa
      if ( Df0 >= 0 ) return std::nullopt;

      Scalar alpha = alpha0;
      m_x_new.resize(x.size());
    
      // Calcola i bound di Goldstein
      Scalar armijo_bound    = f0 + m_c1 * alpha * Df0;
      Scalar goldstein_bound = f0 + (1 - m_c1) * alpha * Df0;
    
      for (size_t k = 0; k < m_max_iters; ++k) {
        m_x_new.noalias() = x + alpha * d;
        Scalar f_new = callback(m_x_new, nullptr); // Solo valore funzione, non gradiente
      
        // DEBUG: stampa per debugging
        // fmt::print("Goldstein iter {}: alpha={}, f_new={}, Armijo_bound={}, Goldstein_bound={}\n",
        //            k, alpha, f_new, armijo_bound, goldstein_bound);
      
        // Verifica condizioni di Goldstein
        bool satisfies_armijo    = (f_new <= armijo_bound);
        bool satisfies_goldstein = (f_new >= goldstein_bound);
      
        if ( satisfies_armijo && satisfies_goldstein ) {
          // fmt::print("Goldstein SUCCESS: alpha={} accepted\n", alpha);
          return alpha;
        }
      
        if ( !satisfies_armijo ) {
          // Passo troppo lungo - riduci
          alpha *= m_step_reduce;
        } else {
          // Passo troppo corto (soddisfa Armijo ma non Goldstein) - aumenta
          alpha *= m_step_expand;
        }
     
        // Ricalcola i bound con il nuovo alpha
        armijo_bound    = f0 + m_c1 * alpha * Df0;
        goldstein_bound = f0 + (1 - m_c1) * alpha * Df0;
     
        // Controllo per alpha troppo piccolo
        if ( alpha < std::numeric_limits<Scalar>::epsilon() * 10 ) {
          // fmt::print("Goldstein FAILED: step too small\n");
          return std::nullopt;
        }
      }
    
      // fmt::print("Goldstein FAILED: max iterations reached\n");
      return std::nullopt;
    }
  };
  
  /**
   * @brief Hager–Zhang inspired line search with enhancements.
   *
   * This class implements a robust line-search procedure inspired by the
   * Hager–Zhang algorithm. It enforces the Armijo (sufficient decrease)
   * condition and a (modified) Wolfe curvature condition with small
   * tolerances to avoid numerical issues. Zoom uses cubic interpolation
   * with explicit δ/σ safeguards; a bisection fallback is used whenever
   * interpolation yields an out-of-range value or NaN/inf.
   *
   * Notes:
   *  - The callback must have signature: Scalar callback(const Vector &x, Vector *g)
   *    returning f(x) and writing gradient into *g (if g != nullptr).
   *
   * @tparam Scalar float/double/etc.
   */
  template <typename Scalar>
  class HagerZhangLineSearch {
  public:
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
  
    /**
     * @brief Default constructor — parameters can be tuned after construction.
     *
     * @param c1 Armijo constant (typically 1e-4).
     * @param c2 Curvature constant (typically in (0.1, 0.9), e.g. 0.9).
     */
    HagerZhangLineSearch(
      Scalar c1 = Scalar(1e-4),
      Scalar c2 = Scalar(0.9)
    ) :
      m_c1(c1),
      m_c2(c2)
    {}
  
    /** @name Tunable parameters */
    Scalar m_c1        = Scalar(1e-4);   ///< Armijo parameter.
    Scalar m_c2        = Scalar(0.9);    ///< Curvature parameter.
    Scalar m_alpha_max = Scalar(10.0);   ///< Maximum allowed step length.
    Scalar m_epsi      = Scalar(1e-12);  ///< Absolute tolerance for interval collapse (used relatively).
    size_t m_max_iters = 50;             ///< Max outer iterations (extrapolation/zoom cycles).
    size_t m_max_evals = 200;            ///< Max function/gradient evaluations total.
  
    // Hager–Zhang specific interpolation safeguards (now used)
    Scalar m_delta     = Scalar(0.1);    ///< left-side safety fraction for interpolation
    Scalar m_sigma     = Scalar(0.9);    ///< right-side safety fraction for interpolation
    Scalar m_epsilon_k = Scalar(1e-6);   ///< tolerance used to relax curvature checks slightly
  
    // internal mutable state (so operator() can be const)
    mutable size_t  m_evals = 0;
    mutable Scalar m_f0;
    mutable Scalar m_Df0;
    mutable Scalar m_c1_Df0;
    mutable Scalar m_c2_Df0;
  
    mutable Vector m_g_new;
    mutable Vector m_x_new;
  
    /**
     * @brief Perform Hager–Zhang-style line search.
     *
     * @tparam Callback signature: Scalar callback(const Vector &x, Vector *g)
     *         where callback returns f(x) and writes gradient into *g (if non-null).
     *
     * @param f0        initial function value f(x).
     * @param Df0       initial directional derivative g(x)^T d (must be < 0).
     * @param x         current point.
     * @param d         search direction (descent direction).
     * @param callback  function+gradient evaluator.
     * @param alpha0    initial step (default 1).
     *
     * @return optional step length satisfying Armijo and curvature conditions, or std::nullopt.
     */
    template <typename Callback>
    std::optional<Scalar>
    operator()(
      Scalar           f0,
      Scalar           Df0,
      Vector   const & x,
      Vector   const & d,
      Callback const & callback,
      Scalar           alpha0 = Scalar(1)
    ) const {
      // require descent
      if (!(Df0 < Scalar(0))) return std::nullopt;
  
      // initialize internal state
      m_f0     = f0;
      m_Df0    = Df0;
      m_c1_Df0 = m_c1 * Df0;
      m_c2_Df0 = m_c2 * Df0;
  
      m_g_new.resize(x.size());
      m_x_new.resize(x.size());
      m_evals = 0;
  
      // evaluation wrapper (records count)
      auto eval = [&](Scalar a, Scalar &f, Scalar &df) {
        ++m_evals;
        m_x_new.noalias() = x + a * d;
        f = callback(m_x_new, &m_g_new);
        df = m_g_new.dot(d);
      };
  
      // initial values
      Scalar alpha      = min(alpha0, m_alpha_max);
      Scalar alpha_prev = Scalar(0);
      Scalar phi_prev   = f0;
      Scalar der_prev   = Df0;
  
      // Evaluate at initial alpha
      Scalar phi, der;
      eval(alpha, phi, der);
  
      // Minimum step threshold to avoid pointless tiny steps
      const Scalar alpha_min_threshold = std::numeric_limits<Scalar>::epsilon() * Scalar(100);
  
      for (size_t iter = 0; iter < m_max_iters && m_evals < m_max_evals; ++iter) {
  
        // Check Armijo (sufficient decrease) and curvature (Wolfe) conditions
        // Armijo: phi <= f0 + alpha * c1 * Df0
        // Curvature: der >= c2 * Df0
        if ( phi <= (m_f0 + alpha * m_c1_Df0) &&
             der >= (m_c2_Df0 - m_epsilon_k) ) {
          return alpha;
        }
  
        // Bracketing: phi is too large (doesn't satisfy Armijo) or phi is not decreasing
        if ( phi > (m_f0 + alpha * m_c1_Df0) ||
            (iter > 0 && phi >= phi_prev) ) {
          return zoom(alpha_prev, phi_prev, der_prev, alpha, phi, der, eval);
        }
  
        // If derivative is already small in magnitude and satisfies modified curvature:
        // Here we check |der| <= -c2*Df0 + epsilon  equivalently der >= c2*Df0 - eps
        if ( abs(der) <= (abs(m_c2_Df0) + m_epsilon_k) &&
             der >= (m_c2_Df0 - m_epsilon_k)) {
          // Accept alpha if derivative is close to satisfying curvature
          return alpha;
        }
  
        // If derivative non-negative, bracket and zoom
        if (der >= Scalar(0)) {
          return zoom(alpha, phi, der, alpha_prev, phi_prev, der_prev, eval);
        }
  
        // Extrapolation: increase alpha (doubling), but respect maximum allowed step
        Scalar alpha_new = min(Scalar(2) * alpha, m_alpha_max);
  
        // If step growth stalls or becomes too tiny, fail
        if (alpha_new <= alpha_min_threshold) return std::nullopt;
  
        alpha_prev = alpha;
        phi_prev   = phi;
        der_prev   = der;
  
        alpha = alpha_new;
  
        // If we've exhausted max evaluations, break
        if (m_evals >= m_max_evals) break;
  
        eval(alpha, phi, der);
      }
  
      return std::nullopt;
    }
  
  private:
    /**
     * @brief Zoom phase: refine [a_lo, a_hi] to find acceptable alpha.
     *
     * Uses cubic minimization (detail::cubic_minimizer) then enforces safety:
     * a_j ∈ [a_lo + δ*(a_hi-a_lo), a_hi - σ*(a_hi-a_lo)].
     * If cubic minimizer returns out of range / NaN / inf, fallback to bisection.
     *
     * @tparam EvalFunc callback type used to compute phi and derivative at a trial alpha.
     */
    template <typename EvalFunc>
    std::optional<Scalar>
    zoom(
      Scalar a_lo, Scalar phi_lo, Scalar der_lo,
      Scalar a_hi, Scalar phi_hi, Scalar der_hi,
      EvalFunc const & eval
    ) const {
      // relative collapse tolerance
      const auto rel_tol = [&](Scalar a, Scalar b) {
        return abs(a - b) <= m_epsi * (Scalar(1) + abs(a));
      };
  
      for (size_t iter = 0; iter < m_max_iters && m_evals < m_max_evals; ++iter) {
        // Attempt cubic minimization using existing helper (assumed available).
        Scalar a_j = detail::cubic_minimizer<Scalar>(a_lo, phi_lo, der_lo, a_hi, phi_hi, der_hi);
  
        // define safeguarded interval for a_j
        Scalar lo_guard = a_lo + m_delta * (a_hi - a_lo);
        Scalar hi_guard = a_hi - m_sigma * (a_hi - a_lo);
  
        // Fallback to bisection if cubic gives bad value or outside guard interval or non-finite
        if (!(a_j > lo_guard && a_j < hi_guard) || !std::isfinite(a_j)) {
          a_j = (a_lo + a_hi) / 2;
        }
  
        Scalar phi_j, der_j;
        eval(a_j, phi_j, der_j);
  
        // Check Armijo at a_j
        if (phi_j > (m_f0 + a_j * m_c1_Df0) || phi_j >= phi_lo) {
          // a_j is too big — tighten upper bound
          a_hi   = a_j;
          phi_hi = phi_j;
          der_hi = der_j;
        } else {
          // sufficient decrease holds at a_j — check curvature
          if (der_j >= (m_c2_Df0 - m_epsilon_k)) {
            // derivative satisfies curvature condition -> accept a_j
            return a_j;
          }
  
          // If derivative has same sign as (a_hi - a_lo), move hi to lo (re-orient interval)
          if (der_j * (a_hi - a_lo) >= Scalar(0)) {
            a_hi   = a_lo;
            phi_hi = phi_lo;
            der_hi = der_lo;
          }
  
          // move lower bound up
          a_lo   = a_j;
          phi_lo = phi_j;
          der_lo = der_j;
        }
  
        // If bracket collapsed sufficiently (relative check), return best point so far
        if (rel_tol(a_hi, a_lo)) {
          // choose midpoint as a candidate
          Scalar a_mid = (a_lo + a_hi) / 2;
          // evaluate midpoint if we still have eval budget
          if (m_evals < m_max_evals) {
            Scalar phi_mid, der_mid;
            eval(a_mid, phi_mid, der_mid);
            // accept midpoint if curvature satisfied
            if (phi_mid <= (m_f0 + a_mid * m_c1_Df0) && der_mid >= (m_c2_Df0 - m_epsilon_k)) {
              return a_mid;
            }
          }
          // otherwise return the last tested a_j as best-effort
          return a_j;
        }
      }
  
      return std::nullopt;
    }
  };

  /**
   * @class MoreThuenteLineSearch
   * @brief Implements the Moré–Thuente line search algorithm for Strong Wolfe conditions.
   *
   * This line search algorithm attempts to find a step length α along a descent
   * direction d that satisfies:
   *
   * 1. **Armijo condition (sufficient decrease)**
   *    \f[ φ(α) ≤ φ(0) + c_1 α φ'(0) \f]
   *
   * 2. **Curvature condition (Strong Wolfe)**
   *    \f[ |φ'(α)| ≤ c_2 |φ'(0)| \f]
   *
   * The algorithm works in two phases:
   * - **Extrapolation phase:** tries to find an interval where the conditions may hold.
   * - **Zoom phase:** refines the interval using cubic interpolation until a valid α is found.
   *
   * ### Algorithm Outline
   *
   * 1. Evaluate initial trial α0.
   * 2. If Strong Wolfe conditions are satisfied → return α0.
   * 3. Otherwise:
   *    - Update bracketing interval [α_lo, α_hi].
   *    - If bracketed → enter zoom phase (iterative cubic interpolation).
   *    - If not bracketed → extrapolate α (doubling or cubic extrapolation).
   * 4. Continue until:
   *    - A step satisfying Strong Wolfe is found, or
   *    - Maximum iterations / function evaluations reached.
   *
   * ### Notes / Improvements
   *
   * - Alpha_min should be enforced to avoid numerical underflow.
   * - Zoom phase should be separated for clarity and robustness.
   * - Relative tolerance for bracket collapse is recommended.
   * - Cubic extrapolation for unbracketed phase can accelerate convergence in steep directions.
   * - Robust derivative checks (using absolute values) are safer than signed comparisons.
   *
   * @tparam Scalar floating point type (float, double, long double)
   */

  template <typename Scalar>
  class MoreThuenteLineSearch {
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    Scalar m_c1        = 1e-4; // Armijo
    Scalar m_c2        = 0.9;  // Strong Wolfe (Curvature)
    Scalar m_alpha_max = 10.0;
    Scalar m_alpha_min = 1e-12; // Passo minimo accettabile
    Scalar m_epsi      = 1e-12;
    size_t m_max_iters = 50;
    
    mutable Vector m_g_new, m_x_new;

    // --- Stati interni (per chiarezza, gestiti tramite i parametri passati) ---
    // Questi potrebbero essere membri privati, ma li gestiremo nel main operator() per isolamento.

  public:

    template <typename Callback>
    std::optional<Scalar>
    operator()(
      Scalar           f0,
      Scalar           Df0,
      Vector   const & x,
      Vector   const & d,
      Callback const & callback,
      Scalar           alpha0 = 1
    ) const {
      if ( !(Df0 < 0) ) return std::nullopt;

      m_g_new.resize( x.size() );
      m_x_new.resize( x.size() );
      size_t evals = 0;

      auto eval = [&]( Scalar a, Scalar & f, Scalar & df ) {
        ++evals;
        m_x_new.noalias() = x + a * d;
        f  = callback( m_x_new, &m_g_new );
        df = m_g_new.dot(d);
      };
        
      Scalar c1_Df0{ m_c1 * Df0 };
      Scalar c2_Df0{ m_c2 * Df0 };

      // ----------------------------------------------------
      // FASE 1: Inizializzazione e ciclo principale (Bracketing & Zoom)
      // ----------------------------------------------------

      // Punti di lavoro
      Scalar alpha_prev = 0;
      Scalar phi_prev   = f0;
      Scalar der_prev   = Df0;
        
      Scalar alpha_curr = alpha0; // Passo di prova corrente
      Scalar phi_curr, der_curr;
        
      // Intervallo di bracketing iniziale
      Scalar alpha_lo = 0;
      Scalar phi_lo   = f0;
      Scalar der_lo   = Df0;
        
      Scalar alpha_hi = m_alpha_max;
        
      bool bracketed = false;

      // Esegui la prima valutazione
      eval(alpha_curr, phi_curr, der_curr);
        
      for ( size_t k{0}; k < m_max_iters && evals < m_max_iters; ++k ) {
            
        // Check Strong Wolfe al passo corrente (MT è un raffinamento per SW)
        if ( phi_curr <= f0 + alpha_curr * c1_Df0 && abs(der_curr) <= -c2_Df0 ) return alpha_curr; // Trovato passo accettabile

        // --- Bracketing Logic (Aggiornamento dell'intervallo [alpha_lo, alpha_hi]) ---

        if ( phi_curr > f0 + alpha_curr * c1_Df0 || (bracketed && phi_curr >= phi_lo) ) {
          // Caso A: Violata Armijo o funzione non decrescente
          // L'intervallo [alpha_lo, alpha_curr] contiene il minimo.
          alpha_hi = alpha_curr;
          // Mantieni alpha_lo (potrebbe essere alpha_prev o 0)
          bracketed = true;
        } else if ( abs(der_curr) <= -c2_Df0 ) {
          // Caso B: Soddisfatta Armijo ma derivata non abbastanza piatta (Strong Wolfe fallita)
          // Se la derivata è troppo negativa (troppo ripida), dobbiamo andare più avanti.
          // Se der_curr < 0, il minimo è oltre alpha_curr.
          alpha_lo = alpha_curr;
          phi_lo   = phi_curr;
          der_lo   = der_curr;
        } else if ( der_curr >= 0 ) {
          // Caso C: Trovata pendenza positiva, l'intervallo è [alpha_prev, alpha_curr]
          alpha_hi = alpha_curr;
          // Mantieni alpha_lo
          bracketed = true;
        }

        // --- Selezione del nuovo passo ---
        Scalar alpha_new;

        if ( bracketed ) {
          // Zoom phase: usa interpolazione cubica/safeguard tra alpha_lo e alpha_hi
          alpha_new = detail::compute_step(
            alpha_lo,   phi_lo,   der_lo,
            alpha_hi,   phi_curr, der_curr, // Nota: usa phi_curr e der_curr per alpha_hi (il passo "cattivo")
            alpha_prev, phi_prev, der_prev,
            m_alpha_max, m_alpha_max,
            false // is_bracketing = false
          );
        } else {
          // Extrapolation phase: usa interpolazione o raddoppia il passo
                
          // Opzione 1: Raddoppia
          alpha_new = 2 * alpha_curr;

          // Opzione 2: Interpolazione/Estrapolazione (più avanzata, ma rischiosa)
          // alpha_new = detail::compute_step(
          //     alpha_lo, phi_lo, der_lo,
          //     alpha_curr, phi_curr, der_curr,
          //     alpha_prev, phi_prev, der_prev,
          //     m_alpha_max, m_alpha_max,
          //     true // is_bracketing = true
          // );
        }

        // --- Salvaguardia finale ---
        if ( alpha_new <= 0 || alpha_new >= m_alpha_max ) alpha_new = min( 2 * alpha_curr, m_alpha_max);
        if ( abs(alpha_hi - alpha_lo) < m_epsi     ) return alpha_curr; // Convergenza dell'intervallo

        // Aggiorna stato per la prossima iterazione
        alpha_prev = alpha_curr;
        phi_prev   = phi_curr;
        der_prev   = der_curr;

        alpha_curr = alpha_new;
        eval(alpha_curr, phi_curr, der_curr);
      }

      return std::nullopt; // Raggiunto il massimo delle iterazioni
    }
  };

  // ---------------------------------------------------------------------------
  // LBFGSMinimizer: high-level optimizer with line-search and optional box bounds
  // ---------------------------------------------------------------------------

  template <typename Scalar = double>
  class LBFGS_minimizer {
  public:
    using Vector   = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Callback = std::function<Scalar(Vector const &, Vector *)>;

    enum class Status {
      CONVERGED          = 0, // Convergenza raggiunta
      MAX_ITERATIONS     = 1, // Massimo numero di iterazioni raggiunto
      LINE_SEARCH_FAILED = 2, // Line search fallita
      GRADIENT_TOO_SMALL = 3, // Gradiente troppo piccolo
      FAILED             = 4  // Fallimento generico
    };

    // AGGIUNTE: Strutture per raccogliere i risultati
    struct IterationData {
      Status status{Status::FAILED};
      size_t iterations{0};
      Scalar final_gradient_norm{0};
      Scalar final_function_value{0};
      Scalar initial_function_value{0};
      size_t function_evaluations{0};
      size_t gradient_evaluations{0};
      size_t line_search_iterations{0};
    };

    struct Options {
      size_t max_iter        { 200   };
      size_t iter_reset      { 50    };
      size_t m               { 20    };

      Scalar g_tol           { 1e-8  };
      Scalar f_tol           { 1e-12 };
      Scalar x_tol           { 1e-10 };
      
      Scalar step_max        { 1     };
      Scalar sty_min_factor  { 1e-12  }; // Più permissivo
      Scalar very_small_step { 1e-8  };

      bool   use_projection  { false };
      bool   verbose         { false };
    };

  private:

    Scalar        m_epsi{ 1e-15 };
    Options       m_options;
    LBFGS<Scalar> m_LBFGS;
    Vector        m_lower;
    Vector        m_upper;

    // project x into bounds
    void
    project_inplace( Vector & x ) const
    { x = x.cwiseMax(m_lower).cwiseMin(m_upper); }

    Vector
    projected_gradient( Vector const & x, Vector const & g ) const {
      auto freeze = ((x.array() <= m_lower.array()) && (g.array() >= 0)) ||
                    ((x.array() >= m_upper.array()) && (g.array() <= 0));
      return freeze.select(Vector::Zero(g.size()), g);
    }

    void
    projected_gradient_inplace( Vector const & x, Vector & g ) const {
      auto freeze = ((x.array() <= m_lower.array()) && (g.array() >= 0)) ||
                    ((x.array() >= m_upper.array()) && (g.array() <= 0));
      g = freeze.select(Vector::Zero(g.size()), g);
    }

    void
    projected_direction_inplace( Vector const & x, Vector & d ) const {
      auto freeze = ((x.array() <= m_lower.array()) && (d.array() < 0)) ||
                    ((x.array() >= m_upper.array()) && (d.array() > 0));
      d = freeze.select(Vector::Zero(d.size()), d);
    }

    Scalar
    projected_gradient_norm( Vector const & x, Vector const & g ) const {
      auto freeze = ((x.array() <= m_lower.array()) && (g.array() >= 0)) ||
                    ((x.array() >= m_upper.array()) && (g.array() <= 0));
      return freeze.select(Vector::Zero(g.size()), g)
                   .array()
                   .abs()
                   .maxCoeff();
    }

    mutable Vector m_x, m_g, m_p, m_x_new, m_g_new, m_s, m_y;
    mutable size_t m_iter_since_reset{0};

    // AGGIUNTE: Contatori per statistiche
    mutable size_t m_function_evaluations{0};
    mutable size_t m_gradient_evaluations{0};
    mutable size_t m_line_search_iterations{0};

  public:

    LBFGS_minimizer( Options opts = Options() )
    : m_options(opts)
    , m_LBFGS(opts.m)
    {}
    
    Vector const & solution() const { return m_x; }

    void
    set_bounds( Vector const & lower, Vector const & upper ) {
      assert( lower.size() == upper.size() );
      m_lower = lower;
      m_upper = upper;
      m_options.use_projection = true;
    }

    void
    set_bounds( size_t n, Scalar const lower[], Scalar const upper[] ) {
      m_lower.resize( n );
      m_upper.resize( n );
      std::copy_n( lower, n, m_lower.data() );
      std::copy_n( upper, n, m_upper.data() );
      m_options.use_projection = true;
    }

    void
    reset_memory() {
      if ( m_options.verbose ) fmt::print("[LBFGS] Periodic memory reset\n");
      m_LBFGS.clear();
      m_iter_since_reset = 0;
    }

    template <typename Linesearch>
    IterationData
    minimize(
      Vector     const & x0,
      Callback   const & callback,
      Linesearch const & linesearch = MoreThuenteLineSearch<Scalar>()
    ) {
      
      Status status { Status::MAX_ITERATIONS };
      Scalar gnorm  { 0 };

      // AGGIUNTA: Reset contatori
      m_function_evaluations = 0;
      m_gradient_evaluations = 0;
      m_line_search_iterations = 0;

      auto const n{ x0.size() };
      m_x.resize(n);     m_g.resize(n);
      m_x_new.resize(n); m_g_new.resize(n);
      m_p.resize(n);     m_s.resize(n); m_y.resize(n);
      
      if ( m_options.use_projection ) {
        assert( m_lower.size() == n );
        assert( m_upper.size() == n );
      }

      m_x.noalias() = x0;
      if ( m_options.use_projection ) project_inplace(m_x);
      
      Scalar f = callback(m_x, &m_g);
      m_function_evaluations++;
      m_gradient_evaluations++;
      
      Scalar f_prev{ f };
      Vector x_prev    = m_x;
      Scalar step_prev = 0;

      // AGGIUNTA: Salva valore iniziale
      Scalar f_initial = f;
      size_t iteration{ 0 };
      bool converged = false;

      for (; iteration < m_options.max_iter && !converged; ++iteration ) {
        ++m_iter_since_reset;
        
        // Check per stagnazione
        if ( iteration > 0 ) {
          Scalar x_change = (m_x - x_prev).norm();
          if ( x_change < step_prev*m_options.x_tol ) {
            if ( m_options.verbose )
              fmt::print("[LBFGS] Converged by x change: {:.2e} < {:.2e}\n", x_change, m_options.x_tol);
            gnorm  = projected_gradient_norm(m_x, m_g);
            status = Status::CONVERGED;
            converged = true;
            break;
          }
        }
        x_prev = m_x;

        gnorm = projected_gradient_norm(m_x, m_g);
        if ( m_options.verbose )
          fmt::print("[LBFGS] iter={:3d} f={:12.4g} ‖pg‖={:12.4g} mem={}\n", iteration, f, gnorm, m_LBFGS.size());
        
        if ( gnorm <= m_options.g_tol ) {
          status = Status::GRADIENT_TOO_SMALL;
          converged = true;
          break;
        }

        // Reset periodico della memoria ogni iter_reset iterazioni
        if ( m_iter_since_reset > m_options.iter_reset ) reset_memory();

        // compute search direction
        Scalar h0 = m_LBFGS.compute_initial_h0(1);
        m_p = -m_LBFGS.two_loop_recursion(m_g, h0);

        if ( m_options.use_projection ) {
          projected_direction_inplace( m_x, m_p );
          if ( m_p.isZero() ) m_p = -projected_gradient(m_x, m_g);
        }

        // Robust descent direction check
        Scalar pg             = m_p.dot(m_g);
        size_t fallback_count = 0;
        size_t max_fallback   = 3;

        while ((!std::isfinite(pg) || pg >= -m_epsi) && fallback_count < max_fallback) {
          if ( fallback_count == 0 ) {
            m_p = -projected_gradient(m_x, m_g);
          } else if ( fallback_count == 1 ) {
            m_LBFGS.clear();
            h0  = m_LBFGS.compute_initial_h0(1);
            m_p = -m_LBFGS.two_loop_recursion(m_g, h0);
            if ( m_options.use_projection ) projected_direction_inplace(m_x, m_p);
          } else {
            m_p = -m_g / (1 + m_g.norm());
            if ( m_options.use_projection ) projected_direction_inplace(m_x, m_p);
          }
          pg = m_p.dot(m_g);
          fallback_count++;
        }

        if (!std::isfinite(pg) || pg >= -m_epsi) {
          if ( m_options.verbose )
            fmt::print("[LBFGS] Cannot find descent direction, stopping\n");
          gnorm  = projected_gradient_norm(m_x, m_g);
          status = Status::FAILED;
          converged = true;
          break;
        }

        // Line search
        auto step_opt = linesearch( f, pg, m_x, m_p, callback, m_options.step_max );
        
        if (!step_opt.has_value()) {
          if ( m_options.verbose )
            fmt::print("[LBFGS] line search failed, resetting memory\n");
          m_LBFGS.clear();
          // Prova con passo fisso piccolo
          m_x_new.noalias() = m_x + m_options.very_small_step * m_p;
          if ( m_options.use_projection ) project_inplace(m_x_new);
          Scalar f_test = callback(m_x_new, &m_g_new);
          m_function_evaluations++;
          m_gradient_evaluations++;
          
          if (f_test < f) {
            // Accetta comunque il passo
            m_s.noalias() = m_x_new - m_x;
            m_y.noalias() = m_g_new - m_g;
            m_x.swap(m_x_new);
            m_g.swap(m_g_new);
            f = f_test;
            continue;
          } else {
            gnorm  = projected_gradient_norm(m_x, m_g);
            status = Status::LINE_SEARCH_FAILED;
            converged = true;
            break;
          }
        }

        Scalar step = *step_opt; step_prev = step > 1 ? step : 1 ;
        
        // evaluate final new point
        m_x_new.noalias() = m_x + step * m_p;
        if ( m_options.use_projection ) project_inplace(m_x_new);
        Scalar f_new = callback(m_x_new, &m_g_new);
        m_function_evaluations++;
        m_gradient_evaluations++;

        // compute s,y
        m_s.noalias() = m_x_new - m_x;
        m_y.noalias() = m_g_new - m_g;

        // Robust curvature check
        Scalar sty           = m_s.dot(m_y);
        Scalar sty_tolerance = max( m_options.sty_min_factor, m_epsi * m_s.squaredNorm() * m_y.squaredNorm());

        if (sty > sty_tolerance) {
          m_LBFGS.add_correction(m_s, m_y);
        } else {
          if ( m_options.verbose && m_LBFGS.size() > 0 ) fmt::print("[LBFGS] curvature condition failed (s^T y = {:.2e}), skipping update\n", sty);
          // Non cancellare la memoria, procedi senza aggiornare
        }

        // move
        m_x.swap(m_x_new);
        m_g.swap(m_g_new);
        f = f_new;

        // check function change
        Scalar f_change = abs(f - f_prev);
        if ( f_change <= m_options.f_tol ) {
          if ( m_options.verbose ) fmt::print("[LBFGS] Converged by function change: {:.2e} < {:.2e}\n", f_change, m_options.f_tol);
          gnorm  = projected_gradient_norm(m_x, m_g);
          status = Status::CONVERGED;
          converged = true;
          break;
        }
        f_prev = f;
      }
      
      if (!converged) {
        gnorm = projected_gradient_norm(m_x, m_g);
      }

      return IterationData{
        status,
        iteration,
        gnorm,
        f,
        f_initial,
        m_function_evaluations,
        m_gradient_evaluations,
        m_line_search_iterations
      };
    }
  };
}

#endif

#endif

//
// eof: Utils_LBFGS.hh
