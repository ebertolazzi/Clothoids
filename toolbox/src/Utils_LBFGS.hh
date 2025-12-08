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

#include <set>

#include "Utils_fmt.hh"
#include "Utils_nonlinear_linesearch.hh"

namespace Utils
{
  using std::abs;
  using std::max;
  using std::min;
  using std::set;

  namespace LBFGS_utils
  {

    // ===========================================================================
    // COLOR DEFINITIONS (consistent with NelderMead)
    // ===========================================================================

    namespace PrintColors
    {
      constexpr auto HEADER    = fmt::fg( fmt::color::light_blue );
      constexpr auto SUCCESS   = fmt::fg( fmt::color::green );
      constexpr auto WARNING   = fmt::fg( fmt::color::yellow );
      constexpr auto ERROR     = fmt::fg( fmt::color::red );
      constexpr auto INFO      = fmt::fg( fmt::color::cyan );
      constexpr auto ITERATION = fmt::fg( fmt::color::white );
      constexpr auto DETAIL    = fmt::fg( fmt::color::gray );
    }  // namespace PrintColors

    // ===========================================================================
    // UTILITY FUNCTIONS
    // ===========================================================================

    /**
     * @brief Format a vector of indices in a compact representation
     *
     * @tparam T Index type (typically size_t or int)
     * @param indices Vector of indices to format
     * @param max_display Maximum number of indices to display before truncating
     * @return std::string Compact string representation
     */
    template <typename T>
    std::string
    format_index_vector_compact( std::vector<T> const & indices, size_t max_display = 5 )
    {
      if ( indices.empty() ) return "[]";

      std::stringstream ss;
      ss << "[";

      if ( indices.size() <= max_display )
      {
        for ( size_t i = 0; i < indices.size(); ++i )
        {
          if ( i > 0 ) ss << ", ";
          ss << indices[i];
        }
      }
      else
      {
        for ( size_t i = 0; i < max_display - 1; ++i )
        {
          if ( i > 0 ) ss << ", ";
          ss << indices[i];
        }
        ss << ", ..., " << indices.back();
      }

      ss << "]";
      return ss.str();
    }

  }  // namespace LBFGS_utils

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
   * where \f$H_0 = \gamma I\f$ with \f$\gamma = \frac{s_{k-1}^T
   * y_{k-1}}{y_{k-1}^T y_{k-1}}\f$
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
  class LBFGS
  {
  public:
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

  private:
    size_t m_capacity;           ///< Maximum number of correction pairs (m)
    size_t m_dimension;          ///< Problem dimension (n)
    Matrix m_S;                  ///< Storage for s vectors (n × m matrix, each column is s_i)
    Matrix m_Y;                  ///< Storage for y vectors (n × m matrix, each column is y_i)
    Vector m_rho;                ///< Storage for ρ_i = 1/(y_i^T s_i) values
    size_t m_current_size{ 0 };  ///< Current number of stored pairs
    size_t m_oldest_index{ 0 };  ///< Circular buffer index of oldest pair
    size_t m_newest_index{ 0 };  ///< Circular buffer index where next pair will be stored

    /* Robustness parameters */
    bool   m_enable_damping{ true };    ///< Enable Powell-style damping by default
    Scalar m_h0_min{ Scalar( 1e-6 ) };  ///< Minimum allowed initial H0 scaling
    Scalar m_h0_max{ Scalar( 1e6 ) };   ///< Maximum allowed initial H0 scaling

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
     * @param maxCorrections Maximum number of stored correction pairs (typical:
     * 5-20). Larger values provide better Hessian approximation but increase
     * memory usage and computational cost.
     * @param problem_dimension Dimension of the optimization problem (n)
     *
     * @note A typical choice is m=10 for most problems. Use smaller m (3-7) for
     *       very large problems or when memory is constrained.
     */
    explicit LBFGS( size_t maxCorrections = 10, size_t problem_dimension = 0 )
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
    clear()
    {
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
    size_t
    size() const
    {
      return m_current_size;
    }

    /**
     * @brief Return maximum capacity for correction pairs
     * @return Maximum number of pairs (m parameter)
     */
    size_t
    capacity() const
    {
      return m_capacity;
    }

    /**
     * @brief Return problem dimension
     * @return Dimension of vectors in optimization problem (n)
     */
    size_t
    dimension() const
    {
      return m_dimension;
    }

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
    resize( size_t const new_dimension )
    {
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
     * Attempts to add a new correction pair after validating the curvature
     * condition. The pair is accepted only if:
     * \f[
     *   s^T y > \max\{\epsilon_{\text{rel}} \|s\| \|y\|, \epsilon_{\text{abs}}
     * \|s\|^2\}
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
    add_correction( Vector const & s, Vector const & y, Scalar const min_curvature_ratio = 1e-8 )
    {
      assert( s.size() == y.size() );
      if ( static_cast<size_t>( s.size() ) != m_dimension ) resize( s.size() );

      Scalar const sty{ s.dot( y ) };
      Scalar const snrm{ s.norm() };
      Scalar const ynrm{ y.norm() };

      // basic positivity check
      if ( !( sty > 0 ) ) return false;

      // relative and absolute thresholds for curvature test (robust)
      Scalar const eps        = std::numeric_limits<Scalar>::epsilon();
      Scalar const rel_thresh = min_curvature_ratio * snrm * ynrm;
      Scalar const abs_thresh = eps * snrm * snrm;
      Scalar const thresh     = max( rel_thresh, abs_thresh );

      if ( sty <= thresh ) return false;

      // Store the new pair
      m_S.col( m_newest_index ) = s;
      m_Y.col( m_newest_index ) = y;
      m_rho( m_newest_index )   = Scalar( 1.0 ) / sty;

      // Update indices and size
      if ( m_current_size < m_capacity )
      {
        // Buffer not full yet
        ++m_current_size;
        m_newest_index = ( m_newest_index + 1 ) % m_capacity;
      }
      else
      {
        // Buffer full - overwrite oldest
        m_oldest_index = ( m_oldest_index + 1 ) % m_capacity;
        m_newest_index = ( m_newest_index + 1 ) % m_capacity;
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
     * @note If no pairs are stored (size() == 0), returns h0 * g (scaled
     * steepest descent)
     *
     * @see Nocedal (1980) for original two-loop recursion algorithm
     * @see Liu & Nocedal (1989) for practical implementation details
     */
    Vector
    two_loop_recursion( Vector const & g, Scalar h0 ) const
    {
      if ( m_current_size == 0 ) return h0 * g;

      m_q = g;  // Start with q = g

      // =====================================================================
      // First loop: Process pairs from newest to oldest
      // =====================================================================
      size_t idx = ( m_newest_index == 0 ) ? m_capacity - 1 : m_newest_index - 1;

      size_t i{ m_current_size };
      while ( i > 0 )
      {
        --i;
        m_alpha( i ) = m_rho( idx ) * m_S.col( idx ).dot( m_q );
        m_q -= m_alpha( i ) * m_Y.col( idx );

        // Move to previous index in circular buffer
        idx = ( idx == 0 ) ? m_capacity - 1 : idx - 1;
      }

      // =====================================================================
      // Apply initial Hessian approximation H0 = h0 * I
      // =====================================================================
      m_r = h0 * m_q;

      // =====================================================================
      // Second loop: Process pairs from oldest to newest
      // =====================================================================
      idx = m_oldest_index;
      for ( size_t i{ 0 }; i < m_current_size; ++i )
      {
        Scalar beta = m_rho( idx ) * m_Y.col( idx ).dot( m_r );
        m_r += m_S.col( idx ) * ( m_alpha( i ) - beta );

        // Move to next index in circular buffer
        idx = ( idx + 1 ) % m_capacity;
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
     * - Self-scaling property: approximates the curvature along the most recent
     * step
     * - Theoretical foundation: relates to Barzilai-Borwein spectral steplength
     * - Empirical success: works well in practice for many problems
     *
     * The result is clamped to [m_h0_min, m_h0_max] to prevent extreme scaling
     * that could lead to numerical instability.
     *
     * @param default_value Value to return if no pairs are stored
     * (default: 1.0)
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
    compute_initial_h0( Scalar default_value = Scalar( 1.0 ) ) const
    {
      if ( m_current_size == 0 ) return default_value;

      // Get the newest pair (just before newest_index in circular buffer)
      size_t       latest_idx = ( m_newest_index == 0 ) ? m_capacity - 1 : m_newest_index - 1;
      auto const & s          = m_S.col( latest_idx );
      auto const & y          = m_Y.col( latest_idx );

      Scalar sty{ s.dot( y ) };
      Scalar yty{ y.dot( y ) };
      Scalar h0 = default_value;
      if ( yty > 0 ) h0 = sty / yty;

      // Clamp to safe range
      h0 = max( m_h0_min, min( m_h0_max, h0 ) );
      return h0;
    }

    /**
     * @brief Add correction with Powell damping for improved robustness
     *
     * Powell damping (also called modified BFGS) adjusts the gradient
     * difference y to ensure the curvature condition is satisfied, even when
     * the standard pair (s,y) would be rejected.
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
    add_correction_with_damping( LBFGS<Scalar> & lb,
                                 Vector const &  s,
                                 Vector const &  y,
                                 Scalar const    min_curvature_ratio = 1e-8 )
    {
      using Vector            = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
      Scalar const sty        = s.dot( y );
      Scalar const snrm       = s.norm();
      Scalar const ynrm       = y.norm();
      Scalar const eps        = std::numeric_limits<Scalar>::epsilon();
      Scalar const rel_thresh = min_curvature_ratio * snrm * ynrm;
      Scalar const abs_thresh = eps * snrm * snrm;
      Scalar const thresh     = max( rel_thresh, abs_thresh );

      if ( sty > thresh ) return lb.add_correction( s, y, min_curvature_ratio );

      // compute h0 and fallback if invalid
      Scalar h0 = lb.compute_initial_h0( 1 );
      if ( !( h0 > 0 ) ) return false;

      Scalar sBs = s.squaredNorm() / h0;  // s^T (B0 s) with B0 = H0^{-1} = (1/h0) I

      // If denominator would be zero or very small, don't attempt damping
      Scalar denom = sty - sBs;
      Vector y_hat( y );
      if ( abs( denom ) < std::numeric_limits<Scalar>::epsilon() )
      {
        // ambiguous, don't damp
        return lb.add_correction( s, y, min_curvature_ratio );
      }
      else
      {
        // choose theta so that s^T y_hat = thresh (solve theta*(sty - sBs) +
        // sBs = thresh)
        Scalar theta   = std::clamp( ( thresh - sBs ) / denom, 0, 1 );
        y_hat          = theta * y + ( 1 - theta ) * ( s / h0 );
        Scalar sty_hat = s.dot( y_hat );
        if ( !( sty_hat > 0 ) ) return false;
        return lb.add_correction( s, y_hat, min_curvature_ratio );
      }
    }
  };

}  // namespace Utils

#endif

#endif

//
// eof: Utils_LBFGS.hh
//
