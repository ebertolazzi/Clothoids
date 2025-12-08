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
// file: Utils_MADS.hh
//
// Implementation of Mesh Adaptive Direct Search (MADS) for derivative-free
// optimization with bound constraints.
//

#pragma once

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifndef UTILS_MADS_dot_HH
#define UTILS_MADS_dot_HH

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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <optional>
#include <random>
#include <set>
#include <utility>
#include <vector>

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

  /**
   * @brief Comparator for Eigen vectors to use in std::set
   */
  template <typename Scalar>
  struct VectorComparator
  {
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    bool
    operator()( Vector const & a, Vector const & b ) const
    {
      if ( a.size() != b.size() ) return a.size() < b.size();
      for ( int i = 0; i < a.size(); ++i )
      {
        if ( a( i ) != b( i ) ) return a( i ) < b( i );
      }
      return false;
    }
  };

  /**
   * @class MADS_minimizer
   * @brief Minimizer using Mesh Adaptive Direct Search (MADS)
   */
  template <typename Scalar = double>
  class MADS_minimizer
  {
  public:
    using Vector   = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;               ///< Vector type
    using Matrix   = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;  ///< Matrix type
    using Callback = std::function<Scalar( Vector const & )>;                ///< Function callback type

    /**
     * @struct Options
     * @brief Options for MADS minimizer
     */
    struct Options
    {
      size_t max_iter{ 1000 };          ///< Maximum iterations
      size_t max_evaluations{ 10000 };  ///< Maximum function evaluations
      Scalar initial_mesh_size{ 1.0 };  ///< Initial mesh size
      Scalar min_mesh_size{ 1e-6 };     ///< Minimum mesh size
      Scalar mesh_contraction{ 0.5 };   ///< Mesh contraction factor
      Scalar mesh_expansion{ 2.0 };     ///< Mesh expansion factor
      bool   verbose{ false };          ///< Verbose output
      size_t print_every{ 10 };         ///< Print every n iterations

      // Search and poll parameters
      bool        enable_search{ true };  ///< Enable search step
      size_t      search_points{ 100 };   ///< Number of search points
      size_t      poll_directions{ 2 };   ///< Number of poll directions per dimension
      std::string poll_method{ "ltm" };   ///< Poll method: "ltm" (LTMADS) or "ortho" (OrthoMADS)

      // Convergence criteria
      Scalar tol_mesh{ 1e-8 };  ///< Mesh size tolerance
      Scalar tol_f{ 1e-8 };     ///< Function value tolerance
      Scalar tol_x{ 1e-6 };     ///< Step size tolerance
      size_t patience{ 50 };    ///< Iterations without improvement before stop

      // Adaptive parameters
      bool   adaptive_poll{ true };     ///< Adaptive poll directions
      Scalar success_threshold{ 0.1 };  ///< Threshold for successful iteration
    };

    /**
     * @struct Result
     * @brief Minimization result
     */
    struct Result
    {
      Scalar final_f;       ///< Final function value
      Vector final_x;       ///< Final point
      size_t iterations;    ///< Number of iterations performed
      size_t f_eval_count;  ///< Number of function evaluations
      bool   converged;     ///< True if convergence reached
      Scalar mesh_size;     ///< Final mesh size
      string message;       ///< Termination message
    };

    /**
     * @struct IterationData
     * @brief Data for each iteration
     */
    struct IterationData
    {
      size_t iteration;
      Scalar f_value;
      Scalar mesh_size;
      bool   success;
      size_t eval_count;
    };

  private:
    Options m_opts;   ///< Configured options
    Vector  m_lower;  ///< Lower bounds
    Vector  m_upper;  ///< Upper bounds

    // Internal state
    Vector                     m_current_x;   ///< Current best point
    Scalar                     m_current_f;   ///< Current best function value
    Scalar                     m_mesh_size;   ///< Current mesh size
    size_t                     m_eval_count;  ///< Function evaluation counter
    std::vector<IterationData> m_history;     ///< Optimization history

    // Cache for visited points (to avoid re-evaluation)
    std::set<Vector, VectorComparator<Scalar>> m_visited_points;

    // Random number generator
    std::mt19937 m_rng{ std::random_device{}() };

  public:
    /**
     * @brief Constructor
     * @param o Options for the minimizer
     */
    MADS_minimizer( Options const & o = Options() ) : m_opts( o ) {}

    /**
     * @brief Set optimization bounds
     * @param lo Lower bounds
     * @param up Upper bounds
     */
    void
    set_bounds( Vector const & lo, Vector const & up )
    {
      assert( lo.size() == up.size() );
      m_lower = lo;
      m_upper = up;
    }

    void
    set_bounds( size_t n, Scalar const lower[], Scalar const upper[] )
    {
      m_lower.resize( n );
      m_upper.resize( n );
      std::copy_n( lower, n, m_lower.data() );
      std::copy_n( upper, n, m_upper.data() );
    }

    /**
     * @brief Project point onto bounds
     * @param x Point to project (modified in-place)
     */
    void
    project( Vector & x ) const
    {
      x = x.cwiseMax( m_lower ).cwiseMin( m_upper );
    }

  private:
    /**
     * @brief Generate poll directions using LTMADS method
     * @param directions Output vector of directions
     */
    void
    generate_ltmads_directions( std::vector<Vector> & directions )
    {
      size_t n = m_current_x.size();
      directions.clear();

      // Generate positive and negative coordinate directions
      for ( size_t i = 0; i < n; ++i )
      {
        for ( int sign = -1; sign <= 1; sign += 2 )
        {
          Vector dir = Vector::Zero( n );
          dir( i )   = static_cast<Scalar>( sign );
          directions.push_back( dir );
        }
      }

      // Add some random orthogonal directions for diversity
      std::normal_distribution<Scalar> dist( 0.0, 1.0 );

      for ( size_t i = 0; i < m_opts.poll_directions; ++i )
      {
        Vector dir( n );
        for ( size_t j = 0; j < n; ++j ) { dir( j ) = dist( m_rng ); }
        dir.normalize();
        directions.push_back( dir );

        // Also add negative direction
        directions.push_back( -dir );
      }
    }

    /**
     * @brief Generate orthogonal poll directions
     * @param directions Output vector of directions
     */
    void
    generate_orthomads_directions( std::vector<Vector> & directions ) const
    {
      size_t n = m_current_x.size();
      directions.clear();

      // Generate orthogonal basis directions
      Matrix                       Q = Matrix::Random( n, n );
      Eigen::HouseholderQR<Matrix> qr( Q );
      Q = qr.householderQ();

      // Use columns of Q as directions (both positive and negative)
      for ( size_t i = 0; i < n; ++i )
      {
        directions.push_back( Q.col( i ) );
        directions.push_back( -Q.col( i ) );
      }
    }

    /**
     * @brief Generate search points around current point
     * @param search_points Output vector of search points
     */
    void
    generate_search_points( std::vector<Vector> & search_points )
    {
      size_t n = m_current_x.size();
      search_points.clear();

      std::uniform_real_distribution<Scalar> dist( -1.0, 1.0 );

      // Generate random points in the mesh neighborhood
      for ( size_t i = 0; i < m_opts.search_points; ++i )
      {
        Vector point = m_current_x;
        for ( size_t j = 0; j < n; ++j )
        {
          Scalar step = m_mesh_size * dist( m_rng );
          point( j ) += step;
        }
        project( point );

        // Only add if not too close to current point
        if ( ( point - m_current_x ).norm() > m_mesh_size * 0.1 ) { search_points.push_back( point ); }
      }
    }

    /**
     * @brief Evaluate point with caching
     * @param point Point to evaluate
     * @param fun Objective function
     * @return Function value
     */
    Scalar
    evaluate_point( Vector const & point, Callback const & fun )
    {
      // Check if point was already evaluated
      auto it = m_visited_points.find( point );
      if ( it != m_visited_points.end() )
      {
        // In a real implementation, we would cache function values
        // For simplicity, we re-evaluate but still count the evaluation
        m_eval_count++;
        return fun( point );
      }

      m_eval_count++;
      Scalar value = fun( point );
      m_visited_points.insert( point );
      return value;
    }

    /**
     * @brief Perform search step
     * @param fun Objective function
     * @return True if improvement found
     */
    bool
    search_step( Callback const & fun )
    {
      if ( !m_opts.enable_search ) return false;

      std::vector<Vector> search_points;
      generate_search_points( search_points );

      Scalar best_f      = m_current_f;
      Vector best_x      = m_current_x;
      bool   improvement = false;

      for ( auto const & point : search_points )
      {
        if ( m_eval_count >= m_opts.max_evaluations ) break;

        Scalar f_val = evaluate_point( point, fun );

        if ( f_val < best_f - m_opts.success_threshold * abs( best_f ) )
        {
          best_f      = f_val;
          best_x      = point;
          improvement = true;
        }
      }

      if ( improvement )
      {
        m_current_x = best_x;
        m_current_f = best_f;
      }

      return improvement;
    }

    /**
     * @brief Perform poll step
     * @param fun Objective function
     * @return True if improvement found
     */
    bool
    poll_step( Callback const & fun )
    {
      std::vector<Vector> directions;

      // Generate poll directions based on selected method
      if ( m_opts.poll_method == "ortho" ) { generate_orthomads_directions( directions ); }
      else
      {
        generate_ltmads_directions( directions );
      }

      Scalar best_f      = m_current_f;
      Vector best_x      = m_current_x;
      bool   improvement = false;

      // Evaluate all poll points
      for ( auto const & dir : directions )
      {
        if ( m_eval_count >= m_opts.max_evaluations ) break;

        Vector poll_point = m_current_x + m_mesh_size * dir;
        project( poll_point );

        // Skip if too close to current point
        if ( ( poll_point - m_current_x ).norm() < m_mesh_size * 1e-6 ) { continue; }

        Scalar f_val = evaluate_point( poll_point, fun );

        if ( f_val < best_f - m_opts.success_threshold * abs( best_f ) )
        {
          best_f      = f_val;
          best_x      = poll_point;
          improvement = true;
          break;  // Stop at first improvement (opportunistic strategy)
        }
      }

      if ( improvement )
      {
        m_current_x = best_x;
        m_current_f = best_f;
      }

      return improvement;
    }

    /**
     * @brief Update mesh size based on iteration success
     * @param success Whether iteration was successful
     */
    void
    update_mesh_size( bool success )
    {
      if ( success )
      {
        // Successful iteration: expand mesh
        m_mesh_size = min( m_mesh_size * m_opts.mesh_expansion, m_opts.initial_mesh_size );
      }
      else
      {
        // Unsuccessful iteration: contract mesh
        m_mesh_size = max( m_mesh_size * m_opts.mesh_contraction, m_opts.min_mesh_size );
      }
    }

    /**
     * @brief Check convergence criteria
     * @param iteration Current iteration
     * @param no_improvement_count Count of iterations without improvement
     * @return Convergence result
     */
    Result
    check_convergence( size_t iteration, size_t no_improvement_count ) const
    {
      // Mesh size convergence
      if ( m_mesh_size <= m_opts.tol_mesh )
      {
        return { m_current_f, m_current_x, iteration, m_eval_count, true, m_mesh_size, "Mesh size below tolerance" };
      }

      // No improvement convergence
      if ( no_improvement_count >= m_opts.patience )
      {
        return {
          m_current_f, m_current_x, iteration, m_eval_count, true, m_mesh_size, "No improvement for patience iterations"
        };
      }

      // Maximum evaluations reached
      if ( m_eval_count >= m_opts.max_evaluations )
      {
        return {
          m_current_f, m_current_x, iteration, m_eval_count, true, m_mesh_size, "Maximum function evaluations reached"
        };
      }

      // Not converged
      return { m_current_f, m_current_x, iteration, m_eval_count, false, m_mesh_size, "" };
    }

  public:
    /**
     * @brief Perform minimization using MADS
     */
    Result
    minimize( Vector const & x0, Callback const & fun )
    {
      // Initialize
      m_current_x = x0;
      project( m_current_x );
      m_current_f  = fun( x0 );
      m_mesh_size  = m_opts.initial_mesh_size;
      m_eval_count = 1;
      m_history.clear();
      m_visited_points.clear();
      m_visited_points.insert( m_current_x );

      size_t n                    = x0.size();
      size_t no_improvement_count = 0;

      if ( m_opts.verbose )
      {
        fmt::print(
            "[MADS] Starting optimization, dimension: {}, initial f: {}, "
            "initial mesh: {}\n",
            n, m_current_f, m_mesh_size );
      }

      for ( size_t k{ 0 }; k < m_opts.max_iter; ++k )
      {
        bool improvement = false;

        // SEARCH STEP (optional)
        if ( m_opts.enable_search ) { improvement = search_step( fun ); }

        // POLL STEP (always performed)
        if ( !improvement ) { improvement = poll_step( fun ); }

        // Update mesh size based on success
        update_mesh_size( improvement );

        // Update improvement counter
        if ( improvement ) { no_improvement_count = 0; }
        else
        {
          no_improvement_count++;
        }

        // Store iteration data
        m_history.push_back( { k, m_current_f, m_mesh_size, improvement, m_eval_count } );

        // Verbose output
        if ( m_opts.verbose && ( k % m_opts.print_every ) == 0 )
        {
          fmt::print( "[MADS] iter={:<4} f={:<10.4} mesh={:<10.4} improv={} evals={}\n", k, m_current_f, m_mesh_size,
                      improvement, m_eval_count );
        }

        // Check convergence
        Result conv_result = check_convergence( k, no_improvement_count );
        if ( conv_result.converged )
        {
          if ( m_opts.verbose ) { fmt::print( "[MADS] Converged at iteration {}: {}\n", k, conv_result.message ); }
          return conv_result;
        }

        // Check maximum evaluations
        if ( m_eval_count >= m_opts.max_evaluations )
        {
          return {
            m_current_f, m_current_x, k, m_eval_count, true, m_mesh_size, "Maximum function evaluations reached"
          };
        }
      }

      return {
        m_current_f, m_current_x, m_opts.max_iter, m_eval_count, false, m_mesh_size, "Maximum iterations reached"
      };
    }

    /**
     * @brief Get optimization history
     * @return Vector of iteration data
     */
    std::vector<IterationData> const &
    get_history() const
    {
      return m_history;
    }

    /**
     * @brief Get current mesh size
     * @return Current mesh size
     */
    Scalar
    get_mesh_size() const
    {
      return m_mesh_size;
    }

    /**
     * @brief Get number of function evaluations
     * @return Evaluation count
     */
    size_t
    get_eval_count() const
    {
      return m_eval_count;
    }

    /**
     * @brief Get current options
     * @return Current options
     */
    Options const &
    get_options() const
    {
      return m_opts;
    }

    /**
     * @brief Set new options
     * @param opts New options
     */
    void
    set_options( Options const & opts )
    {
      m_opts = opts;
    }

    /**
     * @brief Get current state information
     * @return String with current state info
     */
    string
    get_state_info() const
    {
      return fmt::format( "MADS State: f={}, mesh_size={}, evaluations={}, visited_points={}", m_current_f, m_mesh_size,
                          m_eval_count, m_visited_points.size() );
    }

    /**
     * @brief Reset the optimizer state
     */
    void
    reset()
    {
      m_current_x  = Vector();
      m_current_f  = std::numeric_limits<Scalar>::max();
      m_mesh_size  = m_opts.initial_mesh_size;
      m_eval_count = 0;
      m_history.clear();
      m_visited_points.clear();
    }
  };

}  // namespace Utils

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif

#endif

//
// eof: Utils_MADS.hh
