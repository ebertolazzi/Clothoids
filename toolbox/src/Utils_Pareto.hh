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

/**
 * @file Utils_Pareto.hh
 * @brief Defines the ParetoFront class for maintaining the Pareto set in
 * multi-objective optimization.
 */

#pragma once

#ifndef UTILS_PARETO_dot_HH
#define UTILS_PARETO_dot_HH

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <variant>
#include <vector>

#include "Utils_eigen.hh"

namespace Utils
{

  using std::array;
  using std::size_t;

  /**
   * @brief A data structure for maintaining the Pareto Front of a set of
   * N-dimensional points.
   *
   * The **Pareto Front** (or Pareto Set) is the set of non-dominated solutions
   * in a multi-objective optimization problem. A solution (point) $A$ is
   * considered
   * **non-dominated** if no other solution $B$ in the set is superior (better
   * or equal) in all objectives and strictly superior in at least one
   * objective.
   *
   * This class supports both incremental updates (`insert`) and efficient batch
   * construction (`batch_build`). It uses a soft-delete mechanism (tombstones)
   * and a periodic rebuild strategy (`maybe_rebuild`) for consistent
   * performance during incremental erasures.
   *
   * @tparam T The numeric type for the point coordinates (e.g., int, float,
   * double).
   * @tparam N The dimensionality of the space (number of objectives). Must be >
   * 0.
   * @tparam Payload The type of auxiliary data stored with each point. Defaults
   * to void (std::monostate).
   * @tparam RebuildThreshold The number of deleted entries (`m_tombstones`)
   * that triggers a physical cleanup of the underlying vector.
   *
   * @section concept Dominance Concept
   * The core of the Pareto Front is the dominance relation. Given two points
   * $A$ and $B$:
   * - **A strictly dominates B** if $A$ is non-worse than $B$ in all objectives
   * and strictly better in at least one.
   * - **A weakly dominates B** if $A$ is non-worse than $B$ in all objectives.
   *
   * @section insertion Incremental Insertion (`insert`)
   * The `insert` method guarantees that the front remains correct after adding
   * a new point $P$:
   * 1. Check if $P$ is an exact duplicate (fail).
   * 2. Check if $P$ is weakly dominated by the existing front (fail, $P$ is
   * inferior).
   * 3. **Prune**: Identify and soft-delete all points $E$ in the current front
   * where $P$ strictly dominates $E$.
   * 4. Add $P$ as a new active entry.
   *
   * The use of soft-deletion (`Entry::alive = false`) maintains pointer
   * stability (if entries were accessed by index) but requires periodic cleanup
   * via `maybe_rebuild`.
   *
   * @section batch Batch Construction (`batch_build`)
   * The batch building method is typically more efficient for large, unsorted
   * datasets. It relies on sorting the input based on a scalar measure (L1
   * norm) to quickly filter out many dominated points.
   * 1. **Sort**: Input points are sorted by their L1 norm (sum of absolute
   * values) plus a lexicographical tie-breaker. This heuristic places
   * promising, low-norm (often near the origin) points earlier.
   * 2. **Iterate & Prune**: Iterate through the sorted points. For each
   * candidate $C$: a. Check if $C$ is weakly dominated by the current *partial*
   * front (expensive). b. If not dominated, remove any point in the partial
   * front that is strictly dominated by $C$. c. Add $C$ to the partial front.
   */
  template <typename T, size_t N, typename Payload = void, size_t RebuildThreshold = 128>
  class ParetoFront
  {
  public:
    /// @brief Type alias for an N-dimensional point. $\mathbf{p} \in
    /// \mathbb{R}^N$.
    using Point_type = std::array<T, N>;  // Eigen::Matrix<T,N,1>; // ;
    /// @brief Type alias for the payload, resolved to std::monostate if Payload
    /// is void.
    using Payload_type = std::conditional_t<std::is_same<Payload, void>::value, std::monostate, Payload>;

    /**
     * @brief Structure representing an entry (a point and its metadata) in the
     * front.
     */
    struct Entry
    {
      Point_type   p;              ///< The N-dimensional point coordinates.
      Payload_type payload;        ///< The auxiliary data associated with the point.
      bool         alive{ true };  ///< Status flag: true if active, false if soft-deleted (tombstone).
      size_t       id{ 0 };        ///< A unique identifier. Used for reliable
                                   ///< deletion/payload retrieval.
    };

    /// @brief Type alias for a vector of point-payload pairs, representing the
    /// current front.
    using Front = std::vector<std::pair<Point_type, Payload_type>>;

  private:
    /// @brief Internal storage vector containing all entries (active and
    /// tombstones).
    std::vector<Entry> m_entries;
    /// @brief Counter for soft-deleted entries. When this reaches
    /// `RebuildThreshold`, a physical clean-up occurs.
    size_t m_tombstones{ 0 };
    /// @brief Next unique ID to be assigned to a new entry.
    size_t m_next_id{ 1 };

    /**
     * @brief Checks for exact equality between two N-dimensional points.
     * @param a The first point.
     * @param b The second point.
     * @return true if $a_i = b_i$ for all $i \in [0, N-1]$.
     */
    static bool
    equals( Point_type const & a, Point_type const & b )
    {
      for ( size_t i{ 0 }; i < N; ++i )
        if ( a[i] != b[i] ) return false;
      return true;
    }

    /**
     * @brief Checks if point A strictly dominates point B.
     *
     * $A$ strictly dominates $B$ (i.e., $A \succ B$) if:
     * 1. $\forall i \in [0, N-1]$, $A$ is non-worse than $B$ in objective $i$.
     * 2. $\exists j \in [0, N-1]$, $A$ is strictly better than $B$ in objective
     * $j$.
     *
     * @param a The potential dominating point.
     * @param b The potentially dominated point.
     * @return true if $A \succ B$, false otherwise.
     */
    bool
    dominates( Point_type const & a, Point_type const & b ) const
    {
      bool any_strict{ false };
      for ( size_t i{ 0 }; i < N; ++i )
      {
        if ( a[i] > b[i] ) return false;       // A is worse than B
        if ( a[i] < b[i] ) any_strict = true;  // A is strictly better
      }
      return any_strict;
    }

    /**
     * @brief Checks if point A weakly dominates point B.
     *
     * $A$ weakly dominates $B$ (i.e., $A \succeq B$) if:
     * 1. $\forall i \in [0, N-1]$, $A$ is non-worse than $B$ in objective $i$.
     *
     * @param a The potential weakly dominating point.
     * @param b The potentially weakly dominated point.
     * @return true if $A \succeq B$, false otherwise.
     */
    bool
    dominates_weakly( Point_type const & a, Point_type const & b ) const
    {
      for ( size_t i{ 0 }; i < N; ++i )
        if ( a[i] > b[i] ) return false;
      return true;
    }

    /**
     * @brief Performs a physical clean-up of the internal storage.
     *
     * This method is called by `insert` and `erase_by_id` when the number
     * of soft-deleted entries (`m_tombstones`) reaches `RebuildThreshold`.
     * It copies only the `alive` entries into a new vector, effectively
     * removing the tombstones to maintain spatial locality and iteration
     * performance.
     *
     * \dot
     * digraph maybe_rebuild {
     * node [shape=box];
     * Start [label="Start maybe_rebuild()"];
     * CheckThreshold [label="m_tombstones >= RebuildThreshold?"];
     * Rebuild [label="Create new_entries vector"];
     * CopyAlive [label="Iterate m_entries, copy if e.alive"];
     * Swap [label="m_entries = std::move(new_entries)"];
     * Reset [label="m_tombstones = 0"];
     * End [label="End"];
     *
     * Start -> CheckThreshold;
     * CheckThreshold -> End [label="No"];
     * CheckThreshold -> Rebuild [label="Yes"];
     * Rebuild -> CopyAlive;
     * CopyAlive -> Swap;
     * Swap -> Reset;
     * Reset -> End;
     * }
     * \enddot
     */
    void
    maybe_rebuild()
    {
      if ( m_tombstones >= RebuildThreshold )
      {
        std::vector<Entry> new_entries;
        new_entries.reserve( m_entries.size() - m_tombstones );
        for ( auto const & e : m_entries )
        {
          if ( e.alive ) new_entries.emplace_back( e );
        }
        m_entries    = std::move( new_entries );
        m_tombstones = 0;
      }
    }

    /**
     * @brief Comparison function for strict-weak ordering based on L1 norm.
     *
     * This ordering is crucial for the efficiency of `batch_build`.
     * 1. **Primary Key**: L1 Norm $\sum_{i=0}^{N-1} |p_i|$. Sorting by
     * ascending norm ensures that points closer to the origin (often better
     * candidates in normalized spaces) are processed first.
     * 2. **Tie-Breaker**: Lexicographical ordering.
     *
     * @param a The first point.
     * @param b The second point.
     * @return true if a comes before b in the sorted order.
     */
    static bool
    compare_points( Point_type const & a, Point_type const & b )
    {
      // Calculate L1 norm
      T norm_a = T( 0 );
      T norm_b = T( 0 );
      for ( size_t i{ 0 }; i < N; ++i )
      {
        norm_a += std::abs( a[i] );
        norm_b += std::abs( b[i] );
      }

      // 1. Sort by L1 norm ascending
      if ( norm_a != norm_b ) return norm_a < norm_b;

      // 2. Lexicographical ordering as a tie-breaker
      for ( size_t i{ 0 }; i < N; ++i )
        if ( a[i] != b[i] ) return a[i] < b[i];

      // Points are identical
      return false;
    }

  public:
    static_assert( N > 0, "N must be > 0" );

    /**
     * @brief Constructs a ParetoFront object.
     * @param minimize An array of booleans defining the optimization direction
     * for each dimension.
     */
    explicit ParetoFront() {}

    /// @brief Returns the number of **active** entries in the front.
    size_t
    size() const
    {
      return m_entries.size() - m_tombstones;
    }

    /// @brief Returns the total number of entries (active + tombstones).
    size_t
    raw_size() const
    {
      return m_entries.size();
    }

    /// @brief Checks if the front is empty (size() == 0).
    bool
    empty() const
    {
      return size() == 0;
    }

    /**
     * @brief Iterates over all active entries and applies a function object.
     * @tparam F Type of the function object (e.g., a lambda).
     * @param f The function to apply to each alive Entry. Signature must accept
     * a const Entry&.
     */
    template <typename F>
    void
    for_each_alive( F && f ) const
    {
      for ( auto const & e : m_entries )
        if ( e.alive ) f( e );
    }

    /**
     * @brief Returns a vector of (point, payload) pairs representing the
     * current Pareto Front.
     * @return A vector of Front type.
     */
    Front
    front() const
    {
      Front out;
      out.reserve( size() );
      for ( auto const & e : m_entries )
        if ( e.alive ) out.emplace_back( e.p, e.payload );
      return out;
    }

    /**
     * @brief Checks if the given point $p$ is weakly dominated by any active
     * point in the front.
     *
     * This acts as a quick filter for new points. If $p$ is dominated by the
     * existing front, it cannot be part of the new front, and the expensive
     * dominance checks against $p$ are avoided.
     *
     * @param p The point to check.
     * @return true if $p$ is weakly dominated by an existing point, false
     * otherwise.
     */
    bool
    is_dominated_by_front( Point_type const & p ) const
    {
      for ( auto const & e : m_entries )
        if ( e.alive && dominates_weakly( e.p, p ) ) return true;
      return false;
    }

    /**
     * @brief Attempts to insert a new point $P$ into the Pareto Front.
     *
     * @param p The point to insert.
     * @param payload The auxiliary data for the point.
     * @return A pair: {success_of_insertion, unique_id_of_new_entry}. unique_id
     * is 0 if insertion failed.
     *
     * \dot
     * digraph insert {
     * node [shape=box];
     * Start [label="Start insert(P)"];
     * CheckDuplicate [label="P is exact duplicate?"];
     * CheckDominated [label="P dominated by current front?"];
     * PruneFront [label="Soft-delete all E where P strictly dominates E (P
     * \succ E)"]; InsertP [label="Add P to m_entries, assign ID"]; Rebuild
     * [label="maybe_rebuild()"]; Success [label="Return {true, ID}"]; Fail
     * [label="Return {false, 0}"];
     *
     * Start -> CheckDuplicate;
     * CheckDuplicate -> Fail [label="Yes"];
     * CheckDuplicate -> CheckDominated [label="No"];
     * CheckDominated -> Fail [label="Yes"];
     * CheckDominated -> PruneFront [label="No"];
     * PruneFront -> InsertP;
     * InsertP -> Rebuild;
     * Rebuild -> Success;
     * }
     * \enddot
     */
    std::pair<bool, size_t>
    insert( Point_type const & p, Payload_type const & payload = Payload_type{} )
    {
      // 1. Controlla duplicati esatti
      for ( auto const & e : m_entries )
        if ( e.alive && equals( e.p, p ) ) return { false, 0 };

      // 2. Controlla se il nuovo punto è dominato dal fronte
      if ( is_dominated_by_front( p ) ) return { false, 0 };

      // 3. Elimina tutti i punti strettamente dominati dal nuovo punto
      size_t tombstones_before = m_tombstones;
      for ( auto & e : m_entries )
      {
        if ( e.alive && dominates( p, e.p ) )
        {
          e.alive = false;
          ++m_tombstones;
        }
      }

      // 4. Inserisci nuovo punto
      size_t e_id{ m_next_id++ };
      m_entries.emplace_back( Entry{ p, payload, true, e_id } );

      // 5. Ricostruisci se necessario (solo se abbiamo effettivamente rimosso
      // punti)
      if ( m_tombstones > tombstones_before ) maybe_rebuild();

      return { true, e_id };
    }

    /**
     * @brief Erases a point from the front using its unique ID.
     *
     * Performs a soft-delete and then immediately checks for rebuild to keep
     * the performance profile consistent after removal.
     *
     * @param id The unique identifier of the entry to erase. Must be > 0.
     * @return true if the entry was found and erased, false otherwise.
     */
    bool
    erase_by_id( size_t id )
    {
      if ( id == 0 ) return false;
      for ( auto & e : m_entries )
      {
        if ( e.alive && e.id == id )
        {
          e.alive = false;
          ++m_tombstones;
          maybe_rebuild();
          return true;
        }
      }
      return false;
    }

    /**
     * @brief Finds the ID of the alive point in the front that is closest to
     * the given point $p$.
     *
     * The distance metric used is the **squared Euclidean distance** ($L_2^2$
     * norm):
     * $$d(e, p)^2 = \sum_{i=0}^{N-1} (e.p[i] - p[i])^2$$
     *
     * @param p The reference point.
     * @return A pair: {found\_status, unique\_id\_of\_nearest\_point}.
     * unique\_id is 0 if no alive entry is found.
     */
    std::pair<bool, size_t>
    find_nearest( Point_type const & p ) const
    {
      double bestd  = std::numeric_limits<double>::infinity();
      size_t bestid = 0;
      bool   found  = false;
      for ( auto const & e : m_entries )
      {
        if ( !e.alive ) continue;
        double d = 0.0;
        for ( size_t i{ 0 }; i < N; ++i )
        {
          // Calculates squared Euclidean distance
          double diff = double( e.p[i] ) - double( p[i] );
          d += diff * diff;
        }
        if ( !found || d < bestd )
        {
          bestd  = d;
          bestid = e.id;
          found  = true;
        }
      }
      return { found, bestid };
    }

    /**
     * @brief Builds the Pareto Front from a batch of points, clearing any
     * existing state.
     *
     * This is the preferred method for initial construction or when the entire
     * set of points is available. The algorithm complexity is dominated by the
     * initial sort $O(M \log M)$ where $M$ is the number of points, followed by
     * a sweep $O(M^2)$ in the worst-case, but generally much faster in
     * practice, especially for low-dimensional problems.
     *
     * @param pts A vector of (point, payload) pairs to build the front from.
     *
     * \dot
     * digraph batch_build {
     * node [shape=box];
     * Start [label="Start batch_build(PTS)"];
     * Clear [label="Clear state"];
     * Sort [label="Sort PTS by L1-Norm (compare_points)"];
     * Loop [label="For each Candidate C in sorted PTS"];
     * CheckDom [label="C dominated by current Front F?"];
     * PruneF [label="Remove E \in F where C \succ E"];
     * AddC [label="Add C to Front F"];
     * End [label="m_entries = F; End"];
     *
     * Start -> Clear;
     * Clear -> Sort;
     * Sort -> Loop;
     * Loop -> CheckDom;
     * CheckDom -> Loop [label="Yes"];
     * CheckDom -> PruneF [label="No"];
     * PruneF -> AddC;
     * AddC -> Loop;
     * Loop -> End [label="Done"];
     * }
     * \enddot
     */
    void
    batch_build( Front pts )
    {
      m_entries.clear();
      m_tombstones = 0;
      m_next_id    = 1;

      // 1. Sort points using L1-norm based strict-weak ordering
      std::sort( pts.begin(), pts.end(),
                 []( auto const & a, auto const & b ) { return compare_points( a.first, b.first ); } );

      std::vector<Entry> front;
      for ( auto const & pr : pts )
      {
        Point_type const &   candidate_point   = pr.first;
        Payload_type const & candidate_payload = pr.second;

        bool dominated{ false };

        // 2a. Check if the candidate is dominated by any point in the current
        // partial front (F)
        for ( size_t j{ 0 }; j < front.size(); ++j )
        {
          if ( dominates_weakly( front[j].p, candidate_point ) )
          {
            dominated = true;
            break;
          }
        }
        if ( dominated ) continue;

        // 2b. Remove points from F that are strictly dominated by the candidate
        // (C)
        for ( size_t j{ 0 }; j < front.size(); )
        {
          if ( dominates( candidate_point, front[j].p ) )
          {
            // Note: In batch_build, we use physical erase since the vector
            // 'front' is temporary.
            front.erase( front.begin() + j );
          }
          else
          {
            ++j;
          }
        }

        // 2c. Add the candidate to F
        front.emplace_back( Entry{ candidate_point, candidate_payload, true, m_next_id++ } );
      }
      m_entries = std::move( front );
    }

    /// @brief Clears the entire Pareto Front, resetting internal state,
    /// entries, and counters.
    void
    clear()
    {
      m_entries.clear();
      m_tombstones = 0;
      m_next_id    = 1;
    }

    /**
     * @brief Retrieves the payload associated with a given unique ID.
     * @param id The unique identifier. Must be > 0.
     * @return A pair: {found\_status, payload}. Payload is default-constructed
     * if not found.
     */
    std::pair<bool, Payload_type>
    get_payload_by_id( size_t id ) const
    {
      for ( auto const & e : m_entries )
        if ( e.alive && e.id == id ) return { true, e.payload };
      return { false, Payload_type{} };
    }
  };

}  // namespace Utils

#endif  // !UTILS_PARETO_dot_HH
