/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
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
 |      Enhanced and improved version                                       |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// file: Utils_GG2D.hh
//
#pragma once

#ifndef UTILS_GG2D_dot_HH
#define UTILS_GG2D_dot_HH

#include <algorithm>
#include <cmath>
#include <vector>

#include "Utils_eigen.hh"

namespace Utils
{

  /*!
   * \addtogroup TwoD
   * @{
   */

  /*\
   |   ____       _       _   ____  ____
   |  |  _ \ ___ (_)_ __ | |_|___ \|  _ \
   |  | |_) / _ \| | '_ \| __| __) | | | |
   |  |  __/ (_) | | | | | |_ / __/| |_| |
   |  |_|   \___/|_|_| |_|\__|_____|____/
  \*/

  //!
  //! \brief A class representing a 2D point in space.
  //!
  //! The `Point2D` class extends the `Eigen::Matrix` class to provide a
  //! convenient representation of a point in a two-dimensional space. It
  //! encapsulates operations related to 2D points, allowing easy access to the
  //! x and y coordinates. This class leverages the Eigen library for efficient
  //! matrix and vector computations.
  //!
  //! \tparam Real The data type used for the point coordinates, typically a
  //! floating-point type (e.g., float, double).
  //!
  template <typename Real>
  class Point2D : public Eigen::Matrix<Real, 2, 1>
  {
    using P2D = Eigen::Matrix<Real, 2, 1>;

  public:
    using typename P2D::Matrix;  // Inherit Eigen constructors

    //! Default constructor
    Point2D() = default;

    //! Constructor from coordinates
    Point2D( Real x, Real y )
    {
      this->coeffRef( 0 ) = x;
      this->coeffRef( 1 ) = y;
    }

    //! Constructor from Eigen vector
    explicit Point2D( P2D const & v ) : P2D( v ) {}

    //! Get x-coordinate
    Real
    x() const noexcept
    {
      return this->coeff( 0 );
    }

    //! Get y-coordinate
    Real
    y() const noexcept
    {
      return this->coeff( 1 );
    }

    //! Set x-coordinate
    void
    x( Real val ) noexcept
    {
      this->coeffRef( 0 ) = val;
    }

    //! Set y-coordinate
    void
    y( Real val ) noexcept
    {
      this->coeffRef( 1 ) = val;
    }

    //! Convert to const Eigen matrix
    P2D const &
    to_eigen() const noexcept
    {
      return *static_cast<P2D const *>( this );
    }

    //! Convert to Eigen matrix
    P2D &
    to_eigen() noexcept
    {
      return *static_cast<P2D *>( this );
    }

    //! Distance to another point
    Real
    distance( Point2D<Real> const & other ) const noexcept
    {
      return ( *this - other ).norm();
    }

    //! Squared distance to another point (faster, avoids sqrt)
    Real
    distance_squared( Point2D<Real> const & other ) const noexcept
    {
      return ( *this - other ).squaredNorm();
    }

    //! Dot product with another point
    Real
    dot( Point2D<Real> const & other ) const noexcept
    {
      return this->x() * other.x() + this->y() * other.y();
    }

    //! Cross product with another point (returns scalar z-component)
    Real
    cross( Point2D<Real> const & other ) const noexcept
    {
      return this->x() * other.y() - this->y() * other.x();
    }

    //! Normalize the point (treat as vector)
    Point2D<Real>
    normalized() const
    {
      Real len = this->norm();
      if ( len < machine_eps<Real>() ) { return Point2D<Real>( 0, 0 ); }
      return Point2D<Real>( this->x() / len, this->y() / len );
    }

    //! Rotate point by angle (radians) around origin
    Point2D<Real>
    rotate( Real angle ) const noexcept
    {
      Real c = std::cos( angle );
      Real s = std::sin( angle );
      return Point2D<Real>( c * this->x() - s * this->y(), s * this->x() + c * this->y() );
    }

    //! Rotate point by angle around a center point
    Point2D<Real>
    rotate( Real angle, Point2D<Real> const & center ) const noexcept
    {
      Point2D<Real> translated = *this - center;
      return translated.rotate( angle ) + center;
    }

    //! Scale point relative to origin
    Point2D<Real>
    scale( Real sx, Real sy ) const noexcept
    {
      return Point2D<Real>( this->x() * sx, this->y() * sy );
    }

    //! Scale point relative to a center point
    Point2D<Real>
    scale( Real sx, Real sy, Point2D<Real> const & center ) const noexcept
    {
      Point2D<Real> translated = *this - center;
      return Point2D<Real>( translated.x() * sx, translated.y() * sy ) + center;
    }

    //! Check if point is approximately equal to another point
    bool
    isApprox( Point2D<Real> const & other, Real prec = machine_eps<Real>() ) const noexcept
    {
      return std::abs( this->x() - other.x() ) <= prec && std::abs( this->y() - other.y() ) <= prec;
    }

    //! Stream output operator
    friend std::ostream &
    operator<<( std::ostream & os, Point2D<Real> const & p )
    {
      os << "(" << p.x() << ", " << p.y() << ")";
      return os;
    }
  };

  /*\
   |   ____            ____  ____
   |  | __ )  _____  _|___ \|  _ \
   |  |  _ \ / _ \ \/ / __) | | | |
   |  | |_) | (_) >  < / __/| |_| |
   |  |____/ \___/_/\_\_____|____/
  \*/

  //!
  //! \brief A class representing a 2D axis-aligned bounding box.
  //!
  //! The `Box2D` class defines a 2D axis-aligned bounding box using two points,
  //! \f$ P_{\text{min}} \f$ and \f$ P_{\text{max}} \f$, which represent the
  //! minimum and maximum corners of the box respectively.
  //!
  //! \tparam Real The data type used for the coordinates of the box's corners.
  //!
  template <typename Real>
  class Box2D
  {
    Point2D<Real> m_Pmin;  //!< Minimum corner (bottom-left).
    Point2D<Real> m_Pmax;  //!< Maximum corner (top-right).

  public:
    //! Default constructor (creates empty box at origin)
    Box2D() : m_Pmin( 0, 0 ), m_Pmax( 0, 0 ) {}

    //! Constructor from min and max points
    Box2D( Point2D<Real> const & pmin, Point2D<Real> const & pmax ) noexcept : m_Pmin( pmin ), m_Pmax( pmax ) {}

    //! Constructor from coordinates
    Box2D( Real xmin, Real ymin, Real xmax, Real ymax ) noexcept : m_Pmin( xmin, ymin ), m_Pmax( xmax, ymax ) {}

    //! Constructor from a single point
    explicit Box2D( Point2D<Real> const & p ) noexcept : m_Pmin( p ), m_Pmax( p ) {}

    //! Get minimum corner
    Point2D<Real> const &
    Pmin() const noexcept
    {
      return m_Pmin;
    }

    //! Get maximum corner
    Point2D<Real> const &
    Pmax() const noexcept
    {
      return m_Pmax;
    }

    //! Set minimum corner
    void
    Pmin( Point2D<Real> const & p ) noexcept
    {
      m_Pmin = p;
    }

    //! Set maximum corner
    void
    Pmax( Point2D<Real> const & p ) noexcept
    {
      m_Pmax = p;
    }

    //! Set both corners
    void
    set( Point2D<Real> const & pmin, Point2D<Real> const & pmax ) noexcept
    {
      m_Pmin = pmin;
      m_Pmax = pmax;
    }

    //! Set from coordinates
    void
    set( Real xmin, Real ymin, Real xmax, Real ymax ) noexcept
    {
      m_Pmin.x( xmin );
      m_Pmin.y( ymin );
      m_Pmax.x( xmax );
      m_Pmax.y( ymax );
    }

    //! Get width of box
    Real
    width() const noexcept
    {
      return m_Pmax.x() - m_Pmin.x();
    }

    //! Get height of box
    Real
    height() const noexcept
    {
      return m_Pmax.y() - m_Pmin.y();
    }

    //! Get area of box
    Real
    area() const noexcept
    {
      return width() * height();
    }

    //! Get perimeter of box
    Real
    perimeter() const noexcept
    {
      return 2 * ( width() + height() );
    }

    //! Get center of box
    Point2D<Real>
    center() const noexcept
    {
      return Point2D<Real>( ( m_Pmin.x() + m_Pmax.x() ) * Real( 0.5 ), ( m_Pmin.y() + m_Pmax.y() ) * Real( 0.5 ) );
    }

    //! Get size as a point (width, height)
    Point2D<Real>
    size() const noexcept
    {
      return Point2D<Real>( width(), height() );
    }

    //! Check if box is valid (max >= min)
    bool
    is_valid() const noexcept
    {
      return m_Pmax.x() >= m_Pmin.x() && m_Pmax.y() >= m_Pmin.y();
    }

    //! Check if box is empty (zero area)
    bool
    is_empty() const noexcept
    {
      return width() <= machine_eps<Real>() || height() <= machine_eps<Real>();
    }

    //! Check if point is inside box (inclusive boundaries)
    bool
    contains( Point2D<Real> const & P ) const noexcept
    {
      return P.x() >= m_Pmin.x() && P.x() <= m_Pmax.x() && P.y() >= m_Pmin.y() && P.y() <= m_Pmax.y();
    }

    //! Check if point is strictly inside box (exclusive boundaries)
    bool
    contains_strict( Point2D<Real> const & P ) const noexcept
    {
      return P.x() > m_Pmin.x() && P.x() < m_Pmax.x() && P.y() > m_Pmin.y() && P.y() < m_Pmax.y();
    }

    //! Check if another box is fully inside this box
    bool
    contains( Box2D<Real> const & other ) const noexcept
    {
      return contains( other.m_Pmin ) && contains( other.m_Pmax );
    }

    //! Check if boxes intersect
    bool
    intersects( Box2D<Real> const & other ) const noexcept
    {
      return !( m_Pmax.x() < other.m_Pmin.x() || m_Pmin.x() > other.m_Pmax.x() || m_Pmax.y() < other.m_Pmin.y() ||
                m_Pmin.y() > other.m_Pmax.y() );
    }

    //! Expand box to include point
    void
    expand( Point2D<Real> const & P ) noexcept
    {
      m_Pmin.x( std::min( m_Pmin.x(), P.x() ) );
      m_Pmin.y( std::min( m_Pmin.y(), P.y() ) );
      m_Pmax.x( std::max( m_Pmax.x(), P.x() ) );
      m_Pmax.y( std::max( m_Pmax.y(), P.y() ) );
    }

    //! Expand box by a margin in all directions
    void
    expand( Real margin ) noexcept
    {
      m_Pmin.x( m_Pmin.x() - margin );
      m_Pmin.y( m_Pmin.y() - margin );
      m_Pmax.x( m_Pmax.x() + margin );
      m_Pmax.y( m_Pmax.y() + margin );
    }

    //! Expand box by different margins in x and y
    void
    expand( Real margin_x, Real margin_y ) noexcept
    {
      m_Pmin.x( m_Pmin.x() - margin_x );
      m_Pmin.y( m_Pmin.y() - margin_y );
      m_Pmax.x( m_Pmax.x() + margin_x );
      m_Pmax.y( m_Pmax.y() + margin_y );
    }

    //! Expand box to include another box
    void
    expand( Box2D<Real> const & other ) noexcept
    {
      expand( other.m_Pmin );
      expand( other.m_Pmax );
    }

    //! Compute intersection of two boxes
    Box2D<Real>
    intersection( Box2D<Real> const & other ) const
    {
      Box2D<Real> result;
      result.m_Pmin.x( std::max( m_Pmin.x(), other.m_Pmin.x() ) );
      result.m_Pmin.y( std::max( m_Pmin.y(), other.m_Pmin.y() ) );
      result.m_Pmax.x( std::min( m_Pmax.x(), other.m_Pmax.x() ) );
      result.m_Pmax.y( std::min( m_Pmax.y(), other.m_Pmax.y() ) );

      // Ensure valid box
      if ( result.m_Pmin.x() > result.m_Pmax.x() || result.m_Pmin.y() > result.m_Pmax.y() )
      {
        // Return empty box
        result.m_Pmin = result.m_Pmax = Point2D<Real>( 0, 0 );
      }
      return result;
    }

    //! Scale box by factor around its center
    void
    scale( Real factor ) noexcept
    {
      if ( factor <= 0 ) return;

      Point2D<Real> c           = center();
      Real          half_width  = width() * factor * Real( 0.5 );
      Real          half_height = height() * factor * Real( 0.5 );

      m_Pmin.x( c.x() - half_width );
      m_Pmin.y( c.y() - half_height );
      m_Pmax.x( c.x() + half_width );
      m_Pmax.y( c.y() + half_height );
    }

    //! Scale box by different factors in x and y around its center
    void
    scale( Real factor_x, Real factor_y ) noexcept
    {
      if ( factor_x <= 0 || factor_y <= 0 ) return;

      Point2D<Real> c           = center();
      Real          half_width  = width() * factor_x * Real( 0.5 );
      Real          half_height = height() * factor_y * Real( 0.5 );

      m_Pmin.x( c.x() - half_width );
      m_Pmin.y( c.y() - half_height );
      m_Pmax.x( c.x() + half_width );
      m_Pmax.y( c.y() + half_height );
    }

    //! Translate box by a vector
    void
    translate( Point2D<Real> const & delta ) noexcept
    {
      m_Pmin = m_Pmin + delta;
      m_Pmax = m_Pmax + delta;
    }

    //! Translate box by coordinates
    void
    translate( Real dx, Real dy ) noexcept
    {
      m_Pmin.x( m_Pmin.x() + dx );
      m_Pmin.y( m_Pmin.y() + dy );
      m_Pmax.x( m_Pmax.x() + dx );
      m_Pmax.y( m_Pmax.y() + dy );
    }

    //! Get corner points of the box
    std::array<Point2D<Real>, 4>
    corners() const noexcept
    {
      return { { m_Pmin, Point2D<Real>( m_Pmax.x(), m_Pmin.y() ), m_Pmax, Point2D<Real>( m_Pmin.x(), m_Pmax.y() ) } };
    }

    //! Check if box is approximately equal to another box
    bool
    isApprox( Box2D<Real> const & other, Real prec = machine_eps<Real>() ) const noexcept
    {
      return m_Pmin.isApprox( other.m_Pmin, prec ) && m_Pmax.isApprox( other.m_Pmax, prec );
    }

    //! Stream output operator
    friend std::ostream &
    operator<<( std::ostream & os, Box2D<Real> const & box )
    {
      os << "Box[min=" << box.Pmin() << ", max=" << box.Pmax() << ", size=" << box.size() << "]";
      return os;
    }
  };

  /*\
   |   ____                                  _   ____  ____
   |  / ___|  ___  __ _ _ __ ___   ___ _ __ | |_|___ \|  _ \
   |  \___ \ / _ \/ _` | '_ ` _ \ / _ \ '_ \| __| __) | | | |
   |   ___) |  __/ (_| | | | | | |  __/ | | | |_ / __/| |_| |
   |  |____/ \___|\__, |_| |_| |_|\___|_| |_|\__|_____|____/
   |              |___/
  \*/

  //!
  //! \brief A class representing a 2D line segment defined by two endpoints.
  //!
  //! The `Segment2D` class encapsulates a line segment in a two-dimensional
  //! space, defined by two points, \f$ P_a \f$ and \f$ P_b \f$. It provides
  //! methods for setting the endpoints, projecting points onto the segment,
  //! evaluating points on the segment, and checking for intersections with
  //! other segments.
  //!
  //! \tparam Real The data type used for the coordinates of the segment's
  //! endpoints, typically a floating-point type (e.g., float, double).
  //!
  template <typename Real>
  class Segment2D
  {
    Point2D<Real> m_Pa;  //!< First endpoint of the segment.
    Point2D<Real> m_Pb;  //!< Second endpoint of the segment.

    static constexpr Real DEFAULT_TOLERANCE = Real( 1e-10 );

    //! Helper: compute orthogonal normal vector (rotated 90° CCW)
    Point2D<Real>
    compute_normal() const noexcept
    {
      Point2D<Real> BA = m_Pb - m_Pa;
      Point2D<Real> N;
      N.x( -BA.y() );
      N.y( BA.x() );
      return N;
    }

  public:
    //! Default constructor
    Segment2D() = default;

    //! Copy constructor
    Segment2D( Segment2D<Real> const & S ) noexcept : m_Pa( S.m_Pa ), m_Pb( S.m_Pb ) {}

    //! Constructor from two points
    Segment2D( Point2D<Real> const & A, Point2D<Real> const & B ) noexcept : m_Pa( A ), m_Pb( B ) {}

    //! Constructor from coordinates
    Segment2D( Real ax, Real ay, Real bx, Real by ) noexcept : m_Pa( ax, ay ), m_Pb( bx, by ) {}

    //! Move constructor
    Segment2D( Segment2D<Real> && S ) noexcept : m_Pa( std::move( S.m_Pa ) ), m_Pb( std::move( S.m_Pb ) ) {}

    //! Copy assignment
    Segment2D &
    operator=( Segment2D<Real> const & S ) noexcept
    {
      if ( this != &S )
      {
        m_Pa = S.m_Pa;
        m_Pb = S.m_Pb;
      }
      return *this;
    }

    //! Move assignment
    Segment2D &
    operator=( Segment2D<Real> && S ) noexcept
    {
      if ( this != &S )
      {
        m_Pa = std::move( S.m_Pa );
        m_Pb = std::move( S.m_Pb );
      }
      return *this;
    }

    //! Setup endpoints using Point2D objects
    void
    setup( Point2D<Real> const & A, Point2D<Real> const & B ) noexcept
    {
      m_Pa = A;
      m_Pb = B;
    }

    //! Setup endpoints using arrays
    void
    setup( Real const A[], Real const B[] ) noexcept
    {
      m_Pa.x( A[0] );
      m_Pa.y( A[1] );
      m_Pb.x( B[0] );
      m_Pb.y( B[1] );
    }

    //! Setup endpoints using coordinates
    void
    setup( Real ax, Real ay, Real bx, Real by ) noexcept
    {
      m_Pa.x( ax );
      m_Pa.y( ay );
      m_Pb.x( bx );
      m_Pb.y( by );
    }

    //! Check if segment is valid (not degenerate)
    bool
    is_valid( Real tol = DEFAULT_TOLERANCE ) const noexcept
    {
      return m_Pa.distance_squared( m_Pb ) > tol * tol;
    }

    //! Get segment length
    Real
    length() const noexcept
    {
      return m_Pa.distance( m_Pb );
    }

    //! Get squared segment length (faster)
    Real
    length_squared() const noexcept
    {
      return m_Pa.distance_squared( m_Pb );
    }

    //! Get direction vector (not normalized)
    Point2D<Real>
    direction() const noexcept
    {
      return m_Pb - m_Pa;
    }

    //! Get normalized direction vector
    Point2D<Real>
    direction_normalized() const
    {
      Point2D<Real> dir = direction();
      Real          len = dir.norm();
      if ( len < machine_eps<Real>() ) return Point2D<Real>( 0, 0 );
      return dir / len;
    }

    //! Get normal vector (rotated 90° CCW, not normalized)
    Point2D<Real>
    normal() const noexcept
    {
      return compute_normal();
    }

    //! Get normalized normal vector
    Point2D<Real>
    normal_normalized() const
    {
      Point2D<Real> n   = compute_normal();
      Real          len = n.norm();
      if ( len < machine_eps<Real>() ) { return Point2D<Real>( 0, 0 ); }
      return n / len;
    }

    //!
    //! \brief Projects a point onto the segment.
    //!
    //! \param P The point to be projected onto the segment.
    //! \param s The parameter along the segment where the projection occurs
    //! (output).
    //!
    //! \return `true` if the projection lies within [0,1] (on the segment),
    //! `false` otherwise.
    //!
    bool
    projection( Point2D<Real> const & P, Real & s ) const noexcept
    {
      Point2D<Real> const BA        = m_Pb - m_Pa;
      Real const          length_sq = BA.squaredNorm();

      // Handle degenerate segment
      if ( length_sq < machine_eps<Real>() )
      {
        s = 0;
        return false;
      }

      Real const numerator = ( P - m_Pa ).dot( BA );
      s                    = numerator / length_sq;
      return s >= 0 && s <= 1;
    }

    //!
    //! \brief Projects a point onto the segment and returns the projected
    //! point.
    //!
    //! \param P The point to be projected onto the segment.
    //! \param s The parameter along the segment where the projection occurs
    //! (output).
    //! \param t The signed distance from the segment (output).
    //!
    //! \return The projected point on the segment (clamped to [0,1]).
    //!
    Point2D<Real>
    projection( Point2D<Real> const & P, Real & s, Real & t ) const
    {
      Point2D<Real> const BA        = m_Pb - m_Pa;
      Real const          length_sq = BA.squaredNorm();

      // Handle degenerate segment
      if ( length_sq < machine_eps<Real>() )
      {
        s = 0;
        t = P.distance( m_Pa );
        return m_Pa;
      }

      Real const numerator = ( P - m_Pa ).dot( BA );
      s                    = numerator / length_sq;

      // Clamp s to [0, 1]
      if ( s < 0 )
        s = 0;
      else if ( s > 1 )
        s = 1;

      Point2D<Real> const PP = m_Pa + s * BA;

      // Compute signed distance using normal
      Point2D<Real> N        = compute_normal();
      Real const    norm_len = N.norm();
      if ( norm_len > machine_eps<Real>() ) { t = ( P - PP ).dot( N ) / norm_len; }
      else
      {
        t = 0;
      }

      return PP;
    }

    //!
    //! \brief Evaluates the segment at a given parameter.
    //!
    //! \param s The parameter to evaluate (typically in [0,1]).
    //!
    //! \return The point on the segment at parameter s.
    //!
    Point2D<Real>
    eval( Real s ) const noexcept
    {
      return m_Pa + s * ( m_Pb - m_Pa );
    }

    //!
    //! \brief Evaluates a point offset from the segment.
    //!
    //! \param s The parameter along the segment.
    //! \param t The signed distance perpendicular to the segment.
    //!
    //! \return The point at position s along the segment, offset by distance t.
    //!
    Point2D<Real>
    eval( Real s, Real t ) const
    {
      Point2D<Real> N   = compute_normal();
      Real          len = N.norm();
      if ( len > machine_eps<Real>() ) { N = N / len; }
      return m_Pa + s * ( m_Pb - m_Pa ) + t * N;
    }

    //! Get first endpoint
    Point2D<Real> const &
    Pa() const noexcept
    {
      return m_Pa;
    }

    //! Get second endpoint
    Point2D<Real> const &
    Pb() const noexcept
    {
      return m_Pb;
    }

    //! Set first endpoint
    void
    Pa( Point2D<Real> const & p ) noexcept
    {
      m_Pa = p;
    }

    //! Set second endpoint
    void
    Pb( Point2D<Real> const & p ) noexcept
    {
      m_Pb = p;
    }

    //! Get midpoint of segment
    Point2D<Real>
    midpoint() const noexcept
    {
      return ( m_Pa + m_Pb ) * Real( 0.5 );
    }

    //! Compute bounding box
    void
    bbox( Point2D<Real> & pmin, Point2D<Real> & pmax ) const noexcept
    {
      pmin.to_eigen() = m_Pa.cwiseMin( m_Pb );
      pmax.to_eigen() = m_Pa.cwiseMax( m_Pb );
    }

    //! Get bounding box as Box2D
    Box2D<Real>
    bounding_box() const noexcept
    {
      Point2D<Real> pmin, pmax;
      bbox( pmin, pmax );
      return Box2D<Real>( pmin, pmax );
    }

    //! Distance from point to segment
    Real
    distance( Point2D<Real> const & P ) const
    {
      Real s, t;
      projection( P, s, t );
      return std::abs( t );
    }

    //! Squared distance from point to segment
    Real
    distance_squared( Point2D<Real> const & P ) const
    {
      Real          s, t;
      Point2D<Real> proj = projection( P, s, t );
      return P.distance_squared( proj );
    }

    //! Closest point on segment to given point
    Point2D<Real>
    closest_point( Point2D<Real> const & P ) const
    {
      Real s, t;
      return projection( P, s, t );
    }

    //! Check if point lies on segment (within tolerance)
    bool
    contains( Point2D<Real> const & P, Real tol = DEFAULT_TOLERANCE ) const
    {
      Real s, t;
      projection( P, s, t );
      return std::abs( t ) <= tol && s >= -tol && s <= 1 + tol;
    }

    //!
    //! \brief Checks for intersection with another segment.
    //!
    //! \param S The other segment to check for intersection.
    //! \param s Parameter along current segment at intersection (output).
    //! \param t Parameter along other segment at intersection (output).
    //! \param tol Tolerance for collinearity check.
    //!
    //! \return `true` if segments intersect, `false` otherwise.
    //!
    bool
    intersect( Segment2D<Real> const & S, Real & s, Real & t, Real tol = DEFAULT_TOLERANCE ) const
    {
      s = t = 0;

      Point2D<Real> const D1 = m_Pb - m_Pa;
      Point2D<Real> const D2 = S.m_Pb - S.m_Pa;

      Real const len1 = D1.norm();
      Real const len2 = D2.norm();

      // Check for degenerate segments
      if ( len1 < machine_eps<Real>() || len2 < machine_eps<Real>() ) { return false; }

      // Check if collinear using cross product
      Real const cross            = D1.cross( D2 );
      Real const normalized_cross = cross / ( len1 * len2 );

      if ( std::abs( normalized_cross ) <= tol )
      {
        // Collinear case: check for overlap
        if ( this->projection( S.m_Pa, s ) ) return true;
        if ( this->projection( S.m_Pb, s ) ) return true;
        if ( S.projection( m_Pa, t ) )
        {
          s = 0;
          return true;
        }
        if ( S.projection( m_Pb, t ) )
        {
          s = 1;
          return true;
        }
        s = t = 0;
        return false;
      }

      // Regular intersection case
      // Solve: Pa + s*D1 = Sa + t*D2
      // [D1, -D2] * [s; t] = Sa - Pa
      Eigen::Matrix<Real, 2, 2> M;
      M.col( 0 ) = D1.to_eigen();
      M.col( 1 ) = -D2.to_eigen();

      Point2D<Real> const RHS = S.m_Pa - m_Pa;
      Point2D<Real> const result( M.colPivHouseholderQr().solve( RHS.to_eigen() ) );

      s = result.x();
      t = result.y();

      return s >= 0 && s <= 1 && t >= 0 && t <= 1;
    }

    //! Get intersection point with another segment
    bool
    intersection_point( Segment2D<Real> const & S, Point2D<Real> & point, Real tol = DEFAULT_TOLERANCE ) const
    {
      Real s, t;
      if ( intersect( S, s, t, tol ) )
      {
        point = eval( s );
        return true;
      }
      return false;
    }

    //! Check if segments intersect (simpler interface)
    bool
    intersects( Segment2D<Real> const & S, Real tol = DEFAULT_TOLERANCE ) const
    {
      Real s, t;
      return intersect( S, s, t, tol );
    }

    //! Check if segments are parallel
    bool
    is_parallel( Segment2D<Real> const & S, Real tol = DEFAULT_TOLERANCE ) const
    {
      Point2D<Real> D1    = direction();
      Point2D<Real> D2    = S.direction();
      Real          cross = D1.cross( D2 );
      Real          len1  = D1.norm();
      Real          len2  = D2.norm();
      if ( len1 < machine_eps<Real>() || len2 < machine_eps<Real>() )
      {
        return true;  // Degenerate segments are considered parallel
      }
      return std::abs( cross / ( len1 * len2 ) ) <= tol;
    }

    //! Check if segments are perpendicular
    bool
    is_perpendicular( Segment2D<Real> const & S, Real tol = DEFAULT_TOLERANCE ) const
    {
      Point2D<Real> D1   = direction();
      Point2D<Real> D2   = S.direction();
      Real          dot  = D1.dot( D2 );
      Real          len1 = D1.norm();
      Real          len2 = D2.norm();
      if ( len1 < machine_eps<Real>() || len2 < machine_eps<Real>() )
      {
        return false;  // Degenerate segments are not perpendicular
      }
      return std::abs( dot / ( len1 * len2 ) ) <= tol;
    }

    //! Reverse the segment (swap endpoints)
    void
    reverse() noexcept
    {
      std::swap( m_Pa, m_Pb );
    }

    //! Get reversed segment
    Segment2D<Real>
    reversed() const noexcept
    {
      return Segment2D<Real>( m_Pb, m_Pa );
    }

    //! Split segment at parameter s, returning two segments
    std::pair<Segment2D<Real>, Segment2D<Real>>
    split( Real s ) const
    {
      Point2D<Real> P = eval( s );
      return std::make_pair( Segment2D<Real>( m_Pa, P ), Segment2D<Real>( P, m_Pb ) );
    }

    //! Stream output operator
    friend std::ostream &
    operator<<( std::ostream & os, Segment2D<Real> const & seg )
    {
      os << "Segment[" << seg.Pa() << " -> " << seg.Pb() << "]";
      return os;
    }
  };

  /*\
   |   _____     _                   _      ____  ____
   |  |_   _| __(_) __ _ _ __   __ _| | ___|___ \|  _ \
   |    | || '__| |/ _` | '_ \ / _` | |/ _ \ __) | | | |
   |    | || |  | | (_| | | | | (_| | |  __// __/| |_| |
   |    |_||_|  |_|\__,_|_| |_|\__, |_|\___|_____|____/
   |                           |___/
  \*/

  //!
  //! \brief A class representing a 2D triangle defined by three vertices.
  //!
  //! The `Triangle2D` class defines a 2D triangle using three points,
  //! \f$ P_a \f$, \f$ P_b \f$, and \f$ P_c \f$, which represent the vertices of
  //! the triangle.
  //!
  //! \tparam Real The data type used for the coordinates of the triangle's
  //! vertices.
  //!
  template <typename Real>
  class Triangle2D
  {
    Point2D<Real> m_Pa;  //!< First vertex
    Point2D<Real> m_Pb;  //!< Second vertex
    Point2D<Real> m_Pc;  //!< Third vertex

  public:
    //! Default constructor
    Triangle2D() = default;

    //! Constructor from three points
    Triangle2D( Point2D<Real> const & A, Point2D<Real> const & B, Point2D<Real> const & C ) noexcept
      : m_Pa( A ), m_Pb( B ), m_Pc( C )
    {
    }

    //! Constructor from coordinates
    Triangle2D( Real ax, Real ay, Real bx, Real by, Real cx, Real cy ) noexcept
      : m_Pa( ax, ay ), m_Pb( bx, by ), m_Pc( cx, cy )
    {
    }

    //! Setup triangle from three points
    void
    setup( Point2D<Real> const & A, Point2D<Real> const & B, Point2D<Real> const & C ) noexcept
    {
      m_Pa = A;
      m_Pb = B;
      m_Pc = C;
    }

    //! Setup triangle from coordinates
    void
    setup( Real ax, Real ay, Real bx, Real by, Real cx, Real cy ) noexcept
    {
      m_Pa.x( ax );
      m_Pa.y( ay );
      m_Pb.x( bx );
      m_Pb.y( by );
      m_Pc.x( cx );
      m_Pc.y( cy );
    }

    //! Get first vertex
    Point2D<Real> const &
    Pa() const noexcept
    {
      return m_Pa;
    }

    //! Get second vertex
    Point2D<Real> const &
    Pb() const noexcept
    {
      return m_Pb;
    }

    //! Get third vertex
    Point2D<Real> const &
    Pc() const noexcept
    {
      return m_Pc;
    }

    //! Set first vertex
    void
    Pa( Point2D<Real> const & p ) noexcept
    {
      m_Pa = p;
    }

    //! Set second vertex
    void
    Pb( Point2D<Real> const & p ) noexcept
    {
      m_Pb = p;
    }

    //! Set third vertex
    void
    Pc( Point2D<Real> const & p ) noexcept
    {
      m_Pc = p;
    }

    //! Compute signed area (positive if CCW, negative if CW)
    Real
    signed_area() const noexcept
    {
      return ( ( m_Pb.x() - m_Pa.x() ) * ( m_Pc.y() - m_Pa.y() ) - ( m_Pc.x() - m_Pa.x() ) * ( m_Pb.y() - m_Pa.y() ) ) *
             Real( 0.5 );
    }

    //! Compute area (always positive)
    Real
    area() const noexcept
    {
      return std::abs( signed_area() );
    }

    //! Check orientation (1 for CCW, -1 for CW, 0 for degenerate)
    int
    orientation( Real tol = Real( 1e-10 ) ) const noexcept
    {
      Real area = signed_area();
      if ( area > tol ) return 1;
      if ( area < -tol ) return -1;
      return 0;
    }

    //! Make triangle CCW oriented
    void
    make_ccw() noexcept
    {
      if ( orientation() < 0 ) { std::swap( m_Pb, m_Pc ); }
    }

    //! Check if triangle is valid (non-degenerate)
    bool
    is_valid( Real tol = Real( 1e-10 ) ) const noexcept
    {
      return area() > tol;
    }

    //! Compute centroid (center of mass)
    Point2D<Real>
    centroid() const noexcept
    {
      return Point2D<Real>( ( m_Pa.x() + m_Pb.x() + m_Pc.x() ) / Real( 3 ),
                            ( m_Pa.y() + m_Pb.y() + m_Pc.y() ) / Real( 3 ) );
    }

    //! Compute perimeter
    Real
    perimeter() const noexcept
    {
      return m_Pa.distance( m_Pb ) + m_Pb.distance( m_Pc ) + m_Pc.distance( m_Pa );
    }

    //! Get edge as segment
    Segment2D<Real>
    edge( int i ) const
    {
      i = ( ( i % 3 ) + 3 ) % 3;  // Handle negative indices
      switch ( i )
      {
        case 0:
          return Segment2D<Real>( m_Pa, m_Pb );
        case 1:
          return Segment2D<Real>( m_Pb, m_Pc );
        default:
          return Segment2D<Real>( m_Pc, m_Pa );
      }
    }

    //! Get all edges
    std::array<Segment2D<Real>, 3>
    edges() const noexcept
    {
      return { { Segment2D<Real>( m_Pa, m_Pb ), Segment2D<Real>( m_Pb, m_Pc ), Segment2D<Real>( m_Pc, m_Pa ) } };
    }

    //! Get vertex by index
    Point2D<Real>
    vertex( int i ) const
    {
      i = ( ( i % 3 ) + 3 ) % 3;  // Handle negative indices
      switch ( i )
      {
        case 0:
          return m_Pa;
        case 1:
          return m_Pb;
        default:
          return m_Pc;
      }
    }

    //! Check if point is inside triangle using barycentric coordinates
    bool
    contains( Point2D<Real> const & P, Real tol = Real( 1e-10 ) ) const
    {
      Real const total_area = signed_area();
      if ( std::abs( total_area ) < tol ) return false;

      // Compute signed areas of sub-triangles
      Real const area1 = ( ( P.x() - m_Pb.x() ) * ( m_Pc.y() - m_Pb.y() ) -
                           ( P.y() - m_Pb.y() ) * ( m_Pc.x() - m_Pb.x() ) ) *
                         Real( 0.5 );
      Real const area2 = ( ( P.x() - m_Pc.x() ) * ( m_Pa.y() - m_Pc.y() ) -
                           ( P.y() - m_Pc.y() ) * ( m_Pa.x() - m_Pc.x() ) ) *
                         Real( 0.5 );
      Real const area3 = ( ( P.x() - m_Pa.x() ) * ( m_Pb.y() - m_Pa.y() ) -
                           ( P.y() - m_Pa.y() ) * ( m_Pb.x() - m_Pa.x() ) ) *
                         Real( 0.5 );

      // Barycentric coordinates
      Real const u = area1 / total_area;
      Real const v = area2 / total_area;
      Real const w = area3 / total_area;

      // Check if point is inside (including edges with tolerance)
      return u >= -tol && v >= -tol && w >= -tol;
    }

    //! Check if point is strictly inside triangle (excluding edges)
    bool
    contains_strict( Point2D<Real> const & P, Real tol = Real( 1e-10 ) ) const
    {
      Real const total_area = signed_area();
      if ( std::abs( total_area ) < tol ) return false;

      Real const area1 = ( ( P.x() - m_Pb.x() ) * ( m_Pc.y() - m_Pb.y() ) -
                           ( P.y() - m_Pb.y() ) * ( m_Pc.x() - m_Pb.x() ) ) *
                         Real( 0.5 );
      Real const area2 = ( ( P.x() - m_Pc.x() ) * ( m_Pa.y() - m_Pc.y() ) -
                           ( P.y() - m_Pc.y() ) * ( m_Pa.x() - m_Pc.x() ) ) *
                         Real( 0.5 );
      Real const area3 = ( ( P.x() - m_Pa.x() ) * ( m_Pb.y() - m_Pa.y() ) -
                           ( P.y() - m_Pa.y() ) * ( m_Pb.x() - m_Pa.x() ) ) *
                         Real( 0.5 );

      Real const u = area1 / total_area;
      Real const v = area2 / total_area;
      Real const w = area3 / total_area;

      return u > tol && v > tol && w > tol;
    }

    //! Get barycentric coordinates of point
    bool
    barycentric_coords( Point2D<Real> const & P, Real & u, Real & v, Real & w, Real tol = Real( 1e-10 ) ) const
    {
      Real const total_area = signed_area();
      if ( std::abs( total_area ) < tol ) return false;

      Real const area1 = ( ( P.x() - m_Pb.x() ) * ( m_Pc.y() - m_Pb.y() ) -
                           ( P.y() - m_Pb.y() ) * ( m_Pc.x() - m_Pb.x() ) ) *
                         Real( 0.5 );
      Real const area2 = ( ( P.x() - m_Pc.x() ) * ( m_Pa.y() - m_Pc.y() ) -
                           ( P.y() - m_Pc.y() ) * ( m_Pa.x() - m_Pc.x() ) ) *
                         Real( 0.5 );
      Real const area3 = ( ( P.x() - m_Pa.x() ) * ( m_Pb.y() - m_Pa.y() ) -
                           ( P.y() - m_Pa.y() ) * ( m_Pb.x() - m_Pa.x() ) ) *
                         Real( 0.5 );

      u = area1 / total_area;
      v = area2 / total_area;
      w = area3 / total_area;

      return true;
    }

    //! Compute point from barycentric coordinates
    Point2D<Real>
    from_barycentric( Real u, Real v, Real w ) const noexcept
    {
      return m_Pa * u + m_Pb * v + m_Pc * w;
    }

    //! Compute bounding box
    void
    bbox( Point2D<Real> & pmin, Point2D<Real> & pmax ) const noexcept
    {
      pmin.x( std::min( { m_Pa.x(), m_Pb.x(), m_Pc.x() } ) );
      pmin.y( std::min( { m_Pa.y(), m_Pb.y(), m_Pc.y() } ) );
      pmax.x( std::max( { m_Pa.x(), m_Pb.x(), m_Pc.x() } ) );
      pmax.y( std::max( { m_Pa.y(), m_Pb.y(), m_Pc.y() } ) );
    }

    //! Get bounding box as Box2D
    Box2D<Real>
    bounding_box() const noexcept
    {
      Point2D<Real> pmin, pmax;
      bbox( pmin, pmax );
      return Box2D<Real>( pmin, pmax );
    }

    //! Check if triangle is acute (all angles < 90°)
    bool
    is_acute() const
    {
      // Compute squared side lengths
      Real a2 = m_Pb.distance_squared( m_Pc );
      Real b2 = m_Pc.distance_squared( m_Pa );
      Real c2 = m_Pa.distance_squared( m_Pb );

      // Check if all angles are acute
      return ( b2 + c2 > a2 ) && ( a2 + c2 > b2 ) && ( a2 + b2 > c2 );
    }

    //! Check if triangle is obtuse (one angle > 90°)
    bool
    is_obtuse() const
    {
      // Compute squared side lengths
      Real a2 = m_Pb.distance_squared( m_Pc );
      Real b2 = m_Pc.distance_squared( m_Pa );
      Real c2 = m_Pa.distance_squared( m_Pb );

      // Check if any angle is obtuse
      return ( b2 + c2 < a2 ) || ( a2 + c2 < b2 ) || ( a2 + b2 < c2 );
    }

    //! Check if triangle is right-angled (one angle = 90°)
    bool
    is_right_angled( Real tol = Real( 1e-6 ) ) const
    {
      // Compute squared side lengths
      Real a2 = m_Pb.distance_squared( m_Pc );
      Real b2 = m_Pc.distance_squared( m_Pa );
      Real c2 = m_Pa.distance_squared( m_Pb );

      // Check Pythagorean theorem within tolerance
      return std::abs( b2 + c2 - a2 ) < tol || std::abs( a2 + c2 - b2 ) < tol || std::abs( a2 + b2 - c2 ) < tol;
    }

    //! Check if triangle is equilateral (all sides equal)
    bool
    is_equilateral( Real tol = Real( 1e-6 ) ) const
    {
      Real a = m_Pb.distance( m_Pc );
      Real b = m_Pc.distance( m_Pa );
      Real c = m_Pa.distance( m_Pb );

      return std::abs( a - b ) < tol && std::abs( b - c ) < tol;
    }

    //! Check if triangle is isosceles (at least two sides equal)
    bool
    is_isosceles( Real tol = Real( 1e-6 ) ) const
    {
      Real a = m_Pb.distance( m_Pc );
      Real b = m_Pc.distance( m_Pa );
      Real c = m_Pa.distance( m_Pb );

      return std::abs( a - b ) < tol || std::abs( b - c ) < tol || std::abs( c - a ) < tol;
    }

    //! Compute circumcircle center and radius
    bool
    circumcircle( Point2D<Real> & center, Real & radius ) const
    {
      Real const D = 2 * ( m_Pa.x() * ( m_Pb.y() - m_Pc.y() ) + m_Pb.x() * ( m_Pc.y() - m_Pa.y() ) +
                           m_Pc.x() * ( m_Pa.y() - m_Pb.y() ) );

      if ( std::abs( D ) < machine_eps<Real>() )
      {
        return false;  // Collinear points
      }

      Real const a2 = m_Pa.x() * m_Pa.x() + m_Pa.y() * m_Pa.y();
      Real const b2 = m_Pb.x() * m_Pb.x() + m_Pb.y() * m_Pb.y();
      Real const c2 = m_Pc.x() * m_Pc.x() + m_Pc.y() * m_Pc.y();

      center.x( ( a2 * ( m_Pb.y() - m_Pc.y() ) + b2 * ( m_Pc.y() - m_Pa.y() ) + c2 * ( m_Pa.y() - m_Pb.y() ) ) / D );
      center.y( ( a2 * ( m_Pc.x() - m_Pb.x() ) + b2 * ( m_Pa.x() - m_Pc.x() ) + c2 * ( m_Pb.x() - m_Pa.x() ) ) / D );

      radius = center.distance( m_Pa );
      return true;
    }

    //! Compute incircle center and radius
    bool
    incircle( Point2D<Real> & center, Real & radius ) const
    {
      Real const a         = m_Pb.distance( m_Pc );
      Real const b         = m_Pc.distance( m_Pa );
      Real const c         = m_Pa.distance( m_Pb );
      Real const perimeter = a + b + c;

      if ( perimeter < machine_eps<Real>() )
      {
        return false;  // Degenerate triangle
      }

      center.x( ( a * m_Pa.x() + b * m_Pb.x() + c * m_Pc.x() ) / perimeter );
      center.y( ( a * m_Pa.y() + b * m_Pb.y() + c * m_Pc.y() ) / perimeter );

      Real const area = this->area();
      radius          = 2 * area / perimeter;
      return true;
    }

    //! Translate triangle by a vector
    void
    translate( Point2D<Real> const & delta ) noexcept
    {
      m_Pa = m_Pa + delta;
      m_Pb = m_Pb + delta;
      m_Pc = m_Pc + delta;
    }

    //! Translate triangle by coordinates
    void
    translate( Real dx, Real dy ) noexcept
    {
      m_Pa.x( m_Pa.x() + dx );
      m_Pa.y( m_Pa.y() + dy );
      m_Pb.x( m_Pb.x() + dx );
      m_Pb.y( m_Pb.y() + dy );
      m_Pc.x( m_Pc.x() + dx );
      m_Pc.y( m_Pc.y() + dy );
    }

    //! Rotate triangle by angle (radians) around origin
    void
    rotate( Real angle ) noexcept
    {
      m_Pa = m_Pa.rotate( angle );
      m_Pb = m_Pb.rotate( angle );
      m_Pc = m_Pc.rotate( angle );
    }

    //! Rotate triangle by angle around a center point
    void
    rotate( Real angle, Point2D<Real> const & center ) noexcept
    {
      m_Pa = m_Pa.rotate( angle, center );
      m_Pb = m_Pb.rotate( angle, center );
      m_Pc = m_Pc.rotate( angle, center );
    }

    //! Scale triangle relative to origin
    void
    scale( Real sx, Real sy ) noexcept
    {
      m_Pa = m_Pa.scale( sx, sy );
      m_Pb = m_Pb.scale( sx, sy );
      m_Pc = m_Pc.scale( sx, sy );
    }

    //! Scale triangle relative to a center point
    void
    scale( Real sx, Real sy, Point2D<Real> const & center ) noexcept
    {
      m_Pa = m_Pa.scale( sx, sy, center );
      m_Pb = m_Pb.scale( sx, sy, center );
      m_Pc = m_Pc.scale( sx, sy, center );
    }

    //! Stream output operator
    friend std::ostream &
    operator<<( std::ostream & os, Triangle2D<Real> const & tri )
    {
      os << "Triangle[" << tri.Pa() << ", " << tri.Pb() << ", " << tri.Pc() << "]";
      return os;
    }
  };

  /*\
   |   ____       _                         ____  ____
   |  |  _ \ ___ | |_   _  __ _  ___  _ __ |___ \|  _ \
   |  | |_) / _ \| | | | |/ _` |/ _ \| '_ \  __) | | | |
   |  |  __/ (_) | | |_| | (_| | (_) | | | |/ __/| |_| |
   |  |_|   \___/|_|\__, |\__, |\___/|_| |_|_____|____/
   |                |___/ |___/
  \*/

  //!
  //! \brief A class representing a 2D polygon using a dynamic matrix of
  //! vertices.
  //!
  //! The `Polygon2D` class defines a 2D polygon by storing its vertices in a
  //! dynamic matrix. Each column of the matrix represents a vertex of the
  //! polygon in 2D space.
  //!
  //! \tparam Real The data type used for the coordinates of the polygon's
  //! vertices.
  //!
  template <typename Real>
  class Polygon2D : public Eigen::Matrix<Real, 2, Eigen::Dynamic>
  {
    using BaseMatrix = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

  public:
    using BaseMatrix::Matrix;  // Inherit constructors

    //! Default constructor
    Polygon2D() : BaseMatrix() {}

    //! Constructor with number of vertices
    explicit Polygon2D( int num_vertices ) : BaseMatrix( 2, num_vertices ) {}

    //! Constructor from vector of points
    explicit Polygon2D( std::vector<Point2D<Real>> const & points )
    {
      resize( points.size() );
      for ( size_t i = 0; i < points.size(); ++i ) { set_vertex( static_cast<int>( i ), points[i] ); }
    }

    //! Constructor from initializer list
    Polygon2D( std::initializer_list<Point2D<Real>> points )
    {
      resize( static_cast<int>( points.size() ) );
      int i = 0;
      for ( auto const & p : points ) { set_vertex( i++, p ); }
    }

    //! Get number of vertices
    int
    num_vertices() const noexcept
    {
      return static_cast<int>( this->cols() );
    }

    //! Check if polygon is empty
    bool
    empty() const noexcept
    {
      return num_vertices() == 0;
    }

    //! Resize polygon to have n vertices
    void
    resize( int n )
    {
      BaseMatrix::resize( 2, n );
    }

    //! Clear polygon (remove all vertices)
    void
    clear()
    {
      BaseMatrix::resize( 2, 0 );
    }

    //! Get vertex at index i (with bounds checking)
    Point2D<Real>
    vertex( int i ) const
    {
      i = ( ( i % num_vertices() ) + num_vertices() ) % num_vertices();  // Wrap index
      return Point2D<Real>( this->coeff( 0, i ), this->coeff( 1, i ) );
    }

    //! Set vertex at index i
    void
    set_vertex( int i, Point2D<Real> const & p )
    {
      if ( i >= 0 && i < num_vertices() )
      {
        this->coeffRef( 0, i ) = p.x();
        this->coeffRef( 1, i ) = p.y();
      }
    }

    //! Add vertex to the end of polygon
    void
    add_vertex( Point2D<Real> const & p )
    {
      int n = num_vertices();
      this->conservativeResize( 2, n + 1 );
      set_vertex( n, p );
    }

    //! Add vertex by coordinates
    void
    add_vertex( Real x, Real y )
    {
      add_vertex( Point2D<Real>( x, y ) );
    }

    //! Insert vertex at position i
    void
    insert_vertex( int i, Point2D<Real> const & p )
    {
      if ( i < 0 || i > num_vertices() ) return;

      int n = num_vertices();
      this->conservativeResize( 2, n + 1 );

      // Shift vertices to the right
      for ( int j = n; j > i; --j ) { set_vertex( j, vertex( j - 1 ) ); }

      set_vertex( i, p );
    }

    //! Remove vertex at position i
    void
    remove_vertex( int i )
    {
      if ( i < 0 || i >= num_vertices() || num_vertices() <= 0 ) return;

      int n = num_vertices();
      for ( int j = i; j < n - 1; ++j ) { set_vertex( j, vertex( j + 1 ) ); }

      this->conservativeResize( 2, n - 1 );
    }

    //! Get edge as segment
    Segment2D<Real>
    edge( int i ) const
    {
      return Segment2D<Real>( vertex( i ), vertex( i + 1 ) );
    }

    //! Get all edges
    std::vector<Segment2D<Real>>
    edges() const
    {
      std::vector<Segment2D<Real>> result;
      int                          n = num_vertices();
      if ( n < 2 ) return result;

      result.reserve( n );
      for ( int i = 0; i < n; ++i ) { result.push_back( edge( i ) ); }
      return result;
    }

    //! Compute centroid
    Point2D<Real>
    centroid() const
    {
      int n = num_vertices();
      if ( n == 0 ) return Point2D<Real>( 0, 0 );
      if ( n == 1 ) return vertex( 0 );

      Point2D<Real> c( 0, 0 );
      for ( int i = 0; i < n; ++i )
      {
        c.x( c.x() + vertex( i ).x() );
        c.y( c.y() + vertex( i ).y() );
      }
      return Point2D<Real>( c.x() / n, c.y() / n );
    }

    //! Compute signed area (positive if CCW, negative if CW)
    Real
    signed_area() const noexcept
    {
      int n = num_vertices();
      if ( n < 3 ) return 0;

      Real area = 0;
      for ( int i = 0; i < n; ++i )
      {
        Point2D<Real> const & p1 = vertex( i );
        Point2D<Real> const & p2 = vertex( ( i + 1 ) % n );
        area += p1.x() * p2.y() - p2.x() * p1.y();
      }
      return area * Real( 0.5 );
    }

    //! Compute area (always positive)
    Real
    area() const noexcept
    {
      return std::abs( signed_area() );
    }

    //! Check orientation (1 for CCW, -1 for CW, 0 for degenerate)
    int
    orientation( Real tol = Real( 1e-10 ) ) const noexcept
    {
      Real area = signed_area();
      if ( area > tol ) return 1;
      if ( area < -tol ) return -1;
      return 0;
    }

    //! Make polygon CCW oriented
    void
    make_ccw() noexcept
    {
      if ( orientation() < 0 ) { reverse(); }
    }

    //! Make polygon CW oriented
    void
    make_cw() noexcept
    {
      if ( orientation() > 0 ) { reverse(); }
    }

    //! Reverse vertex order
    void
    reverse() noexcept
    {
      int n = num_vertices();
      for ( int i = 0; i < n / 2; ++i )
      {
        Point2D<Real> temp = vertex( i );
        set_vertex( i, vertex( n - 1 - i ) );
        set_vertex( n - 1 - i, temp );
      }
    }

    //! Compute perimeter
    Real
    perimeter() const noexcept
    {
      int n = num_vertices();
      if ( n < 2 ) return 0;

      Real p = 0;
      for ( int i = 0; i < n; ++i ) { p += vertex( i ).distance( vertex( ( i + 1 ) % n ) ); }
      return p;
    }

    //! Check if polygon is convex
    bool
    is_convex( Real tol = Real( 1e-10 ) ) const
    {
      int n = num_vertices();
      if ( n < 3 ) return false;

      // Check cross product signs
      int sign = 0;
      for ( int i = 0; i < n; ++i )
      {
        Point2D<Real> const & p0 = vertex( i );
        Point2D<Real> const & p1 = vertex( ( i + 1 ) % n );
        Point2D<Real> const & p2 = vertex( ( i + 2 ) % n );

        Real cross = ( p1.x() - p0.x() ) * ( p2.y() - p1.y() ) - ( p1.y() - p0.y() ) * ( p2.x() - p1.x() );

        if ( std::abs( cross ) > tol )
        {
          int cur_sign = ( cross > 0 ) ? 1 : -1;
          if ( sign == 0 ) { sign = cur_sign; }
          else if ( sign != cur_sign ) { return false; }
        }
      }
      return true;
    }

    //! Check if polygon is simple (non-self-intersecting)
    bool
    is_simple() const
    {
      int n = num_vertices();
      if ( n < 3 ) return true;

      // Brute-force check for edge intersections
      auto edges = this->edges();
      for ( int i = 0; i < n; ++i )
      {
        for ( int j = i + 1; j < n; ++j )
        {
          // Skip adjacent edges
          if ( j == i + 1 || ( i == 0 && j == n - 1 ) ) continue;

          if ( edges[i].intersects( edges[j] ) ) { return false; }
        }
      }
      return true;
    }

    //! Check if point is inside polygon using winding number
    bool
    contains( Point2D<Real> const & P ) const
    {
      int n = num_vertices();
      if ( n < 3 ) return false;

      int wn = 0;  // Winding number
      for ( int i = 0; i < n; ++i )
      {
        Point2D<Real> const & Vi = vertex( i );
        Point2D<Real> const & Vj = vertex( ( i + 1 ) % n );

        if ( Vi.y() <= P.y() )
        {
          if ( Vj.y() > P.y() )
          {
            // Upward crossing
            Real is_left = ( Vj.x() - Vi.x() ) * ( P.y() - Vi.y() ) - ( P.x() - Vi.x() ) * ( Vj.y() - Vi.y() );
            if ( is_left > 0 ) { ++wn; }
          }
        }
        else
        {
          if ( Vj.y() <= P.y() )
          {
            // Downward crossing
            Real is_left = ( Vj.x() - Vi.x() ) * ( P.y() - Vi.y() ) - ( P.x() - Vi.x() ) * ( Vj.y() - Vi.y() );
            if ( is_left < 0 ) { --wn; }
          }
        }
      }
      return wn != 0;
    }

    //! Check if point is on polygon boundary
    bool
    on_boundary( Point2D<Real> const & P, Real tol = Real( 1e-6 ) ) const
    {
      int n = num_vertices();
      for ( int i = 0; i < n; ++i )
      {
        Segment2D<Real> edge( vertex( i ), vertex( ( i + 1 ) % n ) );
        if ( edge.contains( P, tol ) ) { return true; }
      }
      return false;
    }

    //! Compute bounding box
    void
    bbox( Point2D<Real> & pmin, Point2D<Real> & pmax ) const
    {
      int n = num_vertices();
      if ( n == 0 )
      {
        pmin = pmax = Point2D<Real>( 0, 0 );
        return;
      }

      pmin = vertex( 0 );
      pmax = vertex( 0 );

      for ( int i = 1; i < n; ++i )
      {
        Point2D<Real> const & p = vertex( i );
        pmin.x( std::min( pmin.x(), p.x() ) );
        pmin.y( std::min( pmin.y(), p.y() ) );
        pmax.x( std::max( pmax.x(), p.x() ) );
        pmax.y( std::max( pmax.y(), p.y() ) );
      }
    }

    //! Get bounding box as Box2D
    Box2D<Real>
    bounding_box() const
    {
      Point2D<Real> pmin, pmax;
      bbox( pmin, pmax );
      return Box2D<Real>( pmin, pmax );
    }

    //! Translate polygon by a vector
    void
    translate( Point2D<Real> const & delta ) noexcept
    {
      for ( int i = 0; i < num_vertices(); ++i )
      {
        Point2D<Real> p = vertex( i );
        p               = p + delta;
        set_vertex( i, p );
      }
    }

    //! Translate polygon by coordinates
    void
    translate( Real dx, Real dy ) noexcept
    {
      translate( Point2D<Real>( dx, dy ) );
    }

    //! Rotate polygon by angle (radians) around origin
    void
    rotate( Real angle ) noexcept
    {
      for ( int i = 0; i < num_vertices(); ++i )
      {
        Point2D<Real> p = vertex( i );
        p               = p.rotate( angle );
        set_vertex( i, p );
      }
    }

    //! Rotate polygon by angle around a center point
    void
    rotate( Real angle, Point2D<Real> const & center ) noexcept
    {
      for ( int i = 0; i < num_vertices(); ++i )
      {
        Point2D<Real> p = vertex( i );
        p               = p.rotate( angle, center );
        set_vertex( i, p );
      }
    }

    //! Scale polygon relative to origin
    void
    scale( Real sx, Real sy ) noexcept
    {
      for ( int i = 0; i < num_vertices(); ++i )
      {
        Point2D<Real> p = vertex( i );
        p               = p.scale( sx, sy );
        set_vertex( i, p );
      }
    }

    //! Scale polygon relative to a center point
    void
    scale( Real sx, Real sy, Point2D<Real> const & center ) noexcept
    {
      for ( int i = 0; i < num_vertices(); ++i )
      {
        Point2D<Real> p = vertex( i );
        p               = p.scale( sx, sy, center );
        set_vertex( i, p );
      }
    }

    //! Get polygon as vector of points
    std::vector<Point2D<Real>>
    to_vector() const
    {
      std::vector<Point2D<Real>> result;
      int                        n = num_vertices();
      result.reserve( n );
      for ( int i = 0; i < n; ++i ) { result.push_back( vertex( i ) ); }
      return result;
    }

    //! Simplify polygon (remove collinear points)
    void
    simplify( Real tol = Real( 1e-6 ) )
    {
      int n = num_vertices();
      if ( n < 3 ) return;

      std::vector<bool> keep( n, true );
      for ( int i = 0; i < n; ++i )
      {
        int prev = ( i - 1 + n ) % n;
        int next = ( i + 1 ) % n;

        Point2D<Real> const & A = vertex( prev );
        Point2D<Real> const & B = vertex( i );
        Point2D<Real> const & C = vertex( next );

        // Check if area of triangle ABC is small
        Real area = std::abs( ( B.x() - A.x() ) * ( C.y() - A.y() ) - ( C.x() - A.x() ) * ( B.y() - A.y() ) ) *
                    Real( 0.5 );

        if ( area < tol ) { keep[i] = false; }
      }

      // Create new polygon
      Polygon2D<Real> new_poly;
      for ( int i = 0; i < n; ++i )
      {
        if ( keep[i] ) { new_poly.add_vertex( vertex( i ) ); }
      }

      // Replace current polygon
      *this = new_poly;
    }

    //! Stream output operator
    friend std::ostream &
    operator<<( std::ostream & os, Polygon2D<Real> const & poly )
    {
      os << "Polygon[";
      int n = poly.num_vertices();
      for ( int i = 0; i < n; ++i )
      {
        os << poly.vertex( i );
        if ( i < n - 1 ) os << ", ";
      }
      os << "]";
      return os;
    }
  };

  /*! @} */

}  // namespace Utils

#endif

//
// eof: Utils_GG2D.hh
//
