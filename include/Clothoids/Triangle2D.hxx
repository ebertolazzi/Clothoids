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
\*--------------------------------------------------------------------------*/

///
/// \file Triangle2D.hxx
/// \brief Definition of the Triangle2D class for bounding box computations
///        in clothoid curve algorithms.
///

namespace G2lib
{

  /*\
   |   _____     _                   _      ____  ____
   |  |_   _| __(_) __ _ _ __   __ _| | ___|___ \|  _ \
   |    | || '__| |/ _` | '_ \ / _` | |/ _ \ __) | | | |
   |    | || |  | | (_| | | | | (_| | |  __// __/| |_| |
   |    |_||_|  |_|\__,_|_| |_|\__, |_|\___|_____|____/
   |                           |___/
  \*/

  //!
  //! \class Triangle2D
  //! \brief Manages a triangle in 2D space for bounding box computations of clothoid curves.
  //!
  //! This class represents a triangle defined by three vertices (P1, P2, P3) in 2D space.
  //! It is used in the bounding box (BB) computations for clothoid curves. Each triangle
  //! stores additional information about the curve segment it represents, including
  //! the parameter range [s0, s1] and the curve index.
  //!
  class Triangle2D
  {
    real_type m_p1[2];   ///< Coordinates of the first vertex (x1, y1)
    real_type m_p2[2];   ///< Coordinates of the second vertex (x2, y2)
    real_type m_p3[2];   ///< Coordinates of the third vertex (x3, y3)
    real_type m_s0;      ///< Starting parameter value of the curve segment represented by this triangle
    real_type m_s1;      ///< Ending parameter value of the curve segment represented by this triangle
    integer   m_icurve;  ///< Index of the curve segment in the original clothoid list

  public:
    ///@name Constructors and Destructor
    ///@{

    //! \brief Copy constructor.
    //! \param[in] t The triangle to copy.
    Triangle2D( Triangle2D const & t ) { *this = t; }

    //! \brief Default constructor. Creates a degenerate triangle at origin.
    Triangle2D()
    {
      m_p1[0]  = 0;
      m_p1[1]  = 0;
      m_p2[0]  = 0;
      m_p2[1]  = 0;
      m_p3[0]  = 0;
      m_p3[1]  = 0;
      m_s0     = 0;
      m_s1     = 0;
      m_icurve = 0;
    }

    //! \brief Construct a triangle from coordinates and curve information.
    //! \param[in] x1 X-coordinate of the first vertex.
    //! \param[in] y1 Y-coordinate of the first vertex.
    //! \param[in] x2 X-coordinate of the second vertex.
    //! \param[in] y2 Y-coordinate of the second vertex.
    //! \param[in] x3 X-coordinate of the third vertex.
    //! \param[in] y3 Y-coordinate of the third vertex.
    //! \param[in] s0 Starting parameter of the curve segment.
    //! \param[in] s1 Ending parameter of the curve segment.
    //! \param[in] icurve Index of the curve segment.
    Triangle2D(
      real_type const x1,
      real_type const y1,
      real_type const x2,
      real_type const y2,
      real_type const x3,
      real_type const y3,
      real_type const s0,
      real_type const s1,
      integer const   icurve )
    {
      m_p1[0]  = x1;
      m_p1[1]  = y1;
      m_p2[0]  = x2;
      m_p2[1]  = y2;
      m_p3[0]  = x3;
      m_p3[1]  = y3;
      m_s0     = s0;
      m_s1     = s1;
      m_icurve = icurve;
    }

    //! \brief Construct a triangle from point arrays and curve information.
    //! \param[in] p1 Array containing coordinates of the first vertex [x1, y1].
    //! \param[in] p2 Array containing coordinates of the second vertex [x2, y2].
    //! \param[in] p3 Array containing coordinates of the third vertex [x3, y3].
    //! \param[in] s0 Starting parameter of the curve segment.
    //! \param[in] s1 Ending parameter of the curve segment.
    //! \param[in] icurve Index of the curve segment.
    Triangle2D(
      real_type const p1[2],
      real_type const p2[2],
      real_type const p3[2],
      real_type const s0,
      real_type const s1,
      integer const   icurve )
    {
      m_p1[0]  = p1[0];
      m_p1[1]  = p1[1];
      m_p2[0]  = p2[0];
      m_p2[1]  = p2[1];
      m_p3[0]  = p3[0];
      m_p3[1]  = p3[1];
      m_s0     = s0;
      m_s1     = s1;
      m_icurve = icurve;
    }

    //! \brief Destructor.
    ~Triangle2D() {}
    ///@}

    ///@name Assignment Operator
    ///@{

    //! \brief Assignment operator.
    //! \param[in] t The triangle to assign.
    //! \return Reference to the assigned triangle.
    Triangle2D const & operator=( Triangle2D const & t )
    {
      m_p1[0]  = t.m_p1[0];
      m_p1[1]  = t.m_p1[1];
      m_p2[0]  = t.m_p2[0];
      m_p2[1]  = t.m_p2[1];
      m_p3[0]  = t.m_p3[0];
      m_p3[1]  = t.m_p3[1];
      m_s0     = t.m_s0;
      m_s1     = t.m_s1;
      m_icurve = t.m_icurve;
      return *this;
    }
    ///@}

    ///@name Building Methods
    ///@{

    //! \brief Build a triangle from point arrays and curve information.
    //! \param[in] p1 Array containing coordinates of the first vertex [x1, y1].
    //! \param[in] p2 Array containing coordinates of the second vertex [x2, y2].
    //! \param[in] p3 Array containing coordinates of the third vertex [x3, y3].
    //! \param[in] s0 Starting parameter of the curve segment.
    //! \param[in] s1 Ending parameter of the curve segment.
    //! \param[in] icurve Index of the curve segment.
    void build(
      real_type const p1[2],
      real_type const p2[2],
      real_type const p3[2],
      real_type const s0,
      real_type const s1,
      integer const   icurve )
    {
      m_p1[0]  = p1[0];
      m_p1[1]  = p1[1];
      m_p2[0]  = p2[0];
      m_p2[1]  = p2[1];
      m_p3[0]  = p3[0];
      m_p3[1]  = p3[1];
      m_s0     = s0;
      m_s1     = s1;
      m_icurve = icurve;
    }

    //! \brief Build a triangle from coordinates and curve information.
    //! \param[in] x1 X-coordinate of the first vertex.
    //! \param[in] y1 Y-coordinate of the first vertex.
    //! \param[in] x2 X-coordinate of the second vertex.
    //! \param[in] y2 Y-coordinate of the second vertex.
    //! \param[in] x3 X-coordinate of the third vertex.
    //! \param[in] y3 Y-coordinate of the third vertex.
    //! \param[in] s0 Starting parameter of the curve segment.
    //! \param[in] s1 Ending parameter of the curve segment.
    //! \param[in] icurve Index of the curve segment.
    void build(
      real_type const x1,
      real_type const y1,
      real_type const x2,
      real_type const y2,
      real_type const x3,
      real_type const y3,
      real_type const s0,
      real_type const s1,
      integer const   icurve )
    {
      m_p1[0]  = x1;
      m_p1[1]  = y1;
      m_p2[0]  = x2;
      m_p2[1]  = y2;
      m_p3[0]  = x3;
      m_p3[1]  = y3;
      m_s0     = s0;
      m_s1     = s1;
      m_icurve = icurve;
    }
    ///@}

    ///@name Access Methods
    ///@{

    //! \brief Get the curve index.
    //! \return The index of the curve segment represented by this triangle.
    integer Icurve() const { return m_icurve; }

    //! \brief Get the x-coordinate of the first vertex.
    //! \return X-coordinate of P1.
    real_type x1() const { return m_p1[0]; }

    //! \brief Get the y-coordinate of the first vertex.
    //! \return Y-coordinate of P1.
    real_type y1() const { return m_p1[1]; }

    //! \brief Get the x-coordinate of the second vertex.
    //! \return X-coordinate of P2.
    real_type x2() const { return m_p2[0]; }

    //! \brief Get the y-coordinate of the second vertex.
    //! \return Y-coordinate of P2.
    real_type y2() const { return m_p2[1]; }

    //! \brief Get the x-coordinate of the third vertex.
    //! \return X-coordinate of P3.
    real_type x3() const { return m_p3[0]; }

    //! \brief Get the y-coordinate of the third vertex.
    //! \return Y-coordinate of P3.
    real_type y3() const { return m_p3[1]; }

    //! \brief Get the starting parameter of the curve segment.
    //! \return The parameter s0.
    real_type S0() const { return m_s0; }

    //! \brief Get the ending parameter of the curve segment.
    //! \return The parameter s1.
    real_type S1() const { return m_s1; }

    //! \brief Get a pointer to the first vertex.
    //! \return Constant pointer to the array containing P1 coordinates.
    real_type const * P1() const { return m_p1; }

    //! \brief Get a pointer to the second vertex.
    //! \return Constant pointer to the array containing P2 coordinates.
    real_type const * P2() const { return m_p2; }

    //! \brief Get a pointer to the third vertex.
    //! \return Constant pointer to the array containing P3 coordinates.
    real_type const * P3() const { return m_p3; }
    ///@}

    ///@name Geometric Transformations
    ///@{

    //! \brief Translate the triangle by (tx, ty).
    //! \param[in] tx Translation in the x-direction.
    //! \param[in] ty Translation in the y-direction.
    void translate( real_type const tx, real_type const ty )
    {
      m_p1[0] += tx;
      m_p2[0] += tx;
      m_p3[0] += tx;
      m_p1[1] += ty;
      m_p2[1] += ty;
      m_p3[1] += ty;
    }

    //! \brief Rotate the triangle around a point (cx, cy) by given angle.
    //! \param[in] angle Rotation angle in radians.
    //! \param[in] cx X-coordinate of the rotation center.
    //! \param[in] cy Y-coordinate of the rotation center.
    void rotate( real_type const angle, real_type const cx, real_type const cy );

    //! \brief Scale the triangle by factor sc (uniform scaling).
    //! \param[in] sc Scaling factor.
    void scale( real_type const sc )
    {
      m_p1[0] *= sc;
      m_p1[1] *= sc;
      m_p2[0] *= sc;
      m_p2[1] *= sc;
      m_p3[0] *= sc;
      m_p3[1] *= sc;
    }
    ///@}

    ///@name Geometric Properties
    ///@{

    //! \brief Compute the bounding box of the triangle.
    //! \param[out] xmin Minimum x-coordinate of the bounding box.
    //! \param[out] ymin Minimum y-coordinate of the bounding box.
    //! \param[out] xmax Maximum x-coordinate of the bounding box.
    //! \param[out] ymax Maximum y-coordinate of the bounding box.
    void bbox( real_type & xmin, real_type & ymin, real_type & xmax, real_type & ymax ) const
    {
      minmax3( m_p1[0], m_p2[0], m_p3[0], xmin, xmax );
      minmax3( m_p1[1], m_p2[1], m_p3[1], ymin, ymax );
    }

    //! \brief Compute the x-coordinate of the barycenter (centroid) of the triangle.
    //! \return X-coordinate of the barycenter.
    real_type baricenter_x() const { return ( m_p1[0] + m_p2[0] + m_p3[0] ) / 3; }

    //! \brief Compute the y-coordinate of the barycenter (centroid) of the triangle.
    //! \return Y-coordinate of the barycenter.
    real_type baricenter_y() const { return ( m_p1[1] + m_p2[1] + m_p3[1] ) / 3; }

    //! \brief Compute the signed area of the triangle.
    //! \return The signed area (positive if vertices are counter-clockwise,
    //!         negative if clockwise, zero if degenerate).
    //!
    //! The area is computed using the shoelace formula:
    //! \f[ \text{Area} = \frac{1}{2} \left| x1(y2 - y3) + x2(y3 - y1) + x3(y1 - y2) \right| \f]
    real_type area() const
    {
      return 0.5 *
             std::abs( ( m_p2[0] - m_p1[0] ) * ( m_p3[1] - m_p1[1] ) - ( m_p3[0] - m_p1[0] ) * ( m_p2[1] - m_p1[1] ) );
    }

    //! \brief Check if this triangle overlaps with another triangle.
    //! \param[in] t The triangle to check for overlap.
    //! \return True if triangles overlap, false otherwise.
    bool overlap( Triangle2D const & t ) const;
    ///@}

    ///@name Orientation and Point Inclusion Tests
    ///@{

    //! \brief Check the orientation of the triangle vertices.
    //! \return +1 if vertices are in counter-clockwise order,
    //!         -1 if vertices are in clockwise order,
    //!          0 if the triangle is degenerate (collinear vertices).
    integer is_counter_clockwise() const { return G2lib::is_counter_clockwise( m_p1, m_p2, m_p3 ); }

    //! \brief Test if a point (x, y) is inside the triangle.
    //! \param[in] x X-coordinate of the point.
    //! \param[in] y Y-coordinate of the point.
    //! \return +1 if the point is strictly inside the triangle,
    //!         -1 if the point is strictly outside,
    //!          0 if the point lies on the border (within numerical tolerance).
    integer is_inside( real_type const x, real_type const y ) const
    {
      real_type const pt[2] = { x, y };
      return is_point_in_triangle( pt, m_p1, m_p2, m_p3 );
    }

    //! \brief Test if a point given as an array is inside the triangle.
    //! \param[in] pt Array containing the point coordinates [x, y].
    //! \return +1 if the point is strictly inside the triangle,
    //!         -1 if the point is strictly outside,
    //!          0 if the point lies on the border (within numerical tolerance).
    integer is_inside( real_type const pt[2] ) const { return is_point_in_triangle( pt, m_p1, m_p2, m_p3 ); }
    ///@}

    ///@name Distance Computations
    ///@{

    //! \brief Compute the minimum distance from a point to the triangle (to vertices or edges).
    //! \param[in] x X-coordinate of the point.
    //! \param[in] y Y-coordinate of the point.
    //! \return The minimum Euclidean distance from the point to the triangle.
    real_type dist_min( real_type const x, real_type const y ) const;

    //! \brief Compute the maximum distance from a point to the triangle (to the farthest vertex).
    //! \param[in] x X-coordinate of the point.
    //! \param[in] y Y-coordinate of the point.
    //! \return The maximum Euclidean distance from the point to the triangle vertices.
    real_type dist_max( real_type const x, real_type const y ) const;
    ///@}

    ///@name Information and Output
    ///@{

    //! \brief Get a string with information about the triangle.
    //! \return A formatted string containing triangle vertices, curve index, and parameter range.
    string info() const;

    //! \brief Write triangle information to an output stream.
    //! \param[in] stream The output stream to write to.
    void info( ostream_type & stream ) const { stream << this->info(); }

    //! \brief Output operator for Triangle2D.
    //! \param[in] stream The output stream.
    //! \param[in] t The triangle to output.
    //! \return Reference to the output stream.
    friend ostream_type & operator<<( ostream_type & stream, Triangle2D const & t )
    {
      t.info( stream );
      return stream;
    }
    ///@}
  };

}  // namespace G2lib

///
/// \endofile Triangle2D.hxx
///
