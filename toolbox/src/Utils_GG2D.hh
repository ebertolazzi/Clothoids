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

//
// file: Utils_GG2D.hh
//
#pragma once

#ifndef UTILS_GG2D_dot_HH
#define UTILS_GG2D_dot_HH

#include "Utils_eigen.hh"

namespace Utils {

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
  //! The `Point2D` class extends the `Eigen::Matrix` class to provide a convenient
  //! representation of a point in a two-dimensional space. It encapsulates
  //! operations related to 2D points, allowing easy access to the x and y coordinates.
  //! This class leverages the Eigen library for efficient matrix and vector computations.
  //!
  //! \tparam Real The data type used for the point coordinates, typically a floating-point type (e.g., float, double).
  //!
  //! **Features**
  //!
  //! - **Coordinate Access**: Provides methods to access the x and y coordinates of the point easily.
  //! - **Eigen Compatibility**: Supports seamless conversion to and from Eigen's matrix representation, allowing
  //!   users to leverage Eigen's functionalities directly with the `Point2D` class.
  //!
  //! **Methods**
  //!
  //! - **Default Constructor**:
  //!   - Initializes a `Point2D` object with default values.
  //!
  //! - **Real x() const**:
  //!   - Returns the x-coordinate of the point.
  //!
  //! - **Real y() const**:
  //!   - Returns the y-coordinate of the point.
  //!
  //! - **P2D const & to_eigen() const**:
  //!   - Returns a constant reference to the underlying Eigen matrix representation of the point.
  //!
  //! - **P2D & to_eigen()**:
  //!   - Returns a reference to the underlying Eigen matrix representation of the point, allowing for modifications.
  //!
  //! **Usage**
  //!
  //! \code
  //! #include <Eigen/Dense>
  //!
  //! // Create a Point2D object
  //! Point2D<double> point;
  //!
  //! // Set coordinates
  //! point.to_eigen() << 3.0, 4.0;
  //!
  //! // Access coordinates
  //! double x = point.x(); // x = 3.0
  //! double y = point.y(); // y = 4.0
  //!
  //! // Convert to Eigen matrix
  //! Eigen::Matrix<double, 2, 1> eigen_point = point.to_eigen();
  //! \endcode
  //!
  template <typename Real>
  class Point2D : public Eigen::Matrix<Real,2,1> {
    using P2D = Eigen::Matrix<Real,2,1>;
  public:
    Point2D() = default;

    Real x() const { return this->coeff(0); }
    Real y() const { return this->coeff(1); }

    P2D const & to_eigen() const { return *static_cast<P2D const *>(this); }
    P2D &       to_eigen()       { return *static_cast<P2D*>(this); }
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
  //! The `Segment2D` class encapsulates a line segment in a two-dimensional space,
  //! defined by two points, \f$ P_a \f$ and \f$ P_b \f$. It provides methods for setting
  //! the endpoints, projecting points onto the segment, evaluating points on the segment,
  //! and checking for intersections with other segments.
  //!
  //! \tparam Real The data type used for the coordinates of the segment's endpoints,
  //! typically a floating-point type (e.g., float, double).
  //!
  //! **Features**
  //!
  //! - **Construction**: Supports various constructors for initializing the segment.
  //! - **Setup**: Methods to set or update the endpoints of the segment.
  //! - **Projection**: Projects a point onto the segment and returns the corresponding scalar parameter.
  //! - **Evaluation**: Computes points on the segment based on parameter values.
  //! - **Bounding Box**: Computes the bounding box of the segment.
  //! - **Intersection**: Checks for intersection with another segment and computes parameters.
  //!
  //! **Methods**
  //!
  //! - **Segment2D()**:
  //!   - Default constructor that initializes an empty segment.
  //!
  //! - **Segment2D(Segment2D<Real> const & S)**:
  //!   - Copy constructor that initializes a segment as a copy of another segment.
  //!
  //! - **Segment2D(Point2D<Real> const & A, Point2D<Real> const & B)**:
  //!   - Constructs a segment from two given points A and B.
  //!
  //! - **void setup(Point2D<Real> const & A, Point2D<Real> const & B)**:
  //!   - Sets the endpoints of the segment to the given points A and B.
  //!
  //! - **void setup(Real const A[], Real const B[])**:
  //!   - Sets the endpoints of the segment using arrays representing the coordinates of A and B.
  //!
  //! - **bool projection(Point2D<Real> const & P, Real & s) const**:
  //!   - Projects the point P onto the segment and returns the corresponding parameter `s`.
  //!   - \return `true` if the projection is successful, `false` otherwise.
  //!
  //! - **Point2D<Real> projection(Point2D<Real> const & P, Real & s, Real & t) const**:
  //!   - Projects the point P onto the segment and returns the projected point.
  //!   - Also returns the parameters `s` and `t` related to the projections on the segments.
  //!
  //! - **Point2D<Real> eval(Real & s) const**:
  //!   - Evaluates the segment at the parameter `s` and returns the corresponding point on the segment.
  //!
  //! - **Point2D<Real> eval(Real & s, Real & t) const**:
  //!   - Evaluates the segment based on parameters `s` and `t` and returns the corresponding point.
  //!
  //! - **Point2D<Real> const & Pa() const**:
  //!   - Returns the first endpoint \f$ P_a \f$ of the segment.
  //!
  //! - **Point2D<Real> const & Pb() const**:
  //!   - Returns the second endpoint \f$ P_b \f$ of the segment.
  //!
  //! - **void bbox(Point2D<Real> & pmin, Point2D<Real> & pmax) const**:
  //!   - Computes the bounding box of the segment and updates `pmin` and `pmax` with the minimum and maximum points.
  //!
  //! - **bool intersect(Segment2D<Real> const & S, Real & s, Real & t) const**:
  //!   - Checks if the current segment intersects with another segment S.
  //!   - If they intersect, returns `true` and updates parameters `s` and `t` with the intersection information.
  //!
  //! **Usage**
  //!
  //! \code
  //! #include "Point2D.h" // Assuming Point2D is defined elsewhere
  //!
  //! // Create two points for the segment
  //! Point2D<double> A, B;
  //! A.to_eigen() << 1.0, 2.0; // Set coordinates for point A
  //! B.to_eigen() << 3.0, 4.0; // Set coordinates for point B
  //!
  //! // Create a Segment2D object
  //! Segment2D<double> segment(A, B);
  //!
  //! // Access endpoints
  //! Point2D<double> start = segment.Pa(); // Get endpoint A
  //! Point2D<double> end = segment.Pb(); // Get endpoint B
  //!
  //! // Project a point onto the segment
  //! Point2D<double> P;
  //! P.to_eigen() << 2.0, 3.0;
  //! Real param;
  //! segment.projection(P, param); // Project P onto the segment
  //!
  //! // Evaluate a point on the segment at parameter s
  //! Real s = 0.5;
  //! Point2D<double> eval_point = segment.eval(s); // Evaluate the midpoint of the segment
  //! \endcode
  //!

  template <typename Real>
  class Segment2D {
    Point2D<Real> m_Pa; //!< First endpoint of the segment.
    Point2D<Real> m_Pb; //!< Second endpoint of the segment.
  public:

    //!
    //! \brief Default constructor for the segment.
    //!
    //! Initializes an empty segment with uninitialized endpoints.
    //!
    Segment2D() = default;
    //~Segment2D() = default;

    //!
    //! \brief Copy constructor for the segment.
    //!
    //! \param S The segment to copy from.
    //!
    Segment2D( Segment2D<Real> const & S )
    : m_Pa(S.m_Pa)
    , m_Pb(S.m_Pb)
    {}

    //!
    //! \brief Constructor to create a segment from two points.
    //!
    //! \param A The first endpoint of the segment.
    //! \param B The second endpoint of the segment.
    //!
    Segment2D( Point2D<Real> const & A, Point2D<Real> const & B )
    : m_Pa(A)
    , m_Pb(B)
    {}

    //!
    //! \brief Setup the endpoints of the segment using Point2D objects.
    //!
    //! \param A The first endpoint of the segment.
    //! \param B The second endpoint of the segment.
    //!
    void
    setup( Point2D<Real> const & A, Point2D<Real> const & B ) {
      m_Pa = A;
      m_Pb = B;
    }

    //!
    //! \brief Setup the endpoints of the segment using arrays.
    //!
    //! \param A Array containing the coordinates of the first endpoint.
    //! \param B Array containing the coordinates of the second endpoint.
    //!
    void
    setup( Real const A[], Real const B[] ) {
      m_Pa.coeffRef(0) = A[0]; m_Pa.coeffRef(1) = A[1];
      m_Pb.coeffRef(0) = B[0]; m_Pb.coeffRef(1) = B[1];
    }

    //Real          signed_distance( Point2D<Real> const & P ) const;

    //!
    //! \brief Projects a point onto the segment.
    //!
    //! \param P The point to be projected onto the segment.
    //! \param s The parameter along the segment where the projection occurs.
    //!
    //! \return `true` if the projection is successful, `false` otherwise.
    //!
    bool projection( Point2D<Real> const & P, Real & s ) const;

    //!
    //! \brief Projects a point onto the segment and returns the projected point.
    //!
    //! \param P The point to be projected onto the segment.
    //! \param s The parameter along the segment where the projection occurs.
    //! \param t The parameter along the other segment (if applicable).
    //!
    //! \return The projected point on the segment.
    //!
    Point2D<Real> projection( Point2D<Real> const & P, Real & s, Real & t ) const;

    //!
    //! \brief Evaluates the segment at a given parameter `s`.
    //!
    //! \param s The parameter to evaluate on the segment.
    //!
    //! \return The point on the segment corresponding to the parameter `s`.
    //!
    Point2D<Real> eval( Real & s ) const;

    //!
    //! \brief Evaluates the segment based on two parameters `s` and `t`.
    //!
    //! \param s The first parameter to evaluate.
    //! \param t The second parameter to evaluate.
    //!
    //! \return The point on the segment corresponding to the parameters `s` and `t`.
    //!
    Point2D<Real> eval( Real & s, Real & t ) const;

    //!
    //! \brief Returns the first endpoint of the segment.
    //!
    //! \return The first endpoint \f$ P_a \f$.
    //!
    Point2D<Real> const & Pa() const { return m_Pa; }


    //!
    //! \brief Returns the second endpoint of the segment.
    //!
    //! \return The second endpoint \f$ P_b \f$.
    //!
    Point2D<Real> const & Pb() const { return m_Pb; }

    //!
    //! \brief Computes the bounding box of the segment.
    //!
    //! \param pmin The minimum point of the bounding box.
    //! \param pmax The maximum point of the bounding box.
    //!
    void bbox( Point2D<Real> & pmin, Point2D<Real> & pmax ) const;

    //!
    //! \brief Checks for intersection with another segment.
    //!
    //! \param S The other segment to check for intersection.
    //! \param s The parameter along the current segment at the intersection point.
    //! \param t The parameter along the other segment at the intersection point.
    //!
    //! \return `true` if the segments intersect, `false` otherwise.
    //!
    bool intersect( Segment2D<Real> const & S, Real & s, Real & t ) const;

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
  //! \f$ P_{\text{min}} \f$ and \f$ P_{\text{max}} \f$, which represent the minimum
  //! and maximum corners of the box respectively. This class provides a structure
  //! to represent and manipulate a bounding box in a two-dimensional space.
  //!
  //! \tparam Real The data type used for the coordinates of the box's corners,
  //! typically a floating-point type (e.g., float, double).
  //!
  template <typename Real>
  class Box2D {
      Point2D<Real> m_Pmin; //!< The minimum corner of the box (bottom-left corner).
      Point2D<Real> m_Pmax; //!< The maximum corner of the box (top-right corner).

  public:
      //!
      //! \brief Default constructor for the bounding box.
      //!
      //! Initializes a new bounding box. The minimum and maximum points are
      //! uninitialized. This constructor can be used to create an empty box
      //! which can later be defined using specific points.
      Box2D() {}
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
  //! \f$ P_a \f$, \f$ P_b \f$, and \f$ P_c \f$, which represent the vertices of the triangle.
  //! This class provides a structure to represent and manipulate a triangle in a
  //! two-dimensional space, allowing for operations such as area calculation,
  //! centroid determination, and intersection checks with other geometric entities.
  //!
  //! \tparam Real The data type used for the coordinates of the triangle's vertices,
  //! typically a floating-point type (e.g., float, double).
  //!
  template <typename Real>
  class Triangle2D {
    Point2D<Real> m_Pa; //!< The first vertex of the triangle.
    Point2D<Real> m_Pb; //!< The second vertex of the triangle.
    Point2D<Real> m_Pc; //!< The third vertex of the triangle.

  public:
    //!
    //! \brief Default constructor for the triangle.
    //!
    //! Initializes a new triangle with uninitialized vertices.
    //! This constructor can be used to create an empty triangle
    //! which can later be defined using specific points.
    Triangle2D() {}
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
  //! \brief A class representing a 2D polygon using a dynamic matrix of vertices.
  //!
  //! The `Polygon2D` class defines a 2D polygon by storing its vertices in a dynamic
  //! matrix. Each column of the matrix represents a vertex of the polygon in 2D space,
  //! allowing for the representation of polygons with varying numbers of vertices.
  //! This class inherits from `Eigen::Matrix` to leverage efficient matrix operations
  //! for polygon manipulations such as transformations and evaluations.
  //!
  //! \tparam Real The data type used for the coordinates of the polygon's vertices,
  //! typically a floating-point type (e.g., float, double).
  //!
  template <typename Real>
  class Polygon2D : public Eigen::Matrix<Real, 2, Eigen::Dynamic> {
  public:
    //!
    //! \brief Default constructor for the polygon.
    //!
    //! Initializes a new polygon with no vertices. The internal matrix
    //! structure is empty, and vertices can be added later through
    //! appropriate methods.
    Polygon2D() {}
  };

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */

  #ifndef UTILS_OS_WINDOWS
  extern template class Point2D<float>;
  extern template class Point2D<double>;

  extern template class Segment2D<float>;
  extern template class Segment2D<double>;

  extern template class Box2D<float>;
  extern template class Box2D<double>;

  extern template class Triangle2D<float>;
  extern template class Triangle2D<double>;

  extern template class Polygon2D<float>;
  extern template class Polygon2D<double>;
  #endif

  /*! @} */

}

#endif

//
// eof: Utils_GG2D.hh
//
