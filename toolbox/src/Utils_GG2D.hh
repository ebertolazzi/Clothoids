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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: Utils_GG2D.hh
///
#pragma once

#ifndef UTILS_GG2D_dot_HH
#define UTILS_GG2D_dot_HH

#include "Utils_eigen.hh"

namespace Utils {

  /*\
   |   ____       _       _   ____  ____
   |  |  _ \ ___ (_)_ __ | |_|___ \|  _ \
   |  | |_) / _ \| | '_ \| __| __) | | | |
   |  |  __/ (_) | | | | | |_ / __/| |_| |
   |  |_|   \___/|_|_| |_|\__|_____|____/
  \*/

  template <typename Real>
  class Point2D : public Eigen::Matrix<Real,2,1> {
    typedef Eigen::Matrix<Real,2,1> P2D;
  public:
    Point2D() = default;
    //~Point2D() = default;

    //Point2D<Real> &
    //operator = ( Point2D<Real> const & RHS ) {
    //  this->to_eigen() = RHS.to_eigen();
    //  return *this;
    //}

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

  template <typename Real>
  class Segment2D {
    Point2D<Real> m_Pa;
    Point2D<Real> m_Pb;
  public:

    Segment2D() = default;
    //~Segment2D() = default;

    Segment2D( Segment2D<Real> const & S )
    : m_Pa(S.m_Pa)
    , m_Pb(S.m_Pb)
    {}

    Segment2D( Point2D<Real> const & A, Point2D<Real> const & B )
    : m_Pa(A)
    , m_Pb(B)
    {}

    void
    setup( Point2D<Real> const & A, Point2D<Real> const & B ) {
      m_Pa = A;
      m_Pb = B;
    }

    void
    setup( Real const * A, Real const * B ) {
      m_Pa.coeffRef(0) = A[0]; m_Pa.coeffRef(1) = A[1];
      m_Pb.coeffRef(0) = B[0]; m_Pb.coeffRef(1) = B[1];
    }

    //Segment2D<Real> const &
    //operator = ( Segment2D<Real> const & RHS ) {
    //  m_Pa = RHS.m_Pa;
    //  m_Pb = RHS.m_Pb;
    //  return *this;
    //}

    //Real          signed_distance( Point2D<Real> const & P ) const;
    bool          projection( Point2D<Real> const & P, Real & s ) const;
    Point2D<Real> projection( Point2D<Real> const & P, Real & s, Real & t ) const;
    Point2D<Real> eval( Real & s ) const;
    Point2D<Real> eval( Real & s, Real & t ) const;

    Point2D<Real> const & Pa() const { return m_Pa; }
    Point2D<Real> const & Pb() const { return m_Pb; }

    void bbox( Point2D<Real> & pmin, Point2D<Real> & pmax ) const;
    bool intersect( Segment2D<Real> const & S, Real & s, Real & t ) const;

  };

  /*\
   |   ____            ____  ____
   |  | __ )  _____  _|___ \|  _ \
   |  |  _ \ / _ \ \/ / __) | | | |
   |  | |_) | (_) >  < / __/| |_| |
   |  |____/ \___/_/\_\_____|____/
  \*/

  template <typename Real>
  class Box2D {
    Point2D<Real> m_Pmin;
    Point2D<Real> m_Pmax;
  public:
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

  template <typename Real>
  class Triangle2D {
    Point2D<Real> m_Pa, m_Pb, m_Pc;
  public:
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

  template <typename Real>
  class Polygon2D : public Eigen::Matrix<Real,2,Eigen::Dynamic> {
  public:
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

}


#endif

///
/// eof: Utils_GG2D.hh
///
