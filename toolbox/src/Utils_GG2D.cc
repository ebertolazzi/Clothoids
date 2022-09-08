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
/// file: Utils_GG2D.cc
///

#include "Utils_GG2D.hh"
#include <cmath>

namespace Utils {

  using std::abs;

  /*\
   |   ____       _       _   ____  ____
   |  |  _ \ ___ (_)_ __ | |_|___ \|  _ \
   |  | |_) / _ \| | '_ \| __| __) | | | |
   |  |  __/ (_) | | | | | |_ / __/| |_| |
   |  |_|   \___/|_|_| |_|\__|_____|____/
  \*/

  /*\
   |   ____                                  _   ____  ____
   |  / ___|  ___  __ _ _ __ ___   ___ _ __ | |_|___ \|  _ \
   |  \___ \ / _ \/ _` | '_ ` _ \ / _ \ '_ \| __| __) | | | |
   |   ___) |  __/ (_| | | | | | |  __/ | | | |_ / __/| |_| |
   |  |____/ \___|\__, |_| |_| |_|\___|_| |_|\__|_____|____/
   |              |___/
  \*/

  template <typename Real>
  bool
  Segment2D<Real>::projection( Point2D<Real> const & P, Real & s ) const {
    //
    // || Pa + t*(Pb-Pa) - P ||^2 --> min
    // ( Pa + t*(Pb-Pa) - P ) . (Pb - Pa) = 0
    // s = ( P - Pa) . ( Pb - Pa) / ( Pb - Pa )^2
    //
    Point2D<Real> PA, BA, RES;
    PA.noalias() = P-m_Pa;
    BA.noalias() = m_Pb-m_Pa;
    s = PA.dot(BA) / BA.dot(BA);
    return s >= 0 && s <= 1;
  }

  template <typename Real>
  Point2D<Real>
  Segment2D<Real>::projection( Point2D<Real> const & P, Real & s, Real & t ) const {
    //
    // || Pa + t*(Pb-Pa) - P ||^2 --> min
    // ( Pa + t*(Pb-Pa) - P ) . (Pb - Pa) = 0
    // s = ( P - Pa) . ( Pb - Pa) / ( Pb - Pa )^2
    //
    Point2D<Real> PA, BA, PP, N;
    PA.noalias() = P-m_Pa;
    BA.noalias() = m_Pb-m_Pa;
    s = PA.dot(BA) / BA.dot(BA);
    if      ( s < 0 ) s = 0;
    else if ( s > 1 ) s = 1;
    PP.noalias() = m_Pa + s * BA;
    //
    // P   = Pa + s*(Pb-Pa) + t*N
    // P.N = (Pa + s*(Pb-Pa)).N + t N.N
    //
    N.coeffRef(0) = -BA.coeff(1);
    N.coeffRef(1) = BA.coeff(0);
    t = (P-PP).dot(N)/N.norm();
    return PP;
  }

  template <typename Real>
  Point2D<Real>
  Segment2D<Real>::eval( Real & s ) const {
    Point2D<Real> RES, BA;
    BA.noalias() = m_Pb-m_Pa;
    RES.noalias() = m_Pa + s*BA;
    return RES;
  }

  template <typename Real>
  Point2D<Real>
  Segment2D<Real>::eval( Real & s, Real & t ) const {
    Point2D<Real> BA, N, RES;
    BA.noalias() = m_Pb-m_Pa;
    N.coeffRef(0) = -BA.coeff(1);
    N.coeffRef(1) = BA.coeff(0);
    N.normalize();
    RES.noalias() = m_Pa + s*BA + t * N;
    return RES;
  }

  template <typename Real>
  void
  Segment2D<Real>::bbox( Point2D<Real> & pmin, Point2D<Real> & pmax ) const {
    pmin.noalias() = m_Pa.cwiseMin(m_Pb);
    pmax.noalias() = m_Pa.cwiseMax(m_Pb);
  }

  template <typename Real>
  bool
  Segment2D<Real>::intersect( Segment2D<Real> const & S, Real & s, Real & t ) const {
    s = t = 0;
    // check if collinear
    Point2D<Real> D1, D2, RES;
    D1.noalias() = m_Pb - m_Pa;
    D2.noalias() = S.m_Pb - S.m_Pa;
    Real CX = (D1.x()*D2.y()-D1.y()*D2.x())/(D1.norm()*D2.norm());
    if ( abs(CX) <= 10*machine_eps<Real>() ) {
      // collinear point check if overlap
      if ( this->projection( S.m_Pa, s ) ) return true;
      if ( this->projection( S.m_Pb, s ) ) return true;
      s = t = 0;
      return false;
    }
    // regular case, find intersection
    // Pa + s*D1 = A + t*D2
    // [D1,-D2]*[s;t] = A-Pa;
    Eigen::Matrix<Real,2,2> M;
    M.col(0) = D1;
    M.col(1) = -D2;
    //Eigen::ColPivHouseholderQR<Eigen::Matrix<Real,2,2>> QR(M);
    RES.noalias() = M.colPivHouseholderQr().solve(S.m_Pa-m_Pa);
    s = RES.coeff(0);
    t = RES.coeff(1);
    return s >= 0 && s <= 1 && t >= 0 && t <= 1;
  }

  /*\
   |   ____            ____  ____
   |  | __ )  _____  _|___ \|  _ \
   |  |  _ \ / _ \ \/ / __) | | | |
   |  | |_) | (_) >  < / __/| |_| |
   |  |____/ \___/_/\_\_____|____/
  \*/

  /*\
   |   _____     _                   _      ____  ____
   |  |_   _| __(_) __ _ _ __   __ _| | ___|___ \|  _ \
   |    | || '__| |/ _` | '_ \ / _` | |/ _ \ __) | | | |
   |    | || |  | | (_| | | | | (_| | |  __// __/| |_| |
   |    |_||_|  |_|\__,_|_| |_|\__, |_|\___|_____|____/
   |                           |___/
  \*/

  /*\
   |   ____       _                         ____  ____
   |  |  _ \ ___ | |_   _  __ _  ___  _ __ |___ \|  _ \
   |  | |_) / _ \| | | | |/ _` |/ _ \| '_ \  __) | | | |
   |  |  __/ (_) | | |_| | (_| | (_) | | | |/ __/| |_| |
   |  |_|   \___/|_|\__, |\__, |\___/|_| |_|_____|____/
   |                |___/ |___/
  \*/

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */

  template class Point2D<float>;
  template class Point2D<double>;

  template class Segment2D<float>;
  template class Segment2D<double>;

  template class Box2D<float>;
  template class Box2D<double>;

  template class Triangle2D<float>;
  template class Triangle2D<double>;

  template class Polygon2D<float>;
  template class Polygon2D<double>;

}

///
/// eof: Utils_GG2D.cc
///
