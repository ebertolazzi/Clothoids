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

#include "Clothoid.hh"
#include "Biarc.hh"

#include <cmath>
#include <cfloat>
#include <fstream>

#ifdef __GCC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

namespace G2lib {

  using namespace std ;

  /*\
   |   ____ _       _   _           _     _ _     _     _
   |  / ___| | ___ | |_| |__   ___ (_) __| | |   (_)___| |_
   | | |   | |/ _ \| __| '_ \ / _ \| |/ _` | |   | / __| __|
   | | |___| | (_) | |_| | | | (_) | | (_| | |___| \__ \ |_
   |  \____|_|\___/ \__|_| |_|\___/|_|\__,_|_____|_|___/\__|
   |
  \*/

  ClothoidList::~ClothoidList() {
    s0.clear() ;
    clotoidList.clear();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::copy( ClothoidList const & L ) {
    s0.resize( L.s0.size() ) ;
    std::copy( L.s0.begin(), L.s0.end(), s0.begin() ) ;
    clotoidList.resize( L.clotoidList.size() ) ;
    std::copy( L.clotoidList.begin(), L.clotoidList.end(), clotoidList.begin() ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::reserve( indexType n ) {
    s0.reserve(size_t(n+1)) ;
    clotoidList.reserve(size_t(n)) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back( ClothoidCurve const & c ) {
    if ( clotoidList.empty() ) {
      s0.push_back(0) ;
      s0.push_back(c.length()) ;
    } else {
      s0.push_back(s0.back()+c.length()) ;
    }
    clotoidList.push_back(c) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back(
    valueType kappa0, valueType dkappa, valueType L
  ) {
    G2LIB_ASSERT( !clotoidList.empty(),
                  "ClothoidList::push_back_G1(...) empty list!");
    ClothoidCurve c ;
    valueType x0     = clotoidList.back().xEnd() ;
    valueType y0     = clotoidList.back().yEnd() ;
    valueType theta0 = clotoidList.back().thetaEnd() ;
    c.build( x0, y0, theta0, kappa0, dkappa, L ) ;
    push_back( c ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back(
    valueType x0, valueType y0, valueType theta0,
    valueType kappa0, valueType dkappa, valueType L
  ) {
    ClothoidCurve c ;
    c.build( x0, y0, theta0, kappa0, dkappa, L ) ;
    push_back( c ) ;
  }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back_G1(
    valueType x1, valueType y1, valueType theta1
  ) {
    G2LIB_ASSERT( !clotoidList.empty(),
                  "ClothoidList::push_back_G1(...) empty list!");
    ClothoidCurve c ;
    valueType x0     = clotoidList.back().xEnd() ;
    valueType y0     = clotoidList.back().yEnd() ;
    valueType theta0 = clotoidList.back().thetaEnd() ;
    c.build_G1( x0, y0, theta0, x1, y1, theta1 ) ;
    push_back( c ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back_G1(
    valueType x0, valueType y0, valueType theta0,
    valueType x1, valueType y1, valueType theta1
  ) {
    ClothoidCurve c ;
    c.build_G1( x0, y0, theta0, x1, y1, theta1 ) ;
    push_back( c ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidList::build_G1( indexType       n,
                          valueType const x[],
                          valueType const y[] ) {
    init() ;
    reserve( n-1 );
    ClothoidCurve c ;

    G2LIB_ASSERT( n > 1, "ClothoidList::build_G1, at least 2 points are necessary" ) ;

    if ( n == 2 ) {

      valueType theta = atan2( y[1] - y[0], x[1] - x[0] ) ;
      c.build_G1( x[0], y[0], theta, x[1], y[1], theta ) ;
      push_back(c);

    } else {

      Biarc b ;
      bool ok, ciclic = hypot( x[0]-x[n-1], y[0]-y[n-1] ) < 1e-10 ;
      valueType thetaC(0) ;
      if ( ciclic ) {
        ok = b.build_3P( x[n-2], y[n-2], x[0], y[0], x[1], y[1] ) ;
        G2LIB_ASSERT( ok, "ClothoidList::build_G1, failed" ) ;
        thetaC = b.thetaStar();
      }
      ok = b.build_3P( x[0], y[0], x[1], y[1], x[2], y[2] ) ;
      G2LIB_ASSERT( ok, "ClothoidList::build_G1, failed" ) ;
      valueType theta0 = ciclic ? thetaC : b.thetaBegin0();
      valueType theta1 = b.thetaStar();
      c.build_G1( x[0], y[0], theta0, x[1], y[1], theta1 ) ;
      push_back(c);
      for ( indexType k = 2 ; k < n-1 ; ++k ) {
        theta0 = theta1 ;
        ok = b.build_3P( x[k-1], y[k-1], x[k], y[k], x[k+1], y[k+1] ) ;
        G2LIB_ASSERT( ok, "ClothoidList::build_G1, failed" ) ;
        theta1 = b.thetaStar();
        c.build_G1( x[k-1], y[k-1], theta0, x[k], y[k], theta1 ) ;
        push_back(c);
      }
      theta0 = theta1 ;
      theta1 = ciclic ? thetaC : b.thetaEnd1();
      c.build_G1( x[n-2], y[n-2], theta0, x[n-1], y[n-1], theta1 ) ;
      push_back(c);

    }
    return true ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidList::build_G1( indexType       n,
                          valueType const x[],
                          valueType const y[],
                          valueType const theta[] ) {

    G2LIB_ASSERT( n > 1, "ClothoidList::build_G1, at least 2 points are necessary" ) ;

    init() ;
    reserve( n-1 );
    ClothoidCurve c ;
    for ( indexType k = 1 ; k < n ; ++k ) {
      c.build_G1( x[k-1], y[k-1], theta[k-1], x[k], y[k], theta[k] ) ;
      push_back(c);
    }
    return true ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidList::build_theta( indexType       n,
                             valueType const x[],
                             valueType const y[],
                             valueType       theta[] ) const {
    G2LIB_ASSERT( n > 1, "ClothoidList::build_theta, at least 2 points are necessary" ) ;

    if ( n == 2 ) {
      theta[0] = theta[1] = atan2( y[1] - y[0], x[1] - x[0] ) ;
    } else {
      Biarc b ;
      bool ok, ciclic = hypot( x[0]-x[n-1], y[0]-y[n-1] ) < 1e-10 ;
      if ( ciclic ) {
        ok = b.build_3P( x[n-2], y[n-2], x[0], y[0], x[1], y[1] ) ;
        G2LIB_ASSERT( ok, "ClothoidList::build_theta, failed" ) ;
        theta[0] = theta[n-1] = b.thetaStar();
      }
      for ( indexType k = 1 ; k < n-1 ; ++k ) {
        ok = b.build_3P( x[k-1], y[k-1], x[k], y[k], x[k+1], y[k+1] ) ;
        G2LIB_ASSERT( ok, "ClothoidList::build_theta, failed" ) ;
        theta[k] = b.thetaStar();
        if ( k == 1   && !ciclic ) theta[0]   = b.thetaBegin0();
        if ( k == n-2 && !ciclic ) theta[n-1] = b.thetaEnd1();
      }
    }
    return true ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ClothoidCurve const &
  ClothoidList::get( indexType idx ) const {
    G2LIB_ASSERT( !clotoidList.empty(), "ClothoidList::get( " << idx << " ) empty list" );
    G2LIB_ASSERT( idx >= 0 && idx < indexType(clotoidList.size()),
                  "ClothoidList::get( " << idx <<
                  " ) bad index, must be in [0," <<
                  clotoidList.size()-1 << "]" );
    return clotoidList[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ClothoidCurve const &
  ClothoidList::getAtS( valueType s ) const {
    findAtS(s);
    return get(last_idx) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidList::findAtS( valueType s ) const {
    indexType ns = indexType(clotoidList.size()) ;
    G2LIB_ASSERT( last_idx >= 0 && last_idx < ns,
                  "ClothoidList::findAtS( " << s << ", " << last_idx <<
                  " ) bad index" );
    valueType const * sL = &s0[last_idx] ;
    if ( s < sL[0] ) {
      if ( s > s0.front() ) {
        valueType const * sB = &s0.front() ;
        last_idx = indexType(std::lower_bound( sB, sL, s )-sB) ;
      } else {
        last_idx = 0 ;
      }
    } else if ( s > sL[1] ) {
      if ( s < s0.back() ) {
        valueType const * sE = &s0[ns+1] ; // past to the last
        last_idx += indexType(std::lower_bound( sL, sE, s )-sL) ;
      } else {
        last_idx = ns-1 ;
      }
    } else {
      return true ; // vale intervallo precedente
    }
    if ( s0[last_idx] > s ) --last_idx ; // aggiustamento caso di bordo
    G2LIB_ASSERT( last_idx >= 0 && last_idx < ns,
                  "ClothoidList::findAtS( " << s <<
                  ") last_idx = " << last_idx <<
                  " range [" << s0.front() << ", " << s0.back() << "]" ) ;
    return true ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  ClothoidList::theta( valueType s ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.theta( s - s0[last_idx] ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  ClothoidList::theta_D( valueType s ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.theta_D( s - s0[last_idx] ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  ClothoidList::theta_DD( valueType s ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.theta_DD( s - s0[last_idx] ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  ClothoidList::X( valueType s ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.X( s - s0[last_idx] ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  ClothoidList::Y( valueType s ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.Y( s - s0[last_idx] ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval( valueType   s,
                      valueType & theta,
                      valueType & kappa,
                      valueType & x,
                      valueType & y ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.eval( s - s0[last_idx], theta, kappa, x, y ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval( valueType   s,
                      valueType & x,
                      valueType & y ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.eval( s - s0[last_idx], x, y ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_D( valueType   s,
                        valueType & x_D,
                        valueType & y_D ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.eval_D( s - s0[last_idx], x_D, y_D ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_DD( valueType   s,
                         valueType & x_DD,
                         valueType & y_DD ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.eval_DD( s - s0[last_idx], x_DD, y_DD ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_DDD( valueType   s,
                          valueType & x_DDD,
                          valueType & y_DDD ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.eval_DDD( s - s0[last_idx], x_DDD, y_DDD ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // offset curve
  void
  ClothoidList::eval( valueType   s,
                      valueType   offs,
                      valueType & x,
                      valueType & y ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.eval( s - s0[last_idx], offs, x, y ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_D( valueType   s,
                        valueType   offs,
                        valueType & x_D,
                        valueType & y_D ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.eval_D( s - s0[last_idx], offs, x_D, y_D ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_DD( valueType   s,
                         valueType   offs,
                         valueType & x_DD,
                         valueType & y_DD ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.eval_DD( s - s0[last_idx], offs, x_DD, y_DD ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_DDD( valueType   s,
                          valueType   offs,
                          valueType & x_DDD,
                          valueType & y_DDD ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.eval_DDD( s - s0[last_idx], offs, x_DDD, y_DDD ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::getSTK( valueType s[],
                        valueType theta[],
                        valueType kappa[] ) const {
    std::vector<ClothoidCurve>::const_iterator ic = clotoidList.begin() ;
    indexType k  = 0 ;
    valueType ss = 0 ;
    while ( ic != clotoidList.end() ) {
      s[k]     = ss ;
      theta[k] = ic->thetaBegin();
      kappa[k] = ic->kappaBegin();
      ss       += ic->length() ;
      ++k ;
      ++ic ;
    }
    --ic ;
    s[k]     = ss ;
    theta[k] = ic->thetaEnd();
    kappa[k] = ic->kappaEnd();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::getXY( valueType x[], valueType y[] ) const {
    std::vector<ClothoidCurve>::const_iterator ic = clotoidList.begin() ;
    indexType k  = 0 ;
    while ( ic != clotoidList.end() ) {
      x[k] = ic->xBegin();
      y[k] = ic->yBegin();
      ++k ; ++ic ;
    }
    --ic ;
    x[k] = ic->xEnd();
    y[k] = ic->yEnd();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::getDeltaTheta( valueType deltaTheta[] ) const {
    std::vector<ClothoidCurve>::const_iterator ic = clotoidList.begin() ;
    indexType k = 0 ;
    for ( ++ic ; ic != clotoidList.end() ; ++ic, ++k ) {
      valueType tmp = ic->thetaBegin()-ic[-1].thetaEnd();
      if      ( tmp >  m_pi ) tmp -= m_2pi ;
      else if ( tmp < -m_pi ) tmp += m_2pi ;
      deltaTheta[k] = tmp;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::getDeltaKappa( valueType deltaKappa[] ) const {
    std::vector<ClothoidCurve>::const_iterator ic = clotoidList.begin() ;
    indexType k = 0 ;
    for ( ++ic ; ic != clotoidList.end() ; ++ic, ++k  )
      deltaKappa[k] = ic->kappaBegin()-ic[-1].kappaEnd();

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  valueType
  ClothoidList::closestPoint( valueType   x,
                              valueType   y,
                              valueType & X,
                              valueType & Y,
                              valueType & S ) const {
    G2LIB_ASSERT( !clotoidList.empty(), "ClothoidList::closestPoint, empty list" );
    std::vector<ClothoidCurve>::const_iterator ic = clotoidList.begin() ;
    std::vector<valueType>::const_iterator     is = s0.begin() ;
    valueType DST = ic->closestPoint( x, y, X, Y, S );
    for ( ++ic, ++is ; ic != clotoidList.end() ; ++ic, ++is ) {
      valueType X1, Y1, S1 ;
      valueType DST1 = ic->closestPoint( x, y, X1, Y1, S1 );
      if ( DST1 < DST ) {
        DST = DST1 ;
        X   = X1;
        Y   = Y1;
        S   = *is + S1;
      }
    }
    return DST ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::rotate( valueType angle, valueType cx, valueType cy ) {
    std::vector<ClothoidCurve>::iterator ic = clotoidList.begin() ;
    for ( ; ic != clotoidList.end() ; ++ic ) ic->rotate( angle, cx, cy ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::translate( valueType tx, valueType ty ) {
    std::vector<ClothoidCurve>::iterator ic = clotoidList.begin() ;
    for ( ; ic != clotoidList.end() ; ++ic ) ic->translate( tx, ty ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::changeOrigin( valueType newx0, valueType newy0 ) {
    std::vector<ClothoidCurve>::iterator ic = clotoidList.begin() ;
    for ( ; ic != clotoidList.end() ; ++ic ) {
      ic->changeOrigin( newx0, newy0 ) ;
      newx0 = ic->xEnd() ;
      newy0 = ic->yEnd() ;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::scale( valueType sfactor ) {
    std::vector<ClothoidCurve>::iterator ic = clotoidList.begin() ;
    valueType newx0 = ic->xBegin() ;
    valueType newy0 = ic->yBegin() ;
    for ( ; ic != clotoidList.end() ; ++ic ) {
      ic->scale( sfactor ) ;
      ic->changeOrigin( newx0, newy0 ) ;
      newx0 = ic->xEnd() ;
      newy0 = ic->yEnd() ;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::reverse() {
    std::vector<ClothoidCurve>::iterator ic = clotoidList.begin() ;
    for ( ; ic != clotoidList.end() ; ++ic ) ic->reverse() ;
    std::reverse( clotoidList.begin(), clotoidList.end() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::export_table( std::ostream & stream ) const {
    stream << "x\ty\ttheta0\tkappa0\tdkappa\tL\n" ;
    std::vector<ClothoidCurve>::const_iterator ic = clotoidList.begin() ;
    for ( ; ic != clotoidList.end() ; ++ic )
      stream << ic->xBegin()     << '\t'
             << ic->yBegin()     << '\t'
             << ic->thetaBegin() << '\t'
             << ic->kappaBegin() << '\t'
             << ic->kappa_D()    << '\t'
             << ic->length()     << '\n' ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::export_ruby( std::ostream & stream ) const {
    stream << "data = {\n" ;
    std::vector<ClothoidCurve>::const_iterator ic = clotoidList.begin() ;
    for ( ; ic != clotoidList.end() ; ++ic )
      stream << ic->xBegin()     << '\t'
             << ic->yBegin()     << '\t'
             << ic->thetaBegin() << '\t'
             << ic->kappaBegin() << '\t'
             << ic->kappa_D()    << '\t'
             << ic->length()     << '\n' ;
    stream << "}\n" ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}

// EOF: ClothoidList.cc
