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

  using namespace std;

  /*\
   |   ____ _       _   _           _     _ _     _     _
   |  / ___| | ___ | |_| |__   ___ (_) __| | |   (_)___| |_
   | | |   | |/ _ \| __| '_ \ / _ \| |/ _` | |   | / __| __|
   | | |___| | (_) | |_| | | | (_) | | (_| | |___| \__ \ |_
   |  \____|_|\___/ \__|_| |_|\___/|_|\__,_|_____|_|___/\__|
   |
  \*/

  ClothoidList::~ClothoidList() {
    s0.clear();
    clotoidList.clear();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::copy( ClothoidList const & L ) {
    s0.resize( L.s0.size() );
    std::copy( L.s0.begin(), L.s0.end(), s0.begin() );
    clotoidList.resize( L.clotoidList.size() );
    std::copy( L.clotoidList.begin(), L.clotoidList.end(), clotoidList.begin() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::reserve( int_type n ) {
    s0.reserve(size_t(n+1));
    clotoidList.reserve(size_t(n));
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back( ClothoidCurve const & c ) {
    if ( clotoidList.empty() ) {
      s0.push_back(0);
      s0.push_back(c.length());
    } else {
      s0.push_back(s0.back()+c.length());
    }
    clotoidList.push_back(c);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back(
    real_type kappa0, real_type dkappa, real_type L
  ) {
    G2LIB_ASSERT( !clotoidList.empty(),
                  "ClothoidList::push_back_G1(...) empty list!");
    ClothoidCurve c;
    real_type x0     = clotoidList.back().xEnd();
    real_type y0     = clotoidList.back().yEnd();
    real_type theta0 = clotoidList.back().thetaEnd();
    c.build( x0, y0, theta0, kappa0, dkappa, L );
    push_back( c );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back(
    real_type x0, real_type y0, real_type theta0,
    real_type kappa0, real_type dkappa, real_type L
  ) {
    ClothoidCurve c;
    c.build( x0, y0, theta0, kappa0, dkappa, L );
    push_back( c );
  }
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back_G1(
    real_type x1, real_type y1, real_type theta1
  ) {
    G2LIB_ASSERT( !clotoidList.empty(),
                  "ClothoidList::push_back_G1(...) empty list!");
    ClothoidCurve c;
    real_type x0     = clotoidList.back().xEnd();
    real_type y0     = clotoidList.back().yEnd();
    real_type theta0 = clotoidList.back().thetaEnd();
    c.build_G1( x0, y0, theta0, x1, y1, theta1 );
    push_back( c );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::push_back_G1(
    real_type x0, real_type y0, real_type theta0,
    real_type x1, real_type y1, real_type theta1
  ) {
    ClothoidCurve c;
    c.build_G1( x0, y0, theta0, x1, y1, theta1 );
    push_back( c );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidList::build_G1( int_type        n,
                          real_type const x[],
                          real_type const y[] ) {
    init();
    reserve( n-1 );
    ClothoidCurve c;

    G2LIB_ASSERT( n > 1, "ClothoidList::build_G1, at least 2 points are necessary" );

    if ( n == 2 ) {

      real_type theta = atan2( y[1] - y[0], x[1] - x[0] );
      c.build_G1( x[0], y[0], theta, x[1], y[1], theta );
      push_back(c);

    } else {

      Biarc b;
      bool ok, ciclic = hypot( x[0]-x[n-1], y[0]-y[n-1] ) < 1e-10;
      real_type thetaC(0);
      if ( ciclic ) {
        ok = b.build_3P( x[n-2], y[n-2], x[0], y[0], x[1], y[1] );
        G2LIB_ASSERT( ok, "ClothoidList::build_G1, failed" );
        thetaC = b.thetaStar();
      }
      ok = b.build_3P( x[0], y[0], x[1], y[1], x[2], y[2] );
      G2LIB_ASSERT( ok, "ClothoidList::build_G1, failed" );
      real_type theta0 = ciclic ? thetaC : b.thetaBegin0();
      real_type theta1 = b.thetaStar();
      c.build_G1( x[0], y[0], theta0, x[1], y[1], theta1 );
      push_back(c);
      for ( int_type k = 2; k < n-1; ++k ) {
        theta0 = theta1;
        ok = b.build_3P( x[k-1], y[k-1], x[k], y[k], x[k+1], y[k+1] );
        G2LIB_ASSERT( ok, "ClothoidList::build_G1, failed" );
        theta1 = b.thetaStar();
        c.build_G1( x[k-1], y[k-1], theta0, x[k], y[k], theta1 );
        push_back(c);
      }
      theta0 = theta1;
      theta1 = ciclic ? thetaC : b.thetaEnd1();
      c.build_G1( x[n-2], y[n-2], theta0, x[n-1], y[n-1], theta1 );
      push_back(c);

    }
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidList::build_G1( int_type        n,
                          real_type const x[],
                          real_type const y[],
                          real_type const theta[] ) {

    G2LIB_ASSERT( n > 1, "ClothoidList::build_G1, at least 2 points are necessary" );

    init();
    reserve( n-1 );
    ClothoidCurve c;
    for ( int_type k = 1; k < n; ++k ) {
      c.build_G1( x[k-1], y[k-1], theta[k-1], x[k], y[k], theta[k] );
      push_back(c);
    }
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidList::build_theta( int_type        n,
                             real_type const x[],
                             real_type const y[],
                             real_type       theta[] ) const {
    G2LIB_ASSERT( n > 1, "ClothoidList::build_theta, at least 2 points are necessary" );

    if ( n == 2 ) {
      theta[0] = theta[1] = atan2( y[1] - y[0], x[1] - x[0] );
    } else {
      Biarc b;
      bool ok, ciclic = hypot( x[0]-x[n-1], y[0]-y[n-1] ) < 1e-10;
      if ( ciclic ) {
        ok = b.build_3P( x[n-2], y[n-2], x[0], y[0], x[1], y[1] );
        G2LIB_ASSERT( ok, "ClothoidList::build_theta, failed" );
        theta[0] = theta[n-1] = b.thetaStar();
      }
      for ( int_type k = 1; k < n-1; ++k ) {
        ok = b.build_3P( x[k-1], y[k-1], x[k], y[k], x[k+1], y[k+1] );
        G2LIB_ASSERT( ok, "ClothoidList::build_theta, failed" );
        theta[k] = b.thetaStar();
        if ( k == 1   && !ciclic ) theta[0]   = b.thetaBegin0();
        if ( k == n-2 && !ciclic ) theta[n-1] = b.thetaEnd1();
      }
    }
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ClothoidCurve const &
  ClothoidList::get( int_type idx ) const {
    G2LIB_ASSERT( !clotoidList.empty(), "ClothoidList::get( " << idx << " ) empty list" );
    G2LIB_ASSERT( idx >= 0 && idx < int_type(clotoidList.size()),
                  "ClothoidList::get( " << idx <<
                  " ) bad index, must be in [0," <<
                  clotoidList.size()-1 << "]" );
    return clotoidList[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ClothoidCurve const &
  ClothoidList::getAtS( real_type s ) const {
    findAtS(s);
    return get(last_idx);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  ClothoidList::findAtS( real_type s ) const {
    int_type ns = int_type(clotoidList.size());
    G2LIB_ASSERT( last_idx >= 0 && last_idx < ns,
                  "ClothoidList::findAtS( " << s << ", " << last_idx <<
                  " ) bad index" );
    real_type const * sL = &s0[last_idx];
    if ( s < sL[0] ) {
      if ( s > s0.front() ) {
        real_type const * sB = &s0.front();
        last_idx = int_type(std::lower_bound( sB, sL, s )-sB);
      } else {
        last_idx = 0;
      }
    } else if ( s > sL[1] ) {
      if ( s < s0.back() ) {
        real_type const * sE = &s0[ns+1]; // past to the last
        last_idx += int_type(std::lower_bound( sL, sE, s )-sL);
      } else {
        last_idx = ns-1;
      }
    } else {
      return true; // vale intervallo precedente
    }
    if ( s0[last_idx] > s ) --last_idx; // aggiustamento caso di bordo
    G2LIB_ASSERT( last_idx >= 0 && last_idx < ns,
                  "ClothoidList::findAtS( " << s <<
                  ") last_idx = " << last_idx <<
                  " range [" << s0.front() << ", " << s0.back() << "]" );
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::theta( real_type s ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.theta( s - s0[last_idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::theta_D( real_type s ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.theta_D( s - s0[last_idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::theta_DD( real_type s ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.theta_DD( s - s0[last_idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::X( real_type s ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.X( s - s0[last_idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::Y( real_type s ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.Y( s - s0[last_idx] );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval( real_type   s,
                      real_type & theta,
                      real_type & kappa,
                      real_type & x,
                      real_type & y ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.eval( s - s0[last_idx], theta, kappa, x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval( real_type   s,
                      real_type & x,
                      real_type & y ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.eval( s - s0[last_idx], x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_D( real_type   s,
                        real_type & x_D,
                        real_type & y_D ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.eval_D( s - s0[last_idx], x_D, y_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_DD( real_type   s,
                         real_type & x_DD,
                         real_type & y_DD ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.eval_DD( s - s0[last_idx], x_DD, y_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_DDD( real_type   s,
                          real_type & x_DDD,
                          real_type & y_DDD ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.eval_DDD( s - s0[last_idx], x_DDD, y_DDD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // offset curve
  void
  ClothoidList::eval( real_type   s,
                      real_type   t,
                      real_type & x,
                      real_type & y ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.eval( s - s0[last_idx], t, x, y );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_D( real_type   s,
                        real_type   t,
                        real_type & x_D,
                        real_type & y_D ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.eval_D( s - s0[last_idx], t, x_D, y_D );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_DD( real_type   s,
                         real_type   t,
                         real_type & x_DD,
                         real_type & y_DD ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.eval_DD( s - s0[last_idx], t, x_DD, y_DD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::eval_DDD( real_type   s,
                          real_type   t,
                          real_type & x_DDD,
                          real_type & y_DDD ) const {
    findAtS( s );
    ClothoidCurve const & c = get( last_idx );
    return c.eval_DDD( s - s0[last_idx], t, x_DDD, y_DDD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::getSTK( real_type s[],
                        real_type theta[],
                        real_type kappa[] ) const {
    std::vector<ClothoidCurve>::const_iterator ic = clotoidList.begin();
    int_type  k  = 0;
    real_type ss = 0;
    while ( ic != clotoidList.end() ) {
      s[k]     = ss;
      theta[k] = ic->thetaBegin();
      kappa[k] = ic->kappaBegin();
      ss       += ic->length();
      ++k;
      ++ic;
    }
    --ic;
    s[k]     = ss;
    theta[k] = ic->thetaEnd();
    kappa[k] = ic->kappaEnd();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::getXY( real_type x[], real_type y[] ) const {
    std::vector<ClothoidCurve>::const_iterator ic = clotoidList.begin();
    int_type k  = 0;
    while ( ic != clotoidList.end() ) {
      x[k] = ic->xBegin();
      y[k] = ic->yBegin();
      ++k; ++ic;
    }
    --ic;
    x[k] = ic->xEnd();
    y[k] = ic->yEnd();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::getDeltaTheta( real_type deltaTheta[] ) const {
    std::vector<ClothoidCurve>::const_iterator ic = clotoidList.begin();
    int_type k = 0;
    for ( ++ic; ic != clotoidList.end(); ++ic, ++k ) {
      real_type tmp = ic->thetaBegin()-ic[-1].thetaEnd();
      if      ( tmp >  m_pi ) tmp -= m_2pi;
      else if ( tmp < -m_pi ) tmp += m_2pi;
      deltaTheta[k] = tmp;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::getDeltaKappa( real_type deltaKappa[] ) const {
    std::vector<ClothoidCurve>::const_iterator ic = clotoidList.begin();
    int_type k = 0;
    for ( ++ic; ic != clotoidList.end(); ++ic, ++k  )
      deltaKappa[k] = ic->kappaBegin()-ic[-1].kappaEnd();

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  real_type
  ClothoidList::closestPoint( real_type   x,
                              real_type   y,
                              real_type & X,
                              real_type & Y,
                              real_type & S ) const {
    G2LIB_ASSERT( !clotoidList.empty(), "ClothoidList::closestPoint, empty list" );
    std::vector<ClothoidCurve>::const_iterator ic = clotoidList.begin();
    std::vector<real_type>::const_iterator     is = s0.begin();
    real_type DST = ic->closestPoint( x, y, X, Y, S );
    for ( ++ic, ++is; ic != clotoidList.end(); ++ic, ++is ) {
      real_type X1, Y1, S1;
      real_type DST1 = ic->closestPoint( x, y, X1, Y1, S1 );
      if ( DST1 < DST ) {
        DST = DST1;
        X   = X1;
        Y   = Y1;
        S   = *is + S1;
      }
    }
    return DST;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int_type
  ClothoidList::findST( real_type   x,
                        real_type   y,
                        real_type & s,
                        real_type & t ) const {

    G2LIB_ASSERT( !clotoidList.empty(), "ClothoidList::findST, empty list" );
    std::vector<ClothoidCurve>::const_iterator ic = clotoidList.begin();
    std::vector<real_type>::const_iterator     is = s0.begin();

    s = t = 0;
    int_type  ipos = 0;
    int_type  iseg = 0;
    real_type S, T;
    bool ok = ic->findST( x, y, S, T );
    if ( ok ) {
      s = *is + S;
      t = T;
      iseg = 0;
    }

    for ( ++ic, ++is, ++ipos;
          ic != clotoidList.end();
          ++ic, ++is, ++ipos ) {
      bool ok1 = ic->findST( x, y, S, T );
      if ( ok && ok1 ) ok1 = std::abs(T) < std::abs(t);
      if ( ok1 ) {
        ok   = true;
        s    = *is + S;
        t    = T;
        iseg = ipos;
      }
    }

    return ok ? iseg : -(1+iseg);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int_type
  ClothoidList::findST( int_type    ibegin,
                        int_type    iend,
                        real_type   x,
                        real_type   y,
                        real_type & s,
                        real_type & t ) const {

    G2LIB_ASSERT( !clotoidList.empty(), "ClothoidList::findST, empty list" );
    G2LIB_ASSERT( ibegin >= 0 && ibegin <= iend &&
                  iend < int_type(clotoidList.size()),
                  "ClothoidList::findST( ibegin=" << ibegin << ", iend = " <<
                  iend << " , x, y, s, t ) bad range not in [0," <<
                  clotoidList.size()-1 << "]" );
    s = t = 0;
    int_type iseg = 0;
    bool ok = false;
    for ( int_type k = ibegin; k <= iend; ++k ) {
      ClothoidCurve const & ck = clotoidList[k];
      real_type S, T;
      bool ok1 = ck.findST( x, y, S, T );
      if ( ok && ok1 ) ok1 = std::abs(T) < std::abs(t);
      if ( ok1 ) {
        ok   = true;
        s    = s0[k] + S;
        t    = T;
        iseg = k;
      }
    }
    return ok ? iseg : -(1+iseg);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::rotate( real_type angle, real_type cx, real_type cy ) {
    std::vector<ClothoidCurve>::iterator ic = clotoidList.begin();
    for (; ic != clotoidList.end(); ++ic ) ic->rotate( angle, cx, cy );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::translate( real_type tx, real_type ty ) {
    std::vector<ClothoidCurve>::iterator ic = clotoidList.begin();
    for (; ic != clotoidList.end(); ++ic ) ic->translate( tx, ty );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::changeOrigin( real_type newx0, real_type newy0 ) {
    std::vector<ClothoidCurve>::iterator ic = clotoidList.begin();
    for (; ic != clotoidList.end(); ++ic ) {
      ic->changeOrigin( newx0, newy0 );
      newx0 = ic->xEnd();
      newy0 = ic->yEnd();
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::scale( real_type sfactor ) {
    std::vector<ClothoidCurve>::iterator ic = clotoidList.begin();
    real_type newx0 = ic->xBegin();
    real_type newy0 = ic->yBegin();
    for (; ic != clotoidList.end(); ++ic ) {
      ic->scale( sfactor );
      ic->changeOrigin( newx0, newy0 );
      newx0 = ic->xEnd();
      newy0 = ic->yEnd();
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::reverse() {
    std::vector<ClothoidCurve>::iterator ic = clotoidList.begin();
    for (; ic != clotoidList.end(); ++ic ) ic->reverse();
    std::reverse( clotoidList.begin(), clotoidList.end() );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::intersect( real_type                offs,
                           ClothoidCurve const &    c,
                           real_type                c_offs,
                           std::vector<real_type> & s1,
                           std::vector<real_type> & s2,
                           int_type                 max_iter,
                           real_type                tolerance ) const {
    s1.clear();
    s2.clear();
    std::vector<real_type>::const_iterator iss1, iss2;
    for ( int_type ns = 0; ns < int_type(clotoidList.size()); ++ns ) {
      std::vector<real_type> ss1, ss2;
      clotoidList[ns].intersect( offs, c, c_offs, ss1, ss2, max_iter, tolerance );
      for ( iss1 = ss1.begin();iss1 != ss1.end(); ++iss1 )
        s1.push_back( s0[ns]+(*iss1) );
      for ( iss2 = ss2.begin(); iss2 != ss2.end(); ++iss2 )
        s2.push_back( *iss2 );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::intersect( real_type                offs,
                           ClothoidList const &     CL,
                           real_type                c_offs,
                           std::vector<real_type> & s1,
                           std::vector<real_type> & s2,
                           int_type                 max_iter,
                           real_type                tolerance ) const {
    s1.clear();
    s2.clear();
    std::vector<real_type>::const_iterator iss1, iss2;
    for ( int_type ns = 0; ns < int_type(clotoidList.size()); ++ns ) {
      ClothoidCurve const & C = clotoidList[ns];
      for ( int_type ns1 = 0; ns1 < CL.numSegment(); ++ns1 ) {
        ClothoidCurve const & C1 = CL.clotoidList[ns1];
        std::vector<real_type> ss1, ss2;
        C.intersect( offs, C1, c_offs, ss1, ss2, max_iter, tolerance );
        for ( iss1 = ss1.begin();iss1 != ss1.end(); ++iss1 )
          s1.push_back( s0[ns]+(*iss1) );
        for ( iss2 = ss2.begin(); iss2 != ss2.end(); ++iss2 )
          s2.push_back( CL.s0[ns1]+*iss2 );
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::export_table( ostream_type & stream ) const {
    stream << "x\ty\ttheta0\tkappa0\tdkappa\tL\n";
    std::vector<ClothoidCurve>::const_iterator ic = clotoidList.begin();
    for (; ic != clotoidList.end(); ++ic )
      stream << ic->xBegin()     << '\t'
             << ic->yBegin()     << '\t'
             << ic->thetaBegin() << '\t'
             << ic->kappaBegin() << '\t'
             << ic->dkappa()     << '\t'
             << ic->length()     << '\n';
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  ClothoidList::export_ruby( ostream_type & stream ) const {
    stream << "data = {\n";
    std::vector<ClothoidCurve>::const_iterator ic = clotoidList.begin();
    for (; ic != clotoidList.end(); ++ic )
      stream << ic->xBegin()     << '\t'
             << ic->yBegin()     << '\t'
             << ic->thetaBegin() << '\t'
             << ic->kappaBegin() << '\t'
             << ic->dkappa()     << '\t'
             << ic->length()     << '\n';
    stream << "}\n";
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}

// EOF: ClothoidList.cc
