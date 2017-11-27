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

#include <cmath>
#include <cfloat>
#include <algorithm>

#ifdef __GCC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

namespace Clothoid {

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

  void
  ClothoidList::copy( ClothoidList const & L ) {
    s0.resize( L.s0.size() ) ;
    std::copy( L.s0.begin(), L.s0.end(), s0.begin() ) ;
    clotoidList.resize( L.clotoidList.size() ) ;
    std::copy( L.clotoidList.begin(), L.clotoidList.end(), clotoidList.begin() ) ;
  }

  void
  ClothoidList::reserve( indexType n ) {
    s0.reserve(size_t(n+1)) ;
    clotoidList.reserve(size_t(n)) ;
  }

  void
  ClothoidList::add( ClothoidCurve const & c ) {
    if ( clotoidList.empty() ) {
      s0.push_back(0) ;
      s0.push_back(c.getL()) ;
    } else {
      s0.push_back(s0.back()+c.getL()) ;
    }
    clotoidList.push_back(c) ;
  }

  ClothoidCurve const &
  ClothoidList::get( indexType idx ) const {
    G2LIB_ASSERT( !clotoidList.empty(), "ClothoidList::get( " << idx << " ) empty list" );
    G2LIB_ASSERT( idx >= 0 && idx < indexType(clotoidList.size()),
                  "ClothoidList::get( " << idx << " ) bad index" );
    return clotoidList[idx];
  }

  ClothoidCurve const &
  ClothoidList::getAtS( valueType s, indexType & last_idx ) const {
    findAtS(s,last_idx);
    return get(last_idx) ;
  }

  bool
  ClothoidList::findAtS( valueType s, indexType & last_idx ) const {
    indexType ns = indexType(clotoidList.size()) ;
    G2LIB_ASSERT( last_idx >= 0 && last_idx < ns,
                  "ClothoidList::findAtS( " << s << ", " << last_idx <<
                  " ) bad index" );
    valueType const * sL = &s0[last_idx] ;
    if ( s < sL[0] ) {
      valueType const * sB = &s0.front() ;
      last_idx = indexType(std::lower_bound( sB, sL, s )-sB) ;
    } else if ( s > sL[1] ) {
      valueType const * sE = &s0[ns+1] ;
      last_idx += indexType(std::lower_bound( sL, sE, s )-sL) ;
    } else {
      return true ; // vale intervallo precedente
    }
    if ( s0[last_idx] > s ) --last_idx ; // aggiustamento caso di bordo
    return true ;
  }

  valueType
  ClothoidList::theta( valueType s, indexType & last_idx ) const {
    findAtS( s, last_idx );
    ClothoidCurve const & c = get( last_idx );
    return c.theta( s - s0[last_idx] ) ;
  }

  valueType
  ClothoidList::theta_D( valueType s, indexType & last_idx ) const {
    findAtS( s, last_idx );
    ClothoidCurve const & c = get( last_idx );
    return c.theta_D( s - s0[last_idx] ) ;
  }

  valueType
  ClothoidList::theta_DD( valueType s, indexType & last_idx ) const {
    findAtS( s, last_idx );
    ClothoidCurve const & c = get( last_idx );
    return c.theta_DD( s - s0[last_idx] ) ;
  }

  valueType
  ClothoidList::X( valueType s, indexType & last_idx ) const {
    findAtS( s, last_idx );
    ClothoidCurve const & c = get( last_idx );
    return c.X( s - s0[last_idx] ) ;
  }

  valueType
  ClothoidList::Y( valueType s, indexType & last_idx ) const {
    findAtS( s, last_idx );
    ClothoidCurve const & c = get( last_idx );
    return c.Y( s - s0[last_idx] ) ;
  }

  void
  ClothoidList::eval( valueType   s,
                      indexType & last_idx,
                      valueType & theta,
                      valueType & kappa,
                      valueType & x,
                      valueType & y ) const {
    findAtS( s, last_idx );
    ClothoidCurve const & c = get( last_idx );
    return c.eval( s - s0[last_idx], theta, kappa, x, y ) ;
  }

  void
  ClothoidList::eval( valueType   s,
                      indexType & last_idx,
                      valueType & x,
                      valueType & y ) const {
    findAtS( s, last_idx );
    ClothoidCurve const & c = get( last_idx );
    return c.eval( s - s0[last_idx], x, y ) ;
  }

  void
  ClothoidList::eval_D( valueType   s,
                        indexType & last_idx,
                        valueType & x_D,
                        valueType & y_D ) const {
    findAtS( s, last_idx );
    ClothoidCurve const & c = get( last_idx );
    return c.eval_D( s - s0[last_idx], x_D, y_D ) ;
  }

  void
  ClothoidList::eval_DD( valueType   s,
                         indexType & last_idx,
                         valueType & x_DD,
                         valueType & y_DD ) const {
    findAtS( s, last_idx );
    ClothoidCurve const & c = get( last_idx );
    return c.eval_DD( s - s0[last_idx], x_DD, y_DD ) ;
  }

  void
  ClothoidList::eval_DDD( valueType   s,
                          indexType & last_idx,
                          valueType & x_DDD,
                          valueType & y_DDD ) const {
    findAtS( s, last_idx );
    ClothoidCurve const & c = get( last_idx );
    return c.eval_DDD( s - s0[last_idx], x_DDD, y_DDD ) ;
  }

  // offset curve
  void
  ClothoidList::eval( valueType   s,
                      indexType & last_idx,
                      valueType   offs,
                      valueType & x,
                      valueType & y ) const {
    findAtS( s, last_idx );
    ClothoidCurve const & c = get( last_idx );
    return c.eval( s - s0[last_idx], offs, x, y ) ;
  }

  void
  ClothoidList::eval_D( valueType   s,
                        indexType & last_idx,
                        valueType   offs,
                        valueType & x_D,
                        valueType & y_D ) const {
    findAtS( s, last_idx );
    ClothoidCurve const & c = get( last_idx );
    return c.eval_D( s - s0[last_idx], offs, x_D, y_D ) ;
  }

  void
  ClothoidList::eval_DD( valueType   s,
                         indexType & last_idx,
                         valueType   offs,
                         valueType & x_DD,
                         valueType & y_DD ) const {
    findAtS( s, last_idx );
    ClothoidCurve const & c = get( last_idx );
    return c.eval_DD( s - s0[last_idx], offs, x_DD, y_DD ) ;
  }

  void
  ClothoidList::eval_DDD( valueType   s,
                          indexType & last_idx,
                          valueType   offs,
                          valueType & x_DDD,
                          valueType & y_DDD ) const {
    findAtS( s, last_idx );
    ClothoidCurve const & c = get( last_idx );
    return c.eval_DDD( s - s0[last_idx], offs, x_DDD, y_DDD ) ;
  }

  valueType
  ClothoidList::closestPoint( valueType   x,
                              valueType   y,
                              valueType   ds,
                              valueType & X,
                              valueType & Y,
                              valueType & S ) const {
    G2LIB_ASSERT( !clotoidList.empty(), "ClothoidList::closestPoint, empty list" );
    std::vector<ClothoidCurve>::const_iterator ic = clotoidList.begin() ;
    std::vector<valueType>::const_iterator     is = s0.begin() ;
    valueType DST = ic->closestPoint( x, y, ds, X, Y, S );
    for ( ++ic, ++is ; ic != clotoidList.end() ; ++ic, ++is ) {
      valueType X1, Y1, S1 ;
      valueType DST1 = ic->closestPoint( x, y, ds, X1, Y1, S1 );
      if ( DST1 < DST ) {
        DST = DST1 ;
        X   = X1;
        Y   = Y1;
        S   = *is + S1;
      }
    }
    return DST ;
  }

  void
  ClothoidList::rotate( valueType angle, valueType cx, valueType cy ) {
    std::vector<ClothoidCurve>::iterator ic = clotoidList.begin() ;
    for ( ; ic != clotoidList.end() ; ++ic ) ic->rotate( angle, cx, cy ) ;
  }

  void
  ClothoidList::translate( valueType tx, valueType ty ) {
    std::vector<ClothoidCurve>::iterator ic = clotoidList.begin() ;
    for ( ; ic != clotoidList.end() ; ++ic ) ic->translate( tx, ty ) ;
  }

  void
  ClothoidList::moveOrigin( valueType newx0, valueType newy0 ) {
    std::vector<ClothoidCurve>::iterator ic = clotoidList.begin() ;
    for ( ; ic != clotoidList.end() ; ++ic ) {
      ic->moveOrigin( newx0, newy0 ) ;
      newx0 = ic->Xend() ;
      newy0 = ic->Yend() ;
    }
  }

  void
  ClothoidList::scale( valueType sfactor ) {
    std::vector<ClothoidCurve>::iterator ic = clotoidList.begin() ;
    valueType newx0 = ic->Xbegin() ;
    valueType newy0 = ic->Ybegin() ;
    for ( ; ic != clotoidList.end() ; ++ic ) {
      ic->scale( sfactor ) ;
      ic->moveOrigin( newx0, newy0 ) ;
      newx0 = ic->Xend() ;
      newy0 = ic->Yend() ;
    }
  }

  void
  ClothoidList::reverse() {
    std::vector<ClothoidCurve>::iterator ic = clotoidList.begin() ;
    for ( ; ic != clotoidList.end() ; ++ic ) ic->reverse() ;
    std::reverse( clotoidList.begin(), clotoidList.end() );
  }
}

// EOF: ClothoidList.cc
