/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2018                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Paolo Bevilacqua and Enrico Bertolazzi                              |
 |                                                                          |
 |      (1) Dipartimento di Ingegneria e Scienza dell'Informazione          |
 |      (2) Dipartimento di Ingegneria Industriale                          |
 |                                                                          |
 |      Universita` degli Studi di Trento                                   |
 |      email: paolo.bevilacqua@unitn.it                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: BBox.cc
///

#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

// Workaround for Visual Studio
#ifdef min
  #undef min
#endif

#ifdef max
  #undef max
#endif

#include <algorithm>

namespace G2lib {

  using std::abs;
  using std::min;
  using std::max;

  /*\
   |   ____  ____
   |  | __ )| __ )  _____  __
   |  |  _ \|  _ \ / _ \ \/ /
   |  | |_) | |_) | (_) >  <
   |  |____/|____/ \___/_/\_\
  \*/

  void
  BBox::join( vector<PtrBBox> const & bboxes ) {
    if ( bboxes.empty() ) {
      std::fill_n( m_bbox, 4, 0 );
    } else {
      this->x_min() = Utils::Inf<real_type>();
      this->y_min() = Utils::Inf<real_type>();
      this->x_max() = -Utils::Inf<real_type>();
      this->y_max() = -Utils::Inf<real_type>();

      for ( auto const & it : bboxes ) {
        if ( it->x_min() < this->x_min() ) this->x_min() = it->x_min();
        if ( it->y_min() < this->y_min() ) this->y_min() = it->y_min();
        if ( it->x_max() > this->x_max() ) this->x_max() = it->x_max();
        if ( it->y_max() > this->y_max() ) this->y_max() = it->y_max();
      }
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  real_type
  BBox::distance( real_type x, real_type y ) const {
    /*\
     |
     |   6          7          8
     |       +-------------+
     |       |             |
     |   3   |      4      |   5
     |       |             |
     |       +-------------+
     |   0          1          2
     |
    \*/
    integer icase = 4;
    if      ( x < x_min() ) icase = 3;
    else if ( x > x_max() ) icase = 5;
    if      ( y < y_min() ) icase -= 3;
    else if ( y > y_max() ) icase += 3;
    real_type dst{0};
    switch ( icase ) {
      case 0: dst = hypot( x-x_min(), y-y_min() ); break;
      case 1: dst = y_min()-y;                     break;
      case 2: dst = hypot( x-x_max(), y-y_min() ); break;
      case 3: dst = x_min()-x;                     break;
      case 4:                                      break;
      case 5: dst = x-x_max();                     break;
      case 6: dst = hypot( x-x_min(), y-y_max() ); break;
      case 7: dst = y-y_max();                     break;
      case 8: dst = hypot( x-x_max(), y-y_max() ); break;
    }
    return dst;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  real_type
  BBox::max_distance( real_type x, real_type y ) const {
    real_type dx = max( abs(x-x_min()), abs(x-x_max()) );
    real_type dy = max( abs(y-y_min()), abs(y-y_max()) );
    return hypot(dx,dy);
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  void
  BBox::print( ostream_type & stream ) const {
    fmt::print( stream,
      "BBOX (xmin,ymin,xmax,ymax) = ( {}, {}, {}, {} )\n",
      x_min(), y_min(), x_max(), y_max()
    );
  }
}

///
/// eof: BBox.cc
///
