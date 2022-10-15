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
      m_xmin = m_ymin = m_xmax = m_ymax = 0;
    } else {
      vector<PtrBBox>::const_iterator it = bboxes.begin();

      m_xmin = (*it)->m_xmin;
      m_ymin = (*it)->m_ymin;
      m_xmax = (*it)->m_xmax;
      m_ymax = (*it)->m_ymax;

      for ( ++it; it != bboxes.end(); ++it ) {
        BBox const & currBox = **it;
        if ( currBox.m_xmin < m_xmin ) m_xmin = currBox.m_xmin;
        if ( currBox.m_ymin < m_ymin ) m_ymin = currBox.m_ymin;
        if ( currBox.m_xmax > m_xmax ) m_xmax = currBox.m_xmax;
        if ( currBox.m_ymax > m_ymax ) m_ymax = currBox.m_ymax;
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
    int_type icase = 4;
    if      ( x < m_xmin ) icase = 3;
    else if ( x > m_xmax ) icase = 5;
    if      ( y < m_ymin ) icase -= 3;
    else if ( y > m_ymax ) icase += 3;
    real_type dst = 0;
    switch ( icase ) {
      case 0: dst = hypot( x-m_xmin, y-m_ymin); break;
      case 1: dst = m_ymin-y;                   break;
      case 2: dst = hypot( x-m_xmax, y-m_ymin); break;
      case 3: dst = m_xmin-x;                   break;
      case 4:                                   break;
      case 5: dst = x-m_xmax;                   break;
      case 6: dst = hypot( x-m_xmin, y-m_ymax); break;
      case 7: dst = y-m_ymax;                   break;
      case 8: dst = hypot( x-m_xmax, y-m_ymax); break;
    }
    return dst;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  real_type
  BBox::maxDistance( real_type x, real_type y ) const {
    real_type dx = max( abs(x-m_xmin), abs(x-m_xmax) );
    real_type dy = max( abs(y-m_ymin), abs(y-m_ymax) );
    return hypot(dx,dy);
  }

}

///
/// eof: BBox.cc
///
