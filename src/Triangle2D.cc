/*
 *  Triangle-Triangle Overlap Test Routines
 *  July, 2002
 *  Updated December 2003
 *
 *  This file contains C implementation of algorithms for
 *  performing two and three-dimensional triangle-triangle intersection test
 *  The algorithms and underlying theory are described in
 *
 *  "Fast and Robust Triangle-Triangle Overlap Test Using Orientation Predicates"
 *  P. Guigue - O. Devillers
 *
 *  Journal of Graphics Tools, 8(1), 2003
 *
 *  Several geometric predicates are defined.  Their parameters are all
 *  points.  Each point is an array of two or three double precision
 *  floating point numbers. The geometric predicates implemented in
 *  this file are:
 *
 *    int tri_tri_overlap_test_3d(p1,q1,r1,p2,q2,r2)
 *    int tri_tri_overlap_test_2d(p1,q1,r1,p2,q2,r2)
 *
 *    int tri_tri_intersection_test_3d(p1,q1,r1,p2,q2,r2,
 *                                     coplanar,source,target)
 *
 *       is a version that computes the segment of intersection when
 *       the triangles overlap (and are not coplanar)
 *
 *    each function returns 1 if the triangles (including their
 *    boundary) intersect, otherwise 0
 *
 *
 *  Other information are available from the Web page
 *  http:<i>//www.acm.org/jgt/papers/GuigueDevillers03/
 *
 */
/*
 *
 *  Two dimensional Triangle-Triangle Overlap Test
 *
 */

#include "Triangle2D.hh"
#include <algorithm>
#include <functional>

namespace G2lib {

  template <typename T>
  inline
  T
  orient_2d( T const a[2], T const b[2], T const c[2] ) {
    return (a[0]-c[0]) * (b[1]-c[1]) - (a[1]-c[1]) * (b[0]-c[0]) ;
  }

  template <typename T>
  inline
  bool
  intersection_test_vertex( T const P1[2], T const Q1[2], T const R1[2],
                            T const P2[2], T const Q2[2], T const R2[2] ) {
    if ( orient_2d(R2,P2,Q1) >= 0 ) {
      if ( orient_2d(R2,Q2,Q1) <= 0 ) {
        if ( orient_2d(P1,P2,Q1) > 0 ) {
          return orient_2d(P1,Q2,Q1) <= 0 ;
        } else {
          return orient_2d(P1,P2,R1) >= 0 && orient_2d(Q1,R1,P2) >= 0 ;
        }
      } else {
        return orient_2d(P1,Q2,Q1) <= 0 &&
               orient_2d(R2,Q2,R1) <= 0 &&
               orient_2d(Q1,R1,Q2) >= 0 ;
      }
    } else {
      if ( orient_2d(R2,P2,R1) >= 0 ) {
        if ( orient_2d(Q1,R1,R2) >= 0 ) {
          return orient_2d(P1,P2,R1) >= 0 ;
        } else {
          return orient_2d(Q1,R1,Q2) >= 0 && orient_2d(R2,R1,Q2) >= 0 ;
        }
      } else {
        return false ;
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  bool
  intersection_test_edge( T const P1[2], T const Q1[2], T const R1[2],
                          T const P2[2], T const Q2[2], T const R2[2] ) {
    if ( orient_2d(R2,P2,Q1) >= 0 ) {
      if ( orient_2d(P1,P2,Q1) >= 0 ) {
        return orient_2d(P1,Q1,R2) >= 0 ;
      } else {
        return orient_2d(Q1,R1,P2) >= 0 && orient_2d(R1,P1,P2) >= 0 ;
      }
    } else if ( orient_2d(R2,P2,R1) >= 0 ) {
      return orient_2d(P1,P2,R1) >= 0 &&
             ( orient_2d(P1,R1,R2) >= 0 || orient_2d(Q1,R1,R2) >= 0 ) ;
    }
    return false ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  bool
  tri_tri_intersection_2d( T const p1[2], T const q1[2], T const r1[2],
                           T const p2[2], T const q2[2], T const r2[2] ) {
    if ( orient_2d(p2,q2,p1) >= 0 ) {
      if ( orient_2d(q2,r2,p1) >= 0 ) {
        return orient_2d(r2,p2,p1) >= 0 || intersection_test_edge(p1,q1,r1,p2,q2,r2) ;
      } else {
        if ( orient_2d(r2,p2,p1) >= 0 ) return intersection_test_edge(p1,q1,r1,r2,p2,q2) ;
        else                            return intersection_test_vertex(p1,q1,r1,p2,q2,r2) ;
      }
    } else {
      if ( orient_2d(q2,r2,p1) >= 0 ) {
        if ( orient_2d(r2,p2,p1) >= 0 ) return intersection_test_edge(p1,q1,r1,q2,r2,p2) ;
        else                            return intersection_test_vertex(p1,q1,r1,q2,r2,p2) ;
      } else {
        return intersection_test_vertex(p1,q1,r1,r2,p2,q2) ;
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  bool
  tri_tri_overlap_test_2d( T const p1[2], T const q1[2], T const r1[2],
                           T const p2[2], T const q2[2], T const r2[2] ) {
    if ( orient_2d(p1,q1,r1) < 0 ) {
      if ( orient_2d(p2,q2,r2) < 0 ) return tri_tri_intersection_2d(p1,r1,q1,p2,r2,q2) ;
      else                           return tri_tri_intersection_2d(p1,r1,q1,p2,q2,r2) ;
    } else {
      if ( orient_2d(p2,q2,r2) < 0 ) return tri_tri_intersection_2d(p1,q1,r1,p2,r2,q2) ;
      else                           return tri_tri_intersection_2d(p1,q1,r1,p2,q2,r2) ;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Triangle2D::maxmin3( real_type   a,
                       real_type   b,
                       real_type   c,
                       real_type & vmin,
                       real_type & vmax) const {
    vmin = vmax = a ;
    if ( b < vmin ) vmin = b ;
    else            vmax = b ;
    if ( c < vmin ) vmin = c ;
    else if ( c > vmax ) vmax = c ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Triangle2D::AABBTree::~AABBTree() {
    switch ( numChildren ) {
    case 0:
      if( data.pTriangle != nullptr ) delete data.pTriangle ; data.pTriangle = nullptr ;
      break ;
    case 2: if( data.pChildren[1] != nullptr ) delete data.pChildren[1] ;
    case 1: if( data.pChildren[0] != nullptr ) delete data.pChildren[0] ;
      data.pChildren[1] = data.pChildren[0] = nullptr ;
      break ;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Triangle2D::AABBTree::AABBTree( std::vector<Triangle2D> & triangles ) {
    std::vector<Triangle2D const *> pTriangles ;
    pTriangles.reserve( triangles.size() ) ;
    typename std::vector<Triangle2D>::const_iterator it = triangles.begin() ;
    for ( ; it != triangles.end() ; ++it )
      pTriangles.push_back( &*it ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  class AABBcomparatorX {
    real_type cutPos ;
  public:
    AABBcomparatorX( real_type const & cp ) : cutPos(cp) {}
    bool operator () ( Triangle2D const * pT ) const
    { return pT->baricenterX() < cutPos; }
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  class AABBcomparatorY {
    real_type cutPos ;
  public:
    AABBcomparatorY( real_type const & cp ) : cutPos(cp) {}
    bool operator () ( Triangle2D const * pT ) const
    { return pT->baricenterY() < cutPos; }
  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Triangle2D::AABBTree::AABBTree(
    typename std::vector<Triangle2D const *>::iterator & begin,
    typename std::vector<Triangle2D const *>::iterator & end
  ) {
    numChildren = -1 ;
    xmin = ymin = xmax = ymax = 0 ;
    data.pChildren[0] =
    data.pChildren[1] = nullptr ;

    if ( begin == end ) return;

    (*begin)->bbox(xmin,ymin,xmax,ymax) ;
    if ( end - begin == 1 ) {
      numChildren = 0 ;
      data.pTriangle = new Triangle2D(**begin) ;
      return;
    }

    typename std::vector<Triangle2D const *>::iterator it = begin ;
    for ( ++it ; it != end ; ++it ) {
      real_type xmi, ymi, xma, yma ;
      (*it)->bbox( xmi, ymi, xma, yma ) ;
      if ( xmi < xmin ) xmin = xmi ;
      if ( xma > xmax ) xmax = xma ;
      if ( ymi < ymin ) ymin = ymi ;
      if ( yma > ymax ) ymax = yma ;
    }

    if ( (ymax - ymin) > (xmax - xmin) ) {
      // cut along Y
      AABBcomparatorY comp( (ymax + ymin)/2 ) ;
      it = std::partition( begin, end, comp );
    } else {
      // cut along X
      AABBcomparatorX comp( (xmax + xmin)/2 ) ;
      it = std::partition( begin, end, comp );

    }
    if ( it - begin > 0 ) {
      data.pChildren[0] = new AABBTree(begin, it);
      numChildren = 1 ;
    }
    if ( end - it > 0 ) {
      data.pChildren[numChildren] = new AABBTree(it, end);
      ++numChildren ;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Triangle2D::AABBTree::AABBTree( Triangle2D const & triangle ) {
    numChildren = 0 ;
    triangle.bbox( xmin, ymin, xmax, ymax ) ;
    data.pTriangle = new Triangle2D(triangle) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Triangle2D::AABBTree::AABBTree(
    AABBTree   const * pTree,
    Triangle2D const & triangle
  ) {
    numChildren = 2 ;
    pTree->bbox( xmin, ymin, xmax, ymax ) ;
    real_type xmi, ymi, xma, yma ;
    triangle.bbox( xmi, ymi, xma, yma ) ;
    if ( xmi < xmin ) xmin = xmi ;
    if ( xma > xmax ) xmax = xma ;
    if ( ymi < ymin ) ymin = ymi ;
    if ( yma > ymax ) ymax = yma ;
    data.pChildren[0] = pTree ;
    data.pChildren[1] = new AABBTree(triangle) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Triangle2D::AABBTree::AABBTree(
    AABBTree const * pTreeL,
    AABBTree const * pTreeR
  ) {
    numChildren = 2 ;
    pTreeL->bbox( xmin, ymin, xmax, ymax ) ;
    real_type xmi, ymi, xma, yma ;
    pTreeR->bbox( xmi, ymi, xma, yma ) ;
    if ( xmi < xmin ) xmin = xmi ;
    if ( xma > xmax ) xmax = xma ;
    if ( ymi < ymin ) ymin = ymi ;
    if ( yma > ymax ) ymax = yma ;
    data.pChildren[0] = pTreeL ;
    data.pChildren[1] = pTreeR ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  Triangle2D::AABBTree::overlap( Triangle2D const & triangle ) const {
    real_type xmi, ymi, xma, yma ;
    triangle.bbox( xmi, ymi, xma, yma ) ;
    if ( xmax < xmi ) return false; // a is left of b
    if ( xmin > xma ) return false; // a is right of b
    if ( ymax < ymi ) return false; // a is above b
    if ( ymin > yma ) return false; // a is below b
    // second level check
    switch ( numChildren ) {
    case 0:
      return data.pTriangle->overlap( triangle ) ;
    case 1:
      return data.pChildren[0]->overlap( triangle ) ;
    case 2:
      return data.pChildren[0]->overlap( triangle ) ||
             data.pChildren[1]->overlap( triangle ) ;
    }
    return true; // unused
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  Triangle2D::AABBTree::overlap( AABBTree const * pTree ) const {
    real_type xmi, ymi, xma, yma ;
    pTree->bbox( xmi, ymi, xma, yma ) ;
    if ( xmax < xmi ) return false; // a is left of b
    if ( xmin > xma ) return false; // a is right of b
    if ( ymax < ymi ) return false; // a is above b
    if ( ymin > yma ) return false; // a is below b
    // second level check
    switch ( numChildren ) {
    case 0:
      return pTree->overlap(*data.pTriangle) ;
    case 1:
      return data.pChildren[0]->overlap( pTree ) ;
    case 2:
      return data.pChildren[0]->overlap( pTree ) ||
             data.pChildren[1]->overlap( pTree ) ;
    }
    return true; // unused
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  Triangle2D::intersect( Triangle2D const & t2 ) const {
    return tri_tri_intersection_2d( p1, p2, p3, t2.p1, t2.p2, t2.p3 ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool
  Triangle2D::overlap( Triangle2D const & t2 ) const {
    return tri_tri_overlap_test_2d( p1, p2, p3, t2.p1, t2.p2, t2.p3 ) ;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Triangle2D::distMinMax( real_type   x,
                          real_type   y,
                          real_type & dmin,
                          real_type & dmax ) const {

    int_type in = isInside( x, y ) ;
    if ( in >= 0 ) { dmin = dmax = 0 ; return ; }

    real_type d1 = hypot( x-p1[0], y-p1[1] ) ;
    real_type d2 = hypot( x-p2[0], y-p2[1] ) ;
    real_type d3 = hypot( x-p3[0], y-p3[1] ) ;
    real_type const * P1 = this->p1 ;
    real_type const * P2 = this->p2 ;
    real_type const * P3 = this->p3 ;
    if ( d1 < d2 ) { std::swap( d1, d2 ) ; std::swap( P1, P2 ) ; }
    if ( d2 < d3 ) { std::swap( d2, d3 ) ; std::swap( P2, P3 ) ; }
    if ( d1 < d2 ) { std::swap( d1, d2 ) ; std::swap( P1, P2 ) ; }

    dmax = d1 ;

    real_type dx  = x     - P2[0] ;
    real_type dy  = y     - P2[1] ;
    real_type dx1 = P3[0] - P2[0] ;
    real_type dy1 = P3[1] - P2[1] ;
    real_type tmp = dx * dx1 + dy * dy1 ;
    if ( tmp < 0 ) {
      dmin = d2 ;
    } else {
      real_type tmp2 = hypot(dx,dy)*hypot(dx1,dy1) ;
      if ( tmp > tmp2 ) {
        dmin = d3 ;
      } else {
        real_type S = tmp/tmp2 ;
        real_type X = P2[0] + S*dx1 ;
        real_type Y = P2[1] + S*dy1 ;
        dmin = hypot( x-X, y-Y ) ;
      }
    }
  }

}

///
/// eof: Triangle2D.cc
///
