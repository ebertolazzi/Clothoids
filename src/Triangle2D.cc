/*
 *  Triangle-Triangle Overlap Test Routines
 *  July, 2002
 *  Updated December 2003
 *
 *  This file contains C implementation of algorithms for
 *  performing two and three-dimensional triangle-triangle intersection test
 *  The algorithms and underlying theory are described in
 *
 * "Fast and Robust Triangle-Triangle Overlap Test
 *  Using Orientation Predicates"  P. Guigue - O. Devillers
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

namespace Triangle2D {

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
  
  template <typename T>
  bool
  Triangle2D<T>::intersect( Triangle2D<T> const & t2 ) const {
    return tri_tri_intersection_2d( p1, p2, p3, t2.p1, t2.p2, t2.p3 ) ;
  }

  template <typename T>
  bool
  Triangle2D<T>::overlap( Triangle2D<T> const & t2 ) const {
    return tri_tri_overlap_test_2d( p1, p2, p3, t2.p1, t2.p2, t2.p3 ) ;
  }

  template class Triangle2D<float> ;
  template class Triangle2D<double> ;

}

///
/// eof: Triangle2D.cc
///
