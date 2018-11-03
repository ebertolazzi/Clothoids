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
/// file: AABBtree.hh
///

#pragma once

#ifndef AABBTREE_HH
#define AABBTREE_HH

#include "G2lib.hh"

#include <vector>
#include <iomanip>
#include <utility> // pair

#ifdef G2LIB_USE_CXX11
#include <memory>  // shared_ptr
#endif

namespace G2lib {

  using std::setw;
  using std::vector;
  using std::pair;

  #ifdef G2LIB_USE_CXX11
  using std::make_shared;
  using std::shared_ptr; // promemoria shared_ptr<Foo>(&foo, [](void*){});
  #endif

  class AABBtree;

  /*\
   |   ____  ____
   |  | __ )| __ )  _____  __
   |  |  _ \|  _ \ / _ \ \/ /
   |  | |_) | |_) | (_) >  <
   |  |____/|____/ \___/_/\_\
  \*/
  class BBox {
  public:
    #ifdef G2LIB_USE_CXX11
    typedef shared_ptr<BBox const> PtrBBox;
    #else
    typedef BBox const * PtrBBox;
    #endif

  private:
    real_type xmin, ymin, xmax, ymax;
    int_type  id;
    int_type  ipos;
    BBox();
    BBox( BBox const & );

  public:

    BBox( real_type _xmin,
          real_type _ymin,
          real_type _xmax,
          real_type _ymax,
          int_type  _id,
          int_type  _ipos ){
      this -> xmin = _xmin;
      this -> ymin = _ymin;
      this -> xmax = _xmax;
      this -> ymax = _ymax;
      this -> id   = _id;
      this -> ipos = _ipos;
    }

    BBox( vector<PtrBBox> const & bboxes,
          int_type                _id,
          int_type                _ipos ) {
      this -> id   = _id;
      this -> ipos = _ipos;
      this -> join( bboxes );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    real_type Xmin() const { return xmin; }
    real_type Ymin() const { return ymin; }
    real_type Xmax() const { return xmax; }
    real_type Ymax() const { return ymax; }

    int_type const & Id()   const { return id; }
    int_type const & Ipos() const { return ipos; }

    BBox const &
    operator = ( BBox const & rhs ) {
      this -> xmin = rhs.xmin;
      this -> ymin = rhs.ymin;
      this -> xmax = rhs.xmax;
      this -> ymax = rhs.ymax;
      this -> id   = rhs.id;
      this -> ipos = rhs.ipos;
      return *this;
    }

    bool
    collision( BBox const & box ) const {
      return !( (box.xmin > xmax ) ||
                (box.xmax < xmin ) ||
                (box.ymin > ymax ) ||
                (box.ymax < ymin ) );
    }

    void
    join( vector<PtrBBox> const & bboxes );

    real_type
    distance( real_type x, real_type y ) const;

    real_type
    maxDistance( real_type x, real_type y ) const;

    void
    print( ostream_type & stream ) const {
      stream
        << "BBOX (xmin,ymin,xmax,ymax) = ( " << xmin
        << ", " << ymin << ", " << xmax << ", " << ymax
        << " )\n";
    }

    friend class AABBtree;
  };

  inline
  ostream_type &
  operator << ( ostream_type & stream, BBox const & bb ) {
    bb.print(stream);
    return stream;
  }

  /*\
   |      _        _    ____  ____  _
   |     / \      / \  | __ )| __ )| |_ _ __ ___  ___
   |    / _ \    / _ \ |  _ \|  _ \| __| '__/ _ \/ _ \
   |   / ___ \  / ___ \| |_) | |_) | |_| | |  __/  __/
   |  /_/   \_\/_/   \_\____/|____/ \__|_|  \___|\___|
  \*/

  class AABBtree {
  public:

  #ifdef G2LIB_USE_CXX11
    typedef shared_ptr<BBox const> PtrBBox;
    typedef shared_ptr<AABBtree>   PtrAABB;
  #else
    typedef BBox const *           PtrBBox;
    typedef AABBtree *             PtrAABB;
  #endif

  typedef pair<PtrBBox,PtrBBox> PairPtrBBox;
  typedef vector<PtrBBox>       VecPtrBBox;
  typedef vector<PairPtrBBox>   VecPairPtrBBox;

  private:

    // bbox of the tree
    PtrBBox         pBBox;
    vector<PtrAABB> children;

    AABBtree( AABBtree const & tree );

    /*!
     | Compute the minimum of the maximum distance
     | between a point
    \*/
    static
    real_type
    min_maxdist( real_type        x,
                 real_type        y,
                 AABBtree const & tree,
                 real_type        mmDist );

    /*!
     | Select the candidate which bbox have distance less than mmDist
    \*/
    static
    void
    min_maxdist_select( real_type        x,
                        real_type        y,
                        real_type        mmDist,
                        AABBtree const & tree,
                        VecPtrBBox     & candidateList );

  public:

    AABBtree();
    ~AABBtree();

    // copy contructor (recursive)
    //AABBtree( PtrBBox const & pbox ) {
    //  this -> pBBox = pbox;
    //  children.clear();
    //}

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void clear();

    bool empty() const;

    void
    bbox( real_type & xmin,
          real_type & ymin,
          real_type & xmax,
          real_type & ymax ) const {
      xmin = pBBox->xmin;
      ymin = pBBox->ymin;
      xmax = pBBox->xmax;
      ymax = pBBox->ymax;
    }

    void
    build( vector<PtrBBox> const & bboxes );

    void
    print( ostream_type & stream, int level = 0 ) const;

    template <typename COLLISION_fun>
    bool
    collision( AABBtree const & tree,
               COLLISION_fun    ifun,
               bool             swap_tree = false ) const {

      // check bbox with
      if ( !tree.pBBox->collision(*pBBox) ) return false;

      int icase = (children.empty() ? 0 : 1) +
                  (tree.children.empty()? 0 : 2);

      switch ( icase ) {
      case 0: // both leaf, use GeomPrimitive intersection algorithm
        if ( swap_tree ) return ifun( tree.pBBox, pBBox );
        else             return ifun( pBBox, tree.pBBox );
      case 1: // first is a tree, second is a leaf
        { typename vector<PtrAABB>::const_iterator it;
          for ( it = children.begin(); it != children.end(); ++it )
            if ( tree.collision( **it, ifun, !swap_tree ) )
              return true;
        }
        break;
      case 2: // first leaf, second is a tree
        { typename vector<PtrAABB>::const_iterator it;
          for ( it = tree.children.begin();
                it != tree.children.end(); ++it )
            if ( this->collision( **it, ifun, swap_tree ) )
              return true;
        }
        break;
      case 3: // first is a tree, second is a tree
        { typename vector<PtrAABB>::const_iterator c1;
          typename vector<PtrAABB>::const_iterator c2;
          for ( c1 = children.begin(); c1 != children.end(); ++c1 )
            for ( c2 = tree.children.begin();
                  c2 != tree.children.end(); ++c2 )
              if ( (*c1)->collision( **c2, ifun, swap_tree ) )
                return true;
        }
        break;
      }
      return false;
    }

    void
    intersect(
      AABBtree const & tree,
      VecPairPtrBBox & intersectionList,
      bool             swap_tree = false
    ) const;

    void
    min_distance( real_type    x,
                  real_type    y,
                  VecPtrBBox & candidateList ) const;

  };

}

#endif

///
/// eof: AABBtree.hh
///
