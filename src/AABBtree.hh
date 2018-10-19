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
#include <memory>   // shared_ptr
#include <iomanip>
#include <utility>  // pair

namespace G2lib {

  using std::setw;
  using std::shared_ptr;
  using std::vector;
  using std::pair;
  using std::make_shared;

  template <typename GeomPrimitive> class AABBtree;

  /*\
   |   ____  ____
   |  | __ )| __ )  _____  __
   |  |  _ \|  _ \ / _ \ \/ /
   |  | |_) | |_) | (_) >  <
   |  |____/|____/ \___/_/\_\
  \*/

  template <typename GeomPrimitive>
  class BBox {
  public:

    typedef shared_ptr<BBox>          PtrBBox;
    typedef shared_ptr<GeomPrimitive> PtrGeom;

  private:

    real_type xmin, ymin, xmax, ymax;
    PtrGeom   pGeom;

  public:

    BBox( real_type _xmin = 0,
          real_type _ymin = 0,
          real_type _xmax = 0,
          real_type _ymax = 0,
          PtrGeom   _prim = PtrGeom() ){
      this -> xmin  = _xmin;
      this -> ymin  = _ymin;
      this -> xmax  = _xmax;
      this -> ymax  = _ymax;
      this -> pGeom = _prim;
    }

    BBox( vector<PtrBBox> const & bboxes ) {
      this -> pGeom = PtrGeom();
      this -> join( bboxes );
    }

    BBox( PtrGeom _prim ) {
      this -> pGeom = _prim;
      _prim->bbox( xmin, ymin, xmax, ymax );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    real_type Xmin() const { return xmin; }
    real_type Ymin() const { return ymin; }
    real_type Xmax() const { return xmax; }
    real_type Ymax() const { return ymax; }

    PtrGeom const getPointerGeom() const { return pGeom; }

    BBox<GeomPrimitive> const &
    operator = ( BBox<GeomPrimitive> const & rhs ) {
      this -> xmin = rhs.xmin;
      this -> ymin = rhs.ymin;
      this -> xmax = rhs.xmax;
      this -> ymax = rhs.ymax;
      this -> prim = rhs.prim;
      return *this;
    }

    bool
    intersect( BBox const & box ) const {
      return !( (box.xmin > xmax ) ||
                (box.xmax < xmin ) ||
                (box.ymin > ymax ) ||
                (box.ymax < ymin ) );
    }

    void
    join( vector<PtrBBox> const & bboxes ) {
      if ( bboxes.empty() ) {
        xmin = ymin = xmax = ymax = 0 ;
      } else {
        typename vector<PtrBBox>::const_iterator it = bboxes.begin();
        PtrBBox const box = *it;

        xmin = box->xmin;
        ymin = box->ymin;
        xmax = box->xmax;
        ymax = box->ymax;

        for ( ++it; it != bboxes.end(); ++it ) {
          PtrBBox const currBox = *it;
          if ( currBox->xmin < xmin ) xmin = currBox->xmin;
          if ( currBox->ymin < ymin ) ymin = currBox->ymin;
          if ( currBox->xmax > xmax ) xmax = currBox->xmax;
          if ( currBox->ymax > ymax ) ymax = currBox->ymax;
        }
      }
    }

    void
    print( ostream_type & stream ) const {
      stream
        << "BBOX xmin = " << setw(12) << xmin
        << " ymin = "     << setw(12) << ymin
        << " xmax = "     << setw(12) << xmax
        << " ymax = "     << setw(12) << ymax
        << "\n";
    }

    friend class AABBtree<GeomPrimitive>;
  };

  /*\
   |      _        _    ____  ____  _
   |     / \      / \  | __ )| __ )| |_ _ __ ___  ___
   |    / _ \    / _ \ |  _ \|  _ \| __| '__/ _ \/ _ \
   |   / ___ \  / ___ \| |_) | |_) | |_| | |  __/  __/
   |  /_/   \_\/_/   \_\____/|____/ \__|_|  \___|\___|
  \*/

  template <typename GeomPrimitive>
  class AABBtree {
  public:
    typedef shared_ptr<GeomPrimitive>            PtrGeom;
    typedef shared_ptr<BBox<GeomPrimitive> >     PtrBBox;
    typedef shared_ptr<AABBtree<GeomPrimitive> > PtrAABB;
    typedef pair<PtrBBox,PtrBBox>                PairPtrBBox;
    typedef vector<PairPtrBBox>                  VecPairPtrBBox;

  private:

    // bbox of the tree
    PtrBBox         pBBox;
    vector<PtrAABB> children;

  public:

    AABBtree() {
      pBBox.reset();
      children.clear();
    }

    // copy contructor (recursive)
    AABBtree( PtrBBox const & pbox ) {
      this -> pBBox = pbox;
      children.clear();
    }

    // copy contructor (recursive)
    AABBtree( AABBtree const & tree ) {
      this -> pBBox = tree.pBBox;
      children.reserve(tree.children.size());
      typename PtrAABB::const_iterator t;
      for ( t = tree.children.begin(); t != tree.children.end(); ++t )
        children.push_back( PtrAABB( new AABBtree<GeomPrimitive>(*t) ) );
    }

    virtual
    ~AABBtree() {
      pBBox.reset();
      children.clear();
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    clear() {
      // empty if an old tree present
      pBBox.reset();
      children.clear();
    }

    bool
    empty() const {
      return children.empty() && !pBBox;
    }

    void
    build( vector<PtrBBox> const & bboxes ) {

      if ( bboxes.empty() ) return;

      size_t size = bboxes.size();

      if ( size == 1 ) {
        this -> pBBox = bboxes.front();
        return;
      }

      pBBox = make_shared<BBox<GeomPrimitive> >( bboxes );

      real_type xmin = pBBox -> xmin;
      real_type ymin = pBBox -> ymin;
      real_type xmax = pBBox -> xmax;
      real_type ymax = pBBox -> ymax;

      vector<PtrBBox> posBoxes;
      vector<PtrBBox> negBoxes;

      if ( (ymax - ymin) > (xmax - xmin) ) {
        real_type cutPos = (ymax + ymin)/2;
        typename vector<PtrBBox>::const_iterator it;
        for ( it = bboxes.begin(); it != bboxes.end(); ++it ) {
          real_type ymid = ( (*it) -> ymin + (*it) -> ymax ) / 2;
          if ( ymid > cutPos ) posBoxes.push_back(*it);
          else                 negBoxes.push_back(*it);
        }
      } else {
        real_type cutPos = (xmax + xmin)/2;
        typename vector<PtrBBox>::const_iterator it;
        for ( it = bboxes.begin(); it != bboxes.end(); ++it ) {
          real_type xmid = ( (*it) -> xmin + (*it) -> xmax ) / 2;
          if ( xmid > cutPos ) posBoxes.push_back(*it);
          else                 negBoxes.push_back(*it);
        }
      }

      if ( negBoxes.empty() ) {
        typename vector<PtrBBox>::iterator midIdx;
        midIdx = posBoxes.begin() + posBoxes.size()/2;
        negBoxes.insert( negBoxes.end(), midIdx, posBoxes.end() );
        posBoxes.erase( midIdx, posBoxes.end() );
      } else if ( posBoxes.empty() ) {
        typename vector<PtrBBox>::iterator midIdx;
        midIdx = negBoxes.begin() + negBoxes.size()/2;
        posBoxes.insert( posBoxes.end(), midIdx, negBoxes.end() );
        negBoxes.erase( midIdx, negBoxes.end() );
      }

      PtrAABB neg = make_shared<AABBtree<GeomPrimitive> >();
      PtrAABB pos = make_shared<AABBtree<GeomPrimitive> >();

      neg->build(negBoxes);
      if (!neg->empty()) children.push_back(neg);

      pos->build(posBoxes);
      if (!pos->empty()) children.push_back(pos);
    }

    void
    build( vector<PtrGeom> const & primitives ) {
      vector<PtrBBox> bboxes;
      bboxes.reserve(primitives.size());
      typename vector<PtrGeom>::const_iterator it;
      for ( it = primitives.begin(); it != primitives.end(); ++it )
        bboxes.push_back( make_shared<BBox<GeomPrimitive> >(*it) );
      this->build(bboxes);
    }

    void
    print( ostream_type & stream, int level = 0 ) const {
      if ( empty() ) {
        stream
          << "[EMPTY AABB tree]\n";
      } else {
        stream
          << "BBOX xmin = " << setw(12) << pBBox->xmin
          << " ymin = "     << setw(12) << pBBox->ymin
          << " xmax = "     << setw(12) << pBBox->xmax
          << " ymax = "     << setw(12) << pBBox->ymax
          << " level = "    << level   << "\n";
        typename vector<PtrAABB>::const_iterator it;
        for ( it = children.begin(); it != children.end(); ++it )
          (*it)->print( stream, level+1 );
      }
    }

    bool
    intersects( AABBtree<GeomPrimitive> const & tree ) const {

      // check bbox with
      if ( !tree.pBBox->intersect(*pBBox) ) return false;

      int icase = (children.empty() ? 0 : 1) +
                  (tree.children.empty()? 0 : 2);

      switch ( icase ) {
      case 0: // both leaf, use GeomPrimitive intersection algorithm
        return pBBox->pGeom->intersect( *(tree.pBBox->pGeom) );
      case 1: // first is a tree, second is a leaf
        { typename vector<PtrAABB>::const_iterator it;
          for ( it = children.begin(); it != children.end(); ++it )
            if ( tree.intersects(**it) )
              return true;
        }
        break;
      case 2: // first leaf, second is a tree
        { typename vector<PtrAABB>::const_iterator it;
          for ( it = tree.children.begin(); it != tree.children.end(); ++it )
            if ( this->intersects(**it) )
              return true;
        }
        break;
      case 3: // first is a tree, second is a tree
        { typename vector<PtrAABB>::const_iterator c1, c2;
          for ( c1 = children.begin(); c1 != children.end(); ++c1 )
            for ( c2 = tree.children.begin(); c2 != tree.children.end(); ++c2 )
              if ( (*c1)->intersects(**c2) )
                return true;
        }
        break;
      }
      return false;
    }

    void
    intersects( AABBtree<GeomPrimitive> const & tree,
                VecPairPtrBBox                & intersectionList,
                bool                            swap_pair = false ) const {

      // check bbox with
      if ( !tree.pBBox->intersect(*pBBox) ) return;

      int icase = (children.empty() ? 0 : 1) +
                  (tree.children.empty()? 0 : 2);

      switch ( icase ) {
      case 0: // both leaf
        if ( swap_pair )
          intersectionList.push_back( PairPtrBBox(tree.pBBox,pBBox) );
        else           intersectionList.push_back( PairPtrBBox(pBBox,tree.pBBox) );
        break;
      case 1: // first is a tree, second is a leaf
        { typename vector<PtrAABB>::const_iterator it;
          for ( it = children.begin(); it != children.end(); ++it)
            tree.intersects( **it, intersectionList, !swap_pair );
        }
        break;
      case 2: // first leaf, second is a tree
        { typename vector<PtrAABB>::const_iterator it;
          for ( it = tree.children.begin(); it != tree.children.end(); ++it)
            this->intersects( **it, intersectionList, swap_pair );
        }
        break;
      case 3: // first is a tree, second is a tree
        { typename vector<PtrAABB>::const_iterator c1, c2;
          for ( c1 = children.begin(); c1 != children.end(); ++c1 )
            for ( c2 = tree.children.begin(); c2 != tree.children.end(); ++c2 )
              (*c1)->intersects(**c2, intersectionList, swap_pair );
        }
        break;
      }
    }
  };

}

#endif

///
/// eof: AABBtree.hh
///
