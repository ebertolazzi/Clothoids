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
/// file: Utils_AABBtree.hh
///
#pragma once

#ifndef UTILS_AABB_TREE_dot_HH
#define UTILS_AABB_TREE_dot_HH

#include "Utils.hh"

#include <string>
#include <vector>
#include <set>
#include <map>

namespace Utils {

  using std::string;
  using std::vector;
  using std::set;
  using std::map;

  /*\
   |      _        _    ____  ____  _
   |     / \      / \  | __ )| __ )| |_ _ __ ___  ___
   |    / _ \    / _ \ |  _ \|  _ \| __| '__/ _ \/ _ \
   |   / ___ \  / ___ \| |_) | |_) | |_| | |  __/  __/
   |  /_/   \_\/_/   \_\____/|____/ \__|_|  \___|\___|
  \*/

  template <typename Real>
  class AABBtree {
  public:

    typedef int integer;
    typedef set<integer>     SET;
    typedef map<integer,SET> MAP;

  private:

    Malloc<Real>    m_rmem;
    Malloc<integer> m_imem;

    // AABBtree structure
    integer m_dim            = 0;
    integer m_2dim           = 0;
    integer m_num_objects    = 0;
    integer m_num_tree_nodes = 0;

    integer * m_father    = nullptr;
    integer * m_child     = nullptr;
    integer * m_ptr_nodes = nullptr;
    integer * m_num_nodes = nullptr;
    integer * m_id_nodes  = nullptr;
    integer * m_stack     = nullptr;
    Real    * m_bbox_tree = nullptr;
    Real    * m_bbox_objs = nullptr;

    // parameters
    integer m_max_num_objects_per_node = 16;
    Real    m_bbox_long_edge_ratio     = Real(0.8);
    Real    m_bbox_overlap_tolerance   = Real(0.1);

    // statistic
    mutable integer m_num_check = 0;

    Real max_bbox_distance( Real const * bbox, Real const * pnt ) const;

  public:

    AABBtree() : m_rmem("AABBtree"), m_imem("AABBtree") {}

    AABBtree( AABBtree<Real> const & t );

    void set_max_num_objects_per_node( integer n );
    void set_bbox_long_edge_ratio( Real ratio );
    void set_bbox_overlap_tolerance( Real tol );

    void allocate( integer nbox, integer dim );

    void
    add_bboxes(
      Real const * bb_min, integer ldim0,
      Real const * bb_max, integer ldim1
    );

    void build();

    void
    build(
      Real const * bb_min, integer ldim0,
      Real const * bb_max, integer ldim1,
      integer nbox,
      integer dim
    ) {
      allocate( nbox, dim );
      add_bboxes( bb_min, ldim0, bb_max, ldim1 );
      build();
    }

    void intersect_with_one_point( Real const * pnt, SET & bb_index ) const;
    void intersect_with_one_bbox( Real const * bbox, SET & bb_index ) const;
    void intersect( AABBtree<Real> const & aabb, MAP & bb_index ) const;

    void intersect_with_one_point_and_refine( Real const * pnt, SET & bb_index ) const;
    void intersect_with_one_bbox_and_refine( Real const * bbox, SET & bb_index ) const;
    void intersect_and_refine( AABBtree<Real> const & aabb, MAP & bb_index ) const;

    Real minimum_max_bbox_distance( Real const * pnt ) const;

    integer dim()            const { return m_dim; }
    integer num_objects()    const { return m_num_objects; }
    integer num_tree_nodes() const { return m_num_tree_nodes; }
    integer num_check()      const { return m_num_check; }

    integer num_tree_nodes( integer nmin ) const;
    void
    get_bboxes_of_the_tree(
      Real * bb_min, integer ldim0,
      Real * bb_max, integer ldim1,
      integer nmin
    ) const;

    void
    get_bbox_indexes_of_a_node( integer i_pos, SET & bb_index ) const;

    string info() const;
  };

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */

  #ifndef UTILS_OS_WINDOWS
  extern template class AABBtree<float>;
  extern template class AABBtree<double>;
  #endif

}

#endif

///
/// eof: Utils_AABBtree.hh
///
