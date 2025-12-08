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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi\unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

//
// file: Utils_AABBtree.hh
//


#pragma once

#ifndef UTILS_AABB_TREE_dot_HH
#define UTILS_AABB_TREE_dot_HH

#include <map>
#include <set>
#include <string>
#include <vector>

#include "Utils.hh"

namespace Utils
{

  using std::map;
  using std::set;
  using std::string;
  using std::vector;

  /*\
   |      _        _    ____  ____  _
   |     / \      / \  | __ )| __ )| |_ _ __ ___  ___
   |    / _ \    / _ \ |  _ \|  _ \| __| '__/ _ \/ _ \
   |   / ___ \  / ___ \| |_) | |_) | |_| | |  __/  __/
   |  /_/   \_\/_/   \_\____/|____/ \__|_|  \___|\___|
  \*/

  //!
  //! \brief A class representing an axis-aligned bounding box tree.
  //!
  //! The AABBtree class provides an efficient way to store and query a set
  //! of axis-aligned bounding boxes. It supports various operations such
  //! as building the tree, adding bounding boxes, and performing intersection
  //! tests.
  //!
  //! \tparam Real Type of the real numbers (e.g., float, double).
  //!
  template <typename Real>
  class AABBtree
  {
  public:
    using integer  = int;                     //!< Integer type.
    using AABB_SET = set<integer>;            //!< Set of integers representing bounding box indices.
    using AABB_MAP = map<integer, AABB_SET>;  //!< Map for bounding box indices
                                              //!< and their overlaps.

  private:
    Malloc<Real>    m_rmem{ "AABBtree_real" };
    Malloc<integer> m_imem{ "AABBtree_integer" };

    // AABBtree structure
    integer m_dim{ 0 };
    integer m_2dim{ 0 };
    integer m_num_objects{ 0 };
    integer m_num_tree_nodes{ 0 };

    integer * m_father{ nullptr };     // m_nmax
    integer * m_child{ nullptr };      // m_nmax
    integer * m_ptr_nodes{ nullptr };  // m_nmax
    integer * m_num_nodes{ nullptr };  // m_nmax
    integer * m_id_nodes{ nullptr };   // m_num_objects
    Real *    m_bbox_tree{ nullptr };  // m_nmax*m_2dim
    Real *    m_bbox_objs{ nullptr };  // m_num_objects*m_2dim

    mutable vector<integer> m_stack;

    integer m_nmax{ 0 };

    // parameters
    integer m_max_num_objects_per_node{ 16 };
    Real    m_bbox_long_edge_ratio{ Real( 0.8 ) };
    Real    m_bbox_overlap_tolerance{ Real( 0.1 ) };
    Real    m_bbox_min_size_tolerance{ Real( 0 ) };

    // statistic
    mutable integer m_num_check = 0;

    using OVERLAP_FUN = bool ( * )( Real const bbox1[], Real const bbox2[], integer dim );

    OVERLAP_FUN m_check_overlap{ nullptr };
    OVERLAP_FUN m_check_overlap_with_point{ nullptr };

    static bool overlap1( Real const bbox1[], Real const bbox2[], integer dim );
    static bool overlap2( Real const bbox1[], Real const bbox2[], integer dim );
    static bool overlap3( Real const bbox1[], Real const bbox2[], integer dim );
    static bool overlap4( Real const bbox1[], Real const bbox2[], integer dim );
    static bool overlap5( Real const bbox1[], Real const bbox2[], integer dim );
    static bool overlap6( Real const bbox1[], Real const bbox2[], integer dim );
    static bool overlap7( Real const bbox1[], Real const bbox2[], integer dim );
    static bool overlap8( Real const bbox1[], Real const bbox2[], integer dim );

    static bool pnt_overlap1( Real const pnt[], Real const bbox2[], integer dim );
    static bool pnt_overlap2( Real const pnt[], Real const bbox2[], integer dim );
    static bool pnt_overlap3( Real const pnt[], Real const bbox2[], integer dim );
    static bool pnt_overlap4( Real const pnt[], Real const bbox2[], integer dim );
    static bool pnt_overlap5( Real const pnt[], Real const bbox2[], integer dim );
    static bool pnt_overlap6( Real const pnt[], Real const bbox2[], integer dim );
    static bool pnt_overlap7( Real const pnt[], Real const bbox2[], integer dim );
    static bool pnt_overlap8( Real const pnt[], Real const bbox2[], integer dim );

    static bool check_overlap( Real const bb1[], Real const bb2[], integer dim );
    static bool check_overlap_with_point( Real const bb1[], Real const pnt[], integer dim );

    Real max_bbox_distance( Real const bbox[], Real const pnt[] ) const;

  public:
    //! Default constructor.
    AABBtree() = default;

    //!
    //! \brief Copy constructor for AABBtree.
    //!
    //! \param t Another AABBtree object to copy from.
    //!
    AABBtree( AABBtree<Real> const & t );

    //!
    //! \brief Sets the maximum number of objects per node.
    //!
    //! \param n Maximum number of objects.
    //!
    void set_max_num_objects_per_node( integer n );

    //!
    //! \brief Sets the bounding box long edge ratio.
    //!
    //! \param ratio Long edge ratio.
    //!
    void set_bbox_long_edge_ratio( Real ratio );

    //!
    //! \brief Sets the bounding box overlap tolerance.
    //!
    //! \param tol Overlap tolerance.
    //!
    void set_bbox_overlap_tolerance( Real tol );

    //!
    //! \brief Sets the bounding box minimum size tolerance.
    //!
    //! \param tol Minimum size tolerance.
    //!
    void set_bbox_min_size_tolerance( Real tol );

    //!
    //! \brief Allocates memory for the AABB tree.
    //!
    //! \param nbox Number of bounding boxes.
    //! \param dim Dimension of the bounding boxes.
    //!
    void allocate( integer nbox, integer dim );

    //!
    //! \brief Adds bounding boxes to the tree.
    //!
    //! \param bb_min Minimum corners of bounding boxes.
    //! \param ldim0 Leading dimension for the minimum corners.
    //! \param bb_max Maximum corners of bounding boxes.
    //! \param ldim1 Leading dimension for the maximum corners.
    //!
    void add_bboxes( Real const bb_min[], integer ldim0, Real const bb_max[], integer ldim1 );

    //!
    //! \brief Replaces a bounding box at a specific position.
    //!
    //! \param bbox_min New minimum corner of the bounding box.
    //! \param bbox_max New maximum corner of the bounding box.
    //! \param ipos Index of the bounding box to replace.
    //!
    void replace_bbox( Real const bbox_min[], Real const bbox_max[], integer ipos );

    //!
    //! \brief Builds the AABB tree.
    //!
    void build();

    //!
    //! \brief Builds the AABB tree with specified bounding boxes.
    //!
    //! \param bb_min Minimum corners of bounding boxes.
    //! \param ldim0  Leading dimension for the minimum corners.
    //! \param bb_max Maximum corners of bounding boxes.
    //! \param ldim1  Leading dimension for the maximum corners.
    //! \param nbox   Number of bounding boxes.
    //! \param dim    Dimension of the bounding boxes.
    //!
    void
    build( Real const bb_min[], integer ldim0, Real const bb_max[], integer ldim1, integer nbox, integer dim )
    {
      allocate( nbox, dim );
      add_bboxes( bb_min, ldim0, bb_max, ldim1 );
      build();
    }

    //!
    //! \brief Intersects the tree with a point.
    //!
    //! \param pnt      The point to intersect with.
    //! \param bb_index Set to store indices of intersecting bounding boxes.
    //!
    void intersect_with_one_point( Real const pnt[], AABB_SET & bb_index ) const;

    //!
    //! \brief Intersects the tree with a bounding box.
    //!
    //! \param bbox The bounding box to intersect with.
    //! \param bb_index Set to store indices of intersecting bounding boxes.
    //!
    void intersect_with_one_bbox( Real const bbox[], AABB_SET & bb_index ) const;


    //!
    //! \brief Intersects the tree with another AABB tree.
    //!
    //! \param aabb The other AABB tree to intersect with.
    //! \param bb_index Map to store indices of intersecting bounding boxes.
    //!
    void intersect( AABBtree<Real> const & aabb, AABB_MAP & bb_index ) const;

    //!
    //! \brief Intersects the tree with a point and refines the search.
    //!
    //! \param pnt The point to intersect with.
    //! \param bb_index Set to store indices of intersecting bounding boxes.
    //!
    void intersect_with_one_point_and_refine( Real const pnt[], AABB_SET & bb_index ) const;

    //!
    //! \brief Intersects the tree with a bounding box and refines the search.
    //!
    //! \param bbox The bounding box to intersect with.
    //! \param bb_index Set to store indices of intersecting bounding boxes.
    //!
    void intersect_with_one_bbox_and_refine( Real const bbox[], AABB_SET & bb_index ) const;


    //!
    //! \brief Intersects the tree with another AABB tree and refines the
    //! search.
    //!
    //! \param aabb The other AABB tree to intersect with.
    //! \param bb_index Map to store indices of intersecting bounding boxes.
    //!
    void intersect_and_refine( AABBtree<Real> const & aabb, AABB_MAP & bb_index ) const;

    //!
    //! \brief Finds candidates for minimum distance to a point.
    //!
    //! \param pnt The point to check against.
    //! \param bb_index Set to store indices of candidate bounding boxes.
    //!
    void min_distance_candidates( Real const pnt[], AABB_SET & bb_index ) const;

    //!
    //! \brief Calculates minimum and maximum distance from a point to a
    //! bounding box.
    //!
    //! \param pnt The point.
    //! \param bbox The bounding box.
    //! \param dmin Minimum distance.
    //! \param dmax Maximum distance.
    //!
    void pnt_bbox_minmax( Real const pnt[], Real const bbox[], Real & dmin, Real & dmax ) const;

    //!
    //! \brief Returns the spatial dimension of the bounding boxes.
    //!
    //! \return Dimension of the bounding boxes.
    //!
    integer
    dim() const
    {
      return m_dim;
    }

    //!
    //! \brief Returns the number of objects (bounding boxes).
    //!
    //! \return Number of objects.
    //!
    integer
    num_objects() const
    {
      return m_num_objects;
    }

    //!
    //! \brief Returns the number of tree nodes.
    //!
    //! \return Number of tree nodes.
    //!
    integer
    num_tree_nodes() const
    {
      return m_num_tree_nodes;
    }

    //!
    //! \brief Returns the number of overlap checks performed.
    //!
    //! \return Number of checks.
    //!
    integer
    num_check() const
    {
      return m_num_check;
    }

    //!
    //! \brief Returns the number of tree nodes with at least nmin objects.
    //!
    //! \param nmin Minimum number of objects in the node.
    //! \return Number of nodes.
    //!
    integer num_tree_nodes( integer nmin ) const;

    //!
    //! \brief Gets the bounding box of the root node.
    //!
    //! \param bb_min Array to store the minimum corner.
    //! \param bb_max Array to store the maximum corner.
    //!
    void get_root_bbox( Real bb_min[], Real bb_max[] ) const;

    //!
    //! \brief Gets bounding boxes of the tree.
    //!
    //! \param bb_min Array to store the minimum corners.
    //! \param ldim0 Leading dimension for minimum corners.
    //! \param bb_max Array to store the maximum corners.
    //! \param ldim1 Leading dimension for maximum corners.
    //! \param nmin Minimum number of objects in a node.
    //!
    void get_bboxes_of_the_tree( Real bb_min[], integer ldim0, Real bb_max[], integer ldim1, integer nmin ) const;

    //!
    //! \brief Gets the indices of bounding boxes in a specific tree node.
    //!
    //! \param i_pos Index of the node.
    //! \param bb_index Set to store indices of bounding boxes.
    //!
    void get_bbox_indexes_of_a_node( integer i_pos, AABB_SET & bb_index ) const;

    //!
    //! \brief Returns information about the AABB tree.
    //!
    //! \return A string containing information about the tree.
    //!
    string info() const;
  };

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */

#ifndef UTILS_OS_WINDOWS
  extern template class AABBtree<float>;
  extern template class AABBtree<double>;
#endif

}  // namespace Utils

#endif

//
// eof: Utils_AABBtree.hh
//
