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
// file: Utils_AABBtree.hh (header-only)
//

#pragma once

#ifndef UTILS_AABB_TREE_dot_HH
#define UTILS_AABB_TREE_dot_HH

#include <map>
#include <set>

#include "Utils.hh"
#include "Utils_fmt.hh"

namespace Utils
{
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
  template <typename Real> class AABBtree
  {
  public:
    using integer  = int;                                   //!< Integer type.
    using AABB_SET = std::set<integer>;                     //!< Set of integers representing bounding box indices.
    using AABB_MAP = std::map<integer, std::set<integer>>;  //!< Map for bounding box indices
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

    mutable std::vector<integer> m_stack;

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

    static bool overlap1( Real const bbox1[], Real const bbox2[], integer )
    {
      return bbox1[0] <= bbox2[1] && bbox1[1] >= bbox2[0];
    }

    static bool overlap2( Real const bbox1[], Real const bbox2[], integer )
    {
      return bbox1[0] <= bbox2[2] && bbox1[2] >= bbox2[0] && bbox1[1] <= bbox2[3] && bbox1[3] >= bbox2[1];
    }

    static bool overlap3( Real const bbox1[], Real const bbox2[], integer )
    {
      return bbox1[0] <= bbox2[3] && bbox1[3] >= bbox2[0] && bbox1[1] <= bbox2[4] && bbox1[4] >= bbox2[1] &&
             bbox1[2] <= bbox2[5] && bbox1[5] >= bbox2[2];
    }

    static bool overlap4( Real const bbox1[], Real const bbox2[], integer )
    {
      return bbox1[0] <= bbox2[4] && bbox1[4] >= bbox2[0] && bbox1[1] <= bbox2[5] && bbox1[5] >= bbox2[1] &&
             bbox1[2] <= bbox2[6] && bbox1[6] >= bbox2[2] && bbox1[3] <= bbox2[7] && bbox1[7] >= bbox2[3];
    }

    static bool overlap5( Real const bbox1[], Real const bbox2[], integer )
    {
      return bbox1[0] <= bbox2[5] && bbox1[5] >= bbox2[0] && bbox1[1] <= bbox2[6] && bbox1[6] >= bbox2[1] &&
             bbox1[2] <= bbox2[7] && bbox1[7] >= bbox2[2] && bbox1[3] <= bbox2[8] && bbox1[8] >= bbox2[3] &&
             bbox1[4] <= bbox2[9] && bbox1[9] >= bbox2[4];
    }

    static bool overlap6( Real const bbox1[], Real const bbox2[], integer )
    {
      return bbox1[0] <= bbox2[6] && bbox1[6] >= bbox2[0] && bbox1[1] <= bbox2[7] && bbox1[7] >= bbox2[1] &&
             bbox1[2] <= bbox2[8] && bbox1[8] >= bbox2[2] && bbox1[3] <= bbox2[9] && bbox1[9] >= bbox2[3] &&
             bbox1[4] <= bbox2[10] && bbox1[10] >= bbox2[4] && bbox1[5] <= bbox2[11] && bbox1[11] >= bbox2[5];
    }

    static bool overlap7( Real const bbox1[], Real const bbox2[], integer )
    {
      return bbox1[0] <= bbox2[7] && bbox1[7] >= bbox2[0] && bbox1[1] <= bbox2[8] && bbox1[8] >= bbox2[1] &&
             bbox1[2] <= bbox2[9] && bbox1[9] >= bbox2[2] && bbox1[3] <= bbox2[10] && bbox1[10] >= bbox2[3] &&
             bbox1[4] <= bbox2[11] && bbox1[11] >= bbox2[4] && bbox1[5] <= bbox2[12] && bbox1[12] >= bbox2[5] &&
             bbox1[6] <= bbox2[13] && bbox1[13] >= bbox2[6];
    }

    static bool overlap8( Real const bbox1[], Real const bbox2[], integer )
    {
      return bbox1[0] <= bbox2[8] && bbox1[8] >= bbox2[0] && bbox1[1] <= bbox2[9] && bbox1[9] >= bbox2[1] &&
             bbox1[2] <= bbox2[10] && bbox1[10] >= bbox2[2] && bbox1[3] <= bbox2[11] && bbox1[11] >= bbox2[3] &&
             bbox1[4] <= bbox2[12] && bbox1[12] >= bbox2[4] && bbox1[5] <= bbox2[13] && bbox1[13] >= bbox2[5] &&
             bbox1[6] <= bbox2[14] && bbox1[14] >= bbox2[6] && bbox1[7] <= bbox2[15] && bbox1[15] >= bbox2[7];
    }

    static bool pnt_overlap1( Real const pnt[], Real const bb2[], integer )
    {
      return pnt[0] <= bb2[1] && pnt[0] >= bb2[0];
    }

    static bool pnt_overlap2( Real const pnt[], Real const bb2[], integer )
    {
      return pnt[0] <= bb2[2] && pnt[0] >= bb2[0] && pnt[1] <= bb2[3] && pnt[1] >= bb2[1];
    }

    static bool pnt_overlap3( Real const pnt[], Real const bbox2[], integer )
    {
      return pnt[0] <= bbox2[3] && pnt[0] >= bbox2[0] && pnt[1] <= bbox2[4] && pnt[1] >= bbox2[1] &&
             pnt[2] <= bbox2[5] && pnt[2] >= bbox2[2];
    }

    static bool pnt_overlap4( Real const pnt[], Real const bbox2[], integer )
    {
      return pnt[0] <= bbox2[4] && pnt[0] >= bbox2[0] && pnt[1] <= bbox2[5] && pnt[1] >= bbox2[1] &&
             pnt[2] <= bbox2[6] && pnt[2] >= bbox2[2] && pnt[3] <= bbox2[7] && pnt[3] >= bbox2[3];
    }

    static bool pnt_overlap5( Real const pnt[], Real const bbox2[], integer )
    {
      return pnt[0] <= bbox2[5] && pnt[0] >= bbox2[0] && pnt[1] <= bbox2[6] && pnt[1] >= bbox2[1] &&
             pnt[2] <= bbox2[7] && pnt[2] >= bbox2[2] && pnt[3] <= bbox2[8] && pnt[3] >= bbox2[3] &&
             pnt[4] <= bbox2[9] && pnt[4] >= bbox2[4];
    }

    static bool pnt_overlap6( Real const pnt[], Real const bbox2[], integer )
    {
      return pnt[0] <= bbox2[6] && pnt[0] >= bbox2[0] && pnt[1] <= bbox2[7] && pnt[1] >= bbox2[1] &&
             pnt[2] <= bbox2[8] && pnt[2] >= bbox2[2] && pnt[3] <= bbox2[9] && pnt[3] >= bbox2[3] &&
             pnt[4] <= bbox2[10] && pnt[4] >= bbox2[4] && pnt[5] <= bbox2[11] && pnt[5] >= bbox2[5];
    }

    static bool pnt_overlap7( Real const pnt[], Real const bbox2[], integer )
    {
      return pnt[0] <= bbox2[7] && pnt[0] >= bbox2[0] && pnt[1] <= bbox2[8] && pnt[1] >= bbox2[1] &&
             pnt[2] <= bbox2[9] && pnt[2] >= bbox2[2] && pnt[3] <= bbox2[10] && pnt[3] >= bbox2[3] &&
             pnt[4] <= bbox2[11] && pnt[4] >= bbox2[4] && pnt[5] <= bbox2[12] && pnt[5] >= bbox2[5] &&
             pnt[6] <= bbox2[13] && pnt[6] >= bbox2[6];
    }

    static bool pnt_overlap8( Real const pnt[], Real const bbox2[], integer )
    {
      return pnt[0] <= bbox2[8] && pnt[0] >= bbox2[0] && pnt[1] <= bbox2[9] && pnt[1] >= bbox2[1] &&
             pnt[2] <= bbox2[10] && pnt[2] >= bbox2[2] && pnt[3] <= bbox2[11] && pnt[3] >= bbox2[3] &&
             pnt[4] <= bbox2[12] && pnt[4] >= bbox2[4] && pnt[5] <= bbox2[13] && pnt[5] >= bbox2[5] &&
             pnt[6] <= bbox2[14] && pnt[6] >= bbox2[6] && pnt[7] <= bbox2[15] && pnt[7] >= bbox2[7];
    }

    static bool check_overlap( Real const bb1[], Real const bb2[], integer dim )
    {
      bool    overlap{ false };
      integer k{ dim % 4 };
      switch ( k )
      {
        case 1: overlap = bb1[0] <= bb2[dim] && bb1[dim] >= bb2[0]; break;
        case 2:
          overlap = bb1[0] <= bb2[dim + 0] && bb1[dim + 0] >= bb2[0] && bb1[1] <= bb2[dim + 1] &&
                    bb1[dim + 1] >= bb2[1];
          break;
        case 3:
          overlap = bb1[0] <= bb2[dim + 0] && bb1[dim + 0] >= bb2[0] && bb1[1] <= bb2[dim + 1] &&
                    bb1[dim + 1] >= bb2[1] && bb1[2] <= bb2[dim + 2] && bb1[dim + 2] >= bb2[2];
          break;
      }
      bb1 += k;
      bb2 += k;
      while ( !overlap && k < dim )
      {
        overlap = bb1[0] <= bb2[dim + 0] && bb1[dim + 0] >= bb2[0] && bb1[1] <= bb2[dim + 1] &&
                  bb1[dim + 1] >= bb2[1] && bb1[2] <= bb2[dim + 2] && bb1[dim + 2] >= bb2[2] &&
                  bb1[3] <= bb2[dim + 3] && bb1[dim + 3] >= bb2[3];
        bb1 += 4;
        bb2 += 4;
      }
      return overlap;
    }

    static bool check_overlap_with_point( Real const pnt[], Real const bb2[], integer dim )
    {
      bool    overlap{ false };
      integer k{ dim % 4 };
      switch ( k )
      {
        case 1: overlap = pnt[0] <= bb2[dim] && pnt[0] >= bb2[0]; break;
        case 2:
          overlap = pnt[0] <= bb2[dim + 0] && pnt[0] >= bb2[0] && pnt[1] <= bb2[dim + 1] && pnt[1] >= bb2[1];
          break;
        case 3:
          overlap = pnt[0] <= bb2[dim + 0] && pnt[0] >= bb2[0] && pnt[1] <= bb2[dim + 1] && pnt[1] >= bb2[1] &&
                    pnt[2] <= bb2[dim + 2] && pnt[2] >= bb2[2];
          break;
      }
      pnt += k;
      bb2 += k;
      while ( !overlap && k < dim )
      {
        overlap = pnt[0] <= bb2[dim + 0] && pnt[0] >= bb2[0] && pnt[1] <= bb2[dim + 1] && pnt[1] >= bb2[1] &&
                  pnt[2] <= bb2[dim + 2] && pnt[2] >= bb2[2] && pnt[3] <= bb2[dim + 3] && pnt[3] >= bb2[3];
        pnt += 4;
        bb2 += 4;
      }
      return overlap;
    }

    Real max_bbox_distance( Real const bbox[], Real const pnt[] ) const
    {
      Real res = 0;
      for ( integer i = 0; i < m_dim; ++i )
      {
        Real r1 = pnt[i] - bbox[i];
        Real r2 = pnt[i] - bbox[i + m_dim];
        Real mx = std::max( r1 * r1, r2 * r2 );
        res += mx * mx;
      }
      return std::sqrt( res );
    }

  public:
    //! Default constructor.
    AABBtree() = default;

    //!
    //! \brief Copy constructor for AABBtree.
    //!
    //! \param t Another AABBtree object to copy from.
    //!
    AABBtree( AABBtree const & T )
    {
      allocate( T.m_num_objects, T.m_dim );

      // Copia tutti i dati
      std::copy_n( T.m_father, m_nmax, m_father );
      std::copy_n( T.m_child, m_nmax, m_child );
      std::copy_n( T.m_ptr_nodes, m_nmax, m_ptr_nodes );
      std::copy_n( T.m_num_nodes, m_nmax, m_num_nodes );
      std::copy_n( T.m_id_nodes, m_num_objects, m_id_nodes );
      std::copy_n( T.m_bbox_tree, m_nmax * m_2dim, m_bbox_tree );
      std::copy_n( T.m_bbox_objs, m_num_objects * m_2dim, m_bbox_objs );

      // Copia i parametri
      m_max_num_objects_per_node = T.m_max_num_objects_per_node;
      m_bbox_long_edge_ratio     = T.m_bbox_long_edge_ratio;
      m_bbox_overlap_tolerance   = T.m_bbox_overlap_tolerance;
      m_bbox_min_size_tolerance  = T.m_bbox_min_size_tolerance;  // Aggiunto

      // Copia i puntatori a funzione (l'allocate li ha già impostati, ma li ricopiamo per sicurezza)
      m_check_overlap            = T.m_check_overlap;             // Aggiunto
      m_check_overlap_with_point = T.m_check_overlap_with_point;  // Aggiunto

      // Copia lo stato dell'albero
      m_num_tree_nodes = T.m_num_tree_nodes;

      // m_stack non viene copiato perché è temporaneo per le query
      // m_num_check viene azzerato per le statistiche
      m_num_check = 0;
    }

    //!
    //! \brief Sets the maximum number of objects per node.
    //!
    //! \param n Maximum number of objects.
    //!
    void set_max_num_objects_per_node( integer n )
    {
      UTILS_ASSERT(
        n > 0 && n <= 4096,
        "AABBtree::set_max_num_objects_per_node( nobj = {} )\n"
        "nobj must be > 0 and <= 4096\n",
        n );
      m_max_num_objects_per_node = n;
    }

    //!
    //! \brief Sets the bounding box long edge ratio.
    //!
    //! \param ratio Long edge ratio.
    //!
    void set_bbox_long_edge_ratio( Real ratio )
    {
      UTILS_ASSERT(
        ratio > 0 && ratio < 1,
        "AABBtree::set_bbox_long_edge_ratio( ratio = {} )\n"
        "tol must be > 0 and < 1\n",
        ratio );
      m_bbox_long_edge_ratio = ratio;
    }

    //!
    //! \brief Sets the bounding box overlap tolerance.
    //!
    //! \param tol Overlap tolerance.
    //!
    void set_bbox_overlap_tolerance( Real tol )
    {
      UTILS_ASSERT(
        tol > 0 && tol < 1,
        "AABBtree::set_bbox_overlap_tolerance( tol = {} )\n"
        "tol must be > 0 and < 1\n",
        tol );
      m_bbox_overlap_tolerance = tol;
    }

    //!
    //! \brief Sets the bounding box minimum size tolerance.
    //!
    //! \param tol Minimum size tolerance.
    //!
    void set_bbox_min_size_tolerance( Real tol )
    {
      UTILS_ASSERT(
        tol >= 0,
        "AABBtree::set_bbox_min_size_tolerance( tol = {} )\n"
        "tol must be >= 0\n",
        tol );
      m_bbox_min_size_tolerance = tol;
    }

    //!
    //! \brief Allocates memory for the AABB tree.
    //!
    //! \param nbox Number of bounding boxes.
    //! \param dim Dimension of the bounding boxes.
    //!
    void allocate( integer const nbox, integer const dim )
    {
      UTILS_WARNING(
        dim <= 10,
        "AABBtree::allocate( nbox, dim={} )\n"
        "dim is greather that 10!!!",
        dim );

      switch ( dim )
      {
        case 1:
          m_check_overlap            = overlap1;
          m_check_overlap_with_point = pnt_overlap1;
          break;
        case 2:
          m_check_overlap            = overlap2;
          m_check_overlap_with_point = pnt_overlap2;
          break;
        case 3:
          m_check_overlap            = overlap3;
          m_check_overlap_with_point = pnt_overlap3;
          break;
        case 4:
          m_check_overlap            = overlap4;
          m_check_overlap_with_point = pnt_overlap4;
          break;
        case 5:
          m_check_overlap            = overlap5;
          m_check_overlap_with_point = pnt_overlap5;
          break;
        case 6:
          m_check_overlap            = overlap6;
          m_check_overlap_with_point = pnt_overlap6;
          break;
        case 7:
          m_check_overlap            = overlap7;
          m_check_overlap_with_point = pnt_overlap7;
          break;
        case 8:
          m_check_overlap            = overlap8;
          m_check_overlap_with_point = pnt_overlap8;
          break;
        default:
          m_check_overlap            = check_overlap;
          m_check_overlap_with_point = check_overlap_with_point;
          break;
      }

      m_rmem.free();
      m_imem.free();

      m_dim         = dim;
      m_2dim        = 2 * dim;
      m_num_objects = nbox;
      m_nmax        = 2 * m_num_objects;  // estimate max memory usage

      m_rmem.allocate( ( m_nmax + m_num_objects ) * m_2dim );
      m_imem.allocate( ( 4 * m_nmax + m_num_objects ) );

      m_bbox_tree = m_rmem( m_nmax * m_2dim );
      m_bbox_objs = m_rmem( m_num_objects * m_2dim );

      m_father    = m_imem( m_nmax );
      m_child     = m_imem( m_nmax );
      m_ptr_nodes = m_imem( m_nmax );
      m_num_nodes = m_imem( m_nmax );
      m_id_nodes  = m_imem( m_num_objects );

      // initialize id nodes, will be reordered during the tree build
      for ( integer i{ 0 }; i < m_num_objects; ++i ) m_id_nodes[i] = i;

      // setup root node
      m_father[0]      = -1;
      m_child[0]       = -1;
      m_ptr_nodes[0]   = 0;
      m_num_nodes[0]   = m_num_objects;
      m_num_tree_nodes = 1;
    }

    //!
    //! \brief Adds bounding boxes to the tree.
    //!
    //! \param bb_min Minimum corners of bounding boxes.
    //! \param ldim0 Leading dimension for the minimum corners.
    //! \param bb_max Maximum corners of bounding boxes.
    //! \param ldim1 Leading dimension for the maximum corners.
    //!
    void add_bboxes( Real const bbox_min[], integer ldim0, Real const bbox_max[], integer ldim1 )
    {
      UTILS_ASSERT(
        ldim0 >= m_dim && ldim1 >= m_dim,
        "AABBtree::add_bboxes( bb_min, ldim0={}, bb_max, ldim1={} )\n"
        "must be ldim0, ldim1 >= dim = {}\n",
        ldim0,
        ldim1,
        m_dim );

      Real * bb = m_bbox_objs;
      for ( integer i = 0; i < m_num_objects; ++i )
      {
        for ( integer j = 0; j < m_dim; ++j )
        {
          UTILS_ASSERT( bbox_min[j] <= bbox_max[j], "AABBtree::add_bboxes, bad bbox N.{} max < min", i );
        }
        std::copy_n( bbox_min, m_dim, bb );
        bb += m_dim;
        bbox_min += ldim0;
        std::copy_n( bbox_max, m_dim, bb );
        bb += m_dim;
        bbox_max += ldim1;
      }
    }

    //!
    //! \brief Replaces a bounding box at a specific position.
    //!
    //! \param bbox_min New minimum corner of the bounding box.
    //! \param bbox_max New maximum corner of the bounding box.
    //! \param ipos Index of the bounding box to replace.
    //!
    void replace_bbox( Real const bbox_min[], Real const bbox_max[], integer ipos )
    {
      UTILS_ASSERT(
        ipos >= 0 && ipos < m_num_objects,
        "AABBtree::replace_bbox( bb_min, bb_max, ipos = {})"
        " ipos must be in [0,{})\n",
        ipos,
        m_num_objects );
      Real * bb = m_bbox_objs + ipos * m_2dim;
      for ( integer j = 0; j < m_dim; ++j )
      {
        UTILS_ASSERT( bbox_min[j] <= bbox_max[j], "AABBtree::replace_bbox, bad bbox N.{} max < min", ipos );
      }
      std::copy_n( bbox_min, m_dim, bb );
      std::copy_n( bbox_max, m_dim, bb + m_dim );
    }

    //!
    //! \brief Builds the AABB tree.
    //!
    void build()
    {
      Real otol{ Real( std::pow( m_bbox_overlap_tolerance, m_dim ) ) };

      // root contains all rectangles, build its bbox
      for ( integer j = 0; j < m_dim; ++j )
      {
        Real &       minj = m_bbox_tree[j];
        Real &       maxj = m_bbox_tree[m_dim + j];
        Real const * pmin{ m_bbox_objs + j };
        Real const * pmax{ m_bbox_objs + j + m_dim };
        minj = *pmin;
        maxj = *pmax;
        UTILS_ASSERT0( maxj >= minj, "AABBtree::build, bad bbox N.0 max < min" );
        for ( integer i = 1; i < m_num_objects; ++i )
        {
          pmin += m_2dim;
          pmax += m_2dim;
          UTILS_ASSERT( *pmax >= *pmin, "AABBtree::build, bad bbox N.{} max < min ({} < {})\n", i, *pmax, *pmin );
          if ( minj > *pmin ) minj = *pmin;
          if ( maxj < *pmax ) maxj = *pmax;
        }
      }

      // main loop: divide nodes until all constraints satisfied
      m_stack.clear();
      m_stack.reserve( 2 * m_num_objects + 1 );
      m_stack.emplace_back( 0 );
      m_num_tree_nodes = 1;

      while ( !m_stack.empty() )
      {
        // pop node from stack
        integer id_father = m_stack.back();
        m_stack.pop_back();
        UTILS_ASSERT_DEBUG(
          id_father < m_nmax,
          "AABBtree::build, id_father = {} must be less than m_nmax ={}\n",
          id_father,
          m_nmax );

        // set no childer for the moment
        m_child[id_father] = -1;

        // get rectangles id in parent
        integer const num{ m_num_nodes[id_father] };

        // if few bbox stop splitting
        if ( num < m_max_num_objects_per_node ) continue;

        integer const iptr{ m_ptr_nodes[id_father] };
        integer *     ptr = m_id_nodes + iptr;

        // split plane on longest axis, use euristic
        Real const * father_min = m_bbox_tree + id_father * m_2dim;
        Real const * father_max = father_min + m_dim;

        integer idim{ 0 };
        Real    mx{ father_max[0] - father_min[0] };
        for ( integer i{ 1 }; i < m_dim; ++i )
        {
          Real mx1 = father_max[i] - father_min[i];
          if ( mx < mx1 )
          {
            mx   = mx1;
            idim = i;
          }
        }

        // if too small bbox stop splitting
        if ( mx < m_bbox_min_size_tolerance ) continue;

        Real tol_len{ m_bbox_long_edge_ratio * mx };
        Real sp{ 0 };

        // separate short/long and accumulate short baricenter
        integer n_long{ 0 };
        integer n_short{ 0 };
        while ( n_long + n_short < num )
        {
          integer id = ptr[n_long];
          UTILS_ASSERT_DEBUG(
            id < m_num_objects,
            "AABBtree::build, id = {} must be less than m_num_objects ={}\n",
            id,
            m_num_objects );
          Real const * id_min = m_bbox_objs + id * m_2dim;
          Real const * id_max = id_min + m_dim;
          Real         id_len = id_max[idim] - id_min[idim];
          if ( id_len > tol_len )
          {
            // found long BBOX, increment n_long and update position
            ++n_long;
          }
          else
          {
            // found short BBOX, increment n_short and exchange with bottom
            ++n_short;
            std::swap( ptr[n_long], ptr[num - n_short] );
            sp += id_max[idim] + id_min[idim];
          }
        }

        // if split rectangles do not improve search, stop split at this level
        if ( n_short < 2 ) continue;

        // select the split position: take the mean of the set of
        // (non-"long") rectangle centers along axis idim
        sp /= 2 * n_short;

        // partition based on centers
        integer n_left{ 0 };
        integer n_right{ 0 };

        while ( n_long + n_left + n_right < num )
        {
          integer const id{ ptr[n_long + n_left] };
          Real const *  id_min{ m_bbox_objs + id * m_2dim };
          Real const *  id_max{ id_min + m_dim };
          Real          id_mid{ ( id_max[idim] + id_min[idim] ) / 2 };
          if ( id_mid < sp )
          {
            ++n_left;  // in right position do nothing
          }
          else
          {
            ++n_right;
            std::swap( ptr[n_long + n_left], ptr[num - n_right] );
          }
        }

        // if cannot improve bbox, stop split at this level!
        if ( n_left == 0 || n_right == 0 ) continue;

        // child indexing
        integer id_left{ m_num_tree_nodes + 0 };
        integer id_right{ m_num_tree_nodes + 1 };

        UTILS_ASSERT_DEBUG(
          id_right < m_nmax,
          "AABBtree::build, id_right = {} must be less than m_nmax ={}\n",
          id_right,
          m_nmax );

        // compute bbox of left and right child
        Real * bb_left_min{ m_bbox_tree + id_left * m_2dim };
        Real * bb_left_max{ bb_left_min + m_dim };
        for ( integer i{ 0 }; i < n_left; ++i )
        {
          integer id{ ptr[n_long + i] };
          UTILS_ASSERT_DEBUG(
            id < m_num_objects,
            "AABBtree::build, id = {} must be less than m_num_objects ={}\n",
            id,
            m_num_objects );
          Real const * bb_id_min{ m_bbox_objs + id * m_2dim };
          Real const * bb_id_max{ bb_id_min + m_dim };
          if ( i == 0 ) { std::copy_n( bb_id_min, m_2dim, bb_left_min ); }
          else
          {
            for ( integer j{ 0 }; j < m_dim; ++j )
            {
              if ( bb_left_min[j] > bb_id_min[j] ) bb_left_min[j] = bb_id_min[j];
              if ( bb_left_max[j] < bb_id_max[j] ) bb_left_max[j] = bb_id_max[j];
            }
          }
        }

        Real * bb_right_min = m_bbox_tree + id_right * m_2dim;
        Real * bb_right_max = bb_right_min + m_dim;
        for ( integer i = 0; i < n_right; ++i )
        {
          integer id = ptr[n_long + n_left + i];
          UTILS_ASSERT_DEBUG(
            id < m_num_objects,
            "AABBtree::build, id = {} must be less than m_num_objects ={}\n",
            id,
            m_num_objects );
          Real const * bb_id_min = m_bbox_objs + id * m_2dim;
          Real const * bb_id_max = bb_id_min + m_dim;
          if ( i == 0 ) { std::copy_n( bb_id_min, m_2dim, bb_right_min ); }
          else
          {
            for ( integer j = 0; j < m_dim; ++j )
            {
              if ( bb_right_min[j] > bb_id_min[j] ) bb_right_min[j] = bb_id_min[j];
              if ( bb_right_max[j] < bb_id_max[j] ) bb_right_max[j] = bb_id_max[j];
            }
          }
        }

        // check again if split improve the AABBtree otherwise stop exploration
        if ( n_left < m_max_num_objects_per_node || n_right < m_max_num_objects_per_node )
        {
          // few nodes, check if improve volume
          Real vo{ 1 };
          Real vL{ 1 };
          Real vR{ 1 };
          for ( integer j{ 0 }; j < m_dim; ++j )
          {
            Real Lmin{ bb_left_min[j] };
            Real Lmax{ bb_left_max[j] };
            Real Rmin{ bb_right_min[j] };
            Real Rmax{ bb_right_max[j] };
            vo *= std::max( std::min( Lmax, Rmax ) - std::max( Lmin, Rmin ), Real( 0 ) );
            vL *= Lmax - Lmin;
            vR *= Rmax - Rmin;
          }
          // if do not improve volume, stop split at this level!
          if ( vo > ( vL + vR - vo ) * otol ) continue;
        }

        // push child nodes onto stack
        m_father[id_left]  = id_father;
        m_father[id_right] = id_father;
        m_child[id_father] = id_left;

        m_num_nodes[id_father] = n_long;

        m_ptr_nodes[id_left] = iptr + n_long;
        m_num_nodes[id_left] = n_left;

        m_ptr_nodes[id_right] = iptr + n_long + n_left;
        m_num_nodes[id_right] = n_right;

        m_stack.emplace_back( id_left );
        m_stack.emplace_back( id_right );
        m_num_tree_nodes += 2;
      }
    }

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
    void build( Real const bb_min[], integer ldim0, Real const bb_max[], integer ldim1, integer nbox, integer dim )
    {
      allocate( nbox, dim );
      add_bboxes( bb_min, ldim0, bb_max, ldim1 );
      build();
    }

    //!
    //! \brief Intersects the tree with a point.
    //! Note: This version returns candidates (superset of exact results).
    //!
    //! \param pnt      The point to intersect with.
    //! \param bb_index Set to store indices of intersecting bounding boxes.
    //!
    void intersect_with_one_point( Real const pnt[], AABB_SET & bb_index ) const
    {
      m_num_check = 0;

      // quick return on empty inputs
      if ( m_num_tree_nodes == 0 ) return;

      // descend tree from root
      m_stack.clear();
      m_stack.reserve( 2 * m_num_tree_nodes + 1 );
      m_stack.emplace_back( 0 );
      while ( !m_stack.empty() )
      {
        // pop node from stack
        integer id_father = m_stack.back();
        m_stack.pop_back();

        // get BBOX
        Real const * bb_father = m_bbox_tree + id_father * m_2dim;

        ++m_num_check;

        // if do not overlap skip
        if ( !m_check_overlap_with_point( pnt, bb_father, m_dim ) ) continue;

        // get rectangles id in parent
        this->get_bbox_indexes_of_a_node( id_father, bb_index );

        if ( integer nn{ m_child[id_father] }; nn > 0 )
        {  // root == 0, children > 0
          // push on stack children
          m_stack.emplace_back( nn );
          m_stack.emplace_back( nn + 1 );
        }
      }
    }

    //!
    //! \brief Intersects the tree with a bounding box.
    //! Note: This version returns candidates (superset of exact results).
    //!
    //! \param bbox The bounding box to intersect with.
    //! \param bb_index Set to store indices of intersecting bounding boxes.
    //!
    void intersect_with_one_bbox( Real const bbox[], AABB_SET & bb_index ) const
    {
      m_num_check = 0;

      // quick return on empty inputs
      if ( m_num_tree_nodes == 0 ) return;

      // descend tree from root
      m_stack.clear();
      m_stack.reserve( 2 * m_num_tree_nodes + 1 );
      m_stack.emplace_back( 0 );
      while ( !m_stack.empty() )
      {
        // pop node from stack
        integer id_father = m_stack.back();
        m_stack.pop_back();

        // get BBOX
        Real const * bb_father = m_bbox_tree + id_father * m_2dim;

        ++m_num_check;

        // if do not overlap skip
        if ( !m_check_overlap( bb_father, bbox, m_dim ) ) continue;

        // get rectangles id in parent
        this->get_bbox_indexes_of_a_node( id_father, bb_index );

        if ( integer nn{ m_child[id_father] }; nn > 0 )
        {  // root == 0, children > 0
          // push on stack children
          m_stack.emplace_back( nn );
          m_stack.emplace_back( nn + 1 );
        }
      }
    }

    //!
    //! \brief Intersects the tree with another AABB tree.
    //! Note: This version returns candidates (superset of exact results).
    //!
    //! \param aabb The other AABB tree to intersect with.
    //! \param bb_index Map to store indices of intersecting bounding boxes.
    //!
    void intersect( AABBtree const & aabb, AABB_MAP & bb_index ) const
    {
      UTILS_ASSERT( this->m_bbox_tree != nullptr, "AABBtree::intersect: this tree not built" );
      UTILS_ASSERT( aabb.m_bbox_tree != nullptr, "AABBtree::intersect: aabb tree not built" );

      m_num_check = 0;

      // quick return on empty inputs
      if ( this->m_num_tree_nodes == 0 || aabb.m_num_tree_nodes == 0 ) return;

      // descend tree from root
      m_stack.clear();
      m_stack.reserve( m_num_tree_nodes + aabb.m_num_tree_nodes + 2 );
      m_stack.emplace_back( 0 );
      m_stack.emplace_back( 0 );
      while ( !m_stack.empty() )
      {
        // pop node from stack
        integer root2 = m_stack.back();
        m_stack.pop_back();
        integer sroot1 = m_stack.back();
        m_stack.pop_back();
        integer root1 = sroot1 >= 0 ? sroot1 : -1 - sroot1;

        // check for intersection
        Real const * bb_root1 = this->m_bbox_tree + root1 * m_2dim;
        Real const * bb_root2 = aabb.m_bbox_tree + root2 * m_2dim;

        ++m_num_check;

        // if do not overlap skip
        if ( !m_check_overlap( bb_root1, bb_root2, m_dim ) ) continue;

        // check if there are elements to check
        integer const nn1{ this->m_num_nodes[root1] };
        integer const nn2{ aabb.m_num_nodes[root2] };
        if ( nn1 > 0 && nn2 > 0 )
        {
          // For each object in node root1, add all objects in node root2 as candidates
          integer const * ptr1{ this->m_id_nodes + this->m_ptr_nodes[root1] };
          integer const * ptr2{ aabb.m_id_nodes + aabb.m_ptr_nodes[root2] };
          for ( integer ii{ 0 }; ii < nn1; ++ii )
          {
            integer    s1{ ptr1[ii] };
            AABB_SET & BB{ bb_index[s1] };
            for ( integer jj{ 0 }; jj < nn2; ++jj ) { BB.insert( ptr2[jj] ); }
          }
        }

        integer id_lr1{ sroot1 >= 0 ? m_child[root1] : -1 };
        integer id_lr2{ aabb.m_child[root2] };

        if ( id_lr1 >= 0 && id_lr2 >= 0 )
        {
          // Both trees have children: explore all 4 combinations
          m_stack.emplace_back( id_lr1 );
          m_stack.emplace_back( id_lr2 );

          m_stack.emplace_back( id_lr1 + 1 );
          m_stack.emplace_back( id_lr2 );

          m_stack.emplace_back( id_lr1 );
          m_stack.emplace_back( id_lr2 + 1 );

          m_stack.emplace_back( id_lr1 + 1 );
          m_stack.emplace_back( id_lr2 + 1 );

          if ( nn1 > 0 )
          {
            m_stack.emplace_back( -1 - root1 );
            m_stack.emplace_back( root2 );
          }
        }
        else if ( id_lr1 >= 0 )
        {
          m_stack.emplace_back( id_lr1 );
          m_stack.emplace_back( root2 );
          m_stack.emplace_back( id_lr1 + 1 );
          m_stack.emplace_back( root2 );
          if ( nn1 > 0 )
          {
            m_stack.emplace_back( -1 - root1 );
            m_stack.emplace_back( root2 );
          }
        }
        else if ( id_lr2 >= 0 )
        {
          m_stack.emplace_back( sroot1 );
          m_stack.emplace_back( id_lr2 );
          m_stack.emplace_back( sroot1 );
          m_stack.emplace_back( id_lr2 + 1 );
        }
      }
    }

    //!
    //! \brief Intersects the tree with a point and refines the search.
    //! Note: This version returns exact results.
    //!
    //! \param pnt The point to intersect with.
    //! \param bb_index Set to store indices of intersecting bounding boxes.
    //!
    void intersect_with_one_point_and_refine( Real const pnt[], AABB_SET & bb_index ) const
    {
      UTILS_ASSERT( pnt != nullptr, "AABBtree::intersect_with_one_point: pnt is null" );
      UTILS_ASSERT( m_bbox_tree != nullptr, "AABBtree::intersect_with_one_point: tree not built" );

      m_num_check = 0;

      // quick return on empty inputs
      if ( m_num_tree_nodes == 0 ) return;

      // descend tree from root
      m_stack.clear();
      m_stack.reserve( 2 * m_num_tree_nodes + 1 );
      m_stack.emplace_back( 0 );
      while ( !m_stack.empty() )
      {
        // pop node from stack
        integer id_father = m_stack.back();
        m_stack.pop_back();

        // get BBOX
        Real const * bb_father = m_bbox_tree + id_father * m_2dim;

        ++m_num_check;

        // if do not overlap skip
        if ( !m_check_overlap_with_point( pnt, bb_father, m_dim ) ) continue;

        // refine candidate
        integer const   num{ this->m_num_nodes[id_father] };
        integer const * ptr{ this->m_id_nodes + this->m_ptr_nodes[id_father] };
        for ( integer ii{ 0 }; ii < num; ++ii )
        {
          integer      s{ ptr[ii] };
          Real const * bb_s{ m_bbox_objs + s * m_2dim };
          ++m_num_check;
          if ( m_check_overlap_with_point( pnt, bb_s, m_dim ) ) bb_index.insert( s );
        }

        if ( integer nn{ m_child[id_father] }; nn > 0 )
        {  // root == 0, children > 0
          // push on stack children
          m_stack.emplace_back( nn );
          m_stack.emplace_back( nn + 1 );
        }
      }
    }

    //!
    //! \brief Intersects the tree with a bounding box and refines the search.
    //! Note: This version returns exact results.
    //!
    //! \param bbox The bounding box to intersect with.
    //! \param bb_index Set to store indices of intersecting bounding boxes.
    //!
    void intersect_with_one_bbox_and_refine( Real const bbox[], AABB_SET & bb_index ) const
    {
      UTILS_ASSERT( bbox != nullptr, "AABBtree::intersect_with_one_bbox: bbox is null" );
      UTILS_ASSERT( m_bbox_tree != nullptr, "AABBtree::intersect_with_one_bbox: tree not built" );

      m_num_check = 0;

      // quick return on empty inputs
      if ( m_num_tree_nodes == 0 ) return;

      // descend tree from root
      m_stack.clear();
      m_stack.reserve( 2 * m_num_tree_nodes + 1 );
      m_stack.emplace_back( 0 );
      while ( !m_stack.empty() )
      {
        // pop node from stack
        integer id_father = m_stack.back();
        m_stack.pop_back();

        // get BBOX
        Real const * bb_father = m_bbox_tree + id_father * m_2dim;

        ++m_num_check;

        // if do not overlap skip
        if ( !m_check_overlap( bb_father, bbox, m_dim ) ) continue;

        // refine candidate
        integer const   num{ this->m_num_nodes[id_father] };
        integer const * ptr{ this->m_id_nodes + this->m_ptr_nodes[id_father] };
        for ( integer ii{ 0 }; ii < num; ++ii )
        {
          integer      s{ ptr[ii] };
          Real const * bb_s{ m_bbox_objs + ptr[ii] * m_2dim };
          ++m_num_check;
          if ( m_check_overlap( bb_s, bbox, m_dim ) ) bb_index.insert( s );
        }

        if ( integer nn{ m_child[id_father] }; nn > 0 )
        {  // root == 0, children > 0
          // push on stack children
          m_stack.emplace_back( nn );
          m_stack.emplace_back( nn + 1 );
        }
      }
    }

    //!
    //! \brief Intersects the tree with another AABB tree and refines the
    //! search.
    //! Note: This version returns exact results.
    //!
    //! \param aabb The other AABB tree to intersect with.
    //! \param bb_index Map to store indices of intersecting bounding boxes.
    //!
    void intersect_and_refine( AABBtree const & aabb, AABB_MAP & bb_index ) const
    {
      UTILS_ASSERT( this->m_bbox_tree != nullptr, "AABBtree::intersect_and_refine: this tree not built" );
      UTILS_ASSERT( aabb.m_bbox_tree != nullptr, "AABBtree::intersect_and_refine: aabb tree not built" );

      m_num_check = 0;

      // quick return on empty inputs
      if ( this->m_num_tree_nodes == 0 || aabb.m_num_tree_nodes == 0 ) return;

      // descend tree from root
      m_stack.clear();
      m_stack.reserve( m_num_tree_nodes + aabb.m_num_tree_nodes + 2 );
      m_stack.emplace_back( 0 );
      m_stack.emplace_back( 0 );
      while ( !m_stack.empty() )
      {
        // pop node from stack
        integer root2 = m_stack.back();
        m_stack.pop_back();
        integer sroot1 = m_stack.back();
        m_stack.pop_back();
        integer root1 = sroot1 >= 0 ? sroot1 : -1 - sroot1;

        // check for intersection
        Real const * bb_root1 = this->m_bbox_tree + root1 * m_2dim;
        Real const * bb_root2 = aabb.m_bbox_tree + root2 * m_2dim;

        ++m_num_check;

        // if do not overlap skip
        if ( !m_check_overlap( bb_root1, bb_root2, m_dim ) ) continue;

        // check if there are elements to check
        integer const nn1{ this->m_num_nodes[root1] };
        integer       nn2{ aabb.m_num_nodes[root2] };
        if ( nn1 > 0 && nn2 > 0 )
        {
          // construct list of intersecting candidated
          integer const * ptr1{ this->m_id_nodes + this->m_ptr_nodes[root1] };
          integer const * ptr2{ aabb.m_id_nodes + aabb.m_ptr_nodes[root2] };
          for ( integer ii{ 0 }; ii < nn1; ++ii )
          {
            integer      s1{ ptr1[ii] };
            Real const * bb_s1{ m_bbox_objs + s1 * m_2dim };
            AABB_SET &   BB{ bb_index[s1] };
            for ( integer jj{ 0 }; jj < nn2; ++jj )
            {
              integer      s2{ ptr2[jj] };
              Real const * bb_s2{ aabb.m_bbox_objs + s2 * m_2dim };
              ++m_num_check;
              if ( m_check_overlap( bb_s1, bb_s2, m_dim ) ) BB.insert( s2 );
            }
          }
        }

        integer id_lr1{ sroot1 >= 0 ? m_child[root1] : -1 };
        integer id_lr2{ aabb.m_child[root2] };

        if ( id_lr1 >= 0 && id_lr2 >= 0 )
        {
          // Both trees have children: explore all 4 combinations
          m_stack.emplace_back( id_lr1 );
          m_stack.emplace_back( id_lr2 );

          m_stack.emplace_back( id_lr1 + 1 );
          m_stack.emplace_back( id_lr2 );

          m_stack.emplace_back( id_lr1 );
          m_stack.emplace_back( id_lr2 + 1 );

          m_stack.emplace_back( id_lr1 + 1 );
          m_stack.emplace_back( id_lr2 + 1 );

          if ( nn1 > 0 )
          {
            m_stack.emplace_back( -1 - root1 );
            m_stack.emplace_back( root2 );
          }
        }
        else if ( id_lr1 >= 0 )
        {
          m_stack.emplace_back( id_lr1 );
          m_stack.emplace_back( root2 );
          m_stack.emplace_back( id_lr1 + 1 );
          m_stack.emplace_back( root2 );
          if ( nn1 > 0 )
          {
            m_stack.emplace_back( -1 - root1 );
            m_stack.emplace_back( root2 );
          }
        }
        else if ( id_lr2 >= 0 )
        {
          m_stack.emplace_back( sroot1 );
          m_stack.emplace_back( id_lr2 );
          m_stack.emplace_back( sroot1 );
          m_stack.emplace_back( id_lr2 + 1 );
        }
      }
    }

    //!
    //! \brief Finds candidates for minimum distance to a point.
    //!
    //! \param pnt The point to check against.
    //! \param bb_index Set to store indices of candidate bounding boxes.
    //!
    void min_distance_candidates( Real const pnt[], AABB_SET & bb_index ) const
    {
      Real dst2_min, dst2_max;

      // quick return on empty inputs
      bb_index.clear();
      if ( this->m_num_tree_nodes == 0 ) return;

      Real min_max_distance2 = Utils::Inf<Real>();

      // descend tree from root
      m_stack.clear();
      m_stack.reserve( m_num_tree_nodes + 1 );
      m_stack.emplace_back( 0 );
      while ( !m_stack.empty() )
      {
        // pop node from stack
        integer const id_father{ m_stack.back() };
        m_stack.pop_back();

        // get BBOX
        Real const * father_bbox{ this->m_bbox_tree + id_father * m_2dim };
        this->pnt_bbox_minmax( pnt, father_bbox, dst2_min, dst2_max );

        if ( dst2_min <= min_max_distance2 )
        {
          if ( m_num_nodes[id_father] > 0 && dst2_max < min_max_distance2 ) min_max_distance2 = dst2_max;

          if ( integer nn{ m_child[id_father] }; nn > 0 )
          {  // root == 0, children > 0
            // push on stack childrens
            m_stack.emplace_back( nn );
            m_stack.emplace_back( nn + 1 );
          }
        }
      }

      // descend tree from root
      m_stack.clear();
      m_stack.emplace_back( 0 );
      while ( !m_stack.empty() )
      {
        // pop node from stack
        integer id_father{ m_stack.back() };
        m_stack.pop_back();
        Real const * father_bbox{ this->m_bbox_tree + id_father * m_2dim };
        this->pnt_bbox_minmax( pnt, father_bbox, dst2_min, dst2_max );
        if ( dst2_min <= min_max_distance2 )
        {
          this->get_bbox_indexes_of_a_node( id_father, bb_index );
          if ( integer nn{ m_child[id_father] }; nn > 0 )
          {  // root == 0, children > 0
            // push on stack childrens
            m_stack.emplace_back( nn );
            m_stack.emplace_back( nn + 1 );
          }
        }
      }
    }

    //!
    //! \brief Calculates minimum and maximum distance from a point to a
    //! bounding box.
    //!
    //! \param pnt The point.
    //! \param bbox The bounding box.
    //! \param dmin Minimum distance.
    //! \param dmax Maximum distance.
    //!
    void pnt_bbox_minmax( Real const pnt[], Real const bbox[], Real & dmin, Real & dmax ) const
    {
      Real const * bb_max = bbox + m_dim;
      Real const * bb_min = bbox;
      dmin                = 0;
      dmax                = 0;
      for ( integer i = 0; i < m_dim; ++i )
      {
        // check overlap
        Real pi    = pnt[i];
        Real dpmin = 0;
        Real dpmax = 0;
        Real t1    = pi - bb_max[i];
        Real t2    = bb_min[i] - pi;
        if ( t1 > 0 )
          dpmin = t1;
        else
          dpmax = t1;
        if ( t2 > 0 )
          dpmin = t2;
        else
          dpmax = t2;
        dmin += dpmin * dpmin;
        dmax += dpmax * dpmax;
      }
    }

    //!
    //! \brief Returns the spatial dimension of the bounding boxes.
    //!
    //! \return Dimension of the bounding boxes.
    //!
    integer dim() const { return m_dim; }

    //!
    //! \brief Returns the number of objects (bounding boxes).
    //!
    //! \return Number of objects.
    //!
    integer num_objects() const { return m_num_objects; }

    //!
    //! \brief Returns the number of tree nodes.
    //!
    //! \return Number of tree nodes.
    //!
    integer num_tree_nodes() const { return m_num_tree_nodes; }

    //!
    //! \brief Returns the number of overlap checks performed.
    //!
    //! \return Number of checks.
    //!
    integer num_check() const { return m_num_check; }

    //!
    //! \brief Returns the number of tree nodes with at least nmin objects.
    //!
    //! \param nmin Minimum number of objects in the node.
    //! \return Number of nodes.
    //!
    integer num_tree_nodes( integer const nmin ) const
    {
      integer n{ 0 };
      for ( integer i{ 0 }; i < m_num_tree_nodes; ++i )
        if ( m_num_nodes[i] >= nmin ) ++n;
      return n;
    }

    //!
    //! \brief Gets the bounding box of the root node.
    //!
    //! \param bb_min Array to store the minimum corner.
    //! \param bb_max Array to store the maximum corner.
    //!
    void get_root_bbox( Real bb_min[], Real bb_max[] ) const
    {
      std::copy_n( m_bbox_tree, m_dim, bb_min );
      std::copy_n( m_bbox_tree + m_dim, m_dim, bb_max );
    }

    //!
    //! \brief Gets bounding boxes of the tree.
    //!
    //! \param bb_min Array to store the minimum corners.
    //! \param ldim0 Leading dimension for minimum corners.
    //! \param bb_max Array to store the maximum corners.
    //! \param ldim1 Leading dimension for maximum corners.
    //! \param nmin Minimum number of objects in a node.
    //!
    void get_bboxes_of_the_tree( Real bbox_min[], integer ldim0, Real bbox_max[], integer ldim1, integer nmin ) const
    {
      UTILS_ASSERT(
        ldim0 >= m_dim && ldim1 >= m_dim,
        "AABBtree::get_bboxes_of_the_tree(\n"
        "  bbox_min, ldim0={},\n"
        "  bbox_max, ldim1={},\n"
        "  nmin={} )\n"
        "must be nmin >= 0 and ldim0:1 >= {}\n",
        ldim0,
        ldim1,
        nmin,
        m_dim );

      for ( integer i = 0; i < m_num_tree_nodes; ++i )
      {
        if ( m_num_nodes[i] >= nmin )
        {
          Real const * b_min = m_bbox_tree + i * m_2dim;
          Real const * b_max = b_min + m_dim;
          std::copy_n( b_min, m_dim, bbox_min );
          bbox_min += ldim0;
          std::copy_n( b_max, m_dim, bbox_max );
          bbox_max += ldim1;
        }
      }
    }

    //!
    //! \brief Gets the indices of bounding boxes in a specific tree node.
    //!
    //! \param i_pos Index of the node.
    //! \param bb_index Set to store indices of bounding boxes.
    //!
    void get_bbox_indexes_of_a_node( integer i_pos, AABB_SET & bb_index ) const
    {
      UTILS_ASSERT(
        i_pos >= 0 && i_pos < m_num_tree_nodes,
        "AABBtree::get_bbox_indexes_of_a_node( i_pos={}, bb_index ) "
        "i_pos must be >= 0 and < {}\n",
        i_pos,
        m_num_tree_nodes );
      integer num = m_num_nodes[i_pos];
      integer ptr = m_ptr_nodes[i_pos];
      while ( num-- > 0 ) bb_index.insert( m_id_nodes[ptr++] );
    }

    //!
    //! \brief Returns information about the AABB tree.
    //!
    //! \return A string containing information about the tree.
    //!
    std::string info( string_view indent = "" ) const
    {
      integer nleaf = 0;
      integer nlong = 0;
      for ( integer i = 0; i < m_num_tree_nodes; ++i )
      {
        if ( m_child[i] < 0 )
          ++nleaf;
        else if ( m_num_nodes[i] > 0 )
          ++nlong;
      }
      std::string res;
      res += fmt::format( "{}┌───────── AABB tree info ────────┐\n", indent );
      res += fmt::format( "{}│  Dimension                {:<5} │\n", indent, m_dim );
      res += fmt::format( "{}│  Number of nodes          {:<5} │\n", indent, m_num_tree_nodes );
      res += fmt::format( "{}│  Number of leaf           {:<5} │\n", indent, nleaf );
      res += fmt::format( "{}│  Number of long node      {:<5} │\n", indent, nlong );
      res += fmt::format( "{}│  Number of objects        {:<5} │\n", indent, m_num_objects );
      res += fmt::format( "{}│  max_num_objects_per_node {:<5} │\n", indent, m_max_num_objects_per_node );
      res += fmt::format( "{}│  bbox_long_edge_ratio     {:<5} │\n", indent, m_bbox_long_edge_ratio );
      res += fmt::format( "{}│  bbox_overlap_tolerance   {:<5} │\n", indent, m_bbox_overlap_tolerance );
      res += fmt::format( "{}└─────────────────────────────────┘\n", indent );
      return res;
    }
  };

}  // namespace Utils

#endif

//
// eof: Utils_AABBtree.hh
//
