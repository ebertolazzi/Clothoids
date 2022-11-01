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
/// file: Utils_AABB_tree.cc
///

#include "Utils_AABB_tree.hh"
#include <algorithm>
#include <utility>

namespace Utils {

  using std::max;
  using std::min;
  using std::swap;
  using std::copy_n;

  template <typename Real>
  bool
  AABBtree<Real>::overlap1( Real const bb1[], Real const bb2[], integer ) {
    return bb1[0] <= bb2[1] && bb1[1] >= bb2[0];
  }

  template <typename Real>
  bool
  AABBtree<Real>::overlap2( Real const bb1[], Real const bb2[], integer ) {
    return bb1[0] <= bb2[2] && bb1[2] >= bb2[0] &&
           bb1[1] <= bb2[3] && bb1[3] >= bb2[1];
  }

  template <typename Real>
  bool
  AABBtree<Real>::overlap3( Real const bb1[], Real const bb2[], integer ) {
    return bb1[0] <= bb2[3] && bb1[3] >= bb2[0] &&
           bb1[1] <= bb2[4] && bb1[4] >= bb2[1] &&
           bb1[2] <= bb2[5] && bb1[5] >= bb2[2];
  }

  template <typename Real>
  bool
  AABBtree<Real>::overlap4( Real const bb1[], Real const bb2[], integer ) {
    return bb1[0] <= bb2[4] && bb1[4] >= bb2[0] &&
           bb1[1] <= bb2[5] && bb1[5] >= bb2[1] &&
           bb1[2] <= bb2[6] && bb1[6] >= bb2[2] &&
           bb1[3] <= bb2[7] && bb1[7] >= bb2[3];
  }

  template <typename Real>
  bool
  AABBtree<Real>::overlap5( Real const bb1[], Real const bb2[], integer ) {
    return bb1[0] <= bb2[5] && bb1[5] >= bb2[0] &&
           bb1[1] <= bb2[6] && bb1[6] >= bb2[1] &&
           bb1[2] <= bb2[7] && bb1[7] >= bb2[2] &&
           bb1[3] <= bb2[8] && bb1[8] >= bb2[3] &&
           bb1[4] <= bb2[9] && bb1[9] >= bb2[4];
  }

  template <typename Real>
  bool
  AABBtree<Real>::overlap6( Real const bb1[], Real const bb2[], integer ) {
    return bb1[0] <= bb2[6]  && bb1[6]  >= bb2[0] &&
           bb1[1] <= bb2[7]  && bb1[7]  >= bb2[1] &&
           bb1[2] <= bb2[8]  && bb1[8]  >= bb2[2] &&
           bb1[3] <= bb2[9]  && bb1[9]  >= bb2[3] &&
           bb1[4] <= bb2[10] && bb1[10] >= bb2[4] &&
           bb1[5] <= bb2[11] && bb1[11] >= bb2[5];
  }

  template <typename Real>
  bool
  AABBtree<Real>::overlap7( Real const bb1[], Real const bb2[], integer ) {
    return bb1[0] <= bb2[7]  && bb1[7]  >= bb2[0] &&
           bb1[1] <= bb2[8]  && bb1[8]  >= bb2[1] &&
           bb1[2] <= bb2[9]  && bb1[9]  >= bb2[2] &&
           bb1[3] <= bb2[10] && bb1[10] >= bb2[3] &&
           bb1[4] <= bb2[11] && bb1[11] >= bb2[4] &&
           bb1[5] <= bb2[12] && bb1[12] >= bb2[5] &&
           bb1[6] <= bb2[13] && bb1[14] >= bb2[6];
  }

  template <typename Real>
  bool
  AABBtree<Real>::overlap8( Real const bb1[], Real const bb2[], integer ) {
    return bb1[0] <= bb2[8]  && bb1[8]  >= bb2[0] &&
           bb1[1] <= bb2[9]  && bb1[9]  >= bb2[1] &&
           bb1[2] <= bb2[10] && bb1[10] >= bb2[2] &&
           bb1[3] <= bb2[11] && bb1[11] >= bb2[3] &&
           bb1[4] <= bb2[12] && bb1[12] >= bb2[4] &&
           bb1[5] <= bb2[13] && bb1[13] >= bb2[5] &&
           bb1[6] <= bb2[14] && bb1[14] >= bb2[6] &&
           bb1[7] <= bb2[15] && bb1[15] >= bb2[7];
  }

  template <typename Real>
  bool
  AABBtree<Real>::pnt_overlap1( Real const pnt[], Real const bb2[], integer ) {
    return pnt[0] <= bb2[1] && pnt[0] >= bb2[0];
  }

  template <typename Real>
  bool
  AABBtree<Real>::pnt_overlap2( Real const pnt[], Real const bb2[], integer ) {
    return pnt[0] <= bb2[2] && pnt[0] >= bb2[0] &&
           pnt[1] <= bb2[3] && pnt[1] >= bb2[1];
  }

  template <typename Real>
  bool
  AABBtree<Real>::pnt_overlap3( Real const pnt[], Real const bb2[], integer ) {
    return pnt[0] <= bb2[3] && pnt[0] >= bb2[0] &&
           pnt[1] <= bb2[4] && pnt[1] >= bb2[1] &&
           pnt[2] <= bb2[5] && pnt[2] >= bb2[2];
  }

  template <typename Real>
  bool
  AABBtree<Real>::pnt_overlap4( Real const pnt[], Real const bb2[], integer ) {
    return pnt[0] <= bb2[4] && pnt[0] >= bb2[0] &&
           pnt[1] <= bb2[5] && pnt[1] >= bb2[1] &&
           pnt[2] <= bb2[6] && pnt[2] >= bb2[2] &&
           pnt[3] <= bb2[7] && pnt[3] >= bb2[3];
  }

  template <typename Real>
  bool
  AABBtree<Real>::pnt_overlap5( Real const pnt[], Real const bb2[], integer ) {
    return pnt[0] <= bb2[5] && pnt[0] >= bb2[0] &&
           pnt[1] <= bb2[6] && pnt[1] >= bb2[1] &&
           pnt[2] <= bb2[7] && pnt[2] >= bb2[2] &&
           pnt[3] <= bb2[8] && pnt[3] >= bb2[3] &&
           pnt[4] <= bb2[9] && pnt[4] >= bb2[4];
  }

  template <typename Real>
  bool
  AABBtree<Real>::pnt_overlap6( Real const pnt[], Real const bb2[], integer ) {
    return pnt[0] <= bb2[6]  && pnt[0] >= bb2[0] &&
           pnt[1] <= bb2[7]  && pnt[1] >= bb2[1] &&
           pnt[2] <= bb2[8]  && pnt[2] >= bb2[2] &&
           pnt[3] <= bb2[9]  && pnt[3] >= bb2[3] &&
           pnt[4] <= bb2[10] && pnt[4] >= bb2[4] &&
           pnt[5] <= bb2[11] && pnt[5] >= bb2[5];
  }

  template <typename Real>
  bool
  AABBtree<Real>::pnt_overlap7( Real const pnt[], Real const bb2[], integer ) {
    return pnt[0] <= bb2[7]  && pnt[0] >= bb2[0] &&
           pnt[1] <= bb2[8]  && pnt[1] >= bb2[1] &&
           pnt[2] <= bb2[9]  && pnt[2] >= bb2[2] &&
           pnt[3] <= bb2[10] && pnt[3] >= bb2[3] &&
           pnt[4] <= bb2[11] && pnt[4] >= bb2[4] &&
           pnt[5] <= bb2[12] && pnt[5] >= bb2[5] &&
           pnt[6] <= bb2[13] && pnt[6] >= bb2[6];
  }

  template <typename Real>
  bool
  AABBtree<Real>::pnt_overlap8( Real const pnt[], Real const bb2[], integer ) {
    return pnt[0] <= bb2[8]  && pnt[0] >= bb2[0] &&
           pnt[1] <= bb2[9]  && pnt[1] >= bb2[1] &&
           pnt[2] <= bb2[10] && pnt[2] >= bb2[2] &&
           pnt[3] <= bb2[11] && pnt[3] >= bb2[3] &&
           pnt[4] <= bb2[12] && pnt[4] >= bb2[4] &&
           pnt[5] <= bb2[13] && pnt[5] >= bb2[5] &&
           pnt[6] <= bb2[14] && pnt[6] >= bb2[6] &&
           pnt[7] <= bb2[15] && pnt[7] >= bb2[7];
  }

  template <typename Real>
  bool
  AABBtree<Real>::check_overlap( Real const bb1[], Real const bb2[], integer dim ) {
    bool overlap = false;
    integer k = dim % 4;
    switch ( k ) {
    case 1:
      overlap = bb1[0] <= bb2[dim] && bb1[dim] >= bb2[0];
      break;
    case 2:
      overlap = bb1[0] <= bb2[dim+0] && bb1[dim+0] >= bb2[0] &&
                bb1[1] <= bb2[dim+1] && bb1[dim+1] >= bb2[1];
      break;
    case 3:
      overlap = bb1[0] <= bb2[dim+0] && bb1[dim+0] >= bb2[0] &&
                bb1[1] <= bb2[dim+1] && bb1[dim+1] >= bb2[1] &&
                bb1[2] <= bb2[dim+2] && bb1[dim+2] >= bb2[2];
      break;
    }
    bb1 += k;
    bb2 += k;
    while ( !overlap && k < dim ) {
      overlap = bb1[0] <= bb2[dim+0] && bb1[dim+0] >= bb2[0] &&
                bb1[1] <= bb2[dim+1] && bb1[dim+1] >= bb2[1] &&
                bb1[2] <= bb2[dim+2] && bb1[dim+2] >= bb2[2] &&
                bb1[3] <= bb2[dim+3] && bb1[dim+3] >= bb2[3];
      bb1 += 4;
      bb2 += 4;
    }
    return overlap;
  }

  template <typename Real>
  bool
  AABBtree<Real>::check_overlap_with_point( Real const pnt[], Real const bb2[], integer dim ) {
    bool overlap = false;
    integer k = dim % 4;
    switch ( k ) {
    case 1:
      overlap = pnt[0] <= bb2[dim] && pnt[0] >= bb2[0];
      break;
    case 2:
      overlap = pnt[0] <= bb2[dim+0] && pnt[0] >= bb2[0] &&
                pnt[1] <= bb2[dim+1] && pnt[1] >= bb2[1];
      break;
    case 3:
      overlap = pnt[0] <= bb2[dim+0] && pnt[0] >= bb2[0] &&
                pnt[1] <= bb2[dim+1] && pnt[1] >= bb2[1] &&
                pnt[2] <= bb2[dim+2] && pnt[2] >= bb2[2];
      break;
    }
    pnt += k;
    bb2 += k;
    while ( !overlap && k < dim ) {
      overlap = pnt[0] <= bb2[dim+0] && pnt[0] >= bb2[0] &&
                pnt[1] <= bb2[dim+1] && pnt[1] >= bb2[1] &&
                pnt[2] <= bb2[dim+2] && pnt[2] >= bb2[2] &&
                pnt[3] <= bb2[dim+3] && pnt[3] >= bb2[3];
      pnt += 4;
      bb2 += 4;
    }
    return overlap;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  AABBtree<Real>::AABBtree( AABBtree<Real> const & T )
  {
    allocate( T.m_num_objects, T.m_dim );

    std::copy_n( T.m_father,    m_nmax,               m_father    );
    std::copy_n( T.m_child,     m_nmax,               m_child     );
    std::copy_n( T.m_ptr_nodes, m_nmax,               m_ptr_nodes );
    std::copy_n( T.m_num_nodes, m_nmax,               m_num_nodes );
    std::copy_n( T.m_id_nodes,  m_num_objects,        m_id_nodes  );
    std::copy_n( T.m_bbox_tree, m_nmax*m_2dim,        m_bbox_tree );
    std::copy_n( T.m_bbox_objs, m_num_objects*m_2dim, m_bbox_objs );

    m_max_num_objects_per_node = T.m_max_num_objects_per_node;
    m_bbox_long_edge_ratio     = T.m_bbox_long_edge_ratio;
    m_bbox_overlap_tolerance   = T.m_bbox_overlap_tolerance;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  Real
  AABBtree<Real>::max_bbox_distance( Real const bbox[], Real const pnt[] ) const {
    Real res = 0;
    for ( integer i = 0; i < m_dim; ++i ) {
      Real r1 = pnt[i] - bbox[i];
      Real r2 = pnt[i] - bbox[i+m_dim];
      Real mx = max(r1*r1,r2*r2);
      res += mx*mx;
    }
    return sqrt(res);
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  string
  AABBtree<Real>::info() const {
    integer nleaf = 0;
    integer nlong = 0;
    for ( integer i = 0; i < m_num_tree_nodes; ++i ) {
      if      ( m_child[i] < 0     ) ++nleaf;
      else if ( m_num_nodes[i] > 0 ) ++nlong;
    }
    string res = "-------- AABB tree info --------\n";
    res += fmt::format( "  Dimension                {}\n", m_dim );
    res += fmt::format( "  Number of nodes          {}\n", m_num_tree_nodes );
    res += fmt::format( "  Number of leaf           {}\n", nleaf );
    res += fmt::format( "  Number of long node      {}\n", nlong );
    res += fmt::format( "  Number of objects        {}\n", m_num_objects );
    res += fmt::format( "  max_num_objects_per_node {}\n", m_max_num_objects_per_node );
    res += fmt::format( "  bbox_long_edge_ratio     {}\n", m_bbox_long_edge_ratio );
    res += fmt::format( "  bbox_overlap_tolerance   {}\n", m_bbox_overlap_tolerance );
    res += "--------------------------------\n";
    return res;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  void
  AABBtree<Real>::set_max_num_objects_per_node( integer n ) {
    UTILS_ASSERT(
      n > 0 && n <= 4096,
      "AABBtree::set_max_num_objects_per_node( nobj = {} )\n"
      "nobj must be > 0 and <= 4096\n",
      n
    );
    m_max_num_objects_per_node = n;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  void
  AABBtree<Real>::set_bbox_long_edge_ratio( Real ratio ) {
    UTILS_ASSERT(
      ratio > 0 && ratio < 1,
      "AABBtree::set_bbox_long_edge_ratio( ratio = {} )\n"
      "tol must be > 0 and < 1\n",
      ratio
    );
    m_bbox_long_edge_ratio = ratio;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  void
  AABBtree<Real>::set_bbox_overlap_tolerance( Real tol ) {
    UTILS_ASSERT(
      tol > 0 && tol < 1,
      "AABBtree::set_bbox_overlap_tolerance( tol = {} )\n"
      "tol must be > 0 and < 1\n",
      tol
    );
    m_bbox_overlap_tolerance = tol;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  void
  AABBtree<Real>::set_bbox_min_size_tolerance( Real tol ) {
    UTILS_ASSERT(
      tol >= 0,
      "AABBtree::set_bbox_min_size_tolerance( tol = {} )\n"
      "tol must be >= 0\n",
      tol
    );
    m_bbox_min_size_tolerance = tol;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  void
  AABBtree<Real>::allocate( integer nbox, integer dim ) {

    UTILS_WARNING(
      dim <= 10,
      "AABBtree::allocate( nbox, dim={} )\n"
      "dim is greather that 10!!!",
      dim
    );

    switch ( dim ) {
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
    m_2dim        = 2*dim;
    m_num_objects = nbox;
    m_nmax        = 2*m_num_objects; // estimate max memory usage

    m_rmem.allocate( size_t((m_nmax+m_num_objects)*m_2dim) );
    m_imem.allocate( size_t(4*m_nmax+m_num_objects) );

    m_bbox_tree = m_rmem( size_t(m_nmax*m_2dim) );
    m_bbox_objs = m_rmem( size_t(m_num_objects*m_2dim) );

    m_father    = m_imem( size_t(m_nmax) );
    m_child     = m_imem( size_t(m_nmax) );
    m_ptr_nodes = m_imem( size_t(m_nmax) );
    m_num_nodes = m_imem( size_t(m_nmax) );
    m_id_nodes  = m_imem( size_t(m_num_objects) );

    // initialize id nodes, will be reordered during the tree build
    for ( integer i = 0; i < m_num_objects; ++i ) m_id_nodes[i] = i;

    // setup root node
    m_father[0]    = -1;
    m_child[0]     = -1;
    m_ptr_nodes[0] = 0;
    m_num_nodes[0] = m_num_objects;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  void
  AABBtree<Real>::add_bboxes(
    Real const bbox_min[], integer ldim0,
    Real const bbox_max[], integer ldim1
  ) {

    UTILS_ASSERT(
      ldim0 >= m_dim && ldim1 >= m_dim,
      "AABBtree::add_bboxes( bb_min, ldim0={}, bb_max, ldim1={} )\n"
      "must be ldim0, ldim1 >= dim = {}\n",
      ldim0, ldim1, m_dim
    );

    Real * bb = m_bbox_objs;
    for ( integer i = 0; i < m_num_objects; ++i ) {
      for ( integer j = 0; j < m_dim; ++j ) {
        UTILS_ASSERT(
          bbox_min[j] <= bbox_max[j],
          "AABBtree::add_bboxes, bad bbox N.{} max < min", i
        );
      }
      std::copy_n( bbox_min, m_dim, bb ); bb += m_dim; bbox_min += ldim0;
      std::copy_n( bbox_max, m_dim, bb ); bb += m_dim; bbox_max += ldim1;
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  void
  AABBtree<Real>::replace_bbox(
    Real const bbox_min[],
    Real const bbox_max[],
    integer    ipos
  ) {
    UTILS_ASSERT(
      ipos >= 0 && ipos < m_num_objects,
      "AABBtree::replace_bbox( bb_min, bb_max, ipos = {})"
      " ipos must be in [0,{})\n",
      ipos, m_num_objects
    );
    Real * bb = m_bbox_objs + ipos*m_2dim;
    for ( integer j = 0; j < m_dim; ++j ) {
      UTILS_ASSERT(
        bbox_min[j] <= bbox_max[j],
        "AABBtree::replace_bbox, bad bbox N.{} max < min", ipos
      );
    }
    std::copy_n( bbox_min, m_dim, bb );
    std::copy_n( bbox_max, m_dim, bb+m_dim );
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  void
  AABBtree<Real>::build() {

    Real otol = pow( m_bbox_overlap_tolerance, m_dim );

    // root contains all rectangles, build its bbox
    for ( integer j = 0; j < m_dim; ++j ) {
      Real & minj = m_bbox_tree[j];
      Real & maxj = m_bbox_tree[m_dim+j];
      Real const * pmin = m_bbox_objs+j;
      Real const * pmax = m_bbox_objs+j+m_dim;
      minj = *pmin;
      maxj = *pmax;
      UTILS_ASSERT0(
        maxj >= minj,
        "AABBtree::build, bad bbox N.0 max < min"
      );
      for ( integer i = 1; i < m_num_objects; ++i ) {
        pmin += m_2dim;
        pmax += m_2dim;
        UTILS_ASSERT(
          *pmax >= *pmin,
          "AABBtree::build, bad bbox N.{} max < min ({} < {})\n",
          i, *pmax, *pmin
        );
        if ( minj > *pmin ) minj = *pmin;
        if ( maxj < *pmax ) maxj = *pmax;
      }
    }

    // main loop: divide nodes until all constraints satisfied
    m_stack.clear();
    m_stack.reserve(2*m_num_objects+1);
    m_stack.emplace_back(0);
    m_num_tree_nodes = 1;

    while ( !m_stack.empty() ) {

      // pop node from stack
      integer id_father = m_stack.back(); m_stack.pop_back();
      UTILS_ASSERT_DEBUG(
        id_father < m_nmax,
        "AABBtree::build, id_father = {} must be less than m_nmax ={}\n",
        id_father, m_nmax
      );

      // set no childer for the moment
      m_child[id_father] = -1;

      // get rectangles id in parent
      integer num = m_num_nodes[id_father];

      // if few bbox stop splitting
      if ( num < m_max_num_objects_per_node ) continue;

      integer  iptr = m_ptr_nodes[id_father];
      integer * ptr = m_id_nodes + iptr;

      // split plane on longest axis, use euristic
      Real const * father_min = m_bbox_tree + id_father * m_2dim;
      Real const * father_max = father_min + m_dim;

      integer idim = 0;
      Real    mx   = father_max[0] - father_min[0];
      for ( integer i = 1; i < m_dim; ++i ) {
        Real mx1 = father_max[i] - father_min[i];
        if ( mx < mx1 ) { mx = mx1; idim = i; }
      }

      // if too small bbox stop splitting
      if ( mx < m_bbox_min_size_tolerance ) continue;

      Real tol_len = m_bbox_long_edge_ratio * mx;
      Real sp      = 0;

      // separate short/long and accumulate short baricenter
      integer n_long  = 0;
      integer n_short = 0;
      while ( n_long + n_short < num ) {
        integer id = ptr[n_long];
        UTILS_ASSERT_DEBUG(
          id < m_num_objects,
          "AABBtree::build, id = {} must be less than m_num_objects ={}\n",
          id, m_num_objects
        );
        Real const * id_min = m_bbox_objs + id * m_2dim;
        Real const * id_max = id_min + m_dim;
        Real id_len = id_max[idim] - id_min[idim];
        if ( id_len > tol_len ) {
          // found long BBOX, increment n_long and update position
          ++n_long;
        } else {
          // found short BBOX, increment n_short and exchange with bottom
          ++n_short;
          swap( ptr[n_long], ptr[num-n_short] );
          sp += id_max[idim] + id_min[idim];
        }
      }

      // if split rectangles do not improve search, stop split at this level
      if ( n_short < 2 ) continue;

      // select the split position: take the mean of the set of
      // (non-"long") rectangle centers along axis idim
      sp /= 2*n_short;

      // partition based on centers
      integer n_left  = 0;
      integer n_right = 0;

      while ( n_long + n_left + n_right < num ) {
        integer id = ptr[n_long+n_left];
        Real const * id_min = m_bbox_objs + id * m_2dim;
        Real const * id_max = id_min + m_dim;
        Real id_mid = (id_max[idim] + id_min[idim])/2;
        if ( id_mid < sp ) {
          ++n_left; // in right position do nothing
        } else {
          ++n_right;
          swap( ptr[n_long+n_left], ptr[num-n_right] );
        }
      }

      // if cannot improve bbox, stop split at this level!
      if ( n_left == 0 || n_right == 0 ) continue;

      // child indexing
      integer id_left  = m_num_tree_nodes + 0;
      integer id_right = m_num_tree_nodes + 1;

      UTILS_ASSERT_DEBUG(
        id_right < m_nmax,
        "AABBtree::build, id_right = {} must be less than m_nmax ={}\n",
        id_right, m_nmax
      );

      // compute bbox of left and right child
      Real * bb_left_min = m_bbox_tree + id_left * m_2dim;
      Real * bb_left_max = bb_left_min + m_dim;
      for ( integer i = 0; i < n_left; ++i ) {
        integer id = ptr[n_long+i];
        UTILS_ASSERT_DEBUG(
          id < m_num_objects,
          "AABBtree::build, id = {} must be less than m_num_objects ={}\n",
          id, m_num_objects
        );
        Real const * bb_id_min = m_bbox_objs + id * m_2dim;
        Real const * bb_id_max = bb_id_min + m_dim;
        if ( i == 0 ) {
          copy_n( bb_id_min, m_2dim, bb_left_min );
          //copy_n( bb_id_max, m_dim, left_max );
        } else {
          for ( integer j = 0; j < m_dim; ++j ) {
            if ( bb_left_min[j] > bb_id_min[j] ) bb_left_min[j] = bb_id_min[j];
            if ( bb_left_max[j] < bb_id_max[j] ) bb_left_max[j] = bb_id_max[j];
          }
        }
      }

      Real * bb_right_min = m_bbox_tree + id_right * m_2dim;
      Real * bb_right_max = bb_right_min + m_dim;
      for ( integer i = 0; i < n_right; ++i ) {
        integer id = ptr[n_long+n_left+i];
        UTILS_ASSERT_DEBUG(
          id < m_num_objects,
          "AABBtree::build, id = {} must be less than m_num_objects ={}\n",
          id, m_num_objects
        );
        Real const * bb_id_min = m_bbox_objs + id * m_2dim;
        Real const * bb_id_max = bb_id_min + m_dim;
        if ( i == 0 ) {
          copy_n( bb_id_min, m_2dim, bb_right_min );
          //copy_n( bb_id_max, m_dim, bb_right_max );
        } else {
          for ( integer j = 0; j < m_dim; ++j ) {
            if ( bb_right_min[j] > bb_id_min[j] ) bb_right_min[j] = bb_id_min[j];
            if ( bb_right_max[j] < bb_id_max[j] ) bb_right_max[j] = bb_id_max[j];
          }
        }
      }

      // check again if split improve the AABBtree otherwise stop exploration
      if ( n_left < m_max_num_objects_per_node || n_right < m_max_num_objects_per_node ) {
        // few nodes, check if improve volume
        Real vo = 1;
        Real vL = 1;
        Real vR = 1;
        for ( integer j = 0l; j < m_dim; ++j ) {
          Real Lmin = bb_left_min[j];
          Real Lmax = bb_left_max[j];
          Real Rmin = bb_right_min[j];
          Real Rmax = bb_right_max[j];
          vo *= max(min(Lmax,Rmax) - max(Lmin,Rmin), Real(0));
          vL *= Lmax - Lmin;
          vR *= Rmax - Rmin;
        }
        // if do not improve volume, stop split at this level!
        if ( vo > (vL+vR-vo)*otol ) continue;
      }

      // push child nodes onto stack
      m_father[id_left]  = id_father;
      m_father[id_right] = id_father;
      m_child[id_father] = id_left;

      m_num_nodes[id_father] = n_long;

      m_ptr_nodes[id_left]  = iptr + n_long;
      m_num_nodes[id_left]  = n_left;

      m_ptr_nodes[id_right] = iptr + n_long + n_left;
      m_num_nodes[id_right] = n_right;

      m_stack.emplace_back(id_left);
      m_stack.emplace_back(id_right);
      m_num_tree_nodes += 2;
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  void
  AABBtree<Real>::intersect_with_one_point(
    Real const pnt[],
    AABB_SET & bb_index
  ) const {

    m_num_check = 0;

    // quick return on empty inputs
    if ( m_num_tree_nodes == 0 ) return;

    // descend tree from root
    m_stack.clear();
    m_stack.reserve(2*m_num_tree_nodes+1);
    m_stack.emplace_back(0);
    while ( !m_stack.empty() ) {
      // pop node from stack
      integer id_father = m_stack.back(); m_stack.pop_back();

      // get BBOX
      Real const * bb_father = m_bbox_tree + id_father * m_2dim;

      ++m_num_check;
      bool overlap = m_check_overlap_with_point( pnt, bb_father, m_dim );

      // if do not overlap skip
      if ( !overlap ) continue;

      // get rectangles id in parent
      this->get_bbox_indexes_of_a_node( id_father, bb_index );

      integer nn = m_child[id_father];
      if ( nn > 0 ) { // root == 0, children > 0
        // push on stack children
        m_stack.emplace_back(nn);
        m_stack.emplace_back(nn+1);
      }
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  void
  AABBtree<Real>::intersect_with_one_point_and_refine(
    Real const pnt[],
    AABB_SET & bb_index
  ) const {

    m_num_check = 0;

    // quick return on empty inputs
    if ( m_num_tree_nodes == 0 ) return;

    // descend tree from root
    m_stack.clear();
    m_stack.reserve(2*m_num_tree_nodes+1);
    m_stack.emplace_back(0);
    while ( !m_stack.empty() ) {
      // pop node from stack
      integer id_father = m_stack.back(); m_stack.pop_back();

      // get BBOX
      Real const * bb_father = m_bbox_tree + id_father * m_2dim;

      ++m_num_check;
      bool overlap = m_check_overlap_with_point( pnt, bb_father, m_dim );

      // if do not overlap skip
      if ( !overlap ) continue;

      // refine candidate
      integer num = this->m_num_nodes[id_father];
      integer const * ptr = this->m_id_nodes + this->m_ptr_nodes[id_father];
      for ( integer ii = 0; ii < num; ++ii ) {
        integer s = ptr[ii];
        Real const * bb_s = m_bbox_objs + s * m_2dim;
        ++m_num_check;
        bool olap = m_check_overlap_with_point( pnt, bb_s, m_dim );
        if ( olap ) bb_index.insert(s);
      }

      integer nn = m_child[id_father];
      if ( nn > 0 ) { // root == 0, children > 0
        // push on stack children
        m_stack.emplace_back(nn);
        m_stack.emplace_back(nn+1);
      }
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  void
  AABBtree<Real>::intersect_with_one_bbox(
    Real const bbox[],
    AABB_SET & bb_index
  ) const {
    m_num_check = 0;

    // quick return on empty inputs
    if ( m_num_tree_nodes == 0 ) return;

    // descend tree from root
    m_stack.clear();
    m_stack.reserve(2*m_num_tree_nodes+1);
    m_stack.emplace_back(0);
    while ( !m_stack.empty() ) {
      // pop node from stack
      integer id_father = m_stack.back(); m_stack.pop_back();

      // get BBOX
      Real const * bb_father = m_bbox_tree + id_father * m_2dim;

      ++m_num_check;
      bool overlap = m_check_overlap( bb_father, bbox, m_dim );

      // if do not overlap skip
      if ( !overlap ) continue;

      // get rectangles id in parent
      this->get_bbox_indexes_of_a_node( id_father, bb_index );

      integer nn = m_child[id_father];
      if ( nn > 0 ) { // root == 0, children > 0
        // push on stack children
        m_stack.emplace_back(nn);
        m_stack.emplace_back(nn+1);
      }
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  void
  AABBtree<Real>::intersect_with_one_bbox_and_refine(
    Real const bbox[],
    AABB_SET & bb_index
  ) const {

    m_num_check = 0;

    // quick return on empty inputs
    if ( m_num_tree_nodes == 0 ) return;

    // descend tree from root
    m_stack.clear();
    m_stack.reserve(2*m_num_tree_nodes+1);
    m_stack.emplace_back(0);
    while ( !m_stack.empty() ) {
      // pop node from stack
      integer id_father = m_stack.back(); m_stack.pop_back();

      // get BBOX
      Real const * bb_father = m_bbox_tree + id_father * m_2dim;

      ++m_num_check;
      bool overlap = m_check_overlap( bb_father, bbox, m_dim );

      // if do not overlap skip
      if ( !overlap ) continue;

      // refine candidate
      integer num = this->m_num_nodes[id_father];
      integer const * ptr = this->m_id_nodes + this->m_ptr_nodes[id_father];
      for ( integer ii = 0; ii < num; ++ii ) {
        integer s = ptr[ii];
        Real const * bb_s = m_bbox_objs + ptr[ii] * m_2dim;
        ++m_num_check;
        bool olap = m_check_overlap( bb_s, bbox, m_dim );
        if ( olap ) bb_index.insert(s);
      }

      integer nn = m_child[id_father];
      if ( nn > 0 ) { // root == 0, children > 0
        // push on stack children
        m_stack.emplace_back(nn);
        m_stack.emplace_back(nn+1);
      }
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  void
  AABBtree<Real>::intersect(
    AABBtree<Real> const & aabb,
    AABB_MAP             & bb_index
  ) const {

    m_num_check = 0;

    // quick return on empty inputs
    if ( this->m_num_tree_nodes == 0 || aabb.m_num_tree_nodes == 0 ) return;

    // descend tree from root
    m_stack.clear();
    m_stack.reserve(m_num_tree_nodes+aabb.m_num_tree_nodes+2);
    m_stack.emplace_back(0);
    m_stack.emplace_back(0);
    integer n_stack = 2;
    while ( !m_stack.empty() ) {
      // pop node from stack
      integer root2  = m_stack.back(); m_stack.pop_back();
      integer sroot1 = m_stack.back(); m_stack.pop_back();
      integer root1  = sroot1 >= 0 ? sroot1 : -1-sroot1;

      // check for intersection
      Real const * bb_root1 = this->m_bbox_tree + root1 * m_2dim;
      Real const * bb_root2 = aabb.m_bbox_tree + root2 * m_2dim;

      ++m_num_check;
      bool overlap = m_check_overlap( bb_root1, bb_root2, m_dim );

      // if do not overlap skip
      if ( !overlap ) continue;

      // check if there are elements to check
      integer nn1 = this->m_num_nodes[root1];
      integer nn2 = aabb.m_num_nodes[root2];
      if ( nn1 > 0 && nn2 > 0 ) aabb.get_bbox_indexes_of_a_node( root2, bb_index[root1] );

      integer id_lr1 = sroot1 >= 0 ? m_child[root1] : -1;
      integer id_lr2 = aabb.m_child[root2];

      UTILS_ASSERT_DEBUG(
        n_stack+3 < 2*m_nmax,
        "AABBtree::intersect, n_stack+2 = {} must be less than 2*m_nmax ={}\n",
        n_stack+2, 2*m_nmax
      );

      if ( id_lr1 >= 0 ) {
        m_stack.emplace_back(id_lr1);   m_stack.emplace_back(root2);
        m_stack.emplace_back(id_lr1+1); m_stack.emplace_back(root2);
        if ( nn1 > 0 ) {
          m_stack.emplace_back(-1-root1); m_stack.emplace_back(root2);
        }
      } else if ( id_lr2 >= 0 ) {
        m_stack.emplace_back(sroot1); m_stack.emplace_back(id_lr2);
        m_stack.emplace_back(sroot1); m_stack.emplace_back(id_lr2+1);
      }
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  void
  AABBtree<Real>::intersect_and_refine(
    AABBtree<Real> const & aabb,
    AABB_MAP             & bb_index
  ) const {

    m_num_check = 0;

    // quick return on empty inputs
    if ( this->m_num_tree_nodes == 0 || aabb.m_num_tree_nodes == 0 ) return;

    // descend tree from root
    m_stack.clear();
    m_stack.reserve(m_num_tree_nodes+aabb.m_num_tree_nodes+2);
    m_stack.emplace_back(0);
    m_stack.emplace_back(0);
    while ( !m_stack.empty() ) {
      // pop node from stack
      integer root2  = m_stack.back(); m_stack.pop_back();
      integer sroot1 = m_stack.back(); m_stack.pop_back();
      integer root1  = sroot1 >= 0 ? sroot1 : -1-sroot1;

      // check for intersection
      Real const * bb_root1 = this->m_bbox_tree + root1 * m_2dim;
      Real const * bb_root2 = aabb.m_bbox_tree + root2 * m_2dim;

      ++m_num_check;
      bool overlap = m_check_overlap( bb_root1, bb_root2, m_dim );

      // if do not overlap skip
      if ( !overlap ) continue;

      // check if there are elements to check
      integer nn1 = this->m_num_nodes[root1];
      integer nn2 = aabb.m_num_nodes[root2];
      if ( nn1 > 0 && nn2 > 0 ) {
        // construct list of intersecting candidated
        integer const * ptr1 = this->m_id_nodes + this->m_ptr_nodes[root1];
        integer const * ptr2 = aabb.m_id_nodes + aabb.m_ptr_nodes[root2];
        for ( integer ii = 0; ii < nn1; ++ii ) {
          integer s1 = ptr1[ii];
          Real const * bb_s1 = m_bbox_objs + s1 * m_2dim;
          AABB_SET & BB = bb_index[s1];
          for ( integer jj = 0; jj < nn2; ++jj ) {
            integer s2 = ptr2[jj];
            Real const * bb_s2 = aabb.m_bbox_objs + s2 * m_2dim;
            ++m_num_check;
            bool olap = m_check_overlap( bb_s1, bb_s2, m_dim );
            //if ( olap ) bb_index[s1].insert(s2);
            if ( olap ) BB.insert(s2);
          }
        }
      }

      integer id_lr1 = sroot1 >= 0 ? m_child[root1] : -1;
      integer id_lr2 = aabb.m_child[root2];

      if ( id_lr1 >= 0 ) {
        m_stack.emplace_back(id_lr1);   m_stack.emplace_back(root2);
        m_stack.emplace_back(id_lr1+1); m_stack.emplace_back(root2);
        if ( nn1 > 0 ) {
          m_stack.emplace_back(-1-root1); m_stack.emplace_back(root2);
        }
      } else if ( id_lr2 >= 0 ) {
        m_stack.emplace_back(sroot1); m_stack.emplace_back(id_lr2);
        m_stack.emplace_back(sroot1); m_stack.emplace_back(id_lr2+1);
      }
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  void
  AABBtree<Real>::pnt_bbox_minmax(
    Real const pnt[],
    Real const bbox[],
    Real     & dmin2,
    Real     & dmax2
  ) const {
    Real const * bb_max = bbox+m_dim;
    Real const * bb_min = bbox;
    dmin2 = 0;
    dmax2 = 0;
    for ( integer i = 0; i < m_dim; ++i ) {
      // check overlap
      Real pi    = pnt[i];
      Real dpmin = 0;
      Real dpmax = 0;
      Real t1    = pi-bb_max[i];
      Real t2    = bb_min[i]-pi;
      if ( t1 > 0 ) dpmin = t1; else dpmax = t1;
      if ( t2 > 0 ) dpmin = t2; else dpmax = t2;
      dmin2 += dpmin*dpmin;
      dmax2 += dpmax*dpmax;
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  void
  AABBtree<Real>::min_distance_candidates(
    Real const pnt[],
    AABB_SET & bb_index
  ) const {;

    Real dst2_min, dst2_max;

    // quick return on empty inputs
    bb_index.clear();
    if ( this->m_num_tree_nodes == 0 ) return;

    Real min_max_distance2 = Utils::Inf<Real>();

    // descend tree from root
    m_stack.clear();
    m_stack.reserve(m_num_tree_nodes+1);
    m_stack.emplace_back(0);
    while ( !m_stack.empty() ) {
      // pop node from stack
      integer id_father = m_stack.back(); m_stack.pop_back();

      // get BBOX
      // check for intersection
      Real const * father_bbox = this->m_bbox_tree + id_father * m_2dim;
      this->pnt_bbox_minmax( pnt, father_bbox, dst2_min, dst2_max );

      if ( dst2_min <= min_max_distance2 ) {

        if ( m_num_nodes[id_father] > 0 && dst2_max < min_max_distance2 )
          min_max_distance2 = dst2_max;

        integer nn = m_child[id_father];
        if ( nn > 0 ) { // root == 0, children > 0
          // push on stack childrens
          m_stack.emplace_back(nn);
          m_stack.emplace_back(nn+1);
        }
      }
    }

    // descend tree from root
    m_stack.clear();
    m_stack.emplace_back(0);
    while ( !m_stack.empty() ) {
      // pop node from stack
      integer id_father = m_stack.back(); m_stack.pop_back();
      Real const * father_bbox = this->m_bbox_tree + id_father * m_2dim;
      this->pnt_bbox_minmax( pnt, father_bbox, dst2_min, dst2_max );
      if ( dst2_min <= min_max_distance2 ) {
        this->get_bbox_indexes_of_a_node( id_father, bb_index );
        integer nn = m_child[id_father];
        if ( nn > 0 ) { // root == 0, children > 0
          // push on stack childrens
          m_stack.emplace_back(nn);
          m_stack.emplace_back(nn+1);
        }
      }
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  void
  AABBtree<Real>::get_bbox_indexes_of_a_node(
    integer    i_pos,
    AABB_SET & bb_index
  ) const {
    UTILS_ASSERT(
      i_pos >= 0 && i_pos < m_num_tree_nodes,
      "AABBtree::get_bbox_indexes_of_a_node( i_pos={}, bb_index ) i_pos must be >= 0 and < {}\n",
      i_pos, m_num_tree_nodes
    );
    integer num = m_num_nodes[i_pos];
    integer ptr = m_ptr_nodes[i_pos];
    while ( num-- > 0 ) bb_index.insert( m_id_nodes[ptr++] );
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  typename AABBtree<Real>::integer
  AABBtree<Real>::num_tree_nodes( integer nmin ) const {
    integer n = 0;
    for ( integer i = 0; i < m_num_tree_nodes; ++i )
      if ( m_num_nodes[i] >= nmin ) ++n;
    return n;
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  void
  AABBtree<Real>::get_root_bbox( Real bb_min[], Real bb_max[] ) const {
    std::copy_n( m_bbox_tree,       m_dim, bb_min );
    std::copy_n( m_bbox_tree+m_dim, m_dim, bb_max );
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template <typename Real>
  void
  AABBtree<Real>::get_bboxes_of_the_tree(
    Real bbox_min[], integer ldim0,
    Real bbox_max[], integer ldim1,
    integer nmin
  ) const {
    UTILS_ASSERT(
      ldim0 >= m_dim && ldim1 >= m_dim,
      "AABBtree::get_bboxes_of_the_tree(\n"
      "  bbox_min, ldim0={},\n"
      "  bbox_max, ldim1={},\n"
      "  nmin={} )\n"
      "must be nmin >= 0 and ldim0:1 >= {}\n",
      ldim0, ldim1, nmin, m_dim
    );

    for ( integer i = 0; i < m_num_tree_nodes; ++i ) {
      if ( m_num_nodes[i] >= nmin ) {
        Real const * b_min = m_bbox_tree + i * m_2dim;
        Real const * b_max = b_min + m_dim;
        std::copy_n( b_min, m_dim, bbox_min ); bbox_min += ldim0;
        std::copy_n( b_max, m_dim, bbox_max ); bbox_max += ldim1;
      }
    }
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  template class AABBtree<float>;
  template class AABBtree<double>;

}

///
/// eof: Utils_AABB_tree.cc
///
