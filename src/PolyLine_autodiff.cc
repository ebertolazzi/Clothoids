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
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Clothoids.hh"
#include "Clothoids_fmt.hh"

namespace G2lib
{
#ifdef AUTODIFF_SUPPORT

#define MAP_AUTODIFF_TO_METHOD( NAME )                                  \
  autodiff::dual1st PolyLine::NAME( autodiff::dual1st const & s ) const \
  {                                                                     \
    real_type         sv  = s.val;                                      \
    integer const     idx = find_at_s( sv );                            \
    autodiff::dual1st ss  = s;                                          \
    ss.val                = sv - m_s0[idx];                             \
    return m_polyline_list[size_t( idx )].NAME( ss );                   \
  }                                                                     \
                                                                        \
  autodiff::dual2nd PolyLine::NAME( autodiff::dual2nd const & s ) const \
  {                                                                     \
    real_type         sv  = s.val.val;                                  \
    integer const     idx = find_at_s( sv );                            \
    autodiff::dual2nd ss  = s;                                          \
    ss.val.val            = sv - m_s0[idx];                             \
    return m_polyline_list[size_t( idx )].NAME( ss );                   \
  }

  MAP_AUTODIFF_TO_METHOD( theta );
  MAP_AUTODIFF_TO_METHOD( X );
  MAP_AUTODIFF_TO_METHOD( Y );

#endif

}  // namespace G2lib

// EOF: PolyLine_autodiff.cc
