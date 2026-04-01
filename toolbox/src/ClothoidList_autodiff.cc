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

#define MAP_AUTODIFF_TO_METHOD( NAME )                                      \
  autodiff::dual1st ClothoidList::NAME( autodiff::dual1st const & s ) const \
  {                                                                         \
    real_type             sv  = s.val;                                      \
    integer const         idx = find_at_s( sv );                            \
    ClothoidCurve const & c   = get( idx );                                 \
    autodiff::dual1st     ss  = s;                                          \
    ss.val                    = sv - m_s0[idx];                             \
    return c.NAME( ss );                                                    \
  }                                                                         \
                                                                            \
  autodiff::dual2nd ClothoidList::NAME( autodiff::dual2nd const & s ) const \
  {                                                                         \
    real_type             sv  = s.val.val;                                  \
    integer const         idx = find_at_s( sv );                            \
    ClothoidCurve const & c   = get( idx );                                 \
    autodiff::dual2nd     ss  = s;                                          \
    ss.val.val                = sv - m_s0[idx];                             \
    return c.NAME( ss );                                                    \
  }

#define MAP_AUTODIFF_TO_METHOD2( NAME )                                                           \
  autodiff::dual1st ClothoidList::NAME( autodiff::dual1st const & s, real_type const offs ) const \
  {                                                                                               \
    real_type             sv  = s.val;                                                            \
    integer const         idx = find_at_s( sv );                                                  \
    ClothoidCurve const & c   = get( idx );                                                       \
    autodiff::dual1st     ss  = s;                                                                \
    ss.val                    = sv - m_s0[idx];                                                   \
    return c.NAME( ss, offs );                                                                    \
  }                                                                                               \
                                                                                                  \
  autodiff::dual2nd ClothoidList::NAME( autodiff::dual2nd const & s, real_type const offs ) const \
  {                                                                                               \
    real_type             sv  = s.val.val;                                                        \
    integer const         idx = find_at_s( sv );                                                  \
    ClothoidCurve const & c   = get( idx );                                                       \
    autodiff::dual2nd     ss  = s;                                                                \
    ss.val.val                = sv - m_s0[idx];                                                   \
    return c.NAME( ss, offs );                                                                    \
  }

  MAP_AUTODIFF_TO_METHOD( theta );
  MAP_AUTODIFF_TO_METHOD( tx );
  MAP_AUTODIFF_TO_METHOD( ty );
  MAP_AUTODIFF_TO_METHOD( X );
  MAP_AUTODIFF_TO_METHOD( Y );
  MAP_AUTODIFF_TO_METHOD2( X_ISO );
  MAP_AUTODIFF_TO_METHOD2( Y_ISO );

#endif

}  // namespace G2lib

// EOF: ClothoidList.cc
