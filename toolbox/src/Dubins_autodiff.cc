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
#define MAP_AUTODIFF_TO_METHOD( CLASS, NAME )                            \
  autodiff::dual1st CLASS::NAME( autodiff::dual1st const & s ) const    \
  {                                                                      \
    real_type         V = NAME( s.val );                                 \
    real_type         D = NAME##_D( s.val );                             \
    autodiff::dual1st res;                                               \
    res.val  = V;                                                        \
    res.grad = s.grad * D;                                               \
    return res;                                                          \
  }                                                                      \
                                                                         \
  autodiff::dual2nd CLASS::NAME( autodiff::dual2nd const & s ) const    \
  {                                                                      \
    real_type         V  = NAME( s.val.val );                            \
    real_type         D  = NAME##_D( s.val.val );                        \
    real_type         DD = NAME##_DD( s.val.val );                       \
    autodiff::dual2nd res;                                               \
    res.val.val   = V;                                                   \
    res.val.grad  = D * s.val.grad;                                      \
    res.grad.val  = D * s.grad.val;                                      \
    res.grad.grad = DD * s.val.grad * s.grad.val + D * s.grad.grad;      \
    return res;                                                          \
  }
  
  MAP_AUTODIFF_TO_METHOD( Dubins, theta );
  MAP_AUTODIFF_TO_METHOD( Dubins, X );
  MAP_AUTODIFF_TO_METHOD( Dubins, Y );

  MAP_AUTODIFF_TO_METHOD( Dubins3p, theta );
  MAP_AUTODIFF_TO_METHOD( Dubins3p, X );
  MAP_AUTODIFF_TO_METHOD( Dubins3p, Y );


#endif

}  // namespace G2lib

// EOF: Fresnel.cc
