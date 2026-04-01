/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2026                                                      |
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

// Standard autodiff mapping for methods that have value and derivatives
#define MAP_AUTODIFF_TO_METHOD( NAME )                                   \
  autodiff::dual1st BaseCurve::NAME( autodiff::dual1st const & s ) const \
  {                                                                      \
    real_type         V = NAME( s.val );                                 \
    real_type         D = NAME##_D( s.val );                             \
    autodiff::dual1st res;                                               \
    res.val  = V;                                                        \
    res.grad = s.grad * D;                                               \
    return res;                                                          \
  }                                                                      \
                                                                         \
  autodiff::dual2nd BaseCurve::NAME( autodiff::dual2nd const & s ) const \
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

// Standard autodiff mapping for methods that have value and derivatives
#define MAP_AUTODIFF_TO_METHOD2( NAME )                                                        \
  autodiff::dual1st BaseCurve::NAME( autodiff::dual1st const & s, real_type const offs ) const \
  {                                                                                            \
    real_type         V = NAME( s.val, offs );                                                 \
    real_type         D = NAME##_D( s.val, offs );                                             \
    autodiff::dual1st res;                                                                     \
    res.val  = V;                                                                              \
    res.grad = s.grad * D;                                                                     \
    return res;                                                                                \
  }                                                                                            \
                                                                                               \
  autodiff::dual2nd BaseCurve::NAME( autodiff::dual2nd const & s, real_type const offs ) const \
  {                                                                                            \
    real_type         V  = NAME( s.val.val, offs );                                            \
    real_type         D  = NAME##_D( s.val.val, offs );                                        \
    real_type         DD = NAME##_DD( s.val.val, offs );                                       \
    autodiff::dual2nd res;                                                                     \
    res.val.val   = V;                                                                         \
    res.val.grad  = D * s.val.grad;                                                            \
    res.grad.val  = D * s.grad.val;                                                            \
    res.grad.grad = DD * s.val.grad * s.grad.val + D * s.grad.grad;                            \
    return res;                                                                                \
  }

  // Methods with derivatives available
  MAP_AUTODIFF_TO_METHOD( theta );
  MAP_AUTODIFF_TO_METHOD( kappa );
  MAP_AUTODIFF_TO_METHOD( tx );
  MAP_AUTODIFF_TO_METHOD( ty );
  MAP_AUTODIFF_TO_METHOD( nx_ISO );
  MAP_AUTODIFF_TO_METHOD( ny_ISO );
  MAP_AUTODIFF_TO_METHOD( nx_SAE );
  MAP_AUTODIFF_TO_METHOD( ny_SAE );
  MAP_AUTODIFF_TO_METHOD( X );
  MAP_AUTODIFF_TO_METHOD( Y );

  // ISO/SAE methods (without offset parameter)
  MAP_AUTODIFF_TO_METHOD2( X_ISO );
  MAP_AUTODIFF_TO_METHOD2( Y_ISO );
  MAP_AUTODIFF_TO_METHOD2( X_SAE );
  MAP_AUTODIFF_TO_METHOD2( Y_SAE );

#endif

}  // namespace G2lib

// EOF: Fresnel.cc
